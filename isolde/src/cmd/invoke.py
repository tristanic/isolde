# @Author: Tristan Croll
# @Date:   27-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 27-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Invoke an agent-safe ISOLDE command from a validated JSON arguments object.

Rather than reimplement ChimeraX's argument parser (atom-spec evaluation, enum
coercion, default-filling), we **rebuild the command line text and run() it**.
The agent supplies atom-specs as plain strings, so they pass straight through to
ChimeraX's own parser -- no lossy object round-tripping, and behaviour is
identical to a human typing the command.

The single choke point :func:`invoke_command` also threads the post-condition
*witness* (design decision 5): it snapshots a pre-condition, runs the command,
then measures what actually changed (e.g. RMSD moved) so an agent can tell "did
work" from "no-op succeeded". Witnesses are keyed by category on the command
record, applied here uniformly -- never bolted onto individual commands.
'''
from typing import Any


def build_command_text(record, json_args: dict) -> str:
    '''
    Serialize a validated JSON args object back into a ChimeraX command line
    (without the leading command name).

    required args  -> positional tokens, in CmdDesc order
    everything else (optional + keyword) -> 'name value' pairs (ChimeraX exposes
        optional args by keyword too, so this avoids fragile positional gaps)
    '''
    from chimerax.core.commands import cli
    desc = record.cmd_desc
    required = dict(getattr(desc, '_required', {}))
    optional = dict(getattr(desc, '_optional', {}))
    keyword = dict(getattr(desc, '_keyword', {}))

    tokens = []

    # Required args: positional, in declared order. Missing required args are a
    # schema-validation failure upstream; we still guard here.
    for name, annotation in required.items():
        if name not in json_args:
            raise ValueError('missing required argument %r for %s' % (name, record.name))
        tokens.append(_format_value(cli, annotation, json_args[name], positional=True))

    # Optional + keyword args: by name. Dedupe names already consumed as required.
    seen = set(required)
    for name, annotation in list(optional.items()) + list(keyword.items()):
        if name in seen or name not in json_args:
            continue
        seen.add(name)
        value = json_args[name]
        formatted = _format_value(cli, annotation, value, positional=False)
        if formatted is _OMIT:
            continue                       # e.g. a false NoArg flag
        if formatted is _BARE_FLAG:
            tokens.append(name)            # NoArg present -> just the keyword
        else:
            tokens.append('%s %s' % (name, formatted))

    return ' '.join(tokens)


# sentinels for flag handling
_OMIT = object()
_BARE_FLAG = object()


def _format_value(cli, annotation, value, positional):
    '''Format a single JSON value as a ChimeraX command-line token.'''
    # NoArg: a value-less flag. True -> emit bare keyword; False -> omit entirely.
    if _is_cls(annotation, getattr(cli, 'NoArg', None)):
        return _BARE_FLAG if value else _OMIT

    # Aggregate (ListOf/SetOf) / RepeatOf / plain list -> comma-joined items
    if isinstance(value, (list, tuple)):
        inner = _inner_annotation(cli, annotation)
        return ','.join(_scalar_token(cli, inner, v) for v in value)

    return _scalar_token(cli, annotation, value)


def _scalar_token(cli, annotation, value):
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if isinstance(value, (int, float)):
        return repr(value) if isinstance(value, float) else str(value)
    # strings: atom-specs, enum values, filenames. Quote if necessary so paths
    # and multi-word values survive the parser.
    s = str(value)
    try:
        return cli.quote_if_necessary(s)
    except Exception:
        return s


def _inner_annotation(cli, annotation):
    if isinstance(annotation, getattr(cli, 'Aggregate', ())):
        return annotation.annotation
    return annotation  # RepeatOf loses inner type; treat items as scalars


def _is_cls(annotation, cls):
    if cls is None:
        return False
    a = annotation if isinstance(annotation, type) else type(annotation)
    try:
        return a is cls or issubclass(a, cls)
    except TypeError:
        return False


def invoke_command(session, name: str, json_args: dict, *, enforce_agent_safe=True) -> dict:
    '''
    Invoke ISOLDE command *name* with *json_args*. Returns a structured result:

        {'command': <text run>, 'result': <jsonified return value>,
         'witness': <post-condition witness or None>, 'log': <captured log>}

    Raises PermissionError if the command is not flagged agent_safe (unless
    enforce_agent_safe=False, for trusted internal callers).
    '''
    from .command_registry import get_command
    from chimerax.core.commands import run

    record = get_command(name)
    if record is None:
        raise KeyError('no such ISOLDE command: %r' % name)
    if enforce_agent_safe and not record.agent_safe:
        raise PermissionError('%r is not exposed to agents (not agent_safe)' % name)

    arg_text = build_command_text(record, json_args)
    command_text = name if not arg_text else '%s %s' % (name, arg_text)

    # Human-facing transparency: a remote agent's commands run under a log
    # capture that hides them from the GUI log (see announce_to_user). For
    # *mutating* commands, echo the action so a human sharing the live session
    # sees what the agent is doing; read-only queries stay quiet (no spam).
    announce = _should_announce(record)
    if announce:
        announce_to_user(session, '%s %s' % (_AGENT_TAG, command_text))

    # Pre-condition witness (best-effort; never blocks the command).
    pre = _witness_pre(session, record, json_args)

    result = run(session, command_text, log=False)

    post = _witness_post(session, record, json_args, pre)
    json_result = _jsonify(result)

    # ...and the outcome the witness observed, so the human sees not just what
    # the agent asked for but what actually happened.
    if announce:
        summary = _outcome_summary(post)
        if summary:
            announce_to_user(session, '%s   → %s' % (_AGENT_TAG, summary))

    # Soft philosophy guardrails (advisory; never block). Detect moves that look
    # like they violate an ISOLDE principle (a large unsettled shift, a clean
    # metric mistaken for correctness, ...) and ride them along in the result;
    # high-severity flags also reach the human's GUI log.
    flags = _philosophy_flags(session, record, post, json_result)

    out = {
        'command': command_text,
        'result': json_result,
        'witness': post,
    }
    if flags:
        out['philosophy_flags'] = flags
    return out


# --- human-facing announcements (REST/MCP collaboration transparency) --------

# The agent's commands always log in full to the *remote* (the captured log
# returned in the response). This controls what the *human* watching the shared
# GUI additionally sees. Read-only categories would just spam the log; every
# other category changes the model / simulation / restraints / view and is
# announced so the human can follow (and, if need be, veto) the agent's moves.
_QUIET_CATEGORIES = frozenset({'query', 'validation'})
_AGENT_TAG = '[agent]'


def _should_announce(record) -> bool:
    cat = getattr(record, 'category', None)
    return bool(getattr(record, 'agent_safe', False)) and cat not in _QUIET_CATEGORIES


def announce_to_user(session, msg: str) -> None:
    '''Write *msg* to the human-facing GUI log even while a REST agent call is in
    flight.

    Agent calls execute inside a ``StringPlainTextLog`` whose
    ``excludes_other_logs`` is True: it redirects ChimeraX log output to the
    remote and *excludes* the GUI log, so without this a human sharing the
    session would see nothing of what the agent does. We briefly clear that
    exclusion so this one message broadcasts to every registered log (the GUI
    included); the captured remote log receives a harmless copy. Best-effort --
    never raises into the command path.
    '''
    logger = getattr(session, 'logger', None)
    if logger is None:
        return
    try:
        from chimerax.core.logger import StringPlainTextLog
        toggled = []
        for log in list(getattr(logger, 'logs', ()) or ()):
            if isinstance(log, StringPlainTextLog) and getattr(log, 'excludes_other_logs', False):
                log.excludes_other_logs = False
                toggled.append(log)
        try:
            logger.info(msg)
        finally:
            for log in toggled:
                log.excludes_other_logs = True
    except Exception:
        pass


def _outcome_summary(witness) -> str:
    '''A short human phrase describing what a witness observed, or '' when there
    is nothing worth announcing. Tolerant of every witness shape (see .witness).'''
    if not isinstance(witness, dict):
        return ''
    # sim_control
    if 'simulation_running' in witness and 'cmd' in witness:
        cmd = witness.get('cmd')
        if cmd == 'start':
            return 'simulation started' if witness.get('started_ok') else 'sim start may have failed'
        if cmd == 'stop':
            # 'stopped' may read False transiently (async teardown); say
            # 'stopping' rather than implying failure.
            return 'simulation stopped' if witness.get('stopped') else 'simulation stopping'
        if cmd:
            return 'sim %s' % cmd
    # restraint_change
    if 'restraint_total_delta' in witness:
        td = witness.get('restraint_total_delta') or 0
        ed = witness.get('restraint_enabled_delta') or 0
        if td > 0:
            return '+%d restraints' % td
        if td < 0:
            return '%d restraints removed' % (-td)
        if ed > 0:
            return '%d restraints enabled' % ed
        if ed < 0:
            return '%d restraints released' % (-ed)
        return 'restraints unchanged (parameters may have been retuned)'
    # convergence (async refinement launch)
    if witness.get('launched'):
        return 'refinement launched (runs in the background)'
    # coord_move
    if 'rmsd_moved' in witness or 'atom_count_changed' in witness:
        if witness.get('atom_count_changed'):
            return 'atom count changed'
        rmsd = witness.get('rmsd_moved')
        if rmsd is None:
            return ''
        if rmsd < 1e-4:
            return ('target set (realized over coming frames)'
                    if witness.get('sim_running') else 'no movement')
        mx = witness.get('max_atom_shift')
        if mx is not None:
            return 'moved RMSD %.2f Å (max %.2f Å)' % (rmsd, mx)
        return 'moved RMSD %.2f Å' % rmsd
    return ''


# --- philosophy guardrails (implemented in .agent_philosophy; lazy/tolerant) -

def _philosophy_flags(session, record, post, json_result):
    '''Evaluate soft guardrails for this command and announce high-severity ones.
    Best-effort: a guardrail must never break the command path.'''
    try:
        from . import agent_philosophy as ap
        flags = ap.evaluate_guardrails(
            category=getattr(record, 'category', None),
            witness=post, result=json_result)
        ap.announce_high_severity(session, flags)
        return flags
    except Exception:
        return []


# --- witness hooks (implemented in .witness; imported lazily & tolerantly) ---

def _witness_pre(session, record, json_args):
    if not getattr(record, 'witness', None):
        return None
    try:
        from .witness import run_pre
        return run_pre(session, record, json_args)
    except Exception:
        return None


def _witness_post(session, record, json_args, pre):
    if not getattr(record, 'witness', None):
        return None
    try:
        from .witness import run_post
        return run_post(session, record, json_args, pre)
    except Exception as e:
        return {'witness_error': str(e)}


# --- result coercion to JSON-safe structures ---

def _jsonify(obj) -> Any:
    '''Coerce a command's return value into JSON-serializable form.'''
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    if isinstance(obj, dict):
        return {str(k): _jsonify(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_jsonify(v) for v in obj]
    # ChimeraX model: report its id string + name.
    if hasattr(obj, 'id_string') and hasattr(obj, 'name'):
        return {'id': obj.id_string, 'name': obj.name}
    # Atoms/Residues/Collections: report a compact summary rather than dumping.
    spec = getattr(obj, 'spec', None)
    if spec is not None:
        out = {'spec': spec}
        try:
            out['count'] = len(obj)
        except TypeError:
            pass
        return out
    # numpy arrays / scalars
    try:
        import numpy
        if isinstance(obj, numpy.ndarray):
            return obj.tolist()
        if isinstance(obj, numpy.generic):
            return obj.item()
    except Exception:
        pass
    # last resort: string form, so a result is never un-serializable.
    return str(obj)
