# @Author: Tristan Croll
# @Date:   27-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 27-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Central registry of ISOLDE commands for the agent/MCP control surface.

ISOLDE has no single in-memory registry of its commands: every command is
registered through ChimeraX's :func:`chimerax.core.commands.register`, called
from ~28 sites scattered across ~13 modules, all funnelled through the single
entry point :func:`chimerax.isolde.cmd.cmd.register_isolde`.

This module turns that funnel into a *stays-current* registry. The trick is
that every one of those call sites does a **function-local** ``from
chimerax.core.commands import register`` at call time (never a module-scope
import that we'd miss). So if we temporarily replace
``chimerax.core.commands.register`` with a recording wrapper for the duration
of ``register_isolde``, every command -- present and future, including ones
that land on as-yet-unmerged branches -- is captured automatically with no edit
to the call sites.

Capture is *not* exposure. A captured command is reachable by an agent only
when it is explicitly flagged ``agent_safe`` (see :mod:`.agent_manifest`). New
or half-baked commands are therefore captured-but-invisible by default -- the
opt-in allowlist posture (design decision 3).
'''
from dataclasses import dataclass, field
from typing import Optional, Callable, Any
import contextlib


@dataclass
class CommandRecord:
    '''One captured ISOLDE command and its agent-surface metadata.'''
    name: str                       # full registered name, e.g. 'isolde sim'
    cmd_desc: Any                   # the CmdDesc instance (typed arg descriptors)
    function: Callable              # the registered callable, func(session, **kw)
    synopsis: str = ''              # one-line description (cmd_desc.synopsis)
    category: Optional[str] = None  # grouping for the manifest (e.g. 'simulation')
    agent_safe: bool = False        # OPT-IN: only True commands are agent-reachable
    witness: Optional[str] = None   # post-condition witness id (see .witness)
    long_running: bool = False      # uses the async job pattern (see jobs.py)


# Process-lifetime singleton. Populated by capture_registrations() the first (and
# only) time register_isolde runs; survives for the life of the ChimeraX session.
_ISOLDE_COMMAND_REGISTRY: "dict[str, CommandRecord]" = {}


def get_registry() -> "dict[str, CommandRecord]":
    '''Return the live registry mapping full command name -> CommandRecord.'''
    return _ISOLDE_COMMAND_REGISTRY


def get_command(name: str) -> Optional[CommandRecord]:
    '''Return the record for *name*, or None if not registered.'''
    return _ISOLDE_COMMAND_REGISTRY.get(name)


def agent_safe_commands() -> "dict[str, CommandRecord]":
    '''Return only the commands explicitly flagged agent_safe.'''
    return {n: r for n, r in _ISOLDE_COMMAND_REGISTRY.items() if r.agent_safe}


@contextlib.contextmanager
def capture_registrations():
    '''
    Context manager that records every ``register(...)`` call made while active
    into :data:`_ISOLDE_COMMAND_REGISTRY`, then restores the original
    ``chimerax.core.commands.register``.

    Exception-safe: the original is always restored in the ``finally`` so a
    failure mid-registration never leaves ChimeraX's ``register`` permanently
    monkeypatched. Idempotent with respect to the registry (re-running simply
    overwrites records by name).
    '''
    import chimerax.core.commands as cc
    original_register = cc.register

    def _recording_register(name, cmd_desc=(), function=None, *args, **kwargs):
        # Record only genuine command registrations: a real CmdDesc plus a
        # callable. ChimeraX also allows register() as a decorator (function
        # omitted) and a deferred (cmd_desc, function) tuple form; we let those
        # through untouched rather than recording a half-specified entry.
        try:
            if function is not None and _looks_like_cmd_desc(cmd_desc):
                _record(name, cmd_desc, function)
        except Exception as e:
            # Never let bookkeeping break command registration.
            _warn('failed to record command %r: %s' % (name, e))
        return original_register(name, cmd_desc, function, *args, **kwargs)

    cc.register = _recording_register
    try:
        yield _ISOLDE_COMMAND_REGISTRY
    finally:
        cc.register = original_register


def _record(name, cmd_desc, function):
    '''Build and store a CommandRecord, merging in opt-in agent metadata.'''
    from .agent_manifest import metadata_for
    meta = metadata_for(name)
    _ISOLDE_COMMAND_REGISTRY[name] = CommandRecord(
        name=name,
        cmd_desc=cmd_desc,
        function=function,
        synopsis=(getattr(cmd_desc, 'synopsis', '') or ''),
        category=meta.get('category'),
        agent_safe=meta.get('agent_safe', False),
        witness=meta.get('witness'),
        long_running=meta.get('long_running', False),
    )


def _looks_like_cmd_desc(obj):
    '''True if *obj* is a CmdDesc (has the typed-arg dicts we introspect).'''
    return all(hasattr(obj, attr) for attr in ('_required', '_optional', '_keyword'))


def _warn(msg):
    try:
        from chimerax.core.session import Session  # noqa: F401
        import sys
        print('[isolde.command_registry] %s' % msg, file=sys.stderr)
    except Exception:
        pass
