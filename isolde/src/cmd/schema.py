# @Author: Tristan Croll
# @Date:   27-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 27-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Translate ChimeraX command type descriptors (CmdDesc / Annotation) into
JSON Schema, so each agent-safe ISOLDE command can be advertised as a typed
MCP/agent tool with a machine-checkable input schema.

The mapping is deliberately conservative: anything we don't recognise degrades
to a permissive ``{"type": "string"}`` with a warning, rather than raising --
an unknown annotation should make a command *loosely* typed, never un-buildable.

Atom/residue/structure/map arguments all map to a plain JSON *string* carrying a
ChimeraX atom-spec (e.g. ``#1/A:10-20@CA``). The agent supplies the spec text and
ChimeraX's own parser resolves it at invocation time (see :mod:`.invoke`), so we
never have to evaluate atom-specs ourselves nor round-trip live objects.
'''
from typing import Any


# Scalar annotation classes compared by identity / subclass. Imported lazily in
# _scalar_schema because the exact import path is only valid once ChimeraX core
# is loaded.
def _atomspec_schema(name_hint=''):
    desc = 'A ChimeraX atom specifier string'
    if name_hint:
        desc += ' (%s)' % name_hint
    desc += ", e.g. '#1', '#1/A:10-20', '#1/A:58@CA', 'sel'."
    return {'type': 'string', 'description': desc}


def annotation_to_json_schema(annotation) -> dict:
    '''
    Map a single ChimeraX command Annotation (class or instance) to a JSON
    Schema fragment. Never raises: unknown annotations fall back to a string.
    '''
    from chimerax.core.commands import cli
    # --- composite annotations (instances carrying inner type info) ---
    # Or(A, B, ...) -> anyOf. ChimeraX's Or tries each annotation in turn and
    # accepts the first that parses, so "matches at least one" (anyOf) is correct.
    # oneOf ("matches exactly one") is wrong here: values like the literal "first"
    # satisfy more than one branch (e.g. an EnumOf *and* a string), which a strict
    # JSON-Schema validator (some MCP clients) then rejects.
    if isinstance(annotation, cli.Or):
        subs = [annotation_to_json_schema(a) for a in annotation.annotations]
        return {'anyOf': subs, 'description': _name_of(annotation)}
    # EnumOf(values) / DynamicEnum -> string enum
    if isinstance(annotation, cli.EnumOf):
        return {'type': 'string', 'enum': [str(v) for v in annotation.values]}
    if isinstance(annotation, getattr(cli, 'DynamicEnum', ())):
        try:
            return {'type': 'string', 'enum': [str(v) for v in annotation.values_func()]}
        except Exception:
            return {'type': 'string'}
    # Bounded(anno, min, max) -> inner schema + numeric bounds
    if isinstance(annotation, cli.Bounded):
        inner = annotation_to_json_schema(annotation.anno)
        if annotation.min is not None:
            inner['minimum' if annotation.inclusive else 'exclusiveMinimum'] = annotation.min
        if annotation.max is not None:
            inner['maximum' if annotation.inclusive else 'exclusiveMaximum'] = annotation.max
        return inner
    # ListOf / SetOf (Aggregate) -> array of inner schema
    if isinstance(annotation, cli.Aggregate):
        items = annotation_to_json_schema(annotation.annotation)
        out = {'type': 'array', 'items': items}
        if getattr(annotation, 'min_size', 0):
            out['minItems'] = annotation.min_size
        return out
    # RepeatOf(anno) -> array. NOTE: RepeatOf does not retain its inner
    # annotation as an attribute in this ChimeraX (it only copies .parse), so we
    # cannot recover the item type; advertise a generic array.
    if isinstance(annotation, cli.RepeatOf):
        return {'type': 'array', 'items': {}}

    # --- scalar annotations (usually the class itself, sometimes an instance) ---
    return _scalar_schema(annotation)


def _scalar_schema(annotation) -> dict:
    from chimerax.core.commands import cli
    cls = annotation if isinstance(annotation, type) else type(annotation)

    # Atom-spec family: one subclass check covers AtomsArg, ResiduesArg,
    # AtomicStructuresArg, StructuresArg, StructureArg, ModelArg, MapArg,
    # AtomSpecArg, and ISOLDE's custom IsoldeStructureArg / map handler args.
    atomspec_bases = []
    for nm in ('AtomSpecArg', 'ModelArg', 'ModelsArg', 'ObjectsArg'):
        b = getattr(cli, nm, None)
        if isinstance(b, type):
            atomspec_bases.append(b)
    try:
        from chimerax import atomic
        for nm in ('AtomsArg', 'ResiduesArg', 'AtomicStructureArg',
                   'AtomicStructuresArg', 'StructureArg', 'StructuresArg'):
            b = getattr(atomic, nm, None)
            if isinstance(b, type):
                atomspec_bases.append(b)
    except Exception:
        pass
    if atomspec_bases and issubclass(cls, tuple(atomspec_bases)):
        return _atomspec_schema(_name_of(annotation))

    # Booleans / flags
    if cls is getattr(cli, 'BoolArg', None):
        return {'type': 'boolean'}
    if cls is getattr(cli, 'NoArg', None):
        return {'type': 'boolean', 'description': 'flag (presence sets it true)'}

    # Integers
    int_args = [getattr(cli, n, None) for n in ('IntArg', 'PositiveIntArg', 'NonNegativeIntArg', 'Bounded')]
    if cls is getattr(cli, 'PositiveIntArg', None):
        return {'type': 'integer', 'minimum': 1}
    if cls is getattr(cli, 'NonNegativeIntArg', None):
        return {'type': 'integer', 'minimum': 0}
    if cls is getattr(cli, 'IntArg', None):
        return {'type': 'integer'}

    # Floats
    if cls in (getattr(cli, 'FloatArg', None), getattr(cli, 'NonNegativeFloatArg', None),
               getattr(cli, 'PositiveFloatArg', None)):
        s = {'type': 'number'}
        if cls is getattr(cli, 'NonNegativeFloatArg', None):
            s['minimum'] = 0
        if cls is getattr(cli, 'PositiveFloatArg', None):
            s['exclusiveMinimum'] = 0
        return s

    # Strings / filenames
    if cls in (getattr(cli, 'StringArg', None), getattr(cli, 'FileNameArg', None),
               getattr(cli, 'OpenFileNameArg', None), getattr(cli, 'SaveFileNameArg', None),
               getattr(cli, 'OpenFolderNameArg', None)):
        s = {'type': 'string'}
        if 'File' in cls.__name__ or 'Folder' in cls.__name__:
            s['description'] = 'A filesystem path'
        return s

    # Unknown: degrade to a string, but keep the human description if available.
    name = _name_of(annotation)
    out = {'type': 'string'}
    if name:
        out['description'] = name
    _warn('no JSON-schema mapping for annotation %r (%s); using string'
          % (getattr(cls, '__name__', cls), name))
    return out


def _name_of(annotation) -> str:
    n = getattr(annotation, 'name', None)
    return n if isinstance(n, str) else ''


def command_input_schema(record) -> dict:
    '''
    Build the JSON-Schema object for a CommandRecord's arguments.

    required[]  -> properties + JSON-Schema "required"
    optional[]  -> properties (positional-but-optional; ChimeraX also copies
                   these into _keyword, so we dedupe)
    keyword[]   -> properties
    '''
    desc = record.cmd_desc
    properties: dict = {}
    required: list = []

    def _add(pairs, mark_required=False):
        if not pairs:
            return
        for name, annotation in dict(pairs).items():
            if name in properties:
                continue  # dedupe optional/keyword overlap; first wins
            try:
                properties[name] = annotation_to_json_schema(annotation)
            except Exception as e:
                properties[name] = {'type': 'string'}
                _warn('schema for arg %r failed: %s' % (name, e))
            if mark_required:
                required.append(name)

    _add(getattr(desc, '_required', {}), mark_required=True)
    _add(getattr(desc, '_optional', {}))
    _add(getattr(desc, '_keyword', {}))

    schema = {'type': 'object', 'properties': properties, 'additionalProperties': False}
    if required:
        schema['required'] = required
    return schema


import re

# MCP tool names allow only these characters (SEP-986); spaces and ChimeraX's '~'
# negation prefix are not among them.
_INVALID_TOOL_NAME_CHAR = re.compile(r'[^A-Za-z0-9_.-]')


def tool_name_for_command(command: str) -> str:
    '''Map a ChimeraX command name to an MCP-valid tool name.

    MCP tool names permit only ``[A-Za-z0-9_.-]``. ChimeraX uses a leading ``~``
    for the negated form of a command (e.g. ``isolde ~ignore`` undoes ``isolde
    ignore``); we render it as a distinct ``not`` token so the two never collide,
    collapse whitespace to underscores, then replace any other stray character.
    The reverse mapping is NOT done by string surgery -- consumers invoke via the
    canonical name carried in the tool definition's ``command`` field.
    '''
    name = command.replace('~', ' not ')
    name = '_'.join(name.split())
    return _INVALID_TOOL_NAME_CHAR.sub('_', name)


def command_tool_definition(record) -> dict:
    '''
    Produce the agent/MCP tool definition for a CommandRecord:
        {name, command, description, category, long_running, input_schema}
    The tool ``name`` is an MCP-valid label (see :func:`tool_name_for_command`);
    the canonical ChimeraX command name is preserved verbatim in ``command`` (and
    the schema's ``x-command``) -- always invoke via that, never by un-munging the
    name.
    '''
    description = record.synopsis
    if not description and record.function is not None:
        description = (record.function.__doc__ or '').strip().split('\n')[0]
    schema = command_input_schema(record)
    schema['x-command'] = record.name
    return {
        'name': tool_name_for_command(record.name),
        'command': record.name,
        'description': description or record.name,
        'category': record.category,
        'long_running': record.long_running,
        'input_schema': schema,
    }


def _warn(msg):
    import sys
    print('[isolde.schema] %s' % msg, file=sys.stderr)
