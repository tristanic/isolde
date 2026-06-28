# @Author: Tristan Croll
# @Date:   27-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 27-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Smoke test for the Phase-A agent control surface (command registry + schema +
invoke). Run inside ChimeraX, e.g.

    run_chimerax.bat --nogui --exit --script src/tests/test_agent_surface.py

It triggers ISOLDE command registration, then asserts:
  - the registry auto-discovered the commands;
  - exposure is opt-in (agent_safe is a strict subset of all commands);
  - every agent-safe command produces a well-formed JSON-Schema tool definition;
  - a representative command round-trips through build_command_text.

Prints PASS/FAIL lines and exits non-zero on failure so it can gate a build.
'''
import sys


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def run(session):
    from chimerax.core.commands import run as run_cmd

    # Trigger registration of all ISOLDE commands (lazy: any isolde command does
    # it). 'isolde status' is side-effect-free when ISOLDE isn't started.
    status = run_cmd(session, 'isolde status', log=False)
    print('isolde status ->', status)

    from chimerax.isolde.cmd.command_registry import (
        get_registry, agent_safe_commands)
    from chimerax.isolde.cmd.schema import (
        command_tool_definition, command_input_schema)
    from chimerax.isolde.cmd.invoke import build_command_text

    reg = get_registry()
    safe = agent_safe_commands()
    print('PASS: registry captured %d commands' % len(reg))
    if len(reg) < 15:
        _fail('expected many commands, got %d' % len(reg))

    # Opt-in: agent-safe is a strict subset; unflagged commands stay invisible.
    if not safe:
        _fail('no agent_safe commands -- manifest not applied')
    if len(safe) >= len(reg):
        _fail('every command is agent_safe -- opt-in posture broken')
    print('PASS: %d/%d commands exposed (opt-in subset)' % (len(safe), len(reg)))

    # A known command that should NOT be exposed (no manifest entry).
    if 'isolde tutorial' in safe:
        _fail('isolde tutorial unexpectedly agent_safe')
    print('PASS: unflagged command (isolde tutorial) correctly hidden')

    # Every agent-safe command yields a well-formed tool definition, with an
    # MCP-valid tool name (SEP-986: only [A-Za-z0-9_.-]) that still carries the
    # canonical command verbatim so the MCP server can invoke without un-munging.
    import re
    valid_name = re.compile(r'^[A-Za-z0-9_.-]+$')

    def _has_key(obj, key):
        '''Recursively test whether `key` appears anywhere in a JSON-Schema dict.'''
        if isinstance(obj, dict):
            return key in obj or any(_has_key(v, key) for v in obj.values())
        if isinstance(obj, (list, tuple)):
            return any(_has_key(v, key) for v in obj)
        return False

    n_schemas = 0
    for name, rec in safe.items():
        tool = command_tool_definition(rec)
        s = tool['input_schema']
        if s.get('type') != 'object' or 'properties' not in s:
            _fail('bad schema for %s: %r' % (name, s))
        if not valid_name.match(tool['name']):
            _fail('tool name %r (from command %r) is not MCP-valid'
                  % (tool['name'], name))
        if tool.get('command') != name:
            _fail('tool for %r lost its canonical command field: %r'
                  % (name, tool.get('command')))
        # We standardise Or -> anyOf; oneOf ("exactly one") makes strict MCP
        # clients reject values that satisfy >1 branch (e.g. 'isolde stepto first').
        if _has_key(s, 'oneOf'):
            _fail('schema for %s uses oneOf (should be anyOf): %r' % (name, s))
        n_schemas += 1
    print('PASS: %d tool schemas, all names MCP-valid + carry canonical command'
          % n_schemas)

    # Specifically exercise the '~' negation case if present (the bug this guards).
    from chimerax.isolde.cmd.schema import tool_name_for_command
    if not valid_name.match(tool_name_for_command('isolde ~ignore')):
        _fail('tool_name_for_command did not sanitise a ~ command')
    if tool_name_for_command('isolde ~ignore') == tool_name_for_command('isolde ignore'):
        _fail("'~ignore' and 'ignore' collide after sanitising")
    print('PASS: ~negation commands sanitise to distinct MCP-valid names')

    # Show a few representative tool definitions for eyeballing.
    for name in ('isolde sim', 'isolde validate rama', 'isolde stepto',
                 'isolde adjust bfactors'):
        rec = reg.get(name)
        if rec is None:
            print('  (note: %r not registered)' % name)
            continue
        import json
        print('--- %s ---' % name)
        print(json.dumps(command_tool_definition(rec), indent=2, default=str))

    # Round-trip a command through the text serializer (no execution).
    rec = reg.get('isolde stepto')
    if rec is not None:
        text = build_command_text(rec, {'view_distance': 12.0, 'select': True})
        print('build_command_text(isolde stepto) ->', repr(text))
        if 'view_distance 12' not in text or 'select true' not in text:
            _fail('unexpected serialization: %r' % text)
        print('PASS: command text serialization round-trips')

    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
