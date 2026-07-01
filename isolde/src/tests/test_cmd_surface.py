# @Author: Tristan Croll
# @Date:   30-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 30-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Characterisation test for the stable edges of ISOLDE's command surface: the
read-only ``isolde status`` report and the custom ``IsoldeStructureArg``
argument parser. Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_cmd_surface.py

Prints PASS/FAIL and exits non-zero on failure.
'''
import os

FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1pmx_1.pdb')


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def run(session):
    from chimerax.core.commands import run as run_cmd

    # --- isolde status: read-only, works before "isolde start" ------------
    from chimerax.isolde.cmd.cmd import isolde_status
    status = isolde_status(session)
    for key in ('isolde_started', 'selected_model', 'simulation_running', 'forcefield'):
        if key not in status:
            _fail('isolde status lost key %r' % key)
    if not isinstance(status['isolde_started'], bool):
        _fail('isolde_started should be bool, got %r' % (status['isolde_started'],))
    if not isinstance(status['simulation_running'], bool):
        _fail('simulation_running should be bool, got %r' % (status['simulation_running'],))
    print('PASS: isolde status returns the documented dict shape')

    # --- IsoldeStructureArg: single-structure parse + rejection -----------
    m = run_cmd(session, 'open "%s"' % FIXTURE.replace('\\', '/'), log=False)[0]
    try:
        from chimerax.isolde.cmd.argspec import IsoldeStructureArg
        from chimerax.core.commands import AnnotationError
        parsed, text, rest = IsoldeStructureArg.parse('#%s' % m.id_string, session)
        if parsed is not m:
            _fail('IsoldeStructureArg did not parse to the opened model')
        print('PASS: IsoldeStructureArg parses a single structure spec')

        rejected = False
        try:
            IsoldeStructureArg.parse('#999', session)
        except AnnotationError:
            rejected = True
        if not rejected:
            _fail('IsoldeStructureArg accepted a spec matching no structure')
        print('PASS: IsoldeStructureArg rejects a no-match spec (AnnotationError)')
    finally:
        session.models.close([m])

    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
