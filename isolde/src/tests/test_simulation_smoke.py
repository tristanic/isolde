# @Author: Tristan Croll
# @Date:   30-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 30-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
End-to-end smoke test: ISOLDE can build a simulation System, run a few
minimisation/equilibration rounds with restraints attached, and stop cleanly.
Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_simulation_smoke.py

Asserts the durable behavioural contract -- "model in -> running, stoppable
simulation out" -- that the future ML-forcefield rewrite must preserve. It does
NOT pin energies/coordinates (Langevin dynamics is stochastic); it only checks
the machine turns over. Headless has no render loop, so it pumps the engine's
'new frame' trigger by hand. Forces the CPU platform for determinism. Prints
PASS/FAIL and exits non-zero on failure.
'''
import os
import time
import numpy

FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1pmx_1.pdb')

N_MOBILE_RESIDUES = 8
N_FIXED_RESIDUES = 6


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def _pump(session, n, stop_when=None, dt=0.02):
    for _ in range(n):
        session.triggers.activate_trigger('new frame', None)
        time.sleep(dt)
        if stop_when is not None and stop_when():
            return True
    return False


def run(session):
    from chimerax.core.commands import run as run_cmd
    from chimerax.isolde import session_extensions as sx
    from chimerax.isolde.openmm import openmm_interface, sim_param_mgr
    from chimerax.isolde.openmm.forcefields import ForcefieldMgr

    m = run_cmd(session, 'open "%s"' % FIXTURE.replace('\\', '/'), log=False)[0]
    handler = None
    try:
        residues = m.residues
        ordered = residues[numpy.lexsort((residues.numbers, residues.chain_ids))]
        mobile = ordered[:N_MOBILE_RESIDUES].atoms
        fixed = ordered[N_MOBILE_RESIDUES:N_MOBILE_RESIDUES + N_FIXED_RESIDUES].atoms

        construct = openmm_interface.SimConstruct(m, mobile, fixed)
        params = sim_param_mgr.SimParams()
        try:
            params.set_param('platform', 'CPU')
        except Exception:
            params.platform = 'CPU'
        from chimerax.isolde.openmm.param_provider import ForceFieldParameterisationProvider
        handler = openmm_interface.SimHandler(session, params, construct,
                                              ForceFieldParameterisationProvider(ForcefieldMgr(session)))

        # Exercise the restraint-force wiring too: position-restrain the mobile
        # heavy atoms (a core feature the rewrite must keep working).
        handler.initialize_restraint_forces()
        pr_mgr = sx.get_position_restraint_mgr(m)
        prs = pr_mgr.add_restraints(construct.mobile_atoms)
        prs.targets = prs.atoms.coords
        prs.spring_constants = 1000.0
        prs.enableds = True
        handler.add_position_restraints(prs)

        updates = [0]
        handler.triggers.add_handler(
            'coord update', lambda *a: updates.__setitem__(0, updates[0] + 1))
        start_coords = m.atoms.coords.copy()

        handler.start_sim()
        _pump(session, 800, stop_when=lambda: updates[0] >= 5)
        if not handler.sim_running:
            _fail('simulation did not start / stay running')
        if updates[0] < 5:
            _fail('sim loop did not advance (only %d coord updates)' % updates[0])
        moved = float(numpy.abs(m.atoms.coords - start_coords).max())
        if moved <= 1e-3:
            _fail('no atom moved -- dynamics not actually running')
        print('PASS: simulation started and advanced (%d coord updates, max move %.3f A)'
              % (updates[0], moved))

        handler.stop()
        _pump(session, 300, stop_when=lambda: not handler.sim_running)
        if handler.sim_running:
            _fail('simulation did not stop cleanly')
        print('PASS: simulation stopped cleanly')
    finally:
        try:
            if handler is not None and handler.sim_running:
                handler.stop()
                _pump(session, 200, stop_when=lambda: not handler.sim_running)
        except Exception:
            pass
        session.models.close([m])

    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
