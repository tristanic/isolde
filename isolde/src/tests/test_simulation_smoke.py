# @Author: Tristan Croll
# @Date:   19-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 19-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
End-to-end *smoke* test: can ISOLDE still build a simulation System, run a few
minimisation/equilibration rounds with restraints attached, and stop cleanly?

This asserts the durable behavioural contract -- "model in -> running, stoppable
simulation out" -- that the future ML-forcefield rewrite of simulation startup
must preserve. It deliberately does NOT pin energies or coordinates (Langevin
dynamics is stochastic); it only checks the machine turns over.

It drives the engine through the manager + SimHandler layer rather than the
``isolde sim start`` command, because that command is hard-gated to GUI mode
(cmd.py: "ISOLDE currently requires ChimeraX to be in GUI mode") and cannot run
under ``--nogui``. Headless, there is no render loop, so we pump the engine's
``'new frame'`` trigger by hand to turn its simulation loop.
'''

import os
import time

import numpy
import pytest

pytestmark = pytest.mark.integration

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE_PDB = os.path.join(HERE, '1pmx_1.pdb')


def _pump(session, n, stop_when=None, dt=0.02):
    '''Fire 'new frame' up to n times, sleeping between to let the worker
    thread finish each chunk. Returns early once stop_when() is true.'''
    for _ in range(n):
        session.triggers.activate_trigger('new frame', None)
        time.sleep(dt)
        if stop_when is not None and stop_when():
            return True
    return False


def _force_cpu(params):
    # This machine (and CI) may default to a GPU platform that isn't registered
    # here; CPU is always available and keeps the smoke test deterministic.
    try:
        params.set_param('platform', 'CPU')
    except Exception:
        params.platform = 'CPU'


def test_simulation_runs_end_to_end(session):
    from chimerax.core.commands import run
    from chimerax.geometry import find_close_points
    from chimerax.isolde import session_extensions as sx
    from chimerax.isolde.openmm import openmm_interface, sim_param_mgr
    from chimerax.isolde.openmm.forcefields import ForcefieldMgr

    m = run(session, 'open "{}"'.format(FIXTURE_PDB.replace('\\', '/')))[0]
    handler = None
    try:
        # Mobile: a handful of residues. Fixed: an ~8 A shell of whole residues.
        mobile = m.residues[:8].atoms
        others = m.atoms.subtract(mobile)
        _, shell_idx = find_close_points(mobile.coords, others.coords, 8.0)
        fixed = others[shell_idx].unique_residues.atoms

        construct = openmm_interface.SimConstruct(m, mobile, fixed)
        params = sim_param_mgr.SimParams()
        _force_cpu(params)
        handler = openmm_interface.SimHandler(
            session, params, construct, ForcefieldMgr(session))

        # Exercise the restraint-force wiring too: position-restrain the mobile
        # atoms (a core ISOLDE feature the rewrite must keep working).
        handler.initialize_restraint_forces()
        pr_mgr = sx.get_position_restraint_mgr(m)
        # add_restraints restrains only the heavy atoms, so drive targets from
        # the restraints' own atoms rather than the full mobile selection.
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

        assert handler.sim_running, 'simulation did not start / stay running'
        assert updates[0] >= 5, \
            'sim loop did not advance (only {} coord updates)'.format(updates[0])
        moved = float(numpy.abs(m.atoms.coords - start_coords).max())
        assert moved > 1e-3, 'no atom moved -- dynamics not actually running'

        handler.stop()
        _pump(session, 300, stop_when=lambda: not handler.sim_running)
        assert not handler.sim_running, 'simulation did not stop cleanly'
    finally:
        try:
            if handler is not None and handler.sim_running:
                handler.stop()
                _pump(session, 200, stop_when=lambda: not handler.sim_running)
        except Exception:
            pass
        session.models.close([m])
