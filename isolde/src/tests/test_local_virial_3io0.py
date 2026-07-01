# @Author: Tristan Croll
# @Date:   01-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 01-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Acceptance test (go/no-go) for the local-virial baseline-subtraction signal,
reproducing section 6 of BRIEF_for_ISOLDE_agent.md. Manual/interactive, like
SimTester in test_simulation.py -- ISOLDE has no pytest/CI (see the project's
CLAUDE.md), so this prints rankings for visual inspection rather than asserting
on hardcoded atom indices (which are specific to the original prototype's own
atom ordering/analysis and aren't safe to assume still line up exactly here).

Run inside ChimeraX, with ISOLDE already available as session.isolde:

    from chimerax.isolde.tests.test_local_virial_3io0 import run_test
    before = run_test(session, before=True)
    after = run_test(session, before=False)

The pass criterion (brief section 6): ranking by EXCESS von Mises, the
persistent Arg-CZ entries that dominate the RAW ranking in both `before` and
`after` should drop substantially in the EXCESS ranking, while the cis-peptide
Calpha error (near the top of `before`'s RAW ranking) should stay near the top
of `before`'s EXCESS ranking and disappear from `after`'s. Compare the two
printed rankings by eye.
'''

import os

import numpy as np

from ..openmm import openmm_interface, sim_param_mgr
from ..validation.local_virial import LocalVirialCalculator
from ..validation.reference_stress import ReferenceStressLibrary

_DEMO_DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'demo_data', '3io0')


def run_test(session, before=True, top_n=15):
    '''
    Load ISOLDE's bundled 3io0 crystallographic demo model (before.pdb -- the
    un-fixed model containing the cis-peptide/flip/rotamer errors described in
    the brief -- or after.pdb, the fixed version), minimize once with ISOLDE's
    own default settings without touching the model further, and print atoms
    ranked by both raw and excess (baseline-subtracted) von Mises stress.

    Returns a dict with the atoms and both scalar-measure dicts, for further
    inspection/comparison in the Python shell.
    '''
    from chimerax.open_command.cmd import provider_open
    from chimerax.atomic import Atoms

    fname = 'before.pdb' if before else 'after.pdb'
    model = provider_open(session, [os.path.join(_DEMO_DATA_DIR, fname)])[0]

    all_atoms = model.atoms
    sim_construct = openmm_interface.SimConstruct(model, all_atoms, Atoms(), None)
    sim_params = sim_param_mgr.SimParams()
    sim_params.platform = 'CPU'
    forcefield_mgr = session.isolde.forcefield_mgr
    sim_handler = openmm_interface.SimHandler(session, sim_params, sim_construct, forcefield_mgr)
    sim_handler.initialize_restraint_forces(
        amber_cmap=False, tugging=False, position_restraints=False,
        distance_restraints=False, adaptive_distance_restraints=False,
        dihedral_restraints=False, adaptive_dihedral_restraints=False)
    sim_handler._prepare_sim()
    th = sim_handler._thread_handler
    th.minimize()
    th.finalize_thread()

    state = sim_handler._context.getState(getPositions=True)
    calc = LocalVirialCalculator(sim_handler._system)
    virial = calc.compute(state.getPositions(asNumpy=True))['virial']

    raw = LocalVirialCalculator.scalar_measures(virial)
    lib = ReferenceStressLibrary()
    excess = lib.compute_excess(virial, sim_construct.all_atoms)

    label = 'before' if before else 'after'
    _print_ranking(sim_construct.all_atoms, raw['von_mises'], f'{label}: RAW von Mises', top_n)
    _print_ranking(
        sim_construct.all_atoms, excess['von_mises'], f'{label}: EXCESS von Mises', top_n)
    return dict(atoms=sim_construct.all_atoms, raw=raw, excess=excess)


def _print_ranking(atoms, values, label, top_n):
    order = np.argsort(np.nan_to_num(values, nan=-np.inf))[::-1]
    print(f'\n=== Top {top_n} by {label} ===')
    for rank, i in enumerate(order[:top_n]):
        a = atoms[i]
        r = a.residue
        print(f'  {rank + 1:2d}. idx={i:5d}  {r.name}{r.number}/{r.chain_id}  {a.name:4s}  '
              f'von_mises={values[i]:+8.2f}')
