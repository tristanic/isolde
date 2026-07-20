# Manual GUI harness for per-group soft-core nonbonded coupling.
#
# This is NOT headless-runnable: it needs a live ISOLDE simulation (a real GL
# context; `isolde start` refuses --nogui). The per-group *force math* is covered
# headlessly by test_nb_group_softcore.py; this harness exercises the live
# build-path + SimHandler API (Phases 2-3) in the GUI.
#
# ------------------------------------------------------------------------------
# How to run (in a GUI ChimeraX with ISOLDE installed):
#
#   1. Open a model and start ISOLDE on it:
#          open 3io0                      # or any small model with a ligand/loop
#          isolde start
#
#   2. Enable per-group coupling BEFORE starting a simulation (the coupling table
#      is sized when the soft-core forces are built, during sim start):
#          session.isolde.sim_params.nb_groups_max = 8
#
#   3. Select the atoms you want to explore (e.g. one sidechain or a ligand) and
#      start a simulation covering them:
#          select /A:60 & sidechain
#          isolde sim start sel
#
#   4. Pause the simulation (spacebar, or `isolde sim pause`), then in the
#      ChimeraX Python shell:
#          from chimerax.isolde.tests.test_nb_groups_gui import check
#          check(session)                 # uses the current selection as the group
#
# `check` puts the current selection into group 1, softens its coupling to the
# rest of the model (group 0), reports the nonbonded + GB force-group energies
# before/after, then ramps it back to full — all on the live context with no
# reinitialisation. Watch the ISOLDE log for any "reinitializing context" message
# (there should be NONE), and confirm the simulation stays stable throughout.
# ------------------------------------------------------------------------------

import numpy


def _group_energies(sh):
    '''Potential energy of just the nonbonded force group, in kJ/mol.'''
    from chimerax.isolde.openmm.openmm_interface import NONBONDED_FORCE_GROUP
    c = sh._context
    state = c.getState(getEnergy=True, groups={NONBONDED_FORCE_GROUP})
    from openmm import unit
    return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


def check(session, atoms=None):
    isolde = session.isolde
    sh = isolde.sim_handler
    if sh is None or not sh.sim_running:
        print('nb-group check: no running simulation. Start one first.')
        return
    if not sh.nb_groups_enabled:
        print('nb-group check: per-group coupling is NOT enabled for this sim.\n'
              '  Set `session.isolde.sim_params.nb_groups_max = 8` BEFORE '
              '`isolde sim start`, then restart the simulation.')
        return
    if atoms is None:
        from chimerax.atomic import selected_atoms
        atoms = selected_atoms(session)
    if len(atoms) == 0:
        print('nb-group check: no atoms selected/passed.')
        return

    print(f'nb-group check: {len(atoms)} atoms -> group 1')
    e0 = _group_energies(sh)
    sh.assign_nb_group(atoms, 1)
    e1 = _group_energies(sh)
    print(f'  nonbonded E: baseline {e0:.3f}  after grouping (coupling still 1) {e1:.3f}')
    print('  (should be ~equal: grouping alone, coupling unchanged, changes nothing)')

    for lam in (0.6, 0.3, 0.1):
        sh.set_nb_coupling(1, 0, lam)
        print(f'  set_nb_coupling(1,0,{lam})  ->  nonbonded E = {_group_energies(sh):.3f} '
              f'(get_nb_coupling(1,0)={sh.get_nb_coupling(1,0)})')

    print('  ramping coupling(1,0) back to 1.0 ...')
    for lam in (0.3, 0.6, 1.0):
        sh.set_nb_coupling(1, 0, lam)
    e_end = _group_energies(sh)
    print(f'  nonbonded E after full re-coupling: {e_end:.3f} (baseline was {e1:.3f})')
    print('  PASS if: no "reinitializing context" in the log, sim stayed stable,\n'
          '  softening reduced clash energy, and re-coupling returned near baseline.')
