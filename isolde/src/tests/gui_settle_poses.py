# Manual GUI harness for settle_poses (NOT headless-runnable -- needs a live ISOLDE
# simulation with a real GL context; `isolde start` refuses --nogui).
#
# The pure logic is covered by test_settle_poses.py and the E_pack force-math invariant by
# test_nb_group_softcore.py. This harness exercises the full engine in the GUI: it feeds
# synthetic poses for a chosen residue, settles + scores + commits them, and checks that the
# live edits leave the simulation consistent.
#
# ------------------------------------------------------------------------------
# How to run (in a GUI ChimeraX with ISOLDE + Clipper installed):
#
#   open 3io0                         # or any small model with a ligand/sidechain
#   isolde start
#   select /A:60 & sidechain          # a residue (or ligand) to fit; or select a ligand
#   isolde sim start sel
#   isolde sim pause                  # (or spacebar) -- FAST mode needs a paused sim
#
#   from chimerax.isolde.tests.gui_settle_poses import check, ab, check_live
#   check(session)                    # fast burst on the selected residue
#   ab(session)                       # localisation A/B: E_pack vs whole-construct CORE
#   check_live(session)               # live (watchable) driver -- sim may be running for this
#
# Watch the ISOLDE log for any "reinitializing context" message (there should be NONE -- all
# edits are live). `check` reports a PASS/FAIL summary of the programmatic assertions.
# ------------------------------------------------------------------------------

import numpy as np


def _selected_residue(session, residue):
    if residue is not None:
        return residue
    from chimerax.atomic import selected_residues, selected_atoms
    rs = selected_residues(session)
    if len(rs) == 1:
        return rs[0]
    ats = selected_atoms(session)
    if len(ats):
        rr = ats.unique_residues
        if len(rr) == 1:
            return rr[0]
    print('gui_settle_poses: select exactly one residue (or pass residue=...).')
    return None


def _synthetic_poses(residue):
    '''A few translated poses (in residue.atoms order) to exercise the machinery: a small
    nudge (should settle back / be kept), a medium shift, and a big clash (should be culled
    or rank worst). The current conformation is added by settle_poses itself.'''
    base = residue.atoms.coords
    def shift(dx):
        c = base.copy()
        c[:, 0] += dx
        return c
    return [('nudge_0.4A', shift(0.4)),
            ('shift_1.2A', shift(1.2)),
            ('clash_2.5A', shift(2.5))]


def _env_coords(session, residue):
    '''Coordinates of every simulated atom EXCEPT this residue's -- the environment that a
    correct graft-commit must leave untouched.'''
    sh = session.isolde.sim_handler
    atoms = sh.atoms
    ridx = atoms.indices(residue.atoms)
    mask = np.ones(len(atoms), dtype=bool)
    mask[ridx[ridx >= 0]] = False
    return atoms.coords[mask], mask


def _mdff_enableds_all_true(sh, residue):
    '''True if every MDFF map has this residue's atoms enabled (i.e. map_decoupled restored
    them). True trivially when there are no maps.'''
    for mgr in sh._mdff_managers():
        m = mgr.get_mdff_atoms(residue.atoms)
        if m is not None and len(m) and not bool(np.all(m.enableds)):
            return False
    return True


def check(session, residue=None, map_decouple=True, debug=True):
    '''Fast-burst functional test on a paused sim. Asserts: the environment is untouched by
    the commit, the sim is still running, temperature is restored, and the ligand's MDFF
    terms are re-enabled. Reports the SettleResult and a PASS/FAIL summary.'''
    from chimerax.isolde.refine import settle_poses
    residue = _selected_residue(session, residue)
    if residue is None:
        return
    isolde = session.isolde
    sh = isolde.sim_handler
    if sh is None or not sh.sim_running:
        print('gui_settle_poses.check: no running simulation. "isolde sim start sel" first.')
        return
    if not sh.pause:
        print('gui_settle_poses.check: pause the sim first ("isolde sim pause"); fast mode '
              'needs a paused sim. (Or use check_live for the running driver.)')
        return

    env_before, env_mask = _env_coords(session, residue)
    configured_t = isolde.sim_params.temperature
    configured_t = (configured_t.value_in_unit(configured_t.unit)
                    if hasattr(configured_t, 'unit') else float(configured_t))

    poses = _synthetic_poses(residue)
    print('gui_settle_poses.check: %s, %d synthetic poses, map_decouple=%s'
          % (residue, len(poses), map_decouple))
    result = settle_poses(session, residue, poses, map_decouple=map_decouple, debug=debug)

    # ---- programmatic assertions ----
    fails = []
    def ck(cond, msg):
        print(('  ok: ' if cond else '  FAIL: ') + msg)
        if not cond:
            fails.append(msg)

    env_after = sh.atoms.coords[env_mask]
    ck(np.allclose(env_before, env_after, atol=1e-4),
       'environment atoms unchanged by the commit (graft moved only the residue)')
    ck(sh.sim_running, 'simulation still running after settle_poses')
    live_t = sh.temperature
    ck(abs(float(live_t) - configured_t) < 1e-6,
       'temperature restored to configured (%.3g K, now %.3g K)' % (configured_t, float(live_t)))
    ck(_mdff_enableds_all_true(sh, residue),
       'ligand MDFF terms re-enabled after the map_decouple run')
    ck(result is not None and len(result.energies) >= 1, 'result carries ranked energies')
    if result is not None:
        print('  committed=%r kept_current=%s n_culled=%d margin=%.2f kJ/mol'
              % (result.committed, result.kept_current, result.n_culled,
                 result.accept_margin_kJ))
        print('  energies (best-first): '
              + ', '.join('%s %.2f' % (l, e) for l, e in result.energies))

    print('\n%s (%d assertion failure(s)). Also confirm NO "reinitializing context" in the '
          'log above.' % ('PASS' if not fails else 'CHECK FAILED', len(fails)))
    return result


def ab(session, residue=None):
    '''Localisation A/B: for each synthetic pose, settle it and log the LOCALISED E_pack
    alongside the WHOLE-CONSTRUCT CORE energy. E_pack should cleanly order clean-vs-clashing
    poses, while the whole-construct energy is dominated by O(N^2) environment noise and does
    not. Requires a paused sim; leaves the residue on its original coords afterward.'''
    from chimerax.isolde.refine import settle_poses as _sp
    from chimerax.isolde.openmm.openmm_interface import CORE_FORCE_GROUPS
    residue = _selected_residue(session, residue)
    if residue is None:
        return
    isolde = session.isolde
    sh = isolde.sim_handler
    if sh is None or not sh.sim_running or not sh.pause:
        print('gui_settle_poses.ab: needs a PAUSED running simulation.')
        return

    prep = _sp._prepare(sh, residue, residue.atoms, 1.0)
    ridx, base, groups = prep['ridx'], prep['base'], prep['groups']
    moved = residue.atoms
    poses = [('(current)', base[ridx].copy())] + _synthetic_poses(residue)
    print('gui_settle_poses.ab: %s -- E_pack (localised) vs whole-construct CORE' % residue)
    print('  %-12s %14s %18s' % ('pose', 'E_pack', 'whole-CORE'))
    try:
        for label, pose_coords in poses:
            e_pack, _settled = _sp._settle_one(sh, base, ridx, pose_coords, moved, groups,
                                               _sp.SETTLE_STEPS, True, _sp.SETTLE_LAMBDA)
            whole = sh.potential_energy(CORE_FORCE_GROUPS)
            print('  %-12s %14.2f %18.1f' % (label, e_pack, whole))
    finally:
        _sp._restore_coupling(session, sh, moved)
        sh.push_coords_to_sim(base, immediate=True)   # restore original coords
        sh.atoms.coords = base
    print('  Expect: E_pack ranks clash_2.5A worst and (current)/nudge best; the whole-CORE\n'
          '  column is dominated by environment noise and need not order sensibly.')


def check_live(session, residue=None, map_decouple=True):
    '''Live (watchable) driver: run the poses through the GUI loop so you can watch each one
    settle, then score + commit at the end. The sim may be running OR paused on entry. Prints
    the SettleResult from the on_done callback when it finishes (a few seconds).'''
    from chimerax.isolde.refine import settle_poses
    residue = _selected_residue(session, residue)
    if residue is None:
        return
    sh = session.isolde.sim_handler
    if sh is None or not sh.sim_running:
        print('gui_settle_poses.check_live: no running simulation.')
        return
    poses = _synthetic_poses(residue)

    def _done(result):
        if result is None:
            print('gui_settle_poses.check_live: FAILED (see exception above).')
            return
        print('gui_settle_poses.check_live DONE: committed=%r kept_current=%s n_culled=%d'
              % (result.committed, result.kept_current, result.n_culled))
        print('  energies: ' + ', '.join('%s %.2f' % (l, e) for l, e in result.energies))
        print('  Confirm: temperature restored, live-map recalc resumed, sim still running.')

    print('gui_settle_poses.check_live: watch %s settle %d pose(s)...'
          % (residue, len(poses)))
    settle_poses(session, residue, poses, map_decouple=map_decouple, live=True,
                 on_done=_done, debug=True)
