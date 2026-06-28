'''GUI test: chiral restraints survive a checkpoint save/revert round-trip.

Frame-driven state machine: start ISOLDE -> open 1crn -> start a sim -> ensure
chiral restraints exist on the mobile atoms -> snapshot (enableds, spring
constants, cutoffs) -> isolde.checkpoint() -> perturb all three -> isolde.
revert_to_checkpoint() -> re-read and compare. Writes a JSON verdict to
CHIRAL_CP_TEST_OUT and stops the sim.

Run in the GUI:  run_chimerax.bat --script src/tests/gui_chiral_checkpoint_test.py
'''
import os, json, numpy
from chimerax.core.commands import run
from chimerax.core.triggerset import DEREGISTER

OUT = os.environ.get('CHIRAL_CP_TEST_OUT', 'chiral_cp_result.json')


def go(session):
    st = {'phase': 'await_model', 'count': 0}

    def write(result):
        with open(OUT, 'w') as f:
            json.dump(result, f)
        session.logger.info('chiral checkpoint test result: %s' % result)

    # Breadcrumb: prove go() ran at all, and capture any setup-time failure
    # (e.g. isolde start / model open) instead of hanging silently.
    write({'phase': 'go_started'})
    try:
        run(session, 'isolde start', log=False)
        write({'phase': 'isolde_started'})
        # Sim-ready demo (has hydrogens, no add-H dialog) -- same model the
        # proven gui_tug_atom_to_test uses; the structure loads as #1.2.
        run(session, 'isolde demo crystal_intro', log=False)
        write({'phase': 'model_opened'})
    except Exception:
        import traceback
        write({'error': 'setup', 'traceback': traceback.format_exc()})
        return

    def _mobile(isolde):
        return isolde.sim_manager.sim_construct.mobile_atoms

    def step(*_):
        try:
            return _step()
        except Exception:
            import traceback
            write({'error': 'exception', 'traceback': traceback.format_exc(),
                   'phase': st.get('phase')})
            try:
                run(session, 'isolde sim stop', log=False)
            except Exception:
                pass
            return DEREGISTER

    def _step():
        isolde = getattr(session, 'isolde', None)
        if isolde is None:
            return
        st['count'] += 1
        if st['count'] > 2000:
            write({'error': 'timed out in phase %s' % st['phase']})
            return DEREGISTER

        if st['phase'] == 'await_model':
            from chimerax.atomic import AtomicStructure
            if any(isinstance(m, AtomicStructure) for m in session.models.list()):
                run(session, 'isolde sim start #1.2', log=False)
                st['phase'] = 'await_sim'
            return

        if st['phase'] == 'await_sim':
            sh = getattr(isolde, 'sim_handler', None)
            if isolde.simulation_running and sh is not None and sh.thread_handler is not None:
                from chimerax.isolde import session_extensions as sx
                m = isolde.selected_model
                mobile = _mobile(isolde)
                crm = sx.get_chiral_restraint_mgr(m)
                # Ensure chiral restraints exist on the mobile atoms (create if
                # the sim didn't auto-add them), so we always exercise the path.
                crm.add_restraints_by_atoms(mobile)
                crs = crm.get_restraints_by_atoms(mobile)
                if not len(crs):
                    write({'error': 'no chiral centres / restraints on mobile atoms'})
                    run(session, 'isolde sim stop', log=False)
                    return DEREGISTER
                st['orig_enabled'] = numpy.array(crs.enableds, copy=True)
                st['orig_k'] = numpy.array(crs.spring_constants, copy=True)
                st['orig_cut'] = numpy.array(crs.cutoffs, copy=True)
                st['orig_targets'] = numpy.array(crs.targets, copy=True)
                st['n'] = int(len(crs))
                isolde.checkpoint()
                # Perturb every checkpointed field, INCLUDING targets: flip the
                # handedness of every other centre via flip_target() -- the same
                # mutation 'isolde chiralflip force true' performs.
                crs.enableds = numpy.logical_not(crs.enableds)
                crs.spring_constants = crs.spring_constants + 137.0
                crs.cutoffs = crs.cutoffs + 0.07
                n_flipped = 0
                for i, r in enumerate(crs):
                    if i % 2 == 0:
                        r.flip_target()
                        n_flipped += 1
                st['n_flipped'] = n_flipped
                st['phase'] = 'revert'
                st['count'] = 0
            return

        if st['phase'] == 'revert':
            # Let a few frames pass so the perturbation is fully applied.
            if st['count'] < 5:
                return
            from chimerax.isolde import session_extensions as sx
            m = isolde.selected_model
            mobile = _mobile(isolde)
            crm = sx.get_chiral_restraint_mgr(m)
            isolde.revert_to_checkpoint()
            crs = crm.get_restraints_by_atoms(mobile)
            enabled_ok = bool(numpy.array_equal(crs.enableds, st['orig_enabled']))
            k_ok = bool(numpy.allclose(crs.spring_constants, st['orig_k']))
            cut_ok = bool(numpy.allclose(crs.cutoffs, st['orig_cut']))
            # Handedness restored: targets back to their original sign.
            target_ok = bool(numpy.array_equal(
                numpy.sign(crs.targets), numpy.sign(st['orig_targets'])))
            write({
                'error': None,
                'n_restraints': st['n'],
                'n_flipped': st.get('n_flipped'),
                'enabled_restored': enabled_ok,
                'spring_restored': k_ok,
                'cutoff_restored': cut_ok,
                'target_restored': target_ok,
                'all_restored': enabled_ok and k_ok and cut_ok and target_ok,
            })
            run(session, 'isolde sim stop', log=False)
            return DEREGISTER

    session.triggers.add_handler('new frame', step)


try:
    session  # noqa
except NameError:
    session = None
if session is not None:
    go(session)
