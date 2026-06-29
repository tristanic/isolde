'''GUI test for ISOLDE_mgr.tug_atom_to() (the public method fixed on master).
Frame-driven state machine: load demo -> start sim -> call tug_atom_to on a heavy
atom -> let it settle -> write {error, moved_A, ...} to TUG_TEST_OUT and stop.

Run in the GUI:  run_chimerax.bat --script src/tests/gui_tug_atom_to_test.py
'''
import os, json, numpy
from chimerax.core.commands import run
from chimerax.core.triggerset import DEREGISTER

OUT = os.environ.get('TUG_TEST_OUT', 'tug_test_result.json')


def go(session):
    run(session, 'isolde start', log=False)
    run(session, 'isolde demo crystal_intro', log=False)
    st = {'phase': 'await_model', 'count': 0}

    def write(result):
        with open(OUT, 'w') as f:
            json.dump(result, f)
        session.logger.info('tug_atom_to test result: %s' % result)

    def step(*_):
        isolde = getattr(session, 'isolde', None)
        if isolde is None:
            return
        if st['phase'] == 'await_model':
            from chimerax.atomic import AtomicStructure
            if any(isinstance(m, AtomicStructure) for m in session.models.list()):
                run(session, 'isolde sim start #1.2', log=False)
                st['phase'] = 'await_sim'
            return
        if st['phase'] == 'await_sim':
            sh = getattr(isolde, 'sim_handler', None)
            if isolde.simulation_running and sh is not None and sh.thread_handler is not None:
                m = isolde.selected_model
                cas = m.atoms[m.atoms.names == 'CA']
                a = cas[len(cas) // 2] if len(cas) else m.atoms[m.atoms.element_names != 'H'][0]
                st['atom'] = a
                st['start'] = numpy.array(a.coord, dtype=float)
                target = (st['start'] + numpy.array([3.0, 0.0, 0.0])).tolist()
                st['target'] = target
                try:
                    isolde.tug_atom_to(a, target)
                    st['error'] = None
                except Exception as e:
                    import traceback
                    st['error'] = str(e)
                    st['traceback'] = traceback.format_exc()
                st['phase'] = 'settle'
                st['count'] = 0
            return
        if st['phase'] == 'settle':
            st['count'] += 1
            if st['count'] >= 150:   # ~2-3 s of frames
                a = st['atom']
                now = numpy.array(a.coord, dtype=float)
                tgt = numpy.array(st['target'], dtype=float)
                moved = float(numpy.linalg.norm(now - st['start']))
                to_target_before = float(numpy.linalg.norm(st['start'] - tgt))
                to_target_after = float(numpy.linalg.norm(now - tgt))
                write({'error': st.get('error'), 'traceback': st.get('traceback'),
                       'moved_A': moved,
                       'dist_to_target_before': to_target_before,
                       'dist_to_target_after': to_target_after,
                       'pulled_closer': to_target_after < to_target_before})
                run(session, 'isolde sim stop', log=False)
                return DEREGISTER
            return

    session.triggers.add_handler('new frame', step)


try:
    session  # noqa
except NameError:
    session = None
if session is not None:
    go(session)
