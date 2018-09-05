
def delayed_reaction(triggerset, trigger_name, initiator_func, initiator_args, ready_test_func,
    final_func, final_func_args):
    '''
    Designed to work together with threaded (Python or C++) objects, to
    start a long-running threaded task and then automatically apply the
    result (in a later GUI update) when done. Can also be used to simply
    conveniently call the desired callback once on the next firing of the
    trigger, by setting ready_test_func to None.
    Args:
        triggerset: the triggerset providing the trigger (e.g. session.triggers)
        trigger_name: the name of the trigger (e.g. 'new frame')
        initiator_func: A handle to the function that kicks off the threaded
            process. Should not return anything.
        initiator_args: A tuple of arguments to be applied to initiator_func
        ready_test_func: Should return True when the threaded task is done,
            false otherwise. Set it to None to just run on the next firing
            of the trigger.
        final_func: Task to run once the thread is done.
        final_func_args: A tuple of arguments to be applied to final_func
            (e.g. to tell it what to do with the result)
    '''
    initiator_func(*initiator_args)
    class _cb:
        def __init__(self, triggerset, trigger_name, ready_test_func, final_func, final_func_args):
            self.tf = ready_test_func
            self.ff = final_func
            self.ff_args = final_func_args
            self.handler = triggerset.add_handler(trigger_name, self.callback)
        def callback(self, *_):
            if self.tf is None or self.tf():
                self.ff(*self.ff_args)
                from chimerax.core.triggerset import DEREGISTER
                return DEREGISTER
    cb = _cb(triggerset, trigger_name, ready_test_func, final_func, final_func_args)

class Threaded_Contour_Test:
    def __init__(self, session):
        self.session=session
        from chimerax.map import Volume
        volumes = []
        for m in session.models.list():
            if isinstance(m, Volume):
                volumes.append(m)
        if not len(volumes):
            raise RuntimeError("No Volume instances found!")
        self.stop = False
        self._count=0
        self._accumulated_time = 0

        for v in volumes:
            for s in v.surfaces:
                self._iterate(v, s)

    def _iterate(self, volume, surface):
        from random import random
        m = volume.matrix()
        level = random()*(m.max()-m.min())+m.min()
        vertex_transform = volume.matrix_indices_to_xyz_transform()
        normal_transform = vertex_transform.inverse().transpose().zero_translation()
        det = vertex_transform.determinant()

        from chimerax.clipper.contour_thread import Contour_Thread_Mgr
        cm = Contour_Thread_Mgr()
        delayed_reaction(self.session.triggers, 'new frame',
            cm.start_compute, (m, level, det, vertex_transform, normal_transform, False, True),
            cm.ready,
            self._thread_done_cb, (cm, volume, surface, level))

    def _thread_done_cb(self, thread_mgr, volume, surface, level):
        from time import time
        start = time()
        va, ta, na = thread_mgr.get_result()
        surface._set_surface(va, na, ta, None)
        # va, na, ta, hidden_edges = surface._adjust_surface_geometry(va, na, ta, volume.rendering_options, level)
        # surface._set_surface(va, na, ta, hidden_edges)
        self._accumulated_time+=time()-start
        self._count += 1
        if not self.stop:
            self._iterate(volume, surface)
