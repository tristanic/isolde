
def delayed_reaction(session, initiator_func, initiator_args, ready_test_func,
    final_func, final_func_args):
    '''
    Designed to work together with threaded (Python or C++) objects, to
    start a long-running threaded task and then automatically apply the
    result (in a later GUI update) when done.
    Args:
        initiator_func: A handle to the function that kicks off the threaded
            process. Should not return anything.
        initiator_args: A tuple of arguments to be applied to initiator_func
        ready_test_func: Should return True when the threaded task is done,
            false otherwise
        final_func: Task to run once the thread is done.
        final_func_args: A tuple of arguments to be applied to final_func
            (e.g. to tell it what to do with the result)
    '''
    initiator_func(*initiator_args)
    class _cb:
        def __init__(self, session, ready_test_func, final_func, final_func_args):
            self.session = session
            self.tf = ready_test_func
            self.ff = final_func
            self.ff_args = final_func_args
            self.handler = session.triggers.add_handler('new frame', self.callback)
        def callback(self, *_):
            if self.tf():
                self.ff(*self.ff_args)
                self.session.triggers.remove_handler(self.handler)
    cb = _cb(session, ready_test_func, final_func, final_func_args)

class Delayed_Reaction_Tester:
    def __init__(self, frame_delay):
        self._counter = 0
        self._frame_delay = frame_delay

    def initiator_func(self):
        pass

    def ready_test_func(self):
        self._counter += 1
        if (self._counter >= self._frame_delay):
            return True
        return False

    def final_func(self, to_print):
        print(to_print)
