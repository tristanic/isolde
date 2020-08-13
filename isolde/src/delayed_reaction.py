# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 08-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll

def call_after_n_events(triggerset, trigger_name, n, callback, callback_args):
    class _cb:
        def __init__(self, triggerset, trigger_name, n, cb, cb_args):
            self.count = 0
            self.max_count = n
            self.cb = cb
            self.cb_args = cb_args
            triggerset.add_handler(trigger_name, self.callback)
        def callback(self, *_):
            c = self.count = self.count + 1
            if c >= self.max_count:
                self.cb(*self.cb_args)
                from chimerax.core.triggerset import DEREGISTER
                return DEREGISTER
    _cb(triggerset, trigger_name, n, callback, callback_args)


def delayed_reaction(triggerset, trigger_name, initiator_func, initiator_args, ready_test_func,
    final_func, final_func_args):
    '''
    Designed to work together with threaded (Python or C++) objects, to start a
    long-running threaded task and then automatically apply the result (in a
    later GUI update) when done. Can also be used to simply conveniently call
    the desired callback once on the next firing of the trigger, by setting
    ready_test_func to None.

    Args:

        * triggerset:
            - the :class:`chimerax.core.triggers.TriggerSet` (or other class
              with similar functionality) providing the trigger
              (e.g. `session.triggers`)
        * trigger_name:
            - the name of the trigger (e.g. 'new frame')
        * initiator_func:
            - A handle to the function that kicks off the threaded process.
              Any return value will be ignored.
        * initiator_args:
            - A tuple of arguments to be applied to `initiator_func`
        * ready_test_func:
            - Should return `True` when the threaded task is done, `False`
              otherwise. Set it to `None` to just run on the next firing of the
              trigger.
        * final_func:
            - Task to run once the thread is done.
        * final_func_args:
            - A tuple of arguments to be applied to `final_func` (e.g. to tell
              it what to do with the thread result)

    A simple example using the :class:`Delayed_Reaction_Tester` below to print
    some output to the ChimeraX log after a delay of 100 frames (assuming
    the code is run from the ChimeraX Shell or some other environment where
    `session` is already defined):

    .. code-block:: python

        dt = Delayed_Reaction_Tester()
        delayed_reaction(session.triggers, 'new frame', dt.initiator_func,
            (100,), dt.ready_test_func, dt.final_func, ('Finished!,'))


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

class Delayed_Reaction_Tester:
    def __init__(self):
        pass

    def initiator_func(self, frame_delay):
        self._frame_delay = frame_delay
        self._counter = 0

    def ready_test_func(self):
        self._counter += 1
        if (self._counter >= self._frame_delay):
            return True
        return False

    def final_func(self, to_print):
        print(to_print)
