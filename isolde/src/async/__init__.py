
class AsyncTask:
    def __init__(self, session, initiator_func, initiator_args, initiator_kwargs, ready_func, result_getter, callback):
        self.session = session
        self.initiator_func = initiator_func
        self.initiator_args = initiator_args
        self.initiator_kwargs = initiator_kwargs
        self.ready_func = ready_func
        self.result_getter = result_getter
        self.callback = callback
    
    @property
    def ready(self):
        return self.ready_func()

    def task_done(self):
        results = self.result_getter()
        self.callback(self.session, *results)
    
class AsyncTaskMgr:
    '''
    Manager for running asynchronous (i.e. threaded) tasks in response to atomic changes.
    Tasks will be started on submission to the manager, and joined on firing of the ChimeraX
    atomic 'changes done' trigger. Each task must be packaged as a :py:class:`AsyncTask` object,
    which requires (a) an initiator function (and its args and kwargs); (b) a function returning 
    True/False to indicate when the result is ready; (c) a function to collect the results; and (d) 
    a callback taking and applying the results back to the ChimeraX session.
    '''
    def __init__(self, session):
        self.session = session
        session._isolde_task_mgr = self
        from chimerax.atomic import get_triggers
        t = self._atomic_triggers = get_triggers()
        self._trigger = t.add_handler('changes done', self.changes_done_cb)
        self.tasks = []
    
    def submit_task(self, task):
        task.initiator_func(*task.initiator_args, **task.initiator_kwargs)
        self.tasks.append(task)
    
    def changes_done_cb(self, *_):
        tasks = self.tasks
        while len(tasks):
            for i,task in enumerate(tasks):
                if task.ready:
                    task.task_done()
                    tasks.pop(i)
                    break
            

def get_async_task_mgr(session):
    return getattr(session, '_isolde_task_mgr', AsyncTaskMgr(session))