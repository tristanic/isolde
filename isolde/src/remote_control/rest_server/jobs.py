# @Author: Tristan Croll
# @Date:   27-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 27-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Async job pattern for long-running ISOLDE operations (e.g. brefine / brsr).

Interactive simulation is already non-blocking (driven by the 'new frame'
trigger). Longer operations need a non-blocking **start -> poll -> result**
contract with job ids so the HTTP request never blocks on the UI thread.

This is a general mechanism: ``job_start`` runs *any* agent-safe command in a
background worker, marshaling the actual command execution onto the UI thread
(via ``session.ui.thread_safe``) so it is engine-safe, while immediately
returning a job id to the caller. ``job_poll`` / ``job_result`` report state and
the structured result/witness when done.

Engines that thread themselves internally (the brefine manager spins its own
worker and returns immediately) complete a job as soon as their launch returns;
streaming their live convergence metrics into ``job_poll`` is a follow-up that
needs the command to expose its manager handle (noted in the design doc).
'''
import threading
import uuid


class Job:
    __slots__ = ('id', 'command', 'args', 'state', 'result', 'error', 'started')

    def __init__(self, job_id, command, args):
        self.id = job_id
        self.command = command
        self.args = args
        self.state = 'running'      # running | done | failed
        self.result = None
        self.error = None
        self.started = None         # set by caller (no clock access here)

    def as_dict(self):
        d = {'job_id': self.id, 'command': self.command, 'state': self.state}
        if self.result is not None:
            d['result'] = self.result
        if self.error is not None:
            d['error'] = self.error
        return d


# Process-lifetime registry of jobs.
_JOBS = {}


def _new_job_id():
    return uuid.uuid4().hex[:16]


def job_start(session, name, args=None):
    '''
    Start command *name* (with JSON *args*) as a background job. Returns
    immediately with ``{job_id, state:'running'}``. Poll with ``job_poll``.
    '''
    from chimerax.isolde.cmd.command_registry import get_command
    if args is None:
        args = {}
    # Resolve underscore form and validate exposure up-front so a bad request
    # fails synchronously rather than in the worker.
    canonical = name if get_command(name) is not None else name.replace('_', ' ')
    record = get_command(canonical)
    if record is None:
        return {'error': 'no such command: %r' % name}
    if not record.agent_safe:
        return {'error': '%r is not agent_safe' % canonical}

    job = Job(_new_job_id(), canonical, args)
    _JOBS[job.id] = job

    def _worker(job=job, canonical=canonical, args=args):
        from queue import Queue
        from chimerax.isolde.cmd.invoke import invoke_command
        q = Queue()

        def _on_ui():
            try:
                q.put(('ok', invoke_command(session, canonical, args)))
            except Exception as e:
                import traceback
                q.put(('err', (str(e), traceback.format_exc())))

        # Marshal command execution to the UI thread (engine-safe). In nogui
        # --script this runs directly; in the GUI it hops to the UI thread.
        session.ui.thread_safe(_on_ui)
        kind, payload = q.get()
        if kind == 'ok':
            job.result = payload
            job.state = 'done'
        else:
            job.error = {'error': payload[0], 'traceback': payload[1]}
            job.state = 'failed'

    threading.Thread(target=_worker, name='isolde-job-%s' % job.id,
                     daemon=True).start()
    return {'job_id': job.id, 'state': 'running', 'command': canonical}


def job_poll(session, job_id):
    '''Return current state (and result/error if finished) for *job_id*.'''
    job = _JOBS.get(job_id)
    if job is None:
        return {'error': 'no such job: %r' % job_id}
    return job.as_dict()


def job_result(session, job_id):
    '''
    Return the final result for *job_id*. Does not block: if the job is still
    running, returns its current state so the caller can keep polling.
    '''
    job = _JOBS.get(job_id)
    if job is None:
        return {'error': 'no such job: %r' % job_id}
    return job.as_dict()
