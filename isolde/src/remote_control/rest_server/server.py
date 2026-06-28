from http.server import BaseHTTPRequestHandler
from chimerax.core.tasks import Task
from chimerax.core.logger import PlainTextLog
import json

# The render/image-capture endpoint is fully implemented (see render.py) but
# PARKED: a Qt 6.11 regression makes any offscreen capture break the main
# window's resize-reshaping (reproducible in vanilla ChimeraX 1.13.dev with
# `save`, absent in 1.12 / Qt 6.10). Flip to True once Qt is fixed, and uncomment
# the isolde_render tool in src/mcp/server.py.
ENABLE_RENDER = False


class IsoldeRESTServer(Task):
    '''
    Listen for HTTP/REST requests, and return machine-friendly
    JSON descriptors of results.
    '''

    def __init__(self, *args, auth_token=None, allow_run=False, **kw):
        self.httpd = None
        # Bearer token required on every request (loopback is still reachable by
        # other local processes, including the MCP server itself).
        self.auth_token = auth_token
        # The free-text 'run' method (arbitrary ChimeraX commands) is the one
        # acknowledged arbitrary-code path; off unless explicitly enabled.
        self.allow_run = allow_run
        super().__init__(*args, **kw)
        self._server_methods = {}
        from .. import default_server_methods, thread_safe
        self.standard_functions = dict(default_server_methods)
        if not allow_run:
            self.standard_functions.pop('run', None)
        # Agent transport plane (see agent_methods.py): typed tool manifest,
        # invoke-by-command, model inventory, atom-spec resolution.
        from .agent_methods import (
            agent_tools, agent_invoke, agent_philosophy, list_models, resolve_spec,
            describe_model, problem_zones, restraint_summary, residue_info, map_info,
            set_mdff, bfactor_outliers, view_state, tug)
        from .jobs import job_start, job_poll, job_result
        self.standard_functions.update({
            'agent_tools': agent_tools,
            'agent_philosophy': agent_philosophy,
            'agent_invoke': agent_invoke,
            'list_models': list_models,
            'describe_model': describe_model,
            'resolve_spec': resolve_spec,
            'problem_zones': problem_zones,
            'restraint_summary': restraint_summary,
            'residue_info': residue_info,
            'map_info': map_info,
            'set_mdff': set_mdff,
            'bfactor_outliers': bfactor_outliers,
            'view_state': view_state,
            'tug': tug,
            'job_start': job_start,
            'job_poll': job_poll,
            'job_result': job_result,
        })

        # `render` is parked (Qt 6.11 capture-breaks-resize regression); only
        # registered when ENABLE_RENDER is True. It also manages its own UI-thread
        # marshalling (it must capture inside the redraw loop, on the 'frame drawn'
        # trigger), so it is registered RAW; everything else is thread_safe-wrapped.
        render_fn = None
        if ENABLE_RENDER:
            from .render import render_view
            render_fn = render_view
        for fname, func in self.standard_functions.items():
            self.register_server_method(fname, thread_safe(func))
        if render_fn is not None:
            self.register_server_method('render', render_fn)
        self._server_methods['batch'] = self.batch_run


    SESSION_SAVE = False

    def batch_run(self, session, batch_commands):
        '''
        Run a series of commands in a single call, to avoid communication delays
        and redrawing between calls.

        Args:

            batch_commands: a list, where each entry in the list is:
                ['func_name', [args], {kwargs}]
        '''
        session = self.session
        from queue import Queue
        from chimerax.core.logger import StringPlainTextLog
        q = Queue()
        def f(session, batch_commands, q=q):
            ret_dict = {}
            for cmd in batch_commands:
                with StringPlainTextLog(session.logger) as rest_log:
                    fname, args, kwargs = cmd
                    f_dict = ret_dict[fname] = {}
                    func = self.standard_functions.get(fname, None)
                    if func is None:
                        err_msg = 'Unrecognised function name: {}'.format(fname)
                        ret_dict['error'] = err_msg
                        break
                    f_dict.update(func(session, *args, **kwargs))
                    f_dict['log'] = rest_log.getvalue()
            q.put(ret_dict)
        session.ui.thread_safe(f, session, batch_commands)
        return q.get()

    @property
    def server_address(self):
        if self.httpd is not None:
            return self.httpd.server_address
        return None

    @property
    def port(self):
        return self.server_address[1]

    def run(self, port):
        from http.server import HTTPServer
        import sys
        if port is None:
            # Defaults to any available port
            port = 0

        httpd = self.httpd = HTTPServer(('localhost', port), RESTHandler)
        httpd.chimerax_session = self.session
        httpd.manager = self
        msg = "ISOLDE REST server started on host {} port {}".format(
            *httpd.server_address
        )
        self.session.ui.thread_safe(self.session.logger.info, msg)
        self._restore_handler = self.session.triggers.add_handler('begin restore session', self._session_restore_cb)
        httpd.serve_forever()

    def terminate(self):
        if self.httpd is not None:
            self.httpd.shutdown()
            self.httpd = None
            if hasattr(self, '_restore_handler') and self._restore_handler is not None:
                self.session.triggers.remove_handler(self._restore_handler)
        super().terminate()

    def register_server_method(self, func_name, func):
        '''
        Register a method to be available as a remote command. The method should
        take the ChimeraX session as its first argument, and the remaining
        arguments should be JSON-serialisable types (e.g. (lists of) strings or
        numbers). For automatic documentation on the client side, the function
        arguments should be annotated, e.g.:

        def f(a: 'string', b: 'int', c:'float')

        The return type should be a JSON-serialisable dict.
        '''
        if func_name not in self._server_methods:
            self._server_methods[func_name] = func

    @property
    def server_methods(self):
        return self._server_methods

    def list_server_methods(self):
        import inspect
        from collections import defaultdict
        ret_dict = defaultdict(lambda: dict())
        for func_name, func in self.server_methods.items():
            func_dict = ret_dict[func_name]
            func_dict['docstring'] = inspect.getdoc(func)
            arg_dict = func_dict['args'] = dict()
            kwarg_dict = func_dict['kwargs'] = dict()
            sig = inspect.signature(func)
            for arg_name, ap in sig.parameters.items():
                if arg_name not in ('self', 'session'):
                    default = ap.default
                    if default != ap.empty:
                        arg_props = kwarg_dict[arg_name] = dict()
                        arg_props['default'] = default
                    else:
                        arg_props = arg_dict[arg_name] = dict()
                    annot = ap.annotation
                    if annot == ap.empty:
                        arg_props['type'] = 'unspecified'
                    else:
                        arg_props['type'] = annot
            # Convert args description to tuple to ensure ordering in Python 2.7 clients
            func_dict['args'] = tuple([(arg_name, arg_props['type']) for arg_name, arg_props in func_dict['args'].items()])
        return ret_dict

    def _session_restore_cb(self, *_):
        from .cmd import stop_server
        stop_server(self.session)

    def take_snapshot(self, session, flags):
        data = {
            'port':     self.port,
        }
        return data

    @staticmethod
    def restore_snapshot(session, data):
        from . import cmd
        return cmd.start_server(session, port=data['port'])


class RESTHandler(BaseHTTPRequestHandler):
    '''Process one REST request.'''

    def log_message(self, format, *args):
        # Suppress http.server's per-request access logging (it was being routed
        # to the ChimeraX log on every agent call). Errors are still returned to
        # the client as structured JSON.
        pass

    def _set_headers(self, status_code=200):
        self.send_response(status_code)
        self.send_header('Content-type', 'application/json')
        self.end_headers()

    def _authorized(self):
        '''True if the request carries the correct bearer token.'''
        token = self.server.manager.auth_token
        if not token:
            return True  # no token configured -> open (legacy/testing only)
        import hmac
        header = self.headers.get('Authorization', '') or ''
        return hmac.compare_digest(header.strip(), 'Bearer ' + token)

    def _reject_unauthorized(self):
        self.send_response(401)
        self.send_header('Content-type', 'application/json')
        self.send_header('WWW-Authenticate', 'Bearer')
        self.end_headers()
        self.wfile.write(json.dumps(
            {'error': 'missing or invalid bearer token'}).encode('utf-8'))

    def do_HEAD(self):
        if not self._authorized():
            return self._reject_unauthorized()
        self._set_headers()

    def do_GET(self):
        '''
        Return a JSON dict listing all available methods
        '''
        if not self._authorized():
            return self._reject_unauthorized()
        self._list_methods()

    def _list_methods(self):
        mgr = self.server.manager
        self._set_headers()
        msg = json.dumps(mgr.list_server_methods())
        self.wfile.write(msg.encode('utf-8'))


    def do_POST(self):
        if not self._authorized():
            return self._reject_unauthorized()
        # Parse the content-type without the deprecated/removed cgi module.
        ctype = (self.headers.get('content-type') or '').split(';')[0].strip()

        # refuse to receive non-json requests
        if ctype != 'application/json':
            self.send_response(400)
            self.end_headers()
            return

        l = int(self.headers.get('content-length'))
        request = json.loads(self.rfile.read(l).decode('utf-8'))
        return_dict = {}
        try:
            return_dict = self._run_post_job(request)
        except Exception as e:
            import traceback
            err_dict = {'error': str(e),
                'traceback': traceback.format_exc()
                }
            self._set_headers(400)
            self.wfile.write(json.dumps(err_dict).encode('utf-8'))
            return
        self._set_headers()
        self.wfile.write(json.dumps(return_dict).encode('utf-8'))




    def _run_post_job(self, request_dict):
        try:
            func_name = request_dict['cmd']
        except KeyError:
            err_dict = {'error': 'You must provide a command name with the key "cmd"!'}
            return err_dict
        mgr = self.server.manager
        f = mgr.server_methods.get(func_name, None)
        if f is None:
            err_dict = {'error': 'No registered server method with the name {}'.format(func_name)}
            return err_dict
        args = request_dict.get('args', [])
        kwargs = request_dict.get('kwargs', {})
        return f(mgr.session, *args, **kwargs)
