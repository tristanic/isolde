from http.server import BaseHTTPRequestHandler
from chimerax.core.tasks import Task
from chimerax.core.logger import PlainTextLog
import json



class IsoldeRESTServer(Task):
    '''
    Listen for HTTP/REST requests, and return machine-friendly
    JSON descriptors of results.
    '''

    def __init__(self, *args, **kw):
        self.httpd = None
        super().__init__(*args, **kw)
        self._server_methods = {}
        from .. import default_server_methods, thread_safe
        self.standard_functions = default_server_methods

        for fname, func in self.standard_functions.items():
            self.register_server_method(fname, thread_safe(func))
        self._server_methods['batch'] = self.batch_run


    SESSION_SAVE = True

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

    def _set_headers(self, status_code=200):
        self.send_response(status_code)
        self.send_header('Content-type', 'application/json')
        self.end_headers()

    def do_HEAD(self):
        self._set_headers()

    def do_GET(self):
        '''
        Return a JSON dict listing all available methods
        '''
        self._list_methods()

    def _list_methods(self):
        mgr = self.server.manager
        self._set_headers()
        msg = json.dumps(mgr.list_server_methods())
        self.wfile.write(msg.encode('utf-8'))


    def do_POST(self):
        from cgi import parse_header
        ctype, pdict = parse_header(self.headers.get('content-type'))

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
