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

    SESSION_SAVE = False

    def take_snapshot(self, session, flags):
        pass

    @classmethod
    def restore_snapshot(cls, session, data):
        pass

    @property
    def server_address(self):
        if self.httpd is not None:
            return self.httpd.server_address
        return None

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
        httpd.serve_forever()

    def terminate(self):
        if self.httpd is not None:
            self.httpd.shutdown()
            self.httpd = None
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
        import json
        import inspect
        from collections import defaultdict
        ret_dict = defaultdict(lambda: dict())
        for func_name, func in self.server_methods.items():
            func_dict = ret_dict[func_name]
            func_dict['docstring'] = inspect.getdoc(func)
            arg_dict = func_dict['args'] = dict()
            sig = inspect.signature(func)
            for arg_name, ap in zip(sig.parameters.items()):
                if arg_name != 'session':
                    arg_props = arg_dict[arg_name] = dict()
                    annot = ap.annotation
                    if annot == ap.empty:
                        arg_props['type'] = None
                    else:
                        arg_props['type'] = annot
                    default = ap.default
                    if default == ap.empty:
                        arg_props['required'] = True
                    else:
                        arg_props['required'] = False
                        arg_props['default'] = default
        return ret_dict













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
        pass

    def _list_methods(self):
        mgr = self.server.manager
        self._set_headers()
        self.wfile.write(json.dumps(mgr.list_server_methods()))


    def do_POST(self):
        from cgi import parse_header
        ctype, pdict = cgi.parse_header(self.headers.getheader('content-type'))

        # refuse to receive non-json requests
        if ctype != 'application/json':
            self.send_response(400)
            self.end_headers()
            return

        l = int(self.headers.getheader('content-length'))
        request = json.loads(self.rfile.read(l))
        try:
            return_dict = self._run_post_job(request)
        except Exception as e:
            import traceback
            err_dict = {'error': str(e),
                'traceback': traceback.format_exc()
                }
            self._set_headers(400)
            self.wfile.write(json.dumps(err_dict))
        self._set_headers()
        self.wfile.write(json.dumps(return_dict))




    def _run_post_job(self, request_dict):
        pass
