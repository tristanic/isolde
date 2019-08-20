import asyncio
import websockets
from chimerax.core.tasks import Task


class SocketManager:
    def __init__(self, port, address='localhost', is_server=True):
        self.port = port
        self.address = address
        self.is_server = is_server

    async def _serve(self, websocket, path):
        consumer_task = asyncio.ensure_future(
            self._incoming_handler(websocket, path)
        )
        producer_task = asyncio.ensure_future(
            self._outgoing_handler(websocket, path)
        )
        done, pending = await asyncio.wait(
            [consumer_task, producer_task],
            return_when=asyncio.FIRST_COMPLETED
        )
        for task in pending:
            task.cancel()

    async def _incoming_handler(websocket, path):
        async for message in websocket:
            await self._handle_incoming(message)

    async def _outgoing_handler(websocket, path):
        while True:
            message = await self._send_next()
            await websocket.send(message)


class IsoldeSocketServer(Task):
    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self._server_methods = {}
        from .. import default_server_methods, thread_safe
        self.standard_functions = default_server_methods
        for fname, func in self.standard_functions.items():
            self.register_server_method(fname, thread_safe(func))
        self._server_methods['batch'] = self.batch_run

    def run(self, port, address='localhost'):
        self.port = port
        self.address = address
        loop = self._async_loop = asyncio.new_event_loop()
        stop = self._stop_trigger = loop.create_future()
        start_server = self._server = websockets.serve(self._serve, 'localhost', port, loop=loop)
        loop.run_until_complete(start_server)
        loop.run_forever()

    def terminate(self):
        self.stop()
        super().terminate()

    def stop(self):
        self._server.ws_server.close()
        self.session.ui.thread_safe(self.session.logger.info, 'Shutting down server on {} port {}'.format(self.address, self.port))

    async def _serve(self, websocket, cmd):
        cmd = await websocket.recv()
        import json
        cmd_dict = json.loads(cmd)
        if cmd_dict is None:
            await websocket.send(json.dumps(self.list_server_methods()))
        else:
            result = self._run_cmd(cmd_dict)
            await websocket.send(json.dumps(result))

    def _run_cmd(self, cmd_dict):
        try:
            func_name = cmd_dict['cmd']
        except KeyError:
            err_dict = {'error': 'You must provide a command name with the key "cmd"!'}
            return err_dict
        f = self.server_methods.get(func_name, None)
        if f is None:
            err_dict = {'error': 'No registered server method with the name {}'.format(func_name)}
            return err_dict
        args = cmd_dict.get('args', [])
        kwargs = cmd_dict.get('kwargs', {})
        return f(self.session, *args, **kwargs)

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
