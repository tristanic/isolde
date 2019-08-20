
class IsoldeRESTClient:
    '''
    Client-side code to interface to ISOLDE's REST API for remote control from
    another program. NOTE: this uses a non-secure HTTP connection which a
    malicious user could in theory use to run arbitrary code on your machine.
    Should be used only on your local machine or a local area network, protected
    by a firewall.

    On running the connect() function, the client instance will be populated
    with methods based on those defined on the server. A dict giving the
    available method names, required/optional arguments and docstrings can be
    obtained with get_available_methods().
    '''
    def __init__(self, address, port, timeout=20):
        self._address = address
        self._port = port
        self._timeout = timeout
        self._headers = {'Content-type': 'application/json'}

    def connect(self):
        import http.client
        self._connection = http.client.HTTPConnection('{}:{}'.format(self._address, self._port), timeout=self._timeout)
        self.get_available_methods()

    def _get_available_methods(self):
        if not self.connected:
            raise RuntimeError('Not connected to server!')
        self._connection.request('GET', '')
        method_dict = self._get_result()
        return method_dict

    def get_available_methods(self):
        method_dict = self._get_available_methods()
        self._method_factory(self._get_available_methods())
        return method_dict

    def _method_factory(self, method_dict):
        cls = self.__class__
        for fname, func_def in method_dict.items():
            if not hasattr(cls, fname):
                args = func_def.get('args', [])
                kwargs = func_def.get('kwargs', {})
                f_str = '''
def {}_server_method(self{}):
    """
    {}
    """
    send_dict = dict()
    args = [{}]
    kwargs = {}

    send_dict['cmd'] = "{}"
    send_dict['args'] = args
    send_dict['kwargs'] = kwargs

    import json
    self._connection.request('POST', '', json.dumps(send_dict).encode('utf-8'), self._headers)
    result = self._get_result()
    if 'error' in result.keys():
        raise RuntimeError(result['error'] + '; Server traceback: \\n' + result.get('traceback', 'None provided'))
    return result
'''
                if args and kwargs:
                    arg_format = ', '+', '.join((
                        ', '.join("{}:'{}'".format(arg, argtype) for (arg, argtype) in args),
                        ', '.join(["{}:'{}'={}".format(kw, val['type'], val['default']) for kw, val in kwargs.items()])
                        ))
                elif args:
                    arg_format = ', '+', '.join("{}:'{}'".format(arg, argtype) for (arg, argtype) in args)
                elif kwargs:
                    arg_format = ', '+', '.join(["{}:'{}'={}".format(kw, val['type'], val['default']) for kw, val in kwargs.items()])
                else:
                    arg_format = ''

                f_str = f_str.format(
                    fname,
                    arg_format,
                    func_def['docstring'],
                    ', '.join([arg_desc[0] for arg_desc in args]),
                    "dict( [{}] )".format(', '.join('("{}", {})'.format(kw, kw) for kw in kwargs.keys())),
                    fname,
                )
                # print(f_str)
                exec(f_str, globals())
                eval('setattr(cls, fname, {}_server_method)'.format(fname), globals(), locals())



    @property
    def connected(self):
        return (hasattr(self, '_connection') and self._connection is not None)
    @property
    def server_address(self):
        return self._address

    @property
    def server_port(self):
        return self._port

    def _get_result(self):
        from cgi import parse_header
        result = self._connection.getresponse()
        ctype, pdict = parse_header(result.headers.get('content-type'))
        if ctype != 'application/json':
            err_string = ('Server returned a non-supported type, "{}". Only "{}" '
                'is allowed.').format(ctype, 'application/json')
            raise TypeError(err_string)
        import json
        return json.loads(result.read().decode('utf-8'))

if __name__=='__main__':
    ic = IsoldeRESTClient('localhost', '12365')
    ic.connect()
    print(ic.run(['open 1a0m structureFactors true']))
