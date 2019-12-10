
_server = None

def _get_server():
    global _server
    if _server is None:
        return None
    if _server.httpd is None:
        print('Server httpd is None!')
        _server = None
    return _server

def start_server(session, port=None):
    from chimerax.core.commands import run
    run(session, 'isolde start')
    global _server
    server = _get_server()
    if server is not None:
        session.logger.warning('ISOLDE REST server is already running')
    else:
        from .server import IsoldeRESTServer
        _server = IsoldeRESTServer(session)
        _server.start(port)
    return server

def report_info(session):
    server=_get_server()
    addr = server.server_address if _server else None
    if addr is None:
        session.logger.info('ISOLDE REST server is not running')
    else:
        session.logger.info('ISOLDE REST server is listening on host {} port {}'.format(*addr))

def stop_server(session):
    global _server
    server = _get_server()
    if server is None:
        session.logger.info('ISOLDE REST server is not running')
    else:
        server.terminate()
        _server = None
        session.logger.info('ISOLDE REST server stopped')


def register_isolde_rest_server(logger):
    from chimerax.core.commands import (
        register,
        CmdDesc, IntArg, BoolArg,
    )
    def register_start(logger):
        desc = CmdDesc(
            keyword=[('port', IntArg),
                ],
            synopsis='Start ISOLDE REST server'
        )
        register('isolde remote rest start', desc, start_server, logger=logger)
    def register_report_port(logger):
        desc = CmdDesc(synopsis='Report ISOLDE REST server address and port')
        register('isolde remote rest info', desc, report_info, logger=logger)
    def register_stop_server(logger):
        desc = CmdDesc(synopsis='Stop the ISOLDE REST server')
        register('isolde remote rest stop', desc, stop_server, logger=logger)
    register_start(logger)
    register_report_port(logger)
    register_stop_server(logger)
