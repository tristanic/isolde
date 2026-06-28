
_server = None

def _get_server():
    global _server
    if _server is None:
        return None
    if _server.httpd is None:
        print('Server httpd is None!')
        _server = None
    return _server

def start_server(session, port=None, allow_run=False):
    from chimerax.core.commands import run
    run(session, 'isolde start')
    global _server
    server = _get_server()
    if server is not None:
        session.logger.warning('ISOLDE REST server is already running')
    else:
        import secrets
        token = secrets.token_urlsafe(32)
        from .server import IsoldeRESTServer
        _server = IsoldeRESTServer(session, auth_token=token, allow_run=bool(allow_run))
        _server.start(port)
        _report_token(session, _server)
    return server

def _report_token(session, server):
    log = session.logger
    log.info(
        'ISOLDE REST server authentication token (required as '
        '"Authorization: Bearer <token>" on every request):'
    )
    # Bold so it is easy to copy from the log; localhost-only, so logging is ok.
    log.info('    {}'.format(server.auth_token), is_html=False)
    if server.allow_run:
        log.warning(
            'ISOLDE REST server started with allowRun=true: the free-text "run" '
            'method (arbitrary ChimeraX commands) is ENABLED. Use only for trusted '
            'local automation.'
        )

def report_info(session):
    server=_get_server()
    addr = server.server_address if _server else None
    if addr is None:
        session.logger.info('ISOLDE REST server is not running')
    else:
        session.logger.info('ISOLDE REST server is listening on host {} port {}'.format(*addr))
        session.logger.info('Authentication token: {}'.format(server.auth_token))
        session.logger.info('Free-text "run" method enabled: {}'.format(server.allow_run))

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
                ('allow_run', BoolArg),
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
