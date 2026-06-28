'''GUI bootstrap for the live collaboration test. Launched at ChimeraX GUI
startup (no --exit). Starts ISOLDE, loads the bundled demo, starts the REST
server, and writes {port, token} to a file once the server is bound AND a model
has loaded -- so an external client can connect and immediately find the model.'''
import json
import os

OUTFILE = os.environ.get('ISOLDE_REST_INFO_FILE', 'rest_info.json')


def bootstrap(session):
    from chimerax.core.commands import run
    run(session, 'isolde start')
    run(session, 'isolde demo crystal_intro')
    run(session, 'isolde remote rest start')

    from chimerax.isolde.remote_control.rest_server import cmd as rc
    from chimerax.core.triggerset import DEREGISTER
    from chimerax.atomic import AtomicStructure

    def _ready(*_):
        s = rc._server
        if s is None or s.httpd is None:
            return
        # Wait until the demo's atomic model has actually loaded.
        has_model = any(isinstance(m, AtomicStructure) for m in session.models.list())
        if not has_model:
            return
        with open(OUTFILE, 'w') as f:
            json.dump({'port': s.port, 'token': s.auth_token}, f)
        session.logger.info('REST info written to %s (port %s)' % (OUTFILE, s.port))
        return DEREGISTER

    session.triggers.add_handler('new frame', _ready)


try:
    session  # noqa
except NameError:
    session = None
if session is not None:
    bootstrap(session)
