# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 11-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def register_remote_commands(logger):
    from .xmlrpc.remotecmd import register_isolde_xmlrpc_server
    register_isolde_xmlrpc_server(logger)
    from.rest_server.cmd import register_isolde_rest_server
    register_isolde_rest_server(logger)

from . import server_methods as _sm

def thread_safe(func):
    import functools
    @functools.wraps(func)
    def thread_safe_func(session, *args, **kwargs):
        from queue import Queue
        q = Queue()
        def inner_func(session, *args, q=q, **kwargs):
            from chimerax.core.logger import StringPlainTextLog
            with StringPlainTextLog(session.logger) as rest_log:
                ret = {}
                try:
                    ret.update(func(session, *args, **kwargs))
                except Exception as e:
                    import traceback
                    ret.update({
                        'error':  str(e),
                        'traceback':  traceback.format_exc()
                    })
                ret['log'] = rest_log.getvalue()
            q.put(ret)
        session.ui.thread_safe(inner_func, session, *args, **kwargs)
        return q.get()
    return thread_safe_func

default_server_methods = {
    'run':  _sm.run_chimerax_command,
    'test': _sm.test_command,
    'next_model_id': _sm.next_id,
    'load_model': _sm.load_model,
    'update_model': _sm.update_model_from_file,
    'close_models': _sm.close_models,
    'load_structure_factors': _sm.load_structure_factors,
    'load_map': _sm.load_map,
    'center_on_coord': _sm.center_on_coord,
    'spotlight_radius': _sm.spotlight_radius,
}
