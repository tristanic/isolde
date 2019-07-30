from http.server import BaseHTTPRequestHandler
from chimerax.core.tasks import Task
from chimerax.core.logger import PlainTextLog
import json

import functools
def thread_safe(func):
    @functools.wraps(func)
    def thread_safe_func(session, *args, **kwargs):
        from queue import Queue
        q = Queue()
        def inner_func(session, *args, q=q, **kwargs):
            from chimerax.core.errors import NotABug
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


class IsoldeRESTServer(Task):
    '''
    Listen for HTTP/REST requests, and return machine-friendly
    JSON descriptors of results.
    '''



    def __init__(self, *args, **kw):
        self.httpd = None
        super().__init__(*args, **kw)
        self._server_methods = {}
        self.standard_functions = {
            'run':  run_chimerax_command,
            'test': test_command,
            'next_model_id': next_id,
            'load_model': load_model,
            'close_models': close_models,
            'load_structure_factors': load_structure_factors,
            'load_map': load_map,
            'center_on_coord': center_on_coord,
            'spotlight_radius': spotlight_radius,
        }

        for fname, func in self.standard_functions.items():
            self.register_server_method(fname, thread_safe(func))
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
        return ret_dict

def _model_id_as_tuple(model_id):
    return tuple(int(id) for id in model_id.split('.'))

def run_chimerax_command(session, commands:'string or list of strings'):
    '''
    Run one or more ChimeraX command-line commands, and receive the resulting
    log messages.

    Args:

        commands: either a single string or a list of strings, where each string
        is a complete, executable ChimeraX command, e.g.:

        ["color #1 bychain","color #1 byhetero","cofr center showPivot true"]
    '''
    if not isinstance(commands, list):
        commands = [commands]
    from chimerax.core.commands import run
    try:
        for cmd in commands:
            run(session, cmd, log=False)
    except NotABug as e:
        logger.info(str(e))
    return {}


def test_command(session, arg1:'whatever', arg2:'you', kwarg1:'like'=None):
    '''
    Will simply echo back a dict of the provided arguments.
    '''
    return {'arg1': arg1, 'arg2':arg2, 'kwarg1': kwarg1}

def next_id(session, parent=None):
    '''
    Get the next ID number that will be applied to a model added to the given
    parent.

    Args:

        parent: the ID of an existing model, or None. If None, the return value
                will be the ID of the next independent model to be added (e.g.
                '1', '2', '3', ...). Otherwise it will be the next available
                submodel specifier (e.g. '1.1.1.3').
    '''
    if parent is not None:
        m = session.models.list(model_id=_model_id_as_tuple(parent))
        if not m:
            raise TypeError('No model with ID {} found!'.format(parent))
        m = m[0]
    else:
        m = None
    return {'next id': '.'.join([str(i) for i in session.models.next_id(parent=m)])}

def load_model(session, file_path:'string'):
    '''
    Load a model from a single PDB or mmCIF file. Multi-model files are not
    supported.

    Args:

        file_path: the full path to a single PDB or mmCIF file

    Returns:

        {'manager': model id for the top-level Clipper manager for the model,
         'model': model id for the atomic structure itself.
         }
    '''
    from chimerax.core.commands import open as cxopen
    from chimerax.clipper import get_symmetry_handler

    m = cxopen.open(session, file_path)[0]
    sh = get_symmetry_handler(m, create=True, auto_add_to_session=True)
    return {'manager': sh.id_string, 'model id': m.id_string}

def load_structure_factors(session, file_path:'string', model_id:'string'):
    '''
    Load a set of structure factors in MTZ or CIF format, generate maps and
    associate them with an existing atomic model. Data may be provided as any
    of the following:

        F / sigF
        I / sigI
        F+ / sigF+ / F- / sigF-
        I+ / sigI+ / I- / sigI-
        Free flags
        F / phi

    Only one experimental dataset should be provided, but any number of F/phi
    maps many be provided. If experimental data is provided, three "live" maps
    will be calculated: a standard 2mFo-DFc map; a second 2mFo-DFc map with a
    resolution-dependent sharpening or smoothing B-factor applied (sharpened
    at resolutions worse than 2.5A, smoothed otherwise); and a mFo-DFc map. For
    best results, the experimental reflections should already be corrected for
    anisotropy and any artefacts such as ice rings or beamstop shadows.
    Anomalous data will be automatically merged, and intensities converted to
    amplitudes using the method of Read & McCoy.

    CAUTION: if no free flags are provided, a new set will be automatically
    generated.

    Any of the generated maps may be closed using close_model() on its id, or
    the whole set may be closed at once by closing the map_mgr.

    Args:

        file_path:  the full path to a single .mtz or .cif file
        model_id:   id string (e.g. as returned by load_model()) of the atomic
                    structure to associate the structure factors with, or its
                    top-level Clipper manager.

    Returns:

        {'manager': model id for the top-level Clipper manager for the model/maps,
         'model':   model id for the atomic structure,
         'map_mgr': model id for the manager of all maps associated with the structure,
         'mapset':  model id for the container holding the maps resulting from this call,
         'maps': {
            'map 1 column names': map 1 model id,
            ...
            }
         }
    '''
    from chimerax.clipper.symmetry import Symmetry_Manager
    from chimerax.atomic import AtomicStructure
    for m in session.models.list():
        if m.id_string == model_id:
            if isinstance(m, AtomicStructure):
                from chimerax.clipper import get_symmetry_handler
                sh = get_symmetry_handler(m, create=True, auto_add_to_session=True)
            elif isinstance(m, Symmetry_Manager):
                sh = m
            else:
                err_string = ('Model ID {} has unrecognised type: {}. '
                    'Should be one of AtomicStructure or Symmetry_Manager.').format(
                        model_id, type(m)
                    )
                raise TypeError(err_string)
            break
    else:
        raise TypeError('The requested model ID {} was not found!'.format(model_id))

    mmgr = sh.map_mgr
    xmapset = mmgr.add_xmapset_from_file(file_path)
    return {
        'manager':  sh.id_string,
        'model':    sh.structure.id_string,
        'map_mgr':  mmgr.id_string,
        'mapset':   xmapset.id_string,
        'maps': {x.name: x.id_string for x in xmapset},
    }

def load_map(session, file_path:'string', model_id:'string'):
    '''
    Load a real-space map in any format that ChimeraX recognises. In ISOLDE,
    each map must be associated with a Clipper data manager object.

    Args:

        file_path: the full or relative path to the map file to open
        model_id: may be the id for an atomic model, a map_mgr or a top-level
            Clipper data manager object

    Returns:

        {
            'manager':  the top-level Clipper manager id,
            'map_mgr':  the id of the map manager for this model,
            'map':      the id of the newly-opened map
        }
    '''
    m = session.models.list(model_id = _model_id_as_tuple(model_id))
    if not len(m):
        raise RuntimeError('No model with ID {} found!'.format(model_id))
    m = m[0]
    from chimerax.atomic import AtomicStructure
    from chimerax.clipper.symmetry import Symmetry_Manager
    from chimerax.clipper.maps import Map_Mgr
    if isinstance(m, AtomicStructure):
        from chimerax.atomic import get_map_mgr
        mmgr = get_map_mgr(m)
    elif isinstance(m, Symmetry_Manager):
        mmgr = m.map_mgr
    elif isinstance(m, Map_Mgr):
        mmgr = m
    else:
        raise RuntimeError('Model ID {} has unrecognised type {}. Should be one of [{}]'.format(
            model_id, str(type(m)), ', '.join(('AtomicStructure', 'Symmetry_Manager', 'Map_Mgr'))
        ))
    nxmapset = mmgr.nxmapset
    new_map = nxmapset.add_nxmap_handler_from_file(file_path)
    return {
        'manager': mmgr.crystal_mgr.id_string,
        'map_mgr': mmgr.id_string,
        'map':     new_map.id_string
    }



def center_on_coord(session, coord:'list', radius:'float' = 5.0,
        spotlight:'bool' = True):
    '''
    Focus the view on an (x,y,z) coordinate and set the zoom so that a circle
    of the given radius is visible.

    Args:

        coord: a list of three floats in Angstroms
        radius: a floating-point value in Angstroms
        spotlight: if True, return view to Clipper's Spotlight mode if not
                   already set. If False, don't change the current setting.

    Returns:
        empty dict
    '''
    if spotlight:
        from chimerax.core.commands import run
        run(session, 'clipper spotlight enable true', log=False)
    from chimerax.isolde.view import focus_on_coord
    focus_on_coord(session, coord, radius, True)
    return {}

def spotlight_radius(session, radius:'float', managers:'list'=[]):
    '''
    Adjust the radius of the "spotlight" (the sphere of density around the
    centre of rotation).

    Args:

        radius: the desired radius in Angstroms
        managers: a list of the IDs corresponding to the Clipper top-level
                  managers for the models for which the radius is to be updated.
                  To update the radius for all currently-loaded models/maps,
                  either omit the argument or provide an empty list.

    Returs:
        empty dict
    '''
    from chimerax.core.commands import run
    if len(managers):
        mgr_string = '#{}'.format(','.join(managers))
    else:
        mgr_string = ''

    run(session, 'clipper spotlight {} radius {}'.format(mgr_string, radius), log=False)
    return {}

def close_models(session, models:'string or list'):
    '''
    Close one or more models and/or maps. Behaviour varies depending on the
    types of models specified. If the model is an atomic structure or a Clipper
    top-level manager, the manager and all maps associated with the model will
    be deleted. If a map, just that map will be deleted. If a map manager, all
    maps handled by that manager will be deleted.

    Args:

        models: a single model ID or a list of model IDs.
    '''
    from chimerax.atomic import AtomicStructure
    if not isinstance(models, list):
        models = [models]
    to_close = []
    not_found = []
    for mid in models:
        t_mid = _model_id_as_tuple(mid)
        mlist = session.models.list(model_id = t_mid)
        if not mlist:
            session.logger.warning('Model ID {} not found!'.format(mid))
            not_found.append(mid)
            continue
        m = mlist[0]
        if isinstance(m, AtomicStructure):
            from chimerax.clipper import get_symmetry_handler
            m = get_symmetry_handler(m)
        to_close.append(m)
    ret = {'closed': [m.id_string for m in to_close],
           'not_found': not_found,
           }
    session.models.close(to_close)
    return ret

def map_style(session, map_id:'string', style:'string or None'=None, color:'list or None'=None):
    '''
    Set the display style and colour(s) for a map.

    Args:
        map_id: the model ID for the map to be adjusted
        style: either 'mesh' or 'surface', or None to keep current setting
        color: if the map has a single contour, this should be a list of four
               floats defining [red, green, blue, opacity] in the range 0..1.
               For a map with two contours (i.e. a difference map) provide two
               lists: [[r,g,b,a]]
    '''
    pass


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
