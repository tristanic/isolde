def _model_id_as_tuple(model_id):
    return tuple(int(id) for id in model_id.split('.'))

def _model_from_id(session, model_id, error_if_not_found=True):
        m = session.models.list(model_id=_model_id_as_tuple(model_id))
        if not m:
            if error_if_not_found:
                raise TypeError('No model with ID {} found!'.format(model_id))
            else:
                return None
        return m[0]

def _get_symmetry_handler(model):
    from chimerax.atomic import AtomicStructure
    from chimerax.clipper import get_symmetry_handler
    from chimerax.clipper.symmetry import SymmetryManager
    if isinstance(model, AtomicStructure):
        return get_symmetry_handler(model)
    elif isinstance(model, SymmetryManager):
        return model
    else:
        err_string = ('The ID {} corresponds to an incompatible model type: {}'.format(model.id_string, type(model)))
        raise TypeError(err_string)

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
    from chimerax.open_command.cmd import provider_open
    from chimerax.clipper import get_symmetry_handler

    m = provider_open(session, [file_path])[0]
    sh = get_symmetry_handler(m, create=True, auto_add_to_session=True)
    return {'manager': sh.id_string, 'model id': m.id_string}

def update_model_from_file(session, model_id:'string', file_path:'string'):
    '''
    Replace a previously-opened model with new coordinates from a file.

    Args:

        model_id: the ID of the existing model to be replaced
        file_path: the full path to a single PDB or mmCIF file

    Returns:

        {'manager': model id for the top-level Clipper manager for the model,
         'model': model id for the atomic structure itself.
         }
    '''
    m = _model_from_id(session, model_id)
    sh = _get_symmetry_handler(m)
    new_m = sh.swap_model_from_file(file_path)
    return {'manager': sh.id_string, 'model id':new_m.id_string}



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
    from chimerax.clipper.symmetry import SymmetryManager
    from chimerax.atomic import AtomicStructure
    m = _model_from_id(session, model_id)
    if isinstance(m, AtomicStructure):
        from chimerax.clipper import get_symmetry_handler
        sh = get_symmetry_handler(m, create=True, auto_add_to_session=True)
    elif isinstance(m, SymmetryManager):
        sh = m
    else:
        err_string = ('Model ID {} has unrecognised type: {}. '
            'Should be one of AtomicStructure or SymmetryManager.').format(
                model_id, type(m)
            )
        raise TypeError(err_string)

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
    m = _model_from_id(session, model_id)
    from chimerax.atomic import AtomicStructure
    from chimerax.clipper.symmetry import SymmetryManager
    from chimerax.clipper.maps import MapMgr
    if isinstance(m, AtomicStructure):
        from chimerax.clipper import get_map_mgr
        mmgr = get_map_mgr(m)
    elif isinstance(m, SymmetryManager):
        mmgr = m.map_mgr
    elif isinstance(m, MapMgr):
        mmgr = m
    else:
        raise RuntimeError('Model ID {} has unrecognised type {}. Should be one of [{}]'.format(
            model_id, str(type(m)), ', '.join(('AtomicStructure', 'SymmetryManager', 'MapMgr'))
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
        run(session, 'clipper spotlight', log=False)
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
        m = _model_from_id(session, mid, error_if_not_found=False)
        if not m:
            session.logger.warning('Model ID {} not found!'.format(mid))
            not_found.append(mid)
            continue
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
    (NOT YET IMPLEMENTED) Set the display style and colour(s) for a map.

    Args:
        map_id: the model ID for the map to be adjusted
        style: either 'mesh' or 'surface', or None to keep current setting
        color: if the map has a single contour, this should be a list of four
               floats defining [red, green, blue, opacity] in the range 0..1.
               For a map with two contours (i.e. a difference map) provide two
               lists: [[r,g,b,a]]
    '''
    raise RuntimeError('This method is not yet implemented')
