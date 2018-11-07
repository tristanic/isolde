
def compiled_lib_extension():
    import platform
    pname = platform.system()
    if pname == "Windows":
        return "dll"
    elif pname == "Darwin":
        return "dylib"
    return "so"

def available_cores():
    import os
    return max(os.cpu_count()-2, 1)

def set_to_default_cartoon(session, model = None):
    '''
    Adjust the ribbon representation to provide information without
    getting in the way.
    '''
    try:
        if model is None:
            atoms = None
        else:
            arg = atomspec.AtomSpecArg('thearg')
            atoms = arg.parse('#' + model.id_string, session)[0]
        cartoon.cartoon(session, atoms = atoms, suppress_backbone_display=False)
        cartoon.cartoon_style(session, atoms = atoms, width=0.4, thickness=0.1, arrows_helix=True, arrow_scale = 2)
        cartoon.cartoon_tether(session, structures=atoms, opacity=0)
    except:
        return

def atom_list_from_sel(atom_list):
    '''
    Takes a ChimeraX Atoms object, and creates a Clipper Atoms_list object
    from the relevant atom properties.
    '''
    n = len(atom_list)
    elements = atom_list.element_names.tolist()
    coords = atom_list.coords
    occupancies = atom_list.occupancies
    import numpy
    from math import pi
    u_iso = numpy.sqrt(atom_list.bfactors/(8*pi**2))
    u_aniso = numpy.ones([n,6],numpy.double)*numpy.nan
    u_aniso[atom_list.has_aniso_u] = atom_list.filter(atom_list.has_aniso_u).aniso_u6
    from .clipper_python import Atom_list
    clipper_atom_list = Atom_list(elements, coords, occupancies, u_iso, u_aniso)
    return clipper_atom_list

def _model_volume(model, exclude_hydrogens=True, radius_scale=1):
    from math import pi
    atoms = model.atoms
    if exclude_hydrogens:
        atoms = atoms[atoms.element_names != 'H']
    return sum( 4/3 * pi * (atoms.radii*radius_scale)**3 )

def guess_suitable_contour(volume, model, mask_radius=3, atom_radius_scale = 0.5):
    '''
    Find the contour level that would make the volume inside the contour for the
    region immediately surrounding the model approximately equal to the volume
    of the model's atoms scaled by atom_radius_scale.
    '''
    import numpy
    session = model.session

    from .symmetry import is_crystal_map
    is_xmap = is_crystal_map(volume)
    if is_xmap:
        sh = volume.manager.crystal

        spotlight_mode = sh.spotlight_mode
        # Expand the map to cover the whole model
        sh.isolate_and_cover_selection(model.atoms, focus=False)

    vv = volume.data.voxel_volume()
    mv = _model_volume(model, radius_scale=atom_radius_scale)

    from chimerax.core.geometry import find_close_points
    grid_coords = volume.grid_points(volume.scene_position)

    data = numpy.array(volume.data.matrix(), order='C')
    close_i = find_close_points(model.atoms.scene_coords, grid_coords, mask_radius)[1]
    close_vals = data.ravel()[close_i]


    target_percentile = (1-mv/(vv*len(close_i)))*100
    level = numpy.percentile(close_vals, target_percentile)

    if is_xmap:
        sh.spotlight_mode = spotlight_mode

    return level
