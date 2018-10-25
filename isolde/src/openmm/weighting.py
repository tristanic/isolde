'''
Automatically guess MDFF weighting for a map
'''

def model_volume(model, exclude_hydrogens=True):
    from math import pi
    atoms = model.atoms
    if exclude_hydrogens:
        atoms = atoms[atoms.element_names != 'H']
    return sum( 4/3 * pi * (atoms.radii)**3 )

def guess_mdff_weight(mdff_mgr, percentile=80, constant = 1):
    '''
    Uses a simple heuristic to make a starting guess at the useful weighting of
    a map for MDFF.
    Finds the contour level that would make the enclosed volume (after masking
    to the vicinity of the model) approximately equal to the  volume of the
    model itself. The final weighting is (constant / contour_level)
    '''
    import numpy
    session = mdff_mgr.session
    volume = mdff_mgr.volume
    model = mdff_mgr.model

    from chimerax.clipper.crystal_exp import XmapHandler_Live
    is_xmap = isinstance(volume, XmapHandler_Live)
    if is_xmap:
        sh = volume.manager.crystal

        spotlight_mode = sh.spotlight_mode
        # Expand the map to cover the whole model
        sh.isolate_and_cover_selection(model.atoms, focus=False)

    vs = volume.voxel_size
    voxel_vol = vs[0]*vs[1]*vs[2]
    model_vol = model_volume(model)
    # Consider only grid points within 3A of model atoms
    from chimerax.core.geometry import find_close_points, identity
    # xyz coordinates of each grid point
    grid_coords = volume.grid_points(volume.scene_position)


    data = numpy.array(volume.data.matrix(), order='C')
    close_i = find_close_points(model.atoms.scene_coords, grid_coords, 3)[1]
    close_points = grid_coords[close_i]
    gradients, oob = volume.interpolated_gradients(close_points, identity(),
        out_of_bounds_list=True)
    gradient_mags = numpy.linalg.norm(gradients, axis=1)
    ref_g = numpy.percentile(gradient_mags, percentile)

    if is_xmap:
        sh.spotlight_mode = spotlight_mode

    #return close_grid_vals
    return constant/ref_g
