# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll

def focus_on_coord(session, center, radius = 5.0, clip=True):
    v = session.view
    import numpy
    from chimerax.core.geometry import Bounds
    xyz_min = center-radius
    xyz_max = center+radius
    bounds = Bounds(xyz_min, xyz_max)
    cofr_method = v.center_of_rotation_method
    v.view_all(bounds)
    v.center_of_rotation = center
    v.center_of_rotation_method = cofr_method
    cam = v.camera
    vd = cam.view_direction()
    if clip:
        cp = v.clip_planes
        cp.set_clip_position('near', center - radius*vd, cam)
        cp.set_clip_position('far', center + radius*vd, cam)


def focus_on_selection(session, atoms, pad=5.0, clip = True):
    v = session.view
    pad = 5.0
    bounds = atoms.scene_bounds
    focus_on_coord(bounds,center(), bounds.radius()+pad)
    session.selection.clear()
    atoms.selected=True
