# @Author: Tristan Croll
# @Date:   05-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   Tristan Croll
# @Last modified time: 18-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll


def focus_on_selection(session, view, atoms, clip = True):
    v = view
    pad = 5.0
    bounds = atoms.scene_bounds
    bounds.xyz_min = bounds.xyz_min - pad
    bounds.xyz_max = bounds.xyz_max + pad
    radius = bounds.radius() + pad
    cofr_method = v.center_of_rotation_method
    v.view_all(bounds)
    v.center_of_rotation = center = bounds.center()
    v.center_of_rotation_method = cofr_method
    cam = v.camera
    vd = cam.view_direction()
    if clip:
        cp = v.clip_planes
        cp.set_clip_position('near', center - radius*vd, cam)
        cp.set_clip_position('far', center + radius*vd, cam)
    session.selection.clear()
    atoms.selected=True
