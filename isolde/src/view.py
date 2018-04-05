# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.



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
