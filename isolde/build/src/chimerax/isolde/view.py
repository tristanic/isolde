
def focus_on_selection(session, view, atoms):
    v = view
    pad = 5.0
    bounds = atoms.scene_bounds
    v.center_of_rotation = center = bounds.center()
    bounds.xyz_min = bounds.xyz_min - pad
    bounds.xyz_max = bounds.xyz_max + pad
    radius = bounds.radius() + pad
    v.view_all(bounds)
    cam = v.camera
    vd = cam.view_direction()
    cp = v.clip_planes
    cp.set_clip_position('near', center - radius*vd, cam)
    cp.set_clip_position('far', center + radius*vd, cam)
    session.selection.clear()
    atoms.selected=True    
