# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 25-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll

def focus_on_coord(session, center, radius = 5.0, clip=True, rescale=True):
    '''
    Focus the main view on a coordinate, maintaining the current center of
    rotation method and optionally updating the near and far clipping planes.

    Args:
        * session:
            - the top-level ChimeraX `session` instance
        * center:
            - (x,y,z) coordinates as a Numpy array
        * radius (default=5.0):
            - the zoom level, camera position and (optionally) clipping planes
              will be adjusted to ensure that everything within at least this
              distance from the center is shown.
        * clip (default=True):
            - if `True`, updates the near and far clipping planes.
    '''
    v = session.main_view
    import numpy
    center = numpy.array(center)
    from chimerax.geometry import Bounds
    xyz_min = center-radius
    xyz_max = center+radius
    bounds = Bounds(xyz_min, xyz_max)
    cofr_method = v.center_of_rotation_method
    cam = v.camera
    vd = cam.view_direction()
    if rescale:
        v.view_all(bounds)
    else:
        # Keep the view direction and scale, and just move the camera
        camera_to_cofr = cam.position.origin() - v.center_of_rotation
        from chimerax.geometry import Place
        cam.position = Place(axes=cam.position.axes(), origin=center+camera_to_cofr)

    v.center_of_rotation = center
    v.center_of_rotation_method = cofr_method
    if clip:
        cp = v.clip_planes
        cp.set_clip_position('near', center - radius*vd, v)
        cp.set_clip_position('far', center + radius*vd, v)


def focus_on_selection(session, atoms, pad=1.0, clip = True):
    '''
    Focus the main view on a selecton of atoms, maintaining the current center
    of rotation method and optionally updating the near and far clipping
    planes.

    Args:
        * session:
            - the top-level ChimeraX `session` instance
        * atoms:
            - a :class:`chimerax.Atoms` instance
        * pad (default=1.0):
            - the zoom level, camera position and (optionally) clipping planes
              will be adjusted to ensure that everything within at least this
              distance from the selection is shown.
        * clip (default=True):
            - if `True`, updates the near and far clipping planes.
    '''
    bounds = atoms.scene_bounds
    focus_on_coord(session, bounds.center(), bounds.radius()+pad)
    session.selection.clear()
    atoms.selected=True
    atoms.intra_bonds.selected = True
