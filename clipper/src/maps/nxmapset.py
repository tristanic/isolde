# @Author: Tristan Croll <tic20>
# @Date:   15-Jan-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 08-May-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



import numpy
from .mapset_base import MapSet_Base

class NXmapSet(MapSet_Base):
    '''
    Manages real-space maps. The most important difference between these and
    crystallographic maps is that there is no guarantee that two maps will have
    the same grid (i.e. voxel size and angles).
    '''
    def __init__(self, manager, name):
        super().__init__(manager, name)
        manager.add([self])

    def add_nxmap_handler_from_file(self, filename, is_difference_map=False,
        color=None, style=None, contour=None):
        from chimerax.map.data import open_file
        grid_data = open_file(filename)[0]
        from chimerax.map.volume import volume_from_grid_data
        return self.add_nxmap_handler_from_volume(
            volume_from_grid_data(grid_data, self.session),
            is_difference_map=is_difference_map,
            color=color, style=style, contour=contour
        )


    def add_nxmap_handler_from_volume(self, volume,
        is_difference_map=False,
        color=None, style=None, contour=None):
        h = NXmapHandler(self, volume)
        if self.spotlight_mode:
            corners = _find_box_corners(self.box_center, self.display_radius,
                h.data.xyz_to_ijk_transform)
            h.expand_to_cover_coords(corners, 15)
        self.add([h])
        self.set_nxmap_display_style(h)
        self.crystal_mgr.normalize_scene_positions()
        self.master_map_mgr._reapply_zone()
        return h

    def set_nxmap_display_style(self, nxmap_handler, is_difference_map=False,
        color=None, style=None, contour=None):
        if style is None:
            style='mesh'
        if is_difference_map and color is not None and len(color) != 2:
            err_string = '''
            ERROR: For a difference map you need to define colours for
            both positive and negative contours, as:
            [[r,g,b,a],[r,g,b,a]] in order [positive, negative].
            '''
            raise TypeError(err_string)
        if color is None:
            if is_difference_map:
                color = self.DEFAULT_DIFF_MAP_COLORS
            elif style == 'mesh':
                color = [self.DEFAULT_MESH_MAP_COLOR]
            else:
                color = [self.DEFAULT_SOLID_MAP_COLOR]
        if contour is None:
            from ..util import guess_suitable_contour
            if is_difference_map:
                pcontour = guess_suitable_contour(nxmap_handler, self.structure,
                    atom_radius_scale = 0.8)
                contour = numpy.array([-pcontour, pcontour])
            else:
                contour = numpy.array([guess_suitable_contour(nxmap_handler, self.structure)])
        nxmap_handler.set_parameters(**{'cap_faces': False,
                                  'surface_levels': contour,
                                  'show_outline_box': False,
                                  'surface_colors': color,
                                  'square_mesh': False})
        nxmap_handler.set_display_style(style)



    def expand_to_cover_coords(self, coords, padding):
        for v in self:
            v.expand_to_cover_coords(coords, padding)



from .map_handler_base import MapHandler_Base
class NXmapHandler(MapHandler_Base):
    '''
    Real-space equivalent to XmapHandler_Static. Doesn't actually use any of
    the clipper engine, but provides a unified interface.
    '''
    def __init__(self, mapset, volume, is_difference_map=False):
        '''
        Takes ownership of the data from an existing Volume object.
        The input volume will be closed.
        '''
        self._original_volume = volume
        super().__init__(mapset, volume.name, volume.data,
            is_difference_map=is_difference_map)
        if volume in self.session.models.list():
            self.session.models.remove([volume])


    def _box_changed_cb(self, name, params):
        self.update_mask()

    def _box_moved_cb(self, name, params):
        self.update_mask()

    @property
    def stats(self):
        return self.mean_sd_rms()

    @property
    def sigma(self):
        return self.mean_sd_rms()[1]

    def update_mask(self):
        if not self.display:
            return
        corners = _find_box_corners(self.center, self.display_radius, self.data.xyz_to_ijk_transform)
        self.new_region(ijk_min=corners[0], ijk_max=corners[1], ijk_step=[1,1,1])

    def expand_to_cover_coords(self, coords, padding):
        self.new_region(*self.bounding_region(coords, padding=padding, step=[1,1,1]))

    def delete(self):
        self._original_volume.delete()
        super().delete()

_corners = numpy.array([(x,y,z) for x in (-1,1) for y in (-1,1) for z in (-1,1)])
def _find_box_corners(center, radius, xyz_to_ijk_transform):
    corners = xyz_to_ijk_transform*(center+radius*_corners)
    return (numpy.min(corners, axis=0), numpy.max(corners, axis=0))
