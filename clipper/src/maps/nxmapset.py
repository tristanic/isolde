import numpy
from .mapset_base import MapSet_Base

class NXmapSet(MapSet_Base):
    '''
    Manages real-space maps. The most important difference between these and
    crystallographic maps is that there is no guarantee that two maps will have
    the same grid (i.e. voxel size and angles).
    '''

    def add_nxmap_handler(self, volume,
        is_difference_map=None,
        color=None, style=None, contour=None):
        h = NXmapHandler(self, volume)
        if self.spotlight_mode:
            corners = _find_box_corners(self.box_center, self.display_radius,
                h.data.xyz_to_ijk_transform)
            h.expand_to_cover_coords(corners, 15)

    def expand_to_cover_coords(self, coords, padding):
        for v in self:
            v.expand_to_cover_coords(self, coords, padding)



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
        super().__init__(mapset, volume.name, volume.data,
            is_difference_map=is_difference_map)

    def _xtal_lattice_to_map_lattice_transform(self):
        xtal_step = self.crystal_mgr.voxel_size
        xtal_angles = self.cell.angles
        my_step = self.data.step
        my_angles = numpy.radians(self.data.cell_angles)


    def _box_changed_cb(self, name, params):
        self._needs_update = True
        if not self.display:
            return
        self.update_mask()

    def _box_moved_cb(self, name, params):
        if not self.display:
            return
        self.update_mask()

    def update_mask(self):
        if not self.display:
            return
        corners = _find_box_corners(self.center, self.display_radius, self.data.xyz_to_ijk_transform)
        self.new_region(ijk_min=corners[0], ijk_max=corners[1], ijk_step=[1,1,1])

    def expand_to_cover_coords(self, coords, padding):
        self.new_region(*self.bounding_region(coords, padding=padding, step=[1,1,1]))

_corners = numpy.array([(x,y,z) for x in (-1,1) for y in (-1,1) for z in (-1,1)])
def _find_box_corners(center, radius, xyz_to_ijk_transform):
    corners = xyz_to_ijk_transform*(center+radius*_corners)
    return (numpy.min(corners, axis=0), numpy.max(corners, axis=0))
