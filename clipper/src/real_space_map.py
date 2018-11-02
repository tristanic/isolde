# @Author: Tristan Croll
# @Date:   19-Mar-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   Tristan Croll
# @Last modified time: 18-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



from chimerax.map import Volume
import numpy

class NXmapHandler(Volume):
    '''
    Real-space equivalent to crystal.XmapHandler. Doesn't actually use any of
    the clipper engine, but provides a unified interface.
    '''
    def __init__(self, session, manager, volume, is_difference_map = False):
        '''
        Takes ownership of the data from an existing Volume object.
        The input volume will be closed.
        '''
        self.session = session
        self.manager = manager
        super().__init__(volume.data, session)
        self.name = volume.name

        self._is_difference_map = is_difference_map
        self.initialize_thresholds()

        # self._needs_update = True
        self.show()
        self._box_changed_cb_handler = self.manager.triggers.add_handler(
            'map box changed', self._box_changed_cb
        )
        self._box_moved_cb_handler = self.manager.triggers.add_handler(
            'map box moved', self._box_moved_cb
        )
        session.models.remove([volume])

    @property
    def center(self):
        return self.manager._box_center

    @property
    def radius(self):
        return self.manager.display_radius

    @property
    def is_difference_map(self):
        return self._is_difference_map

    def show(self, *args, **kwargs):
        self.update_mask()
        super().show(*args, **kwargs)

    def _box_changed_cb(self, name, params):
        self._needs_update = True
        if not self.display:
            return
        self.update_mask()

    def _box_moved_cb(self, name, params):
        #self.box_params = params
        if not self.display:
            return
        self.update_mask()

    def update_mask(self):
        if not self.display:
            return
        corners = _find_box_corners(self.center, self.radius, self.data.xyz_to_ijk_transform)
        self.new_region(ijk_min=corners[0], ijk_max=corners[1], ijk_step=[1,1,1])

    def delete(self):
        bh = self._box_shape_changed_cb_handler
        if bh is not None:
            self.manager.triggers.remove_handler(bh)
            self._box_changed_cb_handler = None
        bm = self._box_moved_cb_handler
        if bm is not None:
            self.manager.triggers.remove_handler(bm)
            self._box_moved_cb_handler = None
        super().delete()

_corners = numpy.array([(x,y,z) for x in (-1,1) for y in (-1,1) for z in (-1,1)])
def _find_box_corners(center, radius, xyz_to_ijk_transform):
    corners = xyz_to_ijk_transform*(center+radius*_corners)
    return (numpy.min(corners, axis=0), numpy.max(corners, axis=0))
