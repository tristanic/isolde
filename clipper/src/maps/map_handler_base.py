
from chimerax.map import Volume
import numpy

class MapHandler_Base(Volume):
    '''
    Base class for all Clipper map objects for use in ChimeraX.
    '''
    def __init__(self, mapset, name, data, is_difference_map=False):
        session = mapset.session
        super().__init__(data, session)
        self.name = name
        self._mapset = mapset

        self._is_difference_map = is_difference_map
        self.initialize_thresholds()

        self.show()
        mh = self._mgr_handlers = []
        mh.append(self.manager.triggers.add_handler(
            'map box changed', self._box_changed_cb)
        )
        mh.append(self.manager.triggers.add_handler(
            'map box moved', self._box_moved_cb
        )

    @property
    def mapset(self):
        return self._mapset

    @property
    def manager(self):
        return self.mapset.master_map_mgr

    @property
    def center(self):
        return self.manager.box_center

    @property
    def display_radius(self):
        return self.manager.display_radius

    @property
    def box_params(self):
        return self.manager.box_params

    @property
    def is_difference_map(self):
        return self._is_difference_map

    def _box_changed_cb(self, name, params):
        pass

    def _box_moved_cb(self, name, params):
        pass

    def delete(self):
        for h in self._mgr_handlers:
            try:
                self.manager.triggers.remove_handler(h)
            except:
                continue
        super().delete()

from chimerax.map.volume import VolumeSurface
class FastVolumeSurface(VolumeSurface):
    '''
    Threaded implementation of ChimeraX VolumeSurface, with threading at the
    C++ level. Should become obsolete once something similar is implemented in
    ChimeraX itself.
    '''
    def __init__(self, volume, level, rgba=(1.0, 1.0, 1.0, 1.0)):
        super().__init__(volume, level, rgba)
        self._update_needed = False

    def _postprocess(self, varray, narray, tarray, rendering_options, level):
        ro = rendering_options
        if ro.flip_normals and level < 0:
          from chimerax.surface import invert_vertex_normals
          invert_vertex_normals(narray, tarray)

        # reverse_triangle_vertex_order done in thread

        if ro.subdivide_surface:
            from chimerax.surface import subdivide_triangles
            for i in range(ro.subdivision_levels):
                varray, tarray, narray = subdivide_triangles(varray, tarray, narray)

        if ro.square_mesh:
            from numpy import empty, uint8
            hidden_edges = empty((len(tarray),), uint8)
            from chimerax.map import _map
            _map.principle_plane_edges(varray, tarray, hidden_edges)
        else:
            hidden_edges = None

        if ro.surface_smoothing:
          sf, si = ro.smoothing_factor, ro.smoothing_iterations
          from chimerax.surface import smooth_vertex_positions
          smooth_vertex_positions(varray, tarray, sf, si)
          smooth_vertex_positions(narray, tarray, sf, si)

        # Transforms and normalization done in thread
        return varray, narray, tarray, hidden_edges

    def _use_fast_thread_result(self, show_mesh, rendering_options):
        sct = self._surf_calc_thread
        if sct is not None:
            va, ta, na = sct.get_result()
            va, na, ta, hidden_edges = self._postprocess(va, na, ta, self.volume.rendering_options, self.level)
            self._set_surface(va, na, ta, hidden_edges)
            self._set_appearance(show_mesh, rendering_options)
            self._surf_calc_thread = None
            if self._update_needed:
                # surface properties were changed while the thread was working
                self.update_surface(show_mesh, rendering_options)
            self._update_needed = False
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def update_surface(self, show_mesh, rendering_options):
        sct = self._surf_calc_thread
        if sct is not None and not sct.ready():
            self._update_needed = True
            return

        v = self.volume
        level = self.level
        vertex_transform = v.matrix_indices_to_xyz_transform()
        normal_transform = vertex_transform.zero_translation()
        #normal_transform = vertex_transform.inverse().transpose().zero_translation()
        det = vertex_transform.determinant()

        from ..delayed_reaction import delayed_reaction
        from ..contour_thread import Contour_Thread_Mgr
        sct = self._surf_calc_thread = Contour_Thread_Mgr()
        delayed_reaction(self.volume.session.triggers, 'new frame',
            sct.start_compute, (v.matrix(), level, det, vertex_transform.matrix, normal_transform.matrix, False, True),
            sct.ready,
            self._use_fast_thread_result, (show_mesh, rendering_options))



class XmapHandler_Base(MapHandler_Base):
    '''
    Base class for XmapHandler_Static and XmapHandler_Live
    '''
    def __init__(self, mapset, name, origin, grid_origin, grid_dim,
            is_difference_map = False):
        self.box_params = (origin, grid_origin, dim)
        darray = self._generate_data_array(origin, grid_origin, dim)
        super().__init__(mapset, name, data)
        self._mapset_handlers = []


    @property
    def hklinfo(self):
        return self.manager.hklinfo

    @property
    def spacegroup(self):
        return self.manager.spacegroup

    @property
    def cell(self):
        return self.manager.cell

    @property
    def res(self):
        return self.hklinfo.resolution

    @property
    def grid(self):
        return self.manager.grid

    @property
    def voxel_size(self):
        return self.cell.dim / self.grid.dim

    @property
    def voxel_size_frac(self):
        return 1/ self.grid.dim

    @property
    def unit_cell(self):
        return self.manager.unit_cell

    @property
    def xmap(self):
        raise NotImplementedError('This property must be overridden in the '
            'derived class!')

    @property
    def stats(self):
        '''
        Should return a 2-tuple of (mean, sigma) for the map
        '''
        raise NotImplementedError('This property must be overridden in the '
            'derived class!')

    def delete(self):
        for h in self._mapset_handlers:
            try:
                self.mapset.triggers.remove_handler(h)
            except:
                continue
        super().delete()

    def mean_sd_rms(self):
        '''
        Overrides the standard Volume method to give the overall values
        from the Clipper object.
        '''
        s = self.stats
        return (s[0], s[1], s[1])

    def sigma(self):
        return self.stats[1]

    def add_surface(self, level, rgba=None):
        ses = self.session
        s = FastVolumeSurface(self, level, rgba)
        self._surfaces.append(s)
        if self.id is None:
            self.add([s])
        else:
            ses.models.add([s], parent=self)
        return s

    def force_refill(self):
        self._box_changed_cb(None, (self.box_params, True))

    def show(self, *args, **kwargs):
        if self._needs_update:
            self._swap_volume_data(self.box_params)
        else:
            # Just set the origin and fill the box with the data for
            # the current location
            origin, grid_origin, ignore = self.box_params
            self._fill_volume_data(self._data_fill_target, grid_origin)
        super().show(*args, **kwargs)

    def _box_changed_cb(self, name, params):
        box_params, force_fill = params
        self.box_params = box_params
        self._needs_update = True
        if not self.display and not force_fill:
            # No sense in wasting cycles on this if the volume is hidden.
            # We'll just store the params and apply them when we show the
            # volume.
            # NOTE: this means we need to over-ride show() to ensure
            # it's updated before re-displaying.
            return
        self._swap_volume_data(box_params)
        self._use_thread = True
        self.data.values_changed()

    def _box_moved_cb(self, name, params):
        self.box_params = params
        if not self.display:
            return
        self.data.set_origin(params[0])
        self._fill_volume_data(self._data_fill_target, params[1])
        for s in self.surfaces:
            s._use_thread=True
        self.data.values_changed()

    def _generate_data_array(self, origin, grid_origin, dim):
        data = self._data_fill_target = numpy.empty(dim, numpy.float32)
        order = numpy.array([2,1,0], int)
        darray = Array_Grid_Data(data.transpose(), origin = origin,
            step = self.voxel_size, cell_angles = self.cell.angles_deg)
        return darray

    def _fill_volume_data(self, target, start_grid_coor):
        from .clipper_python import Coord_grid
        xmap = self.xmap
        xmap.export_section_numpy(Coord_grid(start_grid_coor), target)

    def _generate_and_fill_data_array(self, origin, grid_origin, dim):
        darray = self._generate_data_array(origin, grid_origin, dim)
        self._fill_volume_data(self._data_fill_target, grid_origin)
        return darray

    def _swap_volume_data(self, params):
        '''
        Replace this Volume's data array with one of a new shape/size
        Args:
            params:
                A tuple of (new_origin, new_grid_origin, new_dim)
        '''
        new_origin, new_grid_origin, new_dim = params
        darray = self._generate_and_fill_data_array(new_origin, new_grid_origin, new_dim)
        self._box_dimensions = new_dim
        self.replace_data(darray)
        self.new_region(ijk_min=(0,0,0), ijk_max=darray.size, ijk_step=(1,1,1), adjust_step=False)
        self._needs_update=False
