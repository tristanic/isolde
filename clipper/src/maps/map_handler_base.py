# @Author: Tristan Croll <tic20>
# @Date:   15-Jan-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 20-May-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



from chimerax.map import Volume
import numpy

class MapHandler_Base(Volume):
    '''
    Base class for all Clipper map objects for use in ChimeraX.
    '''
    pickable=False
    def __init__(self, mapset, name, data, is_difference_map=False):
        session = mapset.session
        data.change_callbacks.clear()
        super().__init__(session, data)
        self.name = name
        ms = self._mapset = mapset

        self._is_difference_map = is_difference_map
        #self.initialize_thresholds()

        self.show()
        mh = self._mgr_handlers = []
        mh.append(
            (ms,
            ms.triggers.add_handler('map box changed', self._box_changed_cb) )
        )
        mh.append(
            (ms,
            ms.triggers.add_handler('map box moved', self._box_moved_cb) )
        )

        # Cache the mean, sigma and RMS, since it's used all the time in contouring
        self._mean_sd_rms = self._full_mean_sd_rms()
        self.data.add_change_callback(self._update_mean_sd_rms_cb)

    def _full_mean_sd_rms(self):
        from chimerax.map.volume import mean_sd_rms
        return mean_sd_rms(self.data.matrix())

    def mean_sd_rms(self):
        return self._mean_sd_rms

    def _update_mean_sd_rms_cb(self, reason):
        if reason == 'values changed':
            self._mean_sd_rms = self._full_mean_sd_rms()

    @property
    def mapset(self):
        return self._mapset

    @property
    def manager(self):
        return self.mapset.master_map_mgr

    @property
    def crystal_mgr(self):
        return self.manager.crystal_mgr

    @property
    def center(self):
        return self.manager.box_center

    @property
    def display_radius(self):
        return self.manager.spotlight_radius

    @property
    def is_difference_map(self):
        return self._is_difference_map

    def _box_changed_cb(self, name, params):
        pass

    def _box_moved_cb(self, name, params):
        pass

    def delete(self):
        for (mgr, h) in self._mgr_handlers:
            try:
                mgr.triggers.remove_handler(h)
            except:
                continue
        super().delete()

    def add_surface(self, level, rgba=None):
        ses = self.session
        s = FastVolumeSurface(self, level, rgba)
        self._surfaces.append(s)
        if self.id is None:
            self.add([s])
        else:
            ses.models.add([s], parent=self)
        return s


from chimerax.map.volume import VolumeSurface
class FastVolumeSurface(VolumeSurface):
    '''
    Threaded implementation of ChimeraX VolumeSurface, with threading at the
    C++ level. Should become obsolete once something similar is implemented in
    ChimeraX itself.
    '''
    def __init__(self, volume, level, rgba=(1.0, 1.0, 1.0, 1.0)):
        if rgba is None:
            rgba = (1.0,1.0,1.0,1.0)
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

    def _use_fast_thread_result(self, rendering_options):
        sct = self._surf_calc_thread
        if sct is not None:
            va, ta, na = sct.get_result()
            va, na, ta, hidden_edges = self._postprocess(va, na, ta, self.volume.rendering_options, self.level)
            self._set_surface(va, na, ta, hidden_edges)
            self._set_appearance(rendering_options)
            self._surf_calc_thread = None
            if self._update_needed:
                # surface properties were changed while the thread was working
                self.update_surface(rendering_options)
            self._update_needed = False
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def update_surface(self, rendering_options):
        sct = self._surf_calc_thread
        if sct is not None: # and not sct.ready():
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
            sct.start_compute, (numpy.asfortranarray(v.matrix()), level, det, vertex_transform.matrix, normal_transform.matrix, False, True),
            sct.ready,
            self._use_fast_thread_result, (rendering_options,))

class XmapHandler_Base(MapHandler_Base):
    '''
    Base class for XmapHandler_Static and XmapHandler_Live
    '''
    def __init__(self, mapset, name,
            is_difference_map = False):
        self._mapset = mapset
        self._name = name
        bp = mapset.box_params

        darray = self._generate_and_fill_data_array(bp.origin_xyz, bp.origin_grid, bp.dim)
        # darray = self._generate_data_array(*mapset.box_params)
        super().__init__(mapset, name, darray,
            is_difference_map=is_difference_map)


    @property
    def box_params(self):
        return self.mapset.box_params

    @property
    def hklinfo(self):
        return self.mapset.hklinfo

    @property
    def spacegroup(self):
        return self.mapset.spacegroup

    @property
    def cell(self):
        return self.mapset.cell

    @property
    def res(self):
        return self.hklinfo.resolution

    @property
    def grid(self):
        return self.mapset.grid

    @property
    def voxel_size(self):
        return self.cell.dim / self.grid.dim

    @property
    def voxel_size_frac(self):
        return 1/ self.grid.dim

    @property
    def unit_cell(self):
        return self.mapset.unit_cell

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

    def mean_sd_rms(self):
        '''
        Overrides the standard Volume method to give the overall values
        from the Clipper object.
        '''
        s = self.stats
        return (s[0], s[1], s[1])

    @property
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

    def _box_changed_cb(self, *_):
        self._swap_volume_data(self.box_params)
        self._use_thread = True
        self.data.values_changed()

    def _box_moved_cb(self, *_):
        bp = self.box_params
        self.data.set_origin(bp.origin_xyz)
        self._fill_volume_data(self._data_fill_target, bp.origin_grid)
        for s in self.surfaces:
            s._use_thread=True
        self.data.values_changed()

    def _generate_data_array(self, origin, grid_origin, dim):
        data = self._data_fill_target = numpy.empty(dim, numpy.float32)
        order = numpy.array([2,1,0], int)
        from chimerax.map.data import ArrayGridData
        darray = ArrayGridData(data.transpose(), origin = origin,
            step = self.voxel_size, cell_angles = self.cell.angles_deg)
        return darray

    def _fill_volume_data(self, target, start_grid_coor):
        from .. import Coord_grid
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
