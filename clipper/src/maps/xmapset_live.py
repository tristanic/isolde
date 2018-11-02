
def map_potential_recommended_bsharp(resolution):
    '''
    Return a recommended sharpening/smoothing B-factor for a given map to
    optimise its use as an MDFF potential. For now this is a simple linear
    function of resolution, passing through zero at 2.5 Angstroms (smoothing
    below, sharpening above). Improved methods will be developed over time.
    '''
    # smooth by 30 A**2 at 1.5A res; sharpen by 30A**2 at 3.5A res; 0 at 2.5A res
    bsharp_base = 30
    return bsharp_base*resolution-2.5*bsharp_base

def viewing_recommended_bsharp(resolution):
    '''
    For viewing purposes it is also often useful to have a smoothed or
    sharpened visualisation of your map, but the optimal degree of sharpening
    for visualisation is not necessarily the same as that for MDFF.
    '''
    # smooth by 50 A**2 at 1.5A res; sharpen by 50 A**2 at 3.5A res, 0 at 2.5A res
    bsharp_base = 50
    return bsharp_base*resolution-2.5*bsharp_base

from .mapset_base import MapSet_Base

class XmapSet_Live(MapSet_Base):
    '''
    Handles creation, deletion, recalculation and visualisation of
    crystallographic maps based on the current model and a set of observed
    reflection data (FOBS, SIGFOBS).
    '''

    def __init__(self, manager, fsigf_name=None, bsharp_vals=[],
        exclude_free_reflections=False, fill_with_fcalc=False,
        live_update=True, show_r_factors=True):
        '''
        Prepare the C++ :class:`Xtal_mgr` object and create a set of
        crystallographic maps. The standard 2mFo-DFc and mFo-DFc maps will
        always be created, while any number of maps with different
        sharpening/smoothing can be created by providing a list of b-factors in
        bsharp_vals.
        Args:
            manager:
                The master :class:`Map_Mgr`
            fsigf_name:
                The label of the :class:`F_sigF_float` object containing the
                observed amplitudes
            bsharp_vals:
                For each value in this list, a 2mFo-DFc map will be generated
                with the given B_sharp value. A negative B_sharp yields a
                sharpened map, while a positive B_sharp gives smoothing. As a
                rough rule of thumb, a value of B_sharp=-100 tends to be
                reasonable for a 3.8 Angstrom map, and -50 for 3 Angstroms.
                For maps with resolutions better than ~2 Angstroms, it may be
                useful in some circumstances to apply smoothing. This can
                sometimes bring out weaker, low-resolution features hidden by
                high-resolution noise in standard maps.
            exclude_free_reflections:
                If True, observed amplitudes corresponding to the free set will
                not be used in generating the maps. The values used in their
                place will depend on the value of `fill_with_fcalc`. Note that
                this only affects maps generated for viewing - the MDFF potential
                map is always generated with free reflections excluded.
            fill_with_fcalc:
                If `exclude_free_reflections` is False this argument will be
                ignored. Otherwise, if `fill_with_fcalc` is True then the
                excluded amplitudes will be replaced by sigmaa-weighted F(calc).
                If False, the excluded amplitudes will be set to zero.
            display_radius:
                The radius (in Angstroms) of the display sphere used in
                live scrolling mode.
            live_update:
                If True, maps will be automatically recalculated whenever
                coordinates change
            show_r_factors:
                If True, print new R factors to the status bar on each map
                recalculation
        '''
        super().__init__(manager, 'Live crystallographic maps')

        from chimerax.core.triggerset import TriggerSet
        trig = self._triggers = TriggerSet()

        trigger_names = (
            'maps recalculated'
        )
        for t in trigger_names:
            trig.add_trigger(t)

        self._live_update = False
        self._recalc_needed = False
        self._model_changes_handler = None
        self._delayed_recalc_handler = None
        self._show_r_factors = show_r_factors

        from ..util import available_cores
        # The master C++ manager for handling all map operations
        from ..clipper_python.ext import Xtal_thread_mgr
        xm = self._xtal_mgr = Xtal_thread_mgr(crystal.hklinfo,
            crystal.mtzdata.free_flags.data, crystal.grid,
            crystal.mtzdata.experimental_data.datasets[fsigf_name].data,
            num_threads=available_cores())

        from .. import atom_list_from_sel
        xm.init(atom_list_from_sel(self.structure.atoms))
        # Only this map will actually be used as the MDFF potential
        self.add_live_xmap('MDFF potential', 0, is_difference_map=False,
            exclude_missing_reflections=True,
            exclude_free_reflections=True, fill_with_fcalc=True,
            display=False)

        self.add_live_xmap('2mFo-DFc', 0, is_difference_map=False,
            exclude_free_reflections=exclude_free_reflections,
            fill_with_fcalc=fill_with_fcalc,
            display=True)

        self.add_live_xmap('mFo-DFc', 0, is_difference_map=True,
            exclude_free_reflections=exclude_free_reflections,
            fill_with_fcalc=fill_with_fcalc,
            display=True)

        for b in bsharp_vals:
            if abs(b) < 5:
                continue
            elif b < 0:
                name_str = "2mFo-DFc_smooth_{:.0f}".format(-b)
            else:
                name_str = "2mFo-DFc_sharp_{:.0f}".format(b)
            self.add_live_xmap(name_str, b, is_difference_map=False,
                exclude_free_reflections=exclude_free_reflections,
                fill_with_fcalc=fill_with_fcalc
            )

        self.master_map_manager.rezone_needed()

    @property
    def xtal_mgr(self):
        '''
        The core manager object responsible for handling the reflection data and
        map recalculations.
        '''
        return self._xtal_mgr

    @property
    def show_r_factors(self):
        return self._show_r_factors

    @show_r_factors.setter
    def show_r_factors(self, flag):
        self._show_r_factors = flag

    @property
    def live_update(self):
        return self._live_update

    @live_update.setter
    def live_update(self, flag):
        if flag == self._live_update:
            return
        if flag:
            if self._model_changes_handler is None:
                self._model_changes_handler = self.model.triggers.add_handler(
                    'changes', self._model_changed_cb
                )
        else:
            if self._model_changes_handler is not None:
                self.model.triggers.remove_handler(
                    self._model_changes_handler
                )
                self._model_changes_handler = None
        self._live_update = flag

    @property
    def rfree(self):
        return self._xtal_mgr.rfree

    @property
    def rwork(self):
        return self._xtal_mgr.rwork

    @property
    def res(self):
        return self.hklinfo.resolution.limit

    @property
    def resolution(self):
        return self.hklinfo.resolution.limit

    @property
    def voxel_size(self):
        return self.cell.dim / self.grid.dim

    @property
    def voxel_size_frac(self):
        return 1/ self.grid.dim

    @property
    def unit_cell(self):
        return self.crystal_mgr.unit_cell

    @property
    def triggers(self):
        return self._triggers

    def add_live_xmap(self, name, b_sharp,
        is_difference_map=False,
        exclude_missing_reflections=False,
        exclude_free_reflections=True,
        fill_with_fcalc=True,
        color=None, style=None, contour=None, display=True):
        xm = self._xtal_mgr
        xm.add_xmap(name, b_sharp, is_difference_map=is_difference_map,
            exclude_missing_reflections=exclude_missing_reflections,
            exclude_free_reflections=exclude_free_reflections,
            fill_with_fcalc = fill_with_fcalc)
        #xmap = xm.get_xmap_ref(name)
        new_handler = self.add_xmap_handler(name, is_difference_map=is_difference_map,
            color=color, style=style, contour=contour)
        if display:
            new_handler.show()
            self.master_map_mgr.rezone_needed()
        else:
            new_handler.display = False

    def add_xmap_handler(self, name, is_difference_map = None,
                color = None, style = None, contour = None):
        '''
        Add a new XmapHandler_Live based on the given reflections and phases.
        Args:
            name:
                A unique string describing this map
            xmap:
                a Clipper :class:`Xmap_float` object.
            is_difference_map:
                Decides whether this map is to be treated as a difference
                map (with positive and negative contours) or a standard
                map with positive contours only. Leave as None to allow
                this to be determined automatically from the dataset,
                or set it to True or False to force an explicit choice.
            color:
                an array of [r, g, b, alpha] as integers in the range 0-255
                for a positive density map, or
                [[r,g,b,alpha],[r,g,b,alpha]] representing the colours of
                the positive and negative contours respectively if the
                map is a difference map
            style:
                one of 'mesh' or 'surface'
            contour:
                The value(s) (in sigma units) at which to contour the map
                on initial loading. For a standard map this should be a
                single value; for a difference map it should be
                [negative contour, positive contour]
        '''
        if is_difference_map and color is not None and len(color) != 2:
            err_string = '''
            ERROR: For a difference map you need to define colours for
            both positive and negative contours, as:
            [[r,g,b,a],[r,g,b,a]] in order [positive, negative].
            '''
            raise TypeError(err_string)
        new_handler = XmapHandler_Live(self.session, self, name,
            *self.box_params,
            is_difference_map = is_difference_map)
        if style is None:
            style = 'mesh'
        if color is None:
            if is_difference_map:
                color = self.DEFAULT_DIFF_MAP_COLORS
            elif style == 'mesh':
                color = [self.DEFAULT_MESH_MAP_COLOR]
            else:
                color = [self.DEFAULT_SOLID_MAP_COLOR]
        if contour is None:
            if is_difference_map:
                contour = self.STANDARD_DIFFERENCE_MAP_CONTOURS
            else:
                contour = self.STANDARD_LOW_CONTOUR
        elif not hasattr(contour, '__len__'):
                contour = numpy.array([contour])
        else:
            contour = numpy.array(contour)
        contour = contour * new_handler.sigma
        self.add([new_handler])
        new_handler.set_representation(style)
        new_handler.set_parameters(**{'cap_faces': False,
                                  'surface_levels': contour,
                                  'show_outline_box': False,
                                  'surface_colors': color,
                                  'square_mesh': False})
        # new_handler.update_surface()
        return new_handler


    _map_impacting_changes = set ((
        "aniso_u changed",
        "bfactor changed",
        "coord changed",
        "coordset changed",
        "occupancy changed",
    ))
    def _model_changed_cb(self, trigger_name, changes):
        if changes is not None:
            changes = set(changes[1].atom_reasons()).intersection(self._map_impacting_changes)
            if changes:
                self._recalc_needed = True
                if self._delayed_recalc_handler is None:
                    self._delayed_recalc_handler = self.session.triggers.add_handler(
                        'new frame', self._recalculate_maps_if_needed
                    )

    def _recalculate_maps_if_needed(self, *_):
        xm = self._xtal_mgr
        if self._recalc_needed and not xm.thread_running:
            self.recalculate_all_maps(self.model.atoms)
            if self._delayed_recalc_handler is not None:
                self.session.triggers.remove_handler(self._delayed_recalc_handler)
                self._delayed_recalc_handler = None

    def recalculate_all_maps(self, atoms):
        from .. import atom_list_from_sel
        from ..delayed_reaction import delayed_reaction
        xm = self._xtal_mgr
        delayed_reaction(self.session.triggers, 'new frame',
            xm.recalculate_all_maps, [atom_list_from_sel(atoms)],
            xm.ready,
            self._apply_new_maps, []
            )
        self._recalc_needed = False

    def _apply_new_maps(self):
        xm = self._xtal_mgr
        xm.apply_new_maps()
        if self.show_r_factors:
            self.session.logger.status(
                'R-work: {:0.4f}  Rfree: {:0.4f}'.format(
                    xm.rwork, xm.rfree
                )
            )
        self.triggers.activate_trigger('maps recalculated', None)

from .map_handler_base import XmapHandler_Base

class XmapHandler_Live(XmapHandler_Base):
    '''
    An XmapHandler_Live is in effect a resizable window into a periodic
    crystallographic map. The actual map data (a clipper Xmap object) is
    recalculated as needed from the combination of atomic coordinates and
    measured diffraction data, and filled into the XmapHandler_Live.data array
    as needed. Mothods are provided for live recalculation, tracking and filling
    a box around the centre of rotation, and static display of a given region.
    '''
    def __init__(self, mapset, name, origin, grid_origin, grid_dim,
        is_difference_map=False):
        '''
        Args:
            mapset:
                The XmapSet_Live object this belongs to
            name:
                A descriptive name for this map
            origin:
                The (x,y,z) coordinates of the bottom left corner of the
                volume.
            grid_origin:
                The (u,v,w) integer grid coordinates corresponding to
                origin.
            grid_dim:
                The shape of the box in (u,v,w) grid coordinates.
            is_difference_map:
                Is this a difference map?
        '''
        super().__init__(mapset, name, origin, grid_origin, grid_dim,
            is_difference_map=is_difference_map)
        self._mapset_handlers.append(
            self.mapset.triggers.add_handler(
                'maps recalculated', self._map_recalc_cb
                )
        )

    @property
    def xtal_mgr(self):
        return self.mapset.xtal_mgr

    @property
    def xmap(self):
        return self.xtal_mgr.get_xmap_ref(self.name)

    @property
    def stats(self):
        return self.manager._xtal_mgr.get_map_stats(self.name)

    def _map_recalc_cb(self, name, *_):
        for s in self.surfaces:
            s._use_thread=True
        self._fill_volume_data(self._data_fill_target, self.box_params[1])
        self.data.values_changed()
