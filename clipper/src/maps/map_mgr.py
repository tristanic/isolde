from chimerax.core.models import Model

class Surface_Zone:
    '''
    Add this as a property to a Volume object to provide it with the
    necessary information to update its triangle mask after re-contouring.
    '''
    def __init__(self, distance, atoms = None, coords = None):
        '''
        Args:
          distance (float in Angstroms):
            distance from points to which the map will be masked
          atoms:
            Atoms to mask to (coordinates will be updated upon re-masking)
          coords:
            (x,y,z) coordinates to mask to (will not be updated upon
            re-masking).

          Set both atoms and coords to None to disable automatic re-masking.
        '''
        self.update(distance, atoms, coords)

    def update(self, distance, atoms = None, coords = None):
        self.distance = distance
        self.atoms = atoms
        self.coords = coords

    @property
    def all_coords(self):
        if self.atoms is not None:
            import numpy
            if self.coords is not None:
                return numpy.concatenate(self.atoms.coords, self.coords)
            return self.atoms.coords
        return self.coords

def surface_zones(models, points, distance):
    from chimerax.surface import zone
    for m in models:
        for s in m.surfaces:
            spoints = s.position.inverse(is_orthonormal=True) * points
            zone.surface_zone(s, spoints, distance, auto_update=True)


class Map_Mgr(Model):
    '''
    Top-level manager for all maps associated with a model.
    '''
    def __init__(self, crystal_manager, spotlight_radius=12):
        cm = self._mgr = crystal_manager
        super().__init__('Map Manager', cm.session)
        self._live_xmapsets = []
        self._static_xmapsets = []
        self._nxmapsets = []

        from chimerax.core.triggerset import TriggerSet
        trig = self._triggers = TriggerSet()

        trigger_names = (
            # Deprecated
            'map box changed',
            # Ask each MapSet to expand its volumes to cover an arbitrary set
            # of coordinates as efficiently as possible.
            'cover coords',
            # Change the radius of the "spotlight" sphere. It is up to each
            # MapSet to work out how to accommodate it
            'spotlight changed',
            'spotlight moved',    # Just changed the centre of the box
        )
        for t in trigger_names:
            trig.add_trigger(t)

        # Handler for live box update
        self._box_update_handler = None
        # Object storing the parameters required for masking (used after
        # adjusting contours)
        self._surface_zone = Surface_Zone(spotlight_radius, None, None)
        # Is the map box moving with the centre of rotation?

        self._spotlight_center = None

        # Radius of the sphere in which the map will be displayed when
        # in live-scrolling mode
        self.spotlight_radius = spotlight_radius

        if self.spotlight_mode:
            self._start_spotlight_mode()

        self.display=False
        self._rezone_pending = False
        # Apply the surface mask

        mh = self._mgr_handlers = []
        mh.append((cm, cm.triggers.add_handler('mode changed',
            self._spotlight_mode_changed_cb)))


        self.session.triggers.add_handler('frame drawn', self._first_init_cb)
        cm.add([self])

    @property
    def xmapsets(self):
        '''
        Sets of crystallographic maps associated with this model. Each XmapSet
        handles the set of maps derived from a single crystallographic dataset.
        '''
        from .xmapset import XmapSet
        return [m for m in self.child_models() if isinstance(m, XmapSet)]

    @property
    def nxmapset(self):
        '''
        Handler for all real-space (non-crystallographic) maps associated with
        this model.
        '''
        from .nxmapset import NXmapSet
        for m in self.child_models():
            if isinstance(m, NXmapSet):
                return m
        return NXmapSet(self, 'Non-crystallographic maps')

    @property
    def all_xtal_maps(self):
        from .map_handler_base import XmapHandler_Base
        return [m for m in self.all_models() if isinstance(m, XmapHandler_Base)]

    @property
    def all_non_xtal_maps(self):
        from .nxmapset import NXmapHandler
        return [m for m in self.all_models() if isinstance(m, NXmapHandler)]

    @property
    def all_maps(self):
        from chimerax.map import Volume
        return [m for m in self.all_models() if isinstance(m, Volume)]

    @property
    def crystal_mgr(self):
        return self._mgr

    @property
    def box_center(self):
        return self.crystal_mgr.spotlight_center

    @property
    def structure(self):
        return self.crystal_mgr.structure

    @property
    def triggers(self):
        return self._triggers

    @property
    def spacegroup(self):
        return self.crystal_mgr.spacegroup

    @property
    def cell(self):
        return self.crystal_mgr.cell

    @property
    def spotlight_radius(self):
        '''Get/set the radius (in Angstroms) of the live map display sphere.'''
        return self._spotlight_radius

    @spotlight_radius.setter
    def spotlight_radius(self, radius):
        import numpy
        self._spotlight_radius = radius
        center = self.crystal_mgr.spotlight_center
        self.triggers.activate_trigger('spotlight changed',
            (center, radius)
        )
        self._surface_zone.update(radius, coords = numpy.array([center]))
        self._reapply_zone()

    @property
    def box_params(self):
        return (self._box_corner_xyz, self._box_corner_grid, self._box_dimensions)


    @property
    def spotlight_mode(self):
        '''
        Is live map scrolling turned on? Can only be changed via the master
        symmetry manager.
        '''
        return self.crystal_mgr.spotlight_mode

    @spotlight_mode.setter
    def spotlight_mode(self, switch):
        raise NotImplementedError(
            'Mode can only be changed via the master symmetry manager!')

    @property
    def spotlight_center(self):
        '''
        Current (x,y,z) position of the centre of the "spotlight". Read-only.
        '''
        return self.crystal_mgr.spotlight_center

    @property
    def last_covered_selection(self):
        '''
        Last set of coordinates covered by
        `self.crystal_mgr.isolate_and_cover_selection()`. Read-only.
        '''
        return self.crystal_mgr.last_covered_selection

    @property
    def display(self):
        return super().display

    @display.setter
    def display(self, switch):
        if switch:
            if self.spotlight_mode:
                self._start_spotlight_mode()
        Model.display.fset(self, switch)

    def add_xmapset_from_mtz(self, mtzfile, oversampling_rate=1.5):
        from ..clipper_mtz import ReflectionDataContainer
        mtzdata = ReflectionDataContainer(self.session, mtzfile,
            shannon_rate = oversampling_rate)
        cm = self.crystal_mgr
        if not cm.has_symmetry:
            self.session.logger.info('(CLIPPER) NOTE: No symmetry information found '
                'in model. Using symmetry from MTZ file.')
            cm = self.crystal_mgr
            cm.cell = mtzdata.cell
            cm.spacegroup = mtzdata.spacegroup
            cm.grid = mtzdata.grid_sampling

            cm.has_symmetry = True
        elif not self.symmetry_matches(mtzdata):
            raise RuntimeError('Symmetry info from MTZ file does not match '
                'symmetry info from model!')
        from .xmapset import XmapSet
        return XmapSet(self, mtzdata)

    def symmetry_matches(self, xtal_data):
        return (
            xtal_data.cell.equals(self.cell, 1.0)
            and xtal_data.spacegroup.spacegroup_number == self.spacegroup.spacegroup_number
        )


    def _spotlight_mode_changed_cb(self, *_):
        if self.spotlight_mode:
            self._start_spotlight_mode()
        else:
            self._stop_spotlight_mode()

    def _start_spotlight_mode(self):
        self.triggers.activate_trigger('spotlight changed',
            (self.spotlight_center, self.spotlight_radius)
        )
        if self._box_update_handler is None:
            self._box_update_handler = self.crystal_mgr.triggers.add_handler(
                'spotlight moved', self.update_spotlight)
            self.update_spotlight(None, self.spotlight_center)
        from chimerax.core.geometry import Places
        self.positions = Places()

    def _stop_spotlight_mode(self):
        if self._box_update_handler is not None:
            self.crystal_mgr.triggers.remove_handler(self._box_update_handler)
            self._box_update_handler = None

    def cover_box(self, minmax):
        '''
        Set the map box to fill a volume encompassed by the provided minimum
        and maximum xyz coordinates. Automatically turns off live scrolling.
        '''
        self.spotlight_mode = False
        xyz_min, xyz_max = minmax
        center = (xyz_min + xyz_max)/2
        self.triggers.activate_trigger('map box changed',
            (center, xyz_min, xyz_max))

    def cover_coords(self, coords, mask_radius=3, extra_padding=0):
        '''
        Expand all maps to a region large enough to cover the coords, plus
        mask_radius+extra_padding in every direction.
        '''
        self.triggers.activate_trigger('cover coords',
            (coords, mask_radius+extra_padding))
        self._surface_zone.update(mask_radius, coords=coords)
        self._reapply_zone()

    def update_spotlight(self, trigger_name, new_center):
        '''
        Update the position of the "spotlight" to surround the current centre of
        rotation. If this manager is not currently displayed, then filling the
        volumes with data around the new position will be deferred unless
        force is set to True.
        '''
        if self.spotlight_mode:
            import numpy
            self.triggers.activate_trigger('spotlight moved',
                new_center)
            self._surface_zone.update(self.spotlight_radius, coords = numpy.array([new_center]))
            self._reapply_zone()

    def rezone_needed(self):
        if not self._rezone_pending:
            self.session.triggers.add_handler('new frame', self._rezone_once_cb)


    # Callbacks

    def _first_init_cb(self, *_):
        self.display = True
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER


    def _rezone_once_cb(self, *_):
        self._reapply_zone()
        self._rezone_pending=False
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def _reapply_zone(self):
        '''
        Reapply any surface zone applied to the volume after changing box
        position.
        '''
        coords = self._surface_zone.all_coords
        radius = self._surface_zone.distance
        if coords is not None:
            surface_zones(self.all_maps, coords, radius)

    def delete(self):
        self._stop_spotlight_mode()
        for (mgr, h) in self._mgr_handlers:
            try:
                mgr.triggers.remove_handler(h)
            except:
                continue
        super().delete()
