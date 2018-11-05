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
    def __init__(self, crystal_manager):
        self._mgr = crystal_manager
        super().__init__('Map Manager', manager.session)
        self._live_xmapsets = []
        self._static_xmapsets = []
        self._nxmapsets = []

        from chimerax.core.triggerset import TriggerSet
        trig = self._triggers = TriggerSet()

        trigger_names = (
            # Change shape/size of box for map viewing.
            # Data is a tuple of (xyz centre, xyz min, xyz max) denoting the
            # rectangular prism that needs to be covered. It is up to
            # each XmapSet or NXmapHandler to work out the necessary box size
            # on its own grid.
            'map box changed',  # Changed shape of box for map viewing
            'map box moved',    # Just changed the centre of the box
        )
        for t in trigger_names:
            trig.add_trigger(t)

        # Handler for live box update
        self._box_update_handler = None
        # Is the box already initialised?
        self._box_initialized = False
        # Object storing the parameters required for masking (used after
        # adjusting contours)
        self._surface_zone = Surface_Zone(display_radius, None, None)
        # Is the map box moving with the centre of rotation?
        self._spotlight_mode = False

        # Radius of the sphere in which the map will be displayed when
        # in live-scrolling mode
        self.display_radius = display_radius
        # Actual box dimensions in (u,v,w) grid coordinates
        self._box_dimensions = None
        # Centre of the box (used when tracking the centre of rotation)
        self._box_center = None
        # Last grid coordinate of the box centre. We only need to update
        # the map if the new centre maps to a different integer grid point
        self._box_center_grid = None
        # Minimum corner of the box in (x,y,z) coordinates. This must
        # correspond to the grid coordinate in self._box_corner_grid
        self._box_corner_xyz = None
        # Minimum corner of the box in grid coordinates
        self._box_corner_grid = None


        self.spotlight_mode = True

        self._box_initialized = True

        self.display=False
        self._rezone_pending = False
        # Apply the surface mask
        self.session.triggers.add_handler('frame drawn', self._first_init_cb)

    @property
    def xmapsets(self):
        from .xmapset import XmapSet
        return (m for m in self.child_models() if isinstance(m, XmapSet))

    @property
    def nxmapsets(self):
        from .nxmapset import NXmapSet
        return (m for m in self.child_models() if isinstance(m, NXmapSet))

    @property
    def all_maps(self):
        from chimerax.map import Volume
        return (m for m in self.all_models() if isinstance(m, Volume))

    @property
    def crystal_mgr(self):
        return self._mgr

    @property
    def structure(self):
        return self.crystal_mgr.structure

    @property
    def triggers(self):
        return self._triggers

    @property
    def hklinfo(self):
        return self.crystal_mgr.hklinfo

    @property
    def spacegroup(self):
        return self.crystal_mgr.spacegroup

    @property
    def cell(self):
        return self.crystal_mgr.cell

    @property
    def grid(self):
        return self.crystal_mgr.grid

    @property
    def display_radius(self):
        '''Get/set the radius (in Angstroms) of the live map display sphere.'''
        return self._display_radius

    @property
    def box_center(self):
        return self._box_center

    @property
    def box_params(self):
        return (self._box_corner_xyz, self._box_corner_grid, self._box_dimensions)

    @display_radius.setter
    def display_radius(self, radius):
        self._display_radius = radius
        v = self.session.view
        cofr = self._box_center = v.center_of_rotation
        xyz_min = cofr-radius
        xyz_max = cofr+radius

        self.triggers.activate_trigger('map box changed',
            (cofr, xyz_min, xyz_max)
        )
        self._surface_zone.update(radius, coords = numpy.array([cofr]))
        self._reapply_zone()

    @property
    def spotlight_mode(self):
        '''Turn live map scrolling on and off.'''
        return self._spotlight_mode

    @spotlight_mode.setter
    def spotlight_mode(self, switch):
        if switch:
            self.position = Place()
            if not self._spotlight_mode:
                '''
                Set the box dimensions to match the stored radius.
                '''
                self.display_radius = self._display_radius
            self._start_spotlight_mode()
        else:
            self._stop_spotlight_mode()

    @property
    def display(self):
        return super().display

    @display.setter
    def display(self, switch):
        if switch:
            if self.spotlight_mode:
                self._start_spotlight_mode()
        Model.display.fset(self, switch)

    def _start_spotlight_mode(self):
        if self._box_update_handler is None:
            self._box_update_handler = self.crystal_mgr.triggers.add_handler(
                'box center moved', self.update_box)
        if self._box_center is not None:
            self.update_box(None, (self._box_center, self._box_center_grid), force=True)
        self.positions = Places()
        self._spotlight_mode = True

    def _stop_spotlight_mode(self):
        if self._box_update_handler is not None:
            self.crystal.triggers.remove_handler(self._box_update_handler)
            self._box_update_handler = None
        self._spotlight_mode = False

    def all_maps(self):
        '''
        Return all maps (of all types) managed by this object.
        '''
        pass

    def cover_box(self, minmax, force_fill=False):
        '''
        Set the map box to fill a volume encompassed by the provided minimum
        and maximum xyz coordinates. Automatically turns off live scrolling.
        '''
        self.spotlight_mode = False
        xyz_min, xyz_max = minmax
        center = (xyz_min + xyz_max)/2
        self.triggers.activate_trigger('map box changed',
            (center, xyz_min, xyz_max))

    def update_spotlight(self, trigger_name, new_center, force=True):
        '''
        Update the position of the "spotlight" to surround the current centre of
        rotation. If this manager is not currently displayed, then filling the
        volumes with data around the new position will be deferred unless
        force is set to True.
        '''
        self._box_center = new_center
        if not self.visible and not force:
            # Just store the box parameters for when we're re-displayed
            return
        if self.spotlight_mode:
            self.triggers.activate_trigger('map box moved',
                new_center)
            self._surface_zone.update(self.display_radius, coords = numpy.array([center]))
            self._reapply_zone()

    def add_xmapset_live(self):
        '''
        Add a set of live-updating maps calculated from atomic coordinates and
        experimental amplitudes.
        '''
        pass

    def add_xmapset_static(self):
        '''
        Add a set of crystallographic maps based on pre-calculated structure
        factors.
        '''
        pass

    def add_nxmapset(self):
        '''
        Add a handler for real-space (non-crystallographic maps)
        '''
        pass

    def _first_init_cb(self, *_):
        self.display = True
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def rezone_needed(self):
        if not self._rezone_pending:
            self.session.triggers.add_handler('new frame', self._rezone_once_cb)

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
        self.spotlight_mode = False
        self.live_update = False
        super().delete()
