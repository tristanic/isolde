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
    def live_xmapsets(self):
        from .xmapset import XmapSet_Live
        return (m for m in self.child_models() if isinstance(m, XmapSet_Live))

    @property
    def static_xmapsets(self):
        from .xmapset import XmapSet_Static
        return (m for m in self.child_models() if isinstance(m, XmapSet_Static))

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
        from .clipper_python import Coord_orth
        self._box_center_grid = Coord_orth(cofr).coord_frac(self.cell).coord_grid(self.grid)
        dim = self._box_dimensions = \
            2 * calculate_grid_padding(radius, self.grid, self.cell)
        self._box_corner_grid, self._box_corner_xyz = _find_box_corner(
            cofr, self.cell, self.grid, radius)
        self.triggers.activate_trigger('map box changed',
            ((self._box_corner_xyz, self._box_corner_grid, dim), False))
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
        and maximum grid coordinates. Automatically turns off live scrolling.
        '''
        self.spotlight_mode = False
        from .clipper_python import Coord_grid
        cmin = Coord_grid(minmax[0])
        cmin_xyz = cmin.coord_frac(self.grid).coord_orth(self.cell).xyz
        dim = (minmax[1]-minmax[0])
        self.triggers.activate_trigger('map box changed',
            ((cmin_xyz, cmin, dim), force_fill))

    def update_spotlight(self, trigger_name, new_center, force=True):
        '''
        Update the position of the "spotlight" to surround the current centre of
        rotation. If this manager is not currently displayed, then filling the
        volumes with data around the new position will be deferred unless
        force is set to True.
        '''
        center, center_grid = new_center
        self._box_center = center
        self._box_center_grid = center_grid
        if not self.visible and not force:
            # Just store the box parameters for when we're re-displayed
            return
        if self.spotlight_mode:
            box_corner_grid, box_corner_xyz = _find_box_corner(
                center, self.cell, self.grid, self.display_radius)
            self.triggers.activate_trigger('map box moved',
                (box_corner_xyz, box_corner_grid, self._box_dimensions))
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

    def update_box(self, trigger_name, new_center, force=True):
        '''Update the map box to surround the current centre of rotation.'''
        center, center_grid = new_center
        self._box_center = center
        self._box_center_grid = center_grid
        if not self.visible and not force:
            # Just store the box parameters for when we're re-displayed
            return
        if self.spotlight_mode:
            box_corner_grid, box_corner_xyz = _find_box_corner(
                center, self.cell, self.grid, self.display_radius)
            self.triggers.activate_trigger(
                'map box moved',
                (box_corner_xyz, box_corner_grid, self._box_dimensions)
            )
            self._surface_zone.update(
                self.display_radius, coords = numpy.array([center]))
            self._reapply_zone()


def calculate_grid_padding(radius, grid, cell):
    '''
    Calculate the number of grid steps needed on each crystallographic axis
    in order to capture at least radius angstroms in x, y and z.
    '''
    corner_mask = numpy.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,0,0],[1,1,1]])
    corners = (corner_mask * radius).astype(float)
    grid_upper = numpy.zeros([8,3], numpy.int)
    grid_lower = numpy.zeros([8,3], numpy.int)
    from .clipper_python import Coord_orth
    for i, c in enumerate(corners):
        co = Coord_orth(c)
        cm = co.coord_frac(cell).coord_map(grid).uvw
        grid_upper[i,:] = numpy.ceil(cm).astype(int)
        grid_lower[i,:] = numpy.floor(cm).astype(int)
    return grid_upper.max(axis=0) - grid_lower.min(axis=0)


def _find_box_corner(center, cell, grid, radius = 20):
    '''
    Find the bottom corner (i.e. the origin) of a rhombohedral box
    big enough to hold a sphere of the desired radius.
    '''
    from .clipper_python import Coord_frac, Coord_orth, Coord_grid
    radii_frac = Coord_frac(radius/cell.dim)
    center_frac = Coord_orth(center).coord_frac(cell)
    bottom_corner_grid = center_frac.coord_grid(grid) \
                - Coord_grid(calculate_grid_padding(radius, grid, cell))
    bottom_corner_orth = bottom_corner_grid.coord_frac(grid).coord_orth(cell)
    return bottom_corner_grid, bottom_corner_orth.xyz

def _get_bounding_box(coords, padding, grid, cell):
    '''
    Find the minimum and maximum grid coordinates of a box which will
    encompass the given (x,y,z) coordinates plus padding (in Angstroms).
    '''
    from .clipper_python import Util, Coord_grid
    grid_pad = calculate_grid_padding(padding, grid, cell)
    box_bounds_grid = Util.get_minmax_grid(coords, cell, grid)\
                        + numpy.array((-grid_pad, grid_pad))
    box_origin_grid = box_bounds_grid[0]
    box_origin_xyz = Coord_grid(box_origin_grid).coord_frac(grid).coord_orth(cell)
    dim = box_bounds_grid[1] - box_bounds_grid[0]
    return [box_origin_grid, box_origin_xyz, dim]
