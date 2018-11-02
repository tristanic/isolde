import numpy

class MapSet_Base(Model):
    '''
    Base class for XmapSet_Live, XmapSet_Static and NXmapSet. Provides basic
    methods for visualisation, masking etc.
    '''

    # Default contour levels and colours for generic maps. Override in the
    # derived class if you wish

    STANDARD_LOW_CONTOUR = numpy.array([1.5])
    STANDARD_HIGH_CONTOUR = numpy.array([2.0])
    STANDARD_DIFFERENCE_MAP_CONTOURS = numpy.array([-3.0, 3.0])

    DEFAULT_MESH_MAP_COLOR = [0,1.0,1.0,1.0] # Solid cyan
    DEFAULT_SOLID_MAP_COLOR = [0,1.0,1.0,0.4] # Transparent cyan
    DEFAULT_DIFF_MAP_COLORS = [[1.0,0,0,1.0],[0,1.0,0,1.0]] #Solid red and green

    def __init__(self, manager, name):
        super().__init__(name, manager.session)
        self._mgr = manager
        manager.add([self])

    @property
    def master_map_mgr(self):
        return self._mgr

    @property
    def crystal_mgr(self):
        return self.master_map_mgr.crystal_mgr

    @property
    def structure(self):
        return self.crystal_mgr.structure

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
        return self.master_map_mgr.display_radius

    def __getitem__(self, name_or_index):
        '''Get one of the child maps by name or index.'''
        if type(name_or_index) == str:
            for m in self.child_models():
                if m.name == name_or_index:
                    return m
            raise KeyError('No map with that name!')
        else:
            return self.child_models()[name_or_index]

    @property
    def spotlight_mode(self):
        return self.master_map_mgr.spotlight_mode

    @spotlight_mode.setter
    def spotlight_mode(self, switch):
        raise NotImplementedError('Spotlight mode can only be enabled/disabled '
            'via the master map manager!')

_pad_base = numpy.array([-1,1], numpy.int)
class XmapSet_Base(MapSet_Base):
    def expand_to_cover_coords(self, coords, padding):
        from .map_mgr import calculate_grid_padding
        cell = self.cell
        grid = self.grid
        pad = calculate_grid_padding(padding, grid, cell)
        from ..clipper_util import get_minmax_grid
        box_bounds_grid = \
            get_minmax_grid(coords, cell, grid) + _pad_base*pad
        self.set_box_limits(box_bounds_grid, force_fill=True)

    def set_box_limits(self, minmax, force_fill = False):
        '''
        Set the map box to fill a volume encompassed by the provided minimum
        and maximum grid coordinates. Automatically turns off live scrolling.
        '''
        self.live_scrolling = False
        from .clipper_python import Coord_grid
        cmin = Coord_grid(minmax[0])
        cmin_xyz = cmin.coord_frac(self.grid).coord_orth(self.cell).xyz
        dim = (minmax[1]-minmax[0])
        for v in self:
            v._box_changed_cb('map box changed',
            ((cmin_xyz, cmin, dim), force_fill))
