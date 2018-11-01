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

    def add_map_handler(self):
        raise NotImplementedError('Function not valid for base class!')
