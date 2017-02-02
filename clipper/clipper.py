from .lib import clipper_python_core as clipper_core
import numpy

#### Message logging
clipper_messages = clipper_core.ClipperMessageStream()

def log_clipper(func):
    def func_wrapper(*args, **kwargs):
        clipper_messages.clear()
        func(*args, **kwargs)
        message_string = clipper_messages.read_and_clear()
        if message_string:
            session.logger.info("CLIPPER WARNING:")
            session.logger.info(message_string)
    return func_wrapper
    
#################################################################
# To make things a little more user-friendly than the raw SWIG wrapper,
# all the major Clipper classes are sub-classed here with some extra
# documentation, useful Python methods, @property decorators, etc.
# This requires overcoming one small problem: when Clipper itself creates 
# and returns objects, they're returned as the base class rather than
# the sub-class. So, we have to over-ride the __new__ method for each
# base class. 
#####################################

def __newAtom__(cls, *args, **kwargs):
    if cls == clipper_core.clipper.Atom:
        return object.__new__(Atom)
    return object.__new__(cls)        
clipper_core.Atom.__new__ = staticmethod(__newAtom__)

class Atom(clipper_core.Atom):
    '''
    A minimalist atom object containing only the properties required for
    electron density calculations
    '''
    ATOM_NAMES = ['H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',
                  'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar',
                  'K',  'Ca', 'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co',
                  'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
                  'Rb', 'Sr', 'Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
                  'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',  'Xe',
                  'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',
                  'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
                  'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl',
                  'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
                  'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 
                  'H1-',  'Li1+', 'Be2+', 'Cval', 'O1-',  'O2-',  'F1-',
                  'Na1+', 'Mg2+', 'Al3+', 'Siva', 'Si4+', 'Cl1-', 'K1+',
                  'Ca2+', 'Sc3+', 'Ti2+', 'Ti3+', 'Ti4+', 'V2+',  'V3+',
                  'V5+',  'Cr2+', 'Cr3+', 'Mn2+', 'Mn3+', 'Mn4+', 'Fe2+',
                  'Fe3+', 'Co2+', 'Co3+', 'Ni2+', 'Ni3+', 'Cu1+', 'Cu2+',
                  'Zn2+', 'Ga3+', 'Ge4+', 'Br1-', 'Rb1+', 'Sr2+', 'Y3+',
                  'Zr4+', 'Nb3+', 'Nb5+', 'Mo3+', 'Mo5+', 'Mo6+', 'Ru3+',
                  'Ru4+', 'Rh3+', 'Rh4+', 'Pd2+', 'Pd4+', 'Ag1+', 'Ag2+',
                  'Cd2+', 'In3+', 'Sn2+', 'Sn4+', 'Sb3+', 'Sb5+', 'I1-',
                  'Cs1+', 'Ba2+', 'La3+', 'Ce3+', 'Ce4+', 'Pr3+', 'Pr4+',
                  'Nd3+', 'Pm3+', 'Sm3+', 'Eu2+', 'Eu3+', 'Gd3+', 'Tb3+',
                  'Dy3+', 'Ho3+', 'Er3+', 'Tm3+', 'Yb2+', 'Yb3+', 'Lu3+',
                  'Hf4+', 'Ta5+', 'W6+',  'Os4+', 'Ir3+', 'Ir4+', 'Pt2+',
                  'Pt4+', 'Au1+', 'Au3+', 'Hg1+', 'Hg2+', 'Tl1+', 'Tl3+',
                  'Pb2+', 'Pb4+', 'Bi3+', 'Bi5+', 'Ra2+', 'Ac3+', 'Th4+',
                  'U3+',  'U4+',  'U6+',  'Np3+', 'Np4+', 'Np6+', 'Pu3+',
                  'Pu4+', 'Pu6+']  
                        
    def __init__(self):
        clipper_core.Atom.__init__(self)
    
    
    @property
    def element(self):
        '''
        The standard abbreviated element (or elemental ion) name. All valid
        names are listed in Atom.ATOM_NAMES.
        '''
        return super(Atom, self).element()
    
    @element.setter
    def element(self, element_name):
        # Check common atom names first to improve performance
        if element_name not in ('H', 'C', 'N', 'O', 'S'):
            if element_name not in self.ATOM_NAMES:
                raise TypeError('Unrecognised element!')
        super(Atom, self).set_element(element_name)
    
    @property
    def coord(self):
        '''
        Get the coordinates of the atom as a numpy array
        '''
        return self.coord_orth().xyz()
    
    @coord.setter
    def coord(self, coord):
        '''
        Set the coordinates of the atom using a list or numpy array
        '''
        self.set_coord_orth(clipper_core.Coord_orth(*coord))

    @property
    def coord_orth(self):
        '''
        Get the Clipper Coord_orth object associated with this atom
        '''
        return super(Atom, self).coord_orth()
    
    @coord_orth.setter
    def coord_orth(self, coord):
        self.coord = coord
        
    @property
    def occupancy(self):
        return super(Atom, self).occupancy()
    
    @occupancy.setter
    def occupancy(self, occ):
        self.set_occupancy(occ)
    
    @property
    def u_iso(self):
        '''
        Get the isotropic b-factor.
        '''
        return super(Atom, self).u_iso()
    
    @u_iso.setter
    def u_iso(self, u_iso):
        self.set_u_iso(u_iso)
    
    @property
    def b_factor(self):
        '''
        Synonym for u_iso.
        '''
        return self.u_iso
    
    @b_factor.setter
    def b_factor(self, b_factor):
        self.u_iso = b_factor
    
    @property
    def u_aniso_orth(self):
        '''
        Get the anisotropic B-factor matrix as a 6-member array:
        [u11, u22, u33, u12, u13, u23].
        '''
        return super(Atom, self).u_aniso_orth().get_vals()
    
    @property
    def _u_aniso_orth(self):
        '''
        Get the Clipper::U_aniso_orth object
        '''
        return super(Atom, self).u_aniso_orth()
    
    @u_aniso_orth.setter
    def u_aniso_orth(self, u_aniso):
        '''
        Set the anisotropic B-factor matrix using a 6-member array:
        [u11, u22, u33, u12, u13, u23].
        '''
        self.set_u_aniso_orth(*u_aniso)
    
    @property
    def is_null(self):
        return super(Atom, self).is_null()

def __newCoord_orth__(cls, *args, **kwargs):
    if cls == clipper_core.Coord_orth:
        return object.__new__(Coord_orth)
    return object.__new__(cls)        
clipper_core.Coord_orth.__new__ = staticmethod(__newCoord_orth__)

class Coord_orth(clipper_core.Coord_orth):
    '''
    Coordinates in orthographic (x,y,z) space.
    '''
    def __init__(self, xyz):
        if isinstance(xyz, numpy.ndarray):
            # Because SWIG does not correctly typemap numpy.float32
            xyz = xyz.astype(float)
        clipper_core.Coord_orth.__init__(self, *xyz)
    
    @property
    def x(self):
        return super(Coord_orth, self).x()
    
    @property
    def y(self):
        return super(Coord_orth, self).y()
        
    @property
    def z(self):
        return super(Coord_orth, self).z()

    @property
    def xyz(self):
        return super(Coord_orth, self).xyz()

def __newCoord_grid__(cls, *args, **kwargs):
    if cls == clipper_core.Coord_grid:
        return object.__new__(Coord_grid)
    return object.__new__(cls)        
clipper_core.Coord_grid.__new__ = staticmethod(__newCoord_grid__)
        
class Coord_grid(clipper_core.Coord_grid):
    '''
    Integer grid coordinates in crystal space.
    '''
    def __init__(self, uvw):
        clipper_core.Coord_grid.__init__(self, *uvw)
    
    @property
    def u(self):
        return super(Coord_grid, self).u()
    
    @property
    def v(self):
        return super(Coord_grid, self).v()
    
    @property
    def w(self):
        return super(Coord_grid, self).w()
    
    @property
    def uvw(self):
        return super(Coord_grid, self).uvw()

def __newCoord_map__(cls, *args, **kwargs):
    if cls == clipper_core.Coord_map:
        return object.__new__(Coord_map)
    return object.__new__(cls)        
clipper_core.Coord_map.__new__ = staticmethod(__newCoord_map__)
    
class Coord_map(clipper_core.Coord_map):
    '''
    Like Coord_grid, but allowing non-integer values.
    '''
    def __init__(self, uvw):
        if isinstance(uvw, numpy.ndarray):
            # Because SWIG does not correctly typemap numpy.float32
            uvw = uvw.astype(float)
        clipper_core.Coord_map.__init__(self, *uvw)
    
    @property
    def u(self):
        return super(Coord_map, self).u()
    
    @property
    def v(self):
        return super(Coord_map, self).v()
    
    @property
    def w(self):
        return super(Coord_map, self).w()
    
    @property
    def uvw(self):
        return super(Coord_map, self).uvw()

def __newCoord_frac__(cls, *args, **kwargs):
    if cls == clipper_core.Coord_frac:
        return object.__new__(Coord_frac)
    return object.__new__(cls)        
clipper_core.Coord_frac.__new__ = staticmethod(__newCoord_frac__)

class Coord_frac(clipper_core.Coord_frac):
    '''
    Fractional coordinates along unit cell axes (a,b,c), scaled to the 
    range 0..1 on each axis.
    '''
    def __init__(self, uvw):
        if isinstance(uvw, numpy.ndarray):
            # Because SWIG does not correctly typemap numpy.float32
            uvw = uvw.astype(float)
        clipper_core.Coord_frac.__init__(self, *uvw)

    @property
    def u(self):
        return super(Coord_frac, self).u()
    
    @property
    def v(self):
        return super(Coord_frac, self).v()
    
    @property
    def w(self):
        return super(Coord_frac, self).w()
    
    @property
    def uvw(self):
        return super(Coord_frac, self).uvw()

def __newCCP4MTZfile__(cls, *args, **kwargs):
    if cls == clipper_core.CCP4MTZfile:
        return object.__new__(CCP4MTZfile)
    return object.__new__(cls)        
clipper_core.CCP4MTZfile.__new__ = staticmethod(__newCCP4MTZfile__)

class CCP4MTZfile(clipper_core.CCP4MTZfile):
    '''
    MTZ import/export parent class for clipper objects.
    
    This is the import/export class which can be linked to an mtz
    file and used to transfer data into or out of a Clipper data structure.
    
    Note that to access the MTZ file efficiently, data reads and writes
    are deferred until the file is closed.

    MTZ column specification:

    Note that the specification of the MTZ column names is quite
    versatile. The MTZ crystal and dataset must be specified, although
    the wildcard '*' may replace a complete name. Several MTZ columns
    will correspond to a single datalist. This may be handled in two
    ways:
    
    - A simple name. The corresponding MTZ columns will be named
    after the datalist name, a dot, the datalist type, a dot, and a
    type name for the indivudal column,
    i.e. /crystal/dataset/datalist.listtype.coltype This is the
    default Clipper naming convention for MTZ data.
    
    - A comma separated list of MTZ column names enclosed in square
    brackets.  This allows MTZ data from legacy applications to be
    accessed.
    
    Thus, for example, an MTZPATH of
    
        native/CuKa/fsigfdata
    
    expands to MTZ column names of
    
        fsigfdata.F_sigF.F
        fsigfdata.F_sigF.sigF
    
    with a crystal called "native" and a dataset called "CuKa". An MTZPATH of
    
        native/CuKa/[FP,SIGFP]
    
    expands to MTZ column names of
    
        FP
        SIGFP
    
    with a crystal called "native" and a dataset called "CuKa".

   Import/export types:

    For an HKL_data object to be imported or exported, an MTZ_iotype
    for that datatype must exist in the MTZ_iotypes_registry. MTZ_iotypes 
    are defined for all the built-in datatypes. If you need to store a 
    user defined type in an MTZ file, then register that type with the 
    MTZ_iotypes_registry. 

    EXAMPLE: Loading essential crystal information and 
    2Fo-Fc amplitudes/phases from an mtz file
    
    fphidata =  HKL_data_F_phi_double()    
    myhkl = hklinfo()
    mtzin = CCP4MTZfile()
    mtzin.open_read(filename)
    mtzin.import_hkl_info(myhkl)
    mtzin.import_hkl_data(fphidata, '/crystal/dataset/[2FOFCWT, PH2FOFCWT]')
    mtzin.close_read()
    '''
    
    def __init__(self):
        clipper_core.CCP4MTZfile.__init__(self)
        
    @log_clipper
    def import_hkl_data(self, cdata, mtzpath):
        return super(CCP4MTZfile, self).import_hkl_data(cdata, mtzpath)

def __newCIFfile__(cls, *args, **kwargs):
    if cls == clipper_core.CIFfile:
        return object.__new__(CIFfile)
    return object.__new__(cls)        
clipper_core.CIFfile.__new__ = staticmethod(__newCIFfile__)
        
class CIFfile(clipper_core.CIFfile):
    '''
    CIF import/export parent class for clipper objects.
      
    This is the import class which can be linked to an cif data
    file and be used to transfer data into a Clipper
    data structure. 
    It is currently a read-only class.
    '''
    
    def __init__(self):
        clipper_core.CIFfile.__init__(self)
    
    @log_clipper
    def resolution(cell):
        return super(CIFfile, self).resolution(cell)



def __newHKL_info__(cls, *args, **kwargs):
    if cls == clipper_core.HKL_info:
        return object.__new__(HKL_info)
    return object.__new__(cls)        
clipper_core.HKL_info.__new__ = staticmethod(__newHKL_info__)

class HKL_info(clipper_core.HKL_info):
    def __init__(self, session):
        clipper_core.HKL_info.__init__(self)
        self.session = session

def __newHKL_data_F_phi__(cls, *args, **kwargs):
    if cls == clipper_core.HKL_data_F_phi_double:
        return object.__new__(HKL_data_F_phi)
    return object.__new__(cls)        
clipper_core.HKL_data_F_phi_double.__new__ = staticmethod(__newHKL_data_F_phi__)

class HKL_data_F_phi (clipper_core.HKL_data_F_phi_double):
    def __init__(self):
        clipper_core.HKL_data_F_phi_double.__init__(self)



def __newGrid_sampling__(cls, *args, **kwargs):
    if cls == clipper_core.Grid_sampling:
        return object.__new__(Grid_sampling)
    return object.__new__(cls)        
clipper_core.Grid_sampling.__new__ = staticmethod(__newGrid_sampling__)

class Grid_sampling(clipper_core.Grid_sampling):
    '''
    Object defining the grid used to sample points in a 3D map.
    '''
    def __init__(self, spacegroup, cell, resolution, rate = 1.5):
        clipper_core.Grid_sampling.__init__(self, spacegroup, cell, resolution, rate)




def __newXmap_double__(cls, *args, **kwargs):
    if cls == clipper_core.Xmap_double:
        return object.__new__(Xmap)
    return object.__new__(cls)        
clipper_core.Xmap_double.__new__ = staticmethod(__newXmap_double__)
    
class Xmap(clipper_core.Xmap_double):
    '''
    A Clipper crystallographic map generated from reciprocal space data.
    '''
    def __init__(self, session, name, spacegroup, cell, grid_sam):
        clipper_core.Xmap_double.__init__(self, spacegroup, cell, grid_sam)
        self.session = session
        self.name = name
        # Some extra useful variables that aren't directly available from
        # the Clipper API
        c = self.cell()
        g = self.grid_sampling()
        self.grid_samples = g.dim()
        
        ###
        # Basic stats. Can only be set after the map has been filled with an
        # FFT. Returned as a tuple in the below order by self.stats()
        ###
        self._max = None
        self._min = None
        self._mean = None
        self._sigma = None
        self._skewness = None
        self._kurtosis = None
        
        # Default "radius" (actually half-width of the rhombohedron) of
        # the box in which the map will be drawn.
        self.box_radius = None
        # Box dimensions in grid points (depends on resolution)
        self._box_dimensions = None
        # Centre point of the box (will typically be set to the centre of
        # rotation).
        self.box_center = None
        # Vector mapping the box centre to its origin
        self.box_origin_offset = None
        # ChimeraX Volume object to draw the map into
        self.volume = None
        # Array_Grid_Data object held by the Volume object
        self._array_grid_data = None
        # Numpy array to send the map data to
        self._box_data = None
        # session.triggers handler for live box update
        self._box_handler = None
        # Distance any axis of the cofr has to move before triggering an update.
        self._cofr_eps = numpy.min(self.voxel_size()) / 2
        
        # Model object to hold Drawings defining the special positions
        self.special_positions_model = None
       
    def recalculate_stats(self):
        '''
        Force recalculation of map statistics (max, min, mean, sigma, 
        skewness and kurtosis).
        '''
        self._min, self._max, self._mean, \
            self._sigma, self._skewness, self._kurtosis = self.stats()         
    
    @property
    def max(self):
        if self._max is None:
            self.recalculate_stats()
        return self._max

    @property
    def min(self):
        if self._min is None:
            self.recalculate_stats()
        return self._min

    @property
    def mean(self):
        if self._mean is None:
            self.recalculate_stats()
        return self._mean
    
    @property
    def sigma(self):
        if self._sigma is None:
            self.recalculate_stats()
        return self._sigma
    
    @property
    def skewness(self):
        if self._skewness is None:
            self.recalculate_stats()
        return self._skewness

    @property
    def kurtosis(self):
        if self._max is None:
            self.recalculate_stats()
        return self._kurtosis
    
    
    def _box_go_live(self):
        if self._box_handler is None:
            self._box_handler = self.session.triggers.add_handler('new frame', self.update_box)
    
    def _box_go_static(self):
        if self._box_handler is not None:
            self.session.triggers.remove_handler(self._box_handler)
            self._box_handler = None
    
    def initialize_box_display(self, radius = 15, pad = 0):
        '''
        Generate a Volume big enough to hold a sphere of the given radius,
        plus an optional padding of pad voxels on each side. The session
        view will be changed to use an orthographic camera, with the centre
        of rotation updating to always remain at the centre of the screen.
        The volume will automatically track the centre of rotation and 
        update its position and contents to reflect the local density.
        '''
        from chimerax.core.commands import camera, cofr
        camera.camera(self.session, 'ortho')
        cofr.cofr(self.session, 'centerOfView')
        self._box_pad = pad
        self.box_radius = radius
        v = self.session.view
        c = self.cell()
        g = self.grid_sampling()
        self.box_center = v.center_of_rotation
        box_corner_grid, box_corner_xyz = self._find_box_corner(self.box_center, radius, pad)
        self.box_origin_offset = box_corner_xyz - self.box_center
        self._box_dimensions = (numpy.ceil(radius / self.voxel_size() * 2)+pad*2).astype(int)[::-1]
        data = self._box_data = numpy.empty(self._box_dimensions, numpy.double)
                
        from chimerax.core.map.data import Array_Grid_Data
        from chimerax.core.map import Volume
        self._array_grid_data = Array_Grid_Data(data, origin = box_corner_xyz,
                step = self.voxel_size(), cell_angles = c.angles_deg())
        self.volume = Volume(self._array_grid_data, self.session)
        self._fill_volume_data(data, box_corner_grid)
        self.volume.initialize_thresholds()
        self.session.models.add([self.volume])
        contour_val = self.contour = 1.5 * self.sigma
        self.change_contour(contour_val)
        self.volume.set_color([0.5,1.0,0.5,0.6])
        self.volume.show()
        self._box_go_live()

    def change_box_radius(self, radius, pad=0):
        self._box_go_static()
        v = self.session.view
        cofr = v.center_of_rotation
        self.box_radius = radius
        dim = (numpy.ceil(radius / self.voxel_size() * 2)+pad*2).astype(int)[::-1]
        data = numpy.empty(dim, numpy.double)
        box_corner_grid, box_corner_xyz = self._find_box_corner(cofr, radius, pad)
        from chimerax.core.map.data import Array_Grid_Data
        darray = Array_Grid_Data(data, origin = box_corner_xyz,
            step = self.voxel_size(), cell_angles = self.cell().angles_deg())
        self.volume.replace_data(darray)
        self._box_pad = pad
        self._box_dimensions = dim
        self._box_data = data
        self._array_grid_data = darray
        self.update_box(force_update=True)
        self._box_go_live()
        
        
        
    def update_box(self, *_, force_update = False):
        v = self.session.view
        cofr = v.center_of_rotation
        if not force_update:
            if numpy.all(abs(self.box_center - cofr) < self._cofr_eps):
                return
        
        self.box_center = cofr
        box_corner_grid, box_corner_xyz = self._find_box_corner(cofr, self.box_radius, self._box_pad)
        self._array_grid_data.set_origin(box_corner_xyz)
        self._fill_volume_data(self._box_data, box_corner_grid)
        self.volume.update_surface()
        
    def change_contour(self, contour_val):
        from chimerax.core.map import volumecommand
        volumecommand.volume(self.session,[self.volume], level=[[contour_val]], cap_faces = False)
        

    def _fill_volume_data(self, target, start_grid_coor):
        shape = (numpy.array(target.shape)[::-1] - 1).tolist()
        end_grid_coor = start_grid_coor + Coord_grid(shape)
        count = self.export_section_numpy(target, start_grid_coor, end_grid_coor, 'C')
                
    
    def _find_box_corner(self, center, radius = 20, pad = 1):
        '''
        Find the bottom corner (i.e. the origin) of a rhombohedral box
        big enough to hold a sphere of the desired radius, and padded by
        pad voxels in each dimension.
        '''
        cell = self.cell()
        grid = self.grid_sampling()
        voxel_size = self.voxel_size()
        radii_in_voxels = radius / voxel_size
        radii_frac = Coord_frac(radii_in_voxels * self.voxel_size_frac())
        center_ortho = Coord_orth(center)
        center_frac = center_ortho.coord_frac(cell)
        bottom_corner_frac = center_frac - radii_frac
        bottom_corner_grid = bottom_corner_frac.coord_grid(grid) -\
                            Coord_grid([pad,pad,pad])
        bottom_corner_orth = bottom_corner_grid.coord_frac(grid).coord_orth(cell)
        return bottom_corner_grid, bottom_corner_orth.xyz
        
                                    
    def draw_unit_cell_and_special_positions(self, model, offset = None):
        from chimerax.core.models import Model, Drawing
        from chimerax.core.geometry import Place, Places
        from chimerax.core.surface.shapes import sphere_geometry
        import copy
        
        ref = model.bounds().center().astype(float)
        frac_coords = Coord_orth(ref).coord_frac(self.cell()).uvw
        if offset is None:
            offset = numpy.array([0,0,0],int)
        corners_frac = numpy.array([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]],numpy.double) + offset
        corners = []
        for c in corners_frac:
            cf = Coord_frac(c).coord_orth(self.cell())
            corners.append(cf.xyz)
        m = self.special_positions_model
        
        if m is None or m.deleted:
            m = self.special_positions_model = Model('Special Positions',self.session)
            spc = numpy.array(self.special_positions_unit_cell_xyz(offset))
            coords = spc[:,0:3]
            multiplicity = spc[:,3].astype(int)
            sphere = numpy.array(sphere_geometry(80))
            sphere[0]*=0.5
            scale_2fold = numpy.identity(3)
            scale_3fold = numpy.identity(3)* 1.5
            scale_4fold = numpy.identity(3)* 2
            scale_6fold = numpy.identity(3)* 3
            rgba_2fold = numpy.array([255,255,255,255],numpy.int32)
            rgba_3fold = numpy.array([0,255,255,255],numpy.int32)
            rgba_4fold = numpy.array([255,255,0,255],numpy.int32)
            rgba_6fold = numpy.array([255,0,0,255],numpy.int32)
            rgba_corner = numpy.array([255,0,255,255],numpy.int32)
            positions = []
            colors = []
            d = Drawing('points')
            d.vertices, d.normals, d.triangles = sphere
            
            for coord, mult in zip(coords, multiplicity):
                if mult == 2:
                    positions.append(Place(axes=scale_2fold, origin=coord))
                    colors.append(rgba_2fold)
                elif mult == 3:
                    positions.append(Place(axes=scale_3fold, origin=coord))
                    colors.append(rgba_3fold)
                elif mult == 4:
                    positions.append(Place(axes=scale_4fold, origin=coord))
                    colors.append(rgba_4fold)
                elif mult == 6:
                    positions.append(Place(axes=scale_6fold, origin=coord))
                    colors.append(rgba_6fold)
            for c in corners:
                positions.append(Place(axes=scale_6fold, origin=c))
                colors.append(rgba_corner)
                
            d.positions = Places(positions)
            d.colors = numpy.array(colors)
            m.add_drawing(d)
            self.session.models.add([m])
            
        
                
                
            
            
        
        
         
                
        
        
    
    
def box_corners(origin_xyz, size_xyz):
    ret = []
    minmax = [origin_xyz, size_xyz]
    for i in range(2):
        for j in range(2):
            for k in range(2):
                ret.append([minmax[i][0],minmax[j][1], minmax[k][2]])
    return ret
                


def apply_b_factors_to_hydrogens(atom_list):
    '''
    When hydrogens are added to a pre-existing structure, their isotropic
    and/or anisotropic B-factors are not necessarily set. This routine simply
    sets them to match their adjacent heavy atoms. The input atom_list can
    be either the hydrogens alone, or the whole structure.
    '''
    hydrogens = atom_list.filter(atom_list.element_names == 'H')
    # Get heavy atoms by taking advantage of the fact that hydrogens are
    # only ever covalently bonded to a single other atom.
    for h in hydrogens:
        b = h.bonds[0]
        for atom in b.atoms:
            if atom.element_name != 'H':
                h.bfactor = atom.bfactor
                h.aniso_u6 = atom.aniso_u6
                break
                
 
def atom_list_from_sel(atom_list):
    '''
    Takes a ChimeraX Atoms object, and creates a Clipper Atoms_list object
    from the relevant atom properties.
    '''
    n = len(atom_list)
    elements = atom_list.element_names.tolist()
    coords = atom_list.coords
    occupancies = atom_list.occupancy
    u_iso = atom_list.bfactors
    # Have to be careful here. Atoms.aniso_u6 returns None if any atom in
    # the array has no aniso_u6 entry. 
    u_aniso = atom_list.aniso_u6
    if u_aniso is None:
        u_aniso = numpy.zeros([n,6],numpy.float32)
        # FIXME Once new ChimeraX build arrives with Atoms.has_aniso_u entry
        for i, a in enumerate(atom_list):
            if a.aniso_u6 is not None:
                u_aniso[i] = a.aniso_u6
    clipper_atom_list = Atom_list()
    clipper_atom_list.extend_by(n)

    clipper_atom_list.set_elements(elements)
    clipper_atom_list.set_coord_orth(coords)
    clipper_atom_list.set_occupancies(occupancies)
    clipper_atom_list.set_u_isos(u_iso)
    clipper_atom_list.set_u_anisos(u_aniso)
    return clipper_atom_list

def import_Xmap_from_mtz_test(session, filename):
    myhkl = HKL_info(session)
    fphidata =  HKL_data_F_phi()  
    mtzin = CCP4MTZfile()
    mtzin.open_read(filename)
    mtzin.import_hkl_info(myhkl)
    mtzin.import_hkl_data(fphidata, '/crystal/dataset/[2FOFCWT, PH2FOFCWT]')
    mtzin.close_read()
    name = '2FOFCWT'
    mygrid = Grid_sampling(myhkl.spacegroup(), myhkl.cell(), myhkl.resolution())
    mymap = Xmap(session, name, myhkl.spacegroup(), myhkl.cell(), mygrid)
    mymap.fft_from(fphidata)
    return (fphidata, myhkl, mymap)
    
def vol_box(hklinfo, xmap, min_coor, max_coor):
    cell = hklifno.cell()
    grid = xmap.grid_sampling()
    min_ortho = Coord_orth(min_coor)
    max_ortho = Coord_orth(max_coor)
    min_grid = min_ortho.coord_frac(cell).coord_grid(grid)
    max_grid = max_ortho.coord_frac(cell).coord_grid(grid)    
    
def make_unit_cell(model, xmap, draw = True):
    from chimerax.core.geometry import Place, Places
    cell = xmap.cell()
    atoms = model.atoms
    clipper_atoms = atom_list_from_sel(atoms)
    coord = model.bounds().center()
    coord_orth = Coord_orth(coord)
    coord_frac = coord_orth.coord_frac(cell)
    sg = xmap.spacegroup()
    unit_cell_frac_symops = xmap.unit_cell_symops(coord_frac, clipper_atoms)
    if draw:
        uc_places = []
        for op in unit_cell_frac_symops.symops():
            op_orth = op.rtop_orth(cell)
            uc_places.append(Place(matrix=op_orth.matrix()[0:3,:]))
        ucp = Places(uc_places)
        model.positions = ucp
    return unit_cell_frac_symops

def test_pack_box(model, xmap, size = 50):
    from chimerax.core.geometry import Place, Places
    uc = make_unit_cell(model, xmap, draw = False)
    coord = model.bounds().center()
    box_size = (numpy.ones(3, numpy.int32)*size)
    bo = xmap.all_symops_in_box(coord-box_size/2, box_size, uc)
    p = []
    for b in bo:
        p.append(Place(matrix=b.rtop_orth(xmap.cell()).matrix()[0:3,:]))
    P = Places(p)
    model.positions = P
    
    from chimerax.core.models import Drawing, Model
    from chimerax.core.surface.shapes import sphere_geometry
    d = Drawing('box corners')
    m = Model('box', session)
    d.vertices, d.normals, d.triangles = sphere_geometry(80)
    box_min = Coord_orth(coord-box_size/2).coord_frac(xmap.cell()).coord_grid(xmap.grid_sampling())
    box_max = box_min + Coord_grid(*box_size.tolist())
    minmax = [box_min, box_max]
    dp = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                thiscg = Coord_grid([minmax[i].u(), minmax[j].v(), minmax[k].w()])
                dp.append(Place(origin=thiscg.coord_frac(xmap.grid_sampling()).coord_orth(xmap.cell()).xyz))
    d.positions = Places(dp)
    m.add_drawing(d)
    session.models.add([m])
    
    return bo

def pack_box(model, xmap, box_origin_xyz, size = 100):
    box_size_xyz = numpy.ones(3)*size
    box_origin_xyz = numpy.array(box_origin_xyz)
    model_bounds = model.bounds().box_corners()
    bo = xmap.pack_xyz_box(model_bounds, box_origin_xyz, box_size_xyz);
    from chimerax.core.geometry import Place, Places

    p = []
    for b in bo:
        p.append(Place(matrix=b.rtop_orth(xmap.cell()).matrix()[0:3,:]))
    P = Places(p)
    model.positions = P
    
    from chimerax.core.models import Drawing, Model
    from chimerax.core.surface.shapes import sphere_geometry
    d = Drawing('box corners')
    m = Model('box', session)
    d.vertices, d.normals, d.triangles = sphere_geometry(80)
    minmax = [box_origin_xyz, box_origin_xyz+box_size_xyz]
    dp = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                dp.append(Place(origin=[minmax[i][0],minmax[j][1],minmax[k][2]]))
    d.positions = Places(dp)
    m.add_drawing(d)
    session.models.add([m])
    return bo

def draw_box(min_corner, max_corner, name='box'):
    from chimerax.core.geometry import Place, Places
    from chimerax.core.models import Drawing, Model
    from chimerax.core.surface.shapes import sphere_geometry
    d = Drawing('corners')
    m = Model(name, session)
    d.vertices, d.normals, d.triangles = sphere_geometry(80)
    minmax = [min_corner, max_corner]
    dp = []
    base_color = numpy.array([255,255,255,255])
    color_increment = numpy.array([0,-32,-32,0])
    colors = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                dp.append(Place(origin=[minmax[i][0],minmax[j][1],minmax[k][2]]))
                colors.append(base_color-color_increment*(i+j+k));
    d.positions = Places(dp)
    d.colors = colors
    m.add_drawing(d)
    session.models.add([m])
    return d
    
    

        
def draw_asu(xmap):    
    from chimerax.core.geometry import Place, Places
    from chimerax.core.models import Drawing, Model
    from chimerax.core.surface.shapes import sphere_geometry
    d = Drawing('asu corners')
    m = Model('asu box', session)
    d.vertices, d.normals, d.triangles = sphere_geometry(80)
    asu = xmap.grid_asu()
    grid = xmap.grid_sampling()
    cell = xmap.cell()
    minmax = [asu.min().coord_frac(grid).coord_orth(cell).xyz, asu.max().coord_frac(grid).coord_orth(cell).xyz]
    dp = []
    for i in range(2):
        for j in range(2):
            for k in range(2):
                dp.append(Place(origin=[minmax[i][0],minmax[j][1],minmax[k][2]]))
    d.positions = Places(dp)
    m.add_drawing(d)
    session.models.add([m])
    return d

@log_clipper
def test_log_warn():
    clipper_core.warn_test()

@log_clipper
def test_log_except():
    clipper_core.except_test()

@log_clipper    
def test_log_no_warn():
    pass
