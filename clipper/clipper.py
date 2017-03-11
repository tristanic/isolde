from .lib import clipper_python_core as clipper_core
import numpy

#### Message logging

from .clipper_decorators import *
_clipper_messages = MessageStreamSingleton().clipper_messages

########################################################################
# ATOMS AND ATOMIC PROPERTIES
########################################################################

@mappedclass(clipper_core.Atom)
class Atom(clipper_core.Atom):
    '''
    A minimalist atom object containing only the properties required for
    electron density calculations
    '''
    # The complete list of scatterers found in clipper/core/atomsf.cpp.
    # Clipper itself doesn't check incoming atom names for legality, but
    # it's easy to do so here in Python. Used by setter methods in Atom and
    # Atom_list objects.
    ATOM_NAMES = set(
         ['H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',
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
          'Pu4+', 'Pu6+'])

    def __init__(self, element, coord, occ, u_iso, u_aniso = None, allow_unknown = False):
        '''
        __init__(self, element, coord, occ, u_iso, u_aniso) -> Atom

        Create a new Clipper Atom with the given properties
        Args:
            element (string):
                The standard abbreviated element (or elemental ion)
                name. All valid names are listed in Atom.ATOM_NAMES.
            coord ([float*3] or Clipper.Coord_orth object):
                (x,y,z) coordinates of the atom in Angstroms.
            occ (float):
                The fractional occupancy of the atom
            u_iso (float):
                Isotropic B-factor
            u_aniso ([float * 6]):
                Anisotropic B-factor matrix as a 6-member array:
                [u00, u11, u22, u01, u02, u12].
        '''
        self.allow_unknown = allow_unknown
        clipper_core.Atom.__init__(self)
        self.element = element
        self.coord = coord
        self.occupancy = occ
        self.u_iso = u_iso
        self.u_aniso_orth = u_aniso


    #@property
    #def element(self):
        #'''
        #The standard abbreviated element (or elemental ion) name. All valid
        #names are listed in Atom.ATOM_NAMES.
        #'''
        #return super(Atom, self).element()

    #@element.setter
    #def element(self, element_name):
        ## Check common atom names first to improve performance
        #if not self.allow_unknown:
            #if element_name not in ('H', 'C', 'N', 'O', 'S'):
                #if element_name not in self.ATOM_NAMES:
                    #raise TypeError('Unrecognised element!')
        #super(Atom, self).set_element(element_name)

    #@property
    #def coord(self):
        #'''
        #(x,y,z) coordinates of the atom in Angstroms. Can be set from
        #a Python list, numpy array, or a Clipper Coord_orth object.
        #'''
        #return self.coord_orth.xyz

    #@coord.setter
    #def coord(self, coord):
        #if isinstance(coord, Coord_orth):
            #self.set_coord_orth(coord)
        #else:
            #self.set_coord_orth(Coord_orth(coord))

    #@property
    #def coord_orth(self):
        #'''
        #Clipper Coord_orth object associated with this atom. Will return
        #a coord_orth object, but can be set with a simple list of 3
        #(x,y,z) coordinates.
        #'''
        #return super(Atom, self).coord_orth()

    #@coord_orth.setter
    #def coord_orth(self, coord):
        #self.coord = coord

    #@property
    #def occupancy(self):
        #'''Fractional occupancy of this atom'''
        #return super(Atom, self).occupancy()

    #@occupancy.setter
    #def occupancy(self, occ):
        #self.set_occupancy(occ)

    #@property
    #def u_iso(self):
        #'''Isotropic b-factor in square Angstroms.'''
        #return super(Atom, self).u_iso()

    #@u_iso.setter
    #def u_iso(self, u_iso):
        #self.set_u_iso(u_iso)

    #@property
    #def b_factor(self):
        #'''Isotropic b-factor in square Angstroms.'''
        #return self.u_iso

    #@b_factor.setter
    #def b_factor(self, b_factor):
        #self.u_iso = b_factor

    #@property
    #def u_aniso_orth(self):
        #'''
        #Anisotropic B-factor matrix as a 6-member array:
        #[u00, u11, u22, u01, u02, u12].
        #For purely isotropic values, set this to None
        #'''
        #return super(Atom, self).u_aniso_orth()._get_vals()

    #@property
    #def _u_aniso_orth(self):
        #'''Get the Clipper::U_aniso_orth object'''
        #return super(Atom, self).u_aniso_orth()

    #@u_aniso_orth.setter
    #def u_aniso_orth(self, u_aniso):
        #if u_aniso is None:
            #from math import nan
            #self.set_u_aniso_orth(clipper_core.U_aniso_orth(*([nan]*6)))
        #else:
            #if type(u_aniso) == numpy.ndarray:
                #u_aniso = u_aniso.tolist()
            #self.set_u_aniso_orth(clipper_core.U_aniso_orth(*u_aniso))

    #@property
    #def is_null(self):
        #'''Check to see if this atom has been initialised'''
        #return super(Atom, self).is_null()

@mappedclass(clipper_core.Atom_list)
class Atom_list(clipper_core.Atom_list):
    '''
    Generate a Clipper Atom_list object from lists or numpy arrays of data.
    '''
    def __init__(self, elements, coords, occupancies, u_isos, u_anisos, allow_unknown_atoms = False):
        '''
        __init__(self, elements, coords, occupancies, u_isos, u_anisos)
            -> Atom_list

        Arguments:
            elements:    standard elemental abbreviations (e.g. "C", "Na1+" etc.)
                         See clipper.Atom.ATOM_NAMES for a full list of legal atom names
            coords:      (x,y,z) coordinates in Angstroms as an N x 3 array
            occupancies: fractional occupancies (one value per atom)
            u_isos:      isotropic B-factors (one value per atom)
            u_anisos:    the six unique entries in the anisotropic B-factor matrix
                         (u00, u11, u22, u01, u02, u12) as an N x 3 array
        '''
        # If true, an error will be raised if any element name is not in
        # the list of known scatterers.
        self.allow_unknown = allow_unknown_atoms
        clipper_core.Atom_list.__init__(self)
        self.extend_by(len(elements))
        self.elements = elements
        self.coord_orth = coords
        self.occupancies = occupancies
        self.u_isos = u_isos
        self.u_anisos = u_anisos


    @property
    def elements(self):
        '''Ordered list of all element names'''
        return super(Atom_list, self)._get_elements()

    @elements.setter
    def elements(self, elements):
        # Quick check to see if all element names are legal
        if not self.allow_unknown:
            if not set(elements).issubset(Atom.ATOM_NAMES):
                bad_atoms = []
                for el in set(elements):
                    if el not in Atom.ATOM_NAMES:
                        bad_atoms.append(el)
                bad_atoms = set(bad_atoms)
                errstring = '''
                    The following atom names are not recognised by Clipper:
                    {}
                    '''.format(bad_atoms)
                raise TypeError(errstring)
        super(Atom_list, self)._set_elements(elements)

    def _set_coord_orth(self, coords):
        n = len(self)
        array_in = numpy.empty((n, 3), numpy.double)
        array_in[:] = coords
        super(Atom_list, self)._set_coord_orth(array_in)

    def _get_coord_orth(self):
        '''Orthographic (x,y,z) coordinates of all atoms'''
        n = len(self)
        coords = numpy.empty((n,3,), numpy.double)
        super(Atom_list, self)._get_coord_orth(coords)
        return coords

    coord_orth = property(_get_coord_orth, _set_coord_orth)

    def _set_occupancies(self, occ):
        n = len(self)
        array_in = numpy.empty(n, numpy.double)
        array_in[:] = occ
        super(Atom_list, self)._set_occupancies(array_in)

    def _get_occupancies(self):
        '''Fractional occupancies of all atoms'''
        n = len(self)
        occ = numpy.empty(n, numpy.double)
        super(Atom_list, self)._get_occupancies(occ)
        return occ

    occupancies = property(_get_occupancies, _set_occupancies)

    def _set_u_isos(self, u_isos):
        n = len(self)
        array_in = numpy.empty(n, numpy.double)
        array_in[:] = u_isos
        super(Atom_list, self)._set_u_isos(array_in)

    def _get_u_isos(self):
        '''Isotropic B-factors of all atoms'''
        n = len(self)
        uiso = numpy.empty(n, numpy.double)
        super(Atom_list, self)._get_u_isos(uiso)
        return uiso

    u_isos = property(_get_u_isos, _set_u_isos)

    def _set_u_anisos(self, u_anisos):
        n = len(self)
        array_in = numpy.empty((n,6), numpy.double)
        array_in[:] = u_anisos
        super(Atom_list, self)._set_u_anisos(array_in)

    def _get_u_anisos(self):
        '''
        Anisotropic B-factor matrices for all atoms as an nx6 array, in the
        format: n*[u00, u11, u22, u01, u02, u12]. For purely isotropic
        atoms, set all elements in their row to math.nan or numpy.nan.
        '''
        n = len(self)
        uaniso = numpy.empty((n,6), numpy.double)
        super(Atom_list, self)._get_u_anisos(uaniso)
        return uaniso

    u_anisos = property(_get_u_anisos, _set_u_anisos)

    def get_minmax_grid(self, cell, grid_sampling):
        '''
        Get the minimum and maximum grid coordinates of a box that just
        encloses all atoms.
        Args:
            cell:
                A clipper.Cell object
            grid_sampling:
                A clipper.Grid_sampling object
        '''
        return super(Atom_list, self).get_minmax_grid(cell, grid_sampling)

########################################################################
# Real space coordinates
########################################################################

#@format_to_string
#@getters_to_properties(('xyz', '_get_xyz'))
#@mappedclass(clipper_core.Coord_orth)
#class Coord_orth(clipper_core.Coord_orth):
    #'''Coordinates in orthographic (x,y,z) space.'''
    #def __init__(self, xyz):
        #'''
        #__init__(self, xyz) -> Coord_orth

        #Args:
            #xyz ([float * 3]): (x, z, z) coordinates in Angstroms
        #'''
        #if isinstance(xyz, numpy.ndarray):
            ## Because SWIG does not correctly typemap numpy.float32
            #xyz = xyz.astype(float)
        #clipper_core.Coord_orth.__init__(self, *xyz)

    #def __add__(self, other):
        #if isinstance(other, Coord_orth):
            #return super(Coord_orth, self).__add__(other)
        #return self + Coord_orth(other)

    #def __radd__(self, other):
        #return self.__add__(other)

    #def __sub__(self, other):
        #if isinstance(other, Coord_orth):
            #return super(Coord_orth, self).__sub__(other)
        #return self - Coord_orth(other)

    #def __rsub__(self, other):
        #if isinstance(other, Coord_orth):
            #return other.__sub__(self)
        #return Coord_orth(other).__sub__(self)


#@format_to_string
#@getters_to_properties(('uvw','_get_uvw'))
#@mappedclass(clipper_core.Coord_grid)
#class Coord_grid(clipper_core.Coord_grid):
    #'''Integer grid coordinates in crystal space.'''
    #def __init__(self, uvw):
        #'''
        #__init__(self, uvw) -> Coord_grid

        #Args:
            #uvw ([int * 3]): (u, v, w) grid coordinates
        #'''
        ## Not sure why, but wrapped C functions fail when fed expansions
        ## of numpy int32 arrays, so we need to convert these to lists first.
        #if isinstance(uvw, numpy.ndarray):
            #uvw = uvw.tolist()
        #clipper_core.Coord_grid.__init__(self, *uvw)

    #def __add__(self, other):
        #if isinstance(other, Coord_grid):
            #return super(Coord_grid, self).__add__(other)
        #return self + Coord_grid(other)

    #def __radd__(self, other):
        #return self.__add__(other)

    #def __sub__(self, other):
        #if isinstance(other, Coord_grid):
            #return super(Coord_grid, self).__sub__(other)
        #return self - Coord_grid(other)

    #def __rsub__(self, other):
        #if isinstance(other, Coord_grid):
            #return other.__sub__(self)
        #return Coord_grid(other).__sub__(self)


# @format_to_string
# @getters_to_properties(('uvw', '_get_uvw'),'ceil','coord_grid','floor')
# @mappedclass(clipper_core.Coord_map)
# class Coord_map(clipper_core.Coord_map):
#     '''Like Coord_grid, but allowing non-integer values.'''
#     def __init__(self, uvw):
#         '''
#         __init__(self, uvw) -> Coord_map
#
#         Args:
#             uvw ([float * 3]): (u, v, w) grid coordinates
#         '''
#         if isinstance(uvw, numpy.ndarray):
#             # Because SWIG does not correctly typemap numpy.float32
#             uvw = uvw.astype(float)
#         clipper_core.Coord_map.__init__(self, *uvw)
#
#     def __add__(self, other):
#         if isinstance(other, Coord_map):
#             return super(Coord_map, self).__add__(other)
#         return self + Coord_map(other)
#
#     def __radd__(self, other):
#         return self.__add__(other)
#
#     def __sub__(self, other):
#         if isinstance(other, Coord_map):
#             return super(Coord_map, self).__sub__(other)
#         return self - Coord_map(other)
#
#     def __rsub__(self, other):
#         if isinstance(other, Coord_map):
#             return other.__sub__(self)
#         return Coord_map(other).__sub__(self)

# @format_to_string
# @getters_to_properties(('uvw', '_get_uvw'))
# @mappedclass(clipper_core.Coord_frac)
# class Coord_frac(clipper_core.Coord_frac):
#     '''
#     Fractional coordinates along unit cell axes (a,b,c), scaled to the
#     range 0..1 on each axis.
#     '''
#     def __init__(self, uvw):
#         '''
#         __init__(self, uvw) -> Coord_frac
#
#         Args:
#             uvw ([float * 3]): (u, v, w) fractional coordinates
#         '''
#         if isinstance(uvw, numpy.ndarray):
#             # Because SWIG does not correctly typemap numpy.float32
#             uvw = uvw.astype(float)
#         clipper_core.Coord_frac.__init__(self, *uvw)
#
#     def __add__(self, other):
#         if isinstance(other, Coord_frac):
#             return super(Coord_frac, self).__add__(other)
#         return self + Coord_frac(other)
#
#     def __radd__(self, other):
#         return self.__add__(other)
#
#     def __sub__(self, other):
#         if isinstance(other, Coord_frac):
#             return super(Coord_frac, self).__sub__(other)
#         return self - Coord_frac(other)
#
#     def __rsub__(self, other):
#         if isinstance(other, Coord_frac):
#             return other.__sub__(self)
#         return Coord_frac(other).__sub__(self)

########################################################################
# Reciprocal-space coordinates
########################################################################

@format_to_string
@getters_to_properties('us','vs','ws')
@mappedclass(clipper_core.Coord_reci_frac)
class Coord_reci_frac(clipper_core.Coord_reci_frac):
    '''
    Fractional hkl coordinates
    '''
    def __init__(self, uvw):
        clipper_core.Coord_reci_frac.__init__(self, *uvw)

@format_to_string
@getters_to_properties('xs','ys','zs')
@mappedclass(clipper_core.Coord_reci_orth)
class Coord_reci_orth(clipper_core.Coord_reci_orth):
    '''
    orthogonal reciprocal coordinate (length of which is invresolsq)
    '''
    def __init__(self, xyz_star):
        clipper_core.Coord_reci_orth.__init__(self, *xyz_star)


########################################################################
# FILE TYPES
########################################################################
@getters_to_properties('ccp4_spacegroup_number', 'cell', 'column_labels',
                       'column_paths','high_res_limit','history',
                       'hkl_sampling','low_res_limit','resolution',
                       'sort_order','spacegroup')
@mappedclass(clipper_core.CCP4MTZfile)
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
        '''
        __init__(self) -> CCP4MTZfile

        Create an empty CCP4MTZfile object. Call open_read(filename) to
        begin reading a MTZ file.
        '''
        clipper_core.CCP4MTZfile.__init__(self)

    @log_clipper
    def open_read(self, filename):
        '''Open an MTZ file for reading'''
        return super(CCP4MTZfile, self).open_read(filename)

    @log_clipper
    def import_hkl_data(self, cdata, mtzpath):
        '''
        Mark a set of data for import. NOTE: the data will not actually
        be imported until the close_read() method has been called.

        Args:
            cdata:   a Clipper.HKL_data_... object matched to the data types
                     in the mtzpath (e.g. HKL_data_F_phi for amplitudes and
                     phases).
            mtzpath: a string giving the path to the target columns within
                     the MTZ file. Available columns can be listed using
                     column_labels(). To import data spanning multiple
                     columns, use e.g.
                        "/crystal/dataset/[2FOFCWT, PH2FOFCWT]"
        '''
        return super(CCP4MTZfile, self).import_hkl_data(cdata, mtzpath)

    @property
    def spacegroup_confidence(self):
        return super(CCP4MTZfile, self).spacegroup_confidence()

    @spacegroup_confidence.setter
    def spacegroup_confidence(self, confidence):
        super(CCP4MTZfile, self).set_spacegroup_confidence(confidence)

    @property
    def title(self):
        return super(CCP4MTZfile, self).title()

    @title.setter
    def title(self, title):
        super(CCP4MTZfile, self).set_title(title)


@mappedclass(clipper_core.CIFfile)
class CIFfile(clipper_core.CIFfile):
    '''
    CIF import/export parent class for clipper objects.

    This is the import class which can be linked to an cif data
    file and be used to transfer data into a Clipper
    data structure.
    It is currently a read-only class.
    '''

    def __init__(self):
        '''
        __init__(self) -> CIFfile

        Create an empty CIFfile object. Call open_read to begin reading
        a .cif structure factor file.
        '''
        clipper_core.CIFfile.__init__(self)

    @log_clipper
    def open_read(self, filename):
        return super(CIFfile, self).open_read(filename)

    @log_clipper
    def resolution(self, cell):
        return super(CIFfile, self).resolution(cell)


@mappedclass(clipper_core.CCP4MAPfile)
class CCP4MAPfile(clipper_core.CCP4MAPfile):
    pass


########################################################################
# RECIPROCAL SPACE DATA TYPES
########################################################################
'''
Data types holding individual H,K,L data points, corresponding to the
array forms (HKL_data_...) below.
'''
@getters_to_properties('data_names','data_size','vals')
@mappedclass(clipper_core.ABCD_double)
class ABCD(clipper_core.ABCD_double):
    '''
    Hendrickson-Lattman coefficients
    '''
    def __init__(self, abcd):
        clipper_core.ABCD_double.__init__(self, *abcd)

@format_to_string
@getters_to_properties('h','k','l')
@mappedclass(clipper_core.HKL)
class HKL(clipper_core.HKL):
    pass

@getters_to_properties('cov', 'friedel', 'missing', 'sigE_mi', 'sigE_pl')
@mappedclass(clipper_core.E_sigE_double)
class E_sigE(clipper_core.E_sigE_double):
    pass

@getters_to_properties('a', 'b', 'friedel', 'missing', 'norm')
@mappedclass(clipper_core.F_phi_double)
class F_phi(clipper_core.F_phi_double):
    pass

@getters_to_properties('f', 'friedel', 'missing', 'sigf')
@mappedclass(clipper_core.F_sigF_ano_double)
class F_sigF_ano(clipper_core.F_sigF_ano_double):
    pass

@getters_to_properties('cov', 'f_mi', 'f_pl', 'friedel', 'missing', 'sigf_mi', 'sigf_pl')
@mappedclass(clipper_core.F_sigF_double)
class F_sigF(clipper_core.F_sigF_double):
    pass

@getters_to_properties('copy', 'friedel', 'missing')
@mappedclass(clipper_core.Flag)
class Flag(clipper_core.Flag):
    def get_flag(self):
        return super(Flag, self).get_flag()

    def set_flag(self, theFlag):
        return super(Flag, self).set_flag(theFlag)

    val = property(get_flag, set_flag)

@getters_to_properties('copy', 'friedel', 'missing')
@mappedclass(clipper_core.Flag_bool)
class Flag_bool(clipper_core.Flag_bool):
    def get_flag(self):
        return super(Flag, self).get_flag()

    def set_flag(self, theFlag):
        return super(Flag, self).set_flag(theFlag)

    val = property(get_flag, set_flag)

@getters_to_properties('cov', 'friedel', 'missing', 'sigI_mi', 'sigI_pl')
@mappedclass(clipper_core.I_sigI_double)
class I_sigI(clipper_core.I_sigI_double):
    pass

@getters_to_properties('friedel', 'missing')
@mappedclass(clipper_core.Phi_fom_double)
class Phi_fom(clipper_core.Phi_fom_double):
    pass




########################################################################
# HKL DATA ARRAY TYPES
########################################################################


@getters_to_properties('cell','spacegroup','resolution')
@mappedclass(clipper_core.HKL_info)
class HKL_info(clipper_core.HKL_info):
    def __init__(self):
        '''
        __init__(self) -> HKL_info

        Create an empty HKL_info object to store the vital statistics
        (h,k,l coordinates, cell parameters, spacegroup etc.) of a set
        of crystal data. It can be filled by passing it as an argument
        to (CCP4MTZfile or CIFfile).import_hkl_info() after loading a
        structure factor file.
        '''
        clipper_core.HKL_info.__init__(self)

@getters_to_properties('base_cell', 'base_hkl_info', 'cell', 'first', 'first_data',
                       'hkl_info', 'hkl_sampling', 'invresolsq_range',
                       'is_null', 'num_obs', 'resolution', 'spacegroup')
@mappedclass(clipper_core.HKL_data_ABCD_double)
class HKL_data_ABCD(clipper_core.HKL_data_ABCD_double):
    pass

@getters_to_properties('base_cell', 'base_hkl_info', 'cell', 'first', 'first_data',
                       'hkl_info', 'hkl_sampling', 'invresolsq_range',
                       'is_null', 'num_obs', 'resolution', 'spacegroup')
@mappedclass(clipper_core.HKL_data_E_sigE_double)
class HKL_data_E_sigE(clipper_core.HKL_data_E_sigE_double):
    pass

@getters_to_properties('base_cell', 'base_hkl_info', 'cell', 'first', 'first_data',
                       'hkl_info', 'hkl_sampling', 'invresolsq_range',
                       'is_null', 'num_obs', 'resolution', 'spacegroup')
@mappedclass(clipper_core.HKL_data_Flag)
class HKL_data_Flag(clipper_core.HKL_data_Flag):
    def __init__(self):
        '''
        __init__(self) -> HKL_data_Flag

        Create an empty object to store an array of integer flags
        (e.g. for holding of R-free flags). Fill it by passing it to
        (CCP4MTZfile or CIFfile).import_hkl_data together with the
        address of a suitable array.
        '''
        clipper_core.HKL_data_Flag.__init__(self)


@getters_to_properties('base_cell', 'base_hkl_info', 'cell', 'first', 'first_data',
                       'hkl_info', 'hkl_sampling', 'invresolsq_range',
                       'is_null', 'num_obs', 'resolution', 'spacegroup')
@mappedclass(clipper_core.HKL_data_Flag_bool)
class HKL_data_Flag_bool(clipper_core.HKL_data_Flag_bool):
    def __init__(self):
        '''
        __init__(self) -> HKL_data_Flag_bool

        Create an empty object to store an array of bool flags
        (e.g. for holding of R-free flags). Fill it by passing it to
        (CCP4MTZfile or CIFfile).import_hkl_data together with the
        address of a suitable array.
        '''
        clipper_core.HKL_data_Flag_bool.__init__(self)


@getters_to_properties('base_cell', 'base_hkl_info', 'cell', 'first', 'first_data',
                       'hkl_info', 'hkl_sampling', 'invresolsq_range',
                       'is_null', 'num_obs', 'resolution', 'spacegroup')
@mappedclass(clipper_core.HKL_data_F_sigF_double)
class HKL_data_F_sigF(clipper_core.HKL_data_F_sigF_double):
    def __init__(self):
        '''
        __init__(self) -> HKL_data_F_sigF

        Create an empty object to store an array of amplitudes and
        their associated sigmas. Fill it by passing it to
        (CCP4MTZfile or CIFfile).import_hkl_data together with the
        addresses of a suitable pair of arrays.
        '''
        clipper_core.HKL_data_F_sigF_double.__init__(self)

@getters_to_properties('base_cell', 'base_hkl_info', 'cell', 'first', 'first_data',
                       'hkl_info', 'hkl_sampling', 'invresolsq_range',
                       'is_null', 'num_obs', 'resolution', 'spacegroup')
@mappedclass(clipper_core.HKL_data_F_sigF_ano_double)
class HKL_data_F_sigF_ano(clipper_core.HKL_data_F_sigF_ano_double):
    pass


@getters_to_properties('base_cell', 'base_hkl_info', 'cell', 'first', 'first_data',
                       'hkl_info', 'hkl_sampling', 'invresolsq_range',
                       'is_null', 'num_obs', 'resolution', 'spacegroup')
@mappedclass(clipper_core.HKL_data_F_phi_double)
class HKL_data_F_phi (clipper_core.HKL_data_F_phi_double):
    def __init__(self):
        '''
        __init__(self) -> HKL_data_F_phi

        Create an empty object to store an array of amplitudes and their
        associated phases (observed or calculated). Fill it by passing it
        to (CCP4MTZfile or CIFfile).import_hkl_data together with the
        addresses of a suitable pair of arrays.
        '''
        clipper_core.HKL_data_F_phi_double.__init__(self)


@getters_to_properties('base_cell', 'base_hkl_info', 'cell', 'first', 'first_data',
                       'hkl_info', 'hkl_sampling', 'invresolsq_range',
                       'is_null', 'num_obs', 'resolution', 'spacegroup')
@mappedclass(clipper_core.HKL_data_I_sigI_double)
class HKL_data_I_sigI(clipper_core.HKL_data_I_sigI_double):
    def __init__(self):
        '''
        __init__(self) -> HKL_data_I_sigI

        Create an empty object to store an array of intensities and
        their associated sigmas. Fill it by passing it to
        (CCP4MTZfile or CIFfile).import_hkl_data together with the
        addresses of a suitable pair of arrays.
        '''
        clipper_core.HKL_data_I_sigI_double.__init__(self)

@getters_to_properties('base_cell', 'base_hkl_info', 'cell', 'first', 'first_data',
                       'hkl_info', 'hkl_sampling', 'invresolsq_range',
                       'is_null', 'num_obs', 'resolution', 'spacegroup')
@mappedclass(clipper_core.HKL_data_Phi_fom_double)
class HKL_data_Phi_fom(clipper_core.HKL_data_Phi_fom_double):
    pass


########################################################################
# CELLS, SYMMETRY, MAPS
########################################################################

@format_to_string
@getters_to_properties('nu', 'nv', 'nw', 'size')
@mappedclass(clipper_core.Grid)
class Grid(clipper_core.Grid):
    pass

@format_to_string
@mappedclass(clipper_core.Grid_range)
class Grid_range(clipper_core.Grid_range):
    pass

@getters_to_properties('dim')
@mappedclass(clipper_core.Grid_sampling)
class Grid_sampling(clipper_core.Grid_sampling):
    '''
    Object defining the grid used to sample points in a 3D map.
    '''
    def __init__(self, spacegroup, cell, resolution, rate = 1.5):
        '''
        __init__(self, spacegroup, cell, resolution, rate = 1.5) -> Grid_sampling

        Args:
            spacegroup (clipper.Spacegroup)
            cell (clipper.Cell)
            resolution (clipper.Resolution)
            rate: the over-sampling rate defining the spacing of grid
                  points. Leave this as the default unless you know what
                  you're doing.

        The first three arguments can be readily obtained from an HKL_info
        object: HKL_info.spacegroup, HKL_info.cell, HKL_info.resolution.

        '''
        clipper_core.Grid_sampling.__init__(self, spacegroup, cell,
                                            resolution, rate)

#@format_to_string
#@mappedclass(clipper_core.Cell_descr)
class Cell_descr(clipper_core.Cell_descr):
    def __init__(self, abc, angles):
        '''
        __init__(self, abc, angles) -> Cell_descr

        Args:
            abc ([float*3]): cell dimensions in Angstroms
            angles ([float*3]): alpha, beta and gamma angles in degrees
        '''
        clipper_core.Cell_descr.__init__(self, *abc, *angles)

#@format_to_string
#@getters_to_properties('cell_descr', 'matrix_frac', 'matrix_orth',
#                       'metric_real', 'metric_reci', 'volume')
#@mappedclass(clipper_core.Cell)
#class Cell(clipper_core.Cell):
    #'''
    #Define a crystallographic unit cell using the lengths of the three
    #sides a, b, c (in Angstroms) and the three angles alpha, beta, gamma
    #(in degrees).
    #'''
    #def __init__(self, abc, angles):
        #'''
        #__init__(self, abc, angles) -> Cell

        #Args:
            #abc ([float*3]): cell dimensions in Angstroms
            #angles ([float*3]): alpha, beta and gamma angles in degrees
        #'''
        #cell_descr = Cell_descr(abc, angles)
        ##cell_descr = clipper_core.Cell_descr(*abc, *angles)
        #clipper_core.Cell.__init__(self, cell_descr)

    #def __eq__(self, other):
        #return self.equals(other)

    #@property
    #def dim(self):
        #'''Returns the (a,b,c) lengths of the cell axes in Angstroms'''
        #return super(Cell, self).dim()

    #@property
    #def angles(self):
        #'''Returns the cell angles (alpha, beta, gamma) in radians'''
        #return super(Cell, self).angles()

    #@property
    #def angles_deg(self):
        #'''Returns the cell angles (alpha, beta, gamma) in degrees'''
        #return super(Cell, self).angles_deg()

    #@property
    #def recip_dim(self):
        #'''Returns the reciprocal cell lengths (a*, b*, c*) in inverse Angstroms'''
        #return super(Cell, self).recip_dim()

    #@property
    #def recip_angles(self):
        #'''Returns the reciprocal cell angles (alpha*, beta*, gamma*) in radians'''
        #return super(Cell, self).recip_angles()

    #@property
    #def recip_angles_deg(self):
        #'''Returns the reciprocal cell angles (alpha*, beta*, gamma*) in degrees'''
        #return super(Cell, self).recip_angles_deg()


@mappedclass(clipper_core.Unit_Cell)
class Unit_Cell(clipper_core.Unit_Cell):
    '''
    Atom-centric object holding the information necessary to construct a unit
    cell from a given atomic model of one asymmetric unit. Also provides
    functions returning the symmetry operations necessary to pack an
    arbitrary box in 3D space.
    '''
    def __init__(self, ref, atom_list, cell,
                 spacegroup, grid_sampling, padding = 0):
        '''
        __init__(self, ref, atom_list, cell, spacegroup, grid_sampling) -> Unit_Cell

        Note: internal calculations for finding symmetry equivalents are
        run using integer symmetry operations on grid coordinates for
        improved performance. This means that if you change the sampling
        rate of your maps, you will need to create a new Unit_Cell object
        to match.

        Args:
            ref ([float*3]):
                (x,y,z) coordinate of the reference you want to construct
                the unit cell relative to. Set this to the centroid of
                your atomic model for best results.
            atom_list (clipper.Atom_list):
                A Clipper Atom_list object containing your reference model.
            cell (clipper.Cell)
            spacegroup (clipper.Spacegroup)
            grid_sampling (clipper.Grid_Sampling)
            padding (int):
                An optional extra padding (in number of grid steps) to
                add around the reference model when defining the reference
                box for finding symmetry operations. In most cases this
                can be left as zero.
        '''
        ref_frac = Coord_orth(ref).coord_frac(cell)
        clipper_core.Unit_Cell.__init__(self, ref_frac, atom_list, cell,
                                        spacegroup, grid_sampling)

    @property
    def symops(self):
        '''Symops necessary to generate a unit cell from the model asu'''
        return super(Unit_Cell, self).symops()

    def all_symops_in_box(self, box_origin_xyz, box_size_uvw,
                          always_include_identity = False,
                          debug = False, sample_frequency = 2):
        '''
        Get an object defining all the symmetry operations mapping any
        part of your reference atomic model to a given box. Calculations
        are done in grid rather than Cartesian coordinates for performance
        reasons, so the precise shape of the box is dependent on the
        unit cell angles for the given spacegroup.

        Args:
            box_origin_xyz ([float*3]):
                The minimum (x,y,z) coordinates of the box. These will
                be rounded to the nearest grid coordinate.
            box_size_uvw ([int*3]):
                The dimensions of the box in grid coordinates
            always_include_identity(bool):
                If True, the identity symmetry operator will always be
                in the returned list.
        '''
        origin = numpy.empty(3)
        origin[:] = box_origin_xyz
        size = numpy.empty(3, numpy.int32)
        size[:] = box_size_uvw
        return super(Unit_Cell, self).all_symops_in_box(origin, size,
                    always_include_identity, debug, sample_frequency)

@format_to_string
@mappedclass(clipper_core.Isymop)
class Isymop(clipper_core.Isymop):
    pass

@mappedclass(clipper_core.Isymops)
class Isymops(clipper_core.Isymops):
    '''
    Python-friendly array class for holding lists of integer symmetry
    operators. For most purposes, you'll be better off using the fractional
    equivalent, Symops.
    '''
    pass

@mappedclass(clipper_core.RTop_frac)
@getters_to_properties('mat34', 'matrix', 'inverse', 'rotation',
                       'translation')
class RTop_frac(clipper_core.RTop_frac):
    '''
    Rotation/translation operator in fractional coordinates. Returned by
    various clipper functions. It is not generally advisable to create
    these for yourself unless you know exactly what you are doing.
    '''

    @property
    def format(self):
        '''
        Returns a symop-like string representation of this operator
        (e.g. "y, -x+y-1, z-5/6")
        '''
        return super(RTop_frac, self).format_as_symop()


    def __str__(self):
        return self.format

    def __hash__(self):
        '''
        Allows the object itself to be used as a key in dicts.
        '''
        return hash(self.format)

    def __eq__(self, other):
        return type(other) == RTop_frac and hash(self) == hash(other)

    def __mul__(self, coord):
        '''
        Apply this transformation operator to either a Python (u,v,w)
        iterable, or to a Coord_frac.
        '''
        if hasattr(coord, '__iter__'):
            return (self.__mul__(Coord_frac(coord))).uvw
        return super(RTop_frac, self).__mul__(coord)


@format_to_string
@getters_to_properties('matrix','mat34','rotation','translation')
@mappedclass(clipper_core.Symop)
class Symop(clipper_core.Symop):
    def __mul__(self, coord):
        '''
        Apply this transformation operator to either a Python (u,v,w)
        iterable, or to a Coord_orth.
        '''
        if hasattr(coord, '__iter__'):
            return (self.__mul__(Coord_frac(coord))).uvw
        return super(Symop, self).__mul__(coord)


#@mappedclass(clipper_core.Symops)
#class Symops(clipper_core.Symops):
    #'''
    #Python-friendly class for holding arrays of fractional symmetry
    #operators, with fast functions for obtaining the symmetry matrices
    #as numpy arrays, in fractional or orthogonal coordinates. Despite
    #the name, the stored elements are actually RTop_frac objects, since
    #Symop objects are not designed to store unit cell offsets.
    #'''

    #@property
    #def all_matrices34_frac(self):
        #'''
        #Get a (nx3x4) Numpy array summarising all the fractional symmetry
        #operators contained in this object.
        #Each entry in the array is of the form:

        #[[rot00 rot01 rot02 trn0]
         #[rot10 rot11 rot12 trn1]
         #[rot20 rot21 rot22 trn2]]
        #'''
        #n = len(self)
        #ret = numpy.empty([n,3,4],numpy.double)
        #super(Symops, self).all_matrices34_frac(ret)
        #return ret

    #def all_matrices34_orth(self, cell):
        #'''
        #Get a (nx3x4) Numpy array summarising all the orthographic symmetry
        #operators contained in this object. Requires a clipper.Cell as
        #an argument.
        #Each entry in the array is of the form:

        #[[rot00 rot01 rot02 trn0]
         #[rot10 rot11 rot12 trn1]
         #[rot20 rot21 rot22 trn2]]
        #'''
        #n = len(self)
        #ret = numpy.empty([n,3,4],numpy.double)
        #super(Symops, self).all_matrices34_orth(cell, ret)
        #return ret

    #@property
    #def all_matrices44_frac(self):
        #'''
        #Get a (nx4x4) Numpy array summarising all the fractional symmetry
        #operators contained in this object.
        #Each entry in the array is of the form:

        #[[rot00 rot01 rot02 trn0]
         #[rot10 rot11 rot12 trn1]
         #[rot20 rot21 rot22 trn2]
         #[  0     0     0    1  ]]
        #'''
        #n = len(self)
        #ret = numpy.empty([n,4,4],numpy.double)
        #super(Symops, self).all_matrices44_frac(ret)
        #return ret

    #def all_matrices44_orth(self, cell):
        #'''
        #Get a (nx4x4) Numpy array summarising all the orthographic symmetry
        #operators contained in this object. Requires a clipper.Cell as
        #an argument.
        #Each entry in the array is of the form:

        #[[rot00 rot01 rot02 trn0]
         #[rot10 rot11 rot12 trn1]
         #[rot20 rot21 rot22 trn2]
         #[  0     0     0    1  ]]
        #'''
        #n = len(self)
        #ret = numpy.empty([n,4,4],numpy.double)
        #super(Symops, self).all_matrices44_orth(cell, ret)
        #return ret



@getters_to_properties('matrix', 'mat34')
@mappedclass(clipper_core.RTop_orth)
class RTop_orth(clipper_core.RTop_orth):
    '''
    Rotation/translation operator in xyz coordinates. Returned by
    various clipper functions. It is not generally advisable to create
    these for yourself unless you know exactly what you are doing.
    '''
    def __mul__(self, coord):
        '''
        Apply this transformation operator to either a Python (x,y,z)
        iterable, or to a Coord_orth.
        '''
        if hasattr(coord, '__iter__'):
            return (self.__mul__(Coord_orth(coord))).xyz
        return super(RTop_orth, self).__mul__(coord)

@getters_to_properties('operator_grid_orth', 'operator_orth_grid', 'first',
                       'grid', 'is_null')
@mappedclass(clipper_core.NXmap_double)
class NXmap(clipper_core.NXmap_double):
    pass


@getters_to_properties(('grid','grid_sampling'), 'grid_sampling', 'cell',
                       'spacegroup', 'operator_grid_orth', 'operator_orth_grid',
                       'stats', 'voxel_size', 'voxel_size_frac', 'grid_asu',
                       'first', 'first_coord', 'is_null')
@mappedclass(clipper_core.Xmap_double)
class Xmap(clipper_core.Xmap_double):
    '''
    A Clipper crystallographic map generated from reciprocal space data.
    '''
    def __init__(self, spacegroup, cell, grid_sam, name = None, hkldata = None):
        '''
        __init__(self, spacegroup, cell, grid_sam) -> Xmap

        Generate a crystallographic map, and optionally fill it with
        data calculated from hkldata. If the hkldata argument is omitted,
        an empty map container will be created that can be filled with data
        later using the fft_from(fphi) method where fphi is a clipper.HKL_data_F_phi
        object.

        Args:
            spacegroup (clipper.Spacegroup)
            cell (clipper.Cell)
            grid_sam (clipper.Grid_sampling)
            name (string):
                Optionally, you can give the map a unique name for later
                identification
            hkldata (clipper.HKL_data_F_phi)
                Dataset from which the map will be generated
        '''
        clipper_core.Xmap_double.__init__(self, spacegroup, cell, grid_sam)

        # Some extra useful variables that aren't directly available from
        # the Clipper API
        self.name = name

        # Get the (nu, nv, nw) grid sampling as a numpy array
        self.grid_samples = self.grid.dim

        # Set this flag to True if you want this to be treated as a difference
        # map (i.e. visualisation at +/- 3 sigma, different handling in
        # simulations).
        self.is_difference_map = False

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

        if hkldata is not None:
            self.fft_from(hkldata)

    def recalculate_stats(self):
        '''
        Force recalculation of map statistics (max, min, mean, sigma,
        skewness and kurtosis).
        '''
        if self.is_null:
            raise RuntimeError('Map has no data!')
        self._min, self._max, self._mean, \
            self._sigma, self._skewness, self._kurtosis = self.stats


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

    def export_numpy(self):
        '''
        export the map asymmetric unit as a numpy array
        '''
        asu = self.grid_asu
        target = numpy.empty([asu.nu(), asu.nv(), asu.nw()],numpy.double)
        super(Xmap, self).export_numpy(target, 'C')
        return target

########################################################################
# ANALYSIS
########################################################################

@getters_to_properties('invresolsq_limit', 'is_null', 'limit')
@mappedclass(clipper_core.Resolution)
class Resolution(clipper_core.Resolution):
    pass


########################################################################
# VECTORS AND MATRICES
########################################################################

#@getters_to_properties('det', 'inverse', 'is_null', 'transpose', 'as_numpy')
#@format_to_string
#@mappedclass(clipper_core.mat33_float)
#class Mat33_float(clipper_core.mat33_float):
    #pass

########################################################################
# UTILITIES
########################################################################

@mappedclass(clipper_core.Util)
class Util(clipper_core.Util):
    '''
    Utility functions. Every method in this class should be static, and
    the class should never need to be instantiated.
    '''
    @staticmethod
    def get_minmax_grid(coords, cell, grid_sampling):
        '''
        Returns a numpy array containing the minimum and maximum grid
        coordinates bounding the rhombohedral box encompassing the
        given coordinates for the given cell and grid sampling.
        Args:
            coords:
                an n*3 numpy array of (x,y,z) coordinates
            cell:
                a clipper.Cell object
            grid_sampling:
                a clipper.Grid_sampling object
        '''
        return super(Util, Util).get_minmax_grid(coords, cell, grid_sampling)

########################################################################
# TEST FUNCTIONS
########################################################################
@log_clipper
def test_log_warn():
    clipper_core.warn_test()

@log_clipper
def test_log_except():
    clipper_core.except_test()

@log_clipper
def test_log_no_warn():
    pass
