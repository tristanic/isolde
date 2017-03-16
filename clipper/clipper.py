from .lib import clipper_python_core as _clipper
import numpy

_clipper_messages = _clipper.MessageStreamSingleton().clipper_messages

# Classes directly passed through from core
#from .lib.clipper_python_core import Coord_orth, Coord_frac, Coord_grid, \
#    Coord_map, Coord_reci_frac, Coord_reci_orth
from .lib.clipper_python_core import Coord_grid
from .lib.clipper_python_core import CCP4MTZfile, CIFfile, CCP4MAPfile
from .lib.clipper_python_core import ABCD_double as ABCD, HKL, E_sigE_double as E_sigE,\
                     F_phi_double as F_phi, F_sigF_ano_double as F_sigF_ano,\
                     F_sigF_double as F_sigF, Flag, Flag_bool, \
                     I_sigI_double as I_sigI, Phi_fom_double as Phi_fom
from .lib.clipper_python_core import HKL_info, HKL_data_ABCD_double as HKL_data_ABCD, \
                     HKL_data_E_sigE_double as HKL_data_E_sigE, \
                     HKL_data_Flag, HKL_data_Flag_bool, \
                     HKL_data_F_sigF_double as HKL_data_F_sigF, \
                     HKL_data_F_sigF_ano_double as HKL_data_F_sigF_ano, \
                     HKL_data_F_phi_double as HKL_data_F_phi, \
                     HKL_data_I_sigI_double as HKL_data_I_sigI, \
                     HKL_data_Phi_fom_double as HKL_data_Phi_fom
from .lib.clipper_python_core import Grid, Grid_range, Grid_sampling, Cell_descr, \
                     Isymop, RTop_frac, Symop, Symops, RTop_orth, \
                     Resolution
from .lib.clipper_python_core import Util
from .lib.clipper_python_core import warn_test, except_test



'''
Sub-classing of specific classes to provide extra functionality.
NOTE: The extra functions will only be available in classes that you 
create directly, not in those created internally by clipper. As such,
you should generally limit extensions here to new __init__ functions 
that would be difficult to implement in SWIG. Anything else should be
added to the SWIG wrapper itself.
'''
class Unit_Cell(_clipper.Unit_Cell):
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
        _clipper.Unit_Cell.__init__(self, ref_frac, atom_list, cell,
                                        spacegroup, grid_sampling)


class Coord_orth(_clipper.Coord_orth):
    def __init__(self, coords):
        try:
            _clipper.Coord_orth.__init__(self, *coords)
        except:
            if type(coords) == numpy.ndarray:
                _clipper.Coord_orth.__init__(self, *(coords.astype(float)))
            else:
                raise

class Coord_frac(_clipper.Coord_frac):
    def __init__(self, coords):
        try:
            _clipper.Coord_frac.__init__(self, *coords)
        except:
            if type(coords) == numpy.ndarray:
                _clipper.Coord_frac.__init__(self, *(coords.astype(float)))
            else:
                raise

class Coord_reci_orth(_clipper.Coord_reci_orth):
    def __init__(self, coords):
        try:
            _clipper.Coord_reci_orth.__init__(self, *coords)
        except:
            if type(coords) == numpy.ndarray:
                _clipper.Coord_reci_orth.__init__(self, *(coords.astype(float)))
            else:
                raise

class Coord_reci_frac(_clipper.Coord_frac):
    def __init__(self, coords):
        try:
            _clipper.Coord_reci_frac.__init__(self, *coords)
        except:
            if type(coords) == numpy.ndarray:
                _clipper.Coord_reci_frac.__init__(self, *(coords.astype(float)))
            else:
                raise

#class Coord_grid(_clipper.Coord_grid):
    #def __init__(self, coords):
        #try:
            #_clipper.Coord_grid.__init__(self, *coords)
        #except:
            #if type(coords) == numpy.ndarray:
                #_clipper.Coord_grid.__init__(self, *(coords.tolist()))
            #else:
                #raise

class Coord_map(_clipper.Coord_map):
    def __init__(self, coords):
        try:
            _clipper.Coord_map.__init__(self, *coords)
        except:
            if type(coords) == numpy.ndarray:
                _clipper.Coord_map.__init__(self, *(coords.astype(float)))
            else:
                raise

      
class Atom(_clipper.Atom):
    def __init__(self, element, coord, occ, u_iso, u_aniso = None, 
                 allow_unknown = False):
        '''
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
            allow_unknown (bool):
                If false, only element names from the list of known
                scatterers in ATOM_NAMES will be allowed.
        '''
        _clipper.Atom.__init__(self)
        self.allow_unknown = allow_unknown
        self.element = element
        self.coord = coord
        self.occupancy = occ
        self.u_iso = u_iso
        self.u_aniso = u_aniso

class Atom_list(_clipper.Atom_list):
    def __init__(self, elements, coords, occupancies, u_isos, u_anisos,
                 allow_unknown_atoms = False):
        _clipper.Atom_list.__init__(self)        
        self.allow_unknown = allow_unknown_atoms
        self.extend_by(len(elements))
        self.elements = elements
        self.coord_orth = coords
        self.occupancies = occupancies
        self.u_isos = u_isos
        self.u_anisos = u_anisos
        
class Xmap(_clipper.Xmap_double):
    def __init__(self, spacegroup, cell, grid_sampling,
                 name = None, hkldata = None, is_difference_map = None):
        _clipper.Xmap_double.__init__(self, spacegroup, cell, grid_sampling)
        self.name = name
        self.is_difference_map = is_difference_map
        if hkldata is not None:
            self.fft_from(hkldata)
                 
