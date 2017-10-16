from .lib import clipper_python_core as _clipper

# Classes directly passed through from core
from _clipper import Coord_orth, Coord_frac, Coord_grid, Coord_map, \
                     Coord_reci_frac, Coord_reci_orth
from _clipper import CCP4MTZfile, CIFfile, CCP4MAPfile
from _clipper import ABCD_double as ABCD, HKL, E_sigE_double as E_sigE,\
                     F_phi_double as F_phi, F_sigF_ano_double as F_sigF_ano,\
                     F_sigF_double as F_sigF, Flag, Flag_bool, \
                     I_sigI_double as I_sigI, Phi_fom_double as Phi_fom
from _clipper import HKL_info, HKL_data_ABCD_double as HKL_data_ABCD, \
                     HKL_data_E_sigE_double as HKL_data_E_sigE, \
                     HKL_data_Flag, HKL_data_Flag_bool, \
                     HKL_data_F_sigF_double as HKL_data_F_sigF, \
                     HKL_data_F_sigF_ano_double as HKL_data_F_sigF_ano, \
                     HKL_data_F_phi_double as HKL_data_F_phi, \
                     HKL_data_I_sigI_double as HKL_data_I_sigI, \
                     HKL_data_Phi_fom_double as HKL_data_Phi_fom
from _clipper import Grid, Grid_range, Grid_sampling, Cell_descr, \
                     Isymop, RTop_frac, Symop, Symops, RTop_orth, \
                     Resolution
from _clipper import Util
from _clipper import warn_test, except_test



'''
Sub-classing of specific classes to provide extra functionality.
NOTE: The extra functions will only be available in classes that you 
create directly, not in those created internally by clipper. As such,
you should generally limit extensions here to new __init__ functions 
that would be difficult to implement in SWIG. Anything else should be
added to the SWIG wrapper itself.
'''

      
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
        super(Atom, self).__init__(self)
        self.allow_unknown = allow_unknown
        self.element = element
        self.coord = coord
        self.occupancy = occ
        self.u_iso = u_iso
        self.u_aniso_orth = u_aniso
        
