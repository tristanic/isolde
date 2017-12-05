
import os, sys, glob
import numpy
from chimerax.core.state import State
from chimerax.core.atomic import molc
from chimerax.core.atomic.molc import CFunctions, string, cptr, pyobject, \
    set_c_pointer, pointer, size_t
# object map lookups
from chimerax.core.atomic.molobject import _py_atoms, _py_bonds, _atoms, \
                _atom_pair, _atom_or_none, _bonds, _chain, _element, \
                _pseudobonds, _residue, _residues, _rings, _non_null_residues, \
                _residue_or_none, _residues_or_nones, _residues_or_nones, \
                _chains, _atomic_structure, _pseudobond_group, \
                _pseudobond_group_map
# object map functions
from chimerax.core.atomic.molobject import object_map, add_to_object_map

from numpy import uint8, int32, uint32, float64, float32, byte, bool as npy_bool

libdir = os.path.dirname(os.path.abspath(__file__))
libfile = glob.glob(os.path.join(libdir, 'molc.cpython*'))[0]

_c_functions = CFunctions(os.path.splitext(libfile)[0])
c_property = _c_functions.c_property
cvec_property = _c_functions.cvec_property
c_function = _c_functions.c_function
c_array_function = _c_functions.c_array_function


def _dihedrals(p):
    from .molarray import Dihedrals
    return Dihedrals(p)
def _dihedral_or_none(p):
    return object_map(p, Dihedral) if p else None

class Dihedral(State):
    def __init__(self, c_pointer):
        set_c_pointer(self, c_pointer)
    
    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')
    
    def __str__(self):
        return self.name
    
    def reset_state(self):
        pass
    
    name = c_property('dihedral_name', string, doc = 'Name of this dihedral.')
    angle = c_property('dihedral_angle', float32, read_only=True, doc = 'Angle in radians. Read only.')
        
