
import os, sys, glob
import numpy
import ctypes

from chimerax.core.atomic import molc
from chimerax.core.atomic.molc import CFunctions, string, cptr, pyobject, \
    set_c_pointer, pointer, size_t
# object map lookups
from chimerax.core.atomic.molobject import _atoms, \
                _atom_pair, _atom_or_none, _bonds, _chain, _element, \
                _pseudobonds, _residue, _residues, _rings, _non_null_residues, \
                _residue_or_none, _residues_or_nones, _residues_or_nones, \
                _chains, _atomic_structure, _pseudobond_group, \
                _pseudobond_group_map

from numpy import uint8, int32, uint32, float64, float32, byte, bool as npy_bool

from .molobject import Proper_Dihedral_Mgr


libdir = os.path.dirname(os.path.abspath(__file__))
libfile = glob.glob(os.path.join(libdir, 'molc.cpython*'))[0]

_c_functions = CFunctions(os.path.splitext(libfile)[0])
c_property = _c_functions.c_property
cvec_property = _c_functions.cvec_property
c_function = _c_functions.c_function
c_array_function = _c_functions.c_array_function


def dihedral_extensions():
    '''
    Extend ChimeraX AtomicStructure, Residue and Residues classes to 
    include ISOLDE dihedrals
    '''
    from chimerax.core.atomic import AtomicStructure
    
    setattr(AtomicStructure, "_proper_dihedral_mgr", None)
    
    def get_dihedral_mgr(self):
        if self._proper_dihedral_mgr is None:
            self._proper_dihedral_mgr = Proper_Dihedral_Mgr(self)
    
    setattr(AtomicStructure, "proper_dihedral_mgr", property(get_dihedral_mgr))
    
    from chimerax.core.atomic import Residues
    def named_dihedrals(self, name):
        us = self.unique_structures
        for s in us:
            pass
    
    
    
