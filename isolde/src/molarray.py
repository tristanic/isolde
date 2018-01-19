
from numpy import uint8, int32, uint32, float64, float32, uintp, byte, bool as npy_bool, integer, empty, array
from chimerax.core.atomic.molc import string, cptr, pyobject, set_cvec_pointer, pointer, size_t
from chimerax.core.atomic.molarray import Collection
from . import molobject
from .molobject import c_function, c_array_function, cvec_property
#from .molobject import object_map
from .molobject import Proper_Dihedral
import ctypes

from chimerax.core.atomic import Atom, Atoms, Residue, Residues

from chimerax.core.atomic.molarray import _atoms, _atoms_or_nones, \
        _non_null_atoms, _bonds, _pseudobond_groups, _pseudobonds, \
        _elements, _residues, _non_null_residues, _chains, \
        _non_null_chains, _atomic_structures, structure_datas, \
        _atoms_pair, _pseudobond_group_map

def _proper_dihedrals(p):
    return Proper_Dihedrals(p)
def _proper_dihedrals_or_nones(p):
    return [Proper_Dihedral.c_ptr_to_py_inst(ptr) if ptr else None for ptr in p]
def _non_null_proper_dihedrals(p):
    return Proper_Dihedrals(p[p!=0])



_proper_dihedrals_from_atoms = c_array_function('proper_dihedral_from_atoms', args=(cptr, cptr),
                            per_object=False)

def proper_dihedrals_from_atoms(atoms, residues, name):
    n = len(atoms)//4
    pointers = empty(n, cptr)
    names = empty(n, string)
    names[:]=name
    _proper_dihedrals_from_atoms(atoms._c_pointers, len(atoms), pointer(names), residues._c_pointers, pointer(pointers))
    d = Proper_Dihedrals(pointers)
    #d.names = name
    return d

class Proper_Dihedrals(Collection):
    
    def __init__(self, d_pointers = None):
        Collection.__init__(self, d_pointers, Proper_Dihedral, Proper_Dihedrals)
    
    names = cvec_property('proper_dihedral_name', string, read_only=True)
    angles = cvec_property('proper_dihedral_angle', float32, read_only=True)
    residues = cvec_property('proper_dihedral_residue', cptr, astype=_residues, read_only=True,
        doc='Returns a :class:`Residues` giving the parent residue of each dihedral. Read only')
