import numpy
import sys, os, glob
import ctypes

dpath = os.path.dirname(os.path.abspath(__file__))
libfile = glob.glob(os.path.join(dpath, '_symmetry.cpython*'))[0]
_symmetry = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), libfile))


_sym_transforms = _symmetry.atom_and_bond_sym_transforms
_sym_transforms.argtypes = (ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_double),
    ctypes.c_size_t, ctypes.POINTER(ctypes.c_double), ctypes.c_double)
_sym_transforms.restype = ctypes.py_object

def sym_transforms_in_sphere(atoms, transforms, center, cutoff):
    natoms = len(atoms);
    tf = numpy.empty(transforms.shape, numpy.double)
    tf[:] = transforms
    c = numpy.empty(3, numpy.double)
    c[:] = center
    n_tf = len(transforms);
    result = _sym_transforms(atoms._c_pointers, natoms, transforms.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        n_tf, c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), cutoff)
    from chimerax.core.atomic.molarray import _atoms, _bonds
    from chimerax.core.geometry import Places
    natoms = len(result[0])
    atom_coords = result[1].reshape((natoms,3))
    nbonds = len(result[3])
    bond_positions = Places(opengl_array=result[4].reshape((nbonds*2,4,4)))
    return (_atoms(result[0]), atom_coords, result[2], _bonds(result[3]), bond_positions, result[5])
