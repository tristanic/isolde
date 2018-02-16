import os, sys, glob
import numpy
import ctypes
from chimerax.core.state import State
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

from numpy import int8, uint8, int32, uint32, float64, float32, byte, bool as npy_bool

libdir = os.path.dirname(os.path.abspath(__file__))
libfile = glob.glob(os.path.join(libdir, '..', 'openmm.cpython*'))[0]

_c_functions = CFunctions(os.path.splitext(libfile)[0])
c_property = _c_functions.c_property
cvec_property = _c_functions.cvec_property
c_function = _c_functions.c_function
c_array_function = _c_functions.c_array_function

class OpenMM_Thread_Handler:
    def __init__(self, integrator, c_pointer=None):
        cname = type(self).__name__.lower()
        if c_pointer is None:
            new_func = cname + '_new'
            int_ptr = int(integrator.this)
            f = c_function(new_func, args=(ctypes.c_void_p,), ret=ctypes.c_void_p)
            c_pointer = f(int_ptr)
        set_c_pointer(self, c_pointer)
        f = c_function('set_'+cname+'_py_instance', args=(ctypes.c_void_p, ctypes.py_object))
        f(self._c_pointer, self)

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def delete(self):
        c_function('openmm_thread_handler_delete', args=(ctypes.c_void_p,))(self._c_pointer)
        delattr(self, '_c_pointer')

    def step(self, steps):
        f = c_function('openmm_thread_handler_step',
            args=(ctypes.c_void_p, ctypes.c_size_t))
        f(self._c_pointer, steps)

    @property
    def thread_finished(self):
        f = c_function('openmm_thread_handler_thread_finished',
            args=(ctypes.c_void_p,),
            ret=npy_bool)
        return f(self._c_pointer)

    def finalize_thread(self):
        f = c_function('openmm_thread_handler_finalize_thread',
            args=(ctypes.c_void_p,))
        f(self._c_pointer)
