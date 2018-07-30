import os, sys, glob
import numpy
import ctypes

from chimerax.atomic import molc

CFunctions = molc.CFunctions
string = molc.string
cptr = molc.cptr
pyobject = molc.pyobject
set_c_pointer = molc.set_c_pointer
pointer = molc.pointer
size_t = molc.size_t

from chimerax.atomic import ctypes_support as convert

from numpy import int8, uint8, int32, uint32, float64, float32, byte, bool as npy_bool

from .util import compiled_lib_extension
libdir = os.path.dirname(os.path.abspath(__file__))
libfile = os.path.join(libdir, 'libclipper.'+compiled_lib_extension())

_c_functions = CFunctions(os.path.splitext(libfile)[0])
c_property = _c_functions.c_property
cvec_property = _c_functions.cvec_property
c_function = _c_functions.c_function
c_array_function = _c_functions.c_array_function

def _clipper_init():
    c_function('clipper_init', args=())()


class CCP4MTZFile:
    '''
    Handles reading/writing of MTZ files.
    '''
    def __init__(self):
        '''
        Just creates the object.
        '''
        self._c_pointer = c_function('mtzfile_new', args=(), ret=(ctypes.c_void_p))()
        self._mode = None

    def __enter__(self):
        return self

    def delete(self):
        c_function("mtzfile_delete", args=(ctypes.c_void_p,))(self._c_pointer)
        self._c_pointer = None

    def __exit__(self):

        self.delete()

    def open_read(self, filename):
        '''
        Open an existing MTZ file for reading.
        '''
        f = c_function('mtzfile_open_read',
            args=(ctypes.c_void_p, ctypes.c_char_p))
        f(self._c_pointer, filename.encode('utf-8'))

    def close_read(self):
        '''
        Close the MTZ file when done
        '''
        c_function('mtzfile_close_read',args=(c_types.c_void_p))(self._c_pointer)

    def open_write(self, filename):
        '''
        Open or create an MTZ file for writing
        '''
        f = c_function('mtzfile_open_write',
            args=(ctypes.c_void_p, ctypes.c_char_p))
        f(self._c_pointer, filename.encode('utf-8'))

    def close_write(self):
        '''
        Close the MTZ file when done
        '''
        c_function('mtzfile_close_write',args=(c_types.c_void_p))(self._c_pointer)