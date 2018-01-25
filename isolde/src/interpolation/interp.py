
import os, sys, glob
from chimerax.core.atomic.molc import CFunctions, string, cptr, pyobject, \
    set_c_pointer, pointer, size_t
from numpy import uint8, int32, uint32, float64, float32, byte, bool as npy_bool
import ctypes
import numpy

libdir = os.path.dirname(os.path.abspath(__file__))
libfile = glob.glob(os.path.join(libdir, '..', '_nd_interp.cpython*'))[0]



_c_functions = CFunctions(os.path.splitext(libfile)[0])
c_property = _c_functions.c_property
cvec_property = _c_functions.cvec_property
c_function = _c_functions.c_function
c_array_function = _c_functions.c_array_function

NPY_FLOAT = numpy.double
FLOAT_TYPE = ctypes.c_double
SIZE_TYPE = ctypes.c_size_t

UINT32_TYPE = ctypes.c_uint32

#c_double_p = ctypes.POINTER(ctypes.c_double)
C_FLOAT_P = ctypes.POINTER(FLOAT_TYPE)
C_UINT32_P = ctypes.POINTER(UINT32_TYPE)

class RegularGridInterpolator:
    '''
    A C++ implementation of n-dimensional regular grid interpolation,
    interfaced to Python using ctypes. About 5 times faster than 
    the SciPy RegularGridInterpolator for 3D data, and more compatible with 
    threading.
    '''
    
    _new_interp = c_function('rg_interp_new', 
        args=(SIZE_TYPE, C_UINT32_P, C_FLOAT_P, C_FLOAT_P, C_FLOAT_P), ret = ctypes.c_void_p)
    _interpolate = c_function('rg_interpolate',
        args=(ctypes.c_void_p, C_FLOAT_P, SIZE_TYPE, C_FLOAT_P))
    _delete = c_function('rg_interp_delete', args=(ctypes.c_void_p,))
    _dim = c_function('rg_interp_dim', args=(ctypes.c_void_p, ), ret=SIZE_TYPE)
    _min = c_function('rg_interp_min', args=(ctypes.c_void_p, C_FLOAT_P))
    _max = c_function('rg_interp_max', args=(ctypes.c_void_p, C_FLOAT_P))
    _values = c_function('rg_interp_values', args=(ctypes.c_void_p, C_FLOAT_P))
    _axis_lengths = c_function('rg_interp_lengths', args=(ctypes.c_void_p, C_UINT32_P))
    _copy = c_function('rg_interp_copy', args=(ctypes.c_void_p, ))
    def __init__(self, dim, axis_lengths, min_vals, max_vals, grid_data):
        '''
        Args:
            dim:
                An integer specifying the number of dimensions
            axis_lengths:
                A numpy int array of length dim, specifying the number of
                points along each axis
            min_vals:
                A numpy float array of length dim, specifying the minimum
                value of each axis
            max_vals:
                A numpy float array of length dim, specifying the maximum
                value of each axis
            grid_data:
                A n-dimensional numpy float array of the given dimensions,
                containing all the gridded data.
        '''
        axis_lengths = convert_and_sanitize_numpy_array(axis_lengths, numpy.uint32)
        min_vals = convert_and_sanitize_numpy_array(min_vals, NPY_FLOAT)
        max_vals = convert_and_sanitize_numpy_array(max_vals, NPY_FLOAT)
        grid_data = convert_and_sanitize_numpy_array(grid_data, NPY_FLOAT)
        
        self._c_pointer = self._new_interp(dim, axis_lengths.ctypes.data_as(C_UINT32_P), 
            min_vals.ctypes.data_as(C_FLOAT_P), max_vals.ctypes.data_as(C_FLOAT_P), 
            grid_data.ctypes.data_as(C_FLOAT_P))
        
        #~ self._dim = dim
        self._min_vals = min_vals
        self._max_vals = max_vals
        
    @property
    def dim(self):
        return self._dim(self._c_pointer)
    
    @property
    def min(self):
        ret = numpy.empty(self.dim, NPY_FLOAT)
        self._min(self._c_pointer, ret.ctypes.data_as(C_FLOAT_P))
        return ret
        
    @property
    def max(self):
        ret = numpy.empty(self.dim, NPY_FLOAT)
        self._max(self._c_pointer, ret.ctypes.data_as(C_FLOAT_P))
        return ret
    
    @property
    def axis_lengths(self):
        ret = numpy.empty(self.dim, numpy.uint32)
        self._axis_lengths(self._c_pointer, ret.ctypes.data_as(C_UINT32_P))
        return ret
    
    @property
    def values(self):
        dim = self.dim
        lengths = self.axis_lengths
        ret = numpy.empty(lengths, NPY_FLOAT)
        self._values(self._c_pointer, ret.ctypes.data_as(C_FLOAT_P))
        return ret
            
    @property
    def grid(self):
        grid = []
        lengths = self.axis_lengths
        min_vals = self.min
        max_vals = self.max
        for length, minv, maxv in zip(lengths, min_vals, max_vals):
            grid.append(numpy.linspace(minv, maxv, length))
        return tuple(grid)
        
    
    def interpolate(self, data):
        if data.shape[1] != self.dim:
            raise TypeError('Wrong number of dimensions! This is a '\
                           +'{}-dimensional interpolator.'.format(self.dim))
        if data.dtype != NPY_FLOAT or not data.flags.c_contiguous:
            in_data = numpy.empty(data.shape, NPY_FLOAT)
            in_data[:] = data
        else:
            in_data = data
        
        n = in_data.shape[0]
        ret = numpy.empty(n, dtype=NPY_FLOAT)
        self._interpolate(self._c_pointer, in_data.ctypes.data_as(C_FLOAT_P),
                          n, ret.ctypes.data_as(C_FLOAT_P))
        return ret
    
    def __call__(self, data):
        return self.interpolate(data)
    
    def __del__(self):
        self._delete(self._c_pointer)


    def __deepcopy__(self, memo):
        from copy import deepcopy
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        setattr(result, '_c_pointer', self._copy(self._c_pointer))
        return result


def test_interpolator(n):
    import numpy
    dim=n
    axis_lengths=numpy.array([36]*n, numpy.uintp)
    mins = numpy.zeros(n)
    maxs = numpy.ones(n)
    data = numpy.random.rand(*[36]*n)
    from scipy.interpolate import RegularGridInterpolator as ScipyInterp
    axis = numpy.array(range(36))/35
    scrg = ScipyInterp([axis]*n, data)
    
    rg = RegularGridInterpolator(dim, axis_lengths, mins, maxs, data)
    #test_data = numpy.random.rand(n,3)
    return (scrg, rg)
    
def convert_and_sanitize_numpy_array(array, dtype):
    '''
    Convert a numpy array to the specified data type, and ensure its
    contents are C-contiguous in memory.
    '''
    #~ if array.flags.c_contiguous:
        #~ if array.dtype == dtype:
            #~ return array
        #~ return array.as_type(dtype)
    ret = numpy.empty(array.shape, dtype)
    ret[:] = array
    return ret
