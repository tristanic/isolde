
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

c_double_p = ctypes.POINTER(ctypes.c_double)
c_size_t_p = ctypes.POINTER(ctypes.c_size_t)

class RegularGridInterpolator:
    _new_interp = c_function('rg_interp_new', 
        args=(ctypes.c_size_t, c_size_t_p, c_double_p, c_double_p, c_double_p), ret = ctypes.c_void_p)
    _interpolate = c_function('rg_interpolate',
        args=(ctypes.c_void_p, c_double_p, ctypes.c_size_t, c_double_p))
    _delete = c_function('rg_interp_delete', args=(ctypes.c_void_p,))
    _dim = c_function('rg_interp_dim', args=(ctypes.c_void_p, ), ret=ctypes.c_size_t)
    _min = c_function('rg_interp_min', args=(ctypes.c_void_p, c_double_p))
    _max = c_function('rg_interp_max', args=(ctypes.c_void_p, c_double_p))
    _values = c_function('rg_interp_values', args=(ctypes.c_void_p, c_double_p))
    _axis_lengths = c_function('rg_interp_lengths', args=(ctypes.c_void_p, c_size_t_p))
    def __init__(self, dim, axis_lengths, min_vals, max_vals, grid_data):
        '''
        A C++ implementation of n-dimensional regular grid interpolation,
        interfaced to Python using ctypes. About 2-3 times faster than 
        the SciPy RegularGridInterpolator, and more compatible with 
        threading.
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
        
        self._c_pointer = self._new_interp(dim, axis_lengths.ctypes.data_as(c_size_t_p), 
            min_vals.ctypes.data_as(c_double_p), max_vals.ctypes.data_as(c_double_p), 
            grid_data.ctypes.data_as(c_double_p))
        
        #~ self._dim = dim
        self._min_vals = min_vals
        self._max_vals = max_vals
        
    @property
    def dim(self):
        return self._dim(self._c_pointer)
    
    @property
    def min(self):
        ret = numpy.empty(self.dim)
        self._min(self._c_pointer, ret.ctypes.data_as(c_double_p))
        return ret
        
    @property
    def max(self):
        ret = numpy.empty(self.dim)
        self._max(self._c_pointer, ret.ctypes.data_as(c_double_p))
        return ret
    
    @property
    def axis_lengths(self):
        ret = numpy.empty(self.dim, numpy.int)
        self._axis_lengths(self._c_pointer, ret.ctypes.data_as(c_size_t_p))
        return ret
    
    @property
    def values(self):
        dim = self.dim
        lengths = self.axis_lengths
        ret = numpy.empty(lengths, numpy.double)
        self._values(self._c_pointer, ret.ctypes.data_as(c_double_p))
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
            raise TypeError('Wrong number of dimensions! This is a {}-'\
                           +'dimensional interpolator.'.format(self.dim))
        #~ if len(data) == 0:
            #~ return
        if data.dtype != numpy.double or not data.flags.c_contiguous:
            in_data = numpy.empty(data.shape, numpy.double)
            in_data[:] = data
        else:
            in_data = data
        
        n = in_data.shape[0]
        ret = numpy.empty(n, dtype=numpy.double)
        self._interpolate(self._c_pointer, in_data.ctypes.data_as(c_double_p),
                          n, ret.ctypes.data_as(c_double_p))
        return ret
    
    def __call__(self, data):
        return self.interpolate(data)
    
    def __del__(self):
        self._delete(self._c_pointer)


def test_interpolator():
    import numpy
    dim=3
    axis_lengths=numpy.array([74,74,74], numpy.uintp)
    mins = numpy.zeros(3)
    maxs = numpy.ones(3)
    data = numpy.random.rand(74,74,74)
    from scipy.interpolate import RegularGridInterpolator as ScipyInterp
    axis = numpy.array(range(74))/73
    scrg = ScipyInterp((axis,axis,axis), data)
    
    rg = RegularGridInterpolator(dim, axis_lengths, mins, maxs, data)
    #test_data = numpy.random.rand(n,3)
    return (scrg, rg)
    
