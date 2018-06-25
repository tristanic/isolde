# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



import os, sys, glob
from chimerax.atomic import molc
# from chimerax.atomic.molc import CFunctions, string, cptr, pyobject, \
#     set_c_pointer, pointer, size_t

CFunctions = molc.CFunctions
string = molc.string
cptr = molc.cptr
pyobject = molc.pyobject
set_c_pointer = molc.set_c_pointer
pointer = molc.pointer
size_t = molc.size_t

from numpy import uint8, int32, uint32, float64, float32, byte, bool as npy_bool
import ctypes
import numpy

from ..ctypes import convert_and_sanitize_numpy_array
from ..util import compiled_lib_extension
libdir = os.path.dirname(os.path.abspath(__file__))
libfile = os.path.join(libdir, '..', 'lib_nd_interp.'+compiled_lib_extension())



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
    A C++ implementation of n-dimensional regular grid interpolation, interfaced
    to Python using ctypes. Designed as an almost-drop-in replacement for the
    SciPy RegularGridInterpolator, but significantly faster (particularly for
    small numbers of interpolations). Whereas the SciPy interpolator has a fixed
    overhead of about 400 microseconds per call, this implementation reduces
    that to 20 microseconds. In addition, the speed of each interpolation is
    increased approximately 5-fold.

    Achieving this speed does impose a few limitations, however. Whereas the
    SciPy interpolator allows uneven steps along each axis, for this
    implementation the gridded data *must* be evenly-spaced (although it is not
    necessary for each axis to have the same spacing). The initialisation
    function is therefore slightly different to the SciPy implementation to
    reflect this.

    While in theory any number of dimensions is supported, in practice memory
    limitations make it impractical beyond 4-5 dimensions for most situations.
    For example, to handle interpolation over a 360 range with periodic wrapping
    on a 10 degree grid requires 39 points per axis. For three dimensions the
    resulting grid contains about 60k points; four dimensions contains about
    2.3 million; five dimensions 90 million and six dimensions 3.5 billion (that
    is, about 25 GB of RAM in double precision).

    Example usage is as follows (also see :func:`test_interpolator`):

    .. code:: python

        import numpy
        dim = 3
        axis_lengths = [50, 50, 50]
        min_vals = [0, 0, 0]
        max_vals = [1, 1, 1]
        grid_data = numpy.random.rand(50,50,50)
        interp = RegularGridInterpolator(dim, axis_lengths, min_vals, max_vals,
                                         grid_data)
        data = numpy.random.rand(10000, 3)

        results = interp(data)
        # or, if you prefer:
        results = interp.interpolate(data)

    Interpolation is fastest if the input is a NumPy double array - anything
    else is converted internally (with an associated performance penalty) prior
    to calling the C++ function.
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
        Prepare the interpolator for a given n-dimensional grid. Once created,
        the grid elements cannot be modified.

        Args:
            * dim:
                - An integer specifying the number of dimensions
            * axis_lengths:
                - A numpy int array of length dim, specifying the number of
                  points along each axis
            * min_vals:
                - A numpy double array of length dim, specifying the minimum
                  value of each axis
            * max_vals:
                  A numpy double array of length dim, specifying the maximum
                  value of each axis
            * grid_data:
                A n-dimensional numpy double array with dimensions matching
                axis_lengths, containing all the gridded data.
        '''
        axis_lengths = convert_and_sanitize_numpy_array(axis_lengths, numpy.uint32)
        min_vals = convert_and_sanitize_numpy_array(min_vals, NPY_FLOAT)
        max_vals = convert_and_sanitize_numpy_array(max_vals, NPY_FLOAT)
        grid_data = convert_and_sanitize_numpy_array(grid_data, NPY_FLOAT)

        self._c_pointer = self._new_interp(dim, pointer(axis_lengths),
            pointer(min_vals), pointer(max_vals), pointer(grid_data))

        #~ self._dim = dim
        self._min_vals = min_vals
        self._max_vals = max_vals

    @property
    def dim(self):
        '''
        Returns the number of dimensions of the interpolator. Read only.
        '''
        return self._dim(self._c_pointer)

    @property
    def min(self):
        '''
        Returns a NumPy array of length :attr:`dim` providing the minimum limit
        of each axis. Read only.
        '''
        ret = numpy.empty(self.dim, NPY_FLOAT)
        self._min(self._c_pointer, pointer(ret))
        return ret

    @property
    def max(self):
        '''
        Returns a NumPy array of length :attr:`dim` providing the maximum limit
        of each axis. Read only.
        '''
        ret = numpy.empty(self.dim, NPY_FLOAT)
        self._max(self._c_pointer, pointer(ret))
        return ret

    @property
    def axis_lengths(self):
        '''
        Returns a NumPy array of length :attr:`dim` providing the number of
        steps on each axis. Read only.
        '''
        ret = numpy.empty(self.dim, numpy.uint32)
        self._axis_lengths(self._c_pointer, pointer(ret))
        return ret

    @property
    def values(self):
        '''
        Returns the original gridded data values as a :attr:`dim`-dimensional
        NumPy array
        '''
        lengths = self.axis_lengths
        ret = numpy.empty(lengths, NPY_FLOAT)
        self._values(self._c_pointer, pointer(ret))
        return ret

    @property
    def grid(self):
        '''
        Returns a :attr:`dim`-length tuple of NumPy arrays giving the axis
        values at each grid step.
        '''
        grid = []
        lengths = self.axis_lengths
        min_vals = self.min
        max_vals = self.max
        for length, minv, maxv in zip(lengths, min_vals, max_vals):
            grid.append(numpy.linspace(minv, maxv, length))
        return tuple(grid)


    def interpolate(self, data):
        '''
        Returns the interpolated values for a set of (x(1), x(2), ...
        x(:attr:`dim`)) points.

        Args:
            * data:
                - a (n, :attr:`dim`) 2D NumPy array providing the coordinates at
                  which to calculate the interpolated values. For fastest
                  performance the data type of the array should be
                  :attr:`numpy.double`.

        Returns:
            * a 1D NumPy double array
        '''
        if data.shape[1] != self.dim:
            raise TypeError('Wrong number of dimensions! This is a '\
                           +'{}-dimensional interpolator.'.format(self.dim))
        in_data = convert_and_sanitize_numpy_array(data, NPY_FLOAT)
        n = len(in_data)
        ret = numpy.empty(n, dtype=NPY_FLOAT)
        self._interpolate(self._c_pointer, pointer(in_data), n, pointer(ret))
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
    '''
    Creates a :class:`scipy.interpolate.RegularGridInterpolator` and a
    :class:`isolde.RegularGridInterpolator` using identical random data for
    side-by-side comparison. Each dimension will have a random number of steps
    between 10 and 50.

    Args:
        * n:
            - number of dimensions

    Returns:
        * a (SciPy interpolator, ISOLDE interpolator) tuple
    '''
    import numpy
    dim=n
    axis_lengths=numpy.random.randint(10, 50, size=n, dtype=numpy.uintp)
    mins = numpy.zeros(n)
    maxs = numpy.ones(n)
    data = numpy.random.rand(*axis_lengths)
    from scipy.interpolate import RegularGridInterpolator as ScipyInterp
    axes = [numpy.array(range(l))/(l-1) for l in axis_lengths]
    #~ axis = numpy.array(range(36))/35
    scrg = ScipyInterp(axes, data)

    rg = RegularGridInterpolator(dim, axis_lengths, mins, maxs, data)
    #test_data = numpy.random.rand(n,3)
    return (scrg, rg)
