import numpy

def convert_and_sanitize_numpy_array(array, dtype):
    '''
    Ensure a numpy array has the specified data type and is C-contiguous in memory.
    '''
    if array.dtype == dtype and array.data.c_contiguous:
        return array
    ret = numpy.empty(array.shape, dtype)
    ret[:] = array
    return ret
