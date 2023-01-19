# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



import numpy

def convert_and_sanitize_numpy_array(array, dtype):
    '''
    Ensure a numpy array has the specified data type and is C-contiguous in memory.
    '''
    if isinstance(array, numpy.ndarray) and array.dtype == dtype \
                                        and array.data.c_contiguous:
        return array
    ret = numpy.empty(array.shape, dtype)
    ret[:] = array
    return ret
