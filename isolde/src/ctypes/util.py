# @Author: Tristan Croll
# @Date:   13-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   Tristan Croll
# @Last modified time: 18-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



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
