# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



import os
import numpy
from math import degrees, radians, pi
from .. import geometry
from ..constants import defaults
#from .dihedrals import Dihedrals

from time import time

package_directory = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(package_directory, 'molprobity_data')

def get_molprobity_data_dir():
    return DATA_DIR

def generate_scipy_interpolator(file_prefix, wrap_axes = True):
    '''
    (Deprecated - use generate_interpolator() instead)
    Load a MolProbity data set and return a SciPy RegularGridInterpolator
    object.
    '''
    from scipy.interpolate import RegularGridInterpolator
    import numpy, pickle
    infile = None
    # First try to load from a pickle file
    try:
        infile = open(file_prefix+'.pickle', 'rb')
        axis, full_grid = pickle.load(infile)
        infile.close()
        infile = None
    except:
        # If pickle load fails for any reason, fall back to loading from
        # text, then regenerate the pickle file at the end.
        if infile is not None:
            infile.close()
        infile = open(file_prefix+'.data', 'rt', encoding='utf-8')
        # Throw away the first line - we don't need it
        infile.readline()
        # Get number of dimensions
        ndim = int(infile.readline().split()[-1])
        # Throw away the next line - it's just headers
        infile.readline()
        lower_bound = []
        upper_bound= []
        number_of_bins = []

        step_size = []
        first_step = []
        last_step = []
        axis = []

        # Read in the header to get the dimensions and step size for each
        # axis, and initialise the axis arrays
        for i in range(ndim):
            line = infile.readline().split()
            lb = float(line[2])
            lower_bound.append(lb)
            ub = float(line[3])
            upper_bound.append(ub)
            nb = int(line[4])
            number_of_bins.append(nb)

            ss = (ub - lb)/nb
            step_size.append(ss)
            # Values are at the midpoint of each bin
            fs = lb + ss/2
            first_step.append(fs)
            ls = ub - ss/2
            last_step.append(ls)
            axis.append(numpy.linspace(fs,ls,nb))

        infile.close()

        full_grid = numpy.zeros(number_of_bins)

        # Slurp in the actual numerical data as a numpy array
        data = numpy.loadtxt(file_prefix+'.data')

        # Convert each coordinate to an integral number of steps along each
        # axis
        axes = []
        for i in range(ndim):
            axes.append([])

        for i in range(ndim):
            ss = step_size[i]
            fs = first_step[i]
            lb = lower_bound[i]
            axis_vals = data[:,i]
            axes[i]=(((axis_vals - ss/2 - lb) / ss).astype(int))

        full_grid[axes] = data[:,ndim]

        # At this point we should have the full n-dimensional matrix, with
        # all values not present in the text file present as zeros.
        # Now we have to consider periodicity. Since we're just going to
        # be doing linear interpretation, the easiest approach is to simply
        # pad the array on all sides with the values from the opposite extreme
        # of the relevant matrix. This can be handily done with numpy.pad
        if wrap_axes:
            full_grid = numpy.pad(full_grid, 1, 'wrap')

            # ... and we need to extend each of the axes by one step to match
            for i, a in enumerate(axis):
                fs = first_step[i]
                ls = last_step[i]
                ss = step_size[i]
                a = numpy.pad(a, 1, mode='constant')
                a[0] = fs - ss
                a[-1] = ls + ss
                axis[i] = a

        # Replace all zero or negative values with the minimum positive non-zero
        # value, so that we can use logs
        full_grid[full_grid<=0] = numpy.min(full_grid[full_grid > 0])
        # Finally, convert the axes to radians
        from math import pi
        axis = numpy.array(axis)
        axis = axis/180*pi

        # Pickle a tuple containing the axes and grid for fast loading in
        # future runs

        outfile = open(file_prefix+'.pickle', 'w+b')
        pickle.dump((axis, full_grid), outfile)
        outfile.close()

    return RegularGridInterpolator(axis, full_grid, bounds_error = True)


def generate_interpolator(file_prefix, wrap_axes = True):
    from .interpolation.interp import RegularGridInterpolator
    ndim, axis_lengths, min_vals, max_vals, grid_data = generate_interpolator_data(file_prefix, wrap_axes)
    return RegularGridInterpolator(ndim, axis_lengths, min_vals, max_vals, grid_data)

def generate_interpolator_data(file_prefix, wrap_axes = True, regenerate = None):
    '''
    Load a MolProbity data set and format it into a form ready for generation of
    a RegularGridInterpolator object for later fast interpolation of values.
    '''
    import numpy, pickle
    infile = None
    if regenerate is None:
        with open(os.path.join(DATA_DIR, 'regenerate'), 'rt') as f:
            regenerate = int(f.readline()[0])

    if regenerate:
        print('Regenerating contour pickle files')
        import glob
        [os.remove(f) for f in glob.glob(os.path.join(DATA_DIR, '*.pickle'))]
        with open(os.path.join(DATA_DIR, 'regenerate'), 'wt') as f:
            f.write('0')

    # First try to load from a pickle file
    try:
        infile = open(file_prefix+'.pickle', 'r+b')
        ndim, axis_lengths, min_vals, max_vals, grid_data = pickle.load(infile)
        infile.close()
        infile = None
    except:
        # If pickle load fails for any reason, fall back to loading from
        # text, then regenerate the pickle file at the end.
        if infile is not None:
            infile.close()
        infile = open(file_prefix+'.data', 'r')
        # Throw away the first line - we don't need it
        infile.readline()
        # Get number of dimensions
        ndim = int(infile.readline().split()[-1])
        # Throw away the next line - it's just headers
        infile.readline()
        lower_bounds = []
        upper_bounds= []
        axis_lengths = []

        step_sizes = []
        min_vals = []
        max_vals = []
        #~ axes = []

        # Read in the header to get the dimensions and step size for each
        # axis, and initialise the axis arrays
        for i in range(ndim):
            line = infile.readline().split()
            lb = float(line[2])
            lower_bounds.append(lb)
            ub = float(line[3])
            upper_bounds.append(ub)
            nb = int(line[4])
            axis_lengths.append(nb)

            ss = (ub - lb)/nb
            step_sizes.append(ss)
            # Values are at the midpoint of each bin
            fs = lb + ss/2
            min_vals.append(fs)
            ls = ub - ss/2
            max_vals.append(ls)
            #~ axis.append(numpy.linspace(fs,ls,nb))

        infile.close()

        grid_data = numpy.zeros(axis_lengths)

        # Slurp in the actual numerical data as a numpy array
        data = numpy.loadtxt(file_prefix+'.data')

        # Convert each coordinate to an integral number of steps along each
        # axis
        axes = []
        for i in range(ndim):
            axes.append([])

        for i in range(ndim):
            ss = step_sizes[i]
            fs = min_vals[i]
            lb = lower_bounds[i]
            axis_vals = data[:,i]
            axes[i]=(((axis_vals - ss/2 - lb) / ss).astype(int))

        grid_data[axes] = data[:,ndim]

        '''
        At this point we should have the full n-dimensional matrix, with
        all values not present in the text file present as zeros.
        Now we have to consider periodicity. Since we're just going to
        be doing linear interpretation, the easiest approach is to simply
        pad the array on all sides with the values from the opposite extreme
        of the relevant matrix. This can be handily done with numpy.pad
        '''
        if wrap_axes:
            grid_data = numpy.pad(grid_data, 1, 'wrap')

            # ... and we need to extend each of the axes by one step to match
            for i, a in enumerate(axes):
                min_v = min_vals[i]
                max_v = max_vals[i]
                ss = step_sizes[i]
                min_vals[i] = min_v - ss
                max_vals[i] = max_v + ss
                axis_lengths[i] += 2
        # Replace all zero or negative values with the minimum positive non-zero
        # value, so that we can use logs
        grid_data[grid_data<=0] = numpy.min(grid_data[grid_data > 0])
        # Finally, convert the axes to radians
        min_vals = numpy.radians(numpy.array(min_vals, numpy.double))
        max_vals = numpy.radians(numpy.array(max_vals, numpy.double))

        # Pickle a tuple containing the axes and grid for fast loading in
        # future runs

        outfile = open(file_prefix+'.pickle', 'w+b')
        axis_lengths = numpy.array(axis_lengths, numpy.int)


        pickle.dump((ndim, axis_lengths, min_vals, max_vals, grid_data), outfile)
        #~ pickle.dump((axis, full_grid), outfile)
        outfile.close()
    return (ndim, axis_lengths, min_vals, max_vals, grid_data)

from ..param_mgr import Param_Mgr, autodoc, param_properties
from .constants import validation_defaults as _val_defaults
@param_properties
@autodoc
class Validation_Params(Param_Mgr):
    _default_params = {
        'track_ramachandran_status':    (_val_defaults.TRACK_RAMACHANDRAN_STATUS, None),
        'rounds_per_rama_update':       (_val_defaults.ROUNDS_PER_RAMA_UPDATE, None),
        'track_rotamer_status':         (_val_defaults.TRACK_ROTAMER_STATUS, None),
        'rounds_per_rota_update':       (_val_defaults.ROUNDS_PER_ROTA_UPDATE, None),
    }
