# Load a MolProbity data set and return a SciPy RegularGridInterpolator 
# object for later fast interpolation of values
def generate_interpolator(filename, wrap_axes = True):
	from scipy.interpolate import RegularGridInterpolator
	import numpy
	infile = open(filename, 'r')
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
	data = numpy.loadtxt(filename)
	
	# Convert each coordinate to an integral number of steps along each
	# axis
	axes_ij = []
	
	for i in range(ndim):
		ss = step_size[i]
		fs = first_step[i]
		lb = lower_bound[i]
		axis_vals = data[:,i]
		axes_ij.append(((axis_vals - ss/2 - lb) / ss).astype(int))
	
	full_grid[axes_ij] = data[:,ndim]
	
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
			
	# Finally, convert the axes to radians
	from math import pi
	axis = numpy.array(axis)
	axis = axis/180*pi
	
	return RegularGridInterpolator(axis, full_grid, bounds_error = True)
		
		
		
	
	
			
		
	
	
	
		
