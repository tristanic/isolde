import os
from math import degrees, radians
package_directory = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(package_directory, 'molprobity_data')

# Load a MolProbity data set and return a SciPy RegularGridInterpolator 
# object for later fast interpolation of values
def generate_interpolator(file_prefix, wrap_axes = True):
    from scipy.interpolate import RegularGridInterpolator
    import numpy, pickle
    infile = None
    
    # First try to load from a pickle file
    try:
        infile = open(file_prefix+'.pickle', 'r+b')
        axis, full_grid = pickle.load(infile)
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

cis_min = radians(-30)
cis_max = radians(30)
trans_min = radians(150)
trans_max = radians(-150)

    
def omega_type(omega):
    if omega is None:
        return None
    if omega >= trans_min or omega <= trans_max:
        return "trans"
    elif omega >= cis_min and omega <= cis_max:
        return "cis"
    return "twisted"
        
class RamaValidator():
    import numpy
    from math import pi
    # Define cases.
    # name: human-readable name for display purposes.
    # file_name: path to the file containing the contour data.
    # cutoffs: [outlier cutoff, 1.0, allowed cutoff].
    # log_cutoffs: cutoffs on a log scale to define color thresholds
    # color_scale: will hold a ThreeColorScale object to generate colours
    #               for display based on scores and the defined cutoffs.
    # validator: handle to hold the RegularGridInterpolator object generated
    #            at runtime.
    
    # Ordered list of keys for the below cases, since dict objects do not
    # guarantee any particular order.
    case_keys = ('General', 'Glycine', 'IleVal', 'PrePro', 'CisPro', 'TransPro')
    
    cases = {}
    cases['CisPro'] = {
        'name': 'Cis proline',
        'file_prefix': os.path.join(data_dir, 'rama8000-cispro'),
        'cutoffs': [0.002, 1.0, 0.02],
        'log_cutoffs': numpy.log(numpy.array([0.002, 1.0, 0.02])),
        'color_scale':  None,
        'validator': None
        }
    cases['TransPro'] = {
        'name': 'Trans proline',
        'file_prefix': os.path.join(data_dir, 'rama8000-transpro'),
        'cutoffs': [0.001, 1.0, 0.02],
        'log_cutoffs': numpy.log(numpy.array([0.001, 1.0, 0.02])),
        'color_scale':  None,
        'validator': None
        }
    cases['Glycine'] = {
        'name': 'Glycine',
        'file_prefix': os.path.join(data_dir, 'rama8000-gly-sym'),
        'cutoffs': [0.001, 1.0, 0.02],
        'log_cutoffs': numpy.log(numpy.array([0.001, 1.0, 0.02])),
        'color_scale':  None,
        'validator': None
        }
    cases['PrePro'] = {
        'name': 'Preceding proline',
        'file_prefix': os.path.join(data_dir, 'rama8000-prepro-noGP'),
        'cutoffs': [0.001, 1.0, 0.02],
        'log_cutoffs': numpy.log(numpy.array([0.001, 1.0, 0.02])),
        'color_scale':  None,
        'validator': None
        }
    cases['IleVal'] = {
        'name': 'Isoleucine or valine',
        'file_prefix': os.path.join(data_dir, 'rama8000-ileval-nopreP'),
        'cutoffs': [0.001, 1.0, 0.02],
        'log_cutoffs': numpy.log(numpy.array([0.001, 1.0, 0.02])),
        'color_scale':  None,
        'validator': None
        }
    cases['General'] = {
        'name': 'General',
        'file_prefix': os.path.join(data_dir, 'rama8000-general-noGPIVpreP'),
        'cutoffs': [0.0005, 1.0, 0.02],
        'log_cutoffs': numpy.log(numpy.array([0.0005, 1.0, 0.02])),
        'color_scale':  None,
        'validator': None
        }
                
    
    def __init__(self, color_scale = 'PiYG'):
        for key, case in self.cases.items():
            case['validator'] = generate_interpolator(case['file_prefix'])
        self.current_structure = None
        self.rama_scores = None
        self.rama_types = None
        # Locally cached array of current phi and psi for plotting purposes
        self.phipsi = None
        # Locally cached list of residues for looking up from the plot
        self.residues = None
        self.case_arrays = {
            'CisPro': [],
            'TransPro': [],
            'Glycine': [],
            'PrePro': [],
            'IleVal': [],
            'General': []
            }
        self._proline_indices = []
        self.color_scale = color_scale
        
        from . import color
        for key, case in self.cases.items():
            minval, maxval, midval = case['log_cutoffs']
            case['color_scale'] = color.standard_three_color_scale(
                    self.color_scale, minval, maxval, midval)
        
        # Array of color values corresponding to scores
        self.current_colors = [];
        
        
        

    def rama_bins(scores, types):
        score_bins = []
        for score, rt in scores, types:
            if rt is None:
                score_bins.append(None)
                continue
            outlier, favored, allowed = self.cases[rt]['cutoffs']
            if score >= allowed:
                score_bins.append('favored')
            elif score >= outlier:
                score_bins.append('allowed')
            else:
                score_bins.append('outlier')




    # Initialise the validator with a sequence and matching set of initial
    # phi, psi and omega values. A value of None in either of the phi or psi 
    # arrays indicates a chain end/break, which should not be included in
    # the analysis. The residues will be sorted into their various MolProbity
    # cases here.
    def load_structure(self, residues, resnames, phi, psi, omega):
        ores, ophi, opsi, oomega = [residues, phi, psi, omega]
        import numpy
        phipsi = numpy.column_stack([phi,psi])
        # Filter out any residues missing either phi or psi
        none_filter = self.none_filter = numpy.where(
            numpy.invert(numpy.any(numpy.isnan(phipsi), axis = 1)))
        phipsi = self.phipsi = phipsi[none_filter]
        residues = self.residues = numpy.array(residues)[none_filter]
        resnames = self.resnames = resnames[none_filter]
        phi = numpy.array(phi)[none_filter]
        psi = numpy.array(psi)[none_filter]
        omega = numpy.array(omega)[none_filter]

        cases = self.cases
        ca = self.case_arrays
        # Reset the arrays of indices for each case
        for key in ca:
            ca[key] = []
        self._proline_indices = []
        num_all_residues = len(ores)
        num_residues = len(resnames)
        # Terminal residues will retain a score of -1
        self.rama_scores = numpy.array([-1] * num_residues, dtype = 'float32')
        self.rama_types = [None] * num_residues
        self.output_colors = numpy.array([128,128,128,255] * num_all_residues, dtype='ubyte').reshape(num_all_residues,4)
        self.current_colors = self.output_colors[none_filter]


        for i, resname in enumerate(resnames):
            if phi[i] == None or psi[i] == None:
                continue
        
            # Determine the category of the residue to match it with a set
            # of contours. There is a hierarchy to the categories: Pro beats
            # PrePro for example.
            
            if resname == 'PRO':
                # we need to sort the prolines into cis and trans every
                # update, so for convenience we'll store a complete list
                # here.
                self._proline_indices.append(i)
                if omega_type(omega[i]) == 'cis':
                    self.rama_types[i] = 'CisPro'
                    ca['CisPro'].append(i)
                else:
                    self.rama_types[i] = 'TransPro'
                    ca['TransPro'].append(i)
            elif resname == 'GLY':
                self.rama_types[i] = 'Glycine'
                ca['Glycine'].append(i)
            elif i < num_residues-1 and resnames[i+1] == 'PRO':
                self.rama_types[i] = 'PrePro'
                ca['PrePro'].append(i)
            elif resname in ('ILE', 'VAL'):
                self.rama_types[i] = 'IleVal'
                ca['IleVal'].append(i)
            else:
                self.rama_types[i] = 'General'
                ca['General'].append(i)
        # Calculate the current scores
        self.update(ophi,opsi,oomega)
        
    def reset(self):
        '''
        Release the current structure and reset phi, psi etc. to None.
        '''
        self.current_structure = None
        self.rama_scores = None
        self.rama_types = None
        # Locally cached array of current phi and psi for plotting purposes
        self.phipsi = None
        self.case_arrays = {
            'CisPro': [],
            'TransPro': [],
            'Glycine': [],
            'PrePro': [],
            'IleVal': [],
            'General': []
            }
        self.current_colors = []
            
        
    
    def update(self, phi, psi, omega, return_colors = True):
        '''
        Calculate and return Ramachandran scores for the current dihedral
        values
        '''
        # Since prolines may change from cis to trans or vice versa
        # in the course of a simulation, we need to double-check these
        # each time.
        import numpy
        phipsi = numpy.column_stack([phi,psi])
        # since None values crash the interpolator, we'll re-cast them to
        # a value outside the range of the maps
        phipsi = self.phipsi = phipsi[self.none_filter]
        ca = self.case_arrays
        ca['CisPro'] = []
        ca['TransPro'] = []
        for i in self._proline_indices:
            if omega_type(omega[i]) == 'cis':
                self.rama_types[i] = 'CisPro'
                ca['CisPro'].append(i)
            elif omega_type(omega[i]) == 'trans':
                self.rama_types[i] = 'TransPro'
                ca['TransPro'].append(i)
        for key, indices in self.case_arrays.items():
            indices = numpy.array(indices, dtype = 'int')
            case = self.cases[key]
            v = case['validator']
            scores = v(phipsi[indices]).astype(float)
            self.rama_scores[indices] = scores
            if return_colors:
                if len(scores):
                    c = case['color_scale']
                    colors = c.get_colors(numpy.log(scores))
                    self.current_colors[indices] = colors
                
        if return_colors:
            self.output_colors[self.none_filter] = self.current_colors
            return self.rama_scores, self.output_colors
        
        return self.rama_scores
        
class RamaPlot():
    
    def __init__(self, session, container, validator):
        import numpy
        self.session = session
        self.container = container
        self.validator = validator
        self.current_case = None
        
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_qt5agg import (
            FigureCanvasQTAgg as FigureCanvas,
            NavigationToolbar2QT as NavigationToolbar
            )
        
        fig = self.figure = Figure()
        
        axes = self.axes = fig.add_subplot(111)
        
        axes.set_xticks([-120,-60,0,60,120])
        axes.set_yticks([-120,-60,0,60,120])
        #axes.set_xticklabels(axes.get_xticklabels(), rotation=60)
        axes.set_xlim(-180,180)
        axes.set_ylim(-180,180)
        axes.minorticks_on()
        axes.autoscale(enable=False)
        
        # Scatter plot needs to always have at least one point, so we'll
        # define a point off the scale for when there's nothing to plot
        self.default_coords = numpy.array([200,200])
        self.default_logscores = numpy.ones(1)
        self.scatter = None
        
        canvas = self.canvas = FigureCanvas(fig)
        self.resize_cid = self.canvas.mpl_connect('resize_event', self.on_resize)
        container.addWidget(canvas)
        
        self.contours = {}
        self.change_case('General')
            
    def cache_contour_plots(self, key):
        import numpy
        case_data = self.validator.cases[key]['validator']
        grid = numpy.degrees(case_data.grid)
        from operator import itemgetter
        contours = itemgetter(0,2)(self.validator.cases[key]['cutoffs'])
        grid = numpy.degrees(case_data.grid)
        values = numpy.rot90(numpy.fliplr(case_data.values)).astype(float)
        logvalues = self.logvalues = numpy.log(values)
        contour_plot = self.axes.contour(*grid, values, contours)
        pcolor_plot = self.axes.pcolor(*grid, logvalues, cmap = 'BuGn')
        for coll in contour_plot.collections:
           coll.remove()
        pcolor_plot.remove()
        self.contours[key] = (contour_plot, pcolor_plot)
        
        
    def on_resize(self, *_):
        axes = self.axes
        axes.set_xlim(-180,180)
        axes.set_ylim(-180,180)
        axes.set_xticks([-120,-60,0,60,120])
        axes.set_yticks([-120,-60,0,60,120])
        axes.autoscale(enable=False)
        if self.scatter.is_figure_set():
            self.scatter.remove()
        self.figure.set_facecolor('0.5')  
        self.canvas.draw()
        self.background = self.canvas.copy_from_bbox(self.axes.bbox)
        self.axes.add_artist(self.scatter)
        self.canvas.draw()
        self.update_scatter()
    
    def on_pick(self, event):
        ind = event.ind[0]
        res_index = self.validator.case_arrays[self.current_case][ind]
        picked_residue = self.validator.residues[res_index]
        from . import view
        view.focus_on_selection(self.session, self.session.main_view, picked_residue.atoms)
        
    def change_case(self, case_key):
        self.current_case = case_key
        import numpy
        self.axes.clear()
        try:
            contourplots = self.contours[case_key]
        except KeyError:
            self.cache_contour_plots(case_key)
            contourplots = self.contours[case_key]
        for coll in contourplots[0].collections:
            self.axes.add_artist(coll)
        self.axes.add_artist(contourplots[1])
        from math import log
        from operator import itemgetter
        contours = itemgetter(0,2)(self.validator.cases[case_key]['cutoffs'])
        P = self.P_limits = [0, -log(contours[0])]
        self.scatter = self.axes.scatter((200),(200), cmap='bwr', picker = 2.0)
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.on_resize()
        
    def update_scatter(self, *_):
        import numpy
        key = self.current_case
        indices = self.validator.case_arrays[key]
        if len(indices):
            phipsi = numpy.degrees(self.validator.phipsi[indices].astype(float))
            logscores = numpy.log(self.validator.rama_scores[indices])
        else:
            phipsi = self.default_coords
            logscores = self.default_logscores
        if phipsi is not None and len(phipsi):
            self.canvas.restore_region(self.background)
            self.scatter.set_offsets(phipsi)
            self.scatter.set_clim(self.P_limits)
            self.scatter.set_array(-logscores)
        else:
            #Just in case of the unexpected
            self.scatter.set_offsets(self.default_coords)
        self.axes.draw_artist(self.scatter)
        self.canvas.blit(self.axes.bbox)

        
class OmegaValidator():
    def __init__(self, annotation_model):
        from chimerax.core.models import Drawing
        self.current_model = None
        self.omega = None
        self.m = annotation_model
        self.name = 'omega planes'
        self.cis_color = [255, 32, 32, 255]
        self.twisted_color = [255,255,32,255]
        existing_names = [d.name for d in self.m.all_drawings()]
        if self.name not in existing_names:          
            self.master_drawing = Drawing(self.name)
            self.m.add_drawing(self.master_drawing)
        else:
            i = existing_names.index(self.name)
            self.master_drawing = self.m.all_drawings()[i]
            self.master_drawing.remove_all_drawings()
        
        self.drawings = {}
        self.currently_drawn = {}
        
        from math import radians
        # maximum deviation from 0 radians to be annotated as cis
        self.cis_max = radians(30)
        # maximum absolute value of torsion angles annotated as twisted
        self.twisted_max = radians(150)
    
    def load_structure(self, model, omega_list):
        self.current_model = model
        self.clear()
        self.omega = [o for o in omega_list if o != None]
        self.drawings = {}
        from chimerax.core.models import Drawing
        for o in self.omega:
            d = self.drawings[o] = Drawing('omega plane')
            self.master_drawing.add_drawing(d)
        self.currently_drawn = {}
    
    def find_outliers(self):
        cis = []
        twisted = []
        for o in self.omega:
            if abs(o.value)  <= self.cis_max:
                cis.append(o)
            elif abs(o.value) <= self.twisted_max:
                twisted.append(o)
        return cis, twisted
    
    def draw_outliers(self, cis, twisted):
        from . import geometry
        for o, d in self.drawings.items():
            d.set_display(False)
        self.currently_drawn = {}
        
        if len(cis):
            for c in cis:
                d = self.drawings[c]
                d.vertices, d.normals, d.triangles = geometry.dihedral_fill_plane(*c.atoms.coords)
                d.set_color(self.cis_color)
                d.set_display(True)
                self.currently_drawn[c] = d
        
        if len(twisted):
            for t in twisted:
                d = self.drawings[t]
                d.vertices, d.normals, d.triangles = geometry.dihedral_fill_plane(*t.atoms.coords)
                d.set_color(self.twisted_color)
                d.set_display(True)
                self.currently_drawn[t] = d
    
    def update_coords(self):
        from . import geometry
        for o, d in self.currently_drawn.items():
            d.vertices, d.normals, d.triangles = geometry.dihedral_fill_plane(*o.atoms.coords)
        
    def clear(self):
        self.master_drawing.remove_all_drawings()
        self.drawings = {}
        self.currently_drawn = {}
    
                
    
    
            
        
    
    
    
        
