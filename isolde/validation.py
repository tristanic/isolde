import os
from math import degrees, radians
package_directory = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(package_directory, 'molprobity_data')

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

# Master list of Ramachandran case keys. As well as the official ones, we'll
# list N- and C-terminal residues here for convenience.
RAMA_CASES = ['Nter', 'Cter', 'CisPro', 'TransPro', 'Glycine', 'PrePro', 'IleVal', 'General']

RAMA_CASE_DETAILS = {
    'Nter': {
        'name': 'N-terminal residues',
        'file_prefix': None,
        'cutoffs': None
    },
    'Cter': {
        'name': 'C-terminal residues',
        'file_prefix': None,
        'cutoffs': None
    },
    'CisPro': {
        'name': 'Cis-proline residues',
        'file_prefix': os.path.join(DATA_DIR, 'rama8000-cispro'),
        'cutoffs': [0.002, 1.0, 0.02]
    },
    'TransPro': {
        'name': 'Trans-proline residues',
        'file_prefix': os.path.join(DATA_DIR, 'rama8000-transpro'),
        'cutoffs': [0.001, 1.0, 0.02]
    },
    'Glycine': {
        'name': 'Glycine residues',
        'file_prefix': os.path.join(DATA_DIR, 'rama8000-gly-sym'),
        'cutoffs': [0.001, 1.0, 0.02]
    },
    'PrePro': {
        'name': 'Residues preceding proline',
        'file_prefix': os.path.join(DATA_DIR, 'rama8000-prepro-noGP'),
        'cutoffs': [0.001, 1.0, 0.02]
    },
    'IleVal': {
        'name': 'Isoleucine or valine residues',
        'file_prefix': os.path.join(DATA_DIR, 'rama8000-ileval-nopreP'),
        'cutoffs': [0.001, 1.0, 0.02]
    },
    'General': {
        'name': 'General amino acid residues',
        'file_prefix': os.path.join(DATA_DIR, 'rama8000-general-noGPIVpreP'),
        'cutoffs': [0.0005, 1.0, 0.02]
    }
}




def sort_into_rama_cases(counts_for_rama, rama_resnames, omega_vals):
    '''
    Sorts a list of residues into the different Ramachandran cases.
    Arguments:
        counts_for_rama: 1D boolean array determining which residues
                         are valid for Ramachandran analysis (i.e. have
                         both phi and psi dihedrals).
        rama_resnames:   2D array containing the names of the first and 
                         second residue in all psi dihedrals
        omega_values:    Values of the omega dihedrals for sorting into
                         cis and trans proline
    Output:
        a dict containing a numpy array of residue indices corresponding to
        each standard Ramachandran case, plus N- and C-terminal residues since
        we can easily do it here.
    '''
    import numpy
    
    case_arrays = {}
    for case in RAMA_CASES:
        case_arrays[case] = []
    
    rama_case_list = []
     
    for i, [counts, names, oval] in enumerate(zip(counts_for_rama, rama_resnames, omega_vals)):
        name1, name2 = names
        if not counts:
            if name1 is not None:
                case_arrays['Nter'].append(i)
                rama_case_list.append('Nter')
            else:
                case_arrays['Cter'].append(i)
                rama_case_list.append('Cter')
        elif name1 == 'PRO':
            otype = omega_type(oval)
            if otype == 'cis':
                case_arrays['CisPro'].append(i)
                rama_case_list.append('CisPro')
            else:
                case_arrays['TransPro'].append(i)
                rama_case_list.append('TransPro')
        elif name1 == 'GLY':
            case_arrays['Glycine'].append(i)
            rama_case_list.append('Glycine')
        elif name2 == 'PRO':
            case_arrays['PrePro'].append(i)
            rama_case_list.append('PrePro')
        elif name1 in ['ILE', 'VAL']:
            case_arrays['IleVal'].append(i)
            rama_case_list.append('IleVal')
        else:
            case_arrays['General'].append(i)
            rama_case_list.append('General')
    
    # Convert lists to numpy integer arrays
    for key, arr in case_arrays.items():
        case_arrays[key] = numpy.array(arr,numpy.int32)

    return case_arrays, numpy.array(rama_case_list)
    

CIS_MIN = radians(-30)
CIS_MAX = radians(30)
TRANS_MIN = radians(150)
TRANS_MAX = radians(-150)

    
def omega_type(omega):
    if omega is None:
        return None
    if omega >= TRANS_MIN or omega <= TRANS_MAX:
        return "trans"
    elif omega >= CIS_MIN and omega <= CIS_MAX:
        return "cis"
    return "twisted"

class RamaCase:
    '''
    Holds the RegularGridInterpolator object and defining parameters for
    a specific MolProbity Ramachandran case.
    '''
    def __init__(self, file_prefix, cutoffs, color_scale_name):
        self._cutoffs = cutoffs
        self._log_cutoffs = None
        self._file_prefix = file_prefix
        self._interpolator = None
        self.color_scale_name = color_scale_name
        self._color_scale = None
    
    @property
    def interpolator(self):
        if self._interpolator is None:
            self._interpolator = generate_interpolator(self._file_prefix)
        return self._interpolator
    
    @property
    def color_scale(self):
        if self._color_scale is None:
            import numpy
            self._log_cutoffs = numpy.log(numpy.array(self._cutoffs, numpy.float32))
            from . import color            
            self._color_scale = color.standard_three_color_scale(
                self.color_scale_name, *self._log_cutoffs)
        return self._color_scale
    
    @color_scale.setter
    def color_scale(self, scale_name):
        from . import color

        try:
            self._color_scale = color.standard_three_color_scale(
                scale_name, *self._log_cutoffs)
            self.color_scale_name = scale_name
        except:
            import warnings
            errstring = 'Please choose a valid color scale! Available scales are: {}'\
                    .format(color.color_scales.keys())
            warnings.warn(errstring)
    
    @property
    def cutoffs(self):
        return self._cutoffs
    
    @cutoffs.setter
    def cutoffs(self, cutoffs):
        c = cutoffs
        if not hasattr(c, '__len__') or len(c) != 3:
            errstring = 'Please provide a list of three values ordered (smallest, largest, midpoint)'
            raise TypeError(errstring)
        import numpy
        c = numpy.array(c, numpy.float32)
        first_smallest = numpy.all(c[[1,2]] > c[0])
        second_largest = c[1] > c[2]
        if not (first_smallest and second_largest):
            raise TypeError('Cutoffs must be in order [smallest, largest, midpoint]!')
        if numpy.any(c <= 0):
            raise TypeError('Only positive non-zero cutoffs are allowed!')
        self._cutoffs = c
        self._log_cutoffs = numpy.log(c)
        
    


        
class RamaValidator():
    '''
    Encapsulates the MolProbity Ramachandran data with a set of handy
    functions for quick look-up of scores for all residues in a
    Backbone_Dihedrals object.
    '''
    
    import numpy
    from math import pi
    
    def __init__(self, color_scale = 'PiYG'):
        # Cached phi and psi values for plotting
        self.phipsi = None
        # Generate the objects handling the interpolation for each case
        self._rama_cases = {}
        for case in RAMA_CASES:
            cd = RAMA_CASE_DETAILS[case]
            prefix = cd['file_prefix']
            cutoffs = cd['cutoffs']
            if prefix is None:
                # Cases that don't have Ramachandran scores
                self._rama_cases[case] = None
                continue
            self._rama_cases[case] = RamaCase(prefix, cutoffs, color_scale)
            
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

    def get_scores(self, bd, update_colors = False):
        '''
        Update the Ramachandran scores for a Backbone_Dihedrals object bd.
        The scores will be automatically filled into the bd.rama_scores
        array. If update_colors is set to True, the colours of the CA
        atoms in the structure will be changed to reflect their scores.
        We'll also cache the Phi and Psi values here for plotting purposes -
        we don't want to be calculating them twice.
        '''
        import numpy
        residues = bd.residues
        phi = bd.phi_vals
        psi = bd.psi_vals
        phipsi = self.phipsi = numpy.column_stack([phi,psi])
        case_index_dict = bd.rama_cases        
        scores = bd.rama_scores
        CAs = bd.CAs
        if update_colors:
            colors = bd.rama_colors
        
        for case in RAMA_CASES:
            if self._rama_cases[case] is None:
                # We're not interested in these residues
                continue
            else:
                rc = self._rama_cases[case]
                v = rc.interpolator
                indices = case_index_dict[case]
                case_scores = v(phipsi[indices])
                scores[indices] = case_scores
                if update_colors:
                    if len(case_scores):
                        cs = rc.color_scale
                        case_colors = cs.get_colors(numpy.log(case_scores))
                        colors[indices] = case_colors
                        
    @property
    def rama_cases(self):
        return self._rama_cases
    
                        
class RamaPlot():
    
    def __init__(self, session, container, validator):
        import numpy
        self.session = session
        self.container = container
        self.validator = validator
        self.current_case = None
        self._last_bd = None
        
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
        case_data = self.validator.rama_cases[key].interpolator
        grid = numpy.degrees(case_data.grid)
        from operator import itemgetter
        contours = itemgetter(0,2)(RAMA_CASE_DETAILS[key]['cutoffs'])
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
        self.update_scatter(self._last_bd)
    
    def on_pick(self, event):
        ind = event.ind[0]
        res_index = self._last_bd.rama_cases[self.current_case][ind]
        picked_residue = self._last_bd.residues[res_index]
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
        contours = itemgetter(0,2)(RAMA_CASE_DETAILS[case_key]['cutoffs'])
        P = self.P_limits = [0, -log(contours[0])]
        self.scatter = self.axes.scatter((200),(200), cmap='bwr', picker = 2.0)
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.on_resize()
        
    def update_scatter(self, bd = None, force_update = False):
        import numpy
        key = self.current_case
        rv = self.validator
        if bd is not None:
            self._last_bd = bd
            indices = bd.rama_cases[key]
            if len(indices):
                if force_update:
                    rv.get_scores(bd)
                phipsi = numpy.degrees(rv.phipsi[indices].astype(float))
                logscores = numpy.log(bd.rama_scores[indices])
        else:
            self._last_bd = None
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
            if abs(o.value)  <= CIS_MAX:
                cis.append(o)
            elif abs(o.value) <= TRANS_MIN:
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
    
                
    
    
            
        
    
    
    
        
