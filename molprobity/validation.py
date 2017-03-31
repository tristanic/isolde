import os
from math import degrees, radians
package_directory = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(package_directory, 'molprobity_data')

def generate_interpolator(file_prefix, wrap_axes = True):
    '''
    The very first time this is run after installation, it will read in
    a MolProbity data set and return a SciPy RegularGridInterpolator 
    object for fast interpolation of scores for given values. It will 
    also write a pickle file of the RegularGridInterpolator object so 
    that subsequent loading will be much faster. Why do it this way? 
    Mainly for maximum compatibility and human readability - we want the
    text files to be there anyway so others can read and understand them,
    and this way ensures that changes to the pickle format can't 
    break functionality.
    This routine should be able to read MolProbity datasets of any
    dimensionality, but a little more work will be needed for any 
    dataset where only some axes are periodic (if any such set exists).
    '''
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

CIS_MIN = radians(-30)
CIS_MAX = radians(30)
TRANS_MIN = radians(150)
TRANS_MAX = radians(-150)

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
        '''
        The Ramachandran validator handles the task of looking up 
        protein phi/psi values against the probability maps corresponding
        to each residue pair. Once initialised, simply call get_scores()
        with a properly-initialised Backbone_Dihedrals object and the
        optional boolean update_colors as arguments. The latter is used
        for live updating of C-alpha colours for tracking of Ramachandran
        scores in a moving simulation.
        Args:
            color_scale:
                One of the following:
                    'RWB': Red-white-blue
                    'BWR': Blue-white-red
                    'RGB': Red-green-blue
                    'BGR': Blue-green-red
                    'RWG': Red-white-green
                    'GWR': Green-white-red
                    'GYO': Green-yellow-orange
                    'OYG': Orange-yellow-green
                    'PiYG': Pink-yellow-green
                    'GYPi': Green-yellow-pink
                Colours change according to the log of the Ramachandran
                score, so the colours are reversed (last colour 
                corresponds to the highest probability score and vice
                versa). The mid-point of the scale is defined as the
                cut-off between "favoured" and "allowed" - so, for 
                example, the default PiYG scale for the General case 
                will change smoothly from green to yellow for scores
                between 1 and 0.02, and then from yellow to pink for
                scores between 0.02 and 0.0005, with no further change
                for lower scores.
        '''
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
        '''
        Sort a pre-calculated set of scores into their classic "favoured,
        allowed, outlier" bins. Might be better to roll this in as an 
        optional argument to get_scores().
        Args:
            scores:
                An iterable of Ramachandran scores
            types:
                An iterable of case names chosen from those in RAMA_CASES
        '''
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
    '''
    A Ramachandran plot implementation based on MatPlotLib, optimised for
    fast redrawing (e.g. to run live during simulation or playback of 
    trajectories) and with functionality to select and display a residue
    when its corresponding point is clicked. Plotted points are coloured
    by their score in a similar fashion to the colouring of atoms by
    RamaValidator (albeit using the MatPlotLib built-in colour scales).
    '''
    def __init__(self, session, container, validator):
        '''
        Initialise the Ramachandran plot widget.
        Args:
            session:
                The ChimeraX session
            container:
                A QtWidgets object suitable for holding and displaying
                the plot (e.g. a QtWidgets.QWidget).
            validator:
                the RamaValidator object that will be responsible for
                recalculating scores.
        '''
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
        '''
        The pcolor plot (which generates a smoothly-varying background 
        coloured according to the local probability value) is very 
        pretty and meaningful, but is somewhat expensive to calculate. 
        Re-calculating it every time the user switches between cases 
        causes noticeable delays, so we'll cache each one the first time
        it's loaded.
        '''
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
        '''
        Resize and reshape the plot when its host window is resized.
        Could use a little more work to make things prettier.
        '''
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
        '''
        Selects the picked residue, and updates the molecule view to
        zoom and focus on it.
        '''
        ind = event.ind[0]
        res_index = self._last_bd.rama_cases[self.current_case][ind]
        picked_residue = self._last_bd.residues[res_index]
        from . import view
        view.focus_on_selection(self.session, self.session.main_view, picked_residue.atoms)
        
    def change_case(self, case_key):
        '''
        Update the plot to show the chosen Ramachandran case.
        '''
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
        self.scatter = self.axes.scatter((200),(200), picker = 2.0)
        self.scatter.set_cmap('bwr')
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.on_resize()
        
    def update_scatter(self, bd = None, force_update = False):
        '''
        Update the scatter plot with dihedral values from the given 
        Backbone_Dihedrals object, or resets it to an empty plot.
        Args:
            bd:
                A Backbone_Dihedrals object defining the residues you 
                want to plot (or None to clear the plot).
            force_update:
                If True the phi, psi and omega values will be recalculated
                from the atomic coordinates. Otherwise the last values
                stored in bd will be used.
        '''
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
                phipsi = self.default_coords
                logscores = self.default_logscores
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

        
class OmegaValidator:
    '''
    Tracking and highlighting of cis and twisted peptide bonds.
    '''
    def __init__(self, annotation_model):
        '''
        The OmegaValidator is designed for fast identification, tracking
        and highlighting of cis and twisted peptide bonds.
        Args:
            annotation_model:
                A ChimeraX Model object (or subclass) to act as the 
                container for the Drawing object highlighting problematic
                peptides.
        '''
        from chimerax.core.models import Drawing
        self.current_model = None
        self.omega = None
        self.m = annotation_model
        self.name = 'omega planes'
        self.cis_drawing = Drawing('cis planes')
        self.cis_color = [255, 32, 32, 255]
        self.cis_drawing.set_color(self.cis_color)
        self.twisted_drawing = Drawing('twisted planes')
        self.twisted_color = [255,255,32,255]
        self.twisted_drawing.set_color(self.twisted_color)
        
        existing_names = [d.name for d in self.m.all_drawings()]
        if self.name not in existing_names:          
            self.master_drawing = Drawing(self.name)
            self.m.add_drawing(self.master_drawing)
        else:
            i = existing_names.index(self.name)
            self.master_drawing = self.m.all_drawings()[i]
            self.master_drawing.remove_all_drawings()
        self.master_drawing.add_drawing(self.cis_drawing)
        self.master_drawing.add_drawing(self.twisted_drawing)
        
            
    def load_structure(self, model, omega_list):
        self.current_model = model
        self.clear()
        self.omega = [o for o in omega_list if o != None]
    
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
        import numpy
        self.clear()
        
        if len(cis):
            d = self.cis_drawing
            vlist = []
            nlist = []
            tlist = []
            
            for i, c in enumerate(cis):
                tv, tn, tt = geometry.dihedral_fill_plane(*c.atoms.coords)
                vlist.append(tv)
                nlist.append(tn)
                tlist.append(tt + i*len(tv))
            v = numpy.concatenate(vlist)
            n = numpy.concatenate(nlist)
            t = numpy.concatenate(tlist)
            d.vertices, d.normals, d.triangles = (v, n, t)

        
        if len(twisted):
            d = self.twisted_drawing
            vlist = []
            nlist = []
            tlist = []
            
            for i, t in enumerate(twisted):
                tv, tn, tt = geometry.dihedral_fill_plane(*t.atoms.coords)
                vlist.append(tv)
                nlist.append(tn)
                tlist.append(tt + i*len(tv))
            v = numpy.concatenate(vlist)
            n = numpy.concatenate(nlist)
            t = numpy.concatenate(tlist)
            d.vertices, d.normals, d.triangles = (v, n, t)
    
        
    def clear(self):
        for d in (self.cis_drawing, self.twisted_drawing):
            d.vertices, d.normals, d.triangles = (None, None, None)
    
                
    
    
            
        
    
    
    
        
