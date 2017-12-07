'''
All structural validation tasks should be run in (a) thread(s) to reduce their
impact on graphics performance. The interface between ISOLDE and the
validation thread(s) is defined here.
'''
import numpy
import multiprocessing as mp
import ctypes
from math import pi, radians, degrees
from warnings import warn
from time import time
from copy import deepcopy
from chimerax.core.atomic import concatenate, Bonds, Residues

from ..threading.shared_array import TypedMPArray, SharedNumpyArray
from ..constants import defaults, validation_cutoffs
from ..param_mgr import Param_Mgr, autodoc, param_properties
from ..dihedrals import Dihedrals

from . import validation
from . import rama_thread
from . import rota_thread

from .rota_annotation import Rotamer_Annotations

FLOAT_TYPE = defaults.FLOAT_TYPE

def error_cb(e):
    print(e.__traceback__)
    print(e)

def _start_thread(thread_type, validator, params, data, comms):
    '''
    Start a thread for one type of validation.
    Args:
        thread_type:
            The name of the Python class/module defining the thread. Must
            contain a method called _init_thread 
        validator:
            If a pre-initialised object is used for validation, it should
            be declared here. Must be deep-copyable.
        params:
            A Param_Mgr subclass holding all pre-defined parameters.
        data:
            The input data needed to define the task (should not be 
            changed after the thread is started)
        comms:
            Thread-safe container for back-and-forth communication
        change_tracker:
            Manages flags to notify master and slave threads of changes,
            and runs callbacks in the thread as necessary.
    '''
    from multiprocessing.pool import Pool
    thread = Pool(processes=1, initializer=thread_type._init_thread, 
                  initargs=(deepcopy(validator), params, data, comms))
    return thread


######################################
# ROTAMERS
######################################

@param_properties
@autodoc
class RotaParams(Param_Mgr):
    '''
    Default parameters for rotamer validation
    '''
    _default_params = {
        'allowed_cutoff':           (validation_cutoffs.ROTA_ALLOWED_CUTOFF, None),
        'outlier_cutoff':           (validation_cutoffs.ROTA_OUTLIER_CUTOFF, None),
    }

class RotaValidationThreadInterface:
    def __init__(self, session, isolde, rotamers, validator, drawing, params = None):
        '''
        In order to maximise the performance in the main ChimeraX 
        thread, all the rotameric dihedrals, their associated residues, 
        and their CA-CB bonds (for annotation purposes) are flattened 
        into single Dihedrals, Residues and Bonds arrays, respectively. 
        Maps are maintained to keep track of what goes where.
        
        '''
        self.session = session
        self.isolde = isolde
        self.rotamers = rotamers
        self.drawing = drawing
        self._counter = 0
        comms = self.comms_object = rota_thread.RotaComms()
        ct = self.change_tracker = comms['changes']
        if params is None:
            params = self.params = RotaParams()
        
        data = self.data = dict()
        dihedral_map = data['dihedral ranges'] = dict()
        score_map = data['residue ranges'] = dict()
        num_chi = data['num chi'] = dict()
        res_count = data['residue count'] = dict()
        symm = data['symmetric terminal chi'] = dict()
        data['color map'] = deepcopy(drawing.color_map)
                
        start = 0
        r_count = 0
        
        dihedrals = Dihedrals()
        ca_cb_bonds = Bonds()
        residues = Residues()
        
        for resname in validator.keys():
            cur_dihedrals = rotamers.dihedrals(resname)
            n = len(cur_dihedrals)
            n_chi = num_chi[resname] = rotamers.num_chi(resname)
            symm[resname] = rotamers.is_symmetric(resname)
            this_r_count = res_count[resname] = n//n_chi
            dihedral_map[resname] = (start, start+n)
            score_map[resname] = (r_count, r_count+this_r_count)
            ca_cb_bonds = concatenate((ca_cb_bonds, cur_dihedrals[0::n_chi].axial_bonds))
            residues = concatenate((residues, cur_dihedrals[0::n_chi].residues))
            start += n
            dihedrals.append(cur_dihedrals)
            r_count += this_r_count
        
        self._master_dihedral_list = dihedrals
        self._ca_cb_bond_list = ca_cb_bonds
        self._residues = residues
        
        n = len(dihedrals)
        coords = comms['coords'] = SharedNumpyArray(TypedMPArray(
            FLOAT_TYPE, n*4*3)).reshape((n*4,3))
        scores = comms['scores'] = SharedNumpyArray(TypedMPArray(
            FLOAT_TYPE, r_count))
        comms['allowed mask'] = SharedNumpyArray(TypedMPArray(
            ctypes.c_bool, r_count))
        comms['outlier mask'] = SharedNumpyArray(TypedMPArray(
            ctypes.c_bool, r_count))
        comms['colors'] = SharedNumpyArray(TypedMPArray(
            ctypes.c_uint8, r_count*4)).reshape((r_count,4))
        
        self._current_scores = numpy.zeros(r_count, dtype=FLOAT_TYPE)
        self._current_allowed = numpy.zeros(r_count, dtype=numpy.bool)
        self._current_outliers = numpy.zeros(r_count, dtype=numpy.bool)
        self._current_colors = numpy.zeros((r_count,4), dtype=numpy.uint8)
        
        thread = self.thread = _start_thread(rota_thread, validator, 
                    params, data, comms)
        self.thread_result = thread.apply_async(rota_thread._rota_thread, args=(), error_callback=error_cb)
        
        
    def update(self, *_):
        isolde = self.isolde
        comms = self.comms_object
        ct = self.change_tracker
        data = self.data
        dihedrals = self._master_dihedral_list
        bonds = self._ca_cb_bond_list
        changes = ct.changes.value
        ct.clear_outputs()
        if changes & ct.VALIDATION_READY:
            # Only provide new coordinates if it's finished with the last set
            scores = comms['scores']
            am = comms['allowed mask']
            om = comms['outlier mask']
            colors = comms['colors']
            with scores.get_lock(), am.get_lock(), om.get_lock(), colors.get_lock():
                self._current_scores[:]=scores
                self._current_allowed[:]=am
                self._current_outliers[:]=om
                self._current_colors[:]=colors
            self._scores_changed_cb()
        c = self._counter = (self._counter+1)%isolde.params.rounds_per_rota_update
        if c == 0 and changes&ct.WAITING:
            comms.thread_safe_set_array_values('coords', dihedrals.coords)
            ct.register_change(ct.COORDS_READY)
        
    def _scores_changed_cb(self):
        scores = self._current_scores
        a_m = self._current_allowed
        o_m = self._current_outliers
        colors = self._current_colors
        iffy_m = numpy.logical_or(a_m,o_m)
        bonds = self._ca_cb_bond_list
        iffy_b = bonds[iffy_m]
        iffy_scores = scores[iffy_m]
        iffy_colors = colors[iffy_m]
        self.drawing.update_scores(iffy_b, iffy_scores, colors = iffy_colors)
        
    def stop_thread(self):
        ct = self.change_tracker
        ct.register_change(ct.STOP)
        self.thread = None    
        
        
        
########################################
# RAMACHANDRAN
########################################

@param_properties
@autodoc
class RamaParams(Param_Mgr):
    '''
    Default parameters for Ramachandran validation
    '''
    _default_params = {
    }

class RamaValidationThreadInterface:
    def __init__(self, session, isolde, backbone_dihedrals, validator, params = None):
        '''
        Threaded implementation of backbone dihedral validation, 
        handling both Ramachandran and omega (peptide bond cis/trans) 
        geometries and scoring.
        '''
        self.session = session
        self.isolde = isolde
        bd = self.backbone_dihedrals = backbone_dihedrals
        self.validator = validator
        self._counter = 0
        comms = self.comms_object = rama_thread.RamaComms()
        ct = self.change_tracker = comms['changes']
        if params is None:
            params = self.params = RamaParams()
        
            
        irc = bd.rama_cases
        
        rama_cases = data['rama cases'] = dict()
        # Proline residues need to be re-sorted every time in case they
        # flip from cis to trans or vice versa.
        for key, indices in irc.items():
            if key not in ('CisPro', 'TransPro'):
                rama_cases[key] = indices.copy()
        all_pro = rama_cases['AllPro'] = \
            numpy.concatenate((irc['CisPro'], irc['TransPro']))
        cispro_mask = comms['cis-pro mask'] = \
            SharedNumpyArray(TypedMPArray(ctypes.c_bool, len(all_pro)))
        cispro_mask[:] = False
        cispro_mask[0:len(irc['CisPro'])] = True
        
        # While the Backbone_Dihedrals class has methods to quickly get
        # Phi, Psi and Omega dihedral values, we want to push as much of
        # the work as possible into the thread. So, it's better to send
        # through the raw coordinates and calculate the dihedral values
        # there.
        
        phi_coords = bd.phi.coords
        psi_coords = bd.psi.coords
        omega_coords = bd.omega.coords
        
        comms['phi coords'] = SharedNumpyArray(TypedMPArray(
            phi_coords.size, FLOAT_TYPE)).reshape(phi_coords.shape)
        comms['psi coords'] = SharedNumpyArray(TypedMPArray(
            psi_coords.size, FLOAT_TYPE)).reshape(psi_coords.shape)
        comms['omega coords'] = SharedNumpyArray(TypedMPArray(
            omega_coords.size, FLOAT_TYPE)).reshape(omega_coords.shape)
    
        comms['scores'] = SharedNumpyArray(TypedMPArray(
            bd.rama_scores.size, FLOAT_TYPE))
        comms['colors'] = SharedNumpyArray(TypedMPArray(
            bd.rama_colors.size, ctypes.c_uint8)).reshape(bd.rama_colors.shape)
        

    
class ChimeraXValidationInterface:
    '''
    Application-facing interface between ChimeraX and the various 
    validation threads. Creates threads as necessary and keeps track of
    the master validator objects.
    '''    
    
    
    def __init__(self, session, isolde):
        self.session = session
        self.isolde = isolde
        
        ####
        # Rotamers
        ####
        self.rota_validator = validation.RotaValidator()
        self.rota_params = RotaParams()

        ####
        # Ramachandran
        ####
        self.rama_validator = validation.RamaValidator()
        self.rama_params = RamaParams()
    
    @property
    def current_model(self):
        return self.isolde.selected_model
    
    
    def _rota_annotation_model(self):
        m = self.current_model
        if m is not None:
            for cm in m.child_models():
                if isinstance(cm, Rotamer_Annotations):
                    return cm
            rp = self.rota_params
            rm = Rotamer_Annotations(self.session, m, rp.allowed_cutoff, 
                                                     rp.outlier_cutoff)
            m.add([rm])
            return rm
    
    def start_validation_threads(self, rotamers):
        ri = self.rota_interface = RotaValidationThreadInterface(self.session, 
                                self.isolde, rotamers, self.rota_validator,
                                self._rota_annotation_model(), 
                                params = self.rota_params)
                                
        self.isolde._isolde_events.add_event_handler('sim rota update',
                    'completed simulation step', ri.update)
    
    def stop_validation_threads(self):
        self.rota_interface.stop_thread()
        self.isolde._isolde_events.remove_event_handler('sim rota update')
        self.rota_interface = None
        
        
    
    
