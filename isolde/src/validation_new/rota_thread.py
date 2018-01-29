'''
Thread to handle rotamer validation
'''
import multiprocessing as mp
import numpy
import ctypes
import traceback

from time import sleep
from math import pi


from ..threading import TypedMPArray, SharedNumpyArray, ThreadComms
from ..threading import ChangeTracker as ChangeTracker_Base
from .. import geometry

from ..constants import defaults

PAUSE_SLEEP = defaults.THREAD_PAUSE_SLEEP_INTERVAL 
TWO_PI = 2*pi

class ChangeTracker(ChangeTracker_Base):
    def __init__(self):
        #Inputs
        self.COORDS_READY = 1
        
        self.STOP = 1<<15

        #Outputs
        self.WAITING          = 1<<30
        self.VALIDATION_READY = 1<<31

        super().__init__()

class RotaComms(ThreadComms):
    def __init__(self):
        super().__init__()
        self.update({
            'changes':          ChangeTracker(),
            'error':            mp.Queue(),
            'status':           mp.Queue(),
        })

def _init_thread(validator, params, data, comms):
    global _thread_obj
    _thread_obj = RotaThread(validator, params, data, comms)

def _rota_thread():
    ''' Main loop for live rotamer validation. '''
    global _thread_obj
    
    rt = _thread_obj
    comms = rt.comms
    ct = rt.change_tracker
    changes = ct.changes
    status_q = comms['status']
    error_q = comms['error']
    
    while True:
        with changes.get_lock():
            current_changes = changes.value
        ct.clear_inputs()
        if current_changes & ct.STOP:
            status_q.put('Rotamer validation thread terminated on command.')
            break
        
        try:
            rt.main_loop(current_changes)
        except Exception as e:
            tb = traceback.format_exc()
            main_str = e.args[0]
            raise Exception(main_str + tb)
            error_q.put((e,tb))
            sleep(0.1)
            ct.error()
            break




class RotaThread:
    '''
    Designed to be run in a worker thread under the control of the
    multiprocessing module.
    '''
    def __init__(self, validator, params, data, comms):
        '''
        Args:
            validator:
                A {resname: scipy.RegularGridInterpolator} dict mapping
                3-letter residue names to the relevant probability 
                contours
            params:
                a {name: val} dict of any necessary parameters
            data:
                a {name: val} dict of all run-specific data that *doesn't*
                change during the simulation
            comms:
                a RotaComms object containing all the data that *can*
                change during the simulation.
            change_tracker:
                a ChangeTracker object handling back-and-forth communication
                between this thread and the master
        '''
        self.validator = validator
        self.params = params
        self.data = data
        self.comms = comms
        ct = self.change_tracker = comms['changes']
        self.changes = ct.changes
        self.error_queue = comms['error']
        
        coords = comms['coords']
        self.local_coords = numpy.empty(coords.shape, coords.dtype)
        
        self.num_dihedrals = len(coords)//4
        
        scores = comms['scores']
        self.local_scores = numpy.empty(scores.shape, scores.dtype)
                
        score_map = data['residue ranges']
        self.color_map = data['color map']
        
        self.resnames = list(score_map.keys())
    
    def main_loop(self, changes):
        par = self.params
        comms = self.comms
        ct = self.change_tracker
        data = self.data
        validator = self.validator
        
        if not changes&ct.COORDS_READY:
            #Nothing to do. Check back in a little while
            sleep(PAUSE_SLEEP)  
            ct.register_change(ct.WAITING)
            return
        
        in_coords = comms['coords']
        coords = self.local_coords
        with in_coords.get_lock():
            coords[:] = in_coords
        
        dihedrals = geometry.get_dihedrals(coords, self.num_dihedrals)
        # Dihedral values are (-pi..pi), MolProbity tables are (0..2*pi)
        #~ dihedrals[dihedrals<0] += TWO_PI
        
        scores = self.local_scores
        res_counts = data['residue count']
        dihedral_ranges = data['dihedral ranges']
        residue_ranges = data['residue ranges']
        chi_counts = data['num chi']
        symms = data['symmetric terminal chi']
        
        for resname in self.resnames:
            res_count = res_counts[resname]
            if res_count == 0:
                continue
            dihedral_range = dihedral_ranges[resname]
            score_range = residue_ranges[resname]
            n_chi = chi_counts[resname]
            symm = symms[resname]
            chi_vals = dihedrals[dihedral_range[0]:dihedral_range[1]].reshape((res_count, n_chi))
            if symm:
                chi_vals[:,-1][chi_vals[:,-1]>pi] -= pi
            try:
                scores[score_range[0]:score_range[1]] = validator[resname](chi_vals)
            except:
                raise Exception(resname)
        
        outliers = scores < par['outlier_cutoff']
        allowed = numpy.logical_xor(scores < par['allowed_cutoff'], outliers)
        colors = self.color_map.get_colors(numpy.log(scores))
        
        
        scores_out = comms['scores']
        allowed_out = comms['allowed mask']
        outliers_out = comms['outlier mask']
        colors_out = comms['colors']
        
        with scores_out.get_lock(), allowed_out.get_lock(), \
            outliers_out.get_lock(), colors_out.get_lock():
            scores_out[:] = scores
            allowed_out[:] = allowed
            outliers_out[:] = outliers
            colors_out[:] = colors
        
        ct.register_change(ct.VALIDATION_READY)
        ct.register_change(ct.WAITING)
        
        
        
        
    
