'''
Thread to handle peptide backbone (Ramachandran and omega) validation.
'''
import multiprocessing as mp
import numpy
import ctypes
import traceback

from ..threading import TypedMPArray, SharedNumpyArray, ThreadComms
from ..threading import ChangeTracker as ChangeTracker_Base
from . import validation

from ..constants import defaults
PAUSE_SLEEP = defaults.THREAD_PAUSE_SLEEP_INTERVAL 


class ChangeTracker(ChangeTracker_Base):
    def __init__(self):
        #Inputs
        self.COORDS_READY = 1
        
        self.STOP = 1<<15

        #Outputs
        self.WAITING          = 1<<30
        self.VALIDATION_READY = 1<<31

        super().__init__()

class RamaComms(ThreadComms):
    def __init__(self):
        super().__init__(self)
        self.update({
            'changes':          ChangeTracker(),
            'error':            mp.Queue(),
            'status':           mp.Queue(),
        })

def _init_rama_thread(params, data, comms_obj, change_tracker):
    global _rama_thread_obj
    _rama_thread_obj = RamaThread(params, data, comms_obj, change_tracker)

def _rama_thread():
    '''
    Main loop for the Ramachandran validation thread.
    '''
    global _rama_thread_obj
    ro = _rama_thread_obj
    comms = ro.comms
    ct = ro.change_tracker
    changes = ct.changes
    status_q = comms['status']
    error_q = comms['error']
    while(True):
        with changes.get_lock():
            current_changes = changes.value
            ct.clear_inputs()

        if current_changes&ct.STOP:
            status_q.put('Ramachandran thread shut down successfully')
            break

        try:
            ro.main_loop(current_changes)
        except Exception as e:
            tb = traceback.format_exc()
            main_str = e.args[0]
            raise Exception(main_str + tb)
            error_q.put((e, tb))
            sleep(0.1)
            ct.error()
            break

class RamaThread:
    '''
    Designed to be run in a thread to provide continuous Ramachandran and peptide
    omega validation.
    '''

    def __init__(self, params, data, comms, change_tracker):
        self.params = params
        self.data = data
        self.comms = comms
        self.change_tracker = change_tracker

        error_q = comms['error']

        self.rama_validator = validation.RamaValidator()
        self.omega_validator = validation.OmegaValidator()
