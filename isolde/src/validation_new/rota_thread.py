'''
Thread to handle rotamer validation
'''
import multiprocessing as mp
import numpy
import ctypes
import traceback
from ..threading import TypedMPArray, SharedNumpyArray, ThreadComms
from ..threading import ChangeTracker as ChangeTracker_Base

class ChangeTracker(ChangeTracker_Base):
    def __init__(self):
        #Inputs
        self.COORDS_READY = 1

        #Outputs
        self.VALIDATION_READY = 1<<31

        super().__init__()

class RotaComms(ThreadComms):
    super().__init__()
    self.update({
        'changes':          ChangeTracker(),
        'error':            mp.Queue(),
        'status':           mp.Queue(),
    })
