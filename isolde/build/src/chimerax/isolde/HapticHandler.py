#def haptics_start(session):
    #singleton(session).start_haptics()
#from chimerax.core.commands import CmdDesc
#haptics_start_desc = CmdDesc()

#def haptics_stop(session):
    #singleton(session).stop_haptics()
#from chimerax.core.commands import CmdDesc
#haptics_stop_desc = CmdDesc()



#_haptic_handler = None
#def singleton(session):
    #global _haptic_handler
    #if _haptic_handler is None:
        #_haptic_handler = HapticHandler(session)
    #return _haptic_handler

class HapticHandler():
    import ctypes
    _HapticHandler = ctypes.CDLL('./HapticHandler.so')
    
    def __init__ (self, session):
        self.session = session
                
        # Number of connected devices detected by the Chai3D library
        self.numHapticDevices = 0
        
        # Current (x,y,z) position of each device
        self.hapticDevicePosition = []
        
        # Target (x,y,z) position of each device
        self.targetPosition = []
        
        # Spring constant for each device
        self.springConstant = []
        
        # Turn damping on/off for each device
        self.useDamping = []
        
        # Turn force feedback on/off for each device
        self.useFeedback = []
        
        # Is the device attached to something?
        self.deviceInUse = []
        
    
    self._startHaptics = self._HapticHandler.startHaptics
    
    # Initialise all connected haptic devices
    def startHaptics(self):
        self._startHaptics()
    
    self._stopHaptics = self._HapticHandler.stopHaptics
    
    # Stop the haptic device handling process
    def stopHaptics(self):
        self._stopHaptics()
        
    self._getNumDevices = self._HapticHandler.getNumDevices
    self._getNumDevices.restype = ctypes.c_int
    
    # Get the number of detected haptic devices
    def getNumDevices():
        n = self.numHapticDevices = self._getNumDevices()
        return n
    
    self._setSpringConstant = self._HapticHandler.setSpringConstant
    self._setSpringConstant.argtypes = [ctypes.c_int, ctypes.c_double]
    
    # Set the spring constant for device i
    def setSpringConstant(i, k):
        self.springConstant[i] = k
        self._setSpringConstant(i, k)
    
    self._setDamping = self._HapticHandler.setDamping
    self._setDamping.argtypes = [ctypes.c_int, ctypes.c_bool]
    
    # Turn damping on/off for device i
    def setDamping(i, d):
        self.useDamping[i] = d
        self._setDamping(i, d)
    
    self._setTargetPosition = self._HapticHandler.setTargetPosition
    self._setTargetPosition.argtypes = [ctypes.c_int, 
            ctypes.c_double, ctypes.c_double, ctypes.c_double]
    
    # Set the target (x, y, z) position for device i    
    def setTargetPosition(i, x, y, z):
        self.targetPosition[i] = [x, y, z]
        self._setTargetPosition(i, x, y, z)
    
    self._startTugging = self._HapticHandler.startTugging
    self._startTugging.argtypes = [ctypes.c_int]
    
    # Start tugging an object with device i
    def startTugging(i):
        self.deviceInUse[i] = True
        self._startTugging(i)
    
    self._stopTugging = self._HapticHandler.stopTugging
    self._stopTugging.argtypes = [ctypes.c_int]
    
    # Stop tugging an object with device i
    def stopTugging(i):
        self.deviceInUse[i] = False
        self._stopTugging(i)
        
    self._turnOnFeedback = self._HapticHandler.turnOnFeedback
    self._turnOnFeedback.argtypes = [ctypes.c_int]
    
    # Turn on force feedback for device i
    def turnOnFeedback(i):
        self.useFeedback[i] = True
        self._turnOnFeedback(i)
    
    self._turnOffFeedback = self._HapticHandler.turnOffFeedback
    self._turnOffFeedback.argtypes = [ctypes.c_int]
    
    # Turn off force feedback for device i
    def turnOffFeedback(i):
        self.useFeedback[i] = False
        self._turnOffFeedback(i)
        
    self._getPosition = self._HapticHandler.getPosition
    self._getPosition.argtypes = [ctypes.c_int]
    self._getPosition.restype = ctypes.c_double * 3
    
    # Get the current position of device i
    def getPosition(i):
        return self._getPosition(i)
    
    self._getButtonStates = self._HapticHandler.getButtonStates
    self._getButtonStates.argtypes = [ctypes.c_int]
    self._getButtonStates.restype = ctypes.c_bool * 4
    
    # Get the states of the buttons on device i
    def getButtonStates(i):
        return self._getButtonStates(i)
    
        
    
        
    
