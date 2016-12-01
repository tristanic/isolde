
def cone_geometry(radius = 1, height = 1, nc = 10, caps = True, flipped = False):
    '''
    Return vertex, normal vector and triangle arrays for cone geometry
    with specified radius and height. If flipped is true, the base of the
    cone will be at the origin, otherwise the point will be at the origin.
    '''
    from numpy import ones, empty, float32, arange, cos, sin, int32, pi
    vc = nc * 2
    tc = nc
    if caps:
        vc += (nc + 1)
        tc += nc
    varray = empty((vc, 3), float32)
    narray = empty((vc, 3), float32)
    tarray = empty((tc, 3), int32)

    # Compute a circle (which may be used twice if caps is true)
    angles = (2 * pi / nc) * arange(nc)
    import sys
    circle = empty((nc, 2), float32)
    circle[:,0] = cos(angles) * radius
    circle[:,1] = sin(angles) * radius

    # Create cone faces (first nc*2 vertices)
    # The normals are wrong, but let's see how they look
    nc2 = nc * 2
    if not flipped:
        varray[:nc] = (0, 0, 0)      # point of cone (multiple normals)
    else:
        varray[:nc] = (0, 0, height)
    narray[:nc,:2] = circle
    if not flipped:
        narray[:nc,2] = 0
    else:
        narray[:nc,2] = height            
    varray[nc:nc2,:2] = circle      # base of cone
    if not flipped:
        varray[nc:nc2,2] = height
    else:
        varray[nc:nc2,2] = 0
    narray[nc:nc2,:2] = circle      # wrong, but close (enough?)
    if not flipped:
        narray[nc:nc2,2] = height
    else:
        narray[nc:nc2,2] = 0
    tarray[:nc,0] = arange(nc)
    tarray[:nc,1] = (arange(nc) + 1) % nc + nc
    tarray[:nc,2] = arange(nc) + nc

    # Create cone base (last nc+1 vertices)
    if caps:
        if not flipped:
            varray[nc2] = (0, 0, height)
        else:
            varray[nc2] = (0, 0, 0)
        varray[nc2+1:,:2] = circle
        if not flipped:
            varray[nc2+1:,2] = height
        else:
            varray[nc2+1:,2] = 0
        narray[nc2:] = (0, 0, 1)
        tarray[nc:,0] = nc2
        tarray[nc:,1] = (arange(nc) + 1) % nc + nc2 + 1
        tarray[nc:,2] = arange(nc) + nc2 + 1

    return varray, narray, tarray



import ctypes
import numpy
class HapticHandler():
    _HapticHandler = ctypes.CDLL('/home/tic20/chimerax_start/haptic_interface/_HapticHandler.so')
    
    def __init__ (self, session):
        self.session = session
        self.log = session.logger.info
        self._MAX_DEVICES = 16
        
        self._running = False
        
        # Handler for arrow models
        self._arrow_model = [None] * self._MAX_DEVICES
        
        # Event handler for display update
        self._event_handler = None
        
        self._arrow_radius = 0.10
        self._arrow_aspect_ratio = 3
        
        from chimerax.core.geometry import place
        # Arrow rotation matrix
        self._arrow_rotation = place.rotation((0,1,0), -90)
        
        # Mapping of device axes to scene axes and vice versa
        self._axis_order = numpy.array([0, 2, 1])
        
        # Direction of device axes relative to scene axes and vice versa
        self._axis_directions = numpy.array([1, 1, -1])
        
        # Current scaling factor between device and scene coordinates
        self._final_axis_scale = None
        
        
        self._final_axes_and_origin = None
        
        # Axes in the reference frame of the camera
        self._display_axes_and_origin = None
        
        # Scaling factor relating interface coordinates to display coordinates
        self.display_scaling_factor = 110.0
        
        # Scaling factor after taking into account current pixel size.
        # Used to maintain constant arrow size, and to scale spring constants.
        self._final_display_scale = None
                
        # Number of connected devices detected by the Chai3D library
        self.numHapticDevices = 0
        
        # Current (x,y,z) position of each device
        self.hapticDevicePosition = [None] * self._MAX_DEVICES
        
        # Target (x,y,z) position of each device
        self.targetPosition = [None] * self._MAX_DEVICES
        
        # Spring constant for each device, with length units in Angstroms
        self.springConstant = [50] * self._MAX_DEVICES
        
        # Turn damping on/off for each device
        self.useDamping = [True] * self._MAX_DEVICES
        
        # Turn force feedback on/off for each device
        self.useFeedback = [True] * self._MAX_DEVICES
        
        # Is the device attached to something?
        self.deviceInUse = [False] * self._MAX_DEVICES
        
        ####
        # functions defined in _HapticHandler.so
        ####
    
        self._startHaptics = self._HapticHandler.startHaptics
        
        self._stopHaptics = self._HapticHandler.stopHaptics

        self._getNumDevices = self._HapticHandler.getNumDevices
        self._getNumDevices.restype = ctypes.c_int
        
        self._setSpringConstant = self._HapticHandler.setSpringConstant
        self._setSpringConstant.argtypes = [ctypes.c_int, ctypes.c_double]
            
        self._setDamping = self._HapticHandler.setDamping
        self._setDamping.argtypes = [ctypes.c_int, ctypes.c_bool]

        self._setTargetPosition = self._HapticHandler.setTargetPosition
        self._setTargetPosition.argtypes = [ctypes.c_int, 
                ctypes.c_double, ctypes.c_double, ctypes.c_double]

        self._startTugging = self._HapticHandler.startTugging
        self._startTugging.argtypes = [ctypes.c_int]

        self._stopTugging = self._HapticHandler.stopTugging
        self._stopTugging.argtypes = [ctypes.c_int]

        self._turnOnFeedback = self._HapticHandler.turnOnFeedback
        self._turnOnFeedback.argtypes = [ctypes.c_int]
        
        self._turnOffFeedback = self._HapticHandler.turnOffFeedback
        self._turnOffFeedback.argtypes = [ctypes.c_int]
         
        self._getPosition = self._HapticHandler.getPosition
        self._getPosition.argtypes = [ctypes.c_int]
        self._getPosition.restype = ctypes.POINTER(ctypes.c_double * 3)
             
        self._getButtonStates = self._HapticHandler.getButtonStates
        self._getButtonStates.argtypes = [ctypes.c_int]
        self._getButtonStates.restype = ctypes.POINTER(ctypes.c_bool * 4)
         
         
    # Initialise all connected haptic devices
    def startHaptics(self):
        if not self._running:
            self._startHaptics()
            self._event_handler = self.session.triggers.add_handler('new frame', self.on_refresh)
        self._running = True
        self.getNumDevices()
    
    
    # Stop the haptic device handling process
    def stopHaptics(self):
        if self._running:
            self._stopHaptics()
            self.session.triggers.remove_handler(self._event_handler)
            for i in range(self.numHapticDevices):
                a = self._arrow_model[i]
                if a is not None:
                    self.session.models.close([a])
                
        self._running = False
        
    
    # Get the number of detected haptic devices
    def getNumDevices(self):
        if self._running:
            n = self.numHapticDevices = self._getNumDevices()
            return n
    
    
    # Set the spring constant for device i. Length units in Angstroms.
    # Needs to be scaled to the device length units, otherwise feedback
    # strength will change as we zoom.
    def setSpringConstant(self, i, k):
        if self._running:
            self.springConstant[i] = k
            self._setSpringConstant(i, k*self._final_display_scale)
    
    
    # Turn damping on/off for device i
    def setDamping(self, i, d):
        if self._running:
            self.useDamping[i] = d
            self._setDamping(i, d)
    
    
    # Set the target (x, y, z) position for device i, in device or
    # scene coordinates  
    def setTargetPosition(self, i, x, y, z, scene_coords = False):
        if self._running:
            if scene_coords:
                from chimerax.core.geometry import place
                pos_in_scene = numpy.array([x,y,z])
                p1 = place.product([self._final_axes_and_origin.inverse(),pos_in_scene])
                p = (p1/self._final_axis_scale * self._axis_directions)[self._axis_order]
            else:
                p = numpy.array([x,y,z])
            self.targetPosition[i] = p
            self._setTargetPosition(i, *p)
    
    
    # Start tugging an object with device i
    def startTugging(self, i):
        if self._running:
            self.deviceInUse[i] = True
            self._startTugging(i)
    
    
    # Stop tugging an object with device i
    def stopTugging(self, i):
        if self._running:
            self.deviceInUse[i] = False
            self._stopTugging(i)
        
    
    # Turn on force feedback for device i
    def turnOnFeedback(self, i):
        if self._running:
            self.useFeedback[i] = True
            self._turnOnFeedback(i)
    
   
    # Turn off force feedback for device i
    def turnOffFeedback(self, i):
        if self._running:
            self.useFeedback[i] = False
            self._turnOffFeedback(i)
        
    
    # Get the current position of device i
    def getPosition(self, i, scene_coords = False):
        if self._running:
            if scene_coords:
                return self._arrow_model[i].position.origin()
            return [i for i in self._getPosition(i).contents]

    
    
    
    # Get the states of the buttons on device i
    def getButtonStates(self, i):
        if self._running:
            buttons = [b for b in self._getButtonStates(i).contents]
            return buttons
    
    def on_refresh(self, *_):
        display_axes_and_origin, axis_scale = self.get_haptic_reference_frame()
        for i in range(self.numHapticDevices):
            pos = self.getPosition(i)
            self._draw_arrow(i, pos, display_axes_and_origin, axis_scale)
            if self.deviceInUse[i]:
                self._setSpringConstant(i, self.springConstant[i]*self._final_display_scale)

    def get_haptic_reference_frame(self):
        s = self.session
        v = s.main_view
        c = v.camera
        camera_place = c.get_position()
        camera_origin = camera_place.origin()
        camera_view_direction = c.view_direction()
        camera_axes = camera_place.axes()
        window_scale = max(v.window_size) / 1000
        near_far = v.near_far_distances(c, 0)
        origin = (near_far[0] + near_far[1])/2 * camera_view_direction + camera_origin
        from chimerax.core.geometry import place
        p = place.Place(axes = camera_axes, origin = origin)
        final_display_scale = self._final_display_scale = self.display_scaling_factor * v.pixel_size()
        final_axis_scale = self._final_axis_scale = self.display_scaling_factor * window_scale
        arrow_transform = self._final_axes_and_origin = place.product([p, self._arrow_rotation, place.scale(final_display_scale)])
        return arrow_transform, final_axis_scale
        
        
        
        
    def _draw_arrow(self, i, pos, reference_frame, axis_scale = 1.0):
        import numpy
        pos = numpy.array(pos)[self._axis_order] * self._axis_directions
        a = self._arrow_model[i]
        s = self.session
        if a is None or a.deleted:
            from chimerax.core.models import Model
            self._arrow_model[i] = a = Model('Haptic tool ' + str(i), s)
            #from chimerax.core.surface import cone_geometry
            a.vertices, a.normals, a.triangles  = cone_geometry(
                    radius = self._arrow_radius,
                    height = self._arrow_radius * self._arrow_aspect_ratio,
                    nc = 20)
            a.color = (0,255,0,255)
            s.models.add([a])
        # Scale and rotate prototype cone.
        from chimerax.core.geometry import place
        r = reference_frame
        r = place.product([r, place.translation(axis_scale*pos)])
        pixel_size = s.main_view.pixel_size(r.origin())
        if not a.position.same(r, angle_tolerance = 1, shift_tolerance = pixel_size):
            a.position = r

        
        
    
