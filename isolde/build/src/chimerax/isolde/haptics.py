
class HapticTugger():
    
    def __init__(self, session, index, annotations):
        self.session = session
        self._hh = session.HapticHandler
        self.index = index
        
        self.tugging = False
        
        self.tug_atom = None
        
        self._arrow_model = None
        
        # Model object to draw the arrow into
        self._annotations = annotations
        
    def start_tugging(self, atom):
        self.tug_atom = atom
        self.tugging = True
        self._hh.setTargetPosition(self.index, *atom.coord, scene_coords = True)
        self._hh.startTugging(self.index)
    
    # Set the target position of the device to the current atom position, and
    # return the device position as the target for the atom
    def update_target(self):
        atom_xyz = self.tug_atom.coord
        pointer_xyz = self._hh.getPosition(self.index, scene_coords = True)
        self._hh.setTargetPosition(self.index, *atom_xyz, scene_coords = True)
        self._draw_arrow(atom_xyz, pointer_xyz, radius = 0.1)
        return pointer_xyz
    
    def stop_tugging(self):
        self._hh.stopTugging(self.index)
        self.tugging = False
        self.tug_atom = None
        a = self._arrow_model
        if a:
            a.display = False

    def cleanup(self):
        if self.tugging:
            self.stop_tugging()
        self._delete_arrow()
        
    def _draw_arrow(self, xyz1, xyz2, radius = 0.1):
        from . import geometry
        a = self._arrow_model
        if a is None:
            self._arrow_model = a = geometry.simple_arrow(radius = radius,
                color = [0,255,0,255])
            self._annotations.add_drawing(a)
        # Scale and rotate prototype cone.
        geometry.arrow_between_points(a,xyz1,xyz2)
        a.display = True

    def _delete_arrow(self):
        a = self._arrow_model
        if a is not None:
            self._annotations.remove_drawing(a)
        self._arrow_model = None
    
    
        
