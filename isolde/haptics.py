
class HapticTugger():
    
    def __init__(self, session, index):
        self.session = session
        self._hh = session.HapticHandler
        self.index = index
        
        self.tugging = False
        
        self.tug_atom = None
        
        self._arrow_model = None
        
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
        if a and not a.deleted:
            a.display = False

    def cleanup(self):
        self._delete_arrow()
        
    def _draw_arrow(self, xyz1, xyz2, radius = 0.1):
        a = self._arrow_model
        if a is None or a.deleted:
            from chimerax.core.models import Model
            s = self.session
            self._arrow_model = a = Model('Tug arrow', s)
            from chimerax.core.surface import cone_geometry
            a.vertices, a.normals, a.triangles  = cone_geometry()
            a.color = (0,255,0,255)
            s.models.add([a])
        # Scale and rotate prototype cylinder.
        from chimerax.core.atomic import structure
        from numpy import array, float32
        p = structure._bond_cylinder_placements(xyz1.reshape((1,3)),
                                                xyz2.reshape((1,3)),
                                                array([radius],float32))
        a.position = p[0]
        a.display = True

    def _delete_arrow(self):
        a = self._arrow_model
        if a is not None:
            self.session.models.close([a])
        self._arrow_model = None
    
    
        
