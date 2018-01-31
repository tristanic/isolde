# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from chimerax.core.models import Drawing, Model

from chimerax.core.ui import MouseMode

class TugAtomsMode(MouseMode):
    name = 'tug'
    #icon_file = 'tug.png'

    def __init__(self, session, atoms):
        MouseMode.__init__(self, session)
        self._tugging = False
        self._last_xy = None
        self.name = 'ISOLDE_mouse_tug'
        a = self._arrow_model = TugArrow(self.name)
        session.main_view.drawing.add_drawing(a)
        self._picked_atom = None
        self._pull_vector = None
        self._xyz0 = None
        self._xy = None
        # Atomic array to pick from
        self._atoms = atoms
        
        # Variables to be set by the caller
        self.last_tugged_atom = None
        self.already_tugging = False

    
    def __del__(self):
        self.cleanup()
        super().__del__()
    
    def cleanup(self):
        self._delete_arrow()

    @property
    def status(self):
        return self._tugging, self._picked_atom, self._xyz0
        print('Atom coords: ' + str(self._picked_atom.scene_coord) + ' xyz0: ' + str(self._xyz0))
    
    @property
    def tugging(self):
        return self._tugging
    
    @tugging.setter
    def tugging(self, flag):
        if flag:
            self._draw_arrow(self._xyz0, atom_xyz)
            self._arrow_model.display = True
        else:
            self._arrow_model.display = False
        self._tugging = flag

    def mouse_down(self, event):
        MouseMode.mouse_down(self, event)
        x,y = event.position()
        self._xy = (x,y)
        view = self.session.main_view
        from . import picking
        pick = picking.pick_closest_to_line(self.session, x, y, self._atoms, 0.5)
        if pick is not None:
            a = self._picked_atom = pick
            atom_xyz, self._pull_vector = self._pull_direction(x, y)
            self._xyz0 = atom_xyz + self._pull_vector
            self.tugging = True


    def mouse_drag(self, event):
        if not self._tugging:
            return
        self._xy = x,y = event.position()
        atom_xyz, self._pull_vector = self._pull_direction(x, y)
        self._xyz0 = atom_xyz + self._pull_vector
        self._draw_arrow(self._xyz0, atom_xyz)


    def mouse_up(self, event):
        MouseMode.mouse_up(self, event)
        self.tugging = False

    def _pull_direction(self, x, y):
        v = self.session.main_view
        x0,x1 = v.clip_plane_points(x, y)
        axyz = self._picked_atom.coord
        # Project atom onto view ray to get displacement.
        dir = x1 - x0
        da = axyz - x0
        from chimerax.core.geometry import inner_product
        offset = da - (inner_product(da, dir)/inner_product(dir,dir)) * dir
        return axyz, -offset

    def _draw_arrow(self, xyz1, xyz2, radius = 0.1):
        from . import geometry
        a = self._arrow_model
        if a is None:
            s = self.session
            self._arrow_model = a = TugArrow(self.name)
        # Scale and rotate prototype cylinder.
        geometry.arrow_between_points(a, xyz2, xyz1)
        a.display = True

    def _delete_arrow(self):
        a = self._arrow_model
        if a is not None:
            self.session.main_view.drawing.remove_drawing(a)
        self._arrow_model = None




class HapticTugger():
    
    def __init__(self, session, index):
        self.session = session
        self._hh = session.HapticHandler
        self.index = index
        
        self.tugging = False
        
        self.tug_atom = None
        
        a = self._arrow = TugArrow(session, "Haptic Tug Arrow {}".format(self.index))
        v = session.main_view
        v.drawing.add_drawing(a)
                
    def start_tugging(self, atom):
        self.tug_atom = atom
        self.tugging = True
        self._arrow.display = True
        self.update_target()
        #~ self._hh.setTargetPosition(self.index, *atom.coord, scene_coords = True)
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
        a = self._arrow
        if a:
            a.display = False
    
    def __del__(self):
        self.cleanup()
    
    def cleanup(self):
        if self.tugging:
            self.stop_tugging()
        self._delete_arrow()
        
    def _draw_arrow(self, xyz1, xyz2, radius = 0.1):
        from . import geometry
        a = self._arrow
        # Scale and rotate prototype cone.
        geometry.arrow_between_points(a,xyz1,xyz2)

    def _delete_arrow(self):
        a = self._arrow
        if a is not None:
            self.session.main_view.drawing.remove_drawing(a)
        self._arrow = None
    
    
class TugArrow(Drawing):
    def __init__(self, session, name, radius = 0.1):
        self.session = session
        super().__init__(name)
        self._create_geometry(0.1)
        self.color=[0,255,0,255]
        self.no_cofr = True # Don't include in cofr calculation
        self.pickable = False
        self.display = False 
    
    def _create_geometry(self, radius):
        from .geometry import simple_arrow_geometry
        self.vertices, self.normals, self.triangles = simple_arrow_geometry(
            radius=radius)
        
            
