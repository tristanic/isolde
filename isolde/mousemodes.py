# vim: set expandtab shiftwidth=4 softtabstop=4:

# Useful mouse modes for working with ISOLDE

from chimerax.core.ui import MouseMode

class TugAtomsMode(MouseMode):
    name = 'tug'
    #icon_file = 'tug.png'

    def __init__(self, session, atoms, annotations):
        MouseMode.__init__(self, session)
        self._tugging = False
        self._last_xy = None
        self._arrow_model = None
        self._picked_atom = None
        self._pull_vector = None
        self._xyz0 = None
        self._xy = None
        self.name = 'ISOLDE_mouse_tug'
        self._handler = None
        # Atomic array to pick from
        self._atoms = atoms
        # Model to draw the arrow into
        self._annotations = annotations
    
    def cleanup(self):
        self._delete_arrow()
        
    def get_status(self):
        return self._tugging, self._picked_atom, self._xyz0
        print('Atom coords: ' + str(self._picked_atom.scene_coord) + ' xyz0: ' + str(self._xyz0))
        
            
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
            self._tugging = True
        

    def mouse_drag(self, event):
        if not self._tugging:
            return
        self._xy = x,y = event.position()
        atom_xyz, self._pull_vector = self._pull_direction(x, y)
        self._xyz0 = atom_xyz + self._pull_vector
        self._draw_arrow(self._xyz0, atom_xyz)
        
        
    def mouse_up(self, event):
        MouseMode.mouse_up(self, event)
        self._tugging = False
        a = self._arrow_model
        if a:
            a.display = False
                
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
            self._arrow_model = a = geometry.simple_arrow(radius = radius, 
                    color = [0,255,0,255])
            self._annotations.add_drawing(a)
        # Scale and rotate prototype cylinder.
        geometry.arrow_between_points(a, xyz2, xyz1)
        a.display = True

    def _delete_arrow(self):
        a = self._arrow_model
        if a is not None:
            self._annotations.remove_drawing(a)
        self._arrow_model = None

class AtomPicker(MouseMode):
    name = 'select'
    '''
    Pick atom(s) only using the mouse, ignoring other objects. Select
    either from only a defined Atoms selection, or from all available
    AtomicStructure models.
    '''
    def __init__(self, session, atoms = None, all_models = False):
        if atoms is None and not all_models:
            raise TypeError('Must provide either an atom selection or\
                set all_models to True!')
        self._atoms = atoms
        self._add_trigger = None
        self._remove_trigger = None
        self.minimum_drag_pixels = 5
        self.drag_color = (0,255,0,255)
        self._drawn_rectangle = None
        if all_models:
            self.pick_from_any()
        
    def pick_from_any():
        self._choose_from_all_models = True
        self._update_atomic_models()
        self._add_handlers()
    
    def _update_atomic_models():
        self._atoms = None
        for m in self._session.models.list():
            if m.atomspec_has_atoms():
                if self_atoms is None:
                    self._atoms = m.atoms
                else:
                    self._atoms = self._atoms.merge(m.atoms)
    
    def pick_from_selection(atoms):
        from chimerax.core.atomic.molarray import Atoms
        if type(atoms) != Atoms:
            raise TypeError('Please provide an Atoms array as your selection!')
        self._choose_from_all_models = False
        self._atoms = atoms
        self._remove_handlers()
        
    def _add_handlers():
        self._add_trigger = self.session.triggers.add_handler(
                        'add models', self._update_atomic_models)
        self._remove_trigger = self.session.triggers.add_handler(
                        'remove models', self._update_atomic_models)
    
    def _remove_handlers():
        if self._add_trigger is not None:
            self.session.triggers.remove_handler(self._add_trigger)
            self._add_trigger = None
        if self._remove_trigger is not None:
            self.session.triggers.remove_handler(self._remove_trigger)
            self._remove_trigger = None
    
    def cleanup():
        self._remove_handlers()
        
    def mouse_down(self, event):
        MouseMode.mouse_down(self, event)
    
    def mouse_drag(self, event):
        if self._is_drag(event):
            self._undraw_drag_rectangle()
            self._draw_drag_rectangle(event)
    
    def mouse_up(self, event):
        self._undraw_drag_rectangle()
        if self._is_drag(event):
            # Select atoms in rectangle
            mouse_drag_select(self.mouse_down_position, event, self._session, self.view)
        
    def _is_drag(self, event):
        dp = self.mouse_down_position
        if dp is None:
            return False
        dx,dy = dp
        x, y = event.position()
        mp = self.minimum_drag_pixels
        return abs(x-dx) > mp or abs(y-dy) > mp

    def _draw_drag_rectangle(self, event):
        dx,dy = self.mouse_down_position
        x, y = event.position()
        v = self._session.main_view
        w,h = v.window_size
        v.draw_xor_rectangle(dx, h-dy, x, h-y, self.drag_color)
        self._drawn_rectangle = (dx,dy), (x,y)
        
    def _undraw_drag_rectangle(self):
        dr = self._drawn_rectangle
        if dr:
            (dx,dy), (x,y) = dr
            v = self._session.main_view
            w,h = v.window_size
            v.draw_xor_rectangle(dx, h-dy, x, h-y, self.drag_color)
            self._drawn_rectangle = None
    



    
class MouseModeRegistry():
    def __init__(self, session):
        self.session = session
        self._registered_modes = {}
        self._existing_modes = {}
    
    def get_names(self):
        return list(self._registered_modes.keys())
    
    def get_mode(self, name):
        if name in list(self._registered_modes.keys()):
            return self._registered_modes[name][0]
        else:
            raise Exception('Unrecognised mouse mode name!')
        
    def register_mode(self, name, mode, button, modifiers = []):
        mm = self.session.ui.mouse_modes
        existing_bindings = mm.bindings
        for b in existing_bindings:
            if b.exact_match(button, modifiers):
                self._existing_modes[name] = (b, button, modifiers)
            else:
                self._existing_modes[name] = (None)
        
        mm.bind_mouse_mode(button, modifiers, mode)
        self._registered_modes[name] = (mode, button, modifiers)
        
    def remove_mode(self, name):
        mm = self.session.ui.mouse_modes
        if self._existing_modes[name] is not None:
            mode, button, modifiers = self._existing_modes[name]
            mm.bind_mouse_mode(button, modifiers, mode)
            mode, button, modifiers = self._registered_modes[name]
            mode.cleanup()
        else:
            mode, button, modifiers = self._registered_modes[name]
            mm.remove_binding(button, modifiers)
            mode.cleanup()
        self._existing_modes.pop(name)
        self._registered_modes.pop(name)
