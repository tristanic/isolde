# @Author: Tristan Croll <tic20>
# @Date:   20-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 14-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



# Useful mouse modes for working with ISOLDE

from chimerax.mouse_modes import MouseMode, SelectMouseMode

from .tugging import TugAtomsMode

class AtomPicker(SelectMouseMode):

    def __init__(self, session, isolde):
        super().__init__(session)
        self.isolde = isolde

    def _pick_exclude(self, d):
        restrict_model = self._restriction_model
        if not getattr(d, 'pickable', True):
            return True
        if restrict_model is None:
            return False
        from chimerax.core.models import Model
        if isinstance(d, Model):
            if restrict_model in d.all_models():
                return False
            return True
        if restrict_model in d.drawing_lineage:
            return False
        return True

    def _find_restrictions(self):
        isolde = self.isolde
        restrict_model = isolde.selected_model
        if isolde.simulation_running:
            restrict_atoms = isolde.sim_manager.sim_construct.mobile_atoms
        else:
            restrict_atoms = None
        return restrict_model, restrict_atoms

    def mouse_up(self, event):
        self._undraw_drag_rectangle()
        mode = self.mode
        # if event.shift_down() and mode == 'replace':
        #     mode = 'add'
        # elif event.alt_down() and mode == 'replace':
        #     mode = 'subtract'

        if self._is_drag(event):
            # Select objects in rectangle
            self.mouse_drag_select(self.mouse_down_position, event, mode, self.session, self.view)
        elif not self.double_click:
            # Select object under pointer
            self.mouse_select(event, mode, self.session, self.view)
        MouseMode.mouse_up(self, event)

    def mouse_select(self, event, mode, session, view):
        x,y = event.position()
        rm, ra = self._find_restrictions()
        self._restriction_model = rm
        pick = view.first_intercept(x, y, self._pick_exclude)
        if ra is not None:
            if hasattr(pick, 'atom'):
                if pick.atom not in ra:
                    pick = None
                    # delattr(pick, 'atom')
            if hasattr(pick, 'residue') and pick.residue is not None:
                if not any([a in ra for a in pick.residue.atoms]):
                    pick = None
                    # delattr(pick, 'residue')
            if hasattr(pick, 'bond'):
                b = pick.bond
                if not any([a in ra for a in b.atoms]):
                    pick = None
        if hasattr(pick, 'bond') and pick.bond is not None:
            distance = pick.distance
            b = pick.bond
            pick = [pick]
            from chimerax.atomic import PickedAtom
            if b.atoms is not None:
                if ra is not None:
                    pick.extend([PickedAtom(a, distance) for a in b.atoms if a in ra])
                else:
                    pick.extend([PickedAtom(a, distance) for a in b.atoms])

        select_pick(session, pick, mode)

    def mouse_drag_select(self, start_xy, event, mode, session, view):
        from chimerax.atomic import Atoms, Bonds, Residues
        rm, ra = self._find_restrictions()
        self._restriction_model = rm
        sx, sy = start_xy
        x, y = event.position()
        pick = view.rectangle_intercept(sx, sy, x, y, exclude=self._pick_exclude)
        if ra is not None:
            rr = ra.unique_residues
            rb = ra.intra_bonds
            for p in reversed(pick):
                if hasattr(p, 'atoms'):
                    p.atoms = p.atoms.intersect(ra)
                if hasattr(p, 'bonds'):
                    p.bonds = p.bonds.intersect(rb)
                if hasattr(p, 'residues'):
                    p.residues = p.residues.intersect(rr)
        select_pick(session, pick, mode)

    def cleanup(self):
        pass



def select_pick(session, pick, mode = 'replace'):
    sel = session.selection
    from chimerax.core.undo import UndoState
    undo_state = UndoState("select")
    sel.undo_add_selected(undo_state, False)
    if pick is None:
        if mode == 'replace':
            from chimerax.core.commands import run
            run(session, 'select clear')
            session.logger.status('cleared selection')
    else:
        if mode == 'replace':
            sel.clear()
            mode = 'add'
        if isinstance(pick, list):
            for p in pick:
                p.select(mode)
        else:
            pick.select(mode)
    sel.clear_promotion_history()
    sel.undo_add_selected(undo_state, True, old_state=False)
    session.undo.register(undo_state)



class AtomPicker_Old(MouseMode):
    name = 'select'
    '''
    Pick atom(s) only using the mouse, ignoring other objects. Select
    either from only a defined Atoms selection, or from all available
    AtomicStructure models.

    Standard binding is ctrl-click (and drag)

    Modifiers:
        - ctrl-shift: Add to current selection
        - ctrl-alt: Subtract from current selection

        If shift is pressed, alt will be ignored.
    '''
    def __init__(self, session, isolde, all_models = False):
        MouseMode.__init__(self, session)

        self.mode = {'select': 'replace',
                     'select add': 'add',
                     'select subtract': 'subtract'}[self.name]
        self._isolde = isolde
        self._atoms = None
        self._add_trigger = None
        self._remove_trigger = None
        self.minimum_drag_pixels = 5
        self.drag_color = (0,255,0,255)
        self._drawn_rectangle = None
        self._choose_from_all_models = all_models
        if all_models:
            self.pick_from_any()
        else:
            self._isolde_model_handler = isolde.triggers.add_handler(
                'selected model changed', self._isolde_changed_model)
        # While a simulation is running, we only want to pick from the
        # mobile atoms
        self._sim_start_handler = isolde.triggers.add_handler(
            'simulation started', self._on_sim_start)
        self._sim_end_handler = isolde.triggers.add_handler(
            'simulation terminated', self._on_sim_end)

    @property
    def atoms(self):
        if not self._choose_from_all_models:
            if self._atoms is None:
                return self._isolde._selected_model.atoms
            return self._atoms
        from chimerax.atomic import Atoms
        ret = Atoms()
        for m in self.session.models.list():
            if m.atomspec_has_atoms():
                ret = ret.merge(m.atoms)
        return ret

    def _isolde_changed_model(self, trigger_name, m):
        self._atoms = None

    def _on_sim_start(self, *_):
        self._atoms = self._isolde.sim_manager.sim_construct.mobile_atoms

    def _on_sim_end(self, *_):
        if self._choose_from_all_models:
            self.pick_from_any()
        else:
            self._atoms = None

    def pick_from_any(self):
        self._choose_from_all_models = True
        self._update_atomic_models()
        self._add_handlers()

    def _update_atomic_models(self, *_):
        self._atoms = None

    def pick_from_selection(atoms):
        from chimerax.atomic import Atoms
        if type(atoms) != Atoms:
            raise TypeError('Please provide an Atoms array as your selection!')
        self._choose_from_all_models = False
        self._atoms = atoms
        self._remove_handlers()

    def _add_handlers(self):
        self._add_trigger = self.session.triggers.add_handler(
                        'add models', self._update_atomic_models)
        self._remove_trigger = self.session.triggers.add_handler(
                        'remove models', self._update_atomic_models)

    def _remove_handlers(self):
        if self._add_trigger is not None:
            self.session.triggers.remove_handler(self._add_trigger)
            self._add_trigger = None
        if self._remove_trigger is not None:
            self.session.triggers.remove_handler(self._remove_trigger)
            self._remove_trigger = None

    def cleanup(self):
        self._remove_handlers()

    def mouse_down(self, event):
        MouseMode.mouse_down(self, event)

    def mouse_drag(self, event):
        if self._is_drag(event):
            self._undraw_drag_rectangle()
            self._draw_drag_rectangle(event)

    def mouse_up(self, event):
        self._undraw_drag_rectangle()
        # if self._atoms is None:
        #     return
        if self._is_drag(event):
            # Select atoms in rectangle
            self._mouse_drag_select(self.mouse_down_position, event)
        else:
            # Select closest atom to a ray projecting from the pointer
            self._mouse_select(event)

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
        v = self.session.main_view
        w,h = v.window_size
        v.draw_xor_rectangle(dx, h-dy, x, h-y, self.drag_color)
        self._drawn_rectangle = (dx,dy), (x,y)

    def _undraw_drag_rectangle(self):
        dr = self._drawn_rectangle
        if dr:
            (dx,dy), (x,y) = dr
            v = self.session.main_view
            w,h = v.window_size
            v.draw_xor_rectangle(dx, h-dy, x, h-y, self.drag_color)
            self._drawn_rectangle = None


    def _mouse_drag_select(self, start_xy, event):
        session = self.session
        mode = self.mode
        view = self.view
        sx, sy = start_xy
        x, y = event.position()
        # Pick all objects under the rectangle
        pick = view.rectangle_intercept(sx, sy, x, y)
        # ... then weed out the non-atomic ones
        #add_toggle = event.shift_down()
        #sub_toggle = event.alt_down()
        if pick is None:
            if not toggle and not sub_toggle:
                session.selection.clear()
                session.logger.status('cleared selection')
                return
        # If shift is down, add to an existing selection, otherwise clear the
        # old selection
        #if not add_toggle and not sub_toggle:
        if mode == 'replace':
            session.selection.clear()
        for p in pick:
            if self._choose_from_all_models:
                if hasattr(p, 'atoms') or hasattr(p, 'residues'):
                    #if add_toggle or not sub_toggle:
                    if not mode == 'subtract':
                        p.select()
                    else:
                        if hasattr(p, 'atoms'):
                            p.atoms.selected = False
                            p.atoms.intra_bonds.selected=False
                        else:
                            p.residues.atoms.selected = False
            else:
                import numpy
                if hasattr(p, 'atoms'):
                    patoms = p.atoms.filter(self.atoms.indices(p.atoms) != -1)
                    if len(patoms):
                        if not mode == 'subtract':
                            patoms.selected = True
                            patoms.intra_bonds.selected = True
                        else:
                            patoms.selected = False
                            patoms.intra_bonds.selected = False
                elif hasattr(p, 'residues'):
                    presidues = p.residues.filter(self.atoms.unique_residues.indices(p.residues) != -1)
                    if len(presidues):
                        if not mode == 'subtract':
                            presidues.atoms.selected = True
                            presidues.atoms.intra_bonds.selected = True
                        else:
                            presidues.atoms.selected = False
                            presidues.atoms.intra_bonds.selected = False


    def _mouse_select(self, event):
        session = self.session
        mode = self.mode
        from . import picking
        x, y = event.position()
        # Allow atoms within 0.5 Angstroms of the cursor to be picked
        cutoff = 0.5
        picked_atom = picking.pick_closest_to_line(session, x, y, self.atoms,
            cutoff, hydrogens = True)
        # If shift is down, add to an existing selection, otherwise clear the
        # old selection
        add_toggle = event.shift_down()
        sub_toggle = event.alt_down()
        if mode == 'replace':
            session.selection.clear()
        if picked_atom is None:
            if mode == 'replace':
                session.logger.status('cleared selection')
            return
        if not mode == 'subtract':
            picked_atom.selected = True
        else:
            picked_atom.selected = False
            for b in picked_atom.bonds:
                b.selected = False


class AtomPickAdd(AtomPicker):
    name = 'select add'
    icon = None

class AtomPickSubtract(AtomPicker):
    name = 'select subtract'
    icon = None


class MouseModeRegistry():
    def __init__(self, session, isolde):
        self._isolde = isolde
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
                break
            else:
                self._existing_modes[name] = (None)

        mm.bind_mouse_mode(button, modifiers, mode)
        self._registered_modes[name] = (mode, button, modifiers)
        return mode

    def register_all_isolde_modes(self):
        # Button, modifier(s), name, class, args
        # from chimerax.clipper.mousemodes import ZoomMouseMode, ClipPlaneAdjuster
        session = self.session
        isolde = self._isolde
        standard_modes = (
            ('left', ['control',], 'select', AtomPicker, (session, isolde)),
            ('left', ['control','shift'], 'select add', AtomPickAdd, (session, isolde)),
            ('left', ['control','alt'], 'select subtract', AtomPickSubtract, (session, isolde)),
            )
        for m in standard_modes:
            self.register_mode(m[2], m[3](*m[4]), m[0], m[1])

        # zoom_mode = self.register_mode('zoom', ZoomMouseMode(session), 'right', ['shift'])
        # self.register_mode('adjust clip', ClipPlaneAdjuster(session, zoom_mode), 'wheel', ['shift'])


    def remove_mode(self, name):
        mm = self.session.ui.mouse_modes
        if self._existing_modes[name] is not None:
            mode, button, modifiers = self._existing_modes[name]
            mm.bind_mouse_mode(button, modifiers, mode.mode)
            mode, *_ = self._registered_modes[name]
            mode.cleanup()
        else:
            mode, button, modifiers = self._registered_modes[name]
            mm.remove_binding(button, modifiers)
            mode.cleanup()
        self._existing_modes.pop(name)
        self._registered_modes.pop(name)

    def remove_all_modes(self):
        names = list(self.get_names())
        for name in names:
            self.remove_mode(name)
