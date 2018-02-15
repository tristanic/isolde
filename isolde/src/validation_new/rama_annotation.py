from chimerax.core.models import Model

class Rama_Annotator(Model):
    '''
    Doesn't actually draw anything, just colours the C-alphas of the parent
    atomic model according to Ramachandran score. But we'll create it as a model
    for consistency with other restraint/validation schemes (so it can be
    conveniently turned on/off via the Model panel).
    '''
    pickable = False

    def __init__(self, session, atomic_structure):
        Model.__init__(self, 'Ramachandran Validation', session)
        structure = self._atomic_structure = atomic_structure
        structure.add([self])
        from .. import session_extensions
        mgr = self._mgr = session_extensions.get_ramachandran_manager(session)
        self.track_whole_model = True
        self._prepare_display()
        t = structure.triggers
        self._structure_change_handler = t.add_handler('changes', self._update_graphics_if_needed)
        self.update_graphics()

    @property
    def track_whole_model(self):
        return self._track_whole_model

    @track_whole_model.setter
    def track_whole_model(self, flag):
        self._track_whole_model = flag
        if flag:
            res = self._selected_residues = self._atomic_structure.residues
            ramas = self._selected_ramas = self._mgr.get_ramas(res)
            self._visible_ramas = ramas[ramas.visibles]
            self.update_graphics()

    @property
    def display(self):
        return Model.display.fget(self)

    @display.setter
    def display(self, flag):
        cflag = self.display
        Model.display.fset(self, flag)
        if flag and not cflag:
            self._prepare_display()
            self.update_graphics()
        if cflag and not flag:
            self._revert_display()

    def _prepare_display(self):
        from chimerax.core.atomic import Atom
        self._selected_ramas.ca_atoms.draw_modes = Atom.BALL_STYLE

    def _revert_display(self):
        res = self._atomic_structure.residues
        ramas = self._mgr.get_ramas(res)
        cas = ramas.ca_atoms
        cs = cas.residues.atoms[cas.residues.atoms.names == 'C']
        cas.colors = cs.colors
        cas.draw_modes = cs.draw_modes


    def restrict_to_selected_residues(self, residues):
        ''' Restrict validation to a defined set of residues. '''
        us = residues.unique_structures
        if not len(us):
            raise TypeError('No residues selected!')
        if len(us) !=1 or us[0] != self._atomic_structure:
            raise TypeError('All residues must be from the parent model!')
        res = self._selected_residues = residues
        ramas = self._selected_ramas = self._mgr.get_ramas(residues)
        self._visible_ramas = ramas[ramas.visibles]
        self.track_whole_model = False
        self.update_graphics()

    @property
    def color_scale(self):
        return self._mgr.color_scale

    def delete(self):
        self._revert_display()
        h = self._structure_change_handler
        if h is not None:
            self._atomic_structure.triggers.remove_handler(h)
        Model.delete(self)

    def _update_graphics_if_needed(self, trigger_name, changes):
        if not self.visible:
            return
        changes = changes[1]
        update_needed = False
        if (self._track_whole_model):
            '''
            Need to update the set of ramas if atoms are added. Deletions will
            take care of themselves.
            '''
            if len(changes.created_atoms()):
                # Trigger rebuild of rama array and graphics update
                self.track_whole_model = True
                self._prepare_display()
                return
        reasons = changes.atom_reasons()
        if 'coord changed' in reasons:
            update_needed = True
        if 'display changed' in reasons or 'hide changed' in reasons:
            self._prepare_display()
            self._visible_ramas = self._selected_ramas[self._selected_ramas.visibles]
            update_needed = True
        if update_needed:
            self.update_graphics()

    def update_graphics(self, *_):
        self._mgr.color_cas_by_rama_score(self._visible_ramas)
