from chimerax.core.models import Model, Drawing
import numpy

class Rama_Annotator(Model):
    '''
    Doesn't actually draw anything, just colours the C-alphas of the parent
    atomic model according to Ramachandran score. But we'll create it as a model
    for consistency with other restraint/validation schemes (so it can be
    conveniently turned on/off via the Model panel).
    '''
    pickable = False

    def __init__(self, session, atomic_structure, hide_favored = False):
        Model.__init__(self, 'Ramachandran Validation', session)
        structure = self._atomic_structure = atomic_structure
        structure.add([self])
        from .. import molobject
        mgr = self._mgr = molobject.get_ramachandran_manager(session)
        self._ca_radius = 0.5
        self._prepare_drawings()
        # self._prepare_ca_display()
        self._hide_favored = hide_favored
        self.track_whole_model = True
        t = structure.triggers
        self._structure_change_handler = t.add_handler('changes', self._update_graphics_if_needed)

    @property
    def ca_radius(self):
        return self._ca_radius

    @ca_radius.setter
    def ca_radius(self, radius):
        if radius != self._ca_radius:
            self._ca_radius = radius
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
    def hide_favored(self):
        return self._hide_favored

    @hide_favored.setter
    def hide_favored(self, flag):
        if flag != self._hide_favored:
            self._hide_favored = flag
            self.update_graphics()

    @property
    def display(self):
        return Model.display.fget(self)

    @display.setter
    def display(self, flag):
        cflag = self.display
        Model.display.fset(self, flag)
        if flag and not cflag:
            # self._prepare_ca_display()
            self.update_graphics()
        # if cflag and not flag:
        #     self._revert_ca_display()

    def _prepare_drawings(self):
        if not hasattr(self, '_omega_drawing'):
            od = self._omega_drawing = Drawing('cis/twisted omegas')
            od.skip_bounds = True
            self.add_drawing(od)
        if not hasattr(self, '_rama_drawing'):
            rd = self._rama_drawing = Drawing('Ramachandran score indicators')
            rd.skip_bounds = True
            from chimerax.core.surface.shapes import sphere_geometry2
            rd.vertices, rd.normals, rd.triangles = sphere_geometry2(80)
            self.add_drawing(rd)




    # def _prepare_ca_display(self):
    #     from chimerax.core.atomic import Atom
    #     self._selected_ramas.ca_atoms.draw_modes = Atom.BALL_STYLE
    #
    # def _revert_ca_display(self):
    #     res = self._atomic_structure.residues
    #     ramas = self._mgr.get_ramas(res)
    #     cas = ramas.ca_atoms
    #     cs = cas.residues.atoms[cas.residues.atoms.names == 'C']
    #     cas.colors = cs.colors
    #     cas.draw_modes = cs.draw_modes


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
        # self._revert_ca_display()
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
                # self._prepare_ca_display()
                return
        reasons = changes.atom_reasons()
        if 'coord changed' in reasons:
            update_needed = True
        if 'display changed' in reasons or 'hide changed' in reasons:
            # self._prepare_ca_display()
            self._visible_ramas = self._selected_ramas[self._selected_ramas.visibles]
            update_needed = True
        # if 'color changed' in reasons:
        #     update_needed = True
        if update_needed:
            self.update_graphics()

    def update_graphics(self, *_):
        ramas = self._visible_ramas
        od = self._omega_drawing
        rd = self._rama_drawing
        if not len(ramas):
            od.display = False
            rd.display = False
            return
        mgr = self._mgr
        #mgr.color_cas_by_rama_score(ramas, self.hide_favored)
        coords, colors = mgr._ca_positions_and_colors(ramas, self.hide_favored)
        n = len(coords)
        if n > 0:
            xyzr = numpy.empty((n, 4), numpy.double)
            xyzr[:,:3] = coords
            xyzr[:,3] = self.ca_radius
            from chimerax.core.geometry import Places
            rd.positions = Places(shift_and_scale = xyzr)
            rd.colors = colors
            rd.display = True

        v, n, t, c = mgr._draw_cis_and_twisted_omegas(ramas)
        if len(v):
            od.vertices, od.normals, od.triangles, od.vertex_colors = v, n, t, c
            od.display = True
        else:
            od.display = False
