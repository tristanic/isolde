'''
Drawing code for live annotation of rotamers
'''
import numpy
from math import log
from time import time

from chimerax.core.models import Drawing, Model
from chimerax.core.atomic import Bonds
from chimerax.core.geometry import translation, rotation, Places

from ..color import standard_three_color_scale
from ..geometry import exclamation_mark, spiral, bond_cylinder_placements
from ..geometry import scale_transforms

class Rotamer_Annotator(Model):
    ''' Model holding annotations for current allowed or outlier rotamers. '''

    pickable = False

    def __init__(self, session, atomic_structure):
        Model.__init__(self, 'Rotamer Validation', session)

        from .. import molobject
        mgr = self._mgr = molobject.get_rotamer_manager(session)
        structure = self._atomic_structure = atomic_structure
        self._MAX_SCALE = 2 # maximum scale factor for annotation drawings
        self._hide_favored = True
        d = self._drawing = self._rota_indicator()
        self.add_drawing(d)
        structure.add([self])
        self.track_whole_model = True
        t = structure.triggers
        self._structure_change_handler = t.add_handler('changes', self._update_graphics_if_needed)

        self.update_graphics()
        self._update_needed = False

    @property
    def display(self):
        return Model.display.fget(self)

    @display.setter
    def display(self, flag):
        cflag = self.display
        Model.display.fset(self, flag)
        if flag and not cflag:
            self.update_graphics()

    @property
    def track_whole_model(self):
        return self._track_whole_model

    @track_whole_model.setter
    def track_whole_model(self, flag):
        self._track_whole_model = flag
        if flag:
            res = self._selected_residues = self._atomic_structure.residues
            self._selected_rotamers = self._mgr.get_rotamers(res)
            self.update_graphics()

    def restrict_to_selected_residues(self, residues):
        ''' Restrict validation to a defined set of residues. '''
        us = residues.unique_structures
        if not len(us):
            raise TypeError('No residues selected!')
        if len(us) !=1 or us[0] != self._atomic_structure:
            raise TypeError('All residues must be from the parent model!')
        self._selected_residues = residues
        self._selected_rotamers = self._mgr.get_rotamers(residues)
        self.track_whole_model = False
        self.update_graphics()

    @property
    def color_scale(self):
        return self._mgr.color_scale

    @property
    def hide_favored(self):
        ''' Show annotations for favoured rotamers, or just non-favoured/outliers?'''
        return self._hide_favored

    @hide_favored.setter
    def hide_favored(self, flag):
        cflag = self._hide_favored
        self._hide_favored = flag
        if flag != cflag:
            self.update_graphics()

    def delete(self):
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
            Need to update the set of rotamers if atoms are added. Deletions will
            take care of themselves.
            '''
            if len(changes.created_atoms()):
                r = self._selected_residues = self._atomic_structure.residues
                self._selected_rotamers = self._mgr.get_rotamers(r)
                update_needed = True
        if changes.num_deleted_atoms():
            update_needed = True
        reasons = changes.atom_reasons()
        if 'coord changed' in reasons:
            update_needed = True
        if 'display changed' in reasons or 'hide changed' in reasons:
            update_needed = True
        if (update_needed):
            self.update_graphics()

    def update_graphics(self, *_, scale_by_scores = True):
        if not self.visible:
            return
        rots, scales, colors = self._mgr.validate_scale_and_color_rotamers(
            self._selected_rotamers, max_scale=self._MAX_SCALE,
            non_favored_only = self._hide_favored)
        d = self._drawing
        if not len(rots):
            d.display = False
            return
        d.display = True
        bonds = rots.ca_cb_bonds
        transforms = bond_cylinder_placements(bonds)
        if scale_by_scores:
            transforms = Places(place_array=scale_transforms(scales, transforms.array()))
        d.positions = transforms
        d.colors = colors

    def _rota_indicator(self):
        v1, n1, t1 = exclamation_mark(radius=0.1, height=0.5, nc = 8)
        v2, n2, t2 = spiral(major_radius=0.3, minor_radius = 0.05, height=0.4,
                            turn_segments=6, circle_segments=3)
        translation((0,0,-0.15)).move(v2)
        v = numpy.concatenate((v1, v2))
        n = numpy.concatenate((n1, n2))
        t = numpy.concatenate((t1, t2+len(v1)))
        r = rotation((1,0,0),180)
        r.move(v)
        r.move(n)
        translation((0,0.5,0.25)).move(v)
        d = Drawing('rotamer indicator')
        d.skip_bounds = True
        d.vertices, d.normals, d.triangles = v, n, t
        return d


    def _exclamation_mark(self):
        v, n, t = exclamation_mark(radius=0.1, height=0.5, nc = 8)
        rotation((1,0,0),180).move(v)
        translation((0,1,0)).move(v)
        d = Drawing('rotamer cb indicator')
        d.vertices, d.normals, d.triangles = v, n, t
        return d

    def _cb_annotation(self):
        from chimerax.core.surface.shapes import cylinder_geometry
        v, n, t = cylinder_geometry(radius=0.1, height=2, nc=8, caps=True)
        d = Drawing('rotamer cb indicator')
        d.vertices, d.normals, d.triangles = v, n, t
        return d
