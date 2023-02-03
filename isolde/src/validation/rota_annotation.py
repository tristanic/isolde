# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 29-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



import numpy
from math import log
from time import time

from chimerax.core.models import Drawing, Model
from chimerax.atomic import Bonds
from chimerax.geometry import translation, rotation, Places

from ..geometry import exclamation_mark, spiral, bond_cylinder_placements
from ..geometry import scale_transforms

class RotamerAnnotator(Model):
    '''
    Handles the task of real-time validation of rotamers for a single
    :py:class:`chimerax.AtomicStructure` and drawing of 3D indicators of their
    current scores. Designed to be mostly "fire-and-forget":

    .. code-block:: python

        ra = RotamerAnnotator(atomic_structure)

    adds the :py:class`RotamerAnnotator` as a child model to `atomic_structure`
    and will update its validation drawing every time the coordinates change.
    Alternatively:

    .. code-block:: python

        from chimerax.isolde import session_extensions as sx
        ra = sx.get_rota_annotator(atomic_model)

    creates the :py:class:`RotamerAnnotator` if it doesn't exist, or returns
    the existing one if it does.

    Each rotamer score visualised as a 3D exclamation mark surrounded by a
    spiral, which changes colour and grows with outlier severity. By default
    only non-favoured rotamers are flagged.

    Turning off the display of the :py:class:`RotamerAnnotator` model (e.g.
    via the `ChimeraX` Model Panel) temporarily turns off automatic validation,
    which will restart when display is turned back on.
    '''
    pickable = False

    def __init__(self, atomic_structure):
        '''
        Create the validator object, and add it as a child model to the target
        structure.

        Args:
            * atomic_structure:
                - a :py:class:`ChimeraX.AtomicStructure` instance
        '''
        structure = self._atomic_structure = atomic_structure
        session = structure.session
        Model.__init__(self, 'Rotamer Validation', session)

        from .. import molobject
        mgr = self._mgr = molobject.get_rotamer_mgr(session)
        self._MAX_SCALE = 2 # maximum scale factor for annotation drawings
        self._hide_favored = True
        d = self._drawing = self._rota_indicator()
        self.add([d])
        structure.add([self])
        self.track_whole_model = True
        t = structure.triggers
        self._structure_change_handler = t.add_handler('changes', self._update_graphics_if_needed)

        self.update_graphics()
        self._update_needed = False

    @property
    def display(self):
        '''
        Show/hide the validation markup (automatic validation will pause while
        hidden).
        '''
        return Model.display.fget(self)

    @display.setter
    def display(self, flag):
        cflag = self.display
        Model.display.fset(self, flag)
        if flag and not cflag:
            self.update_graphics()

    @property
    def track_whole_model(self):
        '''
        Tell the validator to track/annotate all rotameric residues in the
        model (the default starting state).
        '''
        return self._track_whole_model

    @track_whole_model.setter
    def track_whole_model(self, flag):
        self._track_whole_model = flag
        if flag:
            res = self._selected_residues = self._atomic_structure.residues
            self._selected_rotamers = self._mgr.get_rotamers(res)
            self.update_graphics()

    def restrict_to_selected_residues(self, residues):
        '''
        Restrict validation to a defined set of residues. Use
        :py:attr:`track_whole_model` = True to once again cover all residues.

        Args:
            * residues:
                - A :py:class:`chimerax.Residues` instance
        '''
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
        ''' Returns a 3-tuple of (r,g,b,a) arrays defining the current colour scale.'''
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
        if h is not None and self._atomic_structure is not None:
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
            created = changes.created_atoms()
            if len(created):
                # Only need to update if we've added new non-hydrogen protein atoms
                from chimerax.atomic import Residue
                ur = created[created.element_names!='H'].unique_residues
                if sum(ur.polymer_types==Residue.PT_AMINO):
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
            from chimerax.atomic import get_triggers
            get_triggers().add_handler('changes done', self.update_graphics)

    def update_graphics(self, *_, scale_by_scores = True):
        from chimerax.core.triggerset import DEREGISTER
        if not self.visible:
            return DEREGISTER
        rots, scales, colors = self._mgr.validate_scale_and_color_rotamers(
            self._selected_rotamers, max_scale=self._MAX_SCALE,
            non_favored_only = self._hide_favored)
        d = self._drawing
        if not len(rots):
            d.display = False
            return DEREGISTER
        d.display = True
        bonds = rots.ca_cb_bonds
        transforms = bond_cylinder_placements(bonds)
        if scale_by_scores:
            transforms = Places(place_array=scale_transforms(scales, transforms.array()))
        d.positions = transforms
        d.colors = colors
        return DEREGISTER

    def _rota_indicator(self):
        v1, n1, t1 = exclamation_mark(radius=0.1, height=0.5, nc = 8)
        v2, n2, t2 = spiral(major_radius=0.3, minor_radius = 0.05, height=0.4,
                            turn_segments=6, circle_segments=3)
        translation((0,0,-0.15)).transform_points(v2, in_place=True)
        v = numpy.concatenate((v1, v2))
        n = numpy.concatenate((n1, n2))
        t = numpy.concatenate((t1, t2+len(v1)))
        r = rotation((1,0,0),180)
        r.transform_points(v, in_place=True)
        r.transform_vectors(n, in_place=True)
        translation((0,0.5,0.25)).transform_points(v, in_place=True)
        d = Model('rotamer indicator', self.session)
        d.skip_bounds = True
        d.pickable = False
        d.set_geometry(v, n, t)
        return d


    def _exclamation_mark(self):
        v, n, t = exclamation_mark(radius=0.1, height=0.5, nc = 8)
        flip = rotation((1,0,0), 180)
        flip.transform_points(v, in_place=True)
        flip.transform_vectors(n, in_place=True)
        translation((0,1,0)).transform_points(v, in_place=True)
        d = Drawing('rotamer cb indicator')
        d.set_geometry(v, n, t)
        return d

    def _cb_annotation(self):
        from chimerax.surface.shapes import cylinder_geometry
        v, n, t = cylinder_geometry(radius=0.1, height=2, nc=8, caps=True)
        d = Drawing('rotamer cb indicator')
        d.set_geometry(v, n, t)
        return d


    def take_snapshot(self, session, flags):
        from chimerax.core.models import Model
        data = {
            'model state': Model.take_snapshot(self, session, flags),
            'structure': self._atomic_structure,
        }
        from .. import ISOLDE_STATE_VERSION
        data['version']=ISOLDE_STATE_VERSION
        return data

    @staticmethod
    def restore_snapshot(session, data):
        from chimerax.core.models import Model
        ra = RotamerAnnotator(data['structure'])
        Model.set_state_from_snapshot(ra, session, data['model state'])
        return ra
