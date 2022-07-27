# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 29-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



from chimerax.core.models import Model, Drawing
import numpy

class RamaAnnotator(Model):
    '''
    Handles the task of real-time Ramachandran validation for a single
    :py:class:`chimerax.AtomicStructure` and visualisation of the results.
    Results are displayed as coloured spheres overlaying the alpha carbon
    atoms, shading from (by default) green through yellow to hot pink as
    the Ramachandran score goes from favoured through allowed to outlier.

    Designed to be mostly "fire-and-forget":

    .. code-block:: python

        ra = RamaAnnotator(atomic_structure)

    adds the :py:class`RamaAnnotator` as a child model to `atomic_structure`
    and will update its validation drawing every time the coordinates change.
    Alternatively:

    .. code-block:: python

        from chimerax.isolde import session_extensions as sx
        ra = sx.get_rama_annotator(atomic_model)

    creates the :py:class:`RamaAnnotator` if it doesn't exist, or returns
    the existing one if it does.

    Turning off the display of the :py:class:`RamaAnnotator` model (e.g.
    via the `ChimeraX` Model Panel) temporarily turns off automatic validation,
    which will restart when display is turned back on.
    '''
    pickable = False

    def __init__(self, atomic_structure, hide_favored = False,
        ignore_ribbon_hides = True):
        '''
        Create the validator object, and add it as a child model to the target
        structure.

        Args:
            * atomic_structure:
                - a :py:class:`ChimeraX.AtomicStructure` instance
            * hide_favored:
                - if True, indicators will only appear for non-favored residues.
                  (this can be changed at any time later)
            * ignore_ribbon_hides:
                - if True, the :attr:`chimerax.HIDE_RIBBON` bit will be ignored
                  when deciding when to display a particular annotation. This is
                  needed because the default ChimeraX ribbon representation
                  hides the CA for any residue whose sidechain is un-displayed
                  or missing (meaning that the glycine CA is always hidden).
        '''
        structure = self._atomic_structure = atomic_structure
        session = structure.session
        Model.__init__(self, 'Ramachandran Validation', session)
        structure.add([self])
        self._ignore_ribbon_hides = ignore_ribbon_hides
        self._hide_favored = hide_favored
        from .. import molobject
        mgr = self._mgr = molobject.get_ramachandran_mgr(session)
        self._ca_radius = 0.5
        self._prepare_drawings()
        # self._prepare_ca_display()
        self.track_whole_model = True
        t = structure.triggers
        self._structure_change_handler = t.add_handler('changes', self._update_graphics_if_needed)

    @property
    def ca_radius(self):
        '''
        Sets the radius (in Angstroms) of the sphere overlaying each CA atom.
        '''
        return self._ca_radius

    @ca_radius.setter
    def ca_radius(self, radius):
        if radius != self._ca_radius:
            self._ca_radius = radius
            self.update_graphics()

    @property
    def track_whole_model(self):
        '''
        Tell the validator to track/annotate all protein residues in the
        model (the default starting state).
        '''
        return self._track_whole_model

    @track_whole_model.setter
    def track_whole_model(self, flag):
        self._track_whole_model = flag
        if flag:
            res = self._selected_residues = self._atomic_structure.residues
            ramas = self._selected_ramas = self._mgr.get_ramas(res)
            self._update_visible_ramas()
            self.update_graphics()

    @property
    def hide_favored(self):
        ''' Show annotations for favoured rotamers, or just non-favoured/outliers?'''
        return self._hide_favored

    @hide_favored.setter
    def hide_favored(self, flag):
        if flag != self._hide_favored:
            self._hide_favored = flag
            self.update_graphics()

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
            # self._prepare_ca_display()
            self.update_graphics()
        # if cflag and not flag:
        #     self._revert_ca_display()

    def _update_visible_ramas(self):
        ramas = self._selected_ramas
        if self._ignore_ribbon_hides:
            self._visible_ramas = ramas[ramas.visibles_ignoring_ribbon]
        else:
            self._visible_ramas = ramas[ramas.visibles]


    def _prepare_drawings(self):
        if not hasattr(self, '_omega_drawing'):
            od = self._omega_drawing = Drawing('cis/twisted omegas')
            od.skip_bounds = True
            od.pickable = False
            self.add_drawing(od)
        if not hasattr(self, '_rama_drawing'):
            rd = self._rama_drawing = Drawing('Ramachandran score indicators')
            rd.skip_bounds = True
            rd.pickable = False
            from chimerax.surface.shapes import sphere_geometry2
            rd.set_geometry(*sphere_geometry2(80))
            self.add_drawing(rd)

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
        res = self._selected_residues = residues
        ramas = self._selected_ramas = self._mgr.get_ramas(residues)
        self._update_visible_ramas()
        self.track_whole_model = False
        self.update_graphics()

    @property
    def color_scale(self):
        ''' Returns a 3-tuple of (r,g,b,a) arrays defining the current colour scale.'''
        return self._mgr.color_scale

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
            Need to update the set of ramas if atoms are added. Deletions will
            take care of themselves.
            '''
            created = changes.created_atoms()
            if len(created):
                from chimerax.atomic import Residue
                # Only need to update if the additions are non-hydrogen protein atoms
                ur = created[created.element_names !='H'].unique_residues
                if any(ur.polymer_types == Residue.PT_AMINO):
                    # Trigger rebuild of rama array and graphics update
                    self.track_whole_model = True
                    # self._prepare_ca_display()
                    return
        reasons = changes.atom_reasons()
        if 'coord changed' in reasons:
            update_needed = True
        if 'display changed' in reasons or 'hide changed' in reasons:
            # self._prepare_ca_display()
            self._update_visible_ramas()
            update_needed = True
        if 'selected changed' in reasons:
            update_needed = True
        # if 'color changed' in reasons:
        #     update_needed = True
        if update_needed:
            from chimerax.atomic import get_triggers
            get_triggers().add_handler('changes done', self.update_graphics)

    def update_graphics(self, *_):
        from chimerax.core.triggerset import DEREGISTER
        ramas = self._visible_ramas
        od = self._omega_drawing
        rd = self._rama_drawing
        if not len(ramas):
            od.display = False
            rd.display = False
            return DEREGISTER
        mgr = self._mgr
        #mgr.color_cas_by_rama_score(ramas, self.hide_favored)
        coords, colors, selecteds = mgr._ca_positions_colors_and_selecteds(ramas, self.hide_favored)
        n = len(coords)
        if n > 0:
            xyzr = numpy.empty((n, 4), numpy.float32)
            xyzr[:,:3] = coords
            xyzr[:,3] = self.ca_radius
            from chimerax.geometry import Places
            rd.positions = Places(shift_and_scale = xyzr)
            rd.colors = colors
            rd.highlighted_positions = selecteds

            rd.display = True
        else:
            rd.display = False

        v, n, t, c = mgr._draw_cis_and_twisted_omegas(ramas)
        if len(v):
            od.set_geometry(v, n, t)
            od.vertex_colors = c
            od.display = True
        else:
            od.display = False
        return DEREGISTER

    def take_snapshot(self, session, flags):
        from chimerax.core.models import Model
        data = {
            'model state': Model.take_snapshot(self, session, flags),
            'structure': self._atomic_structure,
            'hide favored': self._hide_favored,
            'ignore ribbon hides': self._ignore_ribbon_hides
        }
        from .. import ISOLDE_STATE_VERSION
        data['version']=ISOLDE_STATE_VERSION
        return data

    @staticmethod
    def restore_snapshot(session, data):
        from chimerax.core.models import Model
        ra = RamaAnnotator(data['structure'], data['hide favored'], data['ignore ribbon hides'])
        Model.set_state_from_snapshot(ra, session, data['model state'])
        return ra
