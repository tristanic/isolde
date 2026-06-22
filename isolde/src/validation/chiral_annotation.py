# @Author: Tristan Croll
# @Date:   20-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 20-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Real-time visual validation of chiral centres, mirroring
:class:`RotamerAnnotator`. Unlike the chiral *restraint* (a silent safety
barrier), this is a *validation* indicator: it lights up only on genuine
outliers -- centres whose geometry is wrong-handed or badly strained -- so the
rare cases that slip past the restraint (or exist in a static, non-simulating
model) are impossible to miss. It reads the centre geometry directly via
:func:`chimerax.isolde.atomic.chirality.chiral_outliers`, so it works with or
without a running simulation.
'''

import numpy
from chimerax.core.models import Model
from chimerax.geometry import (translation, rotation, Places, Place,
    vector_rotation, scale, norm)

from ..geometry import exclamation_mark

# Glyph sizing. The marker is drawn in scene (Angstrom) units, so by default it
# scales with the model. To keep it legible when scanning a zoomed-out model we
# enlarge it as the field of view widens: it stays at its natural size while the
# view is narrow (<= VIEW_WIDTH_MIN_A, so it scales with the model when zoomed in)
# and ramps linearly up to a cap (~ an alanine residue) by the time the view is
# wide (>= VIEW_WIDTH_MAX_A), so it can neither vanish nor balloon. All tunable.
GLYPH_BASE_A = 0.85  # natural glyph extent in Angstrom (~mark height + offset)
MAX_GLYPH_A = 4.0  # size cap (~ an alanine residue)
VIEW_WIDTH_MIN_A = 6.0  # view width at/below which the glyph is its natural size
VIEW_WIDTH_MAX_A = 40.0  # view width at/above which the glyph reaches the cap


class ChiralAnnotator(Model):
    '''
    Draws a 3D marker on every chiral centre whose handedness is wrong or badly
    strained, scaled and coloured by severity. Add as a child model of an
    :class:`AtomicStructure`; toggling its display off pauses validation.
    '''
    pickable = False

    def __init__(self, atomic_structure):
        structure = self._atomic_structure = atomic_structure
        session = structure.session
        Model.__init__(self, 'Chiral Validation', session)
        # Chirality is binary (wrong-handed or not), so every marker is the same
        # fixed size -- no severity scaling (cf. the cis/trans peptide markup).
        self._scale = 1.0
        # Cached outlier set (refreshed only on structure changes) + the pixel
        # size at which the glyph Places were last built (for the camera gate).
        self._current_chirals = None
        self._last_pixel_size = None
        d = self._drawing = self._chiral_indicator()
        self.add([d])
        structure.add([self])
        t = structure.triggers
        self._structure_change_handler = t.add_handler(
            'changes', self._update_graphics_if_needed
        )
        # Rescale the glyphs to a fixed minimum on-screen size as the camera zooms.
        self._graphics_update_handler = session.triggers.add_handler(
            'graphics update', self._rescale_for_camera
        )
        self.update_graphics()

    @property
    def display(self):
        return Model.display.fget(self)

    @display.setter
    def display(self, flag):
        cflag = self.display
        Model.display.fset(self, flag)
        if flag and not cflag:
            self.update_graphics()

    def delete(self):
        h = getattr(self, '_structure_change_handler', None)
        if h is not None and self._atomic_structure is not None \
                and not self._atomic_structure.deleted:
            self._atomic_structure.triggers.remove_handler(h)
        self._structure_change_handler = None
        gh = getattr(self, '_graphics_update_handler', None)
        if gh is not None:
            self.session.triggers.remove_handler(gh)
        self._graphics_update_handler = None
        Model.delete(self)

    def _update_graphics_if_needed(self, trigger_name, changes):
        if not self.visible:
            return
        changes = changes[1]
        update_needed = False
        if len(changes.created_atoms()):
            update_needed = True
        if changes.num_deleted_atoms():
            update_needed = True
        reasons = changes.atom_reasons()
        if 'coord changed' in reasons:
            update_needed = True
        if 'display changed' in reasons or 'hide changed' in reasons:
            update_needed = True
        if 'active_coordset changed' in changes.structure_reasons():
            update_needed = True
        if update_needed:
            from chimerax.atomic import get_triggers
            get_triggers().add_handler('changes done', self.update_graphics)

    def update_graphics(self, *_):
        # Structure-change path: re-detect outliers (the expensive part) and cache
        # them, then (re)build the glyph placements. The per-frame camera rescale
        # reuses the cached set so it never re-runs outlier detection.
        from chimerax.core.triggerset import DEREGISTER
        if not self.visible:
            return DEREGISTER
        from ..atomic.chirality import chiral_outliers
        chirals, oriented, severity = chiral_outliers(
            self.session, self._atomic_structure.atoms
        )
        # Only mark up centres whose atom is actually shown, consistent with
        # ISOLDE's other validation glyphs -- so hiding part of the model hides its
        # markup too. "Shown" is ribbon-aware, exactly as the Rama annotator's
        # Atom visible-ignoring-ribbon test (chiral.cpp):
        #     display() && !(hide() & ~HIDE_RIBBON)
        # i.e. displayed AND not hidden by anything other than the ribbon, so
        # backbone chiral centres keep their markup in the usual ribbon view but a
        # genuinely hidden atom loses it. Refreshed on display/hide changes (see
        # _update_graphics_if_needed). TODO: a per-validator option to keep outlier
        # markup visible even for hidden atoms.
        if len(chirals):
            cas = chirals.chiral_atoms
            shown = cas.displays & ((cas.hides & ~cas.HIDE_RIBBON) == 0)
            chirals = chirals[shown]
        self._current_chirals = chirals
        self._rebuild_places()
        return DEREGISTER

    def _glyph_scale(self, cc):
        # Ramp the glyph size with the field-of-view width: natural size while the
        # view is <= VIEW_WIDTH_MIN_A wide (so it scales with the model when zoomed
        # in), linearly up to MAX_GLYPH_A by VIEW_WIDTH_MAX_A (so it stays legible
        # when zoomed out), capped beyond. view_width() is evaluated at the centre's
        # scene coordinate so it is correct under perspective.
        vw = self.session.main_view.camera.view_width(cc.chiral_atom.scene_coord)
        if vw <= VIEW_WIDTH_MIN_A:
            eff_a = GLYPH_BASE_A
        elif vw >= VIEW_WIDTH_MAX_A:
            eff_a = MAX_GLYPH_A
        else:
            frac = (vw - VIEW_WIDTH_MIN_A) / (VIEW_WIDTH_MAX_A - VIEW_WIDTH_MIN_A)
            eff_a = GLYPH_BASE_A + frac * (MAX_GLYPH_A - GLYPH_BASE_A)
        return self._scale * eff_a / GLYPH_BASE_A

    def _rebuild_places(self):
        d = self._drawing
        chirals = self._current_chirals
        if chirals is None or not len(chirals):
            d.display = False
            return
        # Remember the representative pixel size these placements were built at, so
        # the camera handler can skip rebuilding while the zoom is unchanged.
        self._last_pixel_size = self.session.main_view.pixel_size()
        places = [
            self._marker_place(chirals[i], self._glyph_scale(chirals[i]))
            for i in range(len(chirals))
        ]
        d.positions = Places(places)
        d.display = True

    def _rescale_for_camera(self, *_):
        # Persistent 'graphics update' handler: rebuild the glyph placements only
        # when the zoom (pixel size) has actually changed, so an idle view is free.
        if not self.visible:
            return
        chirals = self._current_chirals
        if chirals is None or not len(chirals):
            return
        ps = self.session.main_view.pixel_size()
        last = self._last_pixel_size
        if last is not None and abs(ps - last) <= 1e-4 * last:
            return
        self._rebuild_places()

    @staticmethod
    def _unit(v, fallback=(0.0, 0.0, 1.0)):
        n = norm(v)
        if n < 1e-6:
            return numpy.asarray(fallback, dtype=float)
        return v / n

    def _perp(self, u):
        # An arbitrary unit vector perpendicular to u.
        a = numpy.array([1.0, 0, 0]) if abs(u[0]) < 0.9 else numpy.array([0, 1.0, 0])
        return self._unit(a - numpy.dot(a, u) * u)

    def _marker_place(self, cc, s):
        # Build the per-centre placement frame. The glyph's local axes map as:
        #   +z -> `out`  (exclamation mark, perpendicular to the smallest
        #                 substituent, pointing out between s2 and s3)
        #   +x -> `u4`   (arrows, along the smallest-substituent axis)
        atoms = cc.atoms
        C = numpy.asarray(atoms[0].coord, dtype=float)
        p2 = numpy.asarray(atoms[2].coord, dtype=float)
        p3 = numpy.asarray(atoms[3].coord, dtype=float)
        s4 = cc.fourth_substituent
        if s4 is not None:
            u4 = self._unit(numpy.asarray(s4.coord, dtype=float) - C)
        else:
            # Unmodelled 4th substituent (e.g. a stripped H): use the implied
            # tetrahedral direction opposite the sum of the other three bonds.
            p1 = numpy.asarray(atoms[1].coord, dtype=float)
            u4 = self._unit(-(self._unit(p1 - C) + self._unit(p2 - C)
                              + self._unit(p3 - C)))
        # `out`: the s2/s3 bisector, made perpendicular to u4.
        bis = self._unit(self._unit(p2 - C) + self._unit(p3 - C))
        out = self._unit(bis - numpy.dot(bis, u4) * u4, fallback=self._perp(u4))
        yax = self._unit(numpy.cross(out, u4), fallback=self._perp(u4))
        M = numpy.empty((3, 4))
        M[:, 0] = u4 * s
        M[:, 1] = yax * s
        M[:, 2] = out * s
        M[:, 3] = C
        return Place(matrix=M)

    def _chiral_indicator(self):
        # Composite glyph in a canonical local frame:
        #   * an exclamation mark along +z (oriented perpendicular to the smallest
        #     substituent), DOT pointing inward toward the chiral atom, shifted so
        #     its innermost point is >= 0.5 A from origin;
        #   * two short arrows along +x and -x (the smallest-substituent axis),
        #     emanating from the midpoint of the mark's barrel - +x = where the
        #     substituent currently is, -x = the opposite side it should swing to.
        from ..geometry import simple_arrow_geometry
        from .constants import rotarama_defaults as _rr
        ev, en, et = exclamation_mark(radius=0.1, height=0.6, nc=8)
        # Flip end-for-end (180 deg about x; a proper rotation, so winding/normals
        # stay correct) so the dot is at the inner end, then offset >= 0.25 A out.
        r180 = rotation((1, 0, 0), 180)
        r180.transform_points(ev, in_place=True)
        r180.transform_vectors(en, in_place=True)
        ev[:, 2] += 0.25 - ev[:, 2].min()
        mid_z = 0.5 * (ev[:, 2].min() + ev[:, 2].max())

        av, an, at = simple_arrow_geometry(radius=0.06, height=0.37, nc=8)
        r_px = rotation((0, 1, 0), 90)    # +z -> +x  (current substituent side)
        pxv = av.copy(); r_px.transform_points(pxv, in_place=True)
        pxn = an.copy(); r_px.transform_vectors(pxn, in_place=True)
        r_nx = rotation((0, 1, 0), -90)   # +z -> -x  (correct/should-go side)
        nxv = av.copy(); r_nx.transform_points(nxv, in_place=True)
        nxn = an.copy(); r_nx.transform_vectors(nxn, in_place=True)
        # Lift both arrows to the midpoint of the mark's barrel.
        shift = translation((0, 0, mid_z))
        shift.transform_points(pxv, in_place=True)
        shift.transform_points(nxv, in_place=True)

        v = numpy.concatenate((ev, pxv, nxv))
        n = numpy.concatenate((en, pxn, nxn))
        t = numpy.concatenate((et, at + len(ev), at + len(ev) + len(pxv)))

        # Fixed two-tone colouring (severity drives size, not colour): the mark and
        # the current-side arrow take the standard outlier colour; the correct/
        # should-go arrow takes the max-favoured colour.
        outlier = numpy.array(_rr.OUTLIER_COLOR, numpy.uint8)
        favored = numpy.array(_rr.MAX_FAVORED_COLOR, numpy.uint8)
        vc = numpy.empty((len(v), 4), numpy.uint8)
        vc[:len(ev) + len(pxv)] = outlier
        vc[len(ev) + len(pxv):] = favored

        d = Model('chiral indicator', self.session)
        d.skip_bounds = True
        d.pickable = False
        d.set_geometry(v, n, t)
        d.vertex_colors = vc
        return d

    def take_snapshot(self, session, flags):
        data = {
            'model state': Model.take_snapshot(self, session, flags),
            'structure': self._atomic_structure,
        }
        from .. import ISOLDE_STATE_VERSION
        data['version'] = ISOLDE_STATE_VERSION
        return data

    @staticmethod
    def restore_snapshot(session, data):
        ca = ChiralAnnotator(data['structure'])
        Model.set_state_from_snapshot(ca, session, data['model state'])
        return ca
