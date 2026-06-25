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
from ._rs_glyphs import letter_R_geometry, letter_S_geometry

# R/S configuration labels are purely informational (not warnings), so unlike the
# outlier glyph they are NOT view-width-scaled: a fixed, unobtrusive height in
# Angstrom, sitting just beside the chiral atom.
LABEL_HEIGHT_A = 0.5
# Default colour mode for the labels (see ChiralAnnotator.set_label_color):
#   'auto'     -- a fixed tone chosen to contrast with the current background
#   'fromatoms'-- each letter takes its chiral atom's colour
#   else a stored uint8 RGBA (a custom ColorArg colour).
DEFAULT_LABEL_COLOR_MODE = 'auto'

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
        # Opt-in R/S configuration labels (camera-facing 3D letters), drawn only on
        # centres the user has flagged via ChiralCenter.label -- independent of the
        # always-on outlier glyph. 'R' and 'S' are distinct geometry, so each gets
        # its own drawing (instanced across its centres via its own Places).
        self._current_labelled = None
        self._current_configs = []
        self._r_centres = []
        self._s_centres = []
        self._last_view_dir = None
        self._last_bg = None
        # Label colour: mode in {'auto', 'fromatoms'} or a uint8 RGBA array.
        self._label_color_mode = DEFAULT_LABEL_COLOR_MODE
        self._r_drawing = self._make_letter_drawing('R', letter_R_geometry())
        self._s_drawing = self._make_letter_drawing('S', letter_S_geometry())
        self.add([d, self._r_drawing, self._s_drawing])
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
        from ..atomic.chirality import chiral_outliers, as_modelled_configs
        from chimerax.isolde.molobject import get_chiral_mgr
        atoms = self._atomic_structure.atoms
        chirals, oriented, severity = chiral_outliers(self.session, atoms)
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
        self._current_chirals = self._filter_shown(chirals)
        # Label layer: every centre the user opted in via ChiralCenter.label,
        # regardless of outlier status. chiral_outliers() has just ensured the
        # definitions exist, so the manager's centre set is complete here.
        mgr = get_chiral_mgr(self.session)
        allc = mgr.get_chirals(atoms, create=True)
        labelled = allc[allc.labels] if len(allc) else allc
        labelled = self._filter_shown(labelled)
        self._current_labelled = labelled
        self._current_configs = (
            as_modelled_configs(self.session, labelled) if len(labelled) else []
        )
        self._rebuild_places()
        return DEREGISTER

    def _filter_shown(self, chirals):
        '''Restrict to centres whose chiral atom is visible, ribbon-aware (see the
        Rama annotator's ``display() && !(hide() & ~HIDE_RIBBON)``).'''
        if chirals is None or not len(chirals):
            return chirals
        cas = chirals.chiral_atoms
        shown = cas.displays & ((cas.hides & ~cas.HIDE_RIBBON) == 0)
        return chirals[shown]

    def _eff_extent_a(self, cc):
        # Effective glyph extent in Angstrom, ramped with the field-of-view width:
        # natural size (GLYPH_BASE_A) while the view is <= VIEW_WIDTH_MIN_A wide (so
        # it scales with the model when zoomed in), linearly up to MAX_GLYPH_A by
        # VIEW_WIDTH_MAX_A (so it stays legible when zoomed out), capped beyond.
        # view_width() is evaluated at the centre's scene coordinate so it is
        # correct under perspective.
        vw = self.session.main_view.camera.view_width(cc.chiral_atom.scene_coord)
        if vw <= VIEW_WIDTH_MIN_A:
            return GLYPH_BASE_A
        if vw >= VIEW_WIDTH_MAX_A:
            return MAX_GLYPH_A
        frac = (vw - VIEW_WIDTH_MIN_A) / (VIEW_WIDTH_MAX_A - VIEW_WIDTH_MIN_A)
        return GLYPH_BASE_A + frac * (MAX_GLYPH_A - GLYPH_BASE_A)

    def _glyph_scale(self, cc):
        # Place-scale for the (GLYPH_BASE_A-sized) exclamation/arrow outlier glyph.
        return self._scale * self._eff_extent_a(cc) / GLYPH_BASE_A

    def _rebuild_places(self):
        # Full rebuild (structure-change path): both layers.
        self._rebuild_outlier_places()
        self._rebuild_label_places()

    def _rebuild_outlier_places(self):
        # Remember the pixel size the (view-width-scaled) outlier glyphs were built
        # at, so the camera handler can skip rebuilding while the zoom is unchanged.
        self._last_pixel_size = self.session.main_view.pixel_size()
        d = self._drawing
        chirals = self._current_chirals
        if chirals is None or not len(chirals):
            d.display = False
            return
        d.positions = Places([
            self._marker_place(chirals[i], self._glyph_scale(chirals[i]))
            for i in range(len(chirals))
        ])
        d.display = True

    def _rebuild_label_places(self):
        # Partition the labelled centres into 'R' and 'S' (their as-modelled config)
        # and build fixed-size, camera-facing placements for each letter drawing.
        # Remember the camera direction these billboards were oriented for.
        self._last_view_dir = self._camera_view_direction()
        labelled = self._current_labelled
        configs = self._current_configs
        r_places = []; s_places = []
        r_centres = []; s_centres = []
        if labelled is not None and len(labelled):
            # Effective drawn radius per chiral atom (sphere/ball/stick-aware), so
            # the letter clears the atom's surface; floored at 0.5 A so it never
            # sits on top of a thin stick (and to keep clear of the Rama CA glyph).
            radii = labelled.chiral_atoms.display_radii(
                self._atomic_structure.ball_scale, 0.2)
            for i in range(len(labelled)):
                cfg = configs[i] if i < len(configs) else None
                if cfg is None:
                    continue
                clearance = max(float(radii[i]), 0.5)
                place = self._label_place(labelled[i], LABEL_HEIGHT_A, clearance)
                if cfg == 'R':
                    r_places.append(place); r_centres.append(labelled[i])
                else:
                    s_places.append(place); s_centres.append(labelled[i])
        self._r_centres = r_centres
        self._s_centres = s_centres
        self._set_label_drawing(self._r_drawing, r_places)
        self._set_label_drawing(self._s_drawing, s_places)
        self._apply_label_colors()

    @staticmethod
    def _set_label_drawing(d, places):
        if places:
            d.positions = Places(places)
            d.display = True
        else:
            d.display = False

    def set_label_color(self, color):
        '''Set the R/S label colour. ``color`` is the string ``'auto'`` (a tone
        contrasting with the background), ``'fromatoms'`` (each letter takes its
        chiral atom's colour), or a uint8 RGBA sequence (a fixed custom colour).'''
        if isinstance(color, str):
            self._label_color_mode = color.lower()
        else:
            self._label_color_mode = numpy.array(color, numpy.uint8)
        self._apply_label_colors()

    def _resolve_uniform_label_color(self):
        mode = self._label_color_mode
        if not isinstance(mode, str):
            return mode  # stored custom RGBA
        # 'auto': a light or dark tone, whichever contrasts with the background.
        bg = self.session.main_view.background_color
        lum = 0.2126 * bg[0] + 0.7152 * bg[1] + 0.0722 * bg[2]
        self._last_bg = tuple(bg)
        return numpy.array((20, 20, 20, 255) if lum > 0.5 else (240, 240, 240, 255),
                           numpy.uint8)

    def _apply_label_colors(self):
        # Colours must be set AFTER positions (they size to len(positions)).
        mode = self._label_color_mode
        if isinstance(mode, str) and mode == 'fromatoms':
            for d, centres in ((self._r_drawing, self._r_centres),
                               (self._s_drawing, self._s_centres)):
                if centres:
                    d.colors = numpy.array(
                        [cc.chiral_atom.color for cc in centres], numpy.uint8)
        else:
            col = self._resolve_uniform_label_color()
            if self._r_drawing.display:
                self._r_drawing.color = col
            if self._s_drawing.display:
                self._s_drawing.color = col

    def _rescale_for_camera(self, *_):
        # Persistent 'graphics update' handler, gated so an idle view does no work:
        #   * zoom (pixel size) changed  -> rebuild the view-width-scaled outlier
        #     glyphs (the fixed-size labels are unaffected by zoom);
        #   * view direction changed     -> re-face the camera-facing R/S letters;
        #   * background changed + 'auto' -> recompute the contrasting label colour.
        if not self.visible:
            return
        have_outliers = self._current_chirals is not None and len(self._current_chirals)
        have_labels = self._current_labelled is not None and len(self._current_labelled)
        if not have_outliers and not have_labels:
            return
        if have_outliers:
            ps = self.session.main_view.pixel_size()
            last = self._last_pixel_size
            if last is None or abs(ps - last) > 1e-4 * last:
                self._rebuild_outlier_places()
        if have_labels:
            vd = self._camera_view_direction()
            lvd = self._last_view_dir
            mode = self._label_color_mode
            if lvd is None or float(numpy.dot(vd, lvd)) < 0.99995:
                self._rebuild_label_places()
            elif isinstance(mode, str) and mode == 'auto' \
                    and tuple(self.session.main_view.background_color) != self._last_bg:
                self._apply_label_colors()

    def _camera_view_direction(self):
        # Unit direction (scene coords) the camera looks along; the rotation sentinel
        # for the camera-facing label layer.
        cam = self.session.main_view.camera.position
        return self._unit(numpy.asarray(cam.transform_vector((0, 0, -1)), dtype=float))

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

    def _make_letter_drawing(self, name, geometry):
        # A child drawing holding one precomputed extruded letter mesh (unit height,
        # centred at the origin), instanced across its centres via .positions.
        v, n, t = geometry
        d = Model('chiral %s label' % name, self.session)
        d.skip_bounds = True
        d.pickable = False
        d.set_geometry(v, n, t)
        # Actual colour is set per-rebuild by _apply_label_colors (mode-dependent).
        d.display = False
        return d

    def _label_place(self, cc, s, clearance):
        # Camera-facing placement for a unit-height letter -- a real depth-tested
        # billboard. The letter's local axes map to the camera's right/up/toward-
        # viewer directions, expressed in the STRUCTURE's local frame (the drawing
        # is a child of the structure, so its Places are in structure coordinates);
        # this keeps the letter readable from any angle while remaining occluded by
        # nearer geometry. Offset along screen-up so the letter's near edge sits
        # `clearance` from the atom centre (so it floats just outside the atom's
        # surface), with a small toward-viewer nudge to avoid coplanar fighting.
        C = numpy.asarray(cc.chiral_atom.coord, dtype=float)
        cam = self.session.main_view.camera.position
        sp_inv = self._atomic_structure.scene_position.inverse()
        right = self._unit(numpy.asarray(
            sp_inv.transform_vector(cam.transform_vector((1, 0, 0))), dtype=float))
        up = self._unit(numpy.asarray(
            sp_inv.transform_vector(cam.transform_vector((0, 1, 0))), dtype=float))
        out = self._unit(numpy.asarray(
            sp_inv.transform_vector(cam.transform_vector((0, 0, 1))), dtype=float))
        pos = C + up * (clearance + 0.5 * s) + out * (0.15 * s)
        M = numpy.empty((3, 4))
        M[:, 0] = right * s
        M[:, 1] = up * s
        M[:, 2] = out * s
        M[:, 3] = pos
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
        # The opt-in R/S `label` flag is live state on each ChiralCenter, but the
        # centres themselves are transient (rebuilt on demand by the manager, not
        # serialized). So the annotator -- which IS session-saved -- persists which
        # centres are labelled by their chiral atoms (Atoms serialize natively) and
        # re-applies the flag on restore.
        from chimerax.isolde.molobject import get_chiral_mgr
        from chimerax.atomic import Atoms
        mgr = get_chiral_mgr(session)
        chirals = mgr.get_chirals(self._atomic_structure.atoms, create=False)
        if chirals is not None and len(chirals):
            labelled_atoms = chirals[chirals.labels].chiral_atoms
        else:
            labelled_atoms = Atoms()
        data = {
            'model state': Model.take_snapshot(self, session, flags),
            'structure': self._atomic_structure,
            'labelled atoms': labelled_atoms,
            'label color mode': self._label_color_mode,
        }
        from .. import ISOLDE_STATE_VERSION
        data['version'] = ISOLDE_STATE_VERSION
        return data

    @staticmethod
    def restore_snapshot(session, data):
        ca = ChiralAnnotator(data['structure'])
        Model.set_state_from_snapshot(ca, session, data['model state'])
        mode = data.get('label color mode', None)
        if mode is not None:
            ca._label_color_mode = mode
        la = data.get('labelled atoms', None)
        if la is not None and len(la):
            from chimerax.isolde.molobject import get_chiral_mgr
            mgr = get_chiral_mgr(session)
            chirals = mgr.get_chirals(la, create=True)
            if len(chirals):
                chirals.labels = True
                ca.update_graphics()
        return ca
