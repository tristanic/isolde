# @Author: Tristan Croll
# @Date:   12-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 12-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Symmetry ghost markup — a first-class capability for markup managers
====================================================================

Each ISOLDE validation/restraint manager (RamaAnnotator, RotamerAnnotator,
PositionRestraintMgr, TuggableAtomsMgr, ...) draws its markup on the real atoms of
a model. When the model is crystallographic, Clipper draws "ghost" copies of the
atoms around the asymmetric unit. This module lets a manager draw its *own* markup
on those ghosts, so a symmetry mate looks like a real atom.

The manager owns a :class:`SymmetryGhostDrawing` (composition, not inheritance,
so it works uniformly across the annotators — which subclass ``Model`` directly —
and the restraint managers — which subclass ``_RestraintMgr``). The helper:

* creates one ghost child :class:`Drawing` per glyph spec, copying the glyph
  geometry from the manager's own real drawing so they match exactly;
* subscribes to Clipper's ``'symmetry display changed'`` trigger (fired once per
  ghost redraw) so it refreshes in lockstep with the ghosts, without polling;
* on refresh, asks Clipper which ghosts are drawn
  (``AtomicSymmetryModel.drawn_ghosts_by_operator()``) and, for each operator,
  asks the manager for its markup on exactly those atoms, transforms it by the
  operator, and draws it (dimmed to match Clipper's ghost grey-out).

The **only** visibility criterion is "is the ghost drawn?" — never the parent
atom's visibility (in spotlight mode Clipper hides distant real atoms while
drawing their nearby ghost, and a ghost can be in a simulation's visible region
while its parent is not). The helper hands the manager only atoms whose ghosts are
drawn; the manager filters solely on its own enabled/validity state.

The per-manager producer (:meth:`symmetry_markup_instances`) is the single seam:
today it returns numpy arrays computed from the manager's existing C++-backed
collection properties, and the ``operator ∘ transform`` composition below is a
cheap ``numpy.einsum``. If the future "parallelise annotator updates via
std::async" work lands, the C++-native managers (position/tuggable/rama) can move
that composition behind a C++ entry point without changing this interface.
'''

import numpy
from chimerax.core.models import Drawing


def _compose_transforms(place, tf):
    '''
    Compose a :class:`Place` operator with a stack of ``(n,3,4)`` instance
    transforms, returning ``operator . transform`` as ``(n,3,4)`` float64
    (``Places(place_array=...)`` requires float64).
    '''
    m = place.matrix
    R = m[:, :3]
    t = m[:, 3]
    tf = numpy.asarray(tf)
    out = numpy.empty((len(tf), 3, 4), numpy.float64)
    out[:, :, :3] = numpy.einsum('ij,njk->nik', R, tf[:, :, :3])
    out[:, :, 3] = tf[:, :, 3] @ R.T + t
    return out


def _as_place_array(value):
    '''
    Normalise a per-instance transform source to an ``(n,3,4)`` array. ISOLDE
    managers expose these as a raw OpenGL ``(n,4,4)`` array (position restraints,
    tuggables) or as an already-built :class:`Places` (distance restraints).
    '''
    from chimerax.geometry import Places
    if isinstance(value, Places):
        return value.array()
    arr = numpy.asarray(value)
    if arr.ndim == 3 and arr.shape[1:] == (4, 4):
        return Places(opengl_array=arr).array()
    return arr


class GhostSpec:
    '''
    Describes one ghost glyph for a manager.

    * ``name``: identifier passed back to ``manager.symmetry_markup_instances``.
    * ``source_drawing``: the manager's own real :class:`Drawing` for this glyph;
      its geometry is copied so the ghost matches, and (for uniform-coloured
      glyphs) its ``.color`` is used.
    * ``kind``: ``'points'`` (``shift_and_scale`` spheres/pins), ``'transforms'``
      (``(n,3,4)`` cylinder/arrow instance transforms), or ``'mesh'`` (bespoke
      per-frame geometry, e.g. the cis/twisted omega "cups" — the manager returns
      ``(vertices, normals, triangles, colors)`` and the operator transforms the
      vertices/normals directly).
    * ``radius_fn``: for ``'points'`` glyphs, a zero-arg callable returning the
      current sphere radius (default 1.0).
    '''
    def __init__(self, name, source_drawing, kind, radius_fn=None):
        self.name = name
        self.source_drawing = source_drawing
        self.kind = kind
        self.radius_fn = radius_fn if radius_fn is not None else (lambda: 1.0)
        self.ghost_drawing = None


class SymmetryGhostDrawing:
    '''
    Draws a manager's markup on the currently-displayed Clipper symmetry ghosts.

    Args:
        * manager: the owning markup manager (a ChimeraX ``Model``). Must
          implement ``symmetry_markup_instances(spec_name, atoms) -> (geometry,
          colors) | None`` (see module docstring).
        * structure: the :class:`AtomicStructure` the manager belongs to.
        * specs: a list of :class:`GhostSpec`.
    '''
    def __init__(self, manager, structure, specs):
        self.manager = manager
        self.structure = structure
        self.session = manager.session
        self.specs = specs
        for spec in specs:
            sd = spec.source_drawing
            gd = Drawing('{} (symmetry)'.format(sd.name))
            gd.pickable = False
            gd.skip_bounds = True
            # 'mesh' glyphs regenerate their geometry every refresh, so there is
            # nothing (and possibly no valid geometry yet) to copy here.
            if spec.kind != 'mesh':
                gd.set_geometry(sd.vertices, sd.normals, sd.triangles)
            gd.display = False
            manager.add_drawing(gd)
            spec.ghost_drawing = gd
        self._sym_handler = None
        self._trigger_handler = None
        self._subscribe()

    # --- Clipper subscription --------------------------------------------
    def _subscribe(self):
        if self._trigger_handler is not None:
            return
        try:
            from chimerax.clipper.symmetry import get_symmetry_handler
            sh = get_symmetry_handler(self.structure, create=False)
        except Exception:
            sh = None
        if sh is None:
            return
        self._sym_handler = sh
        self._trigger_handler = sh.triggers.add_handler(
            'symmetry display changed', self._display_changed_cb)

    def _display_changed_cb(self, _name, asm):
        self.refresh(asm)

    def _get_asm(self):
        sh = self._sym_handler
        if sh is None or sh.deleted:
            # Handler absent or recreated (e.g. after session restore) — (re)resolve.
            self._sym_handler = None
            self._trigger_handler = None
            self._subscribe()
            sh = self._sym_handler
            if sh is None:
                return None
        try:
            return sh.atomic_symmetry_model
        except Exception:
            return None

    # --- rebuild ----------------------------------------------------------
    def _hide(self):
        for spec in self.specs:
            spec.ghost_drawing.display = False

    def refresh(self, asm=None):
        if asm is None:
            asm = self._get_asm()
        if asm is None or not asm.visible:
            self._hide()
            return
        try:
            ops = asm.drawn_ghosts_by_operator()
        except Exception:
            self._hide()
            return
        if not ops:
            self._hide()
            return
        try:
            dim = float(asm.dim_factor)
        except Exception:
            dim = 1.0
        for spec in self.specs:
            # Guard per spec so one misbehaving glyph can't break the others or the
            # manager's real drawing (refresh() runs from update_graphics).
            try:
                if spec.kind == 'mesh':
                    self._refresh_mesh(spec, ops, dim)
                else:
                    self._refresh_instanced(spec, ops, dim)
            except Exception as e:
                spec.ghost_drawing.display = False
                if not getattr(spec, '_reported_error', False):
                    spec._reported_error = True
                    self.session.logger.warning(
                        'Symmetry markup "{}" disabled after error: {}'.format(spec.name, e))

    def _refresh_instanced(self, spec, ops, dim):
        '''Points (shift_and_scale) or transforms (place_array) glyphs.'''
        from chimerax.geometry import Places
        mgr = self.manager
        d = spec.ghost_drawing
        accum = []
        accum_colors = []
        have_colors = False
        for place, atoms in ops:
            res = mgr.symmetry_markup_instances(spec.name, atoms)
            if res is None:
                continue
            geom, colors = res
            if geom is None or not len(geom):
                continue
            if spec.kind == 'points':
                gc = place.transform_points(numpy.asarray(geom))
                xyzr = numpy.ones((len(gc), 4), numpy.float32)
                xyzr[:, :3] = gc
                xyzr[:, 3] = spec.radius_fn()
                accum.append(xyzr)
            else:
                accum.append(_compose_transforms(place, _as_place_array(geom)))
            if colors is not None:
                accum_colors.append(colors)
                have_colors = True
        if not accum:
            d.display = False
            return
        merged = numpy.concatenate(accum)
        if spec.kind == 'points':
            d.positions = Places(shift_and_scale=merged.astype(numpy.float32))
        else:
            d.positions = Places(place_array=merged.astype(numpy.float64))
        if have_colors:
            colors = numpy.concatenate(accum_colors).astype(numpy.float32)
            # Match Clipper's ghost dimming: scale RGB (not alpha).
            colors[:, :3] *= dim
            d.colors = colors.astype(numpy.uint8)
        else:
            uc = spec.source_drawing.color
            if uc is not None:
                uc = numpy.array(uc, numpy.float32)
                uc[:3] *= dim
                d.color = uc.astype(numpy.uint8)
        d.display = True

    def _refresh_mesh(self, spec, ops, dim):
        '''
        Bespoke per-frame geometry (e.g. cis/twisted omega cups). The manager
        returns (vertices, normals, triangles, colors) for the given atoms; the
        operator transforms the vertices (points) and normals (vectors), and the
        per-operator meshes are concatenated (with triangle-index offsets) into the
        single ghost drawing.
        '''
        mgr = self.manager
        d = spec.ghost_drawing
        vs = []
        ns = []
        ts = []
        cs = []
        offset = 0
        have_colors = True
        for place, atoms in ops:
            res = mgr.symmetry_markup_instances(spec.name, atoms)
            if res is None:
                continue
            v, n, t, c = res
            if v is None or not len(v):
                continue
            vs.append(place.transform_points(numpy.asarray(v, numpy.float32)))
            ns.append(place.transform_vectors(numpy.asarray(n, numpy.float32)))
            ts.append(numpy.asarray(t, numpy.int32) + offset)
            offset += len(v)
            if c is None:
                have_colors = False
            else:
                cs.append(numpy.asarray(c))
        if not vs:
            d.display = False
            return
        d.set_geometry(numpy.concatenate(vs).astype(numpy.float32),
                       numpy.concatenate(ns).astype(numpy.float32),
                       numpy.concatenate(ts).astype(numpy.int32))
        if have_colors and cs:
            C = numpy.concatenate(cs).astype(numpy.float32)
            C[:, :3] *= dim
            d.vertex_colors = C.astype(numpy.uint8)
        d.display = True

    def delete(self):
        h = self._trigger_handler
        if h is not None and self._sym_handler is not None and not self._sym_handler.deleted:
            self._sym_handler.triggers.remove_handler(h)
        self._trigger_handler = None
        self._sym_handler = None
