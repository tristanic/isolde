# @Author: Tristan Croll
# @Date:   12-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 12-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Sourcing layer for crystallographic symmetry-equivalent atoms in ISOLDE
simulations.

ChimeraX-Clipper already draws symmetry "ghost" atoms around the asymmetric
unit, but the simulation has historically been blind to them. Recent OpenMM
(>= 8.4) provides :class:`openmm.SymmetrySite` -- a virtual site that ties a
particle's position to a single parent via a rigid-body transform
``r' = R.r + v`` and folds any force on the copy back onto the parent as
``R^T.f``. That primitive lets us represent each symmetry mate as a zero-DOF
copy of its parent atom, so crystal contacts become real, two-way forces during
interactive MD.

This module is purely the *sourcing* layer: it decides which symmetry copies to
build and hands back a lightweight, OpenMM-free description of them. Actually
instantiating the particles/virtual sites and wiring them into the forces is the
job of :class:`chimerax.isolde.openmm.openmm_interface.SimHandler`.

The crystallographic operators come from Clipper via
:func:`chimerax.clipper.sym_realize.crystal_symmetry_for` and
:func:`~chimerax.clipper.sym_realize.sym_select_within`, which return, for a set
of query atoms, every master-model atom whose symmetry image falls within a
cutoff, tagged with the operator (a 3x4 orthogonal-space matrix, in Angstroms)
that maps it there.
'''

import numpy

# Residue-level self-overlap detection.
#
# A residue clashes with its own crystallographic image only under a *proper
# rotation* operator (in the Sohncke space groups of biomolecular crystals, the
# only operators are rotations, screw rotations and translations; a screw
# carries an on-axis atom tens of Angstroms along the axis, and a translation
# obviously cannot self-overlap, so neither can produce a self-clash). We
# therefore (1) use Clipper's affine axis/screw decomposition to test whether an
# operator instance is a pure rotation at all -- skipping screws/translations
# outright -- and only then (2) test the residue's actual image against itself.
#
# The overlap test is a direct heavy-atom clash check (image atom within
# SELF_OVERLAP_CLASH_DISTANCE of *any* original heavy atom). Measuring real
# image-vs-original distances catches every overlap mode -- an atom approaching
# its own image, and two atoms *swapping* across the axis (bonded or not) -- with
# no proxy and no fold-dependent threshold. Two distinct heavy atoms never sit
# this close in a genuine crystal contact (heavy-heavy bottoms out ~2.5 A, even
# for a short H-bond), so the threshold cleanly separates self-overlap from a
# real (even tight) self-symmetry contact. (e.g. GOL /A:308 in 9cwy, P 6 2 2:
# a 2-fold threads the glycerol body, image atoms land 0.5-1.5 A from originals.)
SELF_OVERLAP_CLASH_DISTANCE = 2.0
# |screw translation| (Angstroms) below which an operator counts as a pure
# rotation (self-overlap capable). Real screw translations are fractions of a
# cell repeat -- many Angstroms -- so this is a wide, safe separation.
SCREW_TRANSLATION_TOLERANCE = 0.5


def symmetry_sims_available():
    '''
    True if the installed OpenMM exposes the :class:`SymmetrySite` virtual site
    (OpenMM >= 8.4). Gate the whole feature on this **capability check**, never
    on a parsed version number: the reference build reports its version as
    ``8.4.0.dev-...``, which orders *below* ``8.4.0`` under PEP 440 and would be
    wrongly rejected by a ``>= 8.4.0`` comparison.
    '''
    import openmm
    return hasattr(openmm, 'SymmetrySite')


class SymmetryCopy:
    '''
    Description of a single symmetry-equivalent atom to be added to the
    simulation. Deliberately OpenMM-free: it names the parent :class:`Atom` and
    the index of the operator (into :attr:`SymmetryShell.symmats`) that maps the
    parent to this copy. Coordinates, particle indices and virtual sites are
    resolved later by the SimHandler.
    '''
    __slots__ = ('parent_atom', 'symop_index', 'source_altloc')

    def __init__(self, parent_atom, symop_index, source_altloc):
        self.parent_atom = parent_atom
        self.symop_index = symop_index
        # The alternate-conformer id of the parent at build time. Currently
        # always '' (ISOLDE strips altlocs before simulating), but carried
        # through so the altloc-vs-symmetry-altloc case can be handled later
        # without re-plumbing (see module notes / project plan).
        self.source_altloc = source_altloc


class SymmetryShell:
    '''
    The full set of symmetry copies to instantiate for one simulation, plus the
    crystallographic context needed to place them.

    Attributes:
        * model: the master :class:`AtomicStructure`.
        * crystal_symmetry: the Clipper
          :class:`~chimerax.clipper.sym_realize.CrystalSymmetry`.
        * symmats: an ``(N, 3, 4)`` float64 array of orthogonal-space operators
          in Angstroms (row-major ``[R | t]``); index 0 is always the identity.
        * copies: a list of :class:`SymmetryCopy` (identity and pruned
          self-overlaps already removed).
    '''
    def __init__(self, model, crystal_symmetry, symmats, copies):
        self.model = model
        self.crystal_symmetry = crystal_symmetry
        self.symmats = numpy.asarray(symmats, dtype=numpy.float64)
        self.copies = copies

    def __len__(self):
        return len(self.copies)

    @property
    def parent_atoms(self):
        '''The unique parent :class:`Atoms` referenced by the copies, in a
        stable order.'''
        from chimerax.atomic import Atoms
        seen = {}
        for c in self.copies:
            seen.setdefault(id(c.parent_atom), c.parent_atom)
        return Atoms(list(seen.values()))

    def operator(self, symop_index):
        '''Return ``(R, t)`` for an operator: R a ``(3, 3)`` rotation and t a
        ``(3,)`` translation in Angstroms.'''
        m = self.symmats[symop_index]
        return m[:, :3], m[:, 3]


# ---------------------------------------------------------------------------
# Operator algebra on orthogonal-space 3x4 [R|t] matrices (Angstroms).
#
# The orthogonal matrices bake in the specific lattice translation, so equality
# within tolerance is EXACT operator identity (two operators differing by a whole
# cell are tens of Angstroms apart in t and never compare equal). That lets the
# symmetry group-weight table (below) test operator relationships without any
# fractional/lattice bookkeeping.
# ---------------------------------------------------------------------------

def _op_inverse(m):
    '''Inverse of an orthogonal 3x4 operator (R orthonormal): ``[R^T | -R^T t]``.'''
    R = m[:, :3]
    t = m[:, 3]
    Rinv = R.T
    return numpy.concatenate([Rinv, (-Rinv @ t)[:, None]], axis=1)


def _op_compose(a, b):
    '''Compose two orthogonal 3x4 operators: ``a . b`` (apply b, then a).'''
    Ra, ta = a[:, :3], a[:, 3]
    Rb, tb = b[:, :3], b[:, 3]
    return numpy.concatenate([Ra @ Rb, (Ra @ tb + ta)[:, None]], axis=1)


def _op_equal(a, b, tol=1e-3):
    '''True if two orthogonal 3x4 operators are equal within ``tol`` (Angstroms).'''
    return numpy.allclose(a, b, atol=tol, rtol=0.0)


def symmetry_group_weight_table(ops_by_group, tol=1e-3):
    '''
    Build the flattened row-major ``n*n`` weight table for the symmetry
    ``grouptable`` (an :class:`openmm.Discrete2DFunction`), given the orthogonal
    operator of each group.

    ``ops_by_group`` is a list of ``(3, 4)`` orthogonal ``[R | t]`` matrices
    indexed by group id: index 0 is the identity (the real asymmetric unit);
    ``1..M`` are one per crystallographic operator actually present in the sim.

    Weights make each *unique* crystal contact contribute exactly once (all
    representations of a contact fold the identical force back to the real atoms,
    so the only requirement is that the weights of all present representations of
    a relationship class ``{S, S^-1}`` sum to 1):

        * ``(0, 0) = 1``            -- real <-> real (ASU internal nonbonded).
        * ``(X, X)``, X>0 ``= 0``   -- a copy's internal energy duplicates the
          ASU's; excluded.
        * real <-> copy ``= 1/2``   -- the two mirror legs ``{0,S}``/``{0,S^-1}``
          of an inter-molecular contact, each half-weighted (the established
          double-count fix).
        * copy <-> copy ``{X,Y}`` (relationship ``S = O_X^-1 . O_Y``) ``= the
          residual`` ``(1 - 1/2*present(S) - 1/2*present(S^-1)) / m``, where
          ``present(D)`` is whether operator ``D`` is one of the copy groups and
          ``m`` is how many copy<->copy pairs share the class ``{S, S^-1}``. When
          both mirror legs are present (whole-construct sims) this is 0 -- exactly
          the previous behaviour. When they are absent (a local simulation at a
          multi-way interface whose linking operator ``S`` was never added to the
          shell) it becomes ``1/m``, so the otherwise-dropped copy<->copy contact
          is counted exactly once.
    '''
    ops = [numpy.asarray(o, dtype=numpy.float64) for o in ops_by_group]
    n = len(ops)

    def present(D):
        return any(_op_equal(ops[g], D, tol) for g in range(1, n))

    copy_pairs = [(x, y) for x in range(1, n) for y in range(x + 1, n)]
    pair_S = {p: _op_compose(_op_inverse(ops[p[0]]), ops[p[1]]) for p in copy_pairs}

    def same_class(s1, s2):
        return _op_equal(s1, s2, tol) or _op_equal(s1, _op_inverse(s2), tol)

    w = [[0.0] * n for _ in range(n)]
    w[0][0] = 1.0
    for x in range(1, n):
        w[0][x] = w[x][0] = 0.5
        # w[x][x] stays 0.0
    for p in copy_pairs:
        S = pair_S[p]
        m = sum(1 for q in copy_pairs if same_class(S, pair_S[q]))
        r = 0.5 * present(S) + 0.5 * present(_op_inverse(S))
        val = max(0.0, (1.0 - r) / m)
        x, y = p
        w[x][y] = w[y][x] = val

    # Flatten to match symmetry_group_table / Discrete2DFunction(n, n, ...):
    # values[i + n*j] = weight(group i, group j). The table is symmetric.
    return [w[i][j] for j in range(n) for i in range(n)]


def _operator_is_pure_rotation(symmat):
    '''
    Return ``True`` if the orthogonal-space operator ``symmat`` (a 3x4
    ``[R | t]`` array, Angstroms) is a proper rotation about a fixed axis --
    i.e. capable of mapping a nearby residue onto itself. ``False`` for pure
    translations and for screw rotations (whose along-axis translation carries
    an on-axis atom many Angstroms away, so they cannot self-overlap).

    Uses Clipper's affine axis/screw decomposition
    (:func:`chimerax.clipper.geometry.util.rotation_screw_axis_and_angle_affine`),
    which returns the along-axis screw component of *this* operator instance --
    correctly folding in any lattice translation (a rotation plus a whole-cell
    shift along its axis reads as a large screw and is skipped, as it should be).
    '''
    from chimerax.geometry import Place
    from chimerax.clipper.geometry.util import rotation_screw_axis_and_angle_affine
    place = Place(matrix=numpy.asarray(symmat, dtype=numpy.float64))
    const_point, axis_dir, angle, screw = rotation_screw_axis_and_angle_affine(place)
    if const_point is None:
        return False        # pure translation: no rotation axis
    screw_mag = float(numpy.linalg.norm(screw))
    return screw_mag < SCREW_TRANSLATION_TOLERANCE


def _classify_self_overlap(residue_heavy_coords, image_heavy_coords,
        clash_distance=SELF_OVERLAP_CLASH_DISTANCE):
    '''
    Decide whether a residue's symmetry image (under a pure-rotation operator)
    overlaps the residue itself. Call only for pure-rotation operators (see
    :func:`_operator_is_pure_rotation`); screws/translations cannot self-overlap.

    Returns one of:
        * ``'none'``             -- no heavy image atom comes within
          ``clash_distance`` of any original heavy atom (a normal, possibly-
          contacting, symmetry mate: keep it).
        * ``'special_position'`` -- the image clashes into the residue (some
          heavy image atom lands within ``clash_distance`` of an original heavy
          atom). This is a direct all-pairs test, so it catches an atom
          approaching its own image *and* two atoms swapping across the axis
          (bonded or not) alike, and needs no fold-dependent threshold.

    A future ``'cross_altloc'`` verdict (kept, occupancy-weighted) will be added
    once ISOLDE can simulate alternate conformers; that is the sole reason this
    is a single decision seam rather than an inline test.
    '''
    # residue_heavy_coords: (M, 3); image_heavy_coords: (M, 3). Minimum distance
    # between any image heavy atom and any original heavy atom.
    deltas = image_heavy_coords[:, None, :] - residue_heavy_coords[None, :, :]
    min_dist = numpy.sqrt((deltas * deltas).sum(axis=2)).min()
    if min_dist < clash_distance:
        return 'special_position'
    return 'none'


def build_symmetry_shell(model, mobile_atoms, cutoff, *, logger=None):
    '''
    Work out which crystallographic symmetry copies should participate in a
    simulation of ``mobile_atoms``, and return them as a :class:`SymmetryShell`
    (or ``None`` if the feature does not apply to this model).

    Args:
        * model: the master :class:`AtomicStructure` being simulated.
        * mobile_atoms: the mobile selection; copies are sourced for every
          symmetry image approaching within ``cutoff`` of these.
        * cutoff: search radius in Angstroms. The caller chooses this per the
          two-regime policy (existing fixed-shell distance for local sims;
          ``max(nonbonded, gbsa) + margin`` when the whole model is mobile).
        * logger: optional ChimeraX logger for status messages.

    Residues whose image under a given operator overlaps themselves (special
    positions) are pruned per :func:`_classify_self_overlap`.

    Returns ``None`` (with a logged note) when OpenMM lacks ``SymmetrySite`` or
    the model has no real crystallographic symmetry (a plain P1 box).
    '''
    if not symmetry_sims_available():
        if logger is not None:
            logger.info('ISOLDE: symmetry-aware simulation requested but the '
                'installed OpenMM has no SymmetrySite support (requires '
                '>= 8.4); continuing without symmetry atoms.')
        return None

    from chimerax.clipper.sym_realize import (crystal_symmetry_for,
        sym_select_within)
    cs = crystal_symmetry_for(model)
    if not cs.has_symmetry:
        if logger is not None:
            logger.info('ISOLDE: model has no crystallographic symmetry (P1 / '
                'non-crystallographic); no symmetry atoms added.')
        return None

    found_atoms, symmats, sym_indices, symops = sym_select_within(
        cs.structure, cs.cell, cs.grid, cs.unit_cell, mobile_atoms, cutoff)

    symmats = numpy.asarray(symmats, dtype=numpy.float64)

    # Which operators are proper rotations (capable of self-overlap)? Screws and
    # translations cannot map a residue onto itself, so we never even test them.
    # Classify each operator instance once.
    pure_rotation = {}
    for symop_index in set(int(s) for s in sym_indices):
        if symop_index == 0:
            continue
        pure_rotation[symop_index] = _operator_is_pure_rotation(symmats[symop_index])

    # Cache the residue-level self-overlap verdict per (residue, operator) so we
    # only classify each combination once, then emit per-atom copies.
    overlap_verdict = {}
    n_pruned_residues = 0

    # Residues flagged ``isolde_ignore`` (e.g. via "isolde ignore :UNL") are
    # excluded from sim generation entirely: they carry no MD parameters, so
    # neither the residue nor its crystallographic image can be simulated. We
    # must not source copies for such parents - otherwise the parent would be
    # promoted to the mobile selection and then subtracted straight back out with
    # the excluded residues, leaving the shell referencing an atom that is absent
    # from the simulation construct (which the parent-lookup guard in
    # ``SimHandler.initialize_symmetry_copies`` then reports as a hard error).
    ignored_residue_ids = {id(r) for r in found_atoms.unique_residues
        if getattr(r, 'isolde_ignore', False)}
    n_ignored_copies = 0

    copies = []
    for atom, symop_index in zip(found_atoms, sym_indices):
        symop_index = int(symop_index)
        # Index 0 is the identity (the ASU itself): the "copy" would be the real
        # atom, so skip it.
        if symop_index == 0:
            continue
        residue = atom.residue
        if id(residue) in ignored_residue_ids:
            n_ignored_copies += 1
            continue
        key = (id(residue), symop_index)
        verdict = overlap_verdict.get(key)
        if verdict is None:
            if not pure_rotation[symop_index]:
                # Screw / translation: cannot self-overlap.
                verdict = 'none'
            else:
                r = symmats[symop_index][:, :3]
                t = symmats[symop_index][:, 3]
                res_atoms = residue.atoms
                heavy = res_atoms[res_atoms.element_names != 'H']
                # Fall back to all atoms if the residue is all-hydrogen (never
                # happens for real residues, but keep the array non-empty).
                ref = heavy if len(heavy) else res_atoms
                res_coords = ref.coords
                image_coords = res_coords @ r.T + t
                verdict = _classify_self_overlap(res_coords, image_coords)
            overlap_verdict[key] = verdict
            if verdict != 'none':
                n_pruned_residues += 1
        # In v1 (altlocs stripped before simulation) the only non-'none' verdict
        # is 'special_position', which we prune. The single-seam design lets the
        # future 'cross_altloc' verdict be kept here instead.
        if verdict == 'special_position':
            continue
        # Normalise the "no alternate conformer" marker (ChimeraX may report it
        # as '' or ' ') to '' so the future cross-altloc comparison is clean.
        source_altloc = (getattr(atom, 'alt_loc', '') or '').strip()
        copies.append(SymmetryCopy(atom, symop_index, source_altloc))

    if logger is not None:
        n_parents = len({id(c.parent_atom) for c in copies})
        msg = (f'ISOLDE: sourced {len(copies)} crystallographic symmetry '
            f'copies from {n_parents} parent atoms (spacegroup '
            f'{cs.spacegroup.symbol_hm}).')
        if n_pruned_residues:
            msg += (f' Pruned {n_pruned_residues} residue/operator '
                'combination(s) coincident with a special position.')
        if n_ignored_copies:
            msg += (f' Skipped {n_ignored_copies} copy(ies) whose parent lies in '
                'an isolde-ignored residue.')
        logger.info(msg)

    if not copies:
        return None
    return SymmetryShell(model, cs, symmats, copies)
