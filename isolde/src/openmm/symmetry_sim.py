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
# A residue clashes with its own crystallographic image only under an operator
# that leaves some point fixed: a rotation about an axis, or -- in the
# centrosymmetric space groups of small-molecule crystals -- an inversion
# centre, a mirror plane or a rotoinversion. A screw or glide carries every atom
# a fraction of a cell repeat along its axis or plane, and a translation
# obviously cannot self-overlap, so none of those can produce a self-clash. We
# therefore (1) test whether an operator instance has a fixed point at all --
# skipping screws/glides/translations outright -- and only then (2) test the
# residue's actual image against itself.
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
# |screw/glide translation| (Angstroms) below which an operator counts as having
# a fixed point (self-overlap capable). Real screw and glide translations are
# fractions of a cell repeat -- many Angstroms -- so this is a wide, safe
# separation.
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
    # Vectorised over all copy<->copy pairs: a small-molecule shell holds hundreds
    # of operators (a ~10 A cell searched cutoff-deep), and the pairwise class
    # count is O(n^4) when done one comparison at a time.
    ops = numpy.array([numpy.asarray(o, dtype=numpy.float64) for o in ops_by_group])
    n = len(ops)
    w = numpy.zeros((n, n))
    w[0, 0] = 1.0
    w[0, 1:] = w[1:, 0] = 0.5
    # w[x][x], x>0, stays 0.0
    if n > 2:
        xs, ys = numpy.triu_indices(n - 1, k=1)
        xs += 1
        ys += 1
        R, t = ops[:, :, :3], ops[:, :, 3]
        Rinv = numpy.transpose(R, (0, 2, 1))
        # S = O_x^-1 . O_y and its inverse, as _op_compose(_op_inverse(.), .) does.
        SR = Rinv[xs] @ R[ys]
        St = (Rinv[xs] @ t[ys][:, :, None])[:, :, 0] - (Rinv[xs] @ t[xs][:, :, None])[:, :, 0]
        SRinv = numpy.transpose(SR, (0, 2, 1))
        Sinvt = -(SRinv @ St[:, :, None])[:, :, 0]
        S = numpy.concatenate([SR, St[:, :, None]], axis=2).reshape(-1, 12)
        Sinv = numpy.concatenate([SRinv, Sinvt[:, :, None]], axis=2).reshape(-1, 12)
        copies = ops[1:].reshape(-1, 12)

        def present(D, chunk=1024):
            out = numpy.empty(len(D), bool)
            for i in range(0, len(D), chunk):
                d = numpy.abs(D[i:i+chunk, None, :] - copies[None, :, :]) <= tol
                out[i:i+chunk] = d.all(axis=2).any(axis=1)
            return out

        # Relationship classes {S, S^-1}: assign each pair to the first earlier
        # representative it matches. Distinct operators differ by far more than
        # tol, so this is the same partition as comparing every pair to every other.
        reps = numpy.empty_like(S)
        n_reps = 0
        cls = numpy.empty(len(S), int)
        for p in range(len(S)):
            if n_reps:
                r = reps[:n_reps]
                hit = numpy.flatnonzero(
                    (numpy.abs(r - S[p]) <= tol).all(axis=1)
                    | (numpy.abs(r - Sinv[p]) <= tol).all(axis=1))
                if len(hit):
                    cls[p] = hit[0]
                    continue
            reps[n_reps] = S[p]
            cls[p] = n_reps
            n_reps += 1
        m = numpy.bincount(cls)[cls]
        r = 0.5 * present(S) + 0.5 * present(Sinv)
        val = numpy.maximum(0.0, (1.0 - r) / m)
        w[xs, ys] = w[ys, xs] = val

    # Flatten to match symmetry_group_table / Discrete2DFunction(n, n, ...):
    # values[i + n*j] = weight(group i, group j). The table is symmetric.
    return w.T.ravel().tolist()


def _operator_has_fixed_point(symmat):
    '''
    Return ``True`` if the orthogonal-space operator ``symmat`` (a 3x4
    ``[R | t]`` array, Angstroms) leaves some point in space fixed -- a rotation
    about an axis, a mirror, an inversion centre or a rotoinversion -- and so
    can map a nearby residue onto itself. ``False`` for translations, screw
    rotations and glides, which carry every atom at least their intrinsic
    translation away, so they cannot self-overlap.

    The intrinsic (screw/glide) translation is the part of ``t`` lying in the
    subspace ``R`` leaves fixed: the axis of a rotation, the plane of a mirror,
    nothing at all for an inversion or rotoinversion (which therefore always
    have a fixed point). ``x -> Rx + t`` has a fixed point exactly when it
    vanishes. It is evaluated for *this* operator instance, so a lattice
    translation along an axis or plane reads as a large screw or glide and is
    skipped, as it should be.
    '''
    m = numpy.asarray(symmat, dtype=numpy.float64)
    r, t = m[:, :3], m[:, 3]
    # The fixed subspace of R is the null space of (I - R). Every other singular
    # value of (I - R) is >= 1 for a crystallographic operator (6-fold minimum).
    _, s, vt = numpy.linalg.svd(numpy.identity(3) - r)
    fixed = vt[s < 1e-3]
    intrinsic = fixed.T @ (fixed @ t)
    return float(numpy.linalg.norm(intrinsic)) < SCREW_TRANSLATION_TOLERANCE


def _classify_self_overlap(residue_heavy_coords, image_heavy_coords,
        clash_distance=SELF_OVERLAP_CLASH_DISTANCE):
    '''
    Decide whether a residue's symmetry image overlaps the residue itself. Call
    only for operators with a fixed point (see :func:`_operator_has_fixed_point`);
    screws, glides and translations cannot self-overlap.

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


def _shell_fragments(residues, mobile_atoms):
    '''
    Label the atoms of ``residues`` by the fragment their symmetry copies are
    completed to. Returns ``(atoms, labels)``, with ``labels[i]`` the fragment of
    ``atoms[i]``.

    A residue in ``mobile_atoms`` is split into its connected pieces, using only
    bonds inside the residue. A protein, nucleic-acid or ligand residue is one
    piece, so this is the same as whole-residue completion. A small-molecule
    (corecif) model puts the whole asymmetric unit in one residue, which splits
    into its separate molecules: completing every copy to the whole asymmetric
    unit instead brought 124 atoms per operator for 2226487 (4 x 31-atom
    molecules), 35% more copies than needed.

    A residue outside ``mobile_atoms`` stays whole: if one of its atoms gets a
    copy, the whole residue is promoted to mobile, and every promoted atom then
    needs an image to carry its MDFF term.
    '''
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    atoms = residues.atoms
    n = len(atoms)
    if not n:
        return atoms, numpy.empty(0, dtype=int)
    res_idx = residues.indices(atoms.residues)
    a1, a2 = atoms.intra_bonds.atoms
    i1, i2 = atoms.indices(a1), atoms.indices(a2)
    same = res_idx[i1] == res_idx[i2]
    graph = coo_matrix((numpy.ones(int(same.sum())), (i1[same], i2[same])),
        shape=(n, n))
    n_frag, labels = connected_components(graph, directed=False)
    in_mobile = mobile_atoms.indices(atoms) >= 0
    labels = numpy.where(in_mobile, labels, n_frag + res_idx)
    return atoms, labels


def _complete_fragments(found_atoms, sym_indices, frag_atoms, frag_labels):
    '''
    Expand per-atom symmetry hits (``found_atoms`` under ``sym_indices``) to
    whole fragments (see :func:`_shell_fragments`), per operator. Returns
    ``(atoms, sym_indices, fragment_labels)`` grouped by ascending operator.
    '''
    if not len(found_atoms):
        return found_atoms, sym_indices, numpy.empty(0, dtype=int)
    fi = frag_atoms.indices(found_atoms)
    parts, ops = [], []
    for k in numpy.unique(sym_indices):
        hit = numpy.unique(frag_labels[fi[sym_indices == k]])
        sel = numpy.flatnonzero(numpy.isin(frag_labels, hit))
        parts.append(sel)
        ops.append(numpy.full(len(sel), k, dtype=sym_indices.dtype))
    idx = numpy.concatenate(parts)
    return frag_atoms[idx], numpy.concatenate(ops), frag_labels[idx]


# Two atoms closer than this (Angstroms) are the same atom twice. A molecule completed
# across a special position re-creates its own symmetry images to coordinate precision,
# while distinct atoms never get this close: bonded heavy atoms are >= ~1.2 A apart and
# genuine contacts (even metal coordination, ~1.9 A) are further still.
COINCIDENT_IMAGE_DISTANCE = 0.5


def _duplicate_fragment_images(real_atoms, frag_atoms, frag_labels, symmats, keys,
        tol=COINCIDENT_IMAGE_DISTANCE):
    '''
    Return the subset of ``keys`` -- ``(fragment, operator)`` images, in shell order --
    whose heavy atoms all land on atoms already present: real atoms, or an image kept
    earlier in ``keys``.

    This happens when the model holds a whole molecule that sits on a special position
    (e.g. after ``clipper fragments complete``): with site symmetry ``S``, ``g.M`` and
    ``g.S.M`` are the same neighbouring molecule, so every neighbour would be built twice,
    atom on top of atom. The coincident pairs belong to different operator groups, whose
    relationship ``S`` is absent from the shell (its self-image was pruned), so the group
    weight table gives them a nonzero copy<->copy weight and they interact at r ~ 0 (COD
    2201140, P2_1/n, on an inversion centre: every one of 5616 copies had a twin).
    '''
    from scipy.spatial import cKDTree
    heavy = real_atoms[real_atoms.element_names != 'H']
    pts = [(heavy if len(heavy) else real_atoms).coords]
    owner = [numpy.full(len(pts[0]), -1)]
    for k, (frag, op) in enumerate(keys):
        a = frag_atoms[frag_labels == frag]
        h = a[a.element_names != 'H']
        ref = h if len(h) else a
        pts.append(ref.coords @ symmats[op][:, :3].T + symmats[op][:, 3])
        owner.append(numpy.full(len(ref), k))
    pts = numpy.concatenate(pts)
    owner = numpy.concatenate(owner)
    pairs = cKDTree(pts).query_pairs(tol, output_type='ndarray')
    pairs = pairs[owner[pairs[:, 0]] != owner[pairs[:, 1]]]
    if not len(pairs):
        return set()
    partners = {}
    for a, b in pairs:
        partners.setdefault(a, []).append(b)
        partners.setdefault(b, []).append(a)
    kept = {-1}
    drop = set()
    for k, key in enumerate(keys):
        idx = numpy.flatnonzero(owner == k)
        if all(any(owner[q] in kept for q in partners.get(p, ())) for p in idx):
            drop.add(key)
        else:
            kept.add(k)
    return drop


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

    Copies are completed to whole fragments (see :func:`_shell_fragments`):
    whole residues, except that a mobile residue holding several molecules is
    split into them. Fragments whose image under a given operator overlaps
    themselves (special positions) are pruned per :func:`_classify_self_overlap`.

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
        cs.structure, cs.cell, cs.grid, cs.unit_cell, mobile_atoms, cutoff,
        whole_residues=False)
    frag_atoms, frag_labels = _shell_fragments(found_atoms.unique_residues,
        mobile_atoms)
    found_atoms, sym_indices, found_frags = _complete_fragments(found_atoms,
        sym_indices, frag_atoms, frag_labels)

    symmats = numpy.asarray(symmats, dtype=numpy.float64)

    # Which operators have a fixed point (capable of self-overlap)? Screws,
    # glides and translations cannot map a residue onto itself, so we never even
    # test them. Classify each operator instance once.
    has_fixed_point = {}
    for symop_index in set(int(s) for s in sym_indices):
        if symop_index == 0:
            continue
        has_fixed_point[symop_index] = _operator_has_fixed_point(symmats[symop_index])

    # Cache the self-overlap verdict per (fragment, operator) so we only classify
    # each combination once, then emit per-atom copies. A fragment is the unit
    # copies are completed to (see _shell_fragments): a residue, or one molecule
    # of a residue that holds several, so one molecule on a special position
    # never prunes the others' images.
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

    copies, copy_keys = [], []
    for atom, symop_index, frag in zip(found_atoms, sym_indices, found_frags):
        symop_index = int(symop_index)
        # Index 0 is the identity (the ASU itself): the "copy" would be the real
        # atom, so skip it.
        if symop_index == 0:
            continue
        residue = atom.residue
        if id(residue) in ignored_residue_ids:
            n_ignored_copies += 1
            continue
        key = (int(frag), symop_index)
        verdict = overlap_verdict.get(key)
        if verdict is None:
            if not has_fixed_point[symop_index]:
                # Screw / glide / translation: cannot self-overlap.
                verdict = 'none'
            else:
                r = symmats[symop_index][:, :3]
                t = symmats[symop_index][:, 3]
                res_atoms = frag_atoms[frag_labels == frag]
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
        copy_keys.append(key)

    # Drop images that only re-create atoms already present (a molecule on a special
    # position held whole in the model builds each neighbour twice).
    unique_keys = list(dict.fromkeys(copy_keys))
    duplicates = _duplicate_fragment_images(cs.structure.atoms, frag_atoms, frag_labels,
        symmats, unique_keys)
    if duplicates:
        copies = [c for c, k in zip(copies, copy_keys) if k not in duplicates]

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
        if duplicates:
            msg += (f' Dropped {len(duplicates)} duplicate residue/operator image(s): '
                'the model holds a whole molecule that sits on a special position, '
                'which re-creates each neighbouring molecule twice.')
        logger.info(msg)

    if not copies:
        return None
    return SymmetryShell(model, cs, symmats, copies)
