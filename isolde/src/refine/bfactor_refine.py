# @Author: Tristan Croll <tic20>
# @Date:   10-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 13-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll

'''
Front-end ChimeraX commands (``isolde brefine`` / ``isolde brsr``) driving the
maximum-likelihood isotropic B-factor / occupancy refinement engine that lives
in the sister ChimeraX-Clipper package.  The engine itself (and all
crystallographic machinery) stays in Clipper; these commands only call into its
existing public API.

The engine no longer builds its own restraints: the caller (i.e. this module)
constructs them and passes them to ``launch()`` / ``launch_realspace()`` as
``b_restraints`` (two-sided pairwise) and ``b_target_restraints`` (one-sided)
tuples, expressed as integer indices into the supplied atom array plus explicit
altloc strings.  Restraint construction here is both altloc-aware and
symmetry-aware (see the builder helpers below).
'''

from math import pi as _PI

from chimerax.core.commands import AtomSpecArg

#: Convert an isotropic B-factor (Å²) to U_iso (Å²): U = B / (8·π²).
U_FROM_B = 1.0 / (8.0 * _PI * _PI)

#: Per-type Barron quadratic-core scale ``c`` (U-space, Å²; ≡ the old GM ``σ``).
#: This sets the quadratic-core width / origin curvature (= weight/c²) and is a
#: pure stiffness knob — the tail shape is now controlled independently by ``α``
#: (see below).  Bonded atoms share B closely → tight c (c_B ≈ 5 Å²); non-bonded
#: context contacts can differ more → wider c (c_B ≈ 15 Å²).
_SIGMA_BONDED = 0.063    # c_B ≈ 5 Å²
_SIGMA_CONTEXT = 0.19    # c_B ≈ 15 Å²

#: Distance kernel scale (Å): w(d) = ε²/(d²+ε²) — ≈1 at contact, decays outward.
_DIST_KERNEL_EPS = 3.5

#: Barron robust-loss shape parameter ``α`` (CVPR 2019).  α=2 harmonic, α=1
#: Charbonnier (quadratic core + saturating linear tail — the default), α=0
#: Cauchy, α=−2 reproduces the old Geman-McClure family, α→−∞ Welsch.  Lower α ⇒
#: more saturating/redescending tail ⇒ a pair is allowed to diverge if the data
#: insists.  ``alpha auto`` (per-restraint) maps a [0,1] confidence onto this band:
#: high confidence (rigid bonds, atoms in solid VdW contact) → ALPHA_HIGH; low
#: confidence (rotamer tips, loose non-contacts) → ALPHA_LOW.
ALPHA_DEFAULT = 1.0      # Charbonnier — flat-α default when no value is given
ALPHA_HIGH = 1.0         # confidence = 1 → Charbonnier (keep B-factors tracking)
ALPHA_LOW = -2.0         # confidence = 0 → Geman-McClure (release floppy pairs)

#: Superatom contact-confidence tuning (non-bonded families, ``alpha auto``).  The
#: refined set is heavy-only, but the model carries H (explicit or implicit), so a
#: heavy-atom VdW radius is inflated by H_SHELL to recover the H-group envelope
#: (hydrophobic CHₙ packing); a donor/acceptor pair gets HBOND_BONUS of extra
#: reach (polar / H-mediated contact, inferred from idatm type — no explicit H
#: needed).  Confidence ramps from 1 at solid overlap to 0 once a GAP_FULL gap
#: opens.  Å throughout; first guesses — refine from the log.
_H_SHELL = 0.8
_GAP_FULL = 1.0
_HBOND_BONUS = 0.5

#: Master normalization constants ("fraction of data energy a restraint family may
#: spend"; see build_b_restraints_for_refinement).  These bridge the unitless data
#: energy to the engine's restraint weight scale and are intended to be CALIBRATED
#: ONCE from a logged run: pick η so a user weight of ~1.0 gives gentle-but-active
#: regularization (reproduces the empirically-good raw ≈ 0.001 at 3 Å).  Because
#: the per-run data energy S already tracks dataset size/|Fo|², the calibrated η is
#: then portable across resolutions.
#: Calibrated for Barron's normalized loss with α auto (the runaway-atom artifact
#: is now controlled by the tail shape, not by cranking weights): a user weight of
#: ~1.0 gives well-behaved regularization on a 3 Å test case.  These are 5× the
#: original GM-tuned values, which were a touch too weak under the new loss.
_ETA_XTAL = 5.0e-7       # crystallographic (brefine); Barron/α-auto calibrated, w≈1
_ETA_RS = 50.0           # real-space (brsr); Barron/α-auto calibrated, w≈1


class _StructuresOrSymmetryMgrsArg(AtomSpecArg):
    name = "a models specifier"

    @classmethod
    def parse(cls, text, session):
        from chimerax.clipper.symmetry import SymmetryManager
        from chimerax.atomic import AtomicStructure
        aspec, text, rest = super().parse(text, session)
        models = aspec.evaluate(session).models
        sym_mgrs = [m for m in models if isinstance(m, SymmetryManager)]
        structures = [m for m in models if type(m) == AtomicStructure
                      and m.parent not in sym_mgrs]
        return structures + sym_mgrs, text, rest


class _MapHandlerArg(AtomSpecArg):
    name = "a Clipper map handler"

    @classmethod
    def parse(cls, text, session):
        from chimerax.clipper.maps.map_handler_base import MapHandlerBase
        aspec, text, rest = super().parse(text, session)
        models = aspec.evaluate(session).models
        handlers = [m for m in models if isinstance(m, MapHandlerBase)]
        return handlers, text, rest


# -----------------------------------------------------------------------------
# Restraint construction
# -----------------------------------------------------------------------------

def _normalize_altloc(al):
    '''Map ChimeraX's blank altloc (``' '`` or ``None``) to the engine's ``''``.'''
    if al is None:
        return ''
    return al.strip()


def _distance_weight(d):
    '''Smooth inverse-square distance kernel; works on scalars or numpy arrays.'''
    e2 = _DIST_KERNEL_EPS * _DIST_KERNEL_EPS
    return e2 / (d * d + e2)


# -----------------------------------------------------------------------------
# Per-restraint Barron α ("auto-confidence") heuristics
#
# In ``alpha auto`` mode each restraint's tail shape is set from a [0,1] local
# confidence: rigid bonds / atoms in solid VdW contact → ALPHA_HIGH (commit, keep
# B-factors tracking); rotamer tips / loose non-contacts → ALPHA_LOW (let the data
# pull them apart).  Classification is driven entirely by ChimeraX ``idatm_type``
# + ``bond.in_cycle`` + heavy-atom degree, so it works for ligands / nucleotides /
# modified residues, not just the standard amino acids.
# -----------------------------------------------------------------------------

def _alpha_from_confidence(c):
    '''Map a confidence in [0,1] onto the α band (ALPHA_LOW … ALPHA_HIGH).'''
    if c < 0.0:
        c = 0.0
    elif c > 1.0:
        c = 1.0
    return ALPHA_LOW + c * (ALPHA_HIGH - ALPHA_LOW)


def _heavy_degree(atom):
    '''Number of non-hydrogen atoms covalently bonded to ``atom``.'''
    return sum(1 for n in atom.neighbors if n.element.number > 1)


def _is_sp3(idatm_type, type_info, tetrahedral):
    '''True if the IDATM type is tetrahedral (sp3) — i.e. rotatable single-bond
    character at this atom.'''
    info = type_info.get(idatm_type)
    return info is not None and info.geometry == tetrahedral


def _is_donor(atom, type_info):
    '''H-bond donor test from heavy-atom IDATM type (no explicit H needed); mirrors
    the predicate in chimerax.clashes — an N/O/S that either has room for an
    (implicit) H or already carries one.'''
    if atom.element.number not in (7, 8, 16):   # N, O, S
        return False
    info = type_info.get(atom.idatm_type)
    if info is not None and atom.num_bonds < info.substituents:
        return True
    return any(n.element.number == 1 for n in atom.neighbors)


def _is_acceptor(atom, type_info):
    '''H-bond acceptor test from IDATM type (lone pair available): substituents <
    geometry.  Mirrors chimerax.clashes.'''
    info = type_info.get(atom.idatm_type)
    return info is not None and info.substituents < info.geometry


def _effective_radius(atom, type_info):
    '''Superatom VdW radius: the heavy-atom radius inflated by ``_H_SHELL`` when the
    atom bears hydrogens (explicit or implicit), recovering the H-group envelope so
    a heavy-only contact test still sees hydrophobic CHₙ packing.'''
    r = atom.default_radius
    info = type_info.get(atom.idatm_type)
    n_h = (info.substituents - _heavy_degree(atom)) if info is not None \
        else sum(1 for n in atom.neighbors if n.element.number == 1)
    if n_h > 0:
        r += _H_SHELL
    return r


def _atom_contact_props(atoms):
    '''Precompute (superatom radius, donor flag, acceptor flag) for an ``Atoms``
    array, returned as numpy arrays parallel to ``atoms`` (so callers can index by
    the same position they already hold).  Only needed for ``alpha auto``.'''
    import numpy
    from chimerax.atomic.idatm import type_info
    n = len(atoms)
    radii = numpy.empty(n, dtype=numpy.float64)
    donor = numpy.zeros(n, dtype=bool)
    acceptor = numpy.zeros(n, dtype=bool)
    for i in range(n):
        a = atoms[i]
        radii[i] = _effective_radius(a, type_info)
        donor[i] = _is_donor(a, type_info)
        acceptor[i] = _is_acceptor(a, type_info)
    return radii, donor, acceptor


def _overlap_confidence(r_i, r_j, donor_i, acc_i, donor_j, acc_j, dist):
    '''Contact confidence in [0,1] from superatom-radius overlap, with a donor/
    acceptor reach bonus for polar / H-mediated contacts.  1 at solid overlap,
    ramping to 0 once a ``_GAP_FULL`` gap opens.'''
    reach = r_i + r_j
    if (donor_i and acc_j) or (donor_j and acc_i):
        reach += _HBOND_BONUS
    c = 1.0 + (reach - dist) / _GAP_FULL
    if c < 0.0:
        return 0.0
    return 1.0 if c > 1.0 else c


def _bonded_alpha(bond, a1, a2, alpha, type_info, tetrahedral):
    '''Per-bond α for the bonded family.  Flat passthrough unless ``alpha == 'auto'``,
    in which case a 3-tier confidence: rigid (ring / non-sp3 endpoint, or a rigid
    coordination hub) → high; acyclic sp3–sp3 interior → mid; terminal tip (an
    endpoint with one heavy neighbour, e.g. Lys Cε–Nζ / Ser Cβ–Oγ) → low.'''
    if alpha != 'auto':
        return float(alpha)
    if bond.in_cycle or not (_is_sp3(a1.idatm_type, type_info, tetrahedral)
                             and _is_sp3(a2.idatm_type, type_info, tetrahedral)):
        return _alpha_from_confidence(1.0)
    d1, d2 = _heavy_degree(a1), _heavy_degree(a2)
    if d1 == 1 or d2 == 1:
        # Terminal heavy atom.  Distinguish a rigid coordination hub from a genuine
        # floppy tip by the *parent* (the non-terminal endpoint; both are already
        # tetrahedral here).  A non-carbon centre holding ≥3 heavy ligands — a
        # sulfate S (Sac), phosphate/phosphonate P (Pac), sulfonyl S, … — is a rigid
        # polyhedron whose central atom moves *with* the group; treat it as rigid
        # (persistent Charbonnier) so the GM tail can't abandon a runaway centre
        # (e.g. sulfate S drifting to b_max while its oxygens stay low).  Charbonnier
        # still saturates, so a weakly-bound group can go out-of-step if the data
        # insists.  Carbon-centred tips (hydroxyl, amine, thiol, methyl, branch
        # methyl) keep the floppy GM tier.
        parent, parent_deg = (a2, d2) if d1 == 1 else (a1, d1)
        if parent.element.number != 6 and parent_deg >= 3:
            return _alpha_from_confidence(1.0)
        return _alpha_from_confidence(0.0)
    return _alpha_from_confidence(0.5)


class _ConformerTable:
    '''Flat per-(atom, altloc) table.  One row per conformer.

    ``refined_index`` is the position in the atom array the table was built from
    (for the refined table that is exactly the index space the engine expects).
    Coordinates are whatever was stored (for the context table they are already
    symmetry-transformed).
    '''
    def __init__(self, refined_index, altloc, coord, bfactor, element_number):
        import numpy
        self.refined_index = numpy.asarray(refined_index, dtype=numpy.int32)
        self.altloc = numpy.asarray(list(altloc), dtype='<U4')
        if len(coord):
            self.coord = numpy.asarray(coord, dtype=numpy.float64).reshape((-1, 3))
        else:
            self.coord = numpy.zeros((0, 3), dtype=numpy.float64)
        self.bfactor = numpy.asarray(bfactor, dtype=numpy.float64)
        self.element_number = numpy.asarray(element_number, dtype=numpy.int32)
        self.master_index = None   # set only for context tables (dedup key)

    def __len__(self):
        return len(self.refined_index)


def expand_conformers(atoms):
    '''Expand a ChimeraX ``Atoms`` array into a :class:`_ConformerTable`.

    Single-conformer atoms (the common case) are read vectorized; multi-altloc
    atoms are read one conformer at a time by switching ``alt_loc`` inside
    ``suppress_alt_loc_change_notifications()`` (auto-restores).
    '''
    import numpy
    n = len(atoms)
    elnums = atoms.elements.numbers
    coords = atoms.coords
    bfactors = atoms.bfactors
    nalt = atoms.num_alt_locs
    multi = nalt > 0
    single_idx = numpy.nonzero(~multi)[0]

    refined_index = [int(i) for i in single_idx]
    altloc = [''] * len(single_idx)
    coord = [coords[i] for i in single_idx]
    bfactor = [bfactors[i] for i in single_idx]
    elnum = [int(elnums[i]) for i in single_idx]

    for i in numpy.nonzero(multi)[0]:
        a = atoms[int(i)]
        for al in a.alt_locs:
            with a.suppress_alt_loc_change_notifications():
                a.alt_loc = al
                c = a.coord
                b = a.bfactor
            refined_index.append(int(i))
            altloc.append(_normalize_altloc(al))
            coord.append(c)
            bfactor.append(b)
            elnum.append(int(elnums[i]))
    return _ConformerTable(refined_index, altloc, coord, bfactor, elnum)


def _altlocs_by_atom(table):
    from collections import defaultdict
    d = defaultdict(list)
    for row in range(len(table)):
        d[int(table.refined_index[row])].append(str(table.altloc[row]))
    return d


def _pair_altlocs(als_i, als_j):
    '''Label-match + blank-wildcard altloc pairing for a bonded/contact pair.

    ``''`` (single-conformer) acts as a wildcard partner.  Matching labels pair;
    a label with no partner falls back to the other atom's blank conformer if it
    has one; genuinely orphaned labels are dropped.
    '''
    si, sj = set(als_i), set(als_j)
    if si == {''} and sj == {''}:
        return [('', '')]
    has_blank_i, has_blank_j = ('' in si), ('' in sj)
    labels_i, labels_j = si - {''}, sj - {''}
    pairs = []
    matched_i, matched_j = set(), set()
    for lab in labels_i & labels_j:
        pairs.append((lab, lab))
        matched_i.add(lab)
        matched_j.add(lab)
    if has_blank_j:
        for lab in labels_i - matched_i:
            pairs.append((lab, ''))
    if has_blank_i:
        for lab in labels_j - matched_j:
            pairs.append(('', lab))
    return pairs


def _pair_key(i, ai, j, aj):
    a, b = (i, ai), (j, aj)
    return (a, b) if a <= b else (b, a)


def _bonded_index_pairs(refined_atoms):
    '''Set of ``frozenset({i, j})`` (positions in ``refined_atoms``) for every
    covalent bond fully inside the set.'''
    bonds = refined_atoms.intra_bonds
    if len(bonds) == 0:
        return set()
    a1s, a2s = bonds.atoms
    i1 = refined_atoms.indices(a1s)
    i2 = refined_atoms.indices(a2s)
    out = set()
    for k in range(len(bonds)):
        if i1[k] >= 0 and i2[k] >= 0:
            out.add(frozenset((int(i1[k]), int(i2[k]))))
    return out


def build_bonded_b_restraints(refined_atoms, table, ignore_hydrogens,
                              alpha=ALPHA_DEFAULT):
    '''Two-sided pairwise restraints along covalent bonds (heavy atoms only).

    The ``ks`` column holds **relative** weights (ω = 1 for every bond); the
    orchestrator rescales them to the data-normalized family budget.  ``alpha`` is
    either a float (broadcast) or ``'auto'`` (per-bond 3-tier confidence; see
    :func:`_bonded_alpha`).  Returns a
    ``(atoms1, altlocs1, atoms2, altlocs2, sigmas, ks, alphas)`` tuple, or ``None``.
    '''
    bonds = refined_atoms.intra_bonds
    if len(bonds) == 0:
        return None
    a1s, a2s = bonds.atoms
    idx1 = refined_atoms.indices(a1s)
    idx2 = refined_atoms.indices(a2s)
    el1 = a1s.elements.numbers
    el2 = a2s.elements.numbers
    als_by_atom = _altlocs_by_atom(table)
    auto = (alpha == 'auto')
    if auto:
        from chimerax.atomic.idatm import type_info, tetrahedral
    flat_alpha = None if auto else float(alpha)
    atoms1, altlocs1, atoms2, altlocs2, sigmas, ks, alphas = \
        [], [], [], [], [], [], []
    seen = set()
    for b in range(len(bonds)):
        i, j = int(idx1[b]), int(idx2[b])
        if i < 0 or j < 0:
            continue
        if el1[b] == 1 or el2[b] == 1:   # heavy-only
            continue
        # α depends only on bond geometry, so it is shared across altloc pairs.
        b_alpha = (_bonded_alpha(bonds[b], a1s[b], a2s[b], alpha, type_info,
                                 tetrahedral) if auto else flat_alpha)
        for ai, aj in _pair_altlocs(als_by_atom.get(i, ['']),
                                    als_by_atom.get(j, [''])):
            key = _pair_key(i, ai, j, aj)
            if key in seen:
                continue
            seen.add(key)
            atoms1.append(i)
            altlocs1.append(ai)
            atoms2.append(j)
            altlocs2.append(aj)
            sigmas.append(_SIGMA_BONDED)
            ks.append(1.0)   # relative weight ω; orchestrator rescales
            alphas.append(b_alpha)
    if not atoms1:
        return None
    return (atoms1, altlocs1, atoms2, altlocs2, sigmas, ks, alphas)


def _apply_place(place, coord):
    import numpy
    return place.transform_points(numpy.array([coord], dtype=numpy.float64))[0]


def _close_sets(core_coords, search_coords, symmats, cutoff):
    '''Symmetry-aware narrowing via ``find_close_points_sets``.

    Returns a list (one int array per operator) of indices into ``search_coords``
    whose image under that operator lies within ``cutoff`` of ``core_coords``.
    This is the optimised C++ proximity primitive; it returns *sets* (which
    points have some nearby counterpart), so exact pairing is done afterward on
    the narrowed candidates.
    '''
    import numpy
    from chimerax.geometry import Place, find_close_points_sets
    core = numpy.asarray(core_coords, dtype=numpy.float32)
    search = numpy.asarray(search_coords, dtype=numpy.float32)
    idmat = numpy.asarray(Place().matrix, dtype=numpy.float32)
    tp_target = [(core, idmat)]
    tp_search = [(search, numpy.asarray(m, dtype=numpy.float32)) for m in symmats]
    i1, _i2 = find_close_points_sets(tp_search, tp_target, cutoff)
    return i1


def _expand_context_candidates(context_master, symmats, per_op):
    '''Build a context :class:`_ConformerTable` (symmetry-transformed coords) from
    the per-operator candidate lists returned by :func:`_close_sets`.

    ``per_op[k]`` are indices into ``context_master`` whose image under operator
    ``symmats[k]`` lies near the refined selection.  Each (atom, altloc) becomes a
    row with its transformed coordinate; ``master_index`` (the position in
    ``context_master``) lets the caller dedup a context atom seen under several
    operators to its nearest image.
    '''
    import numpy
    from chimerax.geometry import Place
    places = [Place(matrix=numpy.asarray(m)) for m in symmats]
    alt, crd, bf, el, mi = [], [], [], [], []
    for k, cand in enumerate(per_op):
        if len(cand) == 0:
            continue
        op = places[k]
        for c in cand:
            c = int(c)
            a = context_master[c]
            als = a.alt_locs
            if len(als) == 0:
                crd.append(_apply_place(op, a.coord))
                bf.append(a.bfactor)
                el.append(a.element.number)
                alt.append('')
                mi.append(c)
            else:
                for al in als:
                    with a.suppress_alt_loc_change_notifications():
                        a.alt_loc = al
                        cc = a.coord
                        bb = a.bfactor
                    crd.append(_apply_place(op, cc))
                    bf.append(bb)
                    el.append(a.element.number)
                    alt.append(_normalize_altloc(al))
                    mi.append(c)
    if not crd:
        return None
    table = _ConformerTable([-1] * len(crd), alt, crd, bf, el)
    table.master_index = numpy.asarray(mi, dtype=numpy.int32)
    return table


def build_boundary_b_target_restraints(sym_mgr, refined_atoms, refined_table,
                                       ignore_hydrogens, cutoff,
                                       alpha=ALPHA_DEFAULT):
    '''One-sided restraints pulling each refined atom toward the distance-weighted
    mean B of its (symmetry-aware) context neighbours.  The ``ks`` column holds
    relative weights (ω = 1 per target); the orchestrator rescales them.  ``alpha``
    is a float (broadcast) or ``'auto'`` — in auto mode each target's confidence is
    the *best* superatom-contact confidence over its contributing context
    neighbours (is this atom genuinely packed against context, or loosely
    surrounded?).

    Returns a ``(atoms, altlocs, target_us, sigmas, ks, alphas)`` tuple, or ``None``.
    '''
    if cutoff <= 0:
        return None
    import numpy
    from chimerax.geometry import find_close_points

    try:
        struct = refined_atoms.unique_structures[0]
    except Exception:
        return None
    # Candidate context = structure heavy atoms NOT in the refined set.  Because
    # refined atoms are excluded here, their symmetry images can never become
    # fixed targets (correct: a refined atom's B is still being optimised).
    context_master = struct.atoms.subtract(refined_atoms)
    context_master = context_master[context_master.elements.numbers != 1]
    if len(context_master) == 0:
        return None

    r_heavy = numpy.nonzero(refined_table.element_number != 1)[0]
    if len(r_heavy) == 0:
        return None

    # Symmetry-aware narrowing: which context atoms (under which operator) lie
    # near the refined selection — one optimised find_close_points_sets call.
    symmats = _operators_for(sym_mgr, refined_atoms, cutoff)
    per_op = _close_sets(refined_table.coord[r_heavy], context_master.coords,
                         symmats, cutoff)
    ctx = _expand_context_candidates(context_master, symmats, per_op)
    if ctx is None or len(ctx) == 0:
        return None

    c_coords = ctx.coord
    c_u = ctx.bfactor * U_FROM_B
    c_mi = ctx.master_index
    c_alt = ctx.altloc

    auto = (alpha == 'auto')
    flat_alpha = None if auto else float(alpha)
    if auto:
        r_radii, r_don, r_acc = _atom_contact_props(refined_atoms)
        c_radii, c_don, c_acc = _atom_contact_props(context_master)

    atoms, altlocs, target_us, alphas = [], [], [], []
    for rr in r_heavy:
        rc = refined_table.coord[rr]
        near = find_close_points(numpy.array([rc]), c_coords, cutoff)[1]
        if len(near) == 0:
            continue
        # Dedup each context conformer to its nearest symmetry image.
        best = {}
        for li in near:
            dvec = c_coords[li] - rc
            dist = float((dvec * dvec).sum()) ** 0.5
            key = (int(c_mi[li]), str(c_alt[li]))
            if key not in best or dist < best[key][0]:
                best[key] = (dist, li)
        dists = numpy.array([v[0] for v in best.values()])
        us = numpy.array([c_u[v[1]] for v in best.values()])
        w = _distance_weight(dists)
        ri = int(refined_table.refined_index[rr])
        atoms.append(ri)
        altlocs.append(str(refined_table.altloc[rr]))
        target_us.append(float((w * us).sum() / w.sum()))
        if auto:
            conf = 0.0
            for d, li in best.values():
                cm = int(c_mi[li])
                cc = _overlap_confidence(r_radii[ri], c_radii[cm], r_don[ri],
                                         r_acc[ri], c_don[cm], c_acc[cm], d)
                if cc > conf:
                    conf = cc
            alphas.append(_alpha_from_confidence(conf))
        else:
            alphas.append(flat_alpha)
    if not atoms:
        return None

    sigmas = [_SIGMA_CONTEXT] * len(atoms)
    ks = [1.0] * len(atoms)   # relative weight ω; orchestrator rescales
    return (atoms, altlocs, target_us, sigmas, ks, alphas)


def _operators_for(sym_mgr, atoms, cutoff):
    '''Symmetry operator matrices (3x4) relevant to ``atoms`` within ``cutoff``;
    always includes the identity.  Falls back to identity-only.'''
    import numpy
    from chimerax.geometry import Place
    idmat = numpy.asarray(Place().matrix)
    mats = None
    asm = getattr(sym_mgr, 'atomic_symmetry_model', None)
    if asm is not None:
        try:
            _found, symmats, _si, _so = asm.sym_select_within(
                atoms, cutoff, whole_residues=False)
            mats = [numpy.asarray(m) for m in symmats]
        except Exception:
            mats = None
    if not mats:
        return [idmat]
    if not any(numpy.allclose(m, idmat) for m in mats):
        mats = [idmat] + mats
    return mats


def build_local_context_b_restraints(sym_mgr, refined_atoms, table,
                                     ignore_hydrogens, cutoff, exclude_bonded=True,
                                     alpha=ALPHA_DEFAULT):
    '''Two-sided per-pair restraints between near (non-bonded) refined atoms,
    including crystallographic symmetry contacts.  The ``ks`` column holds
    relative weights ω = w(d)/√(K_i·K_j) (distance kernel ÷ the endpoints' contact
    counts); the orchestrator rescales them.  ``alpha`` is a float (broadcast) or
    ``'auto'`` (per-pair superatom-contact confidence; see
    :func:`_overlap_confidence`).  Returns a 7-tuple or ``None``.

    ``exclude_bonded`` — when the bonded family is also active, covalently-bonded
    pairs are handled there (with their tighter σ), so they are dropped here to
    avoid double restraints; the orchestrator passes ``False`` when bonded
    restraints are off, so those near contacts are still covered locally.

    Cost note: over a whole structure this enumerates all heavy contacts within
    ``cutoff`` (× symmetry operators) — bounded by the cutoff, but not free; it is
    on by default (``wLocal`` 1.0) because at ISOLDE's usual resolutions the local
    network is almost always needed.  Set ``wLocal 0`` to disable.
    '''
    if cutoff <= 0:
        return None
    import numpy
    from chimerax.geometry import Place, find_close_points

    heavy = numpy.nonzero(table.element_number != 1)[0]
    if len(heavy) < 2:
        return None
    H = table.coord[heavy]
    master = table.refined_index[heavy]
    altl = table.altloc[heavy]
    bonded = _bonded_index_pairs(refined_atoms) if exclude_bonded else set()

    symmats = _operators_for(sym_mgr, refined_atoms, cutoff)
    places = [Place(matrix=numpy.asarray(m)) for m in symmats]
    # Symmetry-aware narrowing: per operator, which atoms' images contact the
    # (untransformed) structure — one optimised find_close_points_sets call.
    per_op = _close_sets(H, H, symmats, cutoff)
    best = {}   # (row_i, row_j) sorted -> nearest contact distance
    for k, cand in enumerate(per_op):
        if len(cand) == 0:
            continue
        tc = places[k].transform_points(H[cand])
        for ci in range(len(cand)):
            j = int(cand[ci])
            near = find_close_points(numpy.array([tc[ci]]), H, cutoff)[1]
            for i in near:
                i = int(i)
                mi, mj = int(master[i]), int(master[j])
                if mi == mj:
                    continue
                if frozenset((mi, mj)) in bonded:
                    continue
                dvec = H[i] - tc[ci]
                dist = float((dvec * dvec).sum()) ** 0.5
                key = (i, j) if i < j else (j, i)
                if key not in best or dist < best[key]:
                    best[key] = dist
    if not best:
        return None

    # Contact count per master atom → coordination normalization.
    from collections import defaultdict
    contacts = defaultdict(int)
    for (a, b) in best:
        contacts[int(master[a])] += 1
        contacts[int(master[b])] += 1

    auto = (alpha == 'auto')
    flat_alpha = None if auto else float(alpha)
    if auto:
        radii, donor, acceptor = _atom_contact_props(refined_atoms)
    atoms1, altlocs1, atoms2, altlocs2, sigmas, ks, alphas = \
        [], [], [], [], [], [], []
    for (a, b), dist in best.items():
        mi, mj = int(master[a]), int(master[b])
        omega = float(_distance_weight(dist)) / ((contacts[mi] * contacts[mj]) ** 0.5)
        atoms1.append(mi)
        altlocs1.append(str(altl[a]))
        atoms2.append(mj)
        altlocs2.append(str(altl[b]))
        sigmas.append(_SIGMA_CONTEXT)
        ks.append(omega)   # relative weight ω; orchestrator rescales
        if auto:
            conf = _overlap_confidence(radii[mi], radii[mj], donor[mi], acceptor[mi],
                                       donor[mj], acceptor[mj], dist)
            alphas.append(_alpha_from_confidence(conf))
        else:
            alphas.append(flat_alpha)
    return (atoms1, altlocs1, atoms2, altlocs2, sigmas, ks, alphas)


def _concat_restraint_tuples(t1, t2, n):
    if t1 is None and t2 is None:
        return None
    if t1 is None:
        return t2
    if t2 is None:
        return t1
    return tuple(list(t1[k]) + list(t2[k]) for k in range(n))


def _strip_self_pairs(b_restraints):
    '''Drop any pairwise restraint whose two endpoints are the same atom (same
    refined index).  Atoms frequently sit close to their own symmetry copies, so
    a contact search can otherwise produce a restraint with the same atom on both
    sides; this is the final guarantee that none survive (the local-context
    builder also skips them at build time).  Returns the filtered 6-tuple, or
    ``None`` if nothing remains.'''
    if b_restraints is None:
        return None
    a1, al1, a2, al2, sig, ks, alphas = b_restraints
    keep = [k for k in range(len(a1)) if a1[k] != a2[k]]
    if len(keep) == len(a1):
        return b_restraints
    if not keep:
        return None
    return tuple([seq[k] for k in keep]
                 for seq in (a1, al1, a2, al2, sig, ks, alphas))


def _data_energy_scale(data_source, crystallographic, n_refined_heavy=0):
    '''Estimate the data-residual energy scale ``S`` used to make restraint weights
    portable.  Crystallographic: ``S = ½·Σ(k·|Fc| − m·|Fo|)²`` over finite working
    reflections (``k`` the least-squares Fo/Fc scale), read from the live xmap
    manager's ``f_obs``/``f_calc``/``weights`` via the ``HKL_data.data`` numpy
    export.  Real-space: ``S = ½·Σ(ρt − ⟨ρt⟩)²`` over the target map voxels (signal-
    energy proxy, since ρc is unavailable caller-side).  Returns a positive float;
    falls back to ``1.0`` if the data can't be read (normalization degrades to an
    η-scaled constant).'''
    import numpy
    try:
        if crystallographic:
            xm = getattr(data_source, 'live_xmap_mgr', None)
            if xm is None:
                return 1.0
            fo = numpy.asarray(xm.f_obs.data[1])[:, 0]
            fc = numpy.asarray(xm.f_calc.data[1])[:, 0]
            fom = numpy.asarray(xm.weights.data[1])[:, 1]
            mask = numpy.isfinite(fo) & numpy.isfinite(fc) & numpy.isfinite(fom)
            fo, fc, fom = fo[mask], fc[mask], fom[mask]
            if fo.size == 0:
                return 1.0
            mfo = fom * fo
            denom = float(numpy.sum(fc * fc))
            k = float(numpy.sum(mfo * fc) / denom) if denom > 0 else 1.0
            resid = k * fc - mfo
            S = 0.5 * float(numpy.sum(resid * resid))
            return S if S > 0 else 1.0
        else:
            # Local, fragment-size-portable proxy: per-atom density-contrast
            # energy.  Scaling by the refined heavy-atom count makes S grow with
            # the fragment (matching the engine's cropped data term), so the
            # per-restraint weight (family budget / Sum(omega), Sum(omega) ~ atom
            # count) is independent of how large the refined selection is.  The
            # earlier whole-map sum made small fragments need huge weights.
            from chimerax.map.volume import mean_sd_rms
            _mean, sd, _rms = mean_sd_rms(data_source.data.matrix())
            S = 0.5 * float(n_refined_heavy) * float(sd) * float(sd)
            return S if S > 0 else 1.0
    except Exception:
        return 1.0


def _normalize_family(restraints, w_user, eta, S):
    '''Rescale a family's *relative* ω weights (the ``ks`` column — now second-to-
    last, since the trailing column holds per-restraint α) to the data-normalized
    budget: ``k_eff(r) = w_user·η·S·ω_r / Σω``.  So the family's total weight ≈
    ``w_user·η·S`` regardless of restraint count, and per-restraint strength tracks
    the data energy ``S`` (→ portable).  Returns the tuple, or ``None`` if the
    family is empty/disabled.'''
    if restraints is None or w_user is None or w_user <= 0:
        return None
    cols = list(restraints)
    omega = cols[-2]
    total = float(sum(omega))
    if total <= 0:
        return None
    scale = float(w_user) * float(eta) * float(S) / total
    cols[-2] = [float(o) * scale for o in omega]
    return tuple(cols)


def build_b_restraints_for_refinement(sym_mgr, refined_atoms, context_atoms,
                                      data_source, *,
                                      ignore_hydrogens, crystallographic,
                                      internal_weight=0.0,
                                      boundary_weight=0.0, boundary_cutoff=0.0,
                                      local_weight=0.0, local_cutoff=0.0,
                                      alpha=ALPHA_DEFAULT):
    '''Build the restraint tuples for a refinement run, with per-family data
    normalization (see :func:`_normalize_family`).

    ``alpha`` is the Barron robust-loss shape — a float broadcast to every
    restraint, or ``'auto'`` for the per-restraint confidence heuristic.  It rides
    in the trailing column of each tuple (the engine defaults it to 1 if absent).

    Returns ``(b_restraints, b_target_restraints, diag)``; the first two may be
    ``None``; ``diag`` carries the data energy / budgets for calibration logging.

    - Crystallographic (``brefine``): ``b_restraints`` = bonded + symmetry-aware
      local-context (both two-sided pairwise); ``b_target_restraints`` is always
      ``None`` (the crystallographic ``launch()`` accepts only ``b_restraints``).
    - Real-space (``brsr``): ``b_restraints`` = bonded; ``b_target_restraints`` =
      one-sided distance-weighted-mean boundary restraints toward context.
    '''
    table = expand_conformers(refined_atoms)
    # Only pay for the data-energy scan if at least one family is active.
    any_active = bool(
        (internal_weight and internal_weight > 0)
        or (crystallographic and local_weight and local_weight > 0)
        or (not crystallographic and boundary_weight and boundary_weight > 0))
    n_refined_heavy = int((table.element_number != 1).sum())
    S = (_data_energy_scale(data_source, crystallographic, n_refined_heavy)
         if any_active else 1.0)
    eta = _ETA_XTAL if crystallographic else _ETA_RS
    diag = {'S_data': S, 'eta': eta, 'budgets': {}}

    bonded = None
    if internal_weight and internal_weight > 0:
        bonded = _normalize_family(
            build_bonded_b_restraints(refined_atoms, table, ignore_hydrogens,
                                      alpha=alpha),
            internal_weight, eta, S)
        if bonded is not None:
            diag['budgets']['bonded'] = internal_weight * eta * S

    if crystallographic:
        local = None
        if local_weight and local_weight > 0:
            # Bonded restraints (tighter σ) take precedence: only drop bonded
            # pairs from the local array when the bonded family is itself active.
            bonded_active = bool(internal_weight and internal_weight > 0)
            local = _normalize_family(
                build_local_context_b_restraints(
                    sym_mgr, refined_atoms, table, ignore_hydrogens, local_cutoff,
                    exclude_bonded=bonded_active, alpha=alpha),
                local_weight, eta, S)
            if local is not None:
                diag['budgets']['local'] = local_weight * eta * S
        b_restraints = _strip_self_pairs(_concat_restraint_tuples(bonded, local, 7))
        return b_restraints, None, diag

    boundary = None
    if boundary_weight and boundary_weight > 0:
        boundary = _normalize_family(
            build_boundary_b_target_restraints(
                sym_mgr, refined_atoms, table, ignore_hydrogens, boundary_cutoff,
                alpha=alpha),
            boundary_weight, eta, S)
        if boundary is not None:
            diag['budgets']['boundary'] = boundary_weight * eta * S
    return _strip_self_pairs(bonded), boundary, diag


def _log_restraint_counts(session, command, b_restraints, b_target_restraints,
                          diag, alpha):
    n_pair = len(b_restraints[0]) if b_restraints is not None else 0
    n_targ = len(b_target_restraints[0]) if b_target_restraints is not None else 0
    msg = f'{command}: {n_pair} pairwise B-factor restraint(s)'
    if b_target_restraints is not None:
        msg += f', {n_targ} one-sided target restraint(s)'

    # Weights are now the second-to-last column (the last column holds α).
    def _maxk(t):
        return max(t[-2]) if (t is not None and len(t[-2])) else 0.0
    rep = max(_maxk(b_restraints), _maxk(b_target_restraints))
    a_vals = []
    for t in (b_restraints, b_target_restraints):
        if t is not None and len(t[-1]):
            a_vals.extend(t[-1])
    a_desc = f'{min(a_vals):.3g}..{max(a_vals):.3g}' if a_vals else 'n/a'
    mode = 'auto' if alpha == 'auto' else f'{float(alpha):.3g}'
    budgets = ', '.join(f'{k}={v:.4g}' for k, v in diag['budgets'].items()) or 'none'
    msg += (f'.  [S_data={diag["S_data"]:.4g}, eta={diag["eta"]:.3g}; '
            f'family budgets: {budgets}; max k_eff={rep:.4g}; '
            f'alpha={mode} (range {a_desc})]')
    session.logger.info(msg)


# -----------------------------------------------------------------------------
# Commands
# -----------------------------------------------------------------------------

def isolde_brefine_cmd(session, structure=None,
                       refine_b=True, max_cycles=50, b_min_factor=2.0,
                       b_max=200.0, n_threads=None, rwork_tolerance=0.02,
                       log_level='info', auto_tolerances=True,
                       delta_tolerance_factor=1.0, epsilon_tolerance_factor=1.0,
                       ignore_hydrogens=True,
                       w_internal=1.0,
                       w_local=1.0,
                       c_local=4.0,
                       alpha=ALPHA_DEFAULT):
    from chimerax.core.errors import UserError
    from chimerax.clipper.symmetry import SymmetryManager, get_symmetry_handler
    from chimerax.atomic import AtomicStructure

    if structure is None or len(structure) == 0:
        sym_mgrs = [m for m in session.models if isinstance(m, SymmetryManager)]
        if not sym_mgrs:
            raise UserError('isolde brefine: no Clipper sessions found — '
                            'load structure factors first with "clipper open".')
        if len(sym_mgrs) > 1:
            raise UserError('isolde brefine: multiple Clipper sessions found; '
                            'please specify a model (e.g. "isolde brefine #1").')
        sym_mgr = sym_mgrs[0]
    elif len(structure) > 1:
        raise UserError('isolde brefine: specify a single structure or Clipper model.')
    else:
        item = structure[0]
        if isinstance(item, SymmetryManager):
            sym_mgr = item
        elif isinstance(item, AtomicStructure):
            sym_mgr = get_symmetry_handler(item)
            if sym_mgr is None:
                raise UserError('isolde brefine: no Clipper session is associated '
                                'with the selected structure — load structure '
                                'factors first with "clipper open".')
        else:
            raise UserError(f'isolde brefine: expected an atomic structure or '
                            f'Clipper model, got {type(item).__name__}.')

    atoms = sym_mgr.structure.atoms
    xmapset = next((x for x in sym_mgr.map_mgr.xmapsets
                    if x.live_xmap_mgr is not None), None)
    if xmapset is None:
        raise UserError('isolde brefine: no live crystallographic map set found — '
                        'the structure factors must be loaded with a live '
                        '(auto-updating) map to drive refinement.')

    level_map = {'none': None, 'info': 'info', 'debug': 'debug'}
    if log_level not in level_map:
        raise UserError(f'isolde brefine: logLevel must be none, info, or debug '
                        f'(got "{log_level}").')

    b_restraints, _, diag = build_b_restraints_for_refinement(
        sym_mgr, atoms, None, xmapset, ignore_hydrogens=ignore_hydrogens,
        crystallographic=True, internal_weight=w_internal,
        local_weight=w_local, local_cutoff=c_local, alpha=alpha)
    _log_restraint_counts(session, 'isolde brefine', b_restraints, None, diag, alpha)

    from chimerax.clipper.clipper_python.ext import RefineConfig
    cfg = RefineConfig()
    cfg.refine_b   = refine_b
    cfg.max_cycles = max_cycles
    cfg.b_max      = b_max

    from chimerax.clipper.refine.bfactor_occ import BFactorOccRefineManager
    mgr = BFactorOccRefineManager(
        session, sym_mgr, xmapset, config=cfg,
        b_min_resolution_factor=b_min_factor, n_threads=n_threads,
        rwork_tolerance=rwork_tolerance, log_level=level_map[log_level])
    mgr.launch(atoms, ignore_hydrogens=ignore_hydrogens,
               auto_tolerances=auto_tolerances,
               delta_tolerance_factor=delta_tolerance_factor,
               epsilon_tolerance_factor=epsilon_tolerance_factor,
               b_restraints=b_restraints)


def isolde_brsr_cmd(session, atoms=None, map=None,
                    context_range=5.0, whole_map=False,
                    padding=6.0, taper_width=3.0,
                    max_cycles=50, b_min_factor=2.0, b_max=200.0,
                    w_boundary=1.0,
                    c_boundary=4.0,
                    w_internal=1.0,
                    n_threads=None, log_level='info', resolution=None,
                    auto_tolerances=True, delta_tolerance_factor=1.0,
                    epsilon_tolerance_factor=1.0, ignore_hydrogens=True,
                    alpha=ALPHA_DEFAULT):
    from chimerax.core.errors import UserError
    from chimerax.clipper.symmetry import SymmetryManager, get_symmetry_handler

    if atoms is not None:
        structures = list(atoms.unique_structures)
        if len(structures) != 1:
            raise UserError('isolde brsr: all specified atoms must belong to a '
                            'single structure.')
        sym_mgr = get_symmetry_handler(structures[0])
        if sym_mgr is None:
            raise UserError('isolde brsr: no Clipper session is associated with '
                            'the structure containing the specified atoms.')
        refined_atoms = atoms
    else:
        sym_mgrs = [m for m in session.models if isinstance(m, SymmetryManager)]
        if not sym_mgrs:
            raise UserError('isolde brsr: no Clipper sessions found — '
                            'open a structure with "clipper open" first.')
        if len(sym_mgrs) > 1:
            raise UserError('isolde brsr: multiple Clipper sessions found; '
                            'specify atoms to identify the target structure '
                            '(e.g. "isolde brsr #1").')
        sym_mgr = sym_mgrs[0]
        refined_atoms = sym_mgr.structure.atoms

    context_atoms = None
    if context_range > 0.0:
        import numpy
        from chimerax.geometry import find_close_points
        non_refined = sym_mgr.structure.atoms.subtract(refined_atoms)
        if len(non_refined) > 0:
            i_nr, _ = find_close_points(non_refined.coords,
                                        refined_atoms.coords, context_range)
            if len(i_nr):
                context_atoms = non_refined[numpy.unique(i_nr)]

    from chimerax.clipper.maps.map_handler_base import MapHandlerBase, XmapHandlerBase
    if map is not None:
        if len(map) != 1:
            raise UserError(f'isolde brsr: specify exactly one map handler '
                            f'(got {len(map)}).')
        target_vol = map[0]
    else:
        all_handlers = [m for m in session.models.list()
                        if isinstance(m, MapHandlerBase)]
        if not all_handlers:
            raise UserError('isolde brsr: no Clipper map handlers found in the '
                            'session; specify the target map with "map #N".')
        if len(all_handlers) > 1:
            names = ', '.join(f'{m.id_string} ({m.name})' for m in all_handlers)
            raise UserError(f'isolde brsr: multiple map handlers found ({names}); '
                            'specify the target map with "map #N".')
        target_vol = all_handlers[0]

    if whole_map and isinstance(target_vol, XmapHandlerBase):
        raise UserError('isolde brsr: wholeMap is only valid for non-'
                        'crystallographic (cryo-EM) maps.  Crystallographic live '
                        'maps cover only the current spotlight box and cannot be '
                        'used in whole-map mode.')

    level_map = {'none': None, 'info': 'info', 'debug': 'debug'}
    if log_level not in level_map:
        raise UserError(f'isolde brsr: logLevel must be none, info, or debug '
                        f'(got "{log_level}").')

    b_restraints, b_target_restraints, diag = build_b_restraints_for_refinement(
        sym_mgr, refined_atoms, context_atoms, target_vol,
        ignore_hydrogens=ignore_hydrogens,
        crystallographic=False, internal_weight=w_internal,
        boundary_weight=w_boundary,
        boundary_cutoff=c_boundary, alpha=alpha)
    _log_restraint_counts(session, 'isolde brsr', b_restraints,
                          b_target_restraints, diag, alpha)

    from chimerax.clipper.clipper_python.ext import RefineConfig
    cfg = RefineConfig()
    cfg.max_cycles = max_cycles
    cfg.b_max      = b_max

    from chimerax.clipper.refine.bfactor_occ import BFactorOccRefineManager
    mgr = BFactorOccRefineManager(
        session, sym_mgr, sym_mgr.map_mgr.nxmapset, config=cfg,
        b_min_resolution_factor=b_min_factor, n_threads=n_threads,
        log_level=level_map[log_level])
    mgr.launch_realspace(
        refined_atoms, context_atoms=context_atoms, target_volume=target_vol,
        whole_map=whole_map, padding=padding, taper_width=taper_width,
        ignore_hydrogens=ignore_hydrogens, resolution=resolution,
        auto_tolerances=auto_tolerances,
        delta_tolerance_factor=delta_tolerance_factor,
        epsilon_tolerance_factor=epsilon_tolerance_factor,
        b_restraints=b_restraints,
        b_target_restraints=b_target_restraints)


def register_bfactor_refine_commands(logger):
    from chimerax.core.commands import (
        register, CmdDesc, BoolArg, FloatArg, IntArg, StringArg, Or, EnumOf)
    from chimerax.atomic import AtomsArg

    # Barron loss shape: a float (broadcast) or the literal 'auto' (per-restraint
    # confidence).  EnumOf is tried first so 'auto' parses as the string.
    AlphaArg = Or(EnumOf(['auto']), FloatArg)

    brefine_desc = CmdDesc(
        optional=[('structure', _StructuresOrSymmetryMgrsArg)],
        keyword=[
            ('refine_b', BoolArg), ('max_cycles', IntArg),
            ('b_min_factor', FloatArg), ('b_max', FloatArg),
            ('n_threads', IntArg), ('log_level', StringArg),
            ('rwork_tolerance', FloatArg), ('auto_tolerances', BoolArg),
            ('delta_tolerance_factor', FloatArg),
            ('epsilon_tolerance_factor', FloatArg),
            ('ignore_hydrogens', BoolArg),
            ('w_internal', FloatArg),
            ('w_local', FloatArg),
            ('c_local', FloatArg),
            ('alpha', AlphaArg),
        ],
        synopsis='Refine B-factors against crystallographic data using the '
                 'ML gradient-map approach')
    register('isolde brefine', brefine_desc, isolde_brefine_cmd, logger=logger)

    brsr_desc = CmdDesc(
        optional=[('atoms', AtomsArg)],
        keyword=[
            ('map', _MapHandlerArg), ('context_range', FloatArg),
            ('whole_map', BoolArg), ('padding', FloatArg),
            ('taper_width', FloatArg), ('max_cycles', IntArg),
            ('b_min_factor', FloatArg), ('b_max', FloatArg),
            ('n_threads', IntArg), ('log_level', StringArg),
            ('resolution', FloatArg), ('auto_tolerances', BoolArg),
            ('delta_tolerance_factor', FloatArg),
            ('epsilon_tolerance_factor', FloatArg),
            ('ignore_hydrogens', BoolArg),
            ('w_boundary', FloatArg),
            ('c_boundary', FloatArg),
            ('w_internal', FloatArg),
            ('alpha', AlphaArg),
        ],
        synopsis='Real-space B-factor refinement against a density map')
    register('isolde brsr', brsr_desc, isolde_brsr_cmd, logger=logger)
