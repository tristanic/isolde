# @Author: Tristan Croll <tic20>
# @Date:   10-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 11-Jun-2026
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

#: Per-type Geman-McClure half-max scale σ (U-space, Å²).  Bonded atoms share B
#: closely → tight σ (σ_B ≈ 5 Å²); non-bonded context contacts can differ more →
#: wider σ (σ_B ≈ 15 Å²) so genuine B-gradients/outliers are still pulled in.
_SIGMA_BONDED = 0.063    # σ_B ≈ 5 Å²
_SIGMA_CONTEXT = 0.19    # σ_B ≈ 15 Å²

#: Distance kernel scale (Å): w(d) = ε²/(d²+ε²) — ≈1 at contact, decays outward.
_DIST_KERNEL_EPS = 3.5

#: Master normalization constants ("fraction of data energy a restraint family may
#: spend"; see build_b_restraints_for_refinement).  These bridge the unitless data
#: energy to the engine's restraint weight scale and are intended to be CALIBRATED
#: ONCE from a logged run: pick η so a user weight of ~1.0 gives gentle-but-active
#: regularization (reproduces the empirically-good raw ≈ 0.001 at 3 Å).  Because
#: the per-run data energy S already tracks dataset size/|Fo|², the calibrated η is
#: then portable across resolutions.  First guesses only — refine from the log.
_ETA_XTAL = 1.0e-7       # crystallographic (brefine); calibrated, sweet spot ~ w=1
_ETA_RS = 10.0           # real-space (brsr); calibrated (sweet spot was ~ w=10 at eta=1.0)


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


def build_bonded_b_restraints(refined_atoms, table, ignore_hydrogens):
    '''Two-sided pairwise restraints along covalent bonds (heavy atoms only).

    The ``ks`` column holds **relative** weights (ω = 1 for every bond); the
    orchestrator rescales them to the data-normalized family budget.  Returns a
    ``(atoms1, altlocs1, atoms2, altlocs2, sigmas, ks)`` tuple, or ``None``.
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
    atoms1, altlocs1, atoms2, altlocs2, sigmas, ks = [], [], [], [], [], []
    seen = set()
    for b in range(len(bonds)):
        i, j = int(idx1[b]), int(idx2[b])
        if i < 0 or j < 0:
            continue
        if el1[b] == 1 or el2[b] == 1:   # heavy-only
            continue
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
    if not atoms1:
        return None
    return (atoms1, altlocs1, atoms2, altlocs2, sigmas, ks)


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
                                       ignore_hydrogens, cutoff):
    '''One-sided restraints pulling each refined atom toward the distance-weighted
    mean B of its (symmetry-aware) context neighbours.  The ``ks`` column holds
    relative weights (ω = 1 per target); the orchestrator rescales them.

    Returns a ``(atoms, altlocs, target_us, sigmas, ks)`` tuple, or ``None``.
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

    atoms, altlocs, target_us = [], [], []
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
        atoms.append(int(refined_table.refined_index[rr]))
        altlocs.append(str(refined_table.altloc[rr]))
        target_us.append(float((w * us).sum() / w.sum()))
    if not atoms:
        return None

    sigmas = [_SIGMA_CONTEXT] * len(atoms)
    ks = [1.0] * len(atoms)   # relative weight ω; orchestrator rescales
    return (atoms, altlocs, target_us, sigmas, ks)


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
                                     ignore_hydrogens, cutoff, exclude_bonded=True):
    '''Two-sided per-pair restraints between near (non-bonded) refined atoms,
    including crystallographic symmetry contacts.  The ``ks`` column holds
    relative weights ω = w(d)/√(K_i·K_j) (distance kernel ÷ the endpoints' contact
    counts); the orchestrator rescales them.  Returns a 6-tuple or ``None``.

    ``exclude_bonded`` — when the bonded family is also active, covalently-bonded
    pairs are handled there (with their tighter σ), so they are dropped here to
    avoid double restraints; the orchestrator passes ``False`` when bonded
    restraints are off, so those near contacts are still covered locally.

    Cost note: over a whole structure this enumerates all heavy contacts within
    ``cutoff`` (× symmetry operators), so it is opt-in (default weight 0).
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

    atoms1, altlocs1, atoms2, altlocs2, sigmas, ks = [], [], [], [], [], []
    for (a, b), dist in best.items():
        mi, mj = int(master[a]), int(master[b])
        omega = float(_distance_weight(dist)) / ((contacts[mi] * contacts[mj]) ** 0.5)
        atoms1.append(mi)
        altlocs1.append(str(altl[a]))
        atoms2.append(mj)
        altlocs2.append(str(altl[b]))
        sigmas.append(_SIGMA_CONTEXT)
        ks.append(omega)   # relative weight ω; orchestrator rescales
    return (atoms1, altlocs1, atoms2, altlocs2, sigmas, ks)


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
    a1, al1, a2, al2, sig, ks = b_restraints
    keep = [k for k in range(len(a1)) if a1[k] != a2[k]]
    if len(keep) == len(a1):
        return b_restraints
    if not keep:
        return None
    return tuple([seq[k] for k in keep] for seq in (a1, al1, a2, al2, sig, ks))


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
    '''Rescale a family's *relative* ω weights (the tuple's last column) to the
    data-normalized budget: ``k_eff(r) = w_user·η·S·ω_r / Σω``.  So the family's
    total weight ≈ ``w_user·η·S`` regardless of restraint count, and per-restraint
    strength tracks the data energy ``S`` (→ portable).  Returns the tuple, or
    ``None`` if the family is empty/disabled.'''
    if restraints is None or w_user is None or w_user <= 0:
        return None
    cols = list(restraints)
    omega = cols[-1]
    total = float(sum(omega))
    if total <= 0:
        return None
    scale = float(w_user) * float(eta) * float(S) / total
    cols[-1] = [float(o) * scale for o in omega]
    return tuple(cols)


def build_b_restraints_for_refinement(sym_mgr, refined_atoms, context_atoms,
                                      data_source, *,
                                      ignore_hydrogens, crystallographic,
                                      internal_weight=0.0,
                                      boundary_weight=0.0, boundary_cutoff=0.0,
                                      local_weight=0.0, local_cutoff=0.0):
    '''Build the restraint tuples for a refinement run, with per-family data
    normalization (see :func:`_normalize_family`).

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
            build_bonded_b_restraints(refined_atoms, table, ignore_hydrogens),
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
                    exclude_bonded=bonded_active),
                local_weight, eta, S)
            if local is not None:
                diag['budgets']['local'] = local_weight * eta * S
        b_restraints = _strip_self_pairs(_concat_restraint_tuples(bonded, local, 6))
        return b_restraints, None, diag

    boundary = None
    if boundary_weight and boundary_weight > 0:
        boundary = _normalize_family(
            build_boundary_b_target_restraints(
                sym_mgr, refined_atoms, table, ignore_hydrogens, boundary_cutoff),
            boundary_weight, eta, S)
        if boundary is not None:
            diag['budgets']['boundary'] = boundary_weight * eta * S
    return _strip_self_pairs(bonded), boundary, diag


def _log_restraint_counts(session, command, b_restraints, b_target_restraints, diag):
    n_pair = len(b_restraints[0]) if b_restraints is not None else 0
    n_targ = len(b_target_restraints[0]) if b_target_restraints is not None else 0
    msg = f'{command}: {n_pair} pairwise B-factor restraint(s)'
    if b_target_restraints is not None:
        msg += f', {n_targ} one-sided target restraint(s)'

    def _maxk(t):
        return max(t[-1]) if (t is not None and len(t[-1])) else 0.0
    rep = max(_maxk(b_restraints), _maxk(b_target_restraints))
    budgets = ', '.join(f'{k}={v:.4g}' for k, v in diag['budgets'].items()) or 'none'
    msg += (f'.  [S_data={diag["S_data"]:.4g}, eta={diag["eta"]:.3g}; '
            f'family budgets: {budgets}; max k_eff={rep:.4g}]')
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
                       w_internal=0.0,
                       w_local=0.0,
                       c_local=4.0):
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
        local_weight=w_local, local_cutoff=c_local)
    _log_restraint_counts(session, 'isolde brefine', b_restraints, None, diag)

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
                    epsilon_tolerance_factor=1.0, ignore_hydrogens=True):
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
        boundary_cutoff=c_boundary)
    _log_restraint_counts(session, 'isolde brsr', b_restraints,
                          b_target_restraints, diag)

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
        register, CmdDesc, BoolArg, FloatArg, IntArg, StringArg)
    from chimerax.atomic import AtomsArg

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
        ],
        synopsis='Real-space B-factor refinement against a density map')
    register('isolde brsr', brsr_desc, isolde_brsr_cmd, logger=logger)
