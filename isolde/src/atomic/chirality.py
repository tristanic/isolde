# @Author: Tristan Croll
# @Date:   20-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 20-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Automatic generation of chiral-centre restraint definitions for arbitrary
residues, from RDKit CIP perception of a reference (CCD ideal) structure.

ISOLDE's ``ChiralMgr`` restrains chirality during simulation, but only for
centres pre-defined in the hand-curated ``chirals.json`` dictionary. It has had
no way to *determine* the absolute configuration of an arbitrary ligand. RDKit
supplies exactly that, and its output maps directly onto ISOLDE's existing
machinery:

``ChiralCenter`` is a ``Dihedral(centre, s1, s2, s3)`` improper whose **sign
encodes handedness** (S -> positive, R -> negative), where ``s1, s2, s3`` are the
three highest-CIP-priority substituents, highest first (the lowest-priority
substituent -- almost always the H -- is omitted). RDKit's canonical CIP ranks
reproduce that ordering automatically (e.g. for CYS the S on CB outranks the
carbonyl O, reordering the substituents and flipping the curated sign), so we
compute the same improper from the reference coordinates and hand it to
``ChiralMgr.add_chiral_def``.

The :func:`validate_against_curated` gate regenerates the hand-curated entries
and diffs them, proving the sign/ordering convention before the generator is
trusted on novel ligands.
'''

from math import radians, pi
import numpy

from rdkit import Chem

from . import rdkit_bridge as rb


def _conf_coord(conf, idx):
    p = conf.GetAtomPosition(idx)
    return numpy.array([p.x, p.y, p.z])


def _signed_volume(pts):
    '''Signed chiral volume V=(s1-c).[(s2-c)x(s3-c)] in Angstrom^3, from
    ``[centre, s1, s2, s3]`` coordinates. The vector ordering is identical to
    ``ChiralCenter::chiral_volume()`` (C++) and the OpenMM
    ``ChiralVolumeRestraintForce`` expression -- keep all three in lock-step.'''
    c, s1, s2, s3 = (numpy.asarray(p, dtype=float) for p in pts)
    return float(numpy.dot(s1 - c, numpy.cross(s2 - c, s3 - c)))


def _leaving_atom_flags(session, resname):
    '''``{atom_id: bool}`` from the CCD ``pdbx_leaving_atom_flag`` (via the local
    ChemComp store). Empty if unavailable -- callers then treat every centre as
    intra-residue (the pre-6a behaviour).'''
    try:
        from chimerax import chemcomp
        rec = chemcomp.record(session, resname)
    except Exception:
        return {}
    if not rec:
        return {}
    return {a[0]: (len(a) > 5 and a[5] == 'Y') for a in rec['atoms']}


def chiral_definitions_from_ccd(session, resname):
    '''Generate chiral-centre definitions for a CCD component.

    Returns a dict
    ``{centre_name: ([s1], [s2], [s3], expected_angle_rad, signed_volume,
    externals)}`` -- substituent names as single-element lists (matching
    ``add_chiral_def``'s "list of possible names" convention) and ``externals`` a
    3-list of 0/1 flags.

    Cross-residue (linkage) centres are handled: if one of the three
    highest-priority substituents is a CCD polymerisation *leaving* atom
    (``pdbx_leaving_atom_flag``), it is replaced by the wildcard ``'*'`` and its
    slot flagged external, so the def matches the *linked* residue -- where that
    atom is gone, replaced by the bond to the neighbouring residue (e.g. a
    glycosidic anomeric centre). The sign/volume come from the monomer reference
    geometry (the leaving atom sits where the linkage partner will). Centres with
    two leaving atoms among the top three are skipped: the single-external model
    can't represent them.
    '''
    mol, _status = rb.template_to_rdkit(session, resname)
    if mol is None or mol.GetNumConformers() == 0:
        return {}
    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception:
        return {}
    leaving = _leaving_atom_flags(session, resname)
    conf = mol.GetConformer()
    from chimerax.geometry import dihedral
    defs = {}
    for atom in mol.GetAtoms():
        if not atom.HasProp('_CIPCode') or not atom.HasProp(rb.NAME_PROP):
            continue
        nbrs = list(atom.GetNeighbors())
        if len(nbrs) < 3:
            continue

        def _rank(a):
            return a.GetIntProp('_CIPRank') if a.HasProp('_CIPRank') else -1
        # Highest CIP priority first; drop the lowest-priority substituent.
        nbrs.sort(key=_rank, reverse=True)
        subs = nbrs[:3]
        if any(not n.HasProp(rb.NAME_PROP) for n in subs):
            continue
        s_names = [n.GetProp(rb.NAME_PROP) for n in subs]
        ext_mask = [bool(leaving.get(nm)) for nm in s_names]
        if sum(ext_mask) >= 2:
            # >1 leaving atom among the top-3 substituents: the single external
            # wildcard can't represent this, so skip rather than mis-define.
            continue
        # The C++ Chiral_Def only permits the wildcard in the final (s3) slot, so
        # put any external substituent last (keeping the other two in CIP order)
        # and compute the target improper for THAT ordering -- the sign stays
        # self-consistent because the C++/OpenMM volume uses the same order.
        order = (
            [i for i in range(3) if not ext_mask[i]] + [i for i in range(3) if ext_mask[i]]
        )
        ordered = [subs[i] for i in order]
        pts = [_conf_coord(conf, atom.GetIdx())]
        pts += [_conf_coord(conf, n.GetIdx()) for n in ordered]
        # Computed on the monomer (leaving atom in place); its sign is the
        # handedness the linked partner atom will occupy.
        angle = radians(dihedral(*pts))
        volume = _signed_volume(pts)
        centre_name = atom.GetProp(rb.NAME_PROP)
        names = [s_names[i] for i in order]
        externals = [0, 0, 0]
        if any(ext_mask):
            # The single external (now last) is absent in the linked residue;
            # name it with the cross-residue wildcard.
            names[2] = '*'
            externals = [0, 0, 1]
        defs[centre_name] = ([names[0]], [names[1]], [names[2]], angle, volume, externals)
    return defs


def lookup_name(residue):
    '''The ChemComp identifier under which ``residue``'s chemistry is stored: its
    ``isolde_chemcomp_id`` when set (e.g. a ``LIG01`` residue placed from a
    registered compound whose real id is long/proprietary), else the residue name
    (the usual case -- standard CCD codes are their own identifiers).'''
    return getattr(residue, 'isolde_chemcomp_id', None) or residue.name


def register_chiral_definitions(session, resname, lookup_id=None, override=False):
    '''Register generated chiral definitions for ``resname`` with the session
    ``ChiralMgr`` so they are restrained during simulation. Idempotent: centres
    already defined (e.g. from ``chirals.json``) are skipped unless ``override``.

    Definitions are *registered under* ``resname`` (the in-model residue name the
    running simulation matches on), but the source chemistry is fetched under
    ``lookup_id`` (default ``resname``) -- so a ``LIG01`` residue placed from a
    registered compound gets restraints derived from that compound's geometry.
    Returns the number of new definitions added.'''
    from chimerax.isolde.molobject import get_chiral_mgr
    mgr = get_chiral_mgr(session)
    existing = mgr.chiral_center_dict.get(resname, {})
    added = 0
    for centre, (s1, s2, s3, angle, volume, externals) in \
            chiral_definitions_from_ccd(session, lookup_id or resname).items():
        if centre in existing and not override:
            continue
        try:
            mgr.add_chiral_def(
                resname,
                centre,
                s1,
                s2,
                s3,
                angle,
                externals=externals,
                expected_volume=volume
            )
            added += 1
        except Exception as e:
            session.logger.info(
                f'Chiral definition for {resname} {centre} skipped: {e}')
    return added


# Oriented four-substituent volume (Angstrom^3) below which a centre is flagged
# as a validation outlier. A correct tetrahedral centre has a true volume well
# above this (~8-10 A^3); values near zero are near-planar and negative values
# are inverted, so this small positive threshold flags inverted and severely
# strained centres without false-flagging healthy ones.
CHIRAL_OUTLIER_THRESHOLD = 0.5  # Angstrom^3


def chiral_outliers(session, atoms, ensure_defs=True):
    '''Identify chiral centres among `atoms` whose handedness is wrong or badly
    strained, from the centre geometry alone (works on a static model -- it does
    not require a running simulation / chiral restraints).

    Uses the robust **four-substituent** signed volume
    (:attr:`ChiralCenter.true_chiral_volume`), which -- unlike the restraint's
    centre-plus-three metric -- cannot be fooled by the lowest-priority
    substituent being stuck on the wrong side. The target handedness sign comes
    from :attr:`expected_volume` (for a correct centre sign(V4) == sign(V3)).

    Returns ``(chirals, oriented, severity)``:
      * ``chirals``  -- the outlier :class:`ChiralCenters` (a sub-selection)
      * ``oriented`` -- their oriented true volume (sign(target)*V4; <0 = inverted)
      * ``severity`` -- per-centre severity in [0, 1] (0 at the threshold, 1 once
        the centre is planar or inverted)
    '''
    import numpy
    from chimerax.isolde.molobject import get_chiral_mgr
    if ensure_defs:
        ensure_chiral_definitions(session, atoms.unique_residues)
    mgr = get_chiral_mgr(session)
    chirals = mgr.get_chirals(atoms, create=True)
    empty = numpy.array([], dtype=float)
    if not len(chirals):
        return chirals, empty, empty
    v4 = chirals.true_chiral_volumes
    ev = chirals.expected_volumes
    oriented = numpy.sign(ev) * v4
    mask = oriented < CHIRAL_OUTLIER_THRESHOLD
    severity = numpy.clip(
        (CHIRAL_OUTLIER_THRESHOLD - oriented) / CHIRAL_OUTLIER_THRESHOLD, 0.0, 1.0
    )
    return chirals[mask], oriented[mask], severity[mask]


def _reference_data(session, lookup_id):
    '''``(cip_codes, ideal_coords)`` for a ChemComp component, both from one RDKit
    pass over its reference (CCD ideal) geometry:

      * ``cip_codes``    -- ``{atom_name: 'R'|'S'}`` (true stereocentres only)
      * ``ideal_coords`` -- ``{atom_name: numpy xyz}`` of the ideal conformer

    Works uniformly for standard CCD codes and registered (store-only) ligands
    (the same path :func:`chiral_definitions_from_ccd` uses). Both are fixed for a
    given component, so callers should cache the result per id. Empty dicts if the
    component is unknown.'''
    codes, coords = {}, {}
    mol, _status = rb.template_to_rdkit(session, lookup_id)
    if mol is None or mol.GetNumConformers() == 0:
        return codes, coords
    conf = mol.GetConformer()
    for a in mol.GetAtoms():
        if a.HasProp(rb.NAME_PROP):
            coords[a.GetProp(rb.NAME_PROP)] = _conf_coord(conf, a.GetIdx())
    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception:
        pass
    for a in mol.GetAtoms():
        if a.HasProp(rb.NAME_PROP) and a.HasProp('_CIPCode'):
            c = a.GetProp('_CIPCode')
            if c in ('R', 'S'):
                codes[a.GetProp(rb.NAME_PROP)] = c
    return codes, coords


def reference_cip_codes(session, lookup_id):
    '''Reference absolute configuration ``{atom_name: 'R'|'S'}`` for the named
    ChemComp component (RDKit CIP on the reference geometry). Thin wrapper over
    :func:`_reference_data`.'''
    return _reference_data(session, lookup_id)[0]


def _reference_volume_sign(cc, ideal_coords):
    '''Sign (+1.0/-1.0/0.0) of the reference signed chiral volume for centre
    ``cc``, from ``ideal_coords`` (by atom name) in the centre's own
    ``[centre, s1, s2, s3]`` ordering. A **flip-immune** stand-in for
    :attr:`ChiralCenter.expected_volume` (which ``chiralflip force`` mutates).
    Returns 0.0 if any of the four atoms is absent from the reference (e.g. a
    cross-residue wildcard substituent), so the caller can skip it.'''
    atoms = cc.atoms
    try:
        pts = [ideal_coords[atoms[j].name] for j in range(4)]
    except KeyError:
        return 0.0
    v = _signed_volume(pts)
    return 1.0 if v > 0 else (-1.0 if v < 0 else 0.0)


def as_modelled_configs(session, chirals):
    '''As-modelled absolute configuration (``'R'``/``'S'``/``None``) for each
    centre in ``chirals`` (a :class:`ChiralCenters`).

    Depends ONLY on the live geometry: the as-modelled handedness sign
    (``sign(true_chiral_volume)``) is compared against the centre's **reference**
    handedness sign, computed from the CCD ideal coordinates in the *same*
    substituent ordering (:func:`_reference_volume_sign`); the reference CIP letter
    is then flipped iff the two signs differ.

    This deliberately avoids :attr:`ChiralCenter.expected_volume`: the experts-only
    ``isolde chiralflip ... force true`` *retargets* a centre to its enantiomer by
    negating ``expected_volume``, so a test based on it would wrongly report the
    flipped geometry as still matching its (also-flipped) target. The *as-modelled*
    label must track the physical geometry, not the restraint target.

    Returns ``None`` where no reference is available (unknown component, missing
    CIP code, or a cross-residue wildcard centre); the caller draws nothing there.
    Computed only when the chiral set / geometry changes, not per frame.'''
    n = len(chirals)
    if not n:
        return []
    cas = chirals.chiral_atoms
    residues = cas.residues
    names = cas.names
    v4 = chirals.true_chiral_volumes
    ref_cache = {}
    out = []
    for i in range(n):
        lid = lookup_name(residues[i])
        ref = ref_cache.get(lid)
        if ref is None:
            ref = ref_cache[lid] = _reference_data(session, lid)
        codes, coords = ref
        letter = codes.get(names[i])
        if letter is None:
            out.append(None)
            continue
        ref_sign = _reference_volume_sign(chirals[i], coords)
        if ref_sign == 0.0:
            out.append(None)
            continue
        if ref_sign * v4[i] < 0:   # live handedness opposite the reference -> flipped
            letter = 'S' if letter == 'R' else 'R'
        out.append(letter)
    return out


def ensure_chiral_definitions(session, residues):
    '''For each unique residue name among ``residues``, generate chiral-centre
    definitions from the CCD and register any centre not already defined, so its
    chirality is restrained during simulation.

    This *supplements* the hand-curated dictionary rather than deferring to it:
    :func:`register_chiral_definitions` keeps any curated centre (skips it per
    atom) and only adds the missing ones. That matters for the curated sugars,
    which define only the anomeric carbon -- their ring stereocentres would
    otherwise go unrestrained.

    Defensive by design -- it never raises, since chirality restraints are
    best-effort and must not break simulation setup. Each residue name is
    attempted at most once per session.
    '''
    attempted = getattr(session, '_isolde_chiral_attempted', None)
    if attempted is None:
        attempted = session._isolde_chiral_attempted = set()
    # Map each in-model residue name to the ChemComp identifier its chemistry is
    # stored under (LIG01 -> the registered id; standard codes map to themselves).
    name_to_id = {}
    for r in residues:
        name_to_id.setdefault(r.name, lookup_name(r))
    for name, lookup_id in name_to_id.items():
        if name in attempted:
            continue
        attempted.add(name)
        try:
            n = register_chiral_definitions(session, name, lookup_id=lookup_id)
            if n:
                session.logger.info(
                    f'ISOLDE: generated {n} chiral restraint definition(s) for '
                    f'residue {name}.'
                )
        except Exception as e:
            session.logger.info(f'Chiral generation for {name} skipped: {e}')


def _angle_delta(a, b):
    '''Smallest signed difference a-b wrapped to (-pi, pi].'''
    return (a - b + pi) % (2 * pi) - pi


def validate_against_curated(session, resnames=None, tol=0.15):
    '''Regenerate chiral definitions from the CCD and diff them against the
    hand-curated ``chirals.json``. This is the correctness gate for the sign /
    substituent-ordering convention.

    Returns a list of discrepancy tuples
    ``(resname, centre, issue, generated, curated)``; an empty list means the
    generator reproduces the curated definitions (within ``tol`` radians on the
    angle).
    '''
    import os
    import json
    from chimerax.isolde.molobject import DICT_DIR
    with open(os.path.join(DICT_DIR, 'chirals.json')) as f:
        curated = json.load(f)
    if resnames is None:
        resnames = list(curated.keys())
    issues = []
    for resname in resnames:
        gen = chiral_definitions_from_ccd(session, resname)
        for centre, cdata in curated.get(resname, {}).items():
            cur_names, cur_angle = cdata[0], cdata[1]
            if centre not in gen:
                issues.append((resname, centre, 'missing in generated', None, cdata))
                continue
            gs1, gs2, gs3, gangle, gvol, gext = gen[centre]
            gen_names = [gs1, gs2, gs3]
            if gen_names != cur_names:
                issues.append((resname, centre, 'substituent order differs',
                               gen_names, cur_names))
            if abs(_angle_delta(gangle, cur_angle)) > tol:
                issues.append((resname, centre, 'angle differs',
                               round(gangle, 3), cur_angle))
            # Ordering self-check: the signed volume (computed with the same
            # vector ordering as the C++/OpenMM restraint) must share the sign of
            # the improper angle. A failure here means the volume formula ordering
            # has drifted out of lock-step.
            if gvol != 0.0 and (gvol > 0) != (gangle > 0):
                issues.append((resname, centre, 'volume/angle sign disagree (ordering bug)',
                               round(gvol, 3), round(gangle, 3)))
    return issues
