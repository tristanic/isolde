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


def chiral_definitions_from_ccd(session, resname):
    '''Generate intra-residue chiral-centre definitions for a CCD component.

    Returns a dict
    ``{centre_name: ([s1], [s2], [s3], expected_angle_rad, signed_volume)}`` (the
    substituent names are single-element lists, matching ``add_chiral_def``'s
    "list of possible names" convention). The signed volume is carried for
    diagnostics / the ordering self-check; the registered restraint derives its
    target sign from the angle. Only centres whose substituent atoms are all
    within the residue are produced; cross-residue (external) centres are handled
    separately.
    '''
    mol, _status = rb.template_to_rdkit(session, resname)
    if mol is None or mol.GetNumConformers() == 0:
        return {}
    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception:
        return {}
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
        pts = [_conf_coord(conf, atom.GetIdx())]
        pts += [_conf_coord(conf, n.GetIdx()) for n in subs]
        # Same improper ISOLDE's ChiralCenter computes: dihedral(centre,s1,s2,s3).
        angle = radians(dihedral(*pts))
        volume = _signed_volume(pts)
        centre_name = atom.GetProp(rb.NAME_PROP)
        s_names = [n.GetProp(rb.NAME_PROP) for n in subs]
        defs[centre_name] = ([s_names[0]], [s_names[1]], [s_names[2]], angle, volume)
    return defs


def register_chiral_definitions(session, resname, override=False):
    '''Register generated chiral definitions for ``resname`` with the session
    ``ChiralMgr`` so they are restrained during simulation. Idempotent: centres
    already defined (e.g. from ``chirals.json``) are skipped unless ``override``.
    Returns the number of new definitions added.'''
    from chimerax.isolde.molobject import get_chiral_mgr
    mgr = get_chiral_mgr(session)
    existing = mgr.chiral_center_dict.get(resname, {})
    added = 0
    for centre, (s1, s2, s3, angle, volume) in chiral_definitions_from_ccd(session, resname).items():
        if centre in existing and not override:
            continue
        try:
            mgr.add_chiral_def(resname, centre, s1, s2, s3, angle,
                expected_volume=volume)
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
        (CHIRAL_OUTLIER_THRESHOLD - oriented) / CHIRAL_OUTLIER_THRESHOLD, 0.0, 1.0)
    return chirals[mask], oriented[mask], severity[mask]


def ensure_chiral_definitions(session, residues):
    '''For each unique residue name among ``residues`` that has no chiral
    definition yet (not in the curated dictionary, not already generated this
    session), generate and register one from the CCD so its chirality is
    restrained during simulation.

    Defensive by design -- it never raises, since chirality restraints are
    best-effort and must not break simulation setup. Each novel residue name is
    attempted at most once per session.
    '''
    from chimerax.isolde.molobject import get_chiral_mgr
    mgr = get_chiral_mgr(session)
    attempted = getattr(session, '_isolde_chiral_attempted', None)
    if attempted is None:
        attempted = session._isolde_chiral_attempted = set()
    for name in set(r.name for r in residues):
        if name in attempted or name in mgr.chiral_center_dict:
            continue
        attempted.add(name)
        try:
            n = register_chiral_definitions(session, name)
            if n:
                session.logger.info(
                    f'ISOLDE: generated {n} chiral restraint definition(s) for '
                    f'residue {name}.')
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
            gs1, gs2, gs3, gangle, gvol = gen[centre]
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
