# @Author: Tristan Croll
# @Date:   20-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 20-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
RDKit bridge: robust interconversion between ChimeraX residues/templates and
RDKit molecules, with stereochemistry assigned from coordinates.

A ChimeraX ``Residue`` gives connectivity only -- per atom you get ``name``,
``element`` and ``coord``, plus ``intra_bonds`` -- but **no bond orders and no
formal charges**, and therefore no usable stereochemistry. ISOLDE's home-grown
graph matcher (``mcsplit``) is blind to both bond order and chirality as a
result. RDKit closes the gap, but only once we hand it a fully-specified
molecule. This module is the single place that does that, in two regimes (try
them in this order):

1. **Template match (exact, preferred).** When a reference with explicit
   chemistry exists -- a CCD component or an MD template -- build the RDKit mol
   from *its* connection table (bond orders + charges are explicit there) and map
   it onto the modelled residue by atom name. No perception required.

2. **Perceive from geometry (fallback, novel ligands).** With no reference,
   complete the valences with ChimeraX's IDATM-based ``addh``, estimate the net
   charge, and let RDKit's ``rdDetermineBonds`` perceive bond orders from the 3D
   coordinates.

In **both** regimes stereochemistry is assigned from coordinates
(``Chem.AssignStereochemistryFrom3D``) before hydrogens are stripped, so the
returned molecule carries CIP chiral tags / double-bond stereo -- the
precondition for every chirality-aware operation downstream.

Atom identity round-trips via two RDKit atom properties that survive
``RemoveHs``/``AddHs``/sanitize and perception reordering:

* ``cxName`` -- the ChimeraX atom name (== CCD atom_id where applicable);
* ``cxIdx``  -- a stable integer index into the source residue's atom list.

Build a ``{cxName: rdkit_idx}`` (or ``cxIdx``) map once and you have a stable
bidirectional correspondence even though heavy-atom counts and ordering change
through conversion.

This implementation ports the battle-tested conversion core from the
ChimeraX-ChemSearch bundle (``chem.py``: ``mol_from_ccd``, ``mol_from_geometry``,
``attach_metals``, ``fragment_by_names``; ``tool.py``: ``_perceive_residue``,
``_names_match``, ``_residue_sig``).

**Threading.** Structure/CIF operations must run on the ChimeraX main thread;
only pure-RDKit work (the ``mol_from_*`` builders, sanitize, FMCS) is safe
off-thread.
'''

from rdkit import Chem, RDLogger

# We do our own failure handling/bookkeeping, so silence RDKit's chatter.
RDLogger.DisableLog('rdApp.*')

#: Identity properties carried on RDKit atoms (survive RemoveHs/AddHs/sanitize).
NAME_PROP = 'cxName'
INDEX_PROP = 'cxIdx'
#: Parent-atom index carried on a capped charge-fragment's real atoms (see
#: :func:`build_capped_fragment`); cap atoms do not carry it.
FRAGSRC_PROP = 'fragSrc'

#: CCD ``value_order`` strings -> RDKit bond types. CCD encodes aromatic rings as
#: Kekule single/double, so we trust the explicit order and let RDKit perceive
#: aromaticity; only the rare explicit 'AROM' forces an aromatic bond.
_BOND_ORDER = {
    'SING': Chem.BondType.SINGLE,
    'DOUB': Chem.BondType.DOUBLE,
    'TRIP': Chem.BondType.TRIPLE,
    'QUAD': Chem.BondType.QUADRUPLE,
    'AROM': Chem.BondType.AROMATIC,
    'DELO': Chem.BondType.SINGLE,   # delocalised; sanitization may refine
    'PI': Chem.BondType.SINGLE,
}

#: Fallback candidate net charges for bond-order perception, tried in order when
#: no estimate is supplied: neutral first, then small magnitudes covering most
#: ligands.
_PERCEIVE_CHARGES = (0, -1, 1, -2, 2, -3, 3, -4, 4)


# ---------------------------------------------------------------------------
# Element / helper utilities
# ---------------------------------------------------------------------------

def _norm_element(symbol):
    '''CCD uses upper-case element symbols (e.g. "CL"); RDKit wants "Cl".'''
    s = (symbol or '').strip()
    if not s:
        return 'C'
    return s[0].upper() + s[1:].lower()


def _charge_candidates(preferred):
    '''Charge-search order: the supplied estimate first, then a widening spread
    around it, then the neutral-first defaults -- de-duplicated.'''
    if preferred is None:
        return list(_PERCEIVE_CHARGES)
    order = [preferred, preferred - 1, preferred + 1, preferred - 2,
             preferred + 2] + list(_PERCEIVE_CHARGES)
    out = []
    for c in order:
        if c not in out:
            out.append(c)
    return out


# ---------------------------------------------------------------------------
# Pure-RDKit builders (no ChimeraX dependency; safe off the main thread)
# ---------------------------------------------------------------------------

def mol_from_ccd(atoms, bonds, name='', coords=None, strip_hs=True):
    '''Build an RDKit molecule from explicit CCD-style atom/bond records.

    Args:
        atoms: iterable of ``(atom_id, type_symbol, charge, aromatic_flag)``.
            ``charge`` may be a string like '0'/'-1'; ``aromatic_flag`` is the raw
            CCD 'Y'/'N' (currently unused -- aromaticity is perceived).
        bonds: iterable of ``(atom_id_1, atom_id_2, value_order, aromatic_flag)``.
        name: optional molecule name.
        coords: optional ``{atom_id: (x, y, z)}``; when given, a 3D conformer is
            attached and stereochemistry is perceived from it so chirality-aware
            operations can distinguish enantiomers/diastereomers.
        strip_hs: remove explicit hydrogens from the returned mol (stereo
            survives).

    Returns:
        ``(mol, status)`` where status is one of 'ok', 'relaxed', 'failed'.
        On 'failed', ``mol`` is None. Heavy atoms carry the ``cxName`` property
        (== CCD ``atom_id``).
    '''
    rw = Chem.RWMol()
    idx = {}
    for atom_id, element, charge, arom in atoms:
        a = Chem.Atom(_norm_element(element))
        try:
            c = int(charge)
        except (TypeError, ValueError):
            c = 0
        if c:
            a.SetFormalCharge(c)
        i = rw.AddAtom(a)
        if atom_id:
            rw.GetAtomWithIdx(i).SetProp(NAME_PROP, str(atom_id))
        idx[atom_id] = i

    for a1, a2, order, arom in bonds:
        i1 = idx.get(a1)
        i2 = idx.get(a2)
        if i1 is None or i2 is None:
            continue
        if order == 'AROM':
            rw.AddBond(i1, i2, Chem.BondType.AROMATIC)
            b = rw.GetBondBetweenAtoms(i1, i2)
            b.SetIsAromatic(True)
            rw.GetAtomWithIdx(i1).SetIsAromatic(True)
            rw.GetAtomWithIdx(i2).SetIsAromatic(True)
        else:
            rw.AddBond(i1, i2, _BOND_ORDER.get(order, Chem.BondType.SINGLE))

    mol = rw.GetMol()
    if coords:
        from rdkit.Geometry import Point3D
        conf = Chem.Conformer(mol.GetNumAtoms())
        have_any = False
        for atom_id, i in idx.items():
            xyz = coords.get(atom_id)
            if xyz is not None:
                conf.SetAtomPosition(i, Point3D(*[float(v) for v in xyz]))
                have_any = True
        if have_any:
            conf.Set3D(True)
            mol.AddConformer(conf, assignId=True)
    if name:
        mol.SetProp('_Name', name)
    return _sanitize_ladder(mol, strip_hs)


def _sanitize_ladder(mol, strip_hs=True):
    '''Try full sanitization; fall back to a relaxed pass (no kekulize /
    aromaticity / property assignment); else give up. Some CCD/template
    definitions won't sanitize cleanly (quaternary N, metals).'''
    try:
        Chem.SanitizeMol(mol)
        return _finish(mol, 'ok', strip_hs)
    except Exception:
        pass
    try:
        mol.UpdatePropertyCache(strict=False)
        ops = (Chem.SanitizeFlags.SANITIZE_ALL
               ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
               ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
               ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        problems = Chem.SanitizeMol(mol, sanitizeOps=ops, catchErrors=True)
        if problems == Chem.SanitizeFlags.SANITIZE_NONE:
            return _finish(mol, 'relaxed', strip_hs)
    except Exception:
        pass
    return None, 'failed'


def _finish(mol, status, strip_hs=True):
    # Perceive stereochemistry from the 3D conformer (if any) BEFORE stripping
    # Hs, so chiral tags / double-bond stereo survive into the returned mol.
    if mol.GetNumConformers() > 0:
        try:
            Chem.AssignStereochemistryFrom3D(mol)
        except Exception:
            pass
    if strip_hs:
        try:
            mol = Chem.RemoveHs(mol)
        except Exception:
            pass  # keep Hs rather than fail outright
    return mol, status


def mol_from_geometry(symbols, coords, bonds, names=None, preferred_charge=None):
    '''Perceive an RDKit molecule (bond orders + formal charges) from atoms with
    3D coordinates and single-bond connectivity, via RDKit's ``rdDetermineBonds``.

    IMPORTANT: hydrogens MUST already be present (add them with ChimeraX ``addh``
    first). ``rdDetermineBonds`` completes valences with bond orders/charges
    rather than implicit H, so on a heavy-atom-only skeleton it invents spurious
    multiple bonds/charges. With Hs present its *global* octet optimization places
    double bonds by matching, which is what gets conjugated polyenes to alternate
    correctly instead of collapsing to a cumulene.

    Args:
        symbols: per-atom element symbols.
        coords: per-atom ``(x, y, z)``.
        bonds: ``(i, j)`` index pairs into ``symbols``/``coords``.
        names: optional per-atom names; when given, stored as the ``cxName``
            property so callers can map perceived atoms back to ChimeraX atoms.
        preferred_charge: net-charge estimate (e.g. from ChimeraX's
            ``estimate_net_charge``) tried first; the search widens around it.

    Returns:
        ``(mol_without_Hs, total_charge)``, or ``(None, None)`` if perception
        fails for every candidate charge (e.g. a coordinated metal). Atoms carry
        the integer ``cxIdx`` property (their index into the input) and, if
        ``names`` was given, ``cxName``.
    '''
    try:
        from rdkit.Chem import rdDetermineBonds
    except Exception:
        return None, None
    from rdkit.Geometry import Point3D
    rw = Chem.RWMol()
    for i, s in enumerate(symbols):
        at = Chem.Atom(_norm_element(s))
        at.SetIntProp(INDEX_PROP, i)   # survives sanitize/RemoveHs
        if names is not None and names[i]:
            at.SetProp(NAME_PROP, str(names[i]))
        rw.AddAtom(at)
    conf = Chem.Conformer(len(symbols))
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    for i, j in bonds:
        i, j = int(i), int(j)
        if i != j and rw.GetBondBetweenAtoms(i, j) is None:
            rw.AddBond(i, j, Chem.BondType.SINGLE)
    base = rw.GetMol()
    base.AddConformer(conf, assignId=True)
    for charge in _charge_candidates(preferred_charge):
        try:
            m = Chem.Mol(base)
            rdDetermineBonds.DetermineBondOrders(m, charge=charge)
            Chem.SanitizeMol(m)
            Chem.AssignStereochemistryFrom3D(m)
            return Chem.RemoveHs(m), charge
        except Exception:
            continue
    return None, None


def attach_metals(mol, metals):
    '''Re-attach coordinated metals to a perceived (metal-free) molecule.

    For a metalloligand (e.g. heme) we strip the metal, perceive the organic
    framework cleanly, then add the metal back here. ``mol`` must carry ``cxIdx``
    atom props (from ``mol_from_geometry``). ``metals`` is a list of
    ``(element_symbol, [neighbour_cxIdx, ...])``. Bonds are added DATIVE
    (donor->metal) so they don't inflate the ligand atoms' valences.
    '''
    rw = Chem.RWMol(mol)
    by_idx = {a.GetIntProp(INDEX_PROP): a.GetIdx()
              for a in rw.GetAtoms() if a.HasProp(INDEX_PROP)}
    for sym, neighbours in metals:
        mi = rw.AddAtom(Chem.Atom(_norm_element(sym)))
        for g in neighbours:
            ai = by_idx.get(g)
            if ai is not None:
                rw.AddBond(ai, mi, Chem.BondType.DATIVE)   # donor -> metal
    m = rw.GetMol()
    try:
        Chem.SanitizeMol(m)
    except Exception:
        try:
            m.UpdatePropertyCache(strict=False)
        except Exception:
            pass
    return m


def fragment_by_names(mol, keep_names):
    '''Return a copy of a CCD-built ``mol`` (heavy atoms carry ``cxName``) with
    the heavy atoms NOT in ``keep_names`` removed -- i.e. the as-modelled fragment
    of a partial residue, inheriting the CCD's correct bond orders/charges (no
    perception). Hydrogens are dropped; RDKit fills the freed valences with
    implicit H.'''
    keep = set(keep_names)
    rw = Chem.RWMol(mol)
    drop = []
    for a in rw.GetAtoms():
        if a.GetAtomicNum() == 1:
            drop.append(a.GetIdx())
            continue
        nm = a.GetProp(NAME_PROP) if a.HasProp(NAME_PROP) else None
        if nm is None or nm not in keep:
            drop.append(a.GetIdx())
    for i in sorted(drop, reverse=True):
        rw.RemoveAtom(i)
    m = rw.GetMol()
    try:
        Chem.SanitizeMol(m)
    except Exception:
        try:
            m.UpdatePropertyCache(strict=False)
        except Exception:
            pass
    return m


def rigid_fragments(mol):
    '''Partition ``mol`` into rigid fragments for block placement: the connected
    components left after cutting *rotatable* bonds (acyclic single bonds between
    two non-terminal atoms). Ring systems and their terminal substituents
    (methyls, hydroxyls, halides, ...) stay together in one fragment.

    Returns a list of sets of heavy-atom names (the ``cxName`` prop). Atoms with
    no name (e.g. the hydrogens the CCD mol omits) are not included; callers
    assign each hydrogen to its parent heavy atom's fragment.
    '''
    # RDKit's standard rotatable-bond pattern: single, not-in-ring, between two
    # non-terminal (degree > 1), non-triple-bonded atoms. On a heavy-atom-only mol
    # this keeps terminal heavy groups (a methyl C, a hydroxyl O) with their parent.
    patt = Chem.MolFromSmarts('[!$(*#*)&!D1]-!@[!$(*#*)&!D1]')
    bond_idx = []
    for a1, a2 in mol.GetSubstructMatches(patt):
        b = mol.GetBondBetweenAtoms(a1, a2)
        if b is not None:
            bond_idx.append(b.GetIdx())
    if bond_idx:
        # addDummies=False keeps atom count/order, so frag indices map back to mol.
        groups = Chem.GetMolFrags(Chem.FragmentOnBonds(mol, bond_idx, addDummies=False))
    else:
        groups = [tuple(range(mol.GetNumAtoms()))]
    out = []
    for grp in groups:
        names = set(
            mol.GetAtomWithIdx(i).GetProp(NAME_PROP) for i in grp
            if mol.GetAtomWithIdx(i).HasProp(NAME_PROP)
        )
        if names:
            out.append(names)
    return out


#: A carbon eligible as a charge-fragmentation cut endpoint: neutral, sp3,
#: non-aromatic, and non-terminal (>=2 heavy neighbours, so cutting it splits the
#: molecule rather than lopping off a lone methyl).
def _cuttable_carbon(a):
    if a.GetAtomicNum() != 6 or a.GetIsAromatic() or a.GetFormalCharge() != 0:
        return False
    if a.GetHybridization() != Chem.HybridizationType.SP3:
        return False
    heavy_deg = sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() != 1)
    return heavy_deg >= 2


def _is_polarish(a):
    '''True for an atom whose local electronics we must NOT split across a charge
    fragment: any heteroatom, any formal charge, aromatic, or bearing a non-single
    bond (a carbonyl/alkene carbon). Hydrogens are transparent.'''
    if a.GetAtomicNum() == 1:
        return False
    if a.GetAtomicNum() != 6 or a.GetFormalCharge() != 0 or a.GetIsAromatic():
        return True
    return any(b.GetBondType() != Chem.BondType.SINGLE for b in a.GetBonds())


def safe_cut_bonds(mol, buffer=2):
    '''Bond indices where the molecule may be cut for divide-and-conquer AM1-BCC
    charging without materially perturbing the charges: an acyclic single C-C bond
    between two :func:`_cuttable_carbon`\\ s, both of which are more than ``buffer``
    bonds from any "polarish" atom (:func:`_is_polarish`). This isolates cuts to
    saturated alkyl linkers and keeps every functional group's electronics intact
    within one fragment.'''
    # Multi-source BFS: graph distance from every polarish atom.
    dist = {a.GetIdx(): 0 for a in mol.GetAtoms() if _is_polarish(a)}
    frontier, d = list(dist), 0
    while frontier:
        d += 1
        nxt = []
        for i in frontier:
            for nb in mol.GetAtomWithIdx(i).GetNeighbors():
                j = nb.GetIdx()
                if j not in dist:
                    dist[j] = d
                    nxt.append(j)
        frontier = nxt
    BIG = mol.GetNumAtoms() + 1

    def far(i):
        return dist.get(i, BIG) > buffer

    out = []
    for b in mol.GetBonds():
        if b.GetBondType() != Chem.BondType.SINGLE or b.IsInRing():
            continue
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if not (_cuttable_carbon(a1) and _cuttable_carbon(a2)):
            continue
        if far(a1.GetIdx()) and far(a2.GetIdx()):
            out.append(b.GetIdx())
    return out


def charge_partition(mol, max_heavy=30, buffer=2):
    '''Partition ``mol`` into fragments of at most ~``max_heavy`` heavy atoms for
    divide-and-conquer charging, cutting ONLY at :func:`safe_cut_bonds`. Returns a
    list of sets of parent-atom indices (each set is one fragment, hydrogens
    included). A single-element list means "do not fragment" (nothing safe to cut,
    or the whole molecule already fits).

    Strategy: cut every safe bond, then greedily merge neighbouring pieces back up
    while the merged heavy-atom count stays within ``max_heavy`` -- minimising the
    number of fragments (hence caps and error). **A piece that contains any
    "polarish" atom (aromatic/heteroatom/charged/multiply-bonded) is treated as an
    irreducible core and is never merged into**, so a cheap alkyl tail can never be
    glued back onto the expensive conjugated core it was cut from; only pure-alkyl
    pieces merge with each other. A core piece larger than ``max_heavy`` (a big rigid
    ring system with no internal safe cut) is left as-is.
    '''
    n = mol.GetNumAtoms()
    heavy_total = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() != 1)
    if heavy_total <= max_heavy:
        return [set(range(n))]
    safe = safe_cut_bonds(mol, buffer=buffer)
    if not safe:
        return [set(range(n))]
    frags = Chem.GetMolFrags(Chem.FragmentOnBonds(mol, safe, addDummies=False))
    atom2frag = {}
    size = []
    core_like = []
    for fi, grp in enumerate(frags):
        for a in grp:
            atom2frag[a] = fi
        size.append(sum(1 for a in grp if mol.GetAtomWithIdx(a).GetAtomicNum() != 1))
        core_like.append(any(_is_polarish(mol.GetAtomWithIdx(a)) for a in grp))
    parent = list(range(len(frags)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    edges = []
    for bidx in safe:
        b = mol.GetBondWithIdx(bidx)
        u, v = atom2frag[b.GetBeginAtomIdx()], atom2frag[b.GetEndAtomIdx()]
        if u != v:
            edges.append((u, v))
    merged = True
    while merged:
        merged = False
        for (u, v) in edges:
            ru, rv = find(u), find(v)
            if ru == rv:
                continue
            # Never merge into/onto a core piece: keep the expensive conjugated core
            # minimal and stop it re-absorbing the cheap alkyl tail it was cut from.
            if core_like[ru] or core_like[rv]:
                continue
            if size[ru] + size[rv] <= max_heavy:
                parent[ru] = rv
                size[rv] += size[ru]
                merged = True
    groups = {}
    for fi, grp in enumerate(frags):
        groups.setdefault(find(fi), set()).update(grp)
    return list(groups.values())


def build_capped_fragment(mol, atom_indices):
    '''Build a standalone, capped RDKit fragment of ``mol`` covering ``atom_indices``
    (from :func:`charge_partition`), ready to feed to ANTECHAMBER for local AM1-BCC
    charging. Every kept (real) atom carries the ``fragSrc`` prop = its parent index,
    so charges map straight back; every bond severed by the cut is replaced by a
    **methyl carbon at the removed neighbour's coordinate** (single bond, tagged
    ``CAP``, no ``fragSrc``) so the cut carbon still sees a carbon neighbour -- the
    true local environment of an alkyl chain, minimising the charge perturbation.
    Hydrogens follow the same policy as :func:`super_residue_to_rdkit` (kept if the
    fragment is already protonated; added to caps only). Returns the mol, or ``None``
    if it will not sanitise.'''
    from rdkit.Geometry import Point3D
    keep = set(atom_indices)
    conf_in = mol.GetConformer()
    rw = Chem.RWMol()
    old2new, coords = {}, {}
    for i in sorted(keep):
        a_in = mol.GetAtomWithIdx(i)
        a = Chem.Atom(a_in.GetAtomicNum())
        a.SetFormalCharge(a_in.GetFormalCharge())
        if a_in.GetIsAromatic():
            a.SetIsAromatic(True)
        a.SetIntProp(FRAGSRC_PROP, i)
        ni = rw.AddAtom(a)
        old2new[i] = ni
        p = conf_in.GetAtomPosition(i)
        coords[ni] = (p.x, p.y, p.z)
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in keep and j in keep:
            rw.AddBond(old2new[i], old2new[j], b.GetBondType())
            if b.GetIsAromatic():
                rw.GetBondBetweenAtoms(old2new[i], old2new[j]).SetIsAromatic(True)
    cap_indices = []
    for i in keep:
        for nb in mol.GetAtomWithIdx(i).GetNeighbors():
            j = nb.GetIdx()
            if j in keep or nb.GetAtomicNum() == 1:
                continue
            cap = Chem.Atom(6)
            cap.SetProp('CAP', '1')
            ci = rw.AddAtom(cap)
            p = conf_in.GetAtomPosition(j)
            coords[ci] = (p.x, p.y, p.z)
            rw.AddBond(old2new[i], ci, Chem.BondType.SINGLE)
            cap_indices.append(ci)
    m = rw.GetMol()
    conf = Chem.Conformer(m.GetNumAtoms())
    for idx, xyz in coords.items():
        conf.SetAtomPosition(idx, Point3D(*xyz))
    conf.Set3D(True)
    m.AddConformer(conf, assignId=True)
    m, status = _sanitize_ladder(m, strip_hs=False)
    if m is None:
        return None
    frag_has_h = any(mol.GetAtomWithIdx(i).GetAtomicNum() == 1 for i in keep)
    try:
        if not frag_has_h:
            m = Chem.AddHs(m, addCoords=True)
        elif cap_indices:
            m = Chem.AddHs(m, addCoords=True, onlyOnAtoms=cap_indices)
    except Exception:
        pass
    return m


def mol_to_sdf(mol, conf_id=-1):
    '''Serialise an RDKit ``mol`` (which MUST carry a 3D conformer) to an MDL
    SDF/MOL block, ready to feed to ANTECHAMBER (``-fi sdf``).

    We hand ANTECHAMBER SDF rather than MOL2 because RDKit writes SDF natively and
    correctly -- bond orders straight from ``mol`` (kekulised), the 3D conformer
    embedded -- which avoids hand-authoring SYBYL atom types. SDF carries no atom
    names, but that is irrelevant here: the caller maps ANTECHAMBER's output back
    to source atoms by atom **order** (ANTECHAMBER preserves it 1:1) via the
    ``cxIdx`` map, then reattaches the real names. Kekule bond orders are what
    ANTECHAMBER wants; no charges are written, so it computes AM1-BCC itself
    (``-c bcc``).
    '''
    if mol is None or mol.GetNumConformers() == 0:
        raise ValueError('mol_to_sdf requires a mol with a 3D conformer')
    try:
        return Chem.MolToMolBlock(mol, confId=conf_id, kekulize=True)
    except Exception:
        # A mol that only passed the relaxed sanitize may not kekulise; fall back
        # to aromatic bond types (ANTECHAMBER still perceives from the geometry).
        return Chem.MolToMolBlock(mol, confId=conf_id, kekulize=False)


# ---------------------------------------------------------------------------
# ChimeraX-side helpers (main thread only)
# ---------------------------------------------------------------------------

def _ccd_tables(session, resname):
    '''Return ``(chem_comp_atom, chem_comp_bond)`` CIFTables for a CCD id, or
    ``(None, None)``.

    Tries the template residue's metadata first (no network); ChimeraX's template
    cache does not always retain those tables, so it falls back to fetching the
    standalone CCD CIF (cached locally after the first fetch), which is the route
    ISOLDE's older ``ccd.py`` uses.'''
    name = resname.upper()
    from chimerax import mmcif
    try:
        tmpl = mmcif.find_template_residue(session, name)
        if tmpl is not None:
            at, bt = mmcif.get_mmcif_tables_from_metadata(
                tmpl, ['chem_comp_atom', 'chem_comp_bond'])
            if at:
                return at, bt
    except Exception:
        pass
    try:
        from chimerax.core.fetch import fetch_file
        url = 'https://files.rcsb.org/ligands/download/%s.cif' % name
        path = fetch_file(session, url, 'CCD %s' % name, '%s.cif' % name, 'CCD')
        at, bt = mmcif.get_mmcif_tables(path, ['chem_comp_atom', 'chem_comp_bond'])
        if at:
            return at, bt
    except Exception:
        pass
    return None, None


def template_geometry_mol(template):
    '''Best-effort RDKit molecule perceived from a ``TmplResidue``'s own atoms
    (which carry ideal coordinates and hydrogens), via :func:`mol_from_geometry`.

    A last-resort fallback when explicit CCD bond-order records can't be obtained.
    Returns ``None`` if any atom has an undefined coordinate (ChimeraX templates
    sometimes do -- see ISOLDE's ``copy_ideal_coords_to_exp``), since perception
    needs real geometry.'''
    import math
    atoms = list(template.atoms)
    pos = {a: i for i, a in enumerate(atoms)}
    symbols = [a.element.name for a in atoms]
    names = [a.name for a in atoms]
    coords = []
    for a in atoms:
        c = a.coord
        if c is None or any(math.isnan(float(v)) for v in c):
            return None
        coords.append((float(c[0]), float(c[1]), float(c[2])))
    bset = set()
    for a in atoms:
        for nb in a.neighbors:
            if nb in pos:
                bset.add(frozenset((pos[a], pos[nb])))
    bonds = [tuple(b) for b in bset]
    mol, _charge = mol_from_geometry(symbols, coords, bonds, names=names)
    return mol


def ccd_records(session, resname):
    '''Fetch the explicit chemistry of a CCD component as plain records.

    Reads the ``chem_comp_atom`` and ``chem_comp_bond`` mmCIF tables from the
    ChimeraX template residue's metadata (the template *object* carries no bond
    order, but its source CIF tables do).

    Returns:
        ``(atoms, bonds, ideal_coords)`` or ``None`` if the component is unknown:

        * ``atoms``  -- list of ``(atom_id, type_symbol, charge, aromatic_flag)``
        * ``bonds``  -- list of ``(atom_id_1, atom_id_2, value_order, arom_flag)``
        * ``ideal_coords`` -- ``{atom_id: (x, y, z)}`` from the ideal coordinates
          (falling back to model coordinates per missing atom).
    '''
    # Consult the local ChemComp store first (instant, offline). It returns the
    # same (atoms, bonds, ideal_coords) shape this function produces, so a store
    # hit is a drop-in for the network path below. Defensive: if the bundle is
    # absent or the store errors, fall through to the per-residue network fetch.
    try:
        from chimerax.chemcomp import lookup as _chemcomp_lookup
        rec = _chemcomp_lookup(session, resname)
        if rec is not None:
            return rec
    except Exception:
        pass
    atom_table, bond_table = _ccd_tables(session, resname)
    if not atom_table:
        return None
    atom_rows = atom_table.fields(
        ['atom_id', 'type_symbol', 'charge', 'pdbx_aromatic_flag',
         'pdbx_model_Cartn_x_ideal', 'pdbx_model_Cartn_y_ideal',
         'pdbx_model_Cartn_z_ideal', 'model_Cartn_x', 'model_Cartn_y',
         'model_Cartn_z'],
        allow_missing_fields=True, missing_value='?')
    atoms = []
    ideal_coords = {}
    for (aid, sym, charge, arom, xi, yi, zi, xm, ym, zm) in atom_rows:
        atoms.append((aid, sym, charge, arom))
        xyz = _coord_or_none(xi, yi, zi)
        if xyz is None:
            xyz = _coord_or_none(xm, ym, zm)
        if xyz is not None:
            ideal_coords[aid] = xyz
    bonds = []
    if bond_table is not None:
        bond_rows = bond_table.fields(
            ['atom_id_1', 'atom_id_2', 'value_order', 'pdbx_aromatic_flag'],
            allow_missing_fields=True, missing_value='N')
        bonds = [(a1, a2, order, arom) for (a1, a2, order, arom) in bond_rows]
    return atoms, bonds, ideal_coords


def _coord_or_none(x, y, z):
    try:
        return (float(x), float(y), float(z))
    except (TypeError, ValueError):
        return None


def perceive_residue(session, residue):
    '''Perceive an RDKit molecule for a modelled residue with no usable template.

    Clones the residue's heavy atoms into a throwaway (session-less)
    ``AtomicStructure``, completes the valences with ChimeraX's IDATM-based
    ``addh`` (which gives the added H real coordinates), estimates the net charge
    to seed the bond-order search, then lets RDKit perceive bond orders. Metals
    are detached, the organic framework perceived, and the metal re-attached with
    dative bonds.

    Returns:
        ``(mol, atom_map)`` where ``atom_map`` maps RDKit atom index -> ChimeraX
        ``Atom`` for the heavy atoms, or ``(None, {})`` if perception fails for
        every candidate charge.
    '''
    from chimerax.atomic import AtomicStructure, AtomicStructures, Atoms
    tmp = AtomicStructure(session, name='isolde-rdkit-perceive',
                          auto_style=False)
    try:
        r2 = tmp.new_residue(residue.name, 'A', 1)
        # Map throwaway atoms back to the live residue atoms by name.
        name_to_live = {}
        for a in residue.atoms:
            if a.element.number == 1:
                continue          # drop modelled H; addh re-adds them cleanly
            na = tmp.new_atom(a.name, a.element)
            na.coord = a.coord
            r2.add_atom(na)
            name_to_live[a.name] = a
        for b in residue.atoms.intra_bonds:
            a1, a2 = b.atoms
            if a1.element.number == 1 or a2.element.number == 1:
                continue
            na1 = r2.find_atom(a1.name)
            na2 = r2.find_atom(a2.name)
            if na1 is not None and na2 is not None:
                try:
                    tmp.new_bond(na1, na2)
                except Exception:
                    pass
        try:
            from chimerax.addh.cmd import cmd_addh
            cmd_addh(session, AtomicStructures([tmp]), hbond=False,
                     in_isolation=True)
        except Exception as e:
            session.logger.info('ISOLDE rdkit_bridge: addh failed (%s).' % e)

        # Strip coordinated metals so the organic framework perceives cleanly;
        # remember each metal's element + which non-metal atoms it coordinates.
        metal_atoms = [a for a in tmp.atoms if a.element.is_metal]
        metals = []
        for ma in metal_atoms:
            nbrs = [nb for nb in ma.neighbors
                    if not nb.element.is_metal and nb.element.number != 1]
            metals.append((ma.element.name, nbrs))   # nbrs survive deletion
        if metal_atoms:
            Atoms(metal_atoms).delete()

        # IDATM-based net-charge estimate (needs the added protons; computed on
        # the metal-free part) seeds rdDetermineBonds' charge search.
        est = None
        try:
            from chimerax.add_charge import estimate_net_charge
            est = estimate_net_charge(tmp.atoms)
        except Exception:
            est = None

        # Extract per-atom symbols/coords/names + index bonds.
        atoms = list(tmp.atoms)
        pos = {a: i for i, a in enumerate(atoms)}
        symbols = [a.element.name for a in atoms]
        names = [a.name for a in atoms]
        coords = [tuple(float(v) for v in a.coord) for a in atoms]
        bonds = [(pos[b.atoms[0]], pos[b.atoms[1]]) for b in tmp.bonds
                 if b.atoms[0] in pos and b.atoms[1] in pos]

        mol, _charge = mol_from_geometry(symbols, coords, bonds, names=names,
                                         preferred_charge=est)
        if mol is None:
            return None, {}
        if metals:
            recs = [(sym, [pos[nb] for nb in nbrs if nb in pos])
                    for sym, nbrs in metals]
            mol = attach_metals(mol, recs)
        atom_map = _build_atom_map(mol, name_to_live)
        return mol, atom_map
    finally:
        tmp.delete()


def _build_atom_map(mol, name_to_live):
    '''Map RDKit atom index -> ChimeraX ``Atom`` using the ``cxName`` property.'''
    atom_map = {}
    for a in mol.GetAtoms():
        if not a.HasProp(NAME_PROP):
            continue
        live = name_to_live.get(a.GetProp(NAME_PROP))
        if live is not None:
            atom_map[a.GetIdx()] = live
    return atom_map


def residue_signature(residue):
    '''A chemistry signature for lazy cache validation: SHA1 over sorted
    ``(heavy-atom name, idatm_type)``.

    Keyed on ``idatm_type`` (hybridisation/charge-state) rather than element, so
    it survives hydrogen add/remove and protonation tweaks but is invalidated by
    a true chemistry change (an in-place edit that alters the heavy-atom skeleton
    or its bonding: swapaa, Build Structure, ISOLDE replace-ligand).
    '''
    import hashlib
    h = hashlib.sha1()
    for n, t in sorted((a.name, a.idatm_type) for a in residue.atoms
                       if a.element.number != 1):
        h.update(('|%s:%s' % (n, t)).encode('utf-8'))
    return h.hexdigest()


def names_match(atoms, residue):
    '''True iff every modelled heavy atom is in the reference ``atoms`` records
    with the same element. Succeeds for partial residues (a subset of the
    reference atoms) but fails for a novel ligand that merely reused a CCD code
    (different atom names).'''
    ref = {}
    for (aid, el, charge, arom) in atoms:
        e = _norm_element(el)
        if e != 'H':
            ref[aid] = e
    for a in residue.atoms:
        if a.element.number == 1:
            continue
        if ref.get(a.name) != a.element.name:
            return False
    return True


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def cip_codes_by_name(mol):
    '''Return ``{cxName: CIP code}`` (e.g. 'R'/'S'/'r'/'s') for every atom that
    has an assigned CIP descriptor. Used to detect stereocentre inversion by
    comparing a modelled residue against its reference template.'''
    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception:
        pass
    out = {}
    for a in mol.GetAtoms():
        if a.HasProp(NAME_PROP) and a.HasProp('_CIPCode'):
            out[a.GetProp(NAME_PROP)] = a.GetProp('_CIPCode')
    return out


def spurious_stereocentre_substituents(mol):
    '''Identify CIP-perceived stereocentres that are actually NON-stereogenic.

    The canonical case is a phosphate/diphosphate phosphorus (or a sulfate/
    phosphonate S/P, etc.) whose two non-bridging oxygens are modelled in the CCD
    as one ``P=O`` and one ``P-OH`` / ``P-O(-)``. That bond-order/protonation
    asymmetry makes the two otherwise-identical oxygens CIP-distinct, so both the
    CCD (``pdbx_stereo_config``) and RDKit mark the centre chiral -- but the
    oxygens are resonance-equivalent, so it is not a real stereocentre. Restraining
    or flagging it is meaningless (and the restraint fights the freely-exchanging
    group in simulation).

    Returns ``{centre_cxName: [equivalent terminal cxNames]}`` for every such
    centre (empty if none). A centre is reported only when >= 2 of its neighbours
    are (a) **terminal** -- exactly one heavy neighbour, the centre itself -- and
    (b) the **same element**. That pair can differ *only* in bond order, formal
    charge or H count (there is nothing else attached to tell them apart), i.e.
    exactly the resonance/protonation degrees of freedom we treat as equivalent --
    so the test needs nothing more. It is also safe: a genuine stereocentre cannot
    have two equivalent substituents (that would give it a local mirror plane), so
    this never suppresses a real centre. Mixed-element centres (e.g. an O/S
    phosphorothioate) and bridging substituents are left untouched.'''
    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception:
        return {}

    def _heavy_deg(a):
        return sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() != 1)

    out = {}
    for centre in mol.GetAtoms():
        if not centre.HasProp('_CIPCode') or not centre.HasProp(NAME_PROP):
            continue
        by_element = {}
        for nb in centre.GetNeighbors():
            if nb.GetAtomicNum() == 1 or not nb.HasProp(NAME_PROP):
                continue
            if _heavy_deg(nb) != 1:
                continue  # bridging / substituted -- not an exchangeable terminus
            by_element.setdefault(nb.GetAtomicNum(), []).append(nb)
        for group in by_element.values():
            if len(group) >= 2:
                out[centre.GetProp(NAME_PROP)] = [a.GetProp(NAME_PROP) for a in group]
                break
    return out


def fmcs_index_correspondence(mol_a, mol_b, match_bond_order=True, timeout=10):
    '''Maximum common substructure atom correspondence between two molecules,
    element- and (optionally) bond-order-aware.

    Returns a list of ``(idx_a, idx_b)`` aligned atom-index pairs, or ``[]`` if no
    common substructure is found.

    Chirality is intentionally NOT enforced in the MCS itself: callers compare CIP
    codes on the matched pairs (see :func:`cip_codes_by_name`) and drop the
    inverted centres, so an inverted stereocentre falls *out* of the
    correspondence -- the signal to rebuild it -- rather than silently shrinking
    or blocking the whole match.
    '''
    from rdkit.Chem import rdFMCS
    res = rdFMCS.FindMCS(
        [mol_a, mol_b],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=(rdFMCS.BondCompare.CompareOrder if match_bond_order
                     else rdFMCS.BondCompare.CompareAny),
        ringMatchesRingOnly=True,
        completeRingsOnly=False,
        timeout=timeout)
    if res.numAtoms == 0 or not res.smartsString:
        return []
    patt = Chem.MolFromSmarts(res.smartsString)
    if patt is None:
        return []
    ma = mol_a.GetSubstructMatch(patt)
    mb = mol_b.GetSubstructMatch(patt)
    if not ma or not mb or len(ma) != len(mb):
        return []
    return list(zip(ma, mb))


def template_to_rdkit(session, resname):
    '''RDKit molecule for a CCD component (the full reference), with ``cxName``
    props and stereo from the ideal coordinates. Returns ``(mol, status)`` as
    :func:`mol_from_ccd`, or ``(None, 'failed')`` if the component is unknown.'''
    recs = ccd_records(session, resname)
    if recs is None:
        return None, 'failed'
    atoms, bonds, coords = recs
    return mol_from_ccd(atoms, bonds, name=resname, coords=coords)


def _resolve_base_templates(session, strings):
    '''Resolve user-supplied ``baseTemplates`` strings to reference RDKit mols with
    explicit chemistry, ONCE. Each string is tried as a CCD id / ChemComp-registered
    template first (:func:`template_to_rdkit`), then as a SMILES
    (``Chem.MolFromSmiles``). Unresolvable strings are warned about and skipped.
    Returns a list of reference mols (empty if ``strings`` is falsy).'''
    if not strings:
        return []
    mols = []
    for s in strings:
        s = str(s).strip()
        if not s:
            continue
        ref, kind = None, None
        try:
            m, _status = template_to_rdkit(session, s)
        except Exception:
            m = None
        if m is not None:
            ref, kind = m, 'CCD/registered template'
        else:
            try:
                sm = Chem.MolFromSmiles(s)
            except Exception:
                sm = None
            if sm is not None:
                ref, kind = sm, 'SMILES'
        if ref is None:
            session.logger.warning(
                'baseTemplate %r could not be resolved as a CCD id, a registered '
                'template, or a SMILES string; ignoring it.' % s)
            continue
        session.logger.info('baseTemplate %r resolved as %s (%d heavy atoms).'
            % (s, kind, sum(1 for a in ref.GetAtoms() if a.GetAtomicNum() != 1)))
        mols.append(ref)
    return mols


def apply_base_templates(frag_mol, ref_mols, min_match=3):
    '''Transfer trusted chemistry from reference exemplars onto a perceived fragment.

    For each reference mol, find a connectivity (bond-order-agnostic) MCS
    correspondence to ``frag_mol`` (:func:`fmcs_index_correspondence`) and copy the
    reference's **formal charges** (the load-bearing quantity -- it sets the net
    charge) plus bond orders and aromatic flags onto the matched fragment atoms/bonds.
    Matches smaller than ``min_match`` heavy atoms are ignored (spurious MCS), and
    atoms already set by an earlier (processed-first) exemplar are not overwritten, so
    several partial exemplars can each patch a region. Only bond types / formal
    charges / aromatic flags are mutated -- atoms are never rebuilt -- so the
    ``cxName``/``cxIdx`` identity props survive. Returns the (re-sanitised) mol.'''
    if not ref_mols:
        return frag_mol
    covered = set()
    changed = False
    for ref in ref_mols:
        try:
            pairs = fmcs_index_correspondence(frag_mol, ref, match_bond_order=False)
        except Exception:
            pairs = []
        heavy_pairs = [(fi, ri) for (fi, ri) in pairs
                       if frag_mol.GetAtomWithIdx(fi).GetAtomicNum() != 1]
        if len(heavy_pairs) < min_match:
            continue
        r2f = {ri: fi for (fi, ri) in pairs}
        fresh = [(fi, ri) for (fi, ri) in pairs if fi not in covered]
        if not fresh:
            continue
        for fi, ri in fresh:
            fa, ra = frag_mol.GetAtomWithIdx(fi), ref.GetAtomWithIdx(ri)
            fa.SetFormalCharge(ra.GetFormalCharge())
            fa.SetIsAromatic(ra.GetIsAromatic())
            covered.add(fi)
        for rb_bond in ref.GetBonds():
            fi1 = r2f.get(rb_bond.GetBeginAtomIdx())
            fi2 = r2f.get(rb_bond.GetEndAtomIdx())
            if fi1 is None or fi2 is None:
                continue
            fb = frag_mol.GetBondBetweenAtoms(fi1, fi2)
            if fb is None:
                continue
            fb.SetBondType(rb_bond.GetBondType())
            fb.SetIsAromatic(rb_bond.GetIsAromatic())
        changed = True
    if not changed:
        return frag_mol
    m2, _status = _sanitize_ladder(frag_mol, strip_hs=False)
    return m2 if m2 is not None else frag_mol


def residue_to_rdkit(residue, template=None):
    '''Return an RDKit molecule for ``residue`` with stereochemistry from its
    coordinates and a ``cxName`` property per atom.

    Regime selection (template match preferred, perception fallback):

    * If a CCD reference is available -- ``template`` given as a resname/template,
      or auto-resolved from ``residue.name`` -- and the modelled atom names are a
      subset of it (:func:`names_match`), build from the reference connection
      table (exact bond orders/charges) using the **modelled** coordinates, then
      fragment to the as-modelled atoms. No perception.
    * Otherwise perceive bond orders from geometry (:func:`perceive_residue`).

    Args:
        residue: a :class:`chimerax.atomic.Residue`.
        template: optional CCD resname (str) or template residue to match
            against; defaults to ``residue.name``.

    Returns:
        ``(mol, atom_map)`` where ``atom_map`` maps RDKit atom index -> ChimeraX
        ``Atom``, or ``(None, {})`` on failure.
    '''
    session = residue.structure.session
    resname = template if isinstance(template, str) else None
    if resname is None and template is not None:
        resname = getattr(template, 'name', None)
    if resname is None:
        resname = residue.name

    recs = ccd_records(session, resname)
    if recs is not None:
        atoms, bonds, ccd_coords = recs
        if names_match(atoms, residue):
            # Use modelled coordinates where available so stereo reflects the
            # model; fall back to the ideal coordinates for unmodelled atoms.
            coords = dict(ccd_coords)
            modelled = set()
            for a in residue.atoms:
                if a.element.number == 1:
                    continue
                coords[a.name] = tuple(float(v) for v in a.coord)
                modelled.add(a.name)
            mol, status = mol_from_ccd(atoms, bonds, name=resname,
                                       coords=coords)
            if mol is not None:
                # Restrict to the as-modelled heavy atoms (partial residues).
                if any(a.HasProp(NAME_PROP)
                       and a.GetProp(NAME_PROP) not in modelled
                       for a in mol.GetAtoms() if a.GetAtomicNum() != 1):
                    mol = fragment_by_names(mol, modelled)
                name_to_live = {a.name: a for a in residue.atoms
                                if a.element.number != 1}
                return mol, _build_atom_map(mol, name_to_live)

    # Fallback: perceive from geometry.
    return perceive_residue(session, residue)


def _unit_chemistry_maps(residue, base_template_mols=None):
    '''Per-residue chemistry for a super-residue as
    ``({atom_name: formal_charge}, {frozenset(name1, name2): rdkit BondType},
    {aromatic atom_name, ...})``.

    Sourced from :func:`residue_to_rdkit`, so it uses the CCD reference when one
    matches by name and **otherwise perceives from geometry** -- we cannot assume
    every ligand has a CCD exemplar. On total failure it returns empty maps and
    the caller defaults bonds to single / charges to zero.

    ``base_template_mols`` (resolved reference exemplars) override the CCD/perceived
    chemistry wherever they MCS-match the residue (:func:`apply_base_templates`) -- a
    standard residue simply won't match a ligand exemplar and is left untouched.'''
    charge_by_name = {}
    order_by_pair = {}
    aromatic = set()
    try:
        mol, _amap = residue_to_rdkit(residue)
    except Exception:
        mol = None
    if mol is None:
        return charge_by_name, order_by_pair, aromatic
    if base_template_mols:
        try:
            mol = apply_base_templates(mol, base_template_mols)
        except Exception:
            pass
    for a in mol.GetAtoms():
        if not a.HasProp(NAME_PROP):
            continue
        nm = a.GetProp(NAME_PROP)
        charge_by_name[nm] = a.GetFormalCharge()
        if a.GetIsAromatic():
            aromatic.add(nm)
    for b in mol.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.HasProp(NAME_PROP) and a2.HasProp(NAME_PROP):
            order_by_pair[frozenset((a1.GetProp(NAME_PROP),
                                     a2.GetProp(NAME_PROP)))] = b.GetBondType()
    return charge_by_name, order_by_pair, aromatic


def super_residue_to_rdkit(residues, exclude=None, neutralize_excluded_donors=False,
                           base_templates=None):
    '''Build ONE RDKit molecule spanning a covalent unit -- a ligand plus the
    standard residue(s) it is covalently bonded to -- with the wider chain
    replaced by caps, ready to feed to ANTECHAMBER / AM1-BCC.

    Every bond *among* the given residues is kept (including the covalent
    ligand<->partner link). Every bond *out* of the set (to the rest of the model)
    is replaced by a cap so the fragment has a sensible closed valence and a
    physical local geometry for the semi-empirical charge calculation:

    * a peptide backbone cut on a standard amino acid's ``N`` -> an **ACE**
      (acetyl) cap; on its ``C`` -> an **NME** (N-methyl) cap;
    * any other external bond -> a **methyl** cap.

    Cap heavy atoms are placed at the coordinates of the real neighbour atoms they
    stand in for (so the acetyl/N-methyl geometry comes straight from the chain).
    Hydrogens are added last with coordinates (``AddHs(addCoords=True)``).

    ``exclude`` (a set/iterable of ChimeraX ``Atom``\\ s) are treated as ABSENT: they
    are left out of the molecule entirely and a bond from a unit atom to an excluded
    atom is **not capped** (the unit atom keeps its modelled valence/protonation).
    This is how a coordinated metal is removed before ANTECHAMBER -- GAFF2 has no
    metal type, so the organic framework is typed metal-free and the metal is
    spliced back afterwards. Do NOT rely on this to trim ordinary chain atoms; those
    should be capped (pass them as unit boundaries, not exclusions).

    ``neutralize_excluded_donors`` (metal path only): a donor atom that the metal
    had deprotonated (bonded to an excluded atom AND carrying a negative CCD formal
    charge -- e.g. a chlorin pyrrole N or a Cys thiolate) is instead built NEUTRAL
    and PROTONATED, so AM1-BCC runs on a well-behaved neutral free base rather than a
    bare, uncompensated macrocyclic multi-anion (which makes sqm slow/unstable). The
    added H atoms are tagged ``CAP`` (dropped from the template like any cap), and
    ``info['neutralized']`` records ``(donor_rd_idx, [added_H_rd_idx...], orig_fc)``
    so the caller can fold the H charge back and re-impose ``orig_fc`` -- restoring
    the deprotonated coordinating state while conserving total charge.

    Heavy atoms of the unit carry a **globally unique** integer ``cxIdx`` (the key
    used to map back to ChimeraX -- atom *names* collide across residues) plus
    ``cxName``. Cap atoms and all hydrogens carry a ``CAP`` property and no
    ``cxIdx``, so they drop out of the round-trip automatically.

    Args:
        residues: an ordered iterable of :class:`chimerax.atomic.Residue` forming
            the covalent unit.

    Returns:
        ``(mol, cxidx_to_atom, info)`` where ``mol`` is an RDKit mol with
        hydrogens and a 3D conformer (``None`` on failure), ``cxidx_to_atom`` maps
        ``cxIdx`` -> ChimeraX ``Atom`` for the real unit atoms (heavy + modelled H;
        cap atoms and cap hydrogens are excluded), and ``info``
        is a dict with ``'status'`` ('ok'/'relaxed'/'failed'), ``'net_charge'``
        (integer total formal charge of the capped fragment) and ``'caps'`` (a
        list of cap descriptors).
    '''
    import numpy
    from chimerax.atomic import Residue, Atoms

    residues = list(residues)
    exclude = set(exclude) if exclude else set()
    # ALL atoms of the unit, heavy + modelled H. We keep the modelled hydrogens
    # (rather than dropping and re-perceiving them) so the emitted template's atom
    # graph matches the model residue EXACTLY -- OpenMM matches templates by element
    # + connectivity incl. H count, so a single mis-perceived bond order (common on
    # nitriles / hydrazones / heteroaromatics without a CCD entry) would otherwise
    # change the H count and make the template fail to match. Bond orders (which do
    # NOT affect matching) still come from the CCD/perception, for ANTECHAMBER only.
    unit_atoms = [a for r in residues for a in r.atoms if a not in exclude]
    unit_heavy = [a for a in unit_atoms if a.element.number != 1]
    unit_set = set(unit_atoms)
    if not unit_heavy:
        return None, {}, {'status': 'failed', 'net_charge': 0, 'caps': []}

    # Resolve any user-supplied chemistry exemplars ONCE (CCD id / registered
    # template / SMILES -> reference mols), then let them override the CCD/perceived
    # chemistry per residue wherever they match.
    base_template_mols = []
    if base_templates:
        base_template_mols = _resolve_base_templates(residues[0].structure.session,
                                                     base_templates)
    chem_maps = {r: _unit_chemistry_maps(r, base_template_mols) for r in residues}

    rw = Chem.RWMol()
    rd_of = {}                 # ChimeraX Atom -> rdkit index
    cxidx_to_atom = {}         # global cxIdx -> ChimeraX Atom
    coords = {}                # rdkit index -> (x, y, z)
    neutralized = []           # (ChimeraX donor atom, original negative fc)
    for gi, atom in enumerate(unit_atoms):
        a = Chem.Atom(_norm_element(atom.element.name))
        if atom.element.number != 1:
            charge_by_name, _order, aromatic = chem_maps[atom.residue]
            fc = charge_by_name.get(atom.name, 0)
            # Metal-deprotonated donor: build neutral + protonate (below) so sqm sees
            # a normal molecule; the deprotonation is restored by the caller.
            if (neutralize_excluded_donors and fc < 0
                    and any(nb in exclude for nb in atom.neighbors)):
                neutralized.append((atom, int(fc)))
                fc = 0
            if fc:
                a.SetFormalCharge(int(fc))
            if atom.name in aromatic:
                a.SetIsAromatic(True)
        a.SetIntProp(INDEX_PROP, gi)
        a.SetProp(NAME_PROP, atom.name)
        idx = rw.AddAtom(a)
        rd_of[atom] = idx
        cxidx_to_atom[gi] = atom
        coords[idx] = tuple(float(v) for v in atom.coord)

    # Bonds among unit atoms. Heavy-heavy intra-residue orders come from the
    # CCD/perceived chemistry; bonds to hydrogen and inter-residue links are single
    # (peptide / glycosidic / thioether links are all single bonds).
    for b in Atoms(unit_atoms).intra_bonds:
        a1, a2 = b.atoms
        i1, i2 = rd_of.get(a1), rd_of.get(a2)
        if i1 is None or i2 is None:
            continue
        if (a1.element.number != 1 and a2.element.number != 1
                and a1.residue is a2.residue):
            _c, order_by_pair, _ar = chem_maps[a1.residue]
            bond_type = order_by_pair.get(frozenset((a1.name, a2.name)),
                                          Chem.BondType.SINGLE)
        else:
            bond_type = Chem.BondType.SINGLE
        rw.AddBond(i1, i2, bond_type)
        if bond_type == Chem.BondType.AROMATIC:
            bd = rw.GetBondBetweenAtoms(i1, i2)
            bd.SetIsAromatic(True)
            rw.GetAtomWithIdx(i1).SetIsAromatic(True)
            rw.GetAtomWithIdx(i2).SetIsAromatic(True)

    caps_info = []
    cap_indices = []

    def _add_cap_atom(element, xyz):
        a = Chem.Atom(_norm_element(element))
        a.SetProp('CAP', '1')
        idx = rw.AddAtom(a)
        coords[idx] = tuple(float(v) for v in xyz)
        cap_indices.append(idx)
        return idx

    def _neighbour_coord(atom, name, fallback):
        for nb in atom.neighbors:
            if nb.name == name:
                return tuple(float(v) for v in nb.coord)
        return fallback

    # Cap every bond leaving the unit (heavy atoms only; H never bond out).
    for u in unit_heavy:
        ui = rd_of[u]
        for n in u.neighbors:
            if n.element.number == 1 or n in unit_set:
                continue
            if n in exclude:
                continue        # coordinated metal (spliced back later): leave open
            u_xyz = numpy.array(coords[ui])
            n_xyz = numpy.array([float(v) for v in n.coord])
            is_amino = (u.residue.polymer_type == Residue.PT_AMINO)
            if is_amino and u.name == 'N' and n.name == 'C':
                # N-side peptide cut -> ACE (acetyl): carbonyl C at prev-C, =O at
                # prev-O, methyl at prev-CA.
                cC = _add_cap_atom('C', n_xyz)
                o_xyz = _neighbour_coord(n, 'O', tuple(n_xyz + numpy.array([0.0, 1.2, 0.0])))
                ca_xyz = _neighbour_coord(n, 'CA', tuple(n_xyz + (n_xyz - u_xyz)))
                cO = _add_cap_atom('O', o_xyz)
                cM = _add_cap_atom('C', ca_xyz)
                rw.AddBond(ui, cC, Chem.BondType.SINGLE)
                rw.AddBond(cC, cO, Chem.BondType.DOUBLE)
                rw.AddBond(cC, cM, Chem.BondType.SINGLE)
                caps_info.append({'kind': 'ACE', 'on': (u.residue.name, u.name)})
            elif is_amino and u.name == 'C' and n.name == 'N':
                # C-side peptide cut -> NME (N-methyl): amide N at next-N, methyl
                # at next-CA.
                cN = _add_cap_atom('N', n_xyz)
                ca_xyz = _neighbour_coord(n, 'CA', tuple(n_xyz + (n_xyz - u_xyz)))
                cM = _add_cap_atom('C', ca_xyz)
                rw.AddBond(ui, cN, Chem.BondType.SINGLE)
                rw.AddBond(cN, cM, Chem.BondType.SINGLE)
                caps_info.append({'kind': 'NME', 'on': (u.residue.name, u.name)})
            else:
                # Generic cut -> methyl cap at the removed neighbour's position.
                cM = _add_cap_atom('C', n_xyz)
                rw.AddBond(ui, cM, Chem.BondType.SINGLE)
                caps_info.append({'kind': 'CH3', 'on': (u.residue.name, u.name)})

    # Attach a conformer, sanitise (perceives aromaticity from the Kekule orders),
    # then add hydrogens with coordinates.
    from rdkit.Geometry import Point3D
    m = rw.GetMol()
    conf = Chem.Conformer(m.GetNumAtoms())
    for idx, xyz in coords.items():
        conf.SetAtomPosition(idx, Point3D(*xyz))
    conf.Set3D(True)
    m.AddConformer(conf, assignId=True)

    m, status = _sanitize_ladder(m, strip_hs=False)
    if m is None:
        return None, {}, {'status': 'failed', 'net_charge': 0, 'caps': caps_info}
    # Hydrogen policy: if the model is protonated, keep its hydrogens exactly (add
    # H only to the caps) so the emitted template matches the model residue's graph.
    # Only if the model carries NO hydrogens at all do we protonate the whole thing
    # (ISOLDE requires H before simulating, but this keeps the all-heavy input from
    # producing a bare fragment for ANTECHAMBER).
    unit_has_h = any(a.element.number == 1 for a in unit_atoms)
    # Neutralised metal donors must also be protonated for the charge calc, even when
    # the rest of the model is already protonated (they were deprotonated in situ).
    neutral_rd = [rd_of[a] for (a, _fc) in neutralized]
    try:
        if not unit_has_h:
            m = Chem.AddHs(m, addCoords=True)
        else:
            only = list(cap_indices) + neutral_rd
            if only:
                m = Chem.AddHs(m, addCoords=True, onlyOnAtoms=only)
    except Exception:
        pass
    # Record the added donor-H atoms (tag CAP so they drop from the template) so the
    # caller can fold their charge back and re-impose the donor's original -ve charge.
    neutralized_info = []
    for (atom, fc) in neutralized:
        di = rd_of[atom]
        h_idxs = [nb.GetIdx() for nb in m.GetAtomWithIdx(di).GetNeighbors()
                  if nb.GetAtomicNum() == 1]
        for hi in h_idxs:
            m.GetAtomWithIdx(hi).SetProp('CAP', '1')
        neutralized_info.append((di, h_idxs, fc))
    net_charge = Chem.GetFormalCharge(m)
    return m, cxidx_to_atom, {'status': status, 'net_charge': net_charge,
                              'caps': caps_info, 'neutralized': neutralized_info}


# ---------------------------------------------------------------------------
# RDKit mol / residue -> ChemComp record dicts (ligand registration)
# ---------------------------------------------------------------------------


def _ccd_value_order(bond):
    '''CCD ``value_order`` string for an RDKit bond (default 'SING').'''
    bt = bond.GetBondType()
    return {
        Chem.BondType.SINGLE: 'SING',
        Chem.BondType.DOUBLE: 'DOUB',
        Chem.BondType.TRIPLE: 'TRIP',
        Chem.BondType.QUADRUPLE: 'QUAD',
        Chem.BondType.AROMATIC: 'AROM',
    }.get(bt, 'SING')


def _generated_atom_names(mol):
    '''Per-element sequential names (C1, C2, N1, O1, H1, ...) for a mol whose
    atoms lack a :data:`NAME_PROP` (e.g. one built from SMILES).'''
    counts = {}
    names = []
    for a in mol.GetAtoms():
        el = a.GetSymbol().upper()
        counts[el] = counts.get(el, 0) + 1
        names.append('%s%d' % (el, counts[el]))
    return names


def records_from_mol(mol, comp_id, display_name=None):
    '''Build a ChemComp record dict from an RDKit ``mol`` that carries a 3D
    conformer. Atom ids come from the ``cxName`` property when present, else are
    generated per-element. Returns the dict shape
    :func:`chimerax.chemcomp.register_record` consumes
    (``{'id', 'atoms'[6], 'bonds'[4], 'coords', 'display_name'?}``), or ``None``
    if the mol has no conformer.'''
    if mol is None or mol.GetNumConformers() == 0:
        return None
    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception:
        pass
    conf = mol.GetConformer()
    have_names = all(a.HasProp(NAME_PROP) for a in mol.GetAtoms())
    gen = None if have_names else _generated_atom_names(mol)
    idx_name = {}
    atoms = []
    coords = {}
    for i, a in enumerate(mol.GetAtoms()):
        name = a.GetProp(NAME_PROP) if a.HasProp(NAME_PROP) else gen[i]
        idx_name[a.GetIdx()] = name
        stereo = a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else ''
        atoms.append(
            [
                name,
                a.GetSymbol().upper(),
                str(a.GetFormalCharge()), 'Y' if a.GetIsAromatic() else 'N', stereo, ''
            ]
        )
        p = conf.GetAtomPosition(a.GetIdx())
        coords[name] = [float(p.x), float(p.y), float(p.z)]
    bonds = []
    for b in mol.GetBonds():
        bonds.append(
            [
                idx_name[b.GetBeginAtomIdx()], idx_name[b.GetEndAtomIdx()],
                _ccd_value_order(b), 'Y' if b.GetIsAromatic() else 'N'
            ]
        )
    rec = {'id': comp_id, 'atoms': atoms, 'bonds': bonds, 'coords': coords}
    if display_name:
        rec['display_name'] = display_name
    return rec


def records_from_residue(residue, comp_id=None, display_name=None):
    '''Build a ChemComp record dict from a modelled ``residue``, keeping **all**
    atoms (hydrogens included) with their model coordinates. Bond orders, formal
    charges, aromaticity and CIP stereo are perceived on the heavy-atom skeleton
    via :func:`residue_to_rdkit`; bonds to/within hydrogens default to single.
    Returns the record dict (see :func:`records_from_mol`).'''
    comp_id = comp_id or residue.name
    order = {}  # frozenset({name1, name2}) -> (value_order, arom_flag)
    stereo = {}  # name -> 'R'/'S'
    charges = {}  # name -> int formal charge
    arom_atom = {}  # name -> bool
    try:
        mol, _amap = residue_to_rdkit(residue)
    except Exception:
        mol = None
    if mol is not None:
        try:
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        except Exception:
            pass
        for a in mol.GetAtoms():
            if not a.HasProp(NAME_PROP):
                continue
            nm = a.GetProp(NAME_PROP)
            charges[nm] = a.GetFormalCharge()
            arom_atom[nm] = a.GetIsAromatic()
            if a.HasProp('_CIPCode'):
                stereo[nm] = a.GetProp('_CIPCode')
        for b in mol.GetBonds():
            a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
            if a1.HasProp(NAME_PROP) and a2.HasProp(NAME_PROP):
                key = frozenset((a1.GetProp(NAME_PROP), a2.GetProp(NAME_PROP)))
                order[key] = (_ccd_value_order(b), 'Y' if b.GetIsAromatic() else 'N')
    atoms = []
    coords = {}
    for a in residue.atoms:
        atoms.append(
            [
                a.name,
                a.element.name.upper(),
                str(charges.get(a.name, 0)), 'Y' if arom_atom.get(a.name) else 'N',
                stereo.get(a.name, ''), ''
            ]
        )
        coords[a.name] = [float(c) for c in a.coord]
    bonds = []
    for b in residue.atoms.intra_bonds:
        a1, a2 = b.atoms
        vo, ar = order.get(frozenset((a1.name, a2.name)), ('SING', 'N'))
        bonds.append([a1.name, a2.name, vo, ar])
    rec = {'id': comp_id, 'atoms': atoms, 'bonds': bonds, 'coords': coords}
    if display_name:
        rec['display_name'] = display_name
    return rec


def records_from_smiles(smiles, comp_id, display_name=None, seed=0xf00d):
    '''Build a ChemComp record dict from a SMILES string: parse, add hydrogens,
    embed a 3D conformer (ETKDG), assign stereo from the geometry, and emit the
    record (see :func:`records_from_mol`). Returns ``None`` if SMILES parsing or
    embedding fails.'''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    from rdkit.Chem import AllChem
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass
    try:
        Chem.AssignStereochemistryFrom3D(mol)
    except Exception:
        pass
    return records_from_mol(mol, comp_id, display_name=display_name)
