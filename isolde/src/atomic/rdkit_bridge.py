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
