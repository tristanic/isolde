# @Author: Tristan Croll
# @Date:   13-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 13-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Covalent MD parameterisation for ISOLDE.

Some covalent chemistry cannot be parameterised residue-by-residue because the
bond joining the residues needs force-field parameters that span them, and
ANTECHAMBER only knows GAFF2. This covers three kinds of case, all handled the
same way:

* a covalent **ligand** bonded to a protein/nucleic-acid residue (e.g. a Cys
  warhead);
* a **residue-residue crosslink** between two otherwise-standard residues (e.g.
  the His-Tyr sidechain-sidechain bond of cytochrome c oxidase, /A:240,244 in
  1v54 -- and the covalently-linked *triplet* in some related species);
* an **unusual main-chain or sidechain modification** of a standard residue.

The approach (see the design plan):

1. gather the covalent *unit* -- the connected cluster of residues joined by
   non-standard covalent bonds -- and cap the surrounding chain
   (:func:`detect_covalent_unit` + ``rdkit_bridge.super_residue_to_rdkit``);
2. type/charge the capped super-residue with ANTECHAMBER/GAFF2 + AM1-BCC
   (``run_antechamber`` -- the single swappable typing/charge seam);
3. freeze a standard ff14SB core on the parts that are chemically unchanged and
   push the GAFF<->ff14SB seam down to the modified bonds, drawing the crossing
   parameters from a curated library (``boundary_params``) with a
   parmchk2-by-analogy fallback;
4. emit per-residue OpenMM templates with ``<ExternalBond>`` records and the seam
   parameters, and load them live.

This module owns steps 1-4; the RDKit chemistry lives in
``chimerax.isolde.atomic.rdkit_bridge``.

**What counts as a "standard backbone continuation"** (a cap boundary the unit is
NOT grown across): the canonical intra-polymer linkage between two standard
residues -- a peptide ``C-N`` bond between amino acids, or an ``O3'-P``
phosphodiester between nucleotides. *Every other* inter-residue covalent bond
(sidechain crosslinks, ligand attachments, isopeptide bonds, main-chain
modifications) is "non-standard": it is crossed, pulling the partner residue into
the unit, and later becomes a GAFF<->ff14SB seam. This is deliberately permissive
about main-chain modifications -- they are exactly the hard case worth supporting.
'''

#: Residue names treated as solvent (never grown into as unit members).
_SOLVENT = frozenset(('HOH', 'WAT', 'DOD', 'H2O'))

#: Atom-name pair of a peptide backbone bond (amino acid C -> next N).
_PEPTIDE_LINK = frozenset(('C', 'N'))

#: Atom-name pairs of a nucleotide phosphodiester backbone bond.
_NUCLEIC_LINKS = (frozenset(("O3'", 'P')), frozenset(('O3*', 'P')))


def _is_standard_polymer(residue):
    '''True for an amino acid or nucleotide -- a residue the base force field
    already parameterises as part of a polymer (a frozen-core candidate).'''
    from chimerax.atomic import Residue
    return residue.polymer_type in (Residue.PT_AMINO, Residue.PT_NUCLEIC)


def _is_backbone_continuation(a1, a2):
    '''True when the bond ``a1-a2`` is a standard intra-polymer backbone linkage
    between two standard residues (peptide ``C-N`` or phosphodiester ``O3'-P``).
    Such bonds are the boundaries the covalent unit is capped at, never crossed.'''
    from chimerax.atomic import Residue
    r1, r2 = a1.residue, a2.residue
    if not (_is_standard_polymer(r1) and _is_standard_polymer(r2)):
        return False
    names = frozenset((a1.name, a2.name))
    if r1.polymer_type == Residue.PT_AMINO and r2.polymer_type == Residue.PT_AMINO:
        return names == _PEPTIDE_LINK
    if r1.polymer_type == Residue.PT_NUCLEIC and r2.polymer_type == Residue.PT_NUCLEIC:
        return names in _NUCLEIC_LINKS
    return False


def _inter_residue_bonds(residue):
    '''Yield ``(this_atom, other_atom)`` for every covalent bond from ``residue``
    to a different residue.'''
    for a in residue.atoms:
        for nb in a.neighbors:
            if nb.residue is not residue:
                yield a, nb


class CovalentUnit:
    '''The residues + non-standard link bonds that make up one covalent unit.

    Attributes:
        residues: every ``Residue`` parameterised together (ligands and/or
            crosslinked/modified standard residues).
        links: list of ``(atom_a, atom_b)`` for the non-standard covalent bonds
            internal to the unit -- these become the GAFF<->ff14SB seams.
    '''

    def __init__(self, residues, links):
        self.residues = list(residues)
        self.links = list(links)

    @property
    def standard_residues(self):
        '''Unit residues the base force field parameterises as polymers -- the
        frozen-core candidates for Strategy A.'''
        return [r for r in self.residues if _is_standard_polymer(r)]

    @property
    def nonstandard_residues(self):
        '''Unit residues with no standard polymer parameterisation (ligands, novel
        components).'''
        return [r for r in self.residues
                if not _is_standard_polymer(r) and r.name not in _SOLVENT]

    @property
    def num_heavy_atoms(self):
        return sum(len([a for a in r.atoms if a.element.number != 1])
                   for r in self.residues)

    def __repr__(self):
        return '<CovalentUnit [%s] %d link(s)>' % (
            ', '.join('%s%d' % (r.name, r.number) for r in self.residues),
            len(self.links))


def detect_covalent_unit(residues, max_heavy_atoms=250):
    '''Resolve the covalent unit reachable from the selected residue(s).

    Starting from the selection, the unit is grown across every inter-residue
    covalent bond that is *not* a standard backbone continuation (see
    :func:`_is_backbone_continuation`) -- so ligand attachments, sidechain
    crosslinks and main-chain modifications all pull their partners in, while
    peptide/phosphodiester bonds bound the unit and are capped later. Solvent and
    coordinated (monatomic) metal ions are never grown into.

    Args:
        residues: a single ``Residue`` or an iterable of them (the selection).
        max_heavy_atoms: refuse units larger than this (AM1-BCC is unlikely to
            converge, and the fragment stops being a sensible local model).

    Raises ``chimerax.core.errors.UserError`` when there is no non-standard
    covalent linkage to parameterise (e.g. a free ligand -- use
    ``isolde parameterise`` -- or a plain residue), or the unit is too large.

    Returns:
        a :class:`CovalentUnit`.
    '''
    from chimerax.core.errors import UserError
    from chimerax.atomic import Residue

    if isinstance(residues, Residue):
        seeds = [residues]
    else:
        seeds = [r for r in residues]
    if not seeds:
        raise UserError('No residue selected for covalent parameterisation.')

    unit = list(seeds)
    unit_set = set(unit)
    queue = list(seeds)
    while queue:
        r = queue.pop()
        for a, nb in _inter_residue_bonds(r):
            other = nb.residue
            if other in unit_set:
                continue
            if _is_backbone_continuation(a, nb):
                continue                       # standard linkage -> cap boundary
            if other.name in _SOLVENT:
                continue                       # solvent, not a covalent member
            if other.num_atoms == 1 and other.atoms[0].element.is_metal:
                continue                       # coordinated metal ion, not covalent
            unit_set.add(other)
            unit.append(other)
            queue.append(other)

    # Non-standard covalent bonds internal to the unit (the seams).
    links = []
    seen = set()
    for r in unit:
        for a, nb in _inter_residue_bonds(r):
            if nb.residue in unit_set and not _is_backbone_continuation(a, nb):
                key = frozenset((a, nb))
                if key not in seen:
                    seen.add(key)
                    links.append((a, nb))

    if not links:
        raise UserError(
            'Selected residue(s) have no non-standard covalent linkage to '
            'parameterise. For a free ligand use "isolde parameterise" instead.')

    cu = CovalentUnit(unit, links)
    if cu.num_heavy_atoms > max_heavy_atoms:
        raise UserError(
            'The covalent unit %r has %d heavy atoms, above the limit of %d. '
            'AM1-BCC on a fragment this large is unlikely to converge; parameterise '
            'it externally.' % (cu, cu.num_heavy_atoms, max_heavy_atoms))
    return cu


# ---------------------------------------------------------------------------
# Typing / charge engine (the swappable seam)
# ---------------------------------------------------------------------------

def _parse_mol2_atoms(path):
    '''Read the ``@<TRIPOS>ATOM`` block of a MOL2 file. Returns three lists in
    file order: SYBYL/atom types, charges (float), and ``(x, y, z)`` coordinates.'''
    types, charges, coords = [], [], []
    in_atoms = False
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith('@<TRIPOS>'):
                in_atoms = (s == '@<TRIPOS>ATOM')
                continue
            if in_atoms and s:
                # id name x y z type subst_id subst_name charge
                p = s.split()
                if len(p) < 9:
                    continue
                coords.append((float(p[2]), float(p[3]), float(p[4])))
                types.append(p[5])
                charges.append(float(p[8]))
    return types, charges, coords


def run_antechamber(session, mol, net_charge, workdir=None, charge_method='bcc'):
    '''Type and charge ``mol`` with ANTECHAMBER (GAFF2 + AM1-BCC). This is the
    single **swappable typing/charge seam**: replace this one function to drive a
    different engine (the author's forthcoming whole-model force field, or a
    third-party per-molecule typer) without touching the rest of the pipeline.

    ``mol`` is a fully-built RDKit mol with a 3D conformer (from
    ``rdkit_bridge.super_residue_to_rdkit``). It is written to SDF, run through
    ``antechamber -fi sdf`` then ``parmchk2``, and the results are mapped back onto
    the RDKit atom indices **by coordinate** (robust even if ANTECHAMBER reorders).

    Args:
        net_charge: integer net charge of the whole capped fragment.
        workdir: scratch directory (a temp dir is made and removed if omitted).

    Returns a dict:
        ``types``   -- GAFF2 atom type per RDKit atom index;
        ``charges`` -- AM1-BCC charge per RDKit atom index;
        ``frcmod``  -- path to the parmchk2 frcmod (all-GAFF bonded parameters);
        ``mol2``    -- path to ANTECHAMBER's output mol2;
        ``workdir`` / ``cleanup`` -- scratch dir and whether the caller should remove it.

    Raises ``chimerax.core.errors.UserError`` on non-convergence / tool failure.
    '''
    from chimerax.core.errors import UserError
    # amber_home / amber_bin locate ChimeraX's *bundled* AMBER tools; ANTECHAMBER
    # shells out to sub-programs (bondtype, atomtype, sqm) via system() and needs
    # AMBERHOME set to find them -- exactly as chimerax.add_charge.charge does.
    from chimerax.amber_info import amber_bin, amber_home
    from chimerax.isolde.atomic import rdkit_bridge
    import os, subprocess, tempfile

    cleanup = False
    if workdir is None:
        workdir = tempfile.mkdtemp(prefix='isolde_covalent_')
        cleanup = True

    os.environ['AMBERHOME'] = amber_home

    sdf_path = os.path.join(workdir, 'input.sdf')
    with open(sdf_path, 'w') as f:
        f.write(rdkit_bridge.mol_to_sdf(mol))

    ante_out = os.path.join(workdir, 'ante.out.mol2')
    # Use the bundled binary by path (forward slashes for its cygwin runtime), and
    # -dr n to skip acdoctor (which rejects many valid but unusual fragments).
    # Do NOT pass -pf y: it makes ANTECHAMBER clean up via system("rm -f ...") which
    # has no `rm`/`/tmp` on Windows and fatals (matches chimerax.add_charge, which
    # omits -pf). We delete the whole temp workdir ourselves instead.
    ante_cmd = [amber_bin + '/antechamber', '-i', sdf_path, '-fi', 'sdf',
                '-o', ante_out, '-fo', 'mol2', '-c', charge_method,
                '-nc', str(int(net_charge)), '-at', 'gaff2',
                '-s', '2', '-dr', 'n']
    proc = subprocess.run(ante_cmd, cwd=workdir, capture_output=True, text=True)
    if proc.returncode != 0 or not os.path.isfile(ante_out):
        sqm_out = os.path.join(workdir, 'sqm.out')
        detail = (proc.stderr or proc.stdout or '')[-1500:]
        if os.path.isfile(sqm_out):
            with open(sqm_out) as f:
                detail += '\n--- sqm.out tail ---\n' + f.read()[-800:]
        raise UserError(
            'ANTECHAMBER failed on the covalent unit (net charge %d). This is often '
            'a wrong net-charge estimate or an AM1 convergence problem on a large or '
            'strained fragment. Output tail:\n%s' % (net_charge, detail))

    base_path = os.path.dirname(os.path.abspath(__file__))
    gaff2_parms = os.path.join(base_path, 'gaff2.dat')
    frcmod = os.path.join(workdir, 'unit.frcmod')
    # -a Y writes the molecule's COMPLETE parameter set (not just the derived/
    # missing terms), so a seam whose GAFF combo is a standard one still has a
    # value we can re-key onto the ff14SB boundary type.
    subprocess.run([amber_bin + '/parmchk2', '-s', '2', '-a', 'Y', '-i', ante_out,
                    '-f', 'mol2', '-p', gaff2_parms, '-o', frcmod],
                   cwd=workdir, capture_output=True, text=True)
    if not os.path.isfile(frcmod):
        raise UserError('parmchk2 did not produce a frcmod for the covalent unit.')

    types_f, charges_f, coords_f = _parse_mol2_atoms(ante_out)
    n = mol.GetNumAtoms()
    if len(types_f) != n:
        raise UserError('ANTECHAMBER returned %d atoms; expected %d. Atom count '
                        'mismatch -- cannot map parameters back.' % (len(types_f), n))

    # Map ANTECHAMBER output onto RDKit atom indices by nearest coordinate (robust
    # to any reordering; -c bcc does not move atoms).
    conf = mol.GetConformer()
    types = [None] * n
    charges = [0.0] * n
    used = [False] * n
    for k, (ox, oy, oz) in enumerate(coords_f):
        best_i, best_d = None, None
        for i in range(n):
            if used[i]:
                continue
            p = conf.GetAtomPosition(i)
            d = (p.x - ox) ** 2 + (p.y - oy) ** 2 + (p.z - oz) ** 2
            if best_d is None or d < best_d:
                best_d, best_i = d, i
        if best_i is None or best_d > 0.01:   # > 0.1 Angstrom -> no match
            raise UserError('Could not match an ANTECHAMBER output atom back to the '
                            'input by coordinate; aborting to avoid mis-assignment.')
        used[best_i] = True
        types[best_i] = types_f[k]
        charges[best_i] = charges_f[k]

    return {'types': types, 'charges': charges, 'frcmod': frcmod,
            'mol2': ante_out, 'workdir': workdir, 'cleanup': cleanup}


# ---------------------------------------------------------------------------
# Strategy A: per-atom typing + charge repartitioning
# ---------------------------------------------------------------------------

def _standard_template(forcefield, residue):
    '''The OpenMM residue template ISOLDE uses for a standard residue, or None.
    Tries the residue name directly (covers most amino acids / nucleotides); the
    caller treats a miss as "no frozen core" (everything GAFF-typed).'''
    templates = getattr(forcefield, '_templates', {})
    return templates.get(residue.name)


def _graph_distance_to_links(unit, shell_radius):
    '''BFS over intra-unit bonds from every link atom. Returns the set of atoms
    within ``shell_radius`` bonds of a seam (link) atom -- the "shell" that keeps
    its ff14SB type but takes an AM1-BCC charge.'''
    unit_atoms = set(a for r in unit.residues for a in r.atoms)
    seed = set()
    for a, b in unit.links:
        seed.add(a)
        seed.add(b)
    shell = set(seed)
    frontier = set(seed)
    for _ in range(shell_radius):
        nxt = set()
        for at in frontier:
            for nb in at.neighbors:
                if nb in unit_atoms and nb not in shell:
                    nxt.add(nb)
        shell |= nxt
        frontier = nxt
    return shell


def classify_atoms(unit, shell_radius=1):
    '''Tag every atom of the unit as ``'frozen'`` (unchanged standard atom -> keep
    ff14SB type + charge), ``'shell'`` (standard atom near a seam -> keep ff14SB
    type, take AM1-BCC charge), or ``'ligand'`` (non-standard residue atom -> GAFF
    type + AM1-BCC charge). Hydrogens inherit their heavy parent's tag.'''
    shell = _graph_distance_to_links(unit, shell_radius)
    standard = set(unit.standard_residues)
    tags = {}
    for r in unit.residues:
        for a in r.atoms:
            if a.element.number == 1:
                continue                    # handled after, from parent
            if r not in standard:
                tags[a] = 'ligand'
            elif a in shell:
                tags[a] = 'shell'
            else:
                tags[a] = 'frozen'
    # Hydrogens follow their (single) heavy neighbour.
    for r in unit.residues:
        for a in r.atoms:
            if a.element.number != 1:
                continue
            heavy = next((nb for nb in a.neighbors if nb.element.number != 1), None)
            tags[a] = tags.get(heavy, 'ligand')
    return tags


def _standard_template(forcefield, residue):
    '''The OpenMM residue template ISOLDE uses for ``residue``, or None.

    For a modified residue the modelled name (e.g. ``CYS``) is not always the right
    base: a Cys whose SG is bonded to a ligand is chemically a thioether, whose SG
    should be ff14SB type ``S`` (as in ``CYX``/``CYScyc``), not ``SH``. We reuse
    ISOLDE's own ``cys_type`` to pick that variant. Other residues fall back to the
    modelled name; a miss means "no frozen core" and the caller GAFF-types the whole
    residue (correct, but pushes the seam onto the backbone).'''
    templates = getattr(forcefield, '_templates', {})
    if residue.name == 'CYS':
        try:
            from ..openmm_interface import cys_type
            ct = cys_type(residue)
            if ct and ct in templates:
                return templates[ct]
        except Exception:
            pass
    return templates.get(residue.name)


def _template_atom_info(template):
    '''From an OpenMM residue template return ``(by_name, h_of_heavy)``:
    ``by_name`` = ``{atom_name: (amber_type, charge)}``; ``h_of_heavy`` =
    ``{heavy_name: (h_type, h_charge)}`` (the type/charge to give hydrogens on that
    heavy atom -- ff14SB H on a given heavy atom share these).'''
    by_name = {}
    is_h = {}
    for i, at in enumerate(template.atoms):
        by_name[at.name] = (at.type, at.parameters.get('charge', 0.0))
        is_h[i] = (at.element is not None and at.element.atomic_number == 1)
    h_of_heavy = {}
    for i, at in enumerate(template.atoms):
        if is_h[i]:
            continue
        for j in at.bondedTo:
            if is_h.get(j):
                hj = template.atoms[j]
                h_of_heavy.setdefault(at.name,
                                      (hj.type, hj.parameters.get('charge', 0.0)))
    return by_name, h_of_heavy


def _resolve_heavy(tag, info, name, ante_type, ante_charge):
    '''(type, charge) for a heavy atom under Strategy A.'''
    if tag == 'ligand' or info is None:
        return ante_type, ante_charge
    by_name = info[0]
    if name in by_name:
        ttype, tcharge = by_name[name]
        return (ttype, tcharge) if tag == 'frozen' else (ttype, ante_charge)
    return ante_type, ante_charge          # modified atom absent from the template


def _resolve_h(ptag, info, parent_name, ante_type, ante_charge):
    '''(type, charge) for a hydrogen, following its heavy parent's tag. Keeping the
    ff14SB H *type* on ff14SB-typed heavy atoms is essential -- a GAFF-typed H on an
    ff14SB carbon would create a spurious seam on every C-H bond.'''
    if ptag == 'ligand' or info is None:
        return ante_type, ante_charge
    h_of_heavy = info[1]
    if parent_name in h_of_heavy:
        htype, hcharge = h_of_heavy[parent_name]
        return (htype, hcharge) if ptag == 'frozen' else (htype, ante_charge)
    return ante_type, ante_charge


def _unique_name(seen, base):
    if base not in seen:
        return base
    i = 1
    while '%s%d' % (base, i) in seen:
        i += 1
    return '%s%d' % (base, i)


def assign_types_and_charges(unit, mol, cxidx_to_atom, ante, forcefield,
                             shell_radius=1):
    '''Build per-residue template records from the typed/charged super-residue.

    Applies Strategy A per atom (see :func:`classify_atoms`): frozen standard atoms
    keep their ff14SB type + charge; shell atoms keep the ff14SB type but take the
    AM1-BCC charge; ligand atoms take GAFF type + AM1-BCC charge. Hydrogens follow
    their heavy parent. Bonds and external-bond points are read from the RDKit
    graph: a bond to a cap or to another unit residue makes its endpoint an
    ``<ExternalBond>``.

    Returns ``{residue: {'atoms': [ {name, element, type, charge, tag, rd} ],
    'bonds': [(name1, name2)], 'external': set(names)}}``.
    '''
    from chimerax.isolde.atomic.rdkit_bridge import INDEX_PROP

    tags = classify_atoms(unit, shell_radius)
    tinfo = {}
    for r in unit.standard_residues:
        tmpl = _standard_template(forcefield, r)
        if tmpl is not None:
            tinfo[r] = _template_atom_info(tmpl)

    types, charges = ante['types'], ante['charges']
    rd2res, rd2name = {}, {}
    per_res = {}

    def rec(r):
        return per_res.setdefault(
            r, {'atoms': [], 'bonds': [], 'external': set(), '_names': set()})

    for a in mol.GetAtoms():
        j = a.GetIdx()
        if a.HasProp('CAP'):
            rd2res[j] = None
            continue
        if a.HasProp(INDEX_PROP):
            cx = cxidx_to_atom.get(a.GetIntProp(INDEX_PROP))
            if cx is None:
                rd2res[j] = None
                continue
            r = cx.residue
            tag = tags.get(cx, 'ligand')
            name = cx.name
            el = cx.element.name
            ttype, tcharge = _resolve_heavy(tag, tinfo.get(r), name,
                                            types[j], charges[j])
        else:                                   # hydrogen added by AddHs
            nbrs = a.GetNeighbors()
            hn = nbrs[0] if nbrs else None
            if hn is None or hn.HasProp('CAP') or not hn.HasProp(INDEX_PROP):
                rd2res[j] = None
                continue
            cx = cxidx_to_atom.get(hn.GetIntProp(INDEX_PROP))
            if cx is None:
                rd2res[j] = None
                continue
            r = cx.residue
            tag = tags.get(cx, 'ligand')
            el = 'H'
            name = _unique_name(rec(r)['_names'], 'H' + cx.name)
            ttype, tcharge = _resolve_h(tag, tinfo.get(r), cx.name,
                                        types[j], charges[j])
        d = rec(r)
        d['_names'].add(name)
        rd2res[j], rd2name[j] = r, name
        # 'type' is the final (ff14SB or GAFF) type; 'gaff_type' is what ANTECHAMBER
        # assigned -- kept so the emitter can look seam params up by GAFF analogy.
        d['atoms'].append({'name': name, 'element': el, 'type': ttype,
                           'gaff_type': types[j], 'charge': tcharge,
                           'tag': tag, 'rd': j})

    for b in mol.GetBonds():
        j, k = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        rj, rk = rd2res.get(j), rd2res.get(k)
        if rj is None and rk is None:
            continue
        if rj is None or rk is None:            # residue atom -> cap: external bond
            r, nm = (rk, rd2name.get(k)) if rj is None else (rj, rd2name.get(j))
            if nm:
                rec(r)['external'].add(nm)
            continue
        if rj is rk:
            rec(rj)['bonds'].append((rd2name[j], rd2name[k]))
        else:                                   # inter-residue link: external both
            rec(rj)['external'].add(rd2name[j])
            rec(rk)['external'].add(rd2name[k])

    for d in per_res.values():
        d.pop('_names', None)
    return per_res


def repartition_charges(per_res):
    '''Force each residue's template charges to sum to an integer, adjusting only
    non-frozen atoms so the frozen ff14SB core stays canonical (its external
    peptide bonds keep clean charge groups). The per-residue integer target is the
    nearest integer to the residue's current total (AM1-BCC over the whole unit
    already sums to the requested net charge, so per-residue rounding is small).

    Mutates ``per_res`` in place; returns ``{residue: integer_target}``.
    '''
    targets = {}
    for r, d in per_res.items():
        atoms = d['atoms']
        total = sum(a['charge'] for a in atoms)
        target = int(round(total))
        targets[r] = target
        adjustable = [a for a in atoms if a['tag'] != 'frozen']
        if not adjustable:
            continue                            # frozen-only: already ~integer
        share = (target - total) / len(adjustable)
        for a in adjustable:
            a['charge'] += share
    return targets


# ---------------------------------------------------------------------------
# Template naming + orchestration
# ---------------------------------------------------------------------------

def _make_template_names(unit, forcefield):
    '''A unique OpenMM template name per unit residue that never collides with a
    standard template (so a modified CYS never shadows plain ``CYS``). Applied to
    the specific residue instances via ``isolde_template_name``.'''
    used = set(getattr(forcefield, '_templates', {}).keys())
    names = {}
    for r in unit.residues:
        base = 'MCOV_' + r.name
        tn, i = base, 0
        while tn in used:
            i += 1
            tn = '%s_%d' % (base, i)
        names[r] = tn
        used.add(tn)
    return names


def parameterise_covalent_unit(session, unit, shell_radius=1, net_charge=None):
    '''Parameterise a covalent unit end to end and load it into ISOLDE's live
    forcefield.

    Pipeline: build the capped super-residue (``super_residue_to_rdkit``) -> type +
    charge it with ANTECHAMBER/AM1-BCC (``run_antechamber``) -> Strategy-A per-atom
    typing + integer charge repartition (``assign_types_and_charges`` /
    ``repartition_charges``) -> emit per-residue templates + seam parameters
    (``amber_convert.covalent_to_ffxml``) -> ``loadFile`` and set
    ``isolde_template_name`` on each unit residue so the new templates match only
    those instances.

    Args:
        unit: a :class:`CovalentUnit` (from :func:`detect_covalent_unit`).
        shell_radius: bonds from a seam within which standard atoms take AM1-BCC
            charges (Strategy A).
        net_charge: override the whole-unit net charge fed to AM1-BCC (defaults to
            the sum of formal charges from the RDKit build).

    Returns ``(xml_path, {residue: template_name})``. Raises
    ``chimerax.core.errors.UserError`` on failure.
    '''
    from chimerax.core.errors import UserError
    from chimerax.isolde.atomic.rdkit_bridge import super_residue_to_rdkit
    from .amber_convert import covalent_to_ffxml
    import os, shutil

    if not hasattr(session, 'isolde'):
        raise UserError('Start ISOLDE before parameterising a covalent unit.')

    mol, cxidx_to_atom, info = super_residue_to_rdkit(unit.residues)
    if mol is None:
        raise UserError('Could not build an RDKit model of the covalent unit '
                        '(status: %s).' % info.get('status'))
    nc = int(net_charge) if net_charge is not None else info['net_charge']

    ffmgr = session.isolde.forcefield_mgr
    forcefield = ffmgr[session.isolde.sim_params.forcefield]

    ante = run_antechamber(session, mol, nc)
    try:
        per_res = assign_types_and_charges(unit, mol, cxidx_to_atom, ante,
                                           forcefield, shell_radius=shell_radius)
        repartition_charges(per_res)
        template_names = _make_template_names(unit, forcefield)

        stem = '_'.join(sorted({r.name for r in unit.residues}))
        xml_path = os.path.join(os.getcwd(), '%s_covalent.xml' % stem)
        gaff2_xml = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gaff2.xml')
        covalent_to_ffxml(per_res, template_names, mol, ante['frcmod'], xml_path,
                          gaff2_xml)

        for tn in template_names.values():
            forcefield._templates.pop(tn, None)
        forcefield.loadFile(xml_path)
        for r, tn in template_names.items():
            r.isolde_template_name = tn
    finally:
        if ante.get('cleanup'):
            shutil.rmtree(ante['workdir'], ignore_errors=True)

    session.logger.info(
        'Parameterised covalent unit %r; wrote %s and loaded templates %s. '
        'These apply for this session; reload the ffXML in future sessions with '
        '"Load residue MD definition(s)".'
        % (unit, xml_path, ', '.join(template_names.values())))
    return xml_path, template_names


def parameterise_free_ligand(session, residue, net_charge=None):
    '''Parameterise a free (non-covalent) ligand through the same RDKit pipeline as
    the covalent path, replacing the classic ANTECHAMBER-via-add_charge route.

    The whole machinery reduces cleanly to the single-residue case: the "unit" is
    just the ligand (no partners, so ``super_residue_to_rdkit`` adds no caps), every
    atom is GAFF-typed, and one self-contained template is emitted with collision-
    proof **prefixed** atom types. Unlike a covalent unit -- whose modified residue
    shares a name with a standard one and so needs a per-instance
    ``isolde_template_name`` override -- a free ligand's name is its own, so the
    template is named ``USER_<name>`` and matched to **all copies** by name (the
    ``find_residue_templates`` USER_ pre-pass), exactly like the classic path.

    Returns ``(xml_path, template_name)``. Raises ``UserError`` on failure.
    '''
    from chimerax.core.errors import UserError
    from chimerax.isolde.atomic.rdkit_bridge import super_residue_to_rdkit
    from .amber_convert import covalent_to_ffxml
    import os, shutil

    if not hasattr(session, 'isolde'):
        raise UserError('Start ISOLDE before parameterising a ligand.')

    unit = CovalentUnit([residue], [])          # no links -> no caps
    mol, cxidx_to_atom, info = super_residue_to_rdkit(unit.residues)
    if mol is None:
        raise UserError('Could not build an RDKit model of %s (status: %s).'
                        % (residue.name, info.get('status')))
    nc = int(net_charge) if net_charge is not None else info['net_charge']

    forcefield = session.isolde.forcefield_mgr[session.isolde.sim_params.forcefield]
    tname = 'USER_' + residue.name

    ante = run_antechamber(session, mol, nc)
    try:
        per_res = assign_types_and_charges(unit, mol, cxidx_to_atom, ante, forcefield)
        repartition_charges(per_res)
        template_names = {residue: tname}
        xml_path = os.path.join(os.getcwd(), residue.name + '.xml')
        gaff2_xml = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gaff2.xml')
        covalent_to_ffxml(per_res, template_names, mol, ante['frcmod'], xml_path,
                          gaff2_xml)
        forcefield._templates.pop(tname, None)
        forcefield.loadFile(xml_path)
    finally:
        if ante.get('cleanup'):
            shutil.rmtree(ante['workdir'], ignore_errors=True)

    session.logger.info(
        'Parameterised ligand %s as template %s (wrote %s). Applies to all copies '
        'this session; reload the ffXML in future sessions with "Load residue MD '
        'definition(s)".' % (residue.name, tname, xml_path))
    return xml_path, tname

