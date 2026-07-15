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


#: Above this heavy-atom count, AM1-BCC charges are computed by divide-and-conquer
#: (fragment -> parallel sqm -> stitch) instead of one whole-molecule sqm, whose
#: ~O(N^3) cost otherwise makes large ligands (e.g. a phytyl-tailed chlorophyll)
#: take tens of minutes. Also the target size for alkyl fragments. Kept fairly low:
#: sqm is cubic, so a 30-atom fragment is ~8x cheaper than a 60-atom one, and the
#: conjugated-core fragment (which has no safe internal cut) is unaffected by this
#: cap either way. See ``_bcc_charges_fragmented``.
MAX_FRAGMENT_HEAVY = 30


def _amber_bin():
    '''Locate ChimeraX's bundled AMBER tools and set AMBERHOME (ANTECHAMBER shells
    out to bondtype/atomtype/sqm via system() and needs it, as chimerax.add_charge
    does). Returns the bin directory.'''
    from chimerax.amber_info import amber_bin, amber_home
    import os
    os.environ['AMBERHOME'] = amber_home
    return amber_bin


def _mmff_minimized_copy(mol):
    '''A geometry-relaxed COPY of ``mol`` (MMFF, with UFF fallback for any element),
    so single-point AM1 charges aren't corrupted by clashy/strained input geometry --
    which is common for freshly-built covalent links and just-placed ligands. Atom
    order + identity props are preserved (a copy), and the caller's model geometry is
    untouched. Returns the original ``mol`` on any failure (charges then fall back to
    the input geometry).'''
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:
        return mol
    if mol.GetNumConformers() == 0:
        return mol
    m = Chem.Mol(mol)
    try:
        props = AllChem.MMFFGetMoleculeProperties(m)
        if props is not None and AllChem.MMFFOptimizeMolecule(m, maxIters=1000) != -1:
            return m
    except Exception:
        pass
    try:
        AllChem.UFFOptimizeMolecule(m, maxIters=1000)
        return m
    except Exception:
        return mol


def _antechamber_on_mol(amber_bin, mol, net_charge, workdir, charge_method, label):
    '''Run ANTECHAMBER once on ``mol`` and map the result back onto its atom indices
    by coordinate. Returns ``(types, charges, ante_out_path)``. ``charge_method`` is
    ``'bcc'`` (AM1-BCC, invokes sqm) or ``'gas'`` (Gasteiger, instant -- used only
    for typing). ``label`` names the fragment/unit in error messages.'''
    from chimerax.core.errors import UserError
    from chimerax.isolde.atomic import rdkit_bridge
    import os, subprocess
    # For AM1-BCC, relax the geometry first (cheap MMFF/UFF) so the single-point sqm
    # charges are not corrupted by strained input (e.g. a freshly-built covalent link
    # or just-placed ligand); the coordinate map-back below uses this same conformer.
    work_mol = _mmff_minimized_copy(mol) if charge_method == 'bcc' else mol
    sdf_path = os.path.join(workdir, 'input.sdf')
    with open(sdf_path, 'w') as f:
        f.write(rdkit_bridge.mol_to_sdf(work_mol))
    ante_out = os.path.join(workdir, 'ante.out.mol2')
    # Bundled binary by path (forward slashes for its cygwin runtime); -dr n skips
    # acdoctor (rejects many valid unusual fragments); NO -pf y (its system("rm ...")
    # cleanup fatals on Windows -- we delete the temp workdir ourselves).
    ante_cmd = [amber_bin + '/antechamber', '-i', sdf_path, '-fi', 'sdf',
                '-o', ante_out, '-fo', 'mol2', '-c', charge_method,
                '-nc', str(int(net_charge)), '-at', 'gaff2', '-s', '2', '-dr', 'n']
    if charge_method == 'bcc':
        # sqm geometry-optimises by default (dozens-hundreds of O(N^3) SCF steps),
        # which utterly dominates the cost on a large/conjugated fragment -- a
        # bacteriochlorin core took ~1000 s optimised vs ~1 s here. maxcyc=0 does a
        # single-point AM1 at the INPUT (modelled) geometry, which is both far faster
        # and more appropriate: AM1-BCC charges are bond-local and geometry-
        # insensitive, and we want them at the conformation being fitted, not at a
        # gas-phase AM1 minimum. (Not needed for -c gas, which never runs sqm.)
        ante_cmd += ['-ek', 'maxcyc=0']
    proc = subprocess.run(ante_cmd, cwd=workdir, capture_output=True, text=True)
    if proc.returncode != 0 or not os.path.isfile(ante_out):
        sqm_out = os.path.join(workdir, 'sqm.out')
        detail = (proc.stderr or proc.stdout or '')[-1500:]
        if os.path.isfile(sqm_out):
            with open(sqm_out) as f:
                detail += '\n--- sqm.out tail ---\n' + f.read()[-800:]
        raise UserError(
            'ANTECHAMBER failed on %s (net charge %d). This is often a wrong '
            'net-charge estimate or an AM1 convergence problem on a strained '
            'fragment. Output tail:\n%s' % (label, net_charge, detail))

    types_f, charges_f, coords_f = _parse_mol2_atoms(ante_out)
    n = work_mol.GetNumAtoms()
    if len(types_f) != n:
        raise UserError('ANTECHAMBER returned %d atoms for %s; expected %d. Atom '
                        'count mismatch -- cannot map parameters back.'
                        % (len(types_f), label, n))
    conf = work_mol.GetConformer()
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
            raise UserError('Could not match an ANTECHAMBER output atom of %s back '
                            'to the input by coordinate.' % label)
        used[best_i] = True
        types[best_i] = types_f[k]
        charges[best_i] = charges_f[k]
    return types, charges, ante_out


def _parmchk2(amber_bin, ante_out, workdir):
    '''Run parmchk2 (-a Y = complete parameter set) on an ANTECHAMBER mol2 -> frcmod
    path. Charges are irrelevant here, so this works off a gas-typed mol2 too.'''
    from chimerax.core.errors import UserError
    import os, subprocess
    gaff2_parms = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gaff2.dat')
    frcmod = os.path.join(workdir, 'unit.frcmod')
    subprocess.run([amber_bin + '/parmchk2', '-s', '2', '-a', 'Y', '-i', ante_out,
                    '-f', 'mol2', '-p', gaff2_parms, '-o', frcmod],
                   cwd=workdir, capture_output=True, text=True)
    if not os.path.isfile(frcmod):
        raise UserError('parmchk2 did not produce a frcmod.')
    return frcmod


def _bcc_charges_fragmented(session, amber_bin, mol, partition, base_workdir):
    '''AM1-BCC charges for a large ``mol`` by divide-and-conquer: build a capped
    fragment per ``partition`` group (``rdkit_bridge.build_capped_fragment``), run
    ``antechamber -c bcc`` on each IN PARALLEL (independent small sqm jobs), then
    stitch charges back onto the parent atoms via the ``fragSrc`` prop, discarding
    caps and renormalising each fragment to its integer net charge. Returns a
    per-parent-atom charge list. Cuts sit only in transferable alkyl regions
    (``safe_cut_bonds``), so the stitched charges closely match a whole-molecule
    fit.'''
    from chimerax.core.errors import UserError
    from chimerax.isolde.atomic import rdkit_bridge
    from chimerax.isolde.atomic.rdkit_bridge import FRAGSRC_PROP
    from concurrent.futures import ThreadPoolExecutor
    import os

    jobs = []
    for gi, atom_set in enumerate(partition):
        fmol = rdkit_bridge.build_capped_fragment(mol, atom_set)
        if fmol is None:
            raise UserError('Could not build charge fragment %d of the ligand.' % gi)
        fnc = sum(mol.GetAtomWithIdx(i).GetFormalCharge() for i in atom_set)
        jobs.append((gi, fmol, int(fnc)))

    def _one(job):
        gi, fmol, fnc = job
        wd = os.path.join(base_workdir, 'frag%d' % gi)
        os.makedirs(wd, exist_ok=True)
        _t, fcharges, _o = _antechamber_on_mol(amber_bin, fmol, fnc, wd, 'bcc',
                                               'charge fragment %d' % gi)
        return gi, fmol, fnc, fcharges

    workers = max(1, min(len(jobs), (os.cpu_count() or 2)))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        results = list(ex.map(_one, jobs))

    charges = [0.0] * mol.GetNumAtoms()
    for gi, fmol, fnc, fcharges in results:
        reals, total = [], 0.0
        for a in fmol.GetAtoms():
            if a.HasProp(FRAGSRC_PROP):
                pi = a.GetIntProp(FRAGSRC_PROP)
                reals.append((pi, fcharges[a.GetIdx()]))
                total += fcharges[a.GetIdx()]
        # Renormalise the real atoms to the fragment's integer net charge (the caps,
        # now discarded, carried a small amount). Tiny, since cuts are in ~0-charge
        # alkyl regions.
        delta = (fnc - total) / len(reals) if reals else 0.0
        for pi, q in reals:
            charges[pi] = q + delta

    # Float cleanup so the whole-molecule total is exactly integer (per-fragment
    # renorm already makes it sum to sum(fnc) == the mol's formal charge).
    net = sum(charges)
    resid = (round(net) - net) / len(charges) if charges else 0.0
    return [c + resid for c in charges]


def run_antechamber(session, mol, net_charge, workdir=None, charge_method='bcc'):
    '''Type and charge ``mol`` with ANTECHAMBER (GAFF2 + AM1-BCC). This is the
    single **swappable typing/charge seam**: replace this one function to drive a
    different engine (the author's forthcoming whole-model force field, or a
    third-party per-molecule typer) without touching the rest of the pipeline.

    ``mol`` is a fully-built RDKit mol with a 3D conformer (from
    ``rdkit_bridge.super_residue_to_rdkit``). Typing (GAFF2) and the frcmod come from
    a whole-molecule ANTECHAMBER pass; AM1-BCC charges are the expensive part (sqm is
    ~O(N^3)), so for a large ligand (> ``MAX_FRAGMENT_HEAVY`` heavy atoms with safe
    alkyl cut points) they are computed by **divide-and-conquer** -- fragment at
    transferable bonds, charge the small fragments in parallel, stitch back
    (``_bcc_charges_fragmented``) -- while typing/frcmod use an instant Gasteiger
    (``-c gas``) whole-molecule pass (atom types are charge-method-independent). A
    small ligand takes the original single ``-c bcc`` pass unchanged.

    Args:
        net_charge: integer net charge of the whole (capped) fragment.
        workdir: scratch directory (a temp dir is made and removed if omitted).

    Returns a dict:
        ``types``   -- GAFF2 atom type per RDKit atom index;
        ``charges`` -- AM1-BCC charge per RDKit atom index;
        ``frcmod``  -- path to the parmchk2 frcmod (all-GAFF bonded parameters);
        ``mol2``    -- path to ANTECHAMBER's output mol2;
        ``workdir`` / ``cleanup`` -- scratch dir and whether the caller should remove it.

    Raises ``chimerax.core.errors.UserError`` on non-convergence / tool failure.
    '''
    from chimerax.isolde.atomic import rdkit_bridge
    import os, tempfile

    cleanup = False
    if workdir is None:
        workdir = tempfile.mkdtemp(prefix='isolde_covalent_')
        cleanup = True
    amber_bin = _amber_bin()

    heavy = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() != 1)
    partition = None
    if charge_method == 'bcc' and heavy > MAX_FRAGMENT_HEAVY:
        try:
            partition = rdkit_bridge.charge_partition(mol, MAX_FRAGMENT_HEAVY)
        except Exception:
            partition = None
        if partition is not None and len(partition) <= 1:
            partition = None       # nothing safe to cut -> whole molecule

    if partition is None:
        if charge_method == 'bcc' and heavy > MAX_FRAGMENT_HEAVY:
            session.logger.info(
                'Ligand has %d heavy atoms and no safe alkyl cut points; running '
                'AM1-BCC on the whole molecule (this may be slow).' % heavy)
        types, charges, ante_out = _antechamber_on_mol(
            amber_bin, mol, net_charge, workdir, charge_method, 'the ligand unit')
        frcmod = _parmchk2(amber_bin, ante_out, workdir)
        return {'types': types, 'charges': charges, 'frcmod': frcmod,
                'mol2': ante_out, 'workdir': workdir, 'cleanup': cleanup}

    # Large ligand: instant Gasteiger typing (whole) + fragmented AM1-BCC charges.
    session.logger.info(
        'Large ligand (%d heavy atoms): AM1-BCC charges computed on %d fragments in '
        'parallel (divide-and-conquer) to avoid a slow whole-molecule sqm.'
        % (heavy, len(partition)))
    types, _gas, ante_out = _antechamber_on_mol(
        amber_bin, mol, net_charge, workdir, 'gas', 'the ligand unit (typing)')
    frcmod = _parmchk2(amber_bin, ante_out, workdir)
    charges = _bcc_charges_fragmented(session, amber_bin, mol, partition, workdir)
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

def _run_antechamber_autocharge(session, mol, atoms, perceived_net, label,
                                override=None):
    '''run_antechamber with automatic net-charge selection. Neither charge source is
    universally right -- the RDKit-perceived/CCD value is wrong when the modelled
    protonation differs from the CCD neutral form (deprotonated acid), while
    ChimeraX's ``estimate_net_charge`` is wrong for some exotic ligands -- so we let
    ANTECHAMBER be the arbiter: it fails loudly (odd electron count / SCF
    non-convergence) on a wrong total charge. Try the perceived value first, then the
    H-aware estimate, skipping any candidate whose electron count is odd (guaranteed
    to fail). ``override`` (the user's ``netCharge``) short-circuits everything.

    Returns the run_antechamber result dict. Raises ``UserError`` naming the tried
    candidates if all fail (pointing at the ``netCharge`` override).'''
    from chimerax.core.errors import UserError
    import tempfile, shutil

    if override is not None:
        return run_antechamber(session, mol, int(override))

    est = None
    try:
        from chimerax.add_charge import estimate_net_charge
        atoms.idatm_types                       # force IDATM recompute (rings etc.)
        est = int(estimate_net_charge(atoms))
    except Exception:
        est = None

    total_electrons = sum(a.GetAtomicNum() for a in mol.GetAtoms())
    candidates = []
    for c in (perceived_net, est):              # perceived first, estimate fallback
        if c is None or c in candidates:
            continue
        if (total_electrons - c) % 2 != 0:
            continue                            # odd electrons -> would fail; skip
        candidates.append(c)
    if not candidates:                          # parity check ruled all out: try anyway
        candidates = [c for c in (perceived_net, est) if c is not None]

    last = None
    for i, c in enumerate(candidates):
        wd = tempfile.mkdtemp(prefix='isolde_charge_')
        try:
            ante = run_antechamber(session, mol, c, workdir=wd)
            ante['cleanup'] = True
            if i > 0:
                session.logger.info(
                    'Net charge %s for %s did not parameterise; used %+d instead.'
                    % (candidates[0], label, c))
            return ante
        except Exception as e:
            last = e
            shutil.rmtree(wd, ignore_errors=True)
    raise UserError(
        'Could not parameterise %s at any automatic net-charge estimate (%s); '
        'ANTECHAMBER failed for all. Set the correct net charge explicitly with the '
        'netCharge option. Last error:\n%s'
        % (label, candidates, str(last)[:400]))


def _make_template_names(unit, forcefield, prefix='MCOV_'):
    '''A unique OpenMM template name per unit residue that never collides with a
    standard template (so a modified CYS never shadows plain ``CYS``). Applied to
    the specific residue instances via ``isolde_template_name``.'''
    used = set(getattr(forcefield, '_templates', {}).keys())
    names = {}
    for r in unit.residues:
        base = prefix + r.name
        tn, i = base, 0
        while tn in used:
            i += 1
            tn = '%s_%d' % (base, i)
        names[r] = tn
        used.add(tn)
    return names


def parameterise_covalent_unit(session, unit, shell_radius=1, net_charge=None,
                               base_templates=None):
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

    mol, cxidx_to_atom, info = super_residue_to_rdkit(unit.residues,
                                                      base_templates=base_templates)
    if mol is None:
        raise UserError('Could not build an RDKit model of the covalent unit '
                        '(status: %s).' % info.get('status'))
    from chimerax.atomic import Atoms
    unit_atoms = Atoms([a for r in unit.residues for a in r.atoms])

    ffmgr = session.isolde.forcefield_mgr
    forcefield = ffmgr[session.isolde.sim_params.forcefield]

    ante = _run_antechamber_autocharge(session, mol, unit_atoms, info['net_charge'],
                                       repr(unit), override=net_charge)
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


def parameterise_free_ligand(session, residue, net_charge=None, base_templates=None):
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
    mol, cxidx_to_atom, info = super_residue_to_rdkit(unit.residues,
                                                      base_templates=base_templates)
    if mol is None:
        raise UserError('Could not build an RDKit model of %s (status: %s).'
                        % (residue.name, info.get('status')))
    forcefield = session.isolde.forcefield_mgr[session.isolde.sim_params.forcefield]
    tname = 'USER_' + residue.name

    ante = _run_antechamber_autocharge(session, mol, residue.atoms,
                                       info['net_charge'], residue.name,
                                       override=net_charge)
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


# ---------------------------------------------------------------------------
# Metal sites (hemes, chlorophylls, Zn fingers, ...)
# ---------------------------------------------------------------------------
#
# GAFF2/ANTECHAMBER cannot type a transition metal, so a metal site is
# parameterised as a BONDED template (like the shipped BRYCE_HEM / iron_sulfur.xml
# hand-curated ones, but generated automatically). The organic framework is typed
# metal-free by the covalent pipeline above; the metal is spliced back in as an ion
# (Lennard-Jones from the bundled ion set, charge partly spread onto its donors)
# with soft, empirical metal-donor bonds and coordination angles from
# ``metal_params``. Making the coordination a real force-field bond -- rather than a
# distance restraint -- is essential: OpenMM then drops the metal<->donor 1-2 and
# donor<->donor 1-3 nonbonded exclusions, so a modest bond constant holds the
# geometry instead of fighting the full Lennard-Jones wall a nonbonded ion would
# present at a ~2 A coordination distance.

#: Heavy elements accepted as metal-coordinating donor atoms.
_DONOR_ELEMENTS = frozenset(('N', 'O', 'S', 'Se', 'P', 'F', 'Cl', 'Br', 'I'))

#: Per-donor-element search radius (Angstrom) for geometric coordination detection.
_DONOR_CUTOFF = {'N': 2.6, 'O': 2.6, 'S': 3.0, 'Se': 3.1, 'P': 3.0,
                 'F': 2.4, 'Cl': 3.0, 'Br': 3.1, 'I': 3.3}
_DONOR_CUTOFF_DEFAULT = 2.8


class MetalSite(CovalentUnit):
    '''A metal coordination site: the metal atom(s), their donor atoms, and the
    residues parameterised together (the metalloligand and any coordinating
    residues). Subclasses :class:`CovalentUnit` so the RDKit build, capping and
    Strategy-A typing all apply unchanged. ``links`` is EMPTY: metal coordination is
    not a covalent seam (it is handled by the empirical metal-donor terms), so a
    coordinating standard residue has no seam bond and is fully frozen to ff14SB --
    and in fact excluded from the ANTECHAMBER fragment entirely, since once the metal
    is removed it is disconnected from the metalloligand (see
    :func:`parameterise_metal_site`).'''

    def __init__(self, residues, metals, coordination):
        super().__init__(residues, [])
        self.metals = list(metals)
        self.coordination = list(coordination)   # list of (metal_atom, donor_atom)

    @property
    def num_heavy_atoms(self):
        return sum(len([a for a in r.atoms if a.element.number != 1])
                   for r in self.residues)

    def __repr__(self):
        return '<MetalSite [%s] metals=%s %d bond(s)>' % (
            ', '.join('%s%d' % (r.name, r.number) for r in self.residues),
            '/'.join(sorted({m.element.name for m in self.metals})),
            len(self.coordination))


def _find_donors(metal):
    '''Coordinating donor atoms of ``metal``: its real heavy donor-element
    neighbours, plus any donor-element heavy atom within the per-element cutoff
    (geometric detection, so a site whose coordination bonds were never drawn is
    still found). Solvent donors are skipped (coordinating waters are not modelled
    as bonded in this first version -- they simply drift as ordinary TIP3P).'''
    import numpy
    structure = metal.structure
    mcoord = numpy.array([float(v) for v in metal.coord])
    donors = []
    seen = set()

    def _consider(a):
        if a is metal or a in seen:
            return
        if a.element.number == 1 or a.element.is_metal:
            return
        if a.element.name not in _DONOR_ELEMENTS:
            return
        if a.residue.name in _SOLVENT:
            return
        seen.add(a)
        donors.append(a)

    for nb in metal.neighbors:
        if not nb.element.is_metal:
            _consider(nb)
    atoms = structure.atoms
    coords = atoms.coords
    d2 = ((coords - mcoord) ** 2).sum(axis=1)
    max_cut = max(max(_DONOR_CUTOFF.values()), _DONOR_CUTOFF_DEFAULT)
    within = numpy.where(d2 <= max_cut ** 2)[0]
    for i in within:
        a = atoms[i]
        if a.element.name not in _DONOR_ELEMENTS:
            continue
        cut = _DONOR_CUTOFF.get(a.element.name, _DONOR_CUTOFF_DEFAULT)
        if d2[i] <= cut ** 2:
            _consider(a)
    return donors


def detect_metal_site(residues, max_heavy_atoms=250):
    '''Resolve the metal coordination site reachable from the selected residue(s).

    Collects the metal atom(s) in (or coordinated by) the selection, finds each
    metal's donor atoms (:func:`_find_donors`), and gathers the metalloligand and
    coordinating residues into a :class:`MetalSite`. The wider chain is capped at
    peptide/phosphodiester bonds exactly as for a covalent unit.

    Raises ``chimerax.core.errors.UserError`` when the selection contains no metal
    (use ``isolde parameterise``), when no donors can be found, or when the site is
    too large for AM1-BCC on its organic framework.
    '''
    from chimerax.core.errors import UserError
    from chimerax.atomic import Residue

    if isinstance(residues, Residue):
        seeds = [residues]
    else:
        seeds = list(residues)
    if not seeds:
        raise UserError('No residue selected for metal-site parameterisation.')

    metals = []
    seen_m = set()

    def _add_metal(a):
        if a not in seen_m and a.element.is_metal and a.residue.name not in _SOLVENT:
            seen_m.add(a)
            metals.append(a)

    for r in seeds:
        for a in r.atoms:
            if a.element.is_metal:
                _add_metal(a)                       # metal in the selection
            for nb in a.neighbors:
                if nb.element.is_metal:
                    _add_metal(nb)                  # metal bonded to the selection

    # Geometric: a metal within coordination distance of any selected atom, so that
    # selecting a coordinating residue resolves the site even if the coordination
    # bonds were never drawn in the model.
    import numpy
    structure = seeds[0].structure
    seed_coords = numpy.array([[float(v) for v in a.coord]
                               for r in seeds for a in r.atoms])
    if len(seed_coords):
        max_cut = max(max(_DONOR_CUTOFF.values()), _DONOR_CUTOFF_DEFAULT)
        for a in structure.atoms:
            if not a.element.is_metal or a.residue.name in _SOLVENT or a in seen_m:
                continue
            ac = numpy.array([float(v) for v in a.coord])
            if ((seed_coords - ac) ** 2).sum(axis=1).min() <= max_cut ** 2:
                _add_metal(a)

    if not metals:
        raise UserError(
            'Selected residue(s) contain no metal atom. For an ordinary organic '
            'ligand use "isolde parameterise".')

    coordination = []
    donor_residues = set()
    for m in metals:
        donors = _find_donors(m)
        if not donors:
            continue
        for d in donors:
            coordination.append((m, d))
            donor_residues.add(d.residue)
    if not coordination:
        raise UserError(
            'No coordinating donor atoms (N/O/S/...) were found around the metal '
            'atom(s) %s. Check the coordination geometry, or select the '
            'coordinating residues explicitly.'
            % '/'.join(sorted({m.element.name for m in metals})))

    residue_set = {m.residue for m in metals} | donor_residues
    site = MetalSite(residue_set, metals, coordination)
    if site.num_heavy_atoms > max_heavy_atoms:
        raise UserError(
            'The metal site %r has %d heavy atoms, above the limit of %d. AM1-BCC '
            'on its organic framework is unlikely to converge; parameterise it '
            'externally.' % (site, site.num_heavy_atoms, max_heavy_atoms))
    return site


def _angle_radians(a, b, c):
    '''Angle a-b-c (b is the vertex) in radians, from ChimeraX atom coordinates.'''
    import numpy, math
    va = numpy.array([float(v) for v in a.coord]) - numpy.array([float(v) for v in b.coord])
    vc = numpy.array([float(v) for v in c.coord]) - numpy.array([float(v) for v in b.coord])
    na, nc = numpy.linalg.norm(va), numpy.linalg.norm(vc)
    if na == 0 or nc == 0:
        return 0.0
    cosv = float(numpy.dot(va, vc) / (na * nc))
    cosv = max(-1.0, min(1.0, cosv))
    return math.acos(cosv)


def _distance_nm(a, b):
    '''Distance (nm) between two ChimeraX atoms, from their modelled coordinates.'''
    import numpy
    d = numpy.array([float(v) for v in a.coord]) - numpy.array([float(v) for v in b.coord])
    return float(numpy.linalg.norm(d)) * 0.1


def _cluster_core_atoms(site):
    '''Bridging inorganic *core* atoms of a polynuclear cluster: non-metal heavy atoms
    of the metalloligand whose heavy neighbours are ALL metals -- e.g. the mu-sulfides
    of an iron-sulfur cluster (each bonded only to Fe). These are NOT organic (a Cys SG
    or a heme pyrrole N has a carbon neighbour, so it fails this test and stays an
    external/organic donor), so they cannot be typed by GAFF/AM1-BCC; they are emitted
    as explicit, observed-geometry, uniquely-typed cluster atoms instead. Requires the
    metal-donor coordination bonds to already be real Bonds (parameterise_metal_site
    creates them first).'''
    nonstandard = set(site.nonstandard_residues)
    core = []
    for r in nonstandard:
        for a in r.atoms:
            if a.element.number == 1 or a.element.is_metal:
                continue
            heavy_nb = [nb for nb in a.neighbors if nb.element.number != 1]
            if heavy_nb and all(nb.element.is_metal for nb in heavy_nb):
                core.append(a)
    return core


def _ion_template_for_element(element, ox_state):
    '''The bundled ion TEMPLATE name (via ``metal_name_map``) for a metal of the
    given element/oxidation state, trying the oxidation-tagged CCD name first
    (e.g. ``FE2``) then the bare element. Returns ``None`` if unparameterised.'''
    from .metal_name_map import metal_name_map
    for key in ('%s%d' % (element.upper(), ox_state), element.upper()):
        if key in metal_name_map:
            tn = metal_name_map[key]
            if tn is not None:
                return tn
    return None


#: Bundled ffXML files that carry the monatomic-ion Lennard-Jones parameters.
_ION_LJ_FILES = ('tip3p_HFE_multivalent.xml', 'tip3p_standard.xml')


def _ion_lj(forcefield, ion_template_name):
    '''``(sigma, epsilon, mass, element_symbol)`` for the atom type of a bundled ion
    template. The atom TYPE name is read from the loaded template; the Lennard-Jones
    sigma/epsilon are read straight from the shipped ion ffXML ``<NonbondedForce>``
    (robust across OpenMM versions, and guaranteed consistent with the ion set
    ISOLDE actually simulates). Returns ``None`` if the template/type is absent.'''
    import os
    import xml.etree.ElementTree as ET
    tmpl = getattr(forcefield, '_templates', {}).get(ion_template_name)
    if tmpl is None or not tmpl.atoms:
        return None
    at = tmpl.atoms[0]
    atype = at.type
    element = getattr(at.element, 'symbol', None)
    mass = getattr(getattr(at.element, 'mass', None), '_value', None)
    mass_s = '%.4f' % mass if mass is not None else '0.0'
    base = os.path.dirname(os.path.abspath(__file__))
    for fn in _ION_LJ_FILES:
        path = os.path.join(base, fn)
        if not os.path.isfile(path):
            continue
        nb = ET.parse(path).getroot().find('NonbondedForce')
        if nb is None:
            continue
        for a in nb.findall('Atom'):
            if (a.get('type') or a.get('class')) == atype:
                return (a.get('sigma'), a.get('epsilon'), mass_s, element or 'X')
    return None


def _build_metal_terms(session, site, forcefield, template_names, keep_fraction,
                       core_atoms=None):
    '''Build the explicit metal-coordination emission for :func:`covalent_to_ffxml`
    (see its ``metal_terms`` argument): the metal atoms + their ion Lennard-Jones,
    the metal-donor bonds and donor-metal-donor angles (soft empirical values from
    ``metal_params``, angles snapped to the nearest ideal), and the per-donor charge
    deltas that spread the metal's oxidation-state charge onto its (ligand) donors.

    Donors in a **standard** coordinating residue (His, Asp, ...) are NOT part of the
    emitted templates -- those residues keep their normal ISOLDE templates -- so their
    coordination terms are keyed by the donor's ff14SB atom type
    (:func:`metal_params.coordinating_donor_type`; a coordinating His N is always type
    ``NB``), reported in ``terms['donor_types']``, and they take no charge delta.

    ``core_atoms`` (a polynuclear cluster's bridging sulfides/oxides, from
    :func:`_cluster_core_atoms`) switches on the CLUSTER path: every metal and core atom
    is uniquely typed, the metal-core bonds / all cluster angles (incl. metal-core-metal)
    take the observed geometry, stiff ``AmoebaUreyBradleyForce`` 1-3 terms
    (``terms['urey_bradley']``) hold the cage, and the core atoms are emitted explicitly
    (``terms['core']``: name/etype/charge) with :data:`metal_params.CORE_ATOM_LJ`. Empty
    (mononuclear) -> the original curated, angle-snapped, soft-bond behaviour.'''
    from chimerax.core.errors import UserError
    from itertools import combinations
    from .metal_params import (metal_bond_params, snap_angle, angle_k,
                               metal_charge_split, guess_ox_state,
                               coordinating_donor_type, CORE_ATOM_LJ,
                               UREY_BRADLEY_K)

    nonstandard = set(site.nonstandard_residues)
    core_set = set(core_atoms or [])
    # A polynuclear inorganic cluster (bridging core atoms present) is held by explicit
    # observed-geometry terms + stiff Urey-Bradley 1-3 distances, with every metal/core
    # atom UNIQUELY typed so no two terms compete. A mononuclear site keeps the curated,
    # angle-snapped, physically-soft-bond behaviour (heme/Zn/Mg) unchanged.
    is_cluster = len(core_set) > 0
    terms = {'metals': [], 'core': [], 'atom_types': [], 'lj': [], 'coord_bonds': [],
             'coord_angles': [], 'donor_deltas': [], 'donor_types': {},
             'urey_bradley': []}

    # Per-atom-unique emitted type for a cluster atom (metal / core). Atom names are
    # unique within a residue, so this never collides -- and for a mononuclear metal
    # named by its element it is identical to the old element-keyed type.
    def _etype(atom):
        return template_names[atom.residue] + '_' + atom.name

    # Core atoms are emitted here (not via the organic per_res); seed each with its
    # formal anion charge (bridging sulfide/oxide = -2) and add the metal-charge spread.
    core_charge = {c: (-2.0 if c.element.name in ('S', 'O') else 0.0) for c in core_set}

    for m in site.metals:
        r = m.residue
        elem = m.element.name
        # ChimeraX Atoms carry no reliable formal charge, so the oxidation state
        # comes from metal_params' per-element default (overridable there).
        ox = guess_ox_state(elem, None)
        donors = [d for (mm, d) in site.coordination if mm is m]
        # Only ligand (non-standard) donors -- organic OR core -- carry a redistributed
        # charge; standard coordinating residues keep their untouched ISOLDE templates.
        lig_donors = [d for d in donors if d.residue in nonstandard]
        q_metal, delta = metal_charge_split(ox, len(lig_donors), keep_fraction)
        etype = _etype(m)
        ion_tmpl = _ion_template_for_element(elem, ox)
        lj = _ion_lj(forcefield, ion_tmpl) if ion_tmpl else None
        if lj is None:
            raise UserError(
                'No bundled Lennard-Jones parameters for metal %s (oxidation state '
                '~%d). It is not in ISOLDE\'s ion set, so a metal-site template '
                'cannot be built; parameterise this site externally.' % (elem, ox))
        sigma, epsilon, mass, esym = lj
        terms['metals'].append({'residue': r, 'name': m.name, 'etype': etype,
                                'charge': q_metal})
        terms['atom_types'].append({'name': etype, 'element': esym, 'mass': mass})
        terms['lj'].append({'type': etype, 'sigma': sigma, 'epsilon': epsilon})
        for d in donors:
            is_core = d in core_set
            # A donor from a separate standard residue (the heme-His case) is an
            # axial/dative bond -- physically weaker than an in-ligand bond; held by
            # the angles, not a stiff bond. No-op for pairs without an axial entry.
            axial = d.residue not in nonstandard
            p = metal_bond_params(elem, d.element.name, axial=axial)
            if p is None:
                raise UserError(
                    'No metal-donor distance for %s-%s; cannot build a coordination '
                    'bond. Add it to metal_params.METAL_SITE_PARAMS.'
                    % (elem, d.element.name))
            r0, k = p
            # Cluster-internal (metal-core) bonds freeze the observed distance so the
            # rigid cluster keeps its own geometry; external/mononuclear use curated r0.
            length = _distance_nm(m, d) if is_core else r0
            terms['coord_bonds'].append({'metal': (r, m.name),
                                         'donor': (d.residue, d.name),
                                         'length': length, 'k': k})
            if d.residue in nonstandard:
                if is_core:
                    core_charge[d] += delta         # folded into the emitted core atom
                else:
                    terms['donor_deltas'].append({'residue': d.residue, 'name': d.name,
                                                  'delta': delta})
            else:
                dt = coordinating_donor_type(d.residue.name, d.name, d.element.name)
                if dt is None:
                    raise UserError(
                        'No ff14SB type known for coordinating donor %s/%s; cannot '
                        'build its metal bond.' % (d.residue.name, d.name))
                terms['donor_types'][(d.residue, d.name)] = dt
        # Angles at the metal vertex, over its donor pairs. Cluster: observed (cluster
        # geometry is non-canonical); mononuclear: snapped to the ideal polyhedron.
        for d1, d2 in combinations(donors, 2):
            th0 = _angle_radians(d1, m, d2) if is_cluster \
                else snap_angle(_angle_radians(d1, m, d2))
            terms['coord_angles'].append({'a': (d1.residue, d1.name),
                                          'metal': (r, m.name),
                                          'c': (d2.residue, d2.name),
                                          'angle': th0, 'k': angle_k()})
            # Urey-Bradley across the metal for a core...core pair (the S...S term).
            # Keyed by the full angle triple (i-vertex-k), as OpenMM's UB force expects.
            if d1 in core_set and d2 in core_set:
                terms['urey_bradley'].append((_etype(d1), _etype(m), _etype(d2),
                                              _distance_nm(d1, d2), UREY_BRADLEY_K))

    # Core-vertex angles (metal-core-metal, e.g. Fe-S-Fe) + the metal...metal Urey-
    # Bradley across each bridge. Observed geometry; the vertex is the core atom.
    for c in core_set:
        etype = _etype(c)
        terms['core'].append({'residue': c.residue, 'name': c.name, 'etype': etype,
                              'charge': core_charge[c]})
        lj = CORE_ATOM_LJ.get(c.element.name)
        if lj is None:
            raise UserError(
                'No bundled Lennard-Jones parameters for bridging core element %s; add '
                'it to metal_params.CORE_ATOM_LJ.' % c.element.name)
        terms['atom_types'].append({'name': etype, 'element': c.element.name,
                                    'mass': '%.4f' % c.element.mass})
        terms['lj'].append({'type': etype, 'sigma': '%.12f' % lj[0],
                            'epsilon': '%.7f' % lj[1]})
        metal_nbrs = [m for m in site.metals if m in c.neighbors]
        for m1, m2 in combinations(metal_nbrs, 2):
            terms['coord_angles'].append({'a': (m1.residue, m1.name),
                                          'metal': (c.residue, c.name),   # vertex = core
                                          'c': (m2.residue, m2.name),
                                          'angle': _angle_radians(m1, c, m2),
                                          'k': angle_k()})
            terms['urey_bradley'].append((_etype(m1), _etype(c), _etype(m2),
                                          _distance_nm(m1, m2), UREY_BRADLEY_K))

    # Clean up floating-point drift so the emitted cluster carries an exact integer net
    # charge (the split conserves the formal metal-oxidation + core-anion total).
    if is_cluster and terms['metals']:
        emitted = sum(x['charge'] for x in terms['metals']) \
            + sum(x['charge'] for x in terms['core'])
        terms['metals'][0]['charge'] += round(emitted) - emitted
    return terms


def parameterise_metal_site(session, site, shell_radius=1, net_charge=None,
                            metal_keep_fraction=0.5, base_templates=None):
    '''Parameterise a metal coordination site end to end and load it into ISOLDE's
    live force field.

    Pipeline: build the capped super-residue with the metal(s) EXCLUDED
    (``super_residue_to_rdkit(..., exclude=metals)``) -> type + charge the organic
    framework with ANTECHAMBER/AM1-BCC -> Strategy-A per-atom typing
    (``assign_types_and_charges``) -> splice the metal back as a bonded ion with soft
    empirical coordination terms + redistributed charge (``_build_metal_terms``) ->
    emit per-residue templates + coordination parameters
    (``amber_convert.covalent_to_ffxml`` with ``metal_terms``) -> ``loadFile`` and
    set ``isolde_template_name`` on every site residue.

    Returns ``(xml_path, {residue: template_name})``. Raises
    ``chimerax.core.errors.UserError`` on failure.
    '''
    from chimerax.core.errors import UserError
    from chimerax.isolde.atomic.rdkit_bridge import super_residue_to_rdkit
    from .amber_convert import covalent_to_ffxml
    import os, shutil

    if not hasattr(session, 'isolde'):
        raise UserError('Start ISOLDE before parameterising a metal site.')

    # The coordination must be REAL bonds in the model: create_openmm_topology only
    # bonds atoms.intra_bonds, and it is those topology bonds that make OpenMM drop
    # the metal<->donor 1-2 / donor<->donor 1-3 nonbonded exclusions the bonded model
    # relies on. Create any coordination bond that is not already present.
    added = 0
    for (m, d) in site.coordination:
        if d not in m.neighbors:
            try:
                m.structure.new_bond(m, d)
                added += 1
            except Exception:
                pass
    if added:
        session.logger.info('Added %d metal-coordination bond(s) to the model so '
                            'the coordination is parameterised as bonded.' % added)

    from chimerax.atomic import Atoms
    metals = set(site.metals)
    # Bridging inorganic core atoms (iron-sulfur mu-sulfides, ...) are NOT organic: they
    # bond only to metals, so RDKit would fill their empty valence with spurious H and
    # AM1-BCC would type them nonsensically. Exclude them from the organic build alongside
    # the metals and emit them as explicit, uniquely-typed cluster atoms (see the CLUSTER
    # path in _build_metal_terms). For a PURE cluster (F3S) nothing organic remains, so
    # super_residue_to_rdkit returns mol=None and the AM1-BCC step is skipped entirely.
    core_atoms = _cluster_core_atoms(site)
    excluded = metals | set(core_atoms)
    ffmgr = session.isolde.forcefield_mgr
    forcefield = ffmgr[session.isolde.sim_params.forcefield]

    # Charge ONLY the metalloligand's organic framework. The coordinating standard
    # residues (His etc.) reach it only through the Fe-N bonds we exclude, so once the
    # metal is gone they are disconnected -- they keep their normal ISOLDE templates
    # and are left out of the sqm run entirely (their coordination is emitted as
    # type-keyed terms in _build_metal_terms). A bare-ion site (metal coordinated only
    # by protein, no organic metalloligand) charges nothing -- mol is None.
    ligand_residues = list(site.nonstandard_residues)
    ligand_unit = CovalentUnit(ligand_residues, [])
    # neutralize_excluded_donors: build metal-deprotonated donors (chlorin pyrrole N,
    # thiolate) as neutral + protonated so AM1-BCC runs on a well-behaved neutral free
    # base; the deprotonation is restored below, charge-conserving.
    mol, cxidx_to_atom, info = super_residue_to_rdkit(
        ligand_residues, exclude=excluded, neutralize_excluded_donors=True,
        base_templates=base_templates)

    per_res, ante, frcmod = {}, None, None
    if mol is not None:
        # Net charge is unreliable from the CCD field, so let ANTECHAMBER arbitrate.
        lig_atoms = Atoms([a for r in ligand_residues for a in r.atoms
                           if not a.element.is_metal])
        ante = _run_antechamber_autocharge(session, mol, lig_atoms, info['net_charge'],
                                           repr(site), override=net_charge)
    try:
        if ante is not None:
            per_res = assign_types_and_charges(ligand_unit, mol, cxidx_to_atom, ante,
                                               forcefield, shell_radius=shell_radius)
            frcmod = ante['frcmod']
            # Restore the metal-deprotonated donors: fold each surrogate H's AM1-BCC
            # charge into its donor and re-impose the donor's original negative formal
            # charge (conserves total charge; charges were computed on the easy neutral
            # species and the H atoms are dropped from the template).
            rd_to_spec = {s['rd']: s for d in per_res.values() for s in d['atoms']}
            for (donor_rd, h_rds, orig_fc) in info.get('neutralized', []):
                spec = rd_to_spec.get(donor_rd)
                if spec is None:
                    continue
                spec['charge'] += sum(ante['charges'][h] for h in h_rds) + orig_fc
        # Templates are emitted only for the metalloligand / metal-ion residues; the
        # coordinating standard residues keep their normal ISOLDE templates. The name
        # is deterministic `MMET_<resname>` (not a per-instance override) so a
        # name-based pass in find_residue_templates binds it to ALL copies -- a
        # metalloligand type is the same across instances (cf. the free-ligand USER_).
        template_names = {r: 'MMET_' + r.name for r in ligand_residues}
        metal_terms = _build_metal_terms(session, site, forcefield, template_names,
                                         metal_keep_fraction, core_atoms=core_atoms)

        stem = '_'.join(sorted({r.name for r in ligand_residues})) or 'metal'
        xml_path = os.path.join(os.getcwd(), '%s_metal.xml' % stem)
        gaff2_xml = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gaff2.xml')
        covalent_to_ffxml(per_res, template_names, mol, frcmod, xml_path,
                          gaff2_xml, metal_terms=metal_terms)

        for tn in template_names.values():
            forcefield._templates.pop(tn, None)
        forcefield.loadFile(xml_path)
    finally:
        if ante is not None and ante.get('cleanup'):
            shutil.rmtree(ante['workdir'], ignore_errors=True)

    # Fan out to EVERY other copy of the metalloligand in the model: give each real
    # coordination bonds so the name-matched template binds it too (create_openmm_
    # topology only bonds real Bonds, and the bonded template only graph-matches a
    # copy that actually carries the Fe-donor bonds). Each copy's own coordination is
    # detected geometrically. The representative site's bonds were already created.
    ligand_names = {r.name for r in ligand_residues}
    structure = ligand_residues[0].structure
    handled = set(site.residues)
    n_copies = 0
    for r in structure.residues:
        if r.name not in ligand_names or r in handled:
            continue
        try:
            copy_site = detect_metal_site(r)
        except Exception:
            continue
        handled.update(copy_site.residues)
        made = False
        for (mm, dd) in copy_site.coordination:
            if dd not in mm.neighbors:
                try:
                    mm.structure.new_bond(mm, dd)
                    made = True
                except Exception:
                    pass
        if made:
            n_copies += 1

    session.logger.info(
        'Parameterised metal site %r; wrote %s and loaded templates %s. Applies to all '
        '%d copy/copies of %s this session; reload the ffXML in future sessions with '
        '"Load residue MD definition(s)".'
        % (site, xml_path, ', '.join(template_names.values()), 1 + n_copies,
           '/'.join(sorted(ligand_names))))
    return xml_path, template_names

