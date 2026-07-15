# @Author: Tristan Croll
# @Date:   14-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 14-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Headless test for divide-and-conquer AM1-BCC charging of large ligands:

* ``rdkit_bridge.safe_cut_bonds`` -- cuts land only in transferable alkyl regions
  (buffered from every heteroatom / charge / pi-system);
* ``rdkit_bridge.charge_partition`` -- size-bounded partition; small mols untouched;
* ``rdkit_bridge.build_capped_fragment`` -- valid capped fragments carrying fragSrc;
* **transferability** -- fragmented charges closely match a single whole-molecule
  AM1-BCC run (the correctness crux), and sum to the right integer.

The ANTECHAMBER calls are CLI subprocesses (no GUI), so this runs headless. Run:

    run_chimerax.bat --nogui --exit --script src/tests/test_charge_fragmentation.py

Prints PASS/FAIL and exits non-zero on failure.
'''


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def run(session):
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        print('SKIP: RDKit not available')
        return

    import os
    import shutil
    import tempfile
    from chimerax.isolde.atomic import rdkit_bridge as rb
    from chimerax.isolde.atomic.rdkit_bridge import FRAGSRC_PROP
    from chimerax.isolde.openmm.amberff import covalent as cov

    def _embed(smiles):
        m = Chem.AddHs(Chem.MolFromSmiles(smiles))
        if AllChem.EmbedMolecule(m, randomSeed=1) != 0:
            _fail('could not embed %s' % smiles)
        AllChem.MMFFOptimizeMolecule(m)
        return m

    # --- Part A: safe_cut_bonds only cuts buffered alkyl bonds --------------
    # Amine-...-alkyl-...-ester: cuts must sit in the saturated middle, never within
    # 2 bonds of the N, the ester O's or the carbonyl C.
    m = _embed('NCCCCCCCCCCCCOC(=O)C')
    cuts = rb.safe_cut_bonds(m, buffer=2)
    if not cuts:
        _fail('safe_cut_bonds found no cuts in a long alkyl chain')
    dmat = Chem.GetDistanceMatrix(m)
    polar = [a.GetIdx() for a in m.GetAtoms() if rb._is_polarish(a)]
    for bidx in cuts:
        b = m.GetBondWithIdx(bidx)
        for end in (b.GetBeginAtom(), b.GetEndAtom()):
            if end.GetAtomicNum() != 6 or end.GetIsAromatic():
                _fail('a cut endpoint is not an aliphatic carbon')
            if min(dmat[end.GetIdx()][p] for p in polar) <= 2:
                _fail('a cut endpoint is within the 2-bond buffer of a polar atom')
    print('PASS: safe_cut_bonds cuts only buffered aliphatic C-C bonds')

    # --- Part B: charge_partition size bound + small-mol passthrough --------
    if len(rb.charge_partition(m, max_heavy=60)) != 1:
        _fail('charge_partition should not split a molecule under the threshold')
    parts = rb.charge_partition(m, max_heavy=6)
    if len(parts) < 2:
        _fail('charge_partition did not split a long chain at a small threshold')
    covered = set().union(*parts)
    if covered != set(range(m.GetNumAtoms())):
        _fail('charge_partition does not cover every atom exactly once')
    if sum(len(p) for p in parts) != m.GetNumAtoms():
        _fail('charge_partition fragments overlap')
    for p in parts:
        heavy = sum(1 for i in p if m.GetAtomWithIdx(i).GetAtomicNum() != 1)
        # A fragment may exceed the bound only if it is one rigid piece; here the
        # ester head (a few heavy atoms) is the only non-alkyl piece, so all fit.
        if heavy > 6 + 3:
            _fail('a partition fragment (%d heavy) far exceeds the size bound' % heavy)
    print('PASS: charge_partition bounds fragment size and partitions cleanly')

    # --- Part C: build_capped_fragment ------------------------------------
    frag = rb.build_capped_fragment(m, parts[0])
    if frag is None:
        _fail('build_capped_fragment returned None')
    real = [a for a in frag.GetAtoms() if a.HasProp(FRAGSRC_PROP)]
    caps = [a for a in frag.GetAtoms() if a.HasProp('CAP')]
    if len(real) != len(parts[0]):
        _fail('capped fragment real-atom count %d != partition %d'
              % (len(real), len(parts[0])))
    if not caps:
        _fail('an interior fragment should carry at least one methyl cap')
    if frag.GetNumConformers() == 0:
        _fail('capped fragment has no conformer')
    print('PASS: build_capped_fragment yields a capped, fragSrc-tagged fragment')

    # --- Part D: transferability -- fragmented vs whole-molecule AM1-BCC ---
    # Palmitic acid: small enough that a whole-molecule sqm is quick, long enough to
    # force safe cuts at a small threshold. Neutral (net charge 0).
    acid = _embed('CCCCCCCCCCCCCCCC(=O)O')
    amber_bin = cov._amber_bin()
    whole_dir = tempfile.mkdtemp(prefix='isolde_fragtest_whole_')
    frag_base = tempfile.mkdtemp(prefix='isolde_fragtest_frag_')
    try:
        _t, whole_q, _o = cov._antechamber_on_mol(amber_bin, acid, 0, whole_dir,
                                                  'bcc', 'palmitic acid (whole)')
        partition = rb.charge_partition(acid, max_heavy=8)
        if len(partition) < 2:
            _fail('palmitic acid did not partition at max_heavy=8')
        frag_q = cov._bcc_charges_fragmented(session, amber_bin, acid, partition,
                                             frag_base)
    except Exception as e:
        if 'AMBERHOME' in str(e) or 'No such file' in str(e):
            _fail('ANTECHAMBER toolchain unavailable: %s' % e)
        raise
    finally:
        shutil.rmtree(whole_dir, ignore_errors=True)
        shutil.rmtree(frag_base, ignore_errors=True)

    if abs(sum(frag_q)) > 1e-3:
        _fail('fragmented charges do not sum to the net charge (0): %.4f' % sum(frag_q))
    dev = [abs(w - f) for w, f in zip(whole_q, frag_q)]
    max_dev = max(dev)
    rms = (sum(d * d for d in dev) / len(dev)) ** 0.5
    # Alkyl charges are highly transferable; buffered cuts should keep the whole-vs-
    # fragmented deviation small. Generous bound (this is a characterisation test).
    if max_dev > 0.06:
        _fail('fragmented vs whole-molecule charge deviation too large: max %.4f' % max_dev)
    print('PASS: fragmented AM1-BCC matches whole-molecule (max dev %.4f e, rms %.4f e)'
          % (max_dev, rms))

    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
