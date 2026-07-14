# @Author: Tristan Croll
# @Date:   13-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 13-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Headless characterisation test for the RDKit-side foundation of covalent-ligand
parameterisation (``chimerax.isolde.atomic.rdkit_bridge``):

* ``mol_to_sdf`` -- the SDF writer that authors ANTECHAMBER's input (correct
  bond orders straight from the mol, preserved atom order);
* ``super_residue_to_rdkit`` -- the multi-residue "super-residue" builder that
  caps the wider chain (ACE/NME on peptide cuts) ready for AM1-BCC.

These are the pure-RDKit / structure steps that need no OpenMM simulation, so they
run headlessly. The ANTECHAMBER call, charge repartitioning, ffXML emission and
template matching need a GUI simulation and are verified there. Run inside
ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_covalent_prep.py

Prints PASS/FAIL and exits non-zero on failure.
'''
import os

FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1pmx_1.pdb')


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def run(session):
    from chimerax.core.commands import run as run_cmd

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        # RDKit is a build dependency on the rdkit branch; if it is somehow
        # absent, skip (no ALL PASS sentinel -> run_tests.bat records a skip).
        print('SKIP: RDKit not available')
        return

    from chimerax.isolde.atomic import rdkit_bridge as rb

    # --- Part A: SDF writer preserves order + full chemistry --------------
    mol = Chem.MolFromSmiles('CC(=O)Nc1ccccc1C(=O)[O-]')  # N-acetyl anthranilate
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=1) != 0:
        _fail('could not embed the test molecule')
    sdf = rb.mol_to_sdf(mol)
    back = Chem.MolFromMolBlock(sdf, removeHs=False)
    if back is None:
        _fail('RDKit could not re-parse the SDF we wrote')
    if back.GetNumAtoms() != mol.GetNumAtoms():
        _fail('SDF round-trip atom count %d != %d' % (back.GetNumAtoms(), mol.GetNumAtoms()))
    if [a.GetSymbol() for a in back.GetAtoms()] != [a.GetSymbol() for a in mol.GetAtoms()]:
        _fail('SDF round-trip did not preserve atom order')
    # Canonical SMILES equality proves bond orders, aromaticity and charge survived.
    if Chem.MolToSmiles(Chem.RemoveHs(back)) != Chem.MolToSmiles(Chem.RemoveHs(mol)):
        _fail('SDF round-trip changed the chemistry')
    print('PASS: mol_to_sdf preserves atom order and full chemistry')

    # --- Part B: super_residue_to_rdkit caps the wider chain --------------
    m = run_cmd(session, 'open %s' % FIXTURE)[0]
    try:
        existing = [r for r in m.chains[0].residues if r is not None]
        target = None
        for i in range(2, len(existing) - 2):
            r = existing[i]
            if {'N', 'CA', 'C', 'O'} <= set(r.atoms.names) and len(r.neighbors) >= 2:
                target = r
                break
        if target is None:
            _fail('no suitable internal residue in fixture')

        mol, cxmap, info = rb.super_residue_to_rdkit([target])
        if mol is None:
            _fail('super_residue_to_rdkit returned None (status=%s)' % info.get('status'))
        kinds = sorted(c['kind'] for c in info['caps'])
        if 'ACE' not in kinds or 'NME' not in kinds:
            _fail('expected ACE + NME caps for an internal residue, got %s' % kinds)
        # cxidx_to_atom covers all the residue's atoms (heavy + modelled H); the
        # template keeps the model's hydrogens so it graph-matches the residue.
        if len(cxmap) != target.num_atoms:
            _fail('cxmap size %d != residue atom count %d' % (len(cxmap), target.num_atoms))
        if any(atom.residue is not target for atom in cxmap.values()):
            _fail('a cxmap atom is outside the target residue')
        n_cxidx = sum(1 for a in mol.GetAtoms() if a.HasProp('cxIdx'))
        if n_cxidx != target.num_atoms:
            _fail('cxIdx-tagged atoms %d != residue atoms %d' % (n_cxidx, target.num_atoms))
        if not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
            _fail('no hydrogens added to the capped fragment')
        # The capped fragment must serialise to SDF for ANTECHAMBER.
        sdf2 = rb.mol_to_sdf(mol)
        rt = Chem.MolFromMolBlock(sdf2, removeHs=False)
        if rt is None or rt.GetNumAtoms() != mol.GetNumAtoms():
            _fail('capped-fragment SDF did not round-trip')
        print('PASS: super_residue_to_rdkit caps a single residue (ACE+NME) and serialises')

        # Two-residue unit: interior peptide bond kept, only outer bonds capped.
        pair = None
        for i in range(2, len(existing) - 3):
            a, b = existing[i], existing[i + 1]
            if b in a.neighbors and {'N', 'CA', 'C'} <= set(a.atoms.names) \
                    and {'N', 'CA', 'C'} <= set(b.atoms.names):
                pair = [a, b]
                break
        if pair is None:
            _fail('no bonded residue pair in fixture')
        _mol2u, cxmap2, info2 = rb.super_residue_to_rdkit(pair)
        if _mol2u is None:
            _fail('two-residue super_residue returned None (status=%s)' % info2.get('status'))
        pair_atoms = sum(r.num_atoms for r in pair)
        if len(cxmap2) != pair_atoms:
            _fail('two-residue cxmap %d != residue atoms %d' % (len(cxmap2), pair_atoms))
        kinds2 = sorted(c['kind'] for c in info2['caps'])
        if kinds2.count('ACE') != 1 or kinds2.count('NME') != 1:
            _fail('dipeptide unit should cap outer bonds only; got %s' % kinds2)
        print('PASS: super_residue_to_rdkit keeps the interior peptide bond of a 2-residue unit')

        # --- Part C: detect_covalent_unit -----------------------------
        import numpy
        from chimerax.core.errors import UserError
        from chimerax.isolde.openmm.amberff import covalent as cov

        # (i) A plain protein residue (only backbone bonds) has no unit. Check
        # this while the model is still pristine, before adding synthetic bonds.
        plain = next((r for r in existing if len(r.neighbors) >= 1), None)
        try:
            cov.detect_covalent_unit(plain)
            _fail('detect_covalent_unit accepted a plain protein residue')
        except UserError:
            pass
        print('PASS: detect_covalent_unit rejects a plain (backbone-only) residue')

        # Distinct hosts with a CB, well separated so the synthetic bonds below
        # do not interfere with one another.
        cb_hosts = [r for r in existing if 'CB' in set(r.atoms.names)]
        if len(cb_hosts) < 4:
            _fail('fixture has too few residues with a CB for the test')
        h_lig, h_x1, h_x2, h_mc = cb_hosts[0], cb_hosts[len(cb_hosts) // 3], \
            cb_hosts[2 * len(cb_hosts) // 3], cb_hosts[-1]

        def _bond_new_ligand(resnum, anchor_atom, name='LIG'):
            lg = m.new_residue(name, 'Z', resnum)
            c = m.new_atom('C1', 'C')
            c.coord = numpy.array(anchor_atom.coord) + numpy.array([1.5, 0.0, 0.0])
            lg.add_atom(c)
            m.new_bond(anchor_atom, c)
            return lg, c

        # (ii) Ligand bonded to a sidechain atom -> ligand + partner in the unit.
        lig, c1 = _bond_new_ligand(901, h_lig.find_atom('CB'))
        unit = cov.detect_covalent_unit(lig)
        if set(unit.residues) != {lig, h_lig}:
            _fail('sidechain-linked ligand unit wrong: %r' % unit)
        if h_lig not in unit.standard_residues or lig not in unit.nonstandard_residues:
            _fail('standard/nonstandard partition wrong for ligand case')
        if cov.detect_covalent_unit(h_lig).residues and \
                set(cov.detect_covalent_unit(h_lig).residues) != {lig, h_lig}:
            _fail('selecting the partner did not resolve the same unit')
        print('PASS: detect_covalent_unit resolves a ligand-attachment unit from either side')

        # (iii) Sidechain-sidechain crosslink between TWO STANDARD residues
        # (the His-Tyr / CcO case): no ligand at all.
        m.new_bond(h_x1.find_atom('CB'), h_x2.find_atom('CB'))
        xunit = cov.detect_covalent_unit(h_x1)
        if set(xunit.residues) != {h_x1, h_x2}:
            _fail('crosslink unit wrong: %r' % xunit)
        if xunit.nonstandard_residues:
            _fail('crosslink of two standard residues should have no nonstandard members')
        if len(xunit.links) != 1:
            _fail('crosslink should have exactly one link bond, got %d' % len(xunit.links))
        print('PASS: detect_covalent_unit handles a standard-standard sidechain crosslink')

        # (iv) A main-chain modification is now SUPPORTED (must not be refused).
        ligm, _cm = _bond_new_ligand(904, h_mc.find_atom('N'))
        munit = cov.detect_covalent_unit(ligm)
        if set(munit.residues) != {ligm, h_mc}:
            _fail('main-chain modification unit wrong: %r' % munit)
        print('PASS: detect_covalent_unit supports a main-chain modification')

        # (v) A free (unbonded) ligand has no covalent linkage -> UserError.
        free = m.new_residue('FRE', 'Z', 905)
        fa = m.new_atom('C1', 'C')
        fa.coord = numpy.array([0.0, 0.0, 0.0])
        free.add_atom(fa)
        try:
            cov.detect_covalent_unit(free)
            _fail('detect_covalent_unit accepted a free (unbonded) ligand')
        except UserError:
            pass
        print('PASS: detect_covalent_unit rejects a free ligand (points to isolde parameterise)')

        # --- Part D: typing / charge / ffXML pipeline (ANTECHAMBER faked) ---
        import os as _os
        import tempfile
        from chimerax.isolde.openmm.amberff import boundary_params as bp
        from chimerax.isolde.openmm.amberff import amber_convert as ac
        from chimerax.isolde.atomic import rdkit_bridge as rbmod
        from openmm.app import ForceField as OMMFF

        if bp.gaff_equivalent('CX') != 'c3' or bp.gaff_equivalent('S') != 'ss':
            _fail('gaff_equivalent map wrong')
        if bp.gaff_equivalent('zz_unknown') != 'zz_unknown':
            _fail('gaff_equivalent should pass unknown types through')
        if bp.lookup_seam('S', 'c3') is not None:
            _fail('curated seam registry should start empty')

        pr = {'A': {'atoms': [
            {'name': 'N', 'tag': 'frozen', 'charge': -0.4157},
            {'name': 'SG', 'tag': 'shell', 'charge': -0.31},
            {'name': 'C1', 'tag': 'ligand', 'charge': 0.20}]}}
        cov.repartition_charges(pr)
        tot = sum(a['charge'] for a in pr['A']['atoms'])
        if abs(tot - round(tot)) > 1e-9:
            _fail('repartition_charges did not reach an integer total')
        if abs(pr['A']['atoms'][0]['charge'] - (-0.4157)) > 1e-12:
            _fail('repartition_charges modified a frozen atom')
        print('PASS: gaff_equivalent/lookup_seam + repartition_charges (integer, frozen kept)')

        # Drive assign -> emit -> OpenMM load with a faked ANTECHAMBER result.
        amberdir = _os.path.dirname(ac.__file__)
        ff = OMMFF(_os.path.join(amberdir, 'amberff14SB.xml'),
                   _os.path.join(amberdir, 'gaff2.xml'))
        lig_unit = cov.detect_covalent_unit(lig)
        umol, ucx, uinfo = rbmod.super_residue_to_rdkit(lig_unit.residues)
        if umol is None:
            _fail('super_residue_to_rdkit failed for the pipeline test')
        gbe = {'C': 'c3', 'N': 'n', 'O': 'o', 'S': 'ss', 'H': 'hc', 'P': 'p5'}
        nn = umol.GetNumAtoms()
        ante = {'types': [gbe.get(umol.GetAtomWithIdx(i).GetSymbol(), 'c3') for i in range(nn)],
                'charges': [0.0] * nn, 'frcmod': None}
        prr = cov.assign_types_and_charges(lig_unit, umol, ucx, ante, ff, shell_radius=1)
        cov.repartition_charges(prr)
        nspec = next((s for s in prr[h_lig]['atoms'] if s['name'] == 'N'), None)
        if nspec is None or nspec['type'] != 'N':
            _fail('assign_types_and_charges did not freeze the ff14SB backbone N')
        tnames = cov._make_template_names(lig_unit, ff)
        tmpdir = tempfile.mkdtemp(prefix='isolde_covtest_')
        try:
            xmlp = _os.path.join(tmpdir, 'covtest.xml')
            gaff2_xml = _os.path.join(amberdir, 'gaff2.xml')
            ac.covalent_to_ffxml(prr, tnames, umol, None, xmlp, gaff2_xml)
            ff2 = OMMFF(_os.path.join(amberdir, 'amberff14SB.xml'),
                        _os.path.join(amberdir, 'gaff2.xml'))
            ff2.loadFile(xmlp)
            for tn in tnames.values():
                if tn not in ff2._templates:
                    _fail('template %s not registered after loadFile' % tn)
        finally:
            import shutil
            shutil.rmtree(tmpdir, ignore_errors=True)
        print('PASS: covalent_to_ffxml output loads into OpenMM (%s)' % list(tnames.values()))

        # --- Part E: a real ANTECHAMBER run (exercises the AMBERHOME setup) ---
        # ANTECHAMBER is a CLI subprocess (no GUI needed), so this runs headless.
        # A FREE 2-carbon ligand with no modelled H: since the whole unit is
        # unprotonated, the H policy protonates it (ethane), giving ANTECHAMBER a
        # valid molecule. (A partially-protonated unit is a contract violation and
        # would be left bare -- ISOLDE requires correct hydrogens up front.)
        elig = m.new_residue('LIG', 'Z', 906)
        e1 = m.new_atom('C1', 'C')
        e1.coord = numpy.array([10.0, 10.0, 10.0])
        e2 = m.new_atom('C2', 'C')
        e2.coord = numpy.array([11.5, 10.0, 10.0])
        elig.add_atom(e1)
        elig.add_atom(e2)
        m.new_bond(e1, e2)
        eunit = cov.CovalentUnit([elig], [])
        emol, ecx, einfo = rbmod.super_residue_to_rdkit(eunit.residues)
        if emol is None:
            _fail('super_residue_to_rdkit failed for the ANTECHAMBER test')
        try:
            ante = cov.run_antechamber(session, emol, einfo['net_charge'])
        except Exception as e:
            # Only a genuinely missing toolchain is a skip; anything else is a fail.
            if 'AMBERHOME' in str(e) or 'No such file' in str(e):
                _fail('run_antechamber could not set up ANTECHAMBER: %s' % e)
            raise
        try:
            if any(t is None for t in ante['types']):
                _fail('run_antechamber left some atoms untyped')
            if len(ante['types']) != emol.GetNumAtoms():
                _fail('run_antechamber type count mismatch')
            if not _os.path.isfile(ante['frcmod']):
                _fail('run_antechamber produced no frcmod')
            if abs(sum(ante['charges']) - einfo['net_charge']) > 0.05:
                _fail('AM1-BCC charge sum %.3f far from net charge %d'
                      % (sum(ante['charges']), einfo['net_charge']))
        finally:
            if ante.get('cleanup'):
                shutil.rmtree(ante['workdir'], ignore_errors=True)
        print('PASS: run_antechamber (real AM1-BCC/GAFF2) types + charges the unit')

        # --- Part F: free (non-covalent) ligand goes through the SAME pipeline ---
        # A free ligand reduces to a one-residue unit with no links -> no caps, all
        # atoms GAFF-typed, a self-contained USER_<name> template with prefixed
        # (collision-proof) atom types.
        flig = m.new_residue('FLG', 'Z', 907)
        fa = m.new_atom('C1', 'C'); fa.coord = numpy.array([0.0, 0.0, 0.0])
        fb = m.new_atom('C2', 'C'); fb.coord = numpy.array([1.5, 0.0, 0.0])
        fo = m.new_atom('O1', 'O'); fo.coord = numpy.array([2.2, 1.1, 0.0])
        for a in (fa, fb, fo):
            flig.add_atom(a)
        m.new_bond(fa, fb); m.new_bond(fb, fo)
        funit = cov.CovalentUnit([flig], [])
        fmol, fcx, finfo = rbmod.super_residue_to_rdkit(funit.residues)
        if fmol is None or finfo['caps']:
            _fail('free ligand should build with no caps; got %s' % (finfo.get('caps'),))
        fante = {'types': [gbe.get(fmol.GetAtomWithIdx(i).GetSymbol(), 'c3')
                           for i in range(fmol.GetNumAtoms())],
                 'charges': [0.0] * fmol.GetNumAtoms(), 'frcmod': None}
        fpr = cov.assign_types_and_charges(funit, fmol, fcx, fante, ff)
        if set(s['tag'] for s in fpr[flig]['atoms']) != {'ligand'}:
            _fail('free ligand atoms should all be tagged "ligand"')
        cov.repartition_charges(fpr)
        ftmp = tempfile.mkdtemp(prefix='isolde_freetest_')
        try:
            fxml = _os.path.join(ftmp, 'free.xml')
            ac.covalent_to_ffxml(fpr, {flig: 'USER_FLG'}, fmol, None,
                                 fxml, _os.path.join(amberdir, 'gaff2.xml'))
            ff3 = OMMFF(_os.path.join(amberdir, 'amberff14SB.xml'),
                        _os.path.join(amberdir, 'gaff2.xml'))
            ff3.loadFile(fxml)
            if 'USER_FLG' not in ff3._templates:
                _fail('free-ligand USER_FLG template not registered')
            import xml.etree.ElementTree as _ET
            atypes = [a.get('type') for a in
                      _ET.parse(fxml).getroot().findall('Residues/Residue/Atom')]
            if not all(t.startswith('USER_FLG_') for t in atypes):
                _fail('free-ligand atom types not all prefixed: %s' % atypes)
        finally:
            shutil.rmtree(ftmp, ignore_errors=True)
        print('PASS: free ligand -> self-contained USER_ template with prefixed types')

        # --- Part G: unparameterisable term -> UserError (fail loudly) ---
        # A frcmod with no bond/angle parameters stands in for GAFF failing to cover
        # the chemistry; covalent_to_ffxml must abort rather than emit a template
        # with unrestrained (zeroed) bonds/angles.
        gtmp = tempfile.mkdtemp(prefix='isolde_zerotest_')
        try:
            empty_frcmod = _os.path.join(gtmp, 'empty.frcmod')
            with open(empty_frcmod, 'w') as f:
                f.write('remark goes here\nMASS\n\nBOND\n\nANGLE\n\n'
                        'DIHE\n\nIMPROPER\n\nNONBON\n\n')
            raised = False
            try:
                ac.covalent_to_ffxml(fpr, {flig: 'USER_FLG'}, fmol, empty_frcmod,
                                     _os.path.join(gtmp, 'zero.xml'),
                                     _os.path.join(amberdir, 'gaff2.xml'))
            except UserError:
                raised = True
            if not raised:
                _fail('covalent_to_ffxml did not raise on unparameterisable terms')
        finally:
            shutil.rmtree(gtmp, ignore_errors=True)
        print('PASS: covalent_to_ffxml raises UserError on zeroed/missing bond-angle params')
    finally:
        session.models.close([m])

    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
