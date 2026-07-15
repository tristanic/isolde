# @Author: Tristan Croll
# @Date:   14-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 14-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Headless characterisation test for the metal-site MD parameterisation pieces that
need no OpenMM simulation:

* ``metal_params`` -- the pure empirical-parameter helpers (curated table + UFF
  fallback distances, angle snapping, charge splitting, oxidation-state guess);
* ``covalent.detect_metal_site`` + ``_find_donors`` -- geometric coordination-site
  detection;
* ``rdkit_bridge.super_residue_to_rdkit(exclude=...)`` -- the metal-free organic
  build (the metal is left out, its coordination bond NOT capped);
* ``covalent._build_metal_terms`` + ``amber_convert.covalent_to_ffxml(metal_terms=)``
  -- the bonded metal-site template emission, which must LOAD in OpenMM.

The running simulation (coordination held, no LJ blow-off) needs the GUI and is
verified there. Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_metal_site.py

Prints PASS/FAIL and exits non-zero on failure.
'''
import os


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def run(session):
    try:
        from rdkit import Chem  # noqa: F401
    except ImportError:
        print('SKIP: RDKit not available')
        return

    import math
    import numpy
    import tempfile
    import shutil
    from chimerax.core.errors import UserError
    from chimerax.atomic import AtomicStructure
    from chimerax.isolde.atomic import rdkit_bridge as rb
    from chimerax.isolde.openmm.amberff import covalent as cov
    from chimerax.isolde.openmm.amberff import amber_convert as ac
    from chimerax.isolde.openmm.amberff import metal_params as mp
    from openmm.app import ForceField as OMMFF

    amberdir = os.path.dirname(ac.__file__)

    # --- Part A: metal_params pure helpers --------------------------------
    r0, k = mp.metal_bond_params('Fe', 'N')          # curated
    if not (0.19 < r0 < 0.21 and k > 0):
        _fail('metal_bond_params Fe-N off: %s' % ((r0, k),))
    uff = mp.metal_bond_params('W', 'O')             # not curated -> UFF radii
    if uff is None or not (0.20 < uff[0] < 0.28):
        _fail('metal_bond_params UFF fallback (W-O) off: %s' % (uff,))
    if mp.metal_bond_params('Fe', 'Xx') is not None:
        _fail('metal_bond_params should fail for an unknown donor element')
    # Force constant is DERIVED from the measured stretch frequency (k = mu(2 pi c nu)^2).
    # Fe-N(porphyrin) at 340 cm^-1 must reproduce the vetted ~46000 kJ/mol/nm^2 (BRYCE_HEM).
    k_fe_n = mp.bond_k_from_frequency('Fe', 'N', 340.0)
    if not (44000 < k_fe_n < 48000):
        _fail('bond_k_from_frequency(Fe-N,340) should be ~46000, got %.0f' % k_fe_n)
    # Physical fidelity: the axial His-Fe bond is WEAKER than the equatorial pyrrole bond.
    _, k_eq = mp.metal_bond_params('Fe', 'N')
    _, k_ax = mp.metal_bond_params('Fe', 'N', axial=True)
    if not (k_ax < k_eq):
        _fail('axial Fe-N (%.0f) should be softer than equatorial (%.0f)' % (k_ax, k_eq))
    # Coordination angles are firm (they, not the soft dative bond, hold a lone axial donor).
    if mp.angle_k() < 400:
        _fail('coordination angle k should be firm (>=400), got %.0f' % mp.angle_k())
    # angle snapping to the nearest canonical value
    if abs(mp.snap_angle(math.radians(88.0)) - math.radians(90.0)) > 1e-6:
        _fail('snap_angle(88) should snap to 90')
    if abs(mp.snap_angle(math.radians(172.0)) - math.radians(180.0)) > 1e-6:
        _fail('snap_angle(172) should snap to 180')
    if abs(mp.snap_angle(math.radians(112.0)) - math.radians(109.47)) > 1e-6:
        _fail('snap_angle(112) should snap to tetrahedral')
    # charge split conserves the oxidation-state charge
    qm, dq = mp.metal_charge_split(2, 4, keep_fraction=0.5)
    if abs(qm + 4 * dq - 2) > 1e-9 or abs(qm - 1.0) > 1e-9:
        _fail('metal_charge_split does not conserve charge: %s' % ((qm, dq),))
    if mp.guess_ox_state('Zn') != 2 or mp.guess_ox_state('Fe', 3) != 3:
        _fail('guess_ox_state wrong')
    print('PASS: metal_params (table + UFF fallback, angle snap, charge split)')

    # --- Build a synthetic Zn-bis(amine) site -----------------------------
    # Zn (its own residue) with two N donors ~2.1 A away, each on a small organic
    # residue, so detection + the metal-free build + emission can all be exercised.
    s = AtomicStructure(session, name='metal-test', auto_style=False)
    try:
        znr = s.new_residue('ZN', 'A', 1)
        zn = s.new_atom('ZN', 'Zn')
        zn.coord = numpy.array([0.0, 0.0, 0.0])
        znr.add_atom(zn)

        lig = s.new_residue('LIG', 'A', 2)
        coords = {'N1': (2.10, 0.0, 0.0), 'C1': (3.0, 1.0, 0.0),
                  'C2': (1.0, 3.0, 0.0), 'N2': (0.0, 2.10, 0.0)}
        atoms = {}
        for nm, xyz in coords.items():
            a = s.new_atom(nm, nm[0])
            a.coord = numpy.array(xyz)
            lig.add_atom(a)
            atoms[nm] = a
        s.new_bond(atoms['N1'], atoms['C1'])
        s.new_bond(atoms['C1'], atoms['C2'])
        s.new_bond(atoms['C2'], atoms['N2'])

        # --- Part B: detection ---------------------------------------------
        site = cov.detect_metal_site(lig)
        if set(site.metals) != {zn}:
            _fail('detect_metal_site did not find the Zn')
        if set(site.residues) != {znr, lig}:
            _fail('metal site residues wrong: %r' % site.residues)
        donors = {d for (m, d) in site.coordination}
        if donors != {atoms['N1'], atoms['N2']}:
            _fail('coordination donors wrong: %s' % sorted(a.name for a in donors))
        print('PASS: detect_metal_site finds the metal, donors and site residues')

        # Selecting the metal itself resolves the same site.
        if set(cov.detect_metal_site(znr).metals) != {zn}:
            _fail('selecting the metal did not resolve the site')
        # A metal-free organic residue is refused (points to isolde parameterise).
        og = s.new_residue('ORG', 'A', 9)
        oa = s.new_atom('C1', 'C'); oa.coord = numpy.array([20.0, 20.0, 20.0])
        og.add_atom(oa)
        try:
            cov.detect_metal_site(og)
            _fail('detect_metal_site accepted a metal-free residue')
        except UserError:
            pass
        print('PASS: detect_metal_site rejects a metal-free residue')

        # --- Part C: metal-free build leaves the metal out, uncapped -------
        s.new_bond(zn, atoms['N1'])          # real coordination bonds
        s.new_bond(zn, atoms['N2'])
        mol, cxmap, info = rb.super_residue_to_rdkit(site.residues, exclude={zn})
        if mol is None:
            _fail('metal-free super_residue returned None (status=%s)' % info.get('status'))
        if any(a.GetSymbol() == 'Zn' for a in mol.GetAtoms()):
            _fail('the metal was not excluded from the organic build')
        if info['caps']:
            _fail('a coordination bond was capped (should be left open): %s' % info['caps'])
        if set(cxmap.values()) != set(lig.atoms):
            _fail('cxmap should cover exactly the organic (LIG) atoms')
        print('PASS: super_residue_to_rdkit(exclude=metal) omits the metal and does not cap it')

        # --- Part D: ion LJ + metal-term construction ----------------------
        # A forcefield that includes ISOLDE's bundled ion set, so _ion_lj resolves.
        ff = OMMFF(os.path.join(amberdir, 'amberff14SB.xml'),
                   os.path.join(amberdir, 'tip3p_standard.xml'),
                   os.path.join(amberdir, 'tip3p_HFE_multivalent.xml'),
                   os.path.join(amberdir, 'gaff2.xml'))
        lj = cov._ion_lj(ff, cov._ion_template_for_element('Zn', 2))
        if lj is None:
            _fail('_ion_lj could not read Zn Lennard-Jones from the bundled ion set')

        tnames = cov._make_template_names(site, ff, prefix='MMET_')
        if znr not in tnames or lig not in tnames:
            _fail('_make_template_names did not cover all site residues')
        mt = cov._build_metal_terms(session, site, ff, tnames, 0.5)
        if len(mt['metals']) != 1:
            _fail('expected one metal in metal_terms')
        if len(mt['coord_bonds']) != 2 or len(mt['coord_angles']) != 1:
            _fail('expected 2 coord bonds + 1 angle, got %d/%d'
                  % (len(mt['coord_bonds']), len(mt['coord_angles'])))
        th0 = mt['coord_angles'][0]['angle']
        if th0 not in mp.CANONICAL_ANGLES:
            _fail('coordination angle was not snapped to a canonical value')
        print('PASS: _ion_lj + _build_metal_terms (bonds, snapped angle, charge deltas)')

        # --- Part E: emit a bonded metal template that LOADS in OpenMM -----
        gbe = {'C': 'c3', 'N': 'n', 'O': 'o', 'H': 'hc'}
        nn = mol.GetNumAtoms()
        ante = {'types': [gbe.get(mol.GetAtomWithIdx(i).GetSymbol(), 'c3')
                          for i in range(nn)],
                'charges': [0.0] * nn, 'frcmod': None}
        per_res = cov.assign_types_and_charges(site, mol, cxmap, ante, ff)
        tmpdir = tempfile.mkdtemp(prefix='isolde_metaltest_')
        try:
            xmlp = os.path.join(tmpdir, 'metal.xml')
            gaff2_xml = os.path.join(amberdir, 'gaff2.xml')
            ac.covalent_to_ffxml(per_res, tnames, mol, None, xmlp, gaff2_xml,
                                 metal_terms=mt)
            import xml.etree.ElementTree as ET
            root = ET.parse(xmlp).getroot()
            metal_etype = mt['metals'][0]['etype']
            types = [t.get('name') for t in root.findall('AtomTypes/Type')]
            if metal_etype not in types:
                _fail('metal atom type %s missing from AtomTypes' % metal_etype)
            bonds = root.findall('HarmonicBondForce/Bond')
            if not any(metal_etype in (b.get('type1'), b.get('type2')) for b in bonds):
                _fail('no metal-donor bond emitted')
            angles = root.findall('HarmonicAngleForce/Angle')
            if not any(a.get('type2') == metal_etype for a in angles):
                _fail('no donor-metal-donor angle emitted')
            # The emitted template must load into a fresh OpenMM ForceField.
            ff2 = OMMFF(os.path.join(amberdir, 'amberff14SB.xml'),
                        os.path.join(amberdir, 'gaff2.xml'))
            ff2.loadFile(xmlp)
            for tn in tnames.values():
                if tn not in ff2._templates:
                    _fail('metal-site template %s not registered after loadFile' % tn)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)
        print('PASS: covalent_to_ffxml(metal_terms=) emits a bonded template that loads in OpenMM')

        # --- Part F: template_name_to_ccd_name maps MMET_ back to the CCD (Issue B) ---
        # So the residue-rebuild code can fetch the CCD coordinate template and place
        # the metal from ideal coords (instead of dropping it).
        from chimerax.isolde.openmm.amberff.template_utils import template_name_to_ccd_name
        ccd_name, _desc = template_name_to_ccd_name('MMET_HEC')
        if ccd_name != 'HEC':
            _fail('template_name_to_ccd_name(MMET_HEC) should give HEC, got %r' % ccd_name)
        print('PASS: template_name_to_ccd_name maps MMET_<ccd> back to the CCD id')

        # --- Part G: MMET_ template name-matches ALL copies (Issue A) ---
        from chimerax.isolde.openmm.openmm_interface import find_residue_templates
        from chimerax.atomic import Residues
        # ff2 (from Part E) has an MMET_LIG template loaded. Make a second LIG copy
        # (same name, no coordination) and confirm the name-pass binds both.
        lig2 = s.new_residue('LIG', 'A', 3)
        for nm, xyz in coords.items():
            a = s.new_atom(nm, nm[0])
            a.coord = numpy.array(xyz) + numpy.array([10.0, 0.0, 0.0])
            lig2.add_atom(a)
        tmap = find_residue_templates(Residues([lig, lig2]), ff2, logger=session.logger)
        if tmap.get(0) != 'MMET_LIG' or tmap.get(1) != 'MMET_LIG':
            _fail('MMET_ name-pass did not bind both LIG copies: %s' % dict(tmap))
        print('PASS: MMET_ template name-matches every copy of the metalloligand')

        # --- Part H: macrocycle selector picks the coordinating ring system only ---
        # amber_convert.macrocycle_atom_indices selects the connected component of
        # RING atoms reachable from the metal donors -- the porphyrin/chlorin core --
        # so its planar impropers can be stiffened while every non-ring substituent
        # (methyl, propionate, vinyl) keeps stock GAFF flexibility.
        from rdkit.Chem import MolFromSmiles
        # Pyridine (aromatic N donor) bearing a methyl and a carboxylate: the ring is
        # the coordinating core; both substituents must be excluded.
        pm = MolFromSmiles('Cc1ccncc1C(=O)[O-]')
        if pm is None:
            _fail('could not build the macrocycle-selector test molecule')
        ndx = next(a.GetIdx() for a in pm.GetAtoms() if a.GetSymbol() == 'N')
        core = ac.macrocycle_atom_indices(pm, [ndx])
        ring_idx = set(a.GetIdx() for a in pm.GetAtoms() if a.GetIsAromatic())
        if core != ring_idx:
            _fail('macrocycle selector did not return exactly the aromatic ring: '
                  '%s vs %s' % (sorted(core), sorted(ring_idx)))
        if any(not pm.GetAtomWithIdx(i).IsInRing() for i in core):
            _fail('macrocycle selector included a non-ring atom')
        # A seed that is NOT a ring atom yields nothing (no false macrocycle).
        methyl = next(a.GetIdx() for a in pm.GetAtoms()
                      if a.GetSymbol() == 'C' and not a.IsInRing()
                      and all(nb.GetSymbol() in ('C',) for nb in a.GetNeighbors()))
        if ac.macrocycle_atom_indices(pm, [methyl]):
            _fail('macrocycle selector should return nothing for a non-ring seed')
        # The scale constant must be a real stiffening (>1), else it is a no-op.
        if not (mp.MACROCYCLE_IMPROPER_SCALE > 1.0):
            _fail('MACROCYCLE_IMPROPER_SCALE must be > 1 to stiffen the macrocycle')
        print('PASS: macrocycle_atom_indices selects the coordinating ring core only')
    finally:
        session.models.close([s])

    # --- Part I: polynuclear cluster (Fe2S2) -- core handling + Urey-Bradley ----
    # A synthetic [2Fe-2S] rhombus in ONE residue: 2 Fe + 2 mu-bridging S, each S bonded
    # to both Fe. Exercises the CLUSTER path: bridging sulfides excluded from the organic
    # build, per-atom-unique types, observed-geometry bonds/angles incl. Fe-S-Fe, stiff
    # Urey-Bradley 1-3 terms, and an emitted template that BUILDS a system in OpenMM.
    sc = AtomicStructure(session, name='cluster-test', auto_style=False)
    try:
        clu = sc.new_residue('CLU', 'A', 1)
        cc = {'FE1': (0.0, 0.0, 0.0), 'FE2': (2.70, 0.0, 0.0),
              'S1': (1.35, 1.55, 0.0), 'S2': (1.35, -1.55, 0.0)}
        catoms = {}
        for nm, xyz in cc.items():
            a = sc.new_atom(nm, nm[:2] if nm.startswith('FE') else 'S')
            a.coord = numpy.array(xyz)
            clu.add_atom(a)
            catoms[nm] = a
        for fe in ('FE1', 'FE2'):
            for s in ('S1', 'S2'):
                sc.new_bond(catoms[fe], catoms[s])   # real coordination bonds

        csite = cov.detect_metal_site(clu)
        if set(csite.metals) != {catoms['FE1'], catoms['FE2']}:
            _fail('cluster detection missed an Fe: %r' % csite.metals)
        core = cov._cluster_core_atoms(csite)
        if set(core) != {catoms['S1'], catoms['S2']}:
            _fail('_cluster_core_atoms wrong: %s' % sorted(a.name for a in core))

        # Metal-free build with metals AND core excluded -> nothing organic -> mol None.
        cmol, _cx, _info = rb.super_residue_to_rdkit(
            csite.residues, exclude=set(csite.metals) | set(core))
        if cmol is not None:
            _fail('a pure cluster should leave no organic framework (mol should be None)')

        ctn = {clu: 'MMET_CLU'}
        cmt = cov._build_metal_terms(session, csite, ff, ctn, 0.5, core_atoms=core)
        if len(cmt['core']) != 2 or len(cmt['metals']) != 2:
            _fail('expected 2 metals + 2 core, got %d/%d'
                  % (len(cmt['metals']), len(cmt['core'])))
        etypes = [x['etype'] for x in cmt['metals'] + cmt['core']]
        if len(set(etypes)) != 4:
            _fail('cluster atoms must be UNIQUELY typed, got %s' % etypes)
        if not cmt['urey_bradley']:
            _fail('cluster produced no Urey-Bradley terms')
        # Both Fe-S-Fe (core vertex) and S-Fe-S (metal vertex) angles must be present.
        vfe = {ca['metal'][1] for ca in cmt['coord_angles']}
        if 'S1' not in vfe and 'S2' not in vfe:
            _fail('no metal-core-metal (Fe-S-Fe) angle emitted (vertices: %s)' % vfe)
        # Net charge is an integer (formal metal-ox + core-anion), here 2*2 + 2*(-2) = 0.
        net = sum(x['charge'] for x in cmt['metals'] + cmt['core'])
        if abs(net - round(net)) > 1e-6:
            _fail('cluster net charge not integer: %.4f' % net)

        cdir = tempfile.mkdtemp(prefix='isolde_clustertest_')
        try:
            cxml = os.path.join(cdir, 'cluster.xml')
            ac.covalent_to_ffxml({}, ctn, None, None, cxml,
                                 os.path.join(amberdir, 'gaff2.xml'), metal_terms=cmt)
            import xml.etree.ElementTree as ET
            croot = ET.parse(cxml).getroot()
            tnames = [t.get('name') for t in croot.findall('AtomTypes/Type')]
            if len(tnames) != len(set(tnames)):
                _fail('duplicate <Type> emitted for the cluster: %s' % tnames)
            ubs = croot.findall('AmoebaUreyBradleyForce/UreyBradley')
            if not ubs:
                _fail('no <AmoebaUreyBradleyForce> emitted')
            xml_bonds = len(croot.findall('HarmonicBondForce/Bond'))
            # The emitted template must build a real OpenMM System, and the Urey-Bradley
            # 1-3 terms must be applied. OpenMM realises Urey-Bradley as extra
            # HarmonicBondForce bonds (force.addBond on the 1-3 pair), so confirm the
            # coordination bonds AND every UB term landed.
            from openmm.app import ForceField as _FF, Topology
            from openmm.app.element import Element as _El
            from openmm import HarmonicBondForce as _HBF
            ff3 = _FF(os.path.join(amberdir, 'amberff14SB.xml'),
                      os.path.join(amberdir, 'gaff2.xml'))
            ff3.loadFile(cxml)
            top = Topology()
            ch = top.addChain()
            tr = top.addResidue('MMET_CLU', ch)
            ta = {nm: top.addAtom(nm, _El.getBySymbol('Fe' if nm.startswith('FE') else 'S'), tr)
                  for nm in cc}
            for fe in ('FE1', 'FE2'):
                for s in ('S1', 'S2'):
                    top.addBond(ta[fe], ta[s])
            omm_sys = ff3.createSystem(top)
            n_hb = sum(f.getNumBonds() for f in omm_sys.getForces() if isinstance(f, _HBF))
            if n_hb != xml_bonds + len(ubs):
                _fail('expected %d harmonic bonds (%d coordination + %d Urey-Bradley), '
                      'got %d -- Urey-Bradley terms did not apply'
                      % (xml_bonds + len(ubs), xml_bonds, len(ubs), n_hb))
        finally:
            shutil.rmtree(cdir, ignore_errors=True)
        print('PASS: Fe2S2 cluster -- core excluded, unique types, Urey-Bradley, builds in OpenMM')
    finally:
        session.models.close([sc])

    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
