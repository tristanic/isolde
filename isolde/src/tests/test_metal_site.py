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

    # --- Part I: polynuclear cluster (Fe2S2) -- core, Urey-Bradley, IDEALISED geometry ---
    # A synthetic [2Fe-2S] named 'FES' (in CLUSTER_IDEAL_COORDS) but with DELIBERATELY
    # DISTORTED coordinates. Exercises the CLUSTER path: bridging sulfides excluded from
    # the organic build, per-atom-unique types, stiff Urey-Bradley 1-3 terms, and -- the
    # key check -- that internal geometry is taken from the curated CCD-ideal reference,
    # NOT the (distorted) modelled coordinates. Emitted template must build in OpenMM.
    sc = AtomicStructure(session, name='cluster-test', auto_style=False)
    try:
        clu = sc.new_residue('FES', 'A', 1)          # in CLUSTER_IDEAL_COORDS (no network)
        # Distorted rhombus: Fe-S ~3.1 A (blown up), far from the ideal ~2.2 A.
        cc = {'FE1': (0.0, 0.0, 0.0), 'FE2': (4.30, 0.0, 0.0),
              'S1': (2.15, 2.25, 0.0), 'S2': (2.15, -2.25, 0.0)}
        catoms = {}
        for nm, xyz in cc.items():
            a = sc.new_atom(nm, 'Fe' if nm.startswith('FE') else 'S')
            a.coord = numpy.array(xyz)
            clu.add_atom(a)
            catoms[nm] = a
        for fe in ('FE1', 'FE2'):
            for s in ('S1', 'S2'):
                sc.new_bond(catoms[fe], catoms[s])   # real coordination bonds
        dist_input = numpy.linalg.norm(catoms['FE1'].coord - catoms['S1'].coord)  # ~3.1 A

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

        ctn = {clu: 'MMET_FES'}
        cmt = cov._build_metal_terms(session, csite, ff, ctn, 0.5, core_atoms=core)
        if len(cmt['core']) != 2 or len(cmt['metals']) != 2:
            _fail('expected 2 metals + 2 core, got %d/%d'
                  % (len(cmt['metals']), len(cmt['core'])))
        etypes = [x['etype'] for x in cmt['metals'] + cmt['core']]
        if len(set(etypes)) != 4:
            _fail('cluster atoms must be UNIQUELY typed, got %s' % etypes)
        if not cmt['urey_bradley']:
            _fail('cluster produced no Urey-Bradley terms')
        # KEY: internal Fe-S bonds must take the IDEALISED length (~0.22 nm), NOT the
        # distorted modelled ~0.31 nm.
        internal = [cb for cb in cmt['coord_bonds'] if cb['donor'][1] in ('S1', 'S2')]
        if not internal or not all(0.21 < cb['length'] < 0.23 for cb in internal):
            _fail('internal Fe-S bonds not idealised (input was %.2f A): %s'
                  % (dist_input, [round(cb['length'] * 10, 3) for cb in internal]))
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
            tr = top.addResidue('MMET_FES', ch)
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
        print('PASS: Fe2S2 cluster -- core excluded, unique types, idealised geometry, '
              'Urey-Bradley, builds in OpenMM')
    finally:
        session.models.close([sc])

    # --- Part J: uncurated cluster geometry source is OPT-IN (no silent network) ------
    # An uncurated cluster with no exemplar and fetch OFF must yield no reference map
    # (-> modelled-geometry fallback, no network); a user exemplar supplies it.
    su = AtomicStructure(session, name='uncurated', auto_style=False)
    ex = AtomicStructure(session, name='exemplar', auto_style=False)
    try:
        def _mk(struct):
            r = struct.new_residue('UNC', 'A', 1)   # not in CLUSTER_IDEAL_COORDS
            cu = {'FE1': (0., 0., 0.), 'FE2': (2.7, 0., 0.),
                  'S1': (1.35, 1.55, 0.), 'S2': (1.35, -1.55, 0.)}
            at = {}
            for nm, xyz in cu.items():
                a = struct.new_atom(nm, 'Fe' if nm.startswith('FE') else 'S')
                a.coord = numpy.array(xyz); r.add_atom(a); at[nm] = a
            for fe in ('FE1', 'FE2'):
                for s in ('S1', 'S2'):
                    struct.new_bond(at[fe], at[s])
            return r
        ur = _mk(su)
        _mk(ex)
        usite = cov.detect_metal_site(ur)
        ucore = cov._cluster_core_atoms(usite)
        m_off = cov._cluster_reference_coord_map(session, usite.metals, ucore,
                                                 reference_model=None, fetch_reference=False)
        if m_off:
            _fail('uncurated cluster returned a reference map without opt-in (should be '
                  'empty -> modelled-geometry fallback): %d' % len(m_off))
        m_ex = cov._cluster_reference_coord_map(session, usite.metals, ucore,
                                                reference_model=ex, fetch_reference=False)
        if len(m_ex) != 4:
            _fail('user exemplar did not supply cluster geometry: %d' % len(m_ex))
        print('PASS: uncurated cluster geometry is opt-in (exemplar / fetch), never silent')
    finally:
        session.models.close([su, ex])

    # --- Part K: re-typed cysteine keeps its SG-touching bonded terms -----------------
    # Re-typing SG to a cluster type orphans the ff14SB CB-SG bond + X-CB-SG angles (they
    # were keyed by SG's base type 'SH'). Confirm (a) the ff14SB tables carry those terms,
    # and (b) covalent_to_ffxml emits pre-resolved extra_bonds/extra_angles + the aux CYS
    # template, so nothing through SG is left unrestrained.
    amber_xml = os.path.join(amberdir, 'amberff14SB.xml')
    bonds_tab, angles_tab = cov._amber_bonded_tables(amber_xml)
    if frozenset(('CT', 'SH')) not in bonds_tab:
        _fail('ff14SB CT-SH (Cys CB-SG) bond not found by _amber_bonded_tables')
    if ('CX', 'CT', 'SH') not in angles_tab:
        _fail('ff14SB CX-CT-SH (CA-CB-SG) angle not found by _amber_bonded_tables')
    cym = cov._read_amber_residue(amber_xml, 'CYM')
    if cym is None or not any(a['name'] == 'SG' for a in cym['atoms']):
        _fail('could not read the CYM template for re-templating')

    # Minimal metal_terms exercising the emission of a re-typed cysteine.
    R = object()
    mt = {'metals': [{'residue': R, 'name': 'FE', 'etype': 'MMET_X_FE', 'charge': 1.0}],
          'core': [], 'atom_types': [{'name': 'MMET_X_FE', 'element': 'Fe', 'mass': '55.85'},
                                     {'name': 'MMET_X_SG', 'element': 'S', 'mass': '32.06'}],
          'lj': [{'type': 'MMET_X_FE', 'sigma': '0.2', 'epsilon': '0.01'},
                 {'type': 'MMET_X_SG', 'sigma': '0.356', 'epsilon': '1.046'}],
          'coord_bonds': [], 'coord_angles': [], 'donor_deltas': [], 'donor_types': {},
          'urey_bradley': [],
          'aux_residues': [{'name': 'MMET_X_CYS',
                            'atoms': [{'name': 'CB', 'type': 'CT', 'charge': -0.24},
                                      {'name': 'SG', 'type': 'MMET_X_SG', 'charge': -0.63}],
                            'bonds': [('CB', 'SG')], 'external': ['SG']}],
          'cys_overrides': [],
          'extra_bonds': [('CT', 'MMET_X_SG', 0.181, 100000.0)],
          'extra_angles': [('CX', 'CT', 'MMET_X_SG', 1.9, 400.0),
                           ('CT', 'MMET_X_SG', 'MMET_X_FE', 1.9, 400.0)]}
    kdir = tempfile.mkdtemp(prefix='isolde_cystest_')
    try:
        kxml = os.path.join(kdir, 'cys.xml')
        ac.covalent_to_ffxml({}, {R: 'MMET_X'}, None, None, kxml,
                             os.path.join(amberdir, 'gaff2.xml'), metal_terms=mt)
        import xml.etree.ElementTree as ET
        kroot = ET.parse(kxml).getroot()
        if 'MMET_X_CYS' not in [r.get('name') for r in kroot.findall('Residues/Residue')]:
            _fail('re-typed CYS aux template not emitted')
        cbsg = [b for b in kroot.findall('HarmonicBondForce/Bond')
                if set((b.get('type1'), b.get('type2'))) == {'CT', 'MMET_X_SG'}]
        if not cbsg:
            _fail('CB-SG bond (CT-MMET_X_SG) missing from emitted template')
        sgang = [a for a in kroot.findall('HarmonicAngleForce/Angle')
                 if 'MMET_X_SG' in (a.get('type1'), a.get('type2'), a.get('type3'))]
        if len(sgang) < 2:
            _fail('SG-touching angles missing from emitted template: %d' % len(sgang))
    finally:
        shutil.rmtree(kdir, ignore_errors=True)
    print('PASS: re-typed cysteine SG keeps its bond + angles (nothing left unrestrained)')

    # --- Part L: bare-ion multi-site (zinc fingers) ----------------------------------
    # Donor-aware charge, auto-deprotonation, per-instance binding, and signature-based
    # grouping/collision-avoidance -- the fixes for the (Cys)4 Zn-finger failures. The
    # full parameterise_metal_site needs a live ISOLDE (GUI) session, so here we drive its
    # component helpers directly, using a lightweight stand-in `site` (real ChimeraX atoms,
    # explicit nonstandard set) since a synthetic residue is never perceived as PT_AMINO.

    # metal_charge_transfer (pure): donor-aware, conserved, clamped.
    qm, dl = mp.metal_charge_transfer(2, ['thiolate'] * 4)
    if abs(qm + sum(dl) - 2) > 1e-9:
        _fail('metal_charge_transfer does not conserve charge: %s' % ((qm, dl),))
    if not (0.35 < qm < 0.45):
        _fail('Zn(Cys)4 metal charge should be donor-aware ~+0.4, got %.3f (was +1 with '
              'the blunt 50/50 split)' % qm)
    qm_his, _ = mp.metal_charge_transfer(2, ['imidazole'] * 4)
    if not (qm_his > qm + 0.5):
        _fail('a neutral-N (His) site should keep far more charge on the metal (%.2f) '
              'than a thiolate site (%.2f)' % (qm_his, qm))
    qm_c, dl_c = mp.metal_charge_transfer(2, ['thiolate'] * 8)   # over-strong -> clamp
    if qm_c < -1e-9 or abs(qm_c + sum(dl_c) - 2) > 1e-9:
        _fail('metal_charge_transfer clamp failed (metal driven negative): %.3f' % qm_c)
    print('PASS: metal_charge_transfer -- donor-aware (Zn~+0.4), conserved, clamped')

    class _Site:
        '''Stand-in for a MetalSite with an EXPLICIT nonstandard_residues set (so the
        coordinating cysteines count as standard without needing PT_AMINO perception).'''
        def __init__(self, metals, coordination, nonstandard, residues):
            self.metals = list(metals)
            self.coordination = list(coordination)
            self.nonstandard_residues = list(nonstandard)
            self.residues = list(residues)

    tet = numpy.array([[1., 1., 1.], [1., -1., -1.], [-1., 1., -1.], [-1., -1., 1.]])
    tet /= numpy.linalg.norm(tet[0])

    sL = AtomicStructure(session, name='znfinger-test', auto_style=False)
    try:
        def _zn_cys_site(centre, resbase, n_cys=4, protonate=True):
            '''Build a bare Zn (own residue) + n_cys CYS residues (SG donor, CB heavy
            neighbour, optional HG). Returns a _Site with nonstandard = [ZN only].'''
            znr = sL.new_residue('ZN', 'A', resbase)
            zn = sL.new_atom('ZN', 'Zn'); zn.coord = numpy.array(centre, float)
            znr.add_atom(zn)
            coordination, cys_res = [], []
            for i in range(n_cys):
                d = tet[i]
                cr = sL.new_residue('CYS', 'A', resbase + 1 + i)
                sg = sL.new_atom('SG', 'S'); sg.coord = zn.coord + 2.3 * d
                cb = sL.new_atom('CB', 'C'); cb.coord = zn.coord + 4.1 * d
                cr.add_atom(sg); cr.add_atom(cb); sL.new_bond(sg, cb)
                if protonate:
                    hg = sL.new_atom('HG', 'H')
                    hg.coord = sg.coord + numpy.array([0.1, 0.9, 0.2])
                    cr.add_atom(hg); sL.new_bond(sg, hg)
                coordination.append((zn, sg)); cys_res.append(cr)
            return _Site([zn], coordination, [znr], [znr] + cys_res), znr, cys_res

        siteL, znrL, cysL = _zn_cys_site([0., 0., 0.], 1, protonate=True)

        # Auto-deprotonation: SG-HG present -> not yet a thiolate; strip HG -> thiolate.
        if cov._retempl_donors(siteL):
            _fail('protonated Cys should NOT yet be re-templatable thiolates')
        n_removed = cov._deprotonate_coordinating_donors(session, siteL)
        if n_removed != 4:
            _fail('auto-deprotonation should have removed 4 HG, removed %d' % n_removed)
        retp = cov._retempl_donors(siteL)
        if len(retp) != 4:
            _fail('after deprotonation all 4 SG should be re-templatable thiolates: %d'
                  % len(retp))
        if any(any(nb.element.number == 1 for nb in d.neighbors) for d in retp):
            _fail('a deprotonated SG still carries a hydrogen')
        print('PASS: auto-deprotonation strips coordinating Cys HG -> thiolate donors')

        # Bare-ion gate + signature: identical sites share a signature; a different
        # composition (Cys2His2) gets a DISTINCT one (no template collision).
        if not cov._site_is_bare_ion(siteL):
            _fail('_site_is_bare_ion should be True for a single-atom Zn residue')
        sig1 = cov._metal_site_signature(siteL, [])

        siteL2, _z2, _c2 = _zn_cys_site([20., 0., 0.], 20, protonate=False)
        if cov._metal_site_signature(siteL2, []) != sig1:
            _fail('two identical Zn(Cys)4 sites must share a signature (grouping)')

        # A Zn(Cys)2(His)2 site: 2 thiolates + 2 imidazole N donors.
        znr3 = sL.new_residue('ZN', 'A', 40)
        zn3 = sL.new_atom('ZN', 'Zn'); zn3.coord = numpy.array([40., 0., 0.])
        znr3.add_atom(zn3)
        coord3 = []
        for i in range(2):
            cr = sL.new_residue('CYS', 'A', 41 + i)
            sg = sL.new_atom('SG', 'S'); sg.coord = zn3.coord + 2.3 * tet[i]
            cb = sL.new_atom('CB', 'C'); cb.coord = zn3.coord + 4.1 * tet[i]
            cr.add_atom(sg); cr.add_atom(cb); sL.new_bond(sg, cb)
            coord3.append((zn3, sg))
        for i in range(2, 4):
            hr = sL.new_residue('HIS', 'A', 41 + i)
            ne2 = sL.new_atom('NE2', 'N'); ne2.coord = zn3.coord + 2.1 * tet[i]
            hr.add_atom(ne2)
            coord3.append((zn3, ne2))
        siteL3 = _Site([zn3], coord3, [znr3], [znr3])
        if cov._metal_site_signature(siteL3, []) == sig1:
            _fail('Zn(Cys)4 and Zn(Cys)2His2 must have DISTINCT signatures (no collision)')
        print('PASS: site signature -- identical sites share, different compositions differ')

        # _build_metal_terms on the bare-ion Zn(Cys)4: donor-aware metal charge, sig-scoped
        # SG type, all four cysteines bound per-instance, aux CYS template, and an emitted
        # ffXML whose Zn-SG bond loads in OpenMM.
        tnamesL = {znrL: 'MMET_' + sig1}
        lmt = cov._build_metal_terms(session, siteL, ff, tnamesL, 0.5, sig=sig1)
        if not (0.35 < lmt['metals'][0]['charge'] < 0.45):
            _fail('bare-ion Zn charge should be donor-aware ~+0.4, got %.3f'
                  % lmt['metals'][0]['charge'])
        sgt = 'MMET_%s_SG' % sig1
        if not lmt['donor_types'] or any(v != sgt for v in lmt['donor_types'].values()):
            _fail('SG donor type not signature-scoped: %s' % lmt['donor_types'])
        cyst = 'MMET_%s_CYS' % sig1
        if sorted(t for (_r, t) in lmt['cys_overrides']) != [cyst] * 4:
            _fail('cys_overrides did not bind all 4 cysteines to %s: %s'
                  % (cyst, lmt['cys_overrides']))
        if not lmt['aux_residues']:
            _fail('no aux CYS template built for the re-typed thiolates')

        gaff2_xml = os.path.join(amberdir, 'gaff2.xml')
        ldir = tempfile.mkdtemp(prefix='isolde_znfinger_')
        try:
            lxml = os.path.join(ldir, 'zn.xml')
            ac.covalent_to_ffxml({}, tnamesL, None, None, lxml, gaff2_xml, metal_terms=lmt)
            import xml.etree.ElementTree as ET
            lroot = ET.parse(lxml).getroot()
            metal_et = lmt['metals'][0]['etype']
            znsg = [b for b in lroot.findall('HarmonicBondForce/Bond')
                    if set((b.get('type1'), b.get('type2'))) == {metal_et, sgt}]
            if not znsg:
                _fail('Zn-SG coordination bond (%s-%s) missing from emitted template'
                      % (metal_et, sgt))
            if cyst not in [r.get('name') for r in lroot.findall('Residues/Residue')]:
                _fail('aux CYS template %s not emitted' % cyst)
            ffL = OMMFF(os.path.join(amberdir, 'amberff14SB.xml'),
                        os.path.join(amberdir, 'gaff2.xml'))
            ffL.loadFile(lxml)
            if ('MMET_' + sig1) not in ffL._templates or cyst not in ffL._templates:
                _fail('metal / cys templates not registered after loadFile')
            print('PASS: bare-ion Zn(Cys)4 -- donor-aware charge, per-instance Cys binding, '
                  'Zn-SG bond loads in OpenMM')

            # Per-instance override beats the metal_name_map free-ion fallback (the user's
            # clobber-fix): a ZN carrying isolde_template_name resolves to OUR bonded
            # template, not the built-in +1 free ion -- so a second site can't silently
            # revert to an unrestrained free ion.
            from chimerax.isolde.openmm.openmm_interface import find_residue_templates
            from chimerax.atomic import Residues
            ffP = OMMFF(os.path.join(amberdir, 'amberff14SB.xml'),
                        os.path.join(amberdir, 'gaff2.xml'),
                        os.path.join(amberdir, 'tip3p_HFE_multivalent.xml'))
            ffP.loadFile(lxml)              # registers MMET_<sig> alongside the free-ion ZN
            znrL.isolde_template_name = 'MMET_' + sig1
            tmap = find_residue_templates(Residues([znrL]), ffP, logger=session.logger)
            if tmap.get(0) != 'MMET_' + sig1:
                _fail('per-instance metal override not honored (clobbered by metal_name_map '
                      'free-ion?): %s' % dict(tmap))
        finally:
            shutil.rmtree(ldir, ignore_errors=True)
        print('PASS: per-instance metal override survives the free-ion name-map fallback')
    finally:
        session.models.close([sL])

    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
