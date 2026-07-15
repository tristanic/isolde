# @Author: Tristan Croll
# @Date:   14-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 14-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Headless test for user-supplied chemistry exemplars (``baseTemplates``):

* ``rdkit_bridge._resolve_base_templates`` -- CCD id / SMILES resolution + skip of
  unresolvable strings;
* ``rdkit_bridge.apply_base_templates`` -- transfer of formal charges + bond orders
  from an exemplar onto a mis-perceived fragment, preserving cxName/cxIdx;
* end-to-end: ``super_residue_to_rdkit(..., base_templates=[...])`` overriding a
  wrong CCD-neutral charge with the exemplar's, giving the correct net charge.

Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_base_templates.py

Prints PASS/FAIL and exits non-zero on failure.
'''


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def run(session):
    try:
        from rdkit import Chem
    except ImportError:
        print('SKIP: RDKit not available')
        return

    import numpy
    from chimerax.atomic import AtomicStructure
    from chimerax.isolde.atomic import rdkit_bridge as rb

    # --- Part A: _resolve_base_templates ----------------------------------
    # 'CC(=O)[O-]' has SMILES-only chars (=,(),[]) so it can't be a CCD id -> SMILES
    # path deterministically. (A bare token like 'CCO' is a real CCD id and would
    # resolve CCD-first -- the documented, intended ambiguity.)
    mols = rb._resolve_base_templates(
        session, ['CC(=O)[O-]', 'not_a_real_ccd_or_smiles_zzz'])
    if len(mols) != 1:
        _fail('_resolve_base_templates: expected 1 resolved (SMILES only), got %d' % len(mols))
    if sum(1 for a in mols[0].GetAtoms() if a.GetAtomicNum() != 1) != 4:
        _fail('_resolve_base_templates: acetate SMILES should have 4 heavy atoms')
    print('PASS: _resolve_base_templates resolves SMILES and skips garbage')

    # --- Part B: apply_base_templates corrects charge/order, keeps identity ---
    # Mis-perceived acetate: neutral gem-diol-like skeleton (both C-O single, net 0).
    frag = Chem.MolFromSmiles('CC(O)O')
    for i, a in enumerate(frag.GetAtoms()):
        a.SetProp(rb.NAME_PROP, 'A%d' % i)
        a.SetIntProp(rb.INDEX_PROP, i)
    if Chem.GetFormalCharge(frag) != 0:
        _fail('setup: mis-perceived fragment should be neutral')
    ref = Chem.MolFromSmiles('CC(=O)[O-]')          # correct acetate, net -1
    out = rb.apply_base_templates(frag, [ref])
    if out is None:
        _fail('apply_base_templates returned None')
    if Chem.GetFormalCharge(out) != -1:
        _fail('apply_base_templates did not transfer the -1 charge (got %d)'
              % Chem.GetFormalCharge(out))
    if not any(b.GetBondType() == Chem.BondType.DOUBLE for b in out.GetBonds()):
        _fail('apply_base_templates did not transfer the C=O double bond')
    names = sorted(a.GetProp(rb.NAME_PROP) for a in out.GetAtoms()
                   if a.HasProp(rb.NAME_PROP))
    if names != ['A0', 'A1', 'A2', 'A3']:
        _fail('apply_base_templates did not preserve cxName identity (got %s)' % names)
    print('PASS: apply_base_templates transfers charge/order and preserves identity')

    # A too-small / non-matching exemplar leaves the fragment unchanged.
    frag2 = Chem.MolFromSmiles('CC(O)O')
    for i, a in enumerate(frag2.GetAtoms()):
        a.SetProp(rb.NAME_PROP, 'B%d' % i)
    unchanged = rb.apply_base_templates(frag2, [Chem.MolFromSmiles('N')], min_match=3)
    if Chem.GetFormalCharge(unchanged) != 0:
        _fail('a non-matching exemplar should not change the fragment')
    print('PASS: apply_base_templates ignores non-matching / tiny exemplars')

    # --- Part C: end-to-end override through super_residue_to_rdkit -------
    # Acetic acid (CCD 'ACY') modelled DEPROTONATED (omit the acid H). The CCD-by-name
    # path charges it neutral (net 0, wrong); a base template of the acetate SMILES
    # overrides that to the correct -1.
    rec = rb.ccd_records(session, 'ACY')
    if rec is None:
        print('SKIP-C: CCD component ACY unavailable offline')
    else:
        atoms, bonds, coords = rec
        o_ids = {aid for (aid, el, q, ar) in atoms if el.upper() == 'O'}
        acid_h = None
        elem = {aid: el.upper() for (aid, el, q, ar) in atoms}
        for (a1, a2, order, ar) in bonds:
            if elem.get(a1) == 'H' and a2 in o_ids:
                acid_h = a1
            elif elem.get(a2) == 'H' and a1 in o_ids:
                acid_h = a2

        def _build():
            s = AtomicStructure(session, name='acy', auto_style=False)
            r = s.new_residue('ACY', 'A', 1)
            am = {}
            for (aid, el, q, ar) in atoms:
                if aid == acid_h:
                    continue
                a = s.new_atom(aid, el)
                a.coord = numpy.array(coords.get(aid, (0.0, 0.0, 0.0)))
                r.add_atom(a)
                am[aid] = a
            for (a1, a2, order, ar) in bonds:
                if a1 in am and a2 in am:
                    s.new_bond(am[a1], am[a2])
            return s, r

        s0, r0 = _build()
        try:
            _m, _c, info_plain = rb.super_residue_to_rdkit([r0])
        finally:
            s0.session.models.close([s0])
        s1, r1 = _build()
        try:
            _m2, _c2, info_bt = rb.super_residue_to_rdkit(
                [r1], base_templates=['CC(=O)[O-]'])
        finally:
            s1.session.models.close([s1])

        if info_plain['net_charge'] != 0:
            _fail('deprotonated ACY via CCD path should read net 0, got %d'
                  % info_plain['net_charge'])
        if info_bt['net_charge'] != -1:
            _fail('base template should override ACY to net -1, got %d'
                  % info_bt['net_charge'])
        print('PASS: base template overrides CCD-neutral charge (net 0 -> -1) end-to-end')

    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
