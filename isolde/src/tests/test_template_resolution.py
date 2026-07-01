# @Author: Tristan Croll
# @Date:   30-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 30-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Characterisation test for ISOLDE's forcefield-template resolution -- the path
the future ML parameterisation backend will replace. Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_template_resolution.py

Driven through the read-only ``isolde preflight parameters`` dry-run (resolves
templates without building an OpenMM Context). Asserts the CURRENT behaviour for
the bundled 1pmx fixture against captured baselines, so a future change that
alters template assignment fails loudly. Prints PASS/FAIL lines and exits
non-zero on failure so it can gate a build.
'''
import os

FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1pmx_1.pdb')

# Baseline captured 2026-06 against amber14 + the bundled forcefield XMLs.
EXPECTED_COUNTS = {
    'n_total': 70, 'n_matched': 70, 'n_ambiguous': 0, 'n_unmatched': 0,
    'ready_for_simulation': True,
}
# find_residue_templates only assigns explicit templates to the cysteines (all
# disulfide-bonded -> CYX); standard residues are left to OpenMM's matcher.
EXPECTED_TEMPLATE_MAP = {
    '/A CYS 6': 'CYX', '/A CYS 18': 'CYX', '/A CYS 47': 'CYX',
    '/A CYS 48': 'CYX', '/A CYS 52': 'CYX', '/A CYS 61': 'CYX',
}


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def run(session):
    from chimerax.core.commands import run as run_cmd
    m = run_cmd(session, 'open "%s"' % FIXTURE.replace('\\', '/'), log=False)[0]
    try:
        from chimerax.isolde.validation.cmd import isolde_preflight_parameters

        # --- preflight dict shape + per-class counts -----------------------
        result = isolde_preflight_parameters(session, model=m, forcefield=None)
        for key in ('model', 'forcefield', 'n_total', 'n_matched', 'n_ambiguous',
                    'n_unmatched', 'ready_for_simulation', 'unmatched', 'ambiguous'):
            if key not in result:
                _fail('preflight result lost key %r' % key)
        if (result['n_matched'] + result['n_ambiguous'] + result['n_unmatched']
                != result['n_total']):
            _fail('preflight counts do not sum to n_total: %r' % result)
        for key, expect in EXPECTED_COUNTS.items():
            if result[key] != expect:
                _fail('preflight %s: expected %r, got %r' % (key, expect, result[key]))
        print('PASS: preflight counts match baseline (%d residues, all matched)'
              % result['n_total'])

        # --- explicit per-residue template assignment ----------------------
        from chimerax.atomic import Residues
        from chimerax.isolde.validation.cmd import _get_forcefield
        from chimerax.isolde.openmm.openmm_interface import find_residue_templates
        ff, ligand_db, ff_name = _get_forcefield(session, None)
        residues = Residues(sorted(m.residues,
            key=lambda r: (r.chain_id, r.number, r.insertion_code)))
        tdict = find_residue_templates(residues, ff, ligand_db=ligand_db,
            logger=session.logger)
        mapping = {}
        for i, tname in tdict.items():
            r = residues[i]
            mapping['/{} {} {}'.format(r.chain_id, r.name, r.number)] = tname
        if mapping != EXPECTED_TEMPLATE_MAP:
            _fail('template map changed.\n  expected %r\n  got      %r'
                  % (EXPECTED_TEMPLATE_MAP, mapping))
        print('PASS: per-residue template map matches baseline (%d explicit)'
              % len(mapping))

        # --- unmatched-residue safety net ----------------------------------
        # OpenMM's template system is the safety net today: an unrecognisable
        # residue is flagged unmatched (the rewrite removes this net). Forcing a
        # heavy atom to boron guarantees no template can match by name/topology.
        from chimerax.atomic import Element
        target = None
        for r in m.residues:
            heavy = r.atoms[r.atoms.element_names != 'H']
            if len(heavy):
                target = (r, heavy[0])
                break
        if target is None:
            _fail('fixture has no heavy atoms?')
        r, atom = target
        atom.element = Element.get_element('B')
        result = isolde_preflight_parameters(session, model=m, forcefield=None)
        if result['ready_for_simulation'] is not False:
            _fail('boron-mutated model still reported ready_for_simulation')
        hit = [u for u in result['unmatched']
               if u['number'] == r.number and u['chain_id'] == r.chain_id]
        if not hit:
            _fail('boron-mutated residue %s %s%d was not flagged unmatched'
                  % (r.name, r.chain_id, r.number))
        if not isinstance(hit[0]['candidates_by_name'], list):
            _fail('unmatched residue lost its candidates_by_name list')
        print('PASS: an unrecognisable residue is flagged unmatched with candidates')
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
