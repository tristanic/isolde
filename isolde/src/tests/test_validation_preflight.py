# @Author: Tristan Croll
# @Date:   30-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 30-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Characterisation test for ISOLDE's read-only validation / preflight QC: the
hydrogen preflight and the Ramachandran / rotamer outlier reporting. Run inside
ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_validation_preflight.py

Asserts the CURRENT reports for the bundled 1pmx fixture against captured
baselines. Prints PASS/FAIL and exits non-zero on failure.
'''
import os

FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1pmx_1.pdb')

# Baselines captured 2026-06 for 1pmx_1.pdb.
EXPECTED_HYDROGENS = {
    'status': 'ok', 'hydrogen_count': 513, 'heavy_atom_count': 533,
    'water_count': 0, 'recommend_addh': False,
}
EXPECTED_RAMA = {'n_outliers': 0, 'n_cis': 0, 'n_twisted': 0}
EXPECTED_ROTA = {'n_rotamers': 57, 'n_non_favored': 10}


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def run(session):
    from chimerax.core.commands import run as run_cmd
    from chimerax.isolde import session_extensions as sx
    m = run_cmd(session, 'open "%s"' % FIXTURE.replace('\\', '/'), log=False)[0]
    try:
        # --- hydrogen preflight -------------------------------------------
        from chimerax.isolde.validation.cmd import isolde_preflight_hydrogens
        h = isolde_preflight_hydrogens(session, model=m)
        for key, expect in EXPECTED_HYDROGENS.items():
            if h.get(key) != expect:
                _fail('hydrogen preflight %s: expected %r, got %r' % (key, expect, h.get(key)))
        print('PASS: hydrogen preflight matches baseline (%d H / %d heavy, %s)'
              % (h['hydrogen_count'], h['heavy_atom_count'], h['status']))

        # --- Ramachandran report ------------------------------------------
        rama = sx.get_ramachandran_mgr(session)
        residues = m.residues
        got_rama = {
            'n_outliers': len(rama.outliers(residues)),
            'n_cis': len(rama.cis(residues)),
            'n_twisted': len(rama.twisted(residues)),
        }
        if got_rama != EXPECTED_RAMA:
            _fail('Ramachandran report changed: %r' % got_rama)
        print('PASS: Ramachandran report matches baseline %r' % got_rama)

        # --- rotamer report -----------------------------------------------
        rota = sx.get_rotamer_mgr(session)
        rotamers = rota.get_rotamers(residues)
        non_favored, scores = rota.non_favored_rotamers(rotamers)
        got_rota = {'n_rotamers': len(rotamers), 'n_non_favored': len(non_favored)}
        if got_rota != EXPECTED_ROTA:
            _fail('rotamer report changed: %r' % got_rota)
        print('PASS: rotamer report matches baseline %r' % got_rota)
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
