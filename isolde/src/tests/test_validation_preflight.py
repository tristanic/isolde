# @Author: Tristan Croll
# @Date:   19-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 19-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Characterise ISOLDE's read-only validation / preflight QC as it behaves today:
the hydrogen preflight, and the Ramachandran / rotamer outlier reporting. These
feed the "quality control and fallbacks" the future rewrite must keep honest.
'''

import pytest

pytestmark = pytest.mark.fast


def test_preflight_hydrogens(session, model, golden):
    from chimerax.isolde.validation.cmd import isolde_preflight_hydrogens
    result = isolde_preflight_hydrogens(session, model=model)
    for key in ('status', 'hydrogen_count', 'heavy_atom_count', 'water_count',
                'recommend_addh'):
        assert key in result
    # Counts are deterministic for the fixture; status is derived from them.
    golden('preflight_hydrogens', {
        'status': result['status'],
        'hydrogen_count': result['hydrogen_count'],
        'heavy_atom_count': result['heavy_atom_count'],
        'water_count': result['water_count'],
        'recommend_addh': result['recommend_addh'],
    })


def test_ramachandran_report(session, model, golden):
    from chimerax.isolde import session_extensions as sx
    mgr = sx.get_ramachandran_mgr(session)
    residues = model.residues
    outliers = mgr.outliers(residues)
    cispeps = mgr.cis(residues)
    twisteds = mgr.twisted(residues)
    golden('ramachandran_report', {
        'n_outliers': len(outliers),
        'n_cis': len(cispeps),
        'n_twisted': len(twisteds),
    })


def test_rotamer_report(session, model, golden):
    from chimerax.isolde import session_extensions as sx
    mgr = sx.get_rotamer_mgr(session)
    rotamers = mgr.get_rotamers(model.residues)
    non_favored, scores = mgr.non_favored_rotamers(rotamers)
    golden('rotamer_report', {
        'n_rotamers': len(rotamers),
        'n_non_favored': len(non_favored),
    })
