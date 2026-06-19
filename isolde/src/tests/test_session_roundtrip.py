# @Author: Tristan Croll
# @Date:   19-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 19-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Characterise that ISOLDE restraints survive a ChimeraX session save/restore --
the persistence substrate the future model-resident force-field store will
piggyback on. Restraints must be re-queried from the restored manager (not
re-derived from serialised C++ state).
'''

import os
import tempfile

import numpy
import pytest

from chimerax.isolde.tests.conftest import FIXTURE_PDB, cx_path

pytestmark = pytest.mark.integration


def _ca_atoms(model, n=5):
    cas = model.atoms[model.atoms.names == 'CA']
    return cas[:n]


def test_restraints_survive_session_roundtrip(session):
    from chimerax.core.commands import run
    from chimerax.atomic import AtomicStructure
    from chimerax.isolde import session_extensions as sx
    from chimerax.isolde import ISOLDE_STATE_VERSION

    # A bump to the state version is a deliberate event that should force a
    # review of this baseline.
    assert isinstance(ISOLDE_STATE_VERSION, int) and ISOLDE_STATE_VERSION > 0

    model = run(session, 'open "{}"'.format(cx_path(FIXTURE_PDB)))[0]
    cas = _ca_atoms(model, 5)

    pr_mgr = sx.get_position_restraint_mgr(model)
    prs = pr_mgr.add_restraints(cas)
    prs.targets = cas.coords + 1.0
    prs.spring_constants = numpy.linspace(100.0, 500.0, len(prs))
    prs.enableds = True

    dr_mgr = sx.get_distance_restraint_mgr(model)
    dr = dr_mgr.add_restraint(cas[0], cas[-1])
    dr.target = 12.5
    dr.spring_constant = 1234.0
    dr.enabled = True

    pre = {
        'n_pos': len(prs),
        'pos_targets': numpy.array(prs.targets),
        'pos_k': numpy.array(prs.spring_constants),
        'pos_en': numpy.array(prs.enableds),
        'dr_target': dr.target,
        'dr_k': dr.spring_constant,
    }

    tmp = tempfile.NamedTemporaryFile(suffix='.cxs', delete=False)
    tmp.close()
    try:
        run(session, 'save "{}" format session'.format(cx_path(tmp.name)))
        run(session, 'close session')
        run(session, 'open "{}"'.format(cx_path(tmp.name)))

        models = [m for m in session.models.list() if type(m) == AtomicStructure]
        assert len(models) == 1, 'expected exactly one structure after restore'
        m2 = models[0]
        cas2 = _ca_atoms(m2, 5)

        pr_mgr2 = sx.get_position_restraint_mgr(m2, create=False)
        assert pr_mgr2 is not None, 'position restraint manager was not restored'
        prs2 = pr_mgr2.get_restraints(cas2)
        assert len(prs2) == pre['n_pos']
        assert numpy.allclose(prs2.targets, pre['pos_targets'], atol=1e-3)
        assert numpy.allclose(prs2.spring_constants, pre['pos_k'], atol=1e-3)
        assert (numpy.array(prs2.enableds) == pre['pos_en']).all()

        dr_mgr2 = sx.get_distance_restraint_mgr(m2, create=False)
        assert dr_mgr2 is not None, 'distance restraint manager was not restored'
        dr2 = dr_mgr2.get_restraint(cas2[0], cas2[-1])
        assert dr2 is not None, 'distance restraint was not restored'
        assert abs(dr2.target - pre['dr_target']) < 1e-3
        assert abs(dr2.spring_constant - pre['dr_k']) < 1e-3
        assert dr2.enabled
    finally:
        os.unlink(tmp.name)
