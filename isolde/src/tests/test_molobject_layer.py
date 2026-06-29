# @Author: Tristan Croll
# @Date:   19-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 19-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Characterise the C-backed molecular-object layer that the future model-resident
force-field store will build on:

  * the per-object <-> vectorized-Collection parity (a singular ``c_property``
    and its plural ``cvec_property`` twin must read the same values), and
  * the custom attributes ISOLDE registers on core ChimeraX objects so they
    survive session save/restore.

The parameter store will add new per-object/Collection property pairs (and a new
Angle class) following exactly this pattern, so the invariant must hold today.
'''

import numpy
import pytest

pytestmark = pytest.mark.fast


def _assert_parity(collection, singular_attr, plural_attr):
    '''The plural cvec_property must equal the element-wise singular reads.'''
    assert len(collection) > 0, 'no objects to compare'
    plural = numpy.asarray(getattr(collection, plural_attr))
    singular = numpy.array([getattr(obj, singular_attr) for obj in collection])
    assert numpy.allclose(plural, singular, equal_nan=True), \
        '{} != per-object {}'.format(plural_attr, singular_attr)


def test_chiral_center_parity(session, model):
    from chimerax.isolde.molobject import get_chiral_mgr
    cm = get_chiral_mgr(session)
    chirals = cm.get_chirals(model.atoms, create=True)
    _assert_parity(chirals, 'angle', 'angles')


def test_proper_dihedral_parity(session, model):
    from chimerax.isolde.session_extensions import get_proper_dihedral_mgr
    pdm = get_proper_dihedral_mgr(session)
    phis = pdm.get_dihedrals(model.residues, 'phi', create=True)
    _assert_parity(phis, 'angle', 'angles')


def test_isolde_custom_attrs_round_trip(session, model):
    '''The attributes ISOLDE registers (src/__init__.py) are usable on live
    objects -- registration is what lets them persist in sessions.'''
    r = model.residues[0]
    original = getattr(r, 'isolde_ignore', False)
    try:
        r.isolde_ignore = True
        assert r.isolde_ignore is True
        r.isolde_ignore = False
        assert r.isolde_ignore is False
    finally:
        r.isolde_ignore = bool(original)

    # AtomicStructure-level registered attr is settable/readable too.
    model.isolde_initialized = True
    assert model.isolde_initialized is True
