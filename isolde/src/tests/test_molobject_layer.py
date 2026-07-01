# @Author: Tristan Croll
# @Date:   30-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 30-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Characterisation test for the C-backed molecular-object layer the future
model-resident force-field store will build on. Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_molobject_layer.py

Checks the per-object <-> vectorized-Collection parity (a singular c_property
and its plural cvec_property twin must read the same values -- the invariant new
parameter properties, and a new Angle class, will follow), and that the custom
attributes ISOLDE registers on core ChimeraX objects are usable. Prints PASS/FAIL
and exits non-zero on failure.
'''
import os
import numpy

FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1pmx_1.pdb')


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def _check_parity(collection, singular_attr, plural_attr, label):
    if len(collection) == 0:
        _fail('%s: no objects to compare' % label)
    plural = numpy.asarray(getattr(collection, plural_attr))
    singular = numpy.array([getattr(obj, singular_attr) for obj in collection])
    if not numpy.allclose(plural, singular, equal_nan=True):
        _fail('%s: %s (plural) != per-object %s' % (label, plural_attr, singular_attr))


def run(session):
    from chimerax.core.commands import run as run_cmd
    m = run_cmd(session, 'open "%s"' % FIXTURE.replace('\\', '/'), log=False)[0]
    try:
        # --- chiral centre parity -----------------------------------------
        from chimerax.isolde.molobject import get_chiral_mgr
        chirals = get_chiral_mgr(session).get_chirals(m.atoms, create=True)
        _check_parity(chirals, 'angle', 'angles', 'ChiralCenters')
        print('PASS: ChiralCenters.angles matches per-object .angle (%d centres)'
              % len(chirals))

        # --- proper-dihedral parity ---------------------------------------
        from chimerax.isolde.session_extensions import get_proper_dihedral_mgr
        phis = get_proper_dihedral_mgr(session).get_dihedrals(m.residues, 'phi', create=True)
        _check_parity(phis, 'angle', 'angles', 'ProperDihedrals')
        print('PASS: ProperDihedrals.angles matches per-object .angle (%d phi)' % len(phis))

        # --- custom-attr round-trip ---------------------------------------
        r = m.residues[0]
        original = getattr(r, 'isolde_ignore', False)
        r.isolde_ignore = True
        if r.isolde_ignore is not True:
            _fail('Residue.isolde_ignore did not round-trip True')
        r.isolde_ignore = bool(original)
        m.isolde_initialized = True
        if m.isolde_initialized is not True:
            _fail('AtomicStructure.isolde_initialized did not round-trip True')
        print('PASS: ISOLDE custom attributes are registered and round-trip')
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
