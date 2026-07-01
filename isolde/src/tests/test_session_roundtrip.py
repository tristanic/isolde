# @Author: Tristan Croll
# @Date:   30-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 30-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Characterisation test: ISOLDE restraints survive a ChimeraX session
save/restore -- the persistence substrate the future model-resident
force-field store will build on. Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_session_roundtrip.py

Adds position + distance restraints, saves the session, reopens it, and checks
the restraints are re-queryable from the restored managers with their targets /
spring constants / enabled flags intact. Prints PASS/FAIL and exits non-zero on
failure.
'''
import os
import tempfile
import numpy

FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1pmx_1.pdb')


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def _ca_atoms(m, n=5):
    cas = m.atoms[m.atoms.names == 'CA']
    return cas[:n]


def run(session):
    from chimerax.core.commands import run as run_cmd
    from chimerax.atomic import AtomicStructure
    from chimerax.isolde import session_extensions as sx
    from chimerax.isolde import ISOLDE_STATE_VERSION

    if not (isinstance(ISOLDE_STATE_VERSION, int) and ISOLDE_STATE_VERSION > 0):
        _fail('ISOLDE_STATE_VERSION should be a positive int, got %r' % (ISOLDE_STATE_VERSION,))

    m = run_cmd(session, 'open "%s"' % FIXTURE.replace('\\', '/'), log=False)[0]
    cas = _ca_atoms(m, 5)

    pr_mgr = sx.get_position_restraint_mgr(m)
    prs = pr_mgr.add_restraints(cas)
    prs.targets = cas.coords + 1.0
    prs.spring_constants = numpy.linspace(100.0, 500.0, len(prs))
    prs.enableds = True

    dr_mgr = sx.get_distance_restraint_mgr(m)
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
        path = tmp.name.replace('\\', '/')
        run_cmd(session, 'save "%s" format session' % path, log=False)
        run_cmd(session, 'close session', log=False)
        run_cmd(session, 'open "%s"' % path, log=False)

        models = [mm for mm in session.models.list() if type(mm) == AtomicStructure]
        if len(models) != 1:
            _fail('expected exactly one structure after restore, got %d' % len(models))
        m2 = models[0]
        cas2 = _ca_atoms(m2, 5)

        pr_mgr2 = sx.get_position_restraint_mgr(m2, create=False)
        if pr_mgr2 is None:
            _fail('position restraint manager was not restored')
        prs2 = pr_mgr2.get_restraints(cas2)
        if len(prs2) != pre['n_pos']:
            _fail('position restraint count %d != %d' % (len(prs2), pre['n_pos']))
        if not numpy.allclose(prs2.targets, pre['pos_targets'], atol=1e-3):
            _fail('position targets did not survive restore')
        if not numpy.allclose(prs2.spring_constants, pre['pos_k'], atol=1e-3):
            _fail('position spring constants did not survive restore')
        if not (numpy.array(prs2.enableds) == pre['pos_en']).all():
            _fail('position enabled flags did not survive restore')
        print('PASS: %d position restraints survived save/restore intact' % len(prs2))

        dr_mgr2 = sx.get_distance_restraint_mgr(m2, create=False)
        if dr_mgr2 is None:
            _fail('distance restraint manager was not restored')
        dr2 = dr_mgr2.get_restraint(cas2[0], cas2[-1])
        if dr2 is None:
            _fail('distance restraint was not restored')
        if abs(dr2.target - pre['dr_target']) > 1e-3:
            _fail('distance target %r != %r' % (dr2.target, pre['dr_target']))
        if abs(dr2.spring_constant - pre['dr_k']) > 1e-3:
            _fail('distance spring constant %r != %r' % (dr2.spring_constant, pre['dr_k']))
        if not dr2.enabled:
            _fail('distance restraint lost its enabled flag')
        print('PASS: distance restraint survived save/restore intact')
        print('PASS: ISOLDE_STATE_VERSION = %d' % ISOLDE_STATE_VERSION)
    finally:
        os.unlink(tmp.name)

    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
