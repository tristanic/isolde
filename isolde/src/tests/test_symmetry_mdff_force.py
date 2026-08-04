# @Author: Tristan Croll
# @Date:   16-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 16-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Regression test for the symmetry-aware MDFF map force
(``SymmetryAwareCubicInterpMapForce`` in ``openmm.custom_forces``), which lets a
real atom feel the map through a crystallographic operator baked into the energy
expression. Builds tiny OpenMM contexts on the Reference platform (no GL/GUI, no
Clipper maps) and checks:

  1. an IDENTITY per-bond transform reproduces the base ``CubicInterpMapForce``
     energy AND force exactly (guards the energy-string rewrite + parameter
     packing, and the bit-for-bit non-symmetry guarantee);
  2. a transformed term samples the map at ``S.r`` and folds the force back onto
     the real atom as ``R^T . grad(map)`` (the SymmetrySite fold-back, computed
     directly by the expression's chain rule);
  3. the analytic force matches a finite-difference of the energy.

Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_symmetry_mdff_force.py

Prints PASS/FAIL and exits non-zero on failure.
'''
import numpy


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


# A fixed pseudo-random map + a simple xyz(nm)->ijk transform placing a ~0.x nm
# point mid-grid (well clear of the box edges for the tricubic stencil).
_DIM = 16
_DATA = numpy.random.RandomState(42).rand(_DIM, _DIM, _DIM).astype(numpy.float32)
_TF = numpy.zeros((3, 4))
for _a in range(3):
    _TF[_a, _a] = 10.0
    _TF[_a, 3] = _DIM / 2.0
_IDENT12 = numpy.array([1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0], dtype=numpy.float64)


def _eval_one(force_cls, pos_nm, transform12=None):
    '''Build a 1-particle System with a single map term; return (energy, force).'''
    from openmm import System, Context, VerletIntegrator, Platform, unit
    sysm = System()
    sysm.addParticle(1.0)
    f = force_cls(_DATA.copy(), _TF.copy(), 'x', units='nanometers', map_sigma=1.0)
    idx = numpy.array([0], numpy.int32)
    if transform12 is None:
        f.add_atoms(idx, numpy.array([1.0]), numpy.array([1.0]))
    else:
        f.add_atoms(idx, numpy.array([1.0]), numpy.array([1.0]),
                    transforms=numpy.array([transform12], numpy.float64))
    f.set_global_k(1.0)
    sysm.addForce(f)
    ctx = Context(sysm, VerletIntegrator(1.0),
                  Platform.getPlatformByName('Reference'))
    ctx.setPositions(numpy.asarray(pos_nm, dtype=numpy.float64).reshape(1, 3))
    st = ctx.getState(getEnergy=True, getForces=True)
    e = st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    fr = st.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer)[0]
    return e, numpy.asarray(fr)


def run(session=None):
    from chimerax.isolde.openmm.custom_forces import (
        CubicInterpMapForce, SymmetryAwareCubicInterpMapForce)

    p = numpy.array([0.30, 0.42, 0.55])

    # (1) identity transform == base force (energy AND force), to machine precision.
    e_base, f_base = _eval_one(CubicInterpMapForce, p)
    e_id, f_id = _eval_one(SymmetryAwareCubicInterpMapForce, p, _IDENT12)
    if abs(e_base - e_id) > 1e-6 * (1 + abs(e_base)):
        _fail('identity transform energy %r != base %r' % (e_id, e_base))
    if not numpy.allclose(f_base, f_id, atol=1e-5, rtol=1e-5):
        _fail('identity transform force %r != base %r' % (f_id, f_base))
    print('PASS: identity transform reproduces base force (E rel %.1e, F max %.1e)'
          % (abs(e_base - e_id) / (1 + abs(e_base)),
             numpy.max(numpy.abs(f_base - f_id))))

    # (2) transformed term samples S.p and folds back as R^T . grad.
    theta = 0.7
    c, s = numpy.cos(theta), numpy.sin(theta)
    R = numpy.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])  # rot about z
    t = numpy.array([0.05, -0.03, 0.02])                           # nm
    S12 = numpy.concatenate([R.ravel(), t])
    q = R @ p + t                                                  # ghost position
    e_sym, f_sym = _eval_one(SymmetryAwareCubicInterpMapForce, p, S12)
    e_bq, f_bq = _eval_one(CubicInterpMapForce, q)
    if abs(e_sym - e_bq) > 1e-5 * (1 + abs(e_bq)):
        _fail('transformed term energy %r != base at S.p %r' % (e_sym, e_bq))
    f_expected = R.T @ f_bq
    if not numpy.allclose(f_sym, f_expected, atol=1e-4, rtol=1e-4):
        _fail('transformed term force %r != R^T.grad %r' % (f_sym, f_expected))
    print('PASS: transformed term samples S.r and folds back R^T.grad '
          '(E rel %.1e, F max %.1e)'
          % (abs(e_sym - e_bq) / (1 + abs(e_bq)),
             numpy.max(numpy.abs(f_sym - f_expected))))

    # (3) analytic force == finite-difference on the transformed term.
    h = 1e-5
    fd = numpy.zeros(3)
    for a in range(3):
        pp = p.copy(); pp[a] += h
        pm = p.copy(); pm[a] -= h
        ep, _ = _eval_one(SymmetryAwareCubicInterpMapForce, pp, S12)
        em, _ = _eval_one(SymmetryAwareCubicInterpMapForce, pm, S12)
        fd[a] = -(ep - em) / (2 * h)
    if not numpy.allclose(f_sym, fd, atol=1e-2, rtol=1e-2):
        _fail('analytic force %r != finite-difference %r' % (f_sym, fd))
    print('PASS: analytic force == finite-difference on transformed term (max %.1e)'
          % numpy.max(numpy.abs(f_sym - fd)))

    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
