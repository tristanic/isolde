# @Author: Tristan Croll
# @Date:   29-Oct-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 29-Oct-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Regression test for the variable-timestep integrator set-up
(``SimHandler._prepare_integrator``): the maximum step size, and the tolerance
correction for zero-mass crystallographic symmetry copies.

OpenMM's variable-step integrators pick ``dt = sqrt(tol / err)``, with ``err`` the RMS
of force/mass over *every* particle, and may double the step each step. Two problems
follow, both seen on COD 2224840 (garnet, symmetry-aware, flexible hydrogens):

  * symmetry copies are zero-mass virtual sites, which count in the mean but add
    nothing to it, loosening the effective tolerance by sqrt(n_total / n_real);
  * the step follows the forces, not the stiffness, so near a minimum it grows until
    fully flexible C-H bonds go unstable and the hydrogens "bounce".

ISOLDE keeps the variable integrator because tugging forces would explode a fixed
step; the cap only limits the step from above, so large forces still shrink it.

Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_variable_integrator.py

Prints PASS/FAIL and exits non-zero on failure.
'''
from types import SimpleNamespace
import numpy


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def test_profile_defaults():
    from chimerax.isolde.openmm.sim_param_mgr import SimParams
    from chimerax.isolde.openmm.forcefields import _GARNET_VARIANTS
    sp = SimParams(); sp.forcefield = 'amber14'
    if sp.variable_integrator_max_step != 0.0:
        _fail('AMBER should keep an uncapped variable step, got %r' % sp.variable_integrator_max_step)
    for ff in _GARNET_VARIANTS:
        sp = SimParams(); sp.forcefield = ff
        if sp.variable_integrator_max_step != 0.002 or sp.variable_integrator_tolerance != 1e-5:
            _fail('%s: max step %r, tolerance %r (expected 0.002 ps, 1e-5)'
                  % (ff, sp.variable_integrator_max_step, sp.variable_integrator_tolerance))
    print('PASS: garnet variants cap the step at 2 fs with tolerance 1e-5; AMBER uncapped')


def test_tolerance_scaling():
    from chimerax.isolde.openmm.openmm_interface import symmetry_scaled_tolerance as f
    if f(1e-5, 46, 46) != 1e-5 or f(1e-5, None, 2532) != 1e-5:
        _fail('tolerance changed without symmetry copies')
    got = f(1e-5, 46, 2532)
    if abs(got - 1e-5 * (46 / 2532) ** 0.5) > 1e-18:
        _fail('scaled tolerance %r' % got)
    print('PASS: tolerance scaled by sqrt(n_real/n_total) only when copies are present')


def test_prepare_integrator_wiring():
    import openmm
    from chimerax.isolde.openmm.sim_param_mgr import SimParams
    from chimerax.isolde.openmm.openmm_interface import SimHandler
    system = openmm.System()
    for _ in range(46):
        system.addParticle(12.0)
    for _ in range(2486):
        system.addParticle(0.0)
    logged = []
    fake = SimpleNamespace(temperature=0.0, _system=system, _num_real_atoms=46,
        session=SimpleNamespace(logger=SimpleNamespace(info=logged.append)))
    sp = SimParams(); sp.forcefield = 'garnet-r11'
    integ = SimHandler._prepare_integrator(fake, sp).getIntegrator(0)
    max_step = integ.getMaximumStepSize()
    max_step = max_step.value_in_unit(openmm.unit.picosecond) if hasattr(max_step, 'unit') else max_step
    if abs(max_step - 0.002) > 1e-12:
        _fail('integrator max step %r (expected 0.002 ps)' % max_step)
    want = 1e-5 * (46 / 2532) ** 0.5
    if abs(integ.getErrorTolerance() - want) > 1e-15:
        _fail('integrator tolerance %r (expected %r)' % (integ.getErrorTolerance(), want))
    if not logged:
        _fail('tolerance scaling was not logged')
    print('PASS: _prepare_integrator applies the 2 fs cap and the scaled tolerance')


def _toy(n_vsites):
    '''23 flexible C-H bonds at HMR-3 masses (H 3 Da, C 10 Da), 2% stretched, plus
    zero-mass virtual sites standing in for symmetry copies and an off-by-default tug.'''
    import openmm
    K, R0, NPAIR = 2.84e5, 0.109, 23
    s = openmm.System(); hb = openmm.HarmonicBondForce(); pos = []
    for p in range(NPAIR):
        c = s.addParticle(10.0); h = s.addParticle(3.0); hb.addBond(c, h, R0, K)
        b = numpy.array([p * 0.5, 0.0, 0.0]); pos += [b, b + [R0 * 1.02, 0, 0]]
    s.addForce(hb)
    for v in range(n_vsites):
        i = s.addParticle(0.0); a = 2 * (v % NPAIR)
        s.setVirtualSite(i, openmm.TwoParticleAverageSite(a, a + 1, 0.5, 0.5))
        pos.append(numpy.zeros(3))
    tug = openmm.CustomExternalForce('0.5*ktug*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    tug.addGlobalParameter('ktug', 0.0)
    for n in ('x0', 'y0', 'z0'):
        tug.addPerParticleParameter(n)
    tug.addParticle(0, [0.0, 1.0, 0.0])
    s.addForce(tug)
    return s, numpy.array(pos), R0, NPAIR


def _run(n_vsites, tol, max_step, n_steps=2500, tug_at=None):
    import openmm
    from openmm import unit
    s, pos, R0, NPAIR = _toy(n_vsites)
    integ = openmm.VariableLangevinIntegrator(0 * unit.kelvin, 40 / unit.picosecond, tol)
    if max_step:
        integ.setMaximumStepSize(max_step)
    ctx = openmm.Context(s, integ, openmm.Platform.getPlatformByName('Reference'))
    ctx.setPositions(pos); ctx.computeVirtualSites()
    late_dev, tug_dts = 0.0, []
    for k in range(n_steps):
        if tug_at is not None and k == tug_at:
            ctx.setParameter('ktug', 5e4)
        integ.step(1)
        if tug_at is not None and k >= tug_at:
            tug_dts.append(integ.getStepSize().value_in_unit(unit.femtosecond))
        if k >= n_steps - 500 and k % 10 == 0:
            x = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            d = numpy.linalg.norm(x[1:2 * NPAIR:2] - x[0:2 * NPAIR:2], axis=1)
            late_dev = max(late_dev, float(numpy.abs(d - R0).max()))
    return late_dev, (min(tug_dts) if tug_dts else None)


def test_openmm_behaviour():
    from chimerax.isolde.openmm.openmm_interface import symmetry_scaled_tolerance
    nv = 2486
    bounce, _ = _run(nv, 1e-4, None)
    if bounce < 0.005:
        _fail('uncapped integrator did not bounce (late C-H deviation %.4f nm) - the toy no '
              'longer reproduces the failure, so this test is not exercising anything' % bounce)
    fixed, _ = _run(nv, symmetry_scaled_tolerance(1e-4, 46, 46 + nv), 0.002)
    if fixed > 1e-3:
        _fail('capped + scaled integrator still bounces (late C-H deviation %.4f nm)' % fixed)
    _, dt_uncapped = _run(nv, 1e-5, None, n_steps=2000, tug_at=1000)
    _, dt_capped = _run(nv, 1e-5, 0.002, n_steps=2000, tug_at=1000)
    # The two runs take different step histories before the tug, so compare the shrink, not
    # the exact step: under the tug the step must drop far below the 2 fs cap, as uncapped.
    if dt_capped > 1.0 or abs(dt_capped - dt_uncapped) > 0.1 * dt_uncapped:
        _fail('under a tug the capped step only shrank to %.4f fs (uncapped %.4f fs)'
              % (dt_capped, dt_uncapped))
    print('PASS: flexible C-H with %d virtual sites bounces uncapped (%.4f nm), settles with the '
          'fixes (%.5f nm); under a tug the step still shrinks to %.3f fs either way'
          % (nv, bounce, fixed, dt_capped))


def run(session=None):
    test_profile_defaults()
    test_tolerance_scaling()
    test_prepare_integrator_wiring()
    test_openmm_behaviour()
    print('ALL PASS')


try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
