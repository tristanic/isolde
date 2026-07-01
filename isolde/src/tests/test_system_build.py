# @Author: Tristan Croll
# @Date:   30-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 30-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Characterisation test for the OpenMM ``System`` ISOLDE builds from forcefield
templates -- the equivalence baseline the future "direct System build" must
reproduce. Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_system_build.py

Builds a small SimConstruct + SimHandler (no simulation start), then checks the
force-group layout, particle/term counts, fixed-atom masses, the nonbonded
charge, and a single-point energy on the deterministic Reference platform
against captured baselines. Prints PASS/FAIL and exits non-zero on failure.
'''
import os
import numpy

FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1pmx_1.pdb')

N_MOBILE_RESIDUES = 8
N_FIXED_RESIDUES = 6

# Baseline captured 2026-06 (amber14 + bundled XMLs, Reference platform).
EXPECTED_FORCE_COUNTS = {
    'HarmonicBondForce': 1, 'HarmonicAngleForce': 1,
    'PeriodicTorsionForce': 1, 'NonbondedForce': 1,
}
EXPECTED_FORCE_GROUPS = {
    'HarmonicBondForce': 1, 'HarmonicAngleForce': 2,
    'PeriodicTorsionForce': 3, 'NonbondedForce': 0,
}
EXPECTED_NUM_PARTICLES = 189
EXPECTED_TOTAL_CHARGE = -2.0
EXPECTED_ENERGY_KJ_MOL = 1359.25532406659


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def _build_handler(session, m, provider):
    from chimerax.isolde.openmm.openmm_interface import SimConstruct, SimHandler
    from chimerax.isolde.openmm.sim_param_mgr import SimParams
    residues = m.residues
    ordered = residues[numpy.lexsort((residues.numbers, residues.chain_ids))]
    mobile = ordered[:N_MOBILE_RESIDUES].atoms
    fixed = ordered[N_MOBILE_RESIDUES:N_MOBILE_RESIDUES + N_FIXED_RESIDUES].atoms
    construct = SimConstruct(m, mobile, fixed)
    handler = SimHandler(session, SimParams(), construct, provider)
    return construct, handler


def _check_baseline(system, construct, label):
    import openmm
    from openmm import unit

    # --- particle count + mobile/fixed mass partition -----------------
    all_atoms = construct.all_atoms
    if system.getNumParticles() != len(all_atoms):
        _fail('%s: particle count %d != atom count %d'
              % (label, system.getNumParticles(), len(all_atoms)))
    fixed_idx = set(all_atoms.indices(construct.fixed_atoms).tolist())
    for i in range(system.getNumParticles()):
        mass = system.getParticleMass(i).value_in_unit(unit.dalton)
        if i in fixed_idx and mass != 0:
            _fail('%s: fixed particle %d not pinned (mass %r)' % (label, i, mass))
        if i not in fixed_idx and mass <= 0:
            _fail('%s: mobile particle %d lost its mass' % (label, i))
    print('PASS (%s): %d particles; fixed atoms pinned, mobile atoms massive'
          % (label, system.getNumParticles()))

    # --- force-group layout + term counts -----------------------------
    force_counts, force_groups = {}, {}
    for f in system.getForces():
        name = type(f).__name__
        force_counts[name] = force_counts.get(name, 0) + 1
        force_groups[name] = f.getForceGroup()
    if force_counts != EXPECTED_FORCE_COUNTS:
        _fail('%s: force counts changed: %r' % (label, force_counts))
    if force_groups != EXPECTED_FORCE_GROUPS:
        _fail('%s: force-group layout changed: %r' % (label, force_groups))
    if system.getNumParticles() != EXPECTED_NUM_PARTICLES:
        _fail('%s: num particles %d != baseline %d'
              % (label, system.getNumParticles(), EXPECTED_NUM_PARTICLES))
    print('PASS (%s): force counts + group layout match baseline' % label)

    # --- nonbonded charge + single-point energy (Reference platform) --
    total_charge = 0.0
    for f in system.getForces():
        if isinstance(f, openmm.NonbondedForce):
            for i in range(f.getNumParticles()):
                q, sigma, eps = f.getParticleParameters(i)
                total_charge += q.value_in_unit(unit.elementary_charge)
            break
    if abs(round(total_charge, 4) - EXPECTED_TOTAL_CHARGE) > 1e-3:
        _fail('%s: total charge %r != baseline %r'
              % (label, round(total_charge, 4), EXPECTED_TOTAL_CHARGE))

    box = 50.0 * unit.nanometer
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(box, 0, 0), openmm.Vec3(0, box, 0), openmm.Vec3(0, 0, box))
    platform = openmm.Platform.getPlatformByName('Reference')
    context = openmm.Context(system, openmm.VerletIntegrator(1.0 * unit.femtoseconds),
                             platform)
    context.setPositions((construct.all_atoms.coords / 10.0) * unit.nanometer)
    energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)
    del context
    tol = 1.0 + 1e-3 * abs(EXPECTED_ENERGY_KJ_MOL)
    if abs(energy - EXPECTED_ENERGY_KJ_MOL) > tol:
        _fail('%s: single-point energy %.4f != baseline %.4f (tol %.2f)'
              % (label, energy, EXPECTED_ENERGY_KJ_MOL, tol))
    print('PASS (%s): total charge %.1f and single-point energy %.1f kJ/mol match baseline'
          % (label, total_charge, energy))


def run(session):
    from chimerax.core.commands import run as run_cmd
    from chimerax.isolde.openmm.forcefields import ForcefieldMgr
    from chimerax.isolde.openmm.param_provider import ForceFieldParameterisationProvider
    m = run_cmd(session, 'open "%s"' % FIXTURE.replace('\\', '/'), log=False)[0]
    try:
        provider = ForceFieldParameterisationProvider(ForcefieldMgr(session))

        # --- first build: cold cache, every residue must be a cache miss --
        construct, handler = _build_handler(session, m, provider)
        _check_baseline(handler._system, construct, 'cold cache')
        stats = provider.last_build_stats
        if stats['residues_cache_miss'] != stats['residues_total']:
            _fail('cold-cache build: expected every residue to miss, got %r' % stats)
        print('PASS: cold-cache build resolved all %d residues via template matching'
              % stats['residues_total'])

        # --- second build: warm cache, every residue must be a cache hit,
        #     and the resulting system must reproduce the same baseline -----
        construct2, handler2 = _build_handler(session, m, provider)
        stats2 = provider.last_build_stats
        if stats2['residues_cache_miss'] != 0 or stats2['residues_cache_hit'] != stats2['residues_total']:
            _fail('warm-cache build: expected every residue to hit, got %r' % stats2)
        print('PASS: warm-cache build resolved all %d residues from cache (0 misses)'
              % stats2['residues_total'])
        _check_baseline(handler2._system, construct2, 'warm cache')
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
