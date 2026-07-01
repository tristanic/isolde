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


def _build_handler(session, m):
    from chimerax.isolde.openmm.openmm_interface import SimConstruct, SimHandler
    from chimerax.isolde.openmm.sim_param_mgr import SimParams
    from chimerax.isolde.openmm.forcefields import ForcefieldMgr
    from chimerax.isolde.openmm.param_provider import ForceFieldParameterisationProvider
    residues = m.residues
    ordered = residues[numpy.lexsort((residues.numbers, residues.chain_ids))]
    mobile = ordered[:N_MOBILE_RESIDUES].atoms
    fixed = ordered[N_MOBILE_RESIDUES:N_MOBILE_RESIDUES + N_FIXED_RESIDUES].atoms
    construct = SimConstruct(m, mobile, fixed)
    handler = SimHandler(session, SimParams(), construct,
                         ForceFieldParameterisationProvider(ForcefieldMgr(session)))
    return construct, handler


def run(session):
    import openmm
    from openmm import unit
    from chimerax.core.commands import run as run_cmd
    m = run_cmd(session, 'open "%s"' % FIXTURE.replace('\\', '/'), log=False)[0]
    try:
        construct, handler = _build_handler(session, m)
        system = handler._system

        # --- particle count + mobile/fixed mass partition -----------------
        all_atoms = construct.all_atoms
        if system.getNumParticles() != len(all_atoms):
            _fail('particle count %d != atom count %d'
                  % (system.getNumParticles(), len(all_atoms)))
        fixed_idx = set(all_atoms.indices(construct.fixed_atoms).tolist())
        for i in range(system.getNumParticles()):
            mass = system.getParticleMass(i).value_in_unit(unit.dalton)
            if i in fixed_idx and mass != 0:
                _fail('fixed particle %d not pinned (mass %r)' % (i, mass))
            if i not in fixed_idx and mass <= 0:
                _fail('mobile particle %d lost its mass' % i)
        print('PASS: %d particles; fixed atoms pinned, mobile atoms massive'
              % system.getNumParticles())

        # --- force-group layout + term counts -----------------------------
        force_counts, force_groups = {}, {}
        for f in system.getForces():
            name = type(f).__name__
            force_counts[name] = force_counts.get(name, 0) + 1
            force_groups[name] = f.getForceGroup()
        if force_counts != EXPECTED_FORCE_COUNTS:
            _fail('force counts changed: %r' % force_counts)
        if force_groups != EXPECTED_FORCE_GROUPS:
            _fail('force-group layout changed: %r' % force_groups)
        if system.getNumParticles() != EXPECTED_NUM_PARTICLES:
            _fail('num particles %d != baseline %d'
                  % (system.getNumParticles(), EXPECTED_NUM_PARTICLES))
        print('PASS: force counts + group layout match baseline')

        # --- nonbonded charge + single-point energy (Reference platform) --
        total_charge = 0.0
        for f in system.getForces():
            if isinstance(f, openmm.NonbondedForce):
                for i in range(f.getNumParticles()):
                    q, sigma, eps = f.getParticleParameters(i)
                    total_charge += q.value_in_unit(unit.elementary_charge)
                break
        if abs(round(total_charge, 4) - EXPECTED_TOTAL_CHARGE) > 1e-3:
            _fail('total charge %r != baseline %r'
                  % (round(total_charge, 4), EXPECTED_TOTAL_CHARGE))

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
            _fail('single-point energy %.4f != baseline %.4f (tol %.2f)'
                  % (energy, EXPECTED_ENERGY_KJ_MOL, tol))
        print('PASS: total charge %.1f and single-point energy %.1f kJ/mol match baseline'
              % (total_charge, energy))
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
