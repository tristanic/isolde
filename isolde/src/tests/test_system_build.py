# @Author: Tristan Croll
# @Date:   19-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 19-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll

'''
Characterise the *current* OpenMM ``System`` that ISOLDE builds from forcefield
templates. This is the equivalence baseline the future "direct System build"
must reproduce: force-group layout, particle/term counts, fixed-atom masses,
the nonbonded charge digest, and a single-point energy on a deterministic
CPU-only platform.
'''

import numpy
import pytest

pytestmark = pytest.mark.integration

# A small, deterministic partition of the fixture so the System stays tiny.
N_MOBILE_RESIDUES = 8
N_FIXED_RESIDUES = 6


def _build_handler(session, model):
    from chimerax.isolde.openmm.openmm_interface import SimConstruct, SimHandler
    from chimerax.isolde.openmm.sim_param_mgr import SimParams
    from chimerax.isolde.openmm.forcefields import ForcefieldMgr

    residues = model.residues
    ordered = residues[numpy.lexsort((residues.numbers, residues.chain_ids))]
    mobile = ordered[:N_MOBILE_RESIDUES].atoms
    fixed = ordered[N_MOBILE_RESIDUES:N_MOBILE_RESIDUES + N_FIXED_RESIDUES].atoms

    construct = SimConstruct(model, mobile, fixed)
    handler = SimHandler(session, SimParams(), construct, ForcefieldMgr(session))
    return construct, handler


def test_particle_count_and_fixed_masses(session, model):
    '''Particle count tracks the simulation atoms, and fixed atoms (and only
    those) are pinned with zero mass -- ISOLDE's mobile/fixed partition.'''
    from openmm import unit
    construct, handler = _build_handler(session, model)
    system = handler._system

    all_atoms = construct.all_atoms
    assert system.getNumParticles() == len(all_atoms)

    fixed_idx = set(all_atoms.indices(construct.fixed_atoms).tolist())
    for i in range(system.getNumParticles()):
        mass = system.getParticleMass(i).value_in_unit(unit.dalton)
        if i in fixed_idx:
            assert mass == 0, 'fixed particle {} not pinned'.format(i)
        else:
            assert mass > 0, 'mobile particle {} lost its mass'.format(i)


def test_system_digest(session, model, cpu_platform, golden):
    '''Golden digest of the built System + a single-point potential energy.'''
    import openmm
    from openmm import unit
    construct, handler = _build_handler(session, model)
    system = handler._system

    force_counts = {}
    force_groups = {}
    for f in system.getForces():
        name = type(f).__name__
        force_counts[name] = force_counts.get(name, 0) + 1
        force_groups[name] = f.getForceGroup()

    # Total nonbonded charge (rounded) -- a cheap fingerprint of the parameters.
    total_charge = 0.0
    for f in system.getForces():
        if isinstance(f, openmm.NonbondedForce):
            for i in range(f.getNumParticles()):
                q, sigma, eps = f.getParticleParameters(i)
                total_charge += q.value_in_unit(unit.elementary_charge)
            break

    # Deterministic single-point energy. A large box keeps any periodic
    # nonbonded method well-defined without perturbing a small fragment.
    box = 50.0 * unit.nanometer
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(box, 0, 0), openmm.Vec3(0, box, 0), openmm.Vec3(0, 0, box))
    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    context = openmm.Context(system, integrator, cpu_platform)
    context.setPositions((construct.all_atoms.coords / 10.0) * unit.nanometer)
    energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)
    del context, integrator

    golden('system_build_digest', {
        'num_particles': system.getNumParticles(),
        'force_counts': force_counts,
        'force_groups': force_groups,
        'total_charge': round(total_charge, 4),
        'potential_energy_kj_mol': energy,
    }, rtol=1e-3, atol=1.0)
