# @Author: Tristan Croll
# @Date:   24-Aug-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 24-Aug-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Build a runnable OpenMM ``System`` for an ISOLDE simulation from garnet
parameters, reproducing ISOLDE's real simulation contract (not garnet's vacuum
defaults).

Force *expressions* are imported verbatim from ``garnet_core.openmm_build`` so
the functional form has exactly one source of truth; only the *assembly*
(subsetting to ``all_atoms``, particle-index remapping, cutoff method, force
groups, constraints) is ISOLDE-specific. Every force's group is assigned at
add-time (garnet's torsion/OOP forces are ``CustomCompoundBondForce`` and would
be invisible to the type-keyed ``_set_core_force_groups``).

The returned ``System`` has particles in ``sim_construct.all_atoms`` order and
carries no fixed atoms yet — ``SimHandler.set_fixed_atoms`` zeroes the buffer
masses afterward, exactly as on the AMBER path.
'''
import math

import numpy
from openmm import openmm, unit


def build_garnet_system(session, sim_params, sim_construct, garnet_params):
    '''
    Assemble ``(topology, system)`` for ``sim_construct`` from garnet parameters.

    Args:
        * session: ChimeraX session (for logging).
        * sim_params: the :class:`SimParams` (cutoff method/distance, rigid_bonds,
          rigid_water are honoured; use_gbsa is handled downstream, unchanged).
        * sim_construct: the :class:`SimConstruct` (defines ``all_atoms`` order).
        * garnet_params: a parameterised :class:`GarnetParameters`.
    '''
    # Deferred imports break the openmm_interface <-> garnet import cycle and reuse
    # the exact topology builder + force-group constants the AMBER path uses.
    from ..openmm_interface import (
        create_openmm_topology,
        BOND_FORCE_GROUP, ANGLE_FORCE_GROUP, DIHEDRAL_FORCE_GROUP, NONBONDED_FORCE_GROUP,
    )
    # garnet's canonical force expressions / builders (single source of truth). The 1-4 and
    # gated-1-3 double-exponential CustomBondForces are built by garnet's own make_*_vdw_force,
    # fed the subset bond graph in particle-index space (so the pairs are computed on the SAME
    # graph the exclusions use — no double-counting).
    from types import SimpleNamespace
    from garnet_core.openmm_build import (_DEXP_ENERGY, _OOP_ENERGY,
                                          make_14_vdw_force, make_13_vdw_force)
    from garnet_core.params import TORSION_PERIODICITIES, TORSION_PHASES
    from garnet_core.energy import TORSION_DAMP_THETA_ON

    logger = session.logger
    all_atoms = sim_construct.all_atoms
    subset = garnet_params.subset_for(all_atoms)
    if subset.missing_atoms:
        logger.warning('garnet: {} atom(s) in the simulation had no cached parameters and were '
                       'assigned zero charge/vdW. This usually means the model changed since '
                       'parameterisation.'.format(len(subset.missing_atoms)))

    top, _ = create_openmm_topology(all_atoms, {})

    system = openmm.System()
    for mass in all_atoms.elements.masses:            # amu
        system.addParticle(float(mass))

    # --- Harmonic bonds ---
    if subset.bonds:
        bf = openmm.HarmonicBondForce()
        for i, j, k, r0 in subset.bonds:
            bf.addBond(i, j, r0, k)                   # r0 nm, k kJ/mol/nm^2
        bf.setForceGroup(BOND_FORCE_GROUP)
        system.addForce(bf)

    # --- Harmonic angles ---
    if subset.angles:
        af = openmm.HarmonicAngleForce()
        for i, j, k, k_theta, theta0 in subset.angles:
            af.addAngle(i, j, k, theta0, k_theta)     # theta0 rad, k kJ/mol/rad^2
        af.setForceGroup(ANGLE_FORCE_GROUP)
        system.addForce(af)

    # --- Proper torsions: angle-damped, via CustomCompoundBondForce ---
    # Expression mirrors garnet_core.openmm_build.build_system (TC ratified 2026-08-17):
    # each torsion is damped by a C1 switch on its own two bend angles so the Hessian stays
    # bounded as a bend approaches 180 deg (a user can drag an atom through linear geometry in
    # interactive MD). The switch constants are IMPORTED, so only the string template is local;
    # the single-point cross-check against build_system guards against drift.
    if subset.propers:
        terms = ' + '.join(
            f'k{t + 1}*(1+cos({TORSION_PERIODICITIES[t]}*phi-({TORSION_PHASES[t]:.17g})))'
            for t in range(6))
        expr = (f'({terms})*s1*s2;'
                f's1=uu1*uu1*(3-2*uu1); uu1=min(1,(PI-angle(p1,p2,p3))/(PI-TON));'
                f's2=uu2*uu2*(3-2*uu2); uu2=min(1,(PI-angle(p2,p3,p4))/(PI-TON));'
                f'phi=dihedral(p1,p2,p3,p4);'
                f'PI={math.pi:.17g}; TON={TORSION_DAMP_THETA_ON:.17g}')
        tf = openmm.CustomCompoundBondForce(4, expr)
        for t in range(6):
            tf.addPerBondParameter('k{}'.format(t + 1))
        for i, j, k, l, ks in subset.propers:
            if all(abs(x) < 1e-10 for x in ks):
                continue
            tf.addBond([i, j, k, l], [float(x) for x in ks])
        if tf.getNumBonds():
            tf.setForceGroup(DIHEDRAL_FORCE_GROUP)
            system.addForce(tf)

    # --- Impropers: harmonic OOP if the model predicts oop_k, else periodic improper ---
    if subset.impropers:
        if subset.has_oop:
            oop = openmm.CustomCompoundBondForce(4, _OOP_ENERGY)
            oop.addPerBondParameter('k')
            for kind, (i, j, k, l), payload in subset.impropers:
                oop.addBond([i, j, k, l], [float(payload)])   # (central, o1, o2, o3)
            if oop.getNumBonds():
                oop.setForceGroup(DIHEDRAL_FORCE_GROUP)
                system.addForce(oop)
        else:
            itf = openmm.PeriodicTorsionForce()
            for kind, (c, o1, o2, o3), payload in subset.impropers:
                for t in range(2):
                    kt = float(payload[t])
                    if abs(kt) < 1e-10:
                        continue
                    # released-style periodic improper: central atom THIRD
                    itf.addTorsion(o1, o2, c, o3, TORSION_PERIODICITIES[t], TORSION_PHASES[t], kt)
            if itf.getNumTorsions():
                itf.setForceGroup(DIHEDRAL_FORCE_GROUP)
                system.addForce(itf)

    # Bonded pairs (particle-index space) for exclusions/exceptions == intra_bonds of all_atoms.
    bond_pairs = [(i, j) for i, j, _, _ in subset.bonds]

    # Nonbonded method + cutoff from ISOLDE's contract (CutoffNonPeriodic, ~1.7 nm; no box).
    # sim_params stores the openmm.app method object (what createSystem consumes); translate it
    # to the low-level NonbondedForce / CustomNonbondedForce enums that setNonbondedMethod wants.
    nb_method, cnb_method = _nb_methods(sim_params.nonbonded_cutoff_method)
    cutoff = sim_params.nonbonded_cutoff_distance
    g = subset.globals

    # --- Coulomb (NonbondedForce, LJ zeroed): 1-2/1-3 excluded, 1-4 scaled by coulomb14scale ---
    # GBSA (if enabled) locates THIS force by isinstance(NonbondedForce) and reads charges back out
    # of it; the double-exponential vdW below is a CustomNonbondedForce and is correctly skipped.
    nb = openmm.NonbondedForce()
    nb.setNonbondedMethod(nb_method)
    try:
        nb.setCutoffDistance(cutoff)
    except Exception:
        pass                                          # NoCutoff has no cutoff distance
    # Reaction-field dielectric 1.0: at CutoffNonPeriodic this leaves plain shifted 1/r with no
    # geometry-dependent r^2 term, so barriers are invariant (GBSA later also sets 1.0).
    nb.setReactionFieldDielectric(1.0)
    for i in range(subset.n_particles):
        nb.addParticle(float(subset.charges[i]), 0.0, 0.0)
    if bond_pairs:
        nb.createExceptionsFromBonds(bond_pairs, g['coulomb14scale'], 1.0)
    nb.setForceGroup(NONBONDED_FORCE_GROUP)
    system.addForce(nb)

    # --- Double-exponential vdW (CustomNonbondedForce): pairs within 3 bonds excluded ---
    cnb = openmm.CustomNonbondedForce(_DEXP_ENERGY)
    cnb.setNonbondedMethod(cnb_method)
    try:
        cnb.setCutoffDistance(cutoff)
    except Exception:
        pass
    cnb.addGlobalParameter('alpha', g['dexp_alpha'])
    cnb.addGlobalParameter('beta', g['dexp_beta'])
    cnb.addPerParticleParameter('sigma')
    cnb.addPerParticleParameter('epsilon')
    for i in range(subset.n_particles):
        cnb.addParticle([float(subset.sigmas[i]), float(subset.epsilons[i])])
    if bond_pairs:
        cnb.createExclusionsFromBonds(bond_pairs, 3)
    cnb.setForceGroup(NONBONDED_FORCE_GROUP)
    system.addForce(cnb)

    # --- Scaled 1-4 and gated-1-3 double-exponential vdW (CustomBondForces) ---
    # The dexp CustomNonbondedForce above excludes everything within 3 bonds (1-2/1-3/1-4), so any
    # 1-4 (and gated coordination 1-3) pairs the checkpoint weights must be added back here at
    # reduced weight (a CustomNonbondedForce can express exclusions but not scaled exceptions).
    # "inert means ABSENT": each force is built only when its scale is non-zero, and garnet's
    # builders return None when the topology has no such pairs (e.g. no coordination centre). We
    # feed garnet's own make_*_vdw_force a fake `data` carrying the SUBSET bond graph in particle
    # indices + subset per-particle sigma/epsilon, so the pairs are computed on the same graph the
    # exclusions used (no double-counting) and the functional form stays single-sourced with garnet.
    if bond_pairs:
        vdw_data = SimpleNamespace(bond_is=[i for i, _ in bond_pairs],
                                   bond_js=[j for _, j in bond_pairs])
        vdw_params = {'sigma': subset.sigmas, 'epsilon': subset.epsilons,
                      'dexp_alpha': g['dexp_alpha'], 'dexp_beta': g['dexp_beta']}
        w14v = g.get('vdw14scale', 0.0)
        if w14v != 0.0:
            f14 = make_14_vdw_force(vdw_data, vdw_params, subset.n_particles, w14v, periodic=False)
            if f14 is not None:
                f14.setForceGroup(NONBONDED_FORCE_GROUP)
                system.addForce(f14)
        w13v = g.get('vdw13scale', 0.0)
        if w13v != 0.0:
            zl = [int(z) for z in all_atoms.elements.numbers]
            f13 = make_13_vdw_force(vdw_data, vdw_params, subset.n_particles, w13v,
                                    atomic_numbers=zl, periodic=False)
            if f13 is not None:
                f13.setForceGroup(NONBONDED_FORCE_GROUP)
                system.addForce(f13)

    _add_constraints(system, subset, all_atoms, sim_params, logger)
    return top, system


def _nb_methods(method):
    '''
    Translate the nonbonded method to ``(NonbondedForce enum, CustomNonbondedForce enum)``.

    ISOLDE stores the ``openmm.app`` method object (e.g. ``app.CutoffNonPeriodic``, what
    ``createSystem`` consumes); the low-level ``set*NonbondedMethod`` APIs want the force-class
    enums instead. Also accepts a low-level ``NonbondedForce`` enum (for tests). ``CustomNonbondedForce``
    supports only a hard cutoff, so periodic-Ewald/PME collapse to ``CutoffPeriodic`` for the vdW force.
    '''
    from openmm import app
    NB = openmm.NonbondedForce
    CNB = openmm.CustomNonbondedForce
    app_to_nb = {
        app.NoCutoff: NB.NoCutoff,
        app.CutoffNonPeriodic: NB.CutoffNonPeriodic,
        app.CutoffPeriodic: NB.CutoffPeriodic,
        app.Ewald: NB.Ewald,
        app.PME: NB.PME,
    }
    nb = app_to_nb.get(method)
    if nb is None:
        low_level = (NB.NoCutoff, NB.CutoffNonPeriodic, NB.CutoffPeriodic, NB.Ewald, NB.PME)
        nb = method if method in low_level else NB.CutoffNonPeriodic
    nb_to_cnb = {
        NB.NoCutoff: CNB.NoCutoff,
        NB.CutoffNonPeriodic: CNB.CutoffNonPeriodic,
        NB.CutoffPeriodic: CNB.CutoffPeriodic,
        NB.Ewald: CNB.CutoffPeriodic,
        NB.PME: CNB.CutoffPeriodic,
    }
    return nb, nb_to_cnb.get(nb, CNB.CutoffNonPeriodic)


def _add_constraints(system, subset, all_atoms, sim_params, logger):
    '''
    Add X-H (and, if requested, rigid-water) constraints honouring
    ``sim_params.rigid_bonds`` / ``rigid_water``, using garnet's ``r0`` as the
    constraint length.

    Every X-H bond is intra-residue, and ISOLDE's fixed shell is whole-residue,
    so no constraint spans the fixed/mobile boundary — the trap that makes
    OpenMM reject a constraint between a zero-mass and nonzero-mass particle
    (checked at Context creation, after ``set_fixed_atoms``) cannot arise here.
    '''
    from openmm.app import HBonds, AllBonds, HAngles
    rigid_bonds = getattr(sim_params, 'rigid_bonds', None)
    rigid_water = getattr(sim_params, 'rigid_water', False)
    if rigid_bonds is None and not rigid_water:
        return
    if rigid_bonds == HAngles:
        logger.warning('garnet: rigid_bonds=HAngles is not supported by the interim builder; '
                       'treating as HBonds.')

    elements = all_atoms.elements.numbers
    constrained = set()

    def _constrain(i, j, length_nm):
        key = (i, j) if i < j else (j, i)
        if key in constrained:
            return
        constrained.add(key)
        system.addConstraint(i, j, length_nm * unit.nanometer)

    constrain_all = (rigid_bonds == AllBonds)
    constrain_h = rigid_bonds in (HBonds, HAngles)
    bond_len = {}
    for i, j, k, r0 in subset.bonds:
        bond_len[(i, j) if i < j else (j, i)] = r0
        has_h = (elements[i] == 1 or elements[j] == 1)
        if constrain_all or (constrain_h and has_h):
            _constrain(i, j, r0)

    if not rigid_water:
        return

    # Rigid water: engage SETTLE by giving each 3-atom water (1 O + 2 H) a rigid triangle
    # (O-H, O-H, H-H). H-H is derived from garnet's O-H r0 and H-O-H theta0 via the law of cosines.
    index = {a: i for i, a in enumerate(all_atoms)}
    angle_theta = {(i, j, k): theta0 for i, j, k, _, theta0 in subset.angles}
    n_water = 0
    for res in all_atoms.unique_residues:
        ratoms = res.atoms
        znums = ratoms.elements.numbers
        if len(ratoms) != 3 or sorted(znums) != [1, 1, 8]:
            continue
        o = ratoms[list(znums).index(8)]
        hs = [a for a in ratoms if a.element.number == 1]
        if o not in index or any(h not in index for h in hs):
            continue
        oi = index[o]
        hi = [index[h] for h in hs]
        roh = [bond_len.get((oi, h) if oi < h else (h, oi)) for h in hi]
        theta = (angle_theta.get((hi[0], oi, hi[1]))
                 or angle_theta.get((hi[1], oi, hi[0])))
        if None in roh or theta is None:
            logger.warning('garnet: could not derive rigid-water geometry for residue '
                           '{} {}; leaving it flexible.'.format(res.name, res.number))
            continue
        for h, r in zip(hi, roh):
            _constrain(oi, h, r)
        hh = math.sqrt(max(0.0, roh[0] ** 2 + roh[1] ** 2 - 2 * roh[0] * roh[1] * math.cos(theta)))
        _constrain(hi[0], hi[1], hh)
        n_water += 1
    if n_water:
        logger.info('garnet: applied rigid-water (SETTLE) constraints to {} water(s).'.format(n_water))
