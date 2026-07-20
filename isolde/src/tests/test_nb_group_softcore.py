# Headless unit test for per-group soft-core nonbonded coupling (Phase 1).
#
# Runs on the OpenMM Reference platform; no GUI, no map, no ISOLDE simulation.
# Mirrors tests/test_symmetry_mdff_force.py: build tiny Systems, compare
# single-point energies, assert to ~machine precision. Validates:
#   (a) NBGroupNonbondedSoftcoreForce(n_nb_groups=1) == plain NonbondedSoftcoreForce
#       (byte-for-byte degradation) and the symmetry force rebased onto it is
#       unchanged in its symmetry-degenerate configuration;
#   (b) per-group pair_lambda energies for group pairs (0,0)/(1,0)/(1,1);
#   (c) THE LINCHPIN: set_coupling + updateParametersInContext changes the live
#       energy with no reinitialisation, and setParticleParameters changes group
#       membership live.

import numpy as np
import openmm as mm
from openmm import unit

from chimerax.isolde.openmm.custom_forces import (
    NonbondedSoftcoreForce, NBGroupNonbondedSoftcoreForce,
    SymmetryAwareNonbondedSoftcoreForce,
)

SIG, EPS = 0.30, 0.60          # nm, kJ/mol
R_REP = 0.24                   # nm; r/sigma = 0.8 -> repulsive regime
LAMBDA = 0.95                  # global ceiling
TOL = 1e-9

_failures = []
def _check(cond, msg):
    if cond:
        print("  ok:", msg)
    else:
        print("  FAIL:", msg)
        _failures.append(msg)

def _make_context(force, per_particle):
    '''per_particle: list of param-lists, one per particle.'''
    system = mm.System()
    for _ in per_particle:
        system.addParticle(12.0)
    for p in per_particle:
        force.addParticle(list(p))
    system.addForce(force)
    integ = mm.VerletIntegrator(1.0 * unit.femtosecond)
    context = mm.Context(system, integ, mm.Platform.getPlatformByName("Reference"))
    context.setPositions(np.array(
        [[0.0, 0.0, 0.0], [R_REP, 0.0, 0.0]]) * unit.nanometer)
    return context

def _energy(context):
    return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)

def run(session=None):
    # --- (a) parity: plain vs NBGroup(n_nb_groups=1) --------------------------
    plain = NonbondedSoftcoreForce(nb_lambda=LAMBDA)
    plain.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    e_plain = _energy(_make_context(plain, [[-0.5, SIG, EPS], [-0.5, SIG, EPS]]))

    g1 = NBGroupNonbondedSoftcoreForce(nb_lambda=LAMBDA, n_nb_groups=1)
    g1.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    e_g1 = _energy(_make_context(g1, [[-0.5, SIG, EPS], [-0.5, SIG, EPS]]))
    _check(abs(e_plain - e_g1) < TOL,
           f"NBGroup(n=1) parity with plain: {e_plain:.9f} vs {e_g1:.9f}")

    # --- symmetry force rebased onto NBGroup: degenerate config == plain ------
    sym = SymmetryAwareNonbondedSoftcoreForce(
        symmetry_ngroups=1, n_nb_groups=1, nb_lambda=LAMBDA)
    sym.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    # symmetry_ngroups=1 -> grouptable is 1x1 all-ones; symgroup=0 for both.
    e_sym = _energy(_make_context(sym, [[-0.5, SIG, EPS, 0.0], [-0.5, SIG, EPS, 0.0]]))
    _check(abs(e_plain - e_sym) < TOL,
           f"SymmetryAware(sym=1,nb=1) parity with plain: {e_plain:.9f} vs {e_sym:.9f}")

    # --- (b)+(c) grouped force, charge=0 -> pure LJ (monotonic in coupling) ---
    n = 3
    f = NBGroupNonbondedSoftcoreForce(nb_lambda=LAMBDA, n_nb_groups=n)
    f.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    ctx = _make_context(f, [[0.0, SIG, EPS, 0.0], [0.0, SIG, EPS, 1.0]])  # groups 0,1

    # reference: plain force, charge 0, same geometry -> full-coupling LJ
    ref = NonbondedSoftcoreForce(nb_lambda=LAMBDA)
    ref.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    e_ref = _energy(_make_context(ref, [[0.0, SIG, EPS], [0.0, SIG, EPS]]))

    e_full = _energy(ctx)     # table all 1 -> (0,1) coupling 1 -> == full
    _check(abs(e_full - e_ref) < TOL,
           f"grouped, all-coupled == plain full: {e_full:.9f} vs {e_ref:.9f}")
    _check(e_full > 0.0, f"repulsive regime gives positive clash energy ({e_full:.4f})")

    # THE LINCHPIN: soften (0,1) live via set_coupling + updateParametersInContext
    f.set_coupling(0, 1, 0.2, context=ctx)
    e_soft = _energy(ctx)
    _check(e_soft < e_full and e_soft > 0.0,
           f"live set_coupling(0,1,0.2) softens clash: {e_soft:.4f} < {e_full:.4f}")

    # restore live
    f.set_coupling(0, 1, 1.0, context=ctx)
    _check(abs(_energy(ctx) - e_full) < TOL, "live restore of coupling matches full")

    # membership: soften (0,1), then move particle 1 into group 0 -> (0,0)=1 -> full
    f.set_coupling(0, 1, 0.2, context=ctx)
    _check(_energy(ctx) < e_full, "re-softened before membership change")
    f.setParticleParameters(1, [0.0, SIG, EPS, 0.0])   # particle 1 -> group 0
    f.updateParametersInContext(ctx)
    _check(abs(_energy(ctx) - e_full) < TOL,
           "live membership change (group 1->0) restores full coupling")

    # internal coupling is independent: put both in group 1, soften only (1,1)
    f2 = NBGroupNonbondedSoftcoreForce(nb_lambda=LAMBDA, n_nb_groups=n)
    f2.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    ctx2 = _make_context(f2, [[0.0, SIG, EPS, 1.0], [0.0, SIG, EPS, 1.0]])  # both group 1
    _check(abs(_energy(ctx2) - e_full) < TOL, "group (1,1) at coupling 1 == full")
    f2.set_coupling(1, 1, 0.2, context=ctx2)
    _check(_energy(ctx2) < e_full, "softening (1,1) softens the internal pair")

    # rotafit copy-shadow scheme: with the target's copies in group 2, softening the
    # crystal self-contact (1,2) must NOT touch the target's INTERNAL (1,1) pair, but
    # must soften a (1,2) pair. This is what keeps the sidechain internally rigid while
    # its contact with its own symmetry image is softened.
    f_int = NBGroupNonbondedSoftcoreForce(nb_lambda=LAMBDA, n_nb_groups=3)
    f_int.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    ctx_int = _make_context(f_int, [[0.0, SIG, EPS, 1.0], [0.0, SIG, EPS, 1.0]])  # both grp 1
    e_int = _energy(ctx_int)
    f_int.set_coupling(1, 2, 0.2, context=ctx_int)
    _check(abs(_energy(ctx_int) - e_int) < TOL,
           "softening (1,2) leaves the (1,1) internal pair FULL")
    f_x = NBGroupNonbondedSoftcoreForce(nb_lambda=LAMBDA, n_nb_groups=3)
    f_x.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    ctx_x = _make_context(f_x, [[0.0, SIG, EPS, 1.0], [0.0, SIG, EPS, 2.0]])       # grp 1,2
    e_x = _energy(ctx_x)
    f_x.set_coupling(1, 2, 0.2, context=ctx_x)
    _check(_energy(ctx_x) < e_x, "softening (1,2) softens a (1,2) self-contact pair")

    # ================= GB implicit-solvent per-group (plan A7) ==============
    from openmm import app
    from chimerax.isolde.openmm.custom_forces import (
        SoftCoreGBSAGBnForce, NBGroupSoftCoreGBSAGBnForce)

    def _gb_context(gbforce, charges, groups=None):
        top = app.Topology()
        ch = top.addChain(); res = top.addResidue('UNK', ch)
        for i in range(len(charges)):
            top.addAtom(f'C{i}', app.Element.getBySymbol('C'), res)
        std = gbforce.getStandardParameters(top)          # (n, 2): (or, sr)
        nch = len(charges)
        pp = np.zeros((nch, 3)); pp[:, 0] = charges; pp[:, 1:] = std
        if groups is not None:
            gbforce._nb_groups = list(groups)
        gbforce.addParticles(pp)
        gbforce.finalize()
        system = mm.System()
        for _ in range(nch):
            system.addParticle(12.0)
        system.addForce(gbforce)
        ctx = mm.Context(system, mm.VerletIntegrator(1.0 * unit.femtosecond),
                         mm.Platform.getPlatformByName("Reference"))
        ctx.setPositions(np.array([[0.0, 0.0, 0.0], [R_REP, 0.0, 0.0]]) * unit.nanometer)
        return ctx

    q = [-0.5, -0.5]
    e_gb_plain = _energy(_gb_context(SoftCoreGBSAGBnForce(nb_lambda=LAMBDA), q))
    e_gb1 = _energy(_gb_context(
        NBGroupSoftCoreGBSAGBnForce(nb_lambda=LAMBDA, n_nb_groups=1), q))
    _check(abs(e_gb_plain - e_gb1) < 1e-6,
           f"GB NBGroup(n=1) parity with plain: {e_gb_plain:.6f} vs {e_gb1:.6f}")

    gbn = NBGroupSoftCoreGBSAGBnForce(nb_lambda=LAMBDA, n_nb_groups=3)
    ctxg = _gb_context(gbn, q, groups=[0, 1])
    e_gb_full = _energy(ctxg)
    _check(abs(e_gb_full - e_gb_plain) < 1e-6,
           f"GB grouped, all-coupled == plain: {e_gb_full:.6f} vs {e_gb_plain:.6f}")
    gbn.set_coupling(0, 1, 0.2, context=ctxg)          # linchpin on the GB force
    e_gb_soft = _energy(ctxg)
    _check(abs(e_gb_soft - e_gb_full) > 1e-6,
           f"GB live set_coupling(0,1,0.2) changes energy: {e_gb_soft:.6f} vs {e_gb_full:.6f}")
    gbn.set_coupling(0, 1, 1.0, context=ctxg)
    _check(abs(_energy(ctxg) - e_gb_full) < 1e-6, "GB live restore of coupling matches full")

    # ========== symmetry + nb-group COMPOSITION (Phase 5) ==========
    # Both mechanisms active on one force. Give every particle symgroup 0 ("real") so
    # the symmetry mask grouptable(0,0)=1 is identity -- this isolates the nb-group layer
    # composed *under* the symmetry wrapper, confirming the two masks coexist and the
    # per-group coupling still works. (Full symmetry-copy behaviour is covered by the
    # existing symmetry tests; the nb layer is orthogonal to it.)
    # (SymmetryAwareNonbondedSoftcoreForce is already imported at module scope.)
    from chimerax.isolde.openmm.custom_forces import SymmetrySoftCoreGBSAGBnForce

    # -- nonbonded --
    symnb = SymmetryAwareNonbondedSoftcoreForce(
        symmetry_ngroups=2, n_nb_groups=3, nb_lambda=LAMBDA)
    symnb.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    # per-particle: [charge, sigma, epsilon, nb_group, symgroup]
    ctx_sn = _make_context(symnb, [[0.0, SIG, EPS, 0.0, 0.0],
                                   [0.0, SIG, EPS, 1.0, 0.0]])
    e_sn = _energy(ctx_sn)
    _check(abs(e_sn - e_ref) < TOL,
           f"sym+nb nonbonded (symgroup 0, coupled) == plain full: {e_sn:.6f} vs {e_ref:.6f}")
    symnb.set_coupling(0, 1, 0.2, context=ctx_sn)
    _check(_energy(ctx_sn) < e_sn,
           "sym+nb nonbonded: softening nb coupling(0,1) softens the clash")

    # -- GB --
    def _symgb_context(charges, symgroups, nb_groups, ngroups=3):
        gb = SymmetrySoftCoreGBSAGBnForce(nb_lambda=LAMBDA, n_nb_groups=ngroups)
        gb._symmetry_ngroups = 2
        gb._symmetry_group_table = None          # fixed-weight symmetry table
        gb._symmetry_groups = list(symgroups)
        gb._nb_groups = list(nb_groups)
        top = app.Topology(); chn = top.addChain(); res = top.addResidue('UNK', chn)
        for i in range(len(charges)):
            top.addAtom(f'C{i}', app.Element.getBySymbol('C'), res)
        std = gb.getStandardParameters(top)
        pp = np.zeros((len(charges), 3)); pp[:, 0] = charges; pp[:, 1:] = std
        gb.addParticles(pp); gb.finalize()
        system = mm.System()
        for _ in charges:
            system.addParticle(12.0)
        system.addForce(gb)
        ctx = mm.Context(system, mm.VerletIntegrator(1.0 * unit.femtosecond),
                         mm.Platform.getPlatformByName("Reference"))
        ctx.setPositions(np.array([[0.0, 0.0, 0.0], [R_REP, 0.0, 0.0]]) * unit.nanometer)
        return gb, ctx

    gb_sn, ctx_gbsn = _symgb_context(q, [0, 0], [0, 1])
    e_gbsn = _energy(ctx_gbsn)
    _check(abs(e_gbsn - e_gb_plain) < 1e-6,
           f"sym+nb GB (symgroup 0, coupled) == plain GB: {e_gbsn:.6f} vs {e_gb_plain:.6f}")
    gb_sn.set_coupling(0, 1, 0.2, context=ctx_gbsn)
    _check(abs(_energy(ctx_gbsn) - e_gbsn) > 1e-6,
           "sym+nb GB: softening nb coupling(0,1) changes GB energy")
    # n_nb_groups=1: symmetry GB with the nb layer OFF must be byte-identical to the
    # plain symmetry path (here, all-real == plain GB).
    _gb1, ctx_gb1 = _symgb_context(q, [0, 0], [0, 0], ngroups=1)
    _check(abs(_energy(ctx_gb1) - e_gb_plain) < 1e-6,
           "sym GB n_nb_groups=1 (nb off) == plain GB")

    if _failures:
        print(f"\n{len(_failures)} CHECK(S) FAILED")
        raise SystemExit(1)
    print("\nALL PASS")

if 'session' in dict(globals()) and session is not None:
    run(session)
