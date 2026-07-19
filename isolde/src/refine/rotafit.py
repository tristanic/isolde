# @Author: Tristan Croll
# @Date:   05-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 05-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
isolde rotafit -- automated rotamer selection by *settling*.

Automates the tedious, error-prone middle step of the manual rotamer-fixing
workflow (start sim -> dial through previews -> commit the best): the human no
longer eyeballs rigid, pre-settle previews (which are misleading -- a fixed CA-CB
rotation can put even the correct rotamer well out of the map, and long lists like
Arg's are prevalence-sorted, not similarity-sorted). Instead, each library rotamer
is *settled* in the running (contracted) full-FF simulation and ranked by OpenMM's
own energy. The full force field reliably falls into the correct conformation once
placed in the right basin, so the settled energy -- not the pre-settle look -- is
the honest arbiter.

Requires a running ISOLDE simulation with the target residue mobile (start one with
``isolde sim start`` on the selection first, or this command will start a contracted
one). Drives ``sim_handler._simulation`` directly while paused, on the main thread.
The best rotamer is auto-committed but the simulation is LEFT RUNNING, so a poor
result is one "revert to checkpoint" away (and a live run then refines it further
as ISOLDE updates the map on R-factor improvement).
'''

import numpy as np

SEVERE_OVERLAP = 2.0     # Angstrom: a moved sidechain heavy atom closer than this
                         #   to a non-residue heavy atom => severe clash, cull it.
SETTLE_STEPS = 100       # 0 K dynamics steps per rotamer during the SOFT search
                         #   (soft-core vdW converges faster by stepping than by the
                         #   minimizer). Now that it's fast, settle each trial deeper.
POLISH_STEPS = 150       # extra 0 K steps applied at full soft-core stiffness to each
                         #   polished candidate, so the committed pose is tight and its
                         #   reported energy is honest (0 => skip the polish phase).
POLISH_TOP = 1           # how many of the top (soft-search-ranked) survivors to polish
                         #   at full stiffness and RE-RANK by the polished energy. 1 =
                         #   winner only (cheapest); N = top N; <=0 = all survivors.
                         #   The soft (0.6) ranking and the stiff (0.95) ranking can
                         #   disagree, so polishing more guards against a soft-search
                         #   tie picking the wrong basin -- at N x the polish cost.
RAMP_INCREMENTS = 5      # number of stages over which the polish ramps the soft-core
                         #   lambda from settle_lambda up to full stiffness, instead of
                         #   jumping straight there (which can jolt atoms sitting in a
                         #   soft overlap). polish_steps is split across the stages.
                         #   1 => the old direct jump.
ACCEPT_MARGIN = 1.0      # "first, do no harm" threshold, as a MULTIPLE of ISOLDE's
                         #   MDFF coupling constant (summed global_k) -- NOT an absolute
                         #   energy. The residue's CURRENT conformation competes as a
                         #   candidate; a library rotamer replaces it only if it wins by
                         #   at least (accept_margin x global_k) kJ/mol. Because global_k
                         #   is ISOLDE's sigma-normalised, resolution-calibrated map
                         #   weight (small for high-sigma cryo-EM maps, larger for
                         #   x-ray), the threshold auto-tracks the map's energy scale
                         #   instead of being fixed. 0 => any improvement wins; larger =>
                         #   stickier to the current pose. Falls back to an ABSOLUTE
                         #   kJ/mol value when the sim has no MDFF map (apo).
CURRENT_LABEL = '(current)'   # candidate label for the residue's starting conformation
SETTLE_LAMBDA = 0.6      # soft-core lambda during the settle. ISOLDE's live default
                         #   (~0.95) is stiff enough that a rotamer seeded into an
                         #   overlap can gain enough force to blow up the local system
                         #   before it relaxes. Softening the wall to ~0.6 lets clashy
                         #   seeds slide apart instead of exploding; lower than ~0.5
                         #   and real geometry starts to suffer. Restored after.
PIN_K = 5000.0           # kJ/mol/nm^2: spring constant for the position restraints
                         #   that pin the surrounding model during the settle (ISOLDE's
                         #   own default position-restraint strength). Because softening
                         #   the soft-core lambda softens EVERYTHING's nonbonded, pinning
                         #   the environment lets us pull settle_lambda down harder to
                         #   free a stuck target WITHOUT the shell deforming. Only the
                         #   target residue is left free. Restraints are restored to
                         #   their prior state afterward (no context reinit needed).
PIN_NEAR_CUTOFF = 5.0    # Angstrom: environment atoms within this distance of the
                         #   volume the residue + any rotamer can occupy are pinned at
                         #   PIN_K (they may legitimately need to accommodate the target).
PIN_DISTANT_MULT = 10.0  # atoms FURTHER than PIN_NEAR_CUTOFF are pinned at PIN_K x this
                         #   (an order of magnitude stiffer): they can't accommodate the
                         #   target, so any motion they make is pure ranking NOISE in the
                         #   core nonbonded term. Freeze them hard.

# Per-group soft-core group ids used by rotafit (see _Softener). The target residue
# goes in group 1; everything else stays in group 0. (The coupling table is provisioned
# with more slots by default -- SimParams.nb_groups_max -- but rotafit only needs these
# two, and settles one target at a time.)
ENV_NB_GROUP = 0
TARGET_NB_GROUP = 1


class _Softener:
    '''The "soften the target against its surroundings" knob, with two backends.

    PREFERRED (``use_groups=True``): per-group soft-core coupling. The target
    residue is placed in its own nonbonded group and only its coupling to the rest
    of the model (group 0) is softened; the environment keeps full-strength
    interactions, so it holds its own shape WITHOUT position restraints. ``get`` /
    ``set`` read/write the group-pair coupling ``(TARGET, ENV)``.

    LEGACY (``use_groups=False``): the global ``softcore_lambda`` scalar softens
    every pair, so the environment must be pinned separately (see
    :func:`_pin_environment`). ``get`` / ``set`` read/write ``softcore_lambda``.

    Everything else in rotafit talks only to this object, so the settle/polish/score
    logic is identical in both modes -- only the knob differs.'''
    def __init__(self, sh, use_groups):
        self.sh = sh
        self.use_groups = use_groups

    def get(self):
        if self.use_groups:
            return self.sh.get_nb_coupling(TARGET_NB_GROUP, ENV_NB_GROUP)
        return self.sh.softcore_lambda

    def set(self, value):
        if value is None:
            return
        if self.use_groups:
            self.sh.set_nb_coupling(TARGET_NB_GROUP, ENV_NB_GROUP, value)
        else:
            self.sh.softcore_lambda = value

    def assign_target(self, residue):
        '''Put the target residue in its own group (no-op in legacy mode).'''
        if self.use_groups:
            self.sh.assign_nb_group(residue.atoms, TARGET_NB_GROUP)

    def release_target(self, residue):
        '''Return the target residue to the environment group (no-op in legacy mode).'''
        if self.use_groups:
            self.sh.assign_nb_group(residue.atoms, ENV_NB_GROUP)


def _rotamer_poses(session, residue, base_res_coords):
    '''Enumerate the library rotamers as full residue-atom coordinate arrays
    (heavy + H, in residue.atoms order), by applying each target's chi angles to
    ``base_res_coords`` (the residue's CURRENT simulation coordinates) -- the same
    conformations the preview buttons show. Returns [(name, coords), ...].'''
    from chimerax.isolde import session_extensions as sx
    from chimerax.geometry import rotation, matrix, dihedral
    from math import degrees
    rota_mgr = sx.get_rotamer_mgr(session)
    rot = rota_mgr.get_rotamer(residue)
    if rot is None:
        return []
    idx_of = {a: i for i, a in enumerate(residue.atoms)}
    nchi = rot.num_chi_dihedrals
    poses = []
    for it in range(rot.num_targets):
        t = rot.get_target(it)
        coords = base_res_coords.copy()
        for i in range(nchi):
            chi = rot.chi_dihedrals[i]
            p = np.array([coords[idx_of[a]] if a in idx_of else np.array(a.coord)
                          for a in chi.atoms])
            axis = p[2] - p[1]
            center = matrix.project_to_axis(p[3], axis, p[1])
            cur = dihedral(p[0], p[1], p[2], p[3])
            tf = rotation(axis, degrees(t['Angles'][i]) - cur, center)
            rows = [idx_of[a] for a in rot.moving_atoms(i) if a in idx_of]
            if rows:
                coords[rows] = tf.transform_points(coords[rows])
        poses.append((t.get('Name', '?'), coords))
    return poses


def _push_now(sh, coords):
    '''Push full-construct coords (Angstrom) into the running sim IMMEDIATELY while
    paused, so they survive resume and are seen by the next step(). Uses the public
    push_coords_to_sim(immediate=True) rather than reaching for a private method.'''
    sh.push_coords_to_sim(coords, immediate=True)


def _sim_coords(sh):
    '''Current simulation coordinates (Angstrom, construct/_atoms order).

    Returns ONLY the real-atom rows. In a symmetry-aware simulation the OpenMM
    System carries extra symmetry-copy virtual-site particles after the real atoms
    ([0, n_real) real, [n_real, n_total) copies), so the raw context state is longer
    than sh._atoms; everything in rotafit is keyed to sh._atoms, so slice to it.
    (A no-op when symmetry is off: n_total == n_real.)'''
    from openmm import unit
    st = sh._simulation.context.getState(getPositions=True)
    coords = st.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    return coords[:len(sh._atoms)]


def _score_groups():
    '''Force groups to rank rotamers on: the CORE force field (bonded + nonbonded) plus
    the MDFF map term -- i.e. real physics + density fit. EXCLUDES ISOLDE's
    RESTRAINT_FORCE_GROUP: every restraint (our environment position pins, plus omega/
    phi-psi/distance/chiral restraints) lumps into that one group, and its energy is
    large, pose-independent NOISE for this decision -- the 5000 kJ/mol/nm^2 pin on
    hundreds of atoms especially. Ranking on the total energy let that noise dominate.'''
    from chimerax.isolde.openmm.openmm_interface import (CORE_FORCE_GROUPS,
                                                         MAP_FORCE_GROUP)
    return set(CORE_FORCE_GROUPS) | {MAP_FORCE_GROUP}


def _score_energy(ctx, softener=None, score_lambda=0.0):
    '''Potential energy over the scoring force groups only (core FF + MDFF), in kJ/mol.
    If ``score_lambda`` > 0 (and a ``softener`` is given), evaluate at that softening
    level -- a softer wall means a smaller, less jitter-sensitive nonbonded term --
    regardless of the level the pose was settled at, then restore. In per-group mode
    this softens only the target's coupling to its surroundings; in legacy mode it is
    the global soft-core lambda. A single-point re-evaluation, no stepping.'''
    from openmm import unit
    changed = bool(score_lambda) and softener is not None
    if changed:
        prev = softener.get()
        softener.set(score_lambda)
    e = ctx.getState(getEnergy=True, groups=_score_groups()).getPotentialEnergy(
        ).value_in_unit(unit.kilojoule_per_mole)
    if changed:
        softener.set(prev)
    return e


def _relax(sh, integrator, ctx, steps, minimize):
    '''Advance the sim toward a local minimum from the current coords. ALWAYS steps the
    main integrator at 0 K first: Newtonian (damped) dynamics carries momentum over small
    barriers and SEATS the pose in the local density well. (A pure minimise from a
    freshly-placed rotamer instead quenches into the nearest minimum, which under the
    gentle cryo-EM map is often off-density -- it has no momentum to ride into the well.)
    Then, if ``minimize``, energy-minimise to converge cleanly to that seated minimum --
    deterministic, killing the run-to-run jitter in the ranking energy. So minimize=True
    is a HYBRID: seat with dynamics, converge with minimisation.'''
    from openmm import unit
    ctx.setVelocitiesToTemperature(0 * unit.kelvin)   # clean 0 K descent
    integrator.step(steps)                            # dynamics: fall into the map well
    if minimize:
        sh._simulation.minimizeEnergy()               # converge to the seated minimum


def _energy_breakdown(ctx):
    '''Per-force-group potential energies (kJ/mol) for the scoring groups, for debug --
    lets us see WHICH term (nonbonded, map, ...) separates a suspicious winner from the
    current conformation. A handful of extra getState calls, but they are cheap.'''
    from openmm import unit
    from chimerax.isolde.openmm.openmm_interface import (
        DEFAULT_FORCE_GROUP, BOND_FORCE_GROUP, ANGLE_FORCE_GROUP,
        DIHEDRAL_FORCE_GROUP, NONBONDED_FORCE_GROUP, MAP_FORCE_GROUP)
    groups = [('def', DEFAULT_FORCE_GROUP), ('bond', BOND_FORCE_GROUP),
              ('ang', ANGLE_FORCE_GROUP), ('dih', DIHEDRAL_FORCE_GROUP),
              ('nonbond', NONBONDED_FORCE_GROUP), ('map', MAP_FORCE_GROUP)]
    return [(name, ctx.getState(getEnergy=True, groups={g}).getPotentialEnergy()
             .value_in_unit(unit.kilojoule_per_mole)) for name, g in groups]


def _breakdown_str(ctx):
    return ', '.join('%s %.0f' % (n, e) for n, e in _energy_breakdown(ctx))


def _reachable_cloud(session, sh, residue):
    '''Heavy-atom coordinate cloud spanning the residue's CURRENT position AND every
    library-rotamer position -- the volume any candidate could occupy. Used to split the
    environment into "near" atoms (that may need to accommodate the target -> softer pin)
    and "distant" atoms (frozen hard). Returns an (N,3) array, or None if no rotamers.'''
    base = _sim_coords(sh)
    ridx = sh._atoms.indices(residue.atoms)
    if np.any(ridx < 0):
        return None
    heavy = residue.atoms.element_names != 'H'
    clouds = [base[ridx][heavy]]
    for _name, pcoords in _rotamer_poses(session, residue, base[ridx]):
        clouds.append(pcoords[heavy])
    return np.vstack(clouds)


def _sim_coupling_constant(isolde):
    '''Summed MDFF coupling constant (global_k) over the maps driving the running sim.
    This is ISOLDE's OWN map weight: sigma-normalised (the map force divides energy by
    the map sigma) and calibrated per resolution, so it is small for high-sigma cryo-EM
    maps and larger for x-ray. Scaling rotafit's do-no-harm margin by it makes the
    threshold track the map's energy scale automatically -- exactly the coupling ISOLDE
    already uses to drive the sim, so we don't reinvent a map weight. Returns 0.0 if the
    sim has no MDFF maps (apo).'''
    sm = getattr(isolde, 'sim_manager', None)
    mgrs = getattr(sm, 'mdff_mgrs', None) if sm is not None else None
    if not mgrs:
        return 0.0
    total = 0.0
    for mgr in mgrs.values():
        try:
            total += float(mgr.global_k)
        except Exception:
            pass
    return total


def _pin_environment(prm, mobile_atoms, residue, spring_constant, near_ref_coords,
                     near_cutoff, distant_multiplier):
    '''Enable ISOLDE's always-present (default-disabled) position restraints on every
    MOBILE HEAVY atom EXCEPT the target residue, anchored to their current positions, so
    an aggressive soft-search lambda can't push the surrounding model around -- only the
    target residue is left free to settle.

    TWO-TIER stiffness: atoms within ``near_cutoff`` of ``near_ref_coords`` (the volume
    the residue + any rotamer can reach) are pinned at ``spring_constant`` -- they may
    legitimately need to give a little. Atoms further out are pinned at
    ``spring_constant * distant_multiplier`` (an order of magnitude stiffer): they can't
    accommodate the target, so their motion is only noise in the ranking energy.

    Returns (prs, saved, n_near, n_far) where saved is the prior state for restoration,
    or None if there is nothing to pin. Applied immediately while the sim is paused.'''
    in_res = mobile_atoms.mask(residue.atoms)
    heavy = mobile_atoms.element_names != 'H'
    env = mobile_atoms[heavy & ~in_res]
    if not len(env):
        return None
    prs = prm.get_restraints(env)
    if not len(prs):
        return None
    saved = (prs.enableds.copy(), prs.targets.copy(), prs.spring_constants.copy())
    coords = prs.atoms.coords
    ks = np.full(len(prs), spring_constant * distant_multiplier)
    if near_ref_coords is not None and len(near_ref_coords):
        from chimerax.geometry import find_close_points
        near_idx, _ = find_close_points(coords, near_ref_coords, near_cutoff)
        ks[near_idx] = spring_constant
    n_near = int(np.sum(ks == spring_constant))
    prs.targets = coords                    # pin where they are right now
    prs.spring_constants = ks
    prs.enableds = True
    return prs, saved, n_near, len(prs) - n_near


def _restore_restraints(prs, saved):
    '''Restore position restraints to their pre-pin state (targets/springs first, then
    the enabled flags -- which for the environment is back to disabled by default).'''
    enableds, targets, spring_constants = saved
    prs.targets = targets
    prs.spring_constants = spring_constants
    prs.enableds = enableds


def _commit_best(sh, softener, residue, results, base, polish_steps, polish_top,
                 ramp_increments, saved_lambda, settle_lambda, eff_margin, outcomes,
                 dlog, debug, minimize=False, score_lambda=0.0):
    '''Polish the top candidate(s) at full stiffness, re-rank, apply the do-no-harm
    margin against the current conformation, and commit the winner. Appends the
    per-residue outcome (change/keep) to ``outcomes``.'''
    # Polish the top candidates at full soft-core stiffness and RE-RANK by the polished
    # energy: the soft (search) lambda and the stiff (live) lambda can disagree, so
    # re-ranking guards against a soft-search tie committing the wrong basin. polish_top:
    # N = top N; <=0 = all survivors; 1 = winner only. ALWAYS polish the current
    # conformation too (do-no-harm needs it judged at the same lambda as its rivals).
    candidates = results
    if polish_steps and saved_lambda is not None:
        n_polish = (len(results) if polish_top <= 0
                    else min(polish_top, len(results)))
        polish_set = list(results[:n_polish])
        if not any(cur for _e, _n, _c, cur in polish_set):
            curr = next((t for t in results if t[3]), None)
            if curr is not None:
                polish_set.append(curr)
        ctx = sh._simulation.context
        integrator = sh._main_integrator
        ridx = sh._atoms.indices(residue.atoms)
        polished = []
        for _e, pname, pcoords, pcur in polish_set:
            pe, pc = _ramped_polish(sh, softener, integrator, ctx, ridx, pcoords,
                                    pcoords[ridx], polish_steps, settle_lambda,
                                    saved_lambda, ramp_increments, minimize, score_lambda)
            polished.append((pe, pname, pc, pcur))
            if debug:                          # per-group breakdown to spot the culprit
                dlog('    %-9s [%s]' % (pname, _breakdown_str(ctx)))
        softener.set(settle_lambda)            # back to soft for the next residue
        polished.sort(key=lambda t: t[0])
        candidates = polished
        plist = ', '.join('%s %.0f' % (n, e) for e, n, _c, _cur in polished)
        dlog('  polished %d (lambda %.2f->%.2f x%d): %s'
             % (len(polished), settle_lambda, saved_lambda, ramp_increments, plist))
    # "First, do no harm": keep the CURRENT conformation unless a rival beats it by
    # >= eff_margin (else the winner is just settling noise).
    best_e, best_name, best_coords = candidates[0][:3]
    kept_current = candidates[0][3]
    curr = next((t for t in candidates if t[3]), None)
    if curr is not None and not kept_current and (curr[0] - best_e) < eff_margin:
        best_e, best_name, best_coords = curr[:3]
        kept_current = True
        dlog('  do-no-harm: winner beats current by <%.2f kJ/mol -> keeping current'
             % eff_margin)
    dlog('  -> commit "%s" (E=%.0f kJ/mol)' % (best_name, best_e))
    # Commit ONLY the target residue, grafted onto the ORIGINAL environment (base). The
    # settled/polished poses are whole-construct snapshots whose ENVIRONMENT was perturbed
    # (pinned + settled/minimised) during evaluation; committing that wholesale moves the
    # surroundings too -- invisible for a short 0 K step, but a minimise shifts them
    # enough to look like a bad commit. rotafit only places the target rotamer, so move
    # its atoms and leave everything else exactly as it was. Keeping current => no change
    # at all (commit == base); the live sim resolves any interface afterward.
    ridx = sh._atoms.indices(residue.atoms)
    commit_coords = base.copy()
    if not kept_current:
        commit_coords[ridx] = best_coords[ridx]
    _push_now(sh, commit_coords)               # into the sim (survives resume)
    sh._atoms.coords = commit_coords           # and into the ChimeraX model
    outcomes.append((residue, 'keep' if kept_current else 'change', best_name))


def _severe_clash(pose_coords, residue, moved_names, env_coords):
    '''True if any MOVED heavy sidechain atom in this rotamer pose comes within
    SEVERE_OVERLAP of a fixed environment heavy atom -- an obvious non-starter.'''
    ratoms = residue.atoms
    is_moved = np.isin(ratoms.names, list(moved_names))
    is_heavy = ratoms.element_names != 'H'
    sel = is_moved & is_heavy
    if not np.any(sel) or not len(env_coords):
        return False
    from chimerax.geometry import find_close_points
    close, _ = find_close_points(pose_coords[sel], env_coords, SEVERE_OVERLAP)
    return len(close) > 0


def _start_sim_on(session, isolde, residues):
    '''Start a simulation around the target residue(s) using ISOLDE's STANDARD sim-start
    selection (its normal padding + soft-shell). The generous buffer gives the local
    environment room to relax around a re-fitted rotamer -- important for fixes that need
    the neighbours to accommodate (an over-contracted region can leave no such room).'''
    from chimerax.atomic import Residues
    from chimerax.core.commands import run
    session.selection.clear()
    Residues(residues).atoms.selected = True
    run(session, 'isolde sim start sel')


def _settle_pose(sh, integrator, ctx, ridx, base, rcoords, steps, minimize=False):
    '''One settle: overwrite the residue rows (ridx) of the construct with rcoords, push,
    and relax (minimise or 0 K step) AT THE CURRENT soft-core lambda. Returns
    (score_energy_kJmol, settled_construct_coords). Used for the per-trial soft search.'''
    coords = base.copy()
    coords[ridx] = rcoords
    _push_now(sh, coords)                       # immediate push (survives resume)
    _relax(sh, integrator, ctx, steps, minimize)
    return _score_energy(ctx), _sim_coords(sh)


def _ramped_polish(sh, softener, integrator, ctx, ridx, base, rcoords, total_steps,
                   lambda_from, lambda_to, increments, minimize=False, score_lambda=0.0):
    '''Polish a candidate by RAMPING the softening level from ``lambda_from`` up to
    ``lambda_to`` over ``increments`` stages instead of jumping straight to full
    stiffness. A sudden jump can jolt atoms sitting in a soft overlap that abruptly
    becomes a hard wall; ramping lets the system stiffen adiabatically and follow the
    potential down. ``total_steps`` is split across the stages, relaxing (minimise or
    0 K step) at each. The softening is applied through ``softener`` (per-group coupling
    of the target to its surroundings, or legacy global lambda). The final energy is
    read at ``score_lambda`` if set (else at ``lambda_to``). Returns
    (score_energy_kJmol, settled_construct_coords). ``increments==1`` reduces to a
    single settle at ``lambda_to`` (the old behaviour).'''
    coords = base.copy()
    coords[ridx] = rcoords
    _push_now(sh, coords)                       # push once; stages continue from here
    increments = max(1, int(increments))
    steps_per = max(1, int(round(total_steps / increments)))
    lambdas = np.linspace(lambda_from, lambda_to, increments + 1)[1:]  # end at lambda_to
    for lam in lambdas:
        softener.set(float(lam))
        _relax(sh, integrator, ctx, steps_per, minimize)
    return _score_energy(ctx, softener, score_lambda), _sim_coords(sh)


def _settle_and_rank(isolde, residue, temperature, settle_steps, log, minimize=False):
    '''Settle each (non-culled) library rotamer in the running sim and return the
    ranked results [(energy_kJmol, name, settled_construct_coords), ...] best-first,
    plus the original construct coords. Energies are at the CURRENT (soft, search)
    lambda -- comparable to each other, the honest arbiter between rotamer basins; the
    winner is polished separately at full stiffness. Sim paused, main thread.'''
    session = residue.structure.session
    sh = isolde.sim_handler
    ctx = sh._simulation.context
    integrator = sh._main_integrator           # step THIS, not sh._simulation.step:
                                               # the latter routes through ISOLDE's
                                               # fast-atom surveillance, which can veto
                                               # / interrupt a clashy settle before we
                                               # get to score+discard it. The main
                                               # integrator does pure Newtonian steps.
    atoms = sh._atoms                          # construct-ordered atom array
    ridx = atoms.indices(residue.atoms)        # this residue's rows in the construct
    if np.any(ridx < 0):
        raise ValueError('target residue is not part of the running simulation')
    base = _sim_coords(sh)                     # the SIM's current coords (Angstrom)

    # Fixed environment heavy atoms (everything in the construct except this residue)
    heavy = atoms.element_names != 'H'
    env_mask = heavy.copy()
    env_mask[ridx] = False
    env_coords = base[env_mask]

    poses = _rotamer_poses(session, residue, base[ridx])
    if not poses:
        return None, base

    # rotamer moving-atom names (for the clash pre-filter)
    from chimerax.isolde import session_extensions as sx
    rot = sx.get_rotamer_mgr(session).get_rotamer(residue)
    moved_names = set()
    for i in range(rot.num_chi_dihedrals):
        moved_names.update(rot.moving_atoms(i).names)

    # "First, do no harm": the CURRENT conformation is always a candidate (exempt from
    # the clash cull) so an already-correct fit can win and be kept.
    poses = [(CURRENT_LABEL, base[ridx].copy())] + poses

    sh.temperature = temperature               # canonical setter (routes to the
                                               # main integrator; the CompoundIntegrator
                                               # has no setTemperature of its own)
    results = []
    n_cull = 0
    for name, rcoords in poses:
        is_current = (name == CURRENT_LABEL)
        if not is_current and _severe_clash(rcoords, residue, moved_names, env_coords):
            n_cull += 1
            continue
        e, settled = _settle_pose(sh, integrator, ctx, ridx, base, rcoords,
                                  settle_steps, minimize)
        results.append((e, name, settled, is_current))
    results.sort(key=lambda t: t[0])
    log('isolde rotafit %s: %d rotamers (+current), %d culled (severe clash), %d settled'
        % (residue, len(poses) - 1, n_cull, len(results)))
    return results, base


def rotafit(session, residues=None, temperature=0.0, settle_steps=SETTLE_STEPS,
            polish_steps=POLISH_STEPS, polish_top=POLISH_TOP,
            ramp_increments=RAMP_INCREMENTS, accept_margin=ACCEPT_MARGIN,
            settle_lambda=SETTLE_LAMBDA, pin_environment=True, pin_k=PIN_K,
            pin_near_cutoff=PIN_NEAR_CUTOFF, pin_distant_multiplier=PIN_DISTANT_MULT,
            minimize=True, score_lambda=0.0, allow_multiple=False,
            nb_groups=True, apply=True, debug=False):
    '''Automated rotamer fit: settle every viable library rotamer of the selected
    residue(s) in the simulation and commit the lowest-energy one. If no simulation is
    running, one is started around the target(s) using ISOLDE's standard sim-start
    selection (its normal padding/soft-shell -- the buffer gives the local environment
    room to relax, which some fixes need). Logs a single summary line unless ``debug
    True`` (then the full per-rotamer ranking, polish and do-no-harm decisions are
    logged). By default only a SINGLE residue may be targeted (see ``allowMultiple``).'''
    from chimerax.core.errors import UserError
    from chimerax.atomic import Residue, selected_residues
    isolde = getattr(session, 'isolde', None)
    if isolde is None or not session.ui.is_gui:
        raise UserError('isolde rotafit requires the ISOLDE GUI. Run "isolde start" '
                        'first.')
    if residues is None:
        residues = selected_residues(session)
    if residues is None or len(residues) == 0:
        raise UserError('isolde rotafit: no residues specified or selected.')
    dlog = session.logger.info if debug else (lambda *a, **k: None)
    targets = [r for r in residues if r.polymer_type == Residue.PT_AMINO]
    if not targets:
        raise UserError('isolde rotafit: no protein residues in the selection.')
    if len(targets) > 1 and not allow_multiple:
        raise UserError(
            'isolde rotafit: %d residues selected. Each residue is settled '
            'individually and can take a few seconds, so this is disabled by default to '
            'avoid an accidental whole-model run. To fit them all, repeat the command '
            'with "allowMultiple true".' % len(targets))

    if not isolde.simulation_running:
        dlog('isolde rotafit: no simulation running - starting one around the '
             'target(s)...')
        _start_sim_on(session, isolde, targets)

    sh = isolde.sim_handler
    # Use per-group coupling if requested AND available on this sim. It is available by
    # default (SimParams.nb_groups_max defaults to 4, so every sim is built group-aware);
    # it reports unavailable only if the user disabled it (nb_groups_max=1) or the sim
    # can't fully support it (e.g. symmetry + GBSA -- see SimHandler.nb_groups_enabled),
    # in which case we fall back to the legacy global-lambda + environment-pinning path.
    use_groups = bool(nb_groups) and getattr(sh, 'nb_groups_enabled', False)
    if nb_groups and not use_groups:
        dlog('isolde rotafit: per-group coupling unavailable for this simulation; '
             'falling back to global-lambda + environment pinning.')
    params = dict(temperature=temperature, settle_steps=settle_steps,
                  polish_steps=polish_steps, polish_top=polish_top,
                  ramp_increments=ramp_increments, accept_margin=accept_margin,
                  settle_lambda=settle_lambda, pin_environment=pin_environment,
                  pin_k=pin_k, pin_near_cutoff=pin_near_cutoff,
                  pin_distant_multiplier=pin_distant_multiplier,
                  minimize=minimize, score_lambda=score_lambda,
                  use_groups=use_groups, apply=apply, debug=debug)
    if sh.pause:
        # Already paused => the sim thread is idle => safe to drive it now.
        _run_rotafit(session, isolde, targets, resume_after=False, **params)
    else:
        # `sh.pause = True` only sets a FLAG -- the sim thread keeps running until
        # it next returns. Driving _simulation now would race the worker thread.
        #
        # Defer to the 'sim paused' trigger, NOT 'coord update'. 'coord update' fires
        # *inside* _update_coordinates_and_repeat (openmm_interface.py:1671), BEFORE
        # its pause/repeat branch (:1689). If we resume (sh.pause=False) from there,
        # the pause setter schedules a repeat via _resume(), and then the still-live
        # _update_coordinates_and_repeat frame ALSO schedules one -- two overlapping
        # step loops from one cycle. Each leaked loop independently fires
        # 'sim terminated' on shutdown, hence the string of remove_change_cb
        # tracebacks. 'sim paused' fires from a 'frame drawn' handler (:1695) AFTER
        # that function has returned down the pause branch (which did NOT repeat), so
        # the sim is genuinely quiescent and resuming schedules exactly one loop.
        sh.pause = True
        def _deferred(trigger_name, data, session=session, isolde=isolde,
                      targets=targets, params=params):
            from chimerax.core.triggerset import DEREGISTER
            if isolde.simulation_running:
                try:
                    _run_rotafit(session, isolde, targets, resume_after=True, **params)
                except Exception:
                    session.logger.report_exception()
            return DEREGISTER                  # one-shot
        sh.triggers.add_handler('sim paused', _deferred)
        dlog('isolde rotafit: scheduled - will run at the next safe pause point.')


def _run_rotafit(session, isolde, targets, temperature, settle_steps, polish_steps,
                 polish_top, ramp_increments, accept_margin, settle_lambda,
                 pin_environment, pin_k, pin_near_cutoff, pin_distant_multiplier,
                 minimize, score_lambda, apply, debug, use_groups, resume_after):
    '''Do the actual rotamer settling. MUST be called only when the sim thread is
    idle (paused and a 'sim paused' has fired), so main-thread access to
    _simulation is safe. Verbose per-rotamer logging is gated behind ``debug``; the
    standard path logs a single summary line at the end.'''
    sh = isolde.sim_handler
    log = session.logger.info
    dlog = log if debug else (lambda *a, **k: None)
    isolde.checkpoint()                        # so the user can revert if unhappy
    # Restore the sim's configured temperature after the 0 K settle. Read it from
    # sim_params (the configured target to return to; unchanged by the settle).
    st = isolde.sim_params.temperature
    from openmm import unit
    saved_temp = (st.value_in_unit(unit.kelvin) if hasattr(st, 'value_in_unit')
                  else float(st))
    softener = _Softener(sh, use_groups)
    # Per-group mode: the target's coupling to its surroundings (full = 1.0), with the
    # environment held by its own full-strength force field -- NO pinning. Legacy mode:
    # the global soft-core lambda (None if the context has no such parameter).
    saved_lambda = softener.get()
    if use_groups:
        dlog('isolde rotafit: per-group soft-core coupling (environment unpinned)')
    # Scale the do-no-harm margin by ISOLDE's MDFF coupling constant so it tracks the
    # map's energy scale (x-ray vs high-sigma cryo-EM) instead of being a fixed kJ/mol.
    gk = _sim_coupling_constant(isolde)
    eff_margin = accept_margin * gk if gk > 0 else accept_margin
    dlog('isolde rotafit: MDFF coupling global_k=%.3g -> do-no-harm margin %.2f kJ/mol'
         % (gk, eff_margin))
    # Environment pinning: freeze the surrounding model with position restraints so a
    # low settle_lambda can free the target without deforming the shell.
    # Environment pinning is the LEGACY mechanism (needed because a global soft-core
    # lambda softens everything); per-group mode keeps the environment stiff via its
    # own force field, so pinning is neither used nor set up.
    prm = mobile_atoms = None
    if pin_environment and not use_groups:
        try:
            from chimerax.isolde import session_extensions as sx
            prm = sx.get_position_restraint_mgr(targets[0].structure, True)
            mobile_atoms = isolde.sim_manager.sim_construct.mobile_atoms
        except Exception as e:
            dlog('  (environment pinning unavailable: %s)' % e)
            prm = None
    outcomes = []                              # (residue, kind, name) per target
    try:
        if saved_lambda is not None:
            softener.set(settle_lambda)         # soften (target-vs-environment coupling
                                                # in per-group mode; the global lambda in
                                                # legacy) so a clashy seed relaxes
        for r in targets:
            pin = None
            softener.assign_target(r)           # target -> its own nb group (no-op legacy)
            if prm is not None:
                try:
                    cloud = _reachable_cloud(session, sh, r)
                    pin = _pin_environment(prm, mobile_atoms, r, pin_k, cloud,
                                           pin_near_cutoff, pin_distant_multiplier)
                    if pin is not None:
                        dlog('  pinned environment: %d near @ k=%g, %d distant @ k=%g '
                             'kJ/mol/nm^2' % (pin[2], pin_k, pin[3],
                                              pin_k * pin_distant_multiplier))
                except Exception as e:
                    dlog('  (could not pin environment: %s)' % e)
                    pin = None
            try:
                try:
                    results, base = _settle_and_rank(isolde, r, temperature,
                                                     settle_steps, dlog, minimize)
                except Exception as e:
                    session.logger.warning('isolde rotafit: skipped %s (%s)' % (r, e))
                    outcomes.append((r, 'skip', None))
                    continue
                if not results:
                    _push_now(sh, base)        # nothing viable -> restore
                    outcomes.append((r, 'none', None))
                    continue
                shortlist = ', '.join('%s %.0f' % (n, e)
                                      for e, n, _c, _cur in results[:6])
                dlog('  search ranking @ lambda=%.2f: %s' % (settle_lambda, shortlist))
                if not apply:
                    _push_now(sh, base)         # nothing committed -> restore
                    outcomes.append((r, 'noapply', None))
                    continue
                _commit_best(sh, softener, r, results, base, polish_steps, polish_top,
                             ramp_increments, saved_lambda, settle_lambda, eff_margin,
                             outcomes, dlog, debug, minimize, score_lambda)
            finally:
                softener.release_target(r)     # target back to env group (no-op legacy)
                if pin is not None:
                    _restore_restraints(pin[0], pin[1])
    finally:
        sh.temperature = saved_temp            # don't leave the sim at 0 K
        if saved_lambda is not None:
            softener.set(saved_lambda)         # restore full coupling / soft-core lambda
        if resume_after:
            sh.pause = False                   # resume (we requested the pause)
    log(_summary_line(outcomes))


def _summary_line(outcomes):
    '''One-line summary of a rotafit run (the only log line in the non-debug path).'''
    committed = any(k in ('change', 'keep') for _r, k, _n in outcomes)
    tail = '; sim left running ("isolde sim revert" to undo)' if committed else ''
    if len(outcomes) == 1:
        r, kind, name = outcomes[0]
        msg = {
            'change': lambda: 'fitted rotamer "%s"' % name,
            'keep':   lambda: 'current conformation kept (already best fit)',
            'none':   lambda: 'no rotamers to try',
            'skip':   lambda: 'skipped (see warning)',
            'noapply': lambda: 'apply false - nothing committed',
        }[kind]()
        return 'isolde rotafit %s: %s%s' % (r, msg, tail)
    n_change = sum(1 for _r, k, _n in outcomes if k == 'change')
    n_keep = sum(1 for _r, k, _n in outcomes if k == 'keep')
    n_other = len(outcomes) - n_change - n_keep
    return ('isolde rotafit: %d fitted, %d kept current, %d unchanged of %d residue(s)%s'
            % (n_change, n_keep, n_other, len(outcomes), tail))


def register_rotafit_command(logger):
    from chimerax.core.commands import register, CmdDesc, FloatArg, IntArg, BoolArg
    from chimerax.atomic import ResiduesArg
    desc = CmdDesc(
        optional=[('residues', ResiduesArg)],
        keyword=[('temperature', FloatArg), ('settle_steps', IntArg),
                 ('polish_steps', IntArg), ('polish_top', IntArg),
                 ('ramp_increments', IntArg),
                 ('accept_margin', FloatArg), ('settle_lambda', FloatArg),
                 ('pin_environment', BoolArg), ('pin_k', FloatArg),
                 ('pin_near_cutoff', FloatArg), ('pin_distant_multiplier', FloatArg),
                 ('minimize', BoolArg), ('score_lambda', FloatArg),
                 ('allow_multiple', BoolArg), ('nb_groups', BoolArg),
                 ('apply', BoolArg), ('debug', BoolArg)],
        synopsis='Automated rotamer fit: settle each library rotamer in the running '
                 'simulation and commit the lowest-energy one'
    )
    register('isolde rotafit', desc, rotafit, logger=logger)
