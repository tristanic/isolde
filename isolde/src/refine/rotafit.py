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

# --- scoring mode (see _rank_key / _LocalScorer / _pose_metrics / _choose_winner) --
# 'classic'      : the legacy whole-construct core-FF + MDFF energy (groups 0-5).
#                  Dominated by the O(N^2) ENVIRONMENT nonbonded term, whose value
#                  drifts by hundreds of kJ/mol as the environment settles into
#                  slightly different local minima pose-to-pose -- noise that swamps
#                  the ~tens-of-kJ/mol per-residue signal (the Thr89 failure: a
#                  strained rotamer won by 524 kJ/mol of pure environment jitter).
# 'local'        : rank by the TARGET residue's OWN MDFF map energy only (isolated by
#                  the difference method, so the environment's map contribution
#                  cancels). Atom-local => immune to the environment nonbonded noise.
#                  Fixes Thr89, but on a model-biased crystallographic map it can still
#                  prefer a "peak-parking" rotamer over the correct one (the Ile114
#                  failure: wrong 'tp' fits the 2mFo-DFc sum better than correct 'pt').
# 'local+diff'   : 'local' primary, but among candidates that are TIED on the main map
#                  (within DENSITY_TIE_BAND x global_k of the best), break the tie by a
#                  DIFFERENCE-density (mFo-DFc) score -- the sigma-normalised sum of the
#                  difference map sampled at the moved heavy atoms. The 2mFo-DFc map is
#                  phase-biased toward the CURRENT model (weak where un-modelled atoms
#                  belong); the mFo-DFc map marks exactly that un-modelled density
#                  (positive) and mismodelled density (negative), WITHOUT recomputing
#                  structure factors per candidate. Used ONLY as a tiebreaker (and the
#                  do-no-harm justification), never as an MDFF driver. Degrades to
#                  'local' when no difference map exists (cryo-EM / apo). THE DEFAULT.
SCORE_MODE = 'local+diff'
DENSITY_TIE_BAND = 20.0  # 'local+diff': candidates whose map_local is within this many
                         #   multiples of global_k of the best are treated as tied on
                         #   the main map, and the mFo-DFc tiebreaker decides among them.
DIFF_MARGIN = 1.5        # 'local+diff': a challenger must explain this many sigma more
                         #   difference density (summed over moved heavy atoms) than the
                         #   current conformation before it displaces it (do-no-harm).
SETTLE_LAMBDA = 0.6      # soft-core lambda during the settle. ISOLDE's live default
                         #   (~0.95) is stiff enough that a rotamer seeded into an
                         #   overlap can gain enough force to blow up the local system
                         #   before it relaxes. Softening the wall to ~0.6 lets clashy
                         #   seeds slide apart instead of exploding; lower than ~0.5
                         #   and real geometry starts to suffer. Restored after.

# Per-group soft-core group ids used by rotafit (see _Softener). The environment is
# group 0; the target residue's real atoms are group 1; and -- so that the target's
# genuine internal (real<->real) interactions can stay rigid while its crystal
# self-contact with its own symmetry copies is softened -- the target's SYMMETRY COPIES
# go in a distinct group 2. Needs >= 3 provisioned slots (SimParams.nb_groups_max
# defaults to 4). In a non-symmetry sim group 2 is simply empty.
ENV_NB_GROUP = 0
TARGET_NB_GROUP = 1
TARGET_COPY_NB_GROUP = 2


class _Softener:
    '''The "soften the target against its surroundings" knob (per-group soft-core).

    The target's real atoms go in group 1 and their symmetry copies in group 2, leaving
    the environment (+ its copies) in group 0. ``set(value)`` softens every coupling that
    TOUCHES the target -- vs environment ``(1,0)``/``(2,0)`` and its crystal self-contact
    ``(1,2)``/``(2,2)`` -- while leaving the target's genuine INTERNAL coupling ``(1,1)``
    and the environment's ``(0,0)`` FULL. So the target sidechain stays internally rigid
    and the environment holds its shape WITHOUT position restraints, yet the target can
    slide through a clash against its own symmetry image (which would otherwise explode to
    a NaN). ``get`` reads the representative ``(1,0)`` coupling.

    (This is the same composition as :meth:`SimHandler.soften_nb_selection`, but split so
    the one-time group assignment (:meth:`assign_target`) is separate from the coupling
    level (:meth:`set`), which the polish ramps repeatedly.)'''
    def __init__(self, sh):
        self.sh = sh

    def get(self):
        return self.sh.get_nb_coupling(TARGET_NB_GROUP, ENV_NB_GROUP)

    def set(self, value):
        if value is None:
            return
        # Soften every coupling touching the target -- vs environment and its own crystal
        # image -- but NOT the target's internal (1,1) or the environment's (0,0).
        sh = self.sh
        sh.set_nb_coupling(TARGET_NB_GROUP, ENV_NB_GROUP, value)              # (1,0)
        sh.set_nb_coupling(TARGET_COPY_NB_GROUP, ENV_NB_GROUP, value)         # (2,0)
        sh.set_nb_coupling(TARGET_NB_GROUP, TARGET_COPY_NB_GROUP, value)      # (1,2)
        sh.set_nb_coupling(TARGET_COPY_NB_GROUP, TARGET_COPY_NB_GROUP, value)  # (2,2)

    def assign_target(self, residue):
        '''Target real atoms -> group 1, their symmetry copies -> group 2.'''
        self.sh.assign_nb_group(residue.atoms, TARGET_NB_GROUP,
                                copy_group_id=TARGET_COPY_NB_GROUP)

    def release_target(self, residue):
        '''Return the target (and its copies) to the environment group.'''
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


class _LocalScorer:
    '''Atom-local per-residue scoring, to replace the noise-dominated whole-construct
    energy (see SCORE_MODE). Built ONCE per rotafit run; its per-pose methods then read
    the CURRENT context state. Two signals:
      * ``diff_score(coords)`` -- the sigma-normalised sum of the mFo-DFc DIFFERENCE map
        at the given atoms. This is the tiebreaker that separates the correct rotamer
        from a "peak-parking" one the model-biased main map prefers (see SCORE_MODE
        'local+diff'). None when the model has no difference map.
      * ``target_map_force(residue, ridx)`` -- (rms, coherence) of the per-atom MDFF map
        force on the target. High rms + high coherence (all sidechain atoms pulled the
        same way) = the sidechain is being coherently pulled toward a nearby density
        lobe ("cooperative gradient"); ~0 rms = stranded on the flat map. A DEBUG-only
        diagnostic for the "refuses to move" class (not used for ranking).'''
    def __init__(self, sh, session, model=None):
        self.sh = sh
        self._atoms = sh._atoms
        self._n_real = len(sh._atoms)
        # The mFo-DFc difference map (for the 'local+diff' tiebreaker), if this is a
        # crystallographic model that has one. None for cryo-EM / apo -> tiebreak
        # degrades to plain 'local'.
        self.diffmap = _get_difference_map(model) if model is not None else None
        self._diff_sigma = 1.0
        if self.diffmap is not None:
            try:
                self._diff_sigma = float(self.diffmap.mean_sd_rms()[1]) or 1.0
            except Exception:
                self._diff_sigma = 1.0

    def diff_score(self, coords_local):
        '''Sigma-normalised sum of the mFo-DFc difference map at the given (model-local,
        Angstrom) coordinates -- typically a pose's MOVED heavy atoms. Positive => the
        atoms sit on POSITIVE difference density (density the current model fails to
        explain; the atoms belong here); negative => on negative difference (atoms the
        data does not support). None when there is no difference map. Cheap -- a static
        map interpolation, no structure-factor recalculation.'''
        if self.diffmap is None or coords_local is None or not len(coords_local):
            return None
        try:
            vals = self.diffmap.interpolated_values(np.asarray(coords_local, dtype=float))
            return float(np.sum(vals) / self._diff_sigma)
        except Exception:
            return None

    def target_map_force(self, residue, ridx):
        '''(rms magnitude, coherence in [0,1]) of the per-atom MDFF map force on the
        target. Coherence = |sum of per-atom forces| / sum of |per-atom force|: 1.0 =
        all atoms pulled the same direction (a net translational tug into density),
        ~0 = forces cancel (seated, or incoherently pinned).'''
        from openmm import unit
        from ..openmm.openmm_interface import MAP_FORCE_GROUP
        try:
            st = self.sh._simulation.context.getState(getForces=True,
                                                      groups={MAP_FORCE_GROUP})
            f = st.getForces(asNumpy=True).value_in_unit(
                unit.kilojoule_per_mole / unit.nanometer)[ridx]
            mags = np.linalg.norm(f, axis=1)
            denom = float(np.sum(mags))
            if denom <= 0:
                return 0.0, 0.0
            rms = float(np.sqrt(np.mean(mags ** 2)))
            coher = float(np.linalg.norm(f.sum(axis=0)) / denom)
            return rms, coher
        except Exception:
            return float('nan'), float('nan')


def _get_difference_map(model):
    '''The model's mFo-DFc difference map (a clipper Volume with is_difference_map=True),
    or None if this isn't a crystallographic model with one. Prefers a LIVE difference
    map ('Fo-Fc') over a static one if both are present -- the live map tracks the
    current coordinates. Never raises (returns None on any failure).'''
    if model is None:
        return None
    try:
        from chimerax.clipper.symmetry import get_map_mgr
        mgr = get_map_mgr(model)
        if mgr is None:
            return None
        diffs = [m for m in mgr.all_maps if getattr(m, 'is_difference_map', False)]
        if not diffs:
            return None
        live = [m for m in diffs if 'Fo-Fc' == getattr(m, 'name', '')]
        return live[0] if live else diffs[0]
    except Exception:
        return None


def _map_energy(ctx):
    '''Whole-construct MDFF map energy (kJ/mol) -- the MAP force group total.'''
    from openmm import unit
    from ..openmm.openmm_interface import MAP_FORCE_GROUP
    return ctx.getState(getEnergy=True, groups={MAP_FORCE_GROUP}).getPotentialEnergy(
        ).value_in_unit(unit.kilojoule_per_mole)


def _target_map_energy(sh, isolde, residue):
    '''The target residue's OWN MDFF map energy (kJ/mol), by the DIFFERENCE method:
    read the whole map energy, temporarily disable the target's per-atom map terms
    (across every map), re-read, and subtract. The map term is a pure per-atom sum, so
    the (unchanged) environment contribution cancels exactly and only the target's
    density fit remains. Symmetry-safe: drives the same ``update_mdff_atoms`` path used
    for live edits (which fans a disable out to each atom's symmetry terms too), and it
    resolves immediately because the sim is paused. Returns 0.0 if the sim has no maps.

    NB lower (more negative) = better fit; the target sits in density when this is a
    large negative number, and ~0 when it is out of density.'''
    sm = getattr(isolde, 'sim_manager', None)
    mgrs = getattr(sm, 'mdff_mgrs', None) if sm is not None else None
    if not mgrs:
        return 0.0
    ctx = sh._simulation.context
    e_all = _map_energy(ctx)
    saved = []
    try:
        for mgr in mgrs.values():
            m = mgr.get_mdff_atoms(residue.atoms)
            if m is None or len(m) == 0:
                continue
            saved.append((mgr, m, np.array(m.enableds, copy=True)))
            m.enableds = np.zeros(len(m), dtype=bool)   # disable the target's map terms
            sh.update_mdff_atoms(m, mgr.volume)    # live (paused => pushed at once)
        e_without = _map_energy(ctx)
    finally:
        for mgr, m, en in saved:                   # always restore
            m.enableds = en
            sh.update_mdff_atoms(m, mgr.volume)
    return e_all - e_without


def _pose_metrics(sh, isolde, scorer, residue, ridx, score_mode=SCORE_MODE, debug=False):
    '''Scoring signals for the pose currently loaded in the context (settled state, no
    stepping). Computes ONLY what the given ``score_mode`` needs to rank, PLUS -- only in
    ``debug`` -- the extra diagnostic terms for the side-by-side log. Missing terms are
    filled with NaN / None so the log formatter stays happy.'''
    ctx = sh._simulation.context
    nan = float('nan')
    m = {'whole': nan, 'map_whole': nan, 'diff_score': None,
         'mapforce_rms': nan, 'mapforce_coher': nan}
    # map_local is the primary key for every 'local*' mode -> always needed.
    m['map_local'] = _target_map_energy(sh, isolde, residue)
    if score_mode == 'classic' or debug:
        m['whole'] = _score_energy(ctx)            # rank key for classic; else diagnostic
    if scorer is not None and (score_mode == 'local+diff' or debug):
        # Difference-density over the target's HEAVY atoms (cheap: one map interpolation).
        # The backbone is identical across a residue's rotamers, so it contributes a
        # constant that cancels in every ranking / do-no-harm comparison.
        heavy_rows = ridx[residue.atoms.element_names != 'H']
        m['diff_score'] = scorer.diff_score(_sim_coords(sh)[heavy_rows])
    if debug:                                      # cheap-ish diagnostics for the log only
        m['map_whole'] = _map_energy(ctx)
        if scorer is not None:
            m['mapforce_rms'], m['mapforce_coher'] = scorer.target_map_force(residue, ridx)
    return m


def _pre_settle_metrics(sh, isolde, scorer, residue, ridx):
    '''Cheap metrics on the RIGID (pre-relax) placement: the target's own map fit and
    the (rms, coherence) of the map force on it -- the "approach" signal only.
    Coherence -> 1 means the whole sidechain is
    being tugged the same way (a cooperative pull toward a nearby density lobe); ~0 means
    the per-atom pulls cancel (already in place, or incoherently splayed across peaks).'''
    m = {'pre_map_local': _target_map_energy(sh, isolde, residue)}
    if scorer is not None:
        m['pre_mapf_rms'], m['pre_mapf_coher'] = scorer.target_map_force(residue, ridx)
    else:
        m['pre_mapf_rms'] = m['pre_mapf_coher'] = float('nan')
    return m


def _rank_key(metrics, score_mode):
    '''The scalar a pose is ranked by (smaller = better). 'classic' = whole-construct
    energy; 'local'/'local+diff' = the target's own map fit (the difference-density
    tiebreak in 'local+diff' is applied across the candidate set in _choose_winner, not
    in this per-pose key).'''
    if score_mode == 'classic':
        return metrics['whole']
    return metrics['map_local']


def _metrics_str(metrics):
    '''Compact one-line dump of the scoring signals, for the debug log. Pre-settle
    fields (pre_*, disp) are only present for the soft-search candidates, so read them
    defensively -- the polished dumps omit them.'''
    g = metrics.get
    ds = g('diff_score', None)
    ds_str = ('%.2f' % ds) if ds is not None else 'n/a'
    s = ('map_local %.1f | diff %s | map_whole %.0f | whole %.0f'
         % (metrics['map_local'], ds_str, metrics['map_whole'], metrics['whole']))
    if 'pre_map_local' in metrics:
        s += (' || pre: map %.1f mapF rms %.0f coher %.2f | disp %.2fA'
              % (g('pre_map_local', float('nan')), g('pre_mapf_rms', float('nan')),
                 g('pre_mapf_coher', float('nan')), g('disp', float('nan'))))
    return s


def _density_tie_band(results, gk):
    '''The set of candidates TIED with the best on the main (2mFo-DFc) map: those whose
    map_local is within DENSITY_TIE_BAND x global_k of the best. These are the poses the
    main map cannot separate -- where the mFo-DFc difference tiebreaker earns its keep.'''
    band = DENSITY_TIE_BAND * gk if gk > 0 else DENSITY_TIE_BAND
    best_map = min(t[4]['map_local'] for t in results)
    return band, [t for t in results if t[4]['map_local'] <= best_map + band]


def _choose_winner(candidates, score_mode, eff_margin, gk, dlog):
    '''Pick the pose to commit from the ranked (best-first by map_local) candidates.
    Returns (rank_key, name, coords, kept_current).

    'local+diff': among the map-tied candidates (:func:`_density_tie_band`), the pose
    that best explains the mFo-DFc DIFFERENCE density wins -- this is what separates the
    correct rotamer from a "peak-parking" one the model-biased main map prefers. The
    current conformation is displaced only if the winner explains >= DIFF_MARGIN sigma
    more difference density than it does (do-no-harm, but arbitrated by the data the main
    map is biased against). Other modes: the plain map/energy do-no-harm margin.'''
    best = candidates[0]
    curr = next((t for t in candidates if t[3]), None)
    have_diff = any(t[4].get('diff_score') is not None for t in candidates)
    use_diff = (score_mode == 'local+diff' and have_diff)
    if use_diff:
        band, tie = _density_tie_band(candidates, gk)
        tie = [t for t in tie if t[4].get('diff_score') is not None]
        if len(tie) > 1:
            best = max(tie, key=lambda t: t[4]['diff_score'])
            dlog('  diff tiebreak (map band %.1f kJ/mol): %s -> "%s"'
                 % (band, ', '.join('%s %.2f' % (t[1], t[4]['diff_score']) for t in tie),
                    best[1]))
    best_e, best_name, best_coords = best[:3]
    kept_current = best[3]
    if curr is not None and not kept_current:
        cd = curr[4].get('diff_score')
        bd = best[4].get('diff_score')
        if use_diff and cd is not None and bd is not None:
            gain = bd - cd
            if gain < DIFF_MARGIN:
                best_e, best_name, best_coords, kept_current = (*curr[:3], True)
                dlog('  do-no-harm: winner explains only %.2f sigma more diff density '
                     'than current (< %.2f) -> keeping current' % (gain, DIFF_MARGIN))
        elif (curr[0] - best[0]) < eff_margin:
            best_e, best_name, best_coords, kept_current = (*curr[:3], True)
            dlog('  do-no-harm: winner beats current by <%.2f kJ/mol -> keeping current'
                 % eff_margin)
    return best_e, best_name, best_coords, kept_current


def _commit_best(sh, isolde, softener, scorer, score_mode, residue, results, base,
                 polish_steps, polish_top, ramp_increments, saved_lambda, settle_lambda,
                 eff_margin, gk, outcomes, dlog, debug, minimize=False, score_lambda=0.0):
    '''Polish the top candidate(s) at full stiffness, re-rank on the same ``score_mode``
    key, apply the do-no-harm margin against the current conformation, and commit the
    winner. Appends the per-residue outcome (change/keep) to ``outcomes``.'''
    # Polish the top candidates at full soft-core stiffness and RE-RANK by the polished
    # score: the soft (search) lambda and the stiff (live) lambda can disagree, so
    # re-ranking guards against a soft-search tie committing the wrong basin. polish_top:
    # N = top N; <=0 = all survivors; 1 = winner only. ALWAYS polish the current
    # conformation too (do-no-harm needs it judged at the same lambda as its rivals).
    candidates = results
    if polish_steps and saved_lambda is not None:
        if score_mode == 'local+diff':
            # Polish the map-tied candidates -- but only the DIFFERENCE-density
            # CONTENDERS among them (search diff within DIFF_MARGIN of the best band
            # diff). A band member with low difference density can neither win the diff
            # tiebreak nor displace current, so minimising it to full stiffness is wasted
            # work. (The winner is nearly always NOT the top-ranked-by-map pose, so we do
            # still polish more than one.) current is re-added below for do-no-harm.
            _band, band_members = _density_tie_band(results, gk)
            band_diffs = [t[4].get('diff_score') for t in band_members
                          if t[4].get('diff_score') is not None]
            if band_diffs:
                cut = max(band_diffs) - DIFF_MARGIN
                polish_set = [t for t in band_members if t[3]
                              or (t[4].get('diff_score') is not None
                                  and t[4]['diff_score'] >= cut)]
            else:
                polish_set = list(band_members)
        else:
            n_polish = (len(results) if polish_top <= 0
                        else min(polish_top, len(results)))
            polish_set = list(results[:n_polish])
        if not any(t[3] for t in polish_set):
            curr = next((t for t in results if t[3]), None)
            if curr is not None:
                polish_set.append(curr)
        ctx = sh._simulation.context
        integrator = sh._main_integrator
        ridx = sh._atoms.indices(residue.atoms)
        from time import perf_counter
        t_relax = t_metrics = 0.0
        polished = []
        for _k, pname, pcoords, pcur, _pm in polish_set:
            tr = perf_counter()
            _pe, pc = _ramped_polish(sh, softener, integrator, ctx, ridx, pcoords,
                                     pcoords[ridx], polish_steps, settle_lambda,
                                     saved_lambda, ramp_increments, minimize, score_lambda)
            t_relax += perf_counter() - tr
            tm = perf_counter()
            metrics = _pose_metrics(sh, isolde, scorer, residue, ridx, score_mode, debug)
            t_metrics += perf_counter() - tm
            polished.append((_rank_key(metrics, score_mode), pname, pc, pcur, metrics))
            if debug:                          # every signal, to see what discriminates
                dlog('    %-9s %s' % (pname, _metrics_str(metrics)))
        softener.set(settle_lambda)            # back to soft for the next residue
        if debug:
            dlog('  polish timing: relax %.2fs + metrics %.2fs (%d polished)'
                 % (t_relax, t_metrics, len(polished)))
        polished.sort(key=lambda t: t[0])
        candidates = polished
        plist = ', '.join('%s %.1f' % (n, k) for k, n, _c, _cur, _m in polished)
        dlog('  polished %d (lambda %.2f->%.2f x%d) ranked by %s: %s'
             % (len(polished), settle_lambda, saved_lambda, ramp_increments,
                score_mode, plist))
    # Pick the winner: 'local+diff' breaks map-ties by the mFo-DFc difference score and
    # applies a data-driven do-no-harm; other modes use the plain map/energy margin
    # ("keep current unless a rival beats it by >= eff_margin", else it's settling noise).
    best_e, best_name, best_coords, kept_current = _choose_winner(
        candidates, score_mode, eff_margin, gk, dlog)
    dlog('  -> commit "%s" (rank=%.1f)' % (best_name, best_e))
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
    # Ramp by STEPPING at each increment (cheap, and what makes the stiffening adiabatic);
    # run the expensive minimise ONCE, at the final full-stiffness increment, to converge.
    # (Minimising at every increment was ~5x the cost for no ranking benefit -- the
    # intermediate minima are thrown away.)
    last = len(lambdas) - 1
    for i, lam in enumerate(lambdas):
        softener.set(float(lam))
        _relax(sh, integrator, ctx, steps_per, minimize and i == last)
    return _score_energy(ctx, softener, score_lambda), _sim_coords(sh)


def _settle_and_rank(isolde, residue, temperature, settle_steps, log, minimize=False,
                     scorer=None, score_mode=SCORE_MODE, debug=False):
    '''Settle each (non-culled) library rotamer in the running sim and return the
    ranked results [(rank_key, name, settled_construct_coords, is_current, metrics), ...]
    best-first, plus the original construct coords. ``rank_key`` is the pose's score
    under ``score_mode`` (see _rank_key); ``metrics`` carries the scoring signals for the
    debug log and the difference tiebreak. The winner is polished separately at full
    stiffness. Sim paused, main thread.'''
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
    from time import perf_counter
    t_relax = t_metrics = 0.0
    for name, rcoords in poses:
        is_current = (name == CURRENT_LABEL)
        if not is_current and _severe_clash(rcoords, residue, moved_names, env_coords):
            n_cull += 1
            continue
        coords = base.copy()
        coords[ridx] = rcoords
        _push_now(sh, coords)
        # The rigid-placement "approach" signals (pre-settle map fit + map-force
        # coherence + travel distance) are diagnostics only -- compute them just in
        # debug (the pre-settle map fit uses the difference-trick toggling, which isn't
        # free). A settled endpoint hides the approach anyway.
        pre = _pre_settle_metrics(sh, isolde, scorer, residue, ridx) if debug else {}
        # SEARCH is dynamics-ONLY (no minimise), regardless of the polish `minimize`
        # setting: the search only has to nominate the density tie-band, and the polish
        # phase does the honest minimised ranking. 0 K dynamics seats each pose in its
        # basin well enough to rank for band selection, and the generous band absorbs the
        # extra jitter -- deleting 1 whole-construct minimise PER POSE from the hot path.
        tr = perf_counter()
        _relax(sh, integrator, ctx, settle_steps, minimize=False)
        t_relax += perf_counter() - tr
        settled = _sim_coords(sh)
        tm = perf_counter()
        metrics = _pose_metrics(sh, isolde, scorer, residue, ridx, score_mode, debug)
        t_metrics += perf_counter() - tm
        metrics.update(pre)
        if debug:
            metrics['disp'] = float(np.sqrt(np.mean(
                np.sum((settled[ridx] - rcoords) ** 2, axis=1))))
        results.append((_rank_key(metrics, score_mode), name, settled, is_current,
                        metrics))
    results.sort(key=lambda t: t[0])
    log('isolde rotafit %s: %d rotamers (+current), %d culled (severe clash), %d settled'
        % (residue, len(poses) - 1, n_cull, len(results)))
    if debug:
        log('  search timing: relax %.2fs + metrics %.2fs (%d poses)'
            % (t_relax, t_metrics, len(results)))
    return results, base


def rotafit(session, residues=None, temperature=0.0, settle_steps=SETTLE_STEPS,
            polish_steps=POLISH_STEPS, polish_top=POLISH_TOP,
            ramp_increments=RAMP_INCREMENTS, accept_margin=ACCEPT_MARGIN,
            settle_lambda=SETTLE_LAMBDA,
            minimize=False, score_lambda=0.0, score_mode=SCORE_MODE,
            allow_multiple=False, apply=True, debug=False):
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
    # rotafit drives per-group soft-core coupling (target soft against its surroundings
    # -- including its own symmetry copies -- while the environment stays rigid without
    # pinning). This needs >= 3 group slots (env=0, target=1, target-copy=2); every sim
    # is built with them by default (SimParams.nb_groups_max defaults to 4), and it is
    # fully available unless the user disabled it or the sim can't support it (see
    # SimHandler.nb_groups_enabled).
    if not (getattr(sh, 'nb_groups_enabled', False)
            and getattr(sh, 'nb_groups_count', 1) >= 3):
        raise UserError(
            'isolde rotafit requires per-group soft-core coupling, which is not '
            'available for this simulation. It needs the soft-core nonbonded potential '
            'and SimParams.nb_groups_max >= 3 (the default is 4) set before the '
            'simulation was started. Restart the simulation with the defaults, or '
            'raise nb_groups_max, then retry.')
    params = dict(temperature=temperature, settle_steps=settle_steps,
                  polish_steps=polish_steps, polish_top=polish_top,
                  ramp_increments=ramp_increments, accept_margin=accept_margin,
                  settle_lambda=settle_lambda,
                  minimize=minimize, score_lambda=score_lambda, score_mode=score_mode,
                  apply=apply, debug=debug)
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
                 minimize, score_lambda, score_mode, apply, debug, resume_after):
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
    softener = _Softener(sh)
    # The target's current coupling to its surroundings (full = 1.0 initially); the
    # environment is held by its own full-strength force field -- NO position pinning.
    saved_lambda = softener.get()
    dlog('isolde rotafit: per-group soft-core coupling (environment unpinned)')
    # Per-residue LOCALISED scoring (see SCORE_MODE): the difference map is looked up
    # once from the (shared) target structure -- rotafit is single-residue by default,
    # so one model.
    model = targets[0].structure if targets else None
    scorer = _LocalScorer(sh, session, model=model)
    dlog('isolde rotafit: scoring mode "%s"%s' % (score_mode,
         '' if scorer.diffmap is not None else ' (no difference map: diff tiebreak off)'))
    # Scale the do-no-harm margin by ISOLDE's MDFF coupling constant so it tracks the
    # map's energy scale (x-ray vs high-sigma cryo-EM) instead of being a fixed kJ/mol.
    gk = _sim_coupling_constant(isolde)
    eff_margin = accept_margin * gk if gk > 0 else accept_margin
    dlog('isolde rotafit: MDFF coupling global_k=%.3g -> do-no-harm margin %.2f kJ/mol'
         % (gk, eff_margin))
    outcomes = []                              # (residue, kind, name) per target
    try:
        softener.set(settle_lambda)            # soften the target against its surroundings
                                               # (incl. its symmetry copies) so a clashy
                                               # seed relaxes; the environment stays rigid
        for r in targets:
            softener.assign_target(r)          # target -> group 1, its copies -> group 2
            try:
                from time import perf_counter
                t0 = perf_counter()
                try:
                    results, base = _settle_and_rank(isolde, r, temperature,
                                                     settle_steps, dlog, minimize,
                                                     scorer, score_mode, debug)
                except Exception as e:
                    session.logger.warning('isolde rotafit: skipped %s (%s)' % (r, e))
                    outcomes.append((r, 'skip', None))
                    continue
                t_search = perf_counter() - t0
                if not results:
                    _push_now(sh, base)        # nothing viable -> restore
                    outcomes.append((r, 'none', None))
                    continue
                shortlist = ', '.join('%s %.1f' % (n, k)
                                      for k, n, _c, _cur, _m in results[:6])
                dlog('  search ranking @ lambda=%.2f (by %s): %s'
                     % (settle_lambda, score_mode, shortlist))
                if debug:                       # full signal dump for every candidate
                    for _k, n, _c, _cur, m in results:
                        dlog('    %-9s %s' % (n, _metrics_str(m)))
                if not apply:
                    _push_now(sh, base)         # nothing committed -> restore
                    outcomes.append((r, 'noapply', None))
                    continue
                t1 = perf_counter()
                _commit_best(sh, isolde, softener, scorer, score_mode, r, results, base,
                             polish_steps, polish_top, ramp_increments, saved_lambda,
                             settle_lambda, eff_margin, gk, outcomes, dlog, debug,
                             minimize, score_lambda)
                dlog('  timing: search %.2fs (%d poses), polish+commit %.2fs'
                     % (t_search, len(results), perf_counter() - t1))
            finally:
                softener.release_target(r)     # target (and its copies) back to group 0
    finally:
        sh.temperature = saved_temp            # don't leave the sim at 0 K
        softener.set(saved_lambda)             # restore full coupling before resume
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
    from chimerax.core.commands import (register, CmdDesc, FloatArg, IntArg, BoolArg,
                                        EnumOf)
    from chimerax.atomic import ResiduesArg
    desc = CmdDesc(
        optional=[('residues', ResiduesArg)],
        keyword=[('temperature', FloatArg), ('settle_steps', IntArg),
                 ('polish_steps', IntArg), ('polish_top', IntArg),
                 ('ramp_increments', IntArg),
                 ('accept_margin', FloatArg), ('settle_lambda', FloatArg),
                 ('minimize', BoolArg), ('score_lambda', FloatArg),
                 ('score_mode',
                  EnumOf(('local+diff', 'local', 'classic'))),
                 ('allow_multiple', BoolArg),
                 ('apply', BoolArg), ('debug', BoolArg)],
        synopsis='Automated rotamer fit: settle each library rotamer in the running '
                 'simulation and commit the lowest-energy one'
    )
    register('isolde rotafit', desc, rotafit, logger=logger)
