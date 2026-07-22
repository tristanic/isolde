# @Author: Tristan Croll
# @Date:   20-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 20-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
settle_poses -- settle a set of CALLER-SUPPLIED candidate poses for one residue in
the running ISOLDE simulation, rank them by the pose's own localised packing energy,
and commit the best.

This is the ``isolde rotafit`` engine with the candidate generator swapped from
"enumerate library rotamers" to "the caller injects poses", plus a packing-first
scoring policy. Its first consumer is a downstream structure-based-drug-design bundle
doing fragment elaboration: it modifies a ligand, generates a handful of 3D poses, and
asks ISOLDE which one packs best against the receptor under the live force field.

Public API (function only -- poses are not command-line typeable):

    from chimerax.isolde.refine import settle_poses, SettleResult
    result = settle_poses(session, residue, poses, moved_atoms=..., ...)

Scoring -- the hard-won lesson from the rotafit rework -- is the fragment's OWN
LOCALISED packing energy, never the whole-construct energy. The whole-construct core
energy is dominated by the O(N^2) environment nonbonded term, which drifts by hundreds
of kJ/mol pose-to-pose as the surroundings settle into slightly different minima; that
noise swamps the per-fragment signal. Instead we read

    E_pack = E_core(fragment fully coupled) - E_core(fragment decoupled to eps)

over the full CORE force groups at fixed (settled) coordinates. Everything that does not
change with the fragment's coupling -- env-env nonbonded, fragment-internal, all bonded
terms, implicit-solvent self-energy -- cancels in the difference, leaving just the
fragment<->surroundings interaction. Lower is better; positive means net clash.

The map (when present) is NOT a ligand-scoring term. In fragment elaboration there is no
experimental density for the ligand, so with ``map_decouple=True`` (default) the ligand's
own map terms are switched off for the fit while the receptor stays map-restrained toward
its experimental geometry; the ligand is judged on packing alone.

Two drivers, selected by ``live``:
  * fast burst (default): requires a PAUSED sim; drives the integrator directly on the
    main thread in one synchronous pass and returns a SettleResult.
  * live (``live=True``): runs the poses through ISOLDE's canonical GUI loop so they can
    be watched settling (debug/demo); returns immediately and delivers the SettleResult
    via ``on_done(result)``.
'''

import numpy as np
from contextlib import contextmanager, nullcontext
from dataclasses import dataclass
from typing import Optional, List, Tuple

from .settle_common import severe_clash

# Soft-core coupling used to "switch off" the fragment's interaction with its
# surroundings when isolating E_pack. NOT zero: the soft-core potential rejects
# lambda <= 0 (custom_forces.py set_coupling: "0 < lambda <= 1"), and the residual
# fragment<->env energy at 1e-4 is negligible, so the difference is clean.
NB_DECOUPLE_EPS = 1e-4
# Soft-core coupling of the fragment to its surroundings DURING the settle (not the
# scoring read). A soft wall lets a clashy input pose slide apart instead of exploding;
# the polish then ramps back to full stiffness. Mirrors rotafit's SETTLE_LAMBDA.
SETTLE_LAMBDA = 0.6
SETTLE_STEPS = 100       # 0 K dynamics steps per pose during the soft search
POLISH_STEPS = 150       # 0 K steps applied while ramping a polished pose to full stiffness
RAMP_INCREMENTS = 5      # stages over which the polish ramps SETTLE_LAMBDA -> full
CURRENT_LABEL = '(current)'   # label given to the residue's starting conformation

# Force groups the fragment's coupling toggle actually changes: soft-core vdW is in the
# NONBONDED group, but the nb-group-aware GBSA force is in the DEFAULT group, so read the
# WHOLE CORE set -- everything static cancels in the coupled-minus-decoupled difference.
# (Resolved at call time from the engine module so the group ids stay single-sourced.)


@dataclass
class SettleResult:
    '''Outcome of a :func:`settle_poses` run.

    * committed: label of the pose whose coordinates were committed, or None if nothing
      was applied (``apply=False``) or the run failed.
    * kept_current: True if the do-no-harm margin kept the residue's starting conformation.
    * applied: True if coordinates were written to the model/sim.
    * energies: ``[(label, E_pack_kJ_per_mol), ...]`` best-first (lower = better packing).
    * n_culled: poses dropped by the severe-clash pre-filter before settling.
    * accept_margin_kJ: the effective do-no-harm margin used (kJ/mol).
    * committed_coords: the committed residue-atom coordinates ``(N,3)`` (Angstrom), or None.
    * map_decoupled / score_mode: echoes of the run configuration.
    * metrics: optional extra per-pose diagnostics when ``debug=True``.
    '''
    committed: Optional[str]
    kept_current: bool
    applied: bool
    energies: List[Tuple[str, float]]
    n_culled: int
    accept_margin_kJ: float
    committed_coords: Optional[np.ndarray] = None
    map_decoupled: bool = True
    score_mode: str = 'E_pack'
    metrics: Optional[dict] = None


def settle_poses(session, residue, poses, *, moved_atoms=None,
                 map_decouple=True, soft_group=True,
                 temperature=0.0, settle_steps=SETTLE_STEPS, minimize=False,
                 polish_top=1, accept_margin=1.0, apply=True,
                 live=False, on_done=None, debug=False):
    '''Settle ``poses`` for ``residue`` in the running simulation and commit the best.

    Args:
        * session: the ChimeraX session (ISOLDE must be started, GUI mode).
        * residue: the ``Residue`` whose atoms the poses describe. Must be mobile in the
          running simulation.
        * poses: list of ``(label, coords)``; ``coords`` is an ``(N,3)`` array in
          ``residue.atoms`` order (Angstrom), ``N == len(residue.atoms)``. The residue's
          current conformation is added automatically as ``'(current)'``.
        * moved_atoms: the subset of ``residue.atoms`` the poses actually change; defaults
          to all of ``residue.atoms``. Drives the clash pre-filter and the soft-group.
        * map_decouple: if True (default), switch the ligand's own MDFF map terms off for
          the fit so it is judged on packing alone while the receptor stays map-restrained.
        * soft_group: if True (default), soften the fragment against its surroundings while
          settling (clashy poses relax instead of exploding). Scoring always uses per-group
          coupling regardless, so per-group coupling is required either way.
        * temperature: sim temperature (K) held for the fit; restored afterward. 0 = a
          deterministic 0 K descent (recommended).
        * settle_steps: 0 K dynamics steps per pose during the soft search.
        * minimize: also energy-minimise each polished pose (default False; 0 K settling is
          usually enough and cheaper).
        * polish_top: how many top search survivors (plus current) to re-settle at full
          stiffness and re-rank by the polished E_pack.
        * accept_margin: do-no-harm stickiness. The current conformation is displaced only
          if a challenger's E_pack beats it by at least this margin; scaled by ISOLDE's
          MDFF coupling constant when a map drives the sim, else an absolute kJ/mol value.
        * apply: if True (default), commit the winning pose; if False, score only.
        * live: if False (default), fast synchronous burst -- REQUIRES a paused sim, returns
          a SettleResult. If True, run through the GUI loop (watchable); returns None and
          delivers the SettleResult via ``on_done(result)``.
        * on_done: callback ``on_done(SettleResult)`` for ``live=True``.
        * debug: verbose per-pose logging + a metrics dict on the result.

    Returns a :class:`SettleResult` (``live=False``) or None (``live=True``; result via
    ``on_done``). Raises ``UserError`` on a precondition failure.
    '''
    from chimerax.core.errors import UserError
    isolde = getattr(session, 'isolde', None)
    if isolde is None or not session.ui.is_gui:
        raise UserError('settle_poses requires the ISOLDE GUI. Run "isolde start" first.')
    if residue is None or residue.deleted or len(residue.atoms) == 0:
        raise UserError('settle_poses: residue is missing, deleted, or has no atoms.')

    ratoms = residue.atoms
    n = len(ratoms)
    poses = _validate_poses(poses, n, UserError)
    if moved_atoms is None:
        moved_atoms = ratoms
    elif len(moved_atoms.subtract(ratoms)):
        raise UserError('settle_poses: moved_atoms must be a subset of residue.atoms.')

    if not isolde.simulation_running:
        raise UserError('settle_poses needs a running ISOLDE simulation covering the '
                        'residue. Start one (e.g. "isolde sim start sel") first.')
    sh = isolde.sim_handler

    # The residue must be part of the mobile simulation construct.
    ridx = sh.atoms.indices(ratoms)
    if np.any(ridx < 0):
        raise UserError('settle_poses: the residue is not part of the running simulation '
                        '(all of its atoms must be mobile). Start a simulation that covers '
                        'it, e.g. "isolde sim start sel".')

    # Per-group soft-core coupling is required for the localised E_pack scoring (and the
    # soft settle): env=0, fragment=1, its symmetry copies=2 -> needs >= 3 slots. Every sim
    # provides these by default (SimParams.nb_groups_max=4). No pin fallback.
    if not (getattr(sh, 'nb_groups_enabled', False)
            and getattr(sh, 'nb_groups_count', 1) >= 3):
        raise UserError(
            'settle_poses requires per-group soft-core coupling, which is not available '
            'for this simulation. It needs the soft-core nonbonded potential and '
            'SimParams.nb_groups_max >= 3 (the default is 4) set before the simulation '
            'started. Restart the simulation with the defaults, or raise nb_groups_max, '
            'then retry. (Known gap: crystallographic symmetry + GBSA is not yet '
            'group-aware.)')

    cfg = dict(moved_atoms=moved_atoms, map_decouple=map_decouple, soft_group=soft_group,
               temperature=float(temperature), settle_steps=int(settle_steps),
               minimize=bool(minimize), polish_top=int(polish_top),
               accept_margin=float(accept_margin), apply=bool(apply), debug=bool(debug))

    if live:
        return _run_settle_live(session, isolde, residue, poses, on_done=on_done, **cfg)
    return _run_settle_burst(session, isolde, residue, poses, **cfg)


def _validate_poses(poses, n, UserError):
    '''Coerce/validate the poses list: each entry ``(label, (n,3) coords)``. Returns a list
    of ``(str_label, float64 (n,3) array)``. Raises with the offending label on any bad
    shape so the caller can find it.'''
    if poses is None or len(poses) == 0:
        raise UserError('settle_poses: no poses supplied.')
    out = []
    for i, entry in enumerate(poses):
        try:
            label, coords = entry
        except (TypeError, ValueError):
            raise UserError('settle_poses: pose %d is not a (label, coords) pair.' % i)
        arr = np.asarray(coords, dtype=float)
        if arr.shape != (n, 3):
            raise UserError('settle_poses: pose "%s" has coordinates of shape %s, but the '
                            'residue has %d atoms (expected (%d, 3)).'
                            % (label, arr.shape, n, n))
        out.append((str(label), arr))
    return out


# ------------------------------------------------------------------
# Fit-scope guards (wrap the WHOLE fit, both drivers): force 0 K and suspend live
# crystallographic-map recalculation, restoring both on exit.
# ------------------------------------------------------------------
def _configured_temperature_k(isolde, sh):
    '''The sim's CONFIGURED target temperature in Kelvin (the value to restore to after a
    0 K fit). Read from sim_params (unchanged by the fit); falls back to the live handler
    temperature. Robust to a Quantity or a bare float.'''
    try:
        st = isolde.sim_params.temperature
        from openmm import unit
        return st.value_in_unit(unit.kelvin) if hasattr(st, 'value_in_unit') else float(st)
    except Exception:
        try:
            return float(sh.temperature)
        except Exception:
            return None


@contextmanager
def _fixed_temperature(isolde, sh, value=0.0):
    '''Hold the sim temperature at ``value`` (default 0 K) for the block, restoring the
    configured temperature on exit. A 0 K hold makes the settle a deterministic damped
    descent (no thermostat noise) -- essential for the live driver, which otherwise runs
    at the user's running temperature.'''
    saved = _configured_temperature_k(isolde, sh)
    try:
        sh.temperature = value
        yield
    finally:
        if saved is not None:
            sh.temperature = saved


@contextmanager
def _suspended_live_map_recalc(model):
    '''Suspend live crystallographic-map (2mFo-DFc / mFo-DFc) recalculation for the block,
    restoring each map set's prior state on exit. Without this the canonical loop would
    recompute maps under the moving poses (the receptor's map restraint would chase a
    moving target). Per-set save/restore because the XmapSet.live_update setter can't report
    its prior value. A no-op when there is no live crystallographic map (cryo-EM / apo).'''
    saved = []
    try:
        from chimerax.clipper.symmetry import get_map_mgr
        mgr = get_map_mgr(model)
        if mgr is not None:
            for xs in mgr.xmapsets:
                # Only sets backed by experimental data can recalc; skip the rest.
                if getattr(xs, 'live_xmap_mgr', None) is not None:
                    saved.append((xs, xs.live_update))
                    xs.live_update = False
    except Exception:
        saved = []
    try:
        yield
    finally:
        for xs, prev in saved:
            try:
                xs.live_update = prev
            except Exception:
                pass


# ------------------------------------------------------------------
# Fast synchronous burst driver.
# ------------------------------------------------------------------
def _run_settle_burst(session, isolde, residue, poses, *, moved_atoms, map_decouple,
                      soft_group, temperature, settle_steps, minimize, polish_top,
                      accept_margin, apply, debug):
    '''Fast path: the sim must be PAUSED on entry (sim thread idle => safe to drive
    _main_integrator on the main thread). Runs the whole settle synchronously and returns
    a SettleResult. The sim is LEFT RUNNING (paused) so a poor result is one revert away.'''
    from chimerax.core.errors import UserError
    sh = isolde.sim_handler
    if not sh.pause:
        raise UserError('settle_poses (fast mode) needs the simulation PAUSED. Pause it '
                        'first ("isolde sim pause"), or call with live=True to run it '
                        'through the GUI loop.')
    isolde.checkpoint()   # so the user can revert if unhappy
    log = session.logger.info
    dlog = log if debug else (lambda *a, **k: None)
    # Guard the whole fit: 0 K descent + no live-map recompute under the moving poses.
    tctx = _fixed_temperature(isolde, sh, temperature)
    mctx = _suspended_live_map_recalc(residue.structure)
    # map_decouple: switch the LIGAND's own map terms off for the fit (receptor stays
    # map-restrained). map_decoupled() is a no-op if the ligand carries no map coupling.
    lctx = sh.map_decoupled(residue.atoms) if map_decouple else nullcontext()
    with tctx, mctx, lctx:
        result = _settle_and_commit(
            session, isolde, residue, poses, moved_atoms=moved_atoms,
            map_decouple=map_decouple, soft_group=soft_group, settle_steps=settle_steps,
            minimize=minimize, polish_top=polish_top, accept_margin=accept_margin,
            apply=apply, debug=debug, dlog=dlog)
    log(_summary_line(residue, result))
    return result


def _run_settle_live(session, isolde, residue, poses, *, on_done=None, **cfg):
    '''Live path: run the poses through ISOLDE's canonical GUI loop so they can be watched
    settling on screen, then score + commit at a single pause at the end (reusing the burst
    scoring). Returns None immediately; the SettleResult is delivered via ``on_done(result)``
    when the run completes.'''
    _LiveSettleRunner(session, isolde, residue, poses, cfg, on_done).start()
    return None


class _LiveSettleRunner:
    '''Drives a settle_poses run through ISOLDE's live GUI loop (watchable), for debugging
    and demonstration. Each pose is pushed into the RUNNING (0 K-forced) simulation and left
    to settle over a few rendered GUI updates -- visible on screen -- then its settled coords
    are snapshotted. After the last pose the sim is paused ONCE and the snapshots are scored
    (E_pack) + polished + committed exactly as the fast driver does, so both drivers agree on
    the same inputs. Slower than the burst driver; its purpose is to *show* what happens.

    Scoring is deliberately NOT done per-pose live: E_pack needs two energy reads (fragment
    coupled vs decoupled) at the SAME coordinates, but while the sim runs a coupling edit only
    takes effect a frame later and the coordinates are moving between reads. Doing all scoring
    at one end-of-run pause (fixed coords, immediate coupling toggles) keeps it clean and uses
    exactly one pause transition.

    The whole fit is wrapped (via an ExitStack held across the async run) in the same guards
    as the burst driver: 0 K temperature and suspended live crystallographic-map recalc, both
    restored on completion, plus (when ``map_decouple``) the ligand's map terms switched off.'''
    FRAMES_PER_POSE = 8   # GUI updates spent visibly settling each pose before snapshotting

    def __init__(self, session, isolde, residue, poses, cfg, on_done):
        self.session = session
        self.isolde = isolde
        self.residue = residue
        self.poses = poses
        self.cfg = cfg
        self.on_done = on_done
        self.sh = isolde.sim_handler
        self.dlog = session.logger.info if cfg.get('debug') else (lambda *a, **k: None)
        self._stack = None
        self._prep = None
        self._queue = []          # [(label, pose_coords, is_current), ...]
        self._settled = []        # [(label, settled_construct_coords, is_current), ...]
        self._n_culled = 0
        self._pose_i = -1
        self._cur_label = None
        self._cur_is_current = False
        self._frames = 0
        self._was_paused = False
        self._cu_handler = None

    def start(self):
        '''Enter the guards, build the (clash-filtered) pose queue, make the sim visibly run,
        and hook the 'coord update' state machine.'''
        from contextlib import ExitStack
        isolde, sh, cfg = self.isolde, self.sh, self.cfg
        isolde.checkpoint()                             # so the user can revert
        self._was_paused = sh.pause
        self._stack = ExitStack()
        self._stack.enter_context(_fixed_temperature(isolde, sh, cfg['temperature']))
        self._stack.enter_context(_suspended_live_map_recalc(self.residue.structure))
        if cfg['map_decouple']:
            self._stack.enter_context(sh.map_decoupled(self.residue.atoms))
        self._prep = _prepare(sh, self.residue, cfg['moved_atoms'], cfg['accept_margin'])
        prep = self._prep
        # Work queue: current conformation first (exempt from the clash cull), then the
        # clash-filtered input poses.
        self._queue = [(CURRENT_LABEL, prep['base'][prep['ridx']].copy(), True)]
        for label, pose_coords in self.poses:
            if severe_clash(pose_coords, prep['ratoms'], prep['moved_mask'],
                            prep['env_coords']):
                self._n_culled += 1
                self.dlog('  culled "%s" (severe clash)' % label)
                continue
            self._queue.append((label, pose_coords, False))
        self.dlog('settle_poses (live) %s: %d pose(s) + current; margin %.2f kJ/mol'
                  % (self.residue, len(self.poses), prep['eff_margin']))
        if sh.pause:
            sh.pause = False                            # make it visibly run
        self._cu_handler = sh.triggers.add_handler('coord update', self._on_coord_update)
        self._start_pose(0)

    def _start_pose(self, i):
        '''Push pose ``i`` into the running sim (it settles over the next few frames) and
        soften the fragment against its surroundings so a clashy seed relaxes.'''
        self._pose_i = i
        label, pose_coords, is_current = self._queue[i]
        prep, sh = self._prep, self.sh
        coords = prep['base'].copy()
        coords[prep['ridx']] = pose_coords
        sh.push_coords_to_sim(coords)                   # applied on the next 'coord update'
        sh.soften_nb_selection(prep['moved_atoms'],
                               SETTLE_LAMBDA if self.cfg['soft_group'] else 1.0)
        self._cur_label = label
        self._cur_is_current = is_current
        self._frames = 0

    def _on_coord_update(self, trigger_name, data):
        '''State machine, one step per GUI update while the sim runs: let the current pose
        settle for FRAMES_PER_POSE updates (visible), snapshot its settled coords, then start
        the next pose. After the last pose, pause once and finalize.'''
        from chimerax.core.triggerset import DEREGISTER
        try:
            self._frames += 1
            if self._frames < self.FRAMES_PER_POSE:
                return                                  # keep settling this pose (visible)
            # Pose settled -> snapshot its coords. Reading state on 'coord update' is the safe
            # point (the worker has synced coordinates to the main thread).
            self._settled.append((self._cur_label, self.sh.sim_coords(),
                                  self._cur_is_current))
            self.dlog('  settled "%s" (live)' % self._cur_label)
            if self._pose_i + 1 < len(self._queue):
                self._start_pose(self._pose_i + 1)
                return
            self._cu_handler = None
            self._begin_finalize()                      # all poses done
            return DEREGISTER                           # stop the per-frame state machine
        except Exception:
            self.session.logger.report_exception()
            self._abort()
            return DEREGISTER

    def _begin_finalize(self):
        '''Pause the sim (safely -- deferring to 'sim paused' if it is still running, to avoid
        the overlapping-step-loop hazard) and then score + commit at fixed coords.'''
        sh = self.sh
        if sh.pause:
            self._finalize_now()
            return
        sh.pause = True
        def _deferred(trigger_name, data):
            from chimerax.core.triggerset import DEREGISTER
            try:
                self._finalize_now()
            except Exception:
                self.session.logger.report_exception()
                self._abort()
            return DEREGISTER                           # one-shot
        sh.triggers.add_handler('sim paused', _deferred)

    def _finalize_now(self):
        '''Paused: re-evaluate each settled snapshot's E_pack (the burst scoring), polish +
        commit, restore all guards + coupling, restore the pre-run pause state, deliver.'''
        sh, prep = self.sh, self._prep
        result = None
        try:
            results = []
            for label, coords, is_current in self._settled:
                sh.push_coords_to_sim(coords, immediate=True)   # fix these coords
                e = _pack_energy(sh, prep['moved_atoms'], prep['groups'])
                results.append((e, label, coords, is_current))
                self.dlog('  scored "%s": E_pack %.2f kJ/mol' % (label, e))
            result = _finalize_settle(
                self.session, self.isolde, self.residue, prep, results, self._n_culled,
                polish_top=self.cfg['polish_top'], minimize=self.cfg['minimize'],
                apply=self.cfg['apply'], map_decouple=self.cfg['map_decouple'],
                dlog=self.dlog)
        finally:
            _restore_coupling(self.session, sh, prep['moved_atoms'])
            if self._stack is not None:
                self._stack.close()                     # restore 0 K temp + live-map recalc
                self._stack = None
            if not self._was_paused:
                sh.pause = False                        # resume if the sim was running before
        self.session.logger.info(_summary_line(self.residue, result))
        if self.on_done is not None:
            try:
                self.on_done(result)
            except Exception:
                self.session.logger.report_exception()

    def _abort(self):
        '''Best-effort cleanup after a mid-run failure: restore coupling + guards + pause
        state so the simulation is not left decoupled or frozen at 0 K.'''
        try:
            _restore_coupling(self.session, self.sh, self._prep['moved_atoms'])
        except Exception:
            pass
        try:
            if self._stack is not None:
                self._stack.close()
                self._stack = None
        except Exception:
            pass
        try:
            if not self._was_paused and self.sh.sim_running:
                self.sh.pause = False
        except Exception:
            pass
        if self.on_done is not None:
            try:
                self.on_done(None)
            except Exception:
                pass


# ------------------------------------------------------------------
# Shared settle + rank + polish + commit (driver-agnostic core).
# ------------------------------------------------------------------
def _core_force_groups():
    from ..openmm.openmm_interface import CORE_FORCE_GROUPS
    return set(CORE_FORCE_GROUPS)


def _pack_energy(sh, atoms, groups):
    '''The fragment's OWN localised packing energy (kJ/mol) at the CURRENT (fixed) coords,
    by the coupled-minus-decoupled difference over the CORE force groups:

        E_pack = E_core(atoms fully coupled) - E_core(atoms decoupled to eps)

    Everything that does not depend on the fragment's coupling cancels, leaving just the
    fragment<->surroundings interaction. Lower is better; positive = net clash. Two
    single-point reads; leaves the coupling at eps (the caller re-establishes the settle
    coupling for the next pose, and the finally-block restores full coupling at the end).'''
    sh.soften_nb_selection(atoms, 1.0)                 # fragment fully coupled
    e_coupled = sh.potential_energy(groups)
    sh.soften_nb_selection(atoms, NB_DECOUPLE_EPS)     # fragment decoupled
    e_decoupled = sh.potential_energy(groups)
    return e_coupled - e_decoupled


def _settle_one(sh, base, ridx, pose_coords, moved_atoms, groups, settle_steps,
                soft_group, settle_lambda):
    '''Push one pose (full residue coords -> the construct's residue rows), settle it, and
    return (E_pack, settled_construct_coords). The search settle is dynamics-only (0 K); the
    honest minimised ranking is left to the polish. Softens the fragment during the settle
    when ``soft_group`` so a clashy seed relaxes instead of exploding.'''
    coords = base.copy()
    coords[ridx] = pose_coords
    sh.push_coords_to_sim(coords, immediate=True)
    sh.soften_nb_selection(moved_atoms, settle_lambda if soft_group else 1.0)
    sh.settle(settle_steps, minimize=False)
    settled = sh.sim_coords()
    return _pack_energy(sh, moved_atoms, groups), settled


def _polish_one(sh, base, ridx, settled_res_coords, moved_atoms, groups, steps,
                lam_from, lam_to, increments, minimize):
    '''Re-settle a candidate at FULL stiffness, RAMPING the fragment's coupling from
    ``lam_from`` up to ``lam_to`` over ``increments`` stages (a sudden jump can jolt atoms
    sitting in a soft overlap that abruptly becomes a hard wall). Starts from the ORIGINAL
    environment (base) + this pose's settled residue coords, so every polished candidate is
    judged against the same surroundings. Returns (E_pack, settled_construct_coords).'''
    coords = base.copy()
    coords[ridx] = settled_res_coords
    sh.push_coords_to_sim(coords, immediate=True)
    increments = max(1, int(increments))
    steps_per = max(1, int(round(steps / increments)))
    lambdas = np.linspace(lam_from, lam_to, increments + 1)[1:]   # end at lam_to (full)
    last = len(lambdas) - 1
    for i, lam in enumerate(lambdas):
        sh.soften_nb_selection(moved_atoms, float(lam))
        sh.settle(steps_per, minimize=(minimize and i == last))
    return _pack_energy(sh, moved_atoms, groups), sh.sim_coords()


def _prepare(sh, residue, moved_atoms, accept_margin):
    '''Common setup shared by both drivers: the residue's construct rows, the base coords,
    the fixed environment heavy atoms (for the clash filter), the moved-atom mask, and the
    do-no-harm margin (scaled by ISOLDE's MDFF coupling constant when a map drives the sim,
    else the absolute ``accept_margin``).'''
    ratoms = residue.atoms
    ridx = sh.atoms.indices(ratoms)
    base = sh.sim_coords()                              # the sim's current coords (Angstrom)
    heavy = sh.atoms.element_names != 'H'
    env_mask = heavy.copy()
    env_mask[ridx] = False
    moved_idx = ratoms.indices(moved_atoms)
    moved_mask = np.zeros(len(ratoms), dtype=bool)
    moved_mask[moved_idx[moved_idx >= 0]] = True
    gk = sh.mdff_coupling_constant()
    eff_margin = accept_margin * gk if gk > 0 else accept_margin
    return dict(ratoms=ratoms, groups=_core_force_groups(), ridx=ridx, base=base,
                env_coords=base[env_mask], moved_mask=moved_mask, moved_atoms=moved_atoms,
                gk=gk, eff_margin=eff_margin)


def _restore_coupling(session, sh, moved_atoms):
    '''Restore full coupling and release the fragment back to the environment group, so the
    simulation resumes exactly as a normal one. Run from every driver's cleanup path.'''
    try:
        sh.soften_nb_selection(moved_atoms, 1.0)
        sh.assign_nb_group(sh.atoms, 0)
    except Exception:
        session.logger.report_exception()


def _settle_search(sh, prep, poses, settle_steps, soft_group, dlog):
    '''Clash pre-filter + soft 0 K search: settle each pose (and the current conformation)
    and record its E_pack. Returns ``(results, n_culled)`` where results is
    ``[(E_pack, label, settled_construct_coords, is_current), ...]`` (unsorted). Synchronous
    (the burst driver); the live driver settles across frames instead but produces the same
    tuples for the shared tail.'''
    ridx, base, ratoms = prep['ridx'], prep['base'], prep['ratoms']
    moved_atoms, moved_mask = prep['moved_atoms'], prep['moved_mask']
    env_coords, groups = prep['env_coords'], prep['groups']
    # The current conformation is always a candidate (exempt from the clash cull) so an
    # already-good pose can win and be kept ("first, do no harm").
    all_poses = [(CURRENT_LABEL, base[ridx].copy())] + list(poses)
    results = []
    n_culled = 0
    for label, pose_coords in all_poses:
        is_current = (label == CURRENT_LABEL)
        if not is_current and severe_clash(pose_coords, ratoms, moved_mask, env_coords):
            n_culled += 1
            dlog('  culled "%s" (severe clash)' % label)
            continue
        e_pack, settled = _settle_one(sh, base, ridx, pose_coords, moved_atoms, groups,
                                      settle_steps, soft_group, SETTLE_LAMBDA)
        results.append((e_pack, label, settled, is_current))
        dlog('  settled "%s": E_pack %.2f kJ/mol' % (label, e_pack))
    return results, n_culled


def _finalize_settle(session, isolde, residue, prep, results, n_culled, *, polish_top,
                     minimize, apply, map_decouple, dlog):
    '''Shared tail (both drivers): polish the top-N survivors (+current) at full stiffness,
    re-rank, apply do-no-harm, and graft-commit ONLY ``residue.atoms`` onto the original
    environment. Assumes the sim is PAUSED (single-point reads + immediate coord push).
    Returns a SettleResult. Does NOT restore coupling -- the caller does that in cleanup.'''
    sh = isolde.sim_handler
    ridx, base = prep['ridx'], prep['base']
    moved_atoms, groups, eff_margin = (prep['moved_atoms'], prep['groups'],
                                       prep['eff_margin'])
    results = sorted(results, key=lambda t: t[0])
    # Polish the top survivors (plus current) at full stiffness and re-rank by the polished
    # E_pack: the soft-search order and the full-stiffness order can disagree.
    candidates = _polish(sh, base, ridx, results, moved_atoms, groups, polish_top,
                         minimize, dlog)
    best, kept_current = _choose_winner(candidates, eff_margin, dlog)
    energies = [(lbl, e) for e, lbl, _c, _cur in candidates]
    if best is None:
        return SettleResult(None, False, False, energies, n_culled, eff_margin,
                            map_decoupled=map_decouple)
    best_e, best_label, best_coords, _bc = best
    committed_res = None
    applied = False
    if apply:
        # Commit ONLY the target residue, grafted onto the ORIGINAL environment: the settled
        # poses are whole-construct snapshots whose surroundings were perturbed during
        # evaluation; move the residue's atoms and leave everything else exactly as it was.
        commit = base.copy()
        if not kept_current:
            commit[ridx] = best_coords[ridx]
        sh.push_coords_to_sim(commit, immediate=True)   # into the running (paused) sim
        sh.atoms.coords = commit                        # into the ChimeraX model
        committed_res = commit[ridx].copy()
        applied = True
        dlog('  -> commit "%s" (E_pack %.2f)%s'
             % (best_label, best_e, ' [kept current]' if kept_current else ''))
    else:
        sh.push_coords_to_sim(base, immediate=True)     # apply=False: restore, commit nothing
        dlog('  apply=False -> nothing committed')
    return SettleResult(
        committed=(best_label if applied else None), kept_current=kept_current,
        applied=applied, energies=energies, n_culled=n_culled,
        accept_margin_kJ=eff_margin, committed_coords=committed_res,
        map_decoupled=map_decouple)


def _settle_and_commit(session, isolde, residue, poses, *, moved_atoms, map_decouple,
                       soft_group, settle_steps, minimize, polish_top, accept_margin,
                       apply, debug, dlog):
    '''Burst path body (synchronous): prepare -> soft search -> polish + do-no-harm +
    commit. Assumes the sim is paused and the fit-scope guards are in force. Restores full
    coupling before returning. Returns a SettleResult.'''
    sh = isolde.sim_handler
    prep = _prepare(sh, residue, moved_atoms, accept_margin)
    dlog('settle_poses %s: %d pose(s) + current; MDFF global_k=%.3g -> margin %.2f kJ/mol'
         % (residue, len(poses), prep['gk'], prep['eff_margin']))
    try:
        results, n_culled = _settle_search(sh, prep, poses, settle_steps, soft_group, dlog)
        return _finalize_settle(session, isolde, residue, prep, results, n_culled,
                                polish_top=polish_top, minimize=minimize, apply=apply,
                                map_decouple=map_decouple, dlog=dlog)
    finally:
        _restore_coupling(session, sh, moved_atoms)


def _polish(sh, base, ridx, results, moved_atoms, groups, polish_top, minimize, dlog):
    '''Polish the top ``polish_top`` search survivors (plus the current conformation) at
    full stiffness and RE-RANK by the polished E_pack. Unpolished survivors keep their
    search E_pack. Returns the full candidate list best-first.'''
    if not results:
        return []
    n_polish = len(results) if polish_top <= 0 else min(polish_top, len(results))
    polish_set = list(results[:n_polish])
    # Always polish current too, so do-no-harm judges it at the same stiffness as its rivals.
    if not any(cur for _e, _l, _c, cur in polish_set):
        cur = next((t for t in results if t[3]), None)
        if cur is not None:
            polish_set.append(cur)
    polished_labels = set()
    out = []
    for _e, label, settled, is_current in polish_set:
        e_pol, coords = _polish_one(sh, base, ridx, settled[ridx], moved_atoms, groups,
                                    POLISH_STEPS, SETTLE_LAMBDA, 1.0, RAMP_INCREMENTS,
                                    minimize)
        out.append((e_pol, label, coords, is_current))
        polished_labels.add(label)
        dlog('  polished "%s": E_pack %.2f kJ/mol' % (label, e_pol))
    # Carry through any survivors we didn't polish, on their search score.
    for t in results:
        if t[1] not in polished_labels:
            out.append(t)
    out.sort(key=lambda t: t[0])
    return out


def _choose_winner(candidates, eff_margin, dlog):
    '''Pick the pose to commit from the ranked (best-first) candidates, applying do-no-harm:
    the current conformation is displaced only if the best challenger beats its E_pack by at
    least ``eff_margin`` kJ/mol. Returns (best_tuple_or_None, kept_current).'''
    if not candidates:
        return None, False
    best = candidates[0]
    curr = next((t for t in candidates if t[3]), None)
    if best[3]:                                        # the best pose IS current
        return best, True
    if curr is not None and (curr[0] - best[0]) < eff_margin:
        dlog('  do-no-harm: best beats current by %.2f kJ/mol (< %.2f) -> keeping current'
             % (curr[0] - best[0], eff_margin))
        return curr, True
    return best, False


def _summary_line(residue, result):
    '''One-line summary of a settle_poses run (the standard non-debug log line).'''
    if result is None or not result.energies:
        return 'settle_poses %s: no viable poses' % residue
    tail = '; sim left running ("isolde sim revert" to undo)' if result.applied else ''
    if not result.applied:
        return 'settle_poses %s: scored %d pose(s), nothing committed%s' % (
            residue, len(result.energies), tail)
    if result.kept_current:
        return 'settle_poses %s: kept current conformation (best fit already)%s' % (
            residue, tail)
    best_e = result.energies[0][1]
    return 'settle_poses %s: committed pose "%s" (E_pack %.2f kJ/mol)%s' % (
        residue, result.committed, best_e, tail)
