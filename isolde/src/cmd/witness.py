# @Author: Tristan Croll
# @Date:   27-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 27-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Post-condition *witnesses* for mutating ISOLDE commands (design decision 5).

The single worst hazard for a non-seeing agent is the **silent no-op**: a
command returns success while nothing actually changed, and the agent proceeds
on false success. A witness turns "ok" into a *verifiable* statement of what
changed -- e.g. "heavy-atom RMSD moved 0.43 A" or "the integrator is advancing
with finite energy" -- so the agent (and our tests) can tell real work from a
no-op that merely didn't error.

Witnesses are keyed by an id on the command record (set via the agent manifest,
by category) and applied uniformly in :func:`chimerax.isolde.cmd.invoke.invoke_command`
-- they are never bolted onto individual command implementations.

Each witness is a ``(pre_fn, post_fn)`` pair:
    pre_fn(session, record, json_args)            -> opaque pre-state (or None)
    post_fn(session, record, json_args, pre_state) -> JSON-serializable dict
'''


def run_pre(session, record, json_args):
    pair = _WITNESSES.get(record.witness)
    if pair is None:
        return None
    return pair[0](session, record, json_args)


def run_post(session, record, json_args, pre):
    pair = _WITNESSES.get(record.witness)
    if pair is None:
        return None
    return pair[1](session, record, json_args, pre)


# --- helpers -----------------------------------------------------------------

def _isolde(session):
    return getattr(session, 'isolde', None)


def _selected_atoms(session):
    isolde = _isolde(session)
    if isolde is None:
        return None
    m = getattr(isolde, 'selected_model', None)
    if m is None:
        return None
    return getattr(m, 'atoms', None)


def _thread_handler(session):
    isolde = _isolde(session)
    if isolde is None:
        return None
    sh = getattr(isolde, 'sim_handler', None)
    return getattr(sh, 'thread_handler', None) if sh is not None else None


def _current_energy(session):
    '''Read the current potential energy SAFELY.

    OpenMM contexts are not thread-safe; the integrator runs on a worker thread.
    We read context state only when ``thread_finished()`` reports the worker idle
    (it is kicked exclusively from the UI thread's 'new frame', and this read is
    also on the UI thread, so an idle worker cannot restart mid-read). Otherwise
    we return None rather than risk a concurrent context access.
    '''
    th = _thread_handler(session)
    if th is None:
        return None
    try:
        if not th.thread_finished():
            return None
        from openmm import unit
        from ..openmm.openmm_interface import CORE_FORCE_GROUPS
        state = th.context.getState(getEnergy=True, groups=CORE_FORCE_GROUPS)
        return float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    except Exception:
        return None


def _sim_steps(session):
    '''Monotonic integrator-step count (or None). The cheap, authoritative
    "is the sim advancing" signal -- increasing across polls proves stepping.'''
    isolde = _isolde(session)
    sh = getattr(isolde, 'sim_handler', None) if isolde else None
    return getattr(sh, 'sim_steps', None) if sh is not None else None


# --- coord_move: did atoms actually move? ------------------------------------

def _coord_move_pre(session, record, json_args):
    atoms = _selected_atoms(session)
    if atoms is None:
        return None
    try:
        return {'n': len(atoms), 'coords': atoms.coords.copy()}
    except Exception:
        return None


def _coord_move_post(session, record, json_args, pre):
    atoms = _selected_atoms(session)
    if atoms is None:
        return {'moved': None, 'reason': 'no selected model'}
    out = {'n_atoms': len(atoms)}
    if pre is None or pre.get('coords') is None or pre['n'] != len(atoms):
        # atom count changed (e.g. add/remove) -- movement is self-evident.
        out['rmsd_moved'] = None
        out['atom_count_changed'] = pre is not None and pre.get('n') != len(atoms)
        return out
    try:
        import numpy
        before = pre['coords']
        after = atoms.coords
        d2 = ((after - before) ** 2).sum(axis=1)
        out['rmsd_moved'] = float(numpy.sqrt(d2.mean()))
        out['max_atom_shift'] = float(numpy.sqrt(d2.max())) if len(d2) else 0.0
        out['coord_checksum'] = int(numpy.frombuffer(
            after.tobytes(), dtype=numpy.uint8).sum())
        # Honesty about deferred effects: goal-directed commands (e.g. pepflip /
        # cisflip / restraint changes) applied during a running simulation set a
        # target the integrator realizes over subsequent frames -- so immediate
        # rmsd_moved can be ~0 even though real work was done. Tell the agent to
        # confirm via the coord_checksum changing across 'isolde sim status' polls.
        if bool(getattr(_isolde(session), 'simulation_running', False)):
            out['sim_running'] = True
            if out['rmsd_moved'] < 1e-4:
                out['note'] = ('no immediate movement, but a simulation is running: '
                               'this command may set a target realized over the next '
                               'frames -- poll "isolde sim status", watch sim_steps '
                               'advance, then re-check coordinates to confirm.')
    except Exception as e:
        out['rmsd_moved'] = None
        out['error'] = str(e)
    return out


# --- sim_control: is the simulation actually real and stepping? --------------

def _sim_control_pre(session, record, json_args):
    isolde = _isolde(session)
    return {
        'sim_running_before': bool(getattr(isolde, 'simulation_running', False)),
        'energy_before': _current_energy(session),
    }


def _sim_control_post(session, record, json_args, pre):
    isolde = _isolde(session)
    running = bool(getattr(isolde, 'simulation_running', False))
    sim_mgr = getattr(isolde, 'sim_manager', None)
    energy = _current_energy(session)
    cmd = json_args.get('cmd')
    out = {
        'cmd': cmd,
        'simulation_running': running,
        # The silent-no-op guard: a sim that builds a real OpenMM context has a
        # sim_manager. Energy populates once the worker has run a step (async via
        # 'new frame'), so it may be null in the same tick as start -- not a
        # no-op signal. The authoritative stepping proof is sim_steps increasing
        # across subsequent 'isolde sim status' polls.
        'sim_manager_present': sim_mgr is not None,
        'current_energy_kJ_mol': energy,
        'sim_steps': _sim_steps(session),
    }
    if cmd == 'start':
        out['started_ok'] = bool(running and sim_mgr is not None)
        if not out['started_ok']:
            out['warning'] = ('sim reports running=%s but sim_manager=%s '
                              '-- possible silent no-op' % (running, sim_mgr is not None))
    elif cmd == 'stop':
        # Stop tears the sim down ASYNCHRONOUSLY: the integrator halts and the
        # context is released on a subsequent 'new frame', so simulation_running
        # can still read True in this same tick -- that is NOT a failure (mirror
        # of start, where energy is null until the worker has stepped). Report
        # the request and the current (possibly stale) state, and tell the agent
        # to confirm by polling rather than trusting this immediate read.
        out['stop_requested'] = True
        out['stopped'] = not running
        if running:
            out['note'] = ('stop requested; teardown is asynchronous -- poll '
                           '"isolde sim status" until simulation_running is false.')
    return out


# --- restraint_change: did a restrain/adjust/release command change state? ---

# session_extensions getter name -> short label for the counts dict.
_RESTRAINT_MGR_GETTERS = (
    ('get_distance_restraint_mgr', 'distance'),
    ('get_proper_dihedral_restraint_mgr', 'proper_dihedral'),
    ('get_adaptive_dihedral_restraint_mgr', 'adaptive_dihedral'),
    ('get_position_restraint_mgr', 'position'),
    ('get_chiral_restraint_mgr', 'chiral'),
)


def _restraint_collection(mgr, residues):
    '''Return a restraint collection (with an `enableds` mask) for a manager,
    tolerant of which accessor it has: `all_restraints` (distance mgrs) or
    `get_all_restraints_for_residues` (dihedral/chiral/position mgrs).'''
    coll = getattr(mgr, 'all_restraints', None)
    if coll is not None:
        return coll
    getter = getattr(mgr, 'get_all_restraints_for_residues', None)
    if getter is not None and residues is not None:
        try:
            return getter(residues)
        except Exception:
            return None
    return None


def _count_one_mgr(mgr, residues):
    '''(total, enabled) for a restraint manager. Falls back to num_restraints
    (total only) if no enabled-bearing collection is available.'''
    coll = _restraint_collection(mgr, residues)
    if coll is not None:
        try:
            total = len(coll)
            enabled = int(coll.enableds.sum())
            return total, enabled
        except Exception:
            pass
    try:
        return int(mgr.num_restraints), None
    except Exception:
        return None, None


def _count_all_restraints(session):
    '''Per-manager {total, enabled} restraint counts for the selected model.
    Uses create=False so absent managers are simply skipped (count 0).'''
    isolde = _isolde(session)
    model = getattr(isolde, 'selected_model', None) if isolde else None
    if model is None:
        return None
    residues = getattr(model, 'residues', None)
    from .. import session_extensions as sx
    counts = {}
    for getter_name, label in _RESTRAINT_MGR_GETTERS:
        getter = getattr(sx, getter_name, None)
        if getter is None:
            continue
        try:
            mgr = getter(model, create=False)
        except Exception:
            mgr = None
        if mgr is None:
            counts[label] = {'total': 0, 'enabled': 0}
            continue
        total, enabled = _count_one_mgr(mgr, residues)
        counts[label] = {'total': total, 'enabled': enabled}
    # Adaptive distance restraints can live in several named groups.
    try:
        mgrs = sx.get_all_adaptive_distance_restraint_mgrs(model)
        tot = en = 0
        for mgr in mgrs:
            t, e = _count_one_mgr(mgr, residues)
            tot += (t or 0)
            en += (e or 0)
        counts['adaptive_distance'] = {'total': tot, 'enabled': en}
    except Exception:
        pass
    return counts


def _restraint_change_pre(session, record, json_args):
    return _count_all_restraints(session)


def _restraint_change_post(session, record, json_args, pre):
    post = _count_all_restraints(session)
    if pre is None or post is None:
        return {'restraint_counts': post,
                'note': 'no selected model; cannot diff restraint state'}
    total_deltas = {}
    enabled_deltas = {}
    total_delta = enabled_delta = 0
    for label, after in post.items():
        before = pre.get(label, {})
        dt = (after.get('total') or 0) - (before.get('total') or 0)
        if dt:
            total_deltas[label] = dt
            total_delta += dt
        be, ae = before.get('enabled'), after.get('enabled')
        if be is not None and ae is not None and ae != be:
            enabled_deltas[label] = ae - be
            enabled_delta += ae - be
    out = {'restraint_total_delta': total_delta,
           'restraint_enabled_delta': enabled_delta,
           'by_type_total_delta': total_deltas,
           'by_type_enabled_delta': enabled_deltas,
           'restraint_counts': post}
    if total_delta == 0 and enabled_delta == 0:
        out['note'] = ('no change in restraint count or enabled-count. For adjust* '
                       'commands that only retune parameters this is expected; '
                       'otherwise the command may have been a no-op.')
    elif total_delta == 0 and enabled_delta != 0:
        # release typically disables (enabled=False) rather than deleting.
        out['note'] = ('%d restraints toggled enabled-state (count unchanged) — '
                       'e.g. a release/disable.' % enabled_delta)
    return out


# --- convergence: delegated to the async job layer (jobs.py) -----------------

def _convergence_pre(session, record, json_args):
    return None


def _convergence_post(session, record, json_args, pre):
    # brefine/brsr launch their own worker thread and return immediately, so the
    # command's completion does NOT mean refinement has converged. The agent
    # tracks progress by polling map_info (Rwork/Rfree for crystallographic) and
    # bfactor_outliers, not by this witness.
    return {'launched': True,
            'note': ('refinement runs asynchronously in the background; poll '
                     'map_info (Rwork/Rfree) and bfactor_outliers to track '
                     'progress/convergence rather than relying on job completion.')}


_WITNESSES = {
    'coord_move':       (_coord_move_pre, _coord_move_post),
    'sim_control':      (_sim_control_pre, _sim_control_post),
    'restraint_change': (_restraint_change_pre, _restraint_change_post),
    'convergence':      (_convergence_pre, _convergence_post),
}
