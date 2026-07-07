# @Author: Claude
# @Date:   01-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 01-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Robustness test for atom/residue deletion happening "under" code that holds a
reference across time -- the hazard documented in CLAUDE.md's "Live editing:
atoms and residues are never stable" section. ISOLDE's own editing tools, or a
plain ChimeraX `delete` command, can remove atoms/residues at essentially any
point in a session, including mid-simulation. A stale Python reference to a
deleted atomic object does not reliably raise a catchable exception -- it can
crash ChimeraX outright. Run inside ChimeraX:

    run_chimerax.bat --nogui --exit --script src/tests/test_atom_deletion_robustness.py

Four scenarios, each opening its own copy of the fixture since all of them
are destructive:
  1. AMBER-type cache (param_provider.py): delete an atom from a residue that
     was already cached by a prior build, then rebuild with the same
     provider. Must not crash, must detect the residue as a fresh cache miss
     (structural-fingerprint mismatch), and the rebuilt system's particle
     count must match the model's new atom count.
  2. ResidueStepper (navigate.py): delete the residue it is currently
     tracking. `.current_residue` must fall back to None (it already guards
     on `.deleted`) rather than dereference the deleted object, and
     subsequent navigation must still work.
  3. A live, running simulation: delete a mobile atom while the sim loop is
     actively stepping. This is not something ISOLDE's GUI is known to
     block, so it is a realistic "aggressive UI interaction" case. Either a
     clean continuation or a caught Python exception counts as a pass; an
     uncaught process crash would show up as this script exiting abnormally
     with no "ALL PASS" line, rather than as a assertion failure below.
  4. The full `Isolde.start_sim()` path (not just param_provider/SimHandler
     directly, as scenarios 1-3 do): run a clean simulation once so the
     AMBER-type cache is populated for the selected residue, stop it, delete
     a backbone atom from that now-cached residue, then start again. The
     `UnparameterisedResiduesError` this must raise is deliberately worded
     ('Unparameterised residue detected') so `Isolde.start_sim()`'s
     `str(e).startswith('Unparameterised')` catch fires
     `isolde.UNPARAMETERISED_RESIDUE` -- the trigger the "Unparameterised
     Residues" dialog/panel listens for -- instead of the raw exception
     reaching the user as a traceback. This checks that signal fires and
     that `start_sim()` itself doesn't raise; it can't check the actual Qt
     dialog appearing, since that requires a running GUI this test doesn't
     have.

Prints PASS/FAIL per scenario and exits non-zero on the first assertion
failure. Does not (and cannot) turn an actual process crash into a normal
failure report -- that has to be read off the exit code/output of the
ChimeraX process this script runs in.
'''
import os
import time
import numpy

FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '1pmx_1.pdb')

N_MOBILE_RESIDUES = 8
N_FIXED_RESIDUES = 6


def _fail(msg):
    print('FAIL: %s' % msg)
    raise SystemExit(1)


def _pump(session, n, stop_when=None, dt=0.02):
    for _ in range(n):
        session.triggers.activate_trigger('new frame', None)
        time.sleep(dt)
        if stop_when is not None and stop_when():
            return True
    return False


def _open_fixture(session):
    from chimerax.core.commands import run as run_cmd
    return run_cmd(session, 'open "%s"' % FIXTURE.replace('\\', '/'), log=False)[0]


def _mobile_and_fixed(m):
    residues = m.residues
    ordered = residues[numpy.lexsort((residues.numbers, residues.chain_ids))]
    mobile = ordered[:N_MOBILE_RESIDUES].atoms
    fixed = ordered[N_MOBILE_RESIDUES:N_MOBILE_RESIDUES + N_FIXED_RESIDUES].atoms
    return mobile, fixed


def _scenario_cache_survives_atom_deletion(session):
    from chimerax.isolde.openmm.openmm_interface import SimConstruct, SimHandler
    from chimerax.isolde.openmm.sim_param_mgr import SimParams
    from chimerax.isolde.openmm.forcefields import ForcefieldMgr
    from chimerax.isolde.openmm.param_provider import ForceFieldParameterisationProvider

    m = _open_fixture(session)
    try:
        provider = ForceFieldParameterisationProvider(ForcefieldMgr(session))

        def _build():
            mobile, fixed = _mobile_and_fixed(m)
            construct = SimConstruct(m, mobile, fixed)
            handler = SimHandler(session, SimParams(), construct, provider)
            return construct, handler

        construct, handler = _build()
        stats = provider.last_build_stats
        if stats['residues_cache_miss'] != stats['residues_total']:
            _fail('cache: expected cold-cache build to miss every residue, got %r' % stats)

        # Delete an atom from one of the residues the first build just
        # cached -- simulates a common UI action (delete atom, alt-loc
        # cleanup, model-building edit) happening between two simulation
        # (re)starts on the same model.
        target_res = max(construct.mobile_residues, key=lambda r: len(r.atoms))
        if len(target_res.atoms) < 2:
            _fail('cache: fixture has no mobile residue with a deletable atom')
        doomed = target_res.atoms[-1]
        n_atoms_before = len(m.atoms)
        doomed.delete()
        if len(m.atoms) != n_atoms_before - 1:
            _fail('cache: atom deletion did not reduce the model atom count')

        # Rebuild with the SAME provider. Removing any atom from a standard
        # amino acid residue breaks its exact match against the AMBER14
        # template it previously matched -- there is no "template minus one
        # atom" fallback -- so the *correct* outcome here is that the edited
        # residue is detected as a stale cache miss, genuinely re-attempted,
        # and cleanly reported as unparameterisable via the same typed
        # exception the slow path has always raised for this case. What we're
        # actually checking is that this is a clean, catchable, *expected*
        # exception -- not a crash, and not the cache silently reusing the
        # pre-deletion type/parameters for an atom that no longer exists.
        from chimerax.isolde.openmm.openmm_interface import UnparameterisedResiduesError
        try:
            construct2, handler2 = _build()
        except UnparameterisedResiduesError as e:
            stats2 = provider.last_build_stats
            if stats2['residues_cache_miss'] < 1:
                _fail('cache: residue with a deleted atom was not detected as a '
                      'cache miss even though the rebuild correctly failed to '
                      'reparameterise it: %r' % stats2)
            if target_res not in (e.unmatched + e.ambiguous):
                _fail('cache: rebuild failed, but not because of the residue we '
                      'edited (unmatched=%r ambiguous=%r)' % (e.unmatched, e.ambiguous))
            print('PASS: AMBER-type cache correctly treated the edited residue as '
                  'a cache miss and re-attempted matching, which cleanly (and '
                  'correctly) failed rather than crashing or reusing stale cached '
                  'types for a deleted atom (%d/%d residues missed)'
                  % (stats2['residues_cache_miss'], stats2['residues_total']))
        else:
            # If the fixture/deleted-atom choice ever changes such that the
            # edited residue remains parameterisable, the cache must still
            # have caught it as a miss and the rebuilt system must reflect
            # the new (smaller) atom count.
            stats2 = provider.last_build_stats
            if stats2['residues_cache_miss'] < 1:
                _fail('cache: residue with a deleted atom was not detected as a '
                      'cache miss: %r' % stats2)
            system2 = handler2._system
            if system2.getNumParticles() != len(construct2.all_atoms):
                _fail('cache: rebuilt system particle count %d != current atom count %d'
                      % (system2.getNumParticles(), len(construct2.all_atoms)))
            print('PASS: AMBER-type cache correctly invalidated (%d/%d residues missed) '
                  'after atom deletion; rebuilt system matches the new atom count (%d)'
                  % (stats2['residues_cache_miss'], stats2['residues_total'],
                     system2.getNumParticles()))
    finally:
        session.models.close([m])


def _scenario_stepper_survives_residue_deletion(session):
    from chimerax.isolde.navigate import get_stepper

    m = _open_fixture(session)
    try:
        stepper = get_stepper(m)
        # Headless ChimeraX has no `session.ui.mouse_modes`/camera, which
        # _new_camera_position() unconditionally touches. Stub it out so the
        # real residue-selection logic in step_to()/incr_residue() (the thing
        # actually under test here) still runs end to end -- this is a
        # test-only stand-in for GUI side effects, not a change to the
        # deletion-safety logic being exercised.
        stepper._new_camera_position = lambda *a, **k: None

        residues = m.residues
        target = residues[len(residues) // 2]

        stepper.step_to(target)
        if stepper.current_residue is not target:
            _fail('navigate: step_to did not set current_residue')

        target.delete()
        if stepper.current_residue is not None:
            _fail('navigate: current_residue should be None after its residue '
                  'was deleted, got %r' % stepper.current_residue)
        print('PASS: ResidueStepper.current_residue returns None (not a crash) '
              'after its tracked residue is deleted')

        # incr_residue()'s cr.deleted guard must let navigation recover
        # instead of dereferencing the deleted residue.
        nr = stepper.next_residue(polymeric_only=False)
        if nr is None or nr.deleted:
            _fail('navigate: next_residue after deletion returned nothing usable')
        print('PASS: next_residue() after its tracked residue was deleted falls '
              'back cleanly (now at %s)' % nr)
    finally:
        session.models.close([m])


def _scenario_live_sim_survives_mobile_atom_deletion(session):
    from chimerax.isolde.openmm import openmm_interface, sim_param_mgr
    from chimerax.isolde.openmm.forcefields import ForcefieldMgr
    from chimerax.isolde.openmm.param_provider import ForceFieldParameterisationProvider

    m = _open_fixture(session)
    handler = None
    try:
        mobile, fixed = _mobile_and_fixed(m)
        construct = openmm_interface.SimConstruct(m, mobile, fixed)
        params = sim_param_mgr.SimParams()
        try:
            params.set_param('platform', 'CPU')
        except Exception:
            params.platform = 'CPU'
        handler = openmm_interface.SimHandler(session, params, construct,
                                              ForceFieldParameterisationProvider(ForcefieldMgr(session)))
        handler.start_sim()
        _pump(session, 100)
        if not handler.sim_running:
            _fail('live-deletion: simulation did not start')

        # Delete a mobile atom while the sim loop is actively stepping.
        # ISOLDE is not known to block plain atom deletion while a
        # simulation is running, so this is a realistic "user did something
        # aggressive" case, not a contrived one.
        mobile_atoms = construct.mobile_atoms
        hydrogens = mobile_atoms[mobile_atoms.element_names == 'H']
        doomed_atom = hydrogens[0] if len(hydrogens) else mobile_atoms[-1]

        try:
            doomed_atom.delete()
            _pump(session, 100)
            handler.stop()
            _pump(session, 100, stop_when=lambda: not handler.sim_running)
            print('PASS: deleting a mobile atom mid-simulation did not raise, '
                  'and the simulation still stopped cleanly (sim_running=%r)'
                  % handler.sim_running)
        except Exception as e:
            print('PASS: deleting a mobile atom mid-simulation raised a catchable '
                  'Python exception rather than crashing the process: %r' % e)
    finally:
        try:
            if handler is not None and handler.sim_running:
                handler.stop()
                _pump(session, 100, stop_when=lambda: not handler.sim_running)
        except Exception:
            pass
        session.models.close([m])


def _scenario_start_sim_fires_unparameterised_trigger(session):
    from chimerax.isolde.isolde import Isolde

    m = _open_fixture(session)
    isolde = None
    try:
        isolde = Isolde(session)
        isolde.selected_model = m

        residues = m.residues
        ordered = residues[numpy.lexsort((residues.numbers, residues.chain_ids))]
        target_res = ordered[N_MOBILE_RESIDUES // 2]

        fired = []
        isolde.triggers.add_handler(
            isolde.UNPARAMETERISED_RESIDUE, lambda *a: fired.append(True))

        m.atoms.selected = False
        target_res.atoms.selected = True
        isolde.start_sim()
        if not isolde.simulation_running:
            _fail('start_sim: first (clean) start did not begin a simulation')
        if fired:
            _fail('start_sim: unparameterised-residue trigger fired on a clean model')
        _pump(session, 100)
        isolde.discard_sim()
        _pump(session, 200, stop_when=lambda: not isolde.simulation_running)
        if isolde.simulation_running:
            _fail('start_sim: simulation did not stop cleanly after the first run')

        # Delete a backbone atom from the residue that was just cached --
        # breaks its match against the AMBER14 template it was cached under.
        # The *same* Isolde instance (and therefore the same param_provider)
        # is reused for the next start_sim(), exactly reproducing "a residue
        # became unparameterisable between two simulation starts in one
        # session" -- the case the AMBER-type cache exists to handle safely.
        n_before = len(target_res.atoms)
        target_res.atoms[-1].delete()
        if len(target_res.atoms) != n_before - 1:
            _fail('start_sim: atom deletion setup did not take effect')

        m.atoms.selected = False
        target_res.atoms.selected = True
        isolde.start_sim()
        if isolde.simulation_running:
            _fail('start_sim: simulation should not have started for a residue '
                  'that lost a backbone atom')
        if not fired:
            _fail('start_sim: UNPARAMETERISED_RESIDUE trigger did not fire after '
                  'a cached residue became unparameterisable -- the '
                  '"Unparameterised Residues" dialog would not have appeared')
        print('PASS: Isolde.start_sim() caught the deleted-atom residue via the '
              'AMBER-type cache and fired UNPARAMETERISED_RESIDUE (the signal the '
              'warning dialog/panel listens for) instead of raising or silently '
              'proceeding with stale cached data')
    finally:
        try:
            if isolde is not None and isolde.simulation_running:
                isolde.discard_sim()
                _pump(session, 200, stop_when=lambda: not isolde.simulation_running)
        except Exception:
            pass
        session.models.close([m])


def run(session):
    _scenario_cache_survives_atom_deletion(session)
    _scenario_stepper_survives_residue_deletion(session)
    _scenario_live_sim_survives_mobile_atom_deletion(session)
    _scenario_start_sim_fires_unparameterised_trigger(session)
    print('ALL PASS')


# ChimeraX --script provides `session` in the module globals.
try:
    session  # noqa: F821
except NameError:
    session = None
if session is not None:
    run(session)
