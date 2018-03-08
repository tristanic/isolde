# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.



import numpy
import os, sys
import multiprocessing as mp
import ctypes
from math import pi, radians, degrees, cos
from warnings import warn
from time import time, sleep
from simtk import unit
from simtk.unit import Quantity, Unit
from chimerax.core.atomic import concatenate, Bonds

from ..checkpoint import CheckPoint
from ..threading.shared_array import TypedMPArray, SharedNumpyArray
from . import sim_thread
from .sim_thread import SimComms, SimThread, ChangeTracker
from ..constants import defaults, sim_outcomes, control
from ..param_mgr import Param_Mgr, autodoc, param_properties

OPENMM_LENGTH_UNIT = defaults.OPENMM_LENGTH_UNIT
OPENMM_FORCE_UNIT = defaults.OPENMM_FORCE_UNIT
OPENMM_SPRING_UNIT = defaults.OPENMM_SPRING_UNIT
OPENMM_RADIAL_SPRING_UNIT = defaults.OPENMM_RADIAL_SPRING_UNIT
OPENMM_ENERGY_UNIT = defaults.OPENMM_ENERGY_UNIT
OPENMM_ANGLE_UNIT = defaults.OPENMM_ANGLE_UNIT
OPENMM_TIME_UNIT = defaults.OPENMM_TIME_UNIT
OPENMM_DIPOLE_UNIT = defaults.OPENMM_DIPOLE_UNIT
OPENMM_TEMPERATURE_UNIT = defaults.OPENMM_TEMPERATURE_UNIT
CHIMERAX_LENGTH_UNIT        = defaults.CHIMERAX_LENGTH_UNIT
CHIMERAX_FORCE_UNIT         = defaults.CHIMERAX_FORCE_UNIT
CHIMERAX_SPRING_UNIT        = defaults.CHIMERAX_SPRING_UNIT

SIM_MODE_MIN                = control.SIM_MODE_MIN
SIM_MODE_EQUIL              = control.SIM_MODE_EQUIL
SIM_MODE_UNSTABLE           = control.SIM_MODE_UNSTABLE

CARBON_MASS = defaults.CARBON_MASS
CARBON_ATOMIC_NUMBER = defaults.CARBON_ATOMIC_NUMBER

def error_cb(e):
    print(e.__traceback__)
    print(e)

FLOAT_TYPE = defaults.FLOAT_TYPE

_amber14 = ['amberff14SB.xml','tip3p_standard.xml',
            'tip3p_HFE_multivalent.xml', 'tip3p_IOD_multivalent.xml']

cwd = os.path.dirname(os.path.abspath(__file__))
amber14 = [os.path.join(cwd,'amberff',f) for f in _amber14]

@param_properties
@autodoc
class SimParams(Param_Mgr):
    '''
    Container for all the parameters needed to initialise a simulation
    '''
    _default_params = {
        'restraint_max_force':                  (defaults.MAX_RESTRAINT_FORCE, OPENMM_FORCE_UNIT),
        'distance_restraint_spring_constant':   (defaults.DISTANCE_RESTRAINT_SPRING_CONSTANT, OPENMM_SPRING_UNIT),
        'position_restraint_spring_constant':   (defaults.POSITION_RESTRAINT_SPRING_CONSTANT, OPENMM_SPRING_UNIT),
        'haptic_spring_constant':               (defaults.HAPTIC_SPRING_CONSTANT, OPENMM_SPRING_UNIT),
        'mouse_tug_spring_constant':            (defaults.MOUSE_TUG_SPRING_CONSTANT, OPENMM_SPRING_UNIT),
        'tug_max_force':                        (defaults.MAX_TUG_FORCE, OPENMM_FORCE_UNIT),
        'dihedral_restraint_cutoff_angle':      (defaults.DIHEDRAL_RESTRAINT_CUTOFF, OPENMM_ANGLE_UNIT),
        'rotamer_restraint_cutoff_angle':       (defaults.ROTAMER_RESTRAINT_CUTOFF, OPENMM_ANGLE_UNIT),
        'rotamer_spring_constant':              (defaults.ROTAMER_SPRING_CONSTANT, OPENMM_RADIAL_SPRING_UNIT),
        'peptide_bond_spring_constant':         (defaults.PEPTIDE_SPRING_CONSTANT, OPENMM_RADIAL_SPRING_UNIT),
        'phi_psi_spring_constant':              (defaults.PHI_PSI_SPRING_CONSTANT, OPENMM_RADIAL_SPRING_UNIT),
        'cis_peptide_bond_cutoff_angle':        (defaults.CIS_PEPTIDE_BOND_CUTOFF, OPENMM_ANGLE_UNIT),
        'max_atom_movement_per_step':           (defaults.MAX_ATOM_MOVEMENT_PER_STEP, OPENMM_LENGTH_UNIT),
        'max_allowable_force':                  (defaults.MAX_ALLOWABLE_FORCE, OPENMM_FORCE_UNIT),
        'max_stable_force':                     (defaults.MAX_STABLE_FORCE, OPENMM_FORCE_UNIT),
        'friction_coefficient':                 (defaults.OPENMM_FRICTION, 1/OPENMM_TIME_UNIT),
        'temperature':                          (defaults.TEMPERATURE, OPENMM_TEMPERATURE_UNIT),

        'nonbonded_cutoff_method':              (defaults.OPENMM_NONBONDED_METHOD, None),
        'nonbonded_cutoff_distance':            (defaults.OPENMM_NONBONDED_CUTOFF, OPENMM_LENGTH_UNIT),

        'vacuum_dielectric_correction':         (defaults.VACUUM_DIELECTRIC_CORR, OPENMM_DIPOLE_UNIT),

        'use_gbsa':                             (defaults.USE_GBSA, None),
        'gbsa_cutoff_method':                   (defaults.GBSA_NONBONDED_METHOD, None),
        'gbsa_solvent_dielectric':              (defaults.GBSA_SOLVENT_DIELECTRIC, OPENMM_DIPOLE_UNIT),
        'gbsa_solute_dielectric':               (defaults.GBSA_SOLUTE_DIELECTRIC, OPENMM_DIPOLE_UNIT),
        'gbsa_sa_method':                       (defaults.GBSA_SA_METHOD, None),
        'gbsa_cutoff':                          (defaults.GBSA_CUTOFF, OPENMM_LENGTH_UNIT),
        'gbsa_kappa':                           (defaults.GBSA_KAPPA, 1/OPENMM_LENGTH_UNIT),

        'rigid_bonds':                          (defaults.RIGID_BONDS, None),
        'rigid_water':                          (defaults.RIGID_WATER, None),
        'remove_c_of_m_motion':                 (defaults.REMOVE_C_OF_M_MOTION, None),

        'platform':                             (defaults.OPENMM_PLATFORM, None),
        'forcefield':                           ('amber14', None),
        'integrator':                           (defaults.OPENMM_INTEGRATOR_TYPE, None),
        'variable_integrator_tolerance':        (defaults.OPENMM_VAR_INTEGRATOR_TOL, None),
        'fixed_integrator_timestep':            (defaults.OPENMM_FIXED_INTEGRATOR_TS, None),
        'constraint_tolerance':                 (defaults.OPENMM_CONSTRAINT_TOL, None),
        'sim_steps_per_gui_update':             (defaults.SIM_STEPS_PER_GUI_UPDATE, None),
        'minimization_steps_per_gui_update':    (defaults.MIN_STEPS_PER_GUI_UPDATE, None),
        'simulation_startup_rounds':            (defaults.SIM_STARTUP_ROUNDS, None),
        'maximum_unstable_rounds':              (defaults.MAX_UNSTABLE_ROUNDS, None),
        'minimization_convergence_tol':         (defaults.MIN_CONVERGENCE_FORCE_TOL, OPENMM_FORCE_UNIT),
        'tug_hydrogens':                        (False, None),
        'hydrogens_feel_maps':                  (defaults.HYDROGENS_FEEL_MAPS, None),
        'target_loop_period':                   (defaults.TARGET_LOOP_PERIOD, None),

    }





def start_pool(sim_params, sim_data, comms_object, change_tracker, use_mp = False):
    #~ try:
        #~ from chimerax import app_bin_dir
    #~ except:
        #~ return
    if use_mp:
        from multiprocessing.pool import Pool
    else:
        from multiprocessing.pool import ThreadPool as Pool
    thread_pool = Pool(processes=1, initializer=sim_thread._init_sim_thread,
        initargs=(sim_params, sim_data, comms_object, change_tracker))
    return thread_pool

class ChimeraXSimInterface:
    '''
    Application-facing interface between ChimeraX and a SimThread object.
    Handles interconversion between ChimeraX and generic Numpy types and
    bi-directional communication. A new instance is created for each
    simulation, and should be discarded once the simulation is done.
    '''
    TIMEOUT = defaults.SIM_TIMEOUT # seconds

    def __init__(self, session, isolde, time_loops = False):
        self.session = session
        self.isolde = isolde
        self._sim_thread = None
        self._pool = None
        self.time_loops = time_loops
        self._last_time = None

            # Registry for all the event handlers that should only be
            # active when a simulation is running. NOTE: these will
            # automatically be destroyed on simulation termination.
        self._sim_event_names = {
            'isolde': [],
            'session': []
            }
        self.isolde._isolde_events.add_event_handler('sim_interface_sim_start',
            'simulation started', self._sim_start_cb)
        self.isolde._isolde_events.add_event_handler('cleanup_on_sim_termination',
            'simulation terminated', self.finalize)

    def finalize(self, *_):
        #self.stop_sim()
        self._release_all_sim_events()
        self._disable_mouse_tugging()
        self.isolde._isolde_events.remove_event_handler('sim_interface_sim_start')
        self.isolde._isolde_events.remove_event_handler('cleanup_on_sim_termination')

    def _sim_start_cb(self, *_):
        self._initialize_mouse_tugging()
        if self.isolde.params.track_ramachandran_status:
            self._register_sim_event('ramachandran update', 'isolde',
                                     'completed simulation step',
                                      self.isolde.update_ramachandran)
        if self.isolde.params.track_rotamer_status:
            self._register_sim_event('omega update', 'isolde',
                                     'completed simulation step',
                                     self.isolde.update_omega_check)



    def _initialize_mouse_tugging(self):
        from .. import mousemodes
        isolde = self.isolde
        mt = self._mouse_tugger = mousemodes.TugAtomsMode(
            self.session, self.tuggable_atoms, isolde._annotations)
        isolde._mouse_modes.register_mode(mt.name, mt, 'right', ('alt',))
        self._register_sim_event('mouse tugging', 'session', 'new frame',
                                 self._update_mouse_tugging)

    def _disable_mouse_tugging(self):
        self.isolde._mouse_modes.remove_mode(self._mouse_tugger.name)

    def _update_mouse_tugging(self, *_):
        mtug = self._mouse_tugger
        cur_tug = mtug.already_tugging
        tugging, tug_atom, xyz0 = mtug.status
        if tugging:
            if not cur_tug:
                mtug.last_tugged_atom = tug_atom
                mtug.already_tugging = True
            self.tug_atom_to(tug_atom, xyz0)
        elif cur_tug:
            self.release_tugged_atom(mtug.last_tugged_atom)
            mtug.last_tugged_atom = None
            mtug.already_tugging = False



    @property
    def sim_running(self):
        return self._pool is not None

    def toggle_pause(self, acknowledge = False, timeout = defaults.COMMS_TIMEOUT):
        '''
        Toggle between pause and resume, optionally waiting for an
        acknowledgement that the simulation is actually paused before
        continuing
        '''
        if not self.sim_running:
            return
        comms = self.comms_object
        ct = self.change_tracker
        status_q = comms['status']
        currently_paused = comms['pause'].value
        comms.thread_safe_set_value('pause', not currently_paused)
        ct.register_change(ct.PAUSE_TOGGLE)
        if acknowledge:
            start_time = time()
            elapsed_time = 0
            message_received = False
            while elapsed_time < timeout:
                if not status_q.empty():
                    break
                else:
                    sleep(1e-2)
                    elapsed_time = time() - start_time
            if elapsed_time < timeout:
                if not status_q.empty():
                    m = status_q.get()
                else:
                    raise Exception('This should not happen')
                if not ((currently_paused and m == 'Resumed') or \
                        (not currently_paused and m == 'Paused')):
                    raise RuntimeError('Unexpected message: "{}" received '\
                        +'while attempting to toggle pause state!'.format(m))
            else:
                raise RuntimeError('Timed out waiting for acknowledgement!')
            print(m)
        if currently_paused:
            self.isolde.triggers.activate_trigger('simulation resumed', None)
        else:
            self.isolde.triggers.activate_trigger('simulation paused', None)

    @property
    def paused(self):
        return self.comms_object['pause'].value

    @paused.setter
    def paused(self, flag):
        comms = self.comms_object
        ct = self.change_tracker
        if flag != comms['pause'].value:
            self.toggle_pause()

    def pause_and_acknowledge(self, timeout=defaults.COMMS_TIMEOUT):
        if self.paused:
            return
        self.toggle_pause(acknowledge=True, timeout=timeout)

    def resume_and_acknowledge(self, timeout=defaults.COMMS_TIMEOUT):
        if not self.paused:
            return
        self.toggle_pause(acknowledge=True, timeout=timeout)


    def stop_sim(self, reason, err = None):
        '''
        Try to stop the simulation gracefully, or halt the simulation
        thread on a timeout.
        '''
        if not self.sim_running:
            return
        self.change_tracker.register_change(self.change_tracker.STOP)
        self._pool.terminate()
        self._pool = None
        if reason == sim_outcomes.DISCARD:
            self.starting_checkpoint.revert()
        self.isolde.triggers.activate_trigger('simulation terminated', (reason, err))
        # TODO: Add handler for timeout situations

    def _set_temperature(self, temperature):
        ct = self.change_tracker
        self.comms_object.thread_safe_set_value('temperature', temperature)
        ct.register_change(ct.TEMPERATURE)

    def _get_temperature(self):
        return self.comms_object['temperature'].value

    temperature = property(_get_temperature, _set_temperature)

    @property
    def new_coords_ready(self):
        return self.change_tracker.changes & self.change_tracker.COORDS_READY

    @property
    def sim_failed(self):
        return self.change_tracker.changes & self.change_tracker.ERROR

    @property
    def sim_unstable(self):
        return self.change_tracker.changes & self.change_tracker.UNSTABLE

    @property
    def sim_mode(self):
        mode = self.comms_object['sim mode'].value
        if mode == SIM_MODE_MIN:
            return 'min'
        elif mode == SIM_MODE_EQUIL:
            return 'equil'
        elif mode == SIM_MODE_UNSTABLE:
            return 'unstable!'
        else:
            raise RuntimeError('Simulation is in unrecognised state!')

    @sim_mode.setter
    def sim_mode(self, mode):
        comms = self.comms_object
        if mode not in ('min', 'equil'):
            raise RuntimeError('Simulation mode should be either "min" or "equil"!')
        if mode == 'min':
            comms.thread_safe_set_value('sim mode', SIM_MODE_MIN)
        elif mode == 'equil':
            comms.thread_safe_set_value('sim mode', SIM_MODE_EQUIL)
        self.change_tracker.register_change(self.change_tracker.MODE)


    def tug_atom_to(self, atom, target, spring_constant = None):
        '''
        Tug an atom towards the target (x,y,z) position.
        @ param index:
            The index of the atom in the array of tuggable atoms.
        @param target:
            The (x,y,z) coordinates to tug the atom towards, either in
            Angstroms or as a simtk Quantity.
        @param spring_constant:
            The spring constant to tug with, either in kJ mol-1 Angstrom-2
            or as a simtk Quantity. If no spring constant is given, the
            default mouse tugging spring constant will be used.
        '''
        #index = self.all_atoms.index(atom)
        params = self.sim_params
        comms = self.comms_object
        index = self.tuggable_atoms.index(atom)
        ct = self.change_tracker
        if spring_constant is None:
            spring_constant = params['mouse_tug_spring_constant']
        if type(spring_constant) == Quantity:
            spring_constant = spring_constant.value_in_unit(OPENMM_SPRING_UNIT)
        else:
            spring_constant = \
                spring_constant * CHIMERAX_SPRING_UNIT/OPENMM_SPRING_UNIT
        if type(target) == Quantity:
            target = target.value_in_unit(OPENMM_LENGTH_UNIT)
        else:
            target = target * CHIMERAX_LENGTH_UNIT/OPENMM_LENGTH_UNIT
        target_key = 'tugging targets'
        k_key = 'tugging spring constants'
        target_array = comms[target_key]
        k_array = comms[k_key]
        with target_array.get_lock(), k_array.get_lock():
            target_array[index] = target
            k_array[index] = spring_constant
        ct.register_array_changes('tugging', indices = index)

    def release_tugged_atom(self, atom):
        zeros = numpy.array([0,0,0], numpy.double)
        self.tug_atom_to(atom, zeros, spring_constant = 0)

    def update_dihedral_restraints(self, dihedrals):
        ct = self.change_tracker
        name = dihedrals[0].name
        master_array = self.named_dihedrals[name]
        indices = master_array.indices(dihedrals)
        sim_targets, sim_ks = ct.get_managed_arrays(name)
        with sim_targets.get_lock(), sim_ks.get_lock():
            sim_targets[indices] = dihedrals.targets
            sim_ks[indices] = dihedrals.spring_constants
        ct.register_array_changes(name, indices = indices)



    def update_dihedral_restraint(self, dihedral):
        ct = self.change_tracker
        name = dihedral.name
        master_array = self.named_dihedrals[name]
        index = master_array.index(dihedral)
        sim_targets, sim_ks = ct.get_managed_arrays(name)
        with sim_targets.get_lock(), sim_ks.get_lock():
            sim_targets[index] = dihedral.target
            sim_ks[index] = dihedral.spring_constant
        ct.register_array_changes(name, indices = index)

    def update_rotamer_target(self, rotamer):
        key = 'rotamers'
        comms = self.comms_object
        ct = self.change_tracker
        r_index = self.mobile_residues.index(rotamer.residue)
        master_dict = comms[key]
        restrained_mask = ct.get_managed_arrays(key)[0]
        target_array = master_dict['targets'][r_index]
        k_array = master_dict['spring constants'][r_index]
        with restrained_mask.get_lock(), target_array.get_lock(), k_array.get_lock():
            restrained_mask[r_index] = rotamer.restrained
            target_array[:] = rotamer.target
            k_array[:] = rotamer.spring_constants
        ct.register_array_changes(key, indices = r_index)

    def update_position_restraint(self, position_restraint):
        '''
        Restrain an atom to an (x,y,z) position with the given spring constant.
        @param position_restraint:
            A ChimeraX Position_Restraint object.
        '''
        comms = self.comms_object
        ct = self.change_tracker
        (targets, ks) = ct.get_managed_arrays('position restraints')
        pr = position_restraint
        atom = pr.atom
        index = self._restrainable_atoms.index(atom)
        target = pr.target
        k = pr.spring_constant
        with targets.get_lock(), ks.get_lock():
            targets[index] = target
            ks[index] = k
        ct.register_array_changes('position restraints', indices = index)

    def update_position_restraints(self, position_restraints):
        '''
        Update a set of position restraints in the simulation
        '''
        comms = self.comms_object
        ct = self.change_tracker

        pr = position_restraints
        atoms = pr.atoms
        indices = self._restrainable_atoms.indices(atoms)
        targets = pr.targets
        ks = pr.spring_constants
        sim_targets, sim_ks = ct.get_managed_arrays('position restraints')
        with sim_targets.get_lock(), sim_ks.get_lock():
            sim_targets[indices] = targets
            sim_ks[indices] = ks
        ct.register_array_changes('position restraints', indices = indices)


    def update_distance_restraints(self, distance_restraints):
        '''
        Update a set of distance restraints in the simulation at once.
        Args:
            distance_restraints:
                A Distance_Restraints object
        '''
        name = distance_restraints.name
        master_list = self.distance_restraints_dict[name]
        indices = master_list.indices(distance_restraints)
        ct = self.change_tracker
        sim_targets, sim_ks = ct.get_managed_arrays(name)
        with sim_targets.get_lock(), sim_ks.get_lock():
            sim_targets[indices] = distance_restraints.targets
            sim_ks[indices] = distance_restraints.spring_constants
        ct.register_array_changes(name, indices = indices)

    def update_distance_restraint(self, distance_restraint):
        '''
        Update a distance restraint in the simulation.
        '''
        name = distance_restraint.name
        master_list = self.distance_restraints_dict[name]
        index = master_list.index(distance_restraint)
        ct = self.change_tracker
        sim_targets, sim_ks = ct.get_managed_arrays(name)
        with sim_targets.get_lock(), sim_ks.get_lock():
            sim_targets[index] = distance_restraint.target_distance
            sim_ks[index] = distance_restraint.spring_constant
        ct.register_array_changes(name, indices = index)

    def update_coords(self, atoms):
        '''
        Push new atomic coordinates to the simulation for the given
        atoms. Use with care! This will trigger a minimisation to
        deal with the almost-inevitable clashes, but minimisation can't
        work miracles.
        '''
        ct = self.change_tracker
        sim_coords = ct.get_managed_arrays('coords')[0]
        indices = self.all_atoms.indices(atoms)
        with sim_coords.get_lock():
            sim_coords[indices] = atoms.coords
        ct.register_array_changes('coords', indices=indices)

    def _register_sim_event(self, name, owner, trigger_name, callback):
        if owner == 'isolde':
            registry = self.isolde._isolde_events
        elif owner == 'session':
            registry = self.isolde._event_handler
        registry.add_event_handler(name, trigger_name, callback)
        self._sim_event_names[owner].append(name)

    def _release_sim_event(self, name, owner):
        name_list = self._sim_event_names[owner]
        name_index = name_list.index(name)
        if owner == 'isolde':
            registry = self.isolde._isolde_events
        elif owner == 'session':
            registry = self.isolde._event_handler
        registry.remove_event_handler(name)
        name_list.pop(name_index)

    def _release_all_sim_events(self, *_):
        for owner, name_list in self._sim_event_names.items():
            for name in reversed(name_list):
                self._release_sim_event(name, owner)


    def _sim_init_handler(self, *_):
        '''
        Handler to run on ChimeraX 'new frame' trigger while the simulation
        is initialising. On successful initialisation, hands things over
        to self._sim_loop_handler().
        '''
        if time() - self._init_start_time > self.TIMEOUT:
            self._pool.terminate()
            self._pool = None
            err = TimeoutError('Timed out after {:0.1f} seconds waiting for '\
                               'the simulation thread to start. Unless your '\
                               'simulation construct is extremely large, this '\
                               'indicates a serious error.'.format(
                               self.TIMEOUT))
            self.isolde.triggers.activate_trigger('simulation terminated', (sim_outcomes.TIMEOUT, err))

        comms = self.comms_object
        ct = self.change_tracker
        isolde = self.isolde
        err_q = comms['error']
        status_q = comms['status']
        while not status_q.empty():
            print(status_q.get())

        with ct.changes.get_lock():
            changes = ct.changes.value
            ct.clear_outputs()

        if changes & ct.ERROR:
            err, traceback = err_q.get()
            print(traceback)
            self.stop_sim(sim_outcomes.ERROR, err)
            #TODO: provide the user with a dialog to choose whether to
            #keep or discard the coordinates in this case.

        if changes & ct.INIT_COMPLETE:
            self.session.logger.status('Initialisation complete')
            self._release_sim_event('simulation_initialization', 'session')
            # Tell the thread it's ok to start
            comms['status'].put('Go')
            self.isolde.triggers.activate_trigger('simulation started', None)
            self._register_sim_event('check_sim_on_gui_update',
                    'session', 'new frame', self._sim_loop_handler)



    def _sim_loop_handler(self, *_):
        '''
        Primary handler for the simulation loop. On each ChimeraX
        'new frame' trigger, it checks to see if the simulation has
        finished an iteration, and if so updates the coordinates in
        ChimeraX, gets any other simulation information, and triggers
        the next cycle.
        '''
        session = self.session
        isolde = self.isolde
        comms = self.comms_object
        ct = self.change_tracker
        thread = self._sim_thread
        err_q = comms['error']
        status_q = comms['status']
        self.message_received = False
        self.messages = []
        while not status_q.empty():
            self.message_received = True
            m = status_q.get()
            self.messages.append(m)
            print(m)
        with ct.changes.get_lock():
            changes = ct.changes.value
            ct.clear_outputs()

        if self.step_counter == 10:
            self.sim_mode = 'equil'

        elif self.step_counter < 10 and (changes & ct.MIN_COMPLETE):
            self.sim_mode = 'equil'

        if changes & ct.ERROR:
            if changes & ct.UNSTABLE:
                self.stop_sim(sim_outcomes.UNSTABLE)
                return

            err_q = comms['error']
            if not err_q.empty():
                err, traceback = err_q.get()
            else:
                err_str = '''
                Simulation thread returned an error code but not an
                exception. This should not happen.
                '''
                err = RuntimeError(err_str)
                traceback = ''
            print(traceback)
            self.stop_sim(sim_outcomes.DISCARD, err)


        if changes & ct.UNSTABLE:
            warn_str = 'Simulation is unstable! Attempting to minimise excessive force...'
            session.logger.info(warn_str)

        if changes & ct.COORDS_READY:
            if self.time_loops:
                print('coords ready!')
                if self._last_time is None:
                    self._last_time = time()
                self.loop_time = time()-self._last_time
                self._last_time = time()
            new_coords = comms['coords']
            with new_coords.get_lock():
                self.all_atoms.coords = new_coords
            self.step_counter += 1
            isolde.triggers.activate_trigger('completed simulation step', self.step_counter)


    def start_sim_thread(self, sim_params, all_atoms, tuggable_atoms, fixed_flags,
                         backbone_dihedrals, rotamers, distance_restraints,
                         position_restraints, density_maps = None):
        '''
        Start an OpenMM simulation in a separate Python instance, with
        communication via shared variables. Returns a dict containing all
        variables to be used for communication/control, with descriptive
        names.
            @param  sim_params:
                An instance of SimParams.params defining the basic simulation
                settings
            @param all_atoms:
                The set of all ChimeraX atoms defining the complete simulation
                construct
            @param fixed_flags:
                A numpy boolean array flagging atoms to be fixed in space.
            @param backbone_dihedrals:
                An ISOLDE Backbone_Dihedrals object defining all phi, psi and
                omega dihedrals in the simulation
            @param rotamers:
                An ISOLDE Rotamers object
            @param distance_restraints:
                A {name: object} dict of Distance_Restraints objects.
            @param position_restraints
                An ISOLDE Position_Restraints object
            @param maps
                A dict of ISOLDE map objects
        '''
        if self.sim_running:
            raise RuntimeError('You already have a simulation running!')

        #~ manager = self.manager = mp.Manager()
        #ct = self.change_tracker = ChangeTracker()

        # Container for all data that can be changed by ISOLDE, the simulation
        # thread, or both
        comms = self.comms_object = SimComms()
        ct = self.change_tracker = comms['changes']
        self.distance_restraints_dict = distance_restraints
        self.sim_params = sim_params
        self.temperature = sim_params['temperature'].value_in_unit(OPENMM_TEMPERATURE_UNIT)

        self.tuggable_atoms = tuggable_atoms
        fixed_indices = self.fixed_indices = numpy.where(fixed_flags)[0]
        mobile_indices = self.mobile_indices = numpy.where(numpy.invert(fixed_flags))[0]
        self.all_atoms = all_atoms
        fixed_atoms = self.fixed_atoms = all_atoms[fixed_indices]
        mobile_atoms = self.mobile_atoms = all_atoms[mobile_indices]
        self.mobile_residues = mobile_atoms.unique_residues
        self.mobile_heavy_atoms = mobile_atoms[mobile_atoms.element_numbers != 1]
        residues = all_atoms.residues
        atom_names    = all_atoms.names
        element_names = all_atoms.element_names
        residue_names = residues.names
        residue_nums  = residues.numbers
        chain_ids     = residues.chain_ids
        self.starting_coords = all_atoms.coords
        coords        = self.starting_coords * unit.angstrom

        bonds               = all_atoms.intra_bonds
        bond_atoms          = bonds.atoms
        bond_atom_indices   = tuple([all_atoms.indices(ba) for ba in bond_atoms])

        unique_residues = self.unique_residues = all_atoms.unique_residues

        comms_coords = comms['coords'] = SharedNumpyArray(
            TypedMPArray(
                FLOAT_TYPE, coords.size)).reshape(coords.shape)
        ct.add_managed_arrays(ct.EDIT_COORDS, 'coords', (comms_coords, ), 'coords_changed_cb')

        # Container for data describing the model that will be fixed for
        # the duration of the simulation

        sim_data = self.sim_data = {
            'fixed indices':                    fixed_indices,
            'mobile indices':                   mobile_indices,
            'atom names':                       atom_names,
            'element names':                    element_names,
            'residue names':                    residue_names,
            'residue numbers':                  residue_nums,
            'chain ids':                        chain_ids,
            'coords':                           coords,
            'bonded atom indices':              bond_atom_indices,
            'tuggable indices':                 None,
            'phi atom indices':                 None,
            'psi atom indices':                 None,
            'omega atom indices':               None,
            'phi cmap resnames':                None,
            'phi cmap indices':                 None,
            'psi cmap resnames':                None,
            'psi cmap indices':                 None,
            'rotamer map':                      None,
            #'rotamer targets':                  None,
            'distance restraint keys':          None,
            'distance restraint indices':       dict(),
            'position restraint indices':       None,
            'density map names':                None,
            'density map transforms':           dict(),
            'residue templates':                dict(),
        }


        # Provide pre-defined residue templates for residues OpenMM can't
        # identify on its own.
        self.find_residue_templates()

        self._prepare_tugging_force(tuggable_atoms)

        # Secondary structure and peptide bond backbone restraints
        bd = self.backbone_dihedrals = backbone_dihedrals
        phi, psi, omega = (bd.phi, bd.psi, bd.omega)

        self.named_dihedrals = {
            'phi':      phi,
            'psi':      psi,
            'omega':    omega,
            }


        self._phi_targets, self._phi_ks = self._prepare_dihedral_restraints(phi, 'phi')
        self._psi_targets, self._psi_ks = self._prepare_dihedral_restraints(psi, 'psi')
        self._omega_targets, self._omega_ks = self._prepare_dihedral_restraints(omega, 'omega')

        self._set_all_omega_restraints_from_current_geometry(omega)

        # AMBER CMAP corrections for implicit solvent
        self._define_cmap_residues(phi, psi)

        # Prepare rotamer mappings and apply any current targets
        self._prepare_rotamers(rotamers)

        # Prepare distance restraint mappings and apply any current targets
        self._prepare_distance_restraints(distance_restraints)

        self._prepare_position_restraints(position_restraints)

        if density_maps is not None:
            self._prepare_density_maps(density_maps)

        self._pool = start_pool(sim_params, sim_data, comms, ct,
                        use_mp = self.isolde.params.use_multiprocessing)

        self._init_start_time = time()
        self.step_counter = 0
        self.starting_checkpoint = CheckPoint(self.isolde)
        self.thread_result = self._pool.apply_async(sim_thread._sim_thread, args=(), error_callback = error_cb)
        self._register_sim_event('simulation_initialization', 'session',
                                 'new frame', self._sim_init_handler)

    def _prepare_tugging_force(self, allowed_targets):
        data = self.sim_data
        comms = self.comms_object
        par = self.sim_params
        ct = self.change_tracker
        n_targets = len(allowed_targets)
        atom_indices = data['tuggable indices'] = self.all_atoms.indices(allowed_targets)
        targets = comms['tugging targets'] = SharedNumpyArray(
            TypedMPArray(FLOAT_TYPE, n_targets*3)).reshape((n_targets, 3))
        spring_constants = comms['tugging spring constants'] = SharedNumpyArray(
            TypedMPArray(FLOAT_TYPE, n_targets))
        change_bit = ct.TUG
        ct.add_managed_arrays(change_bit, 'tugging', (targets, spring_constants), 'tugging_cb')

    def _prepare_dihedral_restraints(self, dihedrals, name):
        '''
        Prepare dihedral restraint arrays, and set current targets/spring constants.
        '''
        n_dihe = len(dihedrals)
        all_atoms = self.all_atoms
        dihe_i = numpy.reshape(all_atoms.indices(dihedrals.atoms), (n_dihe, 4))
        atom_key = name + ' atom indices'
        target_key = name + ' targets'
        k_key = name + ' spring constants'
        self.sim_data[atom_key] = dihe_i
        t_array = self.comms_object[target_key] = SharedNumpyArray(TypedMPArray(FLOAT_TYPE, n_dihe))
        k_array = self.comms_object[k_key] = SharedNumpyArray(TypedMPArray(FLOAT_TYPE, n_dihe))
        t_array[:] = (dihedrals.targets * OPENMM_LENGTH_UNIT).value_in_unit(CHIMERAX_LENGTH_UNIT)
        k_array[:] = dihedrals.spring_constants

        change_tracker = self.change_tracker
        change_bit = change_tracker.DIHEDRAL_RESTRAINT
        change_tracker.add_managed_arrays(change_bit, name, (t_array, k_array), 'dihedral_restraint_cb')
        return t_array, k_array

    def _set_all_omega_restraints_from_current_geometry(self, omega, track_changes = True):
        '''
        Use the current peptide bond omega dihedral values to decide
        whether to restrain each bond to cis or trans.
        '''
        cis_offset = self.sim_params['cis_peptide_bond_cutoff_angle']/unit.radians
        o_vals = omega.values
        omega_targets = omega.targets = numpy.logical_or(
                    o_vals > cis_offset, o_vals < -cis_offset).astype(float) * pi
        self._omega_targets[:] = omega_targets
        k = self.sim_params.peptide_bond_spring_constant.value_in_unit(OPENMM_RADIAL_SPRING_UNIT)
        self._omega_ks[:] = k
        omega.spring_constants = k
        if track_changes:
            self.change_tracker.register_array_changes('omega')

    def _define_cmap_residues(self, phi, psi):
        '''
        Define residues subject to AMBER CMAP corrections.
        '''
        all_atoms = self.all_atoms
        phi_has_psi = phi.residues.indices(psi.residues)
        psi_has_phi = psi.residues.indices(phi.residues)
        good_phi = phi[phi_has_psi[phi_has_psi != -1]]
        good_psi = psi[psi_has_phi[psi_has_phi != -1]]
        self.sim_data['phi cmap resnames'] = good_phi.residues.names
        self.sim_data['psi cmap resnames'] = good_psi.residues.names
        n = len(good_psi)
        self.sim_data['phi cmap indices'] = all_atoms.indices(good_phi.atoms).reshape((n,4))
        self.sim_data['psi cmap indices'] = all_atoms.indices(good_psi.atoms).reshape((n,4))

    def _prepare_rotamers(self, rotamers):
        '''
        Prepare rotamer communication arrays and fill them with current
        target values.
        '''
        all_atoms = self.all_atoms
        comms = self.comms_object
        input_rotamer_map = {}
        input_rotamer_master = comms['rotamers'] = {}
        input_rotamer_targets = input_rotamer_master['targets'] = {}
        input_rotamer_ks = input_rotamer_master['spring constants'] = {}
        mobile_res = self.mobile_residues

        restrained_mask = SharedNumpyArray(TypedMPArray(ctypes.c_bool, len(mobile_res)))
        input_rotamer_master['restrained mask'] = restrained_mask

        for i, res in enumerate(mobile_res):
            try:
                rot= rotamers[res]
            except KeyError:
                # Non-rotameric residue
                input_rotamer_map[i] = None
                continue
            dlist = rot.dihedrals
            atoms = dlist.atoms
            indices = all_atoms.indices(atoms).reshape(len(dlist),4)
            input_rotamer_map[i] = indices

            rot_target = rot.target
            target_mp_arr = TypedMPArray(FLOAT_TYPE, rot_target.size)
            target_shared_arr = SharedNumpyArray(target_mp_arr)
            target_shared_arr[:] = rot_target
            input_rotamer_targets[i] = target_shared_arr
            rot_ks = rot.spring_constants
            k_mp_arr = TypedMPArray(FLOAT_TYPE, rot_ks.size)
            k_shared_arr = SharedNumpyArray(k_mp_arr)
            k_shared_arr[:] = rot_ks
            input_rotamer_ks[i] = k_shared_arr
            restrained_mask[i] = rot.restrained

        self.sim_data['rotamer map'] = input_rotamer_map
        change_tracker = self.change_tracker
        change_bit = change_tracker.ROTAMER_RESTRAINT
        change_tracker.add_managed_arrays(change_bit,
            'rotamers', (restrained_mask, ), 'rotamer_restraint_cb')


    def _prepare_distance_restraints(self, distance_restraints):
        '''
        Prepare distance restraint communication arrays and fill them in
        with current targets and spring constants.
        '''
        all_atoms = self.all_atoms
        distance_restraint_keys = []
        sim_data = self.sim_data
        comms = self.comms_object
        self.sim_data['distance restraint keys'] = distance_restraint_keys
        distance_restraint_indices = self.sim_data['distance restraint indices']
        change_tracker = self.change_tracker
        change_bit = change_tracker.DISTANCE_RESTRAINT

        for key, d_r in distance_restraints.items():
            n_restraints = len(d_r)
            distance_restraint_keys.append(key)
            target_key = key + ' targets'
            k_key = key + ' k'
            t_array = comms[target_key] = SharedNumpyArray(TypedMPArray(FLOAT_TYPE, n_restraints))
            k_array = comms[k_key] = SharedNumpyArray(TypedMPArray(FLOAT_TYPE, n_restraints))
            i_array = distance_restraint_indices[key] = []
            for i,r in enumerate(d_r):
                t_array[i] = r.target_distance
                k_array[i] = r.spring_constant
                i_array.append(all_atoms.indices(r.atoms))
            change_tracker.add_managed_arrays(change_bit, key, (t_array, k_array), 'distance_restraint_cb')

    def _prepare_position_restraints(self, position_restraints):
        '''
        One position restraint is defined for each heavy atom, with a
        default target of (0,0,0) and spring constant of 0. Turning one
        on involves setting the target and setting the constant to a
        positive non-zero value.
        '''
        all_atoms = self.all_atoms
        comms = self.comms_object
        data = self.sim_data
        pr_atoms = self._restrainable_atoms = position_restraints.atoms
        pr_indices = all_atoms.indices(pr_atoms)
        data['position restraint indices'] = pr_indices
        n_pr = len(pr_indices)
        pr_ks = SharedNumpyArray(TypedMPArray(FLOAT_TYPE, n_pr))
        comms['position restraint spring constants'] = pr_ks
        pr_ks[:] = position_restraints.spring_constants
        pr_targets = SharedNumpyArray(TypedMPArray(FLOAT_TYPE, n_pr*3)).reshape((n_pr,3))
        comms['position restraint targets'] = pr_targets
        pr_targets[:] = position_restraints.targets
        change_tracker = self.change_tracker
        change_bit = change_tracker.POSITION_RESTRAINT
        change_tracker.add_managed_arrays(
            change_bit, 'position restraints', (pr_targets, pr_ks), 'position_restraint_cb')

    def _prepare_density_maps(self, density_maps, normalize = False):
        '''
        Even though updating these during a simulation is not yet possible,
        we want to leave the option open. For now, all mobile heavy atoms
        will be coupled to all maps. In the future, ideally we want to
        be able to specify precisely which atoms are coupled to which
        maps.
        '''
        if self.sim_params.hydrogens_feel_maps:
            atoms = self.mobile_atoms
        else:
            atoms = self.mobile_heavy_atoms

        # Forces scaled by atomic number. Should really be by number of
        # electrons (to account for ions), but that's not available in
        # ChimeraX. Will have to add a look-up table in ISOLDE.
        default_ks = atoms.element_numbers / CARBON_ATOMIC_NUMBER

        #default_ks = atoms.elements.masses / CARBON_MASS
        atom_indices = self.all_atoms.indices(atoms)
        density_map_names = []
        comms = self.comms_object
        sim_data = self.sim_data
        sim_data['density map atoms'] = atom_indices
        sim_data['density map names'] = density_map_names
        transforms = sim_data['density map transforms']
        change_tracker = self.change_tracker
        coupling_change_bit = change_tracker.MAP_COUPLING
        data_change_bit = change_tracker.MAP_DATA
        if density_maps is not None:
            for key, imap in density_maps.items():
                name = imap.get_name()
                density_map_names.append(name)
                pad = imap.mask_cutoff
                vd, r = imap.crop_to_selection(atoms, pad, normalize)

                tf = r.xyz_to_ijk_transform.matrix
                transforms[name] = tf

                map_data_key = 'density map data - ' + name
                global_k_key = 'density map global ks - ' + name
                atom_k_key = 'density map atom ks - ' + name
                map_data = comms[map_data_key]\
                         = SharedNumpyArray(
                            TypedMPArray(FLOAT_TYPE, vd.size)).reshape(vd.shape)
                map_data[:] = vd
                global_k = comms[global_k_key] = mp.Value(FLOAT_TYPE, imap.coupling_constant)
                #default_ks = numpy.ones(len(self.mobile_heavy_atoms))
                atom_ks = comms[atom_k_key] = SharedNumpyArray(TypedMPArray(FLOAT_TYPE, default_ks))
                change_tracker.add_managed_arrays(
                    coupling_change_bit, name, (atom_ks,), 'density_map_coupling_cb')
                change_tracker.add_managed_arrays(
                    data_change_bit, map_data_key, (map_data, ), 'density_map_data_change_cb')

    def find_residue_templates(self):
        '''
        Allowing OpenMM to ignore missing bonds to adjacent residues is
        necessary, but comes with a drawback: in a small number of cases,
        the topology minus the external bond is also a valid residue,
        and OpenMM can't tell which to apply. The most common example of
        this is cysteine, where we have standard CYS, disulfide bonded
        CYS, and the rare-but-occasionally-valid deprotonated CYS. If
        the extent of our simulation construct is such that one CYS from
        a disulfide-bonded pair is left out, OpenMM needs to be told that
        this is not a deprotonated residue.
        '''
        residues = self.unique_residues
        templates = self.sim_data['residue templates']
        cys_indices = numpy.where(residues.names == 'CYS')[0]
        for c_i in cys_indices:
            templates[c_i] = cys_type(residues[c_i])








def openmm_topology_from_model(model):
    '''
    Take an AtomicStructure model from ChimeraX and return an OpenMM
    topology (e.g. for use with the OpenMM Modeller class).
    '''


    a = model.atoms
    b = model.bonds
    n = len(a)
    r = a.residues
    aname = a.names
    ename = a.element_names
    rname = r.names
    rnum = r.numbers
    cids = r.chain_ids
    from simtk.openmm.app import Topology, Element
    from simtk import unit
    top = Topology()
    cmap = {}
    rmap = {}
    atoms = {}
    for i in range(n):
        cid = cids[i]
        if not cid in cmap:
            cmap[cid] = top.addChain()   # OpenMM chains have no name
        rid = (rname[i], rnum[i], cid)
        if not rid in rmap:
            res = rmap[rid] = top.addResidue(rname[i], cmap[cid])
        element = Element.getBySymbol(ename[i])
        atoms[i] = top.addAtom(aname[i], element,rmap[rid])

    a1, a2 = b.atoms
    for i1, i2 in zip(a.indices(a1), a.indices(a2)):
        if -1 not in [i1, i2]:
            top.addBond(atoms[i1],  atoms[i2])

    return top

def cys_type(residue):
    atoms = residue.atoms
    names = atoms.names
    sulfur_atom = atoms[names == 'SG'][0]
    bonds = Bonds(sulfur_atom.bonds)
    if len(bonds) == 1:
        # Deprotonated
        return 'CYM'
    bonded_atoms = concatenate(bonds.atoms)
    for a in bonded_atoms:
        if a.residue != residue:
            if 'OXT' in names:
                return 'CCYX'
            if 'H1' in names:
                return 'NCYX'
            return 'CYX'
        if a.name == 'HG':
            if 'OXT' in names:
                return 'CCYS'
            if 'H1' in names:
                return 'NCYS'
            return 'CYS'

class available_forcefields():
    '''Main force field files'''
    main_files = (
        amber14,
        ['amber99sbildn.xml',],
        ['amber99sbnmr.xml',],
        ['amber10.xml',],
        ['charmm36.xml',],
        )

    main_file_descriptions = (
        'AMBER14 withi improved backbone & sidechain torsions',
        'AMBER99 with improved backbone & sidechain torsions',
        'AMBER99 with modifications to fit NMR data',
        'AMBER10',
        'CHARMM36',
        )

    # Implicit solvent force field files. The main file and the implicit
    # solvent file must match.
    implicit_solvent_files = [
        None,
        'amber99_obc.xml',
        'amber99_obc.xml',
        'amber10_obc.xml',
        None
        ]

    # Explicit water models
    explicit_water_files = [
        None,
        'tip3pfb.xml',
        'tip4pfb.xml',
        'tip3p.xml',
        ]

    explicit_water_descriptions = [
        None,
        'TIP3P-FB (DOI: 10.1021/jz500737m)',
        'TIP4P-FB (DOI: 10.1021/jz500737m)',
        'Original TIP3P water (not recommended)',
        ]

def get_available_platforms():
        from simtk.openmm import Platform
        platform_names = []
        for i in range(Platform.getNumPlatforms()):
            p = Platform.getPlatform(i)
            name = p.getName()
            platform_names.append(name)
        return platform_names


class GenericMapObject:
    def __init__(self, c3dfunc, per_atom_coupling, coupling_k):
        '''
        A minimal class to hold a Continuous3DFunction and its parameters.
        '''
        self._potential_function = c3dfunc
        self._per_atom_coupling = per_atom_coupling
        self._k = coupling_k

    def get_potential_function(self):
        return self._potential_function

    def get_per_atom_coupling_params(self):
        return self._k

    def per_atom_coupling(self):
        return hasattr(self._k, '__len__') or hasattr(self._k, '__iter__')
