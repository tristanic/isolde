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
from time import time
from simtk import unit
from simtk.unit import Quantity, Unit
from chimerax.core.atomic import concatenate

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


try:
    from chimerax import app_bin_dir
    import os
    spawn.set_executable(os.path.join(app_bin_dir, 'python3.6'))
    #~ mp.set_start_method('spawn')
except:
    # We're not in ChimeraX any more, Toto!
    pass

# Workaround since opening the ChimeraX shell breaks the multiprocessing
# spawn method (even for tools that aren't launched from the console).
# After opening the shell, sys.modules['__main__'] will return:
# <IPython.core.interactiveshell.DummyMod>
# which has no __spec__ attribute.
main = sys.modules['__main__']
if not hasattr(main, '__spec__'):
    main.__spec__ = None


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
            'restraint_max_force':              (defaults.MAX_RESTRAINT_FORCE, OPENMM_FORCE_UNIT),
            'haptic_spring_constant':           (defaults.HAPTIC_SPRING_CONSTANT, OPENMM_SPRING_UNIT),
            'mouse_tug_spring_constant':        (defaults.MOUSE_TUG_SPRING_CONSTANT, OPENMM_SPRING_UNIT),
            'tug_max_force':                    (defaults.MAX_TUG_FORCE, OPENMM_FORCE_UNIT),
            'dihedral_restraint_cutoff_angle':  (defaults.DIHEDRAL_RESTRAINT_CUTOFF, OPENMM_ANGLE_UNIT),
            'rotamer_restraint_cutoff_angle':   (defaults.ROTAMER_RESTRAINT_CUTOFF, OPENMM_ANGLE_UNIT),
            'rotamer_spring_constant':          (defaults.ROTAMER_SPRING_CONSTANT, OPENMM_RADIAL_SPRING_UNIT),
            'peptide_bond_spring_constant':     (defaults.PEPTIDE_SPRING_CONSTANT, OPENMM_RADIAL_SPRING_UNIT),
            'phi_psi_spring_constant':          (defaults.PHI_PSI_SPRING_CONSTANT, OPENMM_RADIAL_SPRING_UNIT),
            'cis_peptide_bond_cutoff_angle':    (defaults.CIS_PEPTIDE_BOND_CUTOFF, OPENMM_ANGLE_UNIT),
            'max_atom_movement_per_step':       (defaults.MAX_ATOM_MOVEMENT_PER_STEP, OPENMM_LENGTH_UNIT),
            'max_allowable_force':              (defaults.MAX_ALLOWABLE_FORCE, OPENMM_FORCE_UNIT),
            'max_stable_force':                 (defaults.MAX_STABLE_FORCE, OPENMM_FORCE_UNIT),
            'friction_coefficient':             (defaults.OPENMM_FRICTION, 1/OPENMM_TIME_UNIT),
            'temperature':                      (defaults.TEMPERATURE, OPENMM_TEMPERATURE_UNIT),

            'nonbonded_cutoff_method':          (defaults.OPENMM_NONBONDED_METHOD, None),
            'nonbonded_cutoff_distance':        (defaults.OPENMM_NONBONDED_CUTOFF, OPENMM_LENGTH_UNIT),

            'gbsa_cutoff_method':               (defaults.GBSA_NONBONDED_METHOD, None),
            'gbsa_solvent_dielectric':          (defaults.GBSA_SOLVENT_DIELECTRIC, OPENMM_DIPOLE_UNIT),
            'gbsa_solute_dielectric':           (defaults.GBSA_SOLUTE_DIELECTRIC, OPENMM_DIPOLE_UNIT),
            'gbsa_sa_method':                   (defaults.GBSA_SA_METHOD, None),
            'gbsa_cutoff':                      (defaults.GBSA_CUTOFF, OPENMM_LENGTH_UNIT),
            'gbsa_kappa':                       (defaults.GBSA_KAPPA, 1/OPENMM_LENGTH_UNIT),

            'rigid_bonds':                      (defaults.RIGID_BONDS, None),
            'rigid_water':                      (defaults.RIGID_WATER, None),
            'remove_c_of_m_motion':             (defaults.REMOVE_C_OF_M_MOTION, None),

            'platform':                         (defaults.OPENMM_PLATFORM, None),
            'forcefield':                       (amber14, None),
            'integrator':                       (defaults.OPENMM_INTEGRATOR_TYPE, None),
            'variable_integrator_tolerance':    (defaults.OPENMM_VAR_INTEGRATOR_TOL, None),
            'fixed_integrator_timestep':        (defaults.OPENMM_FIXED_INTEGRATOR_TS, None),
            'constraint_tolerance':             (defaults.OPENMM_CONSTRAINT_TOL, None),
            'sim_steps_per_gui_update':         (defaults.SIM_STEPS_PER_GUI_UPDATE, None),
            'minimization_steps_per_gui_update':(defaults.MIN_STEPS_PER_GUI_UPDATE, None),
            'simulation_startup_rounds':        (defaults.SIM_STARTUP_ROUNDS, None),
            'maximum_unstable_rounds':          (defaults.MAX_UNSTABLE_ROUNDS, None),
            'minimization_convergence_tol':     (defaults.MIN_CONVERGENCE_FORCE_TOL, OPENMM_FORCE_UNIT),
            'tug_hydrogens':                    (False, None),
            'hydrogens_feel_maps':              (False, None),

            }





def start_pool(sim_params, sim_data, _comms_object, _change_tracker):
    try:
        from chimerax import app_bin_dir
    except:
        return
    thread_pool = mp.Pool(processes=1, initializer=sim_thread._init_sim_thread,
        initargs=(sim_params, sim_data, _comms_object, _change_tracker))
    return thread_pool

class ChimeraXSimInterface:
    '''
    Application-facing interface between ChimeraX and a SimThread object.
    Handles interconversion between ChimeraX and generic Numpy types and
    bi-directional communication. A new instance is created for each
    simulation, and should be discarded once the simulation is done.
    '''
    TIMEOUT = 20 # seconds

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
        self._register_sim_event('ramachandran update', 'isolde',
                                 'completed simulation step',
                                  self.isolde.update_ramachandran)
        self._register_sim_event('omega update', 'isolde',
                                 'completed simulation step',
                                 self.isolde.update_omega_check)



    def _initialize_mouse_tugging(self):
        from .. import mousemodes
        isolde = self.isolde
        mt = self._mouse_tugger = mousemodes.TugAtomsMode(
            self.session, self.tuggable_atoms, isolde._annotations)
        isolde._mouse_modes.register_mode(mt.name, mt, 'right', ('control',))
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
                tug_index = self.tuggable_atoms.index(tug_atom)
                mtug.last_tugged_index = tug_index
                mtug.already_tugging = True
            else:
                tug_index = mtug.last_tugged_index
            self.tug_atom_to(tug_index, xyz0)
        elif cur_tug:
            self.release_tugged_atom(mtug.last_tugged_index)
            mtug.last_tugged_index = None
            mtug.already_tugging = False



    @property
    def sim_running(self):
        return self._pool is not None

    def toggle_pause(self):
        if not self.sim_running:
            return
        comms = self.comms_object
        ct = self.change_tracker
        currently_paused = comms['pause'].value
        comms.thread_safe_set_value('pause', not currently_paused)
        ct.register_change(ct.PAUSE_TOGGLE)
        if currently_paused:
            self.isolde.triggers.activate_trigger('simulation resumed', None)
        else:
            self.isolde.triggers.activate_trigger('simulation paused', None)

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
        if reason == sim_outcomes.DISCARD or reason == sim_outcomes.UNSTABLE:
            self.all_atoms.coords = self.starting_coords
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


    def tug_atom_to(self, index, target, spring_constant = None):
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
        key = 'rotamer targets'
        comms = self.comms_object
        ct = self.change_tracker
        r_index = self.mobile_residues.index(rotamer.residue)
        restrained_mask = ct.get_managed_arrays(key)[0]
        target_array = comms[key][r_index]
        with restrained_mask.get_lock(), target_array.get_lock():
            restrained_mask[r_index] = rotamer.restrained
            target_array[:] = rotamer.target
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


    def release_tugged_atom(self, index):
        zeros = numpy.array([0,0,0], numpy.double)
        self.tug_atom_to(index, zeros, spring_constant = 0)

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
        while status_q.full():
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
        while status_q.full():
            print(status_q.get())
        with ct.changes.get_lock():
            changes = ct.changes.value
            ct.clear_outputs()

        if self.step_counter == 10:
            self.sim_mode = 'equil'

        if self.step_counter < 10 and (changes & ct.MIN_COMPLETE):
            self.sim_mode = 'equil'

        if changes & ct.ERROR:
            if changes & ct.UNSTABLE:
                self.stop_sim(sim_outcomes.UNSTABLE)
                return

            err_q = comms['error']
            if err_q.full():
                err, traceback = err_q.get()
            else:
                err_str = '''
                Simulation thread returned an error code but not an
                exception. This should not happen.
                '''
                err = RuntimeError(err_str)
                traceback = ''
            #~ if thread.is_alive():
                #~ thread.terminate()
            print(traceback)
            self.stop_sim(sim_outcomes.DISCARD, err)

        #~ if not thread.is_alive():
            #~ session.triggers.remove_handler(self._update_handler)
            #~ print(bin(changes))
            #~ raise RuntimeError('Simulation thread terminated unexpectedly!')

        if changes & ct.UNSTABLE:
            warn_str = 'Simulation is unstable! Attempting to minimise excessive force...'
            session.logger.info(warn_str)

        if changes & ct.COORDS_READY:
            if self.time_loops:
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
        ct = self.change_tracker = ChangeTracker()

        # Container for all data that can be changed by ISOLDE, the simulation
        # thread, or both
        comms = self.comms_object = SimComms()
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

        self._pool = start_pool(sim_params, sim_data, comms, ct)

        self._init_start_time = time()
        self.step_counter = 0
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
        input_rotamer_targets = comms['rotamer targets'] = {}
        mobile_res = self.mobile_residues

        restrained_mask = SharedNumpyArray(TypedMPArray(ctypes.c_bool, len(mobile_res)))
        comms['restrained rotamers'] = restrained_mask

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
            restrained_mask[i] = rot.restrained

        self.sim_data['rotamer map'] = input_rotamer_map
        change_tracker = self.change_tracker
        change_bit = change_tracker.ROTAMER_RESTRAINT
        change_tracker.add_managed_arrays(change_bit,
            'rotamer targets', (restrained_mask, ), 'rotamer_restraint_cb')


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
        atoms = self.mobile_heavy_atoms
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
                default_ks = numpy.ones(len(self.mobile_heavy_atoms))
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
    sulfur_atom = residue.atoms.filter(residue.atoms.names == 'SG')[0]
    bonds = sulfur_atom.bonds
    if len(bonds) == 1:
        # Deprotonated
        return 'CYM'
    bonded_atoms = concatenate(bonds.atoms)
    for a in bonded_atoms:
        if a.residue != residue:
            return 'CYX'
        if a.name == 'HG':
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





# class SimHandler():
#
#     def __init__(self, session = None):
#         self.session = session
#         # Forcefield used in this simulation
#         self._forcefield = None
#         # Overall simulation topology
#         self._topology = None
#         # Atoms in simulation topology
#         self._atoms = None
#         # Atoms in ChimeraX construct
#         self._chimerax_atoms = None
#
#         # Dict holding each custom external force object, its global and its
#         # per-particle parameters. Format:
#         # {'name': [object,
#         #           {'global_name': value},
#         #           [per_particle_names],
#         #           [per_particle_default_values],
#         #           {index_in_topology: index_in_force}
#         # }
#         self._custom_external_forces = {}
#
#         # List of IsoldeMap objects registered with this simulation
#         self._maps = []
#         # CustomExternalForce handling haptic interactions
#         self._tugging_force = None
#         # Dict for mapping atom indices in the topology to indices in the tugging force
#         self._tug_force_lookup = {}
#
#
#     def update_restraints_in_context(self, context):
#         self.update_distance_restraints_in_context(context)
#         self.update_position_restraints_in_context(context)
#         self.update_dihedral_restraints_in_context(context)
#         for m in self._maps:
#             pf = m.get_potential_function()
#             pf.update_context_if_needed(context)
#
#
#     ####
#     # Positional Restraints
#     ####
#
#         ##
#         # Before simulation starts
#         ##
#
#     def initialize_position_restraints_force(self, max_force):
#         '''Just create the force object.'''
#         rf = self._position_restraints_force = TopOutRestraintForce(max_force*10)
#         return rf
#
#     def add_position_restraints(self, restraints, sim_construct):
#         '''Add all atoms in a Position_Restraints object to the force.'''
#         for r in restraints:
#             self.add_position_restraint(r, sim_construct)
#
#     def add_position_restraint(self, restraint, sim_construct):
#         '''Add one Position_Restraint object to the force.'''
#         r = restraint
#         atom = r.atom
#         sc = sim_construct
#         index = sc.index(atom)
#         rf = self._position_restraints_force
#         r.sim_handler = self
#         r.sim_force_index = rf.addParticle(index, (r.spring_constant * 100, *(r.target/10)))
#
#         ##
#         # During simulation
#         ##
#
#     def change_position_restraint_parameters(self, restraint):
#         rf = self._position_restraints_force
#         index = self._chimerax_atoms.index(restraint.atom)
#         rf.setParticleParameters(restraint.sim_force_index, index,
#             (restraint.spring_constant*100, *(restraint.target/10)))
#         rf.update_needed = True
#
#     def update_position_restraints_in_context(self, context):
#         rf = self._position_restraints_force
#         if rf.update_needed:
#             rf.updateParametersInContext(context)
#         rf.update_needed = False
#
#         ##
#         # Cleanup on simulation termination
#         ##
#
#     def disconnect_position_restraints_from_sim(self, restraints):
#         for r in restraints:
#             r.sim_force_index = -1
#             r.sim_handler = None
#
#     ####
#     # Distance Restraints
#     ####
#
#         ##
#         # Before simulation starts
#         ##
#
#     def initialize_distance_restraints_force(self, max_force):
#         tf = self._distance_restraints_force = TopOutBondForce(max_force*10)
#         return tf
#
#     def add_distance_restraints(self, restraints, sim_construct):
#         for r in restraints:
#             self.add_distance_restraint(r, sim_construct)
#
#     def add_distance_restraint(self, restraint, sim_construct):
#         r = restraint
#         atoms = r.atoms
#         indices = sim_construct.indices(atoms)
#         tf = self._distance_restraints_force
#         r.sim_handler = self
#         r.sim_force_index = tf.addBond(*indices.tolist(), (r.spring_constant*100, r.target_distance/10))
#
#         ##
#         # During simulation
#         ##
#
#     def change_distance_restraint_parameters(self, restraint):
#         tf = self._distance_restraints_force
#         indices = self._chimerax_atoms.indices(restraint.atoms).tolist()
#         tf.setBondParameters(restraint._sim_force_index, *indices,
#             (restraint.spring_constant*100, restraint.target_distance/10))
#         tf.update_needed = True
#
#     def update_distance_restraints_in_context(self, context):
#         tf = self._distance_restraints_force
#         if tf.update_needed:
#             tf.updateParametersInContext(context)
#         tf.update_needed = False
#
#         ##
#         # Cleanup on simulation termination
#         ##
#
#     def disconnect_distance_restraints_from_sim(self, restraints):
#         for r in restraints:
#             r.sim_force_index = -1
#             r.sim_handler = None
#
#
#     ####
#     # AMBER CMAP corrections
#     ####
#
#         ##
#         # Before simulation starts
#         ##
#
#     def initialize_amber_cmap_force(self):
#         self._amber_cmap_force = AmberCMAPForce()
#
#     def add_amber_cmap_torsions(self, phi_array, psi_array, sim_construct):
#         cf = self._amber_cmap_force
#         sc = sim_construct
#         phi_has_psi = phi_array.residues.indices(psi_array.residues)
#         psi_has_phi = psi_array.residues.indices(phi_array.residues)
#         good_phi = phi_array[phi_has_psi[phi_has_psi != -1]]
#         good_psi = psi_array[psi_has_phi[psi_has_phi != -1]]
#         phi_resnames = good_phi.residues.names
#         psi_resnames = good_psi.residues.names
#         n = len(good_psi)
#         phi_indices = sc.indices(good_phi.atoms).reshape((n,4))
#         psi_indices = sc.indices(good_psi.atoms).reshape((n,4))
#
#         for pname, fname, pi, fi in zip(phi_resnames, psi_resnames,
#                                         phi_indices, psi_indices):
#             assert(pname == fname)
#             cf.addTorsion(pname, pi, fi)
#
#
#
#     ####
#     # Dihedral restraints
#     ####
#
#         ##
#         # Before simulation starts
#         ##
#
#     def initialize_dihedral_restraint_force(self, default_cutoff):
#
#         self._dihedral_restraint_force = FlatBottomTorsionRestraintForce()
#         self.default_torsion_cutoff = default_cutoff
#
#
#     def initialize_dihedral_restraint(self, indices, cutoff = None):
#         #top = self._topology
#         c = (cutoff or self.default_torsion_cutoff)
#         force = self._dihedral_restraint_force
#         index_in_force = force.addTorsion(*indices.tolist(), [0, 0, cos(c)])
#         return index_in_force
#
#         ##
#         # During simulation
#         ##
#
#     def update_dihedral_restraints_in_context(self, context):
#         rf = self._dihedral_restraint_force
#         if rf.update_needed:
#             rf.updateParametersInContext(context)
#         rf.update_needed = False
#
#     def set_dihedral_restraints(self, dihedrals, target, k, degrees = False, cutoffs = None):
#         c = (cutoffs or [self.default_torsion_cutoff]*len(dihedrals))
#         variable_t = hasattr(target, '__iter__')
#         variable_k = hasattr(k, '__iter__')
#         variable_c = hasattr(c, '__iter__')
#         for i, d in enumerate(dihedrals):
#             if variable_t:
#                 t = target[i]
#             else:
#                 t = target
#             if degrees:
#                 t = radians(t)
#             if variable_k:
#                 thisk = k[i]
#             else:
#                 thisk = k
#             if variable_c:
#                 thisc = c[i]
#             else:
#                 thisc = c
#
#             self.update_dihedral_restraint(d.sim_index, target=t, k=thisk, cutoff=thisc)
#
#     def update_dihedral_restraint(self, sim_index, target = None,
#                             k = None, cutoff = None, degrees = False):
#         if target is not None and degrees:
#             target = radians(target)
#         self._dihedral_restraint_force.update_target(sim_index, target=target, k=k, cutoff=cutoff)
#
#     def register_custom_external_force(self, name, force, global_params,
#                                         per_particle_params,
#                                         per_particle_default_vals):
#         self._custom_external_forces[name] = [force, global_params, per_particle_params, per_particle_default_vals, {}]
#
#     def get_custom_external_force_by_name(self, name):
#         return self._custom_external_forces[name]
#
#     def get_all_custom_external_forces(self):
#         return self._custom_external_forces
#
#     def set_custom_external_force_particle_params(self, name, index, params):
#         fparams = self._custom_external_forces[name]
#         force = fparams[0]
#         index_lookup = fparams[4]
#         index_in_force = index_lookup[index]
#         force.setParticleParameters(index_in_force, index, params)
#
#     def register_map(self, map_object):
#         self._maps.append(map_object)
#
#     #TODO: Define the tugging force as a subclass in custom_forces.py
#     def initialize_tugging_force(self):
#         from simtk import openmm as mm
#         potential_equation = '0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)'
#         per_particle_parameters = ['k','x0','y0','z0']
#         per_particle_defaults = [0,0,0,0]
#         global_parameters = None
#         global_defaults = None
#         f = self._tugging_force = mm.CustomExternalForce(potential_equation)
#         if global_parameters is not None:
#             for p, v in zip(global_parameters, g_vals):
#                 f.addGlobalParameter(g, v)
#         if per_particle_parameters is not None:
#             for p in per_particle_parameters:
#                 f.addPerParticleParameter(p)
#         self.register_custom_external_force('tug', f, global_parameters,
#                 per_particle_parameters, per_particle_defaults)
#
#         return f
#
#     def _openmm_topology_and_external_forces(self, sim_data,
#                                             tug_hydrogens = False,
#                                             hydrogens_feel_maps = False):
#         aname   = sim_data['atom names']
#         ename   = sim_data['element names']
#         rname   = sim_data['residue names']
#         rnum    = sim_data['residue numbers']
#         cids    = sim_data['chain ids']
#         coords  = sim_data['coords']
#         bond_i  = sim_data['bonded atom indices']
#
#         fixed_indices = sim_data['fixed indices']
#         fixed_flags = numpy.zeros(len(aname),numpy.bool)
#         fixed_flags[fixed_indices] = True
#
#         from simtk.openmm.app import Topology, Element
#         from simtk import unit
#         top = self._simulation_topology = Topology()
#         cmap = {}
#         rmap = {}
#         atoms = self._atoms = {}
#         for i in range(n):
#             cid = cids[i]
#             if not cid in cmap:
#                 cmap[cid] = top.addChain()   # OpenMM chains have no name
#             rid = (rname[i], rnum[i], cid)
#             rid = (rname[i], rnum[i], cid)
#             if not rid in rmap:
#                 res = rmap[rid] = top.addResidue(rname[i], cmap[cid])
#                 if rname[i] == 'CYS':
#                     ctype = cys_type(r[i])
#                     if ctype != 'CYS':
#                         templates[res] = ctype
#
#             element = Element.getBySymbol(ename[i])
#             atoms[i] = top.addAtom(aname[i], element,rmap[rid])
#
#             if not fixed_flags[i]:
#                 # Register atoms with forces
#                 if ename[i] is not 'H' or (ename[i] is 'H' and tug_hydrogens):
#                     # All CustomExternalForces
#                     for key, ff in self._custom_external_forces.items():
#                         f = ff[0]
#                         per_particle_param_vals = ff[3]
#                         index_map = ff[4]
#                         index_map[i] = f.addParticle(i, per_particle_param_vals)
#
#                 if ename[i] is not 'H' or (ename[i] is 'H' and hydrogens_feel_maps):
#                     # All map forces
#                     for m in self._maps:
#                         self.couple_atom_to_map(i, m)
#
#         for i1, i2 in zip(*bond_i):
#             top.addBond(atoms[i1],  atoms[i2])
#
#         return top, templates
#
#
#
#
#
#     def openmm_topology_and_external_forces(self, sim_construct,
#                                         sim_bonds, fixed_flags,
#                                         tug_hydrogens = False,
#                                         hydrogens_feel_maps = False):
#         '''
#         Prepares the openmm topology and binds atoms to existing force fields.
#         Since looping over all atoms can take a long time, it's best to do
#         topology generation and force binding as a single concerted loop as
#         much as we can. This should be possible for all external-type forces
#         (tugging and map forces). Forces involving two or more atoms (e.g.
#         H-bond or dihedral restraints) will have to be handled in separate
#         loops.
#         Args:
#             sim_construct:
#                 The atoms to be simulated
#             sim_bonds:
#                 The set of all bonds between simulated atoms
#             fixed_flags:
#                 A boolean array indicating which atoms will be fixed. No
#                 need to add these to any custom forces.
#             tug_hydrogens:
#                 Do we want to be able to interactively pull on hydrogens?
#             hydrogens_feel_maps:
#                 Do we want the hydrogens to be pulled into the maps?
#         '''
#         # When we allow missing external bonds, some residues become ambiguous.
#         # In particular, a cysteine with a bare sulphur might be part of a
#         # disulphide bond but have the connecting residue missing, or may be
#         # a deprotonated (negatively-charged) Cys. In such cases, we need to
#         # explicitly tell OpenMM which templates to use.
#         templates = {}
#         a = self._chimerax_atoms = sim_construct
#         n = len(a)
#         r = a.residues
#         aname = a.names
#         ename = a.element_names
#         rname = r.names
#         rnum = r.numbers
#         cids = r.chain_ids
#         from simtk.openmm.app import Topology, Element
#         from simtk import unit
#         top = self._simulation_topology = Topology()
#         cmap = {}
#         rmap = {}
#         atoms = self._atoms = {}
#         for i in range(n):
#             cid = cids[i]
#             if not cid in cmap:
#                 cmap[cid] = top.addChain()   # OpenMM chains have no name
#             rid = (rname[i], rnum[i], cid)
#             rid = (rname[i], rnum[i], cid)
#             if not rid in rmap:
#                 res = rmap[rid] = top.addResidue(rname[i], cmap[cid])
#                 if rname[i] == 'CYS':
#                     ctype = cys_type(r[i])
#                     if ctype != 'CYS':
#                         templates[res] = ctype
#
#             element = Element.getBySymbol(ename[i])
#             atoms[i] = top.addAtom(aname[i], element,rmap[rid])
#
#             if not fixed_flags[i]:
#                 # Register atoms with forces
#                 if ename[i] is not 'H' or (ename[i] is 'H' and tug_hydrogens):
#                     # All CustomExternalForces
#                     for key, ff in self._custom_external_forces.items():
#                         f = ff[0]
#                         per_particle_param_vals = ff[3]
#                         index_map = ff[4]
#                         index_map[i] = f.addParticle(i, per_particle_param_vals)
#
#                 if ename[i] is not 'H' or (ename[i] is 'H' and hydrogens_feel_maps):
#                     # All map forces
#                     for m in self._maps:
#                         self.couple_atom_to_map(i, m)
#
#
#         a1, a2 = sim_bonds.atoms
#         for i1, i2 in zip(a.indices(a1), a.indices(a2)):
#             if -1 not in [i1, i2]:
#                 top.addBond(atoms[i1],  atoms[i2])
#
#         pos = a.coords # in Angstrom (convert to nm for OpenMM)
#         return top, pos, templates
#
#
#     def couple_atom_to_map(self, index, map_object):
#         '''
#         Adds an atom to a map-derived potential field.
#         The argument per_atom_coupling must be either a single value,
#         or an array with one value per atom
#         '''
#         m = map_object
#         map_field = m.get_potential_function()
#         k = m.get_per_atom_coupling_params()
#         if m.per_atom_coupling():
#             k = k[index]
#         map_field.addBond([index],[k])
#
#
#     #######################
#     # OLD VERSIONS
#     #######################
#
#
#     #def continuous3D_from_volume(self, volume):
#         #'''
#         #Takes a volumetric map and uses it to generate an OpenMM
#         #Continuous3DFunction. Returns the function.
#         #'''
#         #import numpy as np
#         #vol_data = volume.data
#         #mincoor = np.array(vol_data.origin)
#         #maxcoor = mincoor + volume.data_origin_and_step()[1]*(np.array(vol_data.size)-1)
#         ##Map data is in Angstroms. Need to convert (x,y,z) positions to nanometres
#         #mincoor = mincoor/10
#         #maxcoor = maxcoor/10
#         ## Continuous3DFunction expects the minimum and maximum coordinates as
#         ## arguments xmin, xmax, ymin, ...
#         #minmax = [val for pair in zip(mincoor, maxcoor) for val in pair]
#         #vol_data_1d = np.ravel(vol_data.matrix(), order = 'C').astype(np.double)
#         #vol_dimensions = (vol_data.size)
#         #print('Volume dimensions: {}; expected number: {}; actual number: {}'\
#                 #.format(vol_dimensions, np.product(vol_dimensions), len(vol_data_1d)))
#         #print('Max: {}, min: {}, nans: {}, infs: {}'.format(
#             #vol_data_1d.max(), vol_data_1d.min(),
#             #np.argwhere(np.isnan(vol_data_1d)),
#             #np.argwhere(np.isinf(vol_data_1d))))
#         #return Continuous3DFunction(*vol_dimensions, vol_data_1d, *minmax)
#
#     #def map_potential_force_field(self, c3d_func, global_k):
#         #'''
#         #Takes a Continuous3DFunction and returns a CustomCompoundBondForce
#         #based on it.
#         #Args:
#             #c3d_func:
#                 #A Continuous3DFunction
#             #global_k:
#                 #An overall global spring constant coupling atoms to the
#                 #map. This can be further adjusted per atom using
#                 #the "individual_k" parameter defined in the
#                 #CustomCompoundBondForce energy function.
#         #'''
#         #from simtk.openmm import CustomCompoundBondForce
#         #f = CustomCompoundBondForce(1,'')
#         #f.addTabulatedFunction(name = 'map_potential', function = c3d_func)
#         #f.addGlobalParameter(name = 'global_k', defaultValue = global_k)
#         #f.addPerBondParameter(name = 'individual_k')
#         #f.setEnergyFunction('-global_k * individual_k * map_potential(x1,y1,z1)')
#         #return f
#
#     #######################
#     # /OLD VERSIONS
#     #######################
#
#     def map_to_force_field(self, imap, atoms,
#                             pad, normalize = True):
#         '''
#         Takes an IsoldeMap object, masks down the map to cover the atoms
#         with the desired padding, and converts it into a LinearInterpMapForce
#         object.
#         '''
#         vd, r = imap.crop_to_selection(atoms, pad, normalize)
#         tf = r.xyz_to_ijk_transform.matrix
#         f = self._map_to_force_field(vd, r, tf, imap.get_coupling_constant())
#         imap.set_potential_function(f)
#
#     def _map_to_force_field(self, volume_data, region, xyz_to_ijk_transform,
#                             coupling_constant):
#         vd = volume_data
#         r = region
#         tf = xyz_to_ijk_transform
#         f = LinearInterpMapForce(vd, tf, units='angstroms')
#         f.set_global_k(coupling_constant)
#         return f
#
#
#     def continuous3D_from_volume(self, vol_data):
#         '''
#         Takes a volumetric map and uses it to generate an OpenMM
#         Continuous3DFunction. Returns the function.
#         '''
#         import numpy as np
#         vol_dimensions = vol_data.shape[::-1]
#         mincoor = np.array([0,0,0], np.double)
#         maxcoor = (np.array(vol_dimensions, np.double) - 1) / 10
#         # Continuous3DFunction expects the minimum and maximum coordinates as
#         # arguments xmin, xmax, ymin, ...
#         minmax = [val for pair in zip(mincoor, maxcoor) for val in pair]
#         vol_data_1d = np.ravel(vol_data, order = 'C')
#         from simtk.openmm.openmm import Continuous3DFunction
#         return Continuous3DFunction(*vol_dimensions, vol_data_1d, *minmax)
#
#     def continuous3D_from_maps(self, master_map_list, keys, atoms, pad, normalize):
#         '''
#         Combine multiple maps into a single Continuous3DFunction by
#         summing their data with applied weights.
#         '''
#         from . import volumetric
#         combined_map = None
#         for key in keys:
#             m = master_map_list[key]
#             vd, r = m.crop_to_selection(atoms, pad, normalize)
#             weighted_map = vd*m.get_coupling_constant()
#             if combined_map is None:
#                 combined_map = weighted_map
#             else:
#                 combined_map += weighted_map
#         c3d = self.continuous3D_from_volume(combined_map)
#         f = self.map_potential_force_field(c3d, 1.0, r)
#         ret = GenericMapObject(f, False, 1.0)
#         return ret
#
#
#     def map_potential_force_field(self, c3d_func, global_k, region):
#         '''
#         Takes a Continuous3DFunction and returns a CustomCompoundBondForce
#         based on it.
#         Args:
#             c3d_func:
#                 A Continuous3DFunction
#             global_k:
#                 An overall global spring constant coupling atoms to the
#                 map. This can be further adjusted per atom using
#                 the "individual_k" parameter defined in the
#                 CustomCompoundBondForce energy function.
#             xyz_to_ijk_transform:
#                 The affine transformation matrix mapping (x,y,z) coordinates
#                 back to (i,j,k) in the c3d_func array
#         '''
#         xyz_to_ijk_transform = region.xyz_to_ijk_transform
#         from simtk.openmm import CustomCompoundBondForce
#         f = CustomCompoundBondForce(1,'')
#         f.addTabulatedFunction(name = 'map_potential', function = c3d_func)
#         f.addGlobalParameter(name = 'global_k', defaultValue = global_k)
#         f.addPerBondParameter(name = 'individual_k')
#         tf = xyz_to_ijk_transform.matrix
#         tf[:,3] /= 10 # OpenMM in nm, ChimeraX in Angstroms
#         i_str = 'x1* {} + y1 * {} + z1 * {} + {}'.format(
#             tf[0][0], tf[0][1], tf[0][2], tf[0][3])
#         j_str = 'x1* {} + y1 * {} + z1 * {} + {}'.format(
#             tf[1][0], tf[1][1], tf[1][2], tf[1][3])
#         k_str = 'x1* {} + y1 * {} + z1 * {} + {}'.format(
#             tf[2][0], tf[2][1], tf[2][2], tf[2][3])
#
#         f.setEnergyFunction('-global_k * individual_k * map_potential({},{},{})'.format(
#         i_str, j_str, k_str))
#         return f
#
#
#     def update_force_in_context(self, force_name, context):
#         force = self._custom_external_forces[force_name][0]
#         force.updateParametersInContext(context)
#
#
# def define_forcefield (forcefield_list):
#     from simtk.openmm.app import ForceField
#     ff = self_forcefield = ForceField(*[f for f in forcefield_list if f is not None])
#     return ff
#
# def initialize_implicit_solvent(system, top):
#     '''Add a Generalised Born Implicit Solvent (GBIS) formulation.'''
#     # Somewhat annoyingly, OpenMM doesn't store atomic charges in a
#     # nice accessible format. So, we have to pull it back out of the
#     # NonbondedForce term.
#     for f in system.getForces():
#         if isinstance(f, NonbondedForce):
#             break
#     charges = []
#     for i in range(f.getNumParticles()):
#         charges.append(f.getParticleParameters(i)[0])
#     gbforce = GBSAForce()
#     params = GBSAForce.getStandardParameters(top)
#     print('charge length: {}, param length: {}, natoms: {}, top_n: {}'.format(
#         len(charges), len(params), f.getNumParticles(), top.getNumAtoms()))
#     i = 0
#     for charge, param in zip(charges, params):
#         gbforce.addParticle([charge, *param])
#         i+= 1
#     print(i)
#     gbforce.finalize()
#     print('GB Force num particles: '.format(gbforce.getNumParticles()))
#     system.addForce(gbforce)
#
#
#
# def create_openmm_system(top, ff, templatedict, force_implicit=False):
#     from simtk.openmm import app
#     from simtk import openmm as mm
#     from simtk import unit
#
#     params = {
#         'nonbondedMethod':      app.CutoffNonPeriodic,
#         'nonbondedCutoff':      1.0*unit.nanometers,
#         'constraints':          app.HBonds,
#         'rigidWater':           True,
#         'removeCMMotion':       False,
#         'residueTemplates':     templatedict,
#         'ignoreExternalBonds':  True,
#         }
#
#
#     try:
#         #~ system = ff.createSystem(top,
#                                 #~ nonbondedMethod = app.CutoffNonPeriodic,
#                                 #~ nonbondedCutoff = 1.0*unit.nanometers,
#                                 #~ constraints = app.HBonds,
#                                 #~ rigidWater = True,
#                                 #~ removeCMMotion = False,
#                                 #~ residueTemplates = templatedict,
#                                 #~ ignoreExternalBonds = True)
#         system = ff.createSystem(top, **params)
#         if force_implicit:
#             print('Setting up implicit solvent...')
#             initialize_implicit_solvent(system, top)
#     except ValueError as e:
#         raise Exception('Missing atoms or parameterisation needed by force field.\n' +
#                               'All heavy atoms and hydrogens with standard names are required.\n' +
#                               str(e))
#     return system
#
#
# def integrator(i_type, temperature, friction, tolerance, timestep):
#     from simtk import openmm as mm
#     if i_type == 'variable':
#         integrator = mm.VariableLangevinIntegrator(temperature, friction, tolerance)
#     elif i_type == 'fixed':
#         integrator = mm.LangevinIntegrator(temperature, friction, timestep)
#     return integrator
#
# def platform(name):
#     from simtk.openmm import Platform
#     return Platform.getPlatformByName(name)
#
#
# def create_sim(topology, system, integrator, platform):
#     from simtk.openmm.app import Simulation
#     return Simulation(topology, system, integrator, platform)
