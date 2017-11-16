import multiprocessing as mp
import numpy
import ctypes
import traceback
from ..threading import TypedMPArray, SharedNumpyArray, ThreadComms
from ..threading import ChangeTracker as ChangeTracker_Base
from .sim_handler import SimHandler
from simtk import unit, openmm
from ..constants import defaults, control
from time import sleep
FLOAT_TYPE = defaults.FLOAT_TYPE

OPENMM_LENGTH_UNIT          = defaults.OPENMM_LENGTH_UNIT
OPENMM_FORCE_UNIT           = defaults.OPENMM_FORCE_UNIT
OPENMM_SPRING_UNIT          = defaults.OPENMM_SPRING_UNIT
OPENMM_RADIAL_SPRING_UNIT   = defaults.OPENMM_RADIAL_SPRING_UNIT
OPENMM_ENERGY_UNIT          = defaults.OPENMM_ENERGY_UNIT
OPENMM_ANGLE_UNIT           = defaults.OPENMM_ANGLE_UNIT
CHIMERAX_LENGTH_UNIT        = defaults.CHIMERAX_LENGTH_UNIT
CHIMERAX_FORCE_UNIT         = defaults.CHIMERAX_FORCE_UNIT
CHIMERAX_SPRING_UNIT        = defaults.CHIMERAX_SPRING_UNIT

SIM_MODE_MIN                = control.SIM_MODE_MIN
SIM_MODE_EQUIL              = control.SIM_MODE_EQUIL
SIM_MODE_UNSTABLE           = control.SIM_MODE_UNSTABLE

PAUSE_SLEEP = defaults.THREAD_PAUSE_SLEEP_INTERVAL 

class ChangeTracker(ChangeTracker_Base):
    def __init__(self):
        #Control Inputs
        self.PAUSE_TOGGLE            = 1
        self.STOP                    = 1<<1
        self.TEMPERATURE             = 1<<2
        self.MODE                    = 1<<3
        self.TUG                     = 1<<4
        self.DIHEDRAL_RESTRAINT      = 1<<5
        self.POSITION_RESTRAINT      = 1<<6
        self.DISTANCE_RESTRAINT      = 1<<7
        self.ROTAMER_RESTRAINT       = 1<<8
        self.MAP_COUPLING            = 1<<9
        self.MAP_DATA                = 1<<10
        self.EDIT_COORDS             = 1<<11

        #Outputs
        self.MIN_COMPLETE = 1<<27
        self.INIT_COMPLETE = 1<<28
        self.COORDS_READY = 1<<29
        self.UNSTABLE = 1<<30
        self.ERROR = 1<<31

        super().__init__()


class SimComms(ThreadComms):
    '''
    Holds all the shared variables/arrays needed to communicate with the
    simulation. Those that are common to all simulations are defined on
    initialisation of the SimComms object, while the remainder will be
    added during initialisation of the simulation. Therefore, a new
    SimComms object is required for each simulation run.
    '''
    def __init__(self):
        super().__init__()
        self.update({
            'changes'                   : ChangeTracker(),
            'error'                     : mp.Queue(),
            'status'                    : mp.Queue(),
            'sim mode'                  : mp.Value(ctypes.c_int, 0),
            'temperature'               : mp.Value(FLOAT_TYPE, 0),
            'pause'                     : mp.Value(ctypes.c_bool, False),
        })

def _init_sim_thread(sim_params, sim_data, sim_comms, change_tracker):
    global _sim_thread_obj
    _sim_thread_obj = SimThread(sim_params, sim_data, sim_comms, change_tracker)

def _sim_thread():
    '''
    Main loop for the simulation thread.
    '''
    global _sim_thread_obj

    so = _sim_thread_obj
    comms = so.comms
    ct = so.change_tracker
    changes = ct.changes
    status_q = comms['status']
    error_q = comms['error']
    # Wait for ISOLDE to give the signal to go ahead:
    status_q.get()
    while True:
        with changes.get_lock():
            current_changes = changes.value
            ct.clear_inputs()
        if current_changes & ct.STOP:
            status_q.put('Simulation terminated on command.')
            break

        try:
            so.main_loop(current_changes)
        except Exception as e:
            tb = traceback.format_exc()
            main_str = e.args[0]
            raise Exception(main_str + tb)
            error_q.put((e, tb))
            sleep(0.1)
            ct.error()
            break


class SimThread:
    '''
    Designed to be run in a worker thread under the control of the
    multiprocessing module. Runs an OpenMM simulation, listening for
    input and returning coordinates and simulation information at each
    iteration.
    '''
    def __init__(self, sim_params, sim_data, sim_comms, change_tracker):
        par = self.sim_params = sim_params
        data = self.sim_data = sim_data
        comms = self.comms = sim_comms
        ct = self.change_tracker = change_tracker
        changes = ct.changes
        force_maps = self.force_index_maps = {
            'tugging':              None,
            'dihedrals':            {},
            'rotamers':             {},
            'distance restraints':  {},
            'position restraints':  None,
            'density maps':         {},
        }
        error_queue = comms['error']
        self.sim_mode = SIM_MODE_MIN
        self.last_max_force = par['max_allowable_force']
        self.pause = False


        # Counters for potential timeout conditions
        self.unstable_counter = 0
        self.MAX_UNSTABLE_ROUNDS = 20

        try:
            sh = self.sim_handler = SimHandler()

            # Just initialize the forces. They're empty at this stage.
            sh.initialize_dihedral_restraint_force(par['dihedral_restraint_cutoff_angle'])
            sh.initialize_amber_cmap_force()
            sh.initialize_distance_restraints_force(par['restraint_max_force'])
            sh.initialize_position_restraints_force(par['restraint_max_force'])
            sh.initialize_tugging_force(par['restraint_max_force'])

            # Tugging force
            self.tugging_force_map = force_maps['tugging']
            self.init_tugging()

            # Backbone dihedral restraints
            self.dihedral_force_maps = force_maps['dihedrals']
            self.init_dihedrals('phi')
            self.init_dihedrals('psi')
            self.init_dihedrals('omega', comms['omega targets'])

            #CMAP corrections
            sh.add_amber_cmap_torsions(sim_data['phi cmap resnames'],
                                        sim_data['psi cmap resnames'],
                                        sim_data['phi cmap indices'],
                                        sim_data['psi cmap indices'])

            #Rotamers
            self.rotamer_force_map = force_maps['rotamers']
            self.init_rotamers(sim_data['rotamer map'],
               comms['rotamers']['targets'],
               comms['rotamers']['spring constants'],
               par['rotamer_restraint_cutoff_angle'] / unit.radians,
               )

            #Distance restraints
            #TODO: improve abstraction here
            self.distance_restraint_force_maps = force_maps['distance restraints'] = {}
            self.init_distance_restraints(data['distance restraint keys'])

            # Position restraints
            self.position_restraints_force_map = force_maps['position restraints']
            self.init_position_restraints(data['position restraint indices'],
                                          comms['position restraint spring constants'],
                                          comms['position restraint targets']
                                          )

            # Density maps
            self.density_map_forces = {}
            self.density_map_force_maps = force_maps['density maps']
            self.init_map_forces(data['density map names'])

            top, templates = sh.create_openmm_topology(data['atom names'],
                                            data['element names'],
                                            data['residue names'],
                                            data['residue numbers'],
                                            data['chain ids'],
                                            data['bonded atom indices'],
                                            data['residue templates']
                                            )

            self.topology = top

            sh.define_forcefield(par['forcefield'])


            system_params = {
                'nonbondedMethod':      par['nonbonded_cutoff_method'],
                'nonbondedCutoff':      par['nonbonded_cutoff_distance'],
                'constraints':          par['rigid_bonds'],
                'rigidWater':           par['rigid_water'],
                'removeCMMotion':       par['remove_c_of_m_motion'],
                'ignoreExternalBonds':  True,
                'residueTemplates':     templates,
                }

            system = sh.create_openmm_system(top, system_params)
            self.system = system

            if par['use_gbsa']:

                gbsa_params = {
                    'solventDielectric':    par['gbsa_solvent_dielectric'],
                    'soluteDielectric':     par['gbsa_solute_dielectric'],
                    'SA':                   par['gbsa_sa_method'],
                    'cutoff':               par['gbsa_cutoff'],
                    'kappa':                par['gbsa_kappa'],
                    'nonbonded_method':     par['gbsa_cutoff_method'],
                    }

                sh.initialize_implicit_solvent(gbsa_params)

            else:
                sh.set_vacuum_permittivity(par['vacuum_dielectric_correction'])

            sh.set_fixed_atoms(data['fixed indices'])
            sh.register_all_forces_with_system()

            t = comms['temperature']
            temperature = t.value

            integrator = par['integrator']
            if integrator == openmm.VariableLangevinIntegrator:
                integrator_params = (
                    temperature,
                    par['friction_coefficient'],
                    par['variable_integrator_tolerance'],
                    )
            elif integrator == openmm.LangevinIntegrator:
                integrator_params = (
                    temperature,
                    par['friction_coefficient'],
                    par['fixed_integrator_timestep'],
                    )
            elif integrator == openmm.MTSIntegrator:
                raise RuntimeError('Multiple timestepping not yet implemented!')
            else:
                raise RuntimeError('Unrecognised or unsupported integrator: {}!'.format(integrator))

            sh.initialize_integrator(integrator, integrator_params)

            sh.set_platform(par['platform'])

            self.sim = sh.create_sim()

            coords = self.current_coords = data['coords']
            sh.set_initial_positions_and_velocities(coords, temperature)

            self.startup = True
            self.sim_startup_counter = 0

            with changes.get_lock():
                changes.value |= ct.INIT_COMPLETE


        except Exception as e:
            tb = traceback.format_exc()
            error_queue.put((e, tb))
            change_tracker.error()

    def main_loop(self, changes):
        par = self.sim_params
        comms = self.comms
        sh = self.sim_handler
        ct = self.change_tracker
        status_q = comms['status']
        sim_mode = self.sim_mode

        # Update parameters for all arrays that have changed
        ct.run_all_necessary_callbacks(self, changes)
        sh.update_restraints_in_context_if_needed()

        # Update coordinates from the last round AFTER checking to see
        # if the master thread has changed any.
        if not self.pause:
            comms.thread_safe_set_array_values('coords', self.current_coords.value_in_unit(CHIMERAX_LENGTH_UNIT))
            ct.register_change(ct.COORDS_READY)


        if changes & ct.PAUSE_TOGGLE:
            p = self.pause = comms['pause'].value
            if p:
                status_q.put('Paused')
            else:
                status_q.put('Resumed')

        if changes & ct.MODE:
            mode = comms['sim mode']
            sim_mode = self.sim_mode = mode.value

        if changes & ct.TEMPERATURE:
            temperature = comms['temperature'].value
            sh.set_temperature(temperature)

        if self.pause:
            sleep(PAUSE_SLEEP)
            return

        if sim_mode == SIM_MODE_EQUIL:
            max_allowed_movement = par['max_atom_movement_per_step']
            coords, fast_indices = sh.equil_step(par['sim_steps_per_gui_update'], max_allowed_movement)
            if fast_indices is not None:
                comms.thread_safe_set_value('sim mode', SIM_MODE_UNSTABLE)
                self.last_max_force = par['max_allowable_force']
                ct.register_change(ct.UNSTABLE)
                out_str = 'Simulation has become unstable! The following '\
                         +'atoms are moving too fast: {}'\
                         .format(fast_indices)
                status_q.put(out_str)

        elif sim_mode == SIM_MODE_MIN:
            max_allowed_force = par['max_allowable_force']
            steps = par['minimization_steps_per_gui_update']
            #coords, max_force, max_index = sh.get_positions_and_max_force(max_allowed_force)
            coords, max_force, max_index = sh.min_step(steps, max_allowed_force)
            if max_index != -1:
                ct.register_change(ct.UNSTABLE)
                out_str = 'Excessive force of {} on atom {}. Trying to '\
                         +'minimise...'
                out_str = out_str.format(max_force.in_units_of(CHIMERAX_FORCE_UNIT), max_index)
                status_q.put(out_str)
            elif abs(self.last_max_force - max_force) < par.minimization_convergence_tol:
                ct.register_change(ct.MIN_COMPLETE)
                self.last_max_force = max_allowed_force
            else:
                self.last_max_force = max_force

        elif sim_mode == SIM_MODE_UNSTABLE:
            stable_force_limit = par['max_stable_force']
            max_allowed_force = par['max_allowable_force']
            steps = par['minimization_steps_per_gui_update']
            coords, max_force, max_index = sh.min_step(steps, stable_force_limit)
            if max_index == -1:
                out_str = 'Forces back within stable limits. Returning to '\
                         +'equilibration.'
                status_q.put(out_str)
                comms.thread_safe_set_value('sim mode', SIM_MODE_EQUIL)
                self.unstable_counter = 0
            elif max_force < max_allowed_force and self.unstable_counter == self.MAX_UNSTABLE_ROUNDS:
                out_str = 'Maximum number of minimisation rounds reached. '\
                        +'Maximum force of {} is problematic, but potentially '\
                        +'stable. Attempting to return to equilibration...'\
                        .format(max_force)
                status_q.put(out_str)
                comms.thread_safe_set_value('sim mode', SIM_MODE_EQUIL)
                self.unstable_counter = 0
            elif self.unstable_counter >= self.MAX_UNSTABLE_ROUNDS:
                ct.register_change(ct.UNSTABLE&ct.ERROR)
            else:
                self.unstable_counter += 1

        self.current_coords = coords


    def coords_changed_cb(self, change_mask, _, arrays):
        '''
        Apply changes in coordinates to the simulation if they've been
        edited externally.
        '''
        sh = self.sim_handler
        new_coords = arrays[0]
        with change_mask.get_lock():
            changed_indices = numpy.where(change_mask)[0]
            change_mask[:] = False
        #OpenMM can only change coordinates all at once, so we'll make sure
        # we're only changing the flagged ones to avoid nasty surprises.
        with new_coords.get_lock():
            self.current_coords[changed_indices] = new_coords[changed_indices]*CHIMERAX_LENGTH_UNIT
        sh.change_coords(self.current_coords, self.sim_params['max_stable_force'])

    def tugging_cb(self, change_mask, _, arrays):
        '''
        Tug or release atoms.
        '''
        sh = self.sim_handler
        m_targets, m_ks = arrays
        with m_targets.get_lock(), change_mask.get_lock(), m_ks.get_lock():
            indices = numpy.where(change_mask)[0]
            change_mask[:] = False
            ts = m_targets[indices].copy()
            ks = m_ks[indices].copy()
        force_map = self.force_index_maps['tugging']
        sh.tug_atoms(force_map[indices], ts, ks)

    def dihedral_restraint_cb(self, change_mask, name, arrays):
        '''
        Adjust dihedral targets and spring constants
        '''
        sh = self.sim_handler
        force_map = self.dihedral_force_maps[name]
        m_ts, m_ks = arrays
        with m_ts.get_lock(), m_ks.get_lock(), change_mask.get_lock():
            indices = numpy.where(change_mask)[0]
            change_mask[:] = False
            ts = m_ts[indices].copy()
            ks = m_ks[indices].copy()
        sh.update_dihedral_restraints(force_map[indices], ts, ks)


    def rotamer_restraint_cb(self, change_mask, key, arrays):
        sh = self.sim_handler
        comms = self.comms
        master_dict = comms[key]
        target_dict = master_dict['targets']
        k_dict = master_dict['spring constants']
        #k = self.sim_params['rotamer_spring_constant'].value_in_unit(OPENMM_RADIAL_SPRING_UNIT)
        restrained_mask = arrays[0]
        with change_mask.get_lock(), restrained_mask.get_lock():
            indices = numpy.where(change_mask)[0]
            change_mask[:] = False
            restrained = restrained_mask[indices].copy()
        force_maps = self.rotamer_force_map
        for i, r in zip(indices, restrained):
            targets = target_dict[i]
            ks = k_dict[i]
            with targets.get_lock(), ks.get_lock():
                sh.update_dihedral_restraints(force_maps[i], targets, ks)

    def distance_restraint_cb(self, change_mask, key, arrays):
        sh = self.sim_handler
        comms = self.comms
        force_map = self.distance_restraint_force_maps[key]
        m_ts, m_ks = arrays
        with change_mask.get_lock(), m_ts.get_lock(), m_ks.get_lock():
            indices = numpy.where(change_mask)[0]
            change_mask[:] = False
            ts = m_ts[indices].copy()*CHIMERAX_LENGTH_UNIT
            ks = m_ks[indices].copy()*CHIMERAX_SPRING_UNIT
        sh.update_distance_restraints(force_map[indices], ts, ks)

    def position_restraint_cb(self, change_mask, key, arrays):
        sh = self.sim_handler
        comms = self.comms
        force_map = self.position_restraints_force_map
        m_ts, m_ks = arrays
        with change_mask.get_lock(), m_ts.get_lock(), m_ks.get_lock():
            indices = numpy.where(change_mask)[0]
            change_mask[:] = False
            ts = m_ts[indices].copy()*CHIMERAX_LENGTH_UNIT
            ks = m_ks[indices].copy()*CHIMERAX_SPRING_UNIT
        sh.update_position_restraints(indices, ts, ks)

    def density_map_coupling_cb(self, change_mask, key, arrays):
        sh = self.sim_handler
        comms = self.comms
        force_map = self.density_map_force_maps[key]
        force = self.density_map_forces[name]
        m_ks = arrays[0]
        with change_mask.get_lock(), m_ks.get_lock():
            indices = numpy.where(change_mask)[0]
            change_mask[:] = False
            ks = m_ks[indices].copy()
        sh.update_density_map_individual_ks(indices, ks)


    def density_map_data_change_cb(self, change_mask, key, arrays):
        raise NotImplementedError('Changing a map in a running simulation is not yet possible.')


    def init_tugging(self):
        data = self.sim_data
        sh = self.sim_handler
        tuggable_indices = data['tuggable indices']
        tfm = self.force_index_maps['tugging'] = numpy.empty(len(tuggable_indices), int)
        sh.couple_atoms_to_tugging_force(tuggable_indices, tfm)

    def init_dihedrals(self, name, target_array = None):
        sh = self.sim_handler
        comms = self.comms
        atom_key = name + ' atom indices'
        target_key = name + ' targets'
        k_key = name + ' spring constants'
        dihedral_atom_indices = self.sim_data[atom_key]
        targets = comms[target_key]
        ks = comms[k_key]
        n_dihe = len(dihedral_atom_indices)
        force_map = numpy.empty(n_dihe, int)
        with targets.get_lock(), ks.get_lock():
            for i, (ai, target, k) in enumerate(zip(dihedral_atom_indices, targets, ks)):
                force_map[i] = sh.initialize_dihedral_restraint(ai, target=target, k=k)
        self.dihedral_force_maps[name] = force_map



    def init_rotamers(self, in_map, targets, ks, cutoff):
        sh = self.sim_handler
        force_map = self.rotamer_force_map
        comms = self.comms
        master_dict = comms['rotamers']
        restrained_mask = master_dict['restrained mask']
        with restrained_mask.get_lock():
            for index, dihedrals in in_map.items():
                if dihedrals is None:
                    continue
                force_indices = force_map[index] = []
                d_targets = targets[index]
                d_ks = ks[index]
                for d_indices in dihedrals:
                    force_indices.append(sh.initialize_dihedral_restraint(d_indices, cutoff))
                with d_targets.get_lock(), d_ks.get_lock():
                    if d_targets is not None:
                        for (fi, di, t, k) in zip(force_indices, dihedrals, d_targets, d_ks):
                            sh.update_dihedral_restraint(fi, target= t, k=k)

    def init_distance_restraints(self, keys):
        comms = self.comms
        data = self.sim_data
        force_map_container = self.distance_restraint_force_maps
        sh = self.sim_handler
        for key in keys:
            t_array = comms[key + ' targets']
            k_array = comms[key + ' k']
            i_array = data['distance restraint indices'][key]
            with t_array.get_lock(), k_array.get_lock():
                force_map = force_map_container[key] = numpy.empty(len(k_array), int)
                sh.add_distance_restraints(i_array, t_array*CHIMERAX_LENGTH_UNIT, k_array*CHIMERAX_SPRING_UNIT, force_map)


    def init_position_restraints(self, atom_indices, ks, targets):
        sh = self.sim_handler
        force_maps = self.force_index_maps
        with ks.get_lock(), targets.get_lock():
            force_map = self.position_restraints_force_map =\
                force_maps['position restraints'] = numpy.empty(len(ks), int)
            sh.add_position_restraints(force_map, atom_indices, ks*CHIMERAX_SPRING_UNIT, targets*CHIMERAX_LENGTH_UNIT)

    def init_map_forces(self, names):
        force_maps = self.density_map_force_maps
        force_dict = self.density_map_forces
        comms = self.comms
        data = self.sim_data
        sh = self.sim_handler
        atom_indices = data['density map atoms']
        density_map_names = data['density map names']
        for name in density_map_names:
            map_data = comms['density map data - ' + name]
            atom_ks = comms['density map atom ks - ' + name]
            global_k = comms['density map global ks - ' + name]
            transform = data['density map transforms'][name]
            force_map = force_maps[name] = numpy.empty(len(atom_indices), int)
            with map_data.get_lock(), global_k.get_lock():
                g_k = global_k.value
                map_force = force_dict[name] = sh.map_to_force_field(
                    map_data, transform, g_k)
            with atom_ks.get_lock():
                sh.couple_atoms_to_map(atom_indices, atom_ks, map_force, force_map)
