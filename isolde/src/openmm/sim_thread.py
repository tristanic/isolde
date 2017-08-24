import multiprocessing as mp
import numpy
import ctypes
import traceback
from ..threading import TypedMPArray, SharedNumpyArray
from .sim_handler import SimHandler
from simtk import unit, openmm
from ..constants import defaults

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

class ChangeTracker:
    '''
    A simple bit-field to handle notifications for I/O between the 
    master session and the simulation thread. The first 16 bits are 
    dedicated to inputs; the last 16 bits to outputs.
    '''
    #Control Inputs
    PAUSE_TOGGLE            = 1
    STOP                    = 1<<1
    TEMPERATURE             = 1<<2
    MODE                    = 1<<3
    TUG                     = 1<<4
    DIHEDRAL_RESTRAINT      = 1<<5
    POSITION_RESTRAINT      = 1<<6
    DISTANCE_RESTRAINT      = 1<<7
    ROTAMER_RESTRAINT       = 1<<8
    MAP_COUPLING            = 1<<9
    MAP_DATA                = 1<<10
    
    #Outputs
    INIT_COMPLETE = 1<<28
    COORDS_READY = 1<<29
    UNSTABLE = 1<<30
    ERROR = 1<<31
    
    ALL_OUTPUTS           = 0xffff0000
    ALL_INPUTS            = 0x0000ffff
    
    def __init__(self):
        self._change_var = mp.Value(ctypes.c_int32, 0)
        self._array_change_keys = dict()
        self._array_change_flags = dict()
        self._array_change_masks = dict()
        self._array_change_bits = dict()
    
    def clear_inputs(self):
        self._change_var.value &= ~self.ALL_INPUTS
    
    def clear_outputs(self):
        self._change_var.value &= ~self.ALL_OUTPUTS
    
    def add_managed_array(self, change_bit_mask, key, array):
        '''
        Registers a shared array with the change tracker. Once registered,
        call register_array_changes(key, [change_mask]) to let the manager
        know that the array has changed, and optionally which elements
        have changed.
            
        @param change_bit_mask:
            The bit (e.g. ChangeTracker.DISTANCE_RESTRAINT) to be flipped
            to notify the top-level change tracker
        @param key:
            The key to the array in the SimComms dict
        @param array:
            The array itself.
        '''
        self._array_change_flags[key] = mp.Value(ctypes.c_bool, False)
        self._array_change_masks[key] = SharedNumpyArray(
                TypedMPArray(ctypes.c_bool, len(array)))
        self._array_change_bits[key] = change_bit_mask
        change_keys = self._array_change_keys
        try:
            key_list = change_keys[change_bit_mask]
        except KeyError:
            key_list = change_keys[change_bit_mask] = list()
        key_list.append(key)
    
    def register_array_changes(self, key, change_mask = None):
        '''
        Let the change tracker know that an array has changed.
        @param key:
            The key to this array in the SimComms dict
        @param change_mask:
            A numpy Boolean array of the same length as the first 
            dimension of the array, True where an element has changed.
            If not provided, it will be assumed that all elements of
            the array have changed.
        '''
        shared_flag = self._array_change_flags[key]
        shared_mask = self._array_change_masks[key]
        with shared_flag.get_lock(), shared_mask.get_lock():
            shared_flag.value = True
            if change_mask is not None:
                shared_mask[:] = change_mask
            else:
                shared_mask[:] = True
        self.changes.value |= self._array_change_bits[key]
        
    
    def clear_array_changes(self, key):
        shared_flag = self._array_change_flags[key]
        shared_mask = self._array_change_masks[key]
        with shared_flag.get_lock(), shared_mask.get_lock():
            shared_flag.value = False
            shared_mask[:] = False
    
    @property
    def changes(self):
        return self._change_var
    
    def error(self):
        self.changes.value |= self.ERROR
    
    def register_change(self, change):
        changes = self.changes
        with changes.get_lock():
            changes.value |= change
    
class SimComms:
    '''
    Holds all the shared variables/arrays needed to communicate with the
    simulation. Those that are common to all simulations are defined on
    initialisation of the SimComms object, while the remainder will be
    added during initialisation of the simulation. Therefore, a new 
    SimComms object is required for each simulation run.
    '''
    def __init__(self):
        manager = self.manager = mp.Manager()
        self._comms_obj = {
            'changes'                   : ChangeTracker(),
            'error'                     : mp.Queue(),
            'status'                    : mp.Queue(),
            'sim mode'                  : mp.Value('i', 0),
            'temperature'               : mp.Value(FLOAT_TYPE, 0),
        }
    
    def __getitem__(self, key):
        return self._comms_obj[key]
    
    def __setitem__(self, key, obj):
        self._comms_obj[key] = obj
    


def _sim_thread(sim_params, sim_data, sim_comms, change_tracker):
    '''
    Initialisation function for the simulation thread.
    '''
    status_q = sim_comms['status']
    error_q = sim_comms['error']
    sim_thread = SimThread(sim_params, sim_data, sim_comms, change_tracker)
    # Wait for ISOLDE to give the signal to go ahead:
    status_q.get()
    changes = change_tracker.changes
    while True:
        with changes.get_lock():
            current_changes = changes.value
            change_tracker.clear_inputs()
        if current_changes & ct.STOP:
            with error_q.get_lock():
                error_q.put((Exception(),'Stopped on command: {}'.format(bin(current_changes))))
            change_tracker.error()
            break
            
        try:
            sim_thread.main_loop(current_changes)
        except Exception as e:
            tb = traceback.format_exc()
            with error_q.get_lock():
                error_q.put((e, tb))
            change_tracker.error()
            break
            

class SimThread:
    '''
    Designed to be run in a worker thread under the control of the
    multiprocessing module. Runs an OpenMM simulation, listening for
    input and returning coordinates and simulation information at each
    iteration.
    '''
    SIM_MODE_MIN = 0
    SIM_MODE_EQUIL = 1
    def __init__(self, sim_params, sim_data, sim_comms, change_tracker):
        par = self.sim_params = sim_params
        data = self.sim_data = sim_data
        comms = self.comms = sim_comms
        ct = self.change_tracker = change_tracker
        changes = ct.changes
        force_maps = self.force_index_maps = {}
        error_queue = comms['error']
        self.sim_mode = self.SIM_MODE_MIN
        try:
            sh = self.sim_handler = SimHandler()
            
            # Just initialize the forces. They're empty at this stage.
            sh.initialize_dihedral_restraint_force(par['dihedral_restraint_cutoff_angle'])
            sh.initialize_amber_cmap_force()
            sh.initialize_distance_restraints_force(par['restraint_max_force'])
            sh.initialize_position_restraints_force(par['restraint_max_force'])
            sh.initialize_tugging_force()
            
            # Tugging force
            mobile_indices = data['mobile indices']
            tfm = self.tug_force_map = numpy.empty(len(mobile_indices), int)            
            sh.couple_atoms_to_tugging_force(mobile_indices, tfm)
            
            # Backbone dihedral restraints
            self.dihedral_force_maps = force_maps['dihedrals'] = {}
            self.init_dihedrals('phi')
            self.init_dihedrals('psi')
            self.init_dihedrals('omega', comms['omega targets'])
            
            #CMAP corrections
            sh.add_amber_cmap_torsions(sim_data['phi cmap resnames'],
                                        sim_data['psi cmap resnames'],
                                        sim_data['phi cmap indices'],
                                        sim_data['psi cmap indices'])
            
            #Rotamers
            self.rotamer_force_map = force_maps['rotamers'] = {}
            self.init_rotamers(sim_data['rotamer map'], 
               comms['rotamer targets'],
               par['rotamer_restraint_cutoff_angle'] / unit.radians,
               par['rotamer_spring_constant'] / (unit.kilojoule_per_mole/unit.radians**2),
               )
            
            #Distance restraints
            #TODO: improve abstraction here
            self.distance_restraint_force_maps = force_maps['distance restraints'] = {}
            self.init_distance_restraints(data['distance restraint keys'], 
                                          sh._distance_restraints_force
                                         )
            
            # Position restraints
            self.position_restraints_force_map = force_maps['position restraints'] = None
            self.init_position_restraints(data['position restraint indices'],
                                          comms['position restraint spring constants'],
                                          comms['position restraint targets']
                                          )
            
            # Density maps
            self.density_map_forces = {}
            self.density_map_force_maps = force_maps['density maps'] = {}
            self.init_map_forces(data['density map names'])
            
            top = sh.create_openmm_topology(data['atom names'],
                                            data['element names'],
                                            data['residue names'],
                                            data['residue numbers'],
                                            data['chain ids'],
                                            data['bonded atom indices']
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
                'residueTemplates':     data['residue templates'],
                }
            
            system = sh.create_openmm_system(top, system_params)
            self.system = system

            gbsa_params = {
                'solventDielectric':    par['gbsa_solvent_dielectric'],
                'soluteDielectric':     par['gbsa_solute_dielectric'],
                'SA':                   par['gbsa_sa_method'],
                'cutoff':               par['gbsa_cutoff'],
                'kappa':                par['gbsa_kappa'],
                'nonbonded_method':     par['gbsa_cutoff_method'],
                }

            sh.initialize_implicit_solvent(gbsa_params)

            
            sh.set_fixed_atoms(data['fixed indices'])
            sh.register_all_forces_with_system()
            
            t = comms['temperature']
            with t.get_lock():
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
            
            coords = data['coords']    
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
        comms = self.comms_object
        sh = self.sim_handler
        ct = self.change_tracker
        status_q = comms['status']
        
        if changes & ct.MODE:
            mode = comms['sim mode']
            with mode.get_lock():
                self.sim_mode = mode.value
        
        sim_mode = self.sim_mode
            
        if sim_mode == self.SIM_MODE_EQUIL:
            max_allowed_movement = par['max_atom_movement_per_step']
            coords, fast_indices = sh.equil_step(par['sim_steps_per_gui_update'], max_allowed_movement)
            if fast_indices:
                ct.register_change(ct.UNSTABLE)
                out_str = 'Simulation has become unstable! The following '\
                         +'atoms are moving too fast: {}'\
                         .format(fast_indices)
                with status_q.get_lock():
                    status_q.put(out_str)
                
            
        elif sim_mode == self.SIM_MODE_MIN:
            max_allowed_force = par['max_allowable_force']
            steps = par['minimization_steps_per_gui_update']
            coords, max_force, max_index = sh.min_step(steps, max_allowed_force)
            if max_index != -1:
                ct.register_change(ct.UNSTABLE)
                out_str = 'Excessive force of {} on atom {}. Trying to '\
                         +'minimise...'\
                         .format(max_force.in_units_of(CHIMERAX_FORCE_UNIT), max_index)
                with status_q.get_lock():
                    status_q.put(out_str)
        
        master_coords = comms['coords']
        with master_coords.get_lock():
            master_coords[:] = coords.value_in_unit(CHIMERAX_LENGTH_UNIT)
        with ct.changes.get_lock():
            ct.changes.value |= ct.COORDS_READY
                
        
    def init_dihedrals(self, name, target_array = None):
        sh = self.sim_handler
        atom_key = name + ' atom indices'
        dihedral_atom_indices = self.sim_data[atom_key]
        n_dihe = len(dihedral_atom_indices)
        force_map = numpy.empty(n_dihe, int)
        for i, ai in enumerate(dihedral_atom_indices):
            force_map[i] = sh.initialize_dihedral_restraint(ai)
        if target_array is not None:
            with target_array.get_lock():
                for si, target in zip(force_map, target_array):
                    sh.update_dihedral_restraint(si, target, 
                        k = self.sim_params['peptide_bond_spring_constant'])
        self.dihedral_force_maps[name] = force_map
    
    def init_rotamers(self, in_map, targets, cutoff, k):
        sh = self.sim_handler
        force_map = self.rotamer_force_map
        for index, dihedrals in in_map.items():
            if dihedrals is None:
                continue
            force_indices = force_map[index] = []
            for d_indices in dihedrals:
                force_indices.append(sh.initialize_dihedral_restraint(d_indices, cutoff))
            d_targets = targets[index]
            with d_targets.get_lock():
                if d_targets is not None:
                    for (fi, di, t) in zip(force_indices, dihedrals, d_targets):
                        sh.update_dihedral_restraint(fi, target= t, k=k)
                    
    def init_distance_restraints(self, keys, force):
        comms = self.comms
        data = self.sim_data
        force_map_container = self.distance_restraint_force_maps
        for key in keys:
            t_array = comms[key + ' targets']
            k_array = comms[key + ' k']
            i_array = data['distance restraint indices'][key]
            with t_array.get_lock(), k_array.get_lock():
                force_map = force_map_container[key] = numpy.empty(len(k_array), int)
                for i, (indices, t, k) in enumerate(zip(i_array, t_array, k_array)):
                    force_map[i] = force.addBond(*indices, (k, t))


    def init_position_restraints(self, atom_indices, ks, targets):
        sh = self.sim_handler
        force_maps = self.force_index_maps
        with ks.get_lock(), targets.get_lock():
            force_map = self.position_restraints_force_map =\
                force_maps['position restraints'] = numpy.empty(len(ks), int)
            sh.add_position_restraints(force_map, atom_indices, ks, targets)
    
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
            force_map = force_maps[name]
            with map_data.get_lock(), global_k.get_lock():
                g_k = global_k.value
                map_force = force_dict[name] = sh._map_to_force_field(
                    map_data, transform, g_k)
            with atom_ks.get_lock():
                sh.couple_atoms_to_map(atom_indices, atom_ks, map_force, force_map)
    
        
