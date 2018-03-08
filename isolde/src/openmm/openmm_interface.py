import os, sys, glob
import numpy
import ctypes
from chimerax.core.state import State
from chimerax.core.atomic import molc
from chimerax.core.atomic.molc import CFunctions, string, cptr, pyobject, \
    set_c_pointer, pointer, size_t
# object map lookups
from chimerax.core.atomic.molobject import _atoms, \
                _atom_pair, _atom_or_none, _bonds, _chain, _element, \
                _pseudobonds, _residue, _residues, _rings, _non_null_residues, \
                _residue_or_none, _residues_or_nones, _residues_or_nones, \
                _chains, _atomic_structure, _pseudobond_group, \
                _pseudobond_group_map

from numpy import int8, uint8, int32, uint32, float64, float32, byte, bool as npy_bool

from simtk import unit, openmm

from ..delayed_reaction import delayed_reaction

libdir = os.path.dirname(os.path.abspath(__file__))
libfile = glob.glob(os.path.join(libdir, '..', 'openmm.cpython*'))[0]

_c_functions = CFunctions(os.path.splitext(libfile)[0])
c_property = _c_functions.c_property
cvec_property = _c_functions.cvec_property
c_function = _c_functions.c_function
c_array_function = _c_functions.c_array_function

class OpenMM_Thread_Handler:
    '''
    A lightweight wrapper class for an OpenMM simulation `Context`, which pushes
    the task of stepping the simulation forward off to a C++ thread so that
    Python performance is not interrupted. Each call to `step()` or `minimize()`
    will create a short-lived thread to run the desired number of steps or
    a round of minimization, respectively. The status of the thread can be
    checked with thread_finished, while the initial and final coordinates can
    be retrieved with `last_coords` and `coords` respectively. If instability
    (overly fast-moving atoms) is detected, the thread will terminate early and
    the `unstable` property will set to True. In such cases it is advisable to
    run one or more minimization rounds. When minimization converges to with
    tolerance, `unstable` will be reset to False.

    Use with care! While nothing prevents you from using the standard single-
    thread OpenMM API alongside this one, it is up to you to ensure that
    no threads are running before making any calls that affect the simulation.
    Additionally, the `OpenMM_Thread_Handler` object *must* be destroyed before
    the Context it is attached to.
    '''
    def __init__(self, context, c_pointer=None):
        cname = type(self).__name__.lower()
        if c_pointer is None:
            new_func = cname + '_new'
            int_ptr = int(context.this)
            f = c_function(new_func, args=(ctypes.c_void_p,), ret=ctypes.c_void_p)
            c_pointer = f(int_ptr)
        set_c_pointer(self, c_pointer)
        f = c_function('set_'+cname+'_py_instance', args=(ctypes.c_void_p, ctypes.py_object))
        f(self._c_pointer, self)
        self.context = context

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def delete(self):
        c_function('openmm_thread_handler_delete', args=(ctypes.c_void_p,))(self._c_pointer)

    def step(self, steps):
        f = c_function('openmm_thread_handler_step',
            args=(ctypes.c_void_p, ctypes.c_size_t))
        f(self._c_pointer, steps)

    def minimize(self):
        '''
        Run an energy minimization on the coordinates. If the minimisation
        converges to within tolerance, unstable will be set to False.
        Don't forget to run reinitialize_velocities() before continuing
        equilibration!
        '''
        f = c_function('openmm_thread_handler_minimize',
            args = (ctypes.c_void_p,))
        f(self._c_pointer)

    def reinitialize_velocities(self):
        if not self.thread_finished:
            raise RuntimeError('Wait for the thread to finish before reinitialising velocities!')
        c = self.context
        i = c.getIntegrator()
        c.setVelocitiesToTemperature(i.getTemperature())

    def reinitialize_context_and_keep_state(self, threaded=True):
        '''
        Reinitialize the Context, keeping the current positions and velocities.
        A reinitialization is typically only required when the number of
        bonds/particles in a Force object changes.
        '''
        if threaded:
            func_name = 'openmm_thread_handler_reinitialize_context_and_keep_state_threaded'
        else:
            func_name = 'openmm_thread_handler_reinitialize_context_and_keep_state'
        f = c_function(func_name,
            args=(ctypes.c_void_p,))
        f(self._c_pointer)

    @property
    def natoms(self):
        f = c_function('openmm_thread_handler_num_atoms',
            args=(ctypes.c_void_p,),
            ret = ctypes.c_size_t)
        return f(self._c_pointer)

    def thread_finished(self):
        f = c_function('openmm_thread_handler_thread_finished',
            args=(ctypes.c_void_p,),
            ret=npy_bool)
        return f(self._c_pointer)

    def finalize_thread(self):
        f = c_function('openmm_thread_handler_finalize_thread',
            args=(ctypes.c_void_p,))
        f(self._c_pointer)

    def unstable(self):
        f = c_function('openmm_thread_handler_unstable',
            args=(ctypes.c_void_p,),
            ret = ctypes.c_bool)
        return f(self._c_pointer)

    @property
    def unstable_atoms(self):
        f = c_function('openmm_thread_handler_unstable_atoms',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p))
        n = self.natoms
        ret = numpy.zeros(n, numpy.bool)
        f(self._c_pointer, n, pointer(ret))
        return ret

    @property
    def last_coords(self):
        f = c_function('openmm_thread_handler_last_coords',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p))
        n = self.natoms
        coords = numpy.empty((n,3), float64)
        f(self._c_pointer, n, pointer(coords))
        return coords

    @property
    def coords(self):
        f = c_function('openmm_thread_handler_current_coords',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p))
        n = self.natoms
        coords = numpy.empty((n,3), float64)
        f(self._c_pointer, n, pointer(coords))
        return coords

    @coords.setter
    def coords(self, coords):
        f = c_function('set_openmm_thread_handler_current_coords',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p))
        n = self.natoms
        f(self._c_pointer, n, pointer(coords))
        self.reinitialize_velocities()

    def _get_min_thread_period(self):
        '''Throttle the simulation to a minimum time period per loop (in ms)'''
        f = c_function('openmm_thread_handler_min_thread_period',
            args=(ctypes.c_void_p,),
            ret=ctypes.c_double)
        return f(self._c_pointer)

    def _set_min_thread_period(self, period):
        f = c_function('set_openmm_thread_handler_min_thread_period',
            args=(ctypes.c_void_p, ctypes.c_double))
        f(self._c_pointer, period)

    min_thread_period = property(_get_min_thread_period, _set_min_thread_period)

class Sim_Construct:
    '''
    Container class defining the ChimeraX atoms (all, mobile and fixed)
    in a simulation.
    '''
    def __init__(self, all_atoms, mobile_atoms, fixed_atoms):
        self.all_atoms = all_atoms
        self.mobile_atoms = mobile_atoms
        self.fixed_atoms = fixed_atoms
        self.residue_templates = find_residue_templates(all_atoms.unique_residues)

class Sim_Handler:
    def __init__(self, session, sim_params, sim_construct, temperature):
        '''
        Prepares the simulation topology parameters and construct, and initialises
        the necessary Force objects to handle restraints. The restraint forces
        must be populated using e.g. add_dihedral_restraints() before initialising
        the context and beginning the simulation.
        '''
        self.session = session
        self._params = sim_params
        self._sim_construct = sim_construct

        self._thread_handler = None

        self._paused = False
        self._sim_running = False

        atoms = self._atoms = sim_construct.all_atoms
        # Forcefield used in this simulation
        from .forcefields import forcefields
        ff = self._forcefield = self.define_forcefield(forcefields[sim_params.forcefield])

        self.all_forces = []

        # Overall simulation topology
        top, residue_templates = self.create_openmm_topology(atoms, sim_construct.residue_templates)
        self._topology = top

        self.temperature = temperature

        system = self._system = self._create_openmm_system(ff, top,
            sim_params, residue_templates)
        self.set_fixed_atoms(sim_construct.fixed_atoms)
        self._thread_handler = None

        # List of IsoldeMap objects registered with this simulation
        self._maps = []
        # CustomExternalForce handling mouse and haptic interactions
        self._tugging_force = None

        self._force_update_pending = False
        self._context_reinit_pending = False

        trigger_names = (
            'coord update',
        )
        from chimerax.core.triggerset import TriggerSet
        t = self.triggers = TriggerSet()
        for name in trigger_names:
            t.add_trigger(name)


    @property
    def atoms(self):
        return self._atoms

    def _create_openmm_system(self, forcefield, top, params, residue_templates):
        system_params = {
            'nonbondedMethod':      params.nonbonded_cutoff_method,
            'nonbondedCutoff':      params.nonbonded_cutoff_distance,
            'constraints':          params.rigid_bonds,
            'rigidWater':           params.rigid_water,
            'removeCMMotion':       params.remove_c_of_m_motion,
            'ignoreExternalBonds':  True,
            'residueTemplates':     residue_templates,
        }
        return forcefield.createSystem(top, **system_params)

    def initialize_custom_forces(self, amber_cmap=True, tugging=True, position_restraints=True,
        distance_restraints=True, dihedral_restraints=True):
        '''
        Create the selected Force objects, and add them to the System. No
        restraints are added at this stage - these should be added using
        (e.g.) add_dihedral_restraints().
        '''
        params = self._params
        if amber_cmap:
            self.initialize_amber_cmap_force()
        if tugging:
            self.initialize_tugging_force(params.restraint_max_force)
        if position_restraints:
            self.initialize_position_restraints_force(params.restraint_max_force)
        if distance_restraints:
            self.initialize_distance_restraints_force(params.restraint_max_force)
        if dihedral_restraints:
            self.initialize_dihedral_restraint_force(params.dihedral_restraint_cutoff_angle)

    def _prepare_sim(self):
        params = self._params
        if params.use_gbsa:
            self.initialize_implicit_solvent(params)
        integrator = self._prepare_integrator(params)
        platform = openmm.Platform.getPlatformByName(params.platform)
        from simtk.openmm import app
        s = self._simulation = app.Simulation(self._topology, self._system,
            integrator, platform)
        c = self._context = s.context
        c.setPositions(0.1*self._atoms.coords)
        c.setVelocitiesToTemperature(self.temperature)
        self._thread_handler = OpenMM_Thread_Handler(c)

    def _prepare_integrator(self, params):
        integrator = params.integrator
        temperature = self.temperature
        if integrator == openmm.VariableLangevinIntegrator:
            integrator_params = (
                temperature,
                params.friction_coefficient,
                params.variable_integrator_tolerance,
                )
        elif integrator == openmm.LangevinIntegrator:
            integrator_params = (
                temperature,
                params.friction_coefficient,
                paramsfixed_integrator_timestep,
                )
        elif integrator == openmm.MTSIntegrator:
            raise RuntimeError('Multiple timestepping not yet implemented!')
        else:
            raise RuntimeError('Unrecognised or unsupported integrator: {}!'.format(integrator))
        _integrator = integrator(*integrator_params)
        _integrator.setConstraintTolerance(params.constraint_tolerance)
        return _integrator

    def start_sim(self):
        if self._sim_running:
            raise RuntimeError('Simulation is already running!')
        self._prepare_sim()
        self._pause = False
        self._stop = False
        self._sim_running = True
        self._minimize_and_go()

    def _minimize_and_go(self):
        th = self.thread_handler
        delayed_reaction(self.session.triggers, 'new frame', th.minimize, [],
            th.thread_finished, self._update_coordinates_and_repeat, [True])

    def _repeat_step(self):
        th = self.thread_handler
        params = self._params
        if self._context_reinit_pending:
            f = self._reinitialize_context
            f_args = []
            final_args = []
        elif th.unstable():
            f = th.minimize
            f_args = []
            final_args = [True]
        else:
            f = th.step
            f_args = (params.sim_steps_per_gui_update,)
            final_args = []
        delayed_reaction(self.session.triggers, 'new frame', f, f_args,
            th.thread_finished, self._update_coordinates_and_repeat, final_args)

    def _update_coordinates_and_repeat(self, reinit_vels = False):
        th = self.thread_handler
        self.atoms.coords = th.coords
        self.triggers.activate_trigger('coord update', None)
        if self._force_update_pending:
            self._update_forces_in_context_if_needed()
        if reinit_vels:
            th.reinitialize_velocities()
        if self._stop:
            self._thread_handler.delete()
            self._thread_handler = None
            self._simulation = None
            self._sim_running = False
            return
        if not self._pause:
            self._repeat_step()

    def push_coords_to_sim(self, coords=None):
        if not self._simulation_running:
            raise TypeError('No simulation running!')
        if not self.pause:
            raise TypeError('Simulation must be paused first!')
        if coords is None:
            coords = self._atoms.coords
        self.thread_handler.coords = coords

    @property
    def thread_handler(self):
        return self._thread_handler

    @property
    def pause(self):
        return self._pause

    @pause.setter
    def pause(self, flag):
        if not self._sim_running:
            raise TypeError('No simulation running!')
        if flag != self._pause:
            self._pause = flag
            if not flag:
                self._update_coordinates_and_repeat()

    def stop(self):
        self._stop = True
        if self.pause:
            self._thread_handler.delete()
            self._thread_handler = None
            self._simulation = None
            self._sim_running = False

    @property
    def sim_running(self):
        return self._sim_running

    def force_update_needed(self):
        if not self.sim_running:
            return
        if self._paused:
            self._update_forces_in_context_if_needed()
        else:
            self._force_update_pending = True


    def _update_forces_in_context_if_needed(self):
        context = self._context
        for f in self.all_forces:
            if f.update_needed:
                f.updateParametersInContext(context)
                f.update_needed = False
        self._force_update_pending = False

    def context_reinit_needed(self):
        if not self.sim_running:
            return
        if self._paused:
            self._reinitialize_context()
        else:
            self._context_reinit_pending = True

    def _reinitialize_context(self):
        th = self._thread_handler
        th.reinitialize_context_and_keep_state()
        self._context_reinit_pending = False

    ####
    # AMBER-specific CMAP backbone torsion corrections
    ####

    def initialize_amber_cmap_force(self):
        from .custom_forces import AmberCMAPForce
        cf = self._amber_cmap_force = AmberCMAPForce()
        self._system.addForce(cf)
        self.all_forces.append(cf)

    def add_amber_cmap_torsions(self, ramas):
        ''' Add CMAP correction terms for AMBER force field. '''
        cf = self._amber_cmap_force
        sc = self._atoms
        valid_ramas = ramas[ramas.valids]
        resnames = valid_ramas.residues.names
        phi_atoms = valid_ramas.phi_dihedrals.atoms
        psi_atoms = valid_ramas.psi_dihedrals.atoms
        phi_indices = numpy.column_stack([sc.indices(atoms) for atoms in phi_atoms])
        psi_indices = numpy.column_stack([sc.indices(atoms) for atoms in psi_atoms])
        cf.add_torsions(resnames, phi_indices, psi_indices)

    ####
    # Dihedral restraints
    ####

        ##
        # Before simulation starts
        ##

    def initialize_dihedral_restraint_force(self, default_cutoff):
        from .custom_forces import FlatBottomTorsionRestraintForce
        df = self._dihedral_restraint_force = FlatBottomTorsionRestraintForce()
        self._system.addForce(df)
        self.all_forces.append(df)

    def add_dihedral_restraints(self, restraints):
        force = self._dihedral_restraint_force
        all_atoms = self._atoms
        dihedral_atoms = restraints.dihedrals.atoms
        atom_indices = [all_atoms.indices(atoms) for atoms in dihedral_atoms]
        restraints.sim_indices = force.add_torsions(atom_indices,
            restraints.enableds, restraints.spring_constants, restraints.targets, restraints.cutoffs)
        self.context_reinit_needed()

    def add_dihedral_restraint(self, restraint):
        force = self._dihedral_restraint_force
        all_atoms = self._atoms
        dihedral_atoms = restraint.dihedral.atoms
        indices = [all_atoms.index(atom) for atom in dihedral_atoms]
        restraint.sim_index = force.add_torsion(*indices,
            (float(restraint.enabled), restraint.spring_constant, restraint.target, cos(restraint.cutoff)))
        self.context_reinit_needed()

        ##
        # During simulation
        ##

    def update_dihedral_restraints(self, restraints):
        force = self._dihedral_restraint_force
        force.update_targets(restraints.sim_indices,
            restraints.enableds, restraints.spring_constants, restraints.targets, restraints.cutoffs)
        self.force_update_needed()

    def update_dihedral_restraint(self, restraint):
        force = self._dihedral_restraint_force
        force.update_target(restraint.sim_index,
            enabled=restraint.enabled, k=restraint.spring_constant,
            target=restraint.target, cutoff=restraint.cutoff)
        self.force_update_needed()
    ####
    # Distance Restraints
    ####

        ##
        # Before simulation starts
        ##

    def initialize_distance_restraints_force(self, max_force):
        from .custom_forces import TopOutBondForce
        tf = self._distance_restraints_force = TopOutBondForce(max_force)
        self._system.addForce(tf)
        self.all_forces.append(tf)
        return tf

    def add_distance_restraints(self, restraints):
        force = self._distance_restraints_force
        all_atoms = self._atoms
        dr_atoms = restraints.atoms
        indices = [all_atoms.indices(atoms) for atoms in dr_atoms]
        restraints.sim_indices = force.add_bonds(indices,
            restraints.enableds, restraints.spring_constants, restraints.targets/10)
        self.context_reinit_needed()

    def add_dihedral_restraint(self, restraint):
        force = self._distance_restraints_force
        all_atoms = self._atoms
        dr_atoms = restraint.atoms
        indices = [all_atoms.index(atom) for atom in dr_atoms]
        restraint.sim_index = force.addBond(*indices,
            (float(restraint.enabled), restraint.spring_constant, restraint.target/10))
        self.context_reinit_needed()

        ##
        # During simulation
        ##

    def update_distance_restraints(self, restraints):
        force = self._distance_restraints_force
        force.update_targets(restraints.sim_indices,
            restraints.enableds, restraints.spring_constants, restraints.targets/10)
        self.force_update_needed()

    def update_distance_restraint(self, restraint):
        force = self._distance_restraints_force
        force.update_target(restraint.sim_index,
            restraint.enabled, k=restraint.spring_constant, target=restraint.target/10)
        self.force_update_needed()

    ####
    # Positional Restraints
    ####

        ##
        # Before simulation starts
        ##

    def initialize_position_restraints_force(self, max_force):
        from .custom_forces import TopOutRestraintForce
        rf = self._position_restraints_force = TopOutRestraintForce(max_force)
        self._system.addForce(rf)
        self.all_forces.append(rf)
        return rf

    def add_position_restraints(self, restraints):
        force = self._position_restraints_force
        all_atoms = self._atoms
        indices = all_atoms.indices(restraints.atoms)
        restraints.sim_indices = force.add_particles(indices,
            restraints.enableds, restraints.spring_constants, restraints.targets/10)
        self.context_reinit_needed()

    def add_position_restraint(self, restraint):
        force = self._position_restraints_force
        index = self._all_atoms.index(restraint.atom)
        target = (restraint.target/10).tolist()
        restraint.sim_index = force.addParticle(index,
            (restraint.enabled, restraint.spring_constant, *target))
        self.context_reinit_needed()

        ##
        # During simulation
        ##

    def update_position_restraints(self, restraints):
        force = self._position_restraints_force
        force.update_targets(restraints.sim_indices,
            restraints.enableds, restraints.spring_constants, restraints.targets/10)
        self.force_update_needed()

    def update_position_restraint(self, restraint):
        force = self._position_restraints_force
        force.update_target(restraint.sim_index,
            restraint.enabled, restraint.spring_constant, restraint.target/10)
        self.force_update_needed()

    ####
    # Tugging force
    ####

        ##
        # Before simulation starts
        ##

    def initialize_tugging_force(self, max_force):
        from .custom_forces import TopOutRestraintForce
        f = self._tugging_force = TopOutRestraintForce(max_force)
        self._system.addForce(f)
        self.all_forces.append(f)

    def add_tuggables(self, tuggables):
        force = self._tugging_force
        all_atoms = self._atoms
        indices = all_atoms.indices(tuggables.atoms)
        tuggables.sim_indices = force.add_particles(indices,
            tuggables.enableds, tuggables.spring_constants, tuggables.targets/10)
        self.context_reinit_needed()

    def add_tuggable(self, tuggable):
        force = self._tugging_force
        index = self._all_atoms.index(tuggable.atom)
        target = (tuggable.target/10).tolist()
        tuggable.sim_index = force.addParticle(index,
            (float(tuggable.enabled), tuggable.spring_constant, *target))
        self.context_reinit_needed()

        ##
        # During simulation
        ##

    def update_tuggable(self, tuggable):
        force = self._tugging_force
        force.update_target(tuggable.sim_index,
            tuggable.enabled, tuggable.spring_constant, tuggable.target/10)
        self.force_update_needed()

    def update_tuggables(self, tuggables):
        force = self._tugging_force
        force.update_targets(tuggables.sim_indices,
            tuggables.enableds, tuggables.spring_constants, tuggables.targets/10)
        self.force_update_needed()


    def set_fixed_atoms(self, fixed_atoms):
        '''
        Fix the desired atoms rigidly in space. NOTE: a fixed atom can not be
        bonded to a mobile atom via a rigid bond. In most cases this means that
        you cannot fix a hydrogen without fixing the heavy atom that it's bonded
        to, and any fixed heavy atom must have all its hydrogens fixed.
        While atoms may be fixed and un-fixed during a simulation, this requires
        a costly re-initialisation of the simulation context. In most cases it's
        best to simply use strong position restraints instead.
        '''
        fixed_indices = self._atoms.indices(fixed_atoms).tolist()
        sys = self._system
        for index in fixed_indices:
            sys.setParticleMass(index, 0)
        self.context_reinit_needed()

    def release_fixed_atoms(self, atoms):
        '''
        Make the desired fixed atoms mobile again. NOTE: a fixed atom can not be
        bonded to a mobile atom via a rigid bond. In most cases this means that
        you cannot fix a hydrogen without fixing the heavy atom that it's bonded
        to, and any fixed heavy atom must have all its hydrogens fixed.
        While atoms may be fixed and un-fixed during a simulation, this requires
        a costly re-initialisation of the simulation context. In most cases it's
        best to simply use strong position restraints instead.
        '''
        indices = self._atoms.indices(fixed_atoms).tolist()
        masses = atoms.elements.masses
        sys = self.system
        for index, mass in zip(indices, masses):
            sys.setParticleMass(index, mass)
        self.context_reinit_needed()




    def define_forcefield(self, forcefield_file_list):
        from simtk.openmm.app import ForceField
        ff = ForceField(*[f for f in forcefield_file_list if f is not None])
        return ff

    def create_openmm_topology(self, atoms, residue_templates):
        '''
        Generate a simulation topology from a set of atoms.
        @param atoms:
            A ChimeraX Atoms object. Residues must be complete, and the atoms
            in each residue should be contiguous in the array.
        @param residue_templates:
            A {residue_index: residue_type} dict for residues whose
            topology is otherwise ambiguous. OpenMM requires a
            {openmm_residue_object: residue_type} dict, so we need to
            do the conversion here.
        '''

        anames   = atoms.names
        n = len(anames)
        enames   = atoms.element_names
        rnames   = atoms.residues.names
        rnums    = atoms.residues.numbers
        cids    = atoms.chain_ids
        bonds = atoms.intra_bonds
        bond_is = [atoms.indices(alist) for alist in bonds.atoms]

        #template_indices = list(residue_templates.keys())
        templates_out = {}
        from simtk.openmm.app import Topology, Element
        top = self.topology = Topology()
        cmap = {}
        rmap = {}
        atoms = {}
        rcount = 0
        for i, (aname, ename, rname, rnum, cid) in enumerate(
                zip(anames, enames, rnames, rnums, cids)):
            if not cid in cmap:
                cmap[cid] = top.addChain()   # OpenMM chains have no name
            rid = (rname, rnum, cid)
            if not rid in rmap:
                res = rmap[rid] = top.addResidue(rname, cmap[cid])
                if rcount in residue_templates.keys():
                    templates_out[res] = residue_templates[rcount]
                rcount += 1


            element = Element.getBySymbol(ename)
            atoms[i] = top.addAtom(aname, element,rmap[rid])

        for i1, i2 in zip(*bond_is):
            top.addBond(atoms[i1],  atoms[i2])

        return top, templates_out

    def initialize_implicit_solvent(self, params):
        '''Add a Generalised Born Implicit Solvent (GBIS) formulation.'''
        # Somewhat annoyingly, OpenMM doesn't store atomic charges in a
        # nice accessible format. So, we have to pull it back out of the
        # NonbondedForce term.

        gbsa_params = {
            'solventDielectric':    params.gbsa_solvent_dielectric,
            'soluteDielectric':     params.gbsa_solute_dielectric,
            'SA':                   params.gbsa_sa_method,
            'cutoff':               params.nonbonded_cutoff_distance,
            'kappa':                params.gbsa_kappa,
            'nonbonded_method':     params.gbsa_cutoff_method,
            }

        top = self.topology
        system = self._system
        from simtk.openmm.openmm import NonbondedForce
        for f in system.getForces():
            if isinstance(f, NonbondedForce):
                break
        charges = []
        for i in range(f.getNumParticles()):
            charges.append(f.getParticleParameters(i)[0])
            from .custom_forces import GBSAForce
        gbforce = self._gbsa_force = GBSAForce(**gbsa_params)
        params = gbforce.getStandardParameters(top)
        for charge, param in zip(charges, params):
            gbforce.addParticle([charge, *param])
        gbforce.finalize()
        system.addForce(gbforce)
        #self.all_forces.append(gbforce)


def find_residue_templates(residues):
    templates = {}
    cys_indices = numpy.where(residues.names == 'CYS')[0]
    for c_i in cys_indices:
        templates[c_i] = cys_type(residues[c_i])
    return templates

def cys_type(residue):
    from chimerax.core.atomic import Bonds, concatenate
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
