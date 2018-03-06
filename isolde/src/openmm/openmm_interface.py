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


class Sim_Handler:
    def __init__(self):
        # Forcefield used in this simulation
        self._forcefield = None
        # Overall simulation topology
        self._topology = None
        # Atoms in simulation topology
        self._atoms = None

        self.all_forces = []

        # Dict holding each custom external force object, its global and its
        # per-particle parameters. Format:
        # {'name': [object,
        #           {'global_name': value},
        #           [per_particle_names],
        #           [per_particle_default_values],
        #           {index_in_topology: index_in_force}
        # }
        self._custom_external_forces = {}

        # List of IsoldeMap objects registered with this simulation
        self._maps = []
        # CustomExternalForce handling mouse and haptic interactions
        self._tugging_force = None

        self._force_update_pending = False

        trigger_names = (
            'coord update',
        )
        from chimerax.core.triggerset import TriggerSet
        t = self.triggers = TriggerSet()
        for name in trigger_names:
            t.add_trigger(name)

        self._thread_handler = None

    def force_update_needed(self):
        if not self._force_update_pending:
            from ..delayed_reaction import delayed_reaction
            delayed_reaction(self.triggers, 'coord update', [], self._thread_handler.thread_finished,
                self._update_forces_in_context_if_needed, [])

    def _update_forces_in_context_if_needed(self):
        context = self.context
        for f in self.all_forces:
            if f.update_needed:
                f.updateParametersInContext(context)
                f.update_needed = False



    # AMBER-specific CMAP backbone torsion corrections

    def initialize_amber_cmap_force(self):
        from .custom_forces import AmberCMAPForce
        cf = self._amber_cmap_force = AmberCMAPForce()
        self.all_forces.append(cf)

    def add_amber_cmap_torsions(self, ramas, sim_construct):
        ''' Add CMAP correction terms for AMBER force field. '''
        cf = self._amber_cmap_force
        sc = sim_construct
        valid_ramas = ramas[ramas.valids]
        resnames = valid_ramas.residues.names
        phi_atoms = valid_ramas.phi_dihedrals.atoms
        psi_atoms = valid_ramas.psi_dihedrals.atoms
        phi_indices = numpy.column_stack([sc.indices(atoms) for atoms in phi_atoms])
        psi_indices = numpy.column_stack([sc.indices(atoms) for atoms in psi_atoms])
        for resname, pi, fi in zip(resnames, phi_indices, psi_indices):
            cf.addTorsion(resname, pi, fi)


    ####
    # Dihedral restraints
    ####

        ##
        # Before simulation starts
        ##

    def initialize_dihedral_restraint_force(self, default_cutoff):
        from .custom_forces import FlatBottomTorsionRestraintForce
        df = self._dihedral_restraint_force = FlatBottomTorsionRestraintForce()
        self.default_torsion_cutoff = default_cutoff
        self.all_forces.append(df)

    def add_dihedral_restraints(self, restraints):
        '''
        Add a set of dihedral restraints to the simulation
        '''
        force = self._dihedral_restraint_force
        all_atoms = self._atoms
        dihedral_atoms = restraints.dihedrals.atoms
        atom_indices = [all_atoms.indices(atoms) for atoms in dihedral_atoms]
        restraints.sim_indices = force.add_torsions(atom_indices,
            restraints.targets, restraints.spring_constants, restraints.cutoffs)

    def add_dihedral_restraint(self, restraint):
        force = self._dihedral_restraint_force
        all_atoms = self._atoms
        dihedral_atoms = restraint.dihedral.atoms
        indices = [all_atoms.index(atom) for atom in dihedral_atoms]
        restraint.sim_index = force.add_torsion(*indices,
            (restraint.target, restraint.spring_constant, cos(restraint.cutoff)))

        ##
        # During simulation
        ##

    def update_dihedral_restraints_in_context(self, context):
        rf = self._dihedral_restraint_force
        if rf.update_needed:
            rf.updateParametersInContext(context)
        rf.update_needed = False

    def update_dihedral_restraints(self, restraints):
        force = self._dihedral_restraint_force
        force.update_targets(restraints.sim_indices, restraints.targets,
            restraints.spring_constants, restraints.cutoffs)

    def update_dihedral_restraint(self, restraint):
        force = self._dihedral_restraint_force
        force.update_target(restraint.sim_index, target=restraint.target,
            k=restraint.spring_constant, cutoff=restraint.cutoff)

    ####
    # Distance Restraints
    ####

        ##
        # Before simulation starts
        ##

    def initialize_distance_restraints_force(self, max_force):
        from .custom_forces import TopOutBondForce
        tf = self._distance_restraints_force = TopOutBondForce(max_force)
        self.all_forces.append(tf)
        return tf

    def add_distance_restraints(self, restraints):
        '''
        Add a set of distance restraints to the simulation
        '''
        force = self._distance_restraints_force
        all_atoms = self._atoms
        dr_atoms = restraints.atoms
        indices = [all_atoms.indices(atoms) for atoms in dr_atoms]
        restraints.sim_indices = force.add_bonds(indices, restraints.targets/10, restraints.spring_constants)

    def add_dihedral_restraint(self, restraint):
        force = self._distance_restraints_force
        all_atoms = self._atoms
        dr_atoms = restraint.atoms
        indices = [all_atoms.index(atom) for atom in dr_atoms]
        restraint.sim_index = force.addBond(*indices, (restraint.spring_constant, restraint.target/10))

        ##
        # During simulation
        ##

    def update_distance_restraints(self, restraints):
        force = self._distance_restraints_force
        force.update_targets(restraints.sim_indices, restraints.targets/10, restraints.spring_constants)

    def update_distance_restraint(self, restraint):
        force = self._distance_restraints_force
        force.update_target(restraint.sim_index, target=restraint.target/10, k=restraint.spring_constant)

    ####
    # Positional Restraints
    ####

        ##
        # Before simulation starts
        ##

    def initialize_position_restraints_force(self, max_force):
        from .custom_forces import TopOutRestraintForce
        rf = self._position_restraints_force = TopOutRestraintForce(max_force)
        self.all_forces.append(rf)
        return rf

    def add_position_restraints(self, restraints):
        force = self._position_restraints_force
        all_atoms = self._atoms
        indices = all_atoms.indices(restraints.atoms)
        restraints.sim_indices = force.add_particles(indices, restraints.targets/10, restraints.spring_constants)

    def add_position_restraint(self, restraint):
        force = self._position_restraints_force
        index = self._all_atoms.index(restraint.atom)
        target = (restraint.target/10).tolist()
        restraint.sim_index = force.addParticle(index, (restraint.spring_constant, *target))

        ##
        # During simulation
        ##

    def update_position_restraints(self, restraints):
        force = self._position_restraints_force
        force.update_targets(restraints.indices, restraints.targets/10, restraints.spring_constants)

    def update_position_restraint(self, restraint):
        force = self._position_restraint_force
        force.update_target(restraint.index, restraint.target/10, restraint.spring_constant)
    
    ####
    # Tugging force
    ####

        ##
        # Before simulation starts
        ##

    def initialize_tugging_force(self, max_force):
        from .custom_forces import TopOutRestraintForce
        f = self._tugging_force = TopOutRestraintForce(max_force)
        self.all_forces.append(f)

    def couple_atoms_to_tugging_force(self, indices, force_map):
        force = self._tugging_force
        for i, index in enumerate(indices):
            force_map[i] = force.addParticle(int(index), (0,0,0,0))

        ##
        # During simulation
        ##

    def tug_atom(self, force_index, target, k):
        tf = self._tugging_force
        tf.update_target(force_index, target, k)

    def tug_atoms(self, force_indices, targets, ks):
        tf = self._tugging_force
        for i, t, k in zip(force_indices, targets, ks):
            tf.update_target(i, t, k)
