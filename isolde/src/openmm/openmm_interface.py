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

from ..constants import defaults

class OpenMM_Thread_Handler:
    '''
    A lightweight wrapper class for a :class:`openmm.Context`, which
    pushes time-consuming tasks off to a C++ thread so that Python performance
    is not interrupted. Each call to :func:`step` or :func:`minimize` will
    create a short-lived thread to run the desired number of steps or a round of
    minimization, respectively. Where necessary (e.g. where new  restraints are
    added to the simulation), the context may be reinitialised  with
    :func:`reinitialize_context_and_keep_state`. The status of the thread can be
    checked with :attr:`thread_finished`, while the initial and final
    coordinates can be retrieved with :attr:`last_coords` and :attr:`coords`
    respectively. Within the thread, the simulation is checked for excessive
    velocities every 10 steps. If instability (overly fast-moving atoms) is
    detected, the thread will terminate early and :attr:`unstable` will be set
    to True. In such cases it is advisable to run one or more minimization
    rounds. When minimization converges to within tolerance, `unstable` will be
    reset to False.

    Use with care! While nothing prevents you from using the standard single-
    thread OpenMM API alongside this one, it is up to you to ensure that no
    threads are running before making any calls that affect the simulation.
    Additionally, the :class:`OpenMM_Thread_Handler` object *must* be destroyed
    before the :class:`openmm.Context` it is attached to.
    '''
    def __init__(self, context, c_pointer=None):
        '''
        Initialise the thread handler.

        Args:
            * context:
                - a :py:class:`openmm.Context` instance
        '''
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
        '''
        Value that can be passed to C++ layer to be used as pointer (Python int)
        '''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def delete(self):
        c_function('openmm_thread_handler_delete', args=(ctypes.c_void_p,))(self._c_pointer)

    def step(self, steps):
        '''
        Advance the simulation integrator by the desired number of steps.

        Args:
            * steps:
                - an integer value
        '''
        f = c_function('openmm_thread_handler_step',
            args=(ctypes.c_void_p, ctypes.c_size_t))
        f(self._c_pointer, steps)

    def minimize(self):
        '''
        Run an energy minimization on the coordinates. If the minimisation
        converges to within tolerance, unstable will be set to False.
        Don't forget to run :func:`reinitialize_velocities` before continuing
        equilibration!
        '''
        f = c_function('openmm_thread_handler_minimize',
            args = (ctypes.c_void_p,))
        f(self._c_pointer)

    def reinitialize_velocities(self):
        '''
        Set the atomic velocities to random values consistent with the current
        temperature. Recommended after any energy minimisation.
        '''
        self.finalize_thread()
        c = self.context
        i = c.getIntegrator()
        c.setVelocitiesToTemperature(i.getTemperature())

    def reinitialize_context_and_keep_state(self):
        '''
        Reinitialize the Context, keeping the current positions and velocities.
        A reinitialization is typically only required when the number of
        bonds/particles in a Force object changes.
        '''
        f = c_function(
            'openmm_thread_handler_reinitialize_context_and_keep_state_threaded',
            args=(ctypes.c_void_p,))
        f(self._c_pointer)

    @property
    def natoms(self):
        '''Number of atoms in the simulation.'''
        f = c_function('openmm_thread_handler_num_atoms',
            args=(ctypes.c_void_p,),
            ret = ctypes.c_size_t)
        return f(self._c_pointer)

    def thread_finished(self):
        '''Has the current thread finished its task?'''
        f = c_function('openmm_thread_handler_thread_finished',
            args=(ctypes.c_void_p,),
            ret=npy_bool)
        return f(self._c_pointer)

    def finalize_thread(self):
        '''
        Wrap up and join the existing thread. Note that if the thread has not
        finished, Python and GUI will hang until it has.
        '''
        f = c_function('openmm_thread_handler_finalize_thread',
            args=(ctypes.c_void_p,))
        f(self._c_pointer)

    def unstable(self):
        '''Returns true if any atoms in the simulation are moving too fast.'''
        f = c_function('openmm_thread_handler_unstable',
            args=(ctypes.c_void_p,),
            ret = ctypes.c_bool)
        return f(self._c_pointer)

    @property
    def unstable_atoms(self):
        '''
        Returns a Boolean mask indicating which atoms are currently moving too
        fast.
        '''
        f = c_function('openmm_thread_handler_unstable_atoms',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p))
        n = self.natoms
        ret = numpy.empty(n, numpy.bool)
        f(self._c_pointer, n, pointer(ret))
        return ret

    @property
    def last_coords(self):
        '''
        Returns the coordinates of the atoms as they were at the start of the
        most recent thread. Switching between these and the final coords can
        be useful for diagnosing instability.
        '''
        f = c_function('openmm_thread_handler_last_coords',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p))
        n = self.natoms
        coords = numpy.empty((n,3), float64)
        f(self._c_pointer, n, pointer(coords))
        return coords

    @property
    def coords(self):
        '''
        Returns the coordinates of the atoms after the most recent thread
        completes. Can also be set, to push edited coordinates back to the
        simulation.
        '''
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
        '''
        Throttle the simulation to a minimum time period per loop (in ms).
        Useful when graphics performance needs to be prioritised over simulation
        performance. The default value is 1.0, which in almost all situations
        means the simulation will not be throttled.
        '''
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
    in a simulation. Also responsible for storing the visualisation state of
    these atoms prior to simulation startup, and reverting it when done.
    '''
    def __init__(self, model, mobile_atoms, fixed_atoms):
        '''
        Prepare the construct. The atoms in each array will be sorted in the
        same order as :attr:`model.residues.atoms`, primarily because OpenMM
        requires atoms to be grouped by residue.

        Args:
            * model:
                - The :py:class:`chimerax.AtomicStructure` containing the atoms
                  in this simulation
            * mobile_atoms:
                - A :py:class:`chimerax.Atoms` instance defining the atoms that
                  are to be mobile
            * fixed_atoms:
                - A :py:class:`chimerax.Atoms` instance defining the atoms that
                  are to be fixed

        NOTE: the simulation will fail if mobile_atoms and fixed_atoms do not
        combine to form a set containing only complete residues. Also note that
        OpenMM does not allow a fixed atom to be rigidly bonded to a mobile
        one (in practice this means any fixed heavy atom must have its hydrogens
        also fixed).
        '''
        self.model = model

        # Sort all the atoms according to their order in the model#
        model_atoms = model.residues.atoms
        if len(mobile_atoms.intersect(fixed_atoms)):
            raise TypeError('Atoms cannot be both fixed and mobile!')
        from chimerax.atomic import concatenate
        all_atoms = concatenate((mobile_atoms, fixed_atoms))
        all_i = model_atoms.indices(all_atoms)
        if -1 in all_i:
            raise TypeError('All atoms must be from the targeted model!')
        all_atoms = self._all_atoms = model_atoms[numpy.sort(all_i)]
        mob_i = model_atoms.indices(mobile_atoms)
        ma = self._mobile_atoms = model_atoms[numpy.sort(mob_i)]
        mr = self._mobile_residues = ma.unique_residues
        fixed_i = model_atoms.indices(fixed_atoms)
        self._fixed_atoms = model_atoms[numpy.sort(fixed_i)]

        self.surroundings = model_atoms.subtract(all_atoms)
        self.residue_templates = find_residue_templates(all_atoms.unique_residues)
        self.store_original_visualisation()

    @property
    def mobile_atoms(self):
        '''
        A :py:class:`chimerax.Atoms` instance containing the mobile selection,
        sorted in the same order as :attr:`model.residues.atoms`.
        '''
        return self._mobile_atoms

    @property
    def mobile_residues(self):
        '''
        A :py:class:`chimerax.Residues` instance containing the mobile residues.
        '''
        return self._mobile_residues

    @property
    def fixed_atoms(self):
        '''
        A :py:class:`chimerax.Atoms` instance containing the fixed selection,
        sorted in the same order as :attr:`model.residues.atoms`.
        '''
        return self._fixed_atoms

    @property
    def all_atoms(self):
        '''
        A :py:class:`chimerax.Atoms` instance containing all atoms in the
        simulation, sorted in the same order as :attr:`model.residues.atoms`.
        '''
        return self._all_atoms

    def store_original_visualisation(self):
        '''
        Store the current visualisation state of the model, so it can be
        reverted once the simulation is done.
        '''
        m = self.model
        atoms = m.atoms
        bonds = m.bonds
        residues = m.residues
        self.spotlight_mode = None
        from chimerax.clipper.symmetry import get_symmetry_handler
        sym = self.symmetry_handler = get_symmetry_handler(m)
        if sym:
            self.spotlight_mode = sym.spotlight_mode
        self._original_atom_states = (
            atoms.colors,
            atoms.draw_modes,
            atoms.displays,
            atoms.radii
            )
        self._original_bond_states = (
            bonds.displays,
            bonds.radii,
        )
        self._original_residue_states = (
            residues.ribbon_displays,
            residues.ribbon_colors
        )


    def revert_visualisation(self):
        '''
        Return the model visualisation to the way it was before the simulation
        began.
        '''
        m = self.model
        atoms = m.atoms
        bonds = m.bonds
        residues = m.residues
        atoms.colors, atoms.draw_modes, atoms.displays, atoms.radii = \
            self._original_atom_states
        bonds.displays, bonds.radii = \
            self._original_bond_states
        residues.ribbon_displays, residues.ribbon_colors = \
            self._original_residue_states

        sym = self.symmetry_handler
        if sym:
            sym.spotlight_mode = self.spotlight_mode



class Sim_Manager:
    '''
    Responsible for creating the :py:class:`Sim_Handler` and managing the
    high-level control of the simulation. Handles all the necessary callbacks
    for automatically updating restraints in the simulation whenever their
    parameters change.
    '''
    def __init__(self, isolde, model, selected_atoms,
        isolde_params, sim_params, expansion_mode = 'extend'):
        '''
        Prepares a simulation according to the following workflow:
            * Expands an initial selection of atoms to complete residues
              according to the rules defined by expansion_mode
            * Finds a shell of residues around this selection to act as the
              fixed context
            * Restricts the live validation managers to focus only on the
              mobile selection (creating the managers as necessary)
            * Finds/creates all restraint managers for the model
            * Expands the fixed atom selection to include any non-mobile
              residues containing atoms participating in distance restraints
              with mobile atoms
            * Creates the :py:class:`Sim_Construct` object
            * Prepares the molecule visualisation for simulation (masking maps
              to the mobile selection, hiding atoms not in the simulation, etc.)
            * Prepares the MDFF managers (NOTE: this *must* be done after the
              preceding step, which ensures that each :py:class:`chimerax.Volume`
              has a region covering the mobile selection with sufficient padding)
            * Prepares all necessary callbacks to update the simulation when
              the parameters of restraints, mdff atom proxies etc. change.
            * Creates the :py:class:`Sim_Handler`
            * Adds all existing restraints and MDFF atom proxies to the
              simulation.

        Args:
            * isolde:
                - the current isolde session
            * model:
                - the :py:class:`chimerax.AtomicStructure` from which the
                  simulation atoms are drawn
            * selected_atoms:
                - a :py:class:`chimerax.Atoms` instance defining the initial
                    selection around which the simulation will be built.
            * isolde_params:
                - a :py:class:`IsoldeParams` instance
            * sim_params:
                - a :py:class:`SimParams` instance
            * expansion_mode:
                - string defining how the initial selection will be expanded.
                  For allowable options, see :func:`expand_mobile_selection`
        '''
        self.isolde = isolde
        self.model = model
        session = self.session = model.session
        self.isolde_params = isolde_params
        self.sim_params = sim_params
        self._revert_to = None
        mobile_atoms = self.expand_mobile_selection(selected_atoms, expansion_mode)
        from ..selections import get_shell_of_residues
        fixed_atoms = get_shell_of_residues(mobile_atoms.unique_residues,
            isolde_params.hard_shell_cutoff_distance).atoms
        self._prepare_validation_managers(mobile_atoms)
        self._prepare_restraint_managers()
        fixed_atoms = self._add_fixed_atoms_from_distance_restraints(mobile_atoms, fixed_atoms)
        sc = self.sim_construct = Sim_Construct(model, mobile_atoms, fixed_atoms)
        self.prepare_sim_visualisation()

        self._prepare_mdff_managers()
        sh = self.sim_handler = Sim_Handler(session, sim_params, sc)
        uh = self._update_handlers = []
        self._initialize_restraints(uh)
        self._initialize_mdff(uh)

    def start_sim(self):
        '''
        Start the simulation running. Should only be run once.
        '''
        sh = self.sim_handler
        sh.start_sim()
        from ..checkpoint import CheckPoint
        self._starting_checkpoint = self._current_checkpoint = CheckPoint(self.isolde)
        sh.triggers.add_handler('sim terminated', self._sim_end_cb)

    def stop_sim(self, revert = None):
        '''
        Stop the simulation and optionally revert to a previous state.

        Args:
            * revert:
                - None: keep the current coordinates
                - 'checkpoint': revert to the last saved checkpoint
                - 'start': revert to the state the model was in prior to
                  starting the simulation
        '''
        self.sim_handler.stop()
        self._revert_to = revert

    def toggle_pause(self):
        '''
        Pause/resume the simulation.
        '''
        self.sim_handler.toggle_pause()

    def checkpoint(self):
        '''
        A :py:class:`CheckPoint` is a snapshot of the current simulation state
        (including the state of all restraints) that can be returned to at any
        time, discarding any changes since the checkpoint was saved.

        :py:class:`Sim_Manager` automatically saves a checkpoint when the
        simulation is started. This method allows the saving of an intermediate
        state - particularly useful before experimenting on an ambiguous
        region of your map. Each call to :func:`checkpoint` overwrites the
        previously-saved one. When ending a simulation you may choose to keep
        the final coordinates or revert to either of the last saved checkpoint
        or the start-of-simulation checkpoint.
        '''
        from ..checkpoint import CheckPoint
        self._current_checkpoint = CheckPoint(self.isolde)

    def revert_to_checkpoint(self):
        '''
        Reverts to the last saved checkpoint. If no checkpoint has been manually
        saved, reverts to the start of the simulation.
        '''
        self._current_checkpoint.revert()

    def _prepare_validation_managers(self, mobile_atoms):
        from .. import session_extensions as sx
        m = self.model
        mobile_res = mobile_atoms.unique_residues
        rama_a = self.rama_annotator = sx.get_rama_annotator(m)
        rota_a = self.rota_annotator = sx.get_rota_annotator(m)
        rama_a.restrict_to_selected_residues(mobile_res)
        rota_a.restrict_to_selected_residues(mobile_res)

    def _prepare_restraint_managers(self):
        from .. import session_extensions as sx
        m = self.model
        self.proper_dihedral_restraint_mgr = sx.get_proper_dihedral_restraint_mgr(m)
        self.position_restraint_mgr = sx.get_position_restraint_mgr(m)
        self.tuggable_atoms_mgr = sx.get_tuggable_atoms_mgr(m)
        self.distance_restraint_mgr = sx.get_distance_restraint_mgr(m)

    def _initialize_restraints(self, update_handlers):
        sh = self.sim_handler
        sc = self.sim_construct
        sim_params = self.sim_params
        uh = update_handlers
        mobile_res = sc.mobile_atoms.unique_residues
        sh.initialize_restraint_forces()
        from .. import session_extensions as sx
        rama_mgr = sx.get_ramachandran_mgr(self.session)
        ramas = rama_mgr.get_ramas(mobile_res)
        ramas = ramas[ramas.valids]
        sh.add_amber_cmap_torsions(ramas)

        pdr_m = self.proper_dihedral_restraint_mgr
        pdrs = pdr_m.add_all_defined_restraints_for_residues(mobile_res)
        if sim_params.restrain_peptide_omegas:
            import numpy
            omega_rs = pdr_m.get_restraints_by_residues_and_name(mobile_res, 'omega')
            omega_angles = omega_rs.dihedrals.angles
            from math import pi
            # Restrain all bonds > cutoff to trans, otherwise cis
            mask = numpy.abs(omega_angles) > \
                sim_params.cis_peptide_bond_cutoff_angle.value_in_unit(unit.radians)
            omega_rs.targets = numpy.ones(len(mask))*pi*mask
            omega_rs.cutoffs = sim_params.dihedral_restraint_cutoff_angle.value_in_unit(unit.radians)
            omega_rs.displays = sim_params.display_omega_restraints
            omega_rs.spring_constants = sim_params.peptide_bond_spring_constant.value_in_unit(unit.kilojoule_per_mole/unit.radians**2)
            omega_rs.enableds = True


        sh.add_dihedral_restraints(pdrs)
        uh.append((pdr_m, pdr_m.triggers.add_handler('changes', self._pdr_changed_cb)))

        dr_m = self.distance_restraint_mgr
        # Pre-create all restraints necessary for secondary structure manipulation
        dr_m.add_ss_restraints(sc.all_atoms.unique_residues)
        drs = dr_m.atoms_restraints(sc.mobile_atoms)
        sh.add_distance_restraints(drs)
        uh.append((dr_m, dr_m.triggers.add_handler('changes', self._dr_changed_cb)))

        pr_m = self.position_restraint_mgr
        prs = pr_m.add_restraints(sc.mobile_atoms)
        sh.add_position_restraints(prs)
        uh.append((pr_m, pr_m.triggers.add_handler('changes', self._pr_changed_cb)))

        ta_m = self.tuggable_atoms_mgr
        tuggables = ta_m.add_tuggables(sc.mobile_atoms)
        uh.append((ta_m, ta_m.triggers.add_handler('changes', self._tug_changed_cb)))
        sh.add_tuggables(tuggables)

    def _prepare_mdff_managers(self):
        from .. import session_extensions as sx
        m = self.model
        mdff_mgr_map = self.mdff_mgrs = {}
        from chimerax.map import Volume
        for v in m.all_models():
            if isinstance(v, Volume):
                mdff_mgr_map[v] = sx.get_mdff_mgr(m, v)
        if len(mdff_mgr_map.keys()):
            isolde_params = self.isolde_params
            from chimerax.clipper.symmetry import get_symmetry_handler
            sym = get_symmetry_handler(m)
            sym.isolate_and_cover_selection(self.sim_construct.mobile_atoms,
                include_surrounding_residues = 0,
                show_context = isolde_params.hard_shell_cutoff_distance,
                mask_radius = isolde_params.standard_map_mask_cutoff,
                extra_padding = 5,
                hide_surrounds = isolde_params.hide_surroundings_during_sim,
                focus = False,
                include_symmetry = True)


    def _initialize_mdff(self, update_handlers):
        sh = self.sim_handler
        sc = self.sim_construct
        sp = self.sim_params
        uh = update_handlers
        mdff_mgrs = self.mdff_mgrs
        sh.initialize_mdff_forces(list(mdff_mgrs.keys()))
        for v, mgr in mdff_mgrs.items():
            mdff_atoms = mgr.add_mdff_atoms(sc.mobile_atoms,
                hydrogens = sp.hydrogens_feel_maps)
            sh.add_mdff_atoms(mdff_atoms, v)
            uh.append((mgr, mgr.triggers.add_handler('changes', self._mdff_changed_cb)))


    def _add_fixed_atoms_from_distance_restraints(self, mobile_atoms, fixed_atoms):
        '''
        If we have a distance restraint where only one atom is mobile, we must
        make sure that the other atom is included in the fixed selection to
        prevent it drifting unrestrained.
        '''
        dr_m = self.distance_restraint_mgr
        drs = dr_m.atoms_restraints(mobile_atoms)
        from chimerax.atomic import concatenate
        dr_atoms = concatenate(drs.atoms, remove_duplicates=True)
        remainder = dr_atoms.subtract(dr_atoms.intersect(mobile_atoms))
        remainder = remainder.subtract(remainder.intersect(fixed_atoms))
        fixed_atoms = concatenate((fixed_atoms, remainder.unique_residues.atoms))
        return fixed_atoms

    def expand_mobile_selection(self, core_atoms, expansion_mode):
        '''
        Expand the set of atoms defined by core_atoms according to a set of
        rules. After the initial rule-based expansion, the selection will be
        further expanded to encompass all residues coming within
        :attr:`IsoldeParams.soft_shell_cutoff_distance` Angstroms of the
        expanded selection. The selection returned will always consist of whole
        residues.

        Args:
            * core_atoms:
                - a :py:class:`chimerax.Atoms` instance
            * expansion_mode:
                - 'extend': each contiguous set of residues is extended
                  backwards and forwards along the chain by the number of
                  residues set by
                  :attr:`IsoldeParams.num_selection_padding_residues`
                - other modes will be added later
        '''
        from .. import selections
        iparams = self.isolde_params
        if expansion_mode == 'extend':
            sel = selections.expand_selection_along_chains(core_atoms,
                iparams.num_selection_padding_residues)
        else:
            raise TypeError('Unrecognised expansion mode!')
        shell = selections.get_shell_of_residues(sel.unique_residues,
            iparams.soft_shell_cutoff_distance).atoms
        from chimerax.atomic import concatenate
        merged_sel = concatenate((sel, shell), remove_duplicates=True)
        return merged_sel

    def prepare_sim_visualisation(self):
        '''
        Set up a simulation-friendly visualisation of the construct.
        '''
        m = self.model
        sc = self.sim_construct
        m.residues.ribbon_displays = False
        fixed_bonds = sc.fixed_atoms.intra_bonds
        fixed_bonds.radii *= self.isolde_params.fixed_bond_radius_ratio
        from ..constants import control
        sc.all_atoms.hides &= ~control.HIDE_ISOLDE


    #######################################
    # CALLBACKS
    #######################################

    def _sim_end_cb(self, *_):
        for mgr, handler in self._update_handlers:
            mgr.triggers.remove_handler(handler)
        self._update_handlers = []
        self._pr_sim_end_cb()
        self._dr_sim_end_cb()
        self._pdr_sim_end_cb()
        self._tug_sim_end_cb()
        self._rama_a_sim_end_cb()
        self._rota_a_sim_end_cb()
        self._mdff_sim_end_cb()
        rt = self._revert_to
        if rt == 'checkpoint':
            self._current_checkpoint.revert(update_sim=False)
        elif rt == 'start':
            self._starting_checkpoint.revert(update_sim=False)
        self.sim_construct.revert_visualisation()

    def _rama_a_sim_end_cb(self, *_):
        self.rama_annotator.track_whole_model = True
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def _rota_a_sim_end_cb(self, *_):
        self.rota_annotator.track_whole_model = True
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def _pr_changed_cb(self, trigger_name, changes):
        mgr, changes = changes
        change_types = list(changes.keys())
        from chimerax.core.atomic import concatenate
        changeds = []
        if 'target changed' in change_types:
            changeds.append(changes['target changed'])
        if 'enabled/disabled' in change_types:
            changeds.append(changes['enabled/disabled'])
        if 'spring constant changed' in change_types:
            changeds.append(changes['spring constant changed'])
        if len(changeds):
            all_changeds = concatenate(changeds, remove_duplicates=True)
            # limit to restraints that are actually in the simulation
            # TODO: might be better to just ignore -1 indices in the update_... functions
            all_changeds = all_changeds[all_changeds.sim_indices != -1]
            self.sim_handler.update_position_restraints(all_changeds)

    def _pr_sim_end_cb(self, *_):
        restraints = self.position_restraint_mgr.get_restraints(self.sim_construct.all_atoms)
        restraints.clear_sim_indices()
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def _dr_changed_cb(self, trigger_name, changes):
        mgr, changes = changes
        change_types = list(changes.keys())
        from chimerax.core.atomic import concatenate
        if 'created' in change_types:
            # avoid double counting
            created = changes['created']
            created = created[created.sim_indices == -1]
            # add only restraints with both atoms in the sim
            all_atoms = self.sim_construct.all_atoms
            indices = numpy.array([all_atoms.indices(atoms) for atoms in created.atoms])
            created = created[numpy.all(indices != -1, axis=0)]
            self.sim_handler.add_distance_restraints(created)
        changeds = []
        if 'target changed' in change_types:
            changeds.append(changes['target changed'])
        if 'enabled/disabled' in change_types:
            changeds.append(changes['enabled/disabled'])
        if 'spring constant changed' in change_types:
            changeds.append(changes['spring constant changed'])
        if len(changeds):
            all_changeds = concatenate(changeds, remove_duplicates=True)
            all_changeds = all_changeds[all_changeds.sim_indices != -1]
            self.sim_handler.update_distance_restraints(all_changeds)


    def _dr_sim_end_cb(self, *_):
        restraints = self.distance_restraint_mgr.intra_restraints(self.sim_construct.all_atoms)
        restraints.clear_sim_indices()
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER


    def _pdr_changed_cb(self, trigger_name, changes):
        mgr, changes = changes
        change_types = list(changes.keys())
        from chimerax.core.atomic import concatenate
        changeds = []
        if 'created' in change_types:
            # avoid double counting
            created = changes['created']
            created = created[created.sim_indices == -1]
            # add only restraints with both atoms in the sim
            all_atoms = self.sim_construct.all_atoms
            indices = numpy.array([all_atoms.indices(atoms) for atoms in created.atoms])
            created = created[numpy.all(indices != -1, axis=0)]
            self.sim_handler.add_dihedral_restraints(created)
        if 'target changed' in change_types:
            changeds.append(changes['target changed'])
        if 'enabled/disabled' in change_types:
            changeds.append(changes['enabled/disabled'])
        if 'spring constant changed' in change_types:
            changeds.append(changes['spring constant changed'])
        if len(changeds):
            all_changeds = concatenate(changeds, remove_duplicates=True)
            all_changeds = all_changeds[all_changeds.sim_indices != -1]
            self.sim_handler.update_dihedral_restraints(all_changeds)


    def _pdr_sim_end_cb(self, *_):
        restraints = self.proper_dihedral_restraint_mgr.get_all_restraints_for_residues(self.sim_construct.all_atoms.unique_residues)
        restraints.clear_sim_indices()
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER


    def _tug_changed_cb(self, trigger_name, changes):
        mgr, changes = changes
        change_types = list(changes.keys())
        from chimerax.core.atomic import concatenate
        changeds = []
        if 'target changed' in change_types:
            changeds.append(changes['target changed'])
        if 'enabled/disabled' in change_types:
            changeds.append(changes['enabled/disabled'])
        if 'spring constant changed' in change_types:
            changeds.append(changes['spring constant changed'])
        if len(changeds):
            all_changeds = concatenate(changeds, remove_duplicates=True)
            all_changeds = all_changeds[all_changeds.sim_indices != -1]
            self.sim_handler.update_tuggables(all_changeds)

    def _tug_sim_end_cb(self, *_):
        tuggables = self.tuggable_atoms_mgr.get_tuggables(self.sim_construct.all_atoms)
        tuggables.clear_sim_indices()
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def _mdff_changed_cb(self, trigger_name, changes):
        mgr, changes = changes
        change_types = list(changes.keys())
        from chimerax.atomic import concatenate
        changeds = []
        if 'enabled/disabled' in change_types:
            changeds.append(changes['enabled/disabled'])
        if 'spring constant changed' in change_types:
            changeds.append(changes['spring constant changed'])
        if len(changeds):
            all_changeds = concatenate(changeds, remove_duplicates=True)
            all_changeds = all_changeds[all_changeds.sim_indices != -1]
            self.sim_handler.update_mdff_atoms(all_changeds, mgr.volume)

    def _mdff_sim_end_cb(self, *_):
        for v, mgr in self.mdff_mgrs.items():
            mdff_atoms = mgr.get_mdff_atoms(self.sim_construct.all_atoms)
            mdff_atoms.clear_sim_indices()
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER




class Sim_Handler:
    '''
    Responsible for creating a :py:class:`openmm.Simulation`, instantiating and
    populating the custom force objects, managing the creation and calling of
    :py:class:`OpenMM_Thread_Handler`, and generally handling all the OpenMM
    side of simulation management.
    '''
    def __init__(self, session, sim_params, sim_construct):
        '''
        Prepares the simulation topology parameters and construct, and
        initialises the necessary Force objects to handle restraints. Most
        restraint forces must be populated using e.g. add_dihedral_restraints()
        before initialising the context and beginning the simulation. While it
        is possible to add new restraints *after* the simulation has started,
        in general this is only advisable in situations where it is impossible
        or impractical to define all possible restraints in advance (since each
        addition requires a costly reinitialisation of the simulation context).
        For example, the performance penalty to having all possible position
        restraints pre-defined (and mostly disabled) is minimal, but it is
        not practical to pre-define distance restraints between all possible
        atom pairs.

        Args:
            * session:
                - the ChimeraX session object
            * sim_params:
                - a :py:class:`SimParams` instance
            * sim_construct:
                - a :py:class:`Sim_Construct` instance
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

        # All custom forces in the simulation
        self.all_forces = []

        # A Volume: LinearInterpMapForce dict covering all MDFF forces
        self.mdff_forces = {}

        # Overall simulation topology
        top, residue_templates = self.create_openmm_topology(atoms, sim_construct.residue_templates)
        self._topology = top

        self._temperature = sim_params.temperature

        system = self._system = self._create_openmm_system(ff, top,
            sim_params, residue_templates)
        self.set_fixed_atoms(sim_construct.fixed_atoms)
        self._thread_handler = None

        # CustomExternalForce handling mouse and haptic interactions
        self._tugging_force = None

        self._force_update_pending = False
        self._context_reinit_pending = False
        self._minimize = False
        self._current_mode = 'min' # One of "min" or "equil"

        trigger_names = (
            'sim started',
            'coord update',
            'sim paused',
            'sim resumed',
            'sim terminated',
        )
        from chimerax.core.triggerset import TriggerSet
        t = self._triggers = TriggerSet()
        for name in trigger_names:
            t.add_trigger(name)

    @property
    def triggers(self):
        '''
        A :py:class:`chimerax.TriggerSet` allowing callbacks to be fired on
        key events. See :func:`triggers.trigger_names` for a list of available
        triggers.
        '''
        return self._triggers

    @property
    def temperature(self):
        ''' Get/set the simulation temperature in Kelvin. '''
        if not self.sim_running:
            return self._temperature
        t = self.sim_handler._simulation.integrator.getTemperature()
        return t.value_in_unit(defaults.OPENMM_TEMPERATURE_UNIT)

    @temperature.setter
    def temperature(self, temperature):
        self._simulation.integrator.setTemperature(temperature)

    @property
    def minimize(self):
        ''' Force the simulation to continue minimizing indefinitely. '''
        return self._minimize

    @minimize.setter
    def minimize(self, flag):
        self._minimize=flag

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

    def initialize_restraint_forces(self, amber_cmap=True, tugging=True, position_restraints=True,
        distance_restraints=True, dihedral_restraints=True):
        '''
        Create the selected Force objects, and add them to the System. No
        restraints are added at this stage - these should be added using
        (e.g.) add_dihedral_restraints().

        Args:
            * amber_cmap:
                - use CMAP corrections to the AMBER forcefield as per
                  http://pubs.acs.org/doi/pdf/10.1021/acs.jctc.5b00662.
                  Should only be used with the combination of AMBER and
                  implicit solvent.
            * tugging:
                - Add interactive tugging force (implemented as a
                  :py:class:`TopOutRestraintForce`)
            * position_restraints:
                - Add position restraints force (implemented as a
                  :py:class:`TopOutRestraintForce`)
            * distance_restraints:
                - Add distance restraints force (implemented as a
                  :py:class:`TopOutBondForce`)
            * dihedral_restraints:
                - Add dihedral restraints force (implemented as a
                  :py:class:FlatBottomTorsionRestraintForce)
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
            self.initialize_dihedral_restraint_force()

    def initialize_mdff_forces(self, volumes):
        '''
        Add a :py:class:`LinearInterpMapForce` for each :py:class:`Volume`
        instance.

        Args:
            * volumes:
                - an iterable of :py:class:`Volume` instances
        '''
        for v in volumes:
            self.initialize_mdff_force(v)

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
                params.fixed_integrator_timestep,
                )
        elif integrator == openmm.MTSIntegrator:
            raise RuntimeError('Multiple timestepping not yet implemented!')
        else:
            raise RuntimeError('Unrecognised or unsupported integrator: {}!'.format(integrator))
        _integrator = integrator(*integrator_params)
        _integrator.setConstraintTolerance(params.constraint_tolerance)
        return _integrator

    def start_sim(self):
        '''
        Start the main simulation loop. Automatically runs a minimisation, then
        switches to equilibration once minimisation is converged.
        '''
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
        elif self.minimize:
            f = th.minimize
            f_args=[]
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
            self.triggers.activate_trigger('sim terminated', None)
            return
        if not self._pause:
            self._repeat_step()

    def push_coords_to_sim(self, coords=None):
        '''
        Change the atomic coordinates within the simulation.

        Args:
            * coords:
                - an array of coordinates in Angstroms, or None. If None, the
                  current coordinates of the simulation construct in ChimeraX
                  will be used.
        '''
        if not self._sim_running:
            raise TypeError('No simulation running!')
        if coords is None:
            coords = self._atoms.coords
        self._pending_coords = coords
        self.triggers.add_handler('coord update', self._push_coords_to_sim)

    def _push_coords_to_sim(self, *_):
        self.thread_handler.coords = self._pending_coords
        self._pending_coords = None
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    @property
    def thread_handler(self):
        '''
        Returns the :py:class:`OpenMM_Thread_Handler` object.
        '''
        return self._thread_handler

    def toggle_pause(self):
        '''Pause the simulation if currently running, or resume if paused'''
        self.pause = not self.pause

    @property
    def pause(self):
        '''Get/set the current pause state.'''
        return self._pause

    @pause.setter
    def pause(self, flag):
        if not self._sim_running:
            raise TypeError('No simulation running!')
        if flag != self._pause:
            self._pause = flag
            if flag:
                self.triggers.activate_trigger('sim paused', None)
            else:
                self.triggers.activate_trigger('sim resumed', None)
                self._update_coordinates_and_repeat()

    def stop(self):
        '''
        Stop the simulation. This will destroy the
        :cpp:class:`OpenMM_Thread_Handler` object, rendering the Python class
        unusable.
        '''
        self._stop = True
        if self.pause:
            self._thread_handler.delete()
            self._thread_handler = None
            self._simulation = None
            self._sim_running = False
            self.triggers.activate_trigger('sim terminated', None)

    @property
    def sim_running(self):
        ''' Is the simulation curently running (i.e. started and not destroyed)? '''
        return self._sim_running

    def force_update_needed(self):
        '''
        This must be called after any changes to force objects to ensure
        the changes are pushed to the simulation context. This happens
        automatically when changes are made through the :py:class:`Sim_Manager`
        API.
        '''
        if not self.sim_running:
            return
        if self._paused:
            self._update_forces_in_context_if_needed()
        else:
            self._force_update_pending = True


    def _update_forces_in_context_if_needed(self):
        if self._context_reinit_pending:
            # defer until the reinit is done
            return
        context = self._context
        for f in self.all_forces:
            if f.update_needed:
                f.updateParametersInContext(context)
                f.update_needed = False
        self._force_update_pending = False

    def context_reinit_needed(self):
        '''
        If the number of particles, bonds etc. in any force object changes, the
        context must be reinitialised. Be aware that while reinitialisation
        happens in the thread (so GUI performance is not affected), it will
        pause the simulation for (typically) a few seconds, so should be saved
        for situations where there is no other option. If the simulation has
        not yet started this call will be ignored; otherwise the context will
        be reinitialised on the next iteration of the main loop.
        '''
        if not self.sim_running:
            return
        # if self._paused:
        #     self._reinitialize_context()
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
        '''
        Add CMAP correction terms for AMBER force field.

        Args:
            * ramas:
                - a :py:class:`Ramas` instance covering all mobile residues
        '''
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

    def initialize_dihedral_restraint_force(self):
        '''
        Just initialise the restraint force. Must be called before the
        simulation starts, and before any restraints are added.
        '''
        from .custom_forces import FlatBottomTorsionRestraintForce
        df = self._dihedral_restraint_force = FlatBottomTorsionRestraintForce()
        self._system.addForce(df)
        self.all_forces.append(df)

    def add_dihedral_restraints(self, restraints):
        '''
        Add a set of dihedral restraints to the simulation. This sets the
        :attr:`sim_index` for each restraint so it knows its own place in the
        simulation for later updating purposes. Automatically calls
        :func:`context_reinit_needed()`.

        Args:
            * restraints:
                - Any python restraints :py:class:`chimerax.Collection` instance
                  based on :cpp:class:`Dihedral_Restraint`. At present this is
                  limited to :py:class:`Proper_Dihedral_Restraints`, but
                  expect others (e.g. chirals) to be added with time.
        '''
        force = self._dihedral_restraint_force
        all_atoms = self._atoms
        dihedral_atoms = restraints.dihedrals.atoms
        atom_indices = [all_atoms.indices(atoms) for atoms in dihedral_atoms]
        restraints.sim_indices = force.add_torsions(atom_indices,
            restraints.enableds, restraints.spring_constants, restraints.targets, restraints.cutoffs)
        self.context_reinit_needed()

    def add_dihedral_restraint(self, restraint):
        '''
        Singular form of :func:`add_dihedral_restraints`.

        Args:
            * restraint:
                - A :py:class:`Proper_Dihedral_Restraint` instance, or any other
                  Python restraint class based on :cpp:class:`Dihedral_Restraint`
        '''
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
        '''
        Update the simulation to match the parameters (target angles,
        cutoffs, spring constants and enabled states) of the given restraints.

        Args:
            * restraints:
                - A :py:class:`Proper_Dihedral_Restraints` instance. Any
                  restraints that are not part of the current simulation will
                  be ignored.
        '''
        force = self._dihedral_restraint_force
        force.update_targets(restraints.sim_indices,
            restraints.enableds, restraints.spring_constants, restraints.targets, restraints.cutoffs)
        self.force_update_needed()

    def update_dihedral_restraint(self, restraint):
        '''
        Update the simulation to match the parameters (target angles,
        cutoffs, spring constants and enabled states) of the given restraint.

        Args:
            * restraint:
                - A :py:class:`Proper_Dihedral_Restraint` instance. If the
                  restraint is not part of the current simulation it will be
                  ignored.
        '''
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
        '''
        Just initialise the force, and set its limiting magnitude. Must be called
        before the simulation starts, and before any restraints are added.

        Args:
            * max_force:
                - the maximum allowable force, in :math:`kJ mol^{-1} nm^{-1}`
        '''
        from .custom_forces import TopOutBondForce
        tf = self._distance_restraints_force = TopOutBondForce(max_force)
        self._system.addForce(tf)
        self.all_forces.append(tf)
        return tf

    def add_distance_restraints(self, restraints):
        '''
        Add a set of distance restraints to the simulation. Sets
        :attr:`sim_index` for each restraint so it knows its place in the
        simulation. Automatically calls :func:`context_reinit_needed`

        Args:
            * restraints:
                - a :py:class:`Distance_Restraints` instance
        '''
        force = self._distance_restraints_force
        all_atoms = self._atoms
        dr_atoms = restraints.atoms
        indices = [all_atoms.indices(atoms) for atoms in dr_atoms]
        restraints.sim_indices = force.add_bonds(indices,
            restraints.enableds, restraints.spring_constants, restraints.targets/10)
        self.context_reinit_needed()

    def add_distance_restraint(self, restraint):
        '''
        Add a single distance restraint to the simulation. Sets
        :attr:`sim_index` for the restraint so it knows its place in the
        simulation. Automatically calls :func:`context_reinit_needed`

        Args:
            * restraint:
                - a :py:class:`Distance_Restraint` instance
        '''
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
        '''
        Update the simulation to reflect the current parameters (target
        distances, spring constants, enabled/disabled states) of the given
        restraints.

        Args:
            * restraints:
                - a :py:class:`Distance_Restraints` instance
        '''
        force = self._distance_restraints_force
        force.update_targets(restraints.sim_indices,
            restraints.enableds, restraints.spring_constants, restraints.targets/10)
        self.force_update_needed()

    def update_distance_restraint(self, restraint):
        '''
        Update the simulation to reflect the current parameters (target
        distance, spring constant, enabled/disabled state) of the given
        restraint.

        Args:
            * restraint:
                - a :py:class:`Distance_Restraint` instance
        '''
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
        '''
        Just initialise the force, and set its limiting magnitude. Must be called
        before the simulation starts, and before any restraints are added.

        Args:
            * max_force:
                - the maximum allowable force, in :math:`kJ mol^{-1} nm^{-1}`
        '''
        from .custom_forces import TopOutRestraintForce
        rf = self._position_restraints_force = TopOutRestraintForce(max_force)
        self._system.addForce(rf)
        self.all_forces.append(rf)
        return rf

    def add_position_restraints(self, restraints):
        '''
        Add a set of position restraints to the simulation. Sets
        :attr:`sim_index` for each restraint so it knows its place in the
        simulation. Automatically calls :func:`context_reinit_needed`

        Args:
            * restraints:
                - a :py:class:`Position_Restraints` instance
        '''
        force = self._position_restraints_force
        all_atoms = self._atoms
        indices = all_atoms.indices(restraints.atoms)
        restraints.sim_indices = force.add_particles(indices,
            restraints.enableds, restraints.spring_constants, restraints.targets/10)
        self.context_reinit_needed()

    def add_position_restraint(self, restraint):
        '''
        Add a single position restraint to the simulation. Sets
        :attr:`sim_index` for the restraint so it knows its place in the
        simulation. Automatically calls :func:`context_reinit_needed`

        Args:
            * restraint:
                - a :py:class:`Position_Restraint` instance
        '''
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
        '''
        Update the simulation to reflect the current parameters (target
        positions, spring constants, enabled/disabled states) of the given
        restraints.

        Args:
            * restraints:
                - a :py:class:`Position_Restraints` instance
        '''
        force = self._position_restraints_force
        force.update_targets(restraints.sim_indices,
            restraints.enableds, restraints.spring_constants, restraints.targets/10)
        self.force_update_needed()

    def update_position_restraint(self, restraint):
        '''
        Update the simulation to reflect the current parameters (target
        position, spring constant, enabled/disabled state) of the given
        restraint.

        Args:
            * restraint:
                - a :py:class:`Position_Restraint` instance
        '''
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
        '''
        Just initialise the force, and set its limiting magnitude. Must be called
        before the simulation starts, and before any restraints are added.

        Args:
            * max_force:
                - the maximum allowable force, in :math:`kJ mol^{-1} nm^{-1}`
        '''
        from .custom_forces import TopOutRestraintForce
        f = self._tugging_force = TopOutRestraintForce(max_force)
        self._system.addForce(f)
        self.all_forces.append(f)

    def add_tuggables(self, tuggables):
        '''
        Add a set of tuggable atom proxies to the simulation. Sets
        :attr:`sim_index` for each tuggable so it knows its place in the
        simulation. Automatically calls :func:`context_reinit_needed`

        Args:
            * tuggables:
                - a :py:class:`Tuggable_Atoms` instance
        '''
        force = self._tugging_force
        all_atoms = self._atoms
        indices = all_atoms.indices(tuggables.atoms)
        tuggables.sim_indices = force.add_particles(indices,
            tuggables.enableds, tuggables.spring_constants, tuggables.targets/10)
        self.context_reinit_needed()

    def add_tuggable(self, tuggable):
        '''
        Add a single tuggable atom proxy to the simulation. Sets
        :attr:`sim_index` so the tuggable knows its place in the simulation.
        Automatically calls :func:`context_reinit_needed`

        Args:
            * tuggable:
                - a :py:class:`Tuggable_Atom` instance
        '''
        force = self._tugging_force
        index = self._all_atoms.index(tuggable.atom)
        target = (tuggable.target/10).tolist()
        tuggable.sim_index = force.addParticle(index,
            (float(tuggable.enabled), tuggable.spring_constant, *target))
        self.context_reinit_needed()

        ##
        # During simulation
        ##

    def update_tuggables(self, tuggables):
        '''
        Update the simulation to reflect the current parameters (target
        positions, spring constants, enabled/disabled states) of the given
        tuggables.

        Args:
            * tuggables:
                - a :py:class:`Tuggable_Atoms` instance
        '''
        force = self._tugging_force
        force.update_targets(tuggables.sim_indices,
            tuggables.enableds, tuggables.spring_constants, tuggables.targets/10)
        self.force_update_needed()

    def update_tuggable(self, tuggable):
        '''
        Update the simulation to reflect the current parameters (target
        positions, spring constants, enabled/disabled states) of the given
        tuggable.

        Args:
            * tuggable:
                - a :py:class:`Tuggable_Atom` instance
        '''
        force = self._tugging_force
        force.update_target(tuggable.sim_index,
            tuggable.enabled, tuggable.spring_constant, tuggable.target/10)
        self.force_update_needed()

    ####
    # MDFF forces
    ####

        ##
        # Before simulation starts
        ##

    def initialize_mdff_force(self, volume):
        '''
        Prepare an MDFF map from a :py:class:`chimerax.Volume` instance.  The
        Volume instance is expected to have a single :attr:`region` applied that
        is big enough to cover the expected range of motion of the MDFF atoms
        that will be coupled to it (and for performance/memory reasons, ideally
        not too much bigger). Must be called before the simulation is started,
        and before any MDFF atoms are added.

        Args:
            * volume:
                - a :py:class:`chimerax.Volume` instance
        '''
        from .custom_forces import LinearInterpMapForce
        v = volume
        region = v.region
        # Ensure that the region ijk step size is [1,1,1]
        v.new_region(ijk_min=region[0], ijk_max=region[1], ijk_step=[1,1,1])
        data = v.region_matrix()
        from chimerax.core.geometry import Place
        tf = v.data.xyz_to_ijk_transform
        # Shift the transform to the origin of the region
        region_tf = Place(axes=tf.axes(), origin = tf.origin() -
            v.data.xyz_to_ijk(v.region_origin_and_step(v.region)[0]))
        f = LinearInterpMapForce(data, region_tf.matrix, units='angstroms')
        self.all_forces.append(f)
        self._system.addForce(f)
        self.mdff_forces[v] = f

    def add_mdff_atoms(self, mdff_atoms, volume):
        '''
        Add a set of MDFF atom proxies to the force associated with the given
        volume. Automatically calls :func:`context_reinit_needed`

        Args:
            * mdff_atoms:
                - a :py:class:`MDFF_Atoms` instance
            * volume:
                - the :py:class:`chimerax.Volume` instance that was used to
                  create the target force.
        '''
        f = self.mdff_forces[volume]
        all_atoms = self._atoms
        indices = all_atoms.indices(mdff_atoms.atoms)
        mdff_atoms.sim_indices = f.add_atoms(indices,
            mdff_atoms.coupling_constants, mdff_atoms.enableds)
        self.context_reinit_needed()

    def add_mdff_atom(self, mdff_atom, volume):
        '''
        Add a singl MDFF atom proxy to the force associated with the given
        volume. Automatically calls :func:`context_reinit_needed`

        Args:
            * mdff_atom:
                - a :py:class:`MDFF_Atom` instance
            * volume:
                - the :py:class:`chimerax.Volume` instance that was used to
                  create the target force.
        '''
        f = self.mdff_forces[volume]
        all_atoms = self._atoms
        index = all_atoms.index(mdff_atom.atom)
        mdff_atom.sim_index = f.addParticle(index,
            (mdff_atom.coupling_constant, float(mdff_atom.enabled)))
        self.context_reinit_needed()

        ##
        # During simulation
        ##

    def set_mdff_global_k(self, volume, k):
        '''
        Set the global coupling constant for the MDFF force associated with
        the given volume. Automatically calls :func:`force_update_needed`.

        Args:
            * volume:
                - the :py:class:`chimerax.Volume` instance that was used to
                  create the target force
            * k:
                - the new coupling constant, in
                  :math:`kJ mol^{-1} (\\text{map density unit})^{-1} nm^3`
        '''
        f = self.mdff_forces[volume]
        f.set_global_k(k)
        self.force_update_needed()

    def update_mdff_atoms(self, mdff_atoms, volume):
        '''
        Update the simulation to reflect the new parameters (individual
        coupling constants, enabled/disabled states) for the given MDFF atom
        proxies.

        Args:
            * mdff_atoms:
                - a :py:class:`MDFF_Atoms` instance
            * volume:
                - the :py:class:`chimerax.Volume` instance that was used to
                  create the target force.
        '''
        f = self.mdff_forces[volume]
        f.update_atoms(mdff_atoms.sim_indices,
            mdff_atoms.coupling_constants, mdff_atoms.enableds)
        self.force_update_needed()

    def update_mdff_atom(self, mdff_atom, volume):
        '''
        Update the simulation to reflect the new parameters (individual
        coupling constants, enabled/disabled states) for the given MDFF atom
        proxy.

        Args:
            * mdff_atom:
                - a :py:class:`MDFF_Atom` instance
            * volume:
                - the :py:class:`chimerax.Volume` instance that was used to
                  create the target force.
        '''
        f = self.mdff_forces[volume]
        f.update_atom(mdff_atom.sim_index,
            mdff_atom.coupling_constant, mdff_atom.enabled)


    def set_fixed_atoms(self, fixed_atoms):
        ''' Fix the desired atoms rigidly in space. NOTE: a fixed atom can not
        be bonded to a mobile atom via a rigid bond. In most cases this means
        that you cannot fix a hydrogen without fixing the heavy atom that it's
        bonded to, and any fixed heavy atom must have all its hydrogens fixed.
        While atoms may in theory be fixed and un-fixed during a simulation,
        ISOLDE wasn't built with this in mind and it requires a costly
        re-initialisation of the simulation context. In most cases it's best to
        simply use strong position restraints wherever you want to
        interactively "fix" atoms.

        Args:
            * fixed_atoms:
                - a :py:class:`chimerax.Atoms` instance
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
        While atoms may in theory be fixed and un-fixed during a simulation,
        ISOLDE wasn't built with this in mind and it requires a costly
        re-initialisation of the simulation context. In most cases it's best to
        simply use strong position restraints wherever you want to
        interactively "fix" atoms.

        Args:
            * atoms:
                - a :py:class:`chimerax.Atoms` instance
        '''
        indices = self._atoms.indices(fixed_atoms).tolist()
        masses = atoms.elements.masses
        sys = self.system
        for index, mass in zip(indices, masses):
            sys.setParticleMass(index, mass)
        self.context_reinit_needed()

    def define_forcefield(self, forcefield_file_list):
        '''
        Create a :py:class:`openmm.ForceField` object from a set of OpenMM
        ffXML files.

        Args:
            * forcefield_file_list:
                - An iterable of file names.
        '''
        from simtk.openmm.app import ForceField
        ff = ForceField(*[f for f in forcefield_file_list if f is not None])
        return ff

    def create_openmm_topology(self, atoms, residue_templates):
        '''
        Generate a simulation topology from a set of atoms.

        Args:
            * atoms:
                - a :py:class:`chimerax.Atoms` instance. Residues must be
                  complete, and the atoms in each residue must be contiguous
                  in the array.
            * residue_templates:
                - A {residue_index: residue_type} dict for residues whose
                  topology is otherwise ambiguous, where residue_type is the
                  name of the residue in the forcefield definition. OpenMM
                  requires a {:py:class:`openmm.Residue`: residue_type} dict, so
                  we need to do the conversion here.
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
        '''
        Add a Generalised Born Implicit Solvent (GBIS) formulation.

        Args:
            * params:
                - a :py:class:`SimParams` instance
        '''
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

def get_available_platforms():
        from simtk.openmm import Platform
        platform_names = []
        for i in range(Platform.getNumPlatforms()):
            p = Platform.getPlatform(i)
            name = p.getName()
            platform_names.append(name)
        return platform_names
