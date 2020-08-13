# @Author: Tristan Croll <tic20>
# @Date:   26-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 11-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



import os, sys, glob
import numpy
import ctypes
from chimerax.core.state import State
from chimerax.atomic import molc
# from chimerax.atomic.molc import CFunctions, string, cptr, pyobject, \
#     set_c_pointer, pointer, size_t

CFunctions = molc.CFunctions
string = molc.string
cptr = molc.cptr
pyobject = molc.pyobject
set_c_pointer = molc.set_c_pointer
pointer = molc.pointer
size_t = molc.size_t
# object map lookups

from numpy import int8, uint8, int32, uint32, float64, float32, byte, bool as npy_bool

from simtk import unit, openmm

from ..delayed_reaction import delayed_reaction
from ..util import compiled_lib_extension
libdir = os.path.dirname(os.path.abspath(__file__))
libfile = os.path.join(libdir, '..', 'libopenmm.'+compiled_lib_extension())

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
    def __init__(self, context, params, c_pointer=None):
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
        self._smoothing = False
        self._last_smooth = False
        self._last_mode = None
        self.params = params

    @property
    def smoothing(self):
        '''
        If true, the displayed coordinates will be a smoothed average of the
        last set of equilibration steps. Note that for large values of
        sim_steps_per_gui_update this can lead to distorted geometry.
        '''
        return self._smoothing

    @smoothing.setter
    def smoothing(self, flag):
        self._smoothing = flag


    def _get_smoothing_alpha(self):
        '''
        ISOLDE uses an exponential smoothing scheme, where
        :attr:`smoothing_alpha` defines the contribution of each new set of
        coordinates to the moving average. Values are limited to the range
        (0.01..1), where 1 indicates no smoothing and 0.01 provides extremely
        strong smoothing. Values outside of this range will be automatically
        clamped. Internally, coordinates are added to the moving  average
        every 10 steps or the number of steps to the next graphics  update,
        whichever is smaller.

        Note that smoothing only affects the *visualisation* of the simulation,
        not the simulation itself. Applying energy minimisation or pausing  the
        simulation fetches the latest instantaneous coordinates and restarts the
        smoothing.
        '''
        f = c_function('openmm_thread_handler_smoothing_alpha',
            args=(ctypes.c_void_p,),
            ret=ctypes.c_double)
        return f(self._c_pointer)

    def _set_smoothing_alpha(self, alpha):
        f = c_function('set_openmm_thread_handler_smoothing_alpha',
            args=(ctypes.c_void_p, ctypes.c_double))
        f(self._c_pointer, alpha)


    smoothing_alpha = property(_get_smoothing_alpha, _set_smoothing_alpha)


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
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_bool))
        f(self._c_pointer, steps, self._smoothing)
        self._last_mode = 'equil'
        self._last_smooth = self._smoothing

    def minimize(self, tolerance=None, max_iterations=None):
        '''
        Run an energy minimization on the coordinates. If the minimisation
        converges to within tolerance, unstable will be set to False.
        Don't forget to run :func:`reinitialize_velocities` before continuing
        equilibration!

        Args:
            * tolerance:
                - Convergence tolerance, in kJ/mol/atom.
            * max_iterations:
                - Maximum number of iterations to run for before returning
                  coordinates. NOTE: minimisation runs best if this number is
                  kept large (at least a few hundred).
        '''
        if tolerance is None:
            tolerance = self.params.minimization_convergence_tol_start
        if max_iterations is None:
            max_iterations=self.params.minimization_max_iterations
        f = c_function('openmm_thread_handler_minimize',
            args = (ctypes.c_void_p, ctypes.c_double, ctypes.c_int))
        f(self._c_pointer, tolerance, max_iterations)
        self._last_mode = 'min'

    @property
    def clashing(self):
        '''
        True if any atom is experiencing an extreme force after energy
        minimisation.
        '''
        f = c_function('openmm_thread_handler_clashing',
            args=(ctypes.c_void_p,),
            ret=ctypes.c_bool)
        return f(self._c_pointer)

    @property
    def minimization_converged(self):
        f = c_function('openmm_thread_handler_converged',
            args=(ctypes.c_void_p,),
            ret=ctypes.c_bool)
        return f(self._c_pointer)

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
        if not self._smoothing or not self._last_smooth or self._last_mode !='equil':
            f = c_function('openmm_thread_handler_current_coords',
                args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p))
        else:
            f = c_function('openmm_thread_handler_smoothed_coords',
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
    def __init__(self, model, mobile_atoms, fixed_atoms, excluded_atoms=None):
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
            * excluded_atoms:
                - A :py:class:`chimerax.Atoms` instance defining any atoms to
                  be excluded from the simulation. This may be set to None.

        NOTE: the simulation will fail if mobile_atoms and fixed_atoms do not
        combine to form a set containing only complete residues. Also note that
        OpenMM does not allow a fixed atom to be rigidly bonded to a mobile
        one (in practice this means any fixed heavy atom must have its hydrogens
        also fixed).
        '''
        self.model = model

        # Chains in OpenMM must be in a single unbroken block
        residues = model.residues
        residues = residues[numpy.lexsort((residues.numbers, residues.chain_ids))]

        # Sort all the atoms according to their order in the model#
        model_atoms = residues.atoms
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
        self._excluded_atoms = excluded_atoms

        self.store_original_visualisation()
        self.surroundings = model_atoms.subtract(all_atoms)
        if excluded_atoms is not None:
            self.surroundings = self.surroundings.subtract(excluded_atoms)
            excluded_atoms[excluded_atoms.element_names=='C'].colors = [50,50,50,255]

    @property
    def mobile_atoms(self):
        '''
        A :py:class:`chimerax.Atoms` instance containing the mobile selection,
        sorted in the same order as :attr:`model.residues.atoms`.
        '''
        return self._mobile_atoms

    @property
    def mobile_heavy_atoms(self):
        '''
        A :py:class:`chimerax.Atoms` instance containing the non-hydrogen mobile
        atoms, sorted in the same order as :attr:`model.residues.atoms`
        '''
        if not hasattr(self, '_mobile_heavy_atoms'):
            ma = self.mobile_atoms
            self._mobile_heavy_atoms = ma[ma.element_names != 'H']
        return self._mobile_heavy_atoms

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

    @property
    def all_residues(self):
        '''
        A :py:class:`chimerax.Residues` instance containing all residues in the
        simulation, sorted in the same order as :attr:`model.residues`
        '''
        return self._all_atoms.unique_residues

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
        # If atoms are added or deleted while a simulation is running, the
        # stored settings will become invalid. In that case, just revert to a
        # default visualisation.
        try:
            atoms.colors, atoms.draw_modes, atoms.displays, atoms.radii = \
                self._original_atom_states
            bonds.displays, bonds.radii = \
                self._original_bond_states
            residues.ribbon_displays, residues.ribbon_colors = \
                self._original_residue_states
        except:
            from ..visualisation import default_atom_visualisation
            default_atom_visualisation(self.model)

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
        isolde_params, sim_params, excluded_residues=None, expansion_mode = 'extend'):
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
            * excluded_residues:
                - optional :py:class:`chimerax.Residues` defining residues that
                  should be excluded from the simulation (typically because
                  there is no MD parameterisation for them). These will remain
                  visible in the model and be considered for map calculations,
                  but will not have any impact whatsoever on simulations. Any
                  atom(s) directly bonded to an excluded residue will be fixed
                  in space along with their attendant hydrogen atoms.
            * expansion_mode:
                - string defining how the initial selection will be expanded.
                  For allowable options, see :func:`expand_mobile_selection`
        '''
        self.isolde = isolde
        self.model = model
        session = self.session = model.session
        logger = session.logger
        self.isolde_params = isolde_params
        self.sim_params = sim_params
        # If changes are made to the model while the simulation is paused, we
        # need to push them to the simulation before resuming
        self._pause_atom_changes_handler = None
        self._revert_to = None
        logger.status('Determining simulation layout')
        mobile_atoms = self.expand_mobile_selection(selected_atoms, expansion_mode)
        from ..selections import get_shell_of_residues
        fixed_atoms = get_shell_of_residues(mobile_atoms.unique_residues,
            isolde_params.hard_shell_cutoff_distance).atoms
        self._prepare_validation_managers(mobile_atoms)
        self._prepare_restraint_managers()
        fixed_atoms = self._add_fixed_atoms_from_distance_restraints(mobile_atoms, fixed_atoms)
        if excluded_residues is not None:
            mobile_atoms, fixed_atoms, excluded_atoms = self._add_fixed_atoms_from_excluded_residues(mobile_atoms, fixed_atoms, excluded_residues)

        sc = self.sim_construct = Sim_Construct(model, mobile_atoms, fixed_atoms, excluded_atoms)
        self.prepare_sim_visualisation()

        logger.status('Preparing simulation handler')
        self._prepare_mdff_managers()
        sh = self.sim_handler = None
        uh = self._update_handlers = []
        try:
            sh = self.sim_handler = Sim_Handler(session, sim_params, sc,
                isolde.forcefield_mgr)
        except Exception as e:
            self._sim_end_cb(None, None)
            # if isinstance(e, RuntimeError):
            #     if 'Unpararameterised' in e.args[0]:
            #         return
            if isinstance(e, ValueError):
                # If it's an error in template handling, parse out the offending
                # residue and tell ISOLDE about it
                self._parse_auto_template_error(e)
                # Return early to avoid further errors. This object will be cleaned
                # up automatically
                raise e
            else:
                # Explicit template provision (e.g. for CYS residues) throws a
                # generic Exception if the template doesn't match the residue
                # topology. Catch and handle that case, then throw it upwards as a
                # ValueError
                if self._parse_explicit_template_error(e):
                    raise ValueError(str(e))
                else:
                    raise e
        sh.triggers.add_handler('sim paused', self._sim_pause_cb)
        sh.triggers.add_handler('sim resumed', self._sim_resume_cb)

        self._initialize_restraints(uh)
        self._initialize_mdff(uh)

    @property
    def sim_running(self):
        '''
        Returns True if the simulation is running, otherwise False. Read only.
        '''
        return self.sim_handler.sim_running

    def start_reporting_performance(self, report_interval=20):
        if hasattr(self, '_performance_tracker') and self._performance_tracker is not None:
            pt = self._performance_tracker
            pt.report_interval = report_interval
            pt.start()
        else:
            pt = self._performance_tracker = Sim_Performance_Tracker(self.sim_handler, report_interval, self.session.logger.status)
            pt.start()

    def stop_reporting_performance(self):
        if hasattr(self, '_performance_tracker') and self._performance_tracker is not None:
            self._performance_tracker.stop()
        self._performance_tracker = None

    def _parse_auto_template_error(self, e):
        err_text = str(e)
        if not err_text.startswith('No template found'):
            # Not a template error. Just raise it
            return
        tokens = err_text.split()
        res_num = int(tokens[5])
        #print("Bad residue number: {}".format(tokens[5]))
        # OpenMM residue numbering starts from 1
        residue = self.sim_construct.all_residues[res_num-1]
        self.isolde._handle_bad_template(residue)

    def _parse_explicit_template_error(self, e):
        err_text = str(e)
        if not err_text.startswith('User-supplied template'):
            return False
        tokens = err_text.split()
        res_num = int(tokens[8])
        #print("Bad residue number: {}".format(tokens[8]))
        # OpenMM residue numbering starts from 1
        residue = self.sim_construct.all_residues[res_num-1]
        self.isolde._handle_bad_template(residue)
        return True


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
        self._revert_to = revert
        self.sim_handler.stop()

    @property
    def pause(self):
        return self.sim_handler.pause

    @pause.setter
    def pause(self, flag):
        self.sim_handler.pause = flag

    def _sim_pause_cb(self, *_):
        '''
        Make sure we capture any changes to the model while paused.
        '''
        if self._pause_atom_changes_handler is None:
            self._pause_atom_changes_handler = self.model.triggers.add_handler(
                'changes', self._atom_changes_while_paused_cb
            )

    def _sim_resume_cb(self, *_):
        if self._pause_atom_changes_handler is not None:
            self.model.triggers.remove_handler(self._pause_atom_changes_handler)
            self._pause_atom_changes_handler = None


    def toggle_pause(self):
        '''
        Pause/resume the simulation.
        '''
        self.pause = not self.pause

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
        rama_a = self.RamaAnnotator = sx.get_RamaAnnotator(m)
        rota_a = self.rota_annotator = sx.get_rota_annotator(m)
        rama_a.restrict_to_selected_residues(mobile_res)
        rota_a.restrict_to_selected_residues(mobile_res)

    def _prepare_restraint_managers(self):
        from .. import session_extensions as sx
        m = self.model
        self.chiral_restraint_mgr = sx.get_chiral_restraint_mgr(m)
        self.proper_dihedral_restraint_mgr = sx.get_proper_dihedral_restraint_mgr(m)
        self.adaptive_dihedral_restraint_mgr = sx.get_adaptive_dihedral_restraint_mgr(m)
        self.position_restraint_mgr = sx.get_position_restraint_mgr(m)
        self.tuggable_atoms_mgr = sx.get_tuggable_atoms_mgr(m, allow_hydrogens=self.sim_params.tug_hydrogens)
        self.distance_restraint_mgr = sx.get_distance_restraint_mgr(m)
        adrms = self.adaptive_distance_restraint_mgrs = sx.get_all_adaptive_distance_restraint_mgrs(m)
        if not len(adrms):
            # Always maintain at least the default manager
            self.adaptive_distance_restraint_mgrs = [sx.get_adaptive_distance_restraint_mgr(m)]

    def _initialize_restraints(self, update_handlers):
        sh = self.sim_handler
        sc = self.sim_construct
        logger = self.session.logger
        logger.status('Initialising restraints...')
        sim_params = self.sim_params
        uh = update_handlers
        mobile_res = sc.mobile_atoms.unique_residues
        amber_cmap = False
        if self.sim_params.forcefield == 'amber14':
            amber_cmap = True
        sh.initialize_restraint_forces(amber_cmap)
        if (amber_cmap):
            from .. import session_extensions as sx
            rama_mgr = sx.get_ramachandran_mgr(self.session)
            ramas = rama_mgr.get_ramas(mobile_res)
            ramas = ramas[ramas.valids]
            # ramas = ramas.restrict_to_sel(self.sim_construct.all_atoms)
            sh.add_amber_cmap_torsions(ramas)

        cr_m = self.chiral_restraint_mgr
        crs = cr_m.add_restraints_by_atoms(sc.mobile_heavy_atoms)
        # crs = crs.restrict_to_sel(sc.all_atoms)
        sh.add_dihedral_restraints(crs)
        uh.append((cr_m, cr_m.triggers.add_handler('changes', self._dihe_r_changed_cb)))

        pdr_m = self.proper_dihedral_restraint_mgr
        pdrs = pdr_m.add_all_defined_restraints_for_residues(mobile_res)
        if sim_params.restrain_peptide_omegas:
            import numpy
            omega_rs = pdr_m.get_restraints_by_residues_and_name(mobile_res, 'omega')
            # Only update targets for restraints which aren't already enabled
            omega_rs = omega_rs[numpy.logical_not(omega_rs.enableds)]
            if len(omega_rs):
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
        uh.append((pdr_m, pdr_m.triggers.add_handler('changes', self._dihe_r_changed_cb)))

        apdr_m = self.adaptive_dihedral_restraint_mgr
        apdrs = apdr_m.get_all_restraints_for_residues(mobile_res)
        sh.add_adaptive_dihedral_restraints(apdrs)
        uh.append((apdr_m, apdr_m.triggers.add_handler('changes', self._adaptive_dihe_r_changed_cb)))

        dr_m = self.distance_restraint_mgr
        # Pre-create all restraints necessary for secondary structure manipulation
        dr_m.add_ss_restraints(sc.all_atoms.unique_residues)
        drs = dr_m.atoms_restraints(sc.mobile_atoms)
        sh.add_distance_restraints(drs)
        uh.append((dr_m, dr_m.triggers.add_handler('changes', self._dr_changed_cb)))

        for adr_m in self.adaptive_distance_restraint_mgrs:
            adrs = adr_m.atoms_restraints(sc.mobile_atoms)
            sh.add_adaptive_distance_restraints(adrs)
            uh.append((adr_m, adr_m.triggers.add_handler('changes', self._adr_changed_cb)))

        pr_m = self.position_restraint_mgr
        prs = pr_m.add_restraints(sc.mobile_atoms)
        sh.add_position_restraints(prs)
        uh.append((pr_m, pr_m.triggers.add_handler('changes', self._pr_changed_cb)))

        ta_m = self.tuggable_atoms_mgr
        tuggables = ta_m.add_tuggables(sc.mobile_atoms)
        uh.append((ta_m, ta_m.triggers.add_handler('changes', self._tug_changed_cb)))
        sh.add_tuggables(tuggables)

    def _prepare_mdff_managers(self):
        '''
        Take experimental density maps and convert to MDFF potential energy
        fields. For crystallographic data, only the map named "MDFF potential"
        will be used, since it is guaranteed to exclude the free reflections.
        '''
        from .. import session_extensions as sx
        isolde_params = self.isolde.params
        m = self.model
        mdff_mgr_map = self.mdff_mgrs = {}
        from chimerax.map import Volume
        from chimerax.clipper.symmetry import get_symmetry_handler
        sh = get_symmetry_handler(m)
        if sh is None:
            return
        for v in sh.map_mgr.all_maps:
            mgr = sx.get_mdff_mgr(m, v, create=False)
            if mgr is not None and mgr.enabled:
                mdff_mgr_map[v] = mgr
        sh.isolate_and_cover_selection(self.sim_construct.mobile_atoms,
            include_surrounding_residues = 0,
            show_context = isolde_params.hard_shell_cutoff_distance,
            mask_radius = isolde_params.standard_map_mask_cutoff,
            extra_padding = 5,
            hide_surrounds = isolde_params.hide_surroundings_during_sim,
            focus = False)


    def _initialize_mdff(self, update_handlers):
        sh = self.sim_handler
        sc = self.sim_construct
        sp = self.sim_params
        uh = update_handlers
        mdff_mgrs = self.mdff_mgrs
        sh.initialize_mdff_forces(list(mdff_mgrs.keys()))
        for v in sh.mdff_forces.keys():
            mgr = mdff_mgrs[v]
            mdff_atoms = mgr.add_mdff_atoms(sc.mobile_atoms,
                hydrogens = sp.hydrogens_feel_maps)
            sh.set_mdff_global_k(v, mgr.global_k)
            sh.add_mdff_atoms(mdff_atoms, v)
            uh.append((mgr, mgr.triggers.add_handler('changes', self._mdff_changed_cb)))
            uh.append((mgr, mgr.triggers.add_handler('global k changed', self._mdff_global_k_change_cb)))


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

    def _add_fixed_atoms_from_excluded_residues(self, mobile_atoms, fixed_atoms, excluded_residues):
        '''
        Filter out any residues that are to be excluded from the simulation, and
        fix any atoms (with their attendant hydrogens) that are directly bonded
        to an excluded residue.
        '''
        from ..molobject import residue_bonded_neighbors
        all_sim_atoms = mobile_atoms.merge(fixed_atoms)
        all_sim_res = mobile_atoms.unique_residues
        sim_excludes = excluded_residues.atoms.intersect(all_sim_atoms)
        sim_exclude_res = sim_excludes.unique_residues
        from chimerax.atomic import Atoms
        extra_fixed = []
        for r in sim_exclude_res:
            neighbors = residue_bonded_neighbors(r)
            for n in neighbors:
                if n in all_sim_res:
                    bonds = r.bonds_between(n)
                    for b in bonds:
                        for a in b.atoms:
                            if a.residue == n:
                                extra_fixed.append(a)
                                # Fix directly attached atoms, atoms directly
                                # attached to those, and associated hydrogens
                                for na in a.neighbors:
                                    if na.residue !=r:
                                        extra_fixed.append(na)
                                        if na.element.name != 'H':
                                            for nb in na.neighbors:
                                                if nb.element.name == 'H':
                                                    extra_fixed.append(nb)
        extra_fixed = Atoms(extra_fixed)
        fixed_atoms = fixed_atoms.merge(extra_fixed)
        fixed_atoms = fixed_atoms.subtract(sim_excludes)
        mobile_atoms = mobile_atoms.subtract(extra_fixed)
        mobile_atoms = mobile_atoms.subtract(sim_excludes)
        return mobile_atoms, fixed_atoms, sim_excludes





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

    def _atom_changes_while_paused_cb(self, trigger_name, changes):
        '''
        If changes are made to the model while the simulation is paused, we need
        to deal with them before resuming. Changes in coordinates need to be
        pushed to the simulation, while addition/deletion of atoms within the
        simulation construct should stop the simulation entirely.
        '''
        changes = changes[1]
        if changes.num_deleted_atoms() or len(changes.created_atoms()):
            self.sim_handler.stop(reason = 'coord length mismatch')
            return
        sh = self.sim_handler
        if 'coord changed' in changes.atom_reasons():
            self.sim_handler.push_coords_to_sim()


    def _sim_end_cb(self, trigger_name, reason):
        if self._pause_atom_changes_handler is not None:
            self.model.triggers.remove_handler(self._pause_atom_changes_handler)
        for mgr, handler in self._update_handlers:
            mgr.triggers.remove_handler(handler)
        self._update_handlers = []
        self._pr_sim_end_cb()
        self._dr_sim_end_cb()
        self._adr_sim_end_cb()
        self._dihe_r_sim_end_cb()
        self._adaptive_dihe_r_sim_end_cb()
        self._tug_sim_end_cb()
        self._rama_a_sim_end_cb()
        self._rota_a_sim_end_cb()
        self._mdff_sim_end_cb()
        if reason == 'coord length mismatch':
            rt = None
        else:
            rt = self._revert_to
        if rt == 'checkpoint':
            self._current_checkpoint.revert()
        elif rt == 'start':
            print('reverting to start')
            self._starting_checkpoint.revert()
        self.sim_construct.revert_visualisation()
        if reason == 'coord length mismatch':
            msg = ('Mismatch between number of simulated atoms and the model. '
                'The most common cause of this is the addition or removal of '
                'atoms while a simulation is running (this should be avoided). '
                'The simulation has been terminated to avoid data corruption, '
                'and the display reverted to a default state.')
            from ..dialog import generic_warning
            generic_warning(msg)

    def _rama_a_sim_end_cb(self, *_):
        self.RamaAnnotator.track_whole_model = True
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def _rota_a_sim_end_cb(self, *_):
        self.rota_annotator.track_whole_model = True
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    _pr_update_reasons = frozenset((
        'target changed', 'enabled/disabled', 'spring constant changed'
    ))

    def _pr_changed_cb(self, trigger_name, changes):
        mgr, changes = changes
        change_types = set(changes.keys())
        from chimerax.atomic import concatenate
        changeds = []
        for reason in self._pr_update_reasons.intersection(change_types):
            changeds.append(changes[reason])
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

    _dr_update_reasons = frozenset(('target changed', 'enabled/disabled', 'spring constant changed'))

    def _dr_changed_cb(self, trigger_name, changes):
        mgr, changes = changes
        change_types = set(changes.keys())
        from chimerax.atomic import concatenate
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
        for reason in self._dr_update_reasons.intersection(change_types):
            changeds.append(changes[reason])
        if len(changeds):
            all_changeds = concatenate(changeds, remove_duplicates=True)
            all_changeds = all_changeds[all_changeds.sim_indices != -1]
            self.sim_handler.update_distance_restraints(all_changeds)

    def _dr_sim_end_cb(self, *_):
        restraints = self.distance_restraint_mgr.intra_restraints(self.sim_construct.all_atoms)
        restraints.clear_sim_indices()
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    _adr_update_reasons = frozenset((
        'target changed', 'enabled/disabled',
        'adaptive restraint constant changed',
        'cutoff changed', 'spring constant changed'
    ))

    def _adr_changed_cb(self, trigger_name, changes):
        mgr, changes = changes
        change_types = set(changes.keys())
        from chimerax.atomic import concatenate
        if 'created' in change_types:
            created = changes['created']
            created = created[created.sim_indices == -1]
            # add only restraints with both atoms in the sim
            all_atoms = self.sim_construct.all_atoms
            indices = numpy.array([all_atoms.indices(atoms) for atoms in created.atoms])
            created = created[numpy.all(indices != -1, axis=0)]
            self.sim_handler.add_adaptive_distance_restraints(created)
        changeds = []
        for reason in self._adr_update_reasons.intersection(change_types):
            changeds.append(changes[reason])
        if len(changeds):
            all_changeds = concatenate(changeds, remove_duplicates=True)
            all_changeds = all_changeds[all_changeds.sim_indices != 1]
            self.sim_handler.update_adaptive_distance_restraints(all_changeds)

    def _adr_sim_end_cb(self, *_):
        for adrm in self.adaptive_distance_restraint_mgrs:
            restraints = adrm.intra_restraints(self.sim_construct.all_atoms)
            restraints.clear_sim_indices()
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def _dihe_r_changed_cb(self, trigger_name, changes):
        '''Used for proper dihedral and chiral restraints.'''
        mgr, changes = changes
        change_types = list(changes.keys())
        from chimerax.atomic import concatenate
        changeds = []
        if 'created' in change_types:
            # avoid double counting
            created = changes['created']
            created = created[created.sim_indices == -1]
            # add only restraints with all atoms in the sim
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


    def _dihe_r_sim_end_cb(self, *_):
        ''' Used for proper dihedral and chiral restraints.'''
        sc = self.sim_construct
        pdrs = self.proper_dihedral_restraint_mgr.get_all_restraints_for_residues(sc.mobile_residues)
        pdrs.clear_sim_indices()
        crs = self.chiral_restraint_mgr.get_restraints_by_atoms(sc.mobile_atoms)
        crs.clear_sim_indices()
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def _adaptive_dihe_r_changed_cb(self, trigger_name, changes):
        '''Used for all forms of dihedral restraints.'''
        mgr, changes = changes
        change_types = list(changes.keys())
        from chimerax.atomic import concatenate
        changeds = []
        if 'created' in change_types:
            # avoid double counting
            created = changes['created']
            created = created[created.sim_indices == -1]
            # add only restraints with all atoms in the sim
            all_atoms = self.sim_construct.all_atoms
            indices = numpy.array([all_atoms.indices(atoms) for atoms in created.atoms])
            created = created[numpy.all(indices != -1, axis=0)]
            self.sim_handler.add_adaptive_dihedral_restraints(created)
        if 'target changed' in change_types:
            changeds.append(changes['target changed'])
        if 'enabled/disabled' in change_types:
            changeds.append(changes['enabled/disabled'])
        if 'spring constant changed' in change_types:
            changeds.append(changes['spring constant changed'])
        if 'adaptive restraint constant changed' in change_types:
            changeds.append(changes['adaptive restraint constant changed'])
        if len(changeds):
            all_changeds = concatenate(changeds, remove_duplicates=True)
            all_changeds = all_changeds[all_changeds.sim_indices != -1]
            self.sim_handler.update_adaptive_dihedral_restraints(all_changeds)


    def _adaptive_dihe_r_sim_end_cb(self, *_):
        ''' Used for all forms of dihedral restraints.'''
        sc = self.sim_construct
        apdrs = self.adaptive_dihedral_restraint_mgr.get_all_restraints_for_residues(sc.mobile_residues)
        apdrs.clear_sim_indices()
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER



    def _tug_changed_cb(self, trigger_name, changes):
        mgr, changes = changes
        change_types = list(changes.keys())
        from chimerax.atomic import concatenate
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

    def _mdff_global_k_change_cb(self, trigger_name, data):
        mgr, k = data
        self.sim_handler.set_mdff_global_k(mgr.volume, k)

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
    def __init__(self, session, sim_params, sim_construct, forcefield_mgr):
        '''
        Prepares the simulation topology parameters and construct, and
        initialises the necessary Force objects to handle restraints. Most
        restraint forces must be populated using e.g.
        :func:`add_dihedral_restraints` before initialising the context and
        beginning the simulation. While it is possible to add new restraints
        *after* the simulation has started, in general this is only advisable in
        situations where it is impossible or impractical to define all possible
        restraints in advance (since each addition requires a costly
        reinitialisation of the simulation context). For example, the
        performance penalty to having all possible position restraints
        pre-defined (and mostly disabled) is minimal, but it is not practical to
        pre-define distance restraints between all possible atom pairs.

        Args:
            * session:
                - the ChimeraX session object
            * sim_params:
                - a :py:class:`SimParams` instance
            * sim_construct:
                - a :py:class:`Sim_Construct` instance
            * forcefield_mgr:
                - a class that behaves as a
                  {name: :py:class:`OpenMM::ForceField`} dict.

        '''
        self.session = session
        self._params = sim_params
        self._sim_construct = sim_construct

        self._thread_handler = None

        self._paused = False
        self._sim_running = False
        self._unstable = True
        self._unstable_counter = 0

        atoms = self._atoms = sim_construct.all_atoms
        # Forcefield used in this simulation
#        from .forcefields import forcefields
        ff = forcefield_mgr[sim_params.forcefield]
        ligand_db = forcefield_mgr.ligand_db(sim_params.forcefield)
#        ff = self._forcefield = self.define_forcefield(forcefields[sim_params.forcefield])

        # All custom forces in the simulation
        self.all_forces = []

        # A Volume: LinearInterpMapForce dict covering all MDFF forces
        self.mdff_forces = {}

        logger = self.session.logger
        # Overall simulation topology
        template_dict = find_residue_templates(sim_construct.all_residues, ff, ligand_db=ligand_db, logger=session.logger)
        top, residue_templates = create_openmm_topology(atoms, template_dict)
        self.topology = top


        self._temperature = sim_params.temperature

        trigger_names = (
            'sim started',
            'clash detected',
            'coord update',
            'sim paused',
            'sim resumed',
            'sim terminated',
        )
        from chimerax.core.triggerset import TriggerSet
        t = self._triggers = TriggerSet()
        for name in trigger_names:
            t.add_trigger(name)

        system = self._system = self._create_openmm_system(ff, top,
            sim_params, residue_templates)



        self.set_fixed_atoms(sim_construct.fixed_atoms)
        self._thread_handler = None
        # CustomExternalForce handling mouse and haptic interactions
        self._tugging_force = None

        self._force_update_pending = False
        self._coord_update_pending = False
        self._coord_push_handler = None

        self._context_reinit_pending = False
        self._minimize = False
        self._current_mode = 'min' # One of "min" or "equil"


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
        t = self._simulation.integrator.getTemperature()
        return t.value_in_unit(defaults.OPENMM_TEMPERATURE_UNIT)

    @temperature.setter
    def temperature(self, temperature):
        self._simulation.integrator.setTemperature(temperature)

    @property
    def smoothing(self):
        if self.thread_handler is not None:
            return self.thread_handler.smoothing
        return self._params.trajectory_smoothing

    @smoothing.setter
    def smoothing(self, flag):
        if self.thread_handler is not None:
            self.thread_handler.smoothing = flag

    smoothing.__doc__ = OpenMM_Thread_Handler.smoothing.__doc__

    @property
    def smoothing_alpha(self):
        if self.thread_handler is not None:
            return self.thread_handler.smoothing_alpha

    @smoothing_alpha.setter
    def smoothing_alpha(self, alpha):
        if self.thread_handler is not None:
            self.thread_handler.smoothing_alpha = alpha

    smoothing_alpha.__doc__ = OpenMM_Thread_Handler.smoothing_alpha.__doc__

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
        residue_to_template, ambiguous, unassigned = forcefield.assignTemplates(
            top, ignoreExternalBonds=True, explicit_templates=residue_templates
        )
        if len(ambiguous) or len(unassigned):
            raise RuntimeError('Unparameterised residue detected')

            # from chimerax.core.errors import UserError
            # err_text = ''
            # if len(ambiguous):
            #     err_text += "The following residues match multiple topologies: \n"
            #     for r, tlist in ambiguous.items():
            #         err_text += "{}{}: ({})\n".format(r.name, r.index, ', '.join([info[0].name for info in tlist]))
            # if len(unassigned):
            #     err_text += "The following residues did not match any template: ({})".format(
            #         ', '.join(['{}{}'.format(r.name, r.index) for r in unassigned])
            #     )
            # raise UserError(err_text)
            #
            # residue_templates = {r: t[0].name for r, t in residue_to_template.items()}

        system_params = {
            'nonbondedMethod':      params.nonbonded_cutoff_method,
            'nonbondedCutoff':      params.nonbonded_cutoff_distance,
            'constraints':          params.rigid_bonds,
            'rigidWater':           params.rigid_water,
            'removeCMMotion':       params.remove_c_of_m_motion,
            'ignoreExternalBonds':  True,
            'residueTemplates':     residue_templates,
        }
        sys = forcefield.createSystem(top, **system_params)
        # self._convert_to_soft_core_potentials(sys)
        return sys

    def _convert_to_soft_core_potentials(self, system):
        from simtk.openmm.openmm import NonbondedForce
        for f in system.getForces():
            if type(f) == NonbondedForce:
                break
        from .custom_forces import NonbondedSoftcoreForce
        sf = NonbondedSoftcoreForce()
        sf.setNonbondedMethod(f.getNonbondedMethod())
        sf.setCutoffDistance(f.getCutoffDistance())
        sf.setSwitchingDistance(f.getSwitchingDistance())
        for i in range(system.getNumParticles()):
            charge, sigma, epsilon = f.getParticleParameters(i)
            sf.addParticle([sigma, epsilon])
            f.setParticleParameters(i, charge, sigma, 0)
        for i in range(f.getNumExceptions()):
            p1, p2, cp, sig, eps = f.getExceptionParameters(i)
            sf.addExclusion(p1, p2)
        system.addForce(sf)


    def initialize_restraint_forces(self, amber_cmap=True, tugging=True, position_restraints=True,
        distance_restraints=True, adaptive_distance_restraints=True,
        dihedral_restraints=True, adaptive_dihedral_restraints=True):
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
            * adaptive_distance_restraints:
                - Add an "adaptive" distance restraints force (implemented as a
                  :py:class:`AdaptiveDistanceRestraintForce`)
            * dihedral_restraints:
                - Add dihedral restraints force (implemented as a
                  :py:class:`FlatBottomTorsionRestraintForce`)
            * adaptive_dihedral_restraints:
                - Add an "adaptive" dihedral restraints force (implemented as a
                  :py:class:`TopOutTorsionForce`)
        '''
        logger = self.session.logger
        params = self._params
        logger.status('Initialising forces...')
        if amber_cmap:
            self.initialize_amber_cmap_force()
        if tugging:
            self.initialize_tugging_force(params.restraint_max_force)
        if position_restraints:
            self.initialize_position_restraints_force(params.restraint_max_force)
        if distance_restraints:
            self.initialize_distance_restraints_force(params.restraint_max_force)
        if adaptive_distance_restraints:
            self.initialize_adaptive_distance_restraints_force()
        if dihedral_restraints:
            self.initialize_dihedral_restraints_force()
        if adaptive_dihedral_restraints:
            self.initialize_adaptive_dihedral_restraints_force()

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
        logger = self.session.logger
        params = self._params
        if params.use_gbsa:
            self.initialize_implicit_solvent(params)
        integrator = self._prepare_integrator(params)
        platform = openmm.Platform.getPlatformByName(params.platform)

        properties = {}
        device_index = params.device_index
        if device_index is not None:
            if params.platform=="CUDA":
                properties['CudaDeviceIndex']=str(device_index)
            elif params.platform=='OpenCL':
                properties['OpenCLDeviceIndex']=str(device_index)


        from simtk.openmm import app
        logger.status('Initialising primary simulation object')
        s = self._simulation = app.Simulation(self.topology, self._system,
            integrator, platform, properties)
        c = self._context = s.context
        c.setPositions(0.1*self._atoms.coords)
        c.setVelocitiesToTemperature(self.temperature)
        self._thread_handler = OpenMM_Thread_Handler(c, params)
        self.smoothing = params.trajectory_smoothing
        self.smoothing_alpha = params.smoothing_alpha
        logger.status('')

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
        self._startup = True
        self._startup_counter = 0
        self._minimize_and_go()

    def find_clashing_atoms(self, max_force = defaults.CLASH_FORCE):
        if not self._sim_running:
            raise RuntimeError('Simulation must be running first!')
        c = self._context
        state = c.getState(getForces=True)
        forces = state.getForces(asNumpy = True)
        import numpy
        force_mags = numpy.linalg.norm(forces, axis=1)
        sort_i = numpy.argsort(force_mags)[::-1]
        sorted_forces = numpy.sort(force_mags)[::-1]
        fmask = (sorted_forces > max_force)
        sc = self._sim_construct
        mmask = numpy.zeros(fmask.shape, numpy.bool)
        mmask[sc.all_atoms.indices(sc.mobile_atoms)] = True

        mask = numpy.logical_and(fmask, mmask[sort_i])
        clashes = sort_i[mask]
        return clashes, sorted_forces[mask]


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
        elif th.unstable() or self._unstable or not th.minimization_converged:
            f = th.minimize
            f_args = []
            final_args = [True]
        elif self.minimize:
            f = th.minimize
            f_args=[params.minimization_convergence_tol_end]
            final_args = [True]
        else:
            f = th.step
            f_args = (params.sim_steps_per_gui_update,)
            final_args = []
        delayed_reaction(self.session.triggers, 'new frame', f, f_args,
            th.thread_finished, self._update_coordinates_and_repeat, final_args)
        self._unstable = False

    def _resume(self):
        if self._force_update_pending:
            self._update_forces_in_context_if_needed()
        if self._coord_update_pending:
            self._push_coords_to_sim()
            self._coord_update_pending=False
        self._repeat_step()


    def _update_coordinates_and_repeat(self, reinit_vels = False):
        if self._startup:
            self._startup_counter += 1
            if self._startup_counter >= self._params.simulation_startup_rounds:
                self._startup = False
        th = self.thread_handler
        try:
            self.atoms.coords = th.coords
        except ValueError:
            self.stop(reason='coord length mismatch')
        self.triggers.activate_trigger('coord update', None)
        if th.clashing:
            self._unstable = True
            if not self._startup:
                if self._unstable_counter >= self._params.maximum_unstable_rounds:
                    self._unstable_counter = 0
                    self._pause=True
                    self.session.triggers.add_handler('frame drawn',
                        self._fire_sim_paused_trigger)
                    self.triggers.activate_trigger('clash detected', self.find_clashing_atoms())
                    return
            self._unstable_counter += 1
        else:
            self._unstable_counter = 0
        if self._force_update_pending:
            self._update_forces_in_context_if_needed()
        if reinit_vels:
            th.reinitialize_velocities()
        if self._stop:
            self._thread_handler.delete()
            self._thread_handler = None
            self._simulation = None
            self._sim_running = False
            self.triggers.activate_trigger('sim terminated', self._stop_reason)
            return
        if not self._pause:
            self._repeat_step()
        else:
            # Need to delay the trigger firing until after the frame is drawn,
            # so that the last coordinate update doesn't falsely register as
            # changes during pause.
            self.session.triggers.add_handler('frame drawn',
                self._fire_sim_paused_trigger)


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
        if self.pause:
            self._coord_update_pending=True
        else:
            if self._coord_push_handler is None:
                self._coord_push_handler = self.triggers.add_handler('coord update', self._push_coords_to_sim)

    def _push_coords_to_sim(self, *_):
        if self._pending_coords is not None:
            self.thread_handler.coords = self._pending_coords
        self._pending_coords = None
        self._coord_push_handler = None
        self._unstable = True
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
                pass
                # Trigger firing delayed until *after* the last coordinate
                # update has been applied
                # self.triggers.activate_trigger('sim paused', None)
            else:
                self.triggers.activate_trigger('sim resumed', None)
                self._resume()

    def _fire_sim_paused_trigger(self, *_):
        self.triggers.activate_trigger('sim paused', None)
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def stop(self, reason=None):
        '''
        Stop the simulation. This will destroy the
        :cpp:class:`OpenMM_Thread_Handler` object, rendering the Python class
        unusable.
        '''
        self._stop = True
        self._stop_reason = reason
        if self.pause:
            self._thread_handler.delete()
            self._thread_handler = None
            self._simulation = None
            self._sim_running = False
            self.triggers.activate_trigger('sim terminated', reason)

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
        phi_filter = numpy.all(phi_indices!=-1, axis=1)
        psi_indices = numpy.column_stack([sc.indices(atoms) for atoms in psi_atoms])
        psi_filter = numpy.all(psi_indices!=-1, axis=1)
        combined_filter = numpy.logical_and(phi_filter, psi_filter)
        resnames = resnames[combined_filter]
        phi_indices = phi_indices[combined_filter]
        psi_indices = psi_indices[combined_filter]
        cf.add_torsions(resnames, phi_indices, psi_indices)

    ####
    # Dihedral restraints
    ####

        ##
        # Before simulation starts
        ##

    def initialize_dihedral_restraints_force(self):
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
                - Either a :py:class:`ProperDihedralRestraints` or a
                  :py:class:`ChiralRestraints`.
        '''
        force = self._dihedral_restraint_force
        all_atoms = self._atoms
        dihedral_atoms = restraints.dihedrals.atoms
        atom_indices = numpy.array([all_atoms.indices(atoms) for atoms in dihedral_atoms])
        # Filter out those which don't have all atoms in the simulation
        ifilter = numpy.all(atom_indices!=-1, axis=0)
        atom_indices = [a[ifilter] for a in atom_indices]
        #atom_indices = atom_indices[ifilter]
        restraints = restraints[ifilter]
        restraints.sim_indices = force.add_torsions(atom_indices,
            restraints.enableds, restraints.spring_constants, restraints.targets, restraints.cutoffs)
        self.context_reinit_needed()

    def add_dihedral_restraint(self, restraint):
        '''
        Singular form of :func:`add_dihedral_restraints`.

        Args:
            * restraint:
                - A :py:class:`ProperDihedralRestraint` or
                  :py:class:`ChiralRestraint` instance.
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
                - A :py:class:`ProperDihedralRestraints` instance. Any
                  restraints that are not part of the current simulation will
                  be ignored.
        '''
        force = self._dihedral_restraint_force
        restraints = restraints[restraints.sim_indices != -1]
        force.update_targets(restraints.sim_indices,
            restraints.enableds, restraints.spring_constants, restraints.targets, restraints.cutoffs)
        self.force_update_needed()

    def update_dihedral_restraint(self, restraint):
        '''
        Update the simulation to match the parameters (target angles,
        cutoffs, spring constants and enabled states) of the given restraint.

        Args:
            * restraint:
                - A :py:class:`ProperDihedralRestraint` instance. If the
                  restraint is not part of the current simulation it will be
                  ignored.
        '''
        force = self._dihedral_restraint_force
        force.update_target(restraint.sim_index,
            enabled=restraint.enabled, k=restraint.spring_constant,
            target=restraint.target, cutoff=restraint.cutoff)
        self.force_update_needed()

    #####
    # Adaptive Dihedral Restraints
    #####

    def initialize_adaptive_dihedral_restraints_force(self):
        '''
        Just initialise the restraint force. Must be called before the
        simulation starts, and before any restraints are added.
        '''
        from .custom_forces import TopOutTorsionForce
        df = self._adaptive_dihedral_restraint_force = TopOutTorsionForce()
        self._system.addForce(df)
        self.all_forces.append(df)

    def add_adaptive_dihedral_restraints(self, restraints):
        '''
        Add a set of adaptive dihedral restraints to the simulation. This sets the
        :attr:`sim_index` for each restraint so it knows its own place in the
        simulation for later updating purposes. Automatically calls
        :func:`context_reinit_needed()`.

        Args:
            * restraints:
                - A :py:class:`AdaptiveDihedralRestraints` object.
        '''
        force = self._adaptive_dihedral_restraint_force
        all_atoms = self._atoms
        dihedral_atoms = restraints.dihedrals.atoms
        atom_indices = numpy.array([all_atoms.indices(atoms) for atoms in dihedral_atoms])
        # Filter out those which don't have all atoms in the simulation
        ifilter = numpy.all(atom_indices!=-1, axis=0)
        atom_indices = [a[ifilter] for a in atom_indices]
        #atom_indices = atom_indices[ifilter]
        restraints = restraints[ifilter]
        restraints.sim_indices = force.add_torsions(atom_indices,
            restraints.enableds, restraints.spring_constants, restraints.targets, restraints.kappas)
        self.context_reinit_needed()

    def add_adaptive_dihedral_restraint(self, restraint):
        '''
        Singular form of :func:`add_dihedral_restraints`.

        Args:
            * restraint:
                - A :py:class:`AdaptiveDihedralRestraint` instance
        '''
        force = self._adaptive_dihedral_restraint_force
        all_atoms = self._atoms
        dihedral_atoms = restraint.dihedral.atoms
        indices = [all_atoms.index(atom) for atom in dihedral_atoms]
        restraint.sim_index = force.add_torsion(*indices,
            (float(restraint.enabled), restraint.spring_constant, restraint.target, cos(restraint.kappa)))
        self.context_reinit_needed()

    def update_adaptive_dihedral_restraints(self, restraints):
        '''
        Update the simulation to match the parameters (target angles,
        kappas, spring constants and enabled states) of the given restraints.

        Args:
            * restraints:
                - A :py:class:`AdaptiveDihedralRestraints` instance. Any
                  restraints that are not part of the current simulation will
                  be ignored.
        '''
        force = self._adaptive_dihedral_restraint_force
        restraints = restraints[restraints.sim_indices != -1]
        force.update_targets(restraints.sim_indices,
            restraints.enableds, restraints.spring_constants, restraints.targets, restraints.kappas)
        self.force_update_needed()

    def update_adaptive_dihedral_restraint(self, restraint):
        '''
        Update the simulation to match the parameters (target angle,
        kappa, spring constant and enabled state) of the given restraint.

        Args:
            * restraint:
                - A :py:class:`AdaptiveDihedralRestraint` instance. If the
                  restraint is not part of the current simulation it will be
                  ignored.
        '''
        force = self._adaptive_dihedral_restraint_force
        force.update_target(restraint.sim_index,
            enabled=restraint.enabled, k=restraint.spring_constant,
            target=restraint.target, kappa=restraint.kappa)
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
                - a :py:class:`DistanceRestraints` instance
        '''
        force = self._distance_restraints_force
        all_atoms = self._atoms
        dr_atoms = restraints.atoms
        indices = numpy.array([all_atoms.indices(atoms) for atoms in dr_atoms])
        ifilter = numpy.all(indices!=-1, axis=0)
        indices = [i[ifilter] for i in indices]
        restraints = restraints[ifilter]
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
                - a :py:class:`DistanceRestraint` instance
        '''
        force = self._distance_restraints_force
        all_atoms = self._atoms
        dr_atoms = restraint.atoms
        indices = [all_atoms.index(atom) for atom in dr_atoms]
        if -1 in indices:
            raise TypeError('At least one atom in this restraint is not in the simulation!')
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
                - a :py:class:`DistanceRestraints` instance
        '''
        force = self._distance_restraints_force
        restraints = restraints[restraints.sim_indices != -1]
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
                - a :py:class:`DistanceRestraint` instance
        '''
        force = self._distance_restraints_force
        force.update_target(restraint.sim_index,
            restraint.enabled, k=restraint.spring_constant, target=restraint.target/10)
        self.force_update_needed()

    ####
    # Adaptive Distance Restraints
    ####

        ##
        # Before simulation starts
        ##

    def initialize_adaptive_distance_restraints_force(self):
        '''
        Just initialise the force. Unlike ISOLDE's other custom forces, this
        does not take a max_force argument. Given that the forces applied by
        this method should be quite weak under almost all circumstances, this
        should not be a problem.
        '''
        from .custom_forces import AdaptiveDistanceRestraintForce
        f = self._adaptive_distance_restraints_force = AdaptiveDistanceRestraintForce()
        self._system.addForce(f)
        self.all_forces.append(f)
        return f

    def add_adaptive_distance_restraints(self, restraints):
        '''
        Add a set of adaptive distance restraints to the simulation. Sets
        :attr:`sim_index` for each restraint so it knows its place in the
        simulation. Automatically calls :func:`context_reinit_needed`.

        Args:
            * restraints:
                - a :py:class:`AdaptiveDistanceRestraints` instance
        '''
        force = self._adaptive_distance_restraints_force
        all_atoms = self._atoms
        dr_atoms = restraints.atoms
        indices = numpy.array([all_atoms.indices(atoms) for atoms in dr_atoms])
        ifilter = numpy.all(indices!=-1, axis=0)
        indices = [i[ifilter] for i in indices]
        restraints = restraints[ifilter]
        restraints.sim_indices = force.add_bonds(indices,
            restraints.enableds, restraints.kappas, restraints.cs/10,
            restraints.targets/10, restraints.tolerances/10, restraints.alphas
        )
        self.context_reinit_needed()

    def add_adaptive_distance_restraint(self, restraint):
        '''
        Add a single adaptive distance restraint to the simulation. Sets
        :attr:`sim_index` for the restraint so it knows its place in the
        simulation. Automatically calls :func:`context_reinit_needed`

        Args:
            * restraint:
                - a :py:class:`AdaptiveDistanceRestraint` instance
        '''
        force = self._adaptive_distance_restraints_force
        all_atoms = self._atoms
        dr_atoms = restraint.atoms
        indices = [all_atoms.index(atom) for atom in dr_atoms]
        if -1 in indices:
            raise TypeError('At least one atom in this restraint is not in the simulation!')
        restraint.sim_index = force.addBond(*indices,
            (float(restraint.enabled), restraint.kappa, restraint.c/10,
            restraint.target/10, restraint.tolerance/10, restraint.alpha))
        self.context_reinit_needed()

        ##
        # During simulation
        ##

    def update_adaptive_distance_restraints(self, restraints):
        '''
        Update the simulation to reflect the current parameters of the given
        restraints.

        Args:
            * restraints:
                - a :py:class:`AdaptiveDistanceRestraints` instance
        '''
        force = self._adaptive_distance_restraints_force
        restraints = restraints[restraints.sim_indices != -1]
        force.update_targets(restraints.sim_indices,
            restraints.enableds, restraints.kappas, restraints.cs/10,
            restraints.targets/10, restraints.tolerances/10, restraints.alphas)
        self.force_update_needed()

    def update_adaptive_distance_restraint(self, restraint):
        '''
        Update the simulation to reflect the current parameters of the given
        restraint.

        Args:
            * restraint:
                - a :py:class:`DistanceRestraint` instance
        '''
        force = self._adaptive_distance_restraints_force
        force.update_target(restraint.sim_index,
            restraint.enabled, restraint.kappa, restraint.c/10,
            restraint.target/10, restraint.tolerance/10, restraint.alpha)
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
                - a :py:class:`PositionRestraints` instance
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
                - a :py:class:`PositionRestraint` instance
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
                - a :py:class:`PositionRestraints` instance
        '''
        force = self._position_restraints_force
        restraints = restraints[restraints.sim_indices != -1]
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
                - a :py:class:`PositionRestraint` instance
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
                - a :py:class:`TuggableAtoms` instance
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
                - a :py:class:`TuggableAtom` instance
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
                - a :py:class:`TuggableAtoms` instance
        '''
        force = self._tugging_force
        tuggables = tuggables[tuggables.sim_indices !=-1]
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
                - a :py:class:`TuggableAtom` instance
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
        from .custom_forces import (CubicInterpMapForce,
            CubicInterpMapForce_Low_Memory
            )
        v = volume
        region = list(v.region)
        # Ensure that the region ijk step size is [1,1,1]
        region[-1] = [1,1,1]
        #v.new_region(ijk_min=region[0], ijk_max=region[1], ijk_step=[1,1,1])
        data = v.region_matrix(region=region)
        if any (dim < 3 for dim in data.shape):
            # Not enough of this map is covering the atoms. Leave it out.
            return
        if not data.data.c_contiguous:
            data_copy = numpy.empty(data.shape, numpy.float32)
            data_copy[:] = data
            data = data_copy
        if data.size < self._params.max_cubic_map_size:
            Map_Force = CubicInterpMapForce
        else:
            print("Map is too large for fast cubic interpolation on the GPU!"\
                  +" Switching to slower, more memory-efficient implementation.")
            Map_Force = CubicInterpMapForce_Low_Memory
        from chimerax.core.geometry import Place
        tf = v.data.xyz_to_ijk_transform
        # Shift the transform to the origin of the region
        region_tf = Place(axes=tf.axes(), origin = tf.origin() -
            v.data.xyz_to_ijk(v.region_origin_and_step(region)[0]))
        # In OpenMM forces, parameters can only be per-particle, or global to
        # the entire context. So if we want a parameter that's constant to all
        # particles in this force, it needs a unique name so it doesn't
        # interfere with other instances of the same force.
        suffix = str(len(self.mdff_forces)+1)
        f = Map_Force(data, region_tf.matrix, suffix, units='angstroms')
        f.setForceGroup(1)
        self.all_forces.append(f)
        self._system.addForce(f)
        self.mdff_forces[v] = f

    def add_mdff_atoms(self, mdff_atoms, volume):
        '''
        Add a set of MDFF atom proxies to the force associated with the given
        volume. Automatically calls :func:`context_reinit_needed`

        Args:
            * mdff_atoms:
                - a :py:class:`MDFFAtoms` instance
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
                - a :py:class:`MDFFAtom` instance
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
        the given volume. NOTE: this will trigger a reinitialisation of the
        simulation context, so use sparingly!

        Args:
            * volume:
                - the :py:class:`chimerax.Volume` instance that was used to
                  create the target force
            * k:
                - the new coupling constant, in
                  :math:`kJ mol^{-1} (\\text{map density unit})^{-1} nm^3`
        '''
        f = self.mdff_forces[volume]
        if self.sim_running:
            context = self.thread_handler.context
            f.set_global_k(k, context = context)
        else:
            f.set_global_k(k)
            self.context_reinit_needed()

    def update_mdff_transform(self, volume):
        f = self.mdff_forces[volume]
        region = volume.region
        tf = volume.data.xyz_to_ijk_transform
        from chimerax.core.geometry import Place
        region_tf = Place(axes=tf.axes(), origin = tf.origin() -
            volume.data.xyz_to_ijk(volume.region_origin_and_step(region)[0]))

        if self.sim_running:
            context = self.thread_handler.context
            f.update_transform(region_tf.matrix, context=context)


    def set_mdff_scale_factor(self, volume, scale_factor):
        '''
        Adjust the dimensions of the mdff map in the simulation. This is
        typically used to optimize the scaling in the course of a single
        simulation. The final scale should be applied back to the original map,
        so that in future simulations the scale factor is 1.0.

        Args:
            * volume:
                - the :py:class:`chimerax.Volume` instance that was used to
                  create the target force
            * scale_factor:
                - the new scale factor (dimensionless). Changing the scale
                  factor by more than 1-2% in a single go is dangerous!
        '''
        f = self.mdff_forces[volume]
        if self.sim_running:
            context = self.thread_handler.context
            f.set_map_scale_factor(scale_factor, context=context)
        else:
            f.set_map_scale_factor(scale_factor)
            self.context_reinit_needed()

    def update_mdff_atoms(self, mdff_atoms, volume):
        '''
        Update the simulation to reflect the new parameters (individual
        coupling constants, enabled/disabled states) for the given MDFF atom
        proxies.

        Args:
            * mdff_atoms:
                - a :py:class:`MDFFAtoms` instance
            * volume:
                - the :py:class:`chimerax.Volume` instance that was used to
                  create the target force.
        '''
        f = self.mdff_forces[volume]
        mdff_atoms = mdff_atoms[mdff_atoms.sim_indices != -1]
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
                - a :py:class:`MDFFAtom` instance
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


    def initialize_implicit_solvent(self, params):
        '''
        Add a Generalised Born Implicit Solvent (GBIS) formulation.

        Args:
            * params:
                - a :py:class:`SimParams` instance
        '''

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
        n = f.getNumParticles()
        params = numpy.empty((n, 3))

        for i in range(f.getNumParticles()):
            params[i,0] = f.getParticleParameters(i)[0].value_in_unit(
                defaults.OPENMM_CHARGE_UNIT
            )
        from .custom_forces import GBSAForce
        gbforce = self._gbsa_force = GBSAForce(**gbsa_params)
        params[:,1:] = gbforce.getStandardParameters(top)
        gbforce.addParticles(params)
        gbforce.finalize()
        system.addForce(gbforce)

class Map_Scale_Optimizer:
    '''
    Tool to optimise the magnification of the map(s) used for MDFF, by finding
    the magnification factor that minimises the conformational energy of the
    atomic model after settling. The reference point for scaling (that is, the
    coordinate that models and map shrink towards or grow away from) will be
    the origin of the *largest* map.
    '''
    def __init__(self, sim_manager, save_timepoints=False):
        import numpy
        self._save = save_timepoints
        sm = self.sim_manager = sim_manager
        self.isolde = sm.isolde
        sh = self.sim_handler = sm.sim_handler
        context = self.context = sh._context
        volumes = self._volumes = list(sh.mdff_forces.keys())
        self._original_scales = {v: numpy.array(v.data.step) for v in volumes}
        self._original_origins = {v: numpy.array(v.data.origin) for v in volumes}
        largest_volume = sorted(volumes, key=lambda v: sum(v.data.size))[-1]
        self._reference_coord = numpy.array(largest_volume.data.origin)

    def optimize_map_scales(self, tolerance=1e-4, step_size=0.01, num_steps=11,
            loops=1):
        '''
        Optimise the dimensional scaling of the MDFF map(s), by minimising the
        conformational energy as a function of scaling.

        Args:
            * tolerance:
                - convergence tolerance at each step (kJ/mol/atom)
            * step_size:
                - fractional change in magnification at each step
            * num_steps:
                - number of steps in the line search *in each direction*
            * loops:
                - number of times to loop over the full search. In general this
                  should be left to 1, and only increased if you wish to check
                  reproducibility of the energies.
        '''
        import numpy
        self.step_size=step_size
        sm = self.sim_manager
        sh = self.sim_handler
        s = self.structure = sm.isolde.selected_model
        self.tolerance = tolerance * len(s.atoms)

        self.context = sh._context
        self._original_temperature = sh.temperature
        self._original_minimize = sh.minimize

        self._search_index = 0
        self._define_search_strategy(step_size, num_steps, loops)
        # self._line_search_down = numpy.linspace(1,1-(num_steps-1)*step_size,num=num_steps)
        # self._line_search_up = numpy.linspace(1,1+(num_steps-1)*step_size, num=num_steps)
        # self._search_down_vals = []
        # self._search_up_vals = []
        e = self._last_energy = self.get_current_energy()

        from ..delayed_reaction import delayed_reaction
        delayed_reaction(sh.triggers, 'coord update', self._start, [], self._energy_converged,
            self._initial_min_and_start, [])

    def _define_search_strategy(self, step_size, num_steps, loops=1):
        '''
        Starting from the initial coordinates and magnification, a "loop" is
        defined as a stepwise line search:

            * increasing magnification (num_steps steps of step_size)
            * decreasing magnification (2*num_steps steps of -step_size)
            * increasing magnification back to the original
              (num_steps steps of step_size)
        '''
        import numpy
        max_mag = 1+(num_steps-1)*step_size
        min_mag = 1-(num_steps-1)*step_size
        steps_up = numpy.linspace(1, max_mag, num=num_steps)
        steps_down = numpy.linspace(max_mag-step_size, min_mag, num=num_steps*2)
        steps_ret = numpy.linspace(min_mag+step_size, 1, num=num_steps-1)
        sv = self._search_vals = numpy.concatenate((steps_up, steps_down, steps_ret))
        sv = self._search_vals = numpy.concatenate([sv]*loops)
        self._energies = numpy.ones(sv.shape, numpy.double)*numpy.nan

    def _start(self):
        sh = self.sim_handler
        sh.minimize=False
        sh.temperature = 0
        # First do a *really* thorough minimisation at the starting magnification
        # to set a good reference point.
        sh.thread_handler.minimize(1e-4, max_iterations=100000)
        sh.pause=False

    def _revert_sim_conditions(self):
        sh = self.sim_handler
        sh.minimize = self._original_minimize
        sh.temperature = self._original_temperature

    def get_current_energy(self):
        from simtk.openmm import unit
        state = self.context.getState(getEnergy=True, groups=set([0]))
        return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    def _energy_converged(self):
        current_e = self.get_current_energy()
        if abs(self._last_energy - current_e) < self.tolerance:
            self._last_energy = current_e
            return True
        self._last_energy = current_e
        return False

    def _initial_min_and_start(self):
        from ..checkpoint import CheckPoint
        self._reference_checkpoint = CheckPoint(self.sim_manager.isolde)
        self._initial_model_center = self.structure.atoms.coords.mean(axis=0)
        self._line_search(self._search_index)

    def _update_all_scales(self, scale):
        import numpy
        sh = self.sim_handler
        rc = self._reference_coord
        for v in self._volumes:
            v.data.set_step(self._original_scales[v]*scale)
            o = self._original_origins[v]
            v.data.set_origin(rc + scale*(o-rc))
            v.update_drawings()
            sh.update_mdff_transform(v)
            # sh.set_mdff_scale_factor(v, scale)
        atoms = self.sim_manager.sim_construct.all_atoms
        c = atoms.coords.mean(axis=0)
        oc = self._initial_model_center
        nc = rc+scale*(oc-rc)
        atoms.coords += (nc-c)
        self.sim_handler.thread_handler.coords = atoms.coords


    def _line_search(self, index):
        search_vals = self._search_vals
        energies = self._energies
        if index >=len(search_vals):
            self._finalize()
            return
        energies[index] = self.get_current_energy()
        scale = search_vals[index]
        if self._save:
            from chimerax.core.commands import save
            filename = 'map_scale_optimize_{:.2f}.pdb'.format(scale)
            save.save(self.sim_manager.isolde.session, filename, [self.sim_manager.isolde.selected_model])
        sh = self.sim_handler
        from ..delayed_reaction import delayed_reaction
        delayed_reaction(sh.triggers, 'coord update',
            self._update_all_scales, [scale],
            self._energy_converged,
            self._line_search, [index+1])

    def _go_smoothly_to_scale(self, current_scale, target_scale):
        from ..delayed_reaction import delayed_reaction
        sh = self.sim_handler
        diff = current_scale-target_scale
        if abs(diff) < 0.01:
            delayed_reaction(sh.triggers, 'coord update',
                self._update_all_scales, [target_scale],
                self._energy_converged,
                self._revert_sim_conditions, [])
            return
        else:
            if diff <0:
                new_scale = current_scale + self.step_size
            else:
                new_scale = current_scale - self.step_size
            delayed_reaction(sh.triggers, 'coord update',
                self._update_all_scales, [new_scale],
                self._energy_converged,
                self._go_smoothly_to_scale, [new_scale, target_scale])


    def _finalize(self):
        import numpy
        x_vals = self._final_scales = self._search_vals
        y_vals = self._final_energies = self._energies
        # x_vals = self._final_scales = list(self._line_search_down[::-1])+list(self._line_search_up[1:])
        # y_vals = self._final_energies = numpy.array(self._search_down_vals[::-1]+self._search_up_vals[1:])
        # from scipy.interpolate import UnivariateSpline
        # spline = self.spline = UnivariateSpline(x_vals, y_vals)
        # from scipy.optimize import minimize
        # final_scale = minimize(spline, 1, bounds=[[min(x_vals), max(x_vals)]]).x[0]
        # min_index = numpy.argmin(y_vals)
        # Fit a quadratic to the data, and select its minimum as the optimum scale
        coeffs = self.quadratic_coeffs = numpy.polyfit(x_vals, y_vals, 2)
        final_scale = -coeffs[1]/(2*coeffs[0])
        # final_scale = x_vals[min_index]
        from ..delayed_reaction import delayed_reaction
        self._update_all_scales(1.0)
        self._reference_checkpoint.revert()
        print("Final map scale factor: {}".format(final_scale))
        self._go_smoothly_to_scale(1.0, final_scale)
        # self._update_all_scales(final_scale)
        # self._revert_sim_conditions()









class Sim_Performance_Tracker:
    c = 0
    def __init__(self, sim_handler, report_interval=50, log_func=None):
        self._ri = report_interval
        self._sh = sim_handler
        self._params = sim_handler._params
        if log_func is not None:
            self._log_func = log_func
        else:
            self._log_func = print
        self._running = False
    def start(self):
        if self._running:
            return
        from time import time
        self._start_time = time()
        self._h = self._sh.triggers.add_handler('coord update', self._update)
        self._running = True

    @property
    def report_interval(self):
        return self._ri

    @report_interval.setter
    def report_interval(self, interval):
        self._ri = interval

    def _update(self, *_):
        self.c = (self.c+1)%self._ri
        if self.c == 0:
            from time import time
            interval = (time()-self._start_time)/self._ri
            self._log_func('Average time per coord update: {:.2f} ms ({:.2f} time steps/s)'.format(interval*1000, 1/interval*self._params.sim_steps_per_gui_update))
            self._start_time = time()
    def stop(self):
        if self._h is not None:
            self._sh.triggers.remove_handler(self._h)
            self._h=None
        self._running = False

def create_openmm_topology(atoms, residue_templates):
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

    m = atoms.unique_structures[0]
    session = m.session
    anames   = atoms.names
    n = len(anames)
    enames   = atoms.element_names
    residues = atoms.residues
    rnames   = residues.names
    rindices = m.residues.indices(residues)
    rnums    = residues.numbers
    insertion_codes = residues.insertion_codes
    cids    = atoms.chain_ids
    bonds = atoms.intra_bonds
    bond_is = [atoms.indices(alist) for alist in bonds.atoms]

    #template_indices = list(residue_templates.keys())
    templates_out = {}
    from simtk.openmm.app import Topology, Element
    top = topology = Topology()
    cmap = {}
    rmap = {}
    unique_check = {}
    atoms = {}
    rcount = 0
    for i, (aname, ename, rname, rindex, rnum, icode, cid) in enumerate(
            zip(anames, enames, rnames, rindices, rnums, insertion_codes, cids)):
        if not cid in cmap:
            cmap[cid] = top.addChain()   # OpenMM chains have no name
        if not rindex in rmap:
            rid = (rnum, icode, cid)
            if rid in unique_check:
                session.logger.warning('Chain {}, residue {}{} specifies more '
                    'than one residue! The simulation can still run, but this '
                    'will probably cause problems later if not rectified by '
                    'renumbering.'.format(cid, rnum, icode))
            unique_check[rid]=rindex
            res = rmap[rindex] = top.addResidue(rname, cmap[cid])
            if rcount in residue_templates.keys():
                templates_out[res] = residue_templates[rcount]
            rcount += 1


        element = Element.getBySymbol(ename)
        atoms[i] = top.addAtom(aname, element,rmap[rindex])

    for i1, i2 in zip(*bond_is):
        top.addBond(atoms[i1],  atoms[i2])

    return top, templates_out


def find_residue_templates(residues, forcefield, ligand_db = None, logger=None):
    '''
    Works out the template name applicable to cysteine residues (since OpenMM
    can't work this out for itself when ignoreExternalBonds is True), and
    looks up the template name for all known parameterised ligands from the
    CCD.
    '''
    template_names = forcefield._templates.keys()
    import numpy
    templates = {}
    cys_indices = numpy.where(residues.names == 'CYS')[0]
    for c_i in cys_indices:
        rtype = cys_type(residues[c_i])
        if rtype is not None:
            templates[c_i] = rtype

    from chimerax.atomic import Residue
    ligands = residues[residues.polymer_types == Residue.PT_NONE]
    ligands = ligands[ligands.names != 'HOH']
    names = numpy.unique(ligands.names)
    from .amberff.glycam import (
        find_glycan_template_name_and_link,
        find_glycan_anchor_name,
        known_sugars
        )
    for name in names:
        user_name = 'USER_{}'.format(name)
        if user_name in template_names:
            indices = numpy.where(residues.names == name)[0]
            for i in indices:
                templates[i] = user_name
            continue

        if name in known_sugars:
            sugars = ligands[ligands.names == name]
            for sugar in sugars:
                tname, prot_res = find_glycan_template_name_and_link(sugar)
                if tname in template_names:
                    templates[residues.index(sugar)] = tname
                if prot_res is not None:
                    templates[residues.index(prot_res)] = find_glycan_anchor_name(prot_res)
            continue

        ccd_name = 'MC_{}'.format(name)
        if ccd_name in template_names:
            indices = numpy.where(residues.names == name)[0]
            for i in indices:
                templates[i] = ccd_name
            continue

        if ligand_db is not None:
            zip, namelist = ligand_db
            if name in namelist:
                from zipfile import ZipFile
                logger.info('Loading residue template for {} from internal database'.format(name))
                with ZipFile(zip) as zf:
                    with zf.open(name+'.xml') as xf:
                        forcefield.loadFile(xf)
                indices = numpy.where(residues.names==name)[0]
                for i in indices:
                    templates[i] = 'MC_{}'.format(name)
                continue

        # warn_str = ('No template with name {} found in the molecular dynamics '
        #     'forcefield. Attempting to match by topology. You may need to '
        #     'provide a custom ligand definition.')
        # logger.warning(warn_str.format(name))

    from .amberff.metal_name_map import metal_name_map
    from chimerax.atomic import Element
    metals = [n.upper() for n in Element.names if Element.get_element(n).is_metal]
    atoms = residues.atoms
    metal_atoms = atoms[numpy.in1d(atoms.names, metals)]
    for a in metal_atoms:
        r = a.residue
        if len(r.atoms) != 1:
            continue
        template_name = metal_name_map.get(r.name, None)
        if template_name is not None:
            templates[residues.index(r)] = template_name


    return templates

def cys_type(residue):
    from chimerax.atomic import Bonds, concatenate
    atoms = residue.atoms
    names = atoms.names
    rneighbors = residue.neighbors
    sulfur_atom = residue.find_atom('SG')
    if sulfur_atom is None:
        return None
    for r in rneighbors:
        if r.name in ('SF4', 'FES'):
            print('Found iron-sulfur cysteine')
            return 'MC_CYF'
    bonds = Bonds(sulfur_atom.bonds)
    if len(bonds) == 1:
        # Deprotonated
        return 'CYM'
    bonded_atoms = concatenate(bonds.atoms)
    for a in bonded_atoms:
        if a.residue != residue:
            if a.name == "SG":
                if 'OXT' in names:
                    return 'CCYX'
                if 'H1' in names:
                    return 'NCYX'
                return 'CYX'
            # Assume metal binding - will eventually need to do something better here
            return 'CYM'
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
