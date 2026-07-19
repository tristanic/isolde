# @Author: Tristan Croll <tic20>
# @Date:   26-Apr-2018
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 17-Sep-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



import numpy


from openmm import unit, openmm

from ..constants import defaults

DEFAULT_FORCE_GROUP=0
BOND_FORCE_GROUP=1
ANGLE_FORCE_GROUP=2
DIHEDRAL_FORCE_GROUP=3
NONBONDED_FORCE_GROUP=4

CORE_FORCE_GROUPS=set((DEFAULT_FORCE_GROUP, BOND_FORCE_GROUP, ANGLE_FORCE_GROUP, DIHEDRAL_FORCE_GROUP, NONBONDED_FORCE_GROUP))
MAP_FORCE_GROUP=5
RESTRAINT_FORCE_GROUP=6



from chimerax.isolde.delayed_reaction import delayed_reaction

from chimerax.isolde._openmm_async import OpenmmThreadHandler as _OpenmmThreadHandlerBase

class OpenmmThreadHandler(_OpenmmThreadHandlerBase):
    def __init__(self, context, params, num_real_atoms=0):
        # num_real_atoms is the count of leading particles backed by real
        # ChimeraX atoms. Any particles beyond it (crystallographic symmetry
        # copies represented as SymmetrySite virtual sites) are managed inside
        # the C++ handler and never exposed here. 0 means "all real" (the
        # default for simulations with no symmetry atoms).
        super().__init__(int(context.this), num_real_atoms)
        self.params = params
        self.context = context
        self._smoothing = params.trajectory_smoothing
        self._last_smooth = False
        self._last_mode = None

    @property
    def smoothing(self):
        '''
        If True, the displayed coordinates will be a smoothed average of the
        last set of equilibration steps. Note that for large values of
        sim_steps_per_gui_update this can lead to distorted geometry.
        '''
        return self._smoothing
    
    @smoothing.setter
    def smoothing(self, flag):
        self._smoothing = flag
    
    def step(self, steps):
        '''
        Advance the simulation integrator by the desired number of steps.
        Args:

            * steps:
                - an integer value
        '''
        super().step(steps)
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
        super().minimize(tolerance, max_iterations)
        self._last_mode = 'min'
    
    def reinitialize_velocities(self):
        '''
        Set the atomic velocities to random values consistent with the current
        temperature. Recommended after any energy minimisation.
        '''
        self.finalize_thread()
        c = self.context
        i = c.getIntegrator()
        c.setVelocitiesToTemperature(i.getIntegrator(0).getTemperature())
    
    @property
    def coords(self):
        '''
        Returns the coordinates of the atoms after the most recent thread
        completes. Can also be set, to push edited coordinates back to the
        simulation.
        '''
        if not self._smoothing or not self._last_smooth or self._last_mode !='equil':
            return self.current_coords
        else:
            return self.smoothed_coords

    @coords.setter
    def coords(self, coords):
        self.current_coords = coords
        self.reinitialize_velocities()
    
        



class SimConstruct:
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

        self.surroundings = model_atoms.subtract(all_atoms)
        if excluded_atoms is not None:
            self.surroundings = self.surroundings.subtract(excluded_atoms)

    @property
    def excluded_atoms(self):
        return self._excluded_atoms

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



            


class SimManager:
    '''
    Responsible for creating the :py:class:`SimHandler` and managing the
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
            * Creates the :py:class:`SimConstruct` object
            * Prepares the molecule visualisation for simulation (masking maps
              to the mobile selection, hiding atoms not in the simulation, etc.)
            * Prepares the MDFF managers (NOTE: this *must* be done after the
              preceding step, which ensures that each :py:class:`chimerax.Volume`
              has a region covering the mobile selection with sufficient padding)
            * Prepares all necessary callbacks to update the simulation when
              the parameters of restraints, mdff atom proxies etc. change.
            * Creates the :py:class:`SimHandler`
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
        # Enlarge the mobile selection to include the parents of any
        # crystallographic symmetry copies, so those copies are slaved to mobile
        # atoms (two-way coupling). Must precede the fixed-shell computation so
        # the shell wraps the newly-mobilised parent regions. Stashes the shell
        # on self._symmetry_shell for the SimHandler.
        mobile_atoms = self._add_symmetry_parents_to_mobile(mobile_atoms)
        from ..selections import get_shell_of_residues
        fixed_atoms = get_shell_of_residues(mobile_atoms.unique_residues,
            isolde_params.hard_shell_cutoff_distance).atoms
        self._prepare_validation_managers(mobile_atoms)
        self._prepare_restraint_managers()
        fixed_atoms = self._add_fixed_atoms_from_distance_restraints(mobile_atoms, fixed_atoms)
        if excluded_residues is not None:
            mobile_atoms, fixed_atoms, excluded_atoms = self._add_fixed_atoms_from_excluded_residues(mobile_atoms, fixed_atoms, excluded_residues)
        if len(mobile_atoms) == 0:
            from chimerax.core.errors import UserError
            raise UserError('Selection leads to a simulation with no mobile atoms!')

        sc = self.sim_construct = SimConstruct(model, mobile_atoms, fixed_atoms, excluded_atoms)
        # Hand the symmetry shell (if any) to the SimHandler via the construct.
        sc.symmetry_shell = self._symmetry_shell
        # The region-of-interest selection that drives map coverage + the
        # covered-multiplicity n_sym (the pre-promotion selection; None for
        # plain / whole-model sims, where the SimHandler falls back to the full
        # mobile set).
        sc.symmetry_coverage_atoms = self._symmetry_coverage_atoms

        logger.status('Preparing simulation handler')
        self._prepare_mdff_managers()
        sh = self.sim_handler = None
        uh = self._update_handlers = []
        try:
            sh = self.sim_handler = SimHandler(session, sim_params, sc,
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
            pt = self._performance_tracker = SimPerformanceTracker(self.sim_handler, report_interval, self.session.logger.status)
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

        :py:class:`SimManager` automatically saves a checkpoint when the
        simulation is started. This method allows the saving of an intermediate
        state - particularly useful before experimenting on an ambiguous
        region of your map. Each call to :func:`checkpoint` overwrites the
        previously-saved one. When ending a simulation you may choose to keep
        the final coordinates or revert to either of the last saved checkpoint
        or the start-of-simulation checkpoint.
        '''
        from ..checkpoint import CheckPoint
        self._current_checkpoint = CheckPoint(self.isolde)
        sh = self.sim_handler
        if sh is not None:
            sh._sim_steps_at_checkpoint = sh.sim_steps

    def revert_to_checkpoint(self):
        '''
        Reverts to the last saved checkpoint. If no checkpoint has been manually
        saved, reverts to the start of the simulation.
        '''
        self._current_checkpoint.revert()


    @property
    def rama_annotator(self):
        from chimerax.isolde import session_extensions as sx
        return sx.get_rama_annotator(self.model)
    
    @property
    def rota_annotator(self):
        from chimerax.isolde import session_extensions as sx
        return sx.get_rota_annotator(self.model)

    def _prepare_validation_managers(self, mobile_atoms):
        from .. import session_extensions as sx
        mobile_res = mobile_atoms.unique_residues
        self.rama_annotator.restrict_to_selected_residues(mobile_res)
        self.rota_annotator.restrict_to_selected_residues(mobile_res)

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
            from chimerax.isolde.restraints.restraint_utils import DEFAULT_ADAPTIVE_RESTRAINT_GROUP_NAME
            self.adaptive_distance_restraint_mgrs = [sx.get_adaptive_distance_restraint_mgr(m, DEFAULT_ADAPTIVE_RESTRAINT_GROUP_NAME)]

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
        # Map coverage tracks the user's region of interest, not the full mobile
        # set: in a symmetry-aware sim the mobile selection is enlarged with
        # distant copy-parent residues, but masking density around those far-away
        # parents is wasteful (memory + performance) and pointless (their density
        # restraint is delivered through their ghost inside this patch). Fall back
        # to the full mobile set for plain sims / whole-model sims.
        coverage_atoms = getattr(self, '_symmetry_coverage_atoms', None)
        if coverage_atoms is None:
            coverage_atoms = self.sim_construct.mobile_atoms
        sh.isolate_and_cover_selection(coverage_atoms,
            include_surrounding_residues = 0,
            show_context = isolde_params.hard_shell_cutoff_distance,
            mask_radius = isolde_params.map_mask_radius,
            extra_padding = 5,
            hide_surrounds = isolde_params.hide_surroundings_during_sim,
            focus = False)
        # When simulating with crystallographic symmetry atoms, extend the map
        # coverage to the symmetry copies' image positions too, so there is
        # density where the (now mobile) ghosts sit. The zone re-masks live as
        # the parents move.
        self._cover_maps_including_symmetry(sh, isolde_params, coverage_atoms)

    def _cover_maps_including_symmetry(self, sym_handler, isolde_params,
            coverage_atoms=None):
        '''
        Extend the map zone to cover the crystallographic symmetry copies that
        actually participate in the simulation, in addition to the real atoms of
        the region of interest. No-op when the simulation is not symmetry-aware
        (so ordinary sims keep their original, tight coverage exactly).

        Coverage is sourced directly from the **simulation shell**
        (:attr:`_symmetry_shell`), not from the Clipper display: the shell copies
        are the atoms that receive OpenMM virtual sites and MDFF force terms, so
        covering exactly them makes the covered map region coincide with the set
        of instances that carry a map term (see
        :func:`initialize_symmetry_copies` / :func:`add_mdff_atoms`, which key the
        covered-multiplicity ``n_sym`` and the term set off the same shell). This
        is also independent of the GUI (works headless).

        ``cover_atoms`` REPLACES the zone set by the preceding
        :func:`isolate_and_cover_selection`, so a single combined call passes the
        real region-of-interest atoms (identity transform) plus each shell copy's
        parent atom under its operator. Clipper masks around each atom's
        post-transform (ghost) position and re-masks live as the parents move, so
        density follows the ghosts (both are driven by the same master coords).
        The single fixed OpenMM MDFF box becomes the AABB of exactly these
        covered positions + padding.
        '''
        shell = getattr(self, '_symmetry_shell', None)
        if shell is None:
            return

        import numpy
        from chimerax.geometry import Places
        from chimerax.atomic import Atoms, concatenate

        # Cover the real atoms of the region of interest (the pre-promotion
        # selection when symmetry enlarged it), not the full mobile set - the
        # distant promoted parents are deliberately left uncovered.
        if coverage_atoms is None:
            coverage_atoms = self.sim_construct.mobile_atoms
        mobile = coverage_atoms[coverage_atoms.element_names != 'H']

        # Dense operator ids: transform 0 = identity (the real ROI atoms); 1..M
        # one per operator actually present among the shell copies (first-seen
        # order). Each copy's parent atom is covered under its operator, so the
        # mask lands at the ghost position R.r + t.
        op_to_place = {}
        place_mats = [numpy.eye(3, 4)]
        copy_parents, copy_place_idx = [], []
        for c in shell.copies:
            symop = c.symop_index
            pi = op_to_place.get(symop)
            if pi is None:
                pi = op_to_place[symop] = len(place_mats)
                place_mats.append(numpy.asarray(shell.symmats[symop],
                    dtype=numpy.float64))
            copy_parents.append(c.parent_atom)
            copy_place_idx.append(pi)

        if not copy_parents:
            # No copies -> nothing extra to cover (original behaviour).
            return

        transforms = Places(place_array=numpy.asarray(place_mats,
            dtype=numpy.float64))
        atoms = concatenate([mobile, Atoms(copy_parents)],
            remove_duplicates=False)
        transform_indices = numpy.concatenate([
            numpy.zeros(len(mobile), dtype=numpy.int32),
            numpy.asarray(copy_place_idx, dtype=numpy.int32)])

        sym_handler.map_mgr.cover_atoms(atoms, transforms=transforms,
            transform_indices=transform_indices,
            mask_radius=isolde_params.map_mask_radius, extra_padding=5)


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
        mgrs = []
        mgrs.append(self.distance_restraint_mgr)
        from chimerax.isolde.restraints import AdaptiveDistanceRestraintMgr
        for m in self.isolde.selected_model.child_models():
            if isinstance(m, AdaptiveDistanceRestraintMgr):
                mgrs.append(m)
        remainders = []
        for dr_m in mgrs:
            drs = dr_m.atoms_restraints(mobile_atoms)
            from chimerax.atomic import concatenate
            dr_atoms = concatenate(drs.atoms, remove_duplicates=True)
            remainder = dr_atoms.subtract(dr_atoms.intersect(mobile_atoms))
            remainder = remainder.subtract(remainder.intersect(fixed_atoms))
            remainders.append(remainder)
        remainder = concatenate(remainders)
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

    def _add_symmetry_parents_to_mobile(self, mobile_atoms):
        '''
        If crystallographic symmetry-aware simulation is enabled and available,
        work out which symmetry copies will participate (stashing the resulting
        shell on ``self._symmetry_shell`` for the SimHandler) and make sure every
        copy's *parent* atom is part of the **mobile** selection.

        A SymmetrySite copy is slaved to its parent, so for a copy to be a live,
        two-way participant (the whole point of the feature) its parent must be
        mobile - promoting parents as *fixed* would turn the copies into static
        walls. Parents whose residue isn't already mobile are therefore pulled
        into the mobile selection as whole residues. Because this enlarges the
        mobile set, it **must run before the fixed shell is computed** so that
        shell wraps the newly-mobilised parent regions.

        Returns the (possibly enlarged) mobile-atom selection. A no-op (returning
        ``mobile_atoms`` unchanged, with ``self._symmetry_shell = None``) when the
        feature is disabled, OpenMM lacks SymmetrySite, or the model has no
        crystallographic symmetry.
        '''
        self._symmetry_shell = None
        # The selection that should drive *map coverage* (and the symmetry
        # display): the user's region of interest, before it is enlarged with
        # distant copy-parent residues below. Covering the promoted parents at
        # their own (often far-away) crystallographic positions masks huge extra
        # map volumes for no benefit - a far parent's density restraint is
        # delivered through its ghost sitting inside this selection's patch, not
        # at the parent itself. None until/unless promotion actually enlarges the
        # selection (a plain sim then falls back to the full mobile set).
        self._symmetry_coverage_atoms = None
        sp = self.sim_params
        # Opt-in feature: default off leaves the simulation symmetry-blind
        # (bit-for-bit the original behaviour). See constants.SYMMETRY_AWARE.
        if not getattr(sp, 'symmetry_aware', False):
            return mobile_atoms
        from .symmetry_sim import symmetry_sims_available, build_symmetry_shell
        if not symmetry_sims_available():
            return mobile_atoms
        model = self.model
        ip = self.isolde_params
        # Two-regime shell radius: for a local simulation match the existing
        # fixed-shell depth (keeping the symmetry boundary at the same fidelity
        # as the rest of the truncated environment); for a whole-model
        # simulation, where the symmetry shell *is* the entire environment, go
        # cutoff-deep so it behaves like an infinite crystal.
        whole_model = (len(mobile_atoms) == model.num_atoms)
        if whole_model:
            cutoff = self._whole_model_symmetry_cutoff()
        else:
            cutoff = ip.hard_shell_cutoff_distance
        shell = build_symmetry_shell(model, mobile_atoms, cutoff,
            logger=self.session.logger)
        self._symmetry_shell = shell
        if shell is None:
            return mobile_atoms
        parents = shell.parent_atoms
        new_parents = parents.subtract(parents.intersect(mobile_atoms))
        if len(new_parents):
            # Remember the region of interest *before* enlargement, so map
            # coverage stays tight around what the user actually selected (plus
            # the ghosts drawn near it) rather than ballooning to the promoted
            # distant parents.
            self._symmetry_coverage_atoms = mobile_atoms
            new_res = new_parents.unique_residues
            mobile_atoms = mobile_atoms.merge(new_res.atoms)
            self.session.logger.info(f'ISOLDE: added {len(new_res)} residue(s) to '
                'the mobile selection as crystallographic symmetry-copy parents '
                '(so their symmetry mates are live, two-way participants).')
        return mobile_atoms

    def _whole_model_symmetry_cutoff(self):
        '''
        Cutoff-deep symmetry-shell radius (Angstroms) for whole-model sims:
        max(nonbonded, GBSA) cutoff plus a margin for atomic motion. Reads the
        live sim params (converting from OpenMM nm Quantities where needed) so it
        tracks the actual cutoffs in use rather than a hardcoded value.
        '''
        sp = self.sim_params
        def _in_angstroms(v):
            try:
                from openmm import unit
                return v.value_in_unit(unit.angstrom)
            except AttributeError:
                # Assume a plain value already in nm (OpenMM's native length).
                return float(v) * 10.0
        nb = _in_angstroms(sp.nonbonded_cutoff_distance)
        cutoff = nb
        if getattr(sp, 'use_gbsa', False):
            cutoff = max(cutoff, _in_angstroms(sp.nonbonded_cutoff_distance))
        # Margin for atomic displacement over the sim's lifetime (the copy set is
        # fixed once built).
        return cutoff + 3.0





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
            sel = selections.expand_selection_to_neighbors(core_atoms,
                iparams.num_selection_padding_residues)
        else:
            raise TypeError('Unrecognised expansion mode!')
        shell = selections.get_shell_of_residues(sel.unique_residues,
            iparams.soft_shell_cutoff_distance).atoms
        from chimerax.atomic import concatenate
        merged_sel = concatenate((sel, shell), remove_duplicates=True)
        return merged_sel


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
        print(f'Sim termination reason: {reason}')
        if reason == SimHandler.REASON_MODEL_DELETED:
            return
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
        if reason == SimHandler.REASON_COORD_LENGTH_MISMATCH:
            rt = None
        else:
            rt = self._revert_to
        if rt == 'checkpoint':
            self._current_checkpoint.revert()
        elif rt == 'start':
            print('reverting to start')
            self._starting_checkpoint.revert()
        if reason == SimHandler.REASON_COORD_LENGTH_MISMATCH:
            msg = ('Mismatch between number of simulated atoms and the model. '
                'The most common cause of this is the addition or removal of '
                'atoms while a simulation is running (this should be avoided). '
                'The simulation has been terminated to avoid data corruption, '
                'and the display reverted to a default state.')
            from ..dialog import generic_warning
            generic_warning(msg)

    def _rama_a_sim_end_cb(self, *_):
        from chimerax.core.triggerset import DEREGISTER
        if self.rama_annotator.deleted:
            return DEREGISTER
        self.rama_annotator.track_whole_model = True
        return DEREGISTER

    def _rota_a_sim_end_cb(self, *_):
        from chimerax.core.triggerset import DEREGISTER
        if self.rota_annotator.deleted:
            return DEREGISTER
        self.rota_annotator.track_whole_model = True
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
        from chimerax.core.triggerset import DEREGISTER
        if self.distance_restraint_mgr.deleted:
            return DEREGISTER
        restraints = self.distance_restraint_mgr.intra_restraints(self.sim_construct.all_atoms)
        restraints.clear_sim_indices()
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
        from chimerax.core.triggerset import DEREGISTER
        for adrm in self.adaptive_distance_restraint_mgrs:
            if adrm.deleted:
                return DEREGISTER
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
        from chimerax.core.triggerset import DEREGISTER
        sc = self.sim_construct
        if not self.proper_dihedral_restraint_mgr.deleted:
            pdrs = self.proper_dihedral_restraint_mgr.get_all_restraints_for_residues(sc.mobile_residues)
            pdrs.clear_sim_indices()
        if not self.chiral_restraint_mgr.deleted:
            crs = self.chiral_restraint_mgr.get_restraints_by_atoms(sc.mobile_atoms)
            crs.clear_sim_indices()
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
        if not self.adaptive_dihedral_restraint_mgr.deleted:
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
            # Pass the unfiltered set: update_mdff_atoms filters sim_index != -1
            # internally for the real-atom self-terms, but must also update the
            # ghost terms of symmetry parents whose OWN self-term was dropped
            # (sim_index == -1) yet whose copies still carry the map coupling.
            self.sim_handler.update_mdff_atoms(all_changeds, mgr.volume)

    def _mdff_global_k_change_cb(self, trigger_name, data):
        mgr, k = data
        self.sim_handler.set_mdff_global_k(mgr.volume, k)

    def _mdff_sim_end_cb(self, *_):
        for v, mgr in self.mdff_mgrs.items():
            if not v.deleted and not mgr.deleted:
                mdff_atoms = mgr.get_mdff_atoms(self.sim_construct.all_atoms)
                mdff_atoms.clear_sim_indices()
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER




def _symmat_to_transform12(symmat):
    '''
    Flatten a Clipper orthogonal-space operator ``symmat`` (a ``(3,4)`` ``[R | t]``
    array with R orthonormal and t in **Angstroms**) into the 12-element per-bond
    transform expected by :class:`SymmetryAwareCubicInterpMapForce` -- row-major R
    (9) followed by t converted to **nanometres** (3). This matches the
    ``SymmetrySite`` construction (which also uses ``t_nm = t * 0.1``) so a term's
    sampled position ``R.x + t_nm`` coincides with the ghost position.
    '''
    m = numpy.asarray(symmat, dtype=numpy.float64)
    return numpy.concatenate([m[:, :3].ravel(), m[:, 3] * 0.1])


class SimHandler:
    '''
    Responsible for creating a :py:class:`openmm.Simulation`, instantiating and
    populating the custom force objects, managing the creation and calling of
    :py:class:`OpenmmThreadHandler`, and generally handling all the OpenMM
    side of simulation management.
    '''
    REASON_COORD_LENGTH_MISMATCH = 'coord length mismatch'
    REASON_MODEL_DELETED = 'model deleted'
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
                - a :py:class:`SimConstruct` instance
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
        # Monotonic count of integrator steps run since the sim started, plus the
        # value captured at the last checkpoint. A cheap, authoritative
        # "is it stepping / how far" signal (see the sim_steps property).
        self._sim_steps = 0
        self._sim_steps_at_checkpoint = 0

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

        # Symmetry-aware MDFF: per-volume ``{real-atom particle index:
        # [(force_index, transform12), ...]}`` recording every map term of each
        # real atom -- its identity term (own position) and one transformed term
        # per crystallographic image -- so live coupling/enabled edits re-apply to
        # all of them. Empty when the simulation has no symmetry atoms.
        self._mdff_symmetry_term_indices = {}
        # Covered multiplicity per real-atom particle index, set alongside the
        # term map when symmetry MDFF terms are built (see
        # _add_mdff_symmetry_terms). None until then.
        self._mdff_n_sym_covered = None

        logger = self.session.logger
        # Overall simulation topology + system. The build can recover once from
        # stale isolde_template_name overrides that no longer match their
        # residue (see _build_system_with_template_recovery).
        top, system = self._build_system_with_template_recovery(
            ff, atoms, sim_construct.all_residues, sim_params, ligand_db, logger)
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

        self._system = system



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
        t = self._main_integrator.getTemperature()
        return t.value_in_unit(defaults.OPENMM_TEMPERATURE_UNIT)

    @temperature.setter
    def temperature(self, temperature):
        self._main_integrator.setTemperature(temperature)

    @property
    def smoothing(self):
        '''
        If True, the displayed coordinates will be a smoothed average of the
        last set of equilibration steps. Note that for large values of
        sim_steps_per_gui_update this can lead to distorted geometry.
        '''
        if self.sim_running:
            return self._smoother.enabled
        return self._params.trajectory_smoothing

    @smoothing.setter
    def smoothing(self, flag):
        if self.sim_running:
            self._smoother.enabled = flag
            self._thread_handler.smoothing = flag
        
    @property
    def smoothing_alpha(self):
        '''
        A value between 0 and 1 setting the strength of trajectory smoothing. 
        Smoothed coordinates are updated as (existing*alpha + new(1-alpha)).
        '''
        if self.sim_running:
            return self._smoother.smoothing_alpha
        return self._params.smoothing_alpha

    @smoothing_alpha.setter
    def smoothing_alpha(self, alpha):
        if self.sim_running:
            self._smoother.smoothing_alpha = alpha

    def update_nonbonded_softcore_parameter(self, name, value):
        param_name_map = {
            'nonbonded_softcore_lambda_minimize': 'softcore_lambda',
            'nonbonded_softcore_lambda_equil': 'softcore_lambda',
            'nonbonded_softcore_alpha': 'softcore_alpha',
        }
        if name == 'nonbonded_softcore_lambda_minimize' and not self.minimize:
            return
        if name == 'nonbonded_softcore_lambda_equil' and self.minimize:
            return
        from .custom_forces import NonbondedSoftcoreForce
        for f in self.all_forces:
            if isinstance(f, NonbondedSoftcoreForce):
                break
        else:
            return
        c = self._context
        if c is None:
            return
        c.setParameter(param_name_map[name], value)     

    @property
    def softcore_lambda(self):
        if self._context is None:
            return None
        return dict(self._context.getParameters()).get('softcore_lambda', None)
    
    @softcore_lambda.setter
    def softcore_lambda(self, val):
        if self._context is not None:
            c = self._context
            if 'softcore_lambda' in c.getParameters().keys():
                c.setParameter('softcore_lambda', val)
        elif self._system is not None:
            from .custom_forces import NonbondedSoftcoreForce
            for f in self.all_forces:
                if isinstance(f, NonbondedSoftcoreForce):
                    break
            else:
                return
            f.set_lambda(val, self._context)

    # ------------------------------------------------------------------
    # Per-group soft-core nonbonded coupling (transient; see the
    # softcore-nb-groups design). Consumed by fitting/docking engines, NOT
    # exposed to users. Usage: call enable_nb_groups() BEFORE the simulation is
    # built (before `isolde sim start`); then assign_nb_group() / set_nb_coupling()
    # freely while it runs. Edits go through the standard force_update_needed()
    # path (applied immediately when paused, on the next frame when running) with
    # no context reinitialisation. Every atom starts in group 0 with an all-ones
    # coupling table, so an enabled-but-unused simulation is physically identical
    # to a plain one.
    # ------------------------------------------------------------------
    def enable_nb_groups(self, max_groups=8):
        '''
        Provision per-group soft-core nonbonded coupling for this simulation,
        sizing the coupling table to ``max_groups`` (default 8). Must be called
        before the soft-core forces are built (before the sim is started).
        No-op (with a warning) if the soft-core potential is disabled or the
        forces are already built.
        '''
        if not self._params.use_softcore_nonbonded_potential:
            self.session.logger.warning('ISOLDE: per-group soft-core coupling '
                'requires the soft-core nonbonded potential; ignoring '
                'enable_nb_groups().')
            return
        if getattr(self, '_nb_softcore_force', None) is not None:
            self.session.logger.warning('ISOLDE: enable_nb_groups() must be called '
                'before the simulation is built; ignoring.')
            return
        self._nb_groups_max = int(max_groups)

    @property
    def nb_groups_enabled(self):
        '''Whether per-group soft-core coupling is active for this simulation.'''
        return getattr(self, '_nb_softcore_force', None) is not None

    def _nb_group_forces(self):
        forces = []
        f = getattr(self, '_nb_softcore_force', None)
        if f is not None:
            forces.append(f)
        g = getattr(self, '_nb_gbsa_force', None)
        if g is not None:
            forces.append(g)
        return forces

    def assign_nb_group(self, atoms, group_id):
        '''
        Put ``atoms`` (a ChimeraX Atoms) into nonbonded group ``group_id``
        (0 <= group_id < max_groups). Applied live to every group-aware force via
        per-particle parameter updates (no reinitialisation). Atoms outside the
        simulation construct are ignored.
        '''
        forces = self._nb_group_forces()
        if not forces:
            raise RuntimeError('per-group soft-core coupling is not enabled for this '
                'simulation (call enable_nb_groups() before starting it)')
        indices = self._atoms.indices(atoms)
        indices = [int(i) for i in indices if i != -1]
        gid = float(group_id)
        for f in forces:
            gi = f._nb_group_index
            for i in indices:
                params = list(f.getParticleParameters(i))
                params[gi] = gid
                f.setParticleParameters(i, params)
            f.update_needed = True
        pg = getattr(self, '_nb_particle_groups', None)
        if pg is not None:
            for i in indices:
                pg[i] = group_id
        self.force_update_needed()

    def set_nb_coupling(self, group_a, group_b, lam):
        '''
        Set the (symmetric) soft-core coupling between two nonbonded groups
        (0 < lam <= 1) on every group-aware force, live.
        '''
        forces = self._nb_group_forces()
        if not forces:
            raise RuntimeError('per-group soft-core coupling is not enabled')
        for f in forces:
            f.set_coupling(group_a, group_b, lam)      # sets update_needed
        self.force_update_needed()

    def get_nb_coupling(self, group_a, group_b):
        '''Current coupling between two nonbonded groups (1.0 if not enabled).'''
        f = getattr(self, '_nb_softcore_force', None)
        if f is None:
            return 1.0
        return f.get_coupling(group_a, group_b)


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

    def _build_system_with_template_recovery(self, ff, atoms, all_residues,
            sim_params, ligand_db, logger):
        '''
        Build the OpenMM topology and system from ``all_residues``.

        If the build fails because one or more residues carry an
        ``isolde_template_name`` override that does not actually match the
        residue (or names a template missing from the forcefield), clear those
        overrides and rebuild **once**. Rebuilding from scratch — rather than
        just re-running template matching — is important: ``find_residue_templates``
        deliberately skips its cysteine / USER_ / ligand logic for residues that
        already carry an override, so the residue only gets a fair shot at
        automatic assignment after the override is gone. If the residue still
        cannot be matched (or the failure was never override-related), the
        UnparameterisedResiduesError propagates and the usual widget path fires.

        Returns ``(topology, system)``.
        '''
        for attempt in (0, 1):
            template_dict = find_residue_templates(all_residues, ff,
                ligand_db=ligand_db, logger=logger,
                clear_failed_overrides=True)
            top, residue_templates = create_openmm_topology(atoms, template_dict)
            try:
                system = self._create_openmm_system(ff, top, sim_params,
                    residue_templates, residues=all_residues)
            except UnparameterisedResiduesError as e:
                offenders = [r for r in (e.unmatched + e.ambiguous)
                    if getattr(r, 'isolde_template_name', None) is not None]
                if attempt == 0 and offenders:
                    for r in offenders:
                        logger.warning(
                            'Residue {} was assigned template "{}" via its '
                            'isolde_template_name override, but that template '
                            'does not match the residue. Clearing the override '
                            'and retrying with automatic template assignment.'
                            .format(r, r.isolde_template_name))
                        r.isolde_template_name = None
                    continue
                raise
            return top, system

    def _create_openmm_system(self, forcefield, top, params, residue_templates,
            residues=None):
        residue_to_template, ambiguous, unassigned = forcefield.assignTemplates(
            top, ignoreExternalBonds=True, explicit_templates=residue_templates
        )
        if len(ambiguous) or len(unassigned):
            # Map the offending OpenMM topology residues back to their ChimeraX
            # residues (OpenMM res.index is the order they were added in
            # create_openmm_topology, which matches `residues`) so the caller
            # can inspect/clear stale isolde_template_name overrides.
            unmatched_cx = []
            ambiguous_cx = []
            if residues is not None:
                unmatched_cx = [residues[r.index] for r in unassigned]
                ambiguous_cx = [residues[r.index] for r in ambiguous]
            raise UnparameterisedResiduesError(unmatched=unmatched_cx,
                ambiguous=ambiguous_cx)

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
        self._set_core_force_groups(sys)
        return sys

    def _set_core_force_groups(self, system):
        from openmm.openmm import HarmonicBondForce, HarmonicAngleForce, PeriodicTorsionForce
        for i, f in enumerate(system.getForces()):
            if isinstance(f, HarmonicBondForce):
                f.setForceGroup(BOND_FORCE_GROUP)
            elif isinstance(f, HarmonicAngleForce):
                f.setForceGroup(ANGLE_FORCE_GROUP)
            elif isinstance(f, PeriodicTorsionForce):
                f.setForceGroup(DIHEDRAL_FORCE_GROUP)

    def _convert_to_soft_core_potentials(self, system):
        from openmm.openmm import NonbondedForce
        for i,f in enumerate(system.getForces()):
            if type(f) == NonbondedForce:
                break
        from .custom_forces import (NonbondedSoftcoreForce,
            NonbondedSoftcoreExceptionForce, NBGroupNonbondedSoftcoreForce,
            SymmetryAwareNonbondedSoftcoreForce)
        p = self._params
        param_dict = {
            'a': p.nonbonded_softcore_a,
            'b': p.nonbonded_softcore_b,
            'c': p.nonbonded_softcore_c,
            'nb_lambda': p.nonbonded_softcore_lambda_minimize,
            'alpha': p.nonbonded_softcore_alpha,
        }
        # Two orthogonal per-particle grouping layers may be active:
        #  * crystallographic symmetry (symgroup + grouptable mask), if copies are
        #    present -- SymmetryAwareMixin; and
        #  * per-group soft-core coupling (nb_group + nb_coupling_table lambda),
        #    if an engine called enable_nb_groups() before the sim was built.
        # The nb-group force is the base for the symmetry-aware one, so any
        # combination composes. Each active layer appends a trailing per-particle
        # group id, base (nb_group) before mixin (symgroup).
        groups = getattr(self, '_symmetry_particle_groups', None)
        nb_max = int(getattr(self, '_nb_groups_max', 1) or 1)
        nb_on = nb_max > 1
        n_particles = system.getNumParticles()
        if nb_on:
            # Every atom starts in group 0 (inert: coupling table is all-ones);
            # the engine assigns groups and couplings live afterwards.
            nb_groups = self._nb_particle_groups = numpy.zeros(n_particles, dtype=int)
        else:
            nb_groups = None
            self._nb_particle_groups = None

        if groups is not None:
            sf = SymmetryAwareNonbondedSoftcoreForce(
                symmetry_ngroups=self._symmetry_ngroups,
                symmetry_group_weights=getattr(self, '_symmetry_group_table', None),
                n_nb_groups=nb_max,
                **param_dict)
        elif nb_on:
            sf = NBGroupNonbondedSoftcoreForce(n_nb_groups=nb_max, **param_dict)
        else:
            sf = NonbondedSoftcoreForce(**param_dict)
        sfb = NonbondedSoftcoreExceptionForce(**param_dict)
        sf.setForceGroup(NONBONDED_FORCE_GROUP)
        sfb.setForceGroup(NONBONDED_FORCE_GROUP)
        sf.setNonbondedMethod(f.getNonbondedMethod())
        sf.setCutoffDistance(f.getCutoffDistance())
        sf.setSwitchingDistance(f.getSwitchingDistance())
        for j in range(n_particles):
            charge, sigma, epsilon = f.getParticleParameters(j)
            pp = [charge, sigma, epsilon]
            if nb_on:
                pp.append(float(nb_groups[j]))
            if groups is not None:
                pp.append(float(groups[j]))
            sf.addParticle(pp)
        for j in range(f.getNumExceptions()):
            p1, p2, cp, sig, eps = f.getExceptionParameters(j)
            sf.addExclusion(p1, p2)
            sfb.addBond(p1, p2, (cp, sig, eps))
        system.removeForce(i)
        system.addForce(sf)
        system.addForce(sfb)
        self.all_forces.extend((sf, sfb))
        if nb_on:
            # Handle for the group manager (nb_group membership + coupling table).
            self._nb_softcore_force = sf
        self._softcore_nb_param_mgr = _SoftCoreNonbondedParamMgr(self, self._params)


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
                  :py:class:`AdaptiveTorsionForce`)
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
        # Crystallographic symmetry copies must be added to the System (and the
        # base NonbondedForce) *before* GBSA and the soft-core conversion, both
        # of which size themselves to the current particle count.
        shell = getattr(self._sim_construct, 'symmetry_shell', None)
        if shell is not None and not hasattr(self, '_symmetry_initialized'):
            self.initialize_symmetry_copies(shell)
            self._symmetry_initialized = True
        if params.use_gbsa and not hasattr(self, '_gbsa_force'):
            self.initialize_implicit_solvent(params)
        if params.use_softcore_nonbonded_potential and not hasattr(self, '_soft_core_initialized'):
            self._convert_to_soft_core_potentials(self._system)
            self._soft_core_initialized=True
        integrator = self._prepare_integrator(params)
        platform = openmm.Platform.getPlatformByName(params.platform)

        properties = {}
        device_index = params.device_index
        if device_index is not None:
            if params.platform=="CUDA":
                properties['CudaDeviceIndex']=str(device_index)
            elif params.platform=='OpenCL':
                properties['OpenCLDeviceIndex']=str(device_index)


        from openmm import app
        logger.status('Initialising primary simulation object')
        s = self._simulation = app.Simulation(self.topology, self._system,
            integrator, platform, properties)
        # Workaround for https://github.com/openmm/openmm/issues/4038 (affecting OpenMM 8.0.0)
        # Can be safely removed once OpenMM is next updated.
        self._smoother.setGlobalVariableByName('reset_smooth', 1.0)
        # End workaround
        c = self._context = s.context
        n_real = len(self._atoms)
        n_total = self._system.getNumParticles()
        if n_total > n_real:
            # Symmetry copies occupy particles [n_real, n_total). Seed them at
            # their parent's position (any valid placeholder), then let OpenMM
            # snap each virtual site to R.r + t before the first step.
            pos = numpy.zeros((n_total, 3))
            pos[:n_real] = 0.1*self._atoms.coords
            pos[n_real:] = pos[self._symmetry_copies['parent_index']]
            c.setPositions(pos)
            c.computeVirtualSites()
        else:
            c.setPositions(0.1*self._atoms.coords)
        c.setVelocitiesToTemperature(self.temperature)
        self._thread_handler = OpenmmThreadHandler(c, params,
            num_real_atoms=n_real)
        self.smoothing = params.trajectory_smoothing
        self.smoothing_alpha = params.smoothing_alpha
        logger.status('')

    def _prepare_integrator(self, params):
        from openmm.openmm import CompoundIntegrator
        from .custom_integrators import VelocityChecker, Smoother
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
        _integrator = self._main_integrator = integrator(*integrator_params)
        _integrator.setConstraintTolerance(params.constraint_tolerance)
        vcheck = VelocityChecker()
        smoother = self._smoother = Smoother(params.smoothing_alpha)
        main_integrator = CompoundIntegrator()
        main_integrator.addIntegrator(_integrator)
        main_integrator.addIntegrator(vcheck)
        main_integrator.addIntegrator(smoother)

        return main_integrator

    def start_sim(self):
        '''
        Start the main simulation loop. Automatically runs a minimisation, then
        switches to equilibration once minimisation is converged.
        '''
        if self._sim_running:
            raise RuntimeError('Simulation is already running!')
        from openmm import OpenMMException
        try:
            self._prepare_sim()
        except OpenMMException as e:
            if self._params.platform=='CUDA':
                self.session.logger.warning(f'Launching using CUDA failed with the below message. Falling back to using OpenCL.\n\n{str(e)}')
                self._params.platform='OpenCL'
                self._prepare_sim()
            else:
                raise

        
        self._pause = False
        self._stop = False
        self._sim_running = True
        self._startup = True
        self._startup_counter = 0
        self._sim_steps = 0
        self._sim_steps_at_checkpoint = 0
        from chimerax.core.models import REMOVE_MODELS
        self._model_deleted_handler = self.session.triggers.add_handler(REMOVE_MODELS, self._model_deleted_cb)
        self._minimize_and_go()

    def find_clashing_atoms(self, max_force = defaults.CLASH_FORCE):
        if not self._sim_running:
            raise RuntimeError('Simulation must be running first!')
        c = self._context
        state = c.getState(getPositions=True, getForces=True)
        forces = state.getForces(asNumpy = True)
        import numpy
        force_mags = numpy.linalg.norm(forces, axis=1)
        sort_i = numpy.argsort(force_mags)[::-1]
        sorted_forces = numpy.sort(force_mags)[::-1]
        fmask = (sorted_forces > max_force)
        sc = self._sim_construct
        mmask = numpy.zeros(fmask.shape, bool)
        mmask[sc.all_atoms.indices(sc.mobile_atoms)] = True

        mask = numpy.logical_and(fmask, mmask[sort_i])
        clashes = sort_i[mask]
        return clashes, sorted_forces[mask]

    def get_excessive_forces(self, force_groups=(BOND_FORCE_GROUP, ANGLE_FORCE_GROUP, DIHEDRAL_FORCE_GROUP, NONBONDED_FORCE_GROUP), threshold=defaults.STRESS_FORCE_THRESHOLD):
        if not self._sim_running:
            raise RuntimeError('Simulation must be running first!')
        force_map = {}
        c = self._context
        for fg in force_groups:
            if not (0 <= fg <= 31):
                raise ValueError('Force groups must be between 0 and 31 inclusive!')
            state = c.getState(getForces=True, groups=set([fg]))
            # Restrict to real atoms: in a symmetry-aware sim the context also holds
            # symmetry-copy virtual-site particles ([n_real, n_total)), whose indices
            # would be meaningless to callers mapping back through self._atoms.
            forces = state.getForces(asNumpy=True)[:len(self._atoms)]
            magnitudes = numpy.linalg.norm(forces, axis=1)
            excessive = numpy.argwhere(magnitudes > threshold)
            force_map[fg] = (excessive, magnitudes[excessive])
        return force_map

    def _minimize_and_go(self):
        th = self.thread_handler
        self._delayed_reaction_handler = delayed_reaction(self.session.triggers, 'new frame', th.minimize, [],
            th.thread_finished, self._update_coordinates_and_repeat, [True])

    def _repeat_step(self):
        th = self.thread_handler
        params = self._params
        if self._context_reinit_pending:
            f = self._reinitialize_context
            f_args = []
            final_args = []
        elif th.unstable() or self._unstable or not th.minimization_converged:
            self.softcore_lambda = params.nonbonded_softcore_lambda_minimize
            f = th.minimize
            f_args = []
            final_args = [True]
        elif self.minimize:
            self.softcore_lambda = params.nonbonded_softcore_lambda_minimize
            f = th.minimize
            f_args=[params.minimization_convergence_tol_end]
            final_args = [True]
        else:
            if th._last_mode == 'min':
                self.softcore_lambda = params.nonbonded_softcore_lambda_equil            
            f = th.step
            f_args = (params.sim_steps_per_gui_update,)
            final_args = []
        self._delayed_reaction_handler = delayed_reaction(self.session.triggers, 'new frame', f, f_args,
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
        if self._sim_construct.model.deleted:
            self.stop(reason=self.REASON_MODEL_DELETED)
        if self._stop:
            self._delayed_reaction_handler.cancel()
            # self._thread_handler.delete()
            self._thread_handler = None
            self._simulation = None
            self._sim_running = False
            self._model_deleted_handler.remove()
            self.triggers.activate_trigger('sim terminated', self._stop_reason)
            return
        if self._startup:
            self._startup_counter += 1
            if self._startup_counter >= self._params.simulation_startup_rounds:
                self._startup = False
        th = self.thread_handler
        try:
            if not (self._stop and self._stop_reason==self.REASON_MODEL_DELETED):
                self.atoms.coords = th.coords
        except ValueError:
            self.stop(reason=self.REASON_COORD_LENGTH_MISMATCH)
        # Count integrator steps completed this cycle (equilibration only;
        # minimisation rounds are not fixed-step).
        if th._last_mode == 'equil':
            self._sim_steps += self._params.sim_steps_per_gui_update
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
        if self._force_update_pending and not self._stop:
            self._update_forces_in_context_if_needed()
        if reinit_vels and not self._stop:
            th.reinitialize_velocities()
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
        Returns the :py:class:`OpenmmThreadHandler` object.
        '''
        return self._thread_handler

    @property
    def sim_steps(self):
        '''
        Total number of integrator (equilibration) steps run since the
        simulation started. Monotonic; a cheap, authoritative signal that the
        simulation is actually advancing (and by how much).
        '''
        return self._sim_steps

    @property
    def steps_since_checkpoint(self):
        '''Integrator steps run since the last saved checkpoint.'''
        return self._sim_steps - self._sim_steps_at_checkpoint

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
        :cpp:class:`OpenmmThreadHandler` object, rendering the Python class
        unusable.
        '''
        self._stop = True
        self._stop_reason = reason
        if self.pause:
            if hasattr(self, '_delayed_reaction_handler'):
                self._delayed_reaction_handler.cancel()
            # self._thread_handler.delete()
            self._thread_handler = None
            self._simulation = None
            self._sim_running = False
            self._model_deleted_handler.remove()
            self.triggers.activate_trigger('sim terminated', reason)

    @property
    def sim_running(self):
        ''' Is the simulation curently running (i.e. started and not destroyed)? '''
        return self._sim_running

    def force_update_needed(self):
        '''
        This must be called after any changes to force objects to ensure
        the changes are pushed to the simulation context. This happens
        automatically when changes are made through the :py:class:`SimManager`
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

    def _model_deleted_cb(self, _, models):
        if self._sim_construct.model in models:
            self.stop(reason=self.REASON_MODEL_DELETED)
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
        df.setForceGroup(RESTRAINT_FORCE_GROUP)
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
        from math import cos
        force = self._dihedral_restraint_force
        all_atoms = self._atoms
        dihedral_atoms = restraint.dihedral.atoms
        indices = [all_atoms.index(atom) for atom in dihedral_atoms]
        restraint.sim_index = force.addTorsion(*indices,
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
        from .custom_forces import AdaptiveTorsionForce
        df = self._adaptive_dihedral_restraint_force = AdaptiveTorsionForce()
        df.setForceGroup(RESTRAINT_FORCE_GROUP)
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
            restraints.enableds, restraints.spring_constants,
            restraints.targets, restraints.kappas, restraints.alphas)
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
        restraint.sim_index = force.addTorsion(*indices,
            float(restraint.enabled), restraint.spring_constant,
            restraint.target, restraint.kappa, restraint.alpha)
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
            restraints.enableds, restraints.spring_constants,
            restraints.targets, restraints.kappas, restraints.alphas)
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
            target=restraint.target, kappa=restraint.kappa,
            alpha=restraint.alpha)
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
        tf.setForceGroup(RESTRAINT_FORCE_GROUP)
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
        f.setForceGroup(RESTRAINT_FORCE_GROUP)
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
        rf.setForceGroup(RESTRAINT_FORCE_GROUP)
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
        f.setForceGroup(RESTRAINT_FORCE_GROUP)
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
        from .custom_forces import CubicInterpMapForce, SymmetryAwareCubicInterpMapForce
        v = volume
        region = list(v.region)
        # Ensure that the region ijk step size is [1,1,1]
        region[-1] = [1,1,1]
        data = v.region_matrix(region=region)
        if any (dim < 3 for dim in data.shape):
            # Not enough of this map is covering the atoms. Leave it out.
            return
        if not data.data.c_contiguous:
            data_copy = numpy.empty(data.shape, numpy.float32)
            data_copy[:] = data
            data = data_copy
        from chimerax.geometry import Place
        tf = v.data.xyz_to_ijk_transform
        # Shift the transform to the origin of the region
        region_tf = Place(axes=tf.axes(), origin = tf.origin() -
            v.data.xyz_to_ijk(v.region_origin_and_step(region)[0]))
        # In OpenMM forces, parameters can only be per-particle, or global to
        # the entire context. So if we want a parameter that's constant to all
        # particles in this force, it needs a unique name so it doesn't
        # interfere with other instances of the same force.
        suffix = str(len(self.mdff_forces)+1)
        # When simulating with crystallographic symmetry, use the symmetry-aware
        # map force: each term carries a per-bond operator so a real atom can feel
        # the map through its symmetry image(s) (sampling at S.r), with the force
        # folded back to the real atom by the expression's own chain rule. Plain
        # sims keep the base force unchanged.
        shell = getattr(self._sim_construct, 'symmetry_shell', None)
        force_cls = (SymmetryAwareCubicInterpMapForce if shell is not None
                     else CubicInterpMapForce)
        f = force_cls(data, region_tf.matrix, suffix, units='angstroms',
                                map_sigma=v.sigma)
        f.setForceGroup(MAP_FORCE_GROUP)
        self.all_forces.append(f)
        self._system.addForce(f)
        self.mdff_forces[v] = f
        from chimerax.clipper.maps import XmapHandler_Live
        if isinstance(v, XmapHandler_Live):
            class RfactorChecker:
                def __init__(self, v):
                    mgr = self.mgr = v.xmap_mgr
                    self.best_r = mgr.rwork
                
                def r_improved(self):
                    r = self.mgr.rwork
                    if r < self.best_r:
                        self.best_r = r
                        return True
                    return False
            def data_changed_cb(reason, v=volume, r=region, checker=RfactorChecker(volume)):
                if reason=='values changed' and checker.r_improved():
                    f.update_map_data(v.region_matrix(region=r))
                    # Keep the sigma normalisation in lockstep with the data the
                    # force actually samples (only swapped in when R-work improves).
                    f.set_map_sigma(v.sigma)
                    self.force_update_needed()
        else:                
            def data_changed_cb(reason, v = volume, r = region):
                if reason == 'values changed':
                    f.update_map_data(v.region_matrix(region=r))
                    f.set_map_sigma(v.sigma)
                    self.force_update_needed()
        v.data.add_change_callback(data_changed_cb)
        def remove_change_cb(*_):
            if v.deleted:
                return
            v.data.remove_change_callback(data_changed_cb)
        self.triggers.add_handler('sim terminated', remove_change_cb)
        

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
        ks = numpy.asarray(mdff_atoms.coupling_constants, dtype=numpy.float64)
        enableds = numpy.asarray(mdff_atoms.enableds, dtype=numpy.float64)
        shell = getattr(self._sim_construct, 'symmetry_shell', None)
        if shell is None:
            mdff_atoms.sim_indices = f.add_atoms(indices, ks, enableds)
        else:
            self._add_mdff_symmetry_terms(f, volume, shell, mdff_atoms,
                indices, ks, enableds)
        self.context_reinit_needed()

    def _add_mdff_symmetry_terms(self, f, volume, shell, mdff_atoms,
            indices, ks, enableds):
        '''
        Symmetry-aware MDFF term creation (:class:`SymmetryAwareCubicInterpMapForce`).

        Every term lives on a **real** atom and carries a per-bond operator; a
        real atom feels the map through both its own position (identity term, only
        where that position is inside the covered box) and each of its
        crystallographic images (transformed terms sampling the ghost position
        ``S.r``). All representations of an atom fold the identical force back to
        it, so the coupling is split by the covered multiplicity
        ``n_sym_covered = (self in ROI) + (# copies)`` and sums to ``k``. No terms
        are placed on copy (virtual-site) particles. Records, per real-atom
        particle index, the ``(force_index, transform12)`` of each of its terms so
        live coupling/enabled edits can re-apply them.
        '''
        from chimerax.atomic import Atoms
        from .custom_forces import SymmetryAwareCubicInterpMapForce as _SA
        all_atoms = self._atoms
        n_real = len(all_atoms)
        # Which real atoms' OWN position is inside the covered map box (the region
        # of interest; full mobile set for whole-model / plain sims).
        cov = getattr(self._sim_construct, 'symmetry_coverage_atoms', None)
        if cov is None:
            cov = self._sim_construct.mobile_atoms
        self_covered = numpy.zeros(n_real, dtype=bool)
        ci = all_atoms.indices(cov)
        self_covered[ci[ci != -1]] = True
        # Per-copy parent particle indices + covered multiplicity per real atom.
        copies = shell.copies
        parent_pidx = all_atoms.indices(Atoms([c.parent_atom for c in copies]))
        counts = numpy.bincount(parent_pidx[parent_pidx != -1], minlength=n_real)
        n_sym_cov = self_covered.astype(numpy.int64) + counts[:n_real]
        # Unscaled (k, enabled) per MDFF real atom, for both term kinds.
        kmap = {int(i): (float(k), float(e))
            for i, k, e in zip(indices, ks, enableds)}
        ident = _SA.IDENTITY_TRANSFORM
        add_i, add_k, add_e, add_tf, owner, is_ident = [], [], [], [], [], []
        # Identity terms: covered MDFF atoms sample the map at their own position.
        for i, k, e in zip(indices, ks, enableds):
            i = int(i)
            # i < 0 (atom not in the construct) should not occur -- mdff_atoms are
            # sourced from the mobile set -- but guard it so a stray -1 can't
            # negative-index into self_covered / n_sym_cov and add a bogus term.
            if i < 0 or not self_covered[i] or n_sym_cov[i] == 0:
                continue
            add_i.append(i); add_k.append(float(k) / n_sym_cov[i])
            add_e.append(float(e)); add_tf.append(ident)
            owner.append(i); is_ident.append(True)
        # Transformed terms: one per shell copy whose parent is an MDFF atom,
        # sampling the ghost position under that copy's operator.
        for c, pidx in zip(copies, parent_pidx):
            pidx = int(pidx)
            if pidx < 0:
                continue
            ke = kmap.get(pidx)
            if ke is None or n_sym_cov[pidx] == 0:
                continue    # parent has no map coupling (e.g. H) -> no ghost term
            k, e = ke
            add_i.append(pidx); add_k.append(k / n_sym_cov[pidx])
            add_e.append(e); add_tf.append(_symmat_to_transform12(
                shell.symmats[c.symop_index]))
            owner.append(pidx); is_ident.append(False)
        sim_indices = numpy.full(len(indices), -1, dtype=numpy.int32)
        d = self._mdff_symmetry_term_indices.setdefault(volume, {})
        if add_i:
            fis = f.add_atoms(
                numpy.array(add_i, dtype=numpy.int32),
                numpy.array(add_k, dtype=numpy.float64),
                numpy.array(add_e, dtype=numpy.float64),
                transforms=numpy.array(add_tf, dtype=numpy.float64))
            pos_in_indices = {int(a): p for p, a in enumerate(indices)}
            for fi, own, tf, isid in zip(fis, owner, add_tf, is_ident):
                d.setdefault(own, []).append((int(fi), tf))
                if isid:
                    sim_indices[pos_in_indices[own]] = int(fi)
        mdff_atoms.sim_indices = sim_indices
        # Covered multiplicity, cached for the live-update path.
        self._mdff_n_sym_covered = n_sym_cov

    def add_mdff_atom(self, mdff_atom, volume):
        '''
        Add a single MDFF atom proxy to the force associated with the given
        volume. Delegates to :func:`add_mdff_atoms` (wrapping the proxy in a
        one-element collection over the same underlying object) so the
        symmetry-aware term construction and the ``sim_index`` write-back are
        handled by one code path. Automatically calls
        :func:`context_reinit_needed`.

        Args:
            * mdff_atom:
                - a :py:class:`MDFFAtom` instance
            * volume:
                - the :py:class:`chimerax.Volume` instance that was used to
                  create the target force.
        '''
        from ..molarray import MDFFAtoms
        self.add_mdff_atoms(MDFFAtoms([mdff_atom]), volume)

        ##
        # During simulation
        ##

    def set_mdff_global_k(self, volume, k):
        '''
        Set the global coupling constant for the MDFF force associated with
        the given volume. 

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
        from chimerax.geometry import Place
        region_tf = Place(axes=tf.axes(), origin = tf.origin() -
            volume.data.xyz_to_ijk(volume.region_origin_and_step(region)[0]))

        if self.sim_running:
            context = self.thread_handler.context
            f.update_transform(region_tf.matrix, context=context)


    def set_mdff_magnification_factor(self, volume, magnification_factor):
        '''
        Adjust the dimensions (coordinate magnification) of the mdff map in the
        simulation. This is typically used to optimize the scaling in the course of
        a single simulation. The final magnification should be applied back to the
        original map, so that in future simulations the factor is 1.0.

        Args:
            * volume:
                - the :py:class:`chimerax.Volume` instance that was used to
                  create the target force
            * magnification_factor:
                - the new magnification factor (dimensionless). Changing it by
                  more than 1-2% in a single go is dangerous!
        '''
        f = self.mdff_forces[volume]
        if self.sim_running:
            context = self.thread_handler.context
            f.set_map_magnification_factor(magnification_factor, context=context)
        else:
            f.set_map_magnification_factor(magnification_factor)
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
        shell = getattr(self._sim_construct, 'symmetry_shell', None)
        if shell is None:
            m = mdff_atoms[mdff_atoms.sim_indices != -1]
            ks = numpy.asarray(m.coupling_constants, dtype=numpy.float64)
            f.update_atoms(m.sim_indices, ks, m.enableds)
        else:
            self._update_mdff_symmetry_terms(f, volume, mdff_atoms)
        self.force_update_needed()

    def _update_mdff_symmetry_terms(self, f, volume, mdff_atoms):
        '''
        Re-apply a coupling/enabled edit to ALL of each changed atom's MDFF terms
        (its identity term and each transformed symmetry term), rescaled by the
        covered multiplicity. Keyed off the real-atom particle index in
        :attr:`_mdff_symmetry_term_indices` (which stores ``(force_index, transform)``
        per term), so it correctly reaches parents whose self-term was dropped
        (they appear only via their transformed terms). The fixed per-term
        transform is re-supplied because the C++ update overwrites all per-bond
        parameters.
        '''
        d = self._mdff_symmetry_term_indices.get(volume)
        if not d:
            return
        nsym = self._mdff_n_sym_covered
        idx = self._atoms.indices(mdff_atoms.atoms)
        u_i, u_k, u_e, u_tf = [], [], [], []
        for i, k, e in zip(idx, mdff_atoms.coupling_constants,
                mdff_atoms.enableds):
            i = int(i)
            if i < 0 or nsym[i] == 0:
                continue
            scaled = float(k) / nsym[i]
            for fi, tf in d.get(i, ()):
                u_i.append(fi); u_k.append(scaled)
                u_e.append(float(e)); u_tf.append(tf)
        if u_i:
            f.update_atoms(
                numpy.array(u_i, dtype=numpy.int32),
                numpy.array(u_k, dtype=numpy.float64),
                numpy.array(u_e, dtype=numpy.float64),
                transforms=numpy.array(u_tf, dtype=numpy.float64))

    def update_mdff_atom(self, mdff_atom, volume):
        '''
        Update the simulation to reflect the new parameters (individual coupling
        constant, enabled/disabled state) for the given MDFF atom proxy. Delegates
        to :func:`update_mdff_atoms` (one-element collection) so the
        symmetry-aware term handling lives in a single path.

        Args:
            * mdff_atom:
                - a :py:class:`MDFFAtom` instance
            * volume:
                - the :py:class:`chimerax.Volume` instance that was used to
                  create the target force.
        '''
        from ..molarray import MDFFAtoms
        self.update_mdff_atoms(MDFFAtoms([mdff_atom]), volume)


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
        indices = self._atoms.indices(atoms).tolist()
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
        from openmm.app import ForceField
        ff = ForceField(*[f for f in forcefield_file_list if f is not None])
        return ff


    def initialize_symmetry_copies(self, shell):
        '''
        Add the crystallographic symmetry copies described by ``shell`` (a
        :py:class:`chimerax.isolde.openmm.symmetry_sim.SymmetryShell`) to the
        OpenMM system as zero-mass :py:class:`openmm.SymmetrySite` virtual
        sites.

        Each copy is slaved to its parent real atom via the operator's
        orthogonal ``r' = R.r + t`` transform; forces on the copy fold back onto
        the parent as ``R^T.f``, so a copy of a *mobile* atom makes the crystal
        contact fully two-way with no extra bookkeeping.

        Must be called **before** :func:`initialize_implicit_solvent` and
        :func:`_convert_to_soft_core_potentials`: every per-particle force must
        cover ``system.getNumParticles()`` particles, and both of those size
        themselves to the particle count at build time, so the copies have to
        exist first.

        Records, for the force wiring and later coordinate/MDFF handling:

            * ``self._num_real_atoms``           -- ``len(self._atoms)``.
            * ``self._symmetry_particle_groups`` -- int array of length
              ``system.getNumParticles()``: 0 for real atoms, a contiguous
              per-operator id (1..M) for copies (used to index ``grouptable``).
            * ``self._symmetry_ngroups``         -- ``M + 1``.
            * ``self._symmetry_copies``          -- dict of parallel arrays
              (parent_index, group, R, t_nm) describing each copy, for
              coordinate seeding and the per-instance MDFF ``/n_sym`` scaling.
        '''
        import openmm
        from openmm import NonbondedForce
        params = self._params
        n_real = self._num_real_atoms = len(self._atoms)
        system = self._system
        if not params.use_softcore_nonbonded_potential:
            # The validated single-force exclusion design lives in the soft-core
            # nonbonded force's group mask; without soft-core there is no clean,
            # scalable way to suppress same-operator copy-copy pairs. Symmetry-
            # aware sims therefore require the soft-core potential in v1.
            self.session.logger.warning('ISOLDE: symmetry-aware simulation '
                'requires the soft-core nonbonded potential; skipping symmetry '
                'atoms for this simulation.')
            return
        # Base NonbondedForce: copies inherit their parent's (charge, sigma,
        # epsilon). The soft-core conversion later reads these straight back off.
        for nb in system.getForces():
            if type(nb) == NonbondedForce:
                break
        else:
            raise RuntimeError('No NonbondedForce found to attach symmetry '
                'copies to!')

        # Resolve each unique parent atom to its real-particle index once.
        parents = shell.parent_atoms
        parent_indices = self._atoms.indices(parents)
        if -1 in parent_indices:
            raise RuntimeError('A symmetry-copy parent atom is missing from the '
                'simulation construct!')
        parent_index_of = {id(a): int(i)
            for a, i in zip(parents, parent_indices)}

        # Dense, contiguous operator ids (grouptable rows): real atoms are group
        # 0; each distinct operator present becomes 1..M in first-seen order.
        op_to_group = {}
        groups = [0] * n_real
        copy_parent, copy_group, copy_R, copy_t = [], [], [], []
        for copy in shell.copies:
            symop = copy.symop_index
            group = op_to_group.get(symop)
            if group is None:
                group = op_to_group[symop] = len(op_to_group) + 1
            pidx = parent_index_of[id(copy.parent_atom)]
            R, t = shell.operator(symop)
            t_nm = t * 0.1                       # Angstrom -> nm
            cidx = system.addParticle(0.0)       # zero mass -> virtual site
            charge, sigma, epsilon = nb.getParticleParameters(pidx)
            nb.addParticle(charge, sigma, epsilon)
            # SymmetrySite(parent, Rrow0, Rrow1, Rrow2, v, useBoxVectors=False):
            # Cartesian mode (no PBC), fed Clipper's orthogonal-space operator.
            system.setVirtualSite(cidx, openmm.SymmetrySite(pidx,
                openmm.Vec3(*R[0]), openmm.Vec3(*R[1]), openmm.Vec3(*R[2]),
                openmm.Vec3(*t_nm), False))
            groups.append(group)
            copy_parent.append(pidx)
            copy_group.append(group)
            copy_R.append(R)
            copy_t.append(t_nm)

        self._symmetry_ngroups = len(op_to_group) + 1
        self._symmetry_particle_groups = numpy.array(groups, dtype=int)
        # Operator-set-aware group weight table (Phase 2c): each group's
        # orthogonal [R|t] (identity for group 0), fed to symmetry_group_weight_table
        # so cross-operator copy<->copy contacts are counted correctly in LOCAL
        # sims (where the fixed copy-copy=0 rule would drop them). Shared by the
        # nonbonded and GBSA forces built later in _prepare_sim.
        ops_by_group = [numpy.concatenate([numpy.eye(3), numpy.zeros((3, 1))], axis=1)]
        ops_by_group += [None] * len(op_to_group)
        for symop, group in op_to_group.items():
            R, t = shell.operator(symop)
            ops_by_group[group] = numpy.concatenate(
                [R, numpy.asarray(t, dtype=float)[:, None]], axis=1)
        from .symmetry_sim import symmetry_group_weight_table
        self._symmetry_group_table = symmetry_group_weight_table(ops_by_group)
        parent_index = numpy.array(copy_parent, dtype=int)
        self._symmetry_copies = {
            'parent_index': parent_index,
            'group':        numpy.array(copy_group, dtype=int),
            'R':            numpy.array(copy_R, dtype=float),   # (ncopy, 3, 3)
            't_nm':         numpy.array(copy_t, dtype=float),   # (ncopy, 3)
        }
        # NB: MDFF map coupling is handled entirely on the real atoms by the
        # symmetry-aware map force (see _add_mdff_symmetry_terms, which computes
        # its own covered multiplicity from the shell); the copies here carry only
        # nonbonded + GBSA, so no per-atom map multiplicity is stored.
        self.session.logger.info(f'ISOLDE: added {len(shell.copies)} symmetry-'
            f'copy virtual site(s) spanning {len(op_to_group)} crystallographic '
            'operator(s) to the simulation.')

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
        from openmm.openmm import NonbondedForce
        for f in system.getForces():
            if isinstance(f, NonbondedForce):
                break
        n = f.getNumParticles()
        pparams = numpy.empty((n, 3))

        for i in range(f.getNumParticles()):
            pparams[i,0] = f.getParticleParameters(i)[0].value_in_unit(
                defaults.OPENMM_CHARGE_UNIT
            )
        # Crystallographic symmetry copies (if any) occupy the particles beyond
        # the real atoms. GBSA cannot exclude particles (every per-particle force
        # must cover all of them), so copies join this force too; their solvation
        # is made symmetry-correct by the masks in SymmetrySoftCoreGBSAGBnForce
        # (see that class). Their charges are already set above (read from the
        # extended NonbondedForce); their (or, sr) are copied from their parent.
        groups = getattr(self, '_symmetry_particle_groups', None)
        nb_max = int(getattr(self, '_nb_groups_max', 1) or 1)
        # Per-group soft-core GB only applies with the soft-core potential.
        nb_on = nb_max > 1 and params.use_softcore_nonbonded_potential
        from .custom_forces import (GBSAForce, SoftCoreGBSAGBnForce,
            NBGroupSoftCoreGBSAGBnForce, SymmetrySoftCoreGBSAGBnForce)
        if params.use_softcore_nonbonded_potential:
            gbsa_params['nb_lambda'] = params.nonbonded_softcore_lambda_minimize
            if groups is not None:
                # Symmetry GB (nb-group layer is inert here until Phase 5).
                gbforce = SymmetrySoftCoreGBSAGBnForce(n_nb_groups=nb_max, **gbsa_params)
            elif nb_on:
                gbforce = NBGroupSoftCoreGBSAGBnForce(n_nb_groups=nb_max, **gbsa_params)
            else:
                gbforce = SoftCoreGBSAGBnForce(**gbsa_params)
        else:
            gbforce = GBSAForce(**gbsa_params)
        self._gbsa_force = gbforce
        std = gbforce.getStandardParameters(top)   # (n_real, 2): (or, sr)
        n_real = len(std)
        pparams[:n_real, 1:] = std
        if groups is not None:
            # Each copy inherits its parent's offset/scaled radius.
            parent_index = self._symmetry_copies['parent_index']
            pparams[n_real:, 1:] = pparams[parent_index, 1:]
            gbforce._symmetry_ngroups = self._symmetry_ngroups
            gbforce._symmetry_groups = groups
            gbforce._symmetry_group_table = getattr(self, '_symmetry_group_table', None)
        elif nb_on:
            # Per-group soft-core GB (no symmetry): every atom starts in group 0;
            # the group manager updates membership + couplings live afterwards.
            gbforce._nb_groups = numpy.zeros(n, dtype=int)
            self._nb_gbsa_force = gbforce
        gbforce.addParticles(pparams)
        gbforce.finalize()
        system.addForce(gbforce)
        # set the base NonbondedForce dielectric to vacuum
        f.setReactionFieldDielectric(1.0)

class _SoftCoreNonbondedParamMgr:
    def __init__(self, sim_handler, param_mgr):
        self.sim_handler = sim_handler
        self.param_mgr = param_mgr
        self._param_changed_handler = param_mgr.triggers.add_handler(
            param_mgr.PARAMETER_CHANGED, self._param_changed_cb
        )
        sim_handler.triggers.add_handler('sim terminated', self._sim_end_cb)

    def _param_changed_cb(self, _, data):
        param, value = data
        if param in (
            'nonbonded_softcore_lambda_minimize',
            'nonbonded_softcore_lambda_equil',
            'nonbonded_softcore_alpha',
        ):
            self.sim_handler.update_nonbonded_softcore_parameter(param, value)
    
    def _sim_end_cb(self, *_):
        self.param_mgr.triggers.remove_handler(self._param_changed_handler)




class MapScaleOptimizer:
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
        from openmm import unit
        state = self.context.getState(getEnergy=True, groups=CORE_FORCE_GROUPS)
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
            # sh.set_mdff_magnification_factor(v, scale)
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
            from chimerax.save_command.cmd import provider_save
            filename = 'map_scale_optimize_{:.2f}.pdb'.format(scale)
            provider_save(self.sim_manager.isolde.session, filename, [self.sim_manager.isolde.selected_model])
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









class SimPerformanceTracker:
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
    from openmm.app import Topology, Element
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


class UnparameterisedResiduesError(RuntimeError):
    '''
    Raised when one or more residues cannot be matched to a unique MD template.

    Carries the offending ChimeraX residues (``unmatched`` and ``ambiguous``)
    so the caller can decide how to recover (e.g. clearing a stale
    ``isolde_template_name`` override and retrying). The default message starts
    with "Unparameterised" so that the existing ``str(e).startswith(...)`` catch
    in :meth:`Isolde.start_sim` still fires the Unparameterised Residues widget
    when this propagates unhandled.
    '''
    def __init__(self, unmatched=None, ambiguous=None,
            message='Unparameterised residue detected'):
        super().__init__(message)
        self.unmatched = list(unmatched) if unmatched else []
        self.ambiguous = list(ambiguous) if ambiguous else []


def find_residue_templates(residues, forcefield, ligand_db = None, logger=None,
        clear_failed_overrides=False):
    '''
    Works out the template name applicable to cysteine residues (since OpenMM
    can't work this out for itself when ignoreExternalBonds is True), and
    looks up the template name for all known parameterised ligands from the
    CCD.

    If ``clear_failed_overrides`` is True, a residue carrying an
    ``isolde_template_name`` override that names a template *not present in the
    forcefield* will have the override reset to ``None`` (rather than merely
    ignored), so it cannot keep re-triggering the same failure. This must only
    be set on the real simulation-build path; read-only callers (the preflight
    command, the GUI display) leave it False so they never mutate the model.
    '''
    template_names = forcefield._templates.keys()
    import numpy
    templates = {}
    registered_template_names = [getattr(r, 'isolde_template_name', None) for r in residues]
    for i, tname in enumerate(registered_template_names):
        if tname is not None:
            if tname in template_names:
                templates[i] = tname
            elif clear_failed_overrides:
                if logger is not None:
                    logger.warning('Residue {} has a registered template name {}, but this is not present in the forcefield. Clearing the override and falling back to automatic template assignment.'.format(residues[i], tname))
                residues[i].isolde_template_name = None
            else:
                if logger is not None:
                    logger.warning('Residue {} has registered template name {}, but this is not present in the forcefield. Ignoring the registered template name.'.format(residues[i], tname))
    cys_indices = numpy.where(residues.names == 'CYS')[0]
    # remove any CYS residues that already have a registered template name
    cys_indices = set(cys_indices) - set(templates.keys())
    for c_i in cys_indices:
        r = residues[c_i]
        rtype = cys_type(r)
        if rtype is not None:
            templates[c_i] = rtype
            # Don't do this yet. Need to add code to clear the template name if the residue's atoms change
            # r.isolde_template_name = rtype

    # Check USER_ templates for ALL residues first, regardless of polymer type.
    # This ensures that custom templates loaded via "Load residue parameters"
    # (which are registered with the USER_ prefix by loadFile's resname_prefix
    # argument) are found even for residues that ChimeraX classifies as polymer
    # (PT_AMINO_ACID / PT_NUCLEIC) rather than PT_NONE.  Without this pre-pass,
    # terminal caps and other non-standard polymer residues fall through to
    # signature-based matching and can collide with built-in templates such as
    # ACEcyc or NHE, producing an "ambiguous template" error.
    for name in numpy.unique(residues.names):
        if name == 'HOH':
            continue
        user_name = 'USER_{}'.format(name)
        if user_name in template_names:
            indices = numpy.where(residues.names == name)[0]
            indices = set(indices) - set(templates.keys())
            for i in indices:
                templates[i] = user_name

    # filter out residues that already have templates assigned
    filter = numpy.ones(len(residues), dtype=bool)
    filter[list(templates.keys())] = False
    filtered_residues = residues[filter]
    from chimerax.atomic import Residue
    ligands = filtered_residues[filtered_residues.polymer_types == Residue.PT_NONE]
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
            # Already handled by the pre-pass above; skip.
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
    metal_atoms = atoms[numpy.isin(atoms.names, metals)]
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
            elif a.element.name == "C":
                # SG bonded to any external carbon -- a thioether.
                # Previously this branch only matched a partner atom
                # literally named "CH3" (the ACEcyc head cap), which
                # missed thioether bonds to any other external carbon
                # (covalent-inhibitor warheads, post-translational
                # modifications, designed bioconjugates, ...).  The
                # CYScyc / CCYScyc templates only depend on SG having
                # one external bond; the partner atom's name does not
                # affect the internal charges, so the broadened match
                # is safe.
                if 'OXT' in names:
                    return 'CCYScyc'
                if 'H1' in names:
                    # No NCYScyc template ships in termods.xml yet, so
                    # return CYM rather than silently mis-parameterise
                    # an N-terminal Cys with an S--C external bond.
                    return 'CYM'
                return 'CYScyc'
            # Assume metal binding - will eventually need to do something better here
            return 'CYM'
        if a.name == 'HG':
            if 'OXT' in names:
                return 'CCYS'
            if 'H1' in names:
                return 'NCYS'
            return 'CYS'

def get_available_platforms():
        from openmm import Platform
        platform_names = []
        for i in range(Platform.getNumPlatforms()):
            p = Platform.getPlatform(i)
            name = p.getName()
            platform_names.append(name)
        return platform_names
