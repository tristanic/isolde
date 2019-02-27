# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



class CheckPoint:
    '''
    Stores all the necessary information (atom positions, restraint
    parameters etc.) to revert a running simulation and the master
    molecule in ChimeraX back to a given state.
    '''
    def __init__(self, isolde):
        '''
        Save a snapshot of ISOLDE's currently running simulation. Raises a
        :class:`TypeError` if no simulation is running. The :class:`CheckPoint`
        only covers the atoms that are mobile in the current simulation, so
        becomes invalid once a new simulation is started.

        Args:
            * isolde
                - the top-level :class:`Isolde` instance
        '''
        if not isolde.simulation_running:
            raise TypeError('Checkpointing is only available when a '\
                        +'simulation is running!')
        session = self.session = isolde.session
        self.isolde = isolde
        structure = self.structure = isolde.selected_model
        sm = self.sim_manager = isolde.sim_manager
        sc = self.sim_construct = sm.sim_construct
        atoms = self.mobile_atoms = sc.mobile_atoms
        heavy_atoms = self.mobile_heavy_atoms = sc.mobile_heavy_atoms
        from . import session_extensions as sx

        pr_m = self.position_restraint_mgr = sx.get_position_restraint_mgr(structure)
        cr_m = self.chiral_restraint_mgr = sx.get_chiral_restraint_mgr(structure)
        pdr_m = self.proper_dihedral_restraint_mgr = sx.get_proper_dihedral_restraint_mgr(structure)
        dr_m = self.distance_restraint_mgr = sx.get_distance_restraint_mgr(structure)
        tug_m = self.tuggable_atoms_mgr = sx.get_tuggable_atoms_mgr(structure)
        mdff_mgrs = self.mdff_mgrs = {}
        from chimerax.map import Volume
        for v in structure.all_models():
            if isinstance(v, Volume):
                mdff_mgrs[v] = sx.get_mdff_mgr(structure, v)

        self.saved_coords = sc.all_atoms.coords

        prs = self.saved_prs = pr_m.get_restraints(atoms)
        self.saved_pr_properties = {
            'targets':              prs.targets,
            'spring_constants':     prs.spring_constants,
            'enableds':             prs.enableds
        }

        chirals = self.saved_chirals = cr_m.get_restraints_by_atoms(heavy_atoms)
        self.saved_cr_properties = {
            'spring_constants':     chirals.spring_constants,
            'enableds':             chirals.enableds,
        }

        pdrs = self.saved_pdrs = pdr_m.get_all_restraints_for_residues(atoms.unique_residues)
        self.saved_pdr_properties = {
            'targets':              pdrs.targets,
            'spring_constants':     pdrs.spring_constants,
            'enableds':             pdrs.enableds,
            'displays':             pdrs.displays,
        }

        drs = self.saved_drs = dr_m.atoms_restraints(atoms)
        self.saved_dr_properties = {
            'targets':              drs.targets,
            'spring_constants':     drs.spring_constants,
            'enableds':             drs.enableds,
        }

        # Don't save states of tuggables. Just disable any existing tugs on reversion

        # TODO: Decide if it makes sense to also checkpoint MDFF atom properties?
        # sm = self.saved_mdff_properties = {}
        # mas = self.saved_mdff_atoms = {}
        # for v, mdff_mgr in mdff_mgrs.items():
        #     props = sm[v] = {}
        #     matoms = mas[v] = mdff_mgr.get_mdff_atoms(atoms)
        #
        #     props['coupling_constants'] = matoms.coupling_constants
        #     props['enableds'] = matoms.enableds

    def revert(self):
        '''
        Revert the master construct and simulation (if applicable) to
        the saved checkpoint state. Can be called after the simulation the
        :class:`CheckPoint` was created in is stopped, but will raise a
        :class:`TypeError` if called after a new simulation is started.
        '''
        sm = self.isolde.sim_manager
        if sm is not None and sm != self.sim_manager:
            raise TypeError('A new simulation has been started since '\
                +'this checkpoint was saved. This checkpoint is no '\
                +'longer valid!')

        # Rather than deleting any restraints that have been added since the
        # checkpoint was saved, it is much more straightforward/less costly to
        # just disable all restraints, then re-enable only those that were
        # enabled at the checkpoint.
        atoms = self.mobile_atoms
        heavy_atoms = self.mobile_heavy_atoms

        pr_m = self.position_restraint_mgr
        saved_prs = self.saved_prs
        all_prs = pr_m.get_restraints(atoms)
        all_prs.enableds = False
        for prop, val in self.saved_pr_properties.items():
            setattr(saved_prs, prop, val)

        cr_m = self.chiral_restraint_mgr
        saved_chirals = self.saved_chirals
        all_chirals = cr_m.get_restraints_by_atoms(heavy_atoms)
        all_chirals.enableds = False
        for prop, val in self.saved_cr_properties.items():
            setattr(saved_chirals, prop, val)

        pdr_m = self.proper_dihedral_restraint_mgr
        saved_pdrs = self.saved_pdrs
        all_pdrs = pdr_m.get_all_restraints_for_residues(atoms.unique_residues)
        all_pdrs.enableds = False
        for prop, val in self.saved_pdr_properties.items():
            setattr(saved_pdrs, prop, val)

        dr_m = self.distance_restraint_mgr
        saved_drs = self.saved_drs
        all_drs = dr_m.atoms_restraints(atoms)
        all_drs.enableds = False
        for prop, val in self.saved_dr_properties.items():
            setattr(saved_drs, prop, val)

        atoms = self.mobile_atoms
        tugs = self.tuggable_atoms_mgr.get_tuggables(atoms)
        tugs.enableds = False


        self.sim_construct.all_atoms.coords = self.saved_coords
        if sm is not None and sm.sim_running:
            sm.sim_handler.push_coords_to_sim(self.saved_coords)
