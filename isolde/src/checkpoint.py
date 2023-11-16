# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 02-Jan-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



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
        self._restraint_data = []
        from .molobject import _RestraintMgr
        for rm in structure.all_models():
            if isinstance(rm, _RestraintMgr):
                self._restraint_data.append(rm.save_checkpoint(atoms))

        self.saved_coords = sc.all_atoms.coords


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

        for rm, data in self._restraint_data:
            if rm.deleted:
                self.session.logger.warning(f'{rm.name} manager has been deleted. Unable to restore its restraints.')
                continue
            rm.restore_checkpoint(data)


        self.sim_construct.all_atoms.coords = self.saved_coords
        if sm is not None and sm.sim_running:
            sm.sim_handler.push_coords_to_sim(self.saved_coords)
