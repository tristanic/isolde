class VolumeTrajectoryTracker:
    def __init__(self, isolde, model, clipper_volume, trajectory_volumes,
                    anneal_start_temp = 300, anneal_end_temp = 0, 
                    anneal_frames = 100, save_file_prefix='trajectory_'):
        self.isolde = isolde
        self.session = isolde.session
        self.model = model
        self.volume = clipper_volume
        self.v_traj = iter(trajectory_volumes)
        self.start_temp = anneal_start_temp
        self.end_temp = anneal_end_temp
        self.anneal_frames = anneal_frames
        self.file_prefix=save_file_prefix
        self._frame_counter = 0
        self._volume_counter = 0
        self._original_map_data = clipper_volume.data.matrix().copy()
        if self.isolde.simulation_running:
            self._sim_start_cb()
        else:
            self.isolde.triggers.add_handler(self.isolde.SIMULATION_STARTED, self._sim_start_cb)
            from chimerax.core.commands import run
            run (self.session, f'isolde sim start #{self.model.id_string}')

    def _sim_start_cb(self, *_):
        sh = self.isolde.sim_handler
        self._sim_update_handler = sh.triggers.add_handler('coord update', self._sim_update_cb)
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def _sim_update_cb(self, *_):
        fc = self._frame_counter
        sh = self.isolde.sim_handler
        sh.temperature = self.start_temp-fc/self.anneal_frames*(self.start_temp-self.end_temp)
        self._frame_counter = (self._frame_counter+1)%self.anneal_frames
        if self._frame_counter == 0:
            from chimerax.core.commands import run
            run(self.session, f'save {self.file_prefix}{self._volume_counter:03}.pdb #{self.model.id_string}')
            try:
                nv = next(self.v_traj)
                self._volume_counter += 1
            except StopIteration:
                # Finished
                self.volume.data.matrix()[:] = self._original_map_data
                self.volume.data.values_changed()
                run(self.session, 'isolde sim stop discardTo start')
                return
            self.volume.data.matrix()[:] = nv.data.matrix()
            self.volume.data.values_changed()





