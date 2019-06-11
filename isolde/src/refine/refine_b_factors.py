# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 11-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll



class B_Factor_Direct_Iso:
    '''
    Simple isotropic B-factor estimation based directly on physical mobility of
    atoms in a MDFF simulation. The model is simulated at a constant 300K until
    sufficient samples have been taken to give a good estimation of the variance
    of each atom's coordinates.
    NOTE: any aniso_u records will be cleared!
    '''
    def __init__(self, isolde, num_samples = 100, base_b_factor=20):
        self.isolde = isolde
        self.session = isolde.session
        m = self.model = isolde.selected_model
        self._num_samples = num_samples

        from copy import deepcopy
        self._original_sim_params = deepcopy(isolde.sim_params)
        isolde.sim_params.temperature = 100
        import numpy

        m.atoms.aniso_u6 = None


        self._count = 0

        from chimerax.clipper.clipper_python import Util

        self._badd = base_b_factor

        self._log_weights = numpy.array([-1, 0, 1], float)
        self._dlog_weight = 1
        self._scores = numpy.array([0,0,0], float)

        self._outer_iteration = 0
        self._max_iterations = 10
        self._best_log_weight = 0
        self._inner_iteration = 0


        m.atoms.selected = True
        isolde.start_sim()
        # Since unparameterised residues are ignored by the simulation, the
        # cohort of mobile atoms may not match the total set of atoms. We only
        # want to update B-factors for the mobile atoms
        sim_construct = isolde.sim_manager.sim_construct
        ma = self._mobile_atoms = sim_construct.mobile_atoms
        self._trajectory = numpy.empty((num_samples, *ma.coords.shape))

        isolde.sim_handler.triggers.add_handler('coord update', self._coord_update_cb)

    def _coord_update_cb(self, *_):
        if self.isolde.simulation_mode == 'min':
            return
        m = self.model
        if self._count < self._num_samples:
            coords = self._mobile_atoms.coords
            self._trajectory[self._count] = coords
            self._count += 1
        else:
            self.isolde.discard_sim(revert_to='start', warn=False)
            self._optimize_u_iso()
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER

    def _optimize_u_iso(self):
        import numpy
        traj = self._trajectory
        variance = traj.var(axis=0).sum(axis=1)
        from chimerax.clipper.clipper_python import Util
        self._u_base = numpy.array([Util.u2b(u) for u in variance]).astype(numpy.float32)
        from chimerax.clipper.symmetry import get_map_mgr
        map_mgr = get_map_mgr(self.model)
        xmapset = self._xmapset = map_mgr.xmapsets[0]
        xmapset.triggers.add_handler('maps recalculated', self._map_update_cb)
        self._map_update_cb()

    def _map_update_cb(self, *_):
        m = self.model
        xmapset = self._xmapset
        scores = self._scores
        logweights = self._log_weights
        inner = self._inner_iteration
        outer = self._outer_iteration
        if (inner > 0):
            scores[inner-1] = xmapset.rwork
        if self._outer_iteration < self._max_iterations:
            if (inner == 3):
                import operator
                min_index, min_val = min(enumerate(scores), key=operator.itemgetter(1))
                scores[:] = 0
                scores[1] = min_val
                best_weight = self._best_log_weight = logweights[min_index]
                logweights[:] = 0
                logweights[1] = best_weight
                self._dlog_weight *= 0.7
                logweights[0] = best_weight-self._dlog_weight
                logweights[2] = best_weight+self._dlog_weight
                self._inner_iteration = 0
                self._outer_iteration += 1
                self._map_update_cb()
                return
            if scores[inner] == 0:
                self._mobile_atoms.bfactors = self._badd+self._u_base*10**logweights[inner]
                self._inner_iteration += 1
            else:
                self._inner_iteration += 1
                self._map_update_cb()
        else:
            final_bfactors = bfactors = self._badd+self._u_base*10**(self._best_log_weight)
            final_bfactors[final_bfactors>999] = 999
            self._mobile_atoms.bfactors = final_bfactors
            self.isolde.sim_params = self._original_sim_params
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER




class B_Factor_Direct_Aniso:
    '''
    Anisotropic B-factor estimation based directly on physical mobility of
    atoms in a MDFF simulation. The model is simulated at a constant 300K until
    sufficient samples have been taken to give a good estimation of the variance
    of each atom's coordinates. Then, for each atom a covariance matrix is
    calculated based on the sampled coords, which is combined with an overall
    weighting factor to give estimated aniso_u values. The weighting factor is
    then optimised to minimise Rwork using a simple Newton method implementation.
    '''
    def __init__(self, isolde, num_samples = 100):
        self.isolde = isolde
        self.session = isolde.session
        m = self.model = isolde.selected_model

        from copy import deepcopy
        self._original_sim_params = deepcopy(isolde.sim_params)
        isolde.sim_params.temperature = 100

        import numpy
        # Prepare a container to hold the sampled coordinates
        self._trajectory = numpy.empty((num_samples, len(m.atoms.coords), 3))
        self._current_frame = 0
        self._num_samples = num_samples

        self._log_weight = 0
        self._dlog_weight = 1
        self._log_weights = numpy.array([0,1,2], float)
        self._scores = numpy.array([0,0,0], float)
        self._outer_iteration = 0
        self._max_iterations = 10
        self._best_log_weight = 0
        self._inner_iteration = 0

        m.atoms.selected = True
        isolde.start_sim()
        isolde.sim_handler.triggers.add_handler('coord update', self._coord_update_cb)


    def _coord_update_cb(self, *_):
        if self.isolde.simulation_mode == 'min':
            return
        if self._current_frame < self._num_samples:
            self._trajectory[self._current_frame] = self.model.atoms.coords
            self._current_frame += 1
        else:
            self.isolde.discard_sim(revert_to='start', warn=False)
            self._optimize_anisou()
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER

    def _optimize_anisou(self):
        import numpy
        from numpy import cov
        natoms = len(self.model.atoms)
        covariances = self._covariances = numpy.empty((natoms,3,3))
        traj = self._trajectory
        for i in range(natoms):
            coords = traj[:,i,:]
            covariances[i] = numpy.abs(cov(coords, rowvar=False))
        #tensors = self._tensors = numpy.empty(covariances.shape)
        # for cv, t in zip(covariances, tensors):
        #     t[:] = numpy.abs(numpy.linalg.inv(cv))

        anisou = self._anisou_base = numpy.empty((natoms, 6))
        anisou[:,0] = covariances[:,0,0]
        anisou[:,1] = covariances[:,1,1]
        anisou[:,2] = covariances[:,2,2]
        anisou[:,3] = covariances[:,0,1]
        anisou[:,4] = covariances[:,0,2]
        anisou[:,5] = covariances[:,1,2]
        # Convert from B to U
        from math import pi
        anisou[:] = numpy.sqrt(anisou/(8*pi**2))
        from chimerax.clipper.symmetry import get_map_mgr
        map_mgr = get_map_mgr(self.model)
        xmapset = self._xmapset = map_mgr.xmapsets[0]
        xmapset.triggers.add_handler('maps recalculated', self._map_update_cb)
        self._map_update_cb()

    def _map_update_cb(self, *_):
        m = self.model
        xmapset = self._xmapset
        scores = self._scores
        logweights = self._log_weights
        inner = self._inner_iteration
        outer = self._outer_iteration
        if (inner > 0):
            scores[inner-1] = xmapset.rwork
        if self._outer_iteration < self._max_iterations:
            if (inner == 3):
                import operator
                min_index, min_val = min(enumerate(scores), key=operator.itemgetter(1))
                scores[:] = 0
                scores[1] = min_val
                best_weight = self._best_log_weight = logweights[min_index]
                logweights[:] = 0
                logweights[1] = best_weight
                self._dlog_weight *= 0.7
                logweights[0] = best_weight-self._dlog_weight
                logweights[2] = best_weight+self._dlog_weight
                self._inner_iteration = 0
                self._outer_iteration += 1
                self._map_update_cb()
                return
            if scores[inner] == 0:
                print("Trying weight = 10**{}".format(logweights[inner]))
                m.atoms.aniso_u6 = self._anisou_base*10**(logweights[inner])
                self._inner_iteration += 1
            else:
                self._inner_iteration += 1
                self._map_update_cb()
        else:
            m.atoms.aniso_u6 = self._anisou_base*10**(self._best_log_weight)
            self.isolde.sim_params = self._original_sim_params
            from chimerax.core.triggerset import DEREGISTER
            return DEREGISTER
