# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



from .sim_param_mgr import SimParams
from .openmm_interface import SimManager, SimHandler, OpenmmThreadHandler, \
                              SimConstruct

class SimSaver:
    def __init__(self, session, sim_handler, file_prefix='trajectory', interval=100):
        self.session = session
        self.sim_handler = sim_handler
        self.file_prefix=file_prefix
        self.interval = interval
        self.handler = sim_handler.triggers.add_handler('coord update', self._sim_coord_update_cb)
        self._counter=0
        self._frame_count = 0
    def _sim_coord_update_cb(self, *_):
        self._counter += 1
        if self._counter % self.interval == 0:
            from chimerax.core.commands import run
            run(self.session, f'save {self.file_prefix}_{self._frame_count:04}.pdb #{self.session.isolde.selected_model.id_string}')
            self._frame_count += 1

class MorphFollower:
    def __init__(self, session, sim_handler, model, target, spring_constant=20, interval=100):
        self.session = session
        self.model = model
        self.sim_handler = sim_handler
        self.target = target
        self.spring_constant = spring_constant
        self.interval = interval
        self._counter = 0
        self.frames = iter(target.coordset_ids)

        from chimerax.isolde import session_extensions as sx
        prm = sx.get_position_restraint_mgr(model)
        prs = prm.add_restraints(model.atoms)
        target_residues = target.residues
        target.active_coordset_id = next(self.frames)
        pr_dict = self._pr_to_target_atom = {}
        for pr in prs:
            a = pr.atom
            r = a.residue
            trc = target_residues[target_residues.chain_ids==r.chain_id]
            tr = trc[trc.numbers==r.number]
            if len(tr) != 1:
                continue
            tr = tr[0]
            ta = tr.find_atom(a.name)
            if ta is None:
                continue
            pr.target = ta.scene_coord
            pr.spring_constant = self.spring_constant
            pr.enabled=True
            pr_dict[pr] = ta
        self.handler = sim_handler.triggers.add_handler('coord update', self._sim_coord_update_cb)
    
    def _sim_coord_update_cb(self, *_):
        self._counter += 1
        if self._counter % self.interval == 0:
            try:
                self.target.active_coordset_id = next(self.frames)
            except StopIteration:
                self.session.logger.info('ISOLDE: morph completed')
                from chimerax.core.triggerset import DEREGISTER
                return DEREGISTER
            for pr, ta in self._pr_to_target_atom.items():
                pr.target = ta.scene_coord
     

            

