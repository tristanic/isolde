import numpy
from random import choice, random
class SimTester:
    def __init__(self, session, model, selected_atoms):
        from .. import session_extensions
        self.session = session
        pd_m = self.proper_dihedral_mgr = session_extensions.get_proper_dihedral_mgr(session)
        rama_m = self.rama_mgr = session_extensions.get_ramachandran_mgr(session)
        rota_m = self.rota_mgr = session_extensions.get_rotamer_mgr(session)

        rama_a = self.rama_annotator = session_extensions.get_rama_annotator(model)
        rota_a = self.rota_annotator = session_extensions.get_rota_annotator(model)

        pdr_m = self.proper_dihedral_restraint_mgr = session_extensions.get_proper_dihedral_restraint_mgr(model)
        pr_m = self.position_restraint_mgr = session_extensions.get_position_restraint_mgr(model)
        ta_m = self.tuggable_atoms_mgr = session_extensions.get_tuggable_atoms_mgr(model)
        dr_m = self.distance_restraints_mgr = session_extensions.get_distance_restraint_mgr(model)

        mobile_atoms = selected_atoms.residues.atoms
        fixed_atoms = get_shell_of_residues(model, mobile_atoms, 8).residues.atoms
        all_sim_atoms = mobile_atoms.merge(fixed_atoms)
        all_atoms = model.residues.atoms
        mob_i = all_atoms.indices(mobile_atoms)
        as_i = all_atoms.indices(all_sim_atoms)
        mob_i.sort()
        as_i.sort()
        mobile_atoms = all_atoms[mob_i]
        all_sim_atoms = all_atoms[as_i]

        mobile_res = mobile_atoms.unique_residues
        rama_a.restrict_to_selected_residues(mobile_res)
        rota_a.restrict_to_selected_residues(mobile_res)

        dihedral_restraints = pdr_m.add_all_defined_restraints_for_residues(mobile_atoms.unique_residues)

        position_restraints = pr_m.add_restraints(mobile_atoms)
        tuggable_atoms = ta_m.add_tuggables(mobile_atoms)

        distance_restraints = []
        nonh = all_sim_atoms[all_sim_atoms.element_names != 'H']
        for i in range(20):
            try:
                a1 = choice(nonh)
                a2 = choice(nonh)
                d = dr_m.add_restraint(a1,a2)
                d.target = numpy.linalg.norm(a2.coord-a1.coord)
                d.spring_constant = random()*5000
                d.enabled = True
                distance_restraints.append(d)
            except:
                continue
        from ..molarray import Distance_Restraints
        distance_restraints = Distance_Restraints(distance_restraints)

        from ..openmm import openmm_interface, sim_param_mgr
        sim_construct = self.sim_construct = openmm_interface.Sim_Construct(all_sim_atoms, mobile_atoms, fixed_atoms)
        sim_params = self.params = sim_param_mgr.SimParams()

        sim_handler = self.sim_handler = openmm_interface.Sim_Handler(session, sim_params, sim_construct, 100.0)

        sim_handler.initialize_custom_forces()

        ramas = rama_m.get_ramas(mobile_atoms.unique_residues)
        ramas = ramas[ramas.valids]
        sim_handler.add_amber_cmap_torsions(ramas)
        sim_handler.add_dihedral_restraints(dihedral_restraints)
        sim_handler.add_distance_restraints(distance_restraints)
        sim_handler.add_position_restraints(position_restraints)
        sim_handler.add_tuggables(tuggable_atoms)

        uh = self._update_handlers = []
        uh.append((pdr_m, pdr_m.triggers.add_handler('changes', self._pdr_changed_cb)))
        uh.append((dr_m, dr_m.triggers.add_handler('changes', self._dr_changed_cb)))
        uh.append((pr_m, pr_m.triggers.add_handler('changes', self._pr_changed_cb)))
        uh.append((ta_m, ta_m.triggers.add_handler('changes', self._tug_changed_cb)))

    def start_sim(self):
        sh = self.sim_handler
        sh.start_sim()
        sh.triggers.add_handler('sim terminated', self._pr_sim_end_cb)
        sh.triggers.add_handler('sim terminated', self._dr_sim_end_cb)
        sh.triggers.add_handler('sim terminated', self._pdr_sim_end_cb)
        sh.triggers.add_handler('sim terminated', self._tug_sim_end_cb)
        sh.triggers.add_handler('sim terminated', self._rama_a_sim_end_cb)
        sh.triggers.add_handler('sim terminated', self._rota_a_sim_end_cb)


    def stop_sim(self):
        for mgr, handler in self._update_handlers:
            mgr.triggers.remove_handler(handler)
        self._update_handlers = []
        self.sim_handler.stop()

    def _rama_a_sim_end_cb(self, *_):
        self.rama_annotator.track_whole_model = True
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def _rota_a_sim_end_cb(self, *_):
        self.rota_annotator.track_whole_model = True
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER

    def _pr_changed_cb(self, trigger_name, changes):
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
            all_changeds = concatenate(changeds)
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
            all_changeds = concatenate(changeds)
            all_changeds = all_changeds[all_changeds.sim_indices != -1]
            self.sim_handler.update_distance_restraints(all_changeds)

    def _dr_sim_end_cb(self, *_):
        restraints = self.distance_restraints_mgr.intra_restraints(self.sim_construct.all_atoms)
        restraints.clear_sim_indices()
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER


    def _pdr_changed_cb(self, trigger_name, changes):
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
            all_changeds = concatenate(changeds)
            all_changeds = all_changeds[all_changeds.sim_indices != -1]
            self.sim_handler.update_dihedral_restraints(all_changeds)

    def _pdr_sim_end_cb(self, *_):
        restraints = self.proper_dihedral_restraint_mgr.get_all_restraints_for_residues(self.sim_construct.all_atoms.unique_residues)
        restraints.clear_sim_indices()
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER


    def _tug_changed_cb(self, trigger_name, changes):
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
            all_changeds = concatenate(changeds)
            all_changeds = all_changeds[all_changeds.sim_indices != -1]
            self.sim_handler.update_tuggables(all_changeds)

    def _tug_sim_end_cb(self, *_):
        tuggables = self.tuggable_atoms_mgr.get_tuggables(self.sim_construct.all_atoms)
        tuggables.clear_sim_indices()
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER






def get_shell_of_residues(model, existing_sel, dist_cutoff):
    from chimerax.core.geometry import find_close_points
    from chimerax.core.atomic import selected_atoms, Atoms, concatenate
    selatoms = existing_sel
    allatoms = model.atoms
    unselected_atoms = allatoms.subtract(selatoms)
    selcoords = selatoms.coords
    unselcoords = unselected_atoms.coords
    ignore, shell_indices = find_close_points(selcoords, unselcoords, dist_cutoff)
    shell_atoms = unselected_atoms[shell_indices].unique_residues.atoms
    return shell_atoms
