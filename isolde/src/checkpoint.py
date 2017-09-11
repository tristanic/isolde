# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


class CheckPoint:
    '''
    Stores all the necessary information (atom positions, restraint 
    parameters etc.) to revert a running simulation and the master 
    molecule in ChimeraX back to a given state.
    '''
    def __init__(self, isolde):
        if not isolde.simulation_running:
            raise TypeError('Checkpointing is only available when a '\
                        +'simulation is running!')
        self.isolde = isolde
        self.sim_interface = isolde._sim_interface
        sc = self.sim_construct = isolde._total_sim_construct
        tm = self.mobile = isolde._total_mobile
        bd = self.backbone_dihedrals = isolde._mobile_backbone_dihedrals
        
        self.residues = tm.unique_residues
        rot = self.rotamers = {r: isolde.rotamers.get(r, None) for r in self.residues}
        dd = self.distance_restraints = isolde._sim_distance_restraint_dict
        pr = self.position_restraints = isolde._sim_pos_restr
        
        self.saved_positions = sc.coords
        
        sd = self.saved_dihedral_target_states = dict()
        for d_type in (bd.phi, bd.psi, bd.omega):
            sd[d_type] = (d_type.targets,
                          d_type.spring_constants
                         )
        
        sr = self.saved_rotamer_target_states = dict()
        for res, r in rot.items():
            if r is not None:
                sr[r] = (r.restrained, r.target)
        
        sdr = self.saved_distance_restraints = dict()
        for dr_name, dr in dd.items():
            sdr[dr] = (dr.targets, dr.spring_constants)
        
        spr = self.saved_position_restraints = (
            pr.targets,
            pr.spring_constants
            )
        
    def revert(self):
        '''
        Revert the master construct and simulation (if applicable) to 
        the saved checkpoint state.
        '''
        si = self.isolde.sim_interface
        if si is not None:
            if si != self.sim_interface:
                raise TypeError('A new simulation has been started since '\
                    +'this checkpoint was saved. This checkpoint is no '\
                    +'longer valid!')
        
        currently_paused = si.paused
        if not currently_paused:
            # Make sure the simulation is not running while we change things
            si.pause_and_acknowledge()
        
        self.sim_construct.coords = self.saved_positions
        si.update_coords(self.sim_construct)
        
        sd = self.saved_dihedral_target_states
        bd = self.backbone_dihedrals
        for d_type in (bd.phi, bd.psi, bd.omega):
            d_type.targets = sd[d_type][0]
            d_type.spring_constants = sd[d_type][1]
            si.update_dihedral_restraints(d_type)
        
        for rot, saved in self.saved_rotamer_target_states.items():
            rot.restrained = saved[0]
            rot.target = saved[1]
            si.update_rotamer_target(rot)
        
        dd = self.distance_restraints
        for dr_name, dr in dd.items():
            sdr = self.saved_distance_restraints[dr]
            dr.targets = sdr[0]
            dr.spring_constants = sdr[1]
            si.update_distance_restraints(dr)
        
        pr = self.position_restraints
        spr = self.saved_position_restraints
        pr.targets = spr[0]
        pr.spring_constants = spr[1]
        si.update_position_restraints(pr)
        
        if not currently_paused:
            si.paused = False
            
        
        
        
