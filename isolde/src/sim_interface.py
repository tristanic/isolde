# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


# vim: set expandtab shiftwidth=4 softtabstop=4:

import numpy
import os
import multiprocessing as mp
import ctypes
from math import pi, radians, degrees, cos
from simtk import unit, openmm
from simtk.unit import Quantity
from simtk.openmm import app
from simtk.openmm.openmm import CustomBondForce, CustomExternalForce, \
                                CustomCompoundBondForce, CustomTorsionForce, \
                                NonbondedForce
from simtk.openmm.openmm import Continuous1DFunction, Continuous3DFunction 
from simtk.openmm.app.internal import customgbforces   
from chimerax.core.atomic import concatenate
from .custom_forces import LinearInterpMapForce, TopOutBondForce, \
                           TopOutRestraintForce, FlatBottomTorsionRestraintForce, \
                           GBSAForce, AmberCMAPForce

from .constants import defaults

class SimParams:
    '''
    Container for all the parameters needed to initialise a simulation
    '''
    _default_params = {
            'dihedral_restraint_cutoff_angle': (defaults.DIHEDRAL_RESTRAINT_CUTOFF, unit.degrees),
            'restraint_max_force': (defaults.MAX_RESTRAINT_FORCE, unit.kilojoule_per_mole/unit.angstrom**2),
            'rotamer_restraint_cutoff_angle': (defaults.ROTAMER_RESTRAINT_CUTOFF, unit.degrees),
            'rotamer_spring_constant': (defaults.ROTAMER_SPRING_CONSTANT, unit.kilojoule_per_mole/unit.radians**2),
            'peptide_bond_spring_constant': (defaults.PEPTIDE_SPRING_CONSTANT, unit.kilojoule_per_mole/unit.radians**2),
            'cis_peptide_bond_cutoff_angle': (defaults.CIS_PEPTIDE_BOND_CUTOFF, unit.degrees),
            'tug_hydrogens': (False, 'bool'),
            'hydrogens_feel_maps': (False, 'bool'),
            }
    def __init__(self, **kw):
        self.params = self._default_params
        valid_keys = self.params.keys()
        for key, val in kw.items():
            try:
                self.params[key]
            except KeyError:
                raise KeyError('Invalid parameter!')
            self.params[key][0] = val
        for key, item in self.params.items():
            if type(item[1]) == unit.unit.Unit:
                self.params[key] = item[0]*item[1]
            else:
                self.params[key] = item[0]
        
    _init_docstring = '''
        Initialise the simulation parameters, specifying any non-default
        values.\n'''
    param_list = ''
    for key, val in _default_params.items():
        param_list += '\t@param {:>} \t ({} {})\n'.format(key, val[0], val[1]).expandtabs(20)
    __init__.__doc__ = _init_docstring+param_list
    
    

def start_sim_thread(sim_params, all_atoms, fixed_flags,  
                     backbone_dihedrals, rotamers, distance_restraints,
                     position_restraints, 
                     density_maps = None):
    '''
    Start an OpenMM simulation in a separate Python instance, with 
    communication via shared variables. Returns a dict containing all
    variables to be used for communication/control, with descriptive
    names.
        @param  sim_params:
            An instance of SimParams.params defining the basic simulation
            settings
        @param all_atoms:
            The set of all ChimeraX atoms defining the complete simulation
            construct
        @param fixed_flags:
            A numpy boolean array flagging atoms to be fixed in space.
        @param backbone_dihedrals:
            An ISOLDE Backbone_Dihedrals object defining all phi, psi and
            omega dihedrals in the simulation
        @param rotamers
            An ISOLDE Rotamers object
        @param distance_restraints
            A {name: object} dict of Distance_Restraints objects.
        @param maps
            A tuple of ISOLDE map objects
    '''
    manager = mp.Manager()
    
    comms_object = {
        'error'                     : mp.Queue(),
        'status'                    : mp.Queue(),
        'initialization complete'   : mp.Value(ctypes.c_bool, False),
        'pause'                     : mp.Value(ctypes.c_bool, False),
        'unstable'                  : mp.Value(ctypes.c_bool, False),
        'sim mode'                  : mp.Value('i', 0),
        'constants'                 : manager.Namespace(),
        'rotamer map'               : manager.dict(),
        'rotamer targets'           : manager.dict(),
        'position restraints map'   : manager.dict(),
    }
    
    
    
    fixed_indices = numpy.where(fixed_flags)
    mobile_indices = numpy.where(numpy.invert(fixed_flags))
    fixed_atoms = all_atoms[fixed_indices]
    mobile_atoms = all_atoms[mobile_indices]
    residues = all_atoms.residues
    atom_names    = all_atoms.names
    element_names = all_atoms.element_names
    residue_names = residues.names
    residue_nums  = residues.numbers
    chain_ids     = residues.chain_ids 
    coords        = all_atoms.coords
    
    bonds               = all_atoms.intra_bonds
    bond_atoms          = bonds.atoms
    bond_atom_indices   = (all_atoms.indices(ba) for ba in bond_atoms)
    
    # Secondary structure and peptide bond backbone restraints

    bd = backbone_dihedrals
    phi = bd.phi
    n_phi = len(phi)
    phi_i = numpy.reshape(all_atoms.indices(phi.atoms), [n_phi,4])
    comms_object['phi sim indices'] = mp.Array('i', n_phi)
    
    
    psi = bd.psi
    n_psi = len(psi)
    psi_i = numpy.reshape(all_atoms.indices(psi.atoms), [n_psi,4])
    comms_object['psi sim indices'] = mp.Array('i', n_psi)
    
    omega = bd.omega
    n_omega = len(omega)
    omega_i = numpy.reshape(all_atoms.indices(omega.atoms), [n_omega,4])
    comms_object['omega sim indices'] = mp.Array('i', n_omega)
    
    cis_offset = sim_params['cis_peptide_bond_cutoff_angle']/unit.radians
    cr = (-cis_offset, cis_offset)
    o_vals = omega.vals
    omega_targets = numpy.logical_or(
                o_vals > cr[1], o_vals < cr[0]).astype(float) * pi
    o_t = comms_object['omega targets'] = mp.Array('f', n_omega)
    o_t[:] = omega_targets

   
    # CMAP corrections only apply to residues that have both phi and
    # psi dihedrals
    phi_has_psi = phi.residues.indices(psi.residues)
    psi_has_phi = psi.residues.indices(phi.residues)
    good_phi = phi[phi_has_psi[phi_has_psi != -1]]
    good_psi = psi[psi_has_phi[psi_has_phi != -1]]
    phi_cmap_resnames = good_phi.residues.names
    psi_cmap_resnames = good_psi.residues.names
    n = len(good_psi)
    phi_cmap_indices = all_atoms.indices(good_phi.atoms).reshape((n,4))
    psi_cmap_indices = all_atoms.indices(good_psi.atoms).reshape((n,4))
    
    # Rotamers are tricky. Each rotamer holds a set of dihedrals which
    # must individually be targeted in the simulation.
    input_rotamer_map = {}
    input_rotamer_targets = comms_object['rotamer targets']
    mobile_res = mobile_atoms.unique_residues
    for i, r in enumerate(mobile_res):
        try:
            rot= rotamers[res]
        except KeyError:
            # Non-rotameric residue
            continue
        dlist = rot.dihedrals
        atoms = dlist.atoms
        indices = all_atoms.indices(atoms).reshape(len(dlist),4)
        input_rotamer_map[i] = indices
        
        if rot.restrained:
            input_rotamer_targets[i] = rot.target
        else:
            input_rotamer_targets[i] = None
    
    
    # Distance restraints. These are handled as a {name: object} dict of
    # different types of restraints (e.g. CA-CA+2, O-N+4, etc.)
    distance_restraint_keys = []
    distance_restraint_indices = {}    
    for key, d_r in distance_restraints.items():
        n_restraints = len(d_r)
        distance_restraint_keys.append(key)
        force_map_key = key + ' map'
        force_map = comms_object[force_map_key] = manager.dict()
        atom_index_key = key + ' atom indices'
        target_key = key + ' targets'
        k_key = key + 'k'
        t_array = comms_object[target_key] = mp.Array('f', n_restraints)
        k_array = comms_object[k_key] = mp.array('f', n_restraints)
        i_array = distance_restraint_indices[atom_index_key] = []
        for i,r in enumerate(d_r):
            t_array[i] = r.target_distance
            k_array[i] = r.spring_constant
            i_array.append(all_atoms.indices(d_r.atoms).tolist())
            
    # Position restraints. One per heavy atom
    pr_atoms = position_restraints.atoms
    pr_indices = all_atoms.indices(pr_atoms)
    n_pr = len(pr_indices)
    pr_ks = comms_object['position restraint spring constants'] = mp.Array('f', n_pr)
    pr_ks[:] = position_restraints.spring_constants
    pr_targets_base = comms_object['position restraint targets'] = mp.Array('f', n_pr*3)
    
    pr_targets = numpy.frombuffer(pr_targets.get_obj()).reshape((n_pr,3))
    pr_targets[:] = position_restraints.targets
    
    # Density maps. Even though updating these during a simulation is not
    # yet possible, we want to leave the option open.
    density_map_names = []
    density_map_transforms = {}
    if density_maps is not None:
        for key, imap in density_maps.items():
            density_map_names.append(imap.get_name())
            volume_data_key = name + ' data'
            transform_key = name + ' transform'
            vd, r = imap.crop_to_selection(atoms, pad, normalize)
            
            tf = r.xyz_to_ijk_transform.matrix
            densit_map_transforms[transform_key] = tf
            
            vd_comms_base = comms_object[volume_data_key] = mp.Array('f', vd.size)
            vd_comms_base[:] = vd.ravel('C')
        
            
    
    
    
    sim_data = {
        'fixed indices':                    fixed_indices,
        'mobile indices':                   mobile_indices,
        'atom names':                       atom_names,
        'element names':                    element_names,
        'residue names':                    residue_names,
        'residue numbers':                  residue_nums,
        'chain ids':                        chain_ids,
        'coords':                           coords,
        'bonded atom indices':              bond_atom_indices,
        'phi atom indices':                 phi_i,
        'psi atom indices':                 psi_i,
        'omega atom indices':               omega_i,
        'phi cmap resnames':                phi_cmap_resnames,
        'phi cmap indices':                 phi_cmap_indices,
        'psi cmap resnames':                psi_cmap_resnames,
        'psi cmap indices':                 psi_cmap_indices,
        'rotamer map':                      input_rotamer_map,
        'rotamer targets':                  input_rotamer_targets,
        'distance restraint keys':          distance_restraint_keys,
        'position restraint indices':       pr_indices,
        'density map names':                density_map_names,
    }
    sim_data.update(distance_restraint_indices)
    sim_data.update(density_map_transforms)

    return comms_object
    
def _sim_thread(sim_params, sim_data, comms_object
                        
                    ):
    '''
    Start a Python multiprocessing thread to hold a running OpenMM 
    simulation. The simulation needs to be entirely self-contained in
    this thread, so inputs should be only generic Python objects (e.g.
    Numpy arrays) and multiprocessing communication variables.
    '''
    try:
        params = sim_params.params
        sh = SimHandler()
        
        
        # Just initialize the forces. They're empty at this stage.
        sh.initialize_dihedral_restraint_force(params['dihedral_restraint_cutoff_angle'][0])
        sh.initialize_amber_cmap_force()
        sh.initialize_distance_restraints_force(params['restraint_max_force'][0])
        sh.initialize_position_restraints_force(params['restraint_max_force'][0])
        sh.initialize_tugging_force()
        
        
        # Backbone dihedral restraints
        phi_sim_i = comms_object['phi sim indices']
        for i, atom_indices in enumerate(sim_data['phi atom indices']):
            phi_sim_i[i] = sh.initialize_dihedral_restraint(atom_indices)
        
        psi_sim_i = comms_object['psi sim indices']
        with psi_sim_i.get_lock():
            for i, atom_indices in enumerate(sim_data['psi atom indices']):
                psi_sim_i[i] = sh.initialize_dihedral_restraint(atom_indices)
        
        omega_sim_i = comms_object['omega sim indices']
        omega_targets = comms_object['omega targets']
        with omega_sim_i.get_lock(), omega_targets.get_lock():
            for i, (atom_indices, target) in enumerate(zip(sim_data['omega atom indices'], omega_targets)):
                si = omega_sim_i[i] = sh.initialize_dihedral_restraint(atom_indices)
                sh.update_dihedral_restraint(si, target = target, 
                            k = sim_params['peptide_bond_spring_constant'])
        
        
        
        
        # CMAP corrections
        sh._add_amber_cmap_torsions(sim_data['phi cmap resnames'],
                                    sim_data['psi cmap resnames'],
                                    sim_data['phi cmap indices'],
                                    sim_data['psi cmap indices'])
        
        # Rotamers
        in_rotamer_map = sim_data['rotamer map']
        rotamer_targets = comms_object['rotamer targets']
        out_rotamer_map = comms_object['rotamer map']
        rotamer_restraint_cutoff = sim_params['rotamer_restraint_cutoff_angle'] / unit.radians
        rotamer_restraint_k = sim_params['rotamer_spring_constant'] / (unit.kilojoule_per_mole/unit.radians**2)
        with rotamer_targets.get_lock(), out_rotamer_map.get_lock():
            for key, indices in in_rotamer_map.items():
                force_indices = out_rotamer_map[key] = []
                for d_indices in indices:
                    force_indices.append(sh.initialize_dihedral_restraint(d_indices, rotamer_restraint_cutoff))
                targets = rotamer_targets[key]
                if targets is not None:
                    for (fi, di, t) in zip(force_indices, indices, targets):
                        sh.update_dihedral_restraint(fi, target = t, 
                                k = sim_params['rotamer_spring_constant']) 
        
        
        # Distance restraints
        distance_restraint_keys = sim_data['distance restraint keys']
        drf = sh._distance_restraints_force    
        for key in distance_restraint_keys:
            force_map = comms_object[key + ' map']
            t_array = comms_object[key + ' targets']
            k_array = comms_object[key + ' k']
            i_array = sim_data[key + ' atom indices']
            with force_map.get_lock(), t_array.get_lock(), k_array.get_lock(): 
                for (indices, t, k) in zip(i_array, t_array, k_array):
                    force_map[i] = drf.addBond(*indices, (k, t))
            
        #Position restraints
        pr_force = sh._position_restraints_force
        pr_force_map = comms_object['position restraints map']
        pr_indices =    sim_data['position restraint indices']
        pr_ks =         comms_object['position restraint spring constants']
        pr_targets_base = comms_object['position restraint targets']
        pr_targets = numpy.frombuffer(pr_targets.get_obj()).reshape((n_pr,3))
        
        
        with pr_ks.get_lock(), pr_targets_base.get_lock():
            for i, (index, k, target) in enumerate(zip(pr_indices, pr_ks, pr_targets)):
                 pr_force_map[i] = pr_force.addParticle(index, (k, target))
        
        
        
        top, templates = sh._openmm_topology_and_external_forces(sim_data,
                            tug_hydrogens = params['tug_hydrogens'],
                            hydrogens_feel_maps = params['hydrogens_feel_maps']
                            )
    except Exception as e:
        comms_object['error'].put(e)
        return
        
    
    



def openmm_topology_from_model(model):
    '''
    Take an AtomicStructure model from ChimeraX and return an OpenMM
    topology (e.g. for use with the OpenMM Modeller class).
    '''
    
    
    a = model.atoms
    b = model.bonds
    n = len(a)
    r = a.residues
    aname = a.names
    ename = a.element_names
    rname = r.names
    rnum = r.numbers
    cids = r.chain_ids
    from simtk.openmm.app import Topology, Element
    from simtk import unit
    top = Topology()
    cmap = {}
    rmap = {}
    atoms = {}
    for i in range(n):
        cid = cids[i]
        if not cid in cmap:
            cmap[cid] = top.addChain()   # OpenMM chains have no name
        rid = (rname[i], rnum[i], cid)
        if not rid in rmap:
            res = rmap[rid] = top.addResidue(rname[i], cmap[cid])
        element = Element.getBySymbol(ename[i])
        atoms[i] = top.addAtom(aname[i], element,rmap[rid])

    a1, a2 = b.atoms
    for i1, i2 in zip(a.indices(a1), a.indices(a2)):
        if -1 not in [i1, i2]:
            top.addBond(atoms[i1],  atoms[i2])
        
    return top

def cys_type(residue):
    sulfur_atom = residue.atoms.filter(residue.atoms.names == 'SG')[0]
    bonds = sulfur_atom.bonds
    if len(bonds) == 1:
        # Deprotonated
        return 'CYM'
    bonded_atoms = concatenate(bonds.atoms)
    for a in bonded_atoms:
        if a.residue != residue:
            return 'CYX'
        if a.name == 'HG':
            return 'CYS'

amber14 = ['amberff14SB.xml','tip3p_standard.xml',
            'tip3p_HFE_multivalent.xml', 'tip3p_IOD_multivalent.xml']

cwd = os.path.dirname(os.path.abspath(__file__))
class available_forcefields():
    '''Main force field files'''
    main_files = (
        [os.path.join(cwd,'amberff',f) for f in amber14],
        ['amber99sbildn.xml',],
        ['amber99sbnmr.xml',],
        ['amber10.xml',],
        ['charmm36.xml',],
        )
    
    main_file_descriptions = (
        'AMBER14 with improved backbone & sidechain torsions',
        'AMBER99 with improved backbone & sidechain torsions',
        'AMBER99 with modifications to fit NMR data',
        'AMBER10',
        'CHARMM36',
        )
    
    # Implicit solvent force field files. The main file and the implicit
    # solvent file must match.
    implicit_solvent_files = [
        None,
        'amber99_obc.xml',
        'amber99_obc.xml',
        'amber10_obc.xml',
        None
        ]
    
    # Explicit water models
    explicit_water_files = [
        None,
        'tip3pfb.xml',
        'tip4pfb.xml',
        'tip3p.xml',
        ]
    
    explicit_water_descriptions = [
        None,
        'TIP3P-FB (DOI: 10.1021/jz500737m)',
        'TIP4P-FB (DOI: 10.1021/jz500737m)',
        'Original TIP3P water (not recommended)',
        ]

def get_available_platforms():
        from simtk.openmm import Platform
        platform_names = []
        for i in range(Platform.getNumPlatforms()):
            p = Platform.getPlatform(i)
            name = p.getName()
            platform_names.append(name)
        return platform_names


        
class GenericMapObject:
    def __init__(self, c3dfunc, per_atom_coupling, coupling_k):
        '''
        A minimal class to hold a Continuous3DFunction and its parameters.
        '''
        self._potential_function = c3dfunc
        self._per_atom_coupling = per_atom_coupling
        self._k = coupling_k
        
    def get_potential_function(self):
        return self._potential_function
    
    def get_per_atom_coupling_params(self):
        return self._k
    
    def per_atom_coupling(self):
        return hasattr(self._k, '__len__') or hasattr(self._k, '__iter__')





class SimHandler():
    
    def __init__(self, session = None):
        self.session = session
        # Forcefield used in this simulation
        self._forcefield = None
        # Overall simulation topology
        self._topology = None
        # Atoms in simulation topology
        self._atoms = None
        # Atoms in ChimeraX construct
        self._chimerax_atoms = None
        
        # Dict holding each custom external force object, its global and its
        # per-particle parameters. Format:
        # {'name': [object, 
        #           {'global_name': value}, 
        #           [per_particle_names],
        #           [per_particle_default_values],
        #           {index_in_topology: index_in_force}
        # } 
        self._custom_external_forces = {}
        
        # List of IsoldeMap objects registered with this simulation
        self._maps = []
        # CustomExternalForce handling haptic interactions
        self._tugging_force = None
        # Dict for mapping atom indices in the topology to indices in the tugging force
        self._tug_force_lookup = {}
        
    
    def update_restraints_in_context(self, context):
        self.update_distance_restraints_in_context(context)
        self.update_position_restraints_in_context(context)
        self.update_dihedral_restraints_in_context(context)
        for m in self._maps:
            pf = m.get_potential_function()
            pf.update_context_if_needed(context)
            
    
    ####
    # Positional Restraints
    ####
    
        ##
        # Before simulation starts
        ##
    
    def initialize_position_restraints_force(self, max_force):
        '''Just create the force object.'''
        rf = self._position_restraints_force = TopOutRestraintForce(max_force*10)
        return rf
    
    def add_position_restraints(self, restraints, sim_construct):
        '''Add all atoms in a Position_Restraints object to the force.'''
        for r in restraints:
            self.add_position_restraint(r, sim_construct)
    
    def add_position_restraint(self, restraint, sim_construct):
        '''Add one Position_Restraint object to the force.'''
        r = restraint
        atom = r.atom
        sc = sim_construct
        index = sc.index(atom)
        rf = self._position_restraints_force
        r.sim_handler = self
        r.sim_force_index = rf.addParticle(index, (r.spring_constant * 100, *(r.target/10)))

        ##
        # During simulation
        ##
    
    def change_position_restraint_parameters(self, restraint):
        rf = self._position_restraints_force
        index = self._chimerax_atoms.index(restraint.atom)
        rf.setParticleParameters(restraint.sim_force_index, index, 
            (restraint.spring_constant*100, *(restraint.target/10)))
        rf.update_needed = True
    
    def update_position_restraints_in_context(self, context):
        rf = self._position_restraints_force
        if rf.update_needed:
            rf.updateParametersInContext(context)
        rf.update_needed = False

        ##
        # Cleanup on simulation termination
        ##

    def disconnect_position_restraints_from_sim(self, restraints):
        for r in restraints:
            r.sim_force_index = -1
            r.sim_handler = None
    
    ####
    # Distance Restraints
    ####
    
        ##
        # Before simulation starts
        ##
        
    def initialize_distance_restraints_force(self, max_force):
        tf = self._distance_restraints_force = TopOutBondForce(max_force*10)
        return tf
    
    def add_distance_restraints(self, restraints, sim_construct):
        for r in restraints:
            self.add_distance_restraint(r, sim_construct)
    
    def add_distance_restraint(self, restraint, sim_construct):
        r = restraint
        atoms = r.atoms
        indices = sim_construct.indices(atoms)
        tf = self._distance_restraints_force
        r.sim_handler = self
        r.sim_force_index = tf.addBond(*indices.tolist(), (r.spring_constant*100, r.target_distance/10))
        
        ##
        # During simulation
        ##

    def change_distance_restraint_parameters(self, restraint):
        tf = self._distance_restraints_force
        indices = self._chimerax_atoms.indices(restraint.atoms).tolist()
        tf.setBondParameters(restraint._sim_force_index, *indices, 
            (restraint.spring_constant*100, restraint.target_distance/10))
        tf.update_needed = True
        
    def update_distance_restraints_in_context(self, context):
        tf = self._distance_restraints_force
        if tf.update_needed:
            tf.updateParametersInContext(context)
        tf.update_needed = False
    
        ##
        # Cleanup on simulation termination
        ##
    
    def disconnect_distance_restraints_from_sim(self, restraints):
        for r in restraints:
            r.sim_force_index = -1
            r.sim_handler = None
    
    
    ####
    # AMBER CMAP corrections
    ####
        
        ##
        # Before simulation starts
        ##
    
    def initialize_amber_cmap_force(self):
        self._amber_cmap_force = AmberCMAPForce()
    
    def add_amber_cmap_torsions(self, phi_array, psi_array, sim_construct):
        cf = self._amber_cmap_force
        sc = sim_construct
        phi_has_psi = phi_array.residues.indices(psi_array.residues)
        psi_has_phi = psi_array.residues.indices(phi_array.residues)
        good_phi = phi_array[phi_has_psi[phi_has_psi != -1]]
        good_psi = psi_array[psi_has_phi[psi_has_phi != -1]]
        phi_resnames = good_phi.residues.names
        psi_resnames = good_psi.residues.names
        n = len(good_psi)
        phi_indices = sc.indices(good_phi.atoms).reshape((n,4))
        psi_indices = sc.indices(good_psi.atoms).reshape((n,4))
        
        for pname, fname, pi, fi in zip(phi_resnames, psi_resnames, 
                                        phi_indices, psi_indices):
            assert(pname == fname)
            cf.addTorsion(pname, pi, fi)
        
    
    
    ####
    # Dihedral restraints
    ####
    
        ##
        # Before simulation starts
        ##
    
    def initialize_dihedral_restraint_force(self, default_cutoff):
        
        self._dihedral_restraint_force = FlatBottomTorsionRestraintForce()
        self.default_torsion_cutoff = default_cutoff

    
    def initialize_dihedral_restraint(self, indices, cutoff = None):
        #top = self._topology
        c = (cutoff or self.default_torsion_cutoff)
        if type(c) == Quantity:
            c = c.value_in_unit(defaults.OPENMM_ANGLE_UNIT)
        force = self._dihedral_restraint_force
        index_in_force = force.addTorsion(*indices.tolist(), [0, 0, cos(c)])
        return index_in_force

        ##
        # During simulation
        ##
    
    def update_dihedral_restraints_in_context(self, context):
        rf = self._dihedral_restraint_force
        if rf.update_needed:
            rf.updateParametersInContext(context)
        rf.update_needed = False
    
    def set_dihedral_restraints(self, dihedrals, target, k, degrees = False, cutoffs = None):
        c = (cutoffs or [self.default_torsion_cutoff]*len(dihedrals))
        variable_t = hasattr(target, '__iter__')
        variable_k = hasattr(k, '__iter__')
        variable_c = hasattr(c, '__iter__')
        for i, d in enumerate(dihedrals):
            if variable_t:
                t = target[i]
            else:
                t = target
            if degrees:
                t = radians(t)
            if variable_k:
                thisk = k[i]
            else:
                thisk = k
            if variable_c:
                thisc = c[i]
            else:
                thisc = c
            
            self.update_dihedral_restraint(d.sim_index, target=t, k=thisk, cutoff=thisc)    
            
    def update_dihedral_restraint(self, sim_index, target = None, 
                            k = None, cutoff = None, degrees = False):
        if target is not None and degrees:
            target = radians(target)
        self._dihedral_restraint_force.update_target(sim_index, target=target, k=k, cutoff=cutoff)
        
    def register_custom_external_force(self, name, force, global_params, 
                                        per_particle_params, 
                                        per_particle_default_vals):
        self._custom_external_forces[name] = [force, global_params, per_particle_params, per_particle_default_vals, {}]
    
    def get_custom_external_force_by_name(self, name):
        return self._custom_external_forces[name]
    
    def get_all_custom_external_forces(self):
        return self._custom_external_forces
    
    def set_custom_external_force_particle_params(self, name, index, params):
        fparams = self._custom_external_forces[name]
        force = fparams[0]
        index_lookup = fparams[4]
        index_in_force = index_lookup[index]
        force.setParticleParameters(index_in_force, index, params)
       
    def register_map(self, map_object):
        self._maps.append(map_object)

    #TODO: Define the tugging force as a subclass in custom_forces.py
    def initialize_tugging_force(self):
        from simtk import openmm as mm
        potential_equation = '0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)'
        per_particle_parameters = ['k','x0','y0','z0']
        per_particle_defaults = [0,0,0,0]
        global_parameters = None
        global_defaults = None
        f = self._tugging_force = mm.CustomExternalForce(potential_equation)
        if global_parameters is not None:
            for p, v in zip(global_parameters, g_vals):
                f.addGlobalParameter(g, v)
        if per_particle_parameters is not None:
            for p in per_particle_parameters:
                f.addPerParticleParameter(p)
        self.register_custom_external_force('tug', f, global_parameters,
                per_particle_parameters, per_particle_defaults)

        return f
    
    def _openmm_topology_and_external_forces(self, sim_data, 
                                            tug_hydrogens = False,
                                            hydrogens_feel_maps = False):
        aname   = sim_data['atom names']
        ename   = sim_data['element names']
        rname   = sim_data['residue names']
        rnum    = sim_data['residue numbers']
        cids    = sim_data['chain ids']
        coords  = sim_data['coords']
        bond_i  = sim_data['bonded atom indices']
        
        fixed_indices = sim_data['fixed indices']
        fixed_flags = numpy.zeros(len(aname),numpy.bool)
        fixed_flags[fixed_indices] = True

        from simtk.openmm.app import Topology, Element
        from simtk import unit
        top = self._simulation_topology = Topology()
        cmap = {}
        rmap = {}
        atoms = self._atoms = {}
        for i in range(n):
            cid = cids[i]
            if not cid in cmap:
                cmap[cid] = top.addChain()   # OpenMM chains have no name
            rid = (rname[i], rnum[i], cid)
            rid = (rname[i], rnum[i], cid)
            if not rid in rmap:
                res = rmap[rid] = top.addResidue(rname[i], cmap[cid])
                if rname[i] == 'CYS':
                    ctype = cys_type(r[i])
                    if ctype != 'CYS':
                        templates[res] = ctype

            element = Element.getBySymbol(ename[i])
            atoms[i] = top.addAtom(aname[i], element,rmap[rid])

            if not fixed_flags[i]:
                # Register atoms with forces
                if ename[i] is not 'H' or (ename[i] is 'H' and tug_hydrogens):
                    # All CustomExternalForces
                    for key, ff in self._custom_external_forces.items():
                        f = ff[0]
                        per_particle_param_vals = ff[3]
                        index_map = ff[4]
                        index_map[i] = f.addParticle(i, per_particle_param_vals)
            
                if ename[i] is not 'H' or (ename[i] is 'H' and hydrogens_feel_maps):
                    # All map forces
                    for m in self._maps:
                        self.couple_atom_to_map(i, m)

        for i1, i2 in zip(*bond_i):
            top.addBond(atoms[i1],  atoms[i2])

        return top, templates


        
        
    
    def openmm_topology_and_external_forces(self, sim_construct,
                                        sim_bonds, fixed_flags,
                                        tug_hydrogens = False,
                                        hydrogens_feel_maps = False):        
        '''
        Prepares the openmm topology and binds atoms to existing force fields.
        Since looping over all atoms can take a long time, it's best to do
        topology generation and force binding as a single concerted loop as
        much as we can. This should be possible for all external-type forces
        (tugging and map forces). Forces involving two or more atoms (e.g.
        H-bond or dihedral restraints) will have to be handled in separate
        loops.
        Args:
            sim_construct:
                The atoms to be simulated
            sim_bonds:
                The set of all bonds between simulated atoms
            fixed_flags:
                A boolean array indicating which atoms will be fixed. No
                need to add these to any custom forces.
            tug_hydrogens:
                Do we want to be able to interactively pull on hydrogens?
            hydrogens_feel_maps:
                Do we want the hydrogens to be pulled into the maps?
        '''
        # When we allow missing external bonds, some residues become ambiguous.
        # In particular, a cysteine with a bare sulphur might be part of a 
        # disulphide bond but have the connecting residue missing, or may be
        # a deprotonated (negatively-charged) Cys. In such cases, we need to 
        # explicitly tell OpenMM which templates to use.
        templates = {}
        a = self._chimerax_atoms = sim_construct
        n = len(a)
        r = a.residues
        aname = a.names
        ename = a.element_names
        rname = r.names
        rnum = r.numbers
        cids = r.chain_ids
        from simtk.openmm.app import Topology, Element
        from simtk import unit
        top = self._simulation_topology = Topology()
        cmap = {}
        rmap = {}
        atoms = self._atoms = {}
        for i in range(n):
            cid = cids[i]
            if not cid in cmap:
                cmap[cid] = top.addChain()   # OpenMM chains have no name
            rid = (rname[i], rnum[i], cid)
            rid = (rname[i], rnum[i], cid)
            if not rid in rmap:
                res = rmap[rid] = top.addResidue(rname[i], cmap[cid])
                if rname[i] == 'CYS':
                    ctype = cys_type(r[i])
                    if ctype != 'CYS':
                        templates[res] = ctype

            element = Element.getBySymbol(ename[i])
            atoms[i] = top.addAtom(aname[i], element,rmap[rid])

            if not fixed_flags[i]:
                # Register atoms with forces
                if ename[i] is not 'H' or (ename[i] is 'H' and tug_hydrogens):
                    # All CustomExternalForces
                    for key, ff in self._custom_external_forces.items():
                        f = ff[0]
                        per_particle_param_vals = ff[3]
                        index_map = ff[4]
                        index_map[i] = f.addParticle(i, per_particle_param_vals)
            
                if ename[i] is not 'H' or (ename[i] is 'H' and hydrogens_feel_maps):
                    # All map forces
                    for m in self._maps:
                        self.couple_atom_to_map(i, m)

        
        a1, a2 = sim_bonds.atoms
        for i1, i2 in zip(a.indices(a1), a.indices(a2)):
            if -1 not in [i1, i2]:
                top.addBond(atoms[i1],  atoms[i2])

        pos = a.coords # in Angstrom (convert to nm for OpenMM)
        return top, pos, templates


    def couple_atom_to_map(self, index, map_object):
        '''
        Adds an atom to a map-derived potential field.
        The argument per_atom_coupling must be either a single value, 
        or an array with one value per atom
        '''
        m = map_object
        map_field = m.get_potential_function()
        k = m.get_per_atom_coupling_params()
        if m.per_atom_coupling():
            k = k[index]
        map_field.addBond([index],[k])


    #######################
    # OLD VERSIONS
    #######################


    #def continuous3D_from_volume(self, volume):
        #'''
        #Takes a volumetric map and uses it to generate an OpenMM 
        #Continuous3DFunction. Returns the function.
        #'''
        #import numpy as np
        #vol_data = volume.data
        #mincoor = np.array(vol_data.origin)
        #maxcoor = mincoor + volume.data_origin_and_step()[1]*(np.array(vol_data.size)-1)
        ##Map data is in Angstroms. Need to convert (x,y,z) positions to nanometres
        #mincoor = mincoor/10
        #maxcoor = maxcoor/10
        ## Continuous3DFunction expects the minimum and maximum coordinates as
        ## arguments xmin, xmax, ymin, ...
        #minmax = [val for pair in zip(mincoor, maxcoor) for val in pair]
        #vol_data_1d = np.ravel(vol_data.matrix(), order = 'C').astype(np.double)
        #vol_dimensions = (vol_data.size)
        #print('Volume dimensions: {}; expected number: {}; actual number: {}'\
                #.format(vol_dimensions, np.product(vol_dimensions), len(vol_data_1d)))
        #print('Max: {}, min: {}, nans: {}, infs: {}'.format(
            #vol_data_1d.max(), vol_data_1d.min(), 
            #np.argwhere(np.isnan(vol_data_1d)),
            #np.argwhere(np.isinf(vol_data_1d))))
        #return Continuous3DFunction(*vol_dimensions, vol_data_1d, *minmax)
        
    #def map_potential_force_field(self, c3d_func, global_k):
        #'''
        #Takes a Continuous3DFunction and returns a CustomCompoundBondForce 
        #based on it.
        #Args:
            #c3d_func:
                #A Continuous3DFunction
            #global_k:
                #An overall global spring constant coupling atoms to the 
                #map. This can be further adjusted per atom using 
                #the "individual_k" parameter defined in the 
                #CustomCompoundBondForce energy function.
        #'''
        #from simtk.openmm import CustomCompoundBondForce
        #f = CustomCompoundBondForce(1,'')
        #f.addTabulatedFunction(name = 'map_potential', function = c3d_func)
        #f.addGlobalParameter(name = 'global_k', defaultValue = global_k)
        #f.addPerBondParameter(name = 'individual_k')
        #f.setEnergyFunction('-global_k * individual_k * map_potential(x1,y1,z1)')
        #return f
    
    #######################
    # /OLD VERSIONS
    #######################

    def map_to_force_field(self, imap, atoms, 
                            pad, normalize = True):
        '''
        Takes an IsoldeMap object, masks down the map to cover the atoms
        with the desired padding, and converts it into a LinearInterpMapForce
        object. 
        '''
        vd, r = imap.crop_to_selection(atoms, pad, normalize)
        tf = r.xyz_to_ijk_transform.matrix
        f = self._map_to_force_field(vd, r, tf, imap.get_coupling_constant())
        imap.set_potential_function(f)
    
    def _map_to_force_field(self, volume_data, region, xyz_to_ijk_transform,
                            coupling_constant):
        vd = volume_data
        r = region
        tf = xyz_to_ijk_transform
        f = LinearInterpMapForce(vd, tf, units='angstroms')
        f.set_global_k(coupling_constant)
        return f
    
    
    def continuous3D_from_volume(self, vol_data):
        '''
        Takes a volumetric map and uses it to generate an OpenMM 
        Continuous3DFunction. Returns the function.
        '''
        import numpy as np
        vol_dimensions = vol_data.shape[::-1]
        mincoor = np.array([0,0,0], np.double)
        maxcoor = (np.array(vol_dimensions, np.double) - 1) / 10
        # Continuous3DFunction expects the minimum and maximum coordinates as
        # arguments xmin, xmax, ymin, ...
        minmax = [val for pair in zip(mincoor, maxcoor) for val in pair]
        vol_data_1d = np.ravel(vol_data, order = 'C')
        from simtk.openmm.openmm import Continuous3DFunction    
        return Continuous3DFunction(*vol_dimensions, vol_data_1d, *minmax)
    
    def continuous3D_from_maps(self, master_map_list, keys, atoms, pad, normalize):
        '''
        Combine multiple maps into a single Continuous3DFunction by 
        summing their data with applied weights. 
        '''
        from . import volumetric
        combined_map = None
        for key in keys:
            m = master_map_list[key]
            vd, r = m.crop_to_selection(atoms, pad, normalize)
            weighted_map = vd*m.get_coupling_constant()
            if combined_map is None:
                combined_map = weighted_map
            else:
                combined_map += weighted_map
        c3d = self.continuous3D_from_volume(combined_map)
        f = self.map_potential_force_field(c3d, 1.0, r)
        ret = GenericMapObject(f, False, 1.0)
        return ret
        

    def map_potential_force_field(self, c3d_func, global_k, region):
        '''
        Takes a Continuous3DFunction and returns a CustomCompoundBondForce 
        based on it.
        Args:
            c3d_func:
                A Continuous3DFunction
            global_k:
                An overall global spring constant coupling atoms to the 
                map. This can be further adjusted per atom using 
                the "individual_k" parameter defined in the 
                CustomCompoundBondForce energy function.
            xyz_to_ijk_transform:
                The affine transformation matrix mapping (x,y,z) coordinates
                back to (i,j,k) in the c3d_func array
        '''
        xyz_to_ijk_transform = region.xyz_to_ijk_transform
        from simtk.openmm import CustomCompoundBondForce
        f = CustomCompoundBondForce(1,'')
        f.addTabulatedFunction(name = 'map_potential', function = c3d_func)
        f.addGlobalParameter(name = 'global_k', defaultValue = global_k)
        f.addPerBondParameter(name = 'individual_k')
        tf = xyz_to_ijk_transform.matrix
        tf[:,3] /= 10 # OpenMM in nm, ChimeraX in Angstroms
        i_str = 'x1* {} + y1 * {} + z1 * {} + {}'.format(
            tf[0][0], tf[0][1], tf[0][2], tf[0][3])
        j_str = 'x1* {} + y1 * {} + z1 * {} + {}'.format(
            tf[1][0], tf[1][1], tf[1][2], tf[1][3])
        k_str = 'x1* {} + y1 * {} + z1 * {} + {}'.format(
            tf[2][0], tf[2][1], tf[2][2], tf[2][3])
        
        f.setEnergyFunction('-global_k * individual_k * map_potential({},{},{})'.format(
        i_str, j_str, k_str))
        return f

    
    def update_force_in_context(self, force_name, context):
        force = self._custom_external_forces[force_name][0]
        force.updateParametersInContext(context)


def define_forcefield (forcefield_list):
    from simtk.openmm.app import ForceField
    ff = self_forcefield = ForceField(*[f for f in forcefield_list if f is not None])
    return ff

def initialize_implicit_solvent(system, top):
    '''Add a Generalised Born Implicit Solvent (GBIS) formulation.'''
    # Somewhat annoyingly, OpenMM doesn't store atomic charges in a 
    # nice accessible format. So, we have to pull it back out of the
    # NonbondedForce term.
    for f in system.getForces():
        if isinstance(f, NonbondedForce):
            break
    charges = []
    for i in range(f.getNumParticles()):
        charges.append(f.getParticleParameters(i)[0])
    gbforce = GBSAForce()
    params = GBSAForce.getStandardParameters(top)
    print('charge length: {}, param length: {}, natoms: {}, top_n: {}'.format(
        len(charges), len(params), f.getNumParticles(), top.getNumAtoms()))
    i = 0
    for charge, param in zip(charges, params):
        gbforce.addParticle([charge, *param])
        i+= 1
    print(i)
    gbforce.finalize()
    print('GB Force num particles: '.format(gbforce.getNumParticles()))
    system.addForce(gbforce)


    
def create_openmm_system(top, ff, templatedict, force_implicit=False):
    from simtk.openmm import app
    from simtk import openmm as mm
    from simtk import unit

    params = {
        'nonbondedMethod':      app.CutoffNonPeriodic,
        'nonbondedCutoff':      1.0*unit.nanometers,
        'constraints':          app.HBonds,
        'rigidWater':           True,
        'removeCMMotion':       False,
        'residueTemplates':     templatedict,
        'ignoreExternalBonds':  True,
        }
    
        
    try:
        #~ system = ff.createSystem(top,
                                #~ nonbondedMethod = app.CutoffNonPeriodic,
                                #~ nonbondedCutoff = 1.0*unit.nanometers,
                                #~ constraints = app.HBonds,
                                #~ rigidWater = True,
                                #~ removeCMMotion = False,
                                #~ residueTemplates = templatedict,
                                #~ ignoreExternalBonds = True)
        system = ff.createSystem(top, **params)
        if force_implicit:
            print('Setting up implicit solvent...')
            initialize_implicit_solvent(system, top)
    except ValueError as e:
        raise Exception('Missing atoms or parameterisation needed by force field.\n' +
                              'All heavy atoms and hydrogens with standard names are required.\n' +
                              str(e))
    return system


def integrator(i_type, temperature, friction, tolerance, timestep):
    from simtk import openmm as mm
    if i_type == mm.VariableLangevinIntegrator:
        return i_type(temperature, friction, tolerance)
    elif i_type == mm.LangevinIntegrator:
        return i_type(temperature, friction, timestep)
    raise TypeError('Unrecognised integrator!')
    
def platform(name):
    from simtk.openmm import Platform
    return Platform.getPlatformByName(name)


def create_sim(topology, system, integrator, platform):
    from simtk.openmm.app import Simulation
    return Simulation(topology, system, integrator, platform)
    

    

