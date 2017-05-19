# vim: set expandtab shiftwidth=4 softtabstop=4:

import numpy
from math import pi
from simtk.openmm.openmm import CustomBondForce, CustomExternalForce
from chimerax.core.atomic import concatenate


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
        if a.name == 'SG':
            return 'CYS'



class available_forcefields():
    '''Main force field files'''
    main_files = [
        'amber99sbildn.xml',
        'amber99sbnmr.xml',
        'amber10.xml'
        ]
    
    main_file_descriptions = [
        'AMBER99 with improved backbone & sidechain torsions',
        'AMBER99 with modifications to fit NMR data',
        'AMBER10'
        ]
    
    # Implicit solvent force field files. The main file and the implicit
    # solvent file must match.
    implicit_solvent_files = [
        'amber99_obc.xml',
        'amber99_obc.xml',
        'amber10_obc.xml'
        ]
    
    # Explicit water models
    explicit_water_files = [
        'tip3pfb.xml',
        'tip4pfb.xml',
        'tip3p.xml'
        ]
    
    explicit_water_descriptions = [
        'TIP3P-FB (DOI: 10.1021/jz500737m)',
        'TIP4P-FB (DOI: 10.1021/jz500737m)',
        'Original TIP3P water (not recommended)'
        ]

def get_available_platforms():
        from simtk.openmm import Platform
        platform_names = []
        for i in range(Platform.getNumPlatforms()):
            p = Platform.getPlatform(i)
            name = p.getName()
            platform_names.append(name)
        return platform_names

class TopOutBondForce(CustomBondForce):
    '''
    Wraps an OpenMM CustomBondForce defined as a standard harmonic potential
    (0.5 * k * (r - r0)^2) with a user-defined fixed maximum cutoff on the
    applied force. This is meant for steering the simulation into new 
    conformations where the starting distance may be far from the target
    bond length, leading to catastrophically large forces with a standard
    harmonic potential.
    '''
    def __init__(self, max_force):
        super().__init__('min(0.5*k*(r-r0)^2, max_force*abs(r-r0))')
        self._max_force = max_force
        self.k_index = self.addPerBondParameter('k')
        self.r0_index = self.addPerBondParameter('r0')
        self.max_force_index = self.addGlobalParameter('max_force', self.max_force)
        self.update_needed = False
    
    @property
    def max_force(self):
        '''Maximum force applied to any given atom, in kJ/mol/nm.'''
        return self._max_force
    
    @max_force.setter
    def max_force(self, force):
        self.setGlobalParameterDefaultValue(self.max_force_index, force)
        self._max_force = force
        self.update_needed = True
    

class TopOutRestraintForce(CustomExternalForce):
    '''
    Wraps an OpenMM CustomExternalForce to restrain atoms to defined positions
    via a standard harmonic potential (0.5 * k * r^2) with a user-defined
    fixed maximum cutoff on the applied force. This is meant for steering 
    the simulation into new conformations where the starting positions
    may be far from the target positions, leading to catastrophically 
    large forces with a standard harmonic potential.
    '''
    def __init__(self, max_force):
        super().__init__('min(0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2), max_force *(abs(x-x0)+abs(y-y0)+abs(z-z0)))')
        self._max_force = max_force
        per_particle_parameters = ['k','x0','y0','z0']
        for p in per_particle_parameters:
            self.addPerParticleParameter(p)
        self.addGlobalParameter('max_force', max_force)
        self.update_needed = False
    
    @property
    def max_force(self):
        '''Maximum force applied to any given atom, in kJ/mol/nm.'''
        return self._max_force
    
    @max_force.setter
    def max_force(self, force):
        self.setGlobalParameterDefaultValue(0, force)
        self._max_force = force
        self.update_needed = True
    
    
    
    
        



class SimHandler():
    
    def __init__(self, session):
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
        
        
        from simtk.openmm.openmm import PeriodicTorsionForce
        
        self._dihedral_restraint_force = PeriodicTorsionForce()
    
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
    # Dihedral restraints
    ####
    
        ##
        # Before simulation starts
        ##
    
    def initialize_dihedral_restraint(self, dihedral, indices):
        #top = self._topology
        force = self._dihedral_restraint_force
        index_in_force = force.addTorsion(*indices.tolist(), 1, 0, 0)
        return index_in_force

        ##
        # During simulation
        ##
    
    def set_dihedral_restraints(self, context, sim_construct, dihedrals, target, k, degrees = False):
        indices = numpy.reshape(sim_construct.indices(dihedrals.atoms), [len(dihedrals),4])
        variable_t = hasattr(target, '__iter__')
        variable_k = hasattr(k, '__iter__')
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
                
            self.set_dihedral_restraint(context, sim_construct, d, indices[i], t, thisk)
        
    def set_dihedral_restraint(self, context, sim_construct, dihedral, indices, target, k, degrees=False):
        from math import pi, radians
        force = self._dihedral_restraint_force
        atoms = dihedral.atoms
        if degrees:
            target = radians(target)
        target += pi
        force.setTorsionParameters(dihedral.sim_index, *indices.tolist(), 1, target, k)
        if context is not None:
            force.updateParametersInContext(context)
    
    
    
        
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

    
    def initialize_tugging_force(self, potential_equation, g_params, g_vals, pp_params):
        from simtk import openmm as mm
        f = self._tugging_force = mm.CustomExternalForce(potential_equation)
        if g_params is not None:
            for p, v in zip(g_params, g_vals):
                f.addGlobalParameter(g, v)
        if pp_params is not None:
            for p in pp_params:
                f.addPerParticleParameter(p)
        return f
    
    
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


    def continuous3D_from_volume(self, volume):
        '''
        Takes a volumetric map and uses it to generate an OpenMM 
        Continuous3DFunction. Returns the function.
        '''
        import numpy as np
        vol_data = volume.data
        mincoor = np.array(vol_data.origin)
        maxcoor = mincoor + volume.data_origin_and_step()[1]*(np.array(vol_data.size)-1)
        #Map data is in Angstroms. Need to convert (x,y,z) positions to nanometres
        mincoor = mincoor/10
        maxcoor = maxcoor/10
        # Continuous3DFunction expects the minimum and maximum coordinates as
        # arguments xmin, xmax, ymin, ...
        minmax = [val for pair in zip(mincoor, maxcoor) for val in pair]
        vol_data_1d = np.ravel(vol_data.matrix(), order = 'C').astype(np.double)
        vol_dimensions = (vol_data.size)
        print('Volume dimensions: {}; expected number: {}; actual number: {}'\
                .format(vol_dimensions, np.product(vol_dimensions), len(vol_data_1d)))
        print('Max: {}, min: {}, nans: {}, infs: {}'.format(
            vol_data_1d.max(), vol_data_1d.min(), 
            np.argwhere(np.isnan(vol_data_1d)),
            np.argwhere(np.isinf(vol_data_1d))))
        from simtk.openmm.openmm import Continuous3DFunction    
        return Continuous3DFunction(*vol_dimensions, vol_data_1d, *minmax)
        
    def map_potential_force_field(self, c3d_func, global_k):
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
        '''
        from simtk.openmm import CustomCompoundBondForce
        f = CustomCompoundBondForce(1,'')
        f.addTabulatedFunction(name = 'map_potential', function = c3d_func)
        f.addGlobalParameter(name = 'global_k', defaultValue = global_k)
        f.addPerBondParameter(name = 'individual_k')
        f.setEnergyFunction('-global_k * individual_k * map_potential(x1,y1,z1)')
        return f
    
    #######################
    # /OLD VERSIONS
    #######################

    #def continuous3D_from_volume(self, volume):
        #'''
        #Takes a volumetric map and uses it to generate an OpenMM 
        #Continuous3DFunction. Returns the function.
        #'''
        #import numpy as np
        #vol_data = volume.data
        #vol_dimensions = vol_data.size
        #mincoor = np.array([0,0,0], np.double)
        #maxcoor = (np.array(vol_dimensions, np.double) - 1) / 10
        ## Continuous3DFunction expects the minimum and maximum coordinates as
        ## arguments xmin, xmax, ymin, ...
        #minmax = [val for pair in zip(mincoor, maxcoor) for val in pair]
        #vol_data_1d = np.ravel(vol_data.matrix(), order = 'C')
        #from simtk.openmm.openmm import Continuous3DFunction    
        #return Continuous3DFunction(*vol_dimensions, vol_data_1d, *minmax)


    #def map_potential_force_field(self, c3d_func, global_k, xyz_to_ijk_transform):
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
            #xyz_to_ijk_transform:
                #The affine transformation matrix mapping (x,y,z) coordinates
                #back to (i,j,k) in the c3d_func array
        #'''
        #from simtk.openmm import CustomCompoundBondForce
        #f = CustomCompoundBondForce(1,'')
        #f.addTabulatedFunction(name = 'map_potential', function = c3d_func)
        #f.addGlobalParameter(name = 'global_k', defaultValue = global_k)
        #f.addPerBondParameter(name = 'individual_k')
        #tf = xyz_to_ijk_transform
        ##tf [0:3, 0:3] *= 10 # OpenMM in nm, ChimeraX in Angstroms
        #tf[:,3] /= 10
        #i_str = 'x1* {} + y1 * {} + z1 * {} + {}'.format(
            #tf[0][0], tf[0][1], tf[0][2], tf[0][3])
        #j_str = 'x1* {} + y1 * {} + z1 * {} + {}'.format(
            #tf[1][0], tf[1][1], tf[1][2], tf[1][3])
        #k_str = 'x1* {} + y1 * {} + z1 * {} + {}'.format(
            #tf[2][0], tf[2][1], tf[2][2], tf[2][3])
        
        #f.setEnergyFunction('-global_k * individual_k * map_potential({},{},{})'.format(
        #i_str, j_str, k_str))
        #return f

    
    def update_force_in_context(self, force_name, context):
        force = self._custom_external_forces[force_name][0]
        force.updateParametersInContext(context)


def define_forcefield (forcefield_list):
    from simtk.openmm.app import ForceField
    ff = self_forcefield = ForceField(*forcefield_list)
    return ff
    
def create_openmm_system(top, ff, templatedict):
    from simtk.openmm import app
    from simtk import openmm as mm
    from simtk import unit
    
    try:
        system = ff.createSystem(top,
                                nonbondedMethod = app.CutoffNonPeriodic,
                                nonbondedCutoff = 1.0*unit.nanometers,
                                constraints = app.HBonds,
                                rigidWater = True,
                                removeCMMotion = False,
                                residueTemplates = templatedict,
                                ignoreExternalBonds = True)
    except ValueError as e:
        raise Exception('Missing atoms or parameterisation needed by force field.\n' +
                              'All heavy atoms and hydrogens with standard names are required.\n' +
                              str(e))
    return system


def integrator(i_type, temperature, friction, tolerance, timestep):
    from simtk import openmm as mm
    if i_type == 'variable':
        integrator = mm.VariableLangevinIntegrator(temperature, friction, tolerance)
    elif i_type == 'fixed':
        integrator = mm.LangevinIntegrator(temperature, friction, timestep)
    return integrator
    
def platform(name):
    from simtk.openmm import Platform
    return Platform.getPlatformByName(name)


def create_sim(topology, system, integrator, platform):
    from simtk.openmm.app import Simulation
    return Simulation(topology, system, integrator, platform)
    

    

