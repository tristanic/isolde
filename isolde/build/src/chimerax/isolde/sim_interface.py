# vim: set expandtab shiftwidth=4 softtabstop=4:

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
            rmap[rid] = top.addResidue(rname[i], cmap[cid])
        element = Element.getBySymbol(ename[i])
        atoms[i] = top.addAtom(aname[i], element,rmap[rid])

    a1, a2 = b.atoms
    for i1, i2 in zip(a.indices(a1), a.indices(a2)):
        if -1 not in [i1, i2]:
            top.addBond(atoms[i1],  atoms[i2])
            
    return top

class available_forcefields():
    # Main force field files
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

class SimHandler():
    
    def __init__(self, session):
        self.session = session
        # Forcefield used in this simulation
        self._forcefield = None
        # Overall simulation topology
        self._topology = None
        # Atoms in simulation topology
        self._atoms = None
        
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
        

        
    
    def initialize_dihedral_restraint(self, dihedral, indices):
        #top = self._topology
        force = self._dihedral_restraint_force
        index_in_force = force.addTorsion(*indices.tolist(), 1, 0, 0)
        return index_in_force
        
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
    
    
    # Prepares the openmm topology and binds atoms to existing force fields.
    # Since looping over all atoms can take a long time, it's best to do
    # topology generation and force binding as a single concerted loop as
    # much as we can. This should be possible for all external-type forces
    # (tugging and map forces). Forces involving two or more atoms (e.g.
    # H-bond or dihedral restraints) will have to be handled in separate
    # loops.
    def openmm_topology_and_external_forces(self, sim_construct,
                                        sim_bonds,
                                        fix_shell_backbones = False,
                                        tug_hydrogens = False,
                                        hydrogens_feel_maps = False,
                                        logging = False,
                                        log = None):        
        a = sim_construct
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
                rmap[rid] = top.addResidue(rname[i], cmap[cid])
            element = Element.getBySymbol(ename[i])
            atoms[i] = top.addAtom(aname[i], element,rmap[rid])
 
            # Register atoms with forces
            if ename is not 'H' or (ename is 'H' and tug_hydrogens):
                # All CustomExternalForces
                for key, ff in self._custom_external_forces.items():
                    f = ff[0]
                    per_particle_param_vals = ff[3]
                    index_map = ff[4]
                    index_map[i] = f.addParticle(i, per_particle_param_vals)
        
            if ename is not 'H' or (ename is 'H' and hydrogens_feel_maps):
                # All map forces
                for m in self._maps:
                    self.couple_atom_to_map(i, m)

        
        a1, a2 = sim_bonds.atoms
        for i1, i2 in zip(a.indices(a1), a.indices(a2)):
            if -1 not in [i1, i2]:
                top.addBond(atoms[i1],  atoms[i2])
        
        from simtk.openmm import Vec3
        pos = a.coords # in Angstrom (convert to nm for OpenMM)
        return top, pos


    # Take the atoms in a topology, and add them to a map-derived potential field.
    # per_atom_coupling must be either a single value, or an array with one value
    # per atom
    def couple_atom_to_map(self, index, map_object):
        m = map_object
        if not m.per_atom_coupling():
            global_coupling = True
            k = m.get_per_atom_coupling_params()
        else:
            global_coupling = False
            per_atom_k = m.get_per_atom_coupling_params()
        # Find the global coupling constant parameter in the Force and set its new value
        map_field = m.get_potential_function()
        if not global_coupling:
            k = per_atom_k[index]
        map_field.addBond([index],[k])





    # Takes a volumetric map and uses it to generate an OpenMM Continuous3DFunction.
    # Returns the function.
    def continuous3D_from_volume(self, volume):
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
        vol_data_1d = np.ravel(vol_data.matrix(), order = 'C')
        vol_dimensions = (vol_data.size)
        from simtk.openmm.openmm import Continuous3DFunction    
        return Continuous3DFunction(*vol_dimensions, vol_data_1d, *minmax)
        
    # Takes a Continuous3DFunction and returns a CustomCompoundBondForce based on it
    def map_potential_force_field(self, c3d_func, global_k):
        from simtk.openmm import CustomCompoundBondForce
        f = CustomCompoundBondForce(1,'')
        f.addTabulatedFunction(name = 'map_potential', function = c3d_func)
        f.addGlobalParameter(name = 'global_k', defaultValue = global_k)
        f.addPerBondParameter(name = 'individual_k')
        f.setEnergyFunction('-global_k * individual_k * map_potential(x1,y1,z1)')
        return f
    
    def update_force_in_context(self, force_name, context):
        force = self._custom_external_forces[force_name][0]
        force.updateParametersInContext(context)


def define_forcefield (forcefield_list):
    from simtk.openmm.app import ForceField
    ff = self_forcefield = ForceField(*forcefield_list)
    return ff
    
def create_openmm_system(top, ff):
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
                                ignoreMissingExternalBonds = True)
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
    

    

