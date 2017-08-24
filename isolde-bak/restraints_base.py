import numpy
from simtk.unit import kilojoule_per_mole, angstrom, nanometer
from chimerax.core.atomic import Atom, Atoms, concatenate

HIDE_ISOLDE = 0x2

class Distance_Restraint:
    '''
    Base class for distance restraints between atoms. Defines the scheme
    for fast look-up of restraints and control within the context of 
    simulations.
    '''
    def __init__(self, atoms, target_distance, spring_constant):
        '''
        Defines a distance restraint between a pair of atoms.
        Args:
            atoms:
                A ChimeraX Atoms object containing exactly two atoms
            target_distance:
                The desired distance between the two atoms in Angstroms
            spring_constant:
                The spring constant for the harmonic restraining potential,
                in kJ/mol/A^3
        '''
        self._atoms = atoms
        self._target_distance = target_distance
        self._spring_constant = spring_constant
        # Placeholder for the index returned when this restraint is added
        # to an OpenMM Force object.
        self._sim_force_index = -1
        # The SimHandler object handling the current simulation
        self._sim_handler = None
    
    @property
    def atoms(self):
        return self._atoms
    
    @property
    def target_distance(self):
        '''
        The target separation of this atom pair in Angstroms.
        '''
        return self._target_distance
    
    @target_distance.setter
    def target_distance(self, distance):
        self._target_distance = distance
        if self._sim_handler is not None:
            self._sim_handler.change_distance_restraint_parameters(self)
    
    @property
    def distance(self):
        '''The current distance between the restrained atoms in Angstroms.'''
        coords = self.atoms.coords
        return numpy.linalg.norm(coords[1]-coords[0])
    
    @property
    def spring_constant(self):
        '''
        Set the spring constant for this restraint, in kJ/mol/A^3
         '''
        return self._spring_constant
    
    @spring_constant.setter
    def spring_constant(self, k):
        self._spring_constant = k
        if self._sim_handler is not None:
            self._sim_handler.change_distance_restraint_parameters(self)
     
    @property
    def sim_handler(self):
        '''The SimHandler object handling the current simulation'''
        return self._sim_handler
    
    @sim_handler.setter
    def sim_handler(self, handler):
        self._sim_handler = handler
    
    @property
    def sim_force_index(self):
        '''The index of this restraint in the OpenMM Force object.'''
        return self._sim_force_index
    
    @sim_force_index.setter
    def sim_force_index(self, index):
        if self.sim_handler is None:
            self._sim_force_index = -1
            raise RuntimeError('This restraint is not associated with a Force!')
        self._sim_force_index = index

class Distance_Restraints:
    '''
    Holds an array of Distance_Restraint objects 
    '''
    def __init__(self, restraints_list):
        '''
        Initialise from an array of Distance_Restraint objects
        '''
        self._restraints = restraints_list
        # ChimeraX Atoms object holding all atoms in this object. The atoms
        # corresponding to restraint i will be found at [2*i: 2*i+2]
        self._atoms = None
        # Current target distances
        self._targets = None
        # Current actual distances
        self._distances = None
        
    
    def __len__(self):
        return len(self._restraints)
    
    def __bool__(self):
        return len(self) > 0
    
    def __iter__(self):
        return iter(self._restraints)
    
    def __getitem__(self, i):
        if not len(self):
            return None
        if isinstance(i,(int, numpy.integer)):
            return self._restraints[i]
        if isinstance(i,(slice)):
            return self.__class__(self._restraints[i])
        if isinstance(i, numpy.ndarray):
            return self.__class__([self._restraints[j] for j in i])
        if isinstance(i, Atom):
            # Find and return all restraints this atom is involved in
            atom_indices = self.atoms.indices(Atoms([i]))
            atom_indices = atom_indices[numpy.where(atom_indices != -1)]
            restraint_indices = atom_indices//2
            return self[restraint_indices]
        raise IndexError('Only integer indices allowed for {}, got {}'
            .format(self.__class__.__name__, str(type(i))))
    
    def index(self, restraint):
        try:
            i = self._restraints.index(restraint)
        except ValueError:
            return -1
        return i
    
    def append(self, r):
        if isinstance(r, Distance_Restraint):
            if self.atoms is None:
                self._restraints = [r]
            else:
                self.restraints.append(r)
                self._atoms = concatenate((self.atoms, r.atoms))
        elif isinstance(r, Distance_Restraints):
            if self.atoms is None:
                self._restraints = r
            else:
                self.restraints.extend(r)
                self._atoms = concatenate((self.atoms, r.atoms))
        else:
            raise TypeError('Can only append a single Distance_Restraint or \
                             a Distance_Restraints object.')
    
    def in_selection(self, sel):
        '''
        Returns a Distance_Restraints object encompassing all distance
        restraints where both atoms are in the given selection.
        Args:
            sel:
                An Atoms object.
        '''
        atom_indices = sel.indices(self.atoms).reshape((len(self.atoms)//2, 2))
        return self[numpy.argwhere(numpy.all(atom_indices != -1, axis = 1)).ravel()]
        
    
    @property
    def atoms(self):
        if not len(self):
            return None
        if self._atoms is None:
            atoms = numpy.empty([len(self), 2], dtype='object')
            for i, r in enumerate(self):
                atoms[i] = r.atoms
            self._atoms = Atoms(numpy.ravel(atoms))
        return self._atoms
    
    @property
    def restraints(self):
        return self._restraints
    
    @property
    def coords(self):
        if self.atoms is None:
            return None
        return self.atoms.coords
    
    @property
    def targets(self):
        return numpy.array([r.target_distance for r in self])
    
    @targets.setter
    def targets(self, targets):
        if len(targets) != len(self):
            raise IndexError('Target array length must equal the number of restraints!')
        for r, t in zip(self, targets):
            r.target_distance = t

    @property
    def all_distances(self):
        '''Returns the current distances between restrained atom pairs, in Angstroms.'''
        coords = numpy.reshape(self.coords, [len(self)//2, 2,3])
        return numpy.linalg.norm(coords[:,1]-coords[:,0],axis=1)
    
    @property
    def spring_constants(self):
        '''Returns the spring constants for all restraints, in kJ/mol/A^3.'''
        return numpy.array([r.spring_constant for r in self.restraints])
    
    
    
    
    
class Position_Restraint:
    '''
    Restrains one atom to an (x,y,z) position in space via a harmonic spring
    with a maximum limit on the applied force.
    '''
    def __init__(self, atom, pseudobond_group, triggers = None):
        '''
        Just prepare the internal data structure and initialise spring 
        constant to zero and target to (0,0,0).
        Args:
            atom:
                The atom managed by this object
            pseudobond_group:
                The PseudobondGroup object that will manage visualisation
            triggers:
                The ISOLDE triggerset
        '''
        self._atom = atom
        self._target = numpy.array([0,0,0], float)
        # Optionally, we can provide an atom representing the target and
        # a pseudobond between the master and target atoms, to provide a
        # visual representation of the restraint.
        self._target_atom = None
        self._pbg = pseudobond_group 
        self._pseudobond = None
        self._spring_constant = 0
        self._sim_force_index = -1
        self._sim_handler = None
        self._triggers = triggers
    
    @property
    def atom(self):
        return self._atom
    
    @property
    def target(self):
        return self._target
    
    @target.setter
    def target(self, xyz):
        self._target = xyz
        if self._target_atom is not None:
            self._target_atom.coord = xyz
        if self._sim_handler is not None:
            self._sim_handler.change_position_restraint_parameters(self)
        
    @property
    def target_indicator(self):
        '''
        Visual representation of the target position. Read-only.
        '''
        return self._target_atom
    
    @target_indicator.setter
    def target_indicator(self, atom):
        self._target_atom = atom
        
    @property
    def target_bond(self):
        '''
        Dashed bond linking the restrained atom to its target. Read-only.
        '''
        return self._pseudobond
    
    @property
    def spring_constant(self):
        return self._spring_constant
    
    @spring_constant.setter
    def spring_constant(self, k):
        self._spring_constant = k
        ta = self._target_atom
        pb = self._pseudobond
        if k == 0:
            if ta is not None:
                ta.display = False
                ta.hide |= HIDE_ISOLDE
            if pb is not None:
                pb.delete()
                self._pseudobond = None
            if self._triggers is not None:
                self._triggers.activate_trigger('position restraint removed', self)
        else:
            if ta is not None:
                ta.display = True
                ta.hide ^= HIDE_ISOLDE
                if pb is None:
                    pb = self._pseudobond = self._pbg.new_pseudobond(self._atom, ta)
                pb.display = True
            if self._triggers is not None:
                self._triggers.activate_trigger('position restraint added', self)
        if self._sim_handler is not None:
            self._sim_handler.change_position_restraint_parameters(self)
    
    @property
    def sim_handler(self):
        '''The SimHandler object handling the current simulation'''
        return self._sim_handler
    
    @sim_handler.setter
    def sim_handler(self, handler):
        self._sim_handler = handler
    
    @property
    def sim_force_index(self):
        '''The index of this restraint in the OpenMM Force object.'''
        return self._sim_force_index
    
    @sim_force_index.setter
    def sim_force_index(self, index):
        if self.sim_handler is None:
            self._sim_force_index = -1
            raise RuntimeError('This restraint is not associated with a Force!')
        self._sim_force_index = index

class Position_Restraints:
    '''Holds an array of Position_Restraint objects.'''
    def __init__(self, session, master_model, restraints_list):
        self.session = session
        self.master_model = master_model
        self._restraints = restraints_list
        # ChimeraX Atoms object holding all atoms in this object.
        self._atoms = None
        # Current target distances
        self._targets = None
        # Current actual distances
        self._distances = None
        # Atoms object for visualising the target positions
        self._target_indicators = None
     
        
    @property
    def atoms(self):
        if not len(self):
            return None
        if self._atoms is None:
            atoms = numpy.empty(len(self), dtype='object')
            for i, r in enumerate(self):
                atoms[i] = r.atom
            self._atoms = Atoms(atoms)
        return self._atoms
    
    @property
    def restraints(self):
        return self._restraints
    
    @property
    def coords(self):
        if self.atoms is None:
            return None
        return self.atoms.coords
    
    @property
    def targets(self):
        if self.target_indicators is not None:
            return self.target_indicators.coords
        return numpy.array([r.target for r in self])
    
    @targets.setter
    def targets(self, targets):
        if len(targets) != len(self):
            raise IndexError('Target array length must equal the number of restraints!')
        for r, t in zip(self, targets):
            r.target = t
    
    @property
    def target_indicators(self):
        
        if self._target_indicators is None:
            targets = [r.target_indicator for r in self if r.target_indicator is not None]
            if not len(targets):
                self._target_indicators = None
            self._target_indicators = Atoms(targets)
        return self._target_indicators
    
    @property
    def all_distances(self):
        '''Returns the current distances between atoms and their targets, in Angstroms.'''
        return numpy.linalg.norm(self.coords - self.targets,axis=1)
    
    @property
    def spring_constants(self):
        '''Returns the spring constants for all restraints, in kJ/mol/A^3.'''
        return numpy.array([r.spring_constant for r in self.restraints])
    
    @property
    def restrained_bond_vectors(self):
        '''
        Returns a (nx1 numpy.int, nx3 numpy.float) tuple where the first
        array holds the indices of all currently active restraints, and
        the second array holds the vectors connecting the atoms to their
        target positions.
        '''
        kraw = self.spring_constants
        indices = numpy.argwhere(kraw != 0).ravel()
        vectors = (self.targets - self.coords)[indices]
        return (indices, vectors)
        
    def release(self):
        '''
        Release all restraints (last target values will be unchanged, but
        spring constants will be set to zero).
        '''
        for r in self:
            r.spring_constant = 0
    
    
    def __len__(self):
        return len(self._restraints)
    
    def __bool__(self):
        return len(self) > 0
    
    def __iter__(self):
        return iter(self._restraints)
    
    def __getitem__(self, i):
        if not len(self):
            return None
        if isinstance(i,(int, numpy.integer)):
            return self._restraints[i]
        if isinstance(i,(slice)):
            return self.__class__(self._restraints[i])
        if isinstance(i, numpy.ndarray):
            if len(i):
                return self.__class__(self.session, self.master_model, [self._restraints[j] for j in i])
            raise IndexError('No indices in array!')
        if isinstance(i, Atom):
            return self[self.atoms.index(i)]
        raise IndexError('Only integer indices allowed for {}, got {}'
            .format(self.__class__.__name__, str(type(i))))
    
    def index(self, restraint):
        try:
            i = self._restraints.index(restraint)
        except ValueError:
            return -1
        return i
    
    def append(self, r):
        if isinstance(r, Position_Restraint):
            if self.atoms is None:
                self._restraints = [r]
            else:
                self.restraints.append(r)
                self._atoms = concatenate((self.atoms, r.atoms))
        elif isinstance(r, Position_Restraints):
            if self.atoms is None:
                self._restraints = r
            else:
                self.restraints.extend(r)
                self._atoms = concatenate((self.atoms, r.atoms))
        else:
            raise TypeError('Can only append a single Distance_Restraint or \
                             a Distance_Restraints object.')
    
    def in_selection(self, sel):
        '''
        Returns a Position_Restraints object encompassing all restraints
        corresponding to atoms in the given selection.
        Args:
            sel:
                An Atoms object.
        '''
        atom_indices = sel.indices(self.atoms)
        return self[numpy.argwhere(atom_indices != -1).ravel()]
    
    
class Restraint_Bond_Drawings:
    '''
    Handles visual representations of positional and distance restraints.
    Yet to be implemented.
    '''
    pass


        
def restraint_bond_geometry():
    '''
    Create a simple prototype dashed bond of length 1.
    '''
    from chimerax.core.surface.shapes import cylinder_geometry
    from chimerax.core.geometry import translation
    from copy import copy
    v0, n0, t0 = cylinder_geometry(0.075, 0.15)
    v = copy(v0)
    n = copy(n0)
    t = copy(t0)
    for i in range(1, 5):
        tz = translation((0,0,0.2*i))
        nv = len(v)
        v = numpy.concatenate((v, tz.moved(v0)))
        n = numpy.concatenate((n, n0))
        t = numpy.concatenate((t, t0 + nv))
    
    return (v, n, t)
     


