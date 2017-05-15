import numpy
from simtk.unit import kilojoule_per_mole, angstrom, nanometer)
from chimerax.core.atomic import Atom, Atoms, concatenate

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
        self.atoms = atoms
        self._target_distance = target_distance
        self._spring_constant = spring_constant
        # Placeholder for the index returned when this restraint is added
        # to an OpenMM Force object.
        self._sim_force_index = None
        # The actual OpenMM Force object this restraint is associated with
        self._sim_force = None
    
    @property
    def target_distance(self):
        '''
        The target separation of this atom pair in Angstroms.
        '''
        return self._target_distance
    
    @target_distance.setter
    def target_distance(self, distance):
        self._target_distance = distance
        if self._sim_force is not None:
            # OpenMM distances are in nm
            self._sim_force.setBondParameters(self._sim_force_index, 
                (self.spring_constant/1000, self.target_distance/10))
            self._sim_force.update_needed = True
    
    @property
    def spring_constant(self):
        '''
        Set the spring constant for this restraint, in kJ/mol/A^3
         '''
        return self._spring_constant
    
    @spring_constant.setter
    def spring_constant(self, k):
        self._spring_constant = k
        if self._sim_force is not None:
            # OpenMM distances are in nm
            self._sim_force.setBondParameters(self._sim_force_index, 
                (self.spring_constant/1000, self.target_distance/10))
            self._sim_force.update_needed = True

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
        if isinstance(i,(int, numpy.integer)):
            return self._restraints[i]
        if isinstance(i,(slice)):
            return Distance_Restraints(self._restraints[i])
        if isinstance(i, numpy.ndarray):
            return Distance_Restraints([self._restraints[j] for j in i])
        if isinstance(i, Atom):
            # Find and return all restraints this atom is involved in
            atom_indices = self.atoms.indices(Atoms([i]))
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
    
    @targets.setter:
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
    
    
    
