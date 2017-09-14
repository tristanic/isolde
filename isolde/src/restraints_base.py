# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


import numpy
from simtk.unit import Quantity
from chimerax.core.atomic import Atom, Atoms, concatenate
from .constants import defaults
CHIMERAX_SPRING_UNIT = defaults.CHIMERAX_SPRING_UNIT
CHIMERAX_LENGTH_UNIT = defaults.CHIMERAX_LENGTH_UNIT

HIDE_ISOLDE = 0x2

class Distance_Restraint:
    '''
    Base class for distance restraints between atoms. Defines the scheme
    for fast look-up of restraints and control within the context of
    simulations.
    '''
    PB_COLOR = numpy.array([138, 43, 226, 255], numpy.uint8) # violet
    def __init__(self, atoms, name, target_distance, spring_constant, 
                 pseudobond_group = None):
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
        self._name = name
        self._atoms = atoms
        self._target_distance = target_distance
        self._spring_constant = spring_constant
        self._pbg = pseudobond_group
        self._pseudobond = None


    @property
    def name(self):
        return self._name

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
        if isinstance(distance, Quantity):
            distance = distance.value_in_unit(CHIMERAX_LENGTH_UNIT)
        self._target_distance = distance


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
        if isinstance(k, Quantity):
            k = k.value_in_unit(CHIMERAX_SPRING_UNIT)
        self._spring_constant = k
        pb = self._pseudobond
        if k == 0:
            if pb is not None:
                pb.delete()
                self._pseudobond = None
        else:
                if pb is None:
                    pb = self._pseudobond = self._pbg.new_pseudobond(*self.atoms)
                    pb.color = self.PB_COLOR
                pb.display = True




class Distance_Restraints:
    '''
    Holds an array of Distance_Restraint objects
    '''
    def __init__(self, session, restraints_list, name = None):
        '''
        Initialise from an array of Distance_Restraint objects
        '''
        self.session = session
        self._name = name
        self._restraints = restraints_list
        # ChimeraX Atoms object holding all atoms in this object. The atoms
        # corresponding to restraint i will be found at [2*i: 2*i+2]
        self._atoms = None
        # Current target distances
        self._targets = None
        # Current actual distances
        self._distances = None


    @property
    def name(self):
        return self._name

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
            return self.__class__(self.session, self._restraints[i], name = self.name)
        if isinstance(i, numpy.ndarray):
            return self.__class__(self.session, [self._restraints[j] for j in i], name = self.name)
        if isinstance(i, Atom):
            # Find and return all restraints this atom is involved in
            return self[Atoms([i,])]
        if isinstance(i, Atoms):
            atom_indices = i.indices(self.atoms)
            atom_indices = numpy.where(atom_indices != -1)[0]
            restraint_indices = numpy.unique(atom_indices//2)
            return self[restraint_indices]
            
        raise IndexError('Only integer indices allowed for {}, got {}'
            .format(self.__class__.__name__, str(type(i))))

    def index(self, restraint):
        try:
            i = self._restraints.index(restraint)
        except ValueError:
            return -1
        return i

    def indices(self, restraints):
        return numpy.array([self.index(r) for r in restraints])

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
        if isinstance(targets, Quantity):
            targets = targets.value_in_unit(CHIMERAX_LENGTH_UNIT)
        if hasattr(targets, '__len__'):
            if len(targets) != len(self):
                raise IndexError('Target array length must equal the number of restraints!')
            for r, t in zip(self, targets):
                r.target_distance = t
        else:
            for r in self:
                r.target_distance = target

    @property
    def all_distances(self):
        '''Returns the current distances between restrained atom pairs, in Angstroms.'''
        coords = numpy.reshape(self.coords, [len(self)//2, 2,3])
        return numpy.linalg.norm(coords[:,1]-coords[:,0],axis=1)

    @property
    def spring_constants(self):
        '''Returns the spring constants for all restraints, in kJ/mol/A^3.'''
        return numpy.array([r.spring_constant for r in self.restraints])

    @spring_constants.setter
    def spring_constants(self, val_or_vals):
        if isinstance(val_or_vals, Quantity):
            val_or_vals = val_or_vals.value_in_unit(CHIMERAX_SPRING_UNIT)
        if hasattr(val_or_vals, '__len__'):
            for r, v in zip(self, val_or_vals):
                r.spring_constant = v
        else:
            for r in self:
                r.spring_constant = val_or_vals



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
        if isinstance(k, Quantity):
            k = k.value_in_unit(CHIMERAX_SPRING_UNIT)
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
                ta.hide &= ~HIDE_ISOLDE
                if pb is None:
                    pb = self._pseudobond = self._pbg.new_pseudobond(self._atom, ta)
                pb.display = True
            if self._triggers is not None:
                self._triggers.activate_trigger('position restraint added', self)


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
    
    @spring_constants.setter
    def spring_constants(self, ks):
        if isinstance(ks, Quantity):
            ks = ks.value_in_unit(CHIMERAX_SPRING_UNIT)
        if len(ks) != len(self):
            raise IndexError('Target array length must equal the number of restraints!')
        for r, k in zip(self, ks):
            r.spring_constant = k

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
            index = self.atoms.index(i)
            if index == -1:
                raise TypeError('Not a restrainable atom!')
            return self._restraints[index]
        raise IndexError('Only integer indices allowed for {}, got {}'
            .format(self.__class__.__name__, str(type(i))))

    def index(self, restraint):
        return self.atoms.index(restraint.atom)

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
