import numpy
from chimerax.core.atomic import Residue, Atoms, AtomicStructure
from .restraints_base import Position_Restraint, Position_Restraints   

MAX_RESTRAINT_FORCE = 10.0 # kJ/mol/A^3 

class Atom_Position_Restraints(Position_Restraints):
    '''
    For restraining of atoms to specific points in space.
    '''
    def __init__(self, atoms_or_list, include_hydrogens = False):
        '''
        Initialises restraint objects for the atoms in the given selection.
        By default, hydrogens are ignored.
        '''
        if type(atoms_or_list) == Atoms:
            if not include_hydrogens:
                atoms = atoms_or_list
                atoms = atoms.filter(atoms.element_names != 'H')
            restraints = []
            for a in atoms:
                restraints.append(Position_Restraint(a))
        else:
            restraints = atoms_or_list
        super().__init__(restraints)
    
    
