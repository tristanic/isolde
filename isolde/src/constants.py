

_defaults = {
    ###
    # Force constants
    ###
    'MAX_RESTRAINT_FORCE':          250,    # kJ/mol/A**2
    'HAPTIC_SPRING_CONSTANT':      2500,    # kJ/mol/A**2
    'PEPTIDE_SPRING_CONSTANT':      500,    # kJ/mol/radian**2
    'PHI_PSI_SPRING_CONSTANT':      250,    # kJ/mol/radian**2 
    'ROTAMER_SPRING_CONSTANT':      1000,   # kJ/mol/radian**2
    
    ###
    # Geometric parameters
    ###
    'DIHEDRAL_RESTRAINT_CUTOFF':    30,     # degrees
    'ROTAMER_RESTRAINT_CUTOFF':     15,     # degrees
    'CIS_PEPTIDE_BOND_CUTOFF':      30,     # degrees
    
    'SELECTION_SEQUENCE_PADDING':   3,      # residues
    'SOFT_SHELL_CUTOFF':            5,      # Angstroms
    'HARD_SHELL_CUTOFF':            8,      # Angstroms
    'FIX_SOFT_SHELL_BACKBONE':      False,  
    
    
}

class _Defaults:
    pass

def _fset_factory():
    def fset(self, value):
        raise TypeError("Can't change value of a constant!")
    return fset

def _fget_factory(val):
    def fget(self):
        return val
    return fget

for (name, val) in _defaults.items():
    setattr(_Defaults, name, property(_fget_factory(val), _fset_factory()))

defaults = _Defaults()

        
    


