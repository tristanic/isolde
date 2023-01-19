# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



from math import pi, radians

'''
Constants are a slightly difficult problem here, in that ChimeraX works
in Angstroms while OpenMM works in nanometres
'''

def _constant_property_factory(key):
    def fget(self):
        return self._constants[key]

    def fset(self, val):
        raise TypeError("Can't change value of a constant!")

    return property(fget, fset)

def _constant_properties(cls):
    for key in cls._constants.keys():
        setattr(cls, key, _constant_property_factory(key))
    return cls


@_constant_properties
class _SS_Helix:
    _constants = {
        'CA_TO_CA_PLUS_TWO_DISTANCE':        5.43,
        'O_TO_N_PLUS_FOUR_DISTANCE':         3.05,
        'PHI_ANGLE':                         radians(-64.0),
        'PSI_ANGLE':                         radians(-41.0),
        'CUTOFF_ANGLE':                      radians(10.0),
    }

@_constant_properties
class _SS_Beta_Parallel:
    _constants = {
        'CA_TO_CA_PLUS_TWO_DISTANCE':        6.81,
        'O_TO_N_PLUS_FOUR_DISTANCE':         11.4,
        'PHI_ANGLE':                         radians(-119.0),
        'PSI_ANGLE':                         radians(113.0),
        'CUTOFF_ANGLE':                      radians(10.0),
    }

@_constant_properties
class _SS_Beta_Antiparallel:
    _constants = {
        'CA_TO_CA_PLUS_TWO_DISTANCE':        6.81,
        'O_TO_N_PLUS_FOUR_DISTANCE':         11.4,
        'PHI_ANGLE':                         radians(-139.0),
        'PSI_ANGLE':                         radians(135.0),
        'CUTOFF_ANGLE':                      radians(10.0),
    }

@_constant_properties
class _SS_Beta_Generic:
    _constants = {
        'CA_TO_CA_PLUS_TWO_DISTANCE':        6.81,
        'O_TO_N_PLUS_FOUR_DISTANCE':         11.4,
        'PHI_ANGLE':                         radians(-129.0),
        'PSI_ANGLE':                         radians(124.0),
        'CUTOFF_ANGLE':                      radians(30.0),
        'H_BOND_LENGTH':                     2.9,
    }

ss_restraints = {
    'helix':                _SS_Helix(),
    'strand':               _SS_Beta_Generic(),
}

@_constant_properties
class _DefaultColors:
    _constants = {
        'PROPER_DIHEDRAL_RESTRAINT_SATISFIED_COLOR': [0,255,0,255], # Bright green
        'PROPER_DIHEDRAL_RESTRAINT_STRAINED_COLOR':     [255,240,50,255], # Yellow
        'PROPER_DIHEDRAL_RESTRAINT_SEVERE_COLOR':     [255,0,100,255], # Hot pink

        'DISTANCE_RESTRAINT_BOND_COLOR':                [168, 255, 230, 255],
        'DISTANCE_RESTRAINT_TARGET_COLOR':              [128, 215, 190, 255],        

        'ADAPTIVE_DIHEDRAL_RESTRAINT_SATISFIED_COLOR':  [0,191,255,255],
        'ADAPTIVE_DIHEDRAL_RESTRAINT_STRAINED_COLOR':   [255,69,0,255],
        'ADAPTIVE_DIHEDRAL_RESTRAINT_SEVERE_COLOR':     [139,0,139,255],

        'ADAPTIVE_DISTANCE_RESTRAINT_SATISFIED_COLOR':  [0, 255, 0, 255],  # Bright green
        'ADAPTIVE_DISTANCE_RESTRAINT_TOO_CLOSE_COLOR':  [204, 204, 0, 255], # Yellow
        'ADAPTIVE_DISTANCE_RESTRAINT_TOO_FAR_COLOR':    [102, 0, 204, 255], # Purple

        'POSITION_RESTRAINT_PIN_COLOR':                 [255,215,0,255],
        'POSITION_RESTRAINT_BOND_COLOR':                [200, 250, 120, 255],

        'TUGGING_ARROW_COLOR':                          [100, 255, 100, 255],

    }

restraint_color_defaults = _DefaultColors()