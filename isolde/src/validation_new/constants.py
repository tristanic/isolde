from math import pi, radians
import ctypes

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
class _Validation_Defaults:
    _constants = {
        'ROTA_ALLOWED_CUTOFF':      0.02,
        'ROTA_OUTLIER_CUTOFF':      0.0005,
        'CISPRO_OUTLIER':    0.002,
        'TRANSPRO_OUTLIER':  0.001,
        'GLYCINE_OUTLIER':   0.001,
        'PREPRO_OUTLIER':    0.001,
        'ILEVAL_OUTLIER':    0.001,
        'GENERAL_OUTLIER':   0.0005,
        'CISPRO_ALLOWED':    0.02,
        'TRANSPRO_ALLOWED':  0.02,
        'GLYCINE_ALLOWED':   0.02,
        'PREPRO_ALLOWED':    0.02,
        'ILEVAL_ALLOWED':    0.02,
        'GENERAL_ALLOWED':   0.02,
        'MAX_FAVORED_COLOR': [0,255,0,255], # Bright green
        'ALLOWED_COLOR':     [255,240,50,255], # Yellow
        'OUTLIER_COLOR':     [255,0,100,255], # Hot pink
        'NA_COLOR':          [100,100,100,255], # Grey
        'TRACK_RAMACHANDRAN_STATUS':    True,
        'ROUNDS_PER_RAMA_UPDATE':       10,
        'TRACK_ROTAMER_STATUS':         True,
        'ROUNDS_PER_ROTA_UPDATE':       10,
    }

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



validation_defaults = _Validation_Defaults()
ss_restraints = {
    'Helix':                _SS_Helix(),
    'Parallel Beta':        _SS_Beta_Parallel(),
    'Antiparallel Beta':    _SS_Beta_Antiparallel(),
}
