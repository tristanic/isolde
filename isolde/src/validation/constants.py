# @Author: Tristan Croll
# @Date:   11-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   Tristan Croll
# @Last modified time: 18-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



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




validation_defaults = _Validation_Defaults()
