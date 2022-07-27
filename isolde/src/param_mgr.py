# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



from openmm import unit
from openmm.unit import Quantity, Unit

def _property_factory(name):
    def fget(self):
        return self._params[name]
    def fset(self, val):
        self.set_param(name, val)
    return property(fget, fset)

def autodoc(cls):
    _init_docstring = '''
        Initialise the simulation parameters, specifying any non-default
        values.\n\nParams:\n\n'''
    _param_list = ''
    for key, val in cls._default_params.items():
        if type(val[0]) == float:
            val0 = '{:0.4f}'.format(val[0])
        else:
            val0 = val[0]
        if val[1] is None:
            val1 = ''
        else:
            val1 = val[1]
        _param_list += '\t* {:>}: \t {} {}\n'.format(key, val0, val1).expandtabs(20)
    cls.__init__.__doc__ = _init_docstring+_param_list
    return cls

def param_properties(cls):
    for key in cls._default_params.keys():
        setattr(cls, key, _property_factory(key))
    return cls

# When subclassing, apply these decorators in this order to generate the
# __init__ docstring and convert each entry into a class property.
#@param_properties
#@autodoc
class Param_Mgr:
    '''
    A manager for important parameters to keep them all in one place, rather than
    spread as variables throughout a large class. On object initialisation, the
    docstring will be automatically generated to list the parameters and their
    default values (and units where applicable). Addition of new parameters after
    initialisation is not allowed. Once initialised, parameter values may be
    changed by either (for example):

    `obj.set_param(name, new_value)`

    or

    `obj.name = new_value`

    Values are retrievable either by:

    `val = obj.name`

    or

    `val = obj[name]`
    '''
    PARAMETER_CHANGED = 'parameter changed'
    # In the derived class, fill _default_params with entries of the form:
    #    'parameter_name': (default_value, units)
    # where units is either a simtk.unit or None. In the former case, the
    # value will be automatically combined with the unit on object initialisation,
    # and changing the value will respect units (or assume the same units if
    # none are provided with the replacement value).
    _default_params = {}
    def __init__(self, **kw):
        from chimerax.core.triggerset import TriggerSet
        self.triggers = TriggerSet()
        self.triggers.add_trigger(self.PARAMETER_CHANGED)
        self._params = {}
        for key, item in self._default_params.items():
            if type(item[1]) == Unit:
                self._params[key] = item[0]*item[1]
            else:
                self._params[key] = item[0]
        for key, val in kw.items():
            self.set_param(key, val)

    def __getitem__(self, key):
        return self._params[key]

    def __setitem__(self, key, val):
        raise KeyError('Parameters should not be set directly! please use the '+
            'set_param() function.')

    def param_names(self):
        return sorted(self._params.keys())

    def __repr__(self):
        return self._params.__repr__()

    def __str__(self):
        return self._params.__str__()

    def set_param(self, key, value):
        '''
        Set the value of a parameter. If the parameter is a numerical
        quantity with a unit, you may choose to set it directly as a
        number, or in an equivalent :class:`simtk.unit.Unit`. For example:

        `set_param('dihedral_restraint_cutoff_angle', pi/6)`

        and

        `set_param('dihedral_restraint_cutoff_angle', 30 * unit.degrees)`

        will yield the same result, but

        `set_param('dihedral_restraint_cutoff_angle', 30 * unit.nanometer)`

        will fail.
        '''
        try:
            units = self._default_params[key][1]
        except KeyError:
            raise KeyError('Unrecognised parameter!')
        if units is not None and type(units) == Unit:
            if type(value) == Quantity:
                self._params[key] = value.in_units_of(units)
            else:
                self._params[key] = value * units
        elif type(value) == Quantity:
            raise TypeError('Tried to set a unitless quantity with units!')
        else:
            self._params[key] = value
        self.triggers.activate_trigger(self.PARAMETER_CHANGED, (key, self._params[key]))
    
    def set_to_default(self, key):
        '''Set one parameter back to the default value.'''
        self.set_param(key, self._defaults[key][0])

    def reset_to_defaults(self):
        '''Reset all parameters to defaults.'''
        for key in self.param_names():
            self.set_to_default(key)
