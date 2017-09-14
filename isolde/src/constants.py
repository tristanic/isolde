from simtk import unit, openmm
from simtk.openmm import app
from math import pi
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
class _Defaults:
    _constants = {
        ###
        # Units
        ###
        'OPENMM_LENGTH_UNIT':         unit.nanometer,
        'OPENMM_FORCE_UNIT':          unit.kilojoule_per_mole/unit.nanometer,
        'OPENMM_SPRING_UNIT':         unit.kilojoule_per_mole/unit.nanometer**2,
        'OPENMM_RADIAL_SPRING_UNIT':  unit.kilojoule_per_mole/unit.radians**2,
        'OPENMM_ENERGY_UNIT':         unit.kilojoule_per_mole,
        'OPENMM_ANGLE_UNIT':          unit.radians,
        'OPENMM_TIME_UNIT':           unit.picoseconds,
        'OPENMM_DIPOLE_UNIT':         unit.debye,
        'OPENMM_TEMPERATURE_UNIT':    unit.kelvin,

        ###
        # Simulation parameters
        ###
        'OPENMM_PLATFORM':            'OpenCL', # 'CPU', 'OpenCL' or 'CUDA'
        'OPENMM_INTEGRATOR_TYPE':     openmm.VariableLangevinIntegrator,
        'OPENMM_NONBONDED_METHOD':    app.CutoffNonPeriodic,
        'OPENMM_NONBONDED_CUTOFF':    1.0, #* unit.nanometers,
        'OPENMM_FRICTION':            5.0, #/ unit.picoseconds,
        'OPENMM_VAR_INTEGRATOR_TOL':  1e-4,
        'OPENMM_CONSTRAINT_TOL':      1e-4,
        'OPENMM_FIXED_INTEGRATOR_TS': 0.001, #* unit.picoseconds,
        'SIM_STEPS_PER_GUI_UPDATE':   20,
        'MIN_STEPS_PER_GUI_UPDATE':   50,
        'SIM_STARTUP_ROUNDS':         10,
        'MAX_UNSTABLE_ROUNDS':        20,
        'TEMPERATURE':                100.0, # * unit.kelvin,
        'USE_GBSA':                   True,
        'GBSA_NONBONDED_METHOD':      openmm.CustomGBForce.CutoffNonPeriodic,
        'GBSA_SOLVENT_DIELECTRIC':    78.5, # * unit.debye,
        'GBSA_SOLUTE_DIELECTRIC':     1.0, # * unit.debye,
        'GBSA_SA_METHOD':             'ACE', # alternative is None
        'GBSA_CUTOFF':                1.0, # *unit.nanometer,
        'GBSA_KAPPA':                 3.0, # /unit.nanometer,
        'VACUUM_DIELECTRIC_CORR':     150, # *unit.debye,
        'RIGID_BONDS':                app.HBonds,
        'RIGID_WATER':                True,
        'REMOVE_C_OF_M_MOTION':       False,
        'MIN_CONVERGENCE_FORCE_TOL':  50.0, # *unit.kilojoule_per_mole/unit.nanometer


        ###
        # Force constants
        ###
        'MAX_RESTRAINT_FORCE':                  25000.0, # * unit.kilojoule_per_mole/unit.nanometer,
        'HAPTIC_SPRING_CONSTANT':                2500.0, # * unit.kilojoule_per_mole/unit.nanometer**2,
        'MOUSE_TUG_SPRING_CONSTANT':            10000.0, # * unit.kilojoule_per_mole/unit.nanometer**2,
        'MAX_TUG_FORCE':                        10000.0, # * unit.kilojoule_per_mole/unit.nanometer,
        'DISTANCE_RESTRAINT_SPRING_CONSTANT':   10000.0, # * unit.kilojoule_per_mole/unit.nanometer**2,
        'POSITION_RESTRAINT_SPRING_CONSTANT':   10000.0, # * unit.kilojoule_per_mole/unit.nanometer**2,
        'PEPTIDE_SPRING_CONSTANT':                500.0, # * unit.kilojoule_per_mole/unit.radians**2,
        'PHI_PSI_SPRING_CONSTANT':                500.0, # * unit.kilojoule_per_mole/unit.radians**2,
        'ROTAMER_SPRING_CONSTANT':               1000.0, # * unit.kilojoule_per_mole/unit.radians**2,
        'STANDARD_MAP_K':                          50.0, # * unit.kilojoule_per_mole*unit.angstrom**3,
        'DIFFERENCE_MAP_K':                        10.0, # * unit.kilojoule_per_mole*unit.angstrom**3,

        ###
        # Numeric limits
        ###
        'MAX_ALLOWABLE_FORCE':        4.0e4, # * unit.kilojoule_per_mole/unit.nanometer,
        'MAX_STABLE_FORCE':            5000, # * unit.kilojoule_per_mole/unit.nanometer,
        'MAX_ATOM_MOVEMENT_PER_STEP':  0.015, # * unit.nanometer,
        'NEARLY_ZERO':                 1e-6,

        ###
        # Geometric parameters
        ###
        'DIHEDRAL_RESTRAINT_CUTOFF':    pi/6,  # * unit.radians, # 30 degrees
        'ROTAMER_RESTRAINT_CUTOFF':     pi/12, # * unit.radians, # 15 degrees
        'CIS_PEPTIDE_BOND_CUTOFF':      pi/6,  # * unit.radians, # 30 degrees
        'TWISTED_PEPTIDE_BOND_DELTA':  pi/6,  # * unit.radians,

        ###
        # Constants specific to ChimeraX - lengths in Angstroms beyond this point!
        ###

        'CHIMERAX_LENGTH_UNIT':         unit.angstrom,
        'CHIMERAX_FORCE_UNIT':          unit.kilojoule_per_mole/unit.angstrom,
        'CHIMERAX_SPRING_UNIT':         unit.kilojoule_per_mole/unit.angstrom**2,

        'SELECTION_SEQUENCE_PADDING':   3,      # residues
        'SOFT_SHELL_CUTOFF':            5,      # Angstroms
        'HARD_SHELL_CUTOFF':            8,      # Angstroms
        'FIX_SOFT_SHELL_BACKBONE':      False,
        'TRACK_RAMACHANDRAN_STATUS':    True,
        'ROUNDS_PER_RAMA_UPDATE':       20,
        'REMASK_MAPS_DURING_SIM':       True,
        'ROUNDS_PER_MAP_REMASK':        50,
        'HIDE_SURROUNDINGS_DURING_SIM': True,
        'USE_HAPTIC_DEVICE_IF_AVAIL':   True,

        'STANDARD_MAP_MASK_RADIUS':   4.0, # Angstroms
        'DIFFERENCE_MAP_MASK_RADIUS': 8.0, # Angstroms

        'COMMS_TIMEOUT':              0.5, # seconds

        ###
        # Types for shared variables
        ###
        'FLOAT_TYPE' :                  ctypes.c_double,
    }

@_constant_properties
class _Sim_Outcomes:
    _constants = {
        'COMMIT':         0,
        'DISCARD':        1,
        'FAIL_TO_START':  2,
        'UNSTABLE':       3,
        'TIMEOUT':        4,
        'ERROR':          5,
        }

@_constant_properties
class _Control:
    _constants = {
        'HIDE_ISOLDE':    0x2,  # For Atom.hide bitmask, dedicated to ISOLDE

        'SIM_MODE_MIN':               0,
        'SIM_MODE_EQUIL':             1,
        'SIM_MODE_UNSTABLE':          2,

    }



defaults = _Defaults()
sim_outcomes = _Sim_Outcomes()
control = _Control()
