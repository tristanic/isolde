# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


# Custom OpenMM forces

import numpy
from math import pi, radians, degrees, cos
from simtk import unit, openmm
from simtk.unit.quantity import Quantity
from simtk.openmm import app
from simtk.openmm.openmm import CustomBondForce, CustomExternalForce, \
                                CustomCompoundBondForce, CustomTorsionForce, \
                                NonbondedForce, CMAPTorsionForce
from simtk.openmm.openmm import Continuous1DFunction, Continuous3DFunction, \
                                Discrete3DFunction
from simtk.openmm.app.internal import customgbforces
from . import amber_cmap
from ..constants import defaults


import os, sys, glob
import ctypes
from chimerax.core.atomic import molc
from chimerax.core.atomic.molc import CFunctions, string, cptr, pyobject, \
    set_c_pointer, pointer, size_t
from numpy import int8, uint8, int32, uint32, float64, float32, byte, bool as npy_bool
libdir = os.path.dirname(os.path.abspath(__file__))
libfile = glob.glob(os.path.join(libdir, '..', 'openmm.cpython*'))[0]

_c_functions = CFunctions(os.path.splitext(libfile)[0])
c_property = _c_functions.c_property
cvec_property = _c_functions.cvec_property
c_function = _c_functions.c_function
c_array_function = _c_functions.c_array_function





NEARLY_ZERO = defaults.NEARLY_ZERO # Spring constants below this value will be treated as zero

OPENMM_LENGTH_UNIT = defaults.OPENMM_LENGTH_UNIT
OPENMM_FORCE_UNIT = defaults.OPENMM_FORCE_UNIT
OPENMM_SPRING_UNIT = defaults.OPENMM_SPRING_UNIT
OPENMM_RADIAL_SPRING_UNIT = defaults.OPENMM_RADIAL_SPRING_UNIT
OPENMM_ENERGY_UNIT = defaults.OPENMM_ENERGY_UNIT
OPENMM_ANGLE_UNIT = defaults.OPENMM_ANGLE_UNIT
OPENMM_TIME_UNIT = defaults.OPENMM_TIME_UNIT
OPENMM_DIPOLE_UNIT = defaults.OPENMM_DIPOLE_UNIT


class AmberCMAPForce(CMAPTorsionForce):
    '''
    CMAP-style corrections to AMBER 12/14 forcefields to give improved
    backbone conformations in implicit solvent.
    Ref: http://pubs.acs.org/doi/pdf/10.1021/acs.jctc.5b00662
    '''
    _map_loader = amber_cmap.CMAPLoader()
    _cmaps = _map_loader.cmaps

    def __init__(self):
        super().__init__()
        for m in self._cmaps:
            self.addMap(m.shape[0], m.flatten())
        self.update_needed = False

    def add_torsions(self, resnames, phi_indices, psi_indices):
        f = c_function('cmaptorsionforce_add_torsions',
            args = (ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_int32),
                ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32),
                ctypes.POINTER(ctypes.c_int32)))
        ml = self._map_loader
        n = len(resnames)
        map_indices = numpy.array([ml[name] for name in resnames], int32)
        phi_i = numpy.empty((n,4), int32)
        phi_i[:] = phi_indices
        psi_i = numpy.empty((n,4), int32)
        psi_i[:] = psi_indices
        ret = numpy.empty(n, int32)
        f(int(self.this), n, pointer(map_indices), pointer(phi_i), pointer(psi_i),
            pointer(ret))
        return ret

    def addTorsion(self, resname, phi_indices, psi_indices):
        map_index = self._map_loader.map_index(resname)
        super().addTorsion(map_index, *phi_indices.tolist(), *psi_indices.tolist())


class LinearInterpMapForce(CustomCompoundBondForce):
    '''
    Converts a map of (i,j,k) data and a (x,y,z)->(i,j,k) transformation
    matrix to a potential field, with trilinear interpolation of values.
    '''
    def __init__(self, data, xyz_to_ijk_transform, units = 'angstroms'):
        '''
        For a given atom at (x,y,z), the map potential will be defined
        as:
            overall_k * individual_k * pot_xyz

        where pot_xyz is calculated by linear interpolation from the
        (i,j,k) map grid after applying xyz_to_ijk_transform.

        Args:
            data:
                The map data as a 3D array with row order (k,j,i)
            xyz_to_ijk_transform:
                A 3x4 transformation matrix mapping (x,y,z) coordinates
                to (i,j,k)
            units:
                The units in which the transformation matrix is defined.
                Either 'angstroms' or 'nanometers'
        '''
        tf = self._transform = xyz_to_ijk_transform
        # Only the 3x3 rotation/scale portion of the transformation
        # matrix changes with length units. The translation is applied
        # after scaling to the map coordinate system.
        map_func = self._discrete3D_from_volume(data)

        energy_func = self._set_energy_function(tf, units)

        super().__init__(1, energy_func)
        self._map_potential_index = self.addTabulatedFunction(
            name = 'map_potential', function = map_func)
        self._global_k_index = self.addGlobalParameter(
            name = 'global_k', defaultValue = 1.0)
        self._individual_k_index = self.addPerBondParameter(
            name = 'individual_k')
        self._enabled_index = self.addPerBondParameter(
            name = 'enabled')
        self.update_needed = False

    def _set_energy_function(self, tf, units):
        if type(tf) == Quantity:
            rot_and_scale = tf[:,0:3].value_in_unit(OPENMM_LENGTH_UNIT)
            tf = tf.value_in_unit(tf.unit)
            tf[:,0:3] = rot_and_scale

        # TODO: Legacy code. Remove when possible.
        elif units == 'angstroms':
            # OpenMM calcs will be in nm
            tf[:,0:3] *= 10
        elif units != 'nanometers':
            raise TypeError('Units must be either "angstroms" or "nanometers"!')

        # Transform xyz to ijk
        tf_strings = ['i = ', 'j = ', 'k = ']
        entries = ('x1 * {}','y1 * {}','z1 * {}','{}')
        for i in range(3):
            vals = tf[i]
            count = 0
            for (entry, val) in zip(entries, vals):
                if count > 0:
                    spacer = ' + '
                else:
                    spacer = ''
                if abs(val) > NEARLY_ZERO:
                    count += 1
                    tf_strings[i] += spacer + entry.format(val)

        i_str, j_str, k_str = tf_strings


        min_str = 'min_i = floor(i); min_j = floor(j); min_k = floor(k)'
        max_str = 'max_i = ceil(i); max_j = ceil(j); max_k = ceil(k)'

        v000_str = 'v000 = map_potential(min_i, min_j, min_k)'
        v100_str = 'v100 = map_potential(max_i, min_j, min_k)'
        v010_str = 'v010 = map_potential(min_i, max_j, min_k)'
        v001_str = 'v001 = map_potential(min_i, min_j, max_k)'
        v101_str = 'v101 = map_potential(max_i, min_j, max_k)'
        v011_str = 'v011 = map_potential(min_i, max_j, max_k)'
        v110_str = 'v110 = map_potential(max_i, max_j, min_k)'
        v111_str = 'v111 = map_potential(max_i, max_j, max_k)'

        val_str = ';'.join((v000_str, v100_str, v010_str, v001_str,
                            v101_str, v011_str, v110_str, v111_str))

        # Shift i, j and k into the range (0..1)
        norm_str = 'ni = i - min_i; nj = j - min_j; nk = k - min_k'

        interp_str = '(v000*(1-ni)*(1-nj)*(1-nk) + v100* ni   *(1-nj)*(1-nk) +\
                      v010*(1-ni)* nj   *(1-nk) + v001*(1-ni)*(1-nj)* nk    +\
                      v101* ni   *(1-nj)* nk    + v011*(1-ni)* nj   * nk    +\
                      v110* ni   * nj   *(1-nk) + v111* ni   * nj   * nk)'

        energy_str = '-global_k * individual_k * {}'.format(interp_str)

        enabled_eqn = 'step(enabled-0.5)'

        funcs = ';'.join(( norm_str, val_str, min_str, max_str,
                        i_str, j_str, k_str))
        final_func = 'select({}, {}, 0); {}'.format(enabled_eqn, energy_str, funcs)

        return final_func

    def _discrete3D_from_volume(self, data):
        dim = data.shape[::-1]
        data_1d = numpy.ravel(data, order = 'C')
        return Discrete3DFunction(*dim, data_1d)

    def set_global_k(self, k):
        if type(k) == Quantity:
            k = k.value_in_unit(OPENMM_SPRING_UNIT)
        self.setGlobalParameterDefaultValue(self._global_k_index, k)
        self.update_needed = True

    #~ def update_volume_data(self, data, xyz_to_ijk_transform):
        #~ if not numpy.allclose(xyz_to_ijk_transform, self._transform):
            #~ raise TypeError('New map must be on the same grid as the original!')
        #~ dim = data.shape[::-1]
        #~ data_1d = numpy.ravel(data, order = 'C')
        #~ f = self.getTabulatedFunction(0)
        #~ f.setFunctionParameters(*dim, data_1d)
        #~ self.update_needed = True

    def add_atoms(self, indices, ks, enableds):
        f = c_function('customcompoundbondforce_add_bonds',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_int32),
                ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int32)))
        n = len(indices)
        ind = numpy.empty(n, int32)
        ind[:] = indices
        params = numpy.empty((n,2), float64)
        params[:,0] = ks
        params[:,1] = enableds
        ret = numpy.empty(n, int32)
        f(int(self.this), n, pointer(ind), pointer(params), pointer(ret))
        return ret

    def update_atoms(self, indices, coupling_constants, enableds):
        f = c_function('customcompoundbondforce_update_bond_parameters',
            args = (ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_int32),
                ctypes.POINTER(ctypes.c_double)))
        n = len(indices)
        ind = numpy.empty(n, int32)
        ind[:] = indices
        params = numpy.empty((n,2), float64)
        params[:,0] = coupling_constants
        params[:,1] = enableds
        f(int(self.this), n, pointer(ind), pointer(params))
        self.update_needed = True

    def update_atom(self, index, coupling_constant, enabled):
        self.update_atoms([index], [coupling_constant], [enabled])

    def update_context_if_needed(self, context):
        if self.update_needed:
            self.updateParametersInContext(context)
            self.update_needed = False

class TopOutBondForce(CustomBondForce):
    '''
    Wraps an OpenMM CustomBondForce defined as a standard harmonic potential
    (0.5 * k * (r - r0)^2) with a user-defined fixed maximum cutoff on the
    applied force. The force on any atom can be switched on (off) by setting
    the 'enabled' parameter to 1 (0). This is meant for steering the simulation
    into new conformations where the starting distance may be far from the target
    bond length, leading to catastrophically large forces with a standard
    harmonic potential.
    Effective energy equation:

    E = 0                       | enabled < 0.5
        max_force *abs(r-r0)    | r - max_force/k > 0
        0.5 * k * (r - r0)^2    | otherwise
    '''
    def __init__(self, max_force):
        linear_eqn = 'max_force * abs(r-r0)'# - 0.5*max_force^2/k'
        quadratic_eqn = '0.5*k*(r-r0)^2'
        transition_eqn = 'step(r - max_force/k)'
        enabled_eqn = 'step(enabled - 0.5)'
        energy_str = 'select({},select({},{},{}),0)'.format(
            enabled_eqn, transition_eqn, linear_eqn, quadratic_eqn)
        super().__init__(energy_str)

        self.enabled_index = self.addPerBondParameter('enabled')
        self.k_index = self.addPerBondParameter('k')
        self.r0_index = self.addPerBondParameter('r0')
        self.max_force_index = self.addGlobalParameter('max_force', 0)
        self.max_force = max_force

        self.update_needed = False

    @property
    def max_force(self):
        '''Maximum force applied to any given atom, in kJ/mol/nm.'''
        return self._max_force

    @max_force.setter
    def max_force(self, force):
        if type(force) == Quantity:
            force = force.value_in_unit(OPENMM_FORCE_UNIT)
        self.setGlobalParameterDefaultValue(self.max_force_index, force)
        self._max_force = force
        self.update_needed = True

    def add_bonds(self, atom_indices, enableds, spring_constants, targets):
        '''
        Add a set of bonds
        @param atom_indices:
            A 2-tuple of arrays giving the indices of the bonded atoms in the
            simulation construct
        @param targets:
            Target distances in nanometers
        @param spring_constants:
            Spring constants in kJ mol-1 nm-2
        '''
        f = c_function('custombondforce_add_bonds',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_int32),
                ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int32)))
        n = len(targets)
        ind = numpy.empty((n,2), int32)
        for i, ai in enumerate(atom_indices):
            ind[:,i] = ai
        params = numpy.empty((n,3), float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2] = targets
        ret = numpy.empty(n, int32)
        f(int(self.this), n, pointer(ind), pointer(params), pointer(ret))
        return ret


    def update_target(self, index, enabled=None, k = None, target=None):
        '''
        Update an existing distance restraint with a new target and/or spring
        constant.
        @param index:
            The index of this restraint in the OpenMM force object
        @param target (default None):
            The new target distance (as a simtk Quantity or in nanometres).
            Leave as None to preserve the existing target distance
        @param k (default None):
            The new spring constant (as a simtk Quantity or in units of
            kJ mol-1 nm-2). Leave as None to preserve the existing spring
            constant.
        '''
        current_params = self.getBondParameters(int(index))
        atom1, atom2 = current_params[0:2]
        new_enabled, new_k, new_target = current_params[2]
        if enabled is not None:
            new_enabled = float(enabled)
        if target is not None:
            if type(target) == Quantity:
                target = target.value_in_unit(OPENMM_LENGTH_UNIT)
            new_target = target
        if k is not None:
            if type(k) == Quantity:
                k = k.value_in_unit(OPENMM_SPRING_UNIT)
            new_k = k
        self.setBondParameters(int(index), atom1, atom2, (new_enabled, new_k, new_target))
        self.update_needed = True

    def update_targets(self, indices, enableds, spring_constants, targets):
        '''
        Update a set of targets all at once using fast C++ code. Fastest if
        the arguments are provided as Numpy arrays, but any iterable will work.
        @param indices:
            The indices of the restraints in the OpenMM force object
        @param targets:
            The new target distances in nanometres
        @param k
            The new spring constants in units of kJ mol-1 nm-2
        '''
        f = c_function('custombondforce_update_bond_parameters',
            args=(ctypes.c_void_p, ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_double)))
        n = len(indices)
        ind = numpy.empty(n, int32)
        ind[:] = indices
        params = numpy.empty((n,3), float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2] = targets
        f(int(self.this), n, pointer(ind), pointer(params))
        self.update_needed = True


class TopOutRestraintForce(CustomExternalForce):
    '''
    Wraps an OpenMM CustomExternalForce to restrain atoms to defined positions
    via a standard harmonic potential (0.5 * k * r^2) with a user-defined
    fixed maximum cutoff on the appli5ed force. This is meant for steering
    the simulation into new conformations where the starting positions
    may be far from the target positions, leading to catastrophically
    large forces with a standard harmonic potential.
    Effective energy equation:

    E = 0                      | enabled < 0.5
        max_force * r          | r - max_force/k > 0
        0.5 * k * r^2          | otherwise

    '''
    def __init__(self, max_force):
        linear_eqn = 'max_force * sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)'
        quadratic_eqn = '0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)'
        transition_eqn = 'step(sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2) - max_force/k)'
        enabled_eqn = 'step(enabled - 0.5)'
        energy_str = 'select({},select({},{},{}),0)'.format(
            enabled_eqn, transition_eqn, linear_eqn, quadratic_eqn)

        super().__init__(energy_str)
        per_particle_parameters = ('enabled','k','x0','y0','z0')
        for p in per_particle_parameters:
            self.addPerParticleParameter(p)
        self.addGlobalParameter('max_force', 0)
        self.max_force = max_force

        self.update_needed = False

    @property
    def max_force(self):
        '''Maximum force applied to any given atom, in kJ/mol/nm.'''
        return self._max_force

    @max_force.setter
    def max_force(self, force):
        if type(force) == Quantity:
            force = force.value_in_unit(OPENMM_FORCE_UNIT)
        self.setGlobalParameterDefaultValue(0, force)
        self._max_force = force
        self.update_needed = True

    def add_particles(self, indices, enableds, spring_constants, targets):
        f = c_function('customexternalforce_add_particles',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_int32),
                ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int32)))
        n = len(indices)
        ind = numpy.empty(n, int32)
        ind[:] = indices
        ret = numpy.empty(n, int32)
        params = numpy.empty((n,5), float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2:] = targets
        f(int(self.this), n, pointer(ind), pointer(params), pointer(ret))
        return ret

    def update_target(self, index, enabled=None, k=None, target=None):
        current_params = self.getParticleParameters(int(index))
        atom_index = current_params[0]
        new_enabled, new_k, new_x, new_y, new_z = current_params[1]
        if enabled is not None:
            new_enabled = float(enabled)
        if target is not None:
            if type(target) == Quantity:
                target = target.value_in_unit(OPENMM_LENGTH_UNIT)
            new_x, new_y, new_z = target
        if k is not None:
            if type(k) == Quantity:
                k = k.value_in_unit(OPENMM_SPRING_UNIT)
            new_k = k
        self.setParticleParameters(int(index), atom_index, (new_enabled, new_k, new_x, new_y, new_z))
        self.update_needed = True

    def update_targets(self, indices, enableds, spring_constants, targets):
        f = c_function('customexternalforce_update_particle_parameters',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_int32),
                ctypes.POINTER(ctypes.c_double)))
        n = len(indices)
        if len(targets) !=n or len(spring_constants) !=n:
            raise TypeError('Parameter array lengths must match number of indices!')
        ind = numpy.empty(n, int32)
        ind[:] = indices
        params = numpy.empty((n,5), float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2:] = targets
        f(int(self.this), n, pointer(ind), pointer(params))
        self.update_needed = True

    def release_restraint(self, index):
        self.update_target(index, enabled=False)


class FlatBottomTorsionRestraintForce(CustomTorsionForce):
    '''
    Wraps an OpenMM CustomTorsionForce to restrain torsion angles while
    allowing free movement within a range (target +/- cutoff). Within
    the cutoff range the potential is constant, while outside it
    is defined as -k * cos (theta-theta0).
    Effective energy function:

    E = 0                    | enabled < 0.5
        -k*cos(cutoff)       | cos(theta-theta0) - cos(cutoff) < 0
        -k*cos(theta-theta0) | otherwise
    '''
    def __init__(self):
        standard_energy = '-k*cos(theta-theta0)'
        flat_energy = '-k*cos_cutoff'
        switch_function = 'step(cos(theta-theta0)-cos_cutoff)'
        enabled_function = 'step(enabled-0.5)'
        complete_function = 'select({},select({},{},{}),0)'.format(
            enabled_function, switch_function, flat_energy, standard_energy)
        super().__init__(complete_function)
        per_bond_parameters = ('enabled', 'k', 'theta0', 'cos_cutoff')
        for p in per_bond_parameters:
            self.addPerTorsionParameter(p)

        self.update_needed = False

    def add_torsions(self, atom_indices, enableds, spring_constants, targets, cutoffs):
        '''
        Add a set of torsion restraints. Returns an array of ints representing
        the indices of the restraints in the force object.
        @param atom_indices:
            A 4-tuple of arrays providing the indices of the dihedral atoms in
            the simulation construct
        @param enableds:
            A boolean array (or any array castable to float) where values > 0.5
            represent enabled
        @param spring_constants:
            Restraint spring constants in kJ mol-1 rad-2
        @param targets:
            Target angles in radians
        @param cutoffs:
            Cutoff angle (below which no restraint is applied) for each restraint
            in radians.
        '''
        f = c_function('customtorsionforce_add_torsions',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_int32),
                ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int32)))
        n = len(targets)
        ind = numpy.empty((n,4), numpy.int32)
        for i, ai in enumerate(atom_indices):
            ind[:,i] = ai
        params = numpy.empty((n,4), float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2] = targets
        params[:,3] = numpy.cos(cutoffs)
        ret = numpy.empty(n, numpy.int32)
        f(int(self.this), n, pointer(ind), pointer(params), pointer(ret))
        return ret


    def update_target(self, index, enabled=None, k = None, target = None, cutoff = None):
        '''
        Change the target, spring constant and/or cutoff angle for the given torsion.
        '''
        # For compatibility with int32
        index = int(index)
        current_params = self.getTorsionParameters(index)
        indices = current_params[0:4]
        new_enabled, new_k, new_theta0, new_cutoff = current_params[4]
        if enabled is not None:
            new_enabled = float(enabled)
        if target is not None:
            if type(target) == Quantity:
                target = target.value_in_unit(OPENMM_ANGLE_UNIT)
            new_theta0 = target
        if k is not None:
            if type(k) == Quantity:
                k = k.value_in_unit(OPENMM_RADIAL_SPRING_UNIT)
            new_k = k
        if cutoff is not None:
            if type(cutoff) == Quantity:
                cutoff = cutoff.value_in_unit(OPENMM_ANGLE_UNIT)
            new_cutoff = cos(cutoff)
        self.setTorsionParameters(index, *indices, (new_enabled, new_k, new_theta0, new_cutoff))
        self.update_needed = True

    def update_targets(self, indices, enableds, spring_constants, targets, cutoffs):
        '''
        Change the target angles, spring constants and cutoff angles for a set
        of dihedral restraints.
        @param indices:
            the indices for each restraint in the force object
        @param enableds:
            a boolean, int or float array where values >0.5 represent enabled
        @param spring_constants:
            the spring constants for the restraints in kJ mol-1 rad-2
        @param targets:
            the target angles in radians
        @param cutoffs:
            the cutoff angles (below which no force is applied) in radians
        '''
        f = c_function('customtorsionforce_update_torsion_parameters',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_int32),
                ctypes.POINTER(ctypes.c_double)))
        n = len(indices)
        ind = numpy.empty(n, int32)
        ind[:] = indices
        params = numpy.empty((n,4), float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2] = targets
        params[:,3] = numpy.cos(cutoffs)
        f(int(self.this), n, pointer(ind), pointer(params))
        self.update_needed = True


class TorsionNCSForce(CustomCompoundBondForce):
    '''
    Provides torsion-angle non-crystallographic symmetry (NCS)
    restraints for a defined number of NCS copies. For a given set of
    NCS-equivalent dihedrals, each dihedral will receive a scaling term
    (defining how strongly it is restrained towards the weighted vector
    mean of all the dihedral angles) and a weighting term (defining how
    much it contributes to the weighted mean angle). Setting the
    weights for all but one dihedral to zero will yield a
    "master-slave" mode where all NCS copies are forced to follow the
    one with the non-zero weight.
    '''
    def __init__(self, num_copies):
        '''
        Initialise a torsion-angle non-crystallographic symmetry (NCS)
        restraint force for a given number of NCS copies.
        '''
        if num_copies < 2:
            raise TypeError('NCS is only applicable if you have multiple '\
                +'equivalent chains!')
        dihedral_def_strings = []
        k_params = []
        w_params = []
        sum_sin_string = 'sum_sin = ('
        sum_cos_string = 'sum_cos = ('
        energy_string = 'global_k *('
        for i in range(num_copies):
            d_particles = []
            for j in range(4):
                d_particles.append('p{}'.format(4*i+j+1))
            d_def = 'd{} = dihedral({})'.format(i+1, ','.join(d_particles))
            dihedral_def_strings.append(d_def)
            d = 'd{}'.format(i+1)
            k = 'k{}'.format(i+1)
            k_params.append(k)
            w = 'w{}'.format(i+1)
            w_params.append(w)

            if i != 0:
                sum_sin_string += '+'
                sum_cos_string += '+'
                energy_string += '+'
            sum_sin_string += '{}*sin({})'.format(w, d)
            sum_cos_string+='{}*cos({})'.format(w, d)
            energy_string += '{}*(cos({})-target)'.format(k, d)

        sum_sin_string += ')'
        sum_cos_string += ')'
        energy_string += ')'

        cos_switch_str = 'cos_switch = step(sum_cos)'
        sin_switch_str = 'sin_switch = step(sum_sin)'
        quadrant_select = \
            'offset = select(cos_switch, 0, select(sin_switch, {}, {}))'\
                .format(str(pi/2), str(-pi/2))

        target_string = 'target = cos(atan(sum_sin/sum_cos) + offset)'

        string_list = [
            energy_string,
            target_string,
            quadrant_select,
            cos_switch_str,
            sin_switch_str,
            sum_sin_string,
            sum_cos_string,
            ]
        string_list.extend(dihedral_def_strings)


        final_eq = ';'.join(string_list)
        super().__init__(num_copies*4, final_eq)
        self._weight_indices = []
        self._k_indices = []
        for w in w_params:
            self._weight_indices.append(self.addPerBondParameter(w))
        for k in k_params:
            self._k_indices.append(self.addPerBondParameter(k))
        self.addGlobalParameter('global_k', 1.0)



class GBSAForce(customgbforces.GBSAGBn2Force):
    def __init__(self, solventDielectric=78.5, soluteDielectric=1,
                SA='ACE', cutoff=1.0, kappa=3.0,
                nonbonded_method = openmm.CustomGBForce.CutoffNonPeriodic):
        '''
        kappa = 3.0/nm --> approx. 0.5M ion concentration at 100K
        '''
        if type(solventDielectric) == Quantity:
            solventDielectric = solventDielectric.value_in_unit(OPENMM_DIPOLE_UNIT)
        if type(soluteDielectric) == Quantity:
            soluteDielectric = soluteDielectric.value_in_unit(OPENMM_DIPOLE_UNIT)
        if type(cutoff) == Quantity:
            cutoff = cutoff.value_in_unit(OPENMM_LENGTH_UNIT)
        if type(kappa) == Quantity:
            kappa = kappa.value_in_unit(1/OPENMM_LENGTH_UNIT)


        super().__init__(solventDielectric=solventDielectric,
                         soluteDielectric=soluteDielectric,
                         SA=SA,
                         cutoff=cutoff,
                         kappa=kappa)

        self.setCutoffDistance(cutoff)

        self.setNonbondedMethod(nonbonded_method)
        self.update_needed = False
