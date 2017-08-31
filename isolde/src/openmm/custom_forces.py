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
        tf = xyz_to_ijk_transform
        # Only the 3x3 rotation/scale portion of the transformation
        # matrix changes with length units. The translation is applied
        # after scaling to the map coordinate system.
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
        map_func = self._discrete3D_from_volume(data)

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

        func = ';'.join((energy_str, norm_str, val_str, min_str, max_str,
                        i_str, j_str, k_str))


        super().__init__(1, func)
        self._map_potential_index = self.addTabulatedFunction(
            name = 'map_potential', function = map_func)
        self._global_k_index = self.addGlobalParameter(
            name = 'global_k', defaultValue = 1.0)
        self._individual_k_index = self.addPerBondParameter(
            name = 'individual_k')
        self.update_needed = False

    def _discrete3D_from_volume(self, data):
        dim = data.shape[::-1]
        data_1d = numpy.ravel(data, order = 'C')
        return Discrete3DFunction(*dim, data_1d)

    def set_global_k(self, k):
        if type(k) == Quantity:
            k = k.value_in_unit(OPENMM_SPRING_UNIT)
        self.setGlobalParameterDefaultValue(self._global_k_index, k)
        self.update_needed = True

    def update_spring_constant(self, index, k):
        params = self.getBondParameters(index)
        atom_i = params[0]
        self.setBondParaeters(index, atom_i, (k,))

    def update_context_if_needed(self, context):
        if self.update_needed:
            self.updateParametersInContext(context)
            self.update_needed = False

class TopOutBondForce(CustomBondForce):
    '''
    Wraps an OpenMM CustomBondForce defined as a standard harmonic potential
    (0.5 * k * (r - r0)^2) with a user-defined fixed maximum cutoff on the
    applied force. This is meant for steering the simulation into new
    conformations where the starting distance may be far from the target
    bond length, leading to catastrophically large forces with a standard
    harmonic potential.
    '''
    def __init__(self, max_force):
        linear_eqn = 'max_force * abs(r-r0)'# - 0.5*max_force^2/k'
        quadratic_eqn = '0.5*k*(r-r0)^2'
        transition_eqn = 'step(r - max_force/k)'
        zero_k_eqn = 'step(min_k - k)'
        energy_str = 'select({},0,select({},{},{}))'.format(
            zero_k_eqn, transition_eqn, linear_eqn, quadratic_eqn)
        super().__init__(energy_str)

        self.k_index = self.addPerBondParameter('k')
        self.r0_index = self.addPerBondParameter('r0')
        self.max_force_index = self.addGlobalParameter('max_force', 0)
        self.max_force = max_force
        self.min_k_index = self.addGlobalParameter('min_k', NEARLY_ZERO)

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

    def update_target(self, index, target=None, k = None):
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
        atom1, atom2 = current_params[0,1]
        new_k, new_target = current_params[2]
        if target is not None:
            if type(target) == Quantity:
                target = target.value_in_unit(OPENMM_LENGTH_UNIT)
            new_target = target
        if k is not None:
            if type(k) == Quantity:
                k = k.value_in_unit(OPENMM_SPRING_UNIT)
            new_k = k
        self.setBondParameters(int(index), atom1, atom2, (new_k, new_t))
        self.update_needed = True




class TopOutRestraintForce(CustomExternalForce):
    '''
    Wraps an OpenMM CustomExternalForce to restrain atoms to defined positions
    via a standard harmonic potential (0.5 * k * r^2) with a user-defined
    fixed maximum cutoff on the appli5ed force. This is meant for steering
    the simulation into new conformations where the starting positions
    may be far from the target positions, leading to catastrophically
    large forces with a standard harmonic potential.
    '''
    def __init__(self, max_force):
        linear_eqn = 'max_force * sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)'
        quadratic_eqn = '0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)'
        transition_eqn = 'step(sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2) - max_force/k)'
        zero_k_eqn = 'step(min_k - k)'
        energy_str = 'select({},0,select({},{},{}))'.format(
            zero_k_eqn, transition_eqn, linear_eqn, quadratic_eqn)

        super().__init__(energy_str)
        per_particle_parameters = ('k','x0','y0','z0')
        for p in per_particle_parameters:
            self.addPerParticleParameter(p)
        self.addGlobalParameter('max_force', 0)
        self.max_force = max_force
        self.addGlobalParameter('min_k', NEARLY_ZERO)

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

    def update_target(self, index, target=None, k = None):
        current_params = self.getParticleParameters(int(index))
        atom_index = current_params[0]
        new_k, new_x, new_y, new_z = current_params[1]
        if target is not None:
            if type(target) == Quantity:
                target = target.value_in_unit(OPENMM_LENGTH_UNIT)
            new_x, new_y, new_z = target
        if k is not None:
            if type(k) == Quantity:
                k = k.value_in_unit(OPENMM_SPRING_UNIT)
            new_k = k
        self.setParticleParameters(int(index), atom_index, (new_k, new_x, new_y, new_z))
        self.update_needed = True

    def release_restraint(self, index):
        self.setParticleParameters(index, (0,0,0,0))


class FlatBottomTorsionRestraintForce(CustomTorsionForce):
    '''
    Wraps an OpenMM CustomTorsionForce to restrain torsion angles while
    allowing free movement within a range (target +/- cutoff). Within
    the cutoff range the potential is constant, while outside it
    is defined as -k * cos (theta-theta0).
    '''
    def __init__(self):
        standard_energy = '-k*cos(theta-theta0)'
        flat_energy = '-k*cos_cutoff'
        switch_function = 'step(cos(theta-theta0)-cos_cutoff)'
        complete_function = 'select({},{},{})'.format(
            switch_function, flat_energy, standard_energy)
        super().__init__(complete_function)
        per_bond_parameters = ('k', 'theta0', 'cos_cutoff')
        for p in per_bond_parameters:
            self.addPerTorsionParameter(p)

        self.update_needed = False

    def update_target(self, index, target = None, k = None, cutoff = None):
        '''
        Change the target, spring constant and/or cutoff angle for the given torsion.
        '''
        # For compatibility with int32
        index = int(index)
        current_params = self.getTorsionParameters(index)
        indices = current_params[0:4]
        new_k, new_theta0, new_cutoff = current_params[4]
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
        self.setTorsionParameters(index, *indices, (new_k, new_theta0, new_cutoff))
        self.update_needed = True


    def set_cutoff_angle(self, i, angle):
        '''
        Set the cut-off angle (below which no force will be applied) in
        radians for one dihedral.
        '''
        if type(angle) == Quantity:
            angle = angle.value_in_unit(OPENMM_ANGLE_UNIT)
        p1, p2, p3, p4, current_params = self.getTorsionParameters(i)
        current_params[2] = angle
        self.setTorsionParameters(i, p1, p2, p3, p4, current_params)
        self.update_needed = True

    def get_cutoff_angle(self, i):
        '''
        Get the cut-off angle in radians for one dihedral.
        '''
        return self.getBondParameters(i)[5][2]*OPENMM_ANGLE_UNIT

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
        
        self.setNonbondedMethod(nonbonded_method)
        self.update_needed = False
