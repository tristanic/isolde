# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 17-Sep-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



# Custom OpenMM forces

import numpy
from math import pi, radians, degrees, cos
from openmm import unit, openmm
from openmm.unit.quantity import Quantity
from openmm import app
from openmm.openmm import CustomBondForce, CustomExternalForce, \
                                CustomCompoundBondForce, CustomTorsionForce, \
                                CMAPTorsionForce, CustomNonbondedForce
from openmm.openmm import Continuous1DFunction, Continuous3DFunction, \
                                Discrete2DFunction, Discrete3DFunction
from openmm.app.internal import customgbforces
from . import amber_cmap
from ..constants import defaults


def _strip_units(x):
    if unit.is_quantity(x):
        return x.value_in_unit_system(unit.md_unit_system)
    return x


from chimerax.isolde import _openmm_force_ext



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
        '''
        Add a set of peptide backbone (phi, psi) pairs to the force.

        Args:
            * resnames:
                - an iterable of upper-case three-letter residue names
            * phi_indices:
                - a (nx4) NumPy integer array giving the indices of the atoms
                  from each phi dihedral in the simulation
            * psi_indices:
                - a (nx4) NumPy integer array giving the indices of the atoms
                  from each psi dihedral in the simulation
        '''

        ml = self._map_loader
        map_indices = numpy.array([ml[name] for name in resnames], numpy.int32)
        return _openmm_force_ext.cmaptorsionforce_add_torsions(int(self.this), map_indices, phi_indices, psi_indices)

    def addTorsion(self, resname, phi_indices, psi_indices):
        '''
        Add a single phi/psi pair to the force.

        Args:
            * resname:
                - the upper-case three-character name of the amino acid residue
            * phi_indices:
                - a NumPy array of four ints giving the indices of the phi
                  dihedral atoms in the simulation
            * psi_indices:
                - a NumPy array of four ints giving the indices of the psi
                  dihedral atoms in the simulation
        '''
        map_index = self._map_loader.map_index(resname)
        super().addTorsion(map_index, *phi_indices.tolist(), *psi_indices.tolist())


class _Map_Force_Base(CustomCompoundBondForce):
    '''
    Base class for :class:`LinearInterpMapForce`,
    :class:`CubicInterpMapForce` and :class:`CubicInterpMapForce_Low_Memory`.
    '''
    def __init__(self, data, xyz_to_ijk_transform, suffix, units = 'angstroms'):
        '''
        For a given atom at (x,y,z), the map potential will be defined
        as:

        .. math::
            global_k * individual_k * pot_xyz

        where `pot_xyz` is calculated by interpolation from the
        (i,j,k) map grid after applying `xyz_to_ijk_transform`.

        Args:
            * data:
                - The map data as a 3D (i,j,k) NumPy array in C-style order
            * xyz_to_ijk_transform:
                - A NumPy 3x4 float array defining the transformation matrix
                  mapping (x,y,z) coordinates to (i,j,k)
            * suffix:
                - In OpenMM, global parameters are global to the *entire
                  context*, not just to the force. To provide a global
                  parameter unique to this instance, the suffix is appended
                  to the base name of the parameter. Should be a unique string.
            * units:
                - The units in which the transformation matrix is defined.
                  Either 'angstroms' or 'nanometers'
        '''
        super().__init__(1, '')
        self._process_transform(xyz_to_ijk_transform, units)
        self._initialize_transform_arguments(suffix)
        map_func = self._openmm_3D_function_from_volume(data)

        global_k_name = self._global_k_name = 'mdff_global_k_{}'.format(suffix)
        scale_factor_name = self._map_scale_factor_name = 'mdff_scale_factor_{}'.format(suffix)

        energy_func = self._set_energy_function(suffix)
        self.setEnergyFunction(energy_func)
        #super().__init__(1, energy_func)
        self._map_potential_index = self.addTabulatedFunction(
            name = 'map_potential', function = map_func)
        self._map_scale_factor_index = self.addGlobalParameter(
            name=scale_factor_name, defaultValue=1.0
        )
        self._global_k_index = self.addGlobalParameter(
            name = global_k_name, defaultValue = 1.0)
        self._individual_k_index = self.addPerBondParameter(
            name = 'individual_k')
        self._enabled_index = self.addPerBondParameter(
            name = 'enabled')
        self.update_needed = False

    @property
    def global_k_name(self):
        return self._global_k_name

    @property
    def map_scale_factor_name(self):
        return self._map_scale_factor_name

    def _process_transform(self, tf, units):
        if type(tf) == Quantity:
            rot_and_scale = tf[:,0:3].value_in_unit(OPENMM_LENGTH_UNIT)
            tf = tf.value_in_unit(tf.unit)
            tf[:,0:3] = rot_and_scale

        # TODO: Legacy code. Remove when possible.
        elif units == 'angstroms':
            # OpenMM calcs will be in nm
            # Only the 3x3 rotation/scale portion of the transformation
            # matrix changes with length units. The translation is applied
            # after scaling to the map coordinate system.
            tf[:,0:3] *= 10
        elif units != 'nanometers':
            raise TypeError('Units must be either "angstroms" or "nanometers"!')
        self._transform = tf

    def _initialize_transform_arguments(self, suffix):
        tf = self._transform
        tfi = self._tf_term_indices = numpy.zeros(tf.shape, int)
        for i in range(3):
            for j in range(3):
                tfi[i][j] = self.addGlobalParameter(
                    name = 'mdff_rot{}{}_{}'.format(i,j, suffix), defaultValue=tf[i,j]
                )
        for j in range(3):
            tfi[j,3] = self.addGlobalParameter(
                name = 'mdff_trn{}_{}'.format(j, suffix), defaultValue=tf[j,3]
            )


    def _openmm_3D_function_from_volume(self, data):
        raise RuntimeError('Cannot instantiate the base class!')

    def _set_energy_function(self, suffix):
        raise RuntimeError('Cannot instantiate the base class!')

    def set_global_k(self, k, context=None):
        '''
        Set the global coupling constant, in units of
        :math:`kJ mol^{-1} (\\text{map density unit})^{-1} nm^3`
        '''
        if context is not None:
            context.setParameter(self._global_k_name, k)
        else:
            self.setGlobalParameterDefaultValue(self._global_k_index, k)
            self.update_needed = True

    def set_map_scale_factor(self, scale_factor, context=None):
        '''
        Re-scale the map dimensions. Value of scale_factor should typically be
        very close to 1.0.
        '''
        if context is not None:
            context.setParameter(self._map_scale_factor_name, scale_factor)
        else:
            self.setGlobalParameterDefaultValue(self._map_scale_factor_index, scale_factor)
            self.update_needed = True

    def update_transform(self, transform, context=None, units='angstroms'):
        '''
        Update the complete map transform.
        '''
        self._process_transform(transform, units=units)
        tf = self._transform
        for i, row in enumerate(tf):
            for j, val in enumerate(row):
                index = self._tf_term_indices[i,j]
                if context is not None:
                    name = self.getGlobalParameterName(int(index))
                    context.setParameter(name, float(val))
                else:
                    self.setGlobalParameterDefaultValue(int(index), float(val))
                    self.update_needed=True

    # def update_transform(self, transform, units):
    #     self._process_transform(transform, units)
    #     tf = self._transform
    #     tfi = self._tf_term_indices
    #     for i in range(3):
    #         for j in range(4):
    #             self.setGlobalParameterDefaultValue(tfi[i][j], tf[i][j])
    #     self.update_needed = True

    def add_atoms(self, indices, ks, enableds):
        '''
        Add a set of atoms to the force, using a fast C++ function. Fastest if
        the inputs are NumPy arrays.

        Args:
            * indices:
                - An int array giving the indices of the atoms in the
                  simulation
            * ks:
                - A float array giving the per-atom scaling constants used
                  to determine the final potential for each atom (dimensionless)
            * enableds:
                - A Boolean array defining which atoms are to be enabled
                  in the force.
        '''
        n = len(indices)
        params = numpy.empty((n,2), numpy.float64)
        params[:,0] = ks
        params[:,1] = enableds
        return _openmm_force_ext.customcompoundbondforce_add_bonds(int(self.this), indices, params)

    def update_atoms(self, indices, ks, enableds):
        '''
        Update the parameters for a set of atoms at once, using a fast C++
        function. Fastest if the inputs are NumPy arrays.

        Args:
            * indices:
                - An int array giving the indices of the atoms in the
                  simulation
            * ks:
                - A float array giving the per-atom scaling constants used
                  to determine the final potential for each atom (dimensionless)
            * enableds:
                - A Boolean array defining which atoms are to be enabled
                  in the force.
        '''
        n = len(indices)
        params = numpy.empty((n,2), numpy.float64)
        params[:,0] = ks
        params[:,1] = enableds
        _openmm_force_ext.customcompoundbondforce_update_bond_parameters(int(self.this), indices, params)
        self.update_needed = True

    def update_atom(self, index, k, enabled):
        '''
        Update the parameters for a single atom in the force.

        Args:
            * index:
                - the integer index of the atom in the force object
            * k:
                - the per-atom scaling constant used to determine the final
                  potential for this atom
            * enabled:
                - a Boolean flag defining whether this atom can feel the MDFF
                  force
        '''
        self.update_atoms([index], [k], [enabled])

    def update_context_if_needed(self, context):
        '''
        If any parameters have changed since the last time this function was
        called, push them to the simulation context. Otherwise, do nothing.

        Args:
            * context:
                - the :class:`openmm.Context` to update.
        '''
        if self.update_needed:
            self.updateParametersInContext(context)
            self.update_needed = False

class CubicInterpMapForce_Old(_Map_Force_Base):
    '''
    (DEPRECATED. Uses the OpenMM :py:class:`Continuous3DFunction`, which pre-calculates
    64 spline parameters per map point. The memory cost of this becomes prohibitive 
    for large maps, and the spline calculation causes slow simulation startup times.
    The newer :py:class:`CubicInterpMapForce` simply uploads the map values to the GPU 
    as a discrete array, and performs the cubic interpolation on the fly.)

    Converts a map of (i,j,k) data and a (x,y,z)->(i,j,k) transformation
    matrix to a potential energy field, with tricubic interpolation of values.
    '''
    def __init__(self, data, xyz_to_ijk_transform, suffix, units = 'angstroms'):
        '''
        For a given atom at (x,y,z), the map potential will be defined
        as:

        .. math::
            global_k * individual_k * pot_xyz

        where `pot_xyz` is calculated by linear interpolation from the
        (i,j,k) map grid after applying `xyz_to_ijk_transform`.

        Args:
            * data:
                - The map data as a 3D (i,j,k) NumPy array in C-style order
            * xyz_to_ijk_transform:
                - A NumPy 3x4 float array defining the transformation matrix
                  mapping (x,y,z) coordinates to (i,j,k)
            * suffix:
                - In OpenMM, global parameters are global to the *entire
                  context*, not just to the force. To provide a global
                  parameter unique to this instance, the suffix is appended
                  to the base name of the parameter. Should be a unique string.
            * units:
                - The units in which the transformation matrix is defined.
                  Either 'angstroms' or 'nanometers'
        '''
        super().__init__(data, xyz_to_ijk_transform, suffix, units=units)

    def _openmm_3D_function_from_volume(self, data):
        dim = data.shape[::-1]
        data_1d = numpy.ravel(data, order = 'C')
        return Continuous3DFunction(*dim, data_1d, 0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1)

    def _set_energy_function(self, suffix):
        tf = self._transform
        # Transform xyz to ijk
        scale_str = ') / mdff_scale_factor_{}'.format(suffix)
        tf_strings = ['i = (', 'j = (', 'k = (']
        entries = ('x1 * {}','y1 * {}','z1 * {}','{}')
        for i in range(3):
            vals = tf[i]
            count = 0
            for j, (entry, val) in enumerate(zip(entries, vals)):
                if count > 0:
                    spacer = ' + '
                else:
                    spacer = ''
                if abs(val) > NEARLY_ZERO:
                    count += 1
                    if j<3:
                        tf_strings[i] += spacer + entry.format('mdff_rot{}{}_{}'.format(i,j, suffix))
                    else:
                        tf_strings[i] += spacer + entry.format('mdff_trn{}_{}'.format(i, suffix) + scale_str)
                elif j==3:
                    tf_strings[i] += scale_str

        funcs = ';'.join(tf_strings)
        enabled_eqn = 'step(enabled-0.5)'
        energy_str = '-{} * individual_k * map_potential(i,j,k)'.format(self._global_k_name)

        final_func = 'select({}, {}, 0); {}'.format(enabled_eqn, energy_str, funcs)
        return final_func

class CubicInterpMapForce(_Map_Force_Base):
    '''
    Creates a MDFF potential from a 3D map of density values.
    '''
    def __init__(self, data, xyz_to_ijk_transform, suffix, units = 'angstroms'):
        '''
        For a given atom at (x,y,z), the map potential will be defined
        as:

        .. math::
            global_k * individual_k * pot_xyz

        where `pot_xyz` is calculated by tricubic interpolation from the
        (i,j,k) map grid after applying `xyz_to_ijk_transform`.

        Args:
            * data:
                - The map data as a 3D (i,j,k) NumPy array in C-style order
            * xyz_to_ijk_transform:
                - A NumPy 3x4 float array defining the transformation matrix
                  mapping (x,y,z) coordinates to (i,j,k)
            * suffix:
                - In OpenMM, global parameters are global to the *entire
                  context*, not just to the force. To provide a global
                  parameter unique to this instance, the suffix is appended
                  to the base name of the parameter. Should be a unique string.
            * units:
                - The units in which the transformation matrix is defined.
                  Either 'angstroms' or 'nanometers'
        '''
        super().__init__(data, xyz_to_ijk_transform, suffix, units=units)

    def _set_energy_function(self, suffix):
        tf = self._transform
        # Transform xyz to ijk
        scale_str = ') / mdff_scale_factor_{}'.format(suffix)
        tf_strings = ['i = (', 'j = (', 'k = (']
        entries = ('x1 * {}','y1 * {}','z1 * {}','{}')
        for i in range(3):
            vals = tf[i]
            count = 0
            for j, (entry, val) in enumerate(zip(entries, vals)):
                if count > 0:
                    spacer = ' + '
                else:
                    spacer = ''
                if abs(val) > NEARLY_ZERO:
                    count += 1
                    if j<3:
                        tf_strings[i] += spacer + entry.format('mdff_rot{}{}_{}'.format(i,j, suffix))
                    else:
                        tf_strings[i] += spacer + entry.format('mdff_trn{}_{}'.format(i, suffix) + scale_str)
                elif j==3:
                    tf_strings[i] += scale_str

        i_str, j_str, k_str = tf_strings

        min_str = 'min_i = floor(i-1); min_j = floor(j-1); min_k = floor(k-1)'

        offset_strings = (
                          'ci0 = 1.0-ci1;'
                          'cj0 = 1.0-cj1;'
                          'ck0 = 1.0-ck1;'
                          'ci1 = i-min_i-1;'
                          'cj1 = j-min_j-1;'
                          'ck1 = k-min_k-1'
                          )

        coeff_strings = ('sc_i0 = -0.5*ci1*ci0*ci0;'
                         'sc_i1 = ci0*( -1.5*ci1*ci1 + ci1 + 1.0);'
                         'sc_i2 = ci1*( -1.5*ci0*ci0 + ci0 + 1.0);'
                         'sc_i3 = -0.5*ci1*ci1*ci0;'
                         'sc_j0 = -0.5*cj1*cj0*cj0;'
                         'sc_j1 = cj0*( -1.5*cj1*cj1 + cj1 + 1.0);'
                         'sc_j2 = cj1*( -1.5*cj0*cj0 + cj0 + 1.0);'
                         'sc_j3 = -0.5*cj1*cj1*cj0;'
                         'sc_k0 = -0.5*ck1*ck0*ck0;'
                         'sc_k1 = ck0*( -1.5*ck1*ck1 + ck1 + 1.0);'
                         'sc_k2 = ck1*( -1.5*ck0*ck0 + ck0 + 1.0);'
                         'sc_k3 = -0.5*ck1*ck1*ck0'
                         )

        sum_strings = []
        for i in range(4):
            sv_strings = []
            for j in range(4):
                str_1 = '(sc_k0 * map_potential(min_i+{}, min_j+{}, min_k))'.format(i,j)
                str_2 = '(sc_k1 * map_potential(min_i+{}, min_j+{}, min_k+1))'.format(i,j)
                str_3 = '(sc_k2 * map_potential(min_i+{}, min_j+{}, min_k+2))'.format(i,j)
                str_4 = '(sc_k3 * map_potential(min_i+{}, min_j+{}, min_k+3))'.format(i,j)
                sv_strings.append('sc_j{} * ({}+{}+{}+{})'.format(
                    j, str_1, str_2, str_3, str_4
                ))
            sum_strings.append('sc_i{} * ({})'.format(
                i, '+'.join(sv_strings)
            ))
        master_sum_str = '+'.join(sum_strings)

        energy_str = '-{} * individual_k * ({})'.format(self._global_k_name, master_sum_str)

        enabled_eqn = 'step(enabled-0.5)'

        funcs = ';'.join((coeff_strings, offset_strings, min_str, i_str, j_str, k_str))
        final_func = 'select({}, {}, 0); {}'.format(enabled_eqn, energy_str, funcs)
        return final_func

    def _openmm_3D_function_from_volume(self, data):
        dim = data.shape[::-1]
        data_1d = numpy.ravel(data, order = 'C')
        return Discrete3DFunction(*dim, data_1d)


class LinearInterpMapForce(_Map_Force_Base):
    '''
    NOTE: This class is deprecated, since there is no conceivable situation in which
    it is superior to either :class:`CubicInterpMapForce` or
    :class:`CubicInterpMapForce_Old`.

    Converts a map of (i,j,k) data and a (x,y,z)->(i,j,k) transformation
    matrix to a potential energy field, with trilinear interpolation of values.
    '''
    def __init__(self, data, xyz_to_ijk_transform, suffix, units = 'angstroms'):
        '''
        For a given atom at (x,y,z), the map potential will be defined
        as:

        .. math::
            global_k * individual_k * pot_xyz

        where `pot_xyz` is calculated by linear interpolation from the
        (i,j,k) map grid after applying `xyz_to_ijk_transform`.

        Args:
            * data:
                - The map data as a 3D (i,j,k) NumPy array in C-style order
            * xyz_to_ijk_transform:
                - A NumPy 3x4 float array defining the transformation matrix
                  mapping (x,y,z) coordinates to (i,j,k)
            * suffix:
                - In OpenMM, global parameters are global to the *entire
                  context*, not just to the force. To provide a global
                  parameter unique to this instance, the suffix is appended
                  to the base name of the parameter. Should be a unique string.
            * units:
                - The units in which the transformation matrix is defined.
                  Either 'angstroms' or 'nanometers'
        '''
        super().__init__(data, xyz_to_ijk_transform, suffix, units=units)

    def _set_energy_function(self, suffix):
        tf = self._transform
        # Transform xyz to ijk
        scale_str = ') / mdff_scale_factor_{}'.format(suffix)
        tf_strings = ['i = (', 'j = (', 'k = (']
        entries = ('x1 * {}','y1 * {}','z1 * {}','{}')
        for i in range(3):
            vals = tf[i]
            count = 0
            for j, (entry, val) in enumerate(zip(entries, vals)):
                if count > 0:
                    spacer = ' + '
                else:
                    spacer = ''
                if abs(val) > NEARLY_ZERO:
                    count += 1
                    if j<3:
                        tf_strings[i] += spacer + entry.format('mdff_rot{}{}_{}'.format(i,j, suffix))
                    else:
                        tf_strings[i] += spacer + entry.format('mdff_trn{}_{}'.format(i, suffix) + scale_str )
                elif j==3:
                    tf_strings[i] += scale_str

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

        energy_str = '-{} * individual_k * {}'.format(self._global_k_name, interp_str)

        enabled_eqn = 'step(enabled-0.5)'

        funcs = ';'.join(( norm_str, val_str, min_str, max_str,
                        i_str, j_str, k_str))
        final_func = 'select({}, {}, 0); {}'.format(enabled_eqn, energy_str, funcs)

        return final_func

    def _openmm_3D_function_from_volume(self, data):
        dim = data.shape[::-1]
        data_1d = numpy.ravel(data, order = 'C')
        return Discrete3DFunction(*dim, data_1d)

class AdaptiveDistanceRestraintForce(CustomBondForce):
    r'''
    A :py:class:`openmm.CustomBondForce` subclass using the generalised adaptive
    loss function described by Jonathan Barron
    (https://arxiv.org/pdf/1701.03077.pdf).

    There are many situations in which collections of distance restraints may
    be approximate, or may contain some subset of restraints that are simply
    incorrect. A simple example is a distance restraint network created based on
    a reference model of your protein (e.g. a non-crystallographic symmetry copy
    or an independent experimental model of the same protein or close
    homologue). In this case the majority of restraints are likely to be
    correct, but any bulk conformational change between reference model and
    target will lead to some subset being wildly incorrect. In such cases it is
    advantageous to use an energy function which flattens out (that is, yields
    low-to-zero net force) not only at the target distance, but also once the
    distance deviates sufficiently from the target. In other words, restraints
    that are close to satisfied will be enforced, but restraints that are
    incompatible with the data will gracefully release.

    Various loss functions meeting these criteria have been proposed and used in
    the past, each with its own specific characteristics. The function applied
    here is a generalisation of the Cauchy/Lorentzian, Geman-McClure,
    Welsch/Leclerc, generalised Charbonnier, Charbonnier/pseudo-Huber/L1-L2 and
    L2 loss functions, using a single "robustness" parameter to tune the rate of
    fall-off of restraint force with distance. This has been further adapted to
    optionally allow for an arbitrary amount of "lee-way" - that is, a range in
    which no bias is applied - either side of the target distance. The complete
    functional form is:


    .. math::
        E = \kappa *
        \begin{cases}
            0, & \text{if}\ enabled < 0.5 \text{ or}\ |r-r_0| < \tau \\
            1/2 (\frac{r-\rho}{c})^2, & \text{if}\ \alpha = 2 \\
            ln(\frac{1}{2} (\frac{r-\rho}{c})^2 + 1), & \text{if}\ \alpha = 0 \\
            \frac{|2-\alpha|}{\alpha} ((\frac{ (\frac{r-\rho}{c})^2 }{|2-\alpha|} + 1)^\frac{\alpha}{2} - 1), & \text{otherwise}
        \end{cases}

    where

    .. math::
        \rho =
        \begin{cases}
            r-\tau, & \text{if}\ (r-r_0) < -\tau \\
            r+\tau, & \text{if}\ (r-r_0) > \tau
        \end{cases}

    :attr:`c` sets the width of the energy well and the distance at which the
    function switches from quadratic.

    :math:`\kappa` adjusts the depth of the well (that is, the absolute
    magnitude of the applied force). Within the central "well" of the function,
    :math:`\frac{\kappa}{c^2}` is equivalent to the spring constant in a
    standard harmonic restraint.

    :math:`\tau` sets the tolerance around the target distance. If
    :math:`|r-r_0| < \tau` no restraining force will be applied.

    :math:`\alpha` is the master parameter setting the "robustness" of the
    restraint. For all values of :math:`\alpha`, the shape of the function within
    the range :math:`-c < r-\rho < c` is essentially parabolic. Outside of this
    region, increasing positive values increase the steepness of the "walls"
    of the restraint (NOTE: since this force does not currently obey the
    :attr:`max_force` imposed on ISOLDE's other custom forces, use this with
    care to avoid instability). Values larger than 2 are inadvisable. A value
    of :math:`\alpha = 2` yields a standard harmonic (quadratic) restraint.
    When :math:`\alpha = 1`, the restraining force transitions from harmonic to
    essentially constant (i.e. linear) between :math:`c < |r-\rho| < 2c`. When
    :math:`\alpha = 0`, the force at large distances is proportional to
    :math:`\frac{1}{|r-\rho|}`. For increasing negative values of :math:`\alpha`
    the applied force falls off faster with distance. When :math:`\alpha = -2`,
    the force falls to 10% of the maximum at approx. :math:`|r-\rho| = 7c`; when
    :math:`\alpha = -10` the force at this point is approx. 0.1% of maximum. For
    very large negative values of :math:`alpha` (beyond about -50) the
    function converges to a form where the applied force is negligible for
    :math:`|r-rho| > 6c`.

    All parameters are settable at the individual restraint level.

    '''
    def __init__(self):

        # Loss functions
        alpha_2_eqn = '1/2 * delta_r_on_c_sq'
        alpha_0_eqn = 'log(1/2*delta_r_on_c_sq + 1)'
        general_eqn = 'abs(2-alpha)/alpha * ((delta_r_on_c_sq/abs(2-alpha) + 1)^(alpha/2) - 1)'


        # Switching functions
        alpha_2_sw = 'delta(alpha-2)'
        alpha_0_sw = 'delta(alpha)'
        enabled_sw = 'step(enabled-0.5)'
        tol_sw = 'step(abs(delta_r)-tau)'

        # Intermediate variable definitions
        delta_r_on_c_sq_def = 'delta_r_on_c_sq = ((r-rho)/c)^2'
        rho_def = 'rho = select(delta_r-tau, r0+tau, r0-tau)'
        delta_r_def = 'delta_r = r-r0'

        energy_fn = 'energy = kappa * select({}, {}, select({}, {}, {}))'.format(
            alpha_2_sw, alpha_2_eqn, alpha_0_sw, alpha_0_eqn, general_eqn
        )

        final_str = 'select({},select({}, energy, 0), 0);{};{};{};{}'.format(
            enabled_sw, tol_sw,
            energy_fn, delta_r_on_c_sq_def, rho_def, delta_r_def
        )
        super().__init__(final_str)
        self.enabled_index = self.addPerBondParameter('enabled')
        self.kappa_index = self.addPerBondParameter('kappa')
        self.c_index = self.addPerBondParameter('c')
        self.r0_index = self.addPerBondParameter('r0')
        self.tau_index = self.addPerBondParameter('tau')
        self.alpha_index = self.addPerBondParameter('alpha')

        self.update_needed = False

    def add_bonds(self, atom_indices, enableds, kappas, cs, targets, tolerances, alphas):
        r'''
        Add a set of bonds to the simulation, using a fast C++ function. Fastest
        if all parameters are supplied as NumPy arrays.

        Args:

            * atom_indices:
                - a 2-tuple of integer arrays giving the indices of the bonded
                  atoms in the simulation construct
            * enableds:
                - a Boolean array defining which restraints are to be active
            * kappas:
                - a float array of energy scaling constants in
                  :math:`kJ mol^{-1}`. For a given restraint,
                  :math:`\frac{\kappa}{c^2}` is equivalent to a harmonic spring
                  constant when the distance is close to the target.
            * cs:
                - a float array setting the "width" of the energy well for each
                  restraint in nanometres. A restraint behaves like a normal
                  harmonic restraint when the current distance is less than
                  :attr:`c` from the target distance.
            * targets:
                - a float array of target distances in nanometres
            * tolerances:
                - a float array in nanometres. When :math:`|r-r_0|` is less than
                  the tolerance, no force will be applied.
            * alphas:
                - a float array of terms setting the "steepness" of each
                  restraint outside of the central well. Values less than one
                  cause the applied force to fall off with increasing distance.
        '''
        n = len(targets)
        ind = numpy.empty((n,2), numpy.int32)
        for i, ai in enumerate(atom_indices):
            ind[:,i] = ai
        params = numpy.empty((n,6), numpy.float64)
        params[:,0] = enableds
        params[:,1] = kappas
        params[:,2] = cs
        params[:,3] = targets
        params[:,4] = tolerances
        params[:,5] = alphas
        return _openmm_force_ext.custombondforce_add_bonds(int(self.this), ind, params)

    def update_target(self, index, enabled=None, kappa=None, c=None,
            target=None, tolerance=None, alpha=None):
        r'''
        Update the parameters for an existing restraint in the simulation.
        Mostly superseded by :func:`update_targets`.

        Args:

            * index:
                - the index of this restraint in the OpenMM force object
            * enabled:
                - Boolean flag defining whether the restraint is to be enabled.
                  None = keep current value.
            * kappa:
                - Energy scaling constant (as a :class:`simtk.Quantity` or in
                  units of :math:`kJ mol^{-1}`). When the distance is close to
                  the target, :math:`\frac{\kappa}{c^2}` is equivalent to a
                  harmonic spring constant. None = keep current value.
            * c:
                - A distance (in nanometres) defining the width of the central
                  "well" of the energy function. None = keep current value.
            * target:
                - the new target distance (as a :class:`simtk.Quantity` or in
                  nanometres). None = keep current value.
            * tolerance:
                - Distance (in nanometres) around the target distance below
                  which no biasing force will be applied. None = keep current
                  value.
            * alpha:
                - Dimensionless value dictating how quickly the energy grows or
                  levels out for large deviations from the target distance.
                  Values less than one cause the applied force to fall off with
                  increasing distance. None = keep current value.
        '''
        current_params = self.getBondParameters(int(index))
        atom1, atom2 = current_params[0:2]
        new_params = list(current_params[2])
        for i, p in enumerate((enabled, kappa, c, target, tolerance, alpha)):
            if p is not None:
                new_params[i] = float(_strip_units(p))
        self.setBondParameters(int(index), atom1, atom2, new_params)
        self.update_needed = True

    def update_targets(self, indices, enableds, kappas, cs, targets, tolerances, alphas):
        r'''
        Update a set of targets all at once using fast C++ code. Fastest if
        the arguments are provided as Numpy arrays, but any iterable will work.

        Args:

            * indices:
                - the indices of the restraints in the OpenMM force object
            * enableds:
                - a Boolean array defining which restraints are to be active
            * kappas:
                - a float array of energy scaling constants in
                  :math:`kJ mol^{-1}`. For a given restraint,
                  :math:`\frac{\kappa}{c^2}` is equivalent to a harmonic spring
                  constant when the distance is close to the target.
            * cs:
                - a float array setting the "width" of the energy well for each
                  restraint in nanometres. A restraint behaves like a normal
                  harmonic restraint when the current distance is less than
                  :attr:`c` from the target distance.
            * targets:
                - a float array of target distances in nanometres
            * tolerances:
                - a float array in nanometres. When :math:`|r-r_0|` is less than
                  the tolerance, no force will be applied.
            * alphas:
                - a float array of terms setting the "steepness" of each
                  restraint outside of the central well. Values less than one
                  cause the applied force to fall off with increasing distance.
        '''
        n = len(indices)
        params = numpy.empty((n,6), numpy.float64)
        params[:,0] = enableds
        params[:,1] = kappas
        params[:,2] = cs
        params[:,3] = targets
        params[:,4] = tolerances
        params[:,5] = alphas
        _openmm_force_ext.custombondforce_update_bond_parameters(int(self.this), indices, params)
        self.update_needed = True



class TopOutBondForce(CustomBondForce):
    r'''
    A :py:class:`openmm.CustomBondForce` subclass defined as a standard harmonic
    potential with a user-defined fixed maximum cutoff on the applied force. Any
    restraint can be switched on (off) by setting the 'enabled' parameter to 1
    (0). This is designed for steering the simulation into new conformations
    where the starting distance may be far from the target bond length, leading
    to catastrophically large forces with a standard harmonic potential. The
    effective energy equation is:

    .. math::

        E =
        \begin{cases}
            0, & \text{if}\ enabled < 0.5 \\
            \text{max\_force} * abs(r-r_0), & \text{if}\ (r-r_0) - \text{max\_force}/k > 0 \\
            0.5 * k * (r - r_0)^2, & \text{otherwise}
        \end{cases}

    '''
    def __init__(self, max_force):
        '''
        Initialise the force object and set the maximum force magnitude.

        Args:
            * max_force:
                - maximum allowable force in :math:`kJ mol^{-1} nm^{-1}`
        '''
        linear_eqn = 'max_force * abs(r-r0)'# - 0.5*max_force^2/k'
        quadratic_eqn = '0.5*k*(r-r0)^2'
        transition_eqn = 'step(r - r0 - max_force/k)'
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
        '''
        Get/set the maximum force to be applied to any given atom, in
        :math:`kJ mol^{-1} nm^{-1}`
        '''
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
        Add a set of bonds to the simulation, using a fast C++ function. Fastest
        if all parameters are supplied as NumPy arrays.

        Args:
            * atom_indices:
                - a 2-tuple of integer arrays giving the indices of the bonded
                  atoms in the simulation construct
            * enableds:
                - a Boolean array defining which restraints are to be active
            * spring_constants:
                - a float array of spring constants in :math:`kJ mol^{-1} nm^{-2}`
            * targets:
                - a float array of target distances in nanometers
        '''
        n = len(targets)
        ind = numpy.empty((n,2), numpy.int32)
        for i, ai in enumerate(atom_indices):
            ind[:,i] = ai
        params = numpy.empty((n,3), numpy.float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2] = targets
        return _openmm_force_ext.custombondforce_add_bonds(int(self.this), ind, params)


    def update_target(self, index, enabled=None, k = None, target=None):
        '''
        Update the parameters for an existing restraint in the simulation.
        Mostly superseded by :func:`update_targets`.

        Args:
            * index:
                - the index of this restraint in the OpenMM force object
            * enabled:
                - Boolean flag defining whether the restraint is to be enabled.
                  None = keep current value.
            * k:
                - The new spring constant (as a :class:`simtk.Quantity` or in
                  units of :math:`kJ mol^{-1} nm^{-2}`). None = keep current
                  value.
            * target:
                - the new target distance (as a :class:`simtk.Quantity` or in
                  nanometres). None = keep current value.
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

        Args:
            * indices:
                - the indices of the restraints in the OpenMM force object
            * enableds:
                - a Boolean array defining which restraints are to be enabled
            * spring_constants:
                - The new spring constants in units of :math:`kJ mol^{-1} nm^{-2}`
            * targets:
                - the new target distances in nanometres
        '''
        n = len(indices)
        params = numpy.empty((n,3), numpy.float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2] = targets
        _openmm_force_ext.custombondforce_update_bond_parameters(int(self.this), indices, params)
        self.update_needed = True


class TopOutRestraintForce(CustomExternalForce):
    r'''
    A :py:class:`openmm.CustomExternalForce` subclass designed to restrain atoms
    to defined positions via a standard harmonic potential with a user-defined
    fixed maximum cutoff on the applied force. Used for position restraints as
    well as for imposing interactive tugging forces, this is designed for
    steering the simulation into new conformations where the starting positions
    may be far from the target positions, leading to catastrophically large
    forces with a standard harmonic potential. The effective energy equation is:

    .. math::

        E =
        \begin{cases}
            0, & \text{if}\ enabled < 0.5 \\
            \text{max\_force} * r, & \text{if}\ r - \text{max\_force}/k > 0 \\
            0.5 * k * r^2, & \text{otherwise}
        \end{cases}

    '''
    def __init__(self, max_force):
        '''
        Initialise the force object and set the maximum force magnitude.

        Args:
            * max_force:
                - maximum allowable force in :math:`kJ mol^{-1} nm^{-1}`
        '''
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
        '''
        Get/set the maximum force applied to any given atom, in
        :math:`kJ mol^{-1} nm^{-1}`.
        '''
        return self._max_force

    @max_force.setter
    def max_force(self, force):
        if type(force) == Quantity:
            force = force.value_in_unit(OPENMM_FORCE_UNIT)
        self.setGlobalParameterDefaultValue(0, force)
        self._max_force = force
        self.update_needed = True

    def add_particles(self, indices, enableds, spring_constants, targets):
        '''
        Add a set of restraints to the simulation, using a fast C++ function.
        Fastest if all parameters are supplied as NumPy arrays.

        Args:
            * atom_indices:
                - integer array giving the indices of the restrained atoms in
                  the simulation construct
            * enableds:
                - a Boolean array defining which restraints are to be active
            * spring_constants:
                - a float array of spring constants in :math:`kJ mol^{-1} nm^{-2}`
            * targets:
                - a (nx3) float array of (x,y,z) target positions in nanometres
        '''
        n = len(indices)
        params = numpy.empty((n,5), numpy.float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2:] = targets
        return _openmm_force_ext.customexternalforce_add_particles(int(self.this), indices, params)

    def update_target(self, index, enabled=None, k=None, target=None):
        '''
        Update a single restraint. This function is mostly superseded by
        :func:`update_targets`.

        Args:
            * index:
                - integer index of the restraint in the force object
            * enabled:
                - enable/disable the restraint. None keeps the current value.
            * k:
                - set the spring constant in :math:`kJ mol^{-1} nm^{-2}`.
                  None keeps the current value.
            * target:
                - set the target (x,y,z) position in nanometres. None keeps
                  the current value.
        '''
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
        '''
        Update a set of targets all at once using fast C++ code. Fastest if
        the arguments are provided as Numpy arrays, but any iterable should work.

        Args:
            * indices:
                - the indices of the restraints in the OpenMM force object
            * enableds:
                - a Boolean array defining which restraints are to be enabled
            * spring_constants:
                - The new spring constants in units of :math:`kJ mol^{-1} nm^{-2}`
            * targets:
                - A (nx3) float array providing the new target (x,y,z) positions
                  in nanometres.
        '''
        n = len(indices)
        if len(enableds)!=n or len(targets) !=n or len(spring_constants) !=n:
            raise TypeError('Parameter array lengths must match number of indices!')
        params = numpy.empty((n,5), numpy.float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2:] = targets
        _openmm_force_ext.customexternalforce_update_particle_parameters(int(self.this), indices, params)
        self.update_needed = True

    def release_restraint(self, index):
        '''
        Disable a single restraint.

        Args:
            * index:
                - the index of the restraint to be disabled.
        '''
        self.update_target(index, enabled=False)


class FlatBottomTorsionRestraintForce(CustomTorsionForce):
    r'''
    A :py:class:`openmm.CustomTorsionForce` subclass designed to restrain
    torsion angles while allowing free movement within a range (target +/-
    cutoff). Within the cutoff range the potential is constant (that is, zero
    force is applied).

    The effective energy function is:

    .. math::

        E =
        \begin{cases}
            0, & \text{if}\ enabled < 0.5 \\
            -k*cos(\theta_\text{cutoff}), & \text{if}\ cos(\theta-\theta_0) - cos(\theta_\text{cutoff}) < 0 \\
            -k*cos(\theta-\theta_0), & \text{otherwise}
        \end{cases}
    '''
    def __init__(self):
        '''
        Initialise the force object. No restraints are added at this stage.
        '''
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
        Add a set of torsion restraints using a fast C++ function. Returns a
        NumPy integer array giving the indices of the restraints in the force
        object. Fastest if the inputs are NumPy arrays, but most iterables
        should work.

        Args:
            * atom_indices:
                - A 4-tuple of arrays providing the indices of the dihedral
                  atoms in the simulation construct
            * enableds:
                - A Boolean array (or any array castable to float) where values
                  > 0.5 represent enabled restraints
            * spring_constants:
                - Restraint spring constants in :math:`kJ mol^{-1} rad^{-2}`
            * targets:
                - Target angles in radians
            * cutoffs:
                - Cutoff angle (below which no force is applied) for each
                  restraint in radians.
        '''
        n = len(targets)
        ind = numpy.empty((n,4), numpy.int32)
        for i, ai in enumerate(atom_indices):
            ind[:,i] = ai
        params = numpy.empty((n,4), numpy.float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2] = targets
        params[:,3] = numpy.cos(cutoffs)
        return _openmm_force_ext.customtorsionforce_add_torsions(int(self.this), ind, params)

    def update_target(self, index, enabled=None, k = None, target = None, cutoff = None):
        '''
        Update a single restraint. This function is mostly superseded by
        :func:`update_targets`.

        Args:
            * index:
                - integer index of the restraint in the force object
            * enabled:
                - enable/disable the restraint. None keeps the current value.
            * k:
                - set the spring constant in :math:`kJ mol^{-1} rad^{-2}`.
                  None keeps the current value.
            * target:
                - set the target angle in radians. None keeps the current value.
            * cutoff:
                - set the cutoff angle in radians. None keeps the current value.
        '''
        # For compatibility with numpy.int32
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
        Update a set of targets all at once using fast C++ code. Fastest if
        the arguments are provided as NumPy arrays, but any iterable should work.

        Args:
            * indices:
                - the indices of the restraints in the OpenMM force object
            * enableds:
                - a Boolean array defining which restraints are to be enabled
            * spring_constants:
                - the new spring constants in units of :math:`kJ mol^{-1} rad^{-2}`
            * targets:
                - the new target angles in radians
            * cutoffs:
                - the new cutoff angles in radians
        '''
        n = len(indices)
        params = numpy.empty((n,4), numpy.float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2] = targets
        params[:,3] = numpy.cos(cutoffs)
        _openmm_force_ext.customtorsionforce_update_torsion_parameters(int(self.this), indices, params)
        self.update_needed = True

class TopOutTorsionForce(CustomTorsionForce):
    r'''
    Torsion-space analogy to the AdaptiveDistanceRestraintForce: often when
    restraining the model to a template (or restraining NCS copies to their
    average) we want to try to ensure that torsions that *truly* differ
    substantially from the target aren't penalised.

    The functional form used here is somewhat analogous to the von Mises
    distribution (an approximation to the spherically-wrapped normal
    distribution), but shouldn't really be considered as a probability
    distribution. Nor should it be really used as an energy potential in
    situations outside of the fairly narrow scope of restraining to templates,
    since it has little physical meaning. It is normalised such that the maximum
    value of the first derivative (i.e. the maximum applied force) is
    independent of :math:`\kappa`, the term setting the width of the "well" in
    which a substantial restraining force will be applied.

    The effective energy function is:

    .. math::

        E =
        \begin{cases}
            0, & \text{if}\ enabled < 0.5 \\
            1-k\frac{ \sqrt{2} e^{\frac{-1}{2}\sqrt{4\kappa^2+1}-\kappa+\frac{1}{2}}
            e^{\kappa(\cos{(\theta-\theta_0)}+1)-1)}}
            {\sqrt{\sqrt{4\kappa^2+1}-1}}, & \text{if}\ \kappa>0 \\
            -k\cos{(\theta-\theta_0)}, & \text{if}\ \kappa=0
        \end{cases}

    For values of :math:`\kappa` greater than about 1, :math:`\frac{1}{\kappa}`
    is approximately equal to the variance of a periodic normal distribution
    centred on :math:`\theta-\theta_0`. As :math:`\kappa` approaches zero, the
    energy term converges to a standard unimodal cosine.
    '''
    def __init__(self):
        default_energy_term = ('-1 + '
                      '-k * sqrt(2)*exp(-1/2*sqrt(4*kappa^2+1)-kappa+1/2)'
                        '* exp(kappa*(cos(theta-theta0)+1)-1)'
                        '/ sqrt(sqrt(4*kappa^2+1)-1)')
        limiting_energy_term = '-k*cos(theta-theta0)'
        energy_switch_term = 'delta(kappa)'
        full_energy_term = 'select({}, {}, {})'.format(
            energy_switch_term, limiting_energy_term, default_energy_term)
        enabled_function = 'step(enabled-0.5)'
        complete_function = 'select({}, {}, 0)'.format(
            enabled_function, full_energy_term)

        super().__init__(complete_function)
        per_bond_parameters = ('enabled', 'k', 'theta0', 'kappa')
        for p in per_bond_parameters:
            self.addPerTorsionParameter(p)

        self.update_needed = False

    def add_torsions(self, atom_indices, enableds, spring_constants, targets, kappas):
        r'''
        Add a set of torsion restraints using a fast C++ function. Returns a
        NumPy integer array giving the indices of the restraints in the force
        object. Fastest if the inputs are NumPy arrays, but most iterables
        should work.

        Args:
            * atom_indices:
                - A 4-tuple of arrays providing the indices of the dihedral
                  atoms in the simulation construct
            * enableds:
                - A Boolean array (or any array castable to float) where values
                  > 0.5 represent enabled restraints
            * spring_constants:
                - Restraint spring constants in :math:`kJ mol^{-1} rad^{-2}`
            * targets:
                - Target angles in radians
            * kappas:
                - Constants defining the width of the restrained well. For
                  :math:`\kappa > 0.5`, :math:`1/\kappa` is approximately equal
                  to the variance of a normal distribution centred on
                  :math:`\theta_0`. The energy approaches
                  :math:`-k*cos(\theta-\theta_0)` as :math:`\kappa` approaches
                  zero. Values less than 0 are not allowed, and will be
                  automatically set to 0.
        '''
        n = len(targets)
        ind = numpy.empty((n,4), numpy.int32)
        for i, ai in enumerate(atom_indices):
            ind[:,i] = ai
        params = numpy.empty((n,4), numpy.float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2] = targets
        params[:,3] = kappas
        params[:,3][params[:,3]<0] = 0
        return _openmm_force_ext.customtorsionforce_add_torsions(int(self.this), ind, params)

    def update_target(self, index, enabled=None, k = None, target = None, kappa = None):
        '''
        Update a single restraint. This function is mostly superseded by
        :func:`update_targets`.

        Args:
            * index:
                - integer index of the restraint in the force object
            * enabled:
                - enable/disable the restraint. None keeps the current value.
            * k:
                - set the spring constant in :math:`kJ mol^{-1} rad^{-2}`.
                  None keeps the current value.
            * target:
                - set the target angle in radians. None keeps the current value.
            * kappa:
                - set the kappa (approximate units of inverse square radians)
        '''
        # For compatibility with numpy.int32
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
        if kappa is not None:
            if type(kappa) == Quantity:
                kappa = kappa.value_in_unit(OPENMM_ANGLE_UNIT)
            new_kappa = kappa
        self.setTorsionParameters(index, *indices, (new_enabled, new_k, new_theta0, new_kappa))
        self.update_needed = True

    def update_targets(self, indices, enableds, spring_constants, targets, kappas):
        '''
        Update a set of targets all at once using fast C++ code. Fastest if
        the arguments are provided as NumPy arrays, but any iterable should work.

        Args:
            * indices:
                - the indices of the restraints in the OpenMM force object
            * enableds:
                - a Boolean array defining which restraints are to be enabled
            * spring_constants:
                - the new spring constants in units of :math:`kJ mol^{-1} rad^{-2}`
            * targets:
                - the new target angles in radians
            * cutoffs:
                - the new kappas in inverse square radians
        '''
        n = len(indices)
        params = numpy.empty((n,4), numpy.float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2] = targets
        params[:,3] = kappas
        _openmm_force_ext.customtorsionforce_update_torsion_parameters(int(self.this), indices, params)
        self.update_needed = True

class AdaptiveTorsionForce(CustomTorsionForce):
    r'''
    Torsion-space analogy to the AdaptiveDistanceRestraintForce: often when
    restraining the model to a template (or restraining NCS copies to their
    average) we want to try to ensure that torsions that *truly* differ
    substantially from the target aren't penalised.

    The functional form of the core energy term used here is the same as that in
    :class:`TopOutTorsionForce`. This is then further modified to add an
    adjustable fall-off outside of the central "well" - that is, where
    :class:`TopOutTorsionForce` settles quickly to a flat potential,
    :class:`AdaptiveTorsionForce` has an adjustable slope in the same region.

    The effective energy function is:

    .. math::

        E_core =
        \begin{cases}
            0, & \text{if}\ enabled < 0.5 \\
            1-\frac{ \sqrt{2} e^{\frac{-1}{2}\sqrt{4\kappa^2+1}-\kappa+\frac{1}{2}}
            e^{\kappa(\cos{(\theta-\theta_0)}+1)-1)}}
            {\sqrt{\sqrt{4\kappa^2+1}-1}}, & \text{if}\ \kappa>0 \\
            -\cos{(\theta-\theta_0)}, & \text{if}\ \kappa=0
        \end{cases}

        E_final = k* (E_core + \alpha*(e^{\sqrt{\alpha}*(E_core-1)}*(1-\cos{\theta-\theta_0}) )

    For values of :math:`\kappa` greater than about 1, :math:`\frac{1}{\kappa}`
    is approximately equal to the variance of a periodic normal distribution
    centred on :math:`\theta-\theta_0`. As :math:`\kappa` approaches zero, the
    energy term converges to a standard unimodal cosine if :math:`\alpha` is
    zero. For most common uses :math:`alpha` should be a value between zero and
    1. Values of :math:`alpha` greater than 1 cause the potential to become
    steeper outside the central well than inside; values less than zero cause
    the energy to *drop* outside the well (that is, torsions that are not inside
    the well will be repelled from it)
    '''
    def __init__(self):
        default_energy_term = ('1 - '
            'sqrt(2)*exp(-1/2*sqrt(4*kappa^2+1)-kappa+1/2)'
            '* exp(kappa*(cos(theta-theta0)+1)-1)'
            '/ sqrt(sqrt(4*kappa^2+1)-1)')
        limiting_energy_term = '-cos(theta-theta0)'
        energy_switch_term = 'delta(kappa)'
        base_energy_term = 'select({}, {}, {})'.format(
            energy_switch_term, limiting_energy_term, default_energy_term)
        full_energy_term = (
            'k * (base_energy + alpha*exp(sqrt(alpha)*(base_energy-1))'
            '*(1-cos(theta-theta0)))').format(base_energy_term)
        enabled_function = 'step(enabled-0.5)'
        complete_function = 'select({}, {}, 0); base_energy={}'.format(
            enabled_function, full_energy_term, base_energy_term)

        super().__init__(complete_function)
        per_bond_parameters = ('enabled', 'k', 'theta0', 'kappa', 'alpha')
        for p in per_bond_parameters:
            self.addPerTorsionParameter(p)

        self.update_needed = False

    def add_torsions(self, atom_indices, enableds, spring_constants, targets, kappas, alphas):
        r'''
        Add a set of torsion restraints using a fast C++ function. Returns a
        NumPy integer array giving the indices of the restraints in the force
        object. Fastest if the inputs are NumPy arrays, but most iterables
        should work.

        Args:
            * atom_indices:
                - A 4-tuple of arrays providing the indices of the dihedral
                  atoms in the simulation construct
            * enableds:
                - A Boolean array (or any array castable to float) where values
                  > 0.5 represent enabled restraints
            * spring_constants:
                - Restraint spring constants in :math:`kJ mol^{-1} rad^{-2}`
            * targets:
                - Target angles in radians
            * kappas:
                - Constants defining the width of the restrained well. For
                  :math:`\kappa > 0.5`, :math:`1/\kappa` is approximately equal
                  to the variance of a normal distribution centred on
                  :math:`\theta_0`. The energy approaches
                  :math:`-k*cos(\theta-\theta_0)` as :math:`\kappa` approaches
                  zero. Values less than 0 are not allowed, and will be
                  automatically set to 0.
            * alphas:
                - Constants defining the rate of fall-off outside the restrained
                  well. For :math:`\alpha = 0` the potential in this region is
                  essentially completely flat. Positive values increase the
                  gradient. Values of :math:`\alpha` greater than 1 lead to a
                  steeper gradient *outside* the well than *inside* - this is
                  not recommended in most circumstances. If :math:`\alpha<0`,
                  the potential will *fall* outside the well (that is, torsions
                  that are not already within the well will be repelled from
                  it). Again, this is not recommended.
        '''
        n = len(targets)
        ind = numpy.empty((n,4), numpy.int32)
        for i, ai in enumerate(atom_indices):
            ind[:,i] = ai
        params = numpy.empty((n,5), numpy.float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2] = targets
        params[:,3] = kappas
        params[:,3][params[:,3]<0] = 0
        params[:,4] = alphas
        return _openmm_force_ext.customtorsionforce_add_torsions(int(self.this), ind, params)

    def update_target(self, index, enabled=None, k = None, target = None,
            kappa = None, alpha = None):
        '''
        Update a single restraint. This function is mostly superseded by
        :func:`update_targets`.

        Args:
            * index:
                - integer index of the restraint in the force object
            * enabled:
                - enable/disable the restraint. None keeps the current value.
            * k:
                - set the spring constant in :math:`kJ mol^{-1} rad^{-2}`.
                  None keeps the current value.
            * target:
                - set the target angle in radians. None keeps the current value.
            * kappa:
                - set the kappa (approximate units of inverse square radians)
            * alpha:
                - set the alpha (steepness of the potential outside the well -
                  typically between 0 and 1)
        '''
        # For compatibility with numpy.int32
        index = int(index)
        current_params = self.getTorsionParameters(index)
        indices = current_params[0:4]
        new_enabled, new_k, new_theta0, new_cutoff, new_alpha = current_params[4]
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
        if kappa is not None:
            if type(kappa) == Quantity:
                kappa = kappa.value_in_unit(OPENMM_ANGLE_UNIT)
            new_kappa = kappa
        if alpha is not None:
            new_alpha = alpha
        self.setTorsionParameters(index, *indices, (new_enabled, new_k, new_theta0, new_kappa, new_alpha))
        self.update_needed = True

    def update_targets(self, indices, enableds, spring_constants, targets, kappas, alphas):
        '''
        Update a set of targets all at once using fast C++ code. Fastest if
        the arguments are provided as NumPy arrays, but any iterable should work.

        Args:
            * indices:
                - the indices of the restraints in the OpenMM force object
            * enableds:
                - a Boolean array defining which restraints are to be enabled
            * spring_constants:
                - the new spring constants in units of :math:`kJ mol^{-1} rad^{-2}`
            * targets:
                - the new target angles in radians
            * cutoffs:
                - the new kappas in inverse square radians
            * alpha:
                - the new alphas (steepness of the potential outside the well -
                  typically between 0 and 1)
        '''
        n = len(indices)
        params = numpy.empty((n,5), numpy.float64)
        params[:,0] = enableds
        params[:,1] = spring_constants
        params[:,2] = targets
        params[:,3] = kappas
        params[:,4] = alphas
        _openmm_force_ext.customtorsionforce_update_torsion_parameters(int(self.this), indices, params)
        self.update_needed=True

class TorsionNCSForce(CustomCompoundBondForce):
    '''
    (WORK IN PROGRESS)
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

e_charge = 1.602176634e-19*unit.coulomb
eps0 = 1e-6*8.8541878128e-12/(unit.AVOGADRO_CONSTANT_NA*e_charge**2)*unit.farad/unit.meter
ONE_ON_4_PI_EPS0 = 1/(4*pi*eps0)*eps0.unit


class NonbondedSoftcoreForce(CustomNonbondedForce):
    '''
    Defines a soft-core Lennard-Jones potential to replace the default version,
    to reduce the generation of out-of-gamut forces during energy minimisation.
    Based on the functional form defined in Pham & Shirts (2011), 
    doi: 10.1063/1.3607597
    '''
    LAMBDA_INDEX = 0
    def __init__(self, a=0.5, b=2, c=6, nb_lambda=0.9, alpha=0.2):
        energy_function = ('lennard_jones + coulombic;'
            'lennard_jones = '
                '4 * epsilon * softcore_lambda^softcore_a * '
                '( lj_base^(12/softcore_c) - lj_base^(6/softcore_c) );'
            'lj_base = '
                '1 / ( softcore_alpha * (1-softcore_lambda)^softcore_b +'
                '(r/sigma)^softcore_c );'
            'sigma = 0.5*(sigma1+sigma2);'
            'epsilon = sqrt(epsilon1*epsilon2);'
            f'coulombic = {ONE_ON_4_PI_EPS0} * charge1 * charge2 * softcore_lambda^softcore_a *'
                '( 1 / ( softcore_alpha*(1-softcore_lambda)^softcore_b + r^softcore_c ) )^(1/softcore_c)' 
            )

        super().__init__(energy_function)
        self.addGlobalParameter('softcore_lambda', nb_lambda)
        self.addGlobalParameter('softcore_a', a)
        self.addGlobalParameter('softcore_b', b)
        self.addGlobalParameter('softcore_c', c)
        self.addGlobalParameter('softcore_alpha', alpha)

        self.addPerParticleParameter('charge')
        self.addPerParticleParameter('sigma')
        self.addPerParticleParameter('epsilon')
        self.update_needed = False
    
    def set_lambda(self, value, context=None):
        if value <=0 or value > 1:
            from chimerax.core.errors import UserError
            raise UserError('Lambda must be in the range 0 < lambda <= 1!')
        if context is not None:
            context.setParameter('softcore_lambda', value)
        self.setGlobalParameterDefaultValue(self.LAMBDA_INDEX, value)

class NonbondedSoftcoreExceptionForce(CustomBondForce):
    def __init__(self, a=0.5, b=2, c=6, nb_lambda=0.9, alpha=0.2):
        energy_function = ('lennard_jones + coulombic;'
            'lennard_jones = '
                '4 * epsilon * softcore_lambda^softcore_a * '
                '( lj_base^(12/softcore_c) - lj_base^(6/softcore_c) );'
            'lj_base = '
                '1 / ( softcore_alpha * (1-softcore_lambda)^softcore_b +'
                '(r/sigma)^softcore_c );'
            f'coulombic = {ONE_ON_4_PI_EPS0} * charge_prod * softcore_lambda^softcore_a *'
                '( 1 / ( softcore_alpha*(1-softcore_lambda)^softcore_b + r^softcore_c ) )^(1/softcore_c)' 
            )
        super().__init__(energy_function)
        self.addGlobalParameter('softcore_lambda', nb_lambda)
        self.addGlobalParameter('softcore_a', a)
        self.addGlobalParameter('softcore_b', b)
        self.addGlobalParameter('softcore_c', c)
        self.addGlobalParameter('softcore_alpha', alpha)
        self.addPerBondParameter('charge_prod')
        self.addPerBondParameter('sigma')
        self.addPerBondParameter('epsilon')
        self.update_needed = False
        




class GBSAForce(customgbforces.GBSAGBnForce):
    '''
    Wrapper around :py:class:`openmm.GBSAGBnForce` which implements the
    generalised Born GB-Neck2 implicit solvent implementation.
    '''
    def __init__(self, solventDielectric=78.5, soluteDielectric=1,
                SA='ACE', cutoff=1.0, kappa=3.0,
                nonbonded_method = openmm.CustomGBForce.CutoffNonPeriodic):
        '''
        Initialise the force object. Defaults are chosen to represent a salt
        concentration of approximately 0.5M at 100K, broadly representative of
        the conditions within typical protein crystals.

        Args:
            * solventDielectric:
                - dielectric constant of solvent regions
            * soluteDielectric:
                - dielectric constant "inside" the macromolecule
            * SA:
                - string choosing the method for determining solvent-accessible
                  surface
            * cutoff:
                - cutoff distance in nanometres (must match the cutoff distance
                  for the other nonbonded forces in the simulation!)
            * kappa:
                - Determines the rate of falloff of the potential with distance.
                  Effectively a proxy for salt concentration, where higher kappa
                  corresponds to higher salt.
            * nonbonded_method:
                - should be left as default in almost all cases.
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

class SoftCoreGBSAGBnForce(customgbforces.GBSAGBnForce):
    '''
    Generalised Born GB-Neck2 implicit solvent implementation, modified to work with the 
    :class:`NonbondedSoftCoreForce`.
    '''
    OFFSET=0.009
    def __init__(self, solventDielectric=78.5, soluteDielectric=1, 
                SA='ACE', cutoff=1.0, kappa=3.0,
                nonbonded_method = openmm.CustomGBForce.CutoffNonPeriodic,
                a=0.5, b=2, c=6, nb_lambda=0.9, alpha=0.2):
        if type(solventDielectric) == Quantity:
            solventDielectric = solventDielectric.value_in_unit(OPENMM_DIPOLE_UNIT)
        if type(soluteDielectric) == Quantity:
            soluteDielectric = soluteDielectric.value_in_unit(OPENMM_DIPOLE_UNIT)
        if type(cutoff) == Quantity:
            cutoff = cutoff.value_in_unit(OPENMM_LENGTH_UNIT)
        if type(kappa) == Quantity:
            kappa = kappa.value_in_unit(1/OPENMM_LENGTH_UNIT)
        self._softcore_params = {'softcore_a':a, 'softcore_b':b, 'softcore_c':c, 'softcore_lambda':nb_lambda, 'softcore_alpha': alpha}
        super().__init__(solventDielectric=solventDielectric,
                         soluteDielectric=soluteDielectric,
                         SA=SA,
                         cutoff=cutoff,
                         kappa=kappa)
        self.setCutoffDistance(cutoff)

        self.setNonbondedMethod(nonbonded_method)
        self.update_needed = False
    
    def _addEnergyTerms(self):
        for name, val in self._softcore_params.items():
            self.addGlobalParameter(name, val)
        self.addPerParticleParameter('charge')
        self.addPerParticleParameter('or') # Offset radius
        self.addPerParticleParameter('sr') # Scaled offset radius
        self.addPerParticleParameter('radindex')

        n = len(self._uniqueRadii)
        m0Table = self._createUniqueTable(customgbforces.m0)
        d0Table = self._createUniqueTable(customgbforces.d0)
        self.addTabulatedFunction("getd0", Discrete2DFunction(n, n, d0Table))
        self.addTabulatedFunction("getm0", Discrete2DFunction(n, n, m0Table))

        self.addComputedValue("I", "softcore_lambda^softcore_a * (Ivdw+neckScale*Ineck);"
                                   "Ineck=step(radius1+radius2+neckCut-r)*getm0(radindex1,radindex2)/(1+100*(r-getd0(radindex1,radindex2))^2+"
                                   "0.3*1000000*(r-getd0(radindex1,radindex2))^6);"
                                   "Ivdw=select(step(r+sr2-or1), 0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r), 0);"
                                   "U=r+sr2;"
                                   "L=max(or1, D);"
                                   "D=abs(r-sr2);"
                                   "radius1=or1+offset; radius2=or2+offset;"
                                   "neckScale=0.361825; neckCut=0.68; offset=0.009", self.ParticlePairNoExclusions)

        self.addComputedValue("B", "1/(1/or-tanh(1.09511284*psi-1.907992938*psi^2+2.50798245*psi^3)/radius);"
                                  "psi=I*or; radius=or+offset; offset=0.009", self.SingleParticle)
        self._createEnergyTerms(self.solventDielectric, self.soluteDielectric, self.SA, self.cutoff, self.kappa, self.OFFSET)

        
    def _createEnergyTerms(self, solventDielectric, soluteDielectric, SA, cutoff, kappa, offset):
        from openmm.app.internal.customgbforces import CustomGBForce
        params = (f'; solventDielectric={solventDielectric:.16g}'
                  f'; soluteDielectric={soluteDielectric:.16g}'
                  f'; kappa={kappa:.16g}; offset={offset:.16g}')
        if cutoff is not None:
            params += f'; cutoff={cutoff:.16g}'
        if kappa > 0:
            self.addEnergyTerm(f'softcore_lambda^softcore_a * -0.5*{ONE_ON_4_PI_EPS0}* (1/soluteDielectric-exp(-kappa*B)/solventDielectric)*charge^2/B' + params,
                CustomGBForce.SingleParticle)
        elif kappa < 0:
            raise ValueError('kappa/ionic strength must be >= 0')
        else:
            self.addEnergyParameterDerivative(f'-0.5*{ONE_ON_4_PI_EPS0}* (1/soluteDielectric-1/solventDielectric)*charge^2/B' + params,
                CustomGBForce.SingleParticle)
        if SA=='ACE':
            self.addEnergyTerm("softcore_lambda^softcore_a * 28.3919551* (radius+0.14)^2*(radius/B)^6; radius=or+offset"+params, CustomGBForce.SingleParticle)
        elif SA is not None:
            raise ValueError('Unknown surface area method: '+SA)
        if cutoff is None:
            if kappa > 0:
                self.addEnergyTerm(f'-{ONE_ON_4_PI_EPS0} * (1/soluteDielectric - exp(-kappa*f)/solventDielectric) * softcore_lambda^(2*softcore_a)*charge1*charge2/f;'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions)
            else:
                self.addEnergyTerm(f'-{ONE_ON_4_PI_EPS0}*(1/soluteDielectric - 1/solventDielectric) * softcore_lambda^(2*softcore_a)*charge1*charge2/f;'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions)
        else:
            if kappa > 0:
                self.addEnergyTerm(f'-{ONE_ON_4_PI_EPS0} *  (1/soluteDielectric - exp(-kappa*f)/solventDielectric) *softcore_lambda^(2*softcore_a) *charge1*charge2* (1/f - 1/cutoff);'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions
                )
            else:
                self.addEnergyTerm(f'-{ONE_ON_4_PI_EPS0}*(1/soluteDielectric - 1/solventDielectric) * softcore_lambda^(2*softcore_a)*charge1*charge2 * (1/f-1/cutoff);'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions
                )



            

            


    



