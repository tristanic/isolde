# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tcroll@altoslabs.com
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
    def __init__(self, data, xyz_to_ijk_transform, suffix, units = 'angstroms', map_sigma = 1.0):
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
        map_func = self._map_func = self._openmm_3D_function_from_volume(data)

        global_k_name = self._global_k_name = 'mdff_global_k_{}'.format(suffix)
        magnification_name = self._map_magnification_factor_name = \
            'mdff_magnification_factor_{}'.format(suffix)
        # Map-magnitude normalisation.  The energy is divided by this so the coupling
        # constant is invariant to overall changes in map scale (e.g. when a live
        # crystallographic map's sigma shifts as the mean B-factor changes).  Set to
        # the map sigma at force creation and updated live (see set_map_sigma); the
        # stored coupling constant is therefore a sigma-normalised weight.
        map_sigma_name = self._map_sigma_name = 'mdff_map_sigma_{}'.format(suffix)

        energy_func = self._set_energy_function(suffix)
        self.setEnergyFunction(energy_func)
        #super().__init__(1, energy_func)
        self._map_potential_index = self.addTabulatedFunction(
            name = 'map_potential', function = map_func)
        self._map_magnification_factor_index = self.addGlobalParameter(
            name=magnification_name, defaultValue=1.0
        )
        self._map_sigma_index = self.addGlobalParameter(
            name=map_sigma_name, defaultValue=map_sigma)
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
    def map_magnification_factor_name(self):
        return self._map_magnification_factor_name

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
    
    def update_map_data(self, data):
        '''
        Update the map data values in the simulation. The dimensions of :var:`data` must be 
        *exactly* the same as those of the data used to start the simulation.
        '''
        dim = data.shape[::-1]
        data_1d = numpy.ravel(data, order = 'C')
        self._map_func.setFunctionParameters(*dim, data_1d)
        self.update_needed = True

    def _set_energy_function(self, suffix):
        raise RuntimeError('Cannot instantiate the base class!')

    def set_global_k(self, k, context=None):
        '''
        Set the global coupling constant.  Because the energy is normalised by the
        map sigma (see :func:`set_map_sigma`), this is a *sigma-normalised* weight:
        the effective coupling is ``global_k / map_sigma`` per unit map density.
        '''
        if context is not None:
            context.setParameter(self._global_k_name, k)
        else:
            self.setGlobalParameterDefaultValue(self._global_k_index, k)
            self.update_needed = True

    def set_map_magnification_factor(self, magnification_factor, context=None):
        '''
        Re-scale the map *dimensions* (coordinate magnification, not magnitude).
        Value should typically be very close to 1.0.
        '''
        if context is not None:
            context.setParameter(self._map_magnification_factor_name, magnification_factor)
        else:
            self.setGlobalParameterDefaultValue(
                self._map_magnification_factor_index, magnification_factor)
            self.update_needed = True

    def set_map_sigma(self, map_sigma, context=None):
        '''
        Set the map-magnitude normalisation factor (the current map sigma) that the
        energy is divided by.  Updating this in lockstep with the map data (see
        :func:`update_map_data`) keeps the effective coupling invariant to overall
        changes in map scale, e.g. as a live crystallographic map's sigma shifts
        when the mean B-factor changes.
        '''
        if context is not None:
            context.setParameter(self._map_sigma_name, map_sigma)
        else:
            self.setGlobalParameterDefaultValue(self._map_sigma_index, map_sigma)
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
    def __init__(self, data, xyz_to_ijk_transform, suffix, units = 'angstroms', map_sigma = 1.0):
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
        super().__init__(data, xyz_to_ijk_transform, suffix, units=units, map_sigma=map_sigma)

    def _openmm_3D_function_from_volume(self, data):
        dim = data.shape[::-1]
        data_1d = numpy.ravel(data, order = 'C')
        return Continuous3DFunction(*dim, data_1d, 0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1)

    def _set_energy_function(self, suffix):
        tf = self._transform
        # Transform xyz to ijk
        scale_str = ') / mdff_magnification_factor_{}'.format(suffix)
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
        energy_str = '-{} * individual_k * map_potential(i,j,k) / {}'.format(
            self._global_k_name, self._map_sigma_name)

        final_func = 'select({}, {}, 0); {}'.format(enabled_eqn, energy_str, funcs)
        return final_func

class CubicInterpMapForce(_Map_Force_Base):
    '''
    Creates a MDFF potential from a 3D map of density values.
    '''
    def __init__(self, data, xyz_to_ijk_transform, suffix, units = 'angstroms', map_sigma = 1.0):
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
        super().__init__(data, xyz_to_ijk_transform, suffix, units=units, map_sigma=map_sigma)

    def _set_energy_function(self, suffix):
        tf = self._transform
        # Transform xyz to ijk
        scale_str = ') / mdff_magnification_factor_{}'.format(suffix)
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

        energy_str = '-{} * individual_k * ({}) / {}'.format(
            self._global_k_name, master_sum_str, self._map_sigma_name)

        enabled_eqn = 'step(enabled-0.5)'

        funcs = ';'.join((coeff_strings, offset_strings, min_str, i_str, j_str, k_str))
        final_func = 'select({}, {}, 0); {}'.format(enabled_eqn, energy_str, funcs)
        return final_func

    def _openmm_3D_function_from_volume(self, data):
        dim = data.shape[::-1]
        data_1d = numpy.ravel(data, order = 'C')
        return Discrete3DFunction(*dim, data_1d)


class SymmetryAwareCubicInterpMapForce(CubicInterpMapForce):
    '''
    A :class:`CubicInterpMapForce` in which each term (bond) carries its own
    orthogonal-space transform ``[R | t]`` (R dimensionless, t in nanometres),
    so the map is sampled at ``xs = R.x + t`` rather than at the atom's own
    position ``x``. With the identity transform ``[I | 0]`` this reproduces the
    base force exactly (to machine precision).

    This lets a *real* atom feel the map through a crystallographic symmetry
    operator: a term with ``S = [R | t]`` samples the density at the atom's ghost
    position ``S.x``, and because ``CustomCompoundBondForce`` differentiates the
    whole (nested) energy expression, the force on the real atom is
    ``R^T . grad(map)`` -- exactly the ``SymmetrySite`` fold-back, computed
    directly on the parent. The map force therefore never needs terms on the copy
    (virtual-site) particles, and correctly handles atoms that enter the map only
    under symmetry (they simply get no identity term).

    Per-bond parameter layout (must match the packing in :func:`add_atoms` /
    :func:`update_atoms`): ``[individual_k, enabled, R (row-major 9), t (3)]``.
    '''
    # Order MUST match the packing in add_atoms/update_atoms (cols 2..13).
    _TRANSFORM_PARAMS = ('mdff_r00', 'mdff_r01', 'mdff_r02',
                         'mdff_r10', 'mdff_r11', 'mdff_r12',
                         'mdff_r20', 'mdff_r21', 'mdff_r22',
                         'mdff_t0', 'mdff_t1', 'mdff_t2')
    N_PER_BOND_PARAMS = 14
    # Identity transform (R = I, t = 0) as the 12 flattened columns.
    IDENTITY_TRANSFORM = numpy.array(
        [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        dtype=numpy.float64)

    def __init__(self, data, xyz_to_ijk_transform, suffix, units='angstroms',
            map_sigma=1.0):
        super().__init__(data, xyz_to_ijk_transform, suffix, units=units,
            map_sigma=map_sigma)
        # individual_k (col 0) and enabled (col 1) were added by the base;
        # append the 12 transform parameters (cols 2..13). The energy string
        # (built in the base __init__ via our overridden _set_energy_function)
        # references these by name; OpenMM only resolves them at context
        # creation, so adding them here -- after setEnergyFunction -- is fine.
        for name in self._TRANSFORM_PARAMS:
            self.addPerBondParameter(name)
        assert self.getNumPerBondParameters() == self.N_PER_BOND_PARAMS, \
            'SymmetryAwareCubicInterpMapForce per-bond parameter count mismatch'

    def _set_energy_function(self, suffix):
        # Start from the base cubic-interpolation energy, then redirect the
        # sampled coordinate from the atom position (x1,y1,z1) to the transformed
        # position (xs,ys,zs). x1/y1/z1 appear ONLY in the i/j/k transform strings
        # (via the 'x1 * {}' entries); every other token (ci1, min_i, sc_i1,
        # map_potential(...), ...) contains no 'x1'/'y1'/'z1' substring, so the
        # scoped rename is safe.
        import re
        base = super()._set_energy_function(suffix)
        base = re.sub(r'([xyz])1', r'\1s', base)
        transform = (
            'xs = mdff_r00*x1 + mdff_r01*y1 + mdff_r02*z1 + mdff_t0;'
            'ys = mdff_r10*x1 + mdff_r11*y1 + mdff_r12*z1 + mdff_t1;'
            'zs = mdff_r20*x1 + mdff_r21*y1 + mdff_r22*z1 + mdff_t2')
        return base + ';' + transform

    def add_atoms(self, indices, ks, enableds, transforms=None):
        '''
        Add a set of terms. ``transforms`` is an ``(n, 12)`` float array of
        flattened ``[R | t_nm]`` operators (row-major R then t in nanometres), or
        ``None`` for all-identity (plain non-symmetry terms).
        '''
        n = len(indices)
        params = numpy.empty((n, self.N_PER_BOND_PARAMS), numpy.float64)
        params[:, 0] = ks
        params[:, 1] = enableds
        if transforms is None:
            params[:, 2:] = self.IDENTITY_TRANSFORM
        else:
            params[:, 2:] = transforms
        return _openmm_force_ext.customcompoundbondforce_add_bonds(
            int(self.this), indices, params)

    def update_atoms(self, indices, ks, enableds, transforms):
        '''
        Update terms in place. The C++ helper overwrites *all* per-bond
        parameters, so the (fixed) per-term transforms must be re-supplied
        (``(n, 12)`` array, same layout as :func:`add_atoms`).
        '''
        n = len(indices)
        params = numpy.empty((n, self.N_PER_BOND_PARAMS), numpy.float64)
        params[:, 0] = ks
        params[:, 1] = enableds
        params[:, 2:] = transforms
        _openmm_force_ext.customcompoundbondforce_update_bond_parameters(
            int(self.this), indices, params)
        self.update_needed = True


class LinearInterpMapForce(_Map_Force_Base):
    '''
    NOTE: This class is deprecated, since there is no conceivable situation in which
    it is superior to either :class:`CubicInterpMapForce` or
    :class:`CubicInterpMapForce_Old`.

    Converts a map of (i,j,k) data and a (x,y,z)->(i,j,k) transformation
    matrix to a potential energy field, with trilinear interpolation of values.
    '''
    def __init__(self, data, xyz_to_ijk_transform, suffix, units = 'angstroms', map_sigma = 1.0):
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
        super().__init__(data, xyz_to_ijk_transform, suffix, units=units, map_sigma=map_sigma)

    def _set_energy_function(self, suffix):
        tf = self._transform
        # Transform xyz to ijk
        scale_str = ') / mdff_magnification_factor_{}'.format(suffix)
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

        energy_str = '-{} * individual_k * ({}) / {}'.format(
            self._global_k_name, interp_str, self._map_sigma_name)

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
    loss function described by `Jonathan Barron`
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
    def __init__(self, a=1, b=2, c=6, nb_lambda=0.9, alpha=0.2):
        '''
        Note: for best simulation performance a, b and c should be integers.
        '''
        super().__init__(self._soft_core_energy(a, b, c))
        self.addGlobalParameter('softcore_lambda', nb_lambda)
        self.addGlobalParameter('softcore_alpha', alpha)

        self.addPerParticleParameter('charge')
        self.addPerParticleParameter('sigma')
        self.addPerParticleParameter('epsilon')
        self.update_needed = False

    @staticmethod
    def _soft_core_energy(a, b, c, lam='softcore_lambda', extra_defs=''):
        '''
        Build the soft-core LJ + Coulomb energy expression. ``lam`` is the name of
        the per-pair coupling variable used throughout; the default
        ``'softcore_lambda'`` (the global parameter) reproduces the plain
        force byte-for-byte. Subclasses that make the coupling per-group pass
        ``lam='pair_lambda'`` and supply ``extra_defs`` defining ``pair_lambda``
        (e.g. from a group-pair table). The energy value is the first clause;
        ``extra_defs`` is appended at the *end*. OpenMM's Lepton parser requires
        an intermediate variable to be defined *after* the expressions that use
        it (matching the ``lennard_jones``-uses-``lj_base`` convention here), so
        it must be appended, never prepended.

        Only the ``lennard_jones``/``lj_base`` block is Lennard-Jones-specific;
        it references ``lam`` and is the single point that a future
        double-exponential vdW form would replace, leaving the coupling
        plumbing untouched.
        '''
        return (
            'lennard_jones + coulombic;'
            'lennard_jones = '
                f'4 * epsilon * {lam}^(1/{a}) * '
                f'( lj_base^(12/{c}) - lj_base^(6/{c}) );'
            'lj_base = '
                f'1 / ( softcore_alpha * (1-{lam})^{b} +'
                f'(r/sigma)^{c} );'
            'sigma = 0.5*(sigma1+sigma2);'
            'epsilon = sqrt(epsilon1*epsilon2);'
            f'coulombic = {ONE_ON_4_PI_EPS0} * charge1 * charge2 * '
                f'( 1 / ( softcore_alpha*(1-{lam})^({b*4}) + r^{c} ) )^(1/{c})'
            + ((';' + extra_defs.rstrip(';')) if extra_defs else '')
            )
    
    def set_lambda(self, value, context=None):
        if value <=0 or value > 1:
            from chimerax.core.errors import UserError
            raise UserError('Lambda must be in the range 0 < lambda <= 1!')
        if context is not None:
            context.setParameter('softcore_lambda', value)
        self.setGlobalParameterDefaultValue(self.LAMBDA_INDEX, value)
    
    @staticmethod
    def potential_values(radii, nb_lambda, a, b, c, alpha, charge=-0.1):
        '''
        Given a Numpy array of radii (in units of nm), returns the corresponding adjusted Lennard-Jones
        and Coulombic potential values corresponding to the given :var:`nb_lambda`, :var:`a`, 
        :var:`b`, :var:`c` and :var:`alpha` for a hypothetical pair of oxygen atoms each with the specified 
        charge, in kJ/mol.
        '''

        sigma = 0.295992190115
        epsilon = 0.87864
        coulombic = ONE_ON_4_PI_EPS0 * charge**2 * (
            1 / (alpha * (1-nb_lambda)**(b*4) + radii**c )
        ) ** (1/c)
        lj_base = 1 / ( alpha * (1-nb_lambda)**b + (radii/sigma)**c )
        lennard_jones = 4 * epsilon * nb_lambda**(1/a) * ( lj_base**(12/c) - lj_base**(6/c))
        return lennard_jones, coulombic




class NonbondedSoftcoreExceptionForce(CustomBondForce):
    def __init__(self, a=1, b=2, c=6, nb_lambda=0.9, alpha=0.2):
        energy_function = ('lennard_jones + coulombic;'
            'lennard_jones = '
                f'4 * epsilon * softcore_lambda^(1/{a}) * '
                f'( lj_base^(12/{c}) - lj_base^(6/{c}) );'
            'lj_base = '
                f'1 / ( softcore_alpha * (1-softcore_lambda)^{b} +'
                f'(r/sigma)^{c} );'
            f'coulombic = {ONE_ON_4_PI_EPS0} * charge_prod * '
                f'( 1 / ( softcore_alpha*(1-softcore_lambda)^({4*b}) + r^{c} ) )^(1/{c})' 
            )
        super().__init__(energy_function)
        self.addGlobalParameter('softcore_lambda', nb_lambda)
        self.addGlobalParameter('softcore_alpha', alpha)
        self.addPerBondParameter('charge_prod')
        self.addPerBondParameter('sigma')
        self.addPerBondParameter('epsilon')
        self.update_needed = False


def symmetry_group_table(n, exclude_intra_operator=True, weights=None):
    '''
    Build the row-major ``n*n`` weight table for the ``grouptable``
    :py:class:`openmm.Discrete2DFunction` shared by the symmetry-aware nonbonded
    and GBSA pairwise terms. Group 0 is the real asymmetric unit; groups
    ``1..n-1`` are one-per-operator symmetry copies.

    If ``weights`` is given (a precomputed flattened ``n*n`` table from
    :func:`chimerax.isolde.openmm.symmetry_sim.symmetry_group_weight_table`), it is
    returned verbatim. That operator-set-aware table is the correct one: the fixed
    copy<->copy = 0 below is only valid when the whole construct is simulated, and
    silently drops cross-operator copy<->copy contacts in local simulations at
    multi-way interfaces. The fixed-weight branch here is kept as a fallback / for
    the ``exclude_intra_operator=False`` demonstration hook.

    Entries are multiplicative weights on the pairwise energy that make each
    *unique* crystallographic contact contribute exactly once while remaining
    fully two-way (validated: ``scratchpad/symsite_double_count.py``):

        * ``[0][0] = 1``            -- real atoms interact with each other.
        * real<->copy = ``1/2``     -- an inter-molecular contact between two
          asymmetric-unit atoms i, j is generated from *both* sides
          (real_i<->copy_j AND real_j<->copy_i, which are the same physical
          contact). Each is half-weighted so their sum, folded back onto the two
          real parents, is 1x the contact -- not the 2x double-count that
          otherwise collapses salt bridges across symmetry interfaces. The
          self-image case (atom near its own symmetry axis) likewise becomes 1x.
        * copy<->copy = ``0``       -- interactions between two symmetry copies
          (any operators) are neighbour<->neighbour contacts that belong to
          other cells' sums; folding them onto the ASU would over-count. Every
          real contact is already captured by a real<->copy term.

    ``exclude_intra_operator=False`` returns the naive all-ones table (every pair
    full weight); it exists only to *demonstrate* the pathology the correct
    weighting prevents -- the inter-molecular double-count and the explosive
    same-operator clash between copies of bonded atoms.
    '''
    if weights is not None:
        if len(weights) != n * n:
            raise ValueError(
                'symmetry group weight table has {} entries, expected {}'.format(
                    len(weights), n * n))
        return list(weights)
    if not exclude_intra_operator:
        return [1.0] * (n * n)

    def val(i, j):
        if i == 0 and j == 0:
            return 1.0          # real <-> real
        if i == 0 or j == 0:
            return 0.5          # real <-> copy (each unique contact once, two-way)
        return 0.0              # copy <-> copy (any operators): excluded
    return [val(i, j) for j in range(n) for i in range(n)]


class NBGroupNonbondedSoftcoreForce(NonbondedSoftcoreForce):
    '''
    Per-group soft-core nonbonded force. Adds a per-particle integer ``nb_group``
    and an ``N x N`` :py:class:`openmm.Discrete2DFunction` coupling matrix
    ``nb_coupling_table``; the effective per-pair coupling used throughout the
    soft-core LJ/Coulomb terms is::

        pair_lambda = min(softcore_lambda, nb_coupling_table(nb_group1, nb_group2))

    i.e. the global ``softcore_lambda`` acts as a ceiling (softest-wins) and the
    table softens the coupling *between* groups. ``nb_coupling_table`` is the
    identity (all 1.0) at construction, so with every atom in group 0 the force
    reduces to the plain :py:class:`NonbondedSoftcoreForce`. ``addParticle`` takes
    ``[charge, sigma, epsilon, nb_group]``.

    ``n_nb_groups == 1`` degrades to the plain force byte-for-byte (no ``nb_group``
    parameter, no table) so a caller can construct this class unconditionally.

    The base class is used for symmetry-aware forces too (see
    :py:class:`SymmetryAwareNonbondedSoftcoreForce`), so per-group soft-core is
    available with or without crystallographic symmetry.
    '''
    def __init__(self, a=1, b=2, c=6, nb_lambda=0.9, alpha=0.2, n_nb_groups=1):
        n = int(n_nb_groups)
        if n <= 1:
            # No groups requested: behave exactly as the plain force.
            super().__init__(a=a, b=b, c=c, nb_lambda=nb_lambda, alpha=alpha)
            self._n_nb_groups = 1
            self._nb_values = None
            self._nb_coupling_table = None
            return
        # Grouped: build the pair_lambda expression. We bypass
        # NonbondedSoftcoreForce.__init__ (which builds the plain expression) and
        # construct the CustomNonbondedForce directly with the group form.
        extra = ('pair_lambda = min(softcore_lambda, '
                 'nb_coupling_table(nb_group1, nb_group2));')
        energy = self._soft_core_energy(a, b, c, lam='pair_lambda', extra_defs=extra)
        CustomNonbondedForce.__init__(self, energy)
        self.addGlobalParameter('softcore_lambda', nb_lambda)
        self.addGlobalParameter('softcore_alpha', alpha)
        self.addPerParticleParameter('charge')
        self.addPerParticleParameter('sigma')
        self.addPerParticleParameter('epsilon')
        self.addPerParticleParameter('nb_group')
        self._nb_group_index = 3       # position of nb_group in the per-particle list
        self._n_nb_groups = n
        self._nb_values = [1.0] * (n * n)      # identity: every pair fully coupled
        self._nb_coupling_table = Discrete2DFunction(n, n, list(self._nb_values))
        self.addTabulatedFunction('nb_coupling_table', self._nb_coupling_table)
        self.update_needed = False

    @property
    def n_nb_groups(self):
        return self._n_nb_groups

    def get_coupling(self, group_a, group_b):
        if self._nb_values is None:
            return 1.0
        return self._nb_values[group_a * self._n_nb_groups + group_b]

    def set_coupling(self, group_a, group_b, lam, context=None):
        '''
        Set the symmetric coupling between two groups (0 < lam <= 1). Flags the
        force for a live parameter update; if ``context`` is given, pushes it
        immediately via ``updateParametersInContext`` (no reinitialisation).
        '''
        if self._nb_coupling_table is None:
            raise ValueError('this force was built with n_nb_groups == 1 (no groups)')
        if lam <= 0 or lam > 1:
            from chimerax.core.errors import UserError
            raise UserError('nb group coupling must be in the range 0 < lambda <= 1!')
        n = self._n_nb_groups
        self._nb_values[group_a * n + group_b] = lam
        self._nb_values[group_b * n + group_a] = lam
        self._nb_coupling_table.setFunctionParameters(n, n, list(self._nb_values))
        self.update_needed = True
        if context is not None:
            self.updateParametersInContext(context)
            self.update_needed = False


class SymmetryAwareMixin:
    '''
    Turns *any* :py:class:`openmm.CustomNonbondedForce` subclass into a
    crystallographic-symmetry-aware one, applied at the **class** level:

        class SymmetryAwareNonbondedSoftcoreForce(
                SymmetryAwareMixin, NonbondedSoftcoreForce):
            pass

    Every particle carries an integer ``symgroup`` per-particle parameter:

        * ``0``      -- a real asymmetric-unit atom.
        * ``1..M``   -- a symmetry copy, one contiguous id per crystallographic
          operator actually present in the simulation.

    The per-pair energy is multiplied by a 0/1 mask ``grouptable(symgroup1,
    symgroup2)`` (an :py:class:`openmm.Discrete2DFunction`) whose entries are:

        * ``[0][0] = 1``            -- real atoms interact with each other.
        * ``[k][k] = 0`` (k > 0)    -- two copies produced by the *same*
          operator do **not** interact. This is the single correctness rule
          established empirically (see the project's symmetry-sim notes): it
          both prevents the explosive clash between copies of bonded atoms and
          stops each molecule's internal strain being counted twice (the copy
          otherwise replicates the ASU's internal energy and folds it straight
          back onto the parents).
        * off-diagonal ``= 1``      -- real<->copy and *different*-operator
          copy<->copy pairs are genuine crystal contacts and interact in full.
          (The 2x force enhancement at self-inverse operators is correct
          constrained-gradient physics and is deliberately left intact.)

    Because the mask is a purely multiplicative factor it composes with any
    pair potential without touching its functional form -- validated to
    reproduce explicit copy-copy exceptions to machine precision for both
    Lennard-Jones and a double-exponential (GARNET-style) potential.

    The base class builds its own energy expression and per-particle
    parameters in ``__init__``; this mixin runs *after* that (via ``super()``,
    while the force is still particle-free), wraps the head of the expression
    and appends the ``symgroup`` parameter + ``grouptable``. ``addParticle``
    then takes the base parameters **plus** the group id as a trailing value.
    '''
    def __init__(self, *args, symmetry_ngroups=2, exclude_intra_operator=True,
            symmetry_group_weights=None, **kwargs):
        super().__init__(*args, **kwargs)
        # Split the base expression on the first ';'. The head is the per-pair
        # energy value; everything after it is the chain of intermediate
        # definitions, which we preserve verbatim.
        head, _, tail = self.getEnergyFunction().partition(';')
        new_expr = (
            f'symmetry_switch*({head});'
            'symmetry_switch=grouptable(symgroup1,symgroup2)'
            + ((';' + tail) if tail else '')
        )
        self.setEnergyFunction(new_expr)
        self.addPerParticleParameter('symgroup')
        n = int(symmetry_ngroups)
        # symmetry_group_weights (if given) is the operator-set-aware table; else
        # fall back to the fixed-weight scheme.
        table = symmetry_group_table(n, exclude_intra_operator,
            weights=symmetry_group_weights)
        self.addTabulatedFunction('grouptable', Discrete2DFunction(n, n, table))


class SymmetryAwareNonbondedSoftcoreForce(SymmetryAwareMixin,
        NBGroupNonbondedSoftcoreForce):
    '''
    Symmetry-aware form of :py:class:`NBGroupNonbondedSoftcoreForce` (which is
    itself a :py:class:`NonbondedSoftcoreForce`). Symmetry-awareness thus layers
    on top of the per-group soft-core capability, so per-group nonbonded coupling
    is available whether or not crystallographic symmetry is in use.

    Constructed with ``symmetry_ngroups = 1 + (distinct operators present)`` and,
    optionally, ``n_nb_groups > 1`` to also enable per-group coupling. With
    ``n_nb_groups == 1`` (the default) the nb-group layer is inert and this is
    byte-identical to the pre-existing symmetry force: ``addParticle`` takes
    ``[charge, sigma, epsilon, symgroup]``. With both active it takes
    ``[charge, sigma, epsilon, nb_group, symgroup]`` (base adds ``nb_group``, the
    mixin appends ``symgroup``). See :py:class:`SymmetryAwareMixin` and
    :py:class:`NBGroupNonbondedSoftcoreForce`.
    '''
    pass




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
                a=1, b=2, c=6, nb_lambda=0.9, alpha=0.2):
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
        self.addGlobalParameter('softcore_lambda', self._softcore_params['softcore_lambda'])
        self.addPerParticleParameter('charge')
        self.addPerParticleParameter('or') # Offset radius
        self.addPerParticleParameter('sr') # Scaled offset radius
        self.addPerParticleParameter('radindex')

        n = len(self._uniqueRadii)
        m0Table = self._createUniqueTable(customgbforces.m0)
        d0Table = self._createUniqueTable(customgbforces.d0)
        self.addTabulatedFunction("getd0", Discrete2DFunction(n, n, d0Table))
        self.addTabulatedFunction("getm0", Discrete2DFunction(n, n, m0Table))
        a = self._softcore_params['softcore_a']

        self.addComputedValue("I", f"softcore_lambda^{a} * (Ivdw+neckScale*Ineck);"
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
        a = self._softcore_params['softcore_a']
        if cutoff is not None:
            params += f'; cutoff={cutoff:.16g}'
        if kappa > 0:
            self.addEnergyTerm(f'softcore_lambda^{a} * -0.5*{ONE_ON_4_PI_EPS0}* (1/soluteDielectric-exp(-kappa*B)/solventDielectric)*charge^2/B' + params,
                CustomGBForce.SingleParticle)
        elif kappa < 0:
            raise ValueError('kappa/ionic strength must be >= 0')
        else:
            self.addEnergyParameterDerivative(f'-0.5*{ONE_ON_4_PI_EPS0}* (1/soluteDielectric-1/solventDielectric)*charge^2/B' + params,
                CustomGBForce.SingleParticle)
        if SA=='ACE':
            self.addEnergyTerm(f"softcore_lambda^{a} * 28.3919551* (radius+0.14)^2*(radius/B)^6; radius=or+offset"+params, CustomGBForce.SingleParticle)
        elif SA is not None:
            raise ValueError('Unknown surface area method: '+SA)
        if cutoff is None:
            if kappa > 0:
                self.addEnergyTerm(f'-{ONE_ON_4_PI_EPS0} * (1/soluteDielectric - exp(-kappa*f)/solventDielectric) * softcore_lambda^({2*a})*charge1*charge2/f;'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions)
            else:
                self.addEnergyTerm(f'-{ONE_ON_4_PI_EPS0}*(1/soluteDielectric - 1/solventDielectric) * softcore_lambda^({2*a})*charge1*charge2/f;'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions)
        else:
            if kappa > 0:
                self.addEnergyTerm(f'-{ONE_ON_4_PI_EPS0} *  (1/soluteDielectric - exp(-kappa*f)/solventDielectric) *softcore_lambda^({2*a}) *charge1*charge2* (1/f - 1/cutoff);'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions
                )
            else:
                self.addEnergyTerm(f'-{ONE_ON_4_PI_EPS0}*(1/soluteDielectric - 1/solventDielectric) * softcore_lambda^({2*a})*charge1*charge2 * (1/f-1/cutoff);'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params, CustomGBForce.ParticlePairNoExclusions
                )


class NBGroupSoftCoreGBSAGBnForce(SoftCoreGBSAGBnForce):
    '''
    Per-group soft-core GB-Neck2 implicit solvent, matching
    :py:class:`NBGroupNonbondedSoftcoreForce`. The descreening integral ``I`` and
    the pairwise polar term scale by the per-group-pair coupling::

        pair_lambda = min(softcore_lambda, nb_coupling_table(nb_group1, nb_group2))

    (``softcore_lambda^a`` -> ``pair_lambda^a`` in ``I``; ``softcore_lambda^(2a)``
    -> ``pair_lambda^(2a)`` in the pairwise term), so a group decoupled from the
    nonbonded force is correspondingly desolvated/screened. The single-particle
    self-energy terms keep the global ``softcore_lambda^a``: a self term has no
    group pair, and its group dependence already flows through the (now
    group-weighted) Born radius ``B``.

    ``n_nb_groups == 1`` degrades to :py:class:`SoftCoreGBSAGBnForce` byte-for-byte.
    Before ``finalize()`` the caller sets ``self._nb_groups`` (per-particle group id,
    in particle order); ``addParticle`` then carries the group id after ``radindex``.
    '''
    def __init__(self, *args, n_nb_groups=1, **kwargs):
        self._n_nb_groups = int(n_nb_groups)
        self._nb_groups = None          # per-particle group ids (set before finalize)
        self._nb_values = None
        self._nb_coupling_table = None
        super().__init__(*args, **kwargs)

    def _addEnergyTerms(self):
        if self._n_nb_groups <= 1:
            return super()._addEnergyTerms()
        a = self._softcore_params['softcore_a']
        self.addGlobalParameter('softcore_lambda', self._softcore_params['softcore_lambda'])
        self.addPerParticleParameter('charge')
        self.addPerParticleParameter('or')          # Offset radius
        self.addPerParticleParameter('sr')          # Scaled offset radius
        self.addPerParticleParameter('radindex')
        self.addPerParticleParameter('nb_group')
        self._nb_group_index = 4       # position of nb_group in the per-particle list

        n = len(self._uniqueRadii)
        m0Table = self._createUniqueTable(customgbforces.m0)
        d0Table = self._createUniqueTable(customgbforces.d0)
        self.addTabulatedFunction("getd0", Discrete2DFunction(n, n, d0Table))
        self.addTabulatedFunction("getm0", Discrete2DFunction(n, n, m0Table))
        ng = self._n_nb_groups
        self._nb_values = [1.0] * (ng * ng)
        self._nb_coupling_table = Discrete2DFunction(ng, ng, list(self._nb_values))
        self.addTabulatedFunction("nb_coupling_table", self._nb_coupling_table)

        # Descreening integral scaled by the per-group-pair coupling (pair_lambda^a).
        self.addComputedValue("I", f"pair_lambda^{a} * (Ivdw+neckScale*Ineck);"
                                   "Ineck=step(radius1+radius2+neckCut-r)*getm0(radindex1,radindex2)/(1+100*(r-getd0(radindex1,radindex2))^2+"
                                   "0.3*1000000*(r-getd0(radindex1,radindex2))^6);"
                                   "Ivdw=select(step(r+sr2-or1), 0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r), 0);"
                                   "U=r+sr2;"
                                   "L=max(or1, D);"
                                   "D=abs(r-sr2);"
                                   "radius1=or1+offset; radius2=or2+offset;"
                                   "neckScale=0.361825; neckCut=0.68; offset=0.009;"
                                   "pair_lambda=min(softcore_lambda, nb_coupling_table(nb_group1,nb_group2))",
                                   self.ParticlePairNoExclusions)

        self.addComputedValue("B", "1/(1/or-tanh(1.09511284*psi-1.907992938*psi^2+2.50798245*psi^3)/radius);"
                                  "psi=I*or; radius=or+offset; offset=0.009", self.SingleParticle)
        self._createEnergyTerms(self.solventDielectric, self.soluteDielectric, self.SA, self.cutoff, self.kappa, self.OFFSET)

    def _createEnergyTerms(self, solventDielectric, soluteDielectric, SA, cutoff, kappa, offset):
        if self._n_nb_groups <= 1:
            return super()._createEnergyTerms(solventDielectric, soluteDielectric, SA, cutoff, kappa, offset)
        from openmm.app.internal.customgbforces import CustomGBForce
        params = (f'; solventDielectric={solventDielectric:.16g}'
                  f'; soluteDielectric={soluteDielectric:.16g}'
                  f'; kappa={kappa:.16g}; offset={offset:.16g}')
        a = self._softcore_params['softcore_a']
        if cutoff is not None:
            params += f'; cutoff={cutoff:.16g}'
        # pair_lambda definition appended (define-after-use) to each pairwise term.
        plam = ';pair_lambda=min(softcore_lambda, nb_coupling_table(nb_group1,nb_group2))'
        # Single-particle self-energy: keep the global softcore_lambda^a (no group pair).
        if kappa > 0:
            self.addEnergyTerm(f'softcore_lambda^{a} * -0.5*{ONE_ON_4_PI_EPS0}* (1/soluteDielectric-exp(-kappa*B)/solventDielectric)*charge^2/B' + params,
                CustomGBForce.SingleParticle)
        elif kappa < 0:
            raise ValueError('kappa/ionic strength must be >= 0')
        else:
            self.addEnergyParameterDerivative(f'-0.5*{ONE_ON_4_PI_EPS0}* (1/soluteDielectric-1/solventDielectric)*charge^2/B' + params,
                CustomGBForce.SingleParticle)
        if SA=='ACE':
            self.addEnergyTerm(f"softcore_lambda^{a} * 28.3919551* (radius+0.14)^2*(radius/B)^6; radius=or+offset"+params, CustomGBForce.SingleParticle)
        elif SA is not None:
            raise ValueError('Unknown surface area method: '+SA)
        # Pairwise polar term: scale by pair_lambda^(2a).
        if cutoff is None:
            if kappa > 0:
                self.addEnergyTerm(f'-{ONE_ON_4_PI_EPS0} * (1/soluteDielectric - exp(-kappa*f)/solventDielectric) * pair_lambda^({2*a})*charge1*charge2/f;'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params+plam, CustomGBForce.ParticlePairNoExclusions)
            else:
                self.addEnergyTerm(f'-{ONE_ON_4_PI_EPS0}*(1/soluteDielectric - 1/solventDielectric) * pair_lambda^({2*a})*charge1*charge2/f;'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params+plam, CustomGBForce.ParticlePairNoExclusions)
        else:
            if kappa > 0:
                self.addEnergyTerm(f'-{ONE_ON_4_PI_EPS0} *  (1/soluteDielectric - exp(-kappa*f)/solventDielectric) *pair_lambda^({2*a}) *charge1*charge2* (1/f - 1/cutoff);'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params+plam, CustomGBForce.ParticlePairNoExclusions
                )
            else:
                self.addEnergyTerm(f'-{ONE_ON_4_PI_EPS0}*(1/soluteDielectric - 1/solventDielectric) * pair_lambda^({2*a})*charge1*charge2 * (1/f-1/cutoff);'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params+plam, CustomGBForce.ParticlePairNoExclusions
                )

    def _addParticles(self):
        if self._n_nb_groups <= 1:
            return super()._addParticles()
        from openmm import CustomGBForce
        for i, p in enumerate(self.parameters):
            radIndex = self._radiusToIndex[p[self.RADIUS_ARG_POSITION]]
            CustomGBForce.addParticle(self, p + [radIndex, float(self._nb_groups[i])])

    def set_coupling(self, group_a, group_b, lam, context=None):
        '''Set symmetric group-pair coupling (0 < lam <= 1); optionally push live.'''
        if self._nb_coupling_table is None:
            raise ValueError('this GB force was built with n_nb_groups == 1 (no groups)')
        if lam <= 0 or lam > 1:
            from chimerax.core.errors import UserError
            raise UserError('nb group coupling must be in the range 0 < lambda <= 1!')
        n = self._n_nb_groups
        self._nb_values[group_a * n + group_b] = lam
        self._nb_values[group_b * n + group_a] = lam
        self._nb_coupling_table.setFunctionParameters(n, n, list(self._nb_values))
        self.update_needed = True
        if context is not None:
            self.updateParametersInContext(context)
            self.update_needed = False


class SymmetrySoftCoreGBSAGBnForce(NBGroupSoftCoreGBSAGBnForce):
    '''
    Crystallographic-symmetry-aware form of :py:class:`SoftCoreGBSAGBnForce`
    (GB-Neck2 implicit solvent, soft-core variant), implementing option (a') of
    the symmetry-sim design.

    Generalised-Born radii are a *collective* property (each particle's ``B`` is
    an integral over every other particle), so implicit solvent cannot be split
    into a separate "copies" force the way pairwise exclusions can -- every
    particle, real or copy, must live in this one force. Symmetry-correctness is
    achieved by masking only the terms that would otherwise double-count a
    parent's own solvation once its copy folds force back onto it:

        * **Descreening integral ``I`` -- left UNMASKED.** A copy self-descreens
          through the normal collective integral, which (by rigid-transform
          invariance) reproduces ``B_copy == B_parent`` to machine precision in
          a complete environment. Real atoms are thereby correctly desolvated by
          their symmetry neighbours, and no explicit Born-radius propagation is
          needed. (Validated: ``scratchpad/gbsa_born_experiment.py``.)
        * **Both single-particle self-energy terms (polar ``charge^2/B`` and the
          nonpolar ACE surface area) -- masked to real atoms** via
          ``step(0.5 - symgroup)``. A copy's self-energy equals its parent's, so
          without this the parent's solvation would be counted twice.
        * **The pairwise polar term -- masked same-operator** via
          ``grouptable(symgroup1, symgroup2)`` (the identical 0/1 table used by
          the nonbonded force: ``[0][0]=1``, ``[k][k]=0`` for ``k>0``, off-
          diagonal ``1``). Real<->copy and cross-operator copy<->copy solvent
          screening (genuine crystal-contact desolvation) is retained.

    Before :func:`finalize`, the caller must set:

        * ``self._symmetry_ngroups`` -- ``1 + (distinct operators present)``.
        * ``self._symmetry_groups``  -- per-particle group id (0 = real, 1..M =
          copy), in particle order and the same length as ``self.parameters``.
    '''

    def _addEnergyTerms(self):
        # When per-group coupling is also active (n_nb_groups > 1) this force carries
        # BOTH group masks: the symmetry mask (grouptable / step(0.5-symgroup)) and the
        # per-group soft-core coupling (pair_lambda from nb_coupling_table). The
        # softening variable ``lam`` is ``pair_lambda`` then, else the plain global
        # ``softcore_lambda`` (byte-identical to the pure-symmetry force).
        nb = self._n_nb_groups > 1
        lam = 'pair_lambda' if nb else 'softcore_lambda'
        # pair_lambda definition, appended (define-after-use) to the pairwise terms.
        plam = (';pair_lambda=min(softcore_lambda, nb_coupling_table(nb_group1,nb_group2))'
                if nb else '')
        self.addGlobalParameter('softcore_lambda',
            self._softcore_params['softcore_lambda'])
        self.addPerParticleParameter('charge')
        self.addPerParticleParameter('or')        # Offset radius
        self.addPerParticleParameter('sr')        # Scaled offset radius
        self.addPerParticleParameter('radindex')
        if nb:
            self.addPerParticleParameter('nb_group')
            self._nb_group_index = 4              # position in the per-particle list
        self.addPerParticleParameter('symgroup')  # 0 = real, 1..M = copy operator

        n = len(self._uniqueRadii)
        m0Table = self._createUniqueTable(customgbforces.m0)
        d0Table = self._createUniqueTable(customgbforces.d0)
        self.addTabulatedFunction("getd0", Discrete2DFunction(n, n, d0Table))
        self.addTabulatedFunction("getm0", Discrete2DFunction(n, n, m0Table))

        # Symmetry mask table, identical to the nonbonded force's.
        ng = int(self._symmetry_ngroups)
        gtable = symmetry_group_table(ng,
            getattr(self, '_symmetry_exclude_intra', True),
            weights=getattr(self, '_symmetry_group_table', None))
        self.addTabulatedFunction("grouptable", Discrete2DFunction(ng, ng, gtable))
        if nb:
            ngg = self._n_nb_groups
            self._nb_values = [1.0] * (ngg * ngg)
            self._nb_coupling_table = Discrete2DFunction(ngg, ngg, list(self._nb_values))
            self.addTabulatedFunction("nb_coupling_table", self._nb_coupling_table)

        a = self._softcore_params['softcore_a']

        # Descreening integral: UNMASKED by symmetry (option a'); scaled by lam, which
        # is per-group pair_lambda when nb-groups are active (so a decoupled group is
        # descreened proportionally), else the global softcore_lambda.
        self.addComputedValue("I", f"{lam}^{a} * (Ivdw+neckScale*Ineck);"
                                   "Ineck=step(radius1+radius2+neckCut-r)*getm0(radindex1,radindex2)/(1+100*(r-getd0(radindex1,radindex2))^2+"
                                   "0.3*1000000*(r-getd0(radindex1,radindex2))^6);"
                                   "Ivdw=select(step(r+sr2-or1), 0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r), 0);"
                                   "U=r+sr2;"
                                   "L=max(or1, D);"
                                   "D=abs(r-sr2);"
                                   "radius1=or1+offset; radius2=or2+offset;"
                                   "neckScale=0.361825; neckCut=0.68; offset=0.009" + plam,
                                   self.ParticlePairNoExclusions)

        # Born radius: a plain function of I, unmasked.
        self.addComputedValue("B", "1/(1/or-tanh(1.09511284*psi-1.907992938*psi^2+2.50798245*psi^3)/radius);"
                                  "psi=I*or; radius=or+offset; offset=0.009", self.SingleParticle)
        self._createEnergyTerms(self.solventDielectric, self.soluteDielectric,
            self.SA, self.cutoff, self.kappa, self.OFFSET)

    def _createEnergyTerms(self, solventDielectric, soluteDielectric, SA, cutoff, kappa, offset):
        from openmm.app.internal.customgbforces import CustomGBForce
        params = (f'; solventDielectric={solventDielectric:.16g}'
                  f'; soluteDielectric={soluteDielectric:.16g}'
                  f'; kappa={kappa:.16g}; offset={offset:.16g}')
        a = self._softcore_params['softcore_a']
        if cutoff is not None:
            params += f'; cutoff={cutoff:.16g}'
        nb = self._n_nb_groups > 1
        lam = 'pair_lambda' if nb else 'softcore_lambda'          # pairwise softening var
        plam = (';pair_lambda=min(softcore_lambda, nb_coupling_table(nb_group1,nb_group2))'
                if nb else '')
        # Polar self-energy: real atoms only (a copy's is identical to its parent's and
        # would be double-counted on fold-back). Self terms have NO group pair, so they
        # stay on the global softcore_lambda even when per-group coupling is active --
        # their group dependence flows through the (now group-weighted) Born radius B.
        if kappa > 0:
            self.addEnergyTerm(f'step(0.5-symgroup) * softcore_lambda^{a} * -0.5*{ONE_ON_4_PI_EPS0}* (1/soluteDielectric-exp(-kappa*B)/solventDielectric)*charge^2/B' + params,
                CustomGBForce.SingleParticle)
        elif kappa < 0:
            raise ValueError('kappa/ionic strength must be >= 0')
        else:
            self.addEnergyParameterDerivative(f'step(0.5-symgroup) * -0.5*{ONE_ON_4_PI_EPS0}* (1/soluteDielectric-1/solventDielectric)*charge^2/B' + params,
                CustomGBForce.SingleParticle)
        # Nonpolar ACE surface-area self-energy: also real atoms only.
        if SA == 'ACE':
            self.addEnergyTerm(f"step(0.5-symgroup) * softcore_lambda^{a} * 28.3919551* (radius+0.14)^2*(radius/B)^6; radius=or+offset"+params, CustomGBForce.SingleParticle)
        elif SA is not None:
            raise ValueError('Unknown surface area method: '+SA)
        # Pairwise polar term: same-operator copy-copy masked out via grouptable AND
        # scaled by the per-group-pair coupling (lam = pair_lambda when active).
        if cutoff is None:
            if kappa > 0:
                self.addEnergyTerm(f'grouptable(symgroup1,symgroup2) * -{ONE_ON_4_PI_EPS0} * (1/soluteDielectric - exp(-kappa*f)/solventDielectric) * {lam}^({2*a})*charge1*charge2/f;'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params+plam, CustomGBForce.ParticlePairNoExclusions)
            else:
                self.addEnergyTerm(f'grouptable(symgroup1,symgroup2) * -{ONE_ON_4_PI_EPS0}*(1/soluteDielectric - 1/solventDielectric) * {lam}^({2*a})*charge1*charge2/f;'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params+plam, CustomGBForce.ParticlePairNoExclusions)
        else:
            if kappa > 0:
                self.addEnergyTerm(f'grouptable(symgroup1,symgroup2) * -{ONE_ON_4_PI_EPS0} *  (1/soluteDielectric - exp(-kappa*f)/solventDielectric) *{lam}^({2*a}) *charge1*charge2* (1/f - 1/cutoff);'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params+plam, CustomGBForce.ParticlePairNoExclusions
                )
            else:
                self.addEnergyTerm(f'grouptable(symgroup1,symgroup2) * -{ONE_ON_4_PI_EPS0}*(1/soluteDielectric - 1/solventDielectric) * {lam}^({2*a})*charge1*charge2 * (1/f-1/cutoff);'
                    "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"+params+plam, CustomGBForce.ParticlePairNoExclusions
                )

    def _addParticles(self):
        # Mirror the base GBSAGBnForce._addParticles, appending the per-particle group
        # ids after radindex in add-order (self.parameters rows are [charge, or, sr]):
        # nb_group (only when per-group coupling is active) then symgroup. Both
        # self._nb_groups and self._symmetry_groups are aligned to particle order.
        from openmm import CustomGBForce
        nb = self._n_nb_groups > 1
        for i, p in enumerate(self.parameters):
            radIndex = self._radiusToIndex[p[self.RADIUS_ARG_POSITION]]
            extra = [radIndex]
            if nb:
                extra.append(float(self._nb_groups[i]))
            extra.append(float(self._symmetry_groups[i]))
            CustomGBForce.addParticle(self, p + extra)

