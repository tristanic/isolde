# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


# Custom OpenMM forces

import numpy
from math import pi, radians, degrees, cos
from simtk import unit, openmm
from simtk.openmm import app
from simtk.openmm.openmm import CustomBondForce, CustomExternalForce, \
                                CustomCompoundBondForce, CustomTorsionForce, \
                                NonbondedForce, CMAPTorsionForce
from simtk.openmm.openmm import Continuous1DFunction, Continuous3DFunction, \
                                Discrete3DFunction
from simtk.openmm.app.internal import customgbforces   
from . import amber_cmap


MIN_K = 0.01 # Spring constants below this value will be treated as zero



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
    
    def addTorsion(self, resname, phi_indices, psi_indices):
        map_index = self._map_loader.map_index(resname)
        super().addTorsion(map_index, *phi_indices.tolist(), *psi_indices.tolist())


class LinearInterpMapForce(CustomCompoundBondForce):
    '''
    Converts a map of (u,v,w) data and a (x,y,z)->(u,v,w) transformation
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
        if units == 'angstroms':
            # OpenMM calcs will be in nm
            tf[:,0:3] *= 10
        elif units != 'nanometers':
            raise TypeError('Units must be either "angstroms" or "nanometers"!')
        map_func = self._discrete3D_from_volume(data)
        
        eps = 1e-6
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
                if abs(val) > eps:
                    count += 1
                    tf_strings[i] += spacer + entry.format(val)
            
        i_str, j_str, k_str = tf_strings
        
        
        
        #~ i_str = 'i = x1* {} + y1 * {} + z1 * {} + {}'.format(
            #~ tf[0][0], tf[0][1], tf[0][2], tf[0][3])
        #~ j_str = 'j = x1* {} + y1 * {} + z1 * {} + {}'.format(
            #~ tf[1][0], tf[1][1], tf[1][2], tf[1][3])
        #~ k_str = 'k = x1* {} + y1 * {} + z1 * {} + {}'.format(
            #~ tf[2][0], tf[2][1], tf[2][2], tf[2][3])
        
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
        self.setGlobalParameterDefaultValue(self._global_k_index, k)
        self.update_needed = True
    
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
        #super().__init__('min(0.5*k*(r-r0)^2, max_force*abs(r-r0))')
        linear_eqn = 'max_force * abs(r-r0)'# - 0.5*max_force^2/k'
        quadratic_eqn = '0.5*k*(r-r0)^2'
        transition_eqn = 'step(r - max_force/k)'
        zero_k_eqn = 'step(min_k - k)'
        energy_str = 'select({},0,select({},{},{}))'.format(
            zero_k_eqn, transition_eqn, linear_eqn, quadratic_eqn)
        #force_str = 'select(' + ','.join((transition_eqn, linear_eqn, quadratic_eqn)) + ')'
        super().__init__(energy_str)
        self._max_force = max_force
        
        self.k_index = self.addPerBondParameter('k')
        self.r0_index = self.addPerBondParameter('r0')
        self.max_force_index = self.addGlobalParameter('max_force', self.max_force)
        self.min_k_index = self.addGlobalParameter('min_k', MIN_K)
        self.update_needed = False
    
    @property
    def max_force(self):
        '''Maximum force applied to any given atom, in kJ/mol/nm.'''
        return self._max_force
    
    @max_force.setter
    def max_force(self, force):
        self.setGlobalParameterDefaultValue(self.max_force_index, force)
        self._max_force = force
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
        #super().__init__('min(0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2), max_force *(abs(x-x0)+abs(y-y0)+abs(z-z0)))')
        linear_eqn = 'max_force * sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)'# - 0.5*max_force^2/k'
        quadratic_eqn = '0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)'
        transition_eqn = 'step(sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2) - max_force/k)'
        zero_k_eqn = 'step(min_k - k)'
        energy_str = 'select({},0,select({},{},{}))'.format(
            zero_k_eqn, transition_eqn, linear_eqn, quadratic_eqn)
        
        #force_str = 'select(' + ','.join((transition_eqn, linear_eqn, quadratic_eqn)) + ')'
        super().__init__(energy_str)
        self._max_force = max_force
        per_particle_parameters = ('k','x0','y0','z0')
        for p in per_particle_parameters:
            self.addPerParticleParameter(p)
        self.addGlobalParameter('max_force', max_force)
        self.addGlobalParameter('min_k', MIN_K)

        self.update_needed = False
    
    @property
    def max_force(self):
        '''Maximum force applied to any given atom, in kJ/mol/nm.'''
        return self._max_force
    
    @max_force.setter
    def max_force(self, force):
        self.setGlobalParameterDefaultValue(0, force)
        self._max_force = force
        self.update_needed = True
    
        

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

        self.update_needed = True

    def update_target(self, index, target = None, k = None, cutoff = None):
        current_params = self.getTorsionParameters(index)
        indices = current_params[0:4]
        new_k, new_theta0, new_cutoff = current_params[4]
        if target is not None:
            new_theta0 = target
        if k is not None:
            new_k = k
        if cutoff is not None:
            new_cutoff = cos(cutoff)
        self.setTorsionParameters(index, *indices, (new_k, new_theta0, new_cutoff))
        self.update_needed = True



#~ class FlatBottomTorsionRestraintForce(CustomTorsionForce):
    #~ '''
    #~ Wraps an OpenMM CustomTorsionForce to restrain torsion angles while 
    #~ allowing free movement within a range (target +/- cutoff). Within 
    #~ the cutoff range the potential is constant, while outside it 
    #~ is defined as -k * cos (theta-theta0). 
    #~ '''
    #~ def __init__(self):
        #~ normalize_fn = '((theta-theta0 + pi) - floor((theta-theta0 + pi) / two_pi) * two_pi - pi)' 
        #~ standard_energy = '0.5*k * dtheta^2'
        #~ flat_energy = '0.5* k * cutoff^2'
        #~ switch_function = 'step(cutoff-dtheta)'
        #~ complete_function = 'select({},{},{});dtheta = {}'.format(
            #~ switch_function, flat_energy, standard_energy, normalize_fn)
        #~ super().__init__(complete_function)
        #~ per_bond_parameters = ('k', 'theta0', 'cutoff')
        #~ for p in per_bond_parameters:
            #~ self.addPerTorsionParameter(p)
        #~ self.addGlobalParameter('pi', pi)
        #~ self.addGlobalParameter('two_pi', 2*pi)
        
        #~ self.update_needed = False
    
    def set_cutoff_angle(self, i, angle):
        '''
        Set the cut-off angle (below which no force will be applied) in 
        radians for one dihedral.
        '''
        p1, p2, p3, p4, current_params = self.getTorsionParameters(i)
        current_params[2] = angle
        self.setTorsionParameters(i, p1, p2, p3, p4, current_params)
        self.update_needed = True
    
    def get_cutoff_angle(self, i):
        '''
        Get the cut-off angle in radians for one dihedral.
        '''
        return self.getBondParameters(i)[5][2]
    
class GBSAForce(customgbforces.GBSAGBn2Force):
    def __init__(self):
        '''
        kappa = 3.0 --> approx. 0.5M ion concentration at 100K
        '''
        super().__init__(solventDielectric=78.5, soluteDielectric=1,
                        SA = 'ACE', cutoff = 1.0, kappa = 3.0)
        self.setNonbondedMethod(openmm.NonbondedForce.CutoffNonPeriodic)

