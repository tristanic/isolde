# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import numpy
from math import degrees, radians, pi
from . import color, geometry
from chimerax.core.geometry import Places, rotation
from chimerax.core.atomic import Bonds
from chimerax.core.models import Drawing
from simtk.unit import Quantity
from .constants import defaults
from time import time


OPENMM_RADIAL_SPRING_UNIT = defaults.OPENMM_RADIAL_SPRING_UNIT

TWO_PI = 2*pi
FLIP_ON_X = rotation([1,0,0], 180).matrix
Z_AXIS = numpy.array([0,0,1],numpy.double)

class Dihedral():
    '''
    Holds information about the dihedral formed by an array of four atoms.
    The atoms must be an Atom array, ordered so that the torsion axis is
    between the middle two atoms. This axis needn't necessarily correspond to
    a physical bond.
    '''
    
    def __init__(self, atoms, residue = None, name = None):
        self.atoms = atoms
        self._axial_bond = None
        cb = atoms[1:3].intra_bonds
        if len(cb):
            self.axis_bond = cb[0]
        else:
            self.axis_bond = None
        self.residue = residue # the Residue object that this dihedral belongs to
        #self.resnames = self.residues.names
        #self.resnums = self.residues.numbers
        self._rama_case = None
        
        self.name = name # The type of dihedral this is (e.g. phi, psi, omega)
        
        self._target = 0.0
        self._spring_constant = 0.0
        
        # Index of the matching CustomTorsionForce in the simulation so
        # we can restrain this dihedral as needed
        self.sim_index = -1
    
    
    @property
    def axial_bond(self):
        if self._axial_bond is None:
            cb = self.atoms[1:3].intra_bonds
            if len(cb):
                self._axial_bond = cb[0]
            else:
                raise TypeError('The axis of this dihedral is not a bond!')
        return self._axial_bond
    
    @property
    def value(self):
        from .geometry import get_dihedral
        return get_dihedral(*self.atoms.coords)
    
    @property
    def target(self):
        '''
        The target angle to push this dihedral towards, in radians.
        '''
        return self._target
    
    @target.setter
    def target(self, val):
        self._target = val
        
    @property
    def spring_constant(self):
        '''
        The spring constant to restrain this dihedral with, in 
        kJ mol-1 rad-2
        '''
        return self._spring_constant
    
    @spring_constant.setter
    def spring_constant(self, k):
        if isinstance(k, Quantity):
            k = k.value_in_unit(OPENMM_RADIAL_SPRING_UNIT)
        self._spring_constant = k
    
    
    @property
    def coords(self):
        return self.atoms.coords
    
    @property
    def names(self):
        return self.atoms.names
    
    @property
    def residues(self):
        # List of four residues corresponding to the atoms in the dihedral
        return self.atoms.residues
    
    @property
    def rama_case(self):
        return self._rama_case
    
    @rama_case.setter
    def rama_case(self, name):
        from . import validation
        if name not in validation.RAMA_CASES:
            errstring = 'Invalid Ramachandran case! Must be one of {}'.format(validation.RAMA_CASES)
            raise Exception(errstring)
        self._rama_case = name

class Dihedrals():
    '''
    Holds an array of Dihedral objects and their statistics, arranged for
    fast look-up and handling.
    '''
    def __init__(self, dlist = None, drawing = None, session = None):
        self.session = session
        self._dihedrals = []
        # ChimeraX Residues object in the same order as the list of dihedrals
        self._residues = None
        # ChimeraX Atoms object holding all dihedral atoms. The atoms 
        # corresponding to dihedral i will be found at [4*i:4*i+4]
        self._atoms = None
        # A ChimeraX Bonds object made up of the axial bond for each 
        # dihedral.
        self._axial_bonds = None
        # Current dihedral values
        self._values = None
        # Array of Ramachandran cases (if applicable to this dihedral type)
        self._rama_cases = None
        if dlist is not None:
            self._dihedrals = dlist
        self._names = None
        # ChimeraX Drawing object to draw restraint annotations into 
        # if desired
        self.colors = None
        self.set_target_drawing(drawing)
        self.color_scale = color.standard_three_color_scale('GYPi', 0, pi/2, radians(30))
        self._restrained_bonds = Bonds()
        self._restrained_dihedrals = None
        self.update_needed = True
    
    def set_target_drawing(self, drawing):
        d = self._drawing = drawing
        if d is not None:
            r = self._ring_drawing = Drawing('rings')
            p = self._post_drawing = Drawing('posts')
            d.add_drawing(r)
            d.add_drawing(p)
            r.vertices, r.normals, r.triangles = geometry.ring_arrow_with_post(0.5, 0.05, 3, 6, 0.25, 0.1, 0.05, 1)
            p.vertices, p.normals, p.triangles = geometry.post_geometry(0.05, 1, caps=False)
            r.positions = Places(places=[])
            p.positions = Places(places=[])
            if self.session is not None:
                from chimerax.core.atomic import get_triggers
                triggers = get_triggers(self.session)
                self._handler = triggers.add_handler('changes', self._update_graphics)

    
    def _update_graphics(self, trigger_name, changes):
        reasons = changes.atom_reasons()
        dchanged = 'display changed' in reasons
        cchanged = 'coord changed' in reasons
        if dchanged or cchanged:
            self.update_graphics()
    
    def update_graphics(self, *_, update_needed = None):
        if update_needed is None:
            update_needed = self.update_needed
        d = self._drawing
        if d is None:
            return
        r = self._ring_drawing
        p = self._post_drawing
        if update_needed:
            self._update_restrained_bonds()
        r_d = self._restrained_dihedrals
        if r_d is None or len(r_d) == 0:
            if r.display or p.display:
                r.display = False
                p.display = False
            #~ r.positions = Places([])
            #~ p.positions = Places([])
            return
        r.display = True
        p.display = True
        print('Displaying dihedral annotations')
        rb = self._restrained_bonds
        showns = rb.showns
        shown_indices = numpy.where(showns)[0]
        r_d = r_d[shown_indices]
        b = rb[showns]
        if len(b) == 0:
            r.positions = Places([])
            p.positions = Places([])
            return
        targets = self._cached_targets[shown_indices]
        offsets = (r_d.values - targets + pi) % TWO_PI - pi 
        flip_mask = offsets<0
        rotations = geometry.rotations(Z_AXIS, offsets)
        shifts = geometry.bond_cylinder_placements(b)
        tf = geometry.flip_rotate_and_shift(flip_mask, FLIP_ON_X,
                rotations, shifts.array())
        r.positions = Places(place_array = tf)
        p.positions = shifts
        colors = self.color_scale.get_colors(numpy.abs(offsets))
        r.colors = colors
        p.colors = colors
    
    
    def _update_restrained_bonds(self):
        restrained_bond_mask = (self.spring_constants > 0)
        indices = numpy.where(restrained_bond_mask)[0]
        self._restrained_bonds = self.axial_bonds[restrained_bond_mask]
        rd = self._restrained_dihedrals = self[indices]
        self._cached_targets = rd.targets
        
        self.update_needed = False
    
    @property
    def axial_bonds(self):
        if self._axial_bonds is None:
            try:
                from chimerax.core.atomic import Bonds
                self._axial_bonds = Bonds([d.axial_bond for d in self])
            except:
                raise TypeError('At least one dihedral in this array '\
                    +'has no axial bond!')
        return self._axial_bonds
    
    def __len__(self):
        return len(self.dihedrals)
    
    def __bool__(self):
        return len(self) > 0
    
    def __iter__(self):
        return iter(self.dihedrals)
    
    def __getitem__(self, i):
        if isinstance(i,(int, numpy.integer)):
            return self.dihedrals[i]
        elif isinstance(i,(slice)):
            return Dihedrals(dlist=self.dihedrals[i])
        elif isinstance(i, numpy.ndarray):
            return Dihedrals([self.dihedrals[j] for j in i])
        else:
            raise IndexError('Only integer indices allowed for %s, got %s'
                % (self.__class__.__name__, str(type(i))))
    
    def by_residue(self, residue):
        from chimerax.core.atomic import Residues
        i = self.residues.index(residue)
        if i == -1:
            return None
        return self[i]
    
    def by_residues(self, residues):
        indices = self.residues.indices(residues)
        indices = indices[indices != -1]
        if len(indices):
            return self[indices]
        return None
   
    def append(self, d):
        from chimerax.core.atomic import concatenate, Residues
        self._axial_bonds = None
        if isinstance(d, Dihedral):
            if self.atoms is None:
                self._dihedrals = [d]
            else:
                self.dihedrals.append(d)
                self._residues = concatenate([self.residues, [d.residue]])
                self._atoms = concatenate([self.atoms, d.atoms])
        elif isinstance(d, Dihedrals):
            if self.atoms is None:
                self._dihedrals = d.dihedrals
                self._residues = d.residues
                self._atoms = d.atoms
            elif len(d):
                self._dihedrals.extend(d.dihedrals)
                self._residues = concatenate([self.residues, d.residues])
                self._atoms = concatenate([self.atoms, d.atoms])
        else:
            raise TypeError('Can only append a single Dihedral or a Dihedrals object.')    
    
    def index(self, d):
        try:
            i = self._dihedrals.index(d)
        except ValueError:
            return -1
        return i
    
    def indices(self, d_list):
        return numpy.array([self.index(d) for d in d_list])
        
    @property
    def atoms(self):
        if not len(self):
            return None
        if self._atoms is None:
            atoms = numpy.empty([len(self),4],dtype='object')
            for i, d in enumerate(self):
                atoms[i] = d.atoms
            from chimerax.core.atomic import Atoms
            self._atoms = Atoms(numpy.ravel(atoms))
        return self._atoms    
    
    @property
    def axis_bonds(self):
        from chimerax.core.atomic import Bonds
        return Bonds([d.axis_bond for d in self])
    
    @property
    def residues(self):
        if self._residues is None:
            residues = [d.residue for d in self.dihedrals] 
            from chimerax.core.atomic import Residues
            self._residues = Residues(residues)            
        return self._residues
        
    
    @property
    def targets(self):
        return numpy.array([d.target for d in self], numpy.float32)
    
    @targets.setter
    def targets(self, t_or_t_array):
        if hasattr(t_or_t_array, '__len__'):
            for d, t in zip(self, t_or_t_array):
                d.target = t
        else:
            for d in self:
                d.target = t_or_t_array
    
    @property
    def names(self):
        if self._names is None:
            self._names = numpy.array([d.name for d in self])
        return self._names
        
                
    
    @property
    def spring_constants(self):
        return numpy.array([d.spring_constant for d in self], numpy.float32)
    
    @spring_constants.setter
    def spring_constants(self, k_or_k_array):
        if isinstance(k_or_k_array, Quantity):
            k_or_k_array = k_or_k_array.value_in_unit(OPENMM_RADIAL_SPRING_UNIT)
        if hasattr(k_or_k_array, '__len__'):
            for d, k in zip(self, k_or_k_array):
                d.spring_constant = k
        else:
            for d in self:
                d.spring_constant = k_or_k_array
    @property
    def coords(self):
        return self.atoms.coords
    
    @property
    def dihedrals(self):
        return self._dihedrals
        
    @property
    def values(self):
        from . import geometry
        return geometry.get_dihedrals(self.coords, len(self))
    
    @property
    def rama_cases(self):
        # Need to check every time in case any have been changed directly
        cases = self._rama_cases = [d.rama_case for d in self]
        return cases
    
    @rama_cases.setter
    def rama_cases(self, cases):
        if not hasattr(cases,'__len__') or len(cases) != len(self):
            raise TypeError('Array size must match number of dihedrals!')
        self._rama_cases = cases
        for d, case in zip(self, cases):
            d.rama_case = case
            
        

class Backbone_Dihedrals():
    '''
    Takes either an atomic model or three Dihedrals objects with pre-defined
    phi, psi and omega dihedrals. If provided with a model, it will find and
    store all phi, psi and omega dihedrals as Dihedrals objects.
    '''
    
    # (phi, psi) angles for common secondary structures
    standard_phi_psi_angles = {
        'helix': (radians(-64.0), radians(-41.0)),
        'parallel beta': (radians(-119), radians(113)),
        'antiparallel beta': (radians(-139.0), radians(135.0)),
        }
        
    
    
    def __init__(self, session, model = None, phi = None, psi = None, omega = None, old = False):
        if model == None and (phi == None or psi == None or omega == None):
            raise TypeError('You must provide either a model or all three of\
                            phi, psi and omega!')
        elif model and (phi or psi or omega):
            raise TypeError('Cannot provide both a model and predefined dihedrals!')
        elif model and not model.atomspec_has_atoms():
            raise TypeError('Please provide a model containing atoms!')
        
        self.session = session
        # It's most convenient to determine and store the Ramachandran case
        # for each residue here, otherwise things start to get messy when
        # working with subsets. We'll also keep a list of all proline indices
        # since we have to double-check their peptide bond state and re-categorise
        # as necessary
        self._rama_cases = None
        # Numpy array holding the current Ramachandran scores
        self._rama_scores = None
        # Colours for C-alpha atoms will be changed according to Ramachandran score
        self._rama_colors = None
        
        self._phi_vals = None
        self._psi_vals = None
        self._omega_vals = None
        
        if model:
            self.residues = model.residues
            # Filter to get only amino acid residues
            from chimerax.core.atomic import Residue
            f = self.residues.polymer_types == Residue.PT_AMINO
            self.residues = self.residues.filter(f)
            self.residues = self.residues[numpy.lexsort((self.residues.numbers, self.residues.chains.chain_ids))]
            self.resnames = self.residues.names
            self.atoms = self.residues.atoms
            # We want to colour C-alphas according to their status, so we'll
            # hold an array of them here.
            self.CAs = self.atoms.filter(self.atoms.names == 'CA')
            nr = self.num_residues = len(self.residues)
            # Empty lists to fill with phi, psi and omega. Once filled, they
            # will be converted to Dihedrals objects
            self.phi = []
            self.psi = []
            self.omega = []
            # We want to keep a single master list of residues, but not all
            # residues will have all three dihedrals. So, we'll keep an array
            # of indices for each type, mapping the dihedral index to the
            # residue index.
            self._phi_indices = []
            self._psi_indices = []
            self._omega_indices = []
            if old:
                self.find_dihedrals_old()
                return
            self.find_dihedrals()
        else:
            self.phi = phi
            self.psi = psi
            self.omega = omega
            
            self.residues = phi.residues.merge(psi.residues).merge(omega.residues)
            residues = self.residues = self.residues[numpy.lexsort((self.residues.numbers, self.residues.chains.chain_ids))]
            self.resnames = residues.names
            atoms = self.residues.atoms
            self.CAs = atoms.filter(atoms.names == 'CA')
            phi_indices = []
            psi_indices = []
            omega_indices = []
            phi_indices = residues.indices(phi.residues)
            psi_indices = residues.indices(psi.residues)
            omega_indices = residues.indices(omega.residues)
            
            self._phi_indices = numpy.array(phi_indices[phi_indices != -1],numpy.int32)
            self._psi_indices = numpy.array(psi_indices[psi_indices != -1],numpy.int32)
            self._omega_indices = numpy.array(omega_indices[omega_indices != -1],numpy.int32)
            
            # Get a list of Ramachandran cases to sort into a dict. We need
            # to merge the arrays for Phi and Psi to make sure we get them
            # all.
            
            rama_cases = numpy.array([None]*len(self.residues))
            rama_cases[self._phi_indices] = phi.rama_cases
            rama_cases[self._psi_indices] = psi.rama_cases
            
            from . import validation
            self._rama_cases = {}
            for key in validation.RAMA_CASES:
                self._rama_cases[key] = numpy.where(rama_cases == key)[0]
                
            
    @property
    def all_vals(self):
        return self.phi_vals, self.psi_vals, self.omega_vals
    
    @property
    def phi_vals(self):
        if self._phi_vals is None:
            self._phi_vals = numpy.ones(len(self.residues),numpy.float32)*numpy.nan
        self._phi_vals[self._phi_indices] = self.phi.values
        return self._phi_vals
    
    @property
    def psi_vals(self):
        if self._psi_vals is None:
            self._psi_vals = numpy.ones(len(self.residues),numpy.float32)*numpy.nan
        self._psi_vals[self._psi_indices] = self.psi.values
        return self._psi_vals
    
    @property    
    def omega_vals(self):
        if self._omega_vals is None:
            self._omega_vals = numpy.ones(len(self.residues),numpy.float32)*numpy.nan
        self._omega_vals[self._omega_indices] = self.omega.values
        return self._omega_vals
    
    @property
    def rama_cases(self):
        self.update_pro_rama_cases(self.omega_vals)
        return self._rama_cases
    
    @property
    def rama_scores(self):
        if self._rama_scores is None:
            self._rama_scores = numpy.ones(len(self.residues),numpy.float32) * -1
        return self._rama_scores
    
    @property
    def rama_colors(self):
        if self._rama_colors is None:
            self._rama_colors = numpy.array([[128,128,128,255]]*len(self.residues),numpy.uint8)
        return self._rama_colors
        
    def by_residue(self, res):
        ''' 
        Return the phi, psi and omega dihedrals for a given residue
        '''
        phi = self.phi.by_residue(res)
        psi = self.psi.by_residue(res)
        omega = self.omega.by_residue(res)
        return phi, psi, omega
    
    def by_residues(self, reslist):
        '''
        Return phi, psi and omega dihedrals for all residues in a 
        Residues array.
        '''
        phi = self.phi.by_residues(reslist)
        psi = self.psi.by_residues(reslist)
        omega = self.omega.by_residues(reslist)
        return phi, psi, omega
    
    def update_pro_rama_cases(self, omega_vals):
        rc = self._rama_cases
        current_trans = rc['TransPro']
        current_cis = rc['CisPro']
        switch_to_cis = []
        switch_to_trans = []
        from . import validation
        for i in current_trans:
            if validation.omega_type(omega_vals[i]) == 'cis':
                switch_to_cis.append(i)
        for i in current_cis:
            if validation.omega_type(omega_vals[i]) == 'trans':
                switch_to_trans.append(i)
        
        if len(switch_to_cis):
            switch_to_cis = numpy.array(switch_to_cis)
            rc['TransPro'] = current_trans[numpy.in1d(current_trans, switch_to_cis, invert=True)]
            rc['CisPro'] = numpy.append(current_cis, switch_to_cis)
            
            switched_residues = self.residues[switch_to_cis]
            for r in switched_residues:
                for d in self.by_residue(r):
                    if d is not None:
                        d.rama_case = 'CisPro'
            
        if len(switch_to_trans):
            switch_to_trans = numpy.array(switch_to_trans)
            rc['CisPro'] = current_cis[numpy.in1d(current_cis, switch_to_trans, invert=True)]
            rc['TransPro'] = numpy.append(current_trans, switch_to_trans)
            
            switched_residues = self.residues[switch_to_trans]
            for r in switched_residues:
                for d in self.by_residue(r):
                    if d is not None:
                        d.rama_case = 'TransPro'
            
         
    def find_dihedrals_old(self):
        '''
        Old, unused version. New code is about fifteen times faster.
        '''
        if len(self.phi) or len(self.psi) or len(self.omega):
            import warnings
            warnings.warn('Backbone dihedrals have already been defined. \
                           If you want to update them, create a new \
                           Backbone_Dihedrals object.')
            return            
        bond_to_last = False
        bond_to_next = False
        last_residue = None
        last_atoms = None
        last_names = None
        from copy import copy
        # Loop through all residues in the selection held by the object,
        # picking out amino acid residues and storing Atom arrays defining
        # their phi, psi and omega atoms where applicable. To avoid 
        # unforseen errors, we'll explicitly check connectivity for the
        # N- and C-termini of every residue.
        for i, r in enumerate(self.residues):
            if not r.polymer_type == r.PT_AMINO:
                continue
            a = r.atoms
            names = a.names
            # Build the lists of dihedral atoms found in this residue
            phi_atoms = []
            for name in ['N','CA', 'C']:
                phi_atoms.append(a.filter(numpy.in1d(a.names, name)))
            # These three atoms are common to both phi and psi, so we can
            # just copy the array
            psi_atoms = copy(phi_atoms)
            omega_atoms = phi_atoms[0:2]
            # Atoms from the previous residue
            # Check to see if this residue is bonded to a previous
            # one.
            N = phi_atoms[0][0]
            N_bonded_atom_list = N.bonds.atoms[0].merge(N.bonds.atoms[1])
            prev_C_l = N_bonded_atom_list.filter(numpy.in1d(N_bonded_atom_list.names, 'C'))
            if len(prev_C_l):
                last_residue = prev_C_l[0].residue
                if not last_residue.PT_AMINO:
                    last_residue = None
                    bond_to_last = False
                else:
                    bond_to_last = True
                    last_atoms = last_residue.atoms
                    last_names = last_atoms.names
            else:
                bond_to_last = False
            if bond_to_last:
                C = last_atoms.filter(numpy.in1d(last_names, 'C'))
                phi_atoms.insert(0, C)
                omega_atoms.insert(0, C)
                CA = last_atoms.filter(numpy.in1d(last_names, 'CA'))
                omega_atoms.insert(0, CA)
            else:
                # N-terminal residues don't have phi or omega dihedrals
                phi_atoms = None
                omega_atoms = None
            
            # Atoms from the next residue
            C = psi_atoms[-1][0]
            C_bonded_atom_list = C.bonds.atoms[0].merge(C.bonds.atoms[1])
            next_N_l = C_bonded_atom_list.filter(numpy.in1d(C_bonded_atom_list.names, 'N'))
            if len(next_N_l):
                next_residue = next_N_l[0].residue
                if next_residue.PT_AMINO:
                    bond_to_next = True
                else:
                    bond_to_next = False
            else:
                bond_to_next = False
                
            if bond_to_next:
                next_atoms = next_residue.atoms
                next_names = next_atoms.names
                CA = next_atoms.filter(numpy.in1d(next_names, 'N'))
                psi_atoms.append(CA)
            else:
                # C-terminal residues have no psi dihedral
                psi_atoms = None
            from chimerax.core.atomic import concatenate, Atoms
            if phi_atoms is not None:
                self.phi.append(Dihedral(concatenate(phi_atoms, Atoms), residue = r))
                self._phi_indices.append(i)
            if psi_atoms is not None:
                self.psi.append(Dihedral(concatenate(psi_atoms, Atoms), residue = r))
                self._psi_indices.append(i)
            if omega_atoms is not None:
                self.omega.append(Dihedral(concatenate(omega_atoms, Atoms), residue = r))
                self._omega_indices.append(i)
        
        # Convert to Dihedrals objects
        self.phi = Dihedrals(self.phi)
        self.psi = Dihedrals(self.psi)
        self.omega = Dihedrals(self.omega)
        
        self._phi_indices = numpy.array(self._phi_indices)
        self._psi_indices = numpy.array(self._psi_indices)
        self._omega_indices = numpy.array(self._omega_indices)            
                    
                        
                
            
    def find_dihedrals(self):
        '''
        Identifies all protein backbone phi, psi and omega dihedrals in
        a model, generates a Dihedrals containing each set, and sorts the
        residues into the MolProbity Ramachandran groupings.
        '''
        if len(self.phi) or len(self.psi) or len(self.omega):
            import warnings
            warnings.warn('Backbone dihedrals have already been defined. \
                           If you want to update them, create a new \
                           Backbone_Dihedrals object.')
            return            
        from chimerax.core.atomic import Residue
        
        self.session.ui.processEvents()
        
        # Get all protein residues and their atoms    
        res = self.residues = self.residues.filter(
                        self.residues.polymer_types == Residue.PT_AMINO)
        resnums = res.numbers
        atoms = res.atoms
        
        # Get all N, CA, C in ordered arrays. We need to hold the CA atoms
        # long-term for visualisation purposes.
        N_atoms = atoms.filter(atoms.names == 'N')
        CA_atoms = self.CAs = atoms.filter(atoms.names == 'CA')
        C_atoms = atoms.filter(atoms.names == 'C')
        
        self.session.ui.processEvents()
        
        # Get all the C-N bonds
        CN_atoms = atoms.filter(numpy.any(numpy.column_stack(
            [atoms.names == 'N', atoms.names == 'C']), axis = 1))
            
        CN_bonds = CN_atoms.intra_bonds
        
        bonded_C = CN_bonds.atoms[0]
        bonded_N = CN_bonds.atoms[1]
        
        bonded_C_indices = atoms.indices(bonded_C)
        bonded_N_indices = atoms.indices(bonded_N)
        # usually this should be the way things work, but occassionally
        # we end up with atoms in the wrong array. So let's catch and
        # sort that out.
        if not numpy.all(bonded_C.names == 'C'):
            bad_indices = numpy.argwhere(bonded_C.names != 'C')
            for i in bad_indices:
                thisN = bonded_C_indices[i]
                bonded_C_indices[i] = bonded_N_indices[i]
                bonded_N_indices[i] = thisN
            bonded_C = atoms[bonded_C_indices]
            bonded_N = atoms[bonded_N_indices]
        bonded_C_indices = C_atoms.indices(bonded_C)
        bonded_N_indices = N_atoms.indices(bonded_N)
        bonded_C_resnames = bonded_C.unique_residues.names
        assert(numpy.all(bonded_N.names == 'N'))
        assert(numpy.all(bonded_C.names == 'C'))
        bonded_N_resnames = bonded_N.unique_residues.names

        # We also need the CA atom from the preceding residue to make up
        # the omega dihedral
        bonded_C_residues = bonded_C.residues
        prev_atoms = bonded_C_residues.atoms
        prev_CA = prev_atoms.filter(prev_atoms.names == 'CA')
        
        self.session.ui.processEvents()
        
        '''
        Build up a 6 * (number of residues) numpy array where each row is
        [CA, C, N, CA, C, N].
        The omega dihedral is entries 0 to 3, phi is 1 to 4, psi is 2 to 5.
        '''
        total_n_res = len(res)
        
        master_array = numpy.array([[None] * 6] * total_n_res)
        master_array[:,2] = N_atoms
        master_array[:,3] = CA_atoms
        master_array[:,4] = C_atoms
        
        master_array[bonded_C_indices+1,1] = bonded_C
        master_array[bonded_C_indices+1,0] = prev_CA
        master_array[bonded_N_indices-1,5] = bonded_N
        
        
        ome_i = self._omega_indices = numpy.where(numpy.all(master_array[:,0:4], axis=1))[0]
        phi_i = self._phi_indices = numpy.where(numpy.all(master_array[:,1:5], axis = 1))[0]
        psi_i = self._psi_indices = numpy.where(numpy.all(master_array[:,2:6], axis = 1))[0]
                
        raw_omegas = master_array[ome_i,0:4]
        raw_phis = master_array[phi_i,1:5]
        raw_psis = master_array[psi_i,2:6]
        
        from chimerax.core.atomic import Atoms
        
        # Create all the Dihedrals arrays
        omega_1d = numpy.ravel(raw_omegas)

        omega_res = res[ome_i]
        omega = []
        for i in range(len(raw_omegas)):
            omega.append(Dihedral(Atoms(omega_1d[4*i:4*i+4]), omega_res[i], name = 'omega'))
        omega = self.omega = Dihedrals(omega)
        
        phi_1d = numpy.ravel(raw_phis)
        phi_res = res[phi_i]
        phi = []
        for i in range(len(raw_phis)):
            phi.append(Dihedral(Atoms(phi_1d[4*i:4*i+4]), phi_res[i], name = 'phi'))
        phi = self.phi = Dihedrals(phi)
        
        psi_1d = numpy.ravel(raw_psis)
        psi_res = res[psi_i]
        psi = []
        for i in range(len(raw_psis)):
            psi.append(Dihedral(Atoms(psi_1d[4*i:4*i+4]), psi_res[i], name = 'psi'))
        psi = self.psi = Dihedrals(psi)
        


        # To determine the MolProbity Ramachandran case for each residue we
        # need to know both its name and the name of the following residue.
        # We also need the value of the omega dihedral to distinguish cis-Pro
        # from trans-Pro. Only residues that have both omega and phi dihedrals
        # count towards Ramachandran statistics.
        
        has_phi = numpy.array([False]*total_n_res)
        has_phi[phi_i] = True
        
        has_psi = numpy.array([False]*total_n_res)
        has_psi[psi_i] = True
        
        counts_for_rama = numpy.all([has_phi, has_psi], axis=0)
          
        first_res = res[psi_i]
        next_res = Atoms(raw_psis[:,3]).unique_residues
        
        # Column 1: The residue to which the Ramachandran score will apply
        # Column 2: The residue to which residue 1's C is bonded
        rama_resnames = numpy.array([[None]*2]*total_n_res)
        rama_resnames[psi_i,0] = first_res.names
        rama_resnames[psi_i,1] = next_res.names
        
        omega_vals = numpy.array([numpy.nan]*total_n_res)
        omega_vals[ome_i] = omega.values
        
        from . import validation
        
        self._rama_cases, rca = validation.sort_into_rama_cases(
                                counts_for_rama, rama_resnames, omega_vals)
        omega_rama = rca[ome_i]
        phi_rama = rca[phi_i]
        psi_rama = rca[psi_i]
        
        omega.rama_cases = omega_rama
        phi.rama_cases = phi_rama
        psi.rama_cases = psi_rama
        
        
        
        
        
        
        
       
        
    
        
        
        
