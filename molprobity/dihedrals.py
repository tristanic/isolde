
class Dihedral():
    '''
    Holds information about the dihedral formed by an array of four atoms.
    The atoms must be an Atoms object, ordered so that the torsion axis is
    between the middle two atoms. This axis needn't necessarily correspond to
    a physical bond.
    '''
    def __init__(self, atoms, residue = None):
        self.atoms = atoms
        self._residue = residue 
        self._rama_case = None # Ramachandran case (if applicable)
        
    
    @property
    def value(self):
        '''The dihedral angle in radians.'''
        from .geometry import get_dihedral
        return get_dihedral(*self.atoms.coords)
    
    @property
    def coords(self):
        return self.atoms.coords
    
    @property
    def names(self):
        return self.atoms.names
    
    @property
    def residue(self):
        '''
        The Residue object that this dihedral belongs to. If the atoms
        in the dihedral span more than one residue, this needs to be 
        set manually.
        '''
        if self._residue is None:
            if len(set(self.residues)) == 1:
                self._residue = self.residues[0]
            else:
                errmsg = '''
                This dihedral involves atoms from more than one residue,
                and therefore cannot be assigned automatically to a 
                specific residue. You can specify the residue a dihedral
                belongs to by either using the 'residue = Residue' argument on
                initialisation, or by '{dihedral}.residue = Residue' on an
                existing object.
                '''
                raise RuntimeError(errmsg)
        return self._residue
    
    @residue.setter
    def residue (self, residue):
        from chimerax.core.atomic.molobject import Residue
        if type(residue) == Residue:
            self._residue = residue
        else:
            raise TypeError('residue must be a ChimeraX Residue object!')
            
    
    @property
    def residues(self):
        '''List of four residues corresponding to the atoms in the dihedral.'''
        return self.atoms.residues
    
    
    @property
    def rama_case(self):
        ''' (Only applicable if this is a phi or psi dihedral)
        The type of Ramachandran distribution this dihedral belongs to
        in terms of the standard MolProbity cases. N-terminal and 
        C-terminal residues are inapplicable to the Ramachandran plot,
        and are sorted into special cases.
        Valid entries:
            'Nter', 'Cter', 'CisPro', 'TransPro', 'Glycine', 'PrePro', 
            'IleVal', 'General'
        '''
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
    def __init__(self, dlist = None):
        self._dihedrals = None
        # ChimeraX Residues object in the same order as the list of dihedrals
        self._residues = None
        # ChimeraX Atoms object holding all dihedral atoms. The atoms 
        # corresponding to dihedral i will be found at [4*i:4*i+4]
        self._atoms = None
        # Current dihedral values
        self._values = None
        # Array of Ramachandran cases (if applicable to this dihedral type)
        self._rama_cases = None
        if dlist is not None:
            self._dihedrals = dlist
    
    def __len__(self):
        return len(self._dihedrals)
    
    def __bool__(self):
        return len(self) > 0
    
    def __iter__(self):
        return iter(self._dihedrals)
    
    def __getitem__(self, i):
        import numpy
        if isinstance(i,(int,numpy.int32)):
            return self._dihedrals[i]
        elif isinstance(i,(slice)):
            return Dihedrals(self._dihedrals[i])
        elif isinstance(i, numpy.ndarray):
            return Dihedrals([self._dihedrals[j] for j in i])
        else:
            raise IndexError('Only integer indices allowed for %s, got %s'
                % (self.__class__.__name__, str(type(i))))
    
    def by_residue(self, residue):
        '''
        Get the dihedral value for a given residue. Argument is a 
        ChimeraX Residue object.
        '''
        from chimerax.core.atomic import Residues
        i = self.residues.index(residue)
        if i == -1:
            return None
        return self[i]
    
    def by_residues(self, residues):
        '''
        Returns a Dihedrals object containing all dihedrals corresponding
        to the residues in the given Residues object.
        '''
        indices = self.residues.indices(residues)
        indices = indices[indices != -1]
        if len(indices):
            return self[indices]
        return None
   
    def append(self, d):
        '''
        Add one or more dihedrals to this object.
        '''
        from chimerax.core.atomic import concatenate, Residues
        if isinstance(d, Dihedral):
            self._dihedrals.append(d)
            self._residues = concatenate([self.residues, [d.residue]])
            self._atoms = concatenate([self.atoms, d.atoms])
        elif isinstance(d, Dihedrals):
            self._dihedrals.extend(d)
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
    
    @property
    def atoms(self):
        '''
        An ordered 1D list of all atoms corresponding to the dihedrals.
        '''
        if self._atoms is None:
            import numpy
            atoms = numpy.empty([len(self),4],dtype='object')
            for i, d in enumerate(self):
                atoms[i] = d.atoms
            from chimerax.core.atomic import Atoms
            self._atoms = Atoms(numpy.ravel(atoms))
        return self._atoms    
    
    @property
    def residues(self):
        '''
        Returns a list of the residues which "own" the dihedrals (one 
        per dihedral).
        '''
        if self._residues is None:
            residues = [d.residue for d in self] 
            from chimerax.core.atomic import Residues
            self._residues = Residues(residues)            
        return self._residues
        
    @property
    def coords(self):
        '''
        Returns an (n*4)*3 numpy array of the current coordinates of the
        dihedral atoms.
        '''
        return self.atoms.coords
            
    @property
    def values(self):
        '''
        Get the current torsion angles (in radians) for all dihedrals as
        a 1D numpy array. C++-accelerated function, about 1 millisecond
        per 1,500 dihedrals.
        '''
        from . import geometry
        return geometry.get_dihedrals(self.coords, len(self))
    
    @property
    def rama_cases(self):
        '''
        MolProbity Ramachandran cases for each dihedral (where applicable).
        Values are managed by the wrapping Backbone_Dihedrals class.
        '''
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
    def __init__(self, session, model = None, phi = None, psi = None, omega = None):
        if model == None and (phi == None or psi == None or omega == None):
            raise TypeError('You must provide either a model or all three of\
                            phi, psi and omega!')
        elif model and (phi or psi or omega):
            raise TypeError('Cannot provide both a model and predefined dihedrals!')
        elif model and not model.atomspec_has_atoms():
            raise TypeError('Please provide a model containing atoms!')
        
        self.session = session
        import numpy
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
            self.residues = self.residues.filter(f
            # ChimeraX provides no guarantee that residue lists will be
            # in the right order (although they usually are), so we need
            # to sort them explicitly by chain ID and residue number.
            self.residues = self.residues\
                [numpy.lexsort(\
                    (self.residues.numbers, self.residues.chains.chain_ids))]
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
            self.find_dihedrals()
        else:
            self.phi = phi
            self.psi = psi
            self.omega = omega
            # Make absolutely sure we've gotten all the residues
            self.residues = phi.residues.merge(psi.residues).merge(omega.residues)
            # ChimeraX provides no guarantee that residue lists will be
            # in the right order (although they usually are), so we need
            # to sort them explicitly by chain ID and residue number.
            self.residues = self.residues[numpy.lexsort((self.residues.numbers, self.residues.chains.chain_ids))]
            self.resnames = self.residues.names
            atoms = self.residues.atoms
            self.CAs = atoms.filter(atoms.names == 'CA')
            phi_indices = []
            psi_indices = []
            omega_indices = []
            phi_indices = self.residues.indices(phi.residues)
            psi_indices = self.residues.indices(psi.residues)
            omega_indices = self.residues.indices(omega.residues)
            # Filter out dihedrals involving residues not in the input  
            # list for that dihedral type.
            # (return -1 in self.residues.indices())
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
        '''
        phi, psi, omega = {Backbone_Dihedrals}.all_vals
        
        Returns all current phi, psi and omega values in radians as three 
        equal-length numpy arrays. Where a residue is missing a dihedral, 
        the corresponding position in the array will be nan.
        '''
        return self.phi_vals, self.psi_vals, self.omega_vals
    
    @property
    def phi_vals(self):
        '''
        Measure and return all phi values (in radians) as a numpy array 
        of length equal to the number of residues. Where a residue has 
        no phi dihedral, the array value will be nan. 
        '''
        if self._phi_vals is None:
            import numpy
            self._phi_vals = numpy.ones(len(self.residues),numpy.float32)*numpy.nan
        self._phi_vals[self._phi_indices] = self.phi.values
        return self._phi_vals
    
    @property
    def psi_vals(self):
        '''
        Measure and return all psi values (in radians) as a numpy array 
        of length equal to the number of residues. Where a residue has 
        no psi dihedral, the array value will be nan. 
        '''
        if self._psi_vals is None:
            import numpy
            self._psi_vals = numpy.ones(len(self.residues),numpy.float32)*numpy.nan
        self._psi_vals[self._psi_indices] = self.psi.values
        return self._psi_vals
    
    @property    
    def omega_vals(self):
        '''
        Measure and return all omega values (in radians) as a numpy array 
        of length equal to the number of residues. Where a residue has 
        no omega dihedral, the array value will be nan. 
        '''
        if self._omega_vals is None:
            import numpy
            self._omega_vals = numpy.ones(len(self.residues),numpy.float32)*numpy.nan
        self._omega_vals[self._omega_indices] = self.omega.values
        return self._omega_vals
    
    @property
    def rama_cases(self):
        '''
        MolProbity Ramachandran cases for each residue. N-terminal and 
        C-terminal residues are inapplicable to the Ramachandran plot,
        and are sorted into special cases. The cases are initialised 
        on creation of the object, and prolines are re-sorted into 
        CisPro and TransPro every time this property is accessed.
        Valid entries:
            'Nter', 'Cter', 'CisPro', 'TransPro', 'Glycine', 'PrePro', 
            'IleVal', 'General'
        '''
        self.update_pro_rama_cases(self.omega_vals)
        return self._rama_cases
    
    @property
    def rama_scores(self):
        '''
        Numpy array holding the current MolProbity probability scores 
        for all residues. This array is passed to the RamaValidator to
        be filled. Until RamaValidator.get_scores() has been called on
        this object, all scores will be -1. Residues which cannot be
        placed on the Ramachandran plot will retain scores of -1.
        '''
        if self._rama_scores is None:
            import numpy
            self._rama_scores = numpy.ones(len(self.residues),numpy.float32) * -1
        return self._rama_scores
    
    @property
    def rama_colors(self):
        '''
        Intended for mapping of Ramachandran scores to atoms for live
        visual feedback. Holds a (n*4) numpy array of integers, where 
        each entry corresponds to (r,g,b,alpha) on a scale of 0-255.
        The colour scale is set by the RamaValidator. 
        Until RamaValidator.get_scores(update_colors = True) is called
        on this object, all colours will remain the default grey. 
        Residues which cannot be placed on the Ramachandran plot will 
        remain grey.
        '''
        if self._rama_colors is None:
            import numpy
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
        import numpy
        phi = self.phi.by_residues(reslist)
        psi = self.psi.by_residues(reslist)
        omega = self.omega.by_residues(reslist)
        return phi, psi, omega
    
    def update_pro_rama_cases(self, omega_vals):
        '''
        Re-sort Proline residues into cis and trans based on their current
        omega dihedrals.
        '''
        import numpy
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
        Old, unused version. New code is about fifteen times faster,
        but this is retained because it's perhaps a little easier to 
        follow.
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
        import numpy
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
        import numpy
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
            
        CN_bonds = CN_atoms.inter_bonds
        
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
            omega.append(Dihedral(Atoms(omega_1d[4*i:4*i+4]), omega_res[i]))
        omega = self.omega = Dihedrals(omega)
        
        phi_1d = numpy.ravel(raw_phis)
        phi_res = res[phi_i]
        phi = []
        for i in range(len(raw_phis)):
            phi.append(Dihedral(Atoms(phi_1d[4*i:4*i+4]), phi_res[i]))
        phi = self.phi = Dihedrals(phi)
        
        psi_1d = numpy.ravel(raw_psis)
        psi_res = res[psi_i]
        psi = []
        for i in range(len(raw_psis)):
            psi.append(Dihedral(Atoms(psi_1d[4*i:4*i+4]), psi_res[i]))
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
        
        
        
        
        
        
        
       
        
    
        
        
        
