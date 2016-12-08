
class Dihedral():
    '''
    Holds information about the dihedral formed by an array of four atoms.
    The atoms must be an Atom array, ordered so that the torsion axis is
    between the middle two atoms. This axis needn't necessarily correspond to
    a physical bond.
    '''
    
    def __init__(self, atoms, residue = None):
        self.atoms = atoms
        self.coords = atoms.coords
        self.names = atoms.names
        self.residue = residue # the Residue object that this dihedral belongs to
        self.residues = atoms.residues
        self.resnames = self.residues.names
        self.resnums = self.residues.numbers
        # Index of the matching CustomTorsionForce in the simulation so
        # we can restrain this dihedral as needed
        self.sim_index = -1
        
    
    @property
    def value(self):
        from .geometry import get_dihedral
        self.coords = self.atoms.coords
        return get_dihedral(*self.coords)

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
        if dlist is not None:
            self._dihedrals = dlist
            residues = [d.residue for d in self._dihedrals] 
            from chimerax.core.atomic import Residues
            self._residues = Residues(residues)
            atoms = [a for d in self.dihedrals for a in d.atoms]
            from chimerax.core.atomic import Atoms
            self._atoms = Atoms(atoms)
    
    def __len__(self):
        return len(self.dihedrals)
    
    def __bool__(self):
        return len(self) > 0
    
    def __iter__(self):
        return iter(self.dihedrals)
    
    def __getitem__(self, i):
        import numpy
        if isinstance(i,(int,numpy.int32)):
            return self.dihedrals[i]
        elif isinstance(i,(slice)):
            return Dihedrals(self.dihedrals[i])
        elif isinstance(i, numpy.ndarray):
            return Dihedrals([self.dihedrals[j] for j in i])
        else:
            raise IndexError('Only integer indices allowed for %s, got %s'
                % (self.__class__.__name__, str(type(i))))
    
    def by_residue(self, residue):
        from chimerax.core.atomic import Residues
        r = Residues([residue])
        i = r.indices(self.residues)[0]
        if i == -1:
            return None
        return self[i]
    
    def by_residues(self, residues):
        indices = residues.indices(self.residues)
        indices = indices[indices != -1]
        if len(indices):
            return self[indices]
        return None
   
    def append(self, d):
        from chimerax.core.atomic import concatenate, Residues
        if isinstance(d, Dihedral):
            self.dihedrals.append(d)
            self.residues = concatenate([self.residues, [d.residue]])
            self.atoms = concatenate([self.atoms, d.atoms])
        elif isinstance(d, Dihedrals):
            self.dihedrals.extend(d)
            self.residues = concatenate([self.residues, d.residues])
            self.atoms = concatenate([self.atoms, d.atoms])
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
        return self._atoms
    
    @property
    def residues(self):
        return self._residues
    
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

class Backbone_Dihedrals():
    '''
    Takes either an atomic model or three Dihedrals objects with pre-defined
    phi, psi and omega dihedrals. If provided with a model, it will find and
    store all phi, psi and omega dihedrals as Dihedrals objects.
    '''
    def __init__(self, model = None, phi = None, psi = None, omega = None):
        if model == None and (phi == None or psi == None or omega == None):
            raise TypeError('You must provide either a model or all three of\
                            phi, psi and omega!')
        elif model and (phi or psi or omega):
            raise TypeError('Cannot provide both a model and predefined dihedrals!')
        elif model and not model.atomspec_has_atoms():
            raise TypeError('Please provide a model containing atoms!')
            
        import numpy
        # It's most convenient to determine and store the Ramachandran case
        # for each residue here, otherwise things start to get messy when
        # working with subsets.
        self.rama_case = []
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
            self.find_dihedrals()
        else:
            self.phi = phi
            self.psi = psi
            self.omega = omega
            
            self.residues = phi.residues.merge(psi.residues).merge(omega.residues)
            self.residues = self.residues[numpy.lexsort((self.residues.numbers, self.residues.chains.chain_ids))]
            self.resnames = self.residues.names
            atoms = self.residues.atoms
            self.CAs = atoms.filter(atoms.names == 'CA')
            phi_indices = []
            psi_indices = []
            omega_indices = []
            phi_indices = phi.residues.indices(self.residues)
            psi_indices = psi.residues.indices(self.residues)
            omega_indices = omega.residues.indices(self.residues)
            
            #for i, r in enumerate(self.residues):
                #if phi.by_residue(r) != None:
                    #phi_indices.append(i)
                #if psi.by_residue(r) != None:
                    #psi_indices.append(i)
                #if omega.by_residue(r) != None:
                    #omega_indices.append(i)
            self._phi_indices = numpy.array(phi_indices[phi_indices != -1],numpy.int32)
            self._psi_indices = numpy.array(psi_indices[psi_indices != -1],numpy.int32)
            self._omega_indices = numpy.array(omega_indices[omega_indices != -1],numpy.int32)
                
            
    @property
    def all_vals(self):
        return self.phi_vals, self.psi_vals, self.omega_vals
    
    @property
    def phi_vals(self):
        import numpy
        vals = numpy.ones(len(self.residues),numpy.float32)*numpy.nan
        vals[self._phi_indices] = self.phi.values
        return vals
    
    @property
    def psi_vals(self):
        import numpy
        vals = numpy.ones(len(self.residues),numpy.float32)*numpy.nan
        vals[self._psi_indices] = self.psi.values
        return vals
    
    @property    
    def omega_vals(self):
        import numpy
        vals = numpy.ones(len(self.residues),numpy.float32)*numpy.nan
        vals[self._omega_indices] = self.omega.values
        return vals

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
         
    def find_dihedrals(self):
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
        if len(self.phi) or len(self.psi) or len(self.omega):
            import warnings
            warnings.warn('Backbone dihedrals have already been defined. \
                           If you want to update them, create a new \
                           Backbone_Dihedrals object.')
            return            
        import numpy
        from chimerax.core.atomic import Residue
            
        res = self.residues.filter(self.residues.polymer_types == Residue.PT_AMINO)
        atoms = res.atoms
        keyatoms = numpy.empty([len(res),3],dtype='object')
        
        for i, name in enumerate(['N','CA', 'C']):
            keyatoms[:,i] = atoms.filter(atoms.names == name)
            
        
        
                
            
        
        
    
        
        
        
