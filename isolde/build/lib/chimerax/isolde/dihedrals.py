
def time_dihedrals(sel, iterations):
    from time import time
    start_time = time()
    dlist = Backbone_Dihedrals(sel)
    end_time = time()
    print('Caching dihedrals for ' + str(len(dlist.residues)) + ' residues took ' + str(end_time - start_time) + ' seconds.')
    
    start_time = time()
    for i in range(iterations):
        for i in range(len(dlist.residues)):
            if dlist.phi[i] is not None:
                dlist.phi[i].get_value()
            if dlist.psi[i] is not None:
                dlist.psi[i].get_value()
            if dlist.omega[i] is not None:
                dlist.omega[i].get_value()
    end_time = time()    
    print('Calculating dihedrals ' + str(iterations) + ' times took ' + str(end_time - start_time) + ' seconds.')



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
        
        self.value = None
        
    def update(self):
        from .geometry import get_dihedral
        self.coords = self.atoms.coords
        self.value = get_dihedral(*self.coords)
    
    def get_value(self):
        self.update()
        return self.value

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
    Given an atomic model, finds all amino acid residues and creates
    Dihedrals objects containing their phi, psi and omega dihedrals.
    '''
    def __init__(self, model):
        if not model.atomspec_has_atoms():
            raise TypeError('Please provide a model containing atoms!')
        self.residues = model.residues
        # Filter to get only amino acid residues
        f = []
        for r in self.residues:
            f.append(r.PT_AMINO)
        if not len(f):
            return
        import numpy
        f = numpy.array(f, dtype = bool)
        self.residues = self.residues.filter(f)
        
        self.residues = self.residues[numpy.lexsort((self.residues.numbers, self.residues.chains.chain_ids))]
        self.resnames = self.residues.names
        self.atoms = self.residues.atoms
        # We want to colour C-alphas according to their status, so we'll
        # hold an array of them here.
        self.CAs = self.atoms.filter(self.atoms.names == 'CA')
        
        nr = self.num_residues = len(self.residues)
        
                
        self.find_dihedrals()
    
    def get_phi_vals(self, update = True):
        if update:
            for i, phi in enumerate(self.phi):
                if phi is not None:
                    self.phi_vals[i] = phi.get_value()
                else:
                    self.phi_vals[i] = None
        return self.phi_vals
    
    def get_psi_vals(self, update = True):
        if update:
            for i, psi in enumerate(self.psi):
                if psi is not None:
                    self.psi_vals[i] = psi.get_value()
                else:
                    self.phi_vals[i] = None
        return self.psi_vals
        
    def get_omega_vals(self, update = True):
        if update:
            for i, omega in enumerate(self.omega):
                if omega is not None:
                    self.omega_vals[i] = omega.get_value()
                else:
                    self.omega_vals[i] = None
        return self.omega_vals

    def lookup_by_residue(self, res):
        ''' 
        Return the phi, psi and omega dihedrals for a given residue
        '''
        i = self.residues.index(res)
        if i == -1:
            return None
        return (self.phi[i], self.psi[i], self.omega[i])
         
    def find_dihedrals(self):
        # Empty Dihedrals objects to fill with phi, psi and omega
        self.phi = []
        self.psi = []
        self.omega = []
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
            if not r.PT_AMINO:
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
            if psi_atoms is not None:
                self.psi.append(Dihedral(concatenate(psi_atoms, Atoms), residue = r))
            if omega_atoms is not None:
                self.omega.append(Dihedral(concatenate(omega_atoms, Atoms), residue = r))
        
        # Convert to Dihedrals objects
        self.phi = Dihedrals(self.phi)
        self.psi = Dihedrals(self.psi)
        self.omega = Dihedrals(self.omega)
                    
                    
                        
                
            
            
        
                
            
        
        
    
        
        
        
