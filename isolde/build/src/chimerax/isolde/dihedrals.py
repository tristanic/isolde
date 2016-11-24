
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

# Holds information about the dihedral formed by an array of four atoms.
# The atoms must be an Atom array, ordered so that the torsion axis is
# between the middle two atoms. This axis needn't necessarily correspond to
# a physical bond.

class Dihedral():
    
    def __init__(self, atoms):
        self.atoms = atoms
        self.coords = atoms.coords
        self.names = atoms.names
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


# Given an atom selection, finds contiguous protein stretches and stores
# arrays containing the phi, psi and omega dihedrals    
class Backbone_Dihedrals():
    def __init__(self, atoms):
        self.residues = atoms.unique_residues
        # Filter to get only amino acid residues
        f = []
        for r in self.residues:
            f.append(r.PT_AMINO)
        import numpy
        f = numpy.array(f, dtype = bool)
        self.residues = self.residues.filter(f)
        
        self.residues = self.residues[numpy.lexsort((self.residues.numbers, self.residues.chains.chain_ids))]
        self.resnames = self.residues.names
        self.atoms = self.residues.atoms
        # We want to colour C-alphas according to their status, so we'll
        # hold an array of them here.
        self.CAs = atoms.filter(atoms.names == 'CA')
        
        nr = self.num_residues = len(self.residues)
        
        # Lists holding the dihedrals for each residue, in the same order
        # as the lists of residues
        self.phi = [None] * nr
        self.phi_vals = [None] * nr
        self.psi = [None] * nr
        self.psi_vals = [None] * nr
        self.omega = [None] * nr
        self.omega_vals = [None] * nr
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

    
    
    def find_dihedrals(self):
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
                phi_atoms.append(a.filter(numpy.in1d(names, [name])))
            # These three atoms are common to both phi and psi, so we can
            # just copy the array
            psi_atoms = copy(phi_atoms)
            omega_atoms = phi_atoms[0:2]
            # Atoms from the previous residue
            # Check to see if this residue is bonded to a previous
            # one.
            N = phi_atoms[0][0]
            N_bonded_atom_list = N.bonds.atoms[0].merge(N.bonds.atoms[1])
            if -1 in N_bonded_atom_list.indices(self.atoms):
                bond_to_last = False
            else:
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
            if -1 in C_bonded_atom_list.indices(self.atoms):
                bond_to_next = False
            else:
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
                self.phi[i] = Dihedral(concatenate(phi_atoms, Atoms))
            if psi_atoms is not None:
                self.psi[i] = Dihedral(concatenate(psi_atoms, Atoms))
            if omega_atoms is not None:
                self.omega[i] = Dihedral(concatenate(omega_atoms, Atoms))
            
                
                    
                    
                        
                
            
            
        
                
            
        
        
    
        
        
        
