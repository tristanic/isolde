
def time_dihedrals(sel, iterations):
    from time import time
    start_time = time()
    dlist = Backbone_Dihedrals(sel)
    end_time = time()
    print('Caching dihedrals for ' + str(len(dlist.residues)) + ' residues took ' + str(end_time - start_time) + ' seconds.')
    
    start_time = time()
    for i in range(iterations):
        for res, phi in dlist.phi.items():
            phi.get_value()
        for res, psi in dlist.psi.items():
            psi.get_value()
        for res, omega in dlist.omega.items():
            omega.get_value()
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
        self.atoms = atoms
        self.residues = atoms.unique_residues
        
        # Dicts holding the dihedrals for each residue, with the residue
        # as the key
        self.phi = {}
        self.psi = {}
        self.omega = {}
        self.find_dihedrals()
        
    
    def find_dihedrals(self):
        bond_to_last = False
        bond_to_next = False
        last_residue = None
        last_atoms = None
        last_names = None
        import numpy
        from copy import copy
        for r in self.residues:
            if not r.PT_AMINO:
                bond_to_last = False
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
            if not bond_to_last:
                # Double-check to see if this residue is bonded to a previous
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
                
            if bond_to_next:
                next_atoms = next_residue.atoms
                next_names = next_atoms.names
                CA = next_atoms.filter(numpy.in1d(next_names, 'CA'))
                psi_atoms.append(CA)
            else:
                # C-terminal residues have no psi dihedral
                psi_atoms = None
            from chimerax.core.atomic import concatenate, Atoms
            if phi_atoms is not None:
                self.phi[r] = Dihedral(concatenate(phi_atoms, Atoms))
            if psi_atoms is not None:
                self.psi[r] = Dihedral(concatenate(psi_atoms, Atoms))
            if omega_atoms is not None:
                self.omega[r] = Dihedral(concatenate(omega_atoms, Atoms))
            
            bond_to_last = bond_to_next
            last_atoms = a
            last_names = names
            #a = next_atoms
            #names = next_names
            
                
                    
                    
                        
                
            
            
        
                
            
        
        
    
        
        
        
