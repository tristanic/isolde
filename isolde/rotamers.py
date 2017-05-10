import os
import numpy
import pickle
from math import degrees, radians
from chimerax.core.geometry import Place, rotation
from chimerax.core.atomic import Residue

from .dihedrals import Dihedral, Dihedrals

package_directory = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(package_directory, 'molprobity_data')

_rot_file_prefix = os.path.join(DATA_DIR, 'penultimate_rotamer_library')

def load_rotamers(file_prefix):
    '''
    Tries to load the rotamer data from a pickle file. If that fails, reads
    from text and creates the pickle file for next time.
    Text file format:
    
    Name Overall alpha beta other χ1 χ2 χ3 χ4
    Arginine ARG 
    χ1 N CA CB CG CD CG CZ HB1 HB2 HD1 HD2 HE HG1 HG2 HH11 HH12 HH21 HH22 NE NH1 NH2
    χ2 CA CB CG CD CD CZ HD1 HD2 HE HG1 HG2 HH11 HH12 HH21 HH22 NE NH1 NH2
    χ3 CB CG CD NE  CZ HD1 HD2 HE HH11 HH12 HH21 HH22 NE NH1 NH2
    χ4 CG CD NE CZ  CZ HE HH11 HH12 HH21 HH22 NH1 NH2
    ptp85° 0.002 0.00 0.01 0.002 62 180 65 85
    ...
    Lysine LYS 
    χ1 N CA CB CG CD CE CG HB1 HB2 HD1 HD2 HE1 HE2 HG1 HG2 HZ1 HZ2 HZ3 NZ
    χ2 CA CB CG CD CD CE HD1 HD2 HE1 HE2 HG1 HG2 HZ1 HZ2 HZ3 NZ
    χ3 CB CG CD CE CE HD1 HD2 HE1 HE2 HZ1 HZ2 HZ3 NZ
    χ4 CG CD CE NZ HE1 HE2 HZ1 HZ2 HZ3 NZ
    ptpt 0.01 0.00 0.02 0.002 62 180 68 180
    ...
    Name Overall alpha beta other χ1 χ2 χ3
    Methionine MET
    ...
    Special
    Name Overall  χ2 χ3 χ2'
    Disulfide CYS
    
    
    
    Lines beginning with 'Name' are used to set the number of rotamers for
    the following residues. Immediately after the line giving the name of
    each residue, there is a set of lines defining its Chi dihedrals:
        Name atom1 atom2 atom3 atom4 {list of atoms controlled by dihedral}
    
    At the end of a file there is a Special section for rotamers that cannot
    be handled by standard code (e.g. those that interact with other 
    residues). At the moment the only special case is disulfide-bonded
    cysteines, but this is likely to grow to include e.g. glycosylated 
    residues.
    '''
    infile = None
    try:
        infile = open(file_prefix+'.pickle', 'r+b')
        rotamers = pickle.load(infile)
        infile.close()
        infile = None
    except:
        if infile is not None:
            infile.close()
        infile = open(file_prefix+'.data', 'r')
        lines = [line.rstrip('\n ').split(' ') for line in infile]
        infile.close()
        rotamers = {}
        special = False
        i = 0
        while i < len(lines):
            l = lines[i]
            if l[0] == 'Name':
                # Number of torsions has changed.
                num_chi_vals = len(l) - 5
                chi_names = l[6:]
                i += 1
                continue
            if len(l) == 2:
                # We're starting on a new residue
                full_name = l[0]
                key = l[1]
                chi_dihedrals = []
                for j in range(num_chi_vals):
                    i += 1
                    l = lines[i]
                    name = l[0]
                    atoms = l[1:5]
                    moving_atoms = l[5:]
                    chi_dihedrals.append(Chi_dihedral(name, atoms, moving_atoms))
                
                rotamers[key] = Residue_rotamers(full_name, chi_dihedrals)
                i += 1
                continue
            if l[0] == 'Special':
                special = True
                i += 1
                continue
            if not special:
                name = l[0]
                probs = [float(p) for p in l[1:5]]
                angles = [float(a) for a in l[5:]]
                r = Rotamer_info(name, *probs, angles)
                rotamers[key].add_rotamer_info(r)
                i += 1
            else:
                # Special cases are in the too-hard basket for now.
                i += 1
        out = open(file_prefix+'.pickle', 'w+b')
        pickle.dump(rotamers, out)
        out.close()
    return rotamers 

            
class Residue_rotamers:
    '''
    Holds all the standard rotamers with their conformation-dependent
    probabilities for one residue, with name look-up and sorting function.
    Also holds the necessary information to describe each dihedral:
        - the atoms making up the given dihedral
        - the atoms which must move on changing the torsion angle
    '''
    def __init__(self, name, chi_dihedrals):
        self.name = name
        self.dihedrals = chi_dihedrals
        self.num_torsions = len(chi_dihedrals)
        self._rotamers = []
        self._by_overall = None
        self._by_alpha = None
        self._by_beta = None
        self._by_coil = None
    
    def add_rotamer_info(self, rotamer_info):
        '''
        Add the information describing one rotamer
        '''
        self._rotamers.append(rotamer_info)
    
    @property
    def rotamers(self):
        return self._rotamers

    @property
    def by_overall_P(self):
        '''
        Return the list of possible rotamers sorted in order of overall 
        probability (highest probability first)
        '''
        if self._by_overall is None:
            self._by_overall = sorted(
                self._rotamers, key = lambda r: r.overall_P, reverse = True)
        return self._by_overall
    
    @property
    def by_alpha_P(self):
        '''
        Return the list of possible rotamers sorted in order of 
        probability of appearing in an alpha helix (highest probability 
        first)
        '''
        if self._by_alpha is None:
            self._by_alpha = sorted(
                self._rotamers, key = lambda r: r.alpha_P, reverse = True)
        return self._by_alpha

    @property
    def by_beta_P(self):
        '''
        Return the list of possible rotamers sorted in order of 
        probability of appearing in a beta strand (highest probability 
        first)
        '''
        if self._by_beta is None:
            self._by_beta = sorted(
                self._rotamers, key = lambda r: r.beta_P, reverse = True)
        return self._by_beta
    
    @property
    def by_coil_P(self):
        '''
        Return the list of possible rotamers sorted in order of 
        probability of appearing in a random coil (highest probability 
        first)
        '''
        if self._by_coil is None:
            self._by_coil = sorted(
                self._rotamers, key = lambda r: r.coil_P, reverse = True)
        return self._by_coil
    
    def by_ss_type(self, ss_type):
        if ss_type == Residue.SS_COIL:
            return self.by_coil_P
        if ss_type == Residue.SS_HELIX:
            return self.by_alpha_P
        return self.by_beta_P

class Chi_dihedral:
    '''
    Describes one dihedral - its atoms, and the atoms that need to move
    when its angle changes.
    '''
    def __init__(self, name, dihedral_atom_names, moving_atom_names):
        self.name = name
        self.atom_names = dihedral_atom_names
        self.moving_names = moving_atom_names
    
class Rotamer_info:
    '''
    Holds the information necessary to describe one rotamer.
    '''
    def __init__(self, name, 
                    overall_P, alpha_P, beta_P, coil_P, chi_angles):
        '''
        Define a rotamer.
        Args:
            name:
                The standard rotamer name (e.g. ptp85°)
            overall_P:
                The probability of finding this rotamer regardless of 
                backbone conformation
            alpha_P:
                The probability of finding this rotamer in an alpha helix
            beta_P:
                The probability of finding this rotamer in a beta strand
            coil_P:
                The probability of finding this rotamer in a random coil
            chi_P:
                The torsion angles for this rotamer in degrees
        '''
        self.name = name
        self.overall_P = overall_P
        self.alpha_P = alpha_P
        self.beta_P = beta_P
        self.coil_P = coil_P
        self._angles = [radians(a) for a in chi_angles]
    
    @property
    def angles(self):
        return self._angles
    
    @property
    def angles_deg(self):
        return [degrees(a) for a in self._angles]
        
class Rotamer:
    '''
    Describes an actual rotamer found in a protein, and provides the 
    methods necessary to query/change it.
    '''
    def __init__(self, residue):
        '''
        Just create the rotamer object. No querying yet.
        Args:
            residue:
                A ChimeraX Residue object
        '''
        self.residue = residue
        self.resname = residue.name
        try:
            self.rotamers = _rotamer_info[self.resname]
        except:
            raise KeyError('No rotamer information available for {}!'\
                            .format(self.resname))
        ratoms = residue.atoms
        anames = ratoms.names
        self._ghost_drawing = None
        self.dihedrals = []
        self._atoms_to_move = []
        for c in self.rotamers.dihedrals:
            indices = []
            for a in c.atom_names:
                indices.append(anames.tolist().index(a))
            self.dihedrals.append(Dihedral(ratoms[numpy.array(indices)], residue))
            
            self._atoms_to_move.append(ratoms[numpy.in1d(anames, c.moving_names)])
        
        self.dihedrals = Dihedrals(self.dihedrals)
        
    
    @property
    def chi_angles(self):
        return self.dihedrals.values
    
    @property
    def chi_angles_deg(self):
        return numpy.degrees(self.dihedrals.values)
    
    @property
    def available_rotamers(self):
        return _rotamer_info[self.residue.name].by_ss_type(self.residue.ss_type)
    
    def change_rotamer(self, target, draw_ghost = True):
        '''
        Move the atoms in this residue to match the rotamer defined by
        the target.
        Args:
            target:
                A Rotamer_info object matching the current residue
        '''
        target_angles = target.angles_deg
        current_angles = self.chi_angles_deg
        difference = target_angles - current_angles
        for angle, chi, moveatoms in zip(difference, self.dihedrals, self._atoms_to_move):
            axis = chi.coords[2] - chi.coords[1]
            center = chi.coords[2]
            tf = rotation(axis, angle, center)
            moveatoms.coords = tf.moved(moveatoms.coords)

        
        
            
    
                
_rotamer_info = load_rotamers(_rot_file_prefix)
        
