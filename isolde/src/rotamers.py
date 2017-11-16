# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


import os
import numpy
import pickle
import contextlib
import itertools
from math import degrees, radians
from chimerax.core.geometry import Place, rotation, matrix
from chimerax.core.atomic import Residue, Atom

from simtk.unit import Quantity

from .dihedrals import Dihedral, Dihedrals
from .constants import defaults
package_directory = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(package_directory, 'molprobity_data')


'''
If the abundance of a rotamer (for a given secondary structure) is less
than RELATIVE_ABUNDANCE_CUTOFF, it will be removed from the list of
valid rotamers for that secondary structure type.
'''
RELATIVE_ABUNDANCE_CUTOFF = 0.0001

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
        key = None
        while i < len(lines):
            l = lines[i]
            if l[0] == 'Special':
                # Special cases are in the too-hard basket for now
                break
            if l[0] == 'Name':
                # Number of torsions has changed.
                num_chi_vals = len(l) - 5
                chi_names = l[6:]
                i += 1
                continue
            if len(l) == 2:
                # We're starting on a new residue
                # Normalise the probabilities for the last one
                if key is not None:
                    rotamers[key]._normalize_probabilities()
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
            if not special:
                name = l[0]
                probs = [float(p) for p in l[1:5]]
                angles = [float(a) for a in l[5:]]
                r = Rotamer_info(name, *probs, angles)
                rotamers[key].add_rotamer_info(r)
                i += 1
        out = open(file_prefix+'.pickle', 'w+b')
        pickle.dump(rotamers, out)
        out.close()
    return rotamers

def all_rotamers_in_selection(session, atoms):
    '''
    Returns a {Residue: Rotamer} dict covering every rotameric residue
    in the input Atoms object.
    '''
    residues = atoms.unique_residues
    ret = {}
    for r in residues:
        try:
            ret[r] = Rotamer(session, r)
        except KeyError:
            continue
    return ret


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

    def _normalize_probabilities(self):
        best_overall_P = self.by_overall_P[0].overall_P
        best_alpha_P = self.by_alpha_P[0].alpha_P
        best_beta_P = self.by_beta_P[0].beta_P
        best_coil_P = self.by_coil_P[0].coil_P
        for r in self._rotamers:
            r.overall_P /= best_overall_P
            r.alpha_P /= best_alpha_P
            r.beta_P /= best_beta_P
            r.coil_P /= best_coil_P


    @property
    def zeros(self):
        '''
        Return a numpy array of zeros, one for each chi dihedral. Used
        to initialise the target array.
        '''
        return numpy.zeros(self.num_torsions, float)

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
            temp = [r for r in self._rotamers if r.overall_P > RELATIVE_ABUNDANCE_CUTOFF]
            self._by_overall = sorted(
                temp, key = lambda r: r.overall_P, reverse = True)

        return self._by_overall

    @property
    def by_alpha_P(self):
        '''
        Return the list of possible rotamers sorted in order of
        probability of appearing in an alpha helix (highest probability
        first)
        '''
        if self._by_alpha is None:
            temp = [r for r in self._rotamers if r.alpha_P > RELATIVE_ABUNDANCE_CUTOFF]
            self._by_alpha = sorted(
                temp, key = lambda r: r.alpha_P, reverse = True)
        return self._by_alpha

    @property
    def by_beta_P(self):
        '''
        Return the list of possible rotamers sorted in order of
        probability of appearing in a beta strand (highest probability
        first)
        '''
        if self._by_beta is None:
            temp = [r for r in self._rotamers if r.beta_P > RELATIVE_ABUNDANCE_CUTOFF]
            self._by_beta = sorted(
                temp, key = lambda r: r.beta_P, reverse = True)
        return self._by_beta

    @property
    def by_coil_P(self):
        '''
        Return the list of possible rotamers sorted in order of
        probability of appearing in a random coil (highest probability
        first)
        '''
        if self._by_coil is None:
            temp = [r for r in self._rotamers if r.coil_P > RELATIVE_ABUNDANCE_CUTOFF]
            self._by_coil = sorted(
                temp, key = lambda r: r.coil_P, reverse = True)
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
        self._angles = numpy.array([radians(a) for a in chi_angles])

    @property
    def angles(self):
        return self._angles

    @property
    def angles_deg(self):
        return [degrees(a) for a in self._angles]

    def relative_abundance(self, residue):
        ss_type = residue.ss_type
        if ss_type == Residue.SS_COIL:
            return self.coil_P
        if ss_type == Residue.SS_HELIX:
            return self.alpha_P
        return self.beta_P

def _update_factory():
    return False

class Rotamers:
    '''
    Holds a set of Rotamer objects, with look-ups to get all rotamers of a
    specific type, or quickly get all dihedrals for the given rotamers
    '''
    def __init__(self, session, residues = None, master_rotamers = None):
        self.residues = residues
        self.session = session
        from collections import defaultdict
        self._needs_update = defaultdict(_update_factory)


        self._rotamers = {}

        # Entry 0: a list of rotamers of this type, in the order they appear
        # in the structure
        # Entry 1: the rotameric dihedrals for these residues, as a single
        # Dihedrals object
        # Entry 2: the number of rotameric dihedrals for this residue type
        # Entry 3: Is the terminal dihedral symmetric?
        self._by_type = {
            'ARG':          [list(), None, 4, False],
            'ASN':          [list(), None, 2, False],
            'ASP':          [list(), None, 2, True],
            'CYS':          [list(), None, 1, False],
            'GLN':          [list(), None, 3, False],
            'GLU':          [list(), None, 3, True],
            'HIS':          [list(), None, 2, False],
            'ILE':          [list(), None, 2, False],
            'LEU':          [list(), None, 2, False],
            'LYS':          [list(), None, 4, False],
            'MET':          [list(), None, 3, False],
            'PHE':          [list(), None, 2, True],
            'PRO':          [list(), None, 1, False],
            'SER':          [list(), None, 1, False],
            'THR':          [list(), None, 1, False],
            'TRP':          [list(), None, 2, False],
            'TYR':          [list(), None, 2, True],
            'VAL':          [list(), None, 1, False]
        }
        for residue in residues:
            try:
                if master_rotamers is None:
                    self[residue] = Rotamer(session, residue)
                else:
                    self[residue] = master_rotamers[residue]
            except KeyError:
                continue
    
    def dihedrals(self, resname):
        '''
        Get a Dihedrals object containing all the dihedrals for residues of this
        type. Dihedrals will be ordered as:
        [Chi1_1, ... ChiN_1, ... Chi1_M, ChiN_M]
        where N is the number of chi dihedrals for this rotamer type, and M is
        the number of residues of this type.
        '''
        t = self._by_type[resname]
        if t[1] is None or self._needs_update[resname]:
            rlist = t[0]
            if len(rlist):
                t[1] = Dihedrals(list(
                    itertools.chain.from_iterable(
                        [r.dihedrals for r in rlist])))
            else:
                t[1] = Dihedrals()
            self._needs_update[resname] = False
        return t[1]

    def num_chi(self, resname):
        return self._by_type[resname][2]

    def is_symmetric(self, resname):
        return self._by_type[resname][3]

    def __setitem__(self, residue, rotamer):
        if type(residue) != Residue:
            raise TypeError('Key must be a ChimeraX Residue object!')
        if type(rotamer) != Rotamer:
            raise TypeError('Value must be a Rotamer object!')
        self._rotamers[residue] = rotamer
        resname = residue.name
        self._by_type[resname][0].append(rotamer)
        self._needs_update[resname] = True


    def __getitem__(self, residue):
        return self._rotamers[residue]

class Rotamer:
    '''
    Describes an actual rotamer found in a protein, and provides the
    methods necessary to query/change it.
    '''

    def __init__(self, session, residue,
                spring_constant = defaults.ROTAMER_SPRING_CONSTANT ):
        '''
        Just create the rotamer object. No querying yet.
        Args:
            residue:
                A ChimeraX Residue object
        '''
        self.session = session
        self._preview_model = None
        self.residue = residue
        self.resname = residue.name
        try:
            self.rotamers = _rotamer_info[self.resname]
        except:
            raise KeyError('No rotamer information available for {}!'\
                            .format(self.resname))
        self.model = residue.atoms.unique_structures[0]
        self._counter = 0
        self._frames_per_rotamer = 50
        self._cycle_handler = None
        self._current_rotamer_index = -1
        self._current_rotamer = None

        # Is this rotamer currently being restrained?
        self._restrained = False
        # Current target conformation (if self.restrained == True)
        self._target = _rotamer_info[self.resname].zeros

        try:
            self._atoms_to_move, self.dihedrals = self.find_atoms_to_move(self.residue)
        except:
            self._atoms_to_move = None
            self.dihedrals = None
        self._default_spring_constant = spring_constant

    @property
    def spring_constants(self):
        return self.dihedrals.spring_constants

    @spring_constants.setter
    def spring_constants(self, k):
        if isinstance(k, Quantity):
            k = k.value_in_unit(defaults.OPENMM_RADIAL_SPRING_UNIT)
        #self._spring_constant = k
        self.dihedrals.spring_constants = k


    @property
    def restrained(self):
        return self._restrained

    @restrained.setter
    def restrained(self, flag):
        curr = self._restrained
        if flag == curr:
            return
        if flag:
            self.dihedrals.spring_constants = self._default_spring_constant
        else:
            self.dihedrals.spring_constants = 0
            self._current_rotamer = None
        self._restrained = flag

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.cleanup()

    def cleanup(self):
        if self._cycle_handler is not None:
            self.session.triggers.remove_handler(self._cycle1_handler)
        self.remove_preview()
        if not self.restrained:
            self._current_rotamer = None

    @property
    def current_target_rotamer(self):
        return self._current_rotamer

    @property
    def target(self):
        return self._target

    @target.setter
    def target(self, target):
        self._target = target
        self.dihedrals.targets = target

    @property
    def chi_angles(self):
        return self.dihedrals.values

    @property
    def chi_angles_deg(self):
        return numpy.degrees(self.dihedrals.values)

    @property
    def available_rotamers(self):
        return _rotamer_info[self.resname].by_ss_type(self.residue.ss_type)

    def find_atoms_to_move(self, residue):
        ratoms = residue.atoms
        anames = ratoms.names
        anames_list = anames.tolist()
        dihedrals = []
        atoms_to_move = []
        for c in self.rotamers.dihedrals:
            indices = []
            for a in c.atom_names:
                try:
                    indices.append(anames_list.index(a))
                except ValueError:
                    warn_string = \
                    'Residue {}-{} on chain {} is missing atom {}! You \
                    can still analyse the structure, but simulation will\
                    be impossible until this is fixed.'.format(
                        residue.name, residue.number, residue.chain_id, a)
                    import warnings
                    warnings.warn(warn_string)
                    return
            dihedrals.append(Dihedral(ratoms[numpy.array(indices)], residue))

            atoms_to_move.append(ratoms[numpy.in1d(anames, c.moving_names)])

        return atoms_to_move, Dihedrals(dihedrals)

    def next_rotamer(self, preview = False):
        rot = self._iter_rotamer(1)
        self.change_rotamer(rot, preview = preview)
        return rot

    def previous_rotamer(self, preview = False):
        rot = self._iter_rotamer(-1)
        self.change_rotamer(rot, preview = preview)
        return rot


    def _iter_rotamer(self, incr):
        self._current_rotamer_index += incr
        self._current_rotamer_index %= len(self.available_rotamers)
        rot = self.available_rotamers[self._current_rotamer_index]
        return rot


    def change_rotamer(self, target, preview = False):
        '''
        Move the atoms in this residue to match the rotamer defined by
        the target.
        Args:
            target:
                A Rotamer_info object matching the current residue
        '''
        target_angles = target.angles_deg
        self._current_rotamer = target

        if preview:
            if self._preview_model is None or self._preview_model.deleted:
                from chimerax.core.commands.split import molecule_from_atoms
                p = self._preview_model = molecule_from_atoms(self.model, self.residue.atoms)
                p.bonds.radii = 0.1
                p.atoms.draw_modes = Atom.STICK_STYLE
                self.model.add([p])
                pres = p.residues[0]
                self._preview_atoms_to_move, self._preview_dihedrals = \
                    self.find_atoms_to_move(pres)

            atoms_to_move = self._preview_atoms_to_move
            dihedrals = self._preview_dihedrals
            current_angles = numpy.degrees(dihedrals.values)
        else:
            atoms_to_move = self._atoms_to_move
            dihedrals = self.dihedrals
            current_angles = self.chi_angles_deg

        difference = target_angles - current_angles

        for angle, chi, moveatoms in zip(difference, dihedrals, atoms_to_move):
            axis = chi.coords[2] - chi.coords[1]
            # center = chi.coords[2]
            center = matrix.project_to_axis(chi.coords[3], axis, chi.coords[1])
            tf = rotation(axis, angle, center)
            moveatoms.coords = tf.moved(moveatoms.coords)

    def cycle(self, frames_per_rotamer = 50):
        if self._cycle_handler is not None:
            return
        self._counter = 0
        self._frames_per_rotamer = frames_per_rotamer
        self._current_rotamer_index = 0
        self._current_rotamer = self.available_rotamers[0]
        self.change_rotamer(self._current_rotamer, preview = True)
        self._cycle_handler = self.session.triggers.add_handler('new frame', self._cycle_iter)

    def _cycle_iter(self, *_):
        self._counter += 1
        self._counter %= self._frames_per_rotamer
        if self._counter == 0:
            next_index = self._current_rotamer_index + 1
            if next_index < len(self.available_rotamers):
                self._current_rotamer = self.available_rotamers[next_index]
                self._current_rotamer_index = next_index
                self.change_rotamer(self._current_rotamer, preview = True)
            else:
                self.session.triggers.remove_handler(self._cycle_handler)
                self._cycle_handler = None

    def remove_preview(self):
        if self._preview_model is not None:
            self.session.models.remove([self._preview_model])
            self._preview_model = None

    def commit_current_preview(self):
        if self._preview_model is None:
            raise RuntimeError('No preview exists!')
        self.residue.atoms.coords = self._preview_model.atoms.coords
        self.remove_preview()

    def current_dihedral_angles(self):
        return self._current_rotamer.angles

_rotamer_info = load_rotamers(_rot_file_prefix)
