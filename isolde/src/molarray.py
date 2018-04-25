# @Author: Tristan Croll
# @Date:   16-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 25-Apr-2018
# @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
# @Copyright: Copyright 2017-2018 Tristan Croll



import numpy
from numpy import uint8, int32, uint32, float64, float32, uintp, byte, bool as npy_bool, integer, empty, array
from chimerax.core.atomic.molc import string, cptr, pyobject, set_cvec_pointer, pointer, size_t
from chimerax.core.atomic.molarray import Collection
from . import molobject
from .molobject import c_function, c_array_function, cvec_property
#from .molobject import object_map
from .molobject import Proper_Dihedral, Rotamer, Rama, Position_Restraint, \
        Tuggable_Atom, MDFF_Atom, Distance_Restraint, Proper_Dihedral_Restraint, \
        Rotamer_Restraint
import ctypes

from chimerax.core.atomic import Atom, Atoms, Residue, Residues

from chimerax.core.atomic.molarray import _atoms, _atoms_or_nones, \
        _bonds, _non_null_atoms, _pseudobond_groups, _pseudobonds, \
        _elements, _residues, _non_null_residues, _chains, \
        _non_null_chains, _atomic_structures, structure_datas, \
        _atoms_pair, _pseudobond_group_map

def _proper_dihedrals(p):
    return Proper_Dihedrals(p)
def _rotamers(p):
    return Rotamers(p)
def _proper_dihedrals_or_nones(p):
    return [Proper_Dihedral.c_ptr_to_py_inst(ptr) if ptr else None for ptr in p]
def _distance_restraints(p):
    return Distance_Restraints(p)
def _non_null_proper_dihedrals(p):
    return Proper_Dihedrals(p[p!=0])
def _atoms_four_tuple(p):
    return tuple((Atoms(p[:,i].copy()) for i in range(4)))
def _proper_dihedral_restraints(p):
    return Proper_Dihedral_Restraints(p)
def _rotamer_restraints(p):
    return Rotamer_Restraints(p)



class _Dihedrals(Collection):
    '''
    Base class for Proper_Dihedrals and Improper_Dihedrals. Do not
    instantiate directly.
    '''

    def __init__(self, c_pointers, single_class, array_class):
        super().__init__(c_pointers, single_class, array_class)

    atoms = cvec_property('dihedral_atoms', cptr, 4, astype=_atoms_four_tuple, read_only=True,
        doc = '''
        Returns a four-tuple of :class:`Atoms` objects. For each dihedral,
        its constituent atoms are in the matching position in the four
        :class:`Atoms` collections. Read only.
        ''')
    angles = cvec_property('dihedral_angle', float64, read_only=True,
        doc='Returns the angle in radians for each dihedral. Read only.')
    names = cvec_property('dihedral_name', string, read_only=True)


class Proper_Dihedrals(_Dihedrals):

    def __init__(self, c_pointers = None):
        super().__init__(c_pointers, Proper_Dihedral, Proper_Dihedrals)

    def delete(self):
        err_string = 'Dihedrals are not in charge of their own creation '\
            +'and deletion. Instead, use '\
            +'session.proper_dihedrals_mgr.delete_dihedrals(dihedrals).'
        raise RuntimeError(err_string)

    residues = cvec_property('dihedral_residue', cptr, astype=_residues, read_only=True,
        doc='Returns a :class:`Residues` giving the parent residue of each dihedral. Read only')
    axial_bonds = cvec_property('proper_dihedral_axial_bond', cptr, astype=_bonds, read_only=True,
        doc='Returns a :class:`Bonds` giving the axial bond for each dihedral. Read-only')

class Ramas(Collection):
    def __init__(self, c_pointers=None):
        super().__init__(c_pointers, Rama, Ramas)

    @property
    def omega_dihedrals(self):
        '''
        Returns a :class:`Proper_Dihedrals` instance pointing to the omega
        (peptide bond) dihedral for each residue. Note that some residues
        will not have omega dihedrals, so the output array may be shorter
        than the :class:`Ramas` instance.
        '''
        f = c_function('rama_omega',
            args = (ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p),
            ret = ctypes.c_size_t)
        n = len(self)
        ret = numpy.empty(n, cptr)
        found = f(self._c_pointers, n, pointer(ret))
        return _proper_dihedrals(ret[:found])

    @property
    def phi_dihedrals(self):
        '''
        Returns a :class:`Proper_Dihedrals` instance pointing to the phi
        dihedral for each residue. Note that some residues will not have phi
        dihedrals, so the output array may be shorter  than the :class:`Ramas`
        instance.
        '''
        f = c_function('rama_phi',
            args = (ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p),
            ret = ctypes.c_size_t)
        n = len(self)
        ret = numpy.empty(n, cptr)
        found = f(self._c_pointers, n, pointer(ret))
        return _proper_dihedrals(ret[:found])

    @property
    def psi_dihedrals(self):
        '''
        Returns a :class:`Proper_Dihedrals` instance pointing to the psi
        dihedral for each residue. Note that some residues will not have psi
        dihedrals, so the output array may be shorter  than the :class:`Ramas`
        instance.
        '''
        f = c_function('rama_psi',
            args = (ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p),
            ret = ctypes.c_size_t)
        n = len(self)
        ret = numpy.empty(n, cptr)
        found = f(self._c_pointers, n, pointer(ret))
        return _proper_dihedrals(ret[:found])

    residues = cvec_property('rama_residue', cptr, astype=_residues, read_only = True,
            doc = 'Returns a :class:`chimerax.Residues` instance giving the residue to which each Rama belongs. Read only.')
    ca_atoms = cvec_property('rama_ca_atom', cptr, astype=_atoms, read_only = True,
            doc = 'Returns a :class:`chimerax.Atoms` instance giving the alpha carbon of each amino acid residue. Read only.')
    valids = cvec_property('rama_is_valid', npy_bool, read_only = True,
            doc = 'True for each residue that has all three of omega, phi and psi. Read only.')
    visibles = cvec_property('rama_visible', npy_bool, read_only = True,
            doc = 'True for each residue whose alpha carbon is visible. Read only.')
    visibles_ignoring_ribbon = cvec_property('rama_only_hidden_by_ribbon', npy_bool, read_only=True,
            doc = 'True if the only thing hiding the alpha carbon is the ribbon display. Read only.')
    scores = cvec_property('rama_score', float64, read_only = True,
            doc = 'The score of each residue on the MolProbity Ramachandran contours. Read only.')
    phipsis = cvec_property('rama_phipsi', float64, 2, read_only = True,
            doc = 'The phi and psi angles for each residue in radians. Read only.')
    angles = cvec_property('rama_omegaphipsi', float64, 3, read_only = True,
            doc = 'The omega, phi and psi angles for each residue in radians. Read only.')
    cases = cvec_property('rama_case', uint8, read_only=True,
            doc = '''Values representing the Ramachandran case for these residues,
                matching the case definitions in :class:`Rama_Mgr.Rama_Case`. Read only.''')

class Rotamers(Collection):
    def __init__(self, c_pointers=None):
        super().__init__(c_pointers, Rotamer, Rotamers)

    residues = cvec_property('rotamer_residue', cptr, astype=_residues, read_only=True,
                doc=':py:class:`chimerax.Residue` this rotamer belongs to. Read only.')
    scores = cvec_property('rotamer_score', float64, read_only=True,
                doc='P-value for the current conformation of this rotamer. Read only.')
    ca_cb_bonds = cvec_property('rotamer_ca_cb_bond', cptr, astype=_bonds, read_only=True,
                doc='The "stem" :py:class:`chimerax.Bond` of this rotamer. Read only.')
    visibles = cvec_property('rotamer_visible', npy_bool, read_only=True,
                doc='True for each rotamer whose CA-CB bond is visible')

class Position_Restraints(Collection):
    def __init__(self, c_pointers=None, single_type = None, poly_type = None):
        if single_type is None:
            single_type = Position_Restraint
        if poly_type is None:
            poly_type = Position_Restraints
        super().__init__(c_pointers, single_type, poly_type)

    @property
    def _bond_cylinder_transforms(self):
        '''Transforms mapping a unit cylinder onto the restraint bonds. Read only.'''
        f = c_function('position_restraint_bond_transform',
            args = (ctypes.c_void_p, ctypes.c_size_t,
                ctypes.POINTER(ctypes.c_double)))
        n = len(self)
        transforms = empty((n,4,4), float64)
        f(self._c_pointers, n, pointer(transforms))
        return transforms

    def clear_sim_indices(self):
        '''
        Run at the end of a simulation to reset :attr:`sim_indices` to -1.
        '''
        f = c_function('position_restraint_clear_sim_index',
            args = (ctypes.c_void_p, ctypes.c_size_t))
        f(self._c_pointers, len(self))

    atoms = cvec_property('position_restraint_atom', cptr, astype=_atoms, read_only=True,
        doc = 'Returns the restrained :py:class:`chimerax.Atoms`. Read-only.')
    targets = cvec_property('position_restraint_target', float64, 3,
        doc = 'Target (x,y,z) positions in Angstroms. Can be written.')
    target_vectors = cvec_property('position_restraint_target_vector', float64, 3, read_only=True,
        doc = 'Returns the vectors ("bonds") connecting each atom to its target. Read only.')
    spring_constants = cvec_property('position_restraint_k', float64,
        doc = 'Restraint spring constants in :math:`kJ mol^{-1} nm^{-2}`. Can be written')
    enableds = cvec_property('position_restraint_enabled', npy_bool,
        doc = 'Enable/disable position restraints with a Numpy boolean array.')
    visibles = cvec_property('position_restraint_visible', npy_bool, read_only=True,
        doc = 'Returns a boolean mask giving the currently visible restraints. Read only.')
    sim_indices = cvec_property('position_restraint_sim_index', int32,
        doc = '''
        Index of each restraint in the associated force object in a running
        simulation. Restraints tha are not part of a simulation have indices of
        -1. Can be set, but only if you know what you are doing.
        ''')

class Tuggable_Atoms(Position_Restraints):
    def __init__(self, c_pointers=None):
        super().__init__(c_pointers, single_type=Tuggable_Atom, poly_type = Tuggable_Atoms)

class MDFF_Atoms(Collection):
    def __init__(self, c_pointers=None):
        super().__init__(c_pointers, MDFF_Atom, MDFF_Atoms)

    def clear_sim_indices(self):
        '''
        Run at the end of a simulation to reset :attr:`sim_indices` to -1.
        '''
        f = c_function('mdff_atom_clear_sim_index',
            args=(ctypes.c_void_p, ctypes.c_size_t))
        f(self._c_pointers, len(self))

    enableds = cvec_property('mdff_atom_enabled', npy_bool,
        doc='Enable/disable MDFF tugging on each atom or get the current states.')
    atoms = cvec_property('mdff_atom_atom', cptr, astype=_atoms, read_only=True,
        doc='Returns the :py:class:`chimerax.Atom`. Read only.')
    coupling_constants = cvec_property('mdff_atom_coupling_constant', float64,
        doc='''
        Per-atom MDFF coupling constants. These are multiplied by the global
        coupling constant for the map when calculating the MDFF potentials.
        Can be set.
        ''')
    sim_indices = cvec_property('mdff_atom_sim_index', int32,
        doc='''
        Index of each atom in the relevant MDFF Force in a running simulation.
        Atoms which are not currently in a simulation have indices equal to -1.
        Can be set, but only if you know what you are doing.
         ''')


class Distance_Restraints(Collection):
    def __init__(self, c_pointers=None):
        super().__init__(c_pointers, Distance_Restraint, Distance_Restraints)

    @property
    def _bond_cylinder_transforms(self):
        '''Transforms mapping a unit cylinder onto the restraint bonds. Read only.'''
        from chimerax.core.geometry import Places
        f = c_function('distance_restraint_bond_transform',
            args = (ctypes.c_void_p, ctypes.c_size_t,
                ctypes.POINTER(ctypes.c_double)))
        n = len(self)
        transforms = empty((n,4,4), float64)
        f(self._c_pointers, n, pointer(transforms))
        return Places(opengl_array=transforms)

    @property
    def _target_transforms(self):
        from chimerax.core.geometry import Places
        f = c_function('distance_restraint_target_transform',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_double)))
        n = len(self)
        transforms=empty((n,4,4), float64)
        f(self._c_pointers, n, pointer(transforms))
        return Places(opengl_array=transforms)

    def clear_sim_indices(self):
        '''
        Run at the end of a simulation to reset sim indices to -1
        '''
        f = c_function('distance_restraint_clear_sim_index',
            args = (ctypes.c_void_p, ctypes.c_size_t))
        f(self._c_pointers, len(self))

    enableds =cvec_property('distance_restraint_enabled', npy_bool,
            doc = 'Enable/disable these restraints or get their current states.')
    visibles = cvec_property('distance_restraint_visible', npy_bool, read_only = True,
            doc = 'Each restraint will be visible if it is enabled and both atoms are visible.')
    atoms = cvec_property('distance_restraint_atoms', cptr, 2, astype=_atoms_pair, read_only=True,
            doc = 'Returns a 2-tuple of :class:`Atoms` containing the restrained atoms. Read only.' )
    targets = cvec_property('distance_restraint_target', float64,
            doc = 'Target distances in Angstroms')
    spring_constants = cvec_property('distance_restraint_k', float64,
            doc = 'Restraint spring constants in :math:`kJ mol^{-1} nm^{-2}`')
    distances = cvec_property('distance_restraint_distance', float64, read_only=True,
            doc = 'Current distances between restrained atoms in Angstroms. Read only.')
    sim_indices = cvec_property('distance_restraint_sim_index', int32,
        doc='''
        Index of each restraint in the relevant MDFF Force in a running
        simulation. Restraints which are not currently in a simulation have
        indices equal to -1. Can be set, but only if you know what you are
        doing.
         ''')

class Proper_Dihedral_Restraints(Collection):
    def __init__(self, c_pointers=None):
        super().__init__(c_pointers, Proper_Dihedral_Restraint, Proper_Dihedral_Restraints)

    def _annotation_transforms(self):
        n = len(self)
        f = c_function('proper_dihedral_restraint_annotation_transform',
            args = (ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_void_p))
        tf1 = numpy.empty((n,4,4), float64)
        tf2 = numpy.empty((n,4,4), float64)
        f(self._c_pointers, n, pointer(tf1), pointer(tf2))
        from chimerax.core.geometry import Places
        return (Places(opengl_array=tf1), Places(opengl_array=tf2))

    def clear_sim_indices(self):
        f = c_function('proper_dihedral_restraint_clear_sim_index',
            args = (ctypes.c_void_p, ctypes.c_size_t))
        f(self._c_pointers, len(self))

    @property
    def atoms(self):
        return self.dihedrals.atoms

    targets = cvec_property('proper_dihedral_restraint_target', float64,
        doc = 'Target angles for each restraint in radians. Can be written.')
    dihedrals = cvec_property('proper_dihedral_restraint_dihedral', cptr, astype=_proper_dihedrals, read_only=True,
        doc = 'The restrained :py:class:`Proper_Dihedrals`. Read only.')
    offsets = cvec_property('proper_dihedral_restraint_offset', float64, read_only = True,
        doc = 'Difference between current and target angles in radians. Read only.')
    cutoffs = cvec_property('proper_dihedral_restraint_cutoff', float64,
        doc = 'Cutoff angle offsets below which no restraint will be applied. Can be set.')
    enableds = cvec_property('proper_dihedral_restraint_enabled', npy_bool,
        doc = 'Enable/disable each restraint or get their current states.')
    displays = cvec_property('proper_dihedral_restraint_display', npy_bool,
        doc = 'Set whether you want each restraint to be displayed when active.')
    visibles = cvec_property('proper_dihedral_restraint_visible', npy_bool, read_only=True,
        doc = 'Is each restraint currently visible? Read-only.')
    spring_constants = cvec_property('proper_dihedral_restraint_k', float64,
        doc = 'Get/set the spring constant for each restraint in :math:`kJ mol^{-1} rad^{-2}`')
    annotation_colors = cvec_property('proper_dihedral_restraint_annotation_color', uint8, 4, read_only=True,
        doc = 'Get the annotation color for each restraint according to the current colormap. Read only.')
    sim_indices = cvec_property('proper_dihedral_restraint_sim_index', int32,
        doc='''
        Index of each restraint in the relevant MDFF Force in a running
        simulation. Restraints which are not currently in a simulation have
        indices equal to -1. Can be set, but only if you know what you are
        doing.
         ''')

class Rotamer_Restraints(Collection):
    def __init__(self, c_pointers=None):
        super().__init__(c_pointers, Rotamer_Restraint, Rotamer_Restraints)

    def set_spring_constant(self, k):
        '''
        Sets the spring constants for all chi dihedrals in all rotamers to k.
        Write-only. To retrieve the current spring constants, use
        :attr:`Rotamer.chi_restraints.spring_constants` (note: this is only
        possible for one rotamer at a time).
        '''
        f = c_function('set_rotamer_restraint_spring_constant',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_double))
        f(self._c_pointers, len(self), k)


    rotamers = cvec_property('rotamer_restraint_rotamer', cptr, astype=_rotamers, read_only=True,
        doc = ':py:class:`Rotamers` to be restrained. Read only.')
    residues = cvec_property('rotamer_restraint_residue', cptr, astype=_residues, read_only=True,
        doc = ':py:class:`chimerax.Residues` to be restrained. Read only.')
    enableds = cvec_property('rotamer_restraint_enabled', npy_bool,
        doc = '''
        Enable/disable chi dihedral restraints. Returns False for any rotamer
        where at least one chi restraint is disabled.
        ''')
    target_index = cvec_property('rotamer_restraint_target_index', int32,
        doc = '''
        Get/set the index of the rotamer target definition giving target angles
        and cutoffs. If no restraint is currently applied, returns the last
        restraint that was applied to this rotamer. If no restraint has ever
        been applied, returns -1.
        ''')
