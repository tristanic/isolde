# @Author: Tristan Croll <tic20>
# @Date:   26-Apr-2018
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 27-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2017-2018 Tristan Croll



import os, sys, glob
import numpy
import ctypes
from enum import IntEnum

from chimerax.core.state import State
from chimerax.core.atomic import molc
from chimerax.core.atomic.molc import CFunctions, string, cptr, pyobject, \
    set_c_pointer, pointer, size_t
# object map lookups
from chimerax.core.atomic.molobject import _atoms, \
                _atom_pair, _atom_or_none, _bonds, _chain, _element, \
                _pseudobonds, _residue, _residues, _rings, _non_null_residues, \
                _residue_or_none, _residues_or_nones, _residues_or_nones, \
                _chains, _atomic_structure, _pseudobond_group, \
                _pseudobond_group_map

from chimerax.core.models import Model, Drawing

from numpy import int8, uint8, int32, uint32, float64, float32, byte, bool as npy_bool
import json

libdir = os.path.dirname(os.path.abspath(__file__))
libfile = glob.glob(os.path.join(libdir, 'molc.cpython*'))[0]
DICT_DIR = os.path.join(libdir, 'dictionaries')

_c_functions = CFunctions(os.path.splitext(libfile)[0])
c_property = _c_functions.c_property
cvec_property = _c_functions.cvec_property
c_function = _c_functions.c_function
c_array_function = _c_functions.c_array_function

def _asptr(arr, dtype):
    return arr.ctypes.data_as(pointer(dtype))

def _chiral_center_or_none(p):
    return Chiral_Center.c_ptr_to_py_inst(p) if p else None
def _chiral_centers(p):
    from .molarray import Chiral_Centers
    return Chiral_Centers(p)
def _proper_dihedrals(p):
    from .molarray import Proper_Dihedrals
    return Proper_Dihedrals(p)
def _ramas(p):
    from .molarray import Ramas
    return Ramas(p)
def _rama_or_none(p):
    return Rama.c_ptr_to_py_inst(p) if p else None
def _rotamers(p):
    from .molarray import Rotamers
    return Rotamers(p)
def _rotamer_or_none(p):
    return Rotamer.c_ptr_to_py_inst(p) if p else None
def _proper_dihedral_or_none(p):
    return Proper_Dihedral.c_ptr_to_py_inst(p) if p else None
def _rotamer_or_none(p):
    return Rotamer.c_ptr_to_py_inst(p) if p else None
def _bond_or_none(p):
    from chimerax.atomic import Bond
    return Bond.c_ptr_to_py_inst(p) if p else None
def _distance_restraint_or_none(p):
    return Distance_Restraint.c_ptr_to_py_inst(p) if p else None
def _distance_restraints(p):
    from .molarray import Distance_Restraints
    return Distance_Restraints(p)
def _position_restraint_or_none(p):
    return Position_Restraint.c_ptr_to_py_inst(p) if p else None
def _position_restraints(p):
    from .molarray import Position_Restraints
    return Position_Restraints(p)
def _tuggable_atom_or_none(p):
    return Tuggable_Atom.c_ptr_to_py_inst(p) if p else None
def _tuggable_atoms(p):
    from .molarray import Tuggable_Atoms
    return Tuggable_Atoms(p)
def _mdff_atom_or_none(p):
    return MDFF_Atom.c_ptr_to_py_inst(p) if p else None
def _mdff_atoms(p):
    from .molarray import MDFF_Atoms
    return MDFF_Atoms(p)
def _rotamer_restraint_or_none(p):
    return Rotamer_Restraint.c_ptr_to_py_inst(p) if p else None
def _rotamer_restraints(p):
    from .molarray import Rotamer_Restraints
    return Rotamer_Restraints(p)
def _pseudobond_or_none(p):
    from chimerax.core.atomic import Pseudobond
    return Pseudobond.c_ptr_to_py_inst(p) if p else None
def _chiral_restraint_or_none(p):
    return Chiral_Restraint.c_ptr_to_py_inst(p) if p else None
def _chiral_restraints(p):
    from .molarray import Chiral_Restraints
    return Chiral_Restraints(p)
def _chiral_restraint_mgr(p):
    return Chiral_Restraint_Mgr.c_ptr_to_existing_py_inst(p) if p else None
def _proper_dihedral_restraint_or_none(p):
    return Proper_Dihedral_Restraint.c_ptr_to_py_inst(p) if p else None
def _proper_dihedral_restraints(p):
    from .molarray import Proper_Dihedral_Restraints
    return Proper_Dihedral_Restraints(p)
def _proper_dihedral_restraint_mgr(p):
    return Proper_Dihedral_Restraint_Mgr.c_ptr_to_existing_py_inst(p) if p else None
def _distance_restraint_mgr(p):
    return Distance_Restraint_Mgr.c_ptr_to_existing_py_inst(p) if p else None
def _position_restraint_mgr(p):
    return Position_Restraint_Mgr.c_ptr_to_existing_py_inst(p) if p else None
def _tuggable_atoms_mgr(p):
    return Tuggable_Atoms_Mgr.c_ptr_to_existing_py_inst(p) if p else None
def _mdff_mgr(p):
    return MDFF_Mgr.c_ptr_to_existing_py_inst(p) if p else None
def _rotamer_restraint_mgr(p):
    return Rotamer_Restraint_Mgr.c_ptr_to_existing_py_inst(p) if p else None

def get_chiral_mgr(session):
    if hasattr(session, 'chiral_mgr') and not session.chiral_mgr.deleted:
        return session.chiral_mgr
    return Chiral_Mgr(session)

def get_proper_dihedral_manager(session):
    if hasattr(session, 'proper_dihedral_mgr') and not session.proper_dihedral_mgr.deleted:
        return session.proper_dihedral_mgr
    return Proper_Dihedral_Mgr(session)

def get_ramachandran_manager(session):
    if hasattr(session, 'rama_mgr') and not session.rama_mgr.deleted:
        return session.rama_mgr
    return Rama_Mgr(session)

def get_rotamer_manager(session):
    if hasattr(session, 'rota_mgr') and not session.rota_mgr.deleted:
        return session.rota_mgr
    return Rota_Mgr(session)

def _get_restraint_change_tracker(session):
    if hasattr(session, 'isolde_changes') and not session.isolde_changes.deleted:
        return session.isolde_changes
    return Restraint_Change_Tracker(session)

class _Dihedral_Mgr:
    '''Base class. Do not instantiate directly.'''
    def __init__(self, session, c_pointer=None):
        cname = type(self).__name__.lower()
        if c_pointer is None:
            new_func = cname + '_new'
            c_pointer = c_function(new_func, ret=ctypes.c_void_p)()
        set_c_pointer(self, c_pointer)
        f = c_function('set_'+cname+'_py_instance', args=(ctypes.c_void_p, ctypes.py_object))
        f(self._c_pointer, self)
        self.session = session

    def __delete__(self):
        self.delete()
        super().__delete__()

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    def delete(self):
        pass

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

class Chiral_Mgr(_Dihedral_Mgr):
    '''
    A session-level singleton managing all chiral centres. Rather than
    instantiating directly, it is best created/retrieved using
    :func:`session_extensions.get_chiral_mgr`.
    '''

    def __init__(self, session, c_pointer=None):
        super().__init__(session, c_pointer=c_pointer)
        if hasattr(session, 'chiral_mgr') and not session.chiral_mgr.deleted:
            raise RuntimeError('Session already has a chiral atoms manager!')
        session.chiral_mgr = self
        self._load_dict()

    def _load_dict(self):
        import json
        with open(os.path.join(DICT_DIR, 'chirals.json'), 'r') as f:
            cdict = self._chiral_dict = json.load(f)

        f = c_function('chiral_mgr_add_chiral_def',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_size_t,ctypes.c_size_t, ctypes.c_size_t,
                ctypes.c_double)
            )
        for resname, res_data in cdict.items():
            for center, cdata in res_data.items():
                self._add_chiral_def(resname, center, *cdata[0], cdata[1])

    @property
    def chiral_center_dict(self):
        return self._chiral_dict

    def _add_chiral_def(self, residue_name, chiral_atom_name, s1_names, s2_names,
        s3_names, expected_angle):
        f = c_function('chiral_mgr_add_chiral_def',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_size_t,ctypes.c_size_t, ctypes.c_size_t,
                ctypes.c_double)
            )
        rn_key = ctypes.py_object()
        rn_key.value = residue_name
        cn_key = ctypes.py_object()
        cn_key.value = chiral_atom_name
        sname_objs = []
        snums = []
        for sname_list in (s1_names, s2_names, s3_names):
            snums.append(len(sname_list))
            sname_objs.append(numpy.array(sname_list, numpy.object))

        f(self._c_pointer, ctypes.byref(rn_key), ctypes.byref(cn_key),
            *[pointer(n) for n in sname_objs], *snums, expected_angle)


    def add_chiral_def(self, residue_name, chiral_atom_name, s1_names, s2_names,
        s3_names, expected_angle):
        '''
        Add a definition to the dictionary of known chiral centres. The
        definition will only be valid for the current session. All instances of
        the chiral centre will be automatically found and tracked.

        Since chiral bonds can occur across residues, substituent atom names
        should be provided as lists of possible names.

        Args:
            * residue_name:
                - the 3-letter name of the residue the centre is found in
            * chiral_atom_name:
                - the name of the central chiral atom
            * s1_names:
                - a list of potential names for the highest priority substituent
            * s2_names:
                - a list of potential names for the second substituent
            * s3_names:
                - a list of potential names for the third substituent
            * expected_angle:
                - the expected dihedral angle (in radians) for the chiral centre
                  at equilibrium. While it is best to provide more accurate
                  measures based on an energy-minimised model of your compound,
                  fairly typical values are 0.6 for a (S) isomer and -0.6 for a
                  (R) isomer. Note that chirality restraints only become active
                  for deviations more than 15 degrees (0.26 radians), from the
                  expected angle, relying on the forcefield parameterisation
                  inside that range.
        '''
        cdict = self._chiral_dict
        try:
            # Raise an error if the entry already exists
            existing = cdict[residue_name][chiral_atom_name]
            raise AttributeError('A defintion already exists for this residue/atom name pair!')
        except KeyError:
            pass
        if residue_name not in cdict.keys():
            cdict[residue_name] = {}
        cdict[residue_name][chiral_atom_name] = [[s1_names, s2_names, s3_names], expected_angle]
        self._add_chiral_def(residue_name, chiral_atom_name, s1_names, s2_names,
            s3_names, expected_angle)

    def get_chiral(self, atom, create=True):
        '''
        Returns a :class:`Chiral_Center` for the given atom if it is a known
        chiral atom. Returns None if the atom is achiral, no definition exists
        for the chiral centre, required substituent atoms are missing, or
        create is False and the c++ object has not previously been created.

        Args:
            * atom:
                - a :class:`chimerax.Atom` instance
            * create:
                - if True, an attempt will be made to create the chiral centre
                  if it doesn't already exist.
        '''
        from chimerax.atomic import Atoms
        a = Atoms([atom])
        c = self.get_chirals(a, create)
        if len(c):
            return c[0]
        return None

    def get_chirals(self, atoms, create=True):
        '''
        Returns a :class:`Chiral_Centers` containing all known chiral centres in
        the given atoms.

        Args:
            * atoms:
                - a :class:`chimerax.Atoms` instance.
            * create:
                - if True, for every atom where a chiral definition exists but
                  no :cpp:class:`Chiral_Center` has been created, an attempt
                  will be made to create it. Otherwise, only previously-created
                  chirals will be returned.
        '''
        f = c_function('chiral_mgr_get_chiral',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
                ctypes.c_bool),
            ret=ctypes.py_object
        )
        n = len(atoms)
        return _chiral_centers(f(self._c_pointer, atoms._c_pointers, n, create))


    @property
    def num_mapped_chiral_centers(self):
        f = c_function('chiral_mgr_num_chirals',
            args=(ctypes.c_void_p,),
            ret=ctypes.c_size_t
        )
        return f(self._c_pointer)

    def delete_chirals(self, chirals):
        '''
        Delete the C++ objects for the given :class:`Chiral_Centers`. Note that
        this only deletes the chirality information, not any atoms/bonds. In
        general it should not be necessary to use this method - creation and
        deletion is fully automatic.
        '''
        f = c_function('chiral_mgr_delete_chiral',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p)
        )
        f(self._c_pointer, len(chirals), chirals._c_pointers)

    def delete(self):
        c_function('chiral_mgr_delete', args=(ctypes.c_void_p,))(self.cpp_pointer)
        delattr(self.session, 'chiral_mgr')



class Proper_Dihedral_Mgr(_Dihedral_Mgr):
    '''
    A session-level singleton managing all proper dihedrals (phi, psi, chi etc.).
    Rather than instantiating directly, it is best created/retrieved using
    :func:`session_extensions.get_proper_dihedral_mgr`.
    '''

    def __init__(self, session, c_pointer=None):
        super().__init__(session, c_pointer=c_pointer)
        if hasattr(session, 'proper_dihedral_mgr') and not session.proper_dihedral_mgr.deleted:
            raise RuntimeError('Session already has a proper dihedral manager!')
        session.proper_dihedral_mgr = self
        self._load_dict()

    def delete(self):
        c_function('proper_dihedral_mgr_delete', args=(ctypes.c_void_p,))(self.cpp_pointer)
        delattr(self.session, 'proper_dihedral_mgr')

    def delete_dihedrals(self, dihedrals):
        '''
        Delete all dihedrals in a :class:`Proper_Dihedrals`. Note that
        this will not affect the constituent atoms in any way, and
        should not actually be necessary in most cases. Dihedrals are
        automatically deleted at the C++ level when their manager or
        any of their constituent atoms are deleted. The deleted dihedrals
        will be automatically re-created if/when needed.
        '''
        f = c_function('proper_dihedral_mgr_delete_dihedral',
                args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p))
        f(self.cpp_pointer, len(dihedrals), dihedrals._c_pointers)

    @property
    def dihedral_dict(self):
        if not hasattr(self, '_dihedral_dict') or self._dihedral_dict is None:
            self._load_dict()
        return self._dihedral_dict

    def _load_dict(self):
        import json
        with open(os.path.join(DICT_DIR, 'named_dihedrals.json'), 'r') as f:
            dd = self._dihedral_dict = json.load(f)
        aa_resnames = dd['aminoacids']
        # Copy the definitions common to all amino acids to each individual
        # entry for easier lookup later
        rd = dd['residues']['protein']
        for aa_key in aa_resnames:
            for bd_key, data in dd['all_protein'].items():
                rd[aa_key][bd_key] = data
        f = c_function('proper_dihedral_mgr_add_dihedral_def',
                        args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                        ctypes.POINTER(ctypes.c_bool)))

        for res_key, res_data in rd.items():
            for d_key, d_data in res_data.items():
                if type(d_data) != list:
                    continue
                externals = numpy.zeros(4, numpy.bool)
                if type(d_data[1]) == list and d_data[1][0] in (0,1):
                    externals = numpy.array(d_data[1]).astype(numpy.bool)
                d_data = d_data[0]
                rk = ctypes.py_object()
                rk.value = res_key
                dk = ctypes.py_object()
                dk.value = d_key
                anames = numpy.array(d_data).astype(string)
                f(self._c_pointer, ctypes.byref(rk), ctypes.byref(dk),
                    pointer(anames), pointer(externals))

    def _reserve(self, n):
        '''
        Pre-allocate n spaces in the main map.
        '''
        f = c_function('proper_dihedral_mgr_reserve_map', args=(ctypes.c_void_p, ctypes.c_size_t))
        f(self._c_pointer, n)

    def create_all_dihedrals(self, residues):
        '''
        Create C++ objects for all known dihedrals in a set of residues. In
        general it is not necessary to do this, and preferable to allow them to
        simply be created as needed.

        Args:
            * residues:
                - A :class:`chimerax.Residues` instance
        '''
        dihedral_dict = self._dihedral_dict
        amino_acid_resnames = dihedral_dict['aminoacids']
        r = residues
        # self._reserve(len(r))
        aa_residues = r[numpy.in1d(r.names, amino_acid_resnames)]
        f = c_function('proper_dihedral_mgr_get_dihedrals',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_size_t, ctypes.c_bool),
            ret=ctypes.py_object)
        self._find_peptide_backbone_dihedrals(dihedral_dict, aa_residues, f)
        self._find_rotameric_dihedrals(dihedral_dict, amino_acid_resnames, aa_residues, f)

    def _find_peptide_backbone_dihedrals(self, dihedral_dict, aa_residues, f):
        for key in dihedral_dict['all_protein'].keys():
            k = ctypes.py_object()
            k.value = key
            f(self._c_pointer, aa_residues._c_pointers, ctypes.byref(k),
                len(aa_residues), True)

    def _find_rotameric_dihedrals(self, dihedral_dict, amino_acid_resnames, all_aa_residues, f):
        res_dict = dihedral_dict['residues']['protein']
        for aa in amino_acid_resnames:
            data = res_dict[aa]
            nchi = data['nchi']
            aa_residues = all_aa_residues[all_aa_residues.names == aa]
            for i in range(nchi):
                key = 'chi'+str(i+1)
                k = ctypes.py_object()
                k.value = key
                f(self._c_pointer, aa_residues._c_pointers, ctypes.byref(k),
                    len(aa_residues), True)

    def get_dihedral(self, residue, name, create=True):
        '''
        Retrieve a :class:`Proper_Dihedral` for the given residue and name,
        or None if no such dihedral exists.

        Args:
            * residue:
                - A :class:`chimerax.Residue` instance
            * name:
                - A string giving the lowercase name of the dihedral (e.g.
                  'phi', 'omega', 'chi1', etc.)
            * create (default: True):
                - If True, if the dihedral does not currently exist an attempt
                  will be made to create it.
        '''
        from chimerax.atomic import Residues
        r = Residues([residue])
        d = self.get_dihedrals(r, name, create)
        if len(d):
            return d[0]
        return None

    def get_dihedrals(self, residues, name, create = True):
        '''
        Returns a :class:`Proper_Dihedrals` providing the named dihedral
        (where it exists) for every residue in residues. The resulting
        array will be in the same order as residues, but may be shorter.

        Args:
            * residues:
                - A :class:`chimerax.Residues` instance
            * name:
                - A string giving the lowercase name of the dihedral (e.g.
                  'phi', 'omega', 'chi1', etc.)
            * create (default = True):
                - If True, if any dihedral does not currently exist an attempt
                  will be made to create it.
        '''
        f = c_function('proper_dihedral_mgr_get_dihedrals', args=(
                        ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                        ctypes.c_size_t, ctypes.c_bool),
                        ret=ctypes.py_object)
        n = len(residues)
        key = ctypes.py_object()
        key.value = name
        return _proper_dihedrals(f(self._c_pointer, residues._c_pointers, ctypes.byref(key), n, create))

    def get_all_dihedrals(self, residues):
        '''
        Returns a :class:`Proper_Dihedrals` containing all dihedrals that have
        been defined for the given residues. Any standard dihedrals (e.g.
        phi, psi, omega) will be created.

        Args:
            * residues:
                - a :class:`chimerax.Residues` instance
        '''
        self.create_all_dihedrals(residues)
        f = c_function('proper_dihedral_mgr_get_residue_dihedrals',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t),
            ret=ctypes.py_object)
        return _proper_dihedrals(f(self._c_pointer, residues._c_pointers, len(residues)))

    @property
    def num_mapped_dihedrals(self):
        '''
        Number of dihedrals currently being managed.
        '''
        f = c_function('proper_dihedral_mgr_num_mapped_dihedrals', args=(ctypes.c_void_p,), ret=ctypes.c_size_t)
        return f(self._c_pointer)

    def __len__(self):
        return self.num_mapped_dihedrals

class Rama_Mgr:
    '''
    Session-level singleton managing the Ramachandran scoring of protein
    residues. Rather than instantiating directly, it is best created/retrieved
    using :func:`session_extensions.get_ramachandran_mgr`.
    '''
    class Rama_Case(IntEnum):
        '''
        Enumerators for the different Ramachandran cases. These match
        an enumerator in the C++ layer, so don't change them unless you
        know *exactly* what you're doing.
        '''
        NONE=0
        CISPRO=1
        TRANSPRO = 2
        GLYCINE=3
        PREPRO=4
        ILEVAL=5
        GENERAL=6

    class Rama_Bin(IntEnum):
        '''
        Enumerator for validation bins. Values match an enumerator in the C++
        layer, so don't change them unless you know *exactly* what you're doing.
        '''
        FAVORED=0
        ALLOWED=1
        OUTLIER=2
        NA=-1

    from .validation.constants import validation_defaults as val_defaults

    RAMA_CASE_DETAILS = {
        Rama_Case.NONE: {
            'name': 'Not applicable',
            'file_prefix': None,
            'cutoffs': None
        },
        Rama_Case.CISPRO: {
            'name': 'Cis-proline residues',
            'file_prefix': 'rama8000-cispro',
            'cutoffs': [val_defaults.CISPRO_OUTLIER, val_defaults.CISPRO_ALLOWED]
        },
        Rama_Case.TRANSPRO: {
            'name': 'Trans-proline residues',
            'file_prefix': 'rama8000-transpro',
            'cutoffs': [val_defaults.TRANSPRO_OUTLIER, val_defaults.TRANSPRO_ALLOWED]
        },
        Rama_Case.GLYCINE: {
            'name': 'Glycine residues',
            'file_prefix': 'rama8000-gly-sym',
            'cutoffs': [val_defaults.GLYCINE_OUTLIER,val_defaults.GLYCINE_ALLOWED]
        },
        Rama_Case.PREPRO: {
            'name': 'Residues preceding proline',
            'file_prefix': 'rama8000-prepro-noGP',
            'cutoffs': [val_defaults.PREPRO_OUTLIER, val_defaults.PREPRO_ALLOWED]
        },
        Rama_Case.ILEVAL: {
            'name': 'Isoleucine or valine residues',
            'file_prefix': 'rama8000-ileval-nopreP',
            'cutoffs': [val_defaults.ILEVAL_OUTLIER, val_defaults.ILEVAL_ALLOWED]
        },
        Rama_Case.GENERAL: {
            'name': 'General amino acid residues',
            'file_prefix': 'rama8000-general-noGPIVpreP',
            'cutoffs': [val_defaults.GENERAL_OUTLIER, val_defaults.GENERAL_ALLOWED]
        }
    }

    def _prepare_all_validators(self):
        from .validation import validation
        data_dir = validation.get_molprobity_data_dir()
        for case, details in self.RAMA_CASE_DETAILS.items():
            file_prefix = details['file_prefix']
            if file_prefix is not None:
                file_prefix = os.path.join(data_dir, file_prefix)
                i_data = validation.generate_interpolator_data(file_prefix, True)
                self._add_interpolator(case, *i_data)

    def __init__(self, session, c_pointer=None):
        if hasattr(session, 'rama_mgr'):
            raise RuntimeError('Session already has a Ramachandran manager!')
        dmgr = self._dihedral_mgr = get_proper_dihedral_manager(session)
        self.session = session
        cname = type(self).__name__.lower()
        if c_pointer is None:
            new_func = cname + '_new'
            c_pointer = c_function(new_func, args=(ctypes.c_void_p,), ret=ctypes.c_void_p)(dmgr._c_pointer)
        set_c_pointer(self, c_pointer)
        f = c_function('set_'+cname+'_py_instance', args=(ctypes.c_void_p, ctypes.py_object))
        f(self._c_pointer, self)
        self._prepare_all_validators()
        self.set_default_cutoffs()
        self.set_default_colors()
        session.rama_mgr = self

    def delete(self):
        c_function('rama_mgr_delete', args=(ctypes.c_void_p,))(self._c_pointer)

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value
        delattr(self.session, 'rama_mgr')

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def set_default_cutoffs(self):
        '''
        Reset the Ramachandran cutoffs to default values.
        '''
        dd = self.RAMA_CASE_DETAILS
        for case, cd in dd.items():
            cutoffs = cd['cutoffs']
            if cutoffs is not None:
                self._set_cutoffs(case, *cutoffs)

    def _set_cutoffs(self, case, outlier, allowed):
        f = c_function('set_rama_mgr_cutoffs',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_double, ctypes.c_double))
        f(self._c_pointer, case, outlier, allowed)

    @property
    def cutoffs(self):
        '''
        Returns a dict giving the allowed and outlier cutoffs for each
        Ramachandran case. Read only.
        '''
        f = c_function('rama_mgr_cutoffs',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_double)))
        cdict = {}
        for case in self.Rama_Case:
            if case != Rama_Case.NONE:
                cutoffs = numpy.empty(2, numpy.double)
                f(self._c_pointer, case, pointer(cutoffs))
                cdict[case] = cutoffs
        return cdict

    def set_default_colors(self):
        '''
        Set the colours for visualisation of scores back to their defaults.
        '''
        from .validation.constants import validation_defaults as val_defaults
        self.set_color_scale(val_defaults.MAX_FAVORED_COLOR, val_defaults.ALLOWED_COLOR,
            val_defaults.OUTLIER_COLOR, val_defaults.NA_COLOR)

    def set_color_scale(self, max_c, mid_c, min_c, na_c):
        '''
        Define a custom colour scale for visualisation of Ramachandran scores.
        All arguments are iterables of four integers providing
        (red, green, blue, alpha) in the range (0..255).

        Args:
            * max_c:
                - colour associated with the maximum (most favourable) score
            * mid_c:
                - colour at the favoured/allowed cutoff
            * min_c:
                - colour at the allowed/outlier cutoff. All scores below the
                  outlier cutoff will have this colour
            * na_c:
                - colour to associate with residues that don't have
                  Ramachandran scores (e.g. N/C termini)
        '''
        f = c_function('rama_mgr_set_color_scale',
            args=(ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint8),
                  ctypes.POINTER(ctypes.c_uint8), ctypes.POINTER(ctypes.c_uint8),
                  ctypes.POINTER(ctypes.c_uint8)))
        maxc = numpy.array(max_c, uint8)
        midc = numpy.array(mid_c, uint8)
        minc = numpy.array(min_c, uint8)
        nac = numpy.array(na_c, uint8)
        for arr in (maxc, midc, minc, nac):
            if len(arr) != 4:
                raise TypeError('Each color should be an array of 4 values in the range (0,255)')
        f(self._c_pointer, pointer(maxc), pointer(midc), pointer(minc), pointer(nac))

    @property
    def color_scale(self):
        '''
        Returns the current colour scale as a 4-tuple of (max, mid, min, n/a)
        '''
        f = c_function('rama_mgr_get_color_scale',
            args=(ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint8),
                ctypes.POINTER(ctypes.c_uint8), ctypes.POINTER(ctypes.c_uint8),
                ctypes.POINTER(ctypes.c_uint8)))
        maxc = numpy.empty(4, uint8)
        midc = numpy.empty(4, uint8)
        minc = numpy.empty(4, uint8)
        nac = numpy.empty(4, uint8)
        f(self._c_pointer, pointer(maxc), pointer(midc), pointer(minc), pointer(nac))
        return (maxc, midc, minc, nac)

    @property
    def dihedral_manager(self):
        '''
        Returns the session :class:`Proper_Dihedral_Mgr` singleton.
        '''
        return self._dihedral_mgr

    def _add_interpolator(self, rama_case, ndim, axis_lengths, min_vals, max_vals, data):
        '''
        Create a RegularGridInterpolator for the given Ramachandran
        contours, and add it to the manager.
        Args:
            * rama_case:
                - An integer corresponding to the Ramachandran case
                  (see Rama_Mgr.Rama_Case for valid values)
            * axis_lengths:
                - A numpy int array giving the number of points along each axis
            * min_vals:
                - A numpy double array giving the minimum value for each axis
            * max_vals:
                - A numpy double array giving the maximum value for each axis
            * data:
                - A 2D numpy array containing the gridded data
        '''
        f = c_function('rama_mgr_add_interpolator',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_size_t,
                ctypes.POINTER(ctypes.c_uint32), ctypes.POINTER(ctypes.c_double),
                ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)))
        axis_lengths = axis_lengths.astype(uint32)
        f(self._c_pointer, rama_case, ndim, pointer(axis_lengths),
            pointer(min_vals), pointer(max_vals), pointer(data))

    def rama_cases(self, residues):
        '''
        Returns an array of integer values corresponding to the Ramachandran
        case enumerator :class:`Rama_Mgr.Rama_Cases`.

        Args:
            * residues:
                - A :class:`chimerax.Residues` instance
        '''
        f = c_function('rama_mgr_rama_case',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_uint8)))
        n = len(residues)
        ret = numpy.empty(n, uint8);
        f(self._c_pointer, residues._c_pointers, n, pointer(ret))
        return ret

    def bin_scores(self, scores, cases):
        '''
        Returns an array of integer values corresponding to the enumerator
        :class:`Rama_Mgr.Rama_Bin` to bin the given Ramachandran scores into
        favoured, allowed, outlier and N/A. The input arrays are the product
        of :func:`validate`.

        Args:
            * scores:
                - an array of floating-point scores
            * cases:
                - a matching integer array defining the Ramachandran cases
        '''
        f = c_function('rama_mgr_bin_scores',
            args=(ctypes.c_void_p, ctypes.POINTER(ctypes.c_double),
                ctypes.POINTER(ctypes.c_uint8), ctypes.c_size_t,
                ctypes.POINTER(ctypes.c_int32)))
        n = len(scores)
        if len(cases) != n:
            raise TypeError('Scores and cases arrays must be the same length!')
        bins = numpy.empty(n, int32)
        f(self._c_pointer, pointer(scores), pointer(cases), n, bins.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)))
        return bins

    def outliers(self, residues):
        '''
        Returns a :class:`chimerax.Residues` instance encompassing the
        subset of input residues that are Ramachandran outliers.

        Args:
            * residues:
                - a :class:`chimerax.Residues` instance
        '''
        scores, cases = self.validate_by_residue(residues)
        bins = self.bin_scores(scores, cases)
        return residues[bins==self.Rama_Bin.OUTLIER]

    def validate(self, residues_or_ramas):
        '''
        Returns an 2-tuple containing an array of Ramachandran scores and an
        array of case enum values for a set of residues or ramas. Residues
        lacking either phi or psi will have a score of -1.

        Args:
            * residues_or_ramas: either a :class:`chimerax.Residues` or
              :class:`Ramas` instance
        '''
        from .molarray import Ramas
        from chimerax.core.atomic import Residues
        if isinstance(residues_or_ramas, Ramas):
            return self._validate(residues_or_ramas)
        elif isinstance(residues_or_ramas, Residues):
            return self._validate_by_residue(residues_or_ramas)
        raise TypeError('input must be either a Ramas or a Residues object!')

    def _validate_by_residue(self, residues):
        '''
        Returns a tuple of (double, uint8) Numpy arrays giving the scores
        and Ramachandran cases for each residue. Non-Ramachandran
        (N- and C-terminal peptide and non-protein) residues will
        receive scores of -1. Case definitions are found in
        :class:`Rama_Mgr`.Rama_Case.
        '''
        f = c_function('rama_mgr_validate_by_residue',
            args=(ctypes.c_void_p, ctypes.c_void_p,ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_uint8)))
        n = len(residues)
        scores = numpy.empty(n, numpy.double)
        cases = numpy.empty(n, uint8)
        f(self._c_pointer, residues._c_pointers, n, pointer(scores), pointer(cases))
        return (scores, cases)

    def _validate(self, ramas):
        f = c_function('rama_mgr_validate',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_uint8)))
        n = len(ramas)
        scores = numpy.empty(n, float64)
        cases = numpy.empty(n, uint8)
        f(self._c_pointer, ramas._c_pointers,  n, pointer(scores), pointer(cases))
        return (scores, cases)

    def get_ramas(self, residues):
        '''
        Returns a :class:`Ramas` with one Rama for each protein residue in the
        list. Non-protein residues will be skipped so the length of the result
        may be different from the input array, but the returned Ramas will be
        in the same order as the protein residues in the input.

        Args:
            * residues:
                - a :class:`chimerax.Residues` instance
        '''
        f = c_function('rama_mgr_get_rama',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p),
            ret=ctypes.c_size_t)
        n = len(residues)
        ramas = numpy.empty(n, cptr)
        found = f(self._c_pointer, residues._c_pointers, n, pointer(ramas))
        return _ramas(ramas[0:found])

    def rama_colors(self, ramas):
        '''
        Returns a nx4 uint8 array giving a color for each rama corresponding
        to the current colormap.

        Args:
            * ramas:
                - a :class:`Ramas` instance
        '''
        f = c_function('rama_mgr_validate_and_color',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
                  ctypes.POINTER(ctypes.c_uint8)))
        n = len(ramas)
        colors = numpy.empty((n,4), uint8)
        f(self._c_pointer, ramas._c_pointers, n, pointer(colors))
        return colors

    def _ca_positions_colors_and_selecteds(self, ramas, hide_favored = False):
        '''
        Provides all the information necessary to draw a live visualisation of
        Ramachandran status overlaying alpha carbons in a single C++ call.

        Args:
            * ramas:
                - a :class: `Ramas` instance
            * hide_favored:
                - if True, only the data for non-favoured residues will be
                  returned.

        Returns:
            * Coordinates corresponding to CA positions
            * Colors for each position
            * A Boolean mask corresponding to selection state of CAs (so
              the corresponding positions in the rama drawing can be set
              accordingly)
        '''
        f = c_function('rama_mgr_ca_positions_and_colors',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
                ctypes.c_bool, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p),
                ret = ctypes.c_size_t)
        n = len(ramas)
        coords = numpy.empty((n,3), float64)
        colors = numpy.empty((n, 4), uint8)
        selecteds = numpy.empty(n, npy_bool)
        count = f(self._c_pointer, ramas._c_pointers, n, hide_favored,
            pointer(coords), pointer(colors), pointer(selecteds))
        return (coords[0:count], colors[0:count], selecteds[0:count])

    def color_cas_by_rama_score(self, ramas, hide_favored = False):
        '''
        Colours the alpha carbon atoms by Ramachandran score.

        Args:
            * ramas:
                - a :class:`Ramas` instance
            * hide_favored:
                - if True, only the data for non-favoured residues will be
                  returned.
        '''
        f = c_function('rama_mgr_validate_and_color_cas',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_bool))
        n = len(ramas)
        f(self._c_pointer, ramas._c_pointers, n, hide_favored)

    def _draw_cis_and_twisted_omegas(self, ramas):
        '''
        Provides the geometry and colour information necessary to provide a
        drawing annotating cis/twisted peptide bonds (by filling in the "cup"
        formed by C-CA-N-C).

        Args:
            * ramas:
                - a :class:`Ramas` instance
        '''
        f = c_function('rama_mgr_draw_cis_and_twisted_omegas',
            args = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
                ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p),
            ret = ctypes.c_size_t)
        n = len(ramas)
        vertices = numpy.empty((n*5,3), float64)
        normals = numpy.empty((n*5,3), float64)
        triangles = numpy.empty((n*3,3), int32)
        colors = numpy.empty((n*5,4), uint8)
        count = f(self._c_pointer, ramas._c_pointers, n, pointer(vertices),
            pointer(normals), pointer(triangles), pointer(colors))
        return vertices[0:count*5], normals[0:count*5], triangles[0:count*3], colors[0:count*5]

    #######
    # Access to the underlying interpolator data
    #######
    def interpolator_dim(self, rama_case):
        '''
        Retrieve the number of dimensions in the :class:`RegularGridInterpolator`
        object for a given Ramachandran case. Should always return 2.

        Args:
            * rama_case:
                - integer value corresponding to the :class:`Rama_Cases` enum
        '''
        f = c_function('rama_mgr_interpolator_dim',
            args=(ctypes.c_void_p, ctypes.c_size_t),
            ret=ctypes.c_size_t)
        return f(self._c_pointer, rama_case)

    def interpolator_axis_lengths(self, rama_case):
        '''
        Retrieve the (phi,psi) axis dimensions of the
        :class:`RegularGridInterpolator` object for a given Ramachandran case.

        Args:
            * rama_case:
                - integer value corresponding to the :class:`Rama_Cases` enum
        '''
        dim = self.interpolator_dim(rama_case)
        f = c_function('rama_mgr_interpolator_axis_lengths',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_uint32)))
        ret = numpy.empty(dim, uint32)
        f(self._c_pointer, rama_case, pointer(ret))
        return ret

    def interpolator_limits(self, rama_case):
        '''
        Returns a (minimum_values, maximum_values) tuple giving the limiting
        values for each axis in the :class:`RegularGridInterpolator` object
        for a given Ramachandran case.

        Args:
            * rama_case:
                - integer value corresponding to the :class:`Rama_Cases` enum
        '''
        dim = self.interpolator_dim(rama_case)
        f = c_function('rama_mgr_interpolator_minmax',
            args=(ctypes.c_void_p, ctypes.c_size_t,
                ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)))
        minvals = numpy.empty(dim, float64)
        maxvals = numpy.empty(dim, float64)
        f(self._c_pointer, rama_case, pointer(minvals), pointer(maxvals))
        return (minvals, maxvals)

    def interpolator_values(self, rama_case):
        '''
        Returns a multidimensional array containing the contour data for a
        given Ramachandran case.

        Args:
            * rama_case:
                - integer value corresponding to the :class:`Rama_Cases` enum
        '''
        shape = self.interpolator_axis_lengths(rama_case)
        f = c_function('rama_mgr_interpolator_values',
            args=(ctypes.c_void_p, ctypes.c_size_t,
                ctypes.POINTER(ctypes.c_double)))
        data = numpy.empty(shape, float64)
        f(self._c_pointer, rama_case, pointer(data))
        return data

    def interpolator_axes(self, rama_case):
        '''
        Convenience function combining :func:`interpolator_axis_lengths` and
        :func:`interpolator_limits` to provide a tuple of arrays giving the
        axis values at each grid point.

        Args:
            * rama_case:
                - integer value corresponding to the :class:`Rama_Cases` enum
        '''
        lengths = self.interpolator_axis_lengths(rama_case)
        minmax = self.interpolator_limits(rama_case)
        axes = [numpy.linspace(minmax[0][i], minmax[1][i], lengths[i]) for i in range(len(lengths))]
        return tuple(axes)

class Rota_Mgr:
    '''
    Session-level singleton managing rotamers and their scoring. Rather than
    instantiating directly, it is best created/retrieved using
    :func:`session_extensions.get_rotamer_mgr`.
    '''
    class Rota_Bin(IntEnum):
        '''
        Enumerator for validation bins. Values match an enumerator in the C++
        layer, so don't change them unless you know *exactly* what you're doing.
        '''
        FAVORED=0
        ALLOWED=1
        OUTLIER=2
        NA=-1

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def delete(self):
        c_function('rota_mgr_delete', args=(ctypes.c_void_p,))(self.cpp_pointer)
        delattr(self.session, 'proper_dihedral_mgr')


    def __init__(self, session, c_pointer=None):
        if hasattr(session, 'rota_mgr'):
            raise RuntimeError('Session already has a Rotamer manager!')
        self._dihedral_mgr = get_proper_dihedral_manager(session)
        cname = type(self).__name__.lower()
        if c_pointer is None:
            new_func = cname + '_new'
            c_pointer = c_function(new_func, args=(ctypes.c_void_p,), ret=ctypes.c_void_p)(self._dihedral_mgr._c_pointer)
        set_c_pointer(self, c_pointer)
        f = c_function('set_'+cname+'_py_instance', args=(ctypes.c_void_p, ctypes.py_object))
        f(self._c_pointer, self)
        self.session = session
        self._prepare_all_validators()
        self._load_rotamer_defs()
        self._load_target_defs()
        self.set_default_cutoffs()
        self.set_default_colors()
        session.rota_mgr = self

    def _load_target_defs(self, degrees=True):
        '''
        Load rotamer targets from their JSON file and store them in order of
        decreasing prevalence for each residue.

        Args:
            * degrees:
                - tells the function whether the values in the JSON file are in
                  degrees or radians
        '''
        func=c_function('rota_mgr_add_target_def',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
                ctypes.c_void_p, ctypes.POINTER(ctypes.c_double),
                ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)))
        with open(os.path.join(DICT_DIR, 'rota_data.json'), 'rt') as f:
            rd = json.load(f)
        for aa, data in rd.items():
            rdata = [(name, d) for (name, d) in data.items()]
            rdata = sorted(rdata, key=lambda d: d[1]['freq'], reverse=True)
            names = numpy.array([d[0] for d in rdata], numpy.object)
            n = len(names)
            frequencies = numpy.array([d[1]['freq'] for d in rdata], float64)
            angles = numpy.concatenate(
                numpy.array([d[1]['angles'] for d in rdata], float64)
            )
            esds = numpy.concatenate(
                numpy.array([d[1]['esds'] for d in rdata], float64)
            )
            if degrees:
                angles=numpy.radians(angles)
                esds=numpy.radians(esds)
            nkey = ctypes.py_object()
            nkey.value = aa
            func(self._c_pointer, ctypes.byref(nkey), n, pointer(names),
                pointer(frequencies), pointer(angles), pointer(esds))

    def get_rota_targets(self, resname):
        '''
        Returns an OrderedDict giving rotamer name, angles, esds and frequencies
        sorted in order of decreasing frequency.

        Args:
            * resname
                - three-letter element_name, in capitals
        '''
        from copy import copy
        rd = self._rota_targets
        if resname in rd.keys():
            return copy(rd[resname])
        return None

    def set_default_colors(self):
        '''
        Reset the colour map for visualisation of rotamer validation back to
        the stored default colours.
        '''
        from .validation.constants import validation_defaults as val_defaults
        self.set_color_scale(val_defaults.MAX_FAVORED_COLOR, val_defaults.ALLOWED_COLOR,
            val_defaults.OUTLIER_COLOR)

    def set_color_scale(self, max_c, mid_c, min_c):
        '''
        Define a custom colour scale for visualisation of rotamer scores.
        All arguments are iterables of four integers providing
        (red, green, blue, alpha) in the range (0..255).

        Args:
            * max_c:
                - colour associated with the maximum (most favourable) score
            * mid_c:
                - colour at the favoured/allowed cutoff
            * min_c:
                - colour at the allowed/outlier cutoff. All scores below the
                  outlier cutoff will have this colour
        '''
        f = c_function('set_rota_mgr_color_scale',
            args=(ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint8),
                  ctypes.POINTER(ctypes.c_uint8), ctypes.POINTER(ctypes.c_uint8)))
        maxc = numpy.array(max_c, uint8)
        midc = numpy.array(mid_c, uint8)
        minc = numpy.array(min_c, uint8)
        for arr in (maxc, midc, minc):
            if len(arr) != 4:
                raise TypeError('Each color should be an array of 4 values in the range (0,255)')
        f(self._c_pointer, pointer(maxc), pointer(midc), pointer(minc))

    @property
    def color_scale(self):
        '''
        Returns the current colour scale as a 3-tuple of (max, mid, min)
        '''
        f = c_function('rota_mgr_color_scale',
            args=(ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint8),
                ctypes.POINTER(ctypes.c_uint8), ctypes.POINTER(ctypes.c_uint8)))
        maxc = numpy.empty(4, uint8)
        midc = numpy.empty(4, uint8)
        minc = numpy.empty(4, uint8)
        f(self._c_pointer, pointer(maxc), pointer(midc), pointer(minc))
        return (maxc, midc, minc)

    def set_default_cutoffs(self):
        '''
        Reset the rotamer P-value cutoffs to default values
        '''
        from .validation.constants import validation_defaults as vc
        self._set_cutoffs(vc.ROTA_ALLOWED_CUTOFF, vc.ROTA_OUTLIER_CUTOFF)

    def _set_cutoffs(self, allowed, outlier):
        f = c_function('rota_mgr_set_cutoffs',
            args=(ctypes.c_void_p, ctypes.c_double, ctypes.c_double))
        f(self._c_pointer, allowed, outlier)

    @property
    def cutoffs(self):
        '''
        Gives the current (allowed, outlier) P-value cutoffs
        '''
        f = c_function('rota_mgr_get_cutoffs',
            args=(ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)))
        cutoffs = numpy.empty(2, numpy.double)
        f(self._c_pointer, pointer(cutoffs))
        return cutoffs

    # TODO: provide a method retrieving these from the C++ layer
    # @property
    # def defined_rotamers(self):
    #     if not hasattr(self, '_defined_rotamer_dict') or self._defined_rotamer_dict is None:
    #         self._load_rotamer_defs()
    #     return self._defined_rotamer_dict

    def _load_defined_rotamers(self):
        with open(os.path.join(DICT_DIR, 'rota_data.json'), 'r') as f:
            self._defined_rotamer_dict = json.load(f)

    def _load_rotamer_defs(self):
        f = c_function('rota_mgr_add_rotamer_def',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
            ctypes.c_size_t, ctypes.c_bool, ctypes.POINTER(ctypes.c_uint8),
            ctypes.POINTER(ctypes.py_object)))
        dd = self._dihedral_mgr.dihedral_dict
        pd = dd['residues']['protein']
        aa_names = dd['aminoacids']
        for aa in aa_names:
            pda = pd[aa]
            nchi = pda['nchi'] # Number of chi dihedrals
            if nchi > 0:
                val_nchi = pda['val_nchi'] # Number of chi dihedrals used for validation
                symm = pda['symm']
                key = ctypes.py_object()
                key.value = aa
                moving_atom_counts = []
                moving_atom_names = []
                for i in range(1, nchi+1):
                    ckey = 'chi'+str(i)
                    ma = pda[ckey][1]
                    moving_atom_counts.append(len(ma))
                    moving_atom_names.extend(ma)
                    mac = numpy.array(moving_atom_counts, numpy.uint8)
                    man = numpy.array(moving_atom_names, numpy.object)
                f(self._c_pointer, ctypes.byref(key), nchi, val_nchi, symm,
                    pointer(mac), pointer(man))

    def _prepare_all_validators(self):
        from .validation import validation
        dmgr = self._dihedral_mgr
        data_dir = validation.get_molprobity_data_dir()
        prefix = os.path.join(data_dir, 'rota8000-')
        for aa in dmgr.dihedral_dict['aminoacids']:
            fname = prefix + aa.lower()
            if not os.path.isfile(fname+'.data') and not os.path.isfile(fname+'.pickle'):
                # Not a rotameric residue
                continue
            idata = validation.generate_interpolator_data(fname, True)
            self._add_interpolator(aa, *idata)

    def _add_interpolator(self, resname, ndim, axis_lengths, min_vals, max_vals, data):
        f = c_function('rota_mgr_add_interpolator',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_uint32), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)))
        key = ctypes.py_object()
        key.value = resname
        axis_lengths = axis_lengths.astype(uint32)
        f(self._c_pointer, ctypes.byref(key), ndim, pointer(axis_lengths),
            pointer(min_vals), pointer(max_vals), pointer(data))

    def get_rotamer(self, residue):
        '''
        Create/retrieve the :class:`Rotamer` object for a given residue

        Args:
            * residue:
                - a :class:`chimerax.Residue` instance
        '''
        from chimerax.core.atomic import Residues
        rots = self.get_rotamers(Residues([residue]))
        if len(rots):
            return rots[0]
        return None

    def get_rotamers(self, residues):
        '''
        Return a :class:`Rotamers` instance containing all valid rotamers in the
        selection.

        Args:
            * residues:
                - a :class:`chimerax.Residues` instance
        '''
        f = c_function('rota_mgr_get_rotamer',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t),
            ret=ctypes.py_object)
        n = len(residues)
        return _rotamers(f(self._c_pointer, residues._c_pointers, n))

    def validate_rotamers(self, rotamers):
        '''
        Returns an array of P-values for the current conformations of the
        given rotamers.

        Args:
            * rotamers:
                - a :class:`Rotamers` instance
        '''
        f = c_function('rota_mgr_validate_rotamer',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_double)))
        n = len(rotamers)
        ret = numpy.empty(n, numpy.double)
        f(self._c_pointer, rotamers._c_pointers, n, pointer(ret))
        return ret

    def validate_residues(self, residues):
        '''
        Returns an array of P-values for the current conformations of all
        rotameric protein residues in the input. Non-rotameric residues will
        receive a score of -1.

        Args:
            * residues:
                - a :class:`chimerax.Residues` instance
        '''
        f = c_function('rota_mgr_validate_residue',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_double)))
        n = len(residues)
        ret = numpy.empty(n, numpy.double)
        f(self._c_pointer, residues._c_pointers, n, pointer(ret))
        return ret

    def non_favored_rotamers(self, rotamers):
        '''
        Returns a 2-tuple containing only the rotamers in non-favoured
        conformations, and their current scores.

        Args:
            * rotamers:
                - a :class:`Rotamers` instance
        '''
        f = c_function('rota_mgr_non_favored',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)),
            ret=ctypes.c_size_t)
        n = len(rotamers)
        ptrs = numpy.empty(n, cptr)
        scores = numpy.empty(n, float64)
        found = f(self._c_pointer, rotamers._c_pointers, n, pointer(ptrs), pointer(scores))
        return (_rotamers(ptrs[0:found]), scores[0:found])

    def validate_scale_and_color_rotamers(self, rotamers, max_scale = 2.0, non_favored_only = True, visible_only = True):
        '''
        Used by :class:`Rotamer_Annotator` for visualising rotamer validation.

        Args:
            * rotamers:
                - a :class:`Rotamers` instance
            * max_scale:
                - size limit when scaling indicators by outlier severity
            * non_favored_only:
                - if True, the return arrays will be limited to only those
                  rotamers outside of favoured conformations
            * visible_only:
                - if True, non-visible residues will be excluded from analysis
                  and return

        Returns:
            * A :class:`Rotamers` instance matching the filtering criteria
            * An array of scale factors (one per rotamer)
            * An array of colours (one per rotamer)
        '''
        f = c_function('rota_mgr_validate_scale_and_color_rotamers',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
                ctypes.c_double, ctypes.c_bool, ctypes.c_bool, ctypes.c_void_p,
                ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_uint8)),
                ret=ctypes.c_size_t)
        n = len(rotamers)
        rot_out = numpy.empty(n, cptr)
        scale_out = numpy.empty(n, float64)
        color_out = numpy.empty((n,4), uint8)
        count = f(self._c_pointer, rotamers._c_pointers, n, max_scale,
            non_favored_only, visible_only, pointer(rot_out), pointer(scale_out), pointer(color_out))
        return (_rotamers(rot_out[0:count]), scale_out[0:count], color_out[0:count])


class Restraint_Change_Tracker:
    '''
    A per-session singleton tracking changes in ISOLDE restraints, and firing
    triggers as necessary. It should very rarely be necessary to work with
    this class directly. Each individual :class:`_Restraint_Mgr` subclass
    instance has its own :class:`TriggerSet`, and :class:`Restraint_Change_Tracker`
    will fire the 'changes' trigger in the appropriate instances whenever
    necessary.
    '''
    _mgr_name_to_class_functions = {
        'Chiral_Restraint_Mgr': (_chiral_restraint_mgr, _chiral_restraints),
        'Proper_Dihedral_Restraint_Mgr': (_proper_dihedral_restraint_mgr, _proper_dihedral_restraints),
        'Position_Restraint_Mgr': (_position_restraint_mgr, _position_restraints),
        'Distance_Restraint_Mgr': (_distance_restraint_mgr, _distance_restraints),
        'Tuggable_Atoms_Mgr': (_tuggable_atoms_mgr, _tuggable_atoms),
        'MDFF_Mgr': (_mdff_mgr, _mdff_atoms),
        'Rotamer_Restraint_Mgr': (_rotamer_restraint_mgr, _rotamer_restraints),
    }
    def __init__(self, session, c_pointer=None):
        cname = 'change_tracker'
        if c_pointer is None:
            if hasattr(session, 'isolde_changes'):
                if not session.isolde_changes.deleted():
                    c_pointer = session.isolde_changes._c_pointer
            if c_pointer is None:
                new_func = cname + '_new'
                c_pointer = c_function(new_func, ret=ctypes.c_void_p)()
        set_c_pointer(self, c_pointer)
        f = c_function('set_'+cname+'_py_instance', args=(ctypes.c_void_p, ctypes.py_object))
        f(self._c_pointer, self)
        self.session = session
        session.isolde_changes = self
        self._update_handler = session.triggers.add_handler('new frame', self._get_and_clear_changes)

    def delete(self):
        self.session.triggers.remove_handler(self._update_handler)
        c_function('change_tracker_delete', args=(ctypes.c_void_p,))(self._c_pointer)

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value
        delattr(self.session, 'rama_mgr')

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def clear(self):
        f = c_function('change_tracker_clear',
            args = (ctypes.c_void_p,))
        f(self._c_pointer)

    def _get_and_clear_changes(self, *_):
        self._get_and_process_changes()
        self.clear()

    def _process_changes(self, changes):
        processed_dict = {}
        for mgr_name, type_dict in changes.items():
            processed_type_dict = processed_dict[mgr_name] = {}
            class_funcs = self._mgr_name_to_class_functions[mgr_name]
            for mgr_ptr, changeds in type_dict.items():
                mgr = class_funcs[0](mgr_ptr)
                processed_changeds = processed_type_dict[mgr] = {}
                for change_type, changed_ptrs in changeds.items():
                    changed_obj = class_funcs[1](changed_ptrs)
                    processed_changeds[change_type] = changed_obj
                mgr.triggers.activate_trigger('changes', (mgr, processed_changeds))
        return processed_dict

    def _get_and_process_changes(self):
        f = c_function('change_tracker_changes',
            args = (ctypes.c_void_p,),
            ret = ctypes.py_object)
        return self._process_changes(f(self._c_pointer))

    @property
    def reason_names(self):
        f= c_function('change_tracker_reason_names',
            args=(ctypes.c_void_p,),
            ret = ctypes.py_object)
        return f(self._c_pointer)

class _Restraint_Mgr(Model):
    '''Base class. Do not instantiate directly.'''
    def __init__(self, name, model, c_pointer=None, c_class_name = None):
        session = model.session
        ct = self._change_tracker = _get_restraint_change_tracker(session)

        if c_class_name is None:
            cname = type(self).__name__.lower()
        else:
            cname = c_class_name
        if c_pointer is None:
            new_func = cname + '_new'
            c_pointer = c_function(
                new_func, args=(ctypes.c_void_p, ctypes.c_void_p),
                          ret=ctypes.c_void_p)(model._c_pointer, ct._c_pointer)
        set_c_pointer(self, c_pointer)
        f = c_function('set_'+cname+'_py_instance', args=(ctypes.c_void_p, ctypes.py_object))
        f(self._c_pointer, self)
        super().__init__(name, session)
        self.pickable = False
        self.model = model
        model.add([self])

    @property
    def triggers(self):
        '''
        A :py:class:`chimerax.TriggerSet` instance. Contains a single trigger,
        'changes' that is fired every time a restraint is created, deleted,
        or changes in any way. To automatically run a method when the restraints
        in a given restraint manager class change, use the following idiom:

        .. code:: python

            def changes_cb(trigger_name, data):
                mgr = data[0]
                changes_dict = data[1]
                reasons = list(changes_dict.keys())
                if 'target changed' in reasons:
                    changed_restraints = changes_dict['target changed']
                    do_something_with(changed_restraints)

            handler = mgr.triggers.add_handler('changes', changes_cb)
            # changes_cb() will be run every time something changes.

            # To stop:
            mgr.trigger.remove_handler(handler)

        '''
        if not hasattr(self, '_triggers') or self._triggers is None:
            from chimerax.core.triggerset import TriggerSet
            t = self._triggers = TriggerSet()
            t.add_trigger('changes')
        return self._triggers

    def delete(self):
        cname = type(self).__name__.lower()
        del_func = cname+'_delete'
        c_function(del_func, args=(ctypes.c_void_p,))(self._c_pointer)
        super().delete()

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    @property
    def display(self):
        return Model.display.fget(self)

    @display.setter
    def display(self, flag):
        cflag = self.display
        Model.display.fset(self, flag)
        if flag and not cflag:
            self.update_graphics()

    def update_graphics(self):
        ''' Should be overridden in derived classes. '''
        pass

class MDFF_Mgr(_Restraint_Mgr):
    '''
    Manages the molecular dynamics flexible fitting (MDFF) steering forces for
    one map. Appears as a Model to provide a consistent API to all the other
    restraint manager classes, but doesn't actually draw anything. Unlike other
    restraint manager classes, it lives as a child to the
    :py:class:`chimerax.Volume` holding the map it manages.

    The preferred way to create/retrieve the MDFF manager for a given
    :py:class:`AtomicStructure` instance m and :py:class:`Volume` instance is:

    .. code:: python

        from chimerax.isolde import session_extensions as sx
        mdff_mgr = sx.get_mdff_mgr(m, v)
    '''
    def __init__(self, model, volume, c_pointer=None):
        name = 'MDFF - ' + volume.name
        super().__init__(name, model, c_pointer=c_pointer, c_class_name = "mdff_mgr")
        self._volume = volume
        # Place as sub-model to the Volume object so deleting the Volume automatically
        # deletes the MDFF manager
        self.volume.add([self])
        if 'global k changed' not in self.triggers.trigger_names():
            self.triggers.add_trigger('global k changed')

    @property
    def volume(self):
        '''
        The Volume object providing the MDFF potential for the managed atoms.
        '''
        return self._volume

    def _get_mdff_atoms(self, atoms, create=False):
        f = c_function('mdff_mgr_get_mdff_atom',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_bool,
                ctypes.c_void_p),
            ret = ctypes.c_size_t)
        n = len(atoms)
        ret = numpy.empty(n, cptr)
        count = f(self._c_pointer, atoms._c_pointers, n, create, pointer(ret))
        return _mdff_atoms(ret[:count])

    def add_mdff_atoms(self, atoms, hydrogens=False):
        '''
        Returns a :py:class:`MDFF_Atoms` instance covering the atoms in the
        input. Only one :cpp:class:`MDFF_Atom` will ever be created for a given
        atom, so it is safe to use this repeatedly on the same set of atoms.
        Use this method when you want to ensure that all atoms in the input will
        be properly coupled to the map in simulations.

        Args:
            * atoms:
                - a :py:class:`chimerax.Atoms` instance
            * hydrogens:
                - should hydrogens be coupled to the map? Setting this to True
                  will substantially increase the amount of computation required
                  and empirically seems to make little difference to the result.
        '''
        if not hydrogens:
            atoms = atoms[atoms.element_names != "H"]
        return self._get_mdff_atoms(atoms, create=True)

    def get_mdff_atoms(self, atoms):
        '''
        Returns a :py:class:`MDFF_Atoms` instance covering the atoms in the
        input. Unlike :func:`add_mdff_atoms`, no new :cpp:class:`MDFF_Atom` will
        be created. Use this method when you've already coupled all the atoms
        you want to this map, and don't want to accidentally add more.

        Args:
            * atoms:
                - a :py:class:`chimerax.Atoms` instance
        '''
        return self._get_mdff_atoms(atoms, create=False)

    def __len__(self):
        f = c_function('mdff_mgr_num_atoms',
            args=(ctypes.c_void_p,),
            ret = ctypes.c_size_t)
        return f(self._c_pointer)

    def _get_global_k(self):
        '''
        Global coupling constant applied to all atoms (scales the individual
        atom coupling constants). Effective units are
        :math:`kJ mol^{-1} (\\text{map density unit})^{-1} nm^3`
        '''
        f = c_function('mdff_mgr_global_k',
            args=(ctypes.c_void_p,),
            ret=ctypes.c_double)
        return f(self._c_pointer)

    def _set_global_k(self, k):
        f = c_function('set_mdff_mgr_global_k',
            args=(ctypes.c_void_p, ctypes.c_double))
        f(self._c_pointer, k)
        self.triggers.activate_trigger('global k changed', (self, k))

    global_k = property(_get_global_k, _set_global_k)


class Position_Restraint_Mgr(_Restraint_Mgr):
    '''
    Manages creation, deletion, mapping and drawing of position restraints for
    a single atomic structure. Appears as a child :py:class:`chimeraX.Model`
    under the :py:class:`chimerax.AtomicStructure` it manages.

    The preferred way to create/retrieve the position restraint manager for a
    given :py:class:`AtomicStructure` instance m is:

    .. code:: python

        from chimerax.isolde import session_extensions as sx
        pr_mgr = sx.get_position_restraint_mgr(m)
    '''
    _DEFAULT_BOND_COLOR = [200, 250, 120, 255]
    _DEFAULT_PIN_COLOR = [255,215,0,255]

    def __init__(self, model, c_pointer=None):
        super().__init__('Position Restraints', model, c_pointer=c_pointer)
        self._prepare_drawings()
        self._model_update_handler = self.model.triggers.add_handler('changes', self._model_changes_cb)
        self._restraint_update_handler = self.triggers.add_handler('changes', self._restraint_changes_cb)

    def delete(self):
        self.model.triggers.remove_handler(self._model_update_handler)
        super().delete()

    def _prepare_drawings(self):
        pd = self._pin_drawing = Drawing('Target pins')
        bd = self._bond_drawing = Drawing('Restraint bonds')
        bd.skip_bounds = True
        self.add_drawing(pd)
        self.add_drawing(bd)
        pd.set_geometry(*self._target_pin_geometry())
        bd.set_geometry(*self._pseudobond_geometry())
        self.set_bond_color(self._DEFAULT_BOND_COLOR)
        self.set_pin_color(self._DEFAULT_PIN_COLOR)
        pd.display = False
        bd.display = False

    def _pseudobond_geometry(self, segments=3):
        from chimerax.core import surface
        return surface.dashed_cylinder_geometry(segments, height=1.0, nc = 6)

    def _target_pin_geometry(self, handle_radius=0.4, pin_radius=0.05, total_height=1.0):
        from .geometry import pin_geometry
        return pin_geometry(handle_radius, pin_radius, total_height)

    def set_pin_color(self, color):
        '''
        Set the colour of the pins used to represent position restraints.

        Args:
            * color:
                - an iterable of four integers (0..255) representing the
                  (r,g,b,a) colour values
        '''
        self._pin_drawing.color = color

    def set_bond_color(self, color):
        '''
        Set the colour of the dashed bond connecting each restrained atom to
        its target pin.

        Args:
            * color:
                - an iterable of four integers (0..255) representing the
                  (r,g,b,a) colour values
        '''
        self._bond_drawing.color = color

    def _restraint_changes_cb(self, trigger_name, changes):
        mgr, changes = changes
        update_bonds = True
        update_targets = False
        change_types = list(changes.keys())
        if 'target changed' in change_types:
            update_targets = True
        if 'enabled/disabled' in change_types:
            update_targets = True
        if 'display changed' in change_types:
            update_targets = True
        self.update_graphics(update_bonds, update_targets)

    def _model_changes_cb(self, trigger_name, changes):
        update_bonds = False
        update_targets = False
        changes = changes[1]
        if changes.num_deleted_atoms() > 0:
            update_bonds = True
            update_targets = True
        atom_reasons = changes.atom_reasons()
        if 'display changed' in atom_reasons or 'hide changed' in atom_reasons:
            update_bonds = True
            update_targets = True
        if 'coord changed' in atom_reasons:
            update_bonds = True
        self.update_graphics(update_bonds, update_targets)

    _show_pd = True
    _show_bd = True

    def update_graphics(self, update_bonds=True, update_targets=True):
        '''
        Update the restraints drawing. Happens automatically every time
        restraints or coordinates change. It should rarely/never be necessary to
        call this manually.
        '''
        if not self.visible:
            return
        pd = self._pin_drawing
        bd = self._bond_drawing
        visibles = self.visible_restraints
        n = len(visibles)
        if n==0:
            pd.display = False
            bd.display = False
            return
        pd.display = self._show_pd
        bd.display = self._show_bd
        if update_bonds:
            self._update_bond_drawing(bd, visibles, n)
        if update_targets:
            self._update_pin_drawing(pd, visibles, n)

    def _update_pin_drawing(self, pd, visibles, n):
        from chimerax.core.geometry import Places
        xyzr = numpy.ones((n,4), numpy.double)
        xyzr[:, :3] = visibles.targets
        pd.positions = Places(shift_and_scale=xyzr)

    def _update_bond_drawing(self, bd, visibles, n):
        from chimerax.core.geometry import Places
        bd.positions = Places(opengl_array = visibles._bond_cylinder_transforms)

    def _get_restraints(self, atoms, create=False):
        '''
        Get restraints for all non-hydrogen atoms in atoms. If create is True,
        any restraints not found will be created. All atoms must be in the
        molecule belonging to this manager!
        '''
        f = c_function('position_restraint_mgr_get_restraint',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_bool, ctypes.c_size_t,
                ctypes.c_void_p),
            ret=ctypes.c_size_t)
        atoms = atoms[atoms.element_names != 'H']
        n = len(atoms)
        ret = numpy.empty(n, cptr)
        num = f(self._c_pointer, atoms._c_pointers, create, n, pointer(ret))
        return _position_restraints(ret[0:num])

    def add_restraints(self, atoms):
        '''
        Returns a :py:class:`Position_Restraints` instance covering the
        non-hydrogen atoms in the input. Hydrogen atoms are unrestrainable by
        design, since their low mass can easily cause instability with strong
        restraints. Only one :cpp:class:`Position_Restraint` will ever be
        created for a given atom, so it is safe to use this repeatedly on the
        same set of atoms. Use this method when you want to ensure that all
        atoms in the input will be restrainable in simulations.

        Args:
            * atoms:
                - a :py:class:`chimerax.Atoms` instance
        '''
        return self._get_restraints(atoms, create=True)

    def add_restraint(self, atom):
        '''
        Singular form of :func:`add_restraints`, returning a
        :py:class:`Position_Restraint` instance (or None if the atom is not
        restrainable).

        Args:
            * atom:
                - a :py:class:`chimerax.Atom` instance
        '''
        from chimerax.core.atomic import Atoms
        r = self._get_restraints(Atoms([atom]), create=True)
        if len(r):
            return r[0]
        return None

    def get_restraints(self, atoms):
        '''
        Returns a :py:class:`Position_Restraints` instance covering the atoms in
        the input. Unlike :func:`add_restraints`, no new
        :cpp:class:`Position_Restraint` will be created. Use this method when
        you've already defined all the restraints you need, and don't want to
        accidentally add more.

        Args:
            * atoms:
                - a :py:class:`chimerax.Atoms` instance
        '''
        return self._get_restraints(atoms)

    def get_restraint(self, atom):
        '''
        Singular form of :func:`get_restraints`, returning a
        :py:class:`Position_Restraint` instance (or None if no restraint exists
        for this atom).

        Args:
            * atom:
                - a :py:class:`chimerax.Atom` instance
        '''
        from chimerax.core.atomic import Atoms
        r = self._get_restraints(Atoms([atom]))
        if len(r):
            return r[0]
        return None


    @property
    def num_restraints(self):
        '''Number of restraints currently managed by this manager.'''
        return c_function('position_restraint_mgr_num_restraints',
            args=(ctypes.c_void_p,), ret=ctypes.c_size_t)(self._c_pointer)

    def __len__(self):
        return self.num_restraints

    @property
    def visible_restraints(self):
        '''
        Returns a :py:class:`Position_Restraints` instance containing only the
        currently visible restraints. A restraint will be visible if it is
        enabled and its atom is visible. Note that it is not possible to hide an
        enabled restraint independent of its atom, other than by hiding the
        entire restraint manager. This is by design: hiding of individual
        restraints will lead to trouble if the user subsequently forgets they
        are there.
        '''
        f = c_function('position_restraint_mgr_visible_restraints',
            args=(ctypes.c_void_p,), ret=ctypes.py_object)
        return _position_restraints(f(self._c_pointer))


class Tuggable_Atoms_Mgr(_Restraint_Mgr):
    '''
    Manages creation, deletion and mapping of :cpp:class:`Position_Restraint`
    objects to act as interactive tugging forces, and drawing of the applied
    force vectors. Appears as a child :py:class:`chimeraX.Model` under the
    :py:class:`chimerax.AtomicStructure` it manages.

    The preferred way to create/retrieve the tuggable atoms manager for a
    given :py:class:`AtomicStructure` instance m is:

    .. code:: python

        from chimerax.isolde import session_extensions as sx
        ta_mgr = sx.get_tuggable_atoms_mgr(m)
    '''
    _DEFAULT_ARROW_COLOR = [100, 255, 100, 255]
    # _NEAR_ATOMS_COLOR = [0,255,255,255]
    def __init__(self, model, c_pointer=None):
        super().__init__('Tuggable atoms', model, c_pointer=c_pointer)
        self._prepare_drawings()
        self._model_update_handler = self.model.triggers.add_handler('changes', self._model_changes_cb)
        self._restraint_update_handler = self.triggers.add_handler('changes', self._restraint_changes_cb)
        # self._show_nearest_atoms = False
        # self._near_atoms_radius = 0.5

    def _prepare_drawings(self):
        ad = self._arrow_drawing = Drawing('Tugging force vectors')
        self.add_drawing(ad)
        ad.set_geometry(*self._arrow_geometry())
        self.set_arrow_color(self._DEFAULT_ARROW_COLOR)

        # nd = self._nearest_atoms_drawing = Drawing('Nearest tuggable atoms')
        # from chimerax.surface.shapes import sphere_geometry2
        # nd.vertices, d.normals, d.triangles = sphere_geometry2(80)
        # self.set_near_atoms_color(self._NEAR_ATOMS_COLOR)

    def set_arrow_color(self, color):
        '''
        Set the colour applied to the arrows used to represent the tugging force
        vector(s).

        Args:
            * color:
                - an (r,g,b,a) iterable of 4 integers in the range (0..255)
        '''
        self._arrow_drawing.color = color

    # def set_near_atoms_color(self, color):
    #     self._nearest_atoms_drawing.color = color

    def _arrow_geometry(self):
        from .geometry import simple_arrow_geometry
        return simple_arrow_geometry(radius=3, centered_on_origin=True)

    def _restraint_changes_cb(self, trigger_name, changes):
        mgr, changes = changes
        self.update_graphics()

    def _model_changes_cb(self, trigger_name, changes):
        update_needed = False
        changes = changes[1]
        if changes.num_deleted_atoms() > 0:
            update_needed = True
        atom_reasons = changes.atom_reasons()
        if 'display changed' in atom_reasons or 'hide changed' in atom_reasons:
            update_needed = True
        if 'coord changed' in atom_reasons:
            update_needed = True
        self.update_graphics()

    def update_graphics(self):
        '''
        Update the force vector drawing. Happens automatically every time
        restraints or coordinates change. It should rarely/never be necessary to
        call this manually.
        '''
        if not self.visible:
            return
        ad = self._arrow_drawing
        visibles = self.visibles
        n = len(visibles)
        if n==0:
            ad.display = False
            return
        ad.display = True
        self._update_arrow_drawing(ad, visibles, n)

    def _update_arrow_drawing(self, ad, visibles, n):
        from chimerax.core.geometry import Places
        ad.positions = Places(opengl_array = visibles._bond_cylinder_transforms)

    def _get_tuggables(self, atoms, create=False):
        f = c_function('tuggable_atoms_mgr_get_restraint',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_bool, ctypes.c_size_t,
                ctypes.c_void_p),
            ret=ctypes.c_size_t)
        atoms = atoms[atoms.element_names != 'H']
        n = len(atoms)
        ret = numpy.empty(n, cptr)
        num = f(self._c_pointer, atoms._c_pointers, create, n, pointer(ret))
        return _tuggable_atoms(ret[0:num])

    def add_tuggables(self, atoms):
        '''
        Returns a :py:class:`Tuggable_Atoms` instance covering the non-hydrogen
        atoms in the input. Hydrogen atoms cannot be tugged by design, since
        their low mass can easily cause instability with strong forces. Only one
        :cpp:class:`Position_Restraint` will ever be created for a given atom,
        so it is safe to use this repeatedly on the same set of atoms. Use this
        method when you want to ensure that all atoms in the input will be
        tuggable in simulations.

        Args:
            * atoms:
                - a :py:class:`chimerax.Atoms` instance
        '''
        return self._get_tuggables(atoms, create=True)

    def add_tuggable(self, atom):
        '''
        Singular form of :func:`add_tuggable`, returning a
        :py:class:`Tuggable_Atom` instance (or None if the atom is not
        tuggable).

        Args:
            * atom:
                - a :py:class:`chimerax.Atom` instance
        '''
        from chimerax.core.atomic import Atoms
        t = self._get_tuggables(Atoms([atom]), create=True)
        if len(t):
            return t[0]
        return None

    def get_tuggables(self, atoms):
        '''
        Returns a :py:class:`Tuggable_Atoms` instance covering the atoms in the
        input. Unlike :func:`add_tuggables`, no new
        :cpp:class:`Position_Restraint` will be created.

        Args:
            * atoms:
                - a :py:class:`chimerax.Atoms` instance
        '''
        return self._get_tuggables(atoms)

    def add_tuggable(self, atom):
        '''
        Singular form of :func:`get_tuggable`, returning a
        :py:class:`Tuggable_Atom` instance (or None if the atom is not
        tuggable).

        Args:
            * atom:
                - a :py:class:`chimerax.Atom` instance
        '''
        from chimerax.core.atomic import Atoms
        t = self._get_tuggables(Atoms([atom]), create=False)
        if len(t):
            return t[0]
        return None


    @property
    def num_tuggables(self):
        '''
        Number of tuggable atoms currently managed by this manager.
        '''
        return c_function('tuggable_atoms_mgr_num_restraints',
            args=(ctypes.c_void_p,), ret=ctypes.c_size_t)(self._c_pointer)

    def __len__(self):
        return self.num_tuggables

    @property
    def visibles(self):
        '''
        Returns a :py:class:`Tuggable_Atoms` instance containing only those
        tuggables with currently-visible force vectors. A vector will be visible
        if and only if the  its associated atom is both visible and currently
        being tugged.
        '''
        f = c_function('tuggable_atoms_mgr_visible_restraints',
            args=(ctypes.c_void_p,), ret=ctypes.py_object)
        return _tuggable_atoms(f(self._c_pointer))

class Distance_Restraint_Mgr(_Restraint_Mgr):
    '''
    Manages distance restraints (Atom pairs with distances and spring constants)
    and their visualisations for a single atomic structure. Appears as a child
    :py:class:`chimerax.Model` under the :py:class:`chimerax.AtomicStructure`
    it manages.

    The preferred way to create/retrieve the distance restraint manager for a
    given :py:class:`AtomicStructure` instance m is:

    .. code:: python

        from chimerax.isolde import session_extensions as sx
        dr_mgr = sx.get_distance_restraint_mgr(m)

    '''
    _DEFAULT_BOND_COLOR = [168, 255, 230, 255]
    _DEFAULT_TARGET_COLOR = [128, 215, 190, 255]
    def __init__(self, model, c_pointer=None):
        '''
        Prepare a distance restraint manager for a given atomic model.
        '''

        session = model.session
        if not hasattr(session, 'isolde_changes') or session.isolde_changes.deleted:
            ct = self._change_tracker = Restraint_Change_Tracker(session)
        else:
            ct = self._change_tracker = session.isolde_changes

        if c_pointer is None:
            f = c_function('distance_restraint_mgr_new',
                args=(ctypes.c_void_p, ctypes.c_void_p,),
                ret=ctypes.c_void_p)
            c_pointer =(f(model._c_pointer, ct._c_pointer))
        super().__init__('Distance Restraints', model, c_pointer)
        self._prepare_drawing()
        self._model_update_handler = self.model.triggers.add_handler('changes', self._model_changes_cb)
        self._restraint_update_handler = self.triggers.add_handler('changes', self._restraint_changes_cb)

    def delete(self):
        self.model.triggers.remove_handler(self._model_update_handler)
        super().delete()

    def _prepare_drawing(self):
        bd = self._bond_drawing = Drawing('Restraint bonds')
        bd.skip_bounds = True
        self.add_drawing(bd)
        bd.set_geometry(*self._pseudobond_geometry())
        self.set_bond_color(self._DEFAULT_BOND_COLOR)
        td = self._target_drawing = Drawing('Target distances')
        td.skip_bounds = True
        td.set_geometry(*self._target_geometry())
        self.add_drawing(td)
        self.set_target_color(self._DEFAULT_TARGET_COLOR)
        bd.display = False

    def set_bond_color(self, color):
        '''
        Set the colour of the cylinder representing the distance restraint.

        Args:
            * color:
                - an iterable of four integers (0..255) representing the
                  (r,g,b,a) colour values
        '''
        self._bond_drawing.color = color

    def set_target_color(self, color):
        '''
        Set the colour of the cylinder representing the target distance.

        Args:
            * color:
                - an iterable of four integers (0..255) representing the
                  (r,g,b,a) colour values
        '''
        self._target_drawing.color = color

    def _target_geometry(self):
        '''
        Length is scaled to the target distance. Radius scales according to
        spring constant.
        '''
        from chimerax.core import surface
        return surface.cylinder_geometry(radius=1.0, height=1.0)

    def _pseudobond_geometry(self):
        '''
        Connects the two restrained atoms. Radius is fixed.
        '''
        from chimerax.core import surface
        return surface.cylinder_geometry(radius = 0.025, height=1.0, caps=False)

    def _restraint_changes_cb(self, trigger_name, changes):
        mgr, changes = changes
        self.update_graphics()

    def _model_changes_cb(self, trigger_name, changes):
        update_needed = False
        changes = changes[1]
        if changes.num_deleted_atoms():
            update_needed = True
        atom_reasons = changes.atom_reasons()
        if 'display changed' in atom_reasons or 'hide changed' in atom_reasons:
            update_needed = True
        if 'coord changed' in atom_reasons:
            update_needed = True
        if update_needed:
            self.update_graphics()

    _show_bd = True
    _show_td = True

    def update_graphics(self):
        '''
        Update the restraints drawing. Happens automatically every time
        restraints or coordinates change. It should rarely/never be necessary to
        call this manually.
        '''
        if not self.visible:
            return
        bd = self._bond_drawing
        td = self._target_drawing
        visibles = self.visible_restraints
        n = len(visibles)
        if n==0:
            bd.display = False
            td.display = False
            return
        bd.display = self._show_bd
        td.display = self._show_td
        self._update_bond_drawing(bd, visibles, n)
        self._update_target_drawing(td, visibles, n)

    def _update_bond_drawing(self, bd, visibles, n):
        bd.positions = visibles._bond_cylinder_transforms

    def _update_target_drawing(self, td, visibles, n):
        td.positions = visibles._target_transforms

    def _get_restraint(self, atom1, atom2, create=False):
        f = c_function('distance_restraint_mgr_get_restraint',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_bool),
            ret = ctypes.c_void_p)
        from chimerax.core.atomic import Atoms
        atoms = Atoms([atom1, atom2])
        return _distance_restraint_or_none(f(self._c_pointer, atoms._c_pointers, create))

    def add_restraint(self, atom1, atom2):
        '''
        Creates and returns a :py:class:`Distance_Restraint` between the given
        atoms if possible, or returns the existing one if it already exists.
        Note that distance restraints can be added to a running simulation, but
        will trigger an expensive (a few seconds) reinitialisation of the
        simulation context.

        Args:
            * atom1, atom2:
                - :py:class:`chimerax.Atom` instances. Both atoms must be in the
                  :py:class:`chimerax.AtomicStructure` belonging to this
                  manager, must not be the same atom or directly bonded to each
                  other, and should not be hydrogens.
        '''
        return self._get_restraint(atom1, atom2, True)

    def get_restraint(self, atom1, atom2):
        '''
        Returns the :py:class:`Distance_Restraint` between the given atoms if it
        exists, otherwise None.

        Args:
            * atom1, atom2:
                - :py:class:`chimerax.Atom` instances. Both atoms must be in the
                  :py:class:`chimerax.AtomicStructure` belonging to this
                  manager.
        '''
        return self._get_restraint(atom1, atom2, False)

    def atom_restraints(self, atom):
        '''
        Returns a :py:class:`Distance_Restraints` instance encompassing all
        distance restraints involving the given atom.

        Args:
            * atom:
                - a :py:class:`chimerax.Atom` instance from the
                  :py:class:`chimerax.AtomicStructure` belonging to this
                  manager.
        '''
        f = c_function('distance_restraint_mgr_atom_restraints',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t),
            ret=ctypes.py_object)
        return _distance_restraints(f(self._c_pointer, atom._c_pointer_ref, 1))

    def atoms_restraints(self, atoms):
        '''
        Returns a :py:class:`Distance_Restraints` instance encompassing all
        distance restraints involving at least one of the atoms. If you want
        only restraints where both atoms are in the set, use
        :func:`intra_restraints` instead.

        Args:
            * atoms:
                - a :py:class:`chimerax.Atoms` instance. All atoms must be part
                  of the :py:class:`chimerax.AtomicStructure` belonging to this
                  manager.
        '''
        f = c_function('distance_restraint_mgr_atom_restraints',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t),
            ret=ctypes.py_object)
        return _distance_restraints(f(self._c_pointer, atoms._c_pointers, len(atoms)))

    def intra_restraints(self, atoms):
        '''
        Returns a :py:class:`Distance_Restraints` instance encompassing all
        distance restraints for which both atoms are in the input set.

        Args:
            * atoms:
                - a :py:class:`chimerax.Atoms` instance. All atoms must be part
                  of the :py:class:`chimerax.AtomicStructure` belonging to this
                  manager.
        '''
        n = len(atoms)
        f = c_function('distance_restraint_mgr_intra_restraints',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t),
            ret = ctypes.py_object)
        return _distance_restraints(f(self._c_pointer, atoms._c_pointers, n))

    def _get_ss_restraints(self, residues, create=False):
        '''
        Returns (optionally creating) distance restraints suitable for restraining
        secondary structure. Result is a tuple of two Distance_Restraints
        objects: O(n)-N(n+4) and CA(n)-CA(n+2).
        '''
        from chimerax.atomic import Residue
        # Reduce to only protein residues
        residues = residues[residues.polymer_types == Residue.PT_AMINO]
        # ensure residues are sorted and unique
        m = self.model
        indices = m.residues.indices(residues)
        if -1 in indices:
            raise TypeError('All residues must be from the model attached to this handler!')
        residues = m.residues[numpy.sort(indices)].unique()
        n = len(residues)
        f = c_function('distance_restraint_mgr_get_ss_restraints',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_bool),
            ret=ctypes.py_object)
        ptrs = f(self._c_pointer, residues._c_pointers, n, create)
        return tuple((_distance_restraints(ptrs[0]), _distance_restraints(ptrs[1])))

    def add_ss_restraints(self, residues):
        '''
        Creates and returns a tuple of two :py:class:`Distance_Restraints`
        instances containing restraints suitable for restraining secondary
        structure. The first, O(n) - N(n+4) defines the classic alpha
        helical H-bonds. The second, CA(n) - CA^(n+2) is useful for
        "stretching out" beta strands as well as cross-bracing helices. Useful
        target differences for restraining specific secondary structures can
        be retrieved from `restraints.constants.ss_restraints`.

        Args:
            * residues:
                - a :py:class:`chimerax.Residues` instance. All residues must
                  be in the model belonging to this manager.
        '''
        return self._get_ss_restraints(residues, create=True)

    def get_ss_restraints(self, residues):
        '''
        Identical to :func:`add_ss_restraints`, but does not add any new
        restraints.
        '''
        return self._get_ss_restraints(residues, create=False)

    @property
    def visible_restraints(self):
        '''
        Returns a :py:class:`Distance_Restraints` instance encompassing all
        currently visible restraints owned by this manager. A restraint will be
        visible if it is enabled and both its atoms are visible.
        '''
        f = c_function('distance_restraint_mgr_visible_restraints',
            args=(ctypes.c_void_p,),
            ret = ctypes.py_object)
        return _distance_restraints(f(self._c_pointer))

    @property
    def all_restraints(self):
        '''
        Returns a :py:class:`Distance_Restraints` instance encompassing all
        restraints owned by this manager.
        '''
        f = c_function('distance_restraint_mgr_all_restraints',
            args=(ctypes.c_void_p,),
            ret = ctypes.py_object)
        return _distance_restraints(f(self._c_pointer))

class Chiral_Restraint_Mgr(_Restraint_Mgr):
    '''
    Manages creation, deletion and mapping of chiral restraints for a single
    atomic structure. Appears as a :py:class:`chimerax.Model` under the
    :py:class:`chimerax.AtomicStructure` it manages.

    Chirality restraints are always on by default, and are
    not drawn. The restraint targets are likewise fixed quantities - if you
    want to switch chirality, you should replace/rename the residue in question
    to your desired isomer.

    The preferred way to create/retrieve the chirality restraint manager for a
    given :py:class:`AtomicStructure` instance m is:

    .. code:: python

        from chimerax.isolde import session_extensions as sx
        chir_mgr = sx.get_chirality_restraint_mgr(m)
    '''
    def __init__(self, model, c_pointer = None):
        super().__init__('Chirality Restraints', model, c_pointer)

    def _get_restraints(self, chirals, create=False):
        n = len(chirals)
        f = c_function('chiral_restraint_mgr_get_restraint',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_bool, ctypes.c_size_t,
                ctypes.c_void_p),
            ret = ctypes.c_size_t)
        ptrs = numpy.empty(n, cptr)
        found = f(self._c_pointer, chirals._c_pointers, create, n, pointer(ptrs))
        return _chiral_restraints(ptrs[0:found])

    def add_restraints(self, chirals):
        '''
        Returns a :py:class:`Chiral_Restraints` instance covering the
        :py:class:`Chiral_Centres` in the input. Only one
        :cpp:class:`Chiral_Restraint` will ever be created for a given chiral
        centre, so it is safe to use this repeatedly on the same set of centres.
        Use this method when you want to ensure that all chiral centres in the
        input will be restrainable in simulations.

        Args:
            * chirals:
                - a :py:class:`Chiral_Centers` instance.
        '''
        return self._get_restraints(chirals, create=True)

    def add_restraint(self, chiral):
        '''
        Singular form of :func:`add_restraints`, returning a
        :py:class:`Chiral_Restraint` instance (or None if things go
        wrong).

        Args:
            * dihedral:
                - a :py:class:`Chiral_Center` instance
        '''
        from .molarray import Chiral_Centers
        result = self._get_restraints(Dihedrals([dihedral]), create=True)
        if len(result):
            return result[0]
        return None

    def get_restraints(self, chirals):
        '''
        Returns a :py:class:`Chiral_Restraints` instance covering the
        :class:`Chiral_Centers` in the input. Unlike :func:`add_restraints`, no
        new :cpp:class:`Chiral_Restraint` will be created. Use this method when
        you've already defined all the restraints you need, and don't want to
        accidentally add more.

        Args:
            * chirals:
                - a :py:class:`Chiral_Centers` instance.
        '''
        return self._get_restraints(chirals)

    def get_restraint(self, chiral):
        '''
        Singular form of :func:`get_restraints`, returning a
        :py:class:`Chiral_Restraint` instance (or None if no restraint
        exists for this chiral centre).

        Args:
            * chiral:
                - a :py:class:`Chiral_Center` instance
        '''
        from .molarray import Chirals
        r = self._get_restraints(Chirals([chiral]), create=False)
        if len(r):
            return r[0]
        return None

    def _get_restraints_by_atoms(self, atoms, create):
        chir_m = get_chiral_mgr(self.session)
        chirals = chir_m.get_chirals(atoms)
        return self._get_restraints(chirals, create=create)

    def add_restraints_by_atoms(self, atoms):
        '''
        Returns a :py:class:`Chiral_Restraints` covering all chiral centers in
        the given atoms.

        Args:
            * atoms:
                - a :py:class:`chimerax.Atoms` instance
        '''
        return self._get_restraints_by_atoms(atoms, create=True)

    def get_restraints_by_atoms(self, atoms):
        '''
        Returns a :py:class:`Chiral_Restraints` covering only those chiral
        centers in  the given atoms for which a restraint has already been
        defined.

        Args:
            * atoms:
                - a :py:class:`chimerax.Atoms` instance
        '''
        return self._get_restraints_by_atoms(atoms, create=False)

    @property
    def num_restraints(self):
        '''Number of chirality restrains currently owned by this manager.'''
        f = c_function('chiral_restraint_mgr_num_restraints',
            args=(ctypes.c_void_p,),
            ret=ctypes.c_size_t)
        return f(self._c_pointer)


class Proper_Dihedral_Restraint_Mgr(_Restraint_Mgr):
    '''
    Manages creation, deletion, mapping and drawing of proper dihedral
    restraints for a single atomic structure. Appears as a child
    :py:class:`chimerax.Model` under the :py:class:`chimerax.AtomicStructure` it
    manages.

    The preferred way to create/retrieve the proper dihedral restraint manager
    for a given :py:class:`AtomicStructure` instance m is:

    .. code:: python

        from chimerax.isolde import session_extensions as sx
        pdr_mgr = sx.get_proper_dihedral_restraint_mgr(m)
    '''
    def __init__(self, model, c_pointer = None):
        super().__init__('Proper Dihedral Restraints', model, c_pointer)
        self.set_default_colors()
        self._update_needed = True
        self._prepare_drawings()
        self._restraint_changes_handler = self.triggers.add_handler('changes', self._restraint_changes_cb)
        self._atom_changes_handler = model.triggers.add_handler('changes', self._model_changes_cb)
        self.update_graphics()

    def delete(self):
        ah = self._atom_changes_handler
        if ah is not None:
            self.model.triggers.remove_handler(ah)
        super().delete()

    def _prepare_drawings(self):
        from . import geometry
        if not hasattr(self, '_ring_drawing') or self._ring_drawing is None:
            ring_d = self._ring_drawing = Drawing('rings')
            ring_d.skip_bounds = True
            self.add_drawing(ring_d)
        if not hasattr(self, '_post_drawing') or self._post_drawing is None:
            post_d = self._post_drawing = Drawing('posts')
            post_d.skip_bounds = True
            self.add_drawing(post_d)
        ring_d = self._ring_drawing
        post_d = self._post_drawing
        ring_d.set_geometry(*geometry.ring_arrow_with_post(0.5, 0.05, 5, 6, 0.25, 0.1, 0.05, 1))
        post_d.set_geometry(*geometry.post_geometry(0.05, 1, caps=True))
        ring_d.display = False
        post_d.display = False

    def _get_restraints(self, dihedrals, create=False):
        n = len(dihedrals)
        f = c_function('proper_dihedral_restraint_mgr_get_restraint',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_bool, ctypes.c_size_t,
                ctypes.c_void_p),
            ret = ctypes.c_size_t)
        ptrs = numpy.empty(n, cptr)
        found = f(self._c_pointer, dihedrals._c_pointers, create, n, pointer(ptrs))
        return _proper_dihedral_restraints(ptrs[0:found])

    def add_restraints(self, dihedrals):
        '''
        Returns a :py:class:`Proper_Dihedral_Restraints` instance covering the
        dihedrals in the input. Only one :cpp:class:`Proper_Dihedral_Restraint`
        will ever be created for a given dihedral, so it is safe to use this
        repeatedly on the same set of dihedrals. Use this method when you want
        to ensure that all dihedrals in the input will be restrainable in
        simulations.

        Args:
            * dihedrals:
                - a :py:class:`Proper_Dihedrals` instance.
        '''
        return self._get_restraints(dihedrals, create=True)

    def add_restraint(self, dihedral):
        '''
        Singular form of :func:`add_restraints`, returning a
        :py:class:`Proper_Dihedral_Restraint` instance (or None if things go
        wrong).

        Args:
            * dihedral:
                - a :py:class:`Proper_Dihedral` instance
        '''
        from .molarray import Dihedrals
        r = self._get_restraints(Dihedrals([dihedral]), create=True)
        if len(r):
            return r[0]
        return None

    def get_restraints(self, dihedrals):
        '''
        Returns a :py:class:`Proper_Dihedral_Restraints` instance covering the
        dihedrals in the input. Unlike :func:`add_restraints`, no new
        :cpp:class:`Proper_Dihedral_Restraint` will be created. Use this method
        when you've already defined all the restraints you need, and don't want
        to accidentally add more.

        Args:
            * dihedrals:
                - a :py:class:`Proper_Dihedrals` instance.
        '''
        return self._get_restraints(dihedrals)

    def get_restraint(self, dihedral):
        '''
        Singular form of :func:`get_restraints`, returning a
        :py:class:`Proper_Dihedral_Restraint` instance (or None if no restraint
        exists for this dihedral).

        Args:
            * dihedral:
                - a :py:class:`Proper_Dihedral` instance
        '''
        from .molarray import Dihedrals
        r = self._get_restraints(Dihedrals([dihedral]), create=False)
        if len(r):
            return r[0]
        return None

    def _get_restraints_by_residues_and_name(self, residues, name, create):
        pdm = get_proper_dihedral_manager(self.session)
        dihedrals = pdm.get_dihedrals(residues, name)
        return self._get_restraints(dihedrals, create=create)

    def get_restraints_by_residues_and_name(self, residues, name):
        '''
        Convenience function wrapping :func:`get_restraints`. Returns a
        :py:class:`Proper_Dihedral_Restraints` instance containing all existing
        restraints for the dihedrals with the given name in residues.

        Args:
            * residues:
                - A :py:class:`chimerax.Residues` instance
            * name:
                - Lowercase name of a known dihedral definition (e.g. 'phi')
        '''
        return self._get_restraints_by_residues_and_name(residues, name, False)

    def add_restraints_by_residues_and_name(self, residues, name):
        '''
        Convenience function wrapping :func:`add_restraints`. Returns a
        :py:class:`Proper_Dihedral_Restraints` instance covering all dihedrals
        with the given name in residues.

        Args:
            * residues:
                - A :py:class:`chimerax.Residues` instance
            * name:
                - Lowercase name of a known dihedral definition (e.g. 'phi')
        '''
        return self._get_restraints_by_residues_and_name(residues, name, True)

    def get_restraint_by_residue_and_name(self, residue, name):
        '''
        Singular form of :func:`get_restraints_by_residues_and_name`. Returns a
        :py:class:`Proper_Dihedral_Restraint` instance for the dihedral with
        the given name in the given residue if it exists, else None.

        Args:
            * residue:
                - A :py:class:`chimerax.Residue` instance
            * name:
                - Lowercase name of a known dihedral definition (e.g. 'phi')
        '''
        from chimerax.core.atomic import Residues
        r = self._get_restraints_by_residues_and_name(Residues([residue]), name, False)
        if len(r):
            return r[0]
        return None

    def add_restraint_by_residue_and_name(self, residue, name):
        '''
        Singular form of :func:`add_restraints_by_residues_and_name`. Returns a
        :py:class:`Proper_Dihedral_Restraint` instance for the dihedral with
        the given name in the given residue if it exists or can be created,
        else None.

        Args:
            * residue:
                - A :py:class:`chimerax.Residue` instance
            * name:
                - Lowercase name of a known dihedral definition (e.g. 'phi')
        '''
        from chimerax.core.atomic import Residues
        r = self._get_restraints_by_residues_and_name(Residues([residue]), name, True)
        if len(r):
            return r[0]
        return None

    def _get_all_restraints_for_residues(self, residues, create):
        pdm = get_proper_dihedral_manager(self.session)
        dihedrals = pdm.get_all_dihedrals(residues)
        return self._get_restraints(dihedrals, create=create)

    def get_all_restraints_for_residues(self, residues):
        '''
        Returns a :py:class:`Proper_Dihedral_Restraints` containing all existing
        dihedral restraints belonging to the given residues.

        Args:
            * residues:
                - a :py:class:`chimerax.Residues` instance
        '''
        return self._get_all_restraints_for_residues(residues, False)

    def add_all_defined_restraints_for_residues(self, residues):
        '''
        Creates and returns all possible dihedral restraints for all named
        proper dihedrals in the given residues as a
        :py:class:`Proper_Dihedral_Restraints` instance.

        Args:
            * residues:
                - a :py:class:`chimerax.Residues` instance
        '''
        return self._get_all_restraints_for_residues(residues, True)

    @property
    def num_restraints(self):
        '''Number of dihedral restrains currently owned by this manager.'''
        f = c_function('proper_dihedral_restraint_mgr_num_restraints',
            args=(ctypes.c_void_p,),
            ret=ctypes.c_size_t)
        return f(self._c_pointer)

    @property
    def visible_restraints(self):
        '''
        Returns a :py:class:`Proper_Dihedral_Restraints` instance containing
        all currently visible restraints. A restraint is visible if it is
        enabled, at least the axial bond of its :class:`Proper_Dihedral` is
        visible, and its :attr:`display` propery is True. Unlike other restraint
        classes the visibility of active restraints *is* controllable, primarily
        because certain restraints (particularly the peptide omega dihedral) are
        best left always on and would become distracting if visible.
        '''
        f = c_function('proper_dihedral_restraint_mgr_visible_restraints',
            args=(ctypes.c_void_p,),
            ret = ctypes.py_object)
        return _proper_dihedral_restraints(f(self._c_pointer))

    def set_default_colors(self):
        '''
        Set the colour scale for validation indicators back to their defaults.
        '''
        from .validation.constants import validation_defaults as val_defaults
        self.set_color_scale(val_defaults.OUTLIER_COLOR, val_defaults.ALLOWED_COLOR,
             val_defaults.MAX_FAVORED_COLOR)

    def set_color_scale(self, max_c, mid_c, min_c):
        '''
        Set a custom colour scale for the validation indicators.

        Args:
            * max_c:
                - The colour applied to the most-favoured rotamers, as an
                  integer (r,g,b,a) array in the range (0..255)
            * mid_c:
                - The colour applied to rotamers at the border between favoured
                  and allowed, as an integer (r,g,b,a) array in the range
                  (0..255)
            * min_c:
                - The colour applied to outlier rotamers, as an integer
                  (r,g,b,a) array in the range (0..255)
        '''
        f = c_function('proper_dihedral_restraint_mgr_set_colors',
            args=(ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint8),
                  ctypes.POINTER(ctypes.c_uint8), ctypes.POINTER(ctypes.c_uint8)))
        maxc = numpy.array(max_c, uint8)
        midc = numpy.array(mid_c, uint8)
        minc = numpy.array(min_c, uint8)
        for arr in (maxc, midc, minc):
            if len(arr) != 4:
                raise TypeError('Each color should be an array of 4 values in the range (0,255)')
        f(self._c_pointer, pointer(maxc), pointer(midc), pointer(minc))

    def _model_changes_cb(self, trigger_name, changes):
        update_needed = False
        changes = changes[1]
        if changes.num_deleted_atoms():
            update_needed = True
        atom_reasons = changes.atom_reasons()
        if 'display changed' in atom_reasons or 'hide changed' in atom_reasons:
            update_needed = True
        if 'coord changed' in atom_reasons:
            update_needed = True
        if update_needed:
            self.update_graphics()

    def _restraint_changes_cb(self, trigger_name, changes):
        mgr, changes = changes
        # For the time being, just update on any trigger
        self.update_graphics()

    def update_graphics(self):
        '''
        Update the restraints drawing. Happens automatically every time
        restraints or coordinates change. It should rarely/never be necessary to
        call this manually.
        '''
        if not self.visible:
            return
        ring_d = self._ring_drawing
        post_d = self._post_drawing
        visibles = self.visible_restraints
        if not len(visibles):
            ring_d.display = False
            post_d.display = False
            return
        ring_d.display = True
        post_d.display = True
        positions = visibles._annotation_transforms()
        colors = visibles.annotation_colors
        ring_d.positions = positions[0]
        post_d.positions = positions[1]
        ring_d.colors = post_d.colors = colors

class Rotamer_Restraint_Mgr(_Restraint_Mgr):
    '''
    Rotamer restraints are a little special in that they primarily exist as
    convenience wrappers around sets of dihedral restraints. As such,
    :py:class:`Rotamer_Restraint_Mgr` does not handle drawing of restraints,
    but it *does* handle drawing of previews when scrolling through alternate
    target rotamer conformations for a given sidechain.
    '''
    def __init__(self, model, c_pointer=None):
        '''Manages rotamer restraints for a single model'''
        session=model.session
        ct = self._change_tracker = _get_restraint_change_tracker(session)
        from . import session_extensions as sx
        pdr_m = self._pdr_m = sx.get_proper_dihedral_restraint_mgr(model)
        rota_m = self._rota_m = sx.get_rotamer_mgr(session)
        if c_pointer is None:
            f = c_function('rotamer_restraint_mgr_new',
                args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p),
                ret=ctypes.c_void_p)
            c_pointer = f(model._c_pointer, ct._c_pointer, pdr_m._c_pointer,
                rota_m._c_pointer)
        set_c_pointer(self, c_pointer)
        f = c_function('set_rotamer_restraint_mgr_py_instance', args=(ctypes.c_void_p, ctypes.py_object))
        f(self._c_pointer, self)
        Model.__init__(self, 'Rotamer Restraints', session)
        self.pickable=False
        self.model = model
        model.add([self])
        self._preview_model = None

    @property
    def num_restraints(self):
        '''
        Number of :cpp:class:`Rotamer_Restraint` objects maintained by this
        manager.
        '''
        f=c_function('rotamer_restraint_mgr_num_restraints',
            args=(ctypes.c_void_p,),
            ret=ctypes.c_size_t)
        return f(self._c_pointer)

    def _get_restraints(self, rotamers, create=False):
        '''
        Get restraints for the given rotamers. If create is True, any restraints
        not found will be created. All rotamers must be in the molecule
        belonging to this manager!
        '''
        f=c_function('rotamer_restraint_mgr_get_restraint',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
                ctypes.c_bool, ctypes.c_void_p),
            ret=ctypes.c_size_t)
        n = len(rotamers)
        ret = numpy.empty(n, cptr)
        num = f(self._c_pointer, rotamers._c_pointers, n, create, pointer(ret))
        return _rotamer_restraints(ret[0:num])

    def add_restraints(self, rotamers_or_residues):
        '''
        Returns a :py:class:`Rotamer_Restraints` instance covering all rotamers
        in the input, creating those that don't yet exist.

        Args:
            * rotamers_or_residues:
                - a :py:class:`Rotamers` or :py:class:`chimerax.Residues` with
                  all elements belonging to the same molecule as this manager
        '''
        from chimerax.core.atomic import Residues
        if isinstance(rotamers_or_residues, Residues):
            rota_mgr = _get_rotamer_manager(session)
            rotamers = rota_mgr.get_rotamers(rotamers_or_residues)
        else:
            rotamers = rotamers_or_residues
        return self._get_restraints(rotamers, True)

    def get_restraints(self, rotamers_or_residues):
        '''
        Returns a :py:class:`Rotamer_Restraints` instance containing all
        previously-created restraints in the input.

        Args:
            * rotamers_or_residues:
                - a :py:class:`Rotamers` or :py:class:`chimerax.Residues` with
                  all elements belonging to the same molecule as this manager
        '''
        from chimerax.core.atomic import Residues
        if isinstance(rotamers_or_residues, Residues):
            rota_mgr = _get_rotamer_manager(session)
            rotamers = rota_mgr.get_rotamers(rotamers_or_residues)
        else:
            rotamers = rotamers_or_residues
        return self._get_restraints(rotamers, False)

    def add_restraint(self, rotamer_or_residue):
        '''
        Returns a :py:class:`Rotamer_Restraint` instance for the given rotamer
        or residue, attempting to create it if it doesn't already exist. If
        creation is impossible (e.g. the residue is non-rotameric or
        incomplete), returns None.

        Args:
            * rotamer_or_residue:
                - a :py:class:`Rotamer` or :py:class:`chimerax.Residue`
                  belonging to the same molecule as this manager
        '''
        from chimerax.core.atomic import Residue, Residues
        if isinstance(rotamer_or_residue, Residue):
            rota_mgr = _get_rotamer_manager(session)
            r = rota_mgr.get_rotamer(r)
        else:
            r = rotamer_or_residue
        if r is None:
            return None
        rr = self._get_restraints(_rotamers([r]), True)
        if not len(rr):
            return None
        return rr[0]

    def get_restraint(self, rotamer_or_residue):
        '''
        Returns a :py:class:`Rotamer_Restraint` instance for the given rotamer
        or residue if the restraint has previously been created, otherwise
        None.

        Args:
            * rotamer_or_residue:
                - a :py:class:`Rotamer` or :py:class:`chimerax.Residue`
                  belonging to the same molecule as this manager
        '''
        from chimerax.core.atomic import Residue, Residues
        if isinstance(rotamer_or_residue, Residue):
            rota_mgr = _get_rotamer_manager(session)
            r = rota_mgr.get_rotamer(r)
        else:
            r = rotamer_or_residue
        if r is None:
            return None
        rr = self._get_restraints(_rotamers([r]), False)
        if not len(rr):
            return None
        return rr[0]

    def _incr_preview(self, rotamer, incr):
        rr = self.add_restraint(rotamer)
        num_targets = rotamer.num_targets
        pm = self._preview_model
        if pm is not None and pm.rotamer == rotamer:
            current_target = pm.target_index
        else:
            current_target = rr.target_index
        if current_target == -1 and incr == -1:
            new_target = num_targets-1
        else:
            new_target = (current_target + incr) % num_targets
        target_def = rotamer.get_target(new_target)
        self._create_preview(rotamer, target_def, new_target)
        return target_def

    def next_preview(self, rotamer):
        '''
        Returns the target definition (a dict providing name, expected
        frequency, chi angles and estimated standard deviations) for the
        conformation next in order of decreasing favourability for this rotamer,
        and draws a preview overlaying the residue. If the current preview is
        the least favoured conformation, wraps back to the most favoured. Only
        one preview may be displayed at a time.

        Args:
            * rotamer:
                - a :py:class:`Rotamer` instance
        '''
        return self._incr_preview(rotamer, 1)

    def prev_preview(self, rotamer):
        '''
        Returns the target definition (a dict providing name, expected
        frequency, chi angles and estimated standard deviations) for the
        conformation next in order of increasing favourability for this rotamer,
        and draws a preview overlaying the residue. If the current preview is
        the most favoured conformation, wraps back to the least favoured. Only
        one preview may be displayed at a time.

        Args:
            * rotamer:
                - a :py:class:`Rotamer` instance
        '''
        return self._incr_preview(rotamer, -1)

    def commit_preview(self, rotamer):
        '''
        Directly set the coordinates of the associated residue to match the
        current preview. NOTE: if a simulation is running this change will be
        ignored unless you immediately run
        :func:`Sim_Handler.push_coords_to_sim`!

        Args:
            * rotamer:
                - a :py:class:`Rotamer` instance
        '''
        pm = self._preview_model
        if pm is None or pm.rotamer != rotamer:
            raise TypeError('Current preview is for a different rotamer!')
        # Clear any existing chi restraints
        rr = self.add_restraint(rotamer)
        rr.enabled=False
        rotamer.residue.atoms.coords = pm.atoms.coords
        self._remove_preview()

    def set_targets(self, rotamer, target_index = None):
        '''
        Set chi dihedral restraints on the rotamer to match the target
        definition. If no preview is currently active, then the index of the
        target definition must be provided. As well as setting the target
        chi angles, the cutoff angle on each restraint will be set to twice the
        estimated standard deviation.

        Args:
            * rotamer:
                - a :py:class:`Rotamer` instance
            * target_index:
                - Mandatory if no preview is current. The index of the target
                  in the list of available target conformations (see
                  :attr:`Rotamer.num_targets` and :func:`Rotamer.get_target`)
        '''
        rr = self.add_restraint(rotamer)
        if target_index is None:
            pm = self._preview_model
            if pm is None or pm.rotamer != rotamer:
                raise TypeError(
                'No target index has been chosen and there is no suitable preview '
               +'to choose it from!')
            else:
               target_index = pm.target_index
        rr.target_index = target_index
        self._remove_preview()



    #TODO: Implement preview as drawing, to do away with the need for a dummy model
    def _create_preview(self, rotamer, target_def, target_index):
        pm = self._preview_model
        if pm is not None and pm.rotamer != rotamer:
            self._remove_preview()
        if self._preview_model is None:
            from chimerax.core.commands.split import molecule_from_atoms
            pm = self._preview_model = molecule_from_atoms(self.model, rotamer.residue.atoms)
            pm.bonds.radii = 0.1
            from chimerax.atomic import Atom
            pm.atoms.draw_modes = Atom.STICK_STYLE
            self.model.add([pm])
            pm.rotamer = rotamer
        pm.target_def = target_def
        pm.target_index = target_index
        target_angles = target_def['Angles']
        current_angles = rotamer.angles
        rot_angles = numpy.degrees(target_angles-current_angles)
        chis = rotamer.chi_dihedrals
        from chimerax.core.geometry import Place, rotation, matrix
        master_atoms = rotamer.residue.atoms
        pm.atoms.coords = master_atoms.coords
        for i in range(rotamer.num_chi_dihedrals):
            ma_indices = master_atoms.indices(rotamer.moving_atoms(i))
            ma = pm.atoms[ma_indices]
            master_chi = chis[i]
            chi = pm.atoms[master_atoms.indices(master_chi.atoms)]
            coords = chi.coords
            axis = coords[2] - coords[1]
            center = matrix.project_to_axis(coords[3], axis, coords[1])
            tf = rotation(axis, rot_angles[i], center)
            ma.coords = tf.moved(ma.coords)

    def _remove_preview(self):
        if self._preview_model is not None and not self._preview_model.deleted:
            self.session.models.remove([self._preview_model])
            self._preview_model = None



class _Dihedral(State):
    '''
    Base class for Proper_Dihedral and Improper_Dihedral classes. Do not
    instantiate directly.
    '''
    def __init__(self, c_pointer):
        set_c_pointer(self, c_pointer)

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def __str__(self):
        return self.name

    def reset_state(self):
        pass

    @property
    def _natoms(self):
        return 4
    #TODO: remove this hack once ChimeraX c_array_property bug is fixed


class Chiral_Center(_Dihedral):
    '''
    A chiral centre is an atom with four chemically distinct substituents, such
    that it cannot be overlaid with its mirror image by rotation and
    translation. Within ISOLDE, chiral centres are treated as improper dihedrals
    (that is, dihedrals for which the axis isn't a real bond). Atom ordering
    within a :class:`Chiral_Center` is defined with the central atom first,
    followed by  the first three substituents in order of decreasing priority by
    standard Cahn-Ingold-Prelog rules. Typical tetrahedral centres will have
    dihedral angles in the vicinity of +/-35 degrees - S positive, R negative.

    :class:`Chiral_Center` instances are not created directly - they are
    generated as  needed by :class:`Chiral_Mgr` based on dictionary definitions.
    '''

    atoms = c_property('chiral_atoms', cptr, '_natoms', astype=_atoms, read_only=True,
        doc = 'Atoms making up this chiral centre. Read only.')

    angle = c_property('chiral_angle', float64, read_only=True,
        doc = 'Angle in radians. Read only.')
    residue = c_property('chiral_residue', cptr, astype=_residue, read_only=True,
        doc = 'Residue this chiral centre belongs to. Read only.')

    expected_angle = c_property('chiral_center_expected_angle', float64, read_only=True,
        doc='The equilibrium angle of the chiral dihedral in its correct isomeric state. Read only.')
    deviation = c_property('chiral_center_deviation', float64, read_only=True,
        doc='The difference between the current dihedral angle and :attr:`expected_angle`. Read only.')
    chiral_atom = c_property('chiral_center_chiral_atom', cptr, astype=_atom_or_none, read_only=True,
        doc='The chiral atom. Read only.')

class Proper_Dihedral(_Dihedral):
    '''
    A Proper_Dihedral is defined as a dihedral in which the four atoms are
    strictly bonded a1-a2-a3-a4.
    '''
    angle = c_property('proper_dihedral_angle', float64, read_only=True,
        doc = 'Angle in radians. Read only.')
    name = c_property('proper_dihedral_name', string, read_only = True,
        doc = 'Name of this dihedral. Read only.')

    atoms = c_property('proper_dihedral_atoms', cptr, '_natoms', astype=_atoms, read_only=True,
        doc = 'Atoms making up this dihedral. Read only.')
    residue = c_property('proper_dihedral_residue', cptr, astype=_residue, read_only=True,
        doc = 'Residue this dihedral belongs to. Read only.')
    axial_bond = c_property('proper_dihedral_axial_bond', cptr, astype=_bond_or_none, read_only=True,
        doc='Bond forming the axis of this dihedral. Read-only')

class Rama(State):
    def __init__(self, c_pointer):
        set_c_pointer(self, c_pointer)

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def __str__(self):
        return "Not implemented"

    def reset_state(self):
        pass

    @property
    def omega_dihedral(self):
        '''
        Returns a :class:`Proper_Dihedral` instance pointing to the omega
        (peptide bond) dihedral for this residue, or None if there is no
        preceding residue.
        '''
        f = c_function('rama_omega',
            args = (ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p),
            ret = ctypes.c_size_t)
        ret = numpy.empty(1, cptr)
        found = f(self._c_pointer_ref, 1, pointer(ret))
        if found:
            return _proper_dihedral_or_none(ret[0])
        return None

    @property
    def phi_dihedral(self):
        '''
        Returns a :class:`Proper_Dihedral` instance pointing to the phi dihedral
        for this residue, or None if there is no preceding residue.
        '''
        f = c_function('rama_phi',
            args = (ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p),
            ret = ctypes.c_size_t)
        ret = numpy.empty(1, cptr)
        found = f(self._c_pointer_ref, 1, pointer(ret))
        if found:
            return _proper_dihedral_or_none(ret[0])
        return None

    @property
    def psi_dihedral(self):
        '''
        Returns a :class:`Proper_Dihedral` instance pointing to the psi dihedral
        for this residue, or None if there is no following residue.
        '''
        f = c_function('rama_psi',
            args = (ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p),
            ret = ctypes.c_size_t)
        ret = numpy.empty(1, cptr)
        found = f(self._c_pointer_ref, 1, pointer(ret))
        if found:
            return _proper_dihedral_or_none(ret[0])
        return None

    residue = c_property('rama_residue', cptr, astype=_residue_or_none, read_only = True,
            doc = 'The :class:`chimerax.Residue` to which this Rama belongs. Read only.')
    ca_atom = c_property('rama_ca_atom', cptr, astype=_atom_or_none, read_only = True,
            doc = 'The alpha carbon :py:class:`chimerax.Atom` of the amino acid residue. Read only.')
    valid = c_property('rama_is_valid', npy_bool, read_only = True,
            doc = 'True if this residue has all three of omega, phi and psi. Read only.')
    visible = c_property('rama_visible', npy_bool, read_only = True,
            doc = 'True if the alpha carbon of this residue is visible. Read only.')
    visible_ignoring_ribbon = c_property('rama_only_hidden_by_ribbon', npy_bool, read_only=True,
            doc = 'True if the only thing hiding the alpha carbon is the ribbon display. Read only.')
    score = c_property('rama_score', float64, read_only = True,
            doc = 'The score of this residue on the MolProbity Ramachandran contours. Read only.')
    phipsi = c_property('rama_phipsi', float64, 2, read_only = True,
            doc = 'The phi and psi angles for this residue in radians. Read only.')
    angles = c_property('rama_omegaphipsi', float64, 3, read_only = True,
            doc = 'The omega, phi and psi angles for this residue in radians. Read only.')
    case = c_property('rama_case', uint8, read_only=True,
            doc = '''A value representing the Ramachandran case for this residue,
                matching the case definitions in :class:`Rama_Mgr.Rama_Case`. Read only.''')

class Rotamer(State):
    def __init__(self, c_pointer):
        set_c_pointer(self, c_pointer)

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def __str__(self):
        return "Not implemented"

    def reset_state(self):
        pass

    @property
    def angles(self):
        '''
        Returns an array giving the current chi angles (chi1, chi2, ...) for
        this rotamer.
        '''
        f = c_function('rotamer_angles', args=(ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)))
        ret = numpy.empty(self.num_chi_dihedrals, numpy.double)
        f(self._c_pointer, pointer(ret))
        return ret


    def get_target(self, index):
        '''
        For each rotamer type, :class:`Rota_Mgr` stores a set of ideal rotamer
        targets, based upon the Ultimate Rotamer Library [Hintze et al]_. These
        are indexed in decreasing order of expected frequency. The number of
        available target definitions is available from :py:attr:`num_targets`.
        The target definition is returned as a dict containing the name,
        expected frequency, ideal angles and their estimated standard deviations.

        .. [Hintze et al] Proteins. 2016 Sep;84(9):1177-89. doi: 10.1002/prot.25039.
        '''
        f = c_function('rotamer_target_def',
            args=(ctypes.c_void_p, ctypes.c_size_t),
            ret=ctypes.py_object)
        return f(self._c_pointer, index)

    def moving_atoms(self, chi_index):
        '''Returns the set of atoms moved by rotating around the given chi dihedral.'''
        f = c_function('rotamer_moving_atoms',
            args=(ctypes.c_void_p, ctypes.c_size_t),
            ret=ctypes.py_object)
        return _atoms(f(self._c_pointer, chi_index))

    @property
    def nearest_target(self):
        '''
        Returns the definition for the nearest ideal rotamer by minimising the
        sum of differences between the current and target angles. Adds an extra
        entry ['Z scores'] to the returned dict giving the Z score for
        abs(angle-target) for each chi dihedral.
        '''
        f = c_function('rotamer_nearest_target',
            args=(ctypes.c_void_p,),
            ret=ctypes.py_object)
        target_index, zscores = f(self._c_pointer)
        tdict = self.get_target(target_index)
        tdict['Z scores'] = zscores
        return tdict

    @property
    def chi_dihedrals(self):
        '''
        Returns a :class:`Proper_Dihedrals` giving the chi dihedrals for this
        rotamer, in order (chi1, chi2, ...)
        '''
        f = c_function('rotamer_chi_dihedrals',
            args=(ctypes.c_void_p, ctypes.c_void_p))
        ret = numpy.empty(self.num_chi_dihedrals, cptr)
        f(self._c_pointer, pointer(ret))
        return _proper_dihedrals(ret)


    residue = c_property('rotamer_residue', cptr, astype=_residue_or_none, read_only=True,
                doc=':py:class:`chimerax.Residue` this rotamer belongs to. Read only.')
    score = c_property('rotamer_score', float64, read_only=True,
                doc='P-value for the current conformation of this rotamer. Read only.')
    ca_cb_bond = c_property('rotamer_ca_cb_bond', cptr, astype=_bond_or_none, read_only=True,
                doc='The "stem"  :py:class:`chimerax.Bond` of this rotamer. Read only.')
    num_chi_dihedrals = c_property('rotamer_num_chi', uint8, read_only=True,
                doc='Number of dihedrals defining this rotamer. Read only.')
    num_targets = c_property('rotamer_num_target_defs', uint32, read_only=True,
                doc='Number of available target conformations. Read only.')
    visible = c_property('rotamer_visible', npy_bool, read_only=True,
                doc='True if the CA-CB bond of the rotamer is visible. Read only.')

class Position_Restraint(State):
    '''
    Defines a user-adjustable and -switchable harmonic spring connecting an
    atom to a fixed point in space.
    '''
    def __init__(self, c_pointer):
        set_c_pointer(self, c_pointer)

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def __str__(self):
        return "Not implemented"

    def reset_state(self):
        pass

    def clear_sim_index(self):
        f = c_function('position_restraint_clear_sim_index',
            args = (ctypes.c_void_p, ctypes.c_size_t))
        f(self._c_pointer_ref, 1)

    atom = c_property('position_restraint_atom', cptr, astype=_atom_or_none, read_only=True,
        doc = 'Returns the restrained :py:class:`chimerax.Atom`. Read-only.')
    target = c_property('position_restraint_target', float64, 3,
        doc = 'Target (x,y,z) position in Angstroms. Can be written.')
    target_vector = c_property('position_restraint_target_vector', float64, 3, read_only=True,
        doc = 'Returns the vector ("bond") connecting the atom to its target. Read only.')
    spring_constant = c_property('position_restraint_k', float64,
        doc = 'Restraint spring constant in :math:`kJ mol^{-1} nm^{-2}`. Can be written')
    enabled = c_property('position_restraint_enabled', npy_bool,
        doc = 'Enable/disable this position restraint.')
    visible = c_property('position_restraint_visible', npy_bool, read_only=True,
        doc = 'Check whether this restraint is currently visible. Read only.')
    sim_index = c_property('position_restraint_sim_index', int32,
        doc = '''
        Index of this restraint in the associated force object in a running
        simulation. Returns -1 if the restraint is not part of any simulation.
        Can be set, but only if you know what you are doing.
        ''')

class Tuggable_Atom(Position_Restraint):
    _c_class_name = 'position_restraint'

class MDFF_Atom(State):
    '''
    A thin wrapper around a :class:`chimerax.Atom` to handle its coupling
    to a molecular dynamics flexible fitting potential derived from a volumetric
    map.
    '''
    def __init__(self, c_pointer):
        set_c_pointer(self, c_pointer)

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def __str__(self):
        return "Not implemented"

    def reset_state(self):
        pass

    def clear_sim_index(self):
        f = c_function('mdff_atom_clear_sim_index',
            args = (ctypes.c_void_p, ctypes.c_size_t))
        f(self._c_pointer_ref, 1)

    enabled = c_property('mdff_atom_enabled', npy_bool,
        doc='Enable/disable MDFF tugging on this atom or get its current state.')
    atom = c_property('mdff_atom_atom', cptr, astype=_atom_or_none, read_only=True,
        doc='Returns the ChimeraX Atom. Read only.')
    coupling_constant = c_property('mdff_atom_coupling_constant', float64,
        doc='''
        Per-atom MDFF coupling constant. This is multiplied by the global
        coupling constant for the map when calculating the MDFF potential.
        Can be set.
        ''')
    sim_index = c_property('mdff_atom_sim_index', int32,
        doc='''
        Index of this atom in the relevant MDFF Force in a running simulation.
        Returns -1 if the atom is not currently in a simulation. Can be set, but
        only if you know what you are doing.
        ''')

class Distance_Restraint(State):
    '''
    Defines a harmonic spring restraining the distance between a pair of atoms.
    '''
    def __init__(self, c_pointer):
        set_c_pointer(self, c_pointer)

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def __str__(self):
        return "Not implemented"

    def reset_state(self):
        pass

    def clear_sim_index(self):
        f = c_function('distance_restraint_clear_sim_index',
            args = (ctypes.c_void_p, ctypes.c_size_t))
        f(self._c_pointer_ref, 1)

    enabled =c_property('distance_restraint_enabled', npy_bool,
            doc = 'Enable/disable this restraint or get its current state.')
    visible = c_property('distance_restraint_visible', npy_bool, read_only = True,
            doc = 'Restraint will be visible if it is enabled and both atoms are visible.')
    atoms = c_property('distance_restraint_atoms', cptr, 2, astype=_atom_pair, read_only=True,
            doc = 'Returns a tuple of :py:class:`chimerax.Atoms` pointing to the pair of restrained atoms. Read only.' )
    target = c_property('distance_restraint_target', float64,
            doc = 'Target distance in Angstroms')
    spring_constant = c_property('distance_restraint_k', float64,
            doc = 'Restraint spring constant in :math:`kJ mol^{-1} nm^{-2}`')
    distance = c_property('distance_restraint_distance', float64, read_only=True,
            doc = 'Current distance between restrained atoms in Angstroms. Read only.')
    sim_index = c_property('distance_restraint_sim_index', int32,
        doc='''
        Index of this restraint in the relevant MDFF Force in a running
        simulation. Returns -1 if the restraint is not currently in a
        simulation. Can be set, but only if you know what you are doing.
        ''')

class Chiral_Restraint(State):
    '''
    Handles the restraint of a single chiral centre in simulations. Unlike other
    restraints, :class:`Chiral_Restraint` instances are enabled by default. In
    addition, :attr:`target` is immutable, and set to the equilibrium angle of
    the chiral improper dihedral (defined in `dictionaries/chirals.json` and
    viewable via :attr:`Chiral_Mgr.chiral_center_dict`).
    '''
    def __init__(self, c_pointer):
        set_c_pointer(self, c_pointer)

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def __str__(self):
        return "Not implemented"

    def reset_state(self):
        pass

    def clear_sim_index(self):
        f = c_function('chiral_restraint_clear_sim_index',
            args = (ctypes.c_void_p, ctypes.c_size_t))
        f(self._c_pointer_ref, 1)

    @property
    def atoms(self):
        return self.dihedral.atoms

    @property
    def chiral_atom(self):
        return self.dihedral.chiral_atom

    target = c_property('chiral_restraint_target', float64, read_only = True,
        doc = 'Target angle for this restraint in radians. Read only.')
    dihedral = c_property('chiral_restraint_chiral_center', cptr, astype=_chiral_center_or_none, read_only=True,
        doc = 'The restrained :py:class:`Chiral_Center`. Read only.')
    offset = c_property('chiral_restraint_offset', float64, read_only = True,
        doc = 'Difference between current and target angle in radians. Read only.')
    cutoff = c_property('chiral_restraint_cutoff', float64,
        doc = 'Cutoff angle offset below which no restraint will be applied. Can be set.')
    enabled = c_property('chiral_restraint_enabled', npy_bool,
        doc = 'Enable/disable this restraint or get its current state.')
    spring_constant = c_property('chiral_restraint_k', float64,
        doc = 'Get/set the spring constant for this restraint in :math:`kJ mol^{-1} rad^{-2}`')
    sim_index = c_property('chiral_restraint_sim_index', int32,
        doc='''
        Index of this restraint in the relevant Force in a running simulation.
        Returns -1 if the restraint is not currently in a simulation. Can be
        set, but only if you know what you are doing.
        ''')


class Proper_Dihedral_Restraint(State):
    def __init__(self, c_pointer):
        set_c_pointer(self, c_pointer)

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def __str__(self):
        return "Not implemented"

    def reset_state(self):
        pass

    def clear_sim_index(self):
        f = c_function('proper_dihedral_restraint_clear_sim_index',
            args = (ctypes.c_void_p, ctypes.c_size_t))
        f(self._c_pointer_ref, 1)

    @property
    def atoms(self):
        return self.dihedral.atoms

    target = c_property('proper_dihedral_restraint_target', float64,
        doc = 'Target angle for this restraint in radians. Can be written.')
    dihedral = c_property('proper_dihedral_restraint_dihedral', cptr, astype=_proper_dihedral_or_none, read_only=True,
        doc = 'The restrained :py:class:`Proper_Dihedral`. Read only.')
    offset = c_property('proper_dihedral_restraint_offset', float64, read_only = True,
        doc = 'Difference between current and target angle in radians. Read only.')
    cutoff = c_property('proper_dihedral_restraint_cutoff', float64,
        doc = 'Cutoff angle offset below which no restraint will be applied. Can be set.')
    enabled = c_property('proper_dihedral_restraint_enabled', npy_bool,
        doc = 'Enable/disable this restraint or get its current state.')
    display = c_property('proper_dihedral_restraint_display', npy_bool,
        doc = 'Set whether you want this restraint to be displayed when active.')
    visible = c_property('proper_dihedral_restraint_visible', npy_bool, read_only=True,
        doc = 'Is this restraint currently visible? Read-only.')
    spring_constant = c_property('proper_dihedral_restraint_k', float64,
        doc = 'Get/set the spring constant for this restraint in :math:`kJ mol^{-1} rad^{-2}`')
    annotation_color = c_property('proper_dihedral_restraint_annotation_color', uint8, 4, read_only=True,
        doc = 'Get the color of the annotation for this restraint according to the current colormap. Read only.')
    sim_index = c_property('proper_dihedral_restraint_sim_index', int32,
        doc='''
        Index of this restraint in the relevant Force in a running simulation.
        Returns -1 if the restraint is not currently in a simulation. Can be
        set, but only if you know what you are doing.
        ''')

class Rotamer_Restraint(State):
    '''
    Rotamer restraints are special in that they only exist as wrappers around
    combinations of :cpp:class:`Proper_Dihedral_Restraint` objects. As such,
    they have no independent presence in simulations. The chi dihedral
    restraints may still be addressed and controlled independently if
    necessary.
    '''
    def __init__(self, c_pointer):
        set_c_pointer(self, c_pointer)

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def __str__(self):
        return "Not implemented"

    def reset_state(self):
        pass

    @property
    def target(self):
        if self.target_index == -1:
            return None
        return self.rotamer.get_target(self.target_index)

    def set_spring_constant(self, k):
        '''
        Sets the spring constants for all chi dihedrals in the rotamer to k.
        Write-only. To retrieve the current spring constants, use
        :attr:`chi_restraints.spring_constants`
        '''
        f = c_function('set_rotamer_restraint_spring_constant',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_double))
        f(self._c_pointer_ref, 1, k)

    @property
    def chi_restraints(self):
        '''
        Returns a :py:class:`Proper_Dihedral_Restraints` covering the chi
        dihedrals in this rotamer.
        '''
        f = c_function('rotamer_restraint_chi_restraints',
            args=(ctypes.c_void_p, ctypes.c_void_p))
        ret = numpy.empty(self.rotamer.num_chi_dihedrals, cptr)
        f(self._c_pointer, pointer(ret))
        return _proper_dihedral_restraints(ret)

    rotamer = c_property('rotamer_restraint_rotamer', cptr, astype=_rotamer_or_none, read_only=True,
        doc = ':py:class:`Rotamer` to be restrained. Read only.')
    residue = c_property('rotamer_restraint_residue', cptr, astype=_residue_or_none, read_only=True,
        doc = ':py:class:`chimerax.Residue` to be restrained. Read only.')
    enabled = c_property('rotamer_restraint_enabled', npy_bool,
        doc = 'Enable/disable chi dihedral restraints. Returns False if any chi restraint is disabled.')
    target_index = c_property('rotamer_restraint_target_index', int32,
        doc = '''
        Get/set the index of the rotamer target definition giving target angles
        and cutoffs. If no restraint is currently applied, returns the last
        restraint that was applied to this rotamer. If no restraint has ever
        been applied, returns -1. Setting the index automatically apples the
        target and cutoff angles for each chi dihedral restraint, according to
        the chosen target. Cutoffs angles are set to twice the estimated
        standard deviations for the target.
        ''')




# tell the C++ layer about class objects whose Python objects can be instantiated directly
# from C++ with just a pointer, and put functions in those classes for getting the instance
# from the pointer (needed by Collections)
for class_obj in (Chiral_Center, Proper_Dihedral, Rama, Rotamer,
        Position_Restraint, Tuggable_Atom, MDFF_Atom, Distance_Restraint,
        Chiral_Restraint, Proper_Dihedral_Restraint, Rotamer_Restraint):
    if hasattr(class_obj, '_c_class_name'):
        cname = class_obj._c_class_name
    else:
        cname = class_obj.__name__.lower()
    func_name = 'set_' + cname + '_pyclass'
    f = c_function(func_name, args = (ctypes.py_object,))
    f(class_obj)

    func_name = cname + '_py_inst'
    class_obj.c_ptr_to_py_inst = lambda ptr, fname=func_name: c_function(fname,
        args = (ctypes.c_void_p,), ret = ctypes.py_object)(ctypes.c_void_p(int(ptr)))
    func_name = cname + '_existing_py_inst'
    class_obj.c_ptr_to_existing_py_inst = lambda ptr, fname=func_name: c_function(fname,
        args = (ctypes.c_void_p,), ret = ctypes.py_object)(ctypes.c_void_p(int(ptr)))

for class_obj in (Proper_Dihedral_Mgr, Chiral_Mgr, Rama_Mgr, Rota_Mgr,
            Position_Restraint_Mgr, Tuggable_Atoms_Mgr, MDFF_Mgr,
            Distance_Restraint_Mgr, Chiral_Restraint_Mgr,
            Proper_Dihedral_Restraint_Mgr, Rotamer_Restraint_Mgr):
    cname = class_obj.__name__.lower()
    func_name = cname + '_py_inst'
    class_obj.c_ptr_to_py_inst = lambda ptr, fname=func_name: c_function(fname,
        args = (ctypes.c_void_p,), ret = ctypes.py_object)(ctypes.c_void_p(int(ptr)))
    func_name = cname + '_existing_py_inst'
    class_obj.c_ptr_to_existing_py_inst = lambda ptr, fname=func_name: c_function(fname,
        args = (ctypes.c_void_p,), ret = ctypes.py_object)(ctypes.c_void_p(int(ptr)))
