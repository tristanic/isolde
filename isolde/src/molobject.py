
import os, sys, glob
import numpy
import ctypes
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

from chimerax.core.models import Model

from numpy import int8, uint8, int32, uint32, float64, float32, byte, bool as npy_bool

import json

libdir = os.path.dirname(os.path.abspath(__file__))
libfile = glob.glob(os.path.join(libdir, 'molc.cpython*'))[0]
DATA_DIR = os.path.join(libdir, 'molprobity_data')
DICT_DIR = os.path.join(libdir, 'dictionaries')

_c_functions = CFunctions(os.path.splitext(libfile)[0])
c_property = _c_functions.c_property
cvec_property = _c_functions.cvec_property
c_function = _c_functions.c_function
c_array_function = _c_functions.c_array_function

def _asptr(arr, dtype):
    return arr.ctypes.data_as(pointer(dtype))

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
def _pseudobond_or_none(p):
    from chimerax.core.atomic import Pseudobond
    return Pseudobond.c_ptr_to_py_inst(p) if p else None
def _proper_dihedral_restraint_or_none(p):
    return Proper_Dihedral_Restraint.c_ptr_to_py_inst(p) if p else None
def _proper_dihedral_restraints(p):
    from .molarray import Proper_Dihedral_Restraints
    return Proper_Dihedral_Restraints(p)

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

    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

class Proper_Dihedral_Mgr(_Dihedral_Mgr):

    def __init__(self, session, c_pointer=None):
        super().__init__(session, c_pointer=c_pointer)
        self._load_dict()
        if hasattr(session, 'proper_dihedral_mgr') and not session.proper_dihedral_mgr.deleted:
            raise RuntimeError('Session already has a proper dihedral manager!')
        session.proper_dihedral_mgr = self

    def __delete__(self):
        self.delete()
        super().__delete__()

    def delete(self):
        c_function('proper_dihedral_mgr_delete', args=(ctypes.c_void_p,))(self.cpp_pointer)
        delattr(self.session, 'proper_dihedral_mgr')

    def delete_dihedrals(self, dihedrals):
        '''
        Delete all dihedrals in a :class:`Proper_Dihedrals`. Note that
        this will not affect the constituent atoms in any way, and
        should not actually be necessary in most cases. Dihedrals are
        automatically deleted at the C++ level when their manager or
        any of their constituent atoms are deleted.
        '''
        f = c_function('proper_dihedral_mgr_delete_dihedral',
                args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p))
        f(self.cpp_pointer, len(dihedrals), dihedrals._c_pointers)

    # def add_dihedrals(self, dihedrals):
    #     f = c_function('proper_dihedral_mgr_add_dihedral',
    #         args = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t),
    #         )
    #     f(self.cpp_pointer, dihedrals._c_pointers, len(dihedrals))

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
                externals = numpy.zeros(4, numpy.bool)
                if type(d_data) == list and type(d_data[1]) == list:
                    externals = numpy.array(d_data[1]).astype(numpy.bool)
                    d_data = d_data[0]
                elif type(d_data) != list:
                    continue
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

    def find_dihedrals(self, model):
        dihedral_dict = self._dihedral_dict
        amino_acid_resnames = dihedral_dict['aminoacids']
        r = model.residues
        self._reserve(len(r))
        aa_residues = r[numpy.in1d(r.names, amino_acid_resnames)]
        f = c_function('proper_dihedral_mgr_new_dihedral',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p))
        self._find_peptide_backbone_dihedrals(dihedral_dict, aa_residues, f)
        self._find_rotameric_dihedrals(dihedral_dict, amino_acid_resnames, aa_residues, f)

    def _find_peptide_backbone_dihedrals(self, dihedral_dict, aa_residues, f):
        for key in dihedral_dict['all_protein'].keys():
            k = ctypes.py_object()
            k.value = key
            f(self._c_pointer, aa_residues._c_pointers, len(aa_residues),
                ctypes.byref(k))

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
                f(self._c_pointer, aa_residues._c_pointers, len(aa_residues),
                    ctypes.byref(k))

    def get_dihedral(self, residue, name, create=True):
        from chimerax.core.atomic import Residues
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
            residues:
                A :class:`Residues` object
            name:
                Name of the desired dihedral
            create (default = True):
                If a dihedral is not found, try to create it.
        '''
        f = c_function('proper_dihedral_mgr_get_dihedrals', args=(
                        ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                        ctypes.c_size_t, ctypes.c_bool),
                        ret=ctypes.py_object)
        n = len(residues)
        key = ctypes.py_object()
        key.value = name
        return _proper_dihedrals(f(self._c_pointer, residues._c_pointers, ctypes.byref(key), n, create))

    @property
    def num_mapped_dihedrals(self):
        f = c_function('proper_dihedral_mgr_num_mapped_dihedrals', args=(ctypes.c_void_p,), ret=ctypes.c_size_t)
        return f(self._c_pointer)


    def __len__(self):
        return self.num_mapped_dihedrals

class Rama_Mgr:
    '''
    Manager for Ramachandran scoring of protein residues. Should have only
    one per session.
    '''
    from enum import IntEnum
    class Rama_Case(IntEnum):
        '''
        Enumerators for the different Ramachandran cases. These match
        an enumerator in the C++ level, so don't change them unless you
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
        FAVORED=0
        ALLOWED=1
        OUTLIER=2
        NA=-1

    from .constants import validation_defaults as val_defaults

    RAMA_CASE_DETAILS = {
        Rama_Case.NONE: {
            'name': 'Not applicable',
            'file_prefix': None,
            'cutoffs': None
        },
        Rama_Case.CISPRO: {
            'name': 'Cis-proline residues',
            'file_prefix': os.path.join(DATA_DIR, 'rama8000-cispro'),
            'cutoffs': [val_defaults.CISPRO_ALLOWED, val_defaults.CISPRO_OUTLIER]
        },
        Rama_Case.TRANSPRO: {
            'name': 'Trans-proline residues',
            'file_prefix': os.path.join(DATA_DIR, 'rama8000-transpro'),
            'cutoffs': [val_defaults.TRANSPRO_ALLOWED, val_defaults.TRANSPRO_OUTLIER]
        },
        Rama_Case.GLYCINE: {
            'name': 'Glycine residues',
            'file_prefix': os.path.join(DATA_DIR, 'rama8000-gly-sym'),
            'cutoffs': [val_defaults.GLYCINE_ALLOWED, val_defaults.GLYCINE_OUTLIER]
        },
        Rama_Case.PREPRO: {
            'name': 'Residues preceding proline',
            'file_prefix': os.path.join(DATA_DIR, 'rama8000-prepro-noGP'),
            'cutoffs': [val_defaults.PREPRO_ALLOWED, val_defaults.PREPRO_OUTLIER]
        },
        Rama_Case.ILEVAL: {
            'name': 'Isoleucine or valine residues',
            'file_prefix': os.path.join(DATA_DIR, 'rama8000-ileval-nopreP'),
            'cutoffs': [val_defaults.ILEVAL_ALLOWED, val_defaults.ILEVAL_OUTLIER]
        },
        Rama_Case.GENERAL: {
            'name': 'General amino acid residues',
            'file_prefix': os.path.join(DATA_DIR, 'rama8000-general-noGPIVpreP'),
            'cutoffs': [val_defaults.GENERAL_ALLOWED, val_defaults.GENERAL_OUTLIER]
        }
    }

    def _prepare_all_validators(self):
        from .validation import generate_interpolator_data
        for case, details in self.RAMA_CASE_DETAILS.items():
            file_prefix = details['file_prefix']
            if file_prefix is not None:
                i_data = generate_interpolator_data(file_prefix, True)
                self._add_interpolator(case, *i_data)

    def __init__(self, session, c_pointer=None):
        if hasattr(session, 'rama_mgr'):
            raise RuntimeError('Session already has a Ramachandran manager!')
        if not hasattr(session, 'proper_dihedral_mgr'):
            raise RuntimeError('Proper_Dihedral_Mgr must be initialised first!')
        self.session = session
        dmgr = self._dihedral_mgr = session.proper_dihedral_mgr
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
        Reset the Ramachandran cutoffs to default values
        '''
        dd = self.RAMA_CASE_DETAILS
        for case, cd in dd.items():
            cutoffs = cd['cutoffs']
            if cutoffs is not None:
                self._set_cutoffs(case, *cutoffs)

    def _set_cutoffs(self, case, allowed, outlier):
        f = c_function('rama_mgr_set_cutoffs',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_double, ctypes.c_double))
        f(self._c_pointer, case, allowed, outlier)

    @property
    def cutoffs(self):
        f = c_function('rama_mgr_get_cutoffs',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_double)))
        cdict = {}
        for case in self.Rama_Case:
            if case != Rama_Case.NONE:
                cutoffs = numpy.empty(2, numpy.double)
                f(self._c_pointer, case, pointer(cutoffs))
                cdict[case] = cutoffs
        return cdict

    def set_default_colors(self):
        from .constants import validation_defaults as val_defaults
        self.set_color_scale(val_defaults.MAX_FAVORED_COLOR, val_defaults.ALLOWED_COLOR,
            val_defaults.OUTLIER_COLOR, val_defaults.NA_COLOR)

    def set_color_scale(self, max_c, mid_c, min_c, na_c):
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
        return self._dihedral_mgr

    def _add_interpolator(self, rama_case, ndim, axis_lengths, min_vals, max_vals, data):
        '''
        Create a RegularGridInterpolator for the given Ramachandran
        contours, and add it to the manager.
        Args:
        rama_case:
            An integer corresponding to the Ramachandran case
            (see Rama_Mgr.Rama_Case for valid values)
        axis_lengths:
            A numpy int array giving the number of points along each axis
        min_vals:
            A numpy double array giving the minimum value for each axis
        max_vals:
            A numpy double array giving the maximum value for each axis
        data:
            A 2D numpy array containing the gridded data
        '''
        f = c_function('rama_mgr_add_interpolator',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.c_size_t,
                ctypes.POINTER(ctypes.c_uint32), ctypes.POINTER(ctypes.c_double),
                ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)))
        axis_lengths = axis_lengths.astype(uint32)
        f(self._c_pointer, rama_case, ndim, pointer(axis_lengths),
            pointer(min_vals), pointer(max_vals), pointer(data))

    # def valid_rama_residues(self, residues):
    #     f = c_function('proper_dihedral_mgr_valid_rama_residues', args=(
    #                     ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
    #                     ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p),
    #                     ret=ctypes.c_size_t)
    #     n = len(residues)
    #     r_ptrs = numpy.empty(n, cptr)
    #     o_ptrs = numpy.empty(n, cptr)
    #     phi_ptrs = numpy.empty(n, cptr)
    #     psi_ptrs = numpy.empty(n, cptr)
    #     num_found = f(self.dihedral_manager.cpp_pointer, residues._c_pointers, n,
    #         pointer(r_ptrs), pointer(o_ptrs), pointer(phi_ptrs), pointer(psi_ptrs))
    #     r = _residues(r_ptrs[0:num_found])
    #     omega = _proper_dihedrals(o_ptrs[0:num_found])
    #     phi = _proper_dihedrals(phi_ptrs[0:num_found])
    #     psi = _proper_dihedrals(psi_ptrs[0:num_found])
    #     return (r, omega, phi, psi)

    def rama_cases(self, residues):
        f = c_function('rama_mgr_rama_case',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_uint8)))
        n = len(residues)
        ret = numpy.empty(n, uint8);
        f(self._c_pointer, residues._c_pointers, n, pointer(ret))
        return ret

    def bin_scores(scores, cases):
        n = len(scores)
        if len(cases) != len(scores):
            raise TypeError('Both arrays must be the same length!')
        f = c_function('rama_mgr_bin_scores',
            args=(ctypes.c_void_p, ctypes.POINTER(ctypes.c_double),
                ctypes.POINTER(ctypes.c_uint8), ctypes.c_size_t,
                ctypes.POINTER(ctypes.c_int8)))
        ret = numpy.empty(n, int8)
        f(self._c_pointer, pointer(scores), pointer(cases), n, ret)
        return ret

    def outliers(self, residues):
        scores, cases = self.validate_by_residue(residues)
        bins = self.bin_scores(scores, cases)
        return residues[bins==self.Rama_Bin.OUTLIER]

    def validate(self, residues_or_ramas):
        '''
        Returns Ramachandran scores for a set of pre-defined valid
        Ramachandran cases. The input to this function is typically
        the output of :class:`Rama_Mgr`.rama_cases(). For a slower
        but more robust method which can handle invalid (non-protein
        and/or terminal) residues, use validate_by_residue().
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
        Returns a nx4 uint8 :class:`Numpy` array giving a color for each
        residue corresponding to the current colormap.
        '''
        f = c_function('rama_mgr_validate_and_color',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
                  ctypes.POINTER(ctypes.c_uint8)))
        n = len(ramas)
        colors = numpy.empty((n,4), uint8)
        f(self._c_pointer, ramas._c_pointers, n, pointer(colors))
        return colors

    def color_cas_by_rama_score(self, ramas):
        f = c_function('rama_mgr_validate_and_color_cas',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,))
        n = len(ramas)
        f(self._c_pointer, ramas._c_pointers, n)

    def bin_scores(self, scores, cases):
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

    #######
    # Access to the underlying interpolator data
    #######
    def interpolator_dim(self, rama_case):
        f = c_function('rama_mgr_interpolator_dim',
            args=(ctypes.c_void_p, ctypes.c_size_t),
            ret=ctypes.c_size_t)
        return f(self._c_pointer, rama_case)

    def interpolator_axis_lengths(self, rama_case):
        dim = self.interpolator_dim(rama_case)
        f = c_function('rama_mgr_interpolator_axis_lengths',
            args=(ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_uint32)))
        ret = numpy.empty(dim, uint32)
        f(self._c_pointer, rama_case, pointer(ret))
        return ret

    def interpolator_limits(self, rama_case):
        dim = self.interpolator_dim(rama_case)
        f = c_function('rama_mgr_interpolator_minmax',
            args=(ctypes.c_void_p, ctypes.c_size_t,
                ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)))
        minvals = numpy.empty(dim, float64)
        maxvals = numpy.empty(dim, float64)
        f(self._c_pointer, rama_case, pointer(minvals), pointer(maxvals))
        return (minvals, maxvals)

    def interpolator_values(self, rama_case):
        shape = self.interpolator_axis_lengths(rama_case)
        f = c_function('rama_mgr_interpolator_values',
            args=(ctypes.c_void_p, ctypes.c_size_t,
                ctypes.POINTER(ctypes.c_double)))
        data = numpy.empty(shape, float64)
        f(self._c_pointer, rama_case, pointer(data))
        return data

    def interpolator_axes(self, rama_case):
        lengths = self.interpolator_axis_lengths(rama_case)
        minmax = self.interpolator_limits(rama_case)
        axes = [numpy.linspace(minmax[0][i], minmax[1][i], lengths[i]) for i in range(len(lengths))]
        return tuple(axes)

class Rota_Mgr:
    from enum import IntEnum
    class Rota_Bin(IntEnum):
        FAVORED=0
        ALLOWED=1
        OUTLIER=2
        NA=-1

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

    def __init__(self, session, c_pointer=None):
        if hasattr(session, 'rota_mgr'):
            raise RuntimeError('Session already has a Rotamer manager!')
        if not hasattr(session, 'proper_dihedral_mgr'):
            raise RuntimeError('Proper_Dihedral_Mgr must be initialised first!')
        cname = type(self).__name__.lower()
        self._dihedral_mgr = session.proper_dihedral_mgr
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

    def _load_target_defs(self):
        '''
        Load rotamer targets from their JSON file and store them in order of
        decreasing prevalence for each residue.
        '''
        with open(os.path.join(DICT_DIR, 'rota_data.json'), 'rt') as f:
            rd = json.load(f)
        from collections import OrderedDict
        ordered_rd = self._rota_targets = dict()
        for aa, data in rd.items():
            rdata = [(name, d) for (name, d) in data.items()]
            rdata = sorted(rdata, key=lambda d: d[1]['freq'], reverse=True)
            r_dict = OrderedDict()
            for r, angles in rdata:
                for k in ('angles', 'esds'):
                    angles[k] = numpy.array(angles[k])
                r_dict[r]=angles
            ordered_rd[aa] = r_dict

    def get_rota_targets(self, resname):
        '''
        Returns an OrderedDict giving rotamer name, angles, esds and frequencies
        sorted in order of decreasing frequency.
        '''
        from copy import copy
        rd = self._rota_targets
        if resname in rd.keys():
            return copy(rd[resname])
        return None

    def nearest_valid_rotamer(self, residue_or_rotamer):
        '''
        Returns the name of the nearest valid rotamer given this residue's
        current conformation, its percentage frequency in the Top8000 dataset,
        and the z_score for the dihedral with the biggest deviation from ideal.
        '''
        from chimerax.core.atomic import Residue
        if type(residue_or_rotamer) == Residue:
            r = self.get_rotamer(residue)
        else:
            r = residue_or_rotamer
        if r is None:
            return None
        angles = numpy.degrees(r.angles)
        targets = self.get_rota_targets(residue.name)
        target_angles = numpy.array([t['angles'] for t in targets.values()])
        differences = angles-target_angles
        differences[differences<-180] += 360
        differences[differences>=180] -= 360
        differences = numpy.abs(differences)
        sums = differences.sum(axis=1)
        min_index = numpy.argmin(sums)
        min_differences = differences[min_index]
        key = list(targets.keys())[min_index]
        rot = targets[key]
        esds = rot['esds']
        z_score = numpy.max(min_differences/esds)
        freq = rot['freq']
        return (key, freq, z_score)

    def set_default_colors(self):
        from .constants import validation_defaults as val_defaults
        self.set_color_scale(val_defaults.MAX_FAVORED_COLOR, val_defaults.ALLOWED_COLOR,
            val_defaults.OUTLIER_COLOR)

    def set_color_scale(self, max_c, mid_c, min_c):
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

    def get_color_scale(self):
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
        Reset the rotamer cutoffs to default values
        '''
        from .constants import validation_defaults as vc
        self._set_cutoffs(vc.ROTA_ALLOWED_CUTOFF, vc.ROTA_OUTLIER_CUTOFF)

    def _set_cutoffs(self, allowed, outlier):
        f = c_function('rota_mgr_set_cutoffs',
            args=(ctypes.c_void_p, ctypes.c_double, ctypes.c_double))
        f(self._c_pointer, allowed, outlier)

    @property
    def cutoffs(self):
        f = c_function('rota_mgr_get_cutoffs',
            args=(ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)))
        cutoffs = numpy.empty(2, numpy.double)
        f(self._c_pointer, pointer(cutoffs))
        return cutoffs

    @property
    def defined_rotamers(self):
        if not hasattr(self, '_defined_rotamer_dict') or self._defined_rotamer_dict is None:
            self._load_rotamer_defs()
        return self._defined_rotamer_dict

    def _load_defined_rotamers(self):
        with open(os.path.join(DICT_DIR, 'rota_data.json'), 'r') as f:
            self._defined_rotamer_dict = json.load(f)

    def _load_rotamer_defs(self):
        f = c_function('rota_mgr_add_rotamer_def',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_bool))
        dd = self._dihedral_mgr.dihedral_dict
        pd = dd['residues']['protein']
        aa_names = dd['aminoacids']
        for aa in aa_names:
            nchi = pd[aa]['nchi']
            if nchi > 0:
                symm = pd[aa]['symm']
                key = ctypes.py_object()
                key.value = aa
                f(self._c_pointer, ctypes.byref(key), nchi, symm)

    def _prepare_all_validators(self):
        from .validation import generate_interpolator_data
        dmgr = self._dihedral_mgr
        prefix = os.path.join(DATA_DIR, 'rota8000-')
        for aa in dmgr.dihedral_dict['aminoacids']:
            fname = prefix + aa.lower()
            if not os.path.isfile(fname+'.data') and not os.path.isfile(fname+'.pickle'):
                # Not a rotameric residue
                continue
            idata = generate_interpolator_data(fname, True)
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
        from chimerax.core.atomic import Residues
        rots = self.get_rotamers(Residues([residue]))
        if len(rots):
            return rots[0]
        return None

    def get_rotamers(self, residues):
        f = c_function('rota_mgr_get_rotamer',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t),
            ret=ctypes.py_object)
        n = len(residues)
        return _rotamers(f(self._c_pointer, residues._c_pointers, n))

    def validate_rotamers(self, rotamers):
        f = c_function('rota_mgr_validate_rotamer',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_double)))
        n = len(rotamers)
        ret = numpy.empty(n, numpy.double)
        f(self._c_pointer, rotamers._c_pointers, n, pointer(ret))
        return ret

    def validate_rotamers_threaded(self, rotamers):
        if self._thread_running():
            raise RuntimeError('Thread already running!')
        init_f = c_function('rota_mgr_validate_rotamer_threaded',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_double)))
        n = len(rotamers)
        ret = self._thread_result = numpy.empty(n, numpy.double)
        from .delayed_reaction import delayed_reaction
        delayed_reaction(self.session, init_f,
            (self._c_pointer, rotamers._c_pointers, n, pointer(ret),),
            self._thread_done, self._get_thread_result, ())

    def _thread_running(self):
        f = c_function('rota_mgr_thread_running',
            args=(ctypes.c_void_p,), ret=ctypes.c_bool)
        return f(self._c_pointer)

    def _thread_done(self):
        f = c_function('rota_mgr_thread_done',
            args=(ctypes.c_void_p,), ret=ctypes.c_bool)
        return f(self._c_pointer)

    def _get_thread_result(self):
        f = c_function('rota_mgr_finalize_thread',
            args=(ctypes.c_void_p,))
        f(self._c_pointer)
        print('thread done!')

    def validate_residues(self, residues):
        f = c_function('rota_mgr_validate_residue',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_double)))
        n = len(residues)
        ret = numpy.empty(n, numpy.double)
        f(self._c_pointer, residues._c_pointers, n, pointer(ret))
        return ret

    def non_favored_rotamers(self, rotamers):
        f = c_function('rota_mgr_non_favored',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)),
            ret=ctypes.c_size_t)
        n = len(rotamers)
        ptrs = numpy.empty(n, cptr)
        scores = numpy.empty(n, float64)
        found = f(self._c_pointer, rotamers._c_pointers, n, pointer(ptrs), pointer(scores))
        print ("Found {} bad rotamers".format(found)) #DELETEME
        return (_rotamers(ptrs[0:found]), scores[0:found])

    def validate_scale_and_color_rotamers(self, rotamers, max_scale = 2.0, non_favored_only = True):
        f = c_function('rota_mgr_validate_scale_and_color_rotamers',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
                ctypes.c_double, ctypes.c_bool, ctypes.c_void_p,
                ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_uint8)),
                ret=ctypes.c_size_t)
        n = len(rotamers)
        rot_out = numpy.empty(n, cptr)
        scale_out = numpy.empty(n, float64)
        color_out = numpy.empty((n,4), uint8)
        count = f(self._c_pointer, rotamers._c_pointers, n, max_scale,
            non_favored_only, pointer(rot_out), pointer(scale_out), pointer(color_out))
        return (_rotamers(rot_out[0:count]), scale_out[0:count], color_out[0:count])


class Restraint_Change_Tracker:
    '''
    A per-session singleton tracking changes in ISOLDE restraints, and firing
    triggers as necessary.
    '''
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


class _Restraint_Mgr(Model):
    '''Base class. Do not instantiate directly.'''
    def __init__(self, name, model, c_pointer=None):
        session = model.session
        if not hasattr(self, '_change_tracker'):
            if not hasattr(session, 'isolde_changes') or session.isolde_changes.deleted:
                ct = self._change_tracker = Restraint_Change_Tracker(session)
            else:
                ct = self._change_tracker = session.isolde_changes

        cname = type(self).__name__.lower()
        if c_pointer is None:
            new_func = cname + '_new'
            c_pointer = c_function(
                new_func, args=(ctypes.c_void_p, ctypes.c_void_p),
                          ret=ctypes.c_void_p)(model._c_pointer, ct._c_pointer)
        set_c_pointer(self, c_pointer)
        f = c_function('set_'+cname+'_py_instance', args=(ctypes.c_void_p, ctypes.py_object))
        f(self._c_pointer, self)
        super().__init__(name+' Manager', session)
        self.model = model
        model.add([self])

    def delete(self):
        cname = type(self).__name__.lower()
        del_func = cname+'_delete'
        c_function(del_func, args=(ctypes.c_void_p,))(self._c_pointer)


    @property
    def cpp_pointer(self):
        '''Value that can be passed to C++ layer to be used as pointer (Python int)'''
        return self._c_pointer.value

    @property
    def deleted(self):
        '''Has the C++ side been deleted?'''
        return not hasattr(self, '_c_pointer')

class Position_Restraint_Mgr_Base(_Restraint_Mgr):
    '''
    Manages position restraints for a single atomic structure.
    '''
    def __init__(self, model, c_pointer=None):
        super().__init__('Position Restraints', model, c_pointer)

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

    def get_restraints(self, atoms):
        return self._get_restraints(atoms)

    def add_restraints(self, atoms):
        return self._get_restraints(atoms, create=True)

    @property
    def num_restraints(self):
        return c_function('position_restraint_mgr_num_restraints',
            args=(ctypes.c_void_p,), ret=ctypes.c_size_t)(self._c_pointer)

    def __len__(self):
        return self.num_restraints

    @property
    def visible_restraints(self):
        f = c_function('position_restraint_mgr_visible_restraints',
            args=(ctypes.c_void_p,), ret=ctypes.py_object)
        return _position_restraints(f(self._c_pointer))


class Distance_Restraint_Mgr_Base(_Restraint_Mgr):
    '''
    Manages distance restraints (Atom pairs with distances and spring constants)
    and their visualisations for a single atomic structure.
    '''
    def __init__(self, model, c_pointer=None):
        '''
        Prepare a distance restraint manager for a given atomic model.
        '''

        # TODO: Defer pseudobondgroup creation to after initialisation to avoid
        # duplicating code in _Restrain_Mgr
        session = model.session
        if not hasattr(session, 'isolde_changes') or session.isolde_changes.deleted:
            ct = self._change_tracker = Restraint_Change_Tracker(session)
        else:
            ct = self._change_tracker = session.isolde_changes

        pbg = self._pbgroup = model.pseudobond_group('ISOLDE distance restraints')
        if c_pointer is None:
            f = c_function('distance_restraint_mgr_new',
                args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,),
                ret=ctypes.c_void_p)
            c_pointer =(f(model._c_pointer, ct._c_pointer, pbg._c_pointer))
        super().__init__('Distance Restraints', model, c_pointer)
        self.add([pbg])

    def add_restraint(self, atom1, atom2):
        f = c_function('distance_restraint_mgr_get_restraint',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_bool),
            ret = ctypes.c_void_p)
        from chimerax.core.atomic import Atoms
        atoms = Atoms([atom1, atom2])
        return _distance_restraint_or_none(f(self._c_pointer, atoms._c_pointers, True))

    def intra_restraints(self, atoms):
        n = len(atoms)
        f = c_function('distance_restraint_mgr_intra_restraints',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t),
            ret = ctypes.py_object)
        return _distance_restraints(f(self._c_pointer, atoms._c_pointers, n))

    # def delete(self):
    #     self.session.models.remove([self._pbgroup])
    #     self._pbgroup.delete()
    #     super().delete(self)


class Proper_Dihedral_Restraint_Mgr(_Restraint_Mgr):
    def __init__(self, model, c_pointer = None):
        super().__init__('Proper Dihedral Restraints', model, c_pointer)
        self.set_default_colors()

    def _get_restraints(self, dihedrals, create=False):
        n = len(dihedrals)
        f = c_function('proper_dihedral_restraint_mgr_get_restraint',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_bool, ctypes.c_size_t,
                ctypes.c_void_p),
            ret = ctypes.c_size_t)
        ptrs = numpy.empty(n, cptr)
        found = f(self._c_pointer, dihedrals._c_pointers, create, n, pointer(ptrs))
        return _proper_dihedral_restraints(ptrs[0:n])

    def add_restraints(self, dihedrals):
        return self._get_restraints(dihedrals, create=True)

    def set_default_colors(self):
        from .constants import validation_defaults as val_defaults
        self.set_color_scale(val_defaults.OUTLIER_COLOR, val_defaults.ALLOWED_COLOR,
             val_defaults.MAX_FAVORED_COLOR)

    def set_color_scale(self, max_c, mid_c, min_c):
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

    angle = c_property('dihedral_angle', float32, read_only=True, doc = 'Angle in radians. Read only.')
    name = c_property('dihedral_name', string, read_only = True, doc = 'Name of this dihedral. Read only.')

class Proper_Dihedral(_Dihedral):
    '''
    A Proper_Dihedral is defined as a dihedral in which the four atoms are
    strictly bonded a1-a2-a3-a4.
    '''
    residue = c_property('proper_dihedral_residue', cptr, astype=_residue, read_only=True, doc = 'Residue this dihedral belongs to. Read only.')
    # target = c_property('proper_dihedral_target', float32,
    #     doc='Target angle in radians. Will be automatically wrapped to (-pi,pi)')
    # spring_constant = c_property('proper_dihedral_spring_constant', float32,
    #     doc='Spring constant for dihedral restraint in kJ/mol/radian**2.')
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

    residue = c_property('rama_residue', cptr, astype=_residue_or_none, read_only = True,
            doc = 'The residue to which this Rama belongs. Read only.')
    ca_atom = c_property('rama_ca_atom', cptr, astype=_atom_or_none, read_only = True,
            doc = 'The alpha carbon of the amino acid residue. Read only.')
    valid = c_property('rama_is_valid', npy_bool, read_only = True,
            doc = 'True if this residue has all three of omega, phi and psi. Read only.')
    score = c_property('rama_score', float64, read_only = True,
            doc = 'The score of this residue on the MolProbity Ramachandran contours. Read only.')
    phipsi = c_property('rama_phipsi', float64, 2, read_only = True,
            doc = 'The phi and psi angles for this residue in radians. Read only.')
    angles = c_property('rama_omegaphipsi', float64, 3, read_only = True,
            doc = 'The omega, phi and psi angles for this residue in radians. Read only.')


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
        f = c_function('rotamer_angles', args=(ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)))
        ret = numpy.empty(self.num_chi_dihedrals, numpy.double)
        f(self._c_pointer, pointer(ret))
        return ret

    residue = c_property('rotamer_residue', cptr, astype=_residue_or_none, read_only=True,
                doc='Residue this rotamer belongs to. Read only.')
    score = c_property('rotamer_score', float32, read_only=True,
                doc='P-value for the current conformation of this rotamer. Read only.')
    ca_cb_bond = c_property('rotamer_ca_cb_bond', cptr, astype=_bond_or_none, read_only=True,
                doc='The "stem" bond of this rotamer. Read only.')
    num_chi_dihedrals = c_property('rotamer_num_chi', uint8, read_only=True,
                doc='Number of dihedrals defining this rotamer')

class Position_Restraint(State):
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
    def _bond_cylinder_transform(self):
        '''Transform mapping a unit cylinder onto the restraint bond. Read only.'''
        f = c_function('position_restraint_bond_transform',
            args = (ctypes.c_void_p, ctypes.c_size_t,
                ctypes.POINTER(ctypes.c_float)))
        n = 1
        transform = empty((n,4,4), float32)
        f(self._c_pointer, n, pointer(transform))
        return transform

    atom = c_property('position_restraint_atom', cptr, astype=_atom_or_none, read_only=True,
        doc = 'Returns the restrained atom. Read-only.')
    target = c_property('position_restraint_target', float64, 3,
        doc = 'Target (x,y,z) position in Angstroms. Can be written.')
    target_vector = c_property('position_restraint_target_vector', float64, 3, read_only=True,
        doc = 'Returns the vector ("bond") connecting the atom to its target. Read only.')
    spring_constant = c_property('position_restraint_k', float64,
        doc = 'Restraint spring constant in kJ mol-1 Angstrom-2. Can be written')
    enabled = c_property('position_restraint_enabled', bool,
        doc = 'Enable/disable this position restraint.')
    visible = c_property('position_restraint_visible', bool, read_only=True,
        doc = 'Check whether this restraint is currently visible. Read only.')


class Distance_Restraint(State):
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

    enabled =c_property('distance_restraint_enabled', npy_bool,
            doc = 'Enable/disable this restraint or get its current state.')
    atoms = c_property('distance_restraint_atoms', cptr, 2, astype=_atom_pair, read_only=True,
            doc = 'Returns the pair of restrained atoms. Read only.' )
    target = c_property('distance_restraint_target', float64,
            doc = 'Target distance in Angstroms')
    spring_constant = c_property('distance_restraint_k', float64,
            doc = 'Restraint spring constant in kJ mol-1 Angstrom-2')
    distance = c_property('distance_restraint_distance', float64, read_only=True,
            doc = 'Current distance between restrained atoms in Angstroms. Read only.')
    pseudobond = c_property('distance_restraint_pbond', cptr, astype=_pseudobond_or_none, read_only=True,
            doc = 'Pseudobond visualisation of the restraint. Read only.')

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



    target = c_property('proper_dihedral_restraint_target', float64,
        doc = 'Target angle for this restraint in radians. Can be written.')
    cutoff = c_property('proper_dihedral_restraint_cutoff', float64,
        doc = 'Cutoff angle below which no restraint will be applied. Can be set.')
    enabled = c_property('proper_dihedral_restraint_enabled', npy_bool,
        doc = 'Enable/disable this restraint or get its current state.')
    display = c_property('proper_dihedral_restraint_display', npy_bool,
        doc = 'Set whether you want this restraint to be displayed when active.')
    visible = c_property('proper_dihedral_restraint_visible', npy_bool, read_only=True,
        doc = 'Is this restraint currently visible? Read-only.')
    spring_constant = c_property('proper_dihedral_restraint_k', float64,
        doc = 'Get/set the spring constant for this restraint in kJ mol-1 rad-2')
    annotation_color = c_property('proper_dihedral_restraint_annotation_color', uint8, 4, read_only=True,
        doc = 'Get the color of the annotation for this restraint according to the current colormap. Read only.')


# tell the C++ layer about class objects whose Python objects can be instantiated directly
# from C++ with just a pointer, and put functions in those classes for getting the instance
# from the pointer (needed by Collections)
for class_obj in (Proper_Dihedral, Rama, Rotamer, Position_Restraint,
                  Distance_Restraint, Proper_Dihedral_Restraint,):
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

for class_obj in (Proper_Dihedral_Mgr, Rama_Mgr, Rota_Mgr, Position_Restraint_Mgr_Base,
            Distance_Restraint_Mgr_Base, Proper_Dihedral_Restraint_Mgr):
    func_name = cname + '_py_inst'
    class_obj.c_ptr_to_py_inst = lambda ptr, fname=func_name: c_function(fname,
        args = (ctypes.c_void_p,), ret = ctypes.py_object)(ctypes.c_void_p(int(ptr)))
    func_name = cname + '_existing_py_inst'
    class_obj.c_ptr_to_existing_py_inst = lambda ptr, fname=func_name: c_function(fname,
        args = (ctypes.c_void_p,), ret = ctypes.py_object)(ctypes.c_void_p(int(ptr)))
