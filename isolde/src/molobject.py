
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

from numpy import uint8, int32, uint32, float64, float32, byte, bool as npy_bool

libdir = os.path.dirname(os.path.abspath(__file__))
libfile = glob.glob(os.path.join(libdir, 'molc.cpython*'))[0]
DATA_DIR = os.path.join(libdir, 'molprobity_data')



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
def _proper_dihedral_or_none(p):
    return Proper_Dihedral.c_ptr_to_py_inst(p) if p else None
def _bond_or_none(p):
    return Bond.c_ptr_to_py_inst(p) if p else None

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
    
    def add_dihedrals(self, dihedrals):
        f = c_function('proper_dihedral_mgr_add_dihedral', 
            args = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t),
            )
        f(self.cpp_pointer, dihedrals._c_pointers, len(dihedrals))
    
    def _load_dict(self):
        import json
        with open(os.path.join(libdir, 'dictionaries', 'named_dihedrals.json'), 'r') as f:
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
                        ctypes.c_size_t, ctypes.c_void_p, ctypes.c_bool), 
                        ret=ctypes.c_size_t)
        n = len(residues)
        key = ctypes.py_object()
        key.value = name
        ptrs  = numpy.empty(n, cptr)
        num_found = f(self._c_pointer, residues._c_pointers, ctypes.byref(key), n, pointer(ptrs), create)
        return _proper_dihedrals(ptrs[0:num_found])
    
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
    
    RAMA_CASE_DETAILS = {
        Rama_Case.NONE: {
            'name': 'Not applicable',
            'file_prefix': None,
            'cutoffs': None
        },
        Rama_Case.CISPRO: {
            'name': 'Cis-proline residues',
            'file_prefix': os.path.join(DATA_DIR, 'rama8000-cispro'),
            'cutoffs': [0.002, 1.0, 0.02]
        },
        Rama_Case.TRANSPRO: {
            'name': 'Trans-proline residues',
            'file_prefix': os.path.join(DATA_DIR, 'rama8000-transpro'),
            'cutoffs': [0.001, 1.0, 0.02]
        },
        Rama_Case.GLYCINE: {
            'name': 'Glycine residues',
            'file_prefix': os.path.join(DATA_DIR, 'rama8000-gly-sym'),
            'cutoffs': [0.001, 1.0, 0.02]
        },
        Rama_Case.PREPRO: {
            'name': 'Residues preceding proline',
            'file_prefix': os.path.join(DATA_DIR, 'rama8000-prepro-noGP'),
            'cutoffs': [0.001, 1.0, 0.02]
        },
        Rama_Case.ILEVAL: {
            'name': 'Isoleucine or valine residues',
            'file_prefix': os.path.join(DATA_DIR, 'rama8000-ileval-nopreP'),
            'cutoffs': [0.001, 1.0, 0.02]
        },
        Rama_Case.GENERAL: {
            'name': 'General amino acid residues',
            'file_prefix': os.path.join(DATA_DIR, 'rama8000-general-noGPIVpreP'),
            'cutoffs': [0.0005, 1.0, 0.02]
        }
    }
    
    def _prepare_all_validators(self):
        from .validation import generate_interpolator_data
        for case, details in self.RAMA_CASE_DETAILS.items():
            file_prefix = details['file_prefix']
            if file_prefix is not None:
                i_data = generate_interpolator_data(file_prefix, True)
                self.add_interpolator(case, *i_data)
    
    def __init__(self, session, c_pointer=None):
        if hasattr(session, 'rama_mgr'):
            raise RuntimeError('Session already has a Ramachandran manager!')
        if not hasattr(session, 'proper_dihedral_mgr'):
            raise RuntimeError('Proper_Dihedral_Mgr must be initialised first!')
        cname = type(self).__name__.lower()
        if c_pointer is None:
            new_func = cname + '_new'
            c_pointer = c_function(new_func, ret=ctypes.c_void_p)()
        set_c_pointer(self, c_pointer)
        f = c_function('set_'+cname+'_py_instance', args=(ctypes.c_void_p, ctypes.py_object))
        f(self._c_pointer, self)
        self.session = session
        self._dihedral_mgr = session.proper_dihedral_mgr     
        self._prepare_all_validators()   
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
    
    @property
    def dihedral_manager(self):
        return self._dihedral_mgr
    
    def add_interpolator(self, rama_case, ndim, axis_lengths, min_vals, max_vals, data):
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




    def valid_rama_residues(self, residues):
        f = c_function('proper_dihedral_mgr_valid_rama_residues', args=(
                        ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t,
                        ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p), 
                        ret=ctypes.c_size_t)
        n = len(residues)
        r_ptrs = numpy.empty(n, cptr)
        o_ptrs = numpy.empty(n, cptr)
        phi_ptrs = numpy.empty(n, cptr)
        psi_ptrs = numpy.empty(n, cptr)
        num_found = f(self.dihedral_manager.cpp_pointer, residues._c_pointers, n, 
            pointer(r_ptrs), pointer(o_ptrs), pointer(phi_ptrs), pointer(psi_ptrs))
        r = _residues(r_ptrs[0:num_found])
        omega = _proper_dihedrals(o_ptrs[0:num_found])
        phi = _proper_dihedrals(phi_ptrs[0:num_found])
        psi = _proper_dihedrals(psi_ptrs[0:num_found])
        return (r, omega, phi, psi)

    def rama_cases(self, omegas, psis):
        f = c_function('rama_mgr_rama_cases', 
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_size_t, ctypes.POINTER(ctypes.c_uint8)))
        n = len(omegas)
        if len(psis) != n:
            raise TypeError('Array lengths must be equal!')
        ret = numpy.empty(n, uint8);
        f(self._c_pointer, omegas._c_pointers,
            psis._c_pointers, n, pointer(ret))
        return ret
    
    def validate(self, residues, omegas, phis, psis, cases = None):
        '''
        Returns Ramachandran scores for all 
        '''
        if cases is None:
            cases = self.rama_cases(omegas, psis)
        n = len(residues)
        for arr in (omegas, phis, psis):
            if len(arr) != n:
                raise TypeError('Array lengths must be equal!')
        return self._validate(residues, omegas, phis, psis, cases, n)

    def _validate(self, residues, omegas, phis, psis, cases, n):
        f = c_function('rama_mgr_validate',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, 
                  ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(ctypes.c_uint8),
                  ctypes.c_size_t, ctypes.POINTER(ctypes.c_double)))
        ret = numpy.empty(n, float64)
        f(self._c_pointer, residues._c_pointers, omegas._c_pointers, 
            phis._c_pointers, psis._c_pointers, pointer(cases), n, pointer(ret))
        return ret
        
        
    


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
    
    residue = c_property('proper_dihedral_residue', cptr, astype=_residue, read_only=True, doc = 'Residue this dihedral belongs to. Read only.')
    target = c_property('proper_dihedral_target', float32,
        doc='Target angle in radians. Will be automatically wrapped to (-pi,pi)')
    spring_constant = c_property('proper_dihedral_spring_constant', float32,
        doc='Spring constant for dihedral restraint in kJ/mol/radian**2.')
    axial_bond = c_property('proper_dihedral_axial_bond', cptr, astype=_bond_or_none, read_only=True,
        doc='Bond forming the axis of this dihedral. Read-only')


# tell the C++ layer about class objects whose Python objects can be instantiated directly
# from C++ with just a pointer, and put functions in those classes for getting the instance
# from the pointer (needed by Collections)
for class_obj in [Proper_Dihedral, ]:
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
    
Proper_Dihedral_Mgr.c_ptr_to_py_inst = lambda ptr: c_function('proper_dihedral_mgr_py_inst',
    args = (ctypes.c_void_p,), ret = ctypes.py_object)(ctypes.c_void_p(int(ptr)))
Proper_Dihedral_Mgr.c_ptr_to_existing_py_inst = lambda ptr: c_function("proper_dihedral_mgr_existing_py_inst",
    args = (ctypes.c_void_p,), ret = ctypes.py_object)(ctypes.c_void_p(int(ptr)))

Rama_Mgr.c_ptr_to_py_inst = lambda ptr: c_function('rama_mgr_py_inst',
    args = (ctypes.c_void_p,), ret = ctypes.py_object)(ctypes.c_void_p(int(ptr)))
Rama_Mgr.c_ptr_to_existing_py_inst = lambda ptr: c_function("rama_mgr_existing_py_inst",
    args = (ctypes.c_void_p,), ret = ctypes.py_object)(ctypes.c_void_p(int(ptr)))

