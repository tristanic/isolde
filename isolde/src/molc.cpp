/*
 * Based upon ChimeraX molc.cpp. Interface between Python and C++ 
 * implementations of objects handling ChimeraX atomic data.
 */
 

#include <Python.h>     // Use PyUnicode_FromString

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>      // use PyArray_*(), NPY_*

#include <atomstruct/Atom.h>
#include <atomstruct/AtomicStructure.h>
#include <atomstruct/Bond.h>
#include <atomstruct/Chain.h>
#include <atomstruct/ChangeTracker.h>
#include <atomstruct/CoordSet.h>
#include <atomstruct/connect.h>
#include <atomstruct/destruct.h>     // Use DestructionObserver
#include <atomstruct/MolResId.h>
#include <atomstruct/PBGroup.h>
#include <atomstruct/polymer.h>
#include <atomstruct/Pseudobond.h>
#include <atomstruct/PBGroup.h>
#include <atomstruct/Residue.h>
#include <atomstruct/RibbonXSection.h>
#include <atomstruct/Ring.h>
#include <atomstruct/seq_assoc.h>
#include <atomstruct/Sequence.h>
#include <arrays/pythonarray.h>           // Use python_voidp_array()
#include <pysupport/convert.h>     // Use cset_of_chars_to_pyset

#include "atomic_cpp/dihedral.h"
#include "atomic_cpp/dihedral_mgr.h"
#include "geometry/geometry.h"
#include "interpolation/nd_interp.h"

#include <functional>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <vector>
#include <cmath>

#include "molc.h"
using namespace atomstruct;
using namespace isolde;

// --------------------------------------------------------------------
// dihedral functions
//

extern "C" EXPORT void set_proper_dihedral_pyclass(PyObject* py_class)
{
    try {
        Dihedral::set_py_class(py_class);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT PyObject* proper_dihedral_py_inst(void* ptr)
{
    Dihedral *d = static_cast<Dihedral*>(ptr);
    try {
        return d->py_instance(true);
    } catch (...) {
        molc_error();
        return nullptr;
    }
}

extern "C" EXPORT PyObject* proper_dihedral_existing_py_inst(void* ptr)
{
    Dihedral *d = static_cast<Dihedral*>(ptr);
    try {
        return d->py_instance(false);
    } catch (...) {
        molc_error();
        return nullptr;
    }
}


extern "C" EXPORT void 
proper_dihedral_from_atoms(void *atoms, size_t n, pyobject_t *names, void *residues, pyobject_t *dihedrals)
{
    // n must be divisible by 4
    try {
        if ((n % 4) != 0) {
            throw std::invalid_argument("Number of atoms must be a multiple of 4");
        }
        Atom **a = static_cast<Atom **>(atoms);
        Residue **r = static_cast<Residue **>(residues);
        for (size_t i=0, j=0; i<n; i+=4, j++) {
            Proper_Dihedral* d = new Proper_Dihedral(a[i], a[i+1], a[i+2], a[i+3], r[j],
                std::string(PyUnicode_AsUTF8(static_cast<PyObject *>(names[j]))));
            //d->set_name(std::string(PyUnicode_AsUTF8(static_cast<PyObject *>(names[j]))));
            dihedrals[j] = d;
        }
    } catch (...) {
        molc_error();
    }
}



extern "C" EXPORT void 
proper_dihedral_angle(void *dihedrals, size_t n, float32_t *angles)
{
    Dihedral **d = static_cast<Dihedral **>(dihedrals);
    error_wrap_array_get(d, n, &Dihedral::angle, angles);
}

extern "C" EXPORT void 
proper_dihedral_name(void *dihedrals, size_t n, pyobject_t *names)
{
    Dihedral **d = static_cast<Dihedral **>(dihedrals);
    try {
        for (size_t i = 0; i<n; ++i)
            names[i] = unicode_from_string(d[i]->name());
    } catch (...) {
        molc_error();
    }
}

//~ extern "C" EXPORT void 
//~ set_dihedral_name(void *dihedrals, size_t n, pyobject_t *names)
//~ {
    //~ Dihedral **d = static_cast<Dihedral **>(dihedrals);
    //~ try {
        //~ for (size_t i=0; i<n; ++i)
            //~ d[i]->set_name(std::string(PyUnicode_AsUTF8(static_cast<PyObject *>(names[i]))));
    //~ } catch (...) {
        //~ molc_error();
    //~ }
//~ }


typedef Dihedral_Mgr<Proper_Dihedral> Proper_Dihedral_Mgr;

extern "C" EXPORT void set_proper_dihedral_mgr_py_instance(void* mgr, PyObject* py_inst)
{
    Proper_Dihedral_Mgr *d = static_cast<Proper_Dihedral_Mgr *>(mgr);
    try {
        d->set_py_instance(py_inst);
    } catch (...) {
        molc_error();
    }
}


extern "C" EXPORT PyObject* proper_dihedral_mgr_py_inst(void* ptr)
{
    Proper_Dihedral_Mgr *d = static_cast<Proper_Dihedral_Mgr *>(ptr);
    try {
        return d->py_instance(true);
    } catch (...) {
        molc_error();
        return nullptr;
    }
}

extern "C" EXPORT PyObject* proper_dihedral_mgr_existing_py_inst(void* ptr)
{
    Proper_Dihedral_Mgr *d = static_cast<Proper_Dihedral_Mgr *>(ptr);
    try {
        return d->py_instance(false);
    } catch (...) {
        molc_error();
        return nullptr;
    }
}






extern "C" EXPORT void*
proper_dihedral_mgr_new()
{
    try {
        Proper_Dihedral_Mgr *mgr = new Proper_Dihedral_Mgr();
        return mgr;
        
    } catch (...) {
        molc_error();
        return nullptr;
    }
}

extern "C" EXPORT void
proper_dihedral_mgr_delete(void *mgr)
{
    try {
        Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
        delete m;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_mgr_add_dihedral(void *mgr, void *dihedrals, size_t n)
{
    try {
        Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
        Proper_Dihedral **d = static_cast<Proper_Dihedral **>(dihedrals);
        for (size_t i=0; i<n; ++i)
            m->add_dihedral(d[i]);
    } catch (...) {
        molc_error();
    }
        
}

extern "C" EXPORT int
proper_dihedral_mgr_get_dihedrals(void *mgr, void *residues, pyobject_t *names, size_t n, pyobject_t *dihedrals)
{
    size_t found=0;
    try {
        Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
        Residue **r = static_cast<Residue **>(residues);
        for (size_t i=0, found=0; i<n; ++i) {
            try {
                dihedrals[found++] = m->get_dihedral(r[i], std::string(PyUnicode_AsUTF8(static_cast<PyObject *>(names[i]))));
            } catch (std::out_of_range) {
                continue;
            }
        }
        return found;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT int
proper_dihedral_mgr_num_dihedrals(void *mgr)
{
    try {
        Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
        return m->num_dihedrals();
    } catch (...) {
        molc_error();
        return 0;
    }
}
        
        

