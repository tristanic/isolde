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

extern "C" EXPORT void proper_dihedral_residue(void *dihedrals, size_t n, pyobject_t *resp)
{
    Dihedral **d = static_cast<Dihedral **>(dihedrals);
    error_wrap_array_get(d, n, &Dihedral::residue, resp);
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
proper_dihedral_mgr_add_dihedral_def(void *mgr, pyobject_t *rname, 
    pyobject_t *dname, pyobject_t *anames, npy_bool *externals)
{
    Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
    try {
        std::string resname(PyUnicode_AsUTF8(static_cast<PyObject *>(rname[0])));
        std::string dihe_name(PyUnicode_AsUTF8(static_cast<PyObject *>(dname[0])));
        std::vector<std::string> atom_names;
        std::vector<bool> externals_bool;
        for (size_t i=0; i<4; ++i) {
            atom_names.push_back(std::string(PyUnicode_AsUTF8(static_cast<PyObject *>(anames[i]))));
            externals_bool.push_back(externals[i]);
        }
        m->add_dihedral_def(resname, dihe_name, atom_names, externals_bool);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_mgr_add_dihedral(void *mgr, void *dihedrals, size_t n)
{
    Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
    Proper_Dihedral **d = static_cast<Proper_Dihedral **>(dihedrals);
    try {
        for (size_t i=0; i<n; ++i)
            m->add_dihedral(d[i]);
    } catch (...) {
        molc_error();
    }
        
}

Proper_Dihedral* _proper_dihedral_mgr_new_dihedral(Proper_Dihedral_Mgr *m,
    Residue *this_residue, const std::string &name, const std::vector<std::string> &anames, 
    const std::vector<bool> &external_atoms, size_t first_internal_atom)
{
    Atom* found_atoms[4];
    Atom* this_atom;
    bool found=false;
    for (auto a: this_residue->atoms()) {
        if (a->name() == anames[first_internal_atom]) {
            found=true;
            found_atoms[first_internal_atom] = a;
            this_atom = a;
            break;
        }
    }
    if (!found) return nullptr;
    
    // Work backwards if necessary
    for (size_t j=first_internal_atom; j>0; j--) {
        found=false;
        for (auto a: this_atom->neighbors()) {
            if (a->name() == anames[j-1]) {
                found=true;
                found_atoms[j-1] = a;
                this_atom = a;
                break;
            }
        }
        if (!found) {
            break;
        }
    }
    if (!found) return nullptr;
    // ...and now work forwards
    this_atom = found_atoms[first_internal_atom];
    for (size_t j=first_internal_atom; j<3; j++) {
        found=false;
        size_t k = j+1;
        for (auto a: this_atom->neighbors()) {
            if ((!external_atoms[k] && a->residue() == this_residue) ||
                 (external_atoms[k] && a->residue() != this_residue)) {
                if (a->name() == anames[k]) {
                    found=true;
                    found_atoms[k] = a;
                    this_atom = a;
                    break;
                }
            }
        }
        if (!found) {
            break; 
        }
        
    }
    if (found) {
        Proper_Dihedral *d = new Proper_Dihedral(found_atoms[0], 
            found_atoms[1], found_atoms[2], found_atoms[3], 
            this_residue, name);
        m->add_dihedral(d);
        return d;
    }
    return nullptr;
                
}
int default_external_atoms[4] = {0,0,0,0};
//! Find the atoms corresponding to a named dihedral for each residue
/*! 
 * The dihedral can optionally span more than one residue. In that case,
 * the optional external_atoms argument should be used with a value of 1 
 * indicating that the corresponding atom should be outside of the 
 * target residue.
 */ 
extern "C" EXPORT void
proper_dihedral_mgr_new_multi_residue_dihedral(void *mgr, void *residues, size_t n, 
    pyobject_t *name, pyobject_t *atom_names, 
    int *external_atoms = default_external_atoms)
{  
    Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
    Residue **r = static_cast<Residue **>(residues);
    size_t first_internal_atom;
    try {
        std::string sname = std::string(PyUnicode_AsUTF8(static_cast<PyObject *>(name[0])));
        std::vector<std::string> anames;
        std::vector<bool> e_bool;
        for (size_t i=0; i<4; ++i) {
            anames.push_back(std::string(PyUnicode_AsUTF8(static_cast<PyObject *>(atom_names[i]))));
            e_bool.push_back((bool)external_atoms[i]);
        }
        for (first_internal_atom=0;; first_internal_atom++) {
            if (first_internal_atom >=4)
                throw std::runtime_error("At least one atom must be inside the target residue!");
            if (!external_atoms[first_internal_atom])
                break;
        }
        for (size_t i=0; i<n; ++i)
            _proper_dihedral_mgr_new_dihedral(m, r[i], sname, anames, e_bool, first_internal_atom);
    } catch (...) {
        molc_error();
    }
}



extern "C" EXPORT void
proper_dihedral_mgr_new_single_residue_dihedral(void *mgr, void *residues, size_t n, 
    pyobject_t *name, pyobject_t *atom_names)
{
    proper_dihedral_mgr_new_multi_residue_dihedral(mgr, residues, n, name, atom_names);
}

extern "C" EXPORT int
proper_dihedral_mgr_get_dihedrals(void *mgr, void *residues, pyobject_t *name, size_t n, pyobject_t *dihedrals, npy_bool create)
{
    Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
    size_t found=0;
    try {
        Proper_Dihedral_Mgr::d_def ddef;
        Residue **r = static_cast<Residue **>(residues);
        std::string dname(PyUnicode_AsUTF8(static_cast<PyObject *>(name[0])));
        for (size_t i=0; i<n; ++i) {
            try {
                dihedrals[found] = m->get_dihedral(r[i], dname);
                found++;
            } catch (std::out_of_range) {
                // Dihedral doesn't exist. Let's see if we can create it.
                try {
                    ddef = m->get_dihedral_def(std::string(r[i]->name()), dname);
                } catch (std::out_of_range) {
                    throw std::invalid_argument("Dihedral name is invalid for this residue!");
                }
                auto external = ddef.second;
                size_t first_internal_atom = 0;
                for (; first_internal_atom < 4; ++first_internal_atom)
                    if (!external[first_internal_atom]) break;
                    
                Proper_Dihedral *d = _proper_dihedral_mgr_new_dihedral(m, r[i], dname, ddef.first, ddef.second, first_internal_atom);
                if (d != nullptr)    
                    dihedrals[found++] = d;                

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
        
extern "C" EXPORT int
proper_dihedral_mgr_num_mapped_dihedrals(void *mgr)
{
    try {
        Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
        return m->num_mapped_dihedrals();
    } catch (...) {
        molc_error();
        return 0;
    }
}
        

