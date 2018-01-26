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
#include "validation_new/rama.h"
#include "validation_new/rota.h"

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

SET_PYTHON_CLASS(proper_dihedral, Proper_Dihedral)
GET_PYTHON_INSTANCES(proper_dihedral, Proper_Dihedral)

/************************************************
 *
 * Generic dihedral functions
 *
 ************************************************/

extern "C" EXPORT void 
dihedral_angle(void *dihedrals, size_t n, float32_t *angles)
{
    Dihedral **d = static_cast<Dihedral **>(dihedrals);
    error_wrap_array_get(d, n, &Dihedral::angle, angles);
}

extern "C" EXPORT void 
dihedral_name(void *dihedrals, size_t n, pyobject_t *names)
{
    Dihedral **d = static_cast<Dihedral **>(dihedrals);
    try {
        for (size_t i = 0; i<n; ++i)
            names[i] = unicode_from_string(d[i]->name());
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
dihedral_atoms(void *dihedrals, size_t n, pyobject_t *atoms)
{
    Dihedral **d = static_cast<Dihedral **>(dihedrals);
    try {
        for (size_t i=0; i<n; ++i) {
            const Dihedral::Atoms &a = d[i]->atoms();
            for (auto ta: a) {
                *atoms++ = ta;
            }
        }
    } catch (...) {
        molc_error();
    }
}

 /**************************************************
  * 
  * Proper_Dihedral functions
  * 
  **************************************************/

extern "C" EXPORT void
proper_dihedral_axial_bond(void *dihedrals, size_t n, pyobject_t *bonds)
{
    Dihedral **d = static_cast<Dihedral **>(dihedrals);
    error_wrap_array_get(d, n, &Dihedral::axial_bond, bonds);
}

extern "C" EXPORT void
proper_dihedral_target(void *dihedrals, size_t n, float32_t *vals)
{
    Dihedral **d = static_cast<Dihedral **>(dihedrals);
    error_wrap_array_get(d, n, &Dihedral::target, vals);
}

extern "C" EXPORT void
set_proper_dihedral_target(void * dihedrals, size_t n, float32_t *vals)
{
    Dihedral **d = static_cast<Dihedral **>(dihedrals);
    error_wrap_array_set(d, n, &Dihedral::set_target, vals);
}
    
extern "C" EXPORT void
proper_dihedral_spring_constant(void *dihedrals, size_t n, float32_t *vals)
{
    Dihedral **d = static_cast<Dihedral **>(dihedrals);
    error_wrap_array_get(d, n, &Dihedral::spring_constant, vals);
}

extern "C" EXPORT void
set_proper_dihedral_spring_constant(void * dihedrals, size_t n, float32_t *vals)
{
    Dihedral **d = static_cast<Dihedral **>(dihedrals);
    error_wrap_array_set(d, n, &Dihedral::set_spring_constant, vals);
}
    
extern "C" EXPORT void proper_dihedral_residue(void *dihedrals, size_t n, pyobject_t *resp)
{
    Dihedral **d = static_cast<Dihedral **>(dihedrals);
    error_wrap_array_get(d, n, &Dihedral::residue, resp);
}



 /*************************************
  * 
  * Proper_Dihedral_Mgr functions
  * 
  *************************************/

SET_PYTHON_INSTANCE(proper_dihedral_mgr, Proper_Dihedral_Mgr)
GET_PYTHON_INSTANCES(proper_dihedral_mgr, Proper_Dihedral_Mgr)

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
    Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_mgr_delete_dihedral(void *mgr, size_t n, void *dihedrals)
{
    Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
    Proper_Dihedral **d = static_cast<Proper_Dihedral **>(dihedrals);
    try {
        std::vector<Proper_Dihedral *> delete_list;
        for (size_t i=0; i<n; ++i) {
            delete_list.push_back(d[i]);
        }
        m->delete_dihedrals(delete_list);
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

extern "C" EXPORT void
proper_dihedral_mgr_reserve_map(void *mgr, size_t n)
{
    Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
    try {
        m->reserve(n);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_mgr_new_dihedral(void *mgr, void*residues, size_t n, pyobject_t *name) 
{
    Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
    Residue **r = static_cast<Residue **>(residues);
    try {
        std::string sname = std::string(PyUnicode_AsUTF8(static_cast<PyObject *>(name[0])));
        for (size_t i=0; i<n; ++i) {
            m->new_dihedral(*r++, sname);
        }
    } catch (...) {
        molc_error();
    }
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
            Proper_Dihedral *d = m->get_dihedral(r[i], dname, (bool)create);
            if (d != nullptr)
                dihedrals[found++] = d;
        }
        return found;
    } catch (...) {
        molc_error();
        return 0;
    }
}
        
extern "C" EXPORT int
proper_dihedral_mgr_num_mapped_dihedrals(void *mgr)
{
    Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
    try {
        return m->num_mapped_dihedrals();
    } catch (...) {
        molc_error();
        return 0;
    }
}

const std::vector<std::string> RAMA_DIHEDRAL_NAMES({"omega", "phi", "psi"}); 
extern "C" EXPORT size_t
proper_dihedral_mgr_valid_rama_residues(void *mgr, void *in_residues, size_t n, 
    pyobject_t *out_residues, pyobject_t *omega, pyobject_t *phi, pyobject_t *psi)
{
    size_t found = 0;
    bool found_all = true;
    Proper_Dihedral* found_dihedrals[3];
    Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
    Residue **r = static_cast<Residue **>(in_residues);
    try {
        for (size_t i=0; i<n; ++i) {
            found_all=true;
            for(size_t j=0; j<3; ++j) {
                Proper_Dihedral *d = m->get_dihedral(r[i], RAMA_DIHEDRAL_NAMES[j], true);
                if (d == nullptr) {
                    found_all = false;
                    break;
                }
                found_dihedrals[j] = d;
            }
            if (found_all) {
                out_residues[found] = r[i];
                omega[found] = found_dihedrals[0];
                phi[found] = found_dihedrals[1];
                psi[found] = found_dihedrals[2];
                found++;
            }
        }
        return found;
    } catch (...) {
        molc_error();
        return 0;
    }
} //proper_dihedral_mgr_valid_rama_residues
                    

/**********************************************************************
 * 
 * Rama_Mgr
 * 
 **********************************************************************/
SET_PYTHON_INSTANCE(rama_mgr, Rama_Mgr)
GET_PYTHON_INSTANCES(rama_mgr, Rama_Mgr)

extern "C" EXPORT void*
rama_mgr_new()
{
    Rama_Mgr *mgr = new Rama_Mgr();
    try {
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
}

extern "C" EXPORT void
rama_mgr_delete(void *mgr)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_mgr_add_interpolator(void *mgr, size_t r_case, size_t dim, 
    uint32_t *n, double *min, double *max, double *data)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    try {
        m->add_interpolator(r_case, dim, n, min, max, data);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_mgr_rama_cases(void *mgr, void *omega, void *psi,
    size_t n, uint8_t *cases)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    Dihedral **o = static_cast<Dihedral **>(omega);
    Dihedral **q = static_cast<Dihedral **>(psi);
    try {
        for (size_t i=0; i<n; ++i) {
            *cases++ = m->rama_case(*o++, *q++);
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_mgr_validate(void *mgr, void *residue, void *omega, void *phi,
    void *psi, uint8_t *r_case, size_t n, double *scores)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    Residue **r = static_cast<Residue **>(residue);
    Dihedral **o = static_cast<Dihedral **>(omega);
    Dihedral **p = static_cast<Dihedral **>(phi);
    Dihedral **q = static_cast<Dihedral **>(psi);
    try {
        m->validate(r, o, p, q, r_case, n, scores);
    } catch (...) {
        molc_error();
    }
}
    
/**************************************************************
 * 
 * Rota_Mgr functions
 * 
 **************************************************************/
SET_PYTHON_INSTANCE(rota_mgr, Rota_Mgr)
GET_PYTHON_INSTANCES(rota_mgr, Rota_Mgr)

extern "C" EXPORT void*
rota_mgr_new()
{
    Rota_Mgr *mgr = new Rota_Mgr();
    try {
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
}
 
 
 /*************************************************************
  * 
  * Rotamer functions
  * 
  *************************************************************/

SET_PYTHON_CLASS(rotamer, Rotamer)
GET_PYTHON_INSTANCES(rotamer, Rotamer)
