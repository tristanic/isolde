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
        //~ Proper_Dihedral_Mgr::d_def ddef;
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
rama_mgr_set_cutoffs(void *mgr, size_t r_case, double allowed, double outlier)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    try {
        m->set_cutoffs(r_case, allowed, outlier);
    } catch (...) {
        molc_error();
    }
}
    
extern "C" EXPORT void
rama_mgr_get_cutoffs(void *mgr, size_t r_case, double* cutoffs)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    try {
        auto c = m->get_cutoffs(r_case);
        cutoffs[0] = c->allowed;
        cutoffs[1] = c->outlier;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_mgr_set_color_scale(void *mgr, uint8_t *max, uint8_t *mid, uint8_t *min, uint8_t *na)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    try {
        m->set_colors(max, mid, min, na);
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

extern "C" EXPORT void
rama_mgr_validate_and_color(void *mgr, void *residue, void *omega, 
    void *phi, void *psi, uint8_t *r_case, size_t n, uint8_t *colors)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    std::vector<double> scores(n);
    rama_mgr_validate(mgr, residue, omega, phi, psi, r_case, n, scores.data());
    try {
        m->color_by_scores(scores.data(), r_case, n, colors);
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
rota_mgr_new(void *dihedral_mgr)
{
    Proper_Dihedral_Mgr *dmgr = static_cast<Proper_Dihedral_Mgr *>(dihedral_mgr);
    Rota_Mgr *mgr = new Rota_Mgr(dmgr);
    try {
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //rota_mgr_new

extern "C" EXPORT void
rota_mgr_add_rotamer_def(void *mgr, pyobject_t *resname, size_t n_chi, npy_bool symmetric)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        std::string rname(PyUnicode_AsUTF8(static_cast<PyObject *>(resname[0])));
        m->add_rotamer_def(rname, n_chi, (bool)symmetric);
    } catch (...) {
        molc_error();
    }
} //rota_mgr_add_rotamer_def

extern "C" EXPORT void
rota_mgr_add_interpolator(void *mgr, pyobject_t *resname, size_t dim,
    uint32_t *n, double *min, double *max, double*data)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        std::string rname(PyUnicode_AsUTF8(static_cast<PyObject *>(resname[0])));
        m->add_interpolator(rname, dim, n, min, max, data);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT size_t
rota_mgr_get_rotamer(void *mgr, void *residue, size_t n, pyobject_t *rotamers)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    Residue **r = static_cast<Residue **>(residue);
    size_t found =0;
    try {
        for (size_t i=0; i<n; ++i) {
            Rotamer* rot = m->get_rotamer(r[i]);
            if (rot != nullptr) {
                rotamers[found++] = rot;
            }
        }
        return found;
    } catch (...) {
        molc_error();
        return 0;
    }
} //rota_mgr_get_rotamer

extern "C" EXPORT void
rota_mgr_validate_rotamer(void *mgr, void *rotamer, size_t n, double *scores)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    try {
        m->validate(r, n, scores);
    } catch (...) {
        molc_error();
    }
} //rota_mgr_validate_rotamer

extern "C" EXPORT void
rota_mgr_validate_residue(void *mgr, void *residue, size_t n, double *scores)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    Residue **r = static_cast<Residue **>(residue);
    try {
        m->validate(r,n,scores);
    } catch (...) {
        molc_error();
    }
} //rotamer_mgr_validate_residue
 /*************************************************************
  *
  * Rotamer functions
  *
  *************************************************************/

SET_PYTHON_CLASS(rotamer, Rotamer)
GET_PYTHON_INSTANCES(rotamer, Rotamer)
extern "C" EXPORT void 
rotamer_score(void *rotamer, size_t n, float32_t *score)
{
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    error_wrap_array_get(r, n, &Rotamer::score, score);
}


extern "C" EXPORT void
rotamer_residue(void *rotamer, size_t n, pyobject_t *residue)
{
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    error_wrap_array_get(r, n, &Rotamer::residue, residue);
}

extern "C" EXPORT void 
rotamer_ca_cb_bond(void *rotamer, size_t n, pyobject_t *bond)
{
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    error_wrap_array_get(r, n, &Rotamer::ca_cb_bond, bond);
}

extern "C" EXPORT void
rotamer_num_chi(void *rotamer, size_t n, uint8_t *nchi)
{
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    error_wrap_array_get(r, n, &Rotamer::n_chi, nchi);
}

extern "C" EXPORT void
rotamer_angles(void *rotamer, double *a) 
{
    Rotamer *r = static_cast<Rotamer *>(rotamer);
    try {
        r->angles(a);
    } catch (...) {
        molc_error();
    }
}
    
    

