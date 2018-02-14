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
#include "restraints_cpp/changetracker.h"
#include "restraints_cpp/position_restraints.h"
#include "restraints_cpp/distance_restraints.h"
#include "restraints_cpp/dihedral_restraints.h"

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
        std::set<Proper_Dihedral *> delete_list;
        for (size_t i=0; i<n; ++i) {
            delete_list.insert(d[i]);
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


extern "C" EXPORT PyObject*
proper_dihedral_mgr_get_dihedrals(void *mgr, void *residues, pyobject_t *name, size_t n, npy_bool create)
{
    Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
    std::vector<Dihedral *> dvec;
    //size_t found=0;
    try {
        //~ Proper_Dihedral_Mgr::d_def ddef;
        Residue **r = static_cast<Residue **>(residues);
        std::string dname(PyUnicode_AsUTF8(static_cast<PyObject *>(name[0])));
        for (size_t i=0; i<n; ++i) {
            Proper_Dihedral *d = m->get_dihedral(r[i], dname, (bool)create);
            if (d != nullptr)
                dvec.push_back(d);
                //dihedrals[found++] = d;
        }
        void **dptr;
        PyObject *da = python_voidp_array(dvec.size(), &dptr);
        size_t i=0;
        for (auto d: dvec) {
            dptr[i++] = d;
        }
        return da;
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

// const std::vector<std::string> RAMA_DIHEDRAL_NAMES({"omega", "phi", "psi"});
// extern "C" EXPORT size_t
// proper_dihedral_mgr_valid_rama_residues(void *mgr, void *in_residues, size_t n,
//     pyobject_t *out_residues, pyobject_t *omega, pyobject_t *phi, pyobject_t *psi)
// {
//     size_t found = 0;
//     bool found_all = true;
//     Proper_Dihedral* found_dihedrals[3];
//     Proper_Dihedral_Mgr *m = static_cast<Proper_Dihedral_Mgr *>(mgr);
//     Residue **r = static_cast<Residue **>(in_residues);
//     try {
//         for (size_t i=0; i<n; ++i) {
//             found_all=true;
//             for(size_t j=0; j<3; ++j) {
//                 Proper_Dihedral *d = m->get_dihedral(r[i], RAMA_DIHEDRAL_NAMES[j], true);
//                 if (d == nullptr) {
//                     found_all = false;
//                     break;
//                 }
//                 found_dihedrals[j] = d;
//             }
//             if (found_all) {
//                 out_residues[found] = r[i];
//                 omega[found] = found_dihedrals[0];
//                 phi[found] = found_dihedrals[1];
//                 psi[found] = found_dihedrals[2];
//                 found++;
//             }
//         }
//         return found;
//     } catch (...) {
//         molc_error();
//         return 0;
//     }
// } //proper_dihedral_mgr_valid_rama_residues


/**********************************************************************
 *
 * Rama_Mgr
 *
 **********************************************************************/
SET_PYTHON_INSTANCE(rama_mgr, Rama_Mgr)
GET_PYTHON_INSTANCES(rama_mgr, Rama_Mgr)

extern "C" EXPORT void*
rama_mgr_new(void *dmgr)
{
    Proper_Dihedral_Mgr *d = static_cast<Proper_Dihedral_Mgr *>(dmgr);
    Rama_Mgr *mgr = new Rama_Mgr(d);
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
rama_mgr_get_color_scale(void *mgr, uint8_t *max, uint8_t *mid, uint8_t *min, uint8_t *na)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    try {
        auto cmap = m->get_colors(1);
        auto &mapped_colors = cmap->mapped_colors();
        auto &na_color = m->default_color();
        for (size_t i=0; i<4; ++i) {
            *min++ = mapped_colors[0].thecolor[i];
            *mid++ = mapped_colors[1].thecolor[i];
            *max++ = mapped_colors[2].thecolor[i];
            *na++ = na_color[i];
        }
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

extern "C" EXPORT size_t
rama_mgr_interpolator_dim(void *mgr, size_t r_case)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    try {
        auto it = m->get_interpolator(r_case);
        return it->dim();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT void
rama_mgr_interpolator_axis_lengths(void *mgr, size_t r_case, uint32_t *ret)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    try {
        auto it = m->get_interpolator(r_case);
        for (auto l: it->length())
            *ret++ = l;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_mgr_interpolator_minmax(void *mgr, size_t r_case, double *minvals, double *maxvals)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    try {
        auto it = m->get_interpolator(r_case);
        for (auto m: it->min())
            *minvals++ = m;
        for (auto m: it->max())
            *maxvals++ = m;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_mgr_interpolator_values(void *mgr, size_t r_case, double *vals)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    try {
        auto it = m->get_interpolator(r_case);
        for (auto d: it->data())
            *vals++ = d;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT size_t
rama_mgr_get_rama(void *mgr, void *residue, size_t n, pyobject_t *ramas)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    Residue **r = static_cast<Residue **>(residue);
    size_t found=0;
    try {
        for (size_t i=0; i<n; ++i) {
            Residue *thisr = *r++;
            if (thisr->polymer_type() != PT_AMINO)
                continue;
            try {
                Rama *ram = m->get_rama(thisr);
                ramas[found++] = ram;
            } catch (std::logic_error) {
                continue;
            }
        }
        return found;
    } catch (...) {
        molc_error();
        return 0;
    }

}

extern "C" EXPORT void
rama_mgr_rama_case(void *mgr, void *residue, size_t n, uint8_t *rcase)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    Residue **r = static_cast<Residue **>(residue);
    try {
        for (size_t i=0; i<n; ++i) {
            *rcase++ = m->rama_case(*r++);
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_mgr_validate_by_residue(void *mgr, void *residue, size_t n, double *score, uint8_t *rcase)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    Residue **r = static_cast<Residue **>(residue);
    try {
        for (size_t i=0; i<n; ++i) {
            *score++ = m->validate(*r);
            *rcase++ = m->rama_case(*r++);
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_mgr_validate(void *mgr, void *rama, size_t n, double *score, uint8_t *rcase)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    Rama **r = static_cast<Rama **>(rama);
    try {
        m->validate(r, n, score, rcase);
    } catch (...) {
        molc_error();
    }
}

//! Provide an array of colors corresponding to Ramachandran scores
extern "C" EXPORT void
rama_mgr_validate_and_color(void *mgr, void *rama, size_t n, uint8_t *colors)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    Rama **r = static_cast<Rama **>(rama);
    std::vector<double> scores(n);
    std::vector<uint8_t> rcases(n);
    try {
        m->validate(r, n, scores.data(), rcases.data());
        m->color_by_scores(scores.data(), rcases.data(), n, colors);
    } catch (...) {
        molc_error();
    }
}

//! Directly apply colors according to Ramachandran scores to CA atoms.
extern "C" EXPORT void
rama_mgr_validate_and_color_cas(void *mgr, void *rama, size_t n)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    Rama **r = static_cast<Rama **>(rama);
    std::vector<double> scores(n);
    std::vector<uint8_t> rcases(n);
    std::vector<uint8_t> colors(n*4);
    try {
        m->validate(r, n, scores.data(), rcases.data());
        m->color_by_scores(scores.data(), rcases.data(), n, colors.data());
        auto cdata = colors.data();
        for (size_t i=0; i<n; ++i) {
            (*r++)->CA_atom()->set_color(*cdata, *(cdata+1), *(cdata+2), *(cdata+3));
            cdata+=4;
        }
    } catch (...) {
        molc_error();
    }
}


extern "C" EXPORT void
rama_mgr_bin_scores(void *mgr, double *score, uint8_t *r_case, size_t n, int32_t *bin)
{
    Rama_Mgr *m = static_cast<Rama_Mgr *>(mgr);
    try
    {
        for(size_t i=0; i<n; ++i) {
            *bin++ = m->bin_score(*score++, *r_case++);
        }
    } catch (...) {
        molc_error();
    }

}

/**************************************************************
 *
 * Rama functions
 *
 **************************************************************/

SET_PYTHON_CLASS(rama, Rama)
GET_PYTHON_INSTANCES(rama, Rama)

extern "C" EXPORT void
rama_ca_atom(void *rama, size_t n, pyobject_t *atom)
{
    Rama **r = static_cast<Rama **>(rama);
    try {
        for (size_t i=0; i<n; ++i) {
            *atom++ = (*r++)->CA_atom();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_residue(void *rama, size_t n, pyobject_t *residuep)
{
    Rama **r = static_cast<Rama **>(rama);
    error_wrap_array_get(r, n, &Rama::residue, residuep);
}


extern "C" EXPORT void
rama_is_valid(void *rama, size_t n, npy_bool *valid)
{
    Rama **r = static_cast<Rama **>(rama);
    try {
        for (size_t i=0; i<n; ++i) {
            *valid++ = (*r++)->is_valid_rama();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_score(void *rama, size_t n, double *score)
{
    Rama **r = static_cast<Rama **>(rama);
    try {
        for (size_t i=0; i<n; ++i) {
            *score++ = (*r++)->score();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_phipsi(void *rama, size_t n, double *angle)
{
    Rama **r = static_cast<Rama **>(rama);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->phipsi(angle);
            angle +=2;
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_omegaphipsi(void *rama, size_t n, double *angles)
{
    Rama **r = static_cast<Rama **>(rama);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->angles(angles);
            angles +=3;
        }
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

extern "C" EXPORT void
rota_mgr_set_cutoffs(void *mgr, double allowed, double outlier)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        m->set_cutoffs(allowed, outlier);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rota_mgr_get_cutoffs(void *mgr, double* cutoffs)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        auto c = m->get_cutoffs();
        cutoffs[0] = c->allowed;
        cutoffs[1] = c->outlier;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_rota_mgr_color_scale(void *mgr, uint8_t *max, uint8_t *mid, uint8_t *min)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        m->set_colors(max, mid, min);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rota_mgr_color_scale(void *mgr, uint8_t *max, uint8_t *mid, uint8_t *min)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        auto cmap = m->get_colors();
        auto &mapped_colors = cmap->mapped_colors();
        for (size_t i=0; i<4; ++i) {
            *min++ = mapped_colors[0].thecolor[i];
            *mid++ = mapped_colors[1].thecolor[i];
            *max++ = mapped_colors[2].thecolor[i];
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT PyObject*
rota_mgr_get_rotamer(void *mgr, void *residue, size_t n)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    Residue **r = static_cast<Residue **>(residue);
    std::vector<Rotamer *> found;
    try {
        for (size_t i=0; i<n; ++i) {
            Rotamer* rot = m->get_rotamer(*r++);
            if (rot != nullptr) {
                found.push_back(rot);
            }
        }
        void ** rptr;
        PyObject *ra = python_voidp_array(found.size(), &rptr);
        size_t i=0;
        for (auto rot: found)
            rptr[i++] = rot;
        return ra;
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

/*******TESTING***********/
extern "C" EXPORT void
rota_mgr_validate_rotamer_threaded(void *mgr, void *rotamer, size_t n, double *scores)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    try {
        m->validate_threaded(r, n, scores);
    } catch (...) {
        molc_error();
    }
} //rota_mgr_validate_rotamer

extern "C" EXPORT npy_bool
rota_mgr_thread_running(void *mgr) {
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        return m->thread_running();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT npy_bool
rota_mgr_thread_done(void *mgr) {
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        return m->thread_done();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT void
rota_mgr_finalize_thread(void *mgr) {
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        m->finalize_thread();
    } catch (...) {
        molc_error();
    }
}


/******END TESTING***********/


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
} //rota_mgr_validate_residue

extern "C" EXPORT size_t
rota_mgr_validate_scale_and_color_rotamers(void *mgr, void *rotamer, size_t n,
    double max_scale, npy_bool non_favored_only,
    pyobject_t *rot_out, double* scale, uint8_t *color_out)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    size_t ret =0;
    try {
        std::vector<double> scores(n);
        auto cutoffs = m->get_cutoffs();
        const auto &log_allowed = cutoffs->log_allowed;
        const auto &allowed = cutoffs->allowed;
        auto log_range = cutoffs->log_outlier-log_allowed;
        m->validate(r,n,scores.data());
        for (auto s: scores) {
            if (non_favored_only)
                if (s>allowed) {
                    r++;
                    continue;
                }
            *rot_out++ = *r++;
            m->color_by_score(&s, 1, color_out);
            color_out +=4;
            auto log_s = log(s);
            auto this_scale = (log_s-log_allowed)/log_range + 1;
            *scale++ = this_scale>max_scale ? max_scale : this_scale;
            ret++;
        }
        return ret;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT void
rota_mgr_color_by_score(void *mgr, double *score, size_t n, uint8_t *color)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        m->color_by_score(score, n, color);
    } catch (...) {
        molc_error();
    }
} //rota_mgr_color_by_score

extern "C" EXPORT size_t
rota_mgr_non_favored(void *mgr, void *rotamer, size_t n, pyobject_t *bad, double *scores)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    size_t found=0;
    std::vector<double> vscores(n);
    int32_t bin;
    int32_t ALLOWED=m->ALLOWED;
    int32_t OUTLIER=m->OUTLIER;
    try {
        m->validate(r, n, vscores.data());
        for(size_t i=0; i<n; ++i) {
            bin = m->bin_score(vscores[i]);
            if (bin==ALLOWED || bin==OUTLIER)
            {
                bad[found] = r[i];
                scores[found++] = vscores[i];
            }
        }
        return found;
    } catch (...) {
        molc_error();
        return 0;
    }
}


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

/*******************************************************
 *
 * Change_Tracker functions
 *
 *******************************************************/
SET_PYTHON_INSTANCE(change_tracker, Change_Tracker)
GET_PYTHON_INSTANCES(change_tracker, Change_Tracker)

extern "C" EXPORT void*
change_tracker_new()
{
    try {
        Change_Tracker *t = new Change_Tracker();
        return t;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT void
change_tracker_delete(void *tracker)
{
    Change_Tracker *t = static_cast<Change_Tracker *>(tracker);
    try {
        delete t;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
change_tracker_clear(void *tracker)
{
    Change_Tracker *t = static_cast<Change_Tracker *>(tracker);
    try {
        t->clear();
    } catch (...) {
        molc_error();
    }
}

static PyObject* changes_as_py_dict(Change_Tracker *t,
    const Change_Tracker::Ptr_Type_to_Set &changes)
{
    PyObject* changes_data = PyDict_New();
    for (const auto &it1: changes)
    {
        const auto& py_classnames = t->get_python_class_names(it1.first);
        PyObject *mgr_type_key = unicode_from_string(py_classnames.first);
        PyObject* mgr_type_dict = PyDict_New();
        for (const auto &it2: it1.second)
        {
            PyObject *mgr_key = PyLong_FromVoidPtr(const_cast<void*>(it2.first));
            PyObject *change_dict = PyDict_New();
            for (const auto &it3: it2.second)
            {
                PyObject *change_key = unicode_from_string(t->reason_string(it3.first));
                void **ptrs;
                PyObject *ptr_array = python_voidp_array(it3.second.size(), &ptrs);
                for (auto ptr: it3.second)
                    (*ptrs++) = const_cast<void*>(ptr);
                PyDict_SetItem(change_dict, change_key, ptr_array);
                Py_DECREF(change_key);
                Py_DECREF(ptr_array);
            }
            PyDict_SetItem(mgr_type_dict, mgr_key, change_dict);
            Py_DECREF(mgr_key);
            Py_DECREF(change_dict);
        }
        PyDict_SetItem(changes_data, mgr_type_key, mgr_type_dict);
        Py_DECREF(mgr_type_key);
        Py_DECREF(mgr_type_dict);
    }
    return changes_data;
}

extern "C" EXPORT PyObject*
change_tracker_changes(void *tracker)
{
    Change_Tracker *t = static_cast<Change_Tracker *>(tracker);
    try {
        return changes_as_py_dict(t, t->get_changes());
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT PyObject*
change_tracker_reason_names(void *tracker)
{
    Change_Tracker *t = static_cast<Change_Tracker *>(tracker);
    try {
        const auto &reasons = t->all_reason_strings();
        size_t size = reasons.size();
        PyObject* ret = PyTuple_New(size);
        size_t i=0;
        for (const auto &r: reasons) {
            PyTuple_SET_ITEM(ret, i++, unicode_from_string(r.second));
        }
        return ret;
    } catch (...) {
        molc_error();
        return 0;
    }
}

/*******************************************************
 *
 * Position_Restraint_Mgr functions
 *
 *******************************************************/
SET_PYTHON_INSTANCE(position_restraint_mgr, Position_Restraint_Mgr)
GET_PYTHON_INSTANCES(position_restraint_mgr, Position_Restraint_Mgr)

extern "C" EXPORT void*
position_restraint_mgr_new(void *structure, void *change_tracker)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    try {
        Position_Restraint_Mgr *mgr = new Position_Restraint_Mgr(s, ct);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //position_restraint_mgr_new

extern "C" EXPORT void
position_restraint_mgr_delete(void *mgr)
{
    Position_Restraint_Mgr *m = static_cast<Position_Restraint_Mgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
} //position_restraint_mgr_delete

extern "C" EXPORT size_t
position_restraint_mgr_get_restraint(void *mgr, void *atom, npy_bool create, size_t n, pyobject_t *restraint)
{
    Position_Restraint_Mgr *m = static_cast<Position_Restraint_Mgr *>(mgr);
    Atom **a = static_cast<Atom **>(atom);
    size_t count = 0;
    try {
        for (size_t i=0; i<n; ++i) {
            Position_Restraint *r = m->get_restraint(*a++, create);
            if (r != nullptr) {
                *restraint++ = r;
                count++;
            }
        }
        return count;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT size_t
position_restraint_mgr_num_restraints(void *mgr)
{
    Position_Restraint_Mgr *m = static_cast<Position_Restraint_Mgr *>(mgr);
    try {
        return m->num_restraints();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT PyObject*
position_restraint_mgr_visible_restraints(void *mgr)
{
    Position_Restraint_Mgr *m = static_cast<Position_Restraint_Mgr *>(mgr);
    try {
        std::vector<Position_Restraint *> visibles = m->visible_restraints();
        void **vptr;
        PyObject *va = python_voidp_array(visibles.size(), &vptr);
        for (auto v: visibles)
            *vptr++ = v;
        return va;
    } catch (...) {
        molc_error();
        return 0;
    }
}



/***************************************************************
 *
 * Position_Restraint functions
 *
 ***************************************************************/
SET_PYTHON_CLASS(position_restraint, Position_Restraint)
GET_PYTHON_INSTANCES(position_restraint, Position_Restraint)

extern "C" EXPORT void
position_restraint_atom(void *restraint, size_t n, pyobject_t *atom)
{
    Position_Restraint **r = static_cast<Position_Restraint **>(restraint);
    error_wrap_array_get(r, n, &Position_Restraint::atom, atom);
}

extern "C" EXPORT void
position_restraint_target(void *restraint, size_t n, double *target)
{
    Position_Restraint **r = static_cast<Position_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->get_target(target);
            target+=3;
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_position_restraint_target(void *restraint, size_t n, double *target)
{
    Position_Restraint **r = static_cast<Position_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_target(target);
            target+=3;
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_position_restraint_k(void *restraint, size_t n, double *k)
{
    Position_Restraint **r = static_cast<Position_Restraint **>(restraint);
    error_wrap_array_set<Position_Restraint, double, double>(r, n, &Position_Restraint::set_k, k);
}

extern "C" EXPORT void
position_restraint_k(void *restraint, size_t n, double *k)
{
    Position_Restraint **r = static_cast<Position_Restraint **>(restraint);
    error_wrap_array_get<Position_Restraint, double, double>(r, n, &Position_Restraint::get_k, k);
}

extern "C" EXPORT void
position_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    Position_Restraint **r = static_cast<Position_Restraint **>(restraint);
    error_wrap_array_get<Position_Restraint, bool, npy_bool>(r, n, &Position_Restraint::enabled, flag);
}

extern "C" EXPORT void
set_position_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    Position_Restraint **r = static_cast<Position_Restraint **>(restraint);
    error_wrap_array_set<Position_Restraint, bool, npy_bool>(r, n, &Position_Restraint::set_enabled, flag);
}

extern "C" EXPORT void
position_restraint_visible(void *restraint, size_t n, npy_bool *flag)
{
    Position_Restraint **r = static_cast<Position_Restraint **>(restraint);
    error_wrap_array_get<Position_Restraint, bool, npy_bool>(r, n, &Position_Restraint::visible, flag);
}

extern "C" EXPORT void
position_restraint_target_vector(void *restraint, size_t n, double *vec)
{
    Position_Restraint **r = static_cast<Position_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->target_vector(vec);
            vec +=3;
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void position_restraint_bond_transform(void *restraint, size_t n, double *transform)
{
    Position_Restraint **r = static_cast<Position_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->bond_cylinder_transform(transform);
            transform += 16;
        }
    } catch (...) {
        molc_error();
    }
}



/*******************************************************
 *
 * Distance_Restraint_Mgr functions
 *
 *******************************************************/

SET_PYTHON_INSTANCE(distance_restraint_mgr, Distance_Restraint_Mgr)
GET_PYTHON_INSTANCES(distance_restraint_mgr, Distance_Restraint_Mgr)

extern "C" EXPORT void*
distance_restraint_mgr_new(void *structure, void *change_tracker)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    try {
        Distance_Restraint_Mgr *mgr = new Distance_Restraint_Mgr(s, ct);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //distance_restraint_mgr_new

extern "C" EXPORT void
distance_restraint_mgr_delete(void *mgr)
{
    Distance_Restraint_Mgr *m = static_cast<Distance_Restraint_Mgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void*
distance_restraint_mgr_get_restraint(void *mgr, void *atoms, bool create)
{
    Distance_Restraint_Mgr *d = static_cast<Distance_Restraint_Mgr *>(mgr);
    Atom **a = static_cast<Atom **>(atoms);
    bool c = (bool)create;
    try {
        return d->get_restraint(*a, *(a+1), c);
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT PyObject*
distance_restraint_mgr_intra_restraints(void *mgr, void *atoms, size_t n)
{
    Distance_Restraint_Mgr *d = static_cast<Distance_Restraint_Mgr *>(mgr);
    Atom **a = static_cast<Atom **>(atoms);
    try {
        std::set<Atom *> aset;
        std::set<Distance_Restraint *> dset;
        for (size_t i=0; i<n; ++i)
            aset.insert(*a++);
        for (auto ta: aset)
        {
            auto &drs = d->get_restraints(ta);
            for (auto dr: drs) {
                auto &datoms = dr->atoms();
                for (auto datom: datoms) {
                    if (datom != ta) {
                        if (aset.find(datom) != aset.end())
                            dset.insert(dr);
                    }
                }
            }
        }
        void **dptr;
        PyObject *da = python_voidp_array(dset.size(), &dptr);
        for (auto dr: dset)
            *dptr++ = dr;
        return da;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT PyObject*
distance_restraint_mgr_visible_restraints(void *mgr)
{
    Distance_Restraint_Mgr *d = static_cast<Distance_Restraint_Mgr *>(mgr);
    try {
        std::vector<Distance_Restraint *> visibles;
        for (auto r: d->all_restraints()) {
            if (r->visible())
                visibles.push_back(r);
        }
        void **rptr;
        PyObject *ra = python_voidp_array(visibles.size(), &rptr);
        for (auto r: visibles)
            *rptr++ = r;
        return ra;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT PyObject*
distance_restraint_mgr_all_restraints(void *mgr)
{
    Distance_Restraint_Mgr *d = static_cast<Distance_Restraint_Mgr *>(mgr);
    try {
        const auto &restraints = d->all_restraints();
        size_t n = restraints.size();
        void **rptr;
        PyObject *ra = python_voidp_array(n, &rptr);
        for (auto r: restraints)
            *rptr++ = r;
        return ra;
    } catch (...) {
        molc_error();
        return 0;
    }
}


/***************************************************************
 *
 * Distance_Restraint functions
 *
 ***************************************************************/
SET_PYTHON_CLASS(distance_restraint, Distance_Restraint)
GET_PYTHON_INSTANCES(distance_restraint, Distance_Restraint)

extern "C" EXPORT void
set_distance_restraint_target(void *restraint, size_t n, double *target)
{
    Distance_Restraint **d = static_cast<Distance_Restraint **>(restraint);
    error_wrap_array_set<Distance_Restraint, double, double>(d, n, &Distance_Restraint::set_target, target);
}

extern "C" EXPORT void
distance_restraint_target(void *restraint, size_t n, double *target)
{
    Distance_Restraint **d = static_cast<Distance_Restraint **>(restraint);
    error_wrap_array_get<Distance_Restraint, double, double>(d, n, &Distance_Restraint::get_target, target);
}

extern "C" EXPORT void
set_distance_restraint_k(void *restraint, size_t n, double *k)
{
    Distance_Restraint **d = static_cast<Distance_Restraint **>(restraint);
    error_wrap_array_set<Distance_Restraint, double, double>(d, n, &Distance_Restraint::set_k, k);
}

extern "C" EXPORT void
distance_restraint_k(void *restraint, size_t n, double *k)
{
    Distance_Restraint **d = static_cast<Distance_Restraint **>(restraint);
    error_wrap_array_get<Distance_Restraint, double, double>(d, n, &Distance_Restraint::get_k, k);
}

extern "C" EXPORT void
distance_restraint_atoms(void *restraint, size_t n, pyobject_t *atoms)
{
    Distance_Restraint **d = static_cast<Distance_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
        {
            auto &a = (*d++)->atoms();
            *atoms++=a[0];
            *atoms++=a[1];
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
distance_restraint_distance(void *restraint, size_t n, double *distance)
{
    Distance_Restraint **d = static_cast<Distance_Restraint **>(restraint);
    error_wrap_array_get<Distance_Restraint, double, double>(d, n, &Distance_Restraint::distance, distance);
}

extern "C" EXPORT void
distance_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    Distance_Restraint **d = static_cast<Distance_Restraint **>(restraint);
    error_wrap_array_get<Distance_Restraint, bool, npy_bool>(d, n, &Distance_Restraint::enabled, flag);
}

extern "C" EXPORT void
set_distance_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    Distance_Restraint **d = static_cast<Distance_Restraint **>(restraint);
    error_wrap_array_set<Distance_Restraint, bool, npy_bool>(d, n, &Distance_Restraint::set_enabled, flag);
}

extern "C" EXPORT void
distance_restraint_visible(void *restraint, size_t n, npy_bool *flag)
{
    Distance_Restraint **d = static_cast<Distance_Restraint **>(restraint);
    error_wrap_array_get<Distance_Restraint, bool, npy_bool>(d, n, &Distance_Restraint::visible, flag);
}

extern "C" EXPORT void
distance_restraint_bond_transform(void *restraint, size_t n, double *rot44)
{
    Distance_Restraint **d = static_cast<Distance_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*d++)->bond_cylinder_transform(rot44);
            rot44+=16;
        }
    } catch (...) {
        molc_error();
    }
}



/***************************************************************
 *
 * Proper_Dihedral_Restraint_Mgr functions
 *
 ***************************************************************/
SET_PYTHON_INSTANCE(proper_dihedral_restraint_mgr, Proper_Dihedral_Restraint_Mgr)
GET_PYTHON_INSTANCES(proper_dihedral_restraint_mgr, Proper_Dihedral_Restraint_Mgr)

extern "C" EXPORT void*
proper_dihedral_restraint_mgr_new(void *structure, void *change_tracker)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    try {
        Proper_Dihedral_Restraint_Mgr *mgr = new Proper_Dihedral_Restraint_Mgr(s, ct);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //distance_restraint_mgr_new

extern "C" EXPORT void
proper_dihedral_restraint_mgr_delete(void *mgr)
{
    Proper_Dihedral_Restraint_Mgr *m = static_cast<Proper_Dihedral_Restraint_Mgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_mgr_set_colors(void *mgr, uint8_t *maxc, uint8_t *midc, uint8_t *minc)
{
    Proper_Dihedral_Restraint_Mgr *m = static_cast<Proper_Dihedral_Restraint_Mgr *>(mgr);
    try {
        m->set_colors(maxc, midc, minc);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT size_t
proper_dihedral_restraint_mgr_num_restraints(void *mgr)
{
    Proper_Dihedral_Restraint_Mgr *m = static_cast<Proper_Dihedral_Restraint_Mgr *>(mgr);
    try {
        return m->num_restraints();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT size_t
proper_dihedral_restraint_mgr_get_restraint(void *mgr, void *dihedral,
        npy_bool create, size_t n, pyobject_t *restraint)
{
    Proper_Dihedral_Restraint_Mgr *m = static_cast<Proper_Dihedral_Restraint_Mgr *>(mgr);
    Proper_Dihedral **d = static_cast<Proper_Dihedral **>(dihedral);
    try {
        size_t found = 0;
        for (size_t i=0; i<n; ++i) {
            auto r = m->get_restraint(*d++, create);
            if (r!=nullptr) {
                *restraint++ = r;
                found++;
            }
        }
        return found;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT PyObject*
proper_dihedral_restraint_mgr_visible_restraints(void *mgr)
{
    Proper_Dihedral_Restraint_Mgr *m = static_cast<Proper_Dihedral_Restraint_Mgr *>(mgr);
    try {
        auto vis = m->visible_restraints();
        void **rptr;
        PyObject *ra = python_voidp_array(vis.size(), &rptr);
        for (auto r: vis) {
            *rptr++ = r;
        }
        return ra;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_mgr_delete_restraint(void *mgr, void *restraint, size_t n)
{
    Proper_Dihedral_Restraint_Mgr *m = static_cast<Proper_Dihedral_Restraint_Mgr *>(mgr);
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        std::set<Proper_Dihedral_Restraint *> to_delete;
        for (size_t i=0; i<n; ++i) {
            to_delete.insert(*r++);
        }
        m->delete_restraints(to_delete);
    } catch (...) {
        molc_error();
    }
}
/***************************************************************
 *
 * Proper_Dihedral_Restraint functions
 *
 ***************************************************************/

SET_PYTHON_CLASS(proper_dihedral_restraint, Proper_Dihedral_Restraint)
GET_PYTHON_INSTANCES(proper_dihedral_restraint, Proper_Dihedral_Restraint)

extern "C" EXPORT void
proper_dihedral_restraint_target(void *restraint, size_t n, double *target)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(target++) = (*r++)->get_target();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_proper_dihedral_restraint_target(void *restraint, size_t n, double *target)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_target(*(target++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_dihedral(void *restraint, size_t n, pyobject_t *dihedral)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *dihedral++ = (*r++)->get_dihedral();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_offset(void *restraint, size_t n, double *offset)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*offset++) = (*r++)->offset();
        }
    } catch (...) {
        molc_error();
    }

}

extern "C" EXPORT void
proper_dihedral_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(flag++) = (*r++)->is_enabled();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_proper_dihedral_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_enabled(*(flag++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_display(void *restraint, size_t n, npy_bool *flag)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(flag++) = (*r++)->get_display();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_proper_dihedral_restraint_display(void *restraint, size_t n, npy_bool *flag)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_display(*(flag++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_visible(void *restraint, size_t n, npy_bool *flag)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(flag++) = (*r++)->visible();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_k(void *restraint, size_t n, double *spring_constant)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(spring_constant++) = (*r++)->get_spring_constant();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_proper_dihedral_restraint_k(void *restraint, size_t n, double *spring_constant)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_spring_constant(*(spring_constant++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_cutoff(void *restraint, size_t n, double *cutoff)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(cutoff++) = (*r++)->get_cutoff();
        }
    } catch (...) {
        molc_error();
    }
}
extern "C" EXPORT void
set_proper_dihedral_restraint_cutoff(void *restraint, size_t n, double *cutoff)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_cutoff(*(cutoff++));
        }
    } catch (...) {
        molc_error();
    }
}


extern "C" EXPORT void
proper_dihedral_restraint_annotation_transform(void *restraint, size_t n, double *tf1, double *tf2)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        double transforms[32];
        for (size_t i=0; i<n; ++i) {
            (*r++)->get_annotation_transform(transforms);
            double *ttf1 = transforms, *ttf2 = transforms+16;
            for (size_t j=0; j<16; ++j) {
                *(tf1++) = *(ttf1++);
                *(tf2++) = *(ttf2++);
            }
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_annotation_color(void *restraint, size_t n, uint8_t *color)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->get_annotation_color(color);
            color +=4;
        }
    } catch (...) {
        molc_error();
    }
}
