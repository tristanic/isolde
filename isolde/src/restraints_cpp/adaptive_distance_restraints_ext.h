/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 09-Apr-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */



#pragma once

#include "adaptive_distance_restraints.h"
#include "../molc.h"

using namespace atomstruct;
using namespace isolde;

/*******************************************************
 *
 * Adaptive_Distance_Restraint_Mgr functions
 *
 *******************************************************/

SET_PYTHON_INSTANCE(adaptive_distance_restraint_mgr, Adaptive_Distance_Restraint_Mgr)
GET_PYTHON_INSTANCES(adaptive_distance_restraint_mgr, Adaptive_Distance_Restraint_Mgr)

extern "C" EXPORT void*
adaptive_distance_restraint_mgr_new(void *structure, void *change_tracker)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    try {
        Adaptive_Distance_Restraint_Mgr *mgr = new Adaptive_Distance_Restraint_Mgr(s, ct);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //distance_restraint_mgr_new

extern "C" EXPORT void
adaptive_distance_restraint_mgr_delete(void *mgr)
{
    Adaptive_Distance_Restraint_Mgr *m = static_cast<Adaptive_Distance_Restraint_Mgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void*
adaptive_distance_restraint_mgr_get_restraint(void *mgr, void *atoms, bool create)
{
    Adaptive_Distance_Restraint_Mgr *d = static_cast<Adaptive_Distance_Restraint_Mgr *>(mgr);
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
adaptive_distance_restraint_mgr_atom_restraints(void *mgr, void *atom, size_t n)
{
    Adaptive_Distance_Restraint_Mgr *d = static_cast<Adaptive_Distance_Restraint_Mgr *>(mgr);
    Atom **a = static_cast<Atom **>(atom);
    try {
        std::set<Adaptive_Distance_Restraint *> dset;
        for (size_t i=0; i<n; ++i) {
            auto &drs = d->get_restraints(*a++);
            dset.insert(drs.begin(), drs.end());
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
adaptive_distance_restraint_mgr_intra_restraints(void *mgr, void *atoms, size_t n)
{
    Adaptive_Distance_Restraint_Mgr *d = static_cast<Adaptive_Distance_Restraint_Mgr *>(mgr);
    Atom **a = static_cast<Atom **>(atoms);
    try {
        std::set<Atom *> aset;
        std::set<Adaptive_Distance_Restraint *> dset;
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

extern "C" EXPORT double
adaptive_distance_restraint_mgr_display_threshold(void *mgr) {
    Adaptive_Distance_Restraint_Mgr *d = static_cast<Adaptive_Distance_Restraint_Mgr *>(mgr);
    try {
        return d->display_threshold();
    } catch (...) {
        molc_error();
        return std::nan("");
    }
}

extern "C" EXPORT void
set_adaptive_distance_restraint_mgr_display_threshold(void *mgr, double t) {
    Adaptive_Distance_Restraint_Mgr *d = static_cast<Adaptive_Distance_Restraint_Mgr *>(mgr);
    try {
        d->set_display_threshold(t);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT PyObject*
adaptive_distance_restraint_mgr_visible_restraints(void *mgr)
{
    Adaptive_Distance_Restraint_Mgr *d = static_cast<Adaptive_Distance_Restraint_Mgr *>(mgr);
    try {
        std::vector<Adaptive_Distance_Restraint *> visibles;
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
adaptive_distance_restraint_mgr_all_restraints(void *mgr)
{
    Adaptive_Distance_Restraint_Mgr *d = static_cast<Adaptive_Distance_Restraint_Mgr *>(mgr);
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

extern "C" EXPORT void
set_adaptive_distance_restraint_mgr_colors(void *mgr, uint8_t *colors)
{
    Adaptive_Distance_Restraint_Mgr *d = static_cast<Adaptive_Distance_Restraint_Mgr *>(mgr);
    try {
        auto maxc = colors;
        auto midc = colors+4;
        auto minc = colors+8;

        d->set_colors(maxc, midc, minc);
    } catch (...) {
        molc_error();
    }
}

/***************************************************************
 *
 * Adaptive_Distance_Restraint functions
 *
 ***************************************************************/
SET_PYTHON_CLASS(adaptive_distance_restraint, Adaptive_Distance_Restraint)
GET_PYTHON_INSTANCES(adaptive_distance_restraint, Adaptive_Distance_Restraint)

extern "C" EXPORT void
set_adaptive_distance_restraint_target(void *restraint, size_t n, double *target)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_set<Adaptive_Distance_Restraint, const double&, double>(d, n, &Adaptive_Distance_Restraint::set_target, target);
}

extern "C" EXPORT void
adaptive_distance_restraint_target(void *restraint, size_t n, double *target)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_get<Adaptive_Distance_Restraint, double, double>(d, n, &Adaptive_Distance_Restraint::get_target, target);
}

extern "C" EXPORT void
set_adaptive_distance_restraint_tolerance(void *restraint, size_t n, double *tol)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_set<Adaptive_Distance_Restraint, const double&, double>(d, n, &Adaptive_Distance_Restraint::set_tolerance, tol);
}

extern "C" EXPORT void
adaptive_distance_restraint_tolerance(void *restraint, size_t n, double *tol)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_get<Adaptive_Distance_Restraint, double, double>(d, n, &Adaptive_Distance_Restraint::get_tolerance, tol);
}

extern "C" EXPORT void
set_adaptive_distance_restraint_kappa(void *restraint, size_t n, double *kappa)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_set<Adaptive_Distance_Restraint, const double&, double>(d, n, &Adaptive_Distance_Restraint::set_kappa, kappa);
}

extern "C" EXPORT void
adaptive_distance_restraint_kappa(void *restraint, size_t n, double *kappa)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_get<Adaptive_Distance_Restraint, double, double>(d, n, &Adaptive_Distance_Restraint::get_kappa, kappa);
}

extern "C" EXPORT void
set_adaptive_distance_restraint_c(void *restraint, size_t n, double *c)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_set<Adaptive_Distance_Restraint, const double&, double>(d, n, &Adaptive_Distance_Restraint::set_c, c);
}

extern "C" EXPORT void
adaptive_distance_restraint_c(void *restraint, size_t n, double *c)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_get<Adaptive_Distance_Restraint, double, double>(d, n, &Adaptive_Distance_Restraint::get_c, c);
}

extern "C" EXPORT void
set_adaptive_distance_restraint_effective_k(void *restraint, size_t n, double *k)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            auto dd = *d++;
            dd->set_kappa(*k++ * pow(dd->get_c(), 2));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_distance_restraint_effective_k(void *restraint, size_t n, double *k)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_get<Adaptive_Distance_Restraint, double, double>(d, n, &Adaptive_Distance_Restraint::effective_spring_constant, k);
}

extern "C" EXPORT void
set_adaptive_distance_restraint_alpha(void *restraint, size_t n, double *alpha)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_set<Adaptive_Distance_Restraint, const double&, double>(d, n, &Adaptive_Distance_Restraint::set_alpha, alpha);
}

extern "C" EXPORT void
adaptive_distance_restraint_alpha(void *restraint, size_t n, double *alpha)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_get<Adaptive_Distance_Restraint, double, double>(d, n, &Adaptive_Distance_Restraint::get_alpha, alpha);
}


extern "C" EXPORT void
adaptive_distance_restraint_sim_index(void *restraint, size_t n, int *index)
{
    Adaptive_Distance_Restraint **r = static_cast<Adaptive_Distance_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            *(index++) = (*r++)->get_sim_index();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_adaptive_distance_restraint_sim_index(void *restraint, size_t n, int *index)
{
    Adaptive_Distance_Restraint **r = static_cast<Adaptive_Distance_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            (*r++)->set_sim_index(*(index++));
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_distance_restraint_clear_sim_index(void *restraint, size_t n)
{
    Adaptive_Distance_Restraint **r = static_cast<Adaptive_Distance_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            (*r++)->clear_sim_index();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_distance_restraint_atoms(void *restraint, size_t n, pyobject_t *atoms)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
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
adaptive_distance_restraint_force_magnitude(void *restraint, size_t n, double *force)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_get<Adaptive_Distance_Restraint, double, double>(d, n, &Adaptive_Distance_Restraint::force_magnitude, force);
}

extern "C" EXPORT void
adaptive_distance_restraint_distance(void *restraint, size_t n, double *distance)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_get<Adaptive_Distance_Restraint, double, double>(d, n, &Adaptive_Distance_Restraint::distance, distance);
}

extern "C" EXPORT void
adaptive_distance_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_get<Adaptive_Distance_Restraint, bool, npy_bool>(d, n, &Adaptive_Distance_Restraint::enabled, flag);
}

extern "C" EXPORT void
set_adaptive_distance_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_set<Adaptive_Distance_Restraint, bool, npy_bool>(d, n, &Adaptive_Distance_Restraint::set_enabled, flag);
}

extern "C" EXPORT void
adaptive_distance_restraint_visible(void *restraint, size_t n, npy_bool *flag)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    error_wrap_array_get<Adaptive_Distance_Restraint, bool, npy_bool>(d, n, &Adaptive_Distance_Restraint::visible, flag);
}

extern "C" EXPORT void
adaptive_distance_restraint_bond_transforms(void *restraint, size_t n, float *rot44_ends, float *rot44_m)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*d++)->bond_transforms(rot44_ends, rot44_m, rot44_ends+16);
            rot44_ends += 32; rot44_m += 16;
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_distance_restraint_color(void *restraint, size_t n, uint8_t *color)
{
    Adaptive_Distance_Restraint **d = static_cast<Adaptive_Distance_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            auto dd = *d++;
            dd->color(color);
            color += 4;
        }
    } catch(...) {
        molc_error();
    }
}
