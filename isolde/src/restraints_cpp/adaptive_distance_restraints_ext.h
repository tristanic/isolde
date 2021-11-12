/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 09-Apr-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#pragma once

#include "adaptive_distance_restraints.h"
#include "../molc.h"

using namespace atomstruct;
using namespace isolde;

/*******************************************************
 *
 * AdaptiveDistanceRestraintMgr functions
 *
 *******************************************************/

SET_PYTHON_INSTANCE(adaptive_distance_restraint_mgr, AdaptiveDistanceRestraintMgr)
GET_PYTHON_INSTANCES(adaptive_distance_restraint_mgr, AdaptiveDistanceRestraintMgr)

extern "C" EXPORT void*
adaptive_distance_restraint_mgr_new(void *structure, void *change_tracker)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    try {
        AdaptiveDistanceRestraintMgr *mgr = new AdaptiveDistanceRestraintMgr(s, ct);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //distance_restraint_mgr_new

extern "C" EXPORT void
adaptive_distance_restraint_mgr_delete(void *mgr)
{
    AdaptiveDistanceRestraintMgr *m = static_cast<AdaptiveDistanceRestraintMgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
}

// extern "C" EXPORT void*
// adaptive_distance_restraint_mgr_get_restraint(void *mgr, void *atoms, bool create)
// {
//     AdaptiveDistanceRestraintMgr *d = static_cast<AdaptiveDistanceRestraintMgr *>(mgr);
//     Atom **a = static_cast<Atom **>(atoms);
//     bool c = (bool)create;
//     try {
//         return d->get_restraint(*a, *(a+1), c);
//     } catch (...) {
//         molc_error();
//         return 0;
//     }
// }

extern "C" EXPORT size_t
adaptive_distance_restraint_mgr_get_restraint(void *mgr, void* atom1, void* atom2, npy_bool create, size_t n, pyobject_t *restraint)
{
    AdaptiveDistanceRestraintMgr *d = static_cast<AdaptiveDistanceRestraintMgr *>(mgr);
    Atom **a1 = static_cast<Atom **>(atom1);
    Atom **a2 = static_cast<Atom **>(atom2);
    size_t count=0;
    try {
        for (size_t i=0; i<n; ++i) {
            AdaptiveDistanceRestraint *r = d->get_restraint(*a1++, *a2++, create);
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


extern "C" EXPORT PyObject*
adaptive_distance_restraint_mgr_atom_restraints(void *mgr, void *atom, size_t n)
{
    AdaptiveDistanceRestraintMgr *d = static_cast<AdaptiveDistanceRestraintMgr *>(mgr);
    Atom **a = static_cast<Atom **>(atom);
    try {
        std::set<AdaptiveDistanceRestraint *> dset;
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
    AdaptiveDistanceRestraintMgr *d = static_cast<AdaptiveDistanceRestraintMgr *>(mgr);
    Atom **a = static_cast<Atom **>(atoms);
    try {
        std::set<Atom *> aset;
        std::set<AdaptiveDistanceRestraint *> dset;
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
    AdaptiveDistanceRestraintMgr *d = static_cast<AdaptiveDistanceRestraintMgr *>(mgr);
    try {
        return d->display_threshold();
    } catch (...) {
        molc_error();
        return std::nan("");
    }
}

extern "C" EXPORT void
set_adaptive_distance_restraint_mgr_display_threshold(void *mgr, double t) {
    AdaptiveDistanceRestraintMgr *d = static_cast<AdaptiveDistanceRestraintMgr *>(mgr);
    try {
        d->set_display_threshold(t);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT PyObject*
adaptive_distance_restraint_mgr_visible_restraints(void *mgr)
{
    AdaptiveDistanceRestraintMgr *d = static_cast<AdaptiveDistanceRestraintMgr *>(mgr);
    try {
        std::vector<AdaptiveDistanceRestraint *> visibles;
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
    AdaptiveDistanceRestraintMgr *d = static_cast<AdaptiveDistanceRestraintMgr *>(mgr);
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
    AdaptiveDistanceRestraintMgr *d = static_cast<AdaptiveDistanceRestraintMgr *>(mgr);
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
 * AdaptiveDistanceRestraint functions
 *
 ***************************************************************/
SET_PYTHON_CLASS(adaptive_distance_restraint, AdaptiveDistanceRestraint)
GET_PYTHON_INSTANCES(adaptive_distance_restraint, AdaptiveDistanceRestraint)

extern "C" EXPORT void
adaptive_distance_restraint_get_manager(void *restraint, size_t n, pyobject_t *mgr)
{
    AdaptiveDistanceRestraint **r = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_get(r, n, &AdaptiveDistanceRestraint::mgr, mgr);
}

extern "C" EXPORT void
set_adaptive_distance_restraint_target(void *restraint, size_t n, double *target)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_set<AdaptiveDistanceRestraint, const double&, double>(d, n, &AdaptiveDistanceRestraint::set_target, target);
}

extern "C" EXPORT void
adaptive_distance_restraint_target(void *restraint, size_t n, double *target)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_get<AdaptiveDistanceRestraint, double, double>(d, n, &AdaptiveDistanceRestraint::get_target, target);
}

extern "C" EXPORT void
set_adaptive_distance_restraint_tolerance(void *restraint, size_t n, double *tol)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_set<AdaptiveDistanceRestraint, const double&, double>(d, n, &AdaptiveDistanceRestraint::set_tolerance, tol);
}

extern "C" EXPORT void
adaptive_distance_restraint_tolerance(void *restraint, size_t n, double *tol)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_get<AdaptiveDistanceRestraint, double, double>(d, n, &AdaptiveDistanceRestraint::get_tolerance, tol);
}

extern "C" EXPORT void
set_adaptive_distance_restraint_kappa(void *restraint, size_t n, double *kappa)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_set<AdaptiveDistanceRestraint, const double&, double>(d, n, &AdaptiveDistanceRestraint::set_kappa, kappa);
}

extern "C" EXPORT void
adaptive_distance_restraint_kappa(void *restraint, size_t n, double *kappa)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_get<AdaptiveDistanceRestraint, double, double>(d, n, &AdaptiveDistanceRestraint::get_kappa, kappa);
}

extern "C" EXPORT void
set_adaptive_distance_restraint_c(void *restraint, size_t n, double *c)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_set<AdaptiveDistanceRestraint, const double&, double>(d, n, &AdaptiveDistanceRestraint::set_c, c);
}

extern "C" EXPORT void
adaptive_distance_restraint_c(void *restraint, size_t n, double *c)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_get<AdaptiveDistanceRestraint, double, double>(d, n, &AdaptiveDistanceRestraint::get_c, c);
}

extern "C" EXPORT void
set_adaptive_distance_restraint_effective_k(void *restraint, size_t n, double *k)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
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
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_get<AdaptiveDistanceRestraint, double, double>(d, n, &AdaptiveDistanceRestraint::effective_spring_constant, k);
}

extern "C" EXPORT void
set_adaptive_distance_restraint_alpha(void *restraint, size_t n, double *alpha)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_set<AdaptiveDistanceRestraint, const double&, double>(d, n, &AdaptiveDistanceRestraint::set_alpha, alpha);
}

extern "C" EXPORT void
adaptive_distance_restraint_alpha(void *restraint, size_t n, double *alpha)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_get<AdaptiveDistanceRestraint, double, double>(d, n, &AdaptiveDistanceRestraint::get_alpha, alpha);
}

extern "C" EXPORT void
set_adaptive_distance_restraint_satisfied_limit(void *restraint, size_t n, double *limit)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_set<AdaptiveDistanceRestraint, double, double>(d, n, &AdaptiveDistanceRestraint::set_satisfied_limit, limit);
}

extern "C" EXPORT void
adaptive_distance_restraint_satisfied_limit(void *restraint, size_t n, double *limit)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_get<AdaptiveDistanceRestraint, double, double>(d, n, &AdaptiveDistanceRestraint::get_satisfied_limit, limit);
}

extern "C" EXPORT void
adaptive_distance_restraint_satisfied(void *restraint, size_t n, npy_bool *satisfied)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_get<AdaptiveDistanceRestraint, bool, npy_bool>(d, n, &AdaptiveDistanceRestraint::satisfied, satisfied);
}

extern "C" EXPORT void
adaptive_distance_restraint_unsatisfied(void *restraint, size_t n, npy_bool *unsatisfied)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            *(unsatisfied++) = !( (*d++)->satisfied() );
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_distance_restraint_center(void *restraint, size_t n, double *coords)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            auto center = (*d++)->center();
            for (size_t j=0; j<3; ++j)
                *coords++ = center[j];
        }
    } catch(...) {
        molc_error();
    }
}


extern "C" EXPORT void
adaptive_distance_restraint_sim_index(void *restraint, size_t n, int *index)
{
    AdaptiveDistanceRestraint **r = static_cast<AdaptiveDistanceRestraint **>(restraint);
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
    AdaptiveDistanceRestraint **r = static_cast<AdaptiveDistanceRestraint **>(restraint);
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
    AdaptiveDistanceRestraint **r = static_cast<AdaptiveDistanceRestraint **>(restraint);
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
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
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
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_get<AdaptiveDistanceRestraint, double, double>(d, n, &AdaptiveDistanceRestraint::force_magnitude, force);
}

extern "C" EXPORT void
adaptive_distance_restraint_distance(void *restraint, size_t n, double *distance)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_get<AdaptiveDistanceRestraint, double, double>(d, n, &AdaptiveDistanceRestraint::distance, distance);
}

extern "C" EXPORT void
adaptive_distance_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_get<AdaptiveDistanceRestraint, bool, npy_bool>(d, n, &AdaptiveDistanceRestraint::enabled, flag);
}

extern "C" EXPORT void
set_adaptive_distance_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_set<AdaptiveDistanceRestraint, bool, npy_bool>(d, n, &AdaptiveDistanceRestraint::set_enabled, flag);
}

extern "C" EXPORT void
adaptive_distance_restraint_visible(void *restraint, size_t n, npy_bool *flag)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
    error_wrap_array_get<AdaptiveDistanceRestraint, bool, npy_bool>(d, n, &AdaptiveDistanceRestraint::visible, flag);
}

extern "C" EXPORT void
adaptive_distance_restraint_bond_transforms(void *restraint, size_t n, float *rot44_ends, float *rot44_m)
{
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
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
    AdaptiveDistanceRestraint **d = static_cast<AdaptiveDistanceRestraint **>(restraint);
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
