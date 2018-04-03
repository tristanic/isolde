#ifndef ISOLDE_DISTANCE_RESTRAINTS_EXT
#define ISOLDE_DISTANCE_RESTRAINTS_EXT

#include "distance_restraints.h"
#include "../molc.h"

using namespace atomstruct;
using namespace isolde;

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
distance_restraint_mgr_atom_restraints(void *mgr, void *atom, size_t n)
{
    Distance_Restraint_Mgr *d = static_cast<Distance_Restraint_Mgr *>(mgr);
    Atom **a = static_cast<Atom **>(atom);
    try {
        std::set<Distance_Restraint *> dset;
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
distance_restraint_sim_index(void *restraint, size_t n, int *index)
{
    Distance_Restraint **r = static_cast<Distance_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            *(index++) = (*r++)->get_sim_index();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_distance_restraint_sim_index(void *restraint, size_t n, int *index)
{
    Distance_Restraint **r = static_cast<Distance_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            (*r++)->set_sim_index(*(index++));
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
distance_restraint_clear_sim_index(void *restraint, size_t n)
{
    Distance_Restraint **r = static_cast<Distance_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            (*r++)->clear_sim_index();
    } catch (...) {
        molc_error();
    }
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

extern "C" EXPORT void
distance_restraint_target_transform(void *restraint, size_t n, double *rot44)
{
    Distance_Restraint **d = static_cast<Distance_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*d++)->target_transform(rot44);
            rot44+=16;
        }
    } catch (...) {
        molc_error();
    }
}


#endif
