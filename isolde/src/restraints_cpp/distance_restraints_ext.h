/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 28-Mar-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#ifndef ISOLDE_DISTANCE_RESTRAINTS_EXT
#define ISOLDE_DISTANCE_RESTRAINTS_EXT

#include "distance_restraints.h"
#include "../molc.h"

using namespace atomstruct;
using namespace isolde;

/*******************************************************
 *
 * DistanceRestraintMgr functions
 *
 *******************************************************/

SET_PYTHON_INSTANCE(distance_restraint_mgr, DistanceRestraintMgr)
GET_PYTHON_INSTANCES(distance_restraint_mgr, DistanceRestraintMgr)

extern "C" EXPORT void*
distance_restraint_mgr_new(void *structure, void *change_tracker)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    try {
        DistanceRestraintMgr *mgr = new DistanceRestraintMgr(s, ct);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //distance_restraint_mgr_new

extern "C" EXPORT void
distance_restraint_mgr_delete(void *mgr)
{
    DistanceRestraintMgr *m = static_cast<DistanceRestraintMgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
}

// extern "C" EXPORT void*
// distance_restraint_mgr_get_restraint(void *mgr, void *atoms, bool create)
// {
//     DistanceRestraintMgr *d = static_cast<DistanceRestraintMgr *>(mgr);
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
distance_restraint_mgr_get_restraint(void *mgr, void* atom1, void* atom2, npy_bool create, size_t n, pyobject_t *restraint)
{
    DistanceRestraintMgr *d = static_cast<DistanceRestraintMgr *>(mgr);
    Atom **a1 = static_cast<Atom **>(atom1);
    Atom **a2 = static_cast<Atom **>(atom2);
    size_t count=0;
    try {
        for (size_t i=0; i<n; ++i) {
            DistanceRestraint *r = d->get_restraint(*a1++, *a2++, create);
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
distance_restraint_mgr_atom_restraints(void *mgr, void *atom, size_t n)
{
    DistanceRestraintMgr *d = static_cast<DistanceRestraintMgr *>(mgr);
    Atom **a = static_cast<Atom **>(atom);
    try {
        std::set<DistanceRestraint *> dset;
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
    DistanceRestraintMgr *d = static_cast<DistanceRestraintMgr *>(mgr);
    Atom **a = static_cast<Atom **>(atoms);
    try {
        std::set<Atom *> aset;
        std::set<DistanceRestraint *> dset;
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
    DistanceRestraintMgr *d = static_cast<DistanceRestraintMgr *>(mgr);
    try {
        std::vector<DistanceRestraint *> visibles;
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
    DistanceRestraintMgr *d = static_cast<DistanceRestraintMgr *>(mgr);
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

//NOTE: Residues are expected to already be sorted in chains
extern "C" EXPORT PyObject*
distance_restraint_mgr_get_ss_restraints(void *mgr, void *residues, size_t n, bool create)
{
    DistanceRestraintMgr *d = static_cast<DistanceRestraintMgr *>(mgr);
    Residue **r = static_cast<Residue **>(residues);
    PyObject* ret = PyTuple_New(2);
    try {
        if (n < 3) {
            throw std::logic_error("Secondary structure restraints require at least three contiguous residues!");
        }
        std::vector<DistanceRestraint *> o_to_n_plus_four;
        std::vector<DistanceRestraint *> ca_to_ca_plus_two;
        for (size_t i=0; i<n-2; ++i)
        {
            Residue* cr = *r++;
            if (cr->polymer_type() != PT_AMINO) {
                continue;
            }
            Atom* cca = cr->find_atom("CA");
            Atom* co = cr->find_atom("O");
            if (cca==nullptr || co==nullptr) {
                continue;
                //throw std::logic_error("Missing backbone atoms detected!");
            }
            Residue* rp1 = *r;
            if (!(cr->connects_to(rp1)) || !(rp1->polymer_type()==PT_AMINO)) {
                continue;
            }
            Residue* rp2 = *(r+1);
            if (!rp1->connects_to(rp2) || !(rp2->polymer_type()==PT_AMINO)) {
                continue;
            }
            Atom* cap2 = rp2->find_atom("CA");
            if (cap2 != nullptr) {
                DistanceRestraint* cad = d->get_restraint(cca, cap2, create);
                if (cad != nullptr)
                    ca_to_ca_plus_two.push_back(cad);
            }
            if (i+4 >= n) continue;
            Residue* rp3 = *(r+2);
            if (!rp2->connects_to(rp3) || !(rp3->polymer_type()==PT_AMINO)) continue;
            Residue* rp4 = *(r+3);
            if (!rp3->connects_to(rp4) || !(rp4->polymer_type()==PT_AMINO)) continue;
            Atom* np4 = rp4->find_atom("N");
            if (np4 != nullptr) {
                DistanceRestraint* on4 = d->get_restraint(co, np4, create);
                if (on4 != nullptr)
                    o_to_n_plus_four.push_back(on4);
            }
        }
        void **onptrs;
        PyObject* on_restr_array = python_voidp_array(o_to_n_plus_four.size(), &onptrs);
        for (const auto &ptr: o_to_n_plus_four)
            *onptrs++ = ptr;

        void **captrs;
        PyObject* ca_restr_array = python_voidp_array(ca_to_ca_plus_two.size(), &captrs);
        for (const auto &ptr: ca_to_ca_plus_two)
            *captrs++ = ptr;

        PyTuple_SET_ITEM(ret, 0, on_restr_array);
        PyTuple_SET_ITEM(ret, 1, ca_restr_array);
        return ret;
    } catch(...) {
        molc_error();
        Py_XDECREF(ret);
        return 0;
    }
}


/***************************************************************
 *
 * DistanceRestraint functions
 *
 ***************************************************************/
SET_PYTHON_CLASS(distance_restraint, DistanceRestraint)
GET_PYTHON_INSTANCES(distance_restraint, DistanceRestraint)

extern "C" EXPORT void
distance_restraint_get_manager(void *restraint, size_t n, pyobject_t *mgr)
{
    DistanceRestraint **r = static_cast<DistanceRestraint **>(restraint);
    error_wrap_array_get(r, n, &DistanceRestraint::mgr, mgr);
}

extern "C" EXPORT void
set_distance_restraint_target(void *restraint, size_t n, double *target)
{
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
    error_wrap_array_set<DistanceRestraint, double, double>(d, n, &DistanceRestraint::set_target, target);
}

extern "C" EXPORT void
distance_restraint_target(void *restraint, size_t n, double *target)
{
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
    error_wrap_array_get<DistanceRestraint, double, double>(d, n, &DistanceRestraint::get_target, target);
}

extern "C" EXPORT void
set_distance_restraint_k(void *restraint, size_t n, double *k)
{
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
    error_wrap_array_set<DistanceRestraint, double, double>(d, n, &DistanceRestraint::set_k, k);
}

extern "C" EXPORT void
distance_restraint_k(void *restraint, size_t n, double *k)
{
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
    error_wrap_array_get<DistanceRestraint, double, double>(d, n, &DistanceRestraint::get_k, k);
}

extern "C" EXPORT void
set_distance_restraint_satisfied_limit(void *restraint, size_t n, double *limit)
{
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
    error_wrap_array_set<DistanceRestraint, double, double>(d, n, &DistanceRestraint::set_satisfied_limit, limit);
}

extern "C" EXPORT void
distance_restraint_satisfied_limit(void *restraint, size_t n, double *limit)
{
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
    error_wrap_array_get<DistanceRestraint, double, double>(d, n, &DistanceRestraint::get_satisfied_limit, limit);
}

extern "C" EXPORT void
distance_restraint_satisfied(void *restraint, size_t n, npy_bool *satisfied)
{
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
    error_wrap_array_get<DistanceRestraint, bool, npy_bool>(d, n, &DistanceRestraint::satisfied, satisfied);
}

extern "C" EXPORT void
distance_restraint_unsatisfied(void *restraint, size_t n, npy_bool *unsatisfied)
{
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            *(unsatisfied++) = !( (*d++)->satisfied() );
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
distance_restraint_center(void *restraint, size_t n, double *coords)
{
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
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
distance_restraint_sim_index(void *restraint, size_t n, int *index)
{
    DistanceRestraint **r = static_cast<DistanceRestraint **>(restraint);
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
    DistanceRestraint **r = static_cast<DistanceRestraint **>(restraint);
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
    DistanceRestraint **r = static_cast<DistanceRestraint **>(restraint);
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
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
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
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
    error_wrap_array_get<DistanceRestraint, double, double>(d, n, &DistanceRestraint::distance, distance);
}

extern "C" EXPORT void
distance_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
    error_wrap_array_get<DistanceRestraint, bool, npy_bool>(d, n, &DistanceRestraint::enabled, flag);
}

extern "C" EXPORT void
set_distance_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
    error_wrap_array_set<DistanceRestraint, bool, npy_bool>(d, n, &DistanceRestraint::set_enabled, flag);
}

extern "C" EXPORT void
distance_restraint_visible(void *restraint, size_t n, npy_bool *flag)
{
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
    error_wrap_array_get<DistanceRestraint, bool, npy_bool>(d, n, &DistanceRestraint::visible, flag);
}

extern "C" EXPORT void
distance_restraint_bond_transform(void *restraint, size_t n, float *rot44)
{
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
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
distance_restraint_target_transform(void *restraint, size_t n, float *rot44)
{
    DistanceRestraint **d = static_cast<DistanceRestraint **>(restraint);
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
