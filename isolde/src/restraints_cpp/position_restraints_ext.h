/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#ifndef ISOLDE_POSITION_RESTRAINTS_EXT
#define ISOLDE_POSITION_RESTRAINTS_EXT

#include "position_restraints.h"
#include "../molc.h"

using namespace atomstruct;
using namespace isolde;

/*******************************************************
 *
 * PositionRestraintMgr functions
 *
 *******************************************************/
SET_PYTHON_INSTANCE(position_restraint_mgr, PositionRestraintMgr)
GET_PYTHON_INSTANCES(position_restraint_mgr, PositionRestraintMgr)

extern "C" EXPORT void*
position_restraint_mgr_new(void *structure, void *change_tracker)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    try {
        PositionRestraintMgr *mgr = new PositionRestraintMgr(s, ct);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //position_restraint_mgr_new

extern "C" EXPORT void
position_restraint_mgr_delete(void *mgr)
{
    PositionRestraintMgr *m = static_cast<PositionRestraintMgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
} //position_restraint_mgr_delete

extern "C" EXPORT size_t
position_restraint_mgr_get_restraint(void *mgr, void *atom, npy_bool create, size_t n, pyobject_t *restraint)
{
    PositionRestraintMgr *m = static_cast<PositionRestraintMgr *>(mgr);
    Atom **a = static_cast<Atom **>(atom);
    size_t count = 0;
    try {
        for (size_t i=0; i<n; ++i) {
            PositionRestraint *r = m->get_restraint(*a++, create);
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
    PositionRestraintMgr *m = static_cast<PositionRestraintMgr *>(mgr);
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
    PositionRestraintMgr *m = static_cast<PositionRestraintMgr *>(mgr);
    try {
        std::vector<PositionRestraint *> visibles = m->visible_restraints();
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

/*******************************************************
 *
 * TuggableAtomsMgr functions
 *
 * NOTE: TuggableAtomsMgr is a simple subclass of
 *       PositionRestraintMgr. The Python TuggableAtom(s)
 *       classes are just re-wrappings of C++ PositionRestraint.
 *
 *******************************************************/
SET_PYTHON_INSTANCE(tuggable_atoms_mgr, TuggableAtomsMgr)
GET_PYTHON_INSTANCES(tuggable_atoms_mgr, TuggableAtomsMgr)

extern "C" EXPORT void*
tuggable_atoms_mgr_new(void *structure, void *change_tracker)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    try {
        TuggableAtomsMgr *mgr = new TuggableAtomsMgr(s, ct);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //tuggable_atoms_mgr_new

extern "C" EXPORT void
tuggable_atoms_mgr_delete(void *mgr)
{
    TuggableAtomsMgr *m = static_cast<TuggableAtomsMgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
} //tuggable_atoms_mgr_delete

extern "C" EXPORT size_t
tuggable_atoms_mgr_get_restraint(void *mgr, void *atom, npy_bool create, size_t n, pyobject_t *restraint)
{
    TuggableAtomsMgr *m = static_cast<TuggableAtomsMgr *>(mgr);
    Atom **a = static_cast<Atom **>(atom);
    size_t count = 0;
    try {
        for (size_t i=0; i<n; ++i) {
            PositionRestraint *r = m->get_restraint(*a++, create);
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
tuggable_atoms_mgr_num_restraints(void *mgr)
{
    TuggableAtomsMgr *m = static_cast<TuggableAtomsMgr *>(mgr);
    try {
        return m->num_restraints();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT PyObject*
tuggable_atoms_mgr_visible_restraints(void *mgr)
{
    TuggableAtomsMgr *m = static_cast<TuggableAtomsMgr *>(mgr);
    try {
        std::vector<PositionRestraint *> visibles = m->visible_restraints();
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
 * PositionRestraint functions
 *
 ***************************************************************/
SET_PYTHON_CLASS(position_restraint, PositionRestraint)
GET_PYTHON_INSTANCES(position_restraint, PositionRestraint)

extern "C" EXPORT void
position_restraint_get_manager(void *restraint, size_t n, pyobject_t *mgr)
{
    PositionRestraint **r = static_cast<PositionRestraint **>(restraint);
    error_wrap_array_get(r, n, &PositionRestraint::mgr, mgr);
}

extern "C" EXPORT void
position_restraint_atom(void *restraint, size_t n, pyobject_t *atom)
{
    PositionRestraint **r = static_cast<PositionRestraint **>(restraint);
    error_wrap_array_get(r, n, &PositionRestraint::atom, atom);
}

extern "C" EXPORT void
position_restraint_target(void *restraint, size_t n, double *target)
{
    PositionRestraint **r = static_cast<PositionRestraint **>(restraint);
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
    PositionRestraint **r = static_cast<PositionRestraint **>(restraint);
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
position_restraint_sim_index(void *restraint, size_t n, int *index)
{
    PositionRestraint **r = static_cast<PositionRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            *(index++) = (*r++)->get_sim_index();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_position_restraint_sim_index(void *restraint, size_t n, int *index)
{
    PositionRestraint **r = static_cast<PositionRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            (*r++)->set_sim_index(*(index++));
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
position_restraint_clear_sim_index(void *restraint, size_t n)
{
    PositionRestraint **r = static_cast<PositionRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            (*r++)->clear_sim_index();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_position_restraint_k(void *restraint, size_t n, double *k)
{
    PositionRestraint **r = static_cast<PositionRestraint **>(restraint);
    error_wrap_array_set<PositionRestraint, double, double>(r, n, &PositionRestraint::set_k, k);
}

extern "C" EXPORT void
position_restraint_k(void *restraint, size_t n, double *k)
{
    PositionRestraint **r = static_cast<PositionRestraint **>(restraint);
    error_wrap_array_get<PositionRestraint, double, double>(r, n, &PositionRestraint::get_k, k);
}

extern "C" EXPORT void
position_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    PositionRestraint **r = static_cast<PositionRestraint **>(restraint);
    error_wrap_array_get<PositionRestraint, bool, npy_bool>(r, n, &PositionRestraint::enabled, flag);
}

extern "C" EXPORT void
set_position_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    PositionRestraint **r = static_cast<PositionRestraint **>(restraint);
    error_wrap_array_set<PositionRestraint, bool, npy_bool>(r, n, &PositionRestraint::set_enabled, flag);
}

extern "C" EXPORT void
position_restraint_visible(void *restraint, size_t n, npy_bool *flag)
{
    PositionRestraint **r = static_cast<PositionRestraint **>(restraint);
    error_wrap_array_get<PositionRestraint, bool, npy_bool>(r, n, &PositionRestraint::visible, flag);
}

extern "C" EXPORT void
position_restraint_target_vector(void *restraint, size_t n, double *vec)
{
    PositionRestraint **r = static_cast<PositionRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->target_vector(vec);
            vec +=3;
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void position_restraint_bond_transform(void *restraint, size_t n, float *transform)
{
    PositionRestraint **r = static_cast<PositionRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->bond_cylinder_transform(transform);
            transform += 16;
        }
    } catch (...) {
        molc_error();
    }
}


#endif
