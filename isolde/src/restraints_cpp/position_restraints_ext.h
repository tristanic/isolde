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

/*******************************************************
 *
 * Tuggable_Atoms_Mgr functions
 *
 * NOTE: Tuggable_Atoms_Mgr is a simple subclass of
 *       Position_Restraint_Mgr. The Python Tuggable_Atom(s)
 *       classes are just re-wrappings of C++ Position_Restraint.
 *
 *******************************************************/
SET_PYTHON_INSTANCE(tuggable_atoms_mgr, Tuggable_Atoms_Mgr)
GET_PYTHON_INSTANCES(tuggable_atoms_mgr, Tuggable_Atoms_Mgr)

extern "C" EXPORT void*
tuggable_atoms_mgr_new(void *structure, void *change_tracker)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    try {
        Tuggable_Atoms_Mgr *mgr = new Tuggable_Atoms_Mgr(s, ct);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //tuggable_atoms_mgr_new

extern "C" EXPORT void
tuggable_atoms_mgr_delete(void *mgr)
{
    Tuggable_Atoms_Mgr *m = static_cast<Tuggable_Atoms_Mgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
} //tuggable_atoms_mgr_delete

extern "C" EXPORT size_t
tuggable_atoms_mgr_get_restraint(void *mgr, void *atom, npy_bool create, size_t n, pyobject_t *restraint)
{
    Tuggable_Atoms_Mgr *m = static_cast<Tuggable_Atoms_Mgr *>(mgr);
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
tuggable_atoms_mgr_num_restraints(void *mgr)
{
    Tuggable_Atoms_Mgr *m = static_cast<Tuggable_Atoms_Mgr *>(mgr);
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
    Tuggable_Atoms_Mgr *m = static_cast<Tuggable_Atoms_Mgr *>(mgr);
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
position_restraint_sim_index(void *restraint, size_t n, int *index)
{
    Position_Restraint **r = static_cast<Position_Restraint **>(restraint);
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
    Position_Restraint **r = static_cast<Position_Restraint **>(restraint);
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
    Position_Restraint **r = static_cast<Position_Restraint **>(restraint);
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

extern "C" EXPORT void position_restraint_bond_transform(void *restraint, size_t n, float *transform)
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


#endif
