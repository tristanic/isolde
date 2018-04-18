/**
 * @Author: Tristan Croll
 * @Date:   03-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   Tristan Croll
 * @Last modified time: 18-Apr-2018
 * @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
 * @Copyright: Copyright 2017-2018 Tristan Croll
 */



#ifndef MDFF_EXT
#define MDFF_EXT

#include "mdff.h"
#include "../molc.h"

using namespace atomstruct;
using namespace isolde;

/*******************************************************
 *
 * MDFF_Mgr functions
 *
 *******************************************************/

SET_PYTHON_INSTANCE(mdff_mgr, MDFF_Mgr)
GET_PYTHON_INSTANCES(mdff_mgr, MDFF_Mgr)

extern "C" EXPORT void*
mdff_mgr_new(void *structure, void *change_tracker)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    try {
        MDFF_Mgr *mgr = new MDFF_Mgr(s, ct);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
}

extern "C" EXPORT void
mdff_mgr_delete(void *mgr)
{
    MDFF_Mgr *m = static_cast<MDFF_Mgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT size_t
mdff_mgr_get_mdff_atom(void *mgr, void *atom, size_t n, npy_bool create, pyobject_t *mdffa)
{
    MDFF_Mgr *m = static_cast<MDFF_Mgr *>(mgr);
    Atom **a = static_cast<Atom **>(atom);
    size_t count =0;
    try {
        for (size_t i=0; i<n; ++i) {
            MDFF_Atom *ma = m->get_mdff_atom(*a++, create);
            if (ma!=nullptr) {
                *mdffa++ = ma;
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
mdff_mgr_num_atoms(void *mgr)
{
    MDFF_Mgr *m = static_cast<MDFF_Mgr *>(mgr);
    try {
        return m->num_atoms();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT double
mdff_mgr_global_k(void *mgr)
{
    MDFF_Mgr *m = static_cast<MDFF_Mgr *>(mgr);
    try {
        return m->global_k();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT void
set_mdff_mgr_global_k(void *mgr, double k)
{
    MDFF_Mgr *m = static_cast<MDFF_Mgr *>(mgr);
    try {
        m->set_global_k(k);
    } catch (...) {
        molc_error();
    }
}

/*******************************************************
 *
 * MDFF_Atom functions
 *
 *******************************************************/

SET_PYTHON_CLASS(mdff_atom, MDFF_Atom)
GET_PYTHON_INSTANCES(mdff_atom, MDFF_Atom)


extern "C" EXPORT void
mdff_atom_atom(void *mdffa, size_t n, pyobject_t *atom)
{
    MDFF_Atom **a = static_cast<MDFF_Atom **>(mdffa);
    error_wrap_array_get(a, n, &MDFF_Atom::atom, atom);
}

extern "C" EXPORT void
mdff_atom_sim_index(void *mdffa, size_t n, int *index)
{
    MDFF_Atom **a = static_cast<MDFF_Atom **>(mdffa);
    try {
        for (size_t i=0; i<n; ++i)
            *(index++) = (*a++)->get_sim_index();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_mdff_atom_sim_index(void *mdffa, size_t n, int *index)
{
    MDFF_Atom **a = static_cast<MDFF_Atom **>(mdffa);
    try {
        for (size_t i=0; i<n; ++i)
            (*a++)->set_sim_index(*(index++));
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
mdff_atom_clear_sim_index(void *mdffa, size_t n)
{
    MDFF_Atom **a = static_cast<MDFF_Atom **>(mdffa);
    try {
        for (size_t i=0; i<n; ++i)
            (*a++)->clear_sim_index();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
mdff_atom_coupling_constant(void *mdffa, size_t n, double *k)
{
    MDFF_Atom **a = static_cast<MDFF_Atom **>(mdffa);
    error_wrap_array_get(a, n, &MDFF_Atom::get_coupling_constant, k);
}

extern "C" EXPORT void
set_mdff_atom_coupling_constant(void *mdffa, size_t n, double *k)
{
    MDFF_Atom **a = static_cast<MDFF_Atom **>(mdffa);
    error_wrap_array_set(a, n, &MDFF_Atom::set_coupling_constant, k);
}

extern "C" EXPORT void
mdff_atom_enabled(void *mdffa, size_t n, npy_bool *flag)
{
    MDFF_Atom **a = static_cast<MDFF_Atom **>(mdffa);
    error_wrap_array_get<MDFF_Atom, bool, npy_bool>(a, n, &MDFF_Atom::enabled, flag);
}

extern "C" EXPORT void
set_mdff_atom_enabled(void *mdffa, size_t n, npy_bool *flag)
{
    MDFF_Atom **a = static_cast<MDFF_Atom **>(mdffa);
    error_wrap_array_set<MDFF_Atom, bool, npy_bool>(a, n, &MDFF_Atom::set_enabled, flag);
}


#endif
