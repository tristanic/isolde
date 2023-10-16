/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#ifndef MDFF_EXT
#define MDFF_EXT

#include "mdff.h"
#include "../molc.h"

using namespace atomstruct;
using namespace isolde;

/*******************************************************
 *
 * MDFFMgr functions
 *
 *******************************************************/

SET_PYTHON_INSTANCE(mdff_mgr, MDFFMgr)
GET_PYTHON_INSTANCES(mdff_mgr, MDFFMgr)

extern "C" EXPORT void*
mdff_mgr_new(void *structure, void *change_tracker)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    try {
        MDFFMgr *mgr = new MDFFMgr(s, ct);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
}

extern "C" EXPORT void
mdff_mgr_delete(void *mgr)
{
    MDFFMgr *m = static_cast<MDFFMgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT size_t
mdff_mgr_get_mdff_atom(void *mgr, void *atom, size_t n, npy_bool create, pyobject_t *mdffa)
{
    MDFFMgr *m = static_cast<MDFFMgr *>(mgr);
    Atom **a = static_cast<Atom **>(atom);
    size_t count =0;
    try {
        for (size_t i=0; i<n; ++i) {
            MDFFAtom *ma = m->get_mdff_atom(*a++, create);
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
    MDFFMgr *m = static_cast<MDFFMgr *>(mgr);
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
    MDFFMgr *m = static_cast<MDFFMgr *>(mgr);
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
    MDFFMgr *m = static_cast<MDFFMgr *>(mgr);
    try {
        m->set_global_k(k);
    } catch (...) {
        molc_error();
    }
}

/*******************************************************
 *
 * MDFFAtom functions
 *
 *******************************************************/

SET_PYTHON_CLASS(mdff_atom, MDFFAtom)
GET_PYTHON_INSTANCES(mdff_atom, MDFFAtom)

extern "C" EXPORT void
mdff_atom_get_manager(void *mdffa, size_t n, pyobject_t *mgr)
{
    MDFFAtom **a = static_cast<MDFFAtom **>(mdffa);
    error_wrap_array_get(a, n, &MDFFAtom::mgr, mgr);
}

extern "C" EXPORT void
mdff_atom_atom(void *mdffa, size_t n, pyobject_t *atom)
{
    MDFFAtom **a = static_cast<MDFFAtom **>(mdffa);
    error_wrap_array_get(a, n, &MDFFAtom::atom, atom);
}

extern "C" EXPORT void
mdff_atom_sim_index(void *mdffa, size_t n, int *index)
{
    MDFFAtom **a = static_cast<MDFFAtom **>(mdffa);
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
    MDFFAtom **a = static_cast<MDFFAtom **>(mdffa);
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
    MDFFAtom **a = static_cast<MDFFAtom **>(mdffa);
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
    MDFFAtom **a = static_cast<MDFFAtom **>(mdffa);
    error_wrap_array_get(a, n, &MDFFAtom::get_coupling_constant, k);
}

extern "C" EXPORT void
set_mdff_atom_coupling_constant(void *mdffa, size_t n, double *k)
{
    MDFFAtom **a = static_cast<MDFFAtom **>(mdffa);
    error_wrap_array_set(a, n, &MDFFAtom::set_coupling_constant, k);
}

extern "C" EXPORT void
mdff_atom_enabled(void *mdffa, size_t n, npy_bool *flag)
{
    MDFFAtom **a = static_cast<MDFFAtom **>(mdffa);
    error_wrap_array_get<MDFFAtom, bool, npy_bool>(a, n, &MDFFAtom::enabled, flag);
}

extern "C" EXPORT void
set_mdff_atom_enabled(void *mdffa, size_t n, npy_bool *flag)
{
    MDFFAtom **a = static_cast<MDFFAtom **>(mdffa);
    error_wrap_array_set<MDFFAtom, bool, npy_bool>(a, n, &MDFFAtom::set_enabled, flag);
}


#endif
