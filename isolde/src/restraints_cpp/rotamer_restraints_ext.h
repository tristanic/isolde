/**
 * @Author: Tristan Croll
 * @Date:   13-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   Tristan Croll
 * @Last modified time: 18-Apr-2018
 * @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
 * @Copyright: Copyright 2017-2018 Tristan Croll
 */



#ifndef ISOLDE_ROTAMER_RESTRAINTS_EXT
#define ISOLDE_ROTAMER_RESTRAINTS_EXT

#include "rotamer_restraints.h"
#include "../molc.h"

using namespace atomstruct;
using namespace isolde;

/***************************************************************
 *
 * Rotamer_Restraint_Mgr functions
 *
 ***************************************************************/
SET_PYTHON_INSTANCE(rotamer_restraint_mgr, Rotamer_Restraint_Mgr)
GET_PYTHON_INSTANCES(rotamer_restraint_mgr, Rotamer_Restraint_Mgr)

extern "C" EXPORT void*
rotamer_restraint_mgr_new(void *structure, void *change_tracker,
    void *proper_dihedral_restraint_mgr, void *rota_mgr)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    Proper_Dihedral_Restraint_Mgr *dmgr
        = static_cast<Proper_Dihedral_Restraint_Mgr *>(proper_dihedral_restraint_mgr);
    Rota_Mgr *rmgr = static_cast<Rota_Mgr *>(rota_mgr);
    try {
        Rotamer_Restraint_Mgr *mgr = new Rotamer_Restraint_Mgr(s, ct, dmgr, rmgr);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //rotamer_restraint_mgr_new

extern "C" EXPORT void
rotamer_restraint_mgr_delete(void *mgr)
{
    Rotamer_Restraint_Mgr *m = static_cast<Rotamer_Restraint_Mgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
} //rotamer_restraint_mgr_delete

extern "C" EXPORT size_t
rotamer_restraint_mgr_num_restraints(void *mgr)
{
    Rotamer_Restraint_Mgr *m = static_cast<Rotamer_Restraint_Mgr *>(mgr);
    try {
        return m->num_restraints();
    } catch (...) {
        molc_error();
        return 0;
    }
} //rotamer_restraint_mgr_num_restraints

extern "C" EXPORT size_t
rotamer_restraint_mgr_get_restraint(void *mgr, void *rotamer, size_t n, bool create, pyobject_t *restraint)
{
    Rotamer_Restraint_Mgr *m = static_cast<Rotamer_Restraint_Mgr *>(mgr);
    Rotamer** r = static_cast<Rotamer **>(rotamer);
    try {
        size_t count=0;
        for (size_t i=0; i<n; ++i) {
            auto rr = m->get_restraint(*r++, create);
            if (rr != nullptr) {
                *restraint++ = rr;
                count++;
            }
        }
        return count;
    } catch (...) {
        molc_error();
        return 0;
    }
} //rotamer_restraint_mgr_get_restraint


/***************************************************************
 *
 * Rotamer_Restraint functions
 *
 ***************************************************************/

SET_PYTHON_CLASS(rotamer_restraint, Rotamer_Restraint)
GET_PYTHON_INSTANCES(rotamer_restraint, Rotamer_Restraint)

extern "C" EXPORT void
rotamer_restraint_rotamer(void *restraint, size_t n, pyobject_t *rotamer)
{
    Rotamer_Restraint** rr = static_cast<Rotamer_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *rotamer++ = (*rr++)->rotamer();
        }
    } catch(...) {
        molc_error();
    }
}

extern "C" EXPORT void
rotamer_restraint_residue(void *restraint, size_t n, pyobject_t *residue)
{
    Rotamer_Restraint** rr = static_cast<Rotamer_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *residue++ = (*rr++)->residue();
        }
    } catch(...) {
        molc_error();
    }
}

extern "C" EXPORT void
rotamer_restraint_chi_restraints(void *restraint, pyobject_t *chi)
{
    Rotamer_Restraint* rr = static_cast<Rotamer_Restraint *>(restraint);
    try {
        const auto& chis = rr->chi_restraints();
        for (auto c: chis) {
            *chi++ = c;
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_rotamer_restraint_spring_constant(void *restraint, size_t n, double k)
{
    Rotamer_Restraint** rr = static_cast<Rotamer_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*rr++)->set_spring_constant(k);
        }
    } catch(...) {
        molc_error();
    }
}

extern "C" EXPORT void
rotamer_restraint_enabled(void *restraint, size_t n, npy_bool *enabled)
{
    Rotamer_Restraint** rr = static_cast<Rotamer_Restraint **>(restraint);
    error_wrap_array_get<Rotamer_Restraint, bool, npy_bool>(rr, n, &Rotamer_Restraint::enabled, enabled);
}

extern "C" EXPORT void
set_rotamer_restraint_enabled(void *restraint, size_t n, npy_bool *enabled)
{
    Rotamer_Restraint** rr = static_cast<Rotamer_Restraint **>(restraint);
    error_wrap_array_set<Rotamer_Restraint, bool, npy_bool>(rr, n, &Rotamer_Restraint::set_enabled, enabled);
}

extern "C" EXPORT void
rotamer_restraint_target_index(void *restraint, size_t n, int32_t *index)
{
    Rotamer_Restraint** rr = static_cast<Rotamer_Restraint **>(restraint);
    error_wrap_array_get<Rotamer_Restraint, int, int32_t>(rr, n, &Rotamer_Restraint::target_index, index);
}

extern "C" EXPORT void
set_rotamer_restraint_target_index(void *restraint, size_t n, int32_t *index)
{
    Rotamer_Restraint** rr = static_cast<Rotamer_Restraint **>(restraint);
    error_wrap_array_set<Rotamer_Restraint, int, int32_t>(rr, n, &Rotamer_Restraint::set_target_index, index);
}



#endif
