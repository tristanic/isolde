/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#ifndef ISOLDE_ROTAMER_RESTRAINTS_EXT
#define ISOLDE_ROTAMER_RESTRAINTS_EXT

#include "rotamer_restraints.h"
#include "../molc.h"

using namespace atomstruct;
using namespace isolde;

/***************************************************************
 *
 * RotamerRestraintMgr functions
 *
 ***************************************************************/
SET_PYTHON_INSTANCE(rotamer_restraint_mgr, RotamerRestraintMgr)
GET_PYTHON_INSTANCES(rotamer_restraint_mgr, RotamerRestraintMgr)

extern "C" EXPORT void*
rotamer_restraint_mgr_new(void *structure, void *change_tracker,
    void *proper_dihedral_restraint_mgr, void *rota_mgr)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    ProperDihedralRestraintMgr *dmgr
        = static_cast<ProperDihedralRestraintMgr *>(proper_dihedral_restraint_mgr);
    RotaMgr *rmgr = static_cast<RotaMgr *>(rota_mgr);
    try {
        RotamerRestraintMgr *mgr = new RotamerRestraintMgr(s, ct, dmgr, rmgr);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //rotamer_restraint_mgr_new

extern "C" EXPORT void
rotamer_restraint_mgr_delete(void *mgr)
{
    RotamerRestraintMgr *m = static_cast<RotamerRestraintMgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
} //rotamer_restraint_mgr_delete

extern "C" EXPORT size_t
rotamer_restraint_mgr_num_restraints(void *mgr)
{
    RotamerRestraintMgr *m = static_cast<RotamerRestraintMgr *>(mgr);
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
    RotamerRestraintMgr *m = static_cast<RotamerRestraintMgr *>(mgr);
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
 * RotamerRestraint functions
 *
 ***************************************************************/

SET_PYTHON_CLASS(rotamer_restraint, RotamerRestraint)
GET_PYTHON_INSTANCES(rotamer_restraint, RotamerRestraint)

extern "C" EXPORT void
rotamer_restraint_get_manager(void *restraint, size_t n, pyobject_t *mgr)
{
    RotamerRestraint **r = static_cast<RotamerRestraint **>(restraint);
    error_wrap_array_get(r, n, &RotamerRestraint::mgr, mgr);
}

extern "C" EXPORT void
rotamer_restraint_rotamer(void *restraint, size_t n, pyobject_t *rotamer)
{
    RotamerRestraint** rr = static_cast<RotamerRestraint **>(restraint);
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
    RotamerRestraint** rr = static_cast<RotamerRestraint **>(restraint);
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
    RotamerRestraint* rr = static_cast<RotamerRestraint *>(restraint);
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
    RotamerRestraint** rr = static_cast<RotamerRestraint **>(restraint);
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
    RotamerRestraint** rr = static_cast<RotamerRestraint **>(restraint);
    error_wrap_array_get<RotamerRestraint, bool, npy_bool>(rr, n, &RotamerRestraint::enabled, enabled);
}

extern "C" EXPORT void
set_rotamer_restraint_enabled(void *restraint, size_t n, npy_bool *enabled)
{
    RotamerRestraint** rr = static_cast<RotamerRestraint **>(restraint);
    error_wrap_array_set<RotamerRestraint, bool, npy_bool>(rr, n, &RotamerRestraint::set_enabled, enabled);
}

extern "C" EXPORT void
rotamer_restraint_target_index(void *restraint, size_t n, int32_t *index)
{
    RotamerRestraint** rr = static_cast<RotamerRestraint **>(restraint);
    error_wrap_array_get<RotamerRestraint, int, int32_t>(rr, n, &RotamerRestraint::target_index, index);
}

extern "C" EXPORT void
set_rotamer_restraint_target_index(void *restraint, size_t n, int32_t *index)
{
    RotamerRestraint** rr = static_cast<RotamerRestraint **>(restraint);
    error_wrap_array_set<RotamerRestraint, int, int32_t>(rr, n, &RotamerRestraint::set_target_index, index);
}



#endif
