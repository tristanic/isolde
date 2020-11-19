/**
 * @Author: Tristan Croll <tic20>
 * @Date:   26-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 17-Sep-2020
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#ifndef ISOLDE_DIHEDRAL_RESTRAINTS_EXT
#define ISOLDE_DIHEDRAL_RESTRAINTS_EXT

#include "dihedral_restraints.h"
#include "../molc.h"

using namespace atomstruct;
using namespace isolde;


/***************************************************************
 *
 * ChiralRestraintMgr functions
 *
 ***************************************************************/
SET_PYTHON_INSTANCE(chiral_restraint_mgr, ChiralRestraintMgr)
GET_PYTHON_INSTANCES(chiral_restraint_mgr, ChiralRestraintMgr)

extern "C" EXPORT void*
chiral_restraint_mgr_new(void *structure, void *change_tracker)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    try {
        ChiralRestraintMgr *mgr = new ChiralRestraintMgr(s, ct);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //chiral_restraint_mgr_new

extern "C" EXPORT void
chiral_restraint_mgr_delete(void *mgr)
{
    ChiralRestraintMgr *m = static_cast<ChiralRestraintMgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT size_t
chiral_restraint_mgr_num_restraints(void *mgr)
{
    ChiralRestraintMgr *m = static_cast<ChiralRestraintMgr *>(mgr);
    try {
        return m->num_restraints();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT size_t
chiral_restraint_mgr_get_restraint(void *mgr, void *chiral,
        npy_bool create, size_t n, pyobject_t *restraint)
{
    ChiralRestraintMgr *m = static_cast<ChiralRestraintMgr *>(mgr);
    ChiralCenter **c = static_cast<ChiralCenter **>(chiral);
    try {
        size_t found = 0;
        for (size_t i=0; i<n; ++i) {
            auto r = m->get_restraint(*c++, create);
            if (r!=nullptr) {
                *restraint++ = r;
                found++;
            }
        }
        return found;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT void
chiral_restraint_mgr_delete_restraint(void *mgr, void *restraint, size_t n)
{
    ChiralRestraintMgr *m = static_cast<ChiralRestraintMgr *>(mgr);
    ChiralRestraint **r = static_cast<ChiralRestraint **>(restraint);
    try {
        std::set<ChiralRestraint *> to_delete;
        for (size_t i=0; i<n; ++i) {
            to_delete.insert(*r++);
        }
        m->delete_restraints(to_delete);
    } catch (...) {
        molc_error();
    }
}


/***************************************************************
 *
 * ProperDihedralRestraintMgr functions
 *
 ***************************************************************/
SET_PYTHON_INSTANCE(proper_dihedral_restraint_mgr, ProperDihedralRestraintMgr)
GET_PYTHON_INSTANCES(proper_dihedral_restraint_mgr, ProperDihedralRestraintMgr)

extern "C" EXPORT void*
proper_dihedral_restraint_mgr_new(void *structure, void *change_tracker)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    try {
        ProperDihedralRestraintMgr *mgr = new ProperDihedralRestraintMgr(s, ct);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //proper_dihedral_restraint_mgr_new

extern "C" EXPORT void
proper_dihedral_restraint_mgr_delete(void *mgr)
{
    ProperDihedralRestraintMgr *m = static_cast<ProperDihedralRestraintMgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_mgr_set_colors(void *mgr, uint8_t *maxc, uint8_t *midc, uint8_t *minc)
{
    ProperDihedralRestraintMgr *m = static_cast<ProperDihedralRestraintMgr *>(mgr);
    try {
        m->set_colors(maxc, midc, minc);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT size_t
proper_dihedral_restraint_mgr_num_restraints(void *mgr)
{
    ProperDihedralRestraintMgr *m = static_cast<ProperDihedralRestraintMgr *>(mgr);
    try {
        return m->num_restraints();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT size_t
proper_dihedral_restraint_mgr_get_restraint(void *mgr, void *dihedral,
        npy_bool create, size_t n, pyobject_t *restraint)
{
    ProperDihedralRestraintMgr *m = static_cast<ProperDihedralRestraintMgr *>(mgr);
    ProperDihedral **d = static_cast<ProperDihedral **>(dihedral);
    try {
        size_t found = 0;
        for (size_t i=0; i<n; ++i) {
            auto r = m->get_restraint(*d++, create);
            if (r!=nullptr) {
                *restraint++ = r;
                found++;
            }
        }
        return found;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT PyObject*
proper_dihedral_restraint_mgr_visible_restraints(void *mgr)
{
    ProperDihedralRestraintMgr *m = static_cast<ProperDihedralRestraintMgr *>(mgr);
    try {
        auto vis = m->visible_restraints();
        void **rptr;
        PyObject *ra = python_voidp_array(vis.size(), &rptr);
        for (auto r: vis) {
            *rptr++ = r;
        }
        return ra;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_mgr_delete_restraint(void *mgr, void *restraint, size_t n)
{
    ProperDihedralRestraintMgr *m = static_cast<ProperDihedralRestraintMgr *>(mgr);
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        std::set<ProperDihedralRestraint *> to_delete;
        for (size_t i=0; i<n; ++i) {
            to_delete.insert(*r++);
        }
        m->delete_restraints(to_delete);
    } catch (...) {
        molc_error();
    }
}

/***************************************************************
 *
 * AdaptiveDihedralRestraintMgr functions
 *
 ***************************************************************/
SET_PYTHON_INSTANCE(adaptive_dihedral_restraint_mgr, AdaptiveDihedralRestraintMgr)
GET_PYTHON_INSTANCES(adaptive_dihedral_restraint_mgr, AdaptiveDihedralRestraintMgr)

extern "C" EXPORT void*
adaptive_dihedral_restraint_mgr_new(void *structure, void *change_tracker)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    try {
        AdaptiveDihedralRestraintMgr *mgr = new AdaptiveDihedralRestraintMgr(s, ct);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //adaptive_dihedral_restraint_mgr_new

extern "C" EXPORT void
adaptive_dihedral_restraint_mgr_delete(void *mgr)
{
    AdaptiveDihedralRestraintMgr *m = static_cast<AdaptiveDihedralRestraintMgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_dihedral_restraint_mgr_set_colors(void *mgr, uint8_t *maxc, uint8_t *midc, uint8_t *minc)
{
    AdaptiveDihedralRestraintMgr *m = static_cast<AdaptiveDihedralRestraintMgr *>(mgr);
    try {
        m->set_colors(maxc, midc, minc);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT size_t
adaptive_dihedral_restraint_mgr_num_restraints(void *mgr)
{
    AdaptiveDihedralRestraintMgr *m = static_cast<AdaptiveDihedralRestraintMgr *>(mgr);
    try {
        return m->num_restraints();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT size_t
adaptive_dihedral_restraint_mgr_get_restraint(void *mgr, void *dihedral,
        npy_bool create, size_t n, pyobject_t *restraint)
{
    AdaptiveDihedralRestraintMgr *m = static_cast<AdaptiveDihedralRestraintMgr *>(mgr);
    ProperDihedral **d = static_cast<ProperDihedral **>(dihedral);
    try {
        size_t found = 0;
        for (size_t i=0; i<n; ++i) {
            auto r = m->get_restraint(*d++, create);
            if (r!=nullptr) {
                *restraint++ = r;
                found++;
            }
        }
        return found;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT PyObject*
adaptive_dihedral_restraint_mgr_visible_restraints(void *mgr)
{
    AdaptiveDihedralRestraintMgr *m = static_cast<AdaptiveDihedralRestraintMgr *>(mgr);
    try {
        auto vis = m->visible_restraints();
        void **rptr;
        PyObject *ra = python_voidp_array(vis.size(), &rptr);
        for (auto r: vis) {
            *rptr++ = r;
        }
        return ra;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT void
adaptive_dihedral_restraint_mgr_delete_restraint(void *mgr, void *restraint, size_t n)
{
    AdaptiveDihedralRestraintMgr *m = static_cast<AdaptiveDihedralRestraintMgr *>(mgr);
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        std::set<AdaptiveDihedralRestraint *> to_delete;
        for (size_t i=0; i<n; ++i) {
            to_delete.insert(*r++);
        }
        m->delete_restraints(to_delete);
    } catch (...) {
        molc_error();
    }
}




/***************************************************************
 *
 * ChiralRestraint functions
 *
 ***************************************************************/

SET_PYTHON_CLASS(chiral_restraint, ChiralRestraint)
GET_PYTHON_INSTANCES(chiral_restraint, ChiralRestraint)

extern "C" EXPORT void
chiral_restraint_get_manager(void *restraint, size_t n, pyobject_t *mgr)
{
    ChiralRestraint **r = static_cast<ChiralRestraint **>(restraint);
    error_wrap_array_get(r, n, &ChiralRestraint::mgr, mgr);
}


extern "C" EXPORT void
chiral_restraint_target(void *restraint, size_t n, double *target)
{
    ChiralRestraint **r = static_cast<ChiralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(target++) = (*r++)->get_target();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
chiral_restraint_chiral_center(void *restraint, size_t n, pyobject_t *center)
{
    ChiralRestraint **r = static_cast<ChiralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *center++ = (*r++)->get_dihedral();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
chiral_restraint_offset(void *restraint, size_t n, double *offset)
{
    ChiralRestraint **r = static_cast<ChiralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*offset++) = (*r++)->offset();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
chiral_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    ChiralRestraint **r = static_cast<ChiralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(flag++) = (*r++)->is_enabled();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_chiral_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    ChiralRestraint **r = static_cast<ChiralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_enabled(*(flag++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
chiral_restraint_k(void *restraint, size_t n, double *spring_constant)
{
    ChiralRestraint **r = static_cast<ChiralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(spring_constant++) = (*r++)->get_spring_constant();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_chiral_restraint_k(void *restraint, size_t n, double *spring_constant)
{
    ChiralRestraint **r = static_cast<ChiralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_spring_constant(*(spring_constant++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
chiral_restraint_sim_index(void *restraint, size_t n, int *index)
{
    ChiralRestraint **r = static_cast<ChiralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            *(index++) = (*r++)->get_sim_index();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_chiral_restraint_sim_index(void *restraint, size_t n, int *index)
{
    ChiralRestraint **r = static_cast<ChiralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            (*r++)->set_sim_index(*(index++));
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
chiral_restraint_clear_sim_index(void *restraint, size_t n)
{
    ChiralRestraint **r = static_cast<ChiralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            (*r++)->clear_sim_index();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
chiral_restraint_cutoff(void *restraint, size_t n, double *cutoff)
{
    ChiralRestraint **r = static_cast<ChiralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(cutoff++) = (*r++)->get_cutoff();
        }
    } catch (...) {
        molc_error();
    }
}
extern "C" EXPORT void
set_chiral_restraint_cutoff(void *restraint, size_t n, double *cutoff)
{
    ChiralRestraint **r = static_cast<ChiralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_cutoff(*(cutoff++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT PyObject*
chiral_all_atoms_in_sel(void* restraint, size_t nr, void* atom, size_t na)
{
    ChiralRestraint **r = static_cast<ChiralRestraint **>(restraint);
    Atom **a = static_cast<Atom **>(atom);
    try {
        std::set<Atom *> search_atoms;
        for (size_t i=0; i<na; ++i)
            search_atoms.insert(*(a++));
        std::vector<ChiralRestraint *> ret;
        for (size_t i=0; i<nr; ++i)
        {
            const auto& atoms = (*r)->get_dihedral()->atoms();
            bool keep=true;
            for (size_t j=0; j<4; ++j)
            {
                if (search_atoms.find(atoms[j]) == search_atoms.end())
                {
                    keep=false;
                    break;
                }
            }
            if (keep)
                ret.push_back(*(r++));
        }
        void **rptr;
        PyObject* ra = python_voidp_array(ret.size(), &rptr);
        for (auto rr: ret)
            *(rptr++) = rr;
        return ra;
    } catch (...) {
        molc_error();
        return 0;
    }
}



/***************************************************************
 *
 * ProperDihedralRestraint functions
 *
 ***************************************************************/

SET_PYTHON_CLASS(proper_dihedral_restraint, ProperDihedralRestraint)
GET_PYTHON_INSTANCES(proper_dihedral_restraint, ProperDihedralRestraint)

extern "C" EXPORT void
proper_dihedral_restraint_get_manager(void *restraint, size_t n, pyobject_t *mgr)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    error_wrap_array_get(r, n, &ProperDihedralRestraint::mgr, mgr);
}

extern "C" EXPORT void
proper_dihedral_restraint_target(void *restraint, size_t n, double *target)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(target++) = (*r++)->get_target();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_proper_dihedral_restraint_target(void *restraint, size_t n, double *target)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_target(*(target++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_dihedral(void *restraint, size_t n, pyobject_t *dihedral)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *dihedral++ = (*r++)->get_dihedral();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_offset(void *restraint, size_t n, double *offset)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*offset++) = (*r++)->offset();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(flag++) = (*r++)->is_enabled();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_proper_dihedral_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_enabled(*(flag++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_display(void *restraint, size_t n, npy_bool *flag)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(flag++) = (*r++)->get_display();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_proper_dihedral_restraint_display(void *restraint, size_t n, npy_bool *flag)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_display(*(flag++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_visible(void *restraint, size_t n, npy_bool *flag)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(flag++) = (*r++)->visible();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_k(void *restraint, size_t n, double *spring_constant)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(spring_constant++) = (*r++)->get_spring_constant();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_proper_dihedral_restraint_k(void *restraint, size_t n, double *spring_constant)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_spring_constant(*(spring_constant++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_sim_index(void *restraint, size_t n, int *index)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            *(index++) = (*r++)->get_sim_index();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_proper_dihedral_restraint_sim_index(void *restraint, size_t n, int *index)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            (*r++)->set_sim_index(*(index++));
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_clear_sim_index(void *restraint, size_t n)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            (*r++)->clear_sim_index();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_cutoff(void *restraint, size_t n, double *cutoff)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(cutoff++) = (*r++)->get_cutoff();
        }
    } catch (...) {
        molc_error();
    }
}
extern "C" EXPORT void
set_proper_dihedral_restraint_cutoff(void *restraint, size_t n, double *cutoff)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_cutoff(*(cutoff++));
        }
    } catch (...) {
        molc_error();
    }
}


extern "C" EXPORT void
proper_dihedral_restraint_annotation_transform(void *restraint, size_t n, float *tf1, float *tf2)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        float transforms[32];
        for (size_t i=0; i<n; ++i) {
            (*r++)->get_annotation_transform(transforms);
            float *ttf1 = transforms, *ttf2 = transforms+16;
            for (size_t j=0; j<16; ++j) {
                *(tf1++) = *(ttf1++);
                *(tf2++) = *(ttf2++);
            }
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_annotation_color(void *restraint, size_t n, uint8_t *color)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->get_annotation_color(color);
            color +=4;
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT PyObject*
proper_dihedral_restraint_all_atoms_in_sel(void* restraint, size_t nr, void* atom, size_t na)
{
    ProperDihedralRestraint **r = static_cast<ProperDihedralRestraint **>(restraint);
    Atom **a = static_cast<Atom **>(atom);
    try {
        std::set<Atom *> search_atoms;
        for (size_t i=0; i<na; ++i)
            search_atoms.insert(*(a++));
        std::vector<ProperDihedralRestraint *> ret;
        for (size_t i=0; i<nr; ++i)
        {
            const auto& atoms = (*r)->get_dihedral()->atoms();
            bool keep=true;
            for (size_t j=0; j<4; ++j)
            {
                if (search_atoms.find(atoms[j]) == search_atoms.end())
                {
                    keep=false;
                    break;
                }
            }
            if (keep)
                ret.push_back(*(r++));
        }
        void **rptr;
        PyObject* ra = python_voidp_array(ret.size(), &rptr);
        for (auto rr: ret)
            *(rptr++) = rr;
        return ra;
    } catch (...) {
        molc_error();
        return 0;
    }
}

/***************************************************************
 *
 * AdaptiveDihedralRestraint functions
 *
 ***************************************************************/

SET_PYTHON_CLASS(adaptive_dihedral_restraint, AdaptiveDihedralRestraint)
GET_PYTHON_INSTANCES(adaptive_dihedral_restraint, AdaptiveDihedralRestraint)

extern "C" EXPORT void
adaptive_dihedral_restraint_get_manager(void *restraint, size_t n, pyobject_t *mgr)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    error_wrap_array_get(r, n, &AdaptiveDihedralRestraint::mgr, mgr);
}

extern "C" EXPORT void
adaptive_dihedral_restraint_target(void *restraint, size_t n, double *target)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(target++) = (*r++)->get_target();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_adaptive_dihedral_restraint_target(void *restraint, size_t n, double *target)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_target(*(target++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_dihedral_restraint_dihedral(void *restraint, size_t n, pyobject_t *dihedral)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *dihedral++ = (*r++)->get_dihedral();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_dihedral_restraint_offset(void *restraint, size_t n, double *offset)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*offset++) = (*r++)->offset();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_dihedral_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(flag++) = (*r++)->is_enabled();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_adaptive_dihedral_restraint_enabled(void *restraint, size_t n, npy_bool *flag)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_enabled(*(flag++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_dihedral_restraint_display(void *restraint, size_t n, npy_bool *flag)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(flag++) = (*r++)->get_display();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_adaptive_dihedral_restraint_display(void *restraint, size_t n, npy_bool *flag)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_display(*(flag++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_dihedral_restraint_visible(void *restraint, size_t n, npy_bool *flag)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(flag++) = (*r++)->visible();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_dihedral_restraint_k(void *restraint, size_t n, double *spring_constant)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(spring_constant++) = (*r++)->get_spring_constant();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_adaptive_dihedral_restraint_k(void *restraint, size_t n, double *spring_constant)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_spring_constant(*(spring_constant++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_dihedral_restraint_sim_index(void *restraint, size_t n, int *index)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            *(index++) = (*r++)->get_sim_index();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_adaptive_dihedral_restraint_sim_index(void *restraint, size_t n, int *index)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            (*r++)->set_sim_index(*(index++));
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_dihedral_restraint_clear_sim_index(void *restraint, size_t n)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i)
            (*r++)->clear_sim_index();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_dihedral_restraint_kappa(void *restraint, size_t n, double *kappa)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(kappa++) = (*r++)->get_kappa();
        }
    } catch (...) {
        molc_error();
    }
}
extern "C" EXPORT void
set_adaptive_dihedral_restraint_kappa(void *restraint, size_t n, double *kappa)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_kappa(*(kappa++));
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_dihedral_restraint_alpha(void *restraint, size_t n, double *alpha)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            *(alpha++) = (*r++)->get_alpha();
        }
    } catch (...) {
        molc_error();
    }
}
extern "C" EXPORT void
set_adaptive_dihedral_restraint_alpha(void *restraint, size_t n, double *alpha)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_alpha(*(alpha++));
        }
    } catch (...) {
        molc_error();
    }
}


extern "C" EXPORT void
adaptive_diheddral_restraint_effective_sd(void *restraint, size_t n, double *sd)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    error_wrap_array_get(r, n, &AdaptiveDihedralRestraint::effective_sdev, sd);
}




extern "C" EXPORT void
adaptive_dihedral_restraint_annotation_transform(void *restraint, size_t n, float *tf1, float *tf2)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        float transforms[32];
        for (size_t i=0; i<n; ++i) {
            (*r++)->get_annotation_transform(transforms);
            float *ttf1 = transforms, *ttf2 = transforms+16;
            for (size_t j=0; j<16; ++j) {
                *(tf1++) = *(ttf1++);
                *(tf2++) = *(ttf2++);
            }
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
adaptive_dihedral_restraint_annotation_color(void *restraint, size_t n, uint8_t *color)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->get_annotation_color(color);
            color +=4;
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT PyObject*
adaptive_dihedral_restraint_all_atoms_in_sel(void* restraint, size_t nr, void* atom, size_t na)
{
    AdaptiveDihedralRestraint **r = static_cast<AdaptiveDihedralRestraint **>(restraint);
    Atom **a = static_cast<Atom **>(atom);
    try {
        std::set<Atom *> search_atoms;
        for (size_t i=0; i<na; ++i)
            search_atoms.insert(*(a++));
        std::vector<AdaptiveDihedralRestraint *> ret;
        for (size_t i=0; i<nr; ++i)
        {
            const auto& atoms = (*r)->get_dihedral()->atoms();
            bool keep=true;
            for (size_t j=0; j<4; ++j)
            {
                if (search_atoms.find(atoms[j]) == search_atoms.end())
                {
                    keep=false;
                    break;
                }
            }
            if (keep)
                ret.push_back(*(r++));
        }
        void **rptr;
        PyObject* ra = python_voidp_array(ret.size(), &rptr);
        for (auto rr: ret)
            *(rptr++) = rr;
        return ra;
    } catch (...) {
        molc_error();
        return 0;
    }
}


#endif
