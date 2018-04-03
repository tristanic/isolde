#ifndef ISOLDE_DIHEDRAL_RESTRAINTS_EXT
#define ISOLDE_DIHEDRAL_RESTRAINTS_EXT

#include "dihedral_restraints.h"
#include "../molc.h"

using namespace atomstruct;
using namespace isolde;


/***************************************************************
 *
 * Proper_Dihedral_Restraint_Mgr functions
 *
 ***************************************************************/
SET_PYTHON_INSTANCE(proper_dihedral_restraint_mgr, Proper_Dihedral_Restraint_Mgr)
GET_PYTHON_INSTANCES(proper_dihedral_restraint_mgr, Proper_Dihedral_Restraint_Mgr)

extern "C" EXPORT void*
proper_dihedral_restraint_mgr_new(void *structure, void *change_tracker)
{
    Structure *s = static_cast<Structure *>(structure);
    isolde::Change_Tracker *ct = static_cast<isolde::Change_Tracker *>(change_tracker);
    try {
        Proper_Dihedral_Restraint_Mgr *mgr = new Proper_Dihedral_Restraint_Mgr(s, ct);
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //proper_dihedral_restraint_mgr_new

extern "C" EXPORT void
proper_dihedral_restraint_mgr_delete(void *mgr)
{
    Proper_Dihedral_Restraint_Mgr *m = static_cast<Proper_Dihedral_Restraint_Mgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
proper_dihedral_restraint_mgr_set_colors(void *mgr, uint8_t *maxc, uint8_t *midc, uint8_t *minc)
{
    Proper_Dihedral_Restraint_Mgr *m = static_cast<Proper_Dihedral_Restraint_Mgr *>(mgr);
    try {
        m->set_colors(maxc, midc, minc);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT size_t
proper_dihedral_restraint_mgr_num_restraints(void *mgr)
{
    Proper_Dihedral_Restraint_Mgr *m = static_cast<Proper_Dihedral_Restraint_Mgr *>(mgr);
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
    Proper_Dihedral_Restraint_Mgr *m = static_cast<Proper_Dihedral_Restraint_Mgr *>(mgr);
    Proper_Dihedral **d = static_cast<Proper_Dihedral **>(dihedral);
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
    Proper_Dihedral_Restraint_Mgr *m = static_cast<Proper_Dihedral_Restraint_Mgr *>(mgr);
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
    Proper_Dihedral_Restraint_Mgr *m = static_cast<Proper_Dihedral_Restraint_Mgr *>(mgr);
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        std::set<Proper_Dihedral_Restraint *> to_delete;
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
 * Proper_Dihedral_Restraint functions
 *
 ***************************************************************/

SET_PYTHON_CLASS(proper_dihedral_restraint, Proper_Dihedral_Restraint)
GET_PYTHON_INSTANCES(proper_dihedral_restraint, Proper_Dihedral_Restraint)

extern "C" EXPORT void
proper_dihedral_restraint_target(void *restraint, size_t n, double *target)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->set_cutoff(*(cutoff++));
        }
    } catch (...) {
        molc_error();
    }
}


extern "C" EXPORT void
proper_dihedral_restraint_annotation_transform(void *restraint, size_t n, double *tf1, double *tf2)
{
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        double transforms[32];
        for (size_t i=0; i<n; ++i) {
            (*r++)->get_annotation_transform(transforms);
            double *ttf1 = transforms, *ttf2 = transforms+16;
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
    Proper_Dihedral_Restraint **r = static_cast<Proper_Dihedral_Restraint **>(restraint);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->get_annotation_color(color);
            color +=4;
        }
    } catch (...) {
        molc_error();
    }
}


#endif
