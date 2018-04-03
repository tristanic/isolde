#ifndef ROTA_EXT
#define ROTA_EXT

#include "../geometry/geometry.h"
#include "../interpolation/nd_interp.h"
#include "rota.h"
#include "../molc.h"

using namespace atomstruct;
using namespace isolde;

/**************************************************************
 *
 * Rota_Mgr functions
 *
 **************************************************************/
SET_PYTHON_INSTANCE(rota_mgr, Rota_Mgr)
GET_PYTHON_INSTANCES(rota_mgr, Rota_Mgr)

extern "C" EXPORT void*
rota_mgr_new(void *dihedral_mgr)
{
    Proper_Dihedral_Mgr *dmgr = static_cast<Proper_Dihedral_Mgr *>(dihedral_mgr);
    Rota_Mgr *mgr = new Rota_Mgr(dmgr);
    try {
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
} //rota_mgr_new

extern "C" EXPORT void
rota_mgr_add_rotamer_def(void *mgr, pyobject_t *resname, size_t n_chi, size_t val_nchi, npy_bool symmetric)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        std::string rname(PyUnicode_AsUTF8(static_cast<PyObject *>(resname[0])));
        m->add_rotamer_def(rname, n_chi, val_nchi, (bool)symmetric);
    } catch (...) {
        molc_error();
    }
} //rota_mgr_add_rotamer_def

extern "C" EXPORT void
rota_mgr_add_target_def(void *mgr, pyobject_t *resname, size_t n,
    pyobject_t *name, double *freq, double *angles, double *esds)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        std::string rname(PyUnicode_AsUTF8(static_cast<PyObject *>(resname[0])));
        auto def = m->get_rotamer_def(rname);
        auto nchi = def->n_chi();
        for (size_t i=0; i<n; ++i)
        {
            std::string tname(PyUnicode_AsUTF8(static_cast<PyObject *>(*name++)));
            def->add_target(tname, *freq++, angles, esds);
            angles += nchi;
            esds += nchi;
        }
        def->sort_targets();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rota_mgr_add_interpolator(void *mgr, pyobject_t *resname, size_t dim,
    uint32_t *n, double *min, double *max, double*data)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        std::string rname(PyUnicode_AsUTF8(static_cast<PyObject *>(resname[0])));
        m->add_interpolator(rname, dim, n, min, max, data);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rota_mgr_set_cutoffs(void *mgr, double allowed, double outlier)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        m->set_cutoffs(allowed, outlier);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rota_mgr_get_cutoffs(void *mgr, double* cutoffs)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        auto c = m->get_cutoffs();
        cutoffs[0] = c->allowed;
        cutoffs[1] = c->outlier;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_rota_mgr_color_scale(void *mgr, uint8_t *max, uint8_t *mid, uint8_t *min)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        m->set_colors(max, mid, min);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rota_mgr_color_scale(void *mgr, uint8_t *max, uint8_t *mid, uint8_t *min)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        auto cmap = m->get_colors();
        auto &mapped_colors = cmap->mapped_colors();
        for (size_t i=0; i<4; ++i) {
            *min++ = mapped_colors[0].thecolor[i];
            *mid++ = mapped_colors[1].thecolor[i];
            *max++ = mapped_colors[2].thecolor[i];
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT PyObject*
rota_mgr_get_rotamer(void *mgr, void *residue, size_t n)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    Residue **r = static_cast<Residue **>(residue);
    std::vector<Rotamer *> found;
    try {
        for (size_t i=0; i<n; ++i) {
            Rotamer* rot = m->get_rotamer(*r++);
            if (rot != nullptr) {
                found.push_back(rot);
            }
        }
        void ** rptr;
        PyObject *ra = python_voidp_array(found.size(), &rptr);
        size_t i=0;
        for (auto rot: found)
            rptr[i++] = rot;
        return ra;
    } catch (...) {
        molc_error();
        return 0;
    }
} //rota_mgr_get_rotamer

extern "C" EXPORT void
rota_mgr_validate_rotamer(void *mgr, void *rotamer, size_t n, double *scores)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    try {
        m->validate(r, n, scores);
    } catch (...) {
        molc_error();
    }
} //rota_mgr_validate_rotamer

/*******TESTING***********/
extern "C" EXPORT void
rota_mgr_validate_rotamer_threaded(void *mgr, void *rotamer, size_t n, double *scores)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    try {
        m->validate_threaded(r, n, scores);
    } catch (...) {
        molc_error();
    }
} //rota_mgr_validate_rotamer

extern "C" EXPORT npy_bool
rota_mgr_thread_running(void *mgr) {
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        return m->thread_running();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT npy_bool
rota_mgr_thread_done(void *mgr) {
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        return m->thread_done();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT void
rota_mgr_finalize_thread(void *mgr) {
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        m->finalize_thread();
    } catch (...) {
        molc_error();
    }
}


/******END TESTING***********/


extern "C" EXPORT void
rota_mgr_validate_residue(void *mgr, void *residue, size_t n, double *scores)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    Residue **r = static_cast<Residue **>(residue);
    try {
        m->validate(r,n,scores);
    } catch (...) {
        molc_error();
    }
} //rota_mgr_validate_residue

extern "C" EXPORT size_t
rota_mgr_validate_scale_and_color_rotamers(void *mgr, void *rotamer, size_t n,
    double max_scale, npy_bool non_favored_only,
    npy_bool visible_only, pyobject_t *rot_out, double* scale, uint8_t *color_out)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    size_t ret =0;
    try {
        std::vector<Rotamer *> visibles;
        if (visible_only) {
            auto rr = r;
            for (size_t i=0; i<n; ++i) {
                if ((*rr)->visible())
                    visibles.push_back(*rr++);
                else rr++;
            }
            n = visibles.size();
            r = visibles.data();
        }
        std::vector<double> scores(n);
        auto cutoffs = m->get_cutoffs();
        const auto &log_allowed = cutoffs->log_allowed;
        const auto &allowed = cutoffs->allowed;
        auto log_range = cutoffs->log_outlier-log_allowed;
        m->validate(r,n,scores.data());
        for (auto &s: scores) {
            if (non_favored_only)
                if (s>allowed) { r++; continue; }
            *rot_out++ = *r++;
            m->color_by_score(&s, 1, color_out);
            color_out +=4;
            auto log_s = log(s);
            auto this_scale = (log_s-log_allowed)/log_range + 1;
            *scale++ = this_scale>max_scale ? max_scale : this_scale;
            ret++;
        }
        return ret;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT void
rota_mgr_color_by_score(void *mgr, double *score, size_t n, uint8_t *color)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    try {
        m->color_by_score(score, n, color);
    } catch (...) {
        molc_error();
    }
} //rota_mgr_color_by_score

extern "C" EXPORT size_t
rota_mgr_non_favored(void *mgr, void *rotamer, size_t n, pyobject_t *bad, double *scores)
{
    Rota_Mgr *m = static_cast<Rota_Mgr *>(mgr);
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    size_t found=0;
    std::vector<double> vscores(n);
    int32_t bin;
    int32_t ALLOWED=m->ALLOWED;
    int32_t OUTLIER=m->OUTLIER;
    try {
        m->validate(r, n, vscores.data());
        for(size_t i=0; i<n; ++i) {
            bin = m->bin_score(vscores[i]);
            if (bin==ALLOWED || bin==OUTLIER)
            {
                bad[found] = r[i];
                scores[found++] = vscores[i];
            }
        }
        return found;
    } catch (...) {
        molc_error();
        return 0;
    }
}


 /*************************************************************
  *
  * Rotamer functions
  *
  *************************************************************/

SET_PYTHON_CLASS(rotamer, Rotamer)
GET_PYTHON_INSTANCES(rotamer, Rotamer)
extern "C" EXPORT void
rotamer_score(void *rotamer, size_t n, double *score)
{
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    error_wrap_array_get(r, n, &Rotamer::score, score);
}

extern "C" EXPORT void
rotamer_visible(void *rotamer, size_t n, npy_bool *visible)
{
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    error_wrap_array_get(r, n, &Rotamer::visible, visible);
}

extern "C" EXPORT void
rotamer_residue(void *rotamer, size_t n, pyobject_t *residue)
{
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    error_wrap_array_get(r, n, &Rotamer::residue, residue);
}

extern "C" EXPORT void
rotamer_ca_cb_bond(void *rotamer, size_t n, pyobject_t *bond)
{
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    error_wrap_array_get(r, n, &Rotamer::ca_cb_bond, bond);
}

extern "C" EXPORT void
rotamer_num_chi(void *rotamer, size_t n, uint8_t *nchi)
{
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    error_wrap_array_get(r, n, &Rotamer::n_chi, nchi);
}

extern "C" EXPORT void
rotamer_angles(void *rotamer, double *a)
{
    Rotamer *r = static_cast<Rotamer *>(rotamer);
    try {
        r->angles(a);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rotamer_num_target_defs(void *rotamer, size_t n, uint32_t *ndefs)
{
    Rotamer **r = static_cast<Rotamer **>(rotamer);
    try {
        error_wrap_array_get(r, n, &Rotamer::num_target_defs, ndefs);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT PyObject*
rotamer_target_def(void *rotamer, size_t i)
{
    Rotamer *r = static_cast<Rotamer *>(rotamer);
    PyObject* target_data = PyDict_New();
    try {
        auto target_def = r->get_target_def(i);
        PyObject *name_key = unicode_from_string(std::string("Name"));
        PyObject *name = unicode_from_string(target_def->name);
        PyDict_SetItem(target_data, name_key, name);
        Py_DECREF(name_key); Py_DECREF(name);

        PyObject *freq_key = unicode_from_string(std::string("Frequency"));
        PyObject *freq = PyFloat_FromDouble(target_def->frequency);
        PyDict_SetItem(target_data, freq_key, freq);
        Py_DECREF(freq_key); Py_DECREF(freq);

        PyObject *angle_key = unicode_from_string(std::string("Angles"));
        double *aptr;
        const auto& angles = target_def->angles;
        PyObject *aa = python_double_array(angles.size(), &aptr);
        for (const auto &a: angles)
            *aptr++ = a;
        PyDict_SetItem(target_data, angle_key, aa);
        Py_DECREF(angle_key); Py_DECREF(aa);

        PyObject *esd_key = unicode_from_string(std::string("Esds"));
        double *eptr;
        const auto& esds = target_def->esds;
        PyObject *ea = python_double_array(esds.size(), &eptr);
        for (const auto &e: esds)
            *eptr++ = e;
        PyDict_SetItem(target_data, esd_key, ea);
        Py_DECREF(esd_key); Py_DECREF(ea);

        return target_data;
    } catch(...) {
        Py_XDECREF(target_data);
        molc_error();
        return 0;
    }
}


#endif
