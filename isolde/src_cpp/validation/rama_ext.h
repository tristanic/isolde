/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#ifndef RAMA_EXT
#define RAMA_EXT

#include "../geometry/geometry.h"
#include "../interpolation/nd_interp.h"
#include "rama.h"
#include "../molc.h"

using namespace atomstruct;
using namespace isolde;

/**********************************************************************
 *
 * RamaMgr
 *
 **********************************************************************/
SET_PYTHON_INSTANCE(rama_mgr, RamaMgr)
GET_PYTHON_INSTANCES(rama_mgr, RamaMgr)

extern "C" EXPORT void*
rama_mgr_new(void *dmgr)
{
    ProperDihedralMgr *d = static_cast<ProperDihedralMgr *>(dmgr);
    RamaMgr *mgr = new RamaMgr(d);
    try {
        return mgr;
    } catch (...) {
        molc_error();
        return nullptr;
    }
}

extern "C" EXPORT void
rama_mgr_delete(void *mgr)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    try {
        delete m;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
set_rama_mgr_cutoffs(void *mgr, size_t r_case, double outlier, double allowed)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    try {
        m->set_cutoffs(r_case, outlier, allowed);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_mgr_cutoffs(void *mgr, size_t r_case, double* cutoffs)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    try {
        auto c = m->get_cutoffs(r_case);
        cutoffs[0] = c->outlier;
        cutoffs[1] = c->allowed;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_mgr_set_color_scale(void *mgr, uint8_t *max, uint8_t *mid, uint8_t *min, uint8_t *na)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    try {
        m->set_colors(max, mid, min, na);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_mgr_get_color_scale(void *mgr, uint8_t *max, uint8_t *mid, uint8_t *min, uint8_t *na)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    try {
        auto cmap = m->get_colors(1);
        auto &mapped_colors = cmap->mapped_colors();
        auto &na_color = m->default_color();
        for (size_t i=0; i<4; ++i) {
            *min++ = (uint8_t)(mapped_colors[0].thecolor[i] * 255.0);
            *mid++ = (uint8_t)(mapped_colors[1].thecolor[i] * 255.0);
            *max++ = (uint8_t)(mapped_colors[2].thecolor[i] * 255.0);
            *na++ = (uint8_t)(na_color[i] * 255.0);
        }
    } catch (...) {
        molc_error();
    }

}

extern "C" EXPORT void
rama_mgr_add_interpolator(void *mgr, size_t r_case, size_t dim,
    uint32_t *n, double *min, double *max, double *data)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    try {
        m->add_interpolator(r_case, dim, n, min, max, data);
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT size_t
rama_mgr_interpolator_dim(void *mgr, size_t r_case)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    try {
        auto it = m->get_interpolator(r_case);
        return it->dim();
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT void
rama_mgr_interpolator_axis_lengths(void *mgr, size_t r_case, uint32_t *ret)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    try {
        auto it = m->get_interpolator(r_case);
        for (auto l: it->length())
            *ret++ = l;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_mgr_interpolator_minmax(void *mgr, size_t r_case, double *minvals, double *maxvals)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    try {
        auto it = m->get_interpolator(r_case);
        for (auto m: it->min())
            *minvals++ = m;
        for (auto m: it->max())
            *maxvals++ = m;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_mgr_interpolator_values(void *mgr, size_t r_case, double *vals)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    try {
        auto it = m->get_interpolator(r_case);
        for (auto d: it->data())
            *vals++ = d;
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT size_t
rama_mgr_get_rama(void *mgr, void *residue, size_t n, pyobject_t *ramas)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    Residue **r = static_cast<Residue **>(residue);
    size_t found=0;
    try {
        for (size_t i=0; i<n; ++i) {
            Residue *thisr = *r++;
            if (thisr->polymer_type() != PT_AMINO)
                continue;
            try {
                Rama *ram = m->get_rama(thisr);
                ramas[found++] = ram;
            } catch (std::logic_error&) {
                continue;
            }
        }
        return found;
    } catch (...) {
        molc_error();
        return 0;
    }

}

extern "C" EXPORT void
rama_mgr_rama_case(void *mgr, void *residue, size_t n, uint8_t *rcase)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    Residue **r = static_cast<Residue **>(residue);
    try {
        for (size_t i=0; i<n; ++i) {
            *rcase++ = m->rama_case(*r++);
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_mgr_validate_by_residue(void *mgr, void *residue, size_t n, double *score, uint8_t *rcase)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    Residue **r = static_cast<Residue **>(residue);
    try {
        for (size_t i=0; i<n; ++i) {
            *score++ = m->validate(*r);
            *rcase++ = m->rama_case(*r++);
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_mgr_validate(void *mgr, void *rama, size_t n, double *score, uint8_t *rcase)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    Rama **r = static_cast<Rama **>(rama);
    try {
        m->validate(r, n, score, rcase);
    } catch (...) {
        molc_error();
    }
}

//! Provide an array of colors corresponding to Ramachandran scores
extern "C" EXPORT void
rama_mgr_validate_and_color(void *mgr, void *rama, size_t n, uint8_t *colors)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    Rama **r = static_cast<Rama **>(rama);
    std::vector<double> scores(n);
    std::vector<uint8_t> rcases(n);
    try {
        m->validate(r, n, scores.data(), rcases.data());
        m->color_by_scores(scores.data(), rcases.data(), n, colors);
    } catch (...) {
        molc_error();
    }
}

//! Provide positions of CA atoms and colors corresponding to Ramachandran scores
extern "C" EXPORT size_t
rama_mgr_ca_positions_and_colors(void *mgr, void *rama, size_t n, npy_bool hide_favored,
        double *coord, uint8_t *color, npy_bool *selected)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    Rama **r = static_cast<Rama **>(rama);
    std::vector<double> scores(n);
    std::vector<uint8_t> rcases(n);
    std::vector<Rama *> non_favored;
    try {
        m->validate(r, n, scores.data(), rcases.data());
        size_t count = 0;
        for (size_t i=0; i<n; ++i)
        {
            if (rcases[i] == m->CASE_NONE) {r++; continue; }
            auto& score = scores[i];
            if (hide_favored)
                if (score > m->get_cutoffs(rcases[i])->allowed) {
                    r++; continue;
                }
            auto ca_atom = (*r)->CA_atom();
            const auto& this_coord = ca_atom->coord();
            for (size_t j=0; j<3; ++j)
                *coord++ = this_coord[j];
            m->color_by_scores(&score, &rcases[i], 1, color);
            color += 4;
            *(selected++) = ca_atom->selected();
            r++;
            count ++;
        }
        return count;
    } catch (...) {
        molc_error();
        return 0;
    }

}

//! Directly apply colors according to Ramachandran scores to CA atoms.
extern "C" EXPORT void
rama_mgr_validate_and_color_cas(void *mgr, void *rama, size_t n, npy_bool hide_favored)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    Rama **r = static_cast<Rama **>(rama);
    std::vector<double> scores(n);
    std::vector<uint8_t> rcases(n);
    std::vector<uint8_t> colors(n*4);
    std::vector<Rama *> non_favored;
    Atom::DrawMode BALL_STYLE = static_cast<Atom::DrawMode>(1); //Atom::DrawMode::Ball;
    try {
        m->validate(r, n, scores.data(), rcases.data());
        if (hide_favored) {
            uint8_t color[4];
            for (size_t i=0; i<n; ++i) {
                if (rcases[i] == m->CASE_NONE) {r++; continue;}
                auto ca = (*r)->CA_atom();
                if (scores[i] > m->get_cutoffs(rcases[i])->allowed) {
                    // Revert to the same display style as the attached CA
                    const auto& neighbors = ca->neighbors();
                    for (auto a: neighbors) {
                        if (a->name() == "C") {
                            ca->set_color(a->color());
                            ca->set_draw_mode(a->draw_mode());
                        }
                    }
                }
                else {
                    ca->set_draw_mode(BALL_STYLE);
                    m->color_by_scores(&scores[i], &rcases[i], 1, color);
                    ca->set_color(color[0], color[1], color[2], color[3]);
                }
                r++;
            }
        } else {
            m->color_by_scores(scores.data(), rcases.data(), n, colors.data());
            auto cdata = colors.data();
            for (size_t i=0; i<n; ++i) {
                auto ca = (*r++)->CA_atom();
                ca->set_color(*cdata, *(cdata+1), *(cdata+2), *(cdata+3));
                ca->set_draw_mode(BALL_STYLE);
                cdata+=4;
            }
        }
    } catch (...) {
        molc_error();
    }
}

const int32_t T_INDICES[9] = {0,1,4,1,2,4,2,3,4};
const size_t D_LENGTH = 12;
const size_t V_LENGTH = 15;
const size_t N_LENGTH = 15;
const size_t T_LENGTH = 9;


extern "C" EXPORT size_t
rama_mgr_draw_cis_and_twisted_omegas(void *mgr, void *rama, size_t n, double *vertices,
    double *normals, int32_t *triangles, uint8_t *colors)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    Rama **r = static_cast<Rama **>(rama);
    try {
        size_t count=0;
        size_t vertex_count = 0;
        double abs_angle;
        for (size_t i=0; i<n; ++i) {
            const auto &resname = (*r)->residue()->name();
            auto omega = (*r++)->omega();

            if (omega == nullptr) continue;

            abs_angle = std::abs(omega->angle());
            if (abs_angle < TWISTED_CUTOFF) {
                const auto &atoms = omega->atoms();
                double *v = vertices;
                for (auto a: atoms) {
                    const auto &coord = a->coord();
                    for (size_t j=0; j<3; ++j) {
                        *v++ = coord[j];
                    }
                }
                for (size_t j=0; j<3; ++j) {
                    *v++ = (vertices[j] + vertices[9+j])/2;
                }
                // Pretend the surface is planar, and assign a single normal
                double v1[3], v2[3];
                for (size_t j=0; j<3; ++j) {
                    v1[j] = vertices[3+j]-vertices[j];
                    v2[j] = vertices[9+j]-vertices[j];
                }
                geometry::cross_product_3D(v1, v2, normals);
                // Copy this normal to the others
                for (size_t j=3; j<N_LENGTH; ++j) {
                    normals[j] = normals[j%3];
                }
                for (size_t j=0; j<T_LENGTH; ++j) {
                    triangles[j] = vertex_count + T_INDICES[j];
                }
                vertex_count += V_LENGTH/3;
                vertices += V_LENGTH;
                normals += N_LENGTH;
                triangles += T_LENGTH;
                count++;
                colors::intcolor the_color;
                if (abs_angle < CIS_CUTOFF) {
                    if (resname == "PRO") {
                        colors::copy_color(m->cis_pro_color(), the_color);
                    } else {
                        colors::copy_color(m->cis_nonpro_color(), the_color);
                    }
                } else {
                    colors::copy_color(m->twisted_color(), the_color);
                }
                for (size_t j=0; j<5; ++j) {
                    for (size_t k=0; k<4; ++k) {
                        *(colors++) = the_color[k];
                    }
                }
            }
        }
        return count;
    } catch (...) {
        molc_error();
        return 0;
    }
}


extern "C" EXPORT void
rama_mgr_bin_scores(void *mgr, double *score, uint8_t *r_case, size_t n, int32_t *bin)
{
    RamaMgr *m = static_cast<RamaMgr *>(mgr);
    try
    {
        for(size_t i=0; i<n; ++i) {
            *bin++ = m->bin_score(*score++, *r_case++);
        }
    } catch (...) {
        molc_error();
    }

}

/**************************************************************
 *
 * Rama functions
 *
 **************************************************************/

SET_PYTHON_CLASS(rama, Rama)
GET_PYTHON_INSTANCES(rama, Rama)

extern "C" EXPORT void
rama_ca_atom(void *rama, size_t n, pyobject_t *atom)
{
    Rama **r = static_cast<Rama **>(rama);
    try {
        for (size_t i=0; i<n; ++i) {
            *atom++ = (*r++)->CA_atom();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_residue(void *rama, size_t n, pyobject_t *residuep)
{
    Rama **r = static_cast<Rama **>(rama);
    error_wrap_array_get(r, n, &Rama::residue, residuep);
}


extern "C" EXPORT void
rama_is_valid(void *rama, size_t n, npy_bool *valid)
{
    Rama **r = static_cast<Rama **>(rama);
    try {
        for (size_t i=0; i<n; ++i) {
            *valid++ = (*r++)->is_valid_rama();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_visible(void *rama, size_t n, npy_bool *visible)
{
    Rama **r = static_cast<Rama **>(rama);
    error_wrap_array_get(r, n, &Rama::visible, visible);
}

extern "C" EXPORT void
rama_only_hidden_by_ribbon(void *rama, size_t n, npy_bool *visible)
{
    Rama **r = static_cast<Rama **>(rama);
    error_wrap_array_get(r, n, &Rama::only_hidden_by_ribbon, visible);
}


extern "C" EXPORT void
rama_score(void *rama, size_t n, double *score)
{
    Rama **r = static_cast<Rama **>(rama);
    try {
        for (size_t i=0; i<n; ++i) {
            *score++ = (*r++)->score();
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_phipsi(void *rama, size_t n, double *angle)
{
    Rama **r = static_cast<Rama **>(rama);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->phipsi(angle);
            angle +=2;
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_omegaphipsi(void *rama, size_t n, double *angles)
{
    Rama **r = static_cast<Rama **>(rama);
    try {
        for (size_t i=0; i<n; ++i) {
            (*r++)->angles(angles);
            angles +=3;
        }
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT size_t
rama_omega(void *rama, size_t n, pyobject_t *dihedralp)
{
    Rama **r = static_cast<Rama **>(rama);
    try {
        size_t found=0;
        for (size_t i=0; i<n; ++i)
        {
            Dihedral* omega = (*r++)->omega();
            if (omega != nullptr) {
                (*dihedralp++) = omega;
                found++;
            }
        }
        return found;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT size_t
rama_phi(void *rama, size_t n, pyobject_t *dihedralp)
{
    Rama **r = static_cast<Rama **>(rama);
    try {
        size_t found=0;
        for (size_t i=0; i<n; ++i)
        {
            Dihedral* phi = (*r++)->phi();
            if (phi != nullptr) {
                (*dihedralp++) = phi;
                found++;
            }
        }
        return found;
    } catch (...) {
        molc_error();
        return 0;
    }
}

extern "C" EXPORT size_t
rama_psi(void *rama, size_t n, pyobject_t *dihedralp)
{
    Rama **r = static_cast<Rama **>(rama);
    try {
        size_t found=0;
        for (size_t i=0; i<n; ++i)
        {
            Dihedral* psi = (*r++)->psi();
            if (psi != nullptr) {
                (*dihedralp++) = psi;
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
rama_case(void *rama, size_t n, uint8_t *rcase)
{
    Rama **r = static_cast<Rama **>(rama);
    try {
        for (size_t i=0; i<n; ++i)
            *(rcase++) = (*r++)->rama_case();
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void
rama_center(void *rama, size_t n, double *coords)
{
    Rama **r = static_cast<Rama **>(rama);
    try {
        for (size_t i=0; i<n; ++i) {
            auto center = (*r++)->center();
            for (size_t j=0; j<3; ++j)
                *coords++ = center[j];
        }
    } catch(...) {
        molc_error();
    }
}

#endif
