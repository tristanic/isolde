#define PYINSTANCE_EXPORT

#include "rama.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::Rama_Mgr>;

namespace isolde 
{

const std::string Rama_Mgr::OMEGA_STR  = "omega";
const std::string Rama_Mgr::PHI_STR  = "phi";
const std::string Rama_Mgr::PSI_STR  = "psi";
const double Rama_Mgr::NONE_VAL = -1.0;


void Rama_Mgr::add_interpolator(size_t r_case, 
    const size_t &dim, uint32_t *n, double *min, double *max, double *data)
{
    _interpolators[r_case] = RegularGridInterpolator<double>(dim, n, min, max, data);
}

void Rama_Mgr::set_colors(uint8_t *max, uint8_t *mid, uint8_t *min, uint8_t *na)
{
    colors::color thecolors[3];
    
    for (size_t i=0; i<4; ++i)
    {
        thecolors[0][i] = ((double) *min++) / 255.0;
        thecolors[1][i] = ((double) *mid++) / 255.0;
        thecolors[2][i] = ((double) *max++) / 255.0;
        _null_color[i] = ((double) *na++) / 255.0;
    }
    for (size_t i=1; i<NUM_RAMA_CASES; ++i) {
        auto cuts = get_cutoffs(i);
        double these_cutoffs[3];
        these_cutoffs[0] = cuts->log_outlier;
        these_cutoffs[1] = cuts->log_allowed;
        these_cutoffs[2] = 0;
        _colors[i] = colors::colormap(these_cutoffs, thecolors, 3);
    }
}


uint8_t Rama_Mgr::rama_case(Dihedral *omega, Dihedral *psi)
{
    uint8_t thecase;
    const ResName &name = omega->atoms()[3]->residue()->name();
    if (name == "PRO") {
        double o_val = omega->angle();
        if (std::abs(o_val) <= _cis_cutoff)
            thecase = CISPRO;
        else
            thecase = TRANSPRO;
    } else if (name == "GLY") {
        thecase = GLYCINE;
    } else if (psi->atoms()[3]->residue()->name() == "PRO") {
        thecase = PREPRO;
    } else if (name == "ILE" || name == "VAL") {
        thecase = ILEVAL;
    } else {
        thecase = GENERAL;
    }
    return thecase;
}    

uint8_t Rama_Mgr::rama_case(Residue *res, Proper_Dihedral_Mgr *dmgr)
{
    Dihedral *omega, *psi;
    try {
        omega = dmgr->get_dihedral(res, OMEGA_STR, true);
        psi = dmgr->get_dihedral(res, PSI_STR, true);
    } catch (std::out_of_range) {
        return CASE_NONE;
    }
    if (omega==nullptr || psi==nullptr)
        return CASE_NONE;
    return rama_case(omega, psi);
}
    
double Rama_Mgr::validate(Residue *residue, Proper_Dihedral_Mgr *dmgr)
{
    auto rcase = rama_case(residue, dmgr);
    if (rcase==CASE_NONE) 
        return NONE_VAL;
        
    Dihedral *phi, *psi;
    phi = dmgr->get_dihedral(residue, PHI_STR, true);
    psi = dmgr->get_dihedral(residue, PSI_STR, true);
    double angles[2];
    angles[0] = phi->angle();
    angles[1] = psi->angle();
    auto &interpolator = _interpolators.at(rcase);
    return interpolator.interpolate(angles);
}

void Rama_Mgr::validate(Residue **residue, Dihedral **omega, Dihedral **phi, 
              Dihedral **psi, uint8_t *r_case, const size_t &n, double *scores)
{
    std::unordered_map<uint8_t, std::vector<size_t>> case_indices;
    uint8_t this_case;
    for (size_t i=0; i<n; ++i) {
        this_case = r_case[i];
        if (this_case==CISPRO || this_case==TRANSPRO) {
            // Prolines need to be checked each time
            if (std::abs(omega[i]->angle()) <= _cis_cutoff) {
                case_indices[CISPRO].push_back(i);
            } else {
                case_indices[TRANSPRO].push_back(i);
            }
        } else {
            case_indices[this_case].push_back(i);
        }
    }
    
    std::vector<size_t> &v = case_indices[0];
    for (auto j: v)
        scores[v[j]] = NONE_VAL;
    
    for (size_t i=1; i<NUM_RAMA_CASES; ++i) {
        auto &interpolator = _interpolators.at(i);
        std::vector<size_t> &v = case_indices[i];
        size_t len = v.size();
        std::vector<double> case_scores(len);
        std::vector<double> angles(len*2);
        for (size_t j=0, k=0; j<len;++j) {
            size_t &ind = v[j];
            angles[k++] = phi[ind]->angle();
            angles[k++] = psi[ind]->angle();
        }
        interpolator.interpolate(angles.data(), len, case_scores.data());
        for (size_t j=0; j<len; ++j) {
            scores[v[j]] = case_scores[j];
        }
    }
}

void Rama_Mgr::color_by_scores(double *scores, uint8_t *r_case, const size_t &n, uint8_t *out)
{
    colors::intcolor default_color;
    colors::color_as_intcolor(_null_color, default_color);
    colors::color this_color;
    for (size_t i=0; i<n; ++i) {
        if (*scores < 0)
        {
            for(size_t j=0; j<4; ++j)
                *out++=default_color[i];
            scores++;
            continue;
        }
        
        auto cmap = get_colors(*r_case++);
        cmap->interpolate(log(*scores++), this_color);
        for (size_t j=0; j<4; ++j) {
            *out++ = (uint8_t)(this_color[j]*255);
        }
    }
} //color_by_scores

int32_t Rama_Mgr::bin_score(const double &score, uint8_t r_case)
{
    if (r_case == CASE_NONE)
        return BIN_NA;
    auto c = get_cutoffs(r_case);
    if (score >= c->allowed)
        return FAVORED;
    if (score < c->outlier)
        if (score > 0)
            return OUTLIER;
        return BIN_NA;
    return ALLOWED;
} //bin_score


} //namespace isolde
