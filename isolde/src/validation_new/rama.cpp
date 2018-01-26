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
        return NONE;
    }
    if (omega==nullptr || psi==nullptr)
        return NONE;
    return rama_case(omega, psi);
}
    
double Rama_Mgr::validate(Residue *residue, Proper_Dihedral_Mgr *dmgr)
{
    auto rcase = rama_case(residue, dmgr);
    if (rcase==NONE) 
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
    std::map<uint8_t, std::vector<size_t>> case_indices;
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




} //namespace isolde
