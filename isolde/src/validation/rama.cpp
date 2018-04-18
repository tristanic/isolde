/**
 * @Author: Tristan Croll
 * @Date:   03-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   Tristan Croll
 * @Last modified time: 18-Apr-2018
 * @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
 * @Copyright: Copyright 2017-2018 Tristan Croll
 */



#define PYINSTANCE_EXPORT

#include "rama.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::Rama>;
template class pyinstance::PythonInstance<isolde::Rama_Mgr>;

namespace isolde
{

const std::string Rama::OMEGA_STR  = "omega";
const std::string Rama::PHI_STR  = "phi";
const std::string Rama::PSI_STR  = "psi";
const AtomName Rama::CA_NAME = "CA";

Rama::Rama(Residue* residue, Proper_Dihedral_Mgr *dmgr, Rama_Mgr *rmgr)
    : _residue(residue), _dmgr(dmgr), _rmgr(rmgr)
{
    if (residue->polymer_type() != PT_AMINO)
        throw std::logic_error(err_msg_not_protein());
    if (no_valid_dihedrals())
        throw std::logic_error(err_msg_no_dihedrals());
    _CA_atom = _residue->find_atom(CA_NAME);
}

Dihedral* Rama::omega()
{
    _omega = (_omega != nullptr) ? _omega : _dmgr->get_dihedral(_residue, OMEGA_STR);
    return _omega;
}
Dihedral* Rama::phi()
{
    _phi = (_phi != nullptr) ? _phi : _dmgr->get_dihedral(_residue, PHI_STR);
    return _phi;
}
Dihedral* Rama::psi()
{
    _psi = (_psi != nullptr) ? _psi : _dmgr->get_dihedral(_residue, PSI_STR);
    return _psi;
}
double Rama::omega_angle()
{
    return (omega()!=nullptr) ? _omega->angle() : NONE_VAL;
}
double Rama::phi_angle()
{
    return (phi()!=nullptr) ? _phi->angle() : NONE_VAL;
}
double Rama::psi_angle()
{
    return (psi()!=nullptr) ? _psi->angle() : NONE_VAL;
}
void Rama::angles(double *angles)
{
    *angles++ = omega_angle();
    *angles++ = phi_angle();
    *angles   = psi_angle();
}
void Rama::phipsi(double *angles)
{
    *angles++ = phi_angle();
    *angles   = psi_angle();
}
bool Rama::is_valid_rama()
{
    return !(omega()==nullptr || phi()==nullptr || psi()==nullptr);
}
bool Rama::no_valid_dihedrals()
{
    return (omega()==nullptr && phi()==nullptr && psi()==nullptr);
}


double Rama::score()
{
    // if (!is_valid_rama())
    //     return NONE_VAL;
    return _rmgr->validate(this);
}

uint8_t Rama::rama_case()
{
    if (!is_valid_rama())
        return _rmgr->CASE_NONE;
    const ResName &name = _residue->name();
    if (name == "PRO") {
        if (std::abs(_omega->angle()) <= CIS_CUTOFF)
            return _rmgr->CISPRO;
        return _rmgr->TRANSPRO;
    }
    if (name == "GLY")
        return _rmgr->GLYCINE;
    if ((_psi->atoms())[3]->residue()->name() == "PRO")
        return _rmgr->PREPRO;
    if (name == "ILE"||name=="VAL")
        return _rmgr->ILEVAL;
    return _rmgr->GENERAL;
}

//! Checks a set of pointers to see if omega, phi and/or psi have been destroyed
/*! Removes any destroyed pointers, and returns true if all three have been
 *  destroyed.
 */
bool Rama::check_for_deleted_dihedrals( const std::set<void *> &destroyed)
{
    if (_omega != nullptr) {
        if (destroyed.find(static_cast<void *>(_omega)) != destroyed.end()) {
            _omega = nullptr;
        }
    }
    if (_phi != nullptr) {
        if (destroyed.find(static_cast<void *>(_phi)) != destroyed.end()) {
            _phi = nullptr;
        }
    }
    if (_psi != nullptr) {
        if (destroyed.find(static_cast<void *>(_psi)) != destroyed.end()) {
            _psi = nullptr;
        }
    }
    return (_omega==nullptr && _phi==nullptr && _psi==nullptr);
}

Rama_Mgr::~Rama_Mgr()
{
    auto du = DestructionUser(this);
    for (auto &it: _residue_to_rama) {
        delete it.second;
    }
}

Rama* Rama_Mgr::get_rama(Residue *res)
{
    auto it = _residue_to_rama.find(res);
    if (it!=_residue_to_rama.end())
        return it->second;
    Rama *r = new Rama(res, _mgr, this);
    _residue_to_rama[res] = r;
    return r;
}

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

uint8_t Rama_Mgr::rama_case(Residue *res)
{
    return get_rama(res)->rama_case();
}

double Rama_Mgr::validate(Rama *r)
{
    auto rcase = r->rama_case();
    if (rcase==CASE_NONE)
        return NO_RAMA_SCORE;
    auto &interpolator = _interpolators.at(rcase);
    double phipsi[2];
    r->phipsi(phipsi);
    return interpolator.interpolate(phipsi);
}

double Rama_Mgr::validate(Residue *r)
{
    Rama *rama = get_rama(r);
    return validate(rama);
}

void Rama_Mgr::validate(Rama **rama, size_t n, double *scores, uint8_t *r_cases)
{
    std::unordered_map<uint8_t, std::pair<std::vector<size_t>, std::vector<double> > > case_map;
    uint8_t this_case;
    double this_phipsi[2];
    for (size_t i=0; i<n; ++i)
    {
        this_case = (*rama)->rama_case();
        *r_cases++ = this_case;
        switch (this_case)
        {
            case CASE_NONE:
                scores[i] = NO_RAMA_SCORE;
                break;
            default:
                (*rama)->phipsi(this_phipsi);
                auto &case_pair = case_map[this_case];
                case_pair.first.push_back(i);
                auto &avec = case_pair.second;
                avec.push_back(this_phipsi[0]);
                avec.push_back(this_phipsi[1]);
        }
        rama++;
    }
    for (auto &cit: case_map) {
        auto &interpolator = _interpolators.at(cit.first);
        auto &cpair = cit.second;
        auto &case_indices = cpair.first;
        auto &case_angles = cpair.second;
        size_t case_n = case_indices.size();
        std::vector<double> case_scores(case_n);
        double *adata = case_angles.data();
        try { //DELETEME
            interpolator.interpolate(adata, case_n, case_scores.data());
        } catch (std::out_of_range) {
            std::cerr << "Bad data on case " << cit.first << ":"; //DELETEME
            for (size_t i=0; i<case_n; i+=2) { //DELETEME
                std::cerr << " " << adata[i] << ", " << adata[i+1] << "; "; //DELETEME
            } std::cerr << std::endl; //DELETEME
        }
        for (size_t i=0; i<case_n; ++i)
            scores[case_indices[i]] = case_scores[i];
    }
}

void Rama_Mgr::color_by_scores(double *score, uint8_t *r_case, size_t n, uint8_t *out)
{
    colors::intcolor default_color;
    colors::color_as_intcolor(_null_color, default_color);
    colors::color this_color;
    for (size_t i=0; i<n; ++i) {
        _color_by_score(*score++, *r_case++, this_color);
        for (size_t j=0; j<4; ++j) {
            *out++ = (uint8_t)(this_color[j]*255.0);
        }
    }
} //color_by_scores

void Rama_Mgr::_color_by_score(const double &score, const uint8_t &r_case, colors::color &color)
{
    if (score < 0 || r_case == CASE_NONE)
    {
        for(size_t i=0; i<4; ++i)
            color[i]=_null_color[i];
        return;
    }
    auto cmap = get_colors(r_case);
    cmap->interpolate(log(score), color);
}

int32_t Rama_Mgr::bin_score(const double &score, uint8_t r_case)
{
    if (r_case == CASE_NONE)
        return BIN_NA;
    auto c = get_cutoffs(r_case);
    if (score >= c->allowed)
        return FAVORED;
    if (score < c->outlier) {
        if (score > 0)
            return OUTLIER;
        return BIN_NA;
    }
    return ALLOWED;
} //bin_score

void Rama_Mgr::delete_ramas(const std::set<Rama *> to_delete)
{
    auto db = DestructionBatcher(this);
    _delete_ramas(to_delete);
}
void Rama_Mgr::_delete_ramas(const std::set<Rama *> to_delete)
{
    for (auto r: to_delete) {
        _residue_to_rama.erase(r->residue());
        delete r;
    }
}

void Rama_Mgr::destructors_done(const std::set<void *>& destroyed)
{
    auto db = DestructionBatcher(this);
    if (destroyed.find(static_cast<void *>(_mgr)) != destroyed.end()) {
        //The dihedral manager is gone. Delete everything.
        delete this;
        return;
    }

    std::set<Rama *> to_delete;
    // We want to delete a Ramachandran case if its CA is gone or if all three
    // of its dihedrals are gone. Otherwise we keep it on as a partial.
    for (auto &it: _residue_to_rama) {
        auto r = it.second;
        if (destroyed.find(static_cast<void *>(r->CA_atom())) != destroyed.end()) {
            to_delete.insert(r);
            continue;
        }
        if (r->check_for_deleted_dihedrals(destroyed)) {
            to_delete.insert(r);
        }
    }
    _delete_ramas(to_delete);
}


} //namespace isolde
