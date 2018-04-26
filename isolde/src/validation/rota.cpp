/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */



#define PYINSTANCE_EXPORT

#include "rota.h"
#include <sstream> //DELETEME
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::Rotamer>;
template class pyinstance::PythonInstance<isolde::Rota_Mgr>;

namespace isolde
{

/**********************************************************
 *
 * Rotamer
 *
 **********************************************************/

void Rota_Def::add_target(const std::string& name, double freq, double* angles, double* esds)
{
    auto target = Rota_Target(name, n_chi(), freq, angles, esds);
    _targets.push_back(target);
}

Rotamer::Rotamer(Residue *res, Rota_Mgr *mgr): _residue(res), _mgr(mgr)
{
    auto rname = res->name();
    _def = mgr->get_rotamer_def(rname);
    auto n_chi = _def->n_chi();
    auto dmgr = mgr->dihedral_mgr();
    static const std::string basename("chi");
    for (size_t i=1; i<=n_chi; ++i) {
        //~ std::string chi_name = basename+std::to_string(i);
        auto d = dmgr->get_dihedral(res, basename+std::to_string(i), true);
        if (d==nullptr) {
            std::cerr << "Missing dihedral " << basename + std::to_string(i) << + " for residue " << res->name() <<std::endl; //DELETEME
            throw std::out_of_range("Rotamer is missing a dihedral!");
        }
        _chi_dihedrals.push_back(d);
    }
}

void Rotamer::angles(std::vector<double> &a) const
{
    angles(a.data());
}

std::vector<double> Rotamer::angles() const
{
    std::vector<double> _angles(_def->n_chi());
    angles(_angles);
    return _angles;
}

void Rotamer::angles(double *a) const
{
    for (size_t i=0; i<n_chi(); ++i) {
        auto aa = _chi_dihedrals[i]->angle();
        *a++ = aa;
    }
    if (is_symmetric()) {
        a--;
        if (*a < 0) {
            *a += M_PI;
        }
    }

}

double Rotamer::score() const
{
    auto interpolator = _mgr->get_interpolator(_residue->name());
    std::vector<double> cur_angles(_def->n_chi());
    angles(cur_angles);
    return interpolator->interpolate(cur_angles.data());
}

/************************************************************
 *
 * Rota_Mgr
 *
 ************************************************************/

Rota_Mgr::~Rota_Mgr()
{
    auto du = DestructionUser(this);
    for (auto &it: _residue_to_rotamer)
        delete it.second;
}

void Rota_Mgr::add_rotamer_def(const std::string &resname, size_t n_chi, size_t val_nchi,
    bool symmetric, const std::vector<std::vector<std::string>>& moving_atom_names)
{
    if (_resname_to_rota_def.find(resname) == _resname_to_rota_def.end()) {
        Rota_Def rdef(n_chi, val_nchi, symmetric, moving_atom_names);
        _resname_to_rota_def[resname] = rdef;
    } else {
        throw std::runtime_error("Rotamer definition alread exists!");
    }
}

void Rota_Mgr::set_colors(uint8_t *max, uint8_t *mid, uint8_t *min)
{
    colors::color thecolors[3];
    for (size_t i=0; i<4; ++i)
    {
        thecolors[0][i] = ((double) *min++) / 255.0;
        thecolors[1][i] = ((double) *mid++) / 255.0;
        thecolors[2][i] = ((double) *max++) / 255.0;
    }
    auto cuts = get_cutoffs();
    double these_cutoffs[3];
    these_cutoffs[0] = cuts->log_outlier;
    these_cutoffs[1] = cuts->log_allowed;
    these_cutoffs[2] = 0;
    _colors = colors::colormap(these_cutoffs, thecolors, 3);
}

Rota_Def* Rota_Mgr::get_rotamer_def(const std::string &resname)
{
    return &(_resname_to_rota_def.at(resname));
}

Rota_Def* Rota_Mgr::get_rotamer_def(const ResName &resname)
{
    return &(_resname_to_rota_def.at(std::string(resname)));
}

void Rota_Mgr::add_interpolator(const std::string &resname, const size_t &dim,
    uint32_t *n, double *min, double *max, double *data)
{
    _interpolators[resname] = RegularGridInterpolator<double>(dim, n, min, max, data);
}

Rotamer* Rota_Mgr::new_rotamer(Residue* residue)
{
    try {
        auto r = new Rotamer(residue, this);
        _residue_to_rotamer[residue] = r;
        return r;
    } catch (...) {
        return nullptr;
    }
}

//! Fetch the rotamer for the current residue
/*!
 * If the desired rotamer is not found, an attempt will be made to
 * create it. NOTE: if the attempt fails, nullptr will be returned.
 */
Rotamer* Rota_Mgr::get_rotamer(Residue* residue)
{
    auto iter = _residue_to_rotamer.find(residue);
    if (iter != _residue_to_rotamer.end()) {
        return iter->second;
    } else {
        return new_rotamer(residue);
    }


    //~ try {
        //~ return _residue_to_rotamer.at(residue);
    //~ } catch (std::out_of_range) {
        //~ return new_rotamer(residue);
    //~ }
}

//! Fast validation of pre-defined rotamers
void Rota_Mgr::validate(Rotamer** rotamers, size_t n, double* scores)
{
    std::map<ResName, std::vector<size_t>> case_indices;
    for (size_t i=0; i<n; ++i) {
        case_indices[rotamers[i]->residue()->name()].push_back(i);
    }
    for (auto &it: case_indices) {
        std::string name = std::string(it.first);
        auto &indices = it.second;

        // This is ever-so-slightly dodgy, but it saves code and it works.
        // The problem is that Proline is a special case: while it has three
        // chi dihedrals, only the first is actually used for validation
        // (since the other two are tightly constrained due to the cyclic
        // sidechain). Rather than introduce a whole lot of extra code for this
        // one special case, we'll get all the chi angles, but overwrite the
        // extras.
        auto rdef = get_rotamer_def(name);
        size_t n_chi = rdef->n_chi();
        size_t val_nchi = rdef->val_nchi();
        size_t n_rot = indices.size();
        size_t n_angles = n_rot*n_chi;

        std::vector<double> chi_angles(n_angles);
        for (size_t i=0; i<n_rot; i++) {
            rotamers[indices[i]]->angles(chi_angles.data()+i*val_nchi);
        }

        auto &interpolator = _interpolators.at(name);
        std::vector<double> cur_scores(n_rot);
        interpolator.interpolate(chi_angles.data(), n_rot, cur_scores.data());

        for (size_t i=0; i<n_rot; ++i) {
            scores[indices[i]] = cur_scores[i];
        }
    }
}

//! Slower, but more robust validation that allows non-rotameric residues.
/*! Residues for which no valid rotamer is available will get a score of -1.
 */
void Rota_Mgr::validate(Residue** residues, size_t n, double* scores)
{
    for (size_t i=0; i<n; ++i) {
        auto res = residues[i];
        auto rot = get_rotamer(res);
        if (rot == nullptr) {
            scores[i] = -1.0;
            continue;
        }
        auto &interpolator = _interpolators.at(std::string(res->name()));
        scores[i] = interpolator.interpolate((rot->angles()));
    }
}

/**********TESTING***********/
void Rota_Mgr::_validate_from_thread(Rotamer **rotamers, size_t n, double* scores)
{
    _thread_running = true;
    _thread_done = false;
    validate(rotamers, n, scores);
    _thread_done = true;
}
void Rota_Mgr::validate_threaded(Rotamer **rotamers, size_t n, double* scores)
{
    _validation_thread = std::thread(&Rota_Mgr::_validate_from_thread, this, rotamers, n, scores);
}
/******END TESTING***********/


void Rota_Mgr::color_by_score(double *score, size_t n, uint8_t *out)
{
    colors::color this_color;
    auto cmap = get_colors();
    for (size_t i=0; i<n; ++i) {
        cmap->interpolate(log(*score++), this_color);
        for(size_t j=0; j<4; ++j) {
            *out++ = (uint8_t)(this_color[j]*255.0);
        }
    }
}

int32_t Rota_Mgr::bin_score(const double &score)
{
    if (score >= _cutoffs.allowed)
        return FAVORED;
    if (score < _cutoffs.outlier) {
        if (score>0)
            return OUTLIER;
        return BIN_NA;
    }
    return ALLOWED;
} //bin_score


void Rota_Mgr::destructors_done(const std::set<void*>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<Rotamer*> to_delete;
    bool del;
    for (auto it=_residue_to_rotamer.begin(); it!=_residue_to_rotamer.end();) {
        auto rot = it->second;
        del = false;
        // If a residue is deleted then any associated dihedrals will be
        // deleted by the Dihedral_Mgr, so we only need to worry about
        // checking for deleted dihedrals and rotamers.
        if (destroyed.find(static_cast<void *>(rot)) != destroyed.end()) {
            del=true;
        } else {
            for (auto d: rot->dihedrals()) {
                if (destroyed.find(static_cast<void *>(d)) != destroyed.end()) {
                    del = true;
                    break;
                }
            }
        }
        if (del) {
            to_delete.insert(rot);
            it = _residue_to_rotamer.erase(it);
        } else {
            ++it;
        }
    }
    for (auto r: to_delete)
        delete r;
}



} //namespace isolde
