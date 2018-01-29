#define PYINSTANCE_EXPORT

#include "rota.h"
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

Rotamer::Rotamer(Residue *res, Rota_Mgr *mgr): _residue(res), _mgr(mgr)
{
    auto rname = res->name();
    _def = mgr->get_rotamer_def(rname);
    auto n_chi = _def->n_chi;
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
    std::vector<double> _angles(_def->n_chi);
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

float32_t Rotamer::score() const
{
    auto interpolator = _mgr->get_interpolator(_residue->name());
    std::vector<double> cur_angles(_def->n_chi);
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

void Rota_Mgr::add_rotamer_def(const std::string &resname, size_t n_chi, bool symmetric)
{
    if (_resname_to_rota_def.find(resname) == _resname_to_rota_def.end()) {
        Rota_Def rdef(n_chi, symmetric);
        _resname_to_rota_def[resname] = rdef;
    } else {
        throw std::runtime_error("Rotamer definition alread exists!");
    }
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
    try {
        return _residue_to_rotamer.at(residue);
    } catch (std::out_of_range) {
        return new_rotamer(residue);
    }
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
        size_t n_chi = get_rotamer_def(name)->n_chi;
        size_t n_rot = indices.size();
        size_t n_angles = n_rot*n_chi;
        std::vector<double> chi_angles(n_angles);
        for (size_t i=0; i<n_rot; i++) {
            rotamers[indices[i]]->angles(chi_angles.data()+i*n_chi);
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
