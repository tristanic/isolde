#define PYINSTANCE_EXPORT
#include "rotamer_restraints.h"
#include "../util.h"
#include <pyinstance/PythonInstance.instantiate.h>


template class pyinstance::PythonInstance<isolde::Rotamer_Restraint_Mgr>;
template class pyinstance::PythonInstance<isolde::Rotamer_Restraint>;

namespace isolde
{


Rotamer_Restraint::Rotamer_Restraint( Rotamer *rotamer, Rotamer_Restraint_Mgr *mgr )
    : _rotamer(rotamer), _mgr(mgr)
{
    const auto& chi_dihedrals = rotamer->dihedrals();
    for (auto d: chi_dihedrals)
    {
        _chi_restraints.push_back(dihedral_restraint_mgr()->get_restraint(d, true));
    }
}

Proper_Dihedral_Restraint_Mgr* Rotamer_Restraint::dihedral_restraint_mgr() const
{
    return _mgr->dihedral_restraint_mgr();
}

void Rotamer_Restraint::set_spring_constant(double k) {
    for (auto r: _chi_restraints)
        r->set_spring_constant(k);
}

void Rotamer_Restraint::set_enabled(bool flag) {
    for (auto r: _chi_restraints)
        r->set_enabled(flag);
}

void Rotamer_Restraint::set_target_index(int t_index) {
    _current_target_def = rotamer()->get_target_def(t_index);
    _current_target_index = t_index;
    const auto& angles = _current_target_def->angles;
    const auto& esds = _current_target_def->esds;
    for (size_t i=0; i<n_chi(); ++i)
    {
        auto r=_chi_restraints[i];
        r->set_target(angles[i]);
        // Leave free movement within two standard deviations
        r->set_cutoff(esds[i]*2.0);
    }
}

bool Rotamer_Restraint::enabled() const {
    for (auto r: _chi_restraints) {
        if (!r->is_enabled()) return false;
    }
    return true;
}


Rotamer_Restraint_Mgr::~Rotamer_Restraint_Mgr()
{
    auto du = DestructionUser(this);
    delete_restraints(all_restraints());
}

Rotamer_Restraint* Rotamer_Restraint_Mgr::new_restraint(Rotamer *rot)
{
    auto it = _restraint_map.find(rot);
    if (it != _restraint_map.end()) {
        throw std::logic_error(error_duplicate());
        return nullptr;
    }
    if (rot->structure() != _structure) {
        throw std::logic_error(error_different_mol());
        return nullptr;
    }
    return _new_restraint(rot);
}

Rotamer_Restraint* Rotamer_Restraint_Mgr::_new_restraint(Rotamer *rot)
{
    Rotamer_Restraint *r = new Rotamer_Restraint(rot, this);
    _restraint_map[rot] = r;
    track_created(r);
    return r;
}

Rotamer_Restraint* Rotamer_Restraint_Mgr::get_restraint(Rotamer *rot, bool create)
{
    auto it = _restraint_map.find(rot);
    if (it != _restraint_map.end())
        return it->second;
    if (create)
        return new_restraint(rot);
    //throw std::logic_error(error_no_restraint());
    return nullptr;
}

std::set<Rotamer_Restraint *> Rotamer_Restraint_Mgr::all_restraints() const
{
    std::set<Rotamer_Restraint *> ret;
    for (const auto &it: _restraint_map)
        ret.insert(it.second);
    return ret;
}

void Rotamer_Restraint_Mgr::delete_restraints(const std::set<Rotamer_Restraint *>& delete_list)
{
    auto db = DestructionBatcher(this);
    _delete_restraints(delete_list);
}

void Rotamer_Restraint_Mgr::_delete_restraints(const std::set<Rotamer_Restraint *>& delete_list)
{
    for (auto d: delete_list)
    {
        auto r = d->rotamer();
        _restraint_map.erase(r);
        delete d;
    }
}

void Rotamer_Restraint_Mgr::destructors_done(const std::set<void *>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<Rotamer_Restraint *> to_delete;
    for (auto it=_restraint_map.begin(); it != _restraint_map.end();) {
        auto r = it->first;
        if (destroyed.find(static_cast<void *>(r)) != destroyed.end()) {
            to_delete.insert(it->second);
            it = _restraint_map.erase(it);
        } else {
            it++;
        }
    }
    _delete_restraints(to_delete);
}

} //namespace isolde
