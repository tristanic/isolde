/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#define PYINSTANCE_EXPORT
#include "rotamer_restraints.h"
#include "../util.h"
#include <pyinstance/PythonInstance.instantiate.h>


template class pyinstance::PythonInstance<isolde::RotamerRestraintMgr>;
template class pyinstance::PythonInstance<isolde::RotamerRestraint>;

namespace isolde
{


RotamerRestraint::RotamerRestraint( Rotamer *rotamer, RotamerRestraintMgr *mgr )
    : _rotamer(rotamer), _mgr(mgr)
{
    const auto& chi_dihedrals = rotamer->dihedrals();
    for (auto d: chi_dihedrals)
    {
        _chi_restraints.push_back(dihedral_restraint_mgr()->get_restraint(d, true));
    }
}

ProperDihedralRestraintMgr* RotamerRestraint::dihedral_restraint_mgr() const
{
    return _mgr->dihedral_restraint_mgr();
}

void RotamerRestraint::set_spring_constant(double k) {
    for (auto r: _chi_restraints)
        r->set_spring_constant(k);
}

void RotamerRestraint::set_enabled(bool flag) {
    for (auto r: _chi_restraints)
        r->set_enabled(flag);
}

void RotamerRestraint::set_target_index(int t_index) {
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

bool RotamerRestraint::enabled() const {
    for (auto r: _chi_restraints) {
        if (!r->is_enabled()) return false;
    }
    return true;
}


RotamerRestraintMgr::~RotamerRestraintMgr()
{
    auto du = DestructionUser(this);
    delete_restraints(all_restraints());
    _restraint_map.clear();
}

RotamerRestraint* RotamerRestraintMgr::new_restraint(Rotamer *rot)
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

RotamerRestraint* RotamerRestraintMgr::_new_restraint(Rotamer *rot)
{
    RotamerRestraint *r = new RotamerRestraint(rot, this);
    _restraint_map[rot] = r;
    track_created(r);
    return r;
}

RotamerRestraint* RotamerRestraintMgr::get_restraint(Rotamer *rot, bool create)
{
    auto it = _restraint_map.find(rot);
    if (it != _restraint_map.end())
        return it->second;
    if (create)
        return new_restraint(rot);
    //throw std::logic_error(error_no_restraint());
    return nullptr;
}

std::set<RotamerRestraint *> RotamerRestraintMgr::all_restraints() const
{
    std::set<RotamerRestraint *> ret;
    for (const auto &it: _restraint_map)
        ret.insert(it.second);
    return ret;
}

void RotamerRestraintMgr::delete_restraints(const std::set<RotamerRestraint *>& delete_list)
{
    auto db = DestructionBatcher(this);
    _delete_restraints(delete_list);
}

void RotamerRestraintMgr::_delete_restraints(const std::set<RotamerRestraint *>& delete_list)
{
    for (auto d: delete_list)
    {
        auto r = d->rotamer();
        _restraint_map.erase(r);
        delete d;
    }
}

void RotamerRestraintMgr::destructors_done(const std::set<void *>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<RotamerRestraint *> to_delete;
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
