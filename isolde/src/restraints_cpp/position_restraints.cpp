/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 02-Apr-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#define PYINSTANCE_EXPORT

#include "position_restraints.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::PositionRestraint>;
template class pyinstance::PythonInstance<isolde::PositionRestraintMgr>;
template class pyinstance::PythonInstance<isolde::TuggableAtomsMgr>;

namespace isolde
{

Change_Tracker* PositionRestraint::change_tracker() const
{
    return _mgr->change_tracker();
}

void PositionRestraint::set_target(const Real &x, const Real &y, const Real &z)
{
    _target[0]=x; _target[1]=y; _target[2]=z;
    mgr()->track_change(this, change_tracker()->REASON_TARGET_CHANGED);
}

void PositionRestraint::set_target(Real *target)
{
    for (size_t i=0; i<3; ++i)
        _target[i] = *(target++);
    mgr()->track_change(this, change_tracker()->REASON_TARGET_CHANGED);
}

void PositionRestraint::set_k(double k)
{
    _spring_constant = k<0 ? 0.0 : ( k > MAX_LINEAR_SPRING_CONSTANT ? MAX_LINEAR_SPRING_CONSTANT : k);
    mgr()->track_change(this, change_tracker()->REASON_SPRING_CONSTANT_CHANGED);
}

void PositionRestraint::set_enabled(bool flag)
{
    if (_enabled != flag)
    {
        _enabled = flag;
        mgr()->track_change(this, change_tracker()->REASON_ENABLED_CHANGED);
    }
}

void PositionRestraint::target_vector(double *vector) const
{
    for (size_t i=0; i<3; ++i)
        *vector++ = _target[i]-_atom->coord()[i];
}

double PositionRestraint::radius() const
{
    return _spring_constant/MAX_LINEAR_SPRING_CONSTANT
        * (LINEAR_RESTRAINT_MAX_RADIUS-LINEAR_RESTRAINT_MIN_RADIUS) + LINEAR_RESTRAINT_MIN_RADIUS;
}

void PositionRestraint::bond_cylinder_transform(float *rot44) const
{
    const Coord &c0 = atom()->coord();
    const Coord &c1 = get_target();
    // float xyz0[3], xyz1[3];
    // for (size_t i=0; i<3; ++i)
    // {
    //     xyz0[i] = c0[i];
    //     xyz1[i] = c1[i];
    // }
    geometry::bond_cylinder_transform_gl<Coord, float>(c0, c1, radius(), 1.0, rot44);
}

PositionRestraint* PositionRestraintMgr_Base::_new_restraint(Atom *atom, const Coord& target)
{
    if (atom->structure() != _atomic_model) {
        throw std::logic_error(error_different_mol());
    }
    // if (atom->element().number() == 1) {
    //     throw std::logic_error(error_hydrogen());
    // }
    PositionRestraint* restraint = new PositionRestraint(atom, target, this);
    _atom_to_restraint[atom] = restraint;
    track_created(restraint);
    return restraint;
}

PositionRestraint* PositionRestraintMgr_Base::_new_restraint(Atom *atom)
{
    return _new_restraint(atom, atom->coord());
}

PositionRestraint* PositionRestraintMgr_Base::get_restraint(Atom *atom, bool create)
{
    auto it = _atom_to_restraint.find(atom);
    if (it != _atom_to_restraint.end())
        return it->second;
    if (create)
        return _new_restraint(atom);
    return nullptr;
}

std::vector<PositionRestraint *> PositionRestraintMgr_Base::visible_restraints() const
{
    std::vector<PositionRestraint *> visibles;
    for (auto &it: _atom_to_restraint)
    {
        auto r = it.second;
        if (r->visible())
            visibles.push_back(r);
    }
    return visibles;
}


void PositionRestraintMgr_Base::delete_restraints(const std::set<PositionRestraint *>& to_delete)
{
    auto db = DestructionBatcher(this);
    _delete_restraints(to_delete);
}

void PositionRestraintMgr_Base::_delete_restraints(const std::set<PositionRestraint *>& to_delete)
{
    for (auto r: to_delete) {
        _atom_to_restraint.erase(r->atom());
        delete r;
    }
}

void PositionRestraintMgr_Base::destructors_done(const std::set<void *>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<PositionRestraint *> to_delete;
    for (auto &it: _atom_to_restraint) {
        auto r = it.second;
        if (destroyed.find(static_cast<void *>(r->atom())) != destroyed.end()) {
            to_delete.insert(r);
        }
    }
    _delete_restraints(to_delete);
}

PositionRestraintMgr_Base::~PositionRestraintMgr_Base()
{
    auto du = DestructionUser(this);
    for (auto &it: _atom_to_restraint) {
        delete it.second;
    }
    _atom_to_restraint.clear();
}

} //namespace isolde
