/**
 * @Author: Tristan Croll
 * @Date:   06-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   Tristan Croll
 * @Last modified time: 18-Apr-2018
 * @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
 * @Copyright: Copyright 2017-2018 Tristan Croll
 */



#define PYINSTANCE_EXPORT

#include "distance_restraints.h"
#include <pyinstance/PythonInstance.instantiate.h>
template class pyinstance::PythonInstance<isolde::Distance_Restraint>;
template class pyinstance::PythonInstance<isolde::Distance_Restraint_Mgr>;


namespace isolde {

Distance_Restraint::Distance_Restraint(Atom *a1, Atom *a2, Distance_Restraint_Mgr *mgr)
    : _mgr(mgr)
{
    for (auto b: a1->bonds())
        for (auto a: b->atoms())
            if (a == a2)
                throw std::logic_error(err_msg_bonded());
    _atoms[0] = a1;
    _atoms[1] = a2;
}

Distance_Restraint::Distance_Restraint(Atom *a1, Atom *a2, Distance_Restraint_Mgr *mgr,
    const double &target, const double &k)
    : Distance_Restraint(a1, a2, mgr)
{
    set_target(target);
    set_k(k);
    set_enabled(false);
}

Change_Tracker* Distance_Restraint::change_tracker() const
{
    return _mgr->change_tracker();
}

void Distance_Restraint::set_target(double target)
{
    _target = target < MIN_DISTANCE_RESTRAINT_TARGET ? MIN_DISTANCE_RESTRAINT_TARGET : target;
    _mgr->track_change(this, change_tracker()->REASON_TARGET_CHANGED);
}

void Distance_Restraint::set_k(double k)
{
    _spring_constant = k<0 ? 0.0 : ( k > MAX_LINEAR_SPRING_CONSTANT ? MAX_LINEAR_SPRING_CONSTANT : k);
    _mgr->track_change(this, change_tracker()->REASON_SPRING_CONSTANT_CHANGED);
}

void Distance_Restraint::set_enabled(bool flag)
{
    if (_enabled != flag) {
        _enabled = flag;
        _mgr->track_change(this, change_tracker()->REASON_ENABLED_CHANGED);
    }
}

bool Distance_Restraint::visible() const
{
    return atoms()[0]->visible() && atoms()[1]->visible() && _enabled;
}

double Distance_Restraint::radius() const
{
    return sqrt(_spring_constant/MAX_LINEAR_SPRING_CONSTANT)
        * (LINEAR_RESTRAINT_MAX_RADIUS-LINEAR_RESTRAINT_MIN_RADIUS) + LINEAR_RESTRAINT_MIN_RADIUS;
}

void Distance_Restraint::target_transform(double *rot44) const
{
    double scale = get_target() / distance();
    _bond_transform(rot44, radius(), scale);
}

void Distance_Restraint::bond_cylinder_transform(double *rot44) const
{
    _bond_transform(rot44, 1.0, 1.0);
}

void Distance_Restraint::_bond_transform(double *rot44, double radius, double length_scale) const
{
    const Coord &c0 = atoms()[0]->coord();
    const Coord &c1 = atoms()[1]->coord();
    double xyz0[3], xyz1[3];
    for (size_t i=0; i<3; ++i)
    {
        xyz0[i] = c0[i];
        xyz1[i] = c1[i];
    }
    geometry::bond_cylinder_transform_gl<double>(xyz0, xyz1, radius, length_scale, rot44);
}


Distance_Restraint* Distance_Restraint_Mgr::new_restraint(Atom *a1, Atom *a2)
{
    auto it = _atom_to_restraints.find(a1);
    if (it != _atom_to_restraints.end())
    {
        auto &dset = it->second;
        for (auto d: dset)
        {
            auto &datoms = d->atoms();
            for (auto a: datoms)
            {
                if (a==a2)
                {
                    throw std::logic_error(error_duplicate());
                    return nullptr;
                }
            }
        }
    }
    return _new_restraint(a1, a2);
}

Distance_Restraint* Distance_Restraint_Mgr::_new_restraint(Atom *a1, Atom *a2)
{
    if (a1->structure()!=_structure || a2->structure()!=_structure)
    {
        throw std::logic_error(error_different_mol());
        return nullptr;
    }
    if (a1 == a2)
    {
        throw std::logic_error(error_same_atom());
        return nullptr;
    }
    Distance_Restraint *d = new Distance_Restraint(a1, a2, this);
    _restraints.insert(d);
    _atom_to_restraints[a1].insert(d);
    _atom_to_restraints[a2].insert(d);
    track_created(d);
    return d;
}

Distance_Restraint* Distance_Restraint_Mgr::get_restraint(Atom *a1, Atom *a2, bool create)
{
    auto it = _atom_to_restraints.find(a1);
    if (it != _atom_to_restraints.end())
    {
        auto &dset = it->second;
        for (auto d: dset)
        {
            auto &datoms = d->atoms();
            for (auto a: datoms)
            {
                if (a == a2)
                    return d;
            }
        }
    }
    if (create)
        return _new_restraint(a1, a2);
    throw std::logic_error(error_no_restraint());
    return nullptr;
}

const std::set<Distance_Restraint *>& Distance_Restraint_Mgr::get_restraints(Atom *a) const
{
    auto it = _atom_to_restraints.find(a);
    if (it != _atom_to_restraints.end())
        return it->second;
    return _null_set;
}


Distance_Restraint_Mgr::~Distance_Restraint_Mgr()
{
    auto du = DestructionUser(this);
    delete_restraints(_restraints);

}

void Distance_Restraint_Mgr::delete_restraints(const std::set<Distance_Restraint *> &delete_list)
{
    auto db = DestructionBatcher(this);
    _delete_restraints(delete_list);
} //delete_restraints

void Distance_Restraint_Mgr::_delete_restraints(const std::set<Distance_Restraint *> &delete_list )
{
    for (auto d: delete_list) {
        _restraints.erase(d);
        for (auto &a: d->atoms()) {
            auto it = _atom_to_restraints.find(a);
            if (it != _atom_to_restraints.end()) {
                auto &dset = it->second;
                dset.erase(d);
                if (!dset.size())
                    _atom_to_restraints.erase(it);
            }
        }
        delete d;
    }
}


void Distance_Restraint_Mgr::destructors_done(const std::set<void *>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<Distance_Restraint *> to_delete;
    // Need to check for deleted atoms and delete their corresponding Distance_Restraints
    for (auto it=_atom_to_restraints.begin(); it != _atom_to_restraints.end();) {
        auto a = it->first;
        if (destroyed.find(static_cast<void *>(a)) != destroyed.end()) {
            auto &dset = it->second;
            for (auto d: dset) {
                to_delete.insert(d);
            }
            it = _atom_to_restraints.erase(it);
        } else {
            it++;
        }
    }
    _delete_restraints(to_delete);
} //destructors_done



} //namespace isolde;
