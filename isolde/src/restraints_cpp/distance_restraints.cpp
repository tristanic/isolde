/**
 * @Author: Tristan Croll <tic20>
 * @Date:   23-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 02-Apr-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#define PYINSTANCE_EXPORT

#include "distance_restraints.h"
#include <pyinstance/PythonInstance.instantiate.h>
template class pyinstance::PythonInstance<isolde::DistanceRestraint>;
template class pyinstance::PythonInstance<isolde::DistanceRestraintMgr>;


namespace isolde {

DistanceRestraint::DistanceRestraint(Atom *a1, Atom *a2, DistanceRestraintMgr_Tmpl<DistanceRestraint> *mgr)
    : _mgr(mgr)
{
    for (auto b: a1->bonds())
        for (auto a: b->atoms())
            if (a == a2)
                throw std::logic_error(err_msg_bonded());
    _atoms[0] = a1;
    _atoms[1] = a2;
}

DistanceRestraint::DistanceRestraint(Atom *a1, Atom *a2, DistanceRestraintMgr_Tmpl<DistanceRestraint> *mgr,
    const double &target, const double &k)
    : DistanceRestraint(a1, a2, mgr)
{
    set_target(target);
    set_k(k);
    set_enabled(false);
}

Change_Tracker* DistanceRestraint::change_tracker() const
{
    return _mgr->change_tracker();
}

void DistanceRestraint::set_target(double target)
{
    _target = target < MIN_DISTANCE_RESTRAINT_TARGET ? MIN_DISTANCE_RESTRAINT_TARGET : target;
    _mgr->track_change(this, change_tracker()->REASON_TARGET_CHANGED);
}

void DistanceRestraint::set_k(double k)
{
    _spring_constant = k<0 ? 0.0 : ( k > MAX_LINEAR_SPRING_CONSTANT ? MAX_LINEAR_SPRING_CONSTANT : k);
    _mgr->track_change(this, change_tracker()->REASON_SPRING_CONSTANT_CHANGED);
}

void DistanceRestraint::set_enabled(bool flag)
{
    if (_enabled != flag) {
        _enabled = flag;
        _mgr->track_change(this, change_tracker()->REASON_ENABLED_CHANGED);
    }
}

bool DistanceRestraint::visible() const
{
    return atoms()[0]->visible() && atoms()[1]->visible() && _enabled;
}

double DistanceRestraint::radius() const
{
    return sqrt(_spring_constant/MAX_LINEAR_SPRING_CONSTANT)
        * (LINEAR_RESTRAINT_MAX_RADIUS-LINEAR_RESTRAINT_MIN_RADIUS) + LINEAR_RESTRAINT_MIN_RADIUS;
}

void DistanceRestraint::target_transform(float *rot44) const
{
    float scale = get_target() / distance();
    _bond_transform(rot44, radius(), scale);
}

void DistanceRestraint::bond_cylinder_transform(float *rot44) const
{
    _bond_transform(rot44, 1.0, 1.0);
}

void DistanceRestraint::_bond_transform(float *rot44, float radius, float length_scale) const
{
    const Coord &c0 = atoms()[0]->coord();
    const Coord &c1 = atoms()[1]->coord();
    // float xyz0[3], xyz1[3];
    // for (size_t i=0; i<3; ++i)
    // {
    //     xyz0[i] = c0[i];
    //     xyz1[i] = c1[i];
    // }
    geometry::bond_cylinder_transform_gl<Coord, float>(c0, c1, radius, length_scale, rot44);
}

template class DistanceRestraintMgr_Tmpl<DistanceRestraint>;

} //namespace isolde;
