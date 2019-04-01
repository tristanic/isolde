/**
 * @Author: Tristan Croll <tic20>
 * @Date:   23-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 01-Apr-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */



#define PYINSTANCE_EXPORT

#include "adaptive_distance_restraints.h"
#include <pyinstance/PythonInstance.instantiate.h>
template class pyinstance::PythonInstance<isolde::Adaptive_Distance_Restraint>;
template class pyinstance::PythonInstance<isolde::Adaptive_Distance_Restraint_Mgr>;


namespace isolde {



Adaptive_Distance_Restraint::Adaptive_Distance_Restraint(Atom *a1, Atom *a2,
    Distance_Restraint_Mgr_Tmpl<Adaptive_Distance_Restraint> *mgr)
    : _mgr(mgr)
{
    for (auto b: a1->bonds())
        for (auto a: b->atoms())
            if (a == a2)
                throw std::logic_error(err_msg_bonded());
    _atoms[0] = a1;
    _atoms[1] = a2;
}

Adaptive_Distance_Restraint::Adaptive_Distance_Restraint(Atom *a1, Atom *a2,
    Distance_Restraint_Mgr_Tmpl<Adaptive_Distance_Restraint> *mgr,
    const double &target, const double &tolerance, const double &kappa,
    const double &c, const double &alpha)
    : Adaptive_Distance_Restraint(a1, a2, mgr)
{
    set_target(target);
    set_tolerance(tolerance);
    set_kappa(kappa);
    set_c(c);
    set_alpha(alpha);
    set_enabled(false);
}

Change_Tracker* Adaptive_Distance_Restraint::change_tracker() const
{
    return _mgr->change_tracker();
}

void Adaptive_Distance_Restraint::set_target(const double &target)
{
    _target = target < MIN_DISTANCE_RESTRAINT_TARGET ? MIN_DISTANCE_RESTRAINT_TARGET : target;
    _mgr->track_change(this, change_tracker()->REASON_TARGET_CHANGED);
}

void Adaptive_Distance_Restraint::set_tolerance(const double &tolerance)
{
    _tolerance = tolerance < 0 ? 0 : tolerance;
    _mgr->track_change(this, change_tracker()->REASON_CUTOFF_CHANGED);
}

void Adaptive_Distance_Restraint::set_kappa(const double &kappa)
{
    _kappa = kappa<0 ? 0.0 : kappa;
    _mgr->track_change(this, change_tracker()->REASON_SPRING_CONSTANT_CHANGED);
}

void Adaptive_Distance_Restraint::set_alpha(const double &alpha)
{
    _alpha = alpha;
    _mgr->track_change(this, change_tracker()->REASON_ADAPTIVE_C_CHANGED);
}

void Adaptive_Distance_Restraint::set_c(const double &c)
{
    _c = c<MIN_C ? MIN_C : c;
    _mgr->track_change(this, change_tracker()->REASON_ADAPTIVE_C_CHANGED);
}

void Adaptive_Distance_Restraint::set_enabled(bool flag)
{
    if (_enabled != flag) {
        _enabled = flag;
        _mgr->track_change(this, change_tracker()->REASON_ENABLED_CHANGED);
    }
}

bool Adaptive_Distance_Restraint::visible() const
{
    return atoms()[0]->visible() && atoms()[1]->visible() && _enabled;
}

double Adaptive_Distance_Restraint::radius() const
{
    auto f = std::min(std::abs(force_magnitude()), SCALING_MAX_FORCE);
    return sqrt(f/SCALING_MAX_FORCE)
        * (LINEAR_RESTRAINT_MAX_RADIUS-LINEAR_RESTRAINT_MIN_RADIUS) + LINEAR_RESTRAINT_MIN_RADIUS;
}

void Adaptive_Distance_Restraint::target_transform(float *rot44) const
{
    float scale = get_target() / distance();
    _bond_transform(rot44, radius(), scale);
}

void Adaptive_Distance_Restraint::bond_cylinder_transform(float *rot44) const
{
    _bond_transform(rot44, 1.0, 1.0);
}

double Adaptive_Distance_Restraint::force_magnitude() const
{
    double r = distance();
    double r_m_r0 = r-_target;
    if (!_enabled || std::abs(r_m_r0) < _tolerance || _kappa == 0)
        return 0;
    double rho;
    if (r_m_r0 < -_tolerance)
        rho = _target - _tolerance;
    else
        rho = _target + _tolerance;
    double r_m_rho = r-rho;
    double c_sq = pow(_c,2);
    double r_m_rho_on_c_squared = (r_m_rho/c_sq);
    if (std::abs(_alpha - 2.0) < EPS )
        return _kappa * r_m_rho_on_c_squared;
    else if (std::abs(_alpha) < EPS )
        return _kappa * 2*r_m_rho / ( 2*c_sq + pow(r_m_rho,2));
    else {
        return
            _kappa * pow(r_m_rho*(pow(r_m_rho,2)/(c_sq*std::abs(2-_alpha)) +1), (_alpha/2-1))/c_sq;
    }
} // force_magnitude

void Adaptive_Distance_Restraint::_bond_transform(float *rot44, float radius, float length_scale) const
{
    const Coord &c0 = atoms()[0]->coord();
    const Coord &c1 = atoms()[1]->coord();
    float xyz0[3], xyz1[3];
    for (size_t i=0; i<3; ++i)
    {
        xyz0[i] = c0[i];
        xyz1[i] = c1[i];
    }
    geometry::bond_cylinder_transform_gl<float>(xyz0, xyz1, radius, length_scale, rot44);
}

template class Distance_Restraint_Mgr_Tmpl<Adaptive_Distance_Restraint>;

} //namespace isolde;
