/**
 * @Author: Tristan Croll <tic20>
 * @Date:   23-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 09-Apr-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#define PYINSTANCE_EXPORT

#include "adaptive_distance_restraints.h"
#include <pyinstance/PythonInstance.instantiate.h>
template class pyinstance::PythonInstance<isolde::AdaptiveDistanceRestraint>;
template class pyinstance::PythonInstance<isolde::AdaptiveDistanceRestraintMgr>;


namespace isolde {



AdaptiveDistanceRestraint::AdaptiveDistanceRestraint(Atom *a1, Atom *a2,
    DistanceRestraintMgr_Tmpl<AdaptiveDistanceRestraint> *mgr)
    : _mgr(mgr)
{
    for (auto b: a1->bonds())
        for (auto a: b->atoms())
            if (a == a2)
                throw std::logic_error(err_msg_bonded());
    _atoms[0] = a1;
    _atoms[1] = a2;
}

AdaptiveDistanceRestraint::AdaptiveDistanceRestraint(Atom *a1, Atom *a2,
    DistanceRestraintMgr_Tmpl<AdaptiveDistanceRestraint> *mgr,
    const double &target, const double &tolerance, const double &kappa,
    const double &c, const double &alpha)
    : AdaptiveDistanceRestraint(a1, a2, mgr)
{
    set_target(target);
    set_tolerance(tolerance);
    set_kappa(kappa);
    set_c(c);
    set_alpha(alpha);
    set_enabled(false);
}

Change_Tracker* AdaptiveDistanceRestraint::change_tracker() const
{
    return _mgr->change_tracker();
}

void AdaptiveDistanceRestraint::set_target(const double &target)
{
    _target = target < MIN_DISTANCE_RESTRAINT_TARGET ? MIN_DISTANCE_RESTRAINT_TARGET : target;
    if (_tolerance > _target)
        set_tolerance(_target);
    _mgr->track_change(this, change_tracker()->REASON_TARGET_CHANGED);
    _update_thresholds();
}

void AdaptiveDistanceRestraint::set_tolerance(const double &tolerance)
{
    _tolerance = tolerance < 0 ? 0 : tolerance;
    _tolerance = _tolerance > _target ? _target : _tolerance;
    _mgr->track_change(this, change_tracker()->REASON_CUTOFF_CHANGED);
    _update_thresholds();
}

void AdaptiveDistanceRestraint::set_kappa(const double &kappa)
{
    _kappa = kappa<0 ? 0.0 : kappa;
    _mgr->track_change(this, change_tracker()->REASON_SPRING_CONSTANT_CHANGED);
}

void AdaptiveDistanceRestraint::set_alpha(const double &alpha)
{
    _alpha = alpha;
    _mgr->track_change(this, change_tracker()->REASON_ADAPTIVE_C_CHANGED);
}

void AdaptiveDistanceRestraint::set_c(const double &c)
{
    _c = c<MIN_C ? MIN_C : c;
    _mgr->track_change(this, change_tracker()->REASON_ADAPTIVE_C_CHANGED);
    _update_thresholds();
}

void AdaptiveDistanceRestraint::set_enabled(bool flag)
{
    if (_enabled != flag) {
        _enabled = flag;
        _mgr->track_change(this, change_tracker()->REASON_ENABLED_CHANGED);
    }
}

double AdaptiveDistanceRestraint::display_threshold() const
{
    return static_cast<AdaptiveDistanceRestraintMgr *>(_mgr)->display_threshold();
}

bool AdaptiveDistanceRestraint::visible() const
{
    bool display_threshold_trigger;
    if (display_threshold() == 0.0)
        display_threshold_trigger = true;
    else if ((std::abs(distance()-_target)-_tolerance)/_c >= display_threshold())
        display_threshold_trigger = true;
    else
        display_threshold_trigger = false;

    return atoms()[0]->visible() && atoms()[1]->visible()
        && _enabled
        && display_threshold_trigger;
}


double AdaptiveDistanceRestraint::radius() const
{
    auto f = std::min(std::abs(force_magnitude()), SCALING_MAX_FORCE);
    return sqrt(f/SCALING_MAX_FORCE)
        * (LINEAR_RESTRAINT_MAX_RADIUS-LINEAR_RESTRAINT_MIN_RADIUS) + LINEAR_RESTRAINT_MIN_RADIUS;
}


double AdaptiveDistanceRestraint::force_magnitude() const
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
    if (std::abs(_alpha) < EPS )
        return _kappa * 2*r_m_rho / ( 2*c_sq + pow(r_m_rho,2));
    return
        _kappa*r_m_rho/c_sq * pow( ( pow(r_m_rho,2) / (c_sq*std::abs(2-_alpha)) +1), (_alpha/2-1));
} // force_magnitude

bool AdaptiveDistanceRestraint::satisfied() const
{
    return (std::abs(distance()-_target)-_tolerance)/_c < _satisfied_limit;
}


colors::variable_colormap*
AdaptiveDistanceRestraint::colormap() const
{ return static_cast<AdaptiveDistanceRestraintMgr*>(_mgr)->colormap();}

void AdaptiveDistanceRestraint::color(uint8_t *color) const
{
    colors::color thecolor;
    colormap()->interpolate(distance(), _thresholds, thecolor);
    for (size_t j=0; j<4; ++j) {
        *(color++) = (uint8_t)(thecolor[j]*255.0);
    }
}

void
AdaptiveDistanceRestraint::bond_transforms(float *rot44_e1, float *rot44_m, float *rot44_e2) const
{
    auto r = radius();
    const Coord &c0 = atoms()[0]->coord();
    const Coord &c1 = atoms()[1]->coord();

    float d = distance();

    auto vec = c1-c0;
    auto end_frac = (d-get_target()-get_tolerance())/(2*d);
    auto end_offset = vec*end_frac;
    auto t1 = c0 + end_offset;
    geometry::bond_cylinder_transform_gl<Coord, float>(t1, c0, r*0.99, 1.0, rot44_e1);
    auto t2 = c1 - end_offset;
    geometry::bond_cylinder_transform_gl<Coord, float>(t1, t2, r, 1.0, rot44_m);
    geometry::bond_cylinder_transform_gl<Coord, float>(t2, c1, r*0.99, 1.0, rot44_e2);
}

template class DistanceRestraintMgr_Tmpl<AdaptiveDistanceRestraint>;

} //namespace isolde;
