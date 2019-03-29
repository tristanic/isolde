/**
 * @Author: Tristan Croll <tic20>
 * @Date:   27-Mar-2019
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 29-Mar-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2019 Tristan Croll
 */



#pragma once

#include <string>
#include "../colors.h"
#include "../constants.h"
#include "../geometry/geometry.h"
#include "distance_restraints.h"
#include "changetracker.h"
#include "sim_restraint_base.h"
#include <atomstruct/destruct.h>
#include <atomstruct/string_types.h>
#include <atomstruct/Atom.h>
#include <atomstruct/Bond.h>
#include <atomstruct/PBGroup.h>
#include <pyinstance/PythonInstance.declare.h>

using namespace atomstruct;

namespace isolde
{

template<class R> class Distance_Restraint_Mgr_Tmpl;

/*! ChimeraX-side management of an "adaptive" distance restraint.
 * A tunable restraint with a flat bottom of user-definable width. Acts like a
 * standard harmonic restraint when close to the target, but can be tuned to
 * taper off (or, for that matter, get stronger) at longer distances. This is
 * useful for (e.g.) large collections of reference-model, NMR or cross-linking
 * restraints where some subset may be uncertain or simply wrong. The functional
 * form is described in https://arxiv.org/pdf/1701.03077.pdf
  */
class Adaptive_Distance_Restraint:
    public pyinstance::PythonInstance<Adaptive_Distance_Restraint>,
    public Sim_Restraint_Base
{
public:
    typedef Atom* Atoms[2];
    Adaptive_Distance_Restraint() : Sim_Restraint_Base() {}
    ~Adaptive_Distance_Restraint() { auto du=DestructionUser(this); }
    // Construct a restraint with target and spring constant initialised to zero.
    Adaptive_Distance_Restraint(Atom *a1, Atom *a2,
            Distance_Restraint_Mgr_Tmpl<Adaptive_Distance_Restraint> *mgr);
    // Construct a restraint with specified target and parameters.
    Adaptive_Distance_Restraint(Atom *a1, Atom *a2,
            Distance_Restraint_Mgr_Tmpl<Adaptive_Distance_Restraint> *mgr,
            const double &target, const double &tolerance, const double &kappa,
            const double &c, const double &alpha);

    // Get/update parameters
    double get_target() const { return _target; }
    void set_target(const double &target);
    double get_tolerance() const { return _tolerance; }
    void set_tolerance(const double &tol);
    double get_kappa() const { return _kappa; }
    void set_kappa(const double &kappa);
    double get_c() const { return _c; }
    void set_c(const double &c);
    double get_alpha() const { return _alpha; }
    void set_alpha (const double &alpha);
    double effective_spring_constant () { return _kappa / pow(_c,2); }
    bool enabled() const { return _enabled; }
    void set_enabled(bool flag);

    // Visualisation functions
    double radius() const;
    void target_transform(float *rot44) const;
    // void tolerance_transforms(float *rot44p, float *rot44m) const;
    void bond_cylinder_transform(float *rot44) const;
    bool visible() const;

    // General monitoring
    const Atoms &atoms() const {return _atoms;}
    double distance() const {return _atoms[0]->coord().distance(_atoms[1]->coord());}
    // Current magnitude of the applied force (for scaling/colouring bonds)
    double force_magnitude() const;
    Structure* structure() const {return _atoms[0]->structure();}
    Change_Tracker *change_tracker() const;

private:
    // to avoid numerical instability, define a large epsilon when switching force
    // calculations to special functions
    double EPS = 1e-4;
    double MIN_C = 1e-2;
    double SCALING_MAX_FORCE = 500.0;
    void _bond_transform(float *rot44, float radius, float length_scale) const;
    Atoms _atoms;
    Distance_Restraint_Mgr_Tmpl<Adaptive_Distance_Restraint> *_mgr;
    double _target = 0;
    double _tolerance = 0;
    double _kappa = 0;
    double _c = MIN_C;
    double _alpha = -2;
    bool _enabled=false;
    const char* err_msg_bonded()
    { return "Can't create a distance restraint between bonded atoms!";}


}; // class Distance_Restraint

class Adaptive_Distance_Restraint_Mgr:
    public Distance_Restraint_Mgr_Tmpl<Adaptive_Distance_Restraint>,
    public pyinstance::PythonInstance<Adaptive_Distance_Restraint_Mgr>
{
public:
    Adaptive_Distance_Restraint_Mgr(Structure *structure, Change_Tracker *change_tracker)
    : Distance_Restraint_Mgr_Tmpl<Adaptive_Distance_Restraint>(
        structure, change_tracker, std::type_index(typeid(this)),
        "Adaptive_Distance_Restraint_Mgr", "Adaptive_Distance_Restraints"
    )
    {}
private:
    std::type_index _mgr_type = std::type_index(typeid(this));
};


} //namespace isolde
