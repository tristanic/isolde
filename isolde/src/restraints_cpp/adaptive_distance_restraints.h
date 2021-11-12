/**
 * @Author: Tristan Croll <tic20>
 * @Date:   27-Mar-2019
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 09-Apr-2019
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

template<class R> class DistanceRestraintMgr_Tmpl;

/*! ChimeraX-side management of an "adaptive" distance restraint.
 * A tunable restraint with a flat bottom of user-definable width. Acts like a
 * standard harmonic restraint when close to the target, but can be tuned to
 * taper off (or, for that matter, get stronger) at longer distances. This is
 * useful for (e.g.) large collections of reference-model, NMR or cross-linking
 * restraints where some subset may be uncertain or simply wrong. The functional
 * form is described in https://arxiv.org/pdf/1701.03077.pdf
  */
class AdaptiveDistanceRestraint:
    public pyinstance::PythonInstance<AdaptiveDistanceRestraint>,
    public Sim_Restraint_Base
{
public:
    typedef Atom* Atoms[2];
    AdaptiveDistanceRestraint() : Sim_Restraint_Base() {}
    ~AdaptiveDistanceRestraint() { auto du=DestructionUser(this); }
    // Construct a restraint with target and spring constant initialised to zero.
    AdaptiveDistanceRestraint(Atom *a1, Atom *a2,
            DistanceRestraintMgr_Tmpl<AdaptiveDistanceRestraint> *mgr);
    // Construct a restraint with specified target and parameters.
    AdaptiveDistanceRestraint(Atom *a1, Atom *a2,
            DistanceRestraintMgr_Tmpl<AdaptiveDistanceRestraint> *mgr,
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
    double effective_spring_constant () const { return _kappa / pow(_c,2); }
    bool enabled() const { return _enabled; }
    void set_enabled(bool flag);

    void set_satisfied_limit(double limit) { _satisfied_limit = limit; }
    double get_satisfied_limit() const { return _satisfied_limit; }

    colors::variable_colormap* colormap() const;
    void color(uint8_t *color) const;

    // Visualisation functions
    double radius() const;
    //! Transforms for a tripartite bond representation
    void bond_transforms(float *rot44_e1, float *rot44_m, float *rot44_e2) const;
    double display_threshold() const;
    bool visible() const;
    Coord center() const { return (_atoms[0]->coord() + _atoms[1]->coord())*0.5; }

    // General monitoring
    const Atoms &atoms() const {return _atoms;}
    double distance() const {return _atoms[0]->coord().distance(_atoms[1]->coord());}
    // Current magnitude of the applied force (for scaling/colouring bonds)
    double force_magnitude() const;
    bool satisfied() const;
    Structure* structure() const {return _atoms[0]->structure();}
    Change_Tracker *change_tracker() const;
    DistanceRestraintMgr_Tmpl<AdaptiveDistanceRestraint> *mgr() const { return _mgr; }

private:
    // to avoid numerical instability, define a large epsilon when switching force
    // calculations to special functions
    double _thresholds[4];
    const double EPS = 1e-4;
    const double MIN_C = 1e-2;
    const double SCALING_MAX_FORCE = 100.0;
    Atoms _atoms;
    DistanceRestraintMgr_Tmpl<AdaptiveDistanceRestraint> *_mgr;
    double _target = 0.1;
    double _tolerance = 0;
    double _kappa = 0;
    double _c = MIN_C;
    double _alpha = -2;
    bool _enabled=false;
    const char* err_msg_bonded()
    { return "Can't create a distance restraint between bonded atoms!";}
    void _update_thresholds()
    {
        _thresholds[0] = _target - _tolerance - 5*_c;
        _thresholds[1] = _target - _tolerance;
        _thresholds[2] = _target + _tolerance;
        _thresholds[3] = _target + _tolerance + 5*_c;
    }
    double _satisfied_limit = 0.5;

}; // class DistanceRestraint

class AdaptiveDistanceRestraintMgr:
    public DistanceRestraintMgr_Tmpl<AdaptiveDistanceRestraint>,
    public pyinstance::PythonInstance<AdaptiveDistanceRestraintMgr>
{
public:
    AdaptiveDistanceRestraintMgr(Structure *structure, Change_Tracker *change_tracker)
    : DistanceRestraintMgr_Tmpl<AdaptiveDistanceRestraint>(
        structure, change_tracker, std::type_index(typeid(this)),
        "AdaptiveDistanceRestraintMgr", "AdaptiveDistanceRestraints"
    )
    {}
    inline void set_colors(uint8_t *maxc, uint8_t *midc, uint8_t *minc)
    {
        colors::color thecolors[4];
        for (size_t i=0; i<4; ++i)
        {
            thecolors[0][i] = ((double) *(minc++)) / 255.0;
            thecolors[1][i] = ((double) *(midc)) / 255.0;
            thecolors[2][i] = ((double) *(midc++)) / 255.0;
            thecolors[3][i] = ((double) *(maxc++)) / 255.0;
        }
        _colormap = colors::variable_colormap(thecolors, 4);
    }
    colors::variable_colormap* colormap() { return &_colormap; }

    double display_threshold() const { return _display_threshold; }
    void set_display_threshold(const double &t) {
        _display_threshold = t > 0 ? t : 0;
    }

private:
    std::type_index _mgr_type = std::type_index(typeid(this));
    colors::variable_colormap _colormap;
    double _display_threshold = 0;
    // A given restraint will be designated unsatisfied if (abs(distance-target)-tolerance)/c >= _satisfied_limit.
};


} //namespace isolde
