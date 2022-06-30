/**
 * @Author: Tristan Croll <tic20>
 * @Date:   26-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 17-Sep-2020
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#ifndef ISOLDE_DIHEDRAL_RESTRAINT
#define ISOLDE_DIHEDRAL_RESTRAINT

#include <cmath>

#include "../constants.h"
#include "../util.h"
#include "../atomic_cpp/dihedral.h"
#include "../atomic_cpp/chiral.h"
#include "../colors.h"
#include "changetracker.h"
#include "sim_restraint_base.h"
#include <atomstruct/destruct.h>
#include <atomstruct/AtomicStructure.h>
#include <atomstruct/Atom.h>
#include <atomstruct/Bond.h>
#include <atomstruct/Coord.h>
#include <atomstruct/Residue.h>
#include <pyinstance/PythonInstance.declare.h>

using namespace atomstruct;
namespace isolde
{

template <class DType, class RType>
class DihedralRestraintMgr_Base;
class ProperDihedralRestraintMgr;
class ChiralRestraintMgr;
class AdaptiveDihedralRestraintMgr;

class Dihedral_Restraint_Change_Mgr
{
public:
    Dihedral_Restraint_Change_Mgr() {}
    ~Dihedral_Restraint_Change_Mgr() {}
    void set_change_tracker(Change_Tracker *ct) { _change_tracker = ct; }
    void set_colors(uint8_t *maxc, uint8_t *midc, uint8_t *minc)
    {
        colors::color thecolors[3];
        for (size_t i=0; i<4; ++i)
        {
            thecolors[0][i] = ((double) *(minc++)) / 255.0;
            thecolors[1][i] = ((double) *(midc++)) / 255.0;
            thecolors[2][i] = ((double) *(maxc++)) / 255.0;
        }
        _colormap = colors::variable_colormap(thecolors, 3);
    }
    colors::variable_colormap* colormap() { return &_colormap; }
    Change_Tracker* change_tracker() const { return _change_tracker; }
    void track_created(const void *r) { change_tracker()->add_created(_mgr_type, _mgr_pointer, r); }
    void track_change(const void *r, int reason) {change_tracker()->add_modified(_mgr_type, _mgr_pointer, r, reason);}

protected:
    std::type_index _mgr_type = std::type_index(typeid(this));
    void *_mgr_pointer = static_cast<void *>(this);
private:
    colors::variable_colormap _colormap;
    Change_Tracker *_change_tracker;

}; //class Dihedral_Restraint_Colormap



template <class DType>
class Dihedral_Restraint_Base: public Sim_Restraint_Base
{
public:
    Dihedral_Restraint_Base() : Sim_Restraint_Base() {}
    ~Dihedral_Restraint_Base() {auto du = DestructionUser(this);}
    Dihedral_Restraint_Base(DType *dihedral, Dihedral_Restraint_Change_Mgr *mgr)
        : _dihedral(dihedral), _mgr(mgr) {}
    DType* get_dihedral() const {return _dihedral;}

    // ChiralRestraint has a fixed target set by its ChiralCenter
    virtual double get_target() const {return _target;}
    /* Setters need to be implemented in the derived classes to ensure the
     * correct pointer is handed to the change tracker.
     */
    virtual void set_target(double target) {}
    virtual void set_enabled(bool flag) {}
    virtual void set_display(bool flag) {}
    virtual void set_spring_constant(const double &k) {}
    virtual void set_cutoff(double cutoff) {
        throw std::runtime_error("This class does not implement a cutoff!");
    }

    void set_satisfied_limit(double limit) { _satisfied_limit = limit; }
    double get_satisfied_limit() const { return _satisfied_limit; }
    bool satisfied() const { return std::abs(offset()) < _satisfied_limit; }

    Coord center() const { return _dihedral->center(); }


    bool get_display() const { return _display; }
    bool visible() const
    {
        return _enabled && _display && _dihedral->visible();
    }
    bool is_enabled() const { return _enabled; }
    //! Set the restraint spring constant in kJ mol-1 rad-1
    double get_spring_constant() const {return _spring_constant;}
    virtual double get_cutoff() const {
        throw std::runtime_error("This class does not implement a cutoff!");
        return 0;
    }
    //! Returns (current angle) - (target angle) in radians
    double offset() const {return util::wrapped_angle(_dihedral->angle()-get_target());}
    //! Returns (current angle) - (target angle) in degrees
    double offset_deg() const {return util::degrees(offset()); }
    //! Get the transform mapping an annotation primitive to the dihedral location
    virtual void get_annotation_transform(float *tf)
    {
        throw std::logic_error(err_msg_not_implemented());
    }
    virtual void get_annotation_color(uint8_t *color)
    {
        throw std::logic_error(err_msg_not_implemented());
    }
    Dihedral_Restraint_Change_Mgr *base_mgr() const { return _mgr; }
    isolde::Change_Tracker* change_tracker() const { return _mgr->change_tracker(); }
    colors::variable_colormap* colormap() const { return _mgr->colormap(); }

protected:
    double _cutoffs[3] {0, 1, M_PI}; // First and last cutoffs fixed, middle changes
    double _target = 0.0;
    double _spring_constant = 500.0;
    double _cutoff = 0.0;
    bool _enabled = false;
    bool _display = true;
    DType *_dihedral;
    double _satisfied_limit = 0.523598; // 30 degrees

private:
    Dihedral_Restraint_Change_Mgr *_mgr;
    const char* err_msg_not_implemented() const
        { return "Not implemented!";}

}; // Dihedral_Restraint_Base

class ChiralRestraint:
    public Dihedral_Restraint_Base<ChiralCenter>,
    public pyinstance::PythonInstance<ChiralRestraint>
{
public:
    ChiralRestraint(ChiralCenter *chiral, Dihedral_Restraint_Change_Mgr *mgr);
    double get_target() const { return _dihedral->expected_angle(); }
    void set_target(double target)
    {
        throw std::logic_error("Chiral restraint targets are immutable!");
    }
    void set_enabled(bool flag)
    {
        if (_enabled != flag)
        {
            _enabled = flag;
           base_mgr()->track_change(this, change_tracker()->REASON_ENABLED_CHANGED);
        }
    }
    void set_display(bool flag)
    {
        if (_display != flag)
        {
            _display = flag;
           base_mgr()->track_change(this, change_tracker()->REASON_DISPLAY_CHANGED);
        }
    }
    void set_spring_constant(const double &k)
    {
        _spring_constant = k<0 ? 0.0 : ( k > MAX_RADIAL_SPRING_CONSTANT ? MAX_RADIAL_SPRING_CONSTANT : k);
       base_mgr()->track_change(this, change_tracker()->REASON_SPRING_CONSTANT_CHANGED);
    }

    // Optional cutoff angle below which no force will be applied
    void set_cutoff(double cutoff)
    {
        _cutoff = cutoff; _cutoffs[1] = cutoff;
       base_mgr()->track_change(this, change_tracker()->REASON_CUTOFF_CHANGED);
    }

    double get_cutoff() const { return _cutoff; }

    // To make error_wrap_array_get automatic template substitution happy
    ChiralRestraintMgr *mgr() const;
}; // class ChiralRestraint



class ProperDihedralRestraintBase:
    public Dihedral_Restraint_Base<ProperDihedral>
{
public:
    ProperDihedralRestraintBase(ProperDihedral *dihedral, Dihedral_Restraint_Change_Mgr *mgr);
    void get_annotation_transform(float *tf);
    virtual void get_annotation_color(uint8_t *color);
    void set_target(double target)
    {
        _target=util::wrapped_angle(target);
       base_mgr()->track_change(this, change_tracker()->REASON_TARGET_CHANGED);
    }
    void set_enabled(bool flag)
    {
        if (_enabled != flag)
        {
            _enabled = flag;
           base_mgr()->track_change(this, change_tracker()->REASON_ENABLED_CHANGED);
        }
    }
    void set_display(bool flag)
    {
        if (_display != flag)
        {
            _display = flag;
           base_mgr()->track_change(this, change_tracker()->REASON_DISPLAY_CHANGED);
        }
    }
    void set_spring_constant(const double &k)
    {
        _spring_constant = k<0 ? 0.0 : ( k > MAX_RADIAL_SPRING_CONSTANT ? MAX_RADIAL_SPRING_CONSTANT : k);
       base_mgr()->track_change(this, change_tracker()->REASON_SPRING_CONSTANT_CHANGED);
    }

protected:
    using Dihedral_Restraint_Base<ProperDihedral>::_cutoff;

private:
}; // ProperDihedralRestraintBase

class ProperDihedralRestraint: public ProperDihedralRestraintBase,
    public pyinstance::PythonInstance<ProperDihedralRestraint>
{
public:
    ProperDihedralRestraint(ProperDihedral *dihedral, Dihedral_Restraint_Change_Mgr *mgr);

    // Optional cutoff angle below which no force will be applied
    void set_cutoff(double cutoff)
    {
        _cutoff = cutoff; _cutoffs[1] = cutoff;
       base_mgr()->track_change(this, change_tracker()->REASON_CUTOFF_CHANGED);
    }

    double get_cutoff() const { return _cutoff; }

    // To make error_wrap_array_get automatic template substitution happy
    ProperDihedralRestraintMgr *mgr() const;

private:

}; // class ProperDihedralRestraint



class AdaptiveDihedralRestraint: public ProperDihedralRestraintBase,
    public pyinstance::PythonInstance<AdaptiveDihedralRestraint>
{
public:
    AdaptiveDihedralRestraint(ProperDihedral *dihedral, Dihedral_Restraint_Change_Mgr *mgr);
    // Cutoff is not used in this class
    void set_kappa(double kappa) {
        _kappa = kappa < KAPPA_MIN ? 0 : kappa;
        _recalculate_cutoffs();
        base_mgr()->track_change(this, change_tracker()->REASON_ADAPTIVE_C_CHANGED);
    }
    double get_kappa() const { return _kappa; }
    void set_alpha(double alpha) {
        alpha = alpha < ALPHA_MIN ? ALPHA_MIN : alpha;
        alpha = alpha > ALPHA_MAX ? ALPHA_MAX : alpha;
        _alpha = alpha;
        base_mgr()->track_change(this, change_tracker()->REASON_ADAPTIVE_C_CHANGED);
    }
    double get_alpha() const { return _alpha; }

    double effective_sdev() const { return atan(sqrt(sqrt(4*_kappa*_kappa+1)-2*_kappa)); }

    double applied_moment() const;


    AdaptiveDihedralRestraintMgr *mgr() const;
private:
    const double DEFAULT_KAPPA=14.59; // gives a standard deviation of about 15 degrees
    const double DEFAULT_FMAX_ANGLE = 0.260; // delta-theta giving peak force for kappa=14.59
    const double KAPPA_MIN = 1e-3; // numerical instability develops if kappa is close to but not exactly zero
    const double DEFAULT_ALPHA=0.0; // Default flat outside well
    const double ALPHA_MIN=-1;
    const double ALPHA_MAX=2;
    double _kappa;
    double _alpha;
    void _recalculate_cutoffs()
    {
        auto fmax_dtheta = effective_sdev();
        _cutoffs[1] = fmax_dtheta;
        _cutoffs[2] = 2*fmax_dtheta;
    }

}; // class AdaptiveDihedralRestraint



template <class DType, class RType>
class DihedralRestraintMgr_Base
    : public DestructionObserver, public Dihedral_Restraint_Change_Mgr
{
public:
    DihedralRestraintMgr_Base() {}
    ~DihedralRestraintMgr_Base();
    DihedralRestraintMgr_Base(Structure *atomic_model, Change_Tracker *change_tracker)
        : _atomic_model(atomic_model)
    {
        set_change_tracker(change_tracker);
    }
    Structure* structure() const { return _atomic_model; }
    RType* get_restraint(DType *dihedral, bool create);
    size_t num_restraints() const { return _dihedral_to_restraint.size(); }
    std::vector<RType *> visible_restraints() const;
    // Change_Tracker* change_tracker() const { return _change_tracker; }
    void delete_restraints(const std::set<RType *>& to_delete);
    virtual void destructors_done(const std::set<void *>& deleted);

protected:
    RType* new_restraint(DType *d);


private:
    std::unordered_map<DType*, RType*> _dihedral_to_restraint;
    Structure* _atomic_model;
    // Change_Tracker* _change_tracker;
    // colors::variable_colormap _colormap;
    const std::string _py_name;
    const std::string _managed_class_py_name;
    const char* error_duplicate() const {
        return "A restraint already exists for this dihedral!";
    }
    const char* error_no_restraint() const {
        return "No restraint exists for this dihedral!";
    }
    RType* _new_restraint(DType *d);
    void _delete_restraints(const std::set<RType *>& to_delete);

};

class ChiralRestraintMgr:
    public DihedralRestraintMgr_Base<ChiralCenter, ChiralRestraint>,
    public pyinstance::PythonInstance<ChiralRestraintMgr>
{
public:
    ChiralRestraintMgr(Structure *s, Change_Tracker *ct)
        :DihedralRestraintMgr_Base<ChiralCenter, ChiralRestraint>(s, ct)
    {
        _mgr_type = std::type_index(typeid(this));
        _mgr_pointer = static_cast<void *>(this);
        change_tracker()->register_mgr(_mgr_type, _py_name, _managed_class_py_name);
    }
private:
    const std::string _py_name = "ChiralRestraintMgr";
    const std::string _managed_class_py_name = "ChiralRestraint";
}; // class ChiralRestraintMgr

class ProperDihedralRestraintMgr:
    public DihedralRestraintMgr_Base<ProperDihedral, ProperDihedralRestraint>,
    public pyinstance::PythonInstance<ProperDihedralRestraintMgr>
{
public:
    ProperDihedralRestraintMgr(Structure *s, Change_Tracker *ct)
        :DihedralRestraintMgr_Base<ProperDihedral, ProperDihedralRestraint>(s, ct)
    {
        _mgr_type = std::type_index(typeid(this));
        _mgr_pointer = static_cast<void *>(this);
        change_tracker()->register_mgr(_mgr_type, _py_name, _managed_class_py_name);
    }

private:
    const std::string _py_name = "ProperDihedralRestraintMgr";
    const std::string _managed_class_py_name = "ProperDihedralRestraint";
}; //ProperDihedralRestraintMgr

class AdaptiveDihedralRestraintMgr:
    public DihedralRestraintMgr_Base<ProperDihedral, AdaptiveDihedralRestraint>,
    public pyinstance::PythonInstance<AdaptiveDihedralRestraintMgr>
{
public:
    AdaptiveDihedralRestraintMgr(Structure *s, Change_Tracker *ct)
        :DihedralRestraintMgr_Base<ProperDihedral, AdaptiveDihedralRestraint>(s, ct)
    {
        _mgr_type = std::type_index(typeid(this));
        _mgr_pointer = static_cast<void *>(this);
        change_tracker()->register_mgr(_mgr_type, _py_name, _managed_class_py_name);
    }

private:
    const std::string _py_name = "AdaptiveDihedralRestraintMgr";
    const std::string _managed_class_py_name = "AdaptiveDihedralRestraint";
}; //ProperDihedralRestraintMgr



}
#endif //ISOLDE_DIHEDRAL_RESTRAINT
