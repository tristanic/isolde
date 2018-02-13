#ifndef ISOLDE_DIHEDRAL_RESTRAINT
#define ISOLDE_DIHEDRAL_RESTRAINT

#include <cmath>

#include "../constants.h"
#include "../util.h"
#include "../atomic_cpp/dihedral.h"
#include "../colors.h"
#include "changetracker.h"
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
class Dihedral_Restraint_Mgr_Base;
class Proper_Dihedral_Restraint_Mgr;

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
protected:

private:
    colors::variable_colormap _colormap;
    Change_Tracker *_change_tracker;

}; //class Dihedral_Restraint_Colormap



template <class DType>
class Dihedral_Restraint_Base
{
public:
    Dihedral_Restraint_Base() {}
    ~Dihedral_Restraint_Base() {auto du = DestructionUser(this);}
    Dihedral_Restraint_Base(DType *dihedral, Dihedral_Restraint_Change_Mgr *mgr)
        : _dihedral(dihedral), _mgr(mgr) {}
    DType* get_dihedral() const {return _dihedral;}
    void set_target(double target) {_target=util::wrapped_angle(target);}
    void set_target_deg(double target) {set_target(util::radians(target));}
    double get_target() const {return _target;}
    double get_target_deg() const { return util::degrees(_target); }
    void enable() { _enabled=true; }
    void disable() { _enabled=false; }
    void set_enabled(bool flag) { _enabled = flag; }
    void set_display(bool flag) { _display = flag; }
    bool get_display() const { return _display; }
    bool visible() const
    {
        return _enabled && _display && _dihedral->visible();
    }
    bool is_enabled() const { return _enabled; }
    //! Set the restraint spring constant in kJ mol-1 rad-1
    void set_spring_constant(const double &k)
    {
        _spring_constant = k<0 ? 0.0 : ( k > MAX_SPRING_CONSTANT ? MAX_SPRING_CONSTANT : k);
    }
    double get_spring_constant() const {return _spring_constant;}
    // Optional cutoff angle below which no force will be applied
    void set_cutoff(double cutoff) { _cutoff = cutoff; _cutoffs[1] = cutoff; }
    double get_cutoff() const { return _cutoff; }
    //! Returns (current angle) - (target angle) in radians
    double offset() const {return util::wrapped_angle(_dihedral->angle()-_target);}
    //! Returns (current angle) - (target angle) in degrees
    double offset_deg() const {return util::degrees(offset()); }
    //! Get the transform mapping an annotation primitive to the dihedral location
    virtual void get_annotation_transform(double *tf)
    {
        throw std::logic_error(err_msg_not_implemented());
    }
    virtual void get_annotation_color(uint8_t *color)
    {
        throw std::logic_error(err_msg_not_implemented());
    }

    isolde::Change_Tracker* change_tracker() const { return _mgr->change_tracker(); }
    colors::variable_colormap* colormap() const { return _mgr->colormap(); }

protected:
    double _cutoffs[3] {0, 1, M_PI}; // First and last cutoffs fixed, middle changes

private:
    DType *_dihedral;
    Dihedral_Restraint_Change_Mgr *_mgr;
    double _target = 0.0;
    double _spring_constant = 0.0;
    double _cutoff = 0.0;
    bool _enabled = false;
    bool _display = true;
    const char* err_msg_not_implemented() const
        { return "Not implemented!";}

}; // Dihedral_Restraint_Base

class Proper_Dihedral_Restraint:
    public Dihedral_Restraint_Base<Proper_Dihedral>,
    public pyinstance::PythonInstance<Proper_Dihedral_Restraint>
{
public:
    Proper_Dihedral_Restraint(Proper_Dihedral *dihedral, Dihedral_Restraint_Change_Mgr *mgr);
    void get_annotation_transform(double *tf);
    virtual void get_annotation_color(uint8_t *color);

private:
}; // Proper_Dihedral_Restraint




template <class DType, class RType>
class Dihedral_Restraint_Mgr_Base
    : public DestructionObserver, public Dihedral_Restraint_Change_Mgr
{
public:
    Dihedral_Restraint_Mgr_Base() {}
    ~Dihedral_Restraint_Mgr_Base();
    Dihedral_Restraint_Mgr_Base(Structure *atomic_model, Change_Tracker *change_tracker)
        : _atomic_model(atomic_model)
    {
        set_change_tracker(change_tracker);
        change_tracker->register_mgr(std::type_index(typeid(this)), _py_name, _managed_class_py_name);
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

class Proper_Dihedral_Restraint_Mgr:
    public Dihedral_Restraint_Mgr_Base<Proper_Dihedral, Proper_Dihedral_Restraint>,
    public pyinstance::PythonInstance<Proper_Dihedral_Restraint_Mgr>
{
public:
    Proper_Dihedral_Restraint_Mgr(Structure *s, Change_Tracker *ct)
        :Dihedral_Restraint_Mgr_Base<Proper_Dihedral, Proper_Dihedral_Restraint>(s, ct)
        {}

private:
    const std::string _py_name = "Proper_Dihedral_Restraint_Mgr";
    const std::string _managed_class_py_name = "Proper_Dihedral_Restraint";
}; //Proper_Dihedral_Restraint_Mgr

}
#endif //ISOLDE_DIHEDRAL_RESTRAINT
