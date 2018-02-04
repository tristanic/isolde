#ifndef ISOLDE_DIHEDRAL_RESTRAINT
#define ISOLDE_DIHEDRAL_RESTRAINT

#include <cmath>

#include "../constants.h"
#include "../util.h"
#include "../atomic_cpp/dihedral.h"
#include "../colors.h"
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

class Dihedral_Restraint_Mgr;

class Dihedral_Restraint_Base
{
public:
    Dihedral_Restraint_Base() {}
    ~Dihedral_Restraint_Base() {auto du = DestructionUser(this);}
    Dihedral_Restraint_Base(Dihedral *dihedral)
        : _dihedral(dihedral) {}
    Dihedral* get_dihedral() const {return _dihedral;}
    void set_target(const double &target) {_target=util::wrapped_angle(target);}
    void set_target_deg(const double &target) {set_target(util::radians(target));}
    const double& get_target() const {return _target;}
    void enable() { _enabled=true; }
    void disable() { _enabled=false; }
    bool is_enabled() { return _enabled; }
    //! Set the restraint spring constant in kJ mol-1 rad-1
    void set_spring_constant(const double &k)
    {
        if (k<0 || std::isnan(k))
            throw std::logic_error(err_msg_negative_k());
        _spring_constant = k;
    }
    const double& get_spring_constant() const {return _spring_constant;}
    //! Returns (current angle) - (target angle) in radians
    double offset() const {return util::wrapped_angle(_dihedral->angle()-_target);}
    //! Returns (current angle) - (target angle) in degrees
    double offset_deg() const {return util::degrees(offset()); }
    //! Get the transform mapping an annotation primitive to the dihedral location
    virtual void get_annotation_transform(double *tf)
    {
        throw std::logic_error(err_msg_not_implemented());
    }
    //! Get the colour for the annotation
    virtual void get_annotation_color(double *color)
    {
        throw std::logic_error(err_msg_not_implemented());
    }

private:
    //! Wrap an angle to (-pi,pi)
    Dihedral *_dihedral;
    double _target = 0.0;
    double _spring_constant = 0.0;
    bool _enabled = false;
    const char* err_msg_negative_k() const
        { return "Spring constant must be greater than or equal to zero!";}
    const char* err_msg_not_implemented() const
        { return "Not implemented!";}

}; // Dihedral_Restraint_Base

class Proper_Dihedral_Restraint:
    public Dihedral_Restraint_Base,
    public pyinstance::PythonInstance<Proper_Dihedral_Restraint>
{

}; // Proper_Dihedral_Restraint

class Dihedral_Restraint_Mgr
    : public DestructionObserver,
      public pyinstance::PythonInstance<Dihedral_Restraint_Mgr>
{
public:
    Dihedral_Restraint_Mgr() {}
    ~Dihedral_Restraint_Mgr();


private:
    std::unordered_map<Proper_Dihedral*, Proper_Dihedral_Restraint*> _proper_dihedral_restraints;

};


}
#endif //ISOLDE_DIHEDRAL_RESTRAINT
