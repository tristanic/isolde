#ifndef ISOLDE_DIHEDRAL_RESTRAINT
#define ISOLDE_DIHEDRAL_RESTRAINT

#include <cmath>

#include "../constants.h"
#include "../atomic_cpp/dihedral.h"
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
    Dihedral_Restraint_Base(Dihedral *dihedral, Dihedral_Restraint_Mgr *mgr)
        : _dihedral(dihedral), _mgr(mgr) {}
    Dihedral* get_dihedral() const {return _dihedral;}
    void set_target(const double &target) {_target=wrapped_angle(target);}
    const double& get_target() const {return _target;}
    void set_spring_constant(const double &k) {_spring_constant = k;}
    const double& get_spring_constant() const {return _spring_constant;}
    double offset() const {return wrapped_angle(_dihedral->angle()-_target);}

private:
    //! Wrap an angle to (-pi,pi)
    double wrapped_angle(const double &angle) const { return remainder(angle, TWO_PI); }
    Dihedral *_dihedral;
    double _target = 0.0;
    double _spring_constant = 0.0;
    Dihedral_Restraint_Mgr *_mgr;

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
