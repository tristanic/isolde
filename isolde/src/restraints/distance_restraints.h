#ifndef ISOLDE_DISTANCE_RESTRAINTS
#define ISOLDE_DISTANCE_RESTRAINTS

#include <string>
#include "../colors.h"
#include <atomstruct/destruct.h>
#include <atomstruct/string_types.h>
#include <atomstruct/Atom.h>
#include <atomstruct/Pseudobond.h>
#include <atomstruct/PBGroup.h>
#include <pyinstance/PythonInstance.declare.h>

using namespace atomstruct;

namespace isolde
{

class Distance_Restraint_Mgr;

class Distance_Restraint: public pyinstance::PythonInstance<Distance_Restraint>
{
public:
    typedef Atom* atoms[2];
    Distance_Restraint() {}
    ~Distance_Restraint() {}
    Distance_Restraint(Atom *a1, Atom *a2, Pseudobond *pbond);
    Distance_Restraint(Atom *a1, Atom *a2, Pseudobond *pbond,
            const double &target, const double &k);
    const double &get_target() const { return _target; }
    void set_target(const double &target) { _target=target; }
    const double &get_k() const { return _k; }
    void set_k(const double &k) { _k = k; }
    const atoms &atoms() const {return _atoms;}
    Pseudobond *get_pbond() const {return _pbond;}

private:
    atoms _atoms;
    Pseudobond *_pbond;
    double _target = 0;
    double _k = 0;


}; // class Distance_Restraint

class Distance_Restraint_Mgr:
    public pyinstance::PythonInstance<Distance_Restraint_Mgr>,
    public DestructionObserver
{
public:
    Distance_Restraint_Mgr() {}
    ~Distance_Restraint_Mgr();
    Distance_Restraint_Mgr(Proxy_PBGroup *pbgroup): _pbgroup(pbgroup) {}

    Distance_Restraint* new_restraint(Atom *a1, Atom *a2);
    void delete_restraints(const std::set<Distance_Restraint *> &to_delete);

    // Atom_Map maps individual atoms to a vector of the Distance_Restrants they belong to
    typedef std::unordered_map<Atom*, std::set<Distance_Restraint *> > Atom_Map;


    virtual void destructors_done(const std::set<void*>& destroyed);
private:
    Atom_Map _atom_to_restraints;
    std::set<Atom *> _mapped_atoms;
    Proxy_PBGroup* _pbgroup;

    const char* error_duplicate() const {
        return "This atom pair already has a distance restraint defined!";
    }
}; // class Distance_Restraint_Mgr


} //namespace isolde


#endif //ISOLDE_DISTANCE_RESTRAINTS
