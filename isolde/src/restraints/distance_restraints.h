#ifndef ISOLDE_DISTANCE_RESTRAINTS
#define ISOLDE_DISTANCE_RESTRAINTS

#include <string>
#include "../colors.h"
#include <atomstruct/destruct.h>
#include <atomstruct/string_types.h>
#include <atomstruct/Atom.h>
#include <atomstruct/Bond.h>
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
    typedef Atom* Atoms[2];
    Distance_Restraint() {}
    ~Distance_Restraint() { auto du=DestructionUser(this); }
    Distance_Restraint(Atom *a1, Atom *a2, Pseudobond *pbond);
    Distance_Restraint(Atom *a1, Atom *a2, Pseudobond *pbond,
            const double &target, const double &k);
    double get_target() const { return _target; }
    void set_target(double target) { _target=target; }
    double get_k() const { return _k; }
    void set_k(double k) { _k = k; }
    const Atoms &atoms() const {return _atoms;}
    Pseudobond *get_pbond() const {return _pbond;}

private:
    Atoms _atoms;
    Pseudobond *_pbond;
    double _target = 0;
    double _k = 0;
    const char* err_msg_bonded()
    { return "Can't create a distance restraint between bonded atoms!";}


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
    Distance_Restraint* get_restraint(Atom *a1, Atom *a2, bool create);

    const std::set<Distance_Restraint *>& get_restraints(Atom *a) const;

    void delete_restraints(const std::set<Distance_Restraint *> &delete_list);

    // Atom_Map maps individual atoms to a set of the Distance_Restrants they belong to
    typedef std::unordered_map<Atom*, std::set<Distance_Restraint *> > Atom_Map;



    virtual void destructors_done(const std::set<void*>& destroyed);
private:
    Distance_Restraint* _new_restraint(Atom *a1, Atom *a2);
    std::set<Distance_Restraint *> _restraints;
    std::set<Distance_Restraint *> _null_set;
    Atom_Map _atom_to_restraints;
    // std::set<Atom *> _mapped_atoms;
    Proxy_PBGroup* _pbgroup;

    const char* error_duplicate() const {
        return "This atom pair already has a distance restraint defined!";
    }
    const char* error_no_restraint() const {
        return "No restraint exists between this pair of atoms!";
    }
}; // class Distance_Restraint_Mgr


} //namespace isolde


#endif //ISOLDE_DISTANCE_RESTRAINTS
