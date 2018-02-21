#ifndef ISOLDE_DISTANCE_RESTRAINTS
#define ISOLDE_DISTANCE_RESTRAINTS

#include <string>
#include "../colors.h"
#include "../constants.h"
#include "../geometry/geometry.h"
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

class Distance_Restraint_Mgr;

class Distance_Restraint:
    public pyinstance::PythonInstance<Distance_Restraint>,
    public Sim_Restraint_Base
{
public:
    typedef Atom* Atoms[2];
    Distance_Restraint() {}
    ~Distance_Restraint() { auto du=DestructionUser(this); }
    Distance_Restraint(Atom *a1, Atom *a2, Distance_Restraint_Mgr *mgr);
    Distance_Restraint(Atom *a1, Atom *a2, Distance_Restraint_Mgr *mgr,
            const double &target, const double &k);
    double get_target() const { return _target; }
    void set_target(double target);
    double get_k() const { return _spring_constant; }
    void set_k(double k);
    bool enabled() const { return _enabled; }
    void set_enabled(bool flag);
    double radius() const;
    void target_transform(double *rot44) const;
    void bond_cylinder_transform(double *rot44) const;
    bool visible() const;
    const Atoms &atoms() const {return _atoms;}
    double distance() const {return _atoms[0]->coord().distance(_atoms[1]->coord());}
    Structure* structure() const {return _atoms[0]->structure();}
    Change_Tracker *change_tracker() const;

private:
    void _bond_transform(double *rot44, double radius, double length_scale) const;
    Atoms _atoms;
    Distance_Restraint_Mgr *_mgr;
    double _target = 0;
    double _spring_constant = 0;
    bool _enabled=false;
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
    Distance_Restraint_Mgr(Structure *structure, Change_Tracker *change_tracker)
        : _structure(structure), _change_tracker(change_tracker)
    {
        change_tracker->register_mgr(_mgr_type, _py_name, _managed_class_py_name);
    }

    Distance_Restraint* new_restraint(Atom *a1, Atom *a2);
    Distance_Restraint* get_restraint(Atom *a1, Atom *a2, bool create);

    const std::set<Distance_Restraint *>& all_restraints() { return _restraints; }
    size_t num_restraints() const { return _restraints.size(); }
    const std::set<Distance_Restraint *>& get_restraints(Atom *a) const;

    void delete_restraints(const std::set<Distance_Restraint *> &delete_list);

    // Atom_Map maps individual atoms to a set of the Distance_Restrants they belong to
    typedef std::unordered_map<Atom*, std::set<Distance_Restraint *> > Atom_Map;
    Structure* structure() const { return _structure; }
    Change_Tracker* change_tracker() const { return _change_tracker; }
    void track_created(const void *r) { change_tracker()->add_created(_mgr_type, this, r); }
    void track_change(const void *r, int reason) {change_tracker()->add_modified(_mgr_type, this, r, reason);}

    virtual void destructors_done(const std::set<void*>& destroyed);
private:
    std::type_index _mgr_type = std::type_index(typeid(this));
    Distance_Restraint* _new_restraint(Atom *a1, Atom *a2);
    std::set<Distance_Restraint *> _restraints;
    std::set<Distance_Restraint *> _null_set;
    Atom_Map _atom_to_restraints;
    // std::set<Atom *> _mapped_atoms;
    Structure* _structure;
    Change_Tracker* _change_tracker;
    const std::string _py_name = "Distance_Restraint_Mgr";
    const std::string _managed_class_py_name = "Distance_Restraints";
    void _delete_restraints(const std::set<Distance_Restraint *>& delete_list);

    const char* error_same_atom() const {
        return "Cannot bond an atom to itself!";
    }
    const char* error_duplicate() const {
        return "This atom pair already has a distance restraint defined!";
    }
    const char* error_no_restraint() const {
        return "No restraint exists between this pair of atoms!";
    }
    const char* error_different_mol() const {
        return "All distance restraints must be in the same model!";
    }
}; // class Distance_Restraint_Mgr


} //namespace isolde


#endif //ISOLDE_DISTANCE_RESTRAINTS
