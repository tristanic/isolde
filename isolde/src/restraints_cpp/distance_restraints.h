/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 28-Mar-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2019 Tristan Croll
 */



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

template<class R> class DistanceRestraintMgr_Tmpl;

/*! ChimeraX-side management of a standard harmonic distance restraint
 * Data management and visualisation of a harmonic distance restraint of the
 * form E = 1/2 * k * (r - r_0)**2.
 */
class DistanceRestraint:
    public pyinstance::PythonInstance<DistanceRestraint>,
    public Sim_Restraint_Base
{
public:
    typedef Atom* Atoms[2];
    DistanceRestraint() : Sim_Restraint_Base() {}
    ~DistanceRestraint() { auto du=DestructionUser(this); }
    // Construct a restraint with target and spring constant initialised to zero.
    DistanceRestraint(Atom *a1, Atom *a2, DistanceRestraintMgr_Tmpl<DistanceRestraint> *mgr);
    // Construct a restraint with specified target and spring constant.
    DistanceRestraint(Atom *a1, Atom *a2, DistanceRestraintMgr_Tmpl<DistanceRestraint> *mgr,
            const double &target, const double &k);

    // Get/update parameters
    double get_target() const { return _target; }
    void set_target(double target);
    double get_k() const { return _spring_constant; }
    void set_k(double k);
    bool enabled() const { return _enabled; }
    void set_enabled(bool flag);
    void set_satisfied_limit(double limit) { _satisfied_limit = limit; }
    double get_satisfied_limit() const { return _satisfied_limit; }

    // Visualisation functions
    double radius() const;
    void target_transform(float *rot44) const;
    void bond_cylinder_transform(float *rot44) const;
    bool visible() const;

    // General monitoring
    const Atoms &atoms() const {return _atoms;}
    Coord center() const { return (_atoms[0]->coord() + _atoms[1]->coord())*0.5; }
    double distance() const {return _atoms[0]->coord().distance(_atoms[1]->coord());}
    bool satisfied() const {return std::abs(distance()-_target) < _satisfied_limit; }
    Structure* structure() const {return _atoms[0]->structure();}
    Change_Tracker *change_tracker() const;
    DistanceRestraintMgr_Tmpl<DistanceRestraint> *mgr() const { return _mgr; }

private:
    void _bond_transform(float *rot44, float radius, float length_scale) const;
    Atoms _atoms;
    DistanceRestraintMgr_Tmpl<DistanceRestraint> *_mgr;
    double _target = 0;
    double _spring_constant = 0;
    bool _enabled=false;
    const char* err_msg_bonded()
    { return "Can't create a distance restraint between bonded atoms!";}
    double _satisfied_limit = 0.2;

}; // class DistanceRestraint

template<class R>
class DistanceRestraintMgr_Tmpl:
    //public pyinstance::PythonInstance<DistanceRestraintMgr_Tmpl<R>>,
    public DestructionObserver
{
public:
    DistanceRestraintMgr_Tmpl<R>() {}
    ~DistanceRestraintMgr_Tmpl<R>();
    DistanceRestraintMgr_Tmpl<R>(Structure *structure, Change_Tracker *change_tracker,
        std::type_index mgr_type, std::string py_name, std::string managed_class_py_name)
        : _structure(structure), _change_tracker(change_tracker),
          _mgr_type(mgr_type), _py_name(py_name),
          _managed_class_py_name(managed_class_py_name)
    {
        change_tracker->register_mgr(_mgr_type, _py_name, _managed_class_py_name);
    }

    R* new_restraint(Atom *a1, Atom *a2);
    R* get_restraint(Atom *a1, Atom *a2, bool create);

    const std::set<R *>& all_restraints() { return _restraints; }
    size_t num_restraints() const { return _restraints.size(); }
    const std::set<R *>& get_restraints(Atom *a) const;

    void delete_restraints(const std::set<R *> &delete_list);

    // Atom_Map maps individual atoms to a set of the DistanceRestraints they belong to
    typedef std::unordered_map<Atom*, std::set<R *> > Atom_Map;
    Structure* structure() const { return _structure; }
    Change_Tracker* change_tracker() const { return _change_tracker; }
    void track_created(const void *r) { change_tracker()->add_created(_mgr_type, this, r); }
    void track_change(const void *r, int reason) {change_tracker()->add_modified(_mgr_type, this, r, reason);}

    virtual void destructors_done(const std::set<void*>& destroyed);
protected:
    Structure* _structure;
    Change_Tracker* _change_tracker;
    std::type_index _mgr_type = std::type_index(typeid(this));
    R* _new_restraint(Atom *a1, Atom *a2);
    std::set<R *> _restraints;
    std::set<R *> _null_set;
    Atom_Map _atom_to_restraints;
    // std::set<Atom *> _mapped_atoms;
    std::string _py_name; // = "DistanceRestraintMgr";
    std::string _managed_class_py_name; // = "DistanceRestraints";
    void _delete_restraints(const std::set<R *>& delete_list);

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
}; // class DistanceRestraintMgr

class DistanceRestraintMgr:
    public DistanceRestraintMgr_Tmpl<DistanceRestraint>,
    public pyinstance::PythonInstance<DistanceRestraintMgr>
{
public:
    DistanceRestraintMgr(Structure *structure, Change_Tracker *change_tracker)
    : DistanceRestraintMgr_Tmpl<DistanceRestraint>(
        structure, change_tracker, std::type_index(typeid(this)),
        "DistanceRestraintMgr", "DistanceRestraints"
    )
    {}
private:
    std::type_index _mgr_type = std::type_index(typeid(this));
};

// TEMPLATE IMPLEMENTATIONS

template <class R>
R* DistanceRestraintMgr_Tmpl<R>::new_restraint(Atom *a1, Atom *a2)
{
    auto it = _atom_to_restraints.find(a1);
    if (it != _atom_to_restraints.end())
    {
        auto &dset = it->second;
        for (auto d: dset)
        {
            auto &datoms = d->atoms();
            for (auto a: datoms)
            {
                if (a==a2)
                {
                    throw std::logic_error(error_duplicate());
                    return nullptr;
                }
            }
        }
    }
    return _new_restraint(a1, a2);
}

template <class R>
R* DistanceRestraintMgr_Tmpl<R>::_new_restraint(Atom *a1, Atom *a2)
{
    if (a1->structure()!=_structure || a2->structure()!=_structure)
    {
        throw std::logic_error(error_different_mol());
        return nullptr;
    }
    if (a1 == a2)
    {
        throw std::logic_error(error_same_atom());
        return nullptr;
    }
    R *d = new R(a1, a2, this);
    _restraints.insert(d);
    _atom_to_restraints[a1].insert(d);
    _atom_to_restraints[a2].insert(d);
    track_created(d);
    return d;
}

template <class R>
R* DistanceRestraintMgr_Tmpl<R>::get_restraint(Atom *a1, Atom *a2, bool create)
{
    auto it = _atom_to_restraints.find(a1);
    if (it != _atom_to_restraints.end())
    {
        auto &dset = it->second;
        for (auto d: dset)
        {
            auto &datoms = d->atoms();
            for (auto a: datoms)
            {
                if (a == a2)
                    return d;
            }
        }
    }
    if (create)
        return _new_restraint(a1, a2);
    //throw std::logic_error(error_no_restraint());
    return nullptr;
}

template <class R>
const std::set<R *>& DistanceRestraintMgr_Tmpl<R>::get_restraints(Atom *a) const
{
    auto it = _atom_to_restraints.find(a);
    if (it != _atom_to_restraints.end())
        return it->second;
    return _null_set;
}

template <class R>
DistanceRestraintMgr_Tmpl<R>::~DistanceRestraintMgr_Tmpl()
{
    auto du = DestructionUser(this);
    for (auto r: _restraints) {
        delete r;
    }
    _atom_to_restraints.clear();
}

template <class R>
void DistanceRestraintMgr_Tmpl<R>::delete_restraints(const std::set<R *> &delete_list)
{
    auto db = DestructionBatcher(this);
    _delete_restraints(delete_list);
} //delete_restraints

template <class R>
void DistanceRestraintMgr_Tmpl<R>::_delete_restraints(const std::set<R *> &delete_list )
{
    for (auto d: delete_list) {
        _restraints.erase(d);
        for (auto &a: d->atoms()) {
            auto it = _atom_to_restraints.find(a);
            if (it != _atom_to_restraints.end()) {
                auto &dset = it->second;
                dset.erase(d);
                if (!dset.size())
                    _atom_to_restraints.erase(it);
            }
        }
        delete d;
    }
}

template <class R>
void DistanceRestraintMgr_Tmpl<R>::destructors_done(const std::set<void *>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<R *> to_delete;
    // Need to check for deleted atoms and delete their corresponding DistanceRestraints
    for (auto it=_atom_to_restraints.begin(); it != _atom_to_restraints.end();) {
        auto a = it->first;
        if (destroyed.find(static_cast<void *>(a)) != destroyed.end()) {
            auto &dset = it->second;
            for (auto d: dset) {
                to_delete.insert(d);
            }
            it = _atom_to_restraints.erase(it);
        } else {
            it++;
        }
    }
    _delete_restraints(to_delete);
} //destructors_done


} //namespace isolde


#endif //ISOLDE_DISTANCE_RESTRAINTS
