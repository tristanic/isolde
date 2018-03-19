#ifndef ISOLDE_MDFF
#define ISOLDE_MDFF

#include <iostream>

#include "../constants.h"
#include "changetracker.h"
#include "sim_restraint_base.h"
#include <atomstruct/destruct.h>
#include <atomstruct/Atom.h>
#include <atomstruct/Coord.h>
#include <pyinstance/PythonInstance.declare.h>

using namespace atomstruct;

namespace isolde
{

class MDFF_Mgr;

class MDFF_Atom:
    public pyinstance::PythonInstance<MDFF_Atom>,
    public Sim_Restraint_Base
{
public:
    MDFF_Atom() : Sim_Restraint_Base () {}
    ~MDFF_Atom() { auto du=DestructionUser(this); }
    MDFF_Atom(Atom* atom, MDFF_Mgr *mgr)
        : _atom(atom), _mgr(mgr)
    {}
    void set_coupling_constant(double k);
    double get_coupling_constant () const { return _coupling_constant; }
    void set_enabled(bool flag);
    bool enabled() const { return _enabled; }
    Atom* atom() const { return _atom; }
    MDFF_Mgr* mgr() const { return _mgr; }
    Change_Tracker* change_tracker() const;


private:
    Atom* _atom;
    MDFF_Mgr *_mgr;
    double _coupling_constant = 0.0;
    bool _enabled = false;

}; //class MDFF_Atom

class MDFF_Mgr:
    public pyinstance::PythonInstance<MDFF_Mgr>,
    public DestructionObserver
{
public:
    MDFF_Mgr() {}
    ~MDFF_Mgr();
    MDFF_Mgr(Structure *atomic_model, Change_Tracker *change_tracker)
        : _atomic_model(atomic_model), _change_tracker(change_tracker)
    {
        _mgr_type = std::type_index(typeid(this));
        change_tracker->register_mgr(_mgr_type, _py_name, _managed_class_py_name);
    }

    Structure* structure() const { return _atomic_model; }
    MDFF_Atom* get_mdff_atom(Atom *atom, bool create);
    size_t num_atoms() const { return _atom_to_mdff.size(); }
    Change_Tracker* change_tracker() const { return _change_tracker; }
    void track_created(const void *r) { change_tracker()->add_created(_mgr_type, this, r); }
    void track_change(const void *r, int reason)
    {
        change_tracker()->add_modified(_mgr_type, this, r, reason);
    }
    double global_k() const { return _global_coupling_constant; }
    void set_global_k(double k) { _global_coupling_constant = k; }

    void delete_mdff_atoms(const std::set<MDFF_Atom *>& to_delete);
    virtual void destructors_done(const std::set<void *>& destroyed);

private:
    const std::string _py_name = "MDFF_Mgr";
    const std::string _managed_class_py_name = "MDFF_Atoms";
    std::type_index _mgr_type = std::type_index(typeid(this));
    double _global_coupling_constant = 1.0;

    Structure* _atomic_model;
    Change_Tracker* _change_tracker;
    std::unordered_map<Atom*, MDFF_Atom*> _atom_to_mdff;
    MDFF_Atom* _new_mdff_atom(Atom *atom);
    const char* error_different_mol() const {
        return "This atom is in the wrong structure!";
    }
    void _delete_mdff_atoms(const std::set<MDFF_Atom *>& to_delete);

}; //class MDFF_Mgr



} //namespace isolde



#endif //ISOLDE_MDFF
