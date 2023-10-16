/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



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

class MDFFMgr;

//! Handles the coupling of an Atom to a MDFF map.
/*! Unlike in the case of restraints, MDFFAtom objects are enabled by default,
 *  and are instantiated with their coupling constants pre-set to the mass of
 *  the atom.
 */
class MDFFAtom:
    public pyinstance::PythonInstance<MDFFAtom>,
    public Sim_Restraint_Base
{
public:
    MDFFAtom() {}
    ~MDFFAtom() { auto du=DestructionUser(this); }
    MDFFAtom(Atom* atom, MDFFMgr *mgr);
    void set_coupling_constant(double k);
    double get_coupling_constant () const { return _coupling_constant; }
    void set_enabled(bool flag);
    bool enabled() const { return _enabled; }
    Atom* atom() const { return _atom; }
    MDFFMgr* mgr() const { return _mgr; }
    Change_Tracker* change_tracker() const;


private:
    Atom* _atom;
    MDFFMgr *_mgr;
    double _coupling_constant = 0.0;
    bool _enabled = true;

}; //class MDFFAtom

class MDFFMgr:
    public pyinstance::PythonInstance<MDFFMgr>,
    public DestructionObserver
{
public:
    MDFFMgr() {}
    ~MDFFMgr();
    MDFFMgr(Structure *atomic_model, Change_Tracker *change_tracker)
        : _atomic_model(atomic_model), _change_tracker(change_tracker)
    {
        _mgr_type = std::type_index(typeid(this));
        change_tracker->register_mgr(_mgr_type, _py_name, _managed_class_py_name);
    }

    Structure* structure() const { return _atomic_model; }
    MDFFAtom* get_mdff_atom(Atom *atom, bool create);
    size_t num_atoms() const { return _atom_to_mdff.size(); }
    Change_Tracker* change_tracker() const { return _change_tracker; }
    void track_created(const void *r) { change_tracker()->add_created(_mgr_type, this, r); }
    void track_change(const void *r, int reason)
    {
        change_tracker()->add_modified(_mgr_type, this, r, reason);
    }
    double global_k() const { return _global_coupling_constant; }
    void set_global_k(double k) { _global_coupling_constant = k; }

    void delete_mdff_atoms(const std::set<MDFFAtom *>& to_delete);
    virtual void destructors_done(const std::set<void *>& destroyed);

private:
    const std::string _py_name = "MDFFMgr";
    const std::string _managed_class_py_name = "MDFFAtoms";
    std::type_index _mgr_type = std::type_index(typeid(this));
    double _global_coupling_constant = 1.0;

    Structure* _atomic_model;
    Change_Tracker* _change_tracker;
    std::unordered_map<Atom*, MDFFAtom*> _atom_to_mdff;
    MDFFAtom* _new_mdff_atom(Atom *atom);
    const char* error_different_mol() const {
        return "This atom is in the wrong structure!";
    }
    void _delete_mdff_atoms(const std::set<MDFFAtom *>& to_delete);

}; //class MDFFMgr



} //namespace isolde



#endif //ISOLDE_MDFF
