/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#ifndef ISOLDE_ROTAMER_RESTRAINT
#define ISOLDE_ROTAMER_RESTRAINT

#include "dihedral_restraints.h"
#include "../validation/rota.h"

using namespace atomstruct;
namespace isolde
{

class Rotamer_Restraint_Mgr;

class Rotamer_Restraint
    : public pyinstance::PythonInstance<Rotamer_Restraint>
{
public:
    Rotamer_Restraint() {}
    ~Rotamer_Restraint() { auto du = DestructionUser(this); }
    Rotamer_Restraint( Rotamer *rotamer, Rotamer_Restraint_Mgr *mgr);

    Rotamer* rotamer() const { return _rotamer; }
    Residue* residue() const { return rotamer()->residue(); }
    const std::vector<Proper_Dihedral_Restraint *>& chi_restraints() const
    {
        return _chi_restraints;
    }
    Proper_Dihedral_Restraint_Mgr* dihedral_restraint_mgr() const;
    size_t n_chi() const { return rotamer()->n_chi(); }

    // Propagate a spring constant to all chi restraints
    void set_spring_constant(double k);
    // Enable/disable restraints on all chi dihedrals
    void set_enabled(bool flag);
    // A rotamer restraint will be considered enabled if all chi restraints are enabled
    bool enabled() const;
    // Set the target angles and cutoffs according to the target definition at t_index
    void set_target_index (int t_index);
    int target_index() const { return _current_target_index; }

private:
    Rotamer* _rotamer;
    Rotamer_Restraint_Mgr* _mgr;
    std::vector<Proper_Dihedral_Restraint *> _chi_restraints;
    int _current_target_index = -1;
    Rota_Target* _current_target_def;


}; // class Rotamer_Restraint

class Rotamer_Restraint_Mgr
    : public DestructionObserver, public pyinstance::PythonInstance<Rotamer_Restraint_Mgr>
{
public:
    Rotamer_Restraint_Mgr() {}
    ~Rotamer_Restraint_Mgr();
    Rotamer_Restraint_Mgr(Structure *structure, Change_Tracker *change_tracker,
        Proper_Dihedral_Restraint_Mgr *dmgr, Rota_Mgr *rmgr)
        : _structure(structure), _change_tracker(change_tracker),
          _dihedral_restraint_mgr(dmgr), _rota_mgr(rmgr)
    {
        change_tracker->register_mgr(_mgr_type, _py_name, _managed_class_py_name);
    }
    Rotamer_Restraint* new_restraint(Rotamer *rot);
    Rotamer_Restraint* get_restraint(Rotamer *rot, bool create);

    std::set<Rotamer_Restraint *> all_restraints() const;
    size_t num_restraints() const { return _restraint_map.size(); }
    void delete_restraints(const std::set<Rotamer_Restraint *>& delete_list);

    Proper_Dihedral_Restraint_Mgr* dihedral_restraint_mgr() const
    {
        return _dihedral_restraint_mgr;
    }

    Structure* structure() const { return _structure; }
    Change_Tracker* change_tracker() const { return _change_tracker; }
    void track_created(const void *r) { change_tracker()->add_created(_mgr_type, this, r); }
    void track_change(const void *r, int reason) { change_tracker()->add_modified(_mgr_type, this, r, reason); }

    Rota_Mgr* rota_mgr() const { return _rota_mgr; }

    virtual void destructors_done(const std::set<void *>& destroyed);

private:
    std::type_index _mgr_type = std::type_index(typeid(this));
    Structure *_structure;
    Change_Tracker *_change_tracker;
    Proper_Dihedral_Restraint_Mgr *_dihedral_restraint_mgr;
    Rota_Mgr *_rota_mgr;

    std::unordered_map<Rotamer*, Rotamer_Restraint*> _restraint_map;
    const std::string _py_name = "Rotamer_Restraint_Mgr";
    const std::string _managed_class_py_name = "Rotamer_Restraints";

    const char* error_duplicate() const {
        return "Restraint already defined for this rotamer!";
    }
    const char* error_different_mol() const {
        return "All rotamer restraints must be in the same model!";
    }

    const char* error_no_restraint() const {
        return "No restraint exists for this rotamer!";
    }


    Rotamer_Restraint* _new_restraint(Rotamer *rot);
    void _delete_restraints(const std::set<Rotamer_Restraint *>& delete_list);
}; //Rotamer_Restraint_Mgr

} //namespace isolde


#endif //ISOLDE_ROTAMER_RESTRAINT
