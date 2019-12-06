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

class RotamerRestraintMgr;

class RotamerRestraint
    : public pyinstance::PythonInstance<RotamerRestraint>
{
public:
    RotamerRestraint() {}
    ~RotamerRestraint() { auto du = DestructionUser(this); }
    RotamerRestraint( Rotamer *rotamer, RotamerRestraintMgr *mgr);

    Rotamer* rotamer() const { return _rotamer; }
    Residue* residue() const { return rotamer()->residue(); }
    const std::vector<ProperDihedralRestraint *>& chi_restraints() const
    {
        return _chi_restraints;
    }
    ProperDihedralRestraintMgr* dihedral_restraint_mgr() const;
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

    RotamerRestraintMgr *mgr() const { return _mgr; }

private:
    Rotamer* _rotamer;
    RotamerRestraintMgr* _mgr;
    std::vector<ProperDihedralRestraint *> _chi_restraints;
    int _current_target_index = -1;
    Rota_Target* _current_target_def;


}; // class RotamerRestraint

class RotamerRestraintMgr
    : public DestructionObserver, public pyinstance::PythonInstance<RotamerRestraintMgr>
{
public:
    RotamerRestraintMgr() {}
    ~RotamerRestraintMgr();
    RotamerRestraintMgr(Structure *structure, Change_Tracker *change_tracker,
        ProperDihedralRestraintMgr *dmgr, RotaMgr *rmgr)
        : _structure(structure), _change_tracker(change_tracker),
          _dihedral_restraint_mgr(dmgr), _rota_mgr(rmgr)
    {
        change_tracker->register_mgr(_mgr_type, _py_name, _managed_class_py_name);
    }
    RotamerRestraint* new_restraint(Rotamer *rot);
    RotamerRestraint* get_restraint(Rotamer *rot, bool create);

    std::set<RotamerRestraint *> all_restraints() const;
    size_t num_restraints() const { return _restraint_map.size(); }
    void delete_restraints(const std::set<RotamerRestraint *>& delete_list);

    ProperDihedralRestraintMgr* dihedral_restraint_mgr() const
    {
        return _dihedral_restraint_mgr;
    }

    Structure* structure() const { return _structure; }
    Change_Tracker* change_tracker() const { return _change_tracker; }
    void track_created(const void *r) { change_tracker()->add_created(_mgr_type, this, r); }
    void track_change(const void *r, int reason) { change_tracker()->add_modified(_mgr_type, this, r, reason); }

    RotaMgr* rota_mgr() const { return _rota_mgr; }

    virtual void destructors_done(const std::set<void *>& destroyed);

private:
    std::type_index _mgr_type = std::type_index(typeid(this));
    Structure *_structure;
    Change_Tracker *_change_tracker;
    ProperDihedralRestraintMgr *_dihedral_restraint_mgr;
    RotaMgr *_rota_mgr;

    std::unordered_map<Rotamer*, RotamerRestraint*> _restraint_map;
    const std::string _py_name = "RotamerRestraintMgr";
    const std::string _managed_class_py_name = "RotamerRestraints";

    const char* error_duplicate() const {
        return "Restraint already defined for this rotamer!";
    }
    const char* error_different_mol() const {
        return "All rotamer restraints must be in the same model!";
    }

    const char* error_no_restraint() const {
        return "No restraint exists for this rotamer!";
    }


    RotamerRestraint* _new_restraint(Rotamer *rot);
    void _delete_restraints(const std::set<RotamerRestraint *>& delete_list);
}; //RotamerRestraintMgr

} //namespace isolde


#endif //ISOLDE_ROTAMER_RESTRAINT
