
#ifndef ISOLDE_DIHEDRAL_MGR
#define ISOLDE_DIHEDRAL_MGR

#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <atomstruct/string_types.h>
#include <atomstruct/destruct.h>
#include <atomstruct/AtomicStructure.h>
#include <atomstruct/Atom.h>
#include <atomstruct/Bond.h>
#include <atomstruct/Coord.h>
#include <atomstruct/Residue.h>
#include <pyinstance/PythonInstance.declare.h>

#include "../geometry/geometry.h"
#include "dihedral.h"

using namespace atomstruct;

namespace isolde
{

//! Top-level manager for handling all dihedrals of a given type for a model.
/*!
 * Implemented for Proper_Dihedral and Improper_Dihedral classes.
 */
template<class DType>
class Dihedral_Mgr: public DestructionObserver, public pyinstance::PythonInstance<Dihedral_Mgr<DType>>
{
public:
    // Rmap maps residue instance and dihedral name to a dihedral instance
    typedef std::unordered_map<std::string, DType*> Dmap;
    typedef std::unordered_map<Residue*, Dmap> Rmap;

    // Atom_Map maps individual atoms to a vector of the dihedral(s) they belong to
    typedef std::unordered_map<Atom*, std::set<DType *> > Atom_Map;

    // Nmap maps residue name and dihedral name to the dihedral definition
    typedef std::pair<std::vector<std::string>, std::vector<bool>> d_def;
    typedef std::unordered_map<std::string, d_def> Amap;
    typedef std::unordered_map<std::string, Amap> Nmap;
    Dihedral_Mgr() {}
    ~Dihedral_Mgr();
    void add_dihedral_def(const std::string &rname, const std::string &dname,
        const std::vector<std::string> &anames, const std::vector<bool> &externals);
    const d_def& get_dihedral_def(const std::string &rname, const std::string &dname) {
        try {
            return _residue_name_map.at(rname).at(dname);
        } catch (std::out_of_range) {
            throw std::out_of_range("No dihedral with that name is registered for this residue type!");
        }
    }

    const d_def& get_dihedral_def(const ResName &rname, const std::string &dname) {
        return get_dihedral_def(std::string(rname), dname);
    }

    //! Attempt to create a new dihedral for the given residue and name
    /*! The <residue name, dihedral name> pair must already exist in the
     *  manager's dihedral definition dict, otherwise an error will be
     *  returned.
     *  If successful, adds the dihedral to the Dihedral_Mgr internal mapping,
     *  and returns a pointer to the dihedral. If any dihedral atoms can't be
     *  found, returns nullptr.
     */
    DType* new_dihedral(Residue *res, const std::string &dname);

    //! Attempt to create a new dihedral for the given residue and parameters
    /*! If successful, adds the dihedral to the Dihedral_Mgr internal mapping,
     *  and returns a pointer to the dihedral. If any dihedral atoms can't be
     *  found, returns nullptr.
     */
    DType* new_dihedral(Residue *res, const std::string &dname,
        const std::vector<std::string> &anames, const std::vector<bool> &external,
        const size_t &first_internal_atom);
    size_t size() const {return _residue_map.size();}
    size_t bucket_count() const {return _residue_map.bucket_count();}
    void reserve(const size_t &n) {_residue_map.reserve(n);}


    void delete_dihedrals(const std::set<DType *> &delete_list);

    size_t num_mapped_dihedrals() const;

    //! Retrieve a dihedral by residue and name
    /*! If the dihedral is not found and create is false, throws
     *  std::out_of_range. If not found and create is true, attempts to
     *  find the dihedral definition by residue name and dihedral name.
     *  If no definition exists for the desired dihedral, throws
     *  std::out_of_range. If the definition exists but one or more atoms
     *  are missing, returns nullptr.
     */
    DType* get_dihedral(Residue *res, const std::string &name, bool create=true);

    //! Get all existing dihedrals belonging to a given residue
    std::vector<DType *> get_dihedrals(Residue *res) const;
    virtual void destructors_done(const std::set<void*>& destroyed);
private:
    //! Add an existing dihedral to the manager.
    /*! NOTE: It's up to you to ensure the same dihedral isn't added twice,
     *  and that you are not over-writing an existing map entry!
     */
    void add_dihedral(DType* d);
    Rmap _residue_map;
    Nmap _residue_name_map;
    Atom_Map _atom_to_dihedral_map;
    std::set<Atom *> _mapped_atoms;
    void _delete_dihedrals(const std::set<DType *> &delete_list);

}; //class Dihedral_Mgr;

typedef Dihedral_Mgr<Proper_Dihedral> Proper_Dihedral_Mgr;

} //namespace isolde

#endif //ISOLDE_DIHEDRAL_MGR
