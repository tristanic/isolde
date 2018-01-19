
#ifndef ISOLDE_DIHEDRAL_MGR
#define ISOLDE_DIHEDRAL_MGR

#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>
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
    
    // Nmap maps residue name and dihedral name to the dihedral definition
    typedef std::pair<std::vector<std::string>, std::vector<bool>> d_def;
    typedef std::unordered_map<std::string, d_def> Amap;
    typedef std::unordered_map<std::string, Amap> Nmap;
    Dihedral_Mgr() {}
    ~Dihedral_Mgr();
    void add_dihedral_def(const std::string &rname, const std::string &dname, 
        const std::vector<std::string> &anames, const std::vector<bool> &externals);
    const d_def& get_dihedral_def(const std::string &rname, const std::string &dname) {
        return _residue_name_map.at(rname).at(dname);
    }
    void add_dihedral(DType* d);
    size_t num_dihedrals() const { return _dihedrals.size(); }
    size_t num_mapped_dihedrals() const;
    DType* get_dihedral(Residue *res, const std::string &name) {return _residue_map.at(res).at(name);}
    virtual void destructors_done(const std::set<void*>& destroyed);
private:
    Rmap _residue_map;
    Nmap _residue_name_map;
    std::vector<DType*> _dihedrals;
    
    
}; //class Dihedral_Mgr;
    
} //namespace isolde

#endif //ISOLDE_DIHEDRAL_MGR

