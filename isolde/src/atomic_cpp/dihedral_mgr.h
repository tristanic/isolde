
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
    typedef std::unordered_map<std::string, DType*> DMap;
    typedef std::unordered_map<Residue*, DMap> Rmap;
    Dihedral_Mgr() {}
    ~Dihedral_Mgr();
    void add_dihedral(DType* d);
    size_t num_dihedrals() const { return _dihedrals.size(); }
    DType* get_dihedral(Residue *res, const std::string &name) {return _residue_map.at(res).at(name);}
    virtual void destructors_done(const std::set<void*>& destroyed);
private:
    Rmap _residue_map;
    std::vector<DType*> _dihedrals;
    
    
}; //class Dihedral_Mgr;
    
} //namespace isolde

#endif //ISOLDE_DIHEDRAL_MGR

