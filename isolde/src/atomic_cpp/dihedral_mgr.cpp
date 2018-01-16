
#include "dihedral_mgr.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class isolde::Dihedral_Mgr<isolde::Proper_Dihedral>;


template class pyinstance::PythonInstance<isolde::Dihedral_Mgr<isolde::Proper_Dihedral>>;

namespace isolde 
{

template <class DType>
Dihedral_Mgr<DType>::~Dihedral_Mgr() 
{
    auto du = DestructionUser(this);
    for (auto d: _dihedrals) {
        delete d;
    }
}

template <class DType>
void Dihedral_Mgr<DType>::add_dihedral(DType* d)
{
    // Add the new Dihedral to the vector if not already stored
    if (std::find(_dihedrals.begin(), _dihedrals.end(), d) == _dihedrals.end()) {
        _dihedrals.push_back(d);
    }
    // Add it to the Residue:name map if it has both a Residue and a name
    try {
        Residue* r = d->residue(); // returns an error if no residue assigned
        const std::string &name = d->name();
        if (name != "") {
            DMap &dm = _residue_map[r];
            dm[name] = d;
        }
    } catch(std::runtime_error) {
        return;
    }

                
}//add_dihedral    



// Need to clear entries when Dihedral or Residue objects are deleted
template <class DType>
void Dihedral_Mgr<DType>::destructors_done(const std::set<void*>& destroyed) {
    for (auto it1 = _residue_map.begin(); it1 != _residue_map.end();) {
        auto dm = it1->second;
        for (auto it2=dm.begin(); it2 != dm.end();) {
            auto dp = it2->second;
            if (destroyed.find(static_cast<void*>(dp)) != destroyed.end())
                it2 = dm.erase(it2);
            else
                it2++;
        }
        auto r = it1->first;
        if (destroyed.find(static_cast<void *>(r)) != destroyed.end())
            it1 = _residue_map.erase(it1);
        else
            it1++;
    }

    for (auto it=_dihedrals.begin(); it != _dihedrals.end();) {
        if (destroyed.find(static_cast<void*>(*it)) != destroyed.end()) {
            it = _dihedrals.erase(it);
        } else {
            ++it;
        }
    }

}
    
    
} //namespace isolde
