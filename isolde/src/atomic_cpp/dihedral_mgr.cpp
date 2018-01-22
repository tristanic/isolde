
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
    _atom_to_dihedral_map.clear();
    for (auto &dm: _residue_map) {
        for (auto r: dm.second)
            delete r.second;
    }
}

template <class DType>
void Dihedral_Mgr<DType>::add_dihedral_def(const std::string &rname, 
    const std::string &dname, const std::vector<std::string> &anames, 
    const std::vector<bool> &externals)
{
    Amap &am = _residue_name_map[rname];
    try {
        d_def &dd = am.at(dname);
        throw std::runtime_error("Dihedral definition already exists!");
    } catch (std::out_of_range) {
        d_def dd(anames, externals);
        am[dname] = dd;
    }
} //add_dihedral_def

template <class DType>
size_t Dihedral_Mgr<DType>::num_mapped_dihedrals() const 
{
    size_t count = 0;
    for (auto rm: _residue_map)
        for (auto dm: rm.second)
            count++;
    return count;
}


//! Add an existing dihedral to the manager.
/*! NOTE: It's up to you to ensure the same dihedral isn't added twice,
 *  and that you are not over-writing an existing map entry!
 */
template <class DType>
void Dihedral_Mgr<DType>::add_dihedral(DType* d)
{
    //~ _dihedrals.push_back(d);
    
    // Atom to dihedral mappings for fast clean-up
    for (auto a: d->atoms()) {
        auto &dvec = _atom_to_dihedral_map[a];
        dvec.push_back(d);
    }
    
    // Add it to the Residue:name map if it has both a Residue and a name
    try {
        Residue* r = d->residue(); // returns an error if no residue assigned
        const std::string &name = d->name();
        if (name != "") {
            Dmap &dm = _residue_map[r];
            dm[name] = d;
        }
    } catch(std::runtime_error) {
        return;
    }

                
}//add_dihedral    



// Need to clear entries when Dihedral or Residue objects are deleted
template <class DType>
void Dihedral_Mgr<DType>::destructors_done(const std::set<void*>& destroyed) 
{
    std::set<DType *> to_delete;
    for (auto it=_atom_to_dihedral_map.begin(); it != _atom_to_dihedral_map.end();) {
        auto a = it->first;
        if (destroyed.find(static_cast<void*>(a)) != destroyed.end()) {
            auto &dvec = it->second;
            for (auto d: dvec) {
                _residue_map[d->residue()].erase(d->name());
                to_delete.insert(d);
            }
            it = _atom_to_dihedral_map.erase(it);
        } else {
            ++it;
        }
    }
    for (auto it = _residue_map.begin(); it != _residue_map.end(); ) {
        auto r = it->first;
        if (destroyed.find(static_cast<void *>(r)) != destroyed.end())
            it = _residue_map.erase(it);
        else
            ++it;
    }
    for (auto it = to_delete.begin(); it != to_delete.end();) {
        delete *it;
        it = to_delete.erase(it);
    }
}
    
    
    
    
    //~ bool delete_this;
    //~ for (auto it=_dihedrals.begin(); it != _dihedrals.end();) {
        //~ delete_this = false;
        //~ auto d = *it;
        //~ for (auto a: d->atoms()) {
            //~ if (destroyed.find(static_cast<void*>(a)) != destroyed.end()) {
                //~ _residue_map[d->residue()].erase(d->name());
                //~ it = _dihedrals.erase(it);
                //~ delete d;
                //~ delete_this = true;
                //~ break;
            //~ }
        //~ }
        //~ if (!delete_this) ++it;
    //~ }
//~ }
    
    
} //namespace isolde
