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
    auto db = DestructionBatcher(this);
    _atom_to_dihedral_map.clear();
    for (auto &dm: _residue_map) {
        for (auto &d: dm.second)
            delete d.second;
    }
} //~Dihedral_Mgr

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
} //num_mapped_dihedrals


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
        auto &dset = _atom_to_dihedral_map[a];
        dset.insert(d);
        _mapped_atoms.insert(a);
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

template <class DType>
void Dihedral_Mgr<DType>::delete_dihedrals(const std::vector<DType *> &delete_list)
{
    auto db = DestructionBatcher(this);
    for (auto d: delete_list) {
        for (auto a: d->atoms() {
            auto &dset = _atom_to_dihedral_map.at(a);
            dset.erase(d);
        }
        _residue_map.at(d->residue()).erase(d->name());
        delete d;
    }
} //delete_dihedrals

// Need to clear entries when Dihedral or Residue objects are deleted
template <class DType> 
void Dihedral_Mgr<DType>::destructors_done(const std::set<void*>& destroyed) 
{
    auto db = DestructionBatcher(this);
    std::set<DType *> to_delete;
    for (auto it=_mapped_atoms.begin(); it != _mapped_atoms.end();) {
        auto a = *it;
        if (destroyed.find(static_cast<void*>(a)) != destroyed.end()) {
            auto &dvec = _atom_to_dihedral_map.at(a);
            for (auto d: dvec) {
                to_delete.insert(d);
            }
            _atom_to_dihedral_map.erase(a);
            it = _mapped_atoms.erase(it);
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
    for (auto d: to_delete)
        _residue_map[d->residue()].erase(d->name());
    for (auto d: to_delete)
        delete d;
} //destructors_done
    
} //namespace isolde
