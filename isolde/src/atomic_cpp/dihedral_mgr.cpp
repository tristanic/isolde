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
} //add_dihedral

//! Attempt to create a new dihedral for the given residue and parameters
/*! If successful, adds the dihedral to the Dihedral_Mgr internal mapping,
 *  and returns a pointer to the dihedral. If any dihedral atoms can't be
 *  found, returns nullptr.
 */
template <class DType>
DType* Dihedral_Mgr<DType>::new_dihedral(Residue *res, const std::string &dname, 
    const std::vector<std::string> &anames, const std::vector<bool> &external,
    const size_t &first_internal_atom)
{
    Atom* found_atoms[4];
    Atom* this_atom;
    bool found=false;
    
    found=false;
    for (auto a: res->atoms()) {
        if (a->name() == anames[first_internal_atom]) {
            found=true;
            found_atoms[first_internal_atom] = a;
            this_atom = a;
            break;
        }
    }
    if (!found) return nullptr;
    
    // Work backwards if necessary
    for (size_t j=first_internal_atom; j>0; j--) {
        found=false;
        for (auto a: this_atom->neighbors()) {
            if (a->name() == anames[j-1] && a->residue() != res) {
                found=true;
                found_atoms[j-1] = a;
                this_atom = a;
                break;
            }
        }
        if (!found) {
            break;
        }
    }
    if (!found) return nullptr;
    
    // ...and now work forwards
    this_atom = found_atoms[first_internal_atom];
    for (size_t j=first_internal_atom; j<3; j++) {
        found=false;
        size_t k = j+1;
        for (auto a: this_atom->neighbors()) {
            if ((!external[k] && a->residue() == res) ||
                 (external[k] && a->residue() != res)) {
                if (a->name() == anames[k]) {
                    found=true;
                    found_atoms[k] = a;
                    this_atom = a;
                    break;
                }
            }
        }
        if (!found) {
            break; 
        }        
    }
    if (found) {
        DType *d = new DType(found_atoms[0], 
            found_atoms[1], found_atoms[2], found_atoms[3], 
            res, dname);
        add_dihedral(d);
        return d;
    }
    return nullptr;
} //new_dihedral

//! Attempt to create a new dihedral for the given residue and name
/*! The <residue name, dihedral name> pair must already exist in the 
 *  manager's dihedral definition dict, otherwise an error will be 
 *  returned. 
 *  If successful, adds the dihedral to the Dihedral_Mgr internal mapping,
 *  and returns a pointer to the dihedral. If any dihedral atoms can't be
 *  found, returns nullptr.
 */
template <class DType>
DType* Dihedral_Mgr<DType>::new_dihedral(Residue *res, const std::string &dname)
{
    const d_def &ddef = get_dihedral_def(res->name(), dname);
    bool found=false;
    const auto &anames = ddef.first;
    const auto &external = ddef.second;
    size_t first_internal_atom = 0;
    for (; first_internal_atom < 4; ++first_internal_atom) {
        if (!external[first_internal_atom]) {
            found=true;
            break;
        }
    }
    if (!found)
        throw std::out_of_range("Unrecognised dihedral name for this residue!");
    
    return new_dihedral(res, dname, anames, external, first_internal_atom);

} //new_dihedral


            
template <class DType>
DType* Dihedral_Mgr<DType>::get_dihedral(Residue *res, const std::string &name, bool create)
{
    d_def ddef;
    try {
        return _residue_map.at(res).at(name);
    } catch (std::out_of_range) {
        if (!create)
            throw std::out_of_range("Dihedral not found!");
        try {
            ddef = get_dihedral_def(res->name(), name);
        } catch (std::out_of_range) {
            throw std::out_of_range("Dihedral name is invalid for this residue!");
        }
            
        DType *d = new_dihedral(res, name);
        return d;
    }
}
        



template <class DType>
void Dihedral_Mgr<DType>::delete_dihedrals(std::vector<DType *> &delete_list)
{
    auto db = DestructionBatcher(this);
    for (auto d: delete_list) {
        for (auto a: d->atoms()) {
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
