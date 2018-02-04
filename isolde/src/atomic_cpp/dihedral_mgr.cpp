#define PYINSTANCE_EXPORT

#include "dihedral_mgr.h"
#include <pyinstance/PythonInstance.instantiate.h>



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
    if (am.find(rname) == am.end()) {
        am[dname] = d_def(anames, externals);
    } else {
        throw std::runtime_error("Dihedral definition already exists!");
    }
} //add_dihedral_def

template <class DType>
size_t Dihedral_Mgr<DType>::num_mapped_dihedrals() const
{
    size_t count = 0;
    for (auto rm: _residue_map)
        count += rm.second.size();
        // for (auto dm: rm.second)
        //     count++;
    return count;
} //num_mapped_dihedrals


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
    auto iter1 = _residue_map.find(res);
    if (iter1 != _residue_map.end()) {
        auto &dmap = iter1->second;
        auto iter2 = dmap.find(name);
        if (iter2 != dmap.end()) {
            return iter2->second;
        }
    }
    if (!create)
        throw std::out_of_range("Dihedral not found!");
    ddef = get_dihedral_def(res->name(), name);
    return new_dihedral(res, name);
}



template <class DType>
void Dihedral_Mgr<DType>::delete_dihedrals(const std::set<DType *> &delete_list)
{
    auto db = DestructionBatcher(this);
    for (auto d: delete_list) {
        for (auto a: d->atoms()) {
            auto &dset = _atom_to_dihedral_map.at(a);
            dset.erase(d);
        }
        auto rit = _residue_map.find(d->residue());
        if (rit != _residue_map.end()) {
            _residue_map.erase(rit);
        }
        delete d;
    }
} //delete_dihedrals

// Need to clear entries when Dihedral or Residue objects are deleted
template <class DType>
void Dihedral_Mgr<DType>::destructors_done(const std::set<void*>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<DType *> to_delete;
    for (auto a: _mapped_atoms) {
        if (destroyed.find(static_cast<void*>(a)) != destroyed.end()) {
            auto &dvec = _atom_to_dihedral_map.at(a);
            for (auto d: dvec) {
                to_delete.insert(d);
            }
        }
    }
    for (auto it = _residue_map.begin(); it != _residue_map.end(); ) {
        auto r = it->first;
        if (destroyed.find(static_cast<void *>(r)) != destroyed.end())
            it = _residue_map.erase(it);
        else
            ++it;
    }

    delete_dihedrals(to_delete);
} //destructors_done

template class Dihedral_Mgr<Proper_Dihedral>;


} //namespace isolde
