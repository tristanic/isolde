#include "distance_restraints.h"
#include <pyinstance/PythonInstance.instantiate.h>

Distance_Restraint::Distance_Restraint(Atom *a1, Atom *a2, Pseudobond *pb):
    _atoms{a1, a2}, _pbond(pb)
{}

Distance_Restraint::Distance_Restraint(Atom *a1, Atom *a2, Pseudobond *pb,
    const double &target, const double &k):
    _atoms{a1, a2}, _pbond(pb), _target(target), _k(k)
{}

Distance_Restraint* Distance_Restraint_Mgr::new_restraint(Atom *a1, Atom *a2)
{
    auto it = _mapped_atoms.find(a1);
    if (it != _mapped_atoms.end())
    {
        auto &dset = it->second;
        for (auto d: dset)
        {
            auto &datoms = d->atoms();
            for (auto a: datoms)
            {
                if (a==a2)
                {
                    throw std::logic_error(error_duplicate());
                    return nullptr;
                }
            }
        }
    }
    Pseudobond* pbond = _pbgroup->new_pseudobond(a1, a2);
    Distance_Restraint *d = new Distance_Restraint(a1, a2, pbond);
    _mapped_atoms.insert(a1);
    _mapped_atoms.insert(a2);
    _atom_to_restaints[a1].insert(d);
    _atom_to_restraints[a2].insert(d);
    return d;
}

void Distance_Restraint_Mgr::delete_restraints(const std::set<Distance_Restraint> &delete_list)
{
    auto db = DestructionBatcher(this);
    for (auto d: delete_list) {
        for (auto a: d->atoms()) {
            auto &dset = _atom_to_restraints.at(a);
            dset.erase(d);
        }
        delete d;
    }
} //delete_restraints


void Distance_Restraint_Mgr::destructors_done(const std::set<void *>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<Distance_Restraint *> to_delete;
    // Need to check for deleted atoms and delete their corresponding Distance_Restraints
    for (auto it=_mapped_atoms.begin(); it != _mapped_atoms.end();) {
        auto a = *it;
        if (destroyed.find(static_cast<void *>(a)) != destroyed.end()) {
            auto &dset = _atom_to_restraints.at(a);
            for (auto d: dset) {
                to_delete.insert(d);
            }
        }
    }
    delete_restraints(to_delete);
} //destructors_done
