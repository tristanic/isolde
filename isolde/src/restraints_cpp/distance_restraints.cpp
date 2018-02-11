#define PYINSTANCE_EXPORT

#include "distance_restraints.h"
#include <pyinstance/PythonInstance.instantiate.h>
template class pyinstance::PythonInstance<isolde::Distance_Restraint>;
template class pyinstance::PythonInstance<isolde::Distance_Restraint_Mgr>;


namespace isolde {
Distance_Restraint::Distance_Restraint(Atom *a1, Atom *a2, Pseudobond *pb):
    _pbond(pb)
{
    for (auto b: a1->bonds())
        for (auto a: b->atoms())
            if (a == a2)
                throw std::logic_error(err_msg_bonded());
    _atoms[0] = a1;
    _atoms[1] = a2;
}

Distance_Restraint::Distance_Restraint(Atom *a1, Atom *a2, Pseudobond *pb,
    const double &target, const double &k)
    : Distance_Restraint(a1, a2, pb)
{
    _target = target;
    _k = k;
    set_enabled(false);
}

Distance_Restraint* Distance_Restraint_Mgr::new_restraint(Atom *a1, Atom *a2)
{
    auto it = _atom_to_restraints.find(a1);
    if (it != _atom_to_restraints.end())
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
    return _new_restraint(a1, a2);
}

Distance_Restraint* Distance_Restraint_Mgr::_new_restraint(Atom *a1, Atom *a2)
{
    if (a1->structure()!=_structure || a2->structure()!=_structure)
    {
        throw std::logic_error(error_different_mol());
        return nullptr;
    }
    if (a1 == a2)
    {
        throw std::logic_error(error_same_atom());
        return nullptr;
    }
    Pseudobond* pbond = _pbgroup->new_pseudobond(a1, a2);
    Distance_Restraint *d = new Distance_Restraint(a1, a2, pbond);
    _restraints.insert(d);
    _atom_to_restraints[a1].insert(d);
    _atom_to_restraints[a2].insert(d);
    return d;
}

Distance_Restraint* Distance_Restraint_Mgr::get_restraint(Atom *a1, Atom *a2, bool create)
{
    auto it = _atom_to_restraints.find(a1);
    if (it != _atom_to_restraints.end())
    {
        auto &dset = it->second;
        for (auto d: dset)
        {
            auto &datoms = d->atoms();
            for (auto a: datoms)
            {
                if (a == a2)
                    return d;
            }
        }
    }
    if (create)
        return _new_restraint(a1, a2);
    throw std::logic_error(error_no_restraint());
    return nullptr;
}

const std::set<Distance_Restraint *>& Distance_Restraint_Mgr::get_restraints(Atom *a) const
{
    auto it = _atom_to_restraints.find(a);
    if (it != _atom_to_restraints.end())
        return it->second;
    return _null_set;
}


Distance_Restraint_Mgr::~Distance_Restraint_Mgr()
{
    auto du = DestructionUser(this);
    delete_restraints(_restraints);

}

void Distance_Restraint_Mgr::delete_restraints(const std::set<Distance_Restraint *> &delete_list)
{
    auto db = DestructionBatcher(this);
    for (auto d: delete_list) {
        _restraints.erase(d);
        for (auto &a: d->atoms()) {
            auto it = _atom_to_restraints.find(a);
            if (it != _atom_to_restraints.end()) {
                auto &dset = it->second;
                dset.erase(d);
            }
        }
        delete d;
    }
} //delete_restraints


void Distance_Restraint_Mgr::destructors_done(const std::set<void *>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<Distance_Restraint *> to_delete;
    // Need to check for deleted atoms and delete their corresponding Distance_Restraints
    for (auto it=_atom_to_restraints.begin(); it != _atom_to_restraints.end();) {
        auto a = it->first;
        if (destroyed.find(static_cast<void *>(a)) != destroyed.end()) {
            auto &dset = it->second;
            for (auto d: dset) {
                to_delete.insert(d);
            }
            it = _atom_to_restraints.erase(it);
        } else {
            it++;
        }
    }
    delete_restraints(to_delete);
} //destructors_done
} //namespace isolde;
