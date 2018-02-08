
#define PYINSTANCE_EXPORT

#include "position_restraints.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::Position_Restraint>;
template class pyinstance::PythonInstance<isolde::Position_Restraint_Mgr>;

namespace isolde
{

void Position_Restraint::target_vector(double *vector) const
{
    for (size_t i=0; i<3; ++i)
        *vector++ = _target[i]-_atom->coord()[i];
}



Position_Restraint* Position_Restraint_Mgr::_new_restraint(Atom *atom, const Coord& target)
{
    if (atom->structure() != _atomic_model) {
        throw std::logic_error(error_different_mol());
    }
    if (atom->element().number() == 1) {
        throw std::logic_error(error_hydrogen());
    }
    Position_Restraint* restraint = new Position_Restraint(atom, target);
    _atom_to_restraint[atom] = restraint;
    return restraint;
}

Position_Restraint* Position_Restraint_Mgr::_new_restraint(Atom *atom)
{
    return _new_restraint(atom, atom->coord());
}

Position_Restraint* Position_Restraint_Mgr::get_restraint(Atom *atom, bool create)
{
    auto it = _atom_to_restraint.find(atom);
    if (it != _atom_to_restraint.end())
        return it->second;
    if (create)
        return _new_restraint(atom);
    return nullptr;
}

void Position_Restraint_Mgr::delete_restraints(const std::set<Position_Restraint *>& to_delete)
{
    auto db = DestructionBatcher(this);
    _delete_restraints(to_delete);
}

void Position_Restraint_Mgr::_delete_restraints(const std::set<Position_Restraint *>& to_delete)
{
    for (auto r: to_delete) {
        _atom_to_restraint.erase(r->atom());
        delete r;
    }
}

void Position_Restraint_Mgr::destructors_done(const std::set<void *>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<Position_Restraint *> to_delete;
    for (auto &it: _atom_to_restraint) {
        auto r = it.second;
        if (destroyed.find(static_cast<void *>(r->atom())) != destroyed.end()) {
            to_delete.insert(r);
        }
    }
    _delete_restraints(to_delete);
}

Position_Restraint_Mgr::~Position_Restraint_Mgr()
{
    auto du = DestructionUser(this);
    for (auto &it: _atom_to_restraint) {
        delete it.second;
    }
}

} //namespace isolde
