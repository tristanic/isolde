#define PYINSTANCE_EXPORT

#include "mdff.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::MDFF_Atom>;
template class pyinstance::PythonInstance<isolde::MDFF_Mgr>;

namespace isolde
{

Change_Tracker* MDFF_Atom::change_tracker() const
{
    return mgr()->change_tracker();
}

void MDFF_Atom::set_coupling_constant(double k)
{
    _coupling_constant = k<0 ? 0.0 : k;
    mgr()->track_change(this, change_tracker()->REASON_SPRING_CONSTANT_CHANGED);
}

void MDFF_Atom::set_enabled(bool flag)
{
    if (_enabled != flag)
    {
        _enabled = flag;
        mgr()->track_change(this, change_tracker()->REASON_ENABLED_CHANGED);
    }
}

MDFF_Atom* MDFF_Mgr::_new_mdff_atom(Atom *atom)
{
    if (atom->structure() != _atomic_model) {
        throw std::logic_error(error_different_mol());
    }
    MDFF_Atom* mdffa = new MDFF_Atom(atom, this);
    _atom_to_mdff[atom] = mdffa;
    track_created(mdffa);
    return mdffa;
}

MDFF_Atom* MDFF_Mgr::get_mdff_atom(Atom *atom, bool create)
{
    auto it = _atom_to_mdff.find(atom);
    if (it != _atom_to_mdff.end())
        return it->second;
    if (create)
        return _new_mdff_atom(atom);
    return nullptr;
}

void MDFF_Mgr::delete_mdff_atoms(const std::set<MDFF_Atom *>& to_delete)
{
    auto db = DestructionBatcher(this);
    _delete_mdff_atoms(to_delete);
}

void MDFF_Mgr::_delete_mdff_atoms(const std::set<MDFF_Atom *>& to_delete)
{
    for (auto a: to_delete) {
        _atom_to_mdff.erase(a->atom());
        delete a;
    }
}

void MDFF_Mgr::destructors_done(const std::set<void *>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<MDFF_Atom *> to_delete;
    for (const auto &it: _atom_to_mdff) {
        auto a = it.second;
        if (destroyed.find(static_cast<void *>(a->atom())) != destroyed.end()) {
            to_delete.insert(a);
        }
    }
    _delete_mdff_atoms(to_delete);
}

MDFF_Mgr::~MDFF_Mgr()
{
    auto du = DestructionUser(this);
    for (auto &it: _atom_to_mdff) {
        delete it.second;
    }
}

} //namespace isolde
