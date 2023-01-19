/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#define PYINSTANCE_EXPORT

#include "mdff.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::MDFFAtom>;
template class pyinstance::PythonInstance<isolde::MDFFMgr>;

namespace isolde
{

MDFFAtom::MDFFAtom(Atom* atom, MDFFMgr *mgr) : _atom(atom), _mgr(mgr)
{
    _coupling_constant = _atom->element().mass();
}

Change_Tracker* MDFFAtom::change_tracker() const
{
    return mgr()->change_tracker();
}

void MDFFAtom::set_coupling_constant(double k)
{
    _coupling_constant = k<0 ? 0.0 : k;
    mgr()->track_change(this, change_tracker()->REASON_SPRING_CONSTANT_CHANGED);
}

void MDFFAtom::set_enabled(bool flag)
{
    if (_enabled != flag)
    {
        _enabled = flag;
        mgr()->track_change(this, change_tracker()->REASON_ENABLED_CHANGED);
    }
}

MDFFAtom* MDFFMgr::_new_mdff_atom(Atom *atom)
{
    if (atom->structure() != _atomic_model) {
        throw std::logic_error(error_different_mol());
    }
    MDFFAtom* mdffa = new MDFFAtom(atom, this);
    _atom_to_mdff[atom] = mdffa;
    track_created(mdffa);
    return mdffa;
}

MDFFAtom* MDFFMgr::get_mdff_atom(Atom *atom, bool create)
{
    auto it = _atom_to_mdff.find(atom);
    if (it != _atom_to_mdff.end())
        return it->second;
    if (create)
        return _new_mdff_atom(atom);
    return nullptr;
}

void MDFFMgr::delete_mdff_atoms(const std::set<MDFFAtom *>& to_delete)
{
    auto db = DestructionBatcher(this);
    _delete_mdff_atoms(to_delete);
}

void MDFFMgr::_delete_mdff_atoms(const std::set<MDFFAtom *>& to_delete)
{
    for (auto a: to_delete) {
        _atom_to_mdff.erase(a->atom());
        delete a;
    }
}

void MDFFMgr::destructors_done(const std::set<void *>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<MDFFAtom *> to_delete;
    for (const auto &it: _atom_to_mdff) {
        auto a = it.first;
        if (destroyed.find(static_cast<void *>(a)) != destroyed.end()) {
            to_delete.insert(it.second);
        }
    }
    _delete_mdff_atoms(to_delete);
}

MDFFMgr::~MDFFMgr()
{
    auto du = DestructionUser(this);
    for (auto &it: _atom_to_mdff) {
        delete it.second;
    }
    _atom_to_mdff.clear();
}

} //namespace isolde
