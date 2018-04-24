/**
 * @Author: Tristan Croll <tic20>
 * @Date:   23-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 24-Apr-2018
 * @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
 * @Copyright: Copyright 2017-2018 Tristan Croll
 */

#define PYINSTANCE_EXPORT

#include "chiral_mgr.h"
#include <set>
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::Chiral_Mgr>;

namespace isolde
{

Chiral_Mgr::~Chiral_Mgr()
{
    auto du = DestructionUser(this);
    for (auto &a: _atom_to_chiral) {
        delete a.second;
    }
} //~Chiral_Mgr

void Chiral_Mgr::add_chiral_def(const std::string& resname,
    const std::string& atom_name,
    const std::vector<std::string>& s1,
    const std::vector<std::string>& s2,
    const std::vector<std::string>& s3,
    double expected_angle)
{
    _defs[resname][atom_name] = Chiral_Def(s1, s2, s3, expected_angle);
}

const Chiral_Def& Chiral_Mgr::get_chiral_def(
    const std::string& resname, const std::string& atom_name)
{
    try {
        return _defs.at(resname).at(atom_name);
    } catch (std::out_of_range) {
        throw std::out_of_range("No chiral is registered for that atom name/residue name combination!");
    }
}

const Chiral_Def& Chiral_Mgr::get_chiral_def(const ResName& resname, const AtomName& atom_name)
{
    return get_chiral_def(std::string(resname), std::string(atom_name));
}


Chiral_Center* Chiral_Mgr::get_chiral(Atom* center, bool create)
{
    auto it = _atom_to_chiral.find(center);
    if (it != _atom_to_chiral.end())
        return it->second;
    if (create)
        return _new_chiral(center);
    return nullptr;
}

Chiral_Center* Chiral_Mgr::_new_chiral(Atom* center)
{
    std::array<Atom*, 3> substituents;
    const auto& neighbors = center->neighbors();
    const auto& def = get_chiral_def(center->residue()->name(), center->name());

    const auto& subnames = def.substituents;
    size_t i = 0;
    for (const auto& slist: subnames)
    {
        for (const auto &s: slist)
        {
            bool found = false;
            for (auto n: neighbors)
            {
                if (n->name() == s)
                {
                    substituents[i] = n;
                    found = true;
                    break;
                }
            }
            if (found) { ++i; continue; }
            return nullptr;
        }
    }
    Chiral_Center* c = new Chiral_Center(
        center, substituents[0], substituents[1], substituents[2], def.expected_angle);
    _add_chiral(c);
    return c;

}

void Chiral_Mgr::_add_chiral(Chiral_Center *c)
{
    _atom_to_chiral[c->chiral_atom()] = c;
}

void Chiral_Mgr::delete_chirals(const std::set<Chiral_Center *>& delete_list)
{
    auto db = DestructionBatcher(this);
    _delete_chirals(delete_list);
}

void Chiral_Mgr::_delete_chirals(const std::set<Chiral_Center *>& delete_list)
{
    for (auto c: delete_list) {
        _atom_to_chiral.erase(c->chiral_atom());
        delete c;
    }
}

void Chiral_Mgr::destructors_done(const std::set<void *>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<Chiral_Center *> to_delete;
    for (const auto &it: _atom_to_chiral)
    {
        auto c = it.second;
        const auto &catoms = c->atoms();
        for (auto a: catoms) {
            if (destroyed.find(static_cast<void *>(a)) != destroyed.end()) {
                to_delete.insert(c);
                break;
            }
        }
    }
    _delete_chirals(to_delete);
}


} // namespace isolde
