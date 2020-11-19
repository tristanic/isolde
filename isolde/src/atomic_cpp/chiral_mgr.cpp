/**
 * @Author: Tristan Croll <tic20>
 * @Date:   25-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 11-Nov-2020
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2016-2019 Tristan Croll
 */



#define PYINSTANCE_EXPORT

#include "chiral_mgr.h"
#include <set>
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::ChiralMgr>;

namespace isolde
{

ChiralMgr::~ChiralMgr()
{
    auto du = DestructionUser(this);
    for (auto &a: _atom_to_chiral) {
        delete a.second;
    }
    _atom_to_chiral.clear();
} //~ChiralMgr

void ChiralMgr::add_chiral_def(const std::string& resname,
    const std::string& atom_name,
    const std::vector<std::string>& s1,
    const std::vector<std::string>& s2,
    const std::vector<std::string>& s3,
    double expected_angle)
{
    _defs[resname][atom_name] = Chiral_Def(s1, s2, s3, expected_angle);
}

void ChiralMgr::add_chiral_def(const std::string& resname,
    const std::string& atom_name,
    const std::vector<std::string>& s1,
    const std::vector<std::string>& s2,
    const std::vector<std::string>& s3,
    double expected_angle,
    const std::vector<bool>& externals)
{
    _defs[resname][atom_name] = Chiral_Def(s1, s2, s3, expected_angle, externals);
}

const Chiral_Def& ChiralMgr::get_chiral_def(
    const std::string& resname, const std::string& atom_name)
{
    try {
        return _defs.at(resname).at(atom_name);
    } catch (std::out_of_range) {
        throw std::out_of_range("No chiral is registered for that atom name/residue name combination!");
    }
}

const Chiral_Def& ChiralMgr::get_chiral_def(const ResName& resname, const AtomName& atom_name)
{
    return get_chiral_def(std::string(resname), std::string(atom_name));
}


ChiralCenter* ChiralMgr::get_chiral(Atom* center, bool create)
{
    // Only atoms with 3 or more non-hydrogen substituents can be chiral centres
    if (center->bonds().size() < 3)
        return nullptr;
    auto it = _atom_to_chiral.find(center);
    if (it != _atom_to_chiral.end())
        return it->second;
    if (create)
        return _new_chiral(center);
    return nullptr;
}

ChiralCenter* ChiralMgr::_new_chiral(Atom* center)
{
    std::array<Atom*, 3> substituents;

    const auto rit = _defs.find(center->residue()->name());
    if ( rit == _defs.end() )
        return nullptr;
    const auto ait = rit->second.find(std::string(center->name()));
    if ( ait == rit->second.end() )
        return nullptr;

    const auto& def = ait->second;


    const auto& neighbors = center->neighbors();
    //const auto& def = get_chiral_def(center->residue()->name(), center->name());

    const auto& subnames = def.substituents;
    size_t i = 0;
    for (const auto& slist: subnames)
    {
        for (const auto &s: slist)
        {
            bool found = false;
            for (auto n: neighbors)
            {
                if ((def.externals[i] && n->residue()==center->residue()) ||
                    (!def.externals[i] && n->residue()!=center->residue()))
                    continue;

                if (n->name() == s)
                {
                    substituents[i] = n;
                    found = true;
                    break;
                }
                else if (s == "*")
                {
                    if (n->element().number()==1) continue;
                    for (size_t j=0; j<i; ++j)
                    {
                        if (n!=substituents[j])
                        {
                            found=true;
                            break;
                        }
                    }
                    if (found)
                    {
                        substituents[i] = n;
                        break;
                    }

                }
            }
            if (found) { ++i; continue; }
            return nullptr;
        }
    }
    ChiralCenter* c = new ChiralCenter(
        center, substituents[0], substituents[1], substituents[2], def.expected_angle);
    _add_chiral(c);
    return c;

}

void ChiralMgr::_add_chiral(ChiralCenter *c)
{
    _atom_to_chiral[c->chiral_atom()] = c;
}

void ChiralMgr::delete_chirals(const std::set<ChiralCenter *>& delete_list)
{
    auto db = DestructionBatcher(this);
    _delete_chirals(delete_list);
}

void ChiralMgr::_delete_chirals(const std::set<ChiralCenter *>& delete_list)
{
    for (auto c: delete_list) {
        _atom_to_chiral.erase(c->chiral_atom());
        delete c;
    }
}

void ChiralMgr::destructors_done(const std::set<void *>& destroyed)
{
    auto db = DestructionBatcher(this);
    std::set<ChiralCenter *> to_delete;
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
        const auto& bonds = c->bonds();
        for (auto b: bonds) {
            if (destroyed.find(static_cast<void *>(b)) != destroyed.end()) {
                to_delete.insert(c);
                break;
            }
        }
    }
    _delete_chirals(to_delete);
}


} // namespace isolde
