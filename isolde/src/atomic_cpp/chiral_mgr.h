/**
 * @Author: Tristan Croll <tic20>
 * @Date:   24-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 08-Nov-2020
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2016-2019 Tristan Croll
 */



#ifndef ISOLDE_CHIRAL_MGR
#define ISOLDE_CHIRAL_MGR

#include <array>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>

#include <atomstruct/string_types.h>
#include <atomstruct/destruct.h>
#include <atomstruct/AtomicStructure.h>
#include <atomstruct/Atom.h>
#include <atomstruct/Bond.h>
#include <atomstruct/Coord.h>
#include <atomstruct/Residue.h>
#include <pyinstance/PythonInstance.declare.h>

#include "chiral.h"

using namespace atomstruct;

namespace isolde
{

struct Chiral_Def
{
    // Since chiral centres can include bonds between residues, we need to allow
    // for multiple different names at each substituent site. This will be
    // particularly important for sugars
    std::array<std::vector<std::string>, 3> substituents;
    std::array<bool, 3> externals = {{false, false, false}};
    double expected_angle;
    Chiral_Def() {} // null constructor
    Chiral_Def(const std::vector<std::string>& s1, const std::vector<std::string>& s2,
        const std::vector<std::string>& s3, double angle)
    {
        for (const auto& s: s1)
            if (s=="*")
                throw std::runtime_error("Wildcards are only allowed for the final atom!");
        for (const auto& s: s2)
            if (s=="*")
                throw std::runtime_error("Wildcards are only allowed for the final atom!");
        substituents[0] = s1;
        substituents[1] = s2;
        substituents[2] = s3;
        expected_angle = angle;
    }
    Chiral_Def(const std::vector<std::string>& s1, const std::vector<std::string>& s2,
        const std::vector<std::string>& s3, double angle, const std::vector<bool>& ext)
        : Chiral_Def(s1, s2, s3, angle)
    {
        for (size_t i=0; i<3; ++i)
            externals[i] = ext[i];
    }
}; // struct Chiral_Def

//! Top-level manager for handling chiral centres
/*!
 * Chirals are too much of a specialisation of Dihedral to use
 * Dihedral_Mgr<ChiralCenter>
 */
class ChiralMgr: public DestructionObserver, public pyinstance::PythonInstance<ChiralMgr>
{

public:
    // maps central atom to its ChiralCenter instance
    typedef std::unordered_map<Atom*, ChiralCenter *> Amap;

    // Definitions are stored as (Residue name):(Atom name):(Expected angle)
    typedef std::unordered_map< std::string, Chiral_Def > Aname_Map;
    typedef std::unordered_map< std::string, Aname_Map > Rname_Map;


    ChiralMgr() {} // constructor
    ~ChiralMgr(); // destructor, deletes all managed ChiralCenter objects

    void add_chiral_def(const std::string& resname, const std::string& atom_name,
        const std::vector<std::string>& s1,
        const std::vector<std::string>& s2,
        const std::vector<std::string>& s3,
        double expected_angle);

    void add_chiral_def(const std::string& resname, const std::string& atom_name,
        const std::vector<std::string>& s1,
        const std::vector<std::string>& s2,
        const std::vector<std::string>& s3,
        double expected_angle,
        const std::vector<bool>& externals);


    const Chiral_Def& get_chiral_def(const std::string& resname, const std::string& atom_name);
    const Chiral_Def& get_chiral_def(const ResName& resname, const AtomName& atom_name);

    //! Retrieve or create a ChiralCenter (if possible) for the given atom
    /*! A ChiralCenter will only be created if a definition exists for the
     *  given atom, and at least three non-hydrogen substituents are present.
     */
    ChiralCenter* get_chiral(Atom* center, bool create=true);

    void delete_chirals(const std::set<ChiralCenter *>& delete_list);

    size_t num_chirals() const { return _atom_to_chiral.size(); }
    virtual void destructors_done(const std::set<void*>& destroyed);

private:
    Rname_Map _defs;
    Amap _atom_to_chiral;
    ChiralCenter* _new_chiral(Atom* center);
    void _add_chiral(ChiralCenter* c);

    void _delete_chirals(const std::set<ChiralCenter *>& delete_list);





}; // class ChiralMgr

}

#endif //ISOLDE_CHIRAL_MGR
