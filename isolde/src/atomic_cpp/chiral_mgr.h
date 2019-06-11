/**
 * @Author: Tristan Croll <tic20>
 * @Date:   24-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 11-Jun-2019
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
    double expected_angle;
    Chiral_Def() {} // null constructor
    Chiral_Def(const std::vector<std::string>& s1, const std::vector<std::string>& s2,
        const std::vector<std::string>& s3, double angle)
    {
        substituents[0] = s1;
        substituents[1] = s2;
        substituents[2] = s3;
        expected_angle = angle;
    }
}; // struct Chiral_Def

//! Top-level manager for handling chiral centres
/*!
 * Chirals are too much of a specialisation of Dihedral to use
 * Dihedral_Mgr<Chiral_Center>
 */
class Chiral_Mgr: public DestructionObserver, public pyinstance::PythonInstance<Chiral_Mgr>
{

public:
    // maps central atom to its Chiral_Center instance
    typedef std::unordered_map<Atom*, Chiral_Center *> Amap;

    // Definitions are stored as (Residue name):(Atom name):(Expected angle)
    typedef std::unordered_map< std::string, Chiral_Def > Aname_Map;
    typedef std::unordered_map< std::string, Aname_Map > Rname_Map;


    Chiral_Mgr() {} // constructor
    ~Chiral_Mgr(); // destructor, deletes all managed Chiral_Center objects

    void add_chiral_def(const std::string& resname, const std::string& atom_name,
        const std::vector<std::string>& s1,
        const std::vector<std::string>& s2,
        const std::vector<std::string>& s3,
        double expected_angle);

    const Chiral_Def& get_chiral_def(const std::string& resname, const std::string& atom_name);
    const Chiral_Def& get_chiral_def(const ResName& resname, const AtomName& atom_name);

    //! Retrieve or create a Chiral_Center (if possible) for the given atom
    /*! A Chiral_Center will only be created if a definition exists for the
     *  given atom, and at least three non-hydrogen substituents are present.
     */
    Chiral_Center* get_chiral(Atom* center, bool create=true);

    void delete_chirals(const std::set<Chiral_Center *>& delete_list);

    size_t num_chirals() const { return _atom_to_chiral.size(); }
    virtual void destructors_done(const std::set<void*>& destroyed);

private:
    Rname_Map _defs;
    Amap _atom_to_chiral;
    Chiral_Center* _new_chiral(Atom* center);
    void _add_chiral(Chiral_Center* c);

    void _delete_chirals(const std::set<Chiral_Center *>& delete_list);





}; // class Chiral_Mgr

}

#endif //ISOLDE_CHIRAL_MGR
