
#ifndef isolde_Dihedral
#define isolde_Dihedral

#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <atomstruct/destruct.h>
#include <atomstruct/AtomicStructure.h>
#include <atomstruct/Atom.h>
#include <atomstruct/Bond.h>
#include <atomstruct/Coord.h>
#include <atomstruct/Residue.h>
#include <pyinstance/PythonInstance.declare.h>

#include "../constants.h"
#include "../util.h"
#include "../geometry/geometry.h"

using namespace atomstruct;


namespace isolde
{
//! Define a dihedral by four atoms.
/*!
 * Atoms must be provided in order, such that the central pair defines
 * the dihedral axis. For the generic base class, the atoms needn't be
 * bonded to each other, but the same atom must not appear more than
 * once.
 *
 * MUST BE ALLOCATED ON THE HEAP (i.e. via Dihedral* d = new Dihedral(...))
 * to work with ChimeraX's automatic clean-up system.
 * If you want the dihedral to be automatically deleted when any of its
 * atoms are deleted (HIGHLY recommended!) then it should be added to a
 * suitable Dihedral_Mgr after creation.
 */
class Dihedral
{
public:
    typedef Atom* Atoms[4];
    typedef Coord Coords[4];
    typedef Bond* Bonds[3];

protected:
    Atoms _atoms;

private:
    Coords _coords;
    Bonds _bonds;
    Residue *_residue;
    const char* err_msg_dup_atom() const
        {return "All atoms must be unique!";}
    const char* err_msg_multi_struct() const
        {return "All atoms must be in the same structure!";}
    std::string _name; // Name of the dihedral (e.g. phi, psi, omega, ...)

public:
    Dihedral(Atom* a1, Atom* a2, Atom* a3, Atom* a4, Residue* owner, std::string name);
    Dihedral() {} // null constructor
    ~Dihedral() { auto du = DestructionUser(this); }
    const Atoms& atoms() const { return _atoms; }
    Structure* structure() const { return atoms()[0]->structure(); }
    Residue* residue() const { return _residue; }
    double angle() const; // return the current dihedral angle in radians
    double angle_deg() const { return util::degrees(angle()); }
    const std::string& name() const { return _name; }
    virtual const Bonds& bonds() const {
        throw std::invalid_argument("Base class Dihedral does not support bonds!");
    }
    virtual Bond* axial_bond() const {
        throw std::invalid_argument("Axial bond is only defined for a Proper_Dihedral!");
    }

    virtual Atom* chiral_center() const {
        throw std::invalid_argument("Chiral center is only defined for a Chiral_Dihedral!");
    }
    //! Returns true only if all four atoms are visible.
    virtual bool visible() const
    {
        for (auto a: _atoms) {
            if (!(a->visible()))
                return false;
        }
        return true;
    }
    const Coords &coords() const;

}; // class Dihedral


//! Define a proper dihedral
/*!
 * Atoms must be provided in order and must all be bonded in strict
 * order atom1--atom2--atom3--atom4.
 */
class Proper_Dihedral: public Dihedral, public pyinstance::PythonInstance<Proper_Dihedral>
{

public:
    typedef Bond* Bonds[3];
    Proper_Dihedral() {} // null constructor
    Proper_Dihedral(Atom* a1, Atom* a2, Atom* a3, Atom* a4, Residue* owner, std::string name);
    const Bonds& bonds() const { return _bonds; }
    Bond* axial_bond() const { return bonds()[1]; }
    //! Returns true if both axial bond atoms are visible.
    bool visible() const
    {
        return _atoms[1]->visible() && _atoms[2]->visible();
    }

private:
    Bonds _bonds;
    const char* err_msg_not_bonded() const
        {return "Atoms must be bonded a1--a2--a3--a4";}

}; // class Proper_Dihedral

} // namespace isolde

#endif // isolde_Dihedral
