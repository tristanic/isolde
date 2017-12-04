
#ifndef isolde_Dihedral
#define isolde_Dihedral

#include <vector>

#include <atomstruct/AtomicStructure.h>
#include <atomstruct/Atom.h>
#include <atomstruct/Bond.h>
#include <atomstruct/Coord.h>
#include <atomstruct/Residue.h>

#include "../geometry/geometry.h"

using namespace atomstruct;

//! Define a dihedral by four atoms.
/*!  
 * Atoms must be provided in order, such that the central pair defines
 * the dihedral axis. For the generic base class, the atoms needn't be
 * bonded to each other, but the same atom must not appear more than
 * once.
 */ 
class Dihedral {

public:
    typedef Atom* Atoms[4];
    
private:
    Atoms _atoms;
    const char* err_msg_dup_atom() const
        {return "All atoms must be unique";}
    const char* err_msg_multi_struct() const
        {return "All atoms must be in the same structure";}
    const char* err_msg_no_residue() const
        {return "This dihedral has not been attached to a residue";}
    std::string _name=""; // Name of the dihedral (e.g. phi, psi, omega, ...)
    // Most dihedrals belong to specific residues by convention, but 
    // we want to leave this optional for this base case.
    Residue* _residue = nullptr; 
    
public:
    Dihedral(Atom* a1, Atom* a2, Atom* a3, Atom* a4);
    ~Dihedral() { auto du = DestructionUser(this); }
    const Atoms& atoms() const {return _atoms;}
    Structure* structure() const {return atoms()[0]->structure();}
    double angle() const; // return the current dihedral angle in radians
    const std::string& name() const { return _name; }
    void set_name(const std::string& name) { _name = name; }
    Residue* residue() const;
    void set_residue(Residue* residue) { _residue = residue; }
    

}; // class Dihedral


//! Define a proper dihedral
/*!
 * Atoms must be provided in order and must all be bonded in strict 
 * order atom1--atom2--atom3--atom4. 
 */ 
class Proper_Dihedral: public Dihedral {

public:
    typedef Bond* Bonds[3];
    Proper_Dihedral(Atom* a1, Atom* a2, Atom* a3, Atom* a4);
    const Bonds& bonds() const { return _bonds; }
    Bond* axial_bond() const { return bonds()[1]; }

private:
    Bonds _bonds;
    const char* err_msg_not_bonded() const
        {return "Atoms must be bonded a1--a2--a3--a4";}

}; // class Proper_Dihedral



#endif // isolde_Dihedral
