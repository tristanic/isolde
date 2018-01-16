
#ifndef isolde_Dihedral
#define isolde_Dihedral

#include <vector>
#include <atomstruct/destruct.h>
#include <atomstruct/AtomicStructure.h>
#include <atomstruct/Atom.h>
#include <atomstruct/Bond.h>
#include <atomstruct/Coord.h>
#include <atomstruct/Residue.h>
#include <pyinstance/PythonInstance.declare.h>

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
 */ 
class Dihedral: public DestructionObserver, public pyinstance::PythonInstance<Dihedral> 
{

public:
    typedef Atom* Atoms[4];
    typedef Coord Coords[4];
    typedef Bond* Bonds[3];
    
private:
    Atoms _atoms;
    Coords _coords;
    Bonds _bonds;
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
    Dihedral() {} // null constructor
    Dihedral(Atom* a1, Atom* a2, Atom* a3, Atom* a4);
    Dihedral(Atom* a1, Atom* a2, Atom* a3, Atom* a4, Residue* owner, std::string name);
    ~Dihedral() { auto du = DestructionUser(this); }
    const Atoms& atoms() const {return _atoms;}
    Structure* structure() const {return atoms()[0]->structure();}
    double angle() const; // return the current dihedral angle in radians
    const std::string& name() const { return _name; }
    void set_name(const std::string& name) { _name = name; }
    Residue* residue() const;
    void set_residue(Residue* residue) { _residue = residue; }
    virtual const Bonds& bonds() const { 
        throw std::invalid_argument("Base class Dihedral does not support bonds!");
    }
    virtual Bond* axial_bond() const {
        throw std::invalid_argument("Axial bond is only defined for a Proper_Dihedral!");
    }    
    
    virtual void destructors_done(const std::set<void*>& destroyed) {
        check_destroyed_atoms(destroyed);
    }
    virtual void check_destroyed_atoms(const std::set<void*>& destroyed);
    const Coords& coords() const;
    
    

}; // class Dihedral


//! Define a proper dihedral
/*!
 * Atoms must be provided in order and must all be bonded in strict 
 * order atom1--atom2--atom3--atom4. 
 */ 
class Proper_Dihedral: public Dihedral 
{

public:
    typedef Bond* Bonds[3];
    Proper_Dihedral() {} // null constructor
    Proper_Dihedral(Atom* a1, Atom* a2, Atom* a3, Atom* a4, Residue* owner, std::string name);
    const Bonds& bonds() const { return _bonds; }
    Bond* axial_bond() const { return bonds()[1]; }

private:
    Bonds _bonds;
    const char* err_msg_not_bonded() const
        {return "Atoms must be bonded a1--a2--a3--a4";}

}; // class Proper_Dihedral

} // namespace isolde

#endif // isolde_Dihedral
