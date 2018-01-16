
#define PYINSTANCE_EXPORT
#include "dihedral.h"
#include <pyinstance/PythonInstance.instantiate.h>


template class pyinstance::PythonInstance<isolde::Dihedral>;

namespace isolde {

Dihedral::Dihedral(Atom* a1, Atom* a2, Atom* a3, Atom* a4)
{
    if (a1==a2||a1==a3||a1==a4||a2==a3||a2==a4||a3==a4)
        throw std::invalid_argument(err_msg_dup_atom());
    
    _atoms[0] = a1;
    _atoms[1] = a2;
    _atoms[2] = a3;
    _atoms[3] = a4;

    
    Structure* as = atoms()[0]->structure();
    for (int i=1; i<4; ++i) {
        if (atoms()[i]->structure() != as)
            throw std::invalid_argument(err_msg_multi_struct());
    }
}

Dihedral::Dihedral(Atom* a1, Atom* a2, Atom* a3, Atom* a4, Residue* owner, std::string name):
    Dihedral(a1, a2, a3, a4) 
{
    _residue = owner;
    _name = name;
}

void copy(const Coord& coord_from, Coord coord_to)
{
    coord_to[0] = coord_from[0];
    coord_to[1] = coord_from[1];
    coord_to[2] = coord_from[2];
    
}

const Dihedral::Coords& 
Dihedral::coords() const
{    
    copy(atoms()[0]->coord(), _coords[0]);
    copy(atoms()[1]->coord(), _coords[1]);
    copy(atoms()[2]->coord(), _coords[2]);
    copy(atoms()[3]->coord(), _coords[3]);
    return _coords;
}

Real 
Dihedral::angle() const
{
    return geometry::dihedral_angle<Coord, Real>(atoms()[0]->coord(), 
        atoms()[1]->coord(), atoms()[2]->coord(), atoms()[3]->coord());
}

Residue* Dihedral::residue() const
{
    if (_residue == nullptr) 
        throw std::runtime_error(err_msg_no_residue());
    return _residue;
}

void
Dihedral::check_destroyed_atoms(const std::set<void*>& destroyed)
{
    for (auto a: atoms()) {
        if (destroyed.find(static_cast<void*>(a)) != destroyed.end()) {
            delete this;
            return;
        }
    }
}



Proper_Dihedral::Proper_Dihedral(Atom* a1, Atom* a2, Atom* a3, Atom* a4, Residue* owner, std::string name)
:   Dihedral(a1, a2, a3, a4, owner, name)
{   
    Atom* ba1 = atoms()[0];
    bool found_bond;
    for (int i=1; i<4; ++i) {
        found_bond = false;
        Atom* ba2 = atoms()[i];
        for (auto bond: ba1->bonds()) {
            if (bond->atoms()[0] == ba2 || bond->atoms()[1] == ba2) {
                _bonds[i-1] = bond;
                found_bond = true;
                break;
            }
        }
        if (!found_bond) {
            throw std::invalid_argument(err_msg_not_bonded());
        }
        ba1 = ba2;
    }
}

} // namespace isolde
