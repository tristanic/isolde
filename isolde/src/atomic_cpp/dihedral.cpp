
#include "dihedral.h"

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

Real Dihedral::angle() const
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



Proper_Dihedral::Proper_Dihedral(Atom* a1, Atom* a2, Atom* a3, Atom* a4)
:   Dihedral(a1, a2, a3, a4)
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
