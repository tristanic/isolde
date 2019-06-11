/**
 * @Author: Tristan Croll <tic20>
 * @Date:   23-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright:2016-2019 Tristan Croll
 */



#define PYINSTANCE_EXPORT
#include "dihedral.h"
#include <pyinstance/PythonInstance.instantiate.h>


template class pyinstance::PythonInstance<isolde::Proper_Dihedral>;

namespace isolde {

Dihedral::Dihedral(Atom* a1, Atom* a2, Atom* a3, Atom* a4, Residue* owner, std::string name)
    :_residue(owner), _name(name)
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

const Dihedral::Coords&
Dihedral::coords() const
{
    util::copy_coord(atoms()[0]->coord(), _coords[0]);
    util::copy_coord(atoms()[1]->coord(), _coords[1]);
    util::copy_coord(atoms()[2]->coord(), _coords[2]);
    util::copy_coord(atoms()[3]->coord(), _coords[3]);
    return _coords;
}

double
Dihedral::angle() const
{
    return geometry::dihedral_angle<Coord, Real>(atoms()[0]->coord(),
        atoms()[1]->coord(), atoms()[2]->coord(), atoms()[3]->coord());
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
