/**
 * @Author: Tristan Croll <tic20>
 * @Date:   23-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 25-Apr-2018
 * @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
 * @Copyright: Copyright 2017-2018 Tristan Croll
 */

#include "dihedral.h"

#ifndef ISOLDE_CHIRAL
#define ISOLDE_CHIRAL

namespace isolde {

//! Define a chiral centre
/*!
 * Atoms should be provided with the chiral atom first, followed by the
 * three highest-priority substituents according to standard R-S (Cahn-Ingold-
 * Prelog) rules:
 *   - substituents with higher atomic number take precedence over those with
 *     lower atomic number
 *   - if two substituents have equal atomic number, proceed further along each
 *     chain until a difference appears. Where there is a difference, the
 *     highest atomic number substituent wins.
 *
 * Following the above rules, a (S) isomer should have a positive dihedral
 * angle. Note that according to the above rules, all chiral amino acids other
 * than CYS are (S) isomers: {CA, N, C, CB} since the sulfur attached to CB is
 * heavier than the oxygen attached to C.
 *
 * A chiral center is always owned by the residue to which the central atom
 * belongs
 */
class Chiral_Center: public Dihedral, public pyinstance::PythonInstance<Chiral_Center>
{

public:
    Chiral_Center() {} // null constructor
    Chiral_Center(Atom* center, Atom* s1, Atom* s2, Atom* s3, double expected_angle);
    // Returns true if the central atom is visible;
    bool visible() const { return _atoms[0]->visible(); }
    double expected_angle() const { return _expected_angle; }
    double deviation() const { return util::wrapped_angle(angle()-expected_angle()); }
    Atom* chiral_atom() const { return _atoms[0]; }

private:
    const char* err_msg_not_bonded() const
        { return "Substituent atoms must all be bonded to the chiral center"; }
    double _expected_angle;
    // const std::string _r = "R";
    // const std::string _s = "S";

}; // Chiral_Center



} // namespace isolde

#endif // ISOLDE_CHIRAL
