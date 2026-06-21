/**
 * @Author: Tristan Croll <tic20>
 * @Date:   25-Apr-2018
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 16-Nov-2020
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2016-2019 Tristan Croll
 */



#include "dihedral.h"
#include "../constants.h"

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
class ChiralCenter: public Dihedral, public pyinstance::PythonInstance<ChiralCenter>
{

public:
    ChiralCenter() {} // null constructor
    ChiralCenter(Atom* center, Atom* s1, Atom* s2, Atom* s3, double expected_angle,
        double expected_volume=NAN_NOT_SET);
    // Returns true if the central atom is visible;
    bool visible() const { return _atoms[0]->visible(); }
    double expected_angle() const { return _expected_angle; }
    double deviation() const { return util::wrapped_angle(angle()-expected_angle()); }
    //! Signed chiral volume V = (s1-c).[(s2-c)x(s3-c)] in Angstrom^3. Its sign is
    //! the handedness; it is zero at the planar inversion barrier and monotonic
    //! across it. The vector ordering MUST match the OpenMM restraint expression.
    double chiral_volume() const;
    //! Target signed volume for the correct isomer (Angstrom^3). Computed by the
    //! generator from the reference (CCD ideal) geometry and supplied at
    //! construction. If not supplied, falls back to sign(expected_angle) times a
    //! nominal ideal-tetrahedral magnitude.
    double expected_volume() const;
    //! Signed volume of the tetrahedron of the FOUR substituents
    //! (s1-s4).[(s2-s4)x(s3-s4)], using the lowest-priority (4th) substituent as
    //! apex. Unlike chiral_volume() (centre + 3 substituents) this cannot be
    //! fooled by the 4th substituent being stuck on the wrong side, so it is the
    //! robust handedness measure used for *validation* (the restraint stays on
    //! chiral_volume()). Falls back to chiral_volume() for 3-coordinate centres.
    double true_chiral_volume() const;
    //! The lowest-priority (4th) substituent, dropped from chiral_volume()'s
    //! improper. nullptr for a 3-coordinate centre. Also the natural axis for
    //! orienting validation markup. Computed from the centre's CURRENT neighbours
    //! on each call (never cached), so it can never dangle when atoms are
    //! deleted/added and always reflects the live bonding.
    Atom* fourth_substituent() const;
    Atom* chiral_atom() const { return _atoms[0]; }
    const Bonds& bonds() const { return _bonds; }

private:
    const char* err_msg_not_bonded() const
        { return "Substituent atoms must all be bonded to the chiral center"; }
    double _expected_angle;
    double _expected_volume;
    // const std::string _r = "R";
    // const std::string _s = "S";

}; // ChiralCenter



} // namespace isolde

#endif // ISOLDE_CHIRAL
