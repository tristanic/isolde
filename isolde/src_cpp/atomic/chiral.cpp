/**
 * @Author: Tristan Croll <tic20>
 * @Date:   23-Apr-2018
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 16-Nov-2020
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2016-2019 Tristan Croll
 */



#define PYINSTANCE_EXPORT
#include "chiral.h"
#include "../constants.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::ChiralCenter>;

namespace isolde
{

ChiralCenter::ChiralCenter(Atom* center, Atom* s1, Atom* s2, Atom* s3, double expected_angle,
    double expected_volume)
    : Dihedral(center, s1, s2, s3, center->residue(), std::string("chiral"))
{
    if (!s1->connects_to(center) || !s2->connects_to(center) || !s3->connects_to(center))
        throw std::invalid_argument(err_msg_not_bonded());
    size_t i=0;
    for (auto bond: center->bonds())
        for (auto a: bond->atoms())
        {
            if (a==s1 || a==s2 || a==s3)
            {
                _bonds[i] = bond;
                ++i;
            }
        }
    _expected_angle = expected_angle;
    _expected_volume = expected_volume;
}

Atom* ChiralCenter::fourth_substituent() const
{
    // Found live from the centre's current neighbours (never cached) so it can't
    // dangle when atoms are deleted, and reflects re-added atoms immediately.
    for (auto nb: _atoms[0]->neighbors())
    {
        if (nb != _atoms[1] && nb != _atoms[2] && nb != _atoms[3])
            return nb;
    }
    return nullptr;
}

double ChiralCenter::chiral_volume() const
{
    // Coordinates are in Angstroms here. Vectors from the chiral centre to its
    // three highest-priority substituents:
    auto c  = _atoms[0]->coord();
    auto s1 = _atoms[1]->coord();
    auto s2 = _atoms[2]->coord();
    auto s3 = _atoms[3]->coord();
    double ax = s1[0]-c[0], ay = s1[1]-c[1], az = s1[2]-c[2];
    double bx = s2[0]-c[0], by = s2[1]-c[1], bz = s2[2]-c[2];
    double cx = s3[0]-c[0], cy = s3[1]-c[1], cz = s3[2]-c[2];
    // Signed volume V = a . (b x c), with a=s1-c, b=s2-c, c=s3-c. This ordering
    // is identical to the OpenMM ChiralVolumeRestraintForce expression; keep them
    // in lock-step or the restraint will drive toward the wrong handedness.
    return ax*(by*cz - bz*cy) + ay*(bz*cx - bx*cz) + az*(bx*cy - by*cx);
}

double ChiralCenter::true_chiral_volume() const
{
    // For a 3-coordinate centre there is no 4th substituent; fall back to the
    // centre-based volume (which is unambiguous there).
    Atom* s4a = fourth_substituent();
    if (s4a == nullptr)
        return chiral_volume();
    auto s4 = s4a->coord();
    auto s1 = _atoms[1]->coord();
    auto s2 = _atoms[2]->coord();
    auto s3 = _atoms[3]->coord();
    // Signed volume of the substituent tetrahedron, apex s4:
    // (s1-s4) . [(s2-s4) x (s3-s4)].
    double ax = s1[0]-s4[0], ay = s1[1]-s4[1], az = s1[2]-s4[2];
    double bx = s2[0]-s4[0], by = s2[1]-s4[1], bz = s2[2]-s4[2];
    double cx = s3[0]-s4[0], cy = s3[1]-s4[1], cz = s3[2]-s4[2];
    return ax*(by*cz - bz*cy) + ay*(bz*cx - bx*cz) + az*(bx*cy - by*cx);
}

void ChiralCenter::flip_expected()
{
    // Negate the stored target handedness (experts-only "force" flip): the centre
    // now expects the mirror isomer. Mutates only this ChiralCenter instance --
    // the shared definition is untouched -- so a model reload restores the
    // original target.
    _expected_angle = -_expected_angle;
    if (!std::isnan(_expected_volume))
        _expected_volume = -_expected_volume;
}

double ChiralCenter::expected_volume() const
{
    // Prefer the real reference volume supplied at construction. If none was
    // given (e.g. a hand-curated definition that carries only an angle), fall
    // back to sign(expected_angle) * nominal magnitude -- the sign is correct
    // because sign(improper dihedral) == sign(signed volume) for this ordering.
    if (!std::isnan(_expected_volume))
        return _expected_volume;
    double s = (_expected_angle >= 0.0) ? 1.0 : -1.0;
    return s * IDEAL_TETRAHEDRAL_CHIRAL_VOLUME;
}


} // namespace isolde
