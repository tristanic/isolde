/**
 * @Author: Tristan Croll <tic20>
 * @Date:   23-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */



#define PYINSTANCE_EXPORT
#include "chiral.h"
#include <pyinstance/PythonInstance.instantiate.h>

template class pyinstance::PythonInstance<isolde::Chiral_Center>;

namespace isolde
{

Chiral_Center::Chiral_Center(Atom* center, Atom* s1, Atom* s2, Atom* s3, double expected_angle)
    : Dihedral(center, s1, s2, s3, center->residue(), std::string("chiral"))
{
    if (!s1->connects_to(center) || !s2->connects_to(center) || !s3->connects_to(center))
        throw std::invalid_argument(err_msg_not_bonded());
    _expected_angle = expected_angle;
}


} // namespace isolde
