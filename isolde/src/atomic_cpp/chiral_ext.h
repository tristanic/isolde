/**
 * @Author: Tristan Croll
 * @Date:   03-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 25-Apr-2018
 * @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
 * @Copyright: Copyright 2017-2018 Tristan Croll
 */



#ifndef CHIRAL_EXT
#define Chiral_EXT

#include "chiral.h"

#include "../molc.h"
using namespace atomstruct;
using namespace isolde;


SET_PYTHON_CLASS(chiral_center, Chiral_Center)
GET_PYTHON_INSTANCES(chiral_center, Chiral_Center)

//NOTE: basic measurement functions are handled by the generic Dihedral class in
//      dihedral_ext.h

extern "C" EXPORT void
chiral_center_expected_angle(void *chirals, size_t n, double *angles)
{
    Chiral_Center **c = static_cast<Chiral_Center **>(chirals);
    error_wrap_array_get(c, n, &Chiral_Center::expected_angle, angles);
}

extern "C" EXPORT void
chiral_center_deviation(void *chirals, size_t n, double *angles)
{
    Chiral_Center **c = static_cast<Chiral_Center **>(chirals);
    error_wrap_array_get(c, n, &Chiral_Center::deviation, angles);
}


extern "C" EXPORT void
chiral_center_chiral_atom(void *chirals, size_t n, pyobject_t *atoms)
{
    Chiral_Center **c = static_cast<Chiral_Center **>(chirals);
    try {
        for (size_t i=0; i<n; ++i) {
            *atoms++ = (*c++)->chiral_atom();
        }
    } catch(...) {
        molc_error();
    }
}



#endif //CHIRAL_EXT
