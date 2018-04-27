/**
 * @Author: Tristan Croll <tic20>
 * @Date:   25-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 27-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */



#ifndef DIHEDRAL_EXT
#define DIHEDRAL_EXT

#include "dihedral.h"

#include "../molc.h"
using namespace atomstruct;
using namespace isolde;

/**************************************************
 *
 * Proper_Dihedral functions
 *
 **************************************************/

SET_PYTHON_CLASS(proper_dihedral, Proper_Dihedral)
GET_PYTHON_INSTANCES(proper_dihedral, Proper_Dihedral)


extern "C" EXPORT void
proper_dihedral_angle(void *dihedrals, size_t n, double *angles)
{
    Proper_Dihedral **d = static_cast<Proper_Dihedral **>(dihedrals);
    Dihedral **dd = reinterpret_cast<Dihedral **>(d);
    error_wrap_array_get(dd, n, &Dihedral::angle, angles);
}

extern "C" EXPORT void
proper_dihedral_name(void *dihedrals, size_t n, pyobject_t *names)
{
    Proper_Dihedral **d = static_cast<Proper_Dihedral **>(dihedrals);
    try {
        for (size_t i = 0; i<n; ++i)
            names[i] = unicode_from_string(d[i]->name());
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void proper_dihedral_residue(void *dihedrals, size_t n, pyobject_t *resp)
{
    Proper_Dihedral **d = static_cast<Proper_Dihedral **>(dihedrals);
    Dihedral **dd = reinterpret_cast<Dihedral **>(d);
    error_wrap_array_get(dd, n, &Dihedral::residue, resp);
}

extern "C" EXPORT void
proper_dihedral_atoms(void *dihedrals, size_t n, pyobject_t *atoms)
{
  Proper_Dihedral **d = static_cast<Proper_Dihedral **>(dihedrals);
  try {
      for (size_t i=0; i<n; ++i) {
          const Proper_Dihedral::Atoms &a = d[i]->atoms();
          for (size_t j=0; j<4; ++j) {
              *atoms++ = a[j];
          }
      }
  } catch (...) {
      molc_error();
  }
}


// extern "C" EXPORT void
// proper_dihedral_atoms(void *dihedrals, size_t n, pyobject_t *atoms)
// {
//   Proper_Dihedral **d = static_cast<Proper_Dihedral **>(dihedrals);
//   try {
//       for (size_t i=0; i<n; ++i) {
//           const Proper_Dihedral::Atoms &a = d[i]->atoms();
//           for (auto ta: a) {
//               *atoms++ = ta;
//           }
//       }
//   } catch (...) {
//       molc_error();
//   }
// }


extern "C" EXPORT void
proper_dihedral_axial_bond(void *dihedrals, size_t n, pyobject_t *bonds)
{
    Proper_Dihedral **d = static_cast<Proper_Dihedral **>(dihedrals);
    Dihedral **dd = reinterpret_cast<Dihedral **>(d);
    error_wrap_array_get(dd, n, &Dihedral::axial_bond, bonds);
}


#endif
