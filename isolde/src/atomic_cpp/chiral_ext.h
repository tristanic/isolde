/**
 * @Author: Tristan Croll <tic20>
 * @Date:   25-Apr-2018
 * @Email:  tcroll@altoslabs.com
 * @Last modified by:   tic20
 * @Last modified time: 11-Jun-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2016-2019 Tristan Croll
 */



#ifndef CHIRAL_EXT
#define CHIRAL_EXT

#include "chiral.h"

#include "../molc.h"
using namespace atomstruct;
using namespace isolde;


SET_PYTHON_CLASS(chiral_center, ChiralCenter)
GET_PYTHON_INSTANCES(chiral_center, ChiralCenter)

extern "C" EXPORT void
chiral_angle(void *chirals, size_t n, double *angles)
{
    ChiralCenter **c = static_cast<ChiralCenter **>(chirals);
    Dihedral **d = reinterpret_cast<Dihedral **>(c);
    error_wrap_array_get(d, n, &Dihedral::angle, angles);
}

extern "C" EXPORT void chiral_residue(void *chirals, size_t n, pyobject_t *resp)
{
    ChiralCenter **c = static_cast<ChiralCenter **>(chirals);
    Dihedral **d = reinterpret_cast<Dihedral **>(c);
    error_wrap_array_get(d, n, &Dihedral::residue, resp);
}

extern "C" EXPORT void
chiral_atoms(void *chirals, size_t n, pyobject_t *atoms)
{
  ChiralCenter **c = static_cast<ChiralCenter **>(chirals);
  try {
      for (size_t i=0; i<n; ++i) {
          const ChiralCenter::Atoms &a = c[i]->atoms();
          for (auto ta: a) {
              *atoms++ = ta;
          }
      }
  } catch (...) {
      molc_error();
  }
}

extern "C" EXPORT void
chiral_center_expected_angle(void *chirals, size_t n, double *angles)
{
    ChiralCenter **c = static_cast<ChiralCenter **>(chirals);
    error_wrap_array_get(c, n, &ChiralCenter::expected_angle, angles);
}

extern "C" EXPORT void
chiral_center_deviation(void *chirals, size_t n, double *angles)
{
    ChiralCenter **c = static_cast<ChiralCenter **>(chirals);
    error_wrap_array_get(c, n, &ChiralCenter::deviation, angles);
}


extern "C" EXPORT void
chiral_center_chiral_atom(void *chirals, size_t n, pyobject_t *atoms)
{
    ChiralCenter **c = static_cast<ChiralCenter **>(chirals);
    try {
        for (size_t i=0; i<n; ++i) {
            *atoms++ = (*c++)->chiral_atom();
        }
    } catch(...) {
        molc_error();
    }
}

extern "C" EXPORT void
chiral_center_center(void *chirals, size_t n, double *coords)
{
    ChiralCenter **c = static_cast<ChiralCenter **>(chirals);
    try {
        for (size_t i=0; i<n; ++i) {
            auto center = (*c++)->center();
            for (size_t j=0; j<3; ++j)
                *coords++ = center[j];
        }
    } catch(...) {
        molc_error();
    }
}


#endif //CHIRAL_EXT
