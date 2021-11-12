/**
 * @Author: Tristan Croll <tic20>
 * @Date:   25-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 11-Jun-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2016-2019 Tristan Croll
 */



#ifndef DIHEDRAL_EXT
#define DIHEDRAL_EXT

#include "dihedral.h"

#include "../molc.h"
using namespace atomstruct;
using namespace isolde;

/**************************************************
 *
 * ProperDihedral functions
 *
 **************************************************/

SET_PYTHON_CLASS(proper_dihedral, ProperDihedral)
GET_PYTHON_INSTANCES(proper_dihedral, ProperDihedral)


extern "C" EXPORT void
proper_dihedral_angle(void *dihedrals, size_t n, double *angles)
{
    ProperDihedral **d = static_cast<ProperDihedral **>(dihedrals);
    Dihedral **dd = reinterpret_cast<Dihedral **>(d);
    error_wrap_array_get(dd, n, &Dihedral::angle, angles);
}

extern "C" EXPORT void
proper_dihedral_name(void *dihedrals, size_t n, pyobject_t *names)
{
    ProperDihedral **d = static_cast<ProperDihedral **>(dihedrals);
    try {
        for (size_t i = 0; i<n; ++i)
            names[i] = unicode_from_string(d[i]->name());
    } catch (...) {
        molc_error();
    }
}

extern "C" EXPORT void proper_dihedral_residue(void *dihedrals, size_t n, pyobject_t *resp)
{
    ProperDihedral **d = static_cast<ProperDihedral **>(dihedrals);
    Dihedral **dd = reinterpret_cast<Dihedral **>(d);
    error_wrap_array_get(dd, n, &Dihedral::residue, resp);
}

extern "C" EXPORT void
proper_dihedral_atoms(void *dihedrals, size_t n, pyobject_t *atoms)
{
  ProperDihedral **d = static_cast<ProperDihedral **>(dihedrals);
  try {
      for (size_t i=0; i<n; ++i) {
          const ProperDihedral::Atoms &a = d[i]->atoms();
          for (size_t j=0; j<4; ++j) {
              *atoms++ = a[j];
          }
      }
  } catch (...) {
      molc_error();
  }
}

extern "C" EXPORT void
proper_dihedral_center(void *dihedrals, size_t n, double *coords)
{
    ProperDihedral **d = static_cast<ProperDihedral **>(dihedrals);
    try {
        for (size_t i=0; i<n; ++i) {
            auto center = (*d++)->center();
            for (size_t j=0; j<3; ++j)
                *coords++ = center[j];
        }
    } catch(...) {
        molc_error();
    }
}


extern "C" EXPORT void
proper_dihedral_axial_bond(void *dihedrals, size_t n, pyobject_t *bonds)
{
    ProperDihedral **d = static_cast<ProperDihedral **>(dihedrals);
    Dihedral **dd = reinterpret_cast<Dihedral **>(d);
    error_wrap_array_get(dd, n, &Dihedral::axial_bond, bonds);
}


#endif
