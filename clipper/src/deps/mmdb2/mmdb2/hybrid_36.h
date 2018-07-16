//  $Id: hybrid_36.h $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//    This library is free software: you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License version 3, modified in accordance with the provisions
//    of the license to address the requirements of UK law.
//
//    You should have received a copy of the modified GNU Lesser
//    General Public License along with this library. If not, copies
//    may be downloaded from http://www.ccp4.ac.uk/ccp4license.php
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//  =================================================================
//

/* If you change the include guards, please be sure to also rename the
   functions below. Otherwise your project will clash with the original
   iotbx declarations and definitions.
 */
#ifndef IOTBX_PDB_HYBRID_36_C_H
#define IOTBX_PDB_HYBRID_36_C_H

#ifdef __cplusplus
extern "C" {
#endif

const char*
hy36encode(unsigned width, int value, char* result);

const char*
hy36decode(unsigned width, const char* s, unsigned s_size, int* result);

#ifdef __cplusplus
}
#endif
#endif /* IOTBX_PDB_HYBRID_36_C_H */
