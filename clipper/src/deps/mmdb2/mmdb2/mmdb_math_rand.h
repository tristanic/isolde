//  $Id: mmdb_math_rand.h $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2013.
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
//    12.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  Rand  <interface>
//       ~~~~~~~~~
//  **** Classes :  RandomNumber ( random number generator )
//       ~~~~~~~~~
//
//   (C) E. Krissinel  1997-2013
//
//  =================================================================
//

#ifndef  __MMDB_MATH_Rand__
#define  __MMDB_MATH_Rand__

#include "mmdb_io_file.h"

namespace mmdb  {

  namespace math  {

    //  -------------------------------------------------------------

    enum RN_MAX_SEED  {
      _RN_MAX_IJ = 31328,
      _RN_MAX_KL = 30081
    };

    DefineClass(RandomNumber);

    class RandomNumber  {
      public :
        RandomNumber ( long IJ=0, long KL=0 );
        void  Init   ( long IJ=0, long KL=0 );
        realtype gauss_rnd(); //!< Gaussian random numbers
        realtype random   (); //!< Uniform [0..1] random number generator
        realtype srandom  (); //!< Uniform [-1..1] random number generator

        void  read  ( io::RFile f );
        void  write ( io::RFile f );

      protected :
        long     I97,J97;
        realtype U[97],C,CD,CM;
        realtype gset;
        long     iset;

    };

  }  // namespace math

}  // namespace mmdb

#endif
