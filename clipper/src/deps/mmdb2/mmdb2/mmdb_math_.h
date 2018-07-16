//  $Id: mmdb_math.h $
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
//    11.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  Math <interface>
//       ~~~~~~~~~
//  **** Functions :   mmdb::math::GetTorsion
//       ~~~~~~~~~~~   mmdb::math::GetAngle
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#ifndef  __MMDB_Math__
#define  __MMDB_Math__

#include "mmdb_mattype.h"

namespace mmdb  {

  namespace math  {

    // ------------------------------------------------------------------

    const realtype NO_TORSION = -MaxReal;

    //  U[0,1,2] = x,y,z
    extern realtype GetTorsion   ( rvector U, rvector W, rvector V );
    extern realtype GetAngle     ( rvector U, rvector V );

    //  Calculates the binomial coefficient n choose m, 0<=n<=500, 0<=m<=n
    extern realtype Combinations ( int n, int m );

    //  Calculates precisely log(1-x) for x<1, including very small x
    extern realtype log1mx  ( realtype x );

    //  Calculates precisely 1 - exp(x) for any x including very small values
    extern realtype expc    ( realtype x );

    inline double exp10  ( double x ) { return exp(x*ln10);   }

    //  Calculates precisely 1-(1-x)**y including very small x and very large y
    extern realtype expc1mx ( realtype x, realtype y );

  }

}


#endif


