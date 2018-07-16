//  $Id: mmdb_math.cpp $
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
//  **** Module  :  Math <implementation>
//       ~~~~~~~~~
//  **** Functions :   mmdb::math::GetTorsion
//       ~~~~~~~~~~~   mmdb::math::GetAngle
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#include <math.h>

#include "mmdb_math_.h"

namespace mmdb  {

  namespace math  {

    // --------------------------------------------------------------

    realtype  GetTorsion ( rvector U, rvector W, rvector V )  {
    //      U     W      V
    //   o<----o----->o----->o
    //
    realtype A[3],B[3],C[3],Wmag,S,T;

      A[0] = U[1]*W[2] - W[1]*U[2];
      A[1] = U[2]*W[0] - W[2]*U[0];
      A[2] = U[0]*W[1] - W[0]*U[1];

      B[0] = V[1]*W[2] - W[1]*V[2];
      B[1] = V[2]*W[0] - W[2]*V[0];
      B[2] = V[0]*W[1] - W[0]*V[1];

      C[0] = A[1]*B[2] - B[1]*A[2];
      C[1] = A[2]*B[0] - B[2]*A[0];
      C[2] = A[0]*B[1] - B[0]*A[1];

      Wmag = sqrt(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]);

      S    = C[0]*W[0] + C[1]*W[1] + C[2]*W[2];
      T    = A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
      T   *= Wmag;

      if ((S==0.0) && (T==0.0))  return NO_TORSION;
                           else  return atan2(S,T);

    }


    realtype GetAngle ( rvector v1, rvector v2 )  {
    realtype l1,l2;

      l1 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
      if (l1==0.0)  l1 = 1.0;
      l2 = v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2];
      if (l2==0.0)  l2 = 1.0;

      return  acos((v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])/sqrt(l1*l2));

    }


    #define  nCombMax  500

    realtype Combinations ( int n, int m )  {
    //  0<=n<=nCombMax,  0<=m<=n
    realtype P[nCombMax+1];
    int      i,j;
      if ((m<0)  || (m>n))    return 0.0;
      if ((m==0) || (m==n))   return 1.0;
      if ((m==1) || (m==n-1)) return realtype(n);
      P[0] = 1.0;
      P[1] = 3.0;
      P[2] = 3.0;
      P[3] = 1.0;
      for (i=4;i<=n;i++)  {
        P[i] = 1.0;
        for (j=i-1;j>0;j--)
          P[j] += P[j-1];
      }
      return P[m];
    }

    realtype log1mx ( realtype x )  {
    //  Calculates precisely log(1-x) for x<1, including
    //  very small x
    realtype z,z1,z2,n;

      if (x>=1.0-10.0*MachEps)  z = -MaxReal;
      else if (fabs(x)>1.0e-8)  z = log(1.0-x);
      else  {
        z1 = x;
        z  = 0.0;
        n  = 1.0;
        do  {
          z2  = z;
          z  -= z1/n;
          z1 *= x;
          n  += 1.0;
        } while (z!=z2);
      }

      return z;

    }

    realtype expc ( realtype x )  {
    //  Calculates precisely 1 - exp(x) for any x including
    //  very small values
    realtype z,z1,z2,n;

      if (x>LnMaxReal)         z = -MaxReal;
      else if (x<-LnMaxReal)   z = 1.0;
      else if (fabs(x)>1.0e-8) z = 1.0 - Exp(x);
      else  {
        z1 = x;
        z  = x;
        n  = 1.0;
        do  {
          z2  = z;
          n  += 1.0;
          z1 *= x/n;
          z  += z1;
        } while (z!=z2);
        z = -z;
      }

      return z;

    }


    realtype expc1mx ( realtype x, realtype y )  {
    //  Calculates precisely 1-(1-x)**y including very small x and
    //  very large y
    realtype z,z1,z2,n,s;

      //  Calculate (1-x)**y as exp(y*log(1-x)).  Get log(1-x) first:
      if (x>1.0e-8)  z = log(1.0-x);
      else  {
        z1 = x;
        z  = 0.0;
        n  = 1.0;
        do  {
          z2  = z;
          z  -= z1/n;
          z1 *= x;
          n  += 1.0;
        } while (z!=z2);
      }

      //  Now calculate 1 - exp(y*log(1-x)) :
      z *= y;
      if (fabs(z)>1.0e-8)  s = 1.0 - exp(z);
      else  {
        z1 = z;
        s  = z;
        n  = 1.0;
        do  {
          z2  = s;
          n  += 1.0;
          z1 *= z/n;
          s  += z1;
        } while (s!=z2);
        s = -s;
      }

      return s;

    }

  }

}
