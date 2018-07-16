//  $Id: mmdb_math_rand.cpp $
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
//  **** Module  :  Rand  <implementation>
//       ~~~~~~~~~
//  **** Classes :  RandomNumber ( random number generator )
//       ~~~~~~~~~
//
//   (C) E. Krissinel  1997-2013
//
//  =================================================================
//

#include <math.h>

#include "mmdb_math_rand.h"


namespace mmdb  {

  namespace math  {

    //  ===================  RandomNumber  ==========================

    RandomNumber::RandomNumber ( long IJ, long KL ) {
      Init ( IJ,KL );
    }

    void  RandomNumber::Init ( long IJ, long KL )  {
    long      i,j,k,l,m, ii,jj;
    realtype  s,t;

      iset = 0;
      gset = 0.0;

      if ((IJ<0) || (IJ>_RN_MAX_IJ) ||
          (KL<0) || (KL>_RN_MAX_KL))  return;

      i = mod(IJ/177,177) + 2;
      j = mod(IJ,177)     + 2;
      k = mod(KL/169,178) + 1;
      l = mod(KL,169);

      for (ii=0;ii<97;ii++)  {
        s = 0.0;
        t = 0.5;
        for (jj=1;jj<=24;jj++)  {
          m = mod(mod(i*j,179)*k,179);
          i = j;
          j = k;
          k = m;
          l = mod(53*l+1,169);
          if (mod(l*m,64)>=32)  s += t;
          t *= 0.5;
        }
        U[ii] = s;
      }

      C  = 362436.0   / 16777216.0;
      CD = 7654321.0  / 16777216.0;
      CM = 16777213.0 / 16777216.0;

      I97 = 96;
      J97 = 32;

    }


    // uniform [0..1] random number generator
    realtype RandomNumber::random()  {
    realtype uni;

      uni = U[I97] - U[J97];
      if (uni<0.0) uni += 1.0;
      U[I97] = uni;
      I97--;
      if (I97<0) I97 = 96;
      J97--;
      if (J97<0) J97 = 96;
      C -= CD;
      if (C<0.0)  C += CM;
      uni -= C;
      if (uni<0.0) uni += 1.0;

      return uni;

    }


    // uniform [-1..1] random number generator
    realtype RandomNumber::srandom()  {
    realtype uni;

      uni = U[I97] - U[J97];
      if (uni<0.0) uni += 1.0;
      U[I97] = uni;
      I97--;
      if (I97<0) I97 = 96;
      J97--;
      if (J97<0) J97 = 96;
      C -= CD;
      if (C<0.0)  C += CM;
      uni -= C;
      if (uni<0.0) uni += 1.0;

      return 2.0*uni - 1.0;

    }

    // gaussian random numbers
    realtype RandomNumber::gauss_rnd()  {
    realtype  v1,v2,r,fac;
      if (iset==0)  {
        do {
          v1 = srandom();
          v2 = srandom();
          r  = v1*v1 + v2*v2;
        } while ((r>=1.0) || (r==0.0));
        fac  = sqrt(-2.0*log(r)/r);
        gset = v1*fac;
        iset = 1;
        return v2*fac;
      } else  {
        iset  = 0;
        return gset;
      }
    }

    void  RandomNumber::write ( io::RFile f )  {
    int Version=1;
      f.WriteFile ( &Version,sizeof(Version) );
      f.WriteFile ( &I97    ,sizeof(I97)     );
      f.WriteFile ( &J97    ,sizeof(J97)     );
      f.WriteFile ( U       ,sizeof(U)       );
      f.WriteFile ( &C      ,sizeof(C)       );
      f.WriteFile ( &CD     ,sizeof(CD)      );
      f.WriteFile ( &CM     ,sizeof(CM)      );
      f.WriteFile ( &gset   ,sizeof(gset)    );
      f.WriteFile ( &iset   ,sizeof(iset)    );
    }

    void  RandomNumber::read ( io::RFile f )  {
    int Version;
      f.ReadFile ( &Version,sizeof(Version) );
      f.ReadFile ( &I97    ,sizeof(I97)     );
      f.ReadFile ( &J97    ,sizeof(J97)     );
      f.ReadFile ( U       ,sizeof(U)       );
      f.ReadFile ( &C      ,sizeof(C)       );
      f.ReadFile ( &CD     ,sizeof(CD)      );
      f.ReadFile ( &CM     ,sizeof(CM)      );
      f.ReadFile ( &gset   ,sizeof(gset)    );
      f.ReadFile ( &iset   ,sizeof(iset)    );
    }


  }  // namespace math

}  // namespace mmdb


/*

static int m1       = 259200;
static int ia1      = 7141;
static int ic1      = 54773;
static realtype rm1 = 1.0/259200.0;

static int m2       = 134456;
static int ia2      = 8121;
static int ic2      = 28411;
static realtype rm2 = 1.0/134456.0;

static int m3       = 243000;
static int ia3      = 4561;
static int ic3      = 51349;

static int ix1 = 0;
static int ix2 = 0;
static int ix3 = 0;

static realtype R[97];

void  randomize ( int iseed )  {
int  j;
  RndInit = True;
  ix1 = mod(ic1-iseed,m1);
  ix1 = mod(ia1*ix1+ic1,m1);
  ix2 = mod(ix1,m2);
  ix1 = mod(ia1*ix1+ic1,m1);
  ix3 = mod(ix1,m3);
  for (j=0;j<97;j++)  {
    ix1  = mod(ia1*ix1+ic1,m1);
    ix2  = mod(ia2*ix2+ic2,m2);
    R[j] = (ix1+ix2*rm2)*rm1;
  }
}

realtype  rand()  {
int      j;
realtype rnd;
  if (!RndInit)  randomize();
  ix1 = mod(ia1*ix1+ic1,m1);
  ix2 = mod(ia2*ix2+ic2,m2);
  ix3 = mod(ia3*ix3+ic3,m3);
  j = 1 + (97*ix3)/m3;
  j = IMax(j,1);
  j = IMin(j,97);
  rnd = R[j-1];
  R[j] = (ix1+ix2*rm2)*rm1;
  return rnd;
}
*/

//  ===========================================================

// End of  Random_N
