//  $Id: mmdb_mattype.h $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2008.
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
//    10.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  MatType_ <interface>
//       ~~~~~~~~~
//  **** Functions :
//       ~~~~~~~~~~~
//               GetString  ( reads substring from a string         )
//               GetStrTer  ( reads substring and put term-ing null )
//               strcpy_n   ( copies not more than n characters     )
//               strcpy_ns  ( like strcpy_ns and pads with spaces   )
//               strcpy_n0  ( like strcpy_n and adds terminating 0  )
//               PadSpaces  ( pads a string with spaces             )
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//


#ifndef  __MMDB_MatType__
#define  __MMDB_MatType__

#include <math.h>

#include "imex.h"

#define  UseDoubleFloat

#ifndef __ClassMacros

# define __ClassMacros

 //  A Class definition macros
# define DefineClass(ClassName)             \
   class MMDB_IMEX ClassName;                         \
   typedef ClassName    * P##ClassName;     \
   typedef ClassName    & R##ClassName;     \
   typedef P##ClassName * PP##ClassName;    \
   typedef P##ClassName & RP##ClassName;

 //  A Structure definition macros
# define DefineStructure(StructureName)             \
   struct MMDB_IMEX StructureName;                            \
   typedef StructureName    * P##StructureName;     \
   typedef StructureName    & R##StructureName;     \
   typedef P##StructureName * PP##StructureName;    \
   typedef P##StructureName & RP##StructureName;

#endif

#define UNUSED_ARGUMENT(x) (void)x

// -----------------------------------------------------

namespace mmdb  {

#ifdef  UseDoubleFloat

  typedef  double       realtype;
  const realtype MinReal  = 2.2250e-307;
  const realtype MaxReal  = 1.7976e+308;
  const realtype fMinReal = 2.2250e-307;
  const realtype fMaxReal = 1.7976e+308;

#else

  typedef  float       realtype;
  const realtype MinReal  = 1.1755e-38;
  const realtype MaxReal  = 3.4020e+38;
  const realtype fMinReal = 1.1755e-38;
  const realtype fMaxReal = 3.4020e+38;

#endif

  typedef   float     shortreal;
  const shortreal MinShortReal = 1.1755e-38f;
  const shortreal MaxShortReal = 3.4020e+38f;

/*
#define   strrchr   LastOccurence
#define   fstrrchr  LastOccurence
#define   strchr    FirstOccurence
#define   fstrchr   FirstOccurence
*/

  typedef   char     *         pstr;
  typedef   const char *       cpstr;
  typedef   unsigned int       word;
  typedef   unsigned char      byte;
  typedef   signed   char      short_int;
//  typedef   byte               Boolean;
  typedef   unsigned int       word2;
  typedef   byte *             byteptr;
  typedef   unsigned long      lword;

  typedef   byte intUniBin      [4];
  typedef   byte shortUniBin    [2];
  typedef   byte longUniBin     [4];
  typedef   byte wordUniBin     [4];
  typedef   byte realUniBin     [10];
  typedef   byte floatUniBin    [5];
  typedef   byte shortrealUniBin[5];

#ifdef _WIN32
  pstr strcasestr ( pstr s1, cpstr s2 );
#endif

#ifdef _MSC_VER
#define   strncasecmp _strnicmp
#define   strcasecmp  _stricmp
#endif

  const int      MaxInt  =  32767;
  const int      MinInt  = -32768;
  const word     MaxWord =  65535L;
  const long int MaxInt4 =  2147483647L;

  //    MinInt4 would have to be defined as  -2147483648,
  // however some compilers do not like that. To be on safe,
  // we define it as -2147483647:
  const long int MinInt4  = -2147483647;
  const lword    MaxWord4 =  4294967295UL;

  const realtype Pi   = 3.141592653589793238462643;
  const realtype Eu   = 2.718281828459045235360287;
  const realtype ln10 = 2.3025850929940456840179915;

  // ***  vectors   X[1..N] :
  typedef   realtype * rvector;
  typedef   int      * ivector;
  typedef   word     * wvector;
  typedef   byte     * bvector;
  typedef   bool     * ovector;
  typedef   long     * lvector;
  typedef   lword    * lwvector;
  typedef   pstr     * psvector;

  // ***  matrices   X[1..N][1..M] :
  typedef   rvector  * rmatrix;
  typedef   ivector  * imatrix;
  typedef   wvector  * wmatrix;
  typedef   bvector  * bmatrix;
  typedef   ovector  * omatrix;
  typedef   lvector  * lmatrix;
  typedef   lwvector * lwmatrix;
  typedef   psvector * psmatrix;

  // ***  matrices   X[1..N][1..M][1..K] :
  typedef   rmatrix  * rmatrix3;
  typedef   imatrix  * imatrix3;
  typedef   wmatrix  * wmatrix3;
  typedef   bmatrix  * bmatrix3;
  typedef   omatrix  * omatrix3;
  typedef   lmatrix  * lmatrix3;
  typedef   lwmatrix * lwmatrix3;
  typedef   psmatrix * psmatrix3;


  // ------------------------------------------------------------

  //  Initialization. Some C++ enviroments do not do call
  // InitMatType() automatically, therefore it is always
  // advisable to call InitMatType() explicitely from the top of
  // main(). It is completely harmless and cheap (although
  // unnecessary) to call InitMatType() multiple times.
  extern MMDB_IMEX bool InitMatType();

  // ------------------------------------------------------------

  MMDB_IMEX inline int mround ( realtype X )  { return (int)floor(X+0.5);   }
  MMDB_IMEX inline int ifloor ( realtype X )  { return (int)floor(X);       }
  MMDB_IMEX inline int Abs    ( int x )       { return ( x >= 0 ? x : -x ); }

  MMDB_IMEX inline void ISwap ( int & x, int & y )
  { int  b = x;  x = y;  y = b; }

  MMDB_IMEX inline void WSwap ( word & x, word & y )
  { word  b = x;  x = y;  y = b; }

  MMDB_IMEX inline void BSwap ( byte & x, byte & y )
  { byte b = x;  x = y;  y = b; }

  MMDB_IMEX inline void OSwap ( bool & x, bool & y )
  { bool b = x;  x = y;  y = b; }

  MMDB_IMEX inline void LSwap ( long & x, long & y )
  { long b = x;  x = y;  y = b; }

  MMDB_IMEX inline void RSwap ( realtype & x, realtype & y )
  { realtype b = x;  x = y;  y = b; }

  MMDB_IMEX inline realtype RMax ( const realtype x1, const realtype x2 )
  { return ( x1 > x2 ? x1 : x2 );  }

  MMDB_IMEX inline long LMax ( const long x1, const long x2 )
  { return ( x1 > x2 ? x1 : x2 );  }

  MMDB_IMEX inline word WMax ( const word x1, const word x2 )
  { return ( x1 > x2 ? x1 : x2 );  }

  MMDB_IMEX inline int IMax ( const int x1,  const int x2  )
  { return ( x1 > x2 ? x1 : x2 );  }

  MMDB_IMEX inline realtype RMin ( const realtype x1, const realtype x2 )
  { return ( x1 < x2 ? x1 : x2 );  }

  MMDB_IMEX inline long LMin ( const long x1, const long x2 )
  { return ( x1 < x2 ? x1 : x2 );  }

  MMDB_IMEX inline word WMin ( const word x1, const word x2 )
  { return ( x1 < x2 ? x1 : x2 );  }

  MMDB_IMEX inline int  IMin ( const int x1,  const int x2  )
  { return ( x1 < x2 ? x1 : x2 );  }

  MMDB_IMEX inline realtype fsign ( const realtype x1,  const realtype x2 )  {
  realtype  ax;
    if (x1>=0.0)  ax = x1;
            else  ax = -x1;
    return ( x2 >= 0.0 ? ax : -ax );
  }


  // ------------------------------------------------------------

  //    Allocated vectors are enumerated as [Shift..Shift+N-1]
  //  rather than [0..N-1] !
  //    Get-functions return  <true>  if memory was allocated;
  //  if allocation attemt fails,  vector is assigned with  NULL

  extern bool GetVectorMemory ( rvector  & V, word N, word Shift=1 );
  extern bool GetVectorMemory ( ivector  & I, word N, word Shift=1 );
  extern bool GetVectorMemory ( wvector  & W, word N, word Shift=1 );
  extern bool GetVectorMemory ( bvector  & B, word N, word Shift=1 );
  extern bool GetVectorMemory ( ovector  & O, word N, word Shift=1 );
  extern bool GetVectorMemory ( lvector  & L, word N, word Shift=1 );
  extern bool GetVectorMemory ( lwvector & L, word N, word Shift=1 );
  extern bool GetVectorMemory ( psvector & P, word N, word Shift=1 );

  //    Shift at deallocation MUST be the same as that at allocation !
  //   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //    Free-functions do nothing if vector has value  NULL (e.g.
  //  after unsuccessful allocation).

  extern void FreeVectorMemory ( rvector  & V, word Shift=1 );
  extern void FreeVectorMemory ( ivector  & I, word Shift=1 );
  extern void FreeVectorMemory ( wvector  & W, word Shift=1 );
  extern void FreeVectorMemory ( bvector  & B, word Shift=1 );
  extern void FreeVectorMemory ( ovector  & O, word Shift=1 );
  extern void FreeVectorMemory ( lvector  & L, word Shift=1 );
  extern void FreeVectorMemory ( lwvector & L, word Shift=1 );
  extern void FreeVectorMemory ( psvector & P, word Shift=1 );

  // -------------------------------------------------------------

  //    Allocated matrices are enumerated as
  //          [ShiftN..ShiftN+N-1, ShiftM..ShiftM+M-1]
  //  rather than [0..N-1,0..M-1] !
  //    Get-functions return  <true>  if memory was allocated;
  //  if allocation attemt fails,  matrix is assigned with  NULL
  //    Free-functions do nothing if matrix has value  NULL (e.g.
  //  after unsuccessful allocation).

  extern bool GetMatrixMemory
       ( rmatrix  & A, word N, word M, word ShiftN=1, word ShiftM=1 );
  extern bool GetMatrixMemory
       ( imatrix  & A, word N, word M, word ShiftN=1, word ShiftM=1 );
  extern bool GetMatrixMemory
       ( wmatrix  & W, word N, word M, word ShiftN=1, word ShiftM=1 );
  extern bool GetMatrixMemory
       ( bmatrix  & B, word N, word M, word ShiftN=1, word ShiftM=1 );
  extern bool GetMatrixMemory
       ( omatrix  & O, word N, word M, word ShiftN=1, word ShiftM=1 );
  extern bool GetMatrixMemory
       ( lmatrix  & L, word N, word M, word ShiftN=1, word ShiftM=1 );
  extern bool GetMatrixMemory
       ( lwmatrix & L, word N, word M, word ShiftN=1, word ShiftM=1 );
  extern bool GetMatrixMemory
       ( psmatrix & P, word N, word M, word ShiftN=1, word ShiftM=1 );

  //    ShiftN and ShiftM at deallocation MUST be the same as those at
  //   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //                          allocation !
  //                         ~~~~~~~~~~~~~

  extern  MMDB_IMEX void FreeMatrixMemory  ( rmatrix  & A,  word N,
                   word ShiftN=1, word ShiftM=1 );
  extern  MMDB_IMEX void FreeMatrixMemory  ( imatrix  & A,  word N,
                   word ShiftN=1, word ShiftM=1 );
  extern  MMDB_IMEX void FreeMatrixMemory  ( wmatrix  & W,  word N,
                   word ShiftN=1, word ShiftM=1 );
  extern  MMDB_IMEX void FreeMatrixMemory  ( bmatrix  & B,  word N,
                   word ShiftN=1, word ShiftM=1 );
  extern  MMDB_IMEX void FreeMatrixMemory  ( omatrix  & O,  word N,
                   word ShiftN=1, word ShiftM=1 );
  extern  MMDB_IMEX void FreeMatrixMemory  ( lmatrix  & L,  word N,
                   word ShiftN=1, word ShiftM=1 );
  extern  MMDB_IMEX void FreeMatrixMemory  ( lwmatrix & L,  word N,
                   word ShiftN=1, word ShiftM=1 );
  extern  MMDB_IMEX void FreeMatrixMemory  ( psmatrix & P,  word N,
                   word ShiftN=1, word ShiftM=1 );


  // -------------------------------------------------------------
  //   3D matrices

  extern bool GetMatrix3Memory
             ( rmatrix3  & A, word N, word M, word K,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );
  extern bool GetMatrix3Memory
             ( imatrix3  & A, word N, word M, word K,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );
  extern bool GetMatrix3Memory
             ( wmatrix3  & A, word N, word M, word K,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );
  extern bool GetMatrix3Memory
             ( bmatrix3  & A, word N, word M, word K,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );
  extern bool GetMatrix3Memory
             ( omatrix3  & A, word N, word M, word K,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );
  extern bool GetMatrix3Memory
             ( lmatrix3  & A, word N, word M, word K,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );
  extern bool GetMatrix3Memory
             ( lwmatrix3 & A, word N, word M, word K,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );
  extern bool GetMatrix3Memory
             ( psmatrix3 & A, word N, word M, word K,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );

  //
  //    ShiftN, ShiftM and ShiftK at deallocation MUST be
  //   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //         the same as those at allocation !
  //        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  extern void FreeMatrix3Memory
                 ( rmatrix3  & A, word N, word M,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );
  extern void FreeMatrix3Memory
                 ( imatrix3  & A, word N, word M,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );
  extern void FreeMatrix3Memory
                 ( wmatrix3  & A, word N, word M,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );
  extern void FreeMatrix3Memory
                 ( bmatrix3  & A, word N, word M,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );
  extern void FreeMatrix3Memory
                 ( omatrix3  & A, word N, word M,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );
  extern void FreeMatrix3Memory
                 ( lmatrix3  & A, word N, word M,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );
  extern void FreeMatrix3Memory
                 ( lwmatrix3 & A, word N, word M,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );
  extern void FreeMatrix3Memory
                 ( psmatrix3 & A, word N, word M,
                   word ShiftN=1, word ShiftM=1, word ShiftK=1 );

  // -------------------------------------------------------------

  extern  realtype  MachEps;
  extern  realtype  floatMachEps;
  extern  realtype  LnMaxReal;
  extern  realtype  LnMinReal;

  extern  MMDB_IMEX realtype  MachinEps     ();
  extern  MMDB_IMEX realtype  floatMachinEps();
  extern  MMDB_IMEX realtype  frac  ( realtype R   );
  extern  MMDB_IMEX long      mod   ( long     x, long     y );
  extern  MMDB_IMEX realtype  Pow   ( realtype X, int      y );
  extern  MMDB_IMEX realtype  Pow1  ( realtype X, realtype Y );
  extern  MMDB_IMEX realtype  Exp   ( realtype X );   // use to avoid catastrophies
  extern  MMDB_IMEX bool      Odd   ( int      i );
  extern  MMDB_IMEX long    HexValL ( cpstr S );
  extern  MMDB_IMEX long    OctValL ( cpstr S );
  extern  MMDB_IMEX long    BinValL ( cpstr S );
  extern  MMDB_IMEX pstr    BinValS ( long L, pstr S  );  // S[sizeof(long)+1] at least

  extern  MMDB_IMEX pstr ParamStr ( pstr D, cpstr S, realtype V, int M=5,
                          cpstr S1=(pstr)"" );
  extern  MMDB_IMEX pstr ParamStr ( pstr D, cpstr S, realtype V, int M,
                          cpstr S1, realtype V2, int M2=5,
                          cpstr S2=(pstr)"" );


  //  ----------  Strings

  //   CreateCopy(..) allocates Dest string and copies the contents of
  // Source into it. If Dest is not NULL prior calling the function,
  // it is attempted to deallocate first.

  extern pstr CreateCopy     ( pstr & Dest, cpstr Source );
  extern pstr CreateCopy_n   ( pstr & Dest, cpstr Source, int n );

  extern pstr CreateConcat   ( pstr & Dest, cpstr Source );
  extern pstr CreateConcat   ( pstr & Dest, cpstr Source1,
                                          cpstr Source2 );
  extern pstr CreateConcat   ( pstr & Dest, cpstr Source1,
                                            cpstr Source2,
                                            cpstr Source3 );
  extern pstr CreateConcat   ( pstr & Dest, cpstr Source1,
                                            cpstr Source2,
                                            cpstr Source3,
                                            cpstr Source4 );
  extern pstr CreateConcat   ( pstr & Dest, cpstr Source1,
                                            cpstr Source2,
                                            cpstr Source3,
                                            cpstr Source4,
                                            cpstr Source5 );

  extern pstr CreateCopCat   ( pstr & Dest, cpstr Source1,
                                            cpstr Source2,
                                            cpstr Source3,
                                            cpstr Source4,
                                            cpstr Source5 );
  extern pstr CreateCopCat   ( pstr & Dest, cpstr Source1,
                                            cpstr Source2,
                                            cpstr Source3,
                                            cpstr Source4 );
  extern pstr CreateCopCat   ( pstr & Dest, cpstr Source1,
                                            cpstr Source2,
                                            cpstr Source3 );
  extern pstr CreateCopCat   ( pstr & Dest, cpstr Source1,
                                            cpstr Source2 );

  extern pstr LastOccurence  ( cpstr S     , char c      );
  extern pstr FirstOccurence ( cpstr S     , char c      );
  extern int  indexOf        ( cpstr S     , char c      );
  extern pstr FirstOccurence ( cpstr S, int Slen,
                               cpstr Q, int Qlen );
  extern int  indexOf        ( cpstr S, int Slen,
                               cpstr Q, int Qlen );

  extern pstr LowerCase ( pstr s );
  extern pstr UpperCase ( pstr s );

  //   GetString(..) copies first M characters of string S into string
  // L, appending the terminating null. If S contains less then M
  // characters, L will be padded with spaces.
  extern void GetString ( pstr L, cpstr S, int M );

  //   GetStrTer(..) copies at least n (or LMax if LMax<n) first symbols
  // of string S into string L, then continues copying until first space
  // or terminating null is found. If the terminating null is met among
  // the first n characters or if SMax<n, the string L will be padded
  // with spaces till the length of minimum of n and LMax and then
  // terminated with the null.
  //   LMax and SMax are the buffer lengths of L and S,
  // respectively. Even if no space is found, the last character
  // in L will be the terminating null.
  extern void GetStrTer ( pstr L, cpstr S, int n, int LMax,
                          int SMax );


  // Version of GetStrTer(..) allowing for spaces in the string.
  //
  //   Copies at least n (or LMax if LMax<n) first symbols of
  // string S into string L, then continues copying until first
  // terminating null is found. If the terminating null
  // is met among the first n characters or if SMax<n, the string
  // L will be padded with spaces till the length of minimum of
  // n and LMax and then terminated with the null.
  //   SMax are buffer lengths of L and S, respectively. The last
  // character in L will be the terminating null.
  extern void GetStrTerWin32File ( pstr L, cpstr S, int n,
                                   int LMax, int SMax );


  //   strcpy_n(..) copies at most n symbols from string s to d,
  // but no more than strlen(s) (s must contain a terminating
  // null). The terminating null IS NEITHER appended OR copied
  // to d.
  extern void strcpy_n  ( pstr d, cpstr s, int n );

  //   strcpy_n1(..) copies at most n last symbols from string s
  // to d, but no more than strlen(s) (s must contain a terminating
  // null). The string in d is aligned to the right and added with
  // spaces at the left, if necessary. The terminating null
  // IS NEITHER appended OR copied to d.
  extern void strcpy_n1 ( pstr d, cpstr s, int n );

  //   Copies at most n symbols from string s to d, but no
  // more than strlen(s) (s must contain a terminating null).
  // The string in d is aligned to the right and added with
  // spaces at the left, if necessary. The terminating null
  // IS NEITHER appended NOR copied to d.
  extern void strcpy_nr ( pstr d, cpstr s, int n );

  //   strcpy_ns(..) copies at most n symbols from string s to d,
  // but no more than strlen(s) (s must contain a terminating
  // null). The terminating null IS NEITHER appended NOR copied
  // to d; rather, d is padded with spaces up to the length of n
  // if strlen(s)<n.
  extern void strcpy_ns ( pstr d, cpstr s, int n );

  //   strcpy_cs(..) copies string s to string d cutting all
  // spaces at the end. Thus, " abcde   " will be copied
  // like " abcde" (terminating null appended).
  //   The function returns d.
  extern pstr strcpy_cs ( pstr d, cpstr s );

  //   strcpy_ncs(..) copies at most n characters from string s
  // to string d cutting all spaces at at the end. Thus, " abcde   "
  // will be copied like " abc" at n=4 and like " abcde" at n>5
  // (terminating null appended).
  //   The function returns d.
  extern pstr strcpy_ncs ( pstr d, cpstr s, int n );

  //   strcpy_css(..) copies string s to string d cutting all
  // spaces at the begining and at the end. Thus, " ab c de  "
  // will be copied like "ab c de" (terminating null appended).
  //   The function returns d.
  extern pstr strcpy_css ( pstr d, cpstr s );

  //   strcpy_ncss(..) copies at most n characters from string s
  // to string d cutting all spaces at the begining and at the end.
  // Thus, " ab c de  " will be copied like "ab" at n=3 (terminating
  // null appended).
  //   The function returns d.
  extern pstr strcpy_ncss ( pstr d, cpstr s, int n );

  //   strcpy_n0(..) copies at most n symbols from string s to d,
  // but no more than strlen(s) (s must contain a terminating
  // null). The terminating null IS appended to d.
  //   The function returns d.
  extern pstr strcpy_n0 ( pstr d, cpstr s, int n );

  //   strlen_des returns the length of a string as if all extra
  // spaces from the latter have been deleted. Extra spaces
  // include all leading and tracing spaces and any sequential
  // spaces when more than one. The string does not change.
  extern int strlen_des ( cpstr s );

  //   strcpy_des copies string s into string d removing all extra
  // spaces from the latter. Extra spaces include all leading and
  // tracing spaces and any sequential spaces when more than one.
  extern pstr strcpy_des ( pstr d, cpstr s );

  //   strcat_des appends string s to string d removing all extra
  // spaces from the latter. Extra spaces include all leading and
  // tracing spaces and any sequential spaces when more than one.
  extern pstr strcat_des ( pstr d, cpstr s );

  //  PadSpaces(..) pads string S with spaces making its length
  // equal to len. The terminating zero is added, so that S should
  // reserve space of a minimum len+1 characters.
  extern void PadSpaces ( pstr S, int len );

  enum SCUTKEY  {
    SCUTKEY_BEGIN  = 0x00000001,
    SCUTKEY_END    = 0x00000002,
    SCUTKEY_BEGEND = 0x00000003
  };

  //   CutSpaces(..) cuts spaces at the begining or end of
  // string S according to the value of CutKey. The function
  // returns S.
  extern pstr CutSpaces ( pstr S, int CutKey );

  //   DelSpaces(..) removes all spaces (or other symbols as
  // specified by 'c') from the string. The string is then
  // shrinked by the number of removed characters. Thus,
  // " as ttt  " becomes "asttt".
  extern pstr DelSpaces ( pstr S, char c=' ' );

  //   EnforceSpaces(..) replaces all unprintable characters,
  // except <CR>, <LF>, <TAB> and some others, for spaces
  extern pstr EnforceSpaces ( pstr S );

  // -------------------------------------------------------------

  ///   This call will produce correct floats in universal binaries but
  /// make them incompatible with old files. Without this call, float
  /// read/write will result in error after 6th digit.
  ///   UniBin read/write of other types (realtype, shortreal, int etc)
  /// is not affected by this call, and to the best of knowledge is
  /// correct (no loss of precision).
  extern void  set_new_float_unibin();
  extern bool  is_new_float_unibin();
  extern void  set_old_float_unibin();

  extern void __modify4();

  extern void int2UniBin       ( int       I,  intUniBin        iUB  );
  extern void short2UniBin     ( short     S,  shortUniBin      sUB  );
  extern void long2UniBin      ( long      L,  longUniBin       lUB  );
  extern void word2UniBin      ( word      W,  wordUniBin       wUB  );
  extern void real2UniBin      ( realtype  R,  realUniBin       rUB  );
  extern void float2UniBin     ( realtype  R,  floatUniBin      fUB  );
  extern void shortreal2UniBin ( shortreal R,  shortrealUniBin  srUB );
  extern void UniBin2int       ( intUniBin        iUB, int       & I );
  extern void UniBin2short     ( shortUniBin      sUB, short     & S );
  extern void UniBin2long      ( longUniBin       lUB, long      & L );
  extern void UniBin2word      ( wordUniBin       wUB, word      & W );
  extern void UniBin2real      ( realUniBin       rUB, realtype  & R );
  extern void UniBin2shortreal ( shortrealUniBin srUB, shortreal & R );
  extern void UniBin2float     ( floatUniBin      fUB, realtype  & R );

  extern void mem_write ( int       I, pstr S, int & l );
  extern void mem_write ( short     I, pstr S, int & l );
  extern void mem_write ( long      I, pstr S, int & l );
  extern void mem_write ( word      W, pstr S, int & l );
  extern void mem_write ( realtype  R, pstr S, int & l );
  extern void mem_write ( shortreal R, pstr S, int & l );
  extern void mem_write ( pstr      L, int len, pstr S, int & l );
  extern void mem_write ( pstr      L, pstr S, int & l );
  extern void mem_write ( bool      B, pstr S, int & l );
  extern void mem_write_byte ( byte B, pstr S, int & l );

  extern void mem_read  ( int       & I, cpstr S, int & l );
  extern void mem_read  ( short     & I, cpstr S, int & l );
  extern void mem_read  ( long      & I, cpstr S, int & l );
  extern void mem_read  ( word      & W, cpstr S, int & l );
  extern void mem_read  ( realtype  & R, cpstr S, int & l );
  extern void mem_read  ( shortreal & R, cpstr S, int & l );
  extern void mem_read  ( pstr        L, int len, cpstr S, int & l );
  extern void mem_read  ( pstr      & L, cpstr S, int & l );
  extern void mem_read  ( bool      & B, cpstr S, int & l );
  extern void mem_read_byte ( byte  & B, cpstr S, int & l );

}

#endif

/* ===================================================  */

