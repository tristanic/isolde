//  $Id: mmdb_utils.cpp $
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
//    27.02.17   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :   MMDBF_Utils <implementation>
//       ~~~~~~~~~
//  **** Project :   MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//
//  **** Classes :   mmdb::ContainerClass ( containered class template )
//       ~~~~~~~~~   mmdb::ContString     ( containered string         )
//                   mmdb::ClassContainer ( container of classes       )
//                   mmdb::AtomPath       ( atom path ID               )
//                   mmdb::QuickSort      ( quick sort of integers     )
//
//  **** Functions : Date9to11  ( DD-MMM-YY   -> DD-MMM-YYYY          )
//       ~~~~~~~~~~~ Date11to9  ( DD-MMM-YYYY -> DD-MMM-YY            )
//                   Date9toCIF ( DD-MMM-YY   -> YYYY-MM-DD           )
//                   Date11toCIF( DD-MMM-YYYY -> YYYY-MM-DD           )
//                   DateCIFto9 ( YYYY-MM-DD  -> DD-MMM-YY            )
//                   DateCIFto11( YYYY-MM-DD  -> DD-MMM-YYYY          )
//                   GetInteger ( reads integer from a string         )
//                   GetReal    ( reads real from a string            )
//                   GetIntIns  ( reads integer and insert code       )
//                   PutInteger ( writes integer into a string        )
//                   PutRealF   ( writes real in F-form into a string )
//                   PutIntIns  ( writes integer and insert code      )
//                   CIFGetInteger ( reads and deletes int from CIF   )
//                   CIFGetReal    ( reads and deletes real from CIF  )
//                   CIFGetString  ( reads and deletes string from CIF)
//                   CIFGetInteger1 (reads and del-s int from CIF loop)
//                   CIFGetReal1    (reads and del-s int from CIF loop)
//                   Mat4Inverse    ( inversion of 4x4 matrices       )
//                   GetErrorDescription (ascii line to an Error_XXXXX)
//                   ParseAtomID    ( parses atom ID line             )
//                   ParseResID     ( parses residue ID line          )
//                   ParseAtomPath  ( parses full atom path           )
//
//   (C) E. Krissinel  2000-2017
//
//  =================================================================
//

#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "mmdb_utils.h"
#include "hybrid_36.h"

namespace mmdb  {

  // ====================== Date functions  =======================

  static cpstr Month[12] = {
    "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
    "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
  };

  static cpstr nMonth[12] = {
    "01", "02", "03", "04", "05", "06",
    "07", "08", "09", "10", "11", "12"
  };

  void  Date9to11 ( cpstr Date9, pstr Date11 )  {
  // converts  DD-MMM-YY to DD-MMM-YYYY
  int i;
    i = 0;
    while ((i<12) && (strncmp(Month[i],&(Date9[3]),3)))  i++;
    if (i<12)  {   // DD-MMM-YY -> DD-MMM-YYYY
      strncpy ( Date11,Date9,7 );
      if (Date9[7]!='0')  strncpy ( &(Date11[7]),"19",2 );
                    else  strncpy ( &(Date11[7]),"20",2 );
      strncpy ( &(Date11[9]),&(Date9[7]),2 );
      Date11[2]  = '-';
      Date11[6]  = '-';
      Date11[11] = char(0);
    } else  {     // DD-MM-YY -> DD-MMM-YYYY
      strncpy ( Date11,Date9,3 );
      i = 0;
      while ((i<12) && (strncmp(nMonth[i],&(Date9[3]),2)))  i++;
      if (i<12)  {
        strncpy ( &(Date11[3]),Month[i],3 );
          if (Date9[6]!='0')  strncpy ( &(Date11[7]),"19",2 );
                    else  strncpy ( &(Date11[7]),"20",2 );
        strncpy ( &(Date11[9]),&(Date9[6]),2 );
        Date11[2]  = '-';
        Date11[6]  = '-';
        Date11[11] = char(0);
      } else
        strcpy ( Date11,"           " );
    }
  }

  void  Date11to9 ( cpstr Date11, pstr Date9 )  {
  // converts DD-MMM-YYYY to DD-MMM-YY
  int i;
    i = 0;
    while ((i<12) && (strncmp(Month[i],&(Date11[3]),3)))  i++;
    if (i<12)  {   // DD-MMM-YYYY -> DD-MMM-YY
      strncpy ( Date9,Date11,7 );
      strncpy ( &(Date9[7]),&(Date11[9]),2 );
      Date9[2] = '-';
      Date9[6] = '-';
    } else  {      // DD-MM-YYYY -> DD-MMM-YY
      strncpy ( Date9,Date11,3 );
      i = 0;
      while ((i<12) && (strncmp(nMonth[i],&(Date11[3]),2)))  i++;
      if (i<12)  {
        strncpy ( &(Date9[3]),Month[i],3 );
        strncpy ( &(Date9[7]),&(Date11[8]),2 );
        Date9[2] = '-';
        Date9[6] = '-';
      } else
        strcpy ( Date9,"         " );
    }
  }

  void  Date9toCIF ( cpstr Date9, pstr DateCIF )  {
  //  DD-MMM-YY -> YYYY-MM-DD             )
  int  i;
    i = 0;
    while ((i<12) && (strncmp(Month[i],&(Date9[3]),3)))  i++;
    if (i<12)  {   //  DD-MMM-YY  -> YYYY-MM-DD
      if (Date9[7]!='0')  strcpy ( DateCIF,"19" );
                    else  strcpy ( DateCIF,"20" );
      strncpy ( &(DateCIF[2]),&(Date9[7]),2 );
      strncpy ( &(DateCIF[5]),nMonth[i],2 );
    } else  {      //  DD-MM-YY  ->  YYYY-MM-DD
      if (Date9[6]!='0')  strcpy ( DateCIF,"19" );
                    else  strcpy ( DateCIF,"20" );
      strncpy ( &(DateCIF[2]),&(Date9[6]),2 );
      strncpy ( &(DateCIF[5]),&(Date9[3]),2 );
    }
    DateCIF[4] = '-';
    DateCIF[7] = '-';
    strncpy ( &(DateCIF[8]),Date9,2 );
    DateCIF[10] = char(0);
  }

  void  Date11toCIF ( cpstr Date11, pstr DateCIF )  {
  //  DD-MMM-YYYY -> YYYY-MM-DD
  int  i;
    i = 0;
    while ((i<12) && (strncmp(Month[i],&(Date11[3]),3)))  i++;
    if (i<12) {
      strncpy ( DateCIF,&(Date11[7]),4 );
      strncpy ( &(DateCIF[5]),nMonth[i],2 );
    } else  {
      strncpy ( DateCIF,&(Date11[6]),4 );
      strncpy ( &(DateCIF[5]),&(Date11[3]),2 );
    }
    DateCIF[4] = '-';
    DateCIF[7] = '-';
    strncpy ( &(DateCIF[8]),Date11,2 );
    DateCIF[10] = char(0);
  }

  void  DateCIFto9 ( cpstr DateCIF, pstr Date9 )  {
  //  YYYY-MM-DD -> DD-MMM-YY
  int  i;
    strncpy ( Date9,&(DateCIF[8]),2 );
    Date9[2] = '-';
    i = 0;
    while ((i<12) && (strncmp(nMonth[i],&(DateCIF[5]),2)))  i++;
    if (i<12) strncpy ( &(Date9[3]),Month[i],3 );
    else  {
      strncpy ( &(Date9[3]),&(DateCIF[5]),2 );
      Date9[5] = 'X';
    }
    Date9[6] = '-';
    strncpy ( &(Date9[7]),&(DateCIF[2]),2 );
  //  DateCIF[9] = char(0);
  }

  void  DateCIFto11 ( cpstr DateCIF, pstr Date11 )  {
  //  YYYY-MM-DD  -> DD-MMM-YYYY
  int  i;
    strncpy ( Date11,&(DateCIF[8]),2 );
    Date11[2] = '-';
    i = 0;
    while ((i<12) && (strncmp(nMonth[i],&(DateCIF[5]),2)))  i++;
    if (i<12) strncpy ( &(Date11[3]),Month[i],3 );
    else  {
      strncpy ( &(Date11[3]),&(DateCIF[5]),2 );
      Date11[5] = 'X';
    }
    Date11[6] = '-';
    strncpy ( &(Date11[7]),DateCIF,4 );
  //  DateCIF[11] = char(0);
  }


  //  =============== Format functions  ===================

  bool GetInteger ( int & N, cpstr S, int M )  {
  //   Returns true if S contains an integer number in its
  // first M characters. This number is returned in N.
  //   The return is false if no integer number may be
  // recognized. In this case, N is assigned MinInt4 value.
  pstr endptr;
  char L[50];
    strncpy ( L,S,M );
    L[M] = char(0);
    N    = mround(strtod(L,&endptr));
    if ((N==0) && (endptr==L))  {
      N = MinInt4;  // no number
      return false;
    } else
      return true;
  }

  bool GetReal ( realtype & R, cpstr S, int M )  {
  //   Returns true if S contains a real number in its
  // first M characters. This number is returned in R.
  //   The return is false if no real number may be
  // recognized. In this case, R is assigned -MaxReal value.
  pstr endptr;
  char L[50];
    strncpy ( L,S,M );
    L[M] = char(0);
    R    = strtod(L,&endptr);
    if ((R==0.0) && (endptr==L))  {
      R = -MaxReal;  // no number
      return false;
    } else
      return true;
  }

  bool  GetIntIns ( int & N, pstr ins, cpstr S, int M )  {
  //   Returns true if S contains an integer number in its
  // first M characters. This number is returned in N. In addition
  // to that, GetIntIns() retrieves the insertion code which may
  // follow the integer and returns it in "ins" (1 character +
  // terminating 0).
  //   The return is false if no integer number may be
  // recognized. In this case, N is assigned MinInt4 value,
  // "ins" just returns (M+1)th symbol of S (+terminating 0).
  pstr endptr;
  char L[50];

    if (S[M]!=' ')  {
      ins[0] = S[M];
      ins[1] = char(0);
    } else
      ins[0] = char(0);

    strncpy ( L,S,M );
    L[M] = char(0);
    if ((M==4) && ((S[0]>='A') || ((S[0]=='-') && (S[1]>='A'))))
      hy36decode ( M,L,M,&N);
    else  {
      endptr = NULL;
      N      = mround(strtod(L,&endptr));
      if ((N==0) && (endptr==L))  {
        N = MinInt4;  // no number
        return false;
      }
    }

    return true;

  }

  void  PutInteger ( pstr S, int N, int M )  {
  //  Integer N is converted into ASCII string of length M
  // and pasted onto first M characters of string S. No
  // terminating zero is added.
  //  If N is set to MinInt4, then first M characters of
  // string S are set to the space character.
  char L[50];
  int  i;
    if (N==MinInt4)
      for (i=0;i<M;i++)
        S[i] = ' ';
    else  {
      sprintf ( L,"%*i",M,N );
      strncpy ( S,L,M );
    }
  }

  void  PutRealF ( pstr S, realtype R, int M, int L )  {
  //  Real R is converted into ASCII string of length M
  // and pasted onto first M characters of string S. No
  // terminating zero is added. The conversion is done
  // according to fixed format FM.L
  //  If R is set to -MaxReal, then first M characters of
  // string S are set to the space character.
  char N[50];
  int  i;
    if (R==-MaxReal)
      for (i=0;i<M;i++)
        S[i] = ' ';
    else  {
      sprintf ( N,"%*.*f",M,L,R );
      strncpy ( S,N,M );
    }
  }

  ERROR_CODE CIFGetIntegerD ( int & I, mmcif::PLoop Loop, cpstr Tag,
                             int defValue )  {
  int        Signal;
  ERROR_CODE RC;
    Signal = 0;
    RC = CIFGetInteger ( I,Loop,Tag,Signal );
    if (RC)
      I = defValue;
    return RC;
  }

  ERROR_CODE CIFGetInteger ( int & I, mmcif::PLoop Loop, cpstr Tag,
                             int & Signal )  {
  pstr F;
  int  RC;
    RC = Loop->GetInteger ( I,Tag,Signal,true );
    if (RC==mmcif::CIFRC_WrongFormat)  {
      F = Loop->GetString ( Tag,Signal,RC );
      if (F) sprintf ( CIFErrorLocation,"loop %s.%s row %i data %s",
                       Loop->GetCategoryName(),Tag,Signal,F );
        else sprintf ( CIFErrorLocation,"loop %s.%s row %i data [NULL]",
                       Loop->GetCategoryName(),Tag,Signal );
      Signal = -Error_UnrecognizedInteger-1;
      return Error_UnrecognizedInteger;
    }
    if (RC==mmcif::CIFRC_WrongIndex)  {
      Signal = -1;
      return Error_NoData;
    }
    if (RC)  {
      F = Loop->GetString ( Tag,Signal,RC );
      if (F) sprintf ( CIFErrorLocation,"loop %s.%s row %i data %s",
                       Loop->GetCategoryName(),Tag,Signal,F );
        else sprintf ( CIFErrorLocation,"loop %s.%s row %i data [NULL]",
                       Loop->GetCategoryName(),Tag,Signal );
      Signal = -Error_NoData-1;
      return Error_NoData;
    }
    return Error_NoError;
  }


  ERROR_CODE CIFGetInteger1 ( int & I, mmcif::PLoop Loop, cpstr Tag,
                              int nrow )  {
  pstr F;
  int  RC;
    RC = Loop->GetInteger ( I,Tag,nrow,true );
    if (RC==mmcif::CIFRC_WrongFormat)  {
      F = Loop->GetString ( Tag,nrow,RC );
      if (F) sprintf ( CIFErrorLocation,"loop %s.%s row %i data %s",
                       Loop->GetCategoryName(),Tag,nrow,F );
        else sprintf ( CIFErrorLocation,"loop %s.%s row %i data [NULL]",
                       Loop->GetCategoryName(),Tag,nrow );
      return Error_UnrecognizedInteger;
    }
    if (RC==mmcif::CIFRC_WrongIndex)
      return Error_NoData;
    if (RC)  {
      F = Loop->GetString ( Tag,nrow,RC );
      if (F) sprintf ( CIFErrorLocation,"loop %s.%s row %i data %s",
                       Loop->GetCategoryName(),Tag,nrow,F );
        else sprintf ( CIFErrorLocation,"loop %s.%s row %i data [NULL]",
                       Loop->GetCategoryName(),Tag,nrow );
      return Error_NoData;
    }
    return Error_NoError;
  }


  ERROR_CODE CIFGetReal ( realtype & R, mmcif::PLoop Loop, cpstr Tag,
                          int & Signal )  {
  pstr F;
  int  RC;
    RC = Loop->GetReal ( R,Tag,Signal,true );
    if (RC==mmcif::CIFRC_WrongFormat)  {
      F = Loop->GetString ( Tag,Signal,RC );
      if (F) sprintf ( CIFErrorLocation,"loop %s.%s row %i data %s",
                       Loop->GetCategoryName(),Tag,Signal,F );
        else sprintf ( CIFErrorLocation,"loop %s.%s row %i data [NULL]",
                       Loop->GetCategoryName(),Tag,Signal );
      Signal = -Error_UnrecognizedReal-1;
      return Error_UnrecognizedReal;
    }
    if (RC==mmcif::CIFRC_WrongIndex)  {
      Signal = -1;
      return Error_NoData;
    }
    if (RC)  {
      F = Loop->GetString ( Tag,Signal,RC );
      if (F) sprintf ( CIFErrorLocation,"loop %s.%s row %i data %s",
                       Loop->GetCategoryName(),Tag,Signal,F );
        else sprintf ( CIFErrorLocation,"loop %s.%s row %i data [NULL]",
                       Loop->GetCategoryName(),Tag,Signal );
      Signal = -Error_NoData-1;
      return Error_NoData;
    }
    return Error_NoError;
  }


  ERROR_CODE CIFGetReal1 ( realtype & R, mmcif::PLoop Loop, cpstr Tag,
                           int nrow )  {
  pstr F;
  int  RC;
    RC = Loop->GetReal ( R,Tag,nrow,true );
    if (RC==mmcif::CIFRC_WrongFormat)  {
      F = Loop->GetString ( Tag,nrow,RC );
      if (F) sprintf ( CIFErrorLocation,"loop %s.%s row %i data %s",
                       Loop->GetCategoryName(),Tag,nrow,F );
        else sprintf ( CIFErrorLocation,"loop %s.%s row %i data [NULL]",
                       Loop->GetCategoryName(),Tag,nrow );
      return Error_UnrecognizedReal;
    }
    if (RC==mmcif::CIFRC_WrongIndex)
      return Error_NoData;
    if (RC)  {
      F = Loop->GetString ( Tag,nrow,RC );
      if (F) sprintf ( CIFErrorLocation,"loop %s.%s row %i data %s",
                       Loop->GetCategoryName(),Tag,nrow,F );
        else sprintf ( CIFErrorLocation,"loop %s.%s row %i data [NULL]",
                       Loop->GetCategoryName(),Tag,nrow );
      return Error_NoData;
    }
    return Error_NoError;
  }


  ERROR_CODE CIFGetString ( pstr S, mmcif::PLoop Loop, cpstr Tag,
                            int row, int SLen, cpstr DefS )  {
  pstr F;
  int RC;
    F = Loop->GetString ( Tag,row,RC );
    if ((!RC) && F)  {
      strncpy ( S,F,SLen-1 );
      Loop->DeleteField ( Tag,row );
      return Error_NoError;
    } else  {
      strcpy ( S,DefS );
      return Error_EmptyCIFLoop;
    }
  }


  ERROR_CODE CIFGetInteger ( int & I, mmcif::PStruct Struct, cpstr Tag,
                             bool Remove )  {
  pstr F;
  int  RC;
    RC = Struct->GetInteger ( I,Tag,Remove );
    if (RC==mmcif::CIFRC_WrongFormat)  {
      F = Struct->GetString ( Tag,RC );
      if (F) sprintf ( CIFErrorLocation,"structure %s.%s data %s",
                       Struct->GetCategoryName(),Tag,F );
        else sprintf ( CIFErrorLocation,"structure %s.%s data [NULL]",
                       Struct->GetCategoryName(),Tag );
      return Error_UnrecognizedInteger;
    }
    if (RC)  {
      F = Struct->GetString ( Tag,RC );
      if (F) sprintf ( CIFErrorLocation,"structure %s.%s data %s",
                       Struct->GetCategoryName(),Tag,F );
        else sprintf ( CIFErrorLocation,"structure %s.%s data [NULL]",
                       Struct->GetCategoryName(),Tag );
      return Error_NoData;
    }
    return Error_NoError;
  }

  ERROR_CODE CIFGetReal ( realtype & R, mmcif::PStruct Struct, cpstr Tag,
                          bool Remove )  {
  pstr F;
  int RC;
    RC = Struct->GetReal ( R,Tag,Remove );
    if (RC==mmcif::CIFRC_WrongFormat)  {
      F = Struct->GetString ( Tag,RC );
      if (F) sprintf ( CIFErrorLocation,"structure %s.%s data %s",
                       Struct->GetCategoryName(),Tag,F );
        else sprintf ( CIFErrorLocation,"structure %s.%s data [NULL]",
                       Struct->GetCategoryName(),Tag );
      return Error_UnrecognizedReal;
    }
    if (RC)  {
      F = Struct->GetString ( Tag,RC );
      if (F) sprintf ( CIFErrorLocation,"structure %s.%s data %s",
                       Struct->GetCategoryName(),Tag,F );
        else sprintf ( CIFErrorLocation,"structure %s.%s data [NULL]",
                       Struct->GetCategoryName(),Tag );
      return Error_NoData;
    }
    return Error_NoError;
  }

  ERROR_CODE CIFGetString ( pstr S, mmcif::PStruct Struct, cpstr Tag,
                            int SLen, cpstr DefS, bool Remove )  {
  pstr F;
  int  RC;
    F = Struct->GetString ( Tag,RC );
    if ((!RC) && F)  {
      strcpy_n0 ( S,F,SLen-1 );
      if (Remove)  Struct->DeleteField ( Tag );
      return Error_NoError;
    } else  {
      strcpy ( S,DefS );
      return Error_EmptyCIFStruct;
    }
  }


  void  PutIntIns ( pstr S, int N, int M, cpstr ins )  {
  //  Integer N is converted into ASCII string of length M
  // and pasted onto first M characters of string S. No
  // terminating zero is added. The insert code ins is put
  // immediately after the integer.
  //  If N is set to MinInt4, then first M+1 characters of
  // string S are set to space, and no insert code are
  // appended.
  char L[50];
  int  i;

    if (N==MinInt4)  {
      for (i=0;i<=M;i++)
        S[i] = ' ';
    } else  {
      if ((M!=4) || ((N>=-999) && (N<=9999)))
           sprintf    ( L,"%*i",M,N );
      else hy36encode ( M,N,L );
      strcpy_n1 ( S,L,M );
      if (ins[0]) S[M] = ins[0];
    }

  }


  void  Mat4Inverse ( const mat44 & A, mat44 & AI )  {
  //       ***  FORMER RBRINV(A,AI)  ***
  //   Function to invert 4*4 matrices (AI=A^{-1})
  mat44    c;
  mat33    x;
  realtype s,s1;
  int      ii,jj,i,i1,j,j1;

  // ---- Get cofactors of 'a' in array 'c'

    s1 = 1.0;
    for (ii=0;ii<4;ii++)  {
      s = s1;
      for (jj=0;jj<4;jj++)  {
        i = -1;
        for (i1=0;i1<4;i1++)
          if (i1!=ii)  {
            i++;
            j = -1;
            for (j1=0;j1<4;j1++)
              if (j1!=jj)  {
                j++;
                x[i][j] = A[i1][j1];
              }
          }
        c[ii][jj] = s*(x[0][0]*(x[1][1]*x[2][2]-x[1][2]*x[2][1]) +
                       x[0][1]*(x[1][2]*x[2][0]-x[1][0]*x[2][2]) +
                       x[0][2]*(x[1][0]*x[2][1]-x[1][1]*x[2][0]));
        s = -s;
      }
      s1 = -s1;
    }

  // ---- Calculate determinant

    s = 0.0;
    for (i=0;i<4;i++)
      s += A[i][0]*c[i][0];

  // ---- Get inverse matrix

    if (s!=0.0)
      for (i=0;i<4;i++)
        for (j=0;j<4;j++)
          AI[i][j] = c[j][i]/s;

  }

  realtype Mat3Inverse ( const mat33 & A, mat33 & AI )  {
  mat33    c,x;
  realtype s;
  int      ii,jj,i,i1,j,j1;

    // Get cofactors of 'a' in array 'c'

    s = 1.0;
    for (ii=0;ii<3;ii++)
      for (jj=0;jj<3;jj++)  {
        i = -1;
        for (i1=0;i1<3;i1++)
          if (i1!=ii)  {
            i++;
            j = -1;
            for (j1=0;j1<3;j1++)
              if (j1!=jj)  {
                j++;
                x[i][j] = A[i1][j1];
              }
          }
        c[ii][jj] = s*(x[0][0]*x[1][1]-x[0][1]*x[1][0]);
        s = -s;
      }

    // Calculate determinant

    s = 0.0;
    for (i=0;i<3;i++)
      s += A[i][0]*c[i][0];

    // Get inverse matrix

    if (s!=0.0)
      for (i=0;i<3;i++)
        for (j=0;j<3;j++)
          AI[i][j] = c[j][i]/s;

    return s;

  }

  void  Mat4Mult ( mat44 & A, const mat44 & B, const mat44 & C )  {
  //  Calculates A=B*C
  int i,j,k;
    for (i=0;i<4;i++)
      for (j=0;j<4;j++)  {
        A[i][j] = 0.0;
        for (k=0;k<4;k++)
          A[i][j] += B[i][k]*C[k][j];
      }
  }

  void  Mat4Div1 ( mat44 & A, const mat44 & B, const mat44 & C )  {
  //  Calculates A=B^{-1}*C
  mat44 B1;
  int   i,j,k;
    B1[0][0] = 1.0; // in order to supress warnings from some
                    // stupid compilers
    Mat4Inverse ( B,B1 );
    for (i=0;i<4;i++)
      for (j=0;j<4;j++)  {
        A[i][j] = 0.0;
        for (k=0;k<4;k++)
          A[i][j] += B1[i][k]*C[k][j];
      }
  }

  void  Mat4Div2 ( mat44 & A, const mat44 & B, const mat44 & C )  {
  //  Calculates A=B*C^{-1}
  mat44 C1;
  int   i,j,k;
    C1[0][0] = 1.0; // in order to supress warnings from some
                    // stupid compilers
    Mat4Inverse ( C,C1 );
    for (i=0;i<4;i++)
      for (j=0;j<4;j++)  {
        A[i][j] = 0.0;
        for (k=0;k<4;k++)
          A[i][j] += B[i][k]*C1[k][j];
      }
  }

  void  Mat4Init ( mat44 & A )  {
  int i,j;
    for (i=0;i<4;i++)  {
      for (j=0;j<4;j++)
        A[i][j] = 0.0;
      A[i][i] = 1.0;
    }
  }

  realtype Mat4RotDet ( const mat44 & T )  {
  //  returns determinant of the rotation part
    return T[0][0]*T[1][1]*T[2][2] +
           T[0][1]*T[1][2]*T[2][0] +
           T[1][0]*T[2][1]*T[0][2] -
           T[0][2]*T[1][1]*T[2][0] -
           T[0][0]*T[1][2]*T[2][1] -
           T[2][2]*T[0][1]*T[1][0];
  }

  bool isMat4Unit ( const mat44 & A, realtype eps, bool rotOnly )  {
  // returns true if A is a unit 4x4 matrix
  int     i,j,k;
  bool B;

    if (rotOnly)  k = 3;
            else  k = 4;

    B = true;
    for (i=0;(i<k) && B;i++)
      for (j=0;(j<k) && B;j++)
        if (i==j)  B = (fabs(1.0-A[i][j])<eps);
             else  B = (fabs(A[i][j])<eps);

    return B;

  }

  void  Mat3Init ( mat33 & A )  {
  int i,j;
    for (i=0;i<3;i++)  {
      for (j=0;j<3;j++)
        A[i][j] = 0.0;
      A[i][i] = 1.0;
    }
  }

  void  Mat4Copy ( const mat44 & A, mat44 & ACopy )  {
  int i,j;
    for (i=0;i<4;i++)
      for (j=0;j<4;j++)
        ACopy[i][j] = A[i][j];
  }

  void  Mat3Copy ( const mat33 & A, mat33 & ACopy )  {
  int i,j;
    for (i=0;i<3;i++)
      for (j=0;j<3;j++)
        ACopy[i][j] = A[i][j];
  }

  bool isMat4Eq ( const mat44 & A, const mat44 & B, realtype eps,
                  bool rotOnly )  {
  // returns true if A is equal to B within precision eps
  int  i,j,k;
  bool Eq;

    if (rotOnly)  k = 3;
            else  k = 4;

    Eq = true;
    for (i=0;(i<k) && Eq;i++)
      for (j=0;(j<k) && Eq;j++)
        Eq = (fabs(A[i][j]-B[i][j])<eps);

    return Eq;

  }


  void TransformXYZ ( const mat44 & T, realtype & X, realtype & Y,
                                                     realtype & Z )  {
  realtype x1,y1,z1;
    x1 = T[0][0]*X + T[0][1]*Y + T[0][2]*Z + T[0][3];
    y1 = T[1][0]*X + T[1][1]*Y + T[1][2]*Z + T[1][3];
    z1 = T[2][0]*X + T[2][1]*Y + T[2][2]*Z + T[2][3];
    X = x1;
    Y = y1;
    Z = z1;
  }

  realtype TransformX ( const mat44 & T, realtype X, realtype Y,
                                                     realtype Z )  {
    return  T[0][0]*X + T[0][1]*Y + T[0][2]*Z + T[0][3];
  }

  realtype TransformY ( const mat44 & T, realtype X, realtype Y,
                                                     realtype Z )  {
    return  T[1][0]*X + T[1][1]*Y + T[1][2]*Z + T[1][3];
  }

  realtype TransformZ ( const mat44 & T, realtype X, realtype Y,
                                                     realtype Z )  {
    return  T[2][0]*X + T[2][1]*Y + T[2][2]*Z + T[2][3];
  }




  char CIFErrorLocation[200] = "no error";

  static cpstr msWrongSection =
    "Wrong section. The sections in PDB file may be put in wrong order.";
  static cpstr msWrongChainID =
    "Wrong chain ID. The input may have changed to another chain.";
  static cpstr msWrongEntryID =
    "Entry ID does not match the header.";

  static cpstr msSEQRES_serNum   =
    "Serial numbers of SEQRES records do not increment by 1.";
  static cpstr msSEQRES_numRes   =
    "Different SEQRES records show different numbers of residues.";
  static cpstr msSEQRES_extraRes =
    "SEQRES records contain more residues than specified.";

  static cpstr msNCSM_Unrecognized =
    "Unrecognized numerical input in MTRIXn.";
  static cpstr msNCSM_AlreadySet   =
    "Duplicate MTRIXn record.";
  static cpstr msNCSM_WrongSerial  =
    "Serial number in MTRIXn record is wrong.";
  static cpstr msNCSM_UnmatchIG    =
    "Different MTRIXn record show different iGiven flag.";

  static cpstr msATOM_Unrecognized =
    "Numerical information in ATOM record is not recognized.";
  static cpstr msATOM_AlreadySet   =
    "Atom is already in the system.";
  static cpstr msATOM_NoResidue    =
    "No residue is found for atom.";
  static cpstr msATOM_Unmatch      =
    "Unmatch in different records for the same atom.";

  static cpstr msCantOpenFile        = "File can not be opened.";
  static cpstr msUnrecognizedInteger =
    "Wrong ASCII format of an integer.";
  static cpstr msWrongModelNo        = "Wrong model number.";
  static cpstr msDuplicatedModel     = "Duplicate model number.";
  static cpstr msNoModel             = "No model defined.";
  static cpstr msForeignFile         =
    "Attempt to read unknown-type file.";
  static cpstr msWrongEdition        =
    "Attempt to read a higher-version file.";

  static cpstr msNoData              = "Expected data field not found.";
  static cpstr msUnrecognizedReal    = "Wrong ASCII format of a real.";
  static cpstr msNotACIFFile         =
    "Not a CIF file ('data_' missing).";
  static cpstr msUnrecognCIFItems    =
    "Unrecognized item(s) in CIF file.";
  static cpstr msMissingCIFField     = "Expected CIF item(s) missing.";
  static cpstr msEmptyCIFLoop        = "Empty CIF loop encountered.";
  static cpstr msUnexpEndOfCIF       = "Unexpected end of CIF file.";
  static cpstr msMissgCIFLoopField   = "Inconsistent CIF loop.";
  static cpstr msNotACIFStructure    =
    "Wrong use of CIF structure (as a loop?).";
  static cpstr msNotACIFLoop         =
    "Wrong use of CIF loop (as a structure?).";

  static cpstr msNoSheetID           = "No Sheet ID on PDB ASCII card.";
  static cpstr msWrongSheetID        = "Wrong Sheet ID.";
  static cpstr msWrongStrandNo       =
    "Wrong Strand number on PDB SHEET card.";

  static cpstr msWrongNumberOfStrands =
    "Wrong number of strands in CIF file.";
  static cpstr msWrongSheetOrder     = "Incomplete _struct_sheet_order.";
  static cpstr msHBondInconsistency  =
    "Inconsistency in _struct_sheet_hbond.";

  static cpstr msEmptyResidueName    =
          "No residue name on PDB ATOM or TER card.";
  static cpstr msDuplicateSeqNum     =
          "Duplicate sequence number and insertion code.";

  static cpstr msEmptyFile           = "Non-existent or empty file.";

  static cpstr msNoLogicalName       = "Logical file name not found.";


  cpstr  GetErrorDescription ( ERROR_CODE ErrorCode )  {

    switch (ErrorCode)  {

      case Error_NoError              :  return "No errors.";

      case Error_WrongSection         :  return msWrongSection;
      case Error_WrongChainID         :  return msWrongChainID;
      case Error_WrongEntryID         :  return msWrongEntryID;

      case Error_SEQRES_serNum        :  return msSEQRES_serNum;
      case Error_SEQRES_numRes        :  return msSEQRES_numRes;
      case Error_SEQRES_extraRes      :  return msSEQRES_extraRes;

      case Error_NCSM_Unrecognized    :  return msNCSM_Unrecognized;
      case Error_NCSM_AlreadySet      :  return msNCSM_AlreadySet;
      case Error_NCSM_WrongSerial     :  return msNCSM_WrongSerial;
      case Error_NCSM_UnmatchIG       :  return msNCSM_UnmatchIG;

      case Error_ATOM_Unrecognized    :  return msATOM_Unrecognized;
      case Error_ATOM_AlreadySet      :  return msATOM_AlreadySet;
      case Error_ATOM_NoResidue       :  return msATOM_NoResidue;
      case Error_ATOM_Unmatch         :  return msATOM_Unmatch;

      case Error_CantOpenFile         :  return msCantOpenFile;
      case Error_UnrecognizedInteger  :  return msUnrecognizedInteger;
      case Error_WrongModelNo         :  return msWrongModelNo;
      case Error_DuplicatedModel      :  return msDuplicatedModel;
      case Error_NoModel              :  return msNoModel;
      case Error_ForeignFile          :  return msForeignFile;
      case Error_WrongEdition         :  return msWrongEdition;

      case Error_NoData               :  return msNoData;
      case Error_UnrecognizedReal     :  return msUnrecognizedReal;
      case Error_NotACIFFile          :  return msNotACIFFile;
      case Error_UnrecognCIFItems     :  return msUnrecognCIFItems;
      case Error_MissingCIFField      :  return msMissingCIFField;
      case Error_EmptyCIFLoop         :  return msEmptyCIFLoop;
      case Error_UnexpEndOfCIF        :  return msUnexpEndOfCIF;
      case Error_MissgCIFLoopField    :  return msMissgCIFLoopField;
      case Error_NotACIFStructure     :  return msNotACIFStructure;
      case Error_NotACIFLoop          :  return msNotACIFLoop;

      case Error_NoSheetID            :  return msNoSheetID;
      case Error_WrongSheetID         :  return msWrongSheetID;
      case Error_WrongStrandNo        :  return msWrongStrandNo;

      case Error_WrongNumberOfStrands :  return msWrongNumberOfStrands;
      case Error_WrongSheetOrder      :  return msWrongSheetOrder;
      case Error_HBondInconsistency   :  return msHBondInconsistency;

      case Error_EmptyResidueName     :  return msEmptyResidueName;
      case Error_DuplicateSeqNum      :  return msDuplicateSeqNum;

      case Error_EmptyFile            :  return msEmptyFile;

      case Error_NoLogicalName        :  return msNoLogicalName;

      default                         :  return "Unknown error.";

    }
  }


  //  ==============  ContainerClass  ====================

  ContainerClass::ContainerClass() : io::Stream()  {
    ContinuationNo = 0;
  }

  ContainerClass::ContainerClass ( io::RPStream Object )
                 : io::Stream(Object)  {
    ContinuationNo = 0;
  }

  bool  ContainerClass::Append ( PContainerClass CC )  {
    return  (CC->ContinuationNo>1);
  }


  //  ===================  ContString  =====================

  ContString::ContString() : ContainerClass()  {
    InitString();
  }

  ContString::ContString ( cpstr S ) : ContainerClass()  {
    InitString();
    ConvertPDBASCII ( S );
  }

  ContString::ContString ( io::RPStream Object )
             : ContainerClass(Object)  {
    InitString();
  }

  ContString::~ContString() {
    if (Line)        delete[] Line;
    if (CIFCategory) delete[] CIFCategory;
    if (CIFTag)      delete[] CIFTag;
  }

  void  ContString::InitString()  {
    Line        = NULL;
    CIFCategory = NULL;
    CIFTag      = NULL;
  }

  ERROR_CODE ContString::ConvertPDBASCII ( cpstr S )  {
    CreateCopy ( Line,S );
    return Error_NoError;
  }

  void  ContString::PDBASCIIDump ( pstr S, int )  {
    if (Line)  strcpy ( S,Line );
         else  strcpy ( S,""   );
  }

  bool ContString::PDBASCIIDump1 ( io::RFile f )  {
    if (Line)  f.WriteLine ( Line );
         else  f.LF();
    return true;
  }

/*
  void ContString::GetCIF1 ( mmcif::PData CIF, ERROR_CODE & Signal,
                             int & pos )  {
  pstr F;
  int  i,RC;
  char c;
    if ((!CIFCategory) || (!CIFTag))  {
      Signal = Error_EmptyCIF;
      return;
    }
    F = CIF->GetString ( CIFCategory,CIFTag,RC );
    if (RC || (!F))  {
      Signal = Error_EmptyCIF;
      return;
    }
    if (Signal>=(int)strlen(F))  {
      CIF->DeleteField ( CIFCategory,CIFTag );
      Signal = Error_EmptyCIF;
      return;
    }
//    i = Signal;
//    while (F[i] && (F[i]!='\n') && (F[i]!='\r'))  i++;
//    if ((Signal==0) && (i==0))  {
//      i++;
//      if (((F[Signal]=='\n') && (F[i]=='\r')) ||
//          ((F[Signal]=='\r') && (F[i]=='\n')))  i++;
//      Signal = i;
//      while (F[i] && (F[i]!='\n') && (F[i]!='\r'))  i++;
//    }
//    c = F[i];
//    F[i] = char(0);
//    CreateCopy ( Line,&(F[Signal]) );
//    if (c)  {
//      F[i]   = c;
//      Signal = i+1;
//      if (((c=='\n') && (F[Signal]=='\r')) ||
//          ((c=='\r') && (F[Signal]=='\n')))  Signal++;
//    } else
//      CIF->DeleteField ( CIFCategory,CIFTag );
    i = pos;
    while (F[i] && (F[i]!='\n') && (F[i]!='\r'))  i++;
    if ((pos==0) && (i==0))  {
      i++;
      if (((F[pos]=='\n') && (F[i]=='\r')) ||
          ((F[pos]=='\r') && (F[i]=='\n')))  i++;
      pos = i;
      while (F[i] && (F[i]!='\n') && (F[i]!='\r'))  i++;
    }
    c = F[i];
    F[i] = char(0);
    CreateCopy ( Line,&(F[pos]) );
    if (c)  {
      F[i] = c;
      pos  = i+1;
      if (((c=='\n') && (F[pos]=='\r')) ||
          ((c=='\r') && (F[pos]=='\n')))  pos++;
    } else
      CIF->DeleteField ( CIFCategory,CIFTag );
  }
*/

  void  ContString::MakeCIF ( mmcif::PData CIF, int N )  {
  pstr S;
    if ((!CIFCategory) || (!CIFTag))  return;
    S = new char[strlen(Line)+5];
    strcpy ( S,"\n" );
    strcat ( S,Line );
    CIF->PutString ( S,CIFCategory,CIFTag,(N!=0) );
    delete[] S;
  }

  bool  ContString::Append ( PContainerClass CC )  {
    if (ContainerClass::Append(CC))  {
      if (!Line)  {
        Line = PContString(CC)->Line;
        PContString(CC)->Line = NULL;
      } else
        CreateConcat ( Line,pstr("\n"),PContString(CC)->Line );
      return true;
    }
    return false;
  }

  void  ContString::Copy ( PContainerClass CString )  {
    CreateCopy ( Line,PContString(CString)->Line );
  }

  void  ContString::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte   ( &Version    );
    f.CreateWrite ( Line        );
    f.CreateWrite ( CIFCategory );
    f.CreateWrite ( CIFTag      );
  }

  void  ContString::read ( io::RFile f )  {
  byte Version;
    f.ReadByte   ( &Version    );
    f.CreateRead ( Line        );
    f.CreateRead ( CIFCategory );
    f.CreateRead ( CIFTag      );
  }

  MakeStreamFunctions(ContString)



  //  ==============  ClassContainer  ====================

  MakeStreamFunctions(ContainerClass)

  ClassContainer::ClassContainer() : io::Stream()  {
    Init();
  }

  ClassContainer::ClassContainer ( io::RPStream Object )
                  : io::Stream(Object)  {
    Init();
  }

  void ClassContainer::Init()  {
    length    = 0;
    Container = NULL;
  }

  ClassContainer::~ClassContainer()  {
    FreeContainer();
  }

  void  ClassContainer::FreeContainer()  {
  int i;
    if (Container)  {
      for (i=0;i<length;i++)
        if (Container[i])
          delete Container[i];
      delete[] Container;
    }
    Container = NULL;
    length    = 0;
  }

  void  ClassContainer::AddData ( PContainerClass Data )  {
  int               i;
  PPContainerClass C1;
    if (!Data)  return;
    if (length>0)  {
      i = length-1;
      while (i>=0)  {
        if (!Container[i])  i--;
        else if (Container[i]->GetClassID()!=Data->GetClassID())  i--;
        else break;
      }
      if (i>=0)  {
        if (Container[i]->Append(Data))  {
          delete Data;
          return;
        }
      }
    }
    C1 = new PContainerClass[length+1];
    for (i=0;i<length;i++)
      C1[i] = Container[i];
    C1[length] = Data;
    if (Container)  delete[] Container;
    Container = C1;
    length++;
  }

  void  ClassContainer::PDBASCIIDump ( io::RFile f )  {
  char S[500];
  int  i,j;
    for (i=0;i<length;i++)
      if (Container[i])  {
        if (!Container[i]->PDBASCIIDump1(f))  {
          Container[i]->PDBASCIIDump ( S,i );
          j = strlen(S);
          while (j<80)  S[j++] = ' ';
          S[80] = char(0);
          f.WriteLine ( S );
        }
      }
  }

  ERROR_CODE ClassContainer::GetCIF ( mmcif::PData CIF, int ClassID )  {
  PContainerClass ContainerClass;
  int             n;
  ERROR_CODE      rc;
    n = -1;
    do  {
      ContainerClass = MakeContainerClass ( ClassID );
      rc = ContainerClass->GetCIF ( CIF,n );
      if (rc==Error_NoError)
        AddData ( ContainerClass );
    } while (rc==Error_NoError);
    delete ContainerClass;
    if (rc==Error_EmptyCIF)
      rc = Error_NoError;
    return rc;
  }

  void  ClassContainer::MakeCIF ( mmcif::PData CIF )  {
  int i;
    for (i=0;i<length;i++)
      if (Container[i])
        Container[i]->MakeCIF ( CIF,i );
  }

  void  ClassContainer::write ( io::RFile f )  {
  int  i,ClassID;
  byte Version=1;
    f.WriteByte ( &Version );
    f.WriteInt  ( &length  );
    for (i=0;i<length;i++)
      if (Container[i])  {
        ClassID = Container[i]->GetClassID();
        f.WriteInt ( &ClassID );
        Container[i]->write ( f );
      } else  {
        ClassID = -1;
        f.WriteInt ( &ClassID );
      }
  }

  PContainerClass ClassContainer::MakeContainerClass ( int ClassID )  {
    if (ClassID==ClassID_String)  return new ContString();
    return new ContainerClass();
  }

  PContainerClass ClassContainer::GetContainerClass (int ContClassNo) {
    if ((ContClassNo<0) || (ContClassNo>=length))  return NULL;
    return Container[ContClassNo];
  }

  void  ClassContainer::Copy ( PClassContainer CContainer )  {
  int i;
    FreeContainer();
    if (CContainer)  {
      length = CContainer->length;
      if (length>0)  {
        Container = new PContainerClass[length];
        for (i=0;i<length;i++)
          if (CContainer->Container[i])  {
            Container[i] = MakeContainerClass (
                              CContainer->Container[i]->GetClassID() );
            Container[i]->Copy ( CContainer->Container[i] );
          } else
            Container[i] = NULL;
      }
    }
  }

  void  ClassContainer::read ( io::RFile f )  {
  int  i,ClassID;
  byte Version;
    FreeContainer();
    f.ReadByte ( &Version );
    f.ReadInt  ( &length  );
    if (length>0)  {
      Container = new PContainerClass[length];
      for (i=0;i<length;i++)  {
        f.ReadInt ( &ClassID );
        if (ClassID>=0)  {
          Container[i] = MakeContainerClass ( ClassID );
          Container[i]->read ( f );
        } else
          Container[i] = NULL;
      }
    }
  }

  MakeStreamFunctions(ClassContainer)



  //  ======================  ID parsers  ==========================


  AtomPath::AtomPath() : io::Stream()  {
    InitAtomPath();
  }

  AtomPath::AtomPath ( cpstr ID ) : io::Stream()  {
    InitAtomPath();
    SetPath ( ID );
  }

  AtomPath::AtomPath  ( io::RPStream Object ) : io::Stream(Object) {
    InitAtomPath();
  }

  AtomPath::~AtomPath() {}

  void AtomPath::InitAtomPath()  {
    modelNo     = 0;
    chainID [0] = char(0);
    seqNum      = MinInt4;
    insCode [0] = char(0);
    resName [0] = char(0);
    atomName[0] = char(0);
    element [0] = char(0);
    altLoc  [0] = char(0);
    isSet       = 0;
  }

  int AtomPath::SetPath ( cpstr ID )  {
  //  1. If ID starts with '/':
  //   /mdl/chn/seq(res).i/atm[elm]:a
  //
  //  2. If ID starts with a letter:
  //        chn/seq(res).i/atm[elm]:a
  //
  //  3. If ID starts with a number:
  //            seq(res).i/atm[elm]:a
  //
  //  4. If ID contains colon ':' then
  //     it may be just
  //                       atm[elm]:a
  //
  //  All spaces are ignored. isSet
  // sets bit for each element present.
  // Any element may be a wildcard '*'.
  // Wildcard for model will set modelNo=0,
  // for sequence number will set
  // seqNum=MinInt4.
  //
  // Returns:
  //   0   <-> Ok
  //   -1  <-> wrong numerical format for model
  //   -2  <-> wrong numerical format for sequence number
  //
  char N[100];
  pstr p,p1;
  int  i,k;

    isSet = 0;  // clear all bits.

    p = pstr(ID);
    while (*p==' ')  p++;

    if (!(*p))  return 0;

    if (*p=='/')  {
      //  model number
      p++;
      i = 0;
      while ((*p) && (*p!='/'))  {
        if (*p!=' ')  N[i++] = *p;
        p++;
      }
      N[i] = char(0);
      if ((!N[0]) || (N[0]=='*'))  modelNo = 0;
      else {
        modelNo = mround(strtod(N,&p1));
        if ((modelNo==0) && (p1==N))  return -1;
      }
      isSet |= APATH_ModelNo;
      if (*p!='/')  return 0;
      p++;
      while (*p==' ')  p++;
    }

    if ((*p<'0') || (*p>'9'))  {
      //  chain ID
      i = 0;
      k = sizeof(ChainID)-1;
      while ((*p) && (*p!='/'))  {
        if ((*p!=' ') && (i<k))  chainID[i++] = *p;
        p++;
      }
      chainID[i] = char(0);
      if (!chainID[0])  {
        chainID[0] = '*';
        chainID[1] = char(0);
      }
      isSet |= APATH_ChainID;
      if (*p!='/')  return 0;
      p++;
      while (*p==' ')  p++;
    }

    if (((*p>='0') && (*p<='9')) || (*p=='-') ||
         (*p=='(') || (*p=='.'))  {
      //  sequence number, residue name and insertion code
      i = 0;
      while ((*p) && (*p!='/'))  {
        if (*p!=' ')  N[i++] = *p;
        p++;
      }
      N[i] = char(0);
      i    = ParseResID ( N,seqNum,insCode,resName );
      if (i==2)  return -2;
      isSet |= APATH_SeqNum | APATH_InsCode | APATH_ResName;
      if (*p!='/')  return 0;
      p++;
      while (*p==' ')  p++;
    }

    if (FirstOccurence(p,':') || FirstOccurence(p,'['))  {
      // atom name, chemical element and alternative location
      i = 0;
      while (*p)  {
        if (*p!=' ')  N[i++] = *p;
        p++;
      }
      N[i] = char(0);
      ParseAtomID ( N,atomName,element,altLoc );
      isSet |= APATH_AtomName | APATH_Element | APATH_AltLoc;
    }

    return 0;

  }

  void AtomPath::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte ( &Version );
    io::Stream::write ( f );
    f.WriteInt  ( &modelNo );
    f.WriteInt  ( &seqNum  );
    f.WriteInt  ( &isSet   );
    f.WriteTerLine ( chainID ,false );
    f.WriteTerLine ( insCode ,false );
    f.WriteTerLine ( resName ,false );
    f.WriteTerLine ( atomName,false );
    f.WriteTerLine ( element ,false );
    f.WriteTerLine ( altLoc  ,false );
  }

  void AtomPath::read ( io::RFile f )  {
  byte Version;
    f.ReadByte ( &Version );
    io::Stream::read ( f );
    f.ReadInt  ( &modelNo );
    f.ReadInt  ( &seqNum  );
    f.ReadInt  ( &isSet   );
    f.ReadTerLine ( chainID ,false );
    f.ReadTerLine ( insCode ,false );
    f.ReadTerLine ( resName ,false );
    f.ReadTerLine ( atomName,false );
    f.ReadTerLine ( element ,false );
    f.ReadTerLine ( altLoc  ,false );
  }


  MakeStreamFunctions(AtomPath)



  //  --------------------------------------------------------

  QuickSort::QuickSort() : io::Stream()  {
    selSortLimit = 15;
    data         = NULL;
    dlen         = 0;
  }

  QuickSort::QuickSort ( io::RPStream Object )
            : io::Stream(Object)  {
    selSortLimit = 15;
    data         = NULL;
    dlen         = 0;
  }

  int QuickSort::Compare ( int i, int j )  {
  // sort by increasing data[i]
    if (((ivector)data)[i]<((ivector)data)[j])  return -1;
    if (((ivector)data)[i]>((ivector)data)[j])  return  1;
    return 0;
  }

  void QuickSort::Swap ( int i, int j )  {
  int b;
    b = ((ivector)data)[i];
    ((ivector)data)[i] = ((ivector)data)[j];
    ((ivector)data)[j] = b;
  }

  void QuickSort::SelectionSort ( int left, int right )  {
  int i,j,imin;
    for (i=left;i<right;i++) {
      imin = i;
      for (j=i+1;j<=right;j++)
        if (Compare(j,imin)<0)  imin = j;
      Swap ( i,imin );
    }
  }

  int QuickSort::Partition ( int left, int right )  {
  int lv = left;
  int lm = left-1;
  int rm = right+1;
    do  {
      do
        rm--;
      while ((rm>0) && (Compare(rm,lv)>0));
      do
        lm++;
      while ((lm<dlen) && (Compare(lm,lv)<0));
      if (lm<rm)  {
        if (lv==lm)  lv = rm;
        else if (lv==rm)  lv = lm;
        Swap ( lm,rm );
      }
    } while (lm<rm);
    return rm;
  }

  void QuickSort::Quicksort ( int left, int right )  {
  int split_pt;
    if (left<(right-selSortLimit)) {
      split_pt = Partition ( left,right );
      Quicksort ( left,split_pt    );
      Quicksort ( split_pt+1,right );
    } else
      SelectionSort ( left,right );
  }

  void QuickSort::Sort ( void * sortdata, int data_len )  {
    data = sortdata;
    dlen = data_len-1;
    if (data)  Quicksort ( 0,data_len-1 );
  }

  //  --------------------------------------------------------

  void  takeWord ( pstr & p, pstr wrd, cpstr ter, int l )  {
  pstr p1;
  int  i;
    p1 = strpbrk ( p,ter );
    if (!p1)
      p1 = p + strlen(p);
    i = 0;
    while ((p!=p1) && (i<l))  {
      wrd[i++] = *p;
      p++;
    }
    if (i>=l)  i = l-1;
    wrd[i] = char(0);
    p      = p1;
  }


  void ParseAtomID ( cpstr ID, AtomName aname, Element elname,
                               AltLoc   aloc )  {
  pstr p;

    p = pstr(ID);
    while (*p==' ')  p++;

    strcpy ( aname ,"*" );
    strcpy ( elname,"*" );
    if (*p)  aloc[0] = char(0);
       else  strcpy ( aloc,"*" );

    takeWord ( p,aname,pstr("[: "),sizeof(AtomName) );

    if (*p=='[')  {
      p++;
      takeWord ( p,elname,pstr("]: "),sizeof(Element) );
      if (*p==']')  p++;
    }

    if (*p==':')  {
      p++;
      takeWord ( p,aloc,pstr(" "),sizeof(AltLoc) );
    }

  }

  int ParseResID ( cpstr ID, int & sn, InsCode inscode,
                                       ResName resname )  {
  int  RC;
  pstr p,p1;
  char N[100];

    RC = 0;

    p = pstr(ID);
    while (*p==' ')  p++;

    sn = ANY_RES;
    strcpy ( inscode,"*" );
    strcpy ( resname,"*" );

    N[0] = char(0);
    takeWord ( p,N,pstr("(./ "),sizeof(N) );
    if ((!N[0]) || (N[0]=='*'))  {
      sn = ANY_RES;
      RC = 1;
    }
    if (!RC)  {
      sn = mround(strtod(N,&p1));
      if (p1==N)  RC = 2;
            else  inscode[0] = char(0);
    }

    if (*p=='(')  {
      p++;
      takeWord ( p,resname,pstr(")./ "),sizeof(ResName) );
      if (*p==')')  p++;
    }

    if (*p=='.')  {
      p++;
      takeWord ( p,inscode,pstr("/ "),sizeof(InsCode) );
    }

    return RC;

  }

  int ParseAtomPath ( cpstr      ID,
                      int &      mdl,
                      ChainID    chn,
                      int &      sn,
                      InsCode    ic,
                      ResName    res,
                      AtomName   atm,
                      Element    elm,
                      AltLoc     aloc,
                      PAtomPath DefPath )  {
  //   /mdl/chn/seq(res).i/atm[elm]:a, may be partial
  char    N[100];
  pstr    p,p1;
  int     i,RC;
  bool wasRes;

    wasRes = false;

    RC = 0;

    p = pstr(ID);
    while (*p==' ')  p++;

    mdl = 0;
    if (*p=='/')  {
      p++;
      N[0] = char(0);
      takeWord ( p,N,pstr("/"),sizeof(N) );
      if ((!N[0]) || (N[0]=='*'))  mdl = 0;
      else {
        mdl = mround(strtod(N,&p1));
        if ((mdl==0) && (p1==N))  return -1;
      }
    } else if (DefPath)  {
      if (DefPath->isSet & APATH_ModelNo)
        mdl = DefPath->modelNo;
    }

    strcpy ( chn,"*" );
    if ((*p<'0') || (*p>'9') || (*p=='/'))  {
      if (*p=='/')  p++;
      p1 = p;
      chn[0] = char(0);
      takeWord ( p,chn,pstr("/"),sizeof(ChainID) );
      if (strpbrk(chn,"(.[:-"))  { // this was not a chain ID!
        if (DefPath)  {
          if (DefPath->isSet & APATH_ChainID)
            strcpy ( chn,DefPath->chainID );
        } else
          strcpy ( chn,"*" );
        p = p1;
      }
    } else if (DefPath)  {
      if (DefPath->isSet & APATH_ChainID)
        strcpy ( chn,DefPath->chainID );
    }

    if (*p=='/')  p++;
    sn = ANY_RES;
    strcpy ( ic ,"*" );
    strcpy ( res,"*" );
    if (((*p>='0') && (*p<='9')) || (*p=='-') ||
         (*p=='(') || (*p=='.'))  {
      wasRes = true;
      N[0] = char(0);
      takeWord ( p,N,pstr("/"),sizeof(N) );
      i = ParseResID ( N,sn,ic,res );
      if (i==2)  return -2;
    } else if (DefPath)  {
      wasRes = (*p=='/');
      if (DefPath->isSet & APATH_SeqNum)
        sn = DefPath->seqNum;
      if (DefPath->isSet & APATH_InsCode)
        strcpy ( ic,DefPath->insCode );
      if (DefPath->isSet & APATH_ResName)
        strcpy ( res,DefPath->resName );
    }

    if (*p=='/')  p++;
    strcpy ( atm ,"*" );
    strcpy ( elm ,"*" );
    strcpy ( aloc,"*" );
    if (wasRes || FirstOccurence(p,':') || FirstOccurence(p,'['))  {
      ParseAtomID ( p,atm,elm,aloc );
    } else if (DefPath)  {
      if (DefPath->isSet & APATH_AtomName)
        strcpy ( atm,DefPath->atomName );
      if (DefPath->isSet & APATH_Element)
        strcpy ( elm,DefPath->element );
      if (DefPath->isSet & APATH_ResName)
        strcpy ( aloc,DefPath->altLoc );
    }

    if (mdl<=0)       RC |= APATH_WC_ModelNo;
    if (chn[0]=='*')  RC |= APATH_WC_ChainID;
    if (sn==ANY_RES)  RC |= APATH_WC_SeqNum;
    if (ic[0]=='*')   RC |= APATH_WC_InsCode;
    if (res[0]=='*')  RC |= APATH_WC_ResName;
    if (atm[0]=='*')  RC |= APATH_WC_AtomName;
    if (elm[0]=='*')  RC |= APATH_WC_Element;
    if (aloc[0]=='*') RC |= APATH_WC_AltLoc;

    if (RC & (APATH_WC_ModelNo  | APATH_WC_ChainID |
              APATH_WC_SeqNum   | APATH_WC_InsCode |
              APATH_WC_AtomName | APATH_WC_AltLoc))
      RC |= APATH_Incomplete;

    return RC;

  }


  int  ParseSelectionPath (
               cpstr   CID,
               int &   iModel,
               pstr    Chains,
               int &   sNum1,
               InsCode ic1,
               int &   sNum2,
               InsCode ic2,
               pstr    RNames,
               pstr    ANames,
               pstr    Elements,
               pstr    altLocs
                          )  {
  int     l,j;
  pstr    p,p1;
  pstr    N;
  int     seqNum [2];
  InsCode insCode[2];
  pstr    ID;
  bool wasModel,wasChain,wasRes,haveNeg;

    l  = IMax(10,strlen(CID))+1;
    ID = new char[l];
    N  = new char[l];

    p  = pstr(CID);
    p1 = ID;
    while (*p)  {
      if (*p!=' ')  {
        *p1 = *p;
        p1++;
      }
      p++;
    }
    *p1 = char(0);

    p = ID;

    iModel = 0;
    strcpy ( Chains,"*" );
    seqNum[0] = ANY_RES;
    seqNum[1] = ANY_RES;
    strcpy ( insCode[0],"*" );
    strcpy ( insCode[1],"*" );
    strcpy ( RNames    ,"*" );
    strcpy ( ANames    ,"*" );
    strcpy ( Elements  ,"*" );
    strcpy ( altLocs   ,"*" );

    wasModel = false;
    wasChain = false;
    wasRes   = false;

    if (*p=='/')  {
      //  CID starts with the slash -- take model number first
      p++;
      N[0] = char(0);
      takeWord ( p,N,pstr("/"),l );
      if ((!N[0]) || (N[0]=='*'))  iModel = 0;
      else {
        iModel = mround(strtod(N,&p1));
        if ((iModel==0) && (p1==N))  return -1;
      }
      if (*p=='/')  p++;
      wasModel = true;
    }

    if ((*p) && (wasModel || (*p<'0') || (*p>'9')))  {
      p1 = p;
      Chains[0] = char(0);
      takeWord ( p,Chains,pstr("/"),l );
      if (strpbrk(Chains,"(.[:-"))  {  // this was not a chain ID!
        strcpy ( Chains,"*" );
        p = p1;
      } else
        wasChain = true;
      if (*p=='/')  p++;
    }

    if ((*p) && (wasChain  || ((*p>='0') && (*p<='9')) || (*p=='-') ||
                 (*p=='(') || (*p=='.')  || (*p=='*')))  {
      j = 0;
      do  {
        // take the sequence number
        haveNeg  = false;
        if (*p=='-') {
          haveNeg = true;
          p++;
        }
        N[0] = char(0);
        takeWord ( p,N,pstr("(.-/"),l );
        if ((!N[0]) || (N[0]=='*'))
          seqNum[j] = ANY_RES;
        else  {
          seqNum[j] = mround(strtod(N,&p1));
          if (p1==N)   return -2;
          if (haveNeg) seqNum[j] = - seqNum[j];
        }
        // take the residue list
        if (*p=='(')  {
          p++;
          takeWord ( p,RNames,pstr(").-/"),l );
          if (*p==')')  p++;
        }
        // take the insertion code
        if (seqNum[j]!=ANY_RES)
          insCode[j][0] = char(0);
        if (*p=='.')  {
          p++;
          takeWord ( p,insCode[j],pstr("-/"),sizeof(InsCode) );
        }
        if (*p=='-')  {
          p++;
          j++;
        } else  {
          if (j==0)  {
            seqNum[1] = seqNum[0];
            strcpy ( insCode[1],insCode[0] );
          }
          j = 10;
        }
      } while (j<2);
      wasRes = true;
    } else
      wasRes = (*p=='/');

    if (*p=='/')  p++;
    if ((*p) && (wasRes || FirstOccurence(p,':') || FirstOccurence(p,'[')))  {
      if (*p)  altLocs[0] = char(0);
      takeWord ( p,ANames,pstr("[:"),l );
      if (!ANames[0])  strcpy ( ANames,"*" );
      if (*p=='[')  {
        p++;
        takeWord ( p,Elements,pstr("]:"),l );
        if (*p==']')  p++;
      }
      if (*p==':')  {
        p++;
        takeWord ( p,altLocs,pstr(" "),l );
      }
    }

    /*
    printf ( "  iModel   = %i\n"
             "  Chains   = '%s'\n"
             "  seqNum1  = %i\n"
             "  insCode1 = '%s'\n"
             "  seqNum2  = %i\n"
             "  insCode2 = '%s'\n"
             "  RNames   = '%s'\n"
             "  ANames   = '%s'\n"
             "  Elements = '%s'\n"
             "  altLocs  = '%s'\n",
             iModel,Chains,seqNum[0],insCode[0],
             seqNum[1],insCode[1],RNames,ANames,
             Elements,altLocs );
    */

    sNum1 = seqNum[0];
    sNum2 = seqNum[1];
    strcpy ( ic1,insCode[0] );
    strcpy ( ic2,insCode[1] );

    delete[] ID;
    delete[] N;

    return 0;

  }


  void MakeSelectionPath (
               pstr       CID,
               int        iModel,
               cpstr      Chains,
               int        sNum1,
               const InsCode ic1,
               int        sNum2,
               const InsCode ic2,
               cpstr      RNames,
               cpstr      ANames,
               cpstr      Elements,
               cpstr      altLocs
                          )  {
  char S[100];
  int  k;

    if (iModel>0)  {
      sprintf ( CID,"/%i",iModel );
      k = 1;
    } else  {
      CID[0] = char(0);
      k = 0;
    }

    if (Chains[0]!='*')  {
      if (k>0)  strcat ( CID,"/"    );
      strcat ( CID,Chains );
      k = 2;
    }

    if ((sNum1!=-MaxInt4) || (ic1[0]!='*'))  {
      if (k>0)  {
        if (k<2)  strcat ( CID,"/*" );
        strcat  ( CID,"/" );
      }
      if (sNum1>-MaxInt4)  sprintf ( S,"%i",sNum1 );
                     else  strcpy  ( S,"*" );
      if (ic1[0]!='*')  {
        strcat ( S,"." );
        strcat ( S,ic1 );
      }
      strcat ( CID,S );

      if ((sNum2!=-MaxInt4) || (ic2[0]!='*'))  {
        strcat ( CID,"-" );
        if (sNum1>-MaxInt4)  sprintf ( S,"%i",sNum2 );
                       else  strcpy  ( S,"*" );
        if (ic2[0]!='*')  {
          strcat ( S,"." );
          strcat ( S,ic2 );
        }
        strcat ( CID,S );
      }

      k = 3;

    }

    if (RNames[0]!='*')  {
      if      (k<1)  strcat ( CID,"("    );
      else if (k<2)  strcat ( CID,"*/*(" );
      else if (k<3)  strcat ( CID,"/*("  );
      strcat ( CID,RNames );
      strcat ( CID,")"    );
      k = 4;
    }

    if (ANames[0]!='*')  {
      if      (k<1)  strcat ( CID,"/*/*/*/" );  // full path
      else if (k<2)  strcat ( CID,"/*/*/"   );  // /mdl + /*/*/
      else if (k<3)  strcat ( CID,"/*/"     );  // /mdl/chn + /*/
      else if (k<4)  strcat ( CID,"/"       );  // /mdl/chn/res + /
      strcat ( CID,ANames );
      strcat ( CID,")"    );
      k = 5;
    }

    if (Elements[0]!='*')  {
      if      (k<1)  strcat ( CID,"["       );
      else if (k<2)  strcat ( CID,"/*/*/*[" );
      else if (k<3)  strcat ( CID,"/*/*["   );
      else if (k<4)  strcat ( CID,"/*["     );
      else if (k<5)  strcat ( CID,"["       );
      strcat ( CID,Elements );
      strcat ( CID,"]"    );
      k = 6;
    }

    if (altLocs[0]!='*')  {
      if      (k<1)  strcat ( CID,":"       );
      else if (k<2)  strcat ( CID,"/*/*/*:" );
      else if (k<3)  strcat ( CID,"/*/*:"   );
      else if (k<4)  strcat ( CID,"/*:"     );
      else if (k<6)  strcat ( CID,":"       );
      strcat ( CID,altLocs );
    }

  }

} // namespace mmdb

