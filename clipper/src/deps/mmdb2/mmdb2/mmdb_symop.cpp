//  $Id: mmdb_symop.h $
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
//  **** Module  :   MMDB_SymOp <implementation>
//       ~~~~~~~~~
//  **** Project :   MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//
//  **** Classes :   mmdb::SymOp  ( symmetry operators )
//       ~~~~~~~~~   mmdb::SymOps ( container of symmetry operators )
//
//   (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "mmdb_symop.h"

namespace mmdb  {

  //  ====================  SymOp  ========================

  SymOp::SymOp() : io::Stream()  {
    InitSymOp();
  }

  SymOp::SymOp ( io::RPStream Object ) : io::Stream(Object)  {
    InitSymOp();
  }

  SymOp::~SymOp()  {
    FreeMemory();
  }

  void SymOp::InitSymOp()  {
  int i,j;
    XYZOp = NULL;
    for (i=0;i<4;i++)  {
      for (j=0;j<4;j++)
        T[i][j] = 0.0;
      T[i][i] = 1.0;
    }
  }

  void SymOp::FreeMemory()  {
    if (XYZOp)  delete[] XYZOp;
    XYZOp = NULL;
  }

  int  SymOp::SetSymOp ( cpstr XYZOperation )  {
  int  i,j;

    CreateCopy ( XYZOp,XYZOperation );
    DelSpaces  ( XYZOp );

    for (i=0;i<4;i++)
      for (j=0;j<4;j++)
        T[i][j] = 0.0;

    i = GetOperation ( 0 );
    if (!i) i = GetOperation ( 1 );
    if (!i) i = GetOperation ( 2 );
    T[3][3] = 1.0;

    return i;

  }

  pstr SymOp::GetSymOp()  {
    if (XYZOp)  return XYZOp;
          else  return pstr("");
  }


  int  SymOp::GetOperation ( int n )  {
  char     L[100];
  pstr     p1,p2;
  int      len;
  realtype V;

    p1 = XYZOp;
    p2 = FirstOccurence ( p1,',' );
    if (!p2)  return SYMOP_WrongSyntax;
    if (n>0)  {
      p1 = p2+1;
      p2 = FirstOccurence ( p1,',' );
      if (!p2)  return SYMOP_WrongSyntax;
    }
    if (n>1)  {
      p1 = p2+1;
      p2 = NULL;
    }

    if (p2)  *p2 = char(0);
    strcpy ( L,p1 );
    if (p2)  *p2 = ',';

    DelSpaces ( L );
    if (!L[0])  return SYMOP_WrongSyntax;
    UpperCase ( L );

    len = strlen ( L );
    T[n][0] = 0.0;
    if (L[0]=='X')  {
      T[n][0] += 1.0;
      L[0] = ' ';
    }
    do  {
      p1 = strstr ( L,"+X" );
      if (p1)  {
        T[n][0] += 1.0;
        strncpy ( p1,"  ",2 );
      }
    } while (p1);
    do  {
      p1 = strstr ( L,"-X" );
      if (p1)  {
        T[n][0] -= 1.0;
        strncpy ( p1,"  ",2 );
      }
    } while (p1);

    T[n][1] = 0.0;
    if (L[0]=='Y')  {
      T[n][1] += 1.0;
      L[0] = ' ';
    }
    do  {
      p1 = strstr ( L,"+Y" );
      if (p1)  {
        T[n][1] += 1.0;
        strncpy ( p1,"  ",2 );
      }
    } while (p1);
    do  {
      p1 = strstr ( L,"-Y" );
      if (p1)  {
        T[n][1] -= 1.0;
        strncpy ( p1,"  ",2 );
      }
    } while (p1);

    T[n][2] = 0.0;
    if (L[0]=='Z')  {
      T[n][2] += 1.0;
      L[0] = ' ';
    }
    do  {
      p1 = strstr ( L,"+Z" );
      if (p1)  {
        T[n][2] += 1.0;
        strncpy ( p1,"  ",2 );
      }
    } while (p1);
    do  {
      p1 = strstr ( L,"-Z" );
      if (p1)  {
        T[n][2] -= 1.0;
        strncpy ( p1,"  ",2 );
      }
    } while (p1);

    DelSpaces ( L );
    if ((int)strlen(L)>=len)  return SYMOP_NotAnOperation;

    // translational part
    p1 = L;
    T[n][3] = strtod ( p1,&p2 );
    if (*p2=='/')  {
      p1 = p2+1;
      V  = strtod ( p1,&p2 );
      if (V==0.0)  return SYMOP_ZeroDenominator;
      T[n][3] /= V;
    }

    return SYMOP_Ok;

  }

  void  MakeSign ( pstr S, realtype V, realtype & AV )  {
  int      l;
    if (V>0.0)  {
      l = strlen ( S );
      if (l>0)  {
        if (S[l-1]!=',')  {
          strcat ( S,"+" );
        }
      }
      AV = V;
    } else if (V<0.0)  {
      strcat ( S,"-" );
      AV = -V;
    } else  {
      AV = V;
      return;
    }
  }


  #define  __eps  1.0e-5

  void  GenTranslation ( pstr S, realtype V )  {
  realtype AV,nAV;
  char     N[50];
  int      n,d;

    if (fabs(V)<=__eps)  return;
    MakeSign ( S,V,AV );

    d = 0;
    n = -1;
    while ((d<=20) && (n<0))  {
      d++;
      nAV = AV*d;
      n   = mround(nAV);
      if (fabs(nAV-n)>__eps)  n = -1;
    }

    if (d<=1)      sprintf ( N,"%i"    ,n   );
    else if (n>=0) sprintf ( N,"%i/%i" ,n,d );
              else sprintf ( N,"%-.10g",AV  );
    strcat ( S,N );

  }

  void  GenTransformation ( pstr S, realtype V, pstr Axis ) {
  realtype AV,nAV;
  char     N[50];
  int      n,d;

    if (fabs(V)<=__eps)  return;
    MakeSign ( S,V,AV );

    if (fabs(AV-1.0)>__eps)  {

      d = 0;
      n = -1;
      while ((d<=20) && (n<0))  {
        d++;
        nAV = AV*d;
        n   = mround(nAV);
        if (fabs(nAV-n)>__eps)  n = -1;
      }

      if (n>=0)  sprintf ( N,"%i/%i*",n,d );
           else  sprintf ( N,"%-.10g*",AV );
      strcat ( S,N );

    }

    strcat ( S,Axis );

  }


  /*
  void  GenTranslation ( pstr S, realtype V )  {
  realtype AV,fAV;
  int      n,d;
  char     N[50];

    if (V==0.0)  return;
    MakeSign ( S,V,AV );

    n = mround(floor(AV+0.00000001));
    fAV = AV-n;

    if (fabs(fAV-0.5)<=__eps)                 { n += 1;  d = 2; }
    else if (fabs(fAV-0.25)<=__eps)           { n += 1;  d = 4; }
    else if (fabs(fAV-0.75)<=__eps)           { n += 3;  d = 4; }
    else if (fabs(fAV-0.33333333333)<=__eps)  { n += 1;  d = 3; }
    else if (fabs(fAV-0.66666666666)<=__eps)  { n += 2;  d = 3; }
    else if (fabs(fAV-0.16666666666)<=__eps)  { n += 1;  d = 6; }
    else if (fabs(fAV-0.83333333333)<=__eps)  { n += 5;  d = 6; }
                                        else  d = 1;

    N[0] = char(0);
    if (d>1)       sprintf  ( N,"%i/%i",n,d );
    else if (n>0)  sprintf  ( N,"%i",n );
             else  ParamStr ( N,pstr(""),AV );
    strcat ( S,N );

  }

  void  GenTransformation ( pstr S, realtype V, pstr Axis ) {
  realtype AV;

    if (V==0.0)  return;
    MakeSign ( S,V,AV );

    if (fabs(AV-0.5)<=__eps)                 strcat   ( S,"1/2*" );
    else if (fabs(AV-0.25)<=__eps)           strcat   ( S,"1/4*" );
    else if (fabs(AV-0.75)<=__eps)           strcat   ( S,"3/4*" );
    else if (fabs(AV-0.33333333333)<=__eps)  strcat   ( S,"1/3*" );
    else if (fabs(AV-0.66666666666)<=__eps)  strcat   ( S,"2/3*" );
    else if (fabs(AV-0.16666666666)<=__eps)  strcat   ( S,"1/6*" );
    else if (fabs(AV-0.83333333333)<=__eps)  strcat   ( S,"5/6*" );
    else if (fabs(AV-1.0)>__eps)             ParamStr ( S,pstr(""),AV,
                                                        10,pstr("*") );

    strcat ( S,Axis );

  }

  */

  bool  SymOp::CompileOpTitle ( pstr S )  {
    return CompileOpTitle ( S,T,true );
  }

  bool  SymOp::CompileOpTitle ( pstr S, mat44 symMat,
                                    bool compare )  {
    S[0] = char(0);
    GenTransformation ( S,symMat[0][0],pstr("X") );
    GenTransformation ( S,symMat[0][1],pstr("Y") );
    GenTransformation ( S,symMat[0][2],pstr("Z") );
    GenTranslation    ( S,symMat[0][3]           );
    strcat ( S,"," );
    GenTransformation ( S,symMat[1][0],pstr("X") );
    GenTransformation ( S,symMat[1][1],pstr("Y") );
    GenTransformation ( S,symMat[1][2],pstr("Z") );
    GenTranslation    ( S,symMat[1][3]           );
    strcat ( S,"," );
    GenTransformation ( S,symMat[2][0],pstr("X") );
    GenTransformation ( S,symMat[2][1],pstr("Y") );
    GenTransformation ( S,symMat[2][2],pstr("Z") );
    GenTranslation    ( S,symMat[2][3]           );
    DelSpaces ( S );
    if ((!compare) || (!strcmp(S,XYZOp)))  return true;
    else  {
      S[0] = char(0);
      GenTranslation    ( S,symMat[0][3]           );
      GenTransformation ( S,symMat[0][0],pstr("X") );
      GenTransformation ( S,symMat[0][1],pstr("Y") );
      GenTransformation ( S,symMat[0][2],pstr("Z") );
      strcat ( S,"," );
      GenTranslation    ( S,symMat[1][3]           );
      GenTransformation ( S,symMat[1][0],pstr("X") );
      GenTransformation ( S,symMat[1][1],pstr("Y") );
      GenTransformation ( S,symMat[1][2],pstr("Z") );
      strcat ( S,"," );
      GenTranslation    ( S,symMat[2][3]           );
      GenTransformation ( S,symMat[2][0],pstr("X") );
      GenTransformation ( S,symMat[2][1],pstr("Y") );
      GenTransformation ( S,symMat[2][2],pstr("Z") );
      DelSpaces ( S );
      if (!strcmp(S,XYZOp))  return true;
    }
    return false;
  }

  void  SymOp::Transform ( realtype & x, realtype & y, realtype & z )  {
  realtype x1,y1,z1;
    x1 = T[0][0]*x + T[0][1]*y + T[0][2]*z + T[0][3];
    y1 = T[1][0]*x + T[1][1]*y + T[1][2]*z + T[1][3];
    z1 = T[2][0]*x + T[2][1]*y + T[2][2]*z + T[2][3];
    x = x1;
    y = y1;
    z = z1;
  }

  void  SymOp::GetTMatrix ( mat44 & TMatrix )  {
  // copies T to TMatrix
  int i,j;
    for (i=0;i<4;i++)
      for (j=0;j<4;j++)
        TMatrix[i][j] = T[i][j];
  }

  void  SymOp::SetTMatrix ( mat44 & TMatrix )  {
  // copies TMatrix to T
  int i,j;
    for (i=0;i<4;i++)
      for (j=0;j<4;j++)
        T[i][j] = TMatrix[i][j];
  }


  void  SymOp::Print()  {
  int i;
    printf ( "  operation '%s'\n",XYZOp );
    for (i=0;i<4;i++)
      printf ( "               %10.3g %10.3g %10.3g  %10.3g\n",
               T[i][0],T[i][1],T[i][2],T[i][3] );
  }

  void  SymOp::Copy ( PSymOp SymOp )  {
  int i,j;
    CreateCopy ( XYZOp,SymOp->XYZOp );
    for (i=0;i<4;i++)
      for (j=0;j<4;j++)
        T[i][j] = SymOp->T[i][j];
  }

  void  SymOp::write ( io::RFile f )  {
  int  i,j;
  byte Version=1;
    f.WriteByte   ( &Version );
    f.CreateWrite ( XYZOp    );
    for (i=0;i<4;i++)
      for (j=0;j<4;j++)
        f.WriteReal ( &(T[i][j]) );
  }

  void  SymOp::read ( io::RFile f )  {
  int  i,j;
  byte Version;
    f.ReadByte   ( &Version );
    f.CreateRead ( XYZOp    );
    for (i=0;i<4;i++)
      for (j=0;j<4;j++)
        f.ReadReal ( &(T[i][j]) );
  }

  MakeStreamFunctions(SymOp)



  //  ====================  SymOps  ========================

  SymOps::SymOps() : io::Stream()  {
    InitSymOps();
  }

  SymOps::SymOps ( io::RPStream Object ) : io::Stream(Object)  {
    InitSymOps();
  }

  SymOps::~SymOps()  {
    FreeMemory();
  }

  void SymOps::InitSymOps()  {
    SpGroup = NULL;
    Nops    = 0;
    symOp   = NULL;
  }

  void SymOps::FreeMemory()  {
  int i;
    if (SpGroup)  delete[] SpGroup;
    SpGroup = NULL;
    if (symOp)  {
      for (i=0;i<Nops;i++)
        if (symOp[i])  delete symOp[i];
      delete[] symOp;
      symOp = NULL;
    }
    Nops = 0;
  }

  #define  symop_file  cpstr("symop.lib")

  int  SymOps::SetGroupSymopLib ( cpstr SpaceGroup,
                                  cpstr symop_lib )  {
  char     S[500];
  char     G[100];
  pstr     p;
  io::File f;
  int      i,RC;

    FreeMemory();

    CreateCopy ( SpGroup,SpaceGroup );

    if (!symop_lib)          p = pstr(symop_file);
    else if (!symop_lib[0])  p = pstr(symop_file);
                       else  p = pstr(symop_lib);
    f.assign ( p,true );
    if (!f.reset(true))  {
      p = getenv ( "SYMOP" );
      if (p)
        strcpy ( S,p );
      else  {
        p = getenv ( "CLIBD" );
        if (p)  {
          strcpy ( S,p );
          if (S[strlen(S)-1]!='/')  strcat ( S,"/" );
          strcat ( S,"symop.lib" );
        } else
          strcpy ( S,"symop.lib" );
      }
      f.assign ( S,true );
      if (!f.reset(true))  return SYMOP_NoLibFile;
    }

    strcpy ( G," '"    );
    strcat ( G,SpGroup );
    strcat ( G,"'"     );
    S[0] = char(0);
    while (!f.FileEnd() && (!strstr(S,G)))
      f.ReadLine ( S,sizeof(S) );
    if (f.FileEnd())  {
      f.shut();
      return SYMOP_UnknownSpaceGroup;
    }

    p = S;
    while (*p==' ')  p++;
    p = FirstOccurence ( p,' ' );
    if (p)  Nops = mround(strtod(p,NULL));
    if (Nops<=0)  return SYMOP_NoSymOps;

    symOp = new PSymOp[Nops];
    RC    = SYMOP_Ok;
    for (i=0;(i<Nops) && (!RC);i++)  {
      f.ReadLine ( S,sizeof(S) );
      symOp[i] = new SymOp();
      RC = symOp[i]->SetSymOp ( S );
    }

    f.shut();

    return RC;

  }


  #define  syminfo_file  cpstr("syminfo.lib")

  int  SymOps::SetGroup ( cpstr SpaceGroup,
                           cpstr syminfo_lib )  {
  io::File f;
  pstr     p;
  psvector lines,lines1;
  char     S[500];
  char     G[100];
  char     O[100];
  mat44    T1,T2,T3;
  int      i,j,k,l,m,n,RC;
  int      nlines,npops,ncops;

    FreeMemory();

    npops = 0;
    ncops = 0;

    CreateCopy ( SpGroup,SpaceGroup );

    if (!syminfo_lib)          p = pstr(syminfo_file);
    else if (!syminfo_lib[0])  p = pstr(syminfo_file);
                         else  p = pstr(syminfo_lib);
    f.assign ( p,true );
    if (!f.reset(true))  {
      p = getenv ( "SYMINFO" );
      if (p)
        strcpy ( S,p );
      else  {
        p = getenv ( "CLIBD" );
        if (p)  {
          strcpy ( S,p );
          if (S[strlen(S)-1]!='/')  strcat ( S,"/" );
          strcat ( S,"syminfo.lib" );
        } else
          strcpy ( S,"syminfo.lib" );
      }
      f.assign ( S,true );
      if (!f.reset(true))  return SYMOP_NoLibFile;
    }


    if (strncasecmp(SpGroup,"Hall:",5))  {
      // normal space group symbol on input
      strcpy ( G," '"    );
      strcat ( G,SpGroup );
      strcat ( G,"'"     );
      S[0] = char(0);
      while (!f.FileEnd() && !(strstr(S,G) && (strstr(S,"symbol xHM") ||
              strstr(S,"symbol old"))))
        f.ReadLine ( S,sizeof(S) );
    } else  {
      // hall descriptor on input
      strcpy ( G," ' " );
      p = &(SpGroup[5]);
      while (*p==' ')  p++;
      strcat ( G,p     );
      strcat ( G,"'"   );
      S[0] = char(0);
      while (!f.FileEnd() && !(strstr(S,G) && strstr(S,"symbol Hall")))
        f.ReadLine ( S,sizeof(S) );
    }
    if (f.FileEnd())  {
      f.shut();
      return SYMOP_UnknownSpaceGroup;
    }

    // found spacegroup, move to symop lines
    while (!f.FileEnd() && (!strstr(S,"symop")))
      f.ReadLine ( S,sizeof(S) );

    nlines = 256;
    GetVectorMemory ( lines,nlines,0 );
    for (i=0;i<nlines;i++)
      lines[i] = NULL;
    n = 0;
    CreateCopy ( lines[n++],S );

    // count primitive operators
    while (!f.FileEnd() && (strstr(S,"symop"))) {
      npops++;
      f.ReadLine ( S,sizeof(S) );
      if (n>=nlines)  {
        nlines += + 256;
        GetVectorMemory ( lines1,nlines,0 );
        for (i=0;i<n;i++)
          lines1[i] = lines[i];
        for (i=n;i<nlines;i++)
          lines1[i] = NULL;
        FreeVectorMemory ( lines,0 );
        lines = lines1;
      }
      CreateCopy ( lines[n++],S );
    }

    // count centering operators
    while (!f.FileEnd() && (strstr(S,"cenop"))) {
      ncops++;
      f.ReadLine ( S,sizeof(S) );
      if (n>=nlines)  {
        nlines += + 256;
        GetVectorMemory ( lines1,nlines,0 );
        for (i=0;i<n;i++)
          lines1[i] = lines[i];
        for (i=n;i<nlines;i++)
          lines1[i] = NULL;
        FreeVectorMemory ( lines,0 );
        lines = lines1;
      }
      CreateCopy ( lines[n++],S );
    }

    Nops = npops*ncops;
    symOp = new PSymOp[Nops];
    RC    = SYMOP_Ok;

    n = 0;  // start second pass here

    // read primitive operators
    for (i=0;(i<npops) && (!RC);i++)  {
      symOp[i] = new SymOp();
      RC = symOp[i]->SetSymOp ( lines[n++]+6 );
    }

    // loop over non-trivial centering operators, and for each loop
    // over primtive operators
    for (i=1;(i<ncops) && (!RC);i++)  {
      n++;  // this also skips the identity operator
      for (j=0;(j<npops) && (!RC);j++)  {
        symOp[i*npops+j] = new SymOp();
        RC = symOp[i*npops+j]->SetSymOp ( lines[n]+6 );
        symOp[i*npops+j]->GetTMatrix(T1);
        symOp[j]->GetTMatrix(T2);
        for (k=0;k<4;k++)
          for (l=0;l<4;l++) {
            T3[k][l] = 0.0;
            for (m=0;m<4;m++)
              T3[k][l] += T1[k][m]*T2[m][l];
          }
        for (k=0;k<3;k++)                  // kdc fix
          T3[k][3] -= floor ( T3[k][3] );  // kdc fix
        symOp[i*npops+j]->CompileOpTitle ( O,T3,false );
        symOp[i*npops+j]->SetSymOp ( O );
      }
    }

    f.shut();

    for (i=0;i<nlines;i++)
      if (lines[i])  delete[] lines[i];
    FreeVectorMemory ( lines,0 );

    return RC;

  }

  /*
  int  SymOps::SetGroup ( cpstr SpaceGroup,
                           cpstr syminfo_lib )  {
  CFile  f;
  pstr   p;
  char   S[500];
  char   G[100];
  char   O[100];
  mat44  T1,T2,T3;
  int    i,j,k,l,m,RC;
  int    npops,ncops;
  long   symop_start;
  //fpos_t symop_start;

    FreeMemory();

    npops = 0;
    ncops = 0;

    CreateCopy ( SpGroup,SpaceGroup );

    if (!syminfo_lib)          p = pstr(syminfo_file);
    else if (!syminfo_lib[0])  p = pstr(syminfo_file);
                         else  p = pstr(syminfo_lib);
    f.assign ( p,true );
    if (!f.reset(true))  {
      p = getenv ( "SYMINFO" );
      if (p)
        strcpy ( S,p );
      else  {
        p = getenv ( "CLIBD" );
        if (p)  {
          strcpy ( S,p );
          if (S[strlen(S)-1]!='/')  strcat ( S,"/" );
          strcat ( S,"syminfo.lib" );
        } else
          strcpy ( S,"syminfo.lib" );
      }
      f.assign ( S,true );
      if (!f.reset(true))  return SYMOP_NoLibFile;
    }

    if (strncasecmp(SpGroup,"Hall:",5))  {
      // normal space group symbol on input
      strcpy ( G," '"    );
      strcat ( G,SpGroup );
      strcat ( G,"'"     );
      S[0] = char(0);
      while (!f.FileEnd() && !(strstr(S,G) && (strstr(S,"symbol xHM") ||
              strstr(S,"symbol old"))))
        f.ReadLine ( S,sizeof(S) );
    } else  {
      // hall descriptor on input
      strcpy ( G," ' " );
      p = &(SpGroup[5]);
      while (*p==' ')  p++;
      strcat ( G,p     );
      strcat ( G,"'"   );
      S[0] = char(0);
      while (!f.FileEnd() && !(strstr(S,G) && strstr(S,"symbol Hall")))
        f.ReadLine ( S,sizeof(S) );
    }
    if (f.FileEnd())  {
      f.shut();
      return SYMOP_UnknownSpaceGroup;
    }

    // found spacegroup, move to symop lines
    while (!f.FileEnd() && (!strstr(S,"symop"))) {
      symop_start = f.Position();
  //#    fgetpos ( f.GetHandle(),&symop_start );
      f.ReadLine ( S,sizeof(S) );
    }
    // count primitive operators
    while (!f.FileEnd() && (strstr(S,"symop"))) {
      npops++;
      f.ReadLine ( S,sizeof(S) );
    }
    // count centering operators
    while (!f.FileEnd() && (strstr(S,"cenop"))) {
      ncops++;
      f.ReadLine ( S,sizeof(S) );
    }
    Nops = npops*ncops;
    f.seek(symop_start);
  //#  fsetpos ( f.GetHandle(),&symop_start );
    SymOp = new PSymOp[Nops];
    RC    = SYMOP_Ok;

    // read primitive operators
    for (i=0;(i<npops) && (!RC);i++)  {
      f.ReadLine ( S,sizeof(S) );
      SymOp[i] = new SymOp();
      RC = SymOp[i]->SetSymOp ( S+6 );
    }

    // skip identity centering operator
    f.ReadLine ( S,sizeof(S) );
    // loop over non-trivial centering operators, and for each loop
    // over primtive operators
    for (i=1;(i<ncops) && (!RC);i++)  {
      f.ReadLine ( S,sizeof(S) );
      for (j=0;(j<npops) && (!RC);j++)  {
        SymOp[i*npops+j] = new SymOp();
        RC = SymOp[i*npops+j]->SetSymOp ( S+6 );

        SymOp[i*npops+j]->GetTMatrix(T1);
        SymOp[j]->GetTMatrix(T2);
        for (k=0;k<4;k++)
          for (l=0;l<4;l++) {
            T3[k][l] = 0.0;
            for (m=0;m<4;m++)
              T3[k][l] += T1[k][m]*T2[m][l];
          }
        for (k=0;k<3;k++)                  // kdc fix
          T3[k][3] -= floor ( T3[k][3] );  // kdc fix
        SymOp[i*npops+j]->CompileOpTitle(O,T3,false);
        SymOp[i*npops+j]->SetSymOp (O);
      }
    }

    f.shut();

    return RC;

  }
  */

  void SymOps::Reset()  {
  // removes all symmetry operations
    FreeMemory();
  }

  int  SymOps::AddSymOp ( cpstr XYZOperation )  {
  // adds a symmetry operation
  PPSymOp symOp1;
  int      i;
    symOp1 = new PSymOp[Nops+1];
    for (i=0;i<Nops;i++)
      symOp1[i] = symOp[i];
    if (symOp) delete[] symOp;
    symOp = symOp1;
    i = Nops;
    symOp[i] = new SymOp();
    Nops++;
    return symOp[i]->SetSymOp ( XYZOperation );
  }

  void SymOps::PutGroupName ( cpstr SpGroupName )  {
    CreateCopy ( SpGroup,SpGroupName );
  }


  int  SymOps::GetNofSymOps()  {
  //  GetNofSymOps()  returns Nops -- the number of symmetry operations
    return Nops;
  }

  pstr SymOps::GetSymOp ( int Nop )  {
    if ((0<=Nop) && (Nop<Nops))  return symOp[Nop]->GetSymOp();
                           else  return pstr("");
  }

  int  SymOps::Transform ( realtype & x, realtype & y, realtype & z,
                            int Nop )  {
  //  Transform(..) transforms the coordinates according to the
  // symmetry operation Nop. The return code is non-zero if
  // Nop is a wrong operation number (must range from 0 to Nops-1).
    if ((Nop<0) || (Nop>=Nops))  return 1;
    if (symOp[Nop])  {
      symOp[Nop]->Transform ( x,y,z );
      return 0;
    } else
      return 2;
  }

  int  SymOps::GetTMatrix ( mat44 & TMatrix, int Nop )  {
  //  GetTMatrix(..) returns the coordinate transformation matrix
  // for the symmetry operation Nop. The return code is non-zero if
  // Nop is a wrong operation number (must range from 0 to Nops-1).
    if ((Nop<0) || (Nop>=Nops))  return 1;
    if (symOp[Nop])  {
      symOp[Nop]->GetTMatrix ( TMatrix );
      return 0;
    } else
      return 2;
  }

  void  SymOps::Print()  {
  int  i;
  char S[200];
    printf ( "  SPACE GROUP  '%s'\n",SpGroup );
    for (i=0;i<Nops;i++)  {
      symOp[i]->Print();
      if (symOp[i]->CompileOpTitle(S))
            printf ( " CHECK STATUS: Ok\n" );
      else  printf ( " CHECK STATUS: Generated '%s'\n",S );
    }
  }


  void  SymOps::Copy ( PSymOps SymOps )  {
  int i;
    FreeMemory();
    CreateCopy ( SpGroup,SymOps->SpGroup );
    Nops = SymOps->Nops;
    if (Nops>0)  {
      symOp = new PSymOp[Nops];
      for (i=0;i<Nops;i++) {
        symOp[i] = new SymOp();
        symOp[i]->Copy ( SymOps->symOp[i] );
      }
    }
  }


  void  SymOps::write ( io::RFile f )  {
  int  i;
  byte Version=1;
    f.WriteByte   ( &Version );
    f.CreateWrite ( SpGroup  );
    f.WriteInt    ( &Nops    );
    for (i=0;i<Nops;i++)
      StreamWrite ( f,symOp[i] );
  }

  void  SymOps::read ( io::RFile f )  {
  int  i;
  byte Version;
    FreeMemory();
    f.ReadByte   ( &Version );
    f.CreateRead ( SpGroup  );
    f.ReadInt    ( &Nops    );
    if (Nops>0)  {
      symOp = new PSymOp[Nops];
      for (i=0;i<Nops;i++) {
        symOp[i] = NULL;
        StreamRead ( f,symOp[i] );
      }
    }
  }


  MakeStreamFunctions(SymOps)

}  // namespace mmdb


/*

void TestSymOps()  {
pstr     p,p1;
int      RC;
char     S[500];
SymOps  SymOps;
CFile    f;

  p = getenv ( "PDB_ROOT" );
  if (p)  {
    strcpy ( S,p );
    strcat ( S,"/lib/" );
  } else
    S[0] = char(0);
  strcat   ( S,"symop.lib" );
  f.assign ( S,true );
  if (!f.reset())  {
    printf ( " +++++ No symop.lib file found.\n" );
    return;
  }

  while (!f.FileEnd())  {
    f.ReadLine ( S,sizeof(S) );
    if (S[0] && (S[0]!=' '))  {
      p = FirstOccurence ( S,'\'' );
      if (p)  {
        p++;
        p1 = FirstOccurence ( p,'\'' );
        if (!p1)  p = NULL;
      }
      if (!p)  {
        printf ( " +++++ Strange line in symop.lib:\n"
                 "%s\n",S );
        return;
      }
      *p1 = char(0);
      RC = SymOps.SetGroup ( p );
      printf ( " =========================================================\n"
               "  RC=%i\n",RC );
      SymOps.Print();
    }
  }

  return;

}

*/
