//  $Id: mmdb_math_align.cpp $
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
//  **** Module  :     Align <implementation>
//       ~~~~~~~~~
//  **** Classes    :  mmdb::math::Alignment  (char strings alignment)
//       ~~~~~~~~~~~~  mmdb::math::Alignment1 (int  vectors alignment)
//
//  (C) E.Krissinel  2000-2013
//
//  =================================================================
//

#include <string.h>
#include <math.h>

#include "mmdb_math_align.h"


namespace mmdb  {

  namespace math  {

    //  =====================   CAligParams   ======================

    AlignParams::AlignParams() : Stream()  {
      InitAlignParams();
    }

    AlignParams::AlignParams ( io::RPStream Object ) :
                  io::Stream ( Object )  {
      InitAlignParams();
    }

    void AlignParams::InitAlignParams()  {
      gapWeight   = -1.0;
      spaceWeight = -1.0;
      equalScore  =  2.0;
      nequalScore = -1.0;
      method      = ALIGN_GLOBAL;
    }

    void AlignParams::write ( io::RFile f )  {
      f.WriteReal ( &gapWeight   );
      f.WriteReal ( &spaceWeight );
      f.WriteReal ( &equalScore  );
      f.WriteReal ( &nequalScore );
      f.WriteInt  ( &method      );
    }

    void AlignParams::read ( io::RFile f )  {
      f.ReadReal ( &gapWeight   );
      f.ReadReal ( &spaceWeight );
      f.ReadReal ( &equalScore  );
      f.ReadReal ( &nequalScore );
      f.ReadInt  ( &method      );
    }

    MakeStreamFunctions(AlignParams)


    //  =====================   Alignment   ======================

    Alignment::Alignment() : io::Stream()  {
      InitAlignment();
    }

    Alignment::Alignment ( io::RPStream Object ) :
              io::Stream ( Object )  {
      InitAlignment();
    }

    Alignment::~Alignment()  {
      FreeMemory();
    }

    void  Alignment::InitAlignment()  {
      Space     = '-';
      SLen      = 0;
      TLen      = 0;
      VT        = NULL;
      ET        = NULL;
      FT        = NULL;
      AlgnS     = NULL;
      AlgnT     = NULL;
      AlignKey  = ALIGN_GLOBAL;
      VAchieved = 0.0;
      SEq       =  2.0;
      SNEq      = -1.0;
      Wg        =  0.0;
      Ws        = -1.0;
    }

    void  Alignment::FreeMemory()  {
      FreeMatrixMemory ( VT,TLen+1,0,0 );
      FreeMatrixMemory ( ET,TLen+1,0,0 );
      FreeMatrixMemory ( FT,TLen+1,0,0 );
      if (AlgnS)  {
        delete[] AlgnS;
        AlgnS = NULL;
      }
      if (AlgnT)  {
        delete[] AlgnT;
        AlgnT = NULL;
      }
      TLen = 0;
      SLen = 0;
    }

    void  Alignment::SetAffineModel ( realtype WGap, realtype WSpace )  {
      Wg = WGap;
      Ws = WSpace;
    }

    void  Alignment::SetScores ( realtype SEqual, realtype SNEqual )  {
      SEq  = SEqual;
      SNEq = SNEqual;
    }

    void  Alignment::Align  ( cpstr S, cpstr T, ALIGN_METHOD Method )  {
    int i,j,i0,j0;

      FreeMemory();

      AlignKey = Method;

      switch (Method)  {

        default             :
        case ALIGN_GLOBAL   : // global pairwise alignment of S and T
                              BuildGATable ( S,T, false,false );
                              VAchieved = VT[TLen][SLen];
                              Backtrace ( S,T,SLen,TLen,false );
                              if ((AlgnS[0]!=Space) && (AlgnT[0]!=Space))
                                VAchieved -= Wg;
                            break;

        case ALIGN_LOCAL    : // local pairwise alignment of S and T
                              BuildLATable ( S,T );
                              VAchieved  = 0.0;
                              i0 = -1;
                              j0 = -1;
                              for (i=0;i<=TLen;i++)
                                for (j=0;j<=SLen;j++)
                                  if (VT[i][j]>VAchieved)  {
                                    VAchieved = VT[i][j];
                                    i0 = i;
                                    j0 = j;
                                  }
                              Backtrace ( S,T,j0,i0,true );
                            break;

        case ALIGN_GLOBLOC  : // global alignment with non-penalized
                              // end gaps in T
                              BuildGATable ( S,T,false,true );
                              VAchieved = -MaxReal;
                              i0 = -1;
                              j0 = -1;
                              for (i=0;i<=TLen;i++)
                                if (VT[i][SLen]>VAchieved)  {
                                  VAchieved = VT[i][SLen];
                                  i0 = i;
                                  j0 = SLen;
                                }
                              Backtrace  ( S,T,j0,i0,false );
                              AdjustEnds ( S,T,j0,i0 );
                            break;

        case ALIGN_FREEENDS : // global alignment with non-penalized
                              // end gaps in both S and T
                              BuildGATable ( S,T,true,true );
                              VAchieved = -MaxReal;
                              i0 = -1;
                              j0 = -1;
                              for (i=0;i<=TLen;i++)
                                if (VT[i][SLen]>VAchieved)  {
                                  VAchieved = VT[i][SLen];
                                  i0 = i;
                                  j0 = SLen;
                                }
                              for (j=0;j<=SLen;j++)
                                if (VT[TLen][j]>VAchieved)  {
                                  VAchieved = VT[TLen][j];
                                  i0 = TLen;
                                  j0 = j;
                                }
                              Backtrace  ( S,T,j0,i0,false );
                              AdjustEnds ( S,T,j0,i0 );

      }

    }


    void  Alignment::BuildGATable ( cpstr S, cpstr T,
                                    bool FreeSEnd, bool FreeTEnd )  {
    int      i,j;
    realtype V1;

      SLen = strlen ( S );
      TLen = strlen ( T );
      GetMatrixMemory ( VT,TLen+1,SLen+1,0,0 );
      GetMatrixMemory ( ET,TLen+1,SLen+1,0,0 );
      GetMatrixMemory ( FT,TLen+1,SLen+1,0,0 );

      //  Base conditions
      if (FreeSEnd || FreeTEnd)  VT[0][0] = RMax(0.0,Wg);
                           else  VT[0][0] = Wg;
      ET[0][0] = VT[0][0];
      FT[0][0] = VT[0][0];

      if (FreeTEnd)
        for (i=1;i<=TLen;i++)  {
          V1       = RMax ( 0.0,VT[i-1][0]+Ws );
          VT[i][0] = V1;
          ET[i][0] = V1;
        }
      else
        for (i=1;i<=TLen;i++)  {
          V1       = VT[i-1][0] + Ws;
          VT[i][0] = V1;
          ET[i][0] = V1;
        }

      if (FreeSEnd)
        for (j=1;j<=SLen;j++)  {
          V1       = RMax ( 0.0,VT[0][j-1]+Ws );
          VT[0][j] = V1;
          FT[0][j] = V1;
        }
      else
        for (j=1;j<=SLen;j++)  {
          V1       = VT[0][j-1] + Ws;
          VT[0][j] = V1;
          FT[0][j] = V1;
        }

      //  Recurrence
      for (i=1;i<=TLen;i++)
        for (j=1;j<=SLen;j++)  {
          V1       = VT[i-1][j-1] + Score(T[i-1],S[j-1]);
          ET[i][j] = RMax ( ET[i][j-1]+Ws,VT[i][j-1]+Wg+Ws );
          FT[i][j] = RMax ( FT[i-1][j]+Ws,VT[i-1][j]+Wg+Ws );
          VT[i][j] = RMax ( RMax(V1,ET[i][j]),FT[i][j] );
        }

      FreeMatrixMemory ( ET,TLen+1,0,0 );
      FreeMatrixMemory ( FT,TLen+1,0,0 );

    //  PrintVT ( S,T );

    }


    void  Alignment::BuildLATable ( cpstr S, cpstr T )  {
    int      i,j;
    realtype V1;

      SLen = strlen ( S );
      TLen = strlen ( T );
      GetMatrixMemory ( VT,TLen+1,SLen+1,0,0 );
      GetMatrixMemory ( ET,TLen+1,SLen+1,0,0 );
      GetMatrixMemory ( FT,TLen+1,SLen+1,0,0 );

      //  Base conditions
      VT[0][0] = RMax ( 0.0,Wg );
      ET[0][0] = VT[0][0];
      FT[0][0] = VT[0][0];
      for (i=1;i<=TLen;i++)  {
        V1       = RMax ( 0.0,VT[i-1][0]+Ws );
        VT[i][0] = V1;
        ET[i][0] = V1;
      }
      for (j=1;j<=SLen;j++)  {
        V1       = RMax ( 0.0,VT[0][j-1]+Ws );
        VT[0][j] = V1;
        FT[0][j] = V1;
      }

      //  Recurrence
      for (i=1;i<=TLen;i++)
        for (j=1;j<=SLen;j++)  {
          V1       = VT[i-1][j-1] + Score(T[i-1],S[j-1]);
          ET[i][j] = RMax ( ET[i][j-1]+Ws,VT[i][j-1]+Wg+Ws );
          FT[i][j] = RMax ( FT[i-1][j]+Ws,VT[i-1][j]+Wg+Ws );
          VT[i][j] = RMax ( RMax(V1,ET[i][j]),RMax(0.0,FT[i][j]) );
        }

      FreeMatrixMemory ( ET,TLen+1,0,0 );
      FreeMatrixMemory ( FT,TLen+1,0,0 );

    //  PrintVT ( S,T );

    }

    void  Alignment::PrintVT ( cpstr S, cpstr T )  {
    int i,j;
      printf ( "\n       " );
      for (j=0;j<=SLen;j++)
        printf ( " %2i",j );
      printf ( " \n           " );
      for (j=1;j<=SLen;j++)
        printf ( " %c ",S[j-1] );
      printf ( " \n\n " );
      for (i=0;i<=TLen;i++)  {
        if (i>0)  printf ( " %2i %c ",i,T[i-1] );
            else  printf ( " %2i   ",i );
        for (j=0;j<=SLen;j++)
          printf ( " %2i",mround(VT[i][j]) );
        printf ( " \n " );
      }
      printf ( " \n" );
    }


    void  Alignment::Backtrace ( cpstr S, cpstr T, int J, int I,
                                 bool  StopAtZero )  {
    int       i,j,k, i1,j1, sk,tk;
    char      C;
    realtype  V,SV,TV;
    bool      Stop;

      //  1. Allocate memory

      if (AlgnS)  delete[] AlgnS;
      if (AlgnT)  delete[] AlgnT;

      i = SLen+TLen+1;
      AlgnS = new char[i];
      AlgnT = new char[i];
      memset ( AlgnS,Space,i );
      memset ( AlgnT,Space,i );

      //  2. Initialize backtracing
      i  = I;   // backtracing
      j  = J;   //    indices
      k  = 0;   // alignment index
      SV = 0.0;  sk = -1;   // alignment indices and leading elements
      TV = 0.0;  tk = -1;   //   for vertical and horizontal sections


      //  3. Backtracing
      Stop = false;
      while ((!Stop) && (i>0) && (j>0))  {

        V = VT[i][j];

        // find next leading element
        if (VT[i][j-1]>VT[i-1][j])  {
          i1 = i;    j1 = j-1;
        } else  {
          i1 = i-1;  j1 = j;
        }
        if (VT[i-1][j-1]>=VT[i1][j1])  {
          i1 = i-1;  j1 = j-1;
        }

    //printf ( "  i=%i  j=%i \n",i,j );

        Stop = StopAtZero && (VT[i1][j1]<=0.0);  // used at local alignment

        // treat horizontal section
        if ((sk<0) || (V>SV))  {
          sk = k;
          SV = V;
        }
        if ((j1!=j) || Stop)  {  // end of horizontal section
          AlgnS[sk] = S[j-1];
          sk = -1;
        }

        // treat vertical section
        if ((tk<0) || (V>TV))  {
          tk = k;
          TV = V;
        }
        if ((i1!=i) || Stop)  {  // end of vertical section
          AlgnT[tk] = T[i-1];
          tk = -1;
        }

        i = i1;
        j = j1;
        k++;

      }

      if (!StopAtZero)  {
        //  4. Finish the last horizontal section
        sk = k;
        while (j>0)  AlgnS[k++] = S[--j];
        //  5. Finish the last vertical section
        while (i>0)  AlgnT[sk++] = T[--i];
        k = IMax ( k,sk );
      }

      //  6. Put the termination character
      AlgnS[k] = char(0);
      AlgnT[k] = char(0);

      //  7. Reverse the strings
      i = 0;
      j = k-1;
      if (StopAtZero)  {
        // should work only for local alignment
        while ((j>0) && ((AlgnS[j]==Space) || (AlgnT[j]==Space)))  j--;
        k = j+1;
        AlgnS[k] = char(0);
        AlgnT[k] = char(0);
      }
      while (j>i)  {
        C = AlgnS[i];  AlgnS[i] = AlgnS[j];  AlgnS[j] = C;
        C = AlgnT[i];  AlgnT[i] = AlgnT[j];  AlgnT[j] = C;
        i++;
        j--;
      }

      //  8. Collapse the alternating spaces
      do  {
        k = 0;
        i = 0;
        while (AlgnS[k])  {
          if ((AlgnS[k]==Space) && (AlgnT[k]==Space))  k++;
          else if ((AlgnS[k]==Space) && (AlgnS[k+1]!=Space) &&
                   (AlgnT[k]!=Space) && (AlgnT[k+1]==Space))  {
            AlgnS[i] = AlgnS[k+1];
            AlgnT[i] = AlgnT[k];
            k++;
          } else if ((AlgnS[k]!=Space) && (AlgnS[k+1]==Space) &&
                     (AlgnT[k]==Space) && (AlgnT[k+1]!=Space))  {
            AlgnS[i] = AlgnS[k];
            AlgnT[i] = AlgnT[k+1];
            k++;
          } else if (i!=k)  {
            AlgnS[i] = AlgnS[k];
            AlgnT[i] = AlgnT[k];
          }
          if (AlgnS[k])  {
            k++;
            i++;
          }
        }
        if (i!=k)  {  // terminating character
          AlgnS[i] = AlgnS[k];
          AlgnT[i] = AlgnT[k];
        }
      } while (k>i);

    }


    void  Alignment::AdjustEnds ( cpstr S, cpstr T, int J, int I )  {
    int si,ti,m;

      if (J<SLen)  strcat ( AlgnS,&(S[J]) );
      if (I<TLen)  strcat ( AlgnT,&(T[I]) );
      si = strlen ( AlgnS );
      ti = strlen ( AlgnT );
      m  = IMax ( si,ti );
      while (si<m)  AlgnS[si++] = Space;
      while (ti<m)  AlgnT[ti++] = Space;
      AlgnS[si] = char(0);
      AlgnT[ti] = char(0);

    /*
    int k,m;

      if (J>I)  {
        k = J-I;
        strcat ( AlgnT,&(T[IMax(0,TLen-k)]) );
        k = strlen ( AlgnS );
        m = strlen ( AlgnT );
        while (k<m)
          AlgnS[k++] = Space;
        AlgnS[k] = char(0);
      } else if (I>J)  {
        k = I-J;
        strcat ( AlgnS,&(S[IMax(0,SLen-k)]) );
        k = strlen ( AlgnT );
        m = strlen ( AlgnS );
        while (k<m)
          AlgnT[k++] = Space;
        AlgnT[k] = char(0);
      }
    */

    }


    realtype Alignment::Score ( char A, char B )  {
      if (A==B)  return SEq;
      if ((A==Space) || (B==Space))  return Ws;
      return SNEq;
    }

    realtype Alignment::GetSimilarity()  {
    realtype s,a;
    int      i,n;

      s = 0.0;
      a = 0.0;
      n = IMin ( strlen(AlgnS),strlen(AlgnT) );

      for (i=0;i<n;i++)
        if ((AlgnS[i]!=Space) || (AlgnT[i]!=Space))  {
          a += RMax ( Score(AlgnS[i],AlgnS[i]),Score(AlgnT[i],AlgnT[i]) );
          s += Score ( AlgnS[i],AlgnT[i] );
        }

      if ((s>0.0) && (a>0.0))  return s/a;
      return 0.0;

    }


    int Alignment::GetNAlign()  {
    int      i,n,ne;
      ne = 0;
      n  = IMin ( strlen(AlgnS),strlen(AlgnT) );
      for (i=0;i<n;i++)
        if ((AlgnT[i]!=Space) && (AlgnS[i]==AlgnT[i]))
          ne++;
      return ne;
    }

    realtype Alignment::GetSeqId()  {
    realtype s;
    int      i,n,ne,ns,nt;

      ne = 0;
      ns = 0;
      nt = 0;
      n  = IMin ( strlen(AlgnS),strlen(AlgnT) );

      for (i=0;i<n;i++)  {
        if (AlgnS[i]!=Space)  ns++;
        if (AlgnT[i]!=Space)  {
          nt++;
          if (AlgnS[i]==AlgnT[i])
            ne++;
        }
      }

      s = IMin ( ns,nt );
      if (s>0.0)  return ne/s;
      return 0.0;

    }


    #define  WrapPeriod  61

    void  Alignment::OutputResults ( io::RFile f, cpstr S, cpstr T )  {
    int   k,l,n;
    char  P[3];

      P[1] = char(0);
      if ((!AlgnS) || (!AlgnT))  {
        f.LF();
        f.WriteLine ( pstr(" NO ALIGNMENT HAS BEEN DONE.") );
        f.shut();
        return;
      }
      f.LF();
      f.WriteLine ( pstr(" ========  INPUT DATA") );
      f.LF();
      f.WriteLine ( pstr(" String S:") );
      f.Write ( pstr(" ") );
      l = 1;
      k = 0;
      while (S[k])  {
        P[0] = S[k++];
        f.Write ( P );
        l++;
        if (l>=WrapPeriod)  {
          f.LF();  f.Write ( pstr(" ") );  l = 1;
        }
      }
      f.LF();
      f.LF();
      f.WriteLine ( pstr(" String T:") );
      f.Write ( pstr(" ") );
      l = 1;
      k = 0;
      while (T[k])  {
        P[0] = T[k++];
        f.Write ( P );
        l++;
        if (l>=WrapPeriod)  {
          f.LF();  f.Write ( pstr(" ") );  l = 1;
        }
      }
      f.LF();
      f.LF();
      f.WriteParameter ( pstr(" Score equal")  ,SEq ,20,10 );
      f.WriteParameter ( pstr(" Score unequal"),SNEq,20,10 );
      f.LF();
      f.WriteParameter ( pstr(" Gap weight")   ,Wg  ,20,10 );
      f.WriteParameter ( pstr(" Space weight") ,Ws  ,20,10 );
      f.LF();
      f.LF();
      f.Write ( pstr(" ========  RESULT OF ") );
      switch (AlignKey)  {
        default             :
        case ALIGN_GLOBAL   : f.Write ( pstr("GLOBAL")       );  break;
        case ALIGN_LOCAL    : f.Write ( pstr("LOCAL")        );  break;
        case ALIGN_GLOBLOC  : f.Write ( pstr("GLOBAL/LOCAL") );  break;
        case ALIGN_FREEENDS : f.Write ( pstr("FREE-ENDS")    );
      }
      f.WriteLine ( pstr(" ALIGNMENT") );
      f.LF();
      if (AlignKey==ALIGN_GLOBLOC)  {
        f.WriteLine ( pstr(" End gaps in T-string were not penalized") );
        f.LF();
      }
      f.WriteParameter ( pstr(" Highest score achieved:"),VAchieved,26,10 );
      f.LF();
      f.WriteLine ( pstr(" Aligned S (upper string) and T (lower string):") );
      f.LF();
      k = 0;
      n = 0;
      l = 1;  f.Write ( pstr(" ") );
      while (AlgnS[k])  {
        P[0] = AlgnS[k++];
        f.Write ( P );
        l++;
        if ((l>=WrapPeriod) || (!AlgnS[k]))  {
          f.LF();  f.Write ( pstr(" ") );  l = 1;
          while (AlgnT[n] && (l<WrapPeriod))  {
            P[0] = AlgnT[n++];
            f.Write ( P );
            l++;
          }
          f.LF(); f.LF(); f.Write ( pstr(" ") );  l = 1;
        }
      }

    }


    //  -----------------  Streaming  -----------------------------

    void  Alignment::write ( io::RFile f )  {
    int Version=1;
      f.WriteFile ( &Version,sizeof(Version) );
      Stream::write ( f );
    }

    void  Alignment::read ( io::RFile f )  {
    int Version;
      f.ReadFile ( &Version,sizeof(Version) );
      Stream::write ( f );
    }



    //  =====================   Alignment1   ======================

    Alignment1::Alignment1() : io::Stream()  {
      InitAlignment1();
    }

    Alignment1::Alignment1 ( io::RPStream Object ) :
                io::Stream ( Object ) {
      InitAlignment1();
    }

    Alignment1::~Alignment1()  {
      FreeMemory();
    }

    void  Alignment1::InitAlignment1()  {
      Space     = 0;
      SLen      = 0;
      TLen      = 0;
      AlgnLen   = 0;
      VT        = NULL;
      ET        = NULL;
      FT        = NULL;
      AlgnS     = NULL;
      AlgnT     = NULL;
      AlignKey  = ALIGN_GLOBAL;
      VAchieved = 0.0;
      SEq       =  2.0;
      SNEq      = -1.0;
      Wg        =  0.0;
      Ws        = -1.0;
    }

    void  Alignment1::FreeMemory()  {
      FreeMatrixMemory ( VT,TLen+1,0,0 );
      FreeMatrixMemory ( ET,TLen+1,0,0 );
      FreeMatrixMemory ( FT,TLen+1,0,0 );
      FreeVectorMemory ( AlgnS,0 );
      FreeVectorMemory ( AlgnT,0 );
      TLen    = 0;
      SLen    = 0;
      AlgnLen = 0;
    }

    void  Alignment1::SetAffineModel ( realtype WGap, realtype WSpace )  {
      Wg = WGap;
      Ws = WSpace;
    }

    void  Alignment1::SetScores ( realtype SEqual, realtype SNEqual )  {
      SEq  = SEqual;
      SNEq = SNEqual;
    }

    void  Alignment1::Align  ( ivector S, int SLength,
                               ivector T, int TLength,
                               ALIGN_METHOD Method )  {
    int i,j,i0,j0;

      FreeMemory();

      SLen = SLength;
      TLen = TLength;

      AlignKey = Method;

      switch (Method)  {

        default             :
        case ALIGN_GLOBAL   : // global pairwise alignment of S and T
                              BuildGATable ( S,T, false,false );
                              VAchieved = VT[TLen][SLen];
                              Backtrace ( S,T,SLen,TLen,false );
                              if ((AlgnS[0]!=Space) && (AlgnT[0]!=Space))
                                VAchieved -= Wg;
                            break;

        case ALIGN_LOCAL    : // local pairwise alignment of S and T
                              BuildLATable ( S,T );
                              VAchieved = 0.0;
                              i0 = -1;
                              j0 = -1;
                              for (i=0;i<=TLen;i++)
                                for (j=0;j<=SLen;j++)
                                  if (VT[i][j]>VAchieved)  {
                                    VAchieved = VT[i][j];
                                    i0 = i;
                                    j0 = j;
                                  }
                              Backtrace ( S,T,j0,i0,true );
                            break;

        case ALIGN_GLOBLOC  : // global alignment with non-penalized
                              // end gaps in T
                              BuildGATable ( S,T,false,true );
                              VAchieved = -MaxReal;
                              i0 = -1;
                              j0 = -1;
                              for (i=0;i<=TLen;i++)
                                if (VT[i][SLen]>VAchieved)  {
                                  VAchieved = VT[i][SLen];
                                  i0 = i;
                                  j0 = SLen;
                                }
                              Backtrace  ( S,T,j0,i0,false );
                              AdjustEnds ( S,T,j0,i0 );
                            break;

        case ALIGN_FREEENDS : // global alignment with non-penalized
                              // end gaps in both S and T
                              BuildGATable ( S,T,true,true );
                              VAchieved = -MaxReal;
                              i0 = -1;
                              j0 = -1;
                              for (i=0;i<=TLen;i++)
                                if (VT[i][SLen]>VAchieved)  {
                                  VAchieved = VT[i][SLen];
                                  i0 = i;
                                  j0 = SLen;
                                }
                              for (j=0;j<=SLen;j++)
                                if (VT[TLen][j]>VAchieved)  {
                                  VAchieved = VT[TLen][j];
                                  i0 = TLen;
                                  j0 = j;
                                }
                              Backtrace  ( S,T,j0,i0,false );
                              AdjustEnds ( S,T,j0,i0 );
      }

    }


    void  Alignment1::BuildGATable ( ivector S, ivector T,
                                     bool FreeSEnd, bool FreeTEnd )  {
    int      i,j;
    realtype V1;

      GetMatrixMemory ( VT,TLen+1,SLen+1,0,0 );
      GetMatrixMemory ( ET,TLen+1,SLen+1,0,0 );
      GetMatrixMemory ( FT,TLen+1,SLen+1,0,0 );

      //  Base conditions
      if (FreeSEnd || FreeTEnd)  VT[0][0] = RMax(0.0,Wg);
                           else  VT[0][0] = Wg;
      ET[0][0] = VT[0][0];
      FT[0][0] = VT[0][0];

      if (FreeTEnd)
        for (i=1;i<=TLen;i++)  {
          V1       = RMax ( 0.0,VT[i-1][0]+Ws );
          VT[i][0] = V1;
          ET[i][0] = V1;
        }
      else
        for (i=1;i<=TLen;i++)  {
          V1       = VT[i-1][0] + Ws;
          VT[i][0] = V1;
          ET[i][0] = V1;
        }

      if (FreeSEnd)
        for (j=1;j<=SLen;j++)  {
          V1       = RMax ( 0.0,VT[0][j-1]+Ws );
          VT[0][j] = V1;
          FT[0][j] = V1;
        }
      else
        for (j=1;j<=SLen;j++)  {
          V1       = VT[0][j-1] + Ws;
          VT[0][j] = V1;
          FT[0][j] = V1;
        }

      //  Recurrence
      for (i=1;i<=TLen;i++)
        for (j=1;j<=SLen;j++)  {
          V1       = VT[i-1][j-1] + Score(T[i-1],S[j-1]);
          ET[i][j] = RMax ( ET[i][j-1]+Ws,VT[i][j-1]+Wg+Ws );
          FT[i][j] = RMax ( FT[i-1][j]+Ws,VT[i-1][j]+Wg+Ws );
          VT[i][j] = RMax ( RMax(V1,ET[i][j]),FT[i][j] );
        }

      FreeMatrixMemory ( ET,TLen+1,0,0 );
      FreeMatrixMemory ( FT,TLen+1,0,0 );

    //  PrintVT ( S,T );

    }


    void  Alignment1::BuildLATable ( ivector S, ivector T )  {
    int      i,j;
    realtype V1;

      GetMatrixMemory ( VT,TLen+1,SLen+1,0,0 );
      GetMatrixMemory ( ET,TLen+1,SLen+1,0,0 );
      GetMatrixMemory ( FT,TLen+1,SLen+1,0,0 );

      //  Base conditions
      VT[0][0] = RMax ( 0.0,Wg );
      ET[0][0] = VT[0][0];
      FT[0][0] = VT[0][0];
      for (i=1;i<=TLen;i++)  {
        V1       = RMax ( 0.0,VT[i-1][0]+Ws );
        VT[i][0] = V1;
        ET[i][0] = V1;
      }
      for (j=1;j<=SLen;j++)  {
        V1       = RMax ( 0.0,VT[0][j-1]+Ws );
        VT[0][j] = V1;
        FT[0][j] = V1;
      }

      //  Recurrence
      for (i=1;i<=TLen;i++)
        for (j=1;j<=SLen;j++)  {
          V1       = VT[i-1][j-1] + Score(T[i-1],S[j-1]);
          ET[i][j] = RMax ( ET[i][j-1]+Ws,VT[i][j-1]+Wg+Ws );
          FT[i][j] = RMax ( FT[i-1][j]+Ws,VT[i-1][j]+Wg+Ws );
          VT[i][j] = RMax ( RMax(V1,ET[i][j]),RMax(0.0,FT[i][j]) );
        }

      FreeMatrixMemory ( ET,TLen+1,0,0 );
      FreeMatrixMemory ( FT,TLen+1,0,0 );

    //  PrintVT ( S,T );

    }

    void  Alignment1::PrintVT ( ivector S, ivector T )  {
    int i,j;
      printf ( "\n       " );
      for (j=0;j<=SLen;j++)
        printf ( " %2i",j );
      printf ( " \n           " );
      for (j=1;j<=SLen;j++)
        printf ( " %3i ",S[j-1] );
      printf ( " \n\n " );
      for (i=0;i<=TLen;i++)  {
        if (i>0)  printf ( " %2i %3i ",i,T[i-1] );
            else  printf ( " %2i   ",i );
        for (j=0;j<=SLen;j++)
          printf ( " %2i",mround(VT[i][j]) );
        printf ( " \n " );
      }
      printf ( " \n" );
    }


    void  Alignment1::Backtrace ( ivector S, ivector T, int J, int I,
                                  bool StopAtZero )  {
    int       i,j,k, i1,j1, sk,tk;
    int       C;
    realtype  V,SV,TV;
    bool      Stop;

      //  1. Allocate memory

      FreeVectorMemory ( AlgnS,0 );
      FreeVectorMemory ( AlgnT,0 );
      AlgnLen = 0;

      k = SLen+TLen+1;
      GetVectorMemory  ( AlgnS,k,0 );
      GetVectorMemory  ( AlgnT,k,0 );
      for (i=0;i<k;i++)  {
        AlgnS[i] = Space;
        AlgnT[i] = Space;
      }

      //  2. Initialize backtracing
      i = I;   // backtracing
      j = J;   //    indices

      k  = 0;   // alignment index
      SV = 0.0;  sk = -1;   // alignment indices and leading elements
      TV = 0.0;  tk = -1;   //   for vertical and horizontal sections


      //  3. Backtracing
      Stop = false;
      while ((!Stop) && (i>0) && (j>0))  {

        V = VT[i][j];

        // find next leading element
        if (VT[i][j-1]>VT[i-1][j])  {
          i1 = i;    j1 = j-1;
        } else  {
          i1 = i-1;  j1 = j;
        }
        if (VT[i-1][j-1]>=VT[i1][j1]) {
          i1 = i-1;  j1 = j-1;
        }

        Stop = StopAtZero && (VT[i1][j1]<=0.0);  // used at local alignment

        // treat horizontal section
        if ((sk<0) || (V>SV))  {
          sk = k;
          SV = V;
        }
        if ((j1!=j) || Stop)  {  // end of horizontal section
          AlgnS[sk] = S[j-1];
          sk = -1;
        }

        // treat vertical section
        if ((tk<0) || (V>TV))  {
          tk = k;
          TV = V;
        }
        if ((i1!=i) || Stop)  {  // end of vertical section
          AlgnT[tk] = T[i-1];
          tk = -1;
        }

        i = i1;
        j = j1;
        k++;

      }

      if (!StopAtZero)  {
        //  4. Finish the last horizontal section
        sk = k;
        while (j>0)  AlgnS[k++] = S[--j];
        //  5. Finish the last vertical section
        while (i>0)  AlgnT[sk++] = T[--i];
        k = IMax ( k,sk );
      }

      //  6. Put the termination character
      AlgnLen  = k;

      //  7. Reverse the strings
      i = 0;
      j = k-1;
      if (StopAtZero)  {
        // should work only for local alignment
        while ((j>0) && ((AlgnS[j]==Space) || (AlgnT[j]==Space)))  j--;
        AlgnLen = j+1;
      }
      while (j>i)  {
        C = AlgnS[i];  AlgnS[i] = AlgnS[j];  AlgnS[j] = C;
        C = AlgnT[i];  AlgnT[i] = AlgnT[j];  AlgnT[j] = C;
        i++;
        j--;
      }

      //  8. Filter out parasite spaces
      k = 0;
      i = 0;
      while (k<AlgnLen) {
        while ((k<AlgnLen) && (AlgnS[k]==Space) && (AlgnT[k]==Space))  k++;
        if (k<AlgnLen) {
          AlgnS[i] = AlgnS[k];
          AlgnT[i] = AlgnT[k];
          k++;
          i++;
        }
      }

      AlgnLen = i;

      //  9. Collapse the alternating spaces
      do  {

        k = 0;
        i = 0;
        while (k<AlgnLen)  {
          if ((AlgnS[k]==Space) && (AlgnT[k]==Space))  k++;
          else if ((k+1<AlgnLen) &&
                   (AlgnS[k]==Space) && (AlgnS[k+1]!=Space) &&
                   (AlgnT[k]!=Space) && (AlgnT[k+1]==Space))  {
            AlgnS[i] = AlgnS[k+1];
            AlgnT[i] = AlgnT[k];
            k++;
          } else if ((k+1<AlgnLen) &&
                     (AlgnS[k]!=Space) && (AlgnS[k+1]==Space) &&
                     (AlgnT[k]==Space) && (AlgnT[k+1]!=Space))  {
            AlgnS[i] = AlgnS[k];
            AlgnT[i] = AlgnT[k+1];
            k++;
          } else if (i!=k)  {
            AlgnS[i] = AlgnS[k];
            AlgnT[i] = AlgnT[k];
          }
          if (k<AlgnLen)  {
            k++;
            i++;
          }
        }

        AlgnLen = i;

      } while (k>i);


    }


    void  Alignment1::AdjustEnds ( ivector S, ivector T, int J, int I )  {
    int is,it;
      is = J;
      it = I;
      while ((is<SLen) || (it<TLen))  {
        if (is<SLen)  AlgnS[AlgnLen] = S[is];
                else  AlgnS[AlgnLen] = Space;
        if (it<TLen)  AlgnT[AlgnLen] = T[it];
                else  AlgnT[AlgnLen] = Space;
        is++;
        it++;
        AlgnLen++;
      }
    }

    realtype Alignment1::Score ( int A, int B )  {
      if (A==B)  {
        if (A==Space)  return 0.0;
                 else  return SEq;
      }
      if ((A==Space) || (B==Space))  return Ws;
      return SNEq;
    }


    realtype Alignment1::GetSimilarity()  {
    realtype s,a;
    int      i;

      s = 0.0;
      a = 0.0;

      for (i=0;i<AlgnLen;i++)
        if ((AlgnS[i]!=Space) || (AlgnT[i]!=Space))  {
          a += RMax ( Score(AlgnS[i],AlgnS[i]),Score(AlgnT[i],AlgnT[i]) );
          s += Score ( AlgnS[i],AlgnT[i] );
        }

      if ((s>0.0) && (a>0.0))  return s/a;
      return 0.0;

    }


    void  Alignment1::OutputResults ( io::RFile  f, ivector S, int lenS,
                                      ivector T, int lenT )  {
    int   k,l,n;
    char  P[10];

      if ((!AlgnS) || (!AlgnT))  {
        f.LF();
        f.WriteLine ( pstr(" NO ALIGNMENT HAS BEEN DONE.") );
        f.shut();
        return;
      }
      f.LF();
      f.WriteLine ( pstr(" ========  INPUT DATA") );
      f.LF();
      f.WriteLine ( pstr(" String S:") );
      f.Write ( pstr(" ") );
      l = 1;
      k = 0;
      while (k<lenS)  {
        sprintf ( P,"%4i ",S[k++] );
        f.Write ( P );
        l += 5;
        if (l>=WrapPeriod)  {
          f.LF();  f.Write ( pstr(" ") );  l = 1;
        }
      }
      f.LF();
      f.LF();
      f.WriteLine ( pstr(" String T:") );
      f.Write ( pstr(" ") );
      l = 1;
      k = 0;
      while (k<lenT)  {
        sprintf ( P,"%4i ",T[k++] );
        f.Write ( P );
        l += 5;
        if (l>=WrapPeriod)  {
          f.LF();  f.Write ( pstr(" ") );  l = 1;
        }
      }
      f.LF();
      f.LF();
      f.WriteParameter ( pstr(" Score equal")  ,SEq ,20,10 );
      f.WriteParameter ( pstr(" Score unequal"),SNEq,20,10 );
      f.LF();
      f.WriteParameter ( pstr(" Gap weight")   ,Wg  ,20,10 );
      f.WriteParameter ( pstr(" Space weight") ,Ws  ,20,10 );
      f.LF();
      f.LF();
      f.Write ( pstr(" ========  RESULT OF ") );
      switch (AlignKey)  {
        default             :
        case ALIGN_GLOBAL   : f.Write ( pstr("GLOBAL")       );  break;
        case ALIGN_LOCAL    : f.Write ( pstr("LOCAL")        );  break;
        case ALIGN_GLOBLOC  : f.Write ( pstr("GLOBAL/LOCAL") );  break;
        case ALIGN_FREEENDS : f.Write ( pstr("FREE-ENDS")    );
      }
      f.WriteLine ( pstr(" ALIGNMENT") );
      f.LF();
      if (AlignKey==ALIGN_GLOBLOC)  {
        f.WriteLine ( pstr(" End gaps in T-string were not penalized") );
        f.LF();
      }
      f.WriteParameter ( pstr(" Highest score achieved:"),
                         VAchieved,26,10 );
      f.LF();
      f.WriteLine ( pstr(" Aligned S (upper string) and T "
                         "(lower string):") );
      f.LF();
      k = 0;
      n = 0;
      l = 1;  f.Write ( pstr(" ") );
      while (k<AlgnLen)  {
        sprintf ( P,"%4i ",AlgnS[k++] );
        f.Write ( P );
        l += 5;
        if ((l>=WrapPeriod) || (!AlgnS[k]))  {
          f.LF();  f.Write ( pstr(" ") );  l = 1;
          while ((n<AlgnLen) && (l<WrapPeriod))  {
            sprintf ( P,"%4i ",AlgnT[n++] );
            f.Write ( P );
            l += 5;
          }
          f.LF(); f.LF(); f.Write ( pstr(" ") );  l = 1;
        }
      }

    }


    //  -----------------  Streaming  -----------------------------

    void  Alignment1::write ( io::RFile f )  {
    int Version=1;
      f.WriteFile ( &Version,sizeof(Version) );
      Stream::write ( f );
    }

    void  Alignment1::read ( io::RFile f )  {
    int Version;
      f.ReadFile ( &Version,sizeof(Version) );
      Stream::write ( f );
    }


  }  // namespace math

}  // namespace mmdb
