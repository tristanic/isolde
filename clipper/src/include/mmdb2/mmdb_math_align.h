//  $Id: mmdb_math_align.h $
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
//  **** Module  :     Align <interface>
//       ~~~~~~~~~
//  **** Classes    :  mmdb::math::Alignment  (char strings alignment)
//       ~~~~~~~~~~~~  mmdb::math::Alignment1 (int  vectors alignment)
//
//  (C) E.Krissinel  2000-2013
//
//  =================================================================
//

#ifndef __MMDB_MATH_Align__
#define __MMDB_MATH_Align__

#include "mmdb_io_stream.h"

namespace mmdb  {

  namespace math  {

    //  =====================   AlignParams   ======================

    DefineClass(AlignParams);
    DefineStreamFunctions(AlignParams);

    class AlignParams : public io::Stream  {

      public :

        realtype  gapWeight,spaceWeight;
        realtype  equalScore,nequalScore;
        int       method;

        AlignParams();
        AlignParams ( io::RPStream Object );

        void write ( io::RFile f );
        void read  ( io::RFile f );

      protected :
        void InitAlignParams();

    };


    //  ======================   Alignment   =======================

    DefineClass(Alignment);

    enum ALIGN_METHOD  {
      ALIGN_GLOBAL   = 0,
      ALIGN_LOCAL    = 1,
      ALIGN_GLOBLOC  = 2,
      ALIGN_FREEENDS = 3
    };

    class  Alignment : public io::Stream  {

      public :

        Alignment  ();
        Alignment  ( io::RPStream Object );
        ~Alignment ();

        void SetAffineModel ( realtype WGap,   realtype WSpace  );
        void SetScores      ( realtype SEqual, realtype SNEqual );

        void Align          ( cpstr S, cpstr T,
                              ALIGN_METHOD Method=ALIGN_GLOBAL );

        inline pstr     GetAlignedS()  {  return AlgnS;      }
        inline pstr     GetAlignedT()  {  return AlgnT;      }
        inline realtype GetScore   ()  {  return VAchieved;  }
        inline char     GetSpace   ()  {  return Space;      }

        realtype GetSimilarity(); // Score-weighted sequence id
        realtype GetSeqId     (); // Primitive sequence id
        int      GetNAlign    (); // number of aligned residues

        virtual void OutputResults ( io::RFile f, cpstr S, cpstr T  );

        void read  ( io::RFile f );
        void write ( io::RFile f );

      protected :

        char     Space;
        int      AlignKey, SLen,TLen;
        rmatrix  VT,ET,FT;
        pstr     AlgnS,AlgnT;
        realtype VAchieved;
        realtype SEq,SNEq, Wg,Ws;

        virtual void  InitAlignment();
        virtual void  FreeMemory   ();
        virtual realtype  Score    ( char A, char B );

        void    BuildGATable ( cpstr S, cpstr T,
                               bool FreeSEnd, bool FreeTEnd );
        void    BuildLATable ( cpstr S, cpstr T );
        void    Backtrace    ( cpstr S, cpstr T, int J, int I,
                               bool StopAtZero );
        void    AdjustEnds   ( cpstr S, cpstr T, int J, int I );
        void    PrintVT      ( cpstr S, cpstr T );

    };



    //  ======================   Alignment1   =======================

    DefineClass(Alignment1);

    class  Alignment1 : public io::Stream  {

      public :

        Alignment1 ();
        Alignment1 ( io::RPStream Object );
        ~Alignment1();

        void SetAffineModel ( realtype WGap,   realtype WSpace  );
        void SetScores      ( realtype SEqual, realtype SNEqual );

        void Align          ( ivector S, int SLength,
                              ivector T, int TLength,
                              ALIGN_METHOD Method=ALIGN_GLOBAL );

        inline ivector  GetAlignedS   ()  { return AlgnS;     }
        inline ivector  GetAlignedT   ()  { return AlgnT;     }
        inline int      GetAlignLength()  { return AlgnLen;   }
        inline realtype GetScore      ()  { return VAchieved; }

        realtype GetSimilarity(); // Score-weighted sequence id

        virtual void OutputResults ( io::RFile f, ivector S, int lenS,
                                                  ivector T, int lenT );

        void read  ( io::RFile f );
        void write ( io::RFile f );

      protected :

        int      Space;
        int      AlignKey, SLen,TLen, AlgnLen;
        rmatrix  VT,ET,FT;
        ivector  AlgnS,AlgnT;
        realtype VAchieved;
        realtype SEq,SNEq, Wg,Ws;

        virtual void  InitAlignment1();
        virtual void  FreeMemory    ();
        virtual realtype  Score     ( int A, int B );

        void    BuildGATable ( ivector S, ivector T,
                               bool FreeSEnds, bool FreeTEnds );
        void    BuildLATable ( ivector S, ivector T );
        void    Backtrace    ( ivector S, ivector T, int J, int I,
                               bool StopAtZero );
        void    AdjustEnds   ( ivector S, ivector T, int J, int I );
        void    PrintVT      ( ivector S, ivector T );

    };

  }  // namespace math

}  // namespace mmdb

#endif
