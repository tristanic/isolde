//  $Id: mmdb_seqsuperpose.h $
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
//    19.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  SeqSuperpose <interface>
//       ~~~~~~~~~
//  **** Classes :  mmdb::SeqSuperpose
//       ~~~~~~~~~
//
//  (C) E.Krissinel  2005-2013
//
//  =================================================================
//

#ifndef  __Seq_Superpose__
#define  __Seq_Superpose__

#include "mmdb_manager.h"
#include "mmdb_math_align.h"

namespace mmdb  {

  //  =================================================================

  enum SEQSP_RC  {
    SEQSP_Ok           =   0,
    SEQSP_IterLimit    = 100,
    SEQSP_SeqThreshold = 101
  };

  DefineClass(SeqSuperpose);

  class SeqSuperpose  {

    public :
      mat44    TMatrix;  // superposes Ca1 over Ca2: |T*Ca1 - Ca2|->min
      realtype Q;        // Q-score
      realtype rmsd;     // rmsd
      realtype seqId;    // sequence identity in structure alignment
      realtype _seqId;   // sequence identity in sequence alignment
      int      Nalign;   // alignment length in structure alignment
      ivector  c1;       // sup-n vector: Ca1[i]->Ca2[c1[i]] if c1[i]>=0
      ivector  c2;       // sup-n vector: Ca2[i]->Ca1[c2[i]] if c2[i]>=0

      SeqSuperpose();
      ~SeqSuperpose();

      //   Given two sets of atoms, Calpha1 and Calpha2, Superpose(...)
      // calculates the rotational-translational matrix TMatrix such
      // that |TMatrix*Calpha1 - Calpha2| is minimal in least-square
      // terms.
      //   In difference of a full-scale SSM, this simplified version
      // uses initial superposition from sequence alignment, hence
      // it should be applied only to similar chains where calculation
      // time is crucial. seqThreshold specifies a threshold of
      // sequence identity (0<=seqThreshold<=1), below which
      // structural alignment is not performed and Superpose(..)
      // returns SEQSP_SeqThreshold.
      //
      //   If keepBricks is set True, then space bricks are not
      // removed in MMDB and may be used in the next call if
      // vector Calpha2 does not change. This saves computation
      // time.
      //
      //   The alignment results return in public fields above:
      //     TMatrix  - transformation matrix (1 if not aligned)
      //     Q        - quality Q-score (-1 if not aligned)
      //     rmsd     - r.m.s.d (MaxReal if not aligned)
      //     seqId    - sequence identity in structure alignment
      //                        (0 if not aligned)
      //     Nalign   - alignment length in structure alignment
      //                        (0 if not aligned)
      //     c1,c2    - atom corrspondences:
      //                Calpha1[i] <=> Calpha2[c1[i]]
      //                Calpha2[i] <=> Calpha1[c2[i]]
      //
      // Upon success, Superpose(...) returns SEQSP_Ok
      //
      int Superpose ( PManager MMDB,
                      PPAtom   Calpha1, int nCalpha1,
                      PPAtom   Calpha2, int nCalpha2,
                      realtype seqThreshold,
                      bool     keepBricks );

    protected :
      math::PAlignment Align;
      PManager    M;        // pointers to
      PPAtom      Ca1,Ca2;    //   the input data
      int         nCa1,nCa2;  // copy chain lengths
      ivector     cn1,cn2;    // temporary contact arrays
      realtype    Rmsd0;      // quality optimization parameter
      realtype    maxContact; // maximal Calpha-pair contact parameter
      PContact    contact;
      int         ncontacts;

      void SeqSuperposeInit();
      void FreeMemory      ();
      realtype  MatchQuality   ( int Nalign, realtype Rmsd,
                                 int nres1,  int nres2  );
      realtype  MatchQuality2  ( int Nalign, realtype dist2,
                                 int nres1,  int nres2  );
      void MakeContacts        ( mat44 & TM, realtype cont_est );
      int  makeStructAlignment ( realtype seqThreshold,
                                 bool     keepBricks );

  };

}  // namespace mmdb

#endif
