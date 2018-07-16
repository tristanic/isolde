//  $Id: mmdb_seqsuperpose.cpp $
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
//  **** Module  :  SeqSuperpose <implementation>
//       ~~~~~~~~~
//  **** Classes :  mmdb::SeqSuperpose
//       ~~~~~~~~~
//
//  (C) E.Krissinel  2005-2013
//
//  =================================================================
//

#include <math.h>
#include <string.h>

#include "mmdb_tables.h"
#include "mmdb_seqsuperpose.h"

namespace mmdb  {

  //  =================================================================

  SeqSuperpose::SeqSuperpose()  {
    SeqSuperposeInit();
  }

  SeqSuperpose::~SeqSuperpose()  {
    FreeMemory();
  }

  void  SeqSuperpose::SeqSuperposeInit()  {
    Align  = NULL;
    Mat4Init ( TMatrix ); // superposes Ca1 over Ca2: |T*Ca1 - Ca2|->min
    Q      = -0.5;        // Q-score
    rmsd   = MaxReal;     // rmsd
    seqId  = MaxReal;     // sequence identity
    _seqId = 0.0;         // sequence identity in sequence alignment
    Nalign = 0;           // alignment length
    c1     = NULL;        // sup-n vector: Ca1[i]->Ca2[c1[i]] if c1[i]>=0
    c2     = NULL;        // sup-n vector: Ca2[i]->Ca1[c2[i]] if c2[i]>=0
    cn1    = NULL;        // temporary contact array #1
    cn2    = NULL;        // temporary contact array #2
    Rmsd0  = 3.0;         // quality optimization parameter
    maxContact = 15.0;    // maximal Calpha-pair contact parameter
    contact    = NULL;
    ncontacts  = 0;
  }

  void  SeqSuperpose::FreeMemory()  {
    if (Align)  {
      delete Align;
      Align = NULL;
    }
    FreeVectorMemory ( c1 ,0 );
    FreeVectorMemory ( c2 ,0 );
    FreeVectorMemory ( cn1,0 );
    FreeVectorMemory ( cn2,0 );
    if (contact)  {
      delete[] contact;
      contact = NULL;
    }
    ncontacts = 0;
  }

  void  makeAAString ( pstr & S, PPAtom C, int nat )  {
  pstr    rname;
  ResName r1;
  int     i,j;
    S = new char[nat+1];
    j = 0;
    for (i=0;i<nat;i++)
      if (C[i])  {
        rname = C[i]->GetResName();
        if (rname)  {
          Get1LetterCode ( rname,r1 );
          S[j++] = r1[0];
        }
      }
    S[j] = char(0);
  }


  realtype  SeqSuperpose::MatchQuality ( int Nalign, realtype Rmsd,
                                         int nres1,  int nres2 )  {
    if (Nalign==0)  return 0.0;
    return MatchQuality2 ( Nalign,Rmsd*Rmsd*Nalign,nres1,nres2 );
  }

  realtype  SeqSuperpose::MatchQuality2 ( int Nalign, realtype dist2,
                                          int nres1,  int nres2 ) {
  realtype  NormN,Na2,NormR;
    NormN = nres1*nres2;
    if (NormN<=0.0) return 0.0;
    Na2   = Nalign*Nalign;
    NormR = dist2/(Nalign*Rmsd0*Rmsd0);
    return  Na2/((1.0+NormR)*NormN);
  }


  void  SeqSuperpose::MakeContacts ( mat44 & TM, realtype cont_est ) {
  //  Find the closest contacts atoms and makes the correspondence
  //  vectors cn1 and cn2
  int  i,j,i1,i2;

    //  1. Find all contacts in the range of 0.0 - cont_est
    if (contact)  {
      delete[] contact;
      contact = NULL;
    }
    ncontacts = 0;
    M->SeekContacts ( Ca2,nCa2,Ca1,nCa1,0.0,cont_est,0,
                      contact,ncontacts,0,&TM,0,
                      BRICK_ON_1 | BRICK_READY );

    //  2. Leave only unique shortest contacts, that is, if Ca1[i]-Ca2[j]
    //     is the shortest contact for atom Ca1[i], it has also to be
    //     the shortest contact for atom Ca2[j].

    if (ncontacts>0)  {

      SortContacts ( contact,ncontacts,CNSORT_DINC );

      for (i=0;i<nCa1;i++)
        cn1[i] = -1;
      for (i=0;i<nCa2;i++)
        cn2[i] = -1;

      j = 0;
      for (i=0;i<ncontacts;i++)  {
        i1 = contact[i].id2;
        i2 = contact[i].id1;
        if ((cn1[i1]<0) && (cn2[i2]<0))  {
          // We only check for unmapped atoms in this version, so that
          // chain misdirection and wide-angle contacts are accepted.
          // However, the method itself is meant to be used only with
          // highly similar chains, so we do not expect difficulties
          // here. Our purpose here is to get maximum performance at
          // high-quality input. See SSM code for more rigorous contact
          // building.
          if (j<i)  contact[j].Copy ( contact[i] );
          // close contact
          cn1[i1] = i2;
          cn2[i2] = i1;
          j++;
        }
      }

      ncontacts = j;

    }

  }

  int  SeqSuperpose::makeStructAlignment ( realtype seqThreshold,
                                           bool  keepBricks )  {
  pstr     S,T;
  mat44    TM;
  realtype dist2,maxRMSD2,Q1,Q0,dist20;
  int      i,i1,i2,nal,rc,iter,iter1;
  char     Space;

    S = Align->GetAlignedS();
    T = Align->GetAlignedT();

    GetVectorMemory ( cn1,nCa1,0 );
    GetVectorMemory ( c1 ,nCa1,0 );
    for (i=0;i<nCa1;i++)  {
      cn1[i] = -1;
      c1 [i] = -1;
    }
    GetVectorMemory ( cn2,nCa2,0 );
    GetVectorMemory ( c2 ,nCa2,0 );
    for (i=0;i<nCa2;i++)  {
      cn2[i] = -1;
      c2 [i] = -1;
    }

    i  = 0;
    i1 = 0;
    i2 = 0;
    Space = Align->GetSpace();
    while (S[i] && (i1<nCa1) && (i2<nCa2))  {
      if ((S[i]==Space) && (T[i]!=Space))      i2++;
      else if ((S[i]!=Space) && (T[i]==Space)) i1++;
      else  {
        if (S[i]==T[i])  {
          cn1[i1] = i2;
          cn2[i2] = i1;
          _seqId += 1.0;
        }
        i1++;
        i2++;
      }
      i++;
    }

    _seqId /= IMax(nCa1,nCa2);

    if (_seqId<seqThreshold)  {
      FreeVectorMemory ( cn1,0 );
      FreeVectorMemory ( cn2,0 );
      return SEQSP_SeqThreshold;
    }

    maxRMSD2 = maxContact*maxContact;

    if ((!keepBricks) || (!M->areBricks()))
      M->MakeBricks ( Ca2,nCa2,1.25*maxContact );

    Q     = -1.0;
    Q0    = -1.0;
    iter  = 0;
    iter1 = 0;

    do  {

      Q = RMax ( Q,Q0 );
      iter++;

      rc = SuperposeAtoms ( TM,Ca1,nCa1,Ca2,cn1 );

      if (rc==SPOSEAT_Ok)  {

        MakeContacts ( TM,maxContact );

        dist2 = 0.0;
        for (i=0;i<ncontacts;i++)  {
          contact[i].dist *= contact[i].dist;
          dist2 += contact[i].dist;
        }

        dist20 = dist2;
        nal    = ncontacts;
        Q0     = RMax ( Q0,MatchQuality2(ncontacts,dist2,nCa1,nCa2) );
        if (ncontacts>0)  {
          i = ncontacts;
          while (i>3)  {
            i--;
            dist2 -= contact[i].dist;
            if (dist2<=i*maxRMSD2)  { // rmsd must be within the limits
              Q1 = MatchQuality2 ( i,dist2,nCa1,nCa2 );
              if (Q1>Q0)  {
                Q0     = Q1;
                nal    = i;
                dist20 = dist2;
              }
            }
          }
          for (i=nal+1;i<ncontacts;i++)  {
            cn1[contact[i].id2] = -1;
            cn2[contact[i].id1] = -1;
          }
          if (Q0>Q)  {
            for (i=0;i<nCa1;i++)
              c1[i] = cn1[i];
            for (i=0;i<nCa2;i++)
              c2[i] = cn2[i];
            Mat4Copy ( TM,TMatrix );
            Nalign = nal;
            rmsd   = dist20;
            iter1  = 0;
          } else
            iter1++;
        }

      }

      if ((!rc) && (iter>100))  rc = SEQSP_IterLimit;

    } while ((rc==SPOSEAT_Ok) && ((Q<Q0) || (iter1>=2)));

    if (Nalign>0)  {
      SuperposeAtoms ( TMatrix,Ca1,nCa1,Ca2,c1 );
      rmsd  = sqrt(rmsd/Nalign);  // rmsd
      seqId = 0.0;
      for (i=0;i<nCa1;i++)
        if (c1[i]>=0)  {
          if (!strcasecmp(Ca1[i]->GetResName(),Ca2[c1[i]]->GetResName()))
            seqId += 1.0;
        }
      seqId = seqId/Nalign;
    } else  {
      rmsd  = MaxReal;
      seqId = 0.0;
    }

    FreeVectorMemory ( cn1,0 );
    FreeVectorMemory ( cn2,0 );

    if (!keepBricks)  M->RemoveBricks();

    return rc;

  }


  int  SeqSuperpose::Superpose ( PManager MMDB,
                                 PPAtom   Calpha1, int nCalpha1,
                                 PPAtom   Calpha2, int nCalpha2,
                                 realtype seqThreshold,
                                 bool     keepBricks )  {
  pstr S,T;

    Mat4Init ( TMatrix ); // superposes Ca1 over Ca2: |T*Ca1 - Ca2|->min
    Q       = 0.0;        // Q-score
    rmsd    = MaxReal;    // rmsd
    seqId   = MaxReal;    // sequence identity in structure alignment
    Nalign  = 0;          // alignment length in structure alignment

    FreeVectorMemory ( c1,0 );
    FreeVectorMemory ( c2,0 );

    _seqId  = IMin(nCalpha1,nCalpha2);
    _seqId /= IMax(nCalpha1,nCalpha2);

    if (_seqId<seqThreshold)
      return SEQSP_SeqThreshold;

    M    = MMDB;
    Ca1  = Calpha1;
    nCa1 = nCalpha1;
    Ca2  = Calpha2;
    nCa2 = nCalpha2;

    makeAAString ( S,Ca1,nCa1 );
    makeAAString ( T,Ca2,nCa2 );

    if (!Align)  Align = new math::Alignment();

    Align->Align ( S,T,math::ALIGN_FREEENDS );

    if (S) delete[] S;
    if (T) delete[] T;

    return makeStructAlignment ( seqThreshold,keepBricks );

  }

}  // namespace mmdb

