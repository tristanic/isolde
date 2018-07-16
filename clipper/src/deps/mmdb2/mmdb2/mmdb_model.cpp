//  $Id: mmdb_model.cpp $
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
//    11.09.15   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  MMDB_Model <implementation>
//       ~~~~~~~~~
//  **** Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::HetCompound  ( description of het compounds  )
//       ~~~~~~~~~  mmdb::HetCompounds (HETNAM, HETSYN, FORMULA records)
//                  mmdb::SSContainer  (container for helixes and turns)
//                  mmdb::Helix        ( helix info                    )
//                  mmdb::Strand       ( strand info                   )
//                  mmdb::Sheet        ( sheet info                    )
//                  mmdb::Sheets       ( container for sheets          )
//                  mmdb::Turn         ( turn info                     )
//                  mmdb::LinkContainer   ( container for link data    )
//                  mmdb::Link            ( link data                  )
//                  mmdb::LinkRContainer  ( container for refmac link  )
//                  mmdb::LinkR           ( link data                  )
//                  mmdb::CisPepContainer ( container for CisPep data  )
//                  mmdb::CisPep          ( CisPep data                )
//                  mmdb::Model        ( PDB model                     )
//
//  Copyright (C) E. Krissinel 2000-2015
//
//  =================================================================
//

#include <string.h>
#include <stdlib.h>

#include "mmdb_model.h"
#include "mmdb_manager.h"
#include "mmdb_cifdefs.h"

namespace mmdb  {

  //  ===================  HetCompound  =========================

  HetCompound::HetCompound ( cpstr HetName ) : io::Stream()  {
    InitHetCompound ( HetName );
  }

  HetCompound::HetCompound ( io::RPStream Object ) : io::Stream(Object)  {
    InitHetCompound ( pstr("---") );
  }

  HetCompound::~HetCompound() {
    FreeMemory();
  }

  void  HetCompound::InitHetCompound ( cpstr HetName )  {
    strcpy_n0 ( hetID,HetName,sizeof(ResName) );
    comment    = NULL;
    nSynonyms  = 0;
    hetSynonym = NULL;
    compNum    = MinInt4;
    wc         = ' ';
    Formula    = NULL;
  }

  void  HetCompound::FreeMemory()  {
  int i;
    if (comment)  {
      delete[] comment;
      comment = NULL;
    }
    if (hetSynonym)  {
      for (i=0;i<nSynonyms;i++)
        if (hetSynonym[i])  delete[] hetSynonym[i];
      delete[] hetSynonym;
      hetSynonym = NULL;
    }
    nSynonyms = 0;
    if (Formula)  {
      delete[] Formula;
      Formula = NULL;
    }
  }

  void  HetCompound::AddKeyWord ( cpstr W, bool Closed )  {
  psvector HS1;
  int      i;
    if (Closed || (!hetSynonym))  {
      // first synonym orthe previous synonym was closed by semicolon
      // -- add a new one
      HS1 = new pstr[nSynonyms+1];
      for (i=0;i<nSynonyms;i++)
        HS1[i] = hetSynonym[i];
      if (hetSynonym)  delete[] hetSynonym;
      hetSynonym = HS1;
      hetSynonym[nSynonyms] = NULL;
      CreateCopy ( hetSynonym[nSynonyms],W );
      nSynonyms++;
    } else  {
      // just add W to the last synonym
      CreateConcat ( hetSynonym[nSynonyms-1],pstr(" "),W );
    }
  }


  void HetCompound::HETNAM_PDBDump ( io::RFile f )  {
  char S[100];
  pstr p1,p2;
  char c;
  int  N,i;

    if (!comment)  return;

    c = ' ';
    N  = 0;
    p1 = comment;
    do  {
      N++;
      if (N==1)  sprintf ( S,"HETNAM     %3s " ,hetID   );
           else  sprintf ( S,"HETNAM  %2i %3s ",N,hetID );
      while (*p1==' ')  p1++;
      p2 = FirstOccurence(p1,'\n');
      if (p2)  {
        c   = *p2;
        *p2 = char(0);
      } else if (strlen(p1)>53)  {
        i = 0;
        while (p1[i] && (i<53) && (p1[i]!=' '))  i++;
        p2  = &(p1[i]);
        c   = *p2;
        *p2 = char(0);
      }
      if (*p1)  {
        strcat      ( S,p1 );
        PadSpaces   ( S,80 );
        f.WriteLine ( S );
      } else
        N--;
      if (p2)  {
        *p2 = c;
        if (c)  p1 = p2+1;
          else  p2 = NULL;
      }
    } while (p2);

  }


  void HetCompound::HETSYN_PDBDump ( io::RFile f )  {
  char S[100];
  pstr p;
  char c;
  int  N,k,i,l;
    if (!hetSynonym)  return;
    N = 0;
    k = 0;
    p = &(hetSynonym[0][0]);
    do  {
      N++;
      if (N==1)  sprintf ( S,"HETSYN     %3s " ,hetID   );
           else  sprintf ( S,"HETSYN  %2i %3s ",N,hetID );
      i = 0;
      do  {
        l = strlen(p)+2;
        if (i+l<54)  {
          strcat ( S,p );
          if (k<nSynonyms-1) strcat ( S,"; " );
          k++;
          i += l;
          if (k<nSynonyms)  p = &(hetSynonym[k][0]);
                      else  i = 60;  // break loop
        } else  {
          if (i==0)  {
            // too long synonym, has to be split over several lines
            i = l-3;
            while (i>51)  {
              i--;
              while ((i>0) && (p[i]!=' '))  i--;
            }
            if (i<2)  i = 51;  // no spaces!
            c    = p[i];
            p[i] = char(0);
            strcat ( S,p );
            p[i] = c;
            p    = &(p[i]);
            while (*p==' ')  p++;
          }
          i = 60;  // break loop
        }
      } while (i<54);
      PadSpaces ( S,80 );
      f.WriteLine ( S );
    } while (k<nSynonyms);
  }


  void HetCompound::FORMUL_PDBDump ( io::RFile f )  {
  char S[100];
  pstr p1,p2;
  char c;
  int  N,i;
    if (!Formula)  return;
    N  = 0;
    p1 = Formula;
    do  {
      N++;
      if (compNum>MinInt4)  {
        if (N==1)  sprintf ( S,"FORMUL  %2i  %3s    " ,compNum,hetID   );
             else  sprintf ( S,"FORMUL  %2i  %3s %2i ",compNum,hetID,N );
      } else  {
        if (N==1)  sprintf ( S,"FORMUL      %3s    " ,hetID   );
             else  sprintf ( S,"FORMUL      %3s %2i ",hetID,N );
      }
      S[18] = wc;
      p2 = FirstOccurence(p1,'\n');
      if (p2)  {
        c   = *p2;
        *p2 = char(0);
      } else if (strlen(p1)>50)  {
        while (*p1==' ')  p1++;
        i = 0;
        while (p1[i] && (i<50) && (p1[i]!=' '))  i++;
        p2  = &(p1[i]);
        c   = *p2;
        *p2 = char(0);
      }
      strcat ( S,p1 );
      if (p2)  {
        *p2 = c;
        p1  = p2+1;
      }
      PadSpaces ( S,80 );
      f.WriteLine ( S );
    } while (p2);
  }


  void  HetCompound::FormComString ( pstr & F )  {
  pstr p;
  int  i;
    if (F)  {
      delete[] F;
      F = NULL;
    }
    if (comment)  {
      CreateCopy ( F,comment );
      i = 0;
      p = comment;
      while (*p)  {
        p++;
        if (*p=='\n')  i = 0;
                 else  i++;
        if (i>68)  {
          F[i] = char(0);
          CreateConcat ( F,pstr("\n"),p );
          i = 0;
        }
      }
    }
  }


  void  HetCompound::FormSynString ( pstr & F )  {
  pstr p;
  char c;
  int  i,k,l;
    if (F)  {
      delete[] F;
      F = NULL;
    }
    if (hetSynonym)  {
      CreateCopy ( F,pstr("  ") );
      k = 0;
      p = &(hetSynonym[0][0]);
      do  {
        l = strlen(p)+2;
        if (l<=60)  {
          if (k<nSynonyms-1)  CreateConcat ( F,p,pstr(";\n  ") );
                        else  CreateConcat ( F,p );
          k++;
          if (k<nSynonyms)  p = &(hetSynonym[k][0]);
        } else  {
          // too long synonym, has to be split over several lines
          i = l-3;
          while (i>60)  {
            i--;
            while ((i>0) && (p[i]!=' '))  i--;
          }
          if (i<2)  i = 60;  // no spaces!
          c    = p[i];
          p[i] = char(0);
          CreateConcat ( F,p,pstr("\n  ") );
          p[i] = c;
          p    = &(p[i]);
          while (*p==' ')  p++;
        }
      } while (k<nSynonyms);
    }
  }

  void  HetCompound::FormForString ( pstr & F )  {
  pstr p;
  int  i;
    if (F)  {
      delete[] F;
      F = NULL;
    }
    if (Formula)  {
      CreateCopy ( F,Formula );
      i = 0;
      p = &(Formula[0]);
      while (*p)  {
        p++;
        if (*p=='\n')  i = 0;
                 else  i++;
        if (i>68)  {
          F[i] = char(0);
          CreateConcat ( F,pstr("\n"),p );
          i = 0;
        }
      }
    }
  }


  void  HetCompound::Copy ( PHetCompound hetCompound )  {
  int i;
    FreeMemory ();
    strcpy     ( hetID  ,hetCompound->hetID   );
    CreateCopy ( comment,hetCompound->comment );
    nSynonyms = hetCompound->nSynonyms;
    if (nSynonyms>0) {
      hetSynonym = new pstr[nSynonyms];
      for (i=0;i<nSynonyms;i++)  {
        hetSynonym[i] = NULL;
        CreateCopy ( hetSynonym[i],hetCompound->hetSynonym[i] );
      }
    }
    compNum = hetCompound->compNum;
    wc      = hetCompound->wc;
    CreateCopy ( Formula,hetCompound->Formula );
  }

  void  HetCompound::write ( io::RFile f )  {
  int  i;
  byte Version=1;
    f.WriteByte    ( &Version    );
    f.WriteTerLine ( hetID,false );
    f.CreateWrite  ( comment     );
    f.WriteInt     ( &nSynonyms  );
    for (i=0;i<nSynonyms;i++)
      f.CreateWrite ( hetSynonym[i] );
    f.WriteInt    ( &compNum       );
    f.WriteFile   ( &wc,sizeof(wc) );
    f.CreateWrite ( Formula        );
  }

  void  HetCompound::read ( io::RFile f )  {
  int  i;
  byte Version;
    FreeMemory();
    f.ReadByte    ( &Version    );
    f.ReadTerLine ( hetID,false );
    f.CreateRead  ( comment     );
    f.ReadInt     ( &nSynonyms  );
    if (nSynonyms>0) {
      hetSynonym = new pstr[nSynonyms];
      for (i=0;i<nSynonyms;i++)  {
        hetSynonym[i] = NULL;
        f.CreateRead ( hetSynonym[i] );
      }
    }
    f.ReadInt    ( &compNum       );
    f.ReadFile   ( &wc,sizeof(wc) );
    f.CreateRead ( Formula        );
  }

  MakeStreamFunctions(HetCompound)


  //  ====================  HetCompounds  =======================


  HetCompounds::HetCompounds() : io::Stream()  {
    InitHetCompounds();
  }

  HetCompounds::HetCompounds ( io::RPStream Object )
              : io::Stream(Object)  {
    InitHetCompounds();
  }

  HetCompounds::~HetCompounds() {
    FreeMemory();
  }

  void  HetCompounds::InitHetCompounds()  {
    nHets       = 0;
    hetCompound = NULL;
    Closed      = false;
  }

  void  HetCompounds::FreeMemory()  {
  int i;
    if (hetCompound)  {
      for (i=0;i<nHets;i++)
        if (hetCompound[i])  delete hetCompound[i];
      delete[] hetCompound;
      hetCompound = NULL;
    }
    nHets = 0;
  }

  void  HetCompounds::ConvertHETNAM ( cpstr S )  {
  ResName hetID;
  char    L[100];
  int     l,i;
    l = strlen(S);
    if (l>12)  {
      strcpy_n0 ( hetID,&(S[11]),3 );
      i = AddHetName ( hetID );
      if (l>15)  {
        if (hetCompound[i]->comment)  strcpy ( L,"\n" );
                                else  L[0] = char(0);
        strcat       ( L,&(S[15])    );
        CutSpaces    ( L,SCUTKEY_END );
        CreateConcat ( hetCompound[i]->comment,L );
      }
    }
  }

  void  HetCompounds::ConvertHETSYN ( cpstr S )  {
  ResName hetID;
  char    L[100];
  int     l,i,j,k;
    l = strlen(S);
    if (l>12)  {
      strcpy_n0 ( hetID,&(S[11]),3 );
      i = AddHetName ( hetID );
      if (l>15)  {
        j = 15;
        do {
          while (S[j]==' ')  j++;
          k = 0;
          if (S[j])  {
            while (S[j] && (S[j]!=';'))
              L[k++] = S[j++];
            L[k--] = char(0);
            while ((k>0) && (L[k]==' '))  L[k--] = char(0);
            if (L[0])  {
              hetCompound[i]->AddKeyWord ( L,Closed );
              Closed = (S[j]==';');
            }
            if (S[j])  j++;
          }
        } while (S[j]);
        /*
        p1 = &(S[15]);
        do  {
          p2 = FirstOccurence ( p1,';' );
          if (p2)  {
            c   = *p2;
            *p2 = char(0);
          }
          strcpy_css ( L,p1 );
          if (L[0])
            hetCompound[i]->AddKeyWord ( L,Closed );
          if (p2) {
            if (L[0]) Closed = true;
            *p2 = c;
            p1 = p2+1;
          } else if (L[0])
            Closed = false;
        } while (p2);
        */
      }
    }
  }

  void  HetCompounds::ConvertFORMUL ( cpstr S )  {
  ResName hetID;
  char    L[100];
  int     l,i;
    l = strlen(S);
    if (l>13)  {
      strcpy_n0 ( hetID,&(S[12]),3 );
      i = AddHetName ( hetID );
      if (l>18) {
        GetInteger ( hetCompound[i]->compNum,&(S[9]),2 );
        hetCompound[i]->wc = S[18];
        if (strlen(S)>19)  {
          if (hetCompound[i]->Formula)  strcpy ( L,"\n" );
                                  else  L[0] = char(0);
          strcat       ( L,&(S[19])    );
          CutSpaces    ( L,SCUTKEY_END );
          CreateConcat ( hetCompound[i]->Formula,L );
        }
      }
    }
  }
  int  HetCompounds::AddHetName ( cpstr H )  {
  PPHetCompound HC1;
  int            i;
    i = 0;
    while (i<nHets)  {
      if (hetCompound[i])  {
        if (!strcmp(hetCompound[i]->hetID,H))  break;
      }
      i++;
    }
    if (i>=nHets)  {
      HC1 = new PHetCompound[nHets+1];
      for (i=0;i<nHets;i++)
        HC1[i] = hetCompound[i];
      if (hetCompound)  delete[] hetCompound;
      hetCompound = HC1;
      hetCompound[nHets] = new HetCompound ( H );
      i = nHets;
      nHets++;
    }
    return i;
  }

  void HetCompounds::PDBASCIIDump ( io::RFile f )  {
  int  i;

    for (i=0;i<nHets;i++)
      if (hetCompound[i])
        hetCompound[i]->HETNAM_PDBDump ( f );

    for (i=0;i<nHets;i++)
      if (hetCompound[i])
        hetCompound[i]->HETSYN_PDBDump ( f );

    for (i=0;i<nHets;i++)
      if (hetCompound[i])
        hetCompound[i]->FORMUL_PDBDump ( f );

  }


  void  HetCompounds::MakeCIF ( mmcif::PData CIF )  {
  mmcif::PLoop Loop;
  pstr        F;
  int         RC;
  int         i;

    if (!hetCompound)  return;

    RC = CIF->AddLoop ( CIFCAT_CHEM_COMP,Loop );
    if (RC!=mmcif::CIFRC_Ok)  {
      Loop->AddLoopTag ( CIFTAG_ID               );
      Loop->AddLoopTag ( CIFTAG_NAME             );
      Loop->AddLoopTag ( CIFTAG_NDB_SYNONYMS     );
      Loop->AddLoopTag ( CIFTAG_NDB_COMPONENT_NO );
      Loop->AddLoopTag ( CIFTAG_FORMULA          );
    }

    F = NULL;
    for (i=0;i<nHets;i++)
      if (hetCompound[i])  {
        Loop->AddString ( hetCompound[i]->hetID );
        hetCompound[i]->FormComString ( F );
        Loop->AddString ( F );
        hetCompound[i]->FormSynString ( F );
        Loop->AddString ( F );
        if (hetCompound[i]->compNum>MinInt4)
              Loop->AddInteger ( hetCompound[i]->compNum );
        else  Loop->AddNoData  ( mmcif::CIF_NODATA_QUESTION );
        hetCompound[i]->FormForString ( F );
        Loop->AddString ( F );
      }

    if (F)  delete[] F;

  }

  ERROR_CODE HetCompounds::GetCIF ( mmcif::PData CIF )  {
  mmcif::PLoop Loop;
  char         L[100];
  ResName      hetID;
  pstr         F,p1,p2;
  char         c;
  int          RC,i,l,k;

    FreeMemory();
    c = char(0);  // only to supress compiler warnings

    Loop = CIF->GetLoop ( CIFCAT_CHEM_COMP );
    if (!Loop)  return Error_NoError;

    l = Loop->GetLoopLength();
    F = NULL;

    for (i=0;i<l;i++)  {
      CIFGetString    ( hetID,Loop,CIFTAG_ID,i,sizeof(hetID),
                        pstr("---") );
      k = AddHetName  ( hetID );
      Loop->GetString ( hetCompound[k]->comment,CIFTAG_NAME,i,true );
      RC = Loop->GetInteger ( hetCompound[k]->compNum,
                                       CIFTAG_NDB_COMPONENT_NO,i,true );
      if (RC)  hetCompound[i]->compNum = MinInt4;
      Loop->GetString ( hetCompound[k]->Formula,CIFTAG_FORMULA,i,true );
      RC = Loop->GetString ( F,CIFTAG_NDB_SYNONYMS,i,true );
      if ((!RC) && F )  {
        p1 = &(F[0]);
        while (*p1)  {
          if (*p1=='\n')  *p1 = ' ';
          p1++;
        }
        p1 = &(F[0]);
        do  {
          p2 = FirstOccurence ( p1,';' );
          if (p2)  {
            c   = *p2;
            *p2 = char(0);
          }
          strcpy_css ( L,p1 );
          hetCompound[i]->AddKeyWord ( L,true );
          if (p2) {
            *p2 = c;
            p1 = p2+1;
          }
        } while (p2);
      }
      hetCompound[i]->wc = ' ';
    }

  //  CIF->DeleteLoop ( CIFCAT_CHEM_COMP );

    if (F)  delete[] F;

    return Error_NoError;

  }

  void  HetCompounds::Copy ( PHetCompounds HetCompounds )  {
  int i;
    FreeMemory();
    nHets = HetCompounds->nHets;
    if (nHets>0)  {
      hetCompound = new PHetCompound[nHets];
      for (i=0;i<nHets;i++)  {
        hetCompound[i] = new HetCompound ( "" );
        hetCompound[i]->Copy ( HetCompounds->hetCompound[i] );
      }
    }
  }

  void  HetCompounds::write ( io::RFile f )  {
  int  i;
  byte Version=1;
    f.WriteByte ( &Version );
    f.WriteInt  ( &nHets   );
    for (i=0;i<nHets;i++)
      hetCompound[i]->write ( f );
  }

  void  HetCompounds::read ( io::RFile f )  {
  int  i;
  byte Version;
    FreeMemory();
    f.ReadByte ( &Version );
    f.ReadInt  ( &nHets   );
    if (nHets>0)  {
      hetCompound = new PHetCompound[nHets];
      for (i=0;i<nHets;i++)  {
        hetCompound[i] = new HetCompound ( "---" );
        hetCompound[i]->read ( f );
      }
    }
  }

  MakeStreamFunctions(HetCompounds)



  //  ====================  SSContainer  =========================

  PContainerClass SSContainer::MakeContainerClass ( int ClassID )  {
    switch (ClassID)  {
      default :
      case ClassID_Template : return
                          ClassContainer::MakeContainerClass(ClassID);
      case ClassID_Helix    : return new Helix();
      case ClassID_Turn     : return new Turn ();
    }
  }

  MakeStreamFunctions(SSContainer)


  //  ================  Helix  ===================

  Helix::Helix() : ContainerClass()  {
    InitHelix();
  }

  Helix::Helix ( cpstr S ) : ContainerClass()  {
    InitHelix();
    ConvertPDBASCII ( S );
  }

  Helix::Helix ( io::RPStream Object ) : ContainerClass(Object)  {
    InitHelix();
  }

  Helix::~Helix() {
    if (comment)  delete[] comment;
  }

  void  Helix::InitHelix()  {

    serNum = 0;                   // serial number
    strcpy ( helixID    ,"---" ); // helix ID
    strcpy ( initResName,"---" ); // name of the helix's initial residue
    strcpy ( initChainID,""    ); // chain ID for the chain
                                  // containing the helix
    initSeqNum = 0;               // sequence number of the initial
                                  //    residue
    strcpy ( initICode  ,""    ); // insertion code of the initial
                                  //    residue
    strcpy ( endResName ,"---" ); // name of the helix's terminal residue
    strcpy ( endChainID ,""    ); // chain ID for the chain
                                  // containing the helix
    endSeqNum  = 0;               // sequence number of the terminal
                                  //    residue
    strcpy ( endICode   ,""    ); // insertion code of the terminal
                                  //    residue
    helixClass = 0;               // helix class
    comment    = NULL;            // comment about the helix
    length     = 0;               // length of the helix

  }

  void  Helix::PDBASCIIDump ( pstr S, int N )  {
  UNUSED_ARGUMENT(N);
  //  makes the ASCII PDB OBSLTE line number N
  //  from the class' data
    strcpy     ( S,"HELIX" );
    PadSpaces  ( S,80 );
    PutInteger ( &(S[7]) ,serNum     ,3  );
    strcpy_n1  ( &(S[11]),helixID    ,3  );
    strcpy_n1  ( &(S[15]),initResName,3  );
    if (initChainID[0])  S[19] = initChainID[0];
    PutIntIns  ( &(S[21]),initSeqNum ,4,initICode );
    strcpy_n1  ( &(S[27]),endResName ,3  );
    if (endChainID[0])   S[31] = endChainID[0];
    PutIntIns  ( &(S[33]),endSeqNum  ,4,endICode  );
    PutInteger ( &(S[38]),helixClass ,2  );
    if (comment)
      strcpy_n ( &(S[40]),comment    ,30 );
    PutInteger ( &(S[71]),length     ,5  );
  }

  void AddStructConfTags ( mmcif::PLoop Loop )  {
    Loop->AddLoopTag ( CIFTAG_CONF_TYPE_ID               );
    Loop->AddLoopTag ( CIFTAG_ID                         );
    Loop->AddLoopTag ( CIFTAG_PDB_ID                     );
    Loop->AddLoopTag ( CIFTAG_BEG_LABEL_COMP_ID          );
    Loop->AddLoopTag ( CIFTAG_BEG_LABEL_ASYM_ID          );
    Loop->AddLoopTag ( CIFTAG_BEG_LABEL_SEQ_ID           );
    Loop->AddLoopTag ( CIFTAG_NDB_BEG_LABEL_INS_CODE_PDB );
    Loop->AddLoopTag ( CIFTAG_END_LABEL_COMP_ID          );
    Loop->AddLoopTag ( CIFTAG_END_LABEL_ASYM_ID          );
    Loop->AddLoopTag ( CIFTAG_END_LABEL_SEQ_ID           );
    Loop->AddLoopTag ( CIFTAG_NDB_END_LABEL_INS_CODE_PDB );
    Loop->AddLoopTag ( CIFTAG_NDB_HELIX_CLASS_PDB        );
    Loop->AddLoopTag ( CIFTAG_DETAILS                    );
    Loop->AddLoopTag ( CIFTAG_NDB_LENGTH                 );
  }

  #define  HelixTypeID  "HELX_P"

  void  Helix::MakeCIF ( mmcif::PData CIF, int N )  {
  UNUSED_ARGUMENT(N);
  mmcif::PLoop Loop;
  int         RC;
    RC = CIF->AddLoop ( CIFCAT_STRUCT_CONF,Loop );
    if (RC!=mmcif::CIFRC_Ok)
      // the category was (re)created, provide tags
      AddStructConfTags ( Loop );
    Loop->AddString  ( pstr(HelixTypeID) );
    Loop->AddInteger ( serNum      );
    Loop->AddString  ( helixID     );
    Loop->AddString  ( initResName );
    Loop->AddString  ( initChainID );
    Loop->AddInteger ( initSeqNum  );
    Loop->AddString  ( initICode,true );
    Loop->AddString  ( endResName  );
    Loop->AddString  ( endChainID  );
    Loop->AddInteger ( endSeqNum   );
    Loop->AddString  ( endICode ,true );
    Loop->AddInteger ( helixClass  );
    Loop->AddString  ( comment     );
    Loop->AddInteger ( length      );
  }

  ERROR_CODE Helix::ConvertPDBASCII ( cpstr S )  {
  char L[100];
    GetInteger  ( serNum     ,&(S[7]) ,3  );
    strcpy_ncss ( helixID    ,&(S[11]),3  );
    strcpy_ncss ( initResName,&(S[15]),3  );
    strcpy_ncss ( initChainID,&(S[19]),1  );
    GetIntIns   ( initSeqNum,initICode,&(S[21]),4  );
    strcpy_ncss ( endResName ,&(S[27]),3  );
    strcpy_ncss ( endChainID ,&(S[31]),1  );
    GetIntIns   ( endSeqNum ,endICode ,&(S[33]),4  );
    GetInteger  ( helixClass ,&(S[38]),2  );
    strcpy_ncss ( L          ,&(S[40]),30 );
    CreateCopy  ( comment    ,L           );
    GetInteger  ( length     ,&(S[71]),5  );
    return Error_NoError;
  }

  ERROR_CODE Helix::GetCIF ( mmcif::PData CIF, int & n )  {
  mmcif::PLoop Loop;
  int          RC,l;
  pstr         F;
  bool         Done;
  ERROR_CODE   rc;

    Loop = CIF->GetLoop ( CIFCAT_STRUCT_CONF );
    if (!Loop)  {
      n = -1;  // signal to finish processing of this structure
      return Error_EmptyCIF;
    }

    l    = Loop->GetLoopLength();
    Done = n>=l;
    while (!Done) {
      F = Loop->GetString ( CIFTAG_CONF_TYPE_ID,n,RC );
      if ((!RC) && F)  Done = (strcmp(F,HelixTypeID)==0);
                 else  Done = false;
      if (!Done)  {
        n++;
        Done = n>=l;
      }
    }

    if (n>=l)  {
      n = -1;  // finish processing of Helix
      return Error_EmptyCIF;
    }

    Loop->DeleteField ( CIFTAG_CONF_TYPE_ID,n );

    rc = CIFGetInteger ( serNum,Loop,CIFTAG_ID,n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    CIFGetString ( helixID    ,Loop,CIFTAG_PDB_ID,
                               n,sizeof(helixID),pstr("   ") );

    CIFGetString ( initResName,Loop,CIFTAG_BEG_LABEL_COMP_ID,
                               n,sizeof(initResName),pstr("   ") );
    CIFGetString ( initChainID,Loop,CIFTAG_BEG_LABEL_ASYM_ID,
                               n,sizeof(initChainID),pstr("") );
    CIFGetString ( initICode  ,Loop,CIFTAG_NDB_BEG_LABEL_INS_CODE_PDB,
                               n,sizeof(initICode),pstr("") );
    if (CIFGetInteger(initSeqNum,Loop,CIFTAG_BEG_LABEL_SEQ_ID,n))
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    CIFGetString ( endResName,Loop,CIFTAG_END_LABEL_COMP_ID,
                              n,sizeof(endResName),pstr("   ") );
    CIFGetString ( endChainID,Loop,CIFTAG_END_LABEL_ASYM_ID,
                              n,sizeof(endChainID),pstr("") );
    CIFGetString ( endICode  ,Loop,CIFTAG_NDB_END_LABEL_INS_CODE_PDB,
                              n,sizeof(endICode),pstr("") );
    rc = CIFGetInteger(endSeqNum,Loop,CIFTAG_END_LABEL_SEQ_ID,n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    rc = CIFGetInteger(helixClass,Loop,CIFTAG_NDB_HELIX_CLASS_PDB,n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    CreateCopy     ( comment,Loop->GetString(CIFTAG_DETAILS,n,RC));
    Loop->DeleteField ( CIFTAG_DETAILS,n );
    rc = CIFGetInteger ( length,Loop,CIFTAG_NDB_LENGTH,n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    n++;

    return Error_NoError;

  }

  void  Helix::Copy ( PContainerClass Helix )  {
    serNum     = PHelix(Helix)->serNum;
    initSeqNum = PHelix(Helix)->initSeqNum;
    endSeqNum  = PHelix(Helix)->endSeqNum;
    helixClass = PHelix(Helix)->helixClass;
    length     = PHelix(Helix)->length;
    strcpy ( helixID    ,PHelix(Helix)->helixID     );
    strcpy ( initResName,PHelix(Helix)->initResName );
    strcpy ( initChainID,PHelix(Helix)->initChainID );
    strcpy ( initICode  ,PHelix(Helix)->initICode   );
    strcpy ( endResName ,PHelix(Helix)->endResName  );
    strcpy ( endChainID ,PHelix(Helix)->endChainID  );
    strcpy ( endICode   ,PHelix(Helix)->endICode    );
    CreateCopy ( comment,PHelix(Helix)->comment );
  }

  void  Helix::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte ( &Version    );
    f.WriteInt  ( &serNum     );
    f.WriteInt  ( &initSeqNum );
    f.WriteInt  ( &endSeqNum  );
    f.WriteInt  ( &helixClass );
    f.WriteInt  ( &length     );
    f.WriteTerLine ( helixID    ,false );
    f.WriteTerLine ( initResName,false );
    f.WriteTerLine ( initChainID,false );
    f.WriteTerLine ( initICode  ,false );
    f.WriteTerLine ( endResName ,false );
    f.WriteTerLine ( endChainID ,false );
    f.WriteTerLine ( endICode   ,false );
    f.CreateWrite ( comment );
  }

  void  Helix::read  ( io::RFile f ) {
  byte Version;
    f.ReadByte ( &Version );
    f.ReadInt  ( &serNum     );
    f.ReadInt  ( &initSeqNum );
    f.ReadInt  ( &endSeqNum  );
    f.ReadInt  ( &helixClass );
    f.ReadInt  ( &length     );
    f.ReadTerLine ( helixID    ,false );
    f.ReadTerLine ( initResName,false );
    f.ReadTerLine ( initChainID,false );
    f.ReadTerLine ( initICode  ,false );
    f.ReadTerLine ( endResName ,false );
    f.ReadTerLine ( endChainID ,false );
    f.ReadTerLine ( endICode   ,false );
    f.CreateRead ( comment );
  }

  MakeStreamFunctions(Helix)



  //  ================  Strand  =====================

  Strand::Strand () : io::Stream()  {
    InitStrand();
  }

  Strand::Strand ( io::RPStream Object ) : io::Stream(Object)  {
    InitStrand();
  }

  Strand::~Strand() {
  }

  void  Strand::InitStrand()  {
    initSeqNum = MinInt4;
    endSeqNum  = MinInt4;
    sense      = 0;
    curResSeq  = MinInt4;
    prevResSeq = MinInt4;
    strandNo   = 0;
    strcpy ( sheetID    ,"sheet_0"  );
    strcpy ( initResName,"   "      );
    strcpy ( initChainID,""         );
    strcpy ( initICode  ,""         );
    strcpy ( endResName ,"   "      );
    strcpy ( endChainID ,""         );
    strcpy ( endICode   ,""         );
    strcpy ( curAtom    ," "        );
    strcpy ( curResName ,"   "      );
    strcpy ( curChainID ,""         );
    strcpy ( curICode   ,""         );
    strcpy ( prevAtom   ," "        );
    strcpy ( prevResName,"   "      );
    strcpy ( prevChainID,""         );
    strcpy ( prevICode  ,""         );
  }

  void  Strand::PDBASCIIDump ( pstr S )  {
  //   Finishes making the ASCII PDB SHEET line number N
  // from the class' data. Making is initiated by Sheet.

    strcpy_n1  ( &(S[17]),initResName,3 );
    if (initChainID[0])  S[21] = initChainID[0];
    PutIntIns  ( &(S[22]),initSeqNum ,4,initICode );

    strcpy_n1  ( &(S[28]),endResName ,3 );
    if (endChainID[0])   S[32] = endChainID[0];
    PutIntIns  ( &(S[33]),endSeqNum  ,4,endICode  );

    PutInteger ( &(S[38]),sense      ,2 );

    strcpy_n1  ( &(S[41]),curAtom    ,4 );
    strcpy_n1  ( &(S[45]),curResName ,3 );
    if (curChainID[0])   S[49] = curChainID[0];
    PutIntIns  ( &(S[50]),curResSeq  ,4,curICode  );

    strcpy_n1  ( &(S[56]),prevAtom   ,4 );
    strcpy_n1  ( &(S[60]),prevResName,3 );
    if (prevChainID[0])  S[64] = prevChainID[0];
    PutIntIns  ( &(S[65]),prevResSeq ,4,prevICode );

  }


  void  Strand::MakeCIF ( mmcif::PData CIF )  {
  mmcif::PLoop Loop;
  int         RC;

    RC = CIF->AddLoop ( CIFCAT_STRUCT_SHEET_RANGE,Loop );
    if (RC!=mmcif::CIFRC_Ok)  {
      // the category was (re)created, provide tags
      Loop->AddLoopTag ( CIFTAG_SHEET_ID                   );
      Loop->AddLoopTag ( CIFTAG_ID                         );
      Loop->AddLoopTag ( CIFTAG_BEG_LABEL_COMP_ID          );
      Loop->AddLoopTag ( CIFTAG_BEG_LABEL_ASYM_ID          );
      Loop->AddLoopTag ( CIFTAG_BEG_LABEL_SEQ_ID           );
      Loop->AddLoopTag ( CIFTAG_NDB_BEG_LABEL_INS_CODE_PDB );
      Loop->AddLoopTag ( CIFTAG_END_LABEL_COMP_ID          );
      Loop->AddLoopTag ( CIFTAG_END_LABEL_ASYM_ID          );
      Loop->AddLoopTag ( CIFTAG_END_LABEL_SEQ_ID           );
      Loop->AddLoopTag ( CIFTAG_NDB_END_LABEL_INS_CODE_PDB );
    }
    Loop->AddString  ( sheetID     );
    Loop->AddInteger ( strandNo    );
    Loop->AddString  ( initResName );
    Loop->AddString  ( initChainID );
    Loop->AddInteger ( initSeqNum  );
    Loop->AddString  ( initICode,true );
    Loop->AddString  ( endResName  );
    Loop->AddString  ( endChainID  );
    Loop->AddInteger ( endSeqNum   );
    Loop->AddString  ( endICode ,true );

  }


  ERROR_CODE Strand::ConvertPDBASCII ( cpstr S )  {

    GetInteger  ( strandNo   ,&(S[7])  ,3 );
    strcpy_ncss ( sheetID    ,&(S[11]) ,3 );

    strcpy_ncss ( initResName,&(S[17]) ,3 );
    strcpy_ncss ( initChainID,&(S[21]) ,1 );
    GetIntIns   ( initSeqNum ,initICode,&(S[22]),4 );

    strcpy_ncss ( endResName ,&(S[28]) ,3 );
    strcpy_ncss ( endChainID ,&(S[32]) ,1 );
    GetIntIns   ( endSeqNum  ,endICode ,&(S[33]),4 );

    GetInteger  ( sense      ,&(S[38]) ,2 );

    GetString   ( curAtom    ,&(S[41]) ,4 );
    strcpy_ncss ( curResName ,&(S[45]) ,3 );
    strcpy_ncss ( curChainID ,&(S[49]) ,1 );
    GetIntIns   ( curResSeq  ,curICode ,&(S[50]),4 );

    GetString   ( prevAtom   ,&(S[56]) ,4 );
    strcpy_ncss ( prevResName,&(S[60]) ,3 );
    strcpy_ncss ( prevChainID,&(S[64]) ,1 );
    GetIntIns   ( prevResSeq ,prevICode,&(S[65]),4 );

    return Error_NoError;

  }

  int  Strand::GetCIF ( mmcif::PData CIF, cpstr sheet_id )  {
  mmcif::PLoop Loop;
  int         RC,l,i,sNo;
  pstr        F;

    Loop = CIF->GetLoop ( CIFCAT_STRUCT_SHEET_RANGE );
    if (Loop)  {
      l = Loop->GetLoopLength();
      i = 0;
      while (i<l)  {
        F = Loop->GetString ( CIFTAG_SHEET_ID,i,RC );
        if (F && (!RC))  {
          if (!strcmp(F,sheet_id))  {
            strcpy ( sheetID,sheet_id );
            if (CIFGetInteger(sNo,Loop,CIFTAG_ID,i))  return i;
            if (sNo==strandNo)  {
              CIFGetString ( initResName,Loop,CIFTAG_BEG_LABEL_COMP_ID,
                             i,sizeof(initResName),pstr("   ") );
              CIFGetString ( initChainID,Loop,CIFTAG_BEG_LABEL_ASYM_ID,
                             i,sizeof(initChainID),pstr("") );
              CIFGetString ( initICode,Loop,
                             CIFTAG_NDB_BEG_LABEL_INS_CODE_PDB,
                             i,sizeof(initICode),pstr("") );
              if (CIFGetInteger(initSeqNum,Loop,
                             CIFTAG_BEG_LABEL_SEQ_ID,i))
                return i;
              CIFGetString ( endResName,Loop,CIFTAG_END_LABEL_COMP_ID,
                             i,sizeof(endResName),pstr("   ") );
              CIFGetString ( endChainID,Loop,CIFTAG_END_LABEL_ASYM_ID,
                             i,sizeof(endChainID),pstr("") );
              CIFGetString ( endICode  ,Loop,
                             CIFTAG_NDB_END_LABEL_INS_CODE_PDB,
                             i,sizeof(endICode),pstr("") );
              if (CIFGetInteger(endSeqNum,Loop,
                             CIFTAG_END_LABEL_SEQ_ID,i))
                return i;
              Loop->DeleteRow ( i );
              i = l+100;  // break the loop
            }
          }
        }
        i++;
      }
    }

    return 0;

  }

  void  Strand::Copy ( PStrand Strand )  {
    initSeqNum = Strand->initSeqNum;
    endSeqNum  = Strand->endSeqNum;
    sense      = Strand->sense;
    curResSeq  = Strand->curResSeq;
    prevResSeq = Strand->prevResSeq;
    strcpy ( initResName,Strand->initResName );
    strcpy ( initChainID,Strand->initChainID );
    strcpy ( initICode  ,Strand->initICode   );
    strcpy ( endResName ,Strand->endResName  );
    strcpy ( endChainID ,Strand->endChainID  );
    strcpy ( endICode   ,Strand->endICode    );
    strcpy ( curAtom    ,Strand->curAtom     );
    strcpy ( curResName ,Strand->curResName  );
    strcpy ( curChainID ,Strand->curChainID  );
    strcpy ( curICode   ,Strand->curICode    );
    strcpy ( prevAtom   ,Strand->prevAtom    );
    strcpy ( prevResName,Strand->prevResName );
    strcpy ( prevChainID,Strand->prevChainID );
    strcpy ( prevICode  ,Strand->prevICode   );
  }

  void  Strand::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte ( &Version    );
    f.WriteInt  ( &initSeqNum );
    f.WriteInt  ( &endSeqNum  );
    f.WriteInt  ( &sense      );
    f.WriteInt  ( &curResSeq  );
    f.WriteInt  ( &prevResSeq );
    f.WriteTerLine ( initResName,false );
    f.WriteTerLine ( initChainID,false );
    f.WriteTerLine ( initICode  ,false );
    f.WriteTerLine ( endResName ,false );
    f.WriteTerLine ( endChainID ,false );
    f.WriteTerLine ( endICode   ,false );
    f.WriteTerLine ( curAtom    ,false );
    f.WriteTerLine ( curResName ,false );
    f.WriteTerLine ( curChainID ,false );
    f.WriteTerLine ( curICode   ,false );
    f.WriteTerLine ( prevAtom   ,false );
    f.WriteTerLine ( prevResName,false );
    f.WriteTerLine ( prevChainID,false );
    f.WriteTerLine ( prevICode  ,false );
  }

  void  Strand::read  ( io::RFile f ) {
  byte Version;
    f.ReadByte ( &Version    );
    f.ReadInt  ( &initSeqNum );
    f.ReadInt  ( &endSeqNum  );
    f.ReadInt  ( &sense      );
    f.ReadInt  ( &curResSeq  );
    f.ReadInt  ( &prevResSeq );
    f.ReadTerLine ( initResName,false );
    f.ReadTerLine ( initChainID,false );
    f.ReadTerLine ( initICode  ,false );
    f.ReadTerLine ( endResName ,false );
    f.ReadTerLine ( endChainID ,false );
    f.ReadTerLine ( endICode   ,false );
    f.ReadTerLine ( curAtom    ,false );
    f.ReadTerLine ( curResName ,false );
    f.ReadTerLine ( curChainID ,false );
    f.ReadTerLine ( curICode   ,false );
    f.ReadTerLine ( prevAtom   ,false );
    f.ReadTerLine ( prevResName,false );
    f.ReadTerLine ( prevChainID,false );
    f.ReadTerLine ( prevICode  ,false );
  }

  MakeStreamFunctions(Strand)




  //  ================  Sheet  ===================

  Sheet::Sheet() : io::Stream()  {
    InitSheet();
  }

  Sheet::Sheet ( io::RPStream Object ) : io::Stream(Object)  {
    InitSheet();
  }

  Sheet::~Sheet()  {
    FreeMemory();
  }

  void  Sheet::InitSheet()  {
    nStrands   = 0;
    sheetID[0] = char(0);
    strand     = NULL;
  }

  void  Sheet::FreeMemory()  {
  int i;
    if (strand)  {
      for (i=0;i<nStrands;i++)
        if (strand[i])  delete strand[i];
      delete[] strand;
      strand = NULL;
    }
    nStrands   = 0;
    sheetID[0] = char(0);
  }

  void  Sheet::PDBASCIIDump ( io::RFile f )  {
  char  S[100];
  int   i;
    if (strand)
      for (i=0;i<nStrands;i++)
        if (strand[i])  {
          strcpy      ( S,"SHEET"           );
          PadSpaces   ( S,80                );
          PutInteger  ( &(S[7]) ,i+1     ,3 );
          strcpy_n1   ( &(S[11]),sheetID ,3 );
          PutInteger  ( &(S[14]),nStrands,2 );
          strand[i]->PDBASCIIDump ( S       );
          f.WriteLine ( S                   );
        }
  }

  void  Sheet::OrderSheet()  {
  int       i,k;
  PPStrand  strand1;
    k = 0;
    for (i=0;i<nStrands;i++)
      if (strand[i])  k++;
    if (k<nStrands)  {
      strand1 = new PStrand[k];
      k = 0;
      for (i=0;i<nStrands;i++)
        if (strand[i])  strand1[k++] = strand[i];
      if (strand)  delete[] strand;
      strand   = strand1;
      nStrands = k;
    }
  }

  void  Sheet::MakeCIF ( mmcif::PData CIF )  {
  mmcif::PLoop Loop;
  int          RC;
  int          i;
  bool         isSense;

    OrderSheet();

    RC = CIF->AddLoop ( CIFCAT_STRUCT_SHEET,Loop );
    if (RC!=mmcif::CIFRC_Ok)  {
      // the category was (re)created, provide tags
      Loop->AddLoopTag ( CIFTAG_SHEET_ID       );
      Loop->AddLoopTag ( CIFTAG_NUMBER_STRANDS );
    }
    Loop->AddString  ( sheetID  );
    Loop->AddInteger ( nStrands );

    for (i=0;i<nStrands;i++)  {
      strand[i]->MakeCIF ( CIF );
      if (strand[i]->sense!=0)  isSense = true;
    }

    if (nStrands>1)  {

      if (isSense)  {
        RC = CIF->AddLoop ( CIFCAT_STRUCT_SHEET_ORDER,Loop );
        if (RC!=mmcif::CIFRC_Ok)  {
          // the category was (re)created, provide tags
          Loop->AddLoopTag ( CIFTAG_SHEET_ID   );
          Loop->AddLoopTag ( CIFTAG_RANGE_ID_1 );
          Loop->AddLoopTag ( CIFTAG_RANGE_ID_2 );
          Loop->AddLoopTag ( CIFTAG_SENSE      );
        }
        for (i=1;i<nStrands;i++)  {
          Loop->AddString  ( sheetID               );
          Loop->AddInteger ( strand[i-1]->strandNo );
          Loop->AddInteger ( strand[i]  ->strandNo );
          if (strand[i]->sense>0)
                Loop->AddString ( pstr("parallel")      );
          else  Loop->AddString ( pstr("anti-parallel") );
        }
      }

      RC = CIF->AddLoop ( CIFCAT_STRUCT_SHEET_HBOND,Loop );
      if (RC!=mmcif::CIFRC_Ok)  {
        // the category was (re)created, provide tags
        Loop->AddLoopTag ( CIFTAG_SHEET_ID                       );
        Loop->AddLoopTag ( CIFTAG_RANGE_ID_1                     );
        Loop->AddLoopTag ( CIFTAG_RANGE_ID_2                     );
        Loop->AddLoopTag ( CIFTAG_RANGE_1_BEG_LABEL_ATOM_ID      );
        Loop->AddLoopTag ( CIFTAG_NDB_RANGE_1_BEG_LABEL_COMP_ID  );
        Loop->AddLoopTag ( CIFTAG_NDB_RANGE_1_BEG_LABEL_ASYM_ID  );
        Loop->AddLoopTag ( CIFTAG_RANGE_1_BEG_LABEL_SEQ_ID       );
        Loop->AddLoopTag ( CIFTAG_NDB_RANGE_1_BEG_LABEL_INS_CODE );
        Loop->AddLoopTag ( CIFTAG_RANGE_1_END_LABEL_ATOM_ID      );
        Loop->AddLoopTag ( CIFTAG_NDB_RANGE_1_END_LABEL_COMP_ID  );
        Loop->AddLoopTag ( CIFTAG_NDB_RANGE_1_END_LABEL_ASYM_ID  );
        Loop->AddLoopTag ( CIFTAG_RANGE_1_END_LABEL_SEQ_ID       );
        Loop->AddLoopTag ( CIFTAG_NDB_RANGE_1_END_LABEL_INS_CODE );
      }
      for (i=1;i<nStrands;i++)  {
        Loop->AddString  ( sheetID                );
        Loop->AddInteger ( strand[i-1]->strandNo  );
        Loop->AddInteger ( strand[i]->strandNo    );
        Loop->AddString  ( strand[i]->curAtom     );
        Loop->AddString  ( strand[i]->curResName  );
        Loop->AddString  ( strand[i]->curChainID  );
        Loop->AddInteger ( strand[i]->curResSeq   );
        Loop->AddString  ( strand[i]->curICode ,true );
        Loop->AddString  ( strand[i]->prevAtom    );
        Loop->AddString  ( strand[i]->prevResName );
        Loop->AddString  ( strand[i]->prevChainID );
        Loop->AddInteger ( strand[i]->prevResSeq  );
        Loop->AddString  ( strand[i]->prevICode,true );
      }
    }

  }


  ERROR_CODE Sheet::ConvertPDBASCII ( cpstr S )  {
  int       i,k,ns;
  SheetID   SID;
  PPStrand  strand1;

    GetInteger  ( k  ,&(S[7]) ,3 );
    strcpy_ncss ( SID,&(S[11]),3 );
    GetInteger  ( ns ,&(S[14]),2 );

  //  if (!SID[0])  return  Error_NoSheetID;
    if (!sheetID[0])  strcpy ( sheetID,SID );
    else if (strcmp(sheetID,SID))
                  return  Error_WrongSheetID;

    if (k<=0)     return  Error_WrongStrandNo;

    ns = IMax(k,ns);
    if (!strand)  {
      strand = new PStrand[ns];
      for (i=0;i<ns;i++)
        strand[i] = NULL;
    } else if (ns>nStrands)  {
      strand1 = new PStrand[ns];
      for (i=0;i<nStrands;i++)
        strand1[i] = strand[i];
      for (i=nStrands;i<ns;i++)
        strand1[i] = NULL;
      if (strand)  delete[] strand;
      strand = strand1;
    }
    nStrands = ns;

    k--;
    if (!strand[k])  strand[k] = new Strand();

    return  strand[k]->ConvertPDBASCII ( S );

  }

  void  Sheet::TryStrand ( int strand_no )  {
  int      i,k;
  PPStrand strand1;
    k = -1;
    for (i=0;(i<nStrands) && (k<0);i++)
      if (strand[i])
        if (strand[i]->strandNo==strand_no)  k = i;
    if (k<0)  {
      strand1 = new PStrand[nStrands+1];
      for (i=0;i<nStrands;i++)
        strand1[i] = strand[i];
      if (strand) delete[] strand;
      strand = strand1;
      strand[nStrands] = new Strand();
      strand[nStrands]->strandNo = strand_no;
      nStrands++;
    }
  }


  void  Sheet::CIFFindStrands ( mmcif::PData CIF, cpstr Category ) {
  // just look for all strands mentioned for the sheet
  mmcif::PLoop Loop;
  pstr         F;
  int          RC,i,l,sNo;
    Loop = CIF->GetLoop ( Category );
    if (Loop)  {
      l = Loop->GetLoopLength();
      for (i=0;i<l;i++)  {
        F = Loop->GetString ( CIFTAG_SHEET_ID,i,RC );
        if (F && (!RC))  {
          if (!strcmp(F,sheetID))  {
            if (!Loop->GetInteger(sNo,CIFTAG_ID,i))
              TryStrand ( sNo );
            if (!Loop->GetInteger(sNo,CIFTAG_RANGE_ID_1,i))
              TryStrand ( sNo );
            if (!Loop->GetInteger(sNo,CIFTAG_RANGE_ID_2,i))
              TryStrand ( sNo );
          }
        }
      }
    }
  }

  int  Sheet::GetStrand ( int strand_no )  {
  int i;
    for (i=0;i<nStrands;i++)
      if (strand[i])  {
        if (strand[i]->strandNo==strand_no)
          return i;
      }
    return -1;
  }

  int Sheet::GetCIF ( mmcif::PData CIF )  {
  mmcif::PLoop Loop;
  int         i,ns,l,k,k2,RC,sNo;
  pstr        F;
  ivector     pair;
  bool     Ok;

    pair = NULL;

    //    First find all strands and create
    // the corresponding classes. The CIF fields
    // are not removed at this stage.

    CIFFindStrands ( CIF,CIFCAT_STRUCT_SHEET_ORDER );
    CIFFindStrands ( CIF,CIFCAT_STRUCT_SHEET_RANGE );
    CIFFindStrands ( CIF,CIFCAT_STRUCT_SHEET_HBOND );

    //  Check number of strands
    Loop = CIF->GetLoop ( CIFCAT_STRUCT_SHEET );
    if (Loop)  {
      l = Loop->GetLoopLength();
      i = 0;
      while (i<l)  {
        F = Loop->GetString ( CIFTAG_SHEET_ID,i,RC );
        if (F && (!RC))  {
          if (!strcmp(F,sheetID))  {
            RC = CIFGetInteger1 ( ns,Loop,CIFTAG_NUMBER_STRANDS,i );
            if ((!RC) && (ns!=nStrands))
              return  Error_WrongNumberOfStrands;
            Loop->DeleteRow ( i );
            i = l+100;  // break loop
          }
        }
        i++;
      }
    }

    //  Read each strand
    RC = 0;
    for (i=0;(i<nStrands) && (!RC);i++)
      RC = strand[i]->GetCIF ( CIF,sheetID );

    if (RC)  return RC;

    if (nStrands>1)  {

      GetVectorMemory ( pair,nStrands,0 );
      for (i=0;i<nStrands;i++)
        pair[i] = -1;

      Loop = CIF->GetLoop ( CIFCAT_STRUCT_SHEET_ORDER );
      if (Loop)  {
        Ok = true;
        l  = Loop->GetLoopLength();
        for (i=0;(i<l) && Ok;i++)  {
          F = Loop->GetString ( CIFTAG_SHEET_ID,i,RC );
          if (F && (!RC))  {
            if (!strcmp(F,sheetID))  {
              if (!Loop->GetInteger(sNo,CIFTAG_RANGE_ID_1,i))  {
                k = GetStrand ( sNo );
                if ((k>=0) &&
                    (!Loop->GetInteger(sNo,CIFTAG_RANGE_ID_2,i)))  {
                  pair[k] = GetStrand ( sNo );
                  if (pair[k]>=0)  {
                    F = Loop->GetString ( CIFTAG_SENSE,i,RC );
                    if (F && (!RC))  {
                      if (!strcasecmp(F,"anti-parallel"))
                        strand[pair[k]]->sense = -1;
                      else if (!strcasecmp(F,"parallel"))
                        strand[pair[k]]->sense =  1;
                    }
                    Loop->DeleteRow ( i );
                  } else
                    Ok = false;
                } else
                  Ok = false;
              } else
                Ok = false;
            }
          }
        }
        if (!Ok)  {
          FreeVectorMemory ( pair,0 );
          return Error_WrongSheetOrder;
        }
      }

      Loop = CIF->GetLoop ( CIFCAT_STRUCT_SHEET_HBOND );
      if (Loop)  {
        Ok = true;
        l  = Loop->GetLoopLength();
        for (i=0;(i<l) && Ok;i++)  {
          F = Loop->GetString ( CIFTAG_SHEET_ID,i,RC );
          if (F && (!RC))  {
            if (!strcmp(F,sheetID))  {
              if (!Loop->GetInteger(sNo,CIFTAG_RANGE_ID_1,i))  {
                k = GetStrand ( sNo );
                if ((k>=0) &&
                    (!Loop->GetInteger(sNo,CIFTAG_RANGE_ID_1,i)))  {
                  k2 = GetStrand ( sNo );
                  if (k2>=0)  {
                    if (pair[k]==k2)  {
                      CIFGetString ( strand[k2]->curAtom,Loop,
                                CIFTAG_RANGE_1_BEG_LABEL_ATOM_ID,
                                i,sizeof(strand[k2]->curAtom),
                                pstr("    ") );
                      CIFGetString ( strand[k2]->curResName,Loop,
                                CIFTAG_NDB_RANGE_1_BEG_LABEL_COMP_ID,
                                i,sizeof(strand[k2]->curResName),
                                pstr("   ") );
                      CIFGetString ( strand[k2]->curChainID,Loop,
                                CIFTAG_NDB_RANGE_1_BEG_LABEL_ASYM_ID,
                                i,sizeof(strand[k2]->curChainID),
                                pstr(" ") );
                      if (CIFGetInteger(strand[k2]->curResSeq,Loop,
                                CIFTAG_RANGE_1_BEG_LABEL_SEQ_ID,i))  {
                        FreeVectorMemory ( pair,0 );
                        return i;
                      }
                      CIFGetString ( strand[k2]->curICode,Loop,
                                CIFTAG_NDB_RANGE_1_BEG_LABEL_INS_CODE,
                                i,sizeof(strand[k2]->curICode),
                                pstr(" ") );
                      CIFGetString ( strand[k2]->prevAtom,Loop,
                                CIFTAG_RANGE_1_END_LABEL_ATOM_ID,
                                i,sizeof(strand[k2]->prevAtom),
                                pstr("    ") );
                      CIFGetString ( strand[k2]->prevResName,Loop,
                                CIFTAG_NDB_RANGE_1_END_LABEL_COMP_ID,
                                i,sizeof(strand[k2]->prevResName),
                                pstr("   ") );
                      CIFGetString ( strand[k2]->prevChainID,Loop,
                                CIFTAG_NDB_RANGE_1_END_LABEL_ASYM_ID,
                                i,sizeof(strand[k2]->prevChainID),
                                pstr(" ") );
                      if (CIFGetInteger(strand[k2]->prevResSeq,Loop,
                                CIFTAG_RANGE_1_END_LABEL_SEQ_ID,i))  {
                        FreeVectorMemory ( pair,0 );
                        return i;
                      }
                      CIFGetString ( strand[k2]->prevICode,Loop,
                                CIFTAG_NDB_RANGE_1_END_LABEL_INS_CODE,
                                i,sizeof(strand[k2]->prevICode),
                                pstr(" ") );
                      Loop->DeleteRow ( i );
                    } else
                        Ok = false;
                  } else
                    Ok = false;
                } else
                  Ok = false;
              } else
                Ok = false;
            }
          }
        }
        if (!Ok)  {
          FreeVectorMemory ( pair,0 );
          return Error_HBondInconsistency;
        }
      }
    }

    FreeVectorMemory ( pair,0 );

    return 0;

  }


  void  Sheet::Copy ( PSheet Sheet )  {
  int i;
    FreeMemory();
    nStrands = Sheet->nStrands;
    if (nStrands>0)  {
      strand = new PStrand[nStrands];
      for (i=0;i<nStrands;i++)
        if (Sheet->strand[i])  {
          strand[i] = new Strand();
          strand[i]->Copy ( Sheet->strand[i] );
        } else
          strand[i] = NULL;
    }
    strcpy ( sheetID,Sheet->sheetID );
  }

  void  Sheet::write ( io::RFile f )  {
  int  i;
  byte Version=1;
    f.WriteByte ( &Version  );
    f.WriteInt  ( &nStrands );
    for (i=0;i<nStrands;i++)
      StreamWrite ( f,strand[i] );
    f.WriteTerLine ( sheetID,false );
  }

  void  Sheet::read  ( io::RFile f ) {
  int  i;
  byte Version;
    FreeMemory();
    f.ReadByte ( &Version  );
    f.ReadInt  ( &nStrands );
    if (nStrands>0)  {
      strand = new PStrand[nStrands];
      for (i=0;i<nStrands;i++)  {
        strand[i] = NULL;
        StreamRead ( f,strand[i] );
      }
    }
    f.ReadTerLine ( sheetID,false );
  }

  MakeStreamFunctions(Sheet)



  //  ====================  Sheets  ============================


  Sheets::Sheets() : io::Stream()  {
    InitSheets();
  }


  Sheets::Sheets ( io::RPStream Object ) : io::Stream ( Object )  {
    InitSheets();
  }


  Sheets::~Sheets()  {
    FreeMemory();
  }


  void  Sheets::InitSheets()  {
    nSheets = 0;
    sheet   = NULL;
  }


  void  Sheets::FreeMemory()  {
  int i;
    if (sheet)  {
      for (i=0;i<nSheets;i++)
        if (sheet[i])  delete sheet[i];
      delete[] sheet;
      sheet = NULL;
    }
    nSheets = 0;
  }


  void  Sheets::PDBASCIIDump ( io::RFile f )  {
  int i;
    if (sheet)
      for (i=0;i<nSheets;i++)
        if (sheet[i])  sheet[i]->PDBASCIIDump ( f );
  }


  void  Sheets::MakeCIF ( mmcif::PData CIF )  {
  int i;
    if (sheet)
      for (i=0;i<nSheets;i++)
        if (sheet[i])  sheet[i]->MakeCIF ( CIF );
  }


  ERROR_CODE  Sheets::ConvertPDBASCII ( cpstr S )  {
  PPSheet  sheet1;
  SheetID  sheetID;
  int      i,k;
    strcpy_ncss ( sheetID,&(S[11]),3 );
    //  if (!sheetID[0]) return  Error_NoSheetID;
    k = -1;
    for (i=0;i<nSheets;i++)
      if (sheet[i])  {
        if (!strcmp(sheetID,sheet[i]->sheetID))  {
          k = i;
          break;
        }
      }
    if (k<0)  {
      sheet1 = new PSheet[nSheets+1];
      for (i=0;i<nSheets;i++)
        sheet1[i] = sheet[i];
      if (sheet) delete[] sheet;
      sheet = sheet1;
      sheet[nSheets] = new Sheet();
      k = nSheets;
      nSheets++;
    }
    return  sheet[k]->ConvertPDBASCII ( S );
  }


  void  Sheets::CIFFindSheets ( mmcif::PData CIF, cpstr Category ) {
  mmcif::PLoop Loop;
  int          RC,i,j,k,l;
  pstr         F;
  PPSheet      sheet1;
    Loop = CIF->GetLoop ( Category );
    if (Loop)  {
      l = Loop->GetLoopLength();
      for (i=0;i<l;i++)  {
        F = Loop->GetString ( CIFTAG_SHEET_ID,i,RC );
        if (F && (!RC))  {
          k = -1;
          j = 0;
          while ((j<nSheets) && (k<0))  {
            if (sheet[j])  {
              if (!strcmp(F,sheet[j]->sheetID))  k = j;
            }
            j++;
          }
          if (k<0)  {
            sheet1 = new PSheet[nSheets+1];
            for (i=0;i<nSheets;i++)
              sheet1[i] = sheet[i];
            if (sheet) delete[] sheet;
            sheet = sheet1;
            sheet[nSheets] = new Sheet();
            strcpy ( sheet[nSheets]->sheetID,F );
            nSheets++;
          }
        }
      }
    }
  }

  int Sheets::GetCIF ( mmcif::PData CIF )  {
  int i,RC;

    FreeMemory();

    //  First find all sheet names and create
    // the corresponding classes. The CIF fields
    // are not removed at this stage.

    CIFFindSheets ( CIF,CIFCAT_STRUCT_SHEET       );
    CIFFindSheets ( CIF,CIFCAT_STRUCT_SHEET_ORDER );
    CIFFindSheets ( CIF,CIFCAT_STRUCT_SHEET_RANGE );
    CIFFindSheets ( CIF,CIFCAT_STRUCT_SHEET_HBOND );

    //  Read each sheet
    i  = 0;
    RC = 0;
    while ((i<nSheets) && (!RC))  {
      RC = sheet[i]->GetCIF ( CIF );
      i++;
    }

    return RC;

  }


  void  Sheets::Copy ( PSheets Sheets )  {
  int i;
    FreeMemory();
    if (Sheets->nSheets>0)  {
      nSheets = Sheets->nSheets;
      sheet = new PSheet[nSheets];
      for (i=0;i<nSheets;i++)
        if (Sheets->sheet[i]) {
          sheet[i] = new Sheet();
          sheet[i]->Copy ( Sheets->sheet[i] );
        } else
          sheet[i] = NULL;
    }
  }


  void  Sheets::write ( io::RFile f )  {
  int  i;
  byte Version=1;
    f.WriteByte ( &Version );
    f.WriteInt  ( &nSheets );
    for (i=0;i<nSheets;i++)
      StreamWrite ( f,sheet[i] );
  }


  void  Sheets::read ( io::RFile f )  {
  int  i;
  byte Version;
    FreeMemory();
    f.ReadByte ( &Version );
    f.ReadInt  ( &nSheets );
    if (nSheets>0)  {
      sheet = new PSheet[nSheets];
      for (i=0;i<nSheets;i++)  {
        sheet[i] = NULL;
        StreamRead ( f,sheet[i] );
      }
    }
  }


  MakeStreamFunctions(Sheets)



  //  ================  Turn  ===================

  Turn::Turn() : ContainerClass()  {
    InitTurn();
  }

  Turn::Turn ( cpstr S ) : ContainerClass()  {
    InitTurn();
    ConvertPDBASCII ( S );
  }

  Turn::Turn ( io::RPStream Object ) : ContainerClass(Object)  {
    InitTurn();
  }

  Turn::~Turn() {
    if (comment)  delete[] comment;
  }

  void  Turn::InitTurn()  {
    serNum = 0;                   // serial number
    strcpy ( turnID     ,"---" ); // turn ID
    strcpy ( initResName,"---" ); // name of the turn's initial residue
    strcpy ( initChainID," "   ); // chain ID for the chain
                                  // containing the turn
    initSeqNum = 0;               // sequence number of the initial
                                  //    residue
    strcpy ( initICode  ," "   ); // insertion code of the initial
                                  //    residue
    strcpy ( endResName ,"---" ); // name of the turn's terminal residue
    strcpy ( endChainID ," "   ); // chain ID for the chain
                                  // containing the turn
    endSeqNum  = 0;               // sequence number of the terminal
                                  //    residue
    strcpy ( endICode   ," "   ); // insertion code of the terminal
                                  //    residue
    comment    = NULL;            // comment about the helix
  }

  void  Turn::PDBASCIIDump ( pstr S, int N )  {
  UNUSED_ARGUMENT(N);
  //  makes the ASCII PDB OBSLTE line number N
  //  from the class' data
    strcpy     ( S,"TURN" );
    PadSpaces  ( S,80 );
    PutInteger ( &(S[7]) ,serNum     ,3  );
    strcpy_n1  ( &(S[11]),turnID     ,3  );
    strcpy_n1  ( &(S[15]),initResName,3  );
    strcpy_n1  ( &(S[19]),initChainID,1  );
    PutIntIns  ( &(S[20]),initSeqNum ,4,initICode );
    strcpy_n1  ( &(S[26]),endResName ,3  );
    strcpy_n1  ( &(S[30]),endChainID ,1  );
    PutIntIns  ( &(S[31]),endSeqNum  ,4,endICode  );
    if (comment)
      strcpy_n ( &(S[40]),comment    ,30 );
  }


  #define  TurnTypeID  "TURN_P"

  void  Turn::MakeCIF ( mmcif::PData CIF, int N )  {
  UNUSED_ARGUMENT(N);
  mmcif::PLoop Loop;
  int         RC;
    RC = CIF->AddLoop ( CIFCAT_STRUCT_CONF,Loop );
    if (RC!=mmcif::CIFRC_Ok)
      // the category was (re)created, provide tags
      AddStructConfTags ( Loop );
    Loop->AddString  ( pstr(TurnTypeID) );
    Loop->AddInteger ( serNum      );
    Loop->AddString  ( turnID      );
    Loop->AddString  ( initResName );
    Loop->AddString  ( initChainID );
    Loop->AddInteger ( initSeqNum  );
    Loop->AddString  ( initICode,true );
    Loop->AddString  ( endResName  );
    Loop->AddString  ( endChainID  );
    Loop->AddInteger ( endSeqNum   );
    Loop->AddString  ( endICode ,true );
    Loop->AddNoData  ( mmcif::CIF_NODATA_QUESTION );
    Loop->AddString  ( comment     );
    Loop->AddNoData  ( mmcif::CIF_NODATA_QUESTION );
  }

  ERROR_CODE Turn::ConvertPDBASCII ( cpstr S )  {
  char L[100];
    GetInteger   ( serNum     ,&(S[7]) ,3  );
    strcpy_ncss  ( turnID     ,&(S[11]),3  );
    strcpy_ncss  ( initResName,&(S[15]),3  );
    strcpy_ncss  ( initChainID,&(S[19]),1  );
    GetIntIns    ( initSeqNum,initICode,&(S[20]),4 );
    strcpy_ncss  ( endResName ,&(S[26]),3  );
    strcpy_ncss  ( endChainID ,&(S[30]),1  );
    GetIntIns    ( endSeqNum ,endICode ,&(S[31]),4 );
    strcpy_ncss  ( L          ,&(S[40]),30 );
    CreateCopy   ( comment    ,L           );
    return Error_NoError;
  }

  ERROR_CODE Turn::GetCIF ( mmcif::PData CIF, int & n )  {
  mmcif::PLoop Loop;
  int          RC,l;
  pstr         F;
  bool         Done;
  ERROR_CODE   rc;

    Loop = CIF->GetLoop ( CIFCAT_STRUCT_CONF );
    if (!Loop)  {
      n = -1;  // signal to finish processing of this structure
      return Error_EmptyCIF;
    }

    l    = Loop->GetLoopLength();
    Done = n>=l;
    while (!Done) {
      F = Loop->GetString ( CIFTAG_CONF_TYPE_ID,n,RC );
      if ((!RC) && F)  Done = (strcmp(F,TurnTypeID)==0);
                 else  Done = false;
      if (!Done)  {
        n++;
        Done = n>=l;
      }
    }

    if (n>=l)  {
      n = -1;  // finish processing of Turn
      return Error_EmptyCIF;
    }

    Loop->DeleteField ( CIFTAG_CONF_TYPE_ID,n );

    rc = CIFGetInteger ( serNum,Loop,CIFTAG_ID,n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    CIFGetString ( turnID,Loop,CIFTAG_PDB_ID,n,
                   sizeof(turnID),pstr("   ") );

    CIFGetString ( initResName,Loop,CIFTAG_BEG_LABEL_COMP_ID,
                               n,sizeof(initResName),pstr("   ") );
    CIFGetString ( initChainID,Loop,CIFTAG_BEG_LABEL_ASYM_ID,
                               n,sizeof(initChainID),pstr(" ") );
    CIFGetString ( initICode  ,Loop,CIFTAG_NDB_BEG_LABEL_INS_CODE_PDB,
                               n,sizeof(initICode),pstr(" ") );
    rc = CIFGetInteger ( initSeqNum,Loop,CIFTAG_BEG_LABEL_SEQ_ID,n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    CIFGetString ( endResName,Loop,CIFTAG_END_LABEL_COMP_ID,
                              n,sizeof(endResName),pstr("   ") );
    CIFGetString ( endChainID,Loop,CIFTAG_END_LABEL_ASYM_ID,
                              n,sizeof(endChainID),pstr(" ") );
    CIFGetString ( endICode  ,Loop,CIFTAG_NDB_END_LABEL_INS_CODE_PDB,
                              n,sizeof(endICode),pstr(" ") );
    rc = CIFGetInteger ( endSeqNum,Loop,CIFTAG_END_LABEL_SEQ_ID,n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    CreateCopy      ( comment,Loop->GetString(CIFTAG_DETAILS,n,RC));
    Loop->DeleteField ( CIFTAG_DETAILS,n );

    n++;

    return Error_NoError;

  }

  void  Turn::Copy ( PContainerClass Turn )  {
    serNum     = PTurn(Turn)->serNum;
    initSeqNum = PTurn(Turn)->initSeqNum;
    endSeqNum  = PTurn(Turn)->endSeqNum;
    strcpy ( turnID     ,PTurn(Turn)->turnID      );
    strcpy ( initResName,PTurn(Turn)->initResName );
    strcpy ( initChainID,PTurn(Turn)->initChainID );
    strcpy ( initICode  ,PTurn(Turn)->initICode   );
    strcpy ( endResName ,PTurn(Turn)->endResName  );
    strcpy ( endChainID ,PTurn(Turn)->endChainID  );
    strcpy ( endICode   ,PTurn(Turn)->endICode    );
    CreateCopy ( comment,PTurn(Turn)->comment );
  }

  void  Turn::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte ( &Version    );
    f.WriteInt  ( &serNum     );
    f.WriteInt  ( &initSeqNum );
    f.WriteInt  ( &endSeqNum  );
    f.WriteTerLine ( turnID     ,false );
    f.WriteTerLine ( initResName,false );
    f.WriteTerLine ( initChainID,false );
    f.WriteTerLine ( initICode  ,false );
    f.WriteTerLine ( endResName ,false );
    f.WriteTerLine ( endChainID ,false );
    f.WriteTerLine ( endICode   ,false );
    f.CreateWrite ( comment );
  }

  void  Turn::read ( io::RFile f )  {
  byte Version;
    f.ReadByte ( &Version );
    f.ReadInt  ( &serNum     );
    f.ReadInt  ( &initSeqNum );
    f.ReadInt  ( &endSeqNum  );
    f.ReadTerLine ( turnID     ,false );
    f.ReadTerLine ( initResName,false );
    f.ReadTerLine ( initChainID,false );
    f.ReadTerLine ( initICode  ,false );
    f.ReadTerLine ( endResName ,false );
    f.ReadTerLine ( endChainID ,false );
    f.ReadTerLine ( endICode   ,false );
    f.CreateRead ( comment );
  }

  MakeStreamFunctions(Turn)


  //  ===================  LinkContainer  ========================

  PContainerClass LinkContainer::MakeContainerClass ( int ClassID )  {
    switch (ClassID)  {
      default :
      case ClassID_Template : return
                          ClassContainer::MakeContainerClass(ClassID);
      case ClassID_Link    : return new Link();
    }
  }

  MakeStreamFunctions(LinkContainer)



  //  ========================  Link  ===========================

  Link::Link() : ContainerClass()  {
    InitLink();
  }

  Link::Link ( cpstr S ) : ContainerClass()  {
    InitLink();
    ConvertPDBASCII ( S );
  }

  Link::Link ( io::RPStream Object ) : ContainerClass(Object)  {
    InitLink();
  }

  Link::~Link() {}

  void  Link::InitLink()  {
    strcpy ( atName1 ,"----" );  // name of 1st linked atom
    strcpy ( aloc1   ," "    );  // alternative location of 1st atom
    strcpy ( resName1,"---"  );  // residue name of 1st linked atom
    strcpy ( chainID1," "    );  // chain ID of 1st linked atom
    seqNum1 = 0;                 // sequence number of 1st linked atom
    strcpy ( insCode1," "    );  // insertion code of 1st linked atom
    strcpy ( atName2 ,"----" );  // name of 2nd linked atom
    strcpy ( aloc2   ," "    );  // alternative location of 2nd atom
    strcpy ( resName2,"---"  );  // residue name of 2nd linked atom
    strcpy ( chainID2," "    );  // chain ID of 2nd linked atom
    seqNum2 = 0;                 // sequence number of 2nd linked atom
    strcpy ( insCode2," "    );  // insertion code of 2nd linked atom
    s1   = 1;  // sym id of 1st atom
    i1   = 5;
    j1   = 5;
    k1   = 5;
    s2   = 1;  // sym id of 2nd atom
    i2   = 5;
    j2   = 5;
    k2   = 5;
    dist = -1.0;  // no distance
  }


  void  Link::PDBASCIIDump ( pstr S, int N )  {
  UNUSED_ARGUMENT(N);
  //  makes the ASCII PDB OBSLTE line number N
  //  from the class' data

    strcpy     ( S,"LINK" );
    PadSpaces  ( S,80 );

    strcpy_n1  ( &(S[12]),atName1 ,4 );
    strcpy_n1  ( &(S[16]),aloc1   ,1 );
    strcpy_n1  ( &(S[17]),resName1,3 );
    strcpy_n1  ( &(S[21]),chainID1,1 );
    PutIntIns  ( &(S[22]),seqNum1 ,4,insCode1 );

    strcpy_n1  ( &(S[42]),atName2 ,4 );
    strcpy_n1  ( &(S[46]),aloc2   ,1 );
    strcpy_n1  ( &(S[47]),resName2,3 );
    strcpy_n1  ( &(S[51]),chainID2,1 );
    PutIntIns  ( &(S[52]),seqNum2 ,4,insCode2 );

    PutInteger ( &(S[59]),s1,3 );
    PutInteger ( &(S[62]),i1,1 );
    PutInteger ( &(S[63]),j1,1 );
    PutInteger ( &(S[64]),k1,1 );

    PutInteger ( &(S[66]),s2,3 );
    PutInteger ( &(S[69]),i2,1 );
    PutInteger ( &(S[70]),j2,1 );
    PutInteger ( &(S[71]),k2,1 );

    if (dist>0.0)
      PutRealF ( &(S[73]),dist,5,3 );

  }


  #define  LinkTypeID  "LINK"

  void AddStructConnTags ( mmcif::PLoop Loop )  {

    Loop->AddLoopTag ( CIFTAG_ID                           );
    Loop->AddLoopTag ( CIFTAG_CONN_TYPE_ID                 );

    Loop->AddLoopTag ( CIFTAG_CONN_PTNR1_AUTH_ATOM_ID      );
    Loop->AddLoopTag ( CIFTAG_CONN_PDBX_PTNR1_AUTH_ALT_ID  );
    Loop->AddLoopTag ( CIFTAG_CONN_PTNR1_AUTH_COMP_ID      );
    Loop->AddLoopTag ( CIFTAG_CONN_PTNR1_AUTH_ASYM_ID      );
    Loop->AddLoopTag ( CIFTAG_CONN_PTNR1_AUTH_SEQ_ID       );
    Loop->AddLoopTag ( CIFTAG_CONN_PDBX_PTNR1_PDB_INS_CODE );

    Loop->AddLoopTag ( CIFTAG_CONN_PTNR2_AUTH_ATOM_ID      );
    Loop->AddLoopTag ( CIFTAG_CONN_PDBX_PTNR2_AUTH_ALT_ID  );
    Loop->AddLoopTag ( CIFTAG_CONN_PTNR2_AUTH_COMP_ID      );
    Loop->AddLoopTag ( CIFTAG_CONN_PTNR2_AUTH_ASYM_ID      );
    Loop->AddLoopTag ( CIFTAG_CONN_PTNR2_AUTH_SEQ_ID       );
    Loop->AddLoopTag ( CIFTAG_CONN_PDBX_PTNR2_PDB_INS_CODE );

    Loop->AddLoopTag ( CIFTAG_CONN_PTNR1_SYMMETRY );
    Loop->AddLoopTag ( CIFTAG_CONN_PTNR2_SYMMETRY );

    Loop->AddLoopTag ( CIFTAG_CONN_DIST );

  }


  void  Link::MakeCIF ( mmcif::PData CIF, int N )  {
  UNUSED_ARGUMENT(N);
  mmcif::PLoop Loop;
  char        S[100];
  int         RC;

    RC = CIF->AddLoop ( CIFCAT_STRUCT_CONN,Loop );
    if (RC!=mmcif::CIFRC_Ok) // the category was (re)created, provide tags
      AddStructConnTags ( Loop );

    Loop->AddString  ( "1"      );  // should be a counter
    Loop->AddString  ( pstr(LinkTypeID) );

    Loop->AddString  ( atName1  );
    Loop->AddString  ( aloc1    );
    Loop->AddString  ( resName1 );
    Loop->AddString  ( chainID1 );
    Loop->AddInteger ( seqNum1  );
    Loop->AddString  ( insCode1 );

    Loop->AddString  ( atName2  );
    Loop->AddString  ( aloc2    );
    Loop->AddString  ( resName2 );
    Loop->AddString  ( chainID2 );
    Loop->AddInteger ( seqNum2  );
    Loop->AddString  ( insCode2 );

    sprintf ( S,"%i%i%i%i",s1,i1,j1,k1 );
    Loop->AddString  ( S );
    sprintf ( S,"%i%i%i%i",s2,i2,j2,k2 );
    Loop->AddString  ( S );

    Loop->AddReal    ( dist     );

  }

  ERROR_CODE Link::ConvertPDBASCII ( cpstr S )  {

    GetString    ( atName1 ,&(S[12]),4 );
    strcpy_ncss  ( aloc1   ,&(S[16]),1 );
    strcpy_ncss  ( resName1,&(S[17]),3 );
    strcpy_ncss  ( chainID1,&(S[21]),1 );
    GetIntIns    ( seqNum1,insCode1,&(S[22]),4 );

    GetString    ( atName2 ,&(S[42]),4 );
    strcpy_ncss  ( aloc2   ,&(S[46]),1 );
    strcpy_ncss  ( resName2,&(S[47]),3 );
    strcpy_ncss  ( chainID2,&(S[51]),1 );
    GetIntIns    ( seqNum2,insCode2,&(S[52]),4 );

    GetInteger   ( s1,&(S[59]),3 );
    GetInteger   ( i1,&(S[62]),1 );
    GetInteger   ( j1,&(S[63]),1 );
    GetInteger   ( k1,&(S[64]),1 );

    GetInteger   ( s2,&(S[66]),3 );
    GetInteger   ( i2,&(S[69]),1 );
    GetInteger   ( j2,&(S[70]),1 );
    GetInteger   ( k2,&(S[71]),1 );

    if (!GetReal(dist,&(S[73]),5))  dist = -1.0;

    return Error_NoError;

  }

  ERROR_CODE Link::GetCIF ( mmcif::PData CIF, int & n )  {
  mmcif::PLoop Loop;
  pstr         F;
  char         S[100];
  int          RC,l;
  bool         Done;
  ERROR_CODE   rc;

    Loop = CIF->GetLoop ( CIFCAT_STRUCT_CONN );
    if (!Loop)  {
      n = -1;  // signal to finish processing of this structure
      return Error_EmptyCIF;
    }

    l    = Loop->GetLoopLength();
    Done = (n>=l);
    while (!Done) {
      F = Loop->GetString ( CIFTAG_CONN_TYPE_ID,n,RC );
      if ((!RC) && F)  Done = (strcmp(F,LinkTypeID)==0);
                 else  Done = false;
      if (!Done)  {
        n++;
        Done = (n>=l);
      }
    }

    if (n>=l)  {
      n = -1;  // finish processing of Turn
      return Error_EmptyCIF;
    }

    Loop->DeleteField ( CIFTAG_CONN_TYPE_ID,n );

  //  CIFGetInteger ( l,Loop,CIFTAG_ID,n );

    CIFGetString ( atName1,Loop,CIFTAG_CONN_PTNR1_AUTH_ATOM_ID,n,
                   sizeof(atName1),pstr("    ") );
    CIFGetString ( aloc1,Loop,CIFTAG_CONN_PDBX_PTNR1_AUTH_ALT_ID,n,
                   sizeof(aloc1),pstr(" ") );
    CIFGetString ( resName1,Loop,CIFTAG_CONN_PTNR1_AUTH_COMP_ID,n,
                   sizeof(resName1),pstr("   ") );
    CIFGetString ( chainID1,Loop,CIFTAG_CONN_PTNR1_AUTH_ASYM_ID,n,
                   sizeof(chainID1),pstr(" ") );
    rc = CIFGetInteger ( seqNum1,Loop,CIFTAG_CONN_PTNR1_AUTH_SEQ_ID,n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    CIFGetString ( insCode1,Loop,CIFTAG_CONN_PDBX_PTNR1_PDB_INS_CODE,
                   n,sizeof(insCode1),pstr(" ") );

    CIFGetString ( atName2,Loop,CIFTAG_CONN_PTNR2_AUTH_ATOM_ID,n,
                   sizeof(atName2),pstr("    ") );
    CIFGetString ( aloc2,Loop,CIFTAG_CONN_PDBX_PTNR2_AUTH_ALT_ID,n,
                   sizeof(aloc2),pstr(" ") );
    CIFGetString ( resName2,Loop,CIFTAG_CONN_PTNR2_AUTH_COMP_ID,n,
                   sizeof(resName2),pstr("   ") );
    CIFGetString ( chainID2,Loop,CIFTAG_CONN_PTNR2_AUTH_ASYM_ID,n,
                   sizeof(chainID2),pstr(" ") );
    rc = CIFGetInteger ( seqNum2,Loop,CIFTAG_CONN_PTNR2_AUTH_SEQ_ID,n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    CIFGetString ( insCode2,Loop,CIFTAG_CONN_PDBX_PTNR2_PDB_INS_CODE,
                   n,sizeof(insCode2),pstr(" ") );

    CIFGetString ( S,Loop,CIFTAG_CONN_PTNR1_SYMMETRY,n,
                   sizeof(S),pstr("") );
    if (S[0])  {
      l  = strlen(S)-1;
      k1 = int(S[l--]) - int('0');
      j1 = int(S[l--]) - int('0');
      i1 = int(S[l--]) - int('0');
      S[l] = char(0);
      s1 = atoi(S);
    }

    CIFGetString ( S,Loop,CIFTAG_CONN_PTNR2_SYMMETRY,n,
                   sizeof(S),pstr("") );
    if (S[0])  {
      l  = strlen(S)-1;
      k2 = int(S[l--]) - int('0');
      j2 = int(S[l--]) - int('0');
      i2 = int(S[l--]) - int('0');
      S[l] = char(0);
      s2 = atoi(S);
    }

    rc = CIFGetReal ( dist,Loop,CIFTAG_CONN_DIST,n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    n++;

    return Error_NoError;

  }

  void  Link::Copy ( PContainerClass Link )  {

    strcpy ( atName1 ,PLink(Link)->atName1  );
    strcpy ( aloc1   ,PLink(Link)->aloc1    );
    strcpy ( resName1,PLink(Link)->resName1 );
    strcpy ( chainID1,PLink(Link)->chainID1 );
    seqNum1 = PLink(Link)->seqNum1;
    strcpy ( insCode1,PLink(Link)->insCode1 );

    strcpy ( atName2 ,PLink(Link)->atName2  );
    strcpy ( aloc2   ,PLink(Link)->aloc2    );
    strcpy ( resName2,PLink(Link)->resName2 );
    strcpy ( chainID2,PLink(Link)->chainID2 );
    seqNum2 = PLink(Link)->seqNum2;
    strcpy ( insCode2,PLink(Link)->insCode2 );

    s1 = PLink(Link)->s1;
    i1 = PLink(Link)->i1;
    j1 = PLink(Link)->j1;
    k1 = PLink(Link)->k1;

    s2 = PLink(Link)->s2;
    i2 = PLink(Link)->i2;
    j2 = PLink(Link)->j2;
    k2 = PLink(Link)->k2;

    dist = PLink(Link)->dist;

  }

  void  Link::write ( io::RFile f )  {
  byte Version=2;

    f.WriteByte ( &Version    );

    f.WriteTerLine ( atName1 ,false );
    f.WriteTerLine ( aloc1   ,false );
    f.WriteTerLine ( resName1,false );
    f.WriteTerLine ( chainID1,false );
    f.WriteInt     ( &seqNum1 );
    f.WriteTerLine ( insCode1,false );

    f.WriteTerLine ( atName2 ,false );
    f.WriteTerLine ( aloc2   ,false );
    f.WriteTerLine ( resName2,false );
    f.WriteTerLine ( chainID2,false );
    f.WriteInt     ( &seqNum2 );
    f.WriteTerLine ( insCode2,false );

    f.WriteInt ( &s1 );
    f.WriteInt ( &i1 );
    f.WriteInt ( &j1 );
    f.WriteInt ( &k1 );

    f.WriteInt ( &s2 );
    f.WriteInt ( &i2 );
    f.WriteInt ( &j2 );
    f.WriteInt ( &k2 );

    f.WriteReal ( &dist );

  }

  void  Link::read ( io::RFile f )  {
  byte Version;

    f.ReadByte ( &Version    );

    f.ReadTerLine ( atName1 ,false );
    f.ReadTerLine ( aloc1   ,false );
    f.ReadTerLine ( resName1,false );
    f.ReadTerLine ( chainID1,false );
    f.ReadInt     ( &seqNum1 );
    f.ReadTerLine ( insCode1,false );

    f.ReadTerLine ( atName2 ,false );
    f.ReadTerLine ( aloc2   ,false );
    f.ReadTerLine ( resName2,false );
    f.ReadTerLine ( chainID2,false );
    f.ReadInt     ( &seqNum2 );
    f.ReadTerLine ( insCode2,false );

    f.ReadInt ( &s1 );
    f.ReadInt ( &i1 );
    f.ReadInt ( &j1 );
    f.ReadInt ( &k1 );

    f.ReadInt ( &s2 );
    f.ReadInt ( &i2 );
    f.ReadInt ( &j2 );
    f.ReadInt ( &k2 );

    if (Version>1)
      f.ReadReal ( &dist );

  }

  MakeStreamFunctions(Link)


  //  ===================  LinkRContainer  =======================

  PContainerClass LinkRContainer::MakeContainerClass ( int ClassID )  {
    switch (ClassID)  {
      default :
      case ClassID_Template : return
                           ClassContainer::MakeContainerClass(ClassID);
      case ClassID_LinkR    : return new LinkR();
    }
  }

  MakeStreamFunctions(LinkRContainer)


  //  ========================  LinkR  ===========================

  LinkR::LinkR() : ContainerClass()  {
    InitLinkR();
  }

  LinkR::LinkR ( cpstr S ) : ContainerClass()  {
    InitLinkR();
    ConvertPDBASCII ( S );
  }

  LinkR::LinkR ( io::RPStream Object ) : ContainerClass(Object)  {
    InitLinkR();
  }

  LinkR::~LinkR() {}

  void  LinkR::InitLinkR()  {
    strcpy ( linkRID ,"----" );  // link name
    strcpy ( atName1 ,"----" );  // name of 1st linked atom
    strcpy ( aloc1   ," "    );  // alternative location of 1st atom
    strcpy ( resName1,"---"  );  // residue name of 1st linked atom
    strcpy ( chainID1," "    );  // chain ID of 1st linked atom
    seqNum1 = 0;                 // sequence number of 1st linked atom
    strcpy ( insCode1," "    );  // insertion code of 1st linked atom
    strcpy ( atName2 ,"----" );  // name of 2nd linked atom
    strcpy ( aloc2   ," "    );  // alternative location of 2nd atom
    strcpy ( resName2,"---"  );  // residue name of 2nd linked atom
    strcpy ( chainID2," "    );  // chain ID of 2nd linked atom
    seqNum2 = 0;                 // sequence number of 2nd linked atom
    strcpy ( insCode2," "    );  // insertion code of 2nd linked atom
    dist    = 0.0;               // link distance
  }

  /*
  LINK             LYS A  27                     PLP A 255                PLPLYS
  LINK             MAN S   3                     MAN S   4                BETA1-4
  LINK        C6  BBEN B   1                O1  BMAF S   2                BEN-MAF
  LINK        OE2 AGLU A 320                C1  AMAF S   2                GLU-MAF
  LINK        OE2  GLU A  67        1.895   ZN   ZN  R   5                GLU-ZN
  LINK        NE2  HIS A  71        2.055   ZN   ZN  R   5                HIS-ZN
  LINK        O    ARG A  69        2.240   NA   NA  R   9                ARG-NA
  012345678901234567890123456789012345678901234567890123456789012345678901234567890
            1         2         3         4         5         6         7
  */

  void  LinkR::PDBASCIIDump ( pstr S, int N )  {
  UNUSED_ARGUMENT(N);
  //  makes the ASCII PDB OBSLTE line number N
  //  from the class' data

    strcpy     ( S,"LINKR" );
    PadSpaces  ( S,80 );

    strcpy_n1  ( &(S[12]),atName1 ,4 );
    strcpy_n1  ( &(S[16]),aloc1   ,1 );
    strcpy_n1  ( &(S[17]),resName1,3 );
    strcpy_n1  ( &(S[21]),chainID1,1 );
    PutIntIns  ( &(S[22]),seqNum1 ,4,insCode1 );

    if (dist>0.0)
      PutRealF ( &(S[32]),dist,7,3 );

    strcpy_n1  ( &(S[42]),atName2 ,4 );
    strcpy_n1  ( &(S[46]),aloc2   ,1 );
    strcpy_n1  ( &(S[47]),resName2,3 );
    strcpy_n1  ( &(S[51]),chainID2,1 );
    PutIntIns  ( &(S[52]),seqNum2 ,4,insCode2 );

    strcpy_ns  ( &(S[72]),linkRID,8 );

  }


  #define  LinkRTypeID  "LINKR"

  void AddStructConnLinkRTags ( mmcif::PLoop Loop )  {

    Loop->AddLoopTag ( CIFTAG_ID                           );
    Loop->AddLoopTag ( CIFTAG_CONN_TYPE_ID                 );

    Loop->AddLoopTag ( CIFTAG_CONN_PTNR1_AUTH_ATOM_ID      );
    Loop->AddLoopTag ( CIFTAG_CONN_PDBX_PTNR1_AUTH_ALT_ID  );
    Loop->AddLoopTag ( CIFTAG_CONN_PTNR1_AUTH_COMP_ID      );
    Loop->AddLoopTag ( CIFTAG_CONN_PTNR1_AUTH_ASYM_ID      );
    Loop->AddLoopTag ( CIFTAG_CONN_PTNR1_AUTH_SEQ_ID       );
    Loop->AddLoopTag ( CIFTAG_CONN_PDBX_PTNR1_PDB_INS_CODE );

    Loop->AddLoopTag ( CIFTAG_CONN_DIST );

    Loop->AddLoopTag ( CIFTAG_CONN_PTNR2_AUTH_ATOM_ID      );
    Loop->AddLoopTag ( CIFTAG_CONN_PDBX_PTNR2_AUTH_ALT_ID  );
    Loop->AddLoopTag ( CIFTAG_CONN_PTNR2_AUTH_COMP_ID      );
    Loop->AddLoopTag ( CIFTAG_CONN_PTNR2_AUTH_ASYM_ID      );
    Loop->AddLoopTag ( CIFTAG_CONN_PTNR2_AUTH_SEQ_ID       );
    Loop->AddLoopTag ( CIFTAG_CONN_PDBX_PTNR2_PDB_INS_CODE );

    Loop->AddLoopTag ( CIFTAG_CONN_NAME );

  }

  void  LinkR::MakeCIF ( mmcif::PData CIF, int N )  {
  UNUSED_ARGUMENT(N);
  mmcif::PLoop Loop;
  int          RC;

    RC = CIF->AddLoop ( CIFCAT_STRUCT_LINKR,Loop );
    if (RC!=mmcif::CIFRC_Ok) // the category was (re)created, provide tags
      AddStructConnLinkRTags ( Loop );

    Loop->AddString  ( "1"      );  // should be a counter
    Loop->AddString  ( pstr(LinkTypeID) );

    Loop->AddString  ( atName1  );
    Loop->AddString  ( aloc1    );
    Loop->AddString  ( resName1 );
    Loop->AddString  ( chainID1 );
    Loop->AddInteger ( seqNum1  );
    Loop->AddString  ( insCode1 );

    Loop->AddReal    ( dist     );

    Loop->AddString  ( atName2  );
    Loop->AddString  ( aloc2    );
    Loop->AddString  ( resName2 );
    Loop->AddString  ( chainID2 );
    Loop->AddInteger ( seqNum2  );
    Loop->AddString  ( insCode2 );

    Loop->AddString  ( linkRID  );

  }

  ERROR_CODE LinkR::ConvertPDBASCII ( cpstr S )  {

    GetString    ( atName1 ,&(S[12]),4 );
    strcpy_ncss  ( aloc1   ,&(S[16]),1 );
    strcpy_ncss  ( resName1,&(S[17]),3 );
    strcpy_ncss  ( chainID1,&(S[21]),1 );
    GetIntIns    ( seqNum1,insCode1,&(S[22]),4 );

    if (!GetReal(dist,&(S[32]),7)) dist = 0.0;

    GetString    ( atName2 ,&(S[42]),4 );
    strcpy_ncss  ( aloc2   ,&(S[46]),1 );
    strcpy_ncss  ( resName2,&(S[47]),3 );
    strcpy_ncss  ( chainID2,&(S[51]),1 );
    GetIntIns    ( seqNum2,insCode2,&(S[52]),4 );

    strcpy_ncss  ( linkRID,&(S[72]),8 );

    return Error_NoError ;

  }

  ERROR_CODE LinkR::GetCIF ( mmcif::PData CIF, int & n )  {
  mmcif::PLoop Loop;
  pstr         F;
  int          RC,l;
  bool         Done;
  ERROR_CODE   rc;

    Loop = CIF->GetLoop ( CIFCAT_STRUCT_CONN );
    if (!Loop)  {
      n = -1;  // signal to finish processing of this structure
      return Error_EmptyCIF;
    }

    l    = Loop->GetLoopLength();
    Done = (n>=l);
    while (!Done) {
      F = Loop->GetString ( CIFTAG_CONN_TYPE_ID,n,RC );
      if ((!RC) && F)  Done = (strcmp(F,LinkTypeID)==0);
                 else  Done = false;
      if (!Done)  {
        n++;
        Done = (n>=l);
      }
    }

    if (n>=l)  {
      n = -1;  // finish processing of Turn
      return Error_EmptyCIF;
    }

    Loop->DeleteField ( CIFTAG_CONN_TYPE_ID,n );

    //  CIFGetInteger ( l,Loop,CIFTAG_ID,n );

    CIFGetString ( atName1,Loop,CIFTAG_CONN_PTNR1_AUTH_ATOM_ID,n,
                   sizeof(atName1),pstr("    ") );
    CIFGetString ( aloc1,Loop,CIFTAG_CONN_PDBX_PTNR1_AUTH_ALT_ID,n,
                   sizeof(aloc1),pstr(" ") );
    CIFGetString ( resName1,Loop,CIFTAG_CONN_PTNR1_AUTH_COMP_ID,n,
                   sizeof(resName1),pstr("   ") );
    CIFGetString ( chainID1,Loop,CIFTAG_CONN_PTNR1_AUTH_ASYM_ID,n,
                   sizeof(chainID1),pstr(" ") );
    rc = CIFGetInteger ( seqNum1,Loop,CIFTAG_CONN_PTNR1_AUTH_SEQ_ID,n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    CIFGetString ( insCode1,Loop,CIFTAG_CONN_PDBX_PTNR1_PDB_INS_CODE,
                   n,sizeof(insCode1),pstr(" ") );

    rc = CIFGetReal ( dist,Loop,CIFTAG_CONN_DIST,n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    CIFGetString ( atName2,Loop,CIFTAG_CONN_PTNR2_AUTH_ATOM_ID,n,
                   sizeof(atName2),pstr("    ") );
    CIFGetString ( aloc2,Loop,CIFTAG_CONN_PDBX_PTNR2_AUTH_ALT_ID,n,
                   sizeof(aloc2),pstr(" ") );
    CIFGetString ( resName2,Loop,CIFTAG_CONN_PTNR2_AUTH_COMP_ID,n,
                   sizeof(resName2),pstr("   ") );
    CIFGetString ( chainID2,Loop,CIFTAG_CONN_PTNR2_AUTH_ASYM_ID,n,
                   sizeof(chainID2),pstr(" ") );
    rc = CIFGetInteger ( seqNum2,Loop,CIFTAG_CONN_PTNR2_AUTH_SEQ_ID,n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    CIFGetString ( insCode2,Loop,CIFTAG_CONN_PDBX_PTNR2_PDB_INS_CODE,
                   n,sizeof(insCode2),pstr(" ") );

    CIFGetString ( linkRID,Loop,CIFTAG_CONN_NAME,n,
                   sizeof(linkRID),pstr(" ") );

    n++;

    return Error_NoError;

  }

  void  LinkR::Copy ( PContainerClass LinkR )  {

    strcpy ( atName1 ,PLinkR(LinkR)->atName1  );
    strcpy ( aloc1   ,PLinkR(LinkR)->aloc1    );
    strcpy ( resName1,PLinkR(LinkR)->resName1 );
    strcpy ( chainID1,PLinkR(LinkR)->chainID1 );
    seqNum1 = PLinkR(LinkR)->seqNum1;
    strcpy ( insCode1,PLinkR(LinkR)->insCode1 );

    dist = PLinkR(LinkR)->dist;

    strcpy ( atName2 ,PLinkR(LinkR)->atName2  );
    strcpy ( aloc2   ,PLinkR(LinkR)->aloc2    );
    strcpy ( resName2,PLinkR(LinkR)->resName2 );
    strcpy ( chainID2,PLinkR(LinkR)->chainID2 );
    seqNum2 = PLinkR(LinkR)->seqNum2;
    strcpy ( insCode2,PLinkR(LinkR)->insCode2 );

    strcpy ( linkRID,PLinkR(LinkR)->linkRID );

  }

  void  LinkR::write ( io::RFile f )  {
  byte Version=1;

    f.WriteByte ( &Version    );

    f.WriteTerLine ( atName1 ,false );
    f.WriteTerLine ( aloc1   ,false );
    f.WriteTerLine ( resName1,false );
    f.WriteTerLine ( chainID1,false );
    f.WriteInt     ( &seqNum1 );
    f.WriteTerLine ( insCode1,false );

    f.WriteReal    ( &dist );

    f.WriteTerLine ( atName2 ,false );
    f.WriteTerLine ( aloc2   ,false );
    f.WriteTerLine ( resName2,false );
    f.WriteTerLine ( chainID2,false );
    f.WriteInt     ( &seqNum2 );
    f.WriteTerLine ( insCode2,false );

    f.WriteTerLine ( linkRID,false );

  }

  void  LinkR::read ( io::RFile f )  {
  byte Version;

    f.ReadByte ( &Version    );

    f.ReadTerLine ( atName1 ,false );
    f.ReadTerLine ( aloc1   ,false );
    f.ReadTerLine ( resName1,false );
    f.ReadTerLine ( chainID1,false );
    f.ReadInt     ( &seqNum1 );
    f.ReadTerLine ( insCode1,false );

    f.ReadReal    ( &dist );

    f.ReadTerLine ( atName2 ,false );
    f.ReadTerLine ( aloc2   ,false );
    f.ReadTerLine ( resName2,false );
    f.ReadTerLine ( chainID2,false );
    f.ReadInt     ( &seqNum2 );
    f.ReadTerLine ( insCode2,false );

    f.ReadTerLine ( linkRID,false );

  }

  MakeStreamFunctions(LinkR)


  //  ===================  CisPepContainer  ======================

  PContainerClass CisPepContainer::MakeContainerClass ( int ClassID )  {
    switch (ClassID)  {
      default :
      case ClassID_Template : return
                          ClassContainer::MakeContainerClass(ClassID);
      case ClassID_CisPep   : return new CisPep();
    }
  }

  MakeStreamFunctions(CisPepContainer)


  //  ========================  CisPep  ==========================

  CisPep::CisPep() : ContainerClass()  {
    InitCisPep();
  }

  CisPep::CisPep ( cpstr S ) : ContainerClass()  {
    InitCisPep();
    ConvertPDBASCII ( S );
  }

  CisPep::CisPep ( io::RPStream Object ) : ContainerClass(Object)  {
    InitCisPep();
  }

  CisPep::~CisPep() {}

  void CisPep::InitCisPep()  {
    serNum  = 1;                //  record serial number
    strcpy ( pep1    ,"---" );  //  residue name
    strcpy ( chainID1," "   );  //  chain identifier 1
    seqNum1 = 0;                //  residue sequence number 1
    strcpy ( icode1  ," "   );  //  insertion code 1
    strcpy ( pep2    ,"---" );  //  residue name 2
    strcpy ( chainID2," "   );  //  chain identifier 2
    seqNum2 = 0;                //  residue sequence number 2
    strcpy ( icode2  ," "   );  //  insertion code 2
    modNum  = 0;                //  model number
    measure = 0.0;              //  measure of the angle in degrees.
  }

  void  CisPep::PDBASCIIDump ( pstr S, int N )  {
  UNUSED_ARGUMENT(N);

    strcpy     ( S,"CISPEP" );
    PadSpaces  ( S,80 );

    PutInteger ( &(S[7]),serNum,3 );

    strcpy_n1  ( &(S[11]),pep1    ,3 );
    strcpy_n1  ( &(S[15]),chainID1,1 );
    PutIntIns  ( &(S[17]),seqNum1 ,4,icode1 );

    strcpy_n1  ( &(S[25]),pep2    ,3 );
    strcpy_n1  ( &(S[29]),chainID2,1 );
    PutIntIns  ( &(S[31]),seqNum2 ,4,icode1 );

    PutInteger ( &(S[43]),modNum,3 );
    PutRealF   ( &(S[53]),measure,6,2 );

  }


  ERROR_CODE CisPep::ConvertPDBASCII ( cpstr S )  {

    GetInteger   ( serNum  ,&(S[7]) ,3 );

    strcpy_ncss  ( pep1    ,&(S[11]),3 );
    strcpy_ncss  ( chainID1,&(S[15]),1 );
    GetIntIns    ( seqNum1,icode1,&(S[17]),4 );

    strcpy_ncss  ( pep2    ,&(S[25]),3 );
    strcpy_ncss  ( chainID2,&(S[29]),1 );
    GetIntIns    ( seqNum2,icode2,&(S[31]),4 );

    GetInteger   ( modNum  ,&(S[43]),3 );
    GetReal      ( measure ,&(S[53]),6 );

    return Error_NoError;

  }


  void  CisPep::Copy ( PContainerClass CisPep )  {

    serNum  = PCisPep(CisPep)->serNum;

    strcpy ( pep1    ,PCisPep(CisPep)->pep1     );
    strcpy ( chainID1,PCisPep(CisPep)->chainID1 );
    seqNum1 = PCisPep(CisPep)->seqNum1;
    strcpy ( icode1  ,PCisPep(CisPep)->icode1   );

    strcpy ( pep2    ,PCisPep(CisPep)->pep2     );
    strcpy ( chainID2,PCisPep(CisPep)->chainID2 );
    seqNum2 = PCisPep(CisPep)->seqNum2;
    strcpy ( icode2  ,PCisPep(CisPep)->icode2   );

    modNum  = PCisPep(CisPep)->modNum;
    measure = PCisPep(CisPep)->measure;

  }

  void  CisPep::write ( io::RFile f )  {
  byte Version=1;

    f.WriteByte    ( &Version   );

    f.WriteInt     ( &serNum );

    f.WriteTerLine ( pep1    ,false );
    f.WriteTerLine ( chainID1,false );
    f.WriteInt     ( &seqNum1 );
    f.WriteTerLine ( icode1  ,false );

    f.WriteTerLine ( pep2    ,false );
    f.WriteTerLine ( chainID2,false );
    f.WriteInt     ( &seqNum2 );
    f.WriteTerLine ( icode2  ,false );

    f.WriteInt     ( &modNum  );
    f.WriteReal    ( &measure );

  }

  void  CisPep::read ( io::RFile f )  {
  byte Version;

    f.ReadByte    ( &Version   );

    f.ReadInt     ( &serNum );

    f.ReadTerLine ( pep1    ,false );
    f.ReadTerLine ( chainID1,false );
    f.ReadInt     ( &seqNum1 );
    f.ReadTerLine ( icode1  ,false );

    f.ReadTerLine ( pep2    ,false );
    f.ReadTerLine ( chainID2,false );
    f.ReadInt     ( &seqNum2 );
    f.ReadTerLine ( icode2  ,false );

    f.ReadInt     ( &modNum  );
    f.ReadReal    ( &measure );

  }

  MakeStreamFunctions(CisPep)



  //  =====================   Model   =======================

  Model::Model() : ProModel() {
    InitModel();
  }

  Model::Model ( PManager MMDBM, int serialNum ) : ProModel() {
    InitModel();
    manager = MMDBM;
    serNum  = serialNum;
  }

  Model::Model ( io::RPStream Object ) : ProModel(Object)  {
    InitModel();
  }

  void  Model::InitModel()  {
    serNum       = 0;
    nChains      = 0;
    nChainsAlloc = 0;
    chain        = NULL;
    manager      = NULL;
    Exclude      = true;
  }

  Model::~Model()  {
    FreeMemory();
    if (manager)  manager->_ExcludeModel ( serNum );
  }

  void Model::FreeMemory()  {

    DeleteAllChains();
    if (chain)  delete[] chain;
    chain        = NULL;
    nChains      = 0;
    nChainsAlloc = 0;

    RemoveSecStructure();
    RemoveHetInfo     ();
    RemoveLinks       ();
    RemoveLinkRs      ();
    RemoveCisPeps     ();

  }


  void  Model::SetMMDBManager ( PManager MMDBM, int serialNum )  {
    manager = MMDBM;
    serNum  = serialNum;
  }

  void  Model::CheckInAtoms()  {
  int i;
    if (manager)
      for (i=0;i<nChains;i++)
        if (chain[i])
          chain[i]->CheckInAtoms();
  }


  int  Model::GetNumberOfAtoms ( bool countTers )  {
  // returns number of atoms in the model
  int i,na;
    na = 0;
    for (i=0;i<nChains;i++)
      if (chain[i])  na += chain[i]->GetNumberOfAtoms ( countTers );
    return na;
  }

  int  Model::GetNumberOfResidues()  {
  // returns number of residues in the model
  PChain   chn;
  int      ic,ir,nr;
    nr = 0;
    for (ic=0;ic<nChains;ic++)  {
      chn = chain[ic];
      if (chn)
        for (ir=0;ir<chn->nResidues;ir++)
          if (chn->residue[ir])  nr++;
    }
    return nr;
  }


  //  ----------------  Extracting chains  --------------------------

  int  Model::GetNumberOfChains()  {
    return nChains;
  }

  PChain Model::GetChain ( int chainNo )  {
    if ((0<=chainNo) && (chainNo<nChains))
          return chain[chainNo];
    else  return NULL;
  }


  void Model::ExpandChainArray ( int nOfChains )  {
  PPChain chain1;
  int     i;
    if (nOfChains>=nChainsAlloc)  {
      nChainsAlloc = nOfChains+10;
      chain1 = new PChain[nChainsAlloc];
      for (i=0;i<nChains;i++)
        chain1[i] = chain[i];
      for (i=nChains;i<nChainsAlloc;i++)
        chain1[i] = NULL;
      if (chain)  delete[] chain;
      chain = chain1;
    }
  }

  PChain Model::GetChainCreate ( const ChainID chID,
                                 bool enforceUniqueChainID )  {
  //   Returns pointer on chain, whose identifier is
  // given in chID. If such a chain is absent in the
  // model, it is created.
  PChain  chn;
  ChainID chainID;
  int     i,k;

    // check if such a chain is already in the model
    chn = NULL;
    if (enforceUniqueChainID)  {
      k = 0;
      for (i=0;i<nChains;i++)
        if (chain[i])  {
          // here we check only first letter as it is kept in all
          // derived names
          if (chID[0]==chain[i]->chainID[0])  {
            chn = chain[i];
            if (chn->GetNumberOfResidues()>0)  k++;
          }
        }
      if (k)          sprintf ( chainID,"%s%i",chID,k-1 );
      else if (!chn)  strcpy  ( chainID,chID ); // chain is absent
                else  return chn;  // the only empty chain
    } else  {
      if (chID[0])  {
        for (i=0;(i<nChains) && (!chn);i++)
          if (chain[i])  {
            if (!strcmp(chID,chain[i]->chainID))
              chn = chain[i]; // it is there; just return the pointer
          }
      } else  {
        for (i=0;(i<nChains) && (!chn);i++)
          if (chain[i])  {
            if (!chain[i]->chainID[0])
              chn = chain[i]; // it is there; just return the pointer
          }
      }
      if (chn)  return chn;
      strcpy ( chainID,chID );
    }

    ExpandChainArray ( nChains );

    // create new chain
    chain[nChains] = newChain();
    chain[nChains]->SetChain ( chainID );
    chain[nChains]->SetModel ( this );
    nChains++;

    return chain[nChains-1];

  }

  PChain Model::CreateChain ( const ChainID chID )  {
  //   CreateChain() creates a new chain with chain ID regardless
  // the presence of same-ID chains in the model. This function
  // was introduced only for compatibility with older CCP4
  // applications and using it in any new developments should be
  // strictly discouraged.

    ExpandChainArray ( nChains );

    // create new chain
    chain[nChains] = newChain();
    chain[nChains]->SetChain ( chID );
    chain[nChains]->SetModel ( this );
    nChains++;

    return chain[nChains-1];

  }


  void  Model::GetChainTable ( PPChain & chainTable,
                               int & NumberOfChains )  {
    chainTable     = chain;
    NumberOfChains = nChains;
  }

  bool Model::GetNewChainID ( ChainID chID, int length )  {
  int     i,k;
  bool found;

    memset ( chID,0,sizeof(ChainID) );
    chID[0] = 'A';

    do  {
      found = false;
      for (i=0;(i<nChains) && (!found);i++)
        if (chain[i])
          found = (!strcmp(chID,chain[i]->chainID));
      if (found)  {
        k = 0;
        while (k<length)
          if (!chID[k])  {
            chID[k] = 'A';
            break;
          } else if (chID[k]<'Z')  {
            chID[k]++;
            break;
          } else  {
            chID[k] = 'A';
            k++;
          }
      } else
        k = 0;
    } while (found && (k<length));

    if (found)  {
      k = strlen(chID);
      while (k<length)
        chID[k++] = 'A';
    }

    return (!found);

  }


  PChain Model::GetChain ( const ChainID chID )  {
  //   Returns pointer on chain, whose identifier is
  // given in chID. If such a chain is absent in the
  // model, returns NULL.
  int     i;
  bool isChainID;
    if (chID)  isChainID = (chID[0]!=char(0));
         else  isChainID = false;
    if (isChainID)  {
      for (i=0;i<nChains;i++)
        if (chain[i])  {
          if (!strcmp(chID,chain[i]->chainID))
            return chain[i]; // it is there; just return the pointer
        }
    } else  {
      for (i=0;i<nChains;i++)
        if (chain[i])  {
          if (!chain[i]->chainID[0])
            return chain[i]; // it is there; just return the pointer
        }
    }
    return NULL;
  }


  //  ------------------  Deleting chains  --------------------------

  int Model::DeleteChain ( int chainNo )  {
    if ((0<=chainNo) && (chainNo<nChains))  {
      if (chain[chainNo])  {
        Exclude = false;
        delete chain[chainNo];
        chain[chainNo] = NULL;
        Exclude = true;
        return 1;
      }
    }
    return 0;
  }

  int Model::DeleteChain ( const ChainID chID )  {
  int i;
    if (chID[0])  {
      for (i=0;i<nChains;i++)
        if (chain[i])  {
          if (!strcmp(chID,chain[i]->chainID))  {
            Exclude  = false;
            delete chain[i];
            chain[i] = NULL;
            Exclude  = true;
            return 1;
          }
        }
    } else  {
      for (i=0;i<nChains;i++)
        if (chain[i])  {
          if (!chain[i]->chainID[0])  {
            Exclude  = false;
            delete chain[i];
            chain[i] = NULL;
            Exclude  = true;
            return 1;
          }
        }
    }
    return 0;
  }


  int Model::DeleteAllChains()  {
  int i,k;
    Exclude = false;
    k = 0;
    for (i=0;i<nChains;i++)
      if (chain[i])  {
        delete chain[i];
        chain[i] = NULL;
        k++;
      }
    nChains = 0;
    Exclude = true;
    return k;
  }

  int Model::DeleteSolventChains()  {
  int i,k;
    Exclude = false;
    k = 0;
    for (i=0;i<nChains;i++)
      if (chain[i])  {
        if (chain[i]->isSolventChain())  {
          delete chain[i];
          chain[i] = NULL;
          k++;
        }
      }
    Exclude = true;
    return k;
  }

  void Model::TrimChainTable()  {
  int i,j;
    Exclude = false;
    j = 0;
    for (i=0;i<nChains;i++)
      if (chain[i])  {
        if (chain[i]->nResidues>0)  {
          if (j<i)  {
            chain[j] = chain[i];
            chain[i] = NULL;
          }
          j++;
        } else  {
          delete chain[i];
          chain[i] = NULL;
        }
      }
    nChains = j;
    Exclude = true;
  }


  int  Model::GetNumberOfResidues ( const ChainID chainID )  {
  PChain chn;
    chn = GetChain ( chainID );
    if (chn)  return chn->nResidues;
    return 0;
  }

  int  Model::GetNumberOfResidues ( int chainNo )  {
    if ((0<=chainNo) && (chainNo<nChains))  {
      if (chain[chainNo])
        return chain[chainNo]->nResidues;
    }
    return 0;
  }

  PResidue Model::GetResidue ( const ChainID chainID, int seqNo,
                                 const InsCode insCode )  {
  PChain chn;
    chn = GetChain ( chainID );
    if (chn)
      return chn->GetResidue ( seqNo,insCode );
    return NULL;
  }

  PResidue Model::GetResidue ( const ChainID chainID, int resNo )  {
  PChain chn;
    chn = GetChain ( chainID );
    if (chn)  {
      if ((0<=resNo) && (resNo<chn->nResidues))
        return chn->residue[resNo];
    }
    return NULL;
  }

  PResidue Model::GetResidue ( int chainNo, int seqNo,
                                 const InsCode insCode )  {
    if ((0<=chainNo) && (chainNo<nChains))  {
      if (chain[chainNo])
        return chain[chainNo]->GetResidue ( seqNo,insCode );
    }
    return NULL;
  }

  PResidue Model::GetResidue ( int chainNo, int resNo )  {
    if ((0<=chainNo) && (chainNo<nChains))  {
      if (chain[chainNo]) {
        if ((0<=resNo) && (resNo<chain[chainNo]->nResidues))
          return chain[chainNo]->residue[resNo];
      }
    }
    return NULL;
  }

  int Model::GetResidueNo ( const ChainID chainID, int seqNo,
                            const InsCode insCode )  {
  PChain chn;
    chn = GetChain ( chainID );
    if (chn)
      return chn->GetResidueNo ( seqNo,insCode );
    return -2;
  }

  int Model::GetResidueNo ( int  chainNo, int seqNo,
                             const InsCode insCode )  {
    if ((0<=chainNo) && (chainNo<nChains))  {
      if (chain[chainNo])
        return chain[chainNo]->GetResidueNo ( seqNo,insCode );
    }
    return -2;
  }


  void Model::GetResidueTable ( PPResidue & resTable,
                                 int & NumberOfResidues )  {
  // resTable has to be NULL or it will be reallocated. The application
  // is responsible for deallocating the resTable (but not of its
  // residues!). This does not apply to other GetResidueTable
  // functions.
  PPChain   chn;
  PPResidue res;
  int       i,j,k,nChns,nResidues;

    if (resTable)  {
      delete[] resTable;
      resTable = NULL;
    }

    NumberOfResidues = 0;
    GetChainTable ( chn,nChns );
    for (i=0;i<nChns;i++)
      if (chn[i])  {
        chn[i]->GetResidueTable ( res,nResidues );
        NumberOfResidues += nResidues;
      }

    if (NumberOfResidues>0)  {
      resTable = new PResidue[NumberOfResidues];
      k = 0;
      GetChainTable ( chn,nChns );
      for (i=0;i<nChns;i++)
        if (chn[i])  {
          chn[i]->GetResidueTable ( res,nResidues );
          for (j=0;j<nResidues;j++)
            if (res[j])  resTable[k++] = res[j];
        }
      NumberOfResidues = k;
    }

  }

  void Model::GetResidueTable ( const ChainID chainID,
                                 PPResidue & resTable,
                                 int & NumberOfResidues )  {
  PChain chn;
    resTable         = NULL;
    NumberOfResidues = 0;
    chn = GetChain ( chainID );
    if (chn)  {
      resTable         = chn->residue;
      NumberOfResidues = chn->nResidues;
    }
  }

  void Model::GetResidueTable ( int chainNo, PPResidue & resTable,
                                 int & NumberOfResidues )  {
    resTable         = NULL;
    NumberOfResidues = 0;
    if ((0<=chainNo) && (chainNo<nChains))  {
      if (chain[chainNo])  {
        resTable         = chain[chainNo]->residue;
        NumberOfResidues = chain[chainNo]->nResidues;
      }
    }
  }


  int Model::DeleteResidue ( const ChainID chainID, int seqNo,
                              const InsCode insCode )  {
  PChain chn;
    chn = GetChain ( chainID );
    if (chn)  return chn->DeleteResidue ( seqNo,insCode );
    return 0;
  }

  int Model::DeleteResidue ( const ChainID chainID, int resNo )  {
  PChain chn;
    chn = GetChain ( chainID );
    if (chn)  return chn->DeleteResidue ( resNo );
    return 0;
  }

  int Model::DeleteResidue ( int  chainNo, int seqNo,
                              const InsCode insCode )  {
    if ((0<=chainNo) && (chainNo<nChains))  {
      if (chain[chainNo])
        return chain[chainNo]->DeleteResidue ( seqNo,insCode );
    }
    return 0;
  }

  int Model::DeleteResidue ( int chainNo, int resNo )  {
    if ((0<=chainNo) && (chainNo<nChains))  {
      if (chain[chainNo])
        return chain[chainNo]->DeleteResidue ( resNo );
    }
    return 0;
  }

  int Model::DeleteAllResidues ( const ChainID chainID )  {
  PChain chn;
    chn = GetChain ( chainID );
    if (chn)  return chn->DeleteAllResidues();
    return 0;
  }

  int Model::DeleteAllResidues ( int chainNo )  {
    if ((0<=chainNo) && (chainNo<nChains))  {
      if (chain[chainNo])
        return chain[chainNo]->DeleteAllResidues();
    }
    return 0;
  }

  int Model::DeleteAllResidues()  {
  int i,k;
    k = 0;
    for (i=0;i<nChains;i++)
      if (chain[i])
        k += chain[i]->DeleteAllResidues();
    return k;
  }


  int Model::DeleteSolvent()  {
  int i,k;
    Exclude = false;
    k = 0;
    for (i=0;i<nChains;i++)
      if (chain[i])  {
        k += chain[i]->DeleteSolvent();
        chain[i]->TrimResidueTable();
        if (chain[i]->nResidues<=0)  {
          delete chain[i];
          chain[i] = NULL;
        }
      }
    Exclude = true;
    return k;
  }


  int Model::AddResidue ( const ChainID chainID, PResidue res )  {
  PChain chn;
    chn = GetChain ( chainID );
    if (chn)  return chn->AddResidue ( res );
    return 0;
  }

  int Model::AddResidue ( int chainNo, PResidue res )  {
    if ((0<=chainNo) && (chainNo<nChains))  {
      if (chain[chainNo])
        return chain[chainNo]->AddResidue ( res );
    }
    return 0;
  }


  int  Model::_ExcludeChain ( const ChainID chainID )  {
  //   _ExcludeChain(..) excludes (but does not dispose!) a chain
  // from the model. Returns 1 if the model gets empty and 0 otherwise.
  int  i,k;

    if (!Exclude)  return 0;

    // find the chain
    k = -1;
    for (i=0;(i<nChains) && (k<0);i++)
      if (!strcmp(chainID,chain[i]->chainID))
        k = i;

    if (k>=0)  {
      for (i=k+1;i<nChains;i++)
        chain[i-1] = chain[i];
      nChains--;
      chain[nChains] = NULL;
    }

    if (nChains<=0)  return 1;
               else  return 0;

  }


  //  --------------------  Sort chains  ----------------------------

  DefineClass(QSortChains)

  class QSortChains : public QuickSort  {
    public :
      QSortChains() : QuickSort() { sKey = 0; }
      int  Compare ( int i, int j );
      void Swap    ( int i, int j );
      void Sort    ( PPChain chain, int nChains, int sortKey );
    private :
      int sKey;
  };

  int QSortChains::Compare ( int i, int j )  {
  int diff;

    diff = strcmp ( (PPChain(data))[i]->GetChainID(),
                    (PPChain(data))[j]->GetChainID() );
    if (diff>0)  diff =  1;
    if (diff<0)  diff = -1;

    if (sKey==SORT_CHAIN_ChainID_Desc) return -diff;

    return diff;

  }

  void QSortChains::Swap ( int i, int j )  {
  PChain chn;
    chn = ((PPChain)data)[i];
    ((PPChain)data)[i] = ((PPChain)data)[j];
    ((PPChain)data)[j] = chn;
  }

  void QSortChains::Sort ( PPChain chain, int nChains, int sortKey )  {
    sKey = sortKey;
    QuickSort::Sort ( &(chain[0]),nChains );
  }

  void Model::SortChains ( int sortKey )  {
  QSortChains SC;
    TrimChainTable();
    SC.Sort ( chain,nChains,sortKey );
  }


  // --------------------  Extracting atoms  -----------------------


  int  Model::GetNumberOfAtoms ( const ChainID chainID, int seqNo,
                                  const InsCode insCode )  {
  PChain   chn;
  PResidue res;
    chn = GetChain ( chainID );
    if (chn)  {
      res = chn->GetResidue ( seqNo,insCode );
      if (res)  return res->nAtoms;
    }
    return 0;
  }

  int  Model::GetNumberOfAtoms ( int chainNo, int seqNo,
                                  const InsCode insCode )  {
  PChain   chn;
  PResidue res;
    chn = GetChain ( chainNo );
    if (chn)  {
      res = chn->GetResidue ( seqNo,insCode );
      if (res)  return res->nAtoms;
    }
    return 0;
  }

  int  Model::GetNumberOfAtoms ( const ChainID chainID, int resNo )  {
  PChain   chn;
  PResidue res;
    chn = GetChain ( chainID );
    if (chn)  {
      if ((0<=resNo) && (resNo<chn->nResidues))  {
        res = chn->residue[resNo];
        if (res)  return res->nAtoms;
      }
    }
    return 0;
  }

  int  Model::GetNumberOfAtoms ( int chainNo, int resNo )  {
  PChain   chn;
  PResidue res;
    if ((0<=chainNo) && (chainNo<nChains))  {
      chn = chain[chainNo];
      if (chn)  {
        if ((0<=resNo) && (resNo<chn->nResidues))  {
          res = chn->residue[resNo];
          if (res)  return res->nAtoms;
        }
      }
    }
    return 0;
  }

  PAtom  Model::GetAtom ( const ChainID  chID,
                            int            seqNo,
                            const InsCode  insCode,
                            const AtomName aname,
                            const Element  elmnt,
                            const AltLoc   aloc
                          )  {
  PChain   chn;
  PResidue res;
    chn = GetChain ( chID );
    if (chn)  {
      res = chn->GetResidue ( seqNo,insCode );
      if (res)
        return res->GetAtom ( aname,elmnt,aloc );
    }
    return NULL;
  }

  PAtom Model::GetAtom ( const ChainID chID,    int seqNo,
                           const InsCode insCode, int   atomNo )  {
  PChain   chn;
  PResidue res;
    chn = GetChain ( chID );
    if (chn)  {
      res = chn->GetResidue ( seqNo,insCode );
      if (res)  {
        if ((0<=atomNo) && (atomNo<res->nAtoms))
          return res->atom[atomNo];
      }
    }
    return NULL;
  }

  PAtom Model::GetAtom ( const ChainID  chID,
                           int            resNo,
                           const AtomName aname,
                           const Element  elmnt,
                           const AltLoc   aloc )  {
  PChain   chn;
  PResidue res;
    chn = GetChain ( chID );
    if (chn)  {
      if ((0<=resNo) && (resNo<chn->nResidues))  {
        res = chn->residue[resNo];
        if (res)
          return res->GetAtom ( aname,elmnt,aloc );
      }
    }
    return NULL;
  }

  PAtom Model::GetAtom ( const ChainID chID, int resNo, int atomNo )  {
  PChain   chn;
  PResidue res;
    chn = GetChain ( chID );
    if (chn)  {
      if ((0<=resNo) && (resNo<chn->nResidues))  {
        res = chn->residue[resNo];
        if (res)  {
          if ((0<=atomNo) && (atomNo<res->nAtoms))
            return res->atom[atomNo];
        }
      }
    }
    return NULL;
  }

  PAtom Model::GetAtom ( int chNo, int seqNo,
                           const InsCode  insCode,
                           const AtomName aname,
                           const Element  elmnt,
                           const AltLoc   aloc )  {
  PResidue res;
    if ((0<=chNo) && (chNo<nChains))  {
      if (chain[chNo])  {
        res = chain[chNo]->GetResidue ( seqNo,insCode );
        if (res)
          return res->GetAtom ( aname,elmnt,aloc );
      }
    }
    return NULL;
  }

  PAtom Model::GetAtom ( int chNo, int seqNo, const InsCode insCode,
                           int atomNo )  {
  PResidue res;
    if ((0<=chNo) && (chNo<nChains))  {
      if (chain[chNo])  {
        res = chain[chNo]->GetResidue ( seqNo,insCode );
        if (res)  {
          if ((0<=atomNo) && (atomNo<res->nAtoms))
            return res->atom[atomNo];
        }
      }
    }
    return NULL;
  }

  PAtom Model::GetAtom ( int chNo, int resNo,
                           const AtomName aname,
                           const Element  elmnt,
                           const AltLoc   aloc )  {
  PResidue res;
    if ((0<=chNo) && (chNo<nChains))  {
      if (chain[chNo])  {
        if ((0<=resNo) && (resNo<chain[chNo]->nResidues))  {
          res = chain[chNo]->residue[resNo];
          if (res)
            return res->GetAtom ( aname,elmnt,aloc );
        }
      }
    }
    return NULL;
  }

  PAtom Model::GetAtom ( int chNo, int resNo, int atomNo )  {
  PResidue res;
    if ((0<=chNo) && (chNo<nChains))  {
      if (chain[chNo])  {
        if ((0<=resNo) && (resNo<chain[chNo]->nResidues))  {
          res = chain[chNo]->residue[resNo];
          if (res)  {
            if ((0<=atomNo) && (atomNo<res->nAtoms))
              return res->atom[atomNo];
          }
        }
      }
    }
    return NULL;
  }


  void Model::GetAtomTable ( const ChainID chainID, int seqNo,
                              const InsCode insCode,
                              PPAtom & atomTable,
                              int & NumberOfAtoms )  {
  PResidue res;
    atomTable     = NULL;
    NumberOfAtoms = 0;
    res = GetResidue ( chainID,seqNo,insCode );
    if (res)  {
      atomTable     = res->atom;
      NumberOfAtoms = res->nAtoms;
    }
  }

  void Model::GetAtomTable ( int chainNo, int seqNo,
                              const InsCode insCode,
                              PPAtom & atomTable,
                              int & NumberOfAtoms )  {
  PResidue res;
    atomTable     = NULL;
    NumberOfAtoms = 0;
    res = GetResidue ( chainNo,seqNo,insCode );
    if (res)  {
      atomTable     = res->atom;
      NumberOfAtoms = res->nAtoms;
    }
  }

  void Model::GetAtomTable ( const ChainID chainID, int resNo,
                              PPAtom & atomTable,
                              int & NumberOfAtoms )  {
  PResidue res;
    atomTable     = NULL;
    NumberOfAtoms = 0;
    res = GetResidue ( chainID,resNo );
    if (res)  {
      atomTable     = res->atom;
      NumberOfAtoms = res->nAtoms;
    }
  }

  void Model::GetAtomTable ( int chainNo, int resNo,
                              PPAtom & atomTable,
                              int & NumberOfAtoms )  {
  PResidue res;
    atomTable     = NULL;
    NumberOfAtoms = 0;
    res = GetResidue ( chainNo,resNo );
    if (res)  {
      atomTable     = res->atom;
      NumberOfAtoms = res->nAtoms;
    }
  }


  void Model::GetAtomTable1 ( const ChainID chainID, int seqNo,
                               const InsCode insCode,
                               PPAtom & atomTable,
                               int & NumberOfAtoms )  {
  PResidue res;
    res = GetResidue ( chainID,seqNo,insCode );
    if (res)
      res->GetAtomTable1 ( atomTable,NumberOfAtoms );
    else  {
      if (atomTable)  delete[] atomTable;
      atomTable     = NULL;
      NumberOfAtoms = 0;
    }
  }

  void Model::GetAtomTable1 ( int chainNo, int seqNo,
                               const InsCode insCode,
                               PPAtom & atomTable,
                               int & NumberOfAtoms )  {
  PResidue res;
    res = GetResidue ( chainNo,seqNo,insCode );
    if (res)
      res->GetAtomTable1 ( atomTable,NumberOfAtoms );
    else  {
      if (atomTable)  delete[] atomTable;
      atomTable     = NULL;
      NumberOfAtoms = 0;
    }
  }

  void Model::GetAtomTable1 ( const ChainID chainID, int resNo,
                               PPAtom & atomTable,
                               int & NumberOfAtoms )  {
  PResidue res;
    res = GetResidue ( chainID,resNo );
    if (res)
      res->GetAtomTable1 ( atomTable,NumberOfAtoms );
    else  {
      if (atomTable)  delete[] atomTable;
      atomTable     = NULL;
      NumberOfAtoms = 0;
    }
  }

  void Model::GetAtomTable1 ( int chainNo, int resNo,
                               PPAtom & atomTable,
                               int & NumberOfAtoms )  {
  PResidue res;
    res = GetResidue ( chainNo,resNo );
    if (res)
      res->GetAtomTable1 ( atomTable,NumberOfAtoms );
    else  {
      if (atomTable)  delete[] atomTable;
      atomTable     = NULL;
      NumberOfAtoms = 0;
    }
  }



  int  Model::DeleteAtom ( const ChainID  chID,
                            int            seqNo,
                            const InsCode  insCode,
                            const AtomName aname,
                            const Element  elmnt,
                            const AltLoc   aloc
                          )  {
  PChain chn;
    chn = GetChain ( chID );
    if (chn)
      return  chn->DeleteAtom ( seqNo,insCode,aname,elmnt,aloc );
    return 0;
  }

  int  Model::DeleteAtom ( const ChainID chID,    int seqNo,
                            const InsCode insCode, int   atomNo )  {
  PChain chn;
    chn = GetChain ( chID );
    if (chn)  return chn->DeleteAtom ( seqNo,insCode,atomNo );
    return 0;
  }

  int  Model::DeleteAtom ( const ChainID  chID,
                            int            resNo,
                            const AtomName aname,
                            const Element  elmnt,
                            const AltLoc   aloc )  {
  PChain chn;
    chn = GetChain ( chID );
    if (chn)  return chn->DeleteAtom ( resNo,aname,elmnt,aloc );
    return 0;
  }

  int  Model::DeleteAtom ( const ChainID chID, int resNo, int atomNo ) {
  PChain chn;
    chn = GetChain ( chID );
    if (chn)  return chn->DeleteAtom ( resNo,atomNo );
    return 0;
  }

  int  Model::DeleteAtom ( int chNo, int seqNo,
                            const InsCode  insCode,
                            const AtomName aname,
                            const Element  elmnt,
                            const AltLoc   aloc )  {
    if ((0<=chNo) && (chNo<nChains))  {
      if (chain[chNo])
        return chain[chNo]->DeleteAtom ( seqNo,insCode,aname,
                                         elmnt,aloc );
    }
    return 0;
  }

  int Model::DeleteAtom ( int chNo, int seqNo, const InsCode insCode,
                           int atomNo )  {
    if ((0<=chNo) && (chNo<nChains))  {
      if (chain[chNo])
        return  chain[chNo]->DeleteAtom ( seqNo,insCode,atomNo );
    }
    return 0;
  }

  int Model::DeleteAtom ( int chNo, int resNo,
                           const AtomName aname,
                           const Element  elmnt,
                           const AltLoc   aloc )  {
    if ((0<=chNo) && (chNo<nChains))  {
      if (chain[chNo])
        return chain[chNo]->DeleteAtom ( resNo,aname,elmnt,aloc );
    }
    return 0;
  }

  int Model::DeleteAtom ( int chNo, int resNo, int atomNo )  {
    if ((0<=chNo) && (chNo<nChains))  {
      if (chain[chNo])
        return chain[chNo]->DeleteAtom ( resNo,atomNo );
    }
    return 0;
  }

  int Model::DeleteAllAtoms ( const ChainID chID, int seqNo,
                               const InsCode insCode )  {
  PChain chn;
    chn = GetChain ( chID );
    if (chn)  return chn->DeleteAllAtoms ( seqNo,insCode );
    return 0;
  }

  int Model::DeleteAllAtoms ( const ChainID chID, int resNo )  {
  PChain chn;
    chn = GetChain ( chID );
    if (chn)  return chn->DeleteAllAtoms ( resNo );
    return 0;
  }

  int Model::DeleteAllAtoms ( const ChainID chID )  {
  PChain chn;
    chn = GetChain ( chID );
    if (chn)  return chn->DeleteAllAtoms();
    return 0;
  }

  int Model::DeleteAllAtoms ( int chNo, int seqNo,
                               const InsCode insCode )  {
    if ((0<=chNo) && (chNo<nChains))  {
      if (chain[chNo])
        return chain[chNo]->DeleteAllAtoms ( seqNo,insCode );
    }
    return 0;
  }

  int Model::DeleteAllAtoms ( int chNo, int resNo )  {
    if ((0<=chNo) && (chNo<nChains))  {
      if (chain[chNo])
        return chain[chNo]->DeleteAllAtoms ( resNo );
    }
    return 0;
  }

  int Model::DeleteAllAtoms ( int chNo )  {
    if ((0<=chNo) && (chNo<nChains))  {
      if (chain[chNo])
        return chain[chNo]->DeleteAllAtoms();
    }
    return 0;
  }

  int Model::DeleteAllAtoms()  {
  int i,k;
    k = 0;
    for (i=0;i<nChains;i++)
      if (chain[i])  k += chain[i]->DeleteAllAtoms();
    return k;
  }

  int Model::DeleteAltLocs()  {
  //  This function leaves only alternative location with maximal
  // occupancy, if those are equal or unspecified, the one with
  // "least" alternative location indicator.
  //  The function returns the number of deleted. All tables remain
  // untrimmed, so that explicit trimming or calling FinishStructEdit()
  // is required.
  int i,n;

    n = 0;
    for (i=0;i<nChains;i++)
      if (chain[i])  n += chain[i]->DeleteAltLocs();

    return n;

  }


  int Model::AddAtom ( const ChainID chID, int seqNo,
                        const InsCode insCode,
                        PAtom atom )  {
  PChain chn;
    chn = GetChain ( chID );
    if (chn)  return chn->AddAtom ( seqNo,insCode,atom );
    return 0;
  }

  int Model::AddAtom ( const ChainID chID, int resNo, PAtom  atom )  {
  PChain chn;
    chn = GetChain ( chID );
    if (chn)  return chn->AddAtom ( resNo,atom );
    return 0;
  }

  int Model::AddAtom ( int chNo, int seqNo, const InsCode insCode,
                        PAtom atom )  {
    if ((0<=chNo) && (chNo<nChains))  {
      if (chain[chNo])
        return chain[chNo]->AddAtom ( seqNo,insCode,atom );
    }
    return 0;
  }

  int Model::AddAtom ( int chNo, int resNo, PAtom atom )  {
    if ((0<=chNo) && (chNo<nChains))  {
      if (chain[chNo])
        return chain[chNo]->AddAtom ( resNo,atom );
    }
    return 0;
  }



  void  Model::GetAtomStatistics ( RAtomStat AS )  {
    AS.Init();
    CalAtomStatistics ( AS );
    AS.Finish();
  }

  void  Model::CalAtomStatistics ( RAtomStat AS )  {
  int i;
    for (i=0;i<nChains;i++)
      if (chain[i])  chain[i]->CalAtomStatistics ( AS );
  }



  ERROR_CODE Model::ConvertPDBString ( pstr PDBString ) {
  //   Interprets PDB records DBREF, SEQADV, SEQRES, MODRES.
  //   Returns zero if the line was converted, otherwise returns a
  // non-negative value of Error_XXXX.
  //   PDBString must be not shorter than 81 characters.
  ChainID    chainID;
  PChain     chn;
  PHelix     helix;
  PTurn      turn;
  PLink      link;
  PLinkR     linkR;
  PCisPep    cispep;
  ERROR_CODE RC;

    //  pad input line with spaces, if necessary
    PadSpaces ( PDBString,80 );

    chainID[0] = char(0);
    chainID[1] = char(0);

    if (!strncmp(PDBString,"DBREF ",6))  {

      if (PDBString[12]!=' ')  chainID[0] = PDBString[12];
      chn = GetChainCreate ( chainID,false );
      return chn->ConvertDBREF ( PDBString );

    } else if (!strncmp(PDBString,"SEQADV",6))  {

      if (PDBString[16]!=' ')  chainID[0] = PDBString[16];
      chn = GetChainCreate ( chainID,false );
      return chn->ConvertSEQADV ( PDBString );

    } else if (!strncmp(PDBString,"SEQRES",6))  {

      if (PDBString[11]!=' ')  chainID[0] = PDBString[11];
      chn = GetChainCreate ( chainID,false );
      return chn->ConvertSEQRES ( PDBString );

    } else if (!strncmp(PDBString,"MODRES",6))  {

      if (PDBString[16]!=' ')  chainID[0] = PDBString[16];
      chn = GetChainCreate ( chainID,false );
      return chn->ConvertMODRES ( PDBString );

    } else if (!strncmp(PDBString,"HET   ",6))  {

      if (PDBString[12]!=' ')  chainID[0] = PDBString[12];
      chn = GetChainCreate ( chainID,false );
      return chn->ConvertHET ( PDBString );

    } else if (!strncmp(PDBString,"HETNAM",6))  {

      hetCompounds.ConvertHETNAM ( PDBString );
      return Error_NoError;

    } else if (!strncmp(PDBString,"HETSYN",6))  {

      hetCompounds.ConvertHETSYN ( PDBString );
      return Error_NoError;

    } else if (!strncmp(PDBString,"FORMUL",6))  {

      hetCompounds.ConvertFORMUL ( PDBString );
      return Error_NoError;

    } else if (!strncmp(PDBString,"HELIX ",6))  {

      helix = new Helix();
      RC    = helix->ConvertPDBASCII(PDBString);
      if (RC==0)  helices.AddData ( helix );
            else  delete helix;
      return RC;

    } else if (!strncmp(PDBString,"SHEET ",6))  {

      return sheets.ConvertPDBASCII ( PDBString );

    } else if (!strncmp(PDBString,"TURN  ",6))  {

      turn = new Turn();
      RC   = turn->ConvertPDBASCII(PDBString);
      if (RC==0)  turns.AddData ( turn );
            else  delete turn;
      return RC;

    } else if (!strncmp(PDBString,"LINK  ",6))  {

      link = new Link();
      RC   = link->ConvertPDBASCII(PDBString);
      if (RC==0)  links.AddData ( link );
            else  delete link;
      return RC;


    } else if (!strncmp(PDBString,"LINKR ",6))  {

      linkR = new LinkR();
      RC   = linkR->ConvertPDBASCII(PDBString);
      if (RC==0)  linkRs.AddData ( linkR );
            else  delete linkR;
      return RC;

    } else if (!strncmp(PDBString,"CISPEP",6))  {

      cispep = new CisPep();
      RC   = cispep->ConvertPDBASCII(PDBString);
      if (RC==0)  cisPeps.AddData ( cispep );
            else  delete cispep;
      return RC;

    } else
      return Error_WrongSection;

  }


  void  Model::PDBASCIIDumpPS ( io::RFile f )  {
  int i;

    for (i=0;i<nChains;i++)
      if (chain[i])
        chain[i]->DBRef.PDBASCIIDump ( f );

    for (i=0;i<nChains;i++)
      if (chain[i])
        chain[i]->seqAdv.PDBASCIIDump ( f );

    for (i=0;i<nChains;i++)
      if (chain[i])
        chain[i]->seqRes.PDBASCIIDump ( f );

    for (i=0;i<nChains;i++)
      if (chain[i])
        chain[i]->modRes.PDBASCIIDump ( f );

    for (i=0;i<nChains;i++)
      if (chain[i])
        chain[i]->Het.PDBASCIIDump ( f );

    hetCompounds.PDBASCIIDump ( f );
    helices     .PDBASCIIDump ( f );
    sheets      .PDBASCIIDump ( f );
    turns       .PDBASCIIDump ( f );
    links       .PDBASCIIDump ( f );
    linkRs      .PDBASCIIDump ( f );

  }

  void  Model::PDBASCIIDumpCP ( io::RFile f )  {
    cisPeps.PDBASCIIDump ( f );
  }

  void  Model::PDBASCIIDump ( io::RFile f )  {
  char  S[100];
  int   i;
  bool  singleModel = true;

    if (manager)
      singleModel = (manager->nModels<=1);

    if (!singleModel)  {
      strcpy      ( S,"MODEL " );
      PadSpaces   ( S,80 );
      PutInteger  ( &(S[10]),serNum,4 );
      f.WriteLine ( S );
    }

    for (i=0;i<nChains;i++)
      if (chain[i])
        chain[i]->PDBASCIIAtomDump ( f );

    if (!singleModel)  {
      strcpy      ( S,"ENDMDL" );
      PadSpaces   ( S,80 );
      f.WriteLine ( S );
    }

  }


  void  Model::MakeAtomCIF ( mmcif::PData CIF )  {
  int  i;
    for (i=0;i<nChains;i++)
      if (chain[i])
        chain[i]->MakeAtomCIF ( CIF );
  }


  void  Model::MakePSCIF ( mmcif::PData CIF )  {
  int  i;

    for (i=0;i<nChains;i++)
      if (chain[i])
        chain[i]->DBRef.MakeCIF ( CIF );

    for (i=0;i<nChains;i++)
      if (chain[i])
        chain[i]->seqAdv.MakeCIF ( CIF );

    for (i=0;i<nChains;i++)
      if (chain[i])
        chain[i]->seqRes.MakeCIF ( CIF );

    for (i=0;i<nChains;i++)
      if (chain[i])
        chain[i]->modRes.MakeCIF ( CIF );

    for (i=0;i<nChains;i++)
      if (chain[i])
        chain[i]->Het.MakeCIF ( CIF );

    hetCompounds.MakeCIF ( CIF );
    helices     .MakeCIF ( CIF );
    sheets      .MakeCIF ( CIF );
    turns       .MakeCIF ( CIF );
    links       .MakeCIF ( CIF );
    linkRs      .MakeCIF ( CIF );

  }

  ERROR_CODE Model::GetCIFPSClass ( mmcif::PData CIF, int ClassID )  {
  ChainContainer  PSClass;
  PChainContainer Dest;
  ERROR_CODE      RC;
  cpstr           chainID;
  PChain          chn;
    PSClass.SetChain ( NULL );
    RC = PSClass.GetCIF ( CIF,ClassID );
    if (RC!=Error_NoError)  return RC;
    chainID = PSClass.Get1stChainID();
    while (chainID)  {
      chn = GetChainCreate ( chainID,false );
      switch (ClassID)  {
        case ClassID_DBReference : Dest = &(chn->DBRef);   break;
        case ClassID_SeqAdv      : Dest = &(chn->seqAdv);  break;
        case ClassID_ModRes      : Dest = &(chn->modRes);  break;
        case ClassID_Het         : Dest = &(chn->Het);     break;
        default                  : Dest = NULL;
      }
      if (Dest)  {
        PSClass.MoveByChainID ( chainID,Dest );
        Dest->SetChain ( chn );
      } else
        printf ( " **** PROGRAM ERROR: wrong call to"
                 " Model::GetCIFPSClass(..)\n" );
      chainID = PSClass.Get1stChainID();
    }
    return Error_NoError;
  }

  ERROR_CODE Model::GetCIF ( mmcif::PData CIF ) {
  SeqRes     seqRes;
  ERROR_CODE RC;
  PChain     chn;

    RC = GetCIFPSClass ( CIF,ClassID_DBReference );
    if (RC!=Error_NoError)  return RC;

    RC = GetCIFPSClass ( CIF,ClassID_SeqAdv );
    if (RC!=Error_NoError)  return RC;

    RC = seqRes.GetCIF ( CIF );
    while (RC==Error_NoError)  {
      chn = GetChainCreate ( seqRes.chainID,false );
      chn->seqRes.Copy ( &seqRes );
      RC  = seqRes.GetCIF ( CIF );
    }

    RC = GetCIFPSClass ( CIF,ClassID_ModRes );
    if (RC!=Error_NoError)  return RC;

    RC = GetCIFPSClass ( CIF,ClassID_Het );
    if (RC!=Error_NoError)  return RC;

    hetCompounds.GetCIF ( CIF );
    helices     .GetCIF ( CIF,ClassID_Helix );
    sheets      .GetCIF ( CIF );
    turns       .GetCIF ( CIF,ClassID_Turn  );
    links       .GetCIF ( CIF,ClassID_Link  );
    linkRs      .GetCIF ( CIF,ClassID_LinkR );

    return RC;

  }

  cpstr  Model::GetEntryID()  {
    if (manager)  return manager->title.idCode;
            else  return pstr("");
  }

  void  Model::SetEntryID ( const IDCode idCode )  {
    if (manager)
      manager->SetEntryID ( idCode );
  }

  int   Model::GetNumberOfAllAtoms()  {
    if (manager)  return manager->nAtoms;
            else  return 0;
  }

  int   Model::GetSerNum()  {
    return serNum;
  }

  PAtom * Model::GetAllAtoms()  {
    if (manager)  return manager->atom;
            else  return NULL;
  }


  cpstr  Model::GetModelID ( pstr modelID )  {
    modelID[0] = char(0);
    sprintf ( modelID,"/%i",serNum );
    return modelID;
  }

  int   Model::GetNumberOfModels()  {
    if (manager)  return manager->nModels;
            else  return 0;
  }


  void  Model::Copy ( PModel model )  {
  //  modify both Model::_copy and Model::Copy methods simultaneously!
  int i;

    FreeMemory();

    if (model)  {

      serNum       = model->serNum;
      nChains      = model->nChains;
      nChainsAlloc = nChains;
      if (nChains>0)  {
        chain = new PChain[nChainsAlloc];
        for (i=0;i<nChains;i++)  {
          if (model->chain[i])  {
            chain[i] = newChain();
            chain[i]->SetModel ( this );
            chain[i]->Copy ( model->chain[i] );
          } else
            chain[i] = NULL;
        }
      }

      hetCompounds.Copy ( &(model->hetCompounds) );
      helices     .Copy ( &(model->helices)      );
      sheets      .Copy ( &(model->sheets)       );
      turns       .Copy ( &(model->turns)        );
      links       .Copy ( &(model->links)        );
      linkRs      .Copy ( &(model->linkRs)       );
      cisPeps     .Copy ( &(model->cisPeps)      );

    }

  }

  void  Model::CopyHets ( PModel model )  {
    if (model)  hetCompounds.Copy ( &(model->hetCompounds) );
  }

  void  Model::CopySecStructure ( PModel model )  {
    if (model)  {
      helices.Copy ( &(model->helices) );
      sheets .Copy ( &(model->sheets)  );
      turns  .Copy ( &(model->turns)   );
    }
  }

  void  Model::CopyLinks ( PModel model )  {
    if (model)links.Copy ( &(model->links) );
  }

  void  Model::CopyLinkRs ( PModel model )  {
    if (model)  linkRs.Copy ( &(model->linkRs) );
  }

  void  Model::CopyCisPeps ( PModel model )  {
    if (model)  cisPeps.Copy ( &(model->cisPeps) );
  }

  void  Model::_copy ( PModel model )  {
  //  modify both Model::_copy and Model::Copy methods simultaneously!
  int i;

    FreeMemory();

    if (model)  {

      serNum       = model->serNum;
      nChains      = model->nChains;
      nChainsAlloc = nChains;
      if (nChains>0)  {
        chain = new PChain[nChainsAlloc];
        for (i=0;i<nChains;i++)  {
          if (model->chain[i])  {
            chain[i] = newChain();
            chain[i]->SetModel ( this );
            chain[i]->_copy ( model->chain[i] );
          } else
            chain[i] = NULL;
        }
      }

      hetCompounds.Copy ( &(model->hetCompounds) );
      helices     .Copy ( &(model->helices)      );
      sheets      .Copy ( &(model->sheets)       );
      turns       .Copy ( &(model->turns)        );
      links       .Copy ( &(model->links)        );
      linkRs      .Copy ( &(model->linkRs)       );
      cisPeps     .Copy ( &(model->cisPeps)      );

    }

  }


  void  Model::_copy ( PModel model, PPAtom atom, int & atom_index ) {
  //  modify both Model::_copy and Model::Copy methods simultaneously!
  //
  //  _copy(PModel,PPAtom,int&) does copy atoms into array 'atom'
  // starting from position atom_index. 'atom' should be able to
  // accept all new atoms - no checks on the length of 'atom'
  // is being made. This function should not be used in applications.
  int i;

    FreeMemory();

    if (model)  {

      serNum       = model->serNum;
      nChains      = model->nChains;
      nChainsAlloc = nChains;
      if (nChains>0)  {
        chain = new PChain[nChainsAlloc];
        for (i=0;i<nChains;i++)  {
          if (model->chain[i])  {
            chain[i] = newChain();
            chain[i]->SetModel ( this );
            chain[i]->_copy ( model->chain[i],atom,atom_index );
          } else
            chain[i] = NULL;
        }
      }

      hetCompounds.Copy ( &(model->hetCompounds) );
      helices     .Copy ( &(model->helices)      );
      sheets      .Copy ( &(model->sheets)       );
      turns       .Copy ( &(model->turns)        );
      links       .Copy ( &(model->links)        );
      linkRs      .Copy ( &(model->linkRs)       );

    }

  }


  int  Model::AddChain ( PChain chn )  {
  //  modify both Model::Copy methods simultaneously!
  //
  //  Copy(PModel,PPAtom,int&) copies atoms into array 'atom'
  // starting from position atom_index. 'atom' should be able to
  // accept all new atoms - no checks on the length of 'atom'
  // is being made. This function should not be used in applications.
  PModel  model1;
  int     i;

    for (i=0;i<nChains;i++)
      if (chain[i]==chn)  return -i;  // this chain is already there

    if (chn)  {

      // get space for new chain
      ExpandChainArray ( nChains );

      if (chn->GetCoordHierarchy())  {
        // The chain is associated with a coordinate hierarchy. It should
        // remain there, therefore we physically copy all its residues
        // and atoms.
        chain[nChains] = newChain();
        chain[nChains]->SetModel ( this );
        if (manager)  {
          // get space for new atoms
          manager->AddAtomArray ( chn->GetNumberOfAtoms(true) );
          chain[nChains]->_copy ( chn,manager->atom,manager->nAtoms );
        } else  {
          for (i=0;i<chn->nResidues;i++)
            chain[nChains]->AddResidue ( chn->residue[i] );
        }
      } else  {
        // The chain is not associated with a coordinate hierarchy. Such
        // unregistered objects are simply taken over, i.e. moved into
        // the new destination (model).
        chain[nChains] = chn;
        // remove chain from its model:
        model1 = chn->GetModel();
        if (model1)
          for (i=0;i<model1->nChains;i++)
            if (model1->chain[i]==chn)  {
              model1->chain[i] = NULL;
              break;
            }
        chain[nChains]->SetModel ( this );
        if (manager)
          chain[nChains]->CheckInAtoms();
      }

      nChains++;

    }

    return nChains;

  }


  void  Model::MoveChain ( PChain & m_chain, PPAtom m_atom,
                            PPAtom  atom, int & atom_index,
                            int  chain_ext )  {
  //   MoveChain(..) adds chain m_chain on the top Chain array.
  // The pointer on chain is then set to NULL (m_chain=NULL).
  // If chain_ext is greater than 0, the moved chain will be
  // forcefully renamed; the new name is composed as the previous
  // one + underscore + chain_ext (e.g. A_1). If thus generated
  // name duplicates any of existing chain IDs, or if chain_ext
  // was set to 0 and there is a duplication of chain IDs, the
  // name is again modified as above, with the extension number
  // generated automatically (this may result in IDs like
  // A_1_10).
  //   m_atom must give pointer to the Atom array, from which
  // the atoms belonging to m_chain, are moved to Atom array
  // given by 'atom', starting from poisition 'atom_index'.
  // 'atom_index' is then automatically updated to the next
  // free position in 'atom'.
  //   Note1: the moved atoms will occupy a continuous range
  // in 'atom' array; no checks on whether the corresponding
  // cells are occupied or not, are performed.
  //   Note2: the 'atom_index' is numbered from 0 on, i.e.
  // it is equal to atom[atom_index]->index-1; atom[]->index
  // is assigned automatically.
  ChainID   chainID;
  int       i,j,k,Ok;
  PPChain   chain1;
  PResidue  crRes;

    if (!m_chain)  return;

    // modify chain ID with the extension given
    if (chain_ext>0)
          sprintf ( chainID,"%s_%i",m_chain->chainID,chain_ext );
    else  strcpy  ( chainID,m_chain->chainID );

    // Choose the chain ID. If a chain with such ID is
    // already present in the model, it will be assigned
    // a new ID 'ID_n', where 'ID' stands for the original
    // chain ID and 'n' is the minimum (integer) number
    // chosen such that 'name_n' represents a new chain ID
    // (in the model).
    k = 0;
    do {
      Ok = true;
      for (i=0;(i<nChains) && (Ok);i++)
        if (chain[i])
          if (!strcmp(chainID,chain[i]->chainID))  Ok = false;
      if (!Ok)  {
        k++;
        if (chain_ext>0)
              sprintf ( chainID,"%s_%i_%i",m_chain->chainID,
                                           chain_ext,k );
        else  sprintf ( chainID,"%s_%i",m_chain->chainID,k );
      }
    } while (!Ok);

    // add chain on the top of Chain array.
    strcpy ( m_chain->chainID,chainID );
    if (nChains>=nChainsAlloc)  {
      nChainsAlloc = nChains+10;
      chain1 = new PChain[nChainsAlloc];
      k = 0;
      for (i=0;i<nChains;i++)
        if (chain[i])  chain1[k++] = chain[i];
      for (i=k;i<nChainsAlloc;i++)
        chain1[i] = NULL;
      if (chain)  delete[] chain;
      chain = chain1;
    }
    chain[nChains] = m_chain;
    chain[nChains]->SetModel ( this );
    nChains++;

    // Move all atoms of the chain. While residues belong
    // atoms belong to the chain's manager class. Therefore
    // they should be moved from one manager to another.
    for (i=0;i<m_chain->nResidues;i++)  {
      crRes = m_chain->residue[i];
      if (crRes)
        for (j=0;j<crRes->nAtoms;j++)
          if (crRes->atom[j])  {
            k = crRes->atom[j]->index-1;
            atom[atom_index] = m_atom[k];
            atom[atom_index]->index = atom_index+1;
            atom_index++;
            m_atom[k] = NULL;  // moved!
          }
    }

    m_chain = NULL;  // moved!

  }

  void Model::GetAIndexRange ( int & i1, int & i2 )  {
  PChain    chn;
  PResidue  res;
  int       ic,ir,ia;
    i1 = MaxInt4;
    i2 = MinInt4;
    for (ic=0;ic<nChains;ic++)  {
      chn = chain[ic];
      if (chn)  {
        for (ir=0;ir<chn->nResidues;ir++)  {
          res = chn->residue[ir];
          if (res)  {
            for (ia=0;ia<res->nAtoms;ia++)
              if (res->atom[ia])  {
                if (res->atom[ia]->index<i1)  i1 = res->atom[ia]->index;
                if (res->atom[ia]->index>i2)  i2 = res->atom[ia]->index;
              }
          }
        }
      }
    }

  }


  void  Model::MaskAtoms ( PMask Mask )  {
  int i;
    for (i=0;i<nChains;i++)
      if (chain[i])  chain[i]->MaskAtoms ( Mask );
  }

  void  Model::MaskResidues ( PMask Mask )  {
  int i;
    for (i=0;i<nChains;i++)
      if (chain[i])  chain[i]->MaskResidues ( Mask );
  }

  void  Model::MaskChains ( PMask Mask )  {
  int i;
    for (i=0;i<nChains;i++)
      if (chain[i])  chain[i]->SetMask ( Mask );
  }

  void  Model::UnmaskAtoms ( PMask Mask )  {
  int i;
    for (i=0;i<nChains;i++)
      if (chain[i])  chain[i]->UnmaskAtoms ( Mask );
  }

  void  Model::UnmaskResidues ( PMask Mask )  {
  int i;
    for (i=0;i<nChains;i++)
      if (chain[i])  chain[i]->UnmaskResidues ( Mask );
  }

  void  Model::UnmaskChains ( PMask Mask )  {
  int i;
    for (i=0;i<nChains;i++)
      if (chain[i])  chain[i]->RemoveMask ( Mask );
  }


  // ------ Getting Secondary Structure Elements

  int  Model::GetNumberOfHelices()  {
    return  helices.Length();
  }

  int  Model::GetNumberOfSheets()  {
    return  sheets.nSheets;
  }

  PHelix  Model::GetHelix ( int serialNum )  {
    return (PHelix)helices.GetContainerClass ( serialNum-1 );
  }

  void  Model::GetSheetID ( int serialNum, SheetID sheetID )  {
    if ((1<=serialNum) && (serialNum<=sheets.nSheets))  {
      if (sheets.sheet[serialNum-1])  {
        strcpy ( sheetID,sheets.sheet[serialNum-1]->sheetID );
        return;
      }
    }
    sheetID[0] = char(0);
  }

  PSheet Model::GetSheet ( int serialNum )  {
    if ((1<=serialNum) && (serialNum<=sheets.nSheets))
          return  sheets.sheet[serialNum-1];
    else  return  NULL;
  }

  PSheet Model::GetSheet ( const SheetID sheetID )  {
  int i;
    for (i=0;i<sheets.nSheets;i++)
      if (sheets.sheet[i])  {
        if (!strcmp(sheets.sheet[i]->sheetID,sheetID))
          return sheets.sheet[i];
      }
    return NULL;
  }

  int  Model::GetNumberOfStrands ( int sheetSerNum )  {
    if ((1<=sheetSerNum) && (sheetSerNum<=sheets.nSheets))  {
      if (sheets.sheet[sheetSerNum-1])
        return  sheets.sheet[sheetSerNum-1]->nStrands;
    }
    return 0;
  }

  int  Model::GetNumberOfStrands ( const SheetID sheetID )  {
  int i;
    for (i=0;i<sheets.nSheets;i++)
      if (sheets.sheet[i])  {
        if (!strcmp(sheets.sheet[i]->sheetID,sheetID))
          return sheets.sheet[i]->nStrands;
      }
    return 0;
  }

  PStrand Model::GetStrand ( int sheetSerNum, int strandSerNum )  {
  PSheet sheet;
    if ((1<=sheetSerNum) && (sheetSerNum<=sheets.nSheets))  {
      sheet = sheets.sheet[sheetSerNum-1];
      if (sheet)  {
        if ((1<=strandSerNum) && (strandSerNum<=sheet->nStrands))
        return  sheet->strand[strandSerNum-1];
      }
    }
    return NULL;
  }

  PStrand Model::GetStrand ( const SheetID sheetID,
                             int strandSerNum )  {
  int    i;
  PSheet sheet;
    for (i=0;i<sheets.nSheets;i++)
      if (sheets.sheet[i])  {
        if (!strcmp(sheets.sheet[i]->sheetID,sheetID))  {
          sheet = sheets.sheet[i];
          if (sheet)  {
            if ((1<=strandSerNum) && (strandSerNum<=sheet->nStrands))
              return  sheet->strand[strandSerNum-1];
          }
        }
      }
    return NULL;
  }

  void  Model::RemoveSecStructure()  {
    helices.FreeContainer();
    sheets .FreeMemory   ();
    turns  .FreeContainer();
  }

  void  Model::RemoveHetInfo()  {
    hetCompounds.FreeMemory();
  }


  int  Model::GetNumberOfLinks()  {
    return  links.Length();
  }

  PLink  Model::GetLink ( int serialNum )  {
    return (PLink)links.GetContainerClass ( serialNum-1 );
  }

  void  Model::RemoveLinks()  {
    links.FreeContainer();
  }

  void  Model::AddLink ( PLink link )  {
    links.AddData ( link );
  }


  int  Model::GetNumberOfLinkRs()  {
    return  linkRs.Length();
  }

  PLinkR  Model::GetLinkR ( int serialNum )  {
    return (PLinkR)linkRs.GetContainerClass ( serialNum-1 );
  }

  void  Model::RemoveLinkRs()  {
    linkRs.FreeContainer();
  }

  void  Model::AddLinkR ( PLinkR linkR )  {
    linkRs.AddData ( linkR );
  }



  int  Model::GetNumberOfCisPeps()  {
    return  cisPeps.Length();
  }

  PCisPep Model::GetCisPep ( int CisPepNum )  {
    return (PCisPep)cisPeps.GetContainerClass ( CisPepNum-1 );
  }

  void  Model::RemoveCisPeps()  {
    cisPeps.FreeContainer();
  }

  void  Model::AddCisPep ( PCisPep cisPep )  {
    cisPeps.AddData ( cisPep );
  }


  void  Model::ApplyTransform ( mat44 & TMatrix )  {
  // transforms all coordinates by multiplying with matrix TMatrix
  int i;
    for (i=0;i<nChains;i++)
      if (chain[i])  chain[i]->ApplyTransform ( TMatrix );
  }

  bool Model::isInSelection ( int selHnd )  {
  PMask  mask;
    if (manager)  {
      mask = PRoot(manager)->GetSelMask ( selHnd );
      if (mask)  return CheckMask ( mask );
    }
    return false;
  }



  // -------  user-defined data handlers

  int  Model::PutUDData ( int UDDhandle, int iudd )  {
    if (UDDhandle & UDRF_MODEL)
          return  UDData::putUDData ( UDDhandle,iudd );
    else  return  UDDATA_WrongUDRType;
  }

  int  Model::PutUDData ( int UDDhandle, realtype rudd )  {
    if (UDDhandle & UDRF_MODEL)
          return  UDData::putUDData ( UDDhandle,rudd );
    else  return  UDDATA_WrongUDRType;
  }

  int  Model::PutUDData ( int UDDhandle, cpstr sudd )  {
    if (UDDhandle & UDRF_MODEL)
          return  UDData::putUDData ( UDDhandle,sudd );
    else  return  UDDATA_WrongUDRType;
  }

  int  Model::GetUDData ( int UDDhandle, int & iudd )  {
    if (UDDhandle & UDRF_MODEL)
          return  UDData::getUDData ( UDDhandle,iudd );
    else  return  UDDATA_WrongUDRType;
  }

  int  Model::GetUDData ( int UDDhandle, realtype & rudd )  {
    if (UDDhandle & UDRF_MODEL)
          return  UDData::getUDData ( UDDhandle,rudd );
    else  return  UDDATA_WrongUDRType;
  }

  int  Model::GetUDData ( int UDDhandle, pstr sudd, int maxLen )  {
    if (UDDhandle & UDRF_MODEL)
          return  UDData::getUDData ( UDDhandle,sudd,maxLen );
    else  return  UDDATA_WrongUDRType;
  }

  int  Model::GetUDData ( int UDDhandle, pstr & sudd )  {
    if (UDDhandle & UDRF_MODEL)
          return  UDData::getUDData ( UDDhandle,sudd );
    else  return  UDDATA_WrongUDRType;
  }


  // -------  calculation of Secondary Structure


  int Model::CalcSecStructure ( bool flagBulge, int aminoSelHnd )  {
  // This function is contributed by Liz Potterton, University of York
  //------------------------------------------------------------------
  // Define a secondary structure type of each amino acid residue in the
  // structure.
  // Procedure:
  // Find all amino acids
  // Find all pairs of amino acids which have inter-Ca distance  < 10.0A
  // Test for hydrogen bonds between the main chain N and O of the close
  // residues and store the information in the hbonds matrix
  // Analyse the info in hbonds matrix to assign secondary structure to
  // secstr vector
  PPResidue Res;
  PPAtom    Ca;
  PChain    chn;
  PContact  contact;
  imatrix   hbonds;
  PPAtom *  hbond_atoms;
  int       nres, ncontacts;
  int       ir1,ir2, irdif;
  int       i,j,k,l;

    // 1a. Get protein residues from selection handle

    if (aminoSelHnd>=0) {

      manager->GetSelIndex(aminoSelHnd,Res,nres);
      if (nres<=0)  return  SSERC_noResidues;

    } else {

      //  1b. Get all protein residues

      nres = 0;
      for (i=0;i<nChains;i++)
        if (chain[i])
          nres += chain[i]->nResidues;

      if (nres<=0)  return  SSERC_noResidues;

      Res  = new PResidue[nres];
      nres = 0;
      for (i=0;i<nChains;i++)  {
        chn = chain[i];
        if (chn)  {
          k = chn->nResidues;
          for (j=0;j<k;j++)
            Res[nres++] = chn->residue[j];
        }
      }


      if (nres<=0)  {
        delete[] Res;
        return   SSERC_noResidues;
      }

   }

    //  2. Get C-alphas of all aminoacids

    Ca = new PAtom[nres];
    k  = 0;
    for (i=0;i<nres;i++)
      if (Res[i])  {
        if (aminoSelHnd>=0 || Res[i]->isAminoacid())  {
          Ca[i] = Res[i]->GetAtom("CA", " C", "*");
          k++;
        } else
          Ca[i] = NULL;
        Res[i]->SSE = SSE_None;
      } else
        Ca[i] = NULL;

    if (k<=0)  {
      delete[] Res;
      delete[] Ca;
      return   SSERC_noAminoacids;
    }


    //  3. Find all close Calphas - i.e. find the contacts between
    //     the two equivalent sets of Ca atoms

    contact   = NULL;
    ncontacts = 0;
    manager->SeekContacts ( Ca,nres, Ca,nres, 2.0,10.0, 2,
                            contact,ncontacts,0 );
    manager->RemoveBricks();
    if (ncontacts<=0)  {
      delete[] Res;
      delete[] Ca;
      if (contact)  delete[] contact;
      return  SSERC_noSSE;
    }


    //  4. Get and initialize memory for analysing the SSE

    GetMatrixMemory ( hbonds,nres,3,0,0 );
    hbond_atoms = new PPAtom[nres];
    for (i=0;i<nres;i++)  {
      hbond_atoms[i] = new PAtom[6];
      for (j=0;j<6;j++) hbond_atoms[i][j] = NULL;
      for (j=0;j<3;j++) hbonds     [i][j] = 0;
    }


    //  5.  Loop over all close (in space) residues - excluding those
    //      that are close in sequence

    for (i=0;i<ncontacts;i++)  {
      ir1   = contact[i].id2;
      ir2   = contact[i].id1;
      irdif = ir1 - ir2;
      if (irdif>2)  {
        //  test if there is donor Hbond from residue ir1
        if (Res[ir1]->isMainchainHBond(Res[ir2]))  {
          k = 0;
          while ((hbonds[ir1][k]!=0) && (k<2))  k++;
          hbonds     [ir1][k]   = -irdif;
          hbond_atoms[ir1][k]   = Res[ir1]->GetAtom ( "N" );
          hbond_atoms[ir1][k+3] = Res[ir2]->GetAtom ( "O" );
        }
        //  test if there is donor Hbond from residue ir2
        if (Res[ir2]->isMainchainHBond(Res[ir1]))  {
          k = 0;
          while ((hbonds[ir2][k]!=0) && (k<2))  k++;
          hbonds     [ir2][k]   = irdif;
          hbond_atoms[ir2][k]   = Res[ir2]->GetAtom ( "N" );
          hbond_atoms[ir2][k+3] = Res[ir1]->GetAtom ( "O" );
        }
      }
    }

    //  6. Assign the turns - if there is bifurcated bond then the 4-turn
    //     takes precedence - read the paper to make sense of this

    for (i=0;i<nres;i++)  {
      k = 0;
      while ((k<=2) && (hbonds[i][k]!=0))  {
        if (hbonds[i][k]==-5)  {
          Res[i-1]->SSE = SSE_5Turn;
          Res[i-2]->SSE = SSE_5Turn;
          Res[i-3]->SSE = SSE_5Turn;
          Res[i-4]->SSE = SSE_5Turn;
        }
        if (hbonds[i][k]==-3)  {
          Res[i-1]->SSE = SSE_3Turn;
          Res[i-2]->SSE = SSE_3Turn;
        }
        k++;
      }
    }
    for (i=0;i<nres;i++)  {
      k = 0;
      while ((k<=2) && (hbonds[i][k]!=0))  {
        if (hbonds[i][k]==-4)  {
          Res[i-1]->SSE = SSE_4Turn;
          Res[i-2]->SSE = SSE_4Turn;
          Res[i-3]->SSE = SSE_4Turn;
        }
        k++;
      }
    }


    //  7. Look for consecutive 4-turns which make alpha helix

    for (i=1;i<nres-3;i++) {
      if (((Res[i  ]->SSE==SSE_Helix) || (Res[i  ]->SSE==SSE_4Turn)) &&
          ((Res[i+1]->SSE==SSE_Helix) || (Res[i+1]->SSE==SSE_4Turn)) &&
          ((Res[i+2]->SSE==SSE_Helix) || (Res[i+2]->SSE==SSE_4Turn)) &&
          ((Res[i+3]->SSE==SSE_Helix) || (Res[i+3]->SSE==SSE_4Turn)))
        for (j=i;j<=i+3;j++)  Res[j]->SSE = SSE_Helix;
    }

    for (i=0;i<nres;i++)  {

      k = 0;
      while ((k<=2) && (hbonds[i][k]!=0))  {

        irdif = hbonds[i][k];
        // Test for 'close' hbond
        j = i + irdif;
        l = 0;
        while ((l<=2) && (hbonds[j][l]!=0))  {
          // Antiparallel strands
          if (hbonds[j][l]==-irdif)  {
            Res[i]->SSE = SSE_Strand;
            Res[j]->SSE = SSE_Strand;
          }
          // Parallel strand
          if (hbonds[j][l]==-irdif-2)  {
            Res[i-1]->SSE = SSE_Strand;
            Res[j  ]->SSE = SSE_Strand;
          }
          // Parallel beta bulge
          if (hbonds[j][l]==-irdif-3)  {
            if (flagBulge) {
              if (Res[i-1]->SSE==SSE_None)  Res[i-1]->SSE = SSE_Bulge;
              if (Res[i-2]->SSE==SSE_None)  Res[i-2]->SSE = SSE_Bulge;
              if (Res[j  ]->SSE==SSE_None)  Res[j  ]->SSE = SSE_Bulge;
            } else  {
              if (Res[i-1]->SSE==SSE_None)  Res[i-1]->SSE = SSE_Strand;
              if (Res[i-2]->SSE==SSE_None)  Res[i-2]->SSE = SSE_Strand;
              if (Res[j  ]->SSE==SSE_None)  Res[j  ]->SSE = SSE_Strand;
            }
          }
          l++;
        }
        // Test for 'wide' hbond
        j = i + hbonds[i][k] + 2;
        if (j<nres)  {
          l = 0;
          while ((l<=2) && (hbonds[j][l]!=0))  {
            // Antiaprallel strands
            if (hbonds[j][l]==-irdif-4)  {
              Res[i-1]->SSE = SSE_Strand;
              Res[j-1]->SSE = SSE_Strand;
            }
            // Parallel strands
            if (hbonds[j][l]==-irdif-2)  {
              Res[i  ]->SSE = SSE_Strand;
          Res[j-1]->SSE = SSE_Strand;
            }
            l++;
          }
        }

        // test for anti-parallel B-bulge between 'close' hbonds
        j = i + hbonds[i][k] - 1;
        if (j>=0)  {
          l = 0;
          while ((l<=2) && (hbonds[j][l]!=0))  {
            if (hbonds[j][l]==-irdif+1)  {
              if (flagBulge)  {
            if (Res[i  ]->SSE==SSE_None)  Res[i  ]->SSE = SSE_Bulge;
            if (Res[j+1]->SSE==SSE_None)  Res[j+1]->SSE = SSE_Bulge;
            if (Res[j  ]->SSE==SSE_None)  Res[j  ]->SSE = SSE_Bulge;
              } else  {
                if (Res[i  ]->SSE==SSE_None)  Res[i  ]->SSE = SSE_Strand;
                if (Res[j+1]->SSE==SSE_None)  Res[j+1]->SSE = SSE_Strand;
                if (Res[j  ]->SSE==SSE_None)  Res[j  ]->SSE = SSE_Strand;
              }
            }
            l++;
          }
        }

        // test for anti-parallel B-bulge between 'wide' hbonds
        j = i + hbonds[i][k] + 3;
        if (j<nres)  {
          l = 0;
          while ((l<=2) && (hbonds[j][l]!=0))  {
            if ((hbonds[j][l]==-irdif+5) && (i>0))  {
              if (flagBulge)  {
                if (Res[i-1]->SSE==SSE_None)  Res[i-1]->SSE = SSE_Bulge;
                if (Res[j-1]->SSE==SSE_None)  Res[j-1]->SSE = SSE_Bulge;
                if (Res[j-2]->SSE==SSE_None)  Res[j-2]->SSE = SSE_Bulge;
              } else  {
                if (Res[i-1]->SSE==SSE_None)  Res[i-1]->SSE = SSE_Strand;
                if (Res[j-1]->SSE==SSE_None)  Res[j-1]->SSE = SSE_Strand;
                if (Res[j-2]->SSE==SSE_None)  Res[j-2]->SSE = SSE_Strand;
              }
            } else if (hbonds[j][l]==-irdif-3)  {
              // and bulge in parallel strand
          if (flagBulge)  {
                if (Res[i  ]->SSE==SSE_None)  Res[i  ]->SSE = SSE_Bulge;
                if (Res[j-1]->SSE==SSE_None)  Res[j-1]->SSE = SSE_Bulge;
                if (Res[j-2]->SSE==SSE_None)  Res[j-2]->SSE = SSE_Bulge;
              }
              else {
                if (Res[i  ]->SSE==SSE_None)  Res[i  ]->SSE = SSE_Strand;
                if (Res[j-1]->SSE==SSE_None)  Res[j-1]->SSE = SSE_Strand;
                if (Res[j-2]->SSE==SSE_None)  Res[j-2]->SSE = SSE_Strand;
              }
            }
            l++;
          }
        }
        k++;

      } // Finish looping over Hbonds for residue (k loop)

    }  // Finish looping over residues ( i loop)


    //  8. Free memory

    if (hbond_atoms)  {
      for (i=0;i<nres;i++)
        if (hbond_atoms[i])  delete[] hbond_atoms[i];
      delete[] hbond_atoms;
    }
    FreeMatrixMemory ( hbonds,nres,0,0 );
    if (contact) delete[] contact;
    if (Res && aminoSelHnd<0) delete[] Res;
    if (Ca)      delete[] Ca;

    return  SSERC_Ok;

  }


  // -------  streaming

  void  Model::write ( io::RFile f )  {
  int  i,k;
  byte Version=4;
  bool compactBinary = false;
  
    PManager M = GetCoordHierarchy();
    if (M)
      compactBinary = M->isCompactBinary();

    f.WriteByte ( &Version       );
    f.WriteBool ( &compactBinary );
    
    f.WriteInt ( &serNum  );
    f.WriteInt ( &nChains );

    for (i=0;i<nChains;i++)  {
      if (chain[i])  k = 1;
               else  k = 0;
      f.WriteInt ( &k );
      if (chain[i]) chain[i]->write ( f );
    }

    if (!compactBinary)  {

      ProModel::write ( f );
  
      hetCompounds.write ( f );
      helices     .write ( f );
      sheets      .write ( f );
      turns       .write ( f );
      links       .write ( f );
      linkRs      .write ( f );
      
    }

  }

  void  Model::read ( io::RFile f )  {
  int  i,k;
  byte Version;
  bool compactBinary;

    FreeMemory();

    f.ReadByte ( &Version       );
    f.ReadBool ( &compactBinary );
    
    f.ReadInt ( &serNum  );
    f.ReadInt ( &nChains );
    nChainsAlloc = nChains;
    if (nChains>0)  {
      chain = new PChain[nChainsAlloc];
      for (i=0;i<nChains;i++)  {
        f.ReadInt ( &k );
        if (k)  {
          chain[i] = newChain();
          chain[i]->SetModel ( this );
          chain[i]->read ( f );
        }
      }
    }

    if (!compactBinary)  {
      
      ProModel::read ( f );
  
      hetCompounds.read ( f );
      helices     .read ( f );
      sheets      .read ( f );
      turns       .read ( f );
      links       .read ( f );
      linkRs      .read ( f );
      
    }

  }

  MakeFactoryFunctions(Model)

}  // namespace mmdb
