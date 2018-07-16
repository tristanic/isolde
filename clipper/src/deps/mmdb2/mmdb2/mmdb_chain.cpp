//  $Id: mmdb_chain.h $
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
//    23.12.15   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  MMDB_Chain <implementation>
//       ~~~~~~~~~
//  **** Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::ProModel     ( a virtue of Model           )
//       ~~~~~~~~~  mmdb::DBReference  ( DBREF  records               )
//             mmdb::ChainContainer ( container of in-chain classes   )
//             mmdb::ContainerChain ( chain containered class template)
//             mmdb::SeqAdv         ( SEQADV records                  )
//             mmdb::SeqRes         ( SEQRES records                  )
//             mmdb::ModRes         ( MODRES records                  )
//             mmdb::HetRec         ( HET    records                  )
//             mmdb::Chain          ( chain class                     )
//
//  Copyright (C) E. Krissinel 2000-2015
//
//  =================================================================
//

#include <string.h>
#include <stdlib.h>

#include "mmdb_chain.h"
#include "mmdb_model.h"
#include "mmdb_manager.h"
#include "mmdb_cifdefs.h"
#include "mmdb_tables.h"

namespace mmdb  {

  //  ==================  ProModel  ======================

  MakeStreamFunctions(ProModel)

  //  ==============  ChainContainer  ====================

  PContainerClass ChainContainer::MakeContainerClass ( int ClassID )  {
    switch (ClassID)  {
      default :
      case ClassID_Template    : return
                                   ClassContainer::MakeContainerClass(ClassID);
      case ClassID_DBReference : return new DBReference ( chain );
      case ClassID_SeqAdv      : return new SeqAdv      ( chain );
      case ClassID_ModRes      : return new ModRes      ( chain );
      case ClassID_Het         : return new HetRec      ( chain );
    }
  }

  void ChainContainer::SetChain ( PChain Chain_Owner )  {
  int i;
    chain = Chain_Owner;
    for (i=0;i<length;i++)
      if (Container[i])
        (void)PContainerChain(Container[i])->SetChain ( chain );
  }

  cpstr ChainContainer::Get1stChainID()  {
  int i;
    i = 0;
    if (Container)  {
      while ((i<length-1) && (!Container[i])) i++;
      if (Container[i])
            return PContainerChain(Container[i])->chainID;
      else  return NULL;
    } else
      return NULL;
  }

  void ChainContainer::MoveByChainID ( const ChainID chainID,
                                       PChainContainer ChainContainer ) {
  int i;
    for (i=0;i<length;i++)
      if (Container[i])  {
        if (!strcmp(PContainerChain(Container[i])->chainID,chainID))  {
          ChainContainer->AddData ( Container[i] );
          Container[i] = NULL;
        }
      }
  }


  MakeStreamFunctions(ChainContainer)


  //  ================  ContainerChain  ===================

  ContainerChain::ContainerChain() : ContainerClass()  {
    chain      = NULL;
    chainID[0] = char(0);
  }

  ContainerChain::ContainerChain ( PChain Chain_Owner)
                : ContainerClass()  {
    chain = Chain_Owner;
    if (chain)  strcpy ( chainID,chain->GetChainID() );
          else  chainID[0] = char(0);
  }

  void ContainerChain::SetChain ( PChain Chain_Owner )  {
    chain = Chain_Owner;
    if (chain)  strcpy ( chainID,chain->GetChainID() );
          else  strcpy ( chainID,"" );
  }

  MakeStreamFunctions(ContainerChain)


  //  ================  DBReference  ===================

  DBReference::DBReference() : ContainerChain()  {
    InitDBReference();
  }

  DBReference::DBReference( PChain Chain_Owner )
             : ContainerChain(Chain_Owner)  {
    InitDBReference();
  }

  DBReference::DBReference ( PChain Chain_Owner, cpstr S )
             : ContainerChain(Chain_Owner)  {
    InitDBReference();
    ConvertPDBASCII ( S );
  }

  DBReference::DBReference ( io::RPStream Object )
             : ContainerChain(Object)  {
    InitDBReference();
  }

  DBReference::~DBReference() {}

  void  DBReference::InitDBReference()  {
    seqBeg = 0;
    strcpy ( insBeg     ,"-"            );
    seqEnd = 0;
    strcpy ( insEnd     ,"-"            );
    strcpy ( database   ,"------"       );
    strcpy ( dbAccession,"--------"     );
    strcpy ( dbIdCode   ,"------------" );
    dbseqBeg = 0;
    strcpy ( dbinsBeg,"-" );
    dbseqEnd = 0;
    strcpy ( dbinsEnd,"-" );
  }

  void  DBReference::PDBASCIIDump ( pstr S, int N )  {
  UNUSED_ARGUMENT(N);
  //  makes the ASCII PDB DBREF line number N
  //  from the class' data
    strcpy ( S,"DBREF" );
    PadSpaces ( S,80 );
    strcpy_n  ( &(S[7]),chain->GetEntryID(),4 );
    if (chain->chainID[0])  S[12] = chain->chainID[0];
    PutIntIns ( &(S[14]),seqBeg,4,insBeg     );
    PutIntIns ( &(S[20]),seqEnd,4,insEnd     );
    strcpy_n  ( &(S[26]),database   ,6       );
    strcpy_n  ( &(S[33]),dbAccession,8       );
    strcpy_n  ( &(S[42]),dbIdCode   ,12      );
    PutIntIns ( &(S[55]),dbseqBeg,5,dbinsBeg );
    PutIntIns ( &(S[62]),dbseqEnd,5,dbinsEnd );
  }

  void  DBReference::MakeCIF ( mmcif::PData CIF, int N )  {
  UNUSED_ARGUMENT(N);
  mmcif::PLoop Loop1,Loop2;
  int          RC1,RC2;

    RC1 = CIF->AddLoop ( CIFCAT_STRUCT_REF_SEQ,Loop1 );
    RC2 = CIF->AddLoop ( CIFCAT_STRUCT_REF    ,Loop2 );

    if ((RC1!=mmcif::CIFRC_Ok) || (RC2!=mmcif::CIFRC_Ok))  {
      // the category was (re)created, provide tags
      Loop1->AddLoopTag ( CIFTAG_NDB_PDB_ID_CODE            );
      Loop1->AddLoopTag ( CIFTAG_NDB_CHAIN_ID               );
      Loop1->AddLoopTag ( CIFTAG_SEQ_ALIGN_BEG              );
      Loop1->AddLoopTag ( CIFTAG_NDB_SEQ_ALIGN_BEG_INS_CODE );
      Loop1->AddLoopTag ( CIFTAG_SEQ_ALIGN_END              );
      Loop1->AddLoopTag ( CIFTAG_NDB_SEQ_ALIGN_END_INS_CODE );
      Loop1->AddLoopTag ( CIFTAG_NDB_DB_ACCESSION           );
      Loop1->AddLoopTag ( CIFTAG_DB_ALIGN_BEG               );
      Loop1->AddLoopTag ( CIFTAG_NDB_DB_ALIGN_BEG_INS_CODE  );
      Loop1->AddLoopTag ( CIFTAG_DB_ALIGN_END               );
      Loop1->AddLoopTag ( CIFTAG_NDB_DB_ALIGN_END_INS_CODE  );
      Loop2->AddLoopTag ( CIFTAG_DB_NAME );
      Loop2->AddLoopTag ( CIFTAG_DB_CODE );
    }

    Loop1->AddString  ( chain->GetEntryID(),true );
    Loop1->AddString  ( chain->chainID     ,true );
    Loop1->AddInteger ( seqBeg                   );
    Loop1->AddString  ( insBeg             ,true );
    Loop1->AddInteger ( seqEnd                   );
    Loop1->AddString  ( insEnd             ,true );
    Loop1->AddString  ( dbAccession        ,true );
    Loop1->AddInteger ( dbseqBeg                 );
    Loop1->AddString  ( dbinsBeg           ,true );
    Loop1->AddInteger ( dbseqEnd                 );
    Loop1->AddString  ( dbinsEnd           ,true );

    Loop2->AddString  ( database,true );
    Loop2->AddString  ( dbIdCode,true );

  }

  ERROR_CODE DBReference::GetCIF ( mmcif::PData CIF, int & n )  {
  //  GetCIF(..) must be always run without reference to Chain,
  //  see CModel::GetCIF(..).
  mmcif::PLoop   Loop1,Loop2;
  mmcif::PStruct Struct2;
  pstr           F;
  int            RC,ref_id1,ref_id2;
  CIF_MODE       CIFMode;
  ERROR_CODE     rc;

    Loop1 = CIF->GetLoop ( CIFCAT_STRUCT_REF_SEQ );

    if (!Loop1)  {
      n = -1;
      return Error_EmptyCIF;
    }

    if (n>=Loop1->GetLoopLength())  {
      n = -1;
      return Error_EmptyCIF;
    }


    //  Determine the ChainID first and store it locally. It will
    // be used by CModel for generating chains and placing the
    // primary structure data BEFORE reading the coordinate section.
    CIFMode = CIF_NDB;
    F = Loop1->GetString ( CIFName(TAG_CHAIN_ID,CIFMode),n,RC );
    if ((RC) || (!F))  {
      CIFMode = CIF_PDBX;
      F = Loop1->GetString ( CIFName(TAG_CHAIN_ID,CIFMode),n,RC );
    }
    if ((!RC) && F)  {
      strcpy_n0 ( chainID,F,sizeof(ChainID)-1 );
      Loop1->DeleteField ( CIFName(TAG_CHAIN_ID,CIFMode),n );
    } else
      strcpy ( chainID,"" );


    rc = CIFGetInteger(seqBeg,Loop1,CIFName(TAG_SEQ_ALIGN_BEG,CIFMode),n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;
    CIFGetString ( insBeg,Loop1,CIFName(TAG_SEQ_ALIGN_BEG_INS_CODE,CIFMode),
                   n,sizeof(InsCode),pstr(" ") );

    rc = CIFGetInteger(seqEnd,Loop1,CIFName(TAG_SEQ_ALIGN_END,CIFMode),n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;
    CIFGetString ( insEnd,Loop1,CIFName(TAG_SEQ_ALIGN_END_INS_CODE,CIFMode),
                   n,sizeof(InsCode),pstr(" ") );
    CIFGetString ( dbAccession,Loop1,CIFName(TAG_DB_ACCESSION,CIFMode),
                   n,sizeof(DBAcCode),pstr("        ") );

    rc = CIFGetInteger(dbseqBeg,Loop1,CIFName(TAG_DB_ALIGN_BEG,CIFMode),n);
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;
    CIFGetString ( dbinsBeg,Loop1,CIFName(TAG_DB_ALIGN_BEG_INS_CODE,CIFMode),
                   n,sizeof(InsCode),pstr(" ") );

    rc = CIFGetInteger(dbseqEnd,Loop1,CIFName(TAG_DB_ALIGN_END,CIFMode),n);
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;
    CIFGetString ( dbinsEnd,Loop1,CIFName(TAG_DB_ALIGN_END_INS_CODE,CIFMode),
                   n,sizeof(InsCode),pstr(" ") );

    Loop2 = CIF->GetLoop ( CIFCAT_STRUCT_REF );
    if (Loop2)  {
      CIFGetString ( database,Loop2,CIFTAG_DB_NAME,n,
                     sizeof(DBName)  ,pstr("      ")       );
      CIFGetString ( dbIdCode,Loop2,CIFTAG_DB_CODE,n,
                     sizeof(DBIdCode),pstr("            ") );
    } else if (CIFMode==CIF_PDBX)  {
      Struct2 = CIF->GetStructure ( CIFCAT_STRUCT_REF );
      if (Struct2 &&
          (!CIFGetInteger(ref_id1,Loop1,CIFTAG_REF_ID,n)) &&
          (!CIFGetInteger(ref_id2,Struct2,CIFTAG_ID,false)))  {
        if (ref_id1==ref_id2)  {
          CIFGetString ( database,Struct2,CIFTAG_DB_NAME,
                         sizeof(DBName)  ,pstr("      ")      ,false );
          CIFGetString ( dbIdCode,Struct2,CIFTAG_DB_CODE,
                         sizeof(DBIdCode),pstr("            "),false );
        }
      }
    }

    n++;

    return Error_NoError;

  }


  ERROR_CODE DBReference::ConvertPDBASCII ( cpstr S )  {
  IDCode idCode;
    if (chain->chainID[0])  {
      if (S[12]!=chain->chainID[0])
        return Error_WrongChainID;
    } else if (S[12]!=' ')  {
      chain->chainID[0] = S[12];
      chain->chainID[1] = char(0);
    } else
      chain->chainID[0] = char(0);
    strcpy ( idCode,chain->GetEntryID() );
    if (idCode[0])  {
      if (strncmp(&(S[7]),idCode,4) && (!ignoreNonCoorPDBErrors))
        return Error_WrongEntryID;
    } else  {
      GetString ( idCode,&(S[7]),4 );
      chain->SetEntryID ( idCode );
    }
    GetIntIns  ( seqBeg,insBeg,&(S[14]),4  );
    GetIntIns  ( seqEnd,insEnd,&(S[20]),4  );
    strcpy_ncs ( database     ,&(S[26]),6  );
    strcpy_ncs ( dbAccession  ,&(S[33]),8  );
    strcpy_ncs ( dbIdCode     ,&(S[42]),12 );
    GetIntIns  ( dbseqBeg,dbinsBeg,&(S[55]),5 );
    GetIntIns  ( dbseqEnd,dbinsEnd,&(S[62]),5 );
    return Error_NoError;
  }

  void  DBReference::Copy ( PContainerClass DBRef )  {

    ContainerChain::Copy ( DBRef );

    seqBeg   = PDBReference(DBRef)->seqBeg;
    seqEnd   = PDBReference(DBRef)->seqEnd;
    dbseqBeg = PDBReference(DBRef)->dbseqBeg;
    dbseqEnd = PDBReference(DBRef)->dbseqEnd;
    strcpy ( insBeg     ,PDBReference(DBRef)->insBeg      );
    strcpy ( insEnd     ,PDBReference(DBRef)->insEnd      );
    strcpy ( database   ,PDBReference(DBRef)->database    );
    strcpy ( dbAccession,PDBReference(DBRef)->dbAccession );
    strcpy ( dbIdCode   ,PDBReference(DBRef)->dbIdCode    );
    strcpy ( dbinsBeg   ,PDBReference(DBRef)->dbinsBeg    );
    strcpy ( dbinsEnd   ,PDBReference(DBRef)->dbinsEnd    );

  }

  void  DBReference::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte ( &Version  );
    f.WriteInt  ( &seqBeg   );
    f.WriteInt  ( &seqEnd   );
    f.WriteInt  ( &dbseqBeg );
    f.WriteInt  ( &dbseqEnd );
    f.WriteTerLine ( insBeg     ,false );
    f.WriteTerLine ( insEnd     ,false );
    f.WriteTerLine ( database   ,false );
    f.WriteTerLine ( dbAccession,false );
    f.WriteTerLine ( dbIdCode   ,false );
    f.WriteTerLine ( dbinsBeg   ,false );
    f.WriteTerLine ( dbinsEnd   ,false );
  }

  void  DBReference::read  ( io::RFile f ) {
  byte Version;
    f.ReadByte ( &Version  );
    f.ReadInt  ( &seqBeg   );
    f.ReadInt  ( &seqEnd   );
    f.ReadInt  ( &dbseqBeg );
    f.ReadInt  ( &dbseqEnd );
    f.ReadTerLine ( insBeg     ,false );
    f.ReadTerLine ( insEnd     ,false );
    f.ReadTerLine ( database   ,false );
    f.ReadTerLine ( dbAccession,false );
    f.ReadTerLine ( dbIdCode   ,false );
    f.ReadTerLine ( dbinsBeg   ,false );
    f.ReadTerLine ( dbinsEnd   ,false );
  }

  MakeStreamFunctions(DBReference)



  //  ================  SeqAdv  ===================

  SeqAdv::SeqAdv() : ContainerChain()  {
    InitSeqAdv();
  }

  SeqAdv::SeqAdv ( PChain Chain_Owner )
         : ContainerChain(Chain_Owner)  {
    InitSeqAdv();
  }

  SeqAdv::SeqAdv ( PChain Chain_Owner, cpstr S )
         : ContainerChain(Chain_Owner)  {
    InitSeqAdv();
    ConvertPDBASCII ( S );
  }

  SeqAdv::SeqAdv ( io::RPStream Object ) : ContainerChain(Object)  {
    InitSeqAdv();
  }

  SeqAdv::~SeqAdv()  {
    if (conflict)  delete[] conflict;
  }

  void  SeqAdv::InitSeqAdv()  {
    strcpy ( resName    ,"---"       );
    seqNum = 0;
    strcpy ( insCode    ,"-"         );
    strcpy ( database   ,"------"    );
    strcpy ( dbAccession,"---------" );
    strcpy ( dbRes      ,"---"       );
    dbSeq = 0;
    conflict = NULL;
    CreateCopy ( conflict,pstr(" ") );
  }

  void  SeqAdv::PDBASCIIDump ( pstr S, int N )  {
  UNUSED_ARGUMENT(N);
  //  makes the ASCII PDB SEQADV line number N
  //  from the class' data
    strcpy     ( S,"SEQADV" );
    PadSpaces  ( S,80 );
    strcpy_n   ( &(S[7]) ,chain->GetEntryID(),4 );
    strcpy_n   ( &(S[12]),resName      ,3 );
    if (chain->chainID[0])  S[16] = chain->chainID[0];
    PutIntIns  ( &(S[18]),seqNum,4,insCode );
    strcpy_n   ( &(S[24]),database   ,4    );
    strcpy_n   ( &(S[29]),dbAccession,9    );
    strcpy_n   ( &(S[39]),dbRes      ,3    );
    PutInteger ( &(S[43]),dbSeq      ,5    );
    strcpy_n   ( &(S[49]),conflict,IMin(strlen(conflict),21) );
  }

  ERROR_CODE SeqAdv::ConvertPDBASCII ( cpstr S )  {
  IDCode idCode;
    if (chain->chainID[0])  {
      if (S[16]!=chain->chainID[0])
        return Error_WrongChainID;
    } else if (S[16]!=' ')  {
      chain->chainID[0] = S[16];
      chain->chainID[1] = char(0);
    } else
      chain->chainID[0] = char(0);
    strcpy ( idCode,chain->GetEntryID() );
    if (idCode[0])  {
      if (strncmp(&(S[7]),idCode,4) && (!ignoreNonCoorPDBErrors))
        return Error_WrongEntryID;
    } else  {
      GetString ( idCode,&(S[7]),4 );
      chain->SetEntryID ( idCode );
    }
    strcpy_ncs ( resName       ,&(S[12]),3 );
    GetIntIns  ( seqNum,insCode,&(S[18]),4 );
    strcpy_ncs ( database      ,&(S[24]),4 );
    strcpy_ncs ( dbAccession   ,&(S[29]),9 );
    strcpy_ncs ( dbRes         ,&(S[39]),3 );
    GetInteger ( dbSeq,&(S[43]),5  );
    CreateCopy ( conflict,&(S[49]) );
    CutSpaces  ( conflict,SCUTKEY_END );
    return Error_NoError;
  }


  void  SeqAdv::MakeCIF ( mmcif::PData CIF, int N )  {
  UNUSED_ARGUMENT(N);
  mmcif::PLoop Loop;
  int          RC;

    RC = CIF->AddLoop ( CIFCAT_STRUCT_REF_SEQ_DIF,Loop );

    if (RC!=mmcif::CIFRC_Ok)  {
      // the category was (re)created, provide tags
      Loop->AddLoopTag ( CIFTAG_NDB_PDB_ID_CODE           );
      Loop->AddLoopTag ( CIFTAG_MON_ID                    );
      Loop->AddLoopTag ( CIFTAG_NDB_PDB_CHAIN_ID          );
      Loop->AddLoopTag ( CIFTAG_SEQ_NUM                   );
      Loop->AddLoopTag ( CIFTAG_NDB_PDB_INS_CODE          );
      Loop->AddLoopTag ( CIFTAG_NDB_SEQ_DB_NAME           );
      Loop->AddLoopTag ( CIFTAG_NDB_SEQ_DB_ACCESSION_CODE );
      Loop->AddLoopTag ( CIFTAG_DB_MON_ID                 );
      Loop->AddLoopTag ( CIFTAG_NDB_SEQ_DB_SEQ_NUM        );
      Loop->AddLoopTag ( CIFTAG_DETAILS                   );
    }

    Loop->AddString  ( chain->GetEntryID(),true );
    Loop->AddString  ( resName            ,true );
    Loop->AddString  ( chain->chainID     ,true );
    Loop->AddInteger ( seqNum                   );
    Loop->AddString  ( insCode            ,true );
    Loop->AddString  ( database           ,true );
    Loop->AddString  ( dbAccession        ,true );
    Loop->AddString  ( dbRes              ,true );
    Loop->AddInteger ( dbSeq                    );
    Loop->AddString  ( conflict           ,true );

  }

  ERROR_CODE SeqAdv::GetCIF ( mmcif::PData CIF, int & n )  {
  //  GetCIF(..) must be always run without reference to Chain,
  //  see CModel::GetCIF(..).
  mmcif::PLoop  Loop;
  pstr          F;
  int           RC;

    Loop = CIF->GetLoop ( CIFCAT_STRUCT_REF_SEQ_DIF );
    if (!Loop)  {
      n = -1;
      return Error_EmptyCIF;
    }

    if (n>=Loop->GetLoopLength())  {
      n = -1;
      return Error_EmptyCIF;
    }

    //  Determine the ChainID first and store it locally. It will
    // be used by CModel for generating chains and placing the
    // primary structure data BEFORE reading the coordinate section.

    F = Loop->GetString ( CIFTAG_NDB_PDB_CHAIN_ID,n,RC );
    if ((!RC) && F)  {
      strcpy_n0 ( chainID,F,sizeof(ChainID)-1 );
      Loop->DeleteField ( CIFTAG_NDB_PDB_CHAIN_ID,n );
    } else
      strcpy ( chainID,"" );

    CIFGetString ( resName,Loop,CIFTAG_MON_ID,n,sizeof(ResName),
                   pstr("UNK") );

    CIFGetIntegerD ( seqNum,Loop,CIFTAG_SEQ_NUM );

    CIFGetString ( insCode,Loop,CIFTAG_NDB_PDB_INS_CODE,
                   n,sizeof(InsCode),pstr(" ") );

    CIFGetString ( database,Loop,CIFTAG_NDB_SEQ_DB_NAME,n,
                   sizeof(DBName),pstr(" ") );

    CIFGetString ( dbAccession,Loop,CIFTAG_NDB_SEQ_DB_ACCESSION_CODE,
                   n,sizeof(DBAcCode),pstr(" ") );

    CIFGetString ( dbRes,Loop,CIFTAG_DB_MON_ID,n,sizeof(ResName),
                   pstr("   ") );

    CIFGetIntegerD ( dbSeq,Loop,CIFTAG_NDB_SEQ_DB_SEQ_NUM );
  //  if (CIFGetInteger1(dbSeq,Loop,CIFTAG_NDB_SEQ_DB_SEQ_NUM,n))
  //    dbSeq = MinInt4;

    F = Loop->GetString ( CIFTAG_DETAILS,n,RC );
    if ((!RC) && F)  {
      CreateCopy ( conflict,F );
      Loop->DeleteField ( CIFTAG_DETAILS,n );
    } else
      CreateCopy ( conflict,pstr(" ") );

    n++;

    return Error_NoError;

  }

  void  SeqAdv::Copy ( PContainerClass SeqAdv )  {

    ContainerClass::Copy ( SeqAdv );

    seqNum = PSeqAdv(SeqAdv)->seqNum;
    dbSeq  = PSeqAdv(SeqAdv)->dbSeq;
    strcpy  ( resName    ,PSeqAdv(SeqAdv)->resName     );
    strcpy  ( insCode    ,PSeqAdv(SeqAdv)->insCode     );
    strcpy  ( database   ,PSeqAdv(SeqAdv)->database    );
    strcpy  ( dbAccession,PSeqAdv(SeqAdv)->dbAccession );
    strcpy  ( dbRes      ,PSeqAdv(SeqAdv)->dbRes       );
    CreateCopy ( conflict,PSeqAdv(SeqAdv)->conflict    );

  }

  void  SeqAdv::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte    ( &Version );
    f.WriteInt     ( &seqNum  );
    f.WriteInt     ( &dbSeq   );
    f.WriteTerLine ( resName    ,false );
    f.WriteTerLine ( insCode    ,false );
    f.WriteTerLine ( database   ,false );
    f.WriteTerLine ( dbAccession,false );
    f.WriteTerLine ( dbRes      ,false );
    f.CreateWrite  ( conflict );
  }

  void  SeqAdv::read  ( io::RFile f ) {
  byte Version;
    f.ReadByte    ( &Version );
    f.ReadInt     ( &seqNum  );
    f.ReadInt     ( &dbSeq   );
    f.ReadTerLine ( resName    ,false );
    f.ReadTerLine ( insCode    ,false );
    f.ReadTerLine ( database   ,false );
    f.ReadTerLine ( dbAccession,false );
    f.ReadTerLine ( dbRes      ,false );
    f.CreateRead  ( conflict );
  }

  MakeStreamFunctions(SeqAdv)



  //  ================  SeqRes  ===================

  SeqRes::SeqRes() : io::Stream()  {
    InitSeqRes();
  }

  SeqRes::SeqRes ( io::RPStream Object ) : io::Stream(Object)  {
    InitSeqRes();
  }

  SeqRes::~SeqRes()  {
    FreeMemory();
  }

  void  SeqRes::SetChain ( PChain Chain_Owner )  {
    chain = Chain_Owner;
    if (chain)  strcpy ( chainID,chain->chainID );
          else  strcpy ( chainID,"" );
  }

  void  SeqRes::InitSeqRes()  {
    chain   = NULL;
    numRes  = -1;
    resName = NULL;
    serNum  = 0;
    strcpy ( chainID,"" );
  }

  void  SeqRes::FreeMemory()  {
    if (resName)  delete[] resName;
    resName = NULL;
    numRes  = -1;
    serNum  = 0;
  }

  void  SeqRes::PDBASCIIDump ( io::RFile f )  {
  //  writes the ASCII PDB SEQRES lines into file f
  char S[100];
  int  i,k,sN;
    if (numRes<0)  return;
    strcpy     ( S,"SEQRES" );
    PadSpaces  ( S,80 );
    if (chain->chainID[0])
      S[11] = chain->chainID[0];
    PutInteger ( &(S[13]),numRes,4 );
    if (resName)  {
      i  = 0;
      sN = 1;
      while (i<numRes)  {
        PutInteger ( &(S[7]),sN,3 );
        k = 19;
        while ((i<numRes) && (k<70))  {
          if (resName[i][0])
                strcpy_n ( &(S[k]),resName[i],3 );
          else  strcpy_n ( &(S[k]),pstr("   "),3 );
          i++;
          k += 4;
        }
        while (k<70)  {
          strcpy_n ( &(S[k]),pstr("   "),3 );
          k += 4;
        }
        f.WriteLine ( S );
        sN++;
      }
    } else  {
      S[9] = '0';
      strcpy_n ( &(S[19]),pstr("UNK"),3 );
      f.WriteLine ( S );
    }
  }

  ERROR_CODE SeqRes::ConvertPDBASCII ( cpstr S )  {
  int i,k,sN,nR;
    if (chain->chainID[0])  {
      if (S[11]!=chain->chainID[0])
        return Error_WrongChainID;
    } else if (S[11]!=' ')  {
      chain->chainID[0] = S[11];
      chain->chainID[1] = char(0);
    } else
      chain->chainID[0] = char(0);
    GetInteger ( sN,&(S[8]) ,3 );
    GetInteger ( nR,&(S[13]),4 );
    if (sN==0)  {
      FreeMemory();
      numRes = nR;
    } else  {
      serNum++;
      if (sN!=serNum)
        return Error_SEQRES_serNum;
      if (sN==1)  {
        FreeMemory();
        resName = new ResName[nR];
        for (i=0;i<nR;i++)
          resName[i][0] = char(0);
        numRes  = nR;
        serNum  = sN;
      } else if (nR!=numRes)
        return Error_SEQRES_numRes;
      i = 0;
      while ((i<nR) && (resName[i][0]))  i++;
      if (i>=nR)
        return Error_SEQRES_extraRes;
      k = 19;
      while ((i<nR) && (k<70))  {
        GetString ( resName[i],&(S[k]),3 );
        if (!strcmp(resName[i],"   "))  resName[i][0] = char(0);
                                  else  i++;
        k += 4;
      }
    }
    return Error_NoError;
  }


  void  SeqRes::MakeCIF ( mmcif::PData CIF )  {
  //  Note that SeqRes only adds sequence to the CIF loop common
  // to all chains. Therefore this loop should be wiped off from
  // CIF structure before putting first sequence into it.
  mmcif::PLoop Loop;
  int          RC,i;

    if (numRes<0)  return;

    RC = CIF->AddLoop ( CIFCAT_NDB_POLY_SEQ_SCHEME,Loop );
    if (RC!=mmcif::CIFRC_Ok)  {
      // the category was (re)created, provide tags
      Loop->AddLoopTag ( CIFTAG_ID     );
      Loop->AddLoopTag ( CIFTAG_MON_ID );
    }

    if (resName)
      for (i=0;i<numRes;i++)  {
        Loop->AddString  ( chain->chainID,true );
        Loop->AddString  ( resName[i]    ,true );
      }
    else
      for (i=0;i<numRes;i++)  {
        Loop->AddString  ( chain->GetEntryID(),true );
        Loop->AddString  ( pstr("UNK")        ,true );
      }

  }

  ERROR_CODE SeqRes::GetCIF ( mmcif::PData CIF )  {
  //   Tries to get sequence from the CIF structure. A sequence
  // for first met chain is extracted and then removed from
  // the CIF structure, so that sequential calls will extract
  // all sequencies. Chain ID is stored locally in chainID;
  // reference to parent chain is neither used nor checked.
  //   Returns 0 if sequence was extracted and 1 otherwise.
  mmcif::PLoop Loop;
  ResName    * rN;
  ChainID      chID;
  pstr         F;
  cpstr        CHAIN_ID;
  int          RC,i,l;
  CIF_MODE     CIFMode;
  bool         isMon;

    FreeMemory();

    CIFMode = CIF_NDB;
    Loop = CIF->GetLoop ( CIFName(CAT_POLY_SEQ_SCHEME,CIFMode) );
    if (!Loop)  {
      CIFMode = CIF_PDBX;
      Loop = CIF->GetLoop ( CIFName(CAT_POLY_SEQ_SCHEME,CIFMode) );
      if (!Loop)  return Error_NoLoop;
    }

    l = Loop->GetLoopLength();
    if (l<=0)  return Error_NoLoop;

    rN         = new ResName[l];
    chainID[0] = char(1);
    numRes     = 0;
    isMon      = false;
    CHAIN_ID   = CIFName(TAG_SEQ_CHAIN_ID,CIFMode);
    for (i=0;i<l;i++)  {
      F = Loop->GetString ( CHAIN_ID,i,RC );
      if (!RC)  {
        if (F)  strcpy ( chID,F );
          else  chID[0] = char(0);
        if (chainID[0]==char(1))  strcpy ( chainID,chID );
        if (!strcmp(chainID,chID))  {
          CIFGetString ( rN[numRes],Loop,CIFTAG_MON_ID,i,
                         sizeof(ResName),pstr("UNK") );
          Loop->DeleteField ( CHAIN_ID,i );
          if (strcmp(rN[numRes],"UNK")) isMon = true;
          numRes++;
        }
      }
    }

    if (numRes==0)  {
      numRes = -1;
      delete[] rN;
      return Error_EmptyCIFLoop;
    }

    if (isMon)  {
      resName = new ResName[numRes];
      for (i=0;i<numRes;i++)
        strcpy ( resName[i],rN[i] );
    }

    delete[] rN;

    return Error_NoError;

  }

  void  SeqRes::Copy ( PSeqRes SeqRes )  {
  int i;

    FreeMemory();

    numRes = SeqRes->numRes;
    serNum = SeqRes->serNum;

    if (SeqRes->resName)  {
      resName = new ResName[numRes];
      for (i=0;i<numRes;i++)
        strcpy ( resName[i],SeqRes->resName[i] );
    }

  }

  void  SeqRes::write ( io::RFile f )  {
  int  i;
  byte Version=1;
    f.WriteByte ( &Version );
    f.WriteInt  ( &numRes  );
    f.WriteInt  ( &serNum  );
    if (resName)  i = 1;
            else  i = 0;
    f.WriteInt ( &i );
    if (resName)
      for (i=0;i<numRes;i++)
        f.WriteTerLine ( resName[i],false );
  }

  void  SeqRes::read  ( io::RFile f ) {
  int  i;
  byte Version;
    FreeMemory();
    f.ReadByte ( &Version );
    f.ReadInt  ( &numRes  );
    f.ReadInt  ( &serNum  );
    f.ReadInt  ( &i       );
    if (i)  {
      resName = new ResName[numRes];
      for (i=0;i<numRes;i++)
        f.ReadTerLine ( resName[i],false );
    }
  }


  MakeStreamFunctions(SeqRes)



  //  ================  ModRes  ===================

  ModRes::ModRes() : ContainerChain()  {
    InitModRes();
  }

  ModRes::ModRes ( PChain Chain_Owner )
        : ContainerChain(Chain_Owner)  {
    InitModRes();
  }

  ModRes::ModRes ( PChain Chain_Owner, cpstr S )
        : ContainerChain(Chain_Owner)  {
    InitModRes();
    ConvertPDBASCII ( S );
  }

  ModRes::ModRes ( io::RPStream Object ) : ContainerChain(Object)  {
    InitModRes();
  }

  ModRes::~ModRes()  {
    if (comment)  delete[] comment;
  }

  void  ModRes::InitModRes()  {
    strcpy     ( resName,"---" );
    seqNum  = 0;
    strcpy     ( insCode,"-"   );
    comment = NULL;
    CreateCopy ( comment,pstr(" ") );
    strcpy     ( stdRes ,"---" );
  }

  void  ModRes::PDBASCIIDump ( pstr S, int N )  {
  UNUSED_ARGUMENT(N);
  //  makes the ASCII PDB MODRES line number N
  //  from the class' data
    strcpy     ( S,"MODRES" );
    PadSpaces  ( S,80 );
    strcpy_n   ( &(S[7]) ,chain->GetEntryID(),4  );
    strcpy_n   ( &(S[12]),resName      ,3  );
    if (chain->chainID[0])  S[16] = chain->chainID[0];
    PutIntIns  ( &(S[18]),seqNum,4,insCode );
    strcpy_n   ( &(S[24]),stdRes       ,3  );
    strcpy_n   ( &(S[29]),comment,IMin(strlen(comment),41) );
  }

  ERROR_CODE ModRes::ConvertPDBASCII ( cpstr S )  {
  IDCode idCode;
    if (chain->chainID[0])  {
      if (S[16]!=chain->chainID[0])
        return Error_WrongChainID;
    } else if (S[16]!=' ')  {
      chain->chainID[0] = S[16];
      chain->chainID[1] = char(0);
    } else
      chain->chainID[0] = char(0);
    strcpy ( idCode,chain->GetEntryID() );
    if (idCode[0])  {
      if (strncmp(&(S[7]),idCode,4) && (!ignoreNonCoorPDBErrors))
        return Error_WrongEntryID;
    } else  {
      GetString ( idCode,&(S[7]),4 );
      chain->SetEntryID ( idCode );
    }
    GetString  ( resName       ,&(S[12]),3 );
    GetIntIns  ( seqNum,insCode,&(S[18]),4 );
    GetString  ( stdRes        ,&(S[24]),3 );
    CreateCopy ( comment       ,&(S[29])   );
    CutSpaces  ( comment,SCUTKEY_END       );
    return Error_NoError;
  }

  void  ModRes::MakeCIF ( mmcif::PData CIF, int N )  {
  UNUSED_ARGUMENT(CIF);
  UNUSED_ARGUMENT(N);
  /*  -- apparently wrong use of _struct_conn, to be revised
  mmcif::PLoop Loop;
  int         RC;

    RC = CIF->AddLoop ( CIFCAT_STRUCT_CONN,Loop );

    if (RC!=mmcif::CIFRC_Ok)  {
      // the category was (re)created, provide tags
      Loop->AddLoopTag ( CIFTAG_CONN_TYPE_ID               );
      Loop->AddLoopTag ( CIFTAG_NDB_PDB_ID                 );
      Loop->AddLoopTag ( CIFTAG_PTNR1_LABEL_COMP_ID        );
      Loop->AddLoopTag ( CIFTAG_PTNR1_LABEL_ASYM_ID        );
      Loop->AddLoopTag ( CIFTAG_PTNR1_LABEL_SEQ_ID         );
      Loop->AddLoopTag ( CIFTAG_NDB_PTNR1_LABEL_INS_CODE   );
      Loop->AddLoopTag ( CIFTAG_NDB_PTNR1_STANDARD_COMP_ID );
      Loop->AddLoopTag ( CIFTAG_DETAILS                    );
    }

    Loop->AddString  ( pstr("MODRES")           );
    Loop->AddString  ( chain->GetEntryID(),true );
    Loop->AddString  ( resName            ,true );
    Loop->AddString  ( chain->chainID     ,true );
    Loop->AddInteger ( seqNum                   );
    Loop->AddString  ( insCode            ,true );
    Loop->AddString  ( stdRes             ,true );
    Loop->AddString  ( comment            ,true );

  */

  }

  ERROR_CODE ModRes::GetCIF ( mmcif::PData CIF, int & n )  {
  UNUSED_ARGUMENT(CIF);
  //  GetCIF(..) must be always run without reference to Chain,
  //  see CModel::GetCIF(..).

  /*  -- apparently wrong use of _struct_conn, to be revised
  mmcif::PLoop   Loop;
  pstr          F;
  int           l,RC;

    Loop = CIF->GetLoop ( CIFCAT_STRUCT_CONN );
    if (!Loop)  {
      n = -1;
      return;
    }

    l = Loop->GetLoopLength();
    while (n<l)  {
      F = Loop->GetString ( CIFTAG_CONN_TYPE_ID,n,RC );
      if ((!RC) && F)  {
        if (!strcmp(F,"MODRES"))  break;
      }
      n++;
    }
    if (n>=l)  {
      n = -1;
      return;
    }

    Loop->DeleteField ( CIFTAG_CONN_TYPE_ID,n );

    //  Determine the ChainID first and store it locally. It will
    // be used by CModel for generating chains and placing the
    // primary structure data BEFORE reading the coordinate section.
    F = Loop->GetString ( CIFTAG_PTNR1_LABEL_ASYM_ID,n,RC );
    if ((!RC) && F)  {
      strcpy_n0 ( chainID,F,sizeof(ChainID)-1 );
      Loop->DeleteField ( CIFTAG_PTNR1_LABEL_ASYM_ID,n );
    } else
      strcpy ( chainID,"" );


    CIFGetString ( resName,Loop,CIFTAG_PTNR1_LABEL_COMP_ID,n,
                   sizeof(ResName),pstr("UNK") );

    if (CIFGetInteger(seqNum,Loop,CIFTAG_PTNR1_LABEL_SEQ_ID,n))
      return;

    CIFGetString ( insCode,Loop,CIFTAG_NDB_PTNR1_LABEL_INS_CODE,
                   n,sizeof(InsCode),pstr(" ") );

    CIFGetString ( stdRes,Loop,CIFTAG_NDB_PTNR1_STANDARD_COMP_ID,n,
                   sizeof(ResName),pstr("UNK") );

    F = Loop->GetString ( CIFTAG_DETAILS,n,RC );
    if ((!RC) && F)  {
      CreateCopy ( comment,F );
      Loop->DeleteField ( CIFTAG_DETAILS,n );
    } else
      CreateCopy ( comment,pstr(" ") );

    n++;

  */

    n = -1;

    return Error_EmptyCIF;

  }

  void  ModRes::Copy ( PContainerClass ModRes )  {
    seqNum = PModRes(ModRes)->seqNum;
    strcpy ( resName,PModRes(ModRes)->resName );
    strcpy ( insCode,PModRes(ModRes)->insCode );
    strcpy ( stdRes ,PModRes(ModRes)->stdRes  );
    CreateCopy ( comment,PModRes(ModRes)->comment );
  }

  void  ModRes::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte    ( &Version );
    f.WriteInt     ( &seqNum  );
    f.WriteTerLine ( resName,false );
    f.WriteTerLine ( insCode,false );
    f.WriteTerLine ( stdRes ,false );
    f.CreateWrite  ( comment  );
  }

  void  ModRes::read  ( io::RFile f ) {
  byte Version;
    f.ReadByte    ( &Version );
    f.ReadInt     ( &seqNum  );
    f.ReadTerLine ( resName,false );
    f.ReadTerLine ( insCode,false );
    f.ReadTerLine ( stdRes ,false );
    f.CreateRead  ( comment  );
  }

  MakeStreamFunctions(ModRes)



  //  ================  HetRec  ======================

  HetRec::HetRec() : ContainerChain()  {
    InitHetRec();
  }

  HetRec::HetRec ( PChain Chain_Owner )
        : ContainerChain(Chain_Owner)  {
    InitHetRec();
  }

  HetRec::HetRec ( PChain Chain_Owner, cpstr S )
        : ContainerChain(Chain_Owner)  {
    InitHetRec();
    ConvertPDBASCII ( S );
  }

  HetRec::HetRec ( io::RPStream Object ) : ContainerChain(Object)  {
    InitHetRec();
  }

  HetRec::~HetRec()  {
    if (comment)  delete[] comment;
  }

  void  HetRec::InitHetRec()  {
    strcpy ( hetID  ,"---" );
    strcpy ( insCode,"-"   );
    seqNum      = 0;
    numHetAtoms = 0;
    comment     = NULL;
    CreateCopy ( comment,pstr(" ") );
  }

  void  HetRec::PDBASCIIDump ( pstr S, int N )  {
  UNUSED_ARGUMENT(N);
  //  makes the ASCII PDB MODRES line number N
  //  from the class' data
    strcpy     ( S,"HET" );
    PadSpaces  ( S,80 );
    strcpy_n   ( &(S[7]) ,hetID,3  );
    if (chain->chainID[0])  S[12] = chain->chainID[0];
    PutIntIns  ( &(S[13]),seqNum,4,insCode );
    PutInteger ( &(S[20]),numHetAtoms,5    );
    strcpy_n   ( &(S[30]),comment,IMin(strlen(comment),40) );
  }

  ERROR_CODE HetRec::ConvertPDBASCII ( cpstr S )  {
    if (chain->chainID[0])  {
      if (S[12]!=chain->chainID[0])
        return Error_WrongChainID;
    } else if (S[12]!=' ')  {
      chain->chainID[0] = S[12];
      chain->chainID[1] = char(0);
    } else
      chain->chainID[0] = char(0);
    GetString  ( hetID         ,&(S[7]) ,3 );
    GetIntIns  ( seqNum,insCode,&(S[13]),4 );
    GetInteger ( numHetAtoms   ,&(S[20]),5 );
    CreateCopy ( comment       ,&(S[30])   );
    CutSpaces  ( comment,SCUTKEY_END       );
    return Error_NoError;
  }

  void  HetRec::MakeCIF ( mmcif::PData CIF, int N )  {
  UNUSED_ARGUMENT(N);
  mmcif::PLoop Loop;
  int         RC;

    RC = CIF->AddLoop ( CIFCAT_NDB_NONSTANDARD_LIST,Loop );

    if (RC!=mmcif::CIFRC_Ok)  {
      // the category was (re)created, provide tags
      Loop->AddLoopTag ( CIFTAG_ID              );
      Loop->AddLoopTag ( CIFTAG_AUTH_ASYM_ID    );
      Loop->AddLoopTag ( CIFTAG_AUTH_SEQ_ID     );
      Loop->AddLoopTag ( CIFTAG_INS_CODE        );
      Loop->AddLoopTag ( CIFTAG_NUMBER_ATOMS_NH );
      Loop->AddLoopTag ( CIFTAG_DETAILS         );
    }

    Loop->AddString  ( hetID         ,true );
    Loop->AddString  ( chain->chainID,true );
    Loop->AddInteger ( seqNum              );
    Loop->AddString  ( insCode       ,true );
    Loop->AddInteger ( numHetAtoms         );
    Loop->AddString  ( comment       ,true );

  }

  ERROR_CODE HetRec::GetCIF ( mmcif::PData CIF, int & n )  {
  //  GetCIF(..) must be always run without reference to Chain,
  //  see CModel::GetCIF(..).
  mmcif::PLoop   Loop;
  pstr           F;
  int            RC;
  ERROR_CODE     rc;

    Loop = CIF->GetLoop ( CIFCAT_NDB_NONSTANDARD_LIST );
    if (!Loop)  {
      n = -1;
      return Error_EmptyCIF;
    }

    if (n>=Loop->GetLoopLength())  {
      n = -1;
      return Error_EmptyCIF;
    }

    //  Determine the ChainID first and store it locally. It will
    // be used by CModel for generating chains and placing the
    // primary structure data BEFORE reading the coordinate section.
    F = Loop->GetString ( CIFTAG_AUTH_ASYM_ID,n,RC );
    if ((!RC) && F)  {
      strcpy_n0 ( chainID,F,sizeof(ChainID)-1 );
      Loop->DeleteField ( CIFTAG_AUTH_ASYM_ID,n );
    } else
      strcpy ( chainID,"" );


    CIFGetString ( hetID,Loop,CIFTAG_ID,n,sizeof(ResName),
                   pstr("UNK") );

    rc = CIFGetInteger ( seqNum,Loop,CIFTAG_AUTH_SEQ_ID,n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    CIFGetString ( insCode,Loop,CIFTAG_INS_CODE,n,sizeof(InsCode),
                   pstr(" ") );

    rc = CIFGetInteger ( numHetAtoms,Loop,CIFTAG_NUMBER_ATOMS_NH,n );
    if (rc==Error_NoData)   return Error_EmptyCIF;
    if (rc!=Error_NoError)  return rc;

    F = Loop->GetString ( CIFTAG_DETAILS,n,RC );
    if ((!RC) && F)  {
      CreateCopy ( comment,F );
      Loop->DeleteField ( CIFTAG_DETAILS,n );
    } else
      CreateCopy ( comment,pstr(" ") );

    n++;

    return Error_NoError;

  }

  void  HetRec::Copy ( PContainerClass Het )  {
    seqNum      = PHetRec(Het)->seqNum;
    numHetAtoms = PHetRec(Het)->numHetAtoms;
    strcpy     ( hetID  ,PHetRec(Het)->hetID   );
    strcpy     ( insCode,PHetRec(Het)->insCode );
    CreateCopy ( comment,PHetRec(Het)->comment );
  }

  void  HetRec::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte    ( &Version      );
    f.WriteInt     ( &seqNum       );
    f.WriteInt     ( &numHetAtoms  );
    f.WriteTerLine ( hetID  ,false );
    f.WriteTerLine ( insCode,false );
    f.CreateWrite  ( comment       );
  }

  void  HetRec::read  ( io::RFile f ) {
  byte Version;
    f.ReadByte    ( &Version      );
    f.ReadInt     ( &seqNum       );
    f.ReadInt     ( &numHetAtoms  );
    f.ReadTerLine ( hetID  ,false );
    f.ReadTerLine ( insCode,false );
    f.CreateRead  ( comment       );
  }

  MakeStreamFunctions(HetRec)



  //  =====================   Chain   =======================

  Chain::Chain() : UDData() {
    InitChain();
    SetChain ( pstr("") );
  }

  Chain::Chain ( PProModel Model, const ChainID chID ) : UDData()  {
    InitChain();
    SetChain ( chID );
    if (Model)  Model->AddChain ( this );
  }

  Chain::Chain ( io::RPStream Object ) : UDData(Object)  {
    InitChain();
    SetChain ( pstr("") );
  }

  void  Chain::InitChain()  {
    nResidues      = 0;
    resLen         = 0;
    residue        = NULL;
    model          = NULL;
    chainID[0]     = char(0);
    prevChainID[0] = char(0);
    nWeights       = 0;
    Weight         = 0.0;
    Exclude        = true;
  }

  void  Chain::SetChain ( const ChainID chID )  {
    strcpy ( chainID,chID );
    if (chID[0]==' ')  chainID[0] = char(0);
    DBRef .SetChain ( this );
    seqAdv.SetChain ( this );
    seqRes.SetChain ( this );
    modRes.SetChain ( this );
    Het   .SetChain ( this );
  }

  void  Chain::SetChainID ( const ChainID chID )  {
    strcpy ( chainID,chID );
    if (chID[0]==' ')  chainID[0] = char(0);
  }

  Chain::~Chain()  {
    FreeMemory();
    if (model)  model->_ExcludeChain ( chainID );
  }

  void  Chain::FreeMemory()  {
    DeleteAllResidues();
    if (residue)  delete[] residue;
    resLen    = 0;
    nResidues = 0;
    residue   = NULL;
    FreeAnnotations();
  }

  void  Chain::FreeAnnotations()  {
    DBRef .FreeContainer();
    seqAdv.FreeContainer();
    seqRes.FreeMemory   ();
    modRes.FreeContainer();
    Het   .FreeContainer();
  }

  void Chain::SetModel ( PProModel Model )  {
    model = Model;
  }

  PManager Chain::GetCoordHierarchy()  {
    if (model)  return model->GetCoordHierarchy();
    return NULL;
  }

  void Chain::CheckInAtoms()  {
  int i;
    if (GetCoordHierarchy())
      for (i=0;i<nResidues;i++)
        if (residue[i])
          residue[i]->CheckInAtoms();
  }

  ERROR_CODE Chain::ConvertDBREF ( cpstr PDBString ) {
  PContainerChain ContainerChain;
  ERROR_CODE      RC;
    ContainerChain = new DBReference(this);
    RC = ContainerChain->ConvertPDBASCII ( PDBString );
    if (RC)  {
      delete ContainerChain;
      return RC;
    }
    DBRef.AddData ( ContainerChain );
    return Error_NoError;
  }

  ERROR_CODE Chain::ConvertSEQADV ( cpstr PDBString ) {
  PContainerChain ContainerChain;
  ERROR_CODE      RC;
    ContainerChain = new SeqAdv(this);
    RC = ContainerChain->ConvertPDBASCII ( PDBString );
    if (RC)  {
      delete ContainerChain;
      return RC;
    }
    seqAdv.AddData ( ContainerChain );
    return Error_NoError;
  }

  ERROR_CODE Chain::ConvertSEQRES ( cpstr PDBString ) {
    return seqRes.ConvertPDBASCII ( PDBString );
  }

  ERROR_CODE Chain::ConvertMODRES ( cpstr PDBString ) {
  PContainerChain ContainerChain;
  ERROR_CODE      RC;
    ContainerChain = new ModRes(this);
    RC = ContainerChain->ConvertPDBASCII ( PDBString );
    if (RC)  {
      delete ContainerChain;
      return RC;
    }
    modRes.AddData ( ContainerChain );
    return Error_NoError;
  }

  ERROR_CODE  Chain::ConvertHET ( cpstr PDBString ) {
  PContainerChain ContainerChain;
  ERROR_CODE      RC;
    ContainerChain = new HetRec(this);
    RC = ContainerChain->ConvertPDBASCII ( PDBString );
    if (RC)  {
      delete ContainerChain;
      return RC;
    }
    Het.AddData ( ContainerChain );
    return Error_NoError;
  }


  void  Chain::PDBASCIIDump ( io::RFile f )  {
  // this function was for test purposes and is not used
  // for normal function of MMDB
    DBRef .PDBASCIIDump ( f );
    seqAdv.PDBASCIIDump ( f );
    seqRes.PDBASCIIDump ( f );
    modRes.PDBASCIIDump ( f );
    Het   .PDBASCIIDump ( f );
  }

  void  Chain::PDBASCIIAtomDump ( io::RFile f )  {
  int i;
    for (i=0;i<nResidues;i++)
      if (residue[i])
        residue[i]->PDBASCIIAtomDump ( f );
  }

  void  Chain::MakeAtomCIF ( mmcif::PData CIF )  {
  int i;
    for (i=0;i<nResidues;i++)
      if (residue[i])
        residue[i]->MakeAtomCIF ( CIF );
  }


  int  Chain::GetNumberOfResidues()  {
    return nResidues;
  }

  PResidue Chain::GetResidue ( int resNo )  {
    if ((0<=resNo) && (resNo<nResidues))
          return residue[resNo];
    else  return NULL;
  }


  PResidue Chain::GetResidueCreate ( const ResName resName,
                                      int           seqNum,
                                      const InsCode insCode,
                                      bool          Enforce )  {
  //   Returns pointer on residue, whose name, sequence number and
  // insert code are given in resName, seqNum and insCode, respectively.
  // If such a residue is absent in the chain, one is created at
  // the end of the chain.
  int i;

    // check if such a residue is already in the chain
    if (insCode[0])  {
      for (i=0;i<nResidues;i++)
        if (residue[i])  {
          if ((seqNum==residue[i]->seqNum) &&
              (!strcmp(insCode,residue[i]->insCode)))  {
            if (!strcmp(resName,residue[i]->name))
              return residue[i]; // it is there; just return the pointer
            else if (!Enforce)
              return NULL;       // duplicate seqNum and insCode!
          }
        }
    } else  {
      for (i=0;i<nResidues;i++)
        if (residue[i])  {
          if ((seqNum==residue[i]->seqNum) &&
              (!residue[i]->insCode[0]))  {
            if (!strcmp(resName,residue[i]->name))
              return residue[i]; // it is there; just return the pointer
            else if (!Enforce)
              return NULL;       // duplicate seqNum and insCode!
          }
        }
    }

    // expand the residue array, if necessary
    if (nResidues>=resLen)
      ExpandResidueArray ( 100 );

    // create new residue
    residue[nResidues] = newResidue();
    residue[nResidues]->SetChain ( this );
    residue[nResidues]->SetResID ( resName,seqNum,insCode );
    residue[nResidues]->index = nResidues;
    nResidues++;

    return residue[nResidues-1];

  }

  void  Chain::ExpandResidueArray ( int inc )  {
  PPResidue Residue1;
  int        i;
    resLen  += inc;
    Residue1 = new PResidue[resLen];
    for (i=0;i<nResidues;i++)
      Residue1[i] = residue[i];
    if (residue) delete[] residue;
    residue = Residue1;
    for (i=nResidues;i<resLen;i++)
      residue[i] = NULL;
  }

  PResidue Chain::GetResidue ( int seqNum, const InsCode insCode )  {
  //   Returns pointer on residue, whose sequence number and
  // insert code are given in seqNum and insCode, respectively.
  // If such a residue is absent in the chain, returns NULL.
  int     i;
  bool isInsCode;
    if (insCode)  isInsCode = insCode[0]!=char(0);
            else  isInsCode = false;
    if (isInsCode)  {
      for (i=0;i<nResidues;i++)
        if (residue[i])  {
          if ((seqNum==residue[i]->seqNum) &&
              (!strcmp(insCode,residue[i]->insCode)))
            return residue[i];
        }
    } else  {
      for (i=0;i<nResidues;i++)
        if (residue[i])  {
          if ((seqNum==residue[i]->seqNum) && (!residue[i]->insCode[0]))
            return residue[i];
        }
    }
    return NULL;
  }

  int  Chain::GetResidueNo ( int seqNum, const InsCode insCode )  {
  //   GetResidueNo(..) returns the residue number in the chain's
  // residues table. Residues are numbered as 0..nres-1 as they appear
  // in the coordinate file.
  //   If residue is not found, the function returns -1.
  int      i;
  bool isInsCode;
    if (insCode)  isInsCode = insCode[0]!=char(0);
            else  isInsCode = false;
    if (isInsCode)  {
      for (i=0;i<nResidues;i++)
        if (residue[i])  {
          if ((seqNum==residue[i]->seqNum) &&
              (!strcmp(insCode,residue[i]->insCode)))
            return i;
        }
    } else  {
      for (i=0;i<nResidues;i++)
        if (residue[i])  {
          if ((seqNum==residue[i]->seqNum) && (!residue[i]->insCode[0]))
            return i;
        }
    }
    return -1;
  }

  void Chain::GetResidueTable ( PPResidue & resTable,
                                int & NumberOfResidues )  {
    resTable         = residue;
    NumberOfResidues = nResidues;
  }

  int  Chain::_ExcludeResidue ( const ResName resName, int seqNum,
                                const InsCode insCode )  {
  //   ExcludeResidue(..) excludes (but does not dispose!) a residue
  // from the chain. Returns 1 if the chain gets empty and 0 otherwise.
  int  i,k;

    if (!Exclude)  return 0;

    // find the residue
    k = -1;
    for (i=0;(i<nResidues) && (k<0);i++)
      if ((seqNum==residue[i]->seqNum)           &&
          (!strcmp(insCode,residue[i]->insCode)) &&
          (!strcmp(resName,residue[i]->name)))
        k = i;

    if (k>=0)  {
      for (i=k+1;i<nResidues;i++)  {
        residue[i-1] = residue[i];
        if (residue[i-1])
          residue[i-1]->index = i-1;
      }
      nResidues--;
      residue[nResidues] = NULL;
    }

    if (nResidues<=0)  return 1;
                 else  return 0;

  }


  void Chain::GetCoordSequence ( pstr & seq )  {
  //   GetCoorSequence(...) returns sequence inferred from list
  // of residues (which may differ from one in the file header).
  // The sequence is returned as a null-terminated string 'seq'.
  // On input, 'seq' should be either NULL or allocated (in which
  // case the original allocation will be released).
  int i,j;

    if (seq) delete[] seq;
    seq = new char[nResidues+1];

    j = 0;
    for (i=0;i<nResidues;i++)
      if (residue[i])
        Get1LetterCode ( residue[i]->GetResName(),seq[j++] );
    seq[j] = char(0);

  }


  //  ------------------  Deleting residues  --------------------------

  int  Chain::DeleteResidue ( int resNo )  {
    if ((0<=resNo) && (resNo<nResidues))  {
      if (residue[resNo])  {
        Exclude = false;
        delete residue[resNo];
        residue[resNo] = NULL;
        Exclude = true;
        return 1;
      }
    }
    return 0;
  }

  int  Chain::DeleteResidue ( int seqNum, const InsCode insCode )  {
  int i;
    if (insCode[0])  {
      for (i=0;i<nResidues;i++)
        if (residue[i])  {
          if ((seqNum==residue[i]->seqNum) &&
              (!strcmp(insCode,residue[i]->insCode)))  {
            Exclude = false;
            delete residue[i];
            residue[i] = NULL;
            Exclude = true;
            return 1;
          }
        }
    } else  {
      for (i=0;i<nResidues;i++)
        if (residue[i])  {
          if ((seqNum==residue[i]->seqNum) && (!residue[i]->insCode[0]))  {
            Exclude = false;
            delete residue[i];
            residue[i] = NULL;
            Exclude = true;
            return 1;
          }
        }
    }
    return 0;
  }


  int  Chain::DeleteAllResidues()  {
  int i,k;
    Exclude = false;
    k = 0;
    for (i=0;i<nResidues;i++)
      if (residue[i])  {
        delete residue[i];
        residue[i] = NULL;
        k++;
      }
    nResidues = 0;
    Exclude = true;
    return k;
  }


  int  Chain::DeleteSolvent()  {
  int i,k;
    Exclude = false;
    k = 0;
    for (i=0;i<nResidues;i++)
      if (residue[i])  {
        if (residue[i]->isSolvent())  {
          delete residue[i];
          residue[i] = NULL;
          k++;
        }
      }
    Exclude = true;
    return k;
  }


  void Chain::TrimResidueTable()  {
  int i,j;
    Exclude = false;
    j = 0;
    for (i=0;i<nResidues;i++)
      if (residue[i])  {
        if (residue[i]->nAtoms>0)  {
          if (j<i)  {
            residue[j] = residue[i];
            residue[j]->index = j;
            residue[i] = NULL;
          }
          j++;
        } else  {
          delete residue[i];
          residue[i] = NULL;
        }
      }
    nResidues = j;
    Exclude   = true;
  }

  int  Chain::AddResidue ( PResidue res )  {
  //  modify both CModel::Copy methods simultaneously!
  //
  //  Copy(PCModel,PPAtom,int&) copies atoms into array 'atom'
  // starting from position atom_index. 'atom' should be able to
  // accept all new atoms - no checks on the length of 'atom'
  // is being made. This function should not be used in applications.
    return InsResidue ( res,nResidues );
  }

  /*
  PCmmdbRoot mmdbRoot;
  PChain    chain1;
  int        i;

    for (i=0;i<nResidues;i++)
      if (residue[i]==res)  return -i;  // this residue is already there

    if (res)  {

      mmdbRoot = PCmmdbRoot(GetCoordHierarchy());

      // get space for new residue
      if (nResidues>=resLen)
        ExpandResidueArray ( 100 );

      if (res->GetCoordHierarchy())  {
        residue[nResidues] = newResidue();
        residue[nResidues]->SetChain ( this );
        residue[nResidues]->SetResID ( res->name,res->seqNum,res->insCode );
        if (mmdbRoot)  {
          // get space for new atoms
          mmdbRoot->AddAtomArray ( res->GetNumberOfAtoms(true) );
          residue[nResidues]->Copy ( res,mmdbRoot->Atom,mmdbRoot->nAtoms );
        } else  {
          for (i=0;i<res->nAtoms;i++)
            residue[nResidues]->AddAtom ( res->atom[i] );
        }
      } else  {
        residue[nResidues] = res;
        chain1 = res->GetChain();
        if (chain1)
          for (i=0;i<chain1->nResidues;i++)
            if (chain1->residue[i]==res)  {
              chain1->residue[i] = NULL;
              break;
            }
        residue[nResidues]->SetChain ( this );
        if (mmdbRoot)
          residue[nResidues]->CheckInAtoms();
      }
      nResidues++;

    }

    return nResidues;

  }
  */

  int  Chain::InsResidue ( PResidue res, int seqNum,
                            const InsCode insCode )  {
    return InsResidue ( res,GetResidueNo(seqNum,insCode) );
  }

  int  Chain::InsResidue ( PResidue res, int pos )  {
  //   Inserts residue res onto position pos of the chain,
  // pos=0..nResidues-1 . Residues pos..nResidues-1 are
  // shifted up the chain.
  //   The function places new atoms on the top of atom
  // index. It is advisable to call
  // CmmdbRoot::PDBCleanup ( PDBCLEAN_INDEX ) after all
  // insertions are done.
  PRoot   mmdbRoot;
  PChain  chain1;
  int     i,pp;

    pp = IMax ( 0,IMin(nResidues,pos) );

    for (i=0;i<nResidues;i++)
      if (residue[i]==res)  return -i;  // this residue is already there

    if (res)  {

      mmdbRoot = PRoot(GetCoordHierarchy());

      // get space for new residue
      if (nResidues>=resLen)
        ExpandResidueArray ( 100 );

      // shift residues to the end of the chain as necessary
      for (i=nResidues;i>pp;i--)
        residue[i] = residue[i-1];

      // insert the new residue
      if (res->GetCoordHierarchy())  {
        residue[pp] = newResidue();
        residue[pp]->SetChain ( this );
        residue[pp]->SetResID ( res->name,res->seqNum,res->insCode );
        if (mmdbRoot)  {
          // get space for new atoms
          mmdbRoot->AddAtomArray ( res->GetNumberOfAtoms(true) );
          residue[pp]->_copy ( res,mmdbRoot->atom,mmdbRoot->nAtoms );
        } else  {
          for (i=0;i<res->nAtoms;i++)
            residue[pp]->AddAtom ( res->atom[i] );
        }
      } else  {
        residue[pp] = res;
        chain1 = res->GetChain();
        if (chain1)
          for (i=0;i<chain1->nResidues;i++)
            if (chain1->residue[i]==res)  {
              chain1->residue[i] = NULL;
              break;
            }
        residue[pp]->SetChain ( this );
        if (mmdbRoot)
          residue[pp]->CheckInAtoms();
      }
      nResidues++;

    }

    return nResidues;

  }


  // --------------------  Extracting atoms  -----------------------

  int Chain::GetNumberOfAtoms ( bool countTers )  {
  int i,na;
    na = 0;
    for (i=0;i<nResidues;i++)
      if (residue[i])  na += residue[i]->GetNumberOfAtoms ( countTers );
    return na;
  }

  int Chain::GetNumberOfAtoms ( int seqNo, const InsCode insCode )  {
  PResidue res;
    res = GetResidue ( seqNo,insCode );
    if (res)  return res->nAtoms;
    return 0;
  }

  int Chain::GetNumberOfAtoms ( int resNo )  {
    if ((0<=resNo) && (resNo<nResidues))  {
      if (residue[resNo])  return residue[resNo]->nAtoms;
    }
    return 0;
  }

  PAtom Chain::GetAtom ( int            seqNo,
                         const InsCode  insCode,
                         const AtomName aname,
                         const Element  elmnt,
                         const AltLoc   aloc )  {
  PResidue res;
    res = GetResidue ( seqNo,insCode );
    if (res) return res->GetAtom ( aname,elmnt,aloc );
    return NULL;
  }

  PAtom Chain::GetAtom ( int seqNo, const InsCode insCode,
                         int atomNo )  {
  PResidue res;
    res = GetResidue ( seqNo,insCode );
    if (res)  {
      if ((0<=atomNo) && (atomNo<res->nAtoms))
        return res->atom[atomNo];
    }
    return NULL;
  }

  PAtom Chain::GetAtom ( int            resNo,
                         const AtomName aname,
                         const Element  elmnt,
                         const AltLoc   aloc )  {
    if ((0<=resNo) && (resNo<nResidues))  {
      if (residue[resNo])
        return residue[resNo]->GetAtom ( aname,elmnt,aloc );
    }
    return NULL;
  }

  PAtom Chain::GetAtom ( int resNo, int atomNo )  {
  PResidue res;
    if ((0<=resNo) && (resNo<nResidues))  {
      res = residue[resNo];
      if (res)  {
        if ((0<=atomNo) && (atomNo<res->nAtoms))
          return res->atom[atomNo];
      }
    }
    return NULL;
  }

  void Chain::GetAtomTable ( int seqNo, const InsCode insCode,
                           PPAtom & atomTable, int & NumberOfAtoms )  {
  PResidue res;
    atomTable     = NULL;
    NumberOfAtoms = 0;
    res = GetResidue ( seqNo,insCode );
    if (res)  {
      atomTable     = res->atom;
      NumberOfAtoms = res->nAtoms;
    }
  }

  void Chain::GetAtomTable ( int resNo, PPAtom & atomTable,
                              int & NumberOfAtoms )  {
  PResidue res;
    atomTable     = NULL;
    NumberOfAtoms = 0;
    if ((0<=resNo) && (resNo<nResidues))  {
      res = residue[resNo];
      if (res)  {
        atomTable     = res->atom;
        NumberOfAtoms = res->nAtoms;
      }
    }
  }


  void Chain::GetAtomTable1 ( int seqNo, const InsCode insCode,
                               PPAtom & atomTable, int & NumberOfAtoms )  {
  PResidue res;
    res = GetResidue ( seqNo,insCode );
    if (res)
      res->GetAtomTable1 ( atomTable,NumberOfAtoms );
    else  {
      if (atomTable)  delete[] atomTable;
      atomTable     = NULL;
      NumberOfAtoms = 0;
    }
  }

  void Chain::GetAtomTable1 ( int resNo, PPAtom & atomTable,
                               int & NumberOfAtoms )  {
  PResidue res;
    if ((0<=resNo) && (resNo<nResidues))
         res = residue[resNo];
    else res = NULL;
    if (res)
      res->GetAtomTable1 ( atomTable,NumberOfAtoms );
    else  {
      if (atomTable)  delete[] atomTable;
      atomTable     = NULL;
      NumberOfAtoms = 0;
    }
  }

  int Chain::DeleteAtom ( int            seqNo,
                           const InsCode  insCode,
                           const AtomName aname,
                           const Element  elmnt,
                           const AltLoc   aloc )  {
  PResidue res;
    res = GetResidue ( seqNo,insCode );
    if (res) return res->DeleteAtom ( aname,elmnt,aloc );
    return 0;
  }

  int Chain::DeleteAtom ( int seqNo, const InsCode insCode,
                           int atomNo )  {
  PResidue res;
    res = GetResidue ( seqNo,insCode );
    if (res) return res->DeleteAtom ( atomNo );
    return 0;
  }

  int Chain::DeleteAtom ( int            resNo,
                           const AtomName aname,
                           const Element  elmnt,
                           const AltLoc   aloc )  {
    if ((0<=resNo) && (resNo<nResidues))  {
      if (residue[resNo])
        return residue[resNo]->DeleteAtom ( aname,elmnt,aloc );
    }
    return 0;
  }

  int Chain::DeleteAtom ( int resNo, int atomNo )  {
    if ((0<=resNo) && (resNo<nResidues))  {
      if (residue[resNo])
        return residue[resNo]->DeleteAtom ( atomNo );
    }
    return 0;
  }


  int Chain::DeleteAllAtoms ( int seqNo, const InsCode insCode )  {
  PResidue res;
    res = GetResidue ( seqNo,insCode );
    if (res) return res->DeleteAllAtoms();
    return 0;
  }

  int Chain::DeleteAllAtoms ( int resNo )  {
    if ((0<=resNo) && (resNo<nResidues))  {
      if (residue[resNo])
        return residue[resNo]->DeleteAllAtoms();
    }
    return 0;
  }

  int Chain::DeleteAllAtoms()  {
  int i,k;
    k = 0;
    for (i=0;i<nResidues;i++)
      if (residue[i])
        k += residue[i]->DeleteAllAtoms();
    return k;
  }

  int Chain::DeleteAltLocs()  {
  //  This function leaves only alternative location with maximal
  // occupancy, if those are equal or unspecified, the one with
  // "least" alternative location indicator.
  //  The function returns the number of deleted. All tables remain
  // untrimmed, so that explicit trimming or calling FinishStructEdit()
  // is required.
  int i,n;

    n = 0;
    for (i=0;i<nResidues;i++)
      if (residue[i])  n += residue[i]->DeleteAltLocs();

    return n;

  }


  int Chain::AddAtom ( int seqNo, const InsCode insCode,
                        PAtom atom )  {
  PResidue res;
    res = GetResidue ( seqNo,insCode );
    if (res) return res->AddAtom ( atom );
    return 0;
  }

  int Chain::AddAtom ( int resNo, PAtom atom )  {
    if ((0<=resNo) && (resNo<nResidues))  {
      if (residue[resNo])
        return residue[resNo]->AddAtom ( atom );
    }
    return 0;
  }


  void  Chain::Copy ( PChain chain )  {
  // modify both Chain::_copy and Chain::Copy methods simultaneously!
  int i;

    FreeMemory();

    if (chain)  {

      CopyAnnotations ( chain );

      nResidues = chain->nResidues;
      resLen    = nResidues;
      if (nResidues>0)  {
        residue = new PResidue[nResidues];
        for (i=0;i<nResidues;i++)  {
          residue[i] = newResidue();
          residue[i]->SetChain ( this );
          residue[i]->Copy ( chain->residue[i] );
        }
      }

    }

  }

  void  Chain::CopyAnnotations ( PChain chain )  {
    if (chain)  {
      strcpy ( chainID    ,chain->chainID     );
      strcpy ( prevChainID,chain->prevChainID );
      DBRef .Copy ( &(chain->DBRef)  );
      seqAdv.Copy ( &(chain->seqAdv) );  //  SEQADV records
      seqRes.Copy ( &(chain->seqRes) );  //  SEQRES data
      modRes.Copy ( &(chain->modRes) );  //  MODRES records
      Het   .Copy ( &(chain->Het)    );  //  HET    records
    }
  }


  void  Chain::_copy ( PChain chain )  {
  // modify both Chain::_copy and Chain::Copy methods simultaneously!
  int i;

    FreeMemory();

    strcpy ( chainID    ,chain->chainID     );
    strcpy ( prevChainID,chain->prevChainID );

    DBRef .Copy ( &(chain->DBRef)  );
    seqAdv.Copy ( &(chain->seqAdv) );  //  SEQADV records
    seqRes.Copy ( &(chain->seqRes) );  //  SEQRES data
    modRes.Copy ( &(chain->modRes) );  //  MODRES records
    Het   .Copy ( &(chain->Het)    );  //  HET    records

    nResidues = chain->nResidues;
    resLen    = nResidues;
    if (nResidues>0)  {
      residue = new PResidue[nResidues];
      for (i=0;i<nResidues;i++)  {
        residue[i] = newResidue();
        residue[i]->SetChain ( this );
        residue[i]->_copy ( chain->residue[i] );
      }
    }

  }

  void  Chain::_copy ( PChain chain, PPAtom atom, int & atom_index )  {
  // modify both Chain::_copy and Chain::Copy methods simultaneously!
  int i;

    FreeMemory();

    strcpy ( chainID    ,chain->chainID     );
    strcpy ( prevChainID,chain->prevChainID );

    DBRef .Copy ( &(chain->DBRef)  );
    seqAdv.Copy ( &(chain->seqAdv) );  //  SEQADV records
    seqRes.Copy ( &(chain->seqRes) );  //  SEQRES data
    modRes.Copy ( &(chain->modRes) );  //  MODRES records
    Het   .Copy ( &(chain->Het)    );  //  HET    records

    nResidues = chain->nResidues;
    resLen    = nResidues;
    if (nResidues>0)  {
      residue = new PResidue[nResidues];
      for (i=0;i<nResidues;i++)
        if (chain->residue[i])  {
          residue[i] = newResidue();
          residue[i]->SetChain ( this );
          residue[i]->_copy ( chain->residue[i],atom,atom_index );
        } else
          residue[i] = NULL;
    }

  }

  /*
  void  Chain::Duplicate ( PChain Chain )  {
  int i;

    FreeMemory();

    strcpy ( chainID    ,chain->chainID     );
    strcpy ( prevChainID,chain->prevChainID );

    DBReference.Copy ( &(chain->DBReference) );
    SeqAdv     .Copy ( &(chain->SeqAdv)      );  //  SEQADV records
    SeqRes     .Copy ( &(chain->SeqRes)      );  //  SEQRES data
    ModRes     .Copy ( &(chain->ModRes)      );  //  MODRES records
    Het        .Copy ( &(chain->Het)         );  //  HET    records

    nResidues = chain->nResidues;
    resLen    = nResidues;
    if (nResidues>0)  {
      Residue = new PResidue[nResidues];
      for (i=0;i<nResidues;i++)  {
        residue[i] = newResidue();
        residue[i]->SetChain ( this );
        residue[i]->Duplicate ( chain->residue[i] );
      }
    }

  }
  */

  cpstr  Chain::GetEntryID()  {
    if (model)  return model->GetEntryID();
          else  return pstr("");
  }

  void  Chain::SetEntryID ( const IDCode idCode )  {
    if (model) model->SetEntryID ( idCode );
  }

  int   Chain::GetModelNum()  {
    if (model)  return model->GetSerNum();
    return 0;
  }

  cpstr  Chain::GetChainID ( pstr ChID )  {
    ChID[0] = char(0);
    if (model)
         sprintf ( ChID,"/%i/",model->GetSerNum() );
    else strcpy  ( ChID,"/-/" );
    strcat ( ChID,chainID );
    return ChID;
  }


  void  Chain::GetAtomStatistics  ( RAtomStat AS )  {
    AS.Init();
    CalAtomStatistics ( AS );
    AS.Finish();
  }

  void  Chain::CalAtomStatistics ( RAtomStat AS )  {
  int i;
    for (i=0;i<nResidues;i++)
      if (residue[i])
        residue[i]->CalAtomStatistics ( AS );
  }

  void  Chain::ApplyTransform ( mat44 & TMatrix )  {
  // transforms all coordinates by multiplying with matrix TMatrix
  int i;
    for (i=0;i<nResidues;i++)
      if (residue[i])  residue[i]->ApplyTransform ( TMatrix );
  }

  bool Chain::isSolventChain()  {
  // returns true if chain contains only solvent molecules
  bool B,P;
  int     i;
    B = true;
    P = false;
    for (i=0;(i<nResidues) && B;i++)
      if (residue[i])  {
        P = true;
        B = residue[i]->isSolvent();
      }
    return (B && P);
  }

  bool Chain::isInSelection ( int selHnd )  {
  PRoot mmdbRoot = (PRoot)GetCoordHierarchy();
  PMask mask;
    if (mmdbRoot)  {
      mask = mmdbRoot->GetSelMask ( selHnd );
      if (mask)  return CheckMask ( mask );
    }
    return false;
  }

  bool Chain::isAminoacidChain()  {
  // returns true if chain contains at least one aminoacid residue
  bool B,P;
  int     i;
    B = false;
    P = false;
    for (i=0;(i<nResidues) && (!B);i++)
      if (residue[i])  {
        P = true;
        B = residue[i]->isAminoacid();
      }
    return (B && P);
  }

  bool Chain::isNucleotideChain()  {
  // returns true if chain contains at least one nucleotide residue
  bool B,P;
  int     i;
    B = false;
    P = false;
    for (i=0;(i<nResidues) && (!B);i++)
      if (residue[i])  {
        P = true;
        B = residue[i]->isNucleotide();
      }
    return (B && P);
  }

  int  Chain::CheckID ( const ChainID chID )  {
    if (chID)  {
      if (!strcmp(chID,chainID))  return 1;
    }
    return 0;
  }

  int  Chain::CheckIDS ( cpstr CID )  {
  ChainID  chn;
  InsCode  inscode;
  ResName  resname;
  AtomName atm;
  Element  elm;
  AltLoc   aloc;
  int      mdl,sn,rc;

    rc = ParseAtomPath ( CID,mdl,chn,sn,inscode,resname,
                         atm,elm,aloc,NULL );
    if (rc>=0)  {
      if (!strcmp(chn,chainID))  return 1;
    }
    return 0;

  }

  int  Chain::GetNumberOfDBRefs()  {
    return  DBRef.Length();
  }

  PDBReference  Chain::GetDBRef ( int dbRefNo )  {
    return  (PDBReference)DBRef.GetContainerClass ( dbRefNo );
  }


  void  Chain::MaskAtoms ( PMask mask )  {
  int i;
    for (i=0;i<nResidues;i++)
      if (residue[i])  residue[i]->MaskAtoms ( mask );
  }

  void  Chain::MaskResidues ( PMask mask )  {
  int i;
    for (i=0;i<nResidues;i++)
      if (residue[i])  residue[i]->SetMask ( mask );
  }

  void  Chain::UnmaskAtoms ( PMask mask )  {
  int i;
    for (i=0;i<nResidues;i++)
      if (residue[i])  residue[i]->UnmaskAtoms ( mask );
  }

  void  Chain::UnmaskResidues ( PMask mask )  {
  int i;
    for (i=0;i<nResidues;i++)
      if (residue[i])  residue[i]->RemoveMask ( mask );
  }




  // -------  user-defined data handlers

  int  Chain::PutUDData ( int UDDhandle, int iudd )  {
    if (UDDhandle & UDRF_CHAIN)
          return  UDData::putUDData ( UDDhandle,iudd );
    else  return  UDDATA_WrongUDRType;
  }

  int  Chain::PutUDData ( int UDDhandle, realtype rudd )  {
    if (UDDhandle & UDRF_CHAIN)
          return  UDData::putUDData ( UDDhandle,rudd );
    else  return  UDDATA_WrongUDRType;
  }

  int  Chain::PutUDData ( int UDDhandle, cpstr sudd )  {
    if (UDDhandle & UDRF_CHAIN)
          return  UDData::putUDData ( UDDhandle,sudd );
    else  return  UDDATA_WrongUDRType;
  }

  int  Chain::GetUDData ( int UDDhandle, int & iudd )  {
    if (UDDhandle & UDRF_CHAIN)
          return  UDData::getUDData ( UDDhandle,iudd );
    else  return  UDDATA_WrongUDRType;
  }

  int  Chain::GetUDData ( int UDDhandle, realtype & rudd )  {
    if (UDDhandle & UDRF_CHAIN)
          return  UDData::getUDData ( UDDhandle,rudd );
    else  return  UDDATA_WrongUDRType;
  }

  int  Chain::GetUDData ( int UDDhandle, pstr sudd, int maxLen )  {
    if (UDDhandle & UDRF_CHAIN)
          return  UDData::getUDData ( UDDhandle,sudd,maxLen );
    else  return  UDDATA_WrongUDRType;
  }

  int  Chain::GetUDData ( int UDDhandle, pstr & sudd )  {
    if (UDDhandle & UDRF_CHAIN)
          return  UDData::getUDData ( UDDhandle,sudd );
    else  return  UDDATA_WrongUDRType;
  }


  //  -------------------------------------------------------------------

  DefineClass(SortResidues)

  class QSortResidues : public QuickSort  {
    public :
      QSortResidues() : QuickSort() {}
      int  Compare ( int i, int j );
      void Swap    ( int i, int j );
      void Sort    ( PPResidue res, int nresidues );
  };

  int QSortResidues::Compare ( int i, int j )  {
  int diff;
    diff = ((PPResidue)data)[i]->seqNum - ((PPResidue)data)[j]->seqNum;
    if (diff==0)
      diff = strcmp( (PPResidue(data))[i]->insCode,
                     (PPResidue(data))[j]->insCode );
    if (diff>0)  return  1;
    if (diff<0)  return -1;
    return 0;
  }

  void QSortResidues::Swap ( int i, int j )  {
  PResidue res;
    res = ((PPResidue)data)[i];
    ((PPResidue)data)[i] = ((PPResidue)data)[j];
    ((PPResidue)data)[j] = res;
  }

  void QSortResidues::Sort ( PPResidue res, int nresidues )  {
    QuickSort::Sort ( &(res[0]),nresidues );
  }

  void  Chain::SortResidues()  {
  QSortResidues SR;
    TrimResidueTable();
    SR.Sort ( residue,nResidues );
  }

  int  Chain::GetNofModResidues()  {
    return modRes.Length();
  }

  PModRes  Chain::GetModResidue ( int modResNo )  {
    return  PModRes(modRes.GetContainerClass(modResNo));
  }

  void  Chain::write ( io::RFile f )  {
  int  i;
  byte Version=2;
  bool compactBinary = false;
  
    PManager M = GetCoordHierarchy();
    if (M)
      compactBinary = M->isCompactBinary();

    f.WriteByte ( &Version       );
    f.WriteBool ( &compactBinary );
    f.WriteTerLine ( chainID,false );
    
    f.WriteInt ( &nResidues );
    for (i=0;i<nResidues;i++)
      residue[i]->write ( f );

    if (!compactBinary)  {

      UDData::write ( f );
  
      f.WriteTerLine ( prevChainID,false );
  
      DBRef .write ( f );  //  Database reference
      seqAdv.write ( f );  //  SEQADV records
      seqRes.write ( f );  //  SEQRES data
      modRes.write ( f );  //  MODRES records
      Het   .write ( f );  //  HET    records
      
    }

  }

  void  Chain::read ( io::RFile f )  {
  //   The Atom array in CmmdbRoot must be already read
  // prior to calling this function!
  int  i;
  byte Version;
  bool compactBinary;

    FreeMemory();

    f.ReadByte ( &Version       );
    f.ReadBool ( &compactBinary );
    f.ReadTerLine ( chainID,false );
    
    SetChain ( chainID );
    
    f.ReadInt ( &nResidues );
    resLen = nResidues;
    if (nResidues>0)  {
      residue = new PResidue[nResidues];
      for (i=0;i<nResidues;i++)  {
        residue[i] = newResidue();
        residue[i]->SetChain ( this );
        residue[i]->read ( f );
      }
    }

    if (!compactBinary)  {

      UDData::read ( f );

      f.ReadTerLine ( prevChainID,false );
  
      DBRef .read ( f );   //  Database reference
      seqAdv.read ( f );   //  SEQADV records
      seqRes.read ( f );   //  SEQRES data
      modRes.read ( f );   //  MODRES records
      Het   .read ( f );   //  HET    records
      
    }

  }


  MakeFactoryFunctions(Chain)

}  // namespace mmdb


// ===================================================================

  /*
void  TestChain() {
//  reads from 'in.chain', writes into
//  'out.chain' and 'abin.chain'
CFile    f;
char     S[81];
PChain  Chain;

  Chain = newChain();

  f.assign ( "in.chain",true );
  if (f.reset()) {
    while (!f.FileEnd()) {
      f.ReadLine ( S,sizeof(S) );
      chain->ConvertPDBString ( S );
    }
    f.shut();
  } else {
    printf ( " Can't open input file 'in.chain' \n" );
    delete Chain;
    return;
  }

  f.assign ( "out.chain",true );
  if (f.rewrite()) {
    chain->PDBASCIIDump ( f );
    f.shut();
  } else {
    printf ( " Can't open output file 'out.chain' \n" );
    delete Chain;
    return;
  }


  f.assign ( "mmdb.chain.bin",false );
  if (f.rewrite()) {
    chain->write ( f );
    f.shut();
  } else {
    printf ( "  Can't open binary chain file for writing.\n" );
    delete Chain;
    return;
  }

  delete Chain;
  printf ( "   Chain deleted.\n" );

  Chain = newChain();
  if (f.reset()) {
    chain->read ( f );
    f.shut();
  } else {
    printf ( "  Can't open binary chain file for reading.\n" );
    delete Chain;
    return;
  }

  f.assign ( "abin.chain",true );
  if (f.rewrite()) {
    chain->PDBASCIIDump ( f );
    f.shut();
  } else
    printf ( " Can't open output file 'abin.chain' \n" );

  delete Chain;

}
  */
