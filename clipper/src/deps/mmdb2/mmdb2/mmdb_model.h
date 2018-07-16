//  $Id: mmdb_model.h $
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
//    10.05.15   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  MMDB_Model <interface>
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

#ifndef __MMDB_Model__
#define __MMDB_Model__

#include "mmdb_io_stream.h"
#include "mmdb_utils.h"
#include "mmdb_chain.h"
#include "mmdb_defs.h"
#include "imex.h"


namespace mmdb  {

  //  ====================  HetCompound  =======================

  DefineClass(HetCompound);
  DefineStreamFunctions(HetCompound);

  class MMDB_IMEX HetCompound : public io::Stream  {

    public :

      ResName  hetID;      // Het identifiers, right-justified
      pstr     comment;
      int      nSynonyms;
      psvector hetSynonym; // synonyms
      int      compNum;    // component number
      char     wc;         // '*' for water, otherwise space
      pstr     Formula;    // formulas

      HetCompound ( cpstr HetName );
      HetCompound ( io::RPStream Object );
      ~HetCompound();

      void  AddKeyWord     ( cpstr W, bool Closed );
      void  HETNAM_PDBDump ( io::RFile f );
      void  HETSYN_PDBDump ( io::RFile f );
      void  FORMUL_PDBDump ( io::RFile f );

      void  FormComString  ( pstr & F );
      void  FormSynString  ( pstr & F );
      void  FormForString  ( pstr & F );

      void  Copy  ( PHetCompound hetCompound );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void  InitHetCompound ( cpstr HetName );
      void  FreeMemory      ();

  };


  //  ====================  SSContainer  ======================

  DefineClass(SSContainer);
  DefineStreamFunctions(SSContainer);

  class MMDB_IMEX SSContainer : public ClassContainer  {

    public :

      SSContainer  () : ClassContainer() {}
      SSContainer  ( io::RPStream Object )
                      : ClassContainer ( Object ) {}
      ~SSContainer () {}

      PContainerClass MakeContainerClass ( int ClassID );

  };


  //  ====================  Helix  ============================

  DefineClass(Helix);
  DefineStreamFunctions(Helix);

  class MMDB_IMEX Helix : public ContainerClass  {

    public :
      int     serNum;      // serial number
      HelixID helixID;     // helix ID
      ResName initResName; // name of the helix's initial residue
      ChainID initChainID; // chain ID for the chain containing the helix
      int     initSeqNum;  // sequence number of the initial residue
      InsCode initICode;   // insertion code of the initial residue
      ResName endResName;  // name of the helix's terminal residue
      ChainID endChainID;  // chain ID for the chain containing the helix
      int     endSeqNum;   // sequence number of the terminal residue
      InsCode endICode;    // insertion code of the terminal residue
      int     helixClass;  // helix class
      pstr    comment;     // comment about the helix
      int     length;      // length of the helix

      Helix ();
      Helix ( cpstr S );
      Helix ( io::RPStream Object );
      ~Helix();

      void       PDBASCIIDump    ( pstr S, int N   );
      void       MakeCIF         ( mmcif::PData CIF, int N );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      ERROR_CODE GetCIF          ( mmcif::PData CIF, int & n );
      CLASS_ID   GetClassID      () { return ClassID_Helix; }

      void  Copy  ( PContainerClass Helix );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitHelix();

  };



  //  ====================  Strand  ============================

  DefineClass(Strand);
  DefineStreamFunctions(Strand);

  class MMDB_IMEX Strand : public io::Stream  {

    public :

      StrandID sheetID;     // sheet ID
      int      strandNo;    // strand number
      ResName  initResName; // name of the strand's initial residue
      ChainID  initChainID; // chain ID of initial residue in the strand
      int      initSeqNum;  // sequence number of the initial residue
      InsCode  initICode;   // insertion code of the initial residue
      ResName  endResName;  // name of the strand's terminal residue
      ChainID  endChainID;  // chain ID of terminal residue in the strand
      int      endSeqNum;   // sequence number of the terminal residue
      InsCode  endICode;    // insertion code of the terminal residue
      int      sense;       // sense of strand with respect to previous
                            //    strand
      AtomName curAtom;     // registration; atom name in current strand
      ResName  curResName;  // registration; residue name in current
                            //    strand
      ChainID  curChainID;  // registration; chain ID in current strand
      int      curResSeq;   // registration; res-e seq numb in current
                            //    strand
      InsCode  curICode;    // registration; ins code in current strand
      AtomName prevAtom;    // registration; atom name in previous strand
      ResName  prevResName; // registration; residue name in previous
                            //    strand
      ChainID  prevChainID; // registration; chain ID in previous strand
      int      prevResSeq;  // registration; res-e seq numb in previous
                            //    strand
      InsCode  prevICode;   // registration; ins code in previous strand

      Strand ();
      Strand ( io::RPStream Object );
      ~Strand();

      void       PDBASCIIDump    ( pstr  S );
      void       MakeCIF         ( mmcif::PData CIF );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      int        GetCIF          ( mmcif::PData CIF, cpstr sheet_id );

      void  Copy  ( PStrand Strand );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitStrand();

  };


  //  ====================  Sheet  ============================

  DefineClass(Sheet);
  DefineStreamFunctions(Sheet);

  class MMDB_IMEX Sheet : public io::Stream  {

    public :
      SheetID  sheetID;   // sheet ID
      int      nStrands;  // number of strands in the sheet
      PPStrand strand;    // array of strands

      Sheet ();
      Sheet ( io::RPStream Object );
      ~Sheet();

      void  FreeMemory();
      void  OrderSheet();

      void       PDBASCIIDump    ( io::RFile f );
      void       MakeCIF         ( mmcif::PData CIF );
      ERROR_CODE ConvertPDBASCII ( cpstr  S    );
      int        GetCIF          ( mmcif::PData CIF );

      void  Copy  ( PSheet sheet );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :
      void  InitSheet      ();
      void  CIFFindStrands ( mmcif::PData CIF, cpstr Category );
      void  TryStrand      ( int strand_no );
      int   GetStrand      ( int strand_no );

  };


  //  ====================  Sheets  ============================

  DefineClass(Sheets);
  DefineStreamFunctions(Sheets);

  class MMDB_IMEX Sheets : public io::Stream  {

    public :
      int     nSheets;
      PPSheet sheet;

      Sheets ();
      Sheets ( io::RPStream Object );
      ~Sheets();

      void  FreeMemory();

      void       PDBASCIIDump    ( io::RFile f );
      void       MakeCIF         ( mmcif::PData CIF );
      ERROR_CODE ConvertPDBASCII ( cpstr  S    );
      int        GetCIF          ( mmcif::PData CIF );

      void  Copy  ( PSheets Sheets );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :
      void  InitSheets    ();
      void  CIFFindSheets ( mmcif::PData CIF, cpstr Category );

  };


  //  ====================  Turn  ============================

  DefineClass(Turn);
  DefineStreamFunctions(Turn);

  class MMDB_IMEX Turn : public ContainerClass  {

    public :
      int     serNum;      // serial number
      TurnID  turnID;      // turn ID
      ResName initResName; // name of the turn's initial residue
      ChainID initChainID; // chain ID for the chain containing the turn
      int     initSeqNum;  // sequence number of the initial residue
      InsCode initICode;   // insertion code of the initial residue
      ResName endResName;  // name of the turn's terminal residue
      ChainID endChainID;  // chain ID for the chain containing the turn
      int     endSeqNum;   // sequence number of the terminal residue
      InsCode endICode;    // insertion code of the terminal residue
      pstr    comment;     // comment about the helix

      Turn ();
      Turn ( cpstr S );
      Turn ( io::RPStream Object );
      ~Turn();

      void       PDBASCIIDump    ( pstr S, int N   );
      void       MakeCIF         ( mmcif::PData CIF, int N );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      ERROR_CODE GetCIF          ( mmcif::PData CIF, int & n );
      CLASS_ID   GetClassID      () { return ClassID_Turn; }

      void  Copy  ( PContainerClass turn );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitTurn();

  };



  //  ====================  HetCompounds  =======================

  DefineClass(HetCompounds);
  DefineStreamFunctions(HetCompounds);

  class MMDB_IMEX HetCompounds : public io::Stream  {

    public :

      int            nHets;
      PPHetCompound  hetCompound;

      HetCompounds ();
      HetCompounds ( io::RPStream Object );
      ~HetCompounds();

      void  FreeMemory    ();

      void  PDBASCIIDump  ( io::RFile f );
      void  ConvertHETNAM ( cpstr S );
      void  ConvertHETSYN ( cpstr S );
      void  ConvertFORMUL ( cpstr S );

      void  MakeCIF       ( mmcif::PData CIF );
      ERROR_CODE GetCIF   ( mmcif::PData CIF );

      void  Copy  ( PHetCompounds hetCompounds );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :
      bool Closed;

      void  InitHetCompounds();
      int   AddHetName      ( cpstr H );

  };


  //  ===================  LinkContainer  =====================

  DefineClass(LinkContainer);
  DefineStreamFunctions(LinkContainer);

  class MMDB_IMEX LinkContainer : public ClassContainer  {

    public :

      LinkContainer  () : ClassContainer() {}
      LinkContainer  ( io::RPStream Object )
                       : ClassContainer ( Object ) {}
      ~LinkContainer () {}

      PContainerClass MakeContainerClass ( int ClassID );

  };


  //  ====================  Link  ============================

  DefineClass(Link);
  DefineStreamFunctions(Link);

  class MMDB_IMEX Link : public ContainerClass  {

    public :
      AtomName atName1;   // name of 1st linked atom
      AltLoc   aloc1;     // alternative location of 1st linked atom
      ResName  resName1;  // residue name of 1st linked atom
      ChainID  chainID1;  // chain ID of 1st linked atom
      int      seqNum1;   // sequence number of 1st linked atom
      InsCode  insCode1;  // insertion code of 1st linked atom
      AtomName atName2;   // name of 2nd linked atom
      AltLoc   aloc2;     // alternative location of 2nd linked atom
      ResName  resName2;  // residue name of 2nd linked atom
      ChainID  chainID2;  // chain ID of 2nd linked atom
      int      seqNum2;   // sequence number of 2nd linked atom
      InsCode  insCode2;  // insertion code of 2nd linked atom
      int      s1,i1,j1,k1;  // sym id of 1st atom
      int      s2,i2,j2,k2;  // sym id of 2nd atom
      realtype dist;      // link distance

      Link ();
      Link ( cpstr S );
      Link ( io::RPStream Object );
      ~Link();

      void       PDBASCIIDump    ( pstr S, int N   );
      void       MakeCIF         ( mmcif::PData CIF, int N );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      ERROR_CODE GetCIF          ( mmcif::PData CIF, int & n );
      CLASS_ID   GetClassID      () { return ClassID_Link; }

      void  Copy  ( PContainerClass link );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitLink();

  };

  //  ===================  LinkRContainer  ====================

  DefineClass(LinkRContainer);
  DefineStreamFunctions(LinkRContainer);

  class MMDB_IMEX LinkRContainer : public ClassContainer  {

    public :

      LinkRContainer  () : ClassContainer() {}
      LinkRContainer  ( io::RPStream Object )
                       : ClassContainer ( Object ) {}
      ~LinkRContainer () {}

      PContainerClass MakeContainerClass ( int ClassID );

  };


  //  ====================  LinkR  ============================

  DefineClass(LinkR);
  DefineStreamFunctions(LinkR);

  /*

  Garib's
  LINK             LYS A  27                     PLP A 255                PLPLYS
  LINK             MAN S   3                     MAN S   4                BETA1-4
  LINK        C6  BBEN B   1                O1  BMAF S   2                BEN-MAF
  LINK        OE2 AGLU A 320                C1  AMAF S   2                GLU-MAF
  LINK        OE2  GLU A  67        1.895   ZN   ZN  R   5                GLU-ZN
  LINK        NE2  HIS A  71        2.055   ZN   ZN  R   5                HIS-ZN
  LINK        O    ARG A  69        2.240   NA   NA  R   9                ARG-NA

  Coot's
  LINKR        O   VAL C 103                NA    NA C 401                VAL-NA
  LINKR        OD1 ASP D  58                NA    NA D 401                ASP-NA
  LINKR        O   ALA D  97                NA    NA D 401                ALA-NA
  LINKR        OG1 THR D  99                NA    NA D 401                THR-NA
  LINKR        O   SER D 101                NA    NA D 401                SER-NA
  LINKR        O   VAL D 103                NA    NA D 401                VAL-NA

  PDB's
  LINK         O   GLY A  49                NA    NA A6001     1555   1555  2.98
  LINK         OG1 THR A  51                NA    NA A6001     1555   1555  2.72
  LINK         OD2 ASP A  66                NA    NA A6001     1555   1555  2.72
  LINK         NE  ARG A  68                NA    NA A6001     1555   1555  2.93

  LINK         NE  ARG A  68                NA    NA A6001     1555   1555  2.93
  LINK         C21 2EG A   7                 C22 2EG B  19     1555   1555  1.56
  */

  class MMDB_IMEX LinkR : public ContainerClass  {

    public :
      LinkRID  linkRID;   // link name
      AtomName atName1;   // name of 1st linked atom
      AltLoc   aloc1;     // alternative location of 1st linked atom
      ResName  resName1;  // residue name of 1st linked atom
      ChainID  chainID1;  // chain ID of 1st linked atom
      int      seqNum1;   // sequence number of 1st linked atom
      InsCode  insCode1;  // insertion code of 1st linked atom
      AtomName atName2;   // name of 2nd linked atom
      AltLoc   aloc2;     // alternative location of 2nd linked atom
      ResName  resName2;  // residue name of 2nd linked atom
      ChainID  chainID2;  // chain ID of 2nd linked atom
      int      seqNum2;   // sequence number of 2nd linked atom
      InsCode  insCode2;  // insertion code of 2nd linked atom
      realtype dist;      // link distance

      LinkR ();
      LinkR ( cpstr S );
      LinkR ( io::RPStream Object );
      ~LinkR();

      void       PDBASCIIDump    ( pstr S, int N   );
      void       MakeCIF         ( mmcif::PData CIF, int N );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      ERROR_CODE GetCIF          ( mmcif::PData CIF, int & n );
      CLASS_ID   GetClassID      () { return ClassID_LinkR; }

      void  Copy  ( PContainerClass LinkR );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitLinkR();

  };



  //  ===================  CisPepContainer  =====================

  DefineClass(CisPepContainer);
  DefineStreamFunctions(CisPepContainer);

  class MMDB_IMEX CisPepContainer : public ClassContainer  {

    public :

      CisPepContainer  () : ClassContainer() {}
      CisPepContainer  ( io::RPStream Object )
                       : ClassContainer ( Object ) {}
      ~CisPepContainer () {}

      PContainerClass MakeContainerClass ( int ClassID );

  };


  //  =====================  CisPep  ===========================

  DefineClass(CisPep);
  DefineStreamFunctions(CisPep);

  class MMDB_IMEX CisPep : public ContainerClass  {

    public :
      int      serNum;   //  record serial number
      ResName  pep1;     //  residue name
      ChainID  chainID1; //  chain identifier 1
      int      seqNum1;  //  residue sequence number 1
      InsCode  icode1;   //  insertion code 1
      ResName  pep2;     //  residue name 2
      ChainID  chainID2; //  chain identifier 2
      int      seqNum2;  //  residue sequence number 2
      InsCode  icode2;   //  insertion code 2
      int      modNum;   //  model number
      realtype measure;  //  measure of the angle in degrees.

      CisPep ();
      CisPep ( cpstr S );
      CisPep ( io::RPStream Object );
      ~CisPep();

      void       PDBASCIIDump    ( pstr S, int N );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      CLASS_ID   GetClassID      () { return ClassID_CisPep; }

      void  Copy  ( PContainerClass cisPep );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitCisPep();

  };



  //  ====================  Model  ===============================

  enum SSE_RC  {
    SSERC_Ok           = 0,
    SSERC_noResidues   = 1,
    SSERC_noAminoacids = 2,
    SSERC_noSSE        = 3
  };

  enum SORT_CHAIN_DIR  {
    SORT_CHAIN_ChainID_Asc  = 0,
    SORT_CHAIN_ChainID_Desc = 1
  };

  DefineFactoryFunctions(Model);

  class MMDB_IMEX Model : public ProModel  {

    friend class Manager;
    friend class BondManager;
    friend class SelManager;
    friend class CoorManager;
    friend class Root;
    friend class Chain;
    friend class Residue;
    friend class Atom;

    public :

      Model ();  // SetMMDBFile() MUST be used after this constructor!
      Model ( PManager MMDBF, int serialNum );
      Model ( io::RPStream Object );
      ~Model();

      void   SetMMDBManager ( PManager MMDBM, int serialNum );
      PManager GetCoordHierarchy() { return manager; }

      //   GetChainCreate() returns pointer on chain, whose identifier
      // is given in chID. If such a chain is absent in the model,
      // it is created. If enforceUniqueChainID is true and chain with
      // the same first letter in chain ID already exists in the model,
      // then the new chain ID will be appended with a serial number
      // in order to keep it unique. The model will contain chains like
      // A, A0, A1, A2, ... in such cases.
      PChain GetChainCreate ( const ChainID chID,
                              bool enforceUniqueChainID );

      //   CreateChain() creates a new chain with chain ID regardless
      // the presence of same-ID chains in the model. This function
      // was introduced only for compatibility with older CCP4
      // applications and using it in any new developments should be
      // strictly discouraged.
      PChain CreateChain    ( const ChainID chID );

      cpstr  GetEntryID ();
      void   SetEntryID ( const IDCode idCode );

      int    GetSerNum  (); // returns the model's serial number

      cpstr  GetModelID ( pstr modelID );  // returns "/mdl"

      int    GetNumberOfModels  (); // returns TOTAL number of models
      int    GetNumberOfAtoms   ( bool countTers ); // returns number
                                                // of atoms in the model
      int    GetNumberOfResidues(); // returns number of residues in
                                    // the model


      //  ----------------  Extracting chains  --------------------------

      int  GetNumberOfChains();  // returns number of chains in the model
      bool GetNewChainID ( ChainID chID, int length=1 );
      //   GetChain() returns pointer on chain, whose identifier
      // is given in chID. If such a chain is absent in the model,
      // returns NULL.
      PChain GetChain ( const ChainID chID );
      PChain GetChain ( int chainNo ); // returns chainNo-th chain
                                        // in the model;
                                        // 0<=chainNo<nChains
      void GetChainTable ( PPChain & chainTable,
                           int & NumberOfChains );

      //  ------------------  Deleting chains  --------------------------

      int  DeleteChain        ( const ChainID chID );
      int  DeleteChain        ( int chainNo );
      int  DeleteAllChains    ();
      int  DeleteSolventChains();
      void TrimChainTable     ();

      //  -------------------  Adding chains  ---------------------------

      int  AddChain ( PChain chn );

      //  --------------------  Sort chains  ----------------------------

      void SortChains ( int sortKey ); // SORT_CHAIN_XXXX

      //  ----------------  Extracting residues  ------------------------

      int GetNumberOfResidues ( const ChainID chainID );
      int GetNumberOfResidues ( int   chainNo );
      PResidue GetResidue ( const ChainID chainID, int seqNo,
                             const InsCode insCode );
      PResidue GetResidue ( const ChainID chainID, int resNo );
      PResidue GetResidue ( int   chainNo, int seqNo,
                             const InsCode insCode );
      PResidue GetResidue ( int   chainNo, int resNo );
      int     GetResidueNo ( const ChainID chainID, int seqNo,
                             const InsCode insCode );
      int     GetResidueNo ( int   chainNo, int seqNo,
                             const InsCode insCode );
      void GetResidueTable ( PPResidue & resTable,
                             int & NumberOfResidues );
      void GetResidueTable ( const ChainID chainID,
                             PPResidue & resTable,
                             int & NumberOfResidues );
      void GetResidueTable ( int   chainNo, PPResidue & resTable,
                             int & NumberOfResidues );

      //  -----------------  Deleting residues  -------------------------

      int DeleteResidue ( const ChainID chainID, int seqNo,
                          const InsCode insCode );
      int DeleteResidue ( const ChainID chainID, int resNo );
      int DeleteResidue ( int   chainNo, int seqNo,
                          const InsCode insCode );
      int DeleteResidue ( int   chainNo, int resNo );
      int DeleteAllResidues ( const ChainID chainID );
      int DeleteAllResidues ( int   chainNo );
      int DeleteSolvent     (); // in difference of DeleteSolventChains,
                                // this will remove all solvent molecules
                                // from the file rather then
                                // solely-solvent chains
      int DeleteAllResidues ();

      //  ------------------  Adding residues  --------------------------

      int AddResidue ( const ChainID chainID, PResidue res );
      int AddResidue ( int   chainNo, PResidue res );

      //  -------------------  Extracting atoms  ------------------------

      int GetNumberOfAllAtoms(); // returns TOTAL number of atoms in all
                                 //    models
      PPAtom    GetAllAtoms (); // returns pointer to Atom array

      int   GetNumberOfAtoms ( const ChainID chainID, int seqNo,
                               const InsCode insCode );
      int   GetNumberOfAtoms ( int   chainNo, int seqNo,
                               const InsCode insCode );
      int   GetNumberOfAtoms ( const ChainID chainID, int resNo );
      int   GetNumberOfAtoms ( int   chainNo, int resNo );

      PAtom GetAtom ( const ChainID  chID,
                      int            seqNo,
                      const InsCode  insCode,
                      const AtomName aname,
                      const Element  elmnt,
                      const AltLoc   aloc );
      PAtom GetAtom ( const ChainID  chID,    int seqNo,
                      const InsCode  insCode, int atomNo );
      PAtom GetAtom ( const ChainID  chID,
                      int            resNo,
                      const AtomName aname,
                      const Element  elmnt,
                      const AltLoc   aloc );
      PAtom GetAtom ( const ChainID  chID,  int resNo, int atomNo );
      PAtom GetAtom ( int chNo,  int seqNo,
                      const InsCode  insCode,
                      const AtomName aname,
                      const Element  elmnt,
                      const AltLoc   aloc );
      PAtom GetAtom ( int chNo,  int seqNo, const InsCode insCode,
                      int atomNo );
      PAtom GetAtom ( int chNo,  int resNo,
                      const AtomName aname,
                      const Element  elmnt,
                      const AltLoc aloc );
      PAtom GetAtom ( int chNo,  int resNo, int atomNo );

      void GetAtomTable ( const ChainID chainID, int seqNo,
                          const InsCode insCode,
                          PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable ( int   chainNo,       int seqNo,
                          const InsCode insCode,
                          PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable ( const ChainID chainID, int resNo,
                          PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable ( int     chainNo,     int resNo,
                          PPAtom & atomTable, int & NumberOfAtoms );

      //   GetAtomTable1(..) returns atom table without TER atoms and
      // without NULL atom pointers. NumberOfAtoms returns the actual
      // number of atom pointers in atomTable.
      //   atomTable is allocated withing the function. If it was
      // not set to NULL before calling the function, the latter will
      // attempt to deallocate it first.
      //   The application is responsible for deleting atomTable,
      // however it must not touch atom pointers, i.e. use simply
      // "delete atomTable;". Never pass atomTable from GetAtomTable(..)
      // into this function, unless you set it to NULL before doing that.
      void GetAtomTable1 ( const ChainID chainID, int seqNo,
                           const InsCode insCode,
                           PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable1 ( int   chainNo, int seqNo,
                           const InsCode insCode,
                           PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable1 ( const ChainID chainID, int resNo,
                           PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable1 ( int     chainNo,     int resNo,
                           PPAtom & atomTable, int & NumberOfAtoms );

      void  GetAtomStatistics ( RAtomStat AS );
      void  CalAtomStatistics ( RAtomStat AS );


      //  --------------------  Deleting atoms  -------------------------

      int DeleteAtom ( const ChainID  chID,
                       int            seqNo,
                       const InsCode  insCode,
                       const AtomName aname,
                       const Element  elmnt,
                       const AltLoc   aloc );
      int DeleteAtom ( const ChainID  chID,    int seqNo,
                       const InsCode  insCode, int atomNo );
      int DeleteAtom ( const ChainID  chID,
                       int            resNo,
                       const AtomName aname,
                       const Element  elmnt,
                       const AltLoc   aloc );
      int DeleteAtom ( const ChainID  chID,  int resNo, int atomNo );
      int DeleteAtom ( int chNo,  int seqNo,
                       const InsCode  insCode,
                       const AtomName aname,
                       const Element  elmnt,
                       const AltLoc   aloc );
      int DeleteAtom ( int chNo,  int seqNo, const InsCode insCode,
                       int atomNo );
      int DeleteAtom ( int chNo,  int resNo,
                       const AtomName aname,
                       const Element  elmnt,
                       const AltLoc   aloc );
      int DeleteAtom ( int chNo,  int resNo, int atomNo );

      int DeleteAllAtoms ( const ChainID chID, int seqNo,
                           const InsCode insCode );
      int DeleteAllAtoms ( const ChainID chID, int resNo );
      int DeleteAllAtoms ( const ChainID chID );
      int DeleteAllAtoms ( int chNo, int seqNo, const InsCode insCode );
      int DeleteAllAtoms ( int chNo, int resNo );
      int DeleteAllAtoms ( int chNo );
      int DeleteAllAtoms ();

      //  DeleteAltLocs() leaves only alternative location with maximal
      // occupancy, if those are equal or unspecified, the one with
      // "least" alternative location indicator.
      //  The function returns the number of deleted. All tables remain
      // untrimmed, so that explicit trimming or calling
      // FinishStructEdit() is required.
      int DeleteAltLocs();


      //  ---------------------  Adding atoms  --------------------------

      int AddAtom ( const ChainID chID, int seqNo,
                    const InsCode insCode, PAtom atom );
      int AddAtom ( const ChainID chID, int resNo, PAtom  atom );
      int AddAtom ( int   chNo, int seqNo, const InsCode insCode,
                    PAtom atom );
      int AddAtom ( int   chNo, int resNo, PAtom  atom );


      //  ---------------------------------------------------------------

      //   ConvertPDBString(..) interprets PDB records DBREF, SEQADV,
      // SEQRES, MODRES.
      //   Returns zero if the line was converted, otherwise returns a
      // non-negative value of Error_XXXX.
      //   PDBString must be not shorter than 81 characters.
      ERROR_CODE ConvertPDBString ( pstr PDBString );

      // PDBASCIIDumpPS(..) makes output of PDB primary structure records
      // excluding cispeps
      void  PDBASCIIDumpPS   ( io::RFile f );

      // PDBASCIIDumpCP(..) makes output of cispep records
      void  PDBASCIIDumpCP   ( io::RFile f );

      // PDBASCIIDump(..) makes output of PDB coordinate (ATOM etc.)
      // records
      void  PDBASCIIDump     ( io::RFile f );

      void  MakeAtomCIF      ( mmcif::PData CIF );
      void  MakePSCIF        ( mmcif::PData CIF );
      ERROR_CODE GetCIF      ( mmcif::PData CIF );

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
      void  MoveChain ( PChain & m_chain, PPAtom m_atom,
                        PPAtom  atom, int & atom_index,
                        int  chain_ext );

      void  GetAIndexRange ( int & i1, int & i2 );

      void  MaskAtoms      ( PMask mask );
      void  MaskResidues   ( PMask mask );
      void  MaskChains     ( PMask mask );
      void  UnmaskAtoms    ( PMask mask );
      void  UnmaskResidues ( PMask mask );
      void  UnmaskChains   ( PMask mask );


      //  ----  Getting Secondary Structure Elements

      int  GetNumberOfHelices ();
      int  GetNumberOfSheets  ();

      PHelix   GetHelix      ( int serialNum ); // 1<=serNum<=NofHelices

      void     GetSheetID    ( int serialNum, SheetID sheetID );
                                                   // '\0' for none

      PSheet   GetSheet      ( int   serialNum ); //1<=serNum<=NofSheets
      PSheet   GetSheet      ( const SheetID sheetID ); // NULL for none
      int  GetNumberOfStrands ( int   sheetSerNum );
      int  GetNumberOfStrands ( const SheetID sheetID );
      PStrand  GetStrand     ( int   sheetSerNum,
                               int strandSerNum );
      PStrand  GetStrand     ( const SheetID sheetID,
                               int strandSerNum );

      inline PSSContainer GetHelices() { return &helices; }
      inline PSheets      GetSheets () { return &sheets;  }

      void  RemoveSecStructure();
      int   CalcSecStructure  ( bool flagBulge=true,
                                int aminoSelHnd=-1 );
  //    int   CalcSecStructure  ( bool flagBulge=true );

      PHetCompounds GetHetInfo() { return &hetCompounds; }
      void  RemoveHetInfo     ();


      //  ----  Working Links

      int   GetNumberOfLinks ();
      PLink          GetLink ( int serialNum ); // 1<=serNum<=NofLinks
      PLinkContainer GetLinks() { return &links; }

      void   RemoveLinks();
      void   AddLink    ( PLink link );

      //  ----  Working Refmac Links

      int    GetNumberOfLinkRs ();
      PLinkR          GetLinkR ( int serialNum ); // 1<=serNum<=NofLinks
      PLinkRContainer GetLinkRs() { return &linkRs; }

      void   RemoveLinkRs();
      void   AddLinkR   ( PLinkR linkR );


      //  ----  Working CisPeps

      int       GetNumberOfCisPeps();
      PCisPep          GetCisPep ( int CisPepNum );
      PCisPepContainer GetCisPeps() { return &cisPeps; }

      void  RemoveCisPeps();
      void  AddCisPep    ( PCisPep cisPep );



      void  ApplyTransform ( mat44 & TMatrix );  // transforms all
                                        // coordinates by multiplying
                                        // with matrix TMatrix

      bool isInSelection ( int selHnd );


      // -------  user-defined data handlers
      int   PutUDData ( int UDDhandle, int      iudd );
      int   PutUDData ( int UDDhandle, realtype rudd );
      int   PutUDData ( int UDDhandle, cpstr    sudd );

      int   GetUDData ( int UDDhandle, int      & iudd );
      int   GetUDData ( int UDDhandle, realtype & rudd );
      int   GetUDData ( int UDDhandle, pstr sudd, int maxLen );
      int   GetUDData ( int UDDhandle, pstr     & sudd );


      void  Copy             ( PModel model );
      void  CopyHets         ( PModel model );
      void  CopySecStructure ( PModel model );
      void  CopyLinks        ( PModel model );
      void  CopyLinkRs       ( PModel model );
      void  CopyCisPeps      ( PModel model );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      int             serNum;       // the model serial number
      PManager        manager;      // pointer to mmdbmanager class

      HetCompounds    hetCompounds; // information on heterocompounds
      SSContainer     helices;      // information on helices
      Sheets          sheets;       // information on sheets
      SSContainer     turns;        // information on turns
      LinkContainer   links;        // information on links
      LinkRContainer  linkRs;       // information on refmac links
      CisPepContainer cisPeps;      // information on cispeps

      int             nChains;      // number of chains
      int             nChainsAlloc; // actual length of Chain[]
      PPChain         chain;        // array of chains

      bool            Exclude;      // used internally

      void  InitModel        ();
      void  FreeMemory       ();
      void  ExpandChainArray ( int nOfChains );
      ERROR_CODE GetCIFPSClass ( mmcif::PData CIF, int ClassID );

      //   _ExcludeChain(..) excludes (but does not dispose!) a chain
      // from the model. Returns 1 if the chain gets empty and 0
      // otherwise.
      int   _ExcludeChain ( const ChainID chainID );

      //  _copy(PModel) does not copy atoms! -- not for use in
      // applications
      void  _copy ( PModel Model );

      //  _copy(PModel,PPAtom,int&) does copy atoms into array 'atom'
      // starting from position atom_index. 'atom' should be able to
      // accept all new atoms - no checks on the length of 'atom'
      // is being made. This function should not be used in applications.
      void  _copy ( PModel Model, PPAtom  atom, int & atom_index );

      void  CheckInAtoms  ();

  };

}  // namespace mmdb

#endif

