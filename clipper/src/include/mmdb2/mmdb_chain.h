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
//  **** Module  :  MMDB_Chain <interface>
//       ~~~~~~~~~
//  **** Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::ProModel     ( abstract Model class         )
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

#ifndef __MMDB_Chain__
#define __MMDB_Chain__

#include "mmdb_io_stream.h"
#include "mmdb_utils.h"
#include "mmdb_atom.h"
#include "mmdb_defs.h"
#include "imex.h"

namespace mmdb  {

  //  ====================  ProModel  ======================

  //    This class is a virtue needed only for defining certain
  // functions of Model, which are used by Chain and
  // Residue

  DefineClass(ProModel);
  DefineStreamFunctions(ProModel);

  DefineClass(Manager);

  class MMDB_IMEX ProModel : public UDData  {

    friend class Chain;

    public :

      ProModel  () : UDData () {}
      ProModel  ( io::RPStream Object ) : UDData ( Object ) {}
      ~ProModel () {}

      virtual cpstr GetEntryID () { return ""; }
      virtual void  SetEntryID ( const IDCode ) {}

      virtual int   AddChain ( PChain ) { return 0; }

      // returns pointer to Root
      virtual PManager GetCoordHierarchy() { return NULL; }

      //  GetNumberOfModels() returns TOTAL number of models
      virtual int GetNumberOfModels() { return 0;    }

      //  GetNumberOfAllAtoms() returns TOTAL number of atoms in
      // all models
      virtual int GetNumberOfAllAtoms() { return 0;    }

      //  returns pointer to the general Atom array
      virtual PPAtom     GetAllAtoms() { return NULL; }

      virtual int  GetSerNum       () { return 0; }

      virtual void ExpandAtomArray ( int )  {}
      virtual void AddAtomArray    ( int )  {}

    protected :

      virtual int  _ExcludeChain ( const ChainID ) { return 0; }

  };



  //  ====================  ChainContainer  ======================

  DefineClass(ChainContainer);
  DefineStreamFunctions(ChainContainer);

  class MMDB_IMEX ChainContainer : public ClassContainer  {

    public :

      ChainContainer  () : ClassContainer () {}
      ChainContainer  ( io::RPStream Object )
                          : ClassContainer ( Object ) {}
      ~ChainContainer () {}

      PContainerClass MakeContainerClass ( int ClassID );

      void  SetChain ( PChain Chain_Owner ); // must be set before using
                                              // the Container

      // special functions used in Model::GetCIF(..)
      cpstr Get1stChainID ();
      void  MoveByChainID ( const ChainID chainID,
                            PChainContainer chainContainer );

    protected :
      PChain chain;

  };


  //  ==================  ContainerChain  =====================

  DefineClass(ContainerChain);
  DefineStreamFunctions(ContainerChain);

  class MMDB_IMEX ContainerChain : public ContainerClass {

    friend class ChainContainer;

    public :

      ContainerChain ();
      ContainerChain ( PChain Chain_Owner  );
      ContainerChain ( io::RPStream Object ) : ContainerClass(Object) {}

      void SetChain   ( PChain Chain_Owner );

    protected :
      PChain  chain;
      ChainID chainID;  // just a copy of Chain->chainID

  };


  //  ==================  DBReference  ========================

  DefineClass(DBReference);
  DefineStreamFunctions(DBReference);

  class MMDB_IMEX DBReference : public ContainerChain  {

    public :

      int      seqBeg;      // initial seq num of the PDB seq-ce segment
      InsCode  insBeg;      // initial ins code of the PDB seq-ce segm-t
      int      seqEnd;      // ending seq number of the PDB seq-ce segm-t
      InsCode  insEnd;      // ending ins code of the PDB seq-ce segment
      DBName   database;    // sequence database name
      DBAcCode dbAccession; // sequence database accession code
      DBIdCode dbIdCode;    // sequence database identification code
      int      dbseqBeg;    // initial seq number of the database segment
      InsCode  dbinsBeg;    // ins code of initial residue of the segment
      int      dbseqEnd;    // ending seq number of the database segment
      InsCode  dbinsEnd;   // ins code of the ending residue of the seg-t

      DBReference ();
      DBReference ( PChain Chain_Owner );
      DBReference ( PChain Chain_Owner, cpstr S );
      DBReference ( io::RPStream Object );
      ~DBReference();

      void       PDBASCIIDump    ( pstr S, int N );
      void       MakeCIF         ( mmcif::PData CIF, int N );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      ERROR_CODE GetCIF          ( mmcif::PData CIF, int & n );
      CLASS_ID   GetClassID      () { return ClassID_DBReference; }

      void  Copy  ( PContainerClass DBRef );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitDBReference();

  };


  //  ====================  SeqAdv  ===========================

  DefineClass(SeqAdv);
  DefineStreamFunctions(SeqAdv);

  class MMDB_IMEX SeqAdv : public ContainerChain  {

    public :

      ResName  resName;     // residue name in conflict
      int      seqNum;      // residue sequence number
      InsCode  insCode;     // residue insertion code
      DBName   database;    // sequence database name
      DBAcCode dbAccession; // sequence database accession code
      ResName  dbRes;       // sequence database residue name
      int      dbSeq;       // sequence database sequence number
      pstr     conflict;    // conflict comment

      SeqAdv ();
      SeqAdv ( PChain Chain_Owner );
      SeqAdv ( PChain Chain_Owner, cpstr S );
      SeqAdv ( io::RPStream Object );
      ~SeqAdv();

      void       PDBASCIIDump    ( pstr S, int N );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      void       MakeCIF         ( mmcif::PData CIF, int N );
      ERROR_CODE GetCIF          ( mmcif::PData CIF, int & n );
      CLASS_ID   GetClassID      () { return ClassID_SeqAdv; }

      void  Copy  ( PContainerClass seqAdv );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitSeqAdv();

  };


  //  ==================  SeqRes  ========================

  DefineClass(SeqRes);
  DefineStreamFunctions(SeqRes);

  class MMDB_IMEX SeqRes : public io::Stream  {

    friend class Model;
    friend class Chain;

    public :

      int       numRes;   // number of residues in the chain
      PResName  resName;  // residue names

      SeqRes ();
      SeqRes ( io::RPStream Object );
      ~SeqRes();

      void       SetChain        ( PChain Chain_Owner );
      void       PDBASCIIDump    ( io::RFile f );
      ERROR_CODE ConvertPDBASCII ( cpstr  S );

      void  MakeCIF     ( mmcif::PData CIF );
      ERROR_CODE GetCIF ( mmcif::PData CIF );

      void  Copy  ( PSeqRes seqRes );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :
      PChain  chain;
      ChainID chainID;
      int     serNum;

      void InitSeqRes();
      void FreeMemory();

  };


  //  ==================  ModRes  ========================

  DefineClass(ModRes);
  DefineStreamFunctions(ModRes);

  class MMDB_IMEX ModRes : public ContainerChain  {

    public :

      ResName  resName;     // residue name used
      int      seqNum;      // residue sequence number
      InsCode  insCode;     // residue insertion code
      ResName  stdRes;      // standard residue name
      pstr     comment;     // description of the residue modification

      ModRes ();
      ModRes ( PChain Chain_Owner );
      ModRes ( PChain Chain_Owner, cpstr S );
      ModRes ( io::RPStream Object );
      ~ModRes();

      void       PDBASCIIDump    ( pstr S, int N );
      void       MakeCIF         ( mmcif::PData CIF, int N );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      ERROR_CODE GetCIF          ( mmcif::PData CIF, int & n );
      CLASS_ID   GetClassID      () { return ClassID_ModRes; }

      void  Copy  ( PContainerClass modRes );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitModRes();

  };


  //  ==================  HetRec  ===========================

  DefineClass(HetRec);
  DefineStreamFunctions(HetRec);

  class MMDB_IMEX HetRec : public ContainerChain  {

    public :

      ResName  hetID;       // Het identifier (right-justified)
      int      seqNum;      // sequence number
      InsCode  insCode;     // insertion code
      int      numHetAtoms; // number of HETATM records for the
                            // group present in the entry
      pstr     comment;     // text describing Het group

      HetRec ();
      HetRec ( PChain Chain_Owner );
      HetRec ( PChain Chain_Owner, cpstr S );
      HetRec ( io::RPStream Object );
      ~HetRec();

      void       PDBASCIIDump    ( pstr S, int N );
      void       MakeCIF         ( mmcif::PData CIF, int N );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      ERROR_CODE GetCIF          ( mmcif::PData CIF, int & n );
      CLASS_ID   GetClassID      () { return ClassID_Het; }

      void  Copy  ( PContainerClass Het );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitHetRec();

  };


  //  =================  Chain  =======================

  DefineFactoryFunctions(Chain);

  class MMDB_IMEX Chain : public UDData  {

    friend class DBReference;
    friend class SeqAdv;
    friend class SeqRes;
    friend class ModRes;
    friend class HetRec;
    friend class Residue;
    friend class Atom;
    friend class Model;
    friend class Root;
    friend class SelManager;
    friend class BondManager;
    friend class CoorManager;
    friend class Manager;

    public :

      ChainContainer DBRef;    // database reference
      ChainContainer seqAdv;   // SEQADV records
      SeqRes         seqRes;   // Sequence residues, SEQRES records
      ChainContainer modRes;   // modification descriptions
      ChainContainer Het;      // non-standard residues descriptions

      Chain ();  // SetModel() MUST be used after this constructor!
      Chain ( PProModel model, const ChainID chID );
      Chain ( io::RPStream Object );
      ~Chain();

      void FreeAnnotations();

      void SetModel ( PProModel    model );
      void SetChain ( const ChainID chID );

      PManager GetCoordHierarchy();   // PRoot

      //   ConvertXXXXX(..) functions do not check for record name
      // and assume that PDBString is at least 81 symbols long
      // (including the terminating null).
      ERROR_CODE ConvertDBREF  ( cpstr PDBString );
      ERROR_CODE ConvertSEQADV ( cpstr PDBString );
      ERROR_CODE ConvertSEQRES ( cpstr PDBString );
      ERROR_CODE ConvertMODRES ( cpstr PDBString );
      ERROR_CODE ConvertHET    ( cpstr PDBString );

      // This function should be used for testing purposes only.
      // A full PDB ASCII dump for all models and chains involved
      // is done by Root class.
      void  PDBASCIIDump     ( io::RFile f );

      void  PDBASCIIAtomDump ( io::RFile f );
      void  MakeAtomCIF      ( mmcif::PData CIF );


      //  -----------------  Extracting residues  -------------------------

      int GetNumberOfResidues(); // returns number of res-s in the chain
      PResidue GetResidue ( int resNo );  // returns resNo-th residue
                                          // in the chain;
                                          // 0<=resNo<nResidues

      //   GetResidue(..) returns pointer on residue, whose sequence
      // number and insert code are given in seqNum and insCode,
      // respectively. If such a residue is absent in the chain,
      // returns NULL.
      PResidue GetResidue ( int seqNum, const InsCode insCode );

      //   GetResidueNo(..) returns the residue number in the chain's
      // residues table. Residues are numbered as 0..nres-1 as they
      // appear in the coordinate file.
      //   If residue is not found, the function returns -1.
      int  GetResidueNo ( int seqNum, const InsCode insCode );

      void GetResidueTable ( PPResidue & resTable,
                             int & NumberOfResidues );

      //   GetResidueCreate(..) returns pointer on residue, whose name,
      // sequence number and insertion code are given by resName, seqNum
      // and insCode, respectively. If such a residue is absent in the
      // chain, one is created at the end of chain.
      //   If a residue with given sequence number and insertion code
      // is present in the chain but has a different name, the function
      // returns NULL unless Enforce is set True. In the latter case,
      // a new residue is still created at the end of chain, but there
      // is no guarantee that any function operating on the sequence
      // number and insertion code will work properly.
      PResidue GetResidueCreate ( const ResName resName, int seqNum,
                                  const InsCode insCode, bool Enforce );

      //   GetCoorSequence(...) returns sequence inferred from list
      // of residues (which may differ from one in the file header).
      // The sequence is returned as a null-terminated string 'seq'.
      // On input, 'seq' should be either NULL or allocated (in which
      // case the original allocation will be released).
      void GetCoordSequence ( pstr & seq );

      //  ------------------  Deleting residues  ----------------------

      int  DeleteResidue ( int resNo ); // returns num of deleted res-s
      int  DeleteResidue ( int seqNum, const InsCode insCode );
      int  DeleteAllResidues();
      int  DeleteSolvent    ();
      void TrimResidueTable ();  // do not forget to call after all dels

      //  -------------------  Adding residues  -----------------------

      //   AddResidue(..) adds residue to the chain, InsResidue inserts
      // the residue on the specified position of the chain (other
      // residues are shifted up to the end of chain). Position in the
      // chain may be specified by a serial number (that is position in
      // the residue table) or by seqNum and insCode of one of the
      // chain's residues (the new residue is then inserted before that
      // one). If the chain is associated with a coordinate hierarchy,
      // and residue 'res' is not, the latter is checked in
      // automatically. If residue 'res' belongs to any coordinate
      // hierarchy (even though that of the residue), it is *copied*
      // rather than simply taken over, and is checked in.
      //   If the chain is not associated with a coordinate hierarchy,
      // all added residues will be checked in automatically once the
      // chain is checked in.
      int  AddResidue ( PResidue res );
      int  InsResidue ( PResidue res, int pos );
      int  InsResidue ( PResidue res, int seqNum, const InsCode insCode );

      //  --------------------  Extracting atoms  ---------------------

      int  GetNumberOfAtoms ( bool countTers );
      int  GetNumberOfAtoms ( int seqNo, const InsCode insCode );
      int  GetNumberOfAtoms ( int resNo );

      PAtom GetAtom ( int            seqNo,
                      const InsCode  insCode,
                      const AtomName aname,
                      const Element  elmnt,
                      const AltLoc   aloc );
      PAtom GetAtom ( int seqNo, const InsCode insCode, int atomNo );
      PAtom GetAtom ( int            resNo,
                      const AtomName aname,
                      const Element  elmnt,
                      const AltLoc   aloc );
      PAtom GetAtom ( int resNo, int atomNo );

      void GetAtomTable ( int seqNo, const InsCode insCode,
                          PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable ( int resNo,
                          PPAtom & atomTable, int & NumberOfAtoms );

      //   GetAtomTable1(..) returns atom table without TER atoms and
      // without NULL atom pointers. NumberOfAtoms returns the actual
      // number of atom pointers in atomTable.
      //   atomTable is allocated withing the function. If it was
      // not set to NULL before calling the function, the latter will
      // attempt to deallocate it first.
      //   The application is responsible for deleting atomTable,
      // however it must not touch atom pointers, i.e. use simply
      // "delete[] atomTable;". Never pass atomTable from
      // GetAtomTable(..) into this function, unless you set it to NULL
      // before doing that.
      void GetAtomTable1 ( int seqNo, const InsCode insCode,
                           PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable1 ( int resNo,
                           PPAtom & atomTable, int & NumberOfAtoms );

      //  ---------------------  Deleting atoms  ----------------------

      int DeleteAtom ( int            seqNo,
                       const InsCode  insCode,
                       const AtomName aname,
                       const Element  elmnt,
                       const AltLoc   aloc );
      int DeleteAtom ( int            seqNo,
                       const InsCode  insCode,
                       int            atomNo );
      int DeleteAtom ( int            resNo,
                       const AtomName aname,
                       const Element  elmnt,
                       const AltLoc   aloc );
      int DeleteAtom ( int resNo, int atomNo );

      int DeleteAllAtoms ( int seqNo, const InsCode insCode );
      int DeleteAllAtoms ( int resNo );
      int DeleteAllAtoms ();

      //  DeleteAltLocs() leaves only alternative location with maximal
      // occupancy, if those are equal or unspecified, the one with
      // "least" alternative location indicator.
      //  The function returns the number of deleted. All tables remain
      // untrimmed, so that explicit trimming or calling
      // FinishStructEdit() is required.
      int DeleteAltLocs();

      //  ----------------------  Adding atoms  -----------------------

      int AddAtom ( int seqNo, const InsCode insCode, PAtom atom );
      int AddAtom ( int resNo, PAtom atom );

      //  -------------------------------------------------------------

      void  ApplyTransform ( mat44 & TMatrix );  // transforms all
                                           // coordinates by multiplying
                                           // with matrix TMatrix

      int    GetModelNum();
      PModel GetModel   ()  { return (PModel)model; }
      cpstr  GetChainID ()  { return chainID;       }
      void   SetChainID ( const ChainID chID );
      cpstr  GetChainID ( pstr  ChID );  // returns /m/c

      void  GetAtomStatistics ( RAtomStat AS );
      void  CalAtomStatistics ( RAtomStat AS );

      int   CheckID    ( const ChainID chID );
      int   CheckIDS   ( cpstr CID  );

      cpstr GetEntryID ();
      void  SetEntryID ( const IDCode idCode );

      int   GetNumberOfDBRefs ();
      PDBReference  GetDBRef ( int dbRefNo );  // 0..nDBRefs-1

      void  MaskAtoms      ( PMask Mask );
      void  MaskResidues   ( PMask Mask );
      void  UnmaskAtoms    ( PMask Mask );
      void  UnmaskResidues ( PMask Mask );

      void  SortResidues   ();

      int     GetNofModResidues();
      PModRes GetModResidue    ( int modResNo );  // 0.. on

      bool   isSolventChain   ();
      bool   isInSelection    ( int selHnd );
      bool   isAminoacidChain ();
      bool   isNucleotideChain();


      // -------  user-defined data handlers
      int   PutUDData ( int UDDhandle, int      iudd );
      int   PutUDData ( int UDDhandle, realtype rudd );
      int   PutUDData ( int UDDhandle, cpstr    sudd );

      int   GetUDData ( int UDDhandle, int      & iudd );
      int   GetUDData ( int UDDhandle, realtype & rudd );
      int   GetUDData ( int UDDhandle, pstr sudd, int maxLen );
      int   GetUDData ( int UDDhandle, pstr     & sudd );

      void  Copy            ( PChain chain );
      void  CopyAnnotations ( PChain chain );

      void  write ( io::RFile f );    // writes header to PDB binary file
      void  read  ( io::RFile f );    // reads header from PDB binary file

    protected :

      ChainID    chainID;     // chain ID
      ChainID    prevChainID; // if chain is renamed, its original
                              // name may be saved here.
      PProModel  model;       // pointer to model class

      int        nWeights;    // used externally for sorting
      realtype   Weight;      //   chains

      int        nResidues;   // number of residues
      PPResidue  residue;     // array of residues

      bool       Exclude;     // used internally

      void  InitChain ();
      void  FreeMemory();

      void  ExpandResidueArray ( int inc );
      //   _ExcludeResidue(..) excludes (but does not dispose!) a residue
      // from the chain. Returns 1 if the chain gets empty and 0
      // otherwise.
      int   _ExcludeResidue ( const ResName resName, int seqNum,
                              const InsCode insCode );
      void  _copy ( PChain chain );
      void  _copy ( PChain chain, PPAtom atom, int & atom_index );
      void  CheckInAtoms();

    private :
      int  resLen;      // length of Residue array

  };


  extern void  TestChain();  //  reads from 'in.chain', writes into
                             //  'out.chain' and 'abin.chain'

}  // namespace mmdb

#endif

