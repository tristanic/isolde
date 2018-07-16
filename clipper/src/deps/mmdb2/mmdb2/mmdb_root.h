//  $Id: mmdb_root.h $
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
//    24.07.15   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  MMDB_Root <interface>
//       ~~~~~~~~~
//       Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::Root
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2000-2015
//
//  =================================================================
//

#ifndef __MMDB_Root__
#define __MMDB_Root__

#include "mmdb_io_file.h"
#include "mmdb_uddata.h"
#include "mmdb_title.h"
#include "mmdb_cryst.h"
#include "mmdb_chain.h"
#include "mmdb_model.h"
#include "mmdb_defs.h"

namespace mmdb  {

  // =======================  Root  ===========================

  // special effect flags
  enum MMDB_READ_FLAG  {
    MMDBF_AutoSerials            = 0x00000001,
    MMDBF_NoCoordRead            = 0x00000002,
    MMDBF_SimRWBROOK             = 0x00000004,
    MMDBF_PrintCIFWarnings       = 0x00000008,
    MMDBF_EnforceSpaces          = 0x00000010,
    MMDBF_IgnoreDuplSeqNum       = 0x00000020,
    MMDBF_IgnoreSegID            = 0x00000040,
    MMDBF_IgnoreElement          = 0x00000080,
    MMDBF_IgnoreCharge           = 0x00000100,
    MMDBF_IgnoreNonCoorPDBErrors = 0x00000200,
    MMDBF_IgnoreUnmatch          = 0x00000400,
    MMDBF_IgnoreBlankLines       = 0x00000800,
    MMDBF_IgnoreHash             = 0x00001000,
    MMDBF_IgnoreRemarks          = 0x00002000,
    MMDBF_AllowDuplChainID       = 0x00004000,
    MMDBF_FixSpaceGroup          = 0x00008000,
    MMDBF_EnforceAtomNames       = 0x00010000,
    MMDBF_EnforceUniqueChainID   = 0x00020000,
    MMDBF_DoNotProcessSpaceGroup = 0x00040000,
    MMDBF_MakeCompactBinary      = 0x00080000
  };

  // MMDBF_EnforceUniqueChainID   will make MMDB to rename chains on
  //         reading a file such as to maintain chains uniquesness. This
  //         is supposed to work only with 1-letter chain IDs and only
  //         if chain names are interchanged in the file. For example,
  //         if file contains a sequence of chains named
  //
  //              A,B, A,B, A,B, A,B, A,B
  //
  //         and this flag is set on, the resulting chain names in
  //         MMDB will be:
  //
  //              A,B, A0,B0, A1,B1, A2,B2, A3,B3
  //

  // file types:
  enum MMDB_FILE_TYPE  {
    MMDB_FILE_Undefined = -1,
    MMDB_FILE_PDB       =  0,
    MMDB_FILE_CIF       =  1,
    MMDB_FILE_Binary    =  2
  };

  // cleanup flags:
  enum PDB_CLEAN_FLAG  {
    PDBCLEAN_ATNAME         = 0x00000001,
    PDBCLEAN_TER            = 0x00000002,
    PDBCLEAN_CHAIN          = 0x00000004,
    PDBCLEAN_CHAIN_STRONG   = 0x00000008,
    PDBCLEAN_ALTCODE        = 0x00000010,
    PDBCLEAN_ALTCODE_STRONG = 0x00000020,
    PDBCLEAN_SERIAL         = 0x00000040,
    PDBCLEAN_SEQNUM         = 0x00000080,
    PDBCLEAN_CHAIN_ORDER    = 0x00000100,
    PDBCLEAN_CHAIN_ORDER_IX = 0x00000200,
    PDBCLEAN_SOLVENT        = 0x00000400,
    PDBCLEAN_INDEX          = 0x00000800,
    PDBCLEAN_ELEMENT        = 0x00001000,
    PDBCLEAN_ELEMENT_STRONG = 0x00002000
  };

  // crystallographic info inquery
  enum MMDB_CRYST_FLAG  {
    CRRDY_NotPrecise       =  0x00000001,
    CRRDY_isTranslation    =  0x00000002,
    CRRDY_NoOrthCode       =  0x00000004,
    CRRDY_Complete         =  0,
    CRRDY_NoTransfMatrices = -1,
    CRRDY_Unchecked        = -2,
    CRRDY_Ambiguous        = -3,
    CRRDY_NoCell           = -4,
    CRRDY_NoSpaceGroup     = -5
  };

  DefineClass(Root);
  DefineStreamFunctions(Root);

  class Root : public UDData  {

    friend class Model;
    friend class Chain;
    friend class Residue;
    friend class Atom;

    public :

      Root ();
      Root ( io::RPStream Object );
      ~Root();

      void  FreeFileMemory();


      //  ---------------  Reading/Writing external files  ---------

      void  SetFlag        ( word Flag );
      void  RemoveFlag     ( word Flag );

      ERROR_CODE ReadPDBASCII   ( cpstr PDBFileName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );
      ERROR_CODE ReadPDBASCII1  ( cpstr PDBLFName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );
      ERROR_CODE ReadPDBASCII   ( io::RFile f );

      ERROR_CODE ReadCIFASCII   ( cpstr CIFFileName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );
      ERROR_CODE ReadCIFASCII1  ( cpstr CIFLFName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );
      ERROR_CODE ReadCIFASCII   ( io::RFile f );
      ERROR_CODE ReadFromCIF    ( mmcif::PData CIFD );

      // adds info from PDB file
      ERROR_CODE AddPDBASCII1   ( cpstr PDBLFName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );
      ERROR_CODE AddPDBASCII    ( cpstr PDBFileName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );

      // auto format recognition
      ERROR_CODE ReadCoorFile   ( cpstr LFName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );
      ERROR_CODE ReadCoorFile1  ( cpstr CFName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );
      ERROR_CODE ReadCoorFile   ( io::RFile f );

      ERROR_CODE WritePDBASCII  ( cpstr PDBFileName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );
      ERROR_CODE WritePDBASCII1 ( cpstr PDBLFName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );
      void       WritePDBASCII  ( io::RFile f );

      ERROR_CODE WriteCIFASCII  ( cpstr CIFFileName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );
      ERROR_CODE WriteCIFASCII1 ( cpstr CIFLFName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );

      ERROR_CODE ReadMMDBF      ( cpstr MMDBRootName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );
      ERROR_CODE ReadMMDBF1     ( cpstr MMDBLFName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );
      ERROR_CODE ReadMMDBF      ( io::File & f );
      ERROR_CODE WriteMMDBF     ( cpstr MMDBRootName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );
      ERROR_CODE WriteMMDBF1    ( cpstr MMDBLFName,
                                  io::GZ_MODE gzipMode=io::GZM_CHECK );
      void       WriteMMDBF     ( io::File & f );

      void  GetInputBuffer ( pstr Line, int & count );

      //  PutPDBString adds a PDB-keyworded string
      // to the existing structure. Note that the string
      // is namely added meaning that it will be the
      // last REMARK, last JRNL, last ATOM etc. string
      // -- always the last one in its group.
      ERROR_CODE PutPDBString ( cpstr PDBString );


      //  PDBCleanup(..) cleans coordinate part to comply with PDB
      // standards and MMDB "expectations":
      //
      //  PDBCLEAN_ATNAME  pads atom names with spaces to form
      //                   4-symbol names
      //  PDBCLEAN_TER     inserts TER cards in the end of each chain
      //  PDBCLEAN_CHAIN   generates 1-character chain ids instead of
      //                   those many-character
      //  PDBCLEAN_CHAIN_STRONG generates 1-character chain ids starting
      //                   from 'A' on for all ids, including the
      //                   single-character ones
      //  PDBCLEAN_ALTCODE generates 1-character alternative codes
      //                   instead of many-character ones
      //  PDBCLEAN_ALTCODE_STRONG generates 1-character alternative codes
      //                   from 'A' on for all codes, including the
      //                   single-character ones
      //  PDBCLEAN_SERIAL  puts serial numbers in due order
      //  PDBCLEAN_SEQNUM  renumbers all residues so that they go
      //                   incrementally-by-one without insertion codes
      //  PDBCLEAN_CHAIN_ORDER puts chains in order of atom's serial
      //                   numbers
      //  PDBCLEAN_CHAIN_ORDER_IX puts chains in order of atom's
      //                   indices internal to MMDB
      //  PDBCLEAN_SOLVENT moves solvent chains at the end of each model
      //
      //  Return codes (as bits):
      //  0                Ok
      //  PDBCLEAN_CHAIN   too many chains for assigning them
      //                   1-letter codes
      //  PDBCLEAN_ATNAME  element names were not available
      //  PDBCLEAN_ALTCODE too many alternative codes encountered.
      //
      word  PDBCleanup ( word CleanKey );

      //     Makes all atoms in chain 'chainID', in all models, as
      //  'Het' atoms if Make is set True, and makes them 'ordinary'
      //  atoms otherwise. 'Ter' is automatically removed when
      //  converting to 'Het' atoms, and is automatically added
      //  when converting to 'ordinary' atoms. This may cause
      //  disorder in serial numbers -- just call
      //  PDBClean(PDBCLEAN_SERIAL) when necessary to fix this.
      void  MakeHetAtoms ( cpstr chainID, bool Make );

      //  ---------------  Working with atoms by serial numbers  ---

      inline PPAtom GetAtomArray   ()  { return atom;   }
      inline int GetAtomArrayLength()  { return atmLen; } // strictly not for
                                               // use in applications!!
      PAtom  GetAtomI ( int index );   // returns Atom[index-1]

      //   PutAtom(..) puts atom with the specified properties
      // into the structure. The current model is used; if no model
      // is set (crModel==NULL), one is created. Coordinates and
      // other parameters of the atom need to be set separately.
      //   The place, at which the atom is put, is determined by
      // index. If index is positive, then it points onto (index-1)th
      // element of the Atom array (the index counts 1,2,... and
      // generally coincides with the atom's serial number). If
      // there is already an atom at this position in the system,
      // the new atom will REPLACE it. The corresponding residues
      // are automatically updated.
      //   If index is null (=0), the new atom will be put on
      // the top of the structure, i.e. it will be put into
      // (index=nAtoms+1)-th position.
      //   If index is negative, then the new atom is INSERTED
      // BEFORE the atom in the (-index)th position. For saving
      // the computational efforts, this WILL NOT cause the
      // recalculation of all atoms' serial numbers according
      // to their actual positions. It will be needed, however,
      // for putting the things in order at a certain point,
      // especially before writing an output ASCII file. NOTE
      // that this ordering is never done automatically.
      //   In a correct PDB file the serial number (serNum) is always
      // equal to its position (index). However here we allow them
      // to be different for easing the management of relations,
      // particularly the connectivity.
      //
      //   Limitation: if PutAtom implies creating new
      // chains/residues, these are always created on the top
      // of existing chains/residues.
      int   PutAtom ( int            index,
                      int            serNum,
                      const AtomName atomName,
                      const ResName  resName,
                      const ChainID  chainID,
                      int            seqNum,
                      const InsCode  insCode,
                      const AltLoc   altLoc,
                      const SegID    segID,
                      const Element  element );

      int   PutAtom (
                 int    index,    // same meaning as above
                 PAtom   A,       // pointer to completed atom class
                 int    serNum=0  // 0 means that the serial number
                                  // will be set equal to index.
                                  // Otherwise the serial number
                                  // is set to the specified
                                  // value
                    );


      //    RemoveAtom(..) removes atom at the specified index
      // in the Atom array. This index is always accessible
      // as Atom[index]->index. If this leaves a residue empty,
      // the residue is removed. If this leaves an empty chain,
      // the chain is removed as well; the same happens to the
      // model.
      void  RemoveAtom ( int index );

      int   FinishStructEdit();

      void  TrimModelTable();

      //  ----------------  Deleting models  -----------------------

      int  DeleteAllModels  ();
      bool GetNewChainID ( int modelNo, ChainID chID, int length=1 );

      //  ---------------  Enquiring -------------------------------
      
      bool isCompactBinary();

      int   CrystReady();
      //    Returns flags:
      // CRRDY_Complete       if crystallographic information is complete
      // CRRDY_NotPrecise     if cryst. inf-n is not precise
      // CRRDY_isTranslation  if cryst. inf-n contains translation
      // CRRDY_NoOrthCode      no orthogonalization code
      //    Fatal:
      // CRRDY_NoTransfMatrices  if transform. matrices were not
      //                         calculated
      // CRRDY_Unchecked         if cryst. inf-n was not checked
      // CRRDY_Ambiguous         if cryst. inf-n is ambiguous
      // CRRDY_NoCell            if cryst. inf-n is unusable
      // CRRDY_NoSpaceGroup      if space group is not set


      bool isCrystInfo   ();  // cell parameters and space group
      bool isCellInfo    ();  // cell param-s a,b,c, alpha,beta,gamma
      bool isSpaceGroup  ();  // space group on CRYST1 card
      bool isTransfMatrix();  // orthogonalizing/fractionalizing
                                 // matrices
      bool isScaleMatrix ();  // SCALEx PDB records
      bool isNCSMatrix   ();  // MTRIXx PDB records
      int  GetNumberOfNCSMatrices();
      int  GetNumberOfNCSMates   ();  // Returns the number of
                                      // NCS mates not given in
                                      // the file (iGiven==0)
      bool GetNCSMatrix  ( int NCSMatrixNo, // 0..N-1
                            mat44 & ncs_m, int & iGiven );

      int  GetNumberOfSymOps ();  // number of symmetry operations
      pstr GetSymOp ( int Nop );  // XYZ symmetry operation name


      //  -------------  User-Defined Data  ------------------------

      int RegisterUDInteger ( UDR_TYPE udr_type, cpstr UDDataID );
      int RegisterUDReal    ( UDR_TYPE udr_type, cpstr UDDataID );
      int RegisterUDString  ( UDR_TYPE udr_type, cpstr UDDataID );
      int GetUDDHandle      ( UDR_TYPE udr_type, cpstr UDDataID );

      //  ----------------------------------------------------------

      void  SetSyminfoLib ( cpstr syminfo_lib );
      pstr  GetSyminfoLib ();
      int   SetSpaceGroup ( cpstr spGroup );
      pstr  GetSpaceGroup ();
      pstr  GetSpaceGroupFix();

      void  GetAtomStatistics ( RAtomStat AS );

      void  SetIgnoreSCALEi ( bool ignoreScalei );

      //  SetCell(..) is for changing cell parameters
      void  SetCell ( realtype cell_a,
                      realtype cell_b,
                      realtype cell_c,
                      realtype cell_alpha,
                      realtype cell_beta,
                      realtype cell_gamma,
                      int      OrthCode=0 );

      //  PutCell(..) is for setting cell parameters
      void  PutCell ( realtype cell_a,
                      realtype cell_b,
                      realtype cell_c,
                      realtype cell_alpha,
                      realtype cell_beta,
                      realtype cell_gamma,
                      int      OrthCode=0 );

      int   GetCell ( realtype & cell_a,
                      realtype & cell_b,
                      realtype & cell_c,
                      realtype & cell_alpha,
                      realtype & cell_beta,
                      realtype & cell_gamma,
                      realtype & vol,
                      int      & OrthCode );

      int  GetRCell ( realtype & cell_as,
                      realtype & cell_bs,
                      realtype & cell_cs,
                      realtype & cell_alphas,
                      realtype & cell_betas,
                      realtype & cell_gammas,
                      realtype & vols,
                      int      & OrthCode );

      void GetROMatrix ( mat44 & RO );

      //  GetTMatrix(..) calculates and returns the coordinate
      //  transformation matrix, which converts orthogonal coordinates
      //  according to the symmetry operation number Nop and places
      //  them into unit cell shifted by cellshift_a a's, cellshift_b
      //  b's and cellshift_c c's.
      //
      //  Return 0 means everything's fine,
      //         1 there's no symmetry operation Nop defined
      //         2 fractionalizing/orthogonalizing matrices were not
      //           calculated
      //         3 cell parameters were not set up.
      int GetTMatrix ( mat44 & TMatrix, int Nop,
                       int cellshift_a, int cellshift_b,
                       int cellshift_c );

      //  GetUCTMatrix(..) calculates and returns the coordinate
      //  transformation matrix, which converts orthogonal coordinates
      //  according to the symmetry operation number Nop. Translation
      //  part of the matrix is being chosen such that point (x,y,z)
      //  has least distance to the center of primary (333) unit cell,
      //  and then it is shifted by cellshift_a a's, cellshift_b b's and
      //  cellshift_c c's.
      //
      //  Return 0 means everything's fine,
      //         1 there's no symmetry operation Nop defined
      //         2 fractionalizing/orthogonalizing matrices were not
      //           calculated
      //         3 cell parameters were not set up.
      int GetUCTMatrix ( mat44 & TMatrix, int Nop,
                         realtype x, realtype y, realtype z,
                         int cellshift_a, int cellshift_b,
                         int cellshift_c );

      //  GetFractMatrix(..) calculates and returns the coordinate
      //  transformation matrix, which converts fractional coordinates
      //  according to the symmetry operation number Nop and places them
      //  into unit cell shifted by cellshift_a a's, cellshift_b b's
      //  and cellshift_c c's.
      //
      //  Return 0 means everything's fine,
      //         1 there's no symmetry operation Nop defined
      //         2 fractionalizing/orthogonalizing matrices were not
      //           calculated
      //         3 cell parameters were not set up.
      int GetFractMatrix ( mat44 & TMatrix, int Nop,
                           int cellshift_a, int cellshift_b,
                           int cellshift_c );


      //  GetSymOpMatrix(..) returns the transformation matrix for
      //  Nop-th symmetry operator in the space group
      //
      //  Return 0 means everything's fine,
      //         1 there's no symmetry operation Nop defined
      //         2 fractionalizing/orthogonalizing matrices were not
      //           calculated
      //         3 cell parameters were not set up.
      //
      int GetSymOpMatrix ( mat44 & TMatrix, int Nop );


      int   AddNCSMatrix    ( mat33 & ncs_m, vect3 & ncs_v, int iGiven );
      int   GenerateNCSMates(); // 1: no NCS matrices, 0: Ok

      pstr  GetEntryID ();
      void  SetEntryID ( const IDCode idCode );

      int   GetNofExpDataRecs();
      pstr  GetExpDataRec ( int recNo );  // 0.. on

      int   GetNofMdlTypeRecs();
      pstr  GetMdlTypeRec ( int recNo );  // 0.. on

      int   GetFileType() { return FType; }

      void  Copy ( PRoot MMDBRoot );

      void  SetCompactBinary();  // leaves only coordinates in binary files

      // -------  user-defined data handlers
      int   PutUDData ( int UDDhandle, int      iudd );
      int   PutUDData ( int UDDhandle, realtype rudd );
      int   PutUDData ( int UDDhandle, cpstr    sudd );

      int   GetUDData ( int UDDhandle, int      & iudd );
      int   GetUDData ( int UDDhandle, realtype & rudd );
      int   GetUDData ( int UDDhandle, pstr sudd, int maxLen );
      int   GetUDData ( int UDDhandle, pstr     & sudd );

      // GetStructureTitle() returns the contents of TITLE record
      // unfolded into single line. If Title is missing, returns
      // contents of COMPND(:MOLECULE). If COMPND is missing, returns
      // HEADER. If Header is missing, returns PDB code. If no PDB
      // code is there, returns "Not available".
      pstr  GetStructureTitle ( pstr & L );

      PCryst GetCrystData()  { return &cryst; }

      PClassContainer GetUnparsedA()  { return &SA; }
      PClassContainer GetUnparsedB()  { return &SB; }
      PClassContainer GetUnparsedC()  { return &SC; }

    protected :

      word       Flags;    // special effect flags
      int        FType;    // type of last file operation:
                           //    -1 : none
                           //     0 : PDB
                           //     1 : CIF
                           //     2 : BIN
                           // encoded as MMDB_FILE_XXXXX above

      Title      title;    // title section
      Cryst      cryst;    // crystallographic information section
      UDRegister udRegister; // register of user-defined data

      int        nModels;  // number of models
      PPModel    model;    // array of models [0..nModels-1]

      int        nAtoms;   // number of atoms
      int        atmLen;   // length of Atom array
      PPAtom     atom;     // array of atoms ordered by serial numbers

      AtomPath   DefPath;  // default coordinate path

      ClassContainer SA;   // string container for unrecognized strings
                           // which are between the title and the
                           // crystallographic sections
      ClassContainer Footnote;  // string container for footnotes
      ClassContainer SB;   // string container for unrecognized strings
                           // which are between the crystallographic and
                           // the coordinate sections
      ClassContainer SC;   // string container for unrecognized strings
                           // following the coordinate section

      //  input buffer
      int          lcount;  // input line counter
      char         S[500];  // read buffer
      mmcif::PData CIF;     // CIF file manager

      PModel     crModel; // current model, used at reading a PDB file
      PChain     crChain; // current chain, used at reading a PDB file
      PResidue   crRes;   // current residue, used at reading a PDB file

      bool       Exclude;            // used internally
      bool       ignoreRemarks;      // used temporarily
      bool       allowDuplChID;      // used temporarily
      bool       enforceUniqueChID;  // used temporarily

      void       InitMMDBRoot    ();
      void       FreeCoordMemory ();
      void       ReadPDBLine     ( io::RFile f, pstr L, int maxlen );
      ERROR_CODE ReadPDBAtom     ( cpstr L );
      ERROR_CODE ReadCIFAtom     ( mmcif::PData CIFD   );
      ERROR_CODE CheckAtomPlace  ( int  index, cpstr L );
      ERROR_CODE CheckAtomPlace  ( int  index, mmcif::PLoop Loop );
      ERROR_CODE SwitchModel     ( cpstr L );
      ERROR_CODE SwitchModel     ( int nM );
      ERROR_CODE AllocateAtom    ( int           index,
                                   const ChainID chainID,
                                   const ChainID label_asym_id,
                                   const ResName resName,
                                   const ResName label_comp_id,
                                   int           seqNum,
                                   int           label_seq_id,
                                   int           label_entity_id,
                                   const InsCode insCode,
                                   bool          Replace );
      void  ExpandAtomArray ( int inc );
      void  AddAtomArray    ( int inc );

      void  ApplyNCSTransform ( int NCSMatrixNo );

      virtual void ResetManager();

      //  ---------------  Stream I/O  -----------------------------
      void  write ( io::RFile f );
      void  read  ( io::RFile f );

      // don't use _ExcludeModel in your applications!
      int   _ExcludeModel ( int serNum );

      int   CheckInAtom   ( int index, PAtom A );
      int   CheckInAtoms  ( int index, PPAtom A, int natms );

      virtual PMask GetSelMask ( int selHnd );

    private :
      int modelCnt;  // used only at reading files

  };



  //  isMMDBBIN will return
  //    -1   if file FName does not exist
  //     0   if file FName is likely a MMDB BIN (binary) file
  //     1   if file FName is not a MMDB BIN (binary) file
  //     2   if file FName is likely a MMDB BIN (binary) file,
  //         but of a wrong edition (i.e. produced by a lower
  //         version of MMDB).
  extern int isMMDBBIN ( cpstr FName, io::GZ_MODE gzipMode=io::GZM_CHECK );
  extern int isMMDBBIN ( io::RFile f );

  //  isPDB will return
  //    -1   if file FName does not exist
  //     0   if file FName is likely a PDB file
  //     1   if file FName is not a PDB file
  extern int isPDB ( cpstr FName, io::GZ_MODE gzipMode=io::GZM_CHECK,
                     bool IgnoreBlankLines=false );
  extern int isPDB ( io::RFile f, bool IgnoreBlankLines=false );

}  // namespace mmdb

#endif

