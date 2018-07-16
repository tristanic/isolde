//  $Id: mmdb_defs.h,v 1.27 2012/01/26 17:52:20 ekr Exp $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2015.
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
//    23.07.17   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :   MMDBF_Defs <interface>
//       ~~~~~~~~~
//  **** Project :   MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//
//      Definition of types, constants and important classes.
//
//  (C) E. Krissinel 2000-2017
//
//  =================================================================
//

#ifndef __MMDB_Defs__
#define __MMDB_Defs__

#include "mmdb_mattype.h"

namespace mmdb  {

  enum  MMDB_VERSION  {
    MAJOR_VERSION = 2,  //!< MMDB major version
    MINOR_VERSION = 0,  //!< MMDB minor version
    MICRO_VERSION = 15  //!< MMDB micro version
  };

  //  =======================  types  =================================

  typedef  char         IDCode  [16];   //!< ID code of PDB entry
  typedef  IDCode *     PIDCode;        //!< pointer to ID code
  typedef  PIDCode &    RPIDCode;       //!< ref. to pointer to ID code
  typedef  char         Date    [12];   //!< date DD-MMM-YYYY
  typedef  char         RecName [7];    //!< name of PDB record

  typedef  char         ChainID [10];   //!< chain ID
  typedef  ChainID *    PChainID;       //!< pointer to chain ID
  typedef  char         InsCode [10];   //!< insertion code
  typedef  char         DBName  [10];   //!< sequence database name
  typedef  char         DBAcCode[20];   //!< seq. database accession code
  typedef  DBAcCode *   PDBAcCode;      //!< pointer to seq. db acc code
  typedef  char         DBIdCode[20];   //!< seq. database ident-n code
  typedef  DBIdCode *   PDBIdCode;      //!< pointer to DBIdCode
  typedef  char         ResName [20];   //!< residue name
  typedef  ResName  *   PResName;       //!< pointer to residue name
  typedef  PResName *   PPResName;      //!< ptr to vector of residue names
  typedef  char         HelixID [20];   //!< helix ID
  typedef  char         StrandID[20];   //!< strand ID
  typedef  char         SheetID [20];   //!< sheet ID
  typedef  char         TurnID  [20];   //!< turn ID
  typedef  char         LinkRID [20];   //!< Refmac link ID

  typedef  char         SymGroup[100];  //!< group of space symmetry
  typedef  realtype     vect3   [3];    //!< vector of 3 real numbers
  typedef  vect3    *   pvect3;         //!< ptr to vector 3
  typedef  pvect3   &   rpvect3;        //!< ref to ptr to vector 3
  typedef  realtype     vect4   [4];    //!< vector of 4 real numbers
  typedef  vect3        mat33   [3];    //!< matrix 3x3 of real numbers

  typedef  vect4        mat44   [4];    //!< matrix 4x4 of real numbers
  typedef  mat44    *   pmat44;         //!< ptr to matrix 4x4
  typedef  mat44    &   rmat44;         //!< ref to matrix 4x4
  typedef  pmat44   *   ppmat44;        //!< ptr to ptr to matrix 4x4
  typedef  pmat44   &   rpmat44;        //!< ref to ptr to matrix 4x4
  typedef  mat33        mat633  [6];    //!< matrix 6x3x3 of real numbers

  typedef  char         AtomName[20];   //!< name of the atom
  typedef  AtomName *   PAtomName;      //!< pointer to atom name
  typedef  char         AltLoc  [20];   //!< alternate location indicator
  typedef  AltLoc   *   PAltLoc;        //!< pointer to alt loc indicator
  typedef  char         SegID   [20];   //!< segment identifier
  typedef  char         Element [10];   //!< chemical element name
  typedef  Element  *   PElement;       //!< ptr to chemical element name
  typedef  char         EnergyType[10]; //!< energy type name
  typedef  EnergyType * PEnergyType;    //!< pointer to energy type name

  // do not forget update this when change the above typedefs:
  typedef  char  maxMMDBName[40];


  //  =====================  constants  ===============================

  //   ANY_RES should be used in selection functions for specifying
  // "any residue" to select. Defined in mmdb_selmngr.cpp
  extern const int ANY_RES;

  //    PRNK_XXXXX are the print keys. PRNK_Silent supresses all print
  // inside mmdb_xxxx unless specifically ordered or catastrophic.
  // PRNK_SimRWBROOK instructs mmdb to issue, whenever possible and
  // necessary, printouts and warnings of RWBROOK (fortran) package.
  enum PRINT_KEY  {
    PRNK_Silent      = 0,
    PRNK_SimRWBROOK  = 1
  };

  //  Error_XXXX may be returned by XX::ConvertPDBString() and GetCIF(..)
  // functions.
  //  Error_WrongSection is returned if the string passed into function
  // does not belong to the corresponding PDB section.


  enum ERROR_CODE  {

    Error_EmptyCIF     =-1,  //!< used as signal at reading CIF files

    Error_NoError      = 0,
    Error_Ok           = 0,
    Error_WrongSection = 1,
    Error_WrongChainID = 2,
    Error_WrongEntryID = 3,

    //  Error_SEQRES_serNum is returned by CSeqRes::ConvertPDBASCII() if
    //  serial numbers of SEQRES records do not increment by 1
    Error_SEQRES_serNum = 4,

    //  Error_SEQRES_numRes is returned by CSeqRes::ConvertPDBASCII() if
    //  SEQRES records show different number of residues
    Error_SEQRES_numRes = 5,

    //  Error_SEQRES_extraRes is returned by CSeqRes::ConvertPDBASCII() if
    //  SEQRES contains more residues than specified
    Error_SEQRES_extraRes = 6,

    Error_NCSM_Unrecognized   = 7,
    Error_NCSM_AlreadySet     = 8,
    Error_NCSM_WrongSerial    = 9,
    Error_NCSM_UnmatchIG      = 10,

    Error_ATOM_Unrecognized   = 11,
    Error_ATOM_AlreadySet     = 12,
    Error_ATOM_NoResidue      = 13,
    Error_ATOM_Unmatch        = 14,

    Error_CantOpenFile        = 15,
    Error_UnrecognizedInteger = 16,
    Error_WrongModelNo        = 17,
    Error_DuplicatedModel     = 18,
    Error_NoModel             = 19,
    Error_ForeignFile         = 20,
    Error_WrongEdition        = 21,

    //  CIF specific
    Error_NotACIFFile         = 22,
    Error_NoData              = 23,
    Error_NoLoop              = 24,
    Error_NoStruct            = 25,
    Error_UnrecognCIFItems    = 26,
    Error_MissingCIFField     = 27,
    Error_EmptyCIFLoop        = 28,
    Error_EmptyCIFStruct      = 29,
    Error_UnexpEndOfCIF       = 30,
    Error_MissgCIFLoopField   = 31,
    Error_NotACIFStructure    = 32,
    Error_NotACIFLoop         = 33,
    Error_UnrecognizedReal    = 34,

    Error_NoSheetID           = 35,
    Error_WrongSheetID        = 36,
    Error_WrongStrandNo       = 37,

    //   Error_WrongNumberOfStrands may be issued when reading
    // sheet data from CIF
    Error_WrongNumberOfStrands = 38,

    //   Error_WrongSheetOrder may be issued when reading
    // sheet data from CIF
    Error_WrongSheetOrder      = 39,

    //   Error_HBondInconsistency may be issued when reading
    // sheet data from CIF
    Error_HBondInconsistency   = 40,

    //   Error_EmptyResidueName is issued when PDB ATOM record
    // does not have a residue name
    Error_EmptyResidueName     = 41,

    //   Error_DuplicateSeqNum is issued when PDB ATOM records
    // show the sequence number and insertion code assigned
    // to more than one residue name
    Error_DuplicateSeqNum      = 42,

    //   Error_NoLogicalName may be returned by file i/o functions
    // if the specified environmental variable for file name
    // is not found.
    Error_NoLogicalName        = 43,

    //   Error_EmptyFile may be returned at reading non-existing
    // coordinate files
    Error_EmptyFile            = 44,

    Error_Unknown              = 45,

    //   Error_CIF_EmptyRow is the event of encountering
    // an empty row in _atom_site loop. It is handled
    // internally and has no effect on API
    Error_CIF_EmptyRow      = 99999,

    Error_GeneralError1     = 10000

  };

  //  ClassID_XXXX are used by container classes for proper
  // creating containered classes when reading from binary file.

  enum CLASS_ID {
    ClassID_Template   ,
    ClassID_String     ,
    ClassID_ObsLine    ,
    ClassID_TitleLine  ,
    ClassID_CAVEAT     ,
    ClassID_Compound   ,
    ClassID_Source     ,
    ClassID_ExpData    ,
    ClassID_MdlType    ,
    ClassID_Author     ,
    ClassID_RevData    ,
    ClassID_Supersede  ,
    ClassID_Journal    ,
    ClassID_Remark     ,
    ClassID_DBReference,
    ClassID_SeqAdv     ,
    ClassID_ModRes     ,
    ClassID_Het        ,
    ClassID_NCSMatrix  ,
    ClassID_TVect      ,
    ClassID_Helix      ,
    ClassID_Turn       ,
    ClassID_Link       ,
    ClassID_LinkR      ,
    ClassID_CisPep
  };


  //  =====================  classes  ===============================

  DefineClass(Atom);
  DefineClass(Residue);
  DefineClass(Chain);
  DefineClass(Model);
  DefineClass(Manager);

}  // namespace mmdb

#endif

