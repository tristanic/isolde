//  $Id: mmdb_rwbrook.h $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2008.
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
//    16.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  MMDB_RWBrook  <interface>
//       ~~~~~~~~~
//  **** Project :  MacroMolecular Data Base (MMDB), "new rwbrook"
//       ~~~~~~~~~
//  **** Functions : mmdb_f_init_     ( initializer                       )
//       ~~~~~~~~~~~ mmdb_f_quit_     ( disposer                          )
//                   autoserials_     ( switch to the autoserials regime  )
//                   setreadcoords_   ( switch for reading coordinates    )
//                   simrwbrook_      ( simulates old RWBROOK printout    )
//                   mmdb_f_openl_    ( associates a unit with a file     )
//                   mmdb_f_open_     ( associates a unit with a file     )
//                   mmdb_f_copy_     ( copies contents of units          )
//                   mmdb_f_delete_   ( deletes part of a unit            )
//                   mmdb_f_settype_  ( changes type of file and r/w mode )
//                   mmdb_f_setname_  ( changes file name                 )
//                   mmdb_f_write_    ( writes a data structure into file )
//                   mmdb_f_close_    ( closes and disposes a data str-re )
//                   mmdb_f_advance_  ( advances the internal pointer     )
//                   mmdb_f_rewd_     ( sets internal pointer on the top  )
//                   mmdb_f_bksp_     ( shifts int-l pointer 1 atom back  )
//                   mmdb_f_atom_     ( reads/writes atom properties      )
//                   mmdb_f_coord_    ( reads/writes atom coordinates     )
//                   mmdb_f_setcell_  ( sets the crystal cell parameters  )
//                   mmdb_f_wbspgrp_  ( sets the space group              )
//                   mmdb_f_rbspgrp_  ( gets the space group              )
//                   mmdb_f_wbcell_   ( sets the crystal cell parameters  )
//                   mmdb_f_rbcell_   ( gets the crystal cell parameters  )
//                   mmdb_f_rbcelln_  ( gets the crystal cell parameters  )
//                   mmdb_f_rbrcel_   ( gets the recipricol cell          )
//                   mmdb_f_rborf_    ( returns or fill transf. matrices  )
//                   mmdb_f_orthmat_  ( calc. standard othogonalisations  )
//                   mmdb_f_cvanisou_ ( converts between cryst-c units    )
//                   mmdb_f_wremark_  ( writes a remark statement         )
//                   mmdb_f_setter
//                   mmdb_f_sethet
//                   mmdb_f_getnofncsmates_
//                   rberrstop_       ( error messenger                   )
//                   rbcheckerr_      ( a simple error messenger          )
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//


#ifndef  __MMDB_RWBrook__
#define  __MMDB_RWBrook__

#include "mmdb_mattype.h"
#include "mmdb_machine_.h"


// ****** mmdb_rwbrook error codes

enum RWB_ERROR  {

  RWBERR_Ok                =  0,
  RWBERR_NoChannel         = -1,
  RWBERR_NoFile            = -2,
  RWBERR_NoLogicalName     = -3,

  RWBERR_CantOpenFile      = -4,
  RWBERR_WrongInteger      = -5,
  RWBERR_WrongModelNo      = -6,
  RWBERR_DuplicatedModel   = -7,
  RWBERR_ForeignFile       = -8,
  RWBERR_WrongEdition      = -9,

  RWBERR_ATOM_Unrecognd    = -10,
  RWBERR_ATOM_AlreadySet   = -11,
  RWBERR_ATOM_NoResidue    = -12,
  RWBERR_ATOM_Unmatch      = -13,

  RWBERR_NoAdvance         = -14,
  RWBERR_EmptyPointer      = -15,
  RWBERR_NoMatrices        = -16,

  RWBERR_NoCoordinates     = -17,

  RWBERR_Disagreement      = -18,
  RWBERR_NoOrthCode        = -19,
  RWBERR_NoCheck           = -20,

  RWBERR_NoCellParams      = -21,

  RWBERR_NotACIFFile       = -22,
  RWBERR_NoData            = -23,
  RWBERR_UnrecognCIFItems  = -24,
  RWBERR_MissingCIFField   = -25,
  RWBERR_EmptyCIFLoop      = -26,
  RWBERR_UnexpEndOfCIF     = -27,
  RWBERR_MissgCIFLoopField = -28,
  RWBERR_NotACIFStructure  = -29,
  RWBERR_NotACIFLoop       = -30,
  RWBERR_WrongReal         = -31,

  RWBERR_WrongChainID      = -32,
  RWBERR_WrongEntryID      = -33,
  RWBERR_SEQRES_serNum     = -34,
  RWBERR_SEQRES_numRes     = -35,
  RWBERR_SEQRES_exraRes    = -36,
  RWBERR_NCSM_Unrecogn     = -37,
  RWBERR_NCSM_AlreadySet   = -38,
  RWBERR_NCSM_WrongSerial  = -39,
  RWBERR_NCSM_UnmatchIG    = -40,
  RWBERR_NoModel           = -41,
  RWBERR_NoSheetID         = -42,
  RWBERR_WrongSheetID      = -43,
  RWBERR_WrongStrandNo     = -44,
  RWBERR_WrongNofStrands   = -45,
  RWBERR_WrongSheetOrder   = -46,
  RWBERR_HBondInconsis     = -47,
  RWBERR_EmptyResidueName  = -48,
  RWBERR_DuplicateSeqNum   = -49,
  RWBERR_GeneralError1     = -50,

  RWBERR_Error1            = -101,
  RWBERR_Error2            = -102,
  RWBERR_Error3            = -103

};

//  ***** mmdb_rwbrook warning flags
//        0x00004000 means "it's a warning"

enum RWB_WARNING  {
  RWBWAR_Warning       = 0x00004000,
  RWBWAR_RewOutput     = 0x00004010,
  RWBWAR_FileTop       = 0x00004020,
  RWBWAR_WrongSerial   = 0x00004040,
  RWBWAR_UnkFormFactor = 0x00004080,
  RWBWAR_AmbFormFactor = 0x00004100,
  RWBWAR_NoOccupancy   = 0x00004200,
  RWBWAR_NoTempFactor  = 0x00004400
};

// ------------------------------------------------------------------

//    mmdb_f_init_() makes a general initialization of the file system.
// It must be called ONLY ONCE from the top of an application.

//  FORTRAN equivalent:   subroutine MMDB_F_Init
//  ~~~~~~~~~~~~~~~~~~~

FORTRAN_SUBR ( MMDB_F_INIT, mmdb_f_init, (),(),() );


// ------------------------------------------------------------------

//    mmdb_f_quit_() disposes the file system. A correct use assumes that
// it will be called before an application quits.

//  FORTRAN equivalent:   subroutine MMDB_F_Quit
//  ~~~~~~~~~~~~~~~~~~~

FORTRAN_SUBR ( MMDB_F_QUIT, mmdb_f_quit, (),(),() );


// ------------------------------------------------------------------

//    autoserials_(..) switches On/Off the regime of automatical
// generation of atom serial numbers at reading from PDB ASCII file.
// The autoserials regime is On if iOnOff parameter is set to
// non-zero, and the regime is turned Off otherwise. The setting
// will last until next call to autoserials_(..)
//
//    When this regime is Off (default state), all atom serial
// numbers are expected to be in strictly incremental order and
// any deviation from this rule will cause end of reading and
// MMDB_F_Open_(..) will issue the RWBERR_AlreadySet error code. If
// this code is then passed to error messengers (rberrstop_(..) or
// rbcheckerr_(..)) the application will stop. It is Ok, however,
// for serial numbers to increment by 2 or more.
//
//    When the autoserials regime is On, MMDB_F_Open_(..) does not pay
// attention to the serial numbers and generates them for each
// atom in strict incremental-by-one. This will work correctly only
// if all atom records ("ATOM"/"HETATM", "SIGATM", "ANISOU" and
// "SIGUIJ") are grouped, for every atom, as they should (precisely,
// "ATOM" or "HETATM" opens the group, then "SIGATM", "ANISOU" and
// "SIGUIJ" should follow until next "ATOM"/"HETATM" is met).

//  FORTRAN equivalent:   subroutine AutoSerials ( iOnOff )
//  ~~~~~~~~~~~~~~~~~~~   integer  iOnOff

FORTRAN_SUBR ( AUTOSERIALS,autoserials,
               ( int * iOnOff ),
               ( int * iOnOff ),
               ( int * iOnOff ) );


// ------------------------------------------------------------------

//    setreadcoords_(..) switches On/Off the reading of atomic
// coordinates when mmdb_f_open_ is called for input. The coordinates
// will be read if iOnOff parameter is set to non-zero, otherwise
// the reading will stop on the coordinate section of PDB file.
// The setting will last until the next call to setreadcoords_(..).
//
//    By default, the coordinates are read.
//
//  FORTRAN equivalent:   subroutine SetReadCoords ( iOnOff )
//  ~~~~~~~~~~~~~~~~~~~   integer  iOnOff

FORTRAN_SUBR ( SETREADCOORDS,setreadcoords,
               ( int * iOnOff ),
               ( int * iOnOff ),
               ( int * iOnOff ) );


// ------------------------------------------------------------------

//    simrwbrook_(..) switches On/Off the regime of exact following
// the old fortran RWBROOK's way of issuing messages and warnings.
//
//    By default, this regime is switched off, which supresses all
// messages from mmdb_rwbrook unless directly ordered or catastrophic.
// Switching this regime on will make the printout of converted
// programs significantly closer to that resulting from the use of
// old fortran RWBROOK package. The setting will last until the
// next call to simrwbrook_(..).
//
//  FORTRAN equivalent:   subroutine SimRWBROOK ( iOnOff )
//  ~~~~~~~~~~~~~~~~~~~   integer  iOnOff

FORTRAN_SUBR ( SIMRWBROOK,simrwbrook,
               ( int * iOnOff ),
               ( int * iOnOff ),
               ( int * iOnOff ) );


// ------------------------------------------------------------------

//    mmdb_f_open_(..) associates a coordinate file with channel number
// given in iUnit. If iUnit was previously associated with another
// or the same file, the file gets complete logical reinitialization
// which means that all previous modifications to the file are lost
// unless stored on disk with mmdb_f_write_(..) or mmdb_f_close_(..).
//
//   If the file is to be opened for input (RWStat is set to
// "INPUT"), all contents of the file is read into memory. It may be
// then modified and written back into a file (same or different).
//
//   If the file is to be opened for output (RWStat is set to
// "OUTPUT"), no file is physically opened, and only empty data
// structure is created in the memory. It may then be added with the
// data and stored in a disk file.
//
//   If FType is set to "PDB" or " ", the physical file is assumed
// to be read or written in the PDB format. "CIF" option is reserved
// for mmCIF files and is not realized at present. "BIN" means
// binary format. Note that both file name and file type may be
// changed before writting the file (see mmdb_f_setname_(..) and
// mmdb_f_settype_(..)).
//
//   mmdb_f_open(..) sets an internal pointer to "before the first"
// atom in the file (therefore it should be advanced to get access
// to the first atom, see mmdb_f_advance1_(..)). This pointer is used
// for getting atomic coordinates and other atomic characteristics
// from the file structure or for storing them into the structure.
// The pointer may be advanced, backspaced or set to a specific
// position in the file structure (see below).
//
//   iRet returns the error code (defined above). Extended
// information on the error may be then obtained through the
// geterror_(..) function immediately after return from
// mmdb_f_open_(..).

//  FORTRAN equivalent:   subroutine MMDB_F_Open ( LName,RWStat,FType,
//  ~~~~~~~~~~~~~~~~~~~                        iUnit,iRet )
//                        character*(*)  LName,RWStat,FType
//                        integer        iUnit,iRet

FORTRAN_SUBR ( MMDB_F_OPEN, mmdb_f_open,
               (    // lengths-at-end list
                mmdb::machine::fpstr LName,  // logical name
                mmdb::machine::fpstr RWStat, // "INPUT" or "OUTPUT"
                mmdb::machine::fpstr FType,  // "PDB", "CIF", "BIN" or " "
                int * iUnit,      // channel number
                int * iRet,       // returns error code
                int   LName_len,  // fortran-hidden length of LName
                int   RWStat_len, // fortran-hidden length of RWStat
                int   FType_len   // fortran-hidden length of FType
               ), ( // lengths-in-structure list
                mmdb::machine::fpstr LName,
                mmdb::machine::fpstr RWStat,
                mmdb::machine::fpstr FType,
                int * iUnit,  int * iRet
               ), ( // lengths-follow list
                mmdb::machine::fpstr LName,   int LName_len,
                mmdb::machine::fpstr RWStat,  int RWStat_len,
                mmdb::machine::fpstr FType,   int FType_len,
                int * iUnit,   int * iRet
               ) );



// ------------------------------------------------------------------

//    mmdb_f_open1_(..) is equivalent to mmdb_f_open_(..) but takes directly
// the file name (FName) instead of logical name (LName).
//
//  FORTRAN equivalent:   subroutine MMDB_F_Open1 ( FName,RWStat,FType,
//  ~~~~~~~~~~~~~~~~~~~                         iUnit,iRet )
//                        character*(*)  FName,RWStat,FType
//                        integer        iUnit,iRet

FORTRAN_SUBR ( MMDB_F_OPEN1, mmdb_f_open1,
               (    // lengths-at-end list
                mmdb::machine::fpstr FName,  // file name
                mmdb::machine::fpstr RWStat, // "INPUT" or "OUTPUT"
                mmdb::machine::fpstr FType,  // "PDB", "CIF", "BIN" or " "
                int * iUnit,      // channel number
                int * iRet,       // returns error code
                int   FName_len,  // fortran-hidden length of FName
                int   RWStat_len, // fortran-hidden length of RWStat
                int   FType_len   // fortran-hidden length of FType
               ), ( // lengths-in-structure list
                mmdb::machine::fpstr FName,
                mmdb::machine::fpstr RWStat,
                mmdb::machine::fpstr FType,
                int * iUnit,  int * iRet
               ), ( // lengths-follow list
                mmdb::machine::fpstr FName,   int FName_len,
                mmdb::machine::fpstr RWStat,  int RWStat_len,
                mmdb::machine::fpstr FType,   int FType_len,
                int * iUnit,   int * iRet
               ) );



// ------------------------------------------------------------------

//    mmdb_f_copy_(..) copies the specified part(s) of iUnit2 into iUnit1.
// All data which contained in the corresponding part(s) of iUnit1
// before the copying, is destroyed.
//
//  FORTRAN equivalent:   subroutine MMDB_F_Copy ( iUnit1,iUnit2,
//  ~~~~~~~~~~~~~~~~~~~                        copyKey,iRet )
//                        integer  iUnit1,iUnit2,copyKey,iRet

FORTRAN_SUBR ( MMDB_F_COPY, mmdb_f_copy,
               (    // lengths-at-end list
                int * iUnit1,    // destination unit
                int * iUnit2,    // source unit
                int * copyKey,   // copy key:
                                 //  = 1  copy all
                                 //  = 2  copy all except coordinates
                                 //  = 3  copy title section only
                                 //  = 4  copy crystallographic
                                 //       section only
                                 //  = 5  copy coordinate section only
                                 // any other value does not do anything
                int * iRet       // return code:
                                 //   =0 if success
                                 //   =RWBERR_NoChannel if a unit
                                 //                   does not exist
                                 //   =RWBERR_NoFile    if a unit
                                 //                   was not opened
               ), ( // lengths-in-structure list
                int * iUnit1,  int * iUnit2,
                int * copyKey, int * iRet
               ), ( // lengths-follow list
                int * iUnit1,  int * iUnit2,
                int * copyKey, int * iRet
               ) );



// ------------------------------------------------------------------

//    mmdb_f_delete_(..) deletes the specified parts of iUnit-th unit.
// The delete keys are exactly the same as copy keys in mmdb_f_copy_(..).
//
//  FORTRAN equivalent:   subroutine MMDB_F_Delete ( iUnit1,iUnit2,
//  ~~~~~~~~~~~~~~~~~~~                        CopyAtoms,iRet )
//                        integer  iUnit1,iUnit2,CopyAtoms,iRet

FORTRAN_SUBR ( MMDB_F_DELETE, mmdb_f_delete,
               (    // lengths-at-end list
                int * iUnit,     // unit number; *iUnit<=0 means
                                 // "the last mentioned unit"
                int * delKey,    // delete key:
                                 //  = 1  delete all
                                 //  = 2  delete all except coordinates
                                 //  = 3  delete title section only
                                 //  = 4  delete crystallographic
                                 //       section only
                                 //  = 5  delete coordinate section only
                                 // any other value does not do anything
                int * iRet       // return code:
                                 //   =0 if success
                                 //   =RWBERR_NoChannel if a unit
                                 //                   does not exist
                                 //   =RWBERR_NoFile    if a unit
                                 //                   was not opened
               ), ( // lengths-in-structure list
                int * iUnit, int * delKey, int * iRet
               ), ( // lengths-follow list
                int * iUnit, int * delKey, int * iRet
               ) );




// ------------------------------------------------------------------

//   mmdb_f_settype_(..) changes the type and/or the read/write mode of
// a unit iUnit, previously initialized with mmdb_f_open_(..). The file
// is not read from or purged onto disk, no data change occurs.
//
//   iRet returns either RWBERR_NoChannel if the unit was not
// previously initialized by mmdb_f_open_(..), or RWBERR_Ok in the case
// of success.

//  FORTRAN equivalent:   subroutine MMDB_F_SetType ( iUnit,FType,
//  ~~~~~~~~~~~~~~~~~~~                           RWState,iRet )
//                        character*(*)  FType,RWState
//                        integer        iUnit,iRet

FORTRAN_SUBR ( MMDB_F_SETTYPE, mmdb_f_settype,
               (    // lengths-at-end list
                int * iUnit,     // unit number
                mmdb::machine::fpstr FType,  // "PDB", "CIF", "BIN" or " "
                mmdb::machine::fpstr RWStat, // "INPUT" or "OUTPUT"
                int * iRet,      // returns -1 if unit not found,
                                 // otherwise 0
                int   FType_len, // fortran-hidden length of FType
                int   RWStat_len // fortran-hidden length of RWStat
               ), ( // lengths-in-structure list
                int * iUnit,  mmdb::machine::fpstr FType,
                mmdb::machine::fpstr RWStat, int * iRet
               ), ( // lengths-follow list
                int * iUnit,
                mmdb::machine::fpstr FType,   int FType_len,
                mmdb::machine::fpstr RWStat,  int RWStat_len,
                int * iRet
               ) );



// ------------------------------------------------------------------

//   mmdb_f_setname_(..) changes the file name for a unit iUnit,
// previously initialized with mmdb_f_open_(..). The file is not
// read from or purged onto disk, no data change occurs.
//
//   iRet returns either RWBERR_NoChannel if the unit was not
// previously initialized by mmdb_f_open_(..), or RWBERR_Ok in the case
// of success.

//  FORTRAN equivalent:   subroutine MMDB_F_SetName ( iUnit,FName,
//  ~~~~~~~~~~~~~~~~~~~                           iRet )
//                        character*(*)  FName
//                        integer        iUnit,iRet

FORTRAN_SUBR ( MMDB_F_SETNAME, mmdb_f_setname,
               (    // lengths-at-end list
                int * iUnit,    // unit number
                mmdb::machine::fpstr FName,    // file name
                int * iRet,     // returns -1 if unit not found,
                                // otherwise 0
                int   FName_len // fortran-hidden length of FName
               ), ( // lengths-in-structure list
                int * iUnit, mmdb::machine::fpstr FName, int * iRet
               ), ( // lengths-follow list
                int * iUnit,
                mmdb::machine::fpstr FName, int FName_len,
                int * iRet
               ) );



// ------------------------------------------------------------------

//    mmdb_f_write_(..) writes content of unit iUnit into a disk file.
// iRet will be set to -1 if the unit was not previously opened with
// mmdb_f_open_(..). If writting was successful, iRet is set to 0,
// otherwise geterror(..) will return an information about the
// error occured.
//
//   Note that you may write even units associated with input in
// call to mmdb_f_open(..). The file type does not change unless
// explicitely changed with mmdb_f_settype_(..).

//  FORTRAN equivalent:   subroutine MMDB_F_Write ( iUnit,iRet )
//  ~~~~~~~~~~~~~~~~~~~   integer  iUnit,iRet

FORTRAN_SUBR ( MMDB_F_WRITE,  mmdb_f_write,
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ) );


// ------------------------------------------------------------------

//   mmdb_f_close_(..) acts as mmdb_f_write_(..) if unit iUnit has been
// associated with output in mmdb_f_open_(..) or in the call to
// mmdb_f_settype(..). After writing the file, the unit iUnit is
// completely disposed.
//   If unit iUnit is associated with input, mmdb_f_close_(..) merely
// disposes it and all the information contained will be lost.

//  FORTRAN equivalent:   subroutine MMDB_F_Close ( iUnit,iRet )
//  ~~~~~~~~~~~~~~~~~~~   integer  iUnit,iRet

FORTRAN_SUBR ( MMDB_F_CLOSE,  mmdb_f_close,
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ) );


// ------------------------------------------------------------------

//   If unit iUnit is associated with input, mmdb_f_advance_(..) sets
// the internal pointer on the next atom in the file. The atom
// properties may then be retrieved using mmdb_f_atom_(..) and
// mmdb_f_coord_(..). If iTer is set to 0, then 'ter' cards are
// completely ignored. If iTer is set to 1, then 'ter' card will
// cause return with iRet=1 with internal pointer left on this 'ter'
// card. iRet=2 mean end of file, and iRet=0 means that the pointer
// was successfully advanced to the next atom.
//
//   If unit iUnit is associated with output, mmdb_f_advance_(..)
// merely advances the pointer. No actual change in the data
// structure or on disk occurs. The new position will be filled with
// atom data after execution of mmdb_f_atom_(..) and/or mmdb_f_coord_(..).
// The pointer WILL NOT be advanced if none of these functions were
// called since last advancement, in which case iRet will return
// RWBERR_NoAdvance. After successfull advancement, iRet will
// return 0.

//  FORTRAN equivalent:   subroutine MMDB_F_Advance ( iUnit,iOut,iTer,iRet )
//  ~~~~~~~~~~~~~~~~~~~   integer  iUnit,iOut,iTer,iRet

//  Relation to the former XYZAdvance fortran subroutione:
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//      subroutine XYZAdvance ( iUnit,iOut,iTer,*,* )
//      integer       iRet
//      character*80  ErrLin
//      call MMDB_F_Advance ( iUnit,iOut,iTer,iRet )
//      if (iRet.eq.1)  return 1
//      if (iRet.eq.2)  return 2
//      if (iRet.eq.RWBERR_NoChannel) then
//        ErrLin = ' ERROR: in MMDB_F_ADVANCE file has not been opened'
//        call CCPErr ( 1,ErrLin )
//      endif
//      return
//      end
//
//  where parameter iOut IS NOT USED.

FORTRAN_SUBR ( MMDB_F_ADVANCE, mmdb_f_advance,
               (   // lengths-at-end list
                int * iUnit, // unit number
                int * iOut,  // output echo file
                int * iTer,  // FLAG =1, return iRet=1 if 'ter' card found
                             //      =0, do not return on 'ter' card
                int * iRet   // =0  if normal return
                             // =1  if return on 'ter' card (iTer=1)
                             // =2  if return on end of file
                             // =3  if return on 'hetatm' card
                             // =RWBERR_NoChannel if unit does not exist
                             // =RWBERR_NoAdvance if pointer was not
                             //                   advanced
               ), ( // lengths-in-structure list
                int * iUnit, int * iOut, int * iTer, int * iRet
               ), ( // lengths-follow list
                int * iUnit, int * iOut, int * iTer, int * iRet
               ) );



// ------------------------------------------------------------------

//   mmdb_f_rewd_(..) sets the internal pointer to the "begining" of
// the data structure associated with unit *iUnit. This means that
// one should the "advance" it with mmdb_f_advance_(..) in order
// to get access to the first atom.
//   iRet returns RWBERR_NoChannel if iUnit-th unit was not
// initialized, RWBWAR_RewOutput if the unit was associated with
// output, and 0 otherwise.

//  FORTRAN equivalent:   subroutine MMDB_F_Rewd ( iUnit,iRet )
//  ~~~~~~~~~~~~~~~~~~~   integer  iUnit,iRet

FORTRAN_SUBR ( MMDB_F_REWD, mmdb_f_rewd,
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ) );



// ------------------------------------------------------------------

//   mmdb_f_bksp_(..) shifts the internal pointer for one atom back in
// the data structure associated with unit *iUnit. This means that
// the combination of mmdb_f_advance1_(..) and mmdb_f_bksp_(..) leaves the
// pointer unchanged.
//   iRet returns RWBERR_NoChannel if iUnit-th unit was not
// initialized, and sets bit RWBWAR_RewOutput if the unit was
// associated with output, RWBWAR_FileTop if the pointer is already
// on the top of the structure, and 0 otherwise.

//  FORTRAN equivalent:   subroutine MMDB_F_Bksp ( iUnit,iRet )
//  ~~~~~~~~~~~~~~~~~~~   integer  iUnit,iRet

FORTRAN_SUBR ( MMDB_F_BKSP, mmdb_f_bksp,
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ) );



// ------------------------------------------------------------------

//   mmdb_f_seek_(..) sets the internal pointer to the specified
// position. If fPos==0, *iRet will return bit RWBWAR_FileTop.
//   If unit iUnit is associated with input, iRet will return 2 if
// fPos is given a value outside the file range, 1 if a 'ter' card
// is met and 3 if a 'hetatm' card is met.
//   iRet returns RWBERR_NoChannel if iUnit-th unit was not
// initialized, and RWBERR_EmptyPointer if fPos-th position in the
// input file is not occupied.

//  FORTRAN equivalent:   subroutine MMDB_F_Seek ( iUnit,fPos,iRet )
//  ~~~~~~~~~~~~~~~~~~~   integer  iUnit,fPos,iRet

FORTRAN_SUBR ( MMDB_F_SEEK, mmdb_f_seek,
               (    // lengths-at-end list
                int * iUnit,  // unit number
                int * fPos,   // position to set
                int * iRet    // return code:
                              //  0   Ok
                              //  1   'ter' card met
                              //  2   end of file
                              //  3   'hetatm' card met
                              //  <0 error:
                              //  RWBERR_NoChannel
                              //      iUnit was not
                              //      initialized
                              //  RWBERR_EmptyPointer
                              //      fPos-th position
               ), ( // lengths-in-structure list
                int * iUnit, int * fPos, int * iRet
               ), ( // lengths-follow list
                int * iUnit, int * fPos, int * iRet
               ) );



// ------------------------------------------------------------------

//   mmdb_f_atom_(..) reads or writes the atom name, residue name,
// chain name and other parameters listed from/into the
// data structure associated with the unit number iUnit (cf.
// mmdb_f_open_(..)). The position in the structure is adjusted
// with the help of mmdb_f_advance_(..), mmdb_f_rewd_(..) or
// mmdb_f_bksp_(..).

//  FORTRAN equivalent:   subroutine MMDB_F_Atom ( iUnit,iSer,AtNam,
//  ~~~~~~~~~~~~~~~~~~~                 ResNam,ChnNam,iResN,ResNo,
//                                      InsCod,AltCod,SegID,IZ,ID,
//                                      iRet )
//                        integer       iUnit,iSer,iResN,IZ,iRet
//                        character*(*) AtNam,ResNam,ChnNam,ResNo
//                        character*(*) InsCod,AltCod,SegID,ID

FORTRAN_SUBR ( MMDB_F_ATOM,  mmdb_f_atom,
           (    // lengths-at-end list
            int * iUnit,    // unit number
            int * iSer,     // atom serial number
            mmdb::machine::fpstr AtNam,    // atom name (left justified)
            mmdb::machine::fpstr ResNam,   // residue name
            mmdb::machine::fpstr ChnNam,   // chain name
            int * iResN,    // residue number as an integer
            mmdb::machine::fpstr ResNo,    // residue number as character (input only)
            mmdb::machine::fpstr InsCod,   // the insertion code
            mmdb::machine::fpstr AltCod,   // the alternate conformation code
            mmdb::machine::fpstr segID,    // segment ID
            int * IZ,       // atomic number (input only, returned as
                            // 7 from ambiguous atoms)
            mmdb::machine::fpstr ID,       // atomic ID related to atomic number
                            // (element symbol right justified), plus
                            // the ionic state +2, +3 etc..
                            //
            int * iRet,     // returns
                            //  RWBERR_NoChannel     if iUnit was not
                            //                       initialized
                            //  RWBERR_EmptyPointer  if atom was not
                            //                       advanced
                            //  RWBERR_Error1        internal error #1
                            //  RWBERR_Error2        internal error #2
                            //  RWBERR_Error3        internal error #3
                            //
                            //  >=0 : success, warning flags:
                            //  RWBWAR_WrongSerial   if serial number
                            //               differs from the position
                            //               number in the file
                            //  RWBWAR_UnkFormFactor unknown formfactor
                            //  RWBWAR_AmbFormFactor ambiguous formfactor
                            //
            int AtNam_len,  // fortran-hidden length of AtNam
            int ResNam_len, // fortran-hidden length of ResNam
            int ChnNam_len, // fortran-hidden length of ChnNam
            int ResNo_len,  // fortran-hidden length of ResNo
            int InsCod_len, // fortran-hidden length of InsCod
            int AltCod_len, // fortran-hidden length of AltCod
            int segID_len,  // fortran-hidden length of SegID
            int ID_len      // fortran-hidden length of ID
           ), ( // lengths-in-structure list
            int * iUnit,  int * iSer,
            mmdb::machine::fpstr AtNam,
            mmdb::machine::fpstr ResNam,
            mmdb::machine::fpstr ChnNam, int * iResN,
            mmdb::machine::fpstr ResNo,
            mmdb::machine::fpstr InsCod,
            mmdb::machine::fpstr AltCod,
            mmdb::machine::fpstr segID, int * IZ,
            mmdb::machine::fpstr ID,
            int * iRet
           ), ( // lengths-follow list
            int * iUnit,  int * iSer,
            mmdb::machine::fpstr AtNam,  int   AtNam_len,
            mmdb::machine::fpstr ResNam, int   ResNam_len,
            mmdb::machine::fpstr ChnNam, int   ChnNam_len,
            int * iResN,
            mmdb::machine::fpstr ResNo,  int   ResNo_len,
            mmdb::machine::fpstr InsCod, int   InsCod_len,
            mmdb::machine::fpstr AltCod, int   AltCod_len,
            mmdb::machine::fpstr segID,  int   segID_len,
            int * IZ,
            mmdb::machine::fpstr ID,     int   ID_len,
            int * iRet
           ) );


// ------------------------------------------------------------------

//    mmdb_f_copyatom_(..) copies atom from current position in
//  channel iUnit1 to current position in channel iUnit2.

//  FORTRAN equivalent:   subroutine MMDB_F_CopyAtom ( iUnit1,
//  ~~~~~~~~~~~~~~~~~~~                 iUnit2,iRet )
//                        integer iUnit1,iUnit2,iRet

FORTRAN_SUBR ( MMDB_F_COPYATOM, mmdb_f_copyatom,
               ( // length-at-end list
    int * iUnit1, // source channel number
                int * iUnit2, // destination channel number
                int * iRet    // returns
                              //  RWBERR_NoChannel     if iUnit was not
                              //                       initialized
                              //  RWBERR_EmptyPointer  if atom was not
                              //                       advanced
                              //  >=0 : success
         ), ( // length-in-structure list
                int * iUnit1, int * iUnit2, int * iRet
               ), ( // length-follow list
                int * iUnit1, int * iUnit2, int * iRet
               ) );


// ------------------------------------------------------------------

//   mmdb_f_coord_(..) reads or writes the atom coordinates, occupancy
// and temperature factor from/into the data structure associated
// with the unit number iUnit (cf. mmdb_f_open_(..)). The position in
// the structure is adjusted with the help of  mmdb_f_advance_(..),
// mmdb_f_rewd_(..) or mmdb_f_bksp_(..).
//   It is important that mmdb_f_coord_(..) was called AFTER
// mmdb_f_atom_(..) if channel iUnit is associated with output
// (otherwise iRet will return RWBERR_EmptyPointer).

//  FORTRAN equivalent:   subroutine MMDB_F_Coord ( iUnit,XFlag,BFlag,
//  ~~~~~~~~~~~~~~~~~~~                 x,y,z,occ,BIso,U,iRet )
//                        integer       iUnit,iRet
//                        character*(*) XFlag,BFlag
//                        real          x,y,z,occ,BIso,U(6)

//  Be sure that real-type parameters of mmdb_f_coord_(..) match those
//  of FORTRAN call. The real type is set with typedef apireal
//  statement in the begining of this file.

FORTRAN_SUBR ( MMDB_F_COORD, mmdb_f_coord,
               (    // lengths-at-end list
                int * iUnit,    // unit number
                mmdb::machine::fpstr XFlag,    // "F" or "O" flag for the fractional
                                // or orthogonal coordinates x,y,z
                                // for output files XFlag may also be
                                // set to "HF" or "HO", where "F" and
                                // "O" have the same meaning as before
                                // and "H" indicates that the atom
                                // should be marked as heteroatom
                mmdb::machine::fpstr BFlag ,   // "F" or "O" flag for temperature
                                // factor in fractional or orthogonal
                                // Us
                mmdb::machine::apireal * x,    // x-coordinate
                mmdb::machine::apireal * y,    // y-coordinate
                mmdb::machine::apireal * z,    // z-coordinate
                mmdb::machine::apireal * occ,  // occupancy
                mmdb::machine::apireal * BIso, // isotropic temperature factor
                mmdb::machine::apireal * U,    // array(6) of the anisotr. t-factor
                int * iRet,     // returns
                                //  RWBERR_NoChannel     if iUnit was not
                                //                       initialized
                                //  RWBERR_EmptyPointer  if atom was not
                                //                       advanced
                                //  RWBERR_NoMatrices    if transformation
                                //                       matrices are
                                //                       undefined
                                //  RWBERR_NoCoordinates if coordinates were
                                //                       not set in the atom
                                //
                                //  >=0 : success, warning flags:
                                //  RWBERR_NoOccupancy   if occupancy was
                                //                       not set in the atom
                                //  RWBERR_NoTempFactor  if temp. factor was
                                //                       not set in the atom
                                //
                int XFlag_len,  // fortran-hidden length of XFlag
                int BFlag_len   // fortran-hidden length of BFlag
               ), ( // lengths-in-structure list
                int * iUnit,
                mmdb::machine::fpstr XFlag,
                mmdb::machine::fpstr BFlag,
                mmdb::machine::apireal * x,
                mmdb::machine::apireal * y,
                mmdb::machine::apireal * z,
                mmdb::machine::apireal * occ,
                mmdb::machine::apireal * BIso,
                mmdb::machine::apireal * U,
                int * iRet
               ), ( // lengths-follow list
                int * iUnit,
                mmdb::machine::fpstr XFlag,   int XFlag_len,
                mmdb::machine::fpstr BFlag,   int BFlag_len,
                mmdb::machine::apireal * x,
                mmdb::machine::apireal * y,
                mmdb::machine::apireal * z,
                mmdb::machine::apireal * occ,
                mmdb::machine::apireal * BIso,
                mmdb::machine::apireal * U,
                int * iRet
               ) );


// ------------------------------------------------------------------

//   mmdb_f_setter_(..) sets the termination flag, so that the current
// atom will be converted into terminator of a chain and appear as
// 'ter' card in the output. The atom should be set by mmdb_f_atom_
// first, but the coordinates (mmdb_f_coord_(..)) do not have to be
// specified.

//  FORTRAN equivalent:   subroutine MMDB_F_SetTer ( iUnit,iRet )
//  ~~~~~~~~~~~~~~~~~~~   integer  iUnit,iRet

FORTRAN_SUBR ( MMDB_F_SETTER, mmdb_f_setter,
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ) );


// ------------------------------------------------------------------

//   mmdb_f_sethet_(..) sets the heteroatom flag, so that the current
// atom will appear as 'hetatm' card in the output. The atom should
// be set by mmdb_f_atom_ first and then mmdb_f_coord_(..) and
// mmdb_f_sethet_(..) should be called to specify its coordinates and
// heteroatom status.

//  FORTRAN equivalent:   subroutine MMDB_F_SetHet ( iUnit,iRet )
//  ~~~~~~~~~~~~~~~~~~~   integer  iUnit,iRet

FORTRAN_SUBR ( MMDB_F_SETHET, mmdb_f_sethet,
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ) );


// ------------------------------------------------------------------

//   mmdb_f_getnofncsmates_(..) returns the number of NCS mates not
//  given in the file (iGiven=0).
//
//  Negative returns N<0  mean an error.
//
//  FORTRAN equivalent:   subroutine MMDB_F_GetNofNCSMates ( iUnit,N )
//  ~~~~~~~~~~~~~~~~~~~   integer  iUnit,N

FORTRAN_SUBR ( MMDB_F_GETNOFNCSMATES, mmdb_f_getnofncsmates,
               ( int * iUnit, int * N ),
               ( int * iUnit, int * N ),
               ( int * iUnit, int * N ) );


// ------------------------------------------------------------------

//   mmdb_f_setcell_(..) sets the crystal cell properties and calculates
// the orthogonal-fractional transformation matrices for unit iUnit.

//  FORTRAN equivalent:   subroutine MMDB_F_SetCell ( iUnit,a,b,c,
//  ~~~~~~~~~~~~~~~~~~~                        alpha,beta,gamma,
//                                             ArgNCode,iRet )
//                        integer       iUnit,ArgNCode,iRet
//                        real          a,b,c,alpha,beta,gamma


//  Relation to the former RBFRAC2 FORTRAN subroutine:
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     subroutine  RBFRAC2 ( a,b,c,alpha,beta,gamma,ArgNCode )
//
//  ** the unit number iUnit and buffer for the return code iRet
//     have to be supplied.

FORTRAN_SUBR ( MMDB_F_SETCELL, mmdb_f_setcell,
               (   //   lengths-at-end list
                int     * iUnit,    // unit number; *iUnit<=0 means
                                    // "the last mentioned unit"
                mmdb::machine::apireal * a,     // cell parameter a, angstroms
                mmdb::machine::apireal * b,     // cell parameter b, angstroms
                mmdb::machine::apireal * c,     // cell parameter c, angstroms
                mmdb::machine::apireal * alpha, // cell parameter alpha, degrees
                mmdb::machine::apireal * beta,  // cell parameter beta,  degrees
                mmdb::machine::apireal * gamma, // cell parameter gamma, degrees
                int     * ArgNCode, // orthogonalization code, 1-6
                int     * iRet      // return code:
                                    //   RWBERR_Ok  - success
                                    //   RWBERR_NoChannel     if unit
                                    //              iUnit was not
                                    //              initialized
                                    //   RWBERR_NoFile        if unit
                                    //              has been disposed
                                    //   RWBERR_Disagreement  if a
                                    //              disagreement in
                                    //              cell parameters
                                    //              was found
                                    //   RWBERR_NoOrthCode    if no
                                    //              orthogonalization
                                    //              code was found
                                    //   RWBERR_NoCheck       if check
                                    //              of cell parameters
                                    //              has failed.
                                    //   The last three returns would
                                    // rather indicate a programming
                                    // error in mmdb_rwbrook.cpp
               ), ( // lengths-in-structure list
                int     * iUnit,
                mmdb::machine::apireal * a,
                mmdb::machine::apireal * b,
                mmdb::machine::apireal * c,
                mmdb::machine::apireal * alpha,
                mmdb::machine::apireal * beta,
                mmdb::machine::apireal * gamma,
                int     * ArgNCode, int     * iRet
               ), ( // lengths-follow list
                int     * iUnit,
                mmdb::machine::apireal * a,
                mmdb::machine::apireal * b,
                mmdb::machine::apireal * c,
                mmdb::machine::apireal * alpha,
                mmdb::machine::apireal * beta,
                mmdb::machine::apireal * gamma,
                int     * ArgNCode, int     * iRet
               )
             );



// ------------------------------------------------------------------

//   mmdb_f_wbspgrp_(..) sets the space group

//  FORTRAN equivalent:   subroutine MMDB_F_WBSpGrp ( iUnit,spGroup,iRet )
//  ~~~~~~~~~~~~~~~~~~~   integer       iUnit,iRet
//                        character*(*) spGroup


//  Relation to the former WBSpGrp FORTRAN subroutine:
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     subroutine  WBSpGrp ( spGroup )
//
//  ** the unit number iUnit and buffer for the return code iRet
//     have to be supplied.

FORTRAN_SUBR ( MMDB_F_WBSPGRP, mmdb_f_wbspgrp,
               (   //   lengths-at-end list
                int * iUnit,    // unit number; *iUnit<=0 means
                                // "the last mentioned unit"
                mmdb::machine::fpstr spGroup,  // space group
                int * iRet,     // return code:
                                //   RWBERR_Ok  - success
                                //   RWBERR_NoChannel     if unit
                                //              iUnit was not
                                //              initialized
                                //   RWBERR_NoFile        if unit
                                //              has been disposed
                int spGroup_len // fortran-hidden length of spGroup
               ), ( // lengths-in-structure list
                int * iUnit, mmdb::machine::fpstr spGroup, int * iRet
               ), ( // lengths-follow list
                int * iUnit, mmdb::machine::fpstr spGroup, int spGroup_len,
                int * iRet
               )
             );



// ------------------------------------------------------------------

//   mmdb_f_rbspgrp_(..) retrieves the space group

//  FORTRAN equivalent:   subroutine MMDB_F_RBSpGrp ( iUnit,spGroup,iRet )
//  ~~~~~~~~~~~~~~~~~~~   integer       iUnit,iRet
//                        character*(*) spGroup


//  Relation to the former RBSpGrp FORTRAN subroutine:
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     subroutine  RBSpGrp ( spGroup )
//
//  ** the unit number iUnit and buffer for the return code iRet
//     have to be supplied.

FORTRAN_SUBR ( MMDB_F_RBSPGRP, mmdb_f_rbspgrp,
               (   //   lengths-at-end list
                int * iUnit,    // unit number; *iUnit<=0 means
                                // "the last mentioned unit"
                mmdb::machine::fpstr spGroup,  // space group
                int * iRet,     // return code:
                                //   RWBERR_Ok  - success
                                //   RWBERR_NoChannel     if unit
                                //              iUnit was not
                                //              initialized
                                //   RWBERR_NoFile        if unit
                                //              has been disposed
                int spGroup_len // fortran-hidden length of spGroup
               ), ( // lengths-in-structure list
                int * iUnit, mmdb::machine::fpstr spGroup, int * iRet
               ), ( // lengths-follow list
                int * iUnit, mmdb::machine::fpstr spGroup, int spGroup_len,
                int * iRet
               )
             );



// ------------------------------------------------------------------

//   mmdb_f_wbcell_(..) sets the crystal cell properties into the
// channel iUnit.

//  FORTRAN equivalent:   subroutine MMDB_F_WBCell ( iUnit,ArgCell,
//  ~~~~~~~~~~~~~~~~~~~                        ArgNCode,iRet )
//                        integer iUnit,ArgNCode,iRet
//                        real    ArgCell(6)
//

//  Relation to the former WBCELL FORTRAN subroutine:
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     subroutine  WBCELL ( iUnit,ArgCell,ArgNCode )
//
//  ** the buffer for the return code iRet has to be supplied

FORTRAN_SUBR ( MMDB_F_WBCELL, mmdb_f_wbcell,
               (    // lengths-at-end list
                int     * iUnit,    // unit number; *iUnit<=0 means
                                    // "the last mentioned unit"
                mmdb::machine::apireal * ArgCell,  // array to accept the cell parameters
                                    // if ArgCell(1) is set to 0, then
                                    // the cell does not change
                int     * ArgNCode, // orthogonalisation code
                                    // if ArgNCode is set to 0, then
                                    // the orthogonalisation matrices
                                    // do not change
                int     * iRet      // return code
                                    //   RWBERR_Ok  - success
                                    //   RWBERR_NoChannel     if unit
                                    //              iUnit was not
                                    //              initialized
                                    //   RWBERR_NoFile        if unit
                                    //              has been disposed
               ), ( // lengths-in-structure list
                int * iUnit,    mmdb::machine::apireal * ArgCell,
                int * ArgNCode, int     * iRet
               ), ( // lengths-follow list
                int * iUnit,    mmdb::machine::apireal * ArgCell,
                int * ArgNCode, int     * iRet
               )
             );



// ------------------------------------------------------------------

//   mmdb_f_rbcell_(..) retrieves the crystal cell properties from the
// channel iUnit.

//  FORTRAN equivalent:   subroutine MMDB_F_RBCell ( iUnit,celld,cvol,iRet )
//  ~~~~~~~~~~~~~~~~~~~   integer       iUnit,iRet
//                        real          celld(6),cvol
//                        character*(*) spGroup

//  Relation to the former RBCELL FORTRAN subroutine:
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     subroutine  RBCELL ( celld,cvol )
//
//  ** the unit number iUnit and buffer for the return code iRet
//     have to be supplied.

FORTRAN_SUBR ( MMDB_F_RBCELL, mmdb_f_rbcell,
               (    // lengths-at-end list
                int     * iUnit,    // unit number; *iUnit<=0 means
                                    // "the last mentioned unit"
                mmdb::machine::apireal * celld,    // array to accept the cell parameters
                mmdb::machine::apireal * cvol,     // returns the cell volume
                int     * iRet      // return code
                        //   RWBERR_Ok  - success
                        //   RWBERR_NoChannel     if unit
                        //              iUnit was not
                        //              initialized
                        //   RWBERR_NoFile        if unit
                        //              has been disposed
                        //   RWBERR_Parameters    if the
                        //              cell parameters
                        //              were not set
                        //   RWBERR_NoOrthCode    if no
                                    //              orthogonalization
                        //              code was found
                        //   RWBERR_NoCheck       if check
                        //              of cell parameters
                        //              has failed.
                                    //   The last three returns would
                        // rather indicate a programming
                        // error in mmdb_rwbrook.cpp
               ), ( // lengths-in-structure list
                int     * iUnit,  mmdb::machine::apireal * celld,
                mmdb::machine::apireal * cvol,   int     * iRet
               ), ( // lengths-follow list
                int     * iUnit,  mmdb::machine::apireal * celld,
                mmdb::machine::apireal * cvol,   int     * iRet
               ) );



// ------------------------------------------------------------------

//   mmdb_f_rbcelln_(..) retrieves the crystal cell properties from the
// channel iUnit.

//  FORTRAN equivalent:   subroutine MMDB_F_RBCellN ( iUnit,celld,cvol,
//  ~~~~~~~~~~~~~~~~~~~                        ArgNCode,iRet )
//                        integer       iUnit,ArgNCode,iRet
//                        real          celld(6),cvol

//  Relation to the former RBCELL FORTRAN subroutine:
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     subroutine  RBCELL ( celld,cvol )
//
//  ** the unit number iUnit, buffer for orthogonalization code
//     ArgNCode and for the return code iRet have to be supplied.

FORTRAN_SUBR ( MMDB_F_RBCELLN, mmdb_f_rbcelln,
               (    // lengths-at-end list
                int     * iUnit,    // unit number; *iUnit<=0 means
                                    // "the last mentioned unit"
                mmdb::machine::apireal * celld,    // array to accept the cell parameters
                mmdb::machine::apireal * cvol,     // returns the cell volume
                int     * ArgNCode, // returns the orthogonalization code, 1-6
                int     * iRet      // return code
                                    //   RWBERR_Ok  - success
                                    //   RWBERR_NoChannel     if unit
                                    //              iUnit was not
                                    //              initialized
                                    //   RWBERR_NoFile        if unit
                                    //              has been disposed
                                    //   RWBERR_Parameters    if the
                                    //              cell parameters
                                    //              were not set
                                    //   RWBERR_NoOrthCode    if no
                                    //              orthogonalization
                                    //              code was found
                                    //   RWBERR_NoCheck       if check
                                    //              of cell parameters
                                    //              has failed.
                                    //   The last three returns would
                                    // rather indicate a programming
                                    // error in mmdb_rwbrook.cpp
               ), ( // lengths-in-structure list
                int * iUnit,
                mmdb::machine::apireal * celld,
                mmdb::machine::apireal * cvol,
                int * ArgNCode, int     * iRet
               ), ( // lengths-follow list
                int * iUnit,
                mmdb::machine::apireal * celld,
                mmdb::machine::apireal * cvol,
                int * ArgNCode, int     * iRet
               )
             );



// ------------------------------------------------------------------

//   mmdb_f_rbrcel_(..) retrieves the reciprocal cell dimensions and
// reciprocal cell volume from the channel iUnit.

//  FORTRAN equivalent:   subroutine MMDB_F_RBRCel ( iUnit,rcell,rvol,
//  ~~~~~~~~~~~~~~~~~~~                        iRet )
//                        integer       iUnit,iRet
//                        real          rcell(6),rvol

//  Relation to the former RBRCEL FORTRAN subroutine:
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     subroutine  RBRCEL ( rcell,rvol )
//
//  ** the unit number iUnit and buffer for the return code iRet
// have to be supplied.

FORTRAN_SUBR ( MMDB_F_RBRCEL, mmdb_f_rbrcel,
               (    // lengths-at-end list
                int     * iUnit,    // unit number
                mmdb::machine::apireal * rcell,    // array to accept the reciprocal
                                    // cell parameters
                mmdb::machine::apireal * rvol,     // returns the reciprocal cell volume
                int     * iRet      // return code
                        //   RWBERR_Ok  - success
                        //   RWBERR_NoChannel     if unit
                        //              iUnit was not
                        //              initialized
                        //   RWBERR_NoFile        if unit
                        //              has been disposed
                        //   RWBERR_Parameters    if the
                        //              cell parameters
                        //              were not set
                        //   RWBERR_NoOrthCode    if no
                                    //              orthogonalization
                        //              code was found
                        //   RWBERR_NoCheck       if check
                        //              of cell parameters
                        //              has failed.
                                    //   The last three returns would
                        // rather indicate a programming
                        // error in mmdb_rwbrook.cpp
               ), ( // lengths-in-structure list
                int * iUnit,
                mmdb::machine::apireal * rcell,
                mmdb::machine::apireal * rvol,
                int * iRet
               ), ( // lengths-follow list
                int * iUnit,
                mmdb::machine::apireal * rcell,
                mmdb::machine::apireal * rvol,
                int * iRet
               ) );


// ------------------------------------------------------------------

//   mmdb_f_rborf_(..) fills or retrieves the fractionalising (RF) and
// orthogonalising (RO) 4x4 matrices, as well as the orthogonalisation
// code (LCode) in/from unit iUnit.
//   If RO[1][1] (fortran notations, or RO[0] in C/C++) is set to 0.0
// then the matrices are retrieved and returned in RF and RO;
// otherwise RF and RO are stored in the unit.

//  FORTRAN equivalent:   subroutine MMDB_F_RBORF ( iUnit,RO,RF,
//  ~~~~~~~~~~~~~~~~~~~                         LCode,iRet )
//                        integer  iUnit,LCode,iRet
//                        real     RO(4,4),RF(4,4)

//  Relation to the former RBRORF2 FORTRAN subroutine:
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     subroutine  RBRORF2 ( RO,RF,LCode )
//
//  ** the unit number iUnit and buffer for the return code iRet
//     have to be supplied.

FORTRAN_SUBR ( MMDB_F_RBORF, mmdb_f_rborf,
               (     // lengths-at-end list
                 int     * iUnit, // unit number; *iUnit<=0 means
                                  // "the last mentioned unit"
                 mmdb::machine::apireal * RO,    // array for orthogonalising matrix
                 mmdb::machine::apireal * RF,    // array for fractionalising matrix
                 int     * LCode, // buffer for orthogonalisation code
                 int     * iRet   // return code:
                                  //   RWBERR_Ok  - success
                                  //   RWBERR_NoChannel     if unit
                                  //              iUnit was not
                                  //              initialized
                                  //   RWBERR_NoFile        if unit
                                  //              has been disposed
                                  //   RWBERR_NoMatrices    if the
                                  //              orthogonalisation
                                  //              matrices were not
                                  //              calculated
                ), ( // lengths-in-structure list
                 int * iUnit,
                 mmdb::machine::apireal * RO,
                 mmdb::machine::apireal * RF,
                 int * LCode, int * iRet
                ), ( // lengths-follow list
                 int * iUnit,
                 mmdb::machine::apireal * RO,
                 mmdb::machine::apireal * RF,
                 int * LCode, int * iRet
                )
             );



// ------------------------------------------------------------------

//   mmdb_f_orthmat_(..) calculates matrices for standard orthogonalisations
// and cell volume.
//   If Cell(1) is greater then zero, the existing cell parameters
// will be substituted. If new cell parameters differ substantially
// from the old ones, the returned value of Vol will be negative.

//  FORTRAN equivalent:   subroutine MMDB_F_OrthMat ( iUnit,Cell,Vol,
//  ~~~~~~~~~~~~~~~~~~~                        RRR,iRet
//                        integer  iUnit,iRet
//                        real     Cell(6),Vol,RRR(3,3,6)

//  Relation to the former RBFRO1 FORTRAN subroutine:
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     subroutine  RBFRO1 ( Cell,Vol,RRR )
//
//  ** the unit number iUnit and buffer for the return code iRet
//     have to be supplied.

FORTRAN_SUBR ( MMDB_F_ORTHMAT, mmdb_f_orthmat,
               (     // lengths-at-end list
                int     * iUnit, // unit number; *iUnit<=0 means
                                 // "the last mentioned unit"
                mmdb::machine::apireal * Cell,  // array of cell parameters:
                                 //  Cell(1) - a   Cell(4) - alpha
                                 //  Cell(2) - b   Cell(5) - beta
                                 //  Cell(3) - c   Cell(6) - gamma
                mmdb::machine::apireal * Vol,   // returns cell volume
                mmdb::machine::apireal * RRR,   // array (3,3,6), returns
                                 // orthogonalisation matrices
                int     * iRet   // return code:
                                 //   RWBERR_Ok  - success
                                 //   RWBERR_NoChannel     if unit
                                 //              iUnit was not
                                 //              initialized
                                 //   RWBERR_NoFile        if unit
                                 //              has been disposed
                                 //   RWBERR_NoMatrices    if the
                                 //              orthogonalisation
                                 //              matrices were not
                                 //              calculated
               ), ( // lengths-in-structure list
                int     * iUnit,
                mmdb::machine::apireal * Cell,
                mmdb::machine::apireal * Vol,
                mmdb::machine::apireal * RRR,   int * iRet
               ), ( // lengths-follow list
                int     * iUnit,
                mmdb::machine::apireal * Cell,
                mmdb::machine::apireal * Vol,
                mmdb::machine::apireal * RRR,   int * iRet
               )
             );



// ------------------------------------------------------------------

//   mmdb_f_cvanisou_(..) converts between crystallographic bs and
// orthogonal Us or the other way round.

//  FORTRAN equivalent:   subroutine MMDB_F_CVAnisou ( iUnit,U,iFlag,iRet )
//  ~~~~~~~~~~~~~~~~~~~   integer  iUnit,iFlag,iRet
//                        real     U(6)

//  Relation to the former CVANISOU FORTRAN subroutine:
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     subroutine  CVANISOU ( U,iFlag )
//
//  ** the unit number iUnit and buffer for the return code iRet
//     have to be supplied.

FORTRAN_SUBR ( MMDB_F_CVANISOU, mmdb_f_cvanisou,
               (     // lengths-at-end list
                int     * iUnit, // unit number; *iUnit<=0 means
                                 // "the last mentioned unit"
                mmdb::machine::apireal * U,     // array of coordinates to convert
                int     * iFlag, // =0: convert from fract. to orthog.
                                 // =1: convert from orthog. to fract.
                int     * iRet   // return code:
                                 //   RWBERR_Ok  - success
                                 //   RWBERR_NoChannel     if unit
                                 //              iUnit was not
                                 //              initialized
                                 //   RWBERR_NoFile        if unit
                                 //              has been disposed
                                 //   RWBERR_NoMatrices    if the
                                 //              orthogonalisation
                                 //              matrices were not
                                 //              calculated
               ), ( // lengths-in-structure list
                int * iUnit,
                mmdb::machine::apireal * U, int * iFlag, int * iRet
               ), ( // lengths-follow list
                int * iUnit,
                mmdb::machine::apireal * U, int * iFlag, int * iRet
               )
             );



// ------------------------------------------------------------------

//   mmdb_f_wremark_(..) writes a remark line into data structure.
// Generally, it puts the line on its place according to a PDB
// keyword which should start the line. The line will be always the
// last one in its group (last remark with given number or without
// it, last JRNL record, last ATOM etc.). If the keyword is not
// recognized, the line is appended after the coordinate section.
//   iRet will return same codes as those in mmdb_f_open1_(..) plus
// additional ones specified below.

//  FORTRAN equivalent:   subroutine MMDB_F_WRemark ( iUnit,Line,iRet )
//  ~~~~~~~~~~~~~~~~~~~   integer       iUnit,iRet
//                        character*(*) Line

//  Relation to the former WRemark FORTRAN subroutine:
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     subroutine  WRemark ( iUnit,Line )
//
//  ** the buffer for return code iRet has to be supplied

FORTRAN_SUBR ( MMDB_F_WREMARK, mmdb_f_wremark,
               (    // lengths-at-end list
                int     * iUnit, // unit number; *iUnit<=0 means
                                 // "the last mentioned unit"
                mmdb::machine::fpstr     Line,  // line to be added
                int     * iRet,  // return code:
                                 //   RWBERR_Ok  - success
                                 //   RWBERR_NoChannel     if unit
                                 //              iUnit was not
                                 //              initialized
                                 //   RWBERR_NoFile        if unit
                                 //              has been disposed
                                 // other return codea are those
                                 // returned by mmdb_f_open1_(..)
                int    Line_len  // fortran-hidden length of Line
               ), ( // lengths-in-structure list
                int * iUnit,
                mmdb::machine::fpstr Line, int * iRet
               ), ( // lengths-follow list
                int * iUnit,
                mmdb::machine::fpstr Line, int Line_len, int *iRet
               )
             );


/*
// ------------------------------------------------------------------

//   rbrinv_(..) takes 4x4 real matrix A and returns its inverse in
// matrix AI.

//  FORTRAN equivalent:   subroutine RBRInv ( A,AI )
//  ~~~~~~~~~~~~~~~~~~~   real  A(4,4),AI(4,4)

FORTRAN_SUBR ( RBRINV, rbrinv,
               ( mmdb::machine::apireal * A, mmdb::machine::apireal * AI ),
               ( mmdb::machine::apireal * A, mmdb::machine::apireal * AI ),
               ( mmdb::machine::apireal * A, mmdb::machine::apireal * AI )
             );

*/
/*

// ------------------------------------------------------------------

//   res3to1_(..) returns the 3-character or 1-character residue
// codes. One of them should be supplied (with the other one set
// blank), the routine returns the other one.

//  FORTRAN equivalent:   subroutine Res3to1 ( ResNm3,resNm1 )
//  ~~~~~~~~~~~~~~~~~~~   character*4 ResNm3
//                        cgaracter*1 ResNm1

FORTRAN_SUBR ( RES3TO1, res3to1,
               (     // lengths-at-end list
                mmdb::machine::fpstr ResNm3,   // 3-char name, 4th char
                                // will be set blank
                mmdb::machine::fpstr ResNm1,   // 1-char name
                int ResNm3_len, // fortran-hidden length of ResNm3
                int ResNm1_len  // fortran-hidden length of ResNm3
               ), ( // lengths-in-structure list
                mmdb::machine::fpstr ResNm3, mmdb::machine::fpstr ResNm1
               ), ( // lengths-follow list
                mmdb::machine::fpstr ResNm3, int ResNm3_len,
                mmdb::machine::fpstr ResNm1, int ResNm1_len
               )
             );

*/

// ------------------------------------------------------------------

//   rberrstop_(..) checks the return code got from one of the above
// functions, and if it indicates an error, it issues the following
// type of message (example)
//
//  *** RWBROOK error: point code unit    function
//  ***                   12   -3    3    MMDB_F_Open
//  *** file   : input.pdb
//  *** reason : cannot open a file
//  *** Execution stopped.
//
// if iStop is set to 0, and one of the following type
//
//  *** RWBROOK error: point code unit    function
//  ***                   12   -3    3    MMDB_F_Open
//  *** file   : input.pdb
//  *** reason : cannot open a file
//  *** continue running, may crash ...
//
// if iStop is not null.
//
//   iPlace (12 in the above samples) should be given a number which
// is unique through an application; it serves to the identifying
// the particular call which caused the problem. The return code
// (-3 in the above samples) is that what is back in the iRet
// parameter to the above functions. If iRet is set to RWBERR_Ok,
// rberrstop_(..) makes no action. If rberrstop_(..) is called
// immediately after a call to an RWBROOK function, e.g.
//
//    call MMDB_F_Open ( FName,RWStat,FType,iUnit,iRet )
//    call RBErrStop   ( 12,iRet,iUnit,0 )
//
// then the name of the misfunctioned call (MMDB_F_Open in the above
// samples) will be identified automatically and correctly.

//  FORTRAN equivalent:   subroutine RBErrStop ( iPlace,iRet,
//  ~~~~~~~~~~~~~~~~~~~                          iUnit ,iStop )
//                        integer  iUnit,iPlace,iRet,iStop

FORTRAN_SUBR ( RBERRSTOP, rberrstop,
               (   //    lengths-at-end list
                int * iPlace, // (unique) identificator inside an application
                int * iRet,   // return code to check
                int * iUnit,  // unit number
                int * iStop   // if 0 then stop if error
               ), ( // lengths-in-structure list
                int * iPlace, int * iRet,
                int * iUnit,  int * iStop
               ), ( // lengths-follow list
                int * iPlace, int * iRet,
                int * iUnit,  int * iStop
               ) );



// ------------------------------------------------------------------

//   rbcheckerr_(..) represents a simplified call to rberrstop_(..).
// It will work properly only if rbcheckerr_(..) is called
// immediately after an API function to be checked:
//
//    call MMDB_F_Open ( FName,RWStat,FType,iUnit,iRet )
//    call RBCheckErr  ( 12,0 )
//

//  FORTRAN equivalent:   subroutine RBCheckErr ( iPlace,iStop )
//  ~~~~~~~~~~~~~~~~~~~   integer  iPlace,iStop

FORTRAN_SUBR ( RBCHECKERR, rbcheckerr,
               (   //    lengths-at-end list
                int * iPlace, // (unique) identificator inside an application
                int * iStop   // if 0 then stop if error
               ), ( // lengths-in-structure list
                int * iPlace, int * iStop
               ), ( // lengths-follow list
                int * iPlace, int * iStop
               ) );


#endif
