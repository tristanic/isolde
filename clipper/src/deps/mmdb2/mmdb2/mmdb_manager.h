//  $Id: mmdb_manager.h $
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
//    15.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  mmdb_manager <interface>
//       ~~~~~~~~~
//       Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::Manager  ( MMDB file manager )
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#ifndef __MMDB_Manager__
#define __MMDB_Manager__

#include "mmdb_bondmngr.h"
#include "imex.h"

namespace mmdb  {

  // =======================  Manager  ===========================

  // copy masks
  enum COPY_MASK  {
    MMDBFCM_None        = 0x00000000,
    MMDBFCM_All         = 0xFFFFFFFF,
    MMDBFCM_Title       = 0x00000001,
    MMDBFCM_TitleKeepBM = 0x00000002,
    MMDBFCM_Cryst       = 0x00000004,
    MMDBFCM_Coord       = 0x00000008,
    MMDBFCM_SecStruct   = 0x00000010,
    MMDBFCM_HetInfo     = 0x00000020,
    MMDBFCM_Links       = 0x00000040,
    MMDBFCM_CisPeps     = 0x00000080,
    MMDBFCM_SA          = 0x00000100,
    MMDBFCM_SB          = 0x00000200,
    MMDBFCM_SC          = 0x00000400,
    MMDBFCM_Footnotes   = 0x00000800,
    MMDBFCM_ChainAnnot  = 0x00001000,
    MMDBFCM_Flags       = 0x00002000,
    MMDBFCM_Buffer      = 0x80000000,
    MMDBFCM_Top         = 0xFFFFFFF7
  };

  DefineStreamFunctions(Manager);

  class MMDB_IMEX Manager : public BondManager  {

    public :

      Manager ();
      Manager ( io::RPStream Object );
      ~Manager();


      //  ---------------  Copying/Deleting  -----------------------

      //   Copy(..) will transfer different sort of information
      // between two MMDB's according to the copy mask given
      // (cf. MMDBFCM_XXXXX values). Note that the copying content
      // replaces the corresponding information (e.g. copying
      // coordinates will replace existing coordinates rather than
      // add to them).
      void  Copy   ( PManager MMDB, COPY_MASK CopyMask );

      //   Delete(..) deletes different sort of information from
      // the MMDB according to the delete mask given.
      void  Delete ( word DelMask );  // DelMask is the same as CopyMask

      PTitleContainer GetRemarks();
      PTitleContainer GetJournal();

      realtype GetResolution(); // -1.0 means no resolution record in file

      int   ParseBiomolecules(); // returns the number of biomolecules,
                                 // -2 for general format error
                                 // -3 for errors in BIOMT records
      int   GetNofBiomolecules();
      void  GetBiomolecules   ( PPBiomolecule & BM, int & nBMs );

      PBiomolecule GetBiomolecule ( int bmNo ); // bmno=0,1,..
                                 // returns NULL if bmNo is incorrect
      PManager MakeBiomolecule ( int bmNo, int modelNo=1 );


    protected :

      //  ---------------  Stream I/O  -----------------------------
      void  write  ( io::RFile f );
      void  read   ( io::RFile f );

  };

}  // namespace mmdb

#endif

