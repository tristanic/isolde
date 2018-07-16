//  $Id: mmdb_uddata.h $
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
//    12.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :   MMDBF_UDData <interface>
//       ~~~~~~~~~
//  **** Project :   MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//
//  **** Classes :   mmdb::UDData ( user-defined data )
//       ~~~~~~~~~
//
//   (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#ifndef __MMDB_UDData__
#define __MMDB_UDData__

#include "mmdb_mask.h"
#include "imex.h"

namespace mmdb  {

  //  =======================  UDRegister  =========================

  enum UDR_TYPE  {
    UDR_ATOM      = 0,
    UDR_RESIDUE   = 1,
    UDR_CHAIN     = 2,
    UDR_MODEL     = 3,
    UDR_HIERARCHY = 4
  };

  enum UDD_FLAG  {
    UDRF_ATOM      = 0x01000000,
    UDRF_RESIDUE   = 0x02000000,
    UDRF_CHAIN     = 0x04000000,
    UDRF_MODEL     = 0x08000000,
    UDRF_HIERARCHY = 0x10000000,
    UDRF_MASK      = 0x00FFFFFF
  };

  DefineClass(UDRegister);
  DefineStreamFunctions(UDRegister);

  class MMDB_IMEX UDRegister : public io::Stream  {

    public :

      UDRegister ();
      UDRegister ( io::RPStream Object );
      ~UDRegister();

      int RegisterUDInteger ( UDR_TYPE udr_type, cpstr UDDataID );
      int RegisterUDReal    ( UDR_TYPE udr_type, cpstr UDDataID );
      int RegisterUDString  ( UDR_TYPE udr_type, cpstr UDDataID );
      int GetUDDHandle      ( UDR_TYPE udr_type, cpstr UDDataID );

      void write ( io::RFile f );
      void read  ( io::RFile f );

    protected :
      int      nIUDR[5],nRUDR[5],nSUDR[5];
      psvector IUDRegister[5];
      psvector RUDRegister[5];
      psvector SUDRegister[5];

      void  InitUDRegister ();
      void  FreeUDRegister ();
      int   RegisterUDData ( psvector & UDRegister,
                             int      & nUDR,
                             cpstr      UDDataID );

  };


  //  ==========================  UDData  ===========================

  enum UDDATA_RC  {
    UDDATA_Ok           =  0,
    UDDATA_WrongHandle  = -1,
    UDDATA_WrongUDRType = -2,
    UDDATA_NoData       = -3
  };

  DefineClass(UDData);
  DefineStreamFunctions(UDData);

  class MMDB_IMEX UDData : public Mask  {

    friend class SelManager;

    public :

      UDData ();
      UDData ( io::RPStream Object );
      ~UDData();

    protected :
      ivector  IUData;
      rvector  RUData;
      psvector SUData;

      void  InitUDData   ();
      void  FreeUDDMemory();
      int   getNofIUData ();
      int   getNofRUData ();
      int   getNofSUData ();
      void  setNofSUData ( int newN );

      int   putUDData ( int UDDhandle, int      iudd );
      int   putUDData ( int UDDhandle, realtype rudd );
      int   putUDData ( int UDDhandle, cpstr    sudd );

      int   getUDData ( int UDDhandle, int      & iudd );
      int   getUDData ( int UDDhandle, realtype & rudd );
      int   getUDData ( int UDDhandle, pstr sudd, int maxLen );
      pstr  getUDData ( int UDDhandle, int * retcode=NULL );
      int   getUDData ( int UDDhandle, pstr     & sudd );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

  };

}  // namespace mmdb

#endif

