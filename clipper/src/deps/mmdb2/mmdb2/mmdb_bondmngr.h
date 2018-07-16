//  $Id: mmdb_bondmngr.h $
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
//  **** Module  :  mmdb_bondmngr <interface>
//       ~~~~~~~~~
//       Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::BondManager ( MMDB bonds maker )
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#ifndef __MMDB_BondMngr__
#define __MMDB_BondMngr__

#include "mmdb_selmngr.h"
#include "imex.h"

namespace mmdb  {

  // =======================  BondManager  =======================

  DefineClass(BondManager);
  DefineStreamFunctions(BondManager);

  class MMDB_IMEX BondManager : public SelManager  {

    public :

      BondManager ();
      BondManager ( io::RPStream Object );
      ~BondManager();

      void  MakeBonds  ( bool calc_only );
      void  RemoveBonds();

    protected :
      void  write ( io::RFile f );
      void  read  ( io::RFile f );

  };

}  // namespace mmdb

#endif

