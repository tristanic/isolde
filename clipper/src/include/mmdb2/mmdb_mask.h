//  $Id: mmdb_mask.h $
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
//  **** Module  :   MMDBF_Mask <interface>
//       ~~~~~~~~~
//  **** Project :   MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//
//  **** Classes :   mmdb::Mask  ( atom selection mask )
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#ifndef __MMDB_Mask__
#define __MMDB_Mask__

#include "mmdb_io_stream.h"
#include "imex.h"

namespace mmdb  {

  //  ==========================  Mask  =============================

  DefineClass(Mask);
  DefineStreamFunctions(Mask);

  class MMDB_IMEX Mask : public io::Stream  {

    public :

      Mask ();
      Mask ( io::RPStream Object );
      ~Mask();

      void SetMaskBit ( int  BitNo );
      void NewMask    ( PPMask Mask, int nMasks );

      void CopyMask   ( PMask Mask );   //  this = Mask
      void SetMask    ( PMask Mask );   //  this = this | Mask
      void RemoveMask ( PMask Mask );   //  this = this & (~Mask)
      void SelMask    ( PMask Mask );   //  this = this & Mask
      void XadMask    ( PMask Mask );   //  this = this ^ Mask
      void ClearMask  ();               //  this = NULL
      void NegMask    ();               //  this = ~this

      bool CheckMask ( PMask Mask );    // true if the bit is on
      bool isMask    ();                // true if any mask bit is on

      inline int getLength() { return mlen; }

      pstr Print ( pstr S ); // returns binary string

      void write ( io::RFile f );
      void read  ( io::RFile f );

    protected :
      int     mlen;
      wvector m;

      void InitMask();
      void Expand  ( int n );

  };

}  // namespace mmdb

#endif

