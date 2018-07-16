/* mmdb_wrapper.cpp: MMDB_WRAPPER wrapper */
//C Copyright (C) 2000-2006 Kevin Cowtan and University of York
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA

#include "clipper_mmdb_wrapper.h"
#include <string.h>


namespace clipper {

namespace mmdb {

void XMMDBCryst::set_spacegroup( const Spacegroup& spacegroup )
{
  String name = spacegroup.descr().symbol_hm();
  strncpy( spaceGroup, name.substr( 0, 29 ).c_str(), 30 );
  SymOps.Reset();
  SymOps.PutGroupName( spaceGroup );
  for ( int i = 0; i < spacegroup.num_symops(); i++ )
    SymOps.AddSymOp( (char *)spacegroup.symop(i).format().c_str() );
  WhatIsSet |= CSET_SpaceGroup;
}

Spacegroup XMMDBManager::spacegroup()
{
  if ( isSpaceGroup() ) {  // get spacegroup from ops
    int nops = GetNumberOfSymOps();
    String ops = "";
    for ( int i = 0; i < nops; i++ ) ops += String( GetSymOp(i) ) + ";";
    return Spacegroup( Spgr_descr( ops, Spacegroup::Symops ) );
  } else {                 // otherwise get spacegroup from name
    String name = String( GetSpaceGroup() ).trim();
    if ( name.find_first_of( "PABCFIR" ) == String::npos ) name = "P 1";
    return Spacegroup( Spgr_descr( name, Spacegroup::HM ) );
  }
}

Cell XMMDBManager::cell()
{
  Cell_descr cd( 0, 0, 0 );
  if ( isCrystInfo() ) {
    cd = Cell_descr( Cryst.a, Cryst.b, Cryst.c,
		     Cryst.alpha, Cryst.beta, Cryst.gamma );
  }
  return Cell( cd );
}

void XMMDBManager::set_spacegroup( const Spacegroup& spacegroup )
{
  // convert to XMMDBCryst so we can set protected symops
  XMMDBCryst crys2;
  crys2.Copy( &Cryst );
  crys2.set_spacegroup( spacegroup );
  Cryst.Copy( &crys2 );
}

void XMMDBManager::set_cell( const Cell& cell )
{
  SetCell( cell.descr().a(), cell.descr().b(), cell.descr().c(),
	   cell.descr().alpha_deg(), cell.descr().beta_deg(),
	   cell.descr().gamma_deg() );
}

} // namespace mmdb

} // namespace clipper
