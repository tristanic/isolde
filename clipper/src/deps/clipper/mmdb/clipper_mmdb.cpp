/* clipper_mmdb.cpp: MMDB wrapper */
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

#include "clipper_mmdb.h"

#include <algorithm>
extern "C" {
#include <string.h>
}

namespace clipper {


// MMDB wrapper types
// atom

String MMDBAtom::id() const
{ return String( name ); }

String MMDBAtom::element() const
{ return String( mmdb::CAtom::element ); }

Coord_orth MMDBAtom::coord_orth() const
{
  if ( !Ter && WhatIsSet & ::mmdb::ASET_Coordinates )
    return Coord_orth( x, y, z );
  else
    return Coord_orth( Coord_orth::null() );
}

ftype MMDBAtom::occupancy() const
{
  if ( !Ter && WhatIsSet & ::mmdb::ASET_Occupancy ) return mmdb::CAtom::occupancy;
  else                                              return Util::nan();
}

ftype MMDBAtom::u_iso() const
{
  if ( !Ter && WhatIsSet & ::mmdb::ASET_tempFactor ) return Util::b2u(tempFactor);
  else                                               return Util::nan();
}

U_aniso_orth MMDBAtom::u_aniso_orth() const
{
  if ( !Ter && WhatIsSet & ::mmdb::ASET_Anis_tFac )
    return U_aniso_orth( u11, u22, u33,
			 u12, u13, u23 );
  else
    return U_aniso_orth( U_aniso_orth::null() );
}

Sig_Coord_orth MMDBAtom::sig_coord_orth() const
{
  if ( !Ter && WhatIsSet & ::mmdb::ASET_CoordSigma )
    return Sig_Coord_orth( sigX, sigY, sigZ );
  else
    return Sig_Coord_orth( Sig_Coord_orth::null() );
}

ftype MMDBAtom::sig_occupancy() const
{
  if ( !Ter && WhatIsSet & ::mmdb::ASET_OccSigma ) return sigOcc;
  else                                     return Util::nan();
}

ftype MMDBAtom::sig_u_iso() const
{
  if ( !Ter && WhatIsSet & ::mmdb::ASET_tFacSigma ) return Util::b2u(sigTemp);
  else                                      return Util::nan();
}

Sig_U_aniso_orth MMDBAtom::sig_u_aniso_orth() const
{
  if ( !Ter && WhatIsSet & ::mmdb::ASET_Anis_tFSigma )
    return Sig_U_aniso_orth( su11, su22, su33,
			     su12, su13, su23 );
  else
    return Sig_U_aniso_orth( Sig_U_aniso_orth::null() );
}

void MMDBAtom::set_id( const String& n )
{ this->SetAtomName( (char *)n.c_str() ); }

void MMDBAtom::set_element( const String& n )
{ SetElementName( (char *)n.c_str() ); }

void MMDBAtom::set_coord_orth( const Coord_orth& v )
{
  WhatIsSet &= ~::mmdb::ASET_Coordinates;
  if ( !v.is_null() ) {
    x = v.x(); y = v.y(); z = v.z();
    WhatIsSet |= ::mmdb::ASET_Coordinates;
  }
}

void MMDBAtom::set_occupancy( const ftype& v )
{ 
  WhatIsSet &= ~::mmdb::ASET_Occupancy;
  if ( !Util::is_nan( v ) ) {
    mmdb::CAtom::occupancy = v;
    WhatIsSet |= ::mmdb::ASET_Occupancy;
  }
}

void MMDBAtom::set_u_iso( const ftype& v )
{ 
  WhatIsSet &= ~::mmdb::ASET_tempFactor;
  if ( !Util::is_nan( v ) ) {
    tempFactor = Util::u2b(v);
    WhatIsSet |= ::mmdb::ASET_tempFactor;
  }
}

void MMDBAtom::set_u_aniso_orth( const U_aniso_orth& v )
{
  WhatIsSet &= ~::mmdb::ASET_Anis_tFac;
  if ( !v.is_null() ) {
    u11 = v(0,0); u22 = v(1,1); u33 = v(2,2);
    u12 = v(0,1); u13 = v(0,2); u23 = v(1,2);
    WhatIsSet |= ::mmdb::ASET_Anis_tFac;
  }
}

void MMDBAtom::set_sig_coord_orth( const Sig_Coord_orth& s )
{
  WhatIsSet &= ~::mmdb::ASET_CoordSigma;
  if ( !s.is_null() ) {
    sigX = s.sigx(); sigY = s.sigy(); sigZ = s.sigz();
    WhatIsSet |= ::mmdb::ASET_CoordSigma;
  }
}

void MMDBAtom::set_sig_occupancy( const ftype& s )
{ 
  WhatIsSet &= ~::mmdb::ASET_OccSigma;
  if ( !Util::is_nan( s ) ) {
    sigOcc = s;
    WhatIsSet |= ::mmdb::ASET_OccSigma;
  }
}

void MMDBAtom::set_sig_u_iso( const ftype& s )
{ 
  WhatIsSet &= ~::mmdb::ASET_tFacSigma;
  if ( !Util::is_nan( s ) ) {
    sigTemp = Util::u2b(s);
    WhatIsSet |= ::mmdb::ASET_tFacSigma;
  }
}

void MMDBAtom::set_sig_u_aniso_orth( const Sig_U_aniso_orth& s )
{
  WhatIsSet &= ~::mmdb::ASET_Anis_tFSigma;
  if ( !s.is_null() ) {
    su11 = s(0,0); su22 = s(1,1); su33 = s(2,2);
    su12 = s(0,1); su13 = s(0,2); su23 = s(1,2);
    WhatIsSet |= ::mmdb::ASET_Anis_tFSigma;
  }
}

/*! \return The atom alternate conformation code. */
String MMDBAtom::altconf() const
{ return String( altLoc ); }

/*! \return The atom serial number. */
int MMDBAtom::serial_num() const
{ return serNum; }

/*! \return The atomic charge. */
String MMDBAtom::charge() const
{ return String( mmdb::CAtom::charge ); }


// residue

String MMDBResidue::type() const
{ return String( name ); }

int MMDBResidue::seqnum() const
{ return seqNum; }

String MMDBResidue::inscode() const
{ return String( insCode ); }

void MMDBResidue::set_type( const String& n )
{ SetResName( (char *)n.c_str() ); }

void MMDBResidue::set_seqnum( const int& n )
{ seqNum = n; }

void MMDBResidue::set_inscode( const String& n )
{ strncpy( insCode, n.c_str(), 10 ); }


// chain

String MMDBChain::id() const
{ return String( ((mmdb::CChain*)this)->GetChainID() ); }

void MMDBChain::set_id( const String& n )
{ SetChainID( (char *)n.c_str() ); }


// model

String MMDBModel::id() const
{ return String( ((mmdb::CModel*)this)->GetEntryID() ); }

void MMDBModel::set_id( const String& n )
{ SetEntryID( (char *)n.c_str() ); }


// MMDBManager methods


/*! For later initialisation: see init() */
MMDBManager::MMDBManager() 
{ ::mmdb::InitMatType(); Message::message( Message_ctor( " [MMDBManager: constructed>" ) ); }

MMDBManager::~MMDBManager()
{ Message::message( Message_dtor( " <MMDBManager: destroyed]" ) ); }

Spacegroup MMDBManager::spacegroup() const
{
  mmdb::CMMDBManager& mmdb = const_cast<MMDBManager&>(*this);
  Spacegroup sg;
  if ( mmdb.isSpaceGroup() ) {  // get spacegroup from ops
    int nops = mmdb.GetNumberOfSymOps();
    String ops = "";
    for ( int i = 0; i < nops; i++ ) ops += String( mmdb.GetSymOp(i) ) + ";";
    Spgr_descr sd( ops, Spacegroup::Symops );
    sg = Spacegroup( sd );
  } else {                 // otherwise get spacegroup from name
    const char* mmdbsg = mmdb.GetSpaceGroup();
    if ( mmdbsg != NULL ) {
      String name = String( mmdbsg ).trim();
      if ( name.find_first_of( "PABCFIR" ) != String::npos ) {
	if ( name.find_first_of( ":" ) != String::npos ) {
	  Spgr_descr sd( name, Spacegroup::XHM );
	  sg = Spacegroup( sd );
	} else {
	  Spgr_descr sd( name, Spacegroup::HM );
	  sg = Spacegroup( sd );
	}
      }
    }
  }
  return sg;
}

Cell MMDBManager::cell() const
{
  mmdb::CMMDBManager& mmdb = const_cast<MMDBManager&>(*this);
  if ( mmdb.isCrystInfo() ) {
    return Cell( Cell_descr( this->cryst.a, this->cryst.b, this->cryst.c,
			     this->cryst.alpha, this->cryst.beta, this->cryst.gamma ) );
  } else {
    return Cell();  // null
  }
}

void MMDBManager::set_spacegroup( const Spacegroup& spacegroup )
{
  // Hall method... (gives bad symbol in file)
  //String hall = "Hall:" + spacegroup.symbol_hall();
  //SetSpaceGroup( (char*)hall.c_str() );
  // HM method...
  //String s = spacegroup.symbol_hm();
  //String x = spacegroup.symbol_hm_ext();
  //if ( x.length() > 0 ) s = x + s.substr(1);
  //SetSpaceGroup( (char*)s.c_str() );
  // xHM method
  SetSpaceGroup( (char*)spacegroup.symbol_xhm().c_str() );
}

void MMDBManager::set_cell( const Cell& cell )
{
  SetCell( cell.a(), cell.b(), cell.c(),
	   cell.alpha_deg(), cell.beta_deg(), cell.gamma_deg() );
}


/*! The atom list is constructed from the MMDB atom pointers by
  copying. This should be fast compared to any task for which Clipper
  uses an atom list.
  \param ppcatom The MMDB array of atom pointers. 
  \param natom The number of atoms. */
MMDBAtom_list::MMDBAtom_list( const mmdb::PPCAtom ppcatom, const int natom )
{
  for ( int i = 0; i < natom; i++ )
    push_back( Atom( *( (const MMDBAtom*)ppcatom[i] ) ) );
}


} // namespace clipper
