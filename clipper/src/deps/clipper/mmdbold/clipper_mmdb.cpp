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

#include <string.h>
#include <set>
#include <algorithm>


namespace clipper {


Message_ctor message_ctor_mmdb( " [MMDB: constructed>" );
Message_dtor message_dtor_mmdb( " <MMDB: destroyed]" );


// MMDB wrapper types
// atom

/*! \return \c true if uninitialised */
bool DBAtom::is_null() const
{ return ptr == NULL; }

/*! \return \c true if not null and not "TER" */
bool DBAtom::is_atom() const
{ if ( is_null() ) return false; return !(ptr->Ter); }

/*! \return The parent residue of this atom (check result.is_null()) */
DBResidue DBAtom::residue() const
{ return DBResidue( ptr->GetResidue() ); }

/*! Set all the hierarchy independent properties of this atom from
  another (DB or non-DB) atom.
  \param a The source atom. */
void DBAtom::set( const DBAtom_base& a )
{
  set_type( a.type() );
  set_element( a.element() );
  set_altconf( a.altconf() );
  set_coord_orth( a.coord_orth() );
  set_occupancy( a.occupancy() );
  set_u_iso( a.u_iso() );
  set_u_aniso_orth( a.u_aniso_orth() );
}

int DBAtom::index() const
{
  DBResidue m = residue();
  int natm;
  for ( natm = 0; natm < m.size(); natm++ )
    if ( ptr == m[natm].pcatom() ) break;
  return natm;
}

String DBAtom::type() const
{ return String( ptr->name ); }

String DBAtom::element() const
{
  if ( is_atom() ) return String( ptr->element );
  return "";
}

String DBAtom::altconf() const
{ return String( ptr->altLoc ); }

Coord_orth DBAtom::coord_orth() const
{
  if ( is_atom() ) if ( ptr->WhatIsSet & ASET_Coordinates )
    return Coord_orth( ptr->x, ptr->y, ptr->z );
  return Coord_orth( Coord_orth::null() );
}

ftype DBAtom::occupancy() const
{
  if ( is_atom() ) if ( ptr->WhatIsSet & ASET_Occupancy )
      return ptr->occupancy;
  return Util::nan();
}

ftype DBAtom::u_iso() const
{
  if ( is_atom() ) if ( ptr->WhatIsSet & ASET_tempFactor )
    return Util::b2u(ptr->tempFactor);
  return Util::nan();
}

U_aniso_orth DBAtom::u_aniso_orth() const
{
  if ( is_atom() ) if ( ptr->WhatIsSet & ASET_Anis_tFac )
    return U_aniso_orth( ptr->u11, ptr->u22, ptr->u33,
			 ptr->u12, ptr->u13, ptr->u23 );
  return U_aniso_orth( U_aniso_orth::null() );
}

Sig_Coord_orth DBAtom::sig_coord_orth() const
{
  if ( ptr->WhatIsSet & ASET_CoordSigma )
    return Sig_Coord_orth( ptr->sigX, ptr->sigY, ptr->sigZ );
  else
    return Sig_Coord_orth( Sig_Coord_orth::null() );
}

ftype DBAtom::sig_occupancy() const
{
  if ( ptr->WhatIsSet & ASET_OccSigma ) return ptr->sigOcc;
  else                                  return Util::nan();
}

ftype DBAtom::sig_u_iso() const
{
  if ( ptr->WhatIsSet & ASET_tFacSigma ) return Util::b2u(ptr->sigTemp);
  else                                   return Util::nan();
}

Sig_U_aniso_orth DBAtom::sig_u_aniso_orth() const
{
  if ( ptr->WhatIsSet & ASET_Anis_tFSigma )
    return Sig_U_aniso_orth( ptr->su11, ptr->su22, ptr->su33,
			ptr->su12, ptr->su13, ptr->su23 );
  else
    return Sig_U_aniso_orth( Sig_U_aniso_orth::null() );
}

void DBAtom::set_type( const String& n )
{ ptr->SetAtomName( (char *)n.c_str() ); }

void DBAtom::set_element( const String& n )
{ ptr->SetElementName( (char *)n.c_str() ); }

void DBAtom::set_altconf( const String& n )
{ /*strncpy( ptr->altLoc, n.c_str(), 20 );*/ }

void DBAtom::set_coord_orth( const Coord_orth& v )
{
  ptr->WhatIsSet &= ~ASET_Coordinates;
  if ( !v.is_null() ) {
    ptr->x = v.x(); ptr->y = v.y(); ptr->z = v.z();
    ptr->WhatIsSet |= ASET_Coordinates;
  }
}

void DBAtom::set_occupancy( const ftype& v )
{ 
  ptr->WhatIsSet &= ~ASET_Occupancy;
  if ( !Util::is_nan( v ) ) {
    ptr->occupancy = v;
    ptr->WhatIsSet |= ASET_Occupancy;
  }
}

void DBAtom::set_u_iso( const ftype& v )
{ 
  ptr->WhatIsSet &= ~ASET_tempFactor;
  if ( !Util::is_nan( v ) ) {
    ptr->tempFactor = Util::u2b(v);
    ptr->WhatIsSet |= ASET_tempFactor;
  }
}

void DBAtom::set_u_aniso_orth( const U_aniso_orth& v )
{
  ptr->WhatIsSet &= ~ASET_Anis_tFac;
  if ( !v.is_null() ) {
    ptr->u11 = v(0,0); ptr->u22 = v(1,1); ptr->u33 = v(2,2);
    ptr->u12 = v(0,1); ptr->u13 = v(0,2); ptr->u23 = v(1,2);
    ptr->WhatIsSet |= ASET_Anis_tFac;
  }
}

void DBAtom::set_sig_coord_orth( const Sig_Coord_orth& s )
{
  ptr->WhatIsSet &= ~ASET_CoordSigma;
  if ( !s.is_null() ) {
    ptr->sigX = s.sigx(); ptr->sigY = s.sigy(); ptr->sigZ = s.sigz();
    ptr->WhatIsSet |= ASET_CoordSigma;
  }
}

void DBAtom::set_sig_occupancy( const ftype& s )
{ 
  ptr->WhatIsSet &= ~ASET_OccSigma;
  if ( !Util::is_nan( s ) ) {
    ptr->sigOcc = s;
    ptr->WhatIsSet |= ASET_OccSigma;
  }
}

void DBAtom::set_sig_u_iso( const ftype& s )
{ 
  ptr->WhatIsSet &= ~ASET_tFacSigma;
  if ( !Util::is_nan( s ) ) {
    ptr->sigTemp = Util::u2b(s);
    ptr->WhatIsSet |= ASET_tFacSigma;
  }
}

void DBAtom::set_sig_u_aniso_orth( const Sig_U_aniso_orth& s )
{
  ptr->WhatIsSet &= ~ASET_Anis_tFSigma;
  if ( !s.is_null() ) {
    ptr->su11 = s(0,0); ptr->su22 = s(1,1); ptr->su33 = s(2,2);
    ptr->su12 = s(0,1); ptr->su13 = s(0,2); ptr->su23 = s(1,2);
    ptr->WhatIsSet |= ASET_Anis_tFSigma;
  }
}

/*! \return The atom serial number. */
int DBAtom::serial_num() const
{ return ptr->serNum; }

/*! \return The atomic charge. */
String DBAtom::charge() const
{ return String( ptr->charge ); }


// residue

bool DBResidue::is_null() const
{ return ptr == NULL; }

DBChain DBResidue::chain() const
{ return DBChain( ptr->GetChain() ); }

void DBResidue::set( const DBResidue_base& r )
{
  set_type( r.type() );
  set_seqnum( r.seqnum() );
  set_inscode( r.inscode() );
}

int DBResidue::index() const
{
  DBChain m = chain();
  int nres;
  for ( nres = 0; nres < m.size(); nres++ )
    if ( ptr == m[nres].pcresidue() ) break;
  return nres;
}

String DBResidue::type() const
{ return String( ptr->GetResName() ); }

int DBResidue::seqnum() const
{ return ptr->GetSeqNum(); }

String DBResidue::inscode() const
{ return String( ptr->GetResName() ); }

void DBResidue::set_type( const String& n )
{ ptr->SetResName( (char *)n.c_str() ); }

void DBResidue::set_seqnum( const int& n )
{ ptr->seqNum = n; }

void DBResidue::set_inscode( const String& n )
{ strncpy( ptr->insCode, n.c_str(), 10 ); }

/*! See MMDB::finalise_edit(). */
DBAtom DBResidue::add_atom( const DBAtom_base& a )
{
  PCAtom p = new CAtom();
  DBAtom atom( p );
  atom.set( a );
  PCMMDBCoorManager pmmdb = (PCMMDBCoorManager)(ptr->GetCoordHierarchy());
  if ( pmmdb->AddAtom( chain().model().index()+1, chain().index(), index(), p )
       <= 0 )
    Message::message( Message_fatal( "DBResidue: attempt to add atom already in residue" ) );
  //((PCMMDBFile)ptr->GetCoordHierarchy())->FinishStructEdit();
  return atom;
}

DBAtom DBResidue::atom( const String& type ) const
{
  DBAtom a;
  for ( int i = 0; i < size(); i++ ) {
    a = (*this)[i];
    if ( a.is_atom() )
      if ( a.type().find( type ) != String::npos ) return a;
  }
  return DBAtom();
}

int DBResidue::size() const
{
  if ( ptr != NULL ) return ptr->GetNumberOfAtoms();
  else               return 0;
}

DBAtom DBResidue::operator[] ( const int& i ) const
{ return DBAtom( ptr->GetAtom(i) ); }


// chain

bool DBChain::is_null() const
{ return ptr == NULL; }

void DBChain::set( const DBChain_base& c )
{ set_id( c.id() ); }

int DBChain::index() const
{
  DBModel m = model();
  int nchn;
  for ( nchn = 0; nchn < m.size(); nchn++ )
    if ( ptr == m[nchn].pcchain() ) break;
  return nchn;
}

String DBChain::id() const
{ return String( ptr->GetChainID() ); }

void DBChain::set_id( const String& n )
{ ptr->SetChainID( (char *)n.c_str() ); }

DBModel DBChain::model() const
{ return DBModel( ptr->GetModel() ); }

/*! See MMDB::finalise_edit(). */
DBResidue DBChain::add_residue( const DBResidue_base& r )
{
  PCResidue p = new CResidue();
  DBResidue res( p );
  res.set( r );
  PCMMDBCoorManager pmmdb = (PCMMDBCoorManager)(ptr->GetCoordHierarchy());
  if ( pmmdb->AddResidue( model().index()+1, index(), p ) <= 0 )
    Message::message( Message_fatal( "DBChain: attempt to add residue already in chain" ) );
  //((PCMMDBFile)ptr->GetCoordHierarchy())->FinishStructEdit();
  return res;
}

DBResidue DBChain::residue( const int& seqnum ) const
{
  DBResidue r0, r1;
  // first try and guess where to find the residue
  if ( size() > 0 ) {
    r0 = (*this)[0];
    if ( !r0.is_null() ) {
      int i = seqnum - r0.seqnum();
      if ( i >= 0 && i < size() ) {
	r1 = (*this)[i];
	if ( !r1.is_null() )
	  if ( r1.seqnum() == seqnum ) return r1;
      }
    }
  }
  // if that fails, do binary slice (FIXME)
  // if that fails, search for it
  for ( int i = 0; i < size(); i++ ) {
    r1 = (*this)[i];
    if ( !r1.is_null() )
      if ( r1.seqnum() == seqnum ) return r1;
  }
  return DBResidue();
}

int DBChain::size() const
{
  if ( ptr != NULL ) return ptr->GetNumberOfResidues();
  else               return 0;
}

DBResidue DBChain::operator[] ( const int& i ) const
{ return DBResidue( ptr->GetResidue(i) ); }


// model

bool DBModel::is_null() const
{ return ptr == NULL; }

void DBModel::set( const DBModel_base& m )
{ set_id( m.id() ); }

int DBModel::index() const
{ return ptr->GetSerNum()-1; }

String DBModel::id() const
{ return String( ptr->GetEntryID() ); }

void DBModel::set_id( const String& n )
{ ptr->SetEntryID( (char *)n.c_str() ); }

DBManager DBModel::manager() const
{ return DBManager( (mmdb::XMMDBManager*)ptr->GetCoordHierarchy() ); }

/*! See MMDB::finalise_edit(). */
DBChain DBModel::add_chain( const DBChain_base& r )
{
  PCChain p = new CChain();
  DBChain chn( p );
  chn.set( r );
  //int nmdl = ptr->GetSerNum();
  PCMMDBCoorManager pmmdb = (PCMMDBCoorManager)(ptr->GetCoordHierarchy());
  if ( pmmdb->AddChain( index()+1, p ) <= 0 )
    Message::message( Message_fatal( "DBModel: attempt to add chain already in model" ) );
  return chn;
}

DBChain DBModel::chain( const String& id ) const
{
  DBChain c;
  for ( int i = 0; i < size(); i++ ) {
    c = (*this)[i];
    if ( !c.is_null() )
      if ( c.id() == id ) return c;
  }
  return DBChain();
}

int DBModel::size() const
{
  if ( ptr != NULL ) return ptr->GetNumberOfChains();
  else               return 0;
}

DBChain DBModel::operator[] ( const int& i ) const
{ return DBChain( ptr->GetChain(i) ); }


// selection types

DBAtom_selection::DBAtom_selection( mmdb::PPCAtom p, int n ) : list( p, p+n )
{}

void DBAtom_selection::add_atom( DBAtom a )
{ list.push_back( a.pcatom() ); }

DBAtom_selection operator& ( const DBAtom_selection& a1, const DBAtom_selection& a2 )
{
  std::set<mmdb::PCAtom> s1( a1.list.begin(), a1.list.end() );
  std::set<mmdb::PCAtom> s2( a2.list.begin(), a2.list.end() );
  std::vector<mmdb::PCAtom> v;
  std::set_intersection( s1.begin(), s1.end(), s2.begin(), s2.end(),
			 std::back_inserter( v ) );
  return DBAtom_selection( &(v[0]), v.size() );
}

DBAtom_selection operator| ( const DBAtom_selection& a1, const DBAtom_selection& a2 )
{
  std::set<mmdb::PCAtom> s1( a1.list.begin(), a1.list.end() );
  std::set<mmdb::PCAtom> s2( a2.list.begin(), a2.list.end() );
  std::vector<mmdb::PCAtom> v;
  std::set_union( s1.begin(), s1.end(), s2.begin(), s2.end(),
		  std::back_inserter( v ) );
  return DBAtom_selection( &(v[0]), v.size() );
}

DBAtom_selection operator^ ( const DBAtom_selection& a1, const DBAtom_selection& a2 )
{
  std::set<mmdb::PCAtom> s1( a1.list.begin(), a1.list.end() );
  std::set<mmdb::PCAtom> s2( a2.list.begin(), a2.list.end() );
  std::vector<mmdb::PCAtom> v;
  std::set_symmetric_difference( s1.begin(), s1.end(), s2.begin(), s2.end(),
				 std::back_inserter( v ) );
  return DBAtom_selection( &(v[0]), v.size() );
}

DBAtom_selection_inv operator! ( const DBAtom_selection& a )
{ return DBAtom_selection_inv( a ); }


DBResidue_selection::DBResidue_selection( mmdb::PPCResidue p, int n ) : list( p, p+n )
{}

void DBResidue_selection::add_residue( DBResidue a )
{ list.push_back( a.pcresidue() ); }

DBResidue_selection operator& ( const DBResidue_selection& a1, const DBResidue_selection& a2 )
{
  std::set<mmdb::PCResidue> s1( a1.list.begin(), a1.list.end() );
  std::set<mmdb::PCResidue> s2( a2.list.begin(), a2.list.end() );
  std::vector<mmdb::PCResidue> v;
  std::set_intersection( s1.begin(), s1.end(), s2.begin(), s2.end(),
			 std::back_inserter( v ) );
  return DBResidue_selection( &(v[0]), v.size() );
}

DBResidue_selection operator| ( const DBResidue_selection& a1, const DBResidue_selection& a2 )
{
  std::set<mmdb::PCResidue> s1( a1.list.begin(), a1.list.end() );
  std::set<mmdb::PCResidue> s2( a2.list.begin(), a2.list.end() );
  std::vector<mmdb::PCResidue> v;
  std::set_union( s1.begin(), s1.end(), s2.begin(), s2.end(),
		  std::back_inserter( v ) );
  return DBResidue_selection( &(v[0]), v.size() );
}

DBResidue_selection operator^ ( const DBResidue_selection& a1, const DBResidue_selection& a2 )
{
  std::set<mmdb::PCResidue> s1( a1.list.begin(), a1.list.end() );
  std::set<mmdb::PCResidue> s2( a2.list.begin(), a2.list.end() );
  std::vector<mmdb::PCResidue> v;
  std::set_symmetric_difference( s1.begin(), s1.end(), s2.begin(), s2.end(),
				 std::back_inserter( v ) );
  return DBResidue_selection( &(v[0]), v.size() );
}

DBResidue_selection_inv operator! ( const DBResidue_selection& a )
{ return DBResidue_selection_inv( a ); }


DBChain_selection::DBChain_selection( mmdb::PPCChain p, int n ) : list( p, p+n )
{}

void DBChain_selection::add_chain( DBChain a )
{ list.push_back( a.pcchain() ); }


DBModel_selection::DBModel_selection( mmdb::PPCModel p, int n ) : list( p, p+n )
{}

void DBModel_selection::add_model( DBModel a )
{ list.push_back( a.pcmodel() ); }



// DBManager methods

/*! Fetch the first (usually the only) model, or the N'th model
  indexed from 1 if an argument is supplied.
  \param i The number of the model to return.
  \return The requested model, or a null model if i does not exist. */
DBModel DBManager::model( const int i )
{ return DBModel( ptr->GetModel(i) ); }

int DBManager::size() const
{
  if ( ptr != NULL ) return ptr->GetNumberOfModels();
  else               return 0;
}

/*! The new model becomes the last in the list. Usually you would just
  use this to add the initial model when starting from a blank db.
  See finalise_edit(). */
DBModel DBManager::add_model( const DBModel_base& m )
{
  PCModel p = new CModel();
  DBModel mod( p );
  mod.set( m );
  if ( ptr->AddModel( p ) <= 0 )
    Message::message( Message_fatal( "DBModel: attempt to add chain already in model" ) );
  return mod;
}

/*! Update the hierarchy after adding or removing objects. This
  invalidates all pointer and all DBModel/DBChain/DBResidue/DBAtom
  objects. */
void DBManager::finalise_edit()
{
  ptr->FinishStructEdit();
}

/*! Fetch the model with a given index, based fromzero. Note that MMDB
  models are indexed from 1, so this method differs by 1 in numbering
  from the model ID you will need for a search function. See
  MMDB::model().
  \param i The index of the model to return.
  \return The requested model, or a null model if i does not exist. */
DBModel DBManager::operator[] ( const int& i ) const
{ return DBModel( ptr->GetModel(i+1) ); }


/*! Select atoms using an MMDB coordinate ID expression.
  \param s The coordinate ID expression.
  \return The atom selection object.
  \par Examples:
  <table>
  <tr><td> \b Coordinate_ID_expression <td> \b Meaning
  <tr><td> \c *                 <td> all atoms in all models/chains/residues
  <tr><td> \c [C]               <td> all carbons in all models/chains/residues
  <tr><td> \c /1///:A           <td> all atoms in alternate location A in all residues of chain without a chain ID, in model 1
  <tr><td> \c 30-120            <td> all atoms in residues 30 to 120 in all models and chains.
  <tr><td> \c A/30.A-120.S      <td> all atoms in residues 30 insertion code A to 120 insertion code S in A-chains of all models.
  <tr><td> \c A/30.A-120.S/[!S] <td> all atoms but sulphur in residues 30 insertion code A to 120 insertion code S in A-chains of all models.
  <tr><td> \c (ALA,SER)         <td> all atoms in residues ALA and SER in all models/chains.
  <tr><td> \c /1/A/(!ALA,SER)/CA[C] <td> all C-alpha atoms in all residues but ALA and SER in model 1 chain A.
  </table> */
DBAtom_selection DBManager::select_atoms( const String& s )
{
  mmdb::PPCAtom psel;
  int hndl, nsel;
  hndl = ptr->NewSelection();
  ptr->Select( hndl, STYPE_ATOM, (char *)s.c_str(), SKEY_NEW );
  ptr->GetSelIndex( hndl, psel, nsel );
  DBAtom_selection result( psel, nsel );
  ptr->DeleteSelection( hndl );
  return result;
}

/*! Select residues using an MMDB coordinate ID expression.
  \param s The coordinate ID expression.
  \return The residue selection object. */
DBResidue_selection DBManager::select_residues( const String& s )
{
  mmdb::PPCResidue psel;
  int hndl, nsel;
  hndl = ptr->NewSelection();
  ptr->Select( hndl, STYPE_RESIDUE, (char *)s.c_str(), SKEY_NEW );
  ptr->GetSelIndex( hndl, psel, nsel );
  DBResidue_selection result( psel, nsel );
  ptr->DeleteSelection( hndl );
  return result;
}

/*! Select chains using an MMDB coordinate ID expression.
  \param s The coordinate ID expression.
  \return The chain selection object. */
DBChain_selection DBManager::select_chains( const String& s )
{
  mmdb::PPCChain psel;
  int hndl, nsel;
  hndl = ptr->NewSelection();
  ptr->Select( hndl, STYPE_CHAIN, (char *)s.c_str(), SKEY_NEW );
  ptr->GetSelIndex( hndl, psel, nsel );
  DBChain_selection result( psel, nsel );
  ptr->DeleteSelection( hndl );
  return result;
}

/*! Select models using an MMDB coordinate ID expression.
  \param s The coordinate ID expression.
  \return The model selection object. */
DBModel_selection DBManager::select_models( const String& s )
{
  mmdb::PPCModel psel;
  int hndl, nsel;
  hndl = ptr->NewSelection();
  ptr->Select( hndl, STYPE_MODEL, (char *)s.c_str(), SKEY_NEW );
  ptr->GetSelIndex( hndl, psel, nsel );
  DBModel_selection result( psel, nsel );
  ptr->DeleteSelection( hndl );
  return result;
}

/*! Select atoms which neighbour the atoms in some list.
  \param s The atom selection containing the list of atoms.
  \param r1 Lower bound on the neighbour distance.
  \param r2 Upper bound on the neighbour distance.
  \return The atom selection object. */
DBAtom_selection DBManager::select_atoms_near( DBAtom_selection& s, const ftype& r1, const ftype& r2 )
{
  mmdb::PPCAtom psel;
  int hndl, nsel;
  hndl = ptr->NewSelection();
  ptr->SelectNeighbours( hndl, STYPE_ATOM, s.ppcatom(), s.size(), r1, r2,
			 SKEY_NEW );
  ptr->GetSelIndex( hndl, psel, nsel );
  DBAtom_selection result( psel, nsel );
  ptr->DeleteSelection( hndl );
  return result;
}

/*! Select residues which neighbour the atoms in some list.
  \param s The atom selection containing the list of atoms.
  \param r1 Lower bound on the neighbour distance.
  \param r2 Upper bound on the neighbour distance.
  \return The atom selection object. */
DBResidue_selection DBManager::select_residues_near( DBAtom_selection& s, const ftype& r1, const ftype& r2 )
{
  mmdb::PPCResidue psel;
  int hndl, nsel;
  hndl = ptr->NewSelection();
  ptr->SelectNeighbours( hndl, STYPE_RESIDUE, s.ppcatom(), s.size(), r1, r2,
			 SKEY_NEW );
  ptr->GetSelIndex( hndl, psel, nsel );
  DBResidue_selection result( psel, nsel );
  ptr->DeleteSelection( hndl );
  return result;
}

/*! Select atoms by a sphere about defined by centre and radius.
  \param c The centre of the sphere.
  \param r The radius of the sphere.
  \return The atom selection object. */
DBAtom_selection DBManager::select_atoms_sphere( const Coord_orth& c, const ftype& r )
{
  mmdb::PPCAtom psel;
  int hndl, nsel;
  hndl = ptr->NewSelection();
  ptr->SelectSphere( hndl, STYPE_ATOM, c.x(), c.y(), c.z(), r, SKEY_NEW );
  ptr->GetSelIndex( hndl, psel, nsel );
  DBAtom_selection result( psel, nsel );
  ptr->DeleteSelection( hndl );
  return result;
}

/*! Select residues by a sphere defined by centre and radius.
  \param c The centre of the sphere.
  \param r The radius of the sphere.
  \return The residue selection object. */
DBResidue_selection DBManager::select_residues_sphere( const Coord_orth& c, const ftype& r )
{
  mmdb::PPCResidue psel;
  int hndl, nsel;
  hndl = ptr->NewSelection();
  ptr->SelectSphere( hndl, STYPE_RESIDUE, c.x(), c.y(), c.z(), r, SKEY_NEW );
  ptr->GetSelIndex( hndl, psel, nsel );
  DBResidue_selection result( psel, nsel );
  ptr->DeleteSelection( hndl );
  return result;
}

/*! Select atoms by a cylinder defined by centres of ends and radius.
  \param c1 The centre of one end.
  \param c2 The centre of the other end.
  \param r The radius of the cylinder.
  \return The atom selection object. */
DBAtom_selection DBManager::select_atoms_cylinder( const Coord_orth& c1, const Coord_orth& c2, const ftype& r )
{
  mmdb::PPCAtom psel;
  int hndl, nsel;
  hndl = ptr->NewSelection();
  ptr->SelectCylinder( hndl, STYPE_ATOM, c1.x(), c1.y(), c1.z(), c2.x(), c2.y(), c2.z(), r, SKEY_NEW );
  ptr->GetSelIndex( hndl, psel, nsel );
  DBAtom_selection result( psel, nsel );
  ptr->DeleteSelection( hndl );
  return result;
}

/*! Atoms are selected whose serial numbers lie between i1 and i2, inclusive.
  If both parameters are omitted. then all atoms are selected.
  \param i1 Lower bound on serial number.
  \param i2 Upper bound on serial number.
  \return The atom selection object. */
DBAtom_selection DBManager::select_atoms_serial( const int i1, const int i2 )
{
  mmdb::PPCAtom psel;
  int hndl, nsel;
  hndl = ptr->NewSelection();
  ptr->SelectAtoms( hndl, i1, i2, SKEY_NEW );
  ptr->GetSelIndex( hndl, psel, nsel );
  DBAtom_selection result( psel, nsel );
  ptr->DeleteSelection( hndl );
  return result;
}


// MMDB methods

/*! For later initialisation: see init() */
MMDB::MMDB() 
{
  InitMatType();
  ptr = new mmdb::XMMDBManager();
  Message::message( message_ctor_mmdb );
}

/*! \param spacegroup \param cell */
MMDB::MMDB( const Spacegroup& spacegroup, const Cell& cell )
{
  InitMatType();
  ptr = new mmdb::XMMDBManager();
  init( spacegroup, cell );
  Message::message( message_ctor_mmdb );
}

/*! WARNING: Unlike copying a DBManager which just creates a new smart pointer, this duplicates the entire MMDB. */
MMDB::MMDB( const MMDB& c )
{
  delete ptr;
  ptr = new mmdb::XMMDBManager();
  init( c.spacegroup(), c.cell() );
  ptr->Copy( c.pcmmdbmanager(), MMDBFCM_All );
  Message::message( message_ctor_mmdb );
}

MMDB::~MMDB()
{
  delete ptr;
  Message::message( message_dtor_mmdb );
}

/*! This can be used to build a model from scratch, or it can also be
  called after a model has been imported to change the spacegroup or
  cell info.
  \note This can be used to set a spacegroup which MMDB itself may not
  recognize. The symops and methods which depend on them should still
  work though.
  \param spacegroup The new spacegroup.
  \param cell The new cell.
 */
void MMDB::init( const Spacegroup& spacegroup, const Cell& cell )
{
  spacegroup_ = spacegroup;
  cell_ = cell;
  ptr->set_cell( cell_ );
  ptr->set_spacegroup( spacegroup_ );
}

bool MMDB::is_null() const
{ return ( spacegroup_.is_null() || cell_.is_null() ); }

const Spacegroup& MMDB::spacegroup() const
{ return spacegroup_; }

const Cell& MMDB::cell() const
{ return cell_; }

/*! The file may be either a PDB or mmCIF file.
  If the spacegroup or cell are not set, they will be taken from the
  file, otherwise the existing values override the ones in the file.
  \param file The filename (or pathname) of the file to read. */
void MMDB::read_file( const String& file )
{
  int err = ptr->ReadCoorFile( (char *)file.c_str() );
  if (err) Message::message( Message_fatal( "MMDB: read_file error: "+file+" : "+String(err) ) );
  // set spacegroup if necessary
  if ( spacegroup_.is_null() ) spacegroup_ = ptr->spacegroup();
  else                         ptr->set_spacegroup( spacegroup_ );
  // set cell if necessary
  if ( cell_.is_null() ) cell_ = ptr->cell();
  else                   ptr->set_cell( cell_ );
}

/*! The output file type will be the same as the file read, otherwise PDB.
  \param file The filename (or pathname) of the file to write.
  \param type 0=PDB, 1=CIF, 2=binary, default=same as input file, or PDB. */
void MMDB::write_file( const String& file, TYPE type )
{
  const TYPE types[3] = { PDB, CIF, Binary };
  int rtype = ptr->GetFileType();
  if ( type == Default && rtype >= 0 && rtype <= 2 ) type = types[ rtype ];
  int err;
  switch ( type ) {
  case Binary:
    err = ptr->WriteMMDBF( (char *)file.c_str() ); break;
  case CIF:
    err = ptr->WriteCIFASCII( (char *)file.c_str() ); break;
  case PDB:
  default:
    err = ptr->WritePDBASCII( (char *)file.c_str() ); break;
  }
  if (err) Message::message( Message_fatal( "MMDB: write_file error: "+file+" : "+String(err) ) ); 
}


void MMDB::debug() const
{
  std::cout << "MMDB object: clipper+mmdb cell, clipper+mmdb spgr\n";
  cell().debug();
  ptr->cell().debug();
  spacegroup().debug();
  ptr->spacegroup().debug();
}


} // namespace clipper
