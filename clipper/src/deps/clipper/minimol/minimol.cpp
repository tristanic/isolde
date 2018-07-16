/* minimol.cpp: atomic model types */
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


#include "minimol.h"
#include "minimol_utils.h"

extern "C" {
#include <string.h>
}


namespace clipper {


Message_ctor message_ctor_mmodel( " [MModel: constructed]" );


// helper functions

namespace MM {
 std::vector<clipper::String> path_split( const clipper::String& s, const clipper::String& sep )
 {
   std::vector<String> splitstr;
   size_t tokbeg = 0, tokend = 0;
   if ( s.find_first_of(sep, 0) == 0 ) tokbeg = 1;
   while (1) {
     tokend = s.find_first_of(sep, tokbeg);
     if (tokend == String::npos) break;
     splitstr.push_back( s.substr(tokbeg, tokend-tokbeg) );
     tokbeg = tokend+1;
   }
   splitstr.push_back( s.substr(tokbeg) );
   return splitstr;
 }
}


// Atom

MAtom::MAtom( const clipper::Atom& atom )
{
 set_element( atom.element() );
 set_coord_orth( atom.coord_orth() );
 set_occupancy( atom.occupancy() );
 set_u_iso( atom.u_iso() );
 set_u_aniso_orth( atom.u_aniso_orth() );
}

void MAtom::set_id( const String& s ) { id_ = id_tidy( s ); }

void MAtom::set_name( const String s, const String altconf )
{
 if ( altconf != "" ) set_id( s + ":" + altconf );
 else                 set_id( s );
}

String MAtom::id_tidy( const String& id )
{
 int pos = id.find( ":" );
 if ( pos == String::npos ) pos = id.length();
 String name( id.substr( 0, pos ) );
 String altc( id.substr( pos ) );
 if ( name.length() < 4 ) {
   name = name + "   ";
   if ( islower( name[1] ) )
     name[1] = toupper( name[1] );
   else
     name = " " + name;
 }
 return name.substr(0,4) + altc;
}

/*! copy from other atom. mode can be MM::COPY_M, COPY_P, COPY_MP,
 COPY_C, COPY_MC, COPY_PC, COPY_MPC, where M means copy members, P
 means copy PropertyMananger properties, and C means copy
 children. Children are copied with the same option. The values
 'MEMBERS', 'PROPERTIES', 'CHILDREN' can also be used. */
MAtom& MAtom::copy( const MAtom& other, const MM::COPY& mode )
{
 if ( mode & MM::COPY_M ) atom() = other.atom();
 if ( mode & MM::COPY_M ) id_ = other.id_;
 if ( mode & MM::COPY_P ) PropertyManager::copy( other );
 return *this;
}

bool MAtom::id_match( const String& id1, const String& id2, const MM::MODE& mode )
{ if ( mode == MM::UNIQUE ) return ( id1 == id2 );
 else return ( id1.substr(0,4) == id2.substr(0,4) ); }


// Monomer

void MMonomer::set_id( const String& s )
{
 id_ = id_tidy( s );
}

void MMonomer::set_type( const String& s )
{ type_ = s; }

void MMonomer::set_seqnum( const int s, const String inscode )
{
 if ( inscode != "" ) set_id( String( s, 4 ) + ":" + inscode );
 else                 set_id( String( s, 4 ) );
}

Atom_list MMonomer::atom_list() const
{
 Atom_list list;
 for ( int a = 0; a < children.size(); a++ )
   list.push_back( Atom( children[a] ) );
 return list;
}

void MMonomer::transform( const RTop_orth rt )
{ for ( int i = 0; i < children.size(); i++ ) children[i].transform( rt ); }

/*! Creates a copy of this monomer containing only the atoms described
 by the selection string. '*' copies all atoms.

 The atom selection must contain an atom ID or a comma separated list
 of atom IDs, or '*' to select all atom. Atom IDs are described in
 s_mm_atom_id.

 The selection string must contain an atom ID or a comma separated
 list of atom IDs. Atom IDs are described in s_mm_atom_id.

 \param sel The selection string.
 \param mode MM::UNIQUE forces an exact match, including alternate conformation code. MM::ANY matches every atom with the right name, ignoring alternate conformation codes.
 \return The selection as a new monomer. */
MMonomer MMonomer::select( const String& sel, const MM::MODE mode ) const
{
 std::vector<String> path = MM::path_split( sel, "/" );
 while ( path.size() < 1 ) path.push_back( "*" );
 MMonomer result;
 result.copy( *this, MM::COPY_MP );
 if ( path[0].trim() == "*" ) {
   for ( int i = 0; i < children.size(); i++ ) result.insert( children[i] );
 } else {
   std::vector<String> list = MM::path_split( path[0], "," );
   for ( int j = 0; j < list.size(); j++ ) {
     String sid = CHILDTYPE::id_tidy( list[j] );
     for ( int i = 0; i < children.size(); i++ )
	if ( CHILDTYPE::id_match( sid, children[i].id(), mode ) )
	  result.insert( children[i] );
   }
 }
 return result;
}

/*! Creates a list of inidices of children matching the given selection string.
 See select().
 \param sel The selection string.
 \param mode MM::UNIQUE or MM::ANY.
 \return A vector of integer indices of the matching children. */
std::vector<int> MMonomer::select_index( const String& sel, const MM::MODE mode ) const
{
 std::vector<int> result;
 if ( sel.trim() == "*" ) {
   for ( int i = 0; i < children.size(); i++ ) result.push_back( i );
 } else {
   std::vector<String> list = MM::path_split( sel, "," );
   for ( int j = 0; j < list.size(); j++ ) {
     String sid = CHILDTYPE::id_tidy( list[j] );
     for ( int i = 0; i < children.size(); i++ )
	if ( CHILDTYPE::id_match( sid, children[i].id(), mode ) )
	  result.push_back( i );
   }
 }
 return result;
}

/*! Lookup atom by ID. If mode=UNIQUE, the alternate conformation code
 must match, otherwise the first atom with the same name is returned.
 \param n The atom ID.  \param mode The search mode.
 \return The atom. */
const MAtom& MMonomer::find( const String& n, const MM::MODE mode ) const
{
 int i = lookup( n, mode );
 if ( i < 0 ) Message::message(Message_fatal("MMonomer: no such atom"));
 return children[i];
}

/*! See MMonomer::find() */
MAtom& MMonomer::find( const String& n, const MM::MODE mode )
{
 int i = lookup( n, mode );
 if ( i < 0 ) Message::message(Message_fatal("MMonomer: no such atom"));
 return children[i];
}

int MMonomer::lookup( const String& str, const MM::MODE& mode ) const
{
 String sid = CHILDTYPE::id_tidy( str );
 for ( int i = 0; i < children.size(); i++ )
   if ( CHILDTYPE::id_match( sid, children[i].id(), mode ) ) return i;
 return -1;
}

void MMonomer::insert( const MAtom& add, int pos )
{
 if ( pos < 0 ) children.push_back( add );
 else children.insert( children.begin() + pos, add );
}

/*! \return atoms which are in both monomers. */
MMonomer operator& ( const MMonomer& m1, const MMonomer& m2 )
{
 MMonomer result;
 result.copy( m1, MM::COPY_MP );
 for ( int i1 = 0; i1 < m1.size(); i1++ )
   for ( int i2 = 0; i2 < m2.size(); i2++ )
     if ( m1[i1].id() == m2[i2].id() ) {
	result.insert( m1[i1] );
	break;
     }
 return result;
}

/*! \return atoms which are in either monomer. */
MMonomer operator| ( const MMonomer& m1, const MMonomer& m2 )
{
 MMonomer result;
 result.copy( m1, MM::COPY_MP );
 int i, j;
 for ( i = 0; i < m1.size(); i++ ) {
   for ( j = 0; j < result.size(); j++ )
     if ( m1[i].id() == result[j].id() ) break;
   if ( j == result.size() )
     result.insert( m1[i] );
 }
 for ( i = 0; i < m2.size(); i++ ) {
   for ( j = 0; j < result.size(); j++ )
     if ( m2[i].id() == result[j].id() ) break;
   if ( j == result.size() )
     result.insert( m2[i] );
 }
 return result;
}

/*! copy from other atom. mode can be MM::COPY_M, COPY_P, COPY_MP,
 COPY_C, COPY_MC, COPY_PC, COPY_MPC, where M means copy members, P
 means copy PropertyMananger properties, and C means copy
 children. Children are copied with the same option. The values
 'MEMBERS', 'PROPERTIES', 'CHILDREN' can also be used. */
MMonomer& MMonomer::copy( const MMonomer& other, const MM::COPY& mode )
{
 if ( mode & MM::COPY_M ) id_ = other.id_;
 if ( mode & MM::COPY_M ) type_ = other.type_;
 if ( mode & MM::COPY_P ) PropertyManager::copy( other );
 if ( mode & MM::COPY_C ) {
   children.resize( other.size() );
   for ( int i = 0; i < size(); i++ ) children[i].copy( other[i], mode );
 }
 return *this;
}

String MMonomer::id_tidy( const String& id )
{
 int pos = id.find( ":" );
 if ( pos == String::npos )
   return String( id.i(), 4 );
 else
   return String( id.i(), 4 ) + id.substr( pos );
}

bool MMonomer::id_match( const String& id1, const String& id2, const MM::MODE& mode )
{ if ( mode == MM::UNIQUE ) return ( id1 == id2 );
 else return ( id1.substr(0,4) == id2.substr(0,4) ); }


// Polymer

void MPolymer::set_id( const String& s ) { id_ = id_tidy( s ); }

Atom_list MPolymer::atom_list() const
{
 Atom_list list;
 for ( int m = 0; m < children.size(); m++ )
   for ( int a = 0; a < children[m].size(); a++ )
     list.push_back( Atom( children[m][a] ) );
 return list;
}

void MPolymer::transform( const RTop_orth rt )
{ for ( int i = 0; i < children.size(); i++ ) children[i].transform( rt ); }

/*! Creates a copy of this polymer containing only the monomers and
 atoms described by the selection string.

 The selection string must be of the form 'X/Y' where X is a monomer
 selection and Y is an atom selection, described under
 MAtom::select(). The monomer selection must contain a monomer ID or
 a comma separated list of monomer IDs, or '*' to select all
 monomers. Monomer IDs are described in s_mm_monomer_id.

 \param sel The selection string.
 \param mode MM::UNIQUE forces an exact match, including insertion code.
             MM::ANY matches any monomer with the right sequence number,
	      ignoring insertion code.
 \return The selection as a new polymer. */
MPolymer MPolymer::select( const String& sel, const MM::MODE mode ) const
{
 std::vector<String> path = MM::path_split( sel, "/" );
 while ( path.size() < 2 ) path.push_back( "*" );
 MPolymer result;
 result.copy( *this, MM::COPY_MP );
 if ( path[0].trim() == "*" ) {
   for ( int i = 0; i < children.size(); i++ )
     result.insert( children[i].select( path[1], mode ) );
 } else {
   std::vector<String> list = MM::path_split( path[0], "," );
   for ( int j = 0; j < list.size(); j++ ) {
     String sid = CHILDTYPE::id_tidy( list[j] );
     for ( int i = 0; i < children.size(); i++ )
	if ( CHILDTYPE::id_match( sid, children[i].id(), mode ) )
	  result.insert( children[i].select( path[1], mode ) );
   }
 }
 return result;
}

/*! Creates a list of inidices of children matching the given selection string.
 See select().
 \param sel The selection string.
 \param mode MM::UNIQUE or MM::ANY.
 \return A vector of integer indices of the matching children. */
std::vector<int> MPolymer::select_index( const String& sel, const MM::MODE mode ) const
{
 std::vector<int> result;
 if ( sel.trim() == "*" ) {
   for ( int i = 0; i < children.size(); i++ )
     result.push_back( i );
 } else {
   std::vector<String> list = MM::path_split( sel, "," );
   for ( int j = 0; j < list.size(); j++ ) {
     String sid = CHILDTYPE::id_tidy( list[j] );
     for ( int i = 0; i < children.size(); i++ )
	if ( CHILDTYPE::id_match( sid, children[i].id(), mode ) )
	  result.push_back( i );
   }
 }
 return result;
}

/*! Lookup monomer by ID. If mode=UNIQUE, the insertion code must match,
 otherwise the first monomer with the same sequence number is returned.
 \param n The monomer ID.  \param mode The search mode.
 \return The monomer. */
const MMonomer& MPolymer::find( const String& n, const MM::MODE mode ) const
{
 int i = lookup( n, mode );
 if ( i < 0 ) Message::message(Message_fatal("MPolymer: no such monomer"));
 return children[i];
}

/*! See MPolymer::find() */
MMonomer& MPolymer::find( const String& n, const MM::MODE mode )
{
 int i = lookup( n, mode );
 if ( i < 0 ) Message::message(Message_fatal("MPolymer: no such monomer"));
 return children[i];
}

int MPolymer::lookup( const String& str, const MM::MODE& mode ) const
{
 String sid = CHILDTYPE::id_tidy( str );
 for ( int i = 0; i < children.size(); i++ )
   if ( CHILDTYPE::id_match( sid, children[i].id(), mode ) ) return i;
 return -1;
}

void MPolymer::insert( const MMonomer& add, int pos )
{
 if ( pos < 0 ) children.push_back( add );
 else children.insert( children.begin() + pos, add );
}

/*! \return monomers and atoms which are in both polymers. */
MPolymer operator& ( const MPolymer& m1, const MPolymer& m2 )
{
 MPolymer result;
 result.copy( m1, MM::COPY_MP );
 for ( int i1 = 0; i1 < m1.size(); i1++ )
   for ( int i2 = 0; i2 < m2.size(); i2++ )
     if ( m1[i1].id() == m2[i2].id() ) {
	result.insert( m1[i1] & m2[i2] );
	break;
     }
 return result;
}

/*! \return monomers and atoms which are in either polymer. */
MPolymer operator| ( const MPolymer& m1, const MPolymer& m2 )
{
 MPolymer result;
 result.copy( m1, MM::COPY_MP );
 int i, j;
 for ( i = 0; i < m1.size(); i++ ) {
   for ( j = 0; j < result.size(); j++ )
     if ( m1[i].id() == result[j].id() ) break;
   if ( j == result.size() )
     result.insert( m1[i] );
   else
     result[j] = result[j] | m1[i];
 }
 for ( i = 0; i < m2.size(); i++ ) {
   for ( j = 0; j < result.size(); j++ )
     if ( m2[i].id() == result[j].id() ) break;
   if ( j == result.size() )
     result.insert( m2[i] );
   else
     result[j] = result[j] | m2[i];
 }
 return result;
}

/*! copy from other atom. mode can be MM::COPY_M, COPY_P, COPY_MP,
 COPY_C, COPY_MC, COPY_PC, COPY_MPC, where M means copy members, P
 means copy PropertyMananger properties, and C means copy
 children. Children are copied with the same option. The values
 'MEMBERS', 'PROPERTIES', 'CHILDREN' can also be used. */
MPolymer& MPolymer::copy( const MPolymer& other, const MM::COPY& mode )
{
 if ( mode & MM::COPY_M ) id_ = other.id_;
 if ( mode & MM::COPY_P ) PropertyManager::copy( other );
 if ( mode & MM::COPY_C ) {
   children.resize( other.size() );
   for ( int i = 0; i < size(); i++ ) children[i].copy( other[i], mode );
 }
 return *this;
}

String MPolymer::id_tidy( const String& id ) { return id; }
bool MPolymer::id_match( const String& id1, const String& id2, const MM::MODE& mode ) { return ( id1 == id2 ); }


// Model

Atom_list MModel::atom_list() const
{
 Atom_list list;
 for ( int p = 0; p < children.size(); p++ )
   for ( int m = 0; m < children[p].size(); m++ )
     for ( int a = 0; a < children[p][m].size(); a++ )
       list.push_back( Atom( children[p][m][a] ) );
 return list;
}

void MModel::transform( const RTop_orth rt )
{ for ( int i = 0; i < children.size(); i++ ) children[i].transform( rt ); }

/*! Creates a copy of this model containing only the polymers,
 monomers and atoms described by the selection string.

 The selection string must be of the form 'X/Y/Z' where X is a
 polymer selection, Y is a monomer selection described under
 MMonomer::select(), and Z is an atom selection described under
 MAtom::select(). The polymer selection must contain a polymer ID or
 a comma separated list of polymer IDs, or '*' to select all
 polymers. Polymer IDs are described in s_mm_monomer_id.

 See s_mm_selections for examples.

 \param sel The selection string.
 \param mode No effect.
 \return The selection as a new model. */
MModel MModel::select( const String& sel, const MM::MODE mode ) const
{
 std::vector<String> path = MM::path_split( sel, "/" );
 while ( path.size() < 3 ) path.push_back( "*" );
 MModel result;
 result.copy( *this, MM::COPY_MP );
 if ( path[0].trim() == "*" ) {
   for ( int i = 0; i < children.size(); i++ )
     result.insert( children[i].select( path[1]+"/"+path[2], mode ) );
 } else {
   std::vector<String> list = MM::path_split( path[0], "," );
   for ( int j = 0; j < list.size(); j++ ) {
     String sid = CHILDTYPE::id_tidy( list[j] );
     for ( int i = 0; i < children.size(); i++ )
	if ( CHILDTYPE::id_match( sid, children[i].id(), mode ) )
	  result.insert( children[i].select( path[1]+"/"+path[2], mode ) );
   }
 }
 return result;
}

/*! Creates a list of inidices of children matching the given selection string.
 See select().
 \param sel The selection string.
 \param mode MM::UNIQUE or MM::ANY.
 \return A vector of integer indices of the matching children. */
std::vector<int> MModel::select_index( const String& sel, const MM::MODE mode ) const
{
 std::vector<int> result;
 if ( sel.trim() == "*" ) {
   for ( int i = 0; i < children.size(); i++ )
     result.push_back( i );
 } else {
   std::vector<String> list = MM::path_split( sel, "," );
   for ( int j = 0; j < list.size(); j++ ) {
     String sid = CHILDTYPE::id_tidy( list[j] );
     for ( int i = 0; i < children.size(); i++ )
	if ( CHILDTYPE::id_match( sid, children[i].id(), mode ) )
	  result.push_back( i );
   }
 }
 return result;
}

/*! Lookup polymer by ID. Currently, mode is ignored.
 \param n The monomer ID.  \param mode The search mode.
 \return The polymer. */
const MPolymer& MModel::find( const String& n, const MM::MODE mode ) const
{
 int i = lookup( n, mode );
 if ( i < 0 ) Message::message(Message_fatal("MModel: no such polymer"));
 return children[i];
}

/*! See MModel::find() */
MPolymer& MModel::find( const String& n, const MM::MODE mode )
{
 int i = lookup( n, mode );
 if ( i < 0 ) Message::message(Message_fatal("MModel: no such polymer"));
 return children[i];
}

int MModel::lookup( const String& str, const MM::MODE& mode ) const
{
 String sid = CHILDTYPE::id_tidy( str );
 for ( int i = 0; i < children.size(); i++ )
   if ( CHILDTYPE::id_match( sid, children[i].id(), mode ) ) return i;
 return -1;
}

void MModel::insert( const MPolymer& add, int pos )
{
 if ( pos < 0 ) children.push_back( add );
 else children.insert( children.begin() + pos, add );
}

/*! \return polymers, monomers and atoms which are in both models. */
MModel operator& ( const MModel& m1, const MModel& m2 )
{
 MModel result;
 result.copy( m1, MM::COPY_MP );
 for ( int i1 = 0; i1 < m1.size(); i1++ )
   for ( int i2 = 0; i2 < m2.size(); i2++ )
     if ( m1[i1].id() == m2[i2].id() ) {
	result.insert( m1[i1] & m2[i2] );
	break;
     }
 return result;
}

/*! \return polymers, monomers and atoms which are in either model. */
MModel operator| ( const MModel& m1, const MModel& m2 )
{
 MModel result;
 result.copy( m1, MM::COPY_MP );
 int i, j;
 for ( i = 0; i < m1.size(); i++ ) {
   for ( j = 0; j < result.size(); j++ )
     if ( m1[i].id() == result[j].id() ) break;
   if ( j == result.size() )
     result.insert( m1[i] );
   else
     result[j] = result[j] | m1[i];
 }
 for ( i = 0; i < m2.size(); i++ ) {
   for ( j = 0; j < result.size(); j++ )
     if ( m2[i].id() == result[j].id() ) break;
   if ( j == result.size() )
     result.insert( m2[i] );
   else
     result[j] = result[j] | m2[i];
 }
 return result;
}

/*! copy from other atom. mode can be MM::COPY_M, COPY_P, COPY_MP,
 COPY_C, COPY_MC, COPY_PC, COPY_MPC, where M means copy members, P
 means copy PropertyMananger properties, and C means copy
 children. Children are copied with the same option. The values
 'MEMBERS', 'PROPERTIES', 'CHILDREN' can also be used. */
MModel& MModel::copy( const MModel& other, const MM::COPY& mode )
{
 if ( mode & MM::COPY_P ) PropertyManager::copy( other );
 if ( mode & MM::COPY_C ) {
   children.resize( other.size() );
   for ( int i = 0; i < size(); i++ ) children[i].copy( other[i], mode );
 }
 return *this;
}

/*! Return an atom on the basis of the MAtomIndex. Behaviour is
   undefined if the index is null.
 \param index The index of the atom in the heirarchy.
 \return A reference to the MAtom. */
const MAtom& MModel::atom( const MAtomIndex& index ) const
{ return children[index.polymer()][index.monomer()][index.atom()]; }

/*! Return an atom on the basis of the MAtomIndex. Behaviour is
   undefined if the index is null.
 \param index The index of the atom in the heirarchy.
 \return A reference to the MAtom. */
MAtom& MModel::atom( const MAtomIndex& index )
{ return children[index.polymer()][index.monomer()][index.atom()]; }

std::vector<MAtomIndex> MModel::select_atom_index( const String& sel, const MM::MODE mode ) const
{
 std::vector<MAtomIndex> result;
 std::vector<String> path = MM::path_split( sel, "/" );
 while ( path.size() < 3 ) path.push_back( "*" );
 const std::vector<int> mps =
   select_index( path[0], mode );
 for ( int p = 0; p < mps.size(); p++ ) {
   const std::vector<int> mms =
     children[mps[p]].select_index( path[1], mode );
   for ( int m = 0; m < mms.size(); m++ ) {
     const std::vector<int> mas =
	children[mps[p]][mms[m]].select_index( path[2], mode );
     for ( int a = 0; a < mas.size(); a++ ) {
	result.push_back( MAtomIndex( mps[p], mms[m], mas[a] ) );
     }
   }
 }
 return result;
}


// MiniMol

MiniMol::MiniMol()
{ Message::message( message_ctor_mmodel ); }

/*! The object is constructed with no atoms.
 \param spacegroup the spacegroup.
 \param cell the cell. */
MiniMol::MiniMol( const Spacegroup& spacegroup, const Cell& cell )
{
 init( spacegroup, cell );
 Message::message( message_ctor_mmodel );
}

/*! The object is initialised with no atoms.
 \param spacegroup the spacegroup.
 \param cell the cell. */
void MiniMol::init( const Spacegroup& spacegroup, const Cell& cell )
{
 spacegroup_ = spacegroup;
 cell_ = cell;
}

MAtom MiniMol::symmetry_atom( const MAtomIndexSymmetry& index )
{
 MAtom atom = MModel::atom(index);
 atom.set_coord_orth((spacegroup_.symop(index.symmetry())*atom.coord_orth().coord_frac(cell_)).coord_orth(cell_));
 return atom;
}

bool MiniMol::is_null() const
{ return ( spacegroup_.is_null() || cell_.is_null() ); }




// UTILITY FUNCTIONS:
// e.g. protein specific tools.


MMonomer::TYPE MMonomer::default_type_ = MMonomer::Richardson;

/*! A carbonyl oxygen is added to this residue if the supplied residue
 contains an appriate N atom bonded to the C. Otherwise, nothing
 happens.
 \param next The next monomer in the chain.
*/
void MMonomer::protein_mainchain_build_carbonyl_oxygen( const MMonomer& next ) {
 // check for mainchain atoms
 int a1 = lookup( " CA ", MM::ANY );
 int c1 = lookup( " C  ", MM::ANY );
 int n2 = next.lookup( " N  ", MM::ANY );
 if ( a1 < 0 || c1 < 0 || n2 < 0 ) return;
 // get coordinates and check bonding
 const clipper::Coord_orth cc1 = children[c1].coord_orth();
 const clipper::Coord_orth ca1 = children[a1].coord_orth() - cc1;
 const clipper::Coord_orth cn2 = next[n2].coord_orth() - cc1;
 if ( cn2.lengthsq() > 2.2 ) return;
 double uiso = children[c1].u_iso();
 double occ  = children[c1].occupancy();
 // delete any existing O
 int o1 = lookup( " O  ", MM::ANY );
 if ( o1 >= 0 ) children.erase( children.begin()+o1 );
 // add the atom
 const clipper::Vec3<> v0 = ca1.unit();
 const clipper::Vec3<> v1 = ( clipper::Vec3<>::cross( v0, clipper::Vec3<>::cross( v0, cn2 ) ) ).unit();
 MAtom atm = MAtom::null();
 atm.set_id( " O  " );
 atm.set_element( "O" );
 // length 1.24 angle 2.11
 atm.set_coord_orth( cc1 + clipper::Coord_orth( -0.637*v0 + 1.064*v1 ) );
 atm.set_occupancy( occ );
 atm.set_u_iso( uiso );
 insert( atm );
}

/*! A carbonyl oxygen is added to this residue in a default
 conformation, if it contains N, CA, and C atoms. Otherwise, nothing
 happens.
*/
void MMonomer::protein_mainchain_build_carbonyl_oxygen() {
 // check for mainchain atoms
 int n1 = lookup( " N  ", MM::ANY );
 int a1 = lookup( " CA ", MM::ANY );
 int c1 = lookup( " C  ", MM::ANY );
 if ( n1 < 0 || a1 < 0 || c1 < 0 ) return;
 // get coordinates and check bonding
 const clipper::Coord_orth cn1 = children[n1].coord_orth();
 const clipper::Coord_orth ca1 = children[a1].coord_orth();
 const clipper::Coord_orth cc1 = children[c1].coord_orth();
 double uiso = children[c1].u_iso();
 double occ  = children[c1].occupancy();
 // delete any existing O
 int o1 = lookup( " O  ", MM::ANY );
 if ( o1 >= 0 ) children.erase( children.begin()+o1 );
 MAtom atm = MAtom::null();
 atm.set_id( " O  " );
 atm.set_element( "O" );
 // length 1.24 angle 2.09
 atm.set_coord_orth( clipper::Coord_orth( cn1, ca1, cc1, 1.24, 2.09, -0.58 ) );
 atm.set_occupancy( occ );
 atm.set_u_iso( uiso );
 insert( atm );
}

// internal function for fast lookup in rotamer lib
int MMonomer::rotamer_find( String res, int rota, TYPE t ) {
 int rotamer_data_size;
 if ( t == Dunbrack ) rotamer_data_size = data::rotamer_data_dunbrack_size;
 else                 rotamer_data_size = data::rotamer_data_richardson_size;
 const data::Rotamer_data* rotamer_data;
 if ( t == Dunbrack ) rotamer_data = data::rotamer_data_dunbrack;
 else                 rotamer_data = data::rotamer_data_richardson;
 if ( res.length() < 3 ) return 0;
 int p1 = -1;
 int p2 = rotamer_data_size - 1;
 while ( p2 - p1 > 1 ) {
   int p = ( p1 + p2 ) / 2;
   int s = strncmp( res.c_str(), rotamer_data[p].resname, 3 );
   if ( s < 0 || ( s == 0 && rota <= rotamer_data[p].rota ) )
     p2 = p;
   else
     p1 = p;
 }
 if ( strncmp( res.c_str(), rotamer_data[p2].resname, 3 ) == 0 &&
      rota == rotamer_data[p2].rota ) return p2;
 return -1;
}

/*! \return The number of stored rotamers for this residue type.
 0 if unknown. */
int MMonomer::protein_sidechain_number_of_rotamers( TYPE t ) const {
 const data::Rotamer_data* rotamer_data;
 if ( t == Dunbrack ) rotamer_data = data::rotamer_data_dunbrack;
 else                 rotamer_data = data::rotamer_data_richardson;
 int r = rotamer_find( type(), 0, t );
 if ( r < 0 ) return 0;
 return rotamer_data[r].num_rota;
}

/*!
 \param n The number of the rotamer required.
 \return The frequency of the given rotamer. */
ftype MMonomer::protein_sidechain_build_rotamer( const int& n, TYPE t ) {
 const data::Rotamer_data* rotamer_data;
 if ( t == Dunbrack ) rotamer_data = data::rotamer_data_dunbrack;
 else                 rotamer_data = data::rotamer_data_richardson;
 int na = lookup( " CA ", MM::ANY );
 int nc = lookup( " C  ", MM::ANY );
 int nn = lookup( " N  ", MM::ANY );
 if ( na < 0 || nc < 0 || nn < 0 ) return 0.0;
 clipper::Coord_orth ca = children[na].coord_orth();
 clipper::Coord_orth c1 = children[nc].coord_orth() - ca;
 clipper::Coord_orth c2 = children[nn].coord_orth() - ca;
 double uiso = children[na].u_iso();
 double occ  = children[na].occupancy();
 // strip old side chain
 for ( int i = children.size()-1; i >= 0; i-- )
   if ( children[i].name() != " CA " && children[i].name() != " N  " &&
	 children[i].name() != " C  " && children[i].name() != " O  " )
     children.erase( children.begin()+i );
 // get rotamer
 int r = rotamer_find( type(), n, t );
 if ( r < 0 ) return 0.0;
 if ( n >= rotamer_data[r].num_rota ) return -1.0;
 // get rtop from standard orientation
 const clipper::Vec3<> v1( (c1.unit()+c2.unit()).unit() );
 const clipper::Vec3<> v2( clipper::Vec3<>::cross(c1,c2).unit() );
 const clipper::Vec3<> v3( clipper::Vec3<>::cross(v1,v2).unit() );
 const clipper::Mat33<> mr( v1[0], v2[0], v3[0],
			     v1[1], v2[1], v3[1],
			     v1[2], v2[2], v3[2] );
 clipper::RTop_orth rtop( mr, ca );
 // add new atoms
 MAtom atm = MAtom::null();
 for ( int dr = 0; dr < rotamer_data[r].num_atom; dr++ ) {
   int i = r + dr;
   String name = rotamer_data[i].atomname;
   atm.set_id( name );
   name = name.substr( 0, 2 );
   name = name.trim();
   atm.set_element( name );
   atm.set_coord_orth( rtop * Coord_orth( rotamer_data[i].x, 
					   rotamer_data[i].y, 
					   rotamer_data[i].z ) );
   atm.set_occupancy( occ );
   atm.set_u_iso( uiso );
   insert( atm );
 }

 return rotamer_data[r].rota_prob;
}


/*! Test if the C of residue 1 is bonded to the N of residue 2,
 within the distance r.
 \param r1 The first residue.
 \param r2 The second residue.
 \param r The maximum allowed bond length.
 \return true if N and C are present and bonded. */
bool MMonomer::protein_peptide_bond( const MMonomer& m1, const MMonomer& m2, ftype r ) {
 int c1 = m1.lookup( " C  ", MM::ANY );
 int n2 = m2.lookup( " N  ", MM::ANY );
 if ( c1 >= 0 && n2 >= 0 )
   if ( ( m1[c1].coord_orth() - m2[n2].coord_orth() ).lengthsq() < r*r )
     return true;
 return false;
}

/*! Return the Ramachadran angle in radians on -pi...pi.
 To check the result, see clipper::Util::is_nan()
 \param r1 The first residue.
 \param r2 The second residue.
 \return The torsion angle in radians, or NaN if atoms are missing. */
double MMonomer::protein_ramachandran_phi( const MMonomer& m1, const MMonomer& m2 )
{
 ftype result = clipper::Util::nan();
 int index_cx = m1.lookup( " C  ", clipper::MM::ANY );
 int index_n  = m2.lookup( " N  ", clipper::MM::ANY );
 int index_ca = m2.lookup( " CA ", clipper::MM::ANY );
 int index_c  = m2.lookup( " C  ", clipper::MM::ANY );
 // if we have all three atoms, then add residue
 if ( index_cx >= 0 && index_ca >= 0 && index_c >= 0 && index_n >= 0 ) {
   Coord_orth coord_cx = m1[index_cx].coord_orth();
   Coord_orth coord_n  = m2[index_n ].coord_orth();
   Coord_orth coord_ca = m2[index_ca].coord_orth();
   Coord_orth coord_c  = m2[index_c ].coord_orth();
   // ramachandran calc
   result = Coord_orth::torsion( coord_cx, coord_n, coord_ca, coord_c );
 }
 return result;
}

/*! Return the Ramachadran angle in radians on -pi...pi.
 To check the result, see clipper::Util::is_nan()
 \param r1 The first residue.
 \param r2 The second residue.
 \return The torsion angle in radians, or NaN if atoms are missing. */
double MMonomer::protein_ramachandran_psi( const MMonomer& m1, const MMonomer& m2 )
{
 ftype result = clipper::Util::nan();
 int index_n  = m1.lookup( " N  ", clipper::MM::ANY );
 int index_ca = m1.lookup( " CA ", clipper::MM::ANY );
 int index_c  = m1.lookup( " C  ", clipper::MM::ANY );
 int index_nx = m2.lookup( " N  ", clipper::MM::ANY );
 // if we have all three atoms, then add residue
 if ( index_ca >= 0 && index_c >= 0 && index_n >= 0 && index_nx >= 0 ) {
   Coord_orth coord_n  = m1[index_n ].coord_orth();
   Coord_orth coord_ca = m1[index_ca].coord_orth();
   Coord_orth coord_c  = m1[index_c ].coord_orth();
   Coord_orth coord_nx = m2[index_nx].coord_orth();
   // ramachandran calc
   result = Coord_orth::torsion( coord_n, coord_ca, coord_c, coord_nx );
 }
 return result;
}


} // namespace clipper
