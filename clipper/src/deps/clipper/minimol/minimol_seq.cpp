/* minimol_seq.cpp: atomic model types */
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

#include "minimol_seq.h"


namespace clipper {


// PolymerSequence

void MPolymerSequence::set_id( const String& s ) { id_ = id_tidy( s ); }

void MPolymerSequence::set_sequence( const String& s ) { seq_ = s; }

// MoleculeSequence

/*! Lookup monomer by ID. If mode=UNIQUE, the insertion code must match,
  otherwise the first monomer with the same sequence number is returned.
  \param n The monomer ID.  \param mode The search mode.
  \return The monomer. */
const MPolymerSequence& MMoleculeSequence::find( const String& n, const MM::MODE mode ) const
{
  int i = lookup( n, mode );
  if ( i < 0 ) Message::message(Message_fatal("MMolecule: no such monomer"));
  return children[i];
}

/*! See MMolecule::find() */
MPolymerSequence& MMoleculeSequence::find( const String& n, const MM::MODE mode )
{
  int i = lookup( n, mode );
  if ( i < 0 ) Message::message(Message_fatal("MMolecule: no such monomer"));
  return children[i];
}

int MMoleculeSequence::lookup( const String& str, const MM::MODE& mode ) const
{
  String sid = CHILDTYPE::id_tidy( str );
  for ( int i = 0; i < children.size(); i++ )
    if ( CHILDTYPE::id_match( sid, children[i].id(), mode ) ) return i;
  return -1;
}

void MMoleculeSequence::insert( const MPolymerSequence& add, int pos )
{
  if ( pos < 0 ) children.push_back( add );
  else children.insert( children.begin() + pos, add );
}


std::pair<std::vector<int>,std::vector<int> > MSequenceAlign::operator() ( const String& seq1, const String& seq2 ) const
{
  enum dirn { NUL, U, L, UL };  // directions: null, up, left, diag

  // pad sequences at start to allow first symbol to be aligned
  std::string s1 = " " + seq1;
  std::string s2 = " " + seq2;
  int n1 = s1.length();
  int n2 = s2.length();

  // initilize matrices.
  Matrix<ftype32> scores( n1, n2, 0.0 );
  Matrix<char>    dirns ( n1, n2, NUL );

  // now fill first row/col
  if ( type_ != LOCAL ) {
    for ( int i1 = 1; i1 < n1; i1++ ) {
      scores(i1,0) = scores(i1-1,0) + scrgap;
      dirns(i1,0)  = U;
    }
    for ( int i2 = 1; i2 < n2; i2++ ) {
      scores(0,i2) = scores(0,i2-1) + scrgap;
      dirns(0,i2)  = L;
    }
  }

  // fill the rest of the matrix
  for ( int i1 = 1; i1 < n1; i1++ )
    for ( int i2 = 1; i2 < n2; i2++ ) {
      // calc bonus for a match at this position
      ftype32 s = ( s1[i1] == s2[i2] ) ? scrmat : scrmis;
      // calc best score obtainable for this position
      ftype32 sul = scores(i1-1,i2-1) + s;
      ftype32 su  = scores(i1-1,i2  ) + scrgap;
      ftype32 sl  = scores(i1  ,i2-1) + scrgap;
      // and select
      if ( sul >= su && sul >= sl ) {
	scores(i1,i2) = sul;
	dirns(i1,i2)  = UL;
      } else if ( su > sl ) {
	scores(i1,i2) = su;
	dirns(i1,i2)  = U;
      } else {
	scores(i1,i2) = sl;
	dirns(i1,i2)  = L;
      }
      if ( type_ == LOCAL )
	if ( scores(i1,i2) <= 0.0 ) {
	  scores(i1,i2) = 0.0;
	  dirns(i1,i2)  = NUL;
	}
    }

  // now trace backwards to build up the best sequence alignment
  std::vector<int> r1(seq1.length(),-1), r2(seq2.length(),-1);
  int i1, i2;
  if ( type_ == LOCAL ) {  // local match: start with highest scoring position
    i1 = i2 = 0;
    for ( int j1 = 0; j1 < n1; j1++ )
      for ( int j2 = 0; j2 < n2; j2++ )
	if ( scores(j1,j2) > scores(i1,i2) ) {
	  i1 = j1;
	  i2 = j2;
	}
  } else {        // global match: start with last element
    i1 = n1 - 1;
    i2 = n2 - 1;
  }
  while ( dirns(i1,i2) != NUL ) {
    if        ( dirns(i1,i2) == U ) {
      i1 = i1 - 1;
    } else if ( dirns(i1,i2) == L ) {
      i2 = i2 - 1;
    } else {
      i1 = i1 - 1;
      i2 = i2 - 1;
      r1[i1] = i2;
      r2[i2] = i1;
    }
  }
  
  std::pair<std::vector<int>,std::vector<int> > result(r1,r2);
  return result;
}

} // namespace clipper
