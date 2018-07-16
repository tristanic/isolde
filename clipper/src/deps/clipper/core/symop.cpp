/* symop.cpp: fundamental data types for the clipper libraries */
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


#include "symop.h"

#include "coords.h"


namespace clipper {


/*! Construct an RT operator from a string description,
  e.g. 1/2x,z-y+2/3,x
  '*' is optional for multiplication, commas are compulsory. */
RTop_frac::RTop_frac( const String& strop )
{
  std::vector<String> rows = strop.split(",");
  if ( rows.size() != 3 )
    Message::message(Message_fatal("RTop_frac: invalid string:"+strop));

  rot() = Mat33<>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  trn() = Vec3<> (0.0,0.0,0.0);

  std::vector<String> cols;
  String num;
  ftype val;
  int nrow, ncol, nchr, npart;
  for ( nrow = 0; nrow < 3; nrow++ ) {
    const String& row = rows[nrow];
    cols.clear(); cols.push_back("");
    for ( nchr = 0; nchr < row.length(); nchr++ ) {
      const char& c = row[nchr];
      if ( c == '+' || c == '-' ) cols.push_back( "" );
      cols.back() += c;
    }
    for ( npart = 0; npart < cols.size(); npart++ ) {
      const String& col = cols[npart];
      ncol = 3;
      num = "";
      for ( nchr = 0; nchr < col.length(); nchr++ ) {
	const char& c = col[nchr];
	if      ( c == 'x' || c == 'X' ) ncol = 0;
	else if ( c == 'y' || c == 'Y' ) ncol = 1;
      	else if ( c == 'z' || c == 'Z' ) ncol = 2;
	else if ( c == '-' || c == '/' || isdigit(c) ) num += c;
      }
      if ( ncol < 3 ) {
	if ( num == "" || num == "+" || num == "-" ) num += "1";
	val = num.rational();
	rot()( nrow, ncol ) = val;
      } else {
	if ( num == "" || num == "+" || num == "-" ) num += "0";
	val = num.rational();
	trn()[ nrow ] = val;
      }
    }
  }
}

/*! \param cell The cell concerned \return The transformed coordinate. */
RTop_orth RTop_frac::rtop_orth( const Cell& cell ) const
{
  return RTop_orth( RTop<>(cell.matrix_orth()) * (*this) * RTop<>(cell.matrix_frac()) );
}

/*! \return The inverse of the operator. */
RTop_frac RTop_frac::inverse() const
{ return RTop_frac( RTop<>::inverse() ); }

/*! \return The identity operator. */
RTop_frac RTop_frac::identity()
{ return RTop_frac( RTop<>::identity() ); }

/*! \return The null (uninitialised) operator. */
RTop_frac RTop_frac::null()
{ return RTop_frac( RTop<>::null() ); }


/*! Construct a symmetry operator and initialise it to the supplied RTop.
  Translations are rounded to a basis of 48, and put on the range 0..1
  \param mat The RTop to use. */
Symop::Symop( const RTop<>& rt )
{
  // initialise to the supplied matrix
  int i, j;
  for ( i=0; i<3; i++) for ( j=0; j<3; j++)
    rot()(i,j) = rint( rt.rot()(i,j) );
  for ( i=0; i<3; i++)
    trn()[i] = ftype( Util::mod( Util::intr(48.*rt.trn()[i]), 48 ) ) / 48.;
}

/*! Construct a symmetry operator and initialise it to the supplied matrix.
  Translations are rounded to a basis of 48, and put on the range 0..1
  \param mat The 4x4 matrix to use. The [i][3] elements contain the
  translation. */
Symop::Symop( const ftype mat[4][4] )
{
  // initialise to the supplied matrix
  int i, j;
  for ( i=0; i<3; i++) for ( j=0; j<3; j++)
    rot()(i,j) = mat[i][j];
  for ( i=0; i<3; i++)
    trn()[i] = ftype( Util::mod( Util::intr(48.*mat[i][3]), 48 ) ) / 48.;
}

/*! Return formatted representation of the symmetry operator.
  \return The formatted text string, e.g. -x, -y+1/2, z. */
String Symop::format() const
{
  String s, t, xyz="xyz";
  for ( int i = 0; i < 3; i++ ) {
    t = "";
    for ( int j = 0; j < 3; j++ )
      if ( rot()(i,j) != 0.0 ) {
	t += ( rot()(i,j) > 0.0 ) ? "+" : "-";
	if ( Util::intr( fabs( rot()(i,j) ) ) != 1 )
	  t += String::rational( fabs( rot()(i,j) ), 24 );
	t += xyz[j];
      }
    if ( trn()[i] != 0.0 )
      t += String::rational( trn()[i], 24, true );
    s += t.substr( ( t[0] == '+' ) ? 1 : 0 );
    if ( i < 2 ) s+= ", ";
  }
  return s;
}


/*! Integerised symops are more efficient when handling integer
  coordinate types, e.g. HKL, Coord_grid. The rotation parts of the
  integerised symop are general and can be used for any recirpocal
  space data. The translation part is specific to an individual grid.
  \param symop The conventional symop.
  \param grid The specific grid.  */
Isymop::Isymop( const Symop& symop, const Grid& grid )
{
  for ( int i = 0; i < 3; i++ )
    for ( int j = 0; j < 3; j++ )
      rot()(i,j) = Util::intr( symop.rot()(i,j) );
  trn()[0] = Util::intr( grid.nu() * symop.trn()[0] );
  trn()[1] = Util::intr( grid.nv() * symop.trn()[1] );
  trn()[2] = Util::intr( grid.nw() * symop.trn()[2] );
}


/* Construct the encoded form of the symop. */
Symop_code::Symop_code( const Symop& op )
{ init( Isymop( op, Grid(24,24,24) ) ); }

/* Construct the encoded form of the symop from integerised symop with
   a grid (base) of (24,24,24).*/
Symop_code::Symop_code( const Isymop& op )
{ init( op ); }

/* Initialise the encoded form of the symop from integerised symop with
   a grid (base) of (24,24,24).*/
void Symop_code::init( const Isymop& op )
{
  // initialise to the supplied code
  int i, j, fac, code_r, code_t;
  code_r = code_t = 0;
  fac = 1;
  for ( i = 0; i < 3; i++ ) {
    code_t += ( Util::mod( op.trn()[i], 24 ) ) * fac;
    fac *= 24;
  }
  fac = 1;
  for ( i = 0; i < 3; i++ ) for ( j = 0; j < 3; j++ ) {
    code_r += ( Util::mod( op.rot()(i,j) + 1, 3 ) ) * fac;
    fac *= 3;
  }
  // xor to make identity zero
  code_ = ( ( code_r ^ 0x4064 ) << 16 ) + code_t;
}

/*! Construct a symmetry operator and initialise it to the matrix
  encoded in the given int.
  \param code The integer code. */
Symop Symop_code::symop() const
{
  Isymop iop = isymop();
  Symop op;
  for ( int i = 0; i < 3; i++ ) {
    op.rot()(i,0) = ftype( iop.rot()(i,0) );
    op.rot()(i,1) = ftype( iop.rot()(i,1) );
    op.rot()(i,2) = ftype( iop.rot()(i,2) );
    op.trn()[i] = ftype( iop.trn()[i] ) / 24.0;
  }
  return op;
}

/*! Construct an integerised symmetry operator and initialise it to
  the matrix encoded in the given int, with a grid (base) of (24,24,24).
  \param code The integer code. */
Isymop Symop_code::isymop() const
{
  Isymop op;
  // initialise rotation and translation
  int i, j, fac, code_r, code_t;
  code_t = code_trn();
  fac = 1;
  for ( i = 0; i < 3; i++ ) {
    op.trn()[i] = Util::mod( code_t/fac, 24 );
    fac *= 24;
  }
  code_r = ( code_rot() >> 16 ) ^ 0x4064;
  fac = 1;
  for ( i = 0; i < 3; i++ ) for ( j = 0; j < 3; j++ ) {
    op.rot()(i,j) = Util::mod( code_r/fac, 3 ) - 1;
    fac *= 3;
  }
  return op;
}

Symop_code Symop_code::code_rot() const
{ return Symop_code( code_ & 0xffff0000 ); }

Symop_code Symop_code::code_trn() const
{ return Symop_code( code_ & 0x0000ffff ); }


} // namespace clipper
