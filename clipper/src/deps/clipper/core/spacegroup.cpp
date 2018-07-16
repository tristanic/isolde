/* spacegroup.cpp: methods for spacegroup symmetry class */
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


#include "spacegroup.h"
#include "coords.h"
#include "clipper_instance.h"

#include <algorithm>
#include <sstream>


namespace clipper {


namespace data {

  // matrices
  Mat33<> mat_i   ( 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 );
  Mat33<> mat_inv (-1.0, 0.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0,-1.0 );
  Mat33<> mat_2z  (-1.0, 0.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0, 1.0 );
  Mat33<> mat_3z  ( 0.0,-1.0, 0.0, 1.0,-1.0, 0.0, 0.0, 0.0, 1.0 );
  Mat33<> mat_4z  ( 0.0,-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 );
  Mat33<> mat_6z  ( 1.0,-1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 );
  Mat33<> mat_2q  ( 0.0,-1.0, 0.0,-1.0, 0.0, 0.0, 0.0, 0.0,-1.0 );
  Mat33<> mat_2qq ( 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,-1.0 );
  Mat33<> mat_3abc( 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 );

  Vec3<> vec_0( 0.0, 0.0, 0.0 );
  Vec3<> vec_a( 0.5, 0.0, 0.0 );
  Vec3<> vec_b( 0.0, 0.5, 0.0 );
  Vec3<> vec_c( 0.0, 0.0, 0.5 );
  Vec3<> vec_n( 0.5, 0.5, 0.5 );
  Vec3<> vec_u( 0.25, 0.0, 0.0 );
  Vec3<> vec_v( 0.0, 0.25, 0.0 );
  Vec3<> vec_w( 0.0, 0.0, 0.25 );
  Vec3<> vec_d( 0.25, 0.25, 0.25 );

  Vec3<> vec_A( 0.0, 0.5, 0.5 );
  Vec3<> vec_B( 0.5, 0.0, 0.5 );
  Vec3<> vec_C( 0.5, 0.5, 0.0 );
  Vec3<> vec_I( 0.5, 0.5, 0.5 );
  Vec3<> vec_R( 0.666667, 0.333333, 0.333333 );
  Vec3<> vec_S( 0.333333, 0.666667, 0.333333 );
  Vec3<> vec_T( 0.333333, 0.333333, 0.666667 );
  Vec3<> vec_H( 0.666667, 0.333333, 0.0 );
  Vec3<> vec_F1( 0.0, 0.5, 0.5 );
  Vec3<> vec_F2( 0.5, 0.0, 0.5 );

} // namespace data


// spacegroup description object

void Spgr_descr::Symop_codes::init_hall( const String& symb )
{
  // interpret Hall symbol
  using namespace data;

  Symop_codes& ops = (*this);
  std::vector<String> tokens;
  String token, sym, chb;
  char c, inv, lat;
  int i, j, k, l;

  // separate the Hall symbol from the change of basis
  tokens = symb.split("()");
  if ( tokens.size() > 0 ) sym = tokens[0].trim();
  if ( tokens.size() > 1 ) chb = tokens[1].trim();

  // now separate the parts of the Hall symbol
  tokens = sym.split(" \t");

  // first part: lattice and inversion
  inv = lat = ' ';
  token = tokens[0];
  for ( j = 0; j < token.length(); j++ ) {
    c = toupper( token[j] );
    if ( c == '-' ) inv = c;
    else            lat = c;
  }

  // next 1-3 parts: generating matrices
  int nmat = tokens.size()-1;
  std::vector<char> invop(nmat,' '), order(nmat,' '), screw(nmat,' '),
    axis1(nmat,' '), axis2(nmat,' ');
  std::vector<String> trans(nmat);
  for ( i = 0; i < nmat; i++ ) {
    token = tokens[i+1];
    for ( j = 0; j < token.length(); j++ ) {
      c = tolower( token[j] );
      if ( c == '-' ) invop[i] = c;
      else if ( c >= '1' && c <= '6' )
        if ( order[i] == ' ' ) order[i] = c;
        else                   screw[i] = c;
      else if ( c >= 'x' && c <= 'z' ) axis1[i] = c;
      else if ( c >= '"' && c <= '*' ) axis2[i] = c;
      else if ( c >= 'a' && c <= 'w' ) trans[i] += c;
    }
  }

  // now interpret all the default conventions
  // default first axis to z
  if ( nmat >= 1 ) {
    if ( axis1[0] == ' ' ) axis1[0] = 'z';
  }
  // default second axis on basis of first
  if ( nmat >= 2 ) {
    if ( axis1[1] == ' ' ) {
      if ( order[1] == '2' ) {
        if        ( order[0] == '2' || order[0] == '4' ) {
          if ( axis2[1] == ' ' ) axis1[1] = 'x';
        } else if ( order[0] == '3' || order[0] == '6' ) {
          axis1[1] = axis1[0];
          if ( axis2[1] == ' ' ) axis2[1] = '\'';
        }
      } else if ( order[1] == '3' ) {
        if ( order[0] == '2' || order[0] == '4' ) {
          if ( axis2[1] == ' ' ) axis2[1] = '*';
        }
      }
    }
  }
  // default third axis (not much choice here)
  if ( nmat >= 3 ) {
    if ( axis1[2] == ' ' ) {
      if ( order[2] == '3' ) {
        if ( axis2[2] == ' ' ) axis2[2] = '*';
      }
    }
  }

  // now check we have everything
  for ( i = 0; i < nmat; i++ ) {
    // add fake z axes for non-axis ops
    if ( axis1[i] == ' ' ) {
      if ( order[i] == '1' ) axis1[i] = 'z';  // fake axis
      if ( axis2[i] != ' ' ) axis1[i] = 'z';  // fake axis
    }
    // error messages
    if ( axis1[i] == ' ' ) Message::message(
      Message_fatal("Spacegroup: Missing x/y/z in Hall symb:"+tokens[i+1]) );
    if ( order[i] == ' ' ) Message::message(
      Message_fatal("Spacegroup: Missing order in Hall symb:"+tokens[i+1]) );
  }

  // add identity and inversion
  ops.clear();
  ops.push_back( Symop_code::identity() );
  if ( inv == '-' ) ops.push_back(Symop_code(Symop(RTop<>(mat_inv))));

  // now make the generator matrices
  Mat33<> mat, matperm;
  Vec3<> vec;
  for ( i = 0; i < nmat; i++ ) {
    // make matrix part
    mat = mat_i;
    if ( order[i] == '2' && axis2[i] == ' ' ) mat = mat_2z;
    if ( order[i] == '2' && axis2[i] == '\'') mat = mat_2q;
    if ( order[i] == '2' && axis2[i] == '"' ) mat = mat_2qq;
    if ( order[i] == '3' && axis2[i] == ' ' ) mat = mat_3z;
    if ( order[i] == '3' && axis2[i] == '*' ) mat = mat_3abc;
    if ( order[i] == '4' && axis2[i] == ' ' ) mat = mat_4z;
    if ( order[i] == '6' && axis2[i] == ' ' ) mat = mat_6z;
    // inverse (improper)
    if ( invop[i] == '-' ) mat = mat_inv * mat;
    // axis permutation
    if      ( axis1[i] == 'x' ) j = 0;
    else if ( axis1[i] == 'y' ) j = 1;
    else                        j = 2;
    for ( k = 0; k < 3; k++ )
      for ( l = 0; l < 3; l++ )
        matperm( k, l ) = mat( (k+2-j)%3, (l+2-j)%3 );
    // make translation part
    vec = vec_0;
    for ( k = 0; k < trans[i].length(); k++ ) {
      if      ( trans[i][k] == 'a' ) vec = vec + vec_a;
      else if ( trans[i][k] == 'b' ) vec = vec + vec_b;
      else if ( trans[i][k] == 'c' ) vec = vec + vec_c;
      else if ( trans[i][k] == 'n' ) vec = vec + vec_n;
      else if ( trans[i][k] == 'u' ) vec = vec + vec_u;
      else if ( trans[i][k] == 'v' ) vec = vec + vec_v;
      else if ( trans[i][k] == 'w' ) vec = vec + vec_w;
      else if ( trans[i][k] == 'd' ) vec = vec + vec_d;
    }
    // screw translation
    if ( screw[i] != ' ' )
      vec[j] += ftype( screw[i] - '0' ) / ftype( order[i] - '0' );
    // store the matrix
    ops.push_back( Symop_code( Symop( RTop<>( matperm, vec ) ) ) );
  }

  // add lattice centering
  if (lat=='A') ops.push_back(Symop_code(Symop(RTop<>(mat_i,vec_A))));
  if (lat=='B') ops.push_back(Symop_code(Symop(RTop<>(mat_i,vec_B))));
  if (lat=='C') ops.push_back(Symop_code(Symop(RTop<>(mat_i,vec_C))));
  if (lat=='I') ops.push_back(Symop_code(Symop(RTop<>(mat_i,vec_I))));
  if (lat=='R') ops.push_back(Symop_code(Symop(RTop<>(mat_i,vec_R))));
  if (lat=='S') ops.push_back(Symop_code(Symop(RTop<>(mat_i,vec_S))));
  if (lat=='T') ops.push_back(Symop_code(Symop(RTop<>(mat_i,vec_T))));
  if (lat=='H') ops.push_back(Symop_code(Symop(RTop<>(mat_i,vec_H))));
  if (lat=='Q') ops.push_back(Symop_code(Symop(RTop<>(mat_i,vec_S))));
  if (lat=='F') {
    ops.push_back(Symop_code(Symop(RTop<>(mat_i,vec_F1))));
    ops.push_back(Symop_code(Symop(RTop<>(mat_i,vec_F2))));
  }

  // apply the change of basis operator
  RTop_frac cbop( RTop_frac::identity() );
  if ( chb.find( ',' ) != String::npos ) {
    cbop = RTop_frac( chb );
  } else {
    std::vector<String> t = chb.split(" ");
    for ( i = 0; i < Util::min( int(t.size()), 3 ); i++ )
      cbop.trn()[i] = ftype( t[i].i() ) / 12.0;
  }
  RTop_frac cbopi( cbop.inverse() );
  for ( i = 0; i < ops.size(); i++ )
    ops[i] = Symop_code( Symop( cbop*ops[i].symop()*cbopi ) );
}

void Spgr_descr::Symop_codes::init_symops( const String& symb )
{
  Symop_codes& ops = (*this);
  std::vector<String> symops = symb.split(";");
  for ( int i = 0; i < symops.size(); i++ )
    ops.push_back( Symop_code( Symop( RTop_frac( symops[i] ) ) ) );
}

Spgr_descr::Symop_codes Spgr_descr::Symop_codes::expand() const
{
  int i, j, l, size;
  Symop_code k;

  const Symop_codes& generators = (*this);
  Symop_codes ops;
  ops.push_back( Symop_code::identity() );  // identity is compulsory

  // generate all the symops
  do {
    size = ops.size();
    for ( i = 0; i < generators.size(); i++ )
      if ( generators[i] != Symop_code::identity() ) {
	for ( j = 0; j < size; j++ ) {
	  k = Symop_code( Isymop( generators[i].isymop() * ops[j].isymop() ) );
	  for ( l = 0; l < ops.size(); l++ )
	    if ( ops[l] == k ) break;
	  if ( l == ops.size() ) ops.push_back( k );
	}
      }
  } while ( ops.size() > size );
  return ops;
}

Spgr_descr::Symop_codes Spgr_descr::Symop_codes::primitive_noninversion_ops() const
{
  Symop_codes pops = primitive_ops();
  if ( inversion_ops().size() > 1 ) {
    Symop_codes nops;
    for ( int i = 0; i < pops.size(); i++ )
      if ( pops[i].symop().rot().det() > 0.0 )
	nops.push_back( pops[i] );
    pops = nops;
  }
  return pops;
}

Spgr_descr::Symop_codes Spgr_descr::Symop_codes::inversion_ops() const
{
  const Symop_codes& ops = (*this);
  Symop_codes pops;
  Symop_code invop = Symop_code( Symop( RTop<>( data::mat_inv ) ) );
  pops.push_back( Symop_code::identity() );
  int i;
  for ( i = 0; i < ops.size(); i++ )
    if ( ops[i].code_rot() == invop ) { pops.push_back( ops[i] ); break; }
  return pops;
}

Spgr_descr::Symop_codes Spgr_descr::Symop_codes::primitive_ops() const
{
  const Symop_codes& ops = (*this);
  Symop_codes pops;
  int i, j;
  pops.push_back( Symop_code::identity() );
  for ( i = 0; i < ops.size(); i++ ) {
    for ( j = 0; j < pops.size(); j++ )
      if ( ops[i].code_rot() == pops[j].code_rot() ) break;
    if ( j == pops.size() ) pops.push_back( ops[i] );
  }
  return pops;
}

Spgr_descr::Symop_codes Spgr_descr::Symop_codes::centering_ops() const
{
  const Symop_codes& ops = (*this);
  Symop_codes cops;
  for ( int i = 0; i < ops.size(); i++ )
    if ( ops[i].code_rot() == Symop_code::identity() )
      cops.push_back( ops[i] );
  return cops;
}

Spgr_descr::Symop_codes Spgr_descr::Symop_codes::laue_ops() const
{
  const Symop_codes& ops = (*this);
  Symop_codes lops;
  lops.push_back( Symop_code( Symop( RTop<>( data::mat_inv ) ) ) );
  for ( int i = 0; i < ops.size(); i++ )
    lops.push_back( ops[i].code_rot() );
  return lops.expand();
}

Spgr_descr::Symop_codes Spgr_descr::Symop_codes::pgrp_ops() const
{
  const Symop_codes& ops = (*this);
  Symop_codes lops;
  for ( int i = 0; i < ops.size(); i++ )
    lops.push_back( ops[i].code_rot() );
  return lops.expand();
}

Spgr_descr::Symop_codes Spgr_descr::Symop_codes::patterson_ops() const
{
  const Symop_codes& ops = (*this);
  Symop_codes lops;
  lops.push_back( Symop_code( Symop( RTop<>( data::mat_inv ) ) ) );
  for ( int i = 0; i < ops.size(); i++ )
    if ( ops[i].code_rot() == Symop_code::identity() )
      lops.push_back( ops[i] );
    else
      lops.push_back( ops[i].code_rot() );
  return lops.expand();
}

Spgr_descr::Symop_codes Spgr_descr::Symop_codes::generator_ops() const
{
  Symop_codes ops = expand();
  std::sort( ops.begin(), ops.end() );
  Symop_codes cens = ops.centering_ops();
  Symop_codes prms = ops.primitive_ops();
  Symop_codes gens;
  Symop_codes gend = gens.expand();

  // first make the centering generators
  for ( int i = 1; i < cens.size(); i++ ) {
    for ( int j = 0; j < gend.size(); j++ )
      if ( cens[i] == gend[j] ) goto skip1;
    gens.push_back( cens[i] );
    gend = gens.expand();
    if ( gend.size() == cens.size() ) break;  // optional optimisation
  skip1:;
  }
  int ncen = gens.size();

  // now add the rest of the generators
  for ( int i = 1; i < prms.size(); i++ ) {
    for ( int j = 0; j < gend.size(); j++ )
      if ( prms[i] == gend[j] ) goto skip2;
    gens.push_back( prms[i] );
    gend = gens.expand();
    if ( gend.size() == ops.size() ) break;  // optional optimisation
  skip2:;
  }

 back:  // finally remove any redundent ops
  for ( int i = ncen; i < gens.size(); i++ ) {
    Symop_codes genp = gens;
    genp.erase( genp.begin() + i );
    if ( genp.expand().size() == ops.size() ) {
      gens = genp;
      goto back;
    }
  }

  return gens;  // return result
}

Spgr_descr::Symop_codes Spgr_descr::Symop_codes::product( const Symop_codes& ops2 ) const
{
  Symop_codes ops1 = (*this);  // copy first list, implying identity
  int i, j, n1, n2;
  n1 = ops1.size();
  n2 = ops2.size();
  for ( j = 0; j < n2; j++ )
    if ( ops2[j] != Symop_code::identity() )  // skip identity, if present
      for ( i = 0; i < n1; i++ )
	ops1.push_back(Symop_code(Isymop(ops1[i].isymop()*ops2[j].isymop())));
  return ops1;
}

unsigned int Spgr_descr::Symop_codes::hash() const
{
  Symop_codes data = expand();
  std::sort( data.begin(), data.end() );
  unsigned int polynomial = 0x04c11db7;
  unsigned int remainder  = 0xffffffff;
  unsigned int datum;
  for ( int word = 0; word < data.size(); word++ ) {
    datum = data[word];
    remainder ^= datum;
    for ( int bit = 0; bit < 32; bit++ ) {
      if ( remainder & 0x80000000 )
        remainder = (remainder << 1) ^ polynomial;
      else
        remainder = (remainder << 1);
    }
  }
  return remainder;
}


// Spacegroup desciption object

char Spgr_descr::pref_12 = '1';
char Spgr_descr::pref_hr = 'H';


/*! Construct a null description spacegroup. The result is initialised
  to an invalid spacegroup code. */
Spgr_descr::Spgr_descr()
{ hash_ = 0; }

/*! Construct a spacegroup description from a text description, i.e. a
  symbol or operators. This may be one of the following:
  - Hall symbol, e.g. " P 2ac 2ab"
  - H-M symbol, e.g. "P 21 21 21"
  - Number, e.g. "19"
  - List of symmetry operators separated by semicolons, e.g.
  "x,y,z;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2"

  It is best to specify the type of symbol being used, however if this
  parameter is omitted a guess will be made. Unfortunately, Hall and
  H-M symbols may be ambiguous. Any ambiguity may be removed by
  appending parentheses "()" to the end of the Hall symbol, otherwise
  the symbol will be interpreted as an H-M symbol, and a Hall symbol
  if that fails.

  H-M symbols and spacegroup numbers may correspond to 2 different
  entries in international tables. The choice between 2 origin
  settings or hexagonal/rhomohedral settings is made using the
  set_preferred() method.

  \param name The spacegroup symbol or operators.
  \param type The type of symbol: Spacegroup::Symops, Spacegroup::Hall,
  Spacegroup::HM, Spacegroup::XHM, Spacegroup::Number */
Spgr_descr::Spgr_descr( const String& symb, TYPE type )
{
  using clipper::data::sgdata;

  String symbx = symb.trim();

  // try and guess symbol type (don't do this!)
  if ( type == Unknown ) {
    if      ( symbx.find_first_of( "()" ) != String::npos )
      type = Hall;
    else if ( symbx.find_first_of( ":" ) != String::npos )
      type = XHM;
    else if ( symbx.find_first_of( "," ) != String::npos )
      type = Symops;
    else if ( symbx.find_first_of( "ABCFHIPRSTQ" ) == String::npos )
      type = Number;
    // otherwise still unknown: try HM then Hall.
  }

  // now make the ops
  int i;
  Symop_codes ops;

  if ( type == Symops ) {
    ops.init_symops( symbx );
  } else if ( type == Hall ) {
    ops.init_hall( symbx );
  } else if ( type == XHM ) {
    char ext = ' ';
    size_t c = symbx.find_first_of( ":" );
    if ( c != String::npos ) {
      if ( c+1 < symbx.length() ) ext = symbx[c+1];
      symbx = symbx.substr(0,c);
      symbx = symbx.trim();
    }
    for ( i = 0; i < data::sgdata_size; i++ )
      if ( symbx == String( sgdata[i].hm ) && sgdata[i].ext == ext ) break;
    if ( i == data::sgdata_size )
      Message::message( Message_fatal( "Spgr_descr: No such HM symbol" ) );
    ops.init_hall( String( sgdata[i].hall ) );
  } else if ( type == HM ) {
    for ( i = 0; i < data::sgdata_size; i++ )
      if ( symbx == String( sgdata[i].hm ) &&
	   ( sgdata[i].ext == pref_12 || sgdata[i].ext == pref_hr ||
	     sgdata[i].ext == ' ' ) ) break;
    if ( i == data::sgdata_size )
      Message::message( Message_fatal( "Spgr_descr: No such HM symbol" ) );
    ops.init_hall( String( sgdata[i].hall ) );
  } else if ( type == Number ) {
    int num = symbx.i();
    for ( i = 0; i < data::sgdata_size; i++ )
      if ( num == sgdata[i].num &&
	   ( sgdata[i].ext == pref_12 || sgdata[i].ext == pref_hr ||
	     sgdata[i].ext == ' ' ) ) break;
    if ( i == data::sgdata_size )
      Message::message( Message_fatal( "Spgr_descr: No such SG number" ) );
    ops.init_hall( String( sgdata[i].hall ) );
  } else {
    for ( i = 0; i < data::sgdata_size; i++ )  // try H-M then Hall.
      if ( symbx == String( sgdata[i].hm ) &&
	   ( sgdata[i].ext == pref_12 || sgdata[i].ext == pref_hr ||
	     sgdata[i].ext == ' ' ) ) break;
    if ( i != data::sgdata_size )
      ops.init_hall( String( sgdata[i].hall ) );
    else
      ops.init_hall( symbx );
  }

  // store the hash and generators
  hash_ = ops.hash();
  generators_ = ops.generator_ops();
}

/*! See previous constuctor.
  \param num The spacegroup number. */
Spgr_descr::Spgr_descr( const int& num )
{
  using data::sgdata;
  Symop_codes ops;
  int i;
  for ( i = 0; i < data::sgdata_size; i++ )
    if ( num == sgdata[i].num &&
	 ( sgdata[i].ext == pref_12 || sgdata[i].ext == pref_hr ||
	   sgdata[i].ext == ' ' ) ) break;
  if ( i == data::sgdata_size )
    Message::message( Message_fatal( "Spgr_descr: No such SG number" ) );
  ops.init_hall( String( sgdata[i].hall ) );
  ops = ops.expand();

  // store the hash and generators
  hash_ = ops.hash();
  generators_ = ops.generator_ops();
}

/*! This is not normally used, except in conjunction with
   Spgr_desc::generator_ops() to derive one group from another. */
Spgr_descr::Spgr_descr( const Symop_codes& ops )
{
  // store the hash and generators
  hash_ = ops.hash();
  generators_ = ops.generator_ops();
}

/*! The spacegroup number is only available if the spacegroup exists in
   the internal table, see Hall & Grosse-Kunstleve.
   \return The spacegroup number, or 0 if unavailable. */
int Spgr_descr::spacegroup_number() const
{
  for ( int i = 0; i < data::sgdata_size; i++ )
    if ( data::sgdata[i].sghash == hash_ ) return data::sgdata[i].num;
  return 0;
}

/*! The Hall symbol is only available if the spacegroup exists in
   the internal table, see Hall & Grosse-Kunstleve.
   \return The Hall symbol, or "Unknown" if unavailable. */
String Spgr_descr::symbol_hall() const
{
  for ( int i = 0; i < data::sgdata_size; i++ )
    if ( data::sgdata[i].sghash == hash_ ) return String(data::sgdata[i].hall);
  return "Unknown";
}

/*! The H-M symbol is only available if the spacegroup exists in
   the internal table, see Hall & Grosse-Kunstleve.
   \return The H-M symbol, or "Unknown" if unavailable. */
String Spgr_descr::symbol_hm() const
{
  for ( int i = 0; i < data::sgdata_size; i++ )
    if ( data::sgdata[i].sghash == hash_ ) return String(data::sgdata[i].hm);
  return "Unknown";
}

/*! The extended H-M symbol is only available if the spacegroup exists in
   the internal table, see Hall & Grosse-Kunstleve.
   \return The extended H-M symbol, or "Unknown" if unavailable. */
String Spgr_descr::symbol_xhm() const
{
  for ( int i = 0; i < data::sgdata_size; i++ )
    if ( data::sgdata[i].sghash == hash_ ) {
      String xhm( data::sgdata[i].hm );
      if ( data::sgdata[i].ext != ' ' )
	xhm = xhm + " :" + data::sgdata[i].ext;
      return xhm;
    }
  return "Unknown";
}

/*! The extension H-M symbol is only available if the spacegroup exists in
   the internal table, see Hall & Grosse-Kunstleve.
   \return The extension H-M symbol, or "" */
String Spgr_descr::symbol_hm_ext() const
{
  String ext = "";
  for ( int i = 0; i < data::sgdata_size; i++ )
    if ( data::sgdata[i].sghash == hash_ )
      if ( data::sgdata[i].ext != ' ' )
	return ext + data::sgdata[i].ext;
  return ext;
}

/*! Sets the preferred origin or setting for initialising all
  Spgr_descr objects using H-M symbols or Spacegroup numbers. cctbx
  uses origin choice '1' by default, CCP4 uses '2'. Both packages use
  'H' in preference to 'R'. Preferred values are stored for
  both. Defaults are '1' and 'H'.

  CCP4 users may wish to add the following before using H-M codes or numbers.
  \code
  Spgr_descr::set_preferred('2');
  \endcode

  \param c Either '1' or '2', 'H' or 'R'. */
void Spgr_descr::set_preferred( const char& c )
{
  if ( c == '1' || c == '2' ) pref_12 = c;
  if ( c == 'H' || c == 'R' ) pref_hr = c;
}



// Spacegroup cache object

Mutex Spgr_cacheobj::mutex = Mutex();

struct Compare_grid{ bool operator() ( const Vec3<>& c1, const Vec3<>& c2 ) const { return (c1[0]*c1[1]*c1[2]+0.001*c1[1]+0.00001*c1[0]) <= (c2[0]*c2[1]*c2[
2]+0.001*c2[1]+0.00001*c2[0]); } };

bool reci_asu( const Spgr_descr::Symop_codes& ops, data::ASUfn asufn )
{
  HKL hkl, equ;
  int i;

  // make integerised symops
  std::vector<Isymop> symops( ops.size() );
  for ( i = 0; i < ops.size(); i++ ) symops[i] = ops[i].isymop();

  // now test asu for uniqueness and completeness
  for ( hkl.h() = -2; hkl.h() <= 2; hkl.h()++ )
    for ( hkl.k() = -2; hkl.k() <= 2; hkl.k()++ )
      for ( hkl.l() = -2; hkl.l() <= 2; hkl.l()++ ) {
        if ( asufn( hkl.h(), hkl.k(), hkl.l() ) ) {
          for ( i = 0; i < symops.size(); i++ ) {
            equ = hkl.transform( symops[i] );
            if ( equ != hkl && asufn(equ.h(),equ.k(),equ.l()) ) break;
            equ = -equ;
            if ( equ != hkl && asufn(equ.h(),equ.k(),equ.l()) ) break;
          }
          if ( i != symops.size() ) return false;
        } else {
          for ( i = 0; i < symops.size(); i++ ) {
            equ = hkl.transform( symops[i] );
            if ( asufn(equ.h(),equ.k(),equ.l()) ) break;
            equ = -equ;
            if ( asufn(equ.h(),equ.k(),equ.l()) ) break;
          }
          if ( i == symops.size() ) return false;
        }
      }
  return true;
}

Vec3<> real_asu( const Spgr_descr::Symop_codes& ops )
{
  int i, j, sym, nasu;

  // make integerised symops
  std::vector<Isymop> symops( ops.size() );
  for ( i = 0; i < ops.size(); i++ ) symops[i] = ops[i].isymop();

  // classify each grid point by a unique 'ASU number'
  Coord_grid c;
  Grid_sampling cgrid(24,24,24);
  int symmap[13824], tstmap[13824];
  for ( i = 0; i < cgrid.size(); i++ ) symmap[i] = -1;
  nasu = 0;
  for ( c = Coord_grid(0,0,0); !c.last(cgrid); c.next(cgrid) ) {
    i = c.index(cgrid);
    if ( symmap[i] == -1 ) {
      for ( sym = 0; sym < symops.size(); sym++ )
        symmap[ c.transform(symops[sym]).unit(cgrid).index(cgrid) ] = nasu;
      nasu++;
    }
  }

  // identify trigonal/hexagonal groups
  bool trighex = false;
  for ( sym = 0; sym < symops.size(); sym++ )
    for ( i = 0; i < 3; i++ ) {
      int c = 0, r = 0;     // trigonal if any row/col adds to 2
      for ( j = 0; j < 3; j++ ) {
        c += abs( symops[sym].rot()(i,j) );
        r += abs( symops[sym].rot()(j,i) );
        trighex = trighex || ( c == 2 ) || ( r == 2 );
      }
    }

  // now set search ASU search grid, dependent on symmetry
  std::vector<ftype> gridlim;
  const ftype d = 0.0001;
  if ( trighex ) {
    ftype lim[] = { 1./12-d, 1./6-d, 1./6+d, 1./3-d, 1./3+d,
                    1./2-d, 1./2+d, 2./3+d, 1.-d };
    gridlim = std::vector<ftype>( &lim[0], &lim[ sizeof(lim)/sizeof(ftype) ] );
  } else {
    ftype lim[] = { 1./8+d, 1./4-d, 1./4+d, 1./2-d, 1./2+d, 1.-d };
    gridlim = std::vector<ftype>( &lim[0], &lim[ sizeof(lim)/sizeof(ftype) ] );
  }

  // make a sorted list of grids
  std::vector<Vec3<> > grids;
  for ( int gui = 0; gui < gridlim.size(); gui++ )
    for ( int gvi = 0; gvi < gridlim.size(); gvi++ )
      for ( int gwi = 0; gwi < gridlim.size(); gwi++ )
        grids.push_back( Vec3<>(gridlim[gui],gridlim[gvi],gridlim[gwi]) );
  std::sort( grids.begin(), grids.end(), Compare_grid() );

  // now find smallest viable grid
  for ( j = 0; j < grids.size(); j++ ) {
    Grid mgrid = Grid( Util::intc(grids[j][0]*cgrid.nu()),
                       Util::intc(grids[j][1]*cgrid.nv()),
                       Util::intc(grids[j][2]*cgrid.nw()) );
    if ( mgrid.size() >= nasu ) {    // is grid big enough to be ASU?
      for ( i = 0; i < nasu; i++ ) tstmap[i] = 0;  // if so, check for all pts
      for ( c = Coord_grid(0,0,0); !c.last(mgrid); c.next(mgrid) )
        tstmap[ symmap[ c.index(cgrid) ] ] = 1;
      for ( i = 0; i < nasu; i++ ) if ( tstmap[i] == 0 ) break;
      if ( i == nasu ) break;        // found a full asu, so we're done
    }
  }
  return grids[j];
}

Spgr_cacheobj::Spgr_cacheobj( const Key& spgr_cachekey )
{
  spgr_cachekey_ = spgr_cachekey;

  // tidy the ops
  Spgr_descr::Symop_codes ops = spgr_cachekey.generator_ops().expand();
  std::sort( ops.begin(), ops.end() );
  Spgr_descr::Symop_codes pops = ops.primitive_noninversion_ops();
  Spgr_descr::Symop_codes iops = ops.inversion_ops();
  Spgr_descr::Symop_codes cops = ops.centering_ops();
  std::sort( pops.begin(), pops.end() );
  std::sort( iops.begin(), iops.end() );
  std::sort( cops.begin(), cops.end() );
  ops = pops;
  ops = ops.product( iops );
  ops = ops.product( cops );
  nsym  = ops.size();
  nsymn = pops.size();
  nsymi = iops.size();
  nsymc = cops.size();
  nsymp = nsymn*nsymi;

  // consistency check (redundent)
  if ( ops.hash() != spgr_cachekey.hash() )
    Message::message( Message_fatal( "Spgr_cacheobj: symops fail" ) );

  // Laue group
  unsigned int lghash = ops.laue_ops().hash();
  // first guess from the hash
  for ( lgrp = 0; lgrp < data::lgdata_size; lgrp++ )
    if ( data::lgdata[lgrp].lghash == lghash ) break;
  if ( lgrp == data::lgdata_size ) lgrp = 0;
  // if that fails, try all the ASU functions
  if ( !reci_asu( pops, data::lgdata[lgrp].asufn ) ) {
    std::ostringstream s;
    s << "Spacegroup_registry: ASU warning, LGhash=0x";
    s.width( 8 ); s.fill( '0' ); s.flags( s.hex | s.right ); s << lghash;
    Message::message( Message_warn( s.str() ) );
    for ( lgrp = 0; lgrp < data::lgdata_size; lgrp++ )
      if ( reci_asu( pops, data::lgdata[lgrp].asufn ) ) break;
    if ( lgrp == data::lgdata_size )
      Message::message( Message_fatal( "Spacegroup_registry: ASU fail" ) );
  }

  // real ASU
  asu_min_ = Vec3<>( -0.0001, -0.0001, -0.0001 );
  asu_max_ = real_asu( ops );

  // now make real symop lists
  for ( int i = 0; i < ops.size(); i++ ) {
    symops.push_back( ops[i].symop() );
    isymops.push_back( ops[i].isymop() );
  }
}

bool Spgr_cacheobj::matches( const Key& spgr_cachekey ) const
{ return spgr_cachekey_.hash() == spgr_cachekey.hash(); }

String Spgr_cacheobj::format() const
{ return spgr_cachekey_.symbol_hall(); }



// spacegroup object

/*! Construct null or P1 spacegroup. This is faster than the normal
  constructor.
  \param type Spacegroup::Null or Spacegroup::P1 */
Spacegroup::Spacegroup( TYPE type )
{ if ( type == P1 ) init( Spgr_descr( "x,y,z", Spgr_descr::Symops ) ); }

/*! Construct a spacegroup and initialise with a spacegroup description.
  \param spgr_descr The spacegroup description. */
Spacegroup::Spacegroup( const Spgr_descr& spgr_descr )
{ init( spgr_descr ); }

/*! Initialise the spacegroup.
  \param spgr_descr The spacegroup description. */
void Spacegroup::init( const Spgr_descr& spgr_descr )
{
  // init spgr descr
  hash_ = spgr_descr.hash();
  generators_ = spgr_descr.generator_ops();
  // init cache entry
  cacheref = ClipperInstantiator::instance().spacegroup_cache().cache( spgr_descr );
  symops  = &(cacheref.data().symops[0]);
  isymops = &(cacheref.data().isymops[0]);
  nsym  = cacheref.data().nsym;
  nsymn = cacheref.data().nsymn;
  nsymi = cacheref.data().nsymi;
  nsymc = cacheref.data().nsymc;
  nsymp = cacheref.data().nsymp;
  asufn = data::lgdata[cacheref.data().lgrp].asufn;
}

/*! \return true if the object has not been initalised. */
bool Spacegroup::is_null() const
{ return cacheref.is_null(); }

/*! The number of rotational operators parallel to the specified axis
  is returned.
  \param  axis The axis, A, B or C.
  \return The order of the axis. */
int Spacegroup::order_of_symmetry_about_axis( const AXIS axis ) const
{
  int n = 0;
  int a0 = int( axis );
  int a1 = (a0+1)%3;
  int a2 = (a0+2)%3;
  for ( int i = 0; i < nsymp; i++ )
    if ( symops[i].rot().det() > 0.0 )
      if ( fabs( symops[i].rot()(a0,a1) ) + fabs( symops[i].rot()(a1,a0) ) +
	   fabs( symops[i].rot()(a0,a2) ) + fabs( symops[i].rot()(a2,a0) ) +
	   fabs( symops[i].rot()(a0,a0) - 1.0 ) < 1.0e-6 ) n++;
  return n;
}

/*! The reflection class describes the type of a reflection in a given
  spacegroup, including centricity, systematic absence, phase
  restriction, and multiplicity.

  This is a shortcut to constructing an HKL_class from the
  spacegroup and HKL.
  \param hkl The reflection HKL */
HKL_class Spacegroup::hkl_class( const HKL& hkl ) const
{ return HKL_class( *this, hkl ); }

/*! The reciprocal ASU is chosen from one of 47 optimised functions.
  \param hkl The HKL to test.
  \return true if the HKL is in the ASU. */
bool Spacegroup::recip_asu( const HKL& hkl ) const
{ return asufn( hkl.h(), hkl.k(), hkl.l() ); }

int Spacegroup::product_op( const int& s1, int& s2 ) const
{
  Symop mat( symops[s1] * symops[s2] );
  for ( int s3 = 0; s3 < nsym; s3++ )
    if ( mat.equals( symops[s3], 0.001 ) ) return s3;
  Message::message( Message_fatal("Spacegroup: Internal spacegroup error - missing product") ); return -1;
}

int Spacegroup::inverse_op( const int& s1 ) const
{
  for ( int s2 = 0; s2 < nsym; s2++ )
    if ( symops[0].equals( symops[s1] * symops[s2], 0.001 ) ) return s2;
  Message::message( Message_fatal("Spacegroup: Internal spacegroup error - missing inverse") ); return -1;
}

/*! The map ASU is an oblong which contains at least one assymetric
  unit. It is guaranteed to be contained withing the unit box. The
  lower limit is always 0,0,0.
  \return Fractional coordinate of the upper bound of the ASU. */
Coord_frac Spacegroup::asu_max() const
{ return Coord_frac( cacheref.data().asu_max_ ); }

/*! The map ASU is an oblong which contains at least one assymetric
  unit. It is guaranteed to be contained withing the unit box. The
  lower limit is always 0,0,0.
  \return Fractional coordinate of the lower bound of the ASU. */
Coord_frac Spacegroup::asu_min() const
{ return Coord_frac( cacheref.data().asu_min_ ); }

/*! Test if hand-change is possible.
  \return true if a change of hand preserves the spacegroup. */
bool Spacegroup::invariant_under_change_of_hand() const
{
  for ( int k=0; k<nsym; k++ )
    for ( int j=0; j<3; j++ )
      if ( isymops[k].rot()(j,j) == 1 )
	if ( isymops[k].trn()[j] != 0 && isymops[k].trn()[j] != 12 )
	  return false;
  return true;
}

/*! \return The Laue group symbol. i.e. one of
  -1, 2/m, 2/mmm, -3, -3m, 4/m, 4/mmm, 6/m, 6/mmm, m-3, m-3m */
String Spacegroup::symbol_laue() const
{ return String( data::lgdata[cacheref.data().lgrp].lgname ); }

void Spacegroup::debug() const
{
  // print the members of the class
  int k;
  std::cout << spacegroup_number() << " " << nsym << " " << nsymp << " " << symbol_hall() << "\n";
  for (k=0; k<nsym; k++)
    std::cout << k << ": " << symop(k).format() << "\n";
}


} // namespace clipper
