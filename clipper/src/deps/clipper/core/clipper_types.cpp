/* clipper_types.cpp: implementation file for clipper helper functions */
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


#include "clipper_types.h"

#include <sstream>


namespace clipper {


String::String( const char* str, const int l )
{
  std::ostringstream s;
  for ( int i = 0; i < l; i++ )
    s << str[i];
  *this = s.str();
}

String::String( const int i, const int w )
{ std::ostringstream s; s.width(w); s << i; *this = s.str(); }

String::String( const long i, const int w )
{ std::ostringstream s; s.width(w); s << i; *this = s.str(); }

String::String( const float f, const int w, const int p )
{ std::ostringstream s; s.width(w); s.precision(p); s << f; *this = s.str(); }

String::String( const double f, const int w, const int p )
{ std::ostringstream s; s.width(w); s.precision(p); s << f; *this = s.str(); }

std::vector<String> String::split(const String sep) const
{
  std::vector<String> splitstr;
  size_t tokbeg = 0, tokend = 0;
  while (1) {
    tokbeg = find_first_not_of(sep, tokend);
    if (tokbeg == String::npos) return splitstr;
    tokend = find_first_of(sep, tokbeg);
    if (tokend == String::npos) break;
    splitstr.push_back( substr(tokbeg, tokend-tokbeg) );
  }
  splitstr.push_back( substr(tokbeg) );
  return splitstr;
}

String String::trim() const
{
  String trimmed;
  int i, j;
  for ( i = 0; i < length(); i++ )
    if ( !isspace( (*this)[i] ) ) break;
  for ( j = length()-1; j >= 0; j-- )
    if ( !isspace( (*this)[j] ) ) break;
  return substr( i, j-i+1 );
}

String String::tail() const
{ return substr( rfind( '/' ) + 1 ); }

String String::head() const
{ return substr( 0, find( '/' ) ); }

String String::nohead() const
{
  size_t p = find( '/' );
  return ( p == String::npos ) ? "" : substr( p + 1 );
}

String String::notail() const
{
  size_t p = rfind( '/' );
  return ( p == String::npos ) ? "" : substr( 0, p );
}

String String::rational( const double f, const int b, const bool sign )
{
  std::ostringstream s;
  int n = Util::intr( fabs( b * f ) );
  int d = b;
  if (sign) s << ( ( f > 0 ) ? "+" : "-" );
  else      s << ( ( f > 0 ) ? "" : "-" );
  for ( int i = 5; i > 1; i-- ) {
    if ( ( n % i == 0 ) && ( d % i == 0 ) ) { n /= i; d /= i; }
  }
  s << n;
  if ( d != 1 ) s << "/" << d;
  return s.str();;
}

int String::i() const
{ std::istringstream s(*this); int i; s >> i; return i; }

long String::l() const
{ std::istringstream s(*this); long i; s >> i; return i; }

ftype32 String::f32() const
{ std::istringstream s(*this); float f; s >> f; return f; }

ftype64 String::f64() const
{ std::istringstream s(*this); double f; s >> f; return f; }

ftype String::f() const
{ std::istringstream s(*this); ftype f; s >> f; return f; }

ftype String::rational() const
{
  const String& s = (*this);
  int i;
  for ( i = 0; i < s.length(); i++ ) if ( s[i] == '/' ) break;
  if ( i == s.length() ) return ( s.f() );
  return ( String(s.substr(0,i)).f() / String(s.substr(i+1)).f() );
}


// template compilations

template class CLIPPER_IMEX Vec3<ftype32>;
template class CLIPPER_IMEX Mat33<ftype32>;
template class CLIPPER_IMEX Mat33sym<ftype32>;
template class CLIPPER_IMEX Matrix<ftype32>;

template class CLIPPER_IMEX Vec3<ftype64>;
template class CLIPPER_IMEX Mat33<ftype64>;
template class CLIPPER_IMEX Mat33sym<ftype64>;
template class CLIPPER_IMEX Matrix<ftype64>;

} // namespace clipper
