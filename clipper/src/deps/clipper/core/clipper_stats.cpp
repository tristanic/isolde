/* clipper_stats.cpp: implementation file for clipper helper functions */
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


#include "clipper_stats.h"


namespace clipper {


ftype Histogram::sum() const
{
  ftype s = 0.0;
  for ( int i = 0; i < data.size(); i++ ) s += data[i];
  return s;
}

ftype Histogram::y( const ftype& x ) const
{
  ftype xi = indexf( truncate( x ) ) - 0.5;
  int i = Util::intf( xi );
  ftype xf = xi - ftype( i );
  ftype y0 = ( i   >= 0     ) ? data[i  ] : -data.front();
  ftype y1 = ( i+1 < size() ) ? data[i+1] : -data.back();
  return ( (1.0-xf)*y0 + xf*y1 );
}

const Histogram& Histogram::operator += ( const Histogram& h )
{
  if ( data.size() != h.data.size() )
    Message::message( Message_fatal("Histogram: sum of unequal histograms") );
  for ( int i = 0; i < data.size(); i++ ) data[i] += h.data[i];
  return (*this);
}

// generic ordinal

void Generic_ordinal::init( const Range<ftype>& range, const int num_ranges )
{
  nranges = ftype( num_ranges );
  hist.clear();
  hist.resize( num_ranges + 1, 0.0 );
  range_ = range;
}

void Generic_ordinal::init( const std::vector<ftype>& values, const int num_ranges )
{
  Range<ftype> range;
  for ( int i = 0; i < values.size(); i++ ) range.include( values[i] );
  init( range, num_ranges );
  for ( int i = 0; i < values.size(); i++ ) accumulate( values[i] );
  prep_ordinal();
}

ftype Generic_ordinal::ordinal( const ftype& value ) const
{
  ftype r = ( value - range_.min() ) / ( range_.max() - range_.min() );
  r = Util::bound( 0.0, r, 0.99999 ) * nranges;
  int i = int( r );
  r -= floor( r );
  return (1.0-r)*hist[i] + r*hist[i+1];
}

void Generic_ordinal::accumulate( const ftype& value )
{
  ftype r = ( value - range_.min() ) / ( range_.max() - range_.min() );
  r = Util::bound( 0.0, r, 0.99999 ) * nranges;
  int i = int( r );
  hist[i+1] += 1.0;
}

void Generic_ordinal::accumulate( const ftype& value, const ftype& weight )
{
  ftype r = ( value - range_.min() ) / ( range_.max() - range_.min() );
  r = Util::bound( 0.0, r, 0.99999 ) * nranges;
  int i = int( r );
  hist[i+1] += weight;
}

void Generic_ordinal::prep_ordinal()
{
  for ( int i = 1; i < hist.size(); i++ )
    hist[i] += hist[i-1];
  for ( int i = 0; i < hist.size(); i++ )
    hist[i] = hist[i] / hist.back();
}

void Generic_ordinal::invert()
{
  std::vector<ftype> hinv( hist.size() );
  hinv[0]             = range_.min();
  hinv[hist.size()-1] = range_.max();
  for ( int i = 1; i < hist.size()-1; i++ ) {
    ftype target = ftype(i)/nranges;
    ftype guess = 0.5 * ( range_.max() + range_.min() );
    ftype step  = 0.5 * ( range_.max() - range_.min() );
    for ( int j = 0; j < 10; j++ ) {
      if ( ordinal( guess ) > target ) guess -= step;
      else                             guess += step;
      step *= 0.5;
    }
    hinv[i] = Util::bound( range_.min(), guess, range_.max() );
  }
  range_ = Range<ftype>( 0.0, 1.0 );
  hist = hinv;
}


// deprecated fns

void Generic_ordinal::init( const int num_ranges )
{
  nranges = ftype( num_ranges );
  hist.clear();
  hist.resize( num_ranges + 1, 0.0 );
}

void Generic_ordinal::add_pass_1( const ftype& value )
{
  range_.include( value );
}

void Generic_ordinal::add_pass_2( const ftype& value )
{
  accumulate( value );
}


} // namespace clipper
