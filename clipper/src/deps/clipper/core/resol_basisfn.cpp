/* resol_basisfn.cpp: implementation file for resolution basis functions */
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
/* This code is derived from the 'dm' source code */

#include "resol_basisfn.h"

namespace clipper {


// resolution ordinal

void Resolution_ordinal::init( const HKL_info& hklinfo, const ftype& power )
{
  HKL_info::HKL_reference_index ih;
  Range<ftype> range;
  for ( ih = hklinfo.first(); !ih.last(); ih.next() )
    range.include( ih.invresolsq() );
  Generic_ordinal::init( range, 1000 );
  for ( ih = hklinfo.first(); !ih.last(); ih.next() )
    accumulate( ih.invresolsq() );
  prep_ordinal();

  for ( int i = 0; i < hist.size(); i++ )
    hist[i] = pow( hist[i], 1.0/power );
}

void Resolution_ordinal::init( const HKL_data_base& hkldata, const ftype& power )
{
  HKL_info::HKL_reference_index ih;
  Range<ftype> range;
  for ( ih = hkldata.first_data(); !ih.last(); hkldata.next_data(ih) )
    range.include( ih.invresolsq() );
  Generic_ordinal::init( range, 1000 );
  for ( ih = hkldata.first_data(); !ih.last(); hkldata.next_data(ih) )
    accumulate( ih.invresolsq() );
  prep_ordinal();

  for ( int i = 0; i < hist.size(); i++ )
    hist[i] = pow( hist[i], 1.0/power );
}

void Resolution_ordinal::init( const HKL_data_base& hkldata, const Cell& cell, const ftype& power )
{
  HKL_info::HKL_reference_index ih;
  Range<ftype> range;
  for ( ih = hkldata.first_data(); !ih.last(); hkldata.next_data(ih) )
    range.include( ih.hkl().invresolsq( cell ) );
  Generic_ordinal::init( range, 1000 );
  for ( ih = hkldata.first_data(); !ih.last(); hkldata.next_data(ih) )
    accumulate( ih.hkl().invresolsq( cell ) );
  prep_ordinal();

  for ( int i = 0; i < hist.size(); i++ )
    hist[i] = pow( hist[i], 1.0/power );
}


// binner

ftype BasisFn_binner::f_s( const ftype& s, const std::vector<ftype>& params ) const
{
  const int& nbins = num_params();
  // convert reso to int
  const int bin = Util::bound( 0, Util::intf( ftype(nbins) * s_ord.ordinal( s ) ), nbins-1 );
  return params[bin];
}

const BasisFn_base::Fderiv& BasisFn_binner::fderiv_s( const ftype& s, const std::vector<ftype>& params ) const
{
  const int& nbins = num_params();
  for ( int i = 0; i < nbins; i++ ) result().df[i] = 0.0;

  // convert reso to int
  const int bin = Util::bound( 0, Util::intf( ftype(nbins) * s_ord.ordinal( s ) ), nbins-1 );
  // make vector of derivatives
  result().f = params[bin];
  result().df[bin] = 1.0;
  return result();
}


// linear binner

ftype BasisFn_linear::f_s( const ftype& s_, const std::vector<ftype>& params ) const
{
  const int& nbins = num_params();
  const ftype s = ftype(nbins) * s_ord.ordinal( s_ );
  const int i = Util::intf( s );
  const ftype ds = s - ftype(i);
  const int i0 = Util::bound( 0, i  , nbins-1 );
  const int i1 = Util::bound( 0, i+1, nbins-1 );
  return ( params[i0]*(1.0-ds) + params[i1]*(ds) );
}

const BasisFn_base::Fderiv& BasisFn_linear::fderiv_s( const ftype& s_, const std::vector<ftype>& params ) const
{
  const int& nbins = num_params();
  for ( int i = 0; i < nbins; i++ ) result().df[i] = 0.0;

  const ftype s = ftype(nbins) * s_ord.ordinal( s_ );
  const int i = Util::intf( s );
  const ftype ds = s - ftype(i);
  const int i0 = Util::bound( 0, i  , nbins-1 );
  const int i1 = Util::bound( 0, i+1, nbins-1 );
  result().f = params[i0]*(1.0-ds) + params[i1]*(ds);
  result().df[i0] += (1.0-ds);
  result().df[i1] += (ds);
  return result();
}


// smooth binner

ftype BasisFn_spline::f_s( const ftype& s_, const std::vector<ftype>& params ) const
{
  const int& nbins = num_params();
  const ftype s = ftype(nbins) * s_ord.ordinal( s_ );
  const int i = Util::intf( s );
  const ftype ds = s - ftype(i) - 0.5;
  const int i0 = Util::bound( 0, i-1, nbins-1 );
  const int i1 = Util::bound( 0, i  , nbins-1 );
  const int i2 = Util::bound( 0, i+1, nbins-1 );
  return ( params[i0]*0.5*(ds-0.5)*(ds-0.5) + params[i1]*(0.75-ds*ds) + params[i2]*0.5*(ds+0.5)*(ds+0.5) );
}

const BasisFn_base::Fderiv& BasisFn_spline::fderiv_s( const ftype& s_, const std::vector<ftype>& params ) const
{
  const int& nbins = num_params();
  for ( int i = 0; i < nbins; i++ ) result().df[i] = 0.0;

  const ftype s = ftype(nbins) * s_ord.ordinal( s_ );
  const int i = Util::intf( s );
  const ftype ds = s - ftype(i) - 0.5;
  const int i0 = Util::bound( 0, i-1, nbins-1 );
  const int i1 = Util::bound( 0, i  , nbins-1 );
  const int i2 = Util::bound( 0, i+1, nbins-1 );
  result().f = params[i0]*0.5*(ds-0.5)*(ds-0.5) + params[i1]*(0.75-ds*ds) + params[i2]*0.5*(ds+0.5)*(ds+0.5);
  result().df[i0] += 0.5*(ds-0.5)*(ds-0.5);
  result().df[i1] += 0.75-ds*ds;
  result().df[i2] += 0.5*(ds+0.5)*(ds+0.5);
  return result();
}


// gaussian
/*
ftype BasisFn_gaussian::f_s( const ftype& s, const std::vector<ftype>& params ) const
{
  // generate Gaussian
  return exp( - params[1] * s + params[0] );
}
*/

const BasisFn_base::Fderiv& BasisFn_gaussian::fderiv_s( const ftype& s, const std::vector<ftype>& params ) const
{
  ftype f = exp( - params[1] * s + params[0] );
  result().f     = result().df[0]    = result().df2(0,0) = f;
  result().df[1] = result().df2(0,1) = result().df2(1,0) = -s * f;
  result().df2(1,1) = s * s * f;
  return result();
}

ftype BasisFn_gaussian::scale( const std::vector<ftype>& params ) const
{
  return exp( params[0] );
}

ftype BasisFn_gaussian::u_iso( const std::vector<ftype>& params ) const
{
  return params[1] / Util::twopi2();
}

// aniso gaussian
/*
ftype BasisFn_aniso_gaussian::f_coord( const Coord_reci_orth& xs, const std::vector<ftype>& params ) const
{
  // generate Gaussian
  return exp( params[0] - ( xs[0]*xs[0]*params[1] + xs[1]*xs[1]*params[2] + xs[2]*xs[2]*params[3] + 2.0*(xs[0]*xs[1]*params[4] + xs[0]*xs[2]*params[5] + xs[1]*xs[2]*params[6]) ) );
}
*/

const BasisFn_base::Fderiv& BasisFn_aniso_gaussian::fderiv_coord( const Coord_reci_orth& xs, const std::vector<ftype>& params ) const
{
  ftype c[7];
  c[0] = 1.0;
  c[1] = -xs[0]*xs[0];
  c[2] = -xs[1]*xs[1];
  c[3] = -xs[2]*xs[2];
  c[4] = -2.0*xs[0]*xs[1];
  c[5] = -2.0*xs[0]*xs[2];
  c[6] = -2.0*xs[1]*xs[2];
  ftype f = exp( params[0] + c[1]*params[1] + c[2]*params[2] + c[3]*params[3] +
		             c[4]*params[4] + c[5]*params[5] + c[6]*params[6] );
  result().f = f;
  int i, j;
  for ( j = 0; j < 7; j++ ) {
    result().df[j] = c[j]*f;
    for ( i = 0; i < 7; i++ )
      result().df2(i,j) = c[i]*c[j]*f;
  }
  return result();
}


ftype BasisFn_aniso_gaussian::scale( const std::vector<ftype>& params ) const
{
  return exp( params[0] );
}

U_aniso_orth BasisFn_aniso_gaussian::u_aniso_orth( const std::vector<ftype>& params ) const
{
  return U_aniso_orth( params[1]/Util::twopi2(), params[2]/Util::twopi2(),
		       params[3]/Util::twopi2(), params[4]/Util::twopi2(),
		       params[5]/Util::twopi2(), params[6]/Util::twopi2() );
}

// log_gaussian
/*
ftype BasisFn_log_gaussian::f_s( const ftype& s, const std::vector<ftype>& params ) const
{
  // generate Gaussian
  return exp( - params[1] * s + params[0] );
}
*/

const BasisFn_base::Fderiv& BasisFn_log_gaussian::fderiv_s( const ftype& s, const std::vector<ftype>& params ) const
{
  ftype f = - params[1] * s + params[0];
  result().f     = f;
  result().df[0] = 1.0;
  result().df[1] = -s;
  return result();
}

ftype BasisFn_log_gaussian::scale( const std::vector<ftype>& params ) const
{
  return exp( params[0] );
}

ftype BasisFn_log_gaussian::u_iso( const std::vector<ftype>& params ) const
{
  return params[1] / Util::twopi2();
}

// log_aniso gaussian
/*
ftype BasisFn_log_aniso_gaussian::f_coord( const Coord_reci_orth& xs, const std::vector<ftype>& params ) const
{
  // generate Gaussian
  return exp( params[0] - ( xs[0]*xs[0]*params[1] + xs[1]*xs[1]*params[2] + xs[2]*xs[2]*params[3] + 2.0*(xs[0]*xs[1]*params[4] + xs[0]*xs[2]*params[5] + xs[1]*xs[2]*params[6]) ) );
}
*/

const BasisFn_base::Fderiv& BasisFn_log_aniso_gaussian::fderiv_coord( const Coord_reci_orth& xs, const std::vector<ftype>& params ) const
{
  result().df[0] = 1.0;
  result().df[1] = -xs[0]*xs[0];
  result().df[2] = -xs[1]*xs[1];
  result().df[3] = -xs[2]*xs[2];
  result().df[4] = -2.0*xs[0]*xs[1];
  result().df[5] = -2.0*xs[0]*xs[2];
  result().df[6] = -2.0*xs[1]*xs[2];
  result().f   = ( params[0] +
                   result().df[1]*params[1] +
                   result().df[2]*params[2] +
                   result().df[3]*params[3] +
                   result().df[4]*params[4] +
                   result().df[5]*params[5] +
                   result().df[6]*params[6] );
  return result();
}

ftype BasisFn_log_aniso_gaussian::scale( const std::vector<ftype>& params ) const
{
  return exp( params[0] );
}

U_aniso_orth BasisFn_log_aniso_gaussian::u_aniso_orth( const std::vector<ftype>& params ) const
{
  return U_aniso_orth( params[1]/Util::twopi2(), params[2]/Util::twopi2(),
                       params[3]/Util::twopi2(), params[4]/Util::twopi2(),
                       params[5]/Util::twopi2(), params[6]/Util::twopi2() );
}

// expcubic
/*
ftype BasisFn_expcubic::f_s( const ftype& s, const std::vector<ftype>& params ) const
{
  // generate Expcubic
  return exp( ( ( - params[3]*s + params[2] )*s - params[1] )*s + params[0] );
}
*/

const BasisFn_base::Fderiv& BasisFn_expcubic::fderiv_s( const ftype& s, const std::vector<ftype>& params ) const
{
  ftype f = exp( ( ( -params[3]*s + params[2] )*s - params[1] )*s + params[0] );
  result().f =
    result().df[0] =
    result().df2(0,0) = f;
  f *= -s;
  result().df[1] = 
    result().df2(0,1) =
    result().df2(1,0) = f;
  f *= -s;
  result().df[2] = 
    result().df2(0,2) =
    result().df2(1,1) =
    result().df2(2,0) = f;
  f *= -s;
  result().df[3] = 
    result().df2(0,3) =
    result().df2(1,2) =
    result().df2(2,1) =
    result().df2(3,0) = f;
  f *= -s;
  result().df2(1,3) =
    result().df2(2,2) =
    result().df2(3,1) = f;
  f *= -s;
  result().df2(2,3) =
    result().df2(3,2) = f;
  f *= -s;
  result().df2(3,3) = f;
  return result();
}


} // namespace clipper
