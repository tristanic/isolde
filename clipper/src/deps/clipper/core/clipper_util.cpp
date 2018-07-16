/* clipper_util.cpp: implementation file for clipper helper functions */
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


#include "clipper_util.h"


namespace clipper {

float  Util::nanf_;  //!< declare static
double Util::nand_;  //!< declare static
ftype  Util::nan_ ;  //!< declare static
ftype Util::onepi_ = M_PI;        //!< one*pi
ftype Util::twopi_ = 2.0*M_PI;    //!< two*pi
ftype Util::twopi2_ = 2.0*M_PI*M_PI; //!< two*pi*pi
ftype Util::eightpi2_ = 8.0*M_PI*M_PI; //!< two*pi*pi
ftype Util::d2rad_ = M_PI/180.0;  //!< degree-radian conversion
ftype Util::sim_a = 1.639294;  //!< sim fn param
ftype Util::sim_b = 3.553967;  //!< sim fn param
ftype Util::sim_c = 2.228716;  //!< sim fn param
ftype Util::sim_d = 3.524142;  //!< sim fn param
ftype Util::sim_e = 7.107935;  //!< sim fn param
ftype Util::sim_A = -1.28173889903;  //!< invsim fn param
ftype Util::sim_B =  0.69231689903;  //!< invsim fn param
ftype Util::sim_C = -1.33099462667;  //!< invsim fn param
ftype Util::sim_g =  2.13643992379;  //!< invsim fn param
ftype Util::sim_p =  0.04613803811;  //!< invsim fn param
ftype Util::sim_q =  1.82167089029;  //!< invsim fn param
ftype Util::sim_r = -0.74817947490;  //!< invsim fn param


Util::Util()
{
  ((U32*)&nanf_)->i = CLIPPER_NULL_MASK_32;
  ((U64*)&nand_)->i = CLIPPER_NULL_MASK_64;
  if        ( sizeof(ftype) == 4 ) {
    ((U32*)&nan_)->i = CLIPPER_NULL_MASK_32;
  } else if ( sizeof(ftype) == 8 ) {
    ((U64*)&nan_)->i = CLIPPER_NULL_MASK_64;
  } else {
    /* fail on build */;
  }
}

/*! \param x The argument \return I1(x)/I0(x) */
ftype Util::sim( const ftype& x )
{
  if (x >= 0.0) return (((x + sim_a)*x + sim_b)*x)
		  / (((x + sim_c)*x + sim_d)*x + sim_e);
  else          return -(-(-(-x + sim_a)*x + sim_b)*x)
		  / (-(-(-x + sim_c)*x + sim_d)*x + sim_e);
}

/*! \param x I1(y)/I0(y) \return y */
ftype Util::invsim( const ftype& x )
{
  ftype x0 = fabs(x);
  ftype a0 = -7.107935*x0;
  ftype a1 = 3.553967-3.524142*x0;
  ftype a2 = 1.639294-2.228716*x0;
  ftype a3 = 1.0-x0;
  ftype w = a2/(3.0*a3);
  ftype p = a1/(3.0*a3)-w*w;
  ftype q = -w*w*w+0.5*(a1*w-a0)/a3;
  ftype d = sqrt(q*q+p*p*p);
  ftype q1 = q + d;
  ftype q2 = q - d;
  ftype r1 = pow(fabs(q1), 1.0/3.0);
  ftype r2 = pow(fabs(q2), 1.0/3.0);
  if (x >= 0.0) return  (((q1>0.0)? r1 : -r1) + ((q2>0.0)? r2 : -r2) - w);
  else          return -(((q1>0.0)? r1 : -r1) + ((q2>0.0)? r2 : -r2) - w);
}

ftype Util::sim_integ( const ftype& x0 )
{
  const ftype x = fabs(x0);
  const ftype z = (x+sim_p)/sim_q;
  return sim_A*log(x+sim_g) + 0.5*sim_B*log(z*z+1.0) + sim_r*atan(z) + x + 1.0;
}

ftype Util::sim_deriv( const ftype& x )
{
  if (x >= 0.0) return (((((sim_c-sim_a)*x+(2.0*sim_d-2.0*sim_b))*x+(3.0*sim_e+sim_a*sim_d-sim_b*sim_c))*x+(2.0*sim_a*sim_e))*x+(sim_b*sim_e)) / pow( (((x + sim_c)*x + sim_d)*x + sim_e), 2.0 );
  else          return (-(-(-(-(sim_c-sim_a)*x+(2.0*sim_d-2.0*sim_b))*x+(3.0*sim_e+sim_a*sim_d-sim_b*sim_c))*x+(2.0*sim_a*sim_e))*x+(sim_b*sim_e)) / pow( (-(-(-x + sim_c)*x + sim_d)*x + sim_e), 2.0 );
}

ftype Util::sim_deriv_recur( const ftype& x0 )
{
  const ftype x = fabs(x0);
  const ftype m = sim(x);
  if ( x > 1.0e-4 )
    return ( -m/x + ( 1.0 - m*m ) );
  else
    return ( 0.5 - m*m );
}

ftype Util::bessel_i0( const ftype& x0 )
{
  ftype i0 = 0.0, t, x;
  x = fabs(x0);
  t=x/3.75;
  if (t < 1.0) {
    t=t*t;
    i0= ((((((t*0.0045813+0.0360768)*t+0.2659732)*t+
	    1.2067492)*t+3.0899424)*t+3.5156229)*t+1.0);
  } else {
    i0= (1.0/sqrt(x))*((((((((t*0.00392377-0.01647633)*t+
			     0.02635537)*t-0.02057706)*t+0.00916281)*t-
		 0.00157565)*t+0.00225319)*t+0.01328592)*t+0.39894228)*exp(x);
  }
  return i0;
}

/*! \param x Angle in degrees \return Angle in radians */
ftype Util::d2rad( const ftype& x )
{ return x*d2rad_; }

/*! \param x Angle in radians \return Angle in degrees */
ftype Util::rad2d( const ftype& x )
{ return x/d2rad_; }


// template instantiations
template ftype32 Util::max<ftype32>( const ftype32& a, const ftype32& b );
template ftype32 Util::min<ftype32>( const ftype32& a, const ftype32& b );
template ftype32 Util::bound<ftype32>( const ftype32& a, const ftype32& b, const ftype32& c );
template void Util::swap<ftype32>( ftype32& a, ftype32& b );
template void Util::swap<ftype32>( ftype32& a, ftype32& b, ftype32& c );
template ftype32 Util::sqr<ftype32>( const ftype32& a );
template ftype32 Util::isqrt<ftype32>( const ftype32& a );

template ftype64 Util::max<ftype64>( const ftype64& a, const ftype64& b );
template ftype64 Util::min<ftype64>( const ftype64& a, const ftype64& b );
template ftype64 Util::bound<ftype64>( const ftype64& a, const ftype64& b, const ftype64& c );
template void Util::swap<ftype64>( ftype64& a, ftype64& b );
template void Util::swap<ftype64>( ftype64& a, ftype64& b, ftype64& c );
template ftype64 Util::sqr<ftype64>( const ftype64& a );
template ftype64 Util::isqrt<ftype64>( const ftype64& a );

template int Util::max<int>( const int& a, const int& b );
template int Util::min<int>( const int& a, const int& b );
template int Util::bound<int>( const int& a, const int& b, const int& c );
template void Util::swap<int>( int& a, int& b );
template void Util::swap<int>( int& a, int& b, int& c );
template int Util::sqr<int>( const int& a );
template int Util::isqrt<int>( const int& a );


} // namespace clipper
