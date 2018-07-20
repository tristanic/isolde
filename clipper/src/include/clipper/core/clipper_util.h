/*! \file lib/clipper_util.h
    Header file for clipper helper functions
*/
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


#ifndef CLIPPER_UTIL
#define CLIPPER_UTIL


#include "clipper_precision.h"
#include "../imex.h"

namespace clipper
{

  //! Utility class
  /*! This class holds a set of useful static functions and values.
   You should never need to instantiate this class: Refer to members
   using Util:: */
  class CLIPPER_IMEX Util
  {
  private:
    typedef union { uitype32 i; ftype32 f; } U32;
    typedef union { uitype64 i; ftype64 f; } U64;
  public:
    Util();   //!< null constructor
    //! fast Util::nan() value
    static const ftype& nan() { return nan_; }
    //! fast Util::nan() value
    static const float& nanf() { return nanf_; }
    //! fast Util::nan() value
    static const double& nand() { return nand_; }
    //! set null generic value
    template<class T> inline static void set_null( T& t ) { t = T(-2147483648); }
    //! set null floating value - a specific value of NaN used for missings
    inline static void set_null( ftype32& f ) { U32* const u1=(U32* const)&f; const U32* const u2=(const U32* const)&nanf_; u1->i = u2->i; }
    //! set null floating value - a specific value of NaN used for missings
    inline static void set_null( ftype64& f ) { U64* const u1=(U64* const)&f; const U64* const u2=(const U64* const)&nand_; u1->i = u2->i; }
    //! fast test for null floating value - only works if set from Util::null()
    inline static bool is_null( const ftype32& f ) { U32 u1,u2; u1.f = f; u2.f = nanf_; return ( u1.i == u2.i ); }
    //! fast test for null floating value - only works if set from Util::null()
    inline static bool is_null( const ftype64& f ) { U64 u1,u2; u1.f = f; u2.f = nand_; return ( u1.i == u2.i ); }
    //! fast Util::nan() test
    /*! Used for missing entries: THIS DOES NOT DISTINGUISH BETWEEN NAN & INF */
    inline static bool is_nan( const ftype32 f ) { U32 u; u.f = f; return ((u.i&CLIPPER_NAN_MASK_A_32)==CLIPPER_NAN_MASK_A_32); }
    //! fast Util::nan() test
    /*! Used for missing entries: THIS DOES NOT DISTINGUISH BETWEEN NAN & INF */
    inline static bool is_nan( const ftype64 f ) { U64 u; u.f = f; return ((u.i&CLIPPER_NAN_MASK_A_64)==CLIPPER_NAN_MASK_A_64); }
    //! slow general NaN test for compatibility
    /*! Works for all architectures with IEEE arithmetic only */
    inline static bool isnan(const ftype32 f) { U32 u; u.f = f; return ((u.i&CLIPPER_NAN_MASK_A_32)==CLIPPER_NAN_MASK_A_32)&&((u.i&CLIPPER_NAN_MASK_B_32)!=0U); }
    //! slow general NaN test for compatibility
    /*! Works for all architectures with IEEE arithmetic only */
    inline static bool isnan(const ftype64 f) { U64 u; u.f = f; return ((u.i&CLIPPER_NAN_MASK_A_64)==CLIPPER_NAN_MASK_A_64)&&((u.i&CLIPPER_NAN_MASK_B_64)!=0U); }
    //! Sim function: I1(X)/I0(X)
    static ftype sim( const ftype& x );
    //! Inverse Sim function: I1(X)/I0(X)
    static ftype invsim( const ftype& x );
    //! Integral of Sim function: log(I0(X))
    static ftype sim_integ( const ftype& x );
    //! Derivative of Sim function: d/dx( I1(X)/I0(x) )
    static ftype sim_deriv( const ftype& x );
    //! Derivative of Sim function using recurrance: -sim(x)/x + (1 - sim(x)^2)
    static ftype sim_deriv_recur( const ftype& x );
    //! Arc hyperbolic tangent
    static ftype atanh( const ftype& x ) { return log((1.0+x)/(1.0-x))/2.0; }
    //! Modified Bessel function of the first kind
    static ftype bessel_i0( const ftype& x );
    //! Convert isotropic U-value to B-factor
    static ftype u2b( const ftype& x ) { return x * eightpi2_; }
    //! Convert isotropic B-factor to U-value
    static ftype b2u( const ftype& x ) { return x / eightpi2_; }
    //! Convert F+/F- to mean F, with NaN checks
    template<class T> inline static T mean( const T& pl, const T& mi )
      {
	if ( Util::is_nan((float)pl) ) return mi;
	else if (Util::is_nan((float)mi) ) return pl;
	else return 0.5*(pl+mi);
      }
    //! Convert sigF+/sigF-/cov to sig F, with NaN checks
    template<class T> inline static T sig_mean( const T& pl, const T& mi, const T& cov )
      {
	if ( Util::is_nan(pl) ) return mi;
	else if (Util::is_nan(mi) ) return pl;
	else if (Util::is_nan(cov) ) return 0.5*sqrt(pl*pl+mi*mi);
	else return 0.5*sqrt(pl*pl+mi*mi+2*cov);
      }

    //! Truncate-to-integer: int(floor(a))
    inline static int intf( const ftype& a ) { return int( floor( a ) ); }
    //! Truncate-to-integer above: int(ceil(a))
    inline static int intc( const ftype& a ) { return int( ceil( a ) ); }
    //! Round-to-integer: int(round(a))
    inline static int intr( const ftype& a ) { return int( rint( a ) ); }

    //! Corrected mod
    inline static ftype mod( const ftype& a, const ftype& b )
      { ftype c = fmod(a, b); if (c < 0) c+=b; return c;}
    //! Corrected mod
    inline static int mod( const int& a, const int& b )
      { int c = a%b; if (c < 0) c+=b; return c; }
    //! max
    template<class T> inline static T max(const T& a, const T& b)
      { return (a > b) ? a : b; }
    //! min
    template<class T> inline static T min(const T& a, const T& b)
      { return (a < b) ? a : b; }
    //! bound a value by limits
    template<class T> inline static T bound( const T& min, const T& val, const T& max ) { return ( (val < max) ? ( (val > min ) ? val : min ) : max ); }
    //! swap the contents of two objects
    template<class T> inline static void swap( T& a, T& b )
      { T c = a; a = b; b = c; }
    //! swap the contents of two objects, using third as store (for speed)
    template<class T> inline static void swap( T& a, T& b, T& c )
      { c = a; a = b; b = c; }
    //! square
    template<class T> inline static T sqr( const T& a ) { return a*a; }
    //! Integer square root (returns floor of sqrt)
    template<class T> inline static T isqrt( const T& n )
      { return T(floor(sqrt(ftype(n)))); }

    //! pi
    inline static const ftype& pi() { return onepi_; }
    //! 2 pi
    inline static const ftype& twopi() { return twopi_; }
    //! 2 pi squared
    inline static const ftype& twopi2() { return twopi2_; }
    //! 8 pi squared
    inline static const ftype& eightpi2() { return eightpi2_; }
    //! degree-to-radian conversion
    static ftype d2rad( const ftype& x );
    //! degree-to-radian conversion
    static ftype rad2d( const ftype& x );

  private:
    static float  nanf_;  //!< float NaN
    static double nand_;  //!< double NaN
    static ftype  nan_;   //!< ftype nan
    static ftype onepi_;  //!< one*pi
    static ftype twopi_;  //!< two*pi
    static ftype twopi2_; //!< two*pi*pi
    static ftype eightpi2_; //!< eight*pi*pi
    static ftype d2rad_;  //!< degree-radian conversion
    static ftype sim_a;   //!< sim fn param
    static ftype sim_b;   //!< sim fn param
    static ftype sim_c;   //!< sim fn param
    static ftype sim_d;   //!< sim fn param
    static ftype sim_e;   //!< sim fn param
    static ftype sim_A;   //!< invsim fn param
    static ftype sim_B;   //!< invsim fn param
    static ftype sim_C;   //!< invsim fn param
    static ftype sim_g;   //!< invsim fn param
    static ftype sim_p;   //!< invsim fn param
    static ftype sim_q;   //!< invsim fn param
    static ftype sim_r;   //!< invsim fn param
  };

} // namespace clipper

#endif
