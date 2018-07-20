/*! \file lib/map_interp.h
    Generic interpolation methods for crystal and non-crystal maps.
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


#ifndef CLIPPER_MAP_INTERP
#define CLIPPER_MAP_INTERP

#include "derivs.h"
#include "../imex.h"

namespace clipper
{

  //! Wrapper class for zeroth-order (nearest neighbour) interpolation fns
  /*! These can be used through the built-in methods in Xmap/NXmap, or
    passed to methods to allow a choice of interpolation methods, or
    directly by providing the map as an argument. For example:
    \code
    NXmap<float> nxmap;
    Coord_map c;
    ...
    x1 = Interp_nearest<float>::interp( nxmap, c );
    x2 = nxmap.interp<Interp_nearest>( c );
    \endcode
  */
  class CLIPPER_IMEX Interp_nearest
  {
  public:
    template<class M> static bool can_interp( const M& map, const Coord_map& pos );  //!< Test if we can interpolate in map M at coord
    template<class T, class M> static void interp( const M& map, const Coord_map& pos, T& val );  //!< Interpolate map M using type T at coord
    inline static int order() { return 0; }  //!< Order of interpolant
  };

  //! Wrapper class for first-order (linear) interpolation fns
  /*! These can be used through the built-in methods in Xmap/NXmap, or
    passed to methods to allow a choice of interpolation methods, or
    directly by providing the map as an argument. For example:
    \code
    NXmap<float> nxmap;
    Coord_map c;
    float x1, x2;
    ...
    Interp_linear<float>::interp( nxmap, c, x1 );
    x2 = nxmap.interp<Interp_linear>( c );
    \endcode
  */
  class CLIPPER_IMEX Interp_linear
  {
  public:
    template<class M> static bool can_interp( const M& map, const Coord_map& pos );  //!< Test if we can interpolate in map M at coord
    template<class T, class M> static void interp( const M& map, const Coord_map& pos, T& val );  //!< Interpolate map M using type T at coord
    inline static int order() { return 1; }  //!< Order of interpolant
  };

  //! Wrapper class for third-order (cubic) interpolation fns
  /*! These can be used through the built-in methods in Xmap/NXmap, or
    passed to methods to allow a choice of interpolation methods, or
    directly by providing the map as an argument. For example:
    \code
    NXmap<float> nxmap;
    Coord_map c;
    float x1, x2;
    ...
    Interp_cubic::interp( nxmap, c, x1 );
    x2 = nxmap.interp<Interp_cubic>( c );
    \endcode
  */
  class CLIPPER_IMEX Interp_cubic
  {
  public:
    template<class M> static bool can_interp( const M& map, const Coord_map& pos );  //!< Test if we can interpolate in map M at coord
    template<class T, class M> static void interp( const M& map, const Coord_map& pos, T& val );  //!< Interpolate map M using type T at coord
    template<class T, class M> static void interp_grad( const M& map, const Coord_map& pos, T& val, Grad_map<T>& grad );
    template<class T, class M> static void interp_curv( const M& map, const Coord_map& pos, T& val, Grad_map<T>& grad, Curv_map<T>& curv );
    inline static int order() { return 3; }  //!< Order of interpolant
  };



  // template implementations

  /*! The map is queried to see if interpolation is possible at the
    given coord. For a crystallographic map, this is always true. For
    a non-crystallographic map, this depends if the point and enough
    neighbours are in the grid.
    \param map The map on which to perform the calculation.
    \param pos The map coord at which the density is to be calcuated. */
  template<class M> bool Interp_nearest::can_interp( const M& map, const Coord_map& pos )
    { return map.in_map( pos.coord_grid() ); }

  /*! The value of the map at the supplied map coordinate is
    calculated by zeroth order (nearest neighbour) interpolation based
    on 1 point.
    \param map The map on which to perform the calculation.
    \param pos The map coord at which the density is to be calcuated.
    \return The value of the density at that point.
   */
  template<class T, class M> void Interp_nearest::interp( const M& map, const Coord_map& pos, T& val )
    { val = map.get_data( pos.coord_grid() ); }


  /*! The map is queried to see if interpolation is possible at the
    given coord. For a crystallographic map, this is always true. For
    a non-crystallographic map, this depends if the point and enough
    neighbours are in the grid.
    \param map The map on which to perform the calculation.
    \param pos The map coord at which the density is to be calcuated. */
  template<class M> bool Interp_linear::can_interp( const M& map, const Coord_map& pos )
  {
    Coord_grid c( pos.floor() );  // for even order change floor to coord_grid
    c.u() -= order()/2; c.v() -= order()/2; c.w() -= order()/2;
    if ( map.in_map( c ) ) {
      c.u() += order(); c.v() += order(); c.w() += order();
      return map.in_map( c );
    }
    return false;
  }

  /*! The value of the map at the supplied map coordinate is
    calculated by first order (linear) interpolation based on 8
    neighbouring points.
    \param map The map on which to perform the calculation.
    \param pos The map coord at which the density is to be calcuated.
    \return The value of the density at that point.
   */
  template<class T, class M> void Interp_linear::interp( const M& map, const Coord_map& pos, T& val )
  {
    ftype u0 = floor( pos.u() );
    ftype v0 = floor( pos.v() );
    ftype w0 = floor( pos.w() );
    typename M::Map_reference_coord
      r( map, Coord_grid( int(u0), int(v0), int(w0) ) );
    T cu1( pos.u() - u0 );
    T cv1( pos.v() - v0 );
    T cw1( pos.w() - w0 );
    T cu0( 1.0 - cu1 );
    T cv0( 1.0 - cv1 );
    T cw0( 1.0 - cw1 );
    T r00 = cw0 * map[ r ];  // careful with evaluation order
    r00  += cw1 * map[ r.next_w() ];
    T r01 = cw1 * map[ r.next_v() ];
    r01  += cw0 * map[ r.prev_w() ];
    T r11 = cw0 * map[ r.next_u() ];
    r11  += cw1 * map[ r.next_w() ];
    T r10 = cw1 * map[ r.prev_v() ];
    r10  += cw0 * map[ r.prev_w() ];
    val = ( cu0*( cv0*r00 + cv1*r01 ) + cu1*( cv0*r10 + cv1*r11 ) );
  }


  /*! The map is queried to see if interpolation is possible at the
    given coord. For a crystallographic map, this is always true. For
    a non-crystallographic map, this depends if the point and enough
    neighbours are in the grid.
    \param map The map on which to perform the calculation.
    \param pos The map coord at which the density is to be calcuated. */
  template<class M> bool Interp_cubic::can_interp( const M& map, const Coord_map& pos )
  {
    Coord_grid c( pos.floor() );  // for even order change floor to coord_grid
    c.u() -= order()/2; c.v() -= order()/2; c.w() -= order()/2;
    if ( map.in_map( c ) ) {
      c.u() += order(); c.v() += order(); c.w() += order();
      return map.in_map( c );
    }
    return false;
  }

  /*! The value of the map at the supplied map coordinate is
    calculated by third order (cubic) interpolation based on the
    surrounding 64 points.
    \param pos The fractional coord at which the density is to be calcuated.
    \return The value of the density at that point. */
  template<class T, class M> void Interp_cubic::interp( const M& map, const Coord_map& pos, T& val )
  {
    ftype u0 = floor( pos.u() );
    ftype v0 = floor( pos.v() );
    ftype w0 = floor( pos.w() );
    typename M::Map_reference_coord iw, iv,
      iu( map, Coord_grid( int(u0)-1, int(v0)-1, int(w0)-1 ) );
    T su, sv, sw, cu[4], cv[4], cw[4];
    T cu1( pos.u() - u0 );
    T cv1( pos.v() - v0 );
    T cw1( pos.w() - w0 );
    T cu0( 1.0 - cu1 );
    T cv0( 1.0 - cv1 );
    T cw0( 1.0 - cw1 );
    cu[0] = -0.5*cu1*cu0*cu0; // cubic spline coeffs: u
    cu[1] = cu0*( -1.5*cu1*cu1 + cu1 + 1.0 );
    cu[2] = cu1*( -1.5*cu0*cu0 + cu0 + 1.0 );
    cu[3] = -0.5*cu1*cu1*cu0;
    cv[0] = -0.5*cv1*cv0*cv0; // cubic spline coeffs: v
    cv[1] = cv0*( -1.5*cv1*cv1 + cv1 + 1.0 );
    cv[2] = cv1*( -1.5*cv0*cv0 + cv0 + 1.0 );
    cv[3] = -0.5*cv1*cv1*cv0;
    cw[0] = -0.5*cw1*cw0*cw0; // cubic spline coeffs: w
    cw[1] = cw0*( -1.5*cw1*cw1 + cw1 + 1.0 );
    cw[2] = cw1*( -1.5*cw0*cw0 + cw0 + 1.0 );
    cw[3] = -0.5*cw1*cw1*cw0;
    su = 0.0;
    int i, j;
    for ( j = 0; j < 4; j++ ) {
      iv = iu;
      sv = 0.0;
      for ( i = 0; i < 4; i++ ) {
	iw = iv;
	  sw  = cw[0] * T( map[ iw ] );
	  sw += cw[1] * T( map[ iw.next_w() ] );
	  sw += cw[2] * T( map[ iw.next_w() ] );
	  sw += cw[3] * T( map[ iw.next_w() ] );
	sv += cv[i] * sw;
	iv.next_v();
      }
      su += cu[j] * sv;
      iu.next_u();
    }
    val = su;
  }


  /*! The value of the map at the supplied map coordinate and its
    gradient are calculated by third order (cubic) interpolation based
    on the surrounding 64 points.
    \param pos The fractional coord at which the density is to be calcuated.
    \param val The value of the density at that point.
    \param grad The interpolated value as a gradient vector with respect
    to the fractional coordinates (see Cell::coord_orth). */
  template<class T, class M> void Interp_cubic::interp_grad( const M& map, const Coord_map& pos, T& val, Grad_map<T>& grad )
  {
    ftype u0 = floor( pos.u() );
    ftype v0 = floor( pos.v() );
    ftype w0 = floor( pos.w() );
    typename M::Map_reference_coord iw, iv,
      iu( map, Coord_grid( int(u0)-1, int(v0)-1, int(w0)-1 ) );
    T s1, s2, s3, du1, dv1, dv2, dw1, dw2, dw3;
    T cu[4], cv[4], cw[4], gu[4], gv[4], gw[4];
    T cu1( pos.u() - u0 );
    T cv1( pos.v() - v0 );
    T cw1( pos.w() - w0 );
    T cu0( 1.0 - cu1 );
    T cv0( 1.0 - cv1 );
    T cw0( 1.0 - cw1 );
    cu[0] = -0.5*cu1*cu0*cu0; // cubic spline coeffs: u
    cu[1] = cu0*( -1.5*cu1*cu1 + cu1 + 1.0 );
    cu[2] = cu1*( -1.5*cu0*cu0 + cu0 + 1.0 );
    cu[3] = -0.5*cu1*cu1*cu0;
    cv[0] = -0.5*cv1*cv0*cv0; // cubic spline coeffs: v
    cv[1] = cv0*( -1.5*cv1*cv1 + cv1 + 1.0 );
    cv[2] = cv1*( -1.5*cv0*cv0 + cv0 + 1.0 );
    cv[3] = -0.5*cv1*cv1*cv0;
    cw[0] = -0.5*cw1*cw0*cw0; // cubic spline coeffs: w
    cw[1] = cw0*( -1.5*cw1*cw1 + cw1 + 1.0 );
    cw[2] = cw1*( -1.5*cw0*cw0 + cw0 + 1.0 );
    cw[3] = -0.5*cw1*cw1*cw0;
    gu[0] =  cu0*( 1.5*cu1 - 0.5 ); // cubic spline grad coeffs: u
    gu[1] =  cu1*( 4.5*cu1 - 5.0 );
    gu[2] = -cu0*( 4.5*cu0 - 5.0 );
    gu[3] = -cu1*( 1.5*cu0 - 0.5 );
    gv[0] =  cv0*( 1.5*cv1 - 0.5 ); // cubic spline grad coeffs: v
    gv[1] =  cv1*( 4.5*cv1 - 5.0 );
    gv[2] = -cv0*( 4.5*cv0 - 5.0 );
    gv[3] = -cv1*( 1.5*cv0 - 0.5 );
    gw[0] =  cw0*( 1.5*cw1 - 0.5 ); // cubic spline grad coeffs: w
    gw[1] =  cw1*( 4.5*cw1 - 5.0 );
    gw[2] = -cw0*( 4.5*cw0 - 5.0 );
    gw[3] = -cw1*( 1.5*cw0 - 0.5 );
    s1 = du1 = dv1 = dw1 = 0.0;
    int i, j;
    for ( j = 0; j < 4; j++ ) {
      iv = iu;
      s2 = dv2 = dw2 = 0.0;
      for ( i = 0; i < 4; i++ ) {
	iw = iv;
	  s3  = cw[0] * T( map[ iw ] );
	  dw3  = gw[0] * T( map[ iw ] );
	  iw.next_w();
	  s3 += cw[1] * T( map[ iw ] );
	  dw3 += gw[1] * T( map[ iw ] );
	  iw.next_w();
	  s3 += cw[2] * T( map[ iw ] );
	  dw3 += gw[2] * T( map[ iw ] );
	  iw.next_w();
	  s3 += cw[3] * T( map[ iw ] );
	  dw3 += gw[3] * T( map[ iw ] );
	s2 += cv[i] * s3;
	dv2 += gv[i] * s3;
        dw2 += cv[i] * dw3;
	iv.next_v();
      }
      s1 += cu[j] * s2;
      du1 += gu[j] * s2;
      dv1 += cu[j] * dv2;
      dw1 += cu[j] * dw2;
      iu.next_u();
    }
    val = s1;
    grad = Grad_map<T>( du1, dv1, dw1 );
  }


  /*! The value of the map at the supplied map coordinate and its
    gradient are calculated by third order (cubic) interpolation based
    on the surrounding 64 points.
    \param pos The fractional coord at which the density is to be calcuated.
    \param val The value of the density at that point.
    \param grad The interpolated value as a gradient vector with respect
    to the fractional coordinates (see Cell::coord_orth). */
  template<class T, class M> void Interp_cubic::interp_curv( const M& map, const Coord_map& pos, T& val, Grad_map<T>& grad, Curv_map<T>& curv )
  {
    ftype u0 = floor( pos.u() );
    ftype v0 = floor( pos.v() );
    ftype w0 = floor( pos.w() );
    typename M::Map_reference_coord iw, iv,
      iu( map, Coord_grid( int(u0)-1, int(v0)-1, int(w0)-1 ) );
    T s1, s2, s3, du1, dv1, dv2, dw1, dw2, dw3;
    T duv1, duw1, dvw1, dvw2, duu1, dvv1, dvv2, dww1, dww2, dww3;
    T cu[4], cv[4], cw[4], gu[4], gv[4], gw[4], ggu[4], ggv[4], ggw[4];
    T cu1( pos.u() - u0 );
    T cv1( pos.v() - v0 );
    T cw1( pos.w() - w0 );
    T cu0( 1.0 - cu1 );
    T cv0( 1.0 - cv1 );
    T cw0( 1.0 - cw1 );
    cu[0] = -0.5*cu1*cu0*cu0; // cubic spline coeffs: u
    cu[1] = cu0*( -1.5*cu1*cu1 + cu1 + 1.0 );
    cu[2] = cu1*( -1.5*cu0*cu0 + cu0 + 1.0 );
    cu[3] = -0.5*cu1*cu1*cu0;
    cv[0] = -0.5*cv1*cv0*cv0; // cubic spline coeffs: v
    cv[1] = cv0*( -1.5*cv1*cv1 + cv1 + 1.0 );
    cv[2] = cv1*( -1.5*cv0*cv0 + cv0 + 1.0 );
    cv[3] = -0.5*cv1*cv1*cv0;
    cw[0] = -0.5*cw1*cw0*cw0; // cubic spline coeffs: w
    cw[1] = cw0*( -1.5*cw1*cw1 + cw1 + 1.0 );
    cw[2] = cw1*( -1.5*cw0*cw0 + cw0 + 1.0 );
    cw[3] = -0.5*cw1*cw1*cw0;
    gu[0] =  cu0*( 1.5*cu1 - 0.5 ); // cubic spline grad coeffs: u
    gu[1] =  cu1*( 4.5*cu1 - 5.0 );
    gu[2] = -cu0*( 4.5*cu0 - 5.0 );
    gu[3] = -cu1*( 1.5*cu0 - 0.5 );
    gv[0] =  cv0*( 1.5*cv1 - 0.5 ); // cubic spline grad coeffs: v
    gv[1] =  cv1*( 4.5*cv1 - 5.0 );
    gv[2] = -cv0*( 4.5*cv0 - 5.0 );
    gv[3] = -cv1*( 1.5*cv0 - 0.5 );
    gw[0] =  cw0*( 1.5*cw1 - 0.5 ); // cubic spline grad coeffs: w
    gw[1] =  cw1*( 4.5*cw1 - 5.0 );
    gw[2] = -cw0*( 4.5*cw0 - 5.0 );
    gw[3] = -cw1*( 1.5*cw0 - 0.5 );
    ggu[0] =  2.0 - 3.0*cu1; // cubic spline curv coeffs: u
    ggu[1] =  9.0*cu1 - 5.0;
    ggu[2] =  9.0*cu0 - 5.0;
    ggu[3] =  2.0 - 3.0*cu0;
    ggv[0] =  2.0 - 3.0*cv1; // cubic spline curv coeffs: v
    ggv[1] =  9.0*cv1 - 5.0;
    ggv[2] =  9.0*cv0 - 5.0;
    ggv[3] =  2.0 - 3.0*cv0;
    ggw[0] =  2.0 - 3.0*cw1; // cubic spline curv coeffs: w
    ggw[1] =  9.0*cw1 - 5.0;
    ggw[2] =  9.0*cw0 - 5.0;
    ggw[3] =  2.0 - 3.0*cw0;
    s1 = du1 = dv1 = dw1 = duv1 = duw1 = dvw1 = duu1 = dvv1 = dww1 = 0.0;
    int i, j;
    for ( j = 0; j < 4; j++ ) {
      iv = iu;
      s2 = dv2 = dw2 = dvw2 = dvv2 = dww2 = 0.0;
      for ( i = 0; i < 4; i++ ) {
	iw = iv;
	  s3  = cw[0] * T( map[ iw ] );
	  dw3  = gw[0] * T( map[ iw ] );
	  dww3 = ggw[0] * T( map[ iw ] );
	  iw.next_w();
	  s3 += cw[1] * T( map[ iw ] );
	  dw3 += gw[1] * T( map[ iw ] );
	  dww3 += ggw[1] * T( map[ iw ] );
	  iw.next_w();
	  s3 += cw[2] * T( map[ iw ] );
	  dw3 += gw[2] * T( map[ iw ] );
	  dww3 += ggw[2] * T( map[ iw ] );
	  iw.next_w();
	  s3 += cw[3] * T( map[ iw ] );
	  dw3 += gw[3] * T( map[ iw ] );
	  dww3 += ggw[3] * T( map[ iw ] );
	s2 += cv[i] * s3;
	dv2 += gv[i] * s3;
        dw2 += cv[i] * dw3;
	dvw2 += gv[i] * dw3;
	dvv2 += ggv[i] * s3;
        dww2 += cv[i] * dww3;
	iv.next_v();
      }
      s1 += cu[j] * s2;
      du1 += gu[j] * s2;
      dv1 += cu[j] * dv2;
      dw1 += cu[j] * dw2;
      duv1 += gu[j] * dv2;
      duw1 += gu[j] * dw2;
      dvw1 += cu[j] * dvw2;
      duu1 += ggu[j] * s2;
      dvv1 += cu[j] * dvv2;
      dww1 += cu[j] * dww2;
      iu.next_u();
    }
    val = s1;
    grad = Grad_map<T>( du1, dv1, dw1 );
    curv = Curv_map<T>( Mat33<T>( duu1, duv1, duw1,
				  duv1, dvv1, dvw1,
				  duw1, dvw1, dww1 ) );
  }


} // namespace clipper


#endif
