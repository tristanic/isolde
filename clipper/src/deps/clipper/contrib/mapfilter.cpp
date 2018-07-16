/* mapfilter.cpp: Electron density calculation implementation */
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


#include "mapfilter.h"


namespace clipper {


/*! The filter function is constructed with the given filter.
  Scaling may also be optionally applied. This may be absolute, which
  just applies a scale factor to the result, or relative, which scales
  the filter relative to its own integral. Therefore, relative scaling
  with a scale factor of 1.0 gives an output map on the same scale as
  the input map.

  Note that the filter is not stored internally, and so must persist
  as long as the MapFilter is required.
  \param fltr The radial filter to apply.
  \param scale The scale factor to apply (default = 1.0).
  \param type The type of scaling to apply: NONE, Absolute, Relative. */
template<class T> MapFilter_slow<T>::MapFilter_slow( const MapFilterFn_base& fltr, const ftype scale, const TYPE type ) : fltr_( &fltr ), scale_( scale ), type_( type )
{}

template<class T> MapFilter_slow<T>::MapFilter_slow( clipper::Xmap<T>& result, const clipper::Xmap<T>& xmap, MapFilterFn_base& fltr, const ftype scale, const TYPE type ) : fltr_( &fltr ), scale_( scale ), type_( type )
{ (*this)( result, xmap ); }

/*! Apply the filter to a given map.
  \param result The filtered map.
  \param xmap The map to be filtered. */
template<class T> bool MapFilter_slow<T>::operator() ( clipper::Xmap<T>& result, const clipper::Xmap<T>& xmap ) const
{
  const MapFilterFn_base& fltr = *fltr_;
  const Grid_sampling& g = xmap.grid_sampling();
  result.init( xmap.spacegroup(), xmap.cell(), g );

  // first determine the effective radius of the radial function
  const int nrad = 1000;
  const ftype drad = 0.25;
  int i;
  ftype r, sum[nrad];
  for ( i = 0; i < nrad; i++ ) {
    r = drad * ( ftype(i) + 0.5 );
    sum[i] = r*r*fabs(fltr(r));
  }
  for ( i = 1; i < nrad; i++ ) sum[i] += sum[i-1];
  for ( i = 0; i < nrad; i++ ) if ( sum[i] > 0.99*sum[nrad-1] ) break;
  ftype rad = drad * ( ftype(i) + 1.0 );

  // now prepare a map from the filter
  Grid_range gm( xmap.cell(), g, rad );
  NXmap<T> flt( xmap.cell(), g, gm );

  // calculate the filter
  typename NXmap<T>::Map_reference_index in;
  ftype64 f000 = 0.0;
  for ( in = flt.first(); !in.last(); in.next() ) {
    r = fltr( sqrt( in.coord_orth().lengthsq() ) );
    f000 += r;
    flt[in] = r;
  }
  // calc scale factor
  ftype scale = 1.0;
  if ( type_ == Absolute ) scale = scale_;
  if ( type_ == Relative ) scale = scale_ / f000;

  // now calculate the search function
  Coord_grid g0, g1;
  typename Xmap<T>::Map_reference_index pos;
  typename Xmap<T>::Map_reference_coord i0, iu, iv, iw;
  ftype s;
  for ( pos = result.first(); !pos.last(); pos.next() ) {
    g0 = pos.coord() + gm.min();
    g1 = pos.coord() + gm.max();
    i0 = Xmap_base::Map_reference_coord( xmap, g0 );
    s = 0.0;
    for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
      for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
	for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
	  in.set_coord( iw.coord() - g0 );
	  s += flt[in] * xmap[iw];
	}
    result[pos] = scale * s;
  }

  return true;
}


/*! The filter function is constructed with the given filter.
  Scaling may also be optionally applied. This may be absolute, which
  just applies a scale factor to the result, or relative, which scales
  the filter relative to its own integral. Therefore, relative scaling
  with a scale factor of 1.0 gives an output map on the same scale as
  the input map.

  Note that the filter is not stored internally, and so must persist
  as long as the MapFilter is required.
  \param fltr The radial filter to apply.
  \param scale The scale factor to apply (default = 1.0).
  \param type The type of scaling to apply: NONE, Absolute, Relative. */
template<class T> MapFilter_fft<T>::MapFilter_fft( const MapFilterFn_base& fltr, const ftype scale, const TYPE type ) : fltr_( &fltr ), scale_( scale ), type_( type )
{}

template<class T> MapFilter_fft<T>::MapFilter_fft( clipper::Xmap<T>& result, const clipper::Xmap<T>& xmap, MapFilterFn_base& fltr, const ftype scale, const TYPE type ) : fltr_( &fltr ), scale_( scale ), type_( type )
{ (*this)( result, xmap ); }

/*! Apply the filter to a given map.
  \param result The filtered map.
  \param xmap The map to be filtered. */
template<class T> bool MapFilter_fft<T>::operator() ( clipper::Xmap<T>& result, const clipper::Xmap<T>& xmap ) const
{
  const MapFilterFn_base& fltr = *fltr_;
  const Grid_sampling& g = xmap.grid_sampling();

  // first determine the effective radius of the radial function
  const int nrad = 1000;
  const ftype drad = 0.25;
  int i;
  ftype r, sum[nrad];
  for ( i = 0; i < nrad; i++ ) {
    r = drad * ( ftype(i) + 0.5 );
    sum[i] = r*r*fabs(fltr(r));
  }
  for ( i = 1; i < nrad; i++ ) sum[i] += sum[i-1];
  for ( i = 0; i < nrad; i++ ) if ( sum[i] > 0.99*sum[nrad-1] ) break;
  ftype rad = drad * ( ftype(i) + 1.0 );

  // make the fft maps
  FFTmap_p1 map( g );
  FFTmap_p1 flt( g );
  // fill the fft maps
  Coord_grid c, half( g.nu()/2, g.nv()/2, g.nw()/2 );
  typename Xmap<T>::Map_reference_coord i0( xmap, Coord_grid(0,0,0) );
  typename Xmap<T>::Map_reference_coord iu, iv, iw;
  ftype64 f000 = 0.0;
  for ( iu = i0; iu.coord().u() < g.nu(); iu.next_u() )
    for ( iv = iu; iv.coord().v() < g.nv(); iv.next_v() )
      for ( iw = iv; iw.coord().w() < g.nw(); iw.next_w() ) {
	c = (iw.coord() + half).unit(g) - half;
	r = sqrt( c.coord_frac(g).lengthsq(xmap.cell()) );
	map.real_data( iw.coord() ) = xmap[iw];
	if ( r < rad ) {
	  r = fltr(r);
	  f000 += r;
	  flt.real_data( iw.coord() ) = r;
	}
      }
  // calc scale factor
  ftype32 scale = 1.0;
  if ( type_ == Absolute ) scale = scale_;
  if ( type_ == Relative ) scale = scale_ / f000;

  // fft
  flt.fft_x_to_h( xmap.cell().volume() );
  map.fft_x_to_h( xmap.cell().volume() );
  // do filter
  const Grid& gh = map.grid_reci();
  for ( c.u() = 0; c.u() < gh.nu(); c.u()++ )
    for ( c.v() = 0; c.v() < gh.nv(); c.v()++ )
      for ( c.w() = 0; c.w() < gh.nw(); c.w()++ )
	map.cplx_data( c ) = scale * map.cplx_data(c) * flt.cplx_data(c);
  // invert
  map.fft_h_to_x( map.grid_real().size() / pow( xmap.cell().volume(), 2 ) );

  // store
  result.init( xmap.spacegroup(), xmap.cell(), g );
  for ( typename Xmap<T>::Map_reference_index ix = result.first();
	!ix.last(); ix.next() )
    result[ix] = map.real_data( ix.coord() );

  return true;
}

/*! Apply the filter to a given map.
  \param result The filtered map.
  \param nxmap The map to be filtered. */
template<class T> bool MapFilter_fft<T>::operator() ( clipper::NXmap<T>& result, const clipper::NXmap<T>& nxmap ) const
{
  const MapFilterFn_base& fltr = *fltr_;
  const Grid& gnx = nxmap.grid();

  // Now make fft grid compatible with this NXmap grid
  Grid_sampling nuvw( gnx.nu()+2, gnx.nv()+2, gnx.nw()+2 );
  for ( int i = 0; i < 3; i++ ) {  // find most efficient sampling
    int j, l, m, n, nbest = 0;
    ftype t, tbest = 1.0e12;
    for ( n = nuvw[i]; n < 2*nuvw[i]; n++ ) {
      if ( n % 2 == 0 ) {
        l = 0; // sum of factors (approx. log n)
        m = n; // what is left in factorisation
        for ( j = 2; j <= n; j++ )
          while ( m%j == 0 ) { m /= j; l += j; }
        // FFT time O( n * l ) - introduce an extra n^2 for 3D.
        t = pow( ftype(n), 3 ) * ftype(l);
        if ( t < tbest ) { nbest = n; tbest = t; }
      }
    }
    nuvw[i] = nbest;
  }
  const Grid_sampling& g = nuvw;

  // Get grid skew matrix
  nxmap.operator_grid_orth().rot(); Coord_orth::zero();
  RTop<> rtop( nxmap.operator_grid_orth().rot(), Coord_orth::zero() );

  // make the fft maps
  FFTmap_p1 map( g );
  FFTmap_p1 flt( g );

  // fill fft map from nx map
  Coord_grid l( (g.nu()-gnx.nu())/2,(g.nv()-gnx.nv())/2,(g.nw()-gnx.nw())/2 );
  Coord_grid h( l + gnx - Coord_grid(1,1,1) );
  typedef NXmap<float>::Map_reference_index MRI;
  for ( MRI ix = nxmap.first(); !ix.last(); ix.next() )
    map.real_data( ix.coord()+l ) = nxmap[ix];
  // mirror the remaining points
  Coord_grid iw, im;
  for ( iw.u() = 0; iw.u() < g.nu(); iw.u()++ )
    for ( iw.v() = 0; iw.v() < g.nv(); iw.v()++ )
      for ( iw.w() = 0; iw.w() < g.nw(); iw.w()++ ) {
	im = iw;
	if ( im.u() < l.u() ) im.u() = l.u() + ( l.u() - im.u() );
	if ( im.v() < l.v() ) im.v() = l.v() + ( l.v() - im.v() );
	if ( im.w() < l.w() ) im.w() = l.w() + ( l.w() - im.w() );
	if ( im.u() > h.u() ) im.u() = h.u() + ( h.u() - im.u() );
	if ( im.v() > h.v() ) im.v() = h.v() + ( h.v() - im.v() );
	if ( im.w() > h.w() ) im.w() = h.w() + ( h.w() - im.w() );
	map.real_data( iw ) = map.real_data( im );
      }

  // first determine the effective radius of the radial function
  const int nrad = 1000;
  const ftype drad = 0.25;
  int i;
  ftype r, sum[nrad];
  for ( i = 0; i < nrad; i++ ) {
    r = drad * ( ftype(i) + 0.5 );
    sum[i] = r*r*fabs(fltr(r));
  }
  for ( i = 1; i < nrad; i++ ) sum[i] += sum[i-1];
  for ( i = 0; i < nrad; i++ ) if ( sum[i] > 0.99*sum[nrad-1] ) break;
  ftype rad = drad * ( ftype(i) + 1.0 );

  // fill the radial function map
  Coord_grid c, half( g.nu()/2, g.nv()/2, g.nw()/2 );
  ftype64 f000 = 0.0;
  for ( iw.u() = 0; iw.u() < g.nu(); iw.u()++ )
    for ( iw.v() = 0; iw.v() < g.nv(); iw.v()++ )
      for ( iw.w() = 0; iw.w() < g.nw(); iw.w()++ ) {
	c = (iw + half).unit(g) - half;
	r = sqrt( Coord_orth(rtop*c.coord_map()).lengthsq() );
	if ( r < rad ) {
	  r = fltr(r);
	  f000 += r;
	  flt.real_data( iw ) = r;
	}
      }

  // calc scale factor
  ftype32 scale = 1.0;
  if ( type_ == Absolute ) scale = scale_;
  if ( type_ == Relative ) scale = scale_ / f000;

  // fft
  flt.fft_x_to_h( 1.0 );
  map.fft_x_to_h( 1.0 );
  // do filter
  const Grid& gh = map.grid_reci();
  for ( c.u() = 0; c.u() < gh.nu(); c.u()++ )
    for ( c.v() = 0; c.v() < gh.nv(); c.v()++ )
      for ( c.w() = 0; c.w() < gh.nw(); c.w()++ )
	map.cplx_data( c ) = scale * map.cplx_data(c) * flt.cplx_data(c);
  // invert
  map.fft_h_to_x( map.grid_real().size() );

  // store
  result.init( gnx, nxmap.operator_orth_grid() );
  for ( MRI ix = nxmap.first(); !ix.last(); ix.next() )
    result[ix] = map.real_data( ix.coord()+l );

  return true;
}

// compile templates

template class CLIPPER_IMEX MapFilter_slow<ftype32>;
template class CLIPPER_IMEX MapFilter_slow<ftype64>;
template class CLIPPER_IMEX MapFilter_fft<ftype32>;
template class CLIPPER_IMEX MapFilter_fft<ftype64>;


} // namespace clipper
