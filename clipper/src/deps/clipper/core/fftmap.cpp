/* fftmap.cpp: implementation file for P1 fft map */
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


#include "fftmap.h"

#include "hkl_datatypes.h"

#include <config.h>
#ifdef FFTW2_PREFIX_S
# include <srfftw.h>
#else
# include <rfftw.h>
#endif

// compile-time check if fftw above is really single-precision
static float* dummy = (fftw_real*) NULL;


namespace clipper {

Message_fatal message_fftmap_get_real_space_error( "FFTmap: get_real_data in reciprocal space" );
Message_fatal message_fftmap_set_real_space_error( "FFTmap: set_real_data in reciprocal space" );
Message_fatal message_fftmap_get_reci_space_error( "FFTmap: get_recip_data in real space" );
Message_fatal message_fftmap_set_reci_space_error( "FFTmap: set_recip_data in real space" );
Message_ctor message_ctor_fftmap( " [FFTmap: constructed]" );


FFTmap_base::FFTtype FFTmap_p1::default_type_ = FFTmap_base::Estimate;
Mutex FFTmap_base::mutex = Mutex();


/*! For later initialisation: see init() */
FFTmap_p1::FFTmap_p1()
{}

/*! Construct an FFTmap_p1 for a given spacegroup, cell, and grid.  The
  map values are initialised to zero.

  The FFTmap_p1 is initially in neither real nor reciprocal spce, however
  as soon as one of the 'set' methods is called, it will be defined as
  in either real or reciprocal space until the next fft.
  \param grid_sam The grid sampling of the unit cell.
  \param type Can be FFTmap_p1::Measure, FFTmap_p1::Estimate.
  Measure performs slow precalculation (first time only) to get faster FFT. */
FFTmap_p1::FFTmap_p1( const Grid_sampling& grid_sam, const FFTtype type )
{ init( grid_sam, type ); }

/*! Initialise an FFTmap_p1 for a given spacegroup, cell, and grid.  The
  map values are initialised to zero.

  The FFTmap_p1 is initially in neither real nor reciprocal spce, however
  as soon as one of the 'set' methods is called, it will be defined as
  in either real or reciprocal space until the next fft.
  \param grid_sam The grid sampling of the unit cell.
  \param type Can be FFTmap_p1::Measure, FFTmap_p1::Estimate.
  Measure performs slow precalculation (first time only) to get faster FFT. */
void FFTmap_p1::init( const Grid_sampling& grid_sam, const FFTtype type )
{
  grid_sam_ = grid_sam;
  type_ = type;
  if ( type_ == Default ) type_ = default_type();

  // to start out with, we don't know which space we are in
  mode = NONE;

  // allocate data
  grid_reci_ = Grid( grid_sam_.nu(), grid_sam_.nv(), grid_sam_.nw()/2+1 );
  grid_real_ = Grid( grid_reci_.nu(), grid_reci_.nv(), grid_reci_.nw()*2 );
  grid_half_ = Grid( grid_sam_.nu()/2, grid_sam_.nv()/2, grid_sam_.nw()/2 );
  datavec.resize( grid_real_.size(), 0.0 );

  // set pointers to data
  data_r = &datavec[0];
  data_c = (std::complex<ffttype>*)data_r;
}

/*! Reset the space and zero all the data, if necessary. */
void FFTmap_p1::reset()
{
  mode = NONE;

  std::vector<ffttype>::iterator i;
  for ( i = datavec.begin(); i != datavec.end(); i++ ) *i = 0.0;
}

/*! \fn const Grid_sampling& FFTmap_p1::grid_real() const
  \return The grid sampling of the real space grid. */

/*! \fn const Grid& FFTmap_p1::grid_reci() const
  The reciprocal grid is half-length, plus one section, in the w
  direction. The remainder of the grid may be generated by Hermitian
  symmetry. When accessing data with reci_data, the coordinate should
  always be in this grid. Some points in this grid are redundent, see
  FFTmap_p1::uniq_reci().
  \return The reciprocal grid. */

/*! \fn bool FFTmap_p1::uniq_reci( const Coord_grid& c ) const
  The w=0 and w=nw/2 sections contain some duplicated points related
  by a cetre of symmetry. On of these is considered to be significant,
  and the other redundent. This function returns 'true' for the
  significant point.  \note For some calculations it may be quicker to
  set the whole grid than call this function for every coordinate.
  \param c The coordinate to test. Must be in grid_reci().
  \return true if the coordinate is for a significant point. */

/*! The data is transformed from recirocal to real space. If the FFTmap_p1
  is already in real space, no action is taken.
  \param Scale factor to apply
  (normally 1/cell_volume). */
void FFTmap_p1::fft_h_to_x( const ftype& scale )
{
  if ( mode == REAL ) return;
  // scale
  ffttype s = ffttype( scale );
  int n = grid_reci_.size();
  for ( int i = 0; i < n; i++ ) data_c[i] = std::conj( s * data_c[i] );
  // fft
  int flags = ( type_ == Measure ) ?
    ( FFTW_IN_PLACE | FFTW_USE_WISDOM | FFTW_MEASURE ) :
    ( FFTW_IN_PLACE | FFTW_USE_WISDOM | FFTW_ESTIMATE );

  mutex.lock();
  fftwnd_plan plan =
    rfftw3d_create_plan( grid_sam_.nu(), grid_sam_.nv(), grid_sam_.nw(),
			 FFTW_COMPLEX_TO_REAL, flags );
  mutex.unlock();

  rfftwnd_one_complex_to_real(plan, (fftw_complex*)data_c, NULL);

  mutex.lock();
  rfftwnd_destroy_plan(plan);
  mutex.unlock();

  // done
  mode = REAL;
}

/*! The data is transformed from real to recirocal space. If the
  FFTmap_p1 is already in reciproal space, no action is taken.
  \param Scale factor to apply (in addition to 1/N_grid factor)
  (normally cell_volume). */
void FFTmap_p1::fft_x_to_h( const ftype& scale )
{
  if ( mode == RECI ) return;
  // fft
  int flags = ( type_ == Measure ) ?
    ( FFTW_IN_PLACE | FFTW_USE_WISDOM | FFTW_MEASURE ) :
    ( FFTW_IN_PLACE | FFTW_USE_WISDOM | FFTW_ESTIMATE );

  mutex.lock();
  fftwnd_plan plan =
    rfftw3d_create_plan( grid_sam_.nu(), grid_sam_.nv(), grid_sam_.nw(),
			 FFTW_REAL_TO_COMPLEX, flags );
  mutex.unlock();

  rfftwnd_one_real_to_complex(plan, (fftw_real*)data_r, NULL);

  mutex.lock();
  rfftwnd_destroy_plan(plan);
  mutex.unlock();

  // scale
  ffttype s = ffttype( scale ) / grid_sam_.size();
  int n = grid_reci_.size();
  for ( int i = 0; i < n; i++ ) data_c[i] = std::conj( s * data_c[i] );
  // done
  mode = RECI;
}

/*! This form returns the data for an HKL. The HKL is converted into a
  grid reference, and the data, or if necessary the conjugate of the
  opposite, is returned.
  \param hkl The HKL of the data. */
std::complex<ffttype> FFTmap_p1::get_hkl( const HKL& hkl ) const
{
  Coord_grid c = Coord_grid( hkl ).unit( grid_sam_ );
  if ( c.w() < grid_reci_.nw() )
    return cplx_data(c);
  else
    return std::conj( cplx_data( Coord_grid(
	  (grid_sam_.nu()-c.u())%grid_sam_.nu(),
	  (grid_sam_.nv()-c.v())%grid_sam_.nv(),
	  (grid_sam_.nw()-c.w())%grid_sam_.nw() ) ) );
}

/*! This form returns the data for an HKL. The HKL is converted into a
  grid reference, and the data, and if necessary the conjugate of the
  opposite, is set.
  \param hkl The HKL of the data. */
void FFTmap_p1::set_hkl( const HKL& hkl, const std::complex<ffttype>& f )
{
  Coord_grid c;
  c = Coord_grid( hkl ).unit( grid_sam_ );
  if ( c.w() < grid_reci_.nw() )
    cplx_data(c) = f;
  c = Coord_grid( -hkl ).unit( grid_sam_ );
  if ( c.w() < grid_reci_.nw() )
    cplx_data(c) = std::conj(f);
}

void FFTmap_p1::debug() const
{
  Coord_grid c;
  int i, j, k; i = j = k = 0;
  for ( c.u() = 0; c.u() < grid_sam_.nu(); c.u()++ )
    for ( c.v() = 0; c.v() < grid_sam_.nv(); c.v()++ )
      for ( c.w() = 0; c.w() < grid_sam_.nw(); c.w()++ ) {
	if ( uniq_reci( c ) && uniq_reci( (-c).unit(grid_sam_) ) ) i++;
	if ( uniq_reci( c ) ) j++;
	k++;
      }
  std::cout << "FFTmap_p1 debug: " << i << "\t" << j << "\t" << k << "\n";
  grid_sam_.debug();
  for ( int i = 0; i < datavec.size(); i++ )
    std::cout << i << " " << datavec[i] << "\n";
}

const FFTmap_p1& FFTmap_p1::copy( const FFTmap_p1& other )
{
  mode = other.mode;
  type_ = other.type_;
  grid_sam_ = other.grid_sam_;
  grid_reci_ = other.grid_reci_;
  grid_real_ = other.grid_real_;
  grid_half_ = other.grid_half_;
  req_kl = other.req_kl;
  req_uv = other.req_uv;
  req_l = other.req_l;
  req_u = other.req_u;
  datavec = other.datavec;
  data_r = &datavec[0];
  data_c = (std::complex<ffttype>*)data_r;
  return *this;
}


/*! For later initialisation: see init() */
FFTmap::FFTmap()
{ Message::message( message_ctor_fftmap ); }

/*! Construct an FFTmap for a given spacegroup, cell, and grid.  The
  map values are initialised to zero.

  The FFTmap is initially in neither real nor reciprocal spce, however
  as soon as one of the 'set' methods is called, it will be defined as
  in either real or reciprocal space until the next fft.
  \param spacegroup The spacegroup.
  \param cell The cell, used for scaling.
  \param grid_sam The grid sampling of the unit cell.
  \param precalc Perform slow precalculation to get faster FFT. (default: no) */
FFTmap::FFTmap( const Spacegroup& spacegroup, const Cell& cell, const Grid_sampling grid_sam, const FFTtype type )
{
  Message::message( message_ctor_fftmap );
  init( spacegroup, cell, grid_sam, type );
}

/*! Initialise an FFTmap for a given spacegroup, cell, and grid.  The
  map values are initialised to zero.

  The FFTmap is initially in neither real nor reciprocal spce, however
  as soon as one of the 'set' methods is called, it will be defined as
  in either real or reciprocal space until the next fft.
  \param spacegroup The spacegroup.
  \param cell The cell, used for scaling.
  \param grid_sam The grid sampling of the unit cell.
  \param precalc Perform slow precalculation to get faster FFT.
  This adds a penalty of about 4s on Linux for the first FFT of any
  grid and direction. Subsequent FFTs will be faster. Set to true for
  programs which will use many FFTs. default: false.*/
void FFTmap::init( const Spacegroup& spacegroup, const Cell& cell, const Grid_sampling grid_sam, const FFTtype type )
{
  FFTmap_p1::init( grid_sam, type );
  spacegroup_ = spacegroup;
  cell_ = cell;

  // Create the intergised symops (assumes legal grid)
  isymop.resize( spacegroup.num_symops() );
  for ( int sym = 0; sym < spacegroup.num_symops(); sym++ )
    isymop[sym] = Isymop( spacegroup.symop(sym), grid_sam );
}

/*! Reset the space and zero all the data, if necessary. */
void FFTmap::reset()
{
  mode = NONE;

  std::vector<ffttype>::iterator i;
  for ( i = datavec.begin(); i != datavec.end(); i++ ) *i = 0.0;
}

/*! The data is transformed from recirocal to real space. A scale
  factor of 1/v (where v is the cell volume) is applied. If the FFTmap
  is already in real space, no action is taken. */
void FFTmap::fft_h_to_x()
{
  if ( mode != RECI ) return;
  FFTmap_p1::fft_h_to_x( 1.0/cell().volume() );
}

/*! The data is transformed from real to recirocal space. A scale
  factor of v/n (where v is the cell volume and n the number of grid
  points) is applied. If the FFTmap is already in reciproal space, no
  action is taken. */
void FFTmap::fft_x_to_h()
{
  if ( mode != REAL ) return;
  FFTmap_p1::fft_x_to_h( cell().volume() );
}

/*! The data value for the given HKL, or the conjugate of its Friedel
  opposite if required, is returned. The symmetry related copies of
  the data are ignored.
  \param rfl The HKL of the data to be returned.
  \param fphi The value, as a magnitude and phase of type \c ffttype */
template<class T> void FFTmap::get_recip_data( const HKL& rfl, datatypes::F_phi<T>& fphi ) const
{
  if ( mode != RECI ) Message::message(message_fftmap_get_reci_space_error);
  fphi = std::complex<T>( FFTmap_p1::get_hkl( rfl ) );
}

/*! The data value for the given HKL, or the conjugate of its Friedel
  opposite if required, is set. All the symmetry related copies of
  the data, and any Friedel copies in the zero section, are also set.
  \param rfl The HKL of the data to be set.
  \param fphi The value, as a magnitude and phase of type \c ffttype */
template<class T> void FFTmap::set_recip_data( const HKL& rfl, const datatypes::F_phi<T>& fphi )
{
  // check space
  if ( mode != RECI ) {
    if ( mode == NONE ) mode = RECI;
    else Message::message(message_fftmap_set_reci_space_error);
  }

  // store all sym copies of reflection
  // could use the F_phi type to do sym stuff, but this is faster
  T phi = fphi.phi();
  FFTmap_p1::set_hkl( rfl, std::complex<ffttype>( fphi.f() * cos( phi ), fphi.f() * sin( phi ) ) );
  for ( int sym = 1; sym < spacegroup_.num_primops(); sym++ ) {
    phi = fphi.phi() + rfl.sym_phase_shift( spacegroup_.symop(sym) );
    FFTmap_p1::set_hkl( rfl.transform( isymop[sym] ), std::complex<ffttype>( fphi.f() * cos( phi ), fphi.f() * sin( phi ) ) );
  }
}

/*! The data value for the given grid coordinate is returned. Symmetry
  related copies are ignored.
  \param c The coordinate of the data to be returned.  
  \param datum The value of the data. */
template<class T> void FFTmap::get_real_data( const Coord_grid& c, T& datum ) const
{
  if ( mode != REAL ) Message::message(message_fftmap_get_real_space_error);
  datum = real_data(c.unit(grid_real()));
}

/*! The data value for the given grid coordinate is set. All the
  symmetry related copies of the data are also set.
  \param c The coordinate of the data to be set.  
  \param datum The value of the data. */
template<class T> void FFTmap::set_real_data( const Coord_grid& c, const T& datum )
{
  // check space
  if ( mode != REAL ) {
    if ( mode == NONE ) mode = REAL;
    else Message::message(message_fftmap_set_real_space_error);
  }

  // set all sym copies of map value
  real_data( c.unit(grid_sam_) ) = ffttype( datum );
  for ( int sym = 1; sym < isymop.size(); sym++ )
    real_data( c.transform(isymop[sym]).unit(grid_sam_) ) = ffttype( datum );
}

/*! \fn datatypes::F_phi<ffttype> FFTmap::get_recip_data( const HKL& rfl ) const
  No error is produced if the space is wrong.
  \param rfl The HKL of the data to be returned.
  \return The value, as magnitude and phase of type \c ffttype */
datatypes::F_phi<ffttype> FFTmap::get_recip_data( const HKL& rfl ) const
{ return datatypes::F_phi<ffttype>( get_hkl(rfl) ); }

/*! \fn const ffttype& FFTmap::get_real_data( const Coord_grid& c ) const
  No error is produced if the space is wrong.
  \param c The grid coordinate of the data to be returned.
  \return The value, as type \c ffttype */

// template instantiations

template void FFTmap::get_recip_data<ftype32>( const HKL& rfl, datatypes::F_phi<ftype32>& fphi ) const;
template void FFTmap::set_recip_data<ftype32>( const HKL& rfl, const datatypes::F_phi<ftype32>& fphi );
template void FFTmap::get_real_data<ftype32>( const Coord_grid& c, ftype32& datum ) const;
template void FFTmap::set_real_data<ftype32>( const Coord_grid& c, const ftype32& datum );
template void FFTmap::get_recip_data<ftype64>( const HKL& rfl, datatypes::F_phi<ftype64>& fphi ) const;
template void FFTmap::set_recip_data<ftype64>( const HKL& rfl, const datatypes::F_phi<ftype64>& fphi );
template void FFTmap::get_real_data<ftype64>( const Coord_grid& c, ftype64& datum ) const;
template void FFTmap::set_real_data<ftype64>( const Coord_grid& c, const ftype64& datum );

} // namespace clipper
