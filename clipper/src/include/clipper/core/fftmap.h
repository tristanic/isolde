/*! \file lib/fftmap.h
    Header file for P1 fft map
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


#ifndef CLIPPER_FFTMAP
#define CLIPPER_FFTMAP


#include "coords.h"
#include "../imex.h"
#include <complex>

namespace clipper
{
  // fft type
  typedef float ffttype;

  // forward definition
  namespace datatypes {
    template<class T> class F_phi;
  }

  // base class for FFT classes
  class CLIPPER_IMEX FFTmap_base {
  public:
    enum FFTtype { Default, Measure, Estimate };  //!< optimisation options
  protected:
    static Mutex mutex;                 //!< Thread safety
  };

  //! FFTmap_p1: low level P1 map used for calculating FFTs
  /*! This is a pure real P1 map, with an extra section in reciprocal
    space to allow generation of the full set of resiprocal space
    magnitudes. Access is by Coord_grid in both spaces, and indices
    must be non-negative and in range. The first and last sections
    along the half-length direction only have half the elements
    stored, the contents of the other half is ignored.
  */
  class CLIPPER_IMEX FFTmap_p1 : public FFTmap_base
  {
  public:
    //! Null constructor
    FFTmap_p1();
    //! Copy constructor
    FFTmap_p1( const FFTmap_p1& other ) { copy(other); }
    //! Constructor: takes grid
    FFTmap_p1( const Grid_sampling& grid_sam, const FFTtype type = Default );
    //! Assignment operator
    const FFTmap_p1& operator =(const FFTmap_p1& other) { return copy(other); }
    //! initialiser: takes grid
    void init( const Grid_sampling& grid_sam, const FFTtype type = Default );
    //! Reset
    void reset();

    //! Return real space grid.
    const Grid_sampling& grid_real() const { return grid_sam_; }
    //! Return reciprocal space grid (i.e. half real grid + 1 section).
    const Grid& grid_reci() const { return grid_reci_; }
    //! Test whether a coordinate is in the valid part of the recip. grid.
    bool uniq_reci( const Coord_grid& c ) const { return ( (c.w()>0 && c.w()<grid_half_.nw()) || (c.w()<=grid_half_.nw() && ( (c.v()>0 && c.v()<grid_half_.nv()) || (c.v()<=grid_half_.nv() && c.u()<=grid_half_.nu()) ) ) ); }

    //! Transform to real space
    void fft_h_to_x( const ftype& scale );
    //! Transform to reciprocal space
    void fft_x_to_h( const ftype& scale );

    //! get reciprocal space data: slow form with hemisphere check
    std::complex<ffttype> get_hkl( const HKL& hkl ) const;
    //! set reciprocal space data: slow form with hemisphere check
    void set_hkl( const HKL& hkl, const std::complex<ffttype>& f );
    //! get reciprocal space data
    const std::complex<ffttype>& cplx_data( const Coord_grid& c ) const
      { return data_c[ grid_reci_.index( c ) ]; }
    //! set reciprocal space data
    std::complex<ffttype>& cplx_data( const Coord_grid& c )
      { return data_c[ grid_reci_.index( c ) ]; }
    //! get real space data
    const ffttype& real_data( const Coord_grid& c ) const
      { return data_r[ grid_real_.index( c ) ]; }
    //! set real space data
    ffttype& real_data( const Coord_grid& c )
      { return data_r[ grid_real_.index( c ) ]; }

    //! set/get default optimisation type
    static FFTtype& default_type() { return default_type_; }

    void debug() const;

  protected:
    const FFTmap_p1& copy( const FFTmap_p1& other );  //!< copy function

    enum FFTmode { NONE, RECI, REAL, OTHER };  //!< space enumeration

    FFTmode mode;                       //!< real or reciprocal space?
    FFTtype type_;                      //!< optimisation options
    Grid_sampling grid_sam_;            //!< unit cell grid
    Grid grid_reci_;                    //!< reciprocal space grid
    Grid grid_real_;                    //!< real space grid
    Grid grid_half_;                    //!< half grid (for marking unique)

    Matrix<char> req_kl, req_uv;        //!< reci section lookup
    std::vector<char> req_l, req_u;     //!< real section lookup

    std::vector<ffttype> datavec;       //!< vector for the data
    ffttype* data_r;                    //!< pointer to real data
    std::complex<ffttype>* data_c;      //!< pointer to complex data

    static FFTtype default_type_;       //!< default optimisation options
  };


  //! FFTmap: P1 map with symmetry used for calculating FFTs
  /*! The FFTmap is represented in P1 in memory. However, it also has
    a spacegroup, and the contained data remains consistent with this
    spacegroup at all times. It has three states - unassigned,
    real-space, and reciprocal space. In real space it contains real
    map data. In reciprocal space it holds a hemisphere of complex
    structure factors, with the Friedels duplicated on the zero
    section.

    The user should be able to ignore all the issues of spacegroup
    symmetry, Friedel opposites, and storage order.
  */
  class CLIPPER_IMEX FFTmap : private FFTmap_p1
  {
  public:
    //! Null constructor
    FFTmap();
    //! Constructor: takes spacegroup, cell, grid
    FFTmap( const Spacegroup& spacegroup, const Cell& cell, const Grid_sampling grid_sam, const FFTtype type = Default );
    //! initialiser
    void init( const Spacegroup& spacegroup, const Cell& cell, const Grid_sampling grid_sam, const FFTtype type = Default );
    //! Reset
    void reset();

    //! get the cell
    const Cell& cell() const { return cell_; }
    //! get the spacegroup
    const Spacegroup& spacegroup() const { return spacegroup_; }
    //! get the cell grid
    const Grid_sampling& grid_sampling() const { return FFTmap_p1::grid_real(); }
    //! Transform to real space
    void fft_h_to_x();
    //! Transform to reciprocal space
    void fft_x_to_h();

    //! get reciprocal space data
    template<class T> void get_recip_data( const HKL& rfl, datatypes::F_phi<T>& fphi ) const;
    //! set reciprocal space data
    template<class T> void set_recip_data( const HKL& rfl, const datatypes::F_phi<T>& fphi );

    //! get real space data
    template<class T> void get_real_data( const Coord_grid& c, T& datum ) const;
    //! set real space data
    template<class T> void set_real_data( const Coord_grid& c, const T& datum );

    //! get reciprocal space data (No error checking)
    datatypes::F_phi<ffttype> get_recip_data( const HKL& rfl ) const;
    //! get real space data (No error checking)
    const ffttype& get_real_data( const Coord_grid& c ) const
      { return real_data(c.unit(grid_real())); }

    //! calculate map-like object from reflection-like object
    template<class H, class X> void fft_rfl_to_map( const H& h, X& x );
    //! calculate reflection-like object from map-like object
    template<class H, class X> void fft_map_to_rfl( const X& x, H& h );

    void debug() const;

  protected:
    Cell cell_;                         //!< unit cell
    Spacegroup spacegroup_;             //!< spacegroup
    std::vector<Isymop> isymop;         //!< Integerised symops
  };


  // template implementations


  /*! Fill this FFTmap object from a reflection object, transform it,
    and fill the given map object from the FFTmap. This will work for
    any reflection data object which implements a HKL_reference_index,
    and every map data object which implements a Map_reference_index.

    For the results to be sensible, the spacegroup, cell and grids
    should match. (The map will be zeroed if necessary).

    \param h The source reflection data object.
    \param x The target map object.
  */
  template<class H, class X> void FFTmap::fft_rfl_to_map( const H& h, X& x )
  {
    // zero the map
    reset();

    // copy from reflection data
    typename H::HKL_reference_index ih;
    for ( ih = h.first_data(); !ih.last(); h.next_data( ih ) )
      set_recip_data( ih.hkl(), h[ih] );

    // fft
    fft_h_to_x();

    // and copy into the map
    typename X::Map_reference_index ix;
    for ( ix = x.first(); !ix.last(); ix.next() )
      get_real_data( ix.coord(), x[ix] );
  }


  /*! Fill this FFTmap object from a map object, transform it, and
    fill the given reflection object from the FFTmap. This will work
    for any reflection data object which implements a
    HKL_reference_index, and every map data object which implements a
    Map_reference_index.

    For the results to be sensible, the spacegroup, cell and grids
    should match. (The map will be zeroed if necessary).

    \param x The source map object.
    \param h The target reflection data object.
  */
  template<class H, class X> void FFTmap::fft_map_to_rfl( const X& x, H& h )
  {
    // zero the map
    reset();

    // copy into the map
    typename X::Map_reference_index ix;
    for ( ix = x.first(); !ix.last(); ix.next() )
      set_real_data( ix.coord(), x[ix] );

    // fft
    fft_x_to_h();

    // now fill it
    typename H::HKL_reference_index ih;
    for ( ih = h.first(); !ih.last(); ih.next() )
      get_recip_data( ih.hkl(), h[ih] );
  }


} // namespace clipper

#endif
