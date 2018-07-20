/*! \file lib/fftmap_sparse.h
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


#ifndef CLIPPER_FFTMAP_SPARSE
#define CLIPPER_FFTMAP_SPARSE


#include "fftmap.h"
#include "../imex.h"

namespace clipper
{

  //! base type for sparse P1 fft maps
  class CLIPPER_IMEX FFTmap_sparse_p1_base : public FFTmap_base {
  public:
    //! initialiser: takes grid
    void init( const Grid_sampling& grid_sam, const FFTtype type = Default );
    //! Destructor
    ~FFTmap_sparse_p1_base();
    //! get real grid sampling
    const Grid_sampling& grid_real() const { return grid_real_; }
    //! get reciprocal grid
    const Grid& grid_reci() const { return grid_reci_; }

    //! set/get default optimisation type
    static FFTtype& default_type() { return default_type_; }
  protected:
    //! return/create row
    ffttype* map_uv( const int& u, const int& v );
    //! return/create row
    std::complex<ffttype>* map_kl( const int& k, const int& l );

    Grid_sampling grid_real_;           //!< real space grid
    Grid grid_reci_;                    //!< reciprocal space grid
    FFTtype type_;                      //!< optimisation options

    Array2d<std::complex<ffttype>*> row_kl; //!< section map
    Array2d<ffttype*> row_uv;               //!< section map

    static FFTtype default_type_;       //!< default optimisation options
  };

  //! FFTmap_sparse_p1_hx: low level sparse P1 map used for calculating FFTs
  /*! This version computes sparse Hermititan...real FFTs.

    By specifying what parts of the map are needed in advance, it is
    possible to perform highly optimised FFTs, including some of the
    benefits of symmetry. */
  class CLIPPER_IMEX FFTmap_sparse_p1_hx : public FFTmap_sparse_p1_base
  {
  public:
    //! Null constuctor
    FFTmap_sparse_p1_hx();
    //! Constructor: takes grid
    FFTmap_sparse_p1_hx( const Grid_sampling& grid_sam, const FFTtype type = Default );
    //-- void init( const Grid_sampling& grid_sam, const FFTtype type = Default );
    //-- const Grid_sampling& grid_real() const;
    //-- const Grid& grid_reci() const;

    //! set reciprocal space data by hkl
    void set_hkl( const HKL& hkl, const std::complex<ffttype>& f );
    //! set reciprocal space data (internal use)
    std::complex<ffttype>& cplx_data( const Coord_grid& uvw )
      { return map_kl( uvw.v(), uvw.w() )[ uvw.u() ]; }
    //! express need for real space data
    void require_real_data( const Coord_grid& uvw )
      { map_uv( uvw.u(), uvw.v() ); }
    //! get real space data ( uvw must be in grid_real() )
    const ffttype& real_data( const Coord_grid& uvw ) const
      { return row_uv( uvw.u(), uvw.v() )[ uvw.w() ]; }

    //! Transform to real space
    void fft_h_to_x( const ftype& scale );

private:
    void transform_along_hu_(void* planu_ptr, const int& start, const int& end);
    void transform_along_kv_(void* planv_ptr, const int& start, const int& end,
        const ffttype& s, const int& nmax);
    void transform_along_lw_(void* planw_ptr, const int& start, const int& end);
    int layers_per_thread_ = 25;
    std::vector<bool> map_l;
    std::vector<bool> row_u;
  };

  //! FFTmap_sparse_p1_xh: low level sparse P1 map used for calculating FFTs
  /*! This version computes sparse real...Hermititan FFTs.

    By specifying what parts of the map are needed in advance, it is
    possible to perform highly optimised FFTs, including some of the
    benefits of symmetry. */
  class CLIPPER_IMEX FFTmap_sparse_p1_xh : public FFTmap_sparse_p1_base
  {
  public:
    //! Null constuctor
    FFTmap_sparse_p1_xh();
    //! Constructor: takes grid
    FFTmap_sparse_p1_xh( const Grid_sampling& grid_sam, const FFTtype type = Default );
    //-- void init( const Grid_sampling& grid_sam, const FFTtype type = Default );
    //-- const Grid_sampling& grid_real() const;
    //-- const Grid& grid_reci() const;

    //! set real space data ( uvw must be in grid_real() )
    ffttype& real_data( const Coord_grid& uvw )
      { return map_uv( uvw.u(), uvw.v() )[ uvw.w() ]; }
    //! express need for reciprocal space data by hkl
    void require_hkl( const HKL& hkl );
    //! get reciprocal space data by hkl
    const std::complex<ffttype> get_hkl( const HKL& hkl ) const;
    //! express need for reciprocal space data (internal use)
    void require_cplx_data( const Coord_grid& hkl )
      { map_kl( hkl.v(), hkl.w() ); }
    //! get reciprocal space data (internal use)
    const std::complex<ffttype>& cplx_data( const Coord_grid& hkl ) const
      { return row_kl( hkl.v(), hkl.w() )[ hkl.u() ]; }

    //! Transform to real space
    void fft_x_to_h( const ftype& scale );
  };


} // namespace clipper

#endif
