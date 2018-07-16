/*! \file lib/nxmap_operator.h
    Header file for non-crystal map operator
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


#ifndef CLIPPER_NXMAP_OPERATOR
#define CLIPPER_NXMAP_OPERATOR

#include "xmap.h"
#include "nxmap.h"
#include "../imex.h"

namespace clipper
{

  //! NX_operator: non-crystal map operator
  /*! This class holds a reference to a non-crystal map frame from
    somewhere within a crystallographic map frame. In the general
    case, an orthogonal rotation-translation operator is provided
    which maps the orthogonal frame of the crystal space onto the
    orthogonal frame of the NXmap space.

    The object calculates and stores optimised transformations between
    the crystallgoraphic frame (described either in fractional or grid
    coordinates), and the NXmap grid. Fast paths are generated
    automatically if the grids are related.
  */
  class CLIPPER_IMEX NX_operator
  {
  public:
    //! null constructor
    NX_operator();
    //! constructor: from Xmap, NXmap, and operator
    NX_operator( const Xmap_base& xmap, const NXmap_base& nxmap, const RTop_orth& rtop );
    //! constructor: from cell, grid sampling, NXmap, and operator
    NX_operator( const Cell& cell, const Grid_sampling& grid, const NXmap_base& nxmap, const RTop_orth& rtop  );
    //! initialiser:: from Xmap, NXmap, and operator
    void init( const Xmap_base& xmap, const NXmap_base& nxmap, const RTop_orth& rtop  );
    //! initialiser:: from cell, grid sampling, NXmap, and operator
    void init( const Cell& cell, const Grid_sampling& grid, const NXmap_base& nxmap, const RTop_orth& rtop  );

    //! convert xtal frac coord to nxmap map coord
    inline Coord_map coord_map( const Coord_frac& c ) const
      { return Coord_map( xfrac_nxgrid * c ); }
    //! convert nxmap map coord to xtal frac coord
    inline Coord_frac coord_frac( const Coord_map& c ) const
      { return Coord_frac( nxgrid_xfrac * c ); }
    //! get value of nxmap at xmap grid coord using fastest appropriate method
    template<class I, class T, class M> T nxmap_data( const M& nxmap, const Coord_grid& c ) const;
    //! get value of xmap at nxmap grid coord using fastest appropriate method
    template<class I, class T, class M> T xmap_data( const M& xmap, const Coord_grid& c ) const;

    //! test if object has been initialised
    bool is_null() const;

    void debug() const;

  protected:
    RTop<> xfrac_nxgrid;    //!< xtal_cell -> nxmap operator
    RTop<> nxgrid_xfrac;    //!< nxmap -> xtal_cell operator
    RTop<> xgrid_nxgrid;    //!< xtal_grid -> nxmap operator
    RTop<> nxgrid_xgrid;    //!< nxmap -> xtal_grid operator
    RTop<int> xgrid_nxgrid_int;  //!< xtal_grid -> nxmap integer operator
    RTop<int> nxgrid_xgrid_int;  //!< nxmap -> xtal_grid integer operator
    bool x_nx_is_int;      //!< true if int operator exists
    bool x_nx_is_trn;      //!< true if int operator exists and is pure transln
    bool nx_x_is_int;      //!< true if int operator exists
    bool nx_x_is_trn;      //!< true if int operator exists and is pure transln
  };


  //! NXmap_operator: non-crystal map operator referencing a particular NXmap
  /*! This class holds a reference to a non-crystal map object from
    somewhere within a crystallographic map frame. In the general
    case, an orthogonal rotation-translation operator is provided
    which maps the orthogonal frame of the crystal space onto the
    orthogonal frame of the NXmap space.

    The object calculates and stores optimised transformations between
    the crystallgoraphic frame (described either in fractional or grid
    coordinates), and the NXmap grid. Fast paths are generated
    automatically if the grids are related.

    \note This object differes from NX_operator in that it keeps a
    reference to an individual NXmap, which may be used to access that
    object directly.
  */
  template<class T> class NXmap_operator : public NX_operator
  {
  public:
    //! null constructor
    NXmap_operator() {}
    //! constructor: from Xmap, NXmap, and operator
    NXmap_operator( const Xmap_base& xmap, const NXmap<T>& nxmap, const RTop_orth& rtop ) { init( xmap, nxmap, rtop ); }
    //! constructor: from cell, grid sampling, NXmap, and operator
    NXmap_operator( const Cell& cell, const Grid_sampling& grid, const NXmap<T>& nxmap, const RTop_orth& rtop  ) { init( cell, grid, nxmap, rtop ); }
    //! initialiser:: from Xmap, NXmap, and operator
    void init( const Xmap_base& xmap, const NXmap<T>& nxmap, const RTop_orth& rtop  ) { init( xmap.cell(), xmap.grid_sampling(), nxmap, rtop ); }
    //! initialiser:: from cell, grid sampling, NXmap, and operator
    void init( const Cell& cell, const Grid_sampling& grid, const NXmap<T>& nxmap, const RTop_orth& rtop  ) { nxmap_ = &nxmap; NX_operator::init( cell, grid, nxmap, rtop ); }

    //! access NXmap directly from xmap grid coord using fastest method
    template<class I> T nxmap_data( const Coord_grid& c ) const
      { return NX_operator::nxmap_data<I,T>( *nxmap_, c ); }

    //! get the target NXmap of this operator
    const NXmap<T>& nxmap() const { return *nxmap_; }

  private:
    const NXmap<T>* nxmap_;  //!< pointer to the nxmap
  };



  // template implementations

  /*! The density of the non-crystal map at the position corresponding to
    a crystallographic map grid coordinate is returned. If the grids
    match exactly either by pure translation or by rotation+translation,
    then fast paths are used to return the requested density
    directly. Otherwise the supplied interpolation template is used.
    No checking is performed for coordinates outside the NXmap.
    \param nxmap The non-crystal map (NXmap) to be queried.
    \param c The grid coordinate in the crystallographic coordinate frame.
    \return The value of the NXmap at the requested position. */
  template<class I, class T, class M> T NX_operator::nxmap_data( const M& nxmap, const Coord_grid& c ) const
  {
    if ( x_nx_is_trn ) {
      return T( nxmap.get_data( Coord_grid( c + xgrid_nxgrid_int.trn() ) ) );
    } else if ( x_nx_is_int ) {
      return T( nxmap.get_data( Coord_grid( xgrid_nxgrid_int * c ) ) );
    } else {
      T val;
      I::interp( nxmap, Coord_map( xgrid_nxgrid * c.coord_map() ), val );
      return val;
    }
  }


  /*! The density of the crystal map at the position corresponding to
    a non-crystallographic map grid coordinate is returned. If the grids
    match exactly either by pure translation or by rotation+translation,
    then fast paths are used to return the requested density
    directly. Otherwise the supplied interpolation template is used.
    \param xmap The crystal map (Xmap) to be queried.
    \param c The grid coordinate in the crystallographic coordinate frame.
    \return The value of the Xmap at the requested position. */
  template<class I, class T, class M> T NX_operator::xmap_data( const M& xmap, const Coord_grid& c ) const
  {
    if ( nx_x_is_trn ) {
      return T( xmap.get_data( Coord_grid( c + nxgrid_xgrid_int.trn() ) ) );
    } else if ( nx_x_is_int ) {
      return T( xmap.get_data( Coord_grid( nxgrid_xgrid_int * c ) ) );
    } else {
      T val;
      I::interp( xmap, Coord_map( nxgrid_xgrid * c.coord_map() ), val );
      return val;
    }
  }


} // namespace clipper

#endif
