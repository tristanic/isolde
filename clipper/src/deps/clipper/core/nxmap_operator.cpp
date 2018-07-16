/* nxmap_operator.cpp: implementation file for non-crystal map operator */
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


#include "nxmap_operator.h"


namespace clipper {


/*! The object is not initialised, and will return is_null(). */
NX_operator::NX_operator()
{
  xfrac_nxgrid = RTop<>::null();
}

/*! The operator and inverse operator, together with any possible
  optimisations, are constructed to relate the give crystallographic
  and non-crystallographic grid frames, using the supplied orthogonal
  operator.
  \param xmap An Xmap defining the crystal grid frame.
  \param nxmap An NXmap defining the non-crystal grid frame.
  \param rtop The operator relating the orthogonal frame of the NXmap
  onto the orthogonal frame of the Xmap. */
NX_operator::NX_operator( const Xmap_base& xmap, const NXmap_base& nxmap, const RTop_orth& rtop )
{
  init( xmap, nxmap, rtop );
}

/*! The operator and inverse operator, together with any possible
  optimisations, are constructed to relate the give crystallographic
  and non-crystallographic grid frames, using the supplied orthogonal
  operator.
  \param cell The cell defining the crystal grid frame.
  \param grid The grid defining the crystal grid frame.
  \param nxmap An NXmap defining the non-crystal grid frame.
  \param rtop The operator relating the orthogonal frame of the NXmap
  onto the orthogonal frame of the Xmap. */
NX_operator::NX_operator( const Cell& cell, const Grid_sampling& grid, const NXmap_base& nxmap, const RTop_orth& rtop  )
{
  init( cell, grid, nxmap, rtop );
}

/*! The operator and inverse operator, together with any possible
  optimisations, are constructed to relate the give crystallographic
  and non-crystallographic grid frames, using the supplied orthogonal
  operator.
  \param xmap An Xmap defining the crystal grid frame.
  \param nxmap An NXmap defining the non-crystal grid frame.
  \param rtop The operator relating the orthogonal frame of the NXmap
  onto the orthogonal frame of the Xmap. */
void NX_operator::init( const Xmap_base& xmap, const NXmap_base& nxmap, const RTop_orth& rtop )
{
  init( xmap.cell(), xmap.grid_sampling(), nxmap, rtop );
}

/*! The operator and inverse operator, together with any possible
  optimisations, are constructed to relate the give crystallographic
  and non-crystallographic grid frames, using the supplied orthogonal
  operator.
  \param cell The cell defining the crystal grid frame.
  \param grid The grid defining the crystal grid frame.
  \param nxmap An NXmap defining the non-crystal grid frame.
  \param rtop The operator relating the orthogonal frame of the NXmap
  onto the orthogonal frame of the Xmap. */
void NX_operator::init( const Cell& cell, const Grid_sampling& grid, const NXmap_base& nxmap, const RTop_orth& rtop  )
{
  // make op to map nxmap grid -> cell frac coords
  nxgrid_xfrac = RTop_orth(cell.matrix_frac()) * rtop * nxmap.operator_grid_orth();
  // make op for grid coords
  nxgrid_xgrid = RTop_orth(grid.matrix_frac_grid()) * nxgrid_xfrac;
  // make inverse op
  xfrac_nxgrid = nxgrid_xfrac.inverse();
  xgrid_nxgrid = nxgrid_xgrid.inverse();

  // now make nearest int op to cell -> nxmap op
  RTop_orth xgrid_nxgrid_rnd, nxgrid_xgrid_rnd;
  for ( int j = 0; j < 3; j++ ) {
    xgrid_nxgrid_rnd.trn()[j] = rint( xgrid_nxgrid.trn()[j] );
    nxgrid_xgrid_rnd.trn()[j] = rint( nxgrid_xgrid.trn()[j] );
    xgrid_nxgrid_int.trn()[j] = int( xgrid_nxgrid_rnd.trn()[j] );
    nxgrid_xgrid_int.trn()[j] = int( nxgrid_xgrid_rnd.trn()[j] );
    for ( int i = 0; i < 3; i++ ) {
      xgrid_nxgrid_rnd.rot()(i,j) = rint( xgrid_nxgrid.rot()(i,j) );
      nxgrid_xgrid_rnd.rot()(i,j) = rint( nxgrid_xgrid.rot()(i,j) );
      xgrid_nxgrid_int.rot()(i,j) = int( xgrid_nxgrid_rnd.rot()(i,j) );
      nxgrid_xgrid_int.rot()(i,j) = int( nxgrid_xgrid_rnd.rot()(i,j) );
    }
  }

  // check for integer mapping
  x_nx_is_int = xgrid_nxgrid_rnd.equals( xgrid_nxgrid, 0.01 );
  x_nx_is_trn = x_nx_is_int &&
    xgrid_nxgrid_rnd.rot().equals( Mat33<>::identity(), 0.01 );
  nx_x_is_int = nxgrid_xgrid_rnd.equals( nxgrid_xgrid, 0.01 );
  nx_x_is_trn = nx_x_is_int &&
    nxgrid_xgrid_rnd.rot().equals( Mat33<>::identity(), 0.01 );
}


bool NX_operator::is_null() const
{
  return xfrac_nxgrid.is_null();
}


void NX_operator::debug() const
{
  std::cout << " X->NX is int? " << x_nx_is_int << "\n";
  std::cout << " X->NX is trn? " << x_nx_is_trn << "\n";
  std::cout << xgrid_nxgrid.format() << "\n";
  std::cout << " NX->X is int? " << x_nx_is_int << "\n";
  std::cout << " NX->X is trn? " << x_nx_is_trn << "\n";
  std::cout << nxgrid_xgrid.format() << "\n";
}


} // namespace clipper
