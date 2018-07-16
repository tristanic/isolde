/* nxmap.cpp: implementation file for crystal maps */
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


#include "nxmap.h"


namespace clipper {


Message_ctor message_ctor_nxmap( " [NXmap: constructed]" );


/*! For later initialisation: see init() */
NXmap_base::NXmap_base()
{
  rt_orth_grid = RTop_orth( RTop<>::null() );
  Message::message( message_ctor_nxmap );
}


/*! Initialise an NXmap to some rhomboid chosen from within a crystal
  coordinate space, specified by the grid and a transformation from
  orthogonal to grid coordinates.
  \param grid The grid dimensions of the desired map.
  \param rt The rotation translation op from orthogonal to grid coordinates. */
void NXmap_base::init( const Grid& grid, const RTop<>& rt )
{
  // set up grid and orth->grid operator
  grid_ = grid;
  rt_orth_grid = rt;
  rt_grid_orth = rt.inverse();

  // set up grid steps
  du = grid_.index( Coord_grid(1,0,0) );
  dv = grid_.index( Coord_grid(0,1,0) );
  dw = grid_.index( Coord_grid(0,0,1) );
}

/*! Initialise an NXmap to some rhomboid chosen from within a crystal
  grid coordinate space, specified by a cell, sampling and box within
  that grid. This is useful for creating an NXmap which exactly
  matches some subregion of a crystallographic map.
  \param cell Unit cell defining the crystal space.
  \param grid The grid sampling of the given unit cell.
  \param grid_extent The map extent within that cell. */
void NXmap_base::init( const Cell& cell, const Grid_sampling& grid, const Grid_range& grid_extent )
{
  // initialise
  init( grid_extent,
	RTop<>( grid.matrix_frac_grid() * cell.matrix_frac(),
		     -Coord_map( grid_extent.min() ) ) );
}

/*! \return true if the object has not been initalised. */
bool NXmap_base::is_null() const
{ return rt_orth_grid.is_null(); }


} // namespace clipper
