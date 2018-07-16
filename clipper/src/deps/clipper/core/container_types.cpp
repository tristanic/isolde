/* container_types.cpp: class file for basic containers */
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

#include "container_types.h"


namespace clipper {


const Resolution NullResolution;  //<! null instance
const Spacegroup NullSpacegroup;  //<! null instance
const Cell NullCell;  //<! null instance
const Grid_sampling NullGrid_sampling;  //<! null instance


/*! The object is initialised, and children are updated.
  \param resolution_ The value to give to the contained object. */
void CResolution::init( const Resolution& resolution_ )
{
  dynamic_cast<Resolution&>(*this) = resolution_;
  Container::update();
}


/*! The object is initialised, and children are updated.
  \param spacegroup_ The value to give to the contained object. */
void CSpacegroup::init( const Spacegroup& spacegroup_ )
{
  dynamic_cast<Spacegroup&>(*this) = spacegroup_;
  Container::update();
}


/*! The object is initialised, and children are updated.
  \param cell_ The value to give to the contained object. */
void CCell::init( const Cell& cell_ )
{
  dynamic_cast<Cell&>(*this) = cell_;
  Container::update();
}


/*! The top object in a tree is initialised from a known grid.
  \param name The object name.
  \param grid The grid sampling. */
CGrid_sampling::CGrid_sampling( const String name, const Grid_sampling& grid ) : Container( name ), Grid_sampling( grid ), rate_( 1.5 ) {}


/*! The normal form for a child object - spacegroup and cell inherited.
  \param parent The objects parent.
  \param name The object name.
  \param rate The Shannon rate (default 1.5). */
CGrid_sampling::CGrid_sampling( Container& parent, const String name, const ftype rate ) : Container( parent, name ), rate_( rate )
{
  init( NullSpacegroup, NullCell, NullResolution, 0.0 );
}


/*! This is still a child object but is initialised directly.
  \param parent The objects parent.
  \param name The object name.
  \param grid The grid sampling. */
CGrid_sampling::CGrid_sampling( Container& parent, const String name, const Grid_sampling& grid ) : Container( parent, name ), Grid_sampling( grid ), rate_( 1.5 ) {}


/*! The object is initialised if the appropriate parent objects are
  available, and children are updated.
  \param spacegroup The spacegroup.
  \param cell The cell.
  \param resolution The resolution.
  \param rate_ The Shannon rate (If <1 previous value is used, default 1.5). */
void CGrid_sampling::init( const Spacegroup& spacegroup, const Cell& cell, const Resolution& resolution, const ftype rate )
{
  // use supplied values by default
  const Spacegroup* sp = &spacegroup;  // use pointers so we can reassign
  const Cell*       cp = &cell;
  const Resolution* rp = &resolution;
  // otherwise get them from the tree
  if ( sp->is_null() ) sp = parent_of_type_ptr<const Spacegroup>();
  if ( cp->is_null() ) cp = parent_of_type_ptr<const Cell>();
  if ( rp->is_null() ) rp = parent_of_type_ptr<const Resolution>();
  if ( rate >= 1.0 ) rate_ = rate;  // use rate if given
  // initialise
  if ( sp != NULL && cp != NULL && rp != NULL )
    if ( !sp->is_null() && !cp->is_null() && !rp->is_null() )
      Grid_sampling::init( *sp, *cp, *rp, rate_ );
  Container::update();
}


/*! The object is initialised, and children are updated.
  \param grid_sampling_ The value to give to the contained object. */
void CGrid_sampling::init( const Grid_sampling& grid_sampling_ )
{
  dynamic_cast<Grid_sampling&>(*this) = grid_sampling_;
  Container::update();
}


/*! Hierarchical update. If this object is uninitialised, an attempt
  is made to initialise the object using information from its parents
  in the hierarchy. The childen of the object are then updated. */
void CGrid_sampling::update()
{
  if ( Grid_sampling::is_null() )
    init( NullSpacegroup, NullCell, NullResolution, 0.0 );
  else
    Container::update();
}


} // namespace clipper
