/* container_minimol.cpp: atomic model types */
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

#include "container_minimol.h"


namespace clipper {


// CMiniMol

/*! The object is constructed at the given location in the hierarchy.
  An attempt is made to initialise the object using information from
  its parents in the hierarchy.
  \param parent An object in the hierarchy (usually the parent of the
  new object).
  \param name The path from \c parent to the new object (usually just
  the name of the new object). */
CMiniMol::CMiniMol( Container& parent, const String name ) : Container( parent, name )
{
  init( NullSpacegroup, NullCell );
}

/*! An attempt is made to initialise the object using information from
  the supplied parameters, or if they are Null, from its parents in
  the hierarchy.
  \param spacegroup The spacegroup for the model.
  \param name The cell for the model. */
void CMiniMol::init( const Spacegroup& spacegroup, const Cell& cell )
{
  // use supplied values by default
  const Spacegroup* sp = &spacegroup;  // use pointers so we can reassign
  const Cell* cp = &cell;
  // otherwise get them from the tree
  if ( sp->is_null() ) sp = parent_of_type_ptr<const Spacegroup>();
  if ( cp->is_null() ) cp = parent_of_type_ptr<const Cell>();
  // initialise
  if ( sp != NULL && cp != NULL )
    if ( !sp->is_null() && !cp->is_null() )
      MiniMol::init( *sp, *cp );
  Container::update();
}

/*! Hierarchical update. If this object is uninitialised, an attempt
  is made to initialise the object using information from its parents
  in the hierarchy. The childen of the object are then updated. */
void CMiniMol::update()
{
  if ( CMiniMol::is_null() )
    init( NullSpacegroup, NullCell );
  else
    Container::update();
}


} // namespace clipper
