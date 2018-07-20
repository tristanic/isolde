/*! \file lib/container_map.h
    Header file for Container versions of various objects
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


#ifndef CLIPPER_CONTAINER_MAP
#define CLIPPER_CONTAINER_MAP


#include "container_types.h"
#include "nxmap_operator.h"


namespace clipper
{

  //! Crystallographic map container
  /*! CXmap: This is a crystallographic map. */
  template<class T> class CXmap : public Container, public Xmap<T>
  {
  public:
    //! null constructor
    CXmap() {}
    //! constructor: inherit spacegroup, cell and grid
    CXmap( Container& parent, const String name = "" );
    //! initialiser: supply or inherit spacegroup, cell and grid
    void init( const Spacegroup& spacegroup, const Cell& cell, const Grid_sampling& grid_sampling );
    //! hierarchical update
    void update();
  };


  //! Non-Crystallographic map container
  /*! CNXmap: This is a non-crystallographic map. Since it does not
    exist in a crystallographic frame, it does not inherit
    anything. */
  template<class T> class CNXmap : public Container, public NXmap<T>
  {
  public:
    //! null constructor
    CNXmap() {}
    //! constructor:
    CNXmap( Container& parent, const String name = "" ) :
      Container( parent, name ) {}
  };


  //! Non-Crystallographic map operator container
  /*! CNXmap: This is an operator relating a non-crystallographic map
    into a crystallgraphic frame. It can inherit the crystallographic
    cell and grid sampling. */
  template<class T> class CNXmap_operator : public Container, public NXmap_operator<T>
  {
  public:
    //! null constructor
    CNXmap_operator() {}
    //! constructor: do not initialise
    CNXmap_operator( Container& parent, const String name = "" );
    //! constructor: inherit cell and grid
    CNXmap_operator( Container& parent, const String name, const NXmap<T>& nxmap, const RTop_orth& nxop = RTop_orth(RTop<>::null()) );
    //! initialier: supply or inherit cell, grid, NXmap, RTop_orth
    void init( const Cell& cell, const Grid_sampling& grid, const NXmap<T>& nxmap, const RTop_orth& nxop );
    //! hierarchical update
    void update();
  private:
    const NXmap<T>* nxmap_;  //!< Store for later initialisation
    RTop_orth nxop_;   //!< Store for later initialisation
  };


  // template implementations


  /*! The object is constructed at the given location in the hierarchy.
    An attempt is made to initialise the object using information from
    its parents in the hierarchy.
    \param parent An object in the hierarchy (usually the parent of the
    new object).
    \param name The path from \c parent to the new object (usually just
    the name of the new object). */
  template<class T> CXmap<T>::CXmap( Container& parent, const String name ) : Container( parent, name )
  {
    init( NullSpacegroup, NullCell, NullGrid_sampling );
  }


  /*! An attempt is made to initialise the object using information from
    the supplied parameters, or if they are Null, from its parents in
    the hierarchy.
    \param spacegroup The spacegroup for the map.
    \param name The cell for the map.
    \param grid The grid for the map. */
  template<class T> void CXmap<T>::init( const Spacegroup& spacegroup, const Cell& cell, const Grid_sampling& grid_sampling )
  {
    // use supplied values by default
    const Spacegroup* sp = &spacegroup;  // use pointers so we can reassign
    const Cell* cp = &cell;
    const Grid_sampling* gp = &grid_sampling;
    // otherwise get them from the tree
    if ( sp->is_null() ) sp = parent_of_type_ptr<const Spacegroup>();
    if ( cp->is_null() ) cp = parent_of_type_ptr<const Cell>();
    if ( gp->is_null() ) gp = parent_of_type_ptr<const Grid_sampling>();
    // initialise
    if ( sp != NULL && cp != NULL && gp != NULL )
      if ( !sp->is_null() && !cp->is_null() && !gp->is_null() )
	Xmap<T>::init( *sp, *cp, *gp );
    Container::update();
  }


  /*! Hierarchical update. If this object is uninitialised, an attempt
    is made to initialise the object using information from its parents
    in the hierarchy. The childen of the object are then updated. */
  template<class T> void CXmap<T>::update()
  {
    if ( Xmap_base::is_null() )
      init( NullSpacegroup, NullCell, NullGrid_sampling );
    else
      Container::update();
  }


  /*! The object is not initialised.
    \param parent The objects parent.
    \param name The object name. */
  template<class T> CNXmap_operator<T>::CNXmap_operator( Container& parent, const String name ) : Container( parent, name ), nxmap_( NULL ), nxop_( RTop_orth::null() ) {}


  /*! The object is initialised if the appropriate parent objects are
    available, and children are updated.
    \param parent The objects parent.
    \param name The object name.
    \param nxmap The non-crystal map object.
    \param nxop The orth. operator mapping the NXmap into the crystal frame. */
  template<class T> CNXmap_operator<T>::CNXmap_operator( Container& parent, const String name, const NXmap<T>& nxmap, const RTop_orth& nxop ) : Container( parent, name ), nxmap_( &nxmap ), nxop_( nxop )
  {
    init( NullCell, NullGrid_sampling, NXmap<T>(), RTop_orth::null() );
  }


  /*! An attempt is made to initialise the object using information from
    the supplied parameters, or if they are Null, from its parents in
    the hierarchy.
    \param cell The unit cell for the crystallographic frame.
    \param grid The grid sampling for the crystallographic frame.
    \param nxmap The non-crystal map object.
    \param nxop The orth. operator mapping the NXmap into the crystal frame. */
  template<class T> void CNXmap_operator<T>::init( const Cell& cell, const Grid_sampling& grid_sampling, const NXmap<T>& nxmap, const RTop_orth& nxop )
  {
    // use supplied values by default
    const Cell* cp = &cell;
    const Grid_sampling* gp = &grid_sampling;
    // otherwise get them from the tree
    if ( cp->is_null() ) cp = parent_of_type_ptr<const Cell>();
    if ( gp->is_null() ) gp = parent_of_type_ptr<const Grid_sampling>();
    if ( !nxmap.is_null() ) nxmap_ = &nxmap;
    if ( !nxop.is_null()  ) nxop_ = nxop;
    // initialise
    if ( cp != NULL && gp != NULL && nxmap_ != NULL )
      if ( !cp->is_null() && !gp->is_null() && !nxmap_->is_null() )
	NXmap_operator<T>::init( *cp, *gp, *nxmap_, nxop_ );
    Container::update();
  }


  /*! Hierarchical update. If this object is uninitialised, an attempt
    is made to initialise the object using information from its parents
    in the hierarchy. The childen of the object are then updated. */
  template<class T> void CNXmap_operator<T>::update()
  {
    if ( NX_operator::is_null() )
      init( NullCell, NullGrid_sampling, NXmap<T>(), RTop_orth::null() );
    else
      Container::update();
  }


} // namespace clipper

#endif
