/*! \file lib/container_hkl.h
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


#ifndef CLIPPER_CONTAINER_HKL
#define CLIPPER_CONTAINER_HKL


#include "container_types.h"
#include "hkl_data.h"
#include "../imex.h"


namespace clipper
{

  extern CLIPPER_IMEX const HKL_info NullHKL_info;  //<! null instance


  //! HKL list and indexing object container
  /*! CHKL_info: This is the reflection list object for reflection
    data objects to reference in order to identify their data
    entries. */
  class CLIPPER_IMEX CHKL_info : public Container, public HKL_info
  {
  public:
    //! constructor: make null object or top object in tree
    CHKL_info( const String name = "",
	       const Spacegroup& spacegroup = NullSpacegroup,
	       const Cell& cell = NullCell,
	       const Resolution& resolution = NullResolution,
	       const bool& generate = false );
    //! constructor: inherit spacegroup, cell and resolution
    CHKL_info( Container& parent, const String name = "",
	       const bool& generate = false );
    //! initialiser: supply or inherit spacegroup, cell and resolution
    void init( const Spacegroup& spacegroup, const Cell& cell,
	       const Resolution& resolution, const bool& generate = false );

    //! synthesize hkl list and update children
    void generate_hkl_list();    
    //! hierarchical update
    void update();
  private:
    bool generate_;  //!< Store for later initialisation
  };


  //! Reflection data list container
  /*! CHKL_data: This is the list object containing the actual data.
    It must be indexed by a parent HKL list. */
  template<class T> class CHKL_data : public Container, public HKL_data<T>
  {
  public:
    //! null constructor
    CHKL_data() {}
    //! constructor: inherit datalist and cell
    CHKL_data( Container& parent, const String name = "" );
    //! initialiser: supply or inherit hkl list, and cell
    void init( const HKL_info& hkl_info, const Cell& cell );
    //! initialiser: from spacegroup, cell, and HKL_sampling
    void init( const Spacegroup& spacegroup, const Cell& cell, const HKL_sampling& hkl_sampling ) { clipper::Message::message( clipper::Message_fatal("CHKL_data deprecated") ); }
    //! hierarchical update
    void update();

    //! assignment operator: copies the data from another list
    HKL_data<T>& operator =( const HKL_data<T>& other )
      { return ( dynamic_cast<HKL_data<T>&>(*this) = other ); }
    //! assignment operator: assigns a single value to the whole list
    HKL_data<T>& operator =( const T& value )
      { return ( dynamic_cast<HKL_data<T>&>(*this) = value ); }
  };


  // template implementations


  /*! The object is constructed at the given location in the hierarchy.
    An attempt is made to initialise the object using information from
    its parents in the hierarchy.
    \param parent An object in the hierarchy (usually the parent of the
    new object).
    \param name The path from \c parent to the new object (usually just
    the name of the new object). */
  template<class T> CHKL_data<T>::CHKL_data( Container& parent, const String name ) : Container( parent, name )
  {
    init( NullHKL_info, NullCell );
  }


  /*! An attempt is made to initialise the object using information from
    the supplied parameters, or if they are Null, from its parents in
    the hierarchy.
    \param hkl_info The reflection list object for this datalist.
    \param cell The cell object for this datalist. */
  template<class T> void CHKL_data<T>::init( const HKL_info& hkl_info, const Cell& cell )
  {
    // use supplied values by default
    const HKL_info* hp = &hkl_info;
    const Cell* cp = &cell;
    // otherwise get them from the tree
    if ( hp->is_null() ) hp = parent_of_type_ptr<const HKL_info>();
    if ( cp->is_null() ) cp = parent_of_type_ptr<const Cell>();
    // if we have an hkl_info, try and init
    if ( hp != NULL ) {
      // if no cell then inherit
      if ( cp == NULL ) cp = &(hp->cell());
      // now initialise
      if ( !hp->is_null() && !cp->is_null() )
	HKL_data<T>::init( *hp, *cp );
    }
    Container::update();
  }


  /*! Hierarchical update. If this object is uninitialised, an attempt
    is made to initialise the object using information from its parents
    in the hierarchy. The childen of the object are then updated.

    The data list is also synchronized with the parent reflection list. */
  template<class T> void CHKL_data<T>::update()
  {
    if ( HKL_data_base::is_null() )
      init( NullHKL_info, NullCell );
    else
      Container::update();
    HKL_data<T>::update();
  }


} // namespace clipper

#endif
