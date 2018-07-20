/*! \file lib/container_types.h
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


#ifndef CLIPPER_CONTAINER_TYPES
#define CLIPPER_CONTAINER_TYPES


#include "container.h"
#include "../imex.h"

namespace clipper
{

  CLIPPER_IMEX extern const Resolution NullResolution;  //<! null instance
  CLIPPER_IMEX extern const Spacegroup NullSpacegroup;  //<! null instance
  CLIPPER_IMEX extern const Cell NullCell;  //<! null instance
  CLIPPER_IMEX extern const Grid_sampling NullGrid_sampling;  //<! null instance


  //! Resolution container
  /*! CResolution: This has a name and a resolution.  It overrides the
    base resolution for any objects below it. */
  class CLIPPER_IMEX CResolution : public Container, public Resolution
  {
  public:
    //! constructor: make null object or top object in tree
    CResolution( const String name = "",
		 const Resolution& resol = NullResolution ) :
      Container( name ), Resolution( resol ) {}
    //! constructor: make child object
    CResolution( Container& parent, const String name = "",
		 const Resolution& resol = NullResolution ) :
      Container( parent, name ), Resolution( resol ) {}
    //! initialiser: from Resolution
    void init( const Resolution& resolution_ );
  };


  //! Spacegroup container
  /*! CSpacegroup: This has a name and a spacegroup.  It overrides the
    base spacegroup for any objects below it. */
  class CLIPPER_IMEX CSpacegroup : public Container, public Spacegroup
  {
  public:
    //! constructor: make null object or top object in tree
    CSpacegroup( const String name = "",
		 const Spacegroup& spgr = NullSpacegroup ) :
      Container( name ), Spacegroup( spgr ) {}
    //! constructor: make child object
    CSpacegroup( Container& parent, const String name = "",
		 const Spacegroup& spgr = NullSpacegroup ) :
      Container( parent, name ), Spacegroup( spgr ) {}
    //! initialiser: from Spacegroup
    void init( const Spacegroup& spacegroup_ );
  };


  //! CCell container
  /*! CCell: This has a name and a cell.  It overrides the base cell
    for any objects below it. */
  class CLIPPER_IMEX CCell : public Container, public Cell
  {
  public:
    //! constructor: make null object or top object in tree
    CCell( const String name = "",
	   const Cell& cell = NullCell ) :
      Container( name ), Cell( cell ) {}
    //! constructor: make child object
    CCell( Container& parent, const String name = "",
	   const Cell& cell = NullCell ) :
      Container( parent, name ), Cell( cell ) {}
    //! initialiser: from Cell
    void init( const Cell& cell_ );
  };


  //! CGrid_sampling container
  /*! CGrid_sampling: This has a name and a grid sampling It overrides
    the grid sampling for any objects below it. */
  class CLIPPER_IMEX CGrid_sampling : public Container, public Grid_sampling
  {
  public:
    //! constructor: make null object or top object in tree
    CGrid_sampling( const String name = "",
		    const Grid_sampling& grid = NullGrid_sampling );
    //! constructor: make child object
    CGrid_sampling( Container& parent, const String name = "",
		    const ftype rate = 1.5 );
    //! constructor: make child object with explicit value
    CGrid_sampling( Container& parent, const String name,
		    const Grid_sampling& grid );
    //! initialiser: from sampling rate
    void init( const Spacegroup& spacegroup, const Cell& cell,
	       const Resolution& resolution, const ftype rate );
    //! initialiser: from Grid_sampling
    void init( const Grid_sampling& grid_sampling_ );
    //! hierarchical update
    void update();
  private:
    ftype rate_;  //!< Store for later initialisation
  };


} // namespace clipper

#endif
