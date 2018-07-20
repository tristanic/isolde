/*! \file lib/nxmap.h
    Header file for non-crystal map
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


#ifndef CLIPPER_NXMAP
#define CLIPPER_NXMAP


#include "derivs.h"
#include "../imex.h"

namespace clipper
{

  //! NXmap_base: base for non-crystallographic map class
  /*!
    The non-crystallographic map class stores a map of arbitrary data
    type. Unlike an Xmap it is finite in extent and has no
    symmetry. An RT operator provides mapping onto an arbitrary
    orthogonal coordinate frame. Iterators provide efficient access to
    data.

    This base contains everything except the data, which is templated
    in the derived type clipper::NXmap<T>.
  */
  class CLIPPER_IMEX NXmap_base
  {
  public:
    //! test if object has been initialised
    bool is_null() const;

    //! return the grid dimensions for this map
    const Grid& grid() const { return grid_; }
    //! return the orthogonal-to-grid coordinate operator
    const RTop<>& operator_orth_grid() const { return rt_orth_grid; }
    //! return the grid-to-orthogonal coordinate operator
    const RTop<>& operator_grid_orth() const { return rt_grid_orth; }
    //! convert map coordinate to orthogonal
    /*! \param cm The grid coordinate to be converted.
      \return The equivalent orthogonal coordinate. */
    inline Coord_orth coord_orth( const Coord_map& cm ) const
      { return Coord_orth( rt_grid_orth*cm ); }
    //! convert orthogonal coordinate to map
    /*! \param co The orthogonal coordinate to be converted.
      \return The equivalent grid coordinate. */
    inline Coord_map coord_map( const Coord_orth& co ) const
      { return Coord_map ( rt_orth_grid*co ); }

    //! is the given coord available in the map?
    bool in_map( const Coord_grid& pos ) const { return grid_.in_grid( pos ); }
    //! is the given coord available in the map using the given interpolant?
    template<class I> bool in_map( const Coord_map& cm ) const;

    //! get multiplicity of a map grid point (always 1 for NXmap)
    int multiplicity( const Coord_grid& ) const { return 1; }

    //! Map reference base class
    /*! This is a reference to an Map. It forms a base class for
      index-like and coordinate-like Map references. If you write a
      method which will work with either, then specify this instead of
      either of the derived classed. \internal */
    class Map_reference_base
    {
    public:
      //! return the parent NXmap
      inline const NXmap_base& base_nxmap() const { return *map_; }
      //! Get the index into the map data array
      inline const int& index() const { return index_; }
      //! Check for end of map
      inline bool last() const { return ( index_ >= map_->grid_.size() ); }
    protected:
      //! pointer to map for which this Map_reference_index is defined
      const NXmap_base* map_;
      //! integer index into map data array
      int index_;
    };

    //! Map reference with index-like behaviour
    /*! This is a reference to a map coordinate. It behaves like a
      simple index into the map, but can be easily converted into a
      coordinate as and when required. It also implements methods for
      iterating through the map. It is very compact, but coord()
      involves some overhead.

      \note The following methods are inherited from
      Map_reference_base but are documented here for convenience:
      base_nxmap(), index(), last().
    */
    class Map_reference_index : public Map_reference_base
    {
    public:
      //! Null constructor
      Map_reference_index() {}
      //! Constructor: need parent map
      explicit Map_reference_index( const NXmap_base& map )
	{ map_ = &map; index_ = 0; }
      //! Constructor: need parent map and coord
      Map_reference_index( const NXmap_base& map, const Coord_grid& pos )
	{ map_ = &map; index_ = map_->grid_.index( pos ); }
      //! Get current grid coordinate
      inline Coord_grid coord() const
	{ return map_->grid_.deindex(index_); }
      //! Get current value of orthogonal coordinate
      inline const Coord_orth coord_orth() const
	{ return map_->coord_orth( coord().coord_map() ); }
      //! Set current value of coordinate - optimised for nearby coords
      inline Map_reference_index& set_coord( const Coord_grid& pos )
	{ index_ = map_->grid_.index( pos ); return *this; }
      //! Simple increment
      inline Map_reference_index& next() { index_++; return *this; }
      //! Index of neighbouring point
      /* Use for e.g. peak search. Valid for -1 <= du/dv/dw <= 1 only.
	 \param du/dv/dw Coordinate offset. \return Map index. */
      inline int index_offset(const int& du,const int& dv,const int& dw) const {
	return index_ + du*map_->du + dv*map_->dv + dw*map_->dw;
      }
      // inherited functions listed for documentation purposes
      //-- const NXmap_base& base_nxmap() const;
      //-- const int& index() const;
      //-- bool last() const;
    };

    //! Map reference with coordinate-like behaviour
    /*! This is a reference to a map coordinate. It behaves like a
      coordinate, but also stores the index of the corresponding point
      in the map. It also implements methods for iterating through the
      a map. Since the current coordinate is stored, coord() is
      fast. However it required 5 words of storage.

      \note The following methods are inherited from
      Map_reference_base but are documented here for convenience:
      base_nxmap(), index(), last().
    */
    class Map_reference_coord : public NXmap_base::Map_reference_base
    {
    public:
      //! Null constructor
      Map_reference_coord() {}
      //! Constructor: need parent map
      explicit Map_reference_coord( const NXmap_base& map )
	{ map_ = &map; }
      //! Constructor: need parent map and coord
      Map_reference_coord( const NXmap_base& map, const Coord_grid& pos )
	{ map_ = &map; set_coord( pos ); }
      //! Get current value of coordinate
      inline Coord_grid coord() const { return pos_; }
      //! Get current value of orthogonal coordinate
      inline const Coord_orth coord_orth() const
	{ return map_->coord_orth( coord().coord_map() ); }
      //! Set current value of coordinate - optimised for nearby coords
      inline Map_reference_coord& set_coord( const Coord_grid& pos )
	{ pos_ = pos; index_ = map_->grid_.index( pos_ ); return *this; }
      //! Simple increment
      /*! Use of this function resets the stored coordinate and sym */
      inline Map_reference_coord& next() {
	index_++;
	pos_ = map_->grid_.deindex(index_);
	return *this;
      }
      // Increment u,v,w
      inline Map_reference_coord& next_u() { pos_.u()++; index_ += map_->du; return *this; }  //!< increment u
      inline Map_reference_coord& next_v() { pos_.v()++; index_ += map_->dv; return *this; }  //!< increment v
      inline Map_reference_coord& next_w() { pos_.w()++; index_ += map_->dw; return *this; }  //!< increment w
      inline Map_reference_coord& prev_u() { pos_.u()--; index_ -= map_->du; return *this; }  //!< decrement u
      inline Map_reference_coord& prev_v() { pos_.v()--; index_ -= map_->dv; return *this; }  //!< decrement v
      inline Map_reference_coord& prev_w() { pos_.w()--; index_ -= map_->dw; return *this; }  //!< decrement w
      //! Assignment operator from a coord
      inline Map_reference_coord& operator =( const Coord_grid& pos )
	{ return set_coord( pos ); }
      // inherited functions listed for documentation purposes
      //-- const NXmap_base& base_nxmap() const;
      //-- const int& index() const;
      //-- bool last() const;
    protected:
      //! Current coord
      Coord_grid pos_;
    };

    //! return a basic Map_reference_index for this map
    Map_reference_index first() const { return Map_reference_index( *this ); }
    //! return a coord Map_reference_index for this map
    Map_reference_coord first_coord() const { return Map_reference_coord( *this ); }

  protected:
    Grid grid_;               //!< grid for the map
    RTop<> rt_orth_grid; //!< orth->grid operator
    RTop<> rt_grid_orth; //!< grid->orth operator
    int du, dv, dw;           //!< steps for shifts along u,v,w

    //! Null constructor, for later initialisation
    NXmap_base();
    //! initialiser: takes grid and orthogonal->grid coordinate operator 
    void init( const Grid& grid, const RTop<>& rt );
    //! initialiser: takes grid, cell, and fraction limits
    void init( const Cell& cell, const Grid_sampling& grid, const Grid_range& grid_extent );

    friend class NXmap_base::Map_reference_base;
    friend class NXmap_base::Map_reference_index;
    friend class NXmap_base::Map_reference_coord;
  };




  //! NXmap<T>: actual non-crystallographic map class
  /*!
    The non-crystallographic map class stores a map of arbitrary data
    type. Unlike an Xmap it is finite in extent and has no
    symmetry. An RT operator provides mapping onto an arbitrary
    orthogonal coordinate frame. Iterators provide efficient access to
    data.

    This is derived from NXmap_base, and adds the templatised data
    itself and the methods which deal with it.

    \note The following methods are inherited from NXmap_base but are
    documented here for convenience: grid(), coord_orth(),
    coord_grid(), first(), first_coord().
  */
  template<class T> class NXmap : public NXmap_base
  {
  public:
    //! Null constructor, for later initialisation
    NXmap() {}
    //! Constructor: takes grid and orthogonal->grid coordinate operator
    NXmap( const Grid& grid, const RTop<>& rt );
    //! Constructor: takes grid, cell, and extent
    NXmap( const Cell& cell, const Grid_sampling& grid, const Grid_range& grid_extent );
    //! initialiser: takes grid and orthogonal->grid coordinate operator 
    void init( const Grid& grid, const RTop<>& rt );
    //! initialiser: takes grid, cell, and fraction limits
    void init( const Cell& cell, const Grid_sampling& grid, const Grid_range& grid_extent );

    //! get data by Map_reference_index
    inline const T& operator[] (const NXmap_base::Map_reference_index i) const
      { return list[i.index()]; }
    //! set data by Map_reference_index
    inline T& operator[] (const NXmap_base::Map_reference_index i)
      { return list[i.index()]; }

    //! get data by Map_reference_coord
    inline const T& operator[] (const NXmap_base::Map_reference_coord i) const
      { return list[i.index()]; }
    //! set data by Map_reference_coord
    inline T& operator[] (const NXmap_base::Map_reference_coord i)
      { return list[i.index()]; }

    //! get a density value for an arbitrary position
    inline const T& get_data( const Coord_grid& pos ) const
      { return list[ grid_.index( pos ) ]; }
    //! set a density value for an arbitrary position
    inline void set_data( const Coord_grid& pos, const T& val )
      { list[ grid_.index( pos ) ] = val; }

    //! get map value for map coord using supplied interpolator
    template<class I> T interp( const Coord_map& pos ) const;
    //! get map value and grad for map coord using supplied interpolator
    template<class I> void interp_grad( const Coord_map& pos, T& val, Grad_map<T>& grad ) const;
    //! get map value and curv for map coord using supplied interpolator
    template<class I> void interp_curv( const Coord_map& pos, T& val, Grad_map<T>& grad, Curv_map<T>& curv ) const;

    // inherited functions listed for documentation purposes
    //-- const Grid& grid() const;
    //-- const RTop<> operator_orth_grid() const;
    //-- const RTop<> operator_grid_orth() const;
    //-- const Coord_orth coord_orth( const Coord_map&  cg ) const;
    //-- const Coord_map  coord_map ( const Coord_orth& co ) const;
    //-- const Map_reference_index first();
    //-- const Map_reference_coord first_coord();

    //! assignment operator: assigns a single value to the whole map
    const T& operator =( const T& value );
    //! add another map to this one
    const NXmap<T>& operator +=( const NXmap<T>& other );
    //! subtract another map from this one
    const NXmap<T>& operator -=( const NXmap<T>& other );

  private:
    std::vector<T> list;
  };



  // template fucntion definitions

  /*! Note that the higher the order of the interpolant, the more of
    the boundary of the map becomes inaccessible.
    \param cm The coord_map to test.
    \return true if interpolation can be performed at that coordinate. */
  template<class I> bool NXmap_base::in_map( const Coord_map& cm ) const
    { return I::can_interp( *this, cm ); }

  /*! Initialise an NXmap to some rhomboid chosen from within a crystal
    coordinate space, specified by the grid and a transformation from
    orthogonal to grid coordinates.
    \param grid The grid dimensions of the desired map.
    \param rt The rotation/transln op from orthogonal to grid coordinates. */
  template<class T> NXmap<T>::NXmap( const Grid& grid, const RTop<>& rt )
    { init( grid, rt ); }

  /*! Initialise an NXmap to some rhomboid chosen from within a crystal
    grid coordinate space, specified by a cell, sampling and box within
    that grid. This is useful for creating an NXmap which exactly
    matches some subregion of a crystallographic map.
    \param cell Unit cell defining the crystal space.
    \param grid The grid sampling of the given unit cell.
    \param grid_extent The map extent within that cell. */
  template<class T> NXmap<T>::NXmap( const Cell& cell, const Grid_sampling& grid, const Grid_range& grid_extent )
    { init( cell, grid, grid_extent ); }

  /*! Initialise an NXmap to some rhomboid chosen from within a crystal
    coordinate space, specified by the grid and a transformation from
    orthogonal to grid coordinates.
    \param grid The grid dimensions of the desired map.
    \param rt The rotation/transln op from orthogonal to grid coordinates. */
  template<class T> void NXmap<T>::init( const Grid& grid, const RTop<>& rt )
    { NXmap_base::init( grid, rt ); list.resize( grid.size() ); }

  /*! Initialise an NXmap to some rhomboid chosen from within a crystal
    grid coordinate space, specified by a cell, sampling and box within
    that grid. This is useful for creating an NXmap which exactly
    matches some subregion of a crystallographic map.
    \param cell Unit cell defining the crystal space.
    \param grid The grid sampling of the given unit cell.
    \param grid_extent The map extent within that cell. */
  template<class T> void NXmap<T>::init( const Cell& cell, const Grid_sampling& grid, const Grid_range& grid_extent )
    { NXmap_base::init( cell, grid, grid_extent ); list.resize( grid_extent.size() ); }


  /*! The value of the map at the desired non-grid map
    coordinate are calculated using
    the supplied interpolator template.
    \param pos The map coord at which the density is to be calcuated.
    \return The value of the density at that point.
    map coordinates (see Cell::coord_orth). */
  template<class T> template<class I> T NXmap<T>::interp( const Coord_map& pos ) const
  {
    T val;
    I::interp( *this, pos, val );
    return val;
  }


  /*! The value of the map at the desired non-grid map
    coordinate and its gradient are calculated using
    the supplied interpolator template.
    \param pos The map coord at which the density is to be calcuated.
    \param val The value of the density at that point.
    \param grad The interpolated gradient vector with respect to the
    map coordinates (see Cell::coord_orth).
    \param curv The interpolated curvature matrix with respect to the
    map coordinates (see Cell::coord_orth). */
  template<class T> template<class I> void NXmap<T>::interp_grad( const Coord_map& pos, T& val, Grad_map<T>& grad ) const
  {
    I::interp_grad( *this, pos, val, grad );
  }


  /*! The value of the map at the desired non-grid map
    coordinate and its gradient and curvature are calculated using
    the supplied interpolator template.
    \param pos The map coord at which the density is to be calcuated.
    \param val The value of the density at that point.
    \param grad The interpolated gradient vector with respect to the
    map coordinates (see Cell::coord_orth).
    \param curv The interpolated curvature matrix with respect to the
    map coordinates (see Cell::coord_orth). */
  template<class T> template<class I> void NXmap<T>::interp_curv( const Coord_map& pos, T& val, Grad_map<T>& grad, Curv_map<T>& curv ) const
  {
    I::interp_curv( *this, pos, val, grad, curv );
  }


  /*! All values, including missing values, are overwritten by the value.
    \param value The value to which the map is to be set. */
  template<class T> const T& NXmap<T>::operator =( const T& value )
  {
    // copy value into map
    Map_reference_index im;
    for ( im = first(); !im.last(); im.next() ) list[im.index()] = value;
    return value;
  }


  /*! The map grids must match. */
  template<class T> const NXmap<T>& NXmap<T>::operator +=( const NXmap<T>& other )
  {
    if ( grid() != other.grid() )
      Message::message( Message_fatal( "NXmap: map mismatch in +=" ) );
    Map_reference_index im;
    for ( im = first(); !im.last(); im.next() ) list[im.index()] += other[im];
    return (*this);
  }

  /*! The map grids must match. */
  template<class T> const NXmap<T>& NXmap<T>::operator -=( const NXmap<T>& other )
  {
    if ( grid() != other.grid() )
      Message::message( Message_fatal( "NXmap: map mismatch in -=" ) );
    Map_reference_index im;
    for ( im = first(); !im.last(); im.next() ) list[im.index()] -= other[im];
    return (*this);
  }


} // namespace clipper

#endif
