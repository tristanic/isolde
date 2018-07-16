/*! \file lib/xmap.h
    Header file for crystal maps
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


#ifndef CLIPPER_XMAP
#define CLIPPER_XMAP


#include "fftmap.h"
#include "fftmap_sparse.h"
#include "derivs.h"
#include "../imex.h"

namespace clipper
{

  class CLIPPER_IMEX Xmap_cacheobj
  {
  public:
    class Key
    {
    public:
      Key( const Spgr_descr& spgr_descr, const Grid_sampling& grid ) :
	spgr_descr_(spgr_descr), grid_sampling_(grid) {}
      const Spgr_descr& spgr_descr() const { return spgr_descr_; }
      const Grid_sampling& grid_sampling() const { return grid_sampling_; }
    private:
      Spgr_descr spgr_descr_;
      Grid_sampling grid_sampling_;
    };

    Xmap_cacheobj( const Key& xmap_cachekey );  //!< construct entry
    bool matches( const Key& xmap_cachekey ) const; //!< compare entry
    String format() const;  //!< string description
    // data
    Key key;                         //!< key
    Grid_sampling xtl_grid;          //!< grid for the cell
    Grid_range asu_grid;             //!< grid for the ASU
    Grid_range map_grid;             //!< grid for the ASU, plus border
    int nsym;                        // number of symops
    std::vector<unsigned char> asu;  //!< ASU flag array
    std::vector<Isymop> isymop;      //!< Integerised symops
    std::vector<int> du, dv, dw;     //!< symmetry grid shifts to index
    Array2d<unsigned char> symperm;  //!< Perumtation matrix of symops
    Mat33<> mat_grid_orth;           //!< for backward compatibility
    static Mutex mutex;              //!< thread safety
private:
    void find_asu_sym_(const Coord_grid& begin, const Coord_grid& end);
    void find_map_sym_(const Coord_grid& begin, const Coord_grid& end);
    int loops_per_thread_ = 2.5e6;

  };


  //! Xmap_base: base for crystallographic map class
  /*!
    The crystallographic map class stores a map of arbitrary data
    type. Its main difference from a 3-d array is that the data extent
    appears to be infinite, and yet internally only a unique ASU is
    stored. Iterators provide efficient access to data.

    This base contains everything except the data, which is templated
    in the derived type Xmap<T>
  */
  class CLIPPER_IMEX Xmap_base
  {
  public:
    enum FFTtype { Default, Normal, Sparse };  //!< FFT backend selection

    //! test if object has been initialised
    bool is_null() const;

    //! get the cell
    const Cell& cell() const { return cell_; }
    //! get the spacegroup
    const Spacegroup& spacegroup() const { return spacegroup_; }
    //! get the cell grid
    const Grid_sampling& grid_sampling() const { return grid_sam_; }
    //! get the ASU grid
    const Grid_range& grid_asu() const { return cacheref.data().asu_grid; }
    //! map coordinate from index
    /*! \param index The index. \return The corresponding grid coordinate. */
    inline Coord_grid coord_of( const int& index ) const
      { return cacheref.data().map_grid.deindex( index ); }
    //! map index from coordinate
    /*! This does not check symmetry equivalents. \param coord The coordinate.
      \return The index, or -1 if it does not exist. */
    inline int index_of( const Coord_grid& coord ) const {
      if ( cacheref.data().asu_grid.in_grid( coord ) ) {
	const int i = cacheref.data().map_grid.index( coord );
	if ( asu[ i ] == 0 ) return i;
      }
      return -1;
    }
    //! function to pick right cell repeat for any grid coord
    Coord_grid to_map_unit( const Coord_grid& pos ) const
      { return pos.unit( grid_sam_ ); }

    //! return the orthogonal-to-grid coordinate operator (translation is zero)
    const RTop<>& operator_orth_grid() const { return rt_orth_grid; }
    //! return the grid-to-orthogonal coordinate operator (translation is zero)
    const RTop<>& operator_grid_orth() const { return rt_grid_orth; }
    //! convert map coordinate to orthogonal
    /*! \param cm The grid coordinate to be converted.
      \return The equivalent orthogonal coordinate. */
    inline Coord_orth coord_orth( const Coord_map& cm ) const
      { return Coord_orth( rt_grid_orth.rot()*cm ); }
    //! convert orthogonal coordinate to map
    /*! \param co The orthogonal coordinate to be converted.
      \return The equivalent grid coordinate. */
    inline Coord_map coord_map( const Coord_orth& co ) const
      { return Coord_map ( rt_orth_grid.rot()*co ); }

    //! (This method is for compatibility with NXmap - it always returns true)
    bool in_map( const Coord_grid& ) const { return true; }
    //! (This method is for compatibility with NXmap - it always returns true)
    template<class I> bool in_map( const Coord_map& cm ) const { return true; }

    //! get multiplicity of a map grid point
    int multiplicity( const Coord_grid& pos ) const;

    //! Map reference base class
    /*! This is a reference to an Map. It forms a base class for
      index-like and coordinate-like Map references. If you write a
      method which will work with either, then specify this instead of
      either of the derived classed. \internal */
    class CLIPPER_IMEX Map_reference_base
    {
    public:
      //! return the parent Xmap
      inline const Xmap_base& base_xmap() const { return *map_; }
      //! Get the index into the map data array
      inline const int& index() const { return index_; }
      //! Check for end of map
      bool last() const { return ( index_ >= map_->map_grid.size() ); }
    protected:
      //! pointer to map for which this Map_reference_index is defined
      const Xmap_base* map_;
      //! integer index_ into map data array
      int index_;
    };

    //! Map reference with index-like behaviour
    /*! This is a reference to a map coordinate. It behaves like a
      simple index into the map, but can be easily converted into a
      coordinate as and when required. It also implements methods for
      iterating through the unique portion of a map. It is very
      compact, but coord() involves some overhead and loses any
      information concerning symmetry equivelents.

      \note The following methods are inherited from
      Map_reference_base but are documented here for convenience:
      base_xmap(), index(), last().
    */
    class CLIPPER_IMEX Map_reference_index : public Map_reference_base
    {
    public:
      //! Null constructor
      Map_reference_index() {}
      //! Constructor: takes parent map
      explicit Map_reference_index( const Xmap_base& map )
	{ map_ = &map; index_=0; next(); }
      //! Constructor: takes parent map and coord
      Map_reference_index( const Xmap_base& map, const Coord_grid& pos ) { map_ = &map; int sym; map_->find_sym( pos, index_, sym ); }
      //! Get current grid coordinate
      inline Coord_grid coord() const
	{ return map_->map_grid.deindex(index_); }
      //! Get current value of orthogonal coordinate
      inline const Coord_orth coord_orth() const
	{ return Coord_orth( map_->rt_grid_orth.rot() * coord().coord_map() ); }
      //! Set current value of coordinate - optimised for nearby coords
      inline Map_reference_index& set_coord( const Coord_grid& pos )
	{ int sym; map_->find_sym( pos, index_, sym ); return *this; }
      //! Simple increment
      inline Map_reference_index& next() {
	do {
	  index_++; if ( last() ) break;
	} while ( map_->asu[index_] != 0 );
	return *this;
      }
      //! Index of neighbouring point
      /* Use for e.g. peak search. Valid for -1 <= du/dv/dw <= 1 only.
	 \param du/dv/dw Coordinate offset. \return Map index. */
      inline int index_offset(const int& du,const int& dv,const int& dw) const {
	int i = index_ + du*map_->du[0] + dv*map_->dv[0] + dw*map_->dw[0];
	if ( map_->asu[i] != 0 ) { i = map_->map_grid.index( map_->to_map_unit( map_->map_grid.deindex(i).transform( map_->isymop[map_->asu[i]-1] ) ) ); }
	return i;
      }
      // inherited functions listed for documentation purposes
      //-- const Xmap_base& base_xmap() const;
      //-- const int& index() const;
      //-- bool last() const;
    };

    //! Map reference with coordinate-like behaviour
    /*! This is a reference to a map coordinate. It behaves like a
      coordinate, but also stores the index of the corresponding point
      in the map, and the symmetry operator required to get there. It
      also implements methods for iterating through the a map. Since
      the current coordinate and symmetry are stored, coord() is
      fast. However, it requires 1 pointer and 5 words of storage.

      \note The following methods are inherited from
      Map_reference_base but are documented here for convenience:
      base_xmap(), index(), last().
    */
    class CLIPPER_IMEX Map_reference_coord : public Map_reference_base
    {
    public:
      //! Null constructor
      Map_reference_coord() {}
      //! Constructor: takes parent map
      explicit Map_reference_coord( const Xmap_base& map )
	{ map_ = &map; index_ = 0; next(); }
      //! Constructor: takes parent map and coord
      Map_reference_coord( const Xmap_base& map, const Coord_grid& pos ) {
	map_ = &map;
	pos_ = pos;
	map_->find_sym( pos_, index_, sym_ );
      }
      //! Get current value of coordinate
      inline const Coord_grid& coord() const { return pos_; }
      //! Get current value of orthogonal coordinate
      inline const Coord_orth coord_orth() const
	{ return Coord_orth( map_->rt_grid_orth.rot() * coord().coord_map() ); }
      //! Get current symmetry operator
      inline const int& sym() const { return sym_; }
      //! Set current value of coordinate - optimised for nearby coords
      Map_reference_coord& set_coord( const Coord_grid& pos );
      //! Simple increment
      /*! Use of this function resets the stored coordinate and sym */
      inline Map_reference_coord& next() {
	sym_ = 0;
	do {
	  index_++; if ( last() ) break;
	} while ( map_->asu[index_] != 0 );
	pos_ = map_->map_grid.deindex(index_);
	return *this;
      }
      // Increment u,v,w
      inline Map_reference_coord& next_u() { pos_.u()++; index_ += map_->du[sym_]; if (map_->asu[index_] != 0) edge(); return *this; }  //!< increment u
      inline Map_reference_coord& next_v() { pos_.v()++; index_ += map_->dv[sym_]; if (map_->asu[index_] != 0) edge(); return *this; }  //!< increment v
      inline Map_reference_coord& next_w() { pos_.w()++; index_ += map_->dw[sym_]; if (map_->asu[index_] != 0) edge(); return *this; }  //!< increment w
      inline Map_reference_coord& prev_u() { pos_.u()--; index_ -= map_->du[sym_]; if (map_->asu[index_] != 0) edge(); return *this; }  //!< increment u
      inline Map_reference_coord& prev_v() { pos_.v()--; index_ -= map_->dv[sym_]; if (map_->asu[index_] != 0) edge(); return *this; }  //!< decrement v
      inline Map_reference_coord& prev_w() { pos_.w()--; index_ -= map_->dw[sym_]; if (map_->asu[index_] != 0) edge(); return *this; }  //!< decrement w
      //! Assignment operator from a coord
      inline Map_reference_coord& operator =( const Coord_grid& pos )
	{ return set_coord( pos ); }
      // inherited functions listed for documentation purposes
      //-- const Xmap_base& base_xmap() const;
      //-- const int& index() const;
      //-- bool last() const;

    protected:
      //! Current symop
      int sym_;
      //! Current coord
      Coord_grid pos_;

      //! Reset index for a different symop when we hit the map border
      void edge();
    };

    //! return a Map_reference_index for this map
    Map_reference_index first() const { return Map_reference_index( *this ); }
    //! return a Map_reference_coord for this map
    Map_reference_coord first_coord() const { return Map_reference_coord( *this ); }
    //! set/get default backend type
    static FFTtype& default_type() { return default_type_; }
  protected:
    ObjectCache<Xmap_cacheobj>::Reference cacheref;  //!< object cache reference
    const unsigned char* asu;  //!< fast access ptr
    const Isymop* isymop;      //!< fast access ptr
    const int* du;             //!< fast access ptr
    const int* dv;             //!< fast access ptr
    const int* dw;             //!< fast access ptr
    Grid_range asu_grid;       //!< fast access copy
    Grid_range map_grid;       //!< fast access copy
    int nsym;                  //!< fast access copy

    Cell cell_;                    //!< unit cell
    Spacegroup spacegroup_;        //!< spacegroup
    Grid_sampling grid_sam_;       //!< grid for the whole cell

    RTop<> rt_orth_grid;           //!< orth->grid operator
    RTop<> rt_grid_orth;           //!< grid->orth operator

    //! Null constructor, for later initialisation
    Xmap_base();
    //! initialiser
    void init( const Spacegroup& spacegroup, const Cell& cell, const Grid_sampling& grid_sam );
    inline void find_sym( const Coord_grid& base, int& index, int& sym ) const;
    void asu_error( const Coord_grid& pos ) const;

    static FFTtype default_type_;       //!< default backend type

    friend class Xmap_base::Map_reference_base;
    friend class Xmap_base::Map_reference_index;
    friend class Xmap_base::Map_reference_coord;
  };




  //! Xmap<T>: actual crystallographic map class
  /*!
    The crystallographic map class stores a map of arbitrary data
    type. Its main difference from a 3-d array is that the data extent
    appears to be infinite, and yet internally only a unique ASU is
    stored. Iterators provide efficient access to data.

    This is derived from Xmap_base, and adds the templatised data
    itself and the methods which deal with it.

    \note The following methods are inherited from Xmap_base but are
    documented here for convenience: cell(), spacegroup(),
    grid_sampling(), grid_asu(), in_asu(), multiplicity(),
    to_map_unit(), first(), first_coord().
  */
  template<class T> class Xmap : public Xmap_base
  {
  public:
    //! Null constructor, for later initialisation
    Xmap() {}
    //! constructor: from spacegroup, cell, and grid
    Xmap( const Spacegroup& spacegroup, const Cell& cell, const Grid_sampling& grid_sam ) { init( spacegroup, cell, grid_sam ); }
    //! initialiser: from spacegroup, cell, and grid
    void init( const Spacegroup& spacegroup, const Cell& cell, const Grid_sampling& grid_sam ) { Xmap_base::init( spacegroup, cell, grid_sam ); list.resize( cacheref.data().asu.size() ); }

    //! get data by Map_reference_index
    inline const T& operator[] (const Xmap_base::Map_reference_index& ix) const
      { return list[ix.index()]; }
    //! set data by Map_reference_index
    inline T& operator[] (const Xmap_base::Map_reference_index& ix)
      { return list[ix.index()]; }

    //! get data by Map_reference_coord
    inline const T& operator[] (const Xmap_base::Map_reference_coord& ix) const
      { return list[ix.index()]; }
    //! set data by Map_reference_coord
    inline T& operator[] (const Xmap_base::Map_reference_coord& ix)
      { return list[ix.index()]; }

    //! get a density value for an arbitrary position
    const T& get_data( const Coord_grid& pos ) const;
    //! set a density value for an arbitrary position
    void set_data( const Coord_grid& pos, const T& val );
    //! get data by index (not recommended)
    inline const T& get_data( const int& index ) const;
    //! set data by index (not recommended)
    bool set_data( const int& index, const T& val  );

    //! get map value for fractional coord using supplied interpolator
    template<class I> T interp( const Coord_frac& pos ) const;
    //! get map value and grad for fractional coord using supplied interpolator
    template<class I> void interp_grad( const Coord_frac& pos, T& val, Grad_frac<T>& grad ) const;
    //! get map value and curv for fractional coord using supplied interpolator
    template<class I> void interp_curv( const Coord_frac& pos, T& val, Grad_frac<T>& grad, Curv_frac<T>& curv ) const;
    //! get map value for map coord using supplied interpolator
    template<class I> T interp( const Coord_map& pos ) const;
    //! get map value and grad for map coord using supplied interpolator
    template<class I> void interp_grad( const Coord_map& pos, T& val, Grad_map<T>& grad ) const;
    //! get map value and curv for map coord using supplied interpolator
    template<class I> void interp_curv( const Coord_map& pos, T& val, Grad_map<T>& grad, Curv_map<T>& curv ) const;

    //! FFT from reflection list to map
    template<class H> void fft_from( const H& fphidata, const FFTtype type = Default );
    //! FFT from map to reflection list
    template<class H> void fft_to  ( H& fphidata, const FFTtype type = Default ) const;

    // inherited functions listed for documentation purposes
    //-- const Cell& cell() const;
    //-- const Spacegroup& spacegroup() const;
    //-- const Grid_sampling& grid_sampling() const;
    //-- const Grid_range& grid_asu() const;
    //-- inline Coord_grid coord_of( const int& index ) const;
    //-- inline int index_of( const Coord_grid& coord ) const;
    //-- const int multiplicity( const Coord_grid& pos ) const;
    //-- const Coord_grid to_map_unit( const Coord_grid& pos ) const;
    //-- const Map_reference_index first() const;
    //-- const Map_reference_coord first_coord() const;

    //! assignment operator: assigns a single value to the whole map
    const T& operator =( const T& value );
    //! add another map to this one
    const Xmap<T>& operator +=( const Xmap<T>& other );
    //! subtract another map from this one
    const Xmap<T>& operator -=( const Xmap<T>& other );

  private:
    std::vector<T> list;
  };



  // implementations

  void Xmap_base::find_sym( const Coord_grid& base, int& index, int& sym ) const
  {
    // first deal with first symop, and cache for highest performance
    Coord_grid rot = to_map_unit( base );
    if ( asu_grid.in_grid( rot ) ) {
      index = map_grid.index( rot );
      if ( asu[ index ] == 0 ) {
	sym = 0;
      } else {
	sym = asu[ index ] - 1;
	index = map_grid.index( to_map_unit( base.transform(isymop[sym]) ) );
      }
    } else {
      // now deal with other symops
      for ( sym = 1; sym < nsym; sym++ ) {
	rot = to_map_unit( base.transform( isymop[sym] ) );
	if ( asu_grid.in_grid( rot ) ) {
	  index = map_grid.index( rot );
	  if ( asu[ index ] == 0 ) return;
	}
      }
      index = 0;  // redundent, to avoid compiler warnings
      asu_error( base );
    }
    return;
  }


  /*! Accessing the data by coordinate, rather than by
    Map_reference_index or Map_reference_coord, involves a symmetry
    lookup and is therefore slow. Avoid using these methods when you
    can. */
  template<class T> const T& Xmap<T>::get_data( const Coord_grid& pos ) const
  {
    int index, sym;
    find_sym( pos, index, sym );
    return list[ index ];
  }

  /*! Accessing the data by coordinate, rather than by
    Map_reference_index or Map_reference_coord, involves a symmetry
    lookup and is therefore slow. Avoid using these methods when you
    can. */
  template<class T> void Xmap<T>::set_data( const Coord_grid& pos, const T& val )
  {
    int index, sym;
    find_sym( pos, index, sym );
    list[ index ] = val;
  }

  /*! Accessing the data by index, rather than by Map_reference_index
    or Map_reference_coord, is generally to be avoided since the
    indices do not start at zero and do not increase
    contiguously. These methods are only useful when a large number of
    references into a map must be stored, e.g. for sorting into
    density order. */
  template<class T> const T& Xmap<T>::get_data( const int& index ) const
  { return list[index]; }

  /*! Accessing the data by index, rather than by Map_reference_index
    or Map_reference_coord, is generally to be avoided since the
    indices do not start at zero and do not increase
    contiguously. These methods are only useful when a large number of
    references into a map must be stored, e.g. for sorting into
    density order.
    \return true if data was set, i.e. index is valid. */
  template<class T> bool Xmap<T>::set_data( const int& index, const T& val  )
  {
    if ( index >= 0 && index < list.size() )
      if ( asu[index] == 0 ) {
	list[index] = val;
	return true;
      }
    return false;
  }

  /*! The value of the map at the desired non-grid fractional
    coordinate are calculated using the supplied interpolator template.
    \code
    Coord_frac f( u, v, w );
    y = xmap.interp<Interp_cubic>( f );
    \endcode
    \param pos The fractional coord at which the density is to be calcuated.
    \return The value of the density at that point. */
  template<class T> template<class I> T Xmap<T>::interp( const Coord_frac& pos ) const
  {
    T val;
    I::interp( *this, pos.coord_map( grid_sam_ ), val );
    return val;
  }


  /*! The value of the map at the desired non-grid fractional
    coordinate and its gradient are calculated using
    the supplied interpolator template.
    \param pos The fractional coord at which the density is to be calcuated.
    \param val The value of the density at that point.
    \param grad The interpolated gradient vector with respect to the
    fractional coordinates (see Cell::coord_orth). */
  template<class T> template<class I> void Xmap<T>::interp_grad( const Coord_frac& pos, T& val, Grad_frac<T>& grad ) const
  {
    Grad_map<T> g;
    I::interp_grad( *this, pos.coord_map( grid_sam_ ), val, g );
    grad = g.grad_frac( grid_sam_ );
  }


  /*! The value of the map at the desired non-grid fractional
    coordinate and its gradient and curvature are calculated using
    the supplied interpolator template. e.g.
    \param pos The fractional coord at which the density is to be calcuated.
    \param val The value of the density at that point.
    \param grad The interpolated gradient vector with respect to the
    fractional coordinates (see Cell::coord_orth).
    \param curv The interpolated curvature matrix with respect to the
    fractional coordinates (see Cell::coord_orth). */
  template<class T> template<class I> void Xmap<T>::interp_curv( const Coord_frac& pos, T& val, Grad_frac<T>& grad, Curv_frac<T>& curv ) const
  {
    Grad_map<T> g;
    Curv_map<T> c;
    I::interp_curv( *this, pos.coord_map( grid_sam_ ), val, g, c );
    grad = g.grad_frac( grid_sam_ );
    curv = c.curv_frac( grid_sam_ );
  }


  /*! The value of the map at the desired non-grid map
    coordinate are calculated using the supplied interpolator template.
    \code
    Coord_map m( u, v, w );
    y = xmap.interp<Interp_cubic>( m );
    \endcode
    \param pos The map coord at which the density is to be calcuated.
    \return The value of the density at that point. */
  template<class T> template<class I> T Xmap<T>::interp( const Coord_map& pos ) const
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
    map coordinates (see Cell::coord_orth). */
  template<class T> template<class I> void Xmap<T>::interp_grad( const Coord_map& pos, T& val, Grad_map<T>& grad ) const
    { I::interp_grad( *this, pos, val, grad ); }


  /*! The value of the map at the desired non-grid map
    coordinate and its gradient and curvature are calculated using
    the supplied interpolator template. e.g.
    \param pos The map coord at which the density is to be calcuated.
    \param val The value of the density at that point.
    \param grad The interpolated gradient vector with respect to the
    map coordinates (see Cell::coord_orth).
    \param curv The interpolated curvature matrix with respect to the
    map coordinates (see Cell::coord_orth). */
  template<class T> template<class I> void Xmap<T>::interp_curv( const Coord_map& pos, T& val, Grad_map<T>& grad, Curv_map<T>& curv ) const
    { I::interp_curv( *this, pos, val, grad, curv ); }


  /*! An FFT is calculated using the provided reflection list of
    F_phi, and used to fill this map. The reflection list is unchanged.
    \param fphidata The reflection data list to use
  */
  template<class T> template<class H> void Xmap<T>::fft_from( const H& fphidata, const FFTtype type )
  {
    if ( type == Sparse || ( type == Default && default_type() == Sparse ) ) {
    //if ( false ) {
      // make a sparse fftmap
      FFTmap_sparse_p1_hx fftmap( grid_sampling() );
      // copy from reflection data
      typename H::HKL_reference_index ih;
      ffttype f, phi0, phi1;
      int sym;
      for ( ih = fphidata.first_data(); !ih.last(); fphidata.next_data( ih ) ) {
	f = fphidata[ih].f();
	if ( f != 0.0 ) {
	  phi0 = fphidata[ih].phi();
	  const HKL& hkl = ih.hkl();
	  fftmap.set_hkl( hkl,
			  std::complex<ffttype>( f*cos(phi0), f*sin(phi0) ) );
	  for ( sym = 1; sym < spacegroup_.num_primops(); sym++ ) {
	    phi1 = phi0 + hkl.sym_phase_shift( spacegroup_.symop(sym) );
	    fftmap.set_hkl( hkl.transform( isymop[sym] ),
			    std::complex<ffttype>( f*cos(phi1), f*sin(phi1) ) );
	  }
	}
      }
      // require output ASU coords
      for ( Map_reference_index ix = first(); !ix.last(); ix.next() )
	fftmap.require_real_data( ix.coord() );
      // do fft
      fftmap.fft_h_to_x(1.0/cell().volume());
      // fill map ASU
      for ( Map_reference_index ix = first(); !ix.last(); ix.next()) {
	         (*this)[ix] = fftmap.real_data( ix.coord() );
      }

    } else {
      // make a normal fftmap
      FFTmap_p1 fftmap( grid_sampling() );
      // copy from reflection data
      typename H::HKL_reference_index ih;
      ffttype f, phi0, phi1;
      int sym;
      for ( ih = fphidata.first_data(); !ih.last(); fphidata.next_data( ih ) ) {
	f = fphidata[ih].f();
	if ( f != 0.0 ) {
	  phi0 = fphidata[ih].phi();
	  const HKL& hkl = ih.hkl();
	  fftmap.set_hkl( hkl,
			  std::complex<ffttype>( f*cos(phi0), f*sin(phi0) ) );
	  for ( sym = 1; sym < spacegroup_.num_primops(); sym++ ) {
	    phi1 = phi0 + hkl.sym_phase_shift( spacegroup_.symop(sym) );
	    fftmap.set_hkl( hkl.transform( isymop[sym] ),
			    std::complex<ffttype>( f*cos(phi1), f*sin(phi1) ) );
	  }
	}
      }
      // do fft
      fftmap.fft_h_to_x(1.0/cell().volume());
      // fill map ASU
      for ( Map_reference_index ix = first(); !ix.last(); ix.next() )
	(*this)[ix] = fftmap.real_data( ix.coord() );
    }
  }


  /*! The Fourier transform of this map is calculated and used to fill
    a reflection list of F_phi. The map is unchanged.

    Arguably this should be part of hkl_data<F_phi<T>>. But that
    requires writing a specialisation of hkl_data for F_phi. This is
    simpler and imposes less demands on the compiler.
    \param fphidata The reflection data list to set.
  */
  template<class T> template<class H> void Xmap<T>::fft_to  ( H& fphidata, const FFTtype type ) const
  {
    if ( type == Sparse || ( type == Default && default_type() == Sparse ) ) {
      // make a sparse fftmap
      FFTmap_sparse_p1_xh fftmap( grid_sampling() );
      // copy from map data
      ffttype f;
      int sym;
      for ( Map_reference_index ix = first(); !ix.last(); ix.next() ) {
	f = (*this)[ix];
	if ( f != 0.0 ) {
	  fftmap.real_data( ix.coord() ) = f;
	  for ( sym = 1; sym < cacheref.data().nsym; sym++ )
	    fftmap.real_data(
              ix.coord().transform( isymop[sym] ).unit( grid_sam_ ) ) = f;
	}
      }
      // require output ASU coords
      typename H::HKL_reference_index ih;
      for ( ih = fphidata.first(); !ih.last(); ih.next() )
	fftmap.require_hkl( ih.hkl() );
      // do fft
      fftmap.fft_x_to_h(cell().volume());
      // fill data ASU
      for ( ih = fphidata.first(); !ih.last(); ih.next() ) {
	std::complex<ffttype> c = fftmap.get_hkl( ih.hkl() );
	fphidata[ih].f() = std::abs(c);
	fphidata[ih].phi() = std::arg(c);
      }
    } else {
      // make a normal fftmap
      FFTmap_p1 fftmap( grid_sampling() );
      // copy from map data
      ffttype f;
      int sym;
      for ( Map_reference_index ix = first(); !ix.last(); ix.next() ) {
	f = (*this)[ix];
	if ( f != 0.0 ) {
	  fftmap.real_data( ix.coord() ) = f;
	  for ( sym = 1; sym < cacheref.data().nsym; sym++ )
	    fftmap.real_data(
              ix.coord().transform( isymop[sym] ).unit( grid_sam_ ) ) = f;
	}
      }
      // do fft
      fftmap.fft_x_to_h(cell().volume());
      // fill data ASU
      typename H::HKL_reference_index ih;
      for ( ih = fphidata.first(); !ih.last(); ih.next() ) {
	std::complex<ffttype> c = fftmap.get_hkl( ih.hkl() );
	fphidata[ih].f() = std::abs(c);
	fphidata[ih].phi() = std::arg(c);
      }
    }
  }


  /*! All values, including missing values, are overwritten by the value.
    \param value The value to which the map is to be set. */
  template<class T> const T& Xmap<T>::operator =( const T& value )
  {
    // copy value into map
    for ( Map_reference_index im = first(); !im.last(); im.next() )
      list[im.index()] = value;
    return value;
  }


  /*! The map grids and spacegroups must match. */
  template<class T> const Xmap<T>& Xmap<T>::operator +=( const Xmap<T>& other )
  {
    if ( spacegroup().hash() != other.spacegroup().hash() ||
	 grid_sampling() != other.grid_sampling() )
      Message::message( Message_fatal( "Xmap: map mismatch in +=" ) );
    for ( Map_reference_index im = first(); !im.last(); im.next() )
      list[im.index()] += other[im];
    return (*this);
  }

  /*! The map grids and spacegroups must match. */
  template<class T> const Xmap<T>& Xmap<T>::operator -=( const Xmap<T>& other )
  {
    if ( spacegroup().hash() != other.spacegroup().hash() ||
	 grid_sampling() != other.grid_sampling() )
      Message::message( Message_fatal( "Xmap: map mismatch in -=" ) );
    for ( Map_reference_index im = first(); !im.last(); im.next() )
      list[im.index()] -= other[im];
    return (*this);
  }


} // namespace clipper

#endif
