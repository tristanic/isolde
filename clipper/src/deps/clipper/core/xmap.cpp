/* xmap.cpp: implementation file for crystal maps */
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


#include "xmap.h"
#include "clipper_instance.h"
#include <future>

namespace clipper {


Message_ctor message_ctor_xmap( " [Xmap: constructed]" );


Mutex Xmap_cacheobj::mutex = Mutex();

Xmap_cacheobj::Xmap_cacheobj( const Key& xmap_cachekey ) :
  key( xmap_cachekey )
{
  Spacegroup spacegroup_( xmap_cachekey.spgr_descr() );

  xtl_grid = xmap_cachekey.grid_sampling();

  // Get the map grid - must contain ASU
  map_grid = asu_grid =
    Grid_range( xtl_grid, spacegroup_.asu_min(), spacegroup_.asu_max() );
  map_grid.add_border(1);

  // Create the intergised symops (assumes legal grid) &
  // Precalculate shifts to index due to symmetry transformed grid steps
  int sym;
  nsym = spacegroup_.num_symops();
  isymop.resize( nsym );
  du.resize( nsym );
  dv.resize( nsym );
  dw.resize( nsym );
  for ( sym = 0; sym < nsym; sym++ ) {
    Isymop op = Isymop( spacegroup_.symop(sym), xtl_grid );
    isymop[sym] = op;
    du[sym] = map_grid.Grid::index( Coord_grid( op.rot()*Coord_grid(1,0,0) ) );
    dv[sym] = map_grid.Grid::index( Coord_grid( op.rot()*Coord_grid(0,1,0) ) );
    dw[sym] = map_grid.Grid::index( Coord_grid( op.rot()*Coord_grid(0,0,1) ) );
  }

  // now store the symmetry permutations
  symperm.resize( nsym, nsym );
  for ( int s1 = 0; s1 < nsym; s1++ )
    for ( int s2 = 0; s2 < nsym; s2++ )
      symperm( s1, s2 ) = spacegroup_.product_op( s1, s2 );


  // Flag all non-ASU points in the grid with sym number + 1
  Coord_grid base, rot;
  asu.clear();
  asu.resize( map_grid.size(), 255 );   // set all non-asu flags to 255
  find_asu_sym_(asu_grid.min(), asu_grid.max());
  // for ( base = asu_grid.min(); !base.last(asu_grid); base.next(asu_grid) ) {
  //   for ( sym = 1; sym < nsym; sym++ ) {
  //     rot = base.transform(isymop[sym]).unit(xtl_grid);
  //     if ( asu_grid.in_grid( rot ) )
	// if ( asu[ map_grid.index( rot ) ] == 0 ) break;
  //   }
  //   if ( sym == nsym ) asu[ map_grid.index( base ) ] = 0;
  // }
  find_map_sym_(map_grid.min(), map_grid.max());
  // for ( base = map_grid.min(); !base.last(map_grid); base.next(map_grid) )
  //   if ( asu[ map_grid.index( base ) ] == 255 ) {
  //     for ( sym = 0; sym < nsym; sym++ ) {
	// rot = base.transform(isymop[sym]).unit(xtl_grid);
	// if ( asu_grid.in_grid( rot ) )
	//   if ( asu[ map_grid.index( rot ) ] == 0 ) break;
  //     }
  //     asu[ map_grid.index( base ) ] = sym + 1;
  //   }

}

void Xmap_cacheobj::find_asu_sym_(const Coord_grid& begin,
                                  const Coord_grid& end)
{
    auto len = asu_grid.index(end) - asu_grid.index(begin);
    Coord_grid base, rot;
    int sym;
    if (len<loops_per_thread_)
    {
        for ( base = begin; base!=end; base.next(asu_grid) ) {
          for ( sym = 1; sym < nsym; sym++ ) {
            rot = base.transform(isymop[sym]).unit(xtl_grid);
            if ( asu_grid.in_grid( rot ) )
      	if ( asu[ map_grid.index( rot ) ] == 0 ) break;
          }
          if ( sym == nsym ) asu[ map_grid.index( base ) ] = 0;
        }
    } else
    {
        auto mid = asu_grid.deindex((asu_grid.index(begin)+asu_grid.index(end))/2);
        auto handle = std::async(std::launch::async, &Xmap_cacheobj::find_asu_sym_, this, mid, end);
        find_asu_sym_(begin, mid);
        handle.get();
    }
}

void Xmap_cacheobj::find_map_sym_(const Coord_grid& begin,
                                  const Coord_grid& end)
{
    auto len = map_grid.index(end) - map_grid.index(begin);
    Coord_grid base, rot;
    int sym;
    if (len<loops_per_thread_)
    {
        for ( base = begin; base != end; base.next(map_grid) )
        {
          if ( asu[ map_grid.index( base ) ] == 255 ) {
            for ( sym = 0; sym < nsym; sym++ ) {
      	rot = base.transform(isymop[sym]).unit(xtl_grid);
      	if ( asu_grid.in_grid( rot ) )
      	  if ( asu[ map_grid.index( rot ) ] == 0 ) break;
            }
            asu[ map_grid.index( base ) ] = sym + 1;
          }
        }
    } else
    {
        auto mid = map_grid.deindex((map_grid.index(begin)+map_grid.index(end))/2);
        auto handle = std::async(std::launch::async, &Xmap_cacheobj::find_map_sym_, this, mid, end);
        find_map_sym_(begin, mid);
        handle.get();

    }
}

bool Xmap_cacheobj::matches( const Key& xmap_cachekey ) const
{
  return
    key.spgr_descr().hash()  == xmap_cachekey.spgr_descr().hash() &&
    key.grid_sampling().nu() == xmap_cachekey.grid_sampling().nu() &&
    key.grid_sampling().nv() == xmap_cachekey.grid_sampling().nv() &&
    key.grid_sampling().nw() == xmap_cachekey.grid_sampling().nw();
}

String Xmap_cacheobj::format() const
{
  return key.spgr_descr().symbol_hall() + " " + key.grid_sampling().format();
}

Xmap_base::FFTtype Xmap_base::default_type_ = Xmap_base::Sparse;
//Xmap_base::FFTtype Xmap_base::default_type_ = Xmap_base::Normal;

/*! For later initialisation: see init() */
Xmap_base::Xmap_base()
{
  Message::message( message_ctor_xmap );
}

Xmap_base::Map_reference_coord& Xmap_base::Map_reference_coord::set_coord( const Coord_grid& pos )
{
  // use pos_ as a temporary variable to try out the current symop
  pos_ = map_->to_map_unit( pos.transform(map_->isymop[sym_]) );
  if ( map_->asu_grid.in_grid( pos_ ) ) {
    index_ = map_->map_grid.index( pos_ );
    if ( map_->asu[ index_ ] == 0 ) {
      pos_ = pos;
      return *this;
    }
  }
  map_->find_sym( pos, index_, sym_ );  // general case
  pos_ = pos;  // store the unmodified coord
  return *this;
}

void Xmap_base::Map_reference_coord::edge()
{
  int newsym = map_->asu[index_]-1;
  index_ = map_->map_grid.index( map_->to_map_unit( map_->map_grid.deindex(index_).transform( map_->isymop[newsym] ) ) );
  sym_ = map_->cacheref.data().symperm( newsym, sym_ );
}


/*! The Xmap is initialised with a given cell, spacegroup, and grid
  sampling. A unique assymetric unit (ASU) of grid cells is selected
  and will be used to store a unique set of data.

  If any of the parameters have null values, the existing values will
  be unchanged. The object will only be fully initialised once all
  parameters are available.
  \param spacegroup The spacegroup for the map
  \param cell       The cell for the map
  \param grid_sam   The grid sampling for the map,
    i.e. the sampling along each axis for one whole cell */
void Xmap_base::init( const Spacegroup& spacegroup, const Cell& cell, const Grid_sampling& grid_sam )
{
  // set parameters
  spacegroup_ = spacegroup;
  cell_ = cell;
  grid_sam_ = grid_sam;

  // check parameters
  if ( is_null() ) return;

  // get cache ref and fast access copies
  Xmap_cacheobj::Key key( spacegroup, grid_sam );

  cacheref = ClipperInstantiator::instance().xmap_cache().cache( key );
  asu      = &(cacheref.data().asu[0]);
  isymop   = &(cacheref.data().isymop[0]);
  du       = &(cacheref.data().du[0]);
  dv       = &(cacheref.data().dv[0]);
  dw       = &(cacheref.data().dw[0]);
  asu_grid = cacheref.data().asu_grid;
  map_grid = cacheref.data().map_grid;
  nsym     = cacheref.data().nsym;

  // store orthogonal conversions
  rt_grid_orth = RTop<>( cell_.matrix_orth()*grid_sam_.matrix_grid_frac() );
  rt_orth_grid = rt_grid_orth.inverse();

}


/*! \return true if the object has not been initalised. */
bool Xmap_base::is_null() const
{ return ( spacegroup_.is_null() || cell_.is_null() || grid_sam_.is_null() ); }


/*! The multiplicity is the number of times the spacegroup operators
  map a particular grid point onto itself. This is required in order
  to properly weight map statistics so as to get the same result from
  just an ASU as using the whole cell.
  \param pos The coordinate of the grid point.
  \return The multiplicty of the point.
*/
int Xmap_base::multiplicity( const Coord_grid& pos ) const
{
  int mult = 1;
  Coord_grid base = pos.unit(grid_sam_);
  for ( int sym = 1; sym < cacheref.data().nsym; sym++ )
    if ( base.transform(isymop[sym]).unit(grid_sam_) == base ) mult++;
  return mult;
}

void Xmap_base::asu_error( const Coord_grid& pos ) const
{
  std::cerr << "Failure to find grid coordinate " << pos.format() << std::endl;
  std::cerr << "Possible integer overflow or conversion from NaN" << std::endl;
  Message::message( Message_fatal( "Xmap: Internal map ASU error - " +
				   cacheref.data().format() ) );
}


// compile templates

template class Xmap<unsigned char>;
template class Xmap<char>;
template class Xmap<unsigned short>;
template class Xmap<short>;
template class Xmap<unsigned int>;
template class Xmap<int>;
template class Xmap<ftype32>;
template class Xmap<ftype64>;


} // namespace clipper
