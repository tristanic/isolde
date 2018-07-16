/* skeleton.cpp: Skeletonisation implementation */
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

#include "skeleton.h"

#include "../core/map_utils.h"

#include <set>


namespace clipper {


// Skeleton_basic

Skeleton_basic::Neighbours::Neighbours( const clipper::Xmap_base &map, const float min_distsq, const float max_distsq )
{
   /*
   note that thisd2 is a mesure of distance in real space

   So here we look around an arbitary current point and take note if
   the distance between the centre point and points one unit away from
   it (all directions) is less than a certain cuttoff distance.  We
   are left with a list of ndn neighbouring cells that are
   sufficicently close to be considered neighbours.
    
   Typically, ndn will be 18-20 at then end of this block.

   Note that body diagonals have length sqrt(3) ~= 1.73, but we passed
   squared limits.
   */

  clipper::Cell_descr rcd( map.cell().descr() );
  clipper::Cell_descr vcd( 1.0,1.0,1.0, rcd.alpha(), rcd.beta(), rcd.gamma() );
  clipper::Cell vcell( vcd );

  clipper::Coord_grid g0(-1,-1,-1); 
  clipper::Coord_grid g1( 1, 1, 1);
  clipper::Grid_sampling vgrid( 1, 1, 1 );

  clipper::Coord_grid iu, iv, iw;
  float thisd2;
  for ( iu = g0; iu.u() <= g1.u(); iu.u()++ ) {
     for ( iv = iu; iv.v() <= g1.v(); iv.v()++ ) {
        for ( iw = iv; iw.w() <= g1.w(); iw.w()++ ) {
           thisd2 = iw.coord_frac( vgrid ).lengthsq( vcell );
           if (thisd2 > min_distsq && thisd2 < max_distsq) nlist.push_back( iw );
        }
     }
  }
}

Skeleton_basic::NCube::NCube( const int& n )
{
  m = clipper::Grid_range( clipper::Coord_grid( -n, -n, -n ),
			 clipper::Coord_grid(  n,  n,  n ) );
  data.resize( m.size() );
}

bool Skeleton_basic::operator() ( clipper::Xmap<int>& xskl, const clipper::Xmap<float>& xmap ) const
{
  std::vector<int> index;
  clipper::Xmap<float>::Map_reference_index ix;
  
  /* now get the map in sorted order.
     We only sort those points which are to be considered, i.e. non-zero */
  for ( ix = xmap.first(); !ix.last(); ix.next() )
    if ( xskl[ix] > 0 )
      index.push_back( ix.index() );
  Map_index_sort::sort_increasing( xmap, index );

  /* make neighbours
     The neighbours of a point are the other grid points which are
     'near' it.  The exact choice depends on the grid geometry. The
     cutoff is chosen to give 18-20 neighbours. For a cubic grid,
     these are a 3x3x3 cube without the vertices, for a hex grid they
     are a hexagonal cylinder. */
  Skeleton_basic::Neighbours neigh( xmap );

  /* make the skeleton map. This will contain:
      0 for non-skeleton grids (inter ridge spaces)
      1 for untested grids
      1 for skeleton grids (ridges)
     (The untested and skeleton grids can have the same value,
     because the untested ones are still in the index).
     We loop through point in order, starting with the lowest
     and decide whether each one is part of the skeleton. */
  for ( int i = 0; i < index.size(); i++ )
    if ( !isInSkel( xskl, xskl.coord_of(index[i]), neigh, box_ ) )
      xskl.set_data( index[ i ], 0 );

  return true;
}

bool Skeleton_basic::isInSkel( const clipper::Xmap<int>& xskl, const clipper::Coord_grid& c, const Skeleton_basic::Neighbours& neigh, const int& box )
{
  Skeleton_basic::NCube cube( box );  // 1 or 2? Look at the results to decide.

  /* Fill the cube with flags Each non-rejected grid point is given
     its own number. We will the reduce these to 'region numbers', so
     that each connected region takes the value of one of its
     neighbours. (A performace benefit is gaine from allocating these
     in decreasing order). */
  clipper::Coord_grid g0 = c + cube.grid().min();
  clipper::Coord_grid g1 = c + cube.grid().max();
  clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;
  i0 = clipper::Xmap<int>::Map_reference_coord( xskl, g0 );
  int i = cube.grid().size();
  for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
    for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
      for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() )
	cube[ iw.coord() - c ] = ( xskl[iw] > 0 ) ? i-- : 0;
  // the centre cannot be  a link - set to zero
  cube[ clipper::Coord_grid(0,0,0) ] = 0;

  /* The following is a simple optimisation which could be omitted.
     We count the number of neightbours or the centre point which are
     not yet eliminated from the skeleton. If this is more than 14 or
     less than 1, then this pint cannot be part of the skeleton. */
  int nneigh = 0;
  for ( i = 0; i < neigh.size(); i++ )
    if ( cube[neigh[i]] > 0 ) nneigh++;
  if ( nneigh > 14 ) return false;
  if ( nneigh < 1 ) return false;

  /* Now reduce the flags to connected regions.  We consider each
     point in turn, and replace it with the highest value of any of
     its neighbours. This is repeated until no points change. When
     this occurs, each connected region will be set to the largest
     value of any of its initial members. */
  int val, nchg;
  do {
    nchg = 0;
    clipper::Coord_grid u, v, w, x;
    for ( u = cube.grid().min(); u.u() <= cube.grid().max().u(); u.u()++ )
      for ( v = u; v.v() <= cube.grid().max().v(); v.v()++ )
	for ( w = v; w.w() <= cube.grid().max().w(); w.w()++ ) {
	  val = cube[w];
	  if ( val != 0 ) {
	    for ( i = 0; i < neigh.size(); i++ ) {
	      x = w + neigh[i];
	      if ( cube.grid().in_grid( x ) )
		if ( cube[x] > val ) val = cube[x];
	    }
	    if ( val > cube[w] ) {
	      cube[w] = val;
	      nchg++;
	    }
	  }
	}
  } while ( nchg > 0 );

  /* The following code uses an STL set to count the number of different,
     non-zero numbers in the region bordering the centre of the cube */
  std::set<int> uniqnbrs;
  for ( i = 0; i < neigh.size(); i++ ) {
    const int& val = cube[neigh[i]];
    if ( val > 0 ) uniqnbrs.insert( val );
  }

  return (uniqnbrs.size() > 1);
}


// Skeleton_fast

template <class T1, class T2> Skeleton_fast<T1,T2>::Neighbours::Neighbours( const clipper::Xmap_base &map, const float min_distsq, const float max_distsq )
{
   /*
   note that thisd2 is a mesure of distance in real space

   So here we look around an arbitary current point and take note if
   the distance between the centre point and points one unit away from
   it (all directions) is less than a certain cuttoff distance.  We
   are left with a list of ndn neighbouring cells that are
   sufficicently close to be considered neighbours.
    
   Typically, ndn will be 18-20 at then end of this block.

   Note that body diagonals have length sqrt(3) ~= 1.73, but we passed
   squared limits.
   */

  clipper::Cell_descr rcd( map.cell().descr() );
  clipper::Cell_descr vcd( 1.0,1.0,1.0, rcd.alpha(), rcd.beta(), rcd.gamma() );
  clipper::Cell vcell( vcd );

  clipper::Coord_grid g0(-1,-1,-1); 
  clipper::Coord_grid g1( 1, 1, 1);
  clipper::Grid_sampling vgrid( 1, 1, 1 );

  clipper::Coord_grid iu, iv, iw;
  float thisd2;
  for ( iu = g0; iu.u() <= g1.u(); iu.u()++ ) {
     for ( iv = iu; iv.v() <= g1.v(); iv.v()++ ) {
        for ( iw = iv; iw.w() <= g1.w(); iw.w()++ ) {
           thisd2 = iw.coord_frac( vgrid ).lengthsq( vcell );
           if (thisd2 > min_distsq && thisd2 < max_distsq) nlist.push_back( iw );
        }
     }
  }
}


template<class T1, class T2> bool Skeleton_fast<T1,T2>::operator() ( clipper::Xmap<T1>& xskl, const clipper::Xmap<T2>& xmap ) const
{
  std::vector<int> index;
  Xmap_base::Map_reference_index ix;
  
  /* now get the map in sorted order.
     We only sort those points which are to be considered, i.e. non-zero */
  for ( ix = xmap.first(); !ix.last(); ix.next() )
    if ( xskl[ix] > 0 )
      index.push_back( ix.index() );
  Map_index_sort::sort_increasing( xmap, index );

  /* make neighbours
     The neighbours of a point are the other grid points which are
     'near' it.  The exact choice depends on the grid geometry. The
     cutoff is chosen to give 18-20 neighbours. For a cubic grid,
     these are a 3x3x3 cube without the vertices, for a hex grid they
     are a hexagonal cylinder. */
  neigh = Neighbours( xmap );

  /* make the skeleton map. This will contain:
      0 for non-skeleton grids (inter ridge spaces)
      1 for untested grids
      1 for skeleton grids (ridges)
     (The untested and skeleton grids can have the same value,
     because the untested ones are still in the index).
     We loop through point in order, starting with the lowest
     and decide whether each one is part of the skeleton. */
  for ( int i = 0; i < index.size(); i++ )
    if ( !isInSkel( xskl, xskl.coord_of(index[i]) ) )
      xskl.set_data( index[ i ], 0 );

  return true;
}


template<class T1, class T2> bool Skeleton_fast<T1,T2>::isInSkel( const clipper::Xmap<T1>& xskl, const clipper::Coord_grid& c ) const
{
  int dx, dy, dz;

  /* Fill the cube with flags Each non-rejected grid point is given
     its own number. We will the reduce these to 'region numbers', so
     that each connected region takes the value of one of its
     neighbours. (A performace benefit is gaine from allocating these
     in decreasing order). */
  clipper::Xmap_base::Map_reference_index ix( xskl, c );
  for ( dz = 0; dz < 3; dz++ )
    for ( dy = 0; dy < 3; dy++ )
      for ( dx = 0; dx < 3; dx++ )
	cube[dx][dy][dz] = xskl.get_data( ix.index_offset( dx-1, dy-1, dz-1 ) );
  // the centre cannot be  a link - set to zero
  cube[1][1][1] = 0;

  /* The following is a simple optimisation which could be omitted.
     We count the number of neightbours or the centre point which are
     not yet eliminated from the skeleton. If this is more than 14 or
     less than 1, then this pint cannot be part of the skeleton. */
  int nneigh = 0;
  for ( int i = 0; i < neigh.size(); i++ ) {
    dx = neigh[i].u()+1;
    dy = neigh[i].v()+1;
    dz = neigh[i].w()+1;
    if ( cube[dx][dy][dz] > 0 ) nneigh++;
  }
  if ( nneigh > 14 ) return false;
  if ( nneigh < 1 ) return false;

  /* Now alter flags for one connected region. If all flags are
     altered, then this point is not part of the skeleton. */
  for ( int i = 0; i < neigh.size(); i++ ) {
    dx = neigh[i].u()+1;
    dy = neigh[i].v()+1;
    dz = neigh[i].w()+1;
    if ( cube[dx][dy][dz] > 0 ) break;
  }
  flood_cube( dx, dy, dz );
  for ( int i = 0; i < neigh.size(); i++ ) {
    dx = neigh[i].u()+1;
    dy = neigh[i].v()+1;
    dz = neigh[i].w()+1;
    if ( cube[dx][dy][dz] > 0 ) return true;
  }
  return false;
}

template<class T1, class T2> void Skeleton_fast<T1,T2>::flood_cube( const int x, const int y, const int z ) const
{
  cube[x][y][z] = -1;
  for ( int i = 0; i < neigh.size(); i++ ) {
    int dx = x + neigh[i].u();
    if ( dx >= 0 && dx < 3 ) {
      int dy = y + neigh[i].v();
      if ( dy >= 0 && dy < 3 ) {
	int dz = z + neigh[i].w();
	if ( dz >= 0 && dz < 3 )
	  if ( cube[dx][dy][dz] > 0 )
	    flood_cube( dx, dy, dz );
      }
    }
  }
}

// template instantiations

template class Skeleton_fast<char,float>;
template class Skeleton_fast<char,double>;
template class Skeleton_fast<short,float>;
template class Skeleton_fast<short,double>;
template class Skeleton_fast<int,float>;
template class Skeleton_fast<int,double>;


} // namespace clipper
