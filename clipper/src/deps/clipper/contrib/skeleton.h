/*!  \file skeleton.h
  Header file for sample skeletonisation impelementation
  \ingroup g_skel
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


#ifndef CLIPPER_SKELETON
#define CLIPPER_SKELETON


#include "function_object_bases.h"
#include "../imex.h"

namespace clipper {


  //! Simple skeletonisation implementation
  /*! \ingroup g_skel */
  class CLIPPER_IMEX Skeleton_basic : public Skeleton_base<int,float> {
  public:
    // helper classes

    class Neighbours {
    public:
      Neighbours( const clipper::Xmap_base &map, const float min_dist = 0.5, const float max_dist = 2.5 ); 
      clipper::Coord_grid operator[] (int i) const { return nlist[i]; }
      int size() const { return nlist.size(); }
    private:
      std::vector<clipper::Coord_grid> nlist;
    };

    class NCube {
    public:
      NCube( const int& n );
      const int& operator[] ( const clipper::Coord_grid& c ) const
	{ return data[m.index(c)]; }
      int& operator[] ( const clipper::Coord_grid& c )
	{ return data[m.index(c)]; }
      const clipper::Grid_range& grid() const { return m; }
    private:
      clipper::Grid_range m;
      std::vector<int> data;
    };

    // methods
    //! constructor
    Skeleton_basic( const int box = 1 ) : box_(box) {}
    //! constructor: shorthand for constructor+operator
    Skeleton_basic( clipper::Xmap<int>& xskl, const clipper::Xmap<float>& xmap, const int box = 1 ) : box_(box) { (*this)( xskl, xmap ); }
    //! Skeletonise a map
    bool operator() ( clipper::Xmap<int>& xskl, const clipper::Xmap<float>& xmap ) const;

  private:
    static bool isInSkel( const clipper::Xmap<int>& xskl, const clipper::Coord_grid& c, const Skeleton_basic::Neighbours& neigh, const int& box );
    int box_;
  };


  //! Fast template skeletonisation implementation
  /*! \ingroup g_skel */
  template <class T1, class T2> class Skeleton_fast : public Skeleton_base<T1,T2> {
  public:
    // helper classes

    class Neighbours {
    public:
      Neighbours() {}
      Neighbours( const clipper::Xmap_base &map, const float min_dist = 0.5, const float max_dist = 2.5 ); 
      clipper::Coord_grid operator[] (int i) const { return nlist[i]; }
      int size() const { return nlist.size(); }
    private:
      std::vector<clipper::Coord_grid> nlist;
    };

    // methods
    //! constructor
    Skeleton_fast() {}
    //! constructor: shorthand for constructor+operator
    Skeleton_fast( clipper::Xmap<T1>& xskl, const clipper::Xmap<T2>& xmap ) { (*this)( xskl, xmap ); }
    //! Skeletonise a map
    bool operator() ( clipper::Xmap<T1>& xskl, const clipper::Xmap<T2>& xmap ) const;

  private:
    mutable int cube[3][3][3];
    mutable Neighbours neigh;
    bool isInSkel( const clipper::Xmap<T1>& xskl, const clipper::Coord_grid& c ) const;
    void flood_cube( const int x, const int y, const int z ) const;
  };


} // namespace clipper

#endif
