/*! \file lib/map_utils.h
    Header file for map statistics type
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


#ifndef CLIPPER_MAP_UTILS
#define CLIPPER_MAP_UTILS


#include "derivs.h"
#include "../imex.h"

namespace clipper
{

  //! Generic map statistics class
  /*! This class is used to calculate and store the mean and standard
    deviation of a generic map object of scalar types (e.g. Xmap,
    NXmap). If the map contains NaN values, those points are excluded
    for the calculation. In the case of an Xmap, the proper
    multiplicty corrections are applied to give statistics for a whole
    unit cell */
  class CLIPPER_IMEX Map_stats
  { 
  public:
    Map_stats() {}                                //!< null constructor
    template<class M> Map_stats( const M& map );  //!< Constructor: from Xmap
    const ftype& mean() const { return mean_; }        //!< Mean of map
    const ftype& std_dev() const { return std_dev_; }  //!< Std deviation of map
    const ftype& min() const { return min_; }          //!< Minimum of map
    const ftype& max() const { return max_; }          //!< Maximum of map
    const Range<> range() const { return Range<>( min_, max_ ); } //!< Range
  private:
    ftype mean_, std_dev_, min_, max_;
  };


  //! Generic map sorting class
  /*! This class is used to sort a vector of integer indices into a
    map. This includes sorting the whole map to get highest or lowest
    density first, or sorting some subset, e.g. a list of peak
    indices. Integer indices are used because they are the most
    compact way of referencing a unique map location. e.g.
    \code
    clipper::Xmap<float>::Map_reference_index ix;
    std::vector<int> index;
    for ( ix = xmap.first(); !ix.last(); ix.next() )
      index.push_back( ix.index() );
    Map_index_sort::sort_decreasing( xmap, index );
    \endcode
  */
  class CLIPPER_IMEX Map_index_sort
  {
  public:
    //! Sort a list into increasing order
    template<class M> static void sort_increasing( const M& map, std::vector<int>& index );
    //! Sort a list into decreasing order
    template<class M> static void sort_decreasing( const M& map, std::vector<int>& index );
    //! Internal helper class used for sorting
  private:
    template<class M> class Compare_density {
    public:
      Compare_density( const M& m ) { p = &m; }
      bool operator() ( const int& i1, const int& i2 ) const { return p->get_data(i1) < p->get_data(i2); }
    private:
      const M* p;
    };
  };


  // template implementations

  /*! For float and double maps
    \params map The map for which moments are to be calculated. */
  template<class M> Map_stats::Map_stats( const M& map )
  {
    ftype64 w, x, s, sx, sxx;
    s = sx = sxx = 0.0;
    min_ =  1.0e12;
    max_ = -1.0e12;
    for ( typename M::Map_reference_index im = map.first();
	  !im.last(); im.next() ) {
      w = 1.0 / ftype64( map.multiplicity( im.coord() ) );
      x = ftype64( map[im] );
      if ( !Util::is_nan(x) ) {
	s += w;
	sx += w*x;
	sxx += w*x*x;
	if ( x < min_ ) min_ = x;
	if ( x > max_ ) max_ = x;
      }
    }
    sx /= s;
    sxx /= s;
    mean_ = sx;
    std_dev_ = sqrt( sxx - sx*sx );
  }


} // namespace clipper

#endif
