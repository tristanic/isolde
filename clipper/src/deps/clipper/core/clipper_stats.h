/*! \file lib/clipper_stats.h
    Header file for clipper helper functions
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


#ifndef CLIPPER_STATS
#define CLIPPER_STATS


#include "clipper_types.h"
#include "../imex.h"

namespace clipper
{

  //! Range - upper and lower bounds of some type
  template<class T = ftype> class Range
  {
  public:
    //! null constructor
    inline Range() { lmin = 999999999; lmax = -999999999; }
    //! constructor
    inline Range( const T& min, const T& max ) { lmin = min; lmax = max; }
    inline const T& min() const { return lmin; } //!< minimum value
    inline const T& max() const { return lmax; } //!< maximum value
    inline T range() const { return lmax-lmin; } //!< range = max - min
    //! update limits to include a new datum
    inline void include( const T& datum )
      { lmin = (datum<lmin)?datum:lmin; lmax = (datum>lmax)?datum:lmax; }
    //! test if data is within limits ( min <= datum <= max )
    inline bool contains( const T& datum ) const
      { return ( datum >= lmin && datum <= lmax ); }
    //! truncate data to be within range
    inline T truncate( const T& datum ) const
      { return Util::bound( lmin, datum, lmax ); }
  private:
    T lmin, lmax;
  };

  //! Range sampling: discrete sampling of a real range.
  class CLIPPER_IMEX Range_sampling : public Range<ftype>
  {
  public:
    //! null constructor
    inline Range_sampling() : n_(0) {}
    //! constructor: from number of samplings
    inline Range_sampling( const int& n ) : n_(n) {}
    //! constructor: from range and number of samplings
    inline Range_sampling( const Range<ftype>& range, const int& n ) :
      Range<ftype>( range ), n_(n) {}
    //! return fractional posn in counting range from x-value (0..n)
    inline ftype indexf( const ftype& x ) const
      { return ftype(size())*(x-min())/range(); }
    //! return x-value  (0..n) from fractional posn in counting range
    inline ftype x( const ftype& i ) const
      { return range()*i/ftype(size())+min(); }
    //! return nearest index to particular x-value
    inline int index( const ftype& x ) const { return Util::intf( indexf(x) ); }
    //! return nearest index to particular x-value (bounded 0...n-1)
    inline int index_bounded( const ftype& x ) const
      { return Util::bound( 0, Util::intf( indexf(x) ), size()-1 ); }
    //! return x-value corresponding to centre of i'th range
    inline ftype x( const int& i ) const     { return x( ftype(i)+0.5 ); }
    //! return x-value corresponding to bottom of i'th range
    inline ftype x_min( const int& i ) const { return x( ftype(i)     ); }
    //! return x-value corresponding to top of i'th range
    inline ftype x_max( const int& i ) const { return x( ftype(i)+1.0 ); }
    //! return number of samplings in range
    inline int size() const { return n_; }
  private:
    int n_;
  };

  //! General histogram class
  /*! This class is used to accumulate and access a histogram of
    values spread over a specified range. On storing data or
    retrieving by interpolation the range is checked. */
  class CLIPPER_IMEX Histogram : public Range_sampling
  {
  public:
    //! null constructor
    Histogram() {}
    //! constructor: from range and sampling
    Histogram( const Range<ftype>& range, const int& n ) :
      Range_sampling( range, n ), data( n, 0.0 ) {}
    //! add value to histogram (if it is in range)
    void accumulate( const ftype& x )
      { if ( contains(x) ) data[ index_bounded(x) ] += 1.0; }
    //! add specified value to histogram (if it is in range)
    void accumulate( const ftype& x, const ftype& w )
      { if ( contains(x) ) data[ index_bounded(x) ] += w; }
    //! return sum of whole histogram
    ftype sum() const;
    //! return value at index in histogram (Note: no bound check on i)
    inline const ftype& y( const int& i ) const { return data[i]; }
    //! return value at interpolated position in histogram
    ftype y( const ftype& x ) const;
    //! add the contents of two histograms (size must match)
    const Histogram& operator += ( const Histogram& h );
    // inherited functions listed for documentation purposes
    //-- inline ftype x( const int& i ) const;
    //-- inline ftype x_min( const int& i ) const;
    //-- inline ftype x_max( const int& i ) const;
    //-- inline int size() const;
  private:
    std::vector<ftype> data;
  };


  //! Generic ordinal gernerator
  /*! This is a generic fast ordinal calculator. It is supplied with a
    list of values, from which it prepares a cumulative distribution
    function. This may the used to return the approximate fracitonal
    ordinal (in the range 0...1) for any given value from the
    distibution.

    The distibution may be initialised by providing a vector of values
    from the distribution, or by adding the values and calling
    prep_ordinal().

    This distribution may also be inverted. Generation of a value from
    an ordinal may be used for generating random values from a given
    distribution, or for histogram matching. */
  class CLIPPER_IMEX Generic_ordinal
  {
  public:
    //! null constructor
    Generic_ordinal() {}
    //! constructor: from range and sampling
    Generic_ordinal( const Range<ftype>& range, const int& n )
      { init( range, n ); }
    //! initialiser: takes the source range and sampling
    void init( const Range<ftype>& range, const int num_ranges = 1000 );
    //! initialiser: takes the source distibution and a number of bins
    void init( const std::vector<ftype>& values, const int num_ranges = 1000 );
    //! return reflection ordinal
    ftype ordinal( const ftype& value ) const;

    //! accumulate values to build the distribution
    void accumulate( const ftype& value );
    //! accumulate values to build the distribution
    void accumulate( const ftype& value, const ftype& weight );
    //! generate the ordinal histogram
    void prep_ordinal();
    //! invert distribution to get value from ordinal
    void invert();

    //! DEPRECATED: initialiser: takes a number of bins for histogram
    void init( const int num_ranges = 1000 );
    //! DEPRECATED: add a value to the distribution (pass 1 of 2)
    void add_pass_1( const ftype& value );
    //! DEPRECATED: add a value to the distribution (pass 2 of 2)
    void add_pass_2( const ftype& value );
  protected:
    ftype nranges;            //!< number of ranges
    Range<ftype> range_;      //!< resolution range of data
    std::vector<ftype> hist;  //!< histogram of reflections vs resolution
  };



} // namespace clipper

#endif
