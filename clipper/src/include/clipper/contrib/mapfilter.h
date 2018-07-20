/*!  \file mapfilter.h
  Header file for sample map filtering impelementation
  \ingroup g_mapf
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


#ifndef CLIPPER_MAPFILTER
#define CLIPPER_MAPFILTER


#include "function_object_bases.h"
#include "../imex.h"

namespace clipper {

  //! Base class for radial map filters
  class CLIPPER_IMEX MapFilterFn_base {
  public:
    virtual ~MapFilterFn_base() {};
    virtual ftype operator() ( const ftype& radius ) const = 0;
  };

  //! Simple slow convolution-based radial mapfiltering implementation
  /*! This version is of course very slow, and is mainly provided so
    you can test the repcision of the fft version.
    \ingroup g_mapf */
  template<class T> class MapFilter_slow : public MapFilter_base<T> {
  public:
    //! Scaling modes
    enum TYPE { NONE, Absolute, Relative };
    //! constructor
    MapFilter_slow( const MapFilterFn_base& fltr, const ftype scale = 1.0, const TYPE type = NONE );
    //! constructor: shorthand for constructor+operator
    MapFilter_slow( clipper::Xmap<T>& result, const clipper::Xmap<T>& xmap, MapFilterFn_base& fltr, const ftype scale = 1.0, const TYPE type = NONE );
    bool operator() ( clipper::Xmap<T>& result, const clipper::Xmap<T>& xmap ) const;
  private:
    const MapFilterFn_base* fltr_;
    ftype scale_;
    TYPE type_;
  };

  //! Simple fft-based radial mapfiltering implementation
  /*! The FFT method is fast, and also gives good precision.

    The following example demonstrates how to use the MapFilter to
    calculate the local mean and local deviation of an electron
    density map, in 'xmap':
    \code
    // make squared map
    clipper::Xmap<float> xmap2( xmap );
    clipper::Xmap<float>::Map_reference_index ix;
    for ( ix = xmap2.first(); !ix.last(); ix.next() )
      xmap2[ix] = pow( xmap2[ix], 2.0 );

    // now calculate local mean, local mean squared
    clipper::MapFilterFn_step fn( filter_radius );
    clipper::MapFilter_fft<float> fltr( fn, 1.0, clipper::MapFilter_fft<float>::Relative );
    clipper::Xmap<float> lmean, lsigm;
    fltr( lmean, xmap );
    fltr( lsigm, xmap2 );

    // calculate std deviation
    for ( ix = lmean.first(); !ix.last(); ix.next() )
      lsigm[ix] = sqrt( lsigm[ix] - pow( lmean[ix], 2.0 ) );
    \endcode

    This would be a useful step in solvent mask determination, for example.

    \ingroup g_mapf */
  template<class T> class MapFilter_fft : public MapFilter_base<T> {
  public:
    //! Scaling modes
    enum TYPE { NONE, Absolute, Relative };
    //! constructor
    MapFilter_fft( const MapFilterFn_base& fltr, const ftype scale = 1.0, const TYPE type = NONE );
    //! constructor: shorthand for constructor+operator
    MapFilter_fft( clipper::Xmap<T>& result, const clipper::Xmap<T>& xmap, MapFilterFn_base& fltr, const ftype scale = 1.0, const TYPE type = NONE );
    bool operator() ( clipper::Xmap<T>& result, const clipper::Xmap<T>& xmap ) const;
    bool operator() ( clipper::NXmap<T>& result, const clipper::NXmap<T>& nxmap ) const;
  private:
    const MapFilterFn_base* fltr_;
    ftype scale_;
    TYPE type_;
  };





  //! Step-function radial map filter
  /*! \ingroup g_mapf */
  class CLIPPER_IMEX MapFilterFn_step : public MapFilterFn_base {
  public:
    //! constructor: takes radius for step function cutoff
    MapFilterFn_step( const ftype& radius ) : radius_( radius ) {}
    //! destructor
    MapFilterFn_step() {}
    //! evaluate radial step function: 1.0 if inside or 0.0 outside
    ftype operator() ( const ftype& radius ) const
      { return (radius<radius_)?1.0:0.0; }
  private:
    ftype radius_;
  };

  //! Linear-function radial map filter
  /*! \ingroup g_mapf */
  class CLIPPER_IMEX MapFilterFn_linear : public MapFilterFn_base {
  public:
    //! constructor: takes radius for step function cutoff
    MapFilterFn_linear( const ftype& radius ) : radius_( radius ) {}
    //! destructor
    MapFilterFn_linear() {}
    //! evaluate radial step function: 1-r/r0 if inside or 0.0 outside
    ftype operator() ( const ftype& radius ) const
      { return (radius<radius_)?(1.0-radius/radius_):0.0; }
  private:
    ftype radius_;
  };

  //! Quadratic-function radial map filter
  /*! \ingroup g_mapf */
  class CLIPPER_IMEX MapFilterFn_quadratic : public MapFilterFn_base {
  public:
    //! constructor: takes radius for step function cutoff
    MapFilterFn_quadratic( const ftype& radius ) : radius_( radius ) {}
    //! destructor
    MapFilterFn_quadratic() {}
    //! evaluate radial step function: (1-r/r0)^2 if inside or 0.0 outside
    ftype operator() ( const ftype& radius ) const
      { return (radius<radius_)?pow(1.0-radius/radius_,2):0.0; }
  private:
    ftype radius_;
  };

} // namespace clipper

#endif
