/*!  \file convolution_search.h
  Header file for sample Convolution_search impelementation
  \ingroup g_convolution_search
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


#ifndef CLIPPER_CONVOLUTION_SEARCH
#define CLIPPER_CONVOLUTION_SEARCH


#include "function_object_bases.h"
#include "../imex.h"

namespace clipper {


  //! FFT-based Convolution_search implementation
  /*! \ingroup g_convolution_search
    This implementation is currently unoptimised, but much faster then
    the simple implementation. */
  template<class T> class Convolution_search_slow : public Convolution_search_base<T> {
  public:
    //! constructor
    Convolution_search_slow( const Xmap<T>& xmap ) { init( xmap ); }
    //! constructor: shorthand for constructor+operator
    Convolution_search_slow( Xmap<T>& result, const NXmap<T>& srchval, const Xmap<T>& xmap, const NX_operator& nxop ) { init( xmap ); (*this)( result, srchval, nxop ); }
    //! initialiser: initialise with the given target Xmap
    void init( const Xmap<T>& xmap ) { xmp = &xmap; }
    bool operator() ( Xmap<T>& result, const NXmap<T>& srchval, const NX_operator& nxop ) const;
    bool operator() ( Xmap<T>& result, const NXmap<T>& srchval ) const;  //!< \deprecated
  private:
    const Xmap<T>* xmp;
  };


  //! FFT-based Convolution_search implementation
  /*! \ingroup g_convolution_search
    This implementation is currently unoptimised, but much faster then
    the simple implementation. */
  template<class T> class Convolution_search_fft : public Convolution_search_base<T> {
  public:
    //! constructor
    Convolution_search_fft( const Xmap<T>& xmap ) { init( xmap ); }
    //! constructor: shorthand for constructor+operator
    Convolution_search_fft( Xmap<T>& result, const NXmap<T>& srchval, const Xmap<T>& xmap, const NX_operator& nxop ) { init( xmap ); (*this)( result, srchval, nxop ); }
    //! initialiser: initialise with the given target Xmap
    void init( const Xmap<T>& xmap );
    bool operator() ( Xmap<T>& result, const NXmap<T>& srchval, const NX_operator& nxop ) const;
    bool operator() ( Xmap<T>& result, const NXmap<T>& srchval ) const;  //!< \deprecated
  private:
    ftype vol;
    FFTmap_p1 rho1;
  };


} // namespace clipper

#endif
