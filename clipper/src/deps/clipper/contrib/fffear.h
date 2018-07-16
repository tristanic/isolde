/*!  \file fffear.h
  Header file for sample fffear impelementation
  \ingroup g_fffear
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


#ifndef CLIPPER_FFFEAR
#define CLIPPER_FFFEAR


#include "function_object_bases.h"
#include "../imex.h"

namespace clipper {


  //! Simple fffear implementation
  /*! \deprecated
    This is for testing purposes only, as it is very slow.

    The search target and weights are not handled generally in this
    implementation: The NXmaps are assumed to have a grid matching the
    map grid with the origin at the centre. */
  class CLIPPER_IMEX FFFear_slow_basic {
  public:
    //! constructor
    FFFear_slow_basic( const Xmap<float>& xmap ) : xmp( &xmap ) {}
    //! constructor: shorthand for constructor+operator
    FFFear_slow_basic( Xmap<float>& result, const NXmap<float>& srchval, const NXmap<float>& srchwgt, const Xmap<float>& xmap ) : xmp( &xmap ) { (*this)( result, srchval, srchwgt ); }
    bool operator() ( Xmap<float>& result, const NXmap<float>& srchval, const NXmap<float>& srchwgt ) const;
  private:
    const Xmap<float>* xmp;
  };

  //! FFT-based fffear implementation
  /*! \deprecated
    This implementation is currently unoptimised, but much faster then
    the simple implementation.

    The search target and weights are not handled generally in this
    implementation: The NXmaps are assumed to have a grid matching the
    map grid with the origin at the centre. */
  class CLIPPER_IMEX FFFear_fft_basic {
  public:
    //! constructor
    FFFear_fft_basic( const Xmap<float>& xmap ) { init( xmap ); }
    //! constructor: shorthand for constructor+operator
    FFFear_fft_basic( Xmap<float>& result, const NXmap<float>& srchval, const NXmap<float>& srchwgt, const Xmap<float>& xmap ) { init( xmap ); (*this)( result, srchval, srchwgt ); }
    //! initialiser: initialise with the given target Xmap
    void init( const Xmap<float>& xmap );
    bool operator() ( Xmap<float>& result, const NXmap<float>& srchval, const NXmap<float>& srchwgt ) const;
  private:
    ftype vol;
    FFTmap_p1 rho1;
    FFTmap_p1 rho2;
  };


  //! FFT-based fffear implementation
  /*! \ingroup g_fffear
    This implementation is currently unoptimised, but much faster then
    the simple implementation. */
  template<class T> class FFFear_slow : public FFFear_base<T> {
  public:
    //! constructor
    FFFear_slow() {}
    //! constructor
    FFFear_slow( const Xmap<T>& xmap ) { init( xmap ); }
    //! constructor: shorthand for constructor+operator
    FFFear_slow( Xmap<T>& result, const NXmap<T>& srchval, const NXmap<T>& srchwgt, const Xmap<T>& xmap, const NX_operator& nxop ) { init( xmap ); (*this)( result, srchval, srchwgt, nxop ); }
    //! initialiser: initialise with the given target Xmap
    void init( const Xmap<T>& xmap ) { xmp = &xmap; }
    bool operator() ( Xmap<T>& result, const NXmap<T>& srchval, const NXmap<T>& srchwgt, const NX_operator& nxop ) const;
    bool operator() ( Xmap<T>& result, const NXmap<T>& srchval, const NXmap<T>& srchwgt ) const;  //!< \deprecated
  private:
    const Xmap<T>* xmp;
  };


  //! FFT-based fffear implementation
  /*! \ingroup g_fffear
    This implementation is currently unoptimised, but much faster then
    the simple implementation. */
  template<class T> class FFFear_fft : public FFFear_base<T> {
  public:
    enum FFTtype { Default, Normal, Sparse };  //!< FFT backend selection
    //! constructor
    FFFear_fft() {}
    //! constructor
    FFFear_fft( const Xmap<T>& xmap ) { init( xmap ); }
    //! constructor: shorthand for constructor+operator
    FFFear_fft( Xmap<T>& result, const NXmap<T>& srchval, const NXmap<T>& srchwgt, const Xmap<T>& xmap, const NX_operator& nxop ) { init( xmap ); (*this)( result, srchval, srchwgt, nxop ); }
    //! initialiser: initialise with the given target Xmap
    void init( const Xmap<T>& xmap );
    void set_fft_type( FFTtype type );       //! option: fft backend
    void set_resolution( Resolution reso );  //! option: resolution cutoff
    bool operator() ( Xmap<T>& result, const NXmap<T>& srchval, const NXmap<T>& srchwgt, const NX_operator& nxop ) const;  //!< search for given target
    bool operator() ( Xmap<T>& result, const NXmap<T>& srchval, const NXmap<T>& srchwgt, const RTop_orth& rtop ) const;  //!< search for given target
    bool operator() ( Xmap<T>& result, const NXmap<T>& srchval, const NXmap<T>& srchwgt ) const;  //!< \deprecated
  private:
    Cell cell;
    FFTmap_p1 rho1, rho2;
    FFTtype ffttype;
  };


} // namespace clipper

#endif
