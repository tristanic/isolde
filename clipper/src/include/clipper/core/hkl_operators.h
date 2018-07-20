/*! \file lib/hkl_operators.h
    HKL_data operators for the clipper libraries
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


#ifndef CLIPPER_HKL_OPERATORS
#define CLIPPER_HKL_OPERATORS

#include "hkl_datatypes.h"
#include "../imex.h"

namespace clipper
{

  namespace datatypes
  {

    // logical operators
    CLIPPER_IMEX HKL_data<Flag_bool> operator &( const HKL_data_base& d1, const HKL_data_base& d2 );
	CLIPPER_IMEX HKL_data<Flag_bool> operator |( const HKL_data_base& d1, const HKL_data_base& d2 );
	CLIPPER_IMEX HKL_data<Flag_bool> operator ^( const HKL_data_base& d1, const HKL_data_base& d2 );
	CLIPPER_IMEX HKL_data<Flag_bool> operator !( const HKL_data_base& d1 );

	CLIPPER_IMEX HKL_data<Flag_bool> operator ==( const HKL_data<Flag>& d1, const int& n );
	CLIPPER_IMEX HKL_data<Flag_bool> operator !=( const HKL_data<Flag>& d1, const int& n );
	CLIPPER_IMEX HKL_data<Flag_bool> operator >=( const HKL_data<Flag>& d1, const int& n );
	CLIPPER_IMEX HKL_data<Flag_bool> operator <=( const HKL_data<Flag>& d1, const int& n );
	CLIPPER_IMEX HKL_data<Flag_bool> operator >( const HKL_data<Flag>& d1, const int& n );
	CLIPPER_IMEX HKL_data<Flag_bool> operator <( const HKL_data<Flag>& d1, const int& n );

    // individual data operators
    template<class dtype> F_phi<dtype> operator +( const F_phi<dtype>& d1, const F_phi<dtype>& d2 );
    template<class dtype> F_phi<dtype> operator -( const F_phi<dtype>& d1, const F_phi<dtype>& d2 );
    template<class dtype> F_phi<dtype> operator -( const F_phi<dtype>& d1 );
    template<class dtype> ABCD<dtype> operator +( const ABCD<dtype>& d1, const ABCD<dtype>& d2 );

    // data list operators
    template<class dtype> HKL_data<F_phi<dtype> > operator +( const HKL_data<F_phi<dtype> >& d1, const HKL_data<F_phi<dtype> >& d2 );
    template<class dtype> HKL_data<F_phi<dtype> > operator -( const HKL_data<F_phi<dtype> >& d1, const HKL_data<F_phi<dtype> >& d2 );
    template<class dtype> HKL_data<F_phi<dtype> > operator *( const HKL_data<F_phi<dtype> >& d1, const ftype& s );
    template<class dtype> HKL_data<F_phi<dtype> > operator -( const HKL_data<F_phi<dtype> >& d1 );
    template<class dtype> HKL_data<ABCD<dtype> > operator +( const HKL_data<ABCD<dtype> >& d1, const HKL_data<ABCD<dtype> >& d2 );
    template<class dtype> HKL_data<F_phi<dtype> > operator *( const ftype& s, const HKL_data<F_phi<dtype> >& d1 ) { return d1*s; }

  } // namespace datatypes



  //! Log phase probability distribution object
  /*! This object is used to store and manipulate phase
    log-probability distributions. Centrics are handled by two values
    on the phase circle, acentrics by a list of values. The values can
    be indexed like and array. The phase() function returns the phase
    corresponding to the given array index. Conversion to and from
    Hendrickson-Lattman coefficients is provided. The object is
    templatised on the sampling of the phase circle. */
  template<int N> class LogPhaseProb {
  public:
    //! constructor: from HKL class
    LogPhaseProb( const HKL_class& hkl_class );
    //! set HL coeffs
    template<class dtype> void set_abcd( const datatypes::ABCD<dtype>& abcd );
    //! get HL coeffs
    template<class dtype> void get_abcd( datatypes::ABCD<dtype>& abcd ) const;
    //! set phi/fom
    template<class dtype> void set_phi_fom( const datatypes::Phi_fom<dtype>& phifom );
    //! get phi/fom
    template<class dtype> void get_phi_fom( datatypes::Phi_fom<dtype>& phifom ) const;
    //! get log probability
    const ftype& operator[] ( const int& p ) const { return q[p]; }
    //! set log probability
    ftype& operator[] ( const int& p ) { return q[p]; }
    //! return phase associated with index
    ftype phase( const int& p ) const
      { return Util::twopi()*ftype(p*pinc+pmin)/ftype(N); }
    int size() const { return q.size(); }  //!< return num. of phases
    static int sampling() { return N; }    //!< return phase sampling
  private:
    int pmin, pinc;
    std::vector<ftype> q;
  };


} // namespace clipper

#endif
