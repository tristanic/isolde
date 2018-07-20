/*! \file lib/resol_targetfn.h
    Header file for resolution function generator
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


#ifndef CLIPPER_RESOL_TARGETFN
#define CLIPPER_RESOL_TARGETFN

#include "resol_basisfn.h"
#include "hkl_datatypes.h"

namespace clipper {


  //! simple mean |F|<sup>n</sup> target
  /*! This class implements the target function for calculating mean
    |F|<sup>n</sup> as a function of position in reciprocal space. It
    includes the appropriate multiplicity correction, and so can be
    applied to any type with an 'f' member with the same dimensions as
    an |F| or |U| (or an uncorrected |E|).

    \Note This function should not be used to scale F's to E's.
    See TargetFn_scaleEsq. */
  template<class T> class TargetFn_meanFnth : public TargetFn_base
  {
  public:
    //! constructor: takes the datalist against which to calc target, and power
    TargetFn_meanFnth( const HKL_data<T>& hkl_data_, const ftype& n );
    //! return the value and derivatives of the target function
    Rderiv rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const;
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return QUADRATIC; }
  private:
    ftype power;
    const HKL_data<T>* hkl_data;
  };


  //! simple mean |E|<sup>n</sup> target
  /*! This class implements the target function for calculating mean
    |E|<sup>n</sup> as a function of position in reciprocal space. It
    includes the appropriate multiplicity correction, and so can be
    applied to any type with an 'E' member with the same dimensions as
    an |E| (or corrected |F| or |U|).

    \Note This function should not be used to scale F's to E's.
    See TargetFn_scaleEsq. */
  template<class T> class TargetFn_meanEnth : public TargetFn_base
  {
  public:
    //! constructor: takes the datalist against which to calc target, and power
    TargetFn_meanEnth( const HKL_data<T>& hkl_data_, const ftype& n );
    //! return the value and derivatives of the target function
    Rderiv rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const;
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return QUADRATIC; }
  private:
    ftype power;
    const HKL_data<T>* hkl_data;
  };


  //! |F|<sup>2</sup> scaling target
  /*! This class implements the target function for calculating the
    scale factor to scale one set of F's to another. The resulting
    scale is the square of the factor that scales the first set of
    data to match the second. */
  template<class T1, class T2> class TargetFn_scaleF1F2 : public TargetFn_base
  {
  public:
    //! constructor: takes the datalist against which to calc target
    TargetFn_scaleF1F2( const HKL_data<T1>& hkl_data1_, const HKL_data<T2>& hkl_data2_ );
    //! return the value and derivatives of the target function
    Rderiv rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const;
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return QUADRATIC; }
  private:
    const HKL_data<T1>* hkl_data1;
    const HKL_data<T2>* hkl_data2;
  };


  //! log |F|<sup>2</sup> scaling target
  /*! This class implements the target function for calculating the
    scale factor to scale the weighted log of one set of F's to
    another. The resulting scale is the square of the factor that
    scales the first set of data to match the second. The log scaling
    target is used in conjunction with the log-Gaussian basis
    functions for a fast and robust approximation to iso/aniso
    Gaussian scaling.
  */
  template<class T1, class T2> class TargetFn_scaleLogF1F2 : public TargetFn_base
  {
  public:
    //! constructor: takes the datalist against which to calc target
    TargetFn_scaleLogF1F2( const HKL_data<T1>& hkl_data1_, const HKL_data<T2>& hkl_data2_ );
    //! return the value and derivatives of the target function
    Rderiv rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const;
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return QUADRATIC; }
  private:
    const HKL_data<T1>* hkl_data1;
    const HKL_data<T2>* hkl_data2;
  };


  /*! This class implements the target function for calculating the
    scale factor to scale one set of I's to another. The resulting
    scale is the square of the factor that scales the first set of
    data to match the second. */
  template<class T1, class T2> class TargetFn_scaleI1I2 : public TargetFn_base
  {
  public:
    //! constructor: takes the datalist against which to calc target
    TargetFn_scaleI1I2( const HKL_data<T1>& hkl_data1_, const HKL_data<T2>& hkl_data2_ );
    //! return the value and derivatives of the target function
    Rderiv rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const;
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return QUADRATIC; }
  private:
    const HKL_data<T1>* hkl_data1;
    const HKL_data<T2>* hkl_data2;
  };


  //! log |I| scaling target
  /*! This class implements the target function for calculating the
    scale factor to scale the weighted log of one set of I's to
    another. The resulting scale is the square of the factor that
    scales the first set of data to match the second. The log scaling
    target is used in conjunction with the log-Gaussian basis
    functions for a fast and robust approximation to iso/aniso
    Gaussian scaling.
  */
  template<class T1, class T2> class TargetFn_scaleLogI1I2 : public TargetFn_base
  {
  public:
    //! constructor: takes the datalist against which to calc target
    TargetFn_scaleLogI1I2( const HKL_data<T1>& hkl_data1_, const HKL_data<T2>& hkl_data2_ );
    //! return the value and derivatives of the target function
    Rderiv rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const;
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return QUADRATIC; }
  private:
    const HKL_data<T1>* hkl_data1;
    const HKL_data<T2>* hkl_data2;
  };


  //! |E|<sup>2</sup> scaling target
  /*! This class implements the target function for calculating the
    scale factor to normalise to <|E|<sup>2</sup>> = 1. Note that this
    is not the same as dividing by <|E|<sup>2</sup>>, except in a few
    special cases, e.g. a simple resolution bins calculation. The
    resulting targen function is the square of the value by which |E|
    should be multiplied to acheive the correct normalisation. */
  template<class T> class TargetFn_scaleEsq : public TargetFn_base
  {
  public:
    //! constructor: takes the datalist against which to calc target
    TargetFn_scaleEsq( const HKL_data<T>& hkl_data_ );
    //! return the value and derivatives of the target function
    Rderiv rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const;
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return QUADRATIC; }
  private:
    const HKL_data<T>* hkl_data;
  };


  //! \deprecated simple sigma_a target function
  /*! This class implements the target function for calculating sigma_a.
    Required is a datalist containing Eo, Ec.

    It actually refines omegaa = sigmaa/(1-sigmaa^2). This has better
    proerties for refinement. To get sigmaa use
    \code sigmaa = ( sqrt( 4*omegaa^2 + 1 ) - 1 ) / ( 2*omegaa ) \endcode
    This is available as a static function:
    \code sigmaa = targetfn.sigmaa( omegaa ) \endcode

    This version simplifies terms in |Eo|^2 and |Ec|^2 which should
    average out to 1 if the normalisation scheme is consistent with
    the sigmaa calc.

    Convergence is good for calculations using the 'binner' basis
    function, however the smooth basis function have convergence
    problems. This is still under investigation.
  */
  template<class T> class TargetFn_sigmaa_omegaa : public TargetFn_base
  {
  public:
    //! constructor: takes the datalist against which to calc target
    TargetFn_sigmaa_omegaa( const HKL_data<T>& eo, const HKL_data<T>& ec );
    //! return the value and derivatives of the target function
    Rderiv rderiv( const HKL_info::HKL_reference_index& ih, const ftype& omegaa ) const;
    //! convert omegaa to sigmaa
    static ftype sigmaa( const ftype& omegaa )
    {
      ftype omeg = (omegaa>0.05) ? omegaa : (0.05*exp(omegaa/0.05-1.0));
      return 0.5 * ( sqrt( 4.0*omeg*omeg + 1.0 ) - 1.0 ) / omeg;
    }
  private:
    const HKL_data<T>* eo_data;
    const HKL_data<T>* ec_data;
  };


  //! \deprecated simple sigma_a target function
  /*! \par Warning: Convergence of this basis-function can be
    unreliable under some circumstances. Use
    clipper::TargetFn_sigmaa_omegaa instead, except for development
    purposes.

    This class implements the target function for calculating sigma_a.
    Required is a datalist containing Eo, Ec.

    This version simplifies terms in |Eo|^2 and |Ec|^2 which should
    average out to 1 if the normalisation scheme is consistent with
    the sigmaa calc.
  */
  template<class T> class TargetFn_sigmaa : public TargetFn_base
  {
  public:
    //! constructor: takes the datalist against which to calc target
    TargetFn_sigmaa( const HKL_data<T>& eo, const HKL_data<T>& ec );
    //! return the value and derivatives of the target function
    Rderiv rderiv( const HKL_info::HKL_reference_index& ih, const ftype& sigmaa0 ) const;
    //! convert function to sigmaa
    static ftype sigmaa( const ftype& sigm ) { return sigm; }
  private:
    const HKL_data<T>* eo_data;
    const HKL_data<T>* ec_data;
  };



  // implementations for template functions

  // mean F^nth

  template<class T> TargetFn_meanFnth<T>::TargetFn_meanFnth( const HKL_data<T>& hkl_data_, const ftype& n )
  {
    power = n;
    hkl_data = &hkl_data_;
  }

  template<class T> TargetFn_base::Rderiv TargetFn_meanFnth<T>::rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const
  {
    Rderiv result;
    const HKL_data<T>& data = *hkl_data;
    if ( !data[ih].missing() ) {
      ftype d = fh - pow( ftype(data[ih].f()) / sqrt(ih.hkl_class().epsilon()),
			  power );
      result.r = d * d;
      result.dr = 2.0 * d;
      result.dr2 = 2.0;
    } else {
      result.r = result.dr = result.dr2 = 0.0;
    }
    return result;
  }

  // mean E^nth

  template<class T> TargetFn_meanEnth<T>::TargetFn_meanEnth( const HKL_data<T>& hkl_data_, const ftype& n )
  {
    power = n;
    hkl_data = &hkl_data_;
  }

  template<class T> TargetFn_base::Rderiv TargetFn_meanEnth<T>::rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const
  {
    Rderiv result;
    const HKL_data<T>& data = *hkl_data;
    if ( !data[ih].missing() ) {
      ftype d = fh - pow( ftype(data[ih].E()), power );
      result.r = d * d;
      result.dr = 2.0 * d;
      result.dr2 = 2.0;
    } else {
      result.r = result.dr = result.dr2 = 0.0;
    }
    return result;
  }

  // F1-F2 scaling

  template<class T1, class T2> TargetFn_scaleF1F2<T1,T2>::TargetFn_scaleF1F2( const HKL_data<T1>& hkl_data1_, const HKL_data<T2>& hkl_data2_ )
  {
    hkl_data1 = &hkl_data1_;
    hkl_data2 = &hkl_data2_;
  }

  template<class T1, class T2> TargetFn_base::Rderiv TargetFn_scaleF1F2<T1,T2>::rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const
  {
    Rderiv result;
    const T1& ft1 = (*hkl_data1)[ih];
    const T2& ft2 = (*hkl_data2)[ih];
    if ( !ft1.missing() && !ft2.missing() ) {
      const ftype eps = ih.hkl_class().epsilon();
      const ftype f1 = pow( ft1.f(), 2 ) / eps;
      const ftype f2 = pow( ft2.f(), 2 ) / eps;
      const ftype d = fh*f1 - f2;
      result.r = d * d / f1;
      result.dr = 2.0 * d;
      result.dr2 = 2.0 * f1;
    } else {
      result.r = result.dr = result.dr2 = 0.0;
    }
    return result;
  }

  // Log F1-F2 scaling

  template<class T1, class T2> TargetFn_scaleLogF1F2<T1,T2>::TargetFn_scaleLogF1F2( const HKL_data<T1>& hkl_data1_, const HKL_data<T2>& hkl_data2_ )
  {
    hkl_data1 = &hkl_data1_;
    hkl_data2 = &hkl_data2_;
  }

  template<class T1, class T2> TargetFn_base::Rderiv TargetFn_scaleLogF1F2<T1,T2>::rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const
  {
    Rderiv result;
    result.r = result.dr = result.dr2 = 0.0;
    const T1& ft1 = (*hkl_data1)[ih];
    const T2& ft2 = (*hkl_data2)[ih];
    if ( !ft1.missing() && !ft2.missing() )
      if ( ft1.f() > 1.0e-6 && ft2.f() > 1.0e-6 ) {
	const ftype eps = ih.hkl_class().epsilon();
	const ftype f1 = pow( ft1.f(), 2 ) / eps;
	const ftype f2 = pow( ft2.f(), 2 ) / eps;
	const ftype w = 1.0;  // f1 * f2;
	const ftype d = fh + log(f1) - log(f2);
	result.r   =       w * d * d;
	result.dr  = 2.0 * w * d;
	result.dr2 = 2.0 * w;
    }
    return result;
  }

  // I1-I2 scaling

  template<class T1, class T2> TargetFn_scaleI1I2<T1,T2>::TargetFn_scaleI1I2( const HKL_data<T1>& hkl_data1_, const HKL_data<T2>& hkl_data2_ )
  {
    hkl_data1 = &hkl_data1_;
    hkl_data2 = &hkl_data2_;
  }

  template<class T1, class T2> TargetFn_base::Rderiv TargetFn_scaleI1I2<T1,T2>::rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const
  {
    Rderiv result;
    const T1& ft1 = (*hkl_data1)[ih];
    const T2& ft2 = (*hkl_data2)[ih];
    if ( !ft1.missing() && !ft2.missing() ) {
      const ftype eps = ih.hkl_class().epsilon();
      const ftype f1 = ft1.I() / eps;
      const ftype f2 = ft2.I() / eps;
      const ftype d = fh*f1 - f2;
      result.r = d * d / f1;
      result.dr = 2.0 * d;
      result.dr2 = 2.0 * f1;
    } else {
      result.r = result.dr = result.dr2 = 0.0;
    }
    return result;
  }

  // Log I1-I2 scaling

  template<class T1, class T2> TargetFn_scaleLogI1I2<T1,T2>::TargetFn_scaleLogI1I2( const HKL_data<T1>& hkl_data1_, const HKL_data<T2>& hkl_data2_ )
  {
    hkl_data1 = &hkl_data1_;
    hkl_data2 = &hkl_data2_;
  }

  template<class T1, class T2> TargetFn_base::Rderiv TargetFn_scaleLogI1I2<T1,T2>::rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const
  {
    Rderiv result;
    result.r = result.dr = result.dr2 = 0.0;
    const T1& ft1 = (*hkl_data1)[ih];
    const T2& ft2 = (*hkl_data2)[ih];
    if ( !ft1.missing() && !ft2.missing() )
      if ( ft1.I() > 1.0e-6 && ft2.I() > 1.0e-6 ) {
	const ftype eps = ih.hkl_class().epsilon();
	const ftype f1 = ft1.I() / eps;
	const ftype f2 = ft2.I() / eps;
	const ftype w = 1.0;  // f1 * f2;
	const ftype d = fh + log(f1) - log(f2);
	result.r   =       w * d * d;
	result.dr  = 2.0 * w * d;
	result.dr2 = 2.0 * w;
    }
    return result;
  }

  // E^2 scaling

  template<class T> TargetFn_scaleEsq<T>::TargetFn_scaleEsq( const HKL_data<T>& hkl_data_ )
  {
    hkl_data = &hkl_data_;
  }

  template<class T> TargetFn_base::Rderiv TargetFn_scaleEsq<T>::rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const
  {
    Rderiv result;
    const HKL_data<T>& data = *hkl_data;
    const ftype two(2.0);
    if ( !data[ih].missing() ) {
      ftype fsq = pow( ftype(data[ih].E()), two );
      ftype d = fsq * fh - 1.0;
      result.r = d * d / fsq;
      result.dr = two * d;
      result.dr2 = two * fsq;
    } else {
      result.r = result.dr = result.dr2 = 0.0;
    }
    return result;
  }


  // sigmaa (omegaa)

  template<class T> TargetFn_sigmaa_omegaa<T>::TargetFn_sigmaa_omegaa( const HKL_data<T>& eo, const HKL_data<T>& ec )
  {
    eo_data = &eo;
    ec_data = &ec;
  }

  template<class T> TargetFn_base::Rderiv TargetFn_sigmaa_omegaa<T>::rderiv( const HKL_info::HKL_reference_index& ih, const ftype& omegaa ) const
  {
    Rderiv result;
    
    const HKL_data<T>& eodata = *eo_data;
    const HKL_data<T>& ecdata = *ec_data;
    if ( eodata[ih].missing() || ecdata[ih].missing()  ) {
      result.r = result.dr = result.dr2 = 0.0;
    } else {
      ftype eo = eodata[ih].E();
      ftype ec = ecdata[ih].E();
      ftype omeg = (omegaa>0.05) ? omegaa : (0.05*exp(omegaa/0.05-1.0));
      ftype sigmaa = 0.5*(sqrt(4*omeg*omeg+1)-1)/omeg;
      ftype dx = 2.0 * eo * ec;
      ftype x = dx * omeg;
      ftype t0 = 1.0/(1-sigmaa*sigmaa) + 0.5*log((1-sigmaa*sigmaa));
      ftype t1 = sigmaa;
      ftype t2 = pow(1-sigmaa*sigmaa,2)/(1+sigmaa*sigmaa);
      if ( ih.hkl_class().centric() ) {
	result.r = 1.0*t0 - log( cosh( x/2 ) );
	result.dr = 1.0*t1 - dx*0.5*tanh( x/2 );
	result.dr2 = 1.0*t2 - dx*dx*0.25*(1.0 - pow(tanh(x/2),2) );
      } else {
	result.r = 2.0*t0 - Util::sim_integ( x );
	result.dr = 2.0*t1 - dx*Util::sim( x );
	result.dr2 = 2.0*t2 - dx*dx*Util::sim_deriv( x );
      }
      if ( omegaa < 0.05 ) {
	ftype dy = exp( omegaa/0.05 ) / exp( 1.0 );
	ftype dy2 = exp( omegaa/0.05 ) / ( 0.05*exp( 1.0 ) );
	result.dr2 = result.dr*dy2 + result.dr2*dy*dy;
	result.dr = result.dr*dy;
      }
    }
    return result;
  }


  // sigmaa (norm)

  template<class T> TargetFn_sigmaa<T>::TargetFn_sigmaa( const HKL_data<T>& eo, const HKL_data<T>& ec )
  {
    eo_data = &eo;
    ec_data = &ec;
  }

  template<class T> TargetFn_base::Rderiv TargetFn_sigmaa<T>::rderiv( const HKL_info::HKL_reference_index& ih, const ftype& sigmaa0 ) const
  {
    Rderiv result;
    
    const HKL_data<T>& eodata = *eo_data;
    const HKL_data<T>& ecdata = *ec_data;
    if ( eodata[ih].missing() || ecdata[ih].missing()  ) {
      result.r = result.dr = result.dr2 = 0.0;
    } else {
      ftype eo = eodata[ih].E();
      ftype ec = ecdata[ih].E();
      ftype sigmaa = (sigmaa0>0.99)?(0.99):((sigmaa0<0.01)?0.01:sigmaa0); 
      ftype dx = 2.0 * eo * ec;
      ftype x = dx * sigmaa/(1-sigmaa*sigmaa);
      ftype t0 = 1.0/(1-sigmaa*sigmaa) + 0.5*log((1-sigmaa*sigmaa));
      ftype t1 = sigmaa;
      ftype t2 = pow(1-sigmaa*sigmaa,2)/(1+sigmaa*sigmaa);
      if ( ih.hkl_class().centric() ) {
	result.r = 1.0*t0 - log( cosh( x/2 ) );
	result.dr = 1.0*t1 - dx*0.5*tanh( x/2 );
	result.dr2 = 1.0*t2 - dx*dx*0.25*(1.0 - pow(tanh(x/2),2) );
      } else {
	result.r = 2.0*t0 - Util::sim_integ( x );
	result.dr = 2.0*t1 - dx*Util::sim( x );
	result.dr2 = 2.0*t2 - dx*dx*Util::sim_deriv( x );
      }
      ftype ds = (1+sigmaa*sigmaa)/pow(1-sigmaa*sigmaa,2);
      ftype ds2 = 2*sigmaa*(3+sigmaa*sigmaa)/pow(1-sigmaa*sigmaa,3);
      result.dr2 = result.dr*ds2 + result.dr2*ds*ds;
      result.dr = result.dr*ds;
    }
    return result;
  }


} // namespace clipper

#endif
