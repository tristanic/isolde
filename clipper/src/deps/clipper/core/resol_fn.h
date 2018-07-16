/*! \file lib/resol_fn.h
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


#ifndef CLIPPER_RESOL_FN
#define CLIPPER_RESOL_FN

#include "hkl_data.h"
#include "../imex.h"

namespace clipper {

  //! abstract base class for resolution function basis functions
  /*! A basis function must be able to return its value and
    derivatives for any given HKL.

    Optionally, performance can be improved by returning a flag to
    indicate if the value of the basis function for a given reflection
    is linearly dependent on the values of the parameter, and a value
    indicating whether the curvature matrix takes an N-diagonal form.

    \b NOTE: for performance reasons the derivatives are returned as a
    reference to an internal object, so if you store a reference to
    the result (which is also good for performance, it will be
    overwritten on the next call. If this matters, store a copy rather
    than a reference. */
  class CLIPPER_IMEX BasisFn_base
  {
  public:
    //! enumeration of function types: optionally used to improve convergence
    enum FNtype { GENERAL, LINEAR };

    //! object holding the basis function and its first two derivatives
    class Fderiv
    {
    public:
      ftype f;                //!< value of function
      std::vector<ftype> df;  //!< first derivative vector w.r.t params
      Matrix<> df2;      //!< second derivative matrix w.r.t params
      Fderiv() {}             //!< null constructor
      Fderiv(const int& np) : df(np,0.0), df2(np,np,0.0) {}  //<! constructor
    };

    //! null constructor
    BasisFn_base() {}
    //! constructor: takes number of parameters
    BasisFn_base( const int& np ) : np_(np), result_(np) {}
    //! the number of parameters of this basis function
    const int& num_params() const { return np_; }
    //! the value of the resolution function
    virtual ftype f( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const { return fderiv( hkl, cell, params ).f; }
    //! the value of the resolution function and its first two derivatives
    virtual const Fderiv& fderiv( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const = 0;
    //! the type of the function: optionally used to improve convergence
    virtual FNtype type() const;
    //! number of non-zero diagonals in the upper triangle of the curvatures
    virtual int num_diagonals() const;
  protected:
    //! provide write access to result for subclasses
    Fderiv& result() const { return result_; }
    virtual ~BasisFn_base() {}  //!< destructor

  private:
    int np_;                 //!< number of params
    mutable Fderiv result_;  //!< internal cache for result
  };

  //! abstract base class for least-squares resolution function target functions
  /*! A target function must be able to return its value given the
    value of the basis function for all HKL, and its derivative with
    respect the values of the basis function for all HKL.

    Optionally, performance can be improved by returning a flag to
    indicate if the target function is quadratic. */
  class CLIPPER_IMEX TargetFn_base
  {
  public:
    //! object holding the residual function and first two derivatives
    class Rderiv
    {
    public:
      ftype r;    //!< the value of the function
      ftype dr;   //!< first derivative w.r.t basis fn
      ftype dr2;  //!< second derivative w.r.t basis fn
    };

    //! enumeration of function types: optionally used to improve convergence
    enum FNtype { GENERAL, QUADRATIC };
    //! return the value and derivatives of the target function
    /*! If the value of f(h) is invalid, rderiv.r should be set to NaN */
    virtual Rderiv rderiv( const HKL_info::HKL_reference_index& ih, const ftype& fh ) const = 0;
    //! the type of the function: optionally used to improve convergence
    virtual FNtype type() const { return GENERAL; }
    virtual ~TargetFn_base() {}  //!< destructor

    //! test that the residuals, gradients, and curvatures are consistent
    void debug( const HKL_info& hkl_info ) const;
  };


  //! 2nd order resolution function evaluator
  /*! This is an automatic evaluator for arbitrary functions of HKL,
    most commonly used for evaluating a function of resolution (such a
    mean F^2 or sigmaa), although more general tasks including local
    scaling of reflections and anisotropic functions can also be
    handled. This form is for target functions which approach zero
    quadratically, e.g. least-squares targets.

    This version implements a naive Newton-Raphson minimiser, which
    only uses the gradient and curvature of the target function,
    ignoring its value. It is ideal for quadratic targets with linear
    basis functions.

    To evaluate a resolution function, this class must be provided
    with two objects:
     - The basis function (and gradients), which describes the value
     of the function for any reflection given a set of paramters.
     - The target function (and derivatives), which is used to determine
     the values of the basis function parameters.

    For example, the following code may be used to calculate a smooth
    scaling function to scale one set of data to another using an
    anisotropic Gaussian scaling function:
    \code
    // make data object
    clipper::HKL_info hkls;
    clipper::HKL_data<clipper::data32::F_sigF> fsig1( hkls );
    clipper::HKL_data<clipper::data32::F_sigF> fsig2( hkls );
    // INITIALISE THEM HERE!

    // calculate the scaling function
    std::vector<clipper::ftype> params( 7, 0.0 );  // initial parameters
    clipper::BasisFn_aniso_gaussian basisfn;       // aniso gaussian
    clipper::TargetFn_scaleF1F2<clipper::data32::F_sigF,clipper::data32::F_sigF>targetfn( fsig1, fsig2 );
    clipper::ResolutionFn_nonlinear rfn( hkls, basisfn, targetfn, params );

    // now scale the data
    clipper::HKL_info::HKL_reference_index ih;
    for ( ih = hkls.first(); !ih.last(); ih.next() )
      fsig1[ih].scale( sqrt( rfn.f(ih) ) );
    \endcode

    The most useful isotropic resolution function is the
    BasisFn_spline, since it is linear and provides a good fit to most
    data.
   */
  class CLIPPER_IMEX ResolutionFn
  {
  public:
    //! constructor: need reflections, basis fn and target fn.
    ResolutionFn( const HKL_info& hkl_info, const BasisFn_base& basisfn, const TargetFn_base& targetfn, const std::vector<ftype>& params, const ftype damp = 0.0, const bool debug = false );
    //! return the value of the basis function with the current paramters
    inline ftype f( const HKL_info::HKL_reference_index& ih ) const { return basisfn_->f( ih.hkl(), cell_, params_ ); }
    //! return the values of the parameters
    const std::vector<ftype>& params() const;

    //! print the target, gradient, and curvatures with respect to the params
    void debug() const;

  protected:
    const HKL_info* hkl_info_;       //!< reflection list
    const TargetFn_base* targetfn_;  //!< target function
    const BasisFn_base* basisfn_;    //!< basis function
    std::vector<ftype> params_;      //!< basis function parameters
    Cell cell_;                      //!< cell

    //! calculate derivatives of target wrt params \internal
    void calc_derivs( const std::vector<ftype>& params, ftype& r, std::vector<ftype>& drdp, Matrix<>& drdp2 ) const;

    //! null constructor
    ResolutionFn() {}
  };


  //! 2nd order resolution function evaluator
  /*! This is an automatic evaluator for arbitrary functions of HKL,
    most commonly used for evaluating a function of resolution (such a
    mean F^2 or sigmaa), although more general tasks including local
    scaling of reflections and anisotropic functions can also be
    handled. This form is for target functions which approach zero
    quadratically, e.g. least-squares targets.

    \note This version implements a minimiser which uses both
    Newton-Raphson and gradient steps depending on the situation. It
    can be used for non-quadratic targets or non-linear basis
    functions.

    To evaluate a resolution function, this class must be provided
    with two objects:
     - The basis function (and gradients), which describes the value
     of the function for any reflection given a set of paramters.
     - The target function (and derivatives), which is used to determine
     the values of the basis function parameters.
   */
  class CLIPPER_IMEX ResolutionFn_nonlinear: public ResolutionFn
  {
  public:
    //! constructor: need reflections, basis fn and target fn.
    ResolutionFn_nonlinear( const HKL_info& hkl_info, const BasisFn_base& basisfn, const TargetFn_base& targetfn, const std::vector<ftype>& params, const ftype damp = 0.0, const bool debug = false );
  };


} // namespace clipper

#endif
