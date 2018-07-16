/* resol_fn.cpp: implementation file for resolution functions */
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
/* This code is derived from the 'dm' source code */


#include "resol_fn.h"


namespace clipper {


void TargetFn_base::debug( const HKL_info& hkl_info ) const
{
  Rderiv r0, r1, r2;
  ftype d = 0.0001;
  for ( HKL_info::HKL_reference_index ih = hkl_info.first(); !ih.last(); ih.next() ) {
    if ( ih.index()%100 == 0 )
      for ( ftype f = 0.1; f < 10; f*=3.15 ) {
	r0 = rderiv( ih, f - d );
	r1 = rderiv( ih, f );
	r2 = rderiv( ih, f + d );
	std::cout << ih.hkl().format() << " " << f << " val " << r1.r << " grad " << r1.dr << " " << (r2.r-r0.r)/(2*d) << " curv " << r1.dr2 << " " << (r2.dr-r0.dr)/(2*d) << " " << (r2.r+r0.r-2*r1.r)/(d*d) << "\n";
      }
  }
}


/*! Defaults to GENERAL, which will always work. If the basis function
  is linearly dependent on the parameters, override this with a
  function returning LINEAR for improved performance. See the provided
  basis functions for examples.
  \return The function type enumeration. */
BasisFn_base::FNtype BasisFn_base::type() const
{ return GENERAL; }


/*! Defaults to 0, which will always work. If the basis function has
  compact support among the parameters, i.e. the value for any HKL
  depends only on a few parameters, then set this to the number of
  non-zero diagonals in the upper triangle of the matrix, i.e.  1 for
  a diagonal matrix, 2 for a tri-diagonal matrix etc.
  \return The number of non-zero upper diagonals, or zero for full-matrix. */
int BasisFn_base::num_diagonals() const
{ return 0; }


/*! The target and its derivatives are calculated by the chain rule
  and accumulated over all HKLs.
  \param params Current value of parameters.
  \param r The value of the target function.
  \param drdp The derivative of the target function.
  \param drdp2 The curvature of the target function. */
void ResolutionFn::calc_derivs( const std::vector<ftype>& params, ftype& r, std::vector<ftype>& drdp, Matrix<>& drdp2 ) const
{
  HKL_class cls;
  HKL_info::HKL_reference_index ih;
  TargetFn_base::Rderiv rderiv;
  ftype w;
  int j, k;
  int nparams = basisfn_->num_params();
  int ndiag = basisfn_->num_diagonals() - 1;

  // zero the results
  r = 0.0;
  for ( j = 0; j < nparams; j++ ) { 
    drdp[j] = 0.0;
    for ( k = 0; k < nparams; k++ )
      drdp2(j,k) = 0.0;
  }

  // calc dr/dp vector and d2r/dp2 matrix
  for ( ih = hkl_info_->first(); !ih.last(); ih.next() ) {
    cls = ih.hkl_class();
    rderiv = targetfn_->rderiv( ih, f( ih ) );
    const BasisFn_base::Fderiv& fderiv =
      basisfn_->fderiv( ih.hkl(), cell_, params );
    const std::vector<ftype>& dfdp = fderiv.df;
    const Matrix<>& dfdp2 = fderiv.df2;

    // weight for this term
    w = 2.0 / cls.epsilonc();
    // residual, gradient
    r += w * rderiv.r;
    for ( j = 0; j < nparams; j++ )
      drdp[j] += w * rderiv.dr * dfdp[j];
    // curvature: optimised for special cases
    if ( ndiag < 0 ) {
      // general case
      for ( j = 0; j < nparams; j++ )
	for ( k = 0; k <  nparams; k++ )
	  drdp2(j,k) +=
	    w * ( rderiv.dr2*dfdp[j]*dfdp[k] + rderiv.dr*dfdp2(j,k) );
    } else {
      // curvature matrix is (2n-1)-diagonal
      for ( j = 0; j < nparams; j++ )
	for ( k = Util::max(j-ndiag,0); k <= Util::min(j+ndiag,nparams-1); k++ )
	  drdp2(j,k) += w*( rderiv.dr2*dfdp[j]*dfdp[k] + rderiv.dr*dfdp2(j,k) );
    }
  }
}


/*! The constructor performs the full minimisation calculation.
  \param hkl_info HKL_info object which provides the reflection list.
  \param basisfn The basis function used to describe the desired property.
  \param targetfn The target function to be minimised.
  \param params Initial values for the function parameters.
  \param damp_ If > 0.0, shifts are fdamped during early cycles to help
  convergence with difficult bases/target conbinations. */
ResolutionFn::ResolutionFn( const HKL_info& hkl_info, const BasisFn_base& basisfn, const TargetFn_base& targetfn, const std::vector<ftype>& params, const ftype damp, const bool debug )
{
  ftype r, g0, g1, scale;
  int i, k, n;

  hkl_info_ = &hkl_info;
  basisfn_  = &basisfn;
  targetfn_ = &targetfn;
  params_ = params;
  cell_ = hkl_info.cell();

  int nparams = basisfn_->num_params();
  Matrix<> dfdp2( nparams, nparams );
  Matrix<> drdp2( nparams, nparams );
  std::vector<ftype> drdp( nparams ), dfdp( nparams );
  std::vector<ftype> shiftn( nparams );
  params_.resize( nparams );

  // loop for 20 cycles refining the params
  g0 = 1.0e25;
  for ( n = 0; n < 20; n++ ) {
    // calc target fn and derivs
    calc_derivs( params_, r, drdp, drdp2 );
    g1 = 0.0;
    for ( k = 0; k < nparams; k++ ) g1 += drdp[k]*drdp[k];
    g1 = sqrt(g1);

    // stop if gradient starts increasing
    if ( g1 < 1.0e-10 || g1 >= g0 ) break;
    g0 = g1;

    // solve for Newton-Raphson shift
    shiftn = drdp2.solve( drdp );

    // apply shift
    scale = (1.0+ftype(n)) / (1.0+ftype(n)+damp);
    for ( i = 0; i< nparams; i++ ) params_[i] -= scale*shiftn[i];

    if ( debug ) std::cout << " Resolution function cycle " << n << " " << g0 << " " << g1 << " " << scale << "\n";

    // stop if the target is quadratic
    if ( basisfn.type()  == BasisFn_base::LINEAR &&
	 targetfn.type() == TargetFn_base::QUADRATIC ) break;
  }
}


/*! \return The refined basis function parameters */
const std::vector<ftype>& ResolutionFn::params() const
{ return params_; }


/*! The constructor performs the full minimisation calculation.
  \param hkl_info HKL_info object which provides the reflection list.
  \param basisfn The basis function used to describe the desired property.
  \param targetfn The target function to be minimised.
  \param damp_ If > 0.0, shifts are fdamped during early cycles to help
  convergence with difficult bases/target conbinations */
ResolutionFn_nonlinear::ResolutionFn_nonlinear( const HKL_info& hkl_info, const BasisFn_base& basisfn, const TargetFn_base& targetfn, const std::vector<ftype>& params, const ftype damp, const bool debug )
{
  ftype r0, r1, w, scale, g, s, dotprod;
  HKL_class cls;
  HKL_info::HKL_reference_index ih;
  int i, j, k, n;

  hkl_info_ = &hkl_info;
  basisfn_  = &basisfn;
  targetfn_ = &targetfn;
  params_ = params;
  cell_ = hkl_info.cell();

  int nparams = basisfn_->num_params();
  Matrix<> dfdp2( nparams, nparams );
  Matrix<> drdp2( nparams, nparams );
  std::vector<ftype> drdp( nparams ), dfdp( nparams );
  std::vector<ftype> shiftn( nparams ), shiftg( nparams );
  std::vector<ftype> params_old( nparams );
  params_.resize( nparams );

  // loop for 20 cycles refining the params
  for ( n = 0; n < 20; n++ ) {
    params_old = params_;

    // calc target fn and derivs
    calc_derivs( params_, r0, drdp, drdp2 );

    // solve for Newton-Raphson shift
    shiftn = drdp2.solve( drdp );

    // get step sizes and check direction
    g = s = dotprod = 0.0;
    for ( k = 0; k < nparams; k++ ) {
      g += drdp[k]*drdp[k];
      s += shiftn[k]*shiftn[k];
      dotprod += drdp[k]*shiftn[k];
    }
    g = sqrt(g);
    s = sqrt(s);

    // make gradient shift to match NR shift
    for ( k = 0; k < nparams; k++ ) shiftg[k] = drdp[k] * s / g;

    // if NR shift opposes gradient, then ignore the NR shift
    if ( dotprod < 0.0 ) shiftn = shiftg;

    if ( debug ) {
      std::cout << "\nResolution function cycle: " << n << "\n";
      if ( dotprod < 0.0 ) std::cout << " Using scaled grad " << s / g << "\n";
      std::cout << " Gradient " << g << "\n";
      for ( j = 0; j < nparams; j++ ) {
	for ( k = 0; k <  nparams; k++ ) std::cout << " " << drdp2(j,k);
	std::cout << "\n";
      }
      for ( k = 0; k < nparams; k++ ) std::cout << " " << k << " " << params_[k] << " " << drdp[k] << " " << shiftn[k] << "\n";
    }

    // now try the step and if necessary reduce the shift
    scale = (1.0+ftype(n)) / (1.0+ftype(n)+damp);
    for ( j = 0; j < 20; j++ ) {
      for ( i = 0; i< nparams; i++ )
	params_[i] = params_old[i] - scale*shiftn[i];
      r1 = 0.0;
      for ( ih = hkl_info.first(); !ih.last(); ih.next() ) {
	cls = ih.hkl_class();
	w = 2.0 / cls.epsilonc();
	r1 += w * targetfn.rderiv( ih, f( ih ) ).r;
      }
      if ( Util::is_nan(r1) ) {
	scale *= 0.5;
      } else {
	scale *= 0.5;
	if ( (r1-r0) <= 0.0 ) break;
      }
    }
 
    if ( debug ) std::cout << " Resolution function cycle " << n << " " << r0 << " " << r1 << " " << scale << "\n";

    if ( fabs( r1 - r0 ) < 1.0e-10 ) break;
  }
}


/* ResolutionFn( const HKL_info& hkl_info, const BasisFn_base& basisfn, const TargetFn_base& targetfn, const std::vector<ftype>& params, const std::vector<ftype>& targetval, const std::vector<ftype>& targetwgt, const ftype damp = 0.0, const bool debug = false );
  The constructor performs the full minimisation calculation.
  Additional parameters allow the function parameters to be restrained.
  \param hkl_info HKL_info object which provides the reflection list.
  \param basisfn The basis function used to describe the desired property.
  \param targetfn The target function to be minimised.
  \param params Initial values for the function parameters.
  \param targetval Target values for the function parameters.
  \param targetwgt Target sigmas for the function parameters.
  \param damp_ If > 0.0, shifts are fdamped during early cycles to help
  convergence with difficult bases/target conbinations.
ResolutionFn::ResolutionFn( const HKL_info& hkl_info, const BasisFn_base& basisfn, const TargetFn_base& targetfn, const std::vector<ftype>& params, const std::vector<ftype>& targetval, const std::vector<ftype>& targetwgt, const ftype damp, const bool debug )
{
  ftype r, g0, g1, scale;
  int i, k, n;

  hkl_info_ = &hkl_info;
  basisfn_  = &basisfn;
  targetfn_ = &targetfn;
  params_ = params;
  cell_ = hkl_info.cell();

  int nparams = basisfn_->num_params();
  Matrix<> dfdp2( nparams, nparams );
  Matrix<> drdp2( nparams, nparams );
  std::vector<ftype> drdp( nparams ), dfdp( nparams );
  std::vector<ftype> shiftn( nparams );
  params_.resize( nparams );

  // loop for 20 cycles refining the params
  g0 = 1.0e25;
  for ( n = 0; n < 20; n++ ) {
    // calc target fn and derivs
    calc_derivs( params_, r, drdp, drdp2 );
    // add parameter restraint terms
    for ( k = 0; k < nparams; k++ ) {
      ftype d = params[k] - targetval[k];
      r          +=       targetwgt[k] * d * d;
      drdp[k]    += 2.0 * targetwgt[k] * d;
      drdp2(k,k) += 2.0 * targetwgt[k];
    }

    g1 = 0.0;
    for ( k = 0; k < nparams; k++ ) g1 += drdp[k]*drdp[k];
    g1 = sqrt(g1);

    // stop if gradient starts increasing
    if ( g1 < 1.0e-10 || g1 >= g0 ) break;
    g0 = g1;

    // solve for Newton-Raphson shift
    shiftn = drdp2.solve( drdp );

    // apply shift
    scale = (1.0+ftype(n)) / (1.0+ftype(n)+damp);
    for ( i = 0; i< nparams; i++ ) params_[i] -= scale*shiftn[i];

    if ( debug ) std::cout << " Resolution function cycle " << n << " " << g0 << " " << g1 << " " << scale << "\n";

    // stop if the target is quadratic
    if ( basisfn.type()  == BasisFn_base::LINEAR &&
	 targetfn.type() == TargetFn_base::QUADRATIC ) break;
  }
}
*/


} // namespace clipper
