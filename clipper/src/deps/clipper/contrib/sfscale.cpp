/* sfscale.cpp: structure factor anisotropic scaling implementation */
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


#include "sfscale.h"
#include "../core/resol_targetfn.h"
#include "../core/hkl_compute.h"


namespace clipper {

template<class T> bool SFscale_aniso<T>::operator() ( HKL_data<datatypes::F_sigF<T> >& fo, const HKL_data<datatypes::F_phi<T> >& fc )
{
  typedef HKL_info::HKL_reference_index HRI;
  // expand to P1 in order to preserve symmetry
  const HKL_info& hkls = fo.hkl_info();
  Spacegroup spgrp1( Spacegroup::P1 );
  HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
  HKL_data<datatypes::F_sigF<T> > fo1( hkl1 );
  HKL_data<datatypes::F_phi<T> >  fc1( hkl1 );
  for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
    datatypes::F_sigF<T> f = fo[ih.hkl()];
    if ( f.f() >= nsig_ * f.sigf() ) {
      fo1[ih] = f;
      fc1[ih] = fc[ih.hkl()];
    }
  }
  // do the aniso scaling
  std::vector<double> param( 7, 0.0 );
  BasisFn_log_aniso_gaussian bfn;
  TargetFn_scaleLogF1F2<datatypes::F_sigF<T>,datatypes::F_phi<T> >
    tfn( fo1, fc1 );
  ResolutionFn rfn( hkl1, bfn, tfn, param );
  for ( HRI ih = hkls.first(); !ih.last(); ih.next() )
    if ( !fo[ih].missing() )
      fo[ih].scale( exp( 0.5*bfn.f(ih.hkl(),hkls.cell(),rfn.params()) ) );
  u_i = bfn.u_aniso_orth( rfn.params() );
  u_f = 0.5 * u_i;
  return true;
}


template<class T> bool SFscale_aniso<T>::operator() ( HKL_data<datatypes::F_phi<T> >& fc, const HKL_data<datatypes::F_sigF<T> >& fo )
{
  typedef HKL_info::HKL_reference_index HRI;
  // expand to P1 in order to preserve symmetry
  const HKL_info& hkls = fo.hkl_info();
  Spacegroup spgrp1( Spacegroup::P1 );
  HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
  HKL_data<datatypes::F_sigF<T> > fo1( hkl1 );
  HKL_data<datatypes::F_phi<T> >  fc1( hkl1 );
  for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
    datatypes::F_sigF<T> f = fo[ih.hkl()];
    if ( f.f() >= nsig_ * f.sigf() ) {
      fo1[ih] = f;
      fc1[ih] = fc[ih.hkl()];
    }
  }
  // do the aniso scaling
  std::vector<double> param( 7, 0.0 );
  BasisFn_log_aniso_gaussian bfn;
  TargetFn_scaleLogF1F2<datatypes::F_phi<T>,datatypes::F_sigF<T> >
    tfn( fc1, fo1 );
  ResolutionFn rfn( hkl1, bfn, tfn, param );
  for ( HRI ih = hkls.first(); !ih.last(); ih.next() )
    if ( !fc[ih].missing() )
      fc[ih].scale( exp( 0.5*bfn.f(ih.hkl(),hkls.cell(),rfn.params()) ) );
  u_i = bfn.u_aniso_orth( rfn.params() );
  u_f = 0.5 * u_i;
  return true;
}


template<class T> bool SFscale_aniso<T>::operator() ( HKL_data<datatypes::F_sigF<T> >& fo )
{
  typedef datatypes::F_sigF<T> DATA;
  typedef TargetFn_scaleF1F2<DATA,DATA> TGT1;
  typedef TargetFn_scaleLogF1F2<DATA,DATA> TGT2;
  typedef BasisFn_spline SCALETYPE;
  return scale<DATA,TGT1,TGT2,SCALETYPE>( fo, -1.0, 12 );
}


template<class T> bool SFscale_aniso<T>::operator() ( HKL_data<datatypes::F_sigF<T> >& fo, ftype resfilter, const int npar_scl )
{
  typedef datatypes::F_sigF<T> DATA;
  typedef TargetFn_scaleF1F2<DATA,DATA> TGT1;
  typedef TargetFn_scaleLogF1F2<DATA,DATA> TGT2;
  typedef BasisFn_spline SCALETYPE;
  return scale<DATA,TGT1,TGT2,SCALETYPE>( fo, resfilter, npar_scl );
}


template<class T> bool SFscale_aniso<T>::operator() ( HKL_data<datatypes::I_sigI<T> >& Io, ftype resfilter, const int npar_scl )
{
  typedef datatypes::I_sigI<T> DATA;
  typedef TargetFn_scaleI1I2<DATA,DATA> TGT1;
  typedef TargetFn_scaleLogI1I2<DATA,DATA> TGT2;
  typedef BasisFn_spline SCALETYPE;
  return scale<DATA,TGT1,TGT2,SCALETYPE>( Io, resfilter, npar_scl );
}


template<class T> template<class D, class T1, class T2, class S> bool SFscale_aniso<T>::scale( HKL_data<D>& fo, const ftype resfilter, const int npar_scl )
{
  typedef HKL_info::HKL_reference_index HRI;
  // expand to P1 in order to preserve symmetry
  const HKL_info& hkls = fo.hkl_info();
  Spacegroup spgrp1( Spacegroup::P1 );
  HKL_info hkl1( spgrp1, hkls.cell(), hkls.resolution(), true );
  HKL_data<D> fo1( hkl1 ), fs1( hkl1 ), fc1( hkl1 );
  for ( HRI ih = hkl1.first(); !ih.last(); ih.next() ) {
    D f = fo[ih.hkl()];
    if ( obs(f) >= nsig_ * sigobs(f) ) fo1[ih] = f;
  }

  // calc target values
  fc1 = D( 1.0, 1.0 );
  for ( HRI ih = fc1.first(); !ih.last(); ih.next() )
    fc1[ih].scale( sqrt( ih.hkl_class().epsilon() ) );

  // do correction
  u_i = u_f = U_aniso_orth( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
  std::vector<ftype> param( 7, 0.0 );
  BasisFn_log_aniso_gaussian bfn;
  ftype dp;
  for ( int c = 0; c < 2; c++ ) {
    // remove current anistropy estimate
    datatypes::Compute_scale_u_aniso<D> compute_s(1.0,-u_f);
    fs1.compute( fo1, compute_s );

    /*{
    T2 tfn( fs1, fc1 );
    param = std::vector<ftype>( 7, 0.0 );
    ResolutionFn rfn( hkl1, bfn, tfn, param );
    param = rfn.params();
    dp = (param[1]+param[2]+param[3])/3.0;
    param[1] -= dp; param[2] -= dp; param[3] -= dp;
    std::cout << "#" << param[1] << " " << param[2] << " " << param[3] << " " << param[4] << " " << param[5] << " " << param[6] << "\n"; std::cout << "\n";
    }*/

    // and calc E scale
    S basissc( fs1, npar_scl, 1.0 );
    std::vector<ftype> param_fo( basissc.num_params(), 1.0 );
    T1 tfn_fo( fs1, fc1 );
    ResolutionFn rfn_fo( hkl1, basissc, tfn_fo, param_fo );
    param_fo = rfn_fo.params();

    // prescale F to E-like scale
    for ( HRI ih = fs1.first(); !ih.last(); ih.next() ) {
      fs1[ih] = fo1[ih];
      fs1[ih].scale( sqrt( basissc.f_s( ih.invresolsq(), param_fo ) ) );
    }

    /*for ( int i = 0; i < param_fo.size(); i++ ) std::cout << param_fo[i] << " "; std::cout << std::endl;
    {
    S basissc( fs1, npar_scl, 1.0 );
    std::vector<ftype> param_fo( basissc.num_params(), 1.0 );
    T1 tfn_fo( fs1, fc1 );
    ResolutionFn rfn_fo( hkl1, basissc, tfn_fo, param_fo );
    param_fo = rfn_fo.params();
    for ( int i = 0; i < param_fo.size(); i++ ) std::cout << param_fo[i] << " "; std::cout << std::endl;
    }*/

    // if required, weight low res terms towards isotropic
    if ( resfilter > 0.0 ) {
      ftype ssq = clipper::Resolution( resfilter ).invresolsq_limit();
      for ( HRI ih = fs1.first(); !ih.last(); ih.next() )
	if ( !fs1[ih].missing() ) {
	  ftype w = exp( -0.5 * ih.invresolsq() / ssq );
	  fs1[ih] = D( exp( (1.0-w)*log(obs(fs1[ih]))+(w)*log(obs(fc1[ih])) ),
		       sigobs(fs1[ih]) );
	}
    }

    // do the aniso scaling
    T2 tfn( fs1, fc1 );
    param = std::vector<ftype>( 7, 0.0 );
    ResolutionFn rfn( hkl1, bfn, tfn, param );
    param = rfn.params();
    // set trace to zero (i.e. no isotropic correction)
    dp = (param[1]+param[2]+param[3])/3.0;
    param[1] -= dp; param[2] -= dp; param[3] -= dp;
    u_i = bfn.u_aniso_orth( param );
    u_f = 0.5 * u_i;
    //std::cout << c << " | " << param[1] << " " << param[2] << " " << param[3] << " " << param[4] << " " << param[5] << " " << param[6] << "\n"; std::cout << " DP " << dp << "\n";
  }

  // sharpen or smooth as required
  Matrix<ftype> m(3,3);
  m(0,0)=       param[1]; m(1,1)=       param[2]; m(2,2)=       param[3];
  m(0,1)=m(1,0)=param[4]; m(0,2)=m(2,0)=param[5]; m(1,2)=m(2,1)=param[6];
  std::vector<ftype> ev = m.eigen();
  //std::cout << "EIGEN " << param[1] << " " << param[2] << " " << param[3] << " " << ev[0] << " " << ev[1] << " " << ev[2] << std::endl;
  dp = 0.0;
  if ( mode_ == SHARPEN   ) dp = ev[2];
  if ( mode_ == UNSHARPEN ) dp = ev[0];
  param[1] -= dp; param[2] -= dp; param[3] -= dp;
  u_i = bfn.u_aniso_orth( param );
  u_f = 0.5 * u_i;

  // store the results
  datatypes::Compute_scale_u_aniso<D> compute_s(1.0,-u_f);
  fo.compute( fo, compute_s );

  return true;
}


template<class T> const U_aniso_orth& SFscale_aniso<T>::u_aniso_orth( TYPE t ) const
{
  if ( t == I ) return u_i;
  else          return u_f;
}


// compile templates

template class CLIPPER_IMEX SFscale_aniso<ftype32>;
template CLIPPER_IMEX bool SFscale_aniso<ftype32>::scale<datatypes::F_sigF<ftype32>,TargetFn_scaleF1F2<datatypes::F_sigF<ftype32>,datatypes::F_sigF<ftype32> >,TargetFn_scaleLogF1F2<datatypes::F_sigF<ftype32>,datatypes::F_sigF<ftype32> >,BasisFn_binner>( HKL_data<datatypes::F_sigF<ftype32> >& fo, const ftype resfilter, const int npar_scl );
template CLIPPER_IMEX bool SFscale_aniso<ftype32>::scale<datatypes::F_sigF<ftype32>,TargetFn_scaleF1F2<datatypes::F_sigF<ftype32>,datatypes::F_sigF<ftype32> >,TargetFn_scaleLogF1F2<datatypes::F_sigF<ftype32>,datatypes::F_sigF<ftype32> >,BasisFn_linear>( HKL_data<datatypes::F_sigF<ftype32> >& fo, const ftype resfilter, const int npar_scl );
template CLIPPER_IMEX bool SFscale_aniso<ftype32>::scale<datatypes::F_sigF<ftype32>,TargetFn_scaleF1F2<datatypes::F_sigF<ftype32>,datatypes::F_sigF<ftype32> >,TargetFn_scaleLogF1F2<datatypes::F_sigF<ftype32>,datatypes::F_sigF<ftype32> >,BasisFn_spline>( HKL_data<datatypes::F_sigF<ftype32> >& fo, const ftype resfilter, const int npar_scl );
template CLIPPER_IMEX bool SFscale_aniso<ftype32>::scale<datatypes::I_sigI<ftype32>,TargetFn_scaleI1I2<datatypes::I_sigI<ftype32>,datatypes::I_sigI<ftype32> >,TargetFn_scaleLogI1I2<datatypes::I_sigI<ftype32>,datatypes::I_sigI<ftype32> >,BasisFn_binner>( HKL_data<datatypes::I_sigI<ftype32> >& fo, const ftype resfilter, const int npar_scl );
template CLIPPER_IMEX bool SFscale_aniso<ftype32>::scale<datatypes::I_sigI<ftype32>,TargetFn_scaleI1I2<datatypes::I_sigI<ftype32>,datatypes::I_sigI<ftype32> >,TargetFn_scaleLogI1I2<datatypes::I_sigI<ftype32>,datatypes::I_sigI<ftype32> >,BasisFn_linear>( HKL_data<datatypes::I_sigI<ftype32> >& fo, const ftype resfilter, const int npar_scl );
template CLIPPER_IMEX bool SFscale_aniso<ftype32>::scale<datatypes::I_sigI<ftype32>,TargetFn_scaleI1I2<datatypes::I_sigI<ftype32>,datatypes::I_sigI<ftype32> >,TargetFn_scaleLogI1I2<datatypes::I_sigI<ftype32>,datatypes::I_sigI<ftype32> >,BasisFn_spline>( HKL_data<datatypes::I_sigI<ftype32> >& fo, const ftype resfilter, const int npar_scl );

template class CLIPPER_IMEX SFscale_aniso<ftype64>;
template CLIPPER_IMEX bool SFscale_aniso<ftype64>::scale<datatypes::F_sigF<ftype64>,TargetFn_scaleF1F2<datatypes::F_sigF<ftype64>,datatypes::F_sigF<ftype64> >,TargetFn_scaleLogF1F2<datatypes::F_sigF<ftype64>,datatypes::F_sigF<ftype64> >,BasisFn_binner>( HKL_data<datatypes::F_sigF<ftype64> >& fo, const ftype resfilter, const int npar_scl );
template CLIPPER_IMEX bool SFscale_aniso<ftype64>::scale<datatypes::F_sigF<ftype64>,TargetFn_scaleF1F2<datatypes::F_sigF<ftype64>,datatypes::F_sigF<ftype64> >,TargetFn_scaleLogF1F2<datatypes::F_sigF<ftype64>,datatypes::F_sigF<ftype64> >,BasisFn_linear>( HKL_data<datatypes::F_sigF<ftype64> >& fo, const ftype resfilter, const int npar_scl );
template CLIPPER_IMEX bool SFscale_aniso<ftype64>::scale<datatypes::F_sigF<ftype64>,TargetFn_scaleF1F2<datatypes::F_sigF<ftype64>,datatypes::F_sigF<ftype64> >,TargetFn_scaleLogF1F2<datatypes::F_sigF<ftype64>,datatypes::F_sigF<ftype64> >,BasisFn_spline>( HKL_data<datatypes::F_sigF<ftype64> >& fo, const ftype resfilter, const int npar_scl );
template CLIPPER_IMEX bool SFscale_aniso<ftype64>::scale<datatypes::I_sigI<ftype64>,TargetFn_scaleI1I2<datatypes::I_sigI<ftype64>,datatypes::I_sigI<ftype64> >,TargetFn_scaleLogI1I2<datatypes::I_sigI<ftype64>,datatypes::I_sigI<ftype64> >,BasisFn_binner>( HKL_data<datatypes::I_sigI<ftype64> >& fo, const ftype resfilter, const int npar_scl );
template CLIPPER_IMEX bool SFscale_aniso<ftype64>::scale<datatypes::I_sigI<ftype64>,TargetFn_scaleI1I2<datatypes::I_sigI<ftype64>,datatypes::I_sigI<ftype64> >,TargetFn_scaleLogI1I2<datatypes::I_sigI<ftype64>,datatypes::I_sigI<ftype64> >,BasisFn_linear>( HKL_data<datatypes::I_sigI<ftype64> >& fo, const ftype resfilter, const int npar_scl );
template CLIPPER_IMEX bool SFscale_aniso<ftype64>::scale<datatypes::I_sigI<ftype64>,TargetFn_scaleI1I2<datatypes::I_sigI<ftype64>,datatypes::I_sigI<ftype64> >,TargetFn_scaleLogI1I2<datatypes::I_sigI<ftype64>,datatypes::I_sigI<ftype64> >,BasisFn_spline>( HKL_data<datatypes::I_sigI<ftype64> >& fo, const ftype resfilter, const int npar_scl );

} // namespace clipper
