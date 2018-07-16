/* hkl_compute.cpp: fundamental conversion ops for the clipper libraries */
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


#include "hkl_compute.h"


namespace clipper {

namespace datatypes {


template<class dtype> Compute_phifom_from_abcd<dtype>::Compute_phifom_from_abcd()
{
  for ( int i = 0; i < 144; i++ ) {
    costab[i] = cos( Util::d2rad(ftype(5*i)) );
    sintab[i] = sin( Util::d2rad(ftype(5*i)) );
  }
}

template<class dtype> const Phi_fom<dtype> Compute_phifom_from_abcd<dtype>::operator()( const HKL_info::HKL_reference_index& ih, const ABCD<dtype>& abcd ) const
{
  clipper::datatypes::Phi_fom<dtype> phifom;
  if ( !abcd.missing() ) {
    const HKL_class& cls = ih.hkl_class();
    Range<ftype64> qrange( -700.0, 700.0 );
    ftype64 q, pq, scos, ssin, ssum;
    if ( cls.centric() ) {
      // centric case: integrate over 2 allowed phases
      int i = Util::mod( Util::intr( Util::rad2d( cls.allowed() ) / 5 ), 36 );
      q = abcd.a()*costab[i] + abcd.b()*sintab[i];
      pq = exp( qrange.truncate(q) );
      scos = ( pq - 1.0/pq ) * costab[i];
      ssin = ( pq - 1.0/pq ) * sintab[i];
      ssum = ( pq + 1.0/pq );
    } else {
      // acentric case: integrate over phase circle in 5 degree steps
      scos = ssin = ssum = 0.0;
      for ( int i = 0; i < 72; i++ ) {
	q = abcd.a()*costab[i]   + abcd.b()*sintab[i]
	  + abcd.c()*costab[2*i] + abcd.d()*sintab[2*i];
	pq = exp( qrange.truncate(q) );
	scos += pq * costab[i];
	ssin += pq * sintab[i];
	ssum += pq;
      }
    }
    std::complex<ftype64> c( scos/ssum, ssin/ssum );
    phifom.phi() = std::arg(c);
    phifom.fom() = std::abs(c);
  }
  return phifom;
}


template<class dtype> const ABCD<dtype> Compute_abcd_from_phifom<dtype>::operator()( const HKL_info::HKL_reference_index& ih, const Phi_fom<dtype>& phifom ) const
{
  clipper::datatypes::ABCD<dtype> abcd;
  if ( !phifom.missing() ) {
    const HKL_class& cls = ih.hkl_class();
    dtype x = Util::min( phifom.fom(), dtype(0.9999) );
    if ( cls.centric() ) x = Util::atanh( x );
    else                 x = Util::invsim( x );
    abcd.a() = x * cos( phifom.phi() );
    abcd.b() = x * sin( phifom.phi() );
    abcd.c() = abcd.d() = 0.0;
  }
  return abcd;
}


template<class dtype> const F_phi<dtype> Compute_fphi_from_fsigf_phifom<dtype>::operator()( const HKL_info::HKL_reference_index& , const F_sigF<dtype>& fsigf, const Phi_fom<dtype>& phifom ) const
{
  clipper::datatypes::F_phi<dtype> fphi;
  if ( !fsigf.missing() && !phifom.missing() ) {
    fphi.f() = fsigf.f() * phifom.fom();
    fphi.phi() = phifom.phi();
  }
  return fphi;
}


template<class dtype> const E_sigE<dtype> Compute_EsigE_from_FsigF<dtype>::operator()( const HKL_info::HKL_reference_index& ih, const F_sigF<dtype>& fsigf ) const
{
  clipper::datatypes::E_sigE<dtype> esige;
  if ( !fsigf.missing() ) {
    ftype sqrteps = sqrt( ih.hkl_class().epsilon() );
    esige.E() = fsigf.f()/sqrteps;
    esige.sigE() = fsigf.sigf()/sqrteps;
  }
  return esige;
}

template<class dtype> const F_sigF<dtype> Compute_mean_fsigf_from_fsigfano<dtype>::operator()( const HKL_info::HKL_reference_index& , const F_sigF_ano<dtype>& fsigfano ) const
{
  clipper::datatypes::F_sigF<dtype> fsigf;
  if        ( Util::is_nan( fsigfano.f_pl() ) ) {
    fsigf.f() = fsigfano.f_mi();
    fsigf.sigf() = fsigfano.sigf_mi();
  } else if ( Util::is_nan( fsigfano.f_mi() ) ) {
    fsigf.f() = fsigfano.f_pl();
    fsigf.sigf() = fsigfano.sigf_pl();
  } else {
    fsigf.f() = 0.5 * ( fsigfano.f_pl() + fsigfano.f_mi() );
    if ( Util::is_nan(fsigfano.cov()) )
      fsigf.sigf() = 0.5 * sqrt( fsigfano.sigf_pl()*fsigfano.sigf_pl() +
				 fsigfano.sigf_mi()*fsigfano.sigf_mi() );
    else
      fsigf.sigf() = 0.5 * sqrt( fsigfano.sigf_pl()*fsigfano.sigf_pl() +
				 fsigfano.sigf_mi()*fsigfano.sigf_mi() +
				 2.0 * fsigfano.cov() );
  }
  return fsigf;
}

template<class dtype> const F_sigF<dtype> Compute_diff_fsigf_from_fsigfano<dtype>::operator()( const HKL_info::HKL_reference_index& , const F_sigF_ano<dtype>& fsigfano ) const
{
  clipper::datatypes::F_sigF<dtype> fsigf;
  if ( Util::is_nan( fsigfano.f_pl() ) || Util::is_nan( fsigfano.f_mi() ) ) {
    fsigf.f() = fsigf.sigf() = Util::nan();
  } else {
    fsigf.f() = 1.0 * ( fsigfano.f_pl() - fsigfano.f_mi() );
    if ( Util::is_nan(fsigfano.cov()) )
      fsigf.sigf() = 1.0 * sqrt( fsigfano.sigf_pl()*fsigfano.sigf_pl() +
				 fsigfano.sigf_mi()*fsigfano.sigf_mi() );
    else
      fsigf.sigf() = 1.0 * sqrt( fsigfano.sigf_pl()*fsigfano.sigf_pl() +
				 fsigfano.sigf_mi()*fsigfano.sigf_mi() -
				 2.0 * fsigfano.cov() );
  }
  return fsigf;
}

template<class dtype> const F_phi<dtype> Compute_neg_fphi<dtype>::operator()( const HKL_info::HKL_reference_index& , const F_phi<dtype>& fphi1 ) const
{
  clipper::datatypes::F_phi<dtype> fphi;
  if ( !fphi.missing() )
    fphi = clipper::datatypes::F_phi<dtype>( -std::complex<dtype>(fphi1) );
  return fphi;
}


template<class dtype> const F_phi<dtype> Compute_add_fphi<dtype>::operator()( const HKL_info::HKL_reference_index& , const F_phi<dtype>& fphi1, const F_phi<dtype>& fphi2 ) const
{
  clipper::datatypes::F_phi<dtype> fphi;
  if ( !fphi1.missing() && !fphi2.missing() )
    fphi = clipper::datatypes::F_phi<dtype>( std::complex<dtype>(fphi1) + std::complex<dtype>(fphi2) );
  return fphi;
}


template<class dtype> const F_phi<dtype> Compute_sub_fphi<dtype>::operator()( const HKL_info::HKL_reference_index& , const F_phi<dtype>& fphi1, const F_phi<dtype>& fphi2 ) const
{
  clipper::datatypes::F_phi<dtype> fphi;
  if ( !fphi1.missing() && !fphi2.missing() )
    fphi = clipper::datatypes::F_phi<dtype>( std::complex<dtype>(fphi1) - std::complex<dtype>(fphi2) );
  return fphi;
}


template<class dtype> const ABCD<dtype> Compute_add_abcd<dtype>::operator()( const HKL_info::HKL_reference_index& , const ABCD<dtype>& abcd1, const ABCD<dtype>& abcd2 ) const
{
  clipper::datatypes::ABCD<dtype> abcd;
  if ( !abcd1.missing() && !abcd2.missing() ) {
    abcd.a() = abcd1.a() + abcd2.a();
    abcd.b() = abcd1.b() + abcd2.b();
    abcd.c() = abcd1.c() + abcd2.c();
    abcd.d() = abcd1.d() + abcd2.d();
  }
  return abcd;
}


/*! DEPRECATED: This operator applies the scale factor against
  intensities and the U value against magnitudes, which is
  counterintuitive. Compute_scale_u_iso is more intutive.

  Construct conversion operator to scale list according to the formula
  I_new = s.exp( b.|h|^2/2 ) I_old
  or
  F_new^2 = s.exp( b.|h|^2/2 ) F_old^2
  where |h| = invresolsq.
  \param s The intensity scale factor.
  \param u The temperature factor (U-value). */
template<class T> Compute_scale_u<T>::Compute_scale_u( const ftype& s, const ftype& u )
{
  s_ = sqrt(s); u_ = Util::twopi2()*u;
}

template<class T> const T Compute_scale_u<T>::operator()( const HKL_info::HKL_reference_index& ih, T data ) const
{
  if ( !data.missing() )
    data.scale( s_ * exp( u_ * ih.invresolsq() ) );
  return data;
}


/*! Construct conversion operator to scale list according to the formula
  I_new = s^2 .exp( 4 &pi;^2 u.|h|^2 ) I_old
  or
  F_new = s.exp( 2 &pi;^2 u.|h|^2 ) F_old
  where |h|^2 = invresolsq.
  \param s The scale factor.
  \param u The temperature factor (U-value). */
template<class T> Compute_scale_u_iso<T>::Compute_scale_u_iso( const ftype& s, const ftype& u )
{
  s_ = s; u_ = Util::twopi2()*u;
}

template<class T> const T Compute_scale_u_iso<T>::operator()( const HKL_info::HKL_reference_index& ih, T data ) const
{
  if ( !data.missing() )
    data.scale( s_ * exp( u_ * ih.invresolsq() ) );
  return data;
}


/*! Construct conversion operator to scale list according to the formula
  I_new = s^2 .exp( 4 &pi;^2 h^T U h ) I_old
  or
  F_new = s.exp( 2 &pi;^2 h^T U h ) F_old
  \param s The scale factor.
  \param u The temperature factor (U-value). */
template<class T> Compute_scale_u_aniso<T>::Compute_scale_u_aniso( const ftype& s, const U_aniso_orth& u )
{
  s_ = s; u_ = Util::twopi2()*u;
}

template<class T> const T Compute_scale_u_aniso<T>::operator()( const HKL_info::HKL_reference_index& ih, T data ) const
{
  if ( !data.missing() )
    data.scale( s_ * exp( u_.quad_form( ih.hkl().coord_reci_orth( ih.base_hkl_info().cell() ) ) ) );
  return data;
}


template<class dtype, class T> const F_sigF<dtype> Compute_FsigF<dtype,T>::operator()( const HKL_info::HKL_reference_index& , const T& fsigf ) const
{
  clipper::datatypes::F_sigF<dtype> fsigfnew;
  if ( !fsigf.missing() ) {
    fsigfnew.f() = fsigf.f();
    fsigfnew.sigf() = fsigf.sigf();
  }
  return fsigfnew;
}


// compile template types

template class CLIPPER_IMEX Compute_phifom_from_abcd<ftype32>;
template class CLIPPER_IMEX Compute_abcd_from_phifom<ftype32>;
template class CLIPPER_IMEX Compute_fphi_from_fsigf_phifom<ftype32>;
template class CLIPPER_IMEX Compute_EsigE_from_FsigF<ftype32>;
template class CLIPPER_IMEX Compute_mean_fsigf_from_fsigfano<ftype32>;
template class CLIPPER_IMEX Compute_diff_fsigf_from_fsigfano<ftype32>;
template class CLIPPER_IMEX Compute_neg_fphi<ftype32>;
template class CLIPPER_IMEX Compute_add_fphi<ftype32>;
template class CLIPPER_IMEX Compute_sub_fphi<ftype32>;
template class CLIPPER_IMEX Compute_add_abcd<ftype32>;
template class CLIPPER_IMEX Compute_scale_u_iso<datatypes::I_sigI<ftype32> >;
template class CLIPPER_IMEX Compute_scale_u_iso<datatypes::I_sigI_ano<ftype32> >;
template class CLIPPER_IMEX Compute_scale_u_iso<datatypes::F_sigF<ftype32> >;
template class CLIPPER_IMEX Compute_scale_u_iso<datatypes::F_sigF_ano<ftype32> >;
template class CLIPPER_IMEX Compute_scale_u_iso<datatypes::F_phi<ftype32> >;
template class CLIPPER_IMEX Compute_scale_u_aniso<datatypes::I_sigI<ftype32> >;
template class CLIPPER_IMEX Compute_scale_u_aniso<datatypes::I_sigI_ano<ftype32> >;
template class CLIPPER_IMEX Compute_scale_u_aniso<datatypes::F_sigF<ftype32> >;
template class CLIPPER_IMEX Compute_scale_u_aniso<datatypes::F_sigF_ano<ftype32> >;
template class CLIPPER_IMEX Compute_scale_u_aniso<datatypes::F_phi<ftype32> >;
template class CLIPPER_IMEX Compute_FsigF<ftype32, F_sigF<ftype32> >;
template class CLIPPER_IMEX Compute_FsigF<ftype32, F_sigF_ano<ftype32> >;

template class CLIPPER_IMEX Compute_phifom_from_abcd<ftype64>;
template class CLIPPER_IMEX Compute_abcd_from_phifom<ftype64>;
template class CLIPPER_IMEX Compute_fphi_from_fsigf_phifom<ftype64>;
template class CLIPPER_IMEX Compute_EsigE_from_FsigF<ftype64>;
template class CLIPPER_IMEX Compute_mean_fsigf_from_fsigfano<ftype64>;
template class CLIPPER_IMEX Compute_diff_fsigf_from_fsigfano<ftype64>;
template class CLIPPER_IMEX Compute_neg_fphi<ftype64>;
template class CLIPPER_IMEX Compute_add_fphi<ftype64>;
template class CLIPPER_IMEX Compute_sub_fphi<ftype64>;
template class CLIPPER_IMEX Compute_add_abcd<ftype64>;
template class CLIPPER_IMEX Compute_scale_u_iso<datatypes::I_sigI<ftype64> >;
template class CLIPPER_IMEX Compute_scale_u_iso<datatypes::I_sigI_ano<ftype64> >;
template class CLIPPER_IMEX Compute_scale_u_iso<datatypes::F_sigF<ftype64> >;
template class CLIPPER_IMEX Compute_scale_u_iso<datatypes::F_sigF_ano<ftype64> >;
template class CLIPPER_IMEX Compute_scale_u_iso<datatypes::F_phi<ftype64> >;
template class CLIPPER_IMEX Compute_scale_u_aniso<datatypes::I_sigI<ftype64> >;
template class CLIPPER_IMEX Compute_scale_u_aniso<datatypes::I_sigI_ano<ftype64> >;
template class CLIPPER_IMEX Compute_scale_u_aniso<datatypes::F_sigF<ftype64> >;
template class CLIPPER_IMEX Compute_scale_u_aniso<datatypes::F_sigF_ano<ftype64> >;
template class CLIPPER_IMEX Compute_scale_u_aniso<datatypes::F_phi<ftype64> >;
template class CLIPPER_IMEX Compute_FsigF<ftype64, F_sigF<ftype64> >;
template class CLIPPER_IMEX Compute_FsigF<ftype64, F_sigF_ano<ftype64> >;


// DEPRECATED
template class Compute_scale_u<datatypes::I_sigI<ftype32> >;
template class Compute_scale_u<datatypes::I_sigI_ano<ftype32> >;
template class Compute_scale_u<datatypes::F_sigF<ftype32> >;
template class Compute_scale_u<datatypes::F_sigF_ano<ftype32> >;
template class Compute_scale_u<datatypes::F_phi<ftype32> >;
template class Compute_scale_u<datatypes::I_sigI<ftype64> >;
template class Compute_scale_u<datatypes::I_sigI_ano<ftype64> >;
template class Compute_scale_u<datatypes::F_sigF<ftype64> >;
template class Compute_scale_u<datatypes::F_sigF_ano<ftype64> >;
template class Compute_scale_u<datatypes::F_phi<ftype64> >;

} // namespace datatypes

} // namespace clipper
