/* hkl_operators.cpp: HKL_data operators for the clipper libraries */
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


#include "hkl_operators.h"


namespace clipper {

namespace datatypes {

HKL_data<Flag_bool> operator &( const HKL_data_base& d1, const HKL_data_base& d2 )
{
  HKL_data<Flag_bool> result( d1.base_hkl_info(), d1.base_cell() );
  for ( int i = 0; i < d1.base_hkl_info().num_reflections(); i++ )
    result[i].flag() = ( !d1.missing(i) ) && ( !d2.missing(i) );
  return result;
}

HKL_data<Flag_bool> operator |( const HKL_data_base& d1, const HKL_data_base& d2 )
{
  HKL_data<Flag_bool> result( d1.base_hkl_info(), d1.base_cell() );
  for ( int i = 0; i < d1.base_hkl_info().num_reflections(); i++ )
    result[i].flag() = ( !d1.missing(i) ) || ( !d2.missing(i) );
  return result;
}

HKL_data<Flag_bool> operator ^( const HKL_data_base& d1, const HKL_data_base& d2 )
{
  HKL_data<Flag_bool> result( d1.base_hkl_info(), d1.base_cell() );
  for ( int i = 0; i < d1.base_hkl_info().num_reflections(); i++ )
    result[i].flag() = ( !d1.missing(i) ) ^ ( !d2.missing(i) );
  return result;
}

HKL_data<Flag_bool> operator !( const HKL_data_base& d1 )
{
  HKL_data<Flag_bool> result( d1.base_hkl_info(), d1.base_cell() );
  for ( int i = 0; i < d1.base_hkl_info().num_reflections(); i++ )
    result[i].flag() = d1.missing(i);
  return result;
}

HKL_data<Flag_bool> operator ==( const HKL_data<Flag>& d1, const int& n )
{
  HKL_data<Flag_bool> result( d1.base_hkl_info(), d1.base_cell() );
  for ( int i = 0; i < d1.base_hkl_info().num_reflections(); i++ )
    result[i].flag() = ( d1[i].flag() == n );
  return result;
}

HKL_data<Flag_bool> operator !=( const HKL_data<Flag>& d1, const int& n )
{
  HKL_data<Flag_bool> result( d1.base_hkl_info(), d1.base_cell() );
  for ( int i = 0; i < d1.base_hkl_info().num_reflections(); i++ )
    result[i].flag() = ( d1[i].flag() != n );
  return result;
}

HKL_data<Flag_bool> operator >=( const HKL_data<Flag>& d1, const int& n )
{
  HKL_data<Flag_bool> result( d1.base_hkl_info(), d1.base_cell() );
  for ( int i = 0; i < d1.base_hkl_info().num_reflections(); i++ )
    result[i].flag() = ( d1[i].flag() >= n );
  return result;
}

HKL_data<Flag_bool> operator <=( const HKL_data<Flag>& d1, const int& n )
{
  HKL_data<Flag_bool> result( d1.base_hkl_info(), d1.base_cell() );
  for ( int i = 0; i < d1.base_hkl_info().num_reflections(); i++ )
    result[i].flag() = ( d1[i].flag() <= n );
  return result;
}

HKL_data<Flag_bool> operator >( const HKL_data<Flag>& d1, const int& n )
{
  HKL_data<Flag_bool> result( d1.base_hkl_info(), d1.base_cell() );
  for ( int i = 0; i < d1.base_hkl_info().num_reflections(); i++ )
    result[i].flag() = ( d1[i].flag() > n );
  return result;
}

HKL_data<Flag_bool> operator <( const HKL_data<Flag>& d1, const int& n )
{
  HKL_data<Flag_bool> result( d1.base_hkl_info(), d1.base_cell() );
  for ( int i = 0; i < d1.base_hkl_info().num_reflections(); i++ )
    result[i].flag() = ( d1[i].flag() < n );
  return result;
}

template<class dtype> F_phi<dtype> operator +( const F_phi<dtype>& d1, const F_phi<dtype>& d2 )
{
  return clipper::datatypes::F_phi<dtype>( std::complex<dtype>(d1) + std::complex<dtype>(d2) );
}

template<class dtype> F_phi<dtype> operator -( const F_phi<dtype>& d1, const F_phi<dtype>& d2 )
{
  return clipper::datatypes::F_phi<dtype>( std::complex<dtype>(d1) - std::complex<dtype>(d2) );
}

template<class dtype> F_phi<dtype> operator -( const F_phi<dtype>& d1 )
{
  return clipper::datatypes::F_phi<dtype>( -std::complex<dtype>(d1) );
}

template<class dtype> ABCD<dtype> operator +( const ABCD<dtype>& d1, const ABCD<dtype>& d2 )
{
  ABCD<dtype> result;
  result.a() = d1.a() + d2.a();
  result.b() = d1.b() + d2.b();
  result.c() = d1.c() + d2.c();
  result.d() = d1.d() + d2.d();
  return result;
}

template<class dtype> HKL_data<F_phi<dtype> > operator +( const HKL_data<F_phi<dtype> >& d1, const HKL_data<F_phi<dtype> >& d2 )
{
  HKL_data<F_phi<dtype> > result( d1.base_hkl_info(), d1.base_cell() );
  for ( int i = 0; i < d1.base_hkl_info().num_reflections(); i++ )
    if ( !d1.missing(i) && !d2.missing(i) )
      result[i] = clipper::datatypes::F_phi<dtype>( std::complex<dtype>(d1[i]) + std::complex<dtype>(d2[i]) );
  return result;
}

template<class dtype> HKL_data<F_phi<dtype> > operator -( const HKL_data<F_phi<dtype> >& d1, const HKL_data<F_phi<dtype> >& d2 )
{
  HKL_data<F_phi<dtype> > result( d1.base_hkl_info(), d1.base_cell() );
  for ( int i = 0; i < d1.base_hkl_info().num_reflections(); i++ )
    if ( !d1.missing(i) && !d2.missing(i) )
      result[i] = clipper::datatypes::F_phi<dtype>( std::complex<dtype>(d1[i]) - std::complex<dtype>(d2[i]) );
  return result;
}

template<class dtype> HKL_data<F_phi<dtype> > operator *( const HKL_data<F_phi<dtype> >& d1, const ftype& s )
{
  HKL_data<F_phi<dtype> > result( d1.base_hkl_info(), d1.base_cell() );
  for ( int i = 0; i < d1.base_hkl_info().num_reflections(); i++ )
    if ( !d1.missing(i) ) {
      result[i] = d1[i];
      result[i].scale(s);
    }
  return result;
}

template<class dtype> HKL_data<F_phi<dtype> > operator -( const HKL_data<F_phi<dtype> >& d1 )
{
  HKL_data<F_phi<dtype> > result( d1.base_hkl_info(), d1.base_cell() );
  for ( int i = 0; i < d1.base_hkl_info().num_reflections(); i++ )
    if ( !d1.missing(i) )
      result[i] = clipper::datatypes::F_phi<dtype>( -std::complex<dtype>(d1[i]) );
  return result;
}

template<class dtype> HKL_data<ABCD<dtype> > operator +( const HKL_data<ABCD<dtype> >& d1, const HKL_data<ABCD<dtype> >& d2 )
{
  HKL_data<ABCD<dtype> > result( d1.base_hkl_info(), d1.base_cell() );
  for ( int i = 0; i < d1.base_hkl_info().num_reflections(); i++ )
    if ( !d1.missing(i) && !d2.missing(i) ) {
      result[i].a() = d1[i].a() + d2[i].a();
      result[i].b() = d1[i].b() + d2[i].b();
      result[i].c() = d1[i].c() + d2[i].c();
      result[i].d() = d1[i].d() + d2[i].d();
    }
  return result;
}


// compile template types
template F_phi<ftype32> CLIPPER_IMEX operator +(const F_phi<ftype32>& d1, const F_phi<ftype32>& d2);
template F_phi<ftype32> CLIPPER_IMEX operator -(const F_phi<ftype32>& d1, const F_phi<ftype32>& d2);
template F_phi<ftype32> CLIPPER_IMEX operator -(const F_phi<ftype32>& d1);
template ABCD<ftype32> CLIPPER_IMEX operator +(const ABCD<ftype32>& d1, const ABCD<ftype32>& d2);

template F_phi<ftype64> CLIPPER_IMEX operator +(const F_phi<ftype64>& d1, const F_phi<ftype64>& d2);
template F_phi<ftype64> CLIPPER_IMEX operator -(const F_phi<ftype64>& d1, const F_phi<ftype64>& d2);
template F_phi<ftype64> CLIPPER_IMEX operator -(const F_phi<ftype64>& d1);
template ABCD<ftype64> CLIPPER_IMEX operator +(const ABCD<ftype64>& d1, const ABCD<ftype64>& d2);

template HKL_data<F_phi<ftype32> > CLIPPER_IMEX operator +(const HKL_data<F_phi<ftype32> >& d1, const HKL_data<F_phi<ftype32> >& d2);
template HKL_data<F_phi<ftype32> > CLIPPER_IMEX operator -(const HKL_data<F_phi<ftype32> >& d1, const HKL_data<F_phi<ftype32> >& d2);
template HKL_data<F_phi<ftype32> > CLIPPER_IMEX operator *(const HKL_data<F_phi<ftype32> >& d1, const ftype& s);
template HKL_data<F_phi<ftype32> > CLIPPER_IMEX operator *(const ftype& s, const HKL_data<F_phi<ftype32> >& d1);
template HKL_data<F_phi<ftype32> > CLIPPER_IMEX operator -(const HKL_data<F_phi<ftype32> >& d1);
template HKL_data<ABCD<ftype32> > CLIPPER_IMEX operator +(const HKL_data<ABCD<ftype32> >& d1, const HKL_data<ABCD<ftype32> >& d2);

template HKL_data<F_phi<ftype64> > CLIPPER_IMEX operator +(const HKL_data<F_phi<ftype64> >& d1, const HKL_data<F_phi<ftype64> >& d2);
template HKL_data<F_phi<ftype64> > CLIPPER_IMEX operator -(const HKL_data<F_phi<ftype64> >& d1, const HKL_data<F_phi<ftype64> >& d2);
template HKL_data<F_phi<ftype64> > CLIPPER_IMEX operator *(const HKL_data<F_phi<ftype64> >& d1, const ftype& s);
template HKL_data<F_phi<ftype64> > CLIPPER_IMEX operator *(const ftype& s, const HKL_data<F_phi<ftype64> >& d1);
template HKL_data<F_phi<ftype64> > CLIPPER_IMEX operator -(const HKL_data<F_phi<ftype64> >& d1);
template HKL_data<ABCD<ftype64> > CLIPPER_IMEX operator +(const HKL_data<ABCD<ftype64> >& d1, const HKL_data<ABCD<ftype64> >& d2);

} // namespace datatypes


template<int N> LogPhaseProb<N>::LogPhaseProb( const HKL_class& hkl_class )
{
  if ( hkl_class.centric() ) {
    pmin = Util::mod( Util::intr( N*hkl_class.allowed()/Util::twopi() ), N/2 );
    pinc = N/2;
    q.resize( 2, 0.0 );
  } else {
    pmin = 0;
    pinc = 1;
    q.resize( N, 0.0 );
  }
}

template<int N> template<class dtype> void LogPhaseProb<N>::set_abcd( const datatypes::ABCD<dtype>& abcd )
{
  if ( !abcd.missing() ) {
    ftype c, s;
    for ( int p = 0; p < q.size(); p++ ) {
      c = cos( phase(p) );
      s = sin( phase(p) );
      q[p] = abcd.a()*c + abcd.b()*s
	+ abcd.c()*(c*c-s*s) + abcd.d()*(2.0*c*s);
    }
  } else {
    for ( int p = 0; p < q.size(); p++ ) q[p] = 0.0;
  }
}

template<int N> template<class dtype> void LogPhaseProb<N>::get_abcd( datatypes::ABCD<dtype>& abcd ) const
{
  ftype q0, q1, c, s;
  q0 = 0.0;
  for ( int p = 0; p < q.size(); p++ ) q0 += q[p];
  q0 /= double( q.size() );
  abcd.a() = abcd.b() = abcd.c() = abcd.d() = 0.0;
  for ( int p = 0; p < q.size(); p++ ) {
    q1 = ( q[p] - q0 ) / double( (pinc==1) ? (N/2) : 2 );
    c = cos( phase(p) );
    s = sin( phase(p) );
    abcd.a() += q1 * c;
    abcd.b() += q1 * s;
    abcd.c() += q1 * (c*c-s*s);
    abcd.d() += q1 * (2.0*c*s);
  }
}

template<int N> template<class dtype> void LogPhaseProb<N>::set_phi_fom( const datatypes::Phi_fom<dtype>& phifom )
{
  if ( !phifom.missing() ) {
    dtype x = Util::min( phifom.fom(), dtype(0.999999) );
    if ( pinc != 1 ) x = Util::atanh( x );
    else             x = Util::invsim( x );
    for ( int p = 0; p < q.size(); p++ )
      q[p] = x*cos( phase(p) - phifom.phi() );
  } else {
    for ( int p = 0; p < q.size(); p++ ) q[p] = 0.0;
  }
}

template<int N> template<class dtype> void LogPhaseProb<N>::get_phi_fom( datatypes::Phi_fom<dtype>& phifom ) const
{
  Range<ftype64> qrange( -700.0, 700.0 );
  ftype s, a, b, q0, pq;
  s = a = b = q0 = 0.0;
  for ( int p = 0; p < q.size(); p++ ) q0 += q[p];
  q0 /= double( q.size() );
  for ( int p = 0; p < q.size(); p++ ) {
    pq = exp( qrange.truncate(q[p] - q0) );
    s += pq;
    a += pq * cos( phase(p) );
    b += pq * sin( phase(p) );
  }
  std::complex<ftype64> pw( a/s, b/s );
  phifom.phi() = std::arg( pw );
  phifom.fom() = std::abs( pw );
}

template class LogPhaseProb<24>;
template void LogPhaseProb<24>::set_abcd<ftype32>(const datatypes::ABCD<ftype32>& abcd);
template void LogPhaseProb<24>::get_abcd<ftype32>(datatypes::ABCD<ftype32>& abcd) const;
template void LogPhaseProb<24>::set_phi_fom<ftype32>( const datatypes::Phi_fom<ftype32>& phifom );
template void LogPhaseProb<24>::get_phi_fom<ftype32>( datatypes::Phi_fom<ftype32>& phifom ) const;
template void LogPhaseProb<24>::set_abcd<ftype64>(const datatypes::ABCD<ftype64>& abcd);
template void LogPhaseProb<24>::get_abcd<ftype64>(datatypes::ABCD<ftype64>& abcd) const;
template void LogPhaseProb<24>::set_phi_fom<ftype64>( const datatypes::Phi_fom<ftype64>& phifom );
template void LogPhaseProb<24>::get_phi_fom<ftype64>( datatypes::Phi_fom<ftype64>& phifom ) const;
template class LogPhaseProb<72>;
template void LogPhaseProb<72>::set_abcd<ftype32>(const datatypes::ABCD<ftype32>& abcd);
template void LogPhaseProb<72>::get_abcd<ftype32>(datatypes::ABCD<ftype32>& abcd) const;
template void LogPhaseProb<72>::set_phi_fom<ftype32>( const datatypes::Phi_fom<ftype32>& phifom );
template void LogPhaseProb<72>::get_phi_fom<ftype32>( datatypes::Phi_fom<ftype32>& phifom ) const;
template void LogPhaseProb<72>::set_abcd<ftype64>(const datatypes::ABCD<ftype64>& abcd);
template void LogPhaseProb<72>::get_abcd<ftype64>(datatypes::ABCD<ftype64>& abcd) const;
template void LogPhaseProb<72>::set_phi_fom<ftype64>( const datatypes::Phi_fom<ftype64>& phifom );
template void LogPhaseProb<72>::get_phi_fom<ftype64>( datatypes::Phi_fom<ftype64>& phifom ) const;
template class LogPhaseProb<180>;
template void LogPhaseProb<180>::set_abcd<ftype32>(const datatypes::ABCD<ftype32>& abcd);
template void LogPhaseProb<180>::get_abcd<ftype32>(datatypes::ABCD<ftype32>& abcd) const;
template void LogPhaseProb<180>::set_phi_fom<ftype32>( const datatypes::Phi_fom<ftype32>& phifom );
template void LogPhaseProb<180>::get_phi_fom<ftype32>( datatypes::Phi_fom<ftype32>& phifom ) const;
template void LogPhaseProb<180>::set_abcd<ftype64>(const datatypes::ABCD<ftype64>& abcd);
template void LogPhaseProb<180>::get_abcd<ftype64>(datatypes::ABCD<ftype64>& abcd) const;
template void LogPhaseProb<180>::set_phi_fom<ftype64>( const datatypes::Phi_fom<ftype64>& phifom );
template void LogPhaseProb<180>::get_phi_fom<ftype64>( datatypes::Phi_fom<ftype64>& phifom ) const;
template class LogPhaseProb<360>;
template void LogPhaseProb<360>::set_abcd<ftype32>(const datatypes::ABCD<ftype32>& abcd);
template void LogPhaseProb<360>::get_abcd<ftype32>(datatypes::ABCD<ftype32>& abcd) const;
template void LogPhaseProb<360>::set_phi_fom<ftype32>( const datatypes::Phi_fom<ftype32>& phifom );
template void LogPhaseProb<360>::get_phi_fom<ftype32>( datatypes::Phi_fom<ftype32>& phifom ) const;
template void LogPhaseProb<360>::set_abcd<ftype64>(const datatypes::ABCD<ftype64>& abcd);
template void LogPhaseProb<360>::get_abcd<ftype64>(datatypes::ABCD<ftype64>& abcd) const;
template void LogPhaseProb<360>::set_phi_fom<ftype64>( const datatypes::Phi_fom<ftype64>& phifom );
template void LogPhaseProb<360>::get_phi_fom<ftype64>( datatypes::Phi_fom<ftype64>& phifom ) const;


} // namespace clipper
