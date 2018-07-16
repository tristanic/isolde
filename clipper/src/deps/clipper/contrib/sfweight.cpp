/* sfweight.cpp: structure factor weighting implementation */
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


#include "sfweight.h"
#include "../core/hkl_operators.h"
#include "../core/resol_targetfn.h"


namespace clipper {


template<class T> typename SFweight_spline<T>::TargetResult SFweight_spline<T>::TargetFo::operator() ( const HKL_class cls, const datatypes::F_sigF<T>& fo0, const datatypes::ABCD<T>& hl0, const datatypes::F_phi<T>& fc0, const ftype& s, const ftype& w, const std::vector<HLterms>& hlterms )
{
  ftype fo(0.0), fc(0.0), phic(0.0), sfo(1.0);
  if ( !fo0.missing() ) { fo = fo0.f(); sfo = fo0.sigf(); }
  if ( !fc0.missing() ) { fc = fc0.f(); phic = fc0.phi(); }
  const ftype epc = cls.epsilonc();
  const ftype s2 = s*s;
  const ftype fo2 = fo*fo;
  const ftype fc2 = fc*fc;
  const ftype d = 2.0*sfo*sfo + epc*w;
  const ftype d2 = d*d;
  const ftype d3 = d*d2;
  const ftype d4 = d*d3;
  const ftype x = 2.0*fo*fc*s/d;
  ftype i0, di0, ddi0, cf;
  TargetResult r;
  if ( cls.centric() ) {
    i0 = (fabs(x)<10.0) ? log(cosh(x)) : fabs(x)+log(0.5);
    di0 = tanh(x);
    ddi0 = 1.0-pow(tanh(x),2);
    cf = 0.5;
  } else {
    i0 = Util::sim_integ(x);
    di0 = Util::sim(x);
    ddi0 = Util::sim_deriv(x);
    cf = 1.0;
  }
  r.r = cf*log(d) + (fo2+s2*fc2)/d - i0;
  r.ds = 2.0*s*fc2/d - (2.0*fo*fc/d)*di0;
  r.dw = epc*( cf/d - (fo2+s2*fc2)/d2 + (2.0*fo*fc*s/d2)*di0 );
  r.dss = 2.0*fc2/d - (4.0*fo2*fc2/d2)*ddi0;
  r.dww = epc*epc*( -cf/d2 + 2.0*(fo2+s2*fc2)/d3
      	      - (4.0*fo*fc*s/d3)*di0 - (4.0*fo2*fc2*s2/d4)*ddi0 );
  r.dsw = epc*( -2.0*s*fc2/d2 + (2.0*fo*fc/d2)*di0
      	  + (4.0*fo2*fc2*s/d3)*ddi0 );
  abcd = datatypes::ABCD<T>( x*cos(phic), x*sin(phic), 0.0, 0.0 );
  phiw = datatypes::Phi_fom<T>( phic, di0 );
  return r;
}


template<class T> typename SFweight_spline<T>::TargetResult SFweight_spline<T>::TargetHL::operator() ( const HKL_class cls, const datatypes::F_sigF<T>& fo0, const datatypes::ABCD<T>& hl0, const datatypes::F_phi<T>& fc0, const ftype& s, const ftype& w, const std::vector<HLterms>& hlterms )
{
  ftype fo(0.0), fc(0.0), phic(0.0), sfo(1.0);
  ftype a0(0.0), b0(0.0), c0(0.0), d0(0.0);
  if ( !fo0.missing() ) { fo = fo0.f(); sfo = fo0.sigf(); }
  if ( !fc0.missing() ) { fc = fc0.f(); phic = fc0.phi(); }
  if ( !hl0.missing() ) { a0=hl0.a(); b0=hl0.b(); c0=hl0.c(); d0=hl0.d(); }
  const ftype epc = cls.epsilonc();
  const ftype s2 = s*s;
  const ftype fo2 = fo*fo;
  const ftype fc2 = fc*fc;
  const ftype d = 2.0*sfo*sfo + epc*w;
  const ftype d2 = d*d;
  const ftype d3 = d*d2;
  //const ftype d4 = d*d3;
  const ftype epcd = epc/d;
  const ftype sf2 = fo2 + s2*fc2;
  const ftype xs = 2.0*fo*fc/d;
  const ftype cosc = cos(phic);
  const ftype sinc = sin(phic);
  const ftype hl_a = a0 + xs*s*cosc;
  const ftype hl_b = b0 + xs*s*sinc;
  const ftype hl_c = c0;
  const ftype hl_d = d0;
  ftype cf = 1.0;
  int a_zero = 0;
  int a_step = 1;
  if ( cls.centric() ) {
    a_step = hlterms.size()/2;
    a_zero = Util::intr( ftype(hlterms.size())*cls.allowed()/Util::twopi() );
    a_zero = Util::mod( a_zero, a_step );
    cf = 0.5;
  }
  ftype an, asum, asum_ds, asum_dss, asum_dw, asum_dww, asum_dsw, asum_a, asum_b;
  an = asum = asum_ds = asum_dss = asum_dw = asum_dww = asum_dsw = asum_a = asum_b = 0.0;
  ftype q, q1, q2, e;
  ftype qmax = sqrt( hl_a*hl_a + hl_b*hl_b );
  for ( int a = a_zero; a < hlterms.size(); a += a_step ) {
    const HLterms& trig = hlterms[a];
    q = hl_a*trig.cosa + hl_b*trig.sina + hl_c*trig.cos2a + hl_d*trig.sin2a;
    q1 = xs * ( cosc*trig.cosa + sinc*trig.sina );
    q2 = s * q1;
    e = exp( q - qmax );
    an       += 1.0;
    asum     += e;
    asum_ds  += e * q1;
    asum_dss += e * q1*q1;
    asum_dw  += e * ( -q2 )*epcd;
    asum_dww += e * ( 2.0 + q2 )*q2*epcd*epcd;
    asum_a   += e * trig.cosa;
    asum_b   += e * trig.sina;
    // asum_dsw += ?;
  }
  asum_a   /= asum;
  asum_b   /= asum;
  asum     /= an;
  asum_ds  /= an;
  asum_dss /= an;
  asum_dw  /= an;
  asum_dww /= an;
  abcd = datatypes::ABCD<T>( hl_a, hl_b, hl_c, hl_d );
  phiw = datatypes::Phi_fom<T>( atan2( asum_b, asum_a ),
				sqrt( asum_a*asum_a + asum_b*asum_b ) );
  TargetResult r;
  r.r   =           cf*log(d) + sf2/d - log( asum ) - qmax;
  r.ds  =               (2.0*s*fc2)/d - asum_ds/asum;
  r.dw  =           epc*(cf/d-sf2/d2) - asum_dw/asum;
  r.dss =                 (2.0*fc2)/d - asum_dss/asum + Util::sqr(asum_ds/asum);
  r.dww = epc*epc*(-cf/d2+2.0*sf2/d3) - asum_dww/asum + Util::sqr(asum_dw/asum);
  r.dsw = Util::nan();
  return r;
}


template<class T> SFweight_spline<T>::SFweight_spline( HKL_data<datatypes::F_phi<T> >& fb, HKL_data<datatypes::F_phi<T> >& fd, HKL_data<datatypes::Phi_fom<T> >& phiw, const HKL_data<datatypes::F_sigF<T> >& fo, const HKL_data<datatypes::F_phi<T> >& fc, const HKL_data<datatypes::Flag>& usage, const int n_reflns, const int n_params )
{
  init( n_reflns, n_params );
  (*this)( fb, fd, phiw, fo, fc, usage );
}


template<class T> void SFweight_spline<T>::init( const int n_reflns, const int n_params, const int n_phases ) {
  nreflns = n_reflns;
  nparams = n_params;
  hlterms.resize( n_phases );
  for ( int a = 0; a < hlterms.size(); a++ ) {
    ftype phi = ( Util::twopi() * ftype(a)/ftype(hlterms.size()) );
    hlterms[a].cosa  = cos(phi);
    hlterms[a].sina  = sin(phi);
    hlterms[a].cos2a = cos(2.0*phi);
    hlterms[a].sin2a = sin(2.0*phi);
  }
}


template<class T> bool SFweight_spline<T>::operator() ( HKL_data<datatypes::F_phi<T> >& fb, HKL_data<datatypes::F_phi<T> >& fd, HKL_data<datatypes::Phi_fom<T> >& phiw, const HKL_data<datatypes::F_sigF<T> >& fo0, const HKL_data<datatypes::F_phi<T> >& fc0, const HKL_data<datatypes::Flag>& usage )
{
  TargetFo tgt;
  HKL_data<datatypes::ABCD<T> > hl0( fo0.hkl_info() ), hl( fo0.hkl_info() );
  return evaluate( fb, fd, phiw, hl, fo0, hl0, fc0, usage, tgt );
}


template<class T> bool SFweight_spline<T>::operator() ( HKL_data<datatypes::F_phi<T> >& fb, HKL_data<datatypes::F_phi<T> >& fd, HKL_data<datatypes::Phi_fom<T> >& phiw, HKL_data<datatypes::ABCD<T> >& hl, const HKL_data<datatypes::F_sigF<T> >& fo0, const HKL_data<datatypes::ABCD<T> >& hl0, const HKL_data<datatypes::F_phi<T> >& fc0, const HKL_data<datatypes::Flag>& usage )
{
  TargetHL tgt;
  return evaluate( fb, fd, phiw, hl, fo0, hl0, fc0, usage, tgt );
}


template<class T> template<class F> bool SFweight_spline<T>::evaluate( HKL_data<datatypes::F_phi<T> >& fb, HKL_data<datatypes::F_phi<T> >& fd, HKL_data<datatypes::Phi_fom<T> >& phiw, HKL_data<datatypes::ABCD<T> >& hl, const HKL_data<datatypes::F_sigF<T> >& fo0, const HKL_data<datatypes::ABCD<T> >& hl0, const HKL_data<datatypes::F_phi<T> >& fc0, const HKL_data<datatypes::Flag>& usage, F tgt )
{
  const HKL_info& hkls = fo0.base_hkl_info();
  typedef clipper::HKL_info::HKL_reference_index HRI;
  bool status = false;

  // count reflections and determine number of params
  HKL_data<datatypes::Flag_bool> flag(hkls);
  for ( HRI ih = flag.first(); !ih.last(); ih.next() )
    flag[ih].flag() = (!fo0[ih].missing()) && (!fc0[ih].missing()) &&
      (usage[ih].flag()!=SFweight_base<T>::NONE);
  int npar_sig = num_params( flag );
  int npar_scl = npar_sig;
  while ( npar_scl < 12 ) npar_scl *= 2;

  // prepare function
  BasisFn_spline basisfn( flag, npar_sig, 1.0 );
  BasisFn_spline basissc( flag, npar_scl, 1.0 );
  BasisFn_base::Fderiv dsdp, dwdp;
  TargetResult fn;
  scale_fo.resize( hkls.num_reflections() );
  scale_fc.resize( hkls.num_reflections() );
  value_s.resize( hkls.num_reflections() );
  value_w.resize( hkls.num_reflections() );

  // create E's for scaling
  HKL_data<datatypes::E_sigE<T> > eo( hkls ), ec( hkls );
  for ( HRI ih = fo0.first(); !ih.last(); ih.next() ) {
    eo[ih].E() = fo0[ih].f() / sqrt( ih.hkl_class().epsilon() );
    ec[ih].E() = fc0[ih].f() / sqrt( ih.hkl_class().epsilon() );
    eo[ih].sigE() = ec[ih].sigE() = 1.0;
  }
  // calc scale
  std::vector<ftype> param_fo( basissc.num_params(), 1.0 );
  TargetFn_scaleEsq<datatypes::E_sigE<T> > tfn_fo( eo );
  ResolutionFn rfn_fo( hkls, basissc, tfn_fo, param_fo );
  param_fo = rfn_fo.params();
  std::vector<ftype> param_fc( basissc.num_params(), 1.0 );
  TargetFn_scaleEsq<datatypes::E_sigE<T> > tfn_fc( ec );
  ResolutionFn rfn_fc( hkls, basissc, tfn_fc, param_fc );
  param_fc = rfn_fc.params();

  // prescale Fo, Fc
  HKL_data<datatypes::F_sigF<T> > fo = fo0;
  HKL_data<datatypes::F_phi<T> >  fc = fc0;
  for ( HRI ih = fo0.first(); !ih.last(); ih.next() ) {
    scale_fo[ ih.index() ] = sqrt( basissc.f_s( ih.invresolsq(), param_fo ) );
    scale_fc[ ih.index() ] = sqrt( basissc.f_s( ih.invresolsq(), param_fc ) );
    fo[ih].scale( scale_fo[ ih.index() ] );
    fc[ih].scale( scale_fc[ ih.index() ] );
  }

  // make first estimate of s
  param_s = std::vector<ftype>( npar_sig, 1.0 );

  // make first estimate of w
  HKL_data<datatypes::F_sigF<T> > ftmp( hkls );
  for ( HRI ih = flag.first(); !ih.last(); ih.next() ) if ( flag[ih].flag() )
    ftmp[ih].f() = ftmp[ih].sigf() =
      pow( fo[ih].f() - sqrt(basisfn.f_s(ih.invresolsq(),param_s))*fc[ih].f(), 2.0 ) / ih.hkl_class().epsilonc();
  TargetFn_meanFnth<datatypes::F_sigF<T> > target_w( ftmp, 1.0 );
  ResolutionFn rfn2( hkls, basisfn, target_w, param_w );
  param_w = rfn2.params();
  //for ( int i = 0; i < npar_sig; i++ ) std::cout << i << " " << param_s[i] << "    \t" << param_w[i] << "\n";

  // smooth the error term
  for ( int i = 0; i < npar_sig-1; i++ )
    param_w[i] = Util::max( param_w[i], 0.5*param_w[i+1] );
  //for ( int i = 0; i < npar_sig; i++ ) std::cout << i << " " << param_s[i] << "    \t" << param_w[i] << "\n";

  ftype ll, ll_old = 1.0e30;
  // now 25 cycles to refine s and w
  int c = 0, clim = 25;
  for ( c = 0; c < clim; c++ ) {
    std::vector<ftype> grad_s( npar_sig, 0.0 ), shft_s( npar_sig, 0.0 );
    std::vector<ftype> grad_w( npar_sig, 0.0 ), shft_w( npar_sig, 0.0 );
    Matrix<ftype> curv_s( npar_sig, npar_sig, 0.0 );
    Matrix<ftype> curv_w( npar_sig, npar_sig, 0.0 );
    ll = 0.0;

    // build matrices
    for ( HRI ih = flag.first(); !ih.last(); ih.next() ) if ( flag[ih].flag() ) {
      dsdp = basisfn.fderiv_s( ih.invresolsq(), param_s );
      dwdp = basisfn.fderiv_s( ih.invresolsq(), param_w );
      fn = tgt( ih.hkl_class(), fo[ih], hl0[ih], fc[ih], dsdp.f, dwdp.f, hlterms );
      //if ( Util::isnan(fn.r) ) std::cout << ih.hkl().format() << fo[ih].f() << " " << fo[ih].sigf() << " " << fc[ih].f() << " " << fc[ih].missing() << flag[ih].missing() << " " << hl0[ih].a() << " " << hl0[ih].b() << " " << hl0[ih].c() << " " << hl0[ih].d() << " " << dsdp.f << " " << dwdp.f << " \tFo,SIGo,Fc,missing\n";
      ll += fn.r;
      for ( int i = 0; i < npar_sig; i++ ) {
	grad_s[i] += fn.ds * dsdp.df[i];
	grad_w[i] += fn.dw * dwdp.df[i];
	//for ( int j = 0; j < npar_sig; j++ ) {
	for ( int j = Util::max(i-1,0); j <= Util::min(i+1,npar_sig-1); j++ ) {
	  curv_s(i,j) += dsdp.df2(i,j)*fn.ds + dsdp.df[i]*dsdp.df[j]*fn.dss;
	  curv_w(i,j) += dwdp.df2(i,j)*fn.dw + dwdp.df[i]*dwdp.df[j]*fn.dww;
	}
      }
    }
    //std::cout << c << "\t" << ll << "\n";

    if ( ll > ll_old ) break;  // break on divergence

    shft_s = curv_s.solve( grad_s );
    shft_w = curv_w.solve( grad_w );
    for ( int i = 0; i < npar_sig; i++ ) {
      //std::cout << i << "   \t" << param_s[i] << "   \t" << param_w[i] << "   \t" << shft_s[i] << "   \t" << shft_w[i] << "\n";
      // soft buffers to prevent negatives
      param_s[i] -= Util::min( shft_s[i], 0.25*param_s[i] );
      param_w[i] -= Util::min( shft_w[i], 0.25*param_w[i] );
    }

    if ( ll / ll_old > 0.999999 ) { status=true; break; }  // break on convergence
    ll_old = ll;
  }

  // store s,w
  for ( HRI ih = fo0.first(); !ih.last(); ih.next() ) {
    value_s[ih.index()] = basisfn.f_s( ih.invresolsq(), param_s );
    value_w[ih.index()] = basisfn.f_s( ih.invresolsq(), param_w );
  }

  // calculate map coeffs and FOM
  reevaluate( fb, fd, phiw, hl, fo0, hl0, fc0, usage, tgt );

  return status;
}


template<class T> template<class F> bool SFweight_spline<T>::reevaluate( HKL_data<datatypes::F_phi<T> >& fb, HKL_data<datatypes::F_phi<T> >& fd, HKL_data<datatypes::Phi_fom<T> >& phiw, HKL_data<datatypes::ABCD<T> >& hl, const HKL_data<datatypes::F_sigF<T> >& fo0, const HKL_data<datatypes::ABCD<T> >& hl0, const HKL_data<datatypes::F_phi<T> >& fc0, const HKL_data<datatypes::Flag>& usage, F tgt )
{
  typedef clipper::HKL_info::HKL_reference_index HRI;

  // prepare function
  TargetResult fn;
  datatypes::F_sigF<T> fo;
  datatypes::F_phi<T> fc, twomfo, mfo, dfc, fzero(0.0,0.0);

  // calculate map coeffs and FOM
  llw = llf = 0.0;
  for ( HRI ih = fo0.first(); !ih.last(); ih.next() ) {
    // prescale Fo, Fc
    fo = fo0[ih];
    fc = fc0[ih];
    fo.scale( scale_fo[ih.index()] );
    fc.scale( scale_fc[ih.index()] );
    // get params an llk
    const ftype s = value_s[ih.index()];
    const ftype w = value_w[ih.index()];
    fn = tgt( ih.hkl_class(), fo, hl0[ih], fc, s, w, hlterms );
    hl[ih] = tgt.abcd;
    phiw[ih] = tgt.phiw;
    mfo    = datatypes::F_phi<T>( tgt.phiw.fom() * fo.f(), tgt.phiw.phi() );
    twomfo = datatypes::F_phi<T>( 2.0 * mfo.f(), mfo.phi() );
    dfc    = datatypes::F_phi<T>( s * fc.f(), fc.phi() );
    if ( (!fo.missing()) && (!fc.missing()) ) {
      if      ( usage[ih].flag()==SFweight_base<T>::BOTH ) llw += fn.r;
      else if ( usage[ih].flag()==SFweight_base<T>::NONE ) llf += fn.r;
    }
    // deal with all possibilities of missing and non-missing
    if ( !fc.missing() ) {
      if ( !fo.missing() ) {
	fb[ih] = twomfo - dfc;
	fd[ih] =    mfo - dfc;
      } else {
	fb[ih] = dfc;
	fd[ih] = fzero;
      }
    } else {
      if ( !fo.missing() ) {
	fb[ih] = mfo;
	fd[ih] = fzero;
      } else {
	fb[ih] = fzero;
	fd[ih] = fzero;
      }
    }
    // undo scaling on fb, fd
    fb[ih].scale( 1.0/scale_fo[ih.index()] );
    fd[ih].scale( 1.0/scale_fo[ih.index()] );
  }

  return true;
}


template<class T> typename SFweight_spline<T>::TargetResult SFweight_spline<T>::targetfn( const HKL_class cls, const datatypes::F_sigF<T>& fo0, const datatypes::F_phi<T>& fc0, const ftype& s, const ftype& w ) const
{
  const datatypes::ABCD<T> hl0;
  TargetFo tgt;
  return tgt( cls, fo0, hl0, fc0, s, w, hlterms );
}


template<class T> typename SFweight_spline<T>::TargetResult SFweight_spline<T>::targethl( const HKL_class cls, const datatypes::F_sigF<T>& fo0, const datatypes::ABCD<T>& hl0, const datatypes::F_phi<T>& fc0, const ftype& s, const ftype& w ) const
{
  TargetHL tgt;
  return tgt( cls, fo0, hl0, fc0, s, w, hlterms );
}


template<class T> int SFweight_spline<T>::num_params( const HKL_data<datatypes::Flag_bool>& flag ) const
{
  int npar;
  int n_use = flag.num_obs();
  if ( nparams == 0 ) {
    npar = Util::max( n_use / nreflns, 2 );
  } else if ( nreflns == 0 ) {
    npar = nparams;
  } else {
    ftype np1 = ftype(nparams+0.499);
    ftype np2 = ftype(n_use) / ftype(nreflns);
    ftype np = sqrt( np1*np1*np2*np2 / ( np1*np1+np2*np2 ) );
    npar = Util::max( int(np), 2 );
  }
  return npar;
}


template<class T> void SFweight_spline<T>::debug() const
{
  TargetResult r00, r01, r10, r11, rxx;
  Spacegroup p1( Spacegroup::P1 );
  HKL_class cls;
  datatypes::F_sigF<T> fo;
  datatypes::F_phi<T>  fc;
  fo.f() = 10.0; fo.sigf() = 2.0;
  fc.f() = 15.0; fc.phi() = 0.0;
  fo.sigf() = 0.0; //!!!!!
  ftype ds = 0.000001;
  ftype dw = 0.000001;
  datatypes::ABCD<T> hl( 0.0, 0.0, 0.0, 0.0 );
  for ( int h = 0; h < 2; h++ ) {
    cls = HKL_class( p1, HKL( h, 0, 0 ) );
    std::cout << "\nCentric? " << cls.centric() << " epsc " << cls.epsilonc() << "\n";
    for ( ftype w = 10.0; w < 1000.0; w *= 3.0 )
      for ( ftype s = 0.4; s < 2.0; s *= 2.0 ) {

        rxx = targethl( cls, fo, hl, fc, s, w );
        r00 = targetfn( cls, fo, fc, s, w );
        r01 = targetfn( cls, fo, fc, s, w+dw );
        r10 = targetfn( cls, fo, fc, s+ds, w );
        r11 = targetfn( cls, fo, fc, s+ds, w+dw );

        std::cout << w << " " << s << "\t" << r00.r << " " << r01.r << " " << r10.r << " " << r11.r << " " << rxx.r << "\n";
        std::cout << (r10.r-r00.r)/ds << "\t" << r00.ds << "\n";
        std::cout << (r01.r-r00.r)/dw << "\t" << r00.dw << "\n";
        std::cout << (r10.ds-r00.ds)/ds << "\t" << r00.dss << "\n";
        std::cout << (r01.dw-r00.dw)/dw << "\t" << r00.dww << "\n";
        std::cout << (r01.ds-r00.ds)/dw << "\t" << r00.dsw << "\n";
        std::cout << (r10.dw-r00.dw)/ds << "\t" << r00.dsw << "\n";
      }
  }
  for ( int h = 0; h < 2; h++ ) {
    cls = HKL_class( p1, HKL( h, 0, 0 ) );
    std::cout << "\nCentric? " << cls.centric() << " epsc " << cls.epsilonc() << "\n";
    for ( ftype w = 10.0; w < 1000.0; w *= 3.0 )
      for ( ftype s = 0.4; s < 2.0; s *= 2.0 ) {

        rxx = targetfn( cls, fo, fc, s, w );
        r00 = targethl( cls, fo, hl, fc, s, w );
        r01 = targethl( cls, fo, hl, fc, s, w+dw );
        r10 = targethl( cls, fo, hl, fc, s+ds, w );
        r11 = targethl( cls, fo, hl, fc, s+ds, w+dw );

        std::cout << w << " " << s << "\t" << r00.r << " " << r01.r << " " << r10.r << " " << r11.r << " " << rxx.r << "\n";
        std::cout << (r10.r-r00.r)/ds << "\t" << r00.ds << "\n";
        std::cout << (r01.r-r00.r)/dw << "\t" << r00.dw << "\n";
        std::cout << (r10.ds-r00.ds)/ds << "\t" << r00.dss << "\n";
        std::cout << (r01.dw-r00.dw)/dw << "\t" << r00.dww << "\n";
        std::cout << (r01.ds-r00.ds)/dw << "\t" << r00.dsw << "\n";
        std::cout << (r10.dw-r00.dw)/ds << "\t" << r00.dsw << "\n";
      }
  }
}


// compile templates

template class CLIPPER_IMEX SFweight_spline<ftype32>;

template class CLIPPER_IMEX SFweight_spline<ftype64>;


} // namespace clipper
