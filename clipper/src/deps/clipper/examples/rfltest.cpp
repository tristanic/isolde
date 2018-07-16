#include <clipper/clipper.h>
#include <clipper/core/fftmap_sparse.h>

#include <stdlib.h>
#include <iostream>

#define T(X) (fmod(Util::rad2d(X)+72000.0,360.0))


using namespace clipper;
using namespace clipper::data32;

int main()
{
  int sg_num,i,n;

  std::cout << "Enter spacegroup: ";
  std::cin >> sg_num;

  CSpacegroup cspgr("spgr",Spacegroup(Spgr_descr(sg_num)));
  CCell ccell(cspgr,"cell",Cell(Cell_descr(100.0,100.0,100.0,Util::pi()/2.0,Util::pi()/2.0,Util::pi()/2.0)));
  CCell ccell2(cspgr,"cell2",Cell(Cell_descr(99.0,99.0,99.0,Util::pi()/2.0,Util::pi()/2.0,Util::pi()/2.0)));
  cspgr.Spacegroup::debug(); ccell.Cell::debug();
  CHKL_info chkl(ccell,"hkl");

  chkl.init( cspgr, ccell, Resolution(19.99) );
  chkl.generate_hkl_list();
  chkl.spacegroup().debug();

  CHKL_data<F_phi> cfphi(chkl, "myfphi");
  CHKL_data<F_phi> cfphi2(chkl, "myfphi2");
  CHKL_data<F_sigF> cfsig(chkl, "myfsig");

  cspgr.Container::debug();

  std::cout << dynamic_cast<const HKL_data_base*>(chkl.find_path_ptr("myfphi"))->type() << "\n";
  std::cout << dynamic_cast<const HKL_data_base*>(chkl.find_path_ptr("myfsig"))->type() << "\n";

  F_sigF dat1;
  F_phi  dat2;
  dat1.f()=dat1.sigf()=0.0f;
  dat2.f()=dat2.phi()=0.0f;
  n = chkl.spacegroup().num_symops(); std::cout << n << " ";
  n = Util::min( n-1, 1 ); std::cout << n << "\n";
  for (i=0; i<chkl.num_reflections(); i++) {
    dat1.f()=i;
    cfsig[i]=dat1;
    dat2.f()=i; dat2.phi()=(23*i)*Util::twopi()/40.0;
    if ( chkl.hkl_of(i) == HKL(0,0,0) ) dat2.phi() = 0.0;
    cfphi[i]=dat2;
    HKL rfl=-chkl.hkl_of(i).transform(cspgr.symop(n));
    std::cout << cfphi[i].f() << "=" << (cfphi[rfl]).f() << "\n";
  }


  for (int h=-1; h<=1; h++)
    for (int k=-1; k<=1; k++)
      for (int l=-1; l<=1; l++)
	if (!cspgr.hkl_class(HKL(h,k,l)).sys_abs()) {
          std::cout << " (" << h << "," << k << "," << l << ") " << " ";
          std::cout <<  (cfphi[HKL(h,k,l)].f()) << " ";
          std::cout << T(cfphi[-HKL(h,k,l).transform(cspgr.symop(n))].phi()) << "\n";
        }


  // now make an fftmap

  Grid_sampling cgrid(24,24,24);

  FFTmap fftmap( cspgr, ccell, cgrid );
  FFTmap_p1 fftmapp1( cgrid );

  for (i=0; i<chkl.num_reflections(); i++)
    fftmap.set_recip_data( chkl.hkl_of(i), cfphi[i] );

  fftmap.fft_h_to_x();
  std::cout << "done fft\n";
  for ( Coord_grid c(0,0,0); !c.last(cgrid); c.next(cgrid) )
    fftmapp1.real_data(c) = fftmap.get_real_data(c);  // copy to fftmap_p1
  fftmapp1.fft_x_to_h(ccell.volume());
  fftmapp1.fft_h_to_x(1.0/ccell.volume());
  fftmapp1.fft_x_to_h(ccell.volume());

  FFTmap_sparse_p1_hx fftmaps1( cgrid );
  FFTmap_sparse_p1_xh fftmaps2( cgrid );
  for (i=0; i<chkl.num_reflections(); i++)
    fftmaps1.set_hkl( chkl.hkl_of(i), std::complex<ffttype>(cfphi[i]) );
  Grid_range mgrid( Coord_grid(-1,-1,-1), Coord_grid(1,1,1) );
  for ( Coord_grid c = mgrid.min(); !c.last(mgrid); c.next(mgrid) )
    fftmaps1.require_real_data(c.unit(cgrid));
  fftmaps1.fft_h_to_x(1.0/ccell.volume());
  for ( Coord_grid c = mgrid.min(); !c.last(mgrid); c.next(mgrid) )
    std::cout << c.format() << "\t" << fftmaps1.real_data(c.unit(cgrid)) << "\t" << fftmap.get_real_data(c) << "\n";
  fftmapp1.reset();
  for ( Coord_grid c = mgrid.min(); !c.last(mgrid); c.next(mgrid) )
    fftmapp1.real_data(c.unit(cgrid)) = fftmaps1.real_data(c.unit(cgrid));
  for ( Coord_grid c = mgrid.min(); !c.last(mgrid); c.next(mgrid) )
    fftmaps2.real_data(c.unit(cgrid)) = fftmaps1.real_data(c.unit(cgrid));
  for (i=0; i<chkl.num_reflections(); i++)
    fftmaps2.require_hkl( chkl.hkl_of(i) );
  fftmaps2.fft_x_to_h(ccell.volume());
  fftmapp1.fft_x_to_h(ccell.volume());
  for (i=0; i<chkl.num_reflections(); i++)
    std::cout << chkl.hkl_of(i).format() << "\t" << fftmapp1.get_hkl( chkl.hkl_of(i) ) << "\t" << fftmaps2.get_hkl( chkl.hkl_of(i) ) << "\n";

  fftmap.fft_x_to_h();

  Xmap<float> xmap( cspgr, ccell, cgrid );
  xmap.fft_from( cfphi );
  std::cout << "done fft\n";
  xmap.fft_to  ( cfphi2 );

  // for (i=0; i<chkl.num_reflections(); i++) cfphi[i] = fftmapp1.get_hkl( chkl.hkl_of(i) );

  std::cout.precision(3);
  for (int h=-2; h<=2; h++)
    for (int k=-2; k<=2; k++)
      for (int l=-2; l<=2; l++) {
	HKL rfl(h,k,l);
	if ((!cspgr.hkl_class(rfl).sys_abs())&&(!cspgr.hkl_class(rfl).centric())) {
          std::cout << rfl.format() << " ";
          std::cout.width(4); std::cout <<  (cfphi[rfl].f()) << " ";
          std::cout.width(4); std::cout << T(cfphi[rfl].phi()) << ":";
          std::cout.width(4); std::cout <<  (cfphi2[rfl].f()) << " ";
          std::cout.width(4); std::cout << T(cfphi2[rfl].phi()) << ":";
	  std::cout.width(4); std::cout << fftmap.get_recip_data( rfl ).f() << " ";
	  std::cout.width(4); std::cout << T(fftmap.get_recip_data( rfl ).phi()) << "\n";
        }
      }

  for ( double x=0.1; x<50; x*=1.5 )
    std::cout << x << " " << Util::sim(x) << " " << Util::invsim(Util::sim(x)) << " " << Util::sim(-x) << " " << Util::invsim(Util::sim(-x)) << " " << (Util::bessel_i0(x+0.0001)/Util::bessel_i0(x)-1.0)/0.0001 << " " << (Util::sim_integ(x+0.0001)-Util::sim_integ(x))/0.0001 << "\n";
  for ( double x=0.1; x>0.0001; x/=3.1 )
    std::cout << x << " " << Util::sim_deriv(x) << " " << Util::sim(x) << " " << Util::sim_integ(x) << " " << Util::sim_deriv(-x) << " " << Util::sim(-x) << " " << Util::sim_integ(-x) << "\n";


  CHKL_data<Phi_fom> cphifom(chkl, "myphifom");
  CHKL_data<ABCD> cabcd(chkl, "myABCD");
  CHKL_data<Phi_fom> cphifom2(chkl, "myphifom2");
  Phi_fom phifom;
  for (i=0; i<chkl.num_reflections(); i++) {
    if ( chkl.hkl_class(i).centric() )
      phifom.phi() = chkl.hkl_class(i).allowed() + Util::pi()*ftype(i%2);
    else
      phifom.phi() = Util::pi() * float(i%40)/20;
    phifom.fom() = 0.1 * float(i%10) + 0.05;
    cphifom[i] = phifom;
  }
  cabcd.compute( cphifom, Compute_abcd_from_phifom() );
  cphifom2.compute( cabcd, Compute_phifom_from_abcd() );
  HKL_info::HKL_reference_index ih;
  for ( ih = cphifom2.first(); !ih.last(); ih.next() ) {
    LogPhaseProb<24> q(ih.hkl_class());
    q.set_phi_fom( cphifom2[ih] );
    q.get_abcd( cabcd[ih] );
    q.set_abcd( cabcd[ih] );
    q.get_phi_fom( cphifom2[ih] );
  }
  for (i=0; i<chkl.num_reflections(); i++) {
    std::cout << i << chkl.hkl_of(i).format() << " : " << cphifom[i].phi() << " " << cphifom[i].fom() << " : " << cabcd[i].a() << " " << cabcd[i].b() << " : " << cphifom2[i].phi() << " " << cphifom2[i].fom() << "\n";
  }

  cfphi.compute( cfsig, cphifom, Compute_fphi_from_fsigf_phifom() );

  // Now test HKL_sampling
  for ( double a = 20.0; a < 22.0; a+=0.01 ) {
    clipper::Cell cell( clipper::Cell_descr( a, 24.0, 18.0, 60.0, 70.0, 80.0 ) );
    clipper::Resolution reso( 2.0 );
    clipper::HKL_sampling hklsam( cell, reso ); 
    HKL lim = hklsam.hkl_limit();
    for ( int ih = -20; ih <= 20; ih++ )
      for ( int ik = -20; ik <= 20; ik++ )
	for ( int il = -20; il <= 20; il++ ) {
	  HKL hkl(ih,ik,il);
	  if ( hklsam.in_resolution(hkl) ^ ( hkl.invresolsq(cell) < reso.invresolsq_limit() ) ) std::cout << "Err " << hkl.format() << hklsam.in_resolution(hkl) << ( hkl.invresolsq(cell) < reso.invresolsq_limit() ) << "\n";
	  if ( hklsam.in_resolution(hkl) && ( abs(ih) > lim.h() || abs(ik) > lim.k() || abs(il) > lim.l() ) ) std::cout << "Err " << lim.format() << " " << hkl.format() << "\n";
	}
  }

  /*
  // Compare reciprocal cell conventions:
  clipper::Cell cell( clipper::Cell_descr( 20.0, 24.0, 18.0, 60.0, 70.0, 80.0 ) );
  
  // Orthogonalisation matrix (A**-1)
  // Reciprocal cell
  double as = cell.a_star();
  double bs = cell.b_star();
  double cs = cell.c_star();
  double cas = cos(cell.alpha_star());
  double cbs = cos(cell.beta_star() );
  double cgs = cos(cell.gamma_star());
  double sas = sin(cell.alpha_star());
  double sbs = sin(cell.beta_star() );
  double sgs = sin(cell.gamma_star());

  double cc = cell.c();
  double ca = cos(cell.alpha());

  Mat33<> bmat = Mat33<>(as, bs*cgs, cs*cbs,
			 0., bs*sgs, -cs*sbs*ca,
			 0., 0., 1.0/cc);
  std::cout << cell.matrix_frac().format() << "\n";
  std::cout << bmat.format() << "\n";
    for ( int ih = 0; ih <= 5; ih++ )
      for ( int ik = 0; ik <= 5; ik++ )
	for ( int il = 0; il <= 5; il++ ) {
	  HKL hkl(ih,ik,il);
	  std::cout << hkl.format() << "\n";
	  Vec3<> v1 = cell.matrix_frac().transpose()*hkl.coord_reci_frac();
	  Vec3<> v2 = bmat*hkl.coord_reci_frac();
	  std::cout << v1.format() << "\t" << v1*v1 << "\n";
	  std::cout << v2.format() << "\t" << v2*v2 << "\n";
	  std::cout << hkl.invresolsq(cell) << "\n";
	}
  */
}
