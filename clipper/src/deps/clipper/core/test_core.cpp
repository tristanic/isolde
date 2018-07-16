/* test_core.cpp: implementation file for clipper core self-test */
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


#include "test_core.h"
#include "test_data.h"
#include "hkl_compute.h"
#include "hkl_operators.h"
#include "resol_targetfn.h"
#include "atomsf.h"
#include "xmap.h"
#include "rotation.h"

#include <algorithm>


namespace clipper {


bool Test_core::operator() () {
  typedef HKL_info::HKL_reference_index HRI;
  data::Test_data data;
  const HKL_data<data32::F_sigF>& fsig = data.hkl_data_f_sigf();
  const HKL_data<data32::ABCD>&   abcd = data.hkl_data_abcd();
  const Spacegroup spgr = fsig.hkl_info().spacegroup();
  const Cell       cell = fsig.hkl_info().cell();
  const Resolution reso = fsig.hkl_info().resolution();

  // test data sizes
  {
    test( "FTYPE32", sizeof(clipper::ftype32), 4 );
    test( "ITYPE32", sizeof(clipper::itype32), 4 );
    test( "FTYPE64", sizeof(clipper::ftype64), 8 );
    test( "ITYPE64", sizeof(clipper::itype64), 8 );
  }

  // test NaN functions
  {
    ftype64 nandd = Util::nand();
    ftype64 nandf = Util::nanf();
    ftype32 nanfd = Util::nand();
    ftype32 nanff = Util::nanf();
    test( "ISNANDD", Util::isnan( nandd )?1:0, 1 );
    test( "ISNANDF", Util::isnan( nandf )?1:0, 1 );
    test( "ISNANFD", Util::isnan( nanfd )?1:0, 1 );
    test( "ISNANFF", Util::isnan( nanff )?1:0, 1 );
    test( "IS_NANDD", Util::is_nan( nandd )?1:0, 1 );
    test( "IS_NANDF", Util::is_nan( nandf )?1:0, 1 );
    test( "IS_NANFD", Util::is_nan( nanfd )?1:0, 1 );
    test( "IS_NANFF", Util::is_nan( nanff )?1:0, 1 );
    for ( double x = -30.0; x < 30.5; x += 1.0 ) {
      test( "ISNANX",  Util::isnan(exp(x))?1:0, 0 );
      test( "IS_NANX",  Util::is_nan(exp(x))?1:0, 0 );
    }
  }

  // test Sim functions for consistency
  for ( double x=0.1; x<50.0; x*=1.5 ) {
    test( "SIMINV", x, Util::invsim(Util::sim(x)), 0.001*x );
    test( "SIMNEG", Util::sim(x), -Util::sim(-x), 0.001 );
    test( "SIMINT", Util::sim(x), (Util::sim_integ(x+0.0001)-Util::sim_integ(x))/0.0001, 0.001 );
    test( "SIMDRV", Util::sim_deriv(x), (Util::sim(x+0.0001)-Util::sim(x))/0.0001, 0.001 );
    if ( x < 2.0 ) test( "SIM_I0", Util::sim(x), (Util::bessel_i0(x+0.0001)/Util::bessel_i0(x)-1.0)/0.0001, 0.001 );
  }

  // test HL<->abcd conversion
  for ( HRI ih = abcd.first(); !ih.last(); ih.next() ) {
    data32::Compute_phifom_from_abcd cap;
    data32::Compute_abcd_from_phifom cpa;
    data32::Phi_fom phiw1, phiw2;
    data32::ABCD    abcd1, abcd2;
    LogPhaseProb<24> q(ih.hkl_class());
    q.set_abcd( abcd[ih] ); // 1st route
    q.get_phi_fom( phiw1 );
    abcd1 = cpa( ih, phiw1 );
    phiw1 = cap( ih, abcd1 );
    phiw2 = cap( ih, abcd[ih] ); // 2nd route
    q.set_phi_fom( phiw2 );
    q.get_abcd( abcd2 );
    phiw2 = cap( ih, abcd2 );
    test( "HLW", phiw1.fom(), phiw2.fom(), 0.003 );
  }

  // test map calculation
  {
    HKL_data<data32::Phi_fom> pw( fsig.hkl_info() );
    HKL_data<data32::F_phi> fp1( fsig.hkl_info() ), fp2( fsig.hkl_info() );
    pw.compute( abcd, data32::Compute_phifom_from_abcd() );
    fp1.compute( fsig, pw, data32::Compute_fphi_from_fsigf_phifom() );
    Grid_sampling grid( spgr, cell, reso );
    Xmap<float> xmap( spgr, cell, grid );
    xmap.fft_from( fp1 );
    xmap.fft_to( fp2 );
    for ( HRI ih = fp1.first(); !ih.last(); ih.next() ) {
      std::complex<float> ab1(0.0,0.0), ab2(0.0,0.0);
      if ( !fp1[ih].missing() ) ab1 = fp1[ih];
      if ( !fp2[ih].missing() ) ab2 = fp2[ih];
      test( "FFT-A", ab1.real(), ab2.real(), 0.01 );
      test( "FFT-B", ab1.imag(), ab2.imag(), 0.01 );
    }
  }

  // test resolution ordinals
  {
    std::vector<clipper::ftype> resols;
    for ( HRI ih = fsig.first(); !ih.last(); ih.next() )
      resols.push_back( ih.invresolsq() );
    clipper::Generic_ordinal ordinal;
    ordinal.init(resols);
    clipper::Generic_ordinal ordinv = ordinal;
    ordinv.invert();
    clipper::Resolution_ordinal resord;
    resord.init( fsig.hkl_info(), 1.0 );
    std::sort( resols.begin(), resols.end() );
    for ( int i = 0; i < resols.size(); i++ )
      test( "RESORD", resols[i], ordinv.ordinal(resord.ordinal(resols[i])), 0.001 );
  }

  // test resolution functions
  {
    std::vector<ftype> param( 10, 1.0 );
    BasisFn_spline basisfn( fsig, 10 );
    TargetFn_meanFnth<data32::F_sigF> targetfn( fsig, 2.0 );
    ResolutionFn rfn( fsig.hkl_info(), basisfn, targetfn, param );
    test( "RESFN0", rfn.params()[0], 229690.9746, 0.1 );
    test( "RESFN1", rfn.params()[1], 216481.7609, 0.1 );
    test( "RESFN2", rfn.params()[2],  78484.9498, 0.1 );
    test( "RESFN3", rfn.params()[3], 148774.2654, 0.1 );
    test( "RESFN4", rfn.params()[4],  69255.6000, 0.1 );
    test( "RESFN5", rfn.params()[5], 143032.5088, 0.1 );
    test( "RESFN6", rfn.params()[6], 110371.3524, 0.1 );
    test( "RESFN7", rfn.params()[7], 108711.3487, 0.1 );
    test( "RESFN8", rfn.params()[8], 150487.5496, 0.1 );
    test( "RESFN9", rfn.params()[9], 141713.7420, 0.1 );
  }

  // test atom shape function derivatives
  {
    ftype d = 0.001;
    Coord_orth co( 1.0, 2.0, 3.0 );
    AtomShapeFn sf( co, "N", 0.25, 1.0 );
    AtomShapeFn sfx( co+Coord_orth(d,0.0,0.0), "N", 0.25, 1.0 );
    AtomShapeFn sfy( co+Coord_orth(0.0,d,0.0), "N", 0.25, 1.0 );
    AtomShapeFn sfz( co+Coord_orth(0.0,0.0,d), "N", 0.25, 1.0 );
    AtomShapeFn sfo( co, "N", 0.25, 1.0+d );
    AtomShapeFn sfu( co, "N", 0.251, 1.0 );
    AtomShapeFn sfx2( co+Coord_orth(2.0*d,0.0,0.0), "N", 0.25, 1.0 );
    AtomShapeFn sfy2( co+Coord_orth(0.0,2.0*d,0.0), "N", 0.25, 1.0 );
    AtomShapeFn sfz2( co+Coord_orth(0.0,0.0,2.0*d), "N", 0.25, 1.0 );
    AtomShapeFn sfo2( co, "N", 0.25, 1.0+d+d );
    AtomShapeFn sfu2( co, "N", 0.252, 1.0 );
    AtomShapeFn sfxy( co+Coord_orth(d,d,0.0), "N", 0.25, 1.0 );
    AtomShapeFn sfyz( co+Coord_orth(0.0,d,d), "N", 0.25, 1.0 );
    AtomShapeFn sfzx( co+Coord_orth(d,0.0,d), "N", 0.25, 1.0 );
    U_aniso_orth uai( 0.25, 0.25, 0.25, 0.0, 0.0, 0.0 );
    AtomShapeFn sfuai( co, "N", uai, 1.0 );
    std::vector<AtomShapeFn::TYPE> params;
    params.push_back( AtomShapeFn::X );
    params.push_back( AtomShapeFn::Y );
    params.push_back( AtomShapeFn::Z );
    params.push_back( AtomShapeFn::Occ );
    params.push_back( AtomShapeFn::Uiso );
    params.push_back( AtomShapeFn::U11 );
    params.push_back( AtomShapeFn::U22 );
    params.push_back( AtomShapeFn::U33 );
    params.push_back( AtomShapeFn::U12 );
    params.push_back( AtomShapeFn::U13 );
    params.push_back( AtomShapeFn::U23 );
    sf.agarwal_params() = params;
    for ( int i = 0; i < 100; i++ ) {
      Coord_orth c2 = co + Coord_orth( 0.11*(i%5-1.9), 0.13*(i%7-2.8), 0.15*(i%9-3.7) );
      double f;
      std::vector<double> g(11);
      Matrix<double> c(11,11);
      test( "ATOMSF-A", sfuai.rho(c2), sf.rho(c2), 1.0e-8 );
      sf.rho_curv( c2, f, g, c );
      test( "ATOMSF-G", (sfx.rho(c2)-sf.rho(c2))/d, g[0], 0.01 );
      test( "ATOMSF-G", (sfy.rho(c2)-sf.rho(c2))/d, g[1], 0.01 );
      test( "ATOMSF-G", (sfz.rho(c2)-sf.rho(c2))/d, g[2], 0.01 );
      test( "ATOMSF-G", (sfo.rho(c2)-sf.rho(c2))/d, g[3], 0.01 );
      test( "ATOMSF-G", (sfu.rho(c2)-sf.rho(c2))/d, g[4], 0.05 );
      test( "ATOMSF-C", (sfx2.rho(c2)-2*sfx.rho(c2)+sf.rho(c2))/(d*d), c(0,0), 0.1 );
      test( "ATOMSF-C", (sfy2.rho(c2)-2*sfy.rho(c2)+sf.rho(c2))/(d*d), c(1,1), 0.1 );
      test( "ATOMSF-C", (sfz2.rho(c2)-2*sfz.rho(c2)+sf.rho(c2))/(d*d), c(2,2), 0.1 );
      test( "ATOMSF-C", (sfxy.rho(c2)-sfx.rho(c2)-sfy.rho(c2)+sf.rho(c2))/(d*d), c(0,1), 0.1 );
      test( "ATOMSF-C", (sfyz.rho(c2)-sfy.rho(c2)-sfz.rho(c2)+sf.rho(c2))/(d*d), c(1,2), 0.1 );
      test( "ATOMSF-C", (sfzx.rho(c2)-sfz.rho(c2)-sfx.rho(c2)+sf.rho(c2))/(d*d), c(2,0), 0.1 );
    }
    for ( int j = 0; j < 20; j++ ) {
      ftype x = 0.19*(j%5-1.9);
      ftype y = 0.15*(j%7-2.8);
      ftype z = 0.13*(j%9-3.7);
      U_aniso_orth ua  ( x*x+0.2, y*y+0.2, z*z+0.2, y*z, z*x, x*y );
      U_aniso_orth ua00( ua.mat00()+d, ua.mat11(), ua.mat22(),
			 ua.mat01(), ua.mat02(), ua.mat12() );
      U_aniso_orth ua11( ua.mat00(), ua.mat11()+d, ua.mat22(),
			 ua.mat01(), ua.mat02(), ua.mat12() );
      U_aniso_orth ua22( ua.mat00(), ua.mat11(), ua.mat22()+d,
			 ua.mat01(), ua.mat02(), ua.mat12() );
      U_aniso_orth ua01( ua.mat00(), ua.mat11(), ua.mat22(),
			 ua.mat01()+d, ua.mat02(), ua.mat12() );
      U_aniso_orth ua02( ua.mat00(), ua.mat11(), ua.mat22(),
			 ua.mat01(), ua.mat02()+d, ua.mat12() );
      U_aniso_orth ua12( ua.mat00(), ua.mat11(), ua.mat22(),
			 ua.mat01(), ua.mat02(), ua.mat12()+d );
      AtomShapeFn sfua ( co, "N", ua, 1.0 );
      AtomShapeFn sfuax( co+Coord_orth(d,0.0,0.0), "N", ua, 1.0 );
      AtomShapeFn sfuay( co+Coord_orth(0.0,d,0.0), "N", ua, 1.0 );
      AtomShapeFn sfuaz( co+Coord_orth(0.0,0.0,d), "N", ua, 1.0 );
      AtomShapeFn sfuao( co, "N", ua, 1.0+d );
      AtomShapeFn sfua00( co, "N", ua00, 1.0 );
      AtomShapeFn sfua11( co, "N", ua11, 1.0 );
      AtomShapeFn sfua22( co, "N", ua22, 1.0 );
      AtomShapeFn sfua01( co, "N", ua01, 1.0 );
      AtomShapeFn sfua02( co, "N", ua02, 1.0 );
      AtomShapeFn sfua12( co, "N", ua12, 1.0 );
      sfua.agarwal_params() = params;
      for ( int i = 0; i < 50; i++ ) {
	Coord_orth c2 = co + Coord_orth( 0.11*(i%5-1.9), 0.13*(i%7-2.8), 0.15*(i%9-3.7) );
	double f;
	std::vector<double> g(11);
	Matrix<double> c(11,11);
	sfua.rho_grad( c2, f, g );
	test( "ATOMSF-AG", ( sfuax.rho(c2)-sfua.rho(c2))/d,  g[0], 0.01 );
	test( "ATOMSF-AG", ( sfuay.rho(c2)-sfua.rho(c2))/d,  g[1], 0.01 );
	test( "ATOMSF-AG", ( sfuaz.rho(c2)-sfua.rho(c2))/d,  g[2], 0.01 );
	test( "ATOMSF-AG", ( sfuao.rho(c2)-sfua.rho(c2))/d,  g[3], 0.01 );
	test( "ATOMSF-AG", (sfua00.rho(c2)-sfua.rho(c2))/d,  g[5], 0.05 );
	test( "ATOMSF-AG", (sfua11.rho(c2)-sfua.rho(c2))/d,  g[6], 0.05 );
	test( "ATOMSF-AG", (sfua22.rho(c2)-sfua.rho(c2))/d,  g[7], 0.05 );
	test( "ATOMSF-AG", (sfua01.rho(c2)-sfua.rho(c2))/d,  g[8], 0.05 );
	test( "ATOMSF-AG", (sfua02.rho(c2)-sfua.rho(c2))/d,  g[9], 0.05 );
	test( "ATOMSF-AG", (sfua12.rho(c2)-sfua.rho(c2))/d, g[10], 0.05 );
      }
    }
  }

  // test spacegroups
  {
    const char *pgs[] = {"-P 1", "-P 2", "-P 2y", "-P 2x", "-P 2\"", "-P 2y\"", "-P 2x\"", "-P 2'", "-P2y'", "-P 2x'", "-P 2 2", "-P 2 2\"", "-P 2 2\"(y,z,x)", "-P 2 2\"(z,x,y)", "-P3", "-P 3 (y,z,x)", "-P 3 (z,x,y)", "-P 3 (-x,y,z)", "-P 3 (y,z,-x)", "-P 3 (z,-x,y)", "-P 3*", "-P 3* (-x,y,z)", "-P 3* (x,-y,z)", "-P 3* (x,y,-z)", "-P 3 2", "-P 3 2 (y,z,x)", "-P 3 2 (z,x,y)", "-P 3* 2", "-P 3* 2 (-x,y,z)", "-P 3* 2 (x,-y,z)", "-P 3* 2 (-x,-y,z)", "-P 3 2\"", "-P 3 2\"(z,x,y)", "-P 3 2\"(y,z,x)", "-P 3 2\"(-x,y,z)", "-P 3 2\"(z,-x,y)", "-P 3 2\"(y,z,-x)", "-P 4", "-P 4 (y,z,x)", "-P 4 (z,x,y)", "-P 4 2", "-P 4 2 (y,z,x)", "-P 4 2 (z,x,y)", "-P 6", "-P 6 (y,z,x)", "-P 6 (z,x,y)", "-P 6 2", "-P 6 2 (y,z,x)", "-P 6 2 (z,x,y)", "-P 2 2 3", "-P 4 2 3" };
    std::vector<String> hallsymbols;
    for ( int i = 0; i < data::sgdata_size; i++ )
      hallsymbols.push_back( data::sgdata[i].hall );
    for ( int i = 0; i < sizeof(pgs)/sizeof(pgs[0]); i++ )
      hallsymbols.push_back( pgs[i] );
    Cell cellc( Cell_descr( 37, 37, 37 ) );
    Cell cellha( Cell_descr( 37, 37, 37, 120, 90, 90 ) );
    Cell cellhb( Cell_descr( 37, 37, 37, 90, 120, 90 ) );
    Cell cellhc( Cell_descr( 37, 37, 37, 90, 90, 120 ) );
    Cell cellha1( Cell_descr( 37, 37, 37, 60, 90, 90 ) );
    Cell cellhb1( Cell_descr( 37, 37, 37, 90, 60, 90 ) );
    Cell cellhc1( Cell_descr( 37, 37, 37, 90, 90, 60 ) );
    Cell cell;
    String symbol;
    Spacegroup sg;
    for ( int s = 0; s < hallsymbols.size(); s++ ) {
      try {
	symbol = hallsymbols[s];
	sg = Spacegroup( Spgr_descr( symbol, Spgr_descr::Hall ) );
	// identify trigonal/hexagonal groups
	cell = cellc;
	for ( int sym = 1; sym < sg.num_symops(); sym++ ) {
	  if ( ( sg.symop(sym).rot()(1,1) * sg.symop(sym).rot()(1,2) == -1 ) ||
	       ( sg.symop(sym).rot()(2,1) * sg.symop(sym).rot()(2,2) == -1 ) )
	    cell = cellha;
	  if ( ( sg.symop(sym).rot()(0,0) * sg.symop(sym).rot()(0,2) == -1 ) ||
	       ( sg.symop(sym).rot()(2,0) * sg.symop(sym).rot()(2,2) == -1 ) )
	    cell = cellhb;
	  if ( ( sg.symop(sym).rot()(0,0) * sg.symop(sym).rot()(0,1) == -1 ) ||
	       ( sg.symop(sym).rot()(1,0) * sg.symop(sym).rot()(1,1) == -1 ) )
	    cell = cellhc;
	  if ( ( sg.symop(sym).rot()(1,1) * sg.symop(sym).rot()(1,2) == 1 ) ||
	       ( sg.symop(sym).rot()(2,1) * sg.symop(sym).rot()(2,2) == 1 ) )
	    cell = cellha1;
	  if ( ( sg.symop(sym).rot()(0,0) * sg.symop(sym).rot()(0,2) == 1 ) ||
	       ( sg.symop(sym).rot()(2,0) * sg.symop(sym).rot()(2,2) == 1 ) )
	    cell = cellhb1;
	  if ( ( sg.symop(sym).rot()(0,0) * sg.symop(sym).rot()(0,1) == 1 ) ||
	       ( sg.symop(sym).rot()(1,0) * sg.symop(sym).rot()(1,1) == 1 ) )
	    cell = cellhc1;
	}
	for ( int i = 0; i < 100; i++ ) {
	  HKL rfl( i%5-2, i%7-3, i%9-4 );
	  ftype s0 = rfl.invresolsq(cell);
	  for ( int sym = 1; sym < sg.num_symops(); sym++ ) {
	    ftype s1 = rfl.transform(sg.symop(sym)).invresolsq(cell);
	    test( "SG "+symbol, s0, s1, 1.0e-12 );
	  }
	}
	Grid_sampling grid( sg, cell, Resolution( 4.0 ) );
	Xmap<ftype32> xmap( sg, cell, grid );
      } catch ( Message_base ) {
	test( "SG "+symbol, sg.spacegroup_number(), -1 );
      }
    }
  }

  // test rotations
  {
    for ( ftype x = -1.0; x < 1.0; x += 0.02 )
      for ( ftype y = -1.0; y < 1.0; y += 0.02 )
	for ( ftype z = -1.0; z < 1.0; z += 0.02 ) {
	  ftype s = x*x + y*y + z*z;
	  if ( s < 1.0 ) {
	    ftype w = sqrt( 1.0 - s );
	    Rotation rot( w, x, y, z );
	    Rotation rotinv = rot.inverse();
	    Rotation r1( rot.matrix() );
	    Rotation r2( rot.polar_ccp4() );
	    test( "ROT/MAT "+rot.format(), (rotinv*r1).abs_angle(), 0.0, 1.0e-6 );
	    test( "ROT/POL "+rot.format(), (rotinv*r2).abs_angle(), 0.0, 1.0e-6 );
	  }
	}
    /*
    Mat33<> mat1( -0.18332,  0.02511, -0.98273,
		  0.02184, -0.99932, -0.02960,
		 -0.98281, -0.02689,  0.18265 );
    Rotation r1( mat1 );
    Mat33<> mat2( r1.matrix() );
    Rotation r2( mat2 );
    std::cout << mat1.format() << std::endl << r1.format() << std::endl << mat2.format() << std::endl << r2.format() << std::endl << std::endl;
    Rotation r1a( mat1 ); r1a.norm();
    Mat33<> mat2a( r1a.matrix() );
    Rotation r2a( mat2a );
    Mat33<> mat3a( r2a.matrix() );
    Rotation r3a( mat3a );
    std::cout << mat1.format() << mat1.det() << std::endl << r1a.format() << std::endl << mat2a.format() << mat2a.det() << std::endl << r2a.format() << std::endl << mat3a.format() << mat3a.det() << std::endl << r3a.format() << std::endl << std::endl;
    Rotation r3(0.02,0.42,0.02,-0.50);
    std::cout << r3.format() << std::endl << Rotation(r3.matrix()).format() << std::endl;
    */
  }

  return ( error_count == 0 );
}


} // namespace clipper
