/* test_contrib.cpp: implementation file for clipper contrib self-test */
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


#include "test_contrib.h"
#include "mapfilter.h"
#include "convolution_search.h"
#include "fffear.h"
#include "sfcalc.h"
#include "sfcalc_obs.h"
#include "sfweight.h"
#include "sfscale.h"
#include "../core/test_data.h"
#include "../core/hkl_compute.h"


namespace clipper {

namespace data {

  float contrib_vals[] = 
    {62085,0,-2.35103e-05,-537.853,2.04953e-05,271.468,313.756,-1.14655e-05,88.939,-1.97969e-06,-1.41777e-05,137.227,-1020.36,-230.467,258.784,19.4909,-21.0352,168.158,82.8908,82.9012,-3.63488e-05,-106.365,-4.12642,311.855,-182.047,291.716,-473.216,24.3527,58.3351,-90.6147,-72.7842,46.5314,123.021,-21.1225,-156.35,110.062,-277.139,-442.304,-142.18,89.8283,1.87174e-05,-73.6002,271.811,-64.1525,-270.657,-265.162,273.196,779.432,256.913,276.635,226.079,544.524,644.946,-244.673,8.01921e-06,-504.567,4.02743,-451.095,27.369,95.8539,-1049.43,0,332.786,-326.619,-120.613,-233.546,374.933,-125.7,-1.43747e-05,328.855,67.4085,-117.594,-579.966,-250.85,-286.214,421.978,-407.023,-3.5e-12,379.259,20.3695,-546.954,33.0029,638.675,-649.005,70.5812,-146.142,224.169,43.466,302.437,49.827,390.892,-184.531,215.482,-101.874,545.265,474.695,-3.76247e-05,364.171,-5.44969e-05,-1246.74,-2.18374e-05,-499.581,0.99955,5.34325e-06,70.7734,-3.57623e-06,-81.8145,0.965747,629.401,-2.29999e-05,196.403,-7.17707e-06,0.999995,56.9272,-1.26714e-06,-5.00108,5.96175e-07,0.47382,-1.68563e-05,163.153,-2.71878e-06,26.3152,0.803509,-636.356,-143.732,169.765,38.3443,0.981854,110.28,8.30599,-57.1612,-4.30521,0.880984,6.49867,-51.951,11.3566,-90.7857,0.504998,187.239,187.262,61.0691,61.0767,0.844098,-3.9956e-05,-116.921,-1.60269e-05,-46.8987,0.19955,-6.80942,514.623,-1.80122,136.128,0.937416,-158.815,254.487,-13.0074,20.8433,0.826057,-459.028,23.6226,-46.6823,2.40237,0.964024,64.2033,-99.7301,9.28487,-14.4226,0.621669,-25.9256,4.45139,-35.6622,6.12313,0.102692,34.3516,-24.1816,77.2152,-54.3552,0.601134,-203.365,-324.564,11.8764,18.9543,0.961987,-104.664,66.1263,3.25264,-2.055,0.706816,9.51035e-05,-373.964,3.97517e-05,-156.311,0.736061,186.115,-43.9266,-23.964,5.65595,0.917733,-149.75,-146.71,34.3991,33.7007,0.918732,216.55,617.822,-2.34832,-6.6998,0.983464,103.747,111.711,-47.8103,-51.4805,0.933024,138.133,332.702,-27.9613,-67.3465,0.978446,600.51,-227.816,37.9841,-14.4101,0.978923,5.68373e-06,-357.619,-7.21378e-07,16.5032,0.999258,3.50463,-392.538,0.0360139,-4.03376,0.967671,12.9659,45.4102,-3.82199,-13.3857,0.486054,-687.716,5e-12,71.8843,4.5e-12,1,231.469,-227.18,-9.31488,9.14227,0.965045,-76.4861,-148.101,8.49971,16.4581,0.864654,333.513,-111.814,11.3842,-3.81669,0.96192,-3.95762e-06,90.5398,6.19183e-06,-82.0133,0.961299,-29.0476,50.6733,-40.9055,71.3594,0.333013,-455.874,-197.177,-2.65068,-1.14648,0.974159,-180.118,265.555,21.9202,-32.3179,0.979928,-5.86636,-3.24e-12,152.988,1.5e-12,0.987289,81.058,4.35351,-100.564,-5.40116,0.932738,-457.759,27.6209,5.42575,-0.327387,0.976693,586.908,-596.401,17.798,-18.0858,0.990935,199.589,-413.259,72.816,-150.769,0.873015,530.378,102.839,175.409,34.0114,0.946546,235.809,38.8499,-10.2464,-1.68811,0.939609,323.264,-152.605,12.3157,-5.81392,0.957015,165.082,-78.0461,-6.46434,3.05617,0.898444,298.646,259.994,-75.8603,-66.0422,0.98677,-3.71094e-06,35.9183,9.71826e-06,-128.722,0.997217};
  float contrib_tols[] = 
    {3.1,1e-09,0.00038,0.027,9.7e-05,0.014,0.016,6.1e-05,0.0045,8.5e-05,0.00012,0.007,0.051,0.012,0.013,0.001,0.0011,0.0085,0.0042,0.0042,0.00025,0.0058,0.00032,0.016,0.0093,0.015,0.024,0.0016,0.0031,0.0046,0.0038,0.0024,0.0063,0.0012,0.008,0.0056,0.014,0.022,0.0073,0.0047,0.00012,0.0038,0.014,0.0033,0.014,0.013,0.014,0.039,0.013,0.014,0.012,0.027,0.032,0.012,9e-05,0.025,0.00036,0.023,0.0014,0.0049,0.053,0.00048,0.017,0.016,0.0062,0.012,0.019,0.0064,0.00012,0.017,0.0035,0.006,0.029,0.013,0.014,0.021,0.021,0.00018,0.019,0.0012,0.028,0.002,0.032,0.033,0.0036,0.0074,0.011,0.0022,0.015,0.0026,0.02,0.0093,0.011,0.0052,0.027,0.024,0.0002,0.018,0.00089,0.063,0.00036,0.025,5e-05,2.5e-05,0.0036,2.9e-05,0.0042,4.8e-05,0.032,0.00012,0.0099,3.8e-05,5e-05,0.0029,5.4e-05,0.00027,4.1e-06,2.4e-05,0.00015,0.0082,2.4e-05,0.0014,4.1e-05,0.032,0.0076,0.0088,0.0021,4.9e-05,0.0056,0.00044,0.0029,0.00024,4.4e-05,0.00036,0.0027,0.00063,0.0046,2.6e-05,0.0095,0.0096,0.0031,0.0031,4.2e-05,0.00027,0.0063,0.00011,0.0025,1.1e-05,0.00052,0.026,0.00014,0.0069,4.7e-05,0.0082,0.013,0.0007,0.0011,4.1e-05,0.023,0.0015,0.0025,0.00015,4.8e-05,0.0033,0.0051,0.00048,0.00075,3.1e-05,0.0013,0.00026,0.0018,0.00036,5.2e-06,0.0018,0.0012,0.004,0.0028,3e-05,0.01,0.016,0.00073,0.0012,4.8e-05,0.0053,0.0035,0.00019,0.00013,3.6e-05,0.0006,0.019,0.00025,0.008,3.8e-05,0.0094,0.0023,0.0013,0.00031,4.6e-05,0.0077,0.0075,0.0018,0.0018,4.6e-05,0.011,0.031,0.00015,0.00038,4.9e-05,0.0053,0.0058,0.0026,0.0028,4.7e-05,0.007,0.017,0.0015,0.0036,4.9e-05,0.03,0.011,0.0022,0.00082,4.9e-05,6.4e-05,0.018,5.9e-06,0.001,5e-05,0.00032,0.02,4.5e-06,0.00035,4.8e-05,0.00067,0.0023,0.0002,0.00069,2.5e-05,0.035,0.00031,0.0038,3.3e-05,5e-05,0.012,0.012,0.00056,0.00053,4.8e-05,0.0039,0.0075,0.00048,0.00093,4.3e-05,0.017,0.0057,0.00067,0.00024,4.8e-05,3.2e-05,0.0045,2.9e-05,0.0041,4.8e-05,0.0015,0.0026,0.0021,0.0037,1.7e-05,0.023,0.01,0.00021,0.00012,4.9e-05,0.0093,0.013,0.0012,0.0017,4.9e-05,0.00042,2.7e-06,0.0078,6.9e-05,4.9e-05,0.0042,0.00025,0.0052,0.00032,4.7e-05,0.023,0.0017,0.00036,2.7e-05,4.9e-05,0.03,0.03,0.001,0.0011,5e-05,0.01,0.021,0.0037,0.0076,4.4e-05,0.027,0.0053,0.0089,0.0017,4.7e-05,0.012,0.002,0.00066,0.00011,4.7e-05,0.016,0.0077,0.00065,0.00031,4.8e-05,0.0084,0.004,0.00042,0.00019,4.5e-05,0.015,0.013,0.004,0.0034,4.9e-05,1.9e-05,0.0019,4.6e-05,0.0065,5e-05};
  int contrib_size = sizeof(contrib_vals)/sizeof(contrib_vals[0]);

}

bool Test_contrib::operator() () {
  data_val = std::vector<float>( data::contrib_vals, 
				 data::contrib_vals+data::contrib_size );
  data_tol = std::vector<float>( data::contrib_tols, 
				 data::contrib_tols+data::contrib_size );

  typedef HKL_info::HKL_reference_index HRI;
  data::Test_data data;
  const HKL_data<data32::F_sigF>& fsig = data.hkl_data_f_sigf();
  const HKL_data<data32::ABCD>&   abcd = data.hkl_data_abcd();
  const clipper::Atom_list&       xyzb = data.atom_list();
  const Spacegroup spgr = fsig.hkl_info().spacegroup();
  const Cell       cell = fsig.hkl_info().cell();
  const Resolution reso = fsig.hkl_info().resolution();

  typedef HKL_info::HKL_reference_index HRI;
  typedef Xmap<float>::Map_reference_index MRI;

  // test sfcalc objects
  {
    // select spacegroup
    std::vector<String> hallsymbols;
    for ( int i = 0; i < data::sgdata_size; i += 10 )
      hallsymbols.push_back( data::sgdata[i].hall );
    // build model
    Atom_list atoms;
    Atom atom = Atom::null();
    atom.set_occupancy( 1.0 );
    atom.set_u_iso( 0.5 );
    atom.set_element( "C" );
    atom.set_coord_orth( Coord_orth( 12, 8, 5 ) );
    atoms.push_back( atom );
    atom.set_element( "N" );
    atom.set_coord_orth( Coord_orth( 11, 6, 4 ) );
    atoms.push_back( atom );
    atom.set_element( "O" );
    atom.set_coord_orth( Coord_orth( 13, 5, 4 ) );
    atoms.push_back( atom );
    // calc cell
    Cell cellc( Cell_descr( 37, 37, 37 ) );
    Cell cellha( Cell_descr( 37, 37, 37, 120, 90, 90 ) );
    Cell cellhb( Cell_descr( 37, 37, 37, 90, 120, 90 ) );
    Cell cellhc( Cell_descr( 37, 37, 37, 90, 90, 120 ) );
    Cell cg;
    String symbol;
    Spacegroup sg;
    for ( int s = 0; s < hallsymbols.size(); s++ ) {
      try {
	symbol = hallsymbols[s];
	sg = Spacegroup( Spgr_descr( symbol, Spgr_descr::Hall ) );
	// identify trigonal/hexagonal groups
	cg = cellc;
	for ( int sym = 1; sym < sg.num_symops(); sym++ ) {
	  if ( ( sg.symop(sym).rot()(1,1) * sg.symop(sym).rot()(1,2) == -1 ) ||
	       ( sg.symop(sym).rot()(2,1) * sg.symop(sym).rot()(2,2) == -1 ) )
	    cg = cellha;
	  if ( ( sg.symop(sym).rot()(0,0) * sg.symop(sym).rot()(0,2) == -1 ) ||
	       ( sg.symop(sym).rot()(2,0) * sg.symop(sym).rot()(2,2) == -1 ) )
	    cg = cellhb;
	  if ( ( sg.symop(sym).rot()(0,0) * sg.symop(sym).rot()(0,1) == -1 ) ||
	       ( sg.symop(sym).rot()(1,0) * sg.symop(sym).rot()(1,1) == -1 ) )
	    cg = cellhc;
	}
	HKL_info hkl_info( sg, cg, Resolution( 5.0 ), true );
	HKL_data<data32::F_phi> fp1( hkl_info );
	HKL_data<data32::F_phi> fp2( hkl_info );
	SFcalc_iso_sum<float>( fp1, atoms );
	SFcalc_iso_fft<float>( fp2, atoms, 2.5, 2.5, 0.25 );
	// begin extra fft tests
	Grid_sampling gg( sg, cg, hkl_info.resolution() );
	Xmap<float> xg( sg, cg, gg );
	FFTmap fftmap( sg, cg, gg );
	xg.fft_from( fp2, Xmap_base::Normal );
	xg.fft_to( fp2, Xmap_base::Sparse );
	fftmap.fft_rfl_to_map( fp2, xg );
	xg.fft_to( fp2, Xmap_base::Normal );
	xg.fft_from( fp2, Xmap_base::Sparse );
	fftmap.fft_map_to_rfl( xg, fp2 );
	// end extra fft tests
	float tol = 0.005 * fp1[ HKL( 0, 0, 0 ) ].f();
	for ( HRI ih = fp1.first(); !ih.last(); ih.next() ) {
	  std::complex<float> ab1(0.0,0.0), ab2(0.0,0.0);
	  if ( !fp1[ih].missing() ) ab1 = fp1[ih];
	  if ( !fp2[ih].missing() ) ab2 = fp2[ih];
	  test( "SF-A", ab1.real(), ab2.real(), tol );
	  test( "SF-B", ab1.imag(), ab2.imag(), tol );
	}
      } catch ( Message_base ) {
	test( "SFSG "+symbol, sg.spacegroup_number(), -1 );
      }
    }
  }

  // test sfcalc_obs and sfweight
  {
    // sfcalc_obs
    HKL_data<data32::F_phi> fcal( fsig.hkl_info() );
    SFcalc_obs_bulk<float> sfcb;
    sfcb( fcal, fsig, xyzb );

    // sfcalc_obs results
    for ( HRI ih = fcal.first(); !ih.last(); ih.next() )
      if ( ih.index() % 20 == 0 ) {
	std::complex<float> ab( fcal[ih] );
	test( "SFO-A", ab.real() );
	test( "SFO-B", ab.imag() );
      }

    // sfweight
    //SFcalc_iso_fft<float> sfc;
    //sfc( fcal, xyzb );
    HKL_data<data32::F_phi> fb1( fsig.hkl_info() ), fb2( fsig.hkl_info() );
    HKL_data<data32::F_phi> fd1( fsig.hkl_info() ), fd2( fsig.hkl_info() );
    HKL_data<data32::Phi_fom> phiw1(fsig.hkl_info()), phiw2(fsig.hkl_info());
    HKL_data<data32::Flag> flag( fsig.hkl_info() );
    HKL_data<data32::ABCD> abcd( fsig.hkl_info() ), abcd2( fsig.hkl_info() );
    abcd = data32::ABCD( 0.0, 0.0, 0.0, 0.0 );
    for ( HRI ih = flag.first(); !ih.last(); ih.next() )
      if ( !fsig[ih].missing() )
	flag[ih].flag() = SFweight_spline<float>::BOTH;
      else
	flag[ih].flag() = SFweight_spline<float>::NONE;
    clipper::SFweight_spline<float> sfw( 600, 12 );
    bool fl = sfw( fb1, fd1, phiw1, fsig, fcal, flag );
    if ( !fl ) test( "SFW-FAIL", 1, 0 );

    // sfweight results
    for ( HRI ih = phiw1.first(); !ih.last(); ih.next() )
      if ( !fsig[ih].missing() )
	if ( ih.index() % 20 == 0 ) {
	  std::complex<float> ab_b( fb1[ih] );
	  std::complex<float> ab_d( fd1[ih] );
	  test( "SFWB-A", ab_b.real() );
	  test( "SFWB-B", ab_b.imag() );
	  test( "SFWD-A", ab_d.real() );
	  test( "SFWD-B", ab_d.imag() );
	  test( "SFW-W", phiw1[ih].fom() );
	}

    // sfweight-hl results
    HKL_data<data32::F_sigF> fsig0( fsig.hkl_info() );
    for ( HRI ih = fsig.first(); !ih.last(); ih.next() )
      if ( !fsig[ih].missing() )
	fsig0[ih] = data32::F_sigF( fsig[ih].f(), 0.0 );
    clipper::SFweight_spline<float> sfw1( 600, 12 ), sfw2( 600, 12 );
    bool fl1 = sfw1( fb1, fd1, phiw1, fsig0, fcal, flag );
    if ( !fl1 ) test( "SFW-FAIL1", 1, 0 );
    bool fl2 = sfw2( fb2, fd2, phiw2, abcd2, fsig0, abcd, fcal, flag );
    if ( !fl2 ) test( "SFW-FAIL2", 1, 0 );
    std::vector<ftype> params_e1 = sfw1.params_error();
    std::vector<ftype> params_s1 = sfw1.params_scale();
    std::vector<ftype> params_e2 = sfw2.params_error();
    std::vector<ftype> params_s2 = sfw2.params_scale();
    for ( int i = 0; i < params_e1.size(); i++ ) {
      test( "SFWHL-E", params_e1[i], params_e2[i], 0.025 );
      test( "SFWHL-S", params_s1[i], params_s2[i], 0.025 );
    }
    for ( HRI ih = fsig0.first(); !ih.last(); ih.next() )
      if ( !fsig0[ih].missing() ) {
	clipper::SFweight_spline<float>::TargetResult r00, r01, r10, r11, rhl;
	for ( double s = 0.20; s <= 1.01; s += 0.2 ) {
	  for ( double p = 0.20; p <= 1.01; p += 0.2 ) {
	    double w = p * fsig0[ih].f();
	    w = w * w;
	    const double d = 0.000001;
	    r00 = sfw1.targetfn( ih.hkl_class(), fsig0[ih], fcal[ih], s, w );
	    rhl = sfw1.targethl( ih.hkl_class(), fsig0[ih], abcd[ih], fcal[ih], s, w );
	    //std::cout << fsig0[ih].f() << " " << fcal[ih].f() << " " << s << " " << w << " " << r00.r << " " << rhl.r << std::endl;
	    test( "SFW-TGT-CMP", r00.r, rhl.r, 1.0 );
	    r00 = sfw1.targetfn( ih.hkl_class(), fsig0[ih], fcal[ih], s, w );
	    r01 = sfw1.targetfn( ih.hkl_class(), fsig0[ih], fcal[ih], s, w+d );
	    r10 = sfw1.targetfn( ih.hkl_class(), fsig0[ih], fcal[ih], s+d, w );
	    r11 = sfw1.targetfn( ih.hkl_class(), fsig0[ih], fcal[ih], s+d, w+d );
 	    test( "SFW-FN-DW ", (r01.r-r00.r)/d, r00.dw, 0.02 );
 	    test( "SFW-FN-DS ", (r10.r-r00.r)/d, r00.ds, 0.02 );
 	    test( "SFW-FN-DWW", (r01.dw-r00.dw)/d, r00.dww, 0.02 );
 	    test( "SFW-FN-DSS", (r10.ds-r00.ds)/d, r00.dss, 0.02 );
	    r00 = sfw1.targethl( ih.hkl_class(), fsig0[ih], abcd[ih], fcal[ih], s, w );
	    r01 = sfw1.targethl( ih.hkl_class(), fsig0[ih], abcd[ih], fcal[ih], s, w+d );
	    r10 = sfw1.targethl( ih.hkl_class(), fsig0[ih], abcd[ih], fcal[ih], s+d, w );
	    r11 = sfw1.targethl( ih.hkl_class(), fsig0[ih], abcd[ih], fcal[ih], s+d, w+d );
 	    test( "SFW-HL-DW ", (r01.r-r00.r)/d, r00.dw, 0.02 );
 	    test( "SFW-HL-DS ", (r10.r-r00.r)/d, r00.ds, 0.02 );
 	    test( "SFW-HL-DWW", (r01.dw-r00.dw)/d, r00.dww, 0.02 );
 	    test( "SFW-HL-DSS", (r10.ds-r00.ds)/d, r00.dss, 0.02 );
	  }
	}
      }
  }

  // test map filter objects
  {
    HKL_data<data32::Phi_fom> pw( fsig.hkl_info() );
    HKL_data<data32::F_phi> fp( fsig.hkl_info() );
    pw.compute( abcd, data32::Compute_phifom_from_abcd() );
    fp.compute( fsig, pw, data32::Compute_fphi_from_fsigf_phifom() );
    Grid_sampling grid( spgr, cell, reso, 2.5 );
    Xmap<float> xmap( spgr, cell, grid );
    xmap.fft_from( fp );

    MapFilterFn_step step( 2.5 );
    MapFilter_slow<float> fltr1( step, 1.0, MapFilter_slow<float>::Relative );
    MapFilter_fft<float>  fltr2( step, 1.0, MapFilter_fft<float>::Relative );
    Xmap<float> f1, f2;
    fltr1( f1, xmap );
    fltr2( f2, xmap );

    for ( MRI ix = xmap.first(); !ix.last(); ix.next() )
      test( "MAPFILTER", f1[ix], f2[ix], 0.0001 );
  }

  // test the convolution and fffear objects
  {
    HKL_data<data32::Phi_fom> pw( fsig.hkl_info() );
    HKL_data<data32::F_phi> fp( fsig.hkl_info() );
    pw.compute( abcd, data32::Compute_phifom_from_abcd() );
    fp.compute( fsig, pw, data32::Compute_fphi_from_fsigf_phifom() );
    Grid_sampling grid( spgr, cell, reso, 2.5 );
    Xmap<float> xmap( spgr, cell, grid );
    xmap.fft_from( fp );

    Xmap<float> r1( Spacegroup::p1(), xmap.cell(), xmap.grid_sampling() );
    Xmap<float> r2( Spacegroup::p1(), xmap.cell(), xmap.grid_sampling() );
    int irad = 3;
    clipper::Grid_range tg( clipper::Coord_grid(-irad,-irad,-irad),
			    clipper::Coord_grid( irad, irad, irad) );
    NXmap<float> target( xmap.cell(), xmap.grid_sampling(), tg );
    NXmap<float> weight( xmap.cell(), xmap.grid_sampling(), tg );
    target = weight = 0.0;
    for ( Coord_grid c = tg.min(); !c.last(tg); c.next(tg) ) {
      if ( c*c <= 5 ) {
	target.set_data(c-tg.min(),xmap.get_data(c));
	weight.set_data(c-tg.min(),1.0);
      }
    }

    Convolution_search_slow<float> conv1( xmap );
    Convolution_search_fft<float>  conv2( xmap );
    conv1( r1, target );
    conv2( r2, target );
    for ( MRI ix = r1.first(); !ix.last(); ix.next() )
      test( "CONVOL", r1[ix], r2[ix], 0.0001 );

    FFFear_slow<float> srch1( xmap );
    FFFear_fft<float>  srch2( xmap );
    srch1( r1, target, weight );
    srch2( r2, target, weight );

    for ( MRI ix = r1.first(); !ix.last(); ix.next() )
      test( "FFFEAR", r1[ix], r2[ix], 0.0001 );
  }

  // test anisotropic scaling
  {
    // expand to P1
    Spacegroup spgrp1( Spacegroup::P1 );
    HKL_info hkl1( spgrp1, cell, fsig.hkl_info().resolution(), true );
    HKL_data<data32::F_sigF> fs( hkl1 );
    for ( HRI ih = hkl1.first(); !ih.last(); ih.next() )
      fs[ih] = fsig[ih.hkl()];
    // make data objects
    HKL_data<data32::F_phi>  fp( fs.hkl_info() );
    HKL_data<data32::F_sigF> fs1 = fs;
    HKL_data<data32::F_sigF> fs2 = fs;
    U_aniso_orth u_ref( 0.10, 0.13, 0.17, -0.02, 0.03, -0.04 );
    // simulate aniso data
    for ( HRI ih = fs.first(); !ih.last(); ih.next() )
      if ( !fs[ih].missing() ) {
	Coord_reci_orth c = ih.hkl().coord_reci_orth(cell);
	double s = exp( Util::twopi2() * u_ref.quad_form( c ) );
	fs1[ih].scale(s);
	fs2[ih].scale(1.0/s);
	fp[ih] = data32::F_phi( fs[ih].f(), 0.0 );
      }
    // now attempt scaling
    SFscale_aniso<float>::TYPE F = SFscale_aniso<float>::F;
    SFscale_aniso<float> sfscl;
    sfscl( fs1, fp );
    U_aniso_orth u_wrk1 = sfscl.u_aniso_orth(F);
    sfscl( fp, fs2 );
    U_aniso_orth u_wrk2 = sfscl.u_aniso_orth(F);
    test( "ANISO-O-00", u_ref.mat00(), u_wrk1.mat00(), 1.0e-6 );
    test( "ANISO-O-11", u_ref.mat11(), u_wrk1.mat11(), 1.0e-6 );
    test( "ANISO-O-22", u_ref.mat22(), u_wrk1.mat22(), 1.0e-6 );
    test( "ANISO-O-01", u_ref.mat01(), u_wrk1.mat01(), 1.0e-6 );
    test( "ANISO-O-02", u_ref.mat02(), u_wrk1.mat02(), 1.0e-6 );
    test( "ANISO-O-12", u_ref.mat12(), u_wrk1.mat12(), 1.0e-6 );
    test( "ANISO-C-00", u_ref.mat00(), u_wrk2.mat00(), 1.0e-6 );
    test( "ANISO-C-11", u_ref.mat11(), u_wrk2.mat11(), 1.0e-6 );
    test( "ANISO-C-22", u_ref.mat22(), u_wrk2.mat22(), 1.0e-6 );
    test( "ANISO-C-01", u_ref.mat01(), u_wrk2.mat01(), 1.0e-6 );
    test( "ANISO-C-02", u_ref.mat02(), u_wrk2.mat02(), 1.0e-6 );
    test( "ANISO-C-12", u_ref.mat12(), u_wrk2.mat12(), 1.0e-6 );
  }

  return ( error_count == 0 );
}


} // namespace clipper
