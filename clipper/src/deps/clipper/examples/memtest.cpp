#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

#include <iostream>


using namespace clipper;


int main()
{
  const char* hallsymbols[] = {
// the 530 tabulated settings
 "P 1",
 "P 2yb",
 "P 2 2b",
 "C 2c -2",
 "F -2 -2",
 "-P 2ab 2ac",
 "P 4c",
 "P 31",
 "-R 3",
 "-P 3 2\"",
 "P 64",
 "I 2 2 3",
 "-F 4vw 2vw 3",
  };


  // test fft in each spacegroup
  Cell cellc( Cell_descr( 37, 37, 37 ) );
  Cell cellha( Cell_descr( 37, 37, 37, 120, 90, 90 ) );
  Cell cellhb( Cell_descr( 37, 37, 37, 90, 120, 90 ) );
  Cell cellhc( Cell_descr( 37, 37, 37, 90, 90, 120 ) );
  Cell cellha1( Cell_descr( 37, 37, 37, 60, 90, 90 ) );
  Cell cellhb1( Cell_descr( 37, 37, 37, 90, 60, 90 ) );
  Cell cellhc1( Cell_descr( 37, 37, 37, 90, 90, 60 ) );
  Cell cell;
  for ( int s = 0; s < sizeof(hallsymbols)/sizeof(hallsymbols[0]); s++ ) {
    // get spacegroup
    String symbol = hallsymbols[s];
    Spacegroup s1( Spgr_descr( symbol, Spgr_descr::Hall ) );

    // identify trigonal/hexagonal groups
    cell = cellc;
    for ( int sym = 1; sym < s1.num_symops(); sym++ ) {
      if ( ( s1.symop(sym).rot()(1,1) * s1.symop(sym).rot()(1,2) == -1 ) ||
	   ( s1.symop(sym).rot()(2,1) * s1.symop(sym).rot()(2,2) == -1 ) )
	cell = cellha;
      if ( ( s1.symop(sym).rot()(0,0) * s1.symop(sym).rot()(0,2) == -1 ) ||
	   ( s1.symop(sym).rot()(2,0) * s1.symop(sym).rot()(2,2) == -1 ) )
	cell = cellhb;
      if ( ( s1.symop(sym).rot()(0,0) * s1.symop(sym).rot()(0,1) == -1 ) ||
	   ( s1.symop(sym).rot()(1,0) * s1.symop(sym).rot()(1,1) == -1 ) )
	cell = cellhc;
      if ( ( s1.symop(sym).rot()(1,1) * s1.symop(sym).rot()(1,2) == 1 ) ||
	   ( s1.symop(sym).rot()(2,1) * s1.symop(sym).rot()(2,2) == 1 ) )
	cell = cellha1;
      if ( ( s1.symop(sym).rot()(0,0) * s1.symop(sym).rot()(0,2) == 1 ) ||
	   ( s1.symop(sym).rot()(2,0) * s1.symop(sym).rot()(2,2) == 1 ) )
	cell = cellhb1;
      if ( ( s1.symop(sym).rot()(0,0) * s1.symop(sym).rot()(0,1) == 1 ) ||
	   ( s1.symop(sym).rot()(1,0) * s1.symop(sym).rot()(1,1) == 1 ) )
	cell = cellhc1;
    }

    // build model
    Atom_list atoms;
    Atom atom = Atom::null();
    atom.occupancy() = 1.0;
    atom.u_iso() = 0.5;
    atom.element() = "C";
    atom.coord_orth() = Coord_orth( 12, 8, 5 );
    atoms.push_back( atom );
    atom.element() = "N";
    atom.coord_orth() = Coord_orth( 11, 6, 4 );
    atoms.push_back( atom );
    atom.element() = "O";
    atom.coord_orth() = Coord_orth( 13, 5, 4 );
    atoms.push_back( atom );

    // calc structure factors by slow fft
    HKL_info hkl_info( s1, cell, Resolution( 3.0 ), true );
    HKL_data<data32::F_phi> f_phi1( hkl_info );
    HKL_data<data32::F_phi> f_phi2( hkl_info );
    SFcalc_iso_sum<float>( f_phi1, atoms );
    SFcalc_iso_fft<float>( f_phi2, atoms );
    Grid_sampling grid( s1, cell, Resolution( 3.0 ) );
    Xmap<float> xmap( s1, cell, grid );
    FFTmap fftmap( s1, cell, grid );
    xmap.fft_from( f_phi2, Xmap_base::Normal );
    xmap.fft_to( f_phi2, Xmap_base::Sparse );
    fftmap.fft_rfl_to_map( f_phi2, xmap );
    xmap.fft_to( f_phi2, Xmap_base::Normal );
    xmap.fft_from( f_phi2, Xmap_base::Sparse );
    fftmap.fft_map_to_rfl( xmap, f_phi2 );

    int nerr = 0;
    for (int h=-2; h<=2; h++)
      for (int k=-2; k<=2; k++)
	for (int l=-2; l<=2; l++) {
	  HKL rfl(h,k,l);
	  if ( !HKL_class( s1, rfl ).sys_abs() )
	    if ( std::abs( std::complex<float>(fftmap.get_recip_data(rfl) ) -
			   std::complex<float>(f_phi2[rfl]) ) > 1.0e-3 ) {
	      nerr++;
	      std::cout << rfl.format() <<
		"\t: " << fftmap.get_recip_data(rfl).f() <<
		"\t" << fftmap.get_recip_data(rfl).phi() <<
		"\t: " << f_phi2[rfl].f() <<
		"\t" << f_phi2[rfl].phi() <<
		"\t" << f_phi2[rfl].missing() << ":" << HKL_class( s1, rfl ).sys_abs() << "\n";
	    }
	  if ( HKL_class( s1, rfl ).sys_abs() != f_phi2[rfl].missing() ) {
	    nerr++;
	    int sym; bool friedel;
	    int ih = hkl_info.index_of( hkl_info.find_sym( rfl, sym, friedel) );
	    std::cout << rfl.format() << ih << " " << sym << " " << friedel <<
	      "\t" << hkl_info.hkl_of( ih ).format() <<
	      "\t" << f_phi2[rfl].missing() <<
	      ":" << HKL_class( s1, rfl ).sys_abs() << "\n";
	  }
	}
    
    HKL_info::HKL_reference_index ih;
    float tol = 0.002 * f_phi1[ HKL( 0, 0, 0 ) ].f();
    for ( ih = hkl_info.first(); !ih.last(); ih.next() )
      if ( std::abs( std::complex<float>(f_phi1[ih]) -
		     std::complex<float>(f_phi2[ih]) ) > tol ) {
	nerr++;
	std::cout << ih.hkl().format() << 
	  "\t: " << f_phi1[ih].f() << "\t" << f_phi1[ih].phi() <<
	  "\t: " << f_phi2[ih].f() << "\t" << f_phi2[ih].phi() <<
	  "\n";
      }

    // diagnostics
    std::cout << "OK   ";
    String name = symbol + "                  ";
    std::cout << name.substr(0,17) << cell.alpha_deg() << " " << cell.beta_deg() << " " << cell.gamma_deg() << "\t" << s1.asu_max().format() << "   " << nerr << "\n";
  }
}
