#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-cctbx.h>
//#include <clipper/clipper-mmdbold.h>

#include <iostream>


using namespace clipper;


int main()
{
  const char* hallsymbols[] = {
// the 530 tabulated settings
 "P 1", "-P 1", "P 2y", "P 2", "P 2x",
 "P 2yb", "P 2c", "P 2xa", "C 2y", "A 2y",
 "I 2y", "A 2", "B 2", "I 2", "B 2x",
 "C 2x", "I 2x", "P -2y", "P -2", "P -2x",
 "P -2yc", "P -2yac", "P -2ya", "P -2a", "P -2ab",
 "P -2b", "P -2xb", "P -2xbc", "P -2xc", "C -2y",
 "A -2y", "I -2y", "A -2", "B -2", "I -2",
 "B -2x", "C -2x", "I -2x", "C -2yc", "A -2yac",
 "I -2ya", "A -2ya", "C -2ybc", "I -2yc", "A -2a",
 "B -2bc", "I -2b", "B -2b", "A -2ac", "I -2a",
 "B -2xb", "C -2xbc", "I -2xc", "C -2xc", "B -2xbc",
 "I -2xb", "-P 2y", "-P 2", "-P 2x", "-P 2yb",
 "-P 2c", "-P 2xa", "-C 2y", "-A 2y", "-I 2y",
 "-A 2", "-B 2", "-I 2", "-B 2x", "-C 2x",
 "-I 2x", "-P 2yc", "-P 2yac", "-P 2ya", "-P 2a",
 "-P 2ab", "-P 2b", "-P 2xb", "-P 2xbc", "-P 2xc",
 "-P 2ybc", "-P 2yn", "-P 2yab", "-P 2ac", "-P 2n",
 "-P 2bc", "-P 2xab", "-P 2xn", "-P 2xac", "-C 2yc",
 "-A 2yac", "-I 2ya", "-A 2ya", "-C 2ybc", "-I 2yc",
 "-A 2a", "-B 2bc", "-I 2b", "-B 2b", "-A 2ac",
 "-I 2a", "-B 2xb", "-C 2xbc", "-I 2xc", "-C 2xc",
 "-B 2xbc", "-I 2xb", "P 2 2", "P 2c 2", "P 2a 2a",
 "P 2 2b", "P 2 2ab", "P 2bc 2", "P 2ac 2ac", "P 2ac 2ab",
 "C 2c 2", "A 2a 2a", "B 2 2b", "C 2 2", "A 2 2",
 "B 2 2", "F 2 2", "I 2 2", "I 2b 2c", "P 2 -2",
 "P -2 2", "P -2 -2", "P 2c -2", "P 2c -2c", "P -2a 2a",
 "P -2 2a", "P -2 -2b", "P -2b -2", "P 2 -2c", "P -2a 2",
 "P -2b -2b", "P 2 -2a", "P 2 -2b", "P -2b 2", "P -2c 2",
 "P -2c -2c", "P -2a -2a", "P 2c -2ac", "P 2c -2b", "P -2b 2a",
 "P -2ac 2a", "P -2bc -2c", "P -2a -2ab", "P 2 -2bc", "P 2 -2ac",
 "P -2ac 2", "P -2ab 2", "P -2ab -2ab", "P -2bc -2bc", "P 2ac -2",
 "P 2bc -2bc", "P -2ab 2ab", "P -2 2ac", "P -2 -2bc", "P -2ab -2",
 "P 2 -2ab", "P -2bc 2", "P -2ac -2ac", "P 2c -2n", "P 2c -2ab",
 "P -2bc 2a", "P -2n 2a", "P -2n -2ac", "P -2ac -2n", "P 2 -2n",
 "P -2n 2", "P -2n -2n", "C 2 -2", "A -2 2", "B -2 -2",
 "C 2c -2", "C 2c -2c", "A -2a 2a", "A -2 2a", "B -2 -2b",
 "B -2b -2", "C 2 -2c", "A -2a 2", "B -2b -2b", "A 2 -2",
 "B 2 -2", "B -2 2", "C -2 2", "C -2 -2", "A -2 -2",
 "A 2 -2c", "B 2 -2c", "B -2c 2", "C -2b 2", "C -2b -2b",
 "A -2c -2c", "A 2 -2a", "B 2 -2b", "B -2b 2", "C -2c 2",
 "C -2c -2c", "A -2a -2a", "A 2 -2ac", "B 2 -2bc", "B -2bc 2",
 "C -2bc 2", "C -2bc -2bc", "A -2ac -2ac", "F 2 -2", "F -2 2",
 "F -2 -2", "F 2 -2d", "F -2d 2", "F -2d -2d", "I 2 -2",
 "I -2 2", "I -2 -2", "I 2 -2c", "I -2a 2", "I -2b -2b",
 "I 2 -2a", "I 2 -2b", "I -2b 2", "I -2c 2", "I -2c -2c",
 "I -2a -2a", "-P 2 2", "P 2 2 -1n", "-P 2ab 2bc", "-P 2 2c",
 "-P 2a 2", "-P 2b 2b", "P 2 2 -1ab", "-P 2ab 2b", "P 2 2 -1bc",
 "-P 2b 2bc", "P 2 2 -1ac", "-P 2a 2c", "-P 2a 2a", "-P 2b 2",
 "-P 2 2b", "-P 2c 2c", "-P 2c 2", "-P 2 2a", "-P 2a 2bc",
 "-P 2b 2n", "-P 2n 2b", "-P 2ab 2c", "-P 2ab 2n", "-P 2n 2bc",
 "-P 2ac 2", "-P 2bc 2bc", "-P 2ab 2ab", "-P 2 2ac", "-P 2 2bc",
 "-P 2ab 2", "-P 2a 2ac", "-P 2b 2c", "-P 2a 2b", "-P 2ac 2c",
 "-P 2bc 2b", "-P 2b 2ab", "-P 2 2ab", "-P 2bc 2", "-P 2ac 2ac",
 "-P 2ab 2ac", "-P 2ac 2bc", "-P 2bc 2ab", "-P 2c 2b", "-P 2c 2ac",
 "-P 2ac 2a", "-P 2b 2a", "-P 2a 2ab", "-P 2bc 2c", "-P 2 2n",
 "-P 2n 2", "-P 2n 2n", "P 2 2ab -1ab", "-P 2ab 2a", "P 2bc 2 -1bc",
 "-P 2c 2bc", "P 2ac 2ac -1ac", "-P 2c 2a", "-P 2n 2ab", "-P 2n 2c",
 "-P 2a 2n", "-P 2bc 2n", "-P 2ac 2b", "-P 2b 2ac", "-P 2ac 2ab",
 "-P 2bc 2ac", "-P 2ac 2n", "-P 2bc 2a", "-P 2c 2ab", "-P 2n 2ac",
 "-P 2n 2a", "-P 2c 2n", "-C 2c 2", "-C 2c 2c", "-A 2a 2a",
 "-A 2 2a", "-B 2 2b", "-B 2b 2", "-C 2bc 2", "-C 2bc 2bc",
 "-A 2ac 2ac", "-A 2 2ac", "-B 2 2bc", "-B 2bc 2", "-C 2 2",
 "-A 2 2", "-B 2 2", "-C 2 2c", "-A 2a 2", "-B 2b 2b",
 "-C 2b 2", "-C 2b 2b", "-A 2c 2c", "-A 2 2c", "-B 2 2c",
 "-B 2c 2", "C 2 2 -1bc", "-C 2b 2bc", "C 2 2 -1bc", "-C 2b 2c",
 "A 2 2 -1ac", "-A 2a 2c", "A 2 2 -1ac", "-A 2ac 2c", "B 2 2 -1bc",
 "-B 2bc 2b", "B 2 2 -1bc", "-B 2b 2bc", "-F 2 2", "F 2 2 -1d",
 "-F 2uv 2vw", "-I 2 2", "-I 2 2c", "-I 2a 2", "-I 2b 2b",
 "-I 2b 2c", "-I 2a 2b", "-I 2b 2", "-I 2a 2a", "-I 2c 2c",
 "-I 2 2b", "-I 2 2a", "-I 2c 2", "P 4", "P 4w",
 "P 4c", "P 4cw", "I 4", "I 4bw", "P -4",
 "I -4", "-P 4", "-P 4c", "P 4ab -1ab", "-P 4a",
 "P 4n -1n", "-P 4bc", "-I 4", "I 4bw -1bw", "-I 4ad",
 "P 4 2", "P 4ab 2ab", "P 4w 2c", "P 4abw 2nw", "P 4c 2",
 "P 4n 2n", "P 4cw 2c", "P 4nw 2abw", "I 4 2", "I 4bw 2bw",
 "P 4 -2", "P 4 -2ab", "P 4c -2c", "P 4n -2n", "P 4 -2c",
 "P 4 -2n", "P 4c -2", "P 4c -2ab", "I 4 -2", "I 4 -2c",
 "I 4bw -2", "I 4bw -2c", "P -4 2", "P -4 2c", "P -4 2ab",
 "P -4 2n", "P -4 -2", "P -4 -2c", "P -4 -2ab", "P -4 -2n",
 "I -4 -2", "I -4 -2c", "I -4 2", "I -4 2bw", "-P 4 2",
 "-P 4 2c", "P 4 2 -1ab", "-P 4a 2b", "P 4 2 -1n", "-P 4a 2bc",
 "-P 4 2ab", "-P 4 2n", "P 4ab 2ab -1ab", "-P 4a 2a", "P 4ab 2n -1ab",
 "-P 4a 2ac", "-P 4c 2", "-P 4c 2c", "P 4n 2c -1n", "-P 4ac 2b",
 "P 4n 2 -1n", "-P 4ac 2bc", "-P 4c 2ab", "-P 4n 2n", "P 4n 2n -1n",
 "-P 4ac 2a", "P 4n 2ab -1n", "-P 4ac 2ac", "-I 4 2", "-I 4 2c",
 "I 4bw 2bw -1bw", "-I 4bd 2", "I 4bw 2aw -1bw", "-I 4bd 2c", "P 3",
 "P 31", "P 32", "R 3", "P 3*", "-P 3",
 "-R 3", "-P 3*", "P 3 2", "P 3 2\"", "P 31 2c (0 0 1)",
 "P 31 2\"", "P 32 2c (0 0 -1)", "P 32 2\"", "R 3 2\"", "P 3* 2",
 "P 3 -2\"", "P 3 -2", "P 3 -2\"c", "P 3 -2c", "R 3 -2\"",
 "P 3* -2", "R 3 -2\"c", "P 3* -2n", "-P 3 2", "-P 3 2c",
 "-P 3 2\"", "-P 3 2\"c", "-R 3 2\"", "-P 3* 2", "-R 3 2\"c",
 "-P 3* 2n", "P 6", "P 61", "P 65", "P 62",
 "P 64", "P 6c", "P -6", "-P 6", "-P 6c",
 "P 6 2", "P 61 2 (0 0 -1)", "P 65 2 (0 0 1)", "P 62 2c (0 0 1)", "P 64 2c (0 0 -1)",
 "P 6c 2c", "P 6 -2", "P 6 -2c", "P 6c -2", "P 6c -2c",
 "P -6 2", "P -6c 2", "P -6 -2", "P -6c -2c", "-P 6 2",
 "-P 6 2c", "-P 6c 2", "-P 6c 2c", "P 2 2 3", "F 2 2 3",
 "I 2 2 3", "P 2ac 2ab 3", "I 2b 2c 3", "-P 2 2 3", "P 2 2 3 -1n",
 "-P 2ab 2bc 3", "-F 2 2 3", "F 2 2 3 -1d", "-F 2uv 2vw 3", "-I 2 2 3",
 "-P 2ac 2ab 3", "-I 2b 2c 3", "P 4 2 3", "P 4n 2 3", "F 4 2 3",
 "F 4d 2 3", "I 4 2 3", "P 4acd 2ab 3", "P 4bd 2ab 3", "I 4bd 2c 3",
 "P -4 2 3", "F -4 2 3", "I -4 2 3", "P -4n 2 3", "F -4c 2 3",
 "I -4bd 2c 3", "-P 4 2 3", "P 4 2 3 -1n", "-P 4a 2bc 3", "-P 4n 2 3",
 "P 4n 2 3 -1n", "-P 4bc 2bc 3", "-F 4 2 3", "-F 4c 2 3", "F 4d 2 3 -1d",
 "-F 4vw 2vw 3", "F 4d 2 3 -1cd", "-F 4cvw 2vw 3", "-I 4 2 3", "-I 4bd 2c 3",
// the 51 pathological point groups
 "-P 1", "-P 2", "-P 2y", "-P 2x", "-P 2\"", "-P 2y\"", "-P 2x\"", "-P 2'", "-P 2y'", "-P 2x'", "-P 2 2", "-P 2 2\"", "-P 2 2\"(y,z,x)", "-P 2 2\"(z,x,y)", "-P 3", "-P 3 (y,z,x)", "-P 3 (z,x,y)", "-P 3 (-x,y,z)", "-P 3 (y,z,-x)", "-P 3 (z,-x,y)", "-P 3*", "-P 3* (-x,y,z)", "-P 3* (x,-y,z)", "-P 3* (x,y,-z)", "-P 3 2", "-P 3 2 (y,z,x)", "-P 3 2 (z,x,y)", "-P 3* 2", "-P 3* 2 (-x,y,z)", "-P 3* 2 (x,-y,z)", "-P 3* 2 (-x,-y,z)", "-P 3 2\"", "-P 3 2\"(z,x,y)", "-P 3 2\"(y,z,x)", "-P 3 2\"(-x,y,z)", "-P 3 2\"(z,-x,y)", "-P 3 2\"(y,z,-x)", "-P 4", "-P 4 (y,z,x)", "-P 4 (z,x,y)", "-P 4 2", "-P 4 2 (y,z,x)", "-P 4 2 (z,x,y)", "-P 6", "-P 6 (y,z,x)", "-P 6 (z,x,y)", "-P 6 2", "-P 6 2 (y,z,x)", "-P 6 2 (z,x,y)", "-P 2 2 3", "-P 4 2 3",
  };


  // check cells
  Cell c( Cell_descr( 11, 12, 13, 70, 80, 90 ) );
  std::cout << c.format() << "\n" << CCTBX::cell(CCTBX::cell(c)).format() << "\n";

  // check hkls
  HKL hkl( 3, 4, 5 );
  std::cout << hkl.format() << "\t" << CCTBX::hkl(CCTBX::Hkl(hkl)).format() << "\n";
  std::cout << hkl.format() << "\t" << CCTBX::Hkl(CCTBX::hkl(hkl)).format() << "\n";

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
    cctbx::sgtbx::SpaceGroup sgops(hallsymbols[s]);
    Spacegroup s1( Spgr_descr( symbol, Spgr_descr::Hall ) );
    Spacegroup s2 = CCTBX::spacegroup( sgops );
    Spacegroup s3 = CCTBX::spacegroup( CCTBX::spacegroup( s2 ) );

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

    /*
    // build mmdb
    MMDB mmdb( s1, cell );
    clipper::NDBModel model;
    clipper::NDBChain chain;
    clipper::NDBResidue residue;
    clipper::NDBAtom atom;
    model.set_id("A");
    chain.set_id("A");
    residue.set_type("GLY");
    atom = clipper::NDBAtom::null();
    atom.set_occupancy(1.0);
    atom.set_u_iso(0.5);
    clipper::DBModel m1 = mmdb.add_model( model );
    clipper::DBChain c1 = m1.add_chain( chain );
    clipper::DBResidue r1 = c1.add_residue( residue );
    atom.set_element( "C" );
    atom.set_coord_orth( Coord_orth( 12, 8, 5 ) );
    r1.add_atom( atom );
    atom.set_element( "N" );
    atom.set_coord_orth( Coord_orth( 11, 6, 4 ) );
    r1.add_atom( atom );
    atom.set_element( "O" );
    atom.set_coord_orth( Coord_orth( 13, 5, 4 ) );
    r1.add_atom( atom );
    mmdb.finalise_edit();
    clipper::DBAtom_selection atoms = mmdb.select_atoms_serial( 0, 0 );
    */

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
    if ( s1.hash() != s2.hash() || s1.hash() != s3.hash() || nerr > 1 )
      std::cout << "Fail ";
    else
      std::cout << "OK   ";
    String name = symbol + "                  ";
    std::cout << name.substr(0,17) << cell.alpha_deg() << " " << cell.beta_deg() << " " << cell.gamma_deg() << "\t" << s1.asu_max().format() << "   " << nerr << "\n";
  }
}
