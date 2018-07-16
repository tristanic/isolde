#include <clipper/clipper.h>
//#include <csymlib.h>


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

  for ( int s = 0; s < sizeof(hallsymbols)/sizeof(hallsymbols[0]); s++ ) {
    // get spacegroup
    String symbol = hallsymbols[s];
    Spacegroup sg( Spgr_descr( symbol, Spgr_descr::Hall ) );

    Spgr_descr::Symop_codes ops;
    for ( int i = 0; i < sg.num_symops(); i++ )
      ops.push_back( Symop_code( sg.symop(i) ) );

    Spgr_descr::Symop_codes cens = ops.centering_ops();
    Spgr_descr::Symop_codes prms = ops.primitive_ops();
    Spgr_descr::Symop_codes gens;
    Spgr_descr::Symop_codes gend = gens.expand();
    for ( int i = 1; i < cens.size(); i++ ) {
      for ( int j = 0; j < gend.size(); j++ )
	if ( cens[i] == gend[j] ) goto skip1;
      gens.push_back( cens[i] );
      gend = gens.expand();
      if ( gend.size() == cens.size() ) break;  // optional optimisation
    skip1:;
    }
    int ncen = gens.size();

    for ( int i = 1; i < prms.size(); i++ ) {
      for ( int j = 0; j < gend.size(); j++ )
	if ( prms[i] == gend[j] ) goto skip2;
      gens.push_back( prms[i] );
      gend = gens.expand();
      if ( gend.size() == ops.size() ) break;  // optional optimisation
    skip2:;
    }
    
  back:
    for ( int i = ncen; i < gens.size(); i++ ) {
      Spgr_descr::Symop_codes genp = gens;
      genp.erase( genp.begin() + i );
      if ( genp.expand().size() == ops.size() ) {
	gens = genp;
	goto back;
      }
    }

    Spgr_descr::Symop_codes symp;
    for ( int i = ncen; i < gens.size(); i++ ) symp.push_back( gens[i] );

    //std::cout << sg.num_symops() << "\t" << gens.size() << " " << sg.generator_ops().size() << "\t" << (sg.symbol_hm()+"        ").substr(0,12) << "\t: " << (sg.symbol_hall()+"          ").substr(0,12) << "\t\t" << symp.expand().size() << " " << sg.num_primops() << " Inv" << sg.invariant_under_change_of_hand() << "\n";
    //std::cout << (sg.symbol_hm()+"         ").substr(0,12) << " Laue: " << Spacegroup( Spgr_descr( ops.laue_ops() ) ).symbol_hm() << " \tPatt: " << Spacegroup( Spgr_descr( ops.patterson_ops() ) ).symbol_hm() << "\n";
    //for ( int i = 0; i < gens.size(); i++ )
    //std::cout << i << " " << gens[i].symop().format() << "\n";

    /* CCP4 tests:
    CSym::CCP4SPG *spacegroup;
    spacegroup = CSym::ccp4spg_load_by_spgname( sg.symbol_hall().c_str() );
    free(spacegroup);
    int h,k,l;
    for (h=0;h<=2;h++)
      for (k=-2;k<=2;k++)
	for (l=-2;l<=2;l++) {
	  HKL rfl(h,k,l);
	  std::cout << rfl.format() << " " << ccp4spg_is_centric(spacegroup,h,k,l) << " " << " " << "\n";
	} */

    // check inversion ops
    std::cout << (sg.symbol_hm()+"         ").substr(0,12)
	      << sg.order_of_symmetry_about_axis(Spacegroup::A) 
	      << sg.order_of_symmetry_about_axis(Spacegroup::B) 
	      << sg.order_of_symmetry_about_axis(Spacegroup::C)
	      << "\t";
    std::cout << (((sg.num_primitive_symops()/sg.order_of_symmetry_about_axis(Spacegroup::C))%2)?0:2) << "\t";
    for ( int i = 0; i < sg.num_primops(); i++ )
      std::cout << sg.primitive_symop(i).rot().det();
    std::cout << "\n";
  }
}
