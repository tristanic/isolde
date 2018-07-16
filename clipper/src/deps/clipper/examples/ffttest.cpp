// Clipper-fft utility.
/* (C) 2000 Kevin Cowtan */
// This is more of a demo application than a serious version


#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <iostream>

extern "C" {
#include <time.h>
}


using namespace clipper;
using namespace clipper::data32;



int main()
{
  String filename, cmapname( "out.map" );
  String dataname;
  double rlim = 0.0;
  Grid_sampling grid(24,24,24);

  String line = "";

  // read and parse input lines
  while ( std::getline( std::cin, line ), !std::cin.eof() ) {
    std::vector<String> tokens = line.split(" ");

    // input file
    if ( tokens[0] == "inputfile" ) {
      filename = tokens[1];
    }
    // output file
    if ( tokens[0] == "outputfile" ) {
      cmapname = tokens[1];
    }
    // import a column list
    else if ( tokens[0] == "F_phi" ) {
      dataname = tokens[1];
    }
    // resolution
    else if ( tokens[0] == "resolution" ) {
      rlim = tokens[1].f();
    }
    // grid
    else if ( tokens[0] == "grid" ) {
      grid = Grid_sampling( tokens[1].i(), tokens[2].i(), tokens[3].i() );
    }

  }

  std::cout << filename << "\n" << cmapname << "\n" << dataname << "\n" << rlim << "\n";

  // Make data objects: Need spacegroup, cell, hkls, and data
  CSpacegroup spgr( "base spgr" );
  CCell cell( spgr, "base cell" );
  CResolution reso( cell, "base reso" );
  CHKL_info rfl( reso, "cfft" );
  CHKL_data<F_phi> fphidata( rfl, "" );

  CCP4MTZfile mtzin; mtzin.open_read( filename );   // open new file
  spgr.init( mtzin.spacegroup() );              // get info from file
  cell.init( mtzin.cell() );
  reso.init( Resolution( Util::max( rlim, mtzin.resolution().limit() ) ) );
  rfl.generate_hkl_list();
  mtzin.import_chkl_data( fphidata, dataname ); // read data
  mtzin.close_read();

  std::cout << "Number of reflections: " << rfl.num_reflections() << "\n";

  // make map
  Xmap<float> xmap( rfl.spacegroup(), rfl.cell(), grid );
  Xmap<float> xmap1( rfl.spacegroup(), rfl.cell(), grid );
  Xmap<float> xmap2( rfl.spacegroup(), rfl.cell(), grid );

  // fft
  xmap.fft_from( fphidata );

  // write map
  CCP4MAPfile mapout; mapout.open_write( cmapname );
  mapout.export_xmap( xmap );
  mapout.close_write();

  // now bench the fft methods
  int t0 = time( NULL );
  for ( int i = 0; i < 20; i++ )
    xmap1.fft_from( fphidata, Xmap_base::Normal );
  int t1 = time( NULL );
  for ( int i = 0; i < 20; i++ )
    xmap2.fft_from( fphidata, Xmap_base::Sparse );
  int t2 = time( NULL );
  std::cout << "Times: " << t2-t1 << "\t" << t1-t0 << "\n";

  Xmap<float>::Map_reference_index ix;
  for ( ix = xmap1.first(); !ix.last(); ix.next() )
    if ( ix.index() % 100 == 0 ) std::cout << ix.coord().format() << "\t" << xmap1[ix] << "\t" << xmap2[ix] << "\n";

  // do some diagnostics
  Coord_grid c;
  for ( c.w() = 0; c.w() < 2; c.w()++ ) {
    for ( c.v() = 0; c.v() < 24; c.v()++ ) {
      for ( c.u() = 0; c.u() < 24; c.u()++ ) {
	std::cout.width(5);
	std::cout << rint(1000*xmap.get_data(c)) << " ";
      }
      std::cout << "\n";
    }
    std::cout << c.w() << "\n";
  }
  std::cout << "at 012 " << xmap.get_data( Coord_grid( 0, 1, 2 ) ) << "\n";
  std::cout << "at 112 " << xmap.get_data( Coord_grid( 1, 1, 2 ) ) << "\n";
  std::cout << "at 212 " << xmap.get_data( Coord_grid( 2, 1, 2 ) ) << "\n";
  for ( ftype u = 0; u < 2.01; u+=0.1 ) {
    Coord_frac f( u/24.0, 1.0/24.0, 2.0/24.0 );
    float v1; Grad_frac<float> g1;
    float v2; Grad_frac<float> g2;
    xmap.interp_grad<Interp_cubic>(f,v1,g1);
    xmap.interp_grad<Interp_cubic>(f,v2,g2);
    std::cout.precision(5);
    std::cout << u << " " << xmap.interp<Interp_linear>(f) << " " << v1 << " " << xmap.interp<Interp_cubic>(f) << " " << g1[0]/24 << " " << g2[0]/24 << " " << g1[2]/24 << " " << g2[2]/24 << "\n";
  }
  std::cout << "---\n";
  for ( ftype u = 0; u < 2.01; u+=0.1 ) {
    Coord_frac f( u/24.0, 1.0+u/36.0, 2.0+u/48.0 );
    float v1; Grad_frac<float> g1;
    float v2; Grad_frac<float> g2; Curv_frac<float> c2;
    xmap.interp_grad<Interp_cubic>(f,v1,g1);
    xmap.interp_curv<Interp_cubic>(f,v2,g2,c2);
    Grad_orth<float> go = g2.grad_orth(cell);
    std::cout.precision(5);
    float v; Interp_linear::interp(xmap,f.coord_map(grid),v);
    std::cout << u << " " << v << " " << v1 << " " << v2 << " " << g1[0]/24 << " " << g2[0]/24 << " " << g1[1]/24 << " " << g2[1]/24 << "\n";
  }

  // test interpolation
  Coord_map m0( 1.1, 2.2, 3.3 ), m1( 1.11, 2.2, 3.3 ), m2( 1.1, 2.21, 3.3 ), m3( 1.1, 2.2, 3.31 );
  float v0; Grad_map<float> g0, g1, g2, g3; Curv_map<float> c0;
  Interp_cubic::interp_curv<float>( xmap, m0, v0, g0, c0 );
  Interp_cubic::interp_grad<float>( xmap, m1, v0, g1 );
  Interp_cubic::interp_grad<float>( xmap, m2, v0, g2 );
  Interp_cubic::interp_grad<float>( xmap, m3, v0, g3 );
  std:: cout << "Interp curv\n" << c0.format() << "\n" << (100.0f*(g1-g0)).format() << "\n" << (100.0f*(g2-g0)).format() << "\n"<< (100.0f*(g3-g0)).format() << "\n";

  // benchmark the interpolators
  //for ( ftype u = 0; u < 20.01; u+=0.00001 ) {
  for ( ftype u = 0; u < 20.01; u+=0.01 ) {
    Coord_frac f( u/24.0, 1.0+u/36.0, 2.0+u/48.0 );
    //Coord_map g = grid.to_map(f);
    ftype x = xmap.interp<Interp_linear>(f);
    //ftype x = Interp_linear<Xmap,float>::interp(xmap,g);
    //ftype y = Interp_linear<Xmap<float> >::interp(xmap,g);
    ftype y = xmap.interp<Interp_linear>(f);
    if ( fabs( x-y ) > 1.0e-6 ) std::cout << u << " " << x << " " << y << "\n";
  }

  // now try making some grids
  Spacegroup sg(Spgr_descr(76));
  for ( ftype x = 30; x < 66; x+=5 ) {
    Cell c(Cell_descr(x,33,33));
    Grid_sampling g( sg, c, Resolution(2.0), 1.0 );
    std::cout << x << " " << g.format() << "\n";
  }

  // now test the fffear objects
  Xmap<float> r1( Spacegroup::p1(), xmap.cell(), xmap.grid_sampling() );
  Xmap<float> r2( Spacegroup::p1(), xmap.cell(), xmap.grid_sampling() );
  int irad = 2;
  clipper::Grid_range tg( clipper::Coord_grid(-irad,-irad,-irad),
			clipper::Coord_grid(irad,irad,irad) );
  NXmap<float> target( xmap.cell(), xmap.grid_sampling(), tg );
  NXmap<float> weight( xmap.cell(), xmap.grid_sampling(), tg );
  target = weight = 0.0;
  for ( Coord_grid c = tg.min(); !c.last(tg); c.next(tg) ) {
    if ( c*c <= 2 ) {
      target.set_data(c-tg.min(),xmap.get_data(c));
      weight.set_data(c-tg.min(),1.0);
    }
  }

  //FFFear_slow_basic srch1( xmap );
  //FFFear_fft_basic  srch2( xmap );
  FFFear_slow<float> srch1( xmap );
  FFFear_fft<float>  srch2( xmap );
  srch1( r1, target, weight );
  srch2( r2, target, weight );

  for ( ix = r1.first(); !ix.last(); ix.next() )
    if ( fabs(r1[ix] - r2[ix]) > 0.001 )
      std::cout << ix.coord().format() << "  \t" << r1[ix] << "  \t" << r2[ix] << "\n";

  MapFilterFn_step step( 2.5 );
  MapFilter_slow<float> fltr1( step, 1.0, MapFilter_slow<float>::Relative );
  MapFilter_fft<float>  fltr2( step, 1.0, MapFilter_fft<float>::Relative );
  Xmap<float> f1, f2;
  fltr1( f1, xmap );
  fltr2( f2, xmap );

  for ( ix = xmap.first(); !ix.last(); ix.next() )
    if ( fabs(f1[ix] - f2[ix]) > 0.001 )
      std::cout << ix.coord().format() << "  \t" << f1[ix] << "  \t" << f2[ix] << "  \t" << xmap[ix] << "\n";

  // test
  SFweight_spline<float> sfw;
  sfw.debug();
}
