// Clipper-map test
/* (C) 2002 Kevin Cowtan */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
//#include <clipper/clipper-ccp4.h>

#include <iostream>


using namespace clipper;
using namespace clipper::data32;



int main(int argc, char** argv)
{
  Cell cell = Cell( Cell_descr( 10,10,10,90,90,90 ) );
  Grid grid = Grid( 20, 20, 20 );
  Grid_range gmap = Grid_range( cell, grid, 5.0 );
  std::cout << cell.format() << "\t" << gmap.format() << "\n";
  cell = Cell( Cell_descr( 10,10,10,90,120,90 ) );
  gmap = Grid_range ( cell, grid, 5.0 );
  std::cout << cell.format() << "\t" << gmap.format() << "\n";
  grid = Grid( 20, 30, 40 );
  gmap = Grid_range ( cell, grid, 5.0 );
  std::cout << cell.format() << "\t" << gmap.format() << "\n";
  gmap = Grid_range ( cell, grid, 2.5 );
  std::cout << cell.format() << "\t" << gmap.format() << "\n";

  Xmap<float> xmap;
  CCP4MAPfile file;

  file.open_read( argv[1] );
  file.import_xmap( xmap );
  file.close_read();
  
  Coord_grid c;
  for ( c.w() = 0; c.w() < 2; c.w()++ ) {
    std::cout << c.w() << "\n";
    for ( c.v() = 0; c.v() < 24; c.v()++ ) {
      for ( c.u() = 0; c.u() < 24; c.u()++ ) {
	std::cout.width(5);
	std::cout << rint(1000*xmap.get_data(c)) << " ";
      }
      std::cout << "\n";
    }
  }

  file.open_write( "out.xmap" );
  file.export_xmap( xmap );
  file.close_write();

  NXmap<float> nxmap;
  file.open_read( argv[1] );
  file.import_nxmap( nxmap );
  file.close_read();  

  for ( c.v() = 0; c.v() < 2; c.v()++ ) {
    std::cout << c.v() << "\n";
    for ( c.w() = 0; c.w() < 24; c.w()++ ) {
      for ( c.u() = 0; c.u() < 24; c.u()++ ) {
	std::cout.width(5);
	std::cout << rint(1000*nxmap.get_data(c)) << " ";
      }
      std::cout << "\n";
    }
  }

  file.open_write( "out.nxmap" );
  file.export_nxmap( nxmap );
  file.close_write();

  RTop_orth rtop( RTop_orth::identity() );
  Grid_range g( Coord_grid(10,10,10), Coord_grid(20,20,20) );
  nxmap.init( xmap.cell(), xmap.grid_sampling(), g );
  NX_operator nxop( xmap, nxmap, rtop );

  nxop.debug();

  NXmap_base::Map_reference_index ix;
  for ( ix = nxmap.first(); !ix.last(); ix.next() )
    nxmap[ix] = nxop.xmap_data<Interp_linear,float>( xmap, ix.coord() );

  for ( c = g.min(); !c.last( g ); c.next( g ) )
    if ( c.index( g ) % 100 == 0 )
      std::cout << c.format() << " " << xmap.get_data( c ) << " " << nxop.nxmap_data<Interp_linear,float>( nxmap, c ) << "\n";

}
