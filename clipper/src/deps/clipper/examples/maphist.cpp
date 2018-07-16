// Clipper-map histogram
/* (C) 2004 Kevin Cowtan */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>


using namespace clipper;
using namespace clipper::data32;



int main(int argc, char** argv)
{
  Xmap<float> xmap;
  Xmap<float> xmsk;
  CCP4MAPfile file;

  file.open_read( argv[1] );
  file.import_xmap( xmap );
  file.close_read();

  xmsk.init( xmap.spacegroup(), xmap.cell(), xmap.grid_sampling() );
  xmsk = 1.0;
  if ( argc > 2 ) {
    file.open_read( argv[2] );
    file.import_xmap( xmsk );
    file.close_read();
  }

  clipper::Range<double> xrng( -1.0, 1.0 );
  clipper::Histogram hist( xrng, 100 );

  clipper::Xmap_base::Map_reference_index ix;
  float s0 = 0.0, s1 = 0.0;
  for ( ix = xmap.first(); !ix.last(); ix.next() ) {
    if ( xmsk[ix] > 0.5 ) {
      s0 += 1.0;
      s1 += xmap[ix];
    }
  }
  s1 /= s0;
  for ( ix = xmap.first(); !ix.last(); ix.next() ) {
    if ( xmsk[ix] > 0.5 ) {
      hist.accumulate( xmap[ix] - s1, 50.0/s0 );
    }
  }
  for ( int i = 0; i < hist.size(); i++ ) {
    std::cout << hist.y(i) << "\n";
  }
}
