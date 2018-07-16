// Clipper create missing data demo
/* (C) 2003 Kevin Cowtan */
// This is more of a demo application than a serious version
/*
Usage: anisodemo input.mtz input_columns output.mtz fraction
e.g.
./anisodemo ../test/testfile.mtz '/native/peak/[FP,SIGFP]' out.mtz 0.05
*/

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>

extern "C" {
#include <stdlib.h>
}


int main( int argc, char** argv )
{
  // make data object
  clipper::HKL_info hkls;
  clipper::HKL_data<clipper::data32::F_sigF> fsig( hkls );

  // read data from mtz
  clipper::CCP4MTZfile mtz;
  mtz.open_read( argv[1] );
  mtz.import_hkl_info( hkls );
  mtz.import_hkl_data( fsig, argv[2] );
  mtz.close_read();

  // get fraction
  double frac = clipper::String( argv[4] ).f();

  clipper::HKL_info::HKL_reference_index ih;
  for ( ih = hkls.first(); !ih.last(); ih.next() )
    if ( 0.0001*double(random()%10000) < frac )
      fsig[ih] = clipper::data32::F_sigF();

  // write data to mtz
  mtz.open_append( argv[1], argv[3] );
  mtz.export_hkl_data( fsig, "/*/*/missing" );
  mtz.close_append();
}
