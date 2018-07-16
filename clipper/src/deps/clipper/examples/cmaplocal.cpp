// Clipper app to calc local map moments
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include <stdlib.h>


int main( int argc, char** argv )
{
  CCP4Program prog( "cmaplocal", "0.1", "$Date: 2004/05/01" );

  // defaults
  clipper::String ipfile = "NONE";
  clipper::String opfile1 = "NONE";
  clipper::String opfile2 = "NONE";
  double statsrad = -1.0;

  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-mapin" ) {
      if ( ++arg < args.size() ) ipfile = args[arg];
    } else if ( args[arg] == "-mapout-1" ) {
      if ( ++arg < args.size() ) opfile1 = args[arg];
    } else if ( args[arg] == "-mapout-2" ) {
      if ( ++arg < args.size() ) opfile2 = args[arg];
    } else if ( args[arg] == "-radius" ) {
      if ( ++arg < args.size() ) statsrad = clipper::String(args[arg]).f();
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: cmaplocal\n\t-mapin <filename>\n\t-mapout-1 <filename>\n\t-mapout-2 <filename>\n\t-radius <radius>\n";
    exit(1);
  }

  clipper::CCP4MAPfile file;
  clipper::Xmap<float> xmap;
  file.open_read( ipfile );
  file.import_xmap( xmap );
  file.close_read();

  // make squared map
  clipper::Xmap<float> lmom1( xmap );
  clipper::Xmap<float> lmom2( xmap );
  for ( clipper::Xmap<float>::Map_reference_index ix = lmom2.first();
	!ix.last(); ix.next() )
    lmom2[ix] = clipper::Util::sqr( lmom2[ix] );

  // now calculate local mom1, local mom1 squared
  clipper::MapFilterFn_step fn( statsrad );
  clipper::MapFilter_fft<float> fltr( fn, 1.0, clipper::MapFilter_fft<float>::Relative );
  fltr( lmom1, lmom1 );
  fltr( lmom2, lmom2 );

  // calculate std deviation
  for ( clipper::Xmap<float>::Map_reference_index ix = lmom1.first();
	!ix.last(); ix.next() )
    lmom2[ix] = sqrt( lmom2[ix] - clipper::Util::sqr( lmom1[ix] ) );

  // output map
  file.open_write( opfile1 );
  file.export_xmap( lmom1 );
  file.close_write();
  file.open_write( opfile2 );
  file.export_xmap( lmom2 );
  file.close_write();
}
