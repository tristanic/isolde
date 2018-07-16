// Clipper test application
/* (C) 2003 Kevin Cowtan */


#include <clipper/clipper.h>
#include <clipper/clipper-phs.h>
#include <clipper/clipper-cif.h>


using clipper::data32::F_sigF;


int main(int argc, char** argv)
{
  clipper::CIFfile ciffile;
  clipper::PHSfile phsfile;

  // make data objects
  clipper::HKL_info wrk_hkl;
  clipper::HKL_data<F_sigF> wrk_fo( wrk_hkl );

  // read data
  ciffile.open_read( argv[1] );
  ciffile.import_hkl_info( wrk_hkl );
  ciffile.import_hkl_data( wrk_fo );
  ciffile.close_read();

  // write data
  phsfile.open_write( argv[2] );
  phsfile.export_hkl_info( wrk_hkl );
  phsfile.import_hkl_data( wrk_fo );
  phsfile.close_read();
}
