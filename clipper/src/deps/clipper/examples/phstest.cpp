#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-phs.h>

#include <iostream>


using namespace clipper;
using namespace clipper::data32;


int main(int argc, char** argv)
{
  CCP4MTZfile file;

  // import an mtz
  HKL_info mydata;
  HKL_data<F_sigF> myfsig( mydata );
  HKL_data<Phi_fom> myphwt( mydata );
  MTZcrystal xtl;
  MTZdataset set;

  file.open_read( argv[1] );
  file.import_hkl_info( mydata, false );
  file.import_hkl_data( myfsig, set, xtl, "*/*/[FP SIGFP]" );
  file.import_hkl_data( myphwt, set, xtl, "*/*/[PHIB FOM]" );
  file.close_read();

  PHSfile phs;
  phs.open_write( "1.phs" );
  phs.export_hkl_info( mydata );
  phs.export_hkl_data( myfsig );
  phs.export_hkl_data( myphwt );
  phs.close_write();

  HKL_info mydata2( mydata.spacegroup(), mydata.cell(), mydata.resolution() );
  HKL_data<F_sigF> myfsig2( mydata2 );
  HKL_data<Phi_fom> myphwt2( mydata2 );
  phs.open_read( "1.phs" );
  phs.import_hkl_info( mydata2 );
  phs.import_hkl_data( myfsig2 );
  phs.import_hkl_data( myphwt2 );
  phs.close_read();

  phs.open_write( "2.phs" );
  phs.export_hkl_info( mydata2 );
  phs.export_hkl_data( myfsig2 );
  phs.export_hkl_data( myphwt2 );
  phs.close_write();
}
