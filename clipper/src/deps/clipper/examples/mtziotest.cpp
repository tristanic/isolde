#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
//#include <clipper/clipper-ccp4.h>
#include <iostream>



using namespace clipper;
using namespace clipper::data32;


int main(int argc, char** argv)
{
  CCP4MTZfile file;

  HKL_info mydata;
  HKL_data<F_sigF> myfsig(mydata);
  file.open_read(argv[1]);
  file.import_hkl_info( mydata, false );
  file.import_hkl_data( myfsig, argv[2] );
  Spacegroup sp = file.spacegroup();
  file.close_read();

  file.open_append(argv[1],"junk1.mtz");
  file.export_hkl_data( myfsig, "*/*/anom");
  file.close_append();
  file.open_write("junk2.mtz");
  file.export_hkl_info( mydata );
  file.export_hkl_data( myfsig, "*/*/anom");
  file.close_write();

  sp.debug();
}
