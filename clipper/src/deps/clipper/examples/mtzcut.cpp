/* Stupid little application to combine two sets of map coefficients
   by taking the low resolution terms from one set and the high
   resolution terms from another. */

#include "mtz_io.h"

using namespace clipper;
using namespace clipper::data32;


int main()
{
  CCP4MTZfile file;
  MTZcrystal xtal;
  MTZdataset dset;

  HKL_info mydata;
  HKL_data<F_phi> fphi1( mydata );
  HKL_data<F_phi> fphi2( mydata );

  String line, filename_in, filename_out, colname_lo, colname_hi, colname_out;
  float resol = 5.0;
  filename_out = "mtzcut.mtz";
  colname_out = "cut";

  // read and parse input lines
  while ( std::getline( std::cin, line ), !std::cin.eof() ) {
    std::vector<String> tokens = line.split(" ");
    if ( tokens[0] == "inputfile" )  filename_in  = tokens[1];
    if ( tokens[0] == "outputfile" ) filename_out = tokens[1];
    if ( tokens[0] == "F_phi_lo" )   colname_lo   = tokens[1];
    if ( tokens[0] == "F_phi_hi" )   colname_hi   = tokens[1];
    if ( tokens[0] == "F_phi_out" )  colname_out  = tokens[1];
    if ( tokens[0] == "resolution" ) resol        = tokens[1].f();
  }

  file.open_read( filename_in );
  file.import_hkl_info( mydata, false );
  file.import_hkl_data( fphi1, dset, xtal, colname_lo );
  file.import_hkl_data( fphi2, dset, xtal, colname_hi );
  file.close_read();

  HKL_data<F_phi> fphi3( mydata );
  HKL_info::HKL_reference_index ih;
  for ( ih = mydata.first(); !ih.last(); ih.next() )
    fphi3[ih] = (ih.invresolsq() < pow(resol,-2)) ? fphi1[ih] : fphi2[ih];

  file.open_append( filename_in, filename_out );
  file.export_hkl_data( fphi3, dset, xtal, colname_out );
  file.close_append();
}  
