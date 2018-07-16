#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-cns.h>

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

  CNS_HKLfile cns;
  cns.open_write( "1.hkl" );
  cns.export_hkl_info( mydata );
  cns.export_hkl_data( myfsig );
  cns.export_hkl_data( myphwt );
  cns.close_write();

  HKL_info mydata2( mydata.spacegroup(), mydata.cell(), mydata.resolution() );
  HKL_data<F_sigF> myfsig2( mydata2 );
  HKL_data<Phi_fom> myphwt2( mydata2 );
  cns.open_read( "1.hkl" );
  cns.import_hkl_info( mydata2 );
  cns.import_hkl_data( myfsig2 );
  cns.import_hkl_data( myphwt2 );
  cns.close_read();

  cns.open_write( "2.hkl" );
  cns.export_hkl_info( mydata2 );
  cns.export_hkl_data( myfsig2 );
  cns.export_hkl_data( myphwt2 );
  cns.close_write();

  HKL_data<data32::F_phi> fphi( mydata );
  fphi.compute( myfsig, myphwt, data32::Compute_fphi_from_fsigf_phifom() );
  Grid_sampling grid( mydata.spacegroup(), mydata.cell(), Resolution(8.0) );
  Xmap<float> xmap1( mydata.spacegroup(), mydata.cell(), grid );
  Xmap<float> xmap2;
  xmap1.fft_from( fphi );
  CNSMAPfile map( mydata.spacegroup() );
  map.open_write( "1.map" );
  map.export_xmap( xmap1 );
  map.close_write();

  map.open_read( "1.map" );
  map.import_xmap( xmap2 );
  map.close_read();

  map.open_write( "2.map" );
  map.export_xmap( xmap2 );
  map.close_write();

  std::cout << xmap1.spacegroup().symbol_hall() << "\t" << xmap1.grid_sampling().format() << "\n";
  std::cout << xmap2.spacegroup().symbol_hall() << "\t" << xmap2.grid_sampling().format() << "\n";

  Xmap<float>::Map_reference_index ix;
  for ( ix = xmap1.first(); !ix.last(); ix.next() )
    if ( fabs( xmap1[ix] - xmap2[ix] ) > 0.0001 )
      std::cout << ix.coord().format() << "\t" << xmap1[ix] << "\t" << xmap2[ix] << "\n";

  
  // now check cns vs ccp4 files: get the test files from the EDS
  /*
  Xmap<float> x1, x2;
  CCP4MAPfile f1;
  f1.open_read( "1ajr.ccp4" );
  f1.import_xmap( x1 );
  f1.close_read();
  CNSMAPfile  f2( x1.spacegroup() );
  f2.open_read( "1ajr.cns" );
  f2.import_xmap( x2 );
  f2.close_read();
  std::cout << x1.spacegroup().symbol_hall() << "\t" << x1.grid_sampling().format() << "\n";
  std::cout << x2.spacegroup().symbol_hall() << "\t" << x2.grid_sampling().format() << "\n";

  for ( ix = x1.first(); !ix.last(); ix.next() )
    if ( fabs( x1[ix] - x2[ix] ) > 0.0001 )
      std::cout << ix.coord().format() << "\t" << x1[ix] << "\t" << x2[ix] << "\n";
  */
}
