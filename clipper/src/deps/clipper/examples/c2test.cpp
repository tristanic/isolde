#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>

void c1()
{
  clipper::CCP4MTZfile mtzin;
  mtzin.open_read( "my.mtz" );        // open new file
  clipper::Spacegroup   spgr = mtzin.spacegroup();
  clipper::Cell         cell = mtzin.cell();
  clipper::Resolution   reso = mtzin.resolution();
  clipper::HKL_info myhkl( spgr, cell, reso, true );  // init hkls
  clipper::HKL_data<clipper::data32::F_phi> fphidata( myhkl );
  mtzin.import_hkl_data( fphidata, "/*/*/[FCAL,PHICAL]" );
  mtzin.close_read();

  clipper::Grid_sampling mygrid( spgr, cell, reso );  // define grid
  clipper::Xmap<float> mymap( spgr, cell, mygrid );  // define map
  mymap.fft_from( fphidata );                  // generate map

  clipper::CCP4MAPfile mapout;
  mapout.open_write( "1.map" );      // write map
  mapout.export_xmap( mymap );
  mapout.close_write();
}

void c2()
{
  clipper::CCP4MTZfile mtzin;
  mtzin.open_read( "my.mtz" );        // open new file
  clipper::Spacegroup   spgr = mtzin.spacegroup();
  clipper::Cell         cell = mtzin.cell();
  clipper::Resolution   reso = mtzin.resolution();
  clipper::HKL_sampling samp = mtzin.hkl_sampling();
  clipper::HKL_data<clipper::data32::F_phi> fphidata( spgr, cell, samp );
  mtzin.import_hkl_data( fphidata, "/*/*/[FCAL,PHICAL]" );
  mtzin.close_read();

  clipper::Grid_sampling mygrid( spgr, cell, reso );  // define grid
  clipper::Xmap<float> mymap( spgr, cell, mygrid );  // define map
  mymap.fft_from( fphidata );                  // generate map

  clipper::CCP4MAPfile mapout;
  mapout.open_write( "2.map" );      // write map
  mapout.export_xmap( mymap );
  mapout.close_write();
}

void c3()
{
  clipper::CCP4MTZfile mtzin;
  mtzin.open_read( "my.mtz" );        // open new file
  clipper::HKL_data<clipper::data32::F_phi> fphidata;
  mtzin.import_hkl_data( fphidata, "/*/*/[FCAL,PHICAL]" );
  mtzin.close_read();

  clipper::Spacegroup   spgr = fphidata.spacegroup();
  clipper::Cell         cell = fphidata.cell();
  clipper::Resolution   reso = fphidata.hkl_info().resolution();

  clipper::Grid_sampling mygrid( spgr, cell, reso );  // define grid
  clipper::Xmap<float> mymap( spgr, cell, mygrid );  // define map
  mymap.fft_from( fphidata );                  // generate map

  clipper::CCP4MAPfile mapout;
  mapout.open_write( "3.map" );      // write map
  mapout.export_xmap( mymap );
  mapout.close_write();
}

int main()
{
  clipper::data::HDcache.set_mode( clipper::ObjectCache<clipper::HKL_data_cacheobj>::MINMEM );
  std::cout << "1\n"; clipper::data::HDcache.debug();
  c1();
  std::cout << "2\n"; clipper::data::HDcache.debug();
  c2();
  std::cout << "3\n"; clipper::data::HDcache.debug();
  c3();
  std::cout << "4\n"; clipper::data::HDcache.debug();
}
