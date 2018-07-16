// Clipper test application
/* (C) 2003 Kevin Cowtan */


#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>


using clipper::data32::F_sigF;
using clipper::data32::F_phi;


int main(int argc, char** argv)
{
  const int n_scl_param = 6;

  clipper::CCP4MTZfile mtzfile;
  clipper::MTZcrystal xtal;
  clipper::MTZdataset dset;

  // make data objects
  clipper::HKL_info wrk_hkl;
  clipper::HKL_data<F_sigF> wrk_fo( wrk_hkl );
  clipper::HKL_data<F_phi> wrk_fc( wrk_hkl );

  // read data
  mtzfile.open_read( argv[1] );
  mtzfile.import_hkl_info( wrk_hkl );
  mtzfile.import_hkl_data( wrk_fo, dset, xtal, argv[2] );
  mtzfile.import_hkl_data( wrk_fc, dset, xtal, argv[3] );
  for ( int i = 0; i < mtzfile.column_labels().size(); i++ ) std::cout << i << mtzfile.column_labels()[i] << "\n";
  mtzfile.close_read();

  mtzfile.spacegroup().debug();
  std::cout << mtzfile.cell().format() << "\n";
  std::cout << mtzfile.resolution().limit() << "\n";
  wrk_hkl.debug();
  for ( clipper::HKL_info::HKL_reference_index ih = wrk_hkl.first(); !ih.last(); ih.next() ) std::cout << ih.hkl().format() << "\t" << wrk_fo[ih].f() << "\t" << wrk_fo[ih].sigf() << "\t" << wrk_fc[ih].f() << "\t" << wrk_fc[ih].phi() << "\n";

  // scale and difference data
  std::vector<clipper::ftype> params( n_scl_param, 1.0 );
  clipper::BasisFn_spline wrk_basis( wrk_hkl, n_scl_param, 2.0 );
  clipper::TargetFn_scaleF1F2<F_phi,F_sigF> wrk_target( wrk_fc, wrk_fo );
  clipper::ResolutionFn wrk_scale( wrk_hkl, wrk_basis, wrk_target, params );
  clipper::HKL_info::HKL_reference_index ih;
  for ( ih = wrk_hkl.first(); !ih.last(); ih.next() )
    if ( !wrk_fc[ih].missing() ) {
      wrk_fc[ih].scale( sqrt( wrk_scale.f(ih) ) );  // scale
      wrk_fc[ih].f() -= wrk_fo[ih].f();             // difference
    }

  // calculate map
  clipper::Grid_sampling wrk_grid( wrk_hkl.spacegroup(), wrk_hkl.cell(),
				   wrk_hkl.resolution() );
  clipper::Xmap<float> wrk_map( wrk_hkl.spacegroup(), wrk_hkl.cell(),
				wrk_grid );
  wrk_map.fft_from( wrk_fc );

  // output the map
  clipper::CCP4MAPfile mapfile;
  mapfile.open_write( "diff.map" );
  mapfile.export_xmap( wrk_map );
  mapfile.close_write();
}
