#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
//#include <clipper/clipper-ccp4.h>
#include <iostream>



using namespace clipper;
using namespace clipper::data32;


int main()
{
  CCP4MTZfile file;

  // import an mtz
  CSpacegroup myspgr( "" );
  CCell mycell( myspgr, "" );
  CResolution myreso( mycell, "" );
  CHKL_info mydata( myreso, "GerE native and MAD.." );
  CHKL_data<F_sigF> myfsig( mydata );
  CHKL_data<F_sigF_ano> myanom( mydata );

  file.open_read("testfile.mtz");
  myspgr.init( file.spacegroup() );
  mycell.init( file.cell() );
  myreso.init( Resolution( 10.0 ) );

  file.import_hkl_info( mydata, false );
  //myfsig.HKL_data<F_sigF>::debug();
  //myanom.HKL_data<F_sigF_ano>::debug();
  new CMTZcrystal(mydata, "xtal", MTZcrystal("xtal", "proj", mycell ));
  new CMTZdataset(mydata, "xtal/dset", MTZdataset("dset", 1.76) );
  file.import_chkl_data( myfsig, "*/native/[FP SIGFP]" );
  file.import_chkl_data( myanom, "*/*/[F(+)SEinfl SIGF(+)SEinfl F(-)SEinfl SIGF(-)SEinfl]", "xtal/dset/myanom");
  std::vector<String> v1 = file.column_paths();   for ( int i = 0; i < v1.size(); i++ ) std::cout << "Column: " << i << " " << v1[i] << "\n";
  std::vector<String> v2 = file.assigned_paths(); for ( int i = 0; i < v2.size(); i++ ) std::cout << "Import: " << i << " " << v2[i] << "\n";
  file.close_read();

  new Container(mydata, "newproj");
  new CMTZcrystal(mydata, "newproj/newcryst", MTZcrystal("newcryst", "newproj", Cell(Cell_descr(10.0,20.0,30.0))));
  new CMTZdataset(mydata, "newproj/newcryst/newdset", MTZdataset("newdset", 1.76) );

  mydata.Container::debug();
  std::cout << mycell.format() << "\n";

  for (int i=0; i<mydata.num_reflections(); i+=100) {
    std::cout << i << " " << myfsig[i].f() << " " << myfsig[i].sigf() << " " << myanom[i].f_pl() << " " << myanom[i].f_mi() << " " << myanom[i].cov() << "\n";
  }

  for (HKL_info::HKL_reference_index i = myfsig.first_data(); !i.last(); myfsig.next_data( i )) {
    std::cout << i.index() << " " << myfsig[i].f() << " " << myfsig[i].sigf() << "\n";
  }

  file.open_append("testfile.mtz","junk.mtz");
  file.export_chkl_data( myanom, "Gere:native/native/anom");
  file.close_append();
  file.open_write("junk2.mtz");
  file.export_hkl_info( mydata );
  file.export_chkl_data( myanom, "Gere:native/native/anom");
  file.close_write();

  // now test NaN
  std::cout << "\n----------------------------------------\n";

  std::cout << "NaN: " << sizeof(clipper::ftype32) << "=" << sizeof(clipper::itype32) << "\n";
  std::cout << "NaN: " << sizeof(clipper::ftype64) << "=" << sizeof(clipper::itype64) << "\n";
  std::cout << "Nan: ";
  {float x; Util::set_null(x);std::cout << x << Util::isnan(x) << Util::is_nan(x) << Util::is_null(x);}
  {double x; Util::set_null(x);std::cout << x << Util::isnan(x) << Util::is_nan(x) << Util::is_null(x);}
  {float x = Util::nanf();std::cout << x << Util::isnan(x) << Util::is_nan(x) << Util::is_null(x);}
  {double x = Util::nand();std::cout << x << Util::isnan(x) << Util::is_nan(x) << Util::is_null(x);}
  {float x = Util::nand();std::cout << x << Util::isnan(x) << Util::is_nan(x) << Util::is_null(x);}
  {double x = Util::nanf();std::cout << x << Util::isnan(x) << Util::is_nan(x) << Util::is_null(x);}
  std::cout << "\n";

  // now test the resolution function evaluator
  std::cout << "\n----------------------------------------\n";
  CHKL_info rfl( mycell, "GerE native and MAD..(2)" );
  CHKL_data<F_sigF> fsigdata( rfl, "" );

  // read data to higher resolution
  file.open_read("testfile.mtz");
  rfl.init( file.spacegroup(), file.cell(), Resolution(3.0) );
  file.import_hkl_info( rfl, false );
  file.import_chkl_data( fsigdata, "*/native/[FP SIGFP]" );
  file.close_read();

  Range<ftype> slim = fsigdata.invresolsq_range();
  std::cout << slim.min() << " " << slim.max() << "\n";

  TargetFn_meanFnth<F_sigF> targetfn( fsigdata, 2.0 );
  BasisFn_binner basisfn( fsigdata, 10 );
  std::vector<ftype> p1( 10, 1.0 );
  ResolutionFn rfn1( rfl, basisfn, targetfn, p1 );
  BasisFn_gaussian basisfn2;
  std::vector<ftype> q( 2, 1.0 );
  ResolutionFn_nonlinear rfn2( rfl, basisfn2, targetfn, q );
  BasisFn_expcubic basisfn3;
  std::vector<ftype> q1( 4, 1.0 );
  ResolutionFn_nonlinear rfn3( rfl, basisfn3, targetfn, q1 );
  BasisFn_spline basisfn4( fsigdata, 10 );
  std::vector<ftype> p2( 10, 1.0 );
  ResolutionFn rfn4( rfl, basisfn4, targetfn, p2 );

  CHKL_data<E_sigE> esigdata( fsigdata );
  esigdata.compute( fsigdata, Compute_EsigE_from_FsigF() );
  TargetFn_meanEnth<E_sigE> targetfn_e( esigdata, 2.0 );
  BasisFn_spline basisfn5( fsigdata, 10 );
  ResolutionFn rfn5( rfl, basisfn5, targetfn_e, p2 );  

  for ( HKL_info::HKL_reference_index ih = rfl.first(); !ih.last(); ih.next() ) {
    ftype eps = ih.hkl_class().epsilon();
    if ( fabs( rfn4.f(ih) - rfn5.f(ih) ) > 0.001 ) std::cout << "err\n"; 
    if (ih.index()%100 == 0) std::cout << ih.invresolsq() << " " << ih.hkl().format() << " " << rfn1.f(ih)*eps << " " << rfn2.f(ih)*eps << " " << rfn3.f(ih)*eps << " " << rfn4.f(ih)*eps << "\n";
  }
}
