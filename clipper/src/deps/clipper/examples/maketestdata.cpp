// Clipper app to make test data
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/core/test_data.h>
#include <clipper/core/test_core.h>
#include <clipper/contrib/test_contrib.h>

#include <fstream>
extern "C" {
#include <stdlib.h>
}


int main( int argc, char** argv )
{
  if ( argc == 1 ) {

    // do self tests
    std::cout << "Test core:\n";
    clipper::Test_core test_core;
    bool result_core = test_core();
    if ( result_core ) std::cout << "OK\n";
    std::cout << "Test contrib:\n";
    clipper::Test_contrib test_contrib;
    bool result_contrib = test_contrib();
    if ( result_contrib ) std::cout << "OK\n";
    // done self tests

    // error exit code
    if ( !( result_core && result_contrib ) ) exit(1);

  } else {

    // make data for self-tests

    // make data objects for reflection data
    clipper::CCP4MTZfile mtzin;
    clipper::HKL_info hkls;
    clipper::HKL_data<clipper::data32::F_sigF> fsig( hkls );
    clipper::HKL_data<clipper::data32::ABCD>   abcd( hkls );
    typedef clipper::HKL_data_base::HKL_reference_index HRI;

    // open file
    mtzin.open_read( argv[1] );
    clipper::Spacegroup spgr = mtzin.spacegroup();
    clipper::Cell       cell = mtzin.cell();
    clipper::Resolution reso( 5.0 );
    hkls.init( spgr, cell, reso, true );
    mtzin.import_hkl_data( fsig, "/*/*/[FNAT,SIGFNAT]" );
    mtzin.import_hkl_data( abcd, "/*/*/[HLA,HLB,HLC,HLD]" );
    mtzin.close_read();

    // print data
    for ( HRI ih = hkls.first(); !ih.last(); ih.next() ) {
      if ( !fsig[ih].missing() && !abcd[ih].missing() ) {
	int h = ih.hkl().h();
	int k = ih.hkl().k();
	int l = ih.hkl().l();
	double f = rint(fsig[ih].f()/0.1)*0.1;
	double sigf = rint(fsig[ih].sigf()/0.1)*0.1;
	double a = rint(abcd[ih].a()/0.01)*0.01;
	double b = rint(abcd[ih].b()/0.01)*0.01;
	double c = rint(abcd[ih].c()/0.01)*0.01;
	double d = rint(abcd[ih].d()/0.01)*0.01;
	std::cout << "{" << h << "," << k << "," << l << ","
		  << f << "," << sigf << ","
		  << a << "," << b << "," << c << "," << d << "},\n";
      }
    }

    // make data objects for model data
    clipper::MMDBManager mmdb;
    mmdb.ReadPDBASCII( argv[2] );
    // get a list of all the non-solvent atoms
    clipper::mmdb::PPCAtom psel;
    int hndl, nsel;
    hndl = mmdb.NewSelection();
    mmdb.Select( hndl, ::mmdb::STYPE_ATOM, "(!WAT,H2O,HOH)", ::mmdb::SKEY_NEW );
    mmdb.GetSelIndex( hndl, psel, nsel );
    clipper::MMDBAtom_list atoms( psel, nsel );
    mmdb.DeleteSelection( hndl );

    // print data
    for ( int i = 0; i < atoms.size(); i++ ) {
      atoms[i].set_element( atoms[i].element().trim() );
      double x = rint(atoms[i].coord_orth().x()/0.001)*0.001;
      double y = rint(atoms[i].coord_orth().y()/0.001)*0.001;
      double z = rint(atoms[i].coord_orth().z()/0.001)*0.001;
      double u_iso = rint(atoms[i].u_iso()/0.001)*0.001;
      double occ = rint(atoms[i].occupancy()/0.001)*0.001;
      std::cout << "{\"" << atoms[i].element() << "\","
		<< x << "," << y << "," << z << ","
		<< u_iso << "," << occ << "},\n";
    }

    // check the existing data
    clipper::data::Test_data data;
    for ( HRI ih = hkls.first(); !ih.last(); ih.next() ) {
      clipper::data32::F_sigF fs = data.hkl_data_f_sigf()[ih.hkl()];
      if ( !fs.missing() )
	if ( fabs( fs.f() - fsig[ih].f() ) > 0.05 )
	  std::cerr << "Err: " << ih.hkl().format() << "\n";
    }
    for ( int i = 0; i < data.atom_list().size(); i++ ) {
      if ( atoms[i].element() != data.atom_list()[i].element() ||
	   fabs(atoms[i].u_iso()    -data.atom_list()[i].u_iso()    ) > 0.001 ||
	   fabs(atoms[i].occupancy()-data.atom_list()[i].occupancy()) > 0.001 ||
	   ( atoms[i].coord_orth() -
	     data.atom_list()[i].coord_orth() ).lengthsq() > 0.001 ) {
	std::cerr << "Err: " << data.atom_list()[i].element() << "\t" <<
	  data.atom_list()[i].coord_orth().format() << "\n";
      }
    }

    // make tables for big result lists
    clipper::Test_contrib test_contrib;
    std::fstream contrib_data( "test_contrib.dat", std::fstream::out );
    test_contrib.set_stream( contrib_data );
    test_contrib();
    contrib_data.close();

    // done make test data

  }
}
