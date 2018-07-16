// Clipper app to perform structure factor calculation
/* Copyright 2003-2004 Kevin Cowtan & University of York all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-mmdb.h>
#include "ccp4-extras.h"


int main( int argc, char** argv )
{
  CCP4program prog( "csfcalc", "0.1", "$Date: 2004/06/01" );

  // defaults
  clipper::String ippdb = "NONE";
  clipper::String ipfile = "NONE";
  clipper::String ipcolfo = "NONE";
  clipper::String ipcolfree = "NONE";
  clipper::String opfile = "sfcalc.mtz";
  clipper::String opcol = "sfcalc";
  bool bulk = true;
  int freeflag = 0;
  int n_refln = 1000;
  int n_param = 20;
  int verbose = 0;

  // command input
  CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if ( args[arg] == "-pdbin" ) {
      if ( ++arg < args.size() ) ippdb = args[arg];
    } else if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipfile = args[arg];
    } else if ( args[arg] == "-colin-fo" ) {
      if ( ++arg < args.size() ) ipcolfo = args[arg];
    } else if ( args[arg] == "-colin-free" ) {
      if ( ++arg < args.size() ) ipcolfree = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opfile = args[arg];
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg];
    } else if ( args[arg] == "-free-flag" ) {
      if ( ++arg < args.size() ) freeflag = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-num-reflns" ) {
      if ( ++arg < args.size() ) n_refln = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-num-params" ) {
      if ( ++arg < args.size() ) n_param = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-no-bulk" ) {
      bulk = false;
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "Unrecognized:\t" << args[arg] << "\n";
      args.clear();
    }
  }
  if ( args.size() <= 1 ) {
    std::cout << "Usage: csfcalc\n\t-pdbin <filename>\n\t-mtzin <filename>\n\t-colin-fo <colpath>\n\t-colin-free <colpath>\n\t-mtzout <filename>\n\t-colout <colpath>\n\t-free-flag <free set>\n\t-num-reflns <reflns per spline param>\n\t-num-params <spline params>\n\t-no-bulk\nStructure factor calculation with bulk solvent correction.\n";
    exit(1);
  }

  // make data objects
  clipper::CCP4MTZfile mtzin, mtzout;
  clipper::MTZcrystal cxtl;
  clipper::HKL_info hkls;
  double bulkfrc, bulkscl;
  typedef clipper::HKL_data_base::HKL_reference_index HRI;

  // open file
  mtzin.open_read( ipfile );
  mtzin.import_hkl_info( hkls );
  mtzin.import_crystal( cxtl, ipcolfo );
  clipper::HKL_data<clipper::data32::F_sigF> fo( hkls, cxtl );
  clipper::HKL_data<clipper::data32::Flag> free( hkls, cxtl );
  mtzin.import_hkl_data( fo, ipcolfo );
  if ( ipcolfree != "NONE" ) mtzin.import_hkl_data( free, ipcolfree );
  if ( opcol[0] != '/' ) opcol = mtzin.assigned_paths()[0].notail()+"/"+opcol;
  mtzin.close_read();

  // atomic model
  clipper::MMDBManager mmdb;
  mmdb.SetFlag( MMDBF_AutoSerials | MMDBF_IgnoreDuplSeqNum );
  mmdb.ReadPDBASCII( (char*)ippdb.c_str() );

  // get a list of all the atoms
  clipper::mmdb::PPCAtom psel;
  int hndl, nsel;
  hndl = mmdb.NewSelection();
  mmdb.SelectAtoms( hndl, 0, 0, SKEY_NEW );
  mmdb.GetSelIndex( hndl, psel, nsel );
  clipper::MMDBAtom_list atoms( psel, nsel );
  mmdb.DeleteSelection( hndl );

  // calculate structure factors
  clipper::HKL_data<clipper::data32::F_phi> fc( hkls, cxtl );
  if ( bulk ) {
    clipper::SFcalc_obs_bulk<float> sfcb;
    sfcb( fc, fo, atoms );
    bulkfrc = sfcb.bulk_frac();
    bulkscl = sfcb.bulk_scale();
  } else {
    clipper::SFcalc_aniso_fft<float> sfc;
    sfc( fc, atoms );
    bulkfrc = bulkscl = 0.0;
  }

  // now do sigmaa calc
  clipper::HKL_data<clipper::data32::F_phi> fb( hkls, cxtl ), fd( hkls, cxtl );
  clipper::HKL_data<clipper::data32::Phi_fom> phiw( hkls, cxtl );
  clipper::HKL_data<clipper::data32::Flag> flag( hkls, cxtl );
  for ( HRI ih = flag.first(); !ih.last(); ih.next() )
    if ( !fo[ih].missing() && (free[ih].missing()||free[ih].flag()==freeflag) )
      flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
    else
      flag[ih].flag() = clipper::SFweight_spline<float>::NONE;

  // do sigmaa calc
  clipper::SFweight_spline<float> sfw( n_refln, n_param );
  bool fl = sfw( fb, fd, phiw, fo, fc, flag );

  // calc map
  clipper::Grid_sampling grid( hkls.spacegroup(), cxtl, hkls.resolution() );
  clipper::Xmap<float> xmap( hkls.spacegroup(), cxtl, grid );
  xmap.fft_from( fd );

  // write difference map
  clipper::CCP4MAPfile mapout;
  mapout.open_write( "reftest.map" );
  mapout.export_xmap( xmap );
  mapout.close_write();

  // now loop over atoms and calculate parameter gradients and curvatures
  const double radius = 2.5;
  clipper::Grid_range gd( xmap.cell(), xmap.grid_sampling(), radius );
  clipper::Coord_grid g0, g1;
  std::vector<clipper::AtomShapeFn::TYPE> params;
  params.push_back( clipper::AtomShapeFn::X );
  params.push_back( clipper::AtomShapeFn::Y );
  params.push_back( clipper::AtomShapeFn::Z );
  params.push_back( clipper::AtomShapeFn::Occ );
  params.push_back( clipper::AtomShapeFn::Uiso );
  std::vector<double>                   func( atoms.size() );
  std::vector<std::vector<double> >     grad( atoms.size() );
  std::vector<clipper::Matrix<double> > curv( atoms.size() );
  double f, rho;
  std::vector<double> g(5);
  clipper::Matrix<double> c(5,5);
  clipper::Xmap<float>::Map_reference_coord i0, iu, iv, iw;
  for ( int a = 0; a < atoms.size(); a++ ) {
    func[a] = 0.0;
    grad[a].resize( 5, 0.0 );
    curv[a].resize( 5, 5, 0.0 );
    clipper::AtomShapeFn sf( atoms[a].coord_orth(), atoms[a].element(),
			     atoms[a].u_iso(), atoms[a].occupancy() );
    sf.agarwal_params() = params;
    g0 = xmap.coord_map( atoms[a].coord_orth() ).coord_grid() + gd.min();
    g1 = xmap.coord_map( atoms[a].coord_orth() ).coord_grid() + gd.max();
    i0 = clipper::Xmap<float>::Map_reference_coord( xmap, g0 );
    for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
      for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
        for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() )
	  if ( (iw.coord_orth()-atoms[a].coord_orth()).lengthsq() <
	       radius*radius ) {
	    rho = xmap[iw];
	    sf.rho_curv( iw.coord_orth(), f, g, c );
	    func[a] += f * rho;
	    for ( int i = 0; i < 5; i++ )
	      grad[a][i] += g[i] * rho;
	    for ( int i = 0; i < 5; i++ )
	      for ( int j = 0; j < 5; j++ )
		curv[a](i,j) += c(i,j) * rho;
	  }
  }

  // now calculate shifts to parameters
  for ( int a = 0; a < atoms.size(); a++ ) {
    // not all cross terms available, so calc xyz and U shifts independently
    std::vector<double> gxyz(3), sxyz(3);
    clipper::Matrix<double> cxyz(3,3);
    for ( int i = 0; i < 3; i++ )
      gxyz[i] = grad[a][i];
    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )
	cxyz(i,j) = curv[a](i,j);
    sxyz = cxyz.solve(gxyz);
    double gu = grad[a][4];
    double cu = curv[a](4,4);
    double su = gu/cu;
    std::cout << "grad " << a << " [" << atoms[a].element() << "]\t" << gxyz[0] << "  \t" << gxyz[1] << "  \t" << gxyz[2] << "  \t" << gu << "\n";
    std::cout << "Atom " << a << " [" << atoms[a].element() << "]\t" << sxyz[0] << "  \t" << sxyz[1] << "  \t" << sxyz[2] << "  \t" << su << "\n";
  }
}
