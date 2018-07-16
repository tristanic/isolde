// Clipper structure factor calculation demo
/* (C) 2002 Kevin Cowtan */
// This is more of a demo application than a serious version

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-mmdb.h>

#include <iostream>
#include <algorithm>


using namespace clipper;

void sfcalc( HKL_data<datatypes::F_phi<float> >& fphidata, const Atom_list& atoms )
{
  // prepare target map
  const HKL_info& hkls = fphidata.base_hkl_info();
  const Grid_sampling grid( hkls.spacegroup(), hkls.cell(), hkls.resolution() );
  Xmap<float> xmap( hkls.spacegroup(), hkls.cell(), grid );

  // work out how big a box we need to calc density over for each atom
  Grid_range gd( hkls.cell(), grid, 3.0 );
  Xmap<float>::Map_reference_coord i0, iu, iv, iw;
  // loop over atoms
  for ( int i = 0; i < atoms.size(); i++ ) if ( !atoms[i].is_null() ) {
    AtomShapeFn sf( atoms[i] );  // get atom shape fn
    Coord_frac uvw = atoms[i].coord_orth().coord_frac( hkls.cell() );
    Coord_grid g0 = uvw.coord_grid( grid ) + gd.min();
    Coord_grid g1 = uvw.coord_grid( grid ) + gd.max();
    i0 = Xmap<float>::Map_reference_coord( xmap, g0 );
    // sum all map contributions from this atoms
    for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
      for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
        for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() )
          xmap[iw] += sf.rho( iw.coord_orth() );
  }

  // now correct for multiplicity of points on special positions
  for ( Xmap<float>::Map_reference_index ix = xmap.first();
	!ix.last(); ix.next() )
    xmap[ix] *= xmap.multiplicity( ix.coord() );

  // calc structure factors from map by fft
  xmap.fft_to( fphidata );
}


int main( int argc, char** argv )
{
  // make all the data objects
  clipper::CSpacegroup spgr( "base spgr" );
  clipper::CCell cell( spgr, "base cell" );
  clipper::CResolution reso( cell, "base reso", clipper::Resolution(3.0) );
  clipper::CHKL_info hkls( reso, "hkls", true );
  clipper::CHKL_data<clipper::data32::F_phi> fphi1( hkls );
  clipper::CHKL_data<clipper::data32::F_phi> fphi2( hkls );
  clipper::CHKL_data<clipper::data32::F_phi> fphi3( hkls );
  clipper::CHKL_data<clipper::data32::F_phi> fphi4( hkls );
  clipper::CHKL_data<clipper::data32::F_phi> fphi5( hkls );
  clipper::CHKL_data<clipper::data32::F_phi> fphi6( hkls );
  clipper::CGrid_sampling grid( reso );
  clipper::CXmap<float> xmap( grid );

  // atomic model
  clipper::MMDBManager mmdb;
  mmdb.ReadPDBASCII( argv[1] );
  //mmdb.write_file( "out.cif", 1 );

  // get info from mmdb
  spgr.init( mmdb.spacegroup() );
  cell.init( mmdb.cell() );

  // make a second mmdb as a test
  clipper::MMDBManager mmdb2;
  mmdb2.set_spacegroup( spgr ); mmdb2.set_cell( cell );
  mmdb2.cell().debug();         mmdb2.spacegroup().debug();

  // calc Z's
  {
    clipper::AtomSF sf( "C", 0.25 );
    clipper::HKL_info::HKL_reference_index ih;
    for ( ih = hkls.first(); !ih.last(); ih.next() )
      if (ih.index()%100000 == 0) {
	std::cout << ih.hkl().format() << " " <<
	  sf.f_iso(ih.hkl().coord_reci_orth(cell).invresolsq()) << " " <<
	  sf.f_aniso(ih.hkl().coord_reci_orth(cell)) << " " << 
	  ih.hkl().coord_reci_orth(cell).invresolsq() << " " << 
	  sf.f_iso(ih.hkl().coord_reci_orth(cell).invresolsq()) /
	  sf.f_aniso(ih.hkl().coord_reci_orth(cell)) << "\n";
      }
  }


  // debug info
  spgr.clipper::Container::debug();
  spgr.clipper::Spacegroup::debug();
  cell.clipper::Cell::debug();
  hkls.clipper::HKL_info::debug();

  // get a list of all the atoms
  //clipper::DBAtom_selection atoms = mmdb.select_atoms_serial( 0, 0 );
  clipper::mmdb::PPCAtom psel;
  int hndl, nsel;
  hndl = mmdb.NewSelection();
  mmdb.SelectAtoms( hndl, 0, 0, ::mmdb::SKEY_NEW );
  mmdb.GetSelIndex( hndl, psel, nsel );
  clipper::MMDBAtom_list atoms( psel, nsel ); 
  mmdb.DeleteSelection( hndl );

  clipper::SFcalc_iso_sum<float>( fphi1, atoms );
  clipper::SFcalc_iso_fft<float>( fphi2, atoms, 2.5, 2.0, 0.25 );
  clipper::SFcalc_aniso_sum<float>( fphi3, atoms );
  clipper::SFcalc_aniso_fft<float>( fphi4, atoms, 2.5, 2.0, 0.25 );

  /*
  sfcalc( fphi5, atoms );
  for ( HKL_info::HKL_reference_index ih = hkls.first(); !ih.last(); ih.next() ) {
    std::cout << ih.hkl().format() << " " << fphi4[ih].f() << " " << fphi3[ih].f() << "\n";
    std::cout << "                        " << clipper::Util::rad2d(fphi4[ih].phi()) << " " << clipper::Util::rad2d(fphi3[ih].phi()) << "\n";
  }
  */

  /*
  clipper::CCP4MAPfile mapout;
  mapout.open_write( "out.map" );
  mapout.export_xmap( xmap );
  mapout.close_write();
  clipper::CCP4MTZfile mtzout;
  mtzout.open_write( "out.mtz" );
  mtzout.export_hkl_info( hkls );
  mtzout.export_hkl_data( fphi1, clipper::MTZdataset("dset",1.0), clipper::MTZcrystal("xtal","proj",cell), "[FC1 PHIC1]" );
  mtzout.export_hkl_data( fphi2, clipper::MTZdataset("dset",1.0), clipper::MTZcrystal("xtal","proj",cell), "[FC2 PHIC2]" );
  mtzout.export_hkl_data( fphi3, clipper::MTZdataset("dset",1.0), clipper::MTZcrystal("xtal","proj",cell), "[FC3 PHIC3]" );
  mtzout.export_hkl_data( fphi4, clipper::MTZdataset("dset",1.0), clipper::MTZcrystal("xtal","proj",cell), "[FC4 PHIC4]" );
  mtzout.close_write();
  */

  clipper::HKL_info::HKL_reference_index ih;
  clipper::HKL_data<clipper::data32::E_sigE> esigdata( hkls );
  clipper::AtomSF sf( atoms[0].element(), 0.0 );
  for ( ih = hkls.first(); !ih.last(); ih.next() ) {
    esigdata[ih].E() = fphi3[ih].f() / sf.f_iso( ih.invresolsq() );
    esigdata[ih].sigE() = 1.0;
  }

  std::vector<clipper::ftype> p( 7, 0.0 );
  clipper::BasisFn_aniso_gaussian basisfn;
  clipper::TargetFn_meanEnth<clipper::data32::E_sigE> targetfn( esigdata, 2.0 );
  clipper::ResolutionFn_nonlinear rfn( hkls, basisfn, targetfn, p );
  std::cout << basisfn.scale( rfn.params() ) << "\n" << basisfn.u_aniso_orth( rfn.params() ).format() << "\n";

  // check the answers

  for ( ih = hkls.first(); !ih.last(); ih.next() ) {
    std::cout << ih.hkl().format() << " " << fphi1[ih].f() << " " << fphi2[ih].f() << " " << fphi3[ih].f() << " " << fphi4[ih].f() << "\n";
    std::cout << "                        " << clipper::Util::rad2d(fphi1[ih].phi()) << " " << clipper::Util::rad2d(fphi2[ih].phi()) << " " << clipper::Util::rad2d(fphi3[ih].phi()) << " " << clipper::Util::rad2d(fphi4[ih].phi()) << "\n";
  }

  // now test ed calc
  
  clipper::Grid_sampling gsam( spgr, cell, reso, 2.0 );
  clipper::Grid_range gmap( clipper::Coord_grid( 1, 2, 3 ),
			  clipper::Coord_grid( 9, 8, 7 ) );
  clipper::Xmap<float> exmap( spgr, cell, gsam );
  clipper::NXmap<float> enxmap( cell, gsam, gmap );

  clipper::EDcalc_iso<float> ediso;
  clipper::EDcalc_aniso<float> edaniso;

  ediso( enxmap, atoms );

  ediso( exmap, atoms );
  exmap.fft_to( fphi5 );
  edaniso( exmap, atoms );
  exmap.fft_to( fphi6 );

  for ( ih = hkls.first(); !ih.last(); ih.next() )
    if ( std::abs(std::complex<float>(fphi5[ih]) - std::complex<float>(fphi2[ih]) ) > 5.0 ) std::cout << "err" << ih.hkl().format() << " " << fphi5[ih].f() << " " << fphi2[ih].f() << " " << fphi6[ih].f() << " " << fphi4[ih].f() << "\n";

  clipper::CCP4MAPfile iomap;
  iomap.open_write( "map.map" );
  iomap.export_xmap( exmap );
  iomap.close_write();
  clipper::EDcalc_mask<float> edmask;
  edmask( exmap, atoms );
  iomap.open_write( "mask.map" );
  iomap.export_xmap( exmap );
  iomap.close_write();

  /*
  for ( int i = 0; i < atoms.size(); i++ )
    if ( atoms[i].is_atom() )
      std::cout << atoms[i].u_aniso_orth().format() << "\n" << atoms[i].u_aniso_orth().sqrt().format() << "\n" << (atoms[i].u_aniso_orth().sqrt()*atoms[i].u_aniso_orth().sqrt()).format() << "\n\n" ;
  */

  // now test some selections
  /*
  clipper::DBAtom_selection s1 = mmdb.select_atoms( "16-17" );
  clipper::DBAtom_selection s2 = mmdb.select_atoms( "15-20/CA[C]" );
  clipper::DBAtom_selection s3 = s1 & s2;
  clipper::DBAtom_selection s4 = s1 | s2;
  clipper::DBAtom_selection s5 = s1 ^ s2;
  for ( int i = 0; i < s3.size(); i++ )
    std::cout << i << " " << s3[i].residue().seqnum() << "\t" << s3[i].type() << "\n";
  for ( int i = 0; i < s4.size(); i++ )
    std::cout << i << " " << s4[i].residue().seqnum() << "\t" << s4[i].type() << "\n";
  for ( int i = 0; i < s5.size(); i++ )
    std::cout << i << " " << s5[i].residue().seqnum() << "\t" << s5[i].type() << "\n";
  */

  // test AtomSF objects and gradients
  clipper::Coord_orth co( 1.0, 2.0, 3.0 );
  clipper::AtomSF sf1( "N", 0.25, 1.0 );
  clipper::AtomShapeFn sf2( co, "N", 0.25, 1.0 );
  clipper::AtomShapeFn sfx( co+clipper::Coord_orth(0.001,0.0,0.0), "N", 0.25, 1.0 );
  clipper::AtomShapeFn sfy( co+clipper::Coord_orth(0.0,0.001,0.0), "N", 0.25, 1.0 );
  clipper::AtomShapeFn sfz( co+clipper::Coord_orth(0.0,0.0,0.001), "N", 0.25, 1.0 );
  clipper::AtomShapeFn sfo( co, "N", 0.25, 1.001 );
  clipper::AtomShapeFn sfu( co, "N", 0.251, 1.0 );
  clipper::AtomShapeFn sfx2( co+clipper::Coord_orth(0.002,0.0,0.0), "N", 0.25, 1.0 );
  clipper::AtomShapeFn sfy2( co+clipper::Coord_orth(0.0,0.002,0.0), "N", 0.25, 1.0 );
  clipper::AtomShapeFn sfz2( co+clipper::Coord_orth(0.0,0.0,0.002), "N", 0.25, 1.0 );
  clipper::AtomShapeFn sfo2( co, "N", 0.25, 1.002 );
  clipper::AtomShapeFn sfu2( co, "N", 0.252, 1.0 );
  clipper::AtomShapeFn sfxy( co+clipper::Coord_orth(0.001,0.001,0.0), "N", 0.25, 1.0 );
  clipper::AtomShapeFn sfyz( co+clipper::Coord_orth(0.0,0.001,0.001), "N", 0.25, 1.0 );
  clipper::AtomShapeFn sfzx( co+clipper::Coord_orth(0.001,0.0,0.001), "N", 0.25, 1.0 );
  std::vector<clipper::AtomShapeFn::TYPE> params;
  params.push_back( clipper::AtomShapeFn::X );
  params.push_back( clipper::AtomShapeFn::Y );
  params.push_back( clipper::AtomShapeFn::Z );
  params.push_back( clipper::AtomShapeFn::Occ );
  params.push_back( clipper::AtomShapeFn::Uiso );
  sf2.agarwal_params() = params;
  for ( int i = 0; i < 50; i++ ) {
    clipper::Coord_orth c1 = co + clipper::Coord_orth( 0.1*double(i), 0.0, 0.0 );
    clipper::Coord_orth c2 = co + clipper::Coord_orth( 0.0, 0.07*double(i), 0.1*double(i) );
    std::cout << sf1.rho_iso( (c1-co).lengthsq() ) << "\t" << sf2.rho( c1 ) << "\t-\t" << sf1.rho_iso( (c2-co).lengthsq() ) << "\t" << sf2.rho( c2 ) << "\n";
    std::cout << sf1.f_iso( 0.01*double(i) ) << "\t" << sf2.f( 0.01*double(i) ) << "\n";
    double f;
    std::vector<double> g(6);
    Matrix<double> c(6,6);
    sf2.rho_curv( c2, f, g, c );
    std::cout << "G> " << (sfx.rho(c2)-sf2.rho(c2))/0.001 << ":" << g[0] << "\t" << (sfy.rho(c2)-sf2.rho(c2))/0.001 << ":" << g[1] << "\t" << (sfz.rho(c2)-sf2.rho(c2))/0.001 << ":" << g[2] << "\t" << (sfo.rho(c2)-sf2.rho(c2))/0.001 << ":" << g[3] << "\t" << (sfu.rho(c2)-sf2.rho(c2))/0.001 << ":" << g[4] << "\n";
    std::cout << "C> " << (sfx2.rho(c2)-2*sfx.rho(c2)+sf2.rho(c2))/0.000001 << ":" << c(0,0) << "\t" << (sfy2.rho(c2)-2*sfy.rho(c2)+sf2.rho(c2))/0.000001 << ":" << c(1,1) << "\t" << (sfz2.rho(c2)-2*sfz.rho(c2)+sf2.rho(c2))/0.000001 << ":" << c(2,2) << "\t" << (sfo2.rho(c2)-2*sfo.rho(c2)+sf2.rho(c2))/0.000001 << ":" << c(3,3) << "\t" << (sfu2.rho(c2)-2*sfu.rho(c2)+sf2.rho(c2))/0.000001 << ":" << c(4,4) << "\n";
    std::cout << "c> " << (sfxy.rho(c2)-sfx.rho(c2)-sfy.rho(c2)+sf2.rho(c2))/0.000001 << ":" << c(0,1) << "\t" << (sfyz.rho(c2)-sfy.rho(c2)-sfz.rho(c2)+sf2.rho(c2))/0.000001 << ":" << c(1,2) << "\t" << (sfzx.rho(c2)-sfz.rho(c2)-sfx.rho(c2)+sf2.rho(c2))/0.000001 << ":" << c(2,0) << "\n";
  }

  std::cout << "----------------------------------------\n";

  clipper::U_aniso_orth uani( 0.25, 0.3, 0.35, 0.1, 0.05, 0.0 );
  clipper::AtomSF sf3( "N", uani, 1.0 );
  clipper::AtomShapeFn sf4( co, "N", uani, 1.0 );
  for ( double x = -1.0; x <= 1.0; x+= 0.333333 )
    for ( double y = -1.0; y <= 1.0; y+= 0.333333 )
      for ( double z = -1.0; z <= 1.0; z+= 0.333333 ) {
	clipper::Coord_orth c(x,y,z);
	clipper::Coord_reci_orth r(x,y,z);
	std::cout << c.format() << sf3.rho_aniso(c) << "\t" << sf4.rho(co+c) << "\t" << sf3.f_aniso(r) << "\t" << sf4.f(r) << "\n";
      }

  std::cout << "----------------------------------------\n";

  // now test ordinal functions
  std::vector<clipper::ftype> resols;
  for ( ih = hkls.first(); !ih.last(); ih.next() )
    resols.push_back( ih.invresolsq() );
  clipper::Generic_ordinal ordinal;
  ordinal.init(resols);
  clipper::Generic_ordinal ordinv = ordinal;
  ordinv.invert();
  clipper::Resolution_ordinal resord;
  resord.init(hkls,1.0);
  std::sort( resols.begin(), resols.end() );
  for ( int i = 0; i < resols.size(); i++ )
    std::cout << i << " " << resols[i] << " " << ordinal.ordinal(resols[i]) << " " << resord.ordinal(resols[i]) << " " << ordinv.ordinal(ordinal.ordinal(resols[i])) << "\n";

  std::cout << "----------------------------------------\n";

  // now test the Ramachandran plot class
  clipper::Ramachandran rama;
  std::cout << "\nNonGlyPro\n";
  rama.init( clipper::Ramachandran::NonGlyPro );
  for ( int psi = 180; psi >= -180; psi -= 10 ) {
    std::cout << clipper::String(psi,4) << " " ;
    for ( int phi = -180; phi <= 180; phi += 10 ) {
      if      ( rama.favored( clipper::Util::d2rad(phi),
			      clipper::Util::d2rad(psi) ) ) std::cout << "#";
      else if ( rama.allowed( clipper::Util::d2rad(phi),
			      clipper::Util::d2rad(psi) ) ) std::cout << "+";
      else                                                  std::cout << "-";
    }
    std::cout << "\n";
  }
  std::cout << "\nPro\n";
  rama.init( clipper::Ramachandran::Pro );
  for ( int psi = 180; psi >= -180; psi -= 10 ) {
    std::cout << clipper::String(psi,4) << " " ;
    for ( int phi = -180; phi <= 180; phi += 10 ) {
      if      ( rama.favored( clipper::Util::d2rad(phi),
			      clipper::Util::d2rad(psi) ) ) std::cout << "#";
      else if ( rama.allowed( clipper::Util::d2rad(phi),
			      clipper::Util::d2rad(psi) ) ) std::cout << "+";
      else                                                  std::cout << "-";
    }
    std::cout << "\n";
  }
  std::cout << "\nGly\n";
  rama.init( clipper::Ramachandran::Gly );
  for ( int psi = 180; psi >= -180; psi -= 10 ) {
    std::cout << clipper::String(psi,4) << " " ;
    for ( int phi = -180; phi <= 180; phi += 10 ) {
      if      ( rama.favored( clipper::Util::d2rad(phi),
			      clipper::Util::d2rad(psi) ) ) std::cout << "#";
      else if ( rama.allowed( clipper::Util::d2rad(phi),
			      clipper::Util::d2rad(psi) ) ) std::cout << "+";
      else                                                  std::cout << "-";
    }
    std::cout << "\n";
  }
  std::cout << "\nAll\n";
  rama.init( clipper::Ramachandran::All );
  for ( int psi = 180; psi >= -180; psi -= 10 ) {
    std::cout << clipper::String(psi,4) << " " ;
    for ( int phi = -180; phi <= 180; phi += 10 ) {
      if      ( rama.favored( clipper::Util::d2rad(phi),
			      clipper::Util::d2rad(psi) ) ) std::cout << "#";
      else if ( rama.allowed( clipper::Util::d2rad(phi),
			      clipper::Util::d2rad(psi) ) ) std::cout << "+";
      else                                                  std::cout << "-";
    }
    std::cout << "\n";
  }
}
