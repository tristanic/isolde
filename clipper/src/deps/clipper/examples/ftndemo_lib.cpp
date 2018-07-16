#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

extern "C" {

/* Input to this function is spacegroup, cell, atom list, and a list
   of HKLs for which structure factors are to be calculated. Output is
   the F's and phi's for those HKLs. */
void sfcalc_( const int* fspgr, const float* fcell, const int* natom, const float* x, const float* y, const float* z, const float* occ, const float* b, const int* atno, const int* nref, const int* h, const int* k, const int* l, float* fc, float* phic )
{
  char atmnames[][4] = { "H" ,"He","Li","Be","B" ,"C" ,"N" ,"O" ,"F" ,"Ne","Na","Mg","Al","Si","P" ,"S" ,"Cl","Ar","K" ,"Ca","Sc","Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y" ,"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf" };

  // make crystal info
  clipper::Spacegroup spgr = clipper::Spacegroup( clipper::Spgr_descr( *fspgr ) );
  clipper::Cell cell = clipper::Cell( clipper::Cell_descr( fcell[0], fcell[1], fcell[2], fcell[3], fcell[4], fcell[5] ) );
  double slim = 0.0;
  for ( int i = 0; i < *nref; i++ ) slim =
    clipper::Util::max(slim,clipper::HKL(h[i],k[i],l[i]).invresolsq(cell));
  clipper::Resolution reso( 0.999/sqrt(slim) );

  std::cout << " Spacegroup " << spgr.symbol_hall() << " \t " << spgr.symbol_hm() << "\n";
  std::cout << cell.format() << "\n";
  std::cout << " Resolution " << reso.limit() << " A \n";

  // make atom list
  std::vector<clipper::Atom> atomvec;
  for ( int i = 0; i < *natom; i++ ) {
    clipper::Atom atm;
    atm.set_coord_orth( clipper::Coord_orth( x[i], y[i], z[i] ) );
    atm.set_occupancy( occ[i] );
    atm.set_u_iso( clipper::Util::b2u( b[i] ) );    
    atm.set_element( clipper::String( atmnames[ atno[i] ] ) );
    atomvec.push_back( atm );
  }
  clipper::Atom_list atoms( atomvec );

  // make reflection list
  clipper::HKL_info hkls( spgr, cell, reso, true );
  clipper::HKL_data<clipper::data32::F_phi> fphi( hkls );

  // do structure factor calculation
  clipper::SFcalc_iso_fft<float> sfc;
  sfc( fphi, atoms );

  // extract the results
  for ( int i = 0; i < *nref; i++ ) {
    clipper::HKL hkl( h[i], k[i], l[i] );
    clipper::data32::F_phi dat = fphi[hkl];
    if ( !dat.missing() ) {
      fc[i]   = dat.f();
      phic[i] = dat.phi();
    } else {
      fc[i] = phic[i] = 0.0;
    }
  }
}

}
