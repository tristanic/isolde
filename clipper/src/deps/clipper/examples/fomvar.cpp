#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>



using namespace clipper;
using namespace clipper::data32;


int main(int argc, char** argv)
{
  CCP4MTZfile file;
  clipper::HKL_data_base::HKL_reference_index ih;

  double mfom[40][20], nfom[40][20];
  for ( int i = 0; i < 40; i++ ) for ( int j = 0; j < 20; j++ )
    mfom[i][j] = nfom[i][j] = 0.0;

  HKL_info hkls;

  file.open_read(argv[1]);
  file.import_hkl_info( hkls, true );
  file.close_read();
  Resolution_ordinal resord;
  resord.init( hkls, 1.0 );

  HKL_data<ABCD> abcd(hkls);
  HKL_data<Phi_fom> fom(hkls);
  HKL_data<Phi_fom> fom1(hkls);
  HKL_data<Phi_fom> fom2(hkls);
  for ( ih = hkls.first(); !ih.last(); ih.next() ) {
    fom1[ih].phi() = fom2[ih].phi() = 0.0;
    fom1[ih].fom() = fom2[ih].fom() = 0.0;
  }

  double n = 0.0;
  for ( int i = 1; i < argc; i++ ) {
    file.open_read( argv[i] );
    file.import_hkl_data( abcd, "/*/*/sigmaa" );
    file.close_read();
    fom.compute( abcd, clipper::data32::Compute_phifom_from_abcd() );
    for ( ih = hkls.first(); !ih.last(); ih.next() ) 
      if ( !fom[ih].missing() ) {
	fom1[ih].fom() += fom[ih].fom();
	fom2[ih].fom() += fom[ih].fom()*fom[ih].fom();
	int bin = int( 39.999*resord.ordinal( ih.invresolsq() ) );
	mfom[bin][i] += fom[ih].fom();
	nfom[bin][i] += 1.0;
      }
    n += 1.0;
  }

  double s0 = 0.0, s1 = 0.0;
  for ( ih = hkls.first(); !ih.last(); ih.next() )
    if ( !fom1[ih].missing() ) {
      fom1[ih].fom() /= n;
      fom2[ih].fom() /= n;
      fom2[ih].fom() = sqrt( Util::max( fom2[ih].fom() - fom1[ih].fom()*fom1[ih].fom(), 0.0f ) );
      s0 += 1.0;
      s1 += fom2[ih].fom();
    }

  std::cout << s1/s0 << "\n";

  resord.invert();
  for ( int i = 0; i < 40; i++ ) {
    std::cout << resord.ordinal(double(i)/40.0) << "  \t";
    for ( int j = 0; j < 20; j++ )
      if ( nfom[i][j] > 0.0 )
	std::cout << mfom[i][j]/nfom[i][j] << "  \t";
    std::cout << "\n";
  }
}
