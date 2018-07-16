#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>
#include "xmap.cpp"



class IntEntry
{
public:
  IntEntry( const int& key ) : key_(key), str(key) {}
  bool matches( const int& key ) { return key == key_; }
  clipper::String format() const { return str; }
  clipper::String value() const { return str; }
private:
  int key_;
  clipper::String str;
};

clipper::ObjectCache<IntEntry> ocie;



int main()
{
  ocie.set_mode(clipper::ObjectCache<IntEntry>::NORMAL);
  clipper::data::SGcache.set_mode(clipper::ObjectCache<clipper::Spgr_cacheobj>::NORMAL);

  // test object cache
  int k1(1);
  int k2(2);
  int k3(3);
  {
    clipper::ObjectCache<IntEntry>::Reference ie1 = ocie.cache(k1);
    clipper::ObjectCache<IntEntry>::Reference ie2 = ocie.cache(k2);
    clipper::ObjectCache<IntEntry>::Reference ie3 = ocie.cache(k1);
    clipper::ObjectCache<IntEntry>::Reference ie4 = ie3;
    clipper::ObjectCache<IntEntry>::Reference ie5(ie3);
    ocie.debug();
    ie1 = ie2;
    ocie.debug();
  }
  ocie.debug();
  clipper::ObjectCache<IntEntry>::Reference ie1 = ocie.cache(k3);
  ocie.debug();

  std::cout << "----------------------------------------\n";

  const char* hallsymbols[] = {
// the 530 tabulated settings
 "P 1", "-P 1", "P 2y", "P 2ac 2ac", "P 2ac 2ab",
 "C 2c 2", "P 4c", "P 4n 2 -1n", "P 31", "-R 3", "-F 4vw 2vw 3" };

  clipper::data::SGcache.debug();
  clipper::Spacegroup sg3;
  clipper::Spacegroup sg4( clipper::Spacegroup::P1 );
  clipper::Spacegroup sg5( clipper::Spgr_descr( 19 ) );
  clipper::Spacegroup sg6( sg4 );
  clipper::Spacegroup sg7 = sg4;
  std::cout << "100 = " << sg3.is_null() << sg4.is_null() << sg5.is_null() << "\n";
  clipper::data::SGcache.debug();
  sg7 = sg3;
  std::cout << "100 = " << sg7.is_null() << sg4.is_null() << sg5.is_null() << "\n";
  clipper::data::SGcache.debug();

  for ( int s = 0; s < sizeof(hallsymbols)/sizeof(hallsymbols[0]); s++ ) {
    // get spacegroup
    clipper::String symbol = hallsymbols[s];
    clipper::Spacegroup sg( clipper::Spgr_descr( symbol, clipper::Spgr_descr::Hall ) );
    clipper::Spacegroup sg2( clipper::Spgr_descr( symbol, clipper::Spgr_descr::Hall ) );

    int errors = 0;
    for ( int i = 0; i < sg.num_symops(); i++ )
      if ( !sg.symop(i).equals( sg2.symop(i), 1.0e-6 ) ) errors++;

    std::cout << (errors==0?"OK      ":"Error   ") << sg.num_symops() << "\t" << sg.generator_ops().size() << "\t" << (sg.symbol_hm()+"        ").substr(0,12) << "\t: " << (sg.symbol_hall()+"          ").substr(0,12) << "\t" << sg.num_primops() << "\t" << errors << "\n";
    clipper::data::SGcache.debug();
  }
  clipper::data::SGcache.debug();

  std::cout << "----------------------------------------\n";

  clipper::ftype x0 = 10.0;
  for ( clipper::ftype x = x0; x < x0+2.0; x += 0.1 ) {
    clipper::Cell c1( clipper::Cell_descr( x0, x0, x0 ) );
    clipper::Cell c2( clipper::Cell_descr( x, x0, x0 ) );
    //clipper::Cell c3( clipper::Cell_descr( x, x, x ) );
    clipper::Cell c3( clipper::Cell_descr( x0, x0, x0, 90.0*x/x0 ) );
    std::cout << x << "\t" << c1.equals(c2) << c1.equals(c3) << c2.equals(c3) << "\n";
  }

  std::cout << "----------------------------------------\n";

  clipper::Cell cellc( clipper::Cell_descr( 37, 37, 37 ) );
  clipper::Cell cellha( clipper::Cell_descr( 37, 37, 37, 120, 90, 90 ) );
  clipper::Cell cellhb( clipper::Cell_descr( 37, 37, 37, 90, 120, 90 ) );
  clipper::Cell cellhc( clipper::Cell_descr( 37, 37, 37, 90, 90, 120 ) );
  clipper::Cell cellha1( clipper::Cell_descr( 37, 37, 37, 60, 90, 90 ) );
  clipper::Cell cellhb1( clipper::Cell_descr( 37, 37, 37, 90, 60, 90 ) );
  clipper::Cell cellhc1( clipper::Cell_descr( 37, 37, 37, 90, 90, 60 ) );
  clipper::Cell cell;
  for ( int s = 0; s < sizeof(hallsymbols)/sizeof(hallsymbols[0])-1; s++ ) {
    // get spacegroup
    clipper::String symbol = hallsymbols[s];
    clipper::Spacegroup sg( clipper::Spgr_descr( symbol, clipper::Spgr_descr::Hall ) );
    std::cout << "SG " << sg.symbol_hm() << "\n";
    // identify trigonal/hexagonal groups
    cell = cellc;
    for ( int sym = 1; sym < sg.num_symops(); sym++ ) {
      if ( ( sg.symop(sym).rot()(1,1) * sg.symop(sym).rot()(1,2) == -1 ) ||
           ( sg.symop(sym).rot()(2,1) * sg.symop(sym).rot()(2,2) == -1 ) )
        cell = cellha;
      if ( ( sg.symop(sym).rot()(0,0) * sg.symop(sym).rot()(0,2) == -1 ) ||
           ( sg.symop(sym).rot()(2,0) * sg.symop(sym).rot()(2,2) == -1 ) )
        cell = cellhb;
      if ( ( sg.symop(sym).rot()(0,0) * sg.symop(sym).rot()(0,1) == -1 ) ||
           ( sg.symop(sym).rot()(1,0) * sg.symop(sym).rot()(1,1) == -1 ) )
        cell = cellhc;
      if ( ( sg.symop(sym).rot()(1,1) * sg.symop(sym).rot()(1,2) == 1 ) ||
           ( sg.symop(sym).rot()(2,1) * sg.symop(sym).rot()(2,2) == 1 ) )
        cell = cellha1;
      if ( ( sg.symop(sym).rot()(0,0) * sg.symop(sym).rot()(0,2) == 1 ) ||
           ( sg.symop(sym).rot()(2,0) * sg.symop(sym).rot()(2,2) == 1 ) )
        cell = cellhb1;
      if ( ( sg.symop(sym).rot()(0,0) * sg.symop(sym).rot()(0,1) == 1 ) ||
           ( sg.symop(sym).rot()(1,0) * sg.symop(sym).rot()(1,1) == 1 ) )
        cell = cellhc1;
    }
    // now test map
    clipper::Grid_sampling grid(6,6,6);
    clipper::Xmap<float> mymap( sg, cell, grid );
    for ( clipper::Xmap<float>::Map_reference_index it = mymap.first(); !it.last(); it.next() ) {
      clipper::Coord_grid g = it.coord().unit(grid);
      mymap[it] = g.w()+10*(g.v()+10*g.u());
    }
    clipper::Coord_grid x;
    for ( x.u() = 0; x.u() < grid.nu(); x.u()++ ) {
      std::cout << x.u() << "\n";
      for ( x.v() = 0; x.v() < grid.nv(); x.v()++ ) {
	std::cout.width(3); std::cout << x.v() << " ";
	for ( x.w() = 0; x.w() < grid.nw(); x.w()++ ) {
	  std::cout.width(4); std::cout << mymap.get_data( x );
	}
	std::cout << "\n";
      }
    }
  }

  //new clipper::MiniMol( clipper::Spacegroup(clipper::Spgr_descr(19)), clipper::Cell(clipper::Cell_descr(100,100,100)) );
}
