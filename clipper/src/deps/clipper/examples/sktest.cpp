#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

#include <algorithm>

int main(int argc, char** argv)
{
  clipper::Spacegroup spgr( clipper::Spgr_descr( 1 ) );
  clipper::Cell cell( clipper::Cell_descr( 10, 10, 10 ));
  clipper::Grid_sampling grid( 6, 6, 6 );
  clipper::Xmap<float> xmap( spgr, cell, grid );
  clipper::Xmap<int> xskl( spgr, cell, grid );
  std::vector<int> index;
  clipper::Xmap<float>::Map_reference_index ix;

  // first make a shuffled map of value 0....n
  for ( ix = xmap.first(); !ix.last(); ix.next() )
    index.push_back( ix.index() );
  std::random_shuffle( index.begin(), index.end() );
  for ( int i = 0; i < index.size(); i++ )
    xmap.set_data( index[i], float(i) + 0.5 );
  index.clear();

  // print the map
  for ( ix = xmap.first(); !ix.last(); ix.next() )
    std::cout << ix.coord().format() << " " << xmap[ix] << "\n";

  // set up the whole map to be skeletonised
  for ( ix = xskl.first(); !ix.last(); ix.next() ) xskl[ix] = 1;
  // and skeletonise it
  clipper::Skeleton_fast<int,float> sk2;
  sk2( xskl, xmap );
  //clipper::Skeleton_basic sk2(1);
  //sk2( xskl, xmap );

  clipper::Coord_grid u, v, w;
  for ( u = clipper::Coord_grid(0,0,0); u.u() < grid.nu(); u.u()++ ) {
    for ( v = u; v.v() < grid.nv(); v.v()++ ) {
      for ( w = v; w.w() < grid.nw(); w.w()++ )
	std::cout << clipper::String(xmap.get_data( w ),5) << " ";
      std::cout << "\n";
    }
    std::cout << "\n";
  }
  for ( u = clipper::Coord_grid(0,0,0); u.u() < grid.nu(); u.u()++ ) {
    for ( v = u; v.v() < grid.nv(); v.v()++ ) {
      for ( w = v; w.w() < grid.nw(); w.w()++ )
	std::cout << clipper::String(xskl.get_data(w),5) << " ";
      std::cout << "\n";
    }
    std::cout << "\n";
  }


  // test nan
  float f;
  clipper::Util::set_null(f);
  std::cout << f << " " << clipper::Util::is_null(f) << " " << clipper::Util::is_nan(f) << " " << clipper::Util::isnan(f) << "\n";

  // now test euler angles
  {
    clipper::Euler_ccp4 e1( 0.3, 0.6, 0.8 );
    clipper::Euler_ccp4 e2 = clipper::Rotation(e1).euler_ccp4();;
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << clipper::Rotation(e1).format() << "\n"; 
  }
  {
    clipper::Euler<clipper::Rotation::EulerZYZr> e1( 0.3, 0.6, 0.8 );
    clipper::Euler<clipper::Rotation::EulerZYZr> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }

  std::cout << "--------------------\n";
  {
    clipper::Euler< 0> e1( 0.3, 0.6, 0.8 );
    clipper::Euler< 0> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler< 1> e1( 0.3, 0.6, 0.8 );
    clipper::Euler< 1> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler< 2> e1( 0.3, 0.6, 0.8 );
    clipper::Euler< 2> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler< 3> e1( 0.3, 0.6, 0.8 );
    clipper::Euler< 3> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler< 4> e1( 0.3, 0.6, 0.8 );
    clipper::Euler< 4> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler< 5> e1( 0.3, 0.6, 0.8 );
    clipper::Euler< 5> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler< 6> e1( 0.3, 0.6, 0.8 );
    clipper::Euler< 6> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler< 7> e1( 0.3, 0.6, 0.8 );
    clipper::Euler< 7> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler< 8> e1( 0.3, 0.6, 0.8 );
    clipper::Euler< 8> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler< 9> e1( 0.3, 0.6, 0.8 );
    clipper::Euler< 9> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler<10> e1( 0.3, 0.6, 0.8 );
    clipper::Euler<10> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler<11> e1( 0.3, 0.6, 0.8 );
    clipper::Euler<11> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler<12> e1( 0.3, 0.6, 0.8 );
    clipper::Euler<12> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler<13> e1( 0.3, 0.6, 0.8 );
    clipper::Euler<13> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler<14> e1( 0.3, 0.6, 0.8 );
    clipper::Euler<14> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler<15> e1( 0.3, 0.6, 0.8 );
    clipper::Euler<15> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler<16> e1( 0.3, 0.6, 0.8 );
    clipper::Euler<16> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler<17> e1( 0.3, 0.6, 0.8 );
    clipper::Euler<17> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler<18> e1( 0.3, 0.6, 0.8 );
    clipper::Euler<18> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler<19> e1( 0.3, 0.6, 0.8 );
    clipper::Euler<19> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler<20> e1( 0.3, 0.6, 0.8 );
    clipper::Euler<20> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler<21> e1( 0.3, 0.6, 0.8 );
    clipper::Euler<21> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler<22> e1( 0.3, 0.6, 0.8 );
    clipper::Euler<22> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }
  {
    clipper::Euler<23> e1( 0.3, 0.6, 0.8 );
    clipper::Euler<23> e2( e1.rotation() );
    std::cout << e1.alpha() << " " << e2.alpha() << "\t" << e1.beta() << " " << e2.beta() << "\t" << e1.gamma() << " " << e2.gamma() << "\n"; 
    std::cout << e1.format() << "\n" << e1.rotation().format() << "\n"; 
  }

  /*
  std::cout << "--------------------\n";
  std::cout << clipper::Euler< 0>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler< 1>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler< 2>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler< 3>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler< 4>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler< 5>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler< 6>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler< 7>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler< 8>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler< 9>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler<10>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler<11>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler<12>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler<13>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler<14>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler<15>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler<16>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler<17>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler<18>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler<19>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler<20>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler<21>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler<22>( 0.3, 0.6, 0.8 ).format() << "\n";
  std::cout << clipper::Euler<23>( 0.3, 0.6, 0.8 ).format() << "\n";

  for ( double u = 0.0; u < 5.0; u+=0.01 ) {
    clipper::Rotation r( fmod(3*u,1.0)-0.5, fmod(7*u,1.0)-0.5, fmod(13*u,1.0)-0.5, fmod(19*u,1.0)-0.5 );
    r.norm();
    std::cout << clipper::Util::rad2d(r.abs_angle()) << "\t" << r.polar_ccp4().format() << "\n";
  }
  */
}
