#include <clipper/clipper.h>

#include <iostream>


using namespace clipper;



int main()
{
  int sg_num,i,n,h,k,l;

  HKL equ;
  std::cout << "Enter spacegroup: ";
  std::cin >> sg_num;
  Spacegroup mysym=Spacegroup(Spgr_descr(sg_num));

  Spacegroup mycp1(mysym);

  mysym.debug();
  //mycp1.debug();
  //mycp2->debug();

  for (h=-2;h<=2;h++)
    for (k=-2;k<=2;k++)
      for (l=-2;l<=2;l++) {
	HKL rfl(h,k,l);
	if (mysym.recip_asu(rfl)) {
	  // if reflection is in ASU, check no mates are unless they
	  // map to the reflection itself
	  for (i=0;i<mysym.num_symops();i++) {
	    equ=rfl.transform(mysym.symop(i));
	    if (equ!=rfl && mysym.recip_asu(equ))
	      std::cout << "Dup (" << h << " " << k << " " << l << ") " << i << "\n";
	  }
	  for (i=0;i<mysym.num_symops();i++) {
	    equ=-rfl.transform(mysym.symop(i));
	    if (equ!=rfl && mysym.recip_asu(equ))
	      std::cout << "Dup (" << h << " " << k << " " << l << ") " << i << "\n";
	  }
	} else {
	  // if reflection is outside ASU, check at least one mate is inside
	  n=0;
	  for (i=0;i<mysym.num_symops();i++)
	    if (mysym.recip_asu( rfl.transform(mysym.symop(i)))) n++;
	  for (i=0;i<mysym.num_symops();i++)
	    if (mysym.recip_asu(-rfl.transform(mysym.symop(i)))) n++;
	  if (n==0) std::cout << "Abs (" << h << " " << k << " " << l << ") " << i << "\n";
	}
      }

  HKL rfl(0,0,0);
  if (!mysym.recip_asu(rfl)) std::cout << "Error (0,0,0)\n";

  for (i=0;i<9;i++) std::cout << "1 0 " << i << " " <<
		 Util::rad2d(HKL_class(mysym, HKL(1,0,i)).allowed()) << "\n";

  // test map stuff
  Grid_sampling c(8,8,8);
  //Grid_sampling c(16,16,16);

  std::cout << mysym.descr().symbol_hm() << "\n";
  std::cout << mysym.asu_min().u() << " " << mysym.asu_min().v() << " " << mysym.asu_min().w() << "\n";
  std::cout << mysym.asu_max().u() << " " << mysym.asu_max().v() << " " << mysym.asu_max().w() << "\n";

  Xmap<float> mymap( Spacegroup(mysym.descr()), Cell(Cell_descr(30.,40.,50.,90.,90.,90.)), c );
  Xmap<float> mymap2( mymap );

  h=k=0;
  for ( int u = 0; u < c.nu(); u++ ) {
    std::cout << u << "\n";
    for ( int v = 0; v < c.nv(); v++ ) {
      std::cout.width(3); std::cout << v << " ";
      for ( int w = 0; w < c.nw(); w++ ) {
	if (mymap.index_of( mymap.to_map_unit( Coord_grid(u,v,w) ) ) != -1) {
	  h += mysym.num_symops() / mymap.multiplicity( Coord_grid(u,v,w) );
	  k++;
	  std::cout << mymap.multiplicity( Coord_grid(u,v,w) );
	} else {
	  std::cout << "0";
	}
      }
      std::cout << "\n";
    }
  }
  std::cout << c.size() << " " << h << " " << k << " " << mysym.num_symops() << "\n";

  for ( Xmap<float>::Map_reference_index it = mymap.first(); !it.last(); it.next() ) {
    Coord_grid g = it.coord().unit(c);
    mymap2[it] = g.w()+10*(g.v()+10*g.u());
  }
  std::cout << "N " << h << "\n";
  Coord_grid x;
  for ( x.u() = 0; x.u() < c.nu(); x.u()++ ) {
    std::cout << x.u() << "\n";
    for ( x.v() = 0; x.v() < c.nv(); x.v()++ ) {
      std::cout.width(3); std::cout << x.v() << " ";
      for ( x.w() = 0; x.w() < c.nw(); x.w()++ ) {
	std::cout.width(4); std::cout << mymap2.get_data( x );
      }
      std::cout << "\n";
    }
  }

  /*
  for (i = 1; i < mysym.num_symops(); i++) {
    Matrix_symop_grid isymop( mysym.symop(i), Grid_sampling(48,48,48) );
    Coord_grid rot = isymop.inv_sym_coord(isymop.sym_coord(Coord_grid(3,5,7)));
    if (rot.u() != 3 || rot.v() != 5 || rot.w() != 7) std::cout << "Mat inv error";
  }
  */

  //c = Grid_sampling(24,32,40);
  c = Grid_sampling(24,24,24);
  Xmap<float> base( mysym, Cell(Cell_descr(30.,40.,50.,90.,90.,90.)), c );
  Xmap<float> test( mymap );
  Spacegroup p1(Spgr_descr("P1"));
  Xmap<float> test2( p1, Cell(Cell_descr(30.,40.,50.,90.,90.,90.)), c );

  for (i = 0; i<3; i++) {

  for (Xmap<float>::Map_reference_coord it( mymap ); !it.last(); it.next())
    base[it] = it.index();    

  // now test boundary checks
  Xmap<float>::Map_reference_coord jt( mymap );
  //  for (Xmap<float>::iterator it = test.first(); !it.last(); it.next()) {
  //  jt.set_coord( it.coord() );
  for (Xmap<float>::Map_reference_coord it( mymap ); !it.last(); it.next()) {
    jt = it;
    jt.next_u().next_v().next_w();
    test[it] = base[it] + base[jt];
  }

  for (Xmap<float>::Map_reference_coord it( mymap ); !it.last(); it.next()) {
    jt = it;
    jt.next_u().next_v().next_w();
    test[it] -=  base[jt];
    if (fabs(test[it]-base[it]) > 0.0001) std::cout << "Err "+it.coord().format()+"\n"; 
  }

  Xmap<float>::Map_reference_coord kt( test2 );
  for (Xmap<float>::Map_reference_coord it( test2 ); !it.last(); it.next()) {
    kt = it;
    kt.next_u();
    float r;
    r = base.get_data(it.coord());
    test2[kt] = r;
  }

  for (Xmap<float>::Map_reference_coord it( mymap ); !it.last(); it.next()) {
    kt = it;
    kt.next_v().next_w();
    test[it] -= test2[kt];
    if (fabs(test[it]-base[it]) > 0.0001) std::cout << "Err "+it.coord().format()+"\n"; 
  }

  }

  // cache test
  for (Xmap<float>::Map_reference_index iy=base.first(); !iy.last(); iy.next())
    base[iy] = iy.index();
  Xmap_base::Map_reference_coord ix( base );
  for ( i = 0; i < 1000000; i++ )
    ix.next_u();
}
