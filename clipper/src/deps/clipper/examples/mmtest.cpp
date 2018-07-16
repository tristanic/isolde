// Clipper minimol test & demo
/* (C) 2003 Kevin Cowtan */
// This is more of a demo application than a serious version

#include <clipper/clipper.h>
#include <clipper/clipper-minimol.h>


using namespace clipper;

class thing : std::string
{
public:
  thing() { std::cout << "[CtorNull:\n"; }
  thing( const std::string& s ) : std::string(s) { std::cout << "[CtorInit:" << *this << "\n"; }
  thing( const thing& t ) : std::string(t) { std::cout << "[CtorCopy:" << *this << "\n"; }
  ~thing() { std::cout << "]Dtor:" << *this << "\n"; }
};



int main(int argc, char** argv)
{
  // first test the property manager

  PropertyManager p;
  {
    std::cout << "--make thing\n";
    thing it( std::string("goodbye") );
    std::cout << "--make prop\n";
    Property<thing> prop( it );
    std::cout << "--insert prop\n";
    p.set_property( "b", prop );
    std::cout << "--inserted prop\n";

    p.set_property( "a", Property<String>( "hello" ) );
    std::cout << p.exists_property( "a" ) << "\n";
    std::cout << dynamic_cast<const Property<String>& >(p.get_property( "a" )).value() << "\n";

    {
      std::cout << "--copy mgr1\n";
      PropertyManager p2;
      std::cout << "--copy mgr2\n";
      PropertyManager p1(p);
      std::cout << "--copy mgr3\n";
      p2 = p1;
      std::cout << "--copy mgr4\n";
      p2 = p;
      std::cout << "--copy mgr5\n";
    }

    std::cout << "-- vector check\n";
    std::vector<PropertyManager> v(3,p);
    std::cout << "-- done vector check\n";

    std::cout << "--delete prop\n";
    p.delete_property( "b" );
    std::cout << "--deleted prop\n";

  }
  std::cout << "--all gone!\n----------------------------------------\n\n";

  // Now test sequence align

  std::string s1 = "THEQUICKBOWNFOXJUMPSOVETHELAZY";
  std::string s2 = "QUICKBRWNFXJUMPSVERTH3LAZYDOG";
  std::string s3 = "ABCDEFTHELAZYDOGHIJKLMNOPQRSTUVWXYZ";
  clipper::MSequenceAlign align, lalign(clipper::MSequenceAlign::LOCAL);
  std::pair<std::vector<int>,std::vector<int> > result;
  result = align( s1, s2 );
  std::cout << s1 << std::endl;
  for ( int i = 0; i < s1.size(); i++ )
    std::cout << s1[i] << " " << (result.first[i]>=0 ? s2[result.first[i]] : ' ') << " " << i << " " << result.first[i] << std::endl;
  std::cout << s2 << std::endl;
  for ( int i = 0; i < s2.size(); i++ )
    std::cout << s2[i] << " " << (result.second[i]>=0 ? s1[result.second[i]] : ' ') << " " << i << " " << result.second[i] << std::endl;

  s1 = s2; s2 = s3;
  result = align( s1, s2 );
  std::cout << s1 << std::endl;
  for ( int i = 0; i < s1.size(); i++ )
    std::cout << s1[i] << " " << (result.first[i]>=0 ? s2[result.first[i]] : ' ') << " " << i << " " << result.first[i] << std::endl;
  std::cout << s2 << std::endl;
  for ( int i = 0; i < s2.size(); i++ )
    std::cout << s2[i] << " " << (result.second[i]>=0 ? s1[result.second[i]] : ' ') << " " << i << " " << result.second[i] << std::endl;

  result = lalign( s1, s2 );
  std::cout << s1 << std::endl;
  for ( int i = 0; i < s1.size(); i++ )
    std::cout << s1[i] << " " << (result.first[i]>=0 ? s2[result.first[i]] : ' ') << " " << i << " " << result.first[i] << std::endl;
  std::cout << s2 << std::endl;
  for ( int i = 0; i < s2.size(); i++ )
    std::cout << s2[i] << " " << (result.second[i]>=0 ? s1[result.second[i]] : ' ') << " " << i << " " << result.second[i] << std::endl;
  exit(1);

  // Now test minimol

  MiniMol tempmol;
  MMDBfile file;
  file.read_file( argv[1] );
  file.import_minimol( tempmol );
  std::cout << tempmol.spacegroup().symbol_hall() << "\n";
  std::cout << tempmol.cell().format() << "\n";

  MiniMol mol;
  mol.model() = tempmol.model();
  for ( int p = 0; p < mol.size(); p++ )
    for ( int m = 0; m < mol[p].size(); m++ )
      for ( int a = 0; a < mol[p][m].size(); a++ ) {
	std::cout << mol[p].id() << "\t" << mol[p][m].id() << "\t" << mol[p][m].type() << "\t<" << mol[p][m][a].id() << ">";
	std::cout << "\t" << dynamic_cast<const Property<String>&>(mol[p][m][a].get_property("CID")).value();
	if ( mol[p][m][a].exists_property("AltConf") ) std::cout << "\t" << dynamic_cast<const Property<String>&>(mol[p][m][a].get_property("AltConf")).value();
	std::cout << " ANISO:" << mol[p][m][a].u_aniso_orth().is_null();
	std::cout << "\n";
	mol[p][m][a].set_coord_orth( Coord_orth(1.0,2.0,3.0) );
      }

  file.export_minimol( mol );
  file.write_file( "mod.pdb" );

  MMDBfile file2;
  file2.export_minimol( mol );
  file2.write_file( "cpy.pdb" );

  // test non-bonding algorithm
  MAtomNonBond nb( tempmol, 5.0 );
  nb.debug();
  for ( double x = -30.0; x <  50.0; x += 10.0 )
    for ( double y = 320.0; y < 370.0; y += 10.0 )
      for ( double z =  30.0; z <  80.0; z += 10.0 ) {
	Coord_orth c(x,y,z);
	std::cout << c.format() << std::endl;
	std::vector<MAtomIndexSymmetry> atoms = nb(c,5.0);
	for ( int i = 0; i < atoms.size(); i++ ) {
	  Coord_orth d = tempmol[atoms[i].polymer()][atoms[i].monomer()][atoms[i].atom()].coord_orth();
	  Coord_frac f1 = c.coord_frac( tempmol.cell() );
	  Coord_frac f2 = d.coord_frac( tempmol.cell() );
	  double l1 = ( f2.symmetry_copy_near( tempmol.spacegroup(), tempmol.cell(), f1 ) - f1 ).lengthsq( tempmol.cell() );
	  f2 = tempmol.spacegroup().symop(atoms[i].symmetry()) * f2;
	  double l2 = ( f2.lattice_copy_near( f1 ) - f1 ).lengthsq( tempmol.cell() );
	  std::cout << i << " " << l1 << " " << l2 << " " << std::endl;
	}
      }



//   // Test selections and logical ops
//   std::cout << "-----------------------------------------------------\n";
//   MModel mol1 = mol.select( "*/13,14,15,16,17/N,CA,C" );
//   for ( int p = 0; p < mol1.size(); p++ )
//     for ( int m = 0; m < mol1[p].size(); m++ )
//       for ( int a = 0; a < mol1[p][m].size(); a++ )
// 	std::cout << mol1[p].id() << "\t" << mol1[p][m].id() << "\t" << mol1[p][m].type() << "\t" << mol1[p][m][a].id() << "\n";
//   std::cout << "-----------------------------------------------------\n";
//   MModel mol2 = mol.select( "*/15,16/*" );
//   for ( int p = 0; p < mol2.size(); p++ )
//     for ( int m = 0; m < mol2[p].size(); m++ )
//       for ( int a = 0; a < mol2[p][m].size(); a++ )
// 	std::cout << mol2[p].id() << "\t" << mol2[p][m].id() << "\t" << mol2[p][m].type() << "\t" << mol2[p][m][a].id() << "\n";
//   std::cout << "-----------------------------------------------------\n";
//   MModel mol3 = mol1 & mol2;
//   for ( int p = 0; p < mol3.size(); p++ )
//     for ( int m = 0; m < mol3[p].size(); m++ )
//       for ( int a = 0; a < mol3[p][m].size(); a++ )
// 	std::cout << mol3[p].id() << "\t" << mol3[p][m].id() << "\t" << mol3[p][m].type() << "\t" << mol3[p][m][a].id() << "\n";
//   std::cout << "-----------------------------------------------------\n";
//   MModel mol4 = mol1 | mol2;
//   for ( int p = 0; p < mol4.size(); p++ )
//     for ( int m = 0; m < mol4[p].size(); m++ )
//       for ( int a = 0; a < mol4[p][m].size(); a++ )
// 	std::cout << mol4[p].id() << "\t" << mol4[p][m].id() << "\t" << mol4[p][m].type() << "\t" << mol4[p][m][a].id() << "\n";

}
