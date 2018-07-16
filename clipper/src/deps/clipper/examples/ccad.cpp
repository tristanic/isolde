// Clipper-cad utility.
/* (C) 2000 Kevin Cowtan */
// This is more of a demo application than a serious version
/*

Documetation:

This will read several mtz's and combine them into one. It will also
flip reflections about reciprocal space, or generate an hkl list onto
which data from files can be imported. It will do some spacegroup
manipulations, but will not expand from a higher order symmetry to a
lower at this point. (Later).

The application is run by a very simple parser, with operations being
performed in the order commands are entered. There is almost no error
checking at this point.


COMMANDS:

cell <a> <b> <c> <alpha> <beta> <gamma>
  Sets the base crystal cell. If this is omitted, it
  will be taken from the first inputfile.

spacegroup <spacegroupname>
  Sets the base crystal spacegroup. If this is omitted, it
  will be taken from the first inputfile.

resolution <resol>
  Set the resolution limit for import or generation of HKLs.

makehkls
  Generates a list of HKLs to the current resolution limit. If this is
  omitted, the list will be read from the first inputfile.

inputfile <filename>
  Open an mtz file for reading.

import <columnpaths> <hkl_datatype>
  Import the mtz columns or column group given by <columnpaths> to a
  hkl_data of type <hkl_datatype>.
  <hkl_datatype> can be one of:
       I_sigI, F_sigF, F_sigF_ano, F_phi, Phi_fom, ABCD, Flag
  <columnpaths> is of the form:
       /<crystal>/<dataset>/[<column_label_1>,<column_label_2>,...]
  or
       /<crystal>/<dataset>/<groupname>
  Automatic identification of column group types will be implemented
  in the future.

move <oldpath> <newpath>
  More or rename any hkl_data, dataset or crystal

outputfile <filename>
  Write all the accumulated data to an output mtz file


EXAMPLES:

./ccad << eof
inputfile testfile.mtz
import /unnamed_crystal1/native/[FP,SIGFP] F_sigF
outputfile junk.mtz
eof


*/

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <stdlib.h>
#include <iostream>


using namespace clipper;
using namespace clipper::data32;



int main()
{
  Message::set_message_level(1);

  CSpacegroup spgr( "spgr" );
  CCell cell( spgr, "cell" );
  CResolution reso( cell, "reso" );
  CHKL_info rfl( reso, "ccad" );

  CCP4MTZfile mtzfile;

  bool have_file = false;
  bool make_hkls = false;

  String line = "";

  // read and parse input lines
  while ( std::getline( std::cin, line ), !std::cin.eof() ) {
    std::vector<String> tokens = line.split(" ");

    // set cell
    if ( tokens[0] == "cell" ) {
      cell.init( Cell( Cell_descr( tokens[1].f(), tokens[2].f(), tokens[3].f(),
				   tokens[4].f(), tokens[5].f(), tokens[6].f() ) ) );
    }

    // set spacegroup
    if ( tokens[0] == "spacegroup" ) {
      spgr.init( Spacegroup( Spgr_descr( tokens[1] ) ) );
    }

    // set resolution
    if ( tokens[0] == "resolution" ) {
      reso.init( Resolution( ftype(tokens[1].f()) ) );
    }

    // generate hkls
    if ( tokens[0] == "makehkls" ) {
      make_hkls = true;
    }

    // new input file
    if ( tokens[0] == "inputfile" ) {
      if ( have_file ) { 
	mtzfile.close_read();    // close old file
      } else {
	mtzfile.open_read( tokens[1] );   // open new file
	if ( spgr.is_null() ) spgr.init( mtzfile.spacegroup() );
	if ( cell.is_null() ) cell.init( mtzfile.cell() );
	if ( reso.is_null() ) reso.init( mtzfile.resolution() );
	mtzfile.import_hkl_info( rfl, make_hkls );    // read hkls
	have_file = true;
      }
    }

    // import a column list
    else if ( tokens[0] == "import" ) {
      if ( !have_file ) exit(1);

      // read the data into a known hkl_data type
      Container* hkldata;
      if      ( tokens[2] == "I_sigI" )	hkldata = new CHKL_data<I_sigI>( rfl );
      else if ( tokens[2] == "F_sigF" )	hkldata = new CHKL_data<F_sigF>( rfl );
      else if (tokens[2]=="F_sigF_ano") hkldata = new CHKL_data<F_sigF_ano>( rfl );
      else if ( tokens[2] == "F_phi"  ) hkldata = new CHKL_data<F_phi>( rfl );
      else if ( tokens[2] == "Phi_fom") hkldata = new CHKL_data<Phi_fom>( rfl );
      else if ( tokens[2] == "ABCD"   ) hkldata = new CHKL_data<ABCD>( rfl );
      else if ( tokens[2] == "Flag"   ) hkldata = new CHKL_data<Flag>( rfl );
      hkldata->set_destroyed_with_parent();
      mtzfile.import_chkl_data( *hkldata, tokens[1] );
    }

    // move stuff about the tree
    else if ( tokens[0] == "move" ) {
      rfl.find_path_ptr( tokens[1] )->move( tokens[2] );
    }

    // write the data
    else if ( tokens[0] == "outputfile" ) {
      if ( have_file ) mtzfile.close_read();    // close old file

      mtzfile.open_write( tokens[1] );
      mtzfile.export_hkl_info( rfl );
      // search over crystals, datasets, and lists
      for ( int i = 0; i < rfl.num_children(); i++ )
	for ( int j = 0; j < rfl.child(i).num_children(); j++ )
	  for ( int k = 0; k < rfl.child(i).child(j).num_children(); k++ )
	    mtzfile.export_chkl_data( rfl.child(i).child(j).child(k), rfl.child(i).child(j).child(k).path() );
      mtzfile.close_write();
    }

    for (int i=0; i<tokens.size(); i++) std::cout << i << ":" << tokens[i] << "\n";
    rfl.Container::debug();
  }

}
