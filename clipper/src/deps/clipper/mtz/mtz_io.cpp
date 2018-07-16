/* mtz_io.cpp: class file for reflection data  mtz importer */
//C Copyright (C) 2000-2004 Kevin Cowtan and University of York

#include "mtz_io.h"
//L   This code is distributed under the terms and conditions of the
//L   CCP4 Program Suite Licence Agreement as a CCP4 Library.
//L   A copy of the CCP4 licence can be obtained by writing to the
//L   CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.


namespace clipper {


char MTZ_type_registry::names[200][12];
char MTZ_type_registry::types[200][4];
ftype32 MTZ_type_registry::scales[200];
MTZ_type_registry mtz_type_registry;


MTZ_type_registry::MTZ_type_registry()
{
  for ( int j = 0; j < 200; j++ ) names[j][0] = '\0';
  add_type( "I", "J", 1.0 );
  add_type( "sigI", "Q", 1.0 );
  add_type( "I+", "G", 1.0 );
  add_type( "sigI+", "L", 1.0 );
  add_type( "I-", "G", 1.0 );
  add_type( "sigI-", "L", 1.0 );
  add_type( "covI+-", "C", 1.0 );
  add_type( "F", "F", 1.0 );
  add_type( "sigF", "Q", 1.0 );
  add_type( "F+", "G", 1.0 );
  add_type( "sigF+", "L", 1.0 );
  add_type( "F-", "G", 1.0 );
  add_type( "sigF-", "L", 1.0 );
  add_type( "covF+-", "C", 1.0 );
  add_type( "E", "F", 1.0 );
  add_type( "sigE", "Q", 1.0 );
  add_type( "E+", "G", 1.0 );
  add_type( "sigE+", "L", 1.0 );
  add_type( "E-", "G", 1.0 );
  add_type( "sigE-", "L", 1.0 );
  add_type( "covE+-", "C", 1.0 );
  add_type( "A", "A", 1.0 );
  add_type( "B", "A", 1.0 );
  add_type( "C", "A", 1.0 );
  add_type( "D", "A", 1.0 );
  add_type( "phi", "P", Util::rad2d(1.0) );
  add_type( "fom", "W", 1.0 );
  add_type( "flag", "I", 1.0 );
}


void MTZ_type_registry::add_type( const String& name, const String& type, const ftype32& scale )
{
  int i, j;
  for ( j = 0; j < 200; j++ )
    if ( names[j][0] == '\0' ) break;
  if ( j == 200 )
    Message::message( Message_fatal( "MTZ_type_registry: registry full." ) );
  for ( i = 0; i < Util::min( int(name.length()), 11 ); i++ )
    names[j][i] = name[i];
  names[j][i] = '\0';
  for ( i = 0; i < Util::min( int(type.length()), 3 ); i++ )
    types[j][i] = type[i];
  types[j][i] = '\0';
  scales[j] = scale;
}


String MTZ_type_registry::type( const String& name )
{
  int j;
  for ( j = 0; j < 200; j++ )
    if ( String( names[j] ) == name ) break;
  if ( j == 200 )
    Message::message( Message_fatal( "MTZ_type_registry: name not found." ) );
  return String( types[j] );
}


ftype32 MTZ_type_registry::scale( const String& name )
{
  int j;
  for ( j = 0; j < 200; j++ )
    if ( String( names[j] ) == name ) break;
  if ( j == 200 )
    Message::message( Message_fatal( "MTZ_type_registry: name not found." ) );
  return scales[j];
}


/*! Check whether a dummy input column name has been given to suppress
  input or output of some of the data.
  \param path The column path & name.
  \return true if path ends with 'nan' or similar. */ 
bool MTZfile::is_virtual_col( const String path ) const
{
  // test for virtual column name
  String name = path.tail();
  return ( name=="MNF" || name=="NAN" || name=="mnf" || name=="nan" );
}


/*! Turn a clipper path name into a list of MTZ column paths. */
const std::vector<String> MTZfile::mtz_assign( const String assign, const String type, const String ftype, const int f_size ) const
{
  // Translate aassignment to mtz column names
  // ... loop through each list in turn
  std::vector<String> mtz_names( f_size, "MNF" );
  // interpret list name in terms of columns
  if ( assign.find( "[" ) == String::npos ) {
    // name is a single path: add type extenders
    std::vector<String> list = ftype.split(" ");
    for (int i=0; i < list.size(); i++) mtz_names[i] = assign + "." + type + "." + list[i];
  } else {
    // name is a list of mtz columns: extract column names from list
    String pref = assign.substr( 0, assign.find("[") );
    String post = assign.substr( assign.find("[") + 1,
				 assign.find("]") - assign.find("[") - 1 );
    std::vector<String> list = post.split(", ");
    for (int i=0; i < list.size(); i++) mtz_names[i] = pref + list[i];
  }
  return mtz_names;
}


/*! Read spacegroup, cell and resolution from header. */
void MTZfile::fetch_header_info()
{
  // get spacegroup by decoding symops
  String symops;
  char symop[80];
  for ( int i = 0; i < mmtz_io::mmtz_get_num_symops( mtzin ); i++ ) {
    mmtz_io::mmtz_get_symop( mtzin, i, symop );
    symops += String( symop ) + ";";
  }
  spacegroup_.init( Spgr_descr( symops, Spgr_descr::Symops ) );

  // get cell
  //float cp[6];
  //mmtz_io::mmtz_get_cell( mtzin, cp );
  //Cell_descr cd( cp[0], cp[1], cp[2], cp[3], cp[4], cp[5] );
  mmtz_io::mmtz_column col;
  mmtz_io::mmtz_dataset set;
  mmtz_io::mmtz_crystal xtl;
  Cell_descr cd;
  for ( int i = 0; i < mmtz_io::mmtz_num_cols( mtzin ); i++ ) {
    mmtz_io::mmtz_get_column( mtzin, i, &col, &set, &xtl );
    if ( col.label[0] == 'H' && col.label[1] == '\0' ) {
      cd = Cell_descr( xtl.cell[0], xtl.cell[1], xtl.cell[2],
		       xtl.cell[3], xtl.cell[4], xtl.cell[5] );
      goto setcell;
    }
  }
  float cp[6];
  mmtz_io::mmtz_get_cell( mtzin, cp );
  cd = Cell_descr( cp[0], cp[1], cp[2], cp[3], cp[4], cp[5] );
 setcell:
  cell_.init( cd );

  // get resolution: use mtz header. FIXME?
  umtz_io::umtz_hdr* hdr;
  float r1, r2;
  for ( hdr = umtz_io::umtz_first_head( mtzin ); hdr < umtz_io::umtz_last_head( mtzin ); hdr++ )
    if ( umtz_io::umtz_keymatch( hdr->data, "RESO" ) ) {
      String s( hdr->data, 80 );
      r1 = s.split(" ")[1].f();
      r2 = s.split(" ")[2].f();
    }
  resolution_.init( 0.9999/sqrt( Util::max(r1,r2) ) );
}


/*! Constructing an MTZfile does nothing except flag the object as not
  attached to any file for either input or output */
MTZfile::MTZfile()
{
  mode = NONE;
}


/*! Close any files which were left open. This is particularly
  important since to access the MTZ file efficiently, data reads and
  writes are deferred until the file is closed. */
MTZfile::~MTZfile()
{
  switch ( mode ) {
  case READ:
    close_read(); break;
  case WRITE:
    close_write(); break;
  case APPEND:
    close_append(); break;
  }
}


/*! The file is opened for reading. This MTZfile object will remain
  attached to this file until it is closed. Until that occurs, no
  other file may be opened with this object, however another MTZfile
  object could be used to access another file.
  \param filename_in The input filename or pathname. */
void MTZfile::open_read( const String filename_in )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "MTZfile: open_read - File already open" ) );

  // read the mtz
  mtzin = mmtz_io::mmtz_open( filename_in.c_str(), "r" );
  if ( mtzin == NULL )
    Message::message( Message_fatal( "MTZfile: open_read  - Could not read: "+filename_in ) );
  fetch_header_info();

  mode = READ;
}


/*! Close the file after reading. This command also actually fills in
  the data in any HKL_data structures which have been marked for
  import. */
void MTZfile::close_read()
{
  if ( mode != READ )
    Message::message( Message_fatal( "MTZfile: no file open for read" ) );

  // close the mtz:
  // this actually imports the data from the mtz file
  // - we save up all the work to do it on a single pass

  HKL hkl;
  float fdata[1000]; // file data: allow up to 1000 columns on file
  int   idata[1000]; // flag data: allow up to 1000 columns on file
  xtype pdata[100];  // program data: allow up to 100 floats per group type
  int ref, lst, col;

  // update the data lists to ensure consistency
  for (lst=0; lst < list_map_i.size(); lst++) list_map_i[lst].list->update();

  // now import the data
  /*
    Since the mtz sorted with reflection index increasing most slowly,
    it must be read a reflection at a time.
  */
  mmtz_io::mmtz_seek_row( mtzin, 0 );
  for (ref=0; ref < mmtz_io::mmtz_num_rows( mtzin ); ref++) {
    // read reflection
    mmtz_io::mmtz_get_row( mtzin, fdata, idata );
    hkl = HKL( int(fdata[0]), int(fdata[1]), int(fdata[2]) );
    // do each data list
    for (lst=0; lst < list_map_i.size(); lst++) {
      // for each list, grab the relevent columns
      for (col=0; col < list_map_i[lst].col.size(); col++) {
	// set to NaN unless readable and present
	pdata[col] = Util::nan();
	// try and read the data
	if ( list_map_i[lst].col[col] != 0 )
	  if ( idata[ list_map_i[lst].col[col] ] )
	    pdata[col] =
	      xtype(fdata[list_map_i[lst].col[col]] / list_map_i[lst].scl[col]);
      }
      // now set the data
      list_map_i[lst].list->data_import( hkl, pdata );
    }
  }
  // tidy up
  list_map_i.clear();
  mmtz_io::mmtz_close( mtzin );

  mode = NONE;
}


/*! The file is opened for writing. This will be a new file, created
  entirely from data from within the program, rather than by extending
  an existing file. Similar restrictions apply as for open_read().

  In practice the open_append() method is usually preferred.
  \param filename_out The output filename or pathname. */
void MTZfile::open_write( const String filename_out )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "MTZfile: open_write - File already open" ) );

  // open the output mtz
  mtzout = mmtz_io::mmtz_open( filename_out.c_str(), "w" );
  if ( mtzout == NULL )
    Message::message( Message_fatal( "MTZfile: open_write - Could not write: "+filename_out ) );

  mode = WRITE;
}


/*! Close the file after writing. This command also actually writes
  the data reflection list from the HKL_info object and the data from
  any HKL_data objects which have been marked for import. */
void MTZfile::close_write()
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "MTZfile: no file open for write" ) );

  // export the marked list data to an mtz file
  HKL hkl;
  float fdata[1000]; // file data: allow up to 1000 columns on file
  int   idata[1000]; // flag data: allow up to 1000 columns on file
  xtype pdata[100];  // program data: allow up to 100 floats per group type
  int ref, lst, col;

  if ( hkl_ptr == NULL )
    Message::message( Message_fatal( "MTZfile: no refln list exported" ) );

  // write the data
  mmtz_io::mmtz_seek_row( mtzout, 0 );
  for (ref=0; ref < hkl_ptr->num_reflections(); ref++) {
    // read/append reflection
    hkl = hkl_ptr->hkl_of(ref);
    fdata[0] = float( hkl.h() );
    fdata[1] = float( hkl.k() );
    fdata[2] = float( hkl.l() );
    idata[0] = idata[1] = idata[2] = 1;
    // output reflection
    for (lst=0; lst < list_map_o.size(); lst++) {
      // for each list, fetch the data
      list_map_o[lst].list->data_export( hkl, pdata );
      // and copy to the relevent columns
      for (col=0; col < list_map_o[lst].col.size(); col++) {
	// set to mtz->mnf, unless we have a value for it
	if ( !Util::is_nan( pdata[col] ) ) {
	  idata[ list_map_o[lst].col[col] ] = 1;
	  fdata[ list_map_o[lst].col[col] ] =
	    float( pdata[col] * list_map_o[lst].scl[col] );
	} else {
	  idata[ list_map_o[lst].col[col] ] = 0;
	  fdata[ list_map_o[lst].col[col] ] = 0.0;
	}
      }
    }
    // write appended record
    mmtz_io::mmtz_add_row( mtzout, fdata, idata );
  }

  // tidy up
  list_map_o.clear();
  mmtz_io::mmtz_close( mtzout );

  mode = NONE;
}


/*! A file is opened for appending. One file is opened for reading,
  and a second is opened for writing. The second file will contain all
  the information from the first, plus any additional columns exported
  from HKL_data objects.
  \param filename_in The input filename or pathname.
  \param filename_out The output filename or pathname. */
void MTZfile::open_append( const String filename_in, const String filename_out )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "MTZfile: open_append - File already open" ) );

  // open the input and output mtzs
  mtzin  = mmtz_io::mmtz_open( filename_in.c_str(), "r" );
  if ( mtzin == NULL )
    Message::message( Message_fatal( "MTZfile: open_append  - Could not read: "+filename_in ) );
  mtzout = mmtz_io::mmtz_open( filename_out.c_str(), "w" );
  if ( mtzout == NULL )
    Message::message( Message_fatal( "MTZfile: open_append - Could not write: "+filename_out ) );
  fetch_header_info();

  // copy the headers across
  mmtz_io::mmtz_copy_headers( mtzout, mtzin );

  mode = APPEND;
}


/*! Close the files after appending. This command actually copies the
  input file to the output file, adding data from any HKL_data objects
  which have been marked for import. */
void MTZfile::close_append()
{
  if ( mode != APPEND )
    Message::message( Message_fatal( "MTZfile: no file open for append" ) );

  // export the marked list data to an mtz file
  HKL hkl;
  float fdata[1000]; // file data: allow up to 1000 columns on file
  int   idata[1000]; // flag data: allow up to 1000 columns on file
  xtype pdata[100];  // program data: allow up to 100 floats per group type
  int ref, lst, col;

  // copy the data
  mmtz_io::mmtz_seek_row( mtzin, 0 );
  for (ref=0; ref < mmtz_io::mmtz_num_rows( mtzin ); ref++) {
    // read/append reflection
    mmtz_io::mmtz_get_row( mtzin, fdata, idata );
    hkl = HKL( int(fdata[0]), int(fdata[1]), int(fdata[2]) );
    // output reflection
    for (lst=0; lst < list_map_o.size(); lst++) {
      // for each list, fetch the data
      list_map_o[lst].list->data_export( hkl, pdata );
      // and copy to the relevent columns
      for (col=0; col < list_map_o[lst].col.size(); col++) {
	// set to mtz->mnf, unless we have a value for it
	if ( !Util::is_nan( pdata[col] ) ) {
	  idata[ list_map_o[lst].col[col] ] = 1;
	  fdata[ list_map_o[lst].col[col] ] =
	    float( pdata[col] * list_map_o[lst].scl[col] );
	} else {
	  idata[ list_map_o[lst].col[col] ] = 0;
	  fdata[ list_map_o[lst].col[col] ] = 0.0;
	}
      }
    }
    // write appended record
    mmtz_io::mmtz_add_row( mtzout, fdata, idata );
  }

  // tidy up
  list_map_o.clear();
  mmtz_io::mmtz_close( mtzin );
  mmtz_io::mmtz_close( mtzout );

  mode = NONE;
}


/*! Get the spacegroup from the MTZ file. \return The spacegroup. */
const Spacegroup& MTZfile::spacegroup() const
{ return spacegroup_; }


/*! Get the base cell from the MTZ file. \return The cell. */
const Cell& MTZfile::cell() const
{ return cell_; }


/*! Get the resolution limit from the MTZ file. \return The resolution. */
const Resolution& MTZfile::resolution() const
{ return resolution_; }


/*! Import the list of reflection HKLs from an MTZ file into an
  HKL_info object. If the resolution limit of the HKL_info object is
  lower than the limit of the file, any excess reflections will be
  rejected, as will any systematic absences or duplicates.
  \param target The HKL_info object to be initialised. */
void MTZfile::import_hkl_list( HKL_info& target )
{
  if ( mode != READ )
    Message::message( Message_fatal( "MTZfile: no file open for read" ) );

  float fdata[1000]; // file data: allow up to 1000 columns on file
  int   idata[1000]; // flag data: allow up to 1000 columns on file
  HKL hkl;
  std::vector<HKL> hkls;

  // read the reflections from the mtz
  ftype slim = target.resolution().invresolsq_limit();
  // read the reflections
  mmtz_io::mmtz_seek_row( mtzin, 0 );
  for (int ref=0; ref < mmtz_io::mmtz_num_rows( mtzin ); ref++) {
    // read reflection
    mmtz_io::mmtz_get_row( mtzin, fdata, idata );
    hkl = HKL( int(fdata[0]), int(fdata[1]), int(fdata[2]) );
    // check the resolution against the master cell
    bool in_res = hkl.invresolsq(target.cell()) < slim;
    // also check against any hkl_datas which have been assigned
    //for (int lst=0; lst < list_map_i.size(); lst++) in_res = in_res ||
    //        (hkl.invresolsq(list_map_i[lst].list->base_cell()) < slim);
    // store reflection
    if ( in_res ) hkls.push_back( hkl );
  }
  target.add_hkl_list( hkls );
}


/*! Import a complete HKL_info object. The supplied HKL_info object is
  examined, and if any of the parameters (spacegroup, cell, or
  resolution) are unset, then they will be set using values from the
  file. The reflections list will then be generated (the default), or
  imported from the file.

  This method is a shortcut which can generally replace most common
  combinations of calls to import_spacegroup(), import_cell(),
  import_resolution() and import_hkl_list().

  \param target The HKL_info object to be initialised.
  \param generate Generate the list of HKLs rather than importing it from the file. */
void MTZfile::import_hkl_info( HKL_info& target, const bool generate )
{
  // first check if the HKL_info params are already set
  Spacegroup s = target.spacegroup();
  Cell       c = target.cell();
  Resolution r = target.resolution();
  // import any missing params
  if ( s.is_null() ) s = spacegroup_;
  if ( c.is_null() ) c = cell_;
  if ( r.is_null() ) r = resolution_;
  target.init( s, c, r );
  // now make the HKL list
  if ( generate )
    target.generate_hkl_list();
  else
    import_hkl_list( target );
}


/*! Import data from an MTZ file into an HKL_data object. The dataset
  and crystal information from the first corresponding MTZ column are
  also returned.

  An MTZ column type must be present in the MTZ_type_registry for the
  HKL_data type element name concerned.

  This routine does not actually read any data, but rather marks the
  data to be read when the file is closed.

  For container objects import_chkl_data() is preferred.
  \param cdata The HKL_data object into which data is to be imported. 
  \param cset The MTZdataset object to store the MTZ column wavelength info. 
  \param cxtl The Crystal object to store the MTZ column cell info. */
void MTZfile::import_hkl_data( HKL_data_base& cdata, MTZdataset& cset, MTZcrystal& cxtl, const String mtzpath )
{
  if ( mode != READ )
    Message::message( Message_fatal( "MTZfile: no file open for read" ) );

  // make sure the date list is sized
  cdata.update();

  // decode the mtz column info
  String coltype;
  String type  = cdata.type();
  int    ncols = cdata.data_size();
  String names = cdata.data_names();
  std::vector<String> col_names = mtz_assign(mtzpath.tail(),type,names,ncols);
  std::vector<String> dat_names = names.split(" ");

  // now look up the column, dataset and crystal
  column_map_i colmap;
  colmap.col.resize( ncols, 0 );  // mark all columns as null
  colmap.scl.resize( ncols, 1.0 );
  mmtz_io::mmtz_column  fcol;
  mmtz_io::mmtz_dataset fset;
  mmtz_io::mmtz_crystal fxtl;
  String setname = mtzpath.notail().tail();
  String xtlname = mtzpath.notail().notail().tail();
  for ( int mcol = 0; mcol < mmtz_io::mmtz_num_cols( mtzin ); mcol++ ) {
    mmtz_io::mmtz_get_column( mtzin, mcol, &fcol, &fset, &fxtl );
    if ( ( setname == "*" || setname == fset.dname ) &&
         ( xtlname == "*" || xtlname == fxtl.xname ) ) {
      // this col belongs to the right dset/xtal, check the name
      for ( int ccol = 0; ccol < ncols; ccol++ )
        if ( col_names[ccol] == fcol.label ) {
          colmap.col[ccol] = mcol; // store column index
	  colmap.scl[ccol] = MTZ_type_registry::scale( dat_names[ccol] );
	  coltype = MTZ_type_registry::type( dat_names[ccol] );
          if ( coltype[0] != fcol.type[0] )   // check type
	    Message::message( Message_warn( "MTZfile: Mtz column type mismatch: "+String(fcol.label)+" "+String(fcol.type,1)+"-"+coltype.substr(0,1) ) );
        }
    }
  }
  // check for non-virtual unassigned columns
  for ( int ccol = 0; ccol < ncols; ccol++ )
    if ( ( colmap.col[ccol] == 0 ) && !is_virtual_col( col_names[ccol] ) )
      Message::message( Message_fatal( "MTZfile: No such column: "+col_names[ccol] ) );
  // fetch crystal
  for ( int ccol = 0; ccol < ncols; ccol++ )
    if ( colmap.col[ccol] != 0 ) {
      mmtz_io::mmtz_get_column( mtzin, colmap.col[ccol], &fcol, &fset, &fxtl );
      break;
    }

  // copy crystal info
  cxtl = MTZcrystal( fxtl.xname, fxtl.pname,
		     Cell( Cell_descr( fxtl.cell[0], fxtl.cell[1],
				       fxtl.cell[2], fxtl.cell[3],
				       fxtl.cell[4], fxtl.cell[5] ) ) );
  cset = MTZdataset( fset.dname, fset.wavel );
  // store the list<->mtzcol mappings
  colmap.list = &cdata;
  list_map_i.push_back( colmap );
}


/*! Export a complete HKL_info object, including spacegroup, cell, and
  list of reflection HKLs from an HKL_info object to an MTZ file. This
  is compulsory when writing an MTZ file, but forbidden when
  appending, since the HKLs will then come from the input MTZ.
  \param target The HKL_info object to supply the parameters. */
void MTZfile::export_hkl_info( const HKL_info& target )
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "MTZfile: no file open for write" ) );

  // write initial mtz headers
  Cell_descr cp = target.cell().descr();
  float cell[6] = { cp.a(), cp.b(), cp.c(), cp.alpha_deg(), cp.beta_deg(), cp.gamma_deg() };
  mmtz_io::mmtz_init_headers( mtzout, "Created by clipper", cell );

  // write sort order - need to check sort order first
  // mmtz_add_sort_header( mtzout, "HKL" );

  // write spacegroup header
  const Spacegroup& sg = target.spacegroup();
  int   sgnum = sg.descr().spacegroup_number();
  String name = sg.descr().symbol_hm();
  String laue = sg.symbol_laue();
  mmtz_io::mmtz_add_syminf_header( mtzout, sg.num_symops(), sg.num_primops(), name[0], sgnum, name.c_str(), laue.c_str() );

  // write symops
  for ( int i = 0; i < sg.num_symops(); i++ )
    mmtz_io::mmtz_add_symop_header( mtzout, sg.symop(i).format().c_str() );

  // add the h,k,l columns
  mmtz_io::mmtz_column  fcol;
  mmtz_io::mmtz_dataset fset;
  mmtz_io::mmtz_crystal fxtl;
  strcpy( fxtl.xname, "HKL" );
  strcpy( fxtl.pname, "HKL" );
  fxtl.cell[0] = cp.a(); fxtl.cell[3] = cp.alpha_deg();
  fxtl.cell[1] = cp.b(); fxtl.cell[4] = cp.beta_deg();
  fxtl.cell[2] = cp.c(); fxtl.cell[5] = cp.gamma_deg();
  strcpy( fset.dname, "HKL" );
  fset.wavel = 9.999;
  strcpy( fcol.type, "H" );
  strcpy( fcol.label, "H" );
  mmtz_io::mmtz_add_column( mtzout, &fcol, &fset, &fxtl );
  strcpy( fcol.label, "K" );
  mmtz_io::mmtz_add_column( mtzout, &fcol, &fset, &fxtl );
  strcpy( fcol.label, "L" );
  mmtz_io::mmtz_add_column( mtzout, &fcol, &fset, &fxtl );

  // save a pointer to the parent hkl list which is supplying the hkls
  hkl_ptr = &target;
}


/*! Export data from an HKL_data object into an MTZ file. MTZdataset and
  crystal information must be supplied, and will be applied to all
  columns in the output MTZ.

  An MTZ column type must be present in the MTZ_type_registry for the
  HKL_data type element name concerned.

  This routine does not actually write any data, but rather marks the
  data to be written when the file is closed.

  Normally export_chkl_data() is preferred.
  \param cdata The HKL_data object from which data is to be exported. 
  \param cset The MTZdataset object holding the MTZ column wavelength info. 
  \param cxtl The MTZcrystal object holding the MTZ column cell info. */
void MTZfile::export_hkl_data( const HKL_data_base& cdata, const MTZdataset& cset, const MTZcrystal& cxtl, const String mtzpath )
{
  if ( mode != WRITE && mode != APPEND )
    Message::message( Message_fatal( "MTZfile: no file open for write/append" ) );  

  // decode the mtz column info
  String coltype;
  String type  = cdata.type();
  int    ncols = cdata.data_size();
  String names = cdata.data_names();
  std::vector<String> col_names = mtz_assign(mtzpath.tail(),type,names,ncols);
  std::vector<String> dat_names = names.split(" ");

  // now look up the dataset and crystal
  column_map_o colmap;
  colmap.col.resize( ncols, 0 );  // mark all columns as null
  colmap.scl.resize( ncols, 1.0 );
  mmtz_io::mmtz_column  fcol;
  mmtz_io::mmtz_dataset fset;
  mmtz_io::mmtz_crystal fxtl;
  // fill the crystal
  strncpy( fxtl.xname, cxtl.crystal_name().substr( 0, 64 ).c_str(), 64 );
  strncpy( fxtl.pname, cxtl.project_name().substr( 0, 64 ).c_str(), 64 );
  Cell_descr cp = cxtl.descr();
  fxtl.cell[0] = cp.a(); fxtl.cell[3] = cp.alpha_deg();
  fxtl.cell[1] = cp.b(); fxtl.cell[4] = cp.beta_deg();
  fxtl.cell[2] = cp.c(); fxtl.cell[5] = cp.gamma_deg();
  // fill the dataset
  strncpy( fset.dname, cset.dataset_name().substr( 0, 64 ).c_str(), 64 );
  fset.wavel = cset.wavelength();

  // create the mtz columns, and store the mappings
  colmap.list = &cdata;
  for (int c=0; c < ncols; c++) {
    strncpy( fcol.label, col_names[c].substr( 0, 30 ).c_str(), 30 );
    coltype =  MTZ_type_registry::type( dat_names[c] );
    fcol.type[0] = coltype[0];
    fcol.type[1] = '\0';
    colmap.col[c] = mmtz_io::mmtz_add_column( mtzout, &fcol, &fset, &fxtl );
    colmap.scl[c] = MTZ_type_registry::scale( dat_names[c] );
  }

  list_map_o.push_back( colmap );
}


/*! Import data from an MTZ into a CHKL_data object. If they don't
  already exist, then CMTZcrystal and CMTZdataset objects will be created to
  match the MTZ crystal and dataset information for the first MTZ
  column used. These new objects will be children of the parent
  CHKL_info object, and this CHKL_data will be moved to become a child
  of the CMTZdataset object.

  Thus, to import data into a CHKL_data, you must first create a
  CHKL_data anywhere below a parent HKL_info. Then call this method,
  and the object will be moved to a position below the HKL_info
  corresponding to its position in the data hierarchy in the MTZ
  file. The following code imports data, dataset, and crystal from an
  MTZ file:

  \code
  CHKL_info myhkl; // must be given cell, spacegroup, and HKL list.
  ...
  CHKL_data<F_sigF> mydata
  MTZfile file = open_read("in.mtz");
  file.import_chkl_data( mydata, "native/CuKa/[FP,SIGFP]" );
  file.close_read();
  \endcode

  An MTZ column type must be present in the MTZ_type_registry for the
  HKL_data type element name concerned.

  This routine does not actually read any data, but rather marks the
  data to be read when the file is closed.
  \param target The HKL_data object into which data is to be imported. 
  \param mtzpath The MTZ column names, as a path. See \ref MTZpaths
  for details.
  \param path Where to put this in the data hierarchy, overriding the
  MTZ crystal and dataset. */
void MTZfile::import_chkl_data( Container& target, const String mtzpath, const String path )
{
  if ( mode != READ )
    Message::message( Message_fatal( "MTZfile: no file open for read" ) );

  // get this object
  HKL_data_base* hp = dynamic_cast<HKL_data_base*>(&target);
  if ( hp == NULL )
    Message::message( Message_fatal( "MTZfile: import object not HKL_data" ) );

  // and the parent hkl_info
  CHKL_info* chkl = target.parent_of_type_ptr<CHKL_info>();
  if ( chkl == NULL )
    Message::message( Message_fatal( "MTZfile: import HKL_data has no HKL_info" ) );

  // now import the data
  MTZcrystal xtl;
  MTZdataset set;
  import_hkl_data( *hp, set, xtl, mtzpath );

  Container *cxtl, *cset;
  // check for matching crystals
  String xtlname = path.notail().notail().tail();
  if ( xtlname == "" ) xtlname = xtl.crystal_name();
  cxtl = chkl->find_path_ptr( xtlname );
  if ( cxtl == NULL ) {
    cxtl = new CMTZcrystal( *chkl, xtlname, xtl );
    cxtl->set_destroyed_with_parent(); // (use garbage collection for this obj)
  }
  // check for matching datasets
  String setname = path.notail().tail();
  if ( setname == "" ) setname = set.dataset_name();
  cset = cxtl->find_path_ptr( setname );
  if ( cset == NULL ) {
    cset = new CMTZdataset( *cxtl, setname, set );
    cset->set_destroyed_with_parent(); // (use garbage collection for this obj)
  }
  // move the data to the new dataset
  String datname = path.tail();
  if ( datname == "" ) datname = mtzpath.tail();
  target.move( cset->path() + "/" + datname );
}


/*! Export data from a CHKL_data object to an MTZ. The object must
  have a parent Cdataset and CMTZcrystal to provide the MTZ crystal and
  dataset information. The MTZ file will be checked for names matching
  the names of these objects, and the new MTZ columns will be added to
  the corresponding dataset if it exists, otherwise it will be
  created.

  An MTZ column type must be present in the MTZ_type_registry for the
  HKL_data type element name concerned.

  This routine does not actually write any data, but rather marks the
  data to be written when the file is closed.
  \param target The HKL_data object from which data is to be exported.
  \param mtzpath The MTZ column names, as a path. See \ref MTZpaths
  for details. */
void MTZfile::export_chkl_data( Container& target, const String mtzpath )
{
  if ( mode != WRITE && mode != APPEND )
    Message::message( Message_fatal( "MTZfile: no file open for write/append" ) );  

  // get this object
  HKL_data_base* hp = dynamic_cast<HKL_data_base*>( &target );
  if ( hp == NULL )
    Message::message( Message_fatal( "MTZfile: export object not HKL_data" ) );
  // get parent objects
  MTZdataset* dp = target.parent_of_type_ptr<MTZdataset>();
  if ( dp == NULL )
    Message::message( Message_fatal( "MTZfile: HKL_data has no parent MTZdataset" ) );
  MTZcrystal* xp = target.parent_of_type_ptr<MTZcrystal>();
  if ( xp == NULL )
    Message::message( Message_fatal( "MTZfile: HKL_data has no parent MTZcrystal" ) );

  export_hkl_data( *hp, *dp, *xp, mtzpath );
}


} // namespace clipper
