/* cns_hkl_io.cpp: class file for reflection data  cns_hkl importer */
//C Copyright (C) 2000-2006 Kevin Cowtan and University of York
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA

#include "cns_hkl_io.h"

extern "C" {
#include <stdio.h>
#include <string.h>
}


namespace clipper {


clipper::String cnstok( FILE* f ) {
  String s;
  char c = fgetc( f );
  while ( c > '\0' && ! ( c > ' ' && c != '=' ) ) {
    c = fgetc( f );
  }
  while ( c > '\0' &&   ( c > ' ' && c != '=' ) ) {
    s += toupper(c);
    c = fgetc( f );
  }
  return s;
}


clipper::String cnsrmk( FILE* f ) {
  String s;
  char c = fgetc( f );
  while ( c >= ' ' ) {
    s += toupper(c);
    c = fgetc( f );
  }
  return s;
}


/*! Constructing an CNS_HKLfile does nothing except flag the object as not
  attached to any file for either input or output */
CNS_HKLfile::CNS_HKLfile()
{
  mode = NONE;
}


/*! Close any files which were left open. This is particularly
  important since to access the CNS_HKL file efficiently, data reads and
  writes are deferred until the file is closed. */
CNS_HKLfile::~CNS_HKLfile()
{
  switch ( mode ) {
  case READ:
    close_read(); break;
  case WRITE:
    close_write(); break;
  case NONE:
    break;
  }
}


/*! The file is opened for reading. This CNS_HKLfile object will remain
  attached to this file until it is closed. Until that occurs, no
  other file may be opened with this object, however another CNS_HKLfile
  object could be used to access another file.
  \param filename_in The input filename or pathname. */
void CNS_HKLfile::open_read( const String filename_in )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "CNS_HKLfile: open_read - File already open" ) );

  // open the cns_hkl
  f_sigf_i = NULL; phi_wt_i = NULL; f_phi_i = NULL; abcd_i = NULL; flag_i = NULL;
  filename = filename_in;

  FILE* cns_hkl = fopen( filename.c_str(), "r" );
  if ( cns_hkl == NULL )
    Message::message( Message_fatal( "CNS_HKLfile: open_read  - Could not read: "+filename ) );
  mode = READ;

  // read the headers... up to first NREF
  String s, t, symops;
  ftype cell[6] = { 0.0, 0.0, 0.0, 90.0, 90.0, 90.0 };
  while ( (s=cnstok(cns_hkl)).length() > 0 ) {
    t = (s+"   ").substr(0,4);
    if ( t == "NREF" ) break;
    if ( t == "REMA" ) {
      s = cnsrmk(cns_hkl);
      std::vector<String> toks = s.split( " =" );
      if ( toks[0] == "SYMOP" ) {
	if ( toks.size() > 1 ) symops += toks[1].split(" ();")[0] + ";";
      }
      if ( toks.size() >= 12 ) {
	if ( toks[0] == "A" && toks[2] == "B" && toks[4] == "C" ) {
	  cell[0] = toks[1].f();
	  cell[1] = toks[3].f();
	  cell[2] = toks[5].f();
	}
      }
      if ( toks.size() >= 12 ) {
	if ( toks[6] == "ALPHA" && toks[8] == "BETA" && toks[10] == "GAMMA" ) {
	  cell[3] = toks[7].f();
	  cell[4] = toks[9].f();
	  cell[5] = toks[11].f();
	}
      }
    }
  }
  if ( cell[0]*cell[1]*cell[2] > 0.0 ) {
    Cell_descr cd( cell[0], cell[1], cell[2], cell[3], cell[4], cell[5] );
    cell_ = Cell( cd );
    resolution_ = resolution( cell_ );
    hkl_sampling_ = HKL_sampling( cell_, resolution_ );
    std::cerr << cell_.format() << resolution_.limit() << std::endl;
  }
  if ( symops.length() > 0 ) {
    Spgr_descr sd( symops );
    spacegroup_ = Spacegroup( sd );
    std::cerr << spacegroup_.symbol_xhm() << std::endl;
  }

  fclose( cns_hkl );
}


/*! Close the file after reading. This command also actually fills in
  the data in any HKL_data structures which have been marked for
  import. */
void CNS_HKLfile::close_read()
{
  if ( mode != READ )
    Message::message( Message_fatal( "CNS_HKLfile: no file open for read" ) );

  // make sure the data list is sized
  if ( f_sigf_i != NULL ) f_sigf_i->update();
  if ( phi_wt_i != NULL ) phi_wt_i->update();
  if ( f_phi_i  != NULL ) f_phi_i->update();
  if ( abcd_i   != NULL ) abcd_i->update();
  if ( flag_i   != NULL ) flag_i->update();
  for ( int i = 0; i < fphis_i.size(); i++ ) fphis_i[i].first->update();

  xtype fo[2], pw[2], fc[2], hl[4], fl[1];
  xtype fphis[256][2];  // limit of 256 fphi column groups per CNS file
  HKL hkl;

  FILE* cns_hkl = fopen( filename.c_str(), "r" );
  if ( cns_hkl == NULL )
    Message::message( Message_fatal( "CNS_HKLfile: import_hkl_data  - Could not read: "+filename ) );
  // read the data from the CNS_HKL
  String s, t;
  while ( (s=cnstok(cns_hkl)).length() > 0 ) {
    t = (s+"   ").substr(0,4);
    if ( t == "INDE" ) break;
  }
  while ( s.length() > 0 ) {
    hkl.h() = cnstok(cns_hkl).i();
    hkl.k() = cnstok(cns_hkl).i();
    hkl.l() = cnstok(cns_hkl).i();
    fo[0]=fo[1]=pw[0]=pw[1]=fc[0]=fc[1]=hl[0]=hl[1]=hl[2]=hl[3]=fl[0]=0.0;
    for ( int i = 0; i < fphis_i.size(); i++ ) fphis[i][0]=fphis[i][1]=0.0;
    while ( 1 ) {
      s=cnstok(cns_hkl);
      t = (s+"   ").substr(0,4);
      if ( s.length() == 0 || t == "INDE" ) {
	break;
      } else if ( t == "FOBS" ) {
	fo[0] = cnstok(cns_hkl).f();
	pw[0] = Util::d2rad(cnstok(cns_hkl).f());
      } else if ( t == "SIGM" ) {
	fo[1] = cnstok(cns_hkl).f();
      } else if ( t == "FOM " ) {
	pw[1] = cnstok(cns_hkl).f();
      } else if ( t == "FCAL" ) {
	fc[0] = cnstok(cns_hkl).f();
	fc[1] = Util::d2rad(cnstok(cns_hkl).f());
      } else if ( t == "ABCD" ) {
	hl[0] = cnstok(cns_hkl).f();
	hl[1] = cnstok(cns_hkl).f();
	hl[2] = cnstok(cns_hkl).f();
	hl[3] = cnstok(cns_hkl).f();
      } else if ( t == "TEST" ) {
	fl[0] = cnstok(cns_hkl).f();
      } else {
	for ( int i = 0; i < fphis_i.size(); i++ )
	  if ( s == fphis_i[i].second ) {
	    fphis[i][0] = cnstok(cns_hkl).f();
	    fphis[i][1] = Util::d2rad(cnstok(cns_hkl).f());
	  }
      }
    }
    if ( f_sigf_i != NULL ) f_sigf_i->data_import( hkl, fo );
    if ( phi_wt_i != NULL ) phi_wt_i->data_import( hkl, pw );
    if ( f_phi_i  != NULL ) f_phi_i->data_import( hkl, fc );
    if ( abcd_i   != NULL ) abcd_i->data_import( hkl, hl );
    if ( flag_i   != NULL ) flag_i->data_import( hkl, fl );
    for ( int i = 0; i < fphis_i.size(); i++ )
      fphis_i[i].first->data_import( hkl, fphis[i] );
  }
  fclose( cns_hkl );

  mode = NONE;
}


/*! The file is opened for writing. This will be a new file, created
  entirely from data from within the program, rather than by extending
  an existing file. Similar restrictions apply as for open_read().

  \param filename_out The output filename or pathname. */
void CNS_HKLfile::open_write( const String filename_out )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "CNS_HKLfile: open_write - File already open" ) );

  // open the output cns_hkl
  hkl_ptr = NULL;
  f_sigf_o = NULL; phi_wt_o = NULL; f_phi_o = NULL; abcd_o = NULL; flag_o = NULL;
  filename = filename_out;

  FILE* cns_hkl = fopen( filename.c_str(), "w" );
  if ( cns_hkl == NULL )
    Message::message( Message_fatal( "CNS_HKLfile: open_write - Could not write: "+filename ) );
  fclose( cns_hkl );

  mode = WRITE;
}


/*! Close the file after writing. This command also actually writes
  the data reflection list from the HKL_info object and the data from
  any HKL_data objects which have been marked for import. */
void CNS_HKLfile::close_write()
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "CNS_HKLfile: close_write - no file open for write" ) );

  // export the marked list data to an cns_hkl file
  if ( hkl_ptr == NULL )
    Message::message( Message_fatal( "CNS_HKLfile: close_write - no refln list exported" ) );
  const HKL_info& hklinf = *hkl_ptr;

  HKL hkl;
  xtype x[4];
  float f1, f2, f3, f4;

  FILE* cns_hkl = fopen( filename.c_str(), "w" );
  if ( cns_hkl == NULL )
    Message::message( Message_fatal( "CNS_HKLfile: close_write - Could not write: "+filename ) );
  fprintf( cns_hkl, "NREF=%i\n", hklinf.num_reflections() );
  HKL_info::HKL_reference_index ih;
  for ( ih = hklinf.first(); !ih.last(); ih.next() ) {
    f1 = f2 = f3 = f4 = 0.0;  // default sigf to 0 in case missing
    hkl = ih.hkl();
    fprintf( cns_hkl, "INDE %i %i %i", hkl.h(), hkl.k(), hkl.l() );
    if ( f_sigf_o != NULL ) {
      f_sigf_o->data_export( hkl, x );
      for ( int i = 0; i < 4; i++ ) if ( Util::is_nan(x[i]) ) x[i] = 0.0;
      f1 = float( x[0] );
      f2 = float( x[1] );
    }
    if ( phi_wt_o != NULL ) {
      phi_wt_o->data_export( hkl, x );
      for ( int i = 0; i < 4; i++ ) if ( Util::is_nan(x[i]) ) x[i] = 0.0;
      x[0] = Util::rad2d(x[0]);
      f3 = float( x[0] );
      f4 = float( x[1] );
    }
    fprintf( cns_hkl, " FOBS=%.3f %.3f SIGM=%.3f FOM=%.3f",f1,f3,f2,f4 );
    if ( f_phi_o != NULL ) {
      f_phi_o->data_export( hkl, x );
      for ( int i = 0; i < 4; i++ ) if ( Util::is_nan(x[i]) ) x[i] = 0.0;
      x[1] = Util::rad2d(x[1]);
      fprintf( cns_hkl, " FCAL=%.3f %.3f",float(x[0]),float(x[1]) );
    }
    if ( abcd_o != NULL ) {
      abcd_o->data_export( hkl, x );
      for ( int i = 0; i < 4; i++ ) if ( Util::is_nan(x[i]) ) x[i] = 0.0;
      fprintf( cns_hkl, " HLA=%.1f HLB=%.1f HLC=%.1f HLD=%.1f",float(x[0]),float(x[1]),float(x[2]),float(x[3]) );
    }
    if ( flag_o != NULL ) {
      abcd_o->data_export( hkl, x );
      for ( int i = 0; i < 4; i++ ) if ( Util::is_nan(x[i]) ) x[i] = 0.0;
      fprintf( cns_hkl, " TEST=%i",Util::intr(x[0]) );
    }
    fprintf( cns_hkl, "\n" );
  }
  fclose( cns_hkl );

  mode = NONE;
}


/*! Get the spacegroup from the MTZ file. \return The spacegroup. */
const Spacegroup& CNS_HKLfile::spacegroup() const
{ return spacegroup_; }

/*! Get the base cell from the MTZ file. \return The cell. */
const Cell& CNS_HKLfile::cell() const
{ return cell_; }

/*! Get the resolution limit from the MTZ file. \return The resolution. */
const Resolution& CNS_HKLfile::resolution() const
{ return resolution_; }

/*! Get the HKL sampling from the MTZ file. \return The hkl_sampling. */
const HKL_sampling& CNS_HKLfile::hkl_sampling() const
{ return hkl_sampling_; }

/*! Get the resolution limit from the CNS_HKL file.
  Since a CNS_HKL file does not contain cell information, a Cell object
  must be supplied, which will be used to determine the resultion.
  The result is the resolution determined by the most extreme
  reflection in the file.
  \return The resolution. */
Resolution CNS_HKLfile::resolution( const Cell& cell ) const
{
  if ( mode != READ )
    Message::message( Message_fatal( "CNS_HKLfile: resolution - no file open for read" ) );

  FILE* cns_hkl = fopen( filename.c_str(), "r" );
  if ( cns_hkl == NULL )
    Message::message( Message_fatal( "CNS_HKLfile: resolution - Could not read: "+filename ) );

  // read the reflections from the cns_hkl
  HKL hkl;
  ftype slim = 0.0;
  String s, t;
  while ( (s=cnstok(cns_hkl)).length() > 0 ) {
    t = (s+"   ").substr(0,4);
    if ( t == "INDE" ) {
      hkl.h() = cnstok(cns_hkl).i();
      hkl.k() = cnstok(cns_hkl).i();
      hkl.l() = cnstok(cns_hkl).i();
      slim = Util::max( slim, hkl.invresolsq(cell) );
    }
  }
  fclose( cns_hkl );

  return Resolution( 1.0/sqrt(slim) );
}


/*! Import the list of reflection HKLs from an CNS_HKL file into an
  HKL_info object.

  For most CNS files, there is no cell or spacegroup information: In
  this case that data must already be present in the HKL_info object.
  If the resolution limit of the HKL_info object is lower than the
  limit of the file, any excess reflections will be rejected, as will
  any systematic absences or duplicates.

  For recent CNS files with the cell and symops in the remarks fields,
  the current state of the HKL_info object is discarded and replaced
  with the information from the file.

  \param target The HKL_info object to be initialised. */
void CNS_HKLfile::import_hkl_info( HKL_info& target )
{
  if ( mode != READ )
    Message::message( Message_fatal( "CNS_HKLfile: import_hkl_info - no file open for read" ) );

  FILE* cns_hkl = fopen( filename.c_str(), "r" );
  if ( cns_hkl == NULL )
    Message::message( Message_fatal( "CNS_HKLfile: import_hkl_info - Could not read: "+filename ) );

  if ( !cell_.is_null() && !spacegroup_.is_null() && !resolution_.is_null() )
    target.init( spacegroup_, cell_, resolution_, false );
  ftype slim = target.resolution().invresolsq_limit();

  // read the reflections from the cns_hkl
  std::vector<HKL> hkls;
  HKL hkl;
  String s, t;
  while ( (s=cnstok(cns_hkl)).length() > 0 ) {
    t = (s+"   ").substr(0,4);
    if ( t == "INDE" ) {
      hkl.h() = cnstok(cns_hkl).i();
      hkl.k() = cnstok(cns_hkl).i();
      hkl.l() = cnstok(cns_hkl).i();
      if ( hkl.invresolsq(target.cell()) < slim ) hkls.push_back( hkl );
    }
  }
  fclose( cns_hkl );

  target.add_hkl_list( hkls );
}


/*! Import data from an CNS_HKL file into an HKL_data object.

  This routine does not actually read any data, but rather marks the
  data to be read when the file is closed.

  The data to be read (F_sigF or Phi_fom) will be selected based on
  the type of the HKL_data object.

  \param cdata The HKL_data object into which data is to be imported. */
void CNS_HKLfile::import_hkl_data( HKL_data_base& cdata )
{
  if ( mode != READ )
    Message::message( Message_fatal( "CNS_HKLfile: import_hkl_data - no file open for read" ) );

  if ( cdata.is_null() ) cdata.init( spacegroup_, cell_, hkl_sampling_ );

  if      ( cdata.type() == data32::F_sigF::type()  ) f_sigf_i = &cdata;
  else if ( cdata.type() == data32::Phi_fom::type() ) phi_wt_i = &cdata;
  else if ( cdata.type() == data32::F_phi::type() )   f_phi_i  = &cdata;
  else if ( cdata.type() == data32::ABCD::type() )    abcd_i   = &cdata;
  else if ( cdata.type() == data32::Flag::type() )    flag_i   = &cdata;
  else
    Message::message( Message_fatal( "CNS_HKLfile: import_hkl_data - data must be F_sigF/Phi_fom/F_phi/ABCD/Flag" ) );
}


/*! Import data from an CNS_HKL file into an HKL_data object.

  This routine does not actually read any data, but rather marks the
  data to be read when the file is closed.

  The data to be read must be of type F_phi. A column name is given.

  \param cdata The HKL_data object into which data is to be imported.
  \param name  The column name of the data to import. */
void CNS_HKLfile::import_hkl_data( HKL_data_base& cdata, const String name )
{
  if ( mode != READ )
    Message::message( Message_fatal( "CNS_HKLfile: import_hkl_data - no file open for read" ) );
  if ( cdata.is_null() ) cdata.init( spacegroup_, cell_, hkl_sampling_ );
  if ( cdata.type() == data32::F_phi::type() )
    fphis_i.push_back( std::pair<HKL_data_base*,String>(&cdata,name) );
  else
    Message::message( Message_fatal( "CNS_HKLfile: import_hkl_data - data must be F_phi" ) );
}


/*! Export the list of reflection HKLs to an CNS_HKL file from an
  HKL_info object. */
void CNS_HKLfile::export_hkl_info( const HKL_info& target )
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "CNS_HKLfile: export_hkl_info - no file open for write" ) );
  hkl_ptr = &target;
}


/*! Export data from an HKL_data object into an CNS_HKL file.

  This routine does not actually write any data, but rather marks the
  data to be written when the file is closed.

  The data to be read (F_sigF or Phi_fom) will be selected based on
  the type of the HKL_data object.

  \param cdata The HKL_data object from which data is to be exported. */
void CNS_HKLfile::export_hkl_data( const HKL_data_base& cdata )
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "CNS_HKLfile: export_hkl_data - no file open for write" ) );

  if      ( cdata.type() == data32::F_sigF::type()  ) f_sigf_o = &cdata;
  else if ( cdata.type() == data32::Phi_fom::type() ) phi_wt_o = &cdata;
  else if ( cdata.type() == data32::F_phi::type() )   f_phi_o  = &cdata;
  else if ( cdata.type() == data32::ABCD::type() )    abcd_o   = &cdata;
  else if ( cdata.type() == data32::Flag::type() )    flag_o   = &cdata;
  else
    Message::message( Message_fatal( "CNS_HKLfile: export_hkl_data - data must be F_sigF/Phi_fom/F_phi/ABCD/Flag" ) );
}


} // namespace clipper
