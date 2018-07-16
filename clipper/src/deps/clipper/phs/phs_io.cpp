/* phs_io.cpp: class file for reflection data  phs importer */
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

#include "phs_io.h"

extern "C" {
#include <stdio.h>
}


namespace clipper {


/*! Constructing an PHSfile does nothing except flag the object as not
  attached to any file for either input or output */
PHSfile::PHSfile()
{
  mode = NONE;
}


/*! Close any files which were left open. This is particularly
  important since to access the PHS file efficiently, data reads and
  writes are deferred until the file is closed. */
PHSfile::~PHSfile()
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


/*! The file is opened for reading. This PHSfile object will remain
  attached to this file until it is closed. Until that occurs, no
  other file may be opened with this object, however another PHSfile
  object could be used to access another file.
  \param filename_in The input filename or pathname. */
void PHSfile::open_read( const String filename_in )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "PHSfile: open_read - File already open" ) );

  // open the phs
  f_sigf_i = NULL; phi_wt_i = NULL;
  filename = filename_in;

  FILE* phs = fopen( filename.c_str(), "r" );
  if ( phs == NULL )
    Message::message( Message_fatal( "PHSfile: open_read  - Could not read: "+filename ) );
  fclose( phs );

  mode = READ;
}


/*! Close the file after reading. This command also actually fills in
  the data in any HKL_data structures which have been marked for
  import. */
void PHSfile::close_read()
{
  if ( mode != READ )
    Message::message( Message_fatal( "PHSfile: no file open for read" ) );

  // make sure the data list is sized
  if ( f_sigf_i != NULL ) f_sigf_i->update();
  if ( phi_wt_i != NULL ) phi_wt_i->update();

  int h, k, l;
  HKL hkl;
  xtype x1[2], x2[2];
  float f1, f2, f3, f4;

  char line[240];
  FILE* phs = fopen( filename.c_str(), "r" );
  if ( phs == NULL )
    Message::message( Message_fatal( "PHSfile: import_hkl_data  - Could not read: "+filename ) );
  // read the data from the PHS
  while ( fgets( line, 240, phs ) != NULL ) {
    f1 = f2 = f3 = f4 = 0.0;  // default sigf to 0 in case missing
    sscanf( line, " %i %i %i %f %f %f %f", &h, &k, &l, &f1, &f2, &f3, &f4 );
    hkl = HKL( h, k, l );
    x1[0] = xtype(f1);
    x1[1] = xtype(f4);
    x2[0] = xtype(Util::d2rad(ftype(f3)));
    x2[1] = xtype(f2);
    if ( f_sigf_i != NULL ) f_sigf_i->data_import( hkl, x1 );
    if ( phi_wt_i != NULL ) phi_wt_i->data_import( hkl, x2 );
  }
  fclose( phs );

  mode = NONE;
}


/*! The file is opened for writing. This will be a new file, created
  entirely from data from within the program, rather than by extending
  an existing file. Similar restrictions apply as for open_read().

  \param filename_out The output filename or pathname. */
void PHSfile::open_write( const String filename_out )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "PHSfile: open_write - File already open" ) );

  // open the output phs
  hkl_ptr = NULL;
  f_sigf_o = NULL; phi_wt_o = NULL;
  filename = filename_out;

  FILE* phs = fopen( filename.c_str(), "w" );
  if ( phs == NULL )
    Message::message( Message_fatal( "PHSfile: open_write - Could not write: "+filename ) );
  fclose( phs );

  mode = WRITE;
}


/*! Close the file after writing. This command also actually writes
  the data reflection list from the HKL_info object and the data from
  any HKL_data objects which have been marked for import. */
void PHSfile::close_write()
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "PHSfile: close_write - no file open for write" ) );

  // export the marked list data to an phs file
  if ( hkl_ptr == NULL )
    Message::message( Message_fatal( "PHSfile: close_write - no refln list exported" ) );
  const HKL_info& hklinf = *hkl_ptr;

  HKL hkl;
  xtype x1[2], x2[2];
  float f1, f2, f3, f4;
  f1 = f2 = f3 = f4 = 0.0;

  FILE* phs = fopen( filename.c_str(), "w" );
  if ( phs == NULL )
    Message::message( Message_fatal( "PHSfile: close_write - Could not write: "+filename ) );
  HKL_info::HKL_reference_index ih;
  for ( ih = hklinf.first(); !ih.last(); ih.next() ) {
    hkl = ih.hkl();
    if ( f_sigf_o != NULL ) f_sigf_o->data_export( hkl, x1 );
    if ( phi_wt_o != NULL ) phi_wt_o->data_export( hkl, x2 );
    f1 = float(x1[0]);
    f4 = float(x1[1]);
    f3 = float(Util::rad2d(ftype(x2[0])));
    f2 = float(x2[1]);
    fprintf( phs, "%6i %6i %6i %11.3f %11.3f %11.3f %11.3f\n",
	     hkl.h(), hkl.k(), hkl.l(), f1, f2, f3, f4 );
  }
  fclose( phs );

  mode = NONE;
}


/*! Get the resolution limit from the PHS file.
  Since a PHS file does not contain cell information, a Cell object
  must be supplied, which will be used to determine the resultion.
  The result is the resolution determined by the most extreme
  reflection in the file.
  \return The resolution. */
Resolution PHSfile::resolution( const Cell& cell ) const
{
  if ( mode != READ )
    Message::message( Message_fatal( "PHSfile: resolution - no file open for read" ) );

  HKL hkl;
  int h, k, l;
  char line[240];
  FILE* phs = fopen( filename.c_str(), "r" );
  if ( phs == NULL )
    Message::message( Message_fatal( "PHSfile: resolution - Could not read: "+filename ) );
  // read the reflections from the phs
  ftype slim = 0.0;
  while ( fgets( line, 240, phs ) != NULL ) {
    sscanf( line, " %i %i %i", &h, &k, &l );
    hkl = HKL( h, k, l );
    slim = Util::max( slim, hkl.invresolsq(cell) );
  }
  fclose( phs );

  return Resolution( 1.0/sqrt(slim) );
}


/*! Import the list of reflection HKLs from an PHS file into an
  HKL_info object. If the resolution limit of the HKL_info object is
  lower than the limit of the file, any excess reflections will be
  rejected, as will any systematic absences or duplicates.
  \param target The HKL_info object to be initialised. */
void PHSfile::import_hkl_info( HKL_info& target )
{
  if ( mode != READ )
    Message::message( Message_fatal( "PHSfile: import_hkl_info - no file open for read" ) );

  std::vector<HKL> hkls;
  HKL hkl;
  int h, k, l;
  char line[240];
  FILE* phs = fopen( filename.c_str(), "r" );
  if ( phs == NULL )
    Message::message( Message_fatal( "PHSfile: import_hkl_info - Could not read: "+filename ) );
  // read the reflections from the phs
  ftype slim = target.resolution().invresolsq_limit();
  while ( fgets( line, 240, phs ) != NULL ) {
    sscanf( line, " %i %i %i", &h, &k, &l );
    hkl = HKL( h, k, l );
    if ( hkl.invresolsq(target.cell()) < slim ) hkls.push_back( hkl );
  }
  fclose( phs );

  target.add_hkl_list( hkls );
}


/*! Import data from an PHS file into an HKL_data object.

  This routine does not actually read any data, but rather marks the
  data to be read when the file is closed.

  The data to be read (F_sigF or Phi_fom) will be selected based on
  the type of the HKL_data object.

  \param cdata The HKL_data object into which data is to be imported. */
void PHSfile::import_hkl_data( HKL_data_base& cdata )
{
  if ( mode != READ )
    Message::message( Message_fatal( "PHSfile: import_hkl_data - no file open for read" ) );

  if      ( cdata.type() == data32::F_sigF::type()  ) f_sigf_i = &cdata;
  else if ( cdata.type() == data32::Phi_fom::type() ) phi_wt_i = &cdata;
  else
    Message::message( Message_fatal( "PHSfile: import_hkl_data - data must be F_sigF or Phi_fom" ) );  
}


/*! Export the list of reflection HKLs to an PHS file from an
  HKL_info object. */
void PHSfile::export_hkl_info( const HKL_info& target )
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "PHSfile: export_hkl_info - no file open for write" ) );
  hkl_ptr = &target;
}


/*! Export data from an HKL_data object into an PHS file.

  This routine does not actually write any data, but rather marks the
  data to be written when the file is closed.

  The data to be read (F_sigF or Phi_fom) will be selected based on
  the type of the HKL_data object.

  \param cdata The HKL_data object from which data is to be exported. */
void PHSfile::export_hkl_data( const HKL_data_base& cdata )
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "PHSfile: export_hkl_data - no file open for write" ) );  

  if      ( cdata.type() == data32::F_sigF::type()  ) f_sigf_o = &cdata;
  else if ( cdata.type() == data32::Phi_fom::type() ) phi_wt_o = &cdata;
  else
    Message::message( Message_fatal( "PHSfile: export_hkl_data - data must be F_sigF or Phi_fom" ) );  
}


} // namespace clipper
