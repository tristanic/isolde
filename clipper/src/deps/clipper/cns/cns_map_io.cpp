/* map_io.cpp: class file for reflection data  map importer */
//c Copyright (C) 2006 Jon A. Christopher, Kevin Cowtan and University of York
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
//L  You should have received a copy of the CNS licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CNS Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA


#include "cns_map_io.h"

#include <stdio.h>
#include <stdlib.h>

namespace clipper {

  namespace data {
    const char* cns_sg_data[] = { "Null","P 1","-P 1","P 2y","P 2yb","C 2y","P -2y","P -2yc","C -2y","C -2yc","-P 2y","-P 2yb","-C 2y","-P 2yc","-P 2ybc","-C 2yc","P 2 2","P 2c 2","P 2 2ab","P 2ac 2ab","C 2c 2","C 2 2","F 2 2","I 2 2","I 2b 2c","P 2 -2","P 2c -2","P 2 -2c","P 2 -2a","P 2c -2ac","P 2 -2bc","P 2ac -2","P 2 -2ab","P 2c -2n","P 2 -2n","C 2 -2","C 2c -2","C 2 -2c","A 2 -2","A 2 -2c","A 2 -2a","A 2 -2ac","F 2 -2","F 2 -2d","I 2 -2","I 2 -2c","I 2 -2a","-P 2 2","P 2 2 -1n","-P 2 2c","P 2 2 -1ab","-P 2a 2a","-P 2a 2bc","-P 2ac 2","-P 2a 2ac","-P 2 2ab","-P 2ab 2ac","-P 2c 2b","-P 2 2n","P 2 2ab -1ab","-P 2n 2ab","-P 2ac 2ab","-P 2ac 2n","-C 2c 2","-C 2bc 2","-C 2 2","-C 2 2c","-C 2b 2","C 2 2 -1bc","-F 2 2","F 2 2 -1d","-I 2 2","-I 2 2c","-I 2b 2c","-I 2b 2","P 4","P 4w","P 4c","P 4cw","I 4","I 4bw","P -4","I -4","-P 4","-P 4c","P 4ab -1ab","P 4n -1n","-I 4","I 4bw -1bw","P 4 2","P 4ab 2ab","P 4w 2c","P 4abw 2nw","P 4c 2","P 4n 2n","P 4cw 2c","P 4nw 2abw","I 4 2","I 4bw 2bw","P 4 -2","P 4 -2ab","P 4c -2c","P 4n -2n","P 4 -2c","P 4 -2n","P 4c -2","P 4c -2ab","I 4 -2","I 4 -2c","I 4bw -2","I 4bw -2c","P -4 2","P -4 2c","P -4 2ab","P -4 2n","P -4 -2","P -4 -2c","P -4 -2ab","P -4 -2n","I -4 -2","I -4 -2c","I -4 2","I -4 2bw","-P 4 2","-P 4 2c","P 4 2 -1ab","P 4 2 -1n","-P 4 2ab","-P 4 2n","P 4ab 2ab -1ab","P 4ab 2n -1ab","-P 4c 2","-P 4c 2c","P 4n 2c -1n","P 4n 2 -1n","-P 4c 2ab","-P 4n 2n","P 4n 2n -1n","P 4n 2ab -1n","-I 4 2","-I 4 2c","I 4bw 2bw -1bw","I 4bw 2aw -1bw","P 3","P 31","P 32","R 3","-P 3","-R 3","P 3 2","P 3 2\"","P 31 2c (0 0 1)","P 31 2\"","P 32 2c (0 0 -1)","P 32 2\"","R 3 2\"","P 3 -2\"","P 3 -2","P 3 -2\"c","P 3 -2c","R 3 -2\"","R 3 -2\"c","-P 3 2","-P 3 2c","-P 3 2\"","-P 3 2\"c","-R 3 2\"","-R 3 2\"c","P 6","P 61","P 65","P 62","P 64","P 6c","P -6","-P 6","-P 6c","P 6 2","P 61 2 (0 0 -1)","P 65 2 (0 0 1)","P 62 2c (0 0 1)","P 64 2c (0 0 -1)","P 6c 2c","P 6 -2","P 6 -2c","P 6c -2","P 6c -2c","P -6 2","P -6c 2","P -6 -2","P -6c -2c","-P 6 2","-P 6 2c","-P 6c 2","-P 6c 2c","P 2 2 3","F 2 2 3","I 2 2 3","P 2ac 2ab 3","I 2b 2c 3","-P 2 2 3","P 2 2 3 -1n","-F 2 2 3","F 2 2 3 -1d","-I 2 2 3","-P 2ac 2ab 3","-I 2b 2c 3","P 4 2 3","P 4n 2 3","F 4 2 3","F 4d 2 3","I 4 2 3","P 4acd 2ab 3","P 4bd 2ab 3","I 4bd 2c 3","P -4 2 3","F -4 2 3","I -4 2 3","P -4n 2 3","F -4c 2 3","I -4bd 2c 3","-P 4 2 3","-P 4 2 3","-P 4n 2 3","P 4n 2 3 -1n","-F 4 2 3","-F 4c 2 3","F 4d 2 3 -1d","F 4d 2 3 -1cd","-I 4 2 3","-I 4bd 2c 3" };
  }  // namespace data


/*! Constructing an CNSMAPfile does nothing except set the spacegroup
  and flag the object as not attached to any file for either input or
  output.
  The default CNS setting is selected by spacegroup number.
  \param sg_num The spacegroup number of the CNS setting. */
CNSMAPfile::CNSMAPfile( unsigned int sg_num )
{
  mode = NONE;
  if ( sg_num == 0 || sg_num > 230 )
    Message::message( Message_fatal( "CNSMAPfile: invalid spacegroup" ) );
  Spgr_descr sd( data::cns_sg_data[sg_num], Spgr_descr::Hall );
  spacegroup_ = Spacegroup( sd );
}


/*! Constructing an CNSMAPfile does nothing except et the spacegroup
  and flag the object as not attached to any file for either input or
  output.
  \param spacegroup The spacegroup. */
CNSMAPfile::CNSMAPfile( Spacegroup spacegroup )
{
  mode = NONE;
  spacegroup_ = spacegroup;
}


/*! Close any files which were left open. This is particularly
  important since to access the MAP file efficiently, data reads and
  writes are deferred until the file is closed. */
CNSMAPfile::~CNSMAPfile()
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


/*! The file is opened for reading. This CNSMAPfile object will remain
  attached to this file until it is closed. Until that occurs, no
  other file may be opened with this object, however another CNSMAPfile
  object could be used to access another file.
  \param filename_in The input filename or pathname. */
void CNSMAPfile::open_read( const String filename_in )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "CNSMAPfile: open_read - File already open" ) );

  filename = filename_in;
  mode = READ;

  FILE* f = fopen( filename.c_str(), "r" );
  if ( f == NULL )
    Message::message( Message_fatal( "CNSMAPfile: open_read  - Could not read: "+filename ) );
  fclose( f );

  mode = READ;
}


/*! Close the file after reading. */
void CNSMAPfile::close_read()
{
  if ( mode != READ )
    Message::message( Message_fatal( "CNSMAPfile: no file open for read" ) );

  mode = NONE;
}


/*! The file is opened for writing. This will be a new file, created
  entirely from data from within the program, rather than by extending
  an existing file. Similar restrictions apply as for open_read().
  \param filename_out The output filename or pathname. */
void CNSMAPfile::open_write( const String filename_out )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "CNSMAPfile: open_write - File already open" ) );

  filename = filename_out;
  mode = WRITE;
}


/*! Close the file after writing. */
void CNSMAPfile::close_write()
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "CNSMAPfile: no file open for write" ) );

  mode = NONE;
}


/*! When writing an NXmap, the cell for the output map must be set
  using set_cell(). Note that the NXmap is rather more general than
  the CNS map, since it can take an arbitrary rotation+skew
  matrix. The resulting map will only be sensible if the NXmap grid
  skew matrix reflects the supplied cell. This is not possible in the
  general case. (Alternatively, opening an equivalent map for read and
  then closing it again will also set the cell).
  \param cell The cell description for the output map file. */
void CNSMAPfile::set_cell( const Cell& cell )
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "CNSMAPfile: no file open for write" ) );
  cell_ = cell;
}


/*! Get the spacegroup from the MTZ file. \return The spacegroup. */
const Spacegroup& CNSMAPfile::spacegroup() const
{ return spacegroup_; }


/*! Get the base cell from the MTZ file. \return The cell. */
const Cell& CNSMAPfile::cell() const
{ return cell_; }


/*! Get the grid sampling from the MTZ file. \return The grid sampling. */
const Grid_sampling& CNSMAPfile::grid_sampling() const
{ return grid_sam_; }


/*! Import a complete Xmap object. The supplied Xmap object is
  examined, and if any of the parameters (spacegroup, cell, or
  grid_sampling) are unset, then they will be set using values from the
  file. The data is the imported from the file.

  If the spacegroups mismatch, the resulting map will obey its
  spacegroup symmetry, but no expansion will be performed if the file
  has a higher symmetry and only holds an asymmetric unit.

  \param xmap The Xmap to be imported.
*/
template<class T> void CNSMAPfile::import_xmap( Xmap<T>& xmap )
{
  if ( mode != READ )
    Message::message( Message_fatal( "CNSMAPfile: no file open for read" ) );

  char buf[512];
  FILE* f = fopen( filename.c_str(), "r" );
  fgets(buf,512,f);;
  unsigned int n_header_lines;
  if (sscanf(buf,"%d",&n_header_lines) != 1)
    fgets(buf,512,f);; // could be blank line first
  if (sscanf(buf,"%d",&n_header_lines) != 1)
    Message::message( Message_fatal( "CNSMAPfile: can't get number of header lines" ) );
  for (unsigned int i=0;i<n_header_lines;++i)
    fgets(buf,512,f);;

  // get dimensions
  fgets(buf,512,f);;
  int nx,x0,xf,ny,y0,yf,nz,z0,zf;
  sscanf(buf,"%d %d %d %d %d %d %d %d %d",
         &nx,&x0,&xf, &ny,&y0,&yf, &nz,&z0,&zf);

  grid_sam_ = Grid_sampling(nx, ny, nz);
  grid_map_ = Grid_range( Coord_grid( x0, y0, z0 ),
                          Coord_grid( xf, yf, zf ) );
    
  // get cell params
  float cell_param[6];
  fgets(buf,512,f);;
  std::string line(buf);
  for (unsigned int i=0;i<6;++i)
    cell_param[i]=atof(line.substr(i*12,12).c_str());
  
  cell_ = Cell( Cell_descr ( cell_param[0], cell_param[1],
                             cell_param[2], cell_param[3],
                             cell_param[4], cell_param[5]) );
  
  xmap.init(spacegroup_, cell_, grid_sam_);

  //read mode line
  fgets(buf,512,f);;
  if (std::string(buf,3)!=std::string("ZYX"))
    Message::message( Message_fatal( "CNSMAPfile: only ZYX mode supported" ) );

  int z,y,x;
  Xmap_base::Map_reference_coord ix( xmap );
  for (z=z0; z <=zf; ++z) 
  {
    fgets(buf,512,f);; // read slice number line

    unsigned int i=6;
    for (y=y0; y<= yf; ++y)
    {
      for (x=x0; x<=xf; ++x)
      {
        if (i==6)
        {
          fgets(buf,512,f);;
          line=std::string(buf);
          i=0;
        }
        ix.set_coord( Coord_grid( x, y, z) );
        xmap[ix]= T (atof(line.substr(i++*12,12).c_str()));
      }
    }
  }
  fclose(f);
}


/*! 
  \param xmap The Xmap to be exported.
*/
template<class T> void CNSMAPfile::export_xmap( const Xmap<T>& xmap )
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "CNSMAPfile: no file open for write" ) );

  FILE* f = fopen( filename.c_str(), "w" );
  fprintf(f,"\n 1\nFrom clipper Xmap\n");

  // get dimensions

  Grid_sampling gs=xmap.grid_sampling();
  Grid_range gr=xmap.grid_asu();

  //FIXME: verify these ranges. should they be max+1?
  int nx=gs.nu(),ny=gs.nv(),nz=gs.nw();
  int x0=gr.min().u(),y0=gr.min().v(),z0=gr.min().w();
  int xf=gr.max().u(),yf=gr.max().v(),zf=gr.max().w();

  fprintf(f,"%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
          nx,x0,xf,
          ny,y0,yf,
          nz,z0,zf);
  
  Cell cell=xmap.cell();
  fprintf(f,"%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",
          cell.a(),         cell.b(),        cell.c(),
          cell.alpha_deg(), cell.beta_deg(), cell.gamma_deg());

  fprintf(f,"ZYX\n");

  int z,y,x;
  Xmap_base::Map_reference_coord ix( xmap );
  for (z=z0; z <=zf; ++z) 
  {
    fprintf(f,"%8d\n",z);

    unsigned int i=0;
    for (y=y0; y<= yf; ++y)
    {
      for (x=x0; x<=xf; ++x)
      {
        ix.set_coord( Coord_grid( x, y, z) );
        fprintf(f,"%12.5E",float(xmap[ix]));
        i++;
        if (i==6)
        {
          fprintf(f,"\n");
          i=0;
        }
      }
    }
    if ( i != 0 ) fprintf(f,"\n");
  }
  fclose(f);
}

template<class T> void CNSMAPfile::import_nxmap( NXmap<T>& nxmap )
{
  if ( mode != READ )
    Message::message( Message_fatal( "CNSMAPfile: no file open for read" ) );

  char buf[512];
  FILE* f = fopen( filename.c_str(), "r" );
  fgets(buf,512,f);;
  unsigned int n_header_lines;
  if (sscanf(buf,"%d",&n_header_lines) != 1)
    fgets(buf,512,f);; // could be blank line first
  if (sscanf(buf,"%d",&n_header_lines) != 1)
    Message::message( Message_fatal( "CNSMAPfile: can't get number of header lines" ) );
  for (unsigned int i=0;i<n_header_lines;++i)
    fgets(buf,512,f);;

  // get dimensions
  fgets(buf,512,f);;
  int nx,x0,xf,ny,y0,yf,nz,z0,zf;
  sscanf(buf,"%d %d %d %d %d %d %d %d %d",
         &nx,&x0,&xf, &ny,&y0,&yf, &nz,&z0,&zf);

  grid_sam_ = Grid_sampling(nx, ny, nz);
  grid_map_ = Grid_range( Coord_grid( x0, y0, z0 ),
                          Coord_grid( xf, yf, zf ) );
    
  // get cell params
  float cell_param[6];
  fgets(buf,512,f);;
  std::string line(buf);
  for (unsigned int i=0;i<6;++i)
    cell_param[i]=atof(line.substr(i*12,12).c_str());
  
  cell_ = Cell( Cell_descr ( cell_param[0], cell_param[1],
                             cell_param[2], cell_param[3],
                             cell_param[4], cell_param[5]) );
  
  nxmap.init(cell_, grid_sam_, grid_map_);

  //read mode line
  fgets(buf,512,f);;
  if (std::string(buf)!=std::string("ZYX"))
    Message::message( Message_fatal( "CNSMAPfile: only ZYX mode supported" ) );

  int z,y,x;
  NXmap_base::Map_reference_coord ix( nxmap );
  for (z=z0; z <=zf; ++z) 
  {
    fgets(buf,512,f);; // read slice number line

    unsigned int i=6;
    for (y=y0; y<= yf; ++y)
    {
      for (x=x0; x<=xf; ++x)
      {
        if (i==6)
        {
          fgets(buf,512,f);;
          line=std::string(buf);
          i=0;
        }
        ix.set_coord( Coord_grid( x, y, z) );
        nxmap[ix]= T (atof(line.substr(i++*12,12).c_str()));
      }
    }
  }
  fclose(f);
}

template<class T> void CNSMAPfile::export_nxmap( const NXmap<T>& nxmap )
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "CNSMAPfile: no file open for write" ) );

  FILE* f = fopen( filename.c_str(), "w" );
  fprintf(f,"\n 1\nFrom clipper NXmap\n");

  // get dimensions
  Grid gs=nxmap.grid();

  /* cribbed from ccp4_map_io.cpp.  Hope I was smart enough */

  /* Because CNS maps don't allow an arbitrary skew matrix, we have
     to jump through hoops to try and fit an NXmap into one. This code
     will work for those cases where it is possible and produce
     garbage otherwise. Don't even try to understand this unless you
     are pretty smart. */
  // cell (set by user, or from previous map)
  float cp[6];
  cp[0] = cell_.descr().a(); cp[3] = cell_.descr().alpha_deg();
  cp[1] = cell_.descr().b(); cp[4] = cell_.descr().beta_deg ();
  cp[2] = cell_.descr().c(); cp[5] = cell_.descr().gamma_deg();
  // grid (calculated to fit with cell provided - assume angles match)
  Coord_frac c0, c1;
  c0 = nxmap.coord_orth( Coord_map(0,0,0) ).coord_frac(cell_);
  c1 = nxmap.coord_orth( Coord_map( nxmap.grid().nu(), nxmap.grid().nv(),
				    nxmap.grid().nw() ) ).coord_frac(cell_);
  grid_sam_ =
    Grid_sampling( Util::intr( ftype(nxmap.grid().nu())/(c1.u()-c0.u()) ),
		   Util::intr( ftype(nxmap.grid().nv())/(c1.v()-c0.v()) ),
		   Util::intr( ftype(nxmap.grid().nw())/(c1.w()-c0.w()) ) );
  Coord_grid g0 = c0.coord_grid(grid_sam_);
  Coord_grid g1 = g0 + Coord_grid(nxmap.grid()) - Coord_grid(1,1,1);

  int x0=g0.u(), y0=g0.v(), z0=g0.w();
  int xf=g1.u(), yf=g1.v(), zf=g1.w();


  fprintf(f,"%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
          gs.nu(), g0.u(), g1.u(),
          gs.nv(), g0.v(), g1.v(),
          gs.nw(), g0.w(), g1.w());
  

  fprintf(f,"%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",
          cell_.a(),         cell_.b(),        cell_.c(),
          cell_.alpha_deg(), cell_.beta_deg(), cell_.gamma_deg());

  fprintf(f,"ZYX\n");

  int z,y,x;
  NXmap_base::Map_reference_coord ix( nxmap );
  for (z=z0; z <=zf; ++z) 
  {
    fprintf(f,"%8d\n",z);

    unsigned int i=0;
    for (y=y0; y<= yf; ++y)
    {
      for (x=x0; x<=xf; ++x)
      {
        ix.set_coord( Coord_grid( x, y, z) );
        fprintf(f,"%12.5E",float(nxmap[ix]));
        i++;
        if (i==6)
        {
          fprintf(f,"\n");
          i=0;
        }
      }
    }
    if ( i != 0 ) fprintf(f,"\n");
  }
  fclose(f);
}


/*! Import a complete Xmap object. The supplied Xmap object is
  examined, and if any of the parameters (spacegroup, cell, or
  grid_sampling) are unset, then they will be set using values from the
  file. The data is the imported from the file.

  If the spacegroups mismatch, the resulting map will obey its
  spacegroup symmetry, but no expansion will be performed if the file
  has a higher symmetry and only holds an asymmetric unit.

  This version performs a check on the symmatry of the inport map
  file. If the import map contains less than a complete ASU, or if it
  contains more than an ASU and symmetry related positions contain
  inconsistent values, then an error is returned.

  This function is only available for float or double maps.

  \param xmap The Xmap to be imported.
  \return The ASU check flag
*/
template<class T> CNSMAPfile::ASUerror CNSMAPfile::import_xmap_check_asu( Xmap<T>& xmap, T missing )
{
  if ( mode != READ )
    Message::message( Message_fatal( "CNSMAPfile: no file open for read" ) );

  char buf[512];
  FILE* f = fopen( filename.c_str(), "r" );
  fgets(buf,512,f);;
  unsigned int n_header_lines;
  if (sscanf(buf,"%d",&n_header_lines) != 1)
    fgets(buf,512,f);; // could be blank line first
  if (sscanf(buf,"%d",&n_header_lines) != 1)
    Message::message( Message_fatal( "CNSMAPfile: can't get number of header lines" ) );
  for (unsigned int i=0;i<n_header_lines;++i)
    fgets(buf,512,f);;

  // get dimensions
  fgets(buf,512,f);;
  int nx,x0,xf,ny,y0,yf,nz,z0,zf;
  sscanf(buf,"%d %d %d %d %d %d %d %d %d",
         &nx,&x0,&xf, &ny,&y0,&yf, &nz,&z0,&zf);

  grid_sam_ = Grid_sampling(nx, ny, nz);
  grid_map_ = Grid_range( Coord_grid( x0, y0, z0 ),
                          Coord_grid( xf, yf, zf ) );
    
  // get cell params
  float cell_param[6];
  fgets(buf,512,f);;
  std::string line(buf);
  for (unsigned int i=0;i<6;++i)
    cell_param[i]=atof(line.substr(i*12,12).c_str());
  
  cell_ = Cell( Cell_descr ( cell_param[0], cell_param[1],
                             cell_param[2], cell_param[3],
                             cell_param[4], cell_param[5]) );
  
  xmap.init(spacegroup_, cell_, grid_sam_);

  // set missing flag
  T flag;
  Util::set_null( flag );
  xmap = flag;
  ftype maxerr = 0.0;

  //read mode line
  fgets(buf,512,f);;
  if (std::string(buf,3)!=std::string("ZYX"))
    Message::message( Message_fatal( "CNSMAPfile: only ZYX mode supported" ) );

  int z,y,x;
  Xmap_base::Map_reference_coord ix( xmap );
  for (z=z0; z <=zf; ++z) 
  {
    fgets(buf,512,f);; // read slice number line

    unsigned int i=6;
    for (y=y0; y<= yf; ++y)
    {
      for (x=x0; x<=xf; ++x)
      {
        if (i==6)
        {
          fgets(buf,512,f);;
          line=std::string(buf);
          i=0;
        }
        ix.set_coord( Coord_grid( x, y, z) );
	T oldval = xmap[ix];
	T newval = T(atof(line.substr(i++*12,12).c_str()));
        if ( !Util::is_nan( oldval ) && !Util::is_nan( newval ) ) {
          maxerr = Util::max( maxerr, fabs( ftype( newval - oldval ) ) );
          newval = Util::max( oldval, newval );
        }
        xmap[ix] = newval;
      }
    }
  }
  fclose(f);

  // ASU checks
  typedef Xmap_base::Map_reference_index MRI;
  ASUerror asuerr = ASUCORRECT;
  ftype s0(0.0), s1(0.0), s2(0.0);
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() ) {
    T val = xmap[ix];
    if ( Util::is_null( val ) ) {
      asuerr = ASUINCOMPLETE;
      xmap[ix] = missing;
    } else {
      s0 += 1.0;
      s1 += val;
      s2 += val*val;
    }
  }
  if ( s0 > 0.5 ) s2 = sqrt( s2*s0 - s1*s1 ) / s0;
  if ( maxerr > 0.01*s2 ) asuerr = ASUINCONSISTENT;

  return asuerr;
}


// instantiate the template functions for all reasonable types
// xmap import
template void CNSMAPfile::import_xmap<unsigned char>( Xmap<unsigned char>& xmap );
template void CNSMAPfile::import_xmap<char>( Xmap<char>& xmap );
template void CNSMAPfile::import_xmap<unsigned short>( Xmap<unsigned short>& xmap );
template void CNSMAPfile::import_xmap<short>( Xmap<short>& xmap );
template void CNSMAPfile::import_xmap<unsigned int>( Xmap<unsigned int>& xmap );
template void CNSMAPfile::import_xmap<int>( Xmap<int>& xmap );
template void CNSMAPfile::import_xmap<ftype32>( Xmap<ftype32>& xmap );
template void CNSMAPfile::import_xmap<ftype64>( Xmap<ftype64>& xmap );
// xmap export
template void CNSMAPfile::export_xmap<unsigned char>( const Xmap<unsigned char>& xmap );
template void CNSMAPfile::export_xmap<char>( const Xmap<char>& xmap );
template void CNSMAPfile::export_xmap<unsigned short>( const Xmap<unsigned short>& xmap );
template void CNSMAPfile::export_xmap<short>( const Xmap<short>& xmap );
template void CNSMAPfile::export_xmap<unsigned int>( const Xmap<unsigned int>& xmap );
template void CNSMAPfile::export_xmap<int>( const Xmap<int>& xmap );
template void CNSMAPfile::export_xmap<ftype32>( const Xmap<ftype32>& xmap );
template void CNSMAPfile::export_xmap<ftype64>( const Xmap<ftype64>& xmap );
// nxmap import
template void CNSMAPfile::import_nxmap<unsigned char>( NXmap<unsigned char>& nxmap );
template void CNSMAPfile::import_nxmap<char>( NXmap<char>& nxmap );
template void CNSMAPfile::import_nxmap<unsigned short>( NXmap<unsigned short>& nxmap );
template void CNSMAPfile::import_nxmap<short>( NXmap<short>& nxmap );
template void CNSMAPfile::import_nxmap<unsigned int>( NXmap<unsigned int>& nxmap );
template void CNSMAPfile::import_nxmap<int>( NXmap<int>& nxmap );
template void CNSMAPfile::import_nxmap<ftype32>( NXmap<ftype32>& nxmap );
template void CNSMAPfile::import_nxmap<ftype64>( NXmap<ftype64>& nxmap );
// nxmap export
template void CNSMAPfile::export_nxmap<unsigned char>( const NXmap<unsigned char>& nxmap );
template void CNSMAPfile::export_nxmap<char>( const NXmap<char>& nxmap );
template void CNSMAPfile::export_nxmap<unsigned short>( const NXmap<unsigned short>& nxmap );
template void CNSMAPfile::export_nxmap<short>( const NXmap<short>& nxmap );
template void CNSMAPfile::export_nxmap<unsigned int>( const NXmap<unsigned int>& nxmap );
template void CNSMAPfile::export_nxmap<int>( const NXmap<int>& nxmap );
template void CNSMAPfile::export_nxmap<ftype32>( const NXmap<ftype32>& nxmap );
template void CNSMAPfile::export_nxmap<ftype64>( const NXmap<ftype64>& nxmap );
// xmap import and check
template CNSMAPfile::ASUerror CNSMAPfile::import_xmap_check_asu<ftype32>( Xmap<ftype32>& xmap, ftype32 missing );
template CNSMAPfile::ASUerror CNSMAPfile::import_xmap_check_asu<ftype64>( Xmap<ftype64>& xmap, ftype64 missing );


} // namespace clipper
