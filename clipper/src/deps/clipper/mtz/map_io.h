/*! \file lib/map_io.h
    Header file for reflection data map importer
*/
//C Copyright (C) 2000-2004 Kevin Cowtan and University of York


//L   This code is distributed under the terms and conditions of the
//L   CCP4 Program Suite Licence Agreement as a CCP4 Library.
//L   A copy of the CCP4 licence can be obtained by writing to the
//L   CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
#ifndef CLIPPER_MAP_IO
#define CLIPPER_MAP_IO


#include "../core/container_map.h"


namespace clipper
{

  //! MAP import/export parent class for clipper objects
  /*! This is the import/export class which can be linked to a CCP4 map
    file and be used to transfer data into or out of a Clipper
    data structure. */
  class MAPfile
  {
   public:
    //! Constructor: does nothing
    MAPfile();
    //! Destructor: close any file that was left open
    ~MAPfile();

    //! Open a file for read access
    void open_read( const String filename_in );
    //! Close a file after reading
    void close_read();
    //! Open a file for read access
    void open_write( const String filename_out );
    //! Close a file after reading
    void close_write();

    //! set cell desription (NXmap write only)
    void set_cell( const Cell& cell );

    //! get file spacegroup
    const Spacegroup& spacegroup() const;
    //! get file cell
    const Cell& cell() const;
    //! get file grid_sampling
    const Grid_sampling& grid_sampling() const;
    //! import data to Xmap
    template<class T> void import_xmap( Xmap<T>& xmap ) const;
    //! export data from Xmap
    template<class T> void export_xmap( const Xmap<T>& xmap );
    //! import data to NXmap
    template<class T> void import_nxmap( NXmap<T>& nxmap ) const;
    //! export data from NXmap
    template<class T> void export_nxmap( const NXmap<T>& nxmap );

  protected:
    enum MAPmode { NONE, READ, WRITE };
    String filename;  //!< filename
    MAPmode mode;     //!< mode

    // header info
    Spacegroup spacegroup_;   //!< map spacegroup
    Cell cell_;               //!< map cell
    Grid_sampling grid_sam_;  //!< cell grid sampling
    Grid_range grid_map_;       //!< map grid extent
  };


} // namespace clipper

#endif
