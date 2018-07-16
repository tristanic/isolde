/*! \file lib/mtz_io.h
    Header file for reflection data mtz importer
*/
//C Copyright (C) 2000-2004 Kevin Cowtan and University of York


//L   This code is distributed under the terms and conditions of the
//L   CCP4 Program Suite Licence Agreement as a CCP4 Library.
//L   A copy of the CCP4 licence can be obtained by writing to the
//L   CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
#ifndef CLIPPER_MTZ_IO
#define CLIPPER_MTZ_IO


#include "mtz_types.h"

#include <mmtzlib.h>


namespace clipper
{

  //! MTZ import/export type registry
  /*! This class acts as a registry of import/export data types, which
    provide a translation from a data element name to an MTZ column
    type.

    The registry is instantiated statically as \c mtz_type_registry
    and initialised with a list of built-in datatypes.
  */

  class MTZ_type_registry
  {
  public:
    //! constructor: initialise the registry with some built-in types
    MTZ_type_registry();
    //! add a new type to the registry
    static void add_type( const String& name, const String& type, const ftype32& scale );
    //! return MTZ column type
    static String type( const String& name );
    //! return MTZ column type
    static ftype32 scale( const String& name );
  private:
    static char names[200][12];
    static char types[200][4];
    static ftype32 scales[200];
  };




  //! MTZ import/export parent class for clipper objects
  /*! This is the import/export class which can be linked to an mtz
    file and be used to transfer data into or out of a Clipper
    data structure.

    Note that to access the MTZ file efficiently, data reads and
    writes are deferred until the file is closed.

    \anchor MTZpaths \par MTZpaths: MTZ column specification

    Note that the specification of the MTZ column names is quite
    versatile. The MTZ crystal and dataset must be specified, although
    the wildcard '*' may replace a complete name. Several MTZ columns
    will correspond to a single datalist. This may be handled in two
    ways:
    -# A simple name. The corresponding MTZ columns will be named
    after the datalist name, a dot, the datalist type, a dot, and a
    type name for the indivudal column,
    i.e. /crystal/dataset/datalist.listtype.coltype This is the
    default Clipper naming convention for MTZ data.
    -# A comma separated list of MTZ column names enclosed in square
    brackets.  This allows MTZ data from legacy applications to be
    accessed. \n
    Thus, for example, an MTZPATH of
    \code
    native/CuKa/fsigfdata
    \endcode
    expands to MTZ column names of
    \code
    fsigfdata.F_sigF.F
    fsigfdata.F_sigF.sigF
    \endcode
    with a crystal called \c native and a dataset called \c CuKa. An MTZPATH of
    \code
    native/CuKa/[FP,SIGFP]
    \endcode
    expands to MTZ column names of
    \code
    FP
    SIGFP
    \endcode
    with a crystal called \c native and a dataset called \c CuKa.

    \archor MTZ_iotypes \par MTZ_iotypes: Import/export types

    For an HKL_data object to be imported or exported, an MTZ_iotype
    for that datatype must exist in the
    MTZ_iotypes_registry. MTZ_iotypes are defined for all the built-in
    datatypes. If you need to store a user defined type in an MTZ
    file, you must first derive a new MTZ_iotype from MTZ_iotype_base,
    and then register a static instance of that type with the
    MTZ_iotypes_registry. */
  class MTZfile
  {
   public:
    //! Constructor: does nothing
    MTZfile();
    //! Destructor: close any file that was left open
    ~MTZfile();

    //! Open a file for read access
    void open_read( const String filename_in );
    //! Close a file after reading
    void close_read();
    //! Open a file for read access
    void open_append( const String filename_in, const String filename_out );
    //! Close a file after reading
    void close_append();
    //! Open a file for read access
    void open_write( const String filename_out );
    //! Close a file after reading
    void close_write();

    //! get file spacegroup
    const Spacegroup& spacegroup() const;
    //! get file cell
    const Cell& cell() const;
    //! get file resolution
    const Resolution& resolution() const;
    //! read the reflection list from the MTZ
    void import_hkl_list( HKL_info& target );
    //! import parameters of HKL_info object from the MTZ
    void import_hkl_info( HKL_info& target, const bool generate = true );
    //! mark a hkl_data for import from MTZ
    void import_hkl_data( HKL_data_base& cdata, MTZdataset& cset, MTZcrystal& cxtl, const String mtzpath );

    //! write the reflection list to the MTZ (along with cell, spacegroup)
    void export_hkl_info( const HKL_info& target );
    //! mark a hkl_data for export to MTZ
    void export_hkl_data( const HKL_data_base& cdata, const MTZdataset& cset, const MTZcrystal& cxtl, const String mtzpath );

    //! mark a chkl_data container for import from MTZ
    void import_chkl_data( Container& target, const String mtzpath, const String path = "" );
    //! mark a chkl_data container for export to MTZ
    void export_chkl_data( Container& target, const String mtzpath );

  private:
    //! mappings from datalists to mtz columns
    struct column_map_i { HKL_data_base* list; std::vector<int> col; std::vector<ftype32> scl; };
    struct column_map_o { const HKL_data_base* list; std::vector<int> col; std::vector<ftype32> scl; };
    enum MTZmode { NONE, READ, WRITE, APPEND };

    //! mtz object
    mmtz_io::mmtzfile mtzin, mtzout;
    //! index from clipper data lists to mtz columns (by name)
    std::vector<column_map_i> list_map_i;
    std::vector<column_map_o> list_map_o;
    //! object which supplies the hkl list (write mode only)
    const HKL_info* hkl_ptr;
    //! mode
    MTZmode mode;

    //! File spacegroup, cell, resolution
    Spacegroup spacegroup_;
    Cell cell_;
    Resolution resolution_;

    // generic methods
    //! Test if a column name represents a dummy column.
    bool is_virtual_col(const String path) const;
    //! Convert a Clipper datalist name to a list of MTZ column names.
    const std::vector<String> mtz_assign(const String assign, const String type, const String ftype, const int f_size ) const;
    //! Shared open code for both read and append.
    void fetch_header_info();
  };


} // namespace clipper

#endif
