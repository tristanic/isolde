/*! \file lib/cns_hkl_io.h
    Header file for reflection data cns_hkl importer
*/
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


#ifndef CLIPPER_CNS_HKL_IO
#define CLIPPER_CNS_HKL_IO


#include "../core/hkl_datatypes.h"
#include "../imex.h"

namespace clipper
{

  //! CNS_HKL import/export parent class for clipper objects
  /*! This is the import/export class which can be linked to an cns_hkl
    file and be used to transfer data into or out of a Clipper
    data structure. */
  class CLIPPER_IMEX CNS_HKLfile
  {
   public:
    //! Constructor: does nothing
    CNS_HKLfile();
    //! Destructor: close any file that was left open
    ~CNS_HKLfile();

    //! Open a file for read access
    void open_read( const String filename_in );
    //! Close a file after reading
    void close_read();
    //! Open a file for write access
    void open_write( const String filename_out );
    //! Close a file after writing
    void close_write();

    //! get file spacegroup
    const Spacegroup& spacegroup() const;
    //! get file cell
    const Cell& cell() const;
    //! get file resolution
    const Resolution& resolution() const;
    //! get file HKL sampling
    const HKL_sampling& hkl_sampling() const;
    //! get file resolution (if no cell in file)
    Resolution resolution( const Cell& cell ) const;

    //! read the reflection list from the CNS_HKL file
    void import_hkl_info( HKL_info& target );
    //! mark a hkl_data for import from CNS_HKL file
    void import_hkl_data( HKL_data_base& cdata );
    //! mark a hkl_data for import from named column of a CNS_HKL file
    void import_hkl_data( HKL_data_base& cdata, const String name );

    //! write the reflection list to the CNS_HKL file
    void export_hkl_info( const HKL_info& target );
    //! mark a hkl_data for export to CNS_HKL file
    void export_hkl_data( const HKL_data_base& cdata );

  private:
    enum CNS_HKLmode { NONE, READ, WRITE };

    CNS_HKLmode mode;                    //!< file mode
    String filename;                     //!< input/output file
    HKL_data_base* f_sigf_i;             //!< input HKL_data object
    HKL_data_base* phi_wt_i;             //!< input HKL_data object
    HKL_data_base* f_phi_i ;             //!< input HKL_data object
    HKL_data_base* abcd_i;               //!< input HKL_data object
    HKL_data_base* flag_i;               //!< input HKL_data object
    const HKL_data_base* f_sigf_o;       //!< output HKL_data object
    const HKL_data_base* phi_wt_o;       //!< output HKL_data object
    const HKL_data_base* f_phi_o;        //!< output HKL_data object
    const HKL_data_base* abcd_o;         //!< output HKL_data object
    const HKL_data_base* flag_o;         //!< output HKL_data object
    //! object which supplies the hkl list (write mode only)
    const HKL_info* hkl_ptr;
    //! input HKL_data object
    std::vector<std::pair<HKL_data_base*,String> > fphis_i;

    //! File spacegroup, cell, resolution
    Spacegroup spacegroup_;
    Cell cell_;
    Resolution resolution_;
    HKL_sampling hkl_sampling_;
  };
} // namespace clipper

#endif
