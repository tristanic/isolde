/*! \file cif/cif_data_io.h
    Header file for reflection data cif importer
*/
//c Copyright (C) 2000-2006 Paul Emsley, Kevin Cowtan and University of York
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


#ifndef CLIPPER_CIF_IO
#define CLIPPER_CIF_IO

#include "../core/hkl_datatypes.h"
#include "../imex.h"

namespace clipper
{

  //! CIF import/export parent class for clipper objects
  /*! This is the import class which can be linked to an cif data
    file and be used to transfer data into a Clipper
    data structure. 
    It is currently a read-only class.
  */
  class CLIPPER_IMEX CIFfile
  {
   public:
    //! Constructor: does nothing
    CIFfile();
    //! Destructor: close any file that was left open
    ~CIFfile();

    //! Open a file for read access
    void open_read( const String filename_in );
    //! Close a file after reading
    void close_read();

    //! get file spacegroup
    const Spacegroup& spacegroup() const;
    //! get file cell
    const Cell& cell() const;
    //! get file resolution
    const Resolution& resolution() const;
    //! get file HKL sampling
    const HKL_sampling& hkl_sampling() const;
    //! get file resolution
    Resolution resolution( const Cell& cell ) const;

    //! read the reflection list from the PHS file
    void import_hkl_info( HKL_info& target );
    //! mark a hkl_data for import from PHS file
    void import_hkl_data( HKL_data_base& cdata );
    //! contains phases predicate
    bool contains_phases_p() const; 


  private:
    enum CIFmode { NONE, READ};

    CIFmode mode;                        //!< file mode
    String filename;                     //!< input/output file
    HKL_data_base* f_sigf_i;             //!< input HKL_data Fs
    HKL_data_base* I_sigI_i;             //!< input HKL_data intensities
    HKL_data_base* f_phi_i;              //!< input HKL_data object (calc sfs)
    HKL_data_base* rfree_i;              //!< input HKL_data rfree flags
    HKL_data_base* phi_fom_i;            //!< input HKL_data phi foms
    HKL_data_base* f_sigf_ano_i;         //!< input anomalous Fs
    HKL_data_base* I_sigI_ano_i;         //!< input anomalous intensities
    HKL_data_base* d_sigd_i;             //!< input anomalous diffs
    HKL_data_base* ABCD_i;               //!< input Hendrickson Lattman coeffs.
    const HKL_data_base* f_sigf_o;       //!< output HKL_data object
    const HKL_data_base* phi_wt_o;       //!< output HKL_data object
    //! object which supplies the hkl list (write mode only)
    const HKL_info* hkl_ptr;
    
    Spacegroup space_group;
    Cell cell_; 
    Resolution resolution_; 
    HKL_sampling hkl_sampling_;
    int set_cell_symm_reso(std::string cif_file_name); 
    int set_cell_symm_reso_by_cif(std::string cif_file_name); 
    int set_cell_symm_reso_by_kludge(std::string cif_file_name); 

    // internal variables, due to the strange and variable nature of
    // mmCIF files, sometime these things are in the file, sometimes
    // not...
    // 
    short int clipper_cell_set_flag; 
    short int clipper_reso_set_flag; 
    short int clipper_symm_set_flag; 

  };

} // namespace clipper

#endif // CLIPPER_CIF_IO


