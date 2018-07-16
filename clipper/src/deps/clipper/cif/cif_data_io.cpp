
/* cif_io.cpp: class file for reflection data  cif importer               */
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


#include "cif_data_io.h"

extern "C" {
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
}
#include <vector>
#include <fstream>

#ifndef  __MMDB_MMCIF__
#include <mmdb2/mmdb_mmcif_.h>
#undef strchr
#endif


namespace clipper {



/*! Constructing an CIFfile does nothing except flag the object as not
  attached to any file for either input or output */
CIFfile::CIFfile()
{
   mode = NONE;

   // The
   // problem was that we were trying to import into f_phi_i in
   // close_read() when f_phi_i had not been initialised by calling
   // import_hkl_data for it (because it comes from a function that
   // will calculate the phases from the model itself).  So what do we
   // do?  We will set flags in the CIFfile object that mean that the
   // HKL_data is not to be read into.  That value is null.  When we
   // call import_hkl_data, those values get modified to be pointers
   // to HKL_data.  We will add tests for non-nullness to close_read
   // before we import into those HKL_data.
   f_sigf_i = NULL;
   f_phi_i  = NULL;
   rfree_i  = NULL;
   d_sigd_i     = NULL;
   ABCD_i       = NULL;
   I_sigI_i     = NULL;
   I_sigI_ano_i = NULL;
   f_sigf_ano_i = NULL;
   phi_fom_i    = NULL;
   
   // set these when we find them....
   clipper_cell_set_flag = 0; 
   clipper_reso_set_flag = 0; 
   clipper_symm_set_flag = 0; 
}

/*! Close any files which were left open. This is particularly
  important since to access the CIF file efficiently, data reads and
  writes are deferred until the file is closed. */
CIFfile::~CIFfile()
{
  switch ( mode ) {
  case READ:
     
     // No way do we want to close_read() here! If try to get
     // import_hkl_info() and that fails with a Message fatal, then we
     // pass the cif reading code, import_hkl_data() and close_read()
     // code, because if import_hkl_data() fails those functions are
     // useless.  So, on distruction of the CIFfile, we don't want to
     // close_read() - which may well crash in data_import() or some
     // such.
     // 
     // close_read();
     break;
//   case WRITE:
//     close_write(); break;
  case NONE:
    break;
  }
}


/*! The file is opened for reading. This CIFfile object will remain
  attached to this file until it is closed. Until that occurs, no
  other file may be opened with this object, however another CIFfile
  object could be used to access another file.
  \param filename_in The input filename or pathname.
*/
void CIFfile::open_read( const String filename_in )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "CIFfile: open_read - File already open" ) );

  // open the cif
  f_sigf_i = NULL;
  f_phi_i  = NULL;
  rfree_i  = NULL;
  
  filename = filename_in;

  FILE* cif = fopen( filename.c_str(), "r" );
  if ( cif == NULL )
    Message::message( Message_fatal( "CIFfile: open_read  - Could not read: "+filename ) );
  fclose( cif );

  mode = READ;
  filename = filename_in;

  set_cell_symm_reso(filename);
  if ( clipper_cell_set_flag && !clipper_reso_set_flag ) {
    resolution_ = resolution( cell_ );
    clipper_reso_set_flag = 1;
  }
  if ( clipper_cell_set_flag && clipper_reso_set_flag ) {
    hkl_sampling_ = clipper::HKL_sampling( cell_, resolution_ );
  }
}

bool
CIFfile::contains_phases_p() const {

   return 0; 
} 


/*! Close the file after reading. This command also actually fills in
  the data in any HKL_data structures which have been marked for
  import.

  Note that we attempt to read in calculated structure factors too.
  These rely on the tags "F_calc" or "F_calc_au" and "phase_calc".  

  http://pdb.rutgers.edu/mmcif/dictionaries/cif_mm.dic/Categories/refln.html

  It is quite usual then for the HKL_data vector of type F_sigF to
  have a different size to the HKL_data for the F_phi's (often this
  will be zero).  I hope that this will not be a problem.

  We test for f_phi_i being non-null before we import data into it -
  CIFfiles can be used with or without reading calculated phases.

  Note to self: how about we make a function of a CIFfile that says
  whether or not it contains phases...

  Note to self: this text need to be properly marked up in doxygen format.
*/

void CIFfile::close_read()
{
  if ( mode != READ )
     Message::message( Message_fatal( "CIFfile: no file open for read" ) );

  // make sure the data list is sized
  if ( f_sigf_i     != NULL ) f_sigf_i->update();
  if (  f_phi_i     != NULL )  f_phi_i->update();
  if (  rfree_i     != NULL )  rfree_i->update();
  if ( d_sigd_i     != NULL ) d_sigd_i->update();
  if (   ABCD_i     != NULL )   ABCD_i->update();
  if ( I_sigI_i     != NULL ) I_sigI_i->update();
  if ( I_sigI_ano_i != NULL ) I_sigI_ano_i->update();
  if ( f_sigf_ano_i != NULL ) f_sigf_ano_i->update();

  int ret_val = 0;
  int n_calc_data = 0; 
  ::mmdb::InitMatType(); 

  // read the data from the CIF file

  // stat cif_file_name.c_str() here, make sure it exists, is readable.

  ::mmdb::mmcif::File ciffile; 
  int ierr = ciffile.ReadMMCIFFile(filename.c_str());
  int ierr_f;
  int ierr_calc;
  int ierr_rfree_flag;
  int ierr_anom_flag;
  int ierr_intensity_flag;
  int ierr_ABCD_flag;
      
  if (ierr!=::mmdb::mmcif::CIFRC_Ok) { 

    std::string mess = "CIFfile: close_read  - Could not read: ";
    mess += filename.c_str();
    mess += " Dirty mmCIF file? ";
    Message::message( Message_fatal( mess )); 
  } else { 

    for(int i=0; i<ciffile.GetNofData(); i++) {  
	 
      ::mmdb::mmcif::PData data = ciffile.GetCIFData(i);
	 
      for (int icat=0; icat<data->GetNumberOfCategories(); icat++) { 

	::mmdb::mmcif::PCategory cat = data->GetCategory(icat);

	std::string cat_name(cat->GetCategoryName()); 

	::mmdb::mmcif::PLoop mmCIFLoop = data->GetLoop(cat_name.c_str());

	if (mmCIFLoop == NULL) { 

	} else {

	  if (cat_name == "_refln") { 

	    // success:
	    ret_val = 1;

	    int h,k,l;
	    ::mmdb::realtype F, sigF, Fcalc, phi, F_plus, sigF_plus, F_minus, sigF_minus, I, sigI;
	    ::mmdb::realtype I_plus, sigI_plus, I_minus, sigI_minus, D, D_sigma;
	    ::mmdb::realtype HLA, HLB, HLC, HLD;
	    xtype x1[2], x_anom[5], x_hl[4];

	    for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {

	      // Note: modification of references
	      // 
	      ierr  = mmCIFLoop->GetInteger(h, "index_h", j); 
	      ierr += mmCIFLoop->GetInteger(k, "index_k", j); 
	      ierr += mmCIFLoop->GetInteger(l, "index_l", j);
	      if ( !ierr ) {

		// Measured data can be given using F_meas or
		// F_meas_au (arbitary units).  So lets try to
		// read a F_meas first then try F_meas_au if that
		// fails.

		ierr_f = mmCIFLoop->GetReal(   F, "F_meas", j);
		ierr_f += mmCIFLoop->GetReal(sigF, "F_meas_sigma", j);
		if (ierr_f) {
		  ierr_f  = mmCIFLoop->GetReal(   F, "F_meas_au", j);
		  ierr_f += mmCIFLoop->GetReal(sigF, "F_meas_sigma_au", j);
		}
		if ( F < -0.9e10 || sigF < -0.9e10 ) ierr_f++;

		ierr_calc = mmCIFLoop->GetReal(Fcalc, "F_calc", j);
		if ( ierr_calc ) {
		  ierr_calc =  mmCIFLoop->GetReal(Fcalc, "F_calc_au", j);
		}
		ierr_calc += mmCIFLoop->GetReal(phi, "phase_calc", j);
		if ( Fcalc < -0.9e10 || phi < -0.9e10 ) ierr_calc++;

		if ( !ierr_f ) {
		  if ( f_sigf_i != NULL ) if ( !f_sigf_i->is_null() ) { 
		    x1[0] = F;
		    x1[1] = sigF; 
		    f_sigf_i->data_import( HKL(h,k,l), x1 );
		  }
		}

		if (! ierr_calc ) {
		   if ( f_phi_i  != NULL )
		      if ( !f_phi_i->is_null() ) { 
			 x1[0] = Fcalc;
			 x1[1] = clipper::Util::d2rad(phi);
			 f_phi_i->data_import( HKL(h,k,l), x1);
			 n_calc_data++;
		      }
		}

		// RFree flag
		if (rfree_i) { 
		   char *s = mmCIFLoop->GetString("status", j, ierr_rfree_flag);
		   if (! ierr_rfree_flag) { 
		      x1[0] = -1; 
		      if (s) {
			 // std::cout << "read status :" << s << ":" << std::endl;
			 // could be also "x"; not observed, just listed
			 if (! strncmp(s, "o", 1)) {
			    x1[0] = 1;
			 } else {
			    if (! strncmp(s, "f", 1)) {
			       x1[0] = 0;
			    }
			 }
		      }
		      // 
		      // Every output reflection gets one of these, even if
		      // it is set to -1.  Is that the correct thing to do?
		      // I have looked at over 200 recent sf mmCIFs and if
		      // the contain conventional data then they all have
		      // the status flag.
		      //
		      // However, SFS from EM data do not
		      //
		      rfree_i->data_import( HKL(h,k,l), x1);
		   }
		}

		
		// Anomalous Fs
		//
		// (only add if all 4 F+ and sigF+ and F- and sigF- are
		// present), that is, don't add an F if a sigma is
		// missing.
		//
		ierr_anom_flag = mmCIFLoop->GetReal(F_plus, "pdbx_F_plus", j);
		if (! ierr_anom_flag) {
		   ierr_anom_flag = mmCIFLoop->GetReal(sigF_plus, "pdbx_F_plus_sigma", j);
		   if (! ierr_anom_flag) {
		      ierr_anom_flag = mmCIFLoop->GetReal(F_minus, "pdbx_F_minus", j);
		      if (! ierr_anom_flag) {
			 ierr_anom_flag = mmCIFLoop->GetReal(sigF_minus, "pdbx_F_minus_sigma", j);
			 if (! ierr_anom_flag) {
			    x_anom[0] = F_plus;
			    x_anom[1] = sigF_plus;
			    x_anom[2] = F_minus;
			    x_anom[3] = sigF_minus;
			    x_anom[4] = 1.0; // no covarience in cif files, hack a value.
			    if (f_sigf_ano_i)
			       f_sigf_ano_i->data_import(HKL(h,k,l), x_anom);
			 }
		      }
		   }
		}

		// Anomalous Differences (on F, presumably)
		// 
		ierr_anom_flag = mmCIFLoop->GetReal(D, "pdbx_anom_difference", j);
		if (! ierr_anom_flag) {
		   ierr_anom_flag = mmCIFLoop->GetReal(D, "pdbx_anom_difference_sigma", j);
		   if (! ierr_anom_flag) {
		      x1[0] = D;
		      x1[1] = D_sigma;
		      if (d_sigd_i)
			 d_sigd_i->data_import(HKL(h,k,l), x1);
		   }
		}


		// Intensities
		//
		ierr_intensity_flag = mmCIFLoop->GetReal(I, "intensity_meas", j);
		if (! ierr_intensity_flag) {
		   ierr_intensity_flag = mmCIFLoop->GetReal(sigI, "intensity_sigma", j);
		   if (! ierr_intensity_flag) {
		      x1[0] = I;
		      x1[1] = sigI;
		      if (I_sigI_i)
			 I_sigI_i->data_import(HKL(h,k,l), x1);
		   }
		}

		// Anomalous Intensities
		//
		ierr_intensity_flag = mmCIFLoop->GetReal(I_plus, "pdbx_I_plus", j);
		if (! ierr_intensity_flag) {
		   ierr_intensity_flag = mmCIFLoop->GetReal(I_plus, "pdbx_I_plus_sigma", j);
		   if (! ierr_intensity_flag) {
		      ierr_intensity_flag = mmCIFLoop->GetReal(I_minus, "pdbx_I_minus", j);
		      if (! ierr_intensity_flag) {
			 ierr_intensity_flag = mmCIFLoop->GetReal(I_minus, "pdbx_I_minus_sigma", j);
			 if (! ierr_intensity_flag) {
			    x_anom[0] = I_plus;
			    x_anom[1] = sigI_plus;
			    x_anom[2] = I_minus;
			    x_anom[3] = sigI_minus;
			    x_anom[4] = 1.0; // no covarience in cif files, hack a value.
			    if (I_sigI_ano_i)
			       I_sigI_ano_i->data_import(HKL(h,k,l), x_anom);
			 }
		      }
		   }
		}

		// Hendrickson Lattman coefficients (not many files have these)
		//
		ierr_ABCD_flag = mmCIFLoop->GetReal(HLA, "pdbx_HLA", j);
		if (! ierr_ABCD_flag) {
		   ierr_ABCD_flag = mmCIFLoop->GetReal(HLB, "pdbx_HLB", j);
		   if (! ierr_ABCD_flag) {
		      ierr_ABCD_flag = mmCIFLoop->GetReal(HLC, "pdbx_HLC", j);
		      if (! ierr_ABCD_flag) {
			 ierr_ABCD_flag = mmCIFLoop->GetReal(HLD, "pdbx_HLD", j);
			 if (! ierr_ABCD_flag) {
			    x_hl[0] = HLA;
			    x_hl[1] = HLB;
			    x_hl[2] = HLC;
			    x_hl[3] = HLD;
			    if (ABCD_i)
			       ABCD_i->data_import(HKL(h,k,l), x_hl);
			 }
		      }
		   }
		} 
	      }
	    }
	  }
	}
      }
    }
  }
  mode = NONE;
}

/*! Get the spacegroup from the MTZ file. \return The spacegroup. */
const Spacegroup& CIFfile::spacegroup() const
{ return space_group; }

/*! Get the base cell from the MTZ file. \return The cell. */
const Cell& CIFfile::cell() const
{ return cell_; }

/*! Get the resolution limit from the MTZ file. \return The resolution. */
const Resolution& CIFfile::resolution() const
{ return resolution_; }

/*! Get the HKL sampling from the MTZ file. \return The hkl_sampling. */
const HKL_sampling& CIFfile::hkl_sampling() const
{ return hkl_sampling_; }

/*! Get the resolution limit from the CIF file.
  Since a CIF file does not contain cell information, a Cell object
  must be supplied, which will be used to determine the resultion.
  The result is the resolution determined by the most extreme
  reflection in the file.
  \return The resolution. */
Resolution CIFfile::resolution( const Cell& cell ) const
{
  if ( mode != READ )
    Message::message( Message_fatal( "CIFfile: resolution - no file open for read" ) );

  HKL hkl;
  int h, k, l;
  //char line[240];
  ftype slim = 0.0;
  FILE* cif = fopen( filename.c_str(), "r" );
  if ( cif == NULL )
    Message::message( Message_fatal( "CIFfile: resolution  - Could not read: "+filename ) );
  ::mmdb::mmcif::File ciffile; 
  int ierr = ciffile.ReadMMCIFFile(filename.c_str()); 
  if (ierr!=::mmdb::mmcif::CIFRC_Ok) { 
    std::string mess = "CIFfile: resolution  - Could not read: ";
    mess += filename.c_str();
    mess += ". Dirty mmCIF file? "; 
    Message::message( Message_warn( mess )); 
  } else { // read the reflections from the phs
    for(int i=0; i<ciffile.GetNofData(); i++) {  
      ::mmdb::mmcif::PData data = ciffile.GetCIFData(i);
      for (int icat=0; icat<data->GetNumberOfCategories(); icat++) { 
	::mmdb::mmcif::PCategory cat = data->GetCategory(icat);
	std::string cat_name(cat->GetCategoryName()); 
	::mmdb::mmcif::PLoop mmCIFLoop = data->GetLoop( (char *) cat_name.c_str() );
	if (mmCIFLoop == NULL) { 
	} else {
	  if (cat_name == "_refln") { 
	    for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) { 
	      ierr  = mmCIFLoop->GetInteger(h, "index_h", j); 
	      ierr += mmCIFLoop->GetInteger(k, "index_k", j); 
	      ierr += mmCIFLoop->GetInteger(l, "index_l", j);
	      if (!ierr) {
		hkl = HKL( h, k, l );
		slim = Util::max( slim, hkl.invresolsq(cell) );
	      }
	    }
	  }
	}
      }
    }
  }
  fclose( cif );

  return Resolution( 1.0/sqrt(slim) );
}


/*! Import the list of reflection HKLs from an CIF file into an
  HKL_info object. At the start of the routine try to determine the
  space group, cell and resolution. If the resolution limit was found
  and if the resolution limit of the HKL_info object is lower than the
  limit of the file, any excess reflections will be rejected, as will
  any systematic absences or duplicates.  If the resolution is not
  found then the resolution will be determined from the input hkl data
  and \parm target will be init()ed before returning from this function.
  \param target The HKL_info object to be initialised. */
void CIFfile::import_hkl_info( HKL_info& target )
{
   std::vector<HKL> hkls;
   //int h, k, l;

  if ( mode != READ )
    Message::message( Message_fatal( "CIFfile: import_hkl_info - no file open for read"+filename ) );

  std::string new_str = filename; 

  int icell = set_cell_symm_reso(new_str);
  
  if (! icell) {
     if (! clipper_cell_set_flag)
	Message::message( Message_fatal( "CIFfile: import_hkl_info - error getting cell "+filename ) );
     
     if (! clipper_symm_set_flag)
	Message::message( Message_fatal( "CIFfile: import_hkl_info - error getting symm "+filename ) );
  } else {
     // we have the cell and symmetry, proceed.

     if (! clipper_reso_set_flag)
	resolution_.init(2.0); // just a dummy value, we will set it to
			      // the right value once we know it after
			      // we have read in all the reflection
			      // hkls.  If clipper_reso_set_flag *was*
			      // set, then we have initialised
			      // resolution already.
       
     // import any missing params
     target.init( space_group, cell_, resolution_ );

     FILE* cif = fopen( filename.c_str(), "r" );
     if ( cif == NULL )
	Message::message( Message_fatal( "CIFfile: import_hkl_info - Could not read: "+filename ) );
     fclose ( cif );
     // read the reflections from the cif
     ftype slim = target.resolution().invresolsq_limit();
     ftype tmp_lim = 0.0; 
  
     ::mmdb::mmcif::File ciffile; 
     int ierr = ciffile.ReadMMCIFFile(filename.c_str()); 
      
     if (ierr!=::mmdb::mmcif::CIFRC_Ok) { 
     
	std::string mess = "CIFfile: import_hkl_data  - Could not read: ";
	mess += filename.c_str();
	mess += ". Dirty mmCIF file? "; 
	Message::message( Message_fatal( mess )); 

     } else { 

	for(int i=0; i<ciffile.GetNofData(); i++) {  
	 
	   ::mmdb::mmcif::PData data = ciffile.GetCIFData(i);
	 
	   for (int icat=0; icat<data->GetNumberOfCategories(); icat++) { 

	      ::mmdb::mmcif::PCategory cat = data->GetCategory(icat);

	      std::string cat_name(cat->GetCategoryName()); 

	      ::mmdb::mmcif::PLoop mmCIFLoop = data->GetLoop(cat_name.c_str() );

	      if (mmCIFLoop) {

		 if (cat_name == "_refln") { 

		    int h,k,l;

		    for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) { 

		       // ("modification of reference")
		       // 
		       ierr  = mmCIFLoop->GetInteger(h, "index_h", j); 
		       ierr += mmCIFLoop->GetInteger(k, "index_k", j); 
		       ierr += mmCIFLoop->GetInteger(l, "index_l", j);

		       if (!ierr) {
			  HKL hkl(h,k,l);
			  if (clipper_reso_set_flag) {
			     if ( hkl.invresolsq(target.cell()) < slim ) {
				hkls.push_back(hkl);
			     }
			  } else {
			     // resolution had not been set
			     hkls.push_back(hkl);
			     if (hkl.invresolsq(target.cell()) > tmp_lim) {
				tmp_lim = hkl.invresolsq(target.cell());
			     }
			  } 
		       }
		    }
		 } 
	      }
	   }
	}
     }
     // ::mmdb::mmcif::File is closed on destruction of ciffile.

     // Now we can initialise target properly, if we OK reading the hkls.
     if (!clipper_reso_set_flag) {
	if (tmp_lim > 0.0) { // the starting value
	   resolution_.init( 1/sqrt(tmp_lim)); // 210170 // just above
	   target.init( space_group, cell_, resolution_ );
	   std::cout << "Resolution limit set to " << resolution_.limit() << std::endl; 
	} else {
	   std::cout << "Disaster couldn't set resolution" << std::endl;
	}
     }
  }
  // Quiet! Silence is the clipper way :)
  // std::cout << "import_hkl_info read " << hkls.size() << " hkls" << std::endl; 
  target.add_hkl_list( hkls );
}



/*! Import data from an CIF file into an HKL_data object.

  This routine does not actually read any data, but rather marks the
  data to be read when the file is closed.

  The data to be read (F_sigF or Phi_fom) will be selected based on
  the type of the HKL_data object.

  \param cdata The HKL_data object into which data is to be imported. */
void CIFfile::import_hkl_data( HKL_data_base& cdata )
{
  if ( mode != READ )
    Message::message( Message_fatal( "CIFfile: import_hkl_data - no file open for read" ) );

  if ( cdata.is_null() ) cdata.init( space_group, cell_, hkl_sampling_ );

  if      ( cdata.type() == data32::F_sigF::type() )        f_sigf_i = &cdata;
  else if ( cdata.type() == data32::F_phi::type()  )         f_phi_i = &cdata;
  else if ( cdata.type() == data32::I_sigI::type() )        I_sigI_i = &cdata;
  else if ( cdata.type() == data32::Flag::type()   )         rfree_i = &cdata;
  else if ( cdata.type() == data32::Phi_fom::type())       phi_fom_i = &cdata;
  else if ( cdata.type() == data32::F_sigF_ano::type()) f_sigf_ano_i = &cdata;
  else if ( cdata.type() == data32::I_sigI_ano::type()) I_sigI_ano_i = &cdata;
  else if ( cdata.type() == data32::ABCD::type())             ABCD_i = &cdata;
  else if ( cdata.type() == data32::D_sigD::type())         d_sigd_i = &cdata;
  else {
     clipper::String m = "CIFfile: import_hkl_data error";
     m += " - data must be F_sigF, F_phi or Flag";
     Message::message( Message_fatal( m ) );
  }
}

// return non-zero if we find a cell and symmetry and can set them.
// 
int 
CIFfile::set_cell_symm_reso_by_cif(std::string cif_file_name) {

   FILE* cif = fopen( cif_file_name.c_str(), "r" );
   if ( cif == NULL )
      Message::message( Message_fatal( "CIFfile: set_cell_symm_reso_by_cif - Could not read: "+cif_file_name ) );
   fclose( cif );
   
   ::mmdb::mmcif::File ciffile; 
   int ierr = ciffile.ReadMMCIFFile(cif_file_name.c_str()); 
   
   if (ierr!=::mmdb::mmcif::CIFRC_Ok) { 
      
      std::string mess = "CIFfile: set_cell_symm_reso_by_cif: Could not read: ";
      mess += filename.c_str();
      mess += " dirty mmCIF file? "; 
      Message::message( Message_fatal( mess )); 
      
   } else { 
      
      for(int i=0; i<ciffile.GetNofData(); i++) {  
	 
	 ::mmdb::mmcif::PData data = ciffile.GetCIFData(i);
	 // data is the same as mmCIF in Eugenes example

	 ::mmdb::mmcif::PStruct mmCIFStruct = data->GetStructure("_cell");

	 char *cat_name_str;
	 std::string cat_name;
	    
	 if (mmCIFStruct != NULL) { 

	    cat_name_str = mmCIFStruct->GetCategoryName();
	    
	    if (cat_name_str == NULL) { 
	       std::cout << "null cat_name_str" << std::endl; 
	    } else {
	       
	       cat_name = cat_name_str; 
	       if (cat_name == "_cell") {
		  ::mmdb::realtype a,b,c,alpha,beta,gamma;
		  int ierr = 0; 
		  ierr += mmCIFStruct->GetReal (a,     "length_a");
		  ierr += mmCIFStruct->GetReal (b,     "length_b");
		  ierr += mmCIFStruct->GetReal (c,     "length_c");
		  ierr += mmCIFStruct->GetReal (alpha, "angle_alpha");
		  ierr += mmCIFStruct->GetReal (beta,  "angle_beta");
		  ierr += mmCIFStruct->GetReal (gamma, "angle_gamma");
		  
		  if (! ierr) {
		     // set clipper cell
		     clipper_cell_set_flag = 1;
		     cell_.init(clipper::Cell_descr(a,b,c,alpha,beta,gamma));
		     // std::cout << "got cell from cif: "
		     // << cell_.format() << std::endl;
		  }
	       }
	    }
	 }

	 // Try reading symmetry construction
	 mmCIFStruct = data->GetStructure("_symmetry");

	 if (mmCIFStruct != NULL) { 
	    cat_name_str = mmCIFStruct->GetCategoryName();
	    if (cat_name_str != NULL) { 
	       cat_name = cat_name_str; 
	       if ( cat_name == "_symmetry") {
		  int ierr; 
		  char *str = mmCIFStruct->GetString("space_group_name_H-M",ierr);
		  if (! ierr) {
		     std::string hmsymm(str);
		     if (str) { 
			hmsymm = str;
			space_group.init(clipper::Spgr_descr(str));
			clipper_symm_set_flag = 1;
// 			std::cout << "INFO space_group from symmetry in cif: "
// 				  << space_group.descr().symbol_hm() << std::endl;
		     }
		  }
	       }
	    }
	 }

	 // Have another try for symmetry: e.g. from shelx cif files:
	 // (such files tell us the symmetry operators, so we can get
	 // the space group from those rather than the name using a
	 // clipper function).
	 for (int icat=0; icat<data->GetNumberOfCategories(); icat++) { 

	    ::mmdb::mmcif::PCategory cat = data->GetCategory(icat);

	    std::string cat_name(cat->GetCategoryName()); 

	    ::mmdb::mmcif::PLoop mmCIFLoop = data->GetLoop(cat_name.c_str() );
	    
	    if (mmCIFLoop) { 
	       if (cat_name == "_symmetry_equiv") {
		  std::string symmetry_ops("");
		  for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {
		     char *str = mmCIFLoop->GetString("pos_as_xyz", j, ierr);
		     symmetry_ops += str;
		     symmetry_ops += " ; ";
		  }
		  if (symmetry_ops != "") {
		     clipper_symm_set_flag = 1;
		     space_group.init(clipper::Spgr_descr(symmetry_ops));
		  }
	       }
	    }
	 }

	 
	 // Reflection meta data:
	 // 
	 mmCIFStruct = data->GetStructure("_reflns");
	 if (mmCIFStruct != NULL) { 
	    cat_name_str = mmCIFStruct->GetCategoryName();
	    if (cat_name_str != NULL) { 
	       cat_name = cat_name_str; 
	       if ( cat_name == "_reflns") {
		  ::mmdb::realtype reso;
		  int ierr = mmCIFStruct->GetReal(reso, "d_resolution_high");
		  if (! ierr ) {
		     clipper_reso_set_flag = 1;
		     resolution_.init(reso);
// 		     std::cout << "got resolution from cif: "
// 			       << resolution_.limit() << std::endl;
		  }
	       }
	    }
	 }
      }
   }
   return
      clipper_symm_set_flag &&
      clipper_cell_set_flag; 
}

// return non-zero if we find a cell and symmetry and can set them.
// 
int 
CIFfile::set_cell_symm_reso(std::string cif_file_name) { 
// return non-zero if we find a cell and symmetry and can set them.

   if (set_cell_symm_reso_by_cif(cif_file_name)) { 
      return 1; 
   } else {
      return set_cell_symm_reso_by_kludge(cif_file_name);
   }
}

int 
CIFfile::set_cell_symm_reso_by_kludge(std::string cif_file_name) { 
   
   // Cif files from the RCSB do not contain cell and symmetry
   // information in cif format (now how totally crap is that? - but I
   // won't go into that now).
   // 
   // However, they do provide us with comment lines that look like a
   // conventional pdb file.  e.g:
   // 
   // #CRYST1  114.300  114.300  155.200  90.00  90.00 120.00 P 61 2 2     12
   // #REMARK 290 SYMMETRY OPERATORS FOR SPACE GROUP: P 61                   
   // 
   // So the plan is the open the file conventionally and grab the
   // cell and symmetry then after we have got to the "#CRYST1" field
   // or find "loop_", we close up.  Notice how we don't grab the
   // space group symbol, since this is fraught with parsing
   // difficulties, instead we grab the symmetry operators and use
   // clipper cleverness to determine the space group.
   // 
   // Return the success status. True is success, 0 failure.
   //
   //
      
   // std::string line; 
   char word[800];
   short int read_symm_flag = 0;   // 2 is good for read
   short int cell_coming_flag = 0; // 1 is good for read
   short int local_found_symm = 0; 
   short int local_found_cell = 0; 
   //short int local_found_reso = 0; 
   std::string clipper_symm_string(""); 
   std::vector<double> cell_bits; 

   std::ifstream from;
   float reso; 
   short int read_reso_flag = 0; 

   from.open(cif_file_name.c_str()); 

   while (from >> &word[0]) { 
 
      if (read_reso_flag == 1) { 

	 // We can do reso immediately.  The others (cell and
	 // symmetry) we set flags and do them at the end.
	 // 
	 char **endptr = new (char *); 
	 reso = strtod(word, endptr); 
	 if (endptr != NULL &&  *endptr != word) { 
	    resolution_.init(reso); 
	    std::cout << " Found reso: " << resolution_.limit() << std::endl; 
	    clipper_reso_set_flag = 1; 
	 }
	 delete endptr;
	 endptr = 0;
	 read_reso_flag = 0; 
      } 

      if (! strncmp(word, "RESOLUTION.", 11)) { 
	 read_reso_flag = 1; 
      }

      if (! strncmp(word, "NNNMMM",6)) { 
	 read_symm_flag = 1; 
      }

      if (! strncmp(word, "OPERATOR",6)) { 
	 read_symm_flag++;
      }

      if (! strncmp(word, "WHERE",5)) { 
	 read_symm_flag =0;
      }

      if (read_symm_flag == 2) { 
	 
	 if (strchr(word,'X') && 
	     strchr(word,'Y') && 
	     strchr(word,'Z') ) { 

	    clipper_symm_string += word;
	    clipper_symm_string += ";";
	    local_found_symm = 1; 
	 }
      }

      // position dependent code (should be before CRYST1 comparison)
      if (cell_coming_flag == 1) { 
	 cell_bits.push_back(atof(word));
	 local_found_cell = 1; 
      } 

      if (! strncmp(word, "#CRYST1", 7)) {
	 cell_coming_flag = 1; 
      }

      if (cell_bits.size() == 6)
	 break; 

      if (! strncmp(word, "_refln.index_h",14)) break; 
   }
      
   from.close(); 

   if (clipper_symm_set_flag == 0) {
      if (local_found_symm) {
	 space_group.init(clipper::Spgr_descr(clipper_symm_string));
	 clipper_symm_set_flag = 1; 
	 std::cout << " Symm: " << space_group.descr().symbol_hm() << std::endl;
      }
   }
   
   if (clipper_cell_set_flag == 0) {
      if (local_found_cell) {
	 if (cell_bits.size() == 6 && clipper_symm_string != "") { 
	    if ( !clipper::Util::is_nan(cell_bits[0]) && !clipper::Util::is_nan(cell_bits[1]) && 
	         !clipper::Util::is_nan(cell_bits[2]) && !clipper::Util::is_nan(cell_bits[3]) && 
		 !clipper::Util::is_nan(cell_bits[4]) && !clipper::Util::is_nan(cell_bits[5]) ) { 
								      
	       cell_.init(clipper::Cell_descr(cell_bits[0], cell_bits[1],
					      cell_bits[2], cell_bits[3],
					      cell_bits[4], cell_bits[5]));
	       std::cout << cell_.format() << std::endl;
	       clipper_cell_set_flag = 1;
	    }
	 }
      }
   }

   return
      clipper_cell_set_flag &&
      clipper_symm_set_flag; 
}
 

} // namespace clipper
