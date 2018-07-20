/*! \file lib/hkl_compute.h
    Data conversion operators for the clipper libraries
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


#ifndef CLIPPER_HKL_COMPUTE
#define CLIPPER_HKL_COMPUTE

#include "hkl_datatypes.h"
#include "../imex.h"

namespace clipper
{

  namespace datatypes
  {
    // Define conversion operators

    //! Compute from ABCD to Phi_fom by phase integration (loses bimodality)
    template<class dtype> class Compute_phifom_from_abcd
    {
    public:
      // constructor: sets up integration tables
      Compute_phifom_from_abcd();
      const Phi_fom<dtype> operator()( const HKL_info::HKL_reference_index& ih, const ABCD<dtype>& abcd ) const;  //!<conv op
    private:
      ftype costab[144], sintab[144];
    };

    //! Compute from Phi_fom to ABCD ( C = D = 0 )
    template<class dtype> class Compute_abcd_from_phifom
    {
    public:
      const ABCD<dtype> operator()( const HKL_info::HKL_reference_index& ih, const Phi_fom<dtype>& phifom ) const;  //!<conv op
    };

    //! Compute from F_sigF+Phi_fom to F_phi
    template<class dtype> class Compute_fphi_from_fsigf_phifom
    {
    public:
      const F_phi<dtype> operator()( const HKL_info::HKL_reference_index& ih, const F_sigF<dtype>& fsigf, const Phi_fom<dtype>& phifom ) const;  //!<conv op
    };

    //! Compute from F_sigF to E_sigE
    template<class dtype> class Compute_EsigE_from_FsigF
    {
    public:
      const E_sigE<dtype> operator()( const HKL_info::HKL_reference_index& ih, const F_sigF<dtype>& fsigf ) const;  //!<conv op
    };

    //! Compute from F_sigF_anom to F_sigF (mean structure factor) 
    template<class dtype> class Compute_mean_fsigf_from_fsigfano
    {
    public:
      const F_sigF<dtype> operator()( const HKL_info::HKL_reference_index& ih, const F_sigF_ano<dtype>& fsigfano ) const;  //!<conv op
    };

    //! Compute from F_sigF_anom to F_sigF (difference structure factor) 
    template<class dtype> class Compute_diff_fsigf_from_fsigfano
    {
    public:
      const F_sigF<dtype> operator()( const HKL_info::HKL_reference_index& ih, const F_sigF_ano<dtype>& fsigfano ) const;  //!<conv op
    };

    //! Negate F_phi (i.e. advance phase by pi)
    template<class dtype> class Compute_neg_fphi
    {
    public:
      const F_phi<dtype> operator()( const HKL_info::HKL_reference_index& ih, const F_phi<dtype>& fphi1 ) const;  //!<conv op
    };

    //! Add two F_phi datalists
    template<class dtype> class Compute_add_fphi
    {
    public:
      const F_phi<dtype> operator()( const HKL_info::HKL_reference_index& ih, const F_phi<dtype>& fphi1, const F_phi<dtype>& fphi2 ) const;  //!<conv op
    };

    //! Subtract two F_phi datalists
    template<class dtype> class Compute_sub_fphi
    {
    public:
      const F_phi<dtype> operator()( const HKL_info::HKL_reference_index& ih, const F_phi<dtype>& fphi1, const F_phi<dtype>& fphi2 ) const;  //!<conv op
    };

    //! Add two ABCD datalists
    template<class dtype> class Compute_add_abcd
    {
    public:
      const ABCD<dtype> operator()( const HKL_info::HKL_reference_index& ih, const ABCD<dtype>& abcd1, const ABCD<dtype>& abcd2 ) const;  //!<conv op
    };

    //! \deprecated Apply scale and U to any scalable datatype
    template<class T> class Compute_scale_u
    {
    public:
      //! constructor: takes scale, U value
      Compute_scale_u( const ftype& s, const ftype& u );
      const T operator()( const HKL_info::HKL_reference_index& ih, T data ) const;  //!<conv op
    private:
      ftype s_, u_;
    };

    //! Apply scale and U to any scalable datatype
    template<class T> class Compute_scale_u_iso
    {
    public:
      //! constructor: takes scale, U value
      Compute_scale_u_iso( const ftype& s, const ftype& u );
      const T operator()( const HKL_info::HKL_reference_index& ih, T data ) const;  //!<conv op
    private:
      ftype s_, u_;
    };

    //! Apply scale and U to any scalable datatype
    template<class T> class Compute_scale_u_aniso
    {
    public:
      //! constructor: takes scale, U value
      Compute_scale_u_aniso( const ftype& s, const U_aniso_orth& u );
      const T operator()( const HKL_info::HKL_reference_index& ih, T data ) const;  //!<conv op
    private:
      ftype s_;
      U_aniso_orth u_;
    };

    //! Compute from F_sigF to F_sigF
    /*! Use this to get F_sigF from F_sigF of a different precision or
      from F_sigF_ano. */
    template<class dtype, class T> class Compute_FsigF
    {
    public:
      const F_sigF<dtype> operator()( const HKL_info::HKL_reference_index& ih, const T& fsigf ) const;  //!<conv op
    };


  }


  namespace data32
  {
    typedef datatypes::Compute_phifom_from_abcd<ftype32> Compute_phifom_from_abcd;  //!< operator
    typedef datatypes::Compute_abcd_from_phifom<ftype32> Compute_abcd_from_phifom;  //!< operator
    typedef datatypes::Compute_fphi_from_fsigf_phifom<ftype32> Compute_fphi_from_fsigf_phifom;  //!< operator
    typedef datatypes::Compute_EsigE_from_FsigF<ftype32> Compute_EsigE_from_FsigF;  //!< operator
    typedef datatypes::Compute_mean_fsigf_from_fsigfano<ftype32> Compute_mean_fsigf_from_fsigfano;  //!< operator
    typedef datatypes::Compute_diff_fsigf_from_fsigfano<ftype32> Compute_diff_fsigf_from_fsigfano;  //!< operator
    typedef datatypes::Compute_neg_fphi<ftype32> Compute_neg_fphi;  //!< operator
    typedef datatypes::Compute_add_fphi<ftype32> Compute_add_fphi;  //!< operator
    typedef datatypes::Compute_sub_fphi<ftype32> Compute_sub_fphi;  //!< operator
    typedef datatypes::Compute_add_abcd<ftype32> Compute_add_abcd;  //!< operator
    typedef datatypes::Compute_scale_u_iso<datatypes::I_sigI<ftype32> > Compute_scale_u_iso_isigi;  //!< operator
    typedef datatypes::Compute_scale_u_iso<datatypes::F_sigF<ftype32> > Compute_scale_u_iso_fsigf;  //!< operator
    typedef datatypes::Compute_scale_u_iso<datatypes::F_sigF_ano<ftype32> > Compute_scale_u_iso_fsigfano;  //!< operator
    typedef datatypes::Compute_scale_u_iso<datatypes::F_phi<ftype32> > Compute_scale_u_iso_fphi;  //!< operator
    typedef datatypes::Compute_scale_u_aniso<datatypes::I_sigI<ftype32> > Compute_scale_u_aniso_isigi;  //!< operator
    typedef datatypes::Compute_scale_u_aniso<datatypes::F_sigF<ftype32> > Compute_scale_u_aniso_fsigf;  //!< operator
    typedef datatypes::Compute_scale_u_aniso<datatypes::F_sigF_ano<ftype32> > Compute_scale_u_aniso_fsigfano;  //!< operator
    typedef datatypes::Compute_scale_u_aniso<datatypes::F_phi<ftype32> > Compute_scale_u_aniso_fphi;  //!< operator
    typedef datatypes::Compute_FsigF<ftype32, datatypes::F_sigF<ftype32> > Compute_FsigF;  //!< operator
    typedef datatypes::Compute_FsigF<ftype32, datatypes::F_sigF_ano<ftype32> > Compute_FsigF_from_ano;  //!< operator

    typedef datatypes::Compute_scale_u<datatypes::I_sigI<ftype32> > Compute_scale_u_isigi;  // DEPRECATED
    typedef datatypes::Compute_scale_u<datatypes::F_sigF<ftype32> > Compute_scale_u_fsigf;  // DEPRECATED
    typedef datatypes::Compute_scale_u<datatypes::F_sigF_ano<ftype32> > Compute_scale_u_fsigfano;  // DEPRECATED
    typedef datatypes::Compute_scale_u<datatypes::F_phi<ftype32> > Compute_scale_u_fphi;  // DEPRECATED
  }

  namespace data64
  {
    typedef datatypes::Compute_phifom_from_abcd<ftype64> Compute_phifom_from_abcd;  //!< operator
    typedef datatypes::Compute_abcd_from_phifom<ftype64> Compute_abcd_from_phifom;  //!< operator
    typedef datatypes::Compute_fphi_from_fsigf_phifom<ftype64> Compute_fphi_from_fsigf_phifom;  //!< operator
    typedef datatypes::Compute_EsigE_from_FsigF<ftype64> Compute_EsigE_from_FsigF;  //!< operator
    typedef datatypes::Compute_mean_fsigf_from_fsigfano<ftype64> Compute_mean_fsigf_from_fsigfano;  //!< operator
    typedef datatypes::Compute_diff_fsigf_from_fsigfano<ftype64> Compute_diff_fsigf_from_fsigfano;  //!< operator
    typedef datatypes::Compute_neg_fphi<ftype64> Compute_neg_fphi;  //!< operator
    typedef datatypes::Compute_add_fphi<ftype64> Compute_add_fphi;  //!< operator
    typedef datatypes::Compute_sub_fphi<ftype64> Compute_sub_fphi;  //!< operator
    typedef datatypes::Compute_add_abcd<ftype64> Compute_add_abcd;  //!< operator
    typedef datatypes::Compute_scale_u_iso<datatypes::I_sigI<ftype64> > Compute_scale_u_iso_isigi;  //!< operator
    typedef datatypes::Compute_scale_u_iso<datatypes::F_sigF<ftype64> > Compute_scale_u_iso_fsigf;  //!< operator
    typedef datatypes::Compute_scale_u_iso<datatypes::F_sigF_ano<ftype64> > Compute_scale_u_iso_fsigfano;  //!< operator
    typedef datatypes::Compute_scale_u_iso<datatypes::F_phi<ftype64> > Compute_scale_u_iso_fphi;  //!< operator
    typedef datatypes::Compute_scale_u_aniso<datatypes::I_sigI<ftype64> > Compute_scale_u_aniso_isigi;  //!< operator
    typedef datatypes::Compute_scale_u_aniso<datatypes::F_sigF<ftype64> > Compute_scale_u_aniso_fsigf;  //!< operator
    typedef datatypes::Compute_scale_u_aniso<datatypes::F_sigF_ano<ftype64> > Compute_scale_u_aniso_fsigfano;  //!< operator
    typedef datatypes::Compute_scale_u_aniso<datatypes::F_phi<ftype64> > Compute_scale_u_aniso_fphi;  //!< operator
    typedef datatypes::Compute_FsigF<ftype64, datatypes::F_sigF<ftype64> > Compute_FsigF;  //!< operator
    typedef datatypes::Compute_FsigF<ftype64, datatypes::F_sigF_ano<ftype64> > Compute_FsigF_from_ano;  //!< operator

    typedef datatypes::Compute_scale_u<datatypes::I_sigI<ftype64> > Compute_scale_u_isigi;  // DEPRECATED
    typedef datatypes::Compute_scale_u<datatypes::F_sigF<ftype64> > Compute_scale_u_fsigf;  // DEPRECATED
    typedef datatypes::Compute_scale_u<datatypes::F_sigF_ano<ftype64> > Compute_scale_u_fsigfano;  // DEPRECATED
    typedef datatypes::Compute_scale_u<datatypes::F_phi<ftype64> > Compute_scale_u_fphi;  // DEPRECATED
  }

} // namespace clipper

#endif
