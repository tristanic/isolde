/*! \file cctbx/clipper_cctbx.h
    Header file for cctbx interface
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


#ifndef CLIPPER_CCTBX
#define CLIPPER_CCTBX


#include "../core/hkl_datatypes.h"

#include <cctbx/sgtbx/space_group.h>
#include <cctbx/hendrickson_lattman.h>


namespace clipper
{
  namespace af = scitbx::af;

  //! cctbx interface class
  /*! This class contains static methods which convert various objects
    between Clipper and cctbx formats. See
    http:://cctbx.sourceforge.net/ */
  class CCTBX
  {
  public:
    //! Convert Clipper HKL to cctbx Miller index
    static cctbx::miller::index<> Hkl( const HKL& hkl );
    //! Convert cctbx Miller index to Clipper HKL
    static HKL Hkl( const cctbx::miller::index<>& hkl );
    //! Convert Clipper HKL to cctbx Miller index
    /*! This verision is inlined and uses a reinterpret cast, and so
      is very fast, but depends on the binary representations used in
      the two packages remaining compatible. */
    static inline const cctbx::miller::index<>& hkl( const HKL& hkl )
      { return reinterpret_cast<const cctbx::miller::index<>&>(hkl); }
    //! Convert cctbx Miller index to Clipper HKL
    /*! This verision is inlined and uses a reinterpret cast, and so
      is very fast, but depends on the binary representations used in
      the two packages remaining compatible. */
    static inline const HKL& hkl( const cctbx::miller::index<>& hkl )
      { return reinterpret_cast<const HKL&>(hkl); }

    //! Convert Clipper cell to cctbx cell
    static cctbx::uctbx::unit_cell cell( const Cell& cell );
    //! Convert cctbx cell to Clipper cell
    static Cell cell( const cctbx::uctbx::unit_cell& cell );

    //! Convert Clipper spacegroup to cctbx spacegroup
    static cctbx::sgtbx::space_group spacegroup( const Spacegroup& spgr );
    //! Convert cctbx spacegroup to Clipper spacegroup
    static Spacegroup spacegroup( const cctbx::sgtbx::space_group& spgr );

    //! Conversion of cctbx unit cell and space group.
    static HKL_info
    as_HKL_info(
      cctbx::uctbx::unit_cell const& unit_cell,
      cctbx::sgtbx::space_group const& space_group,
      double d_min,
      double tolerance=1.e-8)
    {
      return HKL_info(
        spacegroup(space_group),
        cell(unit_cell),
        Resolution(d_min-tolerance));
    }

    //! Conversion of cctbx unit cell and space group.
    static HKL_info
    as_HKL_info(
      cctbx::uctbx::unit_cell const& unit_cell,
      cctbx::sgtbx::space_group const& space_group,
      af::const_ref<cctbx::miller::index<> > const& miller_indices,
      double tolerance=1.e-8)
    {
      double max_d_star_sq = unit_cell.max_d_star_sq(miller_indices);
      CCTBX_ASSERT(max_d_star_sq != 0);
      HKL_info result = as_HKL_info(
        unit_cell, space_group, 1/std::sqrt(max_d_star_sq), tolerance);
      std::vector<HKL> hkl_list;
      hkl_list.reserve(miller_indices.size());
      for(std::size_t i=0;i<miller_indices.size();i++) {
        hkl_list.push_back(Hkl(miller_indices[i]));
      }
      result.add_hkl_list(hkl_list);
      return result;
    }

    //! Conversion of F, sigF.
    static HKL_data<data64::F_sigF>
    as_HKL_data(
      HKL_info& hkl_info,
      af::const_ref<cctbx::miller::index<> > const& miller_indices,
      af::const_ref<double> const& data,
      af::const_ref<double> const& sigmas)
    {
      CCTBX_ASSERT(data.size() == miller_indices.size());
      CCTBX_ASSERT(sigmas.size() == miller_indices.size());
      HKL_data<data64::F_sigF> hkl_data(hkl_info);
      for(std::size_t i=0;i<miller_indices.size();i++) {
        data64::F_sigF datum;
        datum.f() = data[i];
        datum.sigf() = sigmas[i];
        CCTBX_ASSERT(
          hkl_data.set_data(Hkl(miller_indices[i]), datum));
      }
      return hkl_data;
    }

    //! Conversion of complex structure factors.
    static HKL_data<data64::F_phi>
    as_HKL_data(
      HKL_info& hkl_info,
      af::const_ref<cctbx::miller::index<> > const& miller_indices,
      af::const_ref<std::complex<double> > const& data)
    {
      CCTBX_ASSERT(data.size() == miller_indices.size());
      HKL_data<data64::F_phi> hkl_data(hkl_info);
      for(std::size_t i=0;i<miller_indices.size();i++) {
        CCTBX_ASSERT(
          hkl_data.set_data(Hkl(miller_indices[i]), data64::F_phi(data[i])));
      }
      return hkl_data;
    }

    //! Conversion of complex structure factors.
    static af::shared<std::complex<double> >
    extract_complex(
      HKL_data<data64::F_phi> const& hkl_data,
      af::const_ref<cctbx::miller::index<> > const& miller_indices)
    {
      af::shared<std::complex<double> >
        result((af::reserve(miller_indices.size())));
      data64::F_phi datum;
      for(std::size_t i=0;i<miller_indices.size();i++) {
        CCTBX_ASSERT(hkl_data.get_data(CCTBX::Hkl(miller_indices[i]), datum));
        CCTBX_ASSERT(!datum.missing());
        result.push_back(datum);
      }
      return result;
    }

    //! Conversion of Hendrickson-Lattman coefficients.
    static HKL_data<data64::ABCD>
    as_HKL_data(
      HKL_info& hkl_info,
      af::const_ref<cctbx::miller::index<> > const& miller_indices,
      af::const_ref<cctbx::hendrickson_lattman<> > const& data)
    {
      CCTBX_ASSERT(data.size() == miller_indices.size());
      HKL_data<data64::ABCD> hkl_data(hkl_info);
      data64::ABCD abcd;
      for(std::size_t i=0;i<miller_indices.size();i++) {
        abcd.data_import(data[i].begin());
        CCTBX_ASSERT(hkl_data.set_data(Hkl(miller_indices[i]), abcd));
      }
      return hkl_data;
    }

    //! Conversion of Hendrickson-Lattman coefficients.
    static af::shared<cctbx::hendrickson_lattman<> >
    extract_hendrickson_lattman(
      HKL_data<data64::ABCD> const& hkl_data,
      af::const_ref<cctbx::miller::index<> > const& miller_indices)
    {
      af::shared<cctbx::hendrickson_lattman<> >
        result((af::reserve(miller_indices.size())));
      data64::ABCD datum;
      for(std::size_t i=0;i<miller_indices.size();i++) {
        CCTBX_ASSERT(hkl_data.get_data(CCTBX::Hkl(miller_indices[i]), datum));
        CCTBX_ASSERT(!datum.missing());
        result.push_back(cctbx::hendrickson_lattman<>(
          datum.a(), datum.b(), datum.c(), datum.d()));
      }
      return result;
    }

    //! Conversion of centroid phases.
    static af::shared<double>
    extract_centroid_phases(
      HKL_data<data64::Phi_fom> const& hkl_data,
      af::const_ref<cctbx::miller::index<> > const& miller_indices)
    {
      af::shared<double>
        result((af::reserve(miller_indices.size())));
      data64::Phi_fom datum;
      for(std::size_t i=0;i<miller_indices.size();i++) {
        CCTBX_ASSERT(hkl_data.get_data(CCTBX::Hkl(miller_indices[i]), datum));
        CCTBX_ASSERT(!datum.missing());
        result.push_back(datum.phi());
      }
      return result;
    }

    //! Conversion of figures of merit.
    static af::shared<double>
    extract_figures_of_merit(
      HKL_data<data64::Phi_fom> const& hkl_data,
      af::const_ref<cctbx::miller::index<> > const& miller_indices)
    {
      af::shared<double>
        result((af::reserve(miller_indices.size())));
      data64::Phi_fom datum;
      for(std::size_t i=0;i<miller_indices.size();i++) {
        CCTBX_ASSERT(hkl_data.get_data(CCTBX::Hkl(miller_indices[i]), datum));
        CCTBX_ASSERT(!datum.missing());
        result.push_back(datum.fom());
      }
      return result;
    }
  };

} // namespace clipper

#endif
