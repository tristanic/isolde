/*!  \file function_object_bases.h
  Header file for specifying interfaces to function objects
  \ingroup g_funcobj
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


#ifndef CLIPPER_FUNCTION_OBJECT_BASES
#define CLIPPER_FUNCTION_OBJECT_BASES


#include "../core/hkl_datatypes.h"
#include "../core/xmap.h"
#include "../core/nxmap_operator.h"

namespace clipper {

  //! Base class for structure factor calculation methods
  /*! \ingroup g_funcobj */
  template<class T> class SFcalc_base {
  public:
    //! Structure factor calculation function definition
    /*! In the implementations, this function will do the actual
      structure factor calculation.
      \param fphidata Fourier coefficients to be calculated.
      On output this holds the Fourier coefficients corresponding to
      the atomic model.
      \param atoms The atom selection from which the density is to be
      calculated.
      \return true on success. */
    virtual bool operator() ( HKL_data<datatypes::F_phi<T> >& fphidata, const Atom_list& atoms ) const = 0;
    virtual ~SFcalc_base() {}  //!< destructor
  };


  //! Base class for electron density calculation methods
  /*! \ingroup g_funcobj */
  template<class T> class EDcalc_base {
  public:
    //! Electron density calculation function definition
    /*! In the implementations, this function will do the actual
      electron density calculation.
      \param xmap Electron density map to be calculated.
      On output this holds the electron density map corresponding to
      the atomic model.
      \param atoms The atom selection from which the density is to be
      calculated.
      \return true on success. */
    virtual bool operator() ( Xmap<T>& xmap, const Atom_list& atoms ) const = 0;
    virtual bool operator() ( NXmap<T>& nxmap, const Atom_list& atoms ) const = 0;
    virtual ~EDcalc_base() {}  //!< destructor
  };


  //! Base class for structure factor calculation vs observed methods
  /*! \ingroup g_funcobj */
  template<class T> class SFcalc_obs_base {
  public:
    //! Structure factor calculation function definition
    /*! In the implementations, this function will do the actual
      structure factor calculation.
      \param fphidata Fourier coefficients to be calculated.
      On output this holds the Fourier coefficients corresponding to
      the atomic model.
      \param fo Observed data against which to reference the Fcalc,
      e.g. to fit scale and/or bulk solvent.
      \param atoms The atom selection from which the density is to be
      calculated.
      \return true on success. */
    virtual bool operator() ( HKL_data<datatypes::F_phi<T> >& fphi, const HKL_data<datatypes::F_sigF<T> >& fsig, const Atom_list& atoms ) = 0;
    virtual ~SFcalc_obs_base() {}  //!< destructor
  };


  //! Base class for structure factor weighting (sigmaa) methods
  /*! \ingroup g_funcobj */
  template<class T> class SFweight_base {
  public:
    //! Flag values for different reflection purposes
    enum TYPE { NONE, SIGMAA, SCALE, BOTH };
    //! Structure factor weighting (sigmaa) definition
    /*! In the implementations, this function will do the actual
      structure factor weighting calculation
      \param fb Output best map coefficients.
      \param fd Output difference map coefficients.
      \param phiw Output phase and figure-of-merit.
      \param fo Input observed structure factors.
      \param fc Input calculated map coefficients.
      \param flag Input flag indicating what to use this reflection for.
      \return true on success. */
    virtual bool operator() ( HKL_data<datatypes::F_phi<T> >& fb, HKL_data<datatypes::F_phi<T> >& fd, HKL_data<datatypes::Phi_fom<T> >& phiw, const HKL_data<datatypes::F_sigF<T> >& fo, const HKL_data<datatypes::F_phi<T> >& fc, const HKL_data<datatypes::Flag>& usage ) = 0;
    virtual ~SFweight_base() {}  //!< destructor
  };


  //! Base class for structure factor scaling methods
  /*! \ingroup g_funcobj */
  template<class T> class SFscale_base {
  public:
    //! Structure factor scaling definition
    /*! In the implementations, this function will do the actual
      structure factor scaling calculation.
      \param fo Observed structure factors to be scaled.
      \param fc Calculated structure factors to scale against.
      \return true on success. */
    virtual bool operator() ( HKL_data<datatypes::F_sigF<T> >& fo, const HKL_data<datatypes::F_phi<T> >& fc ) = 0;
    //! Structure factor scaling definition
    /*! In the implementations, this function will do the actual
      structure factor scaling calculation.
      \param fc Calculated structure factors to be scaled.
      \param fo Observed structure factors to scale against.
      \return true on success. */
    virtual bool operator() ( HKL_data<datatypes::F_phi<T> >& fc, const HKL_data<datatypes::F_sigF<T> >& fo ) = 0;
    //! Structure factor scaling definition
    /*! In the implementations, this function will do the actual
      structure factor scaling calculation.
      \param fo Observed structure factors to be scaled.
      \return true on success. */
    virtual bool operator() ( HKL_data<datatypes::F_sigF<T> >& fo ) = 0;
    virtual ~SFscale_base() {}  //!< destructor
  };


  //! Base class for map filtering calculation methods
  /*! \ingroup g_funcobj */
  template<class T> class MapFilter_base {
  public:
    //! Map filtering function definition
    /*! In the implementations, this function will do the actual
      map filtering calculation.
      \param result 
      On output this holds the filtered electron density.
      \param xmap Electron density map to be filtered.
      \return true on success. */
    virtual bool operator() ( Xmap<T>& result, const Xmap<T>& xmap ) const = 0;
    virtual ~MapFilter_base() {}  //!< destructor
  };


  //! Base class for convolution search methods
  /*! \ingroup g_funcobj Note: that target map must be supplied in the
   constructor or initialiser. This allows precalculation of some of
   the coefficients required. The () operator performs an actual
   search in that map for a given target. */
  template<class T> class Convolution_search_base {
  public:
    /*! In the implementations, this function will do the actual
      convolution search calculation.
      \param result The convolution search residual map, containing
      the value of the search function for each offset of the search
      model. The search function is a least squares residual, or may
      be a negative log likelihood, therefore it is positive and low
      values represent better fits.
      \param srchval Map of the target density for which to search.
      \param nxop Operator relating the coordinate frames of the search target and the crystallographic map, used to rotate the search target.
      \return true on success. */
    virtual bool operator() ( Xmap<T>& result, const NXmap<T>& srchval, const NX_operator& nxop ) const = 0;
    virtual ~Convolution_search_base() {}  //!< destructor
  };

  //! Base class for Fast Fourier Feature recongition (FFFEAR) methods
  /*! \ingroup g_funcobj Note: that target map must be supplied in the
   constructor or initialiser. This allows precalculation of some of
   the coefficients required. The () operator performs an actual
   search in that map for a given target. */
  template<class T> class FFFear_base {
  public:
    /*! In the implementations, this function will do the actual
      fffear calculation.
      \param result The fffear residual map, containing the value of
      the search function for each offset of the search model. The
      search function is a least squares residual, or may be a
      negative log likelihood, therefore it is positive and low values
      represent better fits.
      \param srchval Map of the target density for which to search.
      \param srchwgt Map of the target weights for each density point.
      \param nxop Operator relating the coordinate frames of the search target and the crystallographic map, used to rotate the search target.
      \return true on success. */
    virtual bool operator() ( Xmap<T>& result, const NXmap<T>& srchval, const NXmap<T>& srchwgt, const NX_operator& nxop ) const = 0;
    virtual ~FFFear_base() {}  //!< destructor
  };

  //! Base class for skeleton calculation methods
  /*! \ingroup g_funcobj */
  template<class T1, class T2> class Skeleton_base {
  public:
    //! Skeletonisation function definition
    /*! In the implementations, this function will do the actual
      skeletonisation.
      \param xskl int/short/char map.
      On input this may hold 1 for any grid point which is to be
      considered for skeletonisation, and 0 for any other point
      (e.g. low density, solvent)
      [This feaure may not be present in all implementations].
      On output this map holds 1 for any grid point which is part of a
      skeleton ridge, and 0 for any point which is not.
      \param xmap float/double map. This map holds the actual electron
      density values for the map.
      \return true on success. */
    virtual bool operator() ( Xmap<T1>& xskl, const Xmap<T2>& xmap ) const = 0;
    virtual ~Skeleton_base() {}  //!< destructor
  };

  //! Base class for origin matching calculation methods
  /*! \ingroup g_funcobj */
  template<class T> class OriginMatch_base {
  public:
    //! Origin matching function definition
    /*! In the implementations, this function will do the actual
      origin matching.
      \param invert True if the phases must be inverted.
      \param shift The coordinate shift required to match the phases.
      \patam fphi1 The first set of map coefficients.
      \patam fphi2 The second set of map coefficients.
      \return true on success. */
    virtual bool operator() ( bool& invert, Coord_frac& shift, const HKL_data<datatypes::F_phi<T> >& fphi1, const HKL_data<datatypes::F_phi<T> >& fphi2 ) const = 0;
    virtual ~OriginMatch_base() {}  //!< destructor
  };

} // namespace clipper

#endif
