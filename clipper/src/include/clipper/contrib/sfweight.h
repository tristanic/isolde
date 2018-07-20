/*!  \file sfweight.h
  Header file for structure factor weighting object
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


#ifndef CLIPPER_SFWEIGHT
#define CLIPPER_SFWEIGHT


#include "function_object_bases.h"
#include "../core/resol_basisfn.h"

#include "../imex.h"

namespace clipper {


  //! Structure factor weighting by sigmaa-related spline method.
  /*! Perform structure factor weighting to obtain likelihood
    weights for structure factors.

    This implementation uses a single list of reflections for both
    scaling and sigmaa, thus the only relevent usage flags are
    NONE/BOTH.

    The number of spline parameters or the number of reflections per
    parameter may be specified. If either is zero, the other takes
    priority. If both are non-zero, a compromise value is used.
    \ingroup g_funcobj */
  template<class T> class SFweight_spline : public SFweight_base<T> {
  public:
    //! constructor
    SFweight_spline( const int n_reflns = 1000, const int n_params = 20, const int n_phases = 24 ) { init( n_reflns, n_params, n_phases ); }
    //! constructor: shorthand for constructor+operator
    SFweight_spline( HKL_data<datatypes::F_phi<T> >& fb, HKL_data<datatypes::F_phi<T> >& fd, HKL_data<datatypes::Phi_fom<T> >& phiw, const HKL_data<datatypes::F_sigF<T> >& fo, const HKL_data<datatypes::F_phi<T> >& fc, const HKL_data<datatypes::Flag>& usage, const int n_reflns = 1000, const int n_params = 20 );
    //! initialise: from parameters
    void init( const int n_reflns = 1000, const int n_params = 20, const int n_phases = 24 );

    // Likelihood weighting and map coefficient computation
    bool operator() ( HKL_data<datatypes::F_phi<T> >& fb, HKL_data<datatypes::F_phi<T> >& fd, HKL_data<datatypes::Phi_fom<T> >& phiw, const HKL_data<datatypes::F_sigF<T> >& fo0, const HKL_data<datatypes::F_phi<T> >& fc0, const HKL_data<datatypes::Flag>& usage );
    // Likelihood weighting and map coefficient computation (MLHL)
    bool operator() ( HKL_data<datatypes::F_phi<T> >& fb, HKL_data<datatypes::F_phi<T> >& fd, HKL_data<datatypes::Phi_fom<T> >& phiw, HKL_data<datatypes::ABCD<T> >& hl, const HKL_data<datatypes::F_sigF<T> >& fo0, const HKL_data<datatypes::ABCD<T> >& hl0, const HKL_data<datatypes::F_phi<T> >& fc0, const HKL_data<datatypes::Flag>& usage );

    template<class F> bool evaluate( HKL_data<datatypes::F_phi<T> >& fb, HKL_data<datatypes::F_phi<T> >& fd, HKL_data<datatypes::Phi_fom<T> >& phiw, HKL_data<datatypes::ABCD<T> >& hl, const HKL_data<datatypes::F_sigF<T> >& fo0, const HKL_data<datatypes::ABCD<T> >& hl0, const HKL_data<datatypes::F_phi<T> >& fc0, const HKL_data<datatypes::Flag>& usage, F tgt );
    template<class F> bool reevaluate( HKL_data<datatypes::F_phi<T> >& fb, HKL_data<datatypes::F_phi<T> >& fd, HKL_data<datatypes::Phi_fom<T> >& phiw, HKL_data<datatypes::ABCD<T> >& hl, const HKL_data<datatypes::F_sigF<T> >& fo0, const HKL_data<datatypes::ABCD<T> >& hl0, const HKL_data<datatypes::F_phi<T> >& fc0, const HKL_data<datatypes::Flag>& usage, F tgt );
    const std::vector<ftype>& params_scale() { return param_s; }
    const std::vector<ftype>& params_error() { return param_w; }
    const double& log_likelihood_work() { return llw; }
    const double& log_likelihood_free() { return llf; }

    void debug() const;

    struct TargetResult { ftype r, ds, dw, dss, dww, dsw; };
    TargetResult targetfn( const HKL_class cls, const datatypes::F_sigF<T>& fo0, const datatypes::F_phi<T>& fc0, const ftype& s, const ftype& w ) const;
    TargetResult targethl( const HKL_class cls, const datatypes::F_sigF<T>& fo0, const datatypes::ABCD<T>& hl0, const datatypes::F_phi<T>& fc0, const ftype& s, const ftype& w ) const;
  private:
    struct HLterms { ftype cosa, sina, cos2a, sin2a; };
    class TargetFo {
    public:
      TargetResult operator() ( const HKL_class cls, const datatypes::F_sigF<T>& fo0, const datatypes::ABCD<T>& hl0, const datatypes::F_phi<T>& fc0, const ftype& s, const ftype& w, const std::vector<HLterms>& hlterms );
      datatypes::ABCD<T>    abcd;
      datatypes::Phi_fom<T> phiw;
    };
    class TargetHL {
    public:
      TargetResult operator() ( const HKL_class cls, const datatypes::F_sigF<T>& fo0, const datatypes::ABCD<T>& hl0, const datatypes::F_phi<T>& fc0, const ftype& s, const ftype& w, const std::vector<HLterms>& hlterms );
      datatypes::ABCD<T>    abcd;
      datatypes::Phi_fom<T> phiw;
    };
    int num_params( const HKL_data<datatypes::Flag_bool>& flag ) const;
    int nparams, nreflns;
    std::vector<ftype> param_s, param_w;
    std::vector<T> scale_fo, scale_fc, value_s, value_w;
    std::vector<HLterms> hlterms;
    double llw, llf;
  };



} // namespace clipper

#endif
