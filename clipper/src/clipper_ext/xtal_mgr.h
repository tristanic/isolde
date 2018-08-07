#pragma once
#include <string>
#include <unordered_map>

#include "imex.h"
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

using namespace clipper;
using namespace clipper::datatypes;
namespace clipper_cx { // ChimeraX extensions to clipper

//! General manager class for handling operations on a set of crystallographic
//  maps and their data
class CLIPPER_CX_IMEX Xtal_mgr_base
{
public:
    Xtal_mgr_base() {} // default constructor
    Xtal_mgr_base(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
        const Grid_sampling& grid_sampling, const HKL_data<F_sigF<ftype32>>& fobs);

    bool fcalc_initialized() const { return fcalc_initialized_; }
    bool coeffs_initialized() const { return coeffs_initialized_; }
    int freeflag() const { return freeflag_; }
    void set_freeflag(int f);

    const HKL_data<F_phi<ftype32>>& fcalc() const
    {
        if (!fcalc_initialized())
            throw std::runtime_error("Coefficients have not been calculated yet!");
        return fcalc_;
    }

    const HKL_data<Flag>& usage() const { return usage_; }
    const HKL_data<F_phi<ftype32>>& base_2fofc() const
    {
        if (!coeffs_initialized())
            throw std::runtime_error("Coefficients have not yet been calculated!");
        return base_2fofc_;
    }
    const HKL_data<F_phi<ftype32>>& base_fofc() const
    {
        if (!coeffs_initialized())
            throw std::runtime_error("Coefficients have not yet been calculated!");
        return base_fofc_;
    }

    const HKL_data<Phi_fom<ftype32>>& weights() const
    {
        if (!coeffs_initialized())
            throw std::runtime_error("Weights have not yet been calculated!");
        return phi_fom_;
    }


    int guess_free_flag_value(const HKL_data<Flag>& flags);

    // Regenerate Fcalc from a set of ChimeraX Atoms
    void generate_fcalc(const Atom_list& atoms);

    // Generate the standard set of map coefficients
    void generate_base_map_coeffs(bool exclude_free_reflections=true,
        bool fill_with_fcalc=true);

    void set_map_free_terms_to_zero();
    void set_map_free_terms_to_dfc();

protected:
    // Safety
    // Have we generated any fcalcs yet?
    bool fcalc_initialized_=false;
    // Have we generated map coefficients yet?
    bool coeffs_initialized_=false;

private:
    int freeflag_ = 0;

    // Basic information
    HKL_info hklinfo_;
    // Original free-r flags from the reflection data file
    HKL_data<Flag> free_flags_;
    // Flags defining which reflections to use in map calculations
    HKL_data<Flag> usage_;
    Grid_sampling grid_sampling_;
    // Observed structure factors
    HKL_data<F_sigF<ftype32>> fobs_;
    // Real-space maps

    // Base map coefficients before filling, sharpening etc.
    HKL_data<F_phi<ftype32>> base_2fofc_;
    HKL_data<F_phi<ftype32>> base_fofc_;

    // Weighting coefficients
    HKL_data<Phi_fom<ftype32>> phi_fom_;

    std::unordered_map<std::string, Xmap<ftype32>> maps_;
    // Map coefficients used to generate each Xmap
    std::unordered_map<std::string, HKL_data<F_phi<ftype32>>> map_coeffs_;
    // Most recent Fcalc values
    HKL_data<F_phi<ftype32>> fcalc_;

    // Calculates sigmaa-weighted maximum-likelihood map coefficients
    SFweight_spline<ftype32> map_calculator_ = SFweight_spline<ftype32>();
}; // class Xmap_mgr

} // namespace clipper_cx
