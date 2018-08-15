#pragma once
#include <string>
#include <unordered_map>
#include <future>
#include <chrono>

#include "imex.h"
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

#include "sfcalc_obs_vdw.h"

using namespace clipper;
using namespace clipper::datatypes;
namespace clipper_cx { // ChimeraX extensions to clipper

// Container for a crystallographic map and its coefficients, including the
// details needed to regenerate it. Just stores everything - calculations are
// handled by Xtal_mgr_base
struct CLIPPER_CX_IMEX Xmap_details
{
public:
    inline Xmap_details(const HKL_data<F_phi<ftype32>>& base_coeffs, const ftype& b_sharp,
        const Grid_sampling& grid_sampling, bool is_difference_map=false,
        bool exclude_free_reflections=true,
        bool fill_with_fcalc=true)
        : base_coeffs_(base_coeffs), b_sharp_(b_sharp),
          is_difference_map_(is_difference_map),
          exclude_freer_(exclude_free_reflections),
          fill_(fill_with_fcalc)
    {
        // initialise empty data objects
        coeffs_ = HKL_data<F_phi<ftype32>>(base_coeffs.hkl_info());
        xmap_ = Xmap<ftype32>(base_coeffs.spacegroup(), base_coeffs.cell(), grid_sampling);
    }


    inline const Xmap<ftype32>& xmap() const { return xmap_; }
    inline Xmap<ftype32>& xmap() { return xmap_; }
    inline const Map_stats& map_stats() const { return map_stats_; }
    inline Map_stats& map_stats() { return map_stats_; }
    inline const HKL_data<F_phi<ftype32>>& coeffs() const { return coeffs_; }
    inline HKL_data<F_phi<ftype32>>& coeffs() { return coeffs_; }
    inline const HKL_data<F_phi<ftype32>>& base_coeffs() { return base_coeffs_; }
    inline const ftype& b_sharp() const { return b_sharp_; }
    inline bool is_difference_map() const { return is_difference_map_; }
    inline bool exclude_free_reflections() const { return exclude_freer_; }
    inline bool fill_with_fcalc() const { return fill_; }


private:
    // Base coefficients before removing free reflections, sharpening etc.
    const HKL_data<F_phi<ftype32>>& base_coeffs_;
    // Final coefficients used to generate xmap_
    ftype b_sharp_=0;
    // If this is a difference map, fill_ is ignored
    bool is_difference_map_ = false;
    // Exclude free reflections from map?
    bool exclude_freer_=true;
    // if exclude_freer_, replace free reflections with DFcalc?
    bool fill_=true;

    HKL_data<F_phi<ftype32>> coeffs_;
    // Sharpening/smoothing B-factor to apply
    // The real-space map
    Xmap<ftype32> xmap_;

    // max, min, mean, std_dev, range
    Map_stats map_stats_;

}; // class Xmap_details



//! General manager class for handling operations on a set of crystallographic
//  maps and their data
class CLIPPER_CX_IMEX Xtal_mgr_base
{
public:
    Xtal_mgr_base() {} // default constructor
    Xtal_mgr_base(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
        const Grid_sampling& grid_sampling, const HKL_data<F_sigF<ftype32>>& fobs);

    inline bool fcalc_initialized() const { return fcalc_initialized_; }
    inline bool coeffs_initialized() const { return coeffs_initialized_; }
    inline int freeflag() const { return freeflag_; }
    void set_freeflag(int f);

    inline const ftype& rwork() { return rwork_; }
    inline const ftype& rfree() { return rfree_; }
    inline const ftype& weighted_rwork() { return w_rwork_; }
    inline const ftype& weighted_rfree() { return w_rfree_; }

    inline const ftype& bulk_frac()
    {
        if (!fcalc_initialized())
            throw std::runtime_error("Coefficients have not been calculated yet!");
        return bulk_solvent_calculator_.bulk_frac();
    }
    inline const ftype& bulk_scale()
    {
        if (!fcalc_initialized())
            throw std::runtime_error("Coefficients have not been calculated yet!");
        return bulk_solvent_calculator_.bulk_scale();
    }

    inline const HKL_data<F_sigF<ftype32>>& fobs() const
    {
        return fobs_;
    }
    inline const HKL_data<F_phi<ftype32>>& fcalc() const
    {
        if (!fcalc_initialized())
            throw std::runtime_error("Coefficients have not been calculated yet!");
        return fcalc_;
    }

    // Fcalc scaled to Fobs
    HKL_data<F_phi<ftype32>> scaled_fcalc() const;

    inline const HKL_data<Flag>& usage() const { return usage_; }
    inline const HKL_data<F_phi<ftype32>>& base_2fofc() const
    {
        if (!coeffs_initialized())
            throw std::runtime_error("Coefficients have not yet been calculated!");
        return base_2fofc_;
    }
    inline const HKL_data<F_phi<ftype32>>& base_fofc() const
    {
        if (!coeffs_initialized())
            throw std::runtime_error("Coefficients have not yet been calculated!");
        return base_fofc_;
    }

    inline const HKL_data<Phi_fom<ftype32>>& weights() const
    {
        if (!coeffs_initialized())
            throw std::runtime_error("Weights have not yet been calculated!");
        return phi_fom_;
    }

    int guess_free_flag_value(const HKL_data<Flag>& flags);

    // Regenerate Fcalc from a set of ChimeraX Atoms
    void generate_fcalc(const Atom_list& atoms);

    // Generate the standard set of map coefficients
    void generate_base_map_coeffs();

    // Calculate Rwork and Rfree. Called automatically by generate_fcalc()
    void calculate_r_factors();

    // Apply in-place B-factor sharpening to a set of map coefficients
    void apply_b_factor_sharpening(HKL_data<F_phi<ftype32>>& coeffs, const ftype& bsharp);

    void add_xmap(const std::string& name, const HKL_data<F_phi<ftype32>>& base_coeffs,
        const ftype& bsharp, bool is_difference_map=false,
        bool exclude_free_reflections=true, bool fill_with_fcalc=true);

    void delete_xmap(const std::string& name) { maps_.erase(name); }

    // Recalculate map coefficients and real-space map for one previously-stored
    // set of parameters (will throw std::out_of_range if nothing has been stored
    // under the given name)
    void recalculate_map(const std::string& name);

    void recalculate_map(Xmap_details& xmd);

    // Generate fresh Fcalc, and regenerate all existing maps
    void recalculate_all(const Atom_list& atoms);

    //
    inline Xmap<ftype32>& get_xmap(const std::string& name) { return maps_.at(name).xmap(); }


protected:
    // Safety
    // Have we generated any fcalcs yet?
    bool fcalc_initialized_=false;
    // Have we generated map coefficients yet?
    bool coeffs_initialized_=false;

private:
    // constants
    //static constexpr ftype ONE_1_ON_4_PI_SQUARED = 1/(4*M_PI*M_PI);
    int freeflag_ = 0;

    // Basic information
    HKL_info hklinfo_;
    Cell cell_;
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

    std::unordered_map<std::string, Xmap_details> maps_;
    // Most recent Fcalc values
    HKL_data<F_phi<ftype32>> fcalc_;

    SFcalc_obs_bulk_vdw<ftype32> bulk_solvent_calculator_ = SFcalc_obs_bulk_vdw<ftype32>();

    // SFcalc_obs_bulk<ftype32> bulk_solvent_calculator_ = SFcalc_obs_bulk<ftype32>();

    // Calculates sigmaa-weighted maximum-likelihood map coefficients
    SFweight_spline<ftype32> map_calculator_ = SFweight_spline<ftype32>();

    // "Straight" or standard rwork/rfree
    ftype rwork_ = 1;
    ftype rfree_ = 1;

    // sigma-weighted rwork/rfree
    ftype w_rwork_ = 1;
    ftype w_rfree_ = 1;

    void set_map_free_terms_to_zero(const HKL_data<F_phi<ftype32>>& source,
        HKL_data<F_phi<ftype32>>& dest);
    void set_map_free_terms_to_dfc(const HKL_data<F_phi<ftype32>>& source,
        HKL_data<F_phi<ftype32>>& dest);


}; // class Xmap_mgr


//! Parallel handler for map calculations
/*! Uses std::async to push the task of (re)generating maps to a separate thread
    from the main program. Further splits up map calculations into individual
    threads if available.
*/
class CLIPPER_CX_IMEX Xtal_mgr_thread
{
public:
    Xtal_mgr_thread(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
        const Grid_sampling& grid_sampling, const HKL_data<F_sigF<ftype32>>& fobs,
        const size_t num_threads = 1);

    inline size_t num_threads() const { return num_threads_; }
    inline void set_num_threads(size_t n) const { num_threads_=std::max(n, 1); }

    bool thread_running() const { return master_thread_result_.valid(); }
    bool ready() const {
        if (!thread_running()) return false;
        // std::future does not yet have an is_ready() function. Current
        // workaround is to call wait_for with zero time and query the result.
        auto result = master_thread_result_.wait_for(
            std::chrono::duration<std::chrono::microseconds>(0))
        return (result == std::future_status::ready);
    }



private:
    Xtal_mgr_base mgr_;
    size_t num_threads_;
    std::future<void> master_thread_result_;
    std::vector<std::future<Xmap_details>> xmap_thread_results_;

}


} // namespace clipper_cx
