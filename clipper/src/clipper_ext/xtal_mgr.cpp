#include "xtal_mgr.h"

#include <set>

namespace clipper_cx
{

Xtal_mgr_base::Xtal_mgr_base(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
    const Grid_sampling& grid_sampling, const HKL_data<F_sigF<ftype32>>& fobs)
    : hklinfo_(hklinfo), free_flags_(free_flags), grid_sampling_(grid_sampling), fobs_(fobs)
{
    fcalc_ = HKL_data<F_phi<ftype32>>(hklinfo);
    base_2fofc_ = HKL_data<F_phi<ftype32>>(hklinfo);
    base_fofc_ = HKL_data<F_phi<ftype32>>(hklinfo);
    phi_fom_ = HKL_data<Phi_fom<ftype32>>(hklinfo);
    usage_ = HKL_data<Flag>(hklinfo);

    set_freeflag(guess_free_flag_value(free_flags));
}


int
Xtal_mgr_base::guess_free_flag_value(const HKL_data<Flag>& flags)
{
    HKL_info::HKL_reference_index ix;
    int f_min = 1e-6, f_max=-1e6;
    std::set<int> flag_vals;
    for (ix=flags.first(); !ix.last(); ix.next())
    {
        const auto &f = flags[ix];
        if (!f.missing())
        {
            if (f.flag() < f_min) f_min = f.flag();
            else if (f.flag() > f_max) f_max = f.flag();
            flag_vals.insert(f.flag());
        }
    }
    if (f_max > 1) /* assume CCP4 */ return 0;
    if (f_min == -1 && flag_vals.size()==2) /* assume SHELX */ return -1;
    else
    {
        // Count the different values, and choose the one that makes most sense
        std::unordered_map<int, int> val_counts;
        for (const auto& v: flag_vals)
            val_counts[v] = 0;
        for (ix=flags.first(); !ix.last(); ix.next())
            val_counts[flags[ix].flag()] += 1;

        std::unordered_map<int, float> val_fracs;
        // Convert to fractions of total reflections
        for (const auto v: flag_vals)
            val_fracs[v] = ((float)val_counts[v]) / flags.hkl_info().num_reflections();

        // Weed out those that don't fit the criteria
        std::set<int> candidates;
        int last_candidate;
        for (auto v: flag_vals)
        {
            if (val_counts[v] < 5000)
                if (val_fracs[v] > 0.005 && val_fracs[v] < 0.1)
                {
                    candidates.insert(v);
                    last_candidate = v;
                }
        }
        if (candidates.size() != 1) {
            throw std::runtime_error("Cannot determine the correct value for the "
            "free reflections in the given FreeR_flag array. Please convert your "
            "reflection file to CCP4 or PHENIX format using sftools or "
            "phenix.reflection_file_editor. ");
            return -1;
        }
        return last_candidate;
    }
} // guess_free_flag_value

void
Xtal_mgr_base::set_freeflag(int f)
{
    freeflag_ = f;
    HKL_info::HKL_reference_index ix;
    for (ix = usage_.first(); !ix.last(); ix.next())
    {
        auto f = free_flags_[ix];
        auto fobs = fobs_[ix];
        if (!f.missing() && !fobs.missing() && !f.flag()==freeflag_)
            usage_[ix].flag() = SFweight_spline<ftype32>::BOTH;
        else
            usage_[ix].flag() = SFweight_spline<ftype32>::NONE;
    }
}

void
Xtal_mgr_base::generate_fcalc(const Atom_list& atoms)
{
    SFcalc_obs_bulk<ftype32>(fcalc_, fobs_, atoms);
    fcalc_initialized_ = true;
}

// Generate the standard set of map coefficients
void
Xtal_mgr_base::generate_base_map_coeffs(bool exclude_free_reflections,
    bool fill_with_fcalc)
{
    if (!fcalc_initialized())
        throw std::runtime_error("No Fcalc values have been calculated! Run "
            " generate_fcalc() on a suitable set of atoms first!");
    map_calculator_(base_2fofc_, base_fofc_, phi_fom_, fobs_, fcalc_, usage_);
    if (exclude_free_reflections)
    {
        if (fill_with_fcalc)
            set_map_free_terms_to_dfc();
        else
            set_map_free_terms_to_zero();
    }

    coeffs_initialized_=true;
}

void
Xtal_mgr_base::set_map_free_terms_to_zero()
{
    if (!coeffs_initialized_())
        throw std::runtime_error("Coefficients have not yet been calculated!");
    HKL_info::HKL_reference_index ix;
    for (ix = base_2fofc_.first(); !ix.last(); ix.next())
    {
        if (usage_[ix].missing() || usage[ix].flag() = SFweight_spline<ftype32>::NONE)
        {
            if (!base_2fofc_[ix].missing())
                base_2fofc_[ix].f() = 0;
            if (!base_fofc_[ix].missing())
                base_fofc_[ix].f() = 0;
        }
    }
}

void
Xtal_mgr_base::set_map_free_terms_to_dfc()
{
    if (!coeffs_initialized_())
        throw std::runtime_error("Coefficients have not yet been calculated!");
    HKL_info::HKL_reference_index ix;
    const auto& param_s = map_calculator_.params_scale();
    // const auto& param_w = map_calculator_.params_error();
    BasisFn_spline basisfn( flag, param_s.size(), 1.0);
    for (ix=base_2fofc_.first(); !ix.last(); ix.next())
    {
        if (usage_[ix].missing() || usage[ix].flag() = SFweight_spline<ftype32>::NONE)
        {
            auto& fpo = base_2fofc_[ix];
            if (!fpo.missing())
            {
                auto s = basisfn.f_s( ix.invresolsq(), param_s);
                fpo.f() = s*fcalc_[ix].f();
            }
            if (!base_fofc_[ix].missing())
                base_fofc_[ix].f() = 0;
        }
    }
}

} //namespace clipper_cx;
