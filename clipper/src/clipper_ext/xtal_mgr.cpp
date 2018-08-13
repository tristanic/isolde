#include "xtal_mgr.h"

#include <set>

namespace clipper_cx
{

Xtal_mgr_base::Xtal_mgr_base(const HKL_info& hklinfo, const HKL_data<Flag>& free_flags,
    const Grid_sampling& grid_sampling, const HKL_data<F_sigF<ftype32>>& fobs)
    : hklinfo_(hklinfo), free_flags_(free_flags), grid_sampling_(grid_sampling), fobs_(fobs)
{
    cell_ = fobs.cell();
    fcalc_ = HKL_data<F_phi<ftype32>>(hklinfo_);
    base_2fofc_ = HKL_data<F_phi<ftype32>>(hklinfo_);
    base_fofc_ = HKL_data<F_phi<ftype32>>(hklinfo_);
    phi_fom_ = HKL_data<Phi_fom<ftype32>>(hklinfo_);
    usage_ = HKL_data<Flag>(hklinfo_);

    set_freeflag(guess_free_flag_value(free_flags));
}


int
Xtal_mgr_base::guess_free_flag_value(const HKL_data<Flag>& flags)
{
    HKL_info::HKL_reference_index ih;
    int f_min = 1e6, f_max=-1e6;
    std::set<int> flag_vals;
    for (ih=flags.first(); !ih.last(); ih.next())
    {
        const auto &f = flags[ih];
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
        for (ih=flags.first(); !ih.last(); ih.next())
            val_counts[flags[ih].flag()] += 1;

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
    HKL_info::HKL_reference_index ih;
    for (ih = usage_.first(); !ih.last(); ih.next())
    {
        auto f = free_flags_[ih];
        auto fobs = fobs_[ih];
        if (!f.missing() && !fobs.missing() && !(f.flag()==freeflag_))
            usage_[ih].flag() = SFweight_spline<ftype32>::BOTH;
        else
            usage_[ih].flag() = SFweight_spline<ftype32>::NONE;
    }
}

void
Xtal_mgr_base::generate_fcalc(const Atom_list& atoms)
{
    bulk_solvent_calculator_(fcalc_, fobs_, atoms);
    fcalc_initialized_ = true;
    calculate_r_factors();
}

// Generate the standard set of map coefficients
void
Xtal_mgr_base::generate_base_map_coeffs()
{
    if (!fcalc_initialized())
        throw std::runtime_error("No Fcalc values have been calculated! Run "
            " generate_fcalc() on a suitable set of atoms first!");
    map_calculator_(base_2fofc_, base_fofc_, phi_fom_, fobs_, fcalc_, usage_);

    coeffs_initialized_=true;
} // generate_base_map_coeffs

HKL_data<F_phi<ftype32>>
Xtal_mgr_base::scaled_fcalc() const
{
    int nparams = 12;
    std::vector<ftype> params(nparams, 1.0);
    HKL_data<F_phi<ftype32>> ret(hklinfo_);
    BasisFn_aniso_gaussian basisfn;
    TargetFn_scaleF1F2<F_phi<ftype32>, F_sigF<ftype32>> targetfn (fcalc_, fobs_);
    ResolutionFn rfn(hklinfo_, basisfn, targetfn, params);
    HKL_info::HKL_reference_index ih;
    for (ih=fcalc_.first(); !ih.last(); ih.next())
    {
        if (!fobs_[ih].missing())
            ret[ih].f() = sqrt(rfn.f(ih))*fcalc_[ih].f();
            ret[ih].phi() = fcalc[ih].phi();
    }
    return ret;

}


void
Xtal_mgr_base::calculate_r_factors()
{
    if (!fcalc_initialized())
        throw std::runtime_error("No Fcalc values have been calculated! Run "
            " generate_fcalc() on a suitable set of atoms first!");

    int nparams = 12;
    std::vector<ftype> params(nparams, 1.0);
    // BasisFn_spline basisfn(hklinfo_, nparams, 1.0);
    BasisFn_aniso_gaussian basisfn;
    TargetFn_scaleF1F2<F_phi<ftype32>, F_sigF<ftype32>> targetfn (fcalc_, fobs_);
    ResolutionFn rfn(hklinfo_, basisfn, targetfn, params);

    HKL_info::HKL_reference_index ih;
    // for standard rwork, rfree
    ftype sum_fwork=0, sum_ffree=0, sum_dwork=0, sum_dfree=0;
    // for sigma-weighted rwork, rfree
    ftype sum_wfwork2=0, sum_wffree2=0, sum_wdwork2=0, sum_wdfree2=0;
    for (ih=fcalc_.first(); !ih.last(); ih.next())
    {
        const auto& fo = fobs_[ih];
        const auto& fc = fcalc_[ih];
        const auto& fflag = free_flags_[ih];
        if (!fo.missing())
        {
            ftype eps = ih.hkl_class().epsilon();
            ftype two_on_eps = 2.0/eps;
            auto scaled_fcalc = sqrt(rfn.f(ih))*fc.f();
            if (fflag.flag()==freeflag_) {
                sum_ffree+=two_on_eps*fo.f();
                sum_dfree+=two_on_eps*std::abs(fo.f()-scaled_fcalc);

                sum_wffree2 += 1/fo.sigf()*pow(two_on_eps*fo.f(), 2);
                sum_wdfree2 += 1/fo.sigf()*pow(two_on_eps*(fo.f()-scaled_fcalc), 2);
            } else {
                sum_fwork+=two_on_eps*fo.f();
                sum_dwork+=two_on_eps*std::abs(fo.f()-scaled_fcalc);

                sum_wfwork2 += 1/fo.sigf()*pow(two_on_eps*fo.f(), 2);
                sum_wdwork2 += 1/fo.sigf()*pow(two_on_eps*(fo.f()-scaled_fcalc), 2);
            }
        }
    }
    rfree_ = sum_dfree/sum_ffree;
    rwork_ = sum_dwork/sum_fwork;

    w_rfree_ = sqrt(sum_wdfree2/sum_wffree2);
    w_rwork_ = sqrt(sum_wdwork2/sum_wfwork2);
} // calculate_r_factors


void
Xtal_mgr_base::apply_b_factor_sharpening(HKL_data<F_phi<ftype32>>& coeffs,
    const ftype& bsharp)
{
    HKL_info::HKL_reference_index ih;
    for (ih=coeffs.first(); !ih.last(); ih.next())
    {
        auto& fphi = coeffs[ih];
        if (!fphi.missing())
            fphi.f()*=exp(-bsharp*ONE_1_ON_4_PI_SQUARED);
    }
}

void
Xtal_mgr_base::add_xmap(const std::string& name, const HKL_data<F_phi<ftype32>>& base_coeffs,
    const ftype& bsharp, bool is_difference_map,
    bool exclude_free_reflections, bool fill_with_fcalc)
{
    maps_.emplace(name, Xmap_details(base_coeffs, bsharp, grid_sampling_, is_difference_map,
                              exclude_free_reflections, fill_with_fcalc));
    recalculate_map(name);
} // add_xmap

void
Xtal_mgr_base::recalculate_map(const std::string& name)
{
    Xmap_details& xmd = maps_.at(name);
    if (xmd.exclude_free_reflections()) {
        if (xmd.fill_with_fcalc() && !xmd.is_difference_map())
            set_map_free_terms_to_dfc(xmd.base_coeffs(), xmd.coeffs());
        else
            set_map_free_terms_to_zero(xmd.base_coeffs(), xmd.coeffs());
    }
    if (xmd.b_sharp() != 0)
        apply_b_factor_sharpening(xmd.coeffs(), xmd.b_sharp());
    xmd.xmap().fft_from(xmd.coeffs());
    xmd.map_stats() = Map_stats(xmd.xmap());
} // recalculate_map


void
Xtal_mgr_base::set_map_free_terms_to_zero(const HKL_data<F_phi<ftype32>>& source,
    HKL_data<F_phi<ftype32>>& dest)
{
    if (!coeffs_initialized())
        throw std::runtime_error("Coefficients have not yet been calculated!");
    HKL_info::HKL_reference_index ih;
    for (ih = source.first(); !ih.last(); ih.next())
    {
        if (usage_[ih].missing() || usage_[ih].flag() == SFweight_spline<ftype32>::NONE)
        {
            if (!source[ih].missing())
                dest[ih].f() = 0;
            else
                dest[ih].set_null();
        } else {
            dest[ih] = source[ih];
        }
    }
}

void
Xtal_mgr_base::set_map_free_terms_to_dfc(const HKL_data<F_phi<ftype32>>& source,
    HKL_data<F_phi<ftype32>>& dest)
{
    if (!coeffs_initialized())
        throw std::runtime_error("Coefficients have not yet been calculated!");
    HKL_info::HKL_reference_index ih;
    HKL_data<Flag_bool> flag(source.hkl_info());
    for (ih=flag.first(); !ih.last(); ih.next())
        flag[ih].flag() = (!source[ih].missing() && usage_[ih].flag()!=SFweight_base<ftype32>::NONE);
    const auto& param_s = map_calculator_.params_scale();
    // const auto& param_w = map_calculator_.params_error();
    BasisFn_spline basisfn( flag, param_s.size(), 1.0);
    for (ih=source.first(); !ih.last(); ih.next())
    {
        const auto& fpo = source[ih];
        if (usage_[ih].missing() || usage_[ih].flag() == SFweight_spline<ftype32>::NONE)
        {
            if (!fpo.missing()) {
                auto s = basisfn.f_s( ih.invresolsq(), param_s);
                dest[ih].f() = s*fcalc_[ih].f();
            } else {
                dest[ih].set_null();
            }
        } else {
            dest[ih] = fpo;
        }
    }
}

} //namespace clipper_cx;
