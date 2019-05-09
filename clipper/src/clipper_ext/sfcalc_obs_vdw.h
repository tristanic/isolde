/**
 * @Author: Tristan Croll <tic20>
 * @Date:   05-Feb-2019
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 09-May-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */



#pragma once

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

#include "imex.h"

using namespace clipper;
using namespace clipper::datatypes;

namespace clipper_cx
{

// Structure factor calculation using bulk solvent taking into account
// van der Waals radii
template<class T>
class SFcalc_obs_bulk_vdw : public SFcalc_obs_base<T>
{
public:
    // default constructor
    SFcalc_obs_bulk_vdw(const int n_params=12, const size_t n_threads=8) : nparams(n_params), nthreads(n_threads) {}
    // constructor: shorthand for constructor+operator
    SFcalc_obs_bulk_vdw(HKL_data<F_phi<T> >& fphi, const HKL_data<F_sigF<T> >& fsig,
        const Atom_list& atoms, const int n_params = 12);
    bool operator() ( HKL_data<datatypes::F_phi<T> >& fphi,
            const HKL_data<datatypes::F_sigF<T> >& fsig, const Atom_list& atoms );
    const ftype& bulk_frac() { return bulkfrc; }
    const ftype& bulk_scale() { return bulkscl; }
private:
    int nparams;
    size_t nthreads;
    double bulkfrc, bulkscl;

}; // class SFcalc_obs_bulk_vdw


} // namespace clipper_cx
