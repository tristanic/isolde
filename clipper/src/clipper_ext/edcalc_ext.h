/**
 * @Author: Tristan Croll <tic20>
 * @Date:   13-Aug-2018
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
#include "vdw.h"

using namespace clipper;
using namespace clipper::datatypes;

namespace clipper_cx
{

//! Bulk solvent mask taking into account van der Waals radii
/*! As per Jiang & Brunger (1994), J Mol Biol 243: 100-115 */
template <class T>
class EDcalc_mask_vdw: public EDcalc_base<T>
{
public:
    // In my empirical tests, these slightly smaller probe/shrink radii seem better than the published ones
    EDcalc_mask_vdw( const ftype grid_radius = 3.0, const ftype probe_radius = 0.6, const ftype shrink_radius = 0.7)
        : grid_radius_(grid_radius), probe_radius_(probe_radius), shrink_radius_(shrink_radius)
    {}
    bool operator() (Xmap<T>& xmap, const Atom_list& atoms) const;
    bool operator() (NXmap<T>& nxmap, const Atom_list& atoms) const;

private:
    const ftype grid_radius_;
    const ftype probe_radius_;
    const ftype shrink_radius_;
}; // class EDcalc_mask_vdw

template<class T> class EDcalc_aniso_thread : public EDcalc_base<T> {
public:
    EDcalc_aniso_thread( const ftype radius = 2.5, const size_t n_threads = 2 ) : radius_(radius), n_threads_(n_threads) {}
    bool operator() (Xmap<T>& xmap, const Atom_list& atoms) const;
    bool operator() ( NXmap<T>& nxmap, const Atom_list& atoms) const {return false;}
private:
    const ftype radius_;
    const size_t n_threads_;
    bool edcalc_xmap_thread_(Xmap<T>& xmap, const Atom_list& atoms, size_t start, size_t end) const;
    bool edcalc_nxmap_thread_(NXmap<T>& nxmap, const Atom_list& atoms, size_t start, size_t end) const {return false;}
};

} // namespace clipper_cx
