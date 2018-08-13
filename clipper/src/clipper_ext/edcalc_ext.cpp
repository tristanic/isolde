#include "edcalc_ext.h"

using namespace clipper;
using namespace clipper::datatypes;

namespace clipper_cx
{

template <class T>
bool EDcalc_mask_vdw<T>::operator() (Xmap<T>& xmap, const Atom_list& atoms ) const
{
    T zero = 0.0;
    xmap = zero; // Set all values to zero
    const Cell& cell = xmap.cell();
    const Grid_sampling& grid = xmap.grid_sampling();

    Coord_orth xyz;
    Coord_grid g0, g1;
    Grid_range gd(cell, grid, grid_radius_);
    ftype atom_radius;
    typename Xmap<T>::Map_reference_coord i0, iu, iv, iw;
    // Step 1: set all grid points within (atom radius + probe radius) to 1.
    for (const Atom &a: atoms) if ( !a.is_null() )
    {
        try {
            atom_radius = clipper_cx::data::vdw_radii.at(a.element().c_str());
        } catch (...) {
            std::stringstream str;
            str << "Could not find van der Waals radius for atom name " << a.element();
            throw std::runtime_error(str.str());
        }
        xyz = a.coord_orth();
        g0 = xmap.coord_map( xyz ).coord_grid() + gd.min();
        g1 = xmap.coord_map( xyz ).coord_grid() + gd.max();
        i0 = typename Xmap<T>::Map_reference_coord( xmap, g0 );
        for (iu = i0; iu.coord().u() < g1.u(); iu.next_u() ) {
            for (iv=iu; iv.coord().v() < g1.v(); iv.next_v() ) {
                for (iw=iv; iw.coord().w() < g1.w(); iw.next_w() ) {
                    if ( (xyz-iw.coord_orth()).lengthsq() < pow(atom_radius+probe_radius_, 2) )
                        xmap[iw] = 1.0;
                }
            }
        }
    }
    // Step 2: "shrink-wrap" the solvent mask back to the model atoms. For every
    // grid point with a non-zero value, check if it is within shrink_radius_ of
    // a grid point with a value of zero. If so, set its value to zero.
    gd = Grid_range(cell, grid, shrink_radius_);
    Xmap<T> unshrunk(xmap);
    typename Xmap<T>::Map_reference_index ix;
    for (ix = unshrunk.first(); !ix.last(); ix.next()) {
        if (unshrunk[ix] != zero)
        {
            xyz = ix.coord_orth();
            g0 = xmap.coord_map( xyz ).coord_grid() + gd.min();
            g1 = xmap.coord_map( xyz ).coord_grid() + gd.max();
            i0 = typename Xmap<T>::Map_reference_coord( xmap, g0 );
            for (iu = i0; iu.coord().u() < g1.u(); iu.next_u() ) {
                for (iv=iu; iv.coord().v() < g1.v(); iv.next_v() ) {
                    for (iw=iv; iw.coord().w() < g1.w(); iw.next_w() ) {
                        if (unshrunk[iw] == zero) {
                            if ( (xyz-iw.coord_orth()).lengthsq() < pow(shrink_radius_, 2) ) {
                                std::cerr << "Filling in value at " << iw.coord().format() << std::endl;
                                xmap[iw] = 0.0;

                            }
                        }
                    }
                }
            }

        }
    }
    return true;
}


template<class T>
bool EDcalc_mask_vdw<T>::operator() (NXmap<T>& nxmap, const Atom_list& atoms) const
{
    throw std::logic_error("Not implemented");
    return false;
}

template class CLIPPER_CX_IMEX EDcalc_mask_vdw<ftype32>;
template class CLIPPER_CX_IMEX EDcalc_mask_vdw<ftype64>;


} // namespace clipper_cx
