/**
 * @Author: Tristan Croll <tic20>
 * @Date:   13-May-2019
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 13-May-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */

//! Bridge between ChimeraX's atomstruct library and Clipper

#include <atomstruct/Atom.h>
#include <clipper/clipper.h>

namespace clipper_cx {

namespace bridge {

template <class ... Types> clipper::Atom cl_atom_from_cx_atom(atomstruct::Atom* cxatom, Types ... args)
{
    clipper::String element = cxatom->element().name();
    const auto& cxcoord = cxatom->coord(args...);
    auto coord = clipper::Coord_orth(cxcoord[0], cxcoord[1], cxcoord[2]);
    auto occ = cxatom->occupancy(args...);
    auto u_iso = clipper::Util::b2u(cxatom->bfactor(args...));
    clipper::U_aniso_orth uani;
    if (cxatom->has_aniso_u(args...))
    {
        const auto& cxau = *(cxatom->aniso_u(args...));
        // ChimeraX C++ layer stores aniso_u in row-major order`
        uani = clipper::U_aniso_orth(cxau[0],cxau[3],cxau[5],cxau[1],cxau[2],cxau[4]);
    } else {
        uani = clipper::U_aniso_orth(clipper::U_aniso_orth::null());
    }
    auto clatom = clipper::Atom();
    clatom.set_element(element);
    clatom.set_coord_orth(coord);
    clatom.set_occupancy(occ);
    clatom.set_u_iso(u_iso);
    clatom.set_u_aniso_orth(uani);
    return clatom;
}

clipper::Atom_list clipper_atoms_from_cx_atoms(atomstruct::Atom** cxatoms, size_t n)
{
    auto al = clipper::Atom_list();
    for (size_t i=0; i<n; ++i)
    {
        auto cxa = cxatoms[i];
        auto altlocs = cxa->alt_locs();
        if (altlocs.size())
        {
            for (const auto& altloc: altlocs)
                al.push_back(cl_atom_from_cx_atom<char>(cxa, altloc));
        } else {
            al.push_back(cl_atom_from_cx_atom<>(cxa));
        }
    }
    return al;
}

} // namespace bridge
} // namespace clipper_cx
