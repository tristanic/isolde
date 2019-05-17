/**
 * @Author: Tristan Croll <tic20>
 * @Date:   13-May-2019
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 17-May-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */

//! Bridge between ChimeraX's atomstruct library and Clipper

#include <atomstruct/Atom.h>
#include <clipper/clipper.h>
#include <chrono>
#include <memory>
#include <atomic>

namespace clipper_cx {

namespace bridge {

// template <class ... Types> clipper::Atom cl_atom_from_cx_atom(atomstruct::Atom* cxatom, Types ... args)
// {
//     clipper::String element = cxatom->element().name();
//     const auto& cxcoord = cxatom->coord(args...);
//     auto coord = clipper::Coord_orth(cxcoord[0], cxcoord[1], cxcoord[2]);
//     auto occ = cxatom->occupancy(args...);
//     auto u_iso = clipper::Util::b2u(cxatom->bfactor(args...));
//     clipper::U_aniso_orth uani;
//     if (cxatom->has_aniso_u(args...))
//     {
//         const auto& cxau = *(cxatom->aniso_u(args...));
//         // ChimeraX C++ layer stores aniso_u in row-major order`
//         uani = clipper::U_aniso_orth(cxau[0],cxau[3],cxau[5],cxau[1],cxau[2],cxau[4]);
//     } else {
//         uani = clipper::U_aniso_orth(clipper::U_aniso_orth::null());
//     }
//     auto clatom = clipper::Atom();
//     clatom.set_element(element);
//     clatom.set_coord_orth(coord);
//     clatom.set_occupancy(occ);
//     clatom.set_u_iso(u_iso);
//     clatom.set_u_aniso_orth(uani);
//     return clatom;
// }

template <class ... Types> clipper::Atom cl_atom_from_cx_atom(atomstruct::Atom* cxatom, Types ... args)
{
    clipper::Atom clatom;
    clatom.set_element(cxatom->element().name());
    const auto& cxcoord = cxatom->coord(args...);
    clatom.set_coord_orth(clipper::Coord_orth(cxcoord[0], cxcoord[1], cxcoord[2]));
    clatom.set_occupancy(cxatom->occupancy(args...));
    clatom.set_u_iso(clipper::Util::b2u(cxatom->bfactor(args...)));

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
    clatom.set_u_aniso_orth(uani);
    return clatom;
}


clipper::Atom_list clipper_atoms_from_cx_atoms(atomstruct::Atom** cxatoms, size_t n)
{
    // auto start_time = std::chrono::steady_clock::now();
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
    // auto end_time = std::chrono::steady_clock::now();
    // std::chrono::duration<double> elapsed = end_time-start_time;
    // std::cout<<"Copying " << n << " atoms from ChimeraX to Clipper took " << elapsed.count() << " seconds." << std::endl;
    return al;
}

clipper::Atom_list clipper_atoms_from_cx_atoms_threaded(atomstruct::Atom** cxatoms, size_t n, size_t n_threads)
{
    // auto start_time = std::chrono::steady_clock::now();
    const size_t min_threaded_size = 10000;
    const size_t min_atoms_per_thread = 4000;
    if (n_threads==1 || n < min_threaded_size)
        return clipper_atoms_from_cx_atoms(cxatoms, n);

    size_t atoms_per_thread = std::max(min_atoms_per_thread, n/n_threads+1);

    std::vector<std::future<clipper::Atom_list>> results;
    size_t start=0, end;
    size_t actual_num_threads=0;
    for (size_t i=0; i<n_threads && end<n; ++i)
    {
        end = std::min(n, start+atoms_per_thread);
        results.push_back(std::async(std::launch::async,
            [](atomstruct::Atom** cxatoms, size_t start, size_t end)
            {
                auto al = clipper::Atom_list();
                for (size_t i=start; i<end; ++i)
                {
                    auto cxa = cxatoms[i];
                    auto altlocs = cxa->alt_locs();
                    if (altlocs.size())
                    {
                        for (const auto& altloc: altlocs)
                        {
                            al.push_back(cl_atom_from_cx_atom<char>(cxa, altloc));
                        }
                    } else {
                        al.push_back(cl_atom_from_cx_atom<>(cxa));
                    }
                }
                return al;
            },
            cxatoms, start, end
        ));
        start += atoms_per_thread;
        actual_num_threads++;
    }
    // Since each altloc counts as a Clipper atom, we don't know exactly how
    // many atoms we'll end up with - and tracking the count with a
    // std::atomic<size_t> is empirically counterproductive. It's still
    // worthwhile to reserve a reasonable amount of memory ahead of time -
    // let's just assume that for *most* models fewer than 5% of atoms will have
    // altlocs.
    auto final_al = clipper::Atom_list();
    final_al.reserve(size_t(n*1.05));
    // The final thread will always have the smallest number of atoms in it.
    // Since Clipper doesn't care about atom order, we should start with that
    // first to save time while the other threads complete.
    for (auto it = results.rbegin(); it != results.rend(); it++)
    // for (auto& r: results)
    {
        auto al = it->get();
        // final_al.reserve(final_al.size() + al.size());
        std::move(std::begin(al), std::end(al), std::back_inserter(final_al));
        //final_al.insert(final_al.end(), al.begin(), al.end());
    }
    // auto end_time = std::chrono::steady_clock::now();
    // std::chrono::duration<double> elapsed = end_time-start_time;
    // std::cout<<"Copying " << n << " atoms from ChimeraX to Clipper using " << actual_num_threads << " threads took " << elapsed.count() << " seconds." << std::endl;
    return final_al;
}

} // namespace bridge
} // namespace clipper_cx
