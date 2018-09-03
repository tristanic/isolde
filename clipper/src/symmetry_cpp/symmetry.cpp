/**
 * @Author: Tristan Croll
 * @Date:   04-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   Tristan Croll
 * @Last modified time: 18-Apr-2018
 * @License: Creative Commons BY-NC-SA 3.0, https://creativecommons.org/licenses/by-nc-sa/3.0/.
 * @Copyright: Copyright 2017-2018 Tristan Croll
 */



#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <iostream>
#include <algorithm>

#include "symmetry.h"

#include <atomstruct/Atom.h>
#include <atomstruct/Bond.h>
#include <atomstruct/Residue.h>
#include <atomstruct/Coord.h>
#include <arrays/pythonarray.h>

#include "../molc.h"

using namespace atomstruct;

void halfbond_cylinder_placement(const Coord &xyz0, const Coord &xyz1, double r, float32_t *m44)
{
    float32_t *m44b = m44 + GL_TF_SIZE;
    float x0 = xyz0[0], y0 = xyz0[1], z0 = xyz0[2], x1 = xyz1[0], y1 = xyz1[1], z1 = xyz1[2];
	float vx = x1-x0, vy = y1-y0, vz = z1-z0;
	float d = sqrtf(vx*vx + vy*vy + vz*vz);
	if (d == 0)
	  { vx = vy = 0 ; vz = 1; }
	else
	  { vx /= d; vy /= d; vz /= d; }

	float c = vz, c1;
	if (c <= -1)
	  c1 = 0;       // Degenerate -z axis case.
	else
	  c1 = 1.0/(1+c);

	float wx = -vy, wy = vx;
	float cx = c1*wx, cy = c1*wy;
	float h = d;

	*m44++ = *m44b++ = r*(cx*wx + c);
	*m44++ = *m44b++ = r*cy*wx;
	*m44++ = *m44b++ = -r*wy;
	*m44++ = *m44b++ = 0;

	*m44++ = *m44b++ = r*cx*wy;
	*m44++ = *m44b++ = r*(cy*wy + c);
	*m44++ = *m44b++ = r*wx;
	*m44++ = *m44b++ = 0;

	*m44++ = *m44b++ = h*wy;
	*m44++ = *m44b++ = -h*wx;
	*m44++ = *m44b++ = h*c;
	*m44++ = *m44b++ = 0;

	*m44++ = .75*x0 + .25*x1;
	*m44++ = .75*y0 + .25*y1;
	*m44++ = .75*z0 + .25*z1;
	*m44++ = 1;

	*m44b++ = .25*x0 + .75*x1;
	*m44b++ = .25*y0 + .75*y1;
	*m44b++ = .25*z0 + .75*z1;
	*m44b++ = 1;
}

// Find all symmetry ribbons where at least one tether projects into the
// sphere defined by center and cutoff. The first transform should always
// be the identity operator.
extern "C" EXPORT PyObject* close_sym_ribbon_transforms(double *tether_coords,
    size_t n_tethers, double *transforms, size_t n_tf, double *center, double cutoff)
{
    try {
        std::vector<size_t> found_tfs;
        double tf_coord[3];
        for (size_t i=1; i<n_tf; ++i)
        {
            double *tf = transforms + i*TF_SIZE;
            double *tc = tether_coords;
            // Add the transform to the list if any tether falls withn the search sphere
            for (size_t j=0; j<n_tethers; ++j) {
                transform_coord(tf, tc, tf_coord);
                if (distance_below_cutoff(tf_coord, center, cutoff))
                {
                    found_tfs.push_back(i);
                    break;
                }
                tc += 3;
            }
        }
        uint8_t *rsym;
        PyObject *sym_array = python_uint8_array(found_tfs.size(), &rsym);
        for (const auto &s: found_tfs)
            *(rsym++) = s;
        return sym_array;
    } catch (...) {
        molc_error();
        return 0;
    }
}

bool only_hidden_by_clipper(Atom *a) {
    return (a->display() && !((a->hide()&a->HIDE_ISOLDE)&&(a->hide()&~(a->HIDE_ISOLDE))));
}

bool only_hidden_by_clipper(Bond *b) {
    const Bond::Atoms &a = b->atoms();
    return b->display() && only_hidden_by_clipper(a[0]) && only_hidden_by_clipper(a[1]);
}

void fill_atom_and_bond_sym_tuple(PyObject *ret,
    const std::vector<Atom *>& primary_atoms, const std::vector<Atom *>& atoms,
    const std::vector<double>& coords, const std::vector<uint8_t>& atom_sym_indices,
    const std::vector<Bond *>& bonds, const std::vector<float32_t>& bond_tfs,
    const std::vector<uint8_t>& bond_sym_indices)
{
    // first tuple item: atoms in the primary asu
    void **paptrs;
    PyObject *primary_atom_ptr_array = python_voidp_array(primary_atoms.size(), &paptrs);
    for (const auto &ptr: primary_atoms)
        *paptrs++ = ptr;
    PyTuple_SET_ITEM(ret, 0, primary_atom_ptr_array);
    // second tuple item: symmetry atoms
    void **aptrs;
    PyObject *atom_ptr_array = python_voidp_array(atoms.size(), &aptrs);
    for (const auto &ptr: atoms)
        *aptrs++ = ptr;
    PyTuple_SET_ITEM(ret, 1, atom_ptr_array);

    // third tuple item: symmetry coords;
    double *acoords;
    PyObject *atom_coord_array = python_double_array(coords.size(), &acoords);
    for (const auto &coord: coords)
        *acoords++ = coord;
    PyTuple_SET_ITEM(ret, 2, atom_coord_array);

    // fourth tuple item: atom symop indices
    uint8_t *asym;
    PyObject *atom_sym_array = python_uint8_array(atom_sym_indices.size(), &asym);
    for (const auto &sym: atom_sym_indices)
        *(asym++) = sym;
    PyTuple_SET_ITEM(ret, 3, atom_sym_array);

    // fifth tuple item: symmetry bonds
    void **bptrs;
    PyObject *bond_ptr_array = python_voidp_array(bonds.size(), &bptrs);
    for (const auto &ptr: bonds)
        *bptrs++ = ptr;
    PyTuple_SET_ITEM(ret, 4, bond_ptr_array);

    // sixth tuple item: bond transforms.
    // This is a bit ugly: the default arrangement for halfbond transforms
    // is to do the first set of halfbonds followed by the second
    // (i.e. AAAAAABBBBBB). But up until now it's been much more convenient
    // to intersperse them (ABABABABAB). So now they need to be remixed into
    // the final arrangement.
    float *btfs;
    size_t n = bond_tfs.size();
    PyObject *bond_transform_array = python_float_array(n, &btfs);
    float *btfs2 = btfs + n/2;
    const float *tf = bond_tfs.data();
    for (size_t i=0; i<n/(GL_TF_SIZE*2); ++i)
    {
        const float *tf2 = tf+GL_TF_SIZE;
        for (size_t j=0; j<GL_TF_SIZE; ++j) {
            *btfs++ = *tf++;
            *btfs2++ = *tf2++;
        }
        tf+=GL_TF_SIZE;
    }
    PyTuple_SET_ITEM(ret, 5, bond_transform_array);

    // seventh tuple item: bond symop indices. One index per halfbond.
    uint8_t *bsym;
    PyObject *bond_sym_array = python_uint8_array(bond_sym_indices.size()*2, &bsym);
    uint8_t *bsym2 = bsym+bond_sym_indices.size();
    for (const auto &sym: bond_sym_indices) {
        *(bsym++) = sym;
        *(bsym2++) = sym;
    }
    PyTuple_SET_ITEM(ret, 6, bond_sym_array);
}

//! Find all whole residues with atoms entering a spherical volume, including symmetry
 /* The first transform must always be the identity operator.
  * Returns a 7-tuple containing:
  *     (1) Atoms from the main model matching the search criteria (these will
  *         be drawn as normal by ChimeraX)
  *     (2) Visible atoms from symmetry copies matching the search criteria
  *     (3) The symmetry coordinates for the atoms in (2)
  *     (4) Indices mapping each atom in (2) to its symmetry operator
  *     (5) All visible bonds between atoms in (2)
  *     (6) The halfbond_cylinder_transforms necessary to draw the bonds in
  *         (5)
  *     (7) Indices mapping each bond in (5) to its symmetry operator.
  */
extern "C" EXPORT PyObject* atom_and_bond_sym_transforms_by_residue(void *residues, size_t nres,
    double *transforms, size_t n_tf, double *center, double cutoff, npy_bool visible_only)
{
    Residue **res = static_cast<Residue **>(residues);
    PyObject *ret = PyTuple_New(7);
    try {
        std::vector<Atom *> primary_atoms;
        std::vector<Atom *> sym_atoms;
        std::vector<double> sym_coords;
        std::vector<uint8_t> atom_sym_indices;
        std::vector<Bond *> sym_bonds;
        std::vector<float32_t> bond_tfs;
        std::vector<uint8_t> bond_sym_indices;

        // Get all the atoms in the search sphere from the primary model first.
        // We don't have to apply a symmetry operator or draw anything for these
        // atoms, but we do need to know what they are so we know what to hide/
        // show.
        Residue **r = res;
        for (size_t i=0; i<nres; ++i)
        {
            for (auto a: (*r)->atoms()) {
                if (distance_below_cutoff(a->coord(), center, cutoff))
                {
                    for (auto aa: (*r)->atoms()) {
                        primary_atoms.push_back(aa);
                    }
                    break;
                }
            }
            r++;
        }
        double tf_coord[3];
        for (size_t i=1; i<n_tf; ++i) {
            r = res;
            double *tf = transforms + i*TF_SIZE;
            std::unordered_map<Atom*, std::array<double, 3>> atom_sym_map;

            for (size_t j = 0; j<nres; ++j) {
                // Find all visible atoms from residues with any atom in the search sphere
                for (auto a: (*r)->atoms()) {
                    transform_coord<double, Coord>(tf, a->coord(), tf_coord);
                    if (distance_below_cutoff(tf_coord, center, cutoff))
                    {
                        // collect all atoms in the residue
                        for (auto aa: (*r)->atoms()) {
                            if (!(visible_only && !only_hidden_by_clipper(aa))) {
                                transform_coord<double, Coord>(tf, aa->coord(), atom_sym_map[aa].data());
                            }
                        }
                        break;
                    }
                }
                r++;
            }
            std::unordered_map<Bond*, std::array<float32_t, GL_TF_SIZE*2> > bond_map;
            for (auto &it: atom_sym_map)
            {
                auto a = it.first;
                const auto &bonds = a->bonds();
                for (auto bb = bonds.begin(); bb!=bonds.end(); ++bb)
                {
                    Bond *b = *bb;
                    if (visible_only && !only_hidden_by_clipper(b)) continue;
                    if (bond_map.find(b) != bond_map.end()) continue;
                    const Bond::Atoms &batoms = b->atoms();
                    auto a1 = atom_sym_map.find(batoms[0]);
                    if (a1 == atom_sym_map.end()) continue;
                    auto a2 = atom_sym_map.find(batoms[1]);
                    if (a2 == atom_sym_map.end()) continue;
                    halfbond_cylinder_placement(
                        Coord(a1->second.data()),
                        Coord(a2->second.data()),
                        b->radius(),
                        bond_map[b].data());
                }
            }
            for (const auto &it: atom_sym_map)
            {
                sym_atoms.push_back(it.first);
                const auto &coord = it.second;
                for (size_t j=0; j<3; ++j) {
                    sym_coords.push_back(coord[j]);
                }
                atom_sym_indices.push_back(i);
            }
            for (const auto &it: bond_map)
            {
                sym_bonds.push_back(it.first);
                const auto &tfs = it.second;
                for (size_t j=0; j<GL_TF_SIZE*2; ++j)
                    bond_tfs.push_back(tfs[j]);
                bond_sym_indices.push_back(i);
            }
        }

        fill_atom_and_bond_sym_tuple(ret, primary_atoms, sym_atoms,
            sym_coords, atom_sym_indices, sym_bonds, bond_tfs, bond_sym_indices);

        return ret;

    } catch (...) {
        Py_XDECREF(ret);
        molc_error();
        return 0;
    }
}

//! Find all atoms in a spherical volume, including symmetry
 /* The first transform must always be the identity operator.
  * Returns a 7-tuple containing:
  *     (1) Atoms from the main model matching the search criteria (these will
  *         be drawn as normal by ChimeraX)
  *     (2) Visible atoms from symmetry copies matching the search criteria
  *     (3) The symmetry coordinates for the atoms in (2)
  *     (4) Indices mapping each atom in (2) to its symmetry operator
  *     (5) All visible bonds between atoms in (2)
  *     (6) The halfbond_cylinder_transforms necessary to draw the bonds in
  *         (5)
  *     (7) Indices mapping each bond in (5) to its symmetry operator.
  */

extern "C" EXPORT PyObject*
atom_and_bond_sym_transforms(void *atoms, size_t natoms, double *transforms,
    size_t n_tf, double *center, double cutoff, npy_bool visible_only)
{
    Atom **a = static_cast<Atom **>(atoms);
    PyObject *ret = PyTuple_New(7);
    try {
        Sym_Close_Points cp = Sym_Close_Points();
        std::vector<double> coords(natoms*3);
        std::vector<Atom *> primary_atoms;
        Atom **aa = a;
        std::vector<Atom *> visible_atoms;
        // Prune to only visible atoms and get their coordinates
        double *c = coords.data();
        for (size_t i=0; i<natoms; ++i)
        {
            auto this_a = *aa++;
            if (distance_below_cutoff(this_a->coord(), center, cutoff))
                primary_atoms.push_back(this_a);
            if (!(visible_only && !only_hidden_by_clipper(this_a))) //(this_a->visible())
            {
                const auto& coord = this_a->coord();
                for (size_t j=0; j<3; ++j)
                    *c++ = coord[j];
                visible_atoms.push_back(this_a);
            }
        }
        natoms = visible_atoms.size();
        coords.resize(natoms*3);
        // Find the coordinates within the cutoff distance for each transform
        find_close_points_sym(center, cutoff, transforms, n_tf, coords.data(), natoms, cp);

        auto &symmap = cp.symmap;

        std::vector<Atom *> ret_atoms;
        std::vector<double> ret_coords;
        std::vector<uint8_t> ret_atom_sym;
        std::vector<Bond *> ret_bonds;
        std::vector<float32_t> ret_bond_tfs;
        std::vector<uint8_t> ret_bond_sym;

        for (size_t i=1; i<n_tf; ++i)
        {
            const auto& indices = symmap[i].first;
            auto& coords = symmap[i].second;
            std::unordered_map<Atom*, Sym_Close_Points::Coord> amap;
            for (size_t j=0; j<indices.size(); ++j)
            {
                copy_coord(coords[j].data(), amap[visible_atoms[indices[j]]].data());
            }
            // Find all the visible bonds for this symmetry operator
            std::unordered_map<Bond*, std::pair<Coord, Coord>> sym_bond_coords;
            for (auto &amap_it: amap)
            {
                auto a = amap_it.first;
                const Atom::Bonds &abonds = a->bonds();
                for (auto b = abonds.begin(); b!= abonds.end(); ++b)
                {
                    Bond *bond = *b;
                    if (visible_only && !only_hidden_by_clipper(bond)) continue;
                    const Bond::Atoms &batoms = bond->atoms();
                    auto bi = sym_bond_coords.find(bond);
                    if (bi != sym_bond_coords.end()) continue;
                    auto a1 = amap.find(batoms[0]);
                    if (a1 == amap.end()) continue;
                    auto a2 = amap.find(batoms[1]);
                    if (a2 == amap.end()) continue;

                    sym_bond_coords.emplace(std::make_pair(bond,
                        std::make_pair(Coord(a1->second.data()), Coord(a2->second.data()))));
                }
            }
            for (const auto &it: amap) {
                ret_atoms.push_back(it.first);
                const auto &coord = it.second;
                for (const auto &c: coord) {
                    ret_coords.push_back(c);
                }
                ret_atom_sym.push_back(i);
            }
            size_t old_bond_tf_size = ret_bond_tfs.size();
            ret_bond_tfs.resize(old_bond_tf_size + GL_TF_SIZE*2*sym_bond_coords.size());
            float32_t *bond_tf = ret_bond_tfs.data() + old_bond_tf_size;
            for (const auto &it: sym_bond_coords) {
                auto b = it.first;
                ret_bonds.push_back(b);
                const auto &cpair = it.second;
                halfbond_cylinder_placement(cpair.first, cpair.second, b->radius(), bond_tf);
                bond_tf += 32;
                ret_bond_sym.push_back(i);
            }
        }

        fill_atom_and_bond_sym_tuple(ret, primary_atoms, ret_atoms, ret_coords,
            ret_atom_sym, ret_bonds, ret_bond_tfs, ret_bond_sym);

        return ret;

    } catch(...) {
        molc_error();
        Py_XDECREF(ret);
        return 0;
    }
}

//! Provide the information necessary to draw a predefined set of symmetry atoms.
 /* The first transform must always be the identity operator.
  * Returns a 7-tuple containing:
  *     (1) Atoms corresponding to the identity symop (these will be drawn as
  *         normal by ChimeraX)
  *     (2) Visible atoms from symmetry copies
  *     (3) The symmetry coordinates for the atoms in (2)
  *     (4) Indices mapping each atom in (2) to its symmetry operator
  *     (5) All visible bonds between atoms in (2)
  *     (6) The halfbond_cylinder_transforms necessary to draw the bonds in
  *         (5)
  *     (7) Indices mapping each bond in (5) to its symmetry operator.
  */
extern "C" EXPORT PyObject*
atom_and_bond_sym_transforms_from_sym_atoms(void *atoms, uint8_t *sym_indices,
    size_t natoms, double *transforms, size_t n_tf, npy_bool visible_only)
{
    Atom** a = static_cast<Atom **>(atoms);
    PyObject *ret = PyTuple_New(7);
    try {
        std::vector<Atom *> primary_atoms;
        std::vector<Atom *> ret_atoms;
        std::vector<double> ret_coords(natoms*3);
        std::vector<uint8_t> ret_atom_sym;
        std::vector<Bond *> ret_bonds;
        std::vector<float32_t> ret_bond_tfs;
        std::vector<uint8_t> ret_bond_sym;

        Atom **aa = a;
        uint8_t *si = sym_indices;
        double *coords = ret_coords.data();
        for (size_t i=0; i<natoms; ++i)
        {
            auto atom = *aa++;
            auto sym_index = *(si++);
            if (sym_index == 0)
                primary_atoms.push_back(atom);
            else if (!(visible_only && !only_hidden_by_clipper(atom)))
            {
                ret_atoms.push_back(atom);
                ret_atom_sym.push_back(sym_index);
                transform_coord(transforms+sym_index*TF_SIZE, atom->coord(), coords);
                coords += 3;
            }
        }
        natoms = ret_atoms.size();
        ret_coords.resize(natoms*3);
        // if (natoms ==0) {
        //     throw std::logic_error("No atoms visible!");
        // }
        if (natoms > 0)
        {
            uint8_t current_sym = ret_atom_sym[0];
            coords = ret_coords.data();
            aa = ret_atoms.data();
            uint8_t *asym = ret_atom_sym.data();
            size_t count = 0;
            for (size_t i=0; i<n_tf && count<natoms; ++i)
            {
                current_sym = *asym;
                std::unordered_map<Atom*, double*> amap;
                while (count < natoms && *(asym++) == current_sym) {
                    amap[*aa++] = coords;
                    coords += 3;
                    count++;
                }
                std::unordered_map<Bond*, std::pair<double*, double*>> sym_bond_coords;
                for (auto &amap_it: amap)
                {
                    auto a = amap_it.first;
                    const Atom::Bonds &abonds = a->bonds();
                    for (auto b = abonds.begin(); b != abonds.end(); ++b)
                    {
                        Bond *bond = *b;
                        if (visible_only && !only_hidden_by_clipper(bond)) continue;
                        const Bond::Atoms &batoms = bond->atoms();
                        auto bi = sym_bond_coords.find(bond);
                        if (bi != sym_bond_coords.end()) continue;
                        auto a1 = amap.find(batoms[0]);
                        if (a1 == amap.end()) continue;
                        auto a2 = amap.find(batoms[1]);
                        if (a2 == amap.end()) continue;

                        sym_bond_coords.emplace(std::make_pair(bond,
                            std::make_pair(a1->second, a2->second)));

                    }
                }
                size_t old_bond_tf_size = ret_bond_tfs.size();
                ret_bond_tfs.resize(old_bond_tf_size + GL_TF_SIZE*2*sym_bond_coords.size());
                float32_t *bond_tf = ret_bond_tfs.data() + old_bond_tf_size;
                for (const auto &it: sym_bond_coords)
                {
                    auto b = it.first;
                    ret_bonds.push_back(b);
                    const auto &cpair = it.second;
                    halfbond_cylinder_placement(Coord(cpair.first), Coord(cpair.second),
                        b->radius(), bond_tf);
                    bond_tf += GL_TF_SIZE*2;
                    ret_bond_sym.push_back(current_sym);
                }
                //current_sym = *asym;
            }
        }
        fill_atom_and_bond_sym_tuple(ret, primary_atoms, ret_atoms, ret_coords,
            ret_atom_sym, ret_bonds, ret_bond_tfs, ret_bond_sym);

        return ret;

    } catch (...) {
        molc_error();
        Py_XDECREF(ret);
        return 0;
    }
}
