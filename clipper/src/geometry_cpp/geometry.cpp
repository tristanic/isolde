#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "geometry.h"

#include <atomstruct/Atom.h>
#include <atomstruct/Bond.h>
#include <atomstruct/Coord.h>
#include <arrays/pythonarray.h>

#include "../molc.h"

using namespace atomstruct;

void halfbond_cylinder_placement(const Coord &xyz0, const Coord &xyz1, double r, float32_t *m44)
{
    float32_t *m44b = m44 + 16;
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


extern "C" EXPORT PyObject* atom_and_bond_sym_transforms(void *atoms, size_t natoms,
        double *transforms, size_t n_tf, double *center, double cutoff)
{
    Atom **a = static_cast<Atom **>(atoms);
    try {



        Sym_Close_Points cp = Sym_Close_Points();
        std::vector<double> coords(natoms*3);
        Atom **aa = a;
        std::vector<Atom *> visible_atoms;
        // Prune to only visible atoms and get their coordinates
        double *c = coords.data();
        for (size_t i=0; i<natoms; ++i)
        {
            auto this_a = *aa++;
            if (this_a->visible())
            {
                const auto& coord = this_a->coord();
                for (size_t j=0; j<3; ++j)
                    *c++ = coord[j];
                visible_atoms.push_back(this_a);
            }
        }
        natoms = visible_atoms.size();

        // Find the coordinates within the cutoff distance for each transform
        find_close_points_sym(center, cutoff, transforms, n_tf, coords.data(), natoms, cp);

        auto &symmap = cp.symmap;

        std::vector<Atom *> ret_atoms;
        std::vector<double> ret_coords;
        std::vector<uint8_t> ret_atom_sym;
        std::vector<Bond *> ret_bonds;
        std::vector<float32_t> ret_bond_tfs;
        std::vector<uint8_t> ret_bond_sym;
        for (size_t i=0; i<n_tf; ++i)
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
                    if (!bond->shown()) continue;
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
                for (const auto &c: coord)
                    ret_coords.push_back(c);
                ret_atom_sym.push_back(i);
            }
            size_t old_bond_tf_size = ret_bond_tfs.size();
            ret_bond_tfs.resize(old_bond_tf_size + 32*sym_bond_coords.size());
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
        PyObject *ret = PyTuple_New(6);

        // first tuple item: symmetry atoms
        void **aptrs;
        PyObject *atom_ptr_array = python_voidp_array(ret_atoms.size(), &aptrs);
        for (const auto &ptr: ret_atoms)
            *aptrs++ = ptr;
        PyTuple_SET_ITEM(ret, 0, atom_ptr_array);

        // second tuple item: symmetry coords;
        double *acoords;
        PyObject *atom_coord_array = python_double_array(ret_coords.size(), &acoords);
        for (const auto &coord: ret_coords)
            *acoords++ = coord;
        PyTuple_SET_ITEM(ret, 1, atom_coord_array);

        // third tuple item: atom symop indices
        uint8_t *asym;
        PyObject *atom_sym_array = python_uint8_array(ret_atom_sym.size(), &asym);
        for (const auto &sym: ret_atom_sym)
            *(asym++) = sym;
        PyTuple_SET_ITEM(ret, 2, atom_sym_array);

        // fourth tuple item: symmetry bonds
        void **bptrs;
        PyObject *bond_ptr_array = python_voidp_array(ret_bonds.size(), &bptrs);
        for (const auto &ptr: ret_bonds)
            *bptrs++ = ptr;
        PyTuple_SET_ITEM(ret, 3, bond_ptr_array);

        // fifth tuple item: bond transforms
        float *btfs;
        PyObject *bond_transform_array = python_float_array(ret_bond_tfs.size(), &btfs);
        for (const auto &tf: ret_bond_tfs)
            *btfs++ = tf;
        PyTuple_SET_ITEM(ret, 4, bond_transform_array);

        // sixth tuple item: bond symop indices
        uint8_t *bsym;
        PyObject *bond_sym_array = python_uint8_array(ret_bond_sym.size(), &bsym);
        for (const auto &sym: ret_bond_sym)
            *(bsym++) = sym;
        PyTuple_SET_ITEM(ret, 5, bond_sym_array);

        return ret;

    } catch(...) {
        molc_error();
        return 0;
    }
}
