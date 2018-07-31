#include <iomanip>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>


#include "type_conversions.h"
#include <clipper/clipper.h>

#include "numpy_helper.h"

namespace py=pybind11;
using namespace clipper;

void init_coords(py::module& m)
{

py::class_<Resolution>(m, "Resolution")
    .def(py::init<>())
    .def(py::init<const ftype&>())
    .def("init", &Resolution::init)
    .def("limit", &Resolution::limit)
    .def_property_readonly("invresolsq_limit", &Resolution::invresolsq_limit)
    .def_property_readonly("is_null", &Resolution::is_null)
    .def("__str__", [](const Resolution& self)
    {
        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << self.limit();
        return std::string("Resolution limit: ") + stream.str() + "Å";
    })
    ;

py::class_<HKL_class>(m, "HKL_class")
    .def(py::init<>())
    .def(py::init<const Spacegroup&, const HKL&>())
    .def_property_readonly("epsilon", &HKL_class::epsilon)
    .def_property_readonly("epsilonc", &HKL_class::epsilonc)
    .def_property_readonly("allowed", &HKL_class::allowed)
    .def_property_readonly("centric", &HKL_class::centric)
    .def_property_readonly("sys_abs", &HKL_class::sys_abs)
    ;

py::class_<RTop_orth, RTop<ftype>>(m, "RTop_orth")
    .def(py::init<>())
    .def(py::init<const RTop<>&>())
    .def(py::init<const Mat33<>&>())
    .def(py::init<const Mat33<>&, const Vec3<>&>())
    .def(py::init<const std::vector<Coord_orth>&, const std::vector<Coord_orth>&>())
    .def(py::init<const std::vector<Coord_orth>&, const std::vector<Coord_orth>&, const std::vector<ftype>&>())
    .def(py::init<const Atom_list&, const Atom_list&>())
    .def("rtop_frac", &RTop_orth::rtop_frac)
    .def_property_readonly("inverse", &RTop_orth::inverse)
    .def("axis_coordinate_near", &RTop_orth::axis_coordinate_near)
    .def_property_readonly("screw_translation", &RTop_orth::screw_translation)
    .def_static("identity", &RTop_orth::identity)
    ;

py::class_<HKL, Vec3<int>>(m, "HKL")
    .def(py::init<>())
    .def(py::init<const Vec3<int>&>())
    .def(py::init<const int&, const int&, const int&>())
    // from numpy array
    .def(py::init([](py::array_t<int> hkl)
    {
        check_numpy_array_shape(hkl, {3}, true);
        int* ptr = (int*)hkl.request().ptr;
        return std::unique_ptr<HKL>(new HKL(ptr[0], ptr[1], ptr[2]));
    }))
    .def_property("h",
        [](const HKL& self){ return self.h(); },
        [](HKL& self, const int& h){ self.h() = h; }
    )
    .def_property("k",
        [](const HKL& self){ return self.h(); },
        [](HKL& self, const int& h){ self.h() = h; }
    )
    .def_property("l",
        [](const HKL& self){ return self.h(); },
        [](HKL& self, const int& h){ self.h() = h; }
    )
    // Get/set to/from numpy
    .def_property("hkl",
        [](const HKL& self) { return array_as_numpy_1d<HKL, int>(self, 3); },
        [](HKL& self, py::array_t<int> hkl) { fill_array_from_numpy_1d(self, 3, hkl); }
    )
    .def("invresolsq", &HKL::invresolsq)
    .def("coord_reci_frac", &HKL::coord_reci_frac)
    .def("coord_reci_orth", &HKL::coord_reci_orth)
    .def("transform", (HKL (HKL::*)(const Symop&) const) &HKL::transform)
    .def("transform", (HKL (HKL::*)(const Isymop&) const) &HKL::transform)
    .def("sym_phase_shift", &HKL::sym_phase_shift)
    .def("__str__", &HKL::format)
    .def("__neg__", [](const HKL& self) { return -self; })
    .def("__add__", [](const HKL& self, const HKL& other){ return self+other; }, py::is_operator())
    .def("__sub__", [](const HKL& self, const HKL& other) { return self-other; }, py::is_operator())
    .def("__mul__", [](const HKL& self, const int& s) { return s*self; }, py::is_operator())
    .def("__rmul__", [](const HKL& self, const int& s) { return s*self; }, py::is_operator())
    .def("__rmul__", [](const HKL& self, const Isymop& op) { return op*self; }, py::is_operator())
    ;

py::class_<Coord_reci_orth, Vec3<>>(m, "Coord_reci_orth")
    .def(py::init<>())
    .def(py::init<const Vec3<>&>())
    .def(py::init<const ftype&, const ftype&, const ftype&>())
    // from numpy array
    .def(py::init([](py::array_t<ftype> xyz_reci)
    {
        check_numpy_array_shape(xyz_reci, {3}, true);
        ftype* ptr = (ftype*)xyz_reci.request().ptr;
        return std::unique_ptr<Coord_reci_orth>(new Coord_reci_orth(ptr[0], ptr[1], ptr[2]));
    }))
    .def_property_readonly("xs", &Coord_reci_orth::xs)
    .def_property_readonly("ys", &Coord_reci_orth::ys)
    .def_property_readonly("zs", &Coord_reci_orth::zs)
    .def_property_readonly("xyz_reci",
        [](const Coord_reci_orth& self) { return array_as_numpy_1d<Coord_reci_orth, ftype>(self, 3); }
    )
    .def_property_readonly("invresolsq", &Coord_reci_orth::invresolsq)
    .def("coord_reci_frac", &Coord_reci_orth::coord_reci_frac)
    .def("transform", &Coord_reci_orth::transform)
    .def("__str__", &Coord_reci_orth::format)
    ;

py::class_<Coord_reci_frac, Vec3<>>(m, "Coord_reci_frac")
    .def(py::init<>())
    .def(py::init<const Vec3<>&>())
    .def(py::init<const ftype&, const ftype&, const ftype&>())
    .def(py::init<const HKL&>())
    // from numpy array
    .def(py::init([](py::array_t<ftype> uvw_reci)
    {
        check_numpy_array_shape(uvw_reci, {3}, true);
        ftype* ptr = (ftype*)uvw_reci.request().ptr;
        return std::unique_ptr<Coord_reci_frac>(new Coord_reci_frac(ptr[0], ptr[1], ptr[2]));
    }))
    .def_property_readonly("hkl", &Coord_reci_frac::hkl)
    .def("invresolsq", &Coord_reci_frac::invresolsq)
    .def_property_readonly("us", &Coord_reci_frac::us)
    .def_property_readonly("vs", &Coord_reci_frac::vs)
    .def_property_readonly("ws", &Coord_reci_frac::ws)
    .def_property_readonly("uvw_reci",
        [](const Coord_reci_frac& self) { return array_as_numpy_1d<Coord_reci_frac, ftype>(self, 3); }
    )
    .def("coord_reci_orth", &Coord_reci_frac::coord_reci_orth)
    .def("transform", &Coord_reci_frac::transform)
    .def("__str__", &Coord_reci_frac::format)
    ;

py::class_<Coord_grid, Vec3<int>>(m, "Coord_grid")
    .def(py::init<>())
    .def(py::init<const Vec3<int>>())
    .def(py::init<const int&, const int&, const int&>())
    .def(py::init<const Grid&, const int&>())
    // from numpy array
    .def(py::init([](py::array_t<int> uvw)
    {
        check_numpy_array_shape(uvw, {3}, true);
        int* ptr = (int*)uvw.request().ptr;
        return std::unique_ptr<Coord_grid>(new Coord_grid(ptr[0], ptr[1], ptr[2]));
    }))
    .def_property("u",
        [](const Coord_grid& self){ return self.u(); },
        [](Coord_grid& self, const int& u){ self.u() = u; }
    )
    .def_property("v",
        [](const Coord_grid& self){ return self.v(); },
        [](Coord_grid& self, const int& v){ self.v() = v; }
    )
    .def_property("w",
        [](const Coord_grid& self){ return self.w(); },
        [](Coord_grid& self, const int& w){ self.w() = w; }
    )
    .def_property("uvw",
        [](const Coord_grid& self) { return array_as_numpy_1d<Coord_grid, int>(self, 3); },
        [](Coord_grid& self, py::array_t<int> vals) { fill_array_from_numpy_1d<Coord_grid, int>(self, 3, vals); }
    )
    .def("coord_map", &Coord_grid::coord_map)
    .def("coord_frac", &Coord_grid::coord_frac)
    .def("transform", &Coord_grid::transform)
    .def("unit", &Coord_grid::unit)
    .def("next", (const Coord_grid& (Coord_grid::*)(const Grid&)) &Coord_grid::next)
    .def("next", (const Coord_grid& (Coord_grid::*)(const Grid_range&)) &Coord_grid::next)
    .def("last", (bool (Coord_grid::*)(const Grid&) const) &Coord_grid::last)
    .def("last", (bool (Coord_grid::*)(const Grid_range&) const) &Coord_grid::last)
    .def("index", &Coord_grid::index)
    .def("deindex", &Coord_grid::deindex)
    .def("__str__", &Coord_grid::format)
    .def("__neg__", [](const Coord_grid& self) { return -self; })
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(int() * py::self)
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(Isymop() * py::self)
    ;

py::class_<Coord_orth, Vec3<>>(m, "Coord_orth")
    .def(py::init<>())
    .def(py::init<const Vec3<>&>())
    .def(py::init<const ftype&, const ftype&, const ftype&>())
    .def(py::init<const Coord_orth&, const Coord_orth&, const Coord_orth&, const ftype&, const ftype&, const ftype&>())
    // from numpy array
    .def(py::init([](py::array_t<ftype> xyz)
    {
        check_numpy_array_shape(xyz, {3}, true);
        ftype* ptr = (ftype*)xyz.request().ptr;
        return std::unique_ptr<Coord_orth>(new Coord_orth(ptr[0], ptr[1], ptr[2]));
    }))
    .def_property_readonly("x", &Coord_orth::x)
    .def_property_readonly("y", &Coord_orth::y)
    .def_property_readonly("z", &Coord_orth::z)
    .def_property_readonly("xyz",
        [](const Coord_orth& self) { return array_as_numpy_1d<Coord_orth, ftype>(self, 3); }
    )
    .def_property_readonly("lengthsq", &Coord_orth::lengthsq)
    .def("coord_frac", &Coord_orth::coord_frac)
    .def("transform", &Coord_orth::transform)
    .def("__str__", &Coord_orth::format)
    .def_static("length", &Coord_orth::length)
    .def_static("angle", &Coord_orth::angle)
    .def_static("torsion", &Coord_orth::torsion)
    .def("__neg__", [](const Coord_orth& self){ return -self; })
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(ftype() * py::self)
    .def(RTop_orth() * py::self)
    ;

py::class_<Coord_frac, Vec3<>>(m, "Coord_frac")
    .def(py::init<>())
    .def(py::init<const Vec3<>&>())
    .def(py::init<const ftype&, const ftype&, const ftype&>())
    // // from numpy array
    .def(py::init([](py::array_t<ftype> uvw)
    {
        check_numpy_array_shape(uvw, {3}, true);
        ftype* ptr = (ftype*)uvw.request().ptr;
        return std::unique_ptr<Coord_frac>(new Coord_frac(ptr[0], ptr[1], ptr[2]));
    }))
    .def_property_readonly("u", &Coord_frac::u)
    .def_property_readonly("v", &Coord_frac::v)
    .def_property_readonly("w", &Coord_frac::w)
    .def_property_readonly("uvw",
        [](const Coord_frac& self) { return array_as_numpy_1d<Coord_frac, ftype>(self, 3); }
    )
    .def("lengthsq", &Coord_frac::lengthsq)
    .def("coord_orth", &Coord_frac::coord_orth)
    .def("coord_map", &Coord_frac::coord_map)
    .def("coord_grid", &Coord_frac::coord_grid)
    .def("transform", &Coord_frac::transform)
    //.def("__rmul__", [](const Coord_frac& self, const RTop_frac& op) { return op*self; })
    .def("lattice_copy_zero", &Coord_frac::lattice_copy_zero)
    .def("lattice_copy_unit", &Coord_frac::lattice_copy_unit)
    .def("lattice_copy_near", &Coord_frac::lattice_copy_near)
    .def("symmetry_copy_near", &Coord_frac::symmetry_copy_near)
    .def("__str__", &Coord_frac::format)
    .def("__neg__", [](const Coord_frac& self) { return -self; })
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(ftype() * py::self)
    .def("__mul__", [](const Coord_frac& self, const ftype& s) { return s * self; })
    .def(RTop_frac() * py::self)
    ;

py::class_<Coord_map, Vec3<>>(m, "Coord_map")
    .def(py::init<>())
    .def(py::init<const Vec3<>&>())
    .def(py::init<const Coord_grid&>())
    .def(py::init<const ftype&, const ftype&, const ftype&>())
    // from numpy array
    .def(py::init([](py::array_t<ftype> uvw)
    {
        check_numpy_array_shape(uvw, {3}, true);
        ftype* ptr = (ftype*)uvw.request().ptr;
        return std::unique_ptr<Coord_map>(new Coord_map(ptr[0], ptr[1], ptr[2]));
    }))
    .def("coord_frac", &Coord_map::coord_frac)
    .def("coord_grid", &Coord_map::coord_grid)
    .def_property_readonly("floor", &Coord_map::floor)
    .def_property_readonly("ceil", &Coord_map::ceil)
    .def_property_readonly("u", &Coord_map::u)
    .def_property_readonly("v", &Coord_map::v)
    .def_property_readonly("w", &Coord_map::w)
    .def_property_readonly("uvw",
        [](const Coord_map& self) { return array_as_numpy_1d<Coord_map, ftype>(self, 3); }
    )
    .def("__str__", &Coord_map::format)
    .def("__neg__", [](const Coord_map& self) { return -self; })
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(ftype() * py::self)
    .def("__rmul__", [](const Coord_map& self, const ftype& s) { return s*self; })
    ;

py::class_<U_aniso_orth, Mat33sym<>>(m, "U_aniso_orth")
    .def(py::init<>())
    .def(py::init<const Mat33sym<>&>())
    .def(py::init<const ftype&>())
    .def(py::init<const ftype&, const ftype&, const ftype&, const ftype&, const ftype&, const ftype&>())
    // from numpy
    .def(py::init([](py::array_t<ftype> vals)
    {
        ftype* ptr = (ftype*)vals.request().ptr;
        return std::unique_ptr<U_aniso_orth>(new U_aniso_orth(ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5]));
    }))
    .def_property_readonly("u_iso", &U_aniso_orth::u_iso)
    .def("u_aniso_frac", &U_aniso_orth::u_aniso_frac)
    .def("transform", &U_aniso_orth::transform)
    .def("__rmul__", [](const U_aniso_orth& self, const RTop_orth& rt) { return self.transform(rt); })
    .def(py::self + py::self)
    .def("__neg__", [](const U_aniso_orth& self) { return -self; })
    .def(ftype() * py::self)
    .def("__rmul__", [](const U_aniso_orth& self, const ftype& s) { return s * self; })
    ;

py::class_<U_aniso_frac, Mat33sym<>>(m, "U_aniso_frac")
    .def(py::init<>())
    .def(py::init<const Mat33sym<>&>())
    .def(py::init<const ftype&, const ftype&, const ftype&, const ftype&, const ftype&, const ftype&>())
    // from numpy
    .def(py::init([](py::array_t<ftype> vals)
    {
        check_numpy_array_shape(vals, {6}, true);
        ftype* ptr = (ftype*)vals.request().ptr;
        return std::unique_ptr<U_aniso_frac>(new U_aniso_frac(ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5]));
    }))
    .def("u_aniso_orth", &U_aniso_frac::u_aniso_orth)
    .def("transform", &U_aniso_frac::transform)
    .def("__rmul__", [](const U_aniso_frac& self, const RTop_frac& rt) { return self.transform(rt); })
    .def(py::self+py::self)
    .def("__neg__", [](const U_aniso_frac& self) { return -self; })
    .def(ftype() * py::self)
    .def("__rmul__", [](const U_aniso_frac& self, const ftype& s) { return s*self; })
    ;

py::class_<Grid, Vec3<int>>(m, "Grid")
    .def(py::init<>())
    .def(py::init<const int&, const int&, const int&>())
    // from numpy
    .def(py::init([](py::array_t<int> vals)
    {
        auto ret = new Grid();
        fill_array_from_numpy_1d<Grid, int>(*ret, 3, vals);
        return std::unique_ptr<Grid>(ret);
    }))
    .def_property_readonly("nu", &Grid::nu)
    .def_property_readonly("nv", &Grid::nv)
    .def_property_readonly("nw", &Grid::nw)
    .def_property_readonly("nuvw",
        [](const Grid& self){ return array_as_numpy_1d<Grid, int>(self, 3); }
    )
    .def_property_readonly("size", &Grid::size)
    .def("in_grid", &Grid::in_grid)
    .def("index", &Grid::index)
    .def("deindex", &Grid::deindex)
    .def("__str__", &Grid::format)
    // Extra useful features for Python
    .def_property_readonly("dim", [](const Grid& self) -> py::array_t<int>
    {
        py::array_t<int> ret(3);
        int* ptr = (int*)ret.request().ptr;
        *ptr++ = self.nu();
        *ptr++ = self.nv();
        *ptr++ = self.nw();
        return ret;
    })
;

py::class_<Grid_sampling, Grid>(m, "Grid_sampling")
    .def(py::init<>())
    .def(py::init<const int&, const int&, const int&>())
    .def(py::init<const Spacegroup&, const Cell&, const Resolution&, const ftype>(),
        py::arg("spacegroup"), py::arg("cell"), py::arg("resolution"), py::arg("rate")=1.5)
    // from numpy
    .def(py::init([](py::array_t<int> vals)
    {
        auto ret = new Grid_sampling();
        fill_array_from_numpy_1d<Grid_sampling, int>(*ret, 3, vals);
        return std::unique_ptr<Grid_sampling>(ret);
    }))
    .def("init", &Grid_sampling::init)
    .def_property_readonly("matrix_grid_frac", &Grid_sampling::matrix_grid_frac)
    .def_property_readonly("matrix_frac_grid", &Grid_sampling::matrix_frac_grid)
    .def_property_readonly("is_null", &Grid_sampling::is_null)
    ;

py::class_<HKL_sampling>(m, "HKL_sampling")
    .def(py::init<>())
    .def(py::init<const Cell&, const Resolution&>())
    .def_property_readonly("hkl_limit", &HKL_sampling::hkl_limit)
    .def("resolution", &HKL_sampling::resolution)
    .def("in_resolution", &HKL_sampling::in_resolution)
    .def_property_readonly("is_null", &HKL_sampling::is_null)
    .def("__str__", &HKL_sampling::format)
    .def(py::self == py::self)
    ;

py::class_<Grid_range, Grid>(m, "Grid_range")
    .def(py::init<>())
    .def(py::init<const Coord_grid&, const Coord_grid&>())
    .def(py::init<const Grid&, const Coord_frac&, const Coord_frac&>())
    .def(py::init<const Cell&, const Grid&, const ftype&>())
    .def_property_readonly("min", &Grid_range::min)
    .def_property_readonly("max", &Grid_range::max)
    .def("add_border", &Grid_range::add_border)
    .def("in_grid", &Grid_range::in_grid)
    .def("index", &Grid_range::index)
    .def("deindex", &Grid_range::deindex)
    ;

py::class_<Atom>(m, "Atom")
    .def(py::init<>())
    .def(py::init<const Atom&>())
    .def_property("element",
        [](const Atom& self) { return self.element(); },
        [](Atom& self, const String& e) { self.set_element(e); }
    )
    .def_property("coord_orth",
        //[](const Atom& self) { return array_as_numpy_1d<Coord_orth, ftype>(self.coord_orth()); }
        &Atom::coord_orth,
        [](Atom& self, py::array_t<ftype> coord)
        {
            auto c = Coord_orth();
            fill_array_from_numpy_1d<Coord_orth, ftype>(c, 3, coord);
            self.set_coord_orth(c);
        }
    )
    .def_property("occupancy", &Atom::occupancy, &Atom::set_occupancy)
    .def_property("u_iso", &Atom::u_iso, &Atom::set_u_iso)
    .def_property("u_aniso_orth",
        &Atom::u_aniso_orth,
        [](Atom& self, py::array_t<ftype> u)
        {
            std::vector<ftype> vals(6);
            fill_array_from_numpy_1d<std::vector<ftype>, ftype>(vals, 6, u);
            auto ua = U_aniso_orth(vals[0], vals[1], vals[2], vals[3], vals[4], vals[5]);
            self.set_u_aniso_orth(ua);
        }
    )
    .def("set_element", &Atom::set_element)
    .def("set_element", [](Atom& self, const std::string& ename) { self.set_element(String(ename)); } )
    .def("set_coord_orth", &Atom::set_coord_orth)
    .def("set_coord_orth", [](Atom& self, py::array_t<ftype> coord)
    {
        auto c = Coord_orth();
        fill_array_from_numpy_1d<Coord_orth, ftype>(c, 3, coord);
        self.set_coord_orth(c);
    })
    .def("set_occupancy", &Atom::set_occupancy)
    .def("set_u_iso", &Atom::set_u_iso)
    .def("set_u_aniso_orth", &Atom::set_u_aniso_orth)
    .def("set_u_aniso_orth", [](Atom& self, py::array_t<ftype> u)
    {
        std::vector<ftype> vals(6);
        fill_array_from_numpy_1d<std::vector<ftype>, ftype>(vals, 6, u);
        auto ua = U_aniso_orth(vals[0], vals[1], vals[2], vals[3], vals[4], vals[5]);
        self.set_u_aniso_orth(ua);
    })
    .def("transform", &Atom::transform)
    .def("__rmul__", [](Atom& self, const RTop_orth& rt) { self.transform(rt); })
    .def_property_readonly("is_null", &Atom::is_null)
    ;


const char* atom_list_docstring_=
"Manages a list of :class:`Atom` objects, with functions to add, remove or modify \
atoms individually or in groups. One big caveat applies: \n\
\n\
    - Atoms are retrieved and inserted as copies, not references. If you \
      retrieve an atom from the list and modify it, you will need to re-insert \
      it for the changes to have any effect. \
";

py::class_<Atom_list>(m, "Atom_list", atom_list_docstring_)
    .def(py::init<>())
    .def(py::init<const std::vector<Atom>&>())
    // from arrays of atom properties
    .def(py::init([](const std::vector<std::string>& elements, py::array_t<ftype> coords,
                     py::array_t<ftype> occupancies, py::array_t<ftype> u_isos, py::array_t<ftype> u_anisos)
                     -> std::unique_ptr<Atom_list>
    {
        int n = elements.size();
        check_numpy_array_shape(coords, {n,3}, true);
        check_numpy_array_shape(u_anisos, {n, 6}, true);
        check_numpy_array_shape(occupancies, {n}, true);
        check_numpy_array_shape(u_isos, {n}, true);
        String e;
        Coord_orth c;
        U_aniso_orth ua;
        ftype* c_ptr = (ftype*)coords.request().ptr;
        ftype* ua_ptr = (ftype*)u_anisos.request().ptr;
        ftype* o_ptr = (ftype*)occupancies.request().ptr;
        ftype* ui_ptr = (ftype*)u_isos.request().ptr;
        Atom_list * al = new Atom_list();
        for (int i=0; i<n; ++i) {
            auto a = Atom();
            a.set_element(String(elements[i]));
            a.set_coord_orth(Coord_orth(*c_ptr, *(c_ptr+1), *(c_ptr+2)));
            c_ptr+=3;
            a.set_u_aniso_orth(U_aniso_orth(*ua_ptr, *(ua_ptr+1), *(ua_ptr+2), *(ua_ptr+3), *(ua_ptr+4), *(ua_ptr+5)));
            ua_ptr+=6;
            a.set_occupancy(*o_ptr++);
            a.set_u_iso(*ui_ptr++);
            al->push_back(a);
        }
        return std::unique_ptr<Atom_list>(al);
    }))
    .def("delete_atom", [](Atom_list& self, const int& i)
    {
        self.erase(self.begin()+i);
    }, "Deletes the atom at the given index and adjusts the array indexing accordingly.")
    .def("__getitem__", [](const Atom_list& self, const int& i) -> Atom
        { return self.at(i); },
        "Returns a copy of the designated atom.")
    .def("__setitem__", [](Atom_list& self, const int&i, const Atom& atom)
        {
            self.at(i) = atom;
        },
        "Replaces the atom at the given index"
    )
    .def("get_elements", [](const Atom_list& self, py::array_t<int> indices) -> std::vector<std::string>
    {
        std::vector<std::string> ret;
        auto buf = indices.request();
        if (buf.ndim !=1)
            throw std::runtime_error("Indices must be a 1D array!");
        int n = buf.shape[0];
        int* iptr = (int*)buf.ptr;
        for (int i=0; i<n; ++i)
            ret.push_back(self.at(*(iptr++)).element());
        return ret;
    },
    "Get the names of the elements at the given indices."
    )
    .def("set_elements", [](Atom_list& self, py::array_t<int> indices, const std::vector<std::string>& elements)
    {
        try {
            check_numpy_array_shape(indices, {(int)elements.size()}, true);
        } catch (std::runtime_error) {
            throw std::runtime_error("Number of elements doesn't match the number of indices!");
        }
        int* ptr = (int*)indices.request().ptr;
        for (size_t i=0; i<elements.size(); ++i)
            self.at(*ptr++).set_element(String(elements[i]));
    },
    "Change the elements at the given indices."
    )
    .def("get_coords", [](const Atom_list& self, py::array_t<int> indices) -> py::array_t<ftype>
    {
        auto ibuf = indices.request();
        if (ibuf.ndim !=1)
            throw std::runtime_error("Indices must be a 1D array!");
        int n = ibuf.shape[0];
        int* iptr = (int*)ibuf.ptr;
        auto ret = py::array_t<ftype>({n,3});
        ftype* rptr = (ftype*)ret.request().ptr;
        for (int i=0; i<n; ++i)
        {
            const auto& coord = self.at(*(iptr++)).coord_orth();
            for (int j=0; j<3; ++j)
                *rptr++ = coord[j];
        }
        return ret;
    },
    "Get the coordinates of the atoms at the given indices in Angstroms."
    )
    .def("set_coords", [](Atom_list& self, py::array_t<int> indices, py::array_t<ftype> coords)
    {
        auto ibuf = indices.request();
        auto cbuf = coords.request();
        if (ibuf.ndim !=1 || cbuf.ndim !=2 || cbuf.shape[0] != ibuf.shape[0] ||
            cbuf.shape[1] != 3)
            throw std::runtime_error("Malformed input arrays! indices must be a 1D "
                "array of length n, and coords must be a n x 3 array.");
        int* iptr = (int*)ibuf.ptr;
        ftype* cptr = (ftype*)cbuf.ptr;
        int n = ibuf.shape[0];
        for (int i=0; i<n; ++i) {
            self.at(*(iptr++)).set_coord_orth(Coord_orth(*cptr, *(cptr+1), *(cptr+2)));
            cptr += 3;
        }
    },
    "Set the coordinates for the given atoms in Angstroms."
    )
    .def("get_u_anisos", [](const Atom_list& self, py::array_t<int> indices) -> py::array_t<ftype>
    {
        auto ibuf = indices.request();
        if (ibuf.ndim !=1)
            throw std::runtime_error("Indices must be a 1D array!");
        int n = ibuf.shape[0];
        int* iptr = (int*)ibuf.ptr;
        auto ret = py::array_t<ftype>({n,6});
        ftype* rptr = (ftype*)ret.request().ptr;
        for (int i=0; i<n; ++i)
        {
            const auto& ua = self.at(*(iptr++)).u_aniso_orth();
            *rptr++ = ua.mat00();
            *rptr++ = ua.mat11();
            *rptr++ = ua.mat22();
            *rptr++ = ua.mat01();
            *rptr++ = ua.mat02();
            *rptr++ = ua.mat12();
        }
        return ret;
    },
    "Get the anisotropic B-factor parameters for the atoms at the given indices."
    )
    .def("set_u_anisos", [](Atom_list& self, py::array_t<int> indices, py::array_t<ftype> uas)
    {
        auto ibuf = indices.request();
        auto ubuf = uas.request();
        if (ibuf.ndim !=1 || ubuf.ndim !=2 || ubuf.shape[0] != ibuf.shape[0] ||
            ubuf.shape[1] != 6)
            throw std::runtime_error("Malformed input arrays! indices must be a 1D "
                "array of length n, and u_anisos must be a n x 6 array.");
        int* iptr = (int*)ibuf.ptr;
        ftype* uptr = (ftype*)ubuf.ptr;
        int n = ibuf.shape[0];
        for(int i=0; i<n; ++i) {
            self.at(*(iptr++)).set_u_aniso_orth(
                U_aniso_orth(*uptr, *(uptr+1), *(uptr+2), *(uptr+3), *(uptr+4), *(uptr+5))
            );
            uptr+=6;
        }
    },
    "Set the anisotropic B-factors for the given atoms."
    )
    .def("get_occupancies", [](const Atom_list& self, py::array_t<int> indices) -> py::array_t<ftype>
    {
        auto ibuf = indices.request();
        if (ibuf.ndim !=1)
            throw std::runtime_error("Indices must be a 1D array!");
        int n = ibuf.shape[0];
        int* iptr = (int*)ibuf.ptr;
        auto ret = py::array_t<ftype>(n);
        ftype* rptr = (ftype*)ret.request().ptr;
        for (int i=0; i<n; ++i)
            *rptr++ = self.at(*(iptr++)).occupancy();
        return ret;
    },
    "Get the fractional occupancies of the atoms at the given indices."
    )
    .def("set_occupancies", [](Atom_list& self, py::array_t<int> indices, py::array_t<ftype> occs)
    {
        auto ibuf = indices.request();
        auto obuf = occs.request();
        if (ibuf.ndim !=1 || obuf.ndim !=1 || obuf.shape[0] != ibuf.shape[0])
            throw std::runtime_error("Index and occupancy arrays must be the same length!");
        int* iptr = (int*)ibuf.ptr;
        ftype* optr = (ftype*)obuf.ptr;
        int n = ibuf.shape[0];
        for (int i=0; i<n; ++i) {
            self.at(*(iptr++)).set_occupancy(*(optr++));
        }
    },
    "Set the fractional occupancies for the given atoms."
    )
    .def("get_u_isos", [](const Atom_list& self, py::array_t<int> indices) -> py::array_t<ftype>
    {
        auto ibuf = indices.request();
        if (ibuf.ndim !=1)
            throw std::runtime_error("Indices must be a 1D array!");
        int n = ibuf.shape[0];
        int* iptr = (int*)ibuf.ptr;
        auto ret = py::array_t<ftype>(n);
        ftype* rptr = (ftype*)ret.request().ptr;
        for (int i=0; i<n; ++i)
            *rptr++ = self.at(*(iptr++)).u_iso();
        return ret;
    },
    "Get the isotropic B-factors of the atoms at the given indices."
    )
    .def("set_u_isos", [](Atom_list& self, py::array_t<int> indices, py::array_t<ftype> u_isos)
    {
        auto ibuf = indices.request();
        auto ubuf = u_isos.request();
        if (ibuf.ndim !=1 || ubuf.ndim !=1 || ubuf.shape[0] != ibuf.shape[0])
            throw std::runtime_error("Index and u_iso arrays must be the same length!");
        int* iptr = (int*)ibuf.ptr;
        ftype* uptr = (ftype*)ubuf.ptr;
        int n = ibuf.shape[0];
        for (int i=0; i<n; ++i) {
            self.at(*(iptr++)).set_u_iso(*(uptr++));
        }
    },
    "Set the occupancies for the given atoms."
    )
    ;







} //init_coords
