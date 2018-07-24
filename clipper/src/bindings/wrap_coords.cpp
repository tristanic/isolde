#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>


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
    .def("invresolsq_limit", &Resolution::invresolsq_limit)
    .def("is_null", &Resolution::is_null)
    .def("__str__", [](const Resolution& self)
    {
        return std::string(String(self.limit())) + " Å";
    })
    ;

py::class_<HKL_class>(m, "HKL_class")
    .def(py::init<>())
    .def(py::init<const Spacegroup&, const HKL&>())
    .def("epsilon", &HKL_class::epsilon)
    .def("epsilonc", &HKL_class::epsilonc)
    .def("allowed", &HKL_class::allowed)
    .def("centric", &HKL_class::centric)
    .def("sys_abs", &HKL_class::sys_abs)
    ;

py::class_<RTop_orth>(m, "RTop_orth")
    .def(py::init<>())
    .def(py::init<const RTop<>&>())
    .def(py::init<const Mat33<>&>())
    .def(py::init<const Mat33<>&, const Vec3<>&>())
    .def(py::init<const std::vector<Coord_orth>&, const std::vector<Coord_orth>&>())
    .def(py::init<const std::vector<Coord_orth>&, const std::vector<Coord_orth>&, const std::vector<ftype>&>())
    .def(py::init<const Atom_list&, const Atom_list&>())
    .def("rtop_frac", &RTop_orth::rtop_frac)
    .def("inverse", &RTop_orth::inverse)
    .def("axis_coordinate_near", &RTop_orth::axis_coordinate_near)
    .def("screw_translation", &RTop_orth::screw_translation)
    .def_static("identity", &RTop_orth::identity)
    //.def_static("null", &RTop_orth::null)
    .def("format", [](const RTop<>& self){ return self.format(); })
    .def("__repr__", [](const RTop<>& self) { return std::string(self.format()); })
    .def("__str__", [](const RTop<>& self) { return std::string(self.format()); })
    // base class methods
    .def("equals", [](const RTop_orth& self, const RTop_orth& other, const ftype& tol) { return self.equals(other, tol); })
    .def_property("rot",
    [](const RTop_orth& self) { return self.rot(); },
    [](RTop_orth& self, const Mat33<ftype>& rot) { self.rot() = rot; }
    )
    .def_property("trn",
         [](const RTop_orth& self) { return self.trn(); },
         [](RTop_orth& self, const Vec3<>& trn) { self.trn() = trn; }
    )
    .def("is_null", [](const RTop_orth& self) { return self.is_null(); })
    ;

py::class_<HKL, Vec3<int>>(m, "HKL")
    .def(py::init<>())
    .def(py::init<const Vec3<int>&>())
    .def(py::init<const int&, const int&, const int&>())
    // from numpy array
    .def(py::init([](py::array_t<int> hkl)
    {
        check_numpy_array_shape(hkl, {3}, true);
        auto buf = hkl.request();
        int* ptr = (int*)buf.ptr;
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
    .def("__repr__", [](const HKL& self) { return std::string(self.format()); })
    .def("__str__", [](const HKL& self) { return std::string(self.format()); })
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
    .def(py::init([](py::array_t<int> xyz_reci)
    {
        check_numpy_array_shape(xyz_reci, {3}, true);
        auto buf = xyz_reci.request();
        int* ptr = (int*)buf.ptr;
        return std::unique_ptr<Coord_reci_orth>(new Coord_reci_orth(ptr[0], ptr[1], ptr[2]));
    }))
    .def_property_readonly("xs", &Coord_reci_orth::xs)
    .def_property_readonly("ys", &Coord_reci_orth::ys)
    .def_property_readonly("zs", &Coord_reci_orth::zs)
    .def_property_readonly("xyz_reci",
        [](const Coord_reci_orth& self) { return array_as_numpy_1d<Coord_reci_orth, ftype>(self, 3); }
    )
    .def("invresolsq", &Coord_reci_orth::invresolsq)
    .def("coord_reci_frac", &Coord_reci_orth::coord_reci_frac)
    .def("transform", &Coord_reci_orth::transform)
    .def("__repr__", [](const Coord_reci_orth& self) { return std::string(self.format()); })
    .def("__str__", [](const Coord_reci_orth& self) { return std::string(self.format()); })
    ;

py::class_<Coord_reci_frac, Vec3<>>(m, "Coord_reci_frac")
    .def(py::init<>())
    .def(py::init<const Vec3<>&>())
    .def(py::init<const ftype&, const ftype&, const ftype&>())
    .def(py::init<const HKL&>())
    // from numpy array
    .def(py::init([](py::array_t<int> uvw_reci)
    {
        check_numpy_array_shape(uvw_reci, {3}, true);
        auto buf = uvw_reci.request();
        int* ptr = (int*)buf.ptr;
        return std::unique_ptr<Coord_reci_frac>(new Coord_reci_frac(ptr[0], ptr[1], ptr[2]));
    }))
    .def("hkl", &Coord_reci_frac::hkl)
    .def("invresolsq", &Coord_reci_frac::invresolsq)
    .def_property_readonly("us", &Coord_reci_frac::us)
    .def_property_readonly("vs", &Coord_reci_frac::vs)
    .def_property_readonly("ws", &Coord_reci_frac::ws)
    .def_property_readonly("uvw_reci",
        [](const Coord_reci_frac& self) { return array_as_numpy_1d<Coord_reci_frac, ftype>(self, 3); }
    )
    .def("coord_reci_orth", &Coord_reci_frac::coord_reci_orth)
    .def("transform", &Coord_reci_frac::transform)
    .def("__repr__", [](const Coord_reci_frac& self) { return std::string(self.format()); })
    .def("__str__", [](const Coord_reci_frac& self) { return std::string(self.format()); })
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
        auto buf = uvw.request();
        int* ptr = (int*)buf.ptr;
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
    .def("__repr__", [](const Coord_grid& self) { return std::string(self.format()); })
    .def("__str__", [](const Coord_grid& self) { return std::string(self.format()); })
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
        auto buf = xyz.request();
        int* ptr = (int*)buf.ptr;
        return std::unique_ptr<Coord_orth>(new Coord_orth(ptr[0], ptr[1], ptr[2]));
    }))
    .def_property_readonly("x", &Coord_orth::x)
    .def_property_readonly("y", &Coord_orth::y)
    .def_property_readonly("z", &Coord_orth::z)
    .def_property_readonly("xyz",
        [](const Coord_orth& self) { return array_as_numpy_1d<Coord_orth, ftype>(self, 3); }
    )
    .def("lengthsq", &Coord_orth::lengthsq)
    .def("coord_frac", &Coord_orth::coord_frac)
    .def("transform", &Coord_orth::transform)
    .def("__repr__", [](const Coord_orth& self) { return std::string(self.format()); })
    .def("__str__", [](const Coord_orth& self) { return std::string(self.format()); })
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
    // from numpy array
    .def(py::init([](py::array_t<ftype> uvw)
    {
        check_numpy_array_shape(uvw, {3}, true);
        auto buf = uvw.request();
        int* ptr = (int*)buf.ptr;
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
    .def("format", &Coord_frac::format)
    .def("__str__", [](const Coord_frac& self) { return self.format().c_str(); })
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
        auto buf = uvw.request();
        int* ptr = (int*)buf.ptr;
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
    .def("format", &Coord_map::format)
    .def("__str__", [](const Coord_map& self) { return self.format().c_str(); })
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
        check_numpy_array_shape(vals, {6}, true);
        auto buf = vals.request();
        ftype* ptr = (ftype*)buf.ptr;
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
        auto buf = vals.request();
        ftype* ptr = (ftype*)buf.ptr;
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
    .def("format", &Grid::format)
    .def("__str__", [](const Grid& self) { return self.format().c_str(); })
    ;

py::class_<Grid_sampling, Grid>(m, "Grid_sampling")
    .def(py::init<>())
    .def(py::init<const int&, const int&, const int&>())
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
    .def("is_null", &Grid_sampling::is_null)
    ;

py::class_<HKL_sampling>(m, "HKL_sampling")
    .def(py::init<>())
    .def(py::init<const Cell&, const Resolution&>())
    .def("hkl_limit", &HKL_sampling::hkl_limit)
    .def("resolution", &HKL_sampling::resolution)
    .def("in_resolution", &HKL_sampling::in_resolution)
    .def("is_null", &HKL_sampling::is_null)
    .def("format", &HKL_sampling::format)
    .def("__str__", [](const HKL_sampling& self) { return self.format().c_str(); })
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
        [](const Atom& self) { return self.element().c_str(); },
        [](Atom& self, const std::string& e) { self.set_element( String(e)); }
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
    .def("is_null", &Atom::is_null)
    ;


py::class_<Atom_list>(m, "Atom_list")
    .def(py::init<>())
    .def(py::init<const std::vector<Atom>&>())
    // from arrays of atom properties
    .def(py::init([](int n, const std::vector<std::string>& elements, py::array_t<ftype> coords,
                     py::array_t<ftype> u_anisos, py::array_t<ftype> occupancies, py::array_t<ftype> u_isos)
                     -> std::unique_ptr<Atom_list>
    {
        if (elements.size() != n) {
            throw std::logic_error("Element array length does not match the number of atoms!");
            return 0;
        }
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

    ;







} //init_coords
