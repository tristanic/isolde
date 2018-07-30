#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "type_conversions.h"
#include <clipper/clipper.h>

namespace py=pybind11;
using namespace clipper;

void init_hkl_info(py::module &m) {

    py::class_<HKL_info>(m, "HKL_info")
        .def(py::init<>())
        .def(py::init<const Spacegroup&, const Cell&, const Resolution&, const bool&>())
        .def("init", (void (HKL_info::*)(const Spacegroup&, const Cell&, const Resolution&, const bool&)) &HKL_info::init)
        .def("init", (void (HKL_info::*)(const Spacegroup&, const Cell&, const HKL_sampling&, const bool&)) &HKL_info::init)
        .def_property_readonly("is_null", &HKL_info::is_null)
        .def_property_readonly("cell", &HKL_info::cell)
        .def_property_readonly("spacegroup", &HKL_info::spacegroup)
        .def_property_readonly("hkl_sampling", &HKL_info::hkl_sampling)
        .def_property_readonly("resolution", &HKL_info::resolution)
        .def("generate_hkl_list", &HKL_info::generate_hkl_list)
        .def("add_hkl_list", &HKL_info::add_hkl_list)
        .def_property_readonly("num_reflections", &HKL_info::num_reflections)
        .def("hkl_of", &HKL_info::hkl_of)
        .def("index_of", &HKL_info::index_of)
        .def("invresolsq", &HKL_info::invresolsq)
        .def_property_readonly("invresolsq_range", &HKL_info::invresolsq_range)
        .def("hkl_class", &HKL_info::hkl_class)
        .def("find_sym", &HKL_info::find_sym)
        .def_property_readonly("first", &HKL_info::first);

    py::class_<HKL_info::HKL_reference_base>(m, "_HKL_reference_base")
        .def_property_readonly("base_hkl_info", &HKL_info::HKL_reference_base::base_hkl_info)
        .def_property_readonly("index", &HKL_info::HKL_reference_base::index)
        .def("invresolsq", (ftype (HKL_info::HKL_reference_base::*)(const HKL_data_base&) const)&HKL_info::HKL_reference_base::invresolsq)
        .def("invresolsq", (ftype (HKL_info::HKL_reference_base::*)() const)&HKL_info::HKL_reference_base::invresolsq)
        .def("last", &HKL_info::HKL_reference_base::last);

    py::class_<HKL_info::HKL_reference_index, HKL_info::HKL_reference_base>(m, "HKL_reference_index")
        .def(py::init<>())
        .def(py::init<const HKL_info&, const int&>())
        .def_property_readonly("hkl", &HKL_info::HKL_reference_index::hkl)
        .def_property_readonly("hkl_class", &HKL_info::HKL_reference_index::hkl_class)
        // avoid creating a new Python object every time we increment
        .def("next", [](HKL_info::HKL_reference_index& self) { self.next(); });

    py::class_<HKL_info::HKL_reference_coord, HKL_info::HKL_reference_base>(m, "HKL_reference_coord")
        .def(py::init<>())
        .def(py::init<const HKL_info&, const HKL&>())
        .def_property("hkl",
            &HKL_info::HKL_reference_coord::hkl,
            [](HKL_info::HKL_reference_coord& self, const HKL& hkl__) { self.set_hkl(hkl__); }
        )
        .def_property_readonly("sym", &HKL_info::HKL_reference_coord::sym)
        .def_property_readonly("friedel", &HKL_info::HKL_reference_coord::friedel)
        // avoid creating new Python objects when incrementing
        .def("next", [](HKL_info::HKL_reference_coord& self) { self.next(); })
        .def("next_h", [](HKL_info::HKL_reference_coord& self) { self.next_h(); })
        .def("next_k", [](HKL_info::HKL_reference_coord& self) { self.next_k(); })
        .def("next_l", [](HKL_info::HKL_reference_coord& self) { self.next_l(); })
        .def("prev_h", [](HKL_info::HKL_reference_coord& self) { self.prev_h(); })
        .def("prev_k", [](HKL_info::HKL_reference_coord& self) { self.prev_k(); })
        .def("prev_l", [](HKL_info::HKL_reference_coord& self) { self.prev_l(); })
        ;



} //init_chkl_info
