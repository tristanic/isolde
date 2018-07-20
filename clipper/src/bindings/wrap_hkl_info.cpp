#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <clipper/clipper.h>

namespace py=pybind11;
using namespace clipper;

void init_hkl_info(py::module &m) {

    py::class_<HKL_info>(m, "HKL_info")
        .def(py::init<>())
        .def(py::init<const Spacegroup&, const Cell&, const Resolution&, const bool&>())
        .def("init", (void (HKL_info::*)(const Spacegroup&, const Cell&, const Resolution&, const bool&)) &HKL_info::init)
        .def("init", (void (HKL_info::*)(const Spacegroup&, const Cell&, const HKL_sampling&, const bool&)) &HKL_info::init)
        .def("is_null", &HKL_info::is_null)
        .def("cell", &HKL_info::cell)
        .def("spacegroup", &HKL_info::spacegroup)
        .def("hkl_sampling", &HKL_info::hkl_sampling)
        .def("resolution", &HKL_info::resolution)
        .def("generate_hkl_list", &HKL_info::generate_hkl_list)
        .def("add_hkl_list", &HKL_info::add_hkl_list)
        .def("num_reflections", &HKL_info::num_reflections)
        .def("hkl_of", &HKL_info::hkl_of)
        .def("index_of", &HKL_info::index_of)
        .def("invresolsq", &HKL_info::invresolsq)
        .def("invresolsq_range", &HKL_info::invresolsq_range)
        .def("hkl_class", &HKL_info::hkl_class)
        .def("find_sym", &HKL_info::find_sym)
        .def("first", &HKL_info::first);

    py::class_<HKL_info::HKL_reference_base>(m, "HKL_reference_base")
        .def("base_hkl_info", &HKL_info::HKL_reference_base::base_hkl_info)
        .def("index", &HKL_info::HKL_reference_base::index)
        .def("invresolsq", (ftype (HKL_info::HKL_reference_base::*)(const HKL_data_base&) const)&HKL_info::HKL_reference_base::invresolsq)
        .def("invresolsq", (ftype (HKL_info::HKL_reference_base::*)() const)&HKL_info::HKL_reference_base::invresolsq)
        .def("last", &HKL_info::HKL_reference_base::last);

    py::class_<HKL_info::HKL_reference_index, HKL_info::HKL_reference_base>(m, "HKL_reference_index")
        .def(py::init<>())
        .def(py::init<const HKL_info&, const int&>())
        .def("hkl", &HKL_info::HKL_reference_index::hkl)
        .def("hkl_class", &HKL_info::HKL_reference_index::hkl_class)
        .def("next", &HKL_info::HKL_reference_index::next);

    py::class_<HKL_info::HKL_reference_coord, HKL_info::HKL_reference_base>(m, "HKL_reference_coord")
        .def(py::init<>())
        .def(py::init<const HKL_info&, const HKL&>())
        .def("hkl", &HKL_info::HKL_reference_coord::hkl)
        .def("sym", &HKL_info::HKL_reference_coord::sym)
        .def("friedel", &HKL_info::HKL_reference_coord::friedel)
        .def("set_hkl", &HKL_info::HKL_reference_coord::set_hkl)
        .def("next", &HKL_info::HKL_reference_coord::next)
        .def("next_h", &HKL_info::HKL_reference_coord::next_h)
        .def("next_k", &HKL_info::HKL_reference_coord::next_k)
        .def("next_l", &HKL_info::HKL_reference_coord::next_l)
        .def("prev_h", &HKL_info::HKL_reference_coord::prev_h)
        .def("prev_k", &HKL_info::HKL_reference_coord::prev_k)
        .def("prev_l", &HKL_info::HKL_reference_coord::prev_l);



} //init_chkl_info
