#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <clipper/clipper.h>


namespace py=pybind11;
using namespace clipper;


void init_hkl_data(py::module &m) {

    py::class_<HKL_data_base, std::unique_ptr<HKL_data_base, py::nodelete>>(m, "HKL_data_base")
        .def("init", (void (HKL_data_base::*)(const HKL_info&, const Cell&)) &HKL_data_base::init)
        .def("init", (void (HKL_data_base::*)(const HKL_data_base&)) &HKL_data_base::init)
        .def("init", (void (HKL_data_base::*)(const Spacegroup&, const Cell&, const HKL_sampling&)) &HKL_data_base::init)
        .def("is_null", &HKL_data_base::is_null)
        .def("base_hkl_info", &HKL_data_base::base_hkl_info)
        .def("base_cell", &HKL_data_base::base_cell)
        .def("spacegroup", &HKL_data_base::spacegroup)
        .def("cell", &HKL_data_base::cell)
        .def("resolution", &HKL_data_base::resolution)
        .def("hkl_sampling", &HKL_data_base::hkl_sampling)
        .def("hkl_info", &HKL_data_base::hkl_info)
        .def("invresolsq", &HKL_data_base::invresolsq)
        .def("invresolsq_range", &HKL_data_base::invresolsq_range)
        .def("num_obs", &HKL_data_base::num_obs)
        .def("update", &HKL_data_base::update)
        .def("type", &HKL_data_base::type)
        .def("missing", &HKL_data_base::missing)
        .def("set_null", &HKL_data_base::set_null)
        .def("data_size", &HKL_data_base::data_size)
        .def("data_names", &HKL_data_base::data_names)
        .def("data_export", &HKL_data_base::data_export)
        .def("data_import", &HKL_data_base::data_import)
        .def("mask", &HKL_data_base::mask)
        .def("first", &HKL_data_base::first)
        .def("first_data", &HKL_data_base::first_data)
        .def("next_data", &HKL_data_base::next_data);

} // init_hkl_data
