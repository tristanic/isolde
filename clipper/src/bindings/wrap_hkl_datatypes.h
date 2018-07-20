#pragma once

#include <pybind11/pybind11.h>

#include <clipper/clipper.h>

namespace py=pybind11;
using namespace clipper;

template <class C>
void declare_HKL_data(py::module &m, const std::string &class_str)
{
    using Class=HKL_data<C>;
    std::string pyclass_name = std::string("HKL_data_") + class_str;
    py::class_<Class>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const HKL_info&>())
        .def(py::init<const HKL_info&, const Cell&>())
        .def(py::init<const Spacegroup&, const Cell&, const HKL_sampling&>())
        .def(py::init<const HKL_data_base&>())
        .def("init", (void (Class::*)(const HKL_info&, const Cell&)) &Class::init)
        .def("init", (void (Class::*)(const Spacegroup&, const Cell&, const HKL_sampling&)) &Class::init)
        .def("init", (void (Class::*)(const HKL_data_base&)) &Class::init)
        .def("update", &Class::update)
        .def("type", &Class::type)
        .def("missing", &Class::missing)
        .def("set_null", &Class::set_null)
        .def("data_size", &Class::data_size)
        .def("data_names", &Class::data_names)
        .def("data_export", &Class::data_export)
        .def("data_import", &Class::data_import)
        .def("mask", &Class::mask)
        .def("__getitem__", [](const Class& self, const HKL_info::HKL_reference_index& i){return self[i];}, py::is_operator())
        .def("__getitem__", [](const Class& self, const HKL_info::HKL_reference_coord& ih) {
            C data;
            if (self.get_data(ih, data))
                return data;
            throw std::out_of_range("No data equivalent to that HKL!");
        }, py::is_operator())
        .def("__setitem__", [](Class& self, const HKL_info::HKL_reference_coord& ih, const C& data) {
            if (!self.set_data(ih, data))
                throw std::out_of_range("No equivalent HKL has been indexed for this dataset!");
        }, py::is_operator())
        .def("__getitem__", [](const Class& self, const int& index) {return self[index];}, py::is_operator())
        .def("__getitem__", [](const Class& self, const HKL& hkl) {
            C data;
            if (self.get_data(hkl, data))
                return data;
            throw std::out_of_range("No data equivalent to that HKL!");
        }, py::is_operator())
        .def("__setitem__", [](Class& self, const HKL& hkl, const C& data) {
            if (!self.set_data(hkl, data))
                throw std::out_of_range("No equivalent HKL has been indexed for this dataset!");
        })
        .def("copy_from", [](Class& self, const Class& other) { self=other; })
        .def("set_all_values_to", [](Class& self, const C& value) { self=value; })
        /*
        Since the base class HKL_data_base has a protected virtual destructor,
        exposing it to Python via PyBind11 requires it (and all derived classes)
        to be wrapped with std::unique_ptr<Class, py::nodelete>, which would
        require the python side code to explicitly handle object deletion. That
        would be a serious pain, so instead we'll hide the base class from
        Python entirely and expose all base class functions using lambdas.
        */
        .def("is_null", [](HKL_data_base& b) { return b.is_null(); })
        .def("base_hkl_info", [](HKL_data_base& b) { return b.base_hkl_info(); })
        .def("base_cell", [](HKL_data_base& b) { return b.base_cell(); })
        .def("spacegroup", [](HKL_data_base& b) { return b.spacegroup(); })
        .def("cell", [](HKL_data_base& b) { return b.cell(); })
        .def("resolution", [](HKL_data_base& b) { return b.resolution(); })
        .def("hkl_sampling", [](HKL_data_base& b) { return b.hkl_sampling(); })
        .def("hkl_info", [](HKL_data_base& b) { return b.hkl_info(); })
        .def("invresolsq", [](HKL_data_base& b, const int& index) { return b.invresolsq(index); })
        .def("invresolsq_range", [](HKL_data_base& b) { return b.invresolsq_range(); })
        .def("num_obs", [](HKL_data_base& b) { return b.num_obs(); })
        .def("first", [](HKL_data_base& b) { return b.first(); })
        .def("first_data", [](HKL_data_base& b) { return b.first_data(); })
        .def("next_data", [](HKL_data_base& b, HKL_info::HKL_reference_index& ih) { return b.next_data(ih); });
}
