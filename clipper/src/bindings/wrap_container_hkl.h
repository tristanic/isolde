#pragma once

#include <pybind11/pybind11.h>

#include <clipper/clipper.h>

namespace py=pybind11;
using namespace clipper;

template <class T>
void declare_CHKL_data(py::module &m, const std::string &class_str)
{
    using Class=CHKL_data<T>;
    std::string pyclass_name = std::string("CHKL_data_" + class_str);
    py::class_<Class, Container, HKL_data<T> >(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<Container&, const String>())
        .def("init", (void (Class::*)(const HKL_info&, const Cell&)) &Class::init)
        .def("init", (void (Class::*)(const Spacegroup&, const Cell&, const HKL_sampling&)) &Class::init)
        .def("update", &Class::update)
        .def("copy_from", [](Class& self, const HKL_data<T>& other) { self = other; })
        .def("set_all_values_to", [](Class& self, const T& value) { self = value; });
}
