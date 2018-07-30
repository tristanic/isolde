#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "type_conversions.h"
#include <clipper/clipper.h>

#include "numpy_helper.h"

namespace py=pybind11;
using namespace clipper;

template<class T>
void declare_grad_orth(py::module& m, const char* dtype)
{
    using Base=Vec3<T>;
    using Class=Grad_orth<T>;
    std::string pyclass_name = std::string("Grad_orth_") + dtype;
    py::class_<Class, Base>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const Base&>())
        .def(py::init<const T&, const T&, const T&>())
        .def("dx", &Class::dx)
        .def("dy", &Class::dy)
        .def("dz", &Class::dz)
        // as numpy
        .def("dxyz", [](const Class& self) { return array_as_numpy_1d<Class, T>(self, 3); })
        .def("grad_frac", &Class::grad_frac)
        .def("format", &Class::format)
        .def("__str__", [](const Class& self) { return self.format().c_str(); })
        ;
}

template<class T>
void declare_grad_frac(py::module& m, const char* dtype)
{
    using Base=Vec3<T>;
    using Class=Grad_frac<T>;
    std::string pyclass_name = std::string("Grad_frac_") + dtype;
    py::class_<Class, Base>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const Base&>())
        .def(py::init<const T&, const T&, const T&>())
        .def("du", &Class::du)
        .def("dv", &Class::dv)
        .def("dw", &Class::dw)
        // as numpy
        .def("duvw", [](const Class& self) { return array_as_numpy_1d<Class, T>(self, 3); })
        .def("grad_orth", &Class::grad_orth)
        .def("grad_map", &Class::grad_map)
        .def("format", &Class::format)
        .def("__str__", [](const Class& self) { return self.format().c_str(); })
        ;
}

template<class T>
void declare_grad_map(py::module& m, const char* dtype)
{
    using Base=Vec3<T>;
    using Class=Grad_map<T>;
    std::string pyclass_name = std::string("Grad_map_") + dtype;
    py::class_<Class, Base>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const Base&>())
        .def(py::init<const T&, const T&, const T&>())
        .def("du", &Class::du)
        .def("dv", &Class::dv)
        .def("dw", &Class::dw)
        // as numpy
        .def("duvw", [](const Class& self) { return array_as_numpy_1d<Class, T>(self, 3); })
        .def("grad_frac", &Class::grad_frac)
        .def("format", &Class::format)
        .def("__str__", [](const Class& self) { return self.format().c_str(); })
        ;
}

template<class T>
void declare_curv_orth(py::module& m, const char* dtype)
{
    using Base=Mat33<T>;
    using Class=Curv_orth<T>;
    std::string pyclass_name = std::string("Curv_orth_") + dtype;
    py::class_<Class, Base>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const Base&>())
        .def("curv_frac", &Class::curv_frac)
        ;
}

template<class T>
void declare_curv_frac(py::module& m, const char* dtype)
{
    using Base=Mat33<T>;
    using Class=Curv_frac<T>;
    std::string pyclass_name = std::string("Curv_frac_") + dtype;
    py::class_<Class, Base>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const Base&>())
        .def("curv_orth", &Class::curv_orth)
        .def("curv_map", &Class::curv_map)
        ;
}

template<class T>
void declare_curv_map(py::module& m, const char* dtype)
{
    using Base=Mat33<T>;
    using Class=Curv_map<T>;
    std::string pyclass_name = std::string("Curv_map_") + dtype;
    py::class_<Class, Base>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const Base&>())
        .def("curv_frac", &Class::curv_frac)
        ;
}



void init_derivs(py::module& m) {
    declare_grad_orth<float>(m, "float");
    declare_grad_frac<float>(m, "float");
    declare_grad_map<float>(m, "float");
    declare_curv_orth<float>(m, "float");
    declare_curv_frac<float>(m, "float");
    declare_curv_map<float>(m, "float");

    declare_grad_orth<double>(m, "double");
    declare_grad_frac<double>(m, "double");
    declare_grad_map<double>(m, "double");
    declare_curv_orth<double>(m, "double");
    declare_curv_frac<double>(m, "double");
    declare_curv_map<double>(m, "double");

} // init_derivs
