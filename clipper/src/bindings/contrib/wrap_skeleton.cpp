#include <pybind11/pybind11.h>

#include "../type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

namespace py=pybind11;
using namespace clipper;
using namespace clipper::datatypes;


template <class T>
void declare_skeleton_base(py::module &m, const char* dtype)
{
    using Class=Skeleton_base<int, T>;
    auto pyclass_name = std::string("_Skeleton_base_") + dtype;
    py::class_<Class, std::unique_ptr<Class, py::nodelete>>(m, pyclass_name.c_str())
        .def("__call__",
            [](const Class& self, Xmap<int>& xskl, const Xmap<T>& xmap)
            { return self(xskl, xmap); })
        ;
} // declare_skeleton_base

template<class T>
void declare_skeleton_fast(py::module &m, const char* dtype)
{
    using Class=Skeleton_fast<int, T>;
    auto pyclass_name = std::string("Skeleton_") + dtype;
    py::class_<Class, Skeleton_base<int, T>>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<Xmap<int>&, const Xmap<T>&>())
        ;
}

template<class T>
void declare_skeleton_neighbours(py::module& m, const char* dtype)
{
    using Class=typename Skeleton_fast<int, T>::Neighbours;
    auto pyclass_name=std::string("Skeleton_neighbours_")+dtype;
    py::class_<Class>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const Xmap_base&, const float, const float>(),
            py::arg("map"), py::arg("min_dist")=0.5, py::arg("max_dist")=2.5)
        .def("__getitem__", [](const Class& self, int i) { return self[i]; })
        .def_property_readonly("size", &Class::size)
        .def("__len__", &Class::size)
        ;
}

void init_skeleton(py::module& m)
{
    declare_skeleton_base<ftype32>(m, "float");
    declare_skeleton_fast<ftype32>(m, "float");
    declare_skeleton_neighbours<ftype32>(m, "float");

    declare_skeleton_base<ftype64>(m, "double");
    declare_skeleton_fast<ftype64>(m, "double");
    declare_skeleton_neighbours<ftype64>(m, "double");

} // init_sfcalc
