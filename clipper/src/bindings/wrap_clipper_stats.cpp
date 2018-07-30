#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

#include "type_conversions.h"
#include <clipper/clipper.h>

#include "numpy_helper.h"

namespace py=pybind11;
using namespace clipper;

template<class T>
void init_range(py::module& m, const std::string& dtype)
{
    using Class=Range<T>;
    auto pyclass_name = std::string("Range_") + dtype;
    py::class_<Class>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const T&, const T&>())
        .def_property_readonly("min", &Class::min)
        .def_property_readonly("max", &Class::max)
        .def_property_readonly("range", &Class::range)
        .def("include", &Class::include)
        .def("contains", &Class::contains)
        .def("truncate", &Class::truncate)
        ;

} // init_range

void init_clipper_stats(py::module &m)
{
init_range<int>(m, "int");
init_range<float>(m, "float");
init_range<double>(m, "double");

py::class_<Range_sampling, Range<ftype>>(m, "Range_sampling")
    .def(py::init<>())
    .def(py::init<const int&>())
    .def(py::init<const Range<ftype>&, const int&>())
    .def("indexf", &Range_sampling::indexf)
    .def("x", (ftype (Range_sampling::*)(const ftype&) const) &Range_sampling::x)
    .def("x", (ftype (Range_sampling::*)(const int&) const) &Range_sampling::x)
    .def("index", &Range_sampling::index)
    .def("index_bounded", &Range_sampling::index_bounded)
    .def("x_min", &Range_sampling::x_min)
    .def("x_max", &Range_sampling::x_max)
    .def_property_readonly("size", &Range_sampling::size)
    ;

py::class_<Histogram, Range_sampling>(m, "Histogram")
    .def(py::init<>())
    .def(py::init<const Range<ftype>&, const int&>())
    .def("accumulate", (void (Histogram::*)(const ftype&)) &Histogram::accumulate)
    .def("accumulate", (void (Histogram::*)(const ftype&, const ftype&)) &Histogram::accumulate)
    .def_property_readonly("sum", &Histogram::sum)
    .def("y", (const ftype& (Histogram::*)(const int&) const) &Histogram::y)
    .def("y", (ftype (Histogram::*)(const ftype&) const) &Histogram::y)
    .def(py::self += py::self)
    //.def("__iadd__", [](Histogram& self, const Histogram& other) { self += other; })
    ;

py::class_<Generic_ordinal>(m, "Generic_ordinal")
    .def(py::init<>())
    .def(py::init<const Range<ftype>&, const int&>())
    .def("init", (void (Generic_ordinal::*)(const Range<ftype>&, const int)) &Generic_ordinal::init)
    .def("init", (void (Generic_ordinal::*)(const std::vector<ftype>&, const int)) &Generic_ordinal::init)
    .def("ordinal", &Generic_ordinal::ordinal)
    .def("accumulate", (void (Generic_ordinal::*)(const ftype&)) &Generic_ordinal::accumulate)
    .def("accumulate", (void (Generic_ordinal::*)(const ftype&, const ftype&)) &Generic_ordinal::accumulate)
    .def("prep_ordinal", &Generic_ordinal::prep_ordinal)
    .def("invert", &Generic_ordinal::invert)
    ;
} // init_clipper_stats
