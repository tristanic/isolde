#include <pybind11/pybind11.h>

#include "type_conversions.h"
#include <clipper/clipper.h>

namespace py=pybind11;
using namespace clipper;

void declare_map_stats(py::module& m)
{
    py::class_<Map_stats>(m, "Map_stats")
        .def(py::init<Xmap<ftype32>>())
        .def(py::init<Xmap<ftype64>>())
        .def(py::init<NXmap<ftype32>>())
        .def(py::init<NXmap<ftype64>>())
        .def_property_readonly("mean", &Map_stats::mean)
        .def_property_readonly("std_dev", &Map_stats::std_dev)
        .def_property_readonly("min", &Map_stats::min)
        .def_property_readonly("max", &Map_stats::max)
        .def_property_readonly("range", &Map_stats::range)
        ;
}



void init_map_utils(py::module& m)
{
    declare_map_stats(m);
} // init_map_utils
