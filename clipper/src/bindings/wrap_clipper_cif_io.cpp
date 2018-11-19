#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/clipper-cif.h>


namespace py=pybind11;
using namespace clipper;

void declare_cif_io(py::module& m)
{
    py::class_<CIFfile> cif_file(m, "CIFfile");
    cif_file
        .def(py::init<>())
        .def("open_read", &CIFfile::open_read)
        .def("close_read", &CIFfile::close_read)
        .def_property_readonly("spacegroup", &CIFfile::spacegroup)
        .def_property_readonly("cell", &CIFfile::cell)
        .def_property_readonly("hkl_sampling", &CIFfile::hkl_sampling)
        .def_property_readonly("resolution", [](const CIFfile& self)
        {
            return self.resolution();
        })
        .def("import_hkl_info", &CIFfile::import_hkl_info)
        .def("import_hkl_data", &CIFfile::import_hkl_data)
        .def_property_readonly("contains_phases_predicate", &CIFfile::contains_phases_p)
        ;
} // declare_cif_file

void init_cif_io(py::module& m)
{
    declare_cif_io(m);
} // init_cif_file
