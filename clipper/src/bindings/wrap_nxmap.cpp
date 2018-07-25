#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>


#include <clipper/clipper.h>
#include "numpy_helper.h"

namespace py=pybind11;
using namespace clipper;

void declare_nxmap_reference_index(py::module& m)
{
    using Class = NXmap_base::Map_reference_index;
    py::class_<Class>(m, "NXmap_reference_index")
        .def(py::init<>())
        .def(py::init<const NXmap_base&>())
        .def(py::init<const NXmap_base&, const Coord_grid&>())
        .def("coord", &Class::coord)
        .def("coord_orth", &Class::coord_orth)
        .def("set_coord", &Class::set_coord)
        .def("next", &Class::next)
        .def("index_offset", &Class::index_offset)
        // base class methods
        .def_property_readonly("base_nxmap", [](const Class& self) { return self.base_nxmap(); })
        .def_property_readonly("index", [](const Class& self) { return self.index(); })
        .def("last", [](const Class& self){ return self.last(); })
        ;
}

void declare_nxmap_reference_coord(py::module& m)
{
    using Class = NXmap_base::Map_reference_coord;
    py::class_<Class>(m, "NXmap_reference_coord")
        .def(py::init<>())
        .def(py::init<const NXmap_base&>())
        .def(py::init<const NXmap_base&, const Coord_grid&>())
        .def("coord", &Class::coord)
        .def("coord_orth", &Class::coord_orth)
        .def("set_coord", &Class::set_coord)
        .def("next", &Class::next)
        .def("next_u", &Class::next_u)
        .def("next_v", &Class::next_v)
        .def("next_w", &Class::next_w)
        .def("prev_u", &Class::prev_u)
        .def("prev_v", &Class::prev_v)
        .def("prev_w", &Class::prev_w)
        // base class methods
        .def_property_readonly("base_nxmap", [](const Class& self) { return self.base_nxmap(); })
        .def_property_readonly("index", [](const Class& self) { return self.index(); })
        .def("last", [](const Class& self){ return self.last(); })
        ;
}

template<class Derived, class T>
void apply_nxmap_base_methods(py::class_<Derived>& pyclass)
{
    pyclass
        .def("is_null", [](const Derived& self) { return self.is_null(); })
        .def("grid", [](const Derived& self) { return self.grid(); })
        .def("operator_orth_grid", [](const Derived& self){ return self.operator_orth_grid(); })
        .def("operator_grid_orth", [](const Derived& self) { return self.operator_grid_orth(); })
        .def("coord_orth", [](const Derived& self, const Coord_map& cm) { return self.coord_orth(cm); })
        .def("coord_map", [](const Derived& self, const Coord_orth& co) { return self.coord_map(co); })
        .def("in_map", [](const Derived& self, const Coord_grid& pos) { return self.in_map(pos); })
        // TODO: handle in_map() for different interpolators
        .def("multiplicity", [](const Derived& self, const Coord_grid&) { return 1; })

        .def("first", [](const Derived& self) { return self.first(); })
        .def("first_coord", [](const Derived& self) { return self.first_coord(); })
        ;
}

template<class T>
void declare_nxmap(py::module& m, const char* dtype)
{
    using Class = NXmap<T>;
    auto pyclass_name = std::string("NXmap_") + dtype;
    py::class_<Class> nxmap(m, pyclass_name.c_str());
    nxmap
        .def(py::init<>())
        .def(py::init<const Grid&, const RTop<>&>())
        .def(py::init<const Cell&, const Grid_sampling&, const Grid_range&>())
        .def("init", (void (Class::*)(const Grid&, const RTop<>&)) &Class::init)
        .def("init", (void (Class::*)(const Cell&, const Grid_sampling&, const Grid_range&)) &Class::init)
        .def("__getitem__", [](const Class& self, const NXmap_base::Map_reference_index i) { return self[i]; })
        .def("__setitem__", [](Class& self, const NXmap_base::Map_reference_index i, const T& val){ self[i] = val; })
        .def("__getitem__", [](const Class& self, const NXmap_base::Map_reference_coord i) { return self[i]; })
        .def("__setitem__", [](Class& self, const NXmap_base::Map_reference_coord i, const T& val) { self[i] = val; })
        .def("get_data", [](const Class& self, const Coord_grid& pos) { return self.get_data(pos); })
        .def("set_data", [](Class& self, const Coord_grid& pos, const T& val) { self.set_data(pos, val); })
        .def("set_all_values_to", [](Class& self, const T& val) { self = val; })
        .def(py::self += py::self)
        .def(py::self -= py::self)
        ;
    apply_nxmap_base_methods<Class, T>(nxmap);
}

void init_nxmap(py::module& m)
{
    declare_nxmap_reference_index(m);
    declare_nxmap_reference_coord(m);
    declare_nxmap<ftype32>(m, "float");
    declare_nxmap<ftype64>(m, "double");
}
