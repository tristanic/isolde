#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "type_conversions.h"
#include <clipper/clipper.h>

namespace py=pybind11;
using namespace clipper;

void init_containers(py::module &m) {
    py::class_<Container>(m, "Container")
        .def(py::init<const String>())
        .def(py::init<Container&, const String&>())
        .def("update", &Container::update)
        .def("path", &Container::path)
        .def("name", &Container::name)
        .def("set_name", &Container::set_name)
        .def("is_destroyed_with_parent", &Container::is_destroyed_with_parent)
        .def("set_destroyed_with_parent", &Container::set_destroyed_with_parent)
        .def("move", &Container::move)
        .def("has_parent", &Container::has_parent)
        // .def("parent", (const Container& (Container::*)() const)&Container::parent)
        .def("parent", (Container& (Container::*)())&Container::parent)
        .def("num_childrent", &Container::num_children)
        // .def("child", (const Container& (Container::*)(const int&) const)&Container::child)
        .def("child", (Container& (Container::*)(const int&))&Container::child)
        // .def("ultimate_parent", (const Container& (Container::*)() const)&Container::ultimate_parent)
        .def("ultimate_parent", (Container& (Container::*)())&Container::ultimate_parent);

    py::class_<CHKL_info, Container, HKL_info>(m, "CHKL_info")
        .def(py::init<const String, const Spacegroup&, const Cell&,
                      const Resolution&, const bool&>())
        .def(py::init<Container&, const String, const bool&>())
        .def("init", &CHKL_info::init)
        .def("generate_hkl_list", &CHKL_info::generate_hkl_list)
        .def("update", &CHKL_info::update);

} // init_containers
