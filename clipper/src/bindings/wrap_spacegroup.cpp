#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "type_conversions.h"
#include <clipper/clipper.h>
namespace py=pybind11;
using namespace clipper;


void declare_spgr_descr(py::module& m)
{
    py::class_<Spgr_descr> spgr_descr(m, "Spgr_descr");
        spgr_descr
        .def(py::init<>())
        .def(py::init<const String&, Spgr_descr::TYPE>())
        .def(py::init<const int&>())
        .def(py::init<const Spgr_descr::Symop_codes&>())
        .def_property_readonly("spacegroup_number", &Spgr_descr::spacegroup_number)
        .def_property_readonly("symbol_hall", &Spgr_descr::symbol_hall)
        .def_property_readonly("symbol_hm", &Spgr_descr::symbol_hm)
        .def_property_readonly("symbol_xhm", &Spgr_descr::symbol_xhm)
        .def_property_readonly("symbol_hm_ext", &Spgr_descr::symbol_hm_ext)
        .def_static("set_preferred", &Spgr_descr::set_preferred)
        .def_property_readonly("generator_ops", &Spgr_descr::generator_ops)
        .def("__hash__", &Spgr_descr::hash)
        ;

    py::enum_<Spgr_descr::TYPE>(spgr_descr, "TYPE")
        .value("Hall", Spgr_descr::TYPE::Hall)
        .value("HM", Spgr_descr::TYPE::HM)
        .value("XHM", Spgr_descr::TYPE::XHM)
        .value("Symops", Spgr_descr::TYPE::Symops)
        .value("Number", Spgr_descr::TYPE::Number)
        .value("Unknown", Spgr_descr::TYPE::Unknown)
        ;
} //init_spgr_descr

void declare_symop_codes(py::module& m)
{
    using Class = Spgr_descr::Symop_codes;
    py::class_<Class>(m, "Symop_codes")
        .def(py::init<>())
        .def("init_hall", &Class::init_hall)
        .def("init_symops", &Class::init_symops)
        .def("expand", &Class::expand)
        .def_property_readonly("primitive_noninversion_ops", &Class::primitive_noninversion_ops)
        .def_property_readonly("inversion_ops", &Class::inversion_ops)
        .def_property_readonly("primitive_ops", &Class::primitive_ops)
        .def_property_readonly("centering_ops", &Class::centering_ops)
        .def_property_readonly("laue_ops", &Class::laue_ops)
        .def_property_readonly("point_group_ops", &Class::pgrp_ops)
        .def_property_readonly("pgrp_ops", &Class::pgrp_ops)
        .def_property_readonly("patterson_ops", &Class::patterson_ops)
        .def_property_readonly("generator_ops", &Class::generator_ops)
        .def("product", &Class::product)
        .def("__hash__", &Class::hash)
        ;
}  //init_symop_codes

void declare_spacegroup(py::module& m)
{
    py::class_<Spacegroup, Spgr_descr> spacegroup(m, "Spacegroup");
    spacegroup
        .def(py::init<>())
        .def(py::init<Spacegroup::TYPE>())
        .def(py::init<const Spgr_descr&>())
        .def("init", &Spacegroup::init)
        .def_property_readonly("is_null", &Spacegroup::is_null)
        .def_property_readonly("descr", &Spacegroup::descr)
        .def_property_readonly("num_symops", &Spacegroup::num_symops)
        .def_property_readonly("num_primitive_symops", &Spacegroup::num_primitive_symops )
        .def_property_readonly("num_centering_symops", &Spacegroup::num_centering_symops )
        .def_property_readonly("num_inversion_symops", &Spacegroup::num_inversion_symops )
        .def_property_readonly("num_primitive_noninversion_symops", &Spacegroup::num_primitive_noninversion_symops)
        .def("symop", &Spacegroup::symop)
        .def("primitive_symop", &Spacegroup::primitive_symop)
        .def("inversion_symop", &Spacegroup::inversion_symop)
        .def("centering_symop", &Spacegroup::centering_symop)
        .def("order_of_symmetry_about_axis", &Spacegroup::order_of_symmetry_about_axis)
        .def("hkl_class", &Spacegroup::hkl_class)
        .def("recip_asu", &Spacegroup::recip_asu)
        .def("product_op", &Spacegroup::product_op)
        .def("inverse_op", &Spacegroup::inverse_op)
        .def_property_readonly("asu_max", &Spacegroup::asu_max)
        .def_property_readonly("asu_min", &Spacegroup::asu_min)
        .def_property_readonly("invariant_under_change_of_hand", &Spacegroup::invariant_under_change_of_hand)
        .def_property_readonly("symbol_laue", [](const Spacegroup& self) { return self.symbol_laue().c_str(); })
        .def_static("p1", &Spacegroup::p1)
        .def_static("null", &Spacegroup::null)
        ;

    py::enum_<Spacegroup::TYPE>(spacegroup, "SGTYPE")
        .value("Null", Spacegroup::TYPE::Null)
        .value("P1", Spacegroup::TYPE::P1)
        ;

    py::enum_<Spacegroup::AXIS>(spacegroup, "AXIS")
        .value("A", Spacegroup::AXIS::A)
        .value("B", Spacegroup::AXIS::B)
        .value("C", Spacegroup::AXIS::C)
        ;

}

void init_spacegroup(py::module& m)
{
    declare_spgr_descr(m);
    declare_symop_codes(m);
    declare_spacegroup(m);
}
