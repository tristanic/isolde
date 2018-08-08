#include <pybind11/pybind11.h>

#include <clipper_ext/xtal_mgr.h>

namespace py=pybind11;
namespace cx=clipper_cx;
using namespace clipper;
using namespace clipper::datatypes;

void declare_xmap_container(py::module& m)
{
    using Class=cx::Xmap_details;
    py::class_<Class>(m, "Xmap_container")
        .def_property_readonly("xmap", (const Xmap<ftype32>& (Class::*)() const) &Class::xmap,
            py::return_value_policy::reference_internal)
        .def_property_readonly("stats", (const Map_stats& (Class::*)() const) &Class::map_stats)
        .def_property_readonly("is_difference_map", &Class::is_difference_map)
        .def_property_readonly("exclude_free_reflections", &Class::exclude_free_reflections)
        .def_property_readonly("fill_with_fcalc", &Class::fill_with_fcalc)
        .def_property_readonly("b_sharp", &Class::b_sharp)
        ;

} // declare_xmap_container

void declare_xtal_mgr(py::module& m)
{
    using Class=cx::Xtal_mgr_base;
    py::class_<Class>(m, "Xtal_mgr")
        .def(py::init<const HKL_info&, const HKL_data<Flag>&, const Grid_sampling&,
            const HKL_data<F_sigF<ftype32>>&>() )
        .def_property_readonly("free_flag", &Class::freeflag)
        .def_property_readonly("f_calc", &Class::fcalc)
        .def_property_readonly("usage_flags", &Class::usage)
        .def_property_readonly("base_2fofc", &Class::base_2fofc)
        .def_property_readonly("base_fofc", &Class::base_fofc)
        .def_property_readonly("weights", &Class::weights)
        .def_property_readonly("rwork", &Class::rwork)
        .def_property_readonly("rfree", &Class::rfree)
        .def_property_readonly("weighted_rwork", &Class::weighted_rwork)
        .def_property_readonly("weighted_rfree", &Class::weighted_rfree)
        .def("generate_fcalc", &Class::generate_fcalc)
        .def("generate_base_map_coeffs", &Class::generate_base_map_coeffs)
        .def("add_xmap", &Class::add_xmap,
            py::arg("name"), py::arg("base_coeffs"), py::arg("b_sharp"),
            py::arg("is_difference_map")=true, py::arg("exclude_free_reflections")=true,
            py::arg("fill_with_fcalc")=true)
        .def("recalculate_map", &Class::recalculate_map)
        // Get a reference to the managed xmap of a given name
        .def("get_xmap", &Class::get_xmap, py::return_value_policy::reference_internal)
        .def("delete_xmap", &Class::delete_xmap)
        ;
}

void init_xtal_mgr(py::module& m)
{
    declare_xmap_container(m);
    declare_xtal_mgr(m);
} // init_xtal_mgr
