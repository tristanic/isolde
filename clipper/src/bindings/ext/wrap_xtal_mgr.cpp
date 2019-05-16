/**
 * @Author: Tristan Croll <tic20>
 * @Date:   14-Sep-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 16-May-2019
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */



#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
        .def_property_readonly("exclude_missing_reflections", &Class::exclude_missing)
        .def_property_readonly("exclude_free_reflections", &Class::exclude_free_reflections)
        .def_property_readonly("fill_with_fcalc", &Class::fill_with_fcalc)
        .def_property_readonly("b_sharp", &Class::b_sharp)
        .def_property_readonly("_coeffs", [](const Class& self) { return self.coeffs(); })
        .def_property_readonly("_base_coeffs", &Class::base_coeffs)
        ;

} // declare_xmap_container

void declare_xtal_mgr(py::module& m)
{
    using Class=cx::Xtal_mgr_base;
    py::class_<Class>(m, "Xtal_mgr")
        .def(py::init<const HKL_info&, const HKL_data<Flag>&, const Grid_sampling&,
            const HKL_data<F_sigF<ftype32>>&>() )
        .def_property_readonly("free_flag", &Class::freeflag)
        .def_property_readonly("flags", &Class::flags)
        .def_property_readonly("f_obs", &Class::fobs)
        .def_property_readonly("f_calc", &Class::fcalc)
        .def_property_readonly("scaled_fcalc", &Class::scaled_fcalc)
        .def_property_readonly("usage_flags", &Class::usage)
        .def_property_readonly("base_2fofc", &Class::base_2fofc)
        .def_property_readonly("base_fofc", &Class::base_fofc)
        .def_property_readonly("weights", &Class::weights)
        .def_property_readonly("rwork", &Class::rwork)
        .def_property_readonly("rfree", &Class::rfree)
        .def_property_readonly("weighted_rwork", &Class::weighted_rwork)
        .def_property_readonly("weighted_rfree", &Class::weighted_rfree)
        .def_property_readonly("bulk_solvent_frac", &Class::bulk_frac)
        .def_property_readonly("bulk_solvent_scale", &Class::bulk_scale)
        .def_property_readonly("scaling_function", &Class::scaling_function)
        .def("generate_fcalc", &Class::generate_fcalc)
        .def("generate_base_map_coeffs", &Class::generate_base_map_coeffs)
        .def("add_xmap", &Class::add_xmap,
            py::arg("name"), py::arg("b_sharp"),
            py::arg("is_difference_map")=true,
            py::arg("exclude_missing_reflections")=false,
            py::arg("exclude_free_reflections")=true,
            py::arg("fill_with_fcalc")=true)
        .def("recalculate_map", (void (Class::*)(const std::string&, size_t)) &Class::recalculate_map,
            py::arg("name"), py::arg("num_threads")=1)
        // Get a reference to the managed xmap of a given name
        .def("get_xmap_ref", &Class::get_xmap, py::return_value_policy::reference_internal)
        .def("get_xmap_copy", &Class::get_xmap)
        .def("delete_xmap", &Class::delete_xmap)
        ;
} // declare_xal_mgr

void declare_xtal_thread_mgr(py::module& m)
{
    using Class=cx::Xtal_thread_mgr;
    py::class_<Class>(m, "Xtal_thread_mgr")
        .def(py::init<const HKL_info&, const HKL_data<Flag>&, const Grid_sampling&,
            const HKL_data<F_sigF<ftype32>>&, const size_t>(),
            py::arg("hkl_info"), py::arg("free_flags"), py::arg("grid_sampling"),
            py::arg("f_obs"), py::arg("num_threads") = 1)
        .def("init", &Class::init)
        .def_property("num_threads", &Class::num_threads, &Class::set_num_threads)
        .def_property_readonly("thread_running", &Class::thread_running)
        .def("ready", &Class::ready)
        .def("recalculate_all_maps", &Class::recalculate_all)
        .def("apply_new_maps", &Class::apply_new_maps)
        .def_property_readonly("free_flag", &Class::freeflag)
        .def_property_readonly("flags", &Class::flags)
        .def_property_readonly("usage_flags", &Class::usage_flags)
        .def_property_readonly("rwork", &Class::rwork)
        .def_property_readonly("rfree", &Class::rfree)
        .def_property_readonly("weighted_rwork", &Class::weighted_rwork)
        .def_property_readonly("weighted_rfree", &Class::weighted_rfree)
        .def_property_readonly("bulk_frac", &Class::bulk_frac)
        .def_property_readonly("bulk_scale", &Class::bulk_scale)
        .def_property_readonly("f_obs", &Class::fobs)
        .def_property_readonly("f_calc", &Class::fcalc)
        .def_property_readonly("scaled_fcalc", &Class::scaled_fcalc)
        .def_property_readonly("base_fofc", &Class::base_fofc)
        .def_property_readonly("base_2fofc", &Class::base_2fofc)
        .def_property_readonly("weights", &Class::weights)
        .def_property_readonly("num_maps", &Class::n_maps)
        .def_property_readonly("map_names", &Class::map_names)
        .def_property_readonly("scaling_function", &Class::scaling_function)
        .def("add_xmap", &Class::add_xmap,
            py::arg("name"), py::arg("b_sharp"),
            py::arg("is_difference_map")=false,
            py::arg("exclude_missing_reflections")=false,
            py::arg("exclude_free_reflections")=true,
            py::arg("fill_with_fcalc")=true)
        .def("delete_xmap", &Class::delete_xmap)
        .def("_xmap_details", &Class::map_details)
        .def("get_xmap_ref", &Class::get_xmap, py::return_value_policy::reference_internal)
        .def("get_xmap_copy", &Class::get_xmap)
        .def("get_map_stats", &Class::get_map_stats)
        ;
}



void init_xtal_mgr(py::module& m)
{
    declare_xmap_container(m);
    declare_xtal_mgr(m);
    declare_xtal_thread_mgr(m);
} // init_xtal_mgr
