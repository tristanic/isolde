#include <pybind11/pybind11.h>

#include <clipper_ext/xtal_mgr.h>

namespace py=pybind11;
namespace cx=clipper_cx;
using namespace clipper;
using namespace clipper::datatypes;

void declare_xtal_mgr(py::module& m)
{
    using Class=cx::Xtal_mgr_base;
    py::class_<Class>(m, "Xtal_mgr_base")
        .def(py::init<const HKL_info&, const HKL_data<Flag>&, const Grid_sampling&,
            const HKL_data<F_sigF<ftype32>>&>() )
        .def_property_readonly("free_flag", &Class::freeflag)
        .def_property_readonly("f_calc", &Class::fcalc)
        .def_property_readonly("usage_flags", &Class::usage)
        .def_property_readonly("base_2fofc", &Class::base_2fofc)
        .def_property_readonly("base_fofc", &Class::base_fofc)
        .def_property_readonly("weights", &Class::weights)
        .def("generate_fcalc", &Class::generate_fcalc)
        .def("generate_base_map_coeffs", &Class::generate_base_map_coeffs)
        ;
}

void init_xtal_mgr(py::module& m)
{
    declare_xtal_mgr(m);
} // init_xtal_mgr
