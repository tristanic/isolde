#include <pybind11/pybind11.h>

#include "../type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

namespace py=pybind11;
using namespace clipper;
using namespace clipper::datatypes;

template <class T>
void declare_sfweight_base(py::module& m, const char* dtype)
{
    using Class=SFweight_base<T>;
    auto pyclass_name = std::string("_SFweight_base_") + dtype;
    py::class_<Class, std::unique_ptr<Class, py::nodelete>>(m, pyclass_name.c_str())
        .def("__call__",
            [](Class& self,
               HKL_data<F_phi<T>>& fb,
               HKL_data<F_phi<T>>& fd,
               HKL_data<Phi_fom<T>>& phiw,
               const HKL_data<F_sigF<T>>& fo,
               const HKL_data<F_phi<T>>& fc,
               const HKL_data<Flag>& usage)
             { return self(fb, fd, phiw, fo, fc, usage); })
         ;
} // declare_sfweight_base

template <class T>
void declare_sfweight_spline(py::module& m, const char* dtype)
{
    using Class = SFweight_spline<T>;
    auto pyclass_name = std::string("SFweight_spline_") + dtype;
    py::class_<Class, SFweight_base<T>> sfweight_spline(m, pyclass_name.c_str());
    sfweight_spline
        .def(py::init<const int, const int, const int>(),
            py::arg("n_reflns")=1000, py::arg("n_params")=20, py::arg("n_phases")=24 )
        .def(py::init<HKL_data<F_phi<T>>&, HKL_data<F_phi<T>>&, HKL_data<Phi_fom<T>>&,
                const HKL_data<F_sigF<T>>&, const HKL_data<F_phi<T>>&, const HKL_data<Flag>&,
                const int, const int>(),
            py::arg("fb"), py::arg("fd"), py::arg("phiw"), py::arg("fo"), py::arg("fc"), py::arg("usage"),
            py::arg("n_reflns")=1000, py::arg("n_params")=20 )
        .def("init", &Class::init, py::arg("n_reflns")=1000, py::arg("n_params")=20, py::arg("n_phases")=24 )
        .def("__call__",
            [](Class& self,
                HKL_data<F_phi<T>>& fb,
                HKL_data<F_phi<T>>& fd,
                HKL_data<Phi_fom<T>>& phiw,
                HKL_data<ABCD<T>>& hl,
                const HKL_data<F_sigF<T>>& fo0,
                const HKL_data<ABCD<T>>& hl0,
                const HKL_data<F_phi<T>>& fc0,
                const HKL_data<Flag>& usage)
            { return self(fb, fd, phiw, hl, fo0, hl0, fc0, usage); })
        .def_property_readonly("params_scale", &Class::params_scale)
        .def_property_readonly("params_error", &Class::params_error)
        .def_property_readonly("log_likelihood_work", &Class::log_likelihood_work)
        .def_property_readonly("log_likelihood_free", &Class::log_likelihood_free)
        .def("targetfn",
            [](const Class& self, const HKL_class cls, const F_sigF<T>& fo0,
                const F_phi<T>& fc0, const ftype& s, const ftype& w)
                -> py::tuple
            {
                auto r = self.targetfn(cls, fo0, fc0, s, w);
                return py::make_tuple(r.r, r.ds, r.dw, r.dss, r.dww, r.dsw);
            })
        .def("targethl",
            [](const Class& self, const HKL_class cls, const F_sigF<T>& fo0,
                const ABCD<T>& hl0, const F_phi<T>& fc0, const ftype& s, const ftype& w)
                -> py::tuple
            {
                auto r = self.targethl(cls, fo0, hl0, fc0, s, w);
                return py::make_tuple(r.r, r.ds, r.dw, r.dss, r.dww, r.dsw);
            })
        ;
}


void init_sfweight(py::module& m)
{
    declare_sfweight_base<ftype32>(m, "float");
    declare_sfweight_spline<ftype32>(m, "float");

    declare_sfweight_base<ftype64>(m, "double");
    declare_sfweight_spline<ftype64>(m, "double");
}
