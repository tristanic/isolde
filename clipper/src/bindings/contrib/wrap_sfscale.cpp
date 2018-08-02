#include <pybind11/pybind11.h>

#include "../type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

namespace py=pybind11;
using namespace clipper;
using namespace clipper::datatypes;

template <class T>
void declare_sfscale_base(py::module &m, const char* dtype)
{
    using Class=SFscale_base<T>;
    auto pyclass_name = std::string("_SFscale_base_") + dtype;
    py::class_<Class, std::unique_ptr<Class, py::nodelete>>(m, pyclass_name.c_str())
        .def("__call__",
            [](Class& self, HKL_data<F_sigF<T>>& fo, HKL_data<F_phi<T>>& fc)
            { return self(fo, fc); })
        .def("__call__",
            [](Class& self, HKL_data<F_phi<T>>& fc, HKL_data<F_sigF<T>>& fo)
            { return self(fc, fo); })
        .def("__call__",
            [](Class& self, HKL_data<F_sigF<T>>& fo)
            { return self(fo); })
        ;
} // declare_sfscale_base

template<class Derived, class Base, class D, class T1, class T2, class S>
void add_scale_specialization(py::class_<Derived, Base>& pyclass, const char* basis_type)
{
    auto fn_name = std::string("scale_") + basis_type;
    pyclass.def(fn_name.c_str(),
    [](Derived& self, HKL_data<D>& fo, const ftype resfilter, const int npar_scl)
    {
        return self.template scale<D, T1, T2, S>(fo, resfilter, npar_scl);
    })
    ;
    // (bool (Derived::*)(HKL_data<D>&, const ftype, const int)) &Derived::scale<D, T1, T2, S>);
}

template<class T>
void declare_sfscale_aniso(py::module& m, const char* dtype)
{
    using Class=SFscale_aniso<T>;
    auto pyclass_name = std::string("SFscale_aniso_") + dtype;
    py::class_<Class, SFscale_base<T>> sfscale_aniso(m, pyclass_name.c_str());

    py::enum_<typename Class::TYPE>(sfscale_aniso, "TYPE")
        .value("F", Class::TYPE::F)
        .value("I", Class::TYPE::I)
        .export_values()
        ;
    py::enum_<typename Class::MODE>(sfscale_aniso, "MODE")
        .value("NORMAL", Class::MODE::NORMAL)
        .value("SHARPEN", Class::MODE::SHARPEN)
        .value("UNSHARPEN", Class::MODE::UNSHARPEN)
        .export_values()
        ;

    sfscale_aniso
        .def(py::init<ftype, typename Class::MODE>(),
            py::arg("nsig")=0.0, py::arg("mode")=Class::NORMAL)
        .def("__call__",
            [](Class& self, HKL_data<F_sigF<T>>& fo, const ftype resfilter, const int npar_scl)
            { return self(fo, resfilter, npar_scl); })
        .def("__call__",
            [](Class& self, HKL_data<I_sigI<T>>& io, const ftype resfilter, const int npar_scl)
            { return self(io, resfilter, npar_scl); })
        .def("u_aniso_orth", (const U_aniso_orth& (Class::*)(typename Class::TYPE) const) &Class::u_aniso_orth)
        ;

    add_scale_specialization<Class, SFscale_base<T>, F_sigF<T>,
                          TargetFn_scaleF1F2<F_sigF<T>, F_sigF<T>>,
                          TargetFn_scaleLogF1F2<F_sigF<T>, F_sigF<T>>,
                          BasisFn_binner>
                          (sfscale_aniso, "binned");
    add_scale_specialization<Class, SFscale_base<T>, F_sigF<T>,
                          TargetFn_scaleF1F2<F_sigF<T>, F_sigF<T>>,
                          TargetFn_scaleLogF1F2<F_sigF<T>, F_sigF<T>>,
                          BasisFn_linear>
                          (sfscale_aniso, "linear");
    add_scale_specialization<Class, SFscale_base<T>, F_sigF<T>,
                          TargetFn_scaleF1F2<F_sigF<T>, F_sigF<T>>,
                          TargetFn_scaleLogF1F2<F_sigF<T>, F_sigF<T>>,
                          BasisFn_spline>
                          (sfscale_aniso, "spline");

    add_scale_specialization<Class, SFscale_base<T>, I_sigI<T>,
                        TargetFn_scaleI1I2<I_sigI<T>, I_sigI<T>>,
                        TargetFn_scaleLogI1I2<I_sigI<T>, I_sigI<T>>,
                        BasisFn_binner>
                        (sfscale_aniso, "binned");
    add_scale_specialization<Class, SFscale_base<T>, I_sigI<T>,
                        TargetFn_scaleI1I2<I_sigI<T>, I_sigI<T>>,
                        TargetFn_scaleLogI1I2<I_sigI<T>, I_sigI<T>>,
                        BasisFn_linear>
                        (sfscale_aniso, "linear");
    add_scale_specialization<Class, SFscale_base<T>, I_sigI<T>,
                        TargetFn_scaleI1I2<I_sigI<T>, I_sigI<T>>,
                        TargetFn_scaleLogI1I2<I_sigI<T>, I_sigI<T>>,
                        BasisFn_spline>
                        (sfscale_aniso, "spline");

}


void init_sfscale(py::module& m)
{
    declare_sfscale_base<ftype32>(m, "float");
    declare_sfscale_aniso<ftype32>(m, "float");


    declare_sfscale_base<ftype64>(m, "double");
    declare_sfscale_aniso<ftype64>(m, "double");
} //init_sfscale
