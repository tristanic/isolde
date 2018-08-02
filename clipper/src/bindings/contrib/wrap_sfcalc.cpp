#include <pybind11/pybind11.h>

#include "../type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

namespace py=pybind11;
using namespace clipper;
using namespace clipper::datatypes;

template <class T>
void declare_sfcalc_base(py::module& m, const char* dtype)
{
    using Class=SFcalc_base<T>;
    auto pyclass_name = std::string("_SFcalc_base_") + dtype;
    py::class_<Class, std::unique_ptr<Class, py::nodelete>>(m, pyclass_name.c_str())
        .def("__call__",
            [](const Class& self,
               HKL_data<F_phi<T>>& fphidata,
               const Atom_list& atoms )
            { return self(fphidata, atoms); })
        ;
} // declare_sfcalc_base

template <class T>
void declare_sfcalc_iso_sum(py::module& m, const char* dtype)
{
    using Class=SFcalc_iso_sum<T>;
    auto pyclass_name = std::string("SFcalc_iso_sum_") + dtype;
    py::class_<Class, SFcalc_base<T>>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<HKL_data<F_phi<T>>&, const Atom_list&>())
        ;
} // declare_sfcalc_iso_sum

template <class T>
void declare_sfcalc_aniso_sum(py::module& m, const char* dtype)
{
    using Class=SFcalc_aniso_sum<T>;
    auto pyclass_name = std::string("SFcalc_aniso_sum_") + dtype;
    py::class_<Class, SFcalc_base<T>>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<HKL_data<F_phi<T>>&, const Atom_list&>())
        ;
} // declare_sfcalc_aniso_sum

template <class T>
void declare_sfcalc_iso_fft(py::module& m, const char* dtype)
{
    using Class=SFcalc_iso_fft<T>;
    auto pyclass_name = std::string("SFcalc_iso_fft_") + dtype;
    py::class_<Class, SFcalc_base<T>>(m, pyclass_name.c_str())
        .def(py::init<const ftype, const ftype, const ftype>(),
            py::arg("radius")=2.5, py::arg("rate")=1.5, py::arg("uadd")=0.0)
        .def(py::init<HKL_data<F_phi<T>>&, const Atom_list&, const ftype, const ftype, const ftype>(),
            py::arg("fphidata_out"), py::arg("atoms"), py::arg("radius")=2.5, py::arg("rate")=1.5, py::arg("uadd")=0.0)
        ;
} // declare_sfcalc_iso_fft

template <class T>
void declare_sfcalc_aniso_fft(py::module& m, const char* dtype)
{
    using Class=SFcalc_aniso_fft<T>;
    auto pyclass_name = std::string("SFcalc_aniso_fft_") + dtype;
    py::class_<Class, SFcalc_base<T>>(m, pyclass_name.c_str())
        .def(py::init<const ftype, const ftype, const ftype>(),
            py::arg("radius")=2.5, py::arg("rate")=1.5, py::arg("uadd")=0.0)
        .def(py::init<HKL_data<F_phi<T>>&, const Atom_list&, const ftype, const ftype, const ftype>(),
            py::arg("fphidata_out"), py::arg("atoms"), py::arg("radius")=2.5, py::arg("rate")=1.5, py::arg("uadd")=0.0)
        ;
} // declare_sfcalc_aniso_fft


void init_sfcalc(py::module& m)
{
    declare_sfcalc_base<ftype32>(m, "float");
    declare_sfcalc_iso_sum<ftype32>(m, "float");
    declare_sfcalc_aniso_sum<ftype32>(m, "float");
    declare_sfcalc_iso_fft<ftype32>(m, "float");
    declare_sfcalc_aniso_fft<ftype32>(m, "float");

    declare_sfcalc_base<ftype64>(m, "double");
    declare_sfcalc_iso_sum<ftype64>(m, "double");
    declare_sfcalc_aniso_sum<ftype64>(m, "double");
    declare_sfcalc_iso_fft<ftype64>(m, "double");
    declare_sfcalc_aniso_fft<ftype64>(m, "double");

} // init_sfcalc
