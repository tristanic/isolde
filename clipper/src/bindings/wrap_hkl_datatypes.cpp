#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

#include "type_conversions.h"
#include <clipper/clipper.h>
#include "numpy_helper.h"

namespace py=pybind11;
using namespace clipper;
using namespace clipper::datatypes;

#define GETSET(CLASS, FUNCNAME) \
[] (const CLASS& self) { return self.FUNCNAME(); },\
[] (CLASS& self, ftype& val) { self.FUNCNAME() = val; }



template<class C>
void catch_null(const C& c)
{
    if (c.is_null()) throw std::length_error("Array is not initialised!");
}

// template<class C, typename F, class dtype, class... Sources,  typename... Args>
// void safe_compute(C& target, const F& func, const Sources&... sources, const Args&... args)
// {
//     catch_null(target);
//     target.compute(sources, func(args));
// }


template<class Derived>
void declare_base_methods(py::class_<Derived> pyclass)
{
pyclass
    .def("set_null", &Derived::set_null)
    .def_property_readonly("type", [](const Derived& self) { return self.type().c_str(); })
    .def("friedel", &Derived::friedel)
    .def("shift_phase", &Derived::shift_phase)
    .def_property_readonly("missing", &Derived::missing)
    //.def_property_readonly("data_size", &Derived::data_size)
    .def("__len__", &Derived::data_size)
    .def_property_readonly("data_names", [](const Derived& self) { return self.data_names().c_str(); })
    // To/from numpy
    .def_property("data",
        [](const Derived& self) -> py::array_t<xtype>
            { return hkl_data_export_numpy<Derived, xtype>(self, self.data_size()); },
        [](Derived& self, py::array_t<xtype> vals)
            { hkl_data_import_numpy<Derived, xtype>(self, self.data_size(), vals); }
    )
    // Python methods common to all
    .def("copy", [](const Derived& self) -> Derived { return self; } )
    ;
} // declare_base_methods

template <class T>
void declare_i_sigi(py::module& m, const char* dtype)
{
    using Class=I_sigI<T>;
    auto pyclass_name = std::string("I_sigI_") + dtype;
    py::class_<Class> isigi(m, pyclass_name.c_str());
    isigi
        .def(py::init<>())
        .def(py::init<const T&, const T&>())
        .def("scale", &Class::scale)
        .def_property("i", GETSET(Class, I))
        .def_property("sigi", GETSET(Class, sigI))
        .def_property_readonly("i_pl", &Class::I_pl)
        .def_property_readonly("sigi_pl", &Class::sigI_pl)
        .def_property_readonly("i_mi", &Class::I_mi)
        .def_property_readonly("sigi_mi", &Class::sigI_mi)
        .def_property_readonly("cov", &Class::cov)
        ;
    declare_base_methods<Class>(isigi);
} // declare_i_sigi

template<class T>
void declare_i_sigi_ano(py::module& m, const char* dtype)
{
    using Class=I_sigI_ano<T>;
    auto pyclass_name = std::string("I_sigI_ano_") + dtype;
    py::class_<Class> isigi_ano(m, pyclass_name.c_str());
    isigi_ano
        .def(py::init<>())
        .def("scale", &Class::scale)
        .def_property("i_pl", GETSET(Class, I_pl))
        .def_property("sigi_pl", GETSET(Class, sigI_pl))
        .def_property("i_mi", GETSET(Class, I_mi))
        .def_property("sigi_mi", GETSET(Class, sigI_mi))
        .def_property("cov", GETSET(Class, cov))
        .def_property_readonly("i", &Class::I)
        .def_property_readonly("sigi", &Class::sigI)
        ;
    declare_base_methods<Class>(isigi_ano);
} // declare_i_sigi_ano

template <class T>
void declare_f_sigf(py::module& m, const char* dtype)
{
    using Class=F_sigF<T>;
    auto pyclass_name = std::string("F_sigF_") + dtype;
    py::class_<Class> fsigf(m, pyclass_name.c_str());
    fsigf
        .def(py::init<>())
        .def(py::init<const T&, const T&>())
        .def("scale", &Class::scale)
        .def_property("f", GETSET(Class, f))
        .def_property("sigf", GETSET(Class, sigf))
        .def_property_readonly("f_pl", &Class::f_pl)
        .def_property_readonly("sigf_pl", &Class::sigf_pl)
        .def_property_readonly("f_mi", &Class::f_mi)
        .def_property_readonly("sigf_mi", &Class::sigf_mi)
        .def_property_readonly("cov", &Class::cov)
        ;
    declare_base_methods<Class>(fsigf);
} // declare_f_sigf

template <class T>
void declare_f_sigf_ano(py::module& m, const char* dtype)
{
    using Class=F_sigF_ano<T>;
    auto pyclass_name=std::string("F_sigF_ano_") + dtype;
    py::class_<Class> fsigf_ano(m, pyclass_name.c_str());
    fsigf_ano
        .def(py::init<>())
        .def("scale", &Class::scale)
        .def_property("f_pl", GETSET(Class, f_pl))
        .def_property("sigf_pl", GETSET(Class, sigf_pl))
        .def_property("f_mi", GETSET(Class, f_mi))
        .def_property("sigf_mi", GETSET(Class, sigf_mi))
        .def_property("cov", GETSET(Class, cov))
        .def_property_readonly("f", &Class::f)
        .def_property_readonly("sigf", &Class::sigf)
        ;
    declare_base_methods<Class>(fsigf_ano);
} // declare_f_sigf_ano

template <class T>
void declare_e_sige(py::module& m, const char* dtype)
{
    using Class=E_sigE<T>;
    auto pyclass_name=std::string("E_sigE_") + dtype;
    py::class_<Class> esige(m, pyclass_name.c_str());
    esige
        .def(py::init<>())
        .def(py::init<const T&, const T&>())
        .def_property("e", GETSET(Class, E))
        .def_property("sige", GETSET(Class, sigE))
        .def_property_readonly("e_pl", &Class::E_pl)
        .def_property_readonly("sige_pl", &Class::sigE_pl)
        .def_property_readonly("e_mi", &Class::E_mi)
        .def_property_readonly("sige_mi", &Class::sigE_mi)
        .def_property_readonly("cov", &Class::cov)
        ;
    declare_base_methods<Class>(esige);
} // declare_e_sige

// E_sigE_ano template is not currently instantiated in the Clipper core
// template<class T>
// void declare_e_sige_ano(py::module& m, const char* dtype)
// {
//     using Class=E_sigE_ano<T>;
//     auto pyclass_name=std::string("E_sigE_ano_") + dtype;
//     py::class_<Class> esige_ano(m, pyclass_name.c_str());
//     esige_ano
//         .def(py::init<>())
//         .def("scale", &Class::scale)
//         .def_property("e_pl", GETSET(Class, E_pl))
//         .def_property("sige_pl", GETSET(Class, sigE_pl))
//         .def_property("e_mi", GETSET(Class, E_mi))
//         .def_property("sige_mi", GETSET(Class, sigE_mi))
//         .def_property("cov", GETSET(Class, cov))
//         .def_property_readonly("e", &Class::E)
//         .def_property_readonly("sige", &Class::sigE)
//         ;
//     declare_base_methods<Class>(esige_ano);
// } // declare_e_sige_ano

template<class T>
void declare_f_phi(py::module &m, const char* dtype)
{
    using Class=F_phi<T>;
    auto pyclass_name=std::string("F_phi_") + dtype;
    py::class_<Class> fphi(m, pyclass_name.c_str());
    fphi
        .def(py::init<>())
        .def(py::init<const T&, const T&>())
        .def(py::init<const std::complex<T>>())
        .def("scale", &Class::scale)
        .def_property("f", GETSET(Class, f))
        .def_property("phi", GETSET(Class, phi))
        .def_property_readonly("a", &Class::a)
        .def_property_readonly("b", &Class::b)
        .def_property_readonly("complex", [](const Class& self) { return std::complex<T>(self); })
        .def("resolve", &Class::resolve)
        .def("norm", &Class::norm)
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(-py::self)
        ;
    declare_base_methods<Class>(fphi);
} //declare_f_phi

template <class T>
void declare_phi_fom(py::module& m, const char* dtype)
{
    using Class=Phi_fom<T>;
    auto pyclass_name=std::string("Phi_fom_") + dtype;
    py::class_<Class> phifom(m, pyclass_name.c_str());
    phifom
        .def(py::init<>())
        .def(py::init<const T&, const T&>())
        .def_property("phi", GETSET(Class, phi))
        .def_property("fom", GETSET(Class, fom))
        ;
    declare_base_methods<Class>(phifom);
} // declare_phi_fom

template<class T>
void declare_abcd(py::module& m, const char* dtype)
{
    using Class=ABCD<T>;
    auto pyclass_name=std::string("ABCD_") + dtype;
    py::class_<Class> abcd(m, pyclass_name.c_str());
    abcd
        .def(py::init<>())
        .def(py::init<const T&, const T&, const T&, const T&>())
        .def_property("a", GETSET(Class, a))
        .def_property("b", GETSET(Class, b))
        .def_property("c", GETSET(Class, c))
        .def_property("d", GETSET(Class, d))
        .def(py::self + py::self)
        ;
    declare_base_methods<Class>(abcd);
} // declare_abcd

void declare_flag(py::module& m)
{
    py::class_<Flag> flag(m, "Flag");
    flag
        .def(py::init<>())
        .def(py::init<const int&>())
        .def_property("flag", GETSET(Flag, flag))
        ;
    declare_base_methods<Flag>(flag);
} // declare_flag

void declare_flag_bool(py::module& m)
{
    py::class_<Flag_bool> flag_bool(m, "Flag_bool");
    flag_bool
        .def(py::init<>())
        .def_property("flag", GETSET(Flag_bool, flag))
        ;
    declare_base_methods<Flag_bool>(flag_bool);
} // declare_flag_bool

template<class T>
void declare_d_sigd(py::module& m, const char* dtype)
{
    using Class=D_sigD<T>;
    auto pyclass_name=std::string("D_sigD_") + dtype;
    py::class_<Class> dsigd(m, pyclass_name.c_str(),
    "Deprecated anomalous difference class for backward compatibility only. Do not use.");
    dsigd
        .def(py::init<>())
        .def(py::init<const T&, const T&>())
        .def("scale", &Class::scale)
        .def_property("d", GETSET(Class, d))
        .def_property("sigd", GETSET(Class, sigd))
        ;
    declare_base_methods<Class>(dsigd);
} // declare_flag_bool


// Actual wrapper initialisation
void init_hkl_datatypes(py::module& m, py::module& m32, py::module& m64)
{
    // Non-floating-point datatypes go in the main module
    declare_flag(m);
    declare_flag_bool(m);
    {
        using namespace clipper::data32;
        // 32-bit floating point datatypes go in the data32 module
        declare_i_sigi<ftype32>(m32, "float");
        declare_i_sigi_ano<ftype32>(m32, "float");
        declare_f_sigf<ftype32>(m32, "float");
        declare_f_sigf_ano<ftype32>(m32, "float");
        declare_e_sige<ftype32>(m32, "float");
        declare_f_phi<ftype32>(m32, "float");
        declare_phi_fom<ftype32>(m32, "float");
        declare_abcd<ftype32>(m32, "float");
        declare_d_sigd<ftype32>(m32, "float");
    }
    {
        using namespace clipper::data64;
        // 64-bit floating point datatypes go in the data64 module
        declare_i_sigi<ftype64>(m64, "double");
        declare_i_sigi_ano<ftype64>(m64, "double");
        declare_f_sigf<ftype64>(m64, "double");
        declare_f_sigf_ano<ftype64>(m64, "double");
        declare_e_sige<ftype64>(m64, "double");
        declare_f_phi<ftype64>(m64, "double");
        declare_phi_fom<ftype64>(m64, "double");
        declare_abcd<ftype64>(m64, "double");
        declare_d_sigd<ftype64>(m64, "double");
    }

}
