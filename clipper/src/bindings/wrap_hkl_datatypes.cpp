#pragma once

#include <pybind11/pybind11.h>

#include <clipper/clipper.h>

namespace py=pybind11;
using namespace clipper;

template <class C>
void declare_HKL_data(py::module &m, const std::string &class_str)
{
    using Class=HKL_data<C>;
    std::string pyclass_name = std::string("HKL_data_") + class_str;
    py::class_<Class>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const HKL_info&>())
        .def(py::init<const HKL_info&, const Cell&>())
        .def(py::init<const Spacegroup&, const Cell&, const HKL_sampling&>())
        .def(py::init<const HKL_data_base&>())
        .def("init", (void (Class::*)(const HKL_info&, const Cell&)) &Class::init)
        .def("init", (void (Class::*)(const Spacegroup&, const Cell&, const HKL_sampling&)) &Class::init)
        .def("init", (void (Class::*)(const HKL_data_base&)) &Class::init)
        .def("update", &Class::update)
        .def("type", &Class::type)
        .def("missing", &Class::missing)
        .def("set_null", &Class::set_null)
        .def("data_size", &Class::data_size)
        .def("data_names", &Class::data_names)
        .def("data_export", &Class::data_export)
        .def("data_import", &Class::data_import)
        .def("mask", &Class::mask)
        .def("__getitem__", [](const Class& self, const HKL_info::HKL_reference_index& i){return self[i];}, py::is_operator())
        .def("__getitem__", [](const Class& self, const HKL_info::HKL_reference_coord& ih) {
            C data;
            if (self.get_data(ih, data))
                return data;
            throw std::out_of_range("No data equivalent to that HKL!");
        }, py::is_operator())
        .def("__setitem__", [](Class& self, const HKL_info::HKL_reference_coord& ih, const C& data) {
            if (!self.set_data(ih, data))
                throw std::out_of_range("No equivalent HKL has been indexed for this dataset!");
        }, py::is_operator())
        .def("__getitem__", [](const Class& self, const int& index) {return self[index];}, py::is_operator())
        .def("__getitem__", [](const Class& self, const HKL& hkl) {
            C data;
            if (self.get_data(hkl, data))
                return data;
            throw std::out_of_range("No data equivalent to that HKL!");
        }, py::is_operator())
        .def("__setitem__", [](Class& self, const HKL& hkl, const C& data) {
            if (!self.set_data(hkl, data))
                throw std::out_of_range("No equivalent HKL has been indexed for this dataset!");
        })
        .def("copy_from", [](Class& self, const Class& other) { self=other; })
        .def("set_all_values_to", [](Class& self, const C& value) { self=value; })
        /*
        Since the base class HKL_data_base has a protected virtual destructor,
        exposing it to Python via PyBind11 requires it (and all derived classes)
        to be wrapped with std::unique_ptr<Class, py::nodelete>, which would
        require the python side code to explicitly handle object deletion. That
        would be a serious pain, so instead we'll hide the base class from
        Python entirely and expose all base class functions using lambdas.
        */
        .def("is_null", [](HKL_data_base& b) { return b.is_null(); })
        .def("base_hkl_info", [](HKL_data_base& b) { return b.base_hkl_info(); })
        .def("base_cell", [](HKL_data_base& b) { return b.base_cell(); })
        .def("spacegroup", [](HKL_data_base& b) { return b.spacegroup(); })
        .def("cell", [](HKL_data_base& b) { return b.cell(); })
        .def("resolution", [](HKL_data_base& b) { return b.resolution(); })
        .def("hkl_sampling", [](HKL_data_base& b) { return b.hkl_sampling(); })
        .def("hkl_info", [](HKL_data_base& b) { return b.hkl_info(); })
        .def("invresolsq", [](HKL_data_base& b, const int& index) { return b.invresolsq(index); })
        .def("invresolsq_range", [](HKL_data_base& b) { return b.invresolsq_range(); })
        .def("num_obs", [](HKL_data_base& b) { return b.num_obs(); })
        .def("first", [](HKL_data_base& b) { return b.first(); })
        .def("first_data", [](HKL_data_base& b) { return b.first_data(); })
        .def("next_data", [](HKL_data_base& b, HKL_info::HKL_reference_index& ih) { return b.next_data(ih); });
}

template <class T>
void declare_CHKL_data(py::module &m, const std::string &class_str)
{
    using Class=CHKL_data<T>;
    std::string pyclass_name = std::string("CHKL_data_" + class_str);
    py::class_<Class, Container, HKL_data<T> >(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<Container&, const String>())
        .def("init", (void (Class::*)(const HKL_info&, const Cell&)) &Class::init)
        .def("init", (void (Class::*)(const Spacegroup&, const Cell&, const HKL_sampling&)) &Class::init)
        .def("update", &Class::update)
        .def("copy_from", [](Class& self, const HKL_data<T>& other) { self = other; })
        .def("set_all_values_to", [](Class& self, const T& value) { self = value; });
}

void init_hkl_datatypes(py::module& m, py::module& m32, py::module& m64)
{
    // Non-floating-point datatypes go in the main module
    declare_HKL_data<clipper::datatypes::Flag>(m, "Flag");
    declare_HKL_data<clipper::datatypes::Flag_bool>(m, "Flag_bool");
    declare_CHKL_data<clipper::datatypes::Flag>(m, "Flag");
    declare_CHKL_data<clipper::datatypes::Flag_bool>(m, "Flag_bool");

    // 32-bit floating point datatypes go in the data32 module
    declare_HKL_data<clipper::data32::I_sigI>(m32, "I_sigI_float");
    declare_HKL_data<clipper::data32::I_sigI_ano>(m32, "I_sigI_anom");
    declare_HKL_data<clipper::data32::F_sigF>(m32, "F_sigF");
    declare_HKL_data<clipper::data32::F_sigF_ano>(m32, "F_sigF_anom");
    declare_HKL_data<clipper::data32::E_sigE>(m32, "E_sigE");
    declare_HKL_data<clipper::data32::F_phi>(m32, "F_phi");
    declare_HKL_data<clipper::data32::Phi_fom>(m32, "Phi_phom");
    declare_HKL_data<clipper::data32::ABCD>(m32, "ABCD");
    declare_HKL_data<clipper::data32::D_sigD>(m32, "D_sigD");

    declare_CHKL_data<clipper::data32::I_sigI>(m32, "I_sigI_float");
    declare_CHKL_data<clipper::data32::I_sigI_ano>(m32, "I_sigI_anom");
    declare_CHKL_data<clipper::data32::F_sigF>(m32, "F_sigF");
    declare_CHKL_data<clipper::data32::F_sigF_ano>(m32, "F_sigF_anom");
    declare_CHKL_data<clipper::data32::E_sigE>(m32, "E_sigE");
    declare_CHKL_data<clipper::data32::F_phi>(m32, "F_phi");
    declare_CHKL_data<clipper::data32::Phi_fom>(m32, "Phi_phom");
    declare_CHKL_data<clipper::data32::ABCD>(m32, "ABCD");
    declare_CHKL_data<clipper::data32::D_sigD>(m32, "D_sigD");

    // 64-bit floating point datatypes go in the data64 module
    declare_HKL_data<clipper::data64::I_sigI>(m64, "I_sigI");
    declare_HKL_data<clipper::data64::I_sigI_ano>(m64, "I_sigI_anom");
    declare_HKL_data<clipper::data64::F_sigF>(m64, "F_sigF");
    declare_HKL_data<clipper::data64::F_sigF_ano>(m64, "F_sigF_anom");
    declare_HKL_data<clipper::data64::E_sigE>(m64, "E_sigE");
    declare_HKL_data<clipper::data64::F_phi>(m64, "F_phi");
    declare_HKL_data<clipper::data64::Phi_fom>(m64, "Phi_phom");
    declare_HKL_data<clipper::data64::ABCD>(m64, "ABCD");
    declare_HKL_data<clipper::data64::D_sigD>(m64, "D_sigD");

    declare_CHKL_data<clipper::data64::I_sigI>(m64, "I_sigI_float");
    declare_CHKL_data<clipper::data64::I_sigI_ano>(m64, "I_sigI_anom");
    declare_CHKL_data<clipper::data64::F_sigF>(m64, "F_sigF");
    declare_CHKL_data<clipper::data64::F_sigF_ano>(m64, "F_sigF_anom");
    declare_CHKL_data<clipper::data64::E_sigE>(m64, "E_sigE");
    declare_CHKL_data<clipper::data64::F_phi>(m64, "F_phi");
    declare_CHKL_data<clipper::data64::Phi_fom>(m64, "Phi_phom");
    declare_CHKL_data<clipper::data64::ABCD>(m64, "ABCD");
    declare_CHKL_data<clipper::data64::D_sigD>(m64, "D_sigD");

}
