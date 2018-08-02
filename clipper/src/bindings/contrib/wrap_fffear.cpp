#include <pybind11/pybind11.h>

#include "../type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

namespace py=pybind11;
using namespace clipper;
using namespace clipper::datatypes;

template <class T>
void declare_fffear_base(py::module &m, const char* dtype)
{
    using Class=FFFear_base<T>;
    auto pyclass_name = std::string("_FFFear_base_") + dtype;
    py::class_<Class, std::unique_ptr<Class, py::nodelete>>(m, pyclass_name.c_str())
        .def("__call__",
            [](const Class& self, Xmap<T>& result,
               const NXmap<T>& srchval, const NXmap<T>& srchwgt, const NX_operator& nxop)
            { return self(result, srchval, srchwgt, nxop); })
        ;
} // declare_fffear_base

template <class T>
void declare_fffear_fft(py::module& m, const char* dtype)
{
    using Class=FFFear_fft<T>;
    auto pyclass_name = std::string("FFFear_fft") + dtype;
    py::class_<Class, FFFear_base<T>>(m, pyclass_name.c_str())
        .def(py::init<const Xmap<T>&>())
        .def(py::init<Xmap<T>&, const NXmap<T>&, const NXmap<T>&, const Xmap<T>&, const NX_operator&>())
        .def("set_fft_type", &Class::set_fft_type)
        .def("set_resolution", &Class::set_resolution)
        .def("__call__",
            [](const Class& self, Xmap<T>& result,
               const NXmap<T>& srchval, const NXmap<T>& srchwgt, const RTop_orth& rtop)
            { return self(result, srchval, srchwgt, rtop); })
        ;
}


void init_fffear(py::module& m)
{
    declare_fffear_base<ftype32>(m, "float");
    declare_fffear_fft<ftype32>(m, "float");


    declare_fffear_base<ftype64>(m, "double");
    declare_fffear_fft<ftype64>(m, "double");

} // init_fffear
