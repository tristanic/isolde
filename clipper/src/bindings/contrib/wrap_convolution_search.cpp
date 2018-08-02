#include <pybind11/pybind11.h>

#include "../type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

namespace py=pybind11;
using namespace clipper;
using namespace clipper::datatypes;

template <class T>
void declare_convolution_search_base(py::module& m, const char* dtype)
{
    using Class=Convolution_search_base<T>;
    auto pyclass_name = std::string("_Convolution_search_base_") + dtype;
    py::class_<Class, std::unique_ptr<Class, py::nodelete>>(m, pyclass_name.c_str())
        .def("__call__",
            [](const Class& self, Xmap<T>& result,
               const NXmap<T>& srchval, const NX_operator& nxop)
            { return self(result, srchval, nxop); })
        ;
} // declare_convolution_search_base


template <class T>
void declare_convolution_search_fft(py::module& m, const char* dtype)
{
    using Class=Convolution_search_fft<T>;
    auto pyclass_name = std::string("Convolution_search_fft_") + dtype;
    py::class_<Class, Convolution_search_base<T>>(m, pyclass_name.c_str())
        .def(py::init<const Xmap<T>&>())
        .def(py::init<Xmap<T>&, const NXmap<T>&, const Xmap<T>&, const NX_operator&>())
        ;
} // declare_convolution_search_fft


void init_convolution_search(py::module& m)
{
    declare_convolution_search_base<ftype32>(m, "float");
    declare_convolution_search_fft<ftype32>(m, "float");

    declare_convolution_search_base<ftype64>(m, "double");
    declare_convolution_search_fft<ftype64>(m, "double");
} // init_convolution_search
