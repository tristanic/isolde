#include <pybind11/pybind11.h>

#include "../type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

namespace py=pybind11;
using namespace clipper;
using namespace clipper::datatypes;

template <class T>
void declare_originmatch_base(py::module& m, const char* dtype)
{
    using Class=OriginMatch_base<T>;
    auto pyclass_name = std::string("_OriginMatch_base_") + dtype;
    py::class_<Class, std::unique_ptr<Class, py::nodelete>>(m, pyclass_name.c_str())
        .def("__call__",
            [](const Class& self, bool& invert, Coord_frac& shift,
               const HKL_data<F_phi<T>>& fphi1, const HKL_data<F_phi<T>>& fphi2)
            { return self(invert, shift, fphi1, fphi2); })
        ;
} // declare_originmatch_base

template <class T>
void declare_originmatch(py::module& m, const char* dtype)
{
    using Class=OriginMatch<T>;
    auto pyclass_name = std::string("OriginMatch_") + dtype;
    py::class_<Class, OriginMatch_base<T>>(m, pyclass_name.c_str())
        .def(py::init<const ftype>(), py::arg("resol_limit") = 0.1)
        .def(py::init<bool&, Coord_frac&, const HKL_data<F_phi<T>>&, const HKL_data<F_phi<T>>&, const ftype>(),
            py::arg("invert"), py::arg("shift"), py::arg("fphi1"), py::arg("fphi2"), py::arg("resol_limit")=0.1)
        ;
} // declare_originmatch


void init_originmatch(py::module& m)
{
    declare_originmatch_base<ftype32>(m, "float");
    declare_originmatch<ftype32>(m, "float");

    declare_originmatch_base<ftype64>(m, "double");
    declare_originmatch<ftype64>(m, "double");

}
