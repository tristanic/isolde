#include <pybind11/pybind11.h>

#include "../type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

namespace py=pybind11;
using namespace clipper;
using namespace clipper::datatypes;






template <class T>
void declare_skeleton_base(py::module &m, const char* dtype)
{
    using Class=Skeleton_base<int, T>;
    auto pyclass_name = std::string("_Skeleton_base_") + dtype;
    py::class_<Class, std::unique_ptr<Class, py::nodelete>>(m, pyclass_name.c_str())
        .def("__call__",
            [](const Class& self, Xmap<int>& xskl, const Xmap<T>& xmap)
            { return self(xskl, xmap); })
        ;
} // declare_skeleton_base


void init_sfcalc_base(py::module& m)
{
    declare_skeleton_base<ftype32>(m, "float");

    declare_skeleton_base<ftype64>(m, "double");

} // init_sfcalc
