#include <pybind11/pybind11.h>

#include "../type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

namespace py=pybind11;
using namespace clipper;
using namespace clipper::datatypes;

template <class T>
void declare_edcalc_base(py::module& m, const char* dtype)
{
    using Class=EDcalc_base<T>;
    auto pyclass_name = std::string("_EDcalc_base_") + dtype;
    py::class_<Class, std::unique_ptr<Class, py::nodelete>>(m, pyclass_name.c_str())
        .def("__call__",
            []( const Class& self, Xmap<T>& xmap, const Atom_list& atoms)
            { return self(xmap, atoms); })
        .def("__call__",
            [](const Class& self, NXmap<T>& nxmap, const Atom_list& atoms)
            { return self(nxmap, atoms); })
        ;
} // declare_edcalc_base

template<class Class, class T>
void declare_edcalc_specialization(py::module& m, const char* pyclass_name)
{
    py::class_<Class, EDcalc_base<T>>(m, pyclass_name)
        .def(py::init<const ftype>())
        ;
} // declare_edcalc_specialization

void init_edcalc(py::module& m)
{
    declare_edcalc_base<ftype32>(m, "float");
    declare_edcalc_specialization<EDcalc_mask<ftype32>, ftype32>(m, "EDcalc_mask_float");
    declare_edcalc_specialization<EDcalc_iso<ftype32>, ftype32>(m, "EDcalc_iso_float");
    declare_edcalc_specialization<EDcalc_aniso<ftype32>, ftype32>(m, "EDcalc_aniso_float");


    declare_edcalc_base<ftype64>(m, "double");
    declare_edcalc_specialization<EDcalc_mask<ftype64>, ftype64>(m, "EDcalc_mask_double");
    declare_edcalc_specialization<EDcalc_iso<ftype64>, ftype64>(m, "EDcalc_iso_double");
    declare_edcalc_specialization<EDcalc_aniso<ftype64>, ftype64>(m, "EDcalc_aniso_double");

}
