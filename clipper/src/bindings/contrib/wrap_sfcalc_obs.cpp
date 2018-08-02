#include <pybind11/pybind11.h>

#include "../type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

namespace py=pybind11;
using namespace clipper;
using namespace clipper::datatypes;

template <class T>
void declare_sfcalc_obs_base(py::module& m, const char* dtype)
{
    using Class=SFcalc_obs_base<T>;
    auto pyclass_name = std::string("_SFcalc_obs_base_") + dtype;
    py::class_<Class, std::unique_ptr<Class, py::nodelete>>(m, pyclass_name.c_str())
        .def("__call__",
            [](Class& self,
               HKL_data<F_phi<T>>& fphi,
               const HKL_data<F_sigF<T>>& fsigf,
               const Atom_list& atoms)
             { return self(fphi, fsigf, atoms); })
        ;
} // declare_sfcalc_obs_base

template<class T>
void declare_sfcalc_obs_bulk(py::module& m, const char* dtype)
{
    using Class=SFcalc_obs_bulk<T>;
    auto pyclass_name = std::string("SFcalc_obs_bulk_") + dtype;
    py::class_<Class, SFcalc_obs_base<T>>(m, pyclass_name.c_str())
        .def(py::init<const int>(), py::arg("n_params") = 12)
        .def(py::init<HKL_data<F_phi<T>>&, const HKL_data<F_sigF<T>>&, const Atom_list&, const int>(),
            py::arg("f_phi_out"), py::arg("fsigf"), py::arg("atoms"), py::arg("n_params") = 12)
        .def_property_readonly("bulk_frac", &Class::bulk_frac)
        .def_property_readonly("bulk_scale", &Class::bulk_scale)
        ;
}

void init_sfcalc_obs(py::module& m)
{
    declare_sfcalc_obs_base<ftype32>(m, "float");
    declare_sfcalc_obs_bulk<ftype32>(m, "float");

    declare_sfcalc_obs_base<ftype64>(m, "double");
    declare_sfcalc_obs_bulk<ftype64>(m, "double");

} // init_sfcalc_obs
