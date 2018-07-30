#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

#include "type_conversions.h"
#include <clipper/clipper.h>

#include "numpy_helper.h"

namespace py=pybind11;
using namespace clipper;

void init_cell(py::module& m)
{

py::class_<Metric_tensor>(m, "Metric_tensor")
    .def(py::init<>())
    .def(py::init<const ftype&, const ftype&,const ftype&, const ftype&,const ftype&, const ftype&>() )
    .def("lengthsq", (ftype (Metric_tensor::*)(const Vec3<>&) const) &Metric_tensor::lengthsq)
    .def("lengthsq", (ftype (Metric_tensor::*)(const Vec3<int>&) const) &Metric_tensor::lengthsq)
    .def("format", &Metric_tensor::format)
    .def("__str__", [](const Metric_tensor& self) { return self.format().c_str(); })
    ;

py::class_<Cell_descr>(m, "Cell_descr")
    .def(py::init<>())
    .def(py::init<const ftype&, const ftype&, const ftype&, const ftype&, const ftype&, const ftype&>())
    .def_property_readonly("a", &Cell_descr::a)
    .def_property_readonly("b", &Cell_descr::b)
    .def_property_readonly("c", &Cell_descr::c)
    .def_property_readonly("dim", [](const Cell_descr& self)
    {
        auto ret = py::array_t<ftype>(3);
        auto buf = ret.request();
        ftype* ptr = (ftype*)buf.ptr;
        ptr[0]=self.a(); ptr[1]=self.b(); ptr[2] = self.c();
        return ret;
    })
    .def_property_readonly("alpha", &Cell_descr::alpha)
    .def_property_readonly("beta", &Cell_descr::beta)
    .def_property_readonly("gamma", &Cell_descr::gamma)
    .def_property_readonly("angles", [](const Cell_descr& self)
    {
        auto ret = py::array_t<ftype>(3);
        auto buf = ret.request();
        ftype* ptr = (ftype*)buf.ptr;
        ptr[0]=self.alpha(); ptr[1]=self.beta(); ptr[2] = self.gamma();
        return ret;
    })
    .def_property_readonly("alpha_deg", &Cell_descr::alpha_deg)
    .def_property_readonly("beta_deg", &Cell_descr::beta_deg)
    .def_property_readonly("gamma_deg", &Cell_descr::gamma_deg)
    .def_property_readonly("angles_deg", [](const Cell_descr& self)
    {
        auto ret = py::array_t<ftype>(3);
        auto buf = ret.request();
        ftype* ptr = (ftype*)buf.ptr;
        ptr[0]=self.alpha_deg(); ptr[1]=self.beta_deg(); ptr[2] = self.gamma_deg();
        return ret;
    })
    .def("format", &Cell_descr::format)
    .def("__str__", [](const Cell_descr& self) { return self.format().c_str(); })
    ;

py::class_<Cell, Cell_descr>(m, "Cell")
    .def(py::init<>())
    .def(py::init<const Cell_descr&>())
    .def("init", &Cell::init)
    .def("is_null", &Cell::is_null)
    .def_property_readonly("a_star", &Cell::a_star)
    .def_property_readonly("b_star", &Cell::b_star)
    .def_property_readonly("c_star", &Cell::c_star)
    .def_property_readonly("recip_dim", [](const Cell& self)
    {
        auto ret = py::array_t<ftype>(3);
        auto buf = ret.request();
        ftype* ptr = (ftype*)buf.ptr;
        ptr[0]=self.a_star(); ptr[1]=self.b_star(); ptr[2] = self.c_star();
        return ret;
    })
    .def_property_readonly("alpha_star", &Cell::alpha_star)
    .def_property_readonly("beta_star", &Cell::beta_star)
    .def_property_readonly("gamma_star", &Cell::gamma_star)
    .def_property_readonly("recip_angles", [](const Cell& self)
    {
        auto ret = py::array_t<ftype>(3);
        auto buf = ret.request();
        ftype* ptr = (ftype*)buf.ptr;
        ptr[0]=self.alpha_star(); ptr[1]=self.beta_star(); ptr[2] = self.gamma_star();
        return ret;
    })
    .def_property_readonly("descr", &Cell::descr)
    .def_property_readonly("volume", &Cell::volume)
    .def("equals", &Cell::equals)
    .def_property_readonly("matrix_orth", &Cell::matrix_orth)
    .def_property_readonly("matrix_frac", &Cell::matrix_frac)
    .def_property_readonly("metric_real", &Cell::metric_real)
    .def_property_readonly("metric_reci", &Cell::metric_reci)
    ;

} // init_cell
