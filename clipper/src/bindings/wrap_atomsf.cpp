#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

#include "type_conversions.h"
#include <clipper/clipper.h>

#include "numpy_helper.h"

namespace py=pybind11;
using namespace clipper;

void init_atomsf(py::module&m)
{
py::class_<AtomShapeFn> atomsf(m, "AtomShapeFn");
atomsf
    .def(py::init<>())
    .def(py::init<const Atom&>())
    .def(py::init<const Coord_orth&, const String&, const ftype, const ftype>())
    .def(py::init<const Coord_orth&, const String&, const U_aniso_orth&, const ftype>())
    .def("init", (void (AtomShapeFn::*)(const Atom&)) &AtomShapeFn::init)
    .def("init", (void (AtomShapeFn::*)(const Coord_orth&, const String&, const ftype, const ftype)) &AtomShapeFn::init)
    .def("init", (void (AtomShapeFn::*)(const Coord_orth&, const String&, const U_aniso_orth&, const ftype)) &AtomShapeFn::init)
    .def("f", (ftype (AtomShapeFn::*)(const Coord_reci_orth&) const) &AtomShapeFn::f)
    .def("f", (ftype (AtomShapeFn::*)(const ftype&) const) &AtomShapeFn::f)
    .def("rho", (ftype (AtomShapeFn::*)(const Coord_orth&) const) &AtomShapeFn::rho)
    .def("rho", (ftype (AtomShapeFn::*)(const ftype&) const) &AtomShapeFn::rho)
    .def("rho_grad", [] (const AtomShapeFn& self, const Coord_orth& xyz)
    {
        ftype rho;
        std::vector<ftype> grad;
        self.rho_grad(xyz, rho, grad);
        return py::make_tuple(rho, grad);
    })
    .def("rho_curv", [] (const AtomShapeFn& self, const Coord_orth& xyz)
    {
        ftype rho;
        std::vector<ftype> grad;
        Matrix<ftype> curv;
        self.rho_curv(xyz, rho, grad, curv);
        return py::make_tuple(rho, grad, curv);
    })
    .def("agarwal_params", &AtomShapeFn::agarwal_params)
    ;


py::enum_<AtomShapeFn::TYPE>(atomsf, "TYPE")
    .value("X", AtomShapeFn::TYPE::X)
    .value("Y", AtomShapeFn::TYPE::Y)
    .value("Z", AtomShapeFn::TYPE::Z)
    .value("Uiso", AtomShapeFn::TYPE::Uiso)
    .value("Occ", AtomShapeFn::TYPE::Occ)
    .value("U11", AtomShapeFn::TYPE::U11)
    .value("U22", AtomShapeFn::TYPE::U22)
    .value("U33", AtomShapeFn::TYPE::U33)
    .value("U12", AtomShapeFn::TYPE::U12)
    .value("U13", AtomShapeFn::TYPE::U13)
    .value("U23", AtomShapeFn::TYPE::U23)
    .export_values();


} // init_atomsf
