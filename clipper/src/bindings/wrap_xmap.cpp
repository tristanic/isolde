#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

#include "type_conversions.h"
#include <clipper/clipper.h>
//#include <clipper/core/map_interp.h>
#include "numpy_helper.h"

namespace py=pybind11;
using namespace clipper;

// Handing these reference coordinates in an efficient way is tricky. On the
// C++ side, the increment/decrement functions return a reference back to
// themselves for convenience. Fine on the C++ side, but a naive wrapping will
// cause the creation of a brand new python object on every iteration - which
// will make things agonizingly slow. A much better option is to wrap all
// increment/decrement functions in lambdas that return void.


void declare_xmap_base(py::module& m)
{
    py::class_<Xmap_base>(m, "_Xmap_base")
        .def_property_readonly("is_null", &Xmap_base::is_null)
        .def_property_readonly("cell", &Xmap_base::cell)
        .def_property_readonly("spacegroup", &Xmap_base::spacegroup)
        .def_property_readonly("grid_sampling", &Xmap_base::grid_sampling)
        .def_property_readonly("grid_asu", &Xmap_base::grid_asu )
        .def("coord_of", &Xmap_base::coord_of)
        .def("index_of", &Xmap_base::index_of)
        .def("to_map_unit", &Xmap_base::to_map_unit)
        .def_property_readonly("operator_orth_grid", &Xmap_base::operator_orth_grid)
        .def_property_readonly("operator_grid_orth", &Xmap_base::operator_grid_orth)
        .def("coord_orth", &Xmap_base::coord_orth)
        .def("coord_map", &Xmap_base::coord_map)
        //.def("in_map", [](const Derived& self, const Coord_grid& cg) { return self.in_map(cg); })
        .def("multiplicity", &Xmap_base::multiplicity)
        .def_property_readonly("first", &Xmap_base::first)
        .def_property_readonly("first_coord", &Xmap_base::first_coord)
        // Let Xmap take control of creation of Map_reference... types
        .def("map_reference_index", [](const Xmap_base& self, const Coord_grid& pos)
        {
            return Xmap_base::Map_reference_index(self, pos);
        })
        .def("map_reference_coord", [](const Xmap_base& self, const Coord_grid& pos)
        {
            return Xmap_base::Map_reference_coord(self, pos);
        })
        ;
} // apply_xmap_base_methods



void declare_xmap_reference_index(py::module &m)
{
    using Class = Xmap_base::Map_reference_index;
    py::class_<Class>(m, "Xmap_reference_index")
        // Let Xmap take control of creation
        // .def(py::init<>())
        // .def(py::init<const Xmap_base&>())
        // .def(py::init<const Xmap_base&, const Coord_grid&>())
        .def_property("coord",
            &Class::coord,
            [](Class& self, const Coord_grid& pos) -> void { self.set_coord(pos); }
        )
        .def_property_readonly("coord_orth", &Class::coord_orth)
        .def("next", [](Class& self) -> void { self.next(); })
        .def("index_offset", &Class::index_offset)
        // base class methods
        .def_property_readonly("base_xmap", [](const Class& self) { return self.base_xmap(); })
        .def_property_readonly("index", [](const Class& self) { return self.index(); })
        .def("last", [](const Class& self) { return self.last(); })
        ;
}

void declare_xmap_reference_coord(py::module &m)
{
    using Class = Xmap_base::Map_reference_coord;
    py::class_<Class>(m, "Xmap_reference_coord")
        // Let Xmap take control of creation
        // .def(py::init<>())
        // .def(py::init<const Xmap_base&>())
        // .def(py::init<const Xmap_base&, const Coord_grid&>())
        .def_property("coord",
            &Class::coord,
            [](Class& self, const Coord_grid& pos) -> void { self.set_coord(pos); }
        )
        .def("coord_orth", &Class::coord_orth)
        .def_property_readonly("sym", &Class::sym)
        .def("next", [](Class& self) -> void { self.next(); })
        .def("next_u", [](Class& self) -> void { self.next_u(); })
        .def("next_v", [](Class& self) -> void { self.next_v(); })
        .def("next_w", [](Class& self) -> void { self.next_w(); })
        .def("prev_u", [](Class& self) -> void { self.prev_u(); })
        .def("prev_v", [](Class& self) -> void { self.prev_v(); })
        .def("prev_w", [](Class& self) -> void { self.prev_w(); })
        // base class methods
        .def_property_readonly("base_xmap", [](const Class& self) { return self.base_xmap(); })
        .def_property_readonly("index", [](const Class& self) { return self.index(); })
        .def("last", [](const Class& self) { return self.last(); })
        ;
}

#include <iostream> //DELETEME
template<typename T>
void numpy_export_core_(const Xmap<T>& xmap, py::array_t<T> target, const Coord_grid& origin)
{
    auto tbuf = target.request();
    T* tptr = (T*)tbuf.ptr;
    int nu, nv, nw;
    nu=tbuf.shape[0]; nv=tbuf.shape[1]; nw=tbuf.shape[2];
    int u,v;
    int maxu, maxv, maxw;
    maxu = origin.u()+nu; maxv=origin.v()+nv; maxw=origin.w()+nw;
    Xmap_base::Map_reference_coord ix(xmap);
    for (u=origin.u(); u<maxu; ++u)
        for (v=origin.v(); v<maxv; ++v)
            for (ix.set_coord(Coord_grid(u, v, origin.w())); ix.coord().w() < maxw; ix.next_w())
                *tptr++ = xmap[ix];
}

template<typename T>
void numpy_import_core_(const Xmap<T>& xmap, py::array_t<T> vals, const Coord_grid& origin)
{
    auto vbuf = vals.request();
    T* vptr = (T*)vbuf.ptr;
    int nu, nv, nw;
    nu=vbuf.shape[0]; nv=vbuf.shape[1]; nw=vbuf.shape[2];
    int u, v, w;
    int maxu, maxv, maxw;
    maxu = origin.u()+nu; maxv=origin.v()+nv; maxw=origin.w()+nw;
    Xmap_base::Map_reference_coord ix(xmap);
    for (u=origin.u(); u<maxu; ++u)
        for (v=origin.v(); v<maxv; ++v)
            for (ix.set_coord(Coord_grid(u,v,origin.w())); ix.coord().w()<maxw; ix.next_w())
                xmap[ix] = *vptr++;
}

template<class C, class T>
void add_xmap_numpy_functions(py::class_<C, Xmap_base>& pyclass)
{
    pyclass
        .def("export_numpy", [](const C& self)
        {
            const auto& g = self.grid_asu();
            auto s = g.max()-g.min();
            auto target = py::array_t<T>({s.u(), s.v(), s.w()});
            numpy_export_core_(self, target, g.min());
            return target;
        })
        .def("export_section_numpy", [](const C& self, const Coord_grid& origin, const Coord_grid& size)
        {
            auto target = py::array_t<T>({size[0], size[1], size[2]});
            numpy_export_core_(self, target, origin);
            return target;
        })
        .def("export_section_numpy", [](const C& self, const Coord_grid& origin, py::array_t<T> target )
        { numpy_export_core_(self, target, origin); })
        // TODO: decide how to handle imports
        //.def("import_numpy", [](C& self, py::array_t<T> vals))
        ;
}


template<class Derived, class T>
void apply_xmap_interpolation_methods(py::class_<Derived, Xmap_base>& pyclass)
{
    pyclass
        .def("interp_linear_frac", [&](const Derived& self, const Coord_frac& pos)
            { return self.template interp<Interp_linear>(pos); })
        .def("interp_cubic_frac", [](const Derived& self, const Coord_frac& pos)
            { return self.template interp<Interp_cubic>(pos); })
        .def("interp_linear_orth", [](const Derived& self, const Coord_orth& xyz )
            { return self.template interp<Interp_linear>(xyz.coord_frac(self.cell())); })
        .def("interp_cubic_orth", [](const Derived& self, const Coord_orth& xyz)
            { return self.template interp<Interp_cubic>(xyz.coord_frac(self.cell())); })
        .def("interp_cubic_grad_frac", [](const Derived& self, const Coord_frac& pos) -> py::tuple
        {
            T val;
            Grad_frac<T> grad;
            self.template interp_grad<Interp_cubic>(pos, val, grad);
            return py::make_tuple(val, grad);
        })
        .def("interp_cubic_curv_frac", [](const Derived& self, const Coord_frac& pos) -> py::tuple
        {
            T val;
            Grad_frac<T> grad;
            Curv_frac<T> curv;
            self.template interp_curv<Interp_cubic>(pos, val, grad, curv);
            return py::make_tuple(val, grad, curv);
        })
        ;
} // apply_xmap_interpolation_methods

template <class Derived, class T, class H>
void apply_xmap_fft_methods(py::class_<Derived, Xmap_base>& pyclass)
{
    pyclass
        .def("fft_from", [](Derived& self, const H& fphidata) { self.fft_from(fphidata); })
        .def("fft_to", [](const Derived& self, H& fphidata) { self.fft_to(fphidata); })
        ;
}

template <class T>
void declare_xmap(py::module& m, const char* dtype)
{
    using MRI = Xmap_base::Map_reference_index;
    using MRC = Xmap_base::Map_reference_coord;
    using Class=Xmap<T>;
    auto pyclass_name = std::string("Xmap_") + dtype;
    py::class_<Class, Xmap_base> xmap(m, pyclass_name.c_str());
    xmap
        .def(py::init<>())
        .def(py::init<const Spacegroup&, const Cell&, const Grid_sampling&>())
        .def("init", &Class::init)
        .def("__getitem__", [](const Class& self, const MRI& ix) { return self[ix]; })
        .def("__setitem__", [](Class& self, const MRI& ix, const T& val) { self[ix] = val; })
        .def("__getitem__", [](const Class& self, const MRC& ix) { return self[ix]; })
        .def("__setitem__", [](Class& self, const MRC& ix, const T& val) { self[ix] = val; })
        //TODO: Would these be better off as __getitem__/__setitem__ as well?
        .def("get_data", [](const Class& self, const Coord_grid& pos) { return self.get_data(pos); })
        .def("set_data", [](Class& self, const Coord_grid& pos, const T& val) { self.set_data(pos, val); })
        .def("get_data", [](const Class& self, const int& index) { return self.get_data(index); })
        .def("set_data", [](Class& self, const int& index, const T& val) { self.set_data(index, val); })
        .def("set_all_values_to", [](Class& self, const T& val) { self = val; })
        .def(py::self += py::self)
        .def(py::self -= py::self)
        ;
    apply_xmap_interpolation_methods<Class, T>(xmap);
    apply_xmap_fft_methods<Class, T, HKL_data<clipper::data32::F_phi>>(xmap);
    apply_xmap_fft_methods<Class, T, HKL_data<clipper::data64::F_phi>>(xmap);
    add_xmap_numpy_functions<Class, T>(xmap);
} // declare_xmap



void init_xmap(py::module &m)
{
    declare_xmap_base(m);
    declare_xmap_reference_coord(m);
    declare_xmap_reference_index(m);
    declare_xmap<ftype32>(m, "float");
    declare_xmap<ftype64>(m, "double");
} // init_xmap
