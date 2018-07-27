#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>


#include <clipper/clipper.h>
#include "numpy_helper.h"

namespace py=pybind11;
using namespace clipper;


template<typename T>
void numpy_export_core_(const NXmap<T>& nxmap,
    py::array_t<T> target,
    const Coord_grid& origin)
{
    auto tbuf = target.request();
    T* tptr = (T*)tbuf.ptr;
    Coord_grid duvw(tbuf.shape[0], tbuf.shape[1], tbuf.shape[2]);
    Coord_grid last = origin + duvw;
    if (!nxmap.in_map(origin))
        throw std::out_of_range("Requested origin is outside the map!");
    if (!nxmap.in_map(last))
    {
        auto g = nxmap.grid();
        last = Coord_grid(g.nu(), g.nv(), g.nw());
    }
    Coord_grid c;
    for (c.u() = origin.u(); c.u() < last.u(); c.u()++)
        for (c.v() = origin.v(); c.v() < last.v(); c.v()++)
            for (c.w() = origin.w(); c.w() < last.w(); c.w()++)
                *tptr++ = nxmap.get_data(c);
}

template<typename T>
void numpy_import_core_(NXmap<T>& nxmap,
    py::array_t<T> vals,
    const Coord_grid& origin)
{
    auto vbuf = vals.request();
    T* vptr = (T*)vbuf.ptr;
    Coord_grid duvw(vbuf.shape[0], vbuf.shape[1], vbuf.shape[2]);
    Coord_grid last = origin + duvw;
    if (!nxmap.in_map(origin))
        throw std::out_of_range("Input origin is outside the map!");
    if (!nxmap.in_map(last))
    {
        auto g = nxmap.grid();
        last = Coord_grid(g.nu(), g.nv(), g.nw());
    }
    Coord_grid c;
    for (c.u() = origin.u(); c.u() < last.u(); c.u()++)
        for (c.v() = origin.v(); c.v() < last.v(); c.v()++)
            for (c.w() = origin.w(); c.w() < last.w(); c.w()++)
                nxmap.set_data(c, *vptr++);
}



template<class C, class T>
void add_nxmap_numpy_functions(py::class_<C>& pyclass)
{
    pyclass
        .def("export_numpy", [](const C& self)
        {
            auto g = self.grid();
            auto target = py::array_t<T>({g[0],g[1],g[2]});
            numpy_export_core_(self, target, Coord_grid(0,0,0));
            return target;
        },
        "Export the whole map to a Numpy array.")
        .def("export_numpy", [](const C& self, py::array_t<T> target)
        { numpy_export_core_(self, target, Coord_grid(0,0,0)); },
        "Export the map into the given numpy array, starting at the origin. "
        "If the target array is smaller than the map, the output will be "
        "truncated to fit.")
        .def("export_fragment_numpy", [](const C& self, const Coord_grid& origin,
            const Coord_grid& size)
        {
            if (!self.in_map(origin+size))
                throw std::out_of_range("Requested data extends beyond the range of the map!");
            auto target = py::array_t<T>({size[0], size[1], size[2]});
            numpy_export_core_(self, target, origin);
            return target;
        },
        "Export a fragment of the map with the given origin and size, as a "
         "numpy array.")
        .def("export_fragment_numpy", [](const C& self, const Coord_grid& origin,
            py::array_t<T> target)
            { numpy_export_core_(self, target, origin); },
        "Export a fragment of the map with the given origin into the given numpy array.")
        .def("import_numpy", [](C& self, const Coord_grid& origin, py::array_t<T> vals)
            { numpy_import_core_(self, vals, origin); },
            "Import data from numpy. Import will start at the origin of the numpy "
            "array, and be written to the map starting at the given origin.")
        ;
}

void declare_nxmap_reference_index(py::module& m)
{
    using Class = NXmap_base::Map_reference_index;
    py::class_<Class>(m, "NXmap_reference_index")
        .def(py::init<>())
        .def(py::init<const NXmap_base&>())
        .def(py::init<const NXmap_base&, const Coord_grid&>())
        .def("coord", &Class::coord)
        .def("coord_orth", &Class::coord_orth)
        .def("set_coord", &Class::set_coord)
        .def("next", &Class::next)
        .def("index_offset", &Class::index_offset)
        // base class methods
        .def_property_readonly("base_nxmap", [](const Class& self) { return self.base_nxmap(); })
        .def_property_readonly("index", [](const Class& self) { return self.index(); })
        .def("last", [](const Class& self){ return self.last(); })
        ;
}

void declare_nxmap_reference_coord(py::module& m)
{
    using Class = NXmap_base::Map_reference_coord;
    py::class_<Class>(m, "NXmap_reference_coord")
        .def(py::init<>())
        .def(py::init<const NXmap_base&>())
        .def(py::init<const NXmap_base&, const Coord_grid&>())
        .def("coord", &Class::coord)
        .def("coord_orth", &Class::coord_orth)
        .def("set_coord", &Class::set_coord)
        .def("next", &Class::next)
        .def("next_u", &Class::next_u)
        .def("next_v", &Class::next_v)
        .def("next_w", &Class::next_w)
        .def("prev_u", &Class::prev_u)
        .def("prev_v", &Class::prev_v)
        .def("prev_w", &Class::prev_w)
        // base class methods
        .def_property_readonly("base_nxmap", [](const Class& self) { return self.base_nxmap(); })
        .def_property_readonly("index", [](const Class& self) { return self.index(); })
        .def("last", [](const Class& self){ return self.last(); })
        ;
}

template<class Derived, class T>
void apply_nxmap_base_methods(py::class_<Derived>& pyclass)
{
    pyclass
        .def("is_null", [](const Derived& self) { return self.is_null(); })
        .def("grid", [](const Derived& self) { return self.grid(); })
        .def("operator_orth_grid", [](const Derived& self){ return self.operator_orth_grid(); })
        .def("operator_grid_orth", [](const Derived& self) { return self.operator_grid_orth(); })
        .def("coord_orth", [](const Derived& self, const Coord_map& cm) { return self.coord_orth(cm); })
        .def("coord_map", [](const Derived& self, const Coord_orth& co) { return self.coord_map(co); })
        .def("in_map", [](const Derived& self, const Coord_grid& pos) { return self.in_map(pos); })
        // TODO: handle in_map() for different interpolators
        .def("multiplicity", [](const Derived& self, const Coord_grid&) { return 1; })

        .def("first", [](const Derived& self) { return self.first(); })
        .def("first_coord", [](const Derived& self) { return self.first_coord(); })
        ;
}

template<class T>
void declare_nxmap(py::module& m, const char* dtype)
{
    using Class = NXmap<T>;
    auto pyclass_name = std::string("NXmap_") + dtype;
    py::class_<Class> nxmap(m, pyclass_name.c_str());
    nxmap
        .def(py::init<>())
        .def(py::init<const Grid&, const RTop<>&>())
        .def(py::init<const Cell&, const Grid_sampling&, const Grid_range&>())
        .def("init", (void (Class::*)(const Grid&, const RTop<>&)) &Class::init)
        .def("init", (void (Class::*)(const Cell&, const Grid_sampling&, const Grid_range&)) &Class::init)
        .def("__getitem__", [](const Class& self, const NXmap_base::Map_reference_index i) { return self[i]; })
        .def("__setitem__", [](Class& self, const NXmap_base::Map_reference_index i, const T& val){ self[i] = val; })
        .def("__getitem__", [](const Class& self, const NXmap_base::Map_reference_coord i) { return self[i]; })
        .def("__setitem__", [](Class& self, const NXmap_base::Map_reference_coord i, const T& val) { self[i] = val; })
        .def("get_data", [](const Class& self, const Coord_grid& pos) { return self.get_data(pos); })
        .def("set_data", [](Class& self, const Coord_grid& pos, const T& val) { self.set_data(pos, val); })
        .def("set_all_values_to", [](Class& self, const T& val) { self = val; })
        .def(py::self += py::self)
        .def(py::self -= py::self)
        ;
    apply_nxmap_base_methods<Class, T>(nxmap);
    add_nxmap_numpy_functions<Class, T>(nxmap);
}

void init_nxmap(py::module& m)
{
    declare_nxmap_reference_index(m);
    declare_nxmap_reference_coord(m);
    declare_nxmap<ftype32>(m, "float");
    declare_nxmap<ftype64>(m, "double");
}
