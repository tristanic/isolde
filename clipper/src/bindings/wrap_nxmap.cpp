#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

#include "type_conversions.h"
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
            auto target = py::array_t<T>({g.nu(),g.nv(),g.nw()});
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

// Handing these reference coordinates in an efficient way is tricky. On the
// C++ side, the increment/decrement functions return a reference back to
// themselves for convenience. Fine on the C++ side, but a naive wrapping will
// cause the creation of a brand new python object on every iteration - which
// will make things agonizingly slow. A much better option is to wrap all
// increment/decrement functions in lambdas that return void.

void declare_nxmap_reference_index(py::module& m)
{
    using Class = NXmap_base::Map_reference_index;
    py::class_<Class>(m, "NXmap_reference_index")
        // Let NXmap take control of creation
        // .def(py::init<>())
        // .def(py::init<const NXmap_base&>())
        // .def(py::init<const NXmap_base&, const Coord_grid&>())
        .def_property("coord",
            &Class::coord,
            [](Class& self, const Coord_grid& pos) -> void { self.set_coord(pos); }
        )
        .def_property_readonly("coord_orth", &Class::coord_orth)
        .def("next", [](Class& self) -> void { self.next(); })
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
        // Let NXmap take control of creation
        // .def(py::init<>())
        // .def(py::init<const NXmap_base&>())
        // .def(py::init<const NXmap_base&, const Coord_grid&>())
        .def_property("coord",
            &Class::coord,
            [](Class& self, const Coord_grid& pos) -> void { self.set_coord(pos); }
        )
        .def_property_readonly("coord_orth", &Class::coord_orth)
        .def("next", [](Class& self) -> void { self.next(); })
        .def("next_u", [](Class& self) -> void { self.next_u(); })
        .def("next_v", [](Class& self) -> void { self.next_v(); })
        .def("next_w", [](Class& self) -> void { self.next_w(); })
        .def("prev_u", [](Class& self) -> void { self.prev_u(); })
        .def("prev_v", [](Class& self) -> void { self.prev_v(); })
        .def("prev_w", [](Class& self) -> void { self.prev_w(); })
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
        // Let NXmap take control of creation of Map_reference... types
        .def("map_reference_index", [](const Derived& self, const Coord_grid& pos)
        {
            return NXmap_base::Map_reference_index(self, pos);
        })
        .def("map_reference_coord", [](const Derived& self, const Coord_grid& pos)
        {
            return NXmap_base::Map_reference_coord(self, pos);
        })
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
        //TODO: Would these be better off as __getitem__/__setitem__ as well?
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
