#include <pybind11/pybind11.h>

#include "../type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>

namespace py=pybind11;
using namespace clipper;
using namespace clipper::datatypes;

template <class T>
void declare_mapfilter_base(py::module& m, const char* dtype)
{
    using Class=MapFilter_base<T>;
    auto pyclass_name = std::string("_MapFilter_base_") + dtype;
    py::class_<Class, std::unique_ptr<Class, py::nodelete>>(m, pyclass_name.c_str())
        .def("__call__",
            [](const Class& self, Xmap<T>& result, Xmap<T>& xmap)
            { return self(result, xmap); })
        ;
} // declare_mapfilter_base

void declare_mapfilterfn_base(py::module& m)
{
    using Class=MapFilterFn_base;
    py::class_<Class, std::unique_ptr<Class, py::nodelete>>(m, "_MapFilterFn_base")
        .def("__call__", [](const Class& self, const ftype& radius) { return self(radius); })
        ;
}

template<class Class>
void declare_mapfilterfn_specialization(py::module& m, const char* pyclass_name)
{
    py::class_<Class, MapFilterFn_base>(m, pyclass_name)
        .def(py::init<const ftype&>())
        ;
}

template<class Class, class T>
void declare_mapfilter_specialization(py::module& m, const char* pyclass_name)
{
    py::class_<Class, MapFilter_base<T>> pyclass(m, pyclass_name);
    py::enum_<typename Class::TYPE>(pyclass, "SCALE_MODE")
        .value("NONE", Class::TYPE::NONE)
        .value("Absolute", Class::TYPE::Absolute)
        .value("Relative", Class::TYPE::Relative)
        .export_values()
        ;
        
    pyclass
        .def(py::init<const MapFilterFn_base&, const ftype, const typename Class::TYPE>(),
            py::arg("filter"), py::arg("scale")=1.0, py::arg("scale_mode") = Class::NONE)
        .def(py::init<Xmap<T>&, const Xmap<T>&, MapFilterFn_base&, const ftype, const typename Class::TYPE>(),
            py::arg("result"), py::arg("xmap"), py::arg("filter"), py::arg("scale")=1.0, py::arg("scale_mode")=Class::NONE)
        ;

}

void init_mapfilter(py::module& m)
{
    declare_mapfilterfn_base(m);
    declare_mapfilterfn_specialization<MapFilterFn_step>(m, "MapFilterFn_step");
    declare_mapfilterfn_specialization<MapFilterFn_linear>(m, "MapFilterFn_linear");
    declare_mapfilterfn_specialization<MapFilterFn_quadratic>(m, "MapFilterFn_quadratic");


    declare_mapfilter_base<ftype32>(m, "float");
    declare_mapfilter_specialization<MapFilter_slow<ftype32>, ftype32>(m, "MapFilter_slow_float");
    declare_mapfilter_specialization<MapFilter_fft<ftype32>, ftype32>(m, "MapFilter_fft_float");

    declare_mapfilter_base<ftype64>(m, "double");
    declare_mapfilter_specialization<MapFilter_slow<ftype64>, ftype64>(m, "MapFilter_slow_double");
    declare_mapfilter_specialization<MapFilter_fft<ftype64>, ftype64>(m, "MapFilter_fft_double");
} // init_mapfilter
