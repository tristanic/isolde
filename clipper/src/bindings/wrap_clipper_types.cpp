#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <clipper/clipper.h>


namespace py=pybind11;
using namespace clipper;

template<class T>
void init_vec3(py::module &m)
{
    using Class=Vec3<T>;
    py::class_<Class>(m, "Vec3")
        .def(py::init<>())
        .def(py::init<const T&, const T&, const T&>())
        .def(py::init<const Vec3<ftype32>&>())
        .def(py::init<const Vec3<ftype64>&>())
        .def(py::init<const Vec3<int>&>())
        .def("equals", &Class::equals)
        .def("__getitem__", [](const Class& self, const int& i){ return self[i]; })
        .def("__setitem__", [](Class& self, const int& i, const T& val) { self[i] = val; })
        .def("unit", &Class::unit)
        .def_static("zero", &Class::zero)
        .def_static("null", &Class::null)
        .def("is_null", &Class::is_null)
        .def_static("dot", &Class::dot)
        .def_static("cross", &Class::cross)
        .def("__str__", &Class::format)
        .def("format", &Class::format)
        .def("__iadd__", &Class::operator+=)
        .def("__isub__", &Class::operator-=)
        .def("__eq__", [](const Class& self, const Class& other){ return self == other; }, py::is_operator())
        .def("__ne__", [](const Class& self, const Class& other){ return self != other; }, py::is_operator())
        .def("__neg__", [](const Class& self){ return -self; })
        .def("__add__", [](const Class& self, const Class& other){ return self + other; }, py::is_operator())
        .def("__sub__", [](const Class& self, const Class& other){ return self + (-other); }, py::is_operator())
        .def("__rmul__", [](const Class& self, const T& s){ return s*self; }, py::is_operator())
        .def("__mul__", [](const Class& self, const T& s){ return self*s; }, py::is_operator())
        ;
} // init_vec3

template<class T>
void init_mat33(py::module &m)
{
    using Class=Mat33<T>;
    py::class_<Class>(m, "Mat33")
        .def(py::init<>())
        .def(py::init<const T&, const T&, const T&, const T&, const T&, const T&, const T&, const T&, const T&>())
        .def(py::init<const Mat33<ftype32>&>())
        .def(py::init<const Mat33<ftype64>&>())
        .def(py::init<const Mat33sym<ftype32>&>())
        .def(py::init<const Mat33sym<ftype64>&>())
        .def("det", &Class::det)
        .def("inverse", &Class::inverse)
        .def("transpose", &Class::transpose)
        .def("equals", &Class::equals)
        .def("get", [](const Class& self, const int& i, const int& j) { return self(i,j); })
        .def("set", [](Class& self, const int& i, const int& j, const T& val){ self(i,j) = val; })
        .def("__repr__", [](const Class& self){ return std::string("Clipper Mat33: ") + std::string(self.format()); })
        .def("__str__", [](const Class& self){ return self.format(); })
        .def_static("identity", &Class::identity)
        .def_static("null", &Class::null)
        .def("is_null", &Class::is_null)
        .def("__mul__", [](const Class& self, const Vec3<T>& v){ return self*v; }, py::is_operator())
        .def("__rmul__", [](const Class& self, const Vec3<T>& v){ return v*self;}, py::is_operator())
        .def("__mul__", [](const Class& self, const Class& other){ return self*other; }, py::is_operator())
        .def("__add__", [](const Class& self, const Class& other){ return self+other; }, py::is_operator())
        .def("__sub__", [](const Class& self, const Class& other){ return self + (-other); }, py::is_operator())
        .def("__neg__", [](const Class& self){ return -self; }, py::is_operator())
        ;
} // init_mat3

template <typename T>
void init_mat33sym(py::module &m)
{
    using Class=Mat33sym<T>;
    py::class_<Class>(m, "Mat33sym")
        .def(py::init<>())
        .def(py::init<const Mat33<ftype32>&>())
        .def(py::init<const Mat33<ftype64>&>())
        .def(py::init<const Mat33sym<ftype32>&>())
        .def(py::init<const Mat33sym<ftype64>&>())
        .def(py::init<const T&, const T&, const T&, const T&, const T&, const T&>())
        .def("__repr__", [](const Class& self){ return std::string("Clipper Mat33sym: ")+ std::string(self.format()); })
        .def("__str__", [](const Class& self){ return std::string(self.format()); })
        .def_static("identity", &Class::identity)
        .def_static("null", &Class::null)
        .def("is_null", &Class::is_null)
        .def("quad_form", &Class::quad_form)
        .def("det", &Class::det)
        .def("sqrt", &Class::sqrt)
        .def("inverse", &Class::inverse)
        .def("mat00", &Class::mat00)
        .def("mat11", &Class::mat11)
        .def("mat12", &Class::mat12)
        .def("mat01", &Class::mat01)
        .def("mat02", &Class::mat02)
        .def("mat12", &Class::mat12)
        .def("get", [](const Class& self, const int& i, const int& j) { return self(i,j); })
        .def("__mul__", [](const Class& self, const Vec3<T>& v){ return self*v; }, py::is_operator())
        .def("__add__", [](const Class& self, const Class& other){ return self + other; }, py::is_operator())
        .def("__sub__", [](const Class& self, const Class& other){ return self + (-other); }, py::is_operator())
        .def("__neg__", [](const Class& self){ return -self; }, py::is_operator())
        ;
} // init_mat33sym

template<typename T>
void init_rtop(py::module &m)
{
    using Class=RTop<T>;
    py::class_<Class>(m, "RTop")
        .def(py::init<>())
        .def(py::init<const Mat33<T>&>())
        .def(py::init<const Mat33<T>&, const Vec3<T>&>())
        .def("inverse", &Class::inverse)
        ;
} // init_rtop

void init_clipper_types(py::module &m, py::module &m32, py::module &m64) {

    py::class_<String>(m, "String")
        .def(py::init<>())
        .def(py::init<const std::string>())
        .def(py::init<const char*>())
        .def(py::init<const char*, const int>())
        .def(py::init<const int, const int>())
        .def(py::init<const long, const int>())
        .def(py::init<const float, const int, const int>())
        .def(py::init<const double, const int, const int>())
        .def("split", &String::split)
        .def("trim", &String::trim)
        .def("tail", &String::tail)
        .def("head", &String::head)
        .def("nohead", &String::nohead)
        .def("notail", &String::notail)
        .def_static("rational", [](const double& f, const int& b, bool sign=false){ return String::rational(f, b, sign); })
        .def("__int__", &String::l)
        .def("__float__", &String::f64)
        .def("__str__", [](const String& self){ return std::string(self); })
        .def("__repr__", [](const String& self) { return std::string(self); });

    init_vec3<ftype32>(m32);
    init_mat33<ftype32>(m32);
    init_mat33sym<ftype32>(m32);

    init_vec3<ftype64>(m64);
    init_mat33<ftype64>(m64);
    init_mat33sym<ftype64>(m64);
}
