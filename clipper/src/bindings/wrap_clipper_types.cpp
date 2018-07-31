#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

#include "type_conversions.h"
#include <clipper/clipper.h>

#include "numpy_helper.h"

namespace py=pybind11;
using namespace clipper;


template<class T>
void init_vec3(py::module &m, const std::string& dtype)
{
    using Class=Vec3<T>;
    auto pyclass_name = std::string("Vec3_")+dtype;
    py::class_<Class>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const T&, const T&, const T&>())
        // From a numpy array
        .def(py::init([](py::array_t<T> vals)
        {
            check_numpy_array_shape(vals, {3}, true);
            T* ptr = (T*)vals.request().ptr;
            return std::unique_ptr<Class>(new Class(ptr[0], ptr[1], ptr[2]));
        }))
        .def(py::init<const Vec3<ftype32>&>())
        .def(py::init<const Vec3<ftype64>&>())
        .def(py::init<const Vec3<int>&>())
        .def("equals", &Class::equals)
        .def("__getitem__", [](const Class& self, const int& i){ return self[i]; })
        .def("as_numpy", [](const Class& self) { return array_as_numpy_1d<Class, T>(self, 3); })
        .def("as_numpy", [](const Class& self, py::array_t<T> target) { array_as_numpy_1d(self, 3, target); })
        // .def("as_numpy", [](const Class& self) -> py::array_t<T>
        //     {
        //         auto result = py::array_t<T>(3);
        //         auto buf=result.request();
        //         T* ptr = (T*)buf.ptr;
        //         for (size_t i=0; i<3; ++i) {
        //             ptr[i] = self[i];
        //         }
        //         return result;
        //     })
        // .def("as_numpy", [](const Class& self, py::array_t<T> result)
        //     {
        //         check_numpy_array_shape(result, {3}, false);
        //         auto buf = result.request();
        //         T* ptr = (T*)buf.ptr;
        //         for(size_t i=0; i<3; ++i)
        //             ptr[i] = self[i];
        //     })
        .def("__setitem__", [](Class& self, const int& i, const T& val) { self[i] = val; })
        .def("from_numpy", [](Class& self, py::array_t<T> vals) { fill_array_from_numpy_1d<Class, T>(self, 3, vals); })
        // .def("from_numpy", [](Class& self, py::array_t<T> vals)
        //     {
        //         check_numpy_array_shape(vals, {3}, true);
        //         auto buf = vals.request();
        //         T* ptr = (T*)buf.ptr;
        //         for (size_t i=0; i<3; ++i)
        //             self[i] = ptr[i];
        //     })
        .def("unit", &Class::unit)
        .def_static("zero", &Class::zero)
        .def_static("null", &Class::null)
        .def("is_null", &Class::is_null)
        .def_static("dot", &Class::dot)
        .def_static("cross", &Class::cross)
        //.def("__repr__", [](const Class& self) { return self.format().c_str(); })
        .def("__str__", [](const Class& self) { return self.format().c_str(); })
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


template<typename T>
void init_mat33(py::module &m, const std::string& dtype)
{
    using Class=Mat33<T>;
    auto pyclass_name = std::string("Mat33_")+dtype;
    py::class_<Class>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const T&, const T&, const T&, const T&, const T&, const T&, const T&, const T&, const T&>())
        // From a 3x3 numpy array
        .def(py::init([](py::array_t<T> vals)
        {
            check_numpy_array_shape(vals, {3,3}, true);
            auto buf = vals.request();
            T* ptr = (T*)buf.ptr;
            return std::unique_ptr<Class>(new Class(
                ptr[0],ptr[1],ptr[2],ptr[3],ptr[4],ptr[5],ptr[6],ptr[7],ptr[8]));
        }))
        .def(py::init<const Mat33<ftype32>&>())
        .def(py::init<const Mat33<ftype64>&>())
        .def(py::init<const Mat33sym<ftype32>&>())
        .def(py::init<const Mat33sym<ftype64>&>())
        .def("det", &Class::det)
        .def("inverse", &Class::inverse)
        .def("transpose", &Class::transpose)
        .def("equals", &Class::equals)
        .def("get", [](const Class& self, const int& i, const int& j) { return self(i,j); })
        .def("as_numpy", [](const Class& self)
        {
            auto result = py::array_t<T>({3,3});
            auto buf = result.request();
            T* ptr = (T*)buf.ptr;
            for (size_t i=0; i<3; ++i)
                for (size_t j=0; j<3; ++j)
                    *ptr++=self(i,j);
            return result;
        })
        .def("as_numpy", [](const Class& self, py::array_t<T> result)
        {
            check_numpy_array_shape(result, {3,3}, false);
            auto buf = result.request();
            T* ptr = (T*)buf.ptr;
            for (size_t i=0; i<3; ++i)
                for (size_t j=0; j<3; ++j)
                    *ptr++ = self(i,j);
        })
        .def("set", [](Class& self, const int& i, const int& j, const T& val){ self(i,j) = val; })
        .def("from_numpy", [](Class& self, py::array_t<T> vals)
        {
            check_numpy_array_shape(vals, {3,3}, true);
            auto buf = vals.request();
            T* ptr = (T*)buf.ptr;
            for (size_t i=0; i<3; ++i)
                for (size_t j=0; j<3; ++j)
                    self(i,j) = *ptr++;
        })
        //.def("__repr__", [](const Class& self) { return self.format().c_str(); })
        .def("__str__", [](const Class& self) { return self.format().c_str(); })
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
} // init_mat33

template <typename T>
void init_mat33sym(py::module &m, const std::string& dtype)
{
    using Class=Mat33sym<T>;
    auto pyclass_name = std::string("Mat33sym_")+dtype;
    py::class_<Class>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const Mat33<ftype32>&>())
        .def(py::init<const Mat33<ftype64>&>())
        .def(py::init<const Mat33sym<ftype32>&>())
        .def(py::init<const Mat33sym<ftype64>&>())
        .def(py::init<const T&, const T&, const T&, const T&, const T&, const T&>())
        // from numpy array
        .def(py::init([](py::array_t<T> vals)
        {
            check_numpy_array_shape(vals, {6}, true);
            auto buf=vals.request();
            T* ptr = (T*)buf.ptr;
            return std::unique_ptr<Class>(new Class(
                ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5]));
        }))
        //.def("__repr__", [](const Class& self){ return self.format().c_str(); })
        .def("__str__", [](const Class& self){ return self.format().c_str(); })
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
        .def("as_numpy", [](const Class& self)
        {
            auto target = py::array_t<T>(6);
            auto buf=target.request();
            T* ptr = (T*)buf.ptr;
            ptr[0]=self.mat00(); ptr[1]=self.mat11(); ptr[2]=self.mat22();
            ptr[3]=self.mat01(); ptr[4]=self.mat02(); ptr[5]=self.mat12();
            return target;
        })
        .def("as_numpy", [](const Class& self, py::array_t<T> target)
        {
            check_numpy_array_shape(target, {6}, false);
            auto buf=target.request();
            T* ptr = (T*)buf.ptr;
            ptr[0]=self.mat00(); ptr[1]=self.mat11(); ptr[2]=self.mat22();
            ptr[3]=self.mat01(); ptr[4]=self.mat02(); ptr[5]=self.mat12();
        })
        // Mat33sym has no setters
        .def("__mul__", [](const Class& self, const Vec3<T>& v){ return self*v; }, py::is_operator())
        .def("__add__", [](const Class& self, const Class& other){ return self + other; }, py::is_operator())
        .def("__sub__", [](const Class& self, const Class& other){ return self + (-other); }, py::is_operator())
        .def("__neg__", [](const Class& self){ return -self; }, py::is_operator())
        ;
} // init_mat33sym


template<typename T>
void init_rtop(py::module &m, const std::string& dtype)
{
    using Class=RTop<T>;
    std::string pyclass_name = std::string("RTop_")+dtype;
    py::class_<Class>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const Mat33<T>&>())
        .def(py::init<const Mat33<T>&, const Vec3<T>&>())
        .def("inverse", &Class::inverse)
        .def("equals", &Class::equals)
        .def_property("rot",
            [](const Class& self) { return self.rot(); },
            [](Class& self, const Mat33<T>& rot) {self.rot() = rot; }
        )
        .def_property("trn",
            [](const Class& self) { return self.trn(); },
            [](Class& self, const Vec3<T>& trn) {self.trn() = trn; }
        )
        .def_static("identity", &Class::identity)
        .def_static("null", &Class::null)
        .def_property_readonly("is_null", &Class::is_null)
        .def("format", &Class::format)
        .def("__str__", &Class::format)
        .def("__mul__", [](const Class& self, const Vec3<T>& v){ return self*v; }, py::is_operator())
        .def("__mul__", [](const Class& self, const Class& other){ return self*other; }, py::is_operator())
        ;
} // init_rtop


template<typename T>
void init_array2d(py::module &m, const std::string& dtype)
{
    using Class=Array2d<T>;
    std::string pyclass_name = std::string("Array2d_") + dtype;
    py::class_<Class>(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const int&, const int&>())
        .def(py::init<const int&, const int&>())
        .def("resize", (void (Class::*)(const int&, const int&)) &Class::resize)
        .def("resize", (void (Class::*)(const int&, const int&, const T&)) &Class::resize)
        .def("size", &Class::size)
        .def("rows", &Class::rows)
        .def("cols", &Class::cols)
        .def("get", [](const Class& self, const int& i1, const int& i2) { return self(i1, i2); })
        .def("set", [](Class& self, const int& i1, const int& i2, const T& val) { self(i1, i2) = val; })
        .def("as_numpy", [](const Class& self)
        {
            int rows = self.rows();
            int cols = self.cols();
            auto target = py::array_t<T>({rows, cols});
            auto buf=target.request();
            T* ptr = (T*)buf.ptr;
            for (int i=0; i<rows; ++i)
                for (int j=0; j<cols; ++j)
                    *ptr++ = self(i,j);
            return target;
        })
        .def("as_numpy", [](const Class& self, py::array_t<T> target)
        {
            int rows = self.rows();
            int cols = self.cols();
            check_numpy_array_shape(target, {rows, cols}, false);
            auto buf=target.request();
            T* ptr = (T*)buf.ptr;
            for (int i=0; i<rows; ++i)
                for (int j=0; j<cols; ++j)
                    *ptr++ = self(i,j);
        })
        .def("from_numpy", [](Class& self, py::array_t<T> vals)
        {
            int rows = self.rows();
            int cols = self.cols();
            check_numpy_array_shape(vals, {rows, cols}, true);
            auto buf=vals.request();
            T* ptr = (T*)buf.ptr;
            for (int i=0; i<rows; ++i)
                for (int j=0; j<cols; ++j)
                    self(i,j) = *ptr++;
        })
        ;
}

template <class T>
void init_matrix(py::module&m, const std::string& dtype)
{
    using Class=Matrix<T>;
    std::string pyclass_name = std::string("Matrix_") + dtype;
    py::class_<Class, Array2d<T> >(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<const int&, const int&>())
        .def(py::init<const int&, const int&, T>())
        .def("solve", &Class::solve)
        .def("eigen", &Class::eigen)
        .def("__mul__", [](const Class& self, const std::vector<T>& v) { return self*v; })
        ;

}

void init_clipper_types(py::module &m, py::module &m32, py::module &m64) {

    // py::class_<String>(m, "String")
    //     .def(py::init<>())
    //     .def(py::init<const std::string>())
    //     .def(py::init<const char*>())
    //     .def(py::init<const char*, const int>())
    //     .def(py::init<const int, const int>())
    //     .def(py::init<const long, const int>())
    //     .def(py::init<const float, const int, const int>())
    //     .def(py::init<const double, const int, const int>())
    //     .def("split", &String::split)
    //     .def("trim", &String::trim)
    //     .def("tail", &String::tail)
    //     .def("head", &String::head)
    //     .def("nohead", &String::nohead)
    //     .def("notail", &String::notail)
    //     .def_static("rational", [](const double& f, const int& b, bool sign=false){ return String::rational(f, b, sign); })
    //     .def("__int__", &String::l)
    //     .def("__float__", &String::f64)
    //     .def("__str__", [](const String& self) { return self.c_str(); })
    //     .def("__repr__", [](const String& self) { return self.c_str(); })
    //     ;

    init_rtop<int>(m, "int");
    init_vec3<int>(m, "int");
    init_mat33<int>(m, "int");
    init_array2d<int>(m, "int");
    init_vec3<ftype64>(m, "double");
    init_mat33<ftype64>(m, "double");
    init_rtop<ftype64>(m, "double");
    init_mat33sym<ftype64>(m, "double");
    init_array2d<ftype64>(m, "double");
    init_matrix<ftype64>(m, "double");

    // TODO: are 32-bit versions of these actually necessary?
    init_vec3<ftype32>(m, "float");
    init_mat33<ftype32>(m, "float");
    init_rtop<ftype32>(m, "float");
    init_mat33sym<ftype32>(m, "float");
    init_array2d<ftype32>(m, "float");
    init_matrix<ftype32>(m, "float");


}
