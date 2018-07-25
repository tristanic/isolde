#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <vector>

namespace py=pybind11;


// check that the given Numpy array matches expected dimensions, and throw an
// error if not. direction is true for incoming, false for outgoing.
template<typename T>
void check_numpy_array_shape(py::array_t<T> target, std::vector<int> dim, bool direction)
{
    auto buf = target.request();
    bool fail=false;
    if ((size_t)buf.ndim != dim.size())
        fail = true;
    else
        for (ssize_t i=0; i<buf.ndim; ++i)
            if (buf.shape[i] != dim[i])
                fail = true;
    if (fail) {
        auto shapet_txt = std::string("( ");
        auto shapeb_txt = std::string("( ");
        for (ssize_t i=0; i<buf.ndim; ++i)
            shapet_txt += std::to_string(buf.shape[i]) + " ";
        for (size_t i=0; i<dim.size(); ++i)
            shapeb_txt += std::to_string(dim[i]) + " ";
        shapet_txt += ")";
        shapeb_txt += ")";

        std::string message;

        auto msg = std::string("Array shape mismatch! ");
        if (direction)
            msg += "Input array shape is ";
        else
            msg += "Target array shape is ";
        msg += shapet_txt;
        msg += ", while expected shape is ";
        msg += shapeb_txt;
        msg += ".";
        throw std::runtime_error(msg.c_str());
    }
}

template<class C, typename T>
py::array_t<T> array_as_numpy_1d(const C& v, int n)
{
    py::array_t<T> ret(n);
    auto buf = ret.request();
    T* ptr = (T*)buf.ptr;
    for (int i=0; i<n; ++i)
        ptr[i] = v[i];
    return ret;
}

template<class C, typename T>
void array_as_numpy_1d(const C& v, int n, py::array_t<T> target)
{
    check_numpy_array_shape(target, {n}, true);
    auto buf = target.request();
    T* ptr = (T*)buf.ptr;
    for (int i=0; i<n; ++i)
        ptr[i] = v[i];
}


template<class C, typename T>
void fill_array_from_numpy_1d(C& v, int n, py::array_t<T> arr)
{
    check_numpy_array_shape(arr, {n}, true);
    auto buf = arr.request();
    T* ptr = (T*)buf.ptr;
    for (int i=0; i<n; ++i)
        v[i] = ptr[i];
}

template<class HKLdtype, class dtype>
py::array_t<dtype> hkl_data_export_numpy(const HKLdtype& self, const int& size)
{
    auto ret = py::array_t<dtype>(size);
    dtype* ptr = (dtype*)ret.request().ptr;
    self.data_export(ptr);
    return ret;
} // hkl_data_export_numpy

template<class HKLdtype, class dtype>
void hkl_data_import_numpy(HKLdtype& self, const int& size, py::array_t<dtype> vals)
{
    check_numpy_array_shape(vals, {size}, true);
    dtype* ptr = (dtype*)vals.request().ptr;
    self.data_import(ptr);
} // hkl_data_import_numpy
