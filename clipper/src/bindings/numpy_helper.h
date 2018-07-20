#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <vector>

namespace py=pybind11;

// check that the given Numpy array matches expected dimensions, and throw an
// error if not. direction is true for incoming, false for outgoing.
template<typename T>
void check_numpy_array_shape(py::array_t<T> target, std::vector<size_t> dim, bool direction)
{
    auto buf = target.request();
    bool fail=false;
    if (buf.ndim != dim.size())
        fail = true;
    else
        for (size_t i=0; i<buf.ndim; ++i)
            if (buf.shape[i] != dim[i])
                fail = true;
    if (fail) {
        auto shapet_txt = std::string("( ");
        auto shapeb_txt = std::string("( ");
        for (size_t i=0; i<buf.ndim; ++i)
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
