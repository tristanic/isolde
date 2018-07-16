#include <pybind11/pybind11.h>

#include <clipper/clipper.h>

namespace py=pybind11;

py::register_exception_translator([](void *p) {
    try {
        throw;
    } catch (const clipper::Message_fatal& e) {
        PyErr_SetString(PyExc_RuntimeError(e.stream().str()));

    }
});
