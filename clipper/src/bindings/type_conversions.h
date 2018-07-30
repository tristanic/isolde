#pragma once
#include <pybind11/pybind11.h>

#include <clipper/clipper.h>

namespace pybind11 { namespace detail {
    template <> struct type_caster<clipper::String> {
    public:
        PYBIND11_TYPE_CASTER(clipper::String, _("String"));

        bool load(handle src, bool) {
            PyObject * source = src.ptr();
            PyObject *tmp = PyUnicode_AsUTF8String(source);
            if (!tmp)
                return false;
            char *buffer;
            ssize_t length;
            if (PYBIND11_BYTES_AS_STRING_AND_SIZE(tmp, &buffer, &length))
                return false;
            value = clipper::String(std::string(buffer, (size_t) length));
            Py_DECREF(tmp);
            return true;
        }

        static handle cast(const clipper::String& src, return_value_policy /* policy */, handle /* parent*/) {
            return PyUnicode_FromString(src.c_str());
        }

    }; //type_caster<clipper::String>

}} // namespace pybind11::detail
