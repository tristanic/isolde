#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

#include "type_conversions.h"
#include <clipper/clipper.h>

#include "numpy_helper.h"
//#include "wrap_rtop_base.h"

namespace py=pybind11;
using namespace clipper;

void init_symop(py::module& m)
{
py::class_<RTop_frac>(m, "RTop_frac")
    .def(py::init<>())
    .def(py::init<const RTop<>&>())
    .def(py::init<const Mat33<>&>())
    .def(py::init<const String&>())
    .def(py::init<const Mat33<>&, const Vec3<>&>())
    .def("rtop_orth", &RTop_frac::rtop_orth)
    .def("inverse", &RTop_frac::inverse)
    .def_static("identity", &RTop_frac::identity)
    // Base class methods
    .def("equals", [](const RTop_frac& self, const RTop_frac& other, const ftype& tol) {return self.equals(other, tol); })
    .def_property("rot",
    [](const RTop_frac& self) { return self.rot(); },
    [](RTop_frac& self, const Mat33<ftype>& rot) { self.rot() = rot; }
    )
    .def_property("trn",
         [](const RTop_frac& self) { return self.trn(); },
         [](RTop_frac& self, const Vec3<ftype>& trn) { self.trn() = trn; }
    )
    .def("is_null", [](const RTop_frac& self) { return self.is_null(); })
    .def("format", [](const RTop_frac& self) { return self.format(); })
    .def("__str__", [](const RTop_frac& self) { return self.format().c_str(); })
    ;

py::class_ <Symop, RTop_frac>(m, "Symop")
    .def(py::init<>())
    .def(py::init<const RTop<>&>())
    // From a 4x4 numpy array
    .def(py::init([](py::array_t<ftype> arr)
    {
        check_numpy_array_shape(arr, {4,4}, true);
        auto buf = arr.request();
        return std::unique_ptr<Symop>(new Symop(*reinterpret_cast<ftype(*)[4][4]>(buf.ptr)));
    }))
    .def("format", &Symop::format)
    .def("__str__", [](const Symop& self) { return self.format().c_str(); })
    ;

py::class_<Isymop>(m, "Isymop")
    .def(py::init<>())
    .def(py::init<const RTop<int>&>())
    .def(py::init<const Symop&, const Grid&>())
    // base class methods
    .def("equals", [](const Isymop& self, const Isymop& other, const int& tol) {return self.equals(other, tol); })
    .def_property("rot",
    [](const Isymop& self) { return self.rot(); },
    [](Isymop& self, const Mat33<int>& rot) { self.rot() = rot; }
    )
    .def_property("trn",
         [](const Isymop& self) { return self.trn(); },
         [](Isymop& self, const Vec3<int>& trn) { self.trn() = trn; }
    )
    .def("is_null", [](const Isymop& self) { return self.is_null(); })
    .def("format", [](const Isymop& self) { return self.format(); })
    .def("__str__", [](const Isymop& self) { return self.format().c_str(); })
    ;

py::class_<Symop_code>(m, "Symop_code")
    .def(py::init<>())
    .def(py::init<const int&>())
    .def(py::init<const Symop&>())
    .def(py::init<const Isymop&>())
    .def("init", &Symop_code::init)
    .def("code_rot", &Symop_code::code_rot)
    .def("code_trn", &Symop_code::code_trn)
    .def("symop", &Symop_code::symop)
    .def("isymop", &Symop_code::isymop)
    .def_static("identity", &Symop_code::identity)
    .def("__int__", [](const Symop_code& self) { return int(self); })
    ;
}
