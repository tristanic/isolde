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
py::class_<RTop_frac, RTop<ftype>>(m, "RTop_frac")
    .def(py::init<>())
    .def(py::init<const RTop<>&>())
    .def(py::init<const Mat33<>&>())
    .def(py::init<const String&>())
    .def(py::init<const Mat33<>&, const Vec3<>&>())
    .def("rtop_orth", &RTop_frac::rtop_orth)
    .def("inverse", &RTop_frac::inverse)
    .def_static("identity", &RTop_frac::identity)
    .def("format_as_symop", [](const RTop_frac& self)
    {
        String s, t, xyz="xyz";
        for (int i=0; i<3; ++i) {
            t = "";
            for (int j=0; j<3; ++j)
            {
                if (self.rot()(i,j) != 0.0 ) {
                    t += (self.rot()(i,j) > 0.0) ? "+" : "-";
                    if ( Util::intr( fabs( self.rot()(i,j))) != 1)
                    {
                        t += String::rational( fabs( self.rot()(i,j)), 24);
                    }
                    t += xyz[j];
                }
            }
            if (self.trn()[i] != 0.0)
            {
                t += String::rational(self.trn()[i], 24, true);
            }
            s += t.substr((t[0] == '+') ? 1 : 0);
            if (i<2) s+=", ";
        }
        return s;
    })
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

py::class_<Isymop, RTop<int>>(m, "Isymop")
    .def(py::init<>())
    .def(py::init<const RTop<int>&>())
    .def(py::init<const Symop&, const Grid&>())
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
