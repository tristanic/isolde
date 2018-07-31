
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

#include "type_conversions.h"
#include "symops.h"

#include "numpy_helper.h"


namespace py=pybind11;
using namespace clipper;


int nrows_from_shape_string(const std::string& shape)
{
    if (shape == "4x4")
        return 4;
    else if (shape == "3x4")
        return 3;
    else
        throw std::logic_error("Shape should be one of \"4x4\" or \"3x4\"!");

}
void export_rot_trn_(const Mat33<ftype>& rot, const Vec3<ftype>& trn, int nrows, ftype* ptr)
{
    const static ftype lastrow[4] = {0,0,0,1};
    for (int j=0; j<nrows; ++j)
    {
        for (int k=0; k<4; ++k)
        {
            if (j==3)
                *ptr++ = lastrow[k];
            else
            {
                if (k==3)
                    *ptr++ = trn[j];
                else
                    *ptr++ = rot(j,k);
            }

        }
    }
}

void all_matrices_frac_(const Symops& self, int nrows, py::array_t<ftype>& target)
{
    int n = self.size();
    check_numpy_array_shape(target, {n, nrows, 4}, false);
    ftype* tptr = (ftype*)target.request().ptr;
    for (int i=0; i<n; ++i)
    {
        const RTop_frac& thisop = self.with_cell_translation(i);
        const Mat33<ftype>& rot = thisop.rot();
        const Vec3<ftype>& trn = thisop.trn();
        export_rot_trn_(rot, trn, nrows, tptr);
        tptr+= nrows*4;
    }
}

void all_matrices_orth_(const Symops& self, const Cell& cell, int nrows, py::array_t<ftype>& target)
{
    int n = self.size();
    check_numpy_array_shape(target, {n, nrows, 4}, false);
    ftype* tptr = (ftype*)target.request().ptr;
    for (int i=0; i<n; ++i)
    {
        const RTop_orth& thisop = self.with_cell_translation(i).rtop_orth(cell);
        const Mat33<ftype>& rot = thisop.rot();
        const Vec3<ftype>& trn = thisop.trn();
        export_rot_trn_(rot, trn, nrows, tptr);
        tptr+=nrows*4;
    }
}


namespace py=pybind11;
using namespace clipper;

void declare_symops(py::module& m)
{
    py::class_<Symops>(m, "Symops")
        .def(py::init<>())
        .def(py::init<std::vector<RTop_frac>&>())
        .def(py::init<std::vector<Symop>&>())
        .def("__getitem__", &Symops::at)
        .def("__setitem__", &Symops::replace)
        .def("with_cell_translation", &Symops::with_cell_translation)
        .def("append", (void (Symops::*)(const Symop& op, const Coord_frac& offset)) &Symops::append)
        .def("append", (void (Symops::*)(const RTop_frac& op)) &Symops::append)
        .def("pop", &Symops::pop)
        .def("__len__", &Symops::size)
        .def("all_matrices_frac",
            [](const Symops& self, const std::string& shape) -> py::array_t<ftype>
            {
                auto nrows = nrows_from_shape_string(shape);
                py::array_t<ftype> target({(int)self.size(), nrows, 4});
                all_matrices_frac_(self, nrows, target);
                return target;
            })
        .def("all_matrices_frac",
            [](const Symops& self, const std::string& shape, py::array_t<ftype> target) -> void
            {
                auto nrows = nrows_from_shape_string(shape);
                all_matrices_frac_(self, nrows, target);
            })
            .def("all_matrices_orth",
                [](const Symops& self, const Cell& cell, const std::string& shape) -> py::array_t<ftype>
                {
                    auto nrows = nrows_from_shape_string(shape);
                    py::array_t<ftype> target({(int)self.size(), nrows, 4});
                    all_matrices_orth_(self, cell, nrows, target);
                    return target;
                })
            .def("all_matrices_orth",
                [](const Symops& self, const Cell& cell, const std::string& shape, py::array_t<ftype> target) -> void
                {
                    auto nrows = nrows_from_shape_string(shape);
                    all_matrices_orth_(self, cell, nrows, target);
                })
        ;
}

void init_symops(py::module& m)
{
    declare_symops(m);
}
