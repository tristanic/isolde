#include <pybind11/pybind11.h>

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-cif.h>
#include <clipper/clipper-cns.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-phs.h>

namespace py=pybind11;
void init_hkl_info(py::module &m);
void init_hkl_data(py::module &m, py::module& m32, py::module& m64);
void init_hkl_datatypes(py::module &m, py::module &m32, py::module &m64);
void init_containers(py::module &m);
void init_clipper_types(py::module &m, py::module &m32, py::module &m64);
void init_coords(py::module &m);
void init_derivs(py::module &m);
void init_symop(py::module& m);
void init_symops(py::module& m);
void init_cell(py::module& m);
void init_unit_cell(py::module& m);
void init_spacegroup(py::module& m);
void init_clipper_stats(py::module& m);
void init_nxmap(py::module& m);
void init_xmap(py::module& m);
void init_map_utils(py::module& m);
void init_clipper_util(py::module& m);
void init_atomsf(py::module& m);


void init_ccp4_mtz_io(py::module& m);


using namespace clipper;


PYBIND11_MODULE(clipper_python, m) {
    m.doc() = "Python wrapper for the Clipper crystallographic library.";

    py::register_exception_translator([](std::exception_ptr p) {
        try {
            if (p) std::rethrow_exception(p);
        } catch (const clipper::Message_fatal& e) {
            PyErr_SetString(PyExc_RuntimeError, e.text().c_str());
        }
    });

    py::module m32 = m.def_submodule("data32", "32-bit data types");
    py::module m64 = m.def_submodule("data64", "64-bit data types");

    init_hkl_info(m);
    init_containers(m);
    init_hkl_data(m, m32, m64);
    init_hkl_datatypes(m, m32, m64);
    init_clipper_types(m, m32, m64);
    init_coords(m);
    init_derivs(m);
    init_symop(m);
    init_symops(m);
    init_cell(m);
    init_unit_cell(m);
    init_spacegroup(m);
    init_clipper_stats(m);

    init_nxmap(m);
    init_xmap(m);
    init_map_utils(m);

    init_atomsf(m);
    init_clipper_util(m);

    init_ccp4_mtz_io(m);


}
