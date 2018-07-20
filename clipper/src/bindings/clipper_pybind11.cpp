#include <pybind11/pybind11.h>

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-cif.h>
#include <clipper/clipper-cns.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-phs.h>

#include "wrap_hkl_datatypes.h"
#include "wrap_container_hkl.h"

namespace py=pybind11;
void init_hkl_info(py::module &m);
//void init_hkl_data(py::module &m);
void init_containers(py::module &m);
void init_clipper_types(py::module &m, py::module &m32, py::module &m64);

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
    //init_hkl_data(m);
    init_clipper_types(m, m32, m64);


    declare_HKL_data<clipper::datatypes::Flag>(m, "Flag");
    declare_HKL_data<clipper::datatypes::Flag_bool>(m, "Flag_bool");
    declare_CHKL_data<clipper::datatypes::Flag>(m, "Flag");
    declare_CHKL_data<clipper::datatypes::Flag_bool>(m, "Flag_bool");

    {
        using namespace clipper::data32;
        declare_HKL_data<I_sigI>(m32, "I_sigI_float");
        declare_HKL_data<I_sigI_ano>(m32, "I_sigI_anom");
        declare_HKL_data<F_sigF>(m32, "F_sigF");
        declare_HKL_data<F_sigF_ano>(m32, "F_sigF_anom");
        declare_HKL_data<E_sigE>(m32, "E_sigE");
        declare_HKL_data<F_phi>(m32, "F_phi");
        declare_HKL_data<Phi_fom>(m32, "Phi_phom");
        declare_HKL_data<ABCD>(m32, "ABCD");
        declare_HKL_data<D_sigD>(m32, "D_sigD");

        declare_CHKL_data<I_sigI>(m32, "I_sigI_float");
        declare_CHKL_data<I_sigI_ano>(m32, "I_sigI_anom");
        declare_CHKL_data<F_sigF>(m32, "F_sigF");
        declare_CHKL_data<F_sigF_ano>(m32, "F_sigF_anom");
        declare_CHKL_data<E_sigE>(m32, "E_sigE");
        declare_CHKL_data<F_phi>(m32, "F_phi");
        declare_CHKL_data<Phi_fom>(m32, "Phi_phom");
        declare_CHKL_data<ABCD>(m32, "ABCD");
        declare_CHKL_data<D_sigD>(m32, "D_sigD");
    }
    {
        using namespace clipper::data64;
        declare_HKL_data<I_sigI>(m64, "I_sigI");
        declare_HKL_data<I_sigI_ano>(m64, "I_sigI_anom");
        declare_HKL_data<F_sigF>(m64, "F_sigF");
        declare_HKL_data<F_sigF_ano>(m64, "F_sigF_anom");
        declare_HKL_data<E_sigE>(m64, "E_sigE");
        declare_HKL_data<F_phi>(m64, "F_phi");
        declare_HKL_data<Phi_fom>(m64, "Phi_phom");
        declare_HKL_data<ABCD>(m64, "ABCD");
        declare_HKL_data<D_sigD>(m64, "D_sigD");

        declare_CHKL_data<I_sigI>(m64, "I_sigI_float");
        declare_CHKL_data<I_sigI_ano>(m64, "I_sigI_anom");
        declare_CHKL_data<F_sigF>(m64, "F_sigF");
        declare_CHKL_data<F_sigF_ano>(m64, "F_sigF_anom");
        declare_CHKL_data<E_sigE>(m64, "E_sigE");
        declare_CHKL_data<F_phi>(m64, "F_phi");
        declare_CHKL_data<Phi_fom>(m64, "Phi_phom");
        declare_CHKL_data<ABCD>(m64, "ABCD");
        declare_CHKL_data<D_sigD>(m64, "D_sigD");
    }

}
