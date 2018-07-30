#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/core/hkl_operators.h>
#include <clipper/core/hkl_compute.h>

#include "numpy_helper.h"

namespace py=pybind11;
using namespace clipper;

/* Hacky workaround because HKL_data_base has a protected destructor, requiring
 * us to use py::nodelete - which leads to an error if we try to instantiate
 * derived classes normally.
*/
template <typename T> struct Deleter { void operator() (T* o) const { delete o; }};

template<class C>
void catch_null(const C& c)
{
    if (c.is_null()) throw std::length_error("Array is not initialised!");
}

template<class C1, class C2>
void catch_mismatched_lengths(const C1& c1, const C2& c2)
{
    if (c1.is_null() || c1.data_size() != c2.data_size())
        throw std::out_of_range("Array sizes must match!");
}

template<class C1, class C2, class C3>
void catch_mismatched_lengths(const C1& c1, const C2& c2, const C3& c3)
{
    if (c1.is_null() || !((c1.data_size() == c2.data_size())
                       && (c2.data_size() == c3.data_size())))
        throw std::out_of_range("Array sizes must match!");
}


void declare_hkl_data_base(py::module &m)
{
    py::class_<HKL_data_base, std::unique_ptr<HKL_data_base, py::nodelete>>(m, "_HKL_data_base")
        .def_property_readonly("is_null", &HKL_data_base::is_null)
        .def_property_readonly("base_hkl_info", &HKL_data_base::base_hkl_info)
        .def_property_readonly("base_cell", &HKL_data_base::base_cell)
        .def_property_readonly("spacegroup", &HKL_data_base::spacegroup)
        .def_property_readonly("cell", &HKL_data_base::cell)
        .def_property_readonly("resolution", &HKL_data_base::resolution)
        .def_property_readonly("hkl_sampling", &HKL_data_base::hkl_sampling)
        .def_property_readonly("hkl_info", &HKL_data_base::hkl_info)
        .def("invresolsq", &HKL_data_base::invresolsq)
        .def("invresolsq_range", &HKL_data_base::invresolsq_range)
        .def_property_readonly("num_obs", &HKL_data_base::num_obs)
        .def_property_readonly("first", &HKL_data_base::first)
        .def_property_readonly("first_data", &HKL_data_base::first_data)
        .def("next_data", [](const HKL_data_base& self, HKL_info::HKL_reference_index& ih) { self.next_data(ih); })
        ;

}

// Common to all HKL datatypes
template <class C>
py::class_<HKL_data<C>> declare_HKL_data(py::module &m, const std::string &class_str,
    const std::string& docstring = "")
{
    using Class=HKL_data<C>;
    std::string pyclass_name = std::string("HKL_data_") + class_str;
    py::class_<Class, /*std::unique_ptr<Class, Deleter<Class>>,*/ HKL_data_base> theclass(m, pyclass_name.c_str(), docstring.c_str());
    theclass
        .def(py::init<>())
        .def(py::init<const HKL_info&>())
        .def(py::init<const HKL_info&, const Cell&>())
        .def(py::init<const Spacegroup&, const Cell&, const HKL_sampling&>())
        .def(py::init<const HKL_data_base&>())
        .def("init", (void (Class::*)(const HKL_info&, const Cell&)) &Class::init)
        .def("init", (void (Class::*)(const Spacegroup&, const Cell&, const HKL_sampling&)) &Class::init)
        .def("init", (void (Class::*)(const HKL_data_base&)) &Class::init)
        .def("update", &Class::update)
        .def("type", &Class::type)
        .def("missing", &Class::missing)
        .def("set_null", &Class::set_null)
        .def("data_size", &Class::data_size)
        .def("data_names", &Class::data_names)
        .def("data_export", &Class::data_export)
        .def("data_import", &Class::data_import)
        .def("mask", &Class::mask)
        .def("__getitem__", [](const Class& self, const HKL_info::HKL_reference_index& i){catch_null(self); return self[i];}, py::is_operator())
        .def("__getitem__", [](const Class& self, const HKL_info::HKL_reference_coord& ih) {
            catch_null(self);
            C data;
            if (self.get_data(ih, data))
                return data;
            throw std::out_of_range("No data equivalent to that HKL!");
        }, py::is_operator())
        .def("__setitem__", [](Class& self, const HKL_info::HKL_reference_coord& ih, const C& data) {
            if (!self.set_data(ih, data))
                throw std::out_of_range("No equivalent HKL has been indexed for this dataset!");
        }, py::is_operator())
        .def("__getitem__", [](const Class& self, const int& index) {catch_null(self); return self[index];}, py::is_operator())
        .def("__getitem__", [](const Class& self, const HKL& hkl) {
            catch_null(self);
            C data;
            if (self.get_data(hkl, data))
                return data;
            throw std::out_of_range("No data equivalent to that HKL!");
        }, py::is_operator())
        .def("__setitem__", [](Class& self, const HKL& hkl, const C& data) {
            if (!self.set_data(hkl, data))
                throw std::out_of_range("No equivalent HKL has been indexed for this dataset!");
        })
        .def("copy_from", [](Class& self, const Class& other) { self=other; })
        .def("set_all_values_to", [](Class& self, const C& value) { self=value; })
        // To/from numpy
        .def_property("data",
        [](const Class& self) -> py::tuple
        {
            auto w = self.data_size();
            auto l = self.num_obs();
            py::array_t<xtype> data({l, w});
            py::array_t<int> hkl({l,3});
            xtype* dptr = (xtype*)data.request().ptr;
            int* hptr = (int*)hkl.request().ptr;
            HKL_info::HKL_reference_index ih;
            for ( ih = self.first(); !ih.last(); ih.next())
            {
                const C& entry = self[ih];
                entry.data_export(dptr);
                dptr+= w;
                const HKL& h = ih.hkl();
                for (size_t i=0; i<3; ++i)
                    *hptr++ = h[i];
            }
            return py::make_tuple(hkl, data);
        },
        [](Class& self, py::array_t<int> hkl, py::array_t<xtype> data)
        {
            auto w = self.data_size();
            auto hbuf = hkl.request();
            auto dbuf = data.request();
            if (hbuf.shape[1] != 3 || dbuf.shape[1] != w || hbuf.shape[0] != dbuf.shape[0])
                throw std::logic_error("Array lengths don't match, or data does not have the expected width!");
            auto l = dbuf.shape[0];
            HKL_info::HKL_reference_index ih;
            for ( ih = self.first_data(); !ih.last(); self.next_data(ih))
                self.set_null(ih.index());
            int* hptr = (int*)hbuf.ptr;
            xtype* dptr = (xtype*)dbuf.ptr;
            HKL h;
            for (int i=0; i<l; ++i) {
                h = HKL(*hptr, *(hptr+1), *(hptr+2)); hptr+=3;
                self.data_import(h, dptr); dptr+= w;
            }
        },
        "Exports the data to/from numpy as a tuple of two Numpy arrays containing "
        "HKL indices and data values respectively. NOTE: setting this way will "
        "destroy all existing data in this array.")
        /*
        Since the base class HKL_data_base has a protected virtual destructor,
        exposing it to Python via PyBind11 requires it (and all derived classes)
        to be wrapped with std::unique_ptr<Class, py::nodelete>, which would
        require the python side code to explicitly handle object deletion. That
        would be a serious pain, so instead we'll hide the base class from
        Python entirely and expose all base class functions using lambdas.
        */
        // .def("is_null", [](const Class& self) { return self.is_null(); })
        // .def("base_hkl_info", [](const Class& self) { return self.base_hkl_info(); })
        // .def("base_cell", [](const Class& self) { return self.base_cell(); })
        // .def("spacegroup", [](const Class& self) { return self.spacegroup(); })
        // .def("cell", [](const Class& self) { return self.cell(); })
        // .def("resolution", [](const Class& self) { return self.resolution(); })
        // .def("hkl_sampling", [](const Class& self) { return self.hkl_sampling(); })
        // .def("hkl_info", [](const Class& self) { return self.hkl_info(); })
        // .def("invresolsq", [](const Class& self, const int& index) { return self.invresolsq(index); })
        // .def("invresolsq_range", [](const Class& self) { return self.invresolsq_range(); })
        // .def("num_obs", [](const Class& self) { return self.num_obs(); })
        // .def("first", [](const Class& self) { return self.first(); })
        // .def("first_data", [](const Class& self) { return self.first_data(); })
        // .def("next_data", [](const Class& self, HKL_info::HKL_reference_index& ih) { return self.next_data(ih); })
        // Comparison operators common to all, from hkl_operators.h
        .def(py::self & py::self)
        .def(py::self | py::self)
        .def(py::self ^ py::self)
        .def(!py::self)
        ;
    return theclass;
} //declare_HKL_data

template<class T>
void declare_hkl_data_i_sigi(py::module& m, const char* dtype)
{
    auto class_str = std::string("I_sigI_") + dtype;
    using namespace clipper::datatypes;
    using Class=HKL_data<I_sigI<T>>;
    auto pyclass = declare_HKL_data<I_sigI<T>>(m, class_str);
    pyclass
        // compute methods from hkl_compute.h
        .def("compute_scale_u_iso_isigi", [](Class& self, const T& scale, const T& u_value, const Class& isigi)
        {
            catch_null(self);
            self.compute( isigi, Compute_scale_u_iso<I_sigI<T>>(scale, u_value));
        })
        .def("compute_scale_u_aniso_isigi", [](Class& self, const T& scale, const U_aniso_orth& u_value, const Class& isigi)
        {
            catch_null(self);
            self.compute(isigi, Compute_scale_u_aniso<I_sigI<T>>(scale, u_value));
        })
        ;
} // declare_hkl_data_isigi

template<class T>
void declare_hkl_data_i_sigi_ano(py::module& m, const char* dtype)
{
    auto class_str = std::string("I_sigI_ano_") + dtype;
    using namespace clipper::datatypes;
    auto pyclass = declare_HKL_data<I_sigI_ano<T>>(m, class_str);
}

template<class T>
void declare_hkl_data_f_sigf(py::module& m, const char* dtype)
{
    auto class_str = std::string("F_sigF_") + dtype;
    using namespace clipper::datatypes;
    using Class=HKL_data<F_sigF<T>>;
    auto pyclass = declare_HKL_data<F_sigF<T>>(m, class_str);
    pyclass
        // compute methods from hkl_compute.h
        .def("compute_mean_from_fano", [](Class& self, const HKL_data<F_sigF_ano<T>>& fano)
        {
            catch_null(self);
            self.compute( fano, Compute_mean_fsigf_from_fsigfano<T>() );
        })
        .def("compute_diff_from_fano", [](Class& self, const HKL_data<F_sigF_ano<T>>& fano)
        {
            catch_null(self);
            self.compute( fano, Compute_diff_fsigf_from_fsigfano<T>() );
        })
        .def("compute_scale_u_iso_fsigf", [](Class& self, const T& scale, const T& u_value, const Class& fsigf)
        {
            catch_null(self);
            self.compute( fsigf, Compute_scale_u_iso<F_sigF<T>>(scale, u_value));
        })
        .def("compute_scale_u_aniso_fsigf", [](Class& self, const T& scale, const U_aniso_orth& u_value, const Class& fsigf)
        {
            catch_null(self);
            self.compute( fsigf, Compute_scale_u_aniso<F_sigF<T>>(scale, u_value));
        })
        ;
} // declare_hkl_data_fsigf

template<class T>
void declare_hkl_data_f_sigf_ano(py::module& m, const char* dtype)
{
    auto class_str = std::string("F_sigF_ano_") + dtype;
    using namespace clipper::datatypes;
    using Class=HKL_data<F_sigF_ano<T>>;
    auto pyclass = declare_HKL_data<F_sigF_ano<T>>(m, class_str);
    pyclass
        // compute methods from hkl_compute.h
        .def("compute_scale_u_iso_fsigfano", [](Class& self, const T& scale, const T& u_value, const Class& fsigfano)
        {
            catch_null(self);
            self.compute( fsigfano, Compute_scale_u_iso<F_sigF_ano<T>>(scale, u_value));
        })
        .def("compute_scale_u_aniso_fsigfano", [](Class& self, const T& scale, const U_aniso_orth& u_value, const Class& fsigfano)
        {
            catch_null(self);
            self.compute( fsigfano, Compute_scale_u_aniso<F_sigF_ano<T>>(scale, u_value));
        })
        ;
} //declare_hkl_data_fsigf_ano

template<class T>
void declare_hkl_data_e_sige(py::module& m, const char* dtype)
{
    auto class_str = std::string("E_sigE_") + dtype;
    using namespace clipper::datatypes;
    using Class=HKL_data<E_sigE<T>>;
    auto pyclass = declare_HKL_data<E_sigE<T>>(m, class_str);
    pyclass
        // compute methods from hkl_compute.h
        .def("compute_from_fsigf", [](Class& self, const HKL_data<F_sigF<T>>& fsigf)
        {
            catch_null(self);
            self.compute( fsigf, Compute_EsigE_from_FsigF<T>() );
        })
        // extra useful methods carried over from SWIG wrappings
        .def("scale_by_sqrt_resolution", [](Class& self, const ResolutionFn& escale)
        {
            catch_null(self);
            for (clipper::HKL_data_base::HKL_reference_index ih = self.first(); !ih.last(); ih.next())
                if ( !self[ih].missing() ) self[ih].scale( sqrt( escale.f(ih) ));
        })
        .def("scale_by_resolution", [](Class& self, const ResolutionFn& escale)
        {
            catch_null(self);
            for (clipper::HKL_data_base::HKL_reference_index ih = self.first(); !ih.last(); ih.next())
                if ( !self[ih].missing() ) self[ih].scale( escale.f(ih) );
        })
        ;
} //declare_hkl_data_e_sige

template<class T>
void declare_hkl_data_f_phi(py::module& m, const char* dtype)
{
    auto class_str = std::string("F_phi_") + dtype;
    using namespace clipper::datatypes;
    using Class=HKL_data<F_phi<T>>;
    auto pyclass = declare_HKL_data<F_phi<T>>(m, class_str);
    pyclass
        .def("compute_neg", [](Class& self, const Class& other)
        {
            //catch_mismatched_lengths(self, other);
            catch_null(self);
            self.compute( other, Compute_neg_fphi<T>());
        })
        .def("compute_add_fphi", [](Class& self, const Class& fphi1, const Class& fphi2)
        {
            catch_null(self);
            self.compute( fphi1, fphi2, Compute_add_fphi<T>());
        })
        .def("compute_sub_fphi", [](Class& self, const Class& fphi1, const Class& fphi2)
        {
            catch_null(self);
            self.compute( fphi1, fphi2, Compute_sub_fphi<T>());
        })
        .def("compute_from_fsigf_phifom", [](Class& self, const HKL_data<F_sigF<T>>& fsigf, const HKL_data<Phi_fom<T>>& phifom)
        {
            catch_null(self);
            self.compute( fsigf, phifom, Compute_fphi_from_fsigf_phifom<T>());
        })
        .def("compute_scale_u_iso_fphi", [](Class& self, const T& scale, const T& u_value, const HKL_data<F_phi<T>>& fphi)
        {
            catch_null(self);
            self.compute( fphi, Compute_scale_u_iso<F_phi<T>>(scale, u_value));
        })
        .def("compute_scale_u_aniso_fphi", [](Class& self, const T& scale, const U_aniso_orth& u_value, const HKL_data<F_phi<T>>& fphi)
        {
            catch_null(self);
            self.compute( fphi, Compute_scale_u_aniso<F_phi<T>>(scale, u_value));
        })
        // from hkl_operators.h
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * ftype32())
        .def(ftype32() * py::self)
        .def(- py::self)
        ;
} // declare_hkl_data_f_phi

template<class T>
void declare_hkl_data_phi_fom(py::module& m, const char* dtype)
{
    auto class_str = std::string("Phi_fom_") + dtype;
    using namespace clipper::datatypes;
    using Class=HKL_data<Phi_fom<T>>;
    auto pyclass = declare_HKL_data<Phi_fom<T>>(m, class_str);
    pyclass
        // compute methods from hkl_compute.h
        .def("compute_from_abcd", [](Class& self, const HKL_data<ABCD<T>>& abcd)
        {
            catch_null(self);
            self.compute( abcd, Compute_phifom_from_abcd<T>() );
        })
        ;
} // declare_hkl_data_phi_fom

template<class T>
void declare_hkl_data_abcd(py::module& m, const char* dtype)
{
    auto class_str = std::string("ABCD_") + dtype;
    using namespace clipper::datatypes;
    using Class=HKL_data<ABCD<T>>;
    auto pyclass = declare_HKL_data<ABCD<T>>(m, class_str);
    pyclass
        // compute method from hkl_compute.h
        .def("compute_from_phi_fom", [](Class& self, const HKL_data<Phi_fom<T>>& phiw)
        {
            catch_null(self);
            self.compute(phiw, Compute_abcd_from_phifom<T>());
        })
        .def("compute_add_abcd", [](Class& self, const Class& abcd1, const Class& abcd2)
        {
            catch_null(self);
            self.compute(abcd1, abcd2, Compute_add_abcd<T>());
        })
        // from hkl_operators.h
        .def(py::self + py::self);
        ;
}

template<class T>
void declare_hkl_data_d_sigd(py::module& m, const char* dtype)
{
    auto class_str = std::string("D_sigD_") + dtype;
    using namespace clipper::datatypes;
    auto pyclass = declare_HKL_data<D_sigD<T>>(m, class_str, std::string("Deprecated. Do not use."));
}


void declare_hkl_data_flag(py::module& m)
{
    using namespace clipper::datatypes;
    auto pyclass = declare_HKL_data<Flag>(m, "Flag");
    pyclass
        // from hkl_operators.h
        .def(py::self == int())
        .def(py::self != int())
        .def(py::self >= int())
        .def(py::self <= int())
        .def(py::self > int())
        .def(py::self < int())
        ;
}

void declare_hkl_data_flag_bool(py::module& m)
{
    using namespace clipper::datatypes;
    auto pyclass = declare_HKL_data<Flag_bool>(m, "Flag_bool");
}


template <class T>
void declare_CHKL_data(py::module &m, const char* class_name, const char* dtype)
{
    using Class=CHKL_data<T>;
    std::string pyclass_name = std::string("CHKL_data_") + class_name + dtype;
    py::class_<Class, /*std::unique_ptr<Class, Deleter<Class>>,*/ Container, HKL_data<T> >(m, pyclass_name.c_str())
        .def(py::init<>())
        .def(py::init<Container&, const String>())
        .def("init", (void (Class::*)(const HKL_info&, const Cell&)) &Class::init)
        .def("init", (void (Class::*)(const Spacegroup&, const Cell&, const HKL_sampling&)) &Class::init)
        .def("update", &Class::update)
        .def("copy_from", [](Class& self, const HKL_data<T>& other) { self = other; })
        .def("set_all_values_to", [](Class& self, const T& value) { self = value; });
}

void init_hkl_data(py::module& m, py::module& m32, py::module& m64)
{
    declare_hkl_data_base(m);
    {
        using namespace clipper::datatypes;
        // Non-floating-point datatypes go in the main module
        declare_hkl_data_flag(m);
        declare_hkl_data_flag_bool(m);

        declare_CHKL_data<clipper::datatypes::Flag>(m, "Flag", "");
        declare_CHKL_data<clipper::datatypes::Flag_bool>(m, "Flag_bool", "");
    }

    {
        using namespace clipper::data32;
        // 32-bit floating point datatypes go in the data32 module
        const char* suffix = "float";
        typedef ftype32 dtype;
        auto module = m32;
        declare_hkl_data_i_sigi<dtype>(module, suffix);
        declare_hkl_data_i_sigi_ano<dtype>(module, suffix);
        declare_hkl_data_f_sigf<dtype>(module, suffix);
        declare_hkl_data_f_sigf_ano<dtype>(module, suffix);
        declare_hkl_data_e_sige<dtype>(module, suffix);
        declare_hkl_data_f_phi<dtype>(module, suffix);
        declare_hkl_data_phi_fom<dtype>(module, suffix);
        declare_hkl_data_abcd<dtype>(module, suffix);
        declare_hkl_data_d_sigd<dtype>(module, suffix);

        declare_CHKL_data<I_sigI>(module, "I_sigI_", suffix);
        declare_CHKL_data<I_sigI_ano>(module, "I_sigI_ano_", suffix);
        declare_CHKL_data<F_sigF>(module, "F_sigF_", suffix);
        declare_CHKL_data<F_sigF_ano>(module, "F_sigF_ano_", suffix);
        declare_CHKL_data<E_sigE>(module, "E_sigE_", suffix);
        declare_CHKL_data<F_phi>(module, "F_phi_", suffix);
        declare_CHKL_data<Phi_fom>(module, "Phi_phom_", suffix);
        declare_CHKL_data<ABCD>(module, "ABCD_", suffix);
        declare_CHKL_data<D_sigD>(module, "D_sigD_", suffix);
    }

    {
        using namespace clipper::data64;
        // 64-bit floating point datatypes go in the data64 module
        const char* suffix = "double";
        typedef ftype64 dtype;
        auto module = m64;
        declare_hkl_data_i_sigi<dtype>(module, suffix);
        declare_hkl_data_i_sigi_ano<dtype>(module, suffix);
        declare_hkl_data_f_sigf<dtype>(module, suffix);
        declare_hkl_data_f_sigf_ano<dtype>(module, suffix);
        declare_hkl_data_e_sige<dtype>(module, suffix);
        declare_hkl_data_f_phi<dtype>(module, suffix);
        declare_hkl_data_phi_fom<dtype>(module, suffix);
        declare_hkl_data_abcd<dtype>(module, suffix);
        declare_hkl_data_d_sigd<dtype>(module, suffix);

        declare_CHKL_data<I_sigI>(module, "I_sigI_", suffix);
        declare_CHKL_data<I_sigI_ano>(module, "I_sigI_ano_", suffix);
        declare_CHKL_data<F_sigF>(module, "F_sigF_", suffix);
        declare_CHKL_data<F_sigF_ano>(module, "F_sigF_ano_", suffix);
        declare_CHKL_data<E_sigE>(module, "E_sigE_", suffix);
        declare_CHKL_data<F_phi>(module, "F_phi_", suffix);
        declare_CHKL_data<Phi_fom>(module, "Phi_phom_", suffix);
        declare_CHKL_data<ABCD>(module, "ABCD_", suffix);
        declare_CHKL_data<D_sigD>(module, "D_sigD_", suffix);
    }


}
