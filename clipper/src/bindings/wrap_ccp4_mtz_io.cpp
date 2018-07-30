#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "type_conversions.h"
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>

namespace py=pybind11;
using namespace clipper;

void declare_ccp4mtzfile(py::module& m)
{
    py::class_<CCP4MTZfile> ccp4mtz(m, "CCP4MTZfile");
    ccp4mtz
        .def(py::init<>())
        .def("open_read", &CCP4MTZfile::open_read)
        // .def("open_read", [](CCP4MTZfile& self, const std::string& fname) { self.open_read(String(fname)); })
        .def("close_read", &CCP4MTZfile::close_read)
        .def("open_append", &CCP4MTZfile::open_append)
        // .def("open_append", [](CCP4MTZfile& self, const std::string& fname) { self.open_append(String(fname)); })
        .def("close_append", &CCP4MTZfile::close_append)
        .def("open_write", &CCP4MTZfile::open_write)
        // .def("open_write", [](CCP4MTZfile& self, const std::string& fname) { self.open_write(String(fname)); })
        .def("close_write", &CCP4MTZfile::close_write)
        .def_property_readonly("spacegroup", &CCP4MTZfile::spacegroup)
        .def_property_readonly("cell", &CCP4MTZfile::cell)
        .def_property_readonly("resolution", &CCP4MTZfile::resolution)
        .def_property_readonly("hkl_sampling", &CCP4MTZfile::hkl_sampling)

        .def("import_hkl_list", &CCP4MTZfile::import_hkl_list)
        .def("import_hkl_info", &CCP4MTZfile::import_hkl_info)
        .def("import_hkl_info", [](CCP4MTZfile& self, HKL_info& target) { self.import_hkl_info(target, true); })
        .def("import_crystal", &CCP4MTZfile::import_crystal)
        .def("import_dataset", &CCP4MTZfile::import_dataset)
        .def("import_hkl_data", (void (CCP4MTZfile::*)(HKL_data_base&, const String)) &CCP4MTZfile::import_hkl_data)
        // .def("import_hkl_data", [](CCP4MTZfile& self, HKL_data_base& cdata, const std::string& mtzpath) { self.import_hkl_data(mtzpath); })

        .def("export_hkl_info", &CCP4MTZfile::export_hkl_info)
        .def("export_crystal", &CCP4MTZfile::export_crystal)
        .def("export_dataset", &CCP4MTZfile::export_dataset)
        .def("export_hkl_data", (void (CCP4MTZfile::*)(const HKL_data_base&, const String)) &CCP4MTZfile::export_hkl_data)

        .def("import_chkl_data", &CCP4MTZfile::import_chkl_data)
        .def("import_chkl_data", [](CCP4MTZfile& self, Container& target, const String& mtzpath) { self.import_chkl_data(target, mtzpath, ""); })
        .def("export_chkl_data", &CCP4MTZfile::export_chkl_data)

        .def_property_readonly("column_paths", &CCP4MTZfile::column_paths)
        .def_property_readonly("assigned_paths", &CCP4MTZfile::assigned_paths)

        .def_property("title",
        [] (const CCP4MTZfile& self) { return self.title().c_str(); },
        [] (CCP4MTZfile& self, const std::string& name) { self.set_title(String(name)); }
        )

        .def_property("history",
        [] (const CCP4MTZfile& self) //-> std::vector<std::string>>
        {
            auto history = self.history();
            std::vector<std::string> hstr;
            for (const auto &h: history)
                hstr.push_back(static_cast<const std::string&>(h));
            return hstr;
        },
        [] (CCP4MTZfile& self, const std::vector<std::string>& hstr)
        {
            std::vector<String> history;
            for (const auto& h: hstr)
                history.push_back(String(h));
            self.set_history(history);
        })
        .def_property_readonly("num_reflections", &CCP4MTZfile::num_reflections)
        .def_property_readonly("sort_order", &CCP4MTZfile::sort_order)
        .def_property_readonly("low_res_limit", &CCP4MTZfile::low_res_limit)
        .def_property_readonly("high_res_limit", &CCP4MTZfile::high_res_limit)
        .def_property_readonly("ccp4_spacegroup_number", &CCP4MTZfile::ccp4_spacegroup_number)
        .def_property("spacegroup_confidence",
            &CCP4MTZfile::spacegroup_confidence,
            &CCP4MTZfile::set_spacegroup_confidence
        )
        .def("set_verbose", &CCP4MTZfile::set_verbose)
        ;
} // declare_ccp4mtzfile





void init_ccp4_mtz_io(py::module& m)
{
    declare_ccp4mtzfile(m);
} // init_ccp4_mtz_io
