
%include "../clipper/ccp4/ccp4_map_io.h"
namespace clipper
{
%extend CCP4MAPfile {
  %template(import_xmap_float) import_xmap<float>;
  %template(import_xmap_double) import_xmap<double>;
  %template(export_xmap_float) export_xmap<float>;
  %template(export_xmap_double) export_xmap<double>;
  %template(import_nxmap_float) import_nxmap<float>;
  %template(import_nxmap_double) import_nxmap<double>;
  %template(export_nxmap_float) export_nxmap<float>;
  %template(export_nxmap_double) export_nxmap<double>;
};
}
%include "../clipper/ccp4/ccp4_mtz_types.h"


%include "../clipper/ccp4/ccp4_mtz_io.h"

namespace clipper
{
%extend CCP4MTZfile
{
  %pythoncode %{
    open_read = log_clipper(open_read)
    import_hkl_data = log_clipper(import_hkl_data)
    
#ifdef PYTHON_PROPERTIES
    ccp4_spacegroup_number = property(ccp4_spacegroup_number)
    cell = property(cell)
    column_labels = property(column_labels)
    column_paths = property(column_paths)
    high_res_limit = property(high_res_limit)
    history = property(history)
    hkl_sampling = property(hkl_sampling)
    low_res_limit = property(low_res_limit)
    resolution = property(resolution)
    sort_order = property(sort_order)
    spacegroup = property(spacegroup)
#endif
  %}
} // extend CCP4MTZfile
} // namespace clipper
