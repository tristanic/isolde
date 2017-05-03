%include "../clipper/core/hkl_info.h"

#ifdef PYTHON_PROPERTIES
namespace clipper
{
%extend HKL_info
{
  %pythoncode %{
    cell = property(cell)
    spacegroup = property(spacegroup)
    resolution = property(resolution)
  %}
} // extend HKL_info
} // namespace clipper
#endif

%include "../clipper/core/hkl_data.h"

namespace clipper
{

class HKL_reference_base
{
  public:
    const HKL_info& base_hkl_info() const;
    const int& index() const;
    ftype invresolsq( const HKL_data_base& hkldata ) const;
    ftype invresolsq() const;
    bool last() const;
}; // class HKL_reference_base
class HKL_reference_index : public HKL_reference_base
{
  public:
    HKL_reference_index();
    HKL_reference_index( const HKL_info& hklinfo_, const int& index );
    const HKL& hkl() const;
    const HKL_class& hkl_class() const;
    HKL_reference_index& next();
}; // class HKL_reference_index
} // namespace clipper

%{
  namespace clipper {
  typedef HKL_info::HKL_reference_base HKL_reference_base;
  typedef HKL_info::HKL_reference_index HKL_reference_index;
  }
%}

namespace clipper
{
%extend HKL {
  HKL __add__(const HKL &h2) { return *self + h2; }
  HKL __sub__(const HKL &h2) { return *self - h2; }
  HKL __neg__() { return -(*self); }
  HKL __mul__(const int& m) { return m * (*self); }
  HKL __rmul__(const int& m) { return m * (*self); }
  bool __eq__ (const HKL& h2) {return *self == h2;}
  // Transforms are handled in the definition of Isymop::__mul__()
} // extend HKL
} // namespace clipper

