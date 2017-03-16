/* Cell::descr() appears to be the only line in cell.h that SWIG balks at -
 * rather than properly handle the implicit cast of Cell back to
 * Cell_descr, it appears to attempt to define a whole new Cell_descr
 * object and then crashes because Cell_descr already exists. Tell
 * SWIG to ignore it, and we can just include cell.h. For Python
 * purposes it's much more sensible to return the dimensions etc. as
 * numpy arrays rather than requiring a separate function call for each
 * side length, angle etc. - so let's add ignore directives for each
 * individual call and add the array calls below.
 */
namespace clipper
{
//%ignore Cell::matrix_orth() const;
//%ignore Cell::matrix_frac() const;
%ignore Cell::descr() const;

#ifdef PYTHON_PROPERTIES
%ignore Cell_descr::a() const;
%ignore Cell_descr::b() const;
%ignore Cell_descr::c() const;
%ignore Cell_descr::alpha() const;
%ignore Cell_descr::beta() const;
%ignore Cell_descr::gamma() const;
%ignore Cell_descr::alpha_deg() const;
%ignore Cell_descr::beta_deg() const;
%ignore Cell_descr::gamma_deg() const;
%ignore Cell::a_star() const;
%ignore Cell::b_star() const;
%ignore Cell::c_star() const;
%ignore Cell::alpha_star() const;
%ignore Cell::beta_star() const;
%ignore Cell::gamma_star() const;
%ignore Cell::is_null() const;
%ignore Cell::equals(const Cell& other, const ftype tol=1.0) const;
#endif

// Redirect format() to Python __str__
%rename(__str__) Cell::format() const;
}




%include "../clipper/core/cell.h"

namespace clipper
{
%extend Cell_descr {
  void dim(double numpy_double_out[3])
  {
    numpy_double_out[0] = self->a();
    numpy_double_out[1] = self->b();
    numpy_double_out[2] = self->c();
  }
  void angles(double numpy_double_out[3])
  {
    numpy_double_out[0] = self->alpha();
    numpy_double_out[1] = self->beta();
    numpy_double_out[2] = self->gamma();
  }
  void angles_deg(double numpy_double_out[3])
  {
    numpy_double_out[0] = self->alpha_deg();
    numpy_double_out[1] = self->beta_deg();
    numpy_double_out[2] = self->gamma_deg();
  }
#ifdef PYTHON_PROPERTIES
  %pythoncode %{
    dim = property(dim)
    angles = property(angles)
    angles_deg = property(angles_deg)
  %}
#endif

} // extend Cell_descr

%extend Cell {
  //Mat33<float> matrix_orth()
  //{
    //Mat33<float> orth;
    //orth(0,0) = (self->matrix_orth())(0,0);
    //orth(0,1) = (self->matrix_orth())(0,1);
    //orth(0,2) = (self->matrix_orth())(0,2);
    //orth(1,0) = (self->matrix_orth())(1,0);
    //orth(1,1) = (self->matrix_orth())(1,1);
    //orth(1,2) = (self->matrix_orth())(1,2);
    //orth(2,0) = (self->matrix_orth())(2,0);
    //orth(2,1) = (self->matrix_orth())(2,1);
    //orth(2,2) = (self->matrix_orth())(2,2);
    //return orth;
  //};
  //Mat33<float> matrix_frac()
  //{
    //Mat33<float> frac;
    //frac(0,0) = (self->matrix_frac())(0,0);
    //frac(0,1) = (self->matrix_frac())(0,1);
    //frac(0,2) = (self->matrix_frac())(0,2);
    //frac(1,0) = (self->matrix_frac())(1,0);
    //frac(1,1) = (self->matrix_frac())(1,1);
    //frac(1,2) = (self->matrix_frac())(1,2);
    //frac(2,0) = (self->matrix_frac())(2,0);
    //frac(2,1) = (self->matrix_frac())(2,1);
    //frac(2,2) = (self->matrix_frac())(2,2);
    //return frac;
  //};
  void dim(double numpy_double_out[3])
  {
    numpy_double_out[0] = self->a();
    numpy_double_out[1] = self->b();
    numpy_double_out[2] = self->c();
  }
  void recip_dim(double numpy_double_out[3])
  {
    numpy_double_out[0] = self->a_star();
    numpy_double_out[1] = self->b_star();
    numpy_double_out[2] = self->c_star();
  }

  void angles(double numpy_double_out[3])
  {
    numpy_double_out[0] = self->alpha();
    numpy_double_out[1] = self->beta();
    numpy_double_out[2] = self->gamma();
  }
  void angles_deg(double numpy_double_out[3])
  {
    numpy_double_out[0] = self->alpha_deg();
    numpy_double_out[1] = self->beta_deg();
    numpy_double_out[2] = self->gamma_deg();
  }
  void recip_angles(double numpy_double_out[3])
  {
    numpy_double_out[0] = self->alpha_star();
    numpy_double_out[1] = self->beta_star();
    numpy_double_out[2] = self->gamma_star();
  }
  void recip_angles_deg(double numpy_double_out[3])
  {
    numpy_double_out[0] = Util::rad2d(self->alpha_star());
    numpy_double_out[1] = Util::rad2d(self->beta_star());
    numpy_double_out[2] = Util::rad2d(self->gamma_star());
  }
  // Functional equivalent to original Cell::descr
  Cell_descr cell_descr()
  {
    return Cell_descr(self->a(), self->b(), self->c(),
                      self->alpha(), self->beta(), self->gamma());
  }
  bool __eq__ (const Cell& other) const {
      return self->equals(other);
  }

#ifdef PYTHON_PROPERTIES
  %pythoncode %{
    matrix_orth = property(matrix_orth)
    matrix_frac = property(matrix_frac)
    dim = property(dim)
    recip_dim = property(recip_dim)
    angles = property(angles)
    angles_deg = property(angles_deg)
    recip_angles = property(recip_angles)
    recip_angles_deg = property(recip_angles)
    metric_real = property(metric_real)
    metric_reci = property(metric_reci)
    volume = property(volume)
    cell_descr = property(cell_descr)
  %}
#endif

}; // extend Cell

} // namespace clipper
