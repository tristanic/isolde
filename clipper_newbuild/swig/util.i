%feature("docstring") clipper::Util::get_minmax_grid "
 Find the minimum and maximum grid coordinates of a box encompassing the
 coordinates in numpy_2d_in given the cell and grid sampling, and return a
 numpy array [min, max]."

namespace clipper
{
%extend Util {
  //! Some useful extensions to the Util class
  /*! Find the minimum and maximum grid coordinates of a box encompassing
  * the coordinates in numpy_2d_in given the cell and grid sampling,
  * and return a numpy array [min, max].
  */
  static void get_minmax_grid(int numpy_int_out[2][3],
    double* numpy_2d_in, int n1, int n2,
    const clipper::Cell& cell, const clipper::Grid_sampling& grid)
  {
    if (n2 != 3) {
      throw std::out_of_range("Input should be an array of 3D coordinates!");
    }
    Coord_grid ref_min =
      (Coord_orth(numpy_2d_in[0], numpy_2d_in[1], numpy_2d_in[2]))
                            .coord_frac(cell).coord_grid(grid);
    Coord_grid ref_max = ref_min;
    for (size_t i = 0; i < n1*n2; i+=n2 ) {
      Coord_grid thiscoord =
        Coord_orth(numpy_2d_in[i], numpy_2d_in[i+1], numpy_2d_in[i+2])
                            .coord_frac(cell).coord_grid(grid);
      for (size_t j = 0; j < 3; j++) {
        if (thiscoord[j] < ref_min[j]) ref_min[j] = thiscoord[j];
        else if (thiscoord[j] > ref_max[j]) ref_max[j] = thiscoord[j];
      }
    }
    for (size_t i = 0; i < 3; i++) {
      numpy_int_out[0][i] = ref_min[i];
      numpy_int_out[1][i] = ref_max[i];
    }
  }

} // extend Util
} // namespace clipper

