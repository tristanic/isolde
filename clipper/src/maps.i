%include "../clipper/core/xmap.h"
%include "../clipper/core/nxmap.h"

//%{
  //namespace clipper {
  //typedef NXmap_base::Map_reference_base NXmap_reference_base;
  //typedef NXmap_base::Map_reference_index NXmap_reference_index;
  //typedef NXmap_base::Map_reference_coord NXmap_reference_coord;
  //typedef Xmap_base::Map_reference_base Xmap_reference_base;
  //typedef Xmap_base::Map_reference_index Xmap_reference_index;
  //typedef Xmap_base::Map_reference_coord Xmap_reference_coord;
  //}

//%}

// numpy support for Xmap and NXmap
// jon is currently editing this code piece


namespace clipper
{
%extend NXmap {
  NXmap ( Grid const & grid, RTop<ftype> const & rtop)
  {
    return new NXmap<T>( grid, rtop );
  }

  int export_numpy ( double *numpy_array, int nu, int nv, int nw, char order = 'F', std::string rot = "xyz")
  {
    std::string orders("FC");
    int oindex = orders.find(order);
    if (oindex == 2) {
      throw std::invalid_argument("Order must be either F (Fortran-style wvu) or C (C-style uvw)");
    }
    int i = 0;
    int top_u, top_v, top_w;

    Coord_grid c;
    Grid map_grid = self->grid();

    int temp;
    if (rot.compare("zyx") == 0) {
      order = orders[(oindex + 1) % 2];
      temp = nw;
      nw = nu;
      nu = temp;
    } else if (rot.compare("xyz") != 0) {
      throw std::invalid_argument ("Rotation must be either \"xyz\" or \"zyx\"!");
    }


    nu > map_grid.nu() ? top_u = map_grid.nu() : top_u = nu;
    nv > map_grid.nv() ? top_v = map_grid.nv() : top_v = nv;
    nw > map_grid.nw() ? top_w = map_grid.nw() : top_w = nw;

    if (order == 'F') {
      for ( c.w() = 0; c.w() < top_w; c.w()++ )
        for ( c.v() = 0; c.v() < top_v; c.v()++ )
          for (  c.u() = 0; c.u() < top_u; c.u()++, i++ )
            numpy_array[i] = self->get_data(c);
    } else { // order == 'C'
     for ( c.u() = 0; c.u() < top_u; c.u()++ )
        for ( c.v() = 0; c.v() < top_v; c.v()++ )
          for (  c.w() = 0; c.w() < top_w; c.w()++, i++ )
            numpy_array[i] = self->get_data(c);
    }
    return i;

  }

  int import_numpy ( double *numpy_3d_in, int nu, int nv, int nw, char order = 'F', std::string rot = "xyz" )
  {
    std::string orders("FC");
    int oindex = orders.find(order);
    if (oindex == 2) {
      throw std::invalid_argument("Order must be either F (Fortran-style wvu) or C (C-style uvw)");
    }
    int i = 0;
    int top_u, top_v, top_w;

    Coord_grid c;
    Grid map_grid = self->grid();

    int temp;
    if (rot.compare("zyx") == 0) {
      order = orders[(oindex + 1) % 2];
      temp = nw;
      nw = nu;
      nu = temp;
    } else if (rot.compare("xyz") != 0) {
      throw std::invalid_argument ("Rotation must be either \"xyz\" or \"zyx\"!");
    }

    nu > map_grid.nu() ? top_u = map_grid.nu() : top_u = nu;
    nv > map_grid.nv() ? top_v = map_grid.nv() : top_v = nv;
    nw > map_grid.nw() ? top_w = map_grid.nw() : top_w = nw;

    if (order == 'F') {
      for ( c.w() = 0; c.w() < top_w; c.w()++ )
        for ( c.v() = 0; c.v() < top_v; c.v()++ )
          for (  c.u() = 0; c.u() < top_u; c.u()++, i++ )
            self->set_data(c, numpy_3d_in[i]);
    } else { // order == 'C'
      for ( c.u() = 0; c.u() < top_u; c.u()++ )
        for ( c.v() = 0; c.v() < top_v; c.v()++ )
          for (  c.w() = 0; c.w() < top_w; c.w()++, i++ )
            self->set_data(c, numpy_3d_in[i]);
    }
    return i;
  }

  RTop<ftype> operator_orth_grid ()
  {
    return self->operator_orth_grid();
  }
  RTop<ftype> operator_grid_orth ()
  {
    return self->operator_grid_orth();
  }

} // extend NXmap




%extend Xmap {
  // TODO: add Python __init__ to add name and is_difference_map properties
  // (has to be done outside of SWIG).
  Xmap<T>::Map_reference_coord get_reference_coord(const Coord_grid& pos)
  {

    Xmap<T>::Map_reference_coord ref(*self, pos);
    return ref;
  }
  
  // reimplemented from Xmap_base
  const Grid_sampling& grid_sampling() const { return self->grid_sampling(); }
  const Cell& cell() const { return self->cell(); }
  const Spacegroup& spacegroup() const { return self->spacegroup(); }
  const Grid_range& grid_asu() const { return self->grid_asu(); }
  bool is_null() const { return self->is_null(); }
  // ~Xmap_base

  std::vector<double> _recalculate_stats ()
  {
    std::vector<double> ret;
    double max, min, mean, sd, skew, kurtosis;
    int n = 0;

    double sum = 0;
    double sum_sq = 0;
    double sum_3rd = 0;
    double sum_4th = 0;

    double rho_sq, rho_3rd, rho_4th;

    Xmap_base::Map_reference_index ix;
    ix = self->first();
    max = min = (*self)[ix];
    for (ix = self->first(); !ix.last(); ix.next()) {
      const double &rho = (*self)[ix];
      if (! clipper::Util::is_nan(rho)) {
        n++;
        if (rho < min) min = rho;
        if (rho > max) max = rho;
        rho_sq = rho*rho;
        rho_3rd = rho_sq*rho;
        rho_4th = rho_3rd*rho;

        sum += rho;
        sum_sq += rho_sq;
        sum_3rd += rho_3rd;
        sum_4th += rho_4th;
      }
    }

    if (n > 0) {
      mean = sum / n;
      double variance = sum_sq / n - mean * mean;
      sd = sqrt (variance);
      skew = sum_3rd / n - 3 * mean * variance - mean * mean * mean;
      double kt = sum_4th
                  - 4 * sum_3rd * mean
                  + 6 * sum_sq  * mean * mean
                  - 4 * sum     * mean * mean * mean
                  + mean * mean * mean * mean * n;
      kurtosis = kt / (n * variance * variance);
      ret.push_back(min);
      ret.push_back(max);
      ret.push_back(mean);
      ret.push_back(sd);
      ret.push_back(skew);
      ret.push_back(kurtosis);
    }
    else {
      throw std::out_of_range("Map has no data!");
    }
    return ret;
  }

  /*! Return list of all  grid points with multiplicity greater than 1
   *  in the asymmetric unit.
   */
  std::vector<std::vector<int> > special_positions ()
  {
    Xmap_base::Map_reference_index ix;
    std::vector<std::vector<int> > ret;
    Coord_grid this_coord;
    ix = self->first();
    size_t mult;
    for (ix = self->first(); !ix.last(); ix.next()) {
      mult = self-> multiplicity(ix.coord());
      if (mult > 1) {
        std::vector<int> this_sp;
        this_coord = ix.coord();

        this_sp.push_back(this_coord.u());
        this_sp.push_back(this_coord.v());
        this_sp.push_back(this_coord.w());
        this_sp.push_back(mult);
        ret.push_back(this_sp);
      }
    }
    return ret;
  }

  //! Return list of all grid points with multiplicity greater than 1 in the unit cell.

  std::vector<std::vector<int> > special_positions_unit_cell_grid (double frac_offset[3])
  {
    size_t count = 0;
    std::vector<std::vector<int> > ret;
    Coord_frac start_frac(frac_offset[0], frac_offset[1], frac_offset[2]);
    Grid_sampling grid = self->grid_sampling();
    Coord_grid start_coord = start_frac.coord_grid(grid);
    Coord_grid this_coord;
    size_t mult;
    bool done = false;
    for (int i = 0; i < grid.nu(); i++) {
      for (int j = 0; j < grid.nv(); j++) {
        for (int k = 0; k < grid.nw(); k++) {
          this_coord = start_coord + clipper::Coord_grid(i,j,k);
          mult = self-> multiplicity(this_coord);
          if (mult > 1) {
            std::vector<int> this_sp;
            this_sp.push_back(this_coord.u());
            this_sp.push_back(this_coord.v());
            this_sp.push_back(this_coord.w());
            this_sp.push_back(mult);
            ret.push_back(this_sp);
          }
        }
      }
    }
    return ret;
  }

  //! Return list of all (x,y,z) points with multiplicity greater than 1 in the unit cell.

  std::vector<std::vector<double> > special_positions_unit_cell_xyz (Unit_Cell uc, double frac_offset[3])
  {
    size_t count = 0;
    Cell cell = self->cell();
    Grid_sampling grid = self->grid_sampling();
    std::vector<std::vector<double> > ret;
    Coord_frac origin = uc.min().coord_frac(grid);
    Coord_frac start_frac = origin + Coord_frac(frac_offset[0], frac_offset[1], frac_offset[2]);
    Coord_grid start_coord = start_frac.coord_grid(grid);
    Coord_grid this_coord;
    size_t mult;
    bool done = false;
    for (int i = 0; i < grid.nu(); i++) {
      for (int j = 0; j < grid.nv(); j++) {
        for (int k = 0; k < grid.nw(); k++) {
          this_coord = start_coord + Coord_grid(i,j,k);
          mult = self-> multiplicity(this_coord);
          if (mult > 1) {
            std::vector<double> this_sp;
            Coord_frac this_frac = this_coord.coord_frac(grid);
            Coord_orth this_coord_orth = this_frac.coord_orth(cell);
            this_sp.push_back(this_coord_orth.x());
            this_sp.push_back(this_coord_orth.y());
            this_sp.push_back(this_coord_orth.z());
            this_sp.push_back(double(mult));
            ret.push_back(this_sp);
          }
        }
      }
    }
    return ret;
  }

  //! Export the whole asymmetric unit as a numpy array

  int export_numpy ( double *numpy_array, int nu, int nv, int nw, char order = 'F', std::string rot = "xyz" )
  {
    std::string orders("FC");
    int oindex = orders.find(order);
    if (oindex == 2) {
      throw std::invalid_argument("Order must be either F (Fortran-style wvu) or C (C-style uvw)");
    }
    int i = 0;
    int top_u, top_v, top_w;

    Coord_grid c;
    Grid map_grid = self->grid_asu();

    int temp;
    if (rot.compare("zyx") == 0) {
      order = orders[(oindex + 1) % 2];
      temp = nw;
      nw = nu;
      nu = temp;
    } else if (rot.compare("xyz") != 0) {
      throw std::invalid_argument ("Rotation must be either \"xyz\" or \"zyx\"!");
    }


    nu > map_grid.nu() ? top_u = map_grid.nu() : top_u = nu;
    nv > map_grid.nv() ? top_v = map_grid.nv() : top_v = nv;
    nw > map_grid.nw() ? top_w = map_grid.nw() : top_w = nw;

    if (order == 'F') {

      for ( c.w() = 0; c.w() < top_w; c.w()++ )
        for ( c.v() = 0; c.v() < top_v; c.v()++ )
          for ( c.u() = 0; c.u() < nu; c.u()++, i++ ) {
            if ( c.u() < map_grid.nu() && c.v() < map_grid.nv() && c.w() < map_grid.nw() )
              numpy_array[i] = self->get_data(c);
            else
              numpy_array[i] = 0.0;
          }
      return i;
    } else { // order == 'C'
      for ( c.u() = 0; c.u() < top_u; c.u()++ )
        for ( c.v() = 0; c.v() < top_v; c.v()++ )
          for ( c.w() = 0; c.w() < nw; c.w()++, i++ ) {
            if ( c.u() < map_grid.nu() && c.v() < map_grid.nv() && c.w() < map_grid.nw() )
              numpy_array[i] = self->get_data(c);
            else
              numpy_array[i] = 0.0;
          }
      return i;
    }

  }

  //! Import the whole asymmetric unit as a numpy array

  int import_numpy ( double *numpy_3d_in, int nu, int nv, int nw, char order = 'F', std::string rot = "xyz" )
  {
    std::string orders("FC");
    int oindex = orders.find(order);
    if (oindex == 2) {
      throw std::invalid_argument("Order must be either F (Fortran-style wvu) or C (C-style uvw)");
    }
    int i = 0;
    int top_u, top_v, top_w;

    Coord_grid c;
    Grid map_grid = self->grid_asu();

    int temp;
    if (rot.compare("zyx") == 0) {
      order = orders[(oindex + 1) % 2];
      temp = nw;
      nw = nu;
      nu = temp;
    } else if (rot.compare("xyz") != 0) {
      throw std::invalid_argument ("Rotation must be either \"xyz\" or \"zyx\"!");
    }

    nu > map_grid.nu() ? top_u = map_grid.nu() : top_u = nu;
    nv > map_grid.nv() ? top_v = map_grid.nv() : top_v = nv;
    nw > map_grid.nw() ? top_w = map_grid.nw() : top_w = nw;

    if (order == 'F') {
      for ( c.w() = 0; c.w() < top_w; c.w()++ )
        for ( c.v() = 0; c.v() < top_v; c.v()++ )
          for (  c.u() = 0; c.u() < top_u; c.u()++, i++ )
            self->set_data(c, numpy_3d_in[i]);
      return i;
    } else { // order == 'C'
      for ( c.u() = 0; c.u() < top_u; c.u()++ )
        for ( c.v() = 0; c.v() < top_v; c.v()++ )
          for (  c.w() = 0; c.w() < top_w; c.w()++, i++ )
            self->set_data(c, numpy_3d_in[i]);
      return i;
    }
  }

  //! Export an arbitrary box as a numpy array

  int _export_section_numpy ( double *numpy_array, int nu, int nv, int nw, Coord_grid& start, char order = 'F', std::string rot = "xyz" )
  {
    std::string orders("FC");
    int oindex = orders.find(order);
    if (oindex == 2) {
      throw std::invalid_argument("Order must be either F (Fortran-style wvu) or C (C-style uvw)");
    }
    int i = 0;
    int w, v, u;
    int maxw, maxv, maxu;
    maxv = start.v() + nv;
    //Some packages choose to store their data in zyx order rather than
    //xyz, so it's nice to provide the option to fill that way here.
    //Swapping z for x is equivalent (and more efficient, code-wise)
    //to swapping the row order and the length of the u and w arrays.
    if (rot.compare("xyz") == 0) {
      maxu = start.u() + nu;
      maxw = start.w() + nw;
    } else if (rot.compare("zyx") == 0) {
      order = orders[(oindex + 1) % 2];
      maxu = start.u() + nw;
      maxw = start.w() + nu;
    } else {
      throw std::invalid_argument ("Rotation must be either \"xyz\" or \"zyx\"!");
    }

    /*
    Coord_grid dim = end - start;
    if (dim.u() * dim.v() * dim.w() > nu*nv*nw) {
      throw std::length_error("Target array is too small to hold the requested data!");
    }
    */

    clipper::Xmap_base::Map_reference_coord ix( *self );

    if (order == 'F') {
      for ( w = start.w(); w < maxw; w++ )
        for ( v = start.v(); v < maxv; v++ )
          for ( ix.set_coord(Coord_grid(start.u(),v,w)); ix.coord().u() < maxu; ix.next_u(), i++ ) {
            numpy_array[i] = (*self)[ix];
          }
      return i;
    } else { // order == 'C'
      for ( u = start.u(); u < maxu; u++ )
        for ( v = start.v(); v < maxv; v++ )
          for ( ix.set_coord(Coord_grid(u,v,start.w())); ix.coord().w() < maxw; ix.next_w(), i++ ) {
            numpy_array[i] = (*self)[ix];
          }
      return i;
    }
  }

  //! Import an arbitrary box as a numpy array

  int import_section_numpy ( double *numpy_3d_in, int nu, int nv, int nw, Coord_grid& start, char order = 'F', std::string rot = "xyz" )
  {
    std::string orders("FC");
    int oindex = orders.find(order);
    if (oindex == 2) {
      throw std::invalid_argument("Order must be either F (Fortran-style wvu) or C (C-style uvw)");
    }
    int i = 0;
    int w, v, u;
    int maxw, maxv, maxu;
    maxv = start.v() + nv;
    //Some packages choose to store their data in zyx order rather than
    //xyz, so it's nice to provide the option to fill that way here.
    //Swapping z for x is equivalent (and more efficient, code-wise)
    //to swapping the row order and the length of the u and w arrays.
    if (rot.compare("xyz") == 0) {
      maxu = start.u() + nu;
      maxw = start.w() + nw;
    } else if (rot.compare("zyx") == 0) {
      order = orders[(oindex + 1) % 2];
      maxu = start.u() + nw;
      maxw = start.w() + nu;
    } else {
      throw std::invalid_argument ("Rotation must be either \"xyz\" or \"zyx\"!");
    }
    /*
    Coord_grid dim = end - start;
    if (dim.u() * dim.v() * dim.w() > nu*nv*nw) {
      throw std::length_error("Target array is too small to hold the requested data!");
    }
    */
    Xmap_base::Map_reference_coord ix( *self );

    if (order == 'F') {
      for ( w = start.w(); w < maxw; w++ )
        for ( v = start.v(); v < maxv; v++ )
          for ( ix.set_coord(Coord_grid(start.u(),v,w)); ix.coord().u() < maxu; ix.next_u(), i++ ) {
            (*self)[ix] = numpy_3d_in[i];
          }
      return i;
    } else { // order == 'C'
      for ( u = start.u(); u < maxu; u++ )
        for ( v = start.v(); v < maxv; v++ )
          for ( ix.set_coord(Coord_grid(u,v,start.w())); ix.coord().w() < maxw; ix.next_w(), i++ ) {
            (*self)[ix] = numpy_3d_in[i];
          }
      return i;
    }
  }

  //! Export volume interpolated onto an orthogonal grid, as a numpy array

  int export_interpolated_box_numpy ( double *numpy_array, int nx, int ny, int nz,
                                      double box_origin_xyz[3], double box_res_xyz[3],
                                      std::string mode = "cubic", char order = 'F',
                                      std::string rot = "xyz")
  {
    std::string orders("FC");
    int oindex = orders.find(order);
    if (oindex == 2) {
      throw std::invalid_argument("Order must be either F (Fortran-style wvu) or C (C-style uvw)");
    }
    int count = 0;
    double x = box_origin_xyz[0];
    double y = box_origin_xyz[1];
    double z = box_origin_xyz[2];
    Coord_orth origin(x, y, z);
    double x_inc = box_res_xyz[0];
    double y_inc = box_res_xyz[1];
    double z_inc = box_res_xyz[2];
    Coord_frac thiscoord;

    const Cell& cell = self->cell();

    if (! (!mode.compare("cubic") || !mode.compare("linear")) ) {
      throw std::invalid_argument ("Interpolator must be either cubic (default) or linear");
    }
    bool mode_cubic = mode.compare("linear");

    int temp;
    if (rot.compare("zyx") == 0) {
      order = orders[(oindex + 1) % 2];
      temp = nx;
      nx = nz;
      nz = temp;
    } else if (!(rot.compare("zyx") == 0)) {
      throw std::invalid_argument ("Rotation must be either \"xyz\" or \"zyx\"!");
    }


    if (order == 'F') {
      for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
          for (int i = 0; i < nx; i++, count++) {
            thiscoord = (origin + Coord_orth(x_inc*i, y_inc*j, z_inc*k)).coord_frac(cell);
            if (mode_cubic) {
              numpy_array[count] = self->interp<Interp_cubic>(thiscoord);
            } else {
              numpy_array[count] = self->interp<Interp_linear>(thiscoord);
            }
          }
        }
      }
    } else { // order == 'C'
      for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
          for (int k = 0; k < nz; k++, count++) {
            thiscoord = (origin + Coord_orth(x_inc*i, y_inc*j, z_inc*k)).coord_frac(cell);
            if (mode_cubic) {
              numpy_array[count] = self->interp<Interp_cubic>(thiscoord);
            } else {
              numpy_array[count] = self->interp<Interp_linear>(thiscoord);
            }
          }
        }
      }
    }
    return count;

  }



  //! Get the length of a voxel along each axis in Angstroms

  void voxel_size(double numpy_double_out[3])
  {
    Grid_sampling g = self->grid_sampling();
    Cell c = self->cell();
    numpy_double_out[0] = c.a() / g.nu();
    numpy_double_out[1] = c.b() / g.nv();
    numpy_double_out[2] = c.c() / g.nw();
  }

  //! Get the length of a voxel along each axis in fractional coordinates

  void voxel_size_frac(double numpy_double_out[3])
  {
    Grid_sampling g = self->grid_sampling();
    numpy_double_out[0] = 1.0 / g.nu();
    numpy_double_out[1] = 1.0 / g.nv();
    numpy_double_out[2] = 1.0 / g.nw();
  }

  //! Return the interpolated density value at a given fractional coordinate

  T interp_cubic_frac_coord(Coord_frac f)
  {
    return self->interp<Interp_cubic>(f);
  }

  //! Return the interpolated density value at a given (x,y,z) coordinate

  T interp_cubic_xyz(double numpy_1d_in[3])
  {
    Coord_frac thecoord = Coord_orth(numpy_1d_in[0], numpy_1d_in[1], numpy_1d_in[2]).coord_frac(self->cell());
    return self->interp<Interp_cubic>(thecoord);
  }

  //! Return the interpolated density value at a given fractional coordinate

  T interp_linear_frac_coord(Coord_frac f)
  {
    return self->interp<Interp_linear>(f);
  }

  //! Return the interpolated density value at a given (x,y,z) coordinate

  T interp_linear_xyz(double numpy_1d_in[3])
  {
    Coord_frac thecoord = Coord_orth(numpy_1d_in[0], numpy_1d_in[1], numpy_1d_in[2]).coord_frac(self->cell());
    return self->interp<Interp_linear>(thecoord);
  }



  RTop<ftype> operator_orth_grid ()
  {
    return self->operator_orth_grid();
  }
  RTop<ftype> operator_grid_orth ()
  {
    return self->operator_grid_orth();
  }

  %pythoncode %{
    __stats = None
    def recalculate_stats(self):
      self.__stats = self._recalculate_stats()
#ifdef PYTHON_PROPERTIES
    @property
#endif
    def stats(self):
      if self.__stats is None:
        self.recalculate_stats()
      return self.__stats

#ifndef PYTHON_PROPERTIES
    def max(self):
      return self.stats()[0]
    def min(self):
      return self.stats()[1]
    def mean(self):
      return self.stats()[2]
    def sigma(self):
      return self.stats()[3]
    def skewness(self):
      return self.stats()[4]
    def kurtosis(self):
      return self.stats()[5]
#else
    @property
    def max(self):
      return self.stats[0]
    @property
    def min(self):
      return self.stats[1]
    @property
    def mean(self):
      return self.stats[2]
    @property
    def sigma(self):
      return self.stats[3]
    @property
    def skewness(self):
      return self.stats[4]
    @property
    def kurtosis(self):
      return self.stats[5]
#endif

    def export_section_numpy(self, start_coord_grid, end_coord_grid = None,
                            target = None, order = 'C', rot = 'xyz'):
      '''
      Export a section of the map into a Numpy array. Required arguments are a
      starting grid coordinate and either a finishing coordinate or a 3D Numpy
      array.
      Args:
        start_coord_grid:
          Minimum corner of the box in grid coordinates (either a Clipper
          Coord_grid or an iterable of 3 ints)
        end_coord_grid (default = None; incompatible with target):
          Maximum corner of the box in grid coordinates. If set, a Numpy array
          of the required size will be automatically created and filled.
        target (default = None; incompatible with end_coord_grid):
          A pre-initialised 3D Numpy array of doubles. If set, the end_grid_coord
          will be automatically calculated.
        order (default = 'C'):
          One of 'F' or 'C'. If 'C', the array will be filled in C-style row-major
          format (the default format of Numpy arrays). If 'F' it will be filled in
          Fortran-style column-major format.
        rot (default = 'xyz'):
          One of 'xyz' or 'zyx'. Some packages choose to store their data in xyz,
          others in zyx.
      '''
      if end_coord_grid is None and target is None:
        raise TypeError('Must provide either an end grid coordinate or a target Numpy array!')
      elif end_coord_grid is not None and target is not None:
        raise TypeError('Cannot specify both an end grid coordinate and a target array!')
      if end_coord_grid is not None:
        import numpy
        # Fill from start coord to end coord inclusive
        array_size = end_coord_grid - start_coord_grid + [1,1,1]
        if type(array_size) == Coord_grid:
#ifdef PYTHON_PROPERTIES
          array_size = array_size.uvw
#else
          array_size = array_size.get_uvw()
#endif
        if rot == 'zyx':
          array_size = array_size[::-1]
        _target = numpy.empty(array_size, numpy.double)
      else:
       _target = target
      result = self._export_section_numpy(_target, start_coord_grid, order, rot)
      if target is not None:
        return result
      return _target

#ifdef PYTHON_PROPERTIES
    
    _name = None
    
    @property
    def name(self):
      return self._name
    
    @name.setter
    def name(self, name):
      self._name = name
    
    _is_difference_map = False
    
    @property
    def is_difference_map(self):
      '''
      If true, this will be treated as a difference map (with positive
      and negative contours displayed, etc.).
      '''
      return self._is_difference_map
    
    @is_difference_map.setter
    def is_difference_map(self, isdiff):
      self._is_difference_map = isdiff
    
    grid = property(grid_sampling)
    grid_sampling = property(grid_sampling)
    cell = property(cell)
    spacegroup = property(spacegroup)
    operator_grid_orth = property(operator_grid_orth)
    operator_orth_grid = property(operator_orth_grid)
    voxel_size = property(voxel_size)
    voxel_size_frac = property(voxel_size_frac)
    grid_asu = property(grid_asu)
    is_null = property(is_null)
    
    @property
    def grid_samples(self):
      return self.grid.dim
        
#endif

  %} //pythoncode

  } // extend Xmap



} // namespace clipper


//%extend Xmap<double> {

  //clipper::Xmap_base::Map_reference_coord map_reference_coord( const Coord_grid& coord )
  //{
    //clipper::Xmap_base::Map_reference_coord Map_reference_coord ret(*self, coord);
    //return ret;
  //}

  //const int& symop_of_reference_coord (const clipper::Xmap_base::Map_reference_coord& coord )
  //{
    //return coord.sym();
  //}
  ///*! Find and add the necessary translations for each symop in order
   //*  to pack the same unit cell as a given fractional coordinate. Return a
   //*  Unit_Cell object.
   //*/

  //clipper::Unit_Cell unit_cell_symops (clipper::Coord_frac ref, const clipper::Atom_list& atoms)
  //{
    //return clipper::Unit_Cell(ref, atoms, self->cell(), self->spacegroup(), self->grid_sampling());
  //}



//}


namespace clipper
{
%extend Xmap<float> {
  void fft_from (const clipper::HKL_data<clipper::data32::F_phi> &fb)
  {
    self->fft_from( fb );
  }

  void fft_to (clipper::HKL_data<clipper::data32::F_phi> &fphidata)
  {
    self->fft_to(fphidata, clipper::Xmap_base::Default);
  }

} // extend Xmap<float>
%extend Xmap<double> {
  void fft_from (const clipper::HKL_data<clipper::data64::F_phi> &fb)
  {
    self->fft_from( fb );
  }

  void fft_to (clipper::HKL_data<clipper::data32::F_phi> &fphidata ) const
  {
    self->fft_to(fphidata, clipper::Xmap_base::Default);
  }

} // extend Xmap<double>

%template(Xmap_float) Xmap<float>;
%template(Xmap_double) Xmap<double>;
%template(Xmap_int) Xmap<int>;
%template(NXmap_float) NXmap<float>;
%template(NXmap_double) NXmap<double>;
%template(NXmap_int) NXmap<int>;
} // namespace clipper
