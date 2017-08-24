namespace clipper
{
%extend RTop_orth
{

  Coord_orth __mul__(const Coord_orth& c) {
    return (*self) * c;
  }

  void mat44 (double numpy_double_out[4][4])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = self->rot()(i,j);
      }
    }
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i][3] = self->trn()[i];
      numpy_double_out[3][i] = 0;
    }
    numpy_double_out[3][3] = 1;
  }
  //! Return the affine transform matrix excluding the final [0,0,0,1] row
  void mat34 (double numpy_double_out[3][4])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = self->rot()(i,j);
      }
    }
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i][3] = self->trn()[i];
    }
  }
  void rotation (double numpy_double_out[3][3])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = self->rot()(i,j);
      }
    }
  }

  void translation (double numpy_double_out[3])
  {
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i] = self->trn()[i];
    }
  }
#ifdef PYTHON_PROPERTIES
  %pythoncode %{
    mat44 = property(mat44)
    mat34 = property(mat34)
    rotation = property(rotation)
    translation = property(translation)
  %}
#endif
} // extend RTop_orth

%extend RTop_frac {

  Coord_frac __mul__(const clipper::Coord_frac& c) {
    return (*self) * c;
  }

  //! Return a full 4x4 affine transform matrix
  void mat44 (double numpy_double_out[4][4])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = self->rot()(i,j);
      }
    }
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i][3] = self->trn()[i];
      numpy_double_out[3][i] = 0;
    }
    numpy_double_out[3][3] = 1;
  }
  //! Return the affine transform matrix excluding the final [0,0,0,1] row
  void mat34 (double numpy_double_out[3][4])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = self->rot()(i,j);
      }
    }
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i][3] = self->trn()[i];
    }
  }

  void rotation (double numpy_double_out[3][3])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = self->rot()(i,j);
      }
    }
  }
  void translation (double numpy_double_out[3])
  {
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i] = self->trn()[i];
    }
  }
  // format() is only defined in the base class, so needs to be
  // re-defined here. Prints out a multi-line representation of the
  // rotation/translation matrix.
  String format() const { return self->format(); }

  // Provides a prettier representation similar to that of bona fide
  // Symop objects - e.g. (x, y+1/2, z)
  String format_as_symop() const
  {
    String s, t, xyz="xyz";
    for ( int i = 0; i < 3; i++ ) {
      t = "";
      for ( int j = 0; j < 3; j++ )
        if ( self->rot()(i,j) != 0.0 ) {
    t += ( self->rot()(i,j) > 0.0 ) ? "+" : "-";
    if ( Util::intr( fabs( self->rot()(i,j) ) ) != 1 )
      t += String::rational( fabs( self->rot()(i,j) ), 24 );
    t += xyz[j];
        }
      if ( self->trn()[i] != 0.0 )
        t += String::rational( self->trn()[i], 24, true );
      s += t.substr( ( t[0] == '+' ) ? 1 : 0 );
      if ( i < 2 ) s+= ", ";
    }
    return s;
  }
  %pythoncode %{
    def __str__(self):
#ifdef PYTHON_PROPERTIES
      return self.format_as_symop
#else
      return self.format_as_symop()
#endif
    def __hash__(self):
      return hash(self.__str__())

    def __eq__(self, other):
      return type(other) == RTop_frac and hash(self) == hash(other)
  %}

#ifdef PYTHON_PROPERTIES
  %pythoncode %{
    mat44 = property(mat44)
    mat34 = property(mat34)
    rotation = property(rotation)
    translation = property(translation)
    format = property(format)
    format_as_symop = property(format_as_symop)
  %}
#endif
} // extend RTop_frac

%extend Symop {

  Coord_frac __mul__(const clipper::Coord_frac& c) {
    return (*self) * c;
  }
  std::string __str__() { return self->format(); }

  //! Return a full 4x4 affine transform matrix
  void mat44 (double numpy_double_out[4][4])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = self->rot()(i,j);
      }
    }
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i][3] = self->trn()[i];
      numpy_double_out[3][i] = 0;
    }
    numpy_double_out[3][3] = 1;
  }
  //! Return the affine transform matrix excluding the final [0,0,0,1] row
  void mat34 (double numpy_double_out[3][4])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = self->rot()(i,j);
      }
    }
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i][3] = self->trn()[i];
    }
  }
  void rotation (double numpy_double_out[3][3])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = self->rot()(i,j);
      }
    }
  }
  void translation (double numpy_double_out[3])
  {
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i] = self->trn()[i];
    }
  }

#ifdef PYTHON_PROPERTIES
  %pythoncode %{
    mat44 = property(mat44)
    mat34 = property(mat34)
    rotation = property(rotation)
    translation = property(translation)
    format = property(format)
  %}
#endif
} // extend Symop

%extend Isymop
{
  clipper::Coord_grid __mul__(const clipper::Coord_grid& c) {
    return (*self) * c;
  }
  clipper::HKL __mul__(const clipper::HKL& c) {
    return (*self) * c;
  }
} // extend Isymop
} // namespace clipper

%template (RTop_float)  clipper::RTop<float>;
%template (RTop_double) clipper::RTop<double>;
