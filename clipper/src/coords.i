namespace clipper {
  /* Strangely, only Coord_grid has both getter and setter accessors
   * for u, v and w. The remainder of the coordinate classes are read-
   * only for us here in SWIG. If we want to implement a setter function
   * to give the appearance of replacing values in-place, it'll have to
   * be done in Python. All the accessors for individual values are
   * renamed to hidden methods here, and a new function is defined below
   * for each class to return (x,y,z) or (u,v,w) as a numpy array.
   */
  %rename("_u_ptr") Coord_grid::u();
  %rename("_v_ptr") Coord_grid::v();
  %rename("_w_ptr") Coord_grid::w();
#ifdef PYTHON_PROPERTIES
  %rename("_u") Coord_grid::u() const;
  %rename("_v") Coord_grid::v() const;
  %rename("_w") Coord_grid::w() const;
  %rename("_u") Coord_frac::u() const;
  %rename("_v") Coord_frac::v() const;
  %rename("_w") Coord_frac::w() const;
  %rename("_x") Coord_orth::x() const;
  %rename("_y") Coord_orth::y() const;
  %rename("_z") Coord_orth::z() const;
  %rename("_u") Coord_map::u() const;
  %rename("_v") Coord_map::v() const;
  %rename("_w") Coord_map::w() const;
  %rename("_nu") Grid_sampling::nu() const;
  %rename("_nv") Grid_sampling::nv() const;
  %rename("_nw") Grid_sampling::nu() const;
  %rename("_get_element") Atom::element() const;
  %rename("_set_element") Atom::set_element;
  %rename("_get_occupancy") Atom::occupancy() const;
  %rename("_set_occupancy") Atom::set_occupancy;
  %rename("_get_u_iso") Atom::u_iso() const;
  %rename("_set_u_iso") Atom::set_u_iso;
  %rename("_get_u_aniso") Atom::u_aniso_orth() const;
  %rename("_set_u_aniso") Atom::set_u_aniso_orth;
  %rename("_get_coord_orth") Atom::coord_orth() const;
  %rename("_set_coord_orth") Atom::set_coord_orth;
#endif
} // namespace clipper

%pythoncode %{
# The complete list of scatterers found in clipper/core/atomsf.cpp.
# Clipper itself doesn't check incoming atom names for legality, but
# it's easy to do so here in Python. Used by setter methods in Atom and
# Atom_list objects.
ATOM_NAMES = set(
     ['H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',
      'Ne', 'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar',
      'K',  'Ca', 'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co',
      'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
      'Rb', 'Sr', 'Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
      'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',  'Xe',
      'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',
      'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
      'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl',
      'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
      'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
      'H1-',  'Li1+', 'Be2+', 'Cval', 'O1-',  'O2-',  'F1-',
      'Na1+', 'Mg2+', 'Al3+', 'Siva', 'Si4+', 'Cl1-', 'K1+',
      'Ca2+', 'Sc3+', 'Ti2+', 'Ti3+', 'Ti4+', 'V2+',  'V3+',
      'V5+',  'Cr2+', 'Cr3+', 'Mn2+', 'Mn3+', 'Mn4+', 'Fe2+',
      'Fe3+', 'Co2+', 'Co3+', 'Ni2+', 'Ni3+', 'Cu1+', 'Cu2+',
      'Zn2+', 'Ga3+', 'Ge4+', 'Br1-', 'Rb1+', 'Sr2+', 'Y3+',
      'Zr4+', 'Nb3+', 'Nb5+', 'Mo3+', 'Mo5+', 'Mo6+', 'Ru3+',
      'Ru4+', 'Rh3+', 'Rh4+', 'Pd2+', 'Pd4+', 'Ag1+', 'Ag2+',
      'Cd2+', 'In3+', 'Sn2+', 'Sn4+', 'Sb3+', 'Sb5+', 'I1-',
      'Cs1+', 'Ba2+', 'La3+', 'Ce3+', 'Ce4+', 'Pr3+', 'Pr4+',
      'Nd3+', 'Pm3+', 'Sm3+', 'Eu2+', 'Eu3+', 'Gd3+', 'Tb3+',
      'Dy3+', 'Ho3+', 'Er3+', 'Tm3+', 'Yb2+', 'Yb3+', 'Lu3+',
      'Hf4+', 'Ta5+', 'W6+',  'Os4+', 'Ir3+', 'Ir4+', 'Pt2+',
      'Pt4+', 'Au1+', 'Au3+', 'Hg1+', 'Hg2+', 'Tl1+', 'Tl3+',
      'Pb2+', 'Pb4+', 'Bi3+', 'Bi5+', 'Ra2+', 'Ac3+', 'Th4+',
      'U3+',  'U4+',  'U6+',  'Np3+', 'Np4+', 'Np6+', 'Pu3+',
      'Pu4+', 'Pu6+'])

  
  
%}

%include "../clipper/core/coords.h"


namespace clipper
{
%extend Coord_grid {
  // Because conversion of Numpy int scalars is a nightmare, let's add
  // an extra constructor that takes an int array directly.
  Coord_grid(long v[3]) {
    Coord_grid *newCoord = new Coord_grid(v[0], v[1], v[2]);
      return newCoord;
  }
  Coord_grid(const long& u, const long& v, const long& w) {
    Coord_grid *newCoord = new Coord_grid(u, v, w);
    return newCoord;
  }

  bool __ne__(const Coord_grid& g2) { return *self != g2; }
  bool __eq__(const Coord_grid& g2) { return *self == g2; }
  Coord_grid __add_base__(const Coord_grid& g2) { return *self + g2; }
  Coord_grid __neg__() { return -(*self); }
  Coord_grid __sub_base__(const Coord_grid& g2) { return *self - g2; }
  Coord_grid __mul__(const int& m) { return m * (*self); }
  Coord_grid __rmul__(const int& m) { return m * (*self); }
  // Transform operator handled in Isymop::__mul__()
  std::string __str__() {return self->format();}
  std::string __repr__() {return "Coord_grid " + self->format();}
  //! Allow addition and subtraction with any Python iterable of 3 numbers
  %pythoncode %{
    @staticmethod
    @safesplat_int
    def _new_coord_grid(u, v, w):
        return Coord_grid(u, v, w)
    
    def __add__(self, other):
      '''
      This __add__ function should allow you to add any iterable of three
      ints to your Coord_grid object, as long as your sum is written 
            Coord_grid + other_iterable
      Note that the reverse:
            other_iterable + Coord_grid
      ... will not work if other_iterable is a Numpy array, since the
      numpy.ndarray __add__() function takes precedence and attempts to
      do:
        [other_iterable[0]+Coord_grid, other_iterable[1]+Coord_grid, ...]
      For this reason, it is probably best to leave the __radd__ function
      unimplemented.
      '''
      try:
        return self.__add_base__(other)
      except:
        return self + self._new_coord_grid(other)

    def __sub__(self, other):
      try:
        return self.__sub_base__(other)
      except:
        return self - self._new_coord_grid(other)

  %}
#ifdef PYTHON_PROPERTIES
  void _get_uvw(int numpy_int_out[3])
#else
  void get_uvw(int numpy_int_out[3])
#endif
  {
    for (int i = 0; i < 3; i++) {
      numpy_int_out[i] = ((*self)[i]);
    }
  }
#ifdef PYTHON_PROPERTIES
  void _set_uvw(long v[3])
#else
  void set_uvw(long v[3])
#endif
  {
    self->u() = v[0];
    self->v() = v[1];
    self->w() = v[2];
  }
#ifdef PYTHON_PROPERTIES
  %pythoncode %{
    uvw = property(_get_uvw, _set_uvw)
  %}
#endif

} // extend Coord_grid

%extend Coord_orth
{
  Coord_orth __add_base__(const Coord_orth &h2) { return *self + h2; }
  Coord_orth __sub_base__(const Coord_orth &h2) { return *self - h2; }
  Coord_orth __neg__() { return -(*self); }
  Coord_orth __mul__ ( const double &f ) { return f * (*self); }
  Coord_orth __rmul__ ( const double &f ) { return f * (*self); }
  // Transforms handled in RTop_orth::__mul__()
  std::string __str__() {return self->format();}
  std::string __repr__() {return "Coord_orth " + self->format();}
  //! Allow addition and subtraction with any Python iterable of 3 numbers
  %pythoncode %{
    @staticmethod
    @safesplat_float
    def _new_coord_orth(x, y, z):
        return Coord_orth(x, y, z)
        
    def __add__(self, other):
      '''
      This __add__ function should allow you to add any iterable of three
      numbers to your Coord_orth object, as long as your sum is written 
            Coord_orth + other_iterable
      Note that the reverse:
            other_iterable + Coord_orth
      ... will not work if other_iterable is a Numpy array, since the
      numpy.ndarray __add__() function takes precedence and attempts to
      do:
        [other_iterable[0]+Coord_orth, other_iterable[1]+Coord_orth, ...]
      For this reason, it is probably best to leave the __radd__ and 
      __rsub__ functions unimplemented.
      '''
      try:
        return self.__add_base__(other)
      except:
        return self + self._new_coord_orth(other)

    def __sub__(self, other):
      try:
        return self.__sub_base__(other)
      except:
        return self - self._new_coord_orth(other)

  %}


#ifdef PYTHON_PROPERTIES
  void _get_xyz(double numpy_double_out[3])
#else
  void get_xyz(double numpy_double_out[3])
#endif
  {
    for (int i = 0; i < 3; i++) {
      numpy_double_out[i] = (*self)[i];
    }
  }
#ifdef PYTHON_PROPERTIES
  %pythoncode %{
    xyz = property(_get_xyz)
  %}
#endif

} // extend Coord_orth

%extend Coord_frac
{
  Coord_frac __add_base__(const Coord_frac &h2) { return *self + h2; }
  Coord_frac __sub_base__(const Coord_frac &h2) { return *self - h2; }
  Coord_frac __neg__() { return -(*self); }
  Coord_frac __mul__ ( const double &f ) { return f * (*self); }
  Coord_frac __rmul__ ( const double &f ) { return f * (*self); }
  // Transforms handled in RTop_frac::__mul__()
  std::string __str__() {return self->format();}
  std::string __repr__() {return "Coord_frac " + self->format();}
#ifdef PYTHON_PROPERTIES
  void _get_uvw(double numpy_double_out[3])
#else
  void get_uvw(double numpy_double_out[3])
#endif
  {
    for (int i = 0; i < 3; i++) {
      numpy_double_out[i] = (*self)[i];
    }
  }
  //! Allow addition and subtraction with any Python iterable of 3 numbers
  %pythoncode %{
    @staticmethod
    @safesplat_float
    def _new_coord_frac(u, v, w):
        return Coord_frac(u, v, w)
        
    def __add__(self, other):
      try:
        return self.__add_base__(other)
      except:
        return self + self._new_coord_frac(other)

    def __sub__(self, other):
      try:
        return self.__sub_base__(other)
      except:
        return self - self._new_coord_frac(other)

  %}


#ifdef PYTHON_PROPERTIES
  void _get_uvw(double numpy_double_out[3])
  {
    for (int i = 0; i < 3; i++) {
      numpy_double_out[i] = (*self)[i];
    }
  }

  %pythoncode %{
    uvw = property(_get_uvw)
  %}
#endif

} // extend Coord_frac

%extend Coord_map
{
  Coord_map __add_base__(const Coord_map &h2) { return *self + h2; }
  Coord_map __sub_base__(const Coord_map &h2) { return *self - h2; }
  Coord_map __neg__() { return -(*self); }
  Coord_map __mul__ ( const double &f ) { return f * (*self); }
  Coord_map __rmul__ ( const double &f ) { return f * (*self); }
  // No transforms for this coordinate type
  std::string __str__() {return self->format();}
  std::string __repr__() {return "Coord_map " + self->format();}
  //! Allow addition and subtraction with any Python iterable of 3 numbers
  %pythoncode %{
    @staticmethod
    @safesplat_float
    def _new_coord_map(u, v, w):
        return Coord_orth(u, v, w)
        
    def __add__(self, other):
      try:
        return self.__add_base__(other)
      except:
        return self + self._new_coord_map(other)

    def __sub__(self, other):
      try:
        return self.__sub_base__(other)
      except:
        return self - self._new_coord_map(other)

  %}


#ifdef PYTHON_PROPERTIES
  void _get_uvw(double numpy_double_out[3])
#else
  void get_uvw(double numpy_double_out[3])
#endif
  {
    for (int i = 0; i < 3; i++) {
      numpy_double_out[i] = (*self)[i];
    }
  }
#ifdef PYTHON_PROPERTIES
  %pythoncode %{
    uvw = property(_get_uvw)
    ceil = property(ceil)
    coord_grid = property(coord_grid)
    floor = property(floor)
  %}
#endif

} // extend Coord_map

%extend HKL_sampling
{
  bool __eq__(const HKL_sampling& hkl2) { return *self == hkl2; }
  bool __ne__(const HKL_sampling& hkl2) { return !(*self == hkl2); }
} // extend HKL_sampling

%extend Grid
{
#ifdef PYTHON_PROPERTIES
  // TODO: Make dim a native numpy out function
  %pythoncode %{
    @property
    def dim(self):
      import numpy
      return numpy.array([self.nu(), self.nv(), self.nw()])
    
    size = property(size)
  %}
#endif
} // extend Grid

%extend Grid_sampling {
  std::string __str__() const { return self->format(); }
  void dim(int numpy_int_out[3])
  {
    numpy_int_out[0] = self->nu();
    numpy_int_out[1] = self->nv();
    numpy_int_out[2] = self->nw();
  }

#ifdef PYTHON_PROPERTIES
  %pythoncode %{
    dim = property(dim)
    format = property(format)
  %}
#endif
} // extend Grid_sampling
} // namespace clipper
