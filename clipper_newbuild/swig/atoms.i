namespace clipper
{
%extend Atom
{
#ifndef PYTHON_PROPERTIES
  std::string element() { return self->element(); }
#endif

#ifdef PYTHON_PROPERTIES
  %pythoncode %{
    _allow_unknown = True
    
    @property
    def allow_unknown(self):
      return self._allow_unknown
    
    @allow_unknown.setter
    def allow_unknown(self, allow):
      self._allow_unknown = allow
    
    def _set_element_with_check(self, element_name):
      if not self.allow_unknown:
        if element_name not in ('H', 'C', 'N', 'O', 'S'):
          if element_name not in ATOM_NAMES:
            raise TypeError('Unrecognised element!')
      self._set_element(element_name)
        
    element = property(_get_element, _set_element_with_check)
    occupancy = property(_get_occupancy, _set_occupancy)
    u_iso = property(_get_u_iso, _set_u_iso)
    _u_aniso = property(_get_u_aniso, _set_u_aniso)
    coord_orth = property(_get_coord_orth, _set_coord_orth)
    is_null = property(is_null)

    def _get_coord_xyz(self):
      ''' Get the orthographic (x,y,z) coordinates as a Numpy array. '''
      return self.coord_orth.xyz
    def _set_coord_xyz(self, coord):
      self.coord_orth = Coord_orth(*coord)
    coord = property(_get_coord_xyz, _set_coord_xyz)

    def _get_u_aniso_vals(self):
      '''
      Anisotropic B-factor matrix as a 6-member array:
      [u00, u11, u22, u01, u02, u12].
      For purely isotropic B-factors, set this to None
      '''
      return self._u_aniso._get_vals()

    def _set_u_aniso_vals(self, u_aniso):
      import numpy
      if u_aniso is None:
        from math import nan
        self._set_u_aniso(U_aniso_orth(*([nan]*6)))
      else:
        try:
          self._set_u_aniso(U_aniso_orth(*u_aniso))
        except:
          if type(u_aniso) == numpy.ndarray:
            self._set_u_aniso(U_aniso_orth(*(u_aniso.astype(float))))
          elif type(u_aniso) == U_aniso_orth:
            self._set_u_aniso(u_aniso)
          else:
            raise TypeError('''
              u_aniso must be None, a list of 6 numbers [u00, u11, u22, u01, u02, u12]\n
              or a clipper U_aniso_orth object.''')

    u_aniso = property(_get_u_aniso_vals, _set_u_aniso_vals)
  %}

#endif
} // extend Atom

%extend Atom_list
/* There is probably a bit more work to do here to make this object really
 * useful for Python. The fundamental problem is that the clipper Atom_list
 * is just a std::vector of Atom objects, meaning that a naive slicing 
 * implementation would return a new Atom_list object filled with copies
 * of the atoms rather than pointers. So Atom_list[0:5].coords = array, while 
 * appearing perfectly legal in Python, would not actually do anything to 
 * the stored atoms.
 * At present the only ways to modify the actual Atoms stored in the array
 * are either one at a time, or all at once. It would be much better to 
 * create a C++ wrapper class that safely stores the actual Atom_list 
 * object, but provides individual/arrays of Atom pointers on indexing 
 * or slicing. We could also then easily create a constructor to initialise 
 * directly from arrays of atom properties.
 *
 */
{
  Atom& __getitem__(int i)
  {
    int array_len = self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    return (*self)[i];
  }

  void __setitem__(int i, Atom atom)
  {
    int array_len = self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    (*self)[i] = atom;
    return;
  }
  size_t __len__()
  {
    return self->size();
  }

  Atom pop(int i)
  {
    size_t array_len = self->size();
    if (i >= array_len || -i > array_len) {
      throw std::out_of_range("");
    }
    Atom ret;
    if (i>=0) ret = (*self)[i];
    else ret = (*self)[array_len - i];
    self -> erase(self -> begin() + i);
    return ret;
  }

  size_t append (Atom& a)
  {
    self->push_back(a);
    return self->size();
  }

  void extend_by(size_t n)
  {
    for (size_t i = 0; i < n; i++ ) {
      self->push_back( Atom() );
    }
  }

#ifdef PYTHON_PROPERTIES
  void _set_elements_base(std::vector<std::string> elements)
#else
  void set_elements(std::vector<std::string> elements)
#endif
  {
    size_t in_len = elements.size();
    size_t my_len = self -> size();
    if (in_len != my_len) {
      std::string errString;
      errString = "Input array length of " + std::to_string(in_len)
                  + " does not match target array length of " + std::to_string(my_len);
      throw std::length_error(errString);
    }
    for (size_t i = 0; i < my_len; i++) {
      (*self)[i].set_element( elements[i] );
    }
  }

#ifdef PYTHON_PROPERTIES
  std::vector< std::string > _get_elements_base()
#else
  std::vector< std::string > get_elements()
#endif
  {
    std::vector < std::string > ret;
    for (size_t i = 0; i < self -> size(); i++) {
      ret.push_back( (*self)[i].element() );
    }
    return ret;
  }

  //! Quickly set all atomic xyz coordinates from an nx3 numpy array
#ifdef PYTHON_PROPERTIES
  void _set_coord_orth_base(double *numpy_2d_in, int n1, int n2)
#else
  void set_coord_orth(double *numpy_2d_in, int n1, int n2)
#endif
  {
    size_t my_len = self -> size();
    if (n1 != my_len) {
      std::string errString = "Input array length of " + std::to_string(n1) +
                              " does not match target array length of " + std::to_string(my_len);
      throw std::length_error(errString);
    }
    if (n2 != 3) {
      throw std::length_error("Coordinates should be in the form of an N x 3 array");
    }
    for (size_t i = 0; i < n1; i++) {
      size_t first = i * n2;
      (*self)[i].set_coord_orth( Coord_orth( numpy_2d_in[first], numpy_2d_in[first+1], numpy_2d_in[first+2] ) );
    }
    return;
  }

  //! Quickly fill an nx3 numpy array with all atomic xyz coordinates
#ifdef PYTHON_PROPERTIES
  void _get_coord_orth_base(double* numpy_array, int n1, int n2)
#else
  void get_coord_orth(double* numpy_array, int n1, int n2)
#endif
  {
    size_t my_len = self -> size();
    if (n1 != my_len) {
      std::string errString = "Target array length of " + std::to_string(n1) +
                              " does not match Atom_list length of " + std::to_string(my_len);
      throw std::length_error(errString);
    }
    if (n2 != 3) {
      throw std::length_error("Coordinates should be in the form of an N x 3 array");
    }
    for (size_t i = 0; i < n1; i++) {
      size_t first = i * n2;
      Coord_orth thiscoord = (*self)[i].coord_orth();
      numpy_array[first] = thiscoord[0];
      numpy_array[first+1] = thiscoord[1];
      numpy_array[first+2] = thiscoord[2];
    }
  }

#ifdef PYTHON_PROPERTIES
  void _set_occupancies_base(double *numpy_1d_in, int n)
#else
  void set_occupancies(double *numpy_1d_in, int n)
#endif
  {
    size_t my_len = self -> size();
    if (n != my_len) {
      std::string errString = "Input array length of " + std::to_string(n) +
                              " does not match target array length of " + std::to_string(my_len);
      throw std::length_error(errString);
    }
    for (size_t i = 0; i < n; i++) {
      (*self)[i].set_occupancy( numpy_1d_in[i] );
    }
    return;
  }

#ifdef PYTHON_PROPERTIES
  void _get_occupancies_base(double *numpy_array, int n)
#else
  void get_occupancies(double *numpy_array, int n)
#endif
  {
    size_t my_len = self -> size();
    if (n != my_len) {
      std::string errString = "Target array length of " + std::to_string(n) +
                              " does not match Atom_list length of " + std::to_string(my_len);
      throw std::length_error(errString);
    }
    for (size_t i = 0; i < n; i++) {
      numpy_array[i] = (*self)[i].occupancy();
    }
  }

#ifdef PYTHON_PROPERTIES
  void _set_u_isos_base(double *numpy_1d_in, int n)
#else
  void set_u_isos(double *numpy_1d_in, int n)
#endif
  {
    size_t my_len = self -> size();
    if (n != my_len) {
      std::string errString = "Input array length of " + std::to_string(n) +
                              " does not match target array length of " + std::to_string(my_len);
      throw std::length_error(errString);
    }
    for (int i = 0; i < n; i++) {
      (*self)[i].set_u_iso( numpy_1d_in[i] );
    }
    return;
  }

#ifdef PYTHON_PROPERTIES
  void _get_u_isos_base(double *numpy_array, int n)
#else
  void get_u_isos(double *numpy_array, int n)
#endif
  {
    size_t my_len = self -> size();
    if (n != my_len) {
      std::string errString = "Target array length of " + std::to_string(n) +
                              " does not match Atom_list length of " + std::to_string(my_len);
      throw std::length_error(errString);
    }
    for (size_t i = 0; i < n; i++) {
      numpy_array[i] = (*self)[i].u_iso();
    }
  }

#ifdef PYTHON_PROPERTIES
  void _set_u_anisos_base(double *numpy_2d_in, int n1, int n2)
#else
  void set_u_anisos(double *numpy_2d_in, int n1, int n2)
#endif
  {
    size_t my_len = self -> size();
    if (n1 != my_len) {
      std::string errString = "Input array length of " + std::to_string(n1) +
                              " does not match target array length of " + std::to_string(my_len);
      throw std::length_error(errString);
    }
    if (n2 != 6) {
      std::string errString = "Input should be in the form of a N x 6 matrix, where"
                              "each row is of the form [u11, u22, u33, u12, u13, u23]";
      throw std::length_error(errString);
    }
    double *ai = numpy_2d_in;
    for (size_t i = 0; i < n1; i++) {
      size_t first = i*n2;
      (*self)[i].set_u_aniso_orth(
        U_aniso_orth(ai[first], ai[first+1], ai[first+2],
                     ai[first+3], ai[first+4], ai[first+5]) );
    }
  }

#ifdef PYTHON_PROPERTIES
  void _get_u_anisos_base(double* numpy_array, int n1, int n2)
#else
  void get_u_anisos(double* numpy_array, int n1, int n2)
#endif
  {
    size_t my_len = self -> size();
    if (n1 != my_len) {
      std::string errString = "Target array length of " + std::to_string(n1) +
                              " does not match Atom_list length of " + std::to_string(my_len);
      throw std::length_error(errString);
    }
    if (n2 != 6) {
      throw std::length_error("Target should be in the form of an N x 6 array");
    }
    for (size_t i = 0; i < n1; i++) {
      size_t first = i * n2;
      U_aniso_orth thisaniso = (*self)[i].u_aniso_orth();
      numpy_array[first] = thisaniso.mat00();
      numpy_array[first+1] = thisaniso.mat11();
      numpy_array[first+2] = thisaniso.mat22();
      numpy_array[first+3] = thisaniso.mat01();
      numpy_array[first+4] = thisaniso.mat02();
      numpy_array[first+5] = thisaniso.mat12();
    }
  }

  void get_minmax_grid(int numpy_int_out[2][3], const Cell& cell, const Grid_sampling& grid)
  {
    /* Find the minimum and maximum grid coordinates of a box encompassing the atoms,
     * and return a numpy array [min, max].
     */
    Coord_grid ref_min = (*self)[0].coord_orth().coord_frac(cell).coord_grid(grid);
    Coord_grid ref_max = ref_min;
    for (Atom_list::const_iterator it = self->begin(); it != self->end(); ++it) {
      Coord_grid thiscoord = it->coord_orth().coord_frac(cell).coord_grid(grid);
      for (size_t i = 0; i < 3; i++) {
        if (thiscoord[i] < ref_min[i]) ref_min[i] = thiscoord[i];
        else if (thiscoord[i] > ref_max[i]) ref_max[i] = thiscoord[i];
      }
    }
    for (size_t i = 0; i < 3; i++) {
      numpy_int_out[0][i] = ref_min[i];
      numpy_int_out[1][i] = ref_max[i];
    }
  }


  %pythoncode %{
    _allow_unknown = True
    
    @property
    def allow_unknown(self):
      return self._allow_unknown
    
    @allow_unknown.setter
    def allow_unknown(self, allow):
      self._allow_unknown = allow
    
#ifdef PYTHON_PROPERTIES

    def _get_elements(self):
      '''Ordered list of all element names'''
      return self._get_elements_base()
    
    def _set_elements(self, elements):
      if not self.allow_unknown:
        if not set(elements).issubset(ATOM_NAMES):
          bad_atoms = []
          for el in set(elements):
            if el not in ATOM_NAMES:
              bad_atoms.append(el)
          bad_atoms = set(bad_atoms)
          errstring = '''
            The following atom names are not recognised by Clipper:
            {}
            '''.format(bad_atoms)
          raise TypeError(errstring)
      self._set_elements_base(elements)
    
    elements = property(_get_elements, _set_elements)
    
          
      
    def _get_coord_orth(self):
      '''Orthographic (x,y,z) coordinates of all atoms'''
      import numpy
      n = len(self)
      coords = numpy.empty((n,3), numpy.double)
      self._get_coord_orth_base(coords)
      return coords
    
    def _set_coord_orth(self, coords):
      import numpy
      n = len(self)
      array_in = numpy.empty((n, 3), numpy.double)
      array_in[:] = coords
      self._set_coord_orth_base(array_in)
    
    coord_orth = property(_get_coord_orth, _set_coord_orth)
    
    def _get_occupancies(self):
      '''Occupancies of all atoms.'''
      import numpy
      n = len(self)
      occ = numpy.empty(n, numpy.double)
      self._get_occupancies_base(occ)
      return occ
    
    def _set_occupancies(self, occ):
      import numpy
      n = len(self)
      array_in = numpy.empty(n, numpy.double)
      array_in[:] = occ
      self._set_occupancies_base(array_in)
    
    occupancies = property(_get_occupancies, _set_occupancies)
    
    def _get_u_isos(self):
      '''Isotropic B-factors of all atoms.'''
      import numpy
      n = len(self)
      uisos = numpy.empty(n, numpy.double)
      self._get_u_isos_base(uisos)
      return uisos
    
    def _set_u_isos(self, uisos):
      import numpy
      n = len(self)
      array_in = numpy.empty(n, numpy.double)
      array_in[:] = uisos
      self._set_u_isos_base(array_in)
    
    u_isos = property(_get_u_isos, _set_u_isos)
    
    def _get_u_anisos(self):
      '''
      Anisotropic B-factor matrices for all atoms as an nx6 array, in the
      format: n*[u00, u11, u22, u01, u02, u12]. For purely isotropic
      atoms, set all elements in their row to math.nan or numpy.nan.
      '''
      import numpy
      n = len(self)
      uaniso = numpy.empty((n,6), numpy.double)
      self._get_u_anisos_base(uaniso)
      return uaniso
    
    def _set_u_anisos(self, u_anisos):
      import numpy
      n = len(self)
      array_in = numpy.empty((n,6), numpy.double)
      array_in[:] = u_anisos
      self._set_u_anisos_base(array_in)
    
    u_anisos = property(_get_u_anisos, _set_u_anisos)

#endif
    
  %}

} // extend Atom_list

%extend U_aniso_orth {
  U_aniso_orth __add__(const U_aniso_orth& u2) { return *self + u2; }
  U_aniso_orth __sub__(const U_aniso_orth& u2) { return *self + (-u2); }
  U_aniso_orth __neg__() { return -(*self); }
  U_aniso_orth __mul__(const double& f) { return f * (*self); }
  U_aniso_orth __rmul__(const double& f) { return f * (*self); }
  void _get_vals(double numpy_double_out[6])
  {
    numpy_double_out[0] = self->mat00();
    numpy_double_out[1] = self->mat11();
    numpy_double_out[2] = self->mat22();
    numpy_double_out[3] = self->mat01();
    numpy_double_out[4] = self->mat02();
    numpy_double_out[5] = self->mat12();
  }
} // extend U_aniso_orth

%extend U_aniso_frac {
  U_aniso_frac __add__(const U_aniso_frac& u2) { return *self + u2; }
  U_aniso_frac __sub__(const U_aniso_frac& u2) { return *self + (-u2); }
  U_aniso_frac __neg__() { return -(*self); }
  U_aniso_frac __mul__(const double& f) { return f * (*self); }
  U_aniso_frac __rmul__(const double& f) { return f * (*self); }
  void _get_vals(double numpy_double_out[6])
  {
    numpy_double_out[0] = self->mat00();
    numpy_double_out[1] = self->mat11();
    numpy_double_out[2] = self->mat22();
    numpy_double_out[3] = self->mat01();
    numpy_double_out[4] = self->mat02();
    numpy_double_out[5] = self->mat12();
  }

} // extend U_aniso_frac


} // namespace clipper
