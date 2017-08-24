%inline %{
  namespace clipper {
  template <class T> struct matrixRowClipper {
    Mat33<T> *mat;
    int row; // Row number
    T __getitem__(int i)
    {
      return (*mat)(row,i);
    };
    void __setitem__(int i, double val)
    {
      (*mat)(row,i) = val;
    };
  };
  }
%}


namespace clipper
{
%extend Vec3 {

  Vec3 ( T v[3] )
  {
    Vec3<T> *new_v = new Vec3<T>( v[0], v[1], v[2] );
    return new_v;
  }

  T __getitem__ (int i)
  {
    i = (i < 0) ? 3 + i : i;
    if (i > 2 || i < 0) {
      throw std::out_of_range("");
    }
    return (*self)[i];
  }

  void __setitem__ (int i, T value)
  {
    i = (i < 0) ? 3 + i : i;
    if (i > 2 || i < 0) {
      throw std::out_of_range("");
    }
    (*self)[i] = value;
  }

  std::string __str__() const { return self-> format(); }

#ifdef PYTHON_PROPERTIES
  %pythoncode %{
    format = property(format)
  %}
#endif
} // extend Vec3




%extend Mat33 {
  matrixRowClipper<T> __getitem__(int i)
  {
    i = (i < 0) ? 3 + i : i;
    matrixRowClipper<T> r;
    if (i > 2 || i < 0) {
      throw std::out_of_range("");
    }
    r.mat = self;
    r.row = i;
    return r;
  };
  void as_numpy(double numpy_double_out[3][3])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = (*self)(i,j);
      }
    }
  }

  Mat33<T> __neg__ ()
  {
    return -(*self);
  };
  Vec3<T> __mul__ (const Vec3<T> &v)
  {
    return ((*self)*v);
  };
  Vec3<T> __rmul__ (const Vec3<T> &v)
  {
    return (v*(*self));
  };
  Mat33<T> __mul__ (const Mat33<T> &m)
  {
    return (*self*m);
  };
  Mat33<T> __rmul__ (const Mat33<T> &m)
  {
    return (m*(*self));
  };
  Mat33<T> __sub__ (const Mat33<T> &m)
  {
    return ((*self)+-m);
  };
  Mat33<T> __add__ (const Mat33<T> &m)
  {
    return ((*self)+m);
  };

  clipper::Coord_orth __mul__ (const clipper::Coord_orth &v_)
  {
    Vec3<T> v;
    v[0] = (T)v_[0];
    v[1] = (T)v_[1];
    v[2] = (T)v_[2];
    Vec3<T> v2 = ((*self)*v);
    Coord_orth c(v2[0],v2[1],v2[2]);
    return c;
  };
  clipper::Coord_frac __mul__ (const clipper::Coord_frac &v_)
  {
    Vec3<T> v;
    v[0] = (T)v_[0];
    v[1] = (T)v_[1];
    v[2] = (T)v_[2];
    Vec3<T> v2 = ((*self)*v);
    Coord_frac c(v2[0],v2[1],v2[2]);
    return c;
  };

#ifdef PYTHON_PROPERTIES
  %pythoncode %{
    as_numpy = property(as_numpy)
  %}
#endif

}; // extend Mat33
} // namespace clipper

%template (mat33_float)  clipper::Mat33<float>;
//%template (mat33_ftype)  clipper::Mat33<ftype>;
%template (mat33_double) clipper::Mat33<double>;
%template (vec3_int)    clipper::Vec3<int>;
%template (vec3_float)  clipper::Vec3<float>;
%template (vec3_double) clipper::Vec3<double>;

//namespace clipper
//{
//%extend mat33_float
//{
  //clipper::Coord_orth __mul__ (const clipper::Coord_orth &v_)
  //{
    //Vec3<float> v;
    //v[0] = v_[0];
    //v[1] = v_[1];
    //v[2] = v_[2];
    //Vec3<float> v2 = ((*self)*v);
    //Coord_orth c(v2[0],v2[1],v2[2]);
    //return c;
  //};
  //clipper::Coord_frac __mul__ (const clipper::Coord_frac &v_)
  //{
    //Vec3<float> v;
    //v[0] = v_[0];
    //v[1] = v_[1];
    //v[2] = v_[2];
    //Vec3<float> v2 = ((*self)*v);
    //Coord_frac c(v2[0],v2[1],v2[2]);
    //return c;
  //};
//}

//}
