%module(directors="1") clipper

%include "std_vector.i"
%include "std_string.i"
//%include "std_complex.i"
%include "exception.i"
%include "std_except.i"
//%include "typemaps.i"

#pragma SWIG nowarn=312,325,361,362,363,389,401,501,505

%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
    import_array();
%}

%apply int { size_t };

//%rename(clipper_Range) clipper::Range;
//%rename(clipper_Batch) clipper::datatypes::Batch;

%rename(SFscale_aniso_TYPE_float) clipper::SFscale_aniso<float>::TYPE;
%rename(SFscale_aniso_MODE_float) clipper::SFscale_aniso<float>::MODE;
%rename(SFscale_aniso_F_float) clipper::SFscale_aniso<float>::F;
%rename(SFscale_aniso_I_float) clipper::SFscale_aniso<float>::I;
%rename(SFscale_aniso_NORMAL_float) clipper::SFscale_aniso<float>::NORMAL;
%rename(SFscale_aniso_SHARPEN_float) clipper::SFscale_aniso<float>::SHARPEN;
%rename(SFscale_aniso_UNSHARPEN_float) clipper::SFscale_aniso<float>::UNSHARPEN;
%rename(MMonomer_TYPE) clipper::MMonomer::TYPE;
%rename(MMDBManager_TYPE) clipper::MMDBManager::TYPE;

/* Definitions for Numpy in/out/in-place arrays */
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *test_numpy_a, int test_numpy_n)};
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(double *numpy_array, int nu, int nv, int nw )};

%apply (double* IN_ARRAY1, int DIM1) {(double *numpy_1d_in, int n)};
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double *numpy_2d_in, int n1, int n2)};
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(double *numpy_3d_in, int n1, int n2, int n3)};
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(double *numpy_3d_in, int nu, int nv, int nw)};
%apply (int ARGOUT_ARRAY1[ANY]) {(int numpy_int_out[3])};
%apply (double ARGOUT_ARRAY1[ANY]) {(double numpy_double_out[3])};
%apply (double ARGOUT_ARRAY1[ANY]) {(double numpy_double_out[6])};
/*

/


%inline 
%{ 
#include <stdio.h>
void ClipperTestPassNumpyArray(double *test_numpy_a, int test_numpy_n){
  printf("Hello world %d\n",test_numpy_n);
}
%}
*/
/* // Test with
>>> import clipper
>>> import numpy
>>> a = numpy.array([1,2,3,4,5,6],numpy.double)
>>> clipper.ClipperTestPassNumpyArray(a)
*/
/* END NUMPY EXAMPLE */


/*             --  Rename overloaded friend operators --              */

%rename(neg_HKL)              operator-  (const HKL&);
%rename(add_HKL)              operator+  (const HKL&, const HKL&);
%rename(subs_HKL)             operator-  (const HKL&, const HKL&);
%rename(product_HKL)          operator*  (const int&, const HKL&);
%rename(transf_HKL)           operator*  (const Isymop&, const HKL&);

%rename(equals_Coord_grid)    operator== (const Coord_grid&, const Coord_grid&);
%rename(notequals_Coord_grid) operator!= (const Coord_grid&, const Coord_grid&); 
%rename(neg_Coord_grid)       operator-  (const Coord_grid&);
%rename(add_Coord_grid)       operator+  (const Coord_grid&, const Coord_grid&);
%rename(subs_Coord_grid)      operator-  (const Coord_grid&, const Coord_grid&);
%rename(product_Coord_grid)   operator*  (const int&, const Coord_grid&);
%rename(transf_Coord_grid)    operator*  (const Isymop&, const Coord_grid&);

%rename(neg_Coord_orth)       operator-  (const Coord_orth&);
%rename(add_Coord_orth)       operator+  (const Coord_orth&, const Coord_orth&);
%rename(subs_Coord_orth)      operator-  (const Coord_orth&, const Coord_orth&);
%rename(product_Coord_orth)   operator*  (const ftype&, const Coord_orth&);
%rename(transf_Coord_orth)    operator*  (const RTop_orth&, const Coord_orth&);

%rename(neg_Coord_frac)       operator-  (const Coord_frac&);
%rename(add_Coord_frac)       operator+  (const Coord_frac&, const Coord_frac&);
%rename(subs_Coord_frac)      operator-  (const Coord_frac&, const Coord_frac&);
%rename(product_Coord_frac)   operator*  (const ftype&, const Coord_frac&);
%rename(transf_Coord_frac)    operator*  (const RTop_frac&, const Coord_frac&);

%rename(neg_Coord_map)        operator-  (const Coord_map&);
%rename(add_Coord_map)        operator+  (const Coord_map&, const Coord_map&);
%rename(subs_Coord_map)       operator-  (const Coord_map&, const Coord_map&);
%rename(product_Coord_map)    operator*  (const ftype&, const Coord_map&);

%rename(neg_U_aniso_orth)     operator-  (const U_aniso_orth&);
%rename(add_U_aniso_orth)     operator+  (const U_aniso_orth&, const U_aniso_orth&);
%rename(product_U_aniso_orth) operator*  (const ftype&, const U_aniso_orth&);

%rename(neg_U_aniso_frac)     operator-  (const U_aniso_frac&);
%rename(add_U_aniso_frac)     operator+  (const U_aniso_frac&, const U_aniso_frac&);
%rename(product_U_aniso_frac) operator*  (const ftype&, const U_aniso_frac&);

%rename(equals_HKL_samp)      operator== (const HKL_sampling&, const HKL_sampling& );

%rename(and_MMonomer)         operator&  (const MMonomer&, const MMonomer&);
%rename(and_MPolymer)         operator&  (const MPolymer&, const MPolymer&);
%rename(and_MModel)           operator&  (const MModel&, const MModel&);

%rename(or_MMonomer)          operator|  (const MMonomer&, const MMonomer&);
%rename(or_MPolymer)          operator|  (const MPolymer&, const MPolymer&);
%rename(or_MModel)            operator|  (const MModel&, const MModel&);

/*                -- End of friend funcion redefinitions --              */



%ignore clipper::CCP4MTZfile::assigned_paths;
%ignore clipper::MMDBManager::write_file(const String&);
%ignore clipper::TargetFn_base::debug;
%ignore clipper::ResolutionFn::debug;


namespace std {
  %template(FloatVector) vector<float>;
  %template(DoubleVector) vector<double>;
  %template(FloatFloatVector) vector<vector<float> >;
  %template(DoubleDoubleVector) vector<vector<double> >;
}

/*
 * Renaming of existing clipper functions must go here, BEFORE
 * the corresponding header is #included.  
*/
namespace clipper {
  %rename("_u_ptr") Coord_grid::u();
  %rename("_v_ptr") Coord_grid::v();
  %rename("_w_ptr") Coord_grid::w();
}



%{
    #include <string>
    #include <limits>
    #include "../clipper/clipper.h"
    #include "../clipper/core/clipper_types.h"
    #include "../clipper/core/hkl_lookup.h"
    #include "../clipper/core/hkl_info.h"
    #include "../clipper/core/xmap.h"
    #include "../clipper/core/nxmap.h"
    #include "../clipper/core/coords.h"
    #include "../clipper/ccp4/ccp4_map_io.h"
    #include "../clipper/ccp4/ccp4_mtz_io.h"
    #include "../clipper/ccp4/ccp4_mtz_types.h"
/*
    #include "../clipper/cns/cns_map_io.h"
    #include "../clipper/cns/cns_hkl_io.h"
*/
    #include "../clipper/core/atomsf.h"
    #include "../clipper/minimol/minimol_seq.h"
    #include "../clipper/minimol/minimol.h"
    #include "../clipper/minimol/minimol_io.h"
    #include "../clipper/minimol/minimol_utils.h"
    #include "../clipper/minimol/minimol_seq.h"
    #include "../clipper/ccp4/ccp4_utils.h"
    #include "../clipper/core/map_utils.h"
    #include "../clipper/core/clipper_stats.h"
    #include "../clipper/mmdb/clipper_mmdb.h"
    #include "../clipper/core/hkl_datatypes.h"
    #include "../clipper/core/hkl_data.h"
    #include "../clipper/core/hkl_compute.h"
    #include "../clipper/cif/cif_data_io.h"
    #include "../clipper/core/ramachandran.h"
    #include "../clipper/contrib/sfcalc_obs.h"
    #include "../clipper/contrib/sfweight.h"
    #include "../clipper/contrib/sfcalc.h"
    #include "../clipper/contrib/edcalc.h"
    #include "../clipper/core/nxmap_operator.h"
    #include "../clipper/core/resol_fn.h"
    #include "../clipper/core/resol_basisfn.h"
    #include "../clipper/core/resol_targetfn.h"
    #include "../clipper/contrib/convolution_search.h"
    #include "../clipper/contrib/mapfilter.h"
    #include "../clipper/contrib/originmatch.h"
    #include "../clipper/contrib/sfscale.h"
    #include "../clipper/core/rotation.h"

    namespace clipper 
    {
        enum errtype {
            NO_ERROR,
            INDEX_OUT_OF_BOUNDS,
            ARRAY_LENGTH_MISMATCH,
            CUSTOM_INDEX_ERROR,
            
            UNSPECIFIED_ERROR
            
        };
        static int myErr = 0; // flag to save error state
        static errtype myErrType = NO_ERROR; // flag to specify error type to raise
        std::string myErrString = ""; // Optional string for custom error messages
        
        void raise_swig_err (errtype err) {
            
            switch(err) {
                case INDEX_OUT_OF_BOUNDS: {
                    SWIG_exception(SWIG_IndexError, "Index out of bounds");
                    break;
                  }
                case ARRAY_LENGTH_MISMATCH: {
                    SWIG_exception(SWIG_IndexError, "Input array length does not match target array length");
                    break;
                  }
                case CUSTOM_INDEX_ERROR: {
                    std::string thisErrString = myErrString;
                    myErrString = "";
                    SWIG_exception(SWIG_IndexError, thisErrString.c_str());
                    break;
                  }
                case UNSPECIFIED_ERROR: {
                    SWIG_exception(SWIG_UnknownError, "Unspecified error");
                    break;
                  }
                
                
            }
            fail: return;
        }
    }

    using namespace clipper;
    #include <string.h>
%}

/* Macros for handling of common exceptions */

#define CATCH_SWIG_EXCEPTION(f) \
  %exception f {\
    assert(!myErr);\
    $action\
    if (myErr) {\
      myErr = 0;\
      errtype thisErrType = myErrType;\
      myErrType = NO_ERROR;\
      raise_swig_err(thisErrType);\
    }\
  }\
  // Blank line required after macro



namespace clipper 
{
    %ignore Vec3<int>;
}

%inline 
%{
    namespace clipper
    {
        typedef ftype64  ftype;
        typedef ftype64  xtype;
        typedef float    ftype32;
        typedef double   ftype64;
    }
    
    std::string ClipperStringAsString(const clipper::String &a) 
    {
        return (std::string)(a);
    }
%}

%ignore clipper::U32;
%ignore clipper::U64;

%director clipper::Container;
%director clipper::TargetFn_base;
%director clipper::HKL_data_base;
%director clipper::BasisFn_base;
%director clipper::OriginMatch_base;
%director clipper::SFscale_base;
%director clipper::SFweight_base;
%director clipper::EDcalc_base;
%director clipper::SFcalc_base;
%director clipper::Util::Vec3;

namespace std 
{
   %template(UnsignedIntVector) vector<unsigned int>;
   %template(IntVector) vector<int>;
   %template(IntIntVector) vector<vector<int> >;
//   %template(DoubleVector) vector<double>;
//   %template(DoubleDoubleVector) vector<vector<double> >;
   %template(ClipperStringVector) vector<clipper::String>;
   %template(HKLVector) vector<clipper::HKL>;
   // Do not know if these 3 will be required.
   //%template(Coord_orthVector) vector<clipper::Coord_orth>;
   //%template(MAtomIndexVector) vector<clipper::MAtomIndex>;                 // Dodgy ...?
   //%template(MAtomIndexSymmetryVector) vector<clipper::MAtomIndexSymmetry>; // Dodgy ...?
   %template(StringVector) vector<string>;
}

%apply std::string { clipper::String }
%apply std::string& { clipper::String& }

%apply std::string { String }
%apply std::string& { String& }
//%apply const std::string& { const String& } 

%typemap(in) (const clipper::String&)
{
%#if PY_MAJOR_VERSION >= 3
   char *my_result;

   if (PyUnicode_Check($input)) {
     PyObject * temp_bytes = PyUnicode_AsEncodedString($input, "ASCII", "strict"); // Owned reference
     if (temp_bytes != NULL) {
        my_result = PyBytes_AS_STRING(temp_bytes); // Borrowed pointer
        my_result = strdup(my_result);
        Py_DECREF(temp_bytes);
     } else {
         std::cout << "Decoding error" << std::endl;
     }
   }

   clipper::String *s = new clipper::String(my_result);
%#else
   std::string ss = PyString_AsString($input);
   clipper::String *s = new clipper::String(ss);
%#endif
   $1 = s;
}

namespace clipper {
  %extend String {
    std::string __str__() {
      return ($self)->c_str();
    }
  }
}


%include "../clipper/core/clipper_util.h"
%include "../clipper/core/clipper_types.h"

%template (RTop_float)  clipper::RTop<float>;
%template (RTop_double) clipper::RTop<double>;
%template (vec3_int)    clipper::Vec3<int>;
%template (vec3_float)  clipper::Vec3<float>;
%template (vec3_double) clipper::Vec3<double>;

%inline %{
namespace clipper {
  struct matrixRowClipper {
    Mat33<float> *mat;
    int row; // Row number
    float __getitem__(int i) {
      return (*mat)(row,i);
    };
    void __setitem__(int i, double val) {
      (*mat)(row,i) = val;
    };
  };
}
%}

%template (mat33_float)  clipper::Mat33<float>;
%template (mat33_ftype)  clipper::Mat33<ftype>;

namespace clipper
{
    %ignore Cell::matrix_orth() const;
    %ignore Cell::matrix_frac() const;
    %rename(is_nan_float) Util::is_nan(const ftype32);
    %rename(is_nan_double) Util::is_nan(const ftype64);
    %rename(is_nan_float_slow) Util::isnan(ftype32);
    %rename(is_nan_double_slow) Util::isnan(ftype64);
    %rename(set_null_float) Util::set_null(ftype32);
    %rename(set_null_double) Util::set_null(ftype64);
    %rename(is_null_float) Util::is_null(ftype32);
    %rename(is_null_double) Util::is_null(ftype64);
}



  /* We are getting a load of warnings whenever SWIG is trying 
     to wrap template base classes that will never be of any use
     in Python, so let's suppress the warnings for these 
  */

%include "../clipper/core/symop.h"


  /* Here we use immutable to deal with const char * vars */
 
%immutable hall;
%immutable hm;
%immutable lgname;

%include "../clipper/core/spacegroup_data.h"
%feature ("flatnested","1");
%include "../clipper/core/spacegroup.h"
%feature ("flatnested","0");

namespace clipper
{

  //! Metric tensor
  /*! The metric tensor is used to determine a distance in real or
    reciprocal space using fraction coordinates or Miller indices. It
    is symmetrical, so only the upper triangle is stored with the
    off-diagonal elements doubled.
  */
  class Metric_tensor
  {
  public:
    //! null constructor
    inline Metric_tensor() {} //!< Null constructor
    //! constructor: takes parameters of normal or inverse cell
    Metric_tensor( const ftype& a, const ftype& b, const ftype& c, const ftype& alph, const ftype& beta, const ftype& gamm );
    //! apply metric to vector
    inline ftype lengthsq( const Vec3<>& v ) const
      { return ( v[0]*(v[0]*m00 + v[1]*m01 + v[2]*m02) +
                 v[1]*(v[1]*m11 + v[2]*m12) + v[2]*(v[2]*m22) ); }
    //! apply metric to int vector
    inline ftype lengthsq( const Vec3<int>& v ) const
      { ftype h = ftype(v[0]); ftype k = ftype(v[1]); ftype l = ftype(v[2]);
        return h*(h*m00 + k*m01 + l*m02) + k*(k*m11 + l*m12) + l*(l*m22); }

    String format() const;  //!< return formatted String representation
  };


  //! cell description (automatically converts to radians)
  /*! The cell description is a compact description of a cell,
    containing just the cell parameters. It is usually used to
    construct a full Cell object, which provides the expected
    functionality.
  */
  class Cell_descr
  {
  public:
    inline Cell_descr() {}  //!< null constructor
    //! constructor: from cell parameters
    Cell_descr( const ftype& a, const ftype& b, const ftype& c,
    const ftype& alpha=90.0f, const ftype& beta=90.0f,
    const ftype& gamma=90.0f );
    inline const ftype& a() const { return a_; } //!< get a
    inline const ftype& b() const { return b_; } //!< get b
    inline const ftype& c() const { return c_; } //!< get c
    inline const ftype& alpha() const { return alpha_; } //!< get alpha
    inline const ftype& beta() const { return beta_; }   //!< get beta
    inline const ftype& gamma() const { return gamma_; } //!< get gamma
    ftype alpha_deg() const; //!< get alpha in degrees
    ftype beta_deg() const;  //!< get alpha in degrees
    ftype gamma_deg() const; //!< get gamma in degrees
    String format() const;   //!< return formatted String representation

  };


  //! Cell object
  /*! The Cell class is the fully functional description of the unit
    cell. In addition to the cell parameters, it stores derived
    information including the cell volume, orthogonalising and
    fractionalising matrices, and the metric tensors.
   */
  class Cell : public Cell_descr
  {
   public:
    //! null constructor: must initialise later
    inline Cell() { vol = 0.0; }
    //! constructor: takes a Cell descriptor
    explicit Cell( const Cell_descr& cell_ ) { init( cell_ ); }
    //! initialiser
    void init( const Cell_descr& cell_ );

    //! test if object has been initialised
    bool is_null() const;

    ftype a_star() const; //!< get a*
    ftype b_star() const; //!< get b*
    ftype c_star() const; //!< get c*
    ftype alpha_star() const; //!< get alpha*
    ftype beta_star()  const; //!< get beta*
    ftype gamma_star() const; //!< get gamma*
    // inherited functions listed for documentation purposes
    //-- const ftype& a() const;
    //-- const ftype& b() const;
    //-- const ftype& c() const;
    //-- const ftype& alpha() const;
    //-- const ftype& beta() const;
    //-- const ftype& gamma() const;
    //-- ftype alpha_deg() const;
    //-- ftype beta_deg() const;
    //-- ftype gamma_deg() const;
    //-- String format() const;

    //! return cell dimensions
    //inline const Cell_descr& descr() const { return (*this); }
    //! return cell volume
    inline const ftype& volume() const { return vol; }
    //! test equality with another cell
    bool equals( const Cell& other, const ftype tol=1.0 ) const;
    //! return orthogonalisation matrix
    inline const Mat33<>& matrix_orth() const { return orthmat; }
    //! return fractionalisation matrix
    inline const Mat33<>& matrix_frac() const { return fracmat; }
    //! return real space metric tensor
    inline const Metric_tensor& metric_real() const { return realmetric; }
    //! return reciprocal space metric tensor   
    inline const Metric_tensor& metric_reci() const { return recimetric; }

    void debug() const;

  };

} // namespace clipper

namespace clipper {
  %extend Cell {
      Mat33<float> matrix_orth() {
      Mat33<float> orth;
      orth(0,0) = (self->matrix_orth())(0,0);
      orth(0,1) = (self->matrix_orth())(0,1);
      orth(0,2) = (self->matrix_orth())(0,2);
      orth(1,0) = (self->matrix_orth())(1,0);
      orth(1,1) = (self->matrix_orth())(1,1);
      orth(1,2) = (self->matrix_orth())(1,2);
      orth(2,0) = (self->matrix_orth())(2,0);
      orth(2,1) = (self->matrix_orth())(2,1);
      orth(2,2) = (self->matrix_orth())(2,2);
      return orth;
    };
      Mat33<float> matrix_frac() {
      Mat33<float> frac;
      frac(0,0) = (self->matrix_frac())(0,0);
      frac(0,1) = (self->matrix_frac())(0,1);
      frac(0,2) = (self->matrix_frac())(0,2);
      frac(1,0) = (self->matrix_frac())(1,0);
      frac(1,1) = (self->matrix_frac())(1,1);
      frac(1,2) = (self->matrix_frac())(1,2);
      frac(2,0) = (self->matrix_frac())(2,0);
      frac(2,1) = (self->matrix_frac())(2,1);
      frac(2,2) = (self->matrix_frac())(2,2);
      return frac;
    };
  };
}


/*
%typemap(in) (const clipper::String&)
{
   std::string ss = PyString_AsString($input);
   clipper::String *s = new clipper::String(ss);
   $1 = s;
}
*/


%include "../clipper/core/coords.h"


namespace clipper {
  CATCH_SWIG_EXCEPTION(Vec3<float>::__getitem__)
  CATCH_SWIG_EXCEPTION(Vec3<float>::__setitem__)
  

  %extend RTop_orth
  {
    RTop_orth ( Mat33<ftype>& mat )
    {
      return new RTop_orth(mat);
    }
  }

  %extend Vec3<float> {
    float __getitem__ (int i) {
      i = (i < 0) ? 3 + i : i;
      if (i > 2 || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return (*self)[0];
      }
      return (*self)[i];
      }

    void __setitem__ (int i, float value) {
      i = (i < 0) ? 3 + i : i;
      if (i > 2 || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return;        
      }
      (*self)[i] = value;
      }
  }


  CATCH_SWIG_EXCEPTION(Mat33<float>::__getitem__)


  %extend Mat33<float> {
    matrixRowClipper __getitem__(int i) {
      i = (i < 0) ? 3 + i : i;
      matrixRowClipper r;
      if (i > 2 || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return r;
      }
      r.mat = self;
      r.row = i;
      return r;
    };
    Mat33<float> __neg__ () {
      return -(*self);
    };
    Vec3<float> __mul__ (const Vec3<float> &v) {
      return ((*self)*v);
    };
    Vec3<float> __rmul__ (const Vec3<float> &v) {
      return ((*self)*v);
    };
    clipper::Coord_orth __mul__ (const clipper::Coord_orth &v_) {
      Vec3<float> v;
      v[0] = v_[0];
      v[1] = v_[1];
      v[2] = v_[2];
      Vec3<float> v2 = ((*self)*v);
      Coord_orth c(v2[0],v2[1],v2[2]);
      return c;
    };
    clipper::Coord_orth __rmul__ (const clipper::Coord_orth &v_) {
      Vec3<float> v;
      v[0] = v_[0];
      v[1] = v_[1];
      v[2] = v_[2];
      Vec3<float> v2 = ((*self)*v);
      Coord_orth c(v2[0],v2[1],v2[2]);
      return c;
    };
    clipper::Coord_frac __mul__ (const clipper::Coord_frac &v_) {
      Vec3<float> v;
      v[0] = v_[0];
      v[1] = v_[1];
      v[2] = v_[2];
      Vec3<float> v2 = ((*self)*v);
      Coord_frac c(v2[0],v2[1],v2[2]);
      return c;
    };
    clipper::Coord_frac __rmul__ (const clipper::Coord_frac &v_) {
      Vec3<float> v;
      v[0] = v_[0];
      v[1] = v_[1];
      v[2] = v_[2];
      Vec3<float> v2 = ((*self)*v);
      Coord_frac c(v2[0],v2[1],v2[2]);
      return c;
    };
    Mat33<float> __mul__ (const Mat33<float> &m) {
      return (*self*m);
    };
    Mat33<float> __rmul__ (const Mat33<float> &m) {
      return (*self*m);
    };
    Mat33<float> __sub__ (const Mat33<float> &m) {
      return ((*self)+-m);
    };
    Mat33<float> __add__ (const Mat33<float> &m) {
      return ((*self)+m);
    };
  };
}

%include "../clipper/core/hkl_info.h"
%include "../clipper/core/hkl_data.h"

namespace clipper {

class HKL_reference_base {
    public:
      const HKL_info& base_hkl_info() const; 
      const int& index() const;
      ftype invresolsq( const HKL_data_base& hkldata ) const;
      ftype invresolsq() const;
      bool last() const;
};
class HKL_reference_index : public HKL_reference_base {
    public:
      HKL_reference_index();
      HKL_reference_index( const HKL_info& hklinfo_, const int& index );
      const HKL& hkl() const;
      const HKL_class& hkl_class() const;
      HKL_reference_index& next();
};
}

%{
  namespace clipper {
    typedef HKL_info::HKL_reference_base HKL_reference_base;
    typedef HKL_info::HKL_reference_index HKL_reference_index;
  }
%}


%include "../clipper/core/xmap.h"
%include "../clipper/core/nxmap.h"


// numpy support for Xmap and NXmap
// jon is currently editing this code piece

namespace clipper
{
    %extend RTop<double>
    {
        std::string __str__( )
        {
            return (*($self)).format();
            fail: return "";
        }
    }

    %extend Mat33<ftype>
    {
        std::string __str__( )
        {
            return (*($self)).format();
            fail: return "";
        }
    }
    
    %extend Vec3<ftype>
    {
        std::string __str__( )
        {
            return (*($self)).format();
            fail: return "";
        }
    }

    %extend Mat33<double>
    {
        std::string __str__( )
        {
            return (*($self)).format();
            fail: return "";
        }
    }
    
    %extend Vec3<double>
    {
        std::string __str__( )
        {
            return (*($self)).format();
            fail: return "";
        }
    }

    %extend NXmap<float>
    {
        NXmap ( clipper::Grid const & grid, clipper::RTop<ftype> const & rtop)
        {
            return new NXmap<float>( grid, rtop );
        }
    
        int export_numpy ( double *numpy_array, int nu, int nv, int nw )
        {
            int i = 0;
            int top_u, top_v, top_w;
        
            clipper::Coord_grid c;
            clipper::Grid map_grid = (*($self)).grid();
        
            nu > map_grid.nu() ? top_u = map_grid.nu() : top_u = nu;
            nv > map_grid.nv() ? top_v = map_grid.nv() : top_v = nv;
            nw > map_grid.nw() ? top_w = map_grid.nw() : top_w = nw;
        
            for ( c.w() = 0; c.w() < top_w; c.w()++ )
                for ( c.v() = 0; c.v() < top_v; c.v()++ )
                    for (  c.u() = 0; c.u() < top_u; c.u()++, i++ )
                        numpy_array[i] = (*($self)).get_data(c);
            return i;
        }
    
        int import_numpy ( double *numpy_3d_in, int nu, int nv, int nw )
        {
            int i = 0;
            int top_u, top_v, top_w;
        
            clipper::Coord_grid c;
            clipper::Grid map_grid = (*($self)).grid();
        
            nu > map_grid.nu() ? top_u = map_grid.nu() : top_u = nu;
            nv > map_grid.nv() ? top_v = map_grid.nv() : top_v = nv;
            nw > map_grid.nw() ? top_w = map_grid.nw() : top_w = nw;
        
            for ( c.w() = 0; c.w() < top_w; c.w()++ )
                for ( c.v() = 0; c.v() < top_v; c.v()++ )
                    for (  c.u() = 0; c.u() < top_u; c.u()++, i++ )
                        (*($self)).set_data(c, numpy_3d_in[i]);
            return i;
        }
    
        RTop<ftype> operator_orth_grid () { return (*($self)).operator_orth_grid(); }
        RTop<ftype> operator_grid_orth () { return (*($self)).operator_grid_orth(); }
    
    }
    
    %extend Xmap<float>
    {
        int export_numpy ( double *numpy_array, int nu, int nv, int nw )
        {
            int i = 0;
            int top_u, top_v, top_w;
        
            clipper::Coord_grid c;
            clipper::Grid map_grid = (*($self)).grid_asu();
        
            nu > map_grid.nu() ? top_u = map_grid.nu() : top_u = nu;
            nv > map_grid.nv() ? top_v = map_grid.nv() : top_v = nv;
            nw > map_grid.nw() ? top_w = map_grid.nw() : top_w = nw;
        
            for ( c.w() = 0; c.w() < top_w; c.w()++ )
                for ( c.v() = 0; c.v() < top_v; c.v()++ )
                    for ( c.u() = 0; c.u() < nu; c.u()++, i++ )
                    {
                        if ( c.u() < map_grid.nu() && c.v() < map_grid.nv() && c.w() < map_grid.nw() )
                            numpy_array[i] = (*($self)).get_data(c);
                        else
                            numpy_array[i] = 0.0;
                    }
            return i;
        }
    
        int import_numpy ( double *numpy_3d_in, int nu, int nv, int nw )
        {
            int i = 0;
            int top_u, top_v, top_w;
        
            clipper::Coord_grid c;
            clipper::Grid map_grid = (*($self)).grid_asu();
        
            nu > map_grid.nu() ? top_u = map_grid.nu() : top_u = nu;
            nv > map_grid.nv() ? top_v = map_grid.nv() : top_v = nv;
            nw > map_grid.nw() ? top_w = map_grid.nw() : top_w = nw;
        
            for ( c.w() = 0; c.w() < top_w; c.w()++ )
                for ( c.v() = 0; c.v() < top_v; c.v()++ )
                    for (  c.u() = 0; c.u() < top_u; c.u()++, i++ )
                        (*($self)).set_data(c, numpy_3d_in[i]);
            return i;
        }
    
        int export_section_numpy ( double *numpy_array, int nu, int nv, int nw, clipper::Coord_grid& start, clipper::Coord_grid& end )
        {
            int i = 0;
            int w, v;
            clipper::Xmap_base::Map_reference_coord ix( (*($self)) );
        
            for ( w = start.w(); w <= end.w(); w++ )
                for ( v = start.v(); v <= end.v(); v++ )
                    for ( ix.set_coord(Coord_grid(start.u(),v,w)); ix.coord().u() <= end.u(); ix.next_u(), i++ )
                    {
                        numpy_array[i] = (*($self)).get_data(ix.coord());
                    }
            return i;
        }
    
        int import_section_numpy ( double *numpy_3d_in, int nu, int nv, int nw, clipper::Coord_grid& start, clipper::Coord_grid& end )
        {
            int i = 0;
            int w, v;
            clipper::Xmap_base::Map_reference_coord ix( (*($self)) );
        
            for ( w = start.w(); w <= end.w(); w++ )
                for ( v = start.v(); v <= end.v(); v++ )
                    for ( ix.set_coord(Coord_grid(start.u(),v,w)); ix.coord().u() <= end.u(); ix.next_u(), i++ )
                    {
                        (*($self)).set_data(ix.coord(), numpy_3d_in[i]);
                    }
            return i;
        }
    
        RTop<ftype> operator_orth_grid () { return (*($self)).operator_orth_grid(); }
        RTop<ftype> operator_grid_orth () { return (*($self)).operator_grid_orth(); }
    
    }
    
    %extend NXmap<double>
    {
        NXmap ( clipper::Grid const & grid, clipper::RTop<double> const & rtop)
        {
            return new NXmap<double>( grid, rtop );
        }
    
        int export_numpy ( double *numpy_array, int nu, int nv, int nw )
        {
            int i = 0;
            int top_u, top_v, top_w;
        
            clipper::Coord_grid c;
            clipper::Grid map_grid = (*($self)).grid();
        
            nu > map_grid.nu() ? top_u = map_grid.nu() : top_u = nu;
            nv > map_grid.nv() ? top_v = map_grid.nv() : top_v = nv;
            nw > map_grid.nw() ? top_w = map_grid.nw() : top_w = nw;
        
            for ( c.w() = 0; c.w() < top_w; c.w()++ )
                for ( c.v() = 0; c.v() < top_v; c.v()++ )
                    for ( c.u() = 0; c.u() < nu; c.u()++, i++ )
                    {
                        if ( c.u() < map_grid.nu() && c.v() < map_grid.nv() && c.w() < map_grid.nw() )
                            numpy_array[i] = (*($self)).get_data(c);
                        else
                            numpy_array[i] = 0.0;
                    }
            return i;
        }
    
        int import_numpy ( double *numpy_3d_in, int nu, int nv, int nw )
        {
            int i = 0;
            int top_u, top_v, top_w;
        
            clipper::Coord_grid c;
            clipper::Grid map_grid = (*($self)).grid();
        
            nu > map_grid.nu() ? top_u = map_grid.nu() : top_u = nu;
            nv > map_grid.nv() ? top_v = map_grid.nv() : top_v = nv;
            nw > map_grid.nw() ? top_w = map_grid.nw() : top_w = nw;
        
            for ( c.w() = 0; c.w() < top_w; c.w()++ )
                for ( c.v() = 0; c.v() < top_v; c.v()++ )
                    for (  c.u() = 0; c.u() < top_u; c.u()++, i++ )
                        (*($self)).set_data(c, numpy_3d_in[i]);
            return i;
        }
            
        RTop<double> operator_orth_grid () { return (*($self)).operator_orth_grid(); }
        RTop<double> operator_grid_orth () { return (*($self)).operator_grid_orth(); }
    }
    
    %extend Xmap<double>
    {
        int export_numpy ( double *numpy_array, int nu, int nv, int nw )
        {
            int i = 0;
            int top_u, top_v, top_w;
        
            clipper::Coord_grid c;
            clipper::Grid map_grid = (*($self)).grid_asu();
        
            for ( c.w() = 0; c.w() < nw; c.w()++ )
                for ( c.v() = 0; c.v() < nv; c.v()++ )
                    for ( c.u() = 0; c.u() < nu; c.u()++, i++ )
                    {
                        if ( c.u() < map_grid.nu() && c.v() < map_grid.nv() && c.w() < map_grid.nw() )
                            numpy_array[i] = (*($self)).get_data(c);
                        else
                            numpy_array[i] = 0.0;
                    }
            return i;
        }
    
        int import_numpy ( double *numpy_3d_in, int nu, int nv, int nw )
        {
            int i = 0;
            int top_u, top_v, top_w;
        
            clipper::Coord_grid c;
            clipper::Grid map_grid = (*($self)).grid_asu();
        
            nu > map_grid.nu() ? top_u = map_grid.nu() : top_u = nu;
            nv > map_grid.nv() ? top_v = map_grid.nv() : top_v = nv;
            nw > map_grid.nw() ? top_w = map_grid.nw() : top_w = nw;
        
            for ( c.w() = 0; c.w() < top_w; c.w()++ )
                for ( c.v() = 0; c.v() < top_v; c.v()++ )
                    for (  c.u() = 0; c.u() < top_u; c.u()++, i++ )
                        (*($self)).set_data(c, numpy_3d_in[i]);
            return i;
        }
    
        int export_section_numpy ( double *numpy_array, int nu, int nv, int nw, clipper::Coord_grid& start, clipper::Coord_grid& end )
        {
            int i = 0;
            int w, v;
            clipper::Xmap_base::Map_reference_coord ix( (*($self)) );
        
            for ( w = start.w(); w <= end.w(); w++ )
            {
                for ( v = start.v(); v <= end.v(); v++ )
                    for ( ix.set_coord(Coord_grid(start.u(),v,w)); ix.coord().u() <= end.u(); ix.next_u(), i++ )
                        numpy_array[i] = (*($self)).get_data(ix.coord());
            }
            return i;
        }
    
        int import_section_numpy ( double *numpy_3d_in, int nu, int nv, int nw, clipper::Coord_grid& start, clipper::Coord_grid& end )
        {
            int i = 0;
            int w, v;
            clipper::Xmap_base::Map_reference_coord ix( (*($self)) );
        
            for ( w = start.w(); w <= end.w(); w++ )
            {
                for ( v = start.v(); v <= end.v(); v++ )
                    for ( ix.set_coord(Coord_grid(start.u(),v,w)); ix.coord().u() <= end.u(); ix.next_u(), i++ )
                    {
                        (*($self)).set_data(ix.coord(), numpy_3d_in[i]);
                    }
            }
            return i;
        }
    
        RTop<double> operator_orth_grid () { return (*($self)).operator_orth_grid(); }
        RTop<double> operator_grid_orth () { return (*($self)).operator_grid_orth(); }
    }
}


namespace clipper {
    class Map_reference_base
    {
    public:
      inline const NXmap_base& base_nxmap() const { return *map_; }
      inline const int& index() const { return index_; }
      inline bool last() const { return ( index_ >= map_->grid_.size() ); }
    protected:
      const NXmap_base* map_;
      int index_;
    };
    
    
    class Map_reference_index : public Map_reference_base
    {
    public:
      Map_reference_index() {}
      explicit Map_reference_index( const NXmap_base& map )
        { map_ = &map; index_ = 0; }
      Map_reference_index( const NXmap_base& map, const Coord_grid& pos )
        { map_ = &map; index_ = map_->grid_.index( pos ); }
      inline Coord_grid coord() const
        { return map_->grid_.deindex(index_); }
      inline const Coord_orth coord_orth() const
        { return map_->coord_orth( coord().coord_map() ); }
      inline Map_reference_index& set_coord( const Coord_grid& pos )
        { index_ = map_->grid_.index( pos ); return *this; }
      inline Map_reference_index& next() { index_++; return *this; }
      inline int index_offset(const int& du,const int& dv,const int& dw) const {
        return index_ + du*map_->du + dv*map_->dv + dw*map_->dw;
      }
    };
    
    class Map_reference_coord : public Map_reference_base
    {
    public:
      Map_reference_coord() {}
      explicit Map_reference_coord( const NXmap_base& map )
        { map_ = &map; }
      Map_reference_coord( const NXmap_base& map, const Coord_grid& pos )
        { map_ = &map; set_coord( pos ); }
      inline Coord_grid coord() const { return pos_; }
      inline const Coord_orth coord_orth() const
        { return map_->coord_orth( coord().coord_map() ); }
      inline Map_reference_coord& set_coord( const Coord_grid& pos )
        { pos_ = pos; index_ = map_->grid_.index( pos_ ); return *this; }
      inline Map_reference_coord& next() {
        index_++;
        pos_ = map_->grid_.deindex(index_);
        return *this;
      }
      inline Map_reference_coord& next_u() { pos_.u()++; index_ += map_->du; return *this; }
      inline Map_reference_coord& next_v() { pos_.v()++; index_ += map_->dv; return *this; }
      inline Map_reference_coord& next_w() { pos_.w()++; index_ += map_->dw; return *this; }
      inline Map_reference_coord& prev_u() { pos_.u()--; index_ -= map_->du; return *this; }
      inline Map_reference_coord& prev_v() { pos_.v()--; index_ -= map_->dv; return *this; }
      inline Map_reference_coord& prev_w() { pos_.w()--; index_ -= map_->dw; return *this; }
      inline Map_reference_coord& operator =( const Coord_grid& pos )
        { return set_coord( pos ); }
    protected:
      Coord_grid pos_;
    };
}

%{
  namespace clipper {
    typedef NXmap_base::Map_reference_base Map_reference_base;
    typedef NXmap_base::Map_reference_index Map_reference_index;
    typedef NXmap_base::Map_reference_coord Map_reference_coord;
  }
%}


namespace clipper {
  %extend Xmap<float>
  {
    void fft_from (const clipper::HKL_data<clipper::data32::F_phi> &fb)
    {
        ($self)->fft_from( fb );
    }
    
    void fft_to (clipper::HKL_data<clipper::data32::F_phi> &fphidata)
    {
        ($self)->fft_to(fphidata, clipper::Xmap_base::Default);
    }
    
  }
  %extend Xmap<double>
  {
    void fft_from (const clipper::HKL_data<clipper::data64::F_phi> &fb)
    {
        ($self)->fft_from( fb );
    }
    
    void fft_to (clipper::HKL_data<clipper::data32::F_phi> &fphidata ) const
    {
        ($self)->fft_to(fphidata, clipper::Xmap_base::Default);
    }
  
  }
  
  %template(Xmap_float) Xmap<float>;
  %template(Xmap_double) Xmap<double>;
  %template(Xmap_int) Xmap<int>;
  %template(NXmap_float) NXmap<float>;
  %template(NXmap_double) NXmap<double>;
  %template(NXmap_int) NXmap<int>;
}

%include "../clipper/ccp4/ccp4_map_io.h"
namespace clipper {
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

%include "../clipper/ccp4/ccp4_utils.h"

/*
%include "../clipper/cns/cns_map_io.h"
%include "../clipper/cns/cns_hkl_io.h"
namespace clipper {
  %extend CNSMAPfile {
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
*/

%include "../clipper/core/clipper_stats.h"

namespace clipper {
  %template(Range_float) Range<float>;
  %template(Range_double) Range<double>;
}

//FIXME  - We do not really want to be producing new objects. Wrapping constructor properly would be preferred.
//         But SWIG does not let me.
namespace clipper {
  %extend Map_stats {
     Map_stats(const Xmap<float> &m){
       return new Map_stats(m);
     };
     Map_stats(const Xmap<double> &m){
       return new Map_stats(m);
     };
     Range<double> range_double() {
       return ($self)->range();
     }
     Range<double> range_float() {
       return ($self)->range();
     }
  };
}

%include "../clipper/core/map_utils.h"


namespace clipper {
  %extend HKL {
    HKL __add__(const HKL &h2){
      HKL ret;
      ret = *($self)+h2;
      return ret;
    }
    HKL __sub__(const HKL &h2){
      HKL ret;
      ret = *($self)-h2;
      return ret;
    }
    HKL __neg__(){
      HKL ret;
      ret = -*($self);
      return ret;
    }
  }
  
/*
%rename(equals_Coord_grid)    operator== (const Coord_grid&, const Coord_grid&);
%rename(notequals_Coord_grid) operator!= (const Coord_grid&, const Coord_grid&); 
%rename(neg_Coord_grid)       operator-  (const Coord_grid&);
%rename(add_Coord_grid)       operator+  (const Coord_grid&, const Coord_grid&);
%rename(subs_Coord_grid)      operator-  (const Coord_grid&, const Coord_grid&);
%rename(product_Coord_grid)   operator*  (const int&, const Coord_grid&);
%rename(transf_Coord_grid)    operator*  (const Isymop&, const Coord_grid&);
*/

  %extend Coord_grid {
    bool __cmp__(const Coord_grid& g2){
      bool ret;
      ret = (*($self) == g2);
      return ret;
    }
    bool __ne__(const Coord_grid& g2){
      bool ret;
      ret = (*($self) != g2);
      return ret;
    }
    Coord_grid __add__(const Coord_grid& g2){
      Coord_grid ret;
      ret = *($self) + g2;
      return ret;
    }
    Coord_grid __neg__( ){
      Coord_grid ret;
      ret = -*($self);
      return ret;
    }
    Coord_grid __sub__(const Coord_grid& g2){
      Coord_grid ret;
      ret = *($self) - g2;
      return ret;
    }
    void uvw(int numpy_int_out[3]) {
      for (int i = 0; i < 3; i++) {
        numpy_int_out[i] = ((*($self))[i]);
      }
    }
      
    /*Coord_grid __mul__(const int& i1){
      Coord_grid ret;
      ret = *($self) * i1;
      return ret;
    }
    Coord_grid __mul__(const Isymop& isym){
      Coord_grid ret;
      ret = *($self) * isym;
      return ret;
    }*/
  }
  
  %extend Coord_orth {
    Coord_orth __add__(const Coord_orth &h2){
      Coord_orth ret;
      ret = *($self)+h2;
      return ret;
    }
    Coord_orth __sub__(const Coord_orth &h2){
      Coord_orth ret;
      ret = *($self)-h2;
      return ret;
    }
    Coord_orth __neg__(){
      Coord_orth ret;
      ret = -*($self);
      return ret;
    }
    void xyz(double numpy_double_out[3]) {
      for (int i = 0; i < 3; i++) {
        numpy_double_out[i] = (*($self))[i];
      }      
    }
    Coord_orth __mul__ ( const float &factor ) 
    {
      return factor * (*self);
    }
    Coord_orth __mul__ ( const RTop_orth& op )
    {
      return op * (*self);
    }
    Coord_orth __rmul__ ( const float &factor ) 
    {
      return factor * (*self);
    }
    Coord_orth __rmul__ ( const RTop_orth& op )
    {
      return op * (*self);
    }

  }
  %extend Coord_frac {
    Coord_frac __add__(const Coord_frac &h2){
      Coord_frac ret;
      ret = *($self)+h2;
      return ret;
    }
    Coord_frac __sub__(const Coord_frac &h2){
      Coord_frac ret;
      ret = *($self)-h2;
      return ret;
    }
    Coord_frac __neg__(){
      Coord_frac ret;
      ret = -*($self);
      return ret;
    }
    void uvw(double numpy_double_out[3]) {
      for (int i = 0; i < 3; i++) {
        numpy_double_out[i] = (*($self))[i];
      }      
    }
    Coord_frac __mul__ ( const float &factor ) 
    {
      return factor * (*self);
    }
    Coord_frac __mul__ ( const RTop_frac& op )
    {
      return op * (*self);
    }
    Coord_frac __rmul__ ( const float &factor ) 
    {
      return factor * (*self);
    }
    Coord_frac __rmul__ ( const RTop_frac& op )
    {
      return op * (*self);
    }
  }

CATCH_SWIG_EXCEPTION(Atom_list::__getitem__)
CATCH_SWIG_EXCEPTION(Atom_list::__setitem__)

CATCH_SWIG_EXCEPTION(Atom_list::set_elements)
CATCH_SWIG_EXCEPTION(Atom_list::set_coord_orth)
CATCH_SWIG_EXCEPTION(Atom_list::set_occupancies)
CATCH_SWIG_EXCEPTION(Atom_list::set_u_isos)
CATCH_SWIG_EXCEPTION(Atom_list::set_u_anisos)



  
  
  %extend Atom_list {
    Atom& __getitem__(int i) {
      int array_len = $self->size(); 
      i = (i < 0) ? array_len + i : i;
      if (i >= array_len || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;        
        return (*($self))[0];
      }
      return (*($self))[i];      
      fail:
        return (*($self))[0];
    }
    void __setitem__(int i, Atom& atom) {
      int array_len = $self->size();
      i = (i < 0) ? array_len + i : i;
      if (i >= array_len || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return;
      }
      (*($self))[i] = atom;
      return;
      fail:
        return;
    }
    size_t __len__() { 
      return ($self)->size();
    }
    
    Atom pop(int i) {
      size_t array_len = $self->size();
      if (i >= array_len || -i > array_len) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return (*($self))[0];
      }
      Atom ret;
      if (i>=0) ret = (*($self))[i];
      else ret = (*($self))[array_len - i];
      $self -> erase(($self) -> begin() + i);
      return ret;
    }
    size_t append (Atom& a) {
      ($self)->push_back(a);
      return ($self)->size();
    }
    void extend_by(size_t n) {
      for (size_t i = 0; i < n; i++ ) {
        ($self)->push_back( Atom() );
      }
    }
    void set_elements(std::vector<std::string> elements) {
      size_t in_len = elements.size();
      size_t my_len = ($self) -> size();
      if (in_len != my_len) {
        myErr = 1;
        myErrType = ARRAY_LENGTH_MISMATCH;
        return;
      }
      for (size_t i = 0; i < my_len; i++) {
        (*($self))[i].set_element( elements[i] );
      }
    }
    void set_coord_orth(double *numpy_2d_in, int n1, int n2) {
      size_t my_len = ($self) -> size();
      if (n1 != my_len) {
        myErr = 1;
        myErrType = ARRAY_LENGTH_MISMATCH;
        return;
      }
      if (n2 != 3) {
        myErr = 1;
        myErrType = CUSTOM_INDEX_ERROR;
        myErrString = "Coordinates should be in the form of an N x 3 array";
        return;
      }
      for (int i = 0; i < n1; i++) {
        int first = i * n2;
          (*($self))[i].set_coord_orth( Coord_orth( numpy_2d_in[first], numpy_2d_in[first+1], numpy_2d_in[first+2] ) );
      }
      return;
    }
    void set_occupancies(double *numpy_1d_in, int n) {
      size_t my_len = ($self) -> size();
      if (n != my_len) {
        myErr = 1;
        myErrType = ARRAY_LENGTH_MISMATCH;
        return;
      }
      for (int i = 0; i < n; i++) {
        (*($self))[i].set_occupancy( numpy_1d_in[i] );
      }
      return;
    }
    void set_u_isos(double *numpy_1d_in, int n) {
      size_t my_len = ($self) -> size();
      if (n != my_len) {
        myErr = 1;
        myErrType = ARRAY_LENGTH_MISMATCH;
        return;
      }
      for (int i = 0; i < n; i++) {
        (*($self))[i].set_u_iso( numpy_1d_in[i] );
      }
      return;
    }
    void set_u_anisos(double *numpy_2d_in, int n1, int n2) {
      size_t my_len = ($self) -> size();
      if (n1 != my_len) {
        myErr = 1;
        myErrType = ARRAY_LENGTH_MISMATCH;
        return;
      }
      if (n2 != 6) {
        myErr = 1;
        myErrType = CUSTOM_INDEX_ERROR;
        myErrString = "Input should be in the form of a N x 6 matrix, where"
                      "each row is of the form [u11, u22, u33, u12, u13, u23]";
        return;
      }
      double *ai = numpy_2d_in;
      for (int i = 0; i < n1; i++) {
        int first = i*n2;
        (*($self))[i].set_u_aniso_orth(
          U_aniso_orth(ai[first], ai[first+1], ai[first+2],
                       ai[first+3], ai[first+4], ai[first+5]) );
      }
      return;
    }
      
    
      
      
      
      


    
    
  }
  
CATCH_SWIG_EXCEPTION(MModel::__getitem__)
CATCH_SWIG_EXCEPTION(MModel::__setitem__)

  %extend MModel {
    MPolymer& __getitem__(int i) { 
      int array_len = $self->size();
      i = (i < 0) ? array_len + i : i;
      if (i >= array_len || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return (*($self))[0];
      }
      return (*($self))[i];
      fail:
        return (*($self))[0];
    }
    void __setitem__(int i, MPolymer& mpol) {
      int array_len = $self->size();
      i = (i < 0) ? array_len + i : i;
      if (i >= array_len || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return;
      }
      (*($self))[i]=mpol;
      return;
      fail:
        return;
    }
    size_t __len__() { 
      return ($self)->size();
    }
  }
  
CATCH_SWIG_EXCEPTION(MPolymer::__getitem__)
CATCH_SWIG_EXCEPTION(MPolymer::__setitem__)

  %extend MPolymer {
    MMonomer& __getitem__(int i) { 
      int array_len = $self->size();
      i = (i < 0) ? array_len + i : i;
      if (i >= array_len || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return (*($self))[0];
      }
      return (*($self))[i];
      fail:
        return (*($self))[0];
    }
    void __setitem__(int i, MMonomer& mmon) {
      int array_len = $self->size();      
      i = (i < 0) ? array_len + i : i;
      if (i >= array_len || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return;
      }
      (*($self))[i]=mmon;
      return;
      fail:
        return;
    }
    size_t __len__() { 
      return ($self)->size();
    }
  }
  
  CATCH_SWIG_EXCEPTION(MMonomer::__getitem__)
  CATCH_SWIG_EXCEPTION(MMonomer::__setitem__)

  %extend MMonomer {
    MAtom& __getitem__(int i) { 
      int array_len = $self->size();            
      i = (i < 0) ? array_len + i : i;
      if (i >= array_len || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return (*($self))[0];
      }
      return (*($self))[i];
      fail:
        return (*($self))[0];
    }
    void __setitem__(int i, MAtom& atom) {
      int array_len = $self->size();                  
      i = (i < 0) ? array_len + i : i;
      if (i >= array_len || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return;
      }
      (*($self))[i]=atom;
      return;
      fail:
        return;
    }
    size_t __len__() { 
      return ($self)->size();
    }
  }
  
  CATCH_SWIG_EXCEPTION(MAtom::set_vals)

  %extend MAtom {
    std::string __str__( ) { 
      return (*($self)).id() + " " + (*($self)).coord_orth().format();
      fail:
        return "";
    }
  }

  %extend U_aniso_orth {
    std::string __str__( ) { 
      return (*($self)).format();
      fail:
        return "";
    }
    void get_vals(double numpy_double_out[6]) {
      numpy_double_out[0] = (*($self)).mat00();
      numpy_double_out[1] = (*($self)).mat11();
      numpy_double_out[2] = (*($self)).mat22();
      numpy_double_out[3] = (*($self)).mat01();
      numpy_double_out[4] = (*($self)).mat02();
      numpy_double_out[5] = (*($self)).mat12();
    }    
  }

}

%include "../clipper/mmdb/clipper_mmdb.h"
%include "../clipper/minimol/minimol.h"
%include "../clipper/minimol/minimol_io.h"
%include "../clipper/minimol/minimol_seq.h"


namespace clipper {
  namespace data64 {
    %extend Flag_bool {
      bool get_flag() { 
        bool theFlag = ($self)->flag();
        return theFlag;
      }
      void set_flag(bool theFlag) { 
        ($self)->flag() = theFlag;
      }
    }
    %extend Flag {
      int get_flag() { 
        int theFlag = ($self)->flag();
        return theFlag;
      }
      void set_flag(int theFlag) { 
        ($self)->flag() = theFlag;
      }
    }
  }
  namespace data32 {
    %extend Flag_bool {
      bool get_flag() { 
        bool theFlag = ($self)->flag();
        return theFlag;
      }
      void set_flag(bool theFlag) { 
        ($self)->flag() = theFlag;
      }
    }
    %extend Flag {
      int get_flag() { 
        int theFlag = ($self)->flag();
        return theFlag;
      }
      void set_flag(int theFlag) { 
        ($self)->flag() = theFlag;
      }
    }
  }
  namespace datatypes {
    %extend Flag_bool {
      bool get_flag() { 
        bool theFlag = ($self)->flag();
        return theFlag;
      }
      void set_flag(bool theFlag) { 
        ($self)->flag() = theFlag;
      }
      clipper::datatypes::Flag_bool copy(){
        clipper::datatypes::Flag_bool ret;
        ret = *($self);
        return ret;
      }
    }
    %extend Flag {
      int get_flag() { 
        int theFlag = ($self)->flag();
        return theFlag;
      }
      void set_flag(int theFlag) { 
        ($self)->flag() = theFlag;
      }
      clipper::datatypes::Flag copy(){
        clipper::datatypes::Flag ret;
        ret = *($self);
        return ret;
      }
    }
  }
}


//%rename (to_complex_float) operator complex<float>();
//%rename (to_complex_double) operator complex<double>();


%include "../clipper/core/hkl_datatypes.h"

namespace clipper 
{
    
    %template(F_sigF_float) clipper::datatypes::F_sigF<float>;
    %template(F_sigF_double) clipper::datatypes::F_sigF<double>;
    %template(HKL_data_F_sigF_float) HKL_data< clipper::datatypes::F_sigF<float> >;
    %template(HKL_data_F_sigF_double) HKL_data< clipper::datatypes::F_sigF<double> >;

    %template(F_sigF_ano_float) clipper::datatypes::F_sigF_ano<float>;
    %template(F_sigF_ano_double) clipper::datatypes::F_sigF_ano<double>;
    %template(HKL_data_F_sigF_ano_float) HKL_data< clipper::data32::F_sigF_ano >;
    %template(HKL_data_F_sigF_ano_double) HKL_data< clipper::data64::F_sigF_ano >;

    %template(I_sigI_float) clipper::datatypes::I_sigI<float>;
    %template(I_sigI_double) clipper::datatypes::I_sigI<double>;
    %template(HKL_data_I_sigI_float) HKL_data< clipper::data32::I_sigI >;
    %template(HKL_data_I_sigI_double) HKL_data< clipper::data64::I_sigI >;

    %template(E_sigE_float) clipper::datatypes::E_sigE<float>;
    %template(E_sigE_double) clipper::datatypes::E_sigE<double>;
    %template(HKL_data_E_sigE_float) HKL_data< clipper::data32::E_sigE >;
    %template(HKL_data_E_sigE_double) HKL_data< clipper::data64::E_sigE >;

    %template(ABCD_float) clipper::datatypes::ABCD<float>;
    %template(ABCD_double) clipper::datatypes::ABCD<double>;
    %template(HKL_data_ABCD_float) HKL_data< clipper::data32::ABCD >;
    %template(HKL_data_ABCD_double) HKL_data< clipper::data64::ABCD >;

    %template(Phi_fom_float) clipper::datatypes::Phi_fom<float>;
    %template(Phi_fom_double) clipper::datatypes::Phi_fom<double>;
    %template(HKL_data_Phi_fom_float) HKL_data< clipper::data32::Phi_fom >;
    %template(HKL_data_Phi_fom_double) HKL_data< clipper::data64::Phi_fom >;

    %template(F_phi_float) clipper::datatypes::F_phi<float>;
    %template(F_phi_double) clipper::datatypes::F_phi<double>;
    %template(HKL_data_F_phi_float) HKL_data< clipper::data32::F_phi >;
    %template(HKL_data_F_phi_double) HKL_data< clipper::data64::F_phi >;

    %template(HKL_data_Flag) HKL_data< clipper::data32::Flag>;
    %template(HKL_data_Flag_bool) HKL_data< clipper::datatypes::Flag_bool>;
    
    CATCH_SWIG_EXCEPTION(HKL_data_Flag::__getitem__)
    CATCH_SWIG_EXCEPTION(HKL_data_Flag_bool::__getitem__)


  %extend HKL_data<clipper::data32::Flag_bool> {
    HKL_data<clipper::datatypes::Flag_bool>  copy(){
      HKL_data<clipper::data32::Flag_bool> ret;
      ret = *($self);
      return ret;
    }
  }

  // Flag works, all others do not (fall over in Python).
  %extend HKL_data<clipper::data32::F_phi> {
    HKL_data<clipper::datatypes::Flag_bool> not_(){
      return !(*($self));
    }
    HKL_data<clipper::datatypes::Flag_bool> __or__(const HKL_data<clipper::datatypes::F_phi<float> > &d1){
      return (*($self)) | d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __xor__(const HKL_data<clipper::datatypes::F_phi<float> > &d1){
      return (*($self)) ^ d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __and__(const HKL_data<clipper::datatypes::F_phi<float> > &d1){
      return (*($self)) & d1;
    }
  }

  %extend HKL_data<clipper::datatypes::F_sigF<float> > {
    HKL_data<clipper::datatypes::Flag_bool> not_(){
      return !(*($self));
    }
    HKL_data<clipper::datatypes::Flag_bool> __or__(const HKL_data<clipper::datatypes::F_sigF<float> > &d1){
      return (*($self)) | d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __xor__(const HKL_data<clipper::datatypes::F_sigF<float> > &d1){
      return (*($self)) ^ d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __and__(const HKL_data<clipper::datatypes::F_sigF<float> > &d1){
      return (*($self)) & d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __or__(const HKL_data<clipper::data32::F_sigF> &d1){
      return (*($self)) | d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __xor__(const HKL_data<clipper::data32::F_sigF> &d1){
      return (*($self)) ^ d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __and__(const HKL_data<clipper::data32::F_sigF> &d1){
      return (*($self)) & d1;
    }
  }

  %extend HKL_data<clipper::data32::F_sigF> {
    HKL_data<clipper::datatypes::Flag_bool> not_(){
      return !(*($self));
    }
    HKL_data<clipper::datatypes::Flag_bool> __or__(const HKL_data<clipper::datatypes::F_sigF<float> > &d1){
      return (*($self)) | d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __xor__(const HKL_data<clipper::datatypes::F_sigF<float> > &d1){
      return (*($self)) ^ d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __and__(const HKL_data<clipper::datatypes::F_sigF<float> > &d1){
      return (*($self)) & d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __or__(const HKL_data<clipper::data32::F_sigF> &d1){
      return (*($self)) | d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __xor__(const HKL_data<clipper::data32::F_sigF> &d1){
      return (*($self)) ^ d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __and__(const HKL_data<clipper::data32::F_sigF> &d1){
      return (*($self)) & d1;
    }
  }

  %extend HKL_data<clipper::data32::F_sigF_ano> {
    HKL_data<clipper::datatypes::Flag_bool> not_(){
      return !(*($self));
    }
    HKL_data<clipper::datatypes::Flag_bool> __or__(const HKL_data<clipper::datatypes::F_sigF_ano<float> > &d1){
      return (*($self)) | d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __xor__(const HKL_data<clipper::datatypes::F_sigF_ano<float> > &d1){
      return (*($self)) ^ d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __and__(const HKL_data<clipper::datatypes::F_sigF_ano<float> > &d1){
      return (*($self)) & d1;
    }
  }

  %extend HKL_data<clipper::data32::I_sigI> {
    HKL_data<clipper::datatypes::Flag_bool> not_(){
      return !(*($self));
    }
    HKL_data<clipper::datatypes::Flag_bool> __or__(const HKL_data<clipper::datatypes::I_sigI<float> > &d1){
      return (*($self)) | d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __xor__(const HKL_data<clipper::datatypes::I_sigI<float> > &d1){
      return (*($self)) ^ d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __and__(const HKL_data<clipper::datatypes::I_sigI<float> > &d1){
      return (*($self)) & d1;
    }
  }

  %extend HKL_data<clipper::data32::E_sigE> {
    HKL_data<clipper::datatypes::Flag_bool> not_(){
      return !(*($self));
    }
    HKL_data<clipper::datatypes::Flag_bool> __or__(const HKL_data<clipper::datatypes::E_sigE<float> > &d1){
      return (*($self)) | d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __xor__(const HKL_data<clipper::datatypes::E_sigE<float> > &d1){
      return (*($self)) ^ d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __and__(const HKL_data<clipper::datatypes::E_sigE<float> > &d1){
      return (*($self)) & d1;
    }
  }

  %extend HKL_data<clipper::data32::ABCD> {
    HKL_data<clipper::datatypes::Flag_bool> not_(){
      return !(*($self));
    }
    HKL_data<clipper::datatypes::Flag_bool> __or__(const HKL_data<clipper::datatypes::ABCD<float> > &d1){
      return (*($self)) | d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __xor__(const HKL_data<clipper::datatypes::ABCD<float> > &d1){
      return (*($self)) ^ d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __and__(const HKL_data<clipper::datatypes::ABCD<float> > &d1){
      return (*($self)) & d1;
    }
  }

  %extend HKL_data<clipper::data32::Phi_fom> {
    HKL_data<clipper::datatypes::Flag_bool> not_(){
      return !(*($self));
    }
    HKL_data<clipper::datatypes::Flag_bool> __or__(const HKL_data<clipper::datatypes::Phi_fom<float> > &d1){
      return (*($self)) | d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __xor__(const HKL_data<clipper::datatypes::Phi_fom<float> > &d1){
      return (*($self)) ^ d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __and__(const HKL_data<clipper::datatypes::Phi_fom<float> > &d1){
      return (*($self)) & d1;
    }
  }

  %extend HKL_data<clipper::data32::Flag> {
    HKL_data<clipper::datatypes::Flag_bool> not_(){
      return !(*($self));
    }
    HKL_data<clipper::datatypes::Flag_bool> __or__(const HKL_data<clipper::data32::Flag> &d1){
      return (*($self)) | d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __xor__(const HKL_data<clipper::data32::Flag> &d1){
      return (*($self)) ^ d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __and__(const HKL_data<clipper::data32::Flag> &d1){
      return (*($self)) & d1;
    }
  }

  %extend HKL_data<clipper::data32::Flag> {
    HKL_data<clipper::datatypes::Flag_bool> __eq__(const int& n){
      return (*($self)) == n;
    }
    HKL_data<clipper::datatypes::Flag_bool> __ne__(const int& n){
      return (*($self)) != n;
    }
    HKL_data<clipper::datatypes::Flag_bool> __ge__(const int& n){
      return (*($self)) >= n;
    }
    HKL_data<clipper::datatypes::Flag_bool> __le__(const int& n){
      return (*($self)) <= n;
    }
    HKL_data<clipper::datatypes::Flag_bool> __gt__(const int& n){
      return (*($self)) > n;
    }
    HKL_data<clipper::datatypes::Flag_bool> __lt__(const int& n){
      return (*($self)) < n;
    }
    HKL_data<clipper::datatypes::Flag>  copy(){
      HKL_data<clipper::data32::Flag> ret;
      ret = *($self);
      return ret;
    }
  }

  %extend datatypes::F_sigF<float> {
    clipper::datatypes::F_sigF<float>  copy(){
      clipper::data32::F_sigF ret;
      ret = *($self);
      return ret;
    }
  }

  %extend datatypes::F_sigF_ano<float> {
    clipper::datatypes::F_sigF_ano<float>  copy(){
      clipper::data32::F_sigF_ano ret;
      ret = *($self);
      return ret;
    }
  }

  %extend datatypes::I_sigI<float> {
    clipper::datatypes::I_sigI<float>  copy(){
      clipper::data32::I_sigI ret;
      ret = *($self);
      return ret;
    }
  }

  %extend datatypes::E_sigE<float> {
    clipper::datatypes::E_sigE<float>  copy(){
      clipper::data32::E_sigE ret;
      ret = *($self);
      return ret;
    }
  }

  %extend datatypes::F_phi<float> {
    clipper::datatypes::F_phi<float>  __add__(const clipper::datatypes::F_phi<float> &h2){
      clipper::data32::F_phi ret;
      ret = *($self)+h2;
      return ret;
    }
    clipper::datatypes::F_phi<float>  __sub__(const clipper::datatypes::F_phi<float> &h2){
      clipper::data32::F_phi ret;
      ret = *($self)-h2;
      return ret;
    }
    clipper::datatypes::F_phi<float>  __neg__(){
      clipper::data32::F_phi ret;
      ret = -*($self);
      return ret;
    }
    clipper::datatypes::F_phi<float>  copy(){
      clipper::data32::F_phi ret;
      ret = *($self);
      return ret;
    }
  }

  %extend datatypes::ABCD<float> {
    clipper::datatypes::ABCD<float>  __add__(const clipper::datatypes::ABCD<float> &h2){
      clipper::data32::ABCD ret;
      ret = *($self)+h2;
      return ret;
    }
    clipper::datatypes::ABCD<float>  copy(){
      clipper::data32::ABCD ret;
      ret = *($self);
      return ret;
    }
  }

  %extend HKL_data<clipper::data32::ABCD> {
    HKL_data<clipper::datatypes::ABCD<float> > __add__(const HKL_data<clipper::datatypes::ABCD<float> > &h2){
      HKL_data<clipper::data32::ABCD> ret;
      ret = *($self)+h2;
      return ret;
    }
    HKL_data<clipper::datatypes::ABCD<float> > copy(){
      HKL_data<clipper::data32::ABCD> ret;
      ret = *($self);
      return ret;
    }
  }

  %extend HKL_data<clipper::data32::F_phi> {
    HKL_data<clipper::datatypes::F_phi<float> >  copy(){
      HKL_data<clipper::data32::F_phi> ret;
      ret = *($self);
      return ret;
    }
    /*
       This would be nice, but what do I do with memo? 
       Python way is:
       memo[id(self)] = result (where result is new class)
       How on earth can I do this in Python?
       But without this method os.deepcopy will never work.
    HKL_data<clipper::datatypes::F_phi<float> >  __deepcopy__(PyObject *memo){
      HKL_data<clipper::data32::F_phi> ret;
      ret = *($self);
      return ret;
    }
    */
    HKL_data<clipper::datatypes::F_phi<float> > __add__(const HKL_data<clipper::datatypes::F_phi<float> > &h2){
      HKL_data<clipper::data32::F_phi> ret;
      ret = *($self)+h2;
      return ret;
    }
    HKL_data<clipper::datatypes::F_phi<float> > __sub__(const HKL_data<clipper::datatypes::F_phi<float> > &h2){
      HKL_data<clipper::datatypes::F_phi<float> > ret;
      ret = *($self)-h2;
      return ret;
    }
    HKL_data<clipper::datatypes::F_phi<float> > __neg__(){
      HKL_data<clipper::datatypes::F_phi<float> > ret;
      ret = -*($self);
      return ret;
    }
    HKL_data<clipper::datatypes::F_phi<float> > __mul__(const float s){
      HKL_data<clipper::datatypes::F_phi<float> > ret;
      ret = *($self)*s;
      return ret;
    }
    HKL_data<clipper::datatypes::F_phi<float> > __rmul__(const float s){
      HKL_data<clipper::datatypes::F_phi<float> > ret;
      ret = *($self)*s;
      return ret;
    }
  }

  %extend HKL_data<clipper::data32::E_sigE> {
    void scaleBySqrtResolution(const clipper::ResolutionFn &escale){
      for ( clipper::HKL_data_base::HKL_reference_index ih = (*($self)).first(); !ih.last(); ih.next() )
        if ( !(*($self))[ih].missing() ) (*($self))[ih].scale( sqrt( escale.f(ih) ) );
    }
    void scaleByResolution(const clipper::ResolutionFn &escale){
      for ( clipper::HKL_data_base::HKL_reference_index ih = (*($self)).first(); !ih.last(); ih.next() )
        if ( !(*($self))[ih].missing() ) (*($self))[ih].scale( escale.f(ih) );
    }
    HKL_data<clipper::datatypes::E_sigE<float> >  copy(){
      HKL_data<clipper::data32::E_sigE> ret;
      ret = *($self);
      return ret;
    }
  }

  %extend HKL_data<clipper::data32::ABCD> {
    void compute_from_phi_fom(const HKL_data< clipper::datatypes::Phi_fom<float> > &phiw){
      ($self)->compute( phiw, clipper::data32::Compute_abcd_from_phifom() );
    }
    void compute_add_abcd(const HKL_data< clipper::datatypes::ABCD<float> > &abcd1, 
                          const HKL_data< clipper::datatypes::ABCD<float> > &abcd2) {
      ($self)->compute( abcd1, abcd2, clipper::data32::Compute_add_abcd() );
    }
  }
  
  %extend HKL_data<clipper::data32::Phi_fom> {
    void compute_from_abcd(const HKL_data< clipper::datatypes::ABCD<float> > &abcd) {
      ($self)->compute( abcd, clipper::data32::Compute_phifom_from_abcd() );
    }
    HKL_data<clipper::datatypes::Phi_fom<float> >  copy(){
      HKL_data<clipper::data32::Phi_fom> ret;
      ret = *($self);
      return ret;
    }
  }

  %extend HKL_data<clipper::data32::F_sigF> {
    void compute_mean_from_fano(const HKL_data< clipper::datatypes::F_sigF_ano<float> > &fano) {
      ($self)->compute( fano, clipper::data32::Compute_mean_fsigf_from_fsigfano() );
    }
    void compute_diff_from_fano(const HKL_data< clipper::datatypes::F_sigF_ano<float> > &fano) {
      ($self)->compute( fano, clipper::data32::Compute_diff_fsigf_from_fsigfano() );
    }
    void compute_scale_u_iso_fsigf(float scale, float u_value, 
                                  const HKL_data< clipper::datatypes::F_sigF<float> > &fsigf ) {
      ($self)->compute( fsigf, clipper::data32::Compute_scale_u_iso_fsigf(scale, u_value) );
    }
    void compute_scale_u_aniso_fsigf(float scale, clipper::U_aniso_orth u_value, 
                                  const HKL_data< clipper::datatypes::F_sigF<float> > &fsigf ) {
      ($self)->compute( fsigf, clipper::data32::Compute_scale_u_aniso_fsigf(scale, u_value) );
    }
    HKL_data<clipper::datatypes::F_sigF<float> >  copy(){
      HKL_data<clipper::data32::F_sigF> ret;
      ret = *($self);
      return ret;
    }
  }

  %extend HKL_data<clipper::data32::F_sigF_ano> {
    void compute_scale_u_iso_fsigfano(float scale, float u_value, 
                                  const HKL_data< clipper::datatypes::F_sigF_ano<float> > &fsigfano ) {
      ($self)->compute( fsigfano, clipper::data32::Compute_scale_u_iso_fsigfano(scale, u_value) );
    }
    void compute_scale_u_aniso_fsigfano(float scale, clipper::U_aniso_orth u_value, 
                                  const HKL_data< clipper::datatypes::F_sigF_ano<float> > &fsigfano ) {
      ($self)->compute( fsigfano, clipper::data32::Compute_scale_u_aniso_fsigfano(scale, u_value) );
    }
    HKL_data<clipper::datatypes::F_sigF_ano<float> >  copy(){
      HKL_data<clipper::data32::F_sigF_ano> ret;
      ret = *($self);
      return ret;
    }
  }

  %extend HKL_data<clipper::data32::I_sigI> {
    void compute_scale_u_iso_isigi(float scale, float u_value, 
                                  const HKL_data< clipper::datatypes::I_sigI<float> > &isigi ) {
      ($self)->compute( isigi, clipper::data32::Compute_scale_u_iso_isigi(scale, u_value) );
    }
    void compute_scale_u_aniso_isigi(float scale, clipper::U_aniso_orth u_value, 
                                  const HKL_data< clipper::datatypes::I_sigI<float> > &isigi ) {
      ($self)->compute( isigi, clipper::data32::Compute_scale_u_aniso_isigi(scale, u_value) );
    }
    HKL_data<clipper::datatypes::I_sigI<float> >  copy(){
      HKL_data<clipper::data32::I_sigI> ret;
      ret = *($self);
      return ret;
    }
  }
        
  

  %extend HKL_data< clipper::datatypes::F_phi<float> > {

    void getDataNumpy(double *test_numpy_a, int test_numpy_n) { 
        int i=0;
        for(clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next(),i++ ){
            if(!((*($self))[ih].missing())){
                std::vector<xtype> thisData((*($self)).data_size());
                (*($self)).data_export(ih.hkl(),&(thisData[0]));
                std::vector<float> thisDataf((*($self)).data_size());
                for(unsigned idat=0;idat<(*($self)).data_size();++idat) {
                    test_numpy_a[i*(*($self)).data_size()+idat] = thisData[idat];
                }
            } else {
                for(unsigned idat=0;idat<(*($self)).data_size();++idat) {
                    test_numpy_a[i*(*($self)).data_size()+idat] = std::numeric_limits<float>::quiet_NaN();
                }
            }
        }
    }
      
      std::vector<std::vector<float> > getData() { 
      std::vector<std::vector<float> > allData;
      for(clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
          if(!((*($self))[ih].missing())){
              std::vector<xtype> thisData((*($self)).data_size());
              (*($self)).data_export(ih.hkl(),&(thisData[0]));
              std::vector<float> thisDataf((*($self)).data_size());
              for(unsigned idat=0;idat<(*($self)).data_size();++idat) {thisDataf[idat] = thisData[idat];}
              allData.push_back(thisDataf);
          }
      }
      return allData;
    }
    clipper::datatypes::F_phi<float>& __getitem__(int i) { 
      int sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      i = (i < 0) ? sz + i : i;
      if (i >= sz || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return (*($self))[0];
      }
      return (*($self))[i];
      fail:
        return (*($self))[0];
    }
    size_t __len__() { 
      size_t sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      return sz;
    }
  }

  %extend HKL_data<clipper::data32::F_phi> {
    void compute_neg(const HKL_data< clipper::datatypes::F_phi<float> > &fphi ) {
      ($self)->compute( fphi, clipper::data32::Compute_neg_fphi() );
    }
    void compute_add_fphi(const HKL_data< clipper::datatypes::F_phi<float> > &fphi1, 
                          const HKL_data< clipper::datatypes::F_phi<float> > &fphi2) {
      ($self)->compute( fphi1, fphi2, clipper::data32::Compute_add_fphi() );
    }
    void compute_sub_fphi(const HKL_data< clipper::datatypes::F_phi<float> > &fphi1, 
                          const HKL_data< clipper::datatypes::F_phi<float> > &fphi2) {
      ($self)->compute( fphi1, fphi2, clipper::data32::Compute_sub_fphi() );
    }
    void compute_from_fsigf_phifom(const HKL_data< clipper::datatypes::F_sigF<float> > &fsigf,
                                   const HKL_data< clipper::datatypes::Phi_fom<float> > &phifom ) {
      ($self)->compute( fsigf, phifom, clipper::data32::Compute_fphi_from_fsigf_phifom() );
    }
    void compute_scale_u_iso_fphi(float scale, float u_value, 
                                  const HKL_data< clipper::datatypes::F_phi<float> > &fphi ) {
      ($self)->compute( fphi, clipper::data32::Compute_scale_u_iso_fphi(scale, u_value) );
    }
    void compute_scale_u_aniso_fphi(float scale, clipper::U_aniso_orth u_value, 
                                  const HKL_data< clipper::datatypes::F_phi<float> > &fphi ) {
      ($self)->compute( fphi, clipper::data32::Compute_scale_u_aniso_fphi(scale, u_value) );
    }

    void getDataNumpy(double *test_numpy_a, int test_numpy_n) { 
        int i=0;
        for(clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next(),i++ ){
            if(!((*($self))[ih].missing())){
                std::vector<xtype> thisData((*($self)).data_size());
                (*($self)).data_export(ih.hkl(),&(thisData[0]));
                std::vector<float> thisDataf((*($self)).data_size());
                for(unsigned idat=0;idat<(*($self)).data_size();++idat) {
                    test_numpy_a[i*(*($self)).data_size()+idat] = thisData[idat];
                }
            } else {
                for(unsigned idat=0;idat<(*($self)).data_size();++idat) {
                    test_numpy_a[i*(*($self)).data_size()+idat] = std::numeric_limits<float>::quiet_NaN();
                }
            }
        }
    }
      
      std::vector<std::vector<float> > getData() { 
      std::vector<std::vector<float> > allData;
      for(clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
          if(!((*($self))[ih].missing())){
              std::vector<xtype> thisData((*($self)).data_size());
              (*($self)).data_export(ih.hkl(),&(thisData[0]));
              std::vector<float> thisDataf((*($self)).data_size());
              for(unsigned idat=0;idat<(*($self)).data_size();++idat) {thisDataf[idat] = thisData[idat];}
              allData.push_back(thisDataf);
          }
      }
      return allData;
    }
    clipper::data32::F_phi& __getitem__(int i) { 
      int sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      i = (i < 0) ? sz + i : i;
      if (i >= sz || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return (*($self))[0];
      }
      return (*($self))[i];
      fail:
        return (*($self))[0];
    }
    size_t __len__() { 
      size_t sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      return sz;
    }
  }
  
  %extend HKL_data<clipper::data32::E_sigE> {
    void compute_from_fsigf(const HKL_data< clipper::datatypes::F_sigF<float> > &fsigf ) {
      ($self)->compute( fsigf, clipper::data32::Compute_EsigE_from_FsigF() );
    }
  }

  %extend HKL_data<clipper::datatypes::Flag_bool> {
    clipper::datatypes::Flag_bool& __getitem__(int i) { 
      int sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      i = (i < 0) ? sz + i : i;
      if (i >= sz || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return (*($self))[0];
      }
      return (*($self))[i];
      fail:
        return (*($self))[0];
    }
    size_t __len__() { 
      size_t sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      return sz;
    }
  }
  %extend HKL_data<clipper::data32::Flag_bool> {
    clipper::data32::Flag_bool& __getitem__(int i) { 
      int sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      i = (i < 0) ? sz + i : i;
      if (i >= sz || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return (*($self))[0];
      }
      return (*($self))[i];
      fail:
        return (*($self))[0];
    }
    size_t __len__() { 
      size_t sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      return sz;
    }
  }
 
  %extend HKL_data<clipper::datatypes::Flag> {
    clipper::datatypes::Flag& __getitem__(int i) { 
      int sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      i = (i < 0) ? sz + i : i;
      if (i >= sz || i < 0) {
        myErr = 1;
        myErrType = INDEX_OUT_OF_BOUNDS;
        return (*($self))[0];
      }
      return (*($self))[i];
      fail:
        return (*($self))[0];
    }
    size_t __len__() { 
      size_t sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      return sz;
    }
  }
  %extend HKL_data<clipper::data32::Flag> {
    clipper::data32::Flag& __getitem__(int i) { 
      int sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      i = (i < 0) ? sz + i : i;
      if (i >= sz || i < 0) {
        myErr = 1;
        return (*($self))[0];
      }
      return (*($self))[i];
      fail:
        return (*($self))[0];
    }
    size_t __len__() { 
      size_t sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      return sz;
    }
  }

  CATCH_SWIG_EXCEPTION(HKL_data< clipper::datatypes::F_sigF<float> >::__getitem__)

  %extend HKL_data< clipper::datatypes::F_sigF<float> > {

    void getDataNumpy(double *test_numpy_a, int test_numpy_n) { 
        int i=0;
        for(clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next(),i++ ){
            if(!((*($self))[ih].missing())){
                std::vector<xtype> thisData((*($self)).data_size());
                (*($self)).data_export(ih.hkl(),&(thisData[0]));
                std::vector<float> thisDataf((*($self)).data_size());
                for(unsigned idat=0;idat<(*($self)).data_size();++idat) {
                    test_numpy_a[i*(*($self)).data_size()+idat] = thisData[idat];
                }
            } else {
                for(unsigned idat=0;idat<(*($self)).data_size();++idat) {
                    test_numpy_a[i*(*($self)).data_size()+idat] = std::numeric_limits<float>::quiet_NaN();
                }
            }
        }
    }
      
      std::vector<std::vector<float> > getData() { 
      std::vector<std::vector<float> > allData;
      for(clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
          if(!((*($self))[ih].missing())){
              std::vector<xtype> thisData((*($self)).data_size());
              (*($self)).data_export(ih.hkl(),&(thisData[0]));
              std::vector<float> thisDataf((*($self)).data_size());
              for(unsigned idat=0;idat<(*($self)).data_size();++idat) {thisDataf[idat] = thisData[idat];}
              allData.push_back(thisDataf);
          }
      }
      return allData;
    }
    clipper::datatypes::F_sigF<float>& __getitem__(int i) { 
      int sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      i = (i < 0) ? sz + i : i;
      if (i >= sz || i < 0) {
        myErr = 1;
        return (*($self))[0];
      }
      return (*($self))[i];
      fail:
        return (*($self))[0];
    }
    size_t __len__() { 
      size_t sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      return sz;
    }
  }

CATCH_SWIG_EXCEPTION(HKL_data< clipper::data32::F_sigF<float> >::__getitem__)

  %extend HKL_data< clipper::data32::F_sigF<float> > {
    clipper::data32::F_sigF<float>& __getitem__(int i) { 
      int sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      i = (i < 0) ? sz + i : i;
      if (i >= sz || i < 0) {
        myErr = 1;
        return (*($self))[0];
      }
      return (*($self))[i];
      fail:
        return (*($self))[0];
    }
    size_t __len__() { 
      size_t sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      return sz;
    }
  }

}

%include "../clipper/core/ramachandran.h"

%include "../clipper/cif/cif_data_io.h"
%include "../clipper/contrib/sfcalc_obs.h"

%include "../clipper/core/hkl_compute.h"
namespace clipper
{
    %template(SFcalc_obs_bulk_float) SFcalc_obs_bulk<float>;    
}

%include "../clipper/contrib/function_object_bases.h"
%include "../clipper/contrib/sfscale.h"
namespace clipper {
  %template(SFweight_base_float) SFweight_base<float>;
  %template(SFcalc_base_float) SFcalc_base<float>;
  %template(EDcalc_base_float) EDcalc_base<float>;
  %template(SFcalc_obs_base_float) SFcalc_obs_base<float>;
  %template(SFscale_base_float) SFscale_base<float>;
  %template(MapFilter_base_float) MapFilter_base<float>;
  %template(Convolution_search_base_float) Convolution_search_base<float>;
  %template(FFFear_base_float) FFFear_base<float>;
  %template(Skeleton_base_float) Skeleton_base<float, float>;
  %template(OriginMatch_base_float) OriginMatch_base<float>;
  %template (SFscale_aniso_float)  SFscale_aniso<float>;
}
%include "../clipper/contrib/sfweight.h"
namespace clipper {
  %template(SFweight_spline_float) SFweight_spline<float>;
}

%include "../clipper/contrib/sfcalc.h"
namespace clipper {
  %template(SFcalc_iso_sum_float) SFcalc_iso_sum<float>;
  %template(SFcalc_aniso_sum_float) SFcalc_aniso_sum<float>;
  %template(SFcalc_iso_fft_float) SFcalc_iso_fft<float>;
  %template(SFcalc_aniso_fft_float) SFcalc_aniso_fft<float>;
}

%{

template <typename T> void CopyOverHKLInfo(const T &d_in, T &d_out,  const clipper::HKL_info &newhkl){
  HKL_info::HKL_reference_index ih;
  for ( ih = newhkl.first(); !ih.last(); ih.next() ) {
    d_out[ih] = d_in[ih.hkl()];
  }
}

void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::F_sigF> &d_in, clipper::HKL_data<clipper::data32::F_sigF> &d_out,  const clipper::HKL_info &newhkl){
    CopyOverHKLInfo(d_in,d_out,newhkl);
}
void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::F_sigF_ano> &d_in, clipper::HKL_data<clipper::data32::F_sigF_ano> &d_out,  const clipper::HKL_info &newhkl){
    CopyOverHKLInfo(d_in,d_out,newhkl);
}
void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::E_sigE> &d_in, clipper::HKL_data<clipper::data32::E_sigE> &d_out,  const clipper::HKL_info &newhkl){
    CopyOverHKLInfo(d_in,d_out,newhkl);
}
void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::ABCD> &d_in, clipper::HKL_data<clipper::data32::ABCD> &d_out,  const clipper::HKL_info &newhkl){
    CopyOverHKLInfo(d_in,d_out,newhkl);
}
void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::Phi_fom> &d_in, clipper::HKL_data<clipper::data32::Phi_fom> &d_out,  const clipper::HKL_info &newhkl){
    CopyOverHKLInfo(d_in,d_out,newhkl);
}
void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::Flag> &d_in, clipper::HKL_data<clipper::data32::Flag> &d_out,  const clipper::HKL_info &newhkl){
    CopyOverHKLInfo(d_in,d_out,newhkl);
}
void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::F_phi> &d_in, clipper::HKL_data<clipper::data32::F_phi> &d_out,  const clipper::HKL_info &newhkl){
    CopyOverHKLInfo(d_in,d_out,newhkl);
}

template <typename T> void CopyIfF_sigFRefNotMissing_float(const T &d_in, T &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref){
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  for (HRI ih = d_ref.first(); !ih.last(); ih.next() ) {
    if (!d_ref[ih].missing()) {
      d_out[ih] = d_in[ih.hkl()];
    } else {
      d_out[ih].set_null();
    }
  }
}

void CopyIfF_sigFRefNotMissingF_sigF_float(const clipper::HKL_data<clipper::data32::F_sigF> &d_in, clipper::HKL_data<clipper::data32::F_sigF> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref){
    CopyIfF_sigFRefNotMissing_float(d_in,d_out,d_ref);
}

void CopyIfF_sigFRefNotMissingF_sigF_ano_float(const clipper::HKL_data<clipper::data32::F_sigF_ano> &d_in, clipper::HKL_data<clipper::data32::F_sigF_ano> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref){
    CopyIfF_sigFRefNotMissing_float(d_in,d_out,d_ref);
}

void CopyIfF_sigFRefNotMissingE_sigE_float(const clipper::HKL_data<clipper::data32::E_sigE> &d_in, clipper::HKL_data<clipper::data32::E_sigE> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref){
    CopyIfF_sigFRefNotMissing_float(d_in,d_out,d_ref);
}

void CopyIfF_sigFRefNotMissingABCD_float(const clipper::HKL_data<clipper::data32::ABCD> &d_in, clipper::HKL_data<clipper::data32::ABCD> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref){
    CopyIfF_sigFRefNotMissing_float(d_in,d_out,d_ref);
}

void CopyIfF_sigFRefNotMissingPhi_fom_float(const clipper::HKL_data<clipper::data32::Phi_fom> &d_in, clipper::HKL_data<clipper::data32::Phi_fom> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref){
    CopyIfF_sigFRefNotMissing_float(d_in,d_out,d_ref);
}

void CopyIfF_sigFRefNotMissingF_phi_float(const clipper::HKL_data<clipper::data32::F_phi> &d_in, clipper::HKL_data<clipper::data32::F_phi> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref){
    CopyIfF_sigFRefNotMissing_float(d_in,d_out,d_ref);
}

void CopyIfF_sigFRefNotMissingFlag_float(const clipper::HKL_data<clipper::data32::Flag> &d_in, clipper::HKL_data<clipper::data32::Flag> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref){
    CopyIfF_sigFRefNotMissing_float(d_in,d_out,d_ref);
}

void PopulateMatchesF_sigF_float(const clipper::HKL_data<clipper::data32::F_sigF> &d_ref, const clipper::HKL_data<clipper::data32::F_sigF> &d, std::vector<clipper::HKL> &matched){
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  for ( HRI ih = d_ref.first(); !ih.last(); ih.next() ) {
    if(d_ref[ih].f()  == d[ih.hkl()].f()  )
       matched.push_back(ih.hkl());
    else if (  d_ref[ih].missing()  && d[ih.hkl()].missing()  )
         matched.push_back(ih.hkl());
    else
      std::cout << ih.hkl().format() << " no match " << clipper::String(  d_ref[ih].f() ) << " " <<  clipper::String( d[ih.hkl()].f() ) << "\n";
  }
}

void PopulateMatchesE_sigE_float(const clipper::HKL_data<clipper::data32::E_sigE> &d_ref, const clipper::HKL_data<clipper::data32::E_sigE> &d, std::vector<clipper::HKL> &matched){
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  for ( HRI ih = d_ref.first(); !ih.last(); ih.next() ) {
    if(d_ref[ih].E()  == d[ih.hkl()].E()  )
       matched.push_back(ih.hkl());
    else if (  d_ref[ih].missing()  && d[ih.hkl()].missing()  )
         matched.push_back(ih.hkl());
    else
      std::cout << ih.hkl().format() << " no match " << clipper::String(  d_ref[ih].E() ) << " " <<  clipper::String( d[ih.hkl()].E() ) << "\n";
  }
}

void PopulateMatchesABCD_float(const clipper::HKL_data<clipper::data32::ABCD> &d_ref, const clipper::HKL_data<clipper::data32::ABCD> &d, std::vector<clipper::HKL> &matched){
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  for ( HRI ih = d_ref.first(); !ih.last(); ih.next() ) {
    if(d_ref[ih].a()  == d[ih.hkl()].a()&&d_ref[ih].b()  == d[ih.hkl()].b()&&d_ref[ih].c()  == d[ih.hkl()].c()&&d_ref[ih].d()  == d[ih.hkl()].d())
       matched.push_back(ih.hkl());
    else if (  d_ref[ih].missing()  && d[ih.hkl()].missing()  )
         matched.push_back(ih.hkl());
    else
      std::cout << ih.hkl().format() << " no match " << clipper::String(  d_ref[ih].a() ) << " " <<  clipper::String( d[ih.hkl()].a() ) << "\n";
  }
}

void PopulateMatchesPhi_fom_float(const clipper::HKL_data<clipper::data32::Phi_fom> &d_ref, const clipper::HKL_data<clipper::data32::Phi_fom> &d, std::vector<clipper::HKL> &matched){
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  for ( HRI ih = d_ref.first(); !ih.last(); ih.next() ) {
    if(d_ref[ih].phi()  == d[ih.hkl()].phi()  )
       matched.push_back(ih.hkl());
    else if (  d_ref[ih].missing()  && d[ih.hkl()].missing()  )
         matched.push_back(ih.hkl());
    else
      std::cout << ih.hkl().format() << " no match " << clipper::String(  d_ref[ih].phi() ) << " " <<  clipper::String( d[ih.hkl()].phi() ) << "\n";
  }
}

void PopulateMatchesF_phi_float(const clipper::HKL_data<clipper::data32::F_phi> &d_ref, const clipper::HKL_data<clipper::data32::F_phi> &d, std::vector<clipper::HKL> &matched){
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  for ( HRI ih = d_ref.first(); !ih.last(); ih.next() ) {
    if(d_ref[ih].f()  == d[ih.hkl()].f()  )
       matched.push_back(ih.hkl());
    else if (  d_ref[ih].missing()  && d[ih.hkl()].missing()  )
         matched.push_back(ih.hkl());
    else
      std::cout << ih.hkl().format() << " no match " << clipper::String(  d_ref[ih].f() ) << " " <<  clipper::String( d[ih.hkl()].f() ) << "\n";
  }
}

void PopulateMatchesFlag_float(const clipper::HKL_data<clipper::data32::Flag> &d_ref, const clipper::HKL_data<clipper::data32::Flag> &d, std::vector<clipper::HKL> &matched){
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  for ( HRI ih = d_ref.first(); !ih.last(); ih.next() ) {
    if(d_ref[ih].flag()  == d[ih.hkl()].flag()  )
       matched.push_back(ih.hkl());
    else if (  d_ref[ih].missing()  && d[ih.hkl()].missing()  )
         matched.push_back(ih.hkl());
    else
      std::cout << ih.hkl().format() << " no match " << clipper::String(  d_ref[ih].flag() ) << " " <<  clipper::String( d[ih.hkl()].flag() ) << "\n";
  }
}

void SetFlagBoth(clipper::HKL_data<clipper::data32::Flag> &flag){
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  for ( HRI ih = flag.first(); !ih.last(); ih.next() )
      flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
}
void SetFlagBothIfMissing(clipper::HKL_data<clipper::data32::Flag> &flag, const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &myfsigf, const clipper::HKL_data< clipper::datatypes::Flag > &status, int freeflag){
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  for ( HRI ih = flag.first(); !ih.last(); ih.next() )
    if ( !myfsigf[ih].missing() && (status[ih].missing()||status[ih].flag()==freeflag) )
      flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
    else
      flag[ih].flag() = clipper::SFweight_spline<float>::NONE;
}
void SetData(const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &F1, const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &F2, const clipper::String &CHECK, const clipper::String &OPS, const clipper::String &ELSE_OPS){
  typedef clipper::HKL_data_base::HKL_reference_index HRI;
  std::vector<clipper::String> ops = OPS.split(",");
  std::vector<clipper::String> else_ops = ELSE_OPS.split(",");
  std::cout << "SetData" << std::endl;
  if(CHECK=="BOTH_PRESENT"){
    for ( HRI ih = F1.first(); !ih.last(); ih.next() ){
      if( !F1[ih].missing() && !F2[ih].missing() ){
        for(unsigned int i=0;i<ops.size();i++){
          std::vector<clipper::String> vs = ops[i].split("=");
          std::cout << vs[0] << std::endl;
        }
      } else {
        std::cout << "else ..." << std::endl;
        for(unsigned int i=0;i<else_ops.size();i++){
          std::vector<clipper::String> vs = else_ops[i].split("=");
          bool isNeg = false;
          clipper::String theArg=vs[1];
          if(vs[1][0]=='-'){
            isNeg = true;
            theArg = vs[1].split("-")[0];
          }
          double val=0.0;
          bool isNull = false;
          if(theArg=="ZERO"){
            val = 0.0;
          } else if(theArg=="NULL"){
            isNull = true;
          } else if(theArg=="2F") {
            val = F2[ih].f();
          } else if(theArg=="2SIGF") {
            val = F2[ih].sigf();
          } else if(theArg=="2F_PL") {
            val = F2[ih].f_pl();
          } else if(theArg=="2SIGF_PL") {
            val = F2[ih].sigf_pl();
          } else if(theArg=="2F_MI") {
            val = F2[ih].f_mi();
          } else if(theArg=="2SIGF_MI") {
            val = F2[ih].sigf_mi();
          } else if(theArg=="2COV") {
            val = F2[ih].cov();
          } else {
            val = theArg.f();
          }
          if(isNeg){
            val = -val;
          }
          std::cout << vs[0] << " " << theArg << " " << val << std::endl;
        }
      }
    }
  }
}
%}
void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::F_sigF> &d_in, clipper::HKL_data<clipper::data32::F_sigF> &d_out,  const clipper::HKL_info &newhkl);
void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::F_sigF_ano> &d_in, clipper::HKL_data<clipper::data32::F_sigF_ano> &d_out,  const clipper::HKL_info &newhkl);
void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::E_sigE> &d_in, clipper::HKL_data<clipper::data32::E_sigE> &d_out,  const clipper::HKL_info &newhkl);
void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::ABCD> &d_in, clipper::HKL_data<clipper::data32::ABCD> &d_out,  const clipper::HKL_info &newhkl);
void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::Phi_fom> &d_in, clipper::HKL_data<clipper::data32::Phi_fom> &d_out,  const clipper::HKL_info &newhkl);
void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::Flag> &d_in, clipper::HKL_data<clipper::data32::Flag> &d_out,  const clipper::HKL_info &newhkl);
void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::F_phi> &d_in, clipper::HKL_data<clipper::data32::F_phi> &d_out,  const clipper::HKL_info &newhkl);
void CopyIfF_sigFRefNotMissingF_sigF_float(const clipper::HKL_data<clipper::data32::F_sigF> &d_in, clipper::HKL_data<clipper::data32::F_sigF> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref);
void CopyIfF_sigFRefNotMissingF_sigF_ano_float(const clipper::HKL_data<clipper::data32::F_sigF_ano> &d_in, clipper::HKL_data<clipper::data32::F_sigF_ano> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref);
void CopyIfF_sigFRefNotMissingE_sigE_float(const clipper::HKL_data<clipper::data32::E_sigE> &d_in, clipper::HKL_data<clipper::data32::E_sigE> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref);
void CopyIfF_sigFRefNotMissingABCD_float(const clipper::HKL_data<clipper::data32::ABCD> &d_in, clipper::HKL_data<clipper::data32::ABCD> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref);
void CopyIfF_sigFRefNotMissingPhi_fom_float(const clipper::HKL_data<clipper::data32::Phi_fom> &d_in, clipper::HKL_data<clipper::data32::Phi_fom> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref);
void CopyIfF_sigFRefNotMissingF_phi_float(const clipper::HKL_data<clipper::data32::F_phi> &d_in, clipper::HKL_data<clipper::data32::F_phi> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref);
void CopyIfF_sigFRefNotMissingFlag_float(const clipper::HKL_data<clipper::data32::Flag> &d_in, clipper::HKL_data<clipper::data32::Flag> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref);

void PopulateMatchesF_sigF_float(const clipper::HKL_data<clipper::data32::F_sigF> &d_ref, const clipper::HKL_data<clipper::data32::F_sigF> &d, std::vector<clipper::HKL> &matched);
void PopulateMatchesE_sigE_float(const clipper::HKL_data<clipper::data32::E_sigE> &d_ref, const clipper::HKL_data<clipper::data32::E_sigE> &d, std::vector<clipper::HKL> &matched);
void PopulateMatchesABCD_float(const clipper::HKL_data<clipper::data32::ABCD> &d_ref, const clipper::HKL_data<clipper::data32::ABCD> &d, std::vector<clipper::HKL> &matched);
void PopulateMatchesPhi_fom_float(const clipper::HKL_data<clipper::data32::Phi_fom> &d_ref, const clipper::HKL_data<clipper::data32::Phi_fom> &d, std::vector<clipper::HKL> &matched);
void PopulateMatchesF_phi_float(const clipper::HKL_data<clipper::data32::F_phi> &d_ref, const clipper::HKL_data<clipper::data32::F_phi> &d, std::vector<clipper::HKL> &matched);
void PopulateMatchesFlag_float(const clipper::HKL_data<clipper::data32::Flag> &d_ref, const clipper::HKL_data<clipper::data32::Flag> &d, std::vector<clipper::HKL> &matched);
void SetFlagBoth(clipper::HKL_data<clipper::data32::Flag> &flag);
void SetFlagBothIfMissing(clipper::HKL_data<clipper::data32::Flag> &flag, const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &myfsigf, const clipper::HKL_data< clipper::datatypes::Flag > &status, int freeflag);
void SetData(const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &F1, const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &F2, const clipper::String &CHECK, const clipper::String &OPS, const clipper::String &ELSE_OPS);

%include "../clipper/contrib/edcalc.h"
namespace clipper {
  %template(EDcalc_mask_float) EDcalc_mask<float>;
  %template(EDcalc_iso_float) EDcalc_iso<float>;
  %template(EDcalc_aniso_float) EDcalc_aniso<float>;
  %extend EDcalc_mask<float>{
    bool compute( Xmap<float>& xmap, const Atom_list& atoms ) const {
      return (*($self))( xmap, atoms );
    }
    bool compute( NXmap<float>& xmap, const Atom_list& atoms ) const {
      return (*($self))( xmap, atoms );
    }
  }
  %extend EDcalc_iso<float>{
    bool compute( Xmap<float>& xmap, const Atom_list& atoms ) const {
      return (*($self))( xmap, atoms );
    }
    bool compute( NXmap<float>& xmap, const Atom_list& atoms ) const {
      return (*($self))( xmap, atoms );
    }
  }
  %extend EDcalc_aniso<float>{
    bool compute( Xmap<float>& xmap, const Atom_list& atoms ) const {
      return (*($self))( xmap, atoms );
    }
    bool compute( NXmap<float>& xmap, const Atom_list& atoms ) const {
      return (*($self))( xmap, atoms );
    }
  }
}

%include "../clipper/core/resol_fn.h"
%include "../clipper/core/resol_basisfn.h"
%include "../clipper/core/resol_targetfn.h"

namespace clipper {
  %template(TargetFn_scaleEsq_E_sigE_T) TargetFn_scaleEsq<clipper::data32::E_sigE>;
  %template(TargetFn_meanFnth_F_phi_T) TargetFn_meanFnth<clipper::data32::F_phi>;
  %template(TargetFn_scaleF1F2_F_sigF_2_T) TargetFn_scaleF1F2<clipper::data32::F_sigF,clipper::data32::F_sigF>;

  %template (TargetFn_scaleLogF1F2_F_sigF_2_T) TargetFn_scaleLogF1F2<data32::F_sigF,data32::F_sigF>;
  %template (TargetFn_scaleI1I2_I_sigI_2_T) TargetFn_scaleI1I2<data32::I_sigI,data32::I_sigI>;
  %template (TargetFn_scaleLogI1I2_I_sigI_2_T) TargetFn_scaleLogI1I2<data32::I_sigI,data32::I_sigI>;
  %template (TargetFn_meanEnth_E_sigE_T) TargetFn_meanEnth<data32::E_sigE>;
  %template (TargetFn_sigmaa_omegaa_E_sigE_T) TargetFn_sigmaa_omegaa<data32::E_sigE>;
  %template (TargetFn_sigmaa_E_sigE_T) TargetFn_sigmaa<data32::E_sigE>;
}

%{
namespace clipper {
  TargetFn_scaleEsq<clipper::data32::E_sigE> TargetFn_scaleEsq_E_sigE(const clipper::HKL_data<clipper::data32::E_sigE>& hkl_data_) {
     TargetFn_scaleEsq<clipper::data32::E_sigE> a(hkl_data_);
     return a;
  }
  TargetFn_meanFnth<clipper::data32::F_phi> TargetFn_meanFnth_F_phi(const clipper::HKL_data<clipper::data32::F_phi>& hkl_data_, float val) {
     TargetFn_meanFnth<clipper::data32::F_phi> a(hkl_data_, val);
     return a;
  }
  TargetFn_scaleF1F2<clipper::data32::F_sigF,clipper::data32::F_sigF> TargetFn_scaleF1F2_F_sigF_2(const clipper::HKL_data<clipper::data32::F_sigF> &F1,const clipper::HKL_data<clipper::data32::F_sigF> &F2) {
     TargetFn_scaleF1F2<clipper::data32::F_sigF,clipper::data32::F_sigF> a( F1, F2 );
     return a;
  }

  TargetFn_scaleLogF1F2<data32::F_sigF,data32::F_sigF> TargetFn_scaleLogF1F2_F_sigF_2( const HKL_data<data32::F_sigF>& hkl_data1_, const HKL_data<data32::F_sigF>& hkl_data2_ ){
     TargetFn_scaleLogF1F2<data32::F_sigF,data32::F_sigF> a(hkl_data1_,hkl_data2_);
     return a;
  }

  TargetFn_scaleI1I2<data32::I_sigI,data32::I_sigI> TargetFn_scaleI1I2_I_sigI_2( const HKL_data<data32::I_sigI>& hkl_data1_, const HKL_data<data32::I_sigI>& hkl_data2_ ){
     TargetFn_scaleI1I2<data32::I_sigI,data32::I_sigI> a(hkl_data1_,hkl_data2_);
     return a;
  }

  TargetFn_scaleLogI1I2<data32::I_sigI,data32::I_sigI> TargetFn_scaleLogI1I2_I_sigI_2( const HKL_data<data32::I_sigI>& hkl_data1_, const HKL_data<data32::I_sigI>& hkl_data2_ ){
     TargetFn_scaleLogI1I2<data32::I_sigI,data32::I_sigI> a(hkl_data1_,hkl_data2_);
     return a;
  }

  TargetFn_meanEnth<data32::E_sigE> TargetFn_meanEnth_E_sigE( const HKL_data<data32::E_sigE>& hkl_data_, const ftype& n ){
     TargetFn_meanEnth<data32::E_sigE> a(hkl_data_,n);
     return a;
  }

  TargetFn_sigmaa_omegaa<data32::E_sigE> TargetFn_sigmaa_omegaa_E_sigE_2( const HKL_data<data32::E_sigE>& eo, const HKL_data<data32::E_sigE>& ec ){
     TargetFn_sigmaa_omegaa<data32::E_sigE> a(eo,ec);
     return a;
  }

  TargetFn_sigmaa<data32::E_sigE> TargetFn_sigmaa_E_sigE_2( const HKL_data<data32::E_sigE>& eo, const HKL_data<data32::E_sigE>& ec ){
     TargetFn_sigmaa<data32::E_sigE> a(eo,ec);
     return a;
  }

}
%}

namespace clipper {
  TargetFn_scaleEsq<clipper::data32::E_sigE> TargetFn_scaleEsq_E_sigE(const clipper::HKL_data<clipper::data32::E_sigE>& hkl_data_);
  TargetFn_meanFnth<clipper::data32::F_phi> TargetFn_meanFnth_F_phi(const clipper::HKL_data<clipper::data32::F_phi>& hkl_data_, float val);
  TargetFn_scaleF1F2<clipper::data32::F_sigF,clipper::data32::F_sigF> TargetFn_scaleF1F2_F_sigF_2(const clipper::HKL_data<clipper::data32::F_sigF> &F1,const clipper::HKL_data<clipper::data32::F_sigF> &F2);

  TargetFn_scaleLogF1F2<data32::F_sigF,data32::F_sigF> TargetFn_scaleLogF1F2_F_sigF_2( const HKL_data<data32::F_sigF>& hkl_data1_, const HKL_data<data32::F_sigF>& hkl_data2_ );
  TargetFn_scaleI1I2<data32::I_sigI,data32::I_sigI> TargetFn_scaleI1I2_I_sigI_2( const HKL_data<data32::I_sigI>& hkl_data1_, const HKL_data<data32::I_sigI>& hkl_data2_ );
  TargetFn_scaleLogI1I2<data32::I_sigI,data32::I_sigI> TargetFn_scaleLogI1I2_I_sigI_2( const HKL_data<data32::I_sigI>& hkl_data1_, const HKL_data<data32::I_sigI>& hkl_data2_ );
  TargetFn_meanEnth<data32::E_sigE> TargetFn_meanEnth_E_sigE( const HKL_data<data32::E_sigE>& hkl_data_, const ftype& n );
  TargetFn_sigmaa_omegaa<data32::E_sigE> TargetFn_sigmaa_omegaa_E_sigE_2( const HKL_data<data32::E_sigE>& eo, const HKL_data<data32::E_sigE>& ec );
  TargetFn_sigmaa<data32::E_sigE> TargetFn_sigmaa_E_sigE_2( const HKL_data<data32::E_sigE>& eo, const HKL_data<data32::E_sigE>& ec );
}

%include "../clipper/core/nxmap_operator.h"
namespace clipper {
  %template(NXmap_operator_float) NXmap_operator<float>;
}

%include "../clipper/contrib/convolution_search.h"
namespace clipper {
  %template(Convolution_search_slow_float_T) Convolution_search_slow<float>;
  %template(Convolution_search_fft_float_T) Convolution_search_fft<float>;
}

%{
namespace clipper {
  Convolution_search_slow<float> Convolution_search_slow_float(const Xmap<float>& xmap) {
     Convolution_search_slow<float> a(xmap);
     return a;
  }
  Convolution_search_slow<float> Convolution_search_slow_float( Xmap<float>& result, const NXmap<float>& srchval, const Xmap<float>& xmap, const NX_operator& nxop ){
     Convolution_search_slow<float> a(result, srchval, xmap, nxop);
     return a;
  }
  Convolution_search_fft<float> Convolution_search_fft_float(const Xmap<float>& xmap) {
     Convolution_search_fft<float> a(xmap);
     return a;
  }
  Convolution_search_fft<float> Convolution_search_fft_float( Xmap<float>& result, const NXmap<float>& srchval, const Xmap<float>& xmap, const NX_operator& nxop ){
     Convolution_search_fft<float> a(result, srchval, xmap, nxop);
     return a;
  }
}
%}

namespace clipper {
  Convolution_search_slow<float> Convolution_search_slow_float(const Xmap<float>& xmap) ;
  Convolution_search_slow<float> Convolution_search_slow_float( Xmap<float>& result, const NXmap<float>& srchval, const Xmap<float>& xmap, const NX_operator& nxop );
  Convolution_search_fft<float> Convolution_search_fft_float(const Xmap<float>& xmap) ;
  Convolution_search_fft<float> Convolution_search_fft_float( Xmap<float>& result, const NXmap<float>& srchval, const Xmap<float>& xmap, const NX_operator& nxop );
}

namespace clipper {
  %extend Convolution_search_slow<float> {
    bool compute( Xmap<float>& res, const NXmap<float>& srchval, const NX_operator& nxop ) const { 
      return (*($self))( res, srchval,  nxop );
    }
  }
  %extend Convolution_search_fft<float> {
    bool compute( Xmap<float>& res, const NXmap<float>& srchval, const NX_operator& nxop ) const { 
      return (*($self))( res, srchval,  nxop );
    }
  }
}

%include "../clipper/contrib/mapfilter.h"
namespace clipper {
  %template(MapFilter_slow_float) MapFilter_slow<float>;
  %template(MapFilter_fft_float) MapFilter_fft<float>;
}

%apply bool *OUTPUT {bool &invert};
%include "../clipper/contrib/originmatch.h"
// FIXME - need a typemap or something for invert return value
namespace clipper {
  %template(OriginMatch_float) OriginMatch<float>;
}

%include "../clipper/core/atomsf.h"
%include "../clipper/core/rotation.h"



