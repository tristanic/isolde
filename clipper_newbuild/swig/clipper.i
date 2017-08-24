%{
#define SWIG_FILE_WITH_INIT
#include <string>
#include <vector>
%}
%include "numpy.i"
%include "pyfragments.i"

%init %{
    import_array();
%}


#define PYTHON_PROPERTIES

%module(directors="1") clipper
%include "std_vector.i"
%include "std_string.i"
//%include "std_complex.i"
%include "exception.i"
%include "std_except.i"
//%include "typemaps.i"

%include "attribute.i"

//%feature("autodoc", "3");
%include "clipper-doc.i"

#pragma SWIG nowarn=312,325,361,362,363,389,401,501,505

%pythoncode %{
  
def safesplat_int(func):
  '''
  C or C++ functions of the form func(int, int, ...) wrapped into Python
  are problematic when (as is quite common) one wants to call them using
  the splat (*) operator - i.e. func(*array). Specifically, this is 
  prone to failure for numpy.int32 and numpy.float32 types due to 
  incorrect typemapping of the scalars. The safesplat_int and 
  safesplat_float decorators are designed to provide a simple workaround
  for these cases. Used as follows:
  @safesplat_int
  def func(arg1, arg2, arg3):
    do stuff
  ... it will try to do the straightforward splat first. If that fails,
  it will assume the argument is a numpy array, and attempt to convert
  that to something that will work. If *that* fails, it will raise a
  TypeError.
  '''
  def func_wrapper(arg_array):
    try:
      return func(*arg_array)
    except:
      try:
        return func(*(arg_array.tolist()))
      except:
        raise NotImplementedError('Input is not a valid array of integers or is the wrong length!')
  return func_wrapper

def safesplat_float(func):
  '''
  C or C++ functions of the form func(int, int, ...) wrapped into Python
  are problematic when (as is quite common) one wants to call them using
  the splat (*) operator - i.e. func(*array). Specifically, this is 
  prone to failure for numpy.int32 and numpy.float32 types due to 
  incorrect typemapping of the scalars. The safesplat_int and 
  safesplat_float decorators are designed to provide a simple workaround
  for these cases. Used as follows:
  @safesplat_float
  def func(arg1, arg2, arg3):
    do stuff
  ... it will try to do the straightforward splat first. If that fails,
  it will assume the argument is a numpy array, and attempt to convert
  that to something that will work. If *that* fails, it will raise a
  TypeError.
  '''
  def func_wrapper(arg_array):
    try:
      return func(*arg_array)
    except:
      try:
        return func(*(arg_array.astype(float)))
      except:
        raise NotImplementedError('Input is not a valid numeric array or is the wrong length!')
  return func_wrapper


%}




// Director for astyle auto-formatting to ignore this section
// *INDENT-OFF*

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
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *numpy_array, int n)}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *numpy_array, int n1, int n2)}
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(double *numpy_array, int n1, int n2, int n3 )};
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(double *numpy_array, int nu, int nv, int nw )};
%apply (double* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3 ) {(double *numpy_array, int nx, int ny, int nz )};

%apply (double* IN_ARRAY1, int DIM1) {(double *numpy_1d_in, int n)};
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double *numpy_2d_in, int n1, int n2)};
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(double *numpy_3d_in, int n1, int n2, int n3)};
%apply (double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(double *numpy_3d_in, int nu, int nv, int nw)};
%apply (int* ARGOUT_ARRAY1, int DIM1) {(int *numpy_int_1d_out, int num)}
%apply (double* ARGOUT_ARRAY1, int DIM1) {(double *numpy_double_1d_out, int num)}
%apply (int ARGOUT_ARRAY1[ANY]) {(int numpy_int_out[ANY])};
%apply (double ARGOUT_ARRAY1[ANY]) {(double numpy_double_out[ANY])};
%apply (int ARGOUT_ARRAY2[ANY][ANY]) {(int numpy_int_out[ANY][ANY])};
%apply (double ARGOUT_ARRAY2[ANY][ANY]) {(double numpy_double_out[ANY][ANY])};
%apply (int ARGOUT_ARRAY3[ANY][ANY][ANY]) {(int numpy_int_out[ANY][ANY][ANY])};
%apply (double ARGOUT_ARRAY3[ANY][ANY][ANY]) {(double numpy_double_out[ANY][ANY][ANY])};

// (u,v,w) offset in fractional coordinates (e.g. for output of a unit cell)
%apply (double IN_ARRAY1[ANY]) {(double frac_offset[3])};

// For making Vec3 objects and their derivatives from Numpy arrays
%apply (float IN_ARRAY1[ANY]) {(float v[3])}
%apply (double IN_ARRAY1[ANY]) {(double v[3])}
%apply (int IN_ARRAY1[ANY]) {(int v[3])}
%apply (long IN_ARRAY1[ANY]) {(long v[3])}



%apply std::string { clipper::String }
%apply std::string& { clipper::String& }


//%apply std::string { String }
//%apply std::string& { String& }


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

// All operators up to here have been mapped to Python functions in the
// relevant objects. Is there any reason to keep these in the API as functions,
// or could they all be switched to ignore directives?

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
  %template(IntVector) vector<int>;
  %template(FloatVector) vector<float>;
  %template(DoubleVector) vector<double>;
  %template(IntIntVector) vector<vector<int> >;
  %template(FloatFloatVector) vector<vector<float> >;
  %template(DoubleDoubleVector) vector<vector<double> >;

  //%template(_string_list) std::vector< std::string >;
}


// *INDENT-ON*

%{
#include "../clipper/core/clipper_message.h"
%}

%exception {
  try
  {
    $action
  } catch (clipper::Message_fatal m)
  {
    SWIG_exception(SWIG_RuntimeError, m.text().c_str() );
    SWIG_fail;
  } catch (std::out_of_range e)
  {
    const char *errString;
    if ( !strcmp(e.what(), "" ) ) {
      errString = "Index out of range!";
    } else {
      errString = e.what();
    }
    SWIG_exception(SWIG_IndexError, errString );
    SWIG_fail;
  } catch (std::length_error e)
  {
    SWIG_exception(SWIG_ValueError, e.what() );
    SWIG_fail;
  } catch (std::invalid_argument e)
  {
    SWIG_exception(SWIG_ValueError, e.what() );
    SWIG_fail;
  } catch (std::exception e)
  {
    SWIG_exception(SWIG_UnknownError, e.what() );
    SWIG_fail;
  } catch (...)
  {
    SWIG_exception(SWIG_UnknownError, "Unknown error" );
    SWIG_fail;
  }

}

/* I'm thinking that Map_reference_coord is probably best left as an
 * internal C++ object hidden from Python.
 *  -- Tristan
 */

/*
 * Some jiggery-pokery to allow SWIG to handle these important nested
 * classes (needed to find the symop of a given coordinate, amongst other
 * things). SWIG won't wrap them directly, so we have to convince it that
 * they're global classes by copying their definitions wholesale from the
 * nxmap.h and xmap.h headers (here, *before* the header files are loaded),
 * and then apply typemaps (later, after the headers are loaded by both the
 * compiler and by SWIG) mapping NXmap_base::Map_reference_... to
 * NXmap_reference_... and Xmap_base::Map_reference_... to Xmap_reference_...
 *
 * Note: since their constructors require the map to be passed by reference,
 * these objects can only be created in functions defined *after* instantiation
 * of the templates - that is, the functions will have to be defined
 * separately for Xmap<double>, Xmap<float> and Xmap<int>. Que sera, sera.
 *
 *
 */
//namespace clipper
//{
//class NXmap_reference_base
//{
  //public:
    //inline const NXmap_base& base_nxmap() const
    //{
      //return *map_;
    //}
    //inline const int& index() const
    //{
      //return index_;
    //}
    //inline bool last() const
    //{
      //return ( index_ >= map_->grid_.size() );
    //}
  //protected:
    //const NXmap_base* map_;
    //int index_;
//};


//class NXmap_reference_index : public NXmap_reference_base
//{
  //public:
    //NXmap_reference_index() {}
    //explicit NXmap_reference_index( const NXmap_base& map )
    //{
      //map_ = &map;
      //index_ = 0;
    //}
    //NXmap_reference_index( const NXmap_base& map, const Coord_grid& pos )
    //{
      //map_ = &map;
      //index_ = map_->grid_.index( pos );
    //}

    //inline Coord_grid coord() const
    //{
      //return map_->grid_.deindex(index_);
    //}
    //inline const Coord_orth coord_orth() const
    //{
      //return map_->coord_orth( coord().coord_map() );
    //}
    //inline NXmap_reference_index& set_coord( const Coord_grid& pos )
    //{
      //index_ = map_->grid_.index( pos );
      //return *this;
    //}
    //inline NXmap_reference_index& next()
    //{
      //index_++;
      //return *this;
    //}
    //inline int index_offset(const int& du,const int& dv,const int& dw) const
    //{
      //return index_ + du*map_->du + dv*map_->dv + dw*map_->dw;
    //}
//};

//class NXmap_reference_coord : public NXmap_reference_base
//{
  //public:
    //NXmap_reference_coord() {}
    //explicit NXmap_reference_coord( const NXmap_base& map )
    //{
      //map_ = &map;
    //}
    //NXmap_reference_coord( const NXmap_base& map, const Coord_grid& pos )
    //{
      //map_ = &map;
      //set_coord( pos );
    //}
    //inline Coord_grid coord() const
    //{
      //return pos_;
    //}
    //inline const Coord_orth coord_orth() const
    //{
      //return map_->coord_orth( coord().coord_map() );
    //}
    //inline NXmap_reference_coord& set_coord( const Coord_grid& pos )
    //{
      //pos_ = pos;
      //index_ = map_->grid_.index( pos_ );
      //return *this;
    //}
    //inline NXmap_reference_coord& next()
    //{
      //index_++;
      //pos_ = map_->grid_.deindex(index_);
      //return *this;
    //}
    //inline NXmap_reference_coord& next_u()
    //{
      //pos_.u()++;
      //index_ += map_->du;
      //return *this;
    //}
    //inline NXmap_reference_coord& next_v()
    //{
      //pos_.v()++;
      //index_ += map_->dv;
      //return *this;
    //}
    //inline NXmap_reference_coord& next_w()
    //{
      //pos_.w()++;
      //index_ += map_->dw;
      //return *this;
    //}
    //inline NXmap_reference_coord& prev_u()
    //{
      //pos_.u()--;
      //index_ -= map_->du;
      //return *this;
    //}
    //inline NXmap_reference_coord& prev_v()
    //{
      //pos_.v()--;
      //index_ -= map_->dv;
      //return *this;
    //}
    //inline NXmap_reference_coord& prev_w()
    //{
      //pos_.w()--;
      //index_ -= map_->dw;
      //return *this;
    //}
    //inline NXmap_reference_coord& operator =( const Coord_grid& pos )
    //{
      //return set_coord( pos );
    //}
  //protected:
    //Coord_grid pos_;
//};

//class Xmap_reference_base
//{
  //public:
    ////! return the parent Xmap
    //inline const Xmap_base& base_xmap() const
    //{
      //return *map_;
    //}
    ////! Get the index into the map data array
    //inline const int& index() const
    //{
      //return index_;
    //}
    ////! Check for end of map
    //bool last() const
    //{
      //return ( index_ >= map_->map_grid.size() );
    //}
  //protected:
    ////! pointer to map for which this Xmap_reference_index is defined
    //const Xmap_base* map_;
    ////! integer index_ into map data array
    //int index_;
//};

////! Map reference with index-like behaviour
///*! This is a reference to a map coordinate. It behaves like a
  //simple index into the map, but can be easily converted into a
  //coordinate as and when required. It also implements methods for
  //iterating through the unique portion of a map. It is very
  //compact, but coord() involves some overhead and loses any
  //information concerning symmetry equivelents.

  //\note The following methods are inherited from
  //Xmap_reference_base but are documented here for convenience:
  //base_xmap(), index(), last().
//*/
//class Xmap_reference_index : public Xmap_reference_base
//{
  //public:
    ////! Null constructor
    //Xmap_reference_index() {}
    ////! Constructor: takes parent map
    //explicit Xmap_reference_index( const Xmap_base& map )
    //{
      //map_ = &map;
      //index_=0;
      //next();
    //}
    ////! Constructor: takes parent map and coord
    //Xmap_reference_index( const Xmap_base& map, const Coord_grid& pos )
    //{
      //map_ = &map;
      //int sym;
      //map_->find_sym( pos, index_, sym );
    //}
    ////! Get current grid coordinate
    //inline Coord_grid coord() const
    //{
      //return map_->map_grid.deindex(index_);
    //}
    ////! Get current value of orthogonal coordinate
    //inline const Coord_orth coord_orth() const
    //{
      //return Coord_orth( map_->rt_grid_orth.rot() * coord().coord_map() );
    //}
    ////! Set current value of coordinate - optimised for nearby coords
    //inline Xmap_reference_index& set_coord( const Coord_grid& pos )
    //{
      //int sym;
      //map_->find_sym( pos, index_, sym );
      //return *this;
    //}
    ////! Simple increment
    //inline Xmap_reference_index& next()
    //{
      //do {
        //index_++;
        //if ( last() ) break;
      //} while ( map_->asu[index_] != 0 );
      //return *this;
    //}
    ////! Index of neighbouring point
    ///* Use for e.g. peak search. Valid for -1 <= du/dv/dw <= 1 only.
    //\param du/dv/dw Coordinate offset. \return Map index. */
    //inline int index_offset(const int& du,const int& dv,const int& dw) const
    //{
      //int i = index_ + du*map_->du[0] + dv*map_->dv[0] + dw*map_->dw[0];
      //if ( map_->asu[i] != 0 ) {
        //i = map_->map_grid.index( map_->to_map_unit( map_->map_grid.deindex(i).transform( map_->isymop[map_->asu[i]-1] ) ) );
      //}
      //return i;
    //}
    //// inherited functions listed for documentation purposes
    ////-- const Xmap_base& base_xmap() const;
    ////-- const int& index() const;
    ////-- bool last() const;
//};

////! Map reference with coordinate-like behaviour
///*! This is a reference to a map coordinate. It behaves like a
  //coordinate, but also stores the index of the corresponding point
  //in the map, and the symmetry operator required to get there. It
  //also implements methods for iterating through the a map. Since
  //the current coordinate and symmetry are stored, coord() is
  //fast. However, it requires 1 pointer and 5 words of storage.

  //\note The following methods are inherited from
  //Xmap_reference_base but are documented here for convenience:
  //base_xmap(), index(), last().
//*/
//class Xmap_reference_coord : public Xmap_reference_base
//{
  //public:
    ////! Null constructor
    //Xmap_reference_coord() {}
    ////! Constructor: takes parent map
    //explicit Xmap_reference_coord( const Xmap_base& map )
    //{
      //map_ = &map;
      //index_ = 0;
      //next();
    //}
    ////! Constructor: takes parent map and coord
    //Xmap_reference_coord( const Xmap_base& map, const Coord_grid& pos )
    //{
      //map_ = &map;
      //pos_ = pos;
      //map_->find_sym( pos_, index_, sym_ );
    //}
    ////! Get current value of coordinate
    //inline const Coord_grid& coord() const
    //{
      //return pos_;
    //}
    ////! Get current value of orthogonal coordinate
    //inline const Coord_orth coord_orth() const
    //{
      //return Coord_orth( map_->rt_grid_orth.rot() * coord().coord_map() );
    //}
    ////! Get current symmetry operator
    //inline const int& sym() const
    //{
      //return sym_;
    //}
    ////! Set current value of coordinate - optimised for nearby coords
    //Xmap_reference_coord& set_coord( const Coord_grid& pos );
    ////! Simple increment
    ///*! Use of this function resets the stored coordinate and sym */
    //inline Xmap_reference_coord& next()
    //{
      //sym_ = 0;
      //do {
        //index_++;
        //if ( last() ) break;
      //} while ( map_->asu[index_] != 0 );
      //pos_ = map_->map_grid.deindex(index_);
      //return *this;
    //}
    //// Increment u,v,w
    //inline Xmap_reference_coord& next_u()
    //{
      //pos_.u()++;  //!< increment u
      //index_ += map_->du[sym_];
      //if (map_->asu[index_] != 0) edge();
      //return *this;
    //}
    //inline Xmap_reference_coord& next_v()
    //{
      //pos_.v()++;  //!< increment v
      //index_ += map_->dv[sym_];
      //if (map_->asu[index_] != 0) edge();
      //return *this;
    //}
    //inline Xmap_reference_coord& next_w()
    //{
      //pos_.w()++;  //!< increment w
      //index_ += map_->dw[sym_];
      //if (map_->asu[index_] != 0) edge();
      //return *this;
    //}
    //inline Xmap_reference_coord& prev_u()
    //{
      //pos_.u()--;  //!< increment u
      //index_ -= map_->du[sym_];
      //if (map_->asu[index_] != 0) edge();
      //return *this;
    //}
    //inline Xmap_reference_coord& prev_v()
    //{
      //pos_.v()--;  //!< decrement v
      //index_ -= map_->dv[sym_];
      //if (map_->asu[index_] != 0) edge();
      //return *this;
    //}
    //inline Xmap_reference_coord& prev_w()
    //{
      //pos_.w()--;  //!< decrement w
      //index_ -= map_->dw[sym_];
      //if (map_->asu[index_] != 0) edge();
      //return *this;
    //}
    ////! Assignment operator from a coord
    //inline Xmap_reference_coord& operator =( const Coord_grid& pos )
    //{
      //return set_coord( pos );
    //}
    //// inherited functions listed for documentation purposes
    ////-- const Xmap_base& base_xmap() const;
    ////-- const int& index() const;
    ////-- bool last() const;

  //protected:
    ////! Current symop
    //int sym_;
    ////! Current coord
    //Coord_grid pos_;

    ////! Reset index for a different symop when we hit the map border
    //void edge();
//};



//} // namespace clipper


%{
#include <sstream>
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
#include "../clipper/clipper-phs.h"
#include "../clipper/phs/phs_io.h"


  using namespace clipper;
#include <string.h>
%}

namespace clipper
{
%ignore Vec3<int>;

}

%include "messagestream.i"

%inline %{
  namespace clipper
  {
    typedef ftype64  ftype;
    typedef ftype64  xtype;
    typedef float    ftype32;
    typedef double   ftype64;    
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
//   %template(IntVector) vector<int>;
//   %template(IntIntVector) vector<vector<int> >;
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



//%typemap(in) (const clipper::String&)
//{
  //%#if PY_MAJOR_VERSION >= 3
  //char *my_result;

//if (PyUnicode_Check($input)) {
    //PyObject * temp_bytes = PyUnicode_AsEncodedString($input, "ASCII", "strict"); // Owned reference
    //if (temp_bytes != NULL) {
      //my_result = PyBytes_AS_STRING(temp_bytes); // Borrowed pointer
      //my_result = strdup(my_result);
      //Py_DECREF(temp_bytes);
    //} else {
      //std::cout << "Decoding error" << std::endl;
    //}
  //}

  //clipper::String *s = new clipper::String(my_result);
  //%#else
    //std::string ss = PyString_AsString($input);
  //clipper::String *s = new clipper::String(ss);
  //%#endif
  //$1 = s;
//}

namespace clipper
{
%extend String {
  std::string __str__()
  {
    return self->c_str();
  }
  std::string __repr__()
  {
    return self->c_str();
  }
}
}

/* Not really sure if these are needed in Python? Would probably be better to
 * rely on the built-in Python equivalents.
 */
namespace clipper
{
%rename(is_nan_float) Util::is_nan(const ftype32);
%rename(is_nan_double) Util::is_nan(const ftype64);
%rename(is_nan_float_slow) Util::isnan(ftype32);
%rename(is_nan_double_slow) Util::isnan(ftype64);
%rename(set_null_float) Util::set_null(ftype32);
%rename(set_null_double) Util::set_null(ftype64);
%rename(is_null_float) Util::is_null(ftype32);
%rename(is_null_double) Util::is_null(ftype64);
}

%include "../clipper/core/clipper_util.h"
%include "../clipper/core/clipper_types.h"


%include "util.i"


/* We are getting a load of warnings whenever SWIG is trying
   to wrap template base classes that will never be of any use
   in Python, so let's suppress the warnings for these
*/

/* Here we use immutable to deal with const char * vars */

%immutable hall;
%immutable hm;
%immutable lgname;

%include "../clipper/core/symop.h"


%include "../clipper/core/spacegroup_data.h"
%feature ("flatnested","1");
%include "../clipper/core/spacegroup.h"
%feature ("flatnested","0");

namespace clipper
{
%extend Spgr_descr
{
  Spgr_descr(std::string descr) 
  {
    String cstr(descr);
    return new Spgr_descr(cstr);
  }
}// extend Spgr_descr
} // namespace clipper

//#ifdef PYTHON_PROPERTIES
//namespace clipper
//{
//%extend Spacegroup
  //{
  //%pythoncode %{
    //spacegroup_number = property(spacegroup_number)
    //symbol_hall = property(symbol_hall)

  //%}
  //}
//}
//#endif

%include "cell.i"


/*
%typemap(in) (const clipper::String&)
{
   std::string ss = PyString_AsString($input);
   clipper::String *s = new clipper::String(ss);
   $1 = s;
}
*/

%include "coords.i"
%include "atoms.i"
%include "symop_arrays.i"
%include "unit_cell.i"

%include "symop.i"

%include "clipper_types.i"


%include "hkl.i"


%include "maps.i"

%include "io.i"

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

%rename("_range_ptr") clipper::Map_stats::range() const;
%include "../clipper/core/clipper_stats.h"
/*
namespace clipper {
  %template(Range_float) Range<float>;
  %template(Range_double) Range<double>;
}
*/
//FIXME  - We do not really want to be producing new objects. Wrapping constructor properly would be preferred.
//         But SWIG does not let me.
// Believe it or not, this is actually memory safe. SWIG automatically
// wraps the below in code that hands responsibility for memory management
// over to Python - when the Python object is deleted, so is the C++ object.
// You can check this in Python by, for example:
// map = clipper.Xmap_float()
// map_stats = clipper.Map_stats()
// map_stats.thisown
// > True
//
// Any SWIG object for which thisown is True is handled by the Python
// garbage collector.
// - TIC
//

namespace clipper
{
%extend Map_stats {
  Map_stats(const Xmap<float> &m)
  {
    return new Map_stats(m);
  };
  Map_stats(const Xmap<double> &m)
  {
    return new Map_stats(m);
  };

  void range(double numpy_double_out[2])
  {
    numpy_double_out[0] = self->range().min();
    numpy_double_out[1] = self->range().max();
  }
}; // extend Map_stats
} // namespace clipper

%include "../clipper/core/map_utils.h"



namespace clipper
{
%extend MModel {
  MPolymer& __getitem__(int i)
  {
    int array_len = self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    return (*self)[i];
  }
  void __setitem__(int i, MPolymer& mpol)
  {
    int array_len = self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    (*self)[i]=mpol;
    return;
  }
  size_t __len__()
  {
    return self->size();
  }
} // extend MModel


%extend MPolymer {
  MMonomer& __getitem__(int i)
  {
    int array_len = self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    return (*self)[i];
  }
  void __setitem__(int i, MMonomer& mmon)
  {
    int array_len = self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    (*self)[i]=mmon;
    return;
  }
  size_t __len__()
  {
    return self->size();
  }
} // extend MPolymer


%extend MMonomer {
  MAtom& __getitem__(int i)
  {
    int array_len = self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    return (*self)[i];
  }
  void __setitem__(int i, MAtom& atom)
  {
    int array_len = self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    (*self)[i]=atom;
    return;
  }
  size_t __len__()
  {
    return self->size();
  }
} // extend MMonomer

%extend MAtom {
  std::string __str__( )
  {
    return self->id() + " " + self->coord_orth().format();
  }
} // extend MAtom

} // namespace clipper

%include "../clipper/mmdb/clipper_mmdb.h"
%include "../clipper/minimol/minimol.h"
%include "../clipper/minimol/minimol_io.h"
%include "../clipper/minimol/minimol_seq.h"

%include "hkl_datatypes.i"

%include "../clipper/core/ramachandran.h"


%include "../clipper/cif/cif_data_io.h"
namespace clipper
{
%extend CIFfile
{
  %pythoncode %{
    open_read = log_clipper(open_read)
  %}
} // extend CIFfile
} // namespace clipper

%include "../clipper/contrib/sfcalc_obs.h"

%include "../clipper/core/hkl_compute.h"
namespace clipper
{
%template(SFcalc_obs_bulk_float) SFcalc_obs_bulk<float>;
} // namespace clipper

%include "../clipper/contrib/function_object_bases.h"
%include "../clipper/contrib/sfscale.h"
namespace clipper
{
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
} // namespace clipper
%include "../clipper/contrib/sfweight.h"
namespace clipper
{
%template(SFweight_spline_float) SFweight_spline<float>;
} // namespace clipper

%include "../clipper/contrib/sfcalc.h"
namespace clipper
{
%template(SFcalc_iso_sum_float) SFcalc_iso_sum<float>;
%template(SFcalc_aniso_sum_float) SFcalc_aniso_sum<float>;
%template(SFcalc_iso_fft_float) SFcalc_iso_fft<float>;
%template(SFcalc_aniso_fft_float) SFcalc_aniso_fft<float>;
} // namespace clipper

%{

  template <typename T> void CopyOverHKLInfo(const T &d_in, T &d_out,  const clipper::HKL_info &newhkl)
  {
    HKL_info::HKL_reference_index ih;
    for ( ih = newhkl.first(); !ih.last(); ih.next() ) {
      d_out[ih] = d_in[ih.hkl()];
    }
  }

  void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::F_sigF> &d_in, clipper::HKL_data<clipper::data32::F_sigF> &d_out,  const clipper::HKL_info &newhkl)
  {
    CopyOverHKLInfo(d_in,d_out,newhkl);
  }
  void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::F_sigF_ano> &d_in, clipper::HKL_data<clipper::data32::F_sigF_ano> &d_out,  const clipper::HKL_info &newhkl)
  {
    CopyOverHKLInfo(d_in,d_out,newhkl);
  }
  void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::E_sigE> &d_in, clipper::HKL_data<clipper::data32::E_sigE> &d_out,  const clipper::HKL_info &newhkl)
  {
    CopyOverHKLInfo(d_in,d_out,newhkl);
  }
  void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::ABCD> &d_in, clipper::HKL_data<clipper::data32::ABCD> &d_out,  const clipper::HKL_info &newhkl)
  {
    CopyOverHKLInfo(d_in,d_out,newhkl);
  }
  void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::Phi_fom> &d_in, clipper::HKL_data<clipper::data32::Phi_fom> &d_out,  const clipper::HKL_info &newhkl)
  {
    CopyOverHKLInfo(d_in,d_out,newhkl);
  }
  void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::Flag> &d_in, clipper::HKL_data<clipper::data32::Flag> &d_out,  const clipper::HKL_info &newhkl)
  {
    CopyOverHKLInfo(d_in,d_out,newhkl);
  }
  void CopyOverHKLInfo(const clipper::HKL_data<clipper::data32::F_phi> &d_in, clipper::HKL_data<clipper::data32::F_phi> &d_out,  const clipper::HKL_info &newhkl)
  {
    CopyOverHKLInfo(d_in,d_out,newhkl);
  }

  template <typename T> void CopyIfF_sigFRefNotMissing_float(const T &d_in, T &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref)
  {
    typedef clipper::HKL_data_base::HKL_reference_index HRI;
    for (HRI ih = d_ref.first(); !ih.last(); ih.next() ) {
      if (!d_ref[ih].missing()) {
        d_out[ih] = d_in[ih.hkl()];
      } else {
        d_out[ih].set_null();
      }
    }
  }

  void CopyIfF_sigFRefNotMissingF_sigF_float(const clipper::HKL_data<clipper::data32::F_sigF> &d_in, clipper::HKL_data<clipper::data32::F_sigF> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref)
  {
    CopyIfF_sigFRefNotMissing_float(d_in,d_out,d_ref);
  }

  void CopyIfF_sigFRefNotMissingF_sigF_ano_float(const clipper::HKL_data<clipper::data32::F_sigF_ano> &d_in, clipper::HKL_data<clipper::data32::F_sigF_ano> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref)
  {
    CopyIfF_sigFRefNotMissing_float(d_in,d_out,d_ref);
  }

  void CopyIfF_sigFRefNotMissingE_sigE_float(const clipper::HKL_data<clipper::data32::E_sigE> &d_in, clipper::HKL_data<clipper::data32::E_sigE> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref)
  {
    CopyIfF_sigFRefNotMissing_float(d_in,d_out,d_ref);
  }

  void CopyIfF_sigFRefNotMissingABCD_float(const clipper::HKL_data<clipper::data32::ABCD> &d_in, clipper::HKL_data<clipper::data32::ABCD> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref)
  {
    CopyIfF_sigFRefNotMissing_float(d_in,d_out,d_ref);
  }

  void CopyIfF_sigFRefNotMissingPhi_fom_float(const clipper::HKL_data<clipper::data32::Phi_fom> &d_in, clipper::HKL_data<clipper::data32::Phi_fom> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref)
  {
    CopyIfF_sigFRefNotMissing_float(d_in,d_out,d_ref);
  }

  void CopyIfF_sigFRefNotMissingF_phi_float(const clipper::HKL_data<clipper::data32::F_phi> &d_in, clipper::HKL_data<clipper::data32::F_phi> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref)
  {
    CopyIfF_sigFRefNotMissing_float(d_in,d_out,d_ref);
  }

  void CopyIfF_sigFRefNotMissingFlag_float(const clipper::HKL_data<clipper::data32::Flag> &d_in, clipper::HKL_data<clipper::data32::Flag> &d_out,  const clipper::HKL_data<clipper::data32::F_sigF> &d_ref)
  {
    CopyIfF_sigFRefNotMissing_float(d_in,d_out,d_ref);
  }

  void PopulateMatchesF_sigF_float(const clipper::HKL_data<clipper::data32::F_sigF> &d_ref, const clipper::HKL_data<clipper::data32::F_sigF> &d, std::vector<clipper::HKL> &matched)
  {
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

  void PopulateMatchesE_sigE_float(const clipper::HKL_data<clipper::data32::E_sigE> &d_ref, const clipper::HKL_data<clipper::data32::E_sigE> &d, std::vector<clipper::HKL> &matched)
  {
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

  void PopulateMatchesABCD_float(const clipper::HKL_data<clipper::data32::ABCD> &d_ref, const clipper::HKL_data<clipper::data32::ABCD> &d, std::vector<clipper::HKL> &matched)
  {
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

  void PopulateMatchesPhi_fom_float(const clipper::HKL_data<clipper::data32::Phi_fom> &d_ref, const clipper::HKL_data<clipper::data32::Phi_fom> &d, std::vector<clipper::HKL> &matched)
  {
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

  void PopulateMatchesF_phi_float(const clipper::HKL_data<clipper::data32::F_phi> &d_ref, const clipper::HKL_data<clipper::data32::F_phi> &d, std::vector<clipper::HKL> &matched)
  {
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

  void PopulateMatchesFlag_float(const clipper::HKL_data<clipper::data32::Flag> &d_ref, const clipper::HKL_data<clipper::data32::Flag> &d, std::vector<clipper::HKL> &matched)
  {
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

  void SetFlagBoth(clipper::HKL_data<clipper::data32::Flag> &flag)
  {
    typedef clipper::HKL_data_base::HKL_reference_index HRI;
    for ( HRI ih = flag.first(); !ih.last(); ih.next() )
      flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
  }
  void SetFlagBothIfMissing(clipper::HKL_data<clipper::data32::Flag> &flag, const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &myfsigf, const clipper::HKL_data< clipper::datatypes::Flag > &status, int freeflag)
  {
    typedef clipper::HKL_data_base::HKL_reference_index HRI;
    for ( HRI ih = flag.first(); !ih.last(); ih.next() )
      if ( !myfsigf[ih].missing() && (status[ih].missing()||status[ih].flag()==freeflag) )
        flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
      else
        flag[ih].flag() = clipper::SFweight_spline<float>::NONE;
  }
  void SetData(const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &F1, const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &F2, const clipper::String &CHECK, const clipper::String &OPS, const clipper::String &ELSE_OPS)
  {
    typedef clipper::HKL_data_base::HKL_reference_index HRI;
    std::vector<clipper::String> ops = OPS.split(",");
    std::vector<clipper::String> else_ops = ELSE_OPS.split(",");
    std::cout << "SetData" << std::endl;
    if(CHECK=="BOTH_PRESENT") {
      for ( HRI ih = F1.first(); !ih.last(); ih.next() ) {
        if( !F1[ih].missing() && !F2[ih].missing() ) {
          for(unsigned int i=0; i<ops.size(); i++) {
            std::vector<clipper::String> vs = ops[i].split("=");
            std::cout << vs[0] << std::endl;
          }
        } else {
          std::cout << "else ..." << std::endl;
          for(unsigned int i=0; i<else_ops.size(); i++) {
            std::vector<clipper::String> vs = else_ops[i].split("=");
            bool isNeg = false;
            clipper::String theArg=vs[1];
            if(vs[1][0]=='-') {
              isNeg = true;
              theArg = vs[1].split("-")[0];
            }
            double val=0.0;
            bool isNull = false;
            if(theArg=="ZERO") {
              val = 0.0;
            } else if(theArg=="NULL") {
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
            if(isNeg) {
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
namespace clipper
{
%template(EDcalc_mask_float) EDcalc_mask<float>;
%template(EDcalc_iso_float) EDcalc_iso<float>;
%template(EDcalc_aniso_float) EDcalc_aniso<float>;
%extend EDcalc_mask<float> {
  bool compute( Xmap<float>& xmap, const Atom_list& atoms ) const {
    return (*self)( xmap, atoms );
  }
  bool compute( NXmap<float>& xmap, const Atom_list& atoms ) const {
    return (*self)( xmap, atoms );
  }
}
%extend EDcalc_iso<float> {
  bool compute( Xmap<float>& xmap, const Atom_list& atoms ) const {
    return (*self)( xmap, atoms );
  }
  bool compute( NXmap<float>& xmap, const Atom_list& atoms ) const {
    return (*self)( xmap, atoms );
  }
}
%extend EDcalc_aniso<float> {
  bool compute( Xmap<float>& xmap, const Atom_list& atoms ) const {
    return (*self)( xmap, atoms );
  }
  bool compute( NXmap<float>& xmap, const Atom_list& atoms ) const {
    return (*self)( xmap, atoms );
  }
}
}

%include "../clipper/core/resol_fn.h"
%include "../clipper/core/resol_basisfn.h"
%include "../clipper/core/resol_targetfn.h"

namespace clipper
{
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
  TargetFn_scaleEsq<clipper::data32::E_sigE> TargetFn_scaleEsq_E_sigE(const clipper::HKL_data<clipper::data32::E_sigE>& hkl_data_)
  {
    TargetFn_scaleEsq<clipper::data32::E_sigE> a(hkl_data_);
    return a;
  }
  TargetFn_meanFnth<clipper::data32::F_phi> TargetFn_meanFnth_F_phi(const clipper::HKL_data<clipper::data32::F_phi>& hkl_data_, float val)
  {
    TargetFn_meanFnth<clipper::data32::F_phi> a(hkl_data_, val);
    return a;
  }
  TargetFn_scaleF1F2<clipper::data32::F_sigF,clipper::data32::F_sigF> TargetFn_scaleF1F2_F_sigF_2(const clipper::HKL_data<clipper::data32::F_sigF> &F1,const clipper::HKL_data<clipper::data32::F_sigF> &F2)
  {
    TargetFn_scaleF1F2<clipper::data32::F_sigF,clipper::data32::F_sigF> a( F1, F2 );
    return a;
  }

  TargetFn_scaleLogF1F2<data32::F_sigF,data32::F_sigF> TargetFn_scaleLogF1F2_F_sigF_2( const HKL_data<data32::F_sigF>& hkl_data1_, const HKL_data<data32::F_sigF>& hkl_data2_ )
  {
    TargetFn_scaleLogF1F2<data32::F_sigF,data32::F_sigF> a(hkl_data1_,hkl_data2_);
    return a;
  }

  TargetFn_scaleI1I2<data32::I_sigI,data32::I_sigI> TargetFn_scaleI1I2_I_sigI_2( const HKL_data<data32::I_sigI>& hkl_data1_, const HKL_data<data32::I_sigI>& hkl_data2_ )
  {
    TargetFn_scaleI1I2<data32::I_sigI,data32::I_sigI> a(hkl_data1_,hkl_data2_);
    return a;
  }

  TargetFn_scaleLogI1I2<data32::I_sigI,data32::I_sigI> TargetFn_scaleLogI1I2_I_sigI_2( const HKL_data<data32::I_sigI>& hkl_data1_, const HKL_data<data32::I_sigI>& hkl_data2_ )
  {
    TargetFn_scaleLogI1I2<data32::I_sigI,data32::I_sigI> a(hkl_data1_,hkl_data2_);
    return a;
  }

  TargetFn_meanEnth<data32::E_sigE> TargetFn_meanEnth_E_sigE( const HKL_data<data32::E_sigE>& hkl_data_, const ftype& n )
  {
    TargetFn_meanEnth<data32::E_sigE> a(hkl_data_,n);
    return a;
  }

  TargetFn_sigmaa_omegaa<data32::E_sigE> TargetFn_sigmaa_omegaa_E_sigE_2( const HKL_data<data32::E_sigE>& eo, const HKL_data<data32::E_sigE>& ec )
  {
    TargetFn_sigmaa_omegaa<data32::E_sigE> a(eo,ec);
    return a;
  }

  TargetFn_sigmaa<data32::E_sigE> TargetFn_sigmaa_E_sigE_2( const HKL_data<data32::E_sigE>& eo, const HKL_data<data32::E_sigE>& ec )
  {
    TargetFn_sigmaa<data32::E_sigE> a(eo,ec);
    return a;
  }

  }
%}

namespace clipper
{
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
namespace clipper
{
%template(NXmap_operator_float) NXmap_operator<float>;
}

%include "../clipper/contrib/convolution_search.h"
namespace clipper
{
%template(Convolution_search_slow_float_T) Convolution_search_slow<float>;
%template(Convolution_search_fft_float_T) Convolution_search_fft<float>;
}

%{
  namespace clipper {
  Convolution_search_slow<float> Convolution_search_slow_float(const Xmap<float>& xmap)
  {
    Convolution_search_slow<float> a(xmap);
    return a;
  }
  Convolution_search_slow<float> Convolution_search_slow_float( Xmap<float>& result, const NXmap<float>& srchval, const Xmap<float>& xmap, const NX_operator& nxop )
  {
    Convolution_search_slow<float> a(result, srchval, xmap, nxop);
    return a;
  }
  Convolution_search_fft<float> Convolution_search_fft_float(const Xmap<float>& xmap)
  {
    Convolution_search_fft<float> a(xmap);
    return a;
  }
  Convolution_search_fft<float> Convolution_search_fft_float( Xmap<float>& result, const NXmap<float>& srchval, const Xmap<float>& xmap, const NX_operator& nxop )
  {
    Convolution_search_fft<float> a(result, srchval, xmap, nxop);
    return a;
  }
  }
%}

namespace clipper
{
Convolution_search_slow<float> Convolution_search_slow_float(const Xmap<float>& xmap) ;
Convolution_search_slow<float> Convolution_search_slow_float( Xmap<float>& result, const NXmap<float>& srchval, const Xmap<float>& xmap, const NX_operator& nxop );
Convolution_search_fft<float> Convolution_search_fft_float(const Xmap<float>& xmap) ;
Convolution_search_fft<float> Convolution_search_fft_float( Xmap<float>& result, const NXmap<float>& srchval, const Xmap<float>& xmap, const NX_operator& nxop );
}

namespace clipper
{
%extend Convolution_search_slow<float> {
  bool compute( Xmap<float>& res, const NXmap<float>& srchval, const NX_operator& nxop ) const {
    return (*self)( res, srchval,  nxop );
  }
}
%extend Convolution_search_fft<float> {
  bool compute( Xmap<float>& res, const NXmap<float>& srchval, const NX_operator& nxop ) const {
    return (*self)( res, srchval,  nxop );
  }
}
}

%include "../clipper/contrib/mapfilter.h"
namespace clipper
{
%template(MapFilter_slow_float) MapFilter_slow<float>;
%template(MapFilter_fft_float) MapFilter_fft<float>;
}

// *INDENT-OFF*
%apply bool *OUTPUT {bool &invert};
// *INDENT-ON*
%include "../clipper/contrib/originmatch.h"
// FIXME - need a typemap or something for invert return value
namespace clipper
{
%template(OriginMatch_float) OriginMatch<float>;
}

%include "../clipper/core/atomsf.h"
%include "../clipper/core/rotation.h"
%include "../clipper/clipper-phs.h"
%include "../clipper/phs/phs_io.h"

