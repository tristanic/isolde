%{
  #include <string>
  #include <vector>
%}


%module(directors="1") clipper
%include "std_vector.i"
%include "std_string.i"
//%include "std_complex.i"
%include "exception.i"
%include "std_except.i"
//%include "typemaps.i"

%feature("autodoc", "3");

#pragma SWIG nowarn=312,325,361,362,363,389,401,501,505

%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
  import_array();
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



%apply std::string { clipper::String }
%apply std::string& { clipper::String& }

%apply std::string { String }
%apply std::string& { String& }


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
namespace clipper
{
class NXmap_reference_base
{
  public:
    inline const NXmap_base& base_nxmap() const
    {
      return *map_;
    }
    inline const int& index() const
    {
      return index_;
    }
    inline bool last() const
    {
      return ( index_ >= map_->grid_.size() );
    }
  protected:
    const NXmap_base* map_;
    int index_;
};


class NXmap_reference_index : public NXmap_reference_base
{
  public:
    NXmap_reference_index() {}
    explicit NXmap_reference_index( const NXmap_base& map )
    {
      map_ = &map;
      index_ = 0;
    }
    NXmap_reference_index( const NXmap_base& map, const Coord_grid& pos )
    {
      map_ = &map;
      index_ = map_->grid_.index( pos );
    }

    inline Coord_grid coord() const
    {
      return map_->grid_.deindex(index_);
    }
    inline const Coord_orth coord_orth() const
    {
      return map_->coord_orth( coord().coord_map() );
    }
    inline NXmap_reference_index& set_coord( const Coord_grid& pos )
    {
      index_ = map_->grid_.index( pos );
      return *this;
    }
    inline NXmap_reference_index& next()
    {
      index_++;
      return *this;
    }
    inline int index_offset(const int& du,const int& dv,const int& dw) const
    {
      return index_ + du*map_->du + dv*map_->dv + dw*map_->dw;
    }
};

class NXmap_reference_coord : public NXmap_reference_base
{
  public:
    NXmap_reference_coord() {}
    explicit NXmap_reference_coord( const NXmap_base& map )
    {
      map_ = &map;
    }
    NXmap_reference_coord( const NXmap_base& map, const Coord_grid& pos )
    {
      map_ = &map;
      set_coord( pos );
    }
    inline Coord_grid coord() const
    {
      return pos_;
    }
    inline const Coord_orth coord_orth() const
    {
      return map_->coord_orth( coord().coord_map() );
    }
    inline NXmap_reference_coord& set_coord( const Coord_grid& pos )
    {
      pos_ = pos;
      index_ = map_->grid_.index( pos_ );
      return *this;
    }
    inline NXmap_reference_coord& next()
    {
      index_++;
      pos_ = map_->grid_.deindex(index_);
      return *this;
    }
    inline NXmap_reference_coord& next_u()
    {
      pos_.u()++;
      index_ += map_->du;
      return *this;
    }
    inline NXmap_reference_coord& next_v()
    {
      pos_.v()++;
      index_ += map_->dv;
      return *this;
    }
    inline NXmap_reference_coord& next_w()
    {
      pos_.w()++;
      index_ += map_->dw;
      return *this;
    }
    inline NXmap_reference_coord& prev_u()
    {
      pos_.u()--;
      index_ -= map_->du;
      return *this;
    }
    inline NXmap_reference_coord& prev_v()
    {
      pos_.v()--;
      index_ -= map_->dv;
      return *this;
    }
    inline NXmap_reference_coord& prev_w()
    {
      pos_.w()--;
      index_ -= map_->dw;
      return *this;
    }
    inline NXmap_reference_coord& operator =( const Coord_grid& pos )
    {
      return set_coord( pos );
    }
  protected:
    Coord_grid pos_;
};

class Xmap_reference_base
{
  public:
    //! return the parent Xmap
    inline const Xmap_base& base_xmap() const
    {
      return *map_;
    }
    //! Get the index into the map data array
    inline const int& index() const
    {
      return index_;
    }
    //! Check for end of map
    bool last() const
    {
      return ( index_ >= map_->map_grid.size() );
    }
  protected:
    //! pointer to map for which this Xmap_reference_index is defined
    const Xmap_base* map_;
    //! integer index_ into map data array
    int index_;
};

//! Map reference with index-like behaviour
/*! This is a reference to a map coordinate. It behaves like a
  simple index into the map, but can be easily converted into a
  coordinate as and when required. It also implements methods for
  iterating through the unique portion of a map. It is very
  compact, but coord() involves some overhead and loses any
  information concerning symmetry equivelents.

  \note The following methods are inherited from
  Xmap_reference_base but are documented here for convenience:
  base_xmap(), index(), last().
*/
class Xmap_reference_index : public Xmap_reference_base
{
  public:
    //! Null constructor
    Xmap_reference_index() {}
    //! Constructor: takes parent map
    explicit Xmap_reference_index( const Xmap_base& map )
    {
      map_ = &map;
      index_=0;
      next();
    }
    //! Constructor: takes parent map and coord
    Xmap_reference_index( const Xmap_base& map, const Coord_grid& pos )
    {
      map_ = &map;
      int sym;
      map_->find_sym( pos, index_, sym );
    }
    //! Get current grid coordinate
    inline Coord_grid coord() const
    {
      return map_->map_grid.deindex(index_);
    }
    //! Get current value of orthogonal coordinate
    inline const Coord_orth coord_orth() const
    {
      return Coord_orth( map_->rt_grid_orth.rot() * coord().coord_map() );
    }
    //! Set current value of coordinate - optimised for nearby coords
    inline Xmap_reference_index& set_coord( const Coord_grid& pos )
    {
      int sym;
      map_->find_sym( pos, index_, sym );
      return *this;
    }
    //! Simple increment
    inline Xmap_reference_index& next()
    {
      do {
        index_++;
        if ( last() ) break;
      } while ( map_->asu[index_] != 0 );
      return *this;
    }
    //! Index of neighbouring point
    /* Use for e.g. peak search. Valid for -1 <= du/dv/dw <= 1 only.
    \param du/dv/dw Coordinate offset. \return Map index. */
    inline int index_offset(const int& du,const int& dv,const int& dw) const
    {
      int i = index_ + du*map_->du[0] + dv*map_->dv[0] + dw*map_->dw[0];
      if ( map_->asu[i] != 0 ) {
        i = map_->map_grid.index( map_->to_map_unit( map_->map_grid.deindex(i).transform( map_->isymop[map_->asu[i]-1] ) ) );
      }
      return i;
    }
    // inherited functions listed for documentation purposes
    //-- const Xmap_base& base_xmap() const;
    //-- const int& index() const;
    //-- bool last() const;
};

//! Map reference with coordinate-like behaviour
/*! This is a reference to a map coordinate. It behaves like a
  coordinate, but also stores the index of the corresponding point
  in the map, and the symmetry operator required to get there. It
  also implements methods for iterating through the a map. Since
  the current coordinate and symmetry are stored, coord() is
  fast. However, it requires 1 pointer and 5 words of storage.

  \note The following methods are inherited from
  Xmap_reference_base but are documented here for convenience:
  base_xmap(), index(), last().
*/
class Xmap_reference_coord : public Xmap_reference_base
{
  public:
    //! Null constructor
    Xmap_reference_coord() {}
    //! Constructor: takes parent map
    explicit Xmap_reference_coord( const Xmap_base& map )
    {
      map_ = &map;
      index_ = 0;
      next();
    }
    //! Constructor: takes parent map and coord
    Xmap_reference_coord( const Xmap_base& map, const Coord_grid& pos )
    {
      map_ = &map;
      pos_ = pos;
      map_->find_sym( pos_, index_, sym_ );
    }
    //! Get current value of coordinate
    inline const Coord_grid& coord() const
    {
      return pos_;
    }
    //! Get current value of orthogonal coordinate
    inline const Coord_orth coord_orth() const
    {
      return Coord_orth( map_->rt_grid_orth.rot() * coord().coord_map() );
    }
    //! Get current symmetry operator
    inline const int& sym() const
    {
      return sym_;
    }
    //! Set current value of coordinate - optimised for nearby coords
    Xmap_reference_coord& set_coord( const Coord_grid& pos );
    //! Simple increment
    /*! Use of this function resets the stored coordinate and sym */
    inline Xmap_reference_coord& next()
    {
      sym_ = 0;
      do {
        index_++;
        if ( last() ) break;
      } while ( map_->asu[index_] != 0 );
      pos_ = map_->map_grid.deindex(index_);
      return *this;
    }
    // Increment u,v,w
    inline Xmap_reference_coord& next_u()
    {
      pos_.u()++;  //!< increment u
      index_ += map_->du[sym_];
      if (map_->asu[index_] != 0) edge();
      return *this;
    }
    inline Xmap_reference_coord& next_v()
    {
      pos_.v()++;  //!< increment v
      index_ += map_->dv[sym_];
      if (map_->asu[index_] != 0) edge();
      return *this;
    }
    inline Xmap_reference_coord& next_w()
    {
      pos_.w()++;  //!< increment w
      index_ += map_->dw[sym_];
      if (map_->asu[index_] != 0) edge();
      return *this;
    }
    inline Xmap_reference_coord& prev_u()
    {
      pos_.u()--;  //!< increment u
      index_ -= map_->du[sym_];
      if (map_->asu[index_] != 0) edge();
      return *this;
    }
    inline Xmap_reference_coord& prev_v()
    {
      pos_.v()--;  //!< decrement v
      index_ -= map_->dv[sym_];
      if (map_->asu[index_] != 0) edge();
      return *this;
    }
    inline Xmap_reference_coord& prev_w()
    {
      pos_.w()--;  //!< decrement w
      index_ -= map_->dw[sym_];
      if (map_->asu[index_] != 0) edge();
      return *this;
    }
    //! Assignment operator from a coord
    inline Xmap_reference_coord& operator =( const Coord_grid& pos )
    {
      return set_coord( pos );
    }
    // inherited functions listed for documentation purposes
    //-- const Xmap_base& base_xmap() const;
    //-- const int& index() const;
    //-- bool last() const;

  protected:
    //! Current symop
    int sym_;
    //! Current coord
    Coord_grid pos_;

    //! Reset index for a different symop when we hit the map border
    void edge();
};



}


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


  using namespace clipper;
#include <string.h>
%}

namespace clipper
{
%ignore Vec3<int>;

}


%inline %{
  namespace clipper
  {
  /*! Clipper messages are directed to an ostream (std::cerr by default).
   * This makes them quite difficult to capture and use sensibly in
   * Python. The solution taken here is to redirect to a stringstream,
   * and apply a simple decorator to key functions in Python to check
   * and report the string buffer as necessary.
   *
   * To make use of this, add code similar to the below to your Python project.
   *
    # You should only ever have a single instance of this object in your project.
    clipper_messages = clipper.ClipperMessageStream()

    def log_clipper(func):
       def func_wrapper(*args, **kwargs):
            clipper_messages.clear()
            func(*args, **kwargs)
            message_string = clipper_messages.read_and_clear()
            if message_string:
                pass # Do whatever you want to do with the messages here
        return func_wrapper
   *
   * For any function with the potential to generate a Clipper warning,
   * add the decorator @log_clipper - i.e.
   *
   * @log_clipper
   * def myfunc(args):
   *    do stuff
   *
   *
   */
  class ClipperMessageStream
  {
    private:
      std::ostringstream stream;
    public:
      ClipperMessageStream()
      {
        this->redirect_clipper_messages();
      }
      ~ClipperMessageStream()
      {
        clipper::Message message;
        message.set_stream(std::cerr);
      }

      std::ostringstream& get_stream ()
      {
        return stream;
      }

      std::ostringstream& redirect_clipper_messages()
      {
        clipper::Message message;
        message.set_stream(stream);
        return stream;
      }

      std::string read_and_clear()
      {
        stream.flush();
        std::string ret = stream.str();
        stream.str("");
        return ret;
      }

      void clear()
      {
        stream.str("");
      }

  };
  };
%}


%inline %{
  namespace clipper
  {

  void warn_test()
  {
    Message::message(Message_warn("Test warning"));
  }

  void except_test()
  {
    throw std::length_error("Test exception");
  }


  typedef ftype64  ftype;
  typedef ftype64  xtype;
  typedef float    ftype32;
  typedef double   ftype64;

  
  std::string ClipperStringAsString(const clipper::String &a)
  {
    return (std::string)(a);
  }
  
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

namespace clipper
{
%extend String {
  std::string __str__()
  {
    return ($self)->c_str();
  }
}
}

%include "../clipper/core/clipper_util.h"
%include "../clipper/core/clipper_types.h"


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

}
}



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




%{
  namespace clipper
  {
  /* Some simple functions for sorting and removing duplicates from
   * vectors of (symop index, Coord_frac) or (symop index, Coord_grid) pairs.
   * These are only used internally for finding symmetry operators mapping
   * the atomic model to a region, and don't need to be wrapped.
   */
  bool compare_int_Coord_frac_pairs (const std::pair<int, Coord_frac > &a, const std::pair<int, Coord_frac > &b)
  {
    if (a.first < b.first) return true;
    if (a.first > b.first) return false;
    for (int i = 0; i < 3; i++) {
      if (a.second[i] > b.second[i]) return false;
      if (a.second[i] < b.second[i]) return true;
    }
    return false;
  }

  bool int_Coord_frac_pairs_equal (const std::pair<int, Coord_frac > &a, const std::pair<int, Coord_frac > &b)
  {
    double eps = 1e-6;
    if (a.first != b.first) return false;
    for (int i=0; i < 3; i++) {
      if (std::abs(a.second[i] - b.second[i]) > eps) return false;
    }
    return true;
  }


  bool compare_int_Coord_grid_pairs (const std::pair<int, Coord_grid > &a, const std::pair<int, Coord_grid > &b)
  {
    if (a.first < b.first) return true;
    if (a.first > b.first) return false;
    for (int i = 0; i < 3; i++) {
      if (a.second[i] > b.second[i]) return false;
      if (a.second[i] < b.second[i]) return true;
    }
    return false;
  }

  bool int_Coord_grid_pairs_equal (const std::pair<int, Coord_grid > &a, const std::pair<int, Coord_grid > &b)
  {

    if (a.first != b.first) return false;
    for (int i=0; i < 3; i++) {
      if (a.second[i] != b.second[i]) return false;
    }
    return true;
  }

  }

%}


%include "../clipper/core/spacegroup_data.h"
%feature ("flatnested","1");
%include "../clipper/core/spacegroup.h"
%feature ("flatnested","0");

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
%ignore Cell::descr() const;
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
}
%include "../clipper/core/cell.h"

namespace clipper
{ 


%extend Cell {
  Mat33<float> matrix_orth()
  {
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
  Mat33<float> matrix_frac()
  {
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
  void dim(double numpy_double_out[3])
  {
    numpy_double_out[0] = $self->a();
    numpy_double_out[1] = $self->b();
    numpy_double_out[2] = $self->c();
  }
  void recip_dim(double numpy_double_out[3])
  {
    numpy_double_out[0] = $self->a_star();
    numpy_double_out[1] = $self->b_star();
    numpy_double_out[2] = $self->c_star();
  }
  
  void angles(double numpy_double_out[3])
  {
    numpy_double_out[0] = $self->alpha();
    numpy_double_out[1] = $self->beta();
    numpy_double_out[2] = $self->gamma();
  }
  void angles_deg(double numpy_double_out[3])
  {
    numpy_double_out[0] = $self->alpha_deg();
    numpy_double_out[1] = $self->beta_deg();
    numpy_double_out[2] = $self->gamma_deg();
  }
  void recip_angles(double numpy_double_out[3])
  {
    numpy_double_out[0] = $self->alpha_star();
    numpy_double_out[1] = $self->beta_star();
    numpy_double_out[2] = $self->gamma_star();
  }
  void recip_angles_deg(double numpy_double_out[3])
  {
    numpy_double_out[0] = Util::rad2d($self->alpha_star());
    numpy_double_out[1] = Util::rad2d($self->beta_star());
    numpy_double_out[2] = Util::rad2d($self->gamma_star());
  }
  // Functional equivalent to Cell::descr
  Cell_descr cell_descr()
  {
    return Cell_descr($self->a(), $self->b(), $self->c(),
                      $self->alpha(), $self->beta(), $self->gamma());
  }

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
%{
  namespace clipper
  {
  template <class T = Coord_frac> T cell_shift (const T& coord)
  {
    T ret;
    for (size_t i = 0; i < 3; i++) {
      ret[i] = floor(coord[i]);
    }
    return ret;
  }

  Coord_grid cell_shift (const Coord_grid& coord, const Grid_sampling& grid)
  {
    Coord_grid ret;
    ret = coord - coord.unit(grid);
    return ret;

  }
  }
%}


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
    
  %rename("_element_ptr") Atom::element() const;
}


%include "../clipper/core/coords.h"


namespace clipper
{
%extend Coord_grid {
  /*
   * This is not right. __cmp__() should return a negative integer if self < g2, 
   * a positive integer if self > g2, and zero if self == g2. I don't see
   * how this can be applied to a 3D coordinate.
  bool __cmp__(const Coord_grid& g2)
  {
    return (*($self) == g2);
  }
  */

  bool __ne__(const Coord_grid& g2) { return *self != g2; }
  bool __eq__(const Coord_grid& g2) { return *self == g2; }
  Coord_grid __add__(const Coord_grid& g2) { return *self + g2; }
  Coord_grid __neg__() { return -(*self); }
  Coord_grid __sub__(const Coord_grid& g2) { return *self - g2; }
  Coord_grid __mul__(const int& m) { return m * (*self); }
  Coord_grid __rmul__(const int& m) { return m * (*self); }
  // Transform operator handled in Isymop::__mul__()
  
  void _get_uvw(int numpy_int_out[3])
  {
    for (int i = 0; i < 3; i++) {
      numpy_int_out[i] = ((*($self))[i]);
    }
  }
  void _set_uvw(int v[3])
  {
    $self->u() = v[0];
    $self->v() = v[1];
    $self->w() = v[2];
  }
}

%extend Coord_orth
{
  Coord_orth __add__(const Coord_orth &h2) { return *self + h2; }
  Coord_orth __sub__(const Coord_orth &h2) { return *self - h2; }
  Coord_orth __neg__() { return -(*self); }
  Coord_orth __mul__ ( const double &f ) { return f * (*self); }
  Coord_orth __rmul__ ( const double &f ) { return f * (*self); }
  // Transforms handled in RTop_orth::__mul__()
  void _get_xyz(double numpy_double_out[3])
  {
    for (int i = 0; i < 3; i++) {
      numpy_double_out[i] = (*self)[i];
    }
  }
}

%extend Coord_frac 
{
  Coord_frac __add__(const Coord_frac &h2) { return *self + h2; }
  Coord_frac __sub__(const Coord_frac &h2) { return *self - h2; }
  Coord_frac __neg__() { return -(*self); }
  Coord_frac __mul__ ( const double &f ) { return f * (*self); }
  Coord_frac __rmul__ ( const double &f ) { return f * (*self); }
  // Transforms handled in RTop_frac::__mul__()
  void _get_uvw(double numpy_double_out[3])
  {
    for (int i = 0; i < 3; i++) {
      numpy_double_out[i] = (*self)[i];
    }
  }
}

%extend Coord_map
{
  Coord_map __add__(const Coord_map &h2) { return *self + h2; }
  Coord_map __sub__(const Coord_map &h2) { return *self - h2; }
  Coord_map __neg__() { return -(*self); }
  Coord_map __mul__ ( const double &f ) { return f * (*self); }
  Coord_map __rmul__ ( const double &f ) { return f * (*self); }
  // No transforms for this coordinate type
  void _get_uvw(double numpy_double_out[3])
  {
    for (int i = 0; i < 3; i++) {
      numpy_double_out[i] = (*($self))[i];
    }
  }
  
}

%extend HKL_sampling
{
  bool __eq__(const HKL_sampling& hkl2) { return *self == hkl2; }
  bool __ne__(const HKL_sampling& hkl2) { return !(*self == hkl2); }
}

%extend Atom 
{
  std::string element()
  {
    return $self->element();
  }
}

%extend Atom_list 
{
  Atom& __getitem__(int i)
  {
    int array_len = $self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    return (*($self))[i];
  }
    
  void __setitem__(int i, Atom& atom)
  {
    int array_len = $self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    (*($self))[i] = atom;
    return;
  }
  size_t __len__()
  {
    return ($self)->size();
  }

  Atom pop(int i)
  {
    size_t array_len = $self->size();
    if (i >= array_len || -i > array_len) {
      throw std::out_of_range("");
    }
    Atom ret;
    if (i>=0) ret = (*($self))[i];
    else ret = (*($self))[array_len - i];
    $self -> erase(($self) -> begin() + i);
    return ret;
  }
  
  size_t append (Atom& a)
  {
    ($self)->push_back(a);
    return ($self)->size();
  }
  
  void extend_by(size_t n)
  {
    for (size_t i = 0; i < n; i++ ) {
      ($self)->push_back( Atom() );
    }
  }
  
  void _set_elements(std::vector<std::string> elements)
  {
    size_t in_len = elements.size();
    size_t my_len = ($self) -> size();
    if (in_len != my_len) {
      std::string errString;
      errString = "Input array length of " + std::to_string(in_len)
                  + " does not match target array length of " + std::to_string(my_len);
      throw std::length_error(errString);
    }
    for (size_t i = 0; i < my_len; i++) {
      (*($self))[i].set_element( elements[i] );
    }
  }
  
  std::vector< std::string > _get_elements()
  {
    std::vector < std::string > ret;
    for (size_t i = 0; i < $self -> size(); i++) {
      ret.push_back( (*($self))[i].element() );
    }
    return ret;
  }
  
  //! Quickly set all atomic xyz coordinates from an nx3 numpy array
  void _set_coord_orth(double *numpy_2d_in, int n1, int n2)
  {
    size_t my_len = ($self) -> size();
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
      (*($self))[i].set_coord_orth( Coord_orth( numpy_2d_in[first], numpy_2d_in[first+1], numpy_2d_in[first+2] ) );
    }
    return;
  }
  
  //! Quickly fill an nx3 numpy array with all atomic xyz coordinates
  void _get_coord_orth(double* numpy_array, int n1, int n2)
  {
    size_t my_len = ($self) -> size();
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
      Coord_orth thiscoord = (*($self))[i].coord_orth();
      numpy_array[first] = thiscoord[0];
      numpy_array[first+1] = thiscoord[1];
      numpy_array[first+2] = thiscoord[2];
    }
  }
    
  void _set_occupancies(double *numpy_1d_in, int n)
  {
    size_t my_len = ($self) -> size();
    if (n != my_len) {
      std::string errString = "Input array length of " + std::to_string(n) +
                              " does not match target array length of " + std::to_string(my_len);
      throw std::length_error(errString);
    }
    for (size_t i = 0; i < n; i++) {
      (*($self))[i].set_occupancy( numpy_1d_in[i] );
    }
    return;
  }
  
  void _get_occupancies(double *numpy_array, int n)
  {
    size_t my_len = ($self) -> size();
    if (n != my_len) {
      std::string errString = "Target array length of " + std::to_string(n) +
                              " does not match Atom_list length of " + std::to_string(my_len);
      throw std::length_error(errString);
    }
    for (size_t i = 0; i < n; i++) {
      numpy_array[i] = (*($self))[i].occupancy();
    }
  }
    
  void _set_u_isos(double *numpy_1d_in, int n)
  {
    size_t my_len = ($self) -> size();
    if (n != my_len) {
      std::string errString = "Input array length of " + std::to_string(n) +
                              " does not match target array length of " + std::to_string(my_len);
      throw std::length_error(errString);
    }
    for (int i = 0; i < n; i++) {
      (*($self))[i].set_u_iso( numpy_1d_in[i] );
    }
    return;
  }
  
  void _get_u_isos(double *numpy_array, int n)
  {
    size_t my_len = ($self) -> size();
    if (n != my_len) {
      std::string errString = "Target array length of " + std::to_string(n) +
                              " does not match Atom_list length of " + std::to_string(my_len);
      throw std::length_error(errString);
    }
    for (size_t i = 0; i < n; i++) {
      numpy_array[i] = (*($self))[i].u_iso();
    }
  }
  
  void _set_u_anisos(double *numpy_2d_in, int n1, int n2)
  {
    size_t my_len = ($self) -> size();
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
      (*($self))[i].set_u_aniso_orth(
        U_aniso_orth(ai[first], ai[first+1], ai[first+2],
                     ai[first+3], ai[first+4], ai[first+5]) );
    }
  }

  void _get_u_anisos(double* numpy_array, int n1, int n2)
  {
    size_t my_len = ($self) -> size();
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
      U_aniso_orth thisaniso = (*($self))[i].u_aniso_orth();
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
    for (Atom_list::const_iterator it = (*self).begin(); it != (*self).end(); ++it) {
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
  

}
} // namespace clipper

%ignore clipper::Symops::size() const;
%rename ("_symops_ptr") clipper::Unit_Cell_clipper::Symops::symops();
%rename("_ref_asu_ptr") clipper::Unit_Cell::ref_asu();

// *INDENT-OFF*

%apply (double IN_ARRAY1[ANY]) {(double box_origin_xyz[3])};
%apply (double IN_ARRAY1[ANY]) {(double box_size_xyz[3])};
%apply (double IN_ARRAY1[ANY]) {(double box_res_xyz[3])};
%apply (int IN_ARRAY1[ANY]) {(int box_size_uvw[3])};
%apply (double IN_ARRAY2[ANY][ANY]) {(double model_bounds_xyz[8][3])};
%apply (double IN_ARRAY1[ANY]) {(double sides_xyz[3])};

// *INDENT-ON*


%inline %{

  namespace clipper
  {
  /*! Convenience class to hold a vector of Symop objects, with
   *  Python getter/setter etc. functions.
   */

  class Symops
  {
    public:
      Symops()
      {
        size_ = 0;
      }
      Symops(std::vector<clipper::RTop_frac > ops)
      {
        symops_ = ops;
        size_ = symops_.size();
      }
      Symops(std::vector<clipper::Symop > ops)
      {
        size_ = 0;
        for (std::vector<clipper::Symop >::iterator it = ops.begin();
        it != ops.end(); ++it) {
          symops_.push_back(RTop_frac(it->rot(), it->trn()));
          unit_translations_.push_back(cell_shift<>(Coord_frac(it->trn())));
          size_ += 1;
        }
      }

      ~Symops() {}

      //! get element
      inline const clipper::RTop_frac operator []( const int& i ) const { return symops_.at(i); }

      const size_t& size() const {
        return size_;
      }

      clipper::RTop_frac __getitem__(int i)
      {
        i = (i < 0) ? size_ + i : i;
        if (i >= size_ || i < 0) {
          throw std::out_of_range("");
        }
        return symops_[i];
      }
      void __setitem__(int i, clipper::RTop_frac op)
      {
        i = (i < 0) ? size_ + i : i;
        if (i >= size_ || i < 0) {
          throw std::out_of_range("");
        }
        symops_[i] = op;
        unit_translations_[i] = cell_shift<>(Coord_frac(op.trn()));
      }
      size_t __len__()
      {
        return size_;
      }
      void append(clipper::RTop_frac op)
      {
        symops_.push_back(op);
        unit_translations_.push_back(cell_shift<>(Coord_frac(op.trn())));
        size_ += 1;
      }
      clipper::RTop_frac pop(int i)
      {
        clipper::RTop_frac ret = symops_[i];
        symops_.erase(symops_.begin() + i);
        unit_translations_.erase(unit_translations_.begin() + i);
        size_ -= 1;
        return ret;
      }
      Coord_frac cell_trans(int i)
      {
        return unit_translations_[i];
      }
      
      //! Return all fractional symops as a single array of 4x4 matrices
      void all_matrices44_frac(double* numpy_array, int n1, int n2, int n3)
      {
        if (size_ != n1) {
          std::string errstring = "Target array length of " + std::to_string(n1) +
                   "does not match Symops array length of " + std::to_string(size_);
          throw std::length_error(errstring);
        }
        if (n2 !=4 || n3 != 4) {
          std::string errstring = "Target should be an nx4x4 numpy array!";
          throw std::length_error(errstring);
        }
        size_t count = 0;
        for (size_t i = 0; i < n1; i++) {
          RTop_frac thisop = symops_[i];
          for (size_t j = 0; j < n2; j++) {
            for (size_t k = 0; k < n3; k++, count++) {
              double thisval;
              // Fill in the bottom row
              if (j == 3) {
                if (k < 3) {
                  numpy_array[count] = 0;
                } else {
                  numpy_array[count] = 1;
                }
                continue;
              }
              // fill in the translations
              if (k == 3) {
                thisval = thisop.trn()[j];
              } else {
                // fill in the rotations
                thisval = thisop.rot()(j, k);
              }
              numpy_array[count] = thisval;
            }
          }
        } 
      }

      //! Return all fractional symops as a single array of 3x4 matrices
      void all_matrices34_frac(double* numpy_array, int n1, int n2, int n3)
      {
        if (size_ != n1) {
          std::string errstring = "Target array length of " + std::to_string(n1) +
                   "does not match Symops array length of " + std::to_string(size_);
          throw std::length_error(errstring);
        }
        if (n2 !=3 || n3 != 4) {
          std::string errstring = "Target should be an nx3x4 numpy array!";
          throw std::length_error(errstring);
        }
        size_t count = 0;
        for (size_t i = 0; i < n1; i++) {
          RTop_frac thisop = symops_[i];
          for (size_t j = 0; j < n2; j++) {
            for (size_t k = 0; k < n3; k++, count++) {
              double thisval;
              // fill in the translations
              if (k == 3) {
                thisval = thisop.trn()[j];
              } else {
                // fill in the rotations
                thisval = thisop.rot()(j, k);
              }
              numpy_array[count] = thisval;
            }
          }
        } 
      }

      //! Return all orthographic symops as a single array of 4x4 matrices
      void all_matrices44_orth(const Cell& cell, double* numpy_array, int n1, int n2, int n3)
      {
        if (size_ != n1) {
          std::string errstring = "Target array length of " + std::to_string(n1) +
                   "does not match Symops array length of " + std::to_string(size_);
          throw std::length_error(errstring);
        }
        if (n2 !=4 || n3 != 4) {
          std::string errstring = "Target should be an nx4x4 numpy array!";
          throw std::length_error(errstring);
        }
        size_t count = 0;
        for (size_t i = 0; i < n1; i++) {
          RTop_orth thisop = symops_[i].rtop_orth(cell);
          for (size_t j = 0; j < n2; j++) {
            for (size_t k = 0; k < n3; k++, count++) {
              double thisval;
              // Fill in the bottom row
              if (j == 3) {
                if (k < 3) {
                  numpy_array[count] = 0;
                } else {
                  numpy_array[count] = 1;
                }
                continue;
              }
              // fill in the translations
              if (k == 3) {
                thisval = thisop.trn()[j];
              } else {
                // fill in the rotations
                thisval = thisop.rot()(j, k);
              }
              numpy_array[count] = thisval;
            }
          }
        } 
      }

      //! Return all orthographic symops as a single array of 3x4 matrices
      void all_matrices34_orth(const Cell& cell, double* numpy_array, int n1, int n2, int n3)
      {
        if (size_ != n1) {
          std::string errstring = "Target array length of " + std::to_string(n1) +
                   "does not match Symops array length of " + std::to_string(size_);
          throw std::length_error(errstring);
        }
        if (n2 !=3 || n3 != 4) {
          std::string errstring = "Target should be an nx3x4 numpy array!";
          throw std::length_error(errstring);
        }
        size_t count = 0;
        for (size_t i = 0; i < n1; i++) {
          RTop_orth thisop = symops_[i].rtop_orth(cell);
          for (size_t j = 0; j < n2; j++) {
            for (size_t k = 0; k < n3; k++, count++) {
              double thisval;
              // fill in the translations
              if (k == 3) {
                thisval = thisop.trn()[j];
              } else {
                // fill in the rotations
                thisval = thisop.rot()(j, k);
              }
              numpy_array[count] = thisval;
            }
          }
        } 
      }


      
    private:
      std::vector<clipper::RTop_frac >  symops_;
      std::vector<Coord_frac> unit_translations_;
      size_t size_;
  };


  /*! Similar to Symops, but for integerised symops. Initialise from
   *  a Symops and a Grid_sampling object.
   */
  class Isymops
  {
    public:
      Isymops()
      {
        size_ = 0;
      }
      Isymops(const Symops& ops, const Grid_sampling& grid)
      {
        for (size_t i = 0; i < ops.size(); i++) {
          const RTop_frac& this_op = ops[i];
          /* There is no constructor for an Isymop from an RTop_frac,
           * so we need to convert to a Symop first. This somewhat
           * annoyingly discards the unit cell translations, so we'll
           * have to add them back to construct the final operator.
           */
          Symop this_symop = Symop(this_op);
          Isymop this_isym(this_symop, grid);
          Coord_grid this_trn = Coord_frac(this_op.trn()).coord_grid(grid);
          this_isym.trn() = this_trn;
          symops_.push_back(this_isym);
          unit_translations_.push_back(cell_shift(this_trn, grid));
        }
        size_ = symops_.size();
      }

      ~Isymops() {}

      //! get element
      inline const clipper::Isymop& operator []( const int& i ) const { return symops_.at(i); }

      const size_t& size() const {
        return size_;
      }

      clipper::Isymop __getitem__(int i)
      {
        i = (i < 0) ? size_ + i : i;
        if (i >= size_ || i < 0) {
          throw std::out_of_range("");
        }
        return symops_[i];
      }
      void __setitem__(int i, clipper::Isymop op, const Grid_sampling& grid)
      {
        i = (i < 0) ? size_ + i : i;
        if (i >= size_ || i < 0) {
          throw std::out_of_range("");
        }
        symops_[i] = op;
        unit_translations_[i] = cell_shift(Coord_grid(op.trn()), grid);
      }
      size_t __len__()
      {
        return size_;
      }
      void append(clipper::Isymop op, const Grid_sampling& grid)
      {
        symops_.push_back(op);
        unit_translations_.push_back(cell_shift(Coord_grid(op.trn()), grid));
        size_ += 1;
      }
      clipper::Isymop pop(int i)
      {
        clipper::Isymop ret = symops_[i];
        symops_.erase(symops_.begin() + i);
        unit_translations_.erase(unit_translations_.begin() + i);
        size_ -= 1;
        return ret;
      }
      Coord_grid cell_trans(int i)
      {
        return unit_translations_[i];
      }
    private:
      std::vector<clipper::Isymop >  symops_;
      std::vector<Coord_grid> unit_translations_;
      size_t size_;
  };


  //! Condensed unit cell description
  /*! Contains a reference fractional coordinate (e.g. corresponding
   *  to the centroid of the modelled asu), a list of symops
   *  with the necessary translations added to ensure all results
   *  cluster in the same unit cell, and a list of pre-computed inverse
   *  symops for convenience. Everything is defined relative to the
   *  given reference coordinate rather than Clipper's internal reference
   *  asu. 
   */

  class Unit_Cell
  {
  public:
    Unit_Cell() {};
    ~Unit_Cell() {};

    //! Find the symops (with translations) necessary to pack a unit cell
    /*! The reference coordinate would typically be the centroid of
     *  the atomic model. Additionally, create a Grid_range object
     *  defining the grid that just encloses all atoms in the asymmetric
     *  unit.
     *
     */
    Unit_Cell(Coord_frac ref, const Atom_list& atoms, const Cell& cell, const Spacegroup& sg, const Grid_sampling& grid, int padding=0)
    {
      ref_ = ref;
      grid_ = grid;
      cell_ = cell;
      sg_ = sg;
      /* Store the nearest grid coordinate to the input reference coordinate.
       * Make the equivalent grid, fractional and map reference coords
       * and store for repeated use.
       */
      grid_ref_ = ref.coord_grid(grid);
      cell_origin_ = grid_ref_ - grid_ref_.unit(grid);
      ref_ = grid_ref_.coord_frac(grid);

      /* Find the minimum and maximum grid coordinates of a box encompassing the atoms,
       * pad it by padding grid units either side, and make a Grid_range object from it.
       */
      Coord_grid ref_min = atoms[0].coord_orth().coord_frac(cell).coord_grid(grid);
      Coord_grid ref_max = ref_min;
      //Coord_grid pad(1,1,1);
      Coord_grid pad(padding,padding,padding);
      for (Atom_list::const_iterator it = atoms.begin(); it != atoms.end(); ++it) {
        const Coord_grid& thiscoord = it->coord_orth().coord_frac(cell).coord_grid(grid);
        for (size_t i = 0; i < 3; i++) {
          if (thiscoord[i] < ref_min[i]) ref_min[i] = thiscoord[i];
          else if (thiscoord[i] > ref_max[i]) ref_max[i] = thiscoord[i];
        }
      }
      //std::cerr << "Reference box corners: " << ref_min.format().c_str() << " " << ref_max.format().c_str() << std::endl;
      //std::cerr << ref_min.coord_frac(grid).coord_orth(cell).format().c_str() << ref_max.coord_frac(grid).coord_orth(cell).format().c_str() << std::endl;

      reference_model_bounds_ = Grid_range( ref_min-pad, ref_max+pad );
      //std::cerr << reference_model_bounds_.min().format().c_str() << reference_model_bounds_.max().format().c_str() << std::endl;
      
      // In order to catch all corner cases, we'll need to search the
      // block of 9 unit cells surrounding our reference model
      Coord_grid this_offset;
      for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
          for (int k = -1; k < 2; k++) {
            this_offset = Coord_grid(grid.nu()*i, grid.nv()*j, grid.nw()*k);
            ref_model_cell_origins_.push_back(cell_origin_ + this_offset);
          }
        }    
      }
      /* Find the side of the reference box covering the shortest
       * distance. When searching a box for relevant symops, we'll
       * sample it in steps of half this size.
       */

      size_t nu, nv, nw;
      nu = reference_model_bounds_.nu();
      nv = reference_model_bounds_.nv();
      nw = reference_model_bounds_.nw();
      min_side_ = nv < nu ? (nw < nv ? nw : nv) : nu;

      for (size_t i = 0; i < sg.num_symops(); i++) {
        clipper::Symop this_symop = sg.symop(i);
        clipper::Coord_frac tc = ref_.transform(this_symop);
        // now we want to find the [int, int, int] shift that brings the
        // transformed coordinates back into the same unit cell as the reference
        // coordinate.
        clipper::RTop_frac new_op(this_symop.rot(), this_symop.trn());

        Coord_frac shift = tc - tc.lattice_copy_unit() - cell_origin_.coord_frac(grid);
        new_op.trn() -= shift;

        symops_.append(new_op);
      }
      for (size_t i = 0; i < symops_.size(); i++) {
        inv_symops_.append(symops_[i].inverse());
      }
      isymops_ = Isymops(symops_, grid);
      inv_isymops_ = Isymops(inv_symops_, grid);

    }
    // Symmetry operators mapping the reference coordinate to each
    // of the other asymmetric units in this unit cell
    const Symops& symops() const { return symops_;} //getter
    Symops& symops()
    {
      return symops_; //setter
    }
    // Symmetry operators mapping asymmetric units back to the
    // reference asu.
    const Symops& inv_symops() const {return inv_symops_;} //getter
    Symops& inv_symops()
    {
      return inv_symops_; //setter
    }
    // Reference coordinate (by default, the centroid of the atomic
    // model)
    const Isymops& isymops() const { return isymops_;} //getter
    Isymops& isymops()
    {
      return isymops_; //setter
    }
    const Isymops& inv_isymops() const {return inv_isymops_;} //getter
    Isymops& inv_isymops()
    {
      return inv_isymops_; //setter
    }

    const Coord_frac& ref() const {return ref_;} //getter
    Coord_frac& ref()
    {
      return ref_; //setter
    }
    // Asymmetric unit containing the reference coordinate
    const Grid_range& ref_box() const {return reference_model_bounds_;} //getter
    Grid_range& ref_box()
    {
      return reference_model_bounds_; //setter
    }
    
    // The length of the shortest side of the reference box
    const size_t& ref_box_min_side() const {return min_side_;}
    
    Coord_grid min() {return cell_origin_;}
    Coord_grid max() {return cell_origin_
            + Coord_grid(grid_.nu()-1, grid_.nv()-1, grid_.nw()-1);}



    //! Find all symops relevant to a given grid coordinate
    /*! Actually finds all inverse symops that map the coordinate back
     *  into our reference box, and appends a std::pair object containing
     *  the symop index and the unit cell translations to the input
     *  pairlist.
     *  Args:
     *    pairlist: vector to append the found symops to
     *    thecoord: input grid coordinate
     */ 
    void find_symops_of_coord(std::vector<std::pair<int, Coord_grid> >& pairlist, const Coord_grid& thecoord)
    {
      Coord_grid shift;
      Coord_grid t_coord;
      
      for (std::vector<Coord_grid>::iterator it = ref_model_cell_origins_.begin();
            it != ref_model_cell_origins_.end(); ++it) {

        // First work out how many unit cells we are away from our reference
        shift = thecoord - thecoord.unit(grid_) - *it;
        /* Now try the inverse symops, and work out which (if any) put us
         * within the bounds of our reference box. NOTE: since our box
         * is the limits of the model (not an ASU), it's entirely possible
         * for the coordinate to be in more than one box - or, indeed, no box
         * at all. Therefore, we should err on the side of caution and test
         * all symops for each point.
        */
        for (size_t i = 0; i < inv_symops_.size(); i++) {
          t_coord = (thecoord.unit(grid_) + *it).transform(inv_isymops_[i]);
          if (reference_model_bounds_.in_grid(t_coord))
            pairlist.push_back(std::pair<int, Coord_grid>(i, shift));
        }
      }
    }

    //! Find all symops mapping the reference asu to positions in a box
    /*! The box is defined by its origin in xyz coordinates, and the
     *  number of grid coordinates along each axis. 
     */
    Symops all_symops_in_box(double box_origin_xyz[3], int box_size_uvw[3], bool always_include_identity = false, bool debug = false, int sample_frequency = 2)
    {
      Symops ret;

      /* We'll do all the work of finding the symops in integer grid
       * coordinates for speed.
       */

      std::vector<std::pair<int, Coord_grid> > ops_and_translations;
      if (always_include_identity) {
        std::pair<int, Coord_grid> identity;
        identity.first = 0;
        identity.second = Coord_grid(0,0,0);
        ops_and_translations.push_back(identity);
      }

      Coord_grid ref_grid = ref_.coord_grid(grid_);
      int num_symops = isymops_.size();

      Coord_grid box_origin = Coord_orth(box_origin_xyz[0], box_origin_xyz[1], box_origin_xyz[2]).coord_frac(cell_).coord_grid(grid_);
      Coord_grid box_max = box_origin + Coord_grid(box_size_uvw[0], box_size_uvw[1], box_size_uvw[2]);

      /* We'll sample the box in steps equal to 1/2 the length of the
       * shortest side of the box encompassing the atomic model, making
       * sure we also capture the faces and corners.
       */
       
      size_t step_size = min_side_/sample_frequency;
      step_size = step_size == 0 ? 1 : step_size;
      if (debug) std::cerr << "Sampling step size: " << step_size << std::endl;
      bool u_done = false, v_done = false, w_done = false;
      int u,v,w;
      Coord_grid thiscg;
      u = box_origin[0];
      while (!u_done) {
        v=box_origin[1];
        v_done = false;
        if (u == box_max[0]) u_done = true;
        while (!v_done) {
          w=box_origin[2];
          w_done = false;
          if (v == box_max[1]) v_done = true;
          while (!w_done) {
            if (w == box_max[2]) w_done = true;
            thiscg = Coord_grid(u,v,w);
            this -> find_symops_of_coord(ops_and_translations, thiscg);
            w += step_size;
            if ( w > box_max[2] ) w = box_max[2];
          }
          v += step_size;
          if ( v > box_max[1] ) v = box_max[1];
        }
        u += step_size;
        if ( u > box_max[0] ) u = box_max[0];
      }
      // Sort the two vectors, and remove any duplicates
      std::sort(ops_and_translations.begin(), ops_and_translations.end(), compare_int_Coord_grid_pairs);
      if(debug) {
        int count = 0;
        for (std::vector<std::pair<int, Coord_grid > >::iterator it = ops_and_translations.begin();
           it != ops_and_translations.end(); ++it, count++) {
             std::cerr << count << "Sorted" << count <<": " << it->first << " " << it->second.format() << std::endl;
        }
      }
          
      ops_and_translations.erase( std::unique( ops_and_translations.begin(),
                                  ops_and_translations.end(),
                                  int_Coord_grid_pairs_equal
                                             ), ops_and_translations.end() );
      if(debug) {
        int count = 0;
        for (std::vector<std::pair<int, Coord_grid > >::iterator it = ops_and_translations.begin();
           it != ops_and_translations.end(); ++it, count++) {
             std::cerr << count << "Unique" << count <<": " << it->first << " " << it->second.format() << std::endl;
        }
      }

      for (std::vector<std::pair<int, Coord_grid > >::iterator it = ops_and_translations.begin();
           it != ops_and_translations.end(); ++it) {
        RTop_frac thisop = symops_[it->first];
        thisop.trn() += (it->second).coord_frac(grid_);
        ret.append(thisop);

      }
      return ret;

    }


  private:
    Cell cell_;
    Spacegroup sg_;
    Grid_sampling grid_;
    Coord_frac ref_;  // Our reference coordinate in fractional...
    Coord_grid grid_ref_; // ... and grid coordinates.
    Symops symops_; // The list of RTops mapping our reference coordinate to equivalent positions in the same unit cell
    Symops inv_symops_; // For convenience, the inverse symops are also saved.
    Isymops isymops_; // Integerised symops for fast operations on grid coordinates
    Isymops inv_isymops_; //... and their inverses
    std::vector<int> symop_map_; // Mapping of the indices of symops_ to Clipper's internal symop indices
    //Grid_range ref_asu_;
    Coord_grid cell_origin_; // The origin of the unit cell containing our reference model
    Grid_range reference_model_bounds_; // Smallest grid range encompassing our atomic model
    size_t min_side_; // The length of the smallest side of the grid box in reference_model_bounds
    /* Because life is never easy, the reference model will almost certainly
     * span across multiple unit cells (unlike the nice, neat map reference
     * asu used by Clipper. Therefore, to find symmetry equivalents we 
     * need to test the main unit cell plus all 26 direct neighbours.
     */
    std::vector<Coord_grid> ref_model_cell_origins_;
  };

  }
%}


namespace clipper
{


%extend RTop_orth {
  
  Coord_orth __mul__(const Coord_orth& c) {
    return (*($self)) * c;
  }
    
  void matrix (double numpy_double_out[4][4])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = ($self)->rot()(i,j);
      }
    }
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i][3] = ($self)->trn()[i];
      numpy_double_out[3][i] = 0;
    }
    numpy_double_out[3][3] = 1;
  }
  //! Return the affine transform matrix excluding the final [0,0,0,1] row
  void mat34 (double numpy_double_out[3][4])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = ($self)->rot()(i,j);
      }
    }
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i][3] = ($self)->trn()[i];
    }
  }
  void rotation (double numpy_double_out[3][3])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = ($self)->rot()(i,j);
      }
    }
  }

  void translation (double numpy_double_out[3])
  {
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i] = ($self)->trn()[i];
    }
  }
}

%extend RTop_frac {
  
  // No idea why, but this only works if the Coord_frac name is fully justified
  clipper::Coord_frac __mul__(const clipper::Coord_frac& c) {
    return (*($self)) * c;
  }
  
  //! Return a full 4x4 affine transform matrix
  void matrix (double numpy_double_out[4][4])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = ($self)->rot()(i,j);
      }
    }
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i][3] = ($self)->trn()[i];
      numpy_double_out[3][i] = 0;
    }
    numpy_double_out[3][3] = 1;
  }
  //! Return the affine transform matrix excluding the final [0,0,0,1] row
  void mat34 (double numpy_double_out[3][4])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = ($self)->rot()(i,j);
      }
    }
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i][3] = ($self)->trn()[i];
    }
  }
  
  void rotation (double numpy_double_out[3][3])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = ($self)->rot()(i,j);
      }
    }
  }
  void translation (double numpy_double_out[3])
  {
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i] = ($self)->trn()[i];
    }
  }
  // format() is only defined in the base class, so needs to be
  // re-defined here.
  String format() const {
    return $self->format();
  }

  String format_as_symop() const
  {
    String s, t, xyz="xyz";
    for ( int i = 0; i < 3; i++ ) {
      t = "";
      for ( int j = 0; j < 3; j++ )
        if ( $self->rot()(i,j) != 0.0 ) {
    t += ( $self->rot()(i,j) > 0.0 ) ? "+" : "-";
    if ( Util::intr( fabs( $self->rot()(i,j) ) ) != 1 )
      t += String::rational( fabs( $self->rot()(i,j) ), 24 );
    t += xyz[j];
        }
      if ( $self->trn()[i] != 0.0 )
        t += String::rational( $self->trn()[i], 24, true );
      s += t.substr( ( t[0] == '+' ) ? 1 : 0 );
      if ( i < 2 ) s+= ", ";
    }
    return s;
  }
  

  


}

%extend Symop {
  //! Return a full 4x4 affine transform matrix
  clipper::Coord_frac __mul__(const clipper::Coord_frac& c) {
    return (*($self)) * c;
  }
  
  void matrix (double numpy_double_out[4][4])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = ($self)->rot()(i,j);
      }
    }
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i][3] = ($self)->trn()[i];
      numpy_double_out[3][i] = 0;
    }
    numpy_double_out[3][3] = 1;
  }
  //! Return the affine transform matrix excluding the final [0,0,0,1] row
  void mat34 (double numpy_double_out[3][4])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = ($self)->rot()(i,j);
      }
    }
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i][3] = ($self)->trn()[i];
    }
  }
  void rotation (double numpy_double_out[3][3])
  {
    for (size_t i = 0; i < 3; i++) {
      for (size_t j = 0; j < 3; j++) {
        numpy_double_out[i][j] = ($self)->rot()(i,j);
      }
    }
  }
  void translation (double numpy_double_out[3])
  {
    for (size_t i = 0; i < 3; i++) {
      numpy_double_out[i] = ($self)->trn()[i];
    }
  }


}

%extend Isymop
{
  clipper::Coord_grid __mul__(const clipper::Coord_grid& c) {
    return (*($self)) * c;
  }
  clipper::HKL __mul__(const clipper::HKL& c) {
    return (*($self)) * c;
  }
}

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
}




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
        numpy_double_out[i][j] = (*($self))(i,j);
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
};
}
%template (mat33_float)  clipper::Mat33<float>;
%template (mat33_ftype)  clipper::Mat33<ftype>;
//%template (mat33_double) clipper::Mat33<double>;

%template (RTop_float)  clipper::RTop<float>;
%template (RTop_double) clipper::RTop<double>;
%template (vec3_int)    clipper::Vec3<int>;
%template (vec3_float)  clipper::Vec3<float>;
%template (vec3_double) clipper::Vec3<double>;

namespace clipper
{
%extend mat33_float
{
  clipper::Coord_orth __mul__ (const clipper::Coord_orth &v_)
  {
    Vec3<float> v;
    v[0] = v_[0];
    v[1] = v_[1];
    v[2] = v_[2];
    Vec3<float> v2 = ((*self)*v);
    Coord_orth c(v2[0],v2[1],v2[2]);
    return c;
  };
  clipper::Coord_frac __mul__ (const clipper::Coord_frac &v_)
  {
    Vec3<float> v;
    v[0] = v_[0];
    v[1] = v_[1];
    v[2] = v_[2];
    Vec3<float> v2 = ((*self)*v);
    Coord_frac c(v2[0],v2[1],v2[2]);
    return c;
  };  
}
  
}


%include "../clipper/core/hkl_info.h"
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
};
class HKL_reference_index : public HKL_reference_base
{
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

%{
  namespace clipper {
  typedef NXmap_base::Map_reference_base NXmap_reference_base;
  typedef NXmap_base::Map_reference_index NXmap_reference_index;
  typedef NXmap_base::Map_reference_coord NXmap_reference_coord;
  typedef Xmap_base::Map_reference_base Xmap_reference_base;
  typedef Xmap_base::Map_reference_index Xmap_reference_index;
  typedef Xmap_base::Map_reference_coord Xmap_reference_coord;
  }

%}

// numpy support for Xmap and NXmap
// jon is currently editing this code piece


namespace clipper
{
%extend Grid_sampling {
  void dim(int numpy_int_out[3])
  {
    numpy_int_out[0] = $self->nu();
    numpy_int_out[1] = $self->nv();
    numpy_int_out[2] = $self->nw();
  }
}


%extend NXmap {
  NXmap ( clipper::Grid const & grid, clipper::RTop<ftype> const & rtop)
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

    clipper::Coord_grid c;
    clipper::Grid map_grid = (*($self)).grid();
    
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
            numpy_array[i] = (*($self)).get_data(c);
    } else if (order == 'C') {
     for ( c.u() = 0; c.u() < top_u; c.u()++ )
        for ( c.v() = 0; c.v() < top_v; c.v()++ )
          for (  c.w() = 0; c.w() < top_w; c.w()++, i++ )
            numpy_array[i] = (*($self)).get_data(c);
    } else {
        throw std::invalid_argument("Order must be either F (Fortran-style wvu) or C (C-style uvw)");
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

    clipper::Coord_grid c;
    clipper::Grid map_grid = (*($self)).grid();

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
            (*($self)).set_data(c, numpy_3d_in[i]);
    } else if (order == 'C') {
      for ( c.u() = 0; c.u() < top_u; c.u()++ )
        for ( c.v() = 0; c.v() < top_v; c.v()++ )
          for (  c.w() = 0; c.w() < top_w; c.w()++, i++ )
            (*($self)).set_data(c, numpy_3d_in[i]);
    } 
    return i;
  }

  RTop<ftype> operator_orth_grid ()
  {
    return (*($self)).operator_orth_grid();
  }
  RTop<ftype> operator_grid_orth ()
  {
    return (*($self)).operator_grid_orth();
  }

}




%extend Xmap {

  clipper::Xmap<T>::Map_reference_coord get_reference_coord(const Coord_grid& pos)
  {

    clipper::Xmap<T>::Map_reference_coord ref((*($self)), pos);
    return ref;
  }

  /*! Return min, max, mean, sigma, skew and kurtosis for this map.
   *  Computationally expensive, so best to save the result on the Python side.
   */
  std::vector<double> stats ()
  {
    std::vector<double> ret;
    double max, min, mean, sd, skew, kurtosis;
    int n = 0;

    double sum = 0;
    double sum_sq = 0;
    double sum_3rd = 0;
    double sum_4th = 0;

    double rho_sq, rho_3rd, rho_4th;

    clipper::Xmap_base::Map_reference_index ix;
    ix = $self->first();
    max = min = (*self)[ix];
    for (ix = $self->first(); !ix.last(); ix.next()) {
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
      return ret;
    }
  }

  /*! Return list of all  grid points with multiplicity greater than 1
   *  in the asymmetric unit.
   */
  std::vector<std::vector<int> > special_positions ()
  {
    clipper::Xmap_base::Map_reference_index ix;
    std::vector<std::vector<int> > ret;
    clipper::Coord_grid this_coord;
    ix = $self->first();
    size_t mult;
    for (ix = $self->first(); !ix.last(); ix.next()) {
      mult = $self-> multiplicity(ix.coord());
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
    clipper::Coord_frac start_frac(frac_offset[0], frac_offset[1], frac_offset[2]);
    clipper::Grid_sampling grid = $self->grid_sampling();
    clipper::Coord_grid start_coord = start_frac.coord_grid(grid);
    clipper::Coord_grid this_coord;
    size_t mult;
    bool done = false;
    for (int i = 0; i < grid.nu(); i++) {
      for (int j = 0; j < grid.nv(); j++) {
        for (int k = 0; k < grid.nw(); k++) {
          this_coord = start_coord + clipper::Coord_grid(i,j,k);
          mult = $self-> multiplicity(this_coord);
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

  std::vector<std::vector<double> > special_positions_unit_cell_xyz (double frac_offset[3])
  {
    size_t count = 0;
    std::vector<std::vector<double> > ret;
    clipper::Coord_frac start_frac(frac_offset[0], frac_offset[1], frac_offset[2]);
    clipper::Grid_sampling grid = $self->grid_sampling();
    clipper::Coord_grid start_coord = start_frac.coord_grid(grid);
    clipper::Coord_grid this_coord;
    clipper::Cell cell = $self->cell();
    size_t mult;
    bool done = false;
    for (int i = 0; i < grid.nu(); i++) {
      for (int j = 0; j < grid.nv(); j++) {
        for (int k = 0; k < grid.nw(); k++) {
          this_coord = start_coord + clipper::Coord_grid(i,j,k);
          mult = $self-> multiplicity(this_coord);
          if (mult > 1) {
            std::vector<double> this_sp;
            clipper::Coord_frac this_frac = this_coord.coord_frac(grid);
            clipper::Coord_orth this_coord_orth = this_frac.coord_orth(cell);
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

    clipper::Coord_grid c;
    clipper::Grid map_grid = (*($self)).grid_asu();
    
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
              numpy_array[i] = (*($self)).get_data(c);
            else
              numpy_array[i] = 0.0;
          }
      return i;
    } else if (order == 'C') {
      for ( c.u() = 0; c.u() < top_u; c.u()++ )
        for ( c.v() = 0; c.v() < top_v; c.v()++ )
          for ( c.w() = 0; c.w() < nw; c.w()++, i++ ) {
            if ( c.u() < map_grid.nu() && c.v() < map_grid.nv() && c.w() < map_grid.nw() )
              numpy_array[i] = (*($self)).get_data(c);
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

    clipper::Coord_grid c;
    clipper::Grid map_grid = (*($self)).grid_asu();

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
            (*($self)).set_data(c, numpy_3d_in[i]);
      return i;
    } else if (order == 'C') {
      for ( c.u() = 0; c.u() < top_u; c.u()++ )
        for ( c.v() = 0; c.v() < top_v; c.v()++ )
          for (  c.w() = 0; c.w() < top_w; c.w()++, i++ )
            (*($self)).set_data(c, numpy_3d_in[i]);
      return i;
    } 
  }
  
  //! Export an arbitrary box as a numpy array

  int export_section_numpy ( double *numpy_array, int nu, int nv, int nw, clipper::Coord_grid& start, char order = 'F', std::string rot = "xyz" )
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
    
    clipper::Xmap_base::Map_reference_coord ix( (*($self)) );

    if (order == 'F') {
      for ( w = start.w(); w < maxw; w++ )
        for ( v = start.v(); v < maxv; v++ )
          for ( ix.set_coord(Coord_grid(start.u(),v,w)); ix.coord().u() < maxu; ix.next_u(), i++ ) {
            numpy_array[i] = (*($self))[ix];
          }
      return i;
    } else if (order == 'C') {
      for ( u = start.u(); u < maxu; u++ )
        for ( v = start.v(); v < maxv; v++ )
          for ( ix.set_coord(Coord_grid(u,v,start.w())); ix.coord().w() < maxw; ix.next_w(), i++ ) {
            numpy_array[i] = (*($self))[ix];
          }
      return i;
    }
  }
  
  //! Import an arbitrary box as a numpy array

  int import_section_numpy ( double *numpy_3d_in, int nu, int nv, int nw, clipper::Coord_grid& start, char order = 'F', std::string rot = "xyz" )
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
    clipper::Xmap_base::Map_reference_coord ix( (*($self)) );

    if (order == 'F') {
      for ( w = start.w(); w < maxw; w++ )
        for ( v = start.v(); v < maxv; v++ )
          for ( ix.set_coord(Coord_grid(start.u(),v,w)); ix.coord().u() < maxu; ix.next_u(), i++ ) {
            (*($self))[ix] = numpy_3d_in[i];
          }
      return i;
    } else if (order == 'C') {
      for ( u = start.u(); u < maxu; u++ )
        for ( v = start.v(); v < maxv; v++ )
          for ( ix.set_coord(Coord_grid(u,v,start.w())); ix.coord().w() < maxw; ix.next_w(), i++ ) {
            (*($self))[ix] = numpy_3d_in[i];
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
    
    const Cell& cell = $self->cell();
    
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
              numpy_array[count] = $self->interp<Interp_cubic>(thiscoord);
            } else {
              numpy_array[count] = $self->interp<Interp_linear>(thiscoord);
            }
          }
        }
      }
    } else if (order == 'C') {
      for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
          for (int k = 0; k < nz; k++, count++) {
            thiscoord = (origin + Coord_orth(x_inc*i, y_inc*j, z_inc*k)).coord_frac(cell);
            if (mode_cubic) {
              numpy_array[count] = $self->interp<Interp_cubic>(thiscoord);
            } else {
              numpy_array[count] = $self->interp<Interp_linear>(thiscoord);
            }
          }
        }
      }
    } else {
      throw std::invalid_argument("Order must be either F (Fortran-style wvu) or C (C-style uvw)");
    }
    return count;
      
  }
  
  
  
  //! Get the length of a voxel along each axis in Angstroms

  void voxel_size(double numpy_double_out[3])
  {
    clipper::Grid_sampling g = $self->grid_sampling();
    clipper::Cell c = $self->cell();
    numpy_double_out[0] = c.a() / g.nu();
    numpy_double_out[1] = c.b() / g.nv();
    numpy_double_out[2] = c.c() / g.nw();
  }
  
  //! Get the length of a voxel along each axis in fractional coordinates

  void voxel_size_frac(double numpy_double_out[3])
  {
    clipper::Grid_sampling g = $self->grid_sampling();
    numpy_double_out[0] = 1.0 / g.nu();
    numpy_double_out[1] = 1.0 / g.nv();
    numpy_double_out[2] = 1.0 / g.nw();
  }
  
  //! Return the interpolated density value at a given fractional coordinate

  T interp_cubic_frac_coord(Coord_frac f)
  {
    return $self->interp<Interp_cubic>(f);
  }
  
  //! Return the interpolated density value at a given (x,y,z) coordinate
  
  T interp_cubic_xyz(double numpy_1d_in[3])
  {
    Coord_frac thecoord = Coord_orth(numpy_1d_in[0], numpy_1d_in[1], numpy_1d_in[2]).coord_frac($self->cell());
    return $self->interp<Interp_cubic>(thecoord);
  }
  
  //! Return the interpolated density value at a given fractional coordinate

  T interp_linear_frac_coord(Coord_frac f)
  {
    return $self->interp<Interp_linear>(f);
  }

  //! Return the interpolated density value at a given (x,y,z) coordinate
  
  T interp_linear_xyz(double numpy_1d_in[3])
  {
    Coord_frac thecoord = Coord_orth(numpy_1d_in[0], numpy_1d_in[1], numpy_1d_in[2]).coord_frac($self->cell());
    return $self->interp<Interp_linear>(thecoord);
  }



  RTop<ftype> operator_orth_grid ()
  {
    return (*($self)).operator_orth_grid();
  }
  RTop<ftype> operator_grid_orth ()
  {
    return (*($self)).operator_grid_orth();
  }

}


%extend Xmap<double> {

  Xmap_reference_coord map_reference_coord( const Coord_grid& coord )
  {
    Xmap_reference_coord ret(*($self), coord);
    return ret;
  }

  const int& symop_of_reference_coord (const Xmap_reference_coord& coord )
  {
    return coord.sym();
  }
  /*! Find and add the necessary translations for each symop in order
   *  to pack the same unit cell as a given fractional coordinate. Return a
   *  Unit_Cell object.
   */

  clipper::Unit_Cell unit_cell_symops (clipper::Coord_frac ref, const clipper::Atom_list& atoms)
  {
    return clipper::Unit_Cell(ref, atoms, $self->cell(), $self->spacegroup(), $self->grid_sampling());
  }



}


}





namespace clipper
{
%extend Xmap<float> {
  void fft_from (const clipper::HKL_data<clipper::data32::F_phi> &fb)
  {
    ($self)->fft_from( fb );
  }

  void fft_to (clipper::HKL_data<clipper::data32::F_phi> &fphidata)
  {
    ($self)->fft_to(fphidata, clipper::Xmap_base::Default);
  }

}
%extend Xmap<double> {
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
namespace clipper
{
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
    numpy_double_out[0] = $self->range().min();
    numpy_double_out[1] = $self->range().max();
  }
};
}

%include "../clipper/core/map_utils.h"


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
}

%extend MModel {
  MPolymer& __getitem__(int i)
  {
    int array_len = $self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    return (*($self))[i];
fail:
    return (*($self))[0];
  }
  void __setitem__(int i, MPolymer& mpol)
  {
    int array_len = $self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    (*($self))[i]=mpol;
    return;
fail:
    return;
  }
  size_t __len__()
  {
    return ($self)->size();
  }
}


%extend MPolymer {
  MMonomer& __getitem__(int i)
  {
    int array_len = $self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    return (*($self))[i];
  }
  void __setitem__(int i, MMonomer& mmon)
  {
    int array_len = $self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    (*($self))[i]=mmon;
    return;
  }
  size_t __len__()
  {
    return ($self)->size();
  }
}


%extend MMonomer {
  MAtom& __getitem__(int i)
  {
    int array_len = $self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    return (*($self))[i];
  }
  void __setitem__(int i, MAtom& atom)
  {
    int array_len = $self->size();
    i = (i < 0) ? array_len + i : i;
    if (i >= array_len || i < 0) {
      throw std::out_of_range("");
    }
    (*($self))[i]=atom;
    return;
  }
  size_t __len__()
  {
    return ($self)->size();
  }
}

%extend MAtom {
  std::string __str__( )
  {
    return (*($self)).id() + " " + (*($self)).coord_orth().format();
  }
}

%extend U_aniso_orth {
  U_aniso_orth __add__(const U_aniso_orth& u2) { return *self + u2; }
  U_aniso_orth __sub__(const U_aniso_orth& u2) { return *self + (-u2); }
  U_aniso_orth __neg__() { return -(*self); }
  U_aniso_orth __mul__(const double& f) { return f * (*self); }
  U_aniso_orth __rmul__(const double& f) { return f * (*self); }
  void _get_vals(double numpy_double_out[6])
  {
    numpy_double_out[0] = (*self).mat00();
    numpy_double_out[1] = (*self).mat11();
    numpy_double_out[2] = (*self).mat22();
    numpy_double_out[3] = (*self).mat01();
    numpy_double_out[4] = (*self).mat02();
    numpy_double_out[5] = (*self).mat12();
  }
  
}
%extend U_aniso_frac {
  U_aniso_frac __add__(const U_aniso_frac& u2) { return *self + u2; }
  U_aniso_frac __sub__(const U_aniso_frac& u2) { return *self + (-u2); }
  U_aniso_frac __neg__() { return -(*self); }
  U_aniso_frac __mul__(const double& f) { return f * (*self); }
  U_aniso_frac __rmul__(const double& f) { return f * (*self); }
  void _get_vals(double numpy_double_out[6])
  {
    numpy_double_out[0] = (*self).mat00();
    numpy_double_out[1] = (*self).mat11();
    numpy_double_out[2] = (*self).mat22();
    numpy_double_out[3] = (*self).mat01();
    numpy_double_out[4] = (*self).mat02();
    numpy_double_out[5] = (*self).mat12();
  }
  
}

}

%include "../clipper/mmdb/clipper_mmdb.h"
%include "../clipper/minimol/minimol.h"
%include "../clipper/minimol/minimol_io.h"
%include "../clipper/minimol/minimol_seq.h"


namespace clipper
{
namespace data64
{
%extend Flag_bool {
  bool get_flag()
  {
    bool theFlag = ($self)->flag();
    return theFlag;
  }
  void set_flag(bool theFlag)
  {
    ($self)->flag() = theFlag;
  }
}
%extend Flag {
  int get_flag()
  {
    int theFlag = ($self)->flag();
    return theFlag;
  }
  void set_flag(int theFlag)
  {
    ($self)->flag() = theFlag;
  }
}
}
namespace data32
{
%extend Flag_bool {
  bool get_flag()
  {
    bool theFlag = ($self)->flag();
    return theFlag;
  }
  void set_flag(bool theFlag)
  {
    ($self)->flag() = theFlag;
  }
}
%extend Flag {
  int get_flag()
  {
    int theFlag = ($self)->flag();
    return theFlag;
  }
  void set_flag(int theFlag)
  {
    ($self)->flag() = theFlag;
  }
}
}
namespace datatypes
{
%extend Flag_bool {
  bool get_flag()
  {
    bool theFlag = ($self)->flag();
    return theFlag;
  }
  void set_flag(bool theFlag)
  {
    ($self)->flag() = theFlag;
  }
  clipper::datatypes::Flag_bool copy()
  {
    clipper::datatypes::Flag_bool ret;
    ret = *($self);
    return ret;
  }
}
%extend Flag {
  int get_flag()
  {
    int theFlag = ($self)->flag();
    return theFlag;
  }
  void set_flag(int theFlag)
  {
    ($self)->flag() = theFlag;
  }
  clipper::datatypes::Flag copy()
  {
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
  %extend HKL_data {
    HKL_data<clipper::datatypes::Flag_bool> not_()
    {
      return !(*($self));
    }
    HKL_data<clipper::datatypes::Flag_bool> __or__(const HKL_data<T> &d1)
    {
      return (*($self)) | d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __xor__(const HKL_data<T> &d1)
    {
      return (*($self)) ^ d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __and__(const HKL_data<T> &d1)
    {
      return (*($self)) & d1;
    }
    
  }
  %extend datatypes::ABCD {
    void vals(double numpy_double_out[4]) {
      numpy_double_out[0] = $self->a();
      numpy_double_out[1] = $self->b();
      numpy_double_out[2] = $self->c();
      numpy_double_out[3] = $self->d();
    }  
  }





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

%extend HKL_data<clipper::data32::Flag_bool> {
  HKL_data<clipper::datatypes::Flag_bool>  copy()
  {
    HKL_data<clipper::data32::Flag_bool> ret;
    ret = *($self);
    return ret;
  }
}

%extend HKL_data<clipper::data32::Flag> {
  HKL_data<clipper::datatypes::Flag_bool> __eq__(const int& n)
  {
    return (*($self)) == n;
  }
  HKL_data<clipper::datatypes::Flag_bool> __ne__(const int& n)
  {
    return (*($self)) != n;
  }
  HKL_data<clipper::datatypes::Flag_bool> __ge__(const int& n)
  {
    return (*($self)) >= n;
  }
  HKL_data<clipper::datatypes::Flag_bool> __le__(const int& n)
  {
    return (*($self)) <= n;
  }
  HKL_data<clipper::datatypes::Flag_bool> __gt__(const int& n)
  {
    return (*($self)) > n;
  }
  HKL_data<clipper::datatypes::Flag_bool> __lt__(const int& n)
  {
    return (*($self)) < n;
  }
  HKL_data<clipper::datatypes::Flag>  copy()
  {
    HKL_data<clipper::data32::Flag> ret;
    ret = *($self);
    return ret;
  }
}

%extend datatypes::F_sigF<float> {
  clipper::datatypes::F_sigF<float>  copy()
  {
    clipper::data32::F_sigF ret;
    ret = *($self);
    return ret;
  }
}

%extend datatypes::F_sigF_ano<float> {
  clipper::datatypes::F_sigF_ano<float>  copy()
  {
    clipper::data32::F_sigF_ano ret;
    ret = *($self);
    return ret;
  }
}

%extend datatypes::I_sigI<float> {
  clipper::datatypes::I_sigI<float>  copy()
  {
    clipper::data32::I_sigI ret;
    ret = *($self);
    return ret;
  }
}

%extend datatypes::E_sigE<float> {
  clipper::datatypes::E_sigE<float>  copy()
  {
    clipper::data32::E_sigE ret;
    ret = *($self);
    return ret;
  }
}

%extend datatypes::F_phi<float> {
  clipper::datatypes::F_phi<float>  __add__(const clipper::datatypes::F_phi<float> &h2)
  {
    clipper::data32::F_phi ret;
    ret = *($self)+h2;
    return ret;
  }
  clipper::datatypes::F_phi<float>  __sub__(const clipper::datatypes::F_phi<float> &h2)
  {
    clipper::data32::F_phi ret;
    ret = *($self)-h2;
    return ret;
  }
  clipper::datatypes::F_phi<float>  __neg__()
  {
    clipper::data32::F_phi ret;
    ret = -*($self);
    return ret;
  }
  clipper::datatypes::F_phi<float>  copy()
  {
    clipper::data32::F_phi ret;
    ret = *($self);
    return ret;
  }
}

%extend datatypes::ABCD<float> {
  clipper::datatypes::ABCD<float>  __add__(const clipper::datatypes::ABCD<float> &h2)
  {
    clipper::data32::ABCD ret;
    ret = *($self)+h2;
    return ret;
  }
  clipper::datatypes::ABCD<float>  copy()
  {
    clipper::data32::ABCD ret;
    ret = *($self);
    return ret;
  }
}

%extend HKL_data<clipper::data32::ABCD> {
  HKL_data<clipper::datatypes::ABCD<float> > __add__(const HKL_data<clipper::datatypes::ABCD<float> > &h2)
  {
    HKL_data<clipper::data32::ABCD> ret;
    ret = *($self)+h2;
    return ret;
  }
  HKL_data<clipper::datatypes::ABCD<float> > copy()
  {
    HKL_data<clipper::data32::ABCD> ret;
    ret = *($self);
    return ret;
  }
}

%extend HKL_data<clipper::data32::F_phi> {
  HKL_data<clipper::datatypes::F_phi<float> >  copy()
  {
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
  HKL_data<clipper::datatypes::F_phi<float> > __add__(const HKL_data<clipper::datatypes::F_phi<float> > &h2)
  {
    HKL_data<clipper::data32::F_phi> ret;
    ret = *($self)+h2;
    return ret;
  }
  HKL_data<clipper::datatypes::F_phi<float> > __sub__(const HKL_data<clipper::datatypes::F_phi<float> > &h2)
  {
    HKL_data<clipper::datatypes::F_phi<float> > ret;
    ret = *($self)-h2;
    return ret;
  }
  HKL_data<clipper::datatypes::F_phi<float> > __neg__()
  {
    HKL_data<clipper::datatypes::F_phi<float> > ret;
    ret = -*($self);
    return ret;
  }
  HKL_data<clipper::datatypes::F_phi<float> > __mul__(const float s)
  {
    HKL_data<clipper::datatypes::F_phi<float> > ret;
    ret = *($self)*s;
    return ret;
  }
  HKL_data<clipper::datatypes::F_phi<float> > __rmul__(const float s)
  {
    HKL_data<clipper::datatypes::F_phi<float> > ret;
    ret = *($self)*s;
    return ret;
  }
}

%extend HKL_data<clipper::data32::E_sigE> {
  void scaleBySqrtResolution(const clipper::ResolutionFn &escale)
  {
    for ( clipper::HKL_data_base::HKL_reference_index ih = (*($self)).first(); !ih.last(); ih.next() )
      if ( !(*($self))[ih].missing() ) (*($self))[ih].scale( sqrt( escale.f(ih) ) );
  }
  void scaleByResolution(const clipper::ResolutionFn &escale)
  {
    for ( clipper::HKL_data_base::HKL_reference_index ih = (*($self)).first(); !ih.last(); ih.next() )
      if ( !(*($self))[ih].missing() ) (*($self))[ih].scale( escale.f(ih) );
  }
  HKL_data<clipper::datatypes::E_sigE<float> >  copy()
  {
    HKL_data<clipper::data32::E_sigE> ret;
    ret = *($self);
    return ret;
  }
}

%extend HKL_data<clipper::data32::ABCD> {
  void compute_from_phi_fom(const HKL_data< clipper::datatypes::Phi_fom<float> > &phiw)
  {
    ($self)->compute( phiw, clipper::data32::Compute_abcd_from_phifom() );
  }
  void compute_add_abcd(const HKL_data< clipper::datatypes::ABCD<float> > &abcd1,
  const HKL_data< clipper::datatypes::ABCD<float> > &abcd2)
  {
    ($self)->compute( abcd1, abcd2, clipper::data32::Compute_add_abcd() );
  }
}

%extend HKL_data<clipper::data32::Phi_fom> {
  void compute_from_abcd(const HKL_data< clipper::datatypes::ABCD<float> > &abcd)
  {
    ($self)->compute( abcd, clipper::data32::Compute_phifom_from_abcd() );
  }
  HKL_data<clipper::datatypes::Phi_fom<float> >  copy()
  {
    HKL_data<clipper::data32::Phi_fom> ret;
    ret = *($self);
    return ret;
  }
}

%extend HKL_data<clipper::data32::F_sigF> {
  void compute_mean_from_fano(const HKL_data< clipper::datatypes::F_sigF_ano<float> > &fano)
  {
    ($self)->compute( fano, clipper::data32::Compute_mean_fsigf_from_fsigfano() );
  }
  void compute_diff_from_fano(const HKL_data< clipper::datatypes::F_sigF_ano<float> > &fano)
  {
    ($self)->compute( fano, clipper::data32::Compute_diff_fsigf_from_fsigfano() );
  }
  void compute_scale_u_iso_fsigf(float scale, float u_value,
  const HKL_data< clipper::datatypes::F_sigF<float> > &fsigf )
  {
    ($self)->compute( fsigf, clipper::data32::Compute_scale_u_iso_fsigf(scale, u_value) );
  }
  void compute_scale_u_aniso_fsigf(float scale, clipper::U_aniso_orth u_value,
  const HKL_data< clipper::datatypes::F_sigF<float> > &fsigf )
  {
    ($self)->compute( fsigf, clipper::data32::Compute_scale_u_aniso_fsigf(scale, u_value) );
  }
  HKL_data<clipper::datatypes::F_sigF<float> >  copy()
  {
    HKL_data<clipper::data32::F_sigF> ret;
    ret = *($self);
    return ret;
  }
}

%extend HKL_data<clipper::data32::F_sigF_ano> {
  void compute_scale_u_iso_fsigfano(float scale, float u_value,
  const HKL_data< clipper::datatypes::F_sigF_ano<float> > &fsigfano )
  {
    ($self)->compute( fsigfano, clipper::data32::Compute_scale_u_iso_fsigfano(scale, u_value) );
  }
  void compute_scale_u_aniso_fsigfano(float scale, clipper::U_aniso_orth u_value,
  const HKL_data< clipper::datatypes::F_sigF_ano<float> > &fsigfano )
  {
    ($self)->compute( fsigfano, clipper::data32::Compute_scale_u_aniso_fsigfano(scale, u_value) );
  }
  HKL_data<clipper::datatypes::F_sigF_ano<float> >  copy()
  {
    HKL_data<clipper::data32::F_sigF_ano> ret;
    ret = *($self);
    return ret;
  }
}

%extend HKL_data<clipper::data32::I_sigI> {
  void compute_scale_u_iso_isigi(float scale, float u_value,
  const HKL_data< clipper::datatypes::I_sigI<float> > &isigi )
  {
    ($self)->compute( isigi, clipper::data32::Compute_scale_u_iso_isigi(scale, u_value) );
  }
  void compute_scale_u_aniso_isigi(float scale, clipper::U_aniso_orth u_value,
  const HKL_data< clipper::datatypes::I_sigI<float> > &isigi )
  {
    ($self)->compute( isigi, clipper::data32::Compute_scale_u_aniso_isigi(scale, u_value) );
  }
  HKL_data<clipper::datatypes::I_sigI<float> >  copy()
  {
    HKL_data<clipper::data32::I_sigI> ret;
    ret = *($self);
    return ret;
  }
}



%extend HKL_data< clipper::datatypes::F_phi<float> > {

  void getDataNumpy(double *test_numpy_a, int test_numpy_n)
  {
    int i=0;
    for(clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next(),i++ ) {
      if(!((*($self))[ih].missing())) {
        std::vector<xtype> thisData((*($self)).data_size());
        (*($self)).data_export(ih.hkl(),&(thisData[0]));
        std::vector<float> thisDataf((*($self)).data_size());
        for(unsigned idat=0; idat<(*($self)).data_size(); ++idat) {
          test_numpy_a[i*(*($self)).data_size()+idat] = thisData[idat];
        }
      } else {
        for(unsigned idat=0; idat<(*($self)).data_size(); ++idat) {
          test_numpy_a[i*(*($self)).data_size()+idat] = std::numeric_limits<float>::quiet_NaN();
        }
      }
    }
  }

  std::vector<std::vector<float> > getData()
  {
    std::vector<std::vector<float> > allData;
    for(clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      if(!((*($self))[ih].missing())) {
        std::vector<xtype> thisData((*($self)).data_size());
        (*($self)).data_export(ih.hkl(),&(thisData[0]));
        std::vector<float> thisDataf((*($self)).data_size());
        for(unsigned idat=0; idat<(*($self)).data_size(); ++idat) {
          thisDataf[idat] = thisData[idat];
        }
        allData.push_back(thisDataf);
      }
    }
    return allData;
  }
  clipper::datatypes::F_phi<float>& __getitem__(int i)
  {
    int sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*($self))[i];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    return sz;
  }
}

%extend HKL_data<clipper::data32::F_phi> {
  void compute_neg(const HKL_data< clipper::datatypes::F_phi<float> > &fphi )
  {
    ($self)->compute( fphi, clipper::data32::Compute_neg_fphi() );
  }
  void compute_add_fphi(const HKL_data< clipper::datatypes::F_phi<float> > &fphi1,
  const HKL_data< clipper::datatypes::F_phi<float> > &fphi2)
  {
    ($self)->compute( fphi1, fphi2, clipper::data32::Compute_add_fphi() );
  }
  void compute_sub_fphi(const HKL_data< clipper::datatypes::F_phi<float> > &fphi1,
  const HKL_data< clipper::datatypes::F_phi<float> > &fphi2)
  {
    ($self)->compute( fphi1, fphi2, clipper::data32::Compute_sub_fphi() );
  }
  void compute_from_fsigf_phifom(const HKL_data< clipper::datatypes::F_sigF<float> > &fsigf,
  const HKL_data< clipper::datatypes::Phi_fom<float> > &phifom )
  {
    ($self)->compute( fsigf, phifom, clipper::data32::Compute_fphi_from_fsigf_phifom() );
  }
  void compute_scale_u_iso_fphi(float scale, float u_value,
  const HKL_data< clipper::datatypes::F_phi<float> > &fphi )
  {
    ($self)->compute( fphi, clipper::data32::Compute_scale_u_iso_fphi(scale, u_value) );
  }
  void compute_scale_u_aniso_fphi(float scale, clipper::U_aniso_orth u_value,
  const HKL_data< clipper::datatypes::F_phi<float> > &fphi )
  {
    ($self)->compute( fphi, clipper::data32::Compute_scale_u_aniso_fphi(scale, u_value) );
  }

  void getDataNumpy(double *test_numpy_a, int test_numpy_n)
  {
    int i=0;
    for(clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next(),i++ ) {
      if(!((*($self))[ih].missing())) {
        std::vector<xtype> thisData((*($self)).data_size());
        (*($self)).data_export(ih.hkl(),&(thisData[0]));
        std::vector<float> thisDataf((*($self)).data_size());
        for(unsigned idat=0; idat<(*($self)).data_size(); ++idat) {
          test_numpy_a[i*(*($self)).data_size()+idat] = thisData[idat];
        }
      } else {
        for(unsigned idat=0; idat<(*($self)).data_size(); ++idat) {
          test_numpy_a[i*(*($self)).data_size()+idat] = std::numeric_limits<float>::quiet_NaN();
        }
      }
    }
  }

  std::vector<std::vector<float> > getData()
  {
    std::vector<std::vector<float> > allData;
    for(clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      if(!((*($self))[ih].missing())) {
        std::vector<xtype> thisData((*($self)).data_size());
        (*($self)).data_export(ih.hkl(),&(thisData[0]));
        std::vector<float> thisDataf((*($self)).data_size());
        for(unsigned idat=0; idat<(*($self)).data_size(); ++idat) {
          thisDataf[idat] = thisData[idat];
        }
        allData.push_back(thisDataf);
      }
    }
    return allData;
  }
  clipper::data32::F_phi& __getitem__(int i)
  {
    int sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*($self))[i];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    return sz;
  }
}

%extend HKL_data<clipper::data32::E_sigE> {
  void compute_from_fsigf(const HKL_data< clipper::datatypes::F_sigF<float> > &fsigf )
  {
    ($self)->compute( fsigf, clipper::data32::Compute_EsigE_from_FsigF() );
  }
}

%extend HKL_data<clipper::datatypes::Flag_bool> {
  clipper::datatypes::Flag_bool& __getitem__(int i)
  {
    int sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*($self))[i];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    return sz;
  }
}
%extend HKL_data<clipper::data32::Flag_bool> {
  clipper::data32::Flag_bool& __getitem__(int i)
  {
    int sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*($self))[i];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    return sz;
  }
}

%extend HKL_data<clipper::datatypes::Flag> {
  clipper::datatypes::Flag& __getitem__(int i)
  {
    int sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*($self))[i];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    return sz;
  }
}
%extend HKL_data<clipper::data32::Flag> {
  clipper::data32::Flag& __getitem__(int i)
  {
    int sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*($self))[i];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    return sz;
  }
}


%extend HKL_data< clipper::datatypes::F_sigF<float> > {

  void getDataNumpy(double *test_numpy_a, int test_numpy_n)
  {
    int i=0;
    for(clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next(),i++ ) {
      if(!((*($self))[ih].missing())) {
        std::vector<xtype> thisData((*($self)).data_size());
        (*($self)).data_export(ih.hkl(),&(thisData[0]));
        std::vector<float> thisDataf((*($self)).data_size());
        for(unsigned idat=0; idat<(*($self)).data_size(); ++idat) {
          test_numpy_a[i*(*($self)).data_size()+idat] = thisData[idat];
        }
      } else {
        for(unsigned idat=0; idat<(*($self)).data_size(); ++idat) {
          test_numpy_a[i*(*($self)).data_size()+idat] = std::numeric_limits<float>::quiet_NaN();
        }
      }
    }
  }

  std::vector<std::vector<float> > getData()
  {
    std::vector<std::vector<float> > allData;
    for(clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      if(!((*($self))[ih].missing())) {
        std::vector<xtype> thisData((*($self)).data_size());
        (*($self)).data_export(ih.hkl(),&(thisData[0]));
        std::vector<float> thisDataf((*($self)).data_size());
        for(unsigned idat=0; idat<(*($self)).data_size(); ++idat) {
          thisDataf[idat] = thisData[idat];
        }
        allData.push_back(thisDataf);
      }
    }
    return allData;
  }
  clipper::datatypes::F_sigF<float>& __getitem__(int i)
  {
    int sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*($self))[i];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    return sz;
  }
}

%extend HKL_data< clipper::data32::F_sigF<float> > {
  clipper::data32::F_sigF<float>& __getitem__(int i)
  {
    int sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*($self))[i];
fail:
    return (*($self))[0];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ) {
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
}
%include "../clipper/contrib/sfweight.h"
namespace clipper
{
%template(SFweight_spline_float) SFweight_spline<float>;
}

%include "../clipper/contrib/sfcalc.h"
namespace clipper
{
%template(SFcalc_iso_sum_float) SFcalc_iso_sum<float>;
%template(SFcalc_aniso_sum_float) SFcalc_aniso_sum<float>;
%template(SFcalc_iso_fft_float) SFcalc_iso_fft<float>;
%template(SFcalc_aniso_fft_float) SFcalc_aniso_fft<float>;
}

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
    return (*($self))( xmap, atoms );
  }
  bool compute( NXmap<float>& xmap, const Atom_list& atoms ) const {
    return (*($self))( xmap, atoms );
  }
}
%extend EDcalc_iso<float> {
  bool compute( Xmap<float>& xmap, const Atom_list& atoms ) const {
    return (*($self))( xmap, atoms );
  }
  bool compute( NXmap<float>& xmap, const Atom_list& atoms ) const {
    return (*($self))( xmap, atoms );
  }
}
%extend EDcalc_aniso<float> {
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



