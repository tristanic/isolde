%module clipper

%include "std_vector.i"
%include "std_string.i"
%include "exception.i"
%include "std_except.i"

//%rename(clipper_Range) clipper::Range;
//%rename(clipper_Batch) clipper::datatypes::Batch;

%rename(MMonomer_TYPE) clipper::MMonomer::TYPE;
%rename(MMDBManager_TYPE) clipper::MMDBManager::TYPE;

%ignore clipper::CCP4MTZfile::assigned_paths;
%ignore clipper::MMDBManager::write_file(const String&);

%{
    #include <string>
    #include "../clipper/core/clipper_types.h"
    #include "../clipper/core/hkl_lookup.h"
    #include "../clipper/core/hkl_info.h"
    #include "../clipper/core/xmap.h"
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
    #include "../clipper/ccp4/ccp4_utils.h"
    #include "../clipper/core/map_utils.h"
    #include "../clipper/core/clipper_stats.h"
    #include "../clipper/mmdb/clipper_mmdb.h"
    #include "../clipper/core/hkl_datatypes.h"
    #include "../clipper/core/hkl_data.h"
    #include "../clipper/core/hkl_compute.h"
    #include "../clipper/cif/cif_data_io.h"
    #include "../clipper/contrib/sfcalc_obs.h"
    #include "../clipper/contrib/sfweight.h"
    #include "../clipper/contrib/sfcalc.h"
    #include "../clipper/contrib/edcalc.h"
    #include "../clipper/core/nxmap_operator.h"
    #include "../clipper/contrib/convolution_search.h"
    #include "../clipper/contrib/mapfilter.h"
    #include "../clipper/contrib/originmatch.h"

    namespace clipper 
    {
        static int myErr = 0; // flag to save error state
    }

    using namespace clipper;
%}


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



namespace std 
{
   %template(UnsignedIntVector) vector<unsigned int>;
   %template(IntVector) vector<int>;
   %template(IntIntVector) vector<vector<int> >;
   %template(DoubleVector) vector<double>;
   %template(DoubleDoubleVector) vector<vector<double> >;
   %template(ClipperStringVector) vector<clipper::String>;
   %template(StringVector) vector<string>;
}

%apply std::string { clipper::String }
%apply std::string& { clipper::String& }

// FIXME - this one breaks compile, but it is necessary
//%apply const std::string& { const clipper::String& } 

%apply std::string { String }
%apply std::string& { String& }
//%apply const std::string& { const String& } 

%typemap(in) (const clipper::String&)
{
   std::string ss = PyString_AsString($input);
   clipper::String s(ss);
   $1 = &(s);
}

%include "../clipper/core/clipper_types.h"
%include "../clipper/core/clipper_util.h"

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


//%template(vec3_int) clipper::Vec3< int >;
//%template(vec3_float) clipper::Vec3< ftype32 >;
//%template(vec3_double) clipper::Vec3< ftype64 >;

%ignore Util::Vec3<>;

%include "../clipper/core/spacegroup.h"

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
      {	ftype h = ftype(v[0]); ftype k = ftype(v[1]); ftype l = ftype(v[2]);
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


%typemap(in) (const clipper::String&)
{
   std::string ss = PyString_AsString($input);
   clipper::String s(ss);
   $1 = &(s);
}

%include "../clipper/core/coords.h"

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
%include "../clipper/mmdb/clipper_mmdb.h"
%include "../clipper/minimol/minimol.h"
%include "../clipper/minimol/minimol_io.h"


namespace clipper {
  %extend Xmap<float> {
    void fft_from_float(const clipper::HKL_data<clipper::data32::F_phi> &fb){
      ($self)->fft_from( fb );
    }
  }
  %extend Xmap<double> {
    void fft_from_double(const clipper::HKL_data<clipper::data64::F_phi> &fb){
      ($self)->fft_from( fb );
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
%include "../clipper/ccp4/ccp4_mtz_io.h"
%include "../clipper/ccp4/ccp4_mtz_types.h"
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

%include "../clipper/core/map_utils.h"

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
  }
  %exception Atom_list::__getitem__ {
    assert(!myErr);
    $action
    if (myErr) {
      myErr = 0; // clear flag for next time
      SWIG_exception(SWIG_IndexError, "Index out of bounds");
    }
  }

  %extend Atom_list {
    Atom __getitem__(size_t i) { 
      if (i >= $self->size()) {
        myErr = 1;
        return (*($self))[0];
      }
      return (*($self))[i];
      fail:
        return (*($self))[0];
    }
    size_t __len__() { 
      return ($self)->size();
    }
  }
  %exception MModel::__getitem__ {
    assert(!myErr);
    $action
    if (myErr) {
      myErr = 0; // clear flag for next time
      SWIG_exception(SWIG_IndexError, "Index out of bounds");
    }
  }

  %extend MModel {
    MPolymer __getitem__(size_t i) { 
      if (i >= $self->size()) {
        myErr = 1;
        return (*($self))[0];
      }
      return (*($self))[i];
      fail:
        return (*($self))[0];
    }
    size_t __len__() { 
      return ($self)->size();
    }
  }
  %exception MPolymer::__getitem__ {
    assert(!myErr);
    $action
    if (myErr) {
      myErr = 0; // clear flag for next time
      SWIG_exception(SWIG_IndexError, "Index out of bounds");
    }
  }

  %extend MPolymer {
    MMonomer __getitem__(size_t i) { 
      if (i >= $self->size()) {
        myErr = 1;
        return (*($self))[0];
      }
      return (*($self))[i];
      fail:
        return (*($self))[0];
    }
    size_t __len__() { 
      return ($self)->size();
    }
  }
  %exception MMonomer::__getitem__ {
    assert(!myErr);
    $action
    if (myErr) {
      myErr = 0; // clear flag for next time
      SWIG_exception(SWIG_IndexError, "Index out of bounds");
    }
  }

  %extend MMonomer {
    Atom __getitem__(size_t i) { 
      if (i >= $self->size()) {
        myErr = 1;
        return (*($self))[0];
      }
      return (*($self))[i];
      fail:
        return (*($self))[0];
    }
    size_t __len__() { 
      return ($self)->size();
    }
  }
}

%include "../clipper/core/hkl_datatypes.h"
namespace clipper {
  namespace data64 {
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
}

namespace clipper 
{
    
    %template(F_sigF_float) clipper::datatypes::F_sigF<float>;
    %template(F_sigF_double) clipper::datatypes::F_sigF<double>;
    %template(HKL_data_F_sigF_float) HKL_data< clipper::data32::F_sigF >;
    %template(HKL_data_F_sigF_double) HKL_data< clipper::data64::F_sigF >;

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
    
    %exception HKL_data_Flag::__getitem__ {
    assert(!myErr);
    $action
    if (myErr) {
      myErr = 0; // clear flag for next time
      SWIG_exception(SWIG_IndexError, "Index out of bounds");
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
  }
  
  %extend HKL_data<clipper::data32::E_sigE> {
    void compute_from_fsigf(const HKL_data< clipper::datatypes::F_sigF<float> > &fsigf ) {
      ($self)->compute( fsigf, clipper::data32::Compute_EsigE_from_FsigF() );
    }
  }
 
    %extend HKL_data<clipper::datatypes::Flag> {
    clipper::datatypes::Flag __getitem__(size_t i) { 
      size_t sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      if (i >= sz) {
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
  %exception HKL_data< clipper::datatypes::F_sigF<float> >::__getitem__ {
    assert(!myErr);
    $action
    if (myErr) {
      myErr = 0; // clear flag for next time
      SWIG_exception(SWIG_IndexError, "Index out of bounds");
    }
  }

  %extend HKL_data< clipper::datatypes::F_sigF<float> > {
    clipper::datatypes::F_sigF<float> __getitem__(size_t i) { 
      size_t sz=0;
      for ( clipper::HKL_data_base::HKL_reference_index ih = ($self)->first(); !ih.last(); ih.next() ){
        sz++;
      }
      if (i >= sz) {
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



%include "../clipper/cif/cif_data_io.h"
%include "../clipper/contrib/sfcalc_obs.h"

%include "../clipper/core/hkl_compute.h"
namespace clipper
{
    %template(SFcalc_obs_bulk_float) SFcalc_obs_bulk<float>;    
}

%include "../clipper/contrib/function_object_bases.h"
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
%}
void SetFlagBoth(clipper::HKL_data<clipper::data32::Flag> &flag);
void SetFlagBothIfMissing(clipper::HKL_data<clipper::data32::Flag> &flag, const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &myfsigf, const clipper::HKL_data< clipper::datatypes::Flag > &status, int freeflag);

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
