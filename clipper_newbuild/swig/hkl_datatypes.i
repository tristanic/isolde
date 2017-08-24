%include "../clipper/core/hkl_datatypes.h"

namespace clipper
{
namespace data64
{
%extend Flag_bool {
  bool get_flag()
  {
    bool theFlag = self->flag();
    return theFlag;
  }
  void set_flag(bool theFlag)
  {
    self->flag() = theFlag;
  }
} // extend Flag_bool
%extend Flag {
  int get_flag()
  {
    int theFlag = self->flag();
    return theFlag;
  }
  void set_flag(int theFlag)
  {
    self->flag() = theFlag;
  }
} // extend Flag
} // namespace data64
namespace data32
{
%extend Flag_bool {
  bool get_flag()
  {
    bool theFlag = self->flag();
    return theFlag;
  }
  void set_flag(bool theFlag)
  {
    self->flag() = theFlag;
  }
} // extend Flag_bool
%extend Flag {
  int get_flag()
  {
    int theFlag = self->flag();
    return theFlag;
  }
  void set_flag(int theFlag)
  {
    self->flag() = theFlag;
  }
} // extend Flag
} // namespace data32
namespace datatypes
{
%extend Flag_bool {
  bool get_flag()
  {
    bool theFlag = self->flag();
    return theFlag;
  }
  void set_flag(bool theFlag)
  {
    self->flag() = theFlag;
  }
  clipper::datatypes::Flag_bool copy()
  {
    clipper::datatypes::Flag_bool ret;
    ret = *self;
    return ret;
  }
} // extend Flag_bool
%extend Flag {
  int get_flag()
  {
    int theFlag = self->flag();
    return theFlag;
  }
  void set_flag(int theFlag)
  {
    self->flag() = theFlag;
  }
  clipper::datatypes::Flag copy()
  {
    clipper::datatypes::Flag ret;
    ret = *self;
    return ret;
  }
} // extend Flag
} // namespace datatypes
} // namespace clipper


//%rename (to_complex_float) operator complex<float>();
//%rename (to_complex_double) operator complex<double>();




namespace clipper
{
  %extend HKL_data {
    HKL_data<clipper::datatypes::Flag_bool> not_()
    {
      return !(*self);
    }
    HKL_data<clipper::datatypes::Flag_bool> __or__(const HKL_data<T> &d1)
    {
      return (*self) | d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __xor__(const HKL_data<T> &d1)
    {
      return (*self) ^ d1;
    }
    HKL_data<clipper::datatypes::Flag_bool> __and__(const HKL_data<T> &d1)
    {
      return (*self) & d1;
    }

  } // extend HKL_data
  %extend datatypes::ABCD {
    void vals(double numpy_double_out[4]) {
      numpy_double_out[0] = self->a();
      numpy_double_out[1] = self->b();
      numpy_double_out[2] = self->c();
      numpy_double_out[3] = self->d();
    }
  } // extend datatypes::ABCD





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
    ret = *self;
    return ret;
  }
} // extend HKL_data<clipper::data32::Flag_bool>

%extend HKL_data<clipper::data32::Flag> {
  HKL_data<clipper::datatypes::Flag_bool> __eq__(const int& n)
  {
    return (*self) == n;
  }
  HKL_data<clipper::datatypes::Flag_bool> __ne__(const int& n)
  {
    return (*self) != n;
  }
  HKL_data<clipper::datatypes::Flag_bool> __ge__(const int& n)
  {
    return (*self) >= n;
  }
  HKL_data<clipper::datatypes::Flag_bool> __le__(const int& n)
  {
    return (*self) <= n;
  }
  HKL_data<clipper::datatypes::Flag_bool> __gt__(const int& n)
  {
    return (*self) > n;
  }
  HKL_data<clipper::datatypes::Flag_bool> __lt__(const int& n)
  {
    return (*self) < n;
  }
  HKL_data<clipper::datatypes::Flag>  copy()
  {
    HKL_data<clipper::data32::Flag> ret;
    ret = *self;
    return ret;
  }
} // extend HKL_data<clipper::data32::Flag>

%extend datatypes::F_sigF<float> {
  clipper::datatypes::F_sigF<float>  copy()
  {
    clipper::data32::F_sigF ret;
    ret = *self;
    return ret;
  }
} // extend datatypes::F_sigF<float>

%extend datatypes::F_sigF_ano<float> {
  clipper::datatypes::F_sigF_ano<float>  copy()
  {
    clipper::data32::F_sigF_ano ret;
    ret = *self;
    return ret;
  }
} // extend datatypes::F_sigF_ano<float>

%extend datatypes::I_sigI<float> {
  clipper::datatypes::I_sigI<float>  copy()
  {
    clipper::data32::I_sigI ret;
    ret = *self;
    return ret;
  }
} // extend datatypes::I_sigI<float>

%extend datatypes::E_sigE<float> {
  clipper::datatypes::E_sigE<float>  copy()
  {
    clipper::data32::E_sigE ret;
    ret = *self;
    return ret;
  }
} // extend datatypes::E_sigE<float>

%extend datatypes::F_phi<float> {
  clipper::datatypes::F_phi<float>  __add__(const clipper::datatypes::F_phi<float> &h2)
  {
    clipper::data32::F_phi ret;
    ret = *self+h2;
    return ret;
  }
  clipper::datatypes::F_phi<float>  __sub__(const clipper::datatypes::F_phi<float> &h2)
  {
    clipper::data32::F_phi ret;
    ret = *self-h2;
    return ret;
  }
  clipper::datatypes::F_phi<float>  __neg__()
  {
    clipper::data32::F_phi ret;
    ret = -*self;
    return ret;
  }
  clipper::datatypes::F_phi<float>  copy()
  {
    clipper::data32::F_phi ret;
    ret = *self;
    return ret;
  }
} // extend datatypes::F_phi<float>

%extend datatypes::ABCD<float> {
  clipper::datatypes::ABCD<float>  __add__(const clipper::datatypes::ABCD<float> &h2)
  {
    clipper::data32::ABCD ret;
    ret = *self+h2;
    return ret;
  }
  clipper::datatypes::ABCD<float>  copy()
  {
    clipper::data32::ABCD ret;
    ret = *self;
    return ret;
  }
} // extend datatypes::ABCD<float>

%extend HKL_data<clipper::data32::ABCD> {
  HKL_data<clipper::datatypes::ABCD<float> > __add__(const HKL_data<clipper::datatypes::ABCD<float> > &h2)
  {
    HKL_data<clipper::data32::ABCD> ret;
    ret = *self+h2;
    return ret;
  }
  HKL_data<clipper::datatypes::ABCD<float> > copy()
  {
    HKL_data<clipper::data32::ABCD> ret;
    ret = *self;
    return ret;
  }
} // extend HKL_data<clipper::data32::ABCD>

%extend HKL_data<clipper::data32::F_phi> {
  HKL_data<clipper::datatypes::F_phi<float> >  copy()
  {
    HKL_data<clipper::data32::F_phi> ret;
    ret = *self;
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
    ret = *self;
    return ret;
  }
  */
  HKL_data<clipper::datatypes::F_phi<float> > __add__(const HKL_data<clipper::datatypes::F_phi<float> > &h2)
  {
    HKL_data<clipper::data32::F_phi> ret;
    ret = *self+h2;
    return ret;
  }
  HKL_data<clipper::datatypes::F_phi<float> > __sub__(const HKL_data<clipper::datatypes::F_phi<float> > &h2)
  {
    HKL_data<clipper::datatypes::F_phi<float> > ret;
    ret = *self-h2;
    return ret;
  }
  HKL_data<clipper::datatypes::F_phi<float> > __neg__()
  {
    HKL_data<clipper::datatypes::F_phi<float> > ret;
    ret = -*self;
    return ret;
  }
  HKL_data<clipper::datatypes::F_phi<float> > __mul__(const float s)
  {
    HKL_data<clipper::datatypes::F_phi<float> > ret;
    ret = *self*s;
    return ret;
  }
  HKL_data<clipper::datatypes::F_phi<float> > __rmul__(const float s)
  {
    HKL_data<clipper::datatypes::F_phi<float> > ret;
    ret = *self*s;
    return ret;
  }
} // extend HKL_data<clipper::data32::F_phi>

%extend HKL_data<clipper::data32::E_sigE> {
  void scaleBySqrtResolution(const clipper::ResolutionFn &escale)
  {
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() )
      if ( !(*self)[ih].missing() ) (*self)[ih].scale( sqrt( escale.f(ih) ) );
  }
  void scaleByResolution(const clipper::ResolutionFn &escale)
  {
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() )
      if ( !(*self)[ih].missing() ) (*self)[ih].scale( escale.f(ih) );
  }
  HKL_data<clipper::datatypes::E_sigE<float> >  copy()
  {
    HKL_data<clipper::data32::E_sigE> ret;
    ret = *self;
    return ret;
  }
} // extend HKL_data<clipper::data32::E_sigE>

%extend HKL_data<clipper::data32::ABCD> {
  void compute_from_phi_fom(const HKL_data< clipper::datatypes::Phi_fom<float> > &phiw)
  {
    self->compute( phiw, clipper::data32::Compute_abcd_from_phifom() );
  }
  void compute_add_abcd(const HKL_data< clipper::datatypes::ABCD<float> > &abcd1,
  const HKL_data< clipper::datatypes::ABCD<float> > &abcd2)
  {
    self->compute( abcd1, abcd2, clipper::data32::Compute_add_abcd() );
  }
} // extend HKL_data<clipper::data32::ABCD>

%extend HKL_data<clipper::data32::Phi_fom> {
  void compute_from_abcd(const HKL_data< clipper::datatypes::ABCD<float> > &abcd)
  {
    self->compute( abcd, clipper::data32::Compute_phifom_from_abcd() );
  }
  HKL_data<clipper::datatypes::Phi_fom<float> >  copy()
  {
    HKL_data<clipper::data32::Phi_fom> ret;
    ret = *self;
    return ret;
  }
} // extend HKL_data<clipper::data32::Phi_fom>

%extend HKL_data<clipper::data32::F_sigF> {
  void compute_mean_from_fano(const HKL_data< clipper::datatypes::F_sigF_ano<float> > &fano)
  {
    self->compute( fano, clipper::data32::Compute_mean_fsigf_from_fsigfano() );
  }
  void compute_diff_from_fano(const HKL_data< clipper::datatypes::F_sigF_ano<float> > &fano)
  {
    self->compute( fano, clipper::data32::Compute_diff_fsigf_from_fsigfano() );
  }
  void compute_scale_u_iso_fsigf(float scale, float u_value,
  const HKL_data< clipper::datatypes::F_sigF<float> > &fsigf )
  {
    self->compute( fsigf, clipper::data32::Compute_scale_u_iso_fsigf(scale, u_value) );
  }
  void compute_scale_u_aniso_fsigf(float scale, clipper::U_aniso_orth u_value,
  const HKL_data< clipper::datatypes::F_sigF<float> > &fsigf )
  {
    self->compute( fsigf, clipper::data32::Compute_scale_u_aniso_fsigf(scale, u_value) );
  }
  HKL_data<clipper::datatypes::F_sigF<float> >  copy()
  {
    HKL_data<clipper::data32::F_sigF> ret;
    ret = *self;
    return ret;
  }
} // extend HKL_data<clipper::data32::F_sigF>

%extend HKL_data<clipper::data32::F_sigF_ano> {
  void compute_scale_u_iso_fsigfano(float scale, float u_value,
  const HKL_data< clipper::datatypes::F_sigF_ano<float> > &fsigfano )
  {
    self->compute( fsigfano, clipper::data32::Compute_scale_u_iso_fsigfano(scale, u_value) );
  }
  void compute_scale_u_aniso_fsigfano(float scale, clipper::U_aniso_orth u_value,
  const HKL_data< clipper::datatypes::F_sigF_ano<float> > &fsigfano )
  {
    self->compute( fsigfano, clipper::data32::Compute_scale_u_aniso_fsigfano(scale, u_value) );
  }
  HKL_data<clipper::datatypes::F_sigF_ano<float> >  copy()
  {
    HKL_data<clipper::data32::F_sigF_ano> ret;
    ret = *self;
    return ret;
  }
} // extend HKL_data<clipper::data32::F_sigF_ano>

%extend HKL_data<clipper::data32::I_sigI> {
  void compute_scale_u_iso_isigi(float scale, float u_value,
  const HKL_data< clipper::datatypes::I_sigI<float> > &isigi )
  {
    self->compute( isigi, clipper::data32::Compute_scale_u_iso_isigi(scale, u_value) );
  }
  void compute_scale_u_aniso_isigi(float scale, clipper::U_aniso_orth u_value,
  const HKL_data< clipper::datatypes::I_sigI<float> > &isigi )
  {
    self->compute( isigi, clipper::data32::Compute_scale_u_aniso_isigi(scale, u_value) );
  }
  HKL_data<clipper::datatypes::I_sigI<float> >  copy()
  {
    HKL_data<clipper::data32::I_sigI> ret;
    ret = *self;
    return ret;
  }
} // extend HKL_data<clipper::data32::I_sigI>



%extend HKL_data< clipper::datatypes::F_phi<float> > {

  void getDataNumpy(double *test_numpy_a, int test_numpy_n)
  {
    int i=0;
    for(clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next(),i++ ) {
      if(!((*self)[ih].missing())) {
        std::vector<xtype> thisData(self->data_size());
        self->data_export(ih.hkl(),&(thisData[0]));
        std::vector<float> thisDataf(self->data_size());
        for(unsigned idat=0; idat<self->data_size(); ++idat) {
          test_numpy_a[i*self->data_size()+idat] = thisData[idat];
        }
      } else {
        for(unsigned idat=0; idat<self->data_size(); ++idat) {
          test_numpy_a[i*self->data_size()+idat] = std::numeric_limits<float>::quiet_NaN();
        }
      }
    }
  }

  std::vector<std::vector<float> > getData()
  {
    std::vector<std::vector<float> > allData;
    for(clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      if(!((*self)[ih].missing())) {
        std::vector<xtype> thisData(self->data_size());
        self->data_export(ih.hkl(),&(thisData[0]));
        std::vector<float> thisDataf(self->data_size());
        for(unsigned idat=0; idat<self->data_size(); ++idat) {
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
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*self)[i];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    return sz;
  }
} // extend HKL_data< clipper::datatypes::F_phi<float> >

%extend HKL_data<clipper::data32::F_phi> {
  void compute_neg(const HKL_data< clipper::datatypes::F_phi<float> > &fphi )
  {
    self->compute( fphi, clipper::data32::Compute_neg_fphi() );
  }
  void compute_add_fphi(const HKL_data< clipper::datatypes::F_phi<float> > &fphi1,
  const HKL_data< clipper::datatypes::F_phi<float> > &fphi2)
  {
    self->compute( fphi1, fphi2, clipper::data32::Compute_add_fphi() );
  }
  void compute_sub_fphi(const HKL_data< clipper::datatypes::F_phi<float> > &fphi1,
  const HKL_data< clipper::datatypes::F_phi<float> > &fphi2)
  {
    self->compute( fphi1, fphi2, clipper::data32::Compute_sub_fphi() );
  }
  void compute_from_fsigf_phifom(const HKL_data< clipper::datatypes::F_sigF<float> > &fsigf,
  const HKL_data< clipper::datatypes::Phi_fom<float> > &phifom )
  {
    self->compute( fsigf, phifom, clipper::data32::Compute_fphi_from_fsigf_phifom() );
  }
  void compute_scale_u_iso_fphi(float scale, float u_value,
  const HKL_data< clipper::datatypes::F_phi<float> > &fphi )
  {
    self->compute( fphi, clipper::data32::Compute_scale_u_iso_fphi(scale, u_value) );
  }
  void compute_scale_u_aniso_fphi(float scale, clipper::U_aniso_orth u_value,
  const HKL_data< clipper::datatypes::F_phi<float> > &fphi )
  {
    self->compute( fphi, clipper::data32::Compute_scale_u_aniso_fphi(scale, u_value) );
  }

  void getDataNumpy(double *test_numpy_a, int test_numpy_n)
  {
    int i=0;
    for(clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next(),i++ ) {
      if(!((*self)[ih].missing())) {
        std::vector<xtype> thisData(self->data_size());
        self->data_export(ih.hkl(),&(thisData[0]));
        std::vector<float> thisDataf(self->data_size());
        for(unsigned idat=0; idat<self->data_size(); ++idat) {
          test_numpy_a[i*self->data_size()+idat] = thisData[idat];
        }
      } else {
        for(unsigned idat=0; idat<self->data_size(); ++idat) {
          test_numpy_a[i*self->data_size()+idat] = std::numeric_limits<float>::quiet_NaN();
        }
      }
    }
  }

  std::vector<std::vector<float> > getData()
  {
    std::vector<std::vector<float> > allData;
    for(clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      if(!((*self)[ih].missing())) {
        std::vector<xtype> thisData(self->data_size());
        self->data_export(ih.hkl(),&(thisData[0]));
        std::vector<float> thisDataf(self->data_size());
        for(unsigned idat=0; idat<self->data_size(); ++idat) {
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
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*self)[i];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    return sz;
  }
} // extend HKL_data<clipper::data32::F_phi>

%extend HKL_data<clipper::data32::E_sigE> {
  void compute_from_fsigf(const HKL_data< clipper::datatypes::F_sigF<float> > &fsigf )
  {
    self->compute( fsigf, clipper::data32::Compute_EsigE_from_FsigF() );
  }
} // extend HKL_data<clipper::data32::E_sigE>

%extend HKL_data<clipper::datatypes::Flag_bool> {
  clipper::datatypes::Flag_bool& __getitem__(int i)
  {
    int sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*self)[i];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    return sz;
  }
} // extend HKL_data<clipper::datatypes::Flag_bool>

%extend HKL_data<clipper::data32::Flag_bool> {
  clipper::data32::Flag_bool& __getitem__(int i)
  {
    int sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*self)[i];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    return sz;
  }
} // extend HKL_data<clipper::data32::Flag_bool>

%extend HKL_data<clipper::datatypes::Flag> {
  clipper::datatypes::Flag& __getitem__(int i)
  {
    int sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*self)[i];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    return sz;
  }
} // extend HKL_data<clipper::datatypes::Flag>

%extend HKL_data<clipper::data32::Flag> {
  clipper::data32::Flag& __getitem__(int i)
  {
    int sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*self)[i];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    return sz;
  }
} // extend HKL_data<clipper::data32::Flag>


%extend HKL_data< clipper::datatypes::F_sigF<float> > {

  void getDataNumpy(double *test_numpy_a, int test_numpy_n)
  {
    int i=0;
    for(clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next(),i++ ) {
      if(!((*self)[ih].missing())) {
        std::vector<xtype> thisData(self->data_size());
        self->data_export(ih.hkl(),&(thisData[0]));
        std::vector<float> thisDataf(self->data_size());
        for(unsigned idat=0; idat<self->data_size(); ++idat) {
          test_numpy_a[i*self->data_size()+idat] = thisData[idat];
        }
      } else {
        for(unsigned idat=0; idat<self->data_size(); ++idat) {
          test_numpy_a[i*self->data_size()+idat] = std::numeric_limits<float>::quiet_NaN();
        }
      }
    }
  }

  std::vector<std::vector<float> > getData()
  {
    std::vector<std::vector<float> > allData;
    for(clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      if(!((*self)[ih].missing())) {
        std::vector<xtype> thisData(self->data_size());
        self->data_export(ih.hkl(),&(thisData[0]));
        std::vector<float> thisDataf(self->data_size());
        for(unsigned idat=0; idat<self->data_size(); ++idat) {
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
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*self)[i];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    return sz;
  }
} // extend HKL_data< clipper::datatypes::F_sigF<float> >

%extend HKL_data< clipper::data32::F_sigF<float> > {
  clipper::data32::F_sigF<float>& __getitem__(int i)
  {
    int sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    i = (i < 0) ? sz + i : i;
    if (i >= sz || i < 0) {
      throw std::out_of_range("");
    }
    return (*self)[i];
fail:
    return (*self)[0];
  }
  size_t __len__()
  {
    size_t sz=0;
    for ( clipper::HKL_data_base::HKL_reference_index ih = self->first(); !ih.last(); ih.next() ) {
      sz++;
    }
    return sz;
  }
} // extend HKL_data< clipper::data32::F_sigF<float> >

} // namespace clipper
