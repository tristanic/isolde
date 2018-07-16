// clipper CNS->MTZ utility
/* (C) 2007 Kevin Cowtan */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-mmdb.h>

using clipper::ftype;
using clipper::xtype;
using clipper::String;
using clipper::Datatype_base;
using clipper::Util;

//! Reflection data type: I
class dataI : private Datatype_base
{
public:
  dataI() { Util::set_null(i_); }
  dataI( const float& i ) : i_(i) {}
  void set_null() { Util::set_null(i_); }
  static String type() { return ""; }
  void friedel() {}
  void shift_phase(const ftype& dphi) {}
  bool missing() const { return (Util::is_nan(i_)); }
  static int data_size() { return 1; }
  static String data_names() { return "I"; }
  void data_export( xtype array[] ) const { array[0] = I(); }
  void data_import( const xtype array[] ) { I() = array[0]; }
  void scale(const ftype& s) { i_ *= s*s; }
  // accessors
  const float& I() const { return i_; }  //<! read access
  float& I() { return i_; }  //<! write access
private:
  float i_;
};

//! Reflection data type: F
class dataF : private Datatype_base
{
public:
  dataF() { Util::set_null(f_); }
  dataF( const float& f ) : f_(f) {}
  void set_null() { Util::set_null(f_); }
  static String type() { return ""; }
  void friedel() {}
  void shift_phase(const ftype& dphi) {}
  bool missing() const { return (Util::is_nan(f_)); }
  static int data_size() { return 1; }
  static String data_names() { return "F"; }
  void data_export( xtype array[] ) const { array[0] = f(); }
  void data_import( const xtype array[] ) { f() = array[0]; }
  void scale(const ftype& s) { f_ *= s; }
  // accessors
  const float& f() const { return f_; }  //<! read access
  float& f() { return f_; }  //<! write access
private:
  float f_;
};

//! Reflection data type: E
class dataE : private Datatype_base
{
public:
  dataE() { Util::set_null(e_); }
  dataE( const float& e ) : e_(e) {}
  void set_null() { Util::set_null(e_); }
  static String type() { return ""; }
  void friedel() {}
  void shift_phase(const ftype& dphi) {}
  bool missing() const { return (Util::is_nan(e_)); }
  static int data_size() { return 1; }
  static String data_names() { return "E"; }
  void data_export( xtype array[] ) const { array[0] = E(); }
  void data_import( const xtype array[] ) { E() = array[0]; }
  void scale(const ftype& s) { e_ *= s; }
  // accessors
  const float& E() const { return e_; }  //<! read access
  float& E() { return e_; }  //<! write access
private:
  float e_;
};

//! Reflection data type: sig
class dataSig : private Datatype_base
{
public:
  dataSig() { Util::set_null(sig_); }
  dataSig( const float& sig ) : sig_(sig) {}
  void set_null() { Util::set_null(sig_); }
  static String type() { return ""; }
  void friedel() {}
  void shift_phase(const ftype& dphi) {}
  bool missing() const { return (Util::is_nan(sig_)); }
  static int data_size() { return 1; }
  static String data_names() { return "sigF"; }
  void data_export( xtype array[] ) const { array[0] = sig(); }
  void data_import( const xtype array[] ) { sig() = array[0]; }
  void scale(const ftype& s) { sig_ *= s; }
  // accessors
  const float& sig() const { return sig_; }  //<! read access
  float& sig() { return sig_; }  //<! write access
private:
  float sig_;
};

//! Reflection data type: phi
class dataPhi : private Datatype_base
{
public:
  dataPhi() { Util::set_null(phi_); }
  dataPhi( const float& phi ) : phi_(phi) {}
  void set_null() { Util::set_null(phi_); }
  static String type() { return ""; }
  void friedel()
    { if (!Util::is_nan(phi_)) phi_=-phi_; }
  void shift_phase(const ftype& dphi)
    { if (!Util::is_nan(phi_)) phi_+=dphi; }
  bool missing() const { return (Util::is_nan(phi_)); }
  static int data_size() { return 1; }
  static String data_names() { return "phi"; }
  void data_export( xtype array[] ) const { array[0] = phi(); }
  void data_import( const xtype array[] ) { phi() = array[0]; }
  void scale(const ftype& s) { phi_ *= s; }
  // accessors
  const float& phi() const { return phi_; }  //<! read access
  float& phi() { return phi_; }  //<! write access
private:
  float phi_;
};

//! Reflection data type: fom
class dataFom : private Datatype_base
{
public:
  dataFom() { Util::set_null(fom_); }
  dataFom( const float& fom ) : fom_(fom) {}
  void set_null() { Util::set_null(fom_); }
  static String type() { return ""; }
  void friedel() {}
  void shift_phase(const ftype& dfom) {}
  bool missing() const { return (Util::is_nan(fom_)); }
  static int data_size() { return 1; }
  static String data_names() { return "fom"; }
  void data_export( xtype array[] ) const { array[0] = fom(); }
  void data_import( const xtype array[] ) { fom() = array[0]; }
  void scale(const ftype& s) { fom_ *= s; }
  // accessors
  const float& fom() const { return fom_; }  //<! read access
  float& fom() { return fom_; }  //<! write access
private:
  float fom_;
};

//! Reflection data type: iano
class dataIano : private Datatype_base
{
public:
  dataIano() { Util::set_null(ipl_); Util::set_null(imi_); }
  dataIano( const float& ipl, const float& imi ) : ipl_(ipl), imi_(imi) {}
  void set_null() { Util::set_null(ipl_); Util::set_null(imi_); }
  static String type() { return ""; }
  void friedel() { float i=ipl_; ipl_=imi_; imi_=i; }
  void shift_phase(const ftype& dphi) {}
  bool missing() const { return (Util::is_nan(ipl_) && Util::is_nan(imi_)); }
  static int data_size() { return 2; }
  static String data_names() { return "I+ I-"; }
  void data_export( xtype array[] ) const
    { array[0] = I_pl(); array[1] = I_mi(); }
  void data_import( const xtype array[] )
    { I_pl() = array[0]; I_mi() = array[1]; }
  void scale(const ftype& s) { ipl_ *= s; imi_ *= s; }
  // accessors
  const float& I_pl() const { return ipl_; }  //<! read access
  const float& I_mi() const { return imi_; }  //<! read access
  float& I_pl() { return ipl_; }  //<! write access
  float& I_mi() { return imi_; }  //<! write access
private:
  float ipl_,imi_;
};

//! Reflection data type: fano
class dataSigIano : private Datatype_base
{
public:
  dataSigIano() { Util::set_null(sigipl_); Util::set_null(sigimi_); }
  dataSigIano( const float& sigipl, const float& sigimi ) : sigipl_(sigipl), sigimi_(sigimi) {}
  void set_null() { Util::set_null(sigipl_); Util::set_null(sigimi_); }
  static String type() { return ""; }
  void friedel() { float sigi=sigipl_; sigipl_=sigimi_; sigimi_=sigi; }
  void shift_phase(const ftype& dphi) {}
  bool missing() const { return (Util::is_nan(sigipl_) && Util::is_nan(sigimi_)); }
  static int data_size() { return 2; }
  static String data_names() { return "sigI+ sigI-"; }
  void data_export( xtype array[] ) const
    { array[0] = sigI_pl(); array[1] = sigI_mi(); }
  void data_import( const xtype array[] )
    { sigI_pl() = array[0]; sigI_mi() = array[1]; }
  void scale(const ftype& s) { sigipl_ *= s; sigimi_ *= s; }
  // accessors
  const float& sigI_pl() const { return sigipl_; }  //<! read access
  const float& sigI_mi() const { return sigimi_; }  //<! read access
  float& sigI_pl() { return sigipl_; }  //<! write access
  float& sigI_mi() { return sigimi_; }  //<! write access
private:
  float sigipl_,sigimi_;
};

//! Reflection data type: fano
class dataFano : private Datatype_base
{
public:
  dataFano() { Util::set_null(fpl_); Util::set_null(fmi_); }
  dataFano( const float& fpl, const float& fmi ) : fpl_(fpl), fmi_(fmi) {}
  void set_null() { Util::set_null(fpl_); Util::set_null(fmi_); }
  static String type() { return ""; }
  void friedel() { float f=fpl_; fpl_=fmi_; fmi_=f; }
  void shift_phase(const ftype& dphi) {}
  bool missing() const { return (Util::is_nan(fpl_) && Util::is_nan(fmi_)); }
  static int data_size() { return 2; }
  static String data_names() { return "F+ F-"; }
  void data_export( xtype array[] ) const
    { array[0] = F_pl(); array[1] = F_mi(); }
  void data_import( const xtype array[] )
    { F_pl() = array[0]; F_mi() = array[1]; }
  void scale(const ftype& s) { fpl_ *= s; fmi_ *= s; }
  // accessors
  const float& F_pl() const { return fpl_; }  //<! read access
  const float& F_mi() const { return fmi_; }  //<! read access
  float& F_pl() { return fpl_; }  //<! write access
  float& F_mi() { return fmi_; }  //<! write access
private:
  float fpl_,fmi_;
};

//! Reflection data type: fano
class dataSigFano : private Datatype_base
{
public:
  dataSigFano() { Util::set_null(sigfpl_); Util::set_null(sigfmi_); }
  dataSigFano( const float& sigfpl, const float& sigfmi ) : sigfpl_(sigfpl), sigfmi_(sigfmi) {}
  void set_null() { Util::set_null(sigfpl_); Util::set_null(sigfmi_); }
  static String type() { return ""; }
  void friedel() { float sigf=sigfpl_; sigfpl_=sigfmi_; sigfmi_=sigf; }
  void shift_phase(const ftype& dphi) {}
  bool missing() const { return (Util::is_nan(sigfpl_) && Util::is_nan(sigfmi_)); }
  static int data_size() { return 2; }
  static String data_names() { return "sigF+ sigF-"; }
  void data_export( xtype array[] ) const
    { array[0] = sigF_pl(); array[1] = sigF_mi(); }
  void data_import( const xtype array[] )
    { sigF_pl() = array[0]; sigF_mi() = array[1]; }
  void scale(const ftype& s) { sigfpl_ *= s; sigfmi_ *= s; }
  // accessors
  const float& sigF_pl() const { return sigfpl_; }  //<! read access
  const float& sigF_mi() const { return sigfmi_; }  //<! read access
  float& sigF_pl() { return sigfpl_; }  //<! write access
  float& sigF_mi() { return sigfmi_; }  //<! write access
private:
  float sigfpl_,sigfmi_;
};

typedef clipper::data32::ABCD dataABCD;
typedef clipper::data32::Flag dataFlag;


class HKLlessthan
{
 public:
  bool operator()( const clipper::HKL& h1, const clipper::HKL& h2 ) const
    { return ( h1.h() < h2.h() ||
	       ( h1.h() == h2.h() && ( h1.k() < h2.k() ||
				       ( h1.k() == h2.k() && h1.l() < h2.l() )
				       ) ) ); }
};


/*
//! Reflection data type: sigI
class dataSigI : private Datatype_base
{
public:
  dataSigI() { Util::set_null(sigi_); }
  dataSigI( const float& sigi ) : sigi_(sigi) {}
  void set_null() { Util::set_null(sigi_); }
  static String type() { return ""; }
  void friedel() {}
  void shift_phase(const ftype& dphi) {}
  bool missing() const { return (Util::is_nan(sigi_)); }
  static int data_size() { return 1; }
  static String data_names() { return "I"; }
  void data_export( xtype array[] ) const { array[0] = sigI(); }
  void data_import( const xtype array[] ) { sigI() = array[0]; }
  void scale(const ftype& s) { sigi_ *= s*s; }
  // accessors
  const float& sigI() const { return sigi_; }  //<! read access
  float& sigI() { return sigi_; }  //<! write access
private:
  float sigi_;
};
*/
/*
//! Reflection data type: sigF
class dataSigF : private Datatype_base
{
public:
  dataSigF() { Util::set_null(sigf_); }
  dataSigF( const float& sigf ) : sigf_(sigf) {}
  void set_null() { Util::set_null(sigf_); }
  static String type() { return ""; }
  void friedel() {}
  void shift_phase(const ftype& dphi) {}
  bool missing() const { return (Util::is_nan(sigf_)); }
  static int data_size() { return 1; }
  static String data_names() { return "I"; }
  void data_export( xtype array[] ) const { array[0] = sigF(); }
  void data_import( const xtype array[] ) { sigF() = array[0]; }
  void scale(const ftype& s) { sigf_ *= s*s; }
  // accessors
  const float& sigF() const { return sigf_; }  //<! read access
  float& sigF() { return sigf_; }  //<! write access
private:
  float sigf_;
};
*/
