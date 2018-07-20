/*! \file lib/hkl_datatypes.h
    Fundamental data types for the clipper libraries
*/
//C Copyright (C) 2000-2006 Kevin Cowtan and University of York
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA


#ifndef CLIPPER_HKL_DATATYPES
#define CLIPPER_HKL_DATATYPES

#include <complex>
#include "hkl_data.h"
#include "../imex.h"

namespace clipper
{

  // Now define some actual datatypes


  namespace datatypes
  {

    //! Reflection data type: I + sigI
    /*! Note that I_sigI also has methods for returning I_pl(),
      sigI_pl(), I_mi, sigI_mi(), so you can use this type in any
      template type where you would use I_sigI_ano. */
    template<class dtype> class I_sigI : private Datatype_base
    {
    public:
      I_sigI() { Util::set_null(I_); Util::set_null(sigI_); }
      I_sigI( const dtype& I, const dtype& sigI ) : I_(I), sigI_(sigI) {}
      void set_null() { Util::set_null(I_); Util::set_null(sigI_); }
      static String type() { return "I_sigI"; }
      void friedel() {}
      void shift_phase(const ftype&) {}
      bool missing() const { return (Util::is_nan(I_) || Util::is_nan(sigI_)); }
      static int data_size() { return 2; }
      static String data_names() { return "I sigI"; }
      void data_export( xtype array[] ) const
	{ array[0] = I(); array[1] = sigI(); }
      void data_import( const xtype array[] )
	{ I() = array[0]; sigI() = array[1]; }
      //! this type is scalable - apply magnitude scale factor
      void scale(const ftype& s) { I_ *= (s*s); sigI_ *= (s*s); }
      // accessors
      const dtype& I() const { return I_; }  //<! read access
      const dtype& sigI() const { return sigI_; }  //<! read access
      dtype& I() { return I_; }  //<! write access
      dtype& sigI() { return sigI_; }  //<! write access
      // anomalous-type accessors
      const dtype& I_pl() const { return I_; }  //<! read access as anom
      const dtype& sigI_pl() const { return sigI_; }  //<! read access as anom
      const dtype& I_mi() const { return I_; }  //<! read access as anom
      const dtype& sigI_mi() const { return sigI_; }  //<! read access as anom
      dtype cov() const { return 1.0; }  //<! read access as anom
    private:
      dtype I_,sigI_;
    };

    //! Reflection data type: I(+) I(+) sigI(+) sigI(-) cov+-
    /*! Note that I_sigI_ano also has methods for returning I(),
      sigI(), so you can use this type in any template type where you
      would use I_sigI. */
    template<class dtype> class I_sigI_ano : private Datatype_base
    { public:
      I_sigI_ano() { set_null(); }
      void set_null() { Util::set_null(I_pl_); Util::set_null(I_mi_); Util::set_null(sigI_pl_); Util::set_null(sigI_mi_); Util::set_null(cov_); }
      static String type() { return "I_sigI_ano"; }
      void friedel() { dtype I=I_pl_; I_pl_=I_mi_; I_mi_=I;
                       I=sigI_pl_; sigI_pl_=sigI_mi_; sigI_mi_=I; }
      void shift_phase(const ftype&) {}
      bool missing() const { return (Util::is_nan(I_pl_) && Util::is_nan(I_mi_)); }
      static int data_size() { return 5; }
      static String data_names() { return "I+ sigI+ I- sigI- covI+-"; }
      void data_export( xtype a[] ) const { a[0] = I_pl(); a[1] = sigI_pl(); a[2] = I_mi(); a[3] = sigI_mi(); a[4] = cov(); }
      void data_import( const xtype a[] ) { I_pl() = a[0]; sigI_pl() = a[1]; I_mi() = a[2]; sigI_mi() = a[3]; cov() = a[4]; }
      //! this type is scalable - apply magnitude scale factor
      void scale(const ftype& s) { I_pl_ *= (s*s); sigI_pl_ *= (s*s); I_mi_ *= (s*s); sigI_mi_ *= (s*s); cov_ *= (s*s); }
      // accessors
      const dtype& I_pl() const { return I_pl_; }  //<! read access
      const dtype& sigI_pl() const { return sigI_pl_; }  //<! read access
      const dtype& I_mi() const { return I_mi_; }  //<! read access
      const dtype& sigI_mi() const { return sigI_mi_; }  //<! read access
      const dtype& cov() const { return cov_; }  //<! read access
      dtype& I_pl() { return I_pl_; }  //<! write access
      dtype& sigI_pl() { return sigI_pl_; }  //<! write access
      dtype& I_mi() { return I_mi_; }  //<! write access
      dtype& sigI_mi() { return sigI_mi_; }  //<! write access
      dtype& cov() { return cov_; }  //<! write access
      // nonanomalous-type accessors
      dtype I() const { return Util::mean(I_pl_,I_mi_); }  //<! read access as simple
      dtype sigI() const { return Util::sig_mean(sigI_pl_,sigI_mi_,cov_); }  //<! read access as simple
    private:
      dtype I_pl_, I_mi_, sigI_pl_, sigI_mi_, cov_;
    };

    //! Reflection data type: F + sigF
    /*! Note that F_sigF also has methods for returning f_pl(),
      sigf_pl(), f_mi, sigf_mi(), so you can use this type in any
      template type where you would use F_sigF_ano. */
    template<class dtype> class F_sigF : private Datatype_base
    {
    public:
      F_sigF() { Util::set_null(f_); Util::set_null(sigf_); }
      F_sigF( const dtype& f, const dtype& sigf ) : f_(f), sigf_(sigf) {}
      void set_null() { Util::set_null(f_); Util::set_null(sigf_); }
      static String type() { return "F_sigF"; }
      void friedel() {}
      void shift_phase(const ftype&) {}
      bool missing() const { return (Util::is_nan(f_) || Util::is_nan(sigf_)); }
      static int data_size() { return 2; }
      static String data_names() { return "F sigF"; }
      void data_export( xtype array[] ) const
	{ array[0] = f(); array[1] = sigf(); }
      void data_import( const xtype array[] )
	{ f() = array[0]; sigf() = array[1]; }
      //! this type is scalable - apply magnitude scale factor
      void scale(const ftype& s) { f_ *= s; sigf_ *= s; }
      // accessors
      const dtype& f() const { return f_; }  //<! read access
      const dtype& sigf() const { return sigf_; }  //<! read access
      dtype& f() { return f_; }  //<! write access
      dtype& sigf() { return sigf_; }  //<! write access
      // anomalous-type accessors
      const dtype& f_pl() const { return f_; }  //<! read access as anom
      const dtype& sigf_pl() const { return sigf_; }  //<! read access as anom
      const dtype& f_mi() const { return f_; }  //<! read access as anom
      const dtype& sigf_mi() const { return sigf_; }  //<! read access as anom
      dtype cov() const { return 1.0; }  //<! read access as anom
    private:
      dtype f_,sigf_;
    };

    //! Reflection data type: F(+) F(+) sigF(+) sigF(-) cov+-
    /*! Note that F_sigF_ano also has methods for returning f(),
      sigf(), so you can use this type in any template type where you
      would use F_sigF. */
    template<class dtype> class F_sigF_ano : private Datatype_base
    { public:
      F_sigF_ano() { set_null(); }
      void set_null() { Util::set_null(f_pl_); Util::set_null(f_mi_); Util::set_null(sigf_pl_); Util::set_null(sigf_mi_); Util::set_null(cov_); }
      static String type() { return "F_sigF_ano"; }
      void friedel() { dtype f=f_pl_; f_pl_=f_mi_; f_mi_=f;
                       f=sigf_pl_; sigf_pl_=sigf_mi_; sigf_mi_=f; }
      void shift_phase(const ftype&) {}
      bool missing() const { return (Util::is_nan(f_pl_) && Util::is_nan(f_mi_)); }
      static int data_size() { return 5; }
      static String data_names() { return "F+ sigF+ F- sigF- covF+-"; }
      void data_export( xtype a[] ) const { a[0] = f_pl(); a[1] = sigf_pl(); a[2] = f_mi(); a[3] = sigf_mi(); a[4] = cov(); }
      void data_import( const xtype a[] ) { f_pl() = a[0]; sigf_pl() = a[1]; f_mi() = a[2]; sigf_mi() = a[3]; cov() = a[4]; }
      //! this type is scalable - apply magnitude scale factor
      void scale(const ftype& s) { f_pl_ *= s; sigf_pl_ *= s; f_mi_ *= s; sigf_mi_ *= s; cov_ *= (s*s); }
      // accessors
      const dtype& f_pl() const { return f_pl_; }  //<! read access
      const dtype& sigf_pl() const { return sigf_pl_; }  //<! read access
      const dtype& f_mi() const { return f_mi_; }  //<! read access
      const dtype& sigf_mi() const { return sigf_mi_; }  //<! read access
      const dtype& cov() const { return cov_; }  //<! read access
      dtype& f_pl() { return f_pl_; }  //<! write access
      dtype& sigf_pl() { return sigf_pl_; }  //<! write access
      dtype& f_mi() { return f_mi_; }  //<! write access
      dtype& sigf_mi() { return sigf_mi_; }  //<! write access
      dtype& cov() { return cov_; }  //<! write access
      // nonanomalous-type accessors
      dtype f() const { return Util::mean(f_pl_,f_mi_); }  //<! read access as simple
      dtype sigf() const { return Util::sig_mean(sigf_pl_,sigf_mi_,cov_); }  //<! read access as simple
    private:
      dtype f_pl_, f_mi_, sigf_pl_, sigf_mi_, cov_;
    };

    //! Reflection data type: E + sigE
    /*! This is not strictly a type for storing E values, but rather a
      type for storing any sturcture factor magnitude-like quantity
      which has already had a symmetry enhancement factor (epsilon)
      removed from it. E's are most commonly stored in this form,
      wheras F's and U's are not. You can compute corrected F's from
      uncorrected F's using:
      \code
      clipper::HKL_data<clipper::data32::F_sigF> fsigf;
      clipper::HKL_data<clipper::data32::E_sigE> esige;
      esige.compute( fsigf, clipper::data32::Compute_EsigE_from_FsigF() );
      \endcode
    */
    template<class dtype> class E_sigE : private Datatype_base
    {
    public:
      E_sigE() { Util::set_null(E_); Util::set_null(sigE_); }
      E_sigE( const dtype& E, const dtype& sigE ) : E_(E), sigE_(sigE) {}
      void set_null() { Util::set_null(E_); Util::set_null(sigE_); }
      static String type() { return "E_sigE"; }
      void friedel() {}
      void shift_phase(const ftype&) {}
      bool missing() const { return (Util::is_nan(E_) || Util::is_nan(sigE_)); }
      static int data_size() { return 2; }
      static String data_names() { return "E sigE"; }
      void data_export( xtype array[] ) const
	{ array[0] = E(); array[1] = sigE(); }
      void data_import( const xtype array[] )
	{ E() = array[0]; sigE() = array[1]; }
      //! this type is scalable - apply magnitude scale factor
      void scale(const ftype& s) { E_ *= s; sigE_ *= s; }
      // accessors
      const dtype& E() const { return E_; }  //<! read access
      const dtype& sigE() const { return sigE_; }  //<! read access
      dtype& E() { return E_; }  //<! write access
      dtype& sigE() { return sigE_; }  //<! write access
      // anomalous-type accessors
      const dtype& E_pl() const { return E_; }  //<! read access as anom
      const dtype& sigE_pl() const { return sigE_; }  //<! read access as anom
      const dtype& E_mi() const { return E_; }  //<! read access as anom
      const dtype& sigE_mi() const { return sigE_; }  //<! read access as anom
      dtype cov() const { return 1.0; }  //<! read access as anom
    private:
      dtype E_,sigE_;
    };

    //! Reflection data type: E(+) E(+) sigE(+) sigE(-) cov+-
    /*! see datatypes::E_sigE */
    template<class dtype> class E_sigE_ano : private Datatype_base
    { public:
      E_sigE_ano() { set_null(); }
      void set_null() { Util::set_null(E_pl_); Util::set_null(E_mi_); Util::set_null(sigE_pl_); Util::set_null(sigE_mi_); Util::set_null(cov_); }
      static String type() { return "E_sigE_ano"; }
      void friedel() { dtype e=E_pl_; E_pl_=E_mi_; E_mi_=e;
                       e=sigE_pl_; sigE_pl_=sigE_mi_; sigE_mi_=e; }
      void shift_phase(const ftype&) {}
      bool missing() const { return (Util::is_nan(E_pl_) && Util::is_nan(E_mi_)); }
      static int data_size() { return 5; }
      static String data_names() { return "E+ sigE+ E- sigE- covE+-"; }
      void data_export( xtype a[] ) const { a[0] = E_pl(); a[1] = sigE_pl(); a[2] = E_mi(); a[3] = sigE_mi(); a[4] = cov(); }
      void data_import( const xtype a[] ) { E_pl() = a[0]; sigE_pl() = a[1]; E_mi() = a[2]; sigE_mi() = a[3]; cov() = a[4]; }
      //! this type is scalable - apply magnitude scale factor
      void scale(const ftype& s) { E_pl_ *= s; sigE_pl_ *= s; E_mi_ *= s; sigE_mi_ *= s; cov_ *= (s*s); }
      // accessors
      const dtype& E_pl() const { return E_pl_; }  //<! read access
      const dtype& sigE_pl() const { return sigE_pl_; }  //<! read access
      const dtype& E_mi() const { return E_mi_; }  //<! read access
      const dtype& sigE_mi() const { return sigE_mi_; }  //<! read access
      const dtype& cov() const { return cov_; }  //<! read access
      dtype& E_pl() { return E_pl_; }  //<! write access
      dtype& sigE_pl() { return sigE_pl_; }  //<! write access
      dtype& E_mi() { return E_mi_; }  //<! write access
      dtype& sigE_mi() { return sigE_mi_; }  //<! write access
      dtype& cov() { return cov_; }  //<! write access
      // nonanomalous-type accessors
      dtype E() const { return Util::mean(E_pl_,E_mi_); }  //<! read access as simple
      dtype sigE() const { return Util::sig_mean(sigE_pl_,sigE_mi_,cov_); }  //<! read access as simple
    private:
      dtype E_pl_, E_mi_, sigE_pl_, sigE_mi_, cov_;
    };

    //! Reflection data type: F + phi model or map coeff (e.g. Fcalc, Fbest)
    template<class dtype> class F_phi : private Datatype_base
    {
    public:
      F_phi() { Util::set_null(f_); Util::set_null(phi_); }
      F_phi( const dtype& f, const dtype& phi ) : f_(f), phi_(phi) {}
      void set_null() { Util::set_null(f_); Util::set_null(phi_); }
      static String type() { return "F_phi"; }
      void friedel()
        { if (!Util::is_nan(phi_)) phi_=-phi_; }
      void shift_phase(const ftype& dphi)
        { if (!Util::is_nan(phi_)) phi_+=dphi; }
      bool missing() const
        { return (Util::is_nan(f_) || Util::is_nan(phi_)); }
      static int data_size() { return 2; }
      static String data_names() { return "F phi"; }
      void data_export( xtype array[] ) const
	{ array[0] = f(); array[1] = phi(); }
      void data_import( const xtype array[] )
	{ f() = array[0]; phi() = array[1]; }
      //! this type is scalable - apply magnitude scale factor
      void scale(const ftype& s) { f_ *= s; }
      // accessors
      const dtype& f() const { return f_; }  //<! read access
      const dtype& phi() const { return phi_; }  //<! read access
      dtype& f() { return f_; }  //<! write access
      dtype& phi() { return phi_; }  //<! write access
      //! read real part
      dtype a() const { return f_ * cos( phi_ ); }
      //! read imag part
      dtype b() const { return f_ * sin( phi_ ); }
      //! convert from complex
      F_phi(const std::complex<dtype> c) { f_=std::abs(c); phi_=std::arg(c); }
      //! convert to complex
      operator std::complex<dtype>() const { return std::polar(f_, phi_); }
      //! resolve along phase direction
      dtype resolve(const dtype phi) { return f_ * cos( phi_ - phi ); }
      //! tidy up so that real part is positive and phase 0...twopi
      const F_phi<dtype>& norm() { if ( f_ < 0.0 ) { f_ = -f_; phi_ += Util::pi(); } phi_ = Util::mod( phi_, Util::twopi() ); return *this; }
    private:
      dtype f_,phi_;
    };

    //! Reflection data type: best phi + fom
    template<class dtype> class Phi_fom : private Datatype_base
    {
    public:
      Phi_fom() { Util::set_null(phi_); Util::set_null(fom_); }
      Phi_fom( const dtype& phi, const dtype& fom ) : phi_(phi), fom_(fom) {}
      void set_null() { Util::set_null(phi_); Util::set_null(fom_); }
      static String type() { return "Phi_fom"; }
      void friedel()
        { if (!Util::is_nan(phi_)) phi_=-phi_; }
      void shift_phase(const ftype& dphi)
        { if (!Util::is_nan(phi_)) phi_+=dphi; }
      bool missing() const
        { return (Util::is_nan(phi_) || Util::is_nan(fom_)); }
      static int data_size() { return 2; }
      static String data_names() { return "phi fom"; }
      void data_export( xtype array[] ) const
	{ array[0] = phi(); array[1] = fom(); }
      void data_import( const xtype array[] )
	{ phi() = array[0]; fom() = array[1]; }
      // accessors
      const dtype& phi() const { return phi_; }  //<! read access
      const dtype& fom() const { return fom_; }  //<! read access
      dtype& phi() { return phi_; }  //<! write access
      dtype& fom() { return fom_; }  //<! write access
    private:
      dtype phi_,fom_;
    };

    //! Reflection data type: Hendrickson-Lattman coeff
    template<class dtype> class ABCD : private Datatype_base
    {
    public:
      ABCD() { Util::set_null(a_); Util::set_null(b_); Util::set_null(c_); Util::set_null(d_); }
      ABCD( const dtype& a, const dtype& b, const dtype& c, const dtype& d ) : a_(a), b_(b), c_(c), d_(d) {}
      void set_null() { Util::set_null(a_); Util::set_null(b_); Util::set_null(c_); Util::set_null(d_); }
      static String type() { return "ABCD"; }
      void friedel() { if ( !missing() ) { b_=-b_; d_=-d_; } }
      void shift_phase(const ftype& dphi)
      {
	if ( !missing() ) {
	  ftype cosd,sind;
	  dtype a1, b1, c1, d1;
	  cosd = cos(dphi);
	  sind = sin(dphi);
	  a1 = a_*cosd - b_*sind;
	  b1 = a_*sind + b_*cosd;
	  cosd = cos(2.0*dphi);
	  sind = sin(2.0*dphi);
	  c1 = c_*cosd - d_*sind;
	  d1 = c_*sind + d_*cosd;
	  a_ = a1; b_ = b1; c_ = c1; d_ = d1;
	}
      }
      bool missing() const { return (Util::is_nan(a_) || Util::is_nan(b_) || Util::is_nan(c_) || Util::is_nan(d_)); }
      static int data_size() { return 4; }
      static String data_names() { return "A B C D"; }
      void data_export( xtype array[] ) const
	{ array[0] = a(); array[1] = b(); array[2] = c(); array[3] = d(); }
      void data_import( const xtype array[] )
	{ a() = array[0]; b() = array[1]; c() = array[2]; d() = array[3]; }
      // accessors
      const dtype& a() const { return a_; }  //<! read access
      const dtype& b() const { return b_; }  //<! read access
      const dtype& c() const { return c_; }  //<! read access
      const dtype& d() const { return d_; }  //<! read access
      dtype& a() { return a_; }  //<! write access
      dtype& b() { return b_; }  //<! write access
      dtype& c() { return c_; }  //<! write access
      dtype& d() { return d_; }  //<! write access
    private:
      dtype a_,b_,c_,d_;
    };

    //! Reflection data type: Free-R flag
    class CLIPPER_IMEX Flag : private Datatype_base
    {
    public:
      Flag() { flag_ = -1; }
      explicit Flag( const int& flag ) : flag_(flag) {}
      void set_null() { flag_ = -1; }
      static String type() { return "Flag"; }
      void friedel() {}
      void shift_phase(const ftype&) {}
      bool missing() const { return (flag_ < 0); }
      static int data_size() { return 1; }
      static String data_names() { return "flag"; }
      void data_export( xtype array[] ) const
	{ array[0] = xtype(flag()); }
      void data_import( const xtype array[] )
	{ flag() = int(array[0]); }
      // accessors
      const int& flag() const { return flag_; }  //<! read access
      int& flag() { return flag_; }  //<! write access
    private:
      int flag_;
    };

    //! Reflection data type: boolean (false = missing)
    class CLIPPER_IMEX Flag_bool : private Datatype_base
    {
    public:
      Flag_bool() : flag_(false) {}
      void set_null() { flag_ = false; }
      static String type() { return "Flag_bool"; }
      void friedel() {}
      void shift_phase(const ftype&) {}
      bool missing() const { return (!flag_); }
      static int data_size() { return 1; }
      static String data_names() { return "flag"; }
      void data_export( xtype array[] ) const
	{ array[0] = xtype(flag()); }
      void data_import( const xtype array[] )
	{ flag() = bool(array[0]); }
      // accessors
      const bool& flag() const { return flag_; }  //<! read access
      bool& flag() { return flag_; }  //<! write access
    private:
      bool flag_;
    };

    //! \deprecated Anomalous difference data type: D + sigD
    /*! Provided for i/o compatibbility with legacy code only. Do not use. */
    template<class dtype> class D_sigD : private Datatype_base
    {
    public:
      D_sigD() { Util::set_null(d_); Util::set_null(sigd_); }
      D_sigD( const dtype& d, const dtype& sigd ) : d_(d), sigd_(sigd) {}
      void set_null() { Util::set_null(d_); Util::set_null(sigd_); }
      static String type() { return "D_sigD"; }
      void friedel() { d_ = -d_; }
      void shift_phase(const ftype&) {}
      bool missing() const { return (Util::is_nan(d_) || Util::is_nan(sigd_)); }
      static int data_size() { return 2; }
      static String data_names() { return "Dano sigDano"; }
      void data_export( xtype array[] ) const
	{ array[0] = d(); array[1] = sigd(); }
      void data_import( const xtype array[] )
	{ d() = array[0]; sigd() = array[1]; }
      //! this type is scalable - apply magnitude scale factor
      void scale(const ftype& s) { d_ *= s; sigd_ *= s; }
      // accessors
      const dtype& d() const { return d_; }  //<! read access
      const dtype& sigd() const { return sigd_; }  //<! read access
      dtype& d() { return d_; }  //<! write access
      dtype& sigd() { return sigd_; }  //<! write access
    private:
      dtype d_,sigd_;
    };

  }


  namespace data32
  {
    typedef clipper::datatypes::I_sigI<ftype32> I_sigI;          //!< datatype
    typedef clipper::datatypes::I_sigI_ano<ftype32> I_sigI_ano;  //!< datatype
    typedef clipper::datatypes::F_sigF<ftype32> F_sigF;          //!< datatype
    typedef clipper::datatypes::F_sigF_ano<ftype32> F_sigF_ano;  //!< datatype
    typedef clipper::datatypes::E_sigE<ftype32> E_sigE;          //!< datatype
    typedef clipper::datatypes::F_phi<ftype32> F_phi;            //!< datatype
    typedef clipper::datatypes::Phi_fom<ftype32> Phi_fom;        //!< datatype
    typedef clipper::datatypes::ABCD<ftype32> ABCD;              //!< datatype
    typedef clipper::datatypes::Flag Flag;                       //!< datatype
    typedef clipper::datatypes::Flag_bool Flag_bool;             //!< datatype
    typedef clipper::datatypes::D_sigD<ftype32> D_sigD;          //!< datatype
  }

  namespace data64
  {
    typedef clipper::datatypes::I_sigI<ftype64> I_sigI;          //!< datatype
    typedef clipper::datatypes::I_sigI_ano<ftype64> I_sigI_ano;  //!< datatype
    typedef clipper::datatypes::F_sigF<ftype64> F_sigF;          //!< datatype
    typedef clipper::datatypes::F_sigF_ano<ftype64> F_sigF_ano;  //!< datatype
    typedef clipper::datatypes::E_sigE<ftype64> E_sigE;          //!< datatype
    typedef clipper::datatypes::F_phi<ftype64> F_phi;            //!< datatype
    typedef clipper::datatypes::Phi_fom<ftype64> Phi_fom;        //!< datatype
    typedef clipper::datatypes::ABCD<ftype64> ABCD;              //!< datatype
    typedef clipper::datatypes::Flag Flag;                       //!< datatype
    typedef clipper::datatypes::Flag_bool Flag_bool;             //!< datatype
    typedef clipper::datatypes::D_sigD<ftype64> D_sigD;          //!< datatype
  }



} // namespace clipper

#endif
