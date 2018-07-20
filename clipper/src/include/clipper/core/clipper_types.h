/*! \file lib/clipper_types.h
    Header file for clipper basic types
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


#ifndef CLIPPER_TYPES
#define CLIPPER_TYPES


#include "clipper_util.h"
#include "clipper_memory.h"

#if defined(_MSC_VER) && _MSC_VER >= 1600
// Workaround to avoid DLL export errors
const std::basic_string<char>::size_type std::basic_string<char>::npos = (std::basic_string<char>::size_type) - 1;
#endif

#include "../imex.h"

namespace clipper
{
  // forward definitions
  template<class T> class Mat33sym;


  //! String extension with simple parsing methods
  /*!
    String extension with primitive 'split' operation for parsing
    and pathname processing operations.
  */

  class CLIPPER_IMEX String : public std::string
  {
  public:
    inline String() : std::string() {}  //!< null constructor
    inline String( const std::string str ) : std::string( str ) {} //!< constructor: from string
    inline String( const char* str ) : std::string( str ) {} //!< constructor: from char*
    String( const char* str, const int l );     //!< constructor: from char*
    String( const int i, const int w = 4 );     //!< constructor: from int
    String( const long i, const int w = 4 );    //!< constructor: from long
    //! constructor: from float
    String( const float f, const int w = 6, const int p = 6 );
    //! constructor: from double
    String( const double f, const int w = 6, const int p = 6 );

    //! String splitter - a very simple parser component
    std::vector<String> split(const String sep) const;
    //! Return copy of string without leading and trailing blanks
    String trim() const;

    //! get trailing path element
    String tail() const;
    //! remove trailing path element
    String head() const;
    //! get leading path element
    String nohead() const;
    //! remove leading path element
    String notail() const;

    //! construct string from rational f using base b
    static String rational( const double f, const int b, const bool sign=false );

    int i() const;     //!< convert to int
    long l() const;    //!< convert to long
    ftype32 f32() const;  //!< convert to float
    ftype64 f64() const;  //!< convert to double
    ftype f() const;      //!< convert to ftype
    ftype rational() const;  //!< convert from rational to ftype
  };
 

  //! 3-vector class
  template<class T = ftype> class Vec3
  {
  public:
    //! null constructor
    inline Vec3() {}
    //! constructor: from individual values
    inline Vec3( const T& v0, const T& v1, const T& v2 )
      { vec[0] = v0; vec[1] = v1; vec[2] = v2; }
    //! constructor: copy/convert
    template<class TT> explicit Vec3( const Vec3<TT>& v )
      { vec[0] = TT(v[0]); vec[1] = TT(v[1]); vec[2] = TT(v[2]); }
    //! test equality
    bool equals( const Vec3<T>& v, const T& tol ) const;
    //! get element
    inline const T& operator []( const int& i ) const { return vec[i]; }
    //! set element
    inline T& operator []( const int& i ) { return vec[i]; }
    //! return unit vector with same direction as this vector
    inline Vec3<T> unit() const
      { return (*this)*T(1.0/sqrt(ftype(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]))); }
    //! return zero vector
    inline static Vec3<T> zero() { return Vec3<T>( 0, 0, 0 ); }
    //! return null vector (only valid for floating point types)
    inline static Vec3<T> null() { Vec3<T> v( 0, 0, 0 ); Util::set_null(v.vec[0]); return v; }
    //! test for null vector
    inline bool is_null() const { return Util::is_nan( (float)vec[0] ); }
    //! Vector dot product (equivalent to *)
    inline static T dot( const Vec3<T>& v1, const Vec3<T>& v2 )
      { return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]); }
    //! Vector cross product
    inline static Vec3<T> cross( const Vec3<T>& v1, const Vec3<T>& v2 )
      { return Vec3<T>(v1[1]*v2[2]-v2[1]*v1[2],
		       v1[2]*v2[0]-v2[2]*v1[0],
		       v1[0]*v2[1]-v2[0]*v1[1]); }
    //! return formatted String representation
    String format() const
      { return "("+String((float)vec[0],10,4)+","+String((float)vec[1],10,4)+","+String((float)vec[2],10,4)+")"; }
    //! add another vector to this one
    inline const Vec3<T>& operator += ( const Vec3<T>& v )
      { vec[0] += v[0]; vec[1] += v[1]; vec[2] += v[2]; return (*this); }
    //! subtract another vector from this one
    inline const Vec3<T>& operator -= ( const Vec3<T>& v )
      { vec[0] -= v[0]; vec[1] -= v[1]; vec[2] -= v[2]; return (*this); }
    //! Vector equality (for floating point types see equals())
    //-- friend int operator == ( const Vec3<T>& v1, const Vec3<T>& v2 );
    //! Vector inequality (for floating point types see equals())
    //-- friend int operator != ( const Vec3<T>& v1, const Vec3<T>& v2 );
    //! Vector negation operator
    //-- friend Vec3<T> operator -( const Vec3<T>& v );
    //! Vector addition operator
    //-- friend Vec3<T> operator +( const Vec3<T>& v1, const Vec3<T> &v2 );
    //! Vector subtraction operator
    //-- friend Vec3<T> operator -( const Vec3<T>& v1, const Vec3<T>& v2 );
    //! Vector scaling operator
    //-- friend Vec3<T> operator *( const T& s, const Vec3<T>& v1 );
    //! Vector scaling operator
    //-- friend Vec3<T> operator *( const Vec3<T>& v1, const T& s );
    //! Vector dot product
    //-- friend T operator *( const Vec3<T>& v1, const Vec3<T>& v2 );
  private:
    T vec[3];
  };
  template<class T> inline int operator == ( const Vec3<T>& v1, const Vec3<T>& v2 ) { return (v1[0]==v2[0] && v1[1]==v2[1] && v1[2]==v2[2]); }
  template<class T> inline int operator != ( const Vec3<T>& v1, const Vec3<T>& v2 ) { return (v1[0]!=v2[0] || v1[1]!=v2[1] || v1[2]!=v2[2]); }
  template<class T> inline Vec3<T> operator -( const Vec3<T>& v )
    { return Vec3<T>( -v[0], -v[1], -v[2] ); }
  template<class T> inline Vec3<T> operator +( const Vec3<T>& v1, const Vec3<T> &v2 ) { return Vec3<T>( v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]); }
  template<class T> inline Vec3<T> operator -( const Vec3<T>& v1, const Vec3<T>& v2 ) { return Vec3<T>( v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]); }
  template<class T> inline Vec3<T> operator *( const T& s, const Vec3<T>& v1 )
    { return Vec3<T>( s*v1[0], s*v1[1], s*v1[2]); }
  template<class T> inline Vec3<T> operator *( const Vec3<T>& v1, const T& s )
    { return Vec3<T>( s*v1[0], s*v1[1], s*v1[2]); }
  template<class T> inline T operator *( const Vec3<T>& v1, const Vec3<T>& v2 )
    { return Vec3<T>::dot(v1,v2); }


  //! 3x3-matrix class
  template<class T = ftype> class Mat33
  {
  public:
    //! null constructor
    inline Mat33() {}
    //! constructor
    inline Mat33( const T& m00, const T& m01, const T& m02, const T& m10, const T& m11, const T& m12, const T& m20, const T& m21, const T& m22 )
      { mat[0][0] = m00; mat[0][1] = m01; mat[0][2] = m02; mat[1][0] = m10; mat[1][1] = m11; mat[1][2] = m12; mat[2][0] = m20; mat[2][1] = m21; mat[2][2] = m22; }
    //! constructor: copy/convert
    template<class TT> explicit Mat33( const Mat33<TT>& m );
    //! constructor: copy/convert from symmetric matrix
    template<class TT> explicit Mat33( const Mat33sym<TT>& m );
    T det() const;               //!< determinant
    Mat33<T> inverse() const;    //!< inverse
    Mat33<T> transpose() const;  //!< transpose
    bool equals( const Mat33<T>& m, const T& tol ) const;  //!< test equality
    inline const T& operator ()( const int& i, const int& j ) const
      { return mat[i][j]; }      //!< get element
    inline T& operator ()( const int& i, const int& j )
      { return mat[i][j]; }      //!< set element
    //! return formatted String representation
    String format() const
      { return "|"+String((float)mat[0][0],10,4)+","+String((float)mat[0][1],10,4)+","+String((float)mat[0][2],10,4)+"|\n|"+String((float)mat[1][0],10,4)+","+String((float)mat[1][1],10,4)+","+String((float)mat[1][2],10,4)+"|\n|"+String((float)mat[2][0],10,4)+","+String((float)mat[2][1],10,4)+","+String((float)mat[2][2],10,4)+"|"; }
    //! return identity matrix
    inline static Mat33<T> identity() { Mat33<T> m; m.mat[0][0] = m.mat[1][1] = m.mat[2][2] = 1; m.mat[0][1] = m.mat[0][2] = m.mat[1][0] = m.mat[1][2] = m.mat[2][0] = m.mat[2][1] = 0; return m; }
    //! return null matrix (only valid for floating point types)
    inline static Mat33<T> null() { Mat33<T> m; m.mat[0][0] = Util::nan(); return m; }
    //! test for null matrix (only valid for floating point types)
    bool inline is_null() const { return Util::is_nan( (float)mat[0][0] ); }
    //! Matrix-vector product
    /*! Assumes a column vector */
    //-- friend Vec3<T> operator *( const Mat33<T>& m, const Vec3<T>& v );
    //! Vector-matrix product
    /*! Assumes a row vector, i.e. equivalent to the matrix-vector
      product with the matrix trasposed */
    //-- friend Vec3<T> operator *( const Vec3<T>& v, const Mat33<T>& m );
    //! Matrix-matrix product
    //-- friend Mat33<T> operator *(const Mat33<T>& m1, const Mat33<T>& m2);
    //! Matrix sum
    //-- friend Mat33<T> operator +(const Mat33<T>& m1, const Mat33<T>& m2);
    //! Unary minus
    //-- friend Mat33<T> operator -(const Mat33<T>& m);
  private:
    T mat[3][3];
  };
  template<class T> inline Vec3<T> operator *( const Mat33<T>& m, const Vec3<T>& v )
    { return Vec3<T>( m(0,0)*v[0]+m(0,1)*v[1]+m(0,2)*v[2],
		      m(1,0)*v[0]+m(1,1)*v[1]+m(1,2)*v[2],
		      m(2,0)*v[0]+m(2,1)*v[1]+m(2,2)*v[2] ); }
  template<class T> inline Vec3<T> operator *( const Vec3<T>& v, const Mat33<T>& m )
    { return Vec3<T>( v[0]*m(0,0)+v[1]*m(1,0)+v[2]*m(2,0),
		      v[0]*m(0,1)+v[1]*m(1,1)+v[2]*m(2,1),
		      v[0]*m(0,2)+v[1]*m(1,2)+v[2]*m(2,2) ); }
  template<class T> inline Mat33<T> operator *( const Mat33<T>& m1, const Mat33<T>& m2 )
    {
      return Mat33<T> ( m1(0,0)*m2(0,0) + m1(0,1)*m2(1,0) + m1(0,2)*m2(2,0),
			m1(0,0)*m2(0,1) + m1(0,1)*m2(1,1) + m1(0,2)*m2(2,1),
			m1(0,0)*m2(0,2) + m1(0,1)*m2(1,2) + m1(0,2)*m2(2,2),
			m1(1,0)*m2(0,0) + m1(1,1)*m2(1,0) + m1(1,2)*m2(2,0),
			m1(1,0)*m2(0,1) + m1(1,1)*m2(1,1) + m1(1,2)*m2(2,1),
			m1(1,0)*m2(0,2) + m1(1,1)*m2(1,2) + m1(1,2)*m2(2,2),
			m1(2,0)*m2(0,0) + m1(2,1)*m2(1,0) + m1(2,2)*m2(2,0),
			m1(2,0)*m2(0,1) + m1(2,1)*m2(1,1) + m1(2,2)*m2(2,1),
			m1(2,0)*m2(0,2) + m1(2,1)*m2(1,2) + m1(2,2)*m2(2,2) );
    }
  template<class T> inline Mat33<T> operator +( const Mat33<T>& m1, const Mat33<T>& m2 )
    { return Mat33<T>( m1(0,0)+m2(0,0), m1(0,1)+m2(0,1), m1(0,2)+m2(0,2),
		       m1(1,0)+m2(1,0), m1(1,1)+m2(1,1), m1(1,2)+m2(1,2),
		       m1(2,0)+m2(2,0), m1(2,1)+m2(2,1), m1(2,2)+m2(2,2) ); }
  template<class T> inline Mat33<T> operator -( const Mat33<T>& m )
    { return Mat33<T>( -m(0,0), -m(0,1), -m(0,2),
		       -m(1,0), -m(1,1), -m(1,2),
		       -m(2,0), -m(2,1), -m(2,2) ); }


  //! Compressed form for 3x3 symmetric matrix class
  template<class T = ftype> class Mat33sym
  {
  public:
    //! null constructor
    inline Mat33sym() {}
    //! constructor: from Mat33 (does not check for symmetry)
    template<class TT> explicit Mat33sym( const Mat33<TT>& m ) :
      m00(m(0,0)), m11(m(1,1)), m22(m(2,2)),
      m01(m(0,1)), m02(m(0,2)), m12(m(1,2)) {}
    //! constructor: from Mat33sym
    template<class TT> explicit Mat33sym( const Mat33sym<TT>& m ) :
      m00(m.mat00()), m11(m.mat11()), m22(m.mat22()),
      m01(m.mat01()), m02(m.mat02()), m12(m.mat12()) {}
    //! constructor: from coefficients
    inline Mat33sym( const T& c00, const T& c11, const T& c22,
		     const T& c01, const T& c02, const T& c12 ) :
      m00(c00), m11(c11), m22(c22), m01(c01), m02(c02), m12(c12) {}
    //! return formatted String representation
    String format() const
      { return "|"+String(m00,10,4)+","+String(m01,10,4)+","+String(m02,10,4)+"|\n|"+String(m01,10,4)+","+String(m11,10,4)+","+String(m12,10,4)+"|\n|"+String(m02,10,4)+","+String(m12,10,4)+","+String(m22,10,4)+"|"; }
    //! return identity matrix
    inline static Mat33sym<T> identity()
      { return Mat33sym<T>( 1, 1, 1, 0, 0, 0 ); }
    //! return null matrix (only valid for floating point types)
    inline static Mat33sym<T> null()
      { return Mat33sym<T>(Util::nan(),0,0,0,0,0); }
    //! test for null matrix (only valid for floating point types)
    inline bool is_null() const { return Util::is_nan( m00 ); }
    //! return quadratic form with vector
    T quad_form( const Vec3<T>& v ) const;
    T det() const;               //!< determinant
    Mat33<T> sqrt() const;       //!< square root
    Mat33sym<T> inverse() const; //!< inverse
    inline const T& mat00() const { return m00; }  //!< element (0,0)
    inline const T& mat11() const { return m11; }  //!< element (1,1)
    inline const T& mat22() const { return m22; }  //!< element (2,2)
    inline const T& mat01() const { return m01; }  //!< element (0,1)
    inline const T& mat02() const { return m02; }  //!< element (0,2)
    inline const T& mat12() const { return m12; }  //!< element (1,2)
    //! access elements as 3x3 matrix (inefficient)
    const T& operator ()( const int& i, const int& j ) const;
    //! Matrix-vector product
    //-- friend Vec3<T> operator *( const Mat33sym<T>& m, const Vec3<T>& v );
    //! Matrix sum
    //-- friend Mat33sym<T> operator +( const Mat33sym<T>& m1, const Mat33sym<T>& m2 );
    //! Unary minus
    //-- friend Mat33sym<T> operator -( const Mat33sym<T>& m );
   private:
    T m00, m11, m22, m01, m02, m12;
  };
  template<class T> inline Vec3<T> operator *( const Mat33sym<T>& m, const Vec3<T>& v )
    { return Vec3<T>( m.mat00()*v[0]+m.mat01()*v[1]+m.mat02()*v[2],
		      m.mat01()*v[0]+m.mat11()*v[1]+m.mat12()*v[2],
		      m.mat02()*v[0]+m.mat12()*v[1]+m.mat22()*v[2] ); }
  template<class T> inline Mat33sym<T> operator +( const Mat33sym<T>& m1, const Mat33sym<T>& m2 )
    { return Mat33sym<T>( m1.mat00()+m2.mat00(), m1.mat11()+m2.mat11(),
			  m1.mat22()+m2.mat22(), m1.mat01()+m2.mat01(),
			  m1.mat02()+m2.mat02(), m1.mat12()+m2.mat12() ); }
  template<class T> inline Mat33sym<T> operator -( const Mat33sym<T>& m )
    { return Mat33sym<T>( -m.mat00(), -m.mat11(), -m.mat22(),
			  -m.mat01(), -m.mat02(), -m.mat12() ); }


  //! Rotation-translation operator
  template<class T = ftype> class RTop
  {
  public:
    //! null constructor
    inline RTop() {}
    //! constructor: from rotation
    inline explicit RTop( const Mat33<T>& r ) : rot_( r ), trn_( Vec3<T>::zero() ) {}
    //! constructor: from rotation and translation
    inline RTop( const Mat33<T>& r, const Vec3<T>& t ) : rot_( r ), trn_( t ) {}
    //! inverse
    RTop<T> inverse() const
      { Mat33<T> minv = rot().inverse(); return RTop<T>(minv, -(minv*trn())); }
    //! test equality with some tolerance
    inline bool equals( const RTop<T>& m, const T& tol ) const
      { return ( rot().equals(m.rot(),tol) && trn().equals(m.trn(),tol) ); }
    inline const Mat33<T>& rot() const { return rot_; }  //!< get rotation
    inline const Vec3<T>&  trn() const { return trn_; }  //!< get translation
    inline Mat33<T>& rot() { return rot_; } //!< set rotation
    inline Vec3<T>&  trn() { return trn_; } //!< set translation
    //! return identity operator
    inline static RTop<T> identity()
      { return RTop<T>( Mat33<T>::identity(), Vec3<T>::zero() ); }
    //! return identity operator
    inline static RTop<T> null()
      { return RTop<T>( Mat33<T>::null(), Vec3<T>::null() ); }
    //! test for null operator
    inline bool is_null() const { return rot_.is_null() || trn_.is_null(); }
    //! return formatted String representation
    String format() const { return rot_.format() + "\n" + trn_.format(); }
    //! apply RTop to vector
    //-- friend Vec3<T> operator *( const RTop<T>& r, const Vec3<T>& v );
    //! RTop product
    //-- friend RTop<T> operator *( const RTop<T>& r1, const RTop<T>& r2 );
  private:
    Mat33<T> rot_; Vec3<T> trn_;
  };
  template<class T> inline Vec3<T> operator *( const RTop<T>& r, const Vec3<T>& v ) { return r.rot()*v + r.trn(); }
  template<class T> inline RTop<T> operator *( const RTop<T>& r1, const RTop<T>& r2 ) { return RTop<T>( r1.rot()*r2.rot(), r1.rot()*r2.trn()+r1.trn() ); }



  //! Simple 2-d array class
  template<class T = ftype> class Array2d
  {
  public:
    //! null constructor
    inline Array2d() { d1_ = d2_ = 0; }
    //! constructor
    inline Array2d( const int& d1, const int& d2 ) { resize( d1, d2 ); }
    //! constructor
    inline Array2d( const int& d1, const int& d2, T val )
      { resize( d1, d2, val ); }
    //! resize
    void inline resize( const int& d1, const int& d2 )
      { d1_ = d1; d2_ = d2; data.resize( size() ); }
    //! resize
    void inline resize( const int& d1, const int& d2, const T& val )
      { d1_ = d1; d2_ = d2; data.resize( size(), val ); }
    inline int size() const { return d1_ * d2_; }    //!< size
    inline const int& rows() const { return d1_; }  //!< number of rows
    inline const int& cols() const { return d2_; }  //!< number of cols
    //! read accessor
    inline const T& operator () ( const int& i1, const int& i2 ) const
      { return data[ i1*d2_ + i2 ]; }
    //! write accessor
    inline T& operator () ( const int& i1, const int& i2 )
      { return data[ i1*d2_ + i2 ]; }
  protected:
    std::vector<T> data;
    int d1_, d2_;
  };


  //! General matrix class: like Array2d but with numerical methods
  template<class T = ftype> class Matrix : public Array2d<T>
  {
  public:
    //! null constructor
    inline Matrix() {}
    //! constructor
    inline Matrix( const int& d1, const int& d2 )
      { Array2d<T>::resize( d1, d2 ); }
    //! constructor
    inline Matrix( const int& d1, const int& d2, T val )
      { Array2d<T>::resize( d1, d2, val ); }
    //! equation solver (square matrices only)
    std::vector<T> solve( const std::vector<T>& b ) const;
    //! eigenvalue calculation (square symmetric matrices only)
    std::vector<T> eigen( const bool sort = true );
    //! Matrix-vector product
    /*! Assumes a column vector */
    //-- friend std::vector<T> operator *( const Matrix<T>& m, const std::vector<T>& v );
  };
  template<class T> std::vector<T> operator *( const Matrix<T>& m, const std::vector<T>& v )
    {
      if ( m.cols() != v.size() )
	Message::message( Message_fatal( "Matrix*vector dimension mismatch" ) );
      std::vector<T> r( m.rows() );
      int i, j; T s;
      for ( i = 0; i < m.rows(); i++ ) {
	s = T(0);
	for ( j = 0; j < m.cols(); j++ ) s += m(i,j) * v[j];
	r[i] = s;
      }
      return r;
    }



  // template implementations

  template<class T> bool Vec3<T>::equals( const Vec3<T>& v, const T& tol ) const
  {
    return ( pow(vec[0]-v[0],T(2)) + pow(vec[1]-v[1],T(2)) + 
	     pow(vec[2]-v[2],T(2)) <= pow(tol,T(2)) );
  }

  template<class T> template<class TT> Mat33<T>::Mat33( const Mat33<TT>& m )
  {
    mat[0][0]=T(m(0,0)); mat[0][1]=T(m(0,1)); mat[0][2]=T(m(0,2));
    mat[1][0]=T(m(1,0)); mat[1][1]=T(m(1,1)); mat[1][2]=T(m(1,2));
    mat[2][0]=T(m(2,0)); mat[2][1]=T(m(2,1)); mat[2][2]=T(m(2,2));
  }

  template<class T> template<class TT> Mat33<T>::Mat33( const Mat33sym<TT>& m )
  {
    mat[0][0]=T(m.mat00());
    mat[1][1]=T(m.mat11());
    mat[2][2]=T(m.mat22());
    mat[0][1]=mat[1][0]=T(m.mat01());
    mat[0][2]=mat[2][0]=T(m.mat02());
    mat[1][2]=mat[2][1]=T(m.mat12());
  }

  template<class T> bool Mat33<T>::equals( const Mat33<T>& m, const T& tol ) const
  {
    return ( pow(mat[0][0]-m(0,0),T(2)) + pow(mat[0][1]-m(0,1),T(2)) +
	     pow(mat[0][2]-m(0,2),T(2)) + pow(mat[1][0]-m(1,0),T(2)) +
	     pow(mat[1][1]-m(1,1),T(2)) + pow(mat[1][2]-m(1,2),T(2)) +
	     pow(mat[2][0]-m(2,0),T(2)) + pow(mat[2][1]-m(2,1),T(2)) +
	     pow(mat[2][2]-m(2,2),T(2)) <= pow(tol,T(2)) );
  }

  template<class T> T Mat33<T>::det() const
  {
    return ( mat[0][0]*(mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1]) + 
	     mat[0][1]*(mat[1][2]*mat[2][0] - mat[1][0]*mat[2][2]) +
	     mat[0][2]*(mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0]) );
  }

  template<class T> Mat33<T> Mat33<T>::inverse() const
  {
    T d = det();
    Mat33<T> inv;
    inv(0,0) = ( mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1] ) / d;
    inv(1,0) = ( mat[1][2]*mat[2][0] - mat[1][0]*mat[2][2] ) / d;
    inv(2,0) = ( mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0] ) / d;
    inv(0,1) = ( mat[2][1]*mat[0][2] - mat[2][2]*mat[0][1] ) / d;
    inv(1,1) = ( mat[2][2]*mat[0][0] - mat[2][0]*mat[0][2] ) / d;
    inv(2,1) = ( mat[2][0]*mat[0][1] - mat[2][1]*mat[0][0] ) / d;
    inv(0,2) = ( mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1] ) / d;
    inv(1,2) = ( mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2] ) / d;
    inv(2,2) = ( mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0] ) / d;
    return inv;
  }

  template<class T> Mat33<T> Mat33<T>::transpose() const
  {
    Mat33<T> t;
    t(0,0) = mat[0][0]; t(0,1) = mat[1][0]; t(0,2) = mat[2][0];
    t(1,0) = mat[0][1]; t(1,1) = mat[1][1]; t(1,2) = mat[2][1];
    t(2,0) = mat[0][2]; t(2,1) = mat[1][2]; t(2,2) = mat[2][2];
    return t;
  }

  template<class T> T Mat33sym<T>::det() const
  {
    return ( m00*(m11*m22 - m12*m12) + m01*(m12*m02 - m01*m22) +
	     m02*(m01*m12 - m11*m02) );
  }

  template<class T> Mat33<T> Mat33sym<T>::sqrt() const
  {
    Mat33<T> half( Mat33sym<T>( 0.5, 0.5, 0.5, 0.0, 0.0, 0.0 ) );
    Mat33<T> target( *this );
    Mat33<T> result( target );
    result(1,0) = result(2,0) = result(2,1) = 0.0;
    for ( int i = 0; i < 10; i++ )
      result = half * ( result.inverse() * target + result );
    return result;
  }

  template<class T> Mat33sym<T> Mat33sym<T>::inverse() const
  {
    T d = det();
    return Mat33sym<T> ( ( m11*m22 - m12*m12 ) / d,
			 ( m22*m00 - m02*m02 ) / d,
			 ( m00*m11 - m01*m01 ) / d,
			 ( m12*m02 - m22*m01 ) / d,
			 ( m01*m12 - m02*m11 ) / d,
			 ( m02*m01 - m00*m12 ) / d );
  }

  template<class T> T Mat33sym<T>::quad_form( const Vec3<T>& v ) const
  {
    return ( v[0]*( v[0]*m00 + 2*(v[1]*m01+v[2]*m02) ) +
	     v[1]*( v[1]*m11 + 2*(v[2]*m12) ) + v[2]*v[2]*m22 );
  }

  template<class T> const T& Mat33sym<T>::operator ()( const int& i, const int& j ) const
  {
    switch (i) {
    case 0:
      switch (j) {
      case 0: return m00;
      case 1: return m01;
      default: return m02;
      }
    case 1:
      switch (j) {
      case 0: return m01;
      case 1: return m11;
      default: return m12;
      }
    default:
      switch (j) {
      case 0: return m02;
      case 1: return m12;
      default: return m22;
      }
    }
  }


// complex matrix methods


/*! Solve the system of linear equations Ax=b for x
  Uses elimination. Only suitable for small systems. */
template<class T> std::vector<T> Matrix<T>::solve( const std::vector<T>& b ) const
{
  if ( Array2d<T>::rows() != Array2d<T>::cols() )
    Message::message( Message_fatal("Matrix.solve() matrix not square") );
  if ( Array2d<T>::rows() != b.size() )
    Message::message( Message_fatal("Matrix.solve() matrix/vector mismatch") );
  const int n = Array2d<T>::rows();

  // solve for X by Gaussian elimination
  T s, pivot;
  int i, j, k;

  Matrix<T> a = *this;
  std::vector<T> x = b;
  for ( i = 0; i < n; i++ ) {
    // pick largest pivot
    j = i;
    for ( k = i+1; k < n; k++ )
      if ( fabs(a(k,i)) > fabs(a(j,i)) ) j = k;
    // swap rows
    for ( k = 0; k < n; k++ )
      Util::swap( a(i,k), a(j,k) );
    Util::swap( x[i], x[j] );
    // perform elimination
    pivot = a(i,i);
    for ( j = 0; j < n; j++ ) {
      if ( j != i ) {
        s = a(j,i) / pivot;
        for ( k = i+1; k < n; k++ ) a(j,k) = a(j,k) - s*a(i,k);
        x[j] = x[j] - s*x[i];
      }
    }
  }
  for ( i = 0; i < n; i++ ) x[i] /= a(i,i);
  return x;
}


/*! Find the Eigenvalues and Eigenvectors of the matrix.
  Uses the Jacobi method. Only suitable for small systems (dimension<20).
  The matrix is replaced by the matrix of eigenvectors (as columns).
  \param sort Sort the eigenvalues and vectors, smallest first. (default=true)
  \return Vector of eigenvalues. */
template<class T> std::vector<T> Matrix<T>::eigen( const bool sort )
{
  if ( Array2d<T>::rows() != Array2d<T>::cols() )
    Message::message( Message_fatal( "Matrix.eigen() matrix not square" ) );
  const int n = Array2d<T>::rows();

  int cyc, j, q, p;
  T spp, spq, t, s, c, theta, tau, h, ap, aq, a_pq;

  Matrix<T>& mat = *this;
  Matrix evec( n, n, 0.0 );
  std::vector<T> eval( n );
  std::vector<T> b( n );
  std::vector<T> z( n );

  // Set evec to identity, eval & b to diagonal, z to 0.
  for( p = 0; p < n; p++ ) {
    evec(p,p) = 1.0;
    eval[p] = b[p] = mat(p,p);
  }

  for ( cyc = 1; cyc <= 50; cyc++ ) {

    // calc sum of diagonal, off-diagonal
    spp = spq = 0.0;
    for ( p=0; p<n-1; p++ ) {
      for ( q=p+1; q<n; q++ )
	spq += fabs( mat(p,q) );
      spp += fabs( mat(p,p) );
    }
    if ( spq <= 1.0e-12 * spp ) break;

    // zero z
    for ( p = 0; p < n; p++ ) z[p] = 0.0;

    // now try and reduce each off-diagonal element in turn
    for( p=0; p<n-1; p++ ) {
      for( q=p+1; q<n; q++ ) {
	a_pq = mat( p, q );
	h = eval[q] - eval[p];
	if ( fabs(a_pq) > 1.0e-12*fabs(h) ) {
	  theta = 0.5*h/a_pq;
	  t = 1.0/(fabs(theta) + sqrt(1.0 + theta*theta));
	  if ( theta < 0.0 ) t = -t;
	} else {
	  t = a_pq/h;
	}

	// calc trig properties
	c   = 1.0/sqrt(1.0+t*t);
	s   = t*c;
	tau = s/(1.0+c);
	h   = t * a_pq;

	// update eigenvalues
	z[p] -= h;
	z[q] += h;
	eval[p] -= h;
	eval[q] += h;

	// rotate the upper diagonal of the matrix
	mat( p, q ) = 0.0;
	for ( j = 0; j < p; j++ ) {
	  ap = mat( j, p );
	  aq = mat( j, q );
	  mat( j, p ) = ap - s * ( aq + ap * tau );
	  mat( j, q ) = aq + s * ( ap - aq * tau );
	}
	for ( j = p+1; j < q; j++ ) {
	  ap = mat( p, j );
	  aq = mat( j, q );
	  mat( p, j ) = ap - s * ( aq + ap * tau );
	  mat( j, q ) = aq + s * ( ap - aq * tau );
	}
	for ( j = q+1; j < n; j++ ) {
	  ap = mat( p, j );
	  aq = mat( q, j );
	  mat( p, j ) = ap - s * ( aq + ap * tau );
	  mat( q, j ) = aq + s * ( ap - aq * tau );
	}
	// apply corresponding rotation to result
	for ( j = 0; j < n; j++ ) {
	  ap = evec( j, p );
	  aq = evec( j, q );
	  evec( j, p ) = ap - s * ( aq + ap * tau );
	  evec( j, q ) = aq + s * ( ap - aq * tau );
	}
      }
    }

    for ( p = 0; p < n; p++ ) {
      b[p] += z[p];
      eval[p] = b[p];
    }
  }

  // sort the eigenvalues
  if ( sort ) {
    for ( p = 0; p < n; p++ ) {
      j = p;        // set j to index of largest remaining eval
      for ( q = p+1; q < n; q++ )
	if ( eval[q] < eval[j] ) j = q;
      Util::swap( eval[p], eval[j] );  // now swap evals, evecs
      for ( q = 0; q < n; q++ )
	Util::swap( evec( q, p ), evec( q, j ) );
    }
  }

  // return eigenvalues
  mat = evec;
  return eval;
}

} // namespace clipper

#endif
