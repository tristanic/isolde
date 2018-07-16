/*! \file lib/rotation.h
    Fundamental types for the clipper libraries
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


#ifndef CLIPPER_ROTATION
#define CLIPPER_ROTATION


#include "clipper_types.h"
#include "../imex.h"

namespace clipper
{

  // forward definition
  class Rotation;


  //! Euler angle class
  /* Rotations are generally handled through the clipper::Rotation
     class. This class only exists for conversion purposes.

     This particular class represents generic Euler angles. The
     convention is selected from the 24 possible conventions according
     to the template parameter. The integer convention code is
     enumerated in the Rotation::EULERtype enumation in the form
     Rotation::EulerZYZr, Rotation::EulerXYZs etc., where the X/Y/Z
     indicates the axes of rotation in order, and the r/s indicates
     static or rotating axes. The type of an Euler class is also given
     as a prefix to the result of format(). */
  template<int T> class Euler {
  public:
    //! constructor: null
    Euler() {}
    //! constructor: from specified angles
    Euler( const ftype& alpha, const ftype& beta,  const ftype& gamma ) :
      alpha_(alpha), beta_(beta), gamma_(gamma) {}
    //! constructor: from rotation
    Euler( const Rotation& rot );
    //! return rotation
    Rotation rotation() const;
    const ftype& alpha() const { return alpha_; }  //!< return alpha
    const ftype& beta()  const { return beta_;  }  //!< return beta
    const ftype& gamma() const { return gamma_; }  //!< return gamma
    String format() const;  //!< return formatted String representation
  private:
    static void params( int& r1, int& r2, int& r3, int& s );
    ftype alpha_, beta_, gamma_;
  };

  //! Euler_ccp4 angle class
  /* Rotations are generally handled through the clipper::Rotation
     class. This class only exists for conversion purposes.

     This particular class represents Euler_ccp4 angles according to the
     CCP4 standard, i.e.
     - Rotation 1 (alpha) about K,
     - Rotation 2 (beta) about the new J,
     - Rotation 3 (gamma) about the new K. */
  class CLIPPER_IMEX Euler_ccp4 {
  public:
    //! constructor: null
    Euler_ccp4() {}
    //! constructor: from specified angles
    Euler_ccp4( const ftype& alpha, const ftype& beta,  const ftype& gamma ) :
      alpha_(alpha), beta_(beta), gamma_(gamma) {}
    const ftype& alpha() const { return alpha_; }  //!< return alpha
    const ftype& beta()  const { return beta_;  }  //!< return beta
    const ftype& gamma() const { return gamma_; }  //!< return gamma
    String format() const;  //!< return formatted String representation
  private:
    ftype alpha_, beta_, gamma_;
  };

  //! Polar_ccp4 angle class
  /* Rotations are generally handled through the clipper::Rotation
     class. This class only exists for conversion purposes.

     This particular class represents Polar_ccp4 angles according to the
     CCP4 standard, i.e.
     - omega gives inclination of rotation axis to K axis,
     - phi gives anticlockwise rotation from I to projection of
     rotation axis onto I-J plane,
     - kappa is the rotation about the rotation axis. */
  class CLIPPER_IMEX Polar_ccp4 {
  public:
    //! null constructor
    Polar_ccp4() {}
    //! constructor: from specified angles
    Polar_ccp4( const ftype& omega, const ftype& phi,  const ftype& kappa ) :
      omega_(omega), phi_(phi), kappa_(kappa) {}
    const ftype& psi() const { return omega_; }    //!< return omega
    const ftype& omega() const { return omega_; }  //!< return omega
    const ftype& phi()   const { return phi_;  }   //!< return phi
    const ftype& kappa() const { return kappa_; }  //!< return kappa
    String format() const;  //!< return formatted String representation
  private:
    ftype omega_, phi_, kappa_;
  };

  //! Rotation class
  /*! This class represents a rotation. The internal representation is
    as a unit quaternion, which is easily combined, inverted, or
    converted to or from other commonly used forms. */
  class CLIPPER_IMEX Rotation {
  public:
    //! null constructor
    Rotation() {}
    //! constructor: from generic Euler
    template<int T> explicit Rotation( const Euler<T>& euler )
      { (*this) = euler.rotation(); }
    //! constructor: from Euler_ccp4
    explicit Rotation( const Euler_ccp4& euler );
    //! constructor: from Polar_ccp4
    explicit Rotation( const Polar_ccp4& polar );
    //! constructor: from Matrix
    explicit Rotation( const Mat33<>& matrix );
    //! constructor: from components
    Rotation( const ftype& w, const ftype& x, const ftype& y,  const ftype& z )
      : w_(w), x_(x), y_(y), z_(z) {}
    const ftype& w() const { return w_; }  //!< return w component
    const ftype& x() const { return x_; }  //!< return x component
    const ftype& y() const { return y_; }  //!< return y component
    const ftype& z() const { return z_; }  //!< return z component
    template<int T> Euler<T> euler() const //!< return Euler angles
      { return Euler<T>( *this ); }
    Euler_ccp4 euler_ccp4() const;  //!< return Euler_ccp4 angles
    Polar_ccp4 polar_ccp4() const;  //!< return Polar_ccp4 angles
    Mat33<> matrix() const;    //!< return 3x3 matrix
    //! normalise this quaternion
    const Rotation& norm();
    //! return absolute rotation angle
    ftype abs_angle() const;
    //! return inverse rotation
    Rotation inverse() const { return Rotation( w_, -x_, -y_, -z_ ); }
    //! return zero rotation
    static Rotation zero() { return Rotation( 1.0, 0.0, 0.0, 0.0 ); }
    //! return null rotation
    static Rotation null() { return Rotation( Util::nan(), 0.0, 0.0, 0.0 ); }
    //! test for null (uninitialised) rotation
    bool is_null() const { return Util::is_nan( w_ ); }
    //! combine two rotations
    friend Rotation operator* ( const Rotation& r1, const Rotation& r2 );
    String format() const;  //!< return formatted String representation
    //! Enumeration of Euler conventions
    enum EULERtype { EulerXYZr,EulerXYZs,EulerXYXr,EulerXYXs,
		     EulerXZXr,EulerXZXs,EulerXZYr,EulerXZYs,
		     EulerYZXr,EulerYZXs,EulerYZYr,EulerYZYs,
		     EulerYXYr,EulerYXYs,EulerYXZr,EulerYXZs,
		     EulerZXYr,EulerZXYs,EulerZXZr,EulerZXZs,
		     EulerZYZr,EulerZYZs,EulerZYXr,EulerZYXs };
  protected:
    ftype w_, x_, y_, z_;
  };


} // namespace clipper

#endif
