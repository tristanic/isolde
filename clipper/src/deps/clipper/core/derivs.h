/*! \file lib/derivs.h
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


#ifndef CLIPPER_DERIVS
#define CLIPPER_DERIVS


#include "coords.h"

namespace clipper
{
  template<class T> class Grad_orth;
  template<class T> class Grad_frac;
  template<class T> class Grad_map;
  template<class T> class Curv_orth;
  template<class T> class Curv_frac;
  template<class T> class Curv_map;


  //! orthogonal (Angstom) gradient, with respect to orthogonal x,y,z
  template<class T> class Grad_orth : public Vec3<T>
  {
  public:
    Grad_orth() {}                 //!< null constructor
    explicit Grad_orth( const Vec3<T>& v ) :
      Vec3<T>( v ) {}          //!< constructor: copy/convert
    Grad_orth( const T& dx, const T& dy, const T& dz ) :
      Vec3<T>( dx, dy, dz ) {} //!< constructor: from d/dx,d/dy,d/dz
    const T& dx() const { return (*this)[0]; }  //!< get d/dx
    const T& dy() const { return (*this)[1]; }  //!< get d/dy
    const T& dz() const { return (*this)[2]; }  //!< get d/dz
    //! orthogonal-fractional derivative conversion
    Grad_frac<T> grad_frac( const Cell& cell ) const;
    String format() const;  //!< return formatted String representation
  };


  //! fractional (cell) gradient, with respect to fractional u,v,w
  template<class T> class Grad_frac : public Vec3<T>
  {
  public:
    Grad_frac() {}                 //!< null constructor
    explicit Grad_frac( const Vec3<T>& v ) :
      Vec3<T>( v ) {}          //!< constructor: copy/convert
    Grad_frac( const T& du, const T& dv, const T& dw ) :
      Vec3<T>( du, dv, dw ) {} //!< constructor: from d/du,d/dv,d/dw
    const T& du() const { return (*this)[0]; }  //!< get d/du
    const T& dv() const { return (*this)[1]; }  //!< get d/dv
    const T& dw() const { return (*this)[2]; }  //!< get d/dw
    //! fractional-orthogonal derivative conversion
    Grad_orth<T> grad_orth( const Cell& cell ) const;
    //! fractional-grid derivative conversion
    Grad_map<T> grad_map( const Grid& g ) const;
    String format() const;  //!< return formatted String representation
  };


  //! map coordinate gradient, with respect to grid u,v,w
  template<class T> class Grad_map : public Vec3<T>
  {
  public:
    Grad_map() {}                 //!< null constructor
    explicit Grad_map( const Vec3<T>& v ) :
      Vec3<T>( v ) {}          //!< constructor: copy/convert
    Grad_map( const T& du, const T& dv, const T& dw ) :
      Vec3<T>( du, dv, dw ) {} //!< constructor: from d/du,d/dv,d/dw
    const T& du() const { return (*this)[0]; }  //!< get d/du
    const T& dv() const { return (*this)[1]; }  //!< get d/dv
    const T& dw() const { return (*this)[2]; }  //!< get d/dw
    //! grid-fractional derivative conversion
    Grad_frac<T> grad_frac( const Grid& g ) const;
    String format() const;  //!< return formatted String representation
  };


  //! orthogonal (Angstom) curvatures, with respect to orthogonal x,y,z
  template<class T> class Curv_orth : public Mat33<T>
  {
  public:
    Curv_orth() {}                 //!< null constructor
    explicit Curv_orth( const Mat33<T>& m ) :
      Mat33<T>( m ) {}         //!< constructor: copy/convert
    //! orthogonal-fractional derivative conversion
    Curv_frac<T> curv_frac( const Cell& cell ) const;
  };


  //! fractional (cell) curvatures, with respect to fractional u,v,w
  template<class T> class Curv_frac : public Mat33<T>
  {
  public:
    Curv_frac() {}                 //!< null constructor
    explicit Curv_frac( const Mat33<T>& m ) :
      Mat33<T>( m ) {}         //!< constructor: copy/convert
    //! fractional-orthogonal derivative conversion
    Curv_orth<T> curv_orth( const Cell& cell ) const;
    //! fractional-grid derivative conversion
    Curv_map<T> curv_map( const Grid& g ) const;
  };


  //! map coordinate curvatures, with respect to grid u,v,w
  template<class T> class Curv_map : public Mat33<T>
  {
  public:
    Curv_map() {}                 //!< null constructor
    explicit Curv_map( const Mat33<T>& m ) :
      Mat33<T>( m ) {}         //!< constructor: copy/convert
    //! grid-fractional derivative conversion
    Curv_frac<T> curv_frac( const Grid& g ) const;
  };



  // template implementations

  /*! The result is an RT operator. This is a redudent representation,
    but is handy for assembling compound operators.
    \return The operator */
  /*! \return The formatted text string */
  template<class T> String Grad_orth<T>::format() const
    { return "d/dx,d/dy,d/dz = ("+String(dx())+","+String(dy())+","+String(dz())+")"; }

  /*! \param cell The cell concerned \return The transformed derivative. */
  template<class T> inline Grad_frac<T> Grad_orth<T>::grad_frac( const Cell& cell ) const
    { return Grad_frac<T>( (*this) * Mat33<T>( cell.matrix_orth() ) ); }


  /*! \return The formatted text string */
  template<class T> String Grad_frac<T>::format() const
    { return "d/du,d/dv,d/dw = ("+String(du())+","+String(dv())+","+String(dw())+")"; }

  /*! \param cell The cell concerned \return The transformed derivative. */
  template<class T> inline Grad_orth<T> Grad_frac<T>::grad_orth( const Cell& cell ) const
    { return Grad_orth<T>( (*this) * Mat33<T>( cell.matrix_frac() ) ); }

  /*! \param g The grid concerned \return The transformed derivative. */
  template<class T> inline Grad_map<T> Grad_frac<T>::grad_map( const Grid& g ) const
    { return Grad_map<T>( du()/g.nu(), dv()/g.nv(), dw()/g.nw() ); }


  /*! \return The formatted text string */
  template<class T> String Grad_map<T>::format() const
    { return "d/du,d/dv,d/dw = ("+String(du())+","+String(dv())+","+String(dw())+")"; }

  /*! \param g The grid concerned \return The transformed derivative. */
  template<class T> inline Grad_frac<T> Grad_map<T>::grad_frac( const Grid& g ) const
    { return Grad_frac<T>( du()*g.nu(), dv()*g.nv(), dw()*g.nw() ); }


  /*! \param cell The cell concerned \return The transformed derivative. */
  template<class T> Curv_frac<T> Curv_orth<T>::curv_frac( const Cell& cell ) const
  {
    Mat33<T> m( cell.matrix_orth() );
    return Curv_frac<T>( m.transpose() * (*this) * m );
  }


  /*! \param cell The cell concerned \return The transformed derivative. */
  template<class T> Curv_orth<T> Curv_frac<T>::curv_orth( const Cell& cell ) const
  {
    Mat33<T> m( cell.matrix_frac() );
    return Curv_orth<T>( m.transpose() * (*this) * m );
  }

  /*! \param g The grid concerned \return The transformed derivative. */
  template<class T> Curv_map<T> Curv_frac<T>::curv_map( const Grid& g ) const
  {
    Curv_map<T> c;
    c(0,0) = (*this)(0,0) / T(g.nu()*g.nu());
    c(0,1) = (*this)(0,1) / T(g.nu()*g.nv());
    c(0,2) = (*this)(0,2) / T(g.nu()*g.nw());
    c(1,0) = (*this)(1,0) / T(g.nv()*g.nu());
    c(1,1) = (*this)(1,1) / T(g.nv()*g.nv());
    c(1,2) = (*this)(1,2) / T(g.nv()*g.nw());
    c(2,0) = (*this)(2,0) / T(g.nw()*g.nu());
    c(2,1) = (*this)(2,1) / T(g.nw()*g.nv());
    c(2,2) = (*this)(2,2) / T(g.nw()*g.nw());
    return c;
  }


  /*! \param g The grid concerned \return The transformed derivative. */
  template<class T> Curv_frac<T> Curv_map<T>::curv_frac( const Grid& g ) const
  {
    Curv_frac<T> c;
    c(0,0) = (*this)(0,0) * T(g.nu()*g.nu());
    c(0,1) = (*this)(0,1) * T(g.nu()*g.nv());
    c(0,2) = (*this)(0,2) * T(g.nu()*g.nw());
    c(1,0) = (*this)(1,0) * T(g.nv()*g.nu());
    c(1,1) = (*this)(1,1) * T(g.nv()*g.nv());
    c(1,2) = (*this)(1,2) * T(g.nv()*g.nw());
    c(2,0) = (*this)(2,0) * T(g.nw()*g.nu());
    c(2,1) = (*this)(2,1) * T(g.nw()*g.nv());
    c(2,2) = (*this)(2,2) * T(g.nw()*g.nw());
    return c;
  }


} // namespace clipper

#endif
