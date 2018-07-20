/*! \file lib/cell.h
    Header file for unit cell class
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
/* This code is derived from the 'dm' source code */


#ifndef CLIPPER_CELL
#define CLIPPER_CELL


#include "clipper_types.h"
#include "../imex.h"

namespace clipper
{


  //! Metric tensor
  /*! The metric tensor is used to determine a distance in real or
    reciprocal space using fraction coordinates or Miller indices. It
    is symmetrical, so only the upper triangle is stored with the
    off-diagonal elements doubled.
  */
  class CLIPPER_IMEX Metric_tensor
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
  private:
    ftype m00, m11, m22, m01, m02, m12;
  };


  //! cell description (automatically converts to radians)
  /*! The cell description is a compact description of a cell,
    containing just the cell parameters. It is usually used to
    construct a full Cell object, which provides the expected
    functionality.
  */
  class CLIPPER_IMEX Cell_descr
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

  protected:
    ftype a_,b_,c_,alpha_,beta_,gamma_;
  };


  //! Cell object
  /*! The Cell class is the fully functional description of the unit
    cell. In addition to the cell parameters, it stores derived
    information including the cell volume, orthogonalising and
    fractionalising matrices, and the metric tensors.
   */
  class CLIPPER_IMEX Cell : public Cell_descr
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
    inline const Cell_descr& descr() const { return (*this); }
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

  private:
    Cell_descr descr_;          //!< unit cell parameters
    ftype vol;                  //!< unit cell volume
    Mat33<> orthmat;       //!< orthogonalisation matrix
    Mat33<> fracmat;       //!< fractionalisation matrix
    Metric_tensor realmetric;   //!< real space metric tensor
    Metric_tensor recimetric;   //!< reciprocal space metric tensor
  };


} // namespace clipper

#endif
