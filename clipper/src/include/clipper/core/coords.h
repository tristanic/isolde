/*! \file lib/coords.h
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


#ifndef CLIPPER_COORDS
#define CLIPPER_COORDS


#include "cell.h"
#include "spacegroup.h"
#include "clipper_stats.h"
#include "../imex.h"

namespace clipper
{
  // forward definitions
  class Grid; class Grid_sampling; class Grid_range;
  class Coord_grid; class Coord_map;
  class Coord_reci_frac; class Coord_reci_orth;
  class Coord_frac; class Coord_orth;
  class U_aniso_frac; class U_aniso_orth;


  //! Resolution in angstroms
  /*! This object represents a resolution limit which will be used for
    all aspects of a calculation.
    This is a base for a donor type. */
  class CLIPPER_IMEX Resolution
  {
  public:
    inline Resolution() : resol(0.0) {}          //!< null constructor
    explicit Resolution( const ftype& resol_ );  //!< constructor: from ftype
    void init( const ftype& resol_ );            //!< initialiser: from ftype
    const ftype& limit() const;                  //!< get resolution limit
    ftype invresolsq_limit() const;  //!< get invresolsq limit
    bool is_null() const;            //!< test if value has been initialised
  private:
    ftype resol;
  };


  //! reflection class
  /*! This describes the type of a reflection in a given spacegroup,
    including centricity, systematic absence, phase restriction, and
    multiplicity. */
  class CLIPPER_IMEX HKL_class
  {
  public:
    //! null constructor
    inline HKL_class() { epsilon_ = 0; allowed_ = 255; }
    //! constructor - from spacegroup and HKL
    HKL_class( const Spacegroup& spgr, const HKL& hkl );
    //! get epsilon
    inline ftype epsilon() const { return ftype(epsilon_); }
    //! get epsilon for acentric, 2x epsilon for centric
    inline ftype epsilonc() const
      { if ( centric() ) return 2.0*ftype(epsilon_);
        else             return ftype(epsilon_); }
    //! get allowed phase
    inline ftype allowed() const { return ftype(allowed_) * (Util::pi()/12.0); }
    inline bool centric() const { return allowed_ != 255; } //!< is centric?
    inline bool sys_abs() const { return epsilon_ == 0; }   //!< is sys abs?
  private:
    unsigned char epsilon_, allowed_;
  };


  //! Orthogonal operator class.
  /*! This class is used for any RT-operator which operates on
    orthogonal coordinates. For a full list of methods, see
    clipper::RTop */
  class CLIPPER_IMEX RTop_orth : public RTop<>
  {
  public:
    //! null constructor
    inline RTop_orth() {}
    //! constructor: copy/convert
    inline explicit RTop_orth( const RTop<>& o ) : RTop<>( o ) {}
    //! constructor: from rotation
    inline explicit RTop_orth( const Mat33<>& r ) : RTop<>( r ) {}
    //! constructor: from rotation and translation
    inline RTop_orth( const Mat33<>& r, const Vec3<>& t ) : RTop<>( r, t ) {}
    //! constructor: from two vectors of Coord_orth
    RTop_orth( const std::vector<Coord_orth>& src, const std::vector<Coord_orth>& tgt );
    //! constructor: from two vectors of Coord_orth
    RTop_orth( const std::vector<Coord_orth>& src, const std::vector<Coord_orth>& tgt, const std::vector<ftype>& wgt );
    //! constructor: from two atom-list type objects
    template<class T> RTop_orth( const T& src, const T& tgt );
    //! orthogonal-fractional conversion
    RTop_frac rtop_frac( const Cell& cell ) const;
    //! inverse operator
    RTop_orth inverse() const;
    //! return point on axis near the specified coordinate
    Coord_orth axis_coordinate_near( const Coord_orth& centre ) const;
    //! return screw translation
    Coord_orth screw_translation() const;
    //! return identity operator
    static RTop_orth identity();
    //! return null (uninitialised) operator
    static RTop_orth null();
  };


  //! reflection 'Miller' index
  class CLIPPER_IMEX HKL : public Vec3<int>
  {
  public:
    inline HKL() {}                //!< null constructor
    inline explicit HKL( const Vec3<int>& v ) :
      Vec3<int>( v ) {}            //!< constructor: copy/convert
    inline HKL( const int& h, const int& k, const int& l ) :
      Vec3<int>( h, k, l ) {}      //!< constructor: from H,K,L
    inline const int& h() const { return (*this)[0]; }  //!< get h
    inline const int& k() const { return (*this)[1]; }  //!< get k
    inline const int& l() const { return (*this)[2]; }  //!< get l
    inline int& h() { return (*this)[0]; }  //!< set h
    inline int& k() { return (*this)[1]; }  //!< set k
    inline int& l() { return (*this)[2]; }  //!< set l
    //! return inverse resolution squared for this reflection in given cell
    inline ftype invresolsq( const Cell& cell ) const;
    //! return fractional reciprocal coordinate (i.e. non-integer HKL)
    inline Coord_reci_frac coord_reci_frac() const;
    //! orthogonal-fractional reciprocal space coordinate conversion
    inline Coord_reci_orth coord_reci_orth( const Cell& cell ) const;
    //! return transformed hkl
    inline HKL transform( const Symop& op ) const;
    //! return transformed hkl
    inline HKL transform( const Isymop& op ) const;
    //! return symmetry phase shift for this HKL under op
    inline ftype sym_phase_shift( const Symop& op ) const;
    String format() const;  //!< return formatted String representation
    friend inline HKL operator -(const HKL& h1)
      { return HKL( -h1.h(), -h1.k(), -h1.l() ); }
    friend inline HKL operator +(const HKL& h1, const HKL& h2)
      { return HKL( h1.h()+h2.h(), h1.k()+h2.k(), h1.l()+h2.l() ); }
    friend inline HKL operator -(const HKL& h1, const HKL& h2)
      { return HKL( h1.h()-h2.h(), h1.k()-h2.k(), h1.l()-h2.l() ); }
    friend inline HKL operator *(const int& s, const HKL& h1)
      { return HKL( s*h1.h(), s*h1.k(), s*h1.l() ); }
    friend inline HKL operator *(const Isymop& op, const HKL& h1)
      { return HKL( h1 * op.rot() ); }
  };


  //! orthogonal reciprocal coordinate (length of which is invresolsq)
  class CLIPPER_IMEX Coord_reci_orth : public Vec3<>
  {
  public:
    inline Coord_reci_orth() {}           //!< null constructor
    inline explicit Coord_reci_orth( const Vec3<>& v ) :
      Vec3<>( v ) {}          //!< constructor: copy/convert
    inline Coord_reci_orth( const ftype& xs, const ftype& ys, const ftype& zs ) : Vec3<>( xs, ys, zs ) {} //!< constructor: from x*,y*,z*
    inline const ftype& xs() const { return (*this)[0]; }  //!< get x*
    inline const ftype& ys() const { return (*this)[1]; }  //!< get y*
    inline const ftype& zs() const { return (*this)[2]; }  //!< get z*
    //! return inverse resolution squared for this coord
    inline ftype invresolsq() const;
    //! orthogonal-fractional reciprocal space coordinate conversion
    inline Coord_reci_frac coord_reci_frac( const Cell& cell ) const;
    //! return transformed coordinate
    inline Coord_reci_orth transform( const RTop_orth& op ) const
      { return Coord_reci_orth( (*this) * op.rot() ); }
    String format() const;  //!< return formatted String representation
  };


  //! fractional reciprocal coordinate (i.e. non-integer hkl)
  class CLIPPER_IMEX Coord_reci_frac : public Vec3<>
  {
  public:
    inline Coord_reci_frac() {}           //!< null constructor
    inline explicit Coord_reci_frac( const Vec3<>& v ) :
      Vec3<>( v ) {}          //!< constructor: copy/convert
    inline Coord_reci_frac( const ftype& us, const ftype& vs, const ftype& ws ) : Vec3<>( us, vs, ws ) {} //!< constructor: from u,v,w
    //! constructor: from HKL
    inline Coord_reci_frac( const HKL& hkl ) :
      Vec3<>( ftype(hkl[0]), ftype(hkl[1]), ftype(hkl[2]) ) {}
    //! round to HKL
    inline HKL hkl() const
      { return HKL( Util::intr(us()), Util::intr(vs()), Util::intr(ws()) ); }
    //! return inverse resolution squared for this reflection in given cell
    inline ftype invresolsq( const Cell& cell ) const;
    inline const ftype& us() const { return (*this)[0]; }  //!< get u*
    inline const ftype& vs() const { return (*this)[1]; }  //!< get v*
    inline const ftype& ws() const { return (*this)[2]; }  //!< get w*
    //! fractional-orthogonal reciprocal space coordinate conversion
    inline Coord_reci_orth coord_reci_orth( const Cell& cell ) const;
    //! return transformed coordinate
    inline Coord_reci_frac transform( const RTop_frac& op ) const
      { return Coord_reci_frac( (*this) * op.rot() ); }
    String format() const;  //!< return formatted String representation
  };


  //! Grid coordinate
  class CLIPPER_IMEX Coord_grid : public Vec3<int>
  {
  public:
    inline Coord_grid() {}  //!< null constructor
    //! constructor: copy/convert
    inline explicit Coord_grid( const Vec3<int> v ) : Vec3<int>( v ) {}
    //! constructor: from u,v,w
    inline Coord_grid( const int& u, const int& v, const int& w ) :
      Vec3<int>(u,v,w) {}
    //! constructor: from a grid and an index in that grid
    inline Coord_grid( const Grid& g, const int& index )
      { deindex( g, index ); }
    inline const int& u() const { return (*this)[0]; }  //!< get u
    inline const int& v() const { return (*this)[1]; }  //!< get v
    inline const int& w() const { return (*this)[2]; }  //!< get w
    inline int& u() { return (*this)[0]; }  //!< set u
    inline int& v() { return (*this)[1]; }  //!< set v
    inline int& w() { return (*this)[2]; }  //!< set w

    //! convert to Coord_map
    inline Coord_map coord_map() const;
    //! convert to Coord_frac using given Grid_sampling
    inline Coord_frac coord_frac( const Grid_sampling& g ) const;
    //! return transformed coordinate
    inline Coord_grid transform( const Isymop& op ) const
      { return op * (*this); }

    //! reduce to unit box: (0..nu-1, 0..nv-1, 0..nw-1)
    inline Coord_grid unit( const Grid_sampling& g ) const;

    //! increment in storage order (see index())
    /*! guaranteed to increment index(g) by 1 */
    inline const Coord_grid& next( const Grid& g );
    //! increment in storage order (see index())
    /*! guaranteed to increment index(g) by 1 */
    inline const Coord_grid& next( const Grid_range& g );
    //! test if done in storage order (see index())
    inline bool last( const Grid& g ) const;
    //! test if done in storage order (see index())
    inline bool last( const Grid_range& g ) const;
    //! grid indexing operator
    inline int index( const Grid& g ) const;
    //! grid deindexing operator
    inline void deindex( const Grid& g, const int& index );
    // grid indexing operator
    //inline int index( const Grid_range& g ) const;
    // grid deindexing operator
    //inline void deindex( const Grid_range& g, const int& index );

    String format() const;  //!< return formatted String representation
    friend inline Coord_grid operator -(const Coord_grid& r1)
      { return ( Coord_grid( -r1.u(), -r1.v(), -r1.w() ) ); }
    friend inline Coord_grid operator +(const Coord_grid& r1, const Coord_grid& r2) { return ( Coord_grid( r1.u()+r2.u(), r1.v()+r2.v(), r1.w()+r2.w() ) ); }
    friend inline Coord_grid operator -(const Coord_grid& r1, const Coord_grid& r2) { return ( Coord_grid( r1.u()-r2.u(), r1.v()-r2.v(), r1.w()-r2.w() ) ); }
    friend inline Coord_grid operator *(const int& s, const Coord_grid& r1)
      { return ( Coord_grid( s*r1.u(), s*r1.v(), s*r1.w() ) ); }
    friend inline int operator == (const Coord_grid& r1,  const Coord_grid& r2)
      { return (r1.u()==r2.u() && r1.v()==r2.v() && r1.w()==r2.w()); }
    friend inline int operator != (const Coord_grid& r1,  const Coord_grid& r2)
      { return (r1.u()!=r2.u() || r1.v()!=r2.v() || r1.w()!=r2.w()); }
    friend inline Coord_grid operator *(const Isymop& op, const Coord_grid& r1)
      { return Coord_grid( op.rot() * r1 + op.trn() ); }
  };


  //! orthogonal (Angstrom) coordinates
  class CLIPPER_IMEX Coord_orth : public Vec3<>
  {
  public:
    inline Coord_orth() {}                //!< null constructor
    inline explicit Coord_orth( const Vec3<>& v ) :
      Vec3<>( v ) {}          //!< constructor: copy/convert
    inline Coord_orth( const ftype& x, const ftype& y, const ftype& z ) :
      Vec3<>( x, y, z ) {}    //!< constructor: from x,y,z
    //! constructor: from 3 coords and bond length, angle, torsion
    Coord_orth( const Coord_orth& x1, const Coord_orth& x2, const Coord_orth& x3, const ftype& length, const ftype& angle, const ftype& torsion );
    inline const ftype& x() const { return (*this)[0]; }  //!< get x
    inline const ftype& y() const { return (*this)[1]; }  //!< get y
    inline const ftype& z() const { return (*this)[2]; }  //!< get z
    //! return square of length of vector in Angstroms
    inline ftype lengthsq() const;
    //! orthogonal-fractional coordinate conversion
    inline Coord_frac coord_frac( const Cell& cell ) const;
    //! return transformed coordinate
    inline Coord_orth transform( const RTop_orth& op ) const
      { return op*(*this); }
    String format() const;  //!< return formatted String representation
    //! Return length of vector between two coord orths
    static ftype length( const Coord_orth& x1, const Coord_orth& x2);
    //! Return angle between three coord orths
    static ftype angle( const Coord_orth& x1, const Coord_orth& x2, 
		        const Coord_orth& x3);
    //! Return torsion between four coord orths
    static ftype torsion( const Coord_orth& x1, const Coord_orth& x2,
			  const Coord_orth& x3, const Coord_orth& x4);
    friend inline Coord_orth operator -(const Coord_orth& x1)
      { return Coord_orth( -x1.x(), -x1.y(), -x1.z() ); }
    friend inline Coord_orth operator +(const Coord_orth& x1, const Coord_orth& x2) { return Coord_orth( x1.x()+x2.x(), x1.y()+x2.y(), x1.z()+x2.z() ); }
    friend inline Coord_orth operator -(const Coord_orth& x1, const Coord_orth& x2) { return Coord_orth( x1.x()-x2.x(), x1.y()-x2.y(), x1.z()-x2.z() ); }
    friend inline Coord_orth operator *(const ftype& s, const Coord_orth& x1)
      { return Coord_orth( s*x1.x(), s*x1.y(), s*x1.z() ); }
    friend inline Coord_orth operator *(const RTop_orth& op, const Coord_orth& x1) { return Coord_orth( op.rot() * x1 + op.trn() ); }
  };


  //! fractional (cell) coordinates
  class CLIPPER_IMEX Coord_frac : public Vec3<>
  {
  public:
    inline Coord_frac() {}                //!< null constructor
    inline explicit Coord_frac( const Vec3<>& v ) :
      Vec3<>( v ) {}          //!< constructor: copy/convert
    inline Coord_frac( const ftype& u, const ftype& v, const ftype& w ) :
      Vec3<>( u, v, w ) {}    //!< constructor: from u,v,w
    inline const ftype& u() const { return (*this)[0]; }  //!< get u
    inline const ftype& v() const { return (*this)[1]; }  //!< get v
    inline const ftype& w() const { return (*this)[2]; }  //!< get w
    //! return square of length of vector in Angstroms
    inline ftype lengthsq( const Cell& cell ) const;
    //! fractional-orthogonal coordinate conversion
    inline Coord_orth coord_orth( const Cell& cell ) const;
    //! fractional-grid coordinate conversion
    inline Coord_map coord_map( const Grid& g ) const;
    //! fractional-grid coordinate conversion
    inline Coord_grid coord_grid( const Grid& g ) const;
    //! return transformed coordinate
    inline Coord_frac transform( const RTop_frac& op ) const
      { return op*(*this); }
    //! return lattice copy nearest origin
    inline Coord_frac lattice_copy_zero() const
      { return Coord_frac(u()-rint(u()),v()-rint(v()),w()-rint(w())); }
    //! return lattice copy in unit box (0...1,0...1,0...1)
    inline Coord_frac lattice_copy_unit() const
      { return Coord_frac(u()-floor(u()),v()-floor(v()),w()-floor(w())); }
    //! return lattice copy near the specified coordinate
    inline Coord_frac lattice_copy_near(const Coord_frac& n) const
      { return (*this-n).lattice_copy_zero()+n; }
    //! return symmetry copy near the specified coordinate
    Coord_frac symmetry_copy_near(const Spacegroup& spgr, const Cell& cell, const Coord_frac& n) const;
    String format() const;  //!< return formatted String representation
    friend inline Coord_frac operator -(const Coord_frac& u1)
      { return Coord_frac( -u1.u(), -u1.v(), -u1.w() ); }
    friend inline Coord_frac operator +(const Coord_frac& u1, const Coord_frac& u2) { return Coord_frac( u1.u()+u2.u(), u1.v()+u2.v(), u1.w()+u2.w() ); }
    friend inline Coord_frac operator -(const Coord_frac& u1, const Coord_frac& u2) { return Coord_frac( u1.u()-u2.u(), u1.v()-u2.v(), u1.w()-u2.w() ); }
    friend inline Coord_frac operator *(const ftype& s, const Coord_frac& u1)
      { return Coord_frac( s*u1.u(), s*u1.v(), s*u1.w() ); }
    friend inline Coord_frac operator *(const RTop_frac& op, const Coord_frac& x1) { return Coord_frac( op.rot() * x1 + op.trn() ); }
  };


  //! map coordinate: this is like Coord_grid, but non-integer
  class CLIPPER_IMEX Coord_map : public Vec3<>
  {
  public:
    inline Coord_map() {}  //!< null constructor
    //! constructor: copy/convert
    inline explicit Coord_map( const Vec3<>& v ) :
      Vec3<>( v ) {}
    //! constructor: from Coord_grid
    inline explicit Coord_map( const Coord_grid& c ) :
      Vec3<>( ftype(c[0]), ftype(c[1]), ftype(c[2]) ) {}
    //! constructor: from u,v,w
    inline Coord_map( const ftype& u, const ftype& v, const ftype& w ) :
      Vec3<>( u, v, w ) {}
    //! grid-fractional coordinate conversion
    inline Coord_frac coord_frac( const Grid& g ) const;
    //! return integer Coord_grid nearest this coordinate
    inline Coord_grid coord_grid() const { return Coord_grid( Util::intr((*this)[0]), Util::intr((*this)[1]), Util::intr((*this)[2]) ); }
    //! return integer Coord_grid below this coordinate
    inline Coord_grid floor() const { return Coord_grid( Util::intf((*this)[0]), Util::intf((*this)[1]), Util::intf((*this)[2]) ); }
    //! return integer Coord_grid above this coordinate
    inline Coord_grid ceil() const { return Coord_grid( Util::intc((*this)[0]), Util::intc((*this)[1]), Util::intc((*this)[2]) ); }
    inline const ftype& u() const { return (*this)[0]; }  //!< get u
    inline const ftype& v() const { return (*this)[1]; }  //!< get v
    inline const ftype& w() const { return (*this)[2]; }  //!< get w
    String format() const;  //!< return formatted String representation
    friend inline Coord_map operator -(const Coord_map& u1)
      { return Coord_map( -u1.u(), -u1.v(), -u1.w() ); }
    friend inline Coord_map operator +(const Coord_map& u1, const Coord_map& u2)
      { return Coord_map( u1.u()+u2.u(), u1.v()+u2.v(), u1.w()+u2.w() ); }
    friend inline Coord_map operator -(const Coord_map& u1, const Coord_map& u2)
      { return Coord_map( u1.u()-u2.u(), u1.v()-u2.v(), u1.w()-u2.w() ); }
    friend inline Coord_map operator *(const ftype& s, const Coord_map& u1)
      { return Coord_map( s*u1.u(), s*u1.v(), s*u1.w() ); }
  };


  //! Anisotropic orthogonal atomic displacement parameters
  /*! These are defined on orthogonal atomic coordinates in
    A<sup>-2</sup>, i.e. they are anisotropic U values. */
  class CLIPPER_IMEX U_aniso_orth : public Mat33sym<>
  {
  public:
    //! null constructor
    inline U_aniso_orth() {};
    //! constructor: from Mat33sym
    inline explicit U_aniso_orth( const Mat33sym<>& m ) : Mat33sym<>(m) {}
    //! constructor: from isotropic U
    inline explicit U_aniso_orth( const ftype& u ) :
      Mat33sym<>( u, u, u, 0.0, 0.0, 0.0 ) {}
    //! constructor: from Uij
    U_aniso_orth( const ftype& u11, const ftype& u22, const ftype& u33,
		  const ftype& u12, const ftype& u13, const ftype& u23 ) :
      Mat33sym<>( u11, u22, u33, u12, u13, u23 ) {}
    //! return nearest isotropic U
    ftype u_iso() const;
    //! orthogonal-fractional conversion
    U_aniso_frac u_aniso_frac( const Cell& cell ) const;
    //! return transformed U_aniso
    U_aniso_orth transform( const RTop_orth& op ) const;
    friend U_aniso_orth operator +(const U_aniso_orth& u1, const U_aniso_orth& u2) { return U_aniso_orth( u1.mat00()+u2.mat00(), u1.mat11()+u2.mat11(), u1.mat22()+u2.mat22(), u1.mat01()+u2.mat01(), u1.mat02()+u2.mat02(), u1.mat12()+u2.mat12() ); }
    friend U_aniso_orth operator -(const U_aniso_orth& u) { return U_aniso_orth( -u.mat00(), -u.mat11(), -u.mat22(), -u.mat01(), -u.mat02(), -u.mat12() ); }
    friend U_aniso_orth operator *(const ftype& s, const U_aniso_orth& u) { return U_aniso_orth( s*u.mat00(), s*u.mat11(), s*u.mat22(), s*u.mat01(), s*u.mat02(), s*u.mat12() ); }
  };


  //! Anisotropic fractional atomic displacement parameters
  /*! These are defined on fractional atomic coordinates in
    A<sup>-2</sup>, i.e. they are anisotropic U values. */
  class CLIPPER_IMEX U_aniso_frac : public Mat33sym<>
  {
  public:
    //! null constructor
    inline U_aniso_frac() {};
    //! constructor: from Mat33sym
    inline explicit U_aniso_frac( const Mat33sym<>& m ) : Mat33sym<>(m) {}
    //! constructor: from Uij
    U_aniso_frac( const ftype& u11, const ftype& u22, const ftype& u33,
		  const ftype& u12, const ftype& u13, const ftype& u23 ) :
      Mat33sym<>( u11, u22, u33, u12, u13, u23 ) {}
    //! fractional-orthogonal conversion
    U_aniso_orth u_aniso_orth( const Cell& cell ) const;
    //! return transformed U_aniso
    U_aniso_frac transform( const RTop_frac& op ) const;
    friend U_aniso_frac operator +(const U_aniso_frac& u1, const U_aniso_frac& u2) { return U_aniso_frac( u1.mat00()+u2.mat00(), u1.mat11()+u2.mat11(), u1.mat22()+u2.mat22(), u1.mat01()+u2.mat01(), u1.mat02()+u2.mat02(), u1.mat12()+u2.mat12() ); }
    friend U_aniso_frac operator -(const U_aniso_frac& u) { return U_aniso_frac( -u.mat00(), -u.mat11(), -u.mat22(), -u.mat01(), -u.mat02(), -u.mat12() ); }
    friend U_aniso_frac operator *(const ftype& s, const U_aniso_frac& u) { return U_aniso_frac( s*u.mat00(), s*u.mat11(), s*u.mat22(), s*u.mat01(), s*u.mat02(), s*u.mat12() ); }
  };


  //! generic grid
  /*! This holds the dimensions of a 3D array, indexed from 0 along
    each dimension. */ 
  class CLIPPER_IMEX Grid : public Vec3<int>
  {
  public:
    inline Grid() {}                     //!< null constructor
    inline Grid( const int& nu, const int& nv, const int& nw ) :
      Vec3<int>( nu, nv, nw ) {}  //!< constructor: from nu,nv,nw
    inline const int& nu() const { return (*this)[0]; } //!< get nu
    inline const int& nv() const { return (*this)[1]; } //!< get nv
    inline const int& nw() const { return (*this)[2]; } //!< get nw
    //! return size of grid array
    inline int size() const { return nu()*nv()*nw(); }
    //! determine if a point is in the grid
    inline bool in_grid( Coord_grid g ) const { return (g.u() >= 0 && g.u() < nu() && g.v() >= 0 && g.v() < nv() && g.w() >= 0 && g.w() < nw()); }

    //! grid indexing operator
    inline int index( const Coord_grid& c ) const { return c.index(*this); }
    //! grid deindexing operator
    inline Coord_grid deindex( const int& index ) const { return Coord_grid( *this, index ); }
    String format() const;  //!< return formatted String representation
    void debug() const;
  };


  //! Grid sampling of a unit cell
  /*! This class represents the grid sampling of a unit cell. It is
  otherwise identical to its parent, clipper::Grid_cell, but has an
  additional constructor which takes a spacegroup, cell and resolution
  and produces an appropriate grid obeying all of the symmetry
  constraints, and using efficient factors for the calculation of
  FFTs.

  \note The following methods are inherited from Grid and Grid_cell
  but are documented here for convenience: nu(), nv(), nw(), size(),
  index(), deindex(), format(), coord_frac(), coord_grid(),
  to_unit().
 */
  class CLIPPER_IMEX Grid_sampling : public Grid
  {
  public:
    //! null constructor
    inline Grid_sampling() : Grid(Grid(0,0,0)) {}
    //! constructor: from nu, nv, nw
    inline Grid_sampling( const int& nu, const int& nv, const int& nw ) :
      Grid( nu, nv, nw ) {}
    //! constructor: from Spacegroup, Cell, Resolution, Shannon rate
    Grid_sampling( const Spacegroup& spacegroup, const Cell& cell, const Resolution& resol, const ftype rate = 1.5 );
    //! initialiser: from Spacegroup, Cell, Resolution, Shannon rate
    void init( const Spacegroup& spacegroup, const Cell& cell, const Resolution& resol, const ftype rate = 1.5 );

    //! return matrix which converts grid to fractional coordinates
    Mat33<> matrix_grid_frac() const;
    //! return matrix which converts fractional to grid coordinates
    Mat33<> matrix_frac_grid() const;

    //! test if object has been initialised
    bool is_null() const;

    // inherited functions listed for documentation purposes
    //-- const int& nu() const;
    //-- const int& nv() const;
    //-- const int& nw() const;
    //-- int size() const;
    //-- int index( const Coord_grid& c ) const;
    //-- Coord_grid deindex( const int& index ) const;
    //-- const String format() const;
  };


  //! HKL sampling of reciprocal space
  /*! The HKL_sampling class uniquely describes a P0 reflection list
    bounded by some resolution limit in reciprocal space. It is
    described in terms of large integers, and so immune from rounding
    errors once the object is constructed. */
  class CLIPPER_IMEX HKL_sampling
  {
  public:
    //! null constructor
    HKL_sampling(); //!< Null constructor
    //! constructor: takes parameters of normal or inverse cell
    HKL_sampling( const Cell& cell, const Resolution& resolution );
    //! return limiting values of H, K, L
    HKL hkl_limit() const;
    //! return approximate resolution given cell
    Resolution resolution( const Cell& cell ) const;
    //! test if a reflection is within the resolution limit
    inline bool in_resolution( const HKL& h ) const
      {	return ( m00*itype64(h.h()*h.h()) + m11*itype64(h.k()*h.k()) +
		 m22*itype64(h.l()*h.l()) + m01*itype64(h.h()*h.k()) +
		 m02*itype64(h.h()*h.l()) + m12*itype64(h.k()*h.l()) )
	  <=   ( sqrt_limit_value*sqrt_limit_value ); }
    //! test if object has been initialised
    bool is_null() const;
    String format() const;  //!< return formatted String representation
    friend inline int operator == (const HKL_sampling& h1, const HKL_sampling& h2)
      { return ( h1.m00==h2.m00 && h1.m11==h2.m11 && h1.m22==h2.m22 &&
		 h1.m01==h2.m01 && h1.m02==h2.m02 && h1.m12==h2.m12 ); }
  private:
    static itype64 sqrt_limit_value;
    itype64 m00, m11, m22, m01, m02, m12;
  };


  //! Grid range class: defines array limits for a grid
  /*! This class is used for describing 3D grids covering an arbitrary
    part of the 3D space, i.e. which do not start from (0,0,0). */
  class CLIPPER_IMEX Grid_range : public Grid
  {
  public:
    //! null constructor
    inline Grid_range() {}
    //! constructor: takes grid limits
    Grid_range( const Coord_grid& min, const Coord_grid& max );
    //! constructor: takes cell grid and fractional limits
    Grid_range( const Grid& grid, const Coord_frac& min, const Coord_frac& max );
    //! constructor: make grid to hold a sphere from cell, grid, radius
    Grid_range( const Cell& cell, const Grid& grid, const ftype& radius );
    //! access grid limits
    const Coord_grid& min() const { return min_; }
    //! access grid limits
    const Coord_grid& max() const { return max_; }
    //! border: increase grid to include given border
    void add_border( const int b );
    //! determine if a point is in the grid
    bool in_grid( Coord_grid g ) const { return (g.u() >= min_.u() && g.u() <= max_.u() && g.v() >= min_.v() && g.v() <= max_.v() && g.w() >= min_.w() && g.w() <= max_.w()); }

    //! grid indexing operator
    int index( const Coord_grid& c ) const { return (c - min_).index(*this); } 
    //! grid deindexing operator
    Coord_grid deindex( const int& index ) const { return Coord_grid( *this, index ) + min_; }
  private:
    Coord_grid min_, max_;
  };
  //! Obsolete form for Grid_range
  typedef Grid_range Grid_map;


  //! Atom class
  /*! This class defines a minimal atom object providing only those
    properties required for an electron density calculation. A
    template constructor allows it to be constructed from any other
    object with appropriate properties. */
  class CLIPPER_IMEX Atom
  {
  public:
    //! null constructor
    Atom() {}
    //! Constructor: from atom-like object
    template<class T> Atom( const T& atom ) : element_(atom.element()), coord_orth_(atom.coord_orth()), u_aniso_orth_(atom.u_aniso_orth()) , occupancy_(atom.occupancy()), u_iso_(atom.u_iso()){}
    //! get atom element name: e.g. "C", "N", "Zn2+"
    const String& element() const { return element_; }
    //! get atom orthogonal (Angstrom) coordinate
    const Coord_orth& coord_orth() const { return coord_orth_; }
    //! get atom occupancy
    const ftype& occupancy() const { return occupancy_; }
    //! get atom orthogonal isotropic U value
    const ftype& u_iso() const { return u_iso_; }
    //! get atom orthogonal anisotropic U value
    const U_aniso_orth& u_aniso_orth() const { return u_aniso_orth_; }
    void set_element( const String& s );             //!< set element
    void set_coord_orth( const Coord_orth& s );      //!< set coord_orth
    void set_occupancy( const ftype& s );            //!< set occupancy
    void set_u_iso( const ftype& s );                //!< set u_iso
    void set_u_aniso_orth( const U_aniso_orth& s );  //!< set u_aniso
    //! apply a rotation-translation operator (RTop) to the atom
    void transform( const RTop_orth rt );
    //! test for null atom: atom is null is coord is null
    bool is_null() const { return coord_orth_.is_null(); }
    //! return null atom
    static Atom null();
  private:
    String element_;
    Coord_orth coord_orth_;
    U_aniso_orth u_aniso_orth_;
    ftype occupancy_, u_iso_;
  };


  //! Atom list class
  /*! This class defines a minimal atom list object providing only
    those properties required for an electron density calculation. It
    is a trivial derivation from std::vector<Atom>. In addition a template
    constructor allows it to be constructed from any other object with
    appropriate properties. */
  class CLIPPER_IMEX Atom_list : public std::vector<Atom>
  {
  public:
    //! null constructor
    Atom_list() {}
    //! constructor: from std::vector<Atom>
    Atom_list( const std::vector<Atom>& list ) : std::vector<Atom>( list ) {}
    //! Constructor: from vector-like list of atom-like objects
    template<class T> Atom_list( const T& list ) { for ( int i = 0; i < list.size(); i++ ) push_back( Atom( list[i] ) ); }
  };


  // some template function definitions

  /*! Construct the operator which relates one atom-list like object
    onto another. The lists must be the same size, and have the
    following properties:
    - a size() method.
    - a [int] operator, with int ranging from 0 to size()-1.
    - the object returned by the [] operator must have a coord_orth() method.
    Suitable objects include a vector of Atom, or an Atom_list. */
  template<class T> RTop_orth::RTop_orth( const T& src, const T& tgt )
  {
    std::vector<Coord_orth> vsrc( src.size() );
    std::vector<Coord_orth> vtgt( tgt.size() );
    for ( int i = 0; i < src.size(); i++ ) vsrc[i] = src[i].coord_orth();
    for ( int i = 0; i < tgt.size(); i++ ) vtgt[i] = tgt[i].coord_orth();
    (*this) = RTop_orth( vsrc, vtgt );
  }

  // some inline function definitions
  /*! Requires integer->ftype->integer transformation.
    \param op The symmetry operator
    \return The transformed coordinate */
  HKL HKL::transform( const Symop& op ) const
    { return Coord_reci_frac(*this).transform(op).hkl(); }
  /*! Optimal version.
    \param op The symmetry operator
    \return The transformed coordinate */
  HKL HKL::transform( const Isymop& op ) const
    { return op*(*this); }
  /*! Get the symmetry phase shift incurred when transforming a
    reflection by this operator.
    \param hkl The reflection HKL to transform.
    \return The phase shift. */
  ftype HKL::sym_phase_shift( const Symop& op ) const
    { return -Util::twopi()*( Coord_reci_frac(*this) * op.trn() ); }

  /*! The grid coordinate is incremented efficiently in a manner which
    is exaclty equivalent to increasing index() by 1 in a zero based grid.
    \param g The grid with which this increment is synchronised. */
  const Coord_grid& Coord_grid::next( const Grid& g )
    { w()++; if ( w() >= g.nw() ) { w() = 0; v()++; if ( v() >= g.nv() ) { v() = 0; u()++; } } return *this; }
  /*! The grid coordinate is incremented efficiently in a manner which
    is exaclty equivalent to increasing index() by 1 in a non-zero based grid.
    \param g The grid with which this increment is synchronised. */
   const Coord_grid& Coord_grid::next( const Grid_range& g )
    { w()++; if ( w() > g.max().w() ) { w() = g.min().w(); v()++; if ( v() > g.max().v() ) { v() = g.min().v(); u()++; } } return *this; }
  /*! Test whether this coordinate has been incremented using next()
    beyond the end of the specified zero based grid.
    \param g The grid concerned. */
  bool Coord_grid::last( const Grid& g ) const
    { return ( u() >= g.nu() ); }
  /*! Test whether this coordinate has been incremented using next()
    beyond the end of the specified non-zero based grid.
    \param g The grid concerned. */
  bool Coord_grid::last( const Grid_range& g ) const
    { return ( u() > g.max().u() ); }
  /*! Return the index in a 1-d array corresponding to this coordinate
    for a zero based grid.
    \param g The grid concerned.
    \return The corresponding index. */
  int Coord_grid::index( const Grid& g ) const
    { return ( u()*g.nv() + v() )*g.nw() + w(); }
  /*! Return the coordinate corresponding to a given index in a zero
    based grid.
    \param g The grid concerned.
    \return The corresponding coordinate. */
  void Coord_grid::deindex( const Grid& g, const int& index )
    { u() = index/(g.nv()*g.nw()); v() = (index/g.nw()) % g.nv(); w() = (index) % g.nw(); }

  /*! \param g The grid concerned \return The transformed coordinate. */
  Coord_grid Coord_grid::unit( const Grid_sampling& g ) const
    { return Coord_grid( Util::mod(u(), g.nu()), Util::mod(v(), g.nv()), Util::mod(w(), g.nw()) ); }
  /*! \return The non-integer coordinate. */
  Coord_map Coord_grid::coord_map() const
    { return Coord_map( *this ); }
  /*! Fractional coordinate is not normalised onto range 0..1
    \param g The grid concerned
    \return The fractional coordinate */
  Coord_frac Coord_grid::coord_frac( const Grid_sampling& g ) const
    { return Coord_frac( ftype(u())/ftype(g.nu()), ftype(v())/ftype(g.nv()), ftype(w())/ftype(g.nw()) ); }

  /*! \note Normally you would get a value through clipper::HKL_info,
    unless you specifically want a value for a different cell. */
  ftype HKL::invresolsq( const Cell& cell ) const
    { return cell.metric_reci().lengthsq( Coord_reci_frac( *this ) ); }
  /*! \return The non-integer coordinate. */
  Coord_reci_frac HKL::coord_reci_frac() const
    { return Coord_reci_frac( *this ); }
  /*! \param cell The cell concerned \return The transformed coordinate. */
  Coord_reci_orth HKL::coord_reci_orth( const Cell& cell ) const
    { return coord_reci_frac().coord_reci_orth( cell ); }
  /*! \return The inverse resolution squared. */
  ftype Coord_reci_orth::invresolsq() const
    { return xs()*xs() + ys()*ys() + zs()*zs(); }
  /*! \param cell The cell concerned \return The transformed coordinate. */
  Coord_reci_frac Coord_reci_orth::coord_reci_frac( const Cell& cell ) const
    { return Coord_reci_frac( (*this) * cell.matrix_orth() ); }
  /*! \param cell The cell concerned \return The inverse resolution squared. */
   ftype Coord_reci_frac::invresolsq( const Cell& cell ) const
    { return cell.metric_reci().lengthsq( *this ); }
  /*! \param cell The cell concerned \return The transformed coordinate. */
  Coord_reci_orth Coord_reci_frac::coord_reci_orth( const Cell& cell ) const
    { return Coord_reci_orth( (*this) * cell.matrix_frac() ); }

  /*! \return The squared length in Angstroms squared */
  ftype Coord_orth::lengthsq() const
    { return x()*x()+y()*y()+z()*z(); }
  /*! \param cell The cell concerned \return The transformed coordinate. */
  Coord_frac Coord_orth::coord_frac( const Cell& cell ) const
    { return Coord_frac( cell.matrix_frac() * (*this) ); }
  /*! \return The squared length in Angstroms squared */
  ftype Coord_frac::lengthsq( const Cell& cell ) const
    { return cell.metric_real().lengthsq( *this ); }
  /*! \param cell The cell concerned \return The transformed coordinate. */
  Coord_orth Coord_frac::coord_orth( const Cell& cell ) const
    { return Coord_orth( cell.matrix_orth() * (*this) ); }
  /*! \param g The grid concerned \return The transformed coordinate. */
  Coord_map Coord_frac::coord_map( const Grid& g ) const
    { return Coord_map( u()*ftype(g.nu()), v()*ftype(g.nv()), w()*ftype(g.nw()) ); }
  /*! \param g The grid concerned \return The transformed coordinate. */
  Coord_grid Coord_frac::coord_grid( const Grid& g ) const
    { return Coord_grid( Util::intr(u()*ftype(g.nu())), Util::intr(v()*ftype(g.nv())), Util::intr(w()*ftype(g.nw())) ); }
  /*! \param g The grid concerned \return The transformed coordinate. */
  Coord_frac Coord_map::coord_frac( const Grid& g ) const
    { return Coord_frac( u()/ftype(g.nu()), v()/ftype(g.nv()), w()/ftype(g.nw()) ); }

} // namespace clipper

#endif
