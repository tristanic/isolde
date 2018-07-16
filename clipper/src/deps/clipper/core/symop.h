/*! \file lib/symop.h
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


#ifndef CLIPPER_SYMOP
#define CLIPPER_SYMOP


#include "clipper_types.h"
#include "../imex.h"

namespace clipper
{
  // forward definitions
  class Cell;
  class Grid;
  class RTop_orth;


  //! Fractional operator class.
  /*! This class is used for any RT-operator which operates on
    fractional coordinates. For a full list of methods, see
    clipper::RTop */
  class CLIPPER_IMEX RTop_frac : public RTop<>
  {
  public:
    //! null constructor
    inline RTop_frac() {}
    //! constructor: copy/convert
    inline explicit RTop_frac( const RTop<>& o ) : RTop<>( o ) {}
    //! constructor: from rotation
    inline explicit RTop_frac( const Mat33<>& r ) : RTop<>( r ) {}
    //! constructor: from string description
    explicit RTop_frac( const String& strop );
    //! constructor: from rotation and translation
    inline RTop_frac( const Mat33<>& r, const Vec3<>& t ) : RTop<>( r, t ) {}
    //! fractional-orthogonal conversion
    RTop_orth rtop_orth( const Cell& cell ) const;
    //! inverse operator
    RTop_frac inverse() const;
    //! return identity operator
    static RTop_frac identity();
    //! return null (uninitialised) operator
    static RTop_frac null();
  };


  //! Crystallographic symmetry operator
  /*! This is identical to a fractional RTop, but has its own class
    since not all fractional RTops are symops. For a full list of
    methods, see clipper::RTop and clipper::RTop_frac */
  class CLIPPER_IMEX Symop : public RTop_frac
  {
  public:
    //! null constructor
    inline Symop() {}
    //! constructor: RTop
    explicit Symop( const RTop<>& rt );
    //! constructor: from 4x4 matrix
    explicit Symop( const ftype mat[4][4] );
    //! return formatted String representation
    String format() const;
  };


  //! Integerised symmetry matrix
  /*! This is used for optimised calculations in real and reciprocal space */
  class CLIPPER_IMEX Isymop : public RTop<int>
  {
  public:
    //! null constructor
    inline Isymop() {}
    //! constructor: RTop
    inline explicit Isymop( const RTop<int>& rt ) : RTop<int>(rt) {}
    //! constructor
    Isymop( const Symop& symop, const Grid& grid );
  };


  //! Compressed encoded symmetry operator
  /*! This is a compresses representation of a crystallographic
    symmetry operator, stored as a single 32-bit integer. It may be
    converted to or from a symop or an int and compared, sorted, etc.
    The following guarantees are made concerning the code:

    - The identity operator has a code of zero.
    - Operators with non-identity rotations will have higher codes
    than operators with identity rotations, for the same translation.
    - Operators with non-zero translations will have higher codes than
    operators with zero translations, for the same rotation. */
  class CLIPPER_IMEX Symop_code
  {
  public:
    //! null constructor
    inline Symop_code() {}
    //! constructor: from int
    inline explicit Symop_code( const int& code ) : code_(code) {}
    //! constructor: from Symop
    explicit Symop_code( const Symop& op );
    //! constructor: from Isymop
    explicit Symop_code( const Isymop& op );
    //! initialiser: from Isymop
    void init( const Isymop& op );
    Symop_code code_rot() const;  //!< return code for rotation part
    Symop_code code_trn() const;  //!< return code for translation part
    Symop symop() const;          //!< convert to symop
    Isymop isymop() const;        //!< convert to integerised symop
    static Symop_code identity() { return Symop_code(0); } //!< identity code 
    //! convert to integer
    inline operator int() const { return code_; }
  private:
    int code_;
  };


} // namespace clipper

#endif
