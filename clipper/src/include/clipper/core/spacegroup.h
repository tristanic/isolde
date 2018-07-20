/*! \file lib/spacegroup.h
    Header file for spacegroup symmetry class
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
/* some code & data converted from CCP4 library symlib.f (various authors) */


#ifndef CLIPPER_SPACEGROUP
#define CLIPPER_SPACEGROUP


#include "symop.h"
#include "spacegroup_data.h"
#include "../imex.h"

namespace clipper {

  // forward definitions
  class HKL;
  class HKL_class;
  class Coord_frac;


  //! spacegroup description
  /*! The spacegroup description is a compact description of a
    spacegroup. It may be initialised from Hall or H-M symbols, a
    string of symops or a number. Internally a hash code is used to
    refer to the spacegroup, so this object is only 32 bits in
    size.

    For more details of spacegroup symbols, see Sydney R. Hall & Ralf
    W. Grosse-Kunstleve 'Concise Space-Group Symbols',
    http://www.kristall.ethz.ch/LFK/software/sginfo/hall_symbols.html
  */
  class CLIPPER_IMEX Spgr_descr
  {
  public:
    enum TYPE { Hall, HM, XHM, Symops, Number, Unknown };
    //! null constructor
    Spgr_descr();
    //! constructor: from symbol or operators.
    explicit Spgr_descr( const String& symb, TYPE type = Unknown );
    //! constructor: from number.
    explicit Spgr_descr( const int& num );
    //! return the spacegroup number
    int spacegroup_number() const;
    //! return the Hall symbol
    String symbol_hall() const;
    //! return the H-M symbol
    String symbol_hm() const;
    //! return the extended H-M symbol
    String symbol_xhm() const;
    //! return the extension H-M symbol
    String symbol_hm_ext() const;
    //! set preferred default spacegroup choice
    static void set_preferred( const char& c );

    //! Vector of symop codes and associated methods
    class CLIPPER_IMEX Symop_codes : public std::vector<Symop_code>
    {
    public:
      //! initialise from Hall symbol
      void init_hall( const String& symb );
      //! initialise from symops
      void init_symops( const String& symb );
      //! expand (incomplete) list of symops
      Symop_codes expand() const;
      //! return primitive non-inversion ops (by computation)
      Symop_codes primitive_noninversion_ops() const;
      //! return inversion ops (by computation)
      Symop_codes inversion_ops() const;
      //! return primitive incl inversion ops (by computation)
      Symop_codes primitive_ops() const;
      //! return lattice centering ops (by computation)
      Symop_codes centering_ops() const;
      //! return Laue ops
      Symop_codes laue_ops() const;
      //! return point group ops
      Symop_codes pgrp_ops() const;
      //! return Patterson ops
      Symop_codes patterson_ops() const;
      //! return minimal list of generator ops
      Symop_codes generator_ops() const;
      //! return product of this (expanded) list by another (expanded) list
      Symop_codes product( const Symop_codes& ops2 ) const;
      //! return hash code of symop list
      unsigned int hash() const;
    };

    //! constructor: from symop list.
    explicit Spgr_descr( const Symop_codes& ops );
    //! return the generators for the spacegroup
    const Symop_codes& generator_ops() const { return generators_; }
    //! return the hash code for the spacegroup \internal
    const unsigned int& hash() const { return hash_; }

  protected:
    unsigned int hash_;       //!< hash code of spacegroup
    Symop_codes generators_;  //!< codes for symop generators

	static char pref_12;
	static char pref_hr;  //!< preferred origin and hex/romb symbols
  };


  // ObjectCache data type
  class CLIPPER_IMEX Spgr_cacheobj
  {
  public:
    typedef Spgr_descr Key;
    Spgr_cacheobj( const Key& spgr_cachekey );  //!< construct entry
    bool matches( const Key& spgr_cachekey ) const; //!< compare entry
    String format() const;  //!< string description
    // data
    Key spgr_cachekey_;                 //!< spacegroup cachekey
    int nsym, nsymn, nsymi, nsymc, nsymp;  //!< number of syms: total, primitive
    int lgrp;                           //!< Laue group number
    std::vector<Symop>  symops;         //!< symmetry operators
    std::vector<Isymop> isymops;        //!< symmetry operators
    Vec3<> asu_min_, asu_max_;          //!< real space ASU
    static Mutex mutex;                 //!< thread safety
  };


  //! Spacegroup object
  /*! The spacegroup object is a full description of a spacegroup,
    including all the most regularly used information in an efficient
    form. It may be initialised from a clipper::Spgr_descr. This
    object.

    For more details of spacegroup symbols, see Sydney R. Hall & Ralf
    W. Grosse-Kunstleve 'Concise Space-Group Symbols',
    http://www.kristall.ethz.ch/LFK/software/sginfo/hall_symbols.html
  */
  class CLIPPER_IMEX Spacegroup : public Spgr_descr
  {
   public:
    //! enumeration for fast construction of Null or P1 spacegroup
    enum TYPE { Null, P1 };
    //! enumeration for cell axes
    enum AXIS { A=0, B=1, C=2 };
    //! null constructor
    Spacegroup() {};
    //! constructor: fast constructor for Null or P1 spacegroup
    explicit Spacegroup( TYPE type );
    //! constructor: from spacegroup description
    explicit Spacegroup( const Spgr_descr& spgr_descr );
    //! initialiser: from spacegroup description
    void init( const Spgr_descr& spgr_descr );

    //! test if object has been initialised
    bool is_null() const;

    // methods
    //! get spacegroup description
    inline const Spgr_descr& descr() const { return (*this); }
    //! get number of symops
    inline const int& num_symops() const { return nsym; }
    //! get number of primitive symops (identical to num_primitive_symops())
    inline const int& num_primops() const { return num_primitive_symops(); }
    //! get number of primitive symops (inc identity and inversion)
    inline const int& num_primitive_symops() const { return nsymp; }
    //! get number of centering symops (inc identity)
    inline const int& num_centering_symops() const { return nsymc; }
    //! get number of inversion symops (inc identity)
    inline const int& num_inversion_symops() const { return nsymi; }
    //! get number of primitive non-inversion symops (inc identity)
    inline const int& num_primitive_noninversion_symops() const { return nsymn;}
    //! get n'th symop
    inline const Symop& symop( const int& sym_no ) const
      { return symops[sym_no]; }
    //! get n'th primitive symop (identical to symop(sym_no))
    inline const Symop& primitive_symop( const int& sym_no ) const
      { return symops[sym_no]; }
    //! get n'th inversion symop (0...1 max)
    inline const Symop& inversion_symop( const int& sym_no ) const
      { return symops[nsymn*sym_no]; }
    //! get n'th centering symop (0...3 max)
    inline const Symop& centering_symop( const int& sym_no ) const
      { return symops[nsymp*sym_no]; }
    //! get the order of rotational symmetry about a given axis
    int order_of_symmetry_about_axis( const AXIS axis ) const;

    //! get 'class' of reflection: multiplicity, allowed phase, absence
    HKL_class hkl_class( const HKL& hkl ) const;
    //! test if hkl is in default reciprocal ASU
    bool recip_asu( const HKL& hkl ) const;

    //! get symop number corresponding to the product of two symops
    int product_op( const int& s1, int& s2 ) const;
    //! get symop number corresponding to the inverse of a symop
    int inverse_op( const int& s ) const;

    //! get map ASU, upper bound
    Coord_frac asu_max() const;
    //! get map ASU, lower bound
    Coord_frac asu_min() const;

    //! test if change of hand preserves spacegroup
    bool invariant_under_change_of_hand() const;

    // inherited functions listed for documentation purposes
    //-- int spacegroup_number() const;
    //-- String symbol_hall() const;
    //-- String symbol_hm() const;
    //! return the Laue group symbol
    String symbol_laue() const;

    //! Return P1 spacegroup
    static Spacegroup p1() { return Spacegroup( P1 ); }
    //! Return null spacegroup
    static Spacegroup null() { return Spacegroup( Null ); }

    void debug() const;

  private:
    ObjectCache<Spgr_cacheobj>::Reference cacheref;  //!< object cache reference
    const Symop* symops;    //!< fast access ptr
    const Isymop* isymops;  //!< fast access ptr
    data::ASUfn asufn;      //!< fast access ptr
    int nsym, nsymn, nsymi, nsymc, nsymp; //!< fast access copies
  };


} // namespace clipper

#endif
