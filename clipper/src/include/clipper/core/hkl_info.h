/*! \file lib/hkl_info.h
    Header file for hkl list class
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


#ifndef CLIPPER_HKL_INFO
#define CLIPPER_HKL_INFO


#include "hkl_lookup.h"
#include "../imex.h"

namespace clipper
{
  class HKL_data_base;


  //! HKL list container and tree root
  /*! This object contains contains a reflection list, and all the
    properties on which such a list depends, i.e. spacegroup, cell,
    resolution. It also keeps a fast reflection lookup list and lookup
    lists for resolutions and reflection classes. */
  class CLIPPER_IMEX HKL_info
  {
   public:
    //! null constructor
    HKL_info();
    //! constructor: Takes spacegroup, cell, and resolution
    HKL_info( const Spacegroup& spacegroup, const Cell& cell, const Resolution& resolution, const bool& generate = false );
    //! initialiser: Takes spacegroup, cell, and resolution
    void init( const Spacegroup& spacegroup, const Cell& cell, const Resolution& resolution, const bool& generate = false );
    //! initialiser: Takes spacegroup, cell, and HKL_sampling
    void init( const Spacegroup& spacegroup, const Cell& cell, const HKL_sampling& hkl_sampling, const bool& generate = true );

    //! test if object has been initialised
    bool is_null() const;

    //! get the cell
    const Cell& cell() const { return cell_; }
    //! get the spacegroup
    const Spacegroup& spacegroup() const { return spacegroup_; }
    //! [CLIPPER2] get HKL_sampling
    const HKL_sampling& hkl_sampling() const { return hkl_sampling_; }
    //! get the resolution
    const Resolution& resolution() const { return resolution_; }

    //! synthesize hkl list
    void generate_hkl_list();
    //! add new reflections to the list
    void add_hkl_list( const std::vector<HKL>& add );

    //! get number of reflections in the object
    inline int num_reflections() const { return hkl.size(); }

    //! reflection hkl from index
    /*! \param index The index. \return The corresponding HKL. */
    inline const HKL& hkl_of( const int& index ) const { return hkl[index]; }
    //! reflection index from hkl
    /*! This does not check symmetry equivalents (see find_sym).
     \param rfl The HKL. \return The index, or -1 if it does not exist. */
    inline int index_of( const HKL& rfl ) const
      { return lookup.index_of( rfl ); }

    //! get reflection resolution using lookup
    inline const ftype32& invresolsq( const int& index ) const
      { return invresolsq_lookup[index]; }
    //! get resolution limits of the list
    inline const Range<ftype>& invresolsq_range() const
      { return invresolsq_range_; }
    //! get reflection class using lookup
    const HKL_class& hkl_class( const int& index ) const
      { return hkl_class_lookup[index]; }
    //! find symop no and friedel to bring an HKL into ASU
    HKL find_sym( const HKL& rfl, int& sym, bool& friedel ) const;


    //! HKL reference base class
    /*! This is a reference to an HKL. It forms a base class for
      index-like and coordinate-like HKL references. If you write a
      method which will work with either, then specify this instead of
      either of the derived classed. \internal */
    class HKL_reference_base
    {
    public:
      //! return the parent HKL_info
      inline const HKL_info& base_hkl_info() const { return *hklinfo; }
      //! return the current index (-1 if invalid)
      inline const int& index() const { return index_; }
      //! return the inv resol sq for the reflection (assumes index valid)
      inline ftype invresolsq( const HKL_data_base& hkldata ) const;
      //! return the inv resol sq for the reflection (assumes index valid)
      inline ftype invresolsq() const
	{ return hklinfo->invresolsq( index_ ); }
      //! test if index has gone past last reflection
      inline bool last() const
	{ return ( index_ >= hklinfo->num_reflections() ); }
    protected:
      const HKL_info* hklinfo;
      int index_;
    };

    //! HKL reference with index-like behaviour
    /*! This is a reference to an HKL. It behaves like a simple index
      into the reflection list, but can be easily converted into an
      HKL as and when required. It also implements methods for
      iterating through a reflection list.

      \note The following methods are inherited from
      HKL_reference_base but are documented here for convenience:
      base_hkl_info(), index(), invresolsq(), last().
    */
    class HKL_reference_index : public HKL_reference_base
    {
    public:
      //! Null constructor
      HKL_reference_index() {}
      //! Constructor: takes parent HKL_info and initial index
      HKL_reference_index( const HKL_info& hklinfo_, const int& index )
	{ hklinfo = &hklinfo_; index_ = index; }
      //! return the current HKL
      inline const HKL& hkl() const { return hklinfo->hkl_of( index_ ); }
      //! return the reflection class for the reflection
      inline const HKL_class& hkl_class() const
	{ return hklinfo->hkl_class( index_ ); }
      //! increment to next reflection
      inline HKL_reference_index& next() { index_++; return *this; }
      // inherited functions listed for documentation purposes
      //-- const HKL_info& base_hkl_info() const;
      //-- const int& index() const;
      //-- const ftype invresolsq() const;
      //-- bool last() const;
    };

    //! HKL reference with coord-like behaviour
    /*! This is a reference to an HKL. It behaves like an HKL, but
      also stores the index of the corresponding reflection in the
      reflection list, if it exists, and the symmetry and friedel
      operators required to get there.

      \note The following methods are inherited from
      HKL_reference_base but are documented here for convenience:
      base_hkl_info(), index(), invresolsq(), last().
    */
    class HKL_reference_coord : public HKL_reference_base
    {
    public:
      //! Null constructor
      HKL_reference_coord() {}
      //! Constructor: takes parent HKL_info and initial HKL
      HKL_reference_coord( const HKL_info& hklinfo_, const HKL& hkl ) {
	hklinfo = &hklinfo_;
	hkl_ = hkl;
	index_ = hklinfo->index_of( hklinfo->find_sym( hkl_, sym_, friedel_ ) );
	if ( index_ < 0 ) Message::message( Message_fatal( "HKL_reference_coord: hkl not found" ) );
      }
      //! return the current HKL
      inline const HKL& hkl() const { return hkl_; }
      //! get current symop number
      inline const int& sym() const { return sym_; }
      //! get current friedel flag
      inline const bool& friedel() const { return friedel_; }
      //! assign from HKL
      /*! The new HKL must exist in the reflection list. The
	calculation is optimised for the case when the new HKL is near
	the old one. */
      inline HKL_reference_coord& set_hkl( const HKL& hkl__ ) {
	hkl_ = hkl__;
	HKL equiv = hkl__.transform(hklinfo->isymop[sym_]);
	if ( friedel_ ) equiv = -equiv;
	index_ = hklinfo->index_of( equiv );
	if ( index_ < 0 ) index_ =
	   hklinfo->index_of( hklinfo->find_sym( hkl_, sym_, friedel_ ) );
	return *this;
      }
      //! increment to next reflection
      inline HKL_reference_coord& next() {
	sym_ = 0; friedel_ = false;
	index_++;
	if ( !last() ) hkl_ = hklinfo->hkl_of( index_ );
	return *this;
      }
      // increment h,k,l
      inline HKL_reference_coord& next_h() { hkl_.h()++; set_hkl( hkl_ ); return *this; }  //!< increment h
      inline HKL_reference_coord& next_k() { hkl_.k()++; set_hkl( hkl_ ); return *this; }  //!< increment k
      inline HKL_reference_coord& next_l() { hkl_.l()++; set_hkl( hkl_ ); return *this; }  //!< increment l
      inline HKL_reference_coord& prev_h() { hkl_.h()--; set_hkl( hkl_ ); return *this; }  //!< decrement h
      inline HKL_reference_coord& prev_k() { hkl_.k()--; set_hkl( hkl_ ); return *this; }  //!< decrement k
      inline HKL_reference_coord& prev_l() { hkl_.l()--; set_hkl( hkl_ ); return *this; }  //!< decrement l
      //! operator assign from HKL
      inline HKL_reference_coord& operator =( const HKL& hkl__ )
	{ return set_hkl( hkl__ ); }
      // inherited functions listed for documentation purposes
      //-- const HKL_info& base_hkl_info() const;
      //-- const int& index() const;
      //-- const ftype invresolsq() const;
      //-- bool last() const;
    protected:
      HKL hkl_;
      int sym_;
      bool friedel_;
    };

    //! return HKL_reference_index pointing to first reflection
    HKL_reference_index first() const { return HKL_reference_index( *this, 0 ); }

    void debug() const;

   protected:
    Spacegroup spacegroup_;             //!< spacegroup
    Cell cell_;                         //!< unit cell
    HKL_sampling hkl_sampling_;         //!< hkl sampling
    Resolution resolution_;             //!< resolution limit
    std::vector<Isymop> isymop;         //!< integer symops

    //! the reflection list
    std::vector<HKL> hkl;
    //! fast epsilon/centricity lookup table
    std::vector<HKL_class> hkl_class_lookup;
    //! fast resolution lookup table
    std::vector<ftype32> invresolsq_lookup;

    //! fast reflection lookup table
    HKL_lookup lookup;
    //! resolution limit of the current reflection list
    Range<ftype> invresolsq_range_;

    // internal methods:
    //! finalise reflection list \internal
    void update_hkl_list();

    friend class HKL_info::HKL_reference_base;
    friend class HKL_info::HKL_reference_index;
    friend class HKL_info::HKL_reference_coord;
  };


} // namespace clipper

#endif
