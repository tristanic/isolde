/* hkl_info.cpp: class file for hkl_info class + children */
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


#include "hkl_info.h"


namespace clipper {


Message_fatal message_recip_asu_error( "HKL_info: find_sym reciprocal space ASU error" );
Message_ctor message_ctor_hkl_info( " [HKL_info: constructed]" );


// methods for 'HKL_info'
// private methods:

/*! Update all the lookup tables to be consistent with the modified
  reflection list */
void HKL_info::update_hkl_list() {
  lookup.init( hkl );
  hkl_class_lookup.resize( num_reflections() );
  invresolsq_lookup.resize( num_reflections() );
  invresolsq_range_ = Range<ftype>();
  for ( int i = 0; i < num_reflections(); i++ ) {
    hkl_class_lookup[i] = spacegroup_.hkl_class( hkl_of( i ) );
    invresolsq_lookup[i] = hkl_of( i ).invresolsq( cell_ );
    invresolsq_range_.include( invresolsq_lookup[i] );
  }
}


// public methods:
// constructors/destructor

HKL_info::HKL_info()
{
  Message::message( message_ctor_hkl_info );
}


/*! Construct and initialise HKL_info object. This updates the
  spacegroup and cell and clears the reflection list. The resolution
  is used as a rejection criterion for reflections - no HKL will be
  stored beyond the given limit. Initially there are no reflections in
  the reflection list: see generate_hkl_list().

  If any of the parameters have null values, the existing values will
  be unchanged. The object will only be fully initialised once all
  parameters are available.
  \param spacegroup The spacegroup.
  \param cell The unit cell.
  \param resolution The resolution limit. */
HKL_info::HKL_info( const Spacegroup& spacegroup, const Cell& cell, const Resolution& resolution, const bool& generate )
{
  init( spacegroup, cell, resolution, generate );
  Message::message( message_ctor_hkl_info );
}


/*! Initialise the HKL_info object. This updates the spacegroup and
  cell and clears the reflection list. The resolution is used as a
  rejection criterion for reflections - no HKL will be stored beyond
  the given limit. Initially there are no reflections in the
  reflection list: see generate_hkl_list().

  If any of the parameters have null values, the existing values will
  be unchanged. The object will only be fully initialised once all
  parameters are available.
  \param spacegroup The spacegroup.
  \param cell The unit cell.
  \param resolution The resolution limit.
  \param generate If true, a reflection list will be generated for an ASU. */
void HKL_info::init( const Spacegroup& spacegroup, const Cell& cell, const Resolution& resolution, const bool& generate )
{
  // set spacegroup and cell
  spacegroup_ = spacegroup;
  cell_ = cell;
  resolution_ = resolution;

  // check parameters
  if ( is_null() ) return;

  // Create the intergised symops (grid irrelevent)
  Grid g( 24, 24, 24 );
  isymop.resize( spacegroup_.num_symops() );
  for ( int sym = 0; sym < spacegroup_.num_symops(); sym++ )
    isymop[sym] = Isymop( spacegroup_.symop(sym), g );

  // reflection lists
  hkl.clear();
  if ( generate ) generate_hkl_list();
  update_hkl_list();
}


/*! Initialise the HKL_info object. This updates the spacegroup and
  cell and clears the reflection list. The HKL_sampling determines
  the reflection list.

  If any of the parameters have null values, the existing values will
  be unchanged. The object will only be fully initialised once all
  parameters are available.
  \param spacegroup The spacegroup.
  \param cell The unit cell.
  \param hkl_sampling The resolution limit.
  \param generate If true, a reflection list will be generated for an ASU. */
void HKL_info::init( const Spacegroup& spacegroup, const Cell& cell, const HKL_sampling& hkl_sampling, const bool& generate )
{
  // set spacegroup and cell
  spacegroup_ = spacegroup;
  cell_ = cell;
  hkl_sampling_ = hkl_sampling;

  // check parameters
  if ( spacegroup_.is_null() || cell_.is_null() || hkl_sampling_.is_null() ) return;

  // estimate resolution
  resolution_ = hkl_sampling.resolution( cell );

  // Create the intergised symops (grid irrelevent)
  Grid g( 24, 24, 24 );
  isymop.resize( spacegroup_.num_symops() );
  for ( int sym = 0; sym < spacegroup_.num_symops(); sym++ )
    isymop[sym] = Isymop( spacegroup_.symop(sym), g );

  // reflection lists
  hkl.clear();
  if ( generate ) {
    std::vector<HKL> add;
    // make a reflection list for the given refln sampling
    HKL lim = hkl_sampling.hkl_limit();
    HKL rfl;
    // make a list of valid reflections
    for (rfl.h()=-lim.h(); rfl.h()<=lim.h(); rfl.h()++)
      for (rfl.k()=-lim.k(); rfl.k()<=lim.k(); rfl.k()++)
	for (rfl.l()=-lim.l(); rfl.l()<=lim.l(); rfl.l()++)
	  if ( spacegroup_.recip_asu(rfl) &&
	       hkl_sampling.in_resolution( rfl ) &&
	       !( spacegroup_.hkl_class(rfl).sys_abs() ) ) add.push_back(rfl);
    // update the reflection data lists
    add_hkl_list( add );
  }

  update_hkl_list();
}


/*! \return true if the object has not been initalised. */
bool HKL_info::is_null() const
{ return ( spacegroup_.is_null() || cell_.is_null() || resolution_.is_null() ); }


/*! Using current cell, spacegroup, resolution. */
void HKL_info::generate_hkl_list()
{
  std::vector<HKL> add;
  // make a reflection list for the given cell, symm and resolution
  /* TO MAKE A BOX TO HOLD A SPHERE IN REAL OR RECIPROCAL SPACE
     In reciprocal space, use the real space cell dimension / resolution
     In real space, use the recip space dimension * radius */
  HKL rfl;
  int hmax = int(cell_.descr().a()/resolution_.limit());
  int kmax = int(cell_.descr().b()/resolution_.limit());
  int lmax = int(cell_.descr().c()/resolution_.limit());
  ftype s_lim = resolution_.invresolsq_limit();
  // make a list of valid reflections
  for (rfl.h()=-hmax; rfl.h()<=hmax; rfl.h()++)
    for (rfl.k()=-kmax; rfl.k()<=kmax; rfl.k()++)
      for (rfl.l()=-lmax; rfl.l()<=lmax; rfl.l()++)
	if ( spacegroup_.recip_asu(rfl) && rfl.invresolsq(cell_) < s_lim &&
	     !( spacegroup_.hkl_class(rfl).sys_abs() ) ) add.push_back(rfl);
  // update the reflection data lists
  hkl.clear();
  add_hkl_list( add );
}


/*! The new HKLs are transformed to the default reciprocal ASU, and
  added to the reflection list. Duplicates and reflections outside the
  resoluution limit are ignored. Then the fast lookup tables for HKL,
  invresolsq, and reflection class are rebuilt.
  \param add The list of new reflections to add. */
void HKL_info::add_hkl_list( const std::vector<HKL>& add ) {
  HKL equiv; int sym; bool friedel;
  for ( int i = 0; i < add.size(); i++ ) {
    if ( add[i].invresolsq( cell_ ) <= resolution_.invresolsq_limit() ) {
      equiv = find_sym( add[i], sym, friedel );
      if ( lookup.index_of( equiv ) < 0 ) hkl.push_back( equiv );
    }
  }
  update_hkl_list();
}


/*! Returns the index of the reflection, the sym no. and Friedel flag.
 \internal */
HKL HKL_info::find_sym( const HKL& rfl, int& sym, bool& friedel ) const
{
  // find the symmetry operator mapping hkl into the reflection
  // list and determine the friedel flag
  HKL equiv;

  // now find the symop which gives a reflection from the list
  for (sym = 0; sym < spacegroup_.num_primops(); sym++) {
    equiv = rfl.transform(isymop[sym]);
    if ( spacegroup_.recip_asu(equiv) ) { friedel = false; return equiv; }
    equiv = -equiv;
    if ( spacegroup_.recip_asu(equiv) ) { friedel = true ; return equiv; }
  }
  Message::message( message_recip_asu_error );
  return equiv;
}


void HKL_info::debug() const
{
  std::cout << "Num reflns " << hkl.size() << "\n";
}


} // namespace clipper
