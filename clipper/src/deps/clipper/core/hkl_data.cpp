/* hkl_data.cpp: class file for reflection data class + children */
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


#include "hkl_data.h"
#include "clipper_instance.h"


namespace clipper {


Message_ctor message_ctor_hkl_data( " [HKL_data: constructed]" );


Mutex HKL_data_cacheobj::mutex = Mutex();

HKL_data_cacheobj::HKL_data_cacheobj( const Key& hkl_data_cachekey ) : key( hkl_data_cachekey )
{
  init( Spacegroup( hkl_data_cachekey.spgr_descr() ),
	Cell( hkl_data_cachekey.cell_descr() ),
	hkl_data_cachekey.hkl_sampling(), true );
}

bool HKL_data_cacheobj::matches( const Key& hkl_data_cachekey ) const
{
  return ( key.spgr_descr().hash() == hkl_data_cachekey.spgr_descr().hash() &&
	   key.hkl_sampling() == hkl_data_cachekey.hkl_sampling() );
}

String HKL_data_cacheobj::format() const
{
  return key.spgr_descr().symbol_hall() + " " + key.hkl_sampling().format();
}


/*! For later initialisation: see init() */
HKL_data_base::HKL_data_base()
{
  parent_hkl_info = NULL;
  parent_cell = NULL;
  Message::message( message_ctor_hkl_data );
}


/*! Initialise the object using a given reflection list and cell.
  \param hkl_info The reflection list object.
  \param cell The unit cell for this datalist. */
void HKL_data_base::init( const HKL_info& hkl_info, const Cell& cell )
{
  parent_hkl_info = &hkl_info;
  parent_cell = &cell;
  cell_matches_parent = cell.equals( hkl_info.cell(), 0.5 );
}

/*! Initialise the object using a given reflection list and cell.
  \param hkl_data Object from which to inherit spacegrpoup, cell, sampling. */
void HKL_data_base::init( const HKL_data_base& hkl_data )
{
  (*this) = hkl_data;
}

/*! Initialise the object using a given spacegroup, cell, and sampling.
  \param spacegroup The spacegroup for this datalist.
  \param cell The unit cell for this datalist.
  \param hkl_sampling The reflection list description. */
void HKL_data_base::init( const Spacegroup& spacegroup, const Cell& cell, const HKL_sampling& hkl_sampling )
{
  // set parameters
  spacegroup_ = spacegroup;
  cell_ = cell;
  hkl_sampling_ = hkl_sampling;

  // check parameters
  if ( spacegroup_.is_null() || cell_.is_null() || hkl_sampling_.is_null() ) return;

  // estimate resolution
  resolution_ = hkl_sampling.resolution( cell );

  // get cache ref
  HKL_data_cacheobj::Key key( spacegroup, cell, hkl_sampling );
  cacheref = ClipperInstantiator::instance().hkl_data_cache().cache( key );

  // store legacy and fast access pointers
  init( cacheref.data(), cell_ );
}

/*! \return true if the object has not been initalised. */
bool HKL_data_base::is_null() const
{
  if ( parent_hkl_info != NULL && parent_cell != NULL )
    return ( parent_hkl_info->is_null() || parent_cell->is_null() );
  else
    return true;
}


/*! Return the resolution of a particular reflection. If the cell of
  this list closely matches (to within 0.5A) the cell of the parent
  list, this is a simple lookup, otherwise a metric calculation is
  required. */
ftype HKL_data_base::invresolsq( const int& index ) const
{
  if ( cell_matches_parent )
    return base_hkl_info().invresolsq(index);
  else
    return base_hkl_info().hkl_of(index).invresolsq(base_cell());
}


/*! \return The high and low resolution limits of the non-missing data. */
Range<ftype> HKL_data_base::invresolsq_range() const
{
  Range<ftype> slim;
  HKL_info::HKL_reference_index ih;
  for ( ih = first_data(); !ih.last(); next_data(ih) )
    slim.include( invresolsq( ih.index() ) );
  return slim;
}


/*! \return The number of non-missing data in the object. */
int HKL_data_base::num_obs() const
{
  int num = 0;
  HKL_info::HKL_reference_index ih;
  for ( ih = first_data(); !ih.last(); next_data(ih) ) num++;
  return num;
}


/*! \return HKL reference to the first data in this object. */
HKL_info::HKL_reference_index HKL_data_base::first() const
{ return HKL_reference_index( *parent_hkl_info, 0 ); } 


/*! \return HKL reference to the first non-missing data in this object. */
HKL_info::HKL_reference_index HKL_data_base::first_data() const
{ HKL_reference_index it( *parent_hkl_info, -1 ); return next_data(it); }


/*! \param ih The HKL reference to increment.
  \return HKL reference to the next non-missing data in this object. */
HKL_info::HKL_reference_index& HKL_data_base::next_data( HKL_info::HKL_reference_index& ih ) const
{
  do {
    ih.next(); if ( ih.last() ) break;
  } while ( missing( ih.index() ) );
  return ih;
}


void HKL_data_base::debug() const
{
  base_hkl_info().debug();
  base_hkl_info().spacegroup().debug();
  base_hkl_info().cell().debug();
  base_cell().debug();
}


} // namespace clipper
