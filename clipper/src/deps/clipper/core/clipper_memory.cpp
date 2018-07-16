/* clipper_memory.cpp: implementation file for clipper helper functions */
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


#include "clipper_memory.h"


namespace clipper {


//Mutex PropertyManager::mutex = Mutex();

/*! Makes copies of all property objects. */
PropertyManager::PropertyManager( const PropertyManager& mgr )
{ copy( mgr ); }

/*! Clears manager then makes copies of all property objects. */
PropertyManager& PropertyManager::operator =( const PropertyManager& mgr )
{ return copy( mgr ); }

/*! Deletes all stored properties. */
PropertyManager::~PropertyManager()
{
  //mutex.lock();
  for ( int i = 0; i < property_.size(); i++ ) delete property_[i].second;
  property_.clear();
  //mutex.unlock();
}

/*! This function is used by the copy constructor and assignement
  operator and is also useful for derived classes. */
PropertyManager& PropertyManager::copy( const PropertyManager& mgr )
{
  //mutex.lock();
  for ( int i = 0; i < property_.size(); i++ ) delete property_[i].second;
  property_.clear();
  for ( int i = 0; i < mgr.property_.size(); i++ )
    property_.push_back( std::pair<std::string,Property_base*>(
      mgr.property_[i].first, mgr.property_[i].second->clone() ) );
  //mutex.unlock();
  return *this;
}

/* \param label The label of the property to be returned.
   \param property The property object.
   \return true on success. */
bool PropertyManager::set_property( const std::string& label, const Property_base& property )
{
  //mutex.lock();
  property_.push_back( std::pair<std::string,Property_base*>( label, property.clone() ) );
  //mutex.unlock();
  return true;
}

/*! \param label The label of the property to be returned.
    \return the property object. */
const Property_base& PropertyManager::get_property( const std::string& label ) const
{
  //mutex.lock();
  const Property_base* result = NULL;
  for ( int i = 0; i < property_.size(); i++ )
    if ( label == property_[i].first ) {
      result = property_[i].second;
      break;
    }
  //mutex.unlock();
  if ( result == NULL )
    Message::message( Message_fatal( "PropertyManager: label not found.\n" ) );
  return *result;
}

/*! \param label The label of the property to be tested.
    \return true on success. */
bool PropertyManager::exists_property( const std::string& label ) const
{
  //mutex.lock();
  bool result = false;
  for ( int i = 0; i < property_.size(); i++ )
    if ( label == property_[i].first ) {
      result = true;
      break;
    }
  //mutex.unlock();
  return result;
}

/* \param label The label of the property to be deleted.
   \return true on success. */
bool PropertyManager::delete_property( const std::string& label )
{
  //mutex.lock();
  bool result = false;
  for ( int i = 0; i < property_.size(); i++ )
    if ( label == property_[i].first ) {
      delete property_[i].second;
      property_.erase( property_.begin() + i );
      result = true;
      break;
    }
  //mutex.unlock();
  return result;
}


// template compilations

template class Property<int>;
template class Property<float>;
template class Property<double>;
template class Property<std::string>;


} // namespace clipper
