/* clipper_instance.cpp: implementation file for clipper helper functions */
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


#include "clipper_instance.h"


namespace clipper {

  // template function definitions for object cache

  template<class T> ObjectCache<T>::Reference::Reference( const Reference& other )
  {
    T::mutex.lock();
    obj_ = other.obj_;
    if ( !is_null() ) obj_->first++;
    T::mutex.unlock();
  }

  template<class T> ObjectCache<T>::Reference::~Reference()
  {
    T::mutex.lock();
    if ( !is_null() ) obj_->first--;
    T::mutex.unlock();
  }

  template<class T> void ObjectCache<T>::Reference::operator =( const Reference& other )
  {
    T::mutex.lock();
    if ( !is_null() ) obj_->first--;
    obj_ = other.obj_;
    if ( !is_null() ) obj_->first++;
    T::mutex.unlock();
  }

  template<class T> ObjectCache<T>::ObjectCache()
    { mode_ = NORMAL; }

  template<class T> ObjectCache<T>::~ObjectCache()
  {
    for ( int i = 0; i < cache_.size(); i++ ) {
      if ( cache_[i]->first != 0 ) {
	std::string num( "0000" );
	num[3] += cache_[i]->first % 10;
	num[2] += cache_[i]->first / 10 % 10;
	num[1] += cache_[i]->first / 100 % 10;
	num[0] += cache_[i]->first / 1000;
	Message::message( Message_warn( "ObjectCache: Leaked "+num+" refs to <"+cache_[i]->second.format()+">" ) );
      }
    }
  }

  template<class T> typename ObjectCache<T>::Reference ObjectCache<T>::cache( const typename T::Key& key )
  {
    T::mutex.lock();
    std::pair<int,T>* ptr = NULL;
    // find existing data
    for ( int i = 0; i < cache_.size(); i++ )
      if ( cache_[i]->second.matches(key) )
        ptr = cache_[i];
    // none found, make new data
    if ( ptr == NULL ) {
      // pick an appropriate place
      if ( mode_ == MINMEM )  // MINMEM: remove unreferenced
        purge();
      if ( mode_ == NORMAL )  // NORMAL: replace unreferenced
        for ( int i = 0; i < cache_.size(); i++ )
          if ( cache_[i]->first == 0 ) {
    	ptr = cache_[i];
    	ptr->second = T(key);
    	break;
          }
      // otherwise add new
      if ( ptr == NULL ) {
        ptr = new std::pair<int,T>( 0, T(key) );
        cache_.push_back( ptr );
      }
    }
    Reference result( ptr );
    T::mutex.unlock();
    return result;  // we have a ref to the new obj, so it is thread safe
  }

  template<class T> void ObjectCache<T>::set_mode( const MODE& mode )
  { mode_ = mode; }

  template<class T> void ObjectCache<T>::purge()
  {
    for ( int i = cache_.size() - 1; i >= 0; i-- )
      if ( cache_[i]->first == 0 ) {
	delete cache_[i];
	cache_.erase( cache_.begin() + i );
      }
  }

  template<class T> void ObjectCache<T>::destroy()
  { cache_.clear(); }

  template<class T> void ObjectCache<T>::debug() const
  {
    for ( int i = 0; i < cache_.size(); i++ )
      std::cout << "Cache pos: " << i << "\t   Refs: " << cache_[i]->first
		<< "\t" << cache_[i]->second.format() << "\n";
  }


  template class CLIPPER_IMEX ObjectCache<Spgr_cacheobj>;
  template class CLIPPER_IMEX ObjectCache<HKL_data_cacheobj>;
  template class CLIPPER_IMEX ObjectCache<Xmap_cacheobj>;


  // clipper instance

  ClipperInstance ClipperInstantiator::inst;

  ClipperInstance::ClipperInstance()
  {
  }

  ClipperInstance::~ClipperInstance()
  {
    xmcache_.purge();
    hdcache_.purge();
    sgcache_.purge();
  }

  void ClipperInstance::destroy()
  {
    xmcache_.destroy();
    hdcache_.destroy();
    sgcache_.destroy();
  }

} // namespace clipper
