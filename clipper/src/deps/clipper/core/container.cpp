/* container.cpp: implementation for base class of reflection heirarchy */
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


#include "container.h"


namespace clipper {


Message_fatal message_parent_of_root( "Container: attempt to access parent of root" );
Message_fatal message_child_out_of_range( "Container: child index out of range" );
Message_fatal message_duplicate_path( "Container: attempt to create duplicate path" );


void Container::add_child( Container& c )
{
  children.push_back( &c );
}

void Container::del_child( Container& c )
{
  for ( int i = 0; i < children.size(); i++ )
    if ( children[i] == &c ) { children.erase( children.begin() + i ); break; }
}

void Container::update()
{
  // recursively update entire tree
  for ( int i = 0; i < children.size(); i++ ) child(i).update();
}

Container::Container( const String name )
{
  destroyed_with_parent = false;
  parent_ = NULL;
  if ( name != "" ) name_ = name;
  else              name_ = "unnamed";

  // log object if required
  if ( Message::message_level() <= Message_ctor::level() )
    Message::message( Message_ctor( "[Container: contructed (root)/"+name_+">" ) );
}

Container::Container( Container& parent, const String& path )
{
  destroyed_with_parent = false;
  parent_ = NULL;
  name_ = path.tail();
  // find desired parent using working root and path provided
  parent_ = parent.find_path_ptr( path.notail() );
  if ( parent_ == NULL )
    Message::message( Message_fatal( "Container: No such path- "+path ) );
  // if no name given, make one up
  if ( name_ == "" )
    for (int i = 1; i < 100; i++ ) {
      name_ = "unnamed" + String(i,2);
      if ( parent_->find_path_ptr( name_ ) == NULL ) break;
    }
  // otherwise check for dulicates
  if ( parent_->find_path_ptr( name_ ) != NULL )
    Message::message( message_duplicate_path );
  parent_->add_child(*this);

  // log object if required
  if ( Message::message_level() <= Message_ctor::level() )
    Message::message( Message_ctor( "[Container: contructed "+parent_->name()+"/"+name_+">" ) );
}

Container::~Container() {
  for ( int i = 0; i < children.size(); i++ ) {
    children[i]->parent_ = NULL;
    if ( children[i]->is_destroyed_with_parent() ) delete children[i];
  }
  if ( parent_ != NULL ) parent_->del_child(*this);

  // log object if required
  if ( Message::message_level() <= Message_dtor::level() )
    Message::message( Message_dtor( "<Container: destroyed "+name_+"]" ) );
}

String Container::path() const {
  if ( has_parent() )
    return parent().path()+"/"+name();
  else
    return "/"+name();
}

String Container::name() const
{ return name_; }

void Container::set_name( const String& name )
{ name_ = name; }

bool Container::is_destroyed_with_parent() const
{ return destroyed_with_parent; }

void Container::set_destroyed_with_parent( const bool d )
{ destroyed_with_parent=d; }

bool Container::has_parent() const
{ return ( parent_ != NULL ); }

const Container& Container::parent() const
{ 
  if ( parent_ == NULL ) Message::message( message_parent_of_root );
  return *parent_;
}

Container& Container::parent()
{ 
  if ( parent_ == NULL ) Message::message( message_parent_of_root );
  return *parent_;
}

int Container::num_children() const
{ return children.size(); }

const Container& Container::child( const int& i ) const
{
  if ( i < 0 || i >= children.size() )
    Message::message( message_child_out_of_range );
  return *children[i];
}

Container& Container::child( const int& i )
{
  if ( i < 0 || i >= children.size() )
    Message::message( message_child_out_of_range );
  return *children[i];
}

const Container& Container::ultimate_parent() const
{
  if ( has_parent() ) return parent().ultimate_parent();
  else return *this;
}

Container& Container::ultimate_parent()
{
  if ( has_parent() ) return parent().ultimate_parent();
  else return *this;
}

Container* Container::parent_ptr()
{
  return parent_;
}

Container* Container::find_path_ptr( const String& path_ )
{
  // find an object in the heirarchy using UNIX path semantics
  String path = path_;

  // check for this object
  if ( path == "" ) return this;

  // check for root-ed path
  if ( path[0] == '/' ) {
    if ( has_parent() ) {
      return parent().find_path_ptr( path );
    } else {
      path = path.nohead();
      if ( path.head() != name() ) return NULL;
      path = path.nohead();
    }
  }

  // check for a child of this object
  for ( int i = 0; i < num_children(); i++ )
    if ( path.head() == "*" || path.head() == child(i).name() )
      return child(i).find_path_ptr( path.nohead() );

  // else fail
  return NULL;
}

void Container::move( const String& path )
{
  // move this object in heirarchy to 'path'
  if ( has_parent() ) parent().del_child( *this );
  parent_ = find_path_ptr( path.notail() );
  if ( parent_ == NULL )
    Message::message( Message_fatal( "Container: No such path- "+path ) );
  parent_->add_child( *this );
  name_ = path.tail();
  update();
}

void Container::debug()
{
  std::cout << path() << "\n";
  for ( int i = 0; i < num_children(); i++ )
    child(i).debug();
}


} // namespace clipper
