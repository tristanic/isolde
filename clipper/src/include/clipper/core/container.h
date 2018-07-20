/*! \file lib/container.h
    Header for generic container object
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


#ifndef CLIPPER_CONTAINER
#define CLIPPER_CONTAINER


#include "coords.h"
#include "../imex.h"

namespace clipper
{


  //! Definition for a generic container Object
  /*! Container is a definition for a generic container object with a
    name, parents, and children. Any object that wants to be part of the
    tree simply subclasses this class. The class also implements search
    and move objects. The tree is navigate using unix-like pathnames. A
    recursive update method can be overridden to update content after
    altering the hierarchy.

    The top container in a tree is created by passing
    Container() as its parent.
  */

  class CLIPPER_IMEX Container
  {
  public:
    //! constructor: make null object or top object in a tree
    Container( const String name = "" );
    //! constructor: from any other member and a relative path
    Container( Container& parent, const String& path );

    //! update: hierarchical content update function
    virtual void update();

    //! get the path of this tree object
    String path() const;
    //! get the name of this tree object
    String name() const;
    //! set the name of this tree object
    void set_name( const String& name );
    //! is this object to be destroyed when parent is destroyed?
    bool is_destroyed_with_parent() const;
    //! set this object to be destroyed when parent is destroyed
    void set_destroyed_with_parent( const bool d=true );
    //! 'move' method moves this object to somewhere else in the hierarchy
    void move( const String& path );

    //! test if this object has a parent
    bool has_parent() const;
    //! get the parent of this object
    const Container& parent() const;
    //! get the parent of this object
    Container& parent();

    //! return number of children
    int num_children() const;
    //! get the i'th child of this object
    const Container& child( const int& i ) const;
    //! get the i'th child of this object
    Container& child( const int& i );

    //! get the ultimate parent of this object - the top of the tree
    const Container& ultimate_parent() const;
    //! get the ultimate parent of this object - the top of the tree
    Container& ultimate_parent();

    //! get the parent of this object (NULL on fail)
    Container* parent_ptr();
    //! search up the tree for a parent of the specified type (NULL on fail)
    template<class T> T* parent_of_type_ptr();
    //! find an object using a directory-like path (NULL on fail)
    Container* find_path_ptr( const String& path );

    //! destructor: virtual
    virtual ~Container();

    void debug();

  private:
    // members
    String name_;
    Container* parent_;
    std::vector<Container*> children;
    bool destroyed_with_parent;

    //! add a child: called by a child on construction
    void add_child( Container& c );
    //! delete a child: called by a child on destruction
    void del_child( Container& c );

    //! override copy constructor to prevent its use
    Container( const Container& ) {}
  };



  // template function definitions

  template<class T> T* Container::parent_of_type_ptr()
  {
    // this can be done recursively with references, but gcc won't compile it.
    Container* p = this;
    while ( p->has_parent() ) {
      p = &(p->parent());
      T* pt = dynamic_cast<T*>(p);
      if ( pt != NULL ) return pt;
    }
    return NULL;
  }


} // namespace clipper

#endif
