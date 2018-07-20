/*! \file lib/clipper_thread.h
    Header file for clipper threading handlers
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


#ifndef CLIPPER_THREAD
#define CLIPPER_THREAD


#include "clipper_message.h"
#include "clipper_sysdep.h"
#include "../imex.h"

namespace clipper
{


#ifndef CLIPPER_DISABLE_THREADS


  //! Mutex class: used for locking and unlocking shared resources.
  /*! Create a mutex for any sharted resource, i.e. non-stack object
    used by a multi-threaded program. The lock and unlock methods lock
    that resource. Recursive locks are not allowed. */
  class CLIPPER_IMEX Mutex {
  public:
    //! constructor: create the mutex
    Mutex()              { CLIPPER_MUTEX_INIT( &mutex ); }
    //! destructor: destroy the mutex
    ~Mutex()             { CLIPPER_MUTEX_FREE( &mutex ); }
    //! lock the mutex
    inline void lock()   { CLIPPER_MUTEX_LOCK( &mutex ); }
    //! unlock the mutex
    inline void unlock() { CLIPPER_MUTEX_UNLK( &mutex ); }
  protected:
    CLIPPER_MUTEX_TYPE mutex;
  };


  //! Thread base class: Override this to create new threads
  /*! To create a thread, override this class. Store data as members
    with accessors to set input and read output. Override the Run()
    method to do the actual work. e.g. the following class implements
    a thread which can sum a list of numbers.
\code
class Thread_test : public Thread_base {
public:
  class Data : public std::vector<int> {
  public:
    Data() {}
    Data( std::vector<int> v ) : std::vector<int>(v) {}
  };

  void Run() {
    sum = 0;
    while ( 1 ) {
        lock();
        int c = current++;
        unlock();
        if ( c >= data_.size() ) break;
        sum += data_[c];
    }
  }
  static void set_data( Data d ) { data_ = d; }
  static Data data_;
  static int current;
  int sum;
};
\endcode
  */

  class CLIPPER_IMEX Thread_base {
  public:
    Thread_base();
    virtual ~Thread_base() {}
    bool run();
    bool join();
    static void lock()   { mutex_global.lock(); }
    static void unlock() { mutex_global.unlock(); }
    int id() const       { return id_; }

  protected:
    virtual void Run() = 0;

  private:
    static CLIPPER_THREAD_RETTYPE Entry( CLIPPER_THREAD_ARGTYPE thisptr );
    static Mutex mutex_global;
    static int   next_id;
    CLIPPER_THREAD_TYPE thread;
    int id_;
  };


#else


  class Mutex {
  public:
    Mutex()       {}
    ~Mutex()      {}
    inline void lock()   {}
    inline void unlock() {}
  };

  class Thread_base {
  public:
    Thread_base() {}
    virtual ~Thread_base() {}
    bool run()           { Run(); return true; }
    bool join()          { return true; }
    static void lock()   {}
    static void unlock() {}
    int id() const       { return -1; }
  protected:
    virtual void Run() { clipper::Message::message( clipper::Message_fatal( "No Run method defined" ) ); }
  };


#endif


} // namespace clipper

#endif
