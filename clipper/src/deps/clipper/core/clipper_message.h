/*! \file lib/clipper_message.h
    Header file for clipper message handler
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


#ifndef CLIPPER_MESSAGE
#define CLIPPER_MESSAGE


#include <string>
#include <iostream>
#include "../imex.h"

namespace clipper
{

  //! Message handler class
  /*! The message handler is a static class which handles messages and
    errors. It has 3 properties:
     - the output stream: to which messages will be directed (default stderr)
     - message level: messages with a level >= this will be output (default 5)
     - fatal level:   messages with a level >= this will be fatal (default 9)

     Levels may be in the range 1-9. They are priorotised as follows:
     - 1-4 messages, never fatal
     - 5-8 warning, may be fatal
     - 9: always fatal.
    The fatal level must be greater than or equal to the message
    level, and greater than or equal to 5.

    A message is any object which implements the following methods:
    \code
    const std::string& text() const;
    int level() const;
    \endcode
    The level method may be static.
    Messages are usually derived from Message_base. */
  class CLIPPER_IMEX Message {
  public:
    Message();  //!< null constuctor
    //! return the current stream
    static std::ostream& stream() { return *stream_; }
    //! return the current message level
    static const int& message_level() { return message_level_; }
    //! return the current fatal error level
    static const int& fatal_level() { return fatal_level_; }
    //! set the output stream
    static void set_stream( std::ostream& stream );
    //! set the current message level
    static void set_message_level( const int& level );
    //! set the current fatal error level
    static void set_fatal_level( const int& level );
    //! pass a message
    template<class T> inline static void message( const T& message )
    {
      if ( message.level() >= Message::message_level() ) {
	Message::stream() << message.text() << "\n";
	if ( message.level() >= Message::fatal_level() ) throw message;
      }
    }
  private:
    static int message_level_;
    static int fatal_level_;
    static std::ostream* stream_;
  };


  //! Base type for messages
  class CLIPPER_IMEX Message_base {
  public:
    const std::string& text() const;
    int level() const;
  protected:
    Message_base() {}
  };

  //! Generic message
  class CLIPPER_IMEX Message_generic : public Message_base {
  public:
    Message_generic( const int& level, const std::string& text ) :
      text_( text ), level_( level ) {}
    const std::string& text() const { return text_; }
    const int& level() const { return level_; }
  private:
    std::string text_;
    int level_;
  };

  //! Fatal message (level = 9)
  class CLIPPER_IMEX Message_fatal : public Message_base {
  public:
    Message_fatal( const std::string& text ) : text_( text ) {}
    const std::string& text() const { return text_; }
    static int level() { return 9; }
  private:
    std::string text_;
  };

  //! Warning message (level = 5)
  class CLIPPER_IMEX Message_warn : public Message_base  {
  public:
    Message_warn( const std::string& text ) : text_( text ) {}
    const std::string& text() const { return text_; }
    static int level() { return 5; }
  private:
    std::string text_;
  };

  //! Info message (level = 1)
  class CLIPPER_IMEX Message_info : public Message_base  {
  public:
    Message_info( const std::string& text ) : text_( text ) {}
    const std::string& text() const { return text_; }
    static int level() { return 1; }
  private:
    std::string text_;
  };

  //! Constructor message (level = 2)
  class CLIPPER_IMEX Message_ctor : public Message_base  {
  public:
    Message_ctor( const std::string& text ) : text_( "+"+text ) {}
    const std::string& text() const { return text_; }
    static int level() { return 2; }
  private:
    std::string text_;
  };

  //! Destructor message (level = 2)
  class CLIPPER_IMEX Message_dtor : public Message_base  {
  public:
    Message_dtor( const std::string& text ) : text_( "-"+text ) {}
    const std::string& text() const { return text_; }
    static int level() { return 2; }
  private:
    std::string text_;
  };


} // namespace clipper

#endif
