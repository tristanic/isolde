//  $Id: mmdb_io_stream.h $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2008.
//
//    This library is free software: you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License version 3, modified in accordance with the provisions
//    of the license to address the requirements of UK law.
//
//    You should have received a copy of the modified GNU Lesser
//    General Public License along with this library. If not, copies
//    may be downloaded from http://www.ccp4.ac.uk/ccp4license.php
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//  =================================================================
//
//    11.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  Stream <interface>
//       ~~~~~~~~~
//  **** Classes :  mmdb::io::Stream ( Basic Stream Class )
//       ~~~~~~~~~
//
//   (C) E. Krissinel 1995-2013
//
//  =================================================================
//

#ifndef __MMDB_IO_Stream__
#define __MMDB_IO_Stream__

#include "mmdb_io_file.h"
#include "imex.h"
//  *******************************************************************

#ifndef __ClassMacros

# define __ClassMacros

 //  A Class definition macros
# define DefineClass(ClassName)             \
   class ClassName;                         \
   typedef ClassName    * P##ClassName;     \
   typedef ClassName    & R##ClassName;     \
   typedef P##ClassName * PP##ClassName;    \
   typedef P##ClassName & RP##ClassName;

 //  A Structure definition macros
# define DefineStructure(StructureName)             \
   struct StructureName;                            \
   typedef StructureName    * P##StructureName;     \
   typedef StructureName    & R##StructureName;     \
   typedef P##StructureName * PP##StructureName;    \
   typedef P##StructureName & RP##StructureName;

#endif


#define  DefineStreamFunctions(ClassName)                              \
  extern MMDB_IMEX void StreamWrite ( mmdb::io::RFile f, RP##ClassName Object ); \
  extern MMDB_IMEX void StreamRead  ( mmdb::io::RFile f, RP##ClassName Object );


#define  MakeStreamFunctions(ClassName)                                \
  void StreamWrite ( mmdb::io::RFile f, RP##ClassName Object )  {      \
    StreamWrite_ ( f,(mmdb::io::RPStream)Object );                     \
  }                                                                    \
  mmdb::io::PStream StreamInit##ClassName ( mmdb::io::RPStream Object ) { \
    return (mmdb::io::PStream)(new ClassName(Object));                 \
  }                                                                    \
  void StreamRead ( mmdb::io::RFile f, RP##ClassName Object )  {       \
    StreamRead_ ( f,(mmdb::io::RPStream)Object,StreamInit##ClassName );\
  }

#define  DefineFactoryFunctions(ClassName)                             \
  typedef P##ClassName      Make##ClassName();                         \
  typedef Make##ClassName * PMake##ClassName;                          \
  typedef P##ClassName  StreamMake##ClassName ( mmdb::io::RPStream Object ); \
  P##ClassName  new##ClassName ();                                     \
  P##ClassName  streamNew##ClassName ( mmdb::io::RPStream Object );    \
  typedef StreamMake##ClassName * PStreamMake##ClassName;              \
  extern MMDB_IMEX void SetMakers##ClassName ( void * defMk, void * streamMk );  \
  extern MMDB_IMEX void StreamWrite ( mmdb::io::RFile f, RP##ClassName Object ); \
  extern MMDB_IMEX void StreamRead  ( mmdb::io::RFile f, RP##ClassName Object );


#define  MakeFactoryFunctions(ClassName)                               \
  static PMake##ClassName       make##ClassName       = NULL;          \
  static PStreamMake##ClassName streamMake##ClassName = NULL;          \
  P##ClassName new##ClassName()  {                                     \
    if (make##ClassName)  return (*make##ClassName)();                 \
                    else  return new ClassName();                      \
  }                                                                    \
  P##ClassName streamNew##ClassName ( mmdb::io::RPStream Object )  {   \
    if (streamMake##ClassName)                                         \
          return (*streamMake##ClassName)(Object);                     \
    else  return new ClassName(Object);                                \
  }                                                                    \
  void SetMakers##ClassName ( void * defMk, void * streamMk )  {       \
    make##ClassName       = PMake##ClassName(defMk);                   \
    streamMake##ClassName = PStreamMake##ClassName(streamMk);          \
  }                                                                    \
  void StreamWrite ( mmdb::io::RFile f, RP##ClassName Object )  {      \
    StreamWrite_ ( f,(mmdb::io::RPStream)Object );                     \
  }                                                                    \
  mmdb::io::PStream StreamInit##ClassName ( mmdb::io::RPStream Object ) { \
    return (mmdb::io::PStream)(streamNew##ClassName(Object));          \
  }                                                                    \
  void StreamRead ( mmdb::io::RFile f, RP##ClassName Object )  {       \
    StreamRead_ ( f,(mmdb::io::RPStream)Object,StreamInit##ClassName ); \
  }

namespace mmdb  {

  namespace io  {

    //  ==========================  Stream  ===========================

    //     Each streamable class should be derived from Stream
    //  and have constructor Class(PStream & Object), which should
    //  initialize all memory of the class, and virtual functions
    //  read(..) and write(..) (see below). Constructor Class(PStream&)
    //  must not touch the Object variable. This constructor is used
    //  only once just before the read(..) function. It is assumed that
    //  read(..)/write(..) functions of the Class provide storage/reading
    //  of  all vital data. Function read(..) must read data in exactly
    //  the same way as function write(..) stores it.
    //     For using Class in streams, three following functions should
    //  be supplied:
    //
    //     1.
    //     void StreamWrite ( File & f, PClass & Object )  {
    //       StreamWrite ( f,(PStream)Object );
    //     }
    //
    //     2.
    //     PStream ClassInit ( PStream & Object )  {
    //       return (PStream)(new Class(Object));
    //     }
    //
    //     3.
    //     void StreamRead ( File & f, PClass & Object )  {
    //       StreamRead_ ( f,(PStream)Object,ClassInit );
    //     }
    //
    //    All these functions are automatically generated by macros
    //  DefineStreamFunctions(Class) -- in the header -- and
    //  MakeStreamFunctions(Class) -- in the implementation body. Note
    //  that macro DefineClass(Class) should always be issued for
    //  streamable classes prior to the stream-making macros. Then
    //  Class may be streamed using functions #1 and #3.
    //    StreamRead will return NULL for Object if it was not in
    //  the stream. If Object existed before calling StreamRead(..)
    //  but was not found in the stream, it will be disposed (NULL
    //  assigned).


    DefineClass(Stream);
    DefineStreamFunctions(Stream);

    class MMDB_IMEX Stream  {
      public :
        Stream            ()           {}
        Stream            ( RPStream ) {}
        virtual ~Stream   ()           {}
        virtual void read  ( RFile )   {}
        virtual void write ( RFile )   {}
    };


    typedef PStream InitStreamObject(RPStream Object);

    extern  void StreamRead_  ( RFile f, RPStream Object,
                                         InitStreamObject Init );

    extern  void StreamWrite_ ( RFile f, RPStream Object );


  }

}

#endif
