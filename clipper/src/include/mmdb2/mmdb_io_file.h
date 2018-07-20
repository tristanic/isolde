//  $Id: mmdb_io_file.h $
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
//    09.10.15   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  file_  <interface>
//       ~~~~~~~~~
//  **** Classes :  mmdb::io::File  - file I/O Support.
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2000-2015
//
//  =================================================================
//

#ifndef  __MMDB_IO_File__
#define  __MMDB_IO_File__

#include <stdio.h>

#include "mmdb_mattype.h"
#include "imex.h"
namespace mmdb  {

  namespace io  {

    //  ========================  File Class  ========================

    enum GZ_MODE  {
      GZM_NONE             = 0,
      GZM_CHECK            = 1,
      GZM_ENFORCE          = 2,
      GZM_ENFORCE_GZIP     = 2,
      GZM_ENFORCE_COMPRESS = 3
    };

    enum FILE_ERROR  {
      FileError_NoMemory              = 110,
      FileError_ShortData             = 111,
      FileError_NoDataFound           = 112,
      FileError_NoColumn              = 113,
      FileError_BadData               = 114,
      FileError_WrongMemoryAllocation = 115
    };

    enum SYSKEY  {
      syskey_unix = 1,
      syskey_win  = 2,
      syskey_all  = 3
    };

    extern const char _dir_sep_c;
    extern cpstr      _dir_sep;

    // ===================  Auxilary Functions  ========================


    extern cpstr GetFPath  ( pstr  FilePath, SYSKEY syskey=syskey_unix );
    extern cpstr GetFName  ( cpstr FilePath, SYSKEY syskey=syskey_unix );
    extern cpstr GetFExt   ( cpstr FilePath );
    extern cpstr ChangeExt ( pstr  FilePath, cpstr newExt,
                                             SYSKEY syskey=syskey_unix );
    extern cpstr FileError        ( int   ErrCode     );
    extern void  RemoveDelimiters ( pstr  S, int SLen );
    extern void  PickOutNumber    ( cpstr S, pstr SV, int SLen, int & j );


    // ==========================  File  ===============================

    DefineClass(File);

    class  MMDB_IMEX File  {

      public :

        File ( word BufSize=4096 );
        virtual  ~File();

        // ---- control functions
        //   FileName allows for "stdin", "stdout" and "stderr" as
        // for standard UNIX streams.
        void  assign       ( cpstr   FileName,
                             bool    Text=false,
                             bool    UniB=false,
                             GZ_MODE gzMode=GZM_NONE );
        //   assign for memory IO
        void  assign       ( word poolSize, word sizeInc, pstr filePool );
        void  takeFilePool ( pstr & filePool, word & fileSize );
        inline void GetFilePool ( pstr & filePool, word & fileSize )  {
          takeFilePool ( filePool,fileSize );
        }

        inline cpstr FileName() { return FName; }
        void  truncate    ( long size ); // call before reset/append
        bool  reset       ( bool ReadOnly=false, int retry=0 );
                                // = true if opened, each retry 1 sec sleep
        bool  erase       ();   // = true if erased
        bool  exists      ();   // = true if exists
        bool  parse       ( cpstr FileName    ); // true if filled
        bool  rename      ( cpstr NewFileName ); // true if renamed
        bool  rewrite     ();    // = true if opened
        bool  append      ();    // = true if opened
        bool  isOpen      ();
        long  Position    ();
        inline long  FileLength() { return FLength; }
        bool  seek        ( long Position );
        bool  FileEnd     ();
        inline bool  Success   () { return IOSuccess; }
        inline void  SetSuccess() { IOSuccess = true; }
        void  flush       ();
        void  shut        ();

        // ---- binary I/O
        word  ReadFile     ( void * Buffer, word Count );
        word  CreateRead   ( pstr & Line );
        word  ReadTerLine  ( pstr  Line, bool longLine=false );
        bool  WriteFile    ( const void * Buffer, word Count );
        bool  CreateWrite  ( cpstr Line );
        bool  WriteTerLine ( cpstr Line, bool longLine=false );

        //  machine-independent binary I/O
        bool  WriteReal   ( realtype * V );
        bool  WriteFloat  ( realtype * V );
        bool  WriteInt    ( int      * I );
        bool  WriteShort  ( short    * S );
        bool  WriteLong   ( long     * L );
        bool  WriteBool   ( bool     * B );
        bool  WriteByte   ( byte     * B );
        bool  WriteWord   ( word     * W );
        bool  ReadReal    ( realtype * V );
        bool  ReadFloat   ( realtype * V );
        bool  ReadInt     ( int      * I );
        bool  ReadShort   ( short    * S );
        bool  ReadLong    ( long     * L );
        bool  ReadBool    ( bool     * B );
        bool  ReadByte    ( byte     * B );
        bool  ReadWord    ( word     * B );
        bool  AddReal     ( realtype * V );
        bool  AddFloat    ( realtype * V );
        bool  AddInt      ( int      * I );
        bool  AddShort    ( short    * S );
        bool  AddLong     ( long     * L );
        bool  AddByte     ( byte     * B );
        bool  AddWord     ( word     * B );

        //  complex data binary I/O
        bool  WriteVector      ( rvector    V, int len,    int Shift );
        bool  WriteVector      ( ivector   iV, int len,    int Shift );
        bool  WriteVector      ( lvector   lV, int len,    int Shift );
        bool  WriteVector      ( bvector    B, int len,    int Shift );
        bool  ReadVector       ( rvector    V, int maxlen, int Shift );
        bool  ReadVector       ( ivector   iV, int maxlen, int Shift );
        bool  ReadVector       ( lvector   lV, int maxlen, int Shift );
        bool  ReadVector       ( bvector    B, int maxlen, int Shift );
        bool  CreateReadVector ( rvector &  V, int & len,  int Shift );
        bool  CreateReadVector ( ivector & iV, int & len,  int Shift );
        bool  CreateReadVector ( lvector & lV, int & len,  int Shift );
        bool  CreateReadVector ( bvector &  B, int & len,  int Shift );
        bool  CreateReadVector ( rvector &  V, int Shift );
        bool  CreateReadVector ( ivector & iV, int Shift );
        bool  CreateReadVector ( lvector & lV, int Shift );
        bool  CreateReadVector ( bvector &  B, int Shift );
        bool  WriteMatrix      ( rmatrix & A,  int N, int M,
                                    int  ShiftN,  int ShiftM );
        bool  CreateReadMatrix ( rmatrix & A,  int ShiftN, int ShiftM );
        bool  CreateReadMatrix ( rmatrix & A,  int & N, int & M,
                                    int ShiftN, int ShiftM );

        /// ---- text I/O
        bool  Write       ( cpstr  Line );     //!< writes without LF
        bool  Write       ( realtype V, int length=10 ); //!< w/o LF
        bool  Write       ( int     iV, int length=5  ); //!< w/o LF
        bool  WriteLine   ( cpstr  Line );     //!< writes and adds LF
        bool  LF          ();                  //!< just adds LF
        word  ReadLine    ( pstr   Line, word MaxLen=255 );
        word  ReadNonBlankLine ( pstr S, word MaxLen=255 );

        ///  complex data text I/O

        // writes with spaces and adds LF
        bool  WriteDataLine  ( realtype X, realtype Y, int length=10 );

        bool  WriteParameter ( cpstr S, realtype X, // writes parameter
                               int ParColumn=40,    // name S and value X
                               int length=10 );     // at column ParColumn
                                                    // and adds LF.

        bool  WriteParameters ( cpstr S, int n_X, // writes parameter
                                rvector X,      // name S and n_X values
                                int ParColumn=40, // X[0..n_X-1] at col
                                int length=10 );  // ParColumn, ads LF.

        bool  ReadParameter  ( pstr S, realtype & X, // reads parameter
                               int ParColumn=40 );   // name S and val X
        bool  ReadParameter  ( pstr S, int & X,
                               int ParColumn=40 );

        bool  ReadParameters ( pstr S, int & n_X,  // reads parameter
                               rvector X,          // name S, counts the
                               int MaxLen=255,     // of values n_X and
                               int ParColumn=40 ); // reads X[0..n_X-1].
                                                  // MaxLen gives sizeof(S)

        //   WriteColumns writes data stored in X, Y and Z in the form
        // of columns, adding a blank line in the end. If Z (or Z and Y)
        // are set to NULL, then only X and Y (or only X) are written.
        //   Shift corresponds to the begining of arrays' enumeration
        // X[Shift..Shift+len-1].
        bool  WriteColumns ( rvector X, rvector Y, rvector Z,
                             int len, int Shift, int MLength );
        bool  WriteColumns ( rvector X, rvector Y,
                             int len, int Shift, int MLength );

        //   ReadColumns reads data stored by WriteColumns. X, Y, and Z
        // must be allocated prior to call.
        //   xCol, yCol and zCol specify the order number of columns
        // (starting from 0) to be read into X, Y and Z, correspondingly.
        // If zCol (or zCol and yCol) < 0 then Z (or Z and Y) are not read.
        //   Shift corresponds to the begining of arrays' enumeration
        // X[Shift..Shift+len-1].
        //   Returns number of lines read.
        int   ReadColumns  ( int maxlen, rvector X, rvector Y, rvector Z,
                             int xCol, int yCol, int zCol, int Shift );
        int   ReadColumns  ( int maxlen, rvector X, rvector Y,
                             int xCol, int yCol, int Shift );

        //   CreateReadColumns reads data stored by WriteColumns. X, Y,
        // and Z must be set to NULL prior to call. They will be allocated
        // within the procedure.
        //   xCol, yCol and zCol specify the order number of columns
        // (starting from 0) to be read into X, Y and Z, correspondingly.
        // If zCol (or zCol and yCol) < 0 then Z (or Z and Y) are not read.
        //   Shift corresponds to the begining of arrays' enumeration
        // X[Shift..Shift+len-1].
        //   Returns number of lines read, errors are reported by
        // ErrorCode().
        int  CreateReadColumns ( rvector & X, rvector & Y, rvector & Z,
                                 int xCol, int yCol, int zCol, int Shift );
        int  CreateReadColumns ( rvector & X, rvector & Y,
                                 int xCol, int yCol, int Shift );

        // ---- miscellaneous
        realtype GetNumber ( cpstr S );
        FILE *   GetHandle () { return hFile; }

      protected :
        word    Buf_Size;
        bool    TextMode,UniBin;
        GZ_MODE gzipMode;
        pstr    IOBuf;
        word    BufCnt,BufLen,BufInc;
        FILE *  hFile;
        bool    EofFile;
        pstr    FName;
        long    FLength;
        bool    IOSuccess;
        int     ErrCode;

        void  FreeBuffer   ();
        void  _ReadColumns ( int & DLen, pstr S, int SLen,
                             rvector X, rvector Y, rvector Z,
                             int xCol, int yCol, int zCol, int Shift );

      private :
        int   gzipIO;
        bool  StdIO,memIO,ownBuf;

    };


    extern void SetGZIPPath     ( pstr gzipPath,     pstr ungzipPath     );
    extern void SetCompressPath ( pstr compressPath, pstr uncompressPath );

    extern bool FileExists      ( cpstr FileName, PFile f=NULL );

  }

}


#endif

