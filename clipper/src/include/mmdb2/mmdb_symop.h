//  $Id: mmdb_symop.h $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2013.
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
//    12.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :   MMDB_SymOp <interface>
//       ~~~~~~~~~
//  **** Project :   MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//
//  **** Classes :   mmdb::SymOp  ( symmetry operators )
//       ~~~~~~~~~   mmdb::SymOps ( container of symmetry operators )
//
//   (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#ifndef __MMDB_SymOp__
#define __MMDB_SymOp__

#include "mmdb_io_stream.h"
#include "mmdb_defs.h"

namespace mmdb  {

  //  ====================  SymOp  ========================

  DefineClass(SymOp);
  DefineStreamFunctions(SymOp);

  class SymOp : public io::Stream  {

    public :

      SymOp ();
      SymOp ( io::RPStream Object );
      ~SymOp();

      int  SetSymOp  ( cpstr XYZOperation );
      pstr GetSymOp  ();

      void Transform ( realtype & x, realtype & y, realtype & z );

      void GetTMatrix ( mat44 & TMatrix );  // copies T to TMatrix
      void SetTMatrix ( mat44 & TMatrix );  // copies TMatrix to T

      bool CompileOpTitle ( pstr S );  // makes XYZOp from matrix T
      bool CompileOpTitle ( pstr S, mat44 symMat, bool compare );
      void Print          ();          // prints operation and matrix

      void Copy  ( PSymOp symOp );

      void write ( io::RFile f );
      void read  ( io::RFile f );

    protected :
      pstr  XYZOp;
      mat44 T;

      void InitSymOp    ();
      void FreeMemory   ();
      int  GetOperation ( int n );

  };


  //  ====================  SymOps  ========================

  enum SYMOP_RC  {
    SYMOP_Ok                =  0,
    SYMOP_NoLibFile         = -1,
    SYMOP_UnknownSpaceGroup = -2,
    SYMOP_NoSymOps          = -3,
    SYMOP_WrongSyntax       = -4,
    SYMOP_NotAnOperation    = -5,
    SYMOP_ZeroDenominator   = -6
  };

  DefineClass(SymOps);
  DefineStreamFunctions(SymOps);

  class SymOps : public io::Stream  {

    public :

      SymOps ();
      SymOps ( io::RPStream Object );
      ~SymOps();

      virtual void FreeMemory();

      int  SetGroupSymopLib ( cpstr SpaceGroup,
                              cpstr symop_lib=NULL );
        // Space Group is taken from symop.lib. Return Code:
        // SYMOP_Ok <=> success

      int  SetGroup ( cpstr SpaceGroup,
                      cpstr syminfo_lib=NULL );
        // Space Group is taken from syminfo.lib. Return Code:
        // SYMOP_Ok <=> success

      void Reset           ();        // removes all symmetry operations
      virtual int AddSymOp ( cpstr XYZOperation ); // adds symmetry
                                                   // operation
      void PutGroupName    ( cpstr SpGroupName  );

      //  GetNofSymOps()  returns Nops -- the number of sym. operations
      int  GetNofSymOps ();
      pstr GetSymOp     ( int Nop );

      //  Transform(..) transforms the coordinates according to the
      // symmetry operation Nop. The return code is non-zero if
      // Nop is a wrong operation number (must range from 0 to Nops-1).
      int  Transform ( realtype & x, realtype & y, realtype & z,
                       int Nop );

      //  GetTMatrix(..) returns the coordinate transformation matrix
      // for the symmetry operation Nop. The return code is non-zero if
      // Nop is a wrong operation number (must range from 0 to Nops-1).
      int  GetTMatrix ( mat44 & TMatrix, int Nop );

      void Print ();

      virtual void Copy ( PSymOps symOps );

      void write ( io::RFile f );
      void read  ( io::RFile f );

    protected :

      pstr    SpGroup;
      int     Nops;
      PPSymOp symOp;

      void InitSymOps();

  };

}  // namespace mmdb

// extern void TestSymOps();

#endif

