//  $Id: mmdb_cryst.h $
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
//  **** Module  :  MMDB_Cryst <interface>
//       ~~~~~~~~~
//  **** Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::CrystContainer ( container for cryst. data )
//       ~~~~~~~~~  mmdb::NCSMatrix  ( non-cryst. symm. matrix class )
//                  mmdb::TVect      ( translation vector class      )
//                  mmdb::Cryst      ( MMDB cryst. section class     )
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#ifndef __MMDB_Cryst__
#define __MMDB_Cryst__

#include "mmdb_io_stream.h"
#include "mmdb_symop.h"
#include "mmdb_defs.h"
#include "mmdb_utils.h"
#include "imex.h"

namespace mmdb  {

  //  ====================  CrystContainer  ======================

  DefineClass(CrystContainer);
  DefineStreamFunctions(CrystContainer);

  class MMDB_IMEX CrystContainer : public ClassContainer  {

    public :

      CrystContainer () : ClassContainer() {}
      CrystContainer ( io::RPStream Object )
                        : ClassContainer ( Object ) {}
      ~CrystContainer() {}

      PContainerClass MakeContainerClass ( int ClassID );

      ERROR_CODE AddMTRIXLine ( cpstr S );

  };


  //  ==================  NCSMatrix  ========================

  enum NCSM_SET  {
    NCSMSET_Matrix1 = 0x00000001,
    NCSMSET_Matrix2 = 0x00000002,
    NCSMSET_Matrix3 = 0x00000004,
    NCSMSET_All     = 0x00000007
  };

  DefineClass(NCSMatrix);
  DefineStreamFunctions(NCSMatrix);

  class MMDB_IMEX NCSMatrix : public ContainerClass  {

    friend class Cryst;

    public :

      int   serNum;   // serial number
      mat33 m;        // non-crystallographic symmetry matrix
      vect3 v;        // translational part of ncs matrix
      int   iGiven;   // iGiven flag (see PDB format)

      NCSMatrix ();
      NCSMatrix ( cpstr S );
      NCSMatrix ( io::RPStream Object );
      ~NCSMatrix();

      bool       PDBASCIIDump1   ( io::RFile f );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      void       MakeCIF         ( mmcif::PData CIF, int N );
      ERROR_CODE GetCIF          ( mmcif::PData CIF, int & n );
      CLASS_ID   GetClassID      () { return ClassID_NCSMatrix; }

      void  SetNCSMatrix    ( int serialNum,
                              mat33 & ncs_m, vect3 & ncs_v,
                              int i_Given );

      void  Copy  ( PContainerClass NCSMatrix );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :
      word  WhatIsSet;  //    mask       field
                        //   0x0001    MTRIX1 was converted
                        //   0x0002    MTRIX2 was converted
                        //   0x0004    MTRIX3 was converted

      void  Init();

  };


  //  ==================  TVect  ========================

  DefineClass(TVect);
  DefineStreamFunctions(TVect);

  class MMDB_IMEX TVect : public ContainerClass  {

    public :

      int   serNum;   // serial number
      vect3 t;        // translation vector
      pstr  comment;  // comment

      TVect ();
      TVect ( cpstr S );
      TVect ( io::RPStream Object );
      ~TVect();

      void       PDBASCIIDump    ( pstr S, int N );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      void       MakeCIF         ( mmcif::PData CIF, int N );
      ERROR_CODE GetCIF          ( mmcif::PData CIF, int & n );
      CLASS_ID   GetClassID      () { return ClassID_TVect; }

      void  Copy  ( PContainerClass TVect );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void  Init();

  };


  //  =================  Cryst  =======================

  DefineClass(Cryst);
  DefineStreamFunctions(Cryst);

  // constants for the CellCheck field
  enum CELL_CHECK  {
    CCHK_Ok           = 0x00000000,
    CCHK_NoCell       = 0x00000001,
    CCHK_Error        = 0x00000002,
    CCHK_Disagreement = 0x00000004,
    CCHK_NoOrthCode   = 0x00000008,
    CCHK_Translations = 0x00000010,
    CCHK_Unchecked    = 0x00001000
  };

  // constants for the WhatIsSet field
  enum CELL_SET  {
    CSET_CellParams1  = 0x00000001,
    CSET_CellParams2  = 0x00000002,
    CSET_CellParams   = 0x00000003,
    CSET_SpaceGroup   = 0x00000004,
    CSET_ZValue       = 0x00000008,
    CSET_CrystCard    = 0x0000000F,
    CSET_OrigMatrix1  = 0x00000010,
    CSET_OrigMatrix2  = 0x00000020,
    CSET_OrigMatrix3  = 0x00000040,
    CSET_OrigMatrix   = 0x00000070,
    CSET_ScaleMatrix1 = 0x00000080,
    CSET_ScaleMatrix2 = 0x00000100,
    CSET_ScaleMatrix3 = 0x00000200,
    CSET_ScaleMatrix  = 0x00000380,
    CSET_Transforms   = 0x00000400,
    CSET_DummyCell    = 0x00001000
  };

  extern cpstr OrthCode[6];

  class MMDB_IMEX Cryst : public io::Stream  {

    friend class Channel;

    public :
      realtype  a,b,c;            // cell parameters
      realtype  alpha,beta,gamma; // cell parameters
      mat44     RO,RF;            // orthogonal-fractional recalculation
                                  //   matrices
      mat44     ROU,RFU;          // ort-frac recalc matrices for
                                  //   anisotr. t-fac
      mat633    RR;               // standard orthogonalizations
      realtype  Vol;              // cell volume
      int       NCode;            // code of orthogonalization matrix
      SymGroup  spaceGroup;       // group of space symmetry as read
                                  //    from file
      SymGroup  spaceGroupFix;    // actually used space group
      int       Z;                // Z-value

      mat33     o;                // orthogonal transformation matrix
      vect3     t;                // translation orthogonal vector
      mat33     s;                // scale matrix
      vect3     u;                // translation part of the scale matrix

      word      CellCheck;        // 0x0000 - Ok
                                  // 0x0001 - no cell stored
                                  // 0x0002 - some error in cell volume
                                  // 0x0004 - disagreement between
                                  //           cell and PDB
                                  // 0x0008 - no orth code derived
                                  // 0x0010 - translations also specified
                                  // 0x1000 - the check was not done
      word   WhatIsSet;           // indicator of the fields set
      bool   ignoreScalei;        // flag to ignore SCALEi cards
      bool   processSG;           // flag to process space group at file
                                  // read
      bool   fixSpaceGroup;       // flag to fix space group at file read

      Cryst ();
      Cryst ( io::RPStream Object );
      ~Cryst();

      void  FreeMemory();
      void  Reset     ();

      //   ConvertPDBString(..) interprets an ASCII PDB line and fills
      // the corresponding data fields. It returns zero if the line was
      // successfully converted, otherwise returns a non-negative value
      // of Error_XXXX.
      //   PDBString must be not shorter than 81 characters.
      ERROR_CODE ConvertPDBString ( pstr PDBString );

      //   RWBROOKReadPrintout() may be invoked after reading PDB file
      // for simulating the old RWBROOK messages and warnings
      void  RWBROOKReadPrintout();

      void  SetCell ( realtype cell_a,
                      realtype cell_b,
                      realtype cell_c,
                      realtype cell_alpha,
                      realtype cell_beta,
                      realtype cell_gamma,
                      int      OrthCode );
      void  PutCell ( realtype cell_a,
                      realtype cell_b,
                      realtype cell_c,
                      realtype cell_alpha,
                      realtype cell_beta,
                      realtype cell_gamma,
                      int      OrthCode );

      void  GetCell ( realtype & cell_a,
                      realtype & cell_b,
                      realtype & cell_c,
                      realtype & cell_alpha,
                      realtype & cell_beta,
                      realtype & cell_gamma,
                      realtype & vol );

      void  GetRCell ( realtype & cell_as,
                       realtype & cell_bs,
                       realtype & cell_cs,
                       realtype & cell_alphas,
                       realtype & cell_betas,
                       realtype & cell_gammas,
                       realtype & vols );

      void  SetSyminfoLib ( cpstr syminfoLib );
      pstr  GetSyminfoLib ();

      int   SetSpaceGroup ( cpstr spGroup );
      pstr  GetSpaceGroup ();
      pstr  GetSpaceGroupFix();

      //   CalcCoordTransforms() should be called once after all data
      // relevant to the crystallographic information, are read and
      // converted. Field CellCheck will then have bits set if there
      // are errors, e.g. bit CCHK_NoCell means that the coordinate
      // transformations cannot be performed.
      void  CalcCoordTransforms();

      // A PDB ASCII dump
      void  PDBASCIIDump ( io::RFile f );

      ERROR_CODE GetCIF  ( mmcif::PData CIF );
      void  MakeCIF ( mmcif::PData CIF );

      bool areMatrices();  // returns True if the orthogonal-to-
                              // fractional and fractional-to-orthogonal
                              // matrices are defined

      //   Frac2Orth(..) and Orth2Frac(..) transform between fractional
      // and orthogonal coordinates, if areMatrices() returns True.
      // If the transformation matrices were not set, the functions just
      // copy the coordinates.  Returns True if the transformation was
      // done; False return means that transformation matrices were not
      // calculated
      bool Frac2Orth (
                realtype x,    realtype y,    realtype z,
                realtype & xx, realtype & yy, realtype & zz );
      bool Orth2Frac (
                realtype x,    realtype y,    realtype z,
                realtype & xx, realtype & yy, realtype & zz );

      //   Below, F and T are transformation matrices in fractional and
      // orthogonal coordinates, respectively.
      bool Frac2Orth ( mat44 & F, mat44 & T );
      bool Orth2Frac ( mat44 & T, mat44 & F );


      //   Cryst2Orth(..) and Orth2Cryst(..) transform between fractional
      // and orthogonal anisotropic temperature factors, if areMatrices()
      // returns True. If the transformation matrices were not set, the
      // functions leave the factors unchanged.
      //   Vector U is composed as follows:
      //      U[0]=u11   U[1]=u22   U[2]=u33
      //      U[3]=u12   U[4]=u13   U[5]=u23
      // Returns True if the transformation was done; False retuen
      // means that transformation matrices were not calculated
      bool  Cryst2Orth ( rvector U );
      bool  Orth2Cryst ( rvector U );

      void  CalcOrthMatrices();  // calculates RR, AC, cella's and Vol

      bool  isNCSMatrix     ();
      bool  isScaleMatrix   ();
      bool  isCellParameters();

      int   GetNumberOfSymOps();
      pstr  GetSymOp ( int Nop );

      int   GetNumberOfNCSMatrices();
      int   GetNumberOfNCSMates   ();  // Returns the number of
                                       // NCS mates not given in
                                       // the file (iGiven==0)

      bool  GetNCSMatrix ( int NCSMatrixNo, mat33 & ncs_m,
                           vect3 & ncs_v );
      bool  GetNCSMatrix ( int NCSMatrixNo, mat44 & ncs_m,
                           int & iGiven ); // no=0..N-1
      int   AddNCSMatrix ( mat33 & ncs_m, vect3 & ncs_v, int iGiven );

      //  GetTMatrix(..) calculates and returns the coordinate
      //  transformation matrix, which converts orthogonal coordinates
      //  according to the symmetry operation number Nop and places
      //  them into unit cell shifted by cellshift_a a's, cellshift_b
      //  b's and cellshift_c c's.
      //
      //  Return 0 means everything's fine,
      //         1 there's no symmetry operation Nop defined
      //         2 fractionalizing/orthogonalizing matrices were not
      //           calculated
      //         3 cell parameters were not set up.
      int   GetTMatrix ( mat44 & TMatrix, int Nop,
                         int cellshift_a, int cellshift_b,
                         int cellshift_c, PSymOps symOpers=NULL );

      //  GetUCTMatrix(..) calculates and returns the coordinate
      //  transformation matrix, which converts orthogonal coordinates
      //  according to the symmetry operation Nop. Translation part
      //  of the matrix is  being chosen such that point (x,y,z) has
      //  least distance to the center of primary (333) unit cell,
      //  and then it is shifted by cellshift_a a's, cellshift_b b's
      //  and cellshift_c c's.
      //
      //  Return 0 means everything's fine,
      //         1 there's no symmetry operation Nop defined
      //         2 fractionalizing/orthogonalizing matrices were not
      //           calculated
      //         3 cell parameters were not set up.
      //
      int   GetUCTMatrix ( mat44 & TMatrix, int Nop,
                           realtype x, realtype y, realtype z,
                           int cellshift_a, int cellshift_b,
                           int cellshift_c, PSymOps symOpers=NULL );

      //  GetFractMatrix(..) calculates and returns the coordinate
      //  transformation matrix, which converts fractional coordinates
      //  according to the symmetry operation number Nop and places them
      //  into unit cell shifted by cellshift_a a's, cellshift_b b's and
      //  cellshift_c c's.
      //
      //  Return 0 means everything's fine,
      //         1 there's no symmetry operation Nop defined
      //         2 fractionalizing/orthogonalizing matrices were not
      //           calculated
      //         3 cell parameters were not set up.
      int   GetFractMatrix ( mat44 & TMatrix, int Nop,
                             int cellshift_a, int cellshift_b,
                             int cellshift_c, PSymOps symOpers=NULL );

      //  GetSymOpMatrix(..) returns the transformation matrix for
      //  Nop-th symmetry operator in the space group
      //
      //  Return 0 means everything's fine,
      //         1 there's no symmetry operation Nop defined
      //         2 fractionalizing/orthogonalizing matrices were not
      //           calculated
      //         3 cell parameters were not set up.
      //
      int   GetSymOpMatrix ( mat44 & TMatrix, int Nop );

      void  Copy  ( PCryst Cryst );

      void  write ( io::RFile f );    // writes header to PDB binary file
      void  read  ( io::RFile f );    // reads header from PDB binary file

    protected :

      CrystContainer ncsMatrix;      // non-cryst. symm. matrices
      CrystContainer tVect;          // translation vectors

      realtype  as,bs,cs;             // calculated 'cell parameters'
      realtype  alphas,betas,gammas;  // calculated 'cell parameters'
      realtype  AC[6];
      realtype  VolChk,VolErr;

      pstr      syminfo_lib;          // path to syminfo.lib
      SymOps    symOps;               // symmetry operations

      void  Init ( bool fullInit );
      int   FixSpaceGroup();

  };

  extern MMDB_IMEX cpstr getOrthCodeName ( int NCode );

}  // namespace mmdb

/*
extern void  TestCryst();  //  reads from 'in.cryst', writes into
                           //  'out.cryst' and 'abin.cryst'
*/


#endif

