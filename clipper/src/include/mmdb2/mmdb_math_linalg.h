//  $Id: mmdb_math_linalg.h $
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
//  **** Module  :  LinAlg  <interface>
//       ~~~~~~~~~
//  **** Project :  MMDB  ( MacroMolecular Data Base )
//       ~~~~~~~~~
//
//  (C) E.Krissinel  2000-2013
//
//  =================================================================
//
//

#ifndef __MMDB_MATH_LinAlg__
#define __MMDB_MATH_LinAlg__

#include "mmdb_mattype.h"

namespace mmdb  {

  namespace math  {

    //  ==========================  Jacobi  =============================


    ///  Diagonalization of symmetric matrices A[1..N][1..N]
    /// by the method of Jacobi.
    extern void  Jacobi ( int     N,     //!< dimension of the matrix
                          rmatrix A,     //!< matrix to diagonalize; the
                                         /// lower triangle, except the
                                         /// diagonal, will remain unchanged
                          rmatrix T,     //!< eigenvectors placed as columns
                          rvector Eigen, //!< vector of eigenvalues, orderd
                                         /// by increasing
                          rvector Aik,   //!< working array
                          int &  Signal  //!< 0 <=> Ok, ItMax <=> iteration
                                         /// limit exchausted.
                        );


    //  A5.5.2  :  Perturbed Cholesky Decomposition
    extern void  PbCholDecomp ( int        N,
                                rvector    HDiag,
                                realtype   MaxOff,
                                realtype   MachEps,
                                rmatrix    L,
                                realtype & MaxAdd );

    //  A3.2.3a  :  Cholesky's   L - Solution  of
    //              L*Y  =  B  ( given  B )
    extern void  LSolve ( int N, rmatrix L, rvector B, rvector Y );

    //  A3.2.3b  :  Cholesky's   LT - Solution  of
    //              LT*X  =  Y  ( given  Y )
    extern void  LTSolve ( int N, rmatrix L, rvector Y, rvector X );

    //  A3.2.3   :  Solution of the equation    L*LT*S = G
    //              by the  Cholesky's  method
    extern void  ChSolve ( int N, rmatrix L, rvector G, rvector S );


    //  ----------------------------------------------------

    extern void  FastInverse (  int N, rmatrix A, ivector J0,
    //#D                          realtype &  Det,
                                int & Signal );
    //
    //      13.09.90  <--  Last Modification Date
    //                    ------------------------
    //
    // ================================================
    //
    //        Fast Inversion of the matrix  A
    //      by the method of  GAUSS - JORDAN  .
    //
    // ------------------------------------------------
    //
    //          Input  parameters  are  :
    //
    //     N   -   dimension of the matrix
    //     A   -   the matrix [1..N][1..N] to be inverted.
    // ------------------------------------------------
    //
    //     J0  -   integer vector [1..N] for temporal storage
    //
    // ------------------------------------------------
    //
    //          Output parameters  are  :
    //
    //     A   -   the inverted matrix
    //     Signal - the error key :
    //            = 0   <=>   O'K
    //             else
    //            degeneration was found, and
    //            the rang of matrix is  Signal-1.
    //
    //        Variable  Det  may return the determinant
    //     of matrix A.  To obtain it, remove all comments
    //     of form  //#D.
    //
    // ================================================


    //  ----------------------------------------------------

    extern void  SVD ( int    NA,  int     M,    int N,
                       rmatrix A,  rmatrix U,    rmatrix V,
                       rvector W,  rvector RV1,
                       bool MatU,  bool   MatV,
                       int & RetCode );
    //
    //      13.12.01  <--  Last Modification Date
    //                    ------------------------
    //
    // ================================================
    //
    //         The    Singular Value Decomposition
    //    of the matrix  A  by the algorithm from
    //      G.Forsait, M.Malkolm, K.Mouler.  Numerical
    //    methods of mathematical calculations //
    //    M., Mir, 1980.
    //
    //         Matrix  A  is represented as
    //
    //         A  =  U * W * VT
    //
    // ------------------------------------------------
    //
    //  All dimensions are indexed from 1 on.
    //
    // ------------------------------------------------
    //
    //         Input  parameters:
    //
    //     NA  -   number of lines in A. NA may be
    //           equal to M or N  only.  If NA=M
    //           then usual SVD will be made. If MA=N
    //           then matrix A is transposed before
    //           the decomposition, and the meaning of
    //           output parameters  U  and  V  is
    //           swapped (U accepts VT and VT accepts U).
    //           In other words, matrix  A  has physical
    //           dimension of  M x N , same as U and V;
    //           however the logical dimension of it
    //           remains that of  N x M .
    //     M   -   number of lines in  U
    //     N   -   number of columns in  U,V and length
    //           of  W,RV1 .  Always provide  M >= N  !
    //     A   -   matrix [1..M][1..N] or [1..N][1..M]
    //           to be decomposed. The matrix does not
    //           change,  and it may coincide with  U  or
    //           V, if NA=M (in which case  A  does change)
    //     MatU -  compute  U , if set True
    //     MatV -  compute  V , if set True
    //     RV1  -  temporary array [1..N].
    //     U    - should be always supplied as an array of
    //            [1..M][1..N], M>=N .
    //     V    - should be suuplied as an array of
    //            [1..N][1..N] if MatV is True .
    //
    // ------------------------------------------------
    //
    //          Output parameters  are  :
    //
    //     W   -   N non-ordered singular values,
    //           if  RetCode=0. If RetCode<>0, the
    //           RetCode+1 ... N -th values are still
    //           valid
    //     U   -   matrix of right singular vectors
    //           (arranged in columns),  corresponding
    //           to the singular values in  W,  if
    //           RetCode=0 and MatU is True.  If MatU
    //           is False, U  is still used as a
    //           temporary array. If RetCode<>0 then
    //           the  RetCode+1 ... N -th vectors
    //           are  valid
    //     V   -   matrix of left singular vectors
    //           (arranged in columns),  corresponding
    //           to the singular values in  W,  if
    //           RetCode=0 and MatV is True. If MatV
    //           is False, V is not used and may be set
    //           to NULL. If RetCode<>0 then the
    //           RetCode+1 ... N -th vectors are valid
    //     RetCode - the error key :
    //            = 0   <=>   O'K
    //              else
    //            = k, if the k-th singular value
    //                 was not computed after 30 iterations.
    //
    // ------------------------------------------------
    //
    //          Key  Variables  are  :
    //
    //     ItnLimit  -  the limit for iterations
    //
    //     This routine does not use any machine-dependent
    //  constants.
    //
    // ================================================
    //
    //

    extern void  OrderSVD ( int M, int N, rmatrix U, rmatrix V,
                            rvector W, bool MatU, bool MatV );


  }
}

#endif
