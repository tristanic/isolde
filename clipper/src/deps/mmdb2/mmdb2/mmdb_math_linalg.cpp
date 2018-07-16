//  $Id: mmdb_math_linalg.cpp $
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
//  **** Module  :  LinAlg  <implementation>
//       ~~~~~~~~~
//  **** Project :  MMDB  ( MacroMolecular Data Base )
//       ~~~~~~~~~
//
//  (C) E.Krissinel  2000-2013
//
//  =================================================================
//
//

#include <stdio.h>
#include <math.h>

#include "mmdb_math_linalg.h"

namespace mmdb  {

  namespace math  {


    //  ==========================  Jacobi  =============================


    void  Jacobi ( int     N,     // dimension of the matrix
                   rmatrix A,     // matrix to diagonalize; the lower
                                  // triangle, except the diagonal,
                                  // will remain unchanged
                   rmatrix T,     // eigenvectors placed as columns
                   rvector Eigen, // vector of eigenvalues, orderd
                                  // by increasing
                   rvector Aik,   // working array
                   int &   Signal // 0 <=> Ok, ItMax <=> iteration limit
                                  // exchausted.
                 )  {
    //  Diagonalization of symmetric matrices  by the method of Jacobi.
    //  Key  variables:
    //     ItMax  -  the maximum available number of iterations
    //     Eps1   -  is used in  SNA  and  CSA  calculations
    //     Eps2   -  is the level of the elimination of the
    //               non-diagonal matrix elements
    //     Eps3   -  the criterion to stop the iterations.
    //               The iterations stop if (1-Sigma1/Sigma2)<=Eps3
    //               where  Sigma1 is the dekart norm of the eigenvalues
    //               at the preceding iteration and  Sigma2 is
    //               the same for the current iteration.

    realtype Eps1,Eps2,Eps3;
    realtype Sigma1,Sigma2,OffDsq, SPQ,CSA,SNA,Q,P, HoldIK,HoldKI;
    int      ItMax;
    int      i,j,k,Iter;

      Eps1  = 6.0e-9;
      Eps2  = 9.0e-12;
      Eps3  = 1.0e-8;
      ItMax = 9999;

      Signal = 0;

      if (N<=1)  {
        T[1][1]  = 1.0;
        Eigen[1] = A[1][1];
        return;
      }

      for (i=1;i<=N;i++)  {
        for (j=1;j<=N;j++)
          T[i][j] = 0.0;
        T[i][i]  = 1.0;
        Eigen[i] = A[i][i];
      }

      Sigma1 = 0.0;
      OffDsq = 0.0;

      //  Sigma1  is the Dekart measure of the diagonal elements
      //  OffDsq  is the Dekart measure of the non-diagonal elements

      for (i=1;i<=N;i++)  {
        Sigma1 += A[i][i]*A[i][i];
        if (i<N)
          for (j=i+1;j<=N;j++)
            OffDsq += A[i][j]*A[i][j];
      }

      if (OffDsq<Eps2*Eps2)  return;

      //  S      = OffDsq*2.0+Sigma1;
      Iter   = 1;
      HoldIK = 1.0;

      while ((Iter<=ItMax) && (HoldIK>Eps3))  {

        for (i=1;i<N;i++)
          for (j=i+1;j<=N;j++)  {

            Q = fabs(A[i][i]-A[j][j]);

            if ((Q<=Eps1) || (fabs(A[i][j])>Eps2))  {
              if (Q>Eps1)  {
                P   = 2.0*A[i][j]*(Q/(A[i][i]-A[j][j]));
                SPQ = sqrt(P*P+Q*Q);
                CSA = sqrt((1.0+Q/SPQ)/2.0);
                SNA = P/(SPQ*CSA*2.0);
              } else  {
                CSA = sqrt(0.5);
                SNA = CSA;
              }
              for (k=1;k<=N;k++)  {
                HoldKI  = T[k][i];
                T[k][i] = HoldKI*CSA + T[k][j]*SNA;
                T[k][j] = HoldKI*SNA - T[k][j]*CSA;
              }

              for (k=i;k<=N;k++)
                if (k<=j)  {
                  Aik[k]  = A[i][k];
                  A[i][k] = CSA*Aik[k] + SNA*A[k][j];
                  if (k==j)  {
                    A[j][k] = SNA*Aik[k] - CSA*A[j][k];
                    Aik[j]  = SNA*Aik[i] - CSA*Aik[j];
                  }
                } else  {
                  HoldIK  = A[i][k];
                  A[i][k] = CSA*HoldIK + SNA*A[j][k];
                  A[j][k] = SNA*HoldIK - CSA*A[j][k];
                }

              for (k=1;k<=j;k++)
                if (k>i)
                  A[k][j] = SNA*Aik[k] - CSA*A[k][j];
                else  {
                  HoldKI  = A[k][i];
                  A[k][i] = CSA*HoldKI + SNA*A[k][j];
                  A[k][j] = SNA*HoldKI - CSA*A[k][j];
                }

            }

          }

        Sigma2 = 0.0;
        for (i=1;i<=N;i++)  {
          Eigen[i] = A[i][i];
          Sigma2  += Eigen[i]*Eigen[i];
        }

        HoldIK = fabs(1.0-Sigma1/Sigma2);
        Sigma1 = Sigma2;
        Iter++;

      }

      if (Iter>ItMax)  Signal = ItMax;

      for (i=1;i<=N;i++)  {
        k = i;
        for (j=i;j<=N;j++)
          if (Eigen[j]<Eigen[k])  k = j;
        if (k!=i)  {
          P = Eigen[k];
          Eigen[k] = Eigen[i];
          Eigen[i] = P;
          for (j=1;j<=N;j++)  {
            P = T[j][k];
            T[j][k] = T[j][i];
            T[j][i] = P;
          }
        }
      }


    }


    // -----------------------------------------------------

    void  PbCholDecomp ( int        N,
                         rvector    HDiag,
                         realtype   MaxOff,
                         realtype   MachEps,
                         rmatrix    L,
                         realtype & MaxAdd )  {

    //  A5.5.2  :  Perturbed Cholesky Decomposition
    // Part of the modular software system from
    // the appendix of the book "Numerical Methods for Unconstrained
    // Optimization and Nonlinear Equations" by Dennis & Schnabel 1983.

    int      i,j,k;
    realtype MinL,MinL2,S,MinLjj,MaxOffl, BB;

      MaxOffl = MaxOff;
      MinL    = sqrt(sqrt(MachEps))*MaxOffl;
      if (MaxOffl==0.0)  {
        for (i=1;i<=N;i++)  {
          BB = fabs(HDiag[i]);
          if (BB>MaxOffl)  MaxOffl = BB;
        }
        MaxOffl = sqrt(MaxOffl);
      }

      MinL2  = sqrt(MachEps)*MaxOffl;
      MaxAdd = 0.0;
      for (j=1;j<=N;j++)  {
        S = 0.0;
        if (j>1)
          for (i=1;i<j;i++)
            S += L[j][i]*L[j][i];
        L[j][j] = HDiag[j] - S;
        MinLjj  = 0.0;
        if (j<N)
          for (i=j+1;i<=N;i++)  {
            S = 0.0;
            if (j>1)
              for (k=1;k<j;k++)
                S += L[i][k]*L[j][k];
            L[i][j] = L[j][i] - S;
            BB = fabs(L[i][j]);
            if (BB>MinLjj)  MinLjj = BB;
          }
        BB = MinLjj/MaxOffl;
        if (BB>MinL)  MinLjj = BB;
                else  MinLjj = MinL;
        if (L[j][j]>MinLjj*MinLjj) L[j][j] = sqrt(L[j][j]);
        else  {
          if (MinL2>MinLjj)  MinLjj = MinL2;
          BB = MinLjj*MinLjj-L[j][j];
          if (BB>MaxAdd)  MaxAdd = BB;
          L[j][j] = MinLjj;
        }
        if (j<N)
          for (i=j+1;i<=N;i++)
            L[i][j] /= L[j][j];
      }

    }



    // -----------------------------------------------------

    void  LSolve ( int N, rmatrix L, rvector B, rvector Y )  {
    //  A3.2.3a  :  Cholesky's   L - Solution  of
    //              L*Y  =  B  ( given  B )
    int  i,j;
      Y[1] = B[1]/L[1][1];
      if (N>1)
        for (i=2;i<=N;i++)  {
          Y[i] = B[i];
          for (j=1;j<i;j++)
            Y[i] -= L[i][j]*Y[j];
          Y[i] /= L[i][i];
        }
    }


    // -----------------------------------------------------

    void  LTSolve ( int N, rmatrix L, rvector Y, rvector X )  {
    //  A3.2.3b  :   Cholesky's   LT - Solution  of
    //               LT*X  =  Y  ( given  Y )
    int  i,j;
      X[N] = Y[N]/L[N][N];
      if (N>1)
        for (i=N-1;i>=1;i--)  {
          X[i] = Y[i];
          for (j=i+1;j<=N;j++)
            X[i] -= L[j][i]*X[j];
          X[i] /= L[i][i];
        }
    }


    // -----------------------------------------------------

    void  ChSolve ( int N, rmatrix L, rvector G, rvector S )  {
    //  A3.2.3  :  Solution of the equation    L*LT*S = G
    //             by the  Cholesky's  method
    //int i;
      LSolve  ( N,L,G,S );
      LTSolve ( N,L,S,S );
    //  for (i=1;i<=N;i++)
    //    S[i] = -S[i];
    }



    //  ----------------------------------------------------

    void  FastInverse (  int N, rmatrix A, ivector J0,
    //#D                     realtype & Det,
                         int & Signal )  {
    //
    //      17.01.91  <--  Last Date of Modification.
    //                    ----------------------------
    //
    // ================================================
    //
    //        Fast Inversion of the matrix  A
    //      by the method of  GAUSS - JOIRDAN  .
    //
    // ------------------------------------------------
    //
    //          Input  parameters  are  :
    //
    //     N   -   dimension of the matrix
    //     A   -   the matrix [1..N][1..N] to be inverted.
    //
    // ------------------------------------------------
    //
    //     J0  -   integer vector [1..N] for temporal storage
    //
    //
    // ------------------------------------------------
    //
    //          Output parameters  are  :
    //
    //     A   -   the inverted matrix
    //     Signal - the error key :
    //            = 0   <=>   O'K
    //              else
    //           degeneration was found, and
    //           the rang of matrix is  Signal-1.
    //
    //        Variable  Det  may return the determinant
    //     of matrix A.  To obtain it, remove all comments
    //     of form  //#D .
    //
    // ------------------------------------------------
    //
    //          Key  Variables  are  :
    //
    //     Eps   -  is the level for the degeneration
    //              detection.  Keep in mind,  that
    //              this routine does not norm the
    //              matrix given,  and thus Eps1
    //              is the  ABSOLUTE  value.
    //
    // ================================================
    //

    realtype  Eps = 1.0e-16;

    int       i,j,k,i0;
    realtype  A0,B;
    rvector   Ai,Ai0;

      Signal = 0;
      if  (N<=1)   {
        if (fabs(A[1][1])<Eps)  {
          Signal = 1;
          return;
        }
        A[1][1] = 1.0/A[1][1];
    //#D   Det      = A[1][1];
        return;
      }

      if (N==2)  {
        A0 = A[1][1];
        B  = A0*A[2][2] - A[1][2]*A[2][1];
    //#D   Det = B;
        if (fabs(B)<Eps)  {
          Signal = 1;
          return;
        }
        A[1][1] = A[2][2]/B;
        A[2][2] = A0/B;
        B       = -B;
        A[1][2] /= B;
        A[2][1] /= B;
        return;
      }

      for (i=1;i<=N;i++)  {
        //  1. Finding of the leading element ( in A0 );
        //     i0  is the number of the leading string
        A0 = 0.0;
        i0 = 0;
        for (j=i;j<=N;j++)  {
          if (fabs(A[j][i])>A0)  {
            A0 = fabs(A[j][i]);
            i0 = j;
          }
        }
        if (A0<Eps)  {
          Signal = i;   //  Degeneration is found here
          return;
        }

        //  2. Swapping the string
        J0[i] = i0;
        B     = 1.0/A[i0][i];
        Ai    = A[i0];
        Ai0   = A[i];
        A[i]  = Ai;
        A[i0] = Ai0;
        for (j=1;j<=N;j++)
          Ai[j] = Ai[j]*B;
        Ai[i] = B;

        //  3. Substracting the strings
        for (j=1;j<=N;j++)
          if (i!=j)  {
            Ai0 = A[j];
            B   = Ai0[i];
            Ai0[i] = 0.0;
            for (k=1;k<=N;k++)
              Ai0[k] = Ai0[k] - B*Ai[k];
           }

    //#D   Det = Det/Ai[i];

      }

      //  4.  Back Swapping the columns
      for (i=N;i>=1;i--)  {
        j = J0[i];
        if (j!=i)  {
    //#D     Det = -Det;
          for (k=1;k<=N;k++)  {
            B       = A[k][i];
            A[k][i] = A[k][j];
            A[k][j] = B;
          }
        }
      }

      return;

    }    //  End of the procedure  FastInverse




    //  ----------------------------------------------------

    realtype Sign ( realtype A, realtype B )  {
      if (B>=0.0)  return  A;
             else  return -A;
    }

    realtype SrX2Y2 ( realtype X, realtype Y )  {
    realtype Ax,Ay;
      Ax = fabs(X);
      Ay = fabs(Y);
      if (Ay>Ax)   return  Ay*sqrt((X*X)/(Y*Y)+1.0);
      if (Ay==Ax)  return  Ax*sqrt(2.0);
      return  Ax*sqrt((Y*Y)/(X*X)+1.0);
    }


    //  ----------------------------------------------------

    void  SVD ( int    NA,  int       M,   int N,
                rmatrix A,  rmatrix   U,   rmatrix V,
                rvector W,  rvector RV1,
                bool MatU,  bool   MatV,
                int & RetCode )  {
    //
    //      13.12.01  <--  Last Modification Date
    //                    ------------------------
    //
    // ================================================
    //
    //         The    Singular Value Decomposition
    //    of the matrix  A  by the algorithm from
    //      G.Forsait, M.Malkolm, K.Mouler.  Numerical
    //    methods of mathematical calculations
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
    int       ItnLimit=300;
    int       i,j,k,l,i1,k1,l1,its,mn,ExitKey;
    realtype  C,G,F,X,S,H,Y,Z,Scale,ANorm,GG;

      l1      = 0;  // this is to keep compiler happy
      RetCode = 0;

      if (U!=A)  {
        if (NA==M)
          for (i=1;i<=M;i++)
            for (j=1;j<=N;j++)
              U[i][j] = A[i][j];
        else
          for (i=1;i<=M;i++)
            for (j=1;j<=N;j++)
              U[i][j] = A[j][i];
      }

      G     = 0.0;
      Scale = 0.0;
      ANorm = 0.0;

      for (i=1;i<=N;i++)  {
        l      = i+1;
        RV1[i] = Scale*G;
        G      = 0.0;
        S      = 0.0;
        Scale  = 0.0;
        if (i<=M)  {
          for (k=i;k<=M;k++)
            Scale += fabs(U[k][i]);
          if (Scale!=0.0)  {
            for (k=i;k<=M;k++)  {
              U[k][i] /= Scale;
              S       += U[k][i]*U[k][i];
            }
            F = U[i][i];
            G = -Sign(sqrt(S),F);
            H = F*G-S;
            U[i][i] = F-G;
            if (i!=N)
              for (j=l;j<=N;j++)  {
                S = 0.0;
                for (k=i;k<=M;k++)
                  S += U[k][i]*U[k][j];
                F = S/H;
                for (k=i;k<=M;k++)
                  U[k][j] += F*U[k][i];
              }
            for (k=i;k<=M;k++)
              U[k][i] *= Scale;
          }
        }

        W[i]  = Scale*G;
        G     = 0.0;
        S     = 0.0;
        Scale = 0.0;

        if ((i<=M) && (i!=N))  {
          for (k=l;k<=N;k++)
            Scale += fabs(U[i][k]);
          if (Scale!=0.0)  {
            for (k=l;k<=N;k++)  {
              U[i][k] /= Scale;
              S       += U[i][k]*U[i][k];
            }
            F = U[i][l];
            G = -Sign(sqrt(S),F);
            H = F*G-S;
            U[i][l] = F-G;
            for (k=l;k<=N;k++)
              RV1[k] = U[i][k]/H;
            if (i!=M)
              for (j=l;j<=M;j++)  {
                S = 0.0;
                for (k=l;k<=N;k++)
                  S += U[j][k]*U[i][k];
                for (k=l;k<=N;k++)
                  U[j][k] += S*RV1[k];
              }
            for (k=l;k<=N;k++)
              U[i][k] *= Scale;
          }
        }

        ANorm = RMax( ANorm,fabs(W[i])+fabs(RV1[i]) );

      }

      //   Accumulation of the right-hand transformations

      if  (MatV)
        for (i=N;i>=1;i--)  {
          if (i!=N)  {
            if (G!=0.0)  {
              for (j=l;j<=N;j++)
                V[j][i] = (U[i][j]/U[i][l]) / G;
              for (j=l;j<=N;j++)  {
                S = 0.0;
                for (k=l;k<=N;k++)
                  S += U[i][k]*V[k][j];
                for (k=l;k<=N;k++)
                  V[k][j] += S*V[k][i];
              }
            }
            for (j=l;j<=N;j++)  {
              V[i][j] = 0.0;
              V[j][i] = 0.0;
            }
          }

          V[i][i] = 1.0;
          G       = RV1[i];
          l       = i;

        }


      //   Accumulation of the left-hand transformations

      if (MatU)  {
        mn = N;
        if (M<N)  mn = M;

        for (i=mn;i>=1;i--)  {
          l = i+1;
          G = W[i];
          if (i!=N)
            for (j=l;j<=N;j++)
              U[i][j] = 0.0;
          if (G!=0.0)  {
            if (i!=mn)
              for (j=l;j<=N;j++)  {
                S = 0.0;
                for (k=l;k<=M;k++)
                  S += U[k][i]*U[k][j];
                F = (S/U[i][i]) / G;
                for (k=i;k<=M;k++)
                  U[k][j] += F*U[k][i];
              }
            for (j=i;j<=M;j++)
              U[j][i] /= G;
          } else
            for (j=i;j<=M;j++)
              U[j][i] = 0.0;

          U[i][i] += 1.0;

        }
      }

      //   Diagonalization of the two-diagonal form.

      for (k=N;k>=1;k--)  {
        k1  = k-1;
        its = 0;

        do  {
          ExitKey  = 0;
          l        = k+1;
          while ((ExitKey==0) && (l>1))  {
            l--;
            l1 = l-1;
            if (fabs(RV1[l])+ANorm==ANorm)   ExitKey=1;
            else if (l1>0)  {
              if (fabs(W[l1])+ANorm==ANorm)  ExitKey=2;
            }
          }

    //      if (ExitKey!=1)  {  <-- this is original statement
          if (ExitKey>1)  {  // <-- prevents from corruption due to l1<1.
                             // This is a rare case as RV1[1] should be
                             // always 0.0 . Apparently this logics is
                             // on the edge of float-point arithmetic,
                             // therefore extra precaution for the case
                             // of l1<1 was found necessary.
            C       = 0.0;
            S       = 1.0;
            ExitKey = 0;
            i       = l;
            while ((ExitKey==0) && (i<=k))  {
              F       =  S*RV1[i];
              RV1[i]  =  C*RV1[i];
              if (fabs(F)+ANorm==ANorm)  ExitKey = 1;
              else  {
                G = W[i];
                H = SrX2Y2(F,G);
                W[i] = H;
                C = G/H;
                S = -F/H;
                if (MatU)
                  for (j=1;j<=M;j++)  {
                    Y         =  U[j][l1];
                    Z         =  U[j][i];
                    U[j][l1]  =  Y*C+Z*S;
                    U[j][i]   =  -Y*S+Z*C;
                  }
                i++;
              }
            }
          }

          //    Convergence  Checking

          Z = W[k];
          if (l!=k)  {
            if (its>=ItnLimit)  {
              RetCode = k;
              return;
            }
            its++;
            X  =  W[l];
            Y  =  W[k1];
            G  =  RV1[k1];
            H  =  RV1[k];
            F  =  ((Y-Z)*(Y+Z) + (G-H)*(G+H)) / ( 2.0*H*Y );
            if (fabs(F)<=1.0)  GG = Sign(sqrt(F*F+1.0),F);
                         else  GG = F*sqrt(1.0+1.0/F/F);
            F  =  ((X-Z)*(X+Z) + H*(Y/(F+GG)-H)) / X;

            //   Next  QR - Transformation

            C  =  1.0;
            S  =  1.0;
            for (i1=l;i1<=k1;i1++)  {
              i = i1+1;
              G = RV1[i];
              Y = W[i];
              H = S*G;
              G = C*G;
              Z = SrX2Y2(F,H);
              RV1[i1] = Z;
              C = F/Z;
              S = H/Z;
              F = X*C+G*S;
              G = -X*S+G*C;
              H = Y*S;
              Y = Y*C;
              if (MatV)
                for (j=1;j<=N;j++)  {
                  X        = V[j][i1];
                  Z        = V[j][i];
                  V[j][i1] = X*C+Z*S;
                  V[j][i]  = -X*S+Z*C;
                }

              Z = SrX2Y2(F,H);
              W[i1] = Z;
              if (Z!=0.0)  {
                C = F/Z;
                S = H/Z;
              }
              F = C*G+S*Y;
              X = -S*G+C*Y;
              if (MatU)
                for (j=1;j<=M;j++)  {
                  Y        = U[j][i1];
                  Z        = U[j][i];
                  U[j][i1] = Y*C+Z*S;
                  U[j][i]  = -Y*S+Z*C;
                }

            }

            RV1[l] = 0.0;
            RV1[k] = F;
            W[k]   = X;

          } else if (Z<0.0)  {

            W[k] = -Z;
            if (MatV)
              for (j=1;j<=N;j++)
                V[j][k] = -V[j][k];
          }

        } while (l!=k);

      }

    }

    //  -----------------------------------------------------

    void  OrderSVD ( int M, int N, rmatrix U, rmatrix V,
                     rvector W, bool MatU, bool MatV )  {

    int       i,k,j;
    realtype  P;

      //  External loop of the re-ordering
      for (i=1;i<N;i++)  {
        k = i;
        P = W[i];

        //  Internal loop :  finding of the index of greatest
        //  singular value over the remaining ones.
        for (j=i+1;j<=N;j++)
          if (W[j]>P)  {
            k = j;
            P = W[j];
          }

        if (k!=i)  {
          //  Swapping the singular value
          W[k] = W[i];
          W[i] = P;
          //  Swapping the U's columns (  if  needed  )
          if (MatU)
            for (j=1;j<=M;j++)  {
              P       = U[j][i];
              U[j][i] = U[j][k];
              U[j][k] = P;
            }
          //  Swapping the V's columns ( if  needed )
          if (MatV)
            for (j=1;j<=N;j++)  {
              P       = V[j][i];
              V[j][i] = V[j][k];
              V[j][k] = P;
            }
        }

      }

    }


    /*

    #ifndef  __STDIO_H
    #include <stdio.h>
    #endif

    int main ( int argc, char ** argv, char ** env )  {
    //  Test Jacobi
    matrix   A,T,A1;
    vector   Eigen,Aik;
    realtype SR;
    int      N,i,j,k,Signal;

      N = 4;

      GetMatrixMemory ( A,N,N,1,1 );
      GetMatrixMemory ( T,N,N,1,1 );
      GetMatrixMemory ( A1,N,N,1,1 );
      GetVectorMemory ( Eigen,N,1 );
      GetVectorMemory ( Aik  ,N,1 );

      k = 1;
      for (i=1;i<=N;i++)
        for (j=i;j<=N;j++)  {
          A[i][j]  = k++;
          A[i][j] *= 1000.0;
          A[j][i]  = A[i][j];
        }

      printf ( "  INITIAL MATRIX:\n" );
      for (i=1;i<=N;i++)  {
        for (j=1;j<=N;j++)
          printf ( "  %10.4f",A[i][j] );
        printf ( "\n" );
      }

      Jacobi (  N,A,T,Eigen,Aik,Signal );

      printf ( "\n  EIGEN VALUES AND EIGEN VECTORS:\n" );
      for (i=1;i<=N;i++)  {
        printf ( "  %10.4f    ",Eigen[i] );
        for (j=1;j<=N;j++)
          printf ( "  %10.4f",T[j][i] );
        printf ( "\n" );
      }
      printf ( "\n       measure: " );
      for (i=1;i<=N;i++)  {
        SR = 0.0;
        for (j=1;j<=N;j++)
          SR += T[j][i]*T[j][i];
        printf ( "  %10.4f",sqrt(SR) );
      }
      printf ( "\n" );

      for (i=1;i<=N;i++)
        for (j=1;j<=N;j++)  {
          A1[i][j] = 0.0;
          for (k=1;k<=N;k++)
            A1[i][j] += T[i][k]*Eigen[k]*T[j][k];
        }

      printf ( "\n  RESTORED INITIAL MATRIX:\n" );
      for (i=1;i<=N;i++)  {
        for (j=1;j<=N;j++)
          printf ( "  %10.4f",A1[j][i] );
        printf ( "\n" );
      }

      FreeMatrixMemory ( A,N,1,1 );
      FreeMatrixMemory ( T,N,1,1 );
      FreeMatrixMemory ( A1,N,1,1 );
      FreeVectorMemory ( Eigen,1 );
      FreeVectorMemory ( Aik  ,1 );


    }

    */

  }

}
