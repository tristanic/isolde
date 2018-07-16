C
C     modlib.f: grab-bag of mathematical functions
C
C     This library is free software: you can redistribute it and/or
C     modify it under the terms of the GNU Lesser General Public License
C     version 3, modified in accordance with the provisions of the 
C     license to address the requirements of UK law.
C 
C     You should have received a copy of the modified GNU Lesser General 
C     Public License along with this library.  If not, copies may be 
C     downloaded from http://www.ccp4.ac.uk/ccp4license.php
C 
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU Lesser General Public License for more details.
C
C  cross.f        dot.f          ea06c.f      ea08c.f     ea09c.f
C  fa01as.f       fa01bs.f       fa01cs.f     fa01ds.f    fm02ad.f
C  icross.f       idot.f         iminv3       match.f     matmul.f
C  matmulnm.f     matmulgen.f
C  matvec.f       mc04b.f        minvn.f      minv3.f     ranmar.f
C  scalev.f       transp.f       unit.f       vdif.f      vset.f
C  vsum.f         zipin.f        zipout.f
C
C Routines added by pjx (June 2000):
C
C     GMPRD multiply two matrices (general)
C     MATMUL4 multiply two 4x4 matrices
C     MATMLN multiply two nxn matrices
C     DETMAT calculate the determinant of 3x3 matrix
C     MATVC4 matrix times vector (matrix is 4x4 array, treated as 3x3,
C            vector is length 3)
C     ML3MAT multiplies three matrices of any size
C     INV44 invert 4x4 matrix
C     MATMLI integer version of 3x3 matrix product MATMUL
C     MATMULTRANS multiply 3x3 matrix by transpose of another 3x3 matrix
C     IMATVEC integer version of MATVEC (post-multiply 3x3 matrix with a
C             vector)
C     TRANSFRM variant of MATVC4, same except that the input vector is
C              overwritten by the output vector
C
C Routines added by mdw (April 2003)
C
C     ZJVX    Compute zero of Bessel function of 1st kind Jv(x)
C     JDVX    Compute Bessel functions of 1st kind Jv(x) and their derivatives
C     JVX     Compute Bessel functions of 1st kind Jv(x)
C     GAMMA   Compute gamma function GA(x)
C     MSTA1   Determine the starting point for backward
C             recurrence such that the magnitude of
C             Jn(x) at that point is about 10^(-MP)
C     MSTA2   Determine the starting point for backward
C             recurrence such that all Jn(x) has MP
C             significant digits
C
C_BEGIN_GMPRD
C
      SUBROUTINE GMPRD(A,B,R,N,M,L)
C     =============================
C
C
C---- Ssp general matrix product
C
C   R(N,L) = A(N,M) * B(M,L)
C
C     .. Scalar Arguments ..
      INTEGER L,M,N
C     ..
C     .. Array Arguments ..
      REAL A(N*M),B(M*L),R(N*L)
C     ..
C     .. Local Scalars ..
      INTEGER I,IB,IK,IR,J,JI,K
C     ..
C
C_END_GMPRD
C
      IR = 0
      IK = -M
      DO 30 K = 1,L
        IK = IK + M
        DO 20 J = 1,N
          IR = IR + 1
          JI = J - N
          IB = IK
          R(IR) = 0.0
          DO 10 I = 1,M
            JI = JI + N
            IB = IB + 1
            R(IR) = A(JI)*B(IB) + R(IR)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      END
C
C_BEGIN_MATMUL4
C
C     =========================
      SUBROUTINE MATMUL4(A,B,C)
C     =========================
C
C     Multiply two 4x4 matrices, A=B*C
C
C
      IMPLICIT NONE
C
C     .. Array Arguments ..
      REAL A(4,4),B(4,4),C(4,4)
C     ..
C     .. Local Scalars ..
      REAL S
      INTEGER I,J,K
C     ..
C_END_MATMUL4
C
      DO 30 I = 1,4
        DO 20 J = 1,4
          S = 0
          DO 10 K = 1,4
            S = B(I,K)*C(K,J) + S
   10     CONTINUE
          A(I,J) = S
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
C
C_BEGIN_MATMLN
C
      subroutine matmln(n,a,b,c)
C     =========================
C
C Multiply two nxn matrices
C  a = b . c
C
      integer n
      real a(n,n),b(n,n),c(n,n)
      integer i,j,k
C
C_END_MATMLN
C
      do 1 i=1,n
         do 2 j=1,n
            a(j,i)=0.
            do 3 k=1,n
               a(j,i)= b(j,k)*c(k,i)+a(j,i)
 3          continue
 2       continue
 1    continue
      return
      end
C
C_BEGIN_DETMAT
C
C     ===========================
      SUBROUTINE DETMAT(RMAT,DET)
C     ============================
C
C---- Calculate determinant DET of 3x3 matrix RMAT
C
C
C     .. Scalar Arguments ..
      REAL DET
C     ..
C     .. Array Arguments ..
      REAL RMAT(3,3)
C     ..
C_END_DETMAT
C
      DET = RMAT(1,1)*RMAT(2,2)*RMAT(3,3) +
     +      RMAT(1,2)*RMAT(2,3)*RMAT(3,1) +
     +      RMAT(1,3)*RMAT(2,1)*RMAT(3,2) -
     +      RMAT(1,3)*RMAT(2,2)*RMAT(3,1) -
     +      RMAT(1,2)*RMAT(2,1)*RMAT(3,3) -
     +      RMAT(1,1)*RMAT(2,3)*RMAT(3,2)
C
      END
C
C
C_BEGIN_MATVC4
C
      SUBROUTINE MATVC4(A,R,B)
C     ========================
C
C Matrix x Vector
C       A(3)  =  R(4x4) . B(3)
C
      REAL A(3), R(4,4), B(3)
C
      REAL AA(4), BB(4)
C
      EXTERNAL GMPRD,VSET
C
C_END_MATVC4
C
      CALL VSET(BB,B)
      BB(4) = 1.0
      CALL GMPRD(R,BB,AA,4,4,1)
      CALL VSET(A,AA)
      RETURN
      END
C
C
C_BEGIN_ML3MAT
C
C 26-Nov-1988       J. W. Pflugrath           Cold Spring Harbor Laboratory
C    Edited to conform to Fortran 77.  Renamed from Multiply_3_matrices to
C    ML3MAT
C
C ==============================================================================
C
C! to multiply three matrices
C
      SUBROUTINE ML3MAT
C          ! input: 1st side of 1st matrix
     1   (P
C          ! input: first matrix
     2   ,A
C          ! input: 2nd side of 1st matrix & 1st side of 2nd matrix
     3   ,Q
C          ! input: second matrix
     4   ,B
C          ! input: 2nd side of 2nd matrix & 1st side of 3rd matrix
     5   ,R
C          ! input: third matrix
     6   ,C
C          ! input: 2nd side of 3rd matrix
     7   ,S
C          ! output: product matrix
     8   ,D)
C
CEE Multiplies three real matrices of any dimensions.  It is not optimised
C   for very large matrices.
C Multiply_3_matrices
C*** this routine is inefficient!
C Multiply_3_matrices
Created: 15-NOV-1985 D.J.Thomas,   MRC Laboratory of Molecular Biology,
C                                  Hills Road, Cambridge, CB2 2QH, England
C
C            ! loop counters
      INTEGER   I,J,K,L
C            ! loop limits
      INTEGER   P,Q,R,S
C           ! first input matrix
      REAL   A  (1:P,1:Q)
C           ! second input matrix
      REAL   B  (1:Q,1:R)
C           ! third input matrix
      REAL   C  (1:R,1:S)
C           ! output matrix
      REAL   D  (1:P,1:S)
C
C_END_ML3MAT
C
      DO 100 L = 1, S
         DO 100 I = 1, P
            D(I,L) = 0.0
            DO 100 K = 1, R
               DO 100 J = 1, Q
C
C                 ! accumulate product matrix D=ABC
C
100               D(I,L) = D(I,L)  +  A(I,J) * B(J,K) * C(K,L)
               CONTINUE
            CONTINUE
         CONTINUE
      CONTINUE
C Multiply_3_matrices
      RETURN
      END
C
C
C_BEGIN_INV44
C
C     ======================
      SUBROUTINE INV44(A,AI)
C     ======================
C
C SUBROUTINE TO INVERT 4*4 MATRICES FOR CONVERSION BETWEEN
C FRACTIONAL AND ORTHOGONAL AXES
C PARAMETERS
C
C           A (I)   4*4 MATRIX TO BE INVERTED
C          AI (O)   INVERSE MATRIX
C
C SPECIFICATION STATEMENTS
C ------------------------
C
C
C GET COFACTORS OF 'A' IN ARRAY 'C'
C ---------------------------------
C
      IMPLICIT NONE
C
C     .. Array Arguments ..
      REAL A(4,4),AI(4,4)
C     ..
C     .. Local Scalars ..
      REAL AM,D
      INTEGER I,I1,II,J,J1,JJ
C     ..
C     .. Local Arrays ..
      REAL C(4,4),X(3,3)
C     ..
C
C_END_INV44
C
      DO 40 II = 1,4
        DO 30 JJ = 1,4
          I = 0
          DO 20 I1 = 1,4
            IF (I1.EQ.II) GO TO 20
            I = I + 1
            J = 0
            DO 10 J1 = 1,4
              IF (J1.EQ.JJ) GO TO 10
              J = J + 1
              X(I,J) = A(I1,J1)
   10       CONTINUE
   20     CONTINUE
          AM = X(1,1)*X(2,2)*X(3,3) - X(1,1)*X(2,3)*X(3,2) +
     +         X(1,2)*X(2,3)*X(3,1) - X(1,2)*X(2,1)*X(3,3) +
     +         X(1,3)*X(2,1)*X(3,2) - X(1,3)*X(2,2)*X(3,1)
          C(II,JJ) = (-1)** (II+JJ)*AM
   30   CONTINUE
   40 CONTINUE
C
C---- calculate determinant
C
      D = 0
      DO 50 I = 1,4
        D = A(I,1)*C(I,1) + D
   50 CONTINUE
C
C---- get inverse matrix
C
      DO 70 I = 1,4
        DO 60 J = 1,4
          AI(I,J) = C(J,I)/D
   60   CONTINUE
   70 CONTINUE
      RETURN
      END
C
C
C_BEGIN_MATMLI
C
      SUBROUTINE MATMLI(A,B,C)
C     ========================
C Integer matrix multiply
      INTEGER A(3,3),B(3,3),C(3,3),I,J,K
C
C_END_MATMULI
C
      DO 1 I=1,3
      DO 2 J=1,3
      A(I,J)=0
      DO 3 K=1,3
      A(I,J)=A(I,J)+B(I,K)*C(K,J)
3     CONTINUE
  2     CONTINUE
 1      CONTINUE
      RETURN
      END
C
C
C_BEGIN_MATMULTRANS
C
      SUBROUTINE MATMULTrans(A,B,C)
C     =============================
C
C  A=B*C(transpose) for 3x3 matrices
C
      IMPLICIT NONE
C
C     ..Array arguments
      REAL A(3,3),B(3,3),C(3,3)
C     ..Local arrays
      REAL CT(3,3)
C
      EXTERNAL MATMUL
C
C_END_MATMULTRANS
C
      CALL TRANSP(CT,C)
      CALL MATMUL(A,B,CT)
      RETURN
      END
C
C
C_BEGIN_IMATVEC
C
      SUBROUTINE IMATVEC(V,A,B)
C     ========================
C
C---- Post-multiply a 3x3 matrix by a vector
C
C       V=AB
C
C     .. Array Arguments ..
      INTEGER A(3,3),B(3),V(3)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,S
C
C_END_IMATVEC
C
      DO 20 I = 1,3
        S = 0
C
C
        DO 10 J = 1,3
          S = A(I,J)*B(J) + S
   10   CONTINUE
C
C
        V(I) = S
   20 CONTINUE
C
C
      END
C
C
C_BEGIN_TRANSFRM
C
      SUBROUTINE TRANSFRM(X,MAT)
C     ==========================
C
C  Transform vector X(3) by quine matrix MAT(4,4)
C  Return transformed vector in X.
C
      IMPLICIT NONE
C
C     ..Array arguments..
      REAL X(3),MAT(4,4)
C     ..Local arrays..
      REAL TMP(3)
C
C_END_TRANSFRM
C
      CALL MATVC4(TMP,MAT,X)
      CALL VSET(X,TMP)
      RETURN
      END
C
C
C
C_BEGIN_CROSS
C
      SUBROUTINE CROSS(A,B,C)
C     =======================
C
C     compute vector product A = B x C
C
C     .. Array Arguments ..
      REAL             A(3),B(3),C(3)
C
C_END_CROSS
C     ..
      A(1) = B(2)*C(3) - C(2)*B(3)
      A(2) = B(3)*C(1) - C(3)*B(1)
      A(3) = B(1)*C(2) - C(1)*B(2)
      END
C
C
C_BEGIN_DOT
C
      REAL FUNCTION DOT(A,B)
C     ======================
C
C     dot product of two vectors
C
C     .. Array Arguments ..
      REAL              A(3),B(3)
C
C_END_DOT
C     ..
      DOT = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
      END
C

C     ******************************************************************
      SUBROUTINE EIGEN_RS_ASC(A, R, N, MV)
C     ******************************************************************
C
C---- SUBROUTINE TO COMPUTE EIGENVALUES & EIGENVECTORS OF A REAL
C---- SYMMETRIC MATRIX, FROM IBM SSP MANUAL (SEE P165).
C---- DESCRIPTION OF PARAMETERS -
C---- A - ORIGINAL MATRIX STORED COLUMNWISE AS UPPER TRIANGLE ONLY,
C---- I.E. "STORAGE MODE" = 1.  EIGENVALUES ARE WRITTEN INTO DIAGONAL
C---- ELEMENTS OF A  I.E.  A(1)  A(3)  A(6)  FOR A 3*3 MATRIX.
C---- R - RESULTANT MATRIX OF EIGENVECTORS STORED COLUMNWISE IN SAME
C---- ORDER AS EIGENVALUES.
C---- N - ORDER OF MATRICES A & R.
C---- MV = 0 TO COMPUTE EIGENVALUES & EIGENVECTORS.
C
C
c      IMPLICIT    NONE

      INTEGER NMAX
      PARAMETER (NMAX=10)

      REAL        A(*), R(*)
      INTEGER     N, MV
C
      INTEGER     IQ, J, I, IJ, IA, IND, L, M, MQ, LQ, LM, LL, MM
      INTEGER     ILQ, IMQ, IM, IL, ILR, IMR
      REAL        ANORM, ANRMX, RANGE, THR, X, Y, SINX, SINX2,
     &            COSX, COSX2, SINCS

C Lapack variables
      LOGICAL LAPACK
      CHARACTER*1        JOBZ, LAPRANGE, UPLO
      INTEGER            INFO, NVECTORS
      REAL               ABSTOL
      INTEGER            ISUPPZ(2*NMAX), IWORK(10*NMAX)
      REAL               WORK(26*NMAX),EVALUES(NMAX),AM(NMAX,NMAX)
C
C-- FOR REAL
C      DATA RANGE/1D-12/
      DATA RANGE/1E-6/

C  Alternative lapack routine - only marginally tested
      LAPACK = .FALSE.
      IF (LAPACK) THEN
       NVECTORS = 0
       IF (N.GT.NMAX) 
     +    CALL CCPERR(1,'s/r EIGEN_RS_ASC: redimension NMAX!')
       IF (MV.EQ.0) JOBZ = 'V'
       LAPRANGE = 'A'
       UPLO = 'U'
       ABSTOL = 0.0
       IA = 0
       DO I = 1,N
         DO J = 1,I
           IA = IA + 1
           AM(J,I) = A(IA)
         ENDDO
       ENDDO
C       CALL SSYEVR(JOBZ, LAPRANGE, UPLO, N, AM, N, 
C     +   1, N, 1, N, ABSTOL, NVECTORS, EVALUES, R, 
C     +   N, ISUPPZ, WORK, 26*N, IWORK, 10*N, INFO)
       IA = 0
       DO I = 1,NVECTORS
         DO J = 1,I
           IA = IA + 1
           IF (J.EQ.I) A(IA) = EVALUES(I)
         ENDDO
       ENDDO
       RETURN
      ENDIF

      IF (MV.EQ.0) THEN
        IQ=-N
        DO    J=1,N
          IQ=IQ+N
          DO     I=1,N
            IJ=IQ+I
            IF (I.EQ.J) THEN
              R(IJ)=1.
            ELSE
              R(IJ)=0.
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
C---- INITIAL AND FINAL NORMS (ANORM & ANRMX)
      IA=0
      ANORM=0.
      DO    I=1,N
        DO    J=1,I
          IA=IA+1
          IF (J.NE.I) ANORM=ANORM+A(IA)**2
        ENDDO
      ENDDO
C
      IF (ANORM.LE.0.) GOTO 165
      ANORM=SQRT(2.*ANORM)
      ANRMX=ANORM*RANGE/N
C
C---- INITIALIZE INDICATORS AND COMPUTE THRESHOLD
      IND=0
      THR=ANORM
45    THR=THR/N
50    L=1
55    LQ=L*(L-1)/2
      LL=L+LQ
      M=L+1
      ILQ=N*(L-1)
C
C---- COMPUTE SIN & COS
60    MQ=M*(M-1)/2
      LM=L+MQ
      IF (A(LM)*A(LM)-THR.LT.0.) GOTO 130
      IND=1
      MM=M+MQ
      X=.5*(A(LL)-A(MM))
      Y=-A(LM)/SQRT(A(LM)**2+X*X)
C---- Protect against rounding error. C. Flensburg 20080307.
      IF (ABS(Y).GT.1.0) Y=SIGN(1.0,Y)
      IF (X.LT.0.) Y=-Y
      SINX=Y/SQRT(2.*(1.+(SQRT(1.-Y*Y))))
      SINX2=SINX**2
      COSX=SQRT(1.-SINX2)
      COSX2=COSX**2
      SINCS=SINX*COSX
C
C---- ROTATE L & M COLUMNS
      IMQ=N*(M-1)
      DO 125 I=1,N
        IQ=I*(I-1)/2
        IF (I.NE.L .AND. I.NE.M) THEN
          IF (I.LT.M) THEN
            IM=I+MQ
          ELSE
            IM=M+IQ
          ENDIF
          IF (I.LT.L) THEN
            IL=I+LQ
          ELSE
            IL=L+IQ
          ENDIF
          X=A(IL)*COSX-A(IM)*SINX
          A(IM)=A(IL)*SINX+A(IM)*COSX
          A(IL)=X
        ENDIF
        IF (MV.EQ.0) THEN
          ILR=ILQ+I
          IMR=IMQ+I
          X=R(ILR)*COSX-R(IMR)*SINX
          R(IMR)=R(ILR)*SINX+R(IMR)*COSX
          R(ILR)=X
        ENDIF
125   CONTINUE
C
      X=2.*A(LM)*SINCS
      Y=A(LL)*COSX2+A(MM)*SINX2-X
      X=A(LL)*SINX2+A(MM)*COSX2+X
      A(LM)=(A(LL)-A(MM))*SINCS+A(LM)*(COSX2-SINX2)
      A(LL)=Y
      A(MM)=X
C
C---- TESTS FOR COMPLETION
C---- TEST FOR M = LAST COLUMN
130   IF (M.NE.N) THEN
        M=M+1
        GOTO 60
      ENDIF
C
C---- TEST FOR L =PENULTIMATE COLUMN
      IF (L.NE.N-1) THEN
        L=L+1
        GOTO55
      ENDIF
      IF (IND.EQ.1) THEN
        IND=0
        GOTO50
      ENDIF
C
C---- COMPARE THRESHOLD WITH FINAL NORM
      IF (THR.GT.ANRMX) GOTO 45
165   RETURN
      END
C
C     The routines ea06c, ea08c, ea09c, fa01as, fa01bs, fa01cs, fa01ds,
C     fm02ad, mc04b, (and possibly others) are from the
C     Harwell Subroutine library.  The conditions on their external use,
C     reproduced from the Harwell manual are:
C     * due acknowledgement is made of the use of subroutines in any
C     research publications resulting from their use.
C     * the subroutines may be modified for use in research applications
C     by external users.  The nature of such modifiactions should be
C     indicated in writing for information to the liaison officer.  At
C     no time however, shall the subroutines or modifications thereof
C     become the property of the external user.
C       The liaison officer for the library's external affairs is listed
C     as: Mr. S. Marlow, Building 8.9, Harwell Laboratory, Didcot, 
C     Oxon OX11 0RA, UK.
C
C_BEGIN_EA06C
C
      SUBROUTINE EA06C(A,VALUE,VECTOR,M,IA,IV,W)
C     ==========================================
C
C**  18/03/70 LAST LIBRARY UPDATE
C
C       ( Calls EA08C(W,W(M1),VALUE,VECTOR,M,IV,W(M+M1))
C       and MC04B(A,W,W(M1),M,IA,W(M+M1)) )
C
C       Given a real MxM symmetric matrix A = {aij} this routine
C       finds all its eigenvalues (lambda)i i=1,2,.....,m  and
C       eigenvectors xj i=1,2,...,m.  i.e. finds the non-trivial
C       solutions of Ax=(lambda)x
C
C       The matrix is reduced to tri-diagonal form by applying
C       Householder transformations. The eigenvalue problem for
C       the reduced problem is then solved using the QR algorithm
C       by calling EA08C.
C
C  Argument list
C  -------------
C
C   IA    (I) (integer) should be set to the first dimension
C			of the array A, i.e. if the allocation
C			for the array A was specified by
C			DIMENSION A(100,50)
C			then IA would be set to 100
C
C   M     (I) (integer) should be set to the order m of the matrix
C
C   IV    (I) (integer) should be set to the first dimension
C			of the 2-dimensional array VECTOR
C
C   VECTOR(IV,M) (O) (real) 2-dimensional array, with first dimension IV,
C		            containing the eigenvectors. The components
C		            of the eigenvector vector(i) corresponding
C		            to the eigenvalue (lambda)i (in VALUE(I))
C		            are placed in VECTOR(J,I) J=1,2,...,M.
C		            The eigenvectors are normalized so that
C		            xT(i)x(i)=1 i=1,2,...,m.
C
C   VALUE(M) (O) (real)  array in which the routine puts
C		         the eigenvalues (lambda)i, i=1,2,...,m.
C		         These are not necessarily in any order.
C
C   W    (I) (real(*)) working array used by the routine for
C		       work space. The dimension must be set
C		       to at least 5*M.
C
C_END_EA06C
C
C
C     .. Scalar Arguments ..
      INTEGER          IA,IV,M
C     ..
C     .. Array Arguments ..
      REAL             A(IA,M),VALUE(M),VECTOR(IV,M),W(*)
C     ..
C     .. Local Scalars ..
      REAL             PP
      INTEGER          I,I1,II,K,L,M1
C     ..
C     .. External Subroutines ..
      EXTERNAL         EA08C,MC04B
C     ..
      M1 = M + 1
      W(1) = A(1,1)
      IF (M-2) 30,10,20
   10 W(2) = A(2,2)
      W(4) = A(2,1)
      GO TO 30

   20 CALL MC04B(A,W,W(M1),M,IA,W(M+M1))
   30 CALL EA08C(W,W(M1),VALUE,VECTOR,M,IV,W(M+M1))
      IF (M.GT.2) THEN
          DO 70 L = 1,M
              DO 60 II = 3,M
                  I = M - II + 1
                  IF (W(M1+I).NE.0) THEN
                      PP = 0.0
                      I1 = I + 1
                      DO 40 K = I1,M
                          PP = A(I,K)*VECTOR(K,L) + PP
   40                     CONTINUE
                      PP = PP/ (A(I,I+1)*W(M1+I))
                      DO 50 K = I1,M
                          VECTOR(K,L) = A(I,K)*PP + VECTOR(K,L)
   50                     CONTINUE
                  END IF

   60             CONTINUE
   70         CONTINUE
      END IF

      END
C
C_BEGIN_EA08C
C
      SUBROUTINE EA08C(A,B,VALUE,VEC,M,IV,W)
C     =====================================
C
C       (Calls EA09C(A,B,W(M+1),M,W))
C
C       This uses QR iteration to find the all the eigenvalues and
C       eigenvectors of the real symmetric tri-diagonal matrix
C       whose diagonal elements are A(i), i=1,M and off-diagonal
C       elements are B(i),i=2,M. The eigenvalues will have unit
C       length.	The array W is used for workspace and must have
C       dimension at least 2*M.	We treat VEC as	if it had
C       dimensions (IV,M).
C
C       First EA09, which uses the QR algorithm, is used to find
C       the eigenvalues; using these as shifts the QR algorithm is
C       again applied but now using the	plane rotations	to generate
C       the eigenvectors. Finally the eigenvalues are refined
C       by taking Rayleigh quotients of the vectors.
C
C    Argument list
C    -------------
C
C       A(M)	(I) (real)    Diagonal elements
C
C       B(M)	(I) (real)    Off-diagonal elements
C
C       IV	(I)  (integer)	should be set to the first dimen-
C                               sion of the 2-dimensional array VEC
C
C       M	(I) (integer) should be	set to the order m of the
C                             matrix
C
C       VALUE(M) (O) (real)   Eigenvalues
C
C       VEC	(O) (real)    Eigenvectors. The	dimensions
C			      should be	set to (IV,M).
C
C       W(*)	(I) (real)    Working array.The	dimension must be
C			      set to at	least 2*M.
C
C_END_EA08C
C
C     .. Scalar Arguments ..
      INTEGER          IV,M
C     ..
C     .. Array Arguments ..
      REAL             A(M),B(M),VALUE(M),VEC(1),W(*)
C     ..
C     .. Local Scalars ..
      REAL             A11,A12,A13,A21,A22,A23,A33,A34,CO,EPS,ROOT,S,SI,
     +                 V1,V2,XAX,XX
      INTEGER          I,II,ITER,J,J1,J2,JI,JK,K,L,N1,N2,N2M1
C     ..
C     .. External Subroutines ..
      EXTERNAL         EA09C
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC        ABS,MIN,SIGN,SQRT
C     ..
C     .. Data statements ..
      DATA             EPS/1.0E-5/,A34/0.0/
C     ..
C
C
      CALL EA09C(A,B,W(M+1),M,W)
C
C---- Set vec to the identity matrix.
C
      DO 20 I = 1,M
          VALUE(I) = A(I)
          W(I) = B(I)
          K = (I-1)*IV + 1
          L = K + M - 1
          DO 10 J = K,L
              VEC(J) = 0.0
   10         CONTINUE
          VEC(K+I-1) = 1.0
   20     CONTINUE
      ITER = 0
      IF (M.NE.1) THEN
          N2 = M
   30     CONTINUE
C
C---- Each qr iteration is performed of rows and columns n1 to n2
C
          DO 40 II = 2,N2
              N1 = 2 + N2 - II
              IF (ABS(W(N1)).LE. (ABS(VALUE(N1-1))+ABS(VALUE(N1)))*
     +            EPS) GO TO 50
   40         CONTINUE
          N1 = 1
   50     IF (N2.NE.N1) THEN
              ROOT = W(M+N2)
              ITER = ITER + 1
              N2M1 = N2 - 1
              A22 = VALUE(N1)
              A12 = A22 - ROOT
              A23 = W(N1+1)
              A13 = A23
              DO 70 I = N1,N2M1
                  A33 = VALUE(I+1)
                  IF (I.NE.N2M1) A34 = W(I+2)
                  S = SIGN(SQRT(A12*A12+A13*A13),A12)
                  SI = A13/S
                  CO = A12/S
                  JK = I*IV + 1
                  J1 = JK - IV
                  J2 = MIN(M,I+ITER) + J1 - 1
                  DO 60 JI = J1,J2
                      V1 = VEC(JI)
                      V2 = VEC(JK)
                      VEC(JI) = V1*CO + V2*SI
                      VEC(JK) = V2*CO - V1*SI
                      JK = JK + 1
   60                 CONTINUE
                  IF (I.NE.N1) W(I) = S
                  A11 = CO*A22 + SI*A23
                  A12 = CO*A23 + SI*A33
                  A13 = SI*A34
                  A21 = CO*A23 - SI*A22
                  A22 = CO*A33 - SI*A23
                  A23 = CO*A34
                  VALUE(I) = A11*CO + A12*SI
                  A12 = -A11*SI + A12*CO
                  W(I+1) = A12
                  A22 = A22*CO - A21*SI
   70             CONTINUE
              VALUE(N2) = A22
              GO TO 30

          ELSE
              N2 = N2 - 1
              IF (N2-1.GT.0) GO TO 30
          END IF
C
C---- Rayleigh quotient
C
          DO 90 J = 1,M
              K = (J-1)*IV
              XX = VEC(K+1)**2
              XAX = A(1)*XX
              DO 80 I = 2,M
                  XX = VEC(K+I)**2 + XX
                  XAX = (B(I)*2.0*VEC(K+I-1)+VEC(K+I)*A(I))*VEC(K+I) +
     +                  XAX
   80             CONTINUE
              VALUE(J) = XAX/XX
   90         CONTINUE
      END IF

      END
C
C
      SUBROUTINE EA09C(A,B,VALUE,M,OFF)
C     =================================
C
C    18/03/70 LAST LIBRARY UPDATE
C
C     .. Scalar Arguments ..
      INTEGER          M
C     ..
C     .. Array Arguments ..
      REAL             A(M),B(M),OFF(M),VALUE(M)
C     ..
C     .. Local Scalars ..
      REAL             A11,A12,A13,A21,A22,A23,A33,A34,BB,CC,CO,EPS,
     +                 ROOT,S,SBB,SI
      INTEGER          I,II,N1,N2,N2M1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC        ABS,SQRT
C     ..
C     .. Data statements ..
      DATA             A34/0.0/,EPS/0.6E-7/
C     ..
C
C
      VALUE(1) = A(1)
      IF (M.NE.1) THEN
          DO 10 I = 2,M
              VALUE(I) = A(I)
              OFF(I) = B(I)
   10         CONTINUE
C
C---- Each qr iteration is performed of rows and columns n1 to n2
C
          N2 = M
   20     CONTINUE
          IF (N2.GT.1) THEN
              DO 30 II = 2,N2
                  N1 = 2 + N2 - II
                  IF (ABS(OFF(N1)).LE. (ABS(VALUE(N1-1))+
     +                ABS(VALUE(N1)))*EPS) GO TO 40
   30             CONTINUE
              N1 = 1
   40         IF (N2.NE.N1) THEN
C
C---- Root  is the eigenvalue of the bottom 2*2 matrix that is nearest
C     to the last matrix element and is used to accelerate the
C     convergence
C
                  BB = (VALUE(N2)-VALUE(N2-1))*0.5
                  CC = OFF(N2)*OFF(N2)
                  SBB = 1.0
                  IF (BB.LT.0.0) SBB = -1.0
                  ROOT = CC/ (SQRT(BB*BB+CC)*SBB+BB) + VALUE(N2)
                  N2M1 = N2 - 1
                  A22 = VALUE(N1)
                  A12 = A22 - ROOT
                  A23 = OFF(N1+1)
                  A13 = A23
                  DO 50 I = N1,N2M1
                      A33 = VALUE(I+1)
                      IF (I.NE.N2M1) A34 = OFF(I+2)
                      S = SQRT(A12*A12+A13*A13)
                      SI = A13/S
                      CO = A12/S
                      IF (I.NE.N1) OFF(I) = S
                      A11 = CO*A22 + SI*A23
                      A12 = CO*A23 + SI*A33
                      A13 = SI*A34
                      A21 = CO*A23 - SI*A22
                      A22 = CO*A33 - SI*A23
                      A23 = CO*A34
                      VALUE(I) = A11*CO + A12*SI
                      A12 = -A11*SI + A12*CO
                      OFF(I+1) = A12
                      A22 = A22*CO - A21*SI
   50                 CONTINUE
                  VALUE(N2) = A22

              ELSE
                  N2 = N2 - 1
              END IF

              GO TO 20

          END IF

      END IF

      END
C
C
C---- Changes put in to make it work on vax (this version has scaling)
C  1. 12 statements for calculation of over/underflow of determinant in ma21
C     replaced by a simple one which will overflow more easily(entry ma21cd)
C  2. changes to mc10ad to replace 370-specific parts....
C     a. u=floati  (6 times)
C     b. new alog16 procedure (twice)
C     c. simpler 16**diag statement (once)
C  3. all double precision replaced by real*8 (not really necessary).
C  4. replace a(n),b(n) by a(1),b(1) in fm02ad to avoid vax array checking.
C
      FUNCTION FA01AS(I)
C     ==================
      DOUBLE PRECISION  G, FA01AS
      INTEGER I
      COMMON/FA01ES/G
      G= 1431655765.D0
C
C
      G=DMOD(G* 9228907.D0,4294967296.D0)
      IF(I.GE.0)FA01AS=G/4294967296.D0
      IF(I.LT.0)FA01AS=2.D0*G/4294967296.D0-1.D0
      RETURN
      END
C
      SUBROUTINE FA01BS(MAX,NRAND)
C     ============================
C
      INTEGER NRAND, MAX
      DOUBLE PRECISION FA01AS
      EXTERNAL FA01AS
      NRAND=INT(FA01AS(1)*FLOAT(MAX))+1
      RETURN
      END
C
      SUBROUTINE FA01CS(IL,IR)
C     ========================
C
      DOUBLE PRECISION  G
      INTEGER IL,IR
      COMMON/FA01ES/G
      G= 1431655765.D0
C
      IL=G/65536.D0
      IR=G-65536.D0*FLOAT(IL)
      RETURN
      END
C
      SUBROUTINE FA01DS(IL,IR)
C     ========================
C
      DOUBLE PRECISION  G
      INTEGER IL,IR
      COMMON/FA01ES/G
      G= 1431655765.D0
C
      G=65536.D0*FLOAT(IL)+FLOAT(IR)
      RETURN
      END
C
C_BEGIN_FM02AD
C
      DOUBLE PRECISION FUNCTION FM02AD(N,A,IA,B,IB)
C     =============================================
C
C     Compute the inner product of two double precision real
C     vectors accumulating the result double precision, when the
C     elements of each vector are stored at some fixed displacement
C     from neighbouring elements. Given vectors A={a(j)},
C     B={b(j)} of length N, evaluates w=a(j)b(j) summed over
C     j=1..N. Can be used to evaluate inner products involving
C     rows of multi-dimensional arrays.
C     It can be used as an alternative to the assembler version,
C     but note that it is likely to be significantly slower in execution.
C
C     Argument list
C     -------------
C
C       N  (I) (integer) The length of the vectors (if N <= 0 FM02AD = 0)
C
C       A  (I) (double precision) The first vector
C
C       IA (I) (integer) Subscript displacement between elements of A
C
C       B  (I) (double precision) The second vector
C
C       IB (I) (integer) Subscript displacement between elements of B
C
C       FM02AD  the result
C
C
C_END_FM02AD
C
      DOUBLE PRECISION  R1,A,B
C
C---- The following statement changed from a(n),b(n) to avoid vax dynamic
C     array check failure.
C
      DIMENSION A(1),B(1)
      INTEGER N,JA,IA,JB,IB,I
C
      R1=0D0
      IF(N.LE.0) GO TO 2
      JA=1
      IF(IA.LT.0) JA=1-(N-1)*IA
      JB=1
      IF(IB.LT.0) JB=1-(N-1)*IB
      I=0
    1 I=I+1
      R1=R1+A(JA)*B(JB)
      JA=JA+IA
      JB=JB+IB
      IF(I.LT.N) GO TO 1
    2 FM02AD=R1
      RETURN
C
      END
C
C
C_BEGIN_ICROSS
C
      SUBROUTINE ICROSS(A,B,C)
C     ========================
C
C    Cross product (integer version)
C
C    A = B x C
C
C     .. Array Arguments ..
      INTEGER           A(3),B(3),C(3)
C
C_END_ICROSS
C     ..
      A(1) = B(2)*C(3) - C(2)*B(3)
      A(2) = B(3)*C(1) - C(3)*B(1)
      A(3) = B(1)*C(2) - C(1)*B(2)
      END
C
C
C_BEGIN_IDOT
C
      INTEGER FUNCTION IDOT(A,B)
C     ==========================
C
C      Dot product (integer version)
C
C      IDOT = A . B
C
C     .. Array Arguments ..
      INTEGER               A(3),B(3)
C
C_END_IDOT
C     ..
      IDOT = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
      END
C
C
C_BEGIN_IMINV3
C
      SUBROUTINE IMINV3(A,B,D)
C     =======================
C
C     Invert a general 3x3 matrix and return determinant in D
C     (integer version)
C
C     A = (B)-1
C
C     .. Scalar Arguments ..
      INTEGER          D
C     ..
C     .. Array Arguments ..
      INTEGER          A(3,3),B(3,3)
C
C_END_IMINV3
C     ..
C     .. Local Scalars ..
      INTEGER          I,J
C     ..
C     .. Local Arrays ..
      INTEGER          C(3,3)
C     ..
C     .. External Functions ..
      INTEGER          IDOT
      EXTERNAL         IDOT
C     ..
C     .. External Subroutines ..
      EXTERNAL         ICROSS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC        ABS
C     ..
      CALL ICROSS(C(1,1),B(1,2),B(1,3))
      CALL ICROSS(C(1,2),B(1,3),B(1,1))
      CALL ICROSS(C(1,3),B(1,1),B(1,2))
      D = IDOT(B(1,1),C(1,1))
C
C---- Test determinant
C
      IF (ABS(D).GT.0) THEN
C
C---- Determinant is non-zero
C
          DO 20 I = 1,3
              DO 10 J = 1,3
                  A(I,J) = C(J,I)/D
   10             CONTINUE
   20         CONTINUE

      ELSE
          D = 0
      END IF

      END
C
C
C_BEGIN_MATMUL
C
      SUBROUTINE MATMUL(A,B,C)
C     ========================
C
C      Multiply two 3x3 matrices
C
C      A = BC
C
C     .. Array Arguments ..
      REAL              A(3,3),B(3,3),C(3,3)
C
C_END_MATMUL
C     ..
C     .. Local Scalars ..
      INTEGER           I,J,K
C     ..
      DO 30 I = 1,3
          DO 20 J = 1,3
            A(I,J) = 0.0
              DO 10 K = 1,3
                  A(I,J) = B(I,K)*C(K,J) + A(I,J)
   10             CONTINUE
   20         CONTINUE
   30     CONTINUE
      END
C
C
C_BEGIN_MATMULNM
C
      SUBROUTINE MATMULNM(N,M,A,B,C)
C     ========================
C
C      Multiply  NxM  MXN matrices
C
C      A = BC
C
C     .. Array Arguments ..
      INTEGER N,M
      REAL              A(N,N),B(N,M),C(M,N)
C
C_END_MATMUL
C     ..
C     .. Local Scalars ..
      INTEGER           I,J,K
C     ..
      DO 30 I = 1,N
          DO 20 J = 1,N
            A(I,J) = 0.0
              DO 10 K = 1,M
                  A(I,J) = B(I,K)*C(K,J) + A(I,J)
   10             CONTINUE
   20         CONTINUE
   30     CONTINUE
      END
C
C
C_BEGIN_MATVEC
C
      SUBROUTINE MATVEC(V,A,B)
C     ========================
C
C      Post-multiply a 3x3 matrix by a vector
C
C      V = AB
C
C     .. Array Arguments ..
      REAL              A(3,3),B(3),V(3)
C
C_END_MATVEC
C     ..
C     .. Local Scalars ..
      REAL              S
      INTEGER           I,J
C     ..
      DO 20 I = 1,3
          S = 0
          DO 10 J = 1,3
              S = A(I,J)*B(J) + S
   10         CONTINUE
          V(I) = S
   20     CONTINUE
      END
C
C
C_BEGIN_MATMULGEN
C
      SUBROUTINE MATMULGEN(Nb,Mbc,Nc,A,B,C)
C     =====================================
C
C      Generalised matrix multiplication subroutine
C      Multiplies a NbxMbc matrix (B) by a  MbcXNc
C      (C) matrix, so that
C
C      A = BC
C
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER           Nb,Mbc,Nc
C     ..
C     .. Array Arguments ..
      REAL              A(Nb,Nc),B(Nb,Mbc),C(Mbc,Nc)
C
C_END_MATMULGEN
C     ..
C     .. Local Scalars ..
      INTEGER           I,J,K
C     ..
      DO 30 J = 1,Nc
          DO 20 I = 1,Nb
            A(I,J) = 0.0
              DO 10 K = 1,Mbc
                  A(I,J) = B(I,K)*C(K,J) + A(I,J)
   10         CONTINUE
   20     CONTINUE
   30   CONTINUE
      END
C
C
C_BEGIN_MC04B
C
      SUBROUTINE MC04B(A,ALPHA,BETA,M,IA,Q)
C     =====================================
C
C     Transforms a real symmetric matrix A={a(i,j)}, i, j=1..IA
C     into a tri-diagonal matrix having the same eigenvalues as A
C     using Householder's method.
C
C     .. Scalar Arguments ..
      INTEGER          IA,M
C     ..
C     .. Array Arguments ..
      REAL             A(IA,1),ALPHA(1),BETA(1),Q(1)
C
C_END_MC04B
C     ..
C     .. Local Scalars ..
      REAL             BIGK,H,PP,PP1,QJ
      INTEGER          I,I1,I2,J,J1,KI,KJ,M1,M2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC        SQRT
C     ..
      ALPHA(1) = A(1,1)
      DO 20 J = 2,M
          J1 = J - 1
          DO 10 I = 1,J1
              A(I,J) = A(J,I)
   10         CONTINUE
          ALPHA(J) = A(J,J)
   20     CONTINUE
      M1 = M - 1
      M2 = M - 2
      DO 110 I = 1,M2
          PP = 0.0
          I1 = I + 1
          DO 30 J = I1,M
              PP = A(I,J)**2 + PP
   30         CONTINUE
          PP1 = SQRT(PP)
          IF (A(I,I+1).LT.0) THEN
              BETA(I+1) = PP1

          ELSE
              BETA(I+1) = -PP1
          END IF

          IF (PP.GT.0) THEN
              H = PP - BETA(I+1)*A(I,I+1)
              A(I,I+1) = A(I,I+1) - BETA(I+1)
              DO 60 KI = I1,M
                  QJ = 0.0
                  DO 40 KJ = I1,KI
                      QJ = A(KJ,KI)*A(I,KJ) + QJ
   40                 CONTINUE
                  IF (KI-M.LT.0) THEN
                      I2 = KI + 1
                      DO 50 KJ = I2,M
                          QJ = A(KI,KJ)*A(I,KJ) + QJ
   50                     CONTINUE
                  END IF

                  Q(KI) = QJ/H
   60             CONTINUE
              BIGK = 0.0
              DO 70 KJ = I1,M
                  BIGK = A(I,KJ)*Q(KJ) + BIGK
   70             CONTINUE
              BIGK = BIGK/ (2.0*H)
              DO 80 KJ = I1,M
                  Q(KJ) = Q(KJ) - A(I,KJ)*BIGK
   80             CONTINUE
              DO 100 KI = I1,M
                  DO 90 KJ = KI,M
                      A(KI,KJ) = A(KI,KJ) - Q(KI)*A(I,KJ) -
     +                           Q(KJ)*A(I,KI)
   90                 CONTINUE
  100             CONTINUE
          END IF

  110     CONTINUE
      DO 120 I = 2,M
          H = ALPHA(I)
          ALPHA(I) = A(I,I)
          A(I,I) = H
  120     CONTINUE
      BETA(M) = A(M-1,M)
      END
C
C
C_BEGIN_MINVN
C
      SUBROUTINE MINVN(A,N,D,L,M)
C     ===========================
C
C
C---- Purpose
C     =======
C
C           invert a matrix
C
C---- Usage
C     ======
C
C           CALL MINVN(A,N,D,L,M)
C
C---- Description of parameters
C     =========================
C
C    A - input matrix, destroyed in computation and replaced by
C        resultant inverse.
C
C    N - order of matrix A
C
C    D - resultant determinant
C
C    L - work vector of length n
C
C    M - work vector of length n
C
C---- Remarks
C     =======
C
C     Matrix a must be a general matrix
C
C---- Subroutines and function subprograms required
C     =============================================
C
C           NONE
C
C---- Method
C     ======
C
C     The standard gauss-jordan method is used. the determinant
C     is also calculated. a determinant of zero indicates that
C     the matrix is singular.
C
C
C---- Note
C     =====
C
C     If a double precision version of this routine is desired, the
C     c in column 1 should be removed from the double precision
C     statement which follows.
C
C     double precision a,d,biga,hold
C
C        the c must also be removed from double precision statements
C        appearing in other routines used in conjunction with this
C        routine.
C
C        The double precision version of this subroutine must also
C        contain double precision fortran functions.  abs in statement
C        10 must be changed to dabs.
C
ccc     REAL*8 D
C
C_END_MINVN
C
C---- Search for largest element
C
C     .. Scalar Arguments ..
      REAL D
      INTEGER N
C     ..
C     .. Array Arguments ..
      REAL A(N*N)
      INTEGER L(N),M(N)
C     ..
C     .. Local Scalars ..
      REAL BIGA,HOLD
      INTEGER I,IJ,IK,IZ,J,JI,JK,JP,JQ,JR,K,KI,KJ,KK,NK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C
C
      D = 1.0
      NK = -N
      DO 90 K = 1,N
        NK = NK + N
        L(K) = K
        M(K) = K
        KK = NK + K
        BIGA = A(KK)
        DO 20 J = K,N
          IZ = (J-1)*N
          DO 10 I = K,N
            IJ = IZ + I
            IF ((ABS(BIGA)-ABS(A(IJ))).LT.0.0) THEN
              BIGA = A(IJ)
              L(K) = I
              M(K) = J
            END IF
   10     CONTINUE
   20   CONTINUE
C
C---- Interchange rows
C
        J = L(K)
        IF ((J-K).GT.0) THEN
          KI = K - N
          DO 30 I = 1,N
            KI = KI + N
            HOLD = -A(KI)
            JI = KI - K + J
            A(KI) = A(JI)
            A(JI) = HOLD
   30     CONTINUE
        END IF
C
C---- Interchange columns
C
        I = M(K)
        IF ((I-K).GT.0) THEN
          JP = (I-1)*N
          DO 40 J = 1,N
            JK = NK + J
            JI = JP + J
            HOLD = -A(JK)
            A(JK) = A(JI)
            A(JI) = HOLD
   40     CONTINUE
        END IF
C
C---- Divide column by minus pivot (value of pivot element is
C     contained in biga)
C
        IF (BIGA.NE.0.0) THEN
          DO 50 I = 1,N
            IF ((I-K).NE.0) THEN
              IK = NK + I
              A(IK) = A(IK)/ (-BIGA)
            END IF
   50     CONTINUE
C
C---- Reduce matrix
C
          DO 70 I = 1,N
            IK = NK + I
            HOLD = A(IK)
            IJ = I - N
            DO 60 J = 1,N
              IJ = IJ + N
              IF ((I-K).NE.0) THEN
                IF ((J-K).NE.0) THEN
                  KJ = IJ - I + K
                  A(IJ) = A(KJ)*HOLD + A(IJ)
                END IF
              END IF
   60       CONTINUE
   70     CONTINUE
C
C---- Divide row by pivot
C
          KJ = K - N
          DO 80 J = 1,N
            KJ = KJ + N
            IF ((J-K).NE.0) A(KJ) = A(KJ)/BIGA
   80     CONTINUE
C
C---- Product of pivots
C
          D = D*BIGA
C
C---- Replace pivot by reciprocal
C
          A(KK) = 1.0/BIGA
        ELSE
          GO TO 130
        END IF
   90 CONTINUE
C
C---- Final row and column interchange
C
      K = N
  100 CONTINUE
      K = (K-1)
      IF (K.GT.0) THEN
        I = L(K)
        IF ((I-K).GT.0) THEN
          JQ = (K-1)*N
          JR = (I-1)*N
          DO 110 J = 1,N
            JK = JQ + J
            HOLD = A(JK)
            JI = JR + J
            A(JK) = -A(JI)
            A(JI) = HOLD
  110     CONTINUE
        END IF
        J = M(K)
        IF ((J-K).GT.0) THEN
          KI = K - N
          DO 120 I = 1,N
            KI = KI + N
            HOLD = A(KI)
            JI = KI - K + J
            A(KI) = -A(JI)
            A(JI) = HOLD
  120     CONTINUE
        END IF
        GO TO 100
      ELSE
        RETURN
      END IF
  130 D = 0.0
C
C
      END
C
C
C_BEGIN_MINV3
C
      SUBROUTINE MINV3(A,B,D)
C     ======================
C
C     Invert a general 3x3 matrix and return determinant in D
C
C     A = (B)-1
C
C     .. Scalar Arguments ..
      REAL            D
C     ..
C     .. Array Arguments ..
      REAL            A(3,3),B(3,3)
C
C_END_MINV3
C     ..
C     .. Local Scalars ..
      INTEGER         I,J
C     ..
C     .. Local Arrays ..
      REAL            C(3,3)
C     ..
C     .. External Functions ..
      REAL            DOT
      EXTERNAL        DOT
C     ..
C     .. External Subroutines ..
      EXTERNAL        CROSS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC       ABS
C     ..
      CALL CROSS(C(1,1),B(1,2),B(1,3))
      CALL CROSS(C(1,2),B(1,3),B(1,1))
      CALL CROSS(C(1,3),B(1,1),B(1,2))
      D = DOT(B(1,1),C(1,1))
C
C---- Test determinant
C
      IF (ABS(D).GT.1.0E-30) THEN
C
C---- Determinant is non-zero
C
          DO 20 I = 1,3
              DO 10 J = 1,3
                  A(I,J) = C(J,I)/D
   10             CONTINUE
   20         CONTINUE

      ELSE
          D = 0.0
      END IF

      END
C
C
C_BEGIN_RANMAR
C
      SUBROUTINE RANMAR(RVEC,LEN)
C     ===========================
C
C     Universal random number generator proposed by Marsaglia and Zaman
C     in report FSU-SCRI-87-50
C     slightly modified by F. James, 1988 to generate a vector
C     of pseudorandom numbers RVEC of length LEN
C     and making the COMMON block include everything needed to
C     specify completely the state of the generator.
C     Transcribed from CERN report DD/88/22.
C     Rather inelegant messing about added by D. Love, Jan. 1989 to
C     make sure initialisation always occurs.
C     *** James says that this is the preferred generator.
C     Gives bit-identical results on all machines with at least
C     24-bit mantissas in the flotaing point representation (i.e.
C     all common 32-bit computers. Fairly fast, satisfies very
C     stringent tests, has very long period and makes it very
C     simple to generate independly disjoint sequences.
C     See also RANECU.
C     The state of the generator may be saved/restored using the
C     whole contents of /RASET1/.
C     Call RANMAR to get a vector, RMARIN to initialise. 
C     
C  Argument list
C  -------------
C
C     VREC (O)                 (REAL)   Random Vector
C
C     LEN  (I)              (INTEGER)   Length of random vector
C
C
C  For ENTRY point RMARIN
C  ----------------------
C
C     Initialisation for RANMAR.  The input values should
C     be in the ranges: 0<=ij<=31328, 0<=kl<=30081
C     This shows the correspondence between the simplified input seeds
C     IJ, KL and the original Marsaglia-Zaman seeds i,j,k,l
C     To get standard values in Marsaglia-Zaman paper,
C     (I=12, J=34, K=56, L=78) put IJ=1802, KL=9373
C
C     IJ   (I)              (INTEGER)   Seed for random number generator
C
C     KL   (I)              (INTEGER)   Seed for randon number generator
C
C_END_RANMAR
C     
C     ..
C     .. Agruments ..
      REAL RVEC(*)
      INTEGER LEN,IJ,KL
C     ..
C     .. Common Variables ..
      REAL C,CD,CM,U
      INTEGER I97,J97
C     ..
C     .. Local Scalars ..
      REAL S,T,UNI
      INTEGER I,II,IVEC,J,JJ,K,L,M
      LOGICAL INITED
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
C     .. Common Blocks ..
      COMMON /RASET1/ U(97),C,CD,CM,I97,J97
C     ..
C     .. Save Statement ..
      SAVE INITED, /RASET1/
C     ..
C     .. Data Statement ..
      DATA INITED /.FALSE./
C
C---- If initialised, fill RVEC and RETURN. If not, do initialisation
C     and return here later.
C
 1    IF (INITED) THEN
        DO 100 IVEC=1,LEN
          UNI=U(I97)-U(J97)
          IF (UNI.LT.0.) UNI=UNI+1.
          U(I97)=UNI
          I97=I97-1
          IF (I97.EQ.0) I97=97
          J97=J97-1
          IF (J97.EQ.0) J97=97
          C=C-CD
          IF (C.LT.0.) C=C+CM
          UNI=UNI-C
          IF (UNI.LT.0.) UNI=UNI+1.
          RVEC(IVEC)=UNI
 100    CONTINUE
        RETURN
      ENDIF
      I=MOD(1802/177,177)+2
      J=MOD(1802,177)+2
      K=MOD(9373/169,178)+1
      L=MOD(9373,169)
C     
      GOTO 10
C
C---- Initialise and return without filling RVEC
C
      ENTRY RMARIN(IJ,KL)
      I=MOD(IJ/177,177)+2
      J=MOD(IJ,177)+2
      K=MOD(KL/169,178)+1
      L=MOD(KL,169)
      INITED=.TRUE.

 10   CONTINUE
      DO 2 II=1,97
        S=0.
        T=.5
        DO 3 JJ=1,24
          M=MOD(MOD(I*J,179)*K,179)
          I=J
          J=K
          K=M
          L=MOD(53*L+1,169)
          IF (MOD(L*M,64).GE.32) S=S+T
          T=0.5*T
 3      CONTINUE
        U(II)=S
 2    CONTINUE
      C=362436./16777216.
      CD=7654321./16777216.
      CM=16777213./16777216.
      I97=97
      J97=33
      IF (.NOT. INITED) THEN
        INITED=.TRUE.
        GOTO 1
      ENDIF

      END
C
C
C_BEGIN_SCALEV
C
      SUBROUTINE SCALEV(A,X,B)
C     ========================
C
C     Scale vector B with scalar X and put result in A
C
C     .. Scalar Arguments ..
      REAL              X
C     ..
C     .. Array Arguments ..
      REAL              A(3),B(3)
C
C_END_SCALEV
C     ..
C     .. Local Scalars ..
      INTEGER           I
C     ..
      DO 10 I = 1,3
          A(I) = B(I)*X
   10     CONTINUE
      END
C
C
C_BEGIN_TRANSP
C
      SUBROUTINE TRANSP(A,B)
C     ======================
C
C---- Transpose a 3x3 matrix
C
C     A = BT
C
C     .. Array Arguments ..
      REAL              A(3,3),B(3,3)
C
C_END_TRANSP
C     ..
C     .. Local Scalars ..
      INTEGER           I,J
C     ..
      DO 20 I = 1,3
          DO 10 J = 1,3
              A(I,J) = B(J,I)
   10         CONTINUE
   20     CONTINUE
      END
C
C
C_BEGIN_UNIT
C
      SUBROUTINE UNIT(V)
C     =================
C
C    Vector V reduced to unit vector
C
C     .. Array Arguments ..
      REAL            V(3)
C
C_END_UNIT
C     ..
C     .. Local Scalars ..
      REAL            VMOD
      INTEGER         I
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC       SQRT
C     ..
      VMOD = V(1)**2 + V(2)**2 + V(3)**2
      VMOD = SQRT(VMOD)
      DO 10 I = 1,3
          V(I) = V(I)/VMOD
   10     CONTINUE
      END
C
C
C_BEGIN_VDIF
C
      SUBROUTINE VDIF(A,B,C)
C     =====================
C
C      Subtract two vectors
C
C      A = B - C
C
C     .. Array Arguments ..
      REAL            A(3),B(3),C(3)
C
C_END_VDIF
C     ..
C     .. Local Scalars ..
      INTEGER         I
C     ..
      DO 10 I = 1,3
          A(I) = B(I) - C(I)
   10     CONTINUE
      END
C
C
C_BEGIN_VSET
C
      SUBROUTINE VSET(A,B)
C     ====================
C
C      Copy a vector from B to A
C
C     .. Array Arguments ..
      REAL            A(3),B(3)
C
C_END_VSET
C     ..
C     .. Local Scalars ..
      INTEGER         I
C     ..
      DO 10 I = 1,3
          A(I) = B(I)
   10     CONTINUE
      END
C
C
C_BEGIN_VSUM
C
      SUBROUTINE VSUM(A,B,C)
C     ======================
C
C     Add two vectors
C
C     A = B + C
C
C     .. Array Arguments ..
      REAL            A(3),B(3),C(3)
C
C_END_VSUM
C     ..
C     .. Local Scalars ..
      INTEGER         I
C     ..
      DO 10 I = 1,3
          A(I) = B(I) + C(I)
   10     CONTINUE
      END
C
C
C_BEGIN_ZIPIN
C
      SUBROUTINE ZIPIN(ID,N,BUF)
C     ==========================
C
C     Fast binary read on unit ID into real array BUF of length N
C
C     .. Scalar Arguments ..
      INTEGER          ID,N
C     ..
C     .. Array Arguments ..
      REAL             BUF(N)
C
C_END_ZIPIN
C     ..
      READ (ID) BUF
      END
C
C
C_BEGIN_ZIPOUT
C
      SUBROUTINE ZIPOUT(ID,N,BUF)
C     ===========================
C
C     Fast binary write to unit ID of real array BUF length N
C
C     .. Scalar Arguments ..
      INTEGER           ID,N
C     ..
C     .. Array Arguments ..
      REAL              BUF(N)
C
C_END_ZIPOUT
C     ..
      WRITE (ID) BUF
      END
C
C The following routines calculate Bessel functions of
C the first kind, their derivatives, and their zeros.
C Subroutine ZJVX is from Ian Tickle.
C Subroutines JDVX, JVX, GAMMA and functions MSTA1, MSTA2
C are from MJYV at http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
C These routines are copyrighted, but permission is given to incorporate these
C routines into programs provided that the copyright is acknowledged. 

      SUBROUTINE ZJVX(V,IZ,BJ,DJ,X)
C
C     =======================================================
C     Purpose: Compute zero of Bessel function Jv(x)
C     Input :  V     --- Order of Jv(x)
C              IZ    --- Index of previous zero of Jv(x)
C              X     --- Value of previous zero of Jv(x)
C     Output:  BJ(N) --- Jv(x) for N = 0...INT(V)
C              DJ(N) --- J'v(x) for N = 0...INT(V)
C              X     --- Value of next zero of Jv(x)
C
      IMPLICIT NONE
      INTEGER IZ,N
      REAL*8 V,VM,X,X0
      REAL*8 BJ(0:*),DJ(0:*)
C
C#### Initial guess at zero: needs previous value if not first.
      N=V
      IF (IZ.EQ.0) THEN
        X=1.99535+.8333883*SQRT(V)+.984584*V
      ELSEIF (N.LE.10) THEN
        X=X+3.11+.0138*V+(.04832+.2804*V)/(IZ+1)**2
      ELSE
        X=X+3.001+.0105*V+(11.52+.48525*V)/(IZ+3)**2
      ENDIF
C
C#### Polish zero by Newton-Cotes iteration.
1     CALL JDVX(V,X,VM,BJ,DJ)
      IF (INT(VM).NE.N) CALL CCPERR(1,'VM != N in ZJVX.')
      X0=X
      X=X-BJ(N)/DJ(N)
      IF (ABS(X-X0).GT.1D-10) GOTO 1
      END
C
C
      SUBROUTINE JDVX(V,X,VM,BJ,DJ)
C
C     =======================================================
C     Purpose: Compute Bessel functions Jv(x) and their derivatives
C     Input :  x --- Argument of Jv(x)
C              v --- Order of Jv(x)
C                    ( v = n+v0, 0 <= v0 < 1, n = 0,1,2,... )
C     Output:  BJ(n) --- Jn+v0(x)
C              DJ(n) --- Jn+v0'(x)
C              VM --- Highest order computed
C     Routines called:
C          (1) GAMMA for computing gamma function
C          (2) MSTA1 and MSTA2 for computing the starting
C              point for backward recurrence
C     =======================================================
C
      IMPLICIT NONE
      INTEGER J,K,K0,L,M,N
      REAL*8 A,A0,BJV0,BJV1,BJVL,CK,CS,F,F0,F1,F2,GA,PI,PX,QX,R,RP,RP2,
     &RQ,SK,V,V0,VG,VL,VM,VV,X,X2,XK
      PARAMETER(PI=3.141592653589793D0, RP2=.63661977236758D0)
      REAL*8 BJ(0:*),DJ(0:*)
      INTEGER MSTA1,MSTA2
C
      X2=X*X
      N=INT(V)
      V0=V-N
      IF (X.LT.1D-100) THEN
        DO K=0,N
          BJ(K)=0D0
          DJ(K)=0D0
        ENDDO
        IF (V0.EQ.0D0) THEN
          BJ(0)=1D0
          DJ(1)=.5D0
        ELSE
          DJ(0)=1D300
        ENDIF
        VM=V
      ELSE
        IF (X.LE.12D0) THEN
          DO L=0,1
            VL=V0+L
            BJVL=1D0
            R=1D0
            DO K=1,40
              R=-.25D0*R*X2/(K*(K+VL))
              BJVL=BJVL+R
              IF (ABS(R).LT.ABS(BJVL)*1D-15) GOTO 20
            ENDDO
20          VG=1D0+VL
            CALL GAMMA(VG,GA)
            A=(.5D0*X)**VL/GA
            IF (L.EQ.0) THEN
              BJV0=BJVL*A
            ELSEIF (L.EQ.1) THEN
              BJV1=BJVL*A
            ENDIF
          ENDDO
        ELSE
          K0=11
          IF (X.GE.35D0) K0=10
          IF (X.GE.50D0) K0=8
          DO J=0,1
            VV=4D0*(J+V0)*(J+V0)
            PX=1D0
            RP=1D0
            DO K=1,K0
              RP=-.78125D-2*RP*(VV-(4*K-3)**2)*(VV-(4*K-1)**2)/
     &        (K*(2*K-1)*X2)
              PX=PX+RP
            ENDDO
            QX=1D0
            RQ=1D0
            DO K=1,K0
              RQ=-.78125D-2*RQ*(VV-(4*K-1)**2)*(VV-(4*K+1)**2)/
     &        (K*(2*K+1)*X2)
              QX=QX+RQ
            ENDDO
            QX=.125D0*(VV-1D0)*QX/X
            XK=X-(.5D0*(J+V0)+.25D0)*PI
            A0=SQRT(RP2/X)
            CK=COS(XK)
            SK=SIN(XK)
            IF (J.EQ.0) THEN
              BJV0=A0*(PX*CK-QX*SK)
            ELSEIF (J.EQ.1) THEN
              BJV1=A0*(PX*CK-QX*SK)
            ENDIF
          ENDDO
        ENDIF
        BJ(0)=BJV0
        BJ(1)=BJV1
        DJ(0)=V0/X*BJ(0)-BJ(1)
        DJ(1)=BJ(0)-(1D0+V0)/X*BJ(1)
        IF (N.GE.2 .AND. N.LE.INT(.9*X)) THEN
          F0=BJV0
          F1=BJV1
          DO K=2,N
            F=2D0*(K+V0-1D0)/X*F1-F0
            BJ(K)=F
            F0=F1
            F1=F
          ENDDO
        ELSEIF (N.GE.2) THEN
          M=MSTA1(X,200)
          IF (M.LT.N) THEN
            N=M
          ELSE
            M=MSTA2(X,N,15)
          ENDIF
          F2=0D0
          F1=1D-100
          DO K=M,0,-1
            F=2D0*(V0+K+1D0)/X*F1-F2
            IF (K.LE.N) BJ(K)=F
            F2=F1
            F1=F
          ENDDO
          IF (ABS(BJV0).GT.ABS(BJV1)) THEN
            CS=BJV0/F
          ELSE
            CS=BJV1/F2
          ENDIF
          DO K=0,N
            BJ(K)=CS*BJ(K)
          ENDDO
        ENDIF
        DO K=2,N
          DJ(K)=BJ(K-1)-(K+V0)/X*BJ(K)
        ENDDO
        VM=N+V0
      ENDIF
      END
C
C
      SUBROUTINE JVX(V,X,VM,BJ)
C
C     =======================================================
C     Purpose: Compute Bessel functions Jv(x)
C     Input :  x --- Argument of Jv(x)
C              v --- Order of Jv(x)
C                    ( v = n+v0, 0 <= v0 < 1, n = 0,1,2,... )
C     Output:  BJ(n) --- Jn+v0(x)
C              VM --- Highest order computed
C     Routines called:
C          (1) GAMMA for computing gamma function
C          (2) MSTA1 and MSTA2 for computing the starting
C              point for backward recurrence
C     =======================================================
C
      IMPLICIT NONE
      INTEGER J,K,K0,L,M,N
      REAL*8 A,A0,BJV0,BJV1,BJVL,CK,CS,F,F0,F1,F2,GA,PI,PX,QX,R,RP,RP2,
     &RQ,SK,V,V0,VG,VL,VM,VV,X,X2,XK
      PARAMETER(PI=3.141592653589793D0, RP2=.63661977236758D0)
      REAL*8 BJ(0:*)
      INTEGER MSTA1,MSTA2
C
      X2=X*X
      N=INT(V)
      V0=V-N
      IF (X.LT.1D-100) THEN
        DO K=0,N
          BJ(K)=0D0
        ENDDO
        IF (V0.EQ.0D0) BJ(0)=1D0
        VM=V
      ELSE
        IF (X.LE.12D0) THEN
          DO L=0,1
            VL=V0+L
            BJVL=1D0
            R=1D0
            DO K=1,40
              R=-.25D0*R*X2/(K*(K+VL))
              BJVL=BJVL+R
              IF (ABS(R).LT.ABS(BJVL)*1D-15) GOTO 20
            ENDDO
20          VG=1D0+VL
            CALL GAMMA(VG,GA)
            A=(.5D0*X)**VL/GA
            IF (L.EQ.0) THEN
              BJV0=BJVL*A
            ELSEIF (L.EQ.1) THEN
              BJV1=BJVL*A
            ENDIF
          ENDDO
        ELSE
          K0=11
          IF (X.GE.35D0) K0=10
          IF (X.GE.50D0) K0=8
          DO J=0,1
            VV=4D0*(J+V0)*(J+V0)
            PX=1D0
            RP=1D0
            DO K=1,K0
              RP=-.78125D-2*RP*(VV-(4*K-3)**2)*(VV-(4*K-1)**2)/
     &        (K*(2*K-1)*X2)
              PX=PX+RP
            ENDDO
            QX=1D0
            RQ=1D0
            DO K=1,K0
              RQ=-.78125D-2*RQ*(VV-(4*K-1)**2)*(VV-(4*K+1)**2)/
     &        (K*(2*K+1)*X2)
              QX=QX+RQ
            ENDDO
            QX=.125D0*(VV-1D0)*QX/X
            XK=X-(.5D0*(J+V0)+.25D0)*PI
            A0=SQRT(RP2/X)
            CK=COS(XK)
            SK=SIN(XK)
            IF (J.EQ.0) THEN
              BJV0=A0*(PX*CK-QX*SK)
            ELSEIF (J.EQ.1) THEN
              BJV1=A0*(PX*CK-QX*SK)
            ENDIF
          ENDDO
        ENDIF
        BJ(0)=BJV0
        BJ(1)=BJV1
        IF (N.GE.2 .AND. N.LE.INT(.9*X)) THEN
          F0=BJV0
          F1=BJV1
          DO K=2,N
            F=2D0*(K+V0-1D0)/X*F1-F0
            BJ(K)=F
            F0=F1
            F1=F
          ENDDO
        ELSEIF (N.GE.2) THEN
          M=MSTA1(X,200)
          IF (M.LT.N) THEN
            N=M
          ELSE
            M=MSTA2(X,N,15)
          ENDIF
          F2=0D0
          F1=1D-100
          DO K=M,0,-1
            F=2D0*(V0+K+1D0)/X*F1-F2
            IF (K.LE.N) BJ(K)=F
            F2=F1
            F1=F
          ENDDO
          IF (ABS(BJV0).GT.ABS(BJV1)) THEN
            CS=BJV0/F
          ELSE
            CS=BJV1/F2
          ENDIF
          DO K=0,N
            BJ(K)=CS*BJ(K)
          ENDDO
        ENDIF
        VM=N+V0
      ENDIF
      END
C
      SUBROUTINE GAMMA(X,GA)
C
C     ==================================================
C     Purpose: Compute gamma function GA(x)
C     Input :  x  --- Argument of GA(x)
C                     ( x is not equal to 0,-1,-2,... )
C     Output:  GA --- GA(x)
C     ==================================================
C
      IMPLICIT NONE
      INTEGER K,M,M1
      REAL*8 GA,GR,PI,R,X,Z
      PARAMETER(PI=3.141592653589793D0)
      REAL*8 G(26)
      DATA G/1D0, .5772156649015329D0,
     &-.6558780715202538D0, -.420026350340952D-1,
     &.1665386113822915D0, -.421977345555443D-1,
     &-.96219715278770D-2, .72189432466630D-2,
     &-.11651675918591D-2, -.2152416741149D-3,
     &.1280502823882D-3, -.201348547807D-4,
     &-.12504934821D-5, .11330272320D-5,
     &-.2056338417D-6, .61160950D-8,
     &.50020075D-8, -.11812746D-8,
     &.1043427D-9, .77823D-11,
     &-.36968D-11, .51D-12,
     &-.206D-13, -.54D-14, .14D-14, .1D-15/
C
      IF (X.EQ.INT(X)) THEN
        IF (X.GT.0D0) THEN
          GA=1D0
          M1=X-1
          DO K=2,M1
            GA=GA*K
          ENDDO
        ELSE
          GA=1D300
        ENDIF
      ELSE
        IF (ABS(X).GT.1D0) THEN
          Z=ABS(X)
          M=INT(Z)
          R=1D0
          DO K=1,M
            R=R*(Z-K)
          ENDDO
          Z=Z-M
        ELSE
          Z=X
        ENDIF
        GR=G(26)
        DO K=25,1,-1
          GR=GR*Z+G(K)
        ENDDO
        GA=1D0/(GR*Z)
        IF (ABS(X).GT.1D0) THEN
          GA=GA*R
          IF (X.LT.0D0) GA=-PI/(X*GA*DSIN(PI*X))
        ENDIF
      ENDIF
      END
C
      INTEGER FUNCTION MSTA1(X,MP)
C
C     ===================================================
C     Purpose: Determine the starting point for backward
C              recurrence such that the magnitude of
C              Jn(x) at that point is about 10^(-MP)
C     Input :  x     --- Argument of Jn(x)
C              MP    --- Value of magnitude
C     Output:  MSTA1 --- Starting point
C     ===================================================
C
      IMPLICIT NONE
      INTEGER IT,MP,N,N0,N1,NN
      REAL*8 A0,F,F0,F1,X
      REAL*8 ENVJ
      ENVJ(N,X)=.5D0*LOG10(6.28D0*N)-N*LOG10(1.36D0*X/N)
C
      A0=ABS(X)
      N0=INT(1.1D0*A0)+1
      F0=ENVJ(N0,A0)-MP
      N1=N0+5
      F1=ENVJ(N1,A0)-MP
      DO IT=1,20
        NN=N1-(N1-N0)/(1D0-F0/F1)
        F=ENVJ(NN,A0)-MP
        IF(ABS(NN-N1).LT.1) GOTO 20
        N0=N1
        F0=F1
        N1=NN
        F1=F
      ENDDO
20    MSTA1=NN
      END
C
      INTEGER FUNCTION MSTA2(X,N,MP)
C
C     ===================================================
C     Purpose: Determine the starting point for backward
C              recurrence such that all Jn(x) has MP
C              significant digits
C     Input :  x  --- Argument of Jn(x)
C              n  --- Order of Jn(x)
C              MP --- Significant digit
C     Output:  MSTA2 --- Starting point
C     ===================================================
C
      IMPLICIT NONE
      INTEGER IT,MP,N,N0,N1,NN
      REAL*8 A0,F,F0,F1,EJN,HMP,OBJ,X
      REAL*8 ENVJ
      ENVJ(N,X)=.5D0*LOG10(6.28D0*N)-N*LOG10(1.36D0*X/N)
C
      A0=ABS(X)
      HMP=.5D0*MP
      EJN=ENVJ(N,A0)
      IF (EJN.LE.HMP) THEN
        OBJ=MP
        N0=INT(1.1D0*A0)
      ELSE
        OBJ=HMP+EJN
        N0=N
      ENDIF
      F0=ENVJ(N0,A0)-OBJ
      N1=N0+5
      F1=ENVJ(N1,A0)-OBJ
      DO IT=1,20
        NN=N1-(N1-N0)/(1D0-F0/F1)
        F=ENVJ(NN,A0)-OBJ
        IF (ABS(NN-N1).LT.1) GOTO 20
        N0=N1
        F0=F1
        N1=NN
        F1=F
      ENDDO
20    MSTA2=NN+10
      END
C
