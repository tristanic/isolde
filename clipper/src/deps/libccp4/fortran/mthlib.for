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
      SUBROUTINE ROTMAT( U, ANG, R)
C     =============================
C
C
      DIMENSION U( 3), R( 3, 3), K( 3, 3), V( 3)
      DATA K / 0, +3, -2,  -3, 0, +1,  +2, -1, 0/
C
C
      DTOR = 3.14159265/180.
      S = SIN( ANG * DTOR)
      C = COS( ANG * DTOR)
      UMAG = SQRT( U(1)*U(1)+U(2)*U(2)+U(3)*U(3))
C
C
      DO 999 I = 1, 3
        V( I) = U( I)/ UMAG   
999   CONTINUE
C
C
      DO 998 I = 1, 3
        DO 997 J = 1, 3
          R( I, J) = V( I) * V( J) * ( 1. - C)
C
C
          IF( I .EQ. J ) THEN
            R( I, J) = R( I, J) + C
          ELSE
            R( I, J) = R( I, J)+
     + ISIGN(1, K( I, J))*S*V( IABS( K( I, J)))
          END IF
C
C
997      CONTINUE
998      CONTINUE
C
C
      RETURN
C
C
      END
C
C
C
      SUBROUTINE MATMULT( A, B, AB)
C     =============================
C
C
      DIMENSION A( 3, 3), B( 3, 3), AB( 3, 3)
C
C
      DO 999 I = 1, 3
        DO 998 J = 1, 3
          SUM = 0.
          DO 997 K = 1, 3
            SUM = SUM + A( I, K) * B( K, J)
997       CONTINUE
          AB( I, J) = SUM
998   CONTINUE
999   CONTINUE
C
C
      RETURN
C
C
      END
C
C
C
      FUNCTION DET( A)
C     ================
C
C
      DIMENSION A( 3, 3)
C
C
      DET = A(1,1)*( A(2,2)*A(3,3) - A(2,3)*A(3,2)) +
     +      A(1,2)*( A(2,3)*A(3,1) - A(2,1)*A(3,3)) +
     +      A(1,3)*( A(2,1)*A(3,2) - A(2,2)*A(3,1))
C
C
      RETURN
C
C
      END
C
C
C
      FUNCTION INCMOD3( I, INC)
C     =========================
C
C
      INCMOD3 = I + INC
      IF ( INCMOD3 .GT. 3) INCMOD3 = INCMOD3 - 3
C
C
      RETURN
C
C
      END
C
C
C
      SUBROUTINE MATINV( A, AI)
C     =========================
C
C
      DIMENSION A( 3, 3), AI( 3, 3)
C
C
      D = DET( A)
C
C
      DO 999 I = 1, 3
        I1 = INCMOD3( I, 1)
        I2 = INCMOD3( I, 2)
        DO 998 J = 1, 3
          J1 = INCMOD3( J, 1)
          J2 = INCMOD3( J, 2)
          AI( J, I) = ( A( I1, J1)*A( I2, J2) - 
     +    A( I1, J2)*A( I2, J1))/D
998    CONTINUE
999    CONTINUE
C
C
      RETURN
C
C
      END
C
C
C
      SUBROUTINE MATCOPY( SOURCE, DEST, IROWS, ICOLS)
C     ===============================================
C
C
      DIMENSION SOURCE( IROWS, ICOLS), DEST( IROWS, ICOLS)
C
C
      DO 999 I = 1, IROWS
        DO 998 J = 1, ICOLS
          DEST( I, J) = SOURCE( I, J)
998      CONTINUE
999      CONTINUE
C
C
      RETURN
C
C
      END
C
C
C
      SUBROUTINE COPY_VECTOR( A, B, NX)
C     =================================
C
C
      DIMENSION A( NX), B( NX)
C
C
      DO 999 I = 1, NX
        B( I) = A( I)
999   CONTINUE
C
C
      RETURN
C
C
      END
C
C
C
      SUBROUTINE SCALE_VECTOR( SCALAR, VECTOR, SCALED_VECTOR, NX)
C     ===========================================================
C
C
      DIMENSION VECTOR( NX), SCALED_VECTOR( NX)
C
C
      DO 999 I = 1, NX
        SCALED_VECTOR( I) = SCALAR * VECTOR( I)
999   CONTINUE
C
C
      RETURN
C
C
      END      
C
C
C
      SUBROUTINE ADDVEC( A, B, AB, NX)
C     ================================
C
C
      DIMENSION A( NX), B( NX), AB( NX)
C
C
      DO 999 I = 1, NX
        AB( I) = A( I) + B( I)
999   CONTINUE
C
C
      RETURN
C
C
      END
C
C
C
      SUBROUTINE VECTOR_DIFF( A, B, A_MINUS_B, NX)
C     ============================================
C
C
      DIMENSION A( NX), B( NX), A_MINUS_B( NX)
C
C
      DO 999 I = 1, NX
        A_MINUS_B( I) = A( I) - B( I)
999   CONTINUE
C
C
      RETURN
C
C
      END
C
C
C
      FUNCTION DOTPRODUCT( A, B, NX)
C     ==============================
C
C
      DIMENSION A( NX), B( NX)
C
C
      SUM = 0.
C
C
      DO 999 I = 1, NX
        SUM = SUM + A( I) * B( I)
999   CONTINUE
C
C
      DOTPRODUCT = SUM
C
C
      RETURN
C
C
      END
C
C
C
      SUBROUTINE CROSSPRODUCT( A, B, AB)
C     ==================================
C
C
      DIMENSION A( 3), B( 3), AB( 3)
C
C
      AB( 1) = A( 2) * B( 3) - A( 3) * B( 2)
      AB( 2) = A( 3) * B( 1) - A( 1) * B( 3)
      AB( 3) = A( 1) * B( 2) - A( 2) * B( 1)
C
C
      RETURN
C
C
      END
C
C
C
      SUBROUTINE EXTRACT_MATRIX( A_SET, NX, NY, NSET, ISET, AMAT)
C     ===========================================================
C
C
      DIMENSION A_SET( NX, NY, NSET), AMAT( NX, NY)
C
C
      DO 999 I = 1, NX
        DO 998 J = 1, NY
          AMAT( I, J) = A_SET( I, J, ISET)
998     CONTINUE
999     CONTINUE
C
C
      RETURN
C
C
      END
C
C
C
      SUBROUTINE EXTRACT_VECTOR( V_SET, NX, NSET, ISET, VEC)
C     ======================================================
C
C
      DIMENSION V_SET( NX, NSET), VEC( NX)
C
C
      DO 999 I = 1, NX
        VEC( I) = V_SET( I, ISET)
999   CONTINUE
C
C
      RETURN
C
C
      END
C
C
C
      SUBROUTINE BUILD_3VECTOR( VEC, X, Y, Z)
C     =======================================
C
C
      DIMENSION VEC( 3)
C
C
      VEC( 1) = X
      VEC( 2) = Y
      VEC( 3) = Z
C
C
      RETURN
C
C
      END
C
C
C      
      SUBROUTINE GEN_TRIGTAB
C     ======================
C
C
      COMMON /TRIGTAB/ AMULTIPLIER, COSTAB( 0:9000), SINTAB( 0:9000)
C
C
      MAXARG = 90. * AMULTIPLIER
C
C
      DO 999 I = 0, MAXARG
        ARG = FLOAT( I) / AMULTIPLIER
        COSTAB( I) = COSD( ARG)
        SINTAB( I) = SIND( ARG)
999     CONTINUE
C
C
      RETURN
C
C
      END
C
C
C
      FUNCTION GET_COSD( ARG )
C     ========================
C
C
      COMMON /TRIGTAB/ AMULTIPLIER, COSTAB( 0:9000), SINTAB( 0:9000)
C
C
      A = ARG
      IF ( A .LT. 0. )  A = -A
10    CONTINUE
C
C
      IF ( A .GE. 360. ) THEN
        A = A - 360.
      GOTO 10
      END IF
C
C
      IF ( A .GT. 180. )  A = 360. - A
C
C
      IF ( A .GT. 90. )  THEN
        A = 180. - A
        ASIGN = -1.
      ELSE
        ASIGN = 1.
      END IF
C
C
      IA = INT( A * AMULTIPLIER + .5 )
      GET_COSD = ASIGN * COSTAB( IA )
C
C
      RETURN
C
C
      END
C
C
C
      FUNCTION GET_SIND( ARG )
C     ========================
C
C
      COMMON /TRIGTAB/ AMULTIPLIER, COSTAB( 0:9000), SINTAB( 0:9000)
C
C
      A = ARG
C
C
      IF ( A .LT. 0. )  THEN
        A = -A
        ASIGN = -1.
      ELSE
        ASIGN = 1.
      END IF
C
10     CONTINUE
C
      IF ( A .GE. 360. ) THEN
        A = A - 360.
        GOTO 10
        END IF
C
C
      IF ( A .GT. 180. )  THEN
        A = 360. - A
        ASIGN = -1. * ASIGN
      END IF
C
C
      IF ( A .GT. 90. )  A = 180. - A
      IA = INT( A * AMULTIPLIER + .5 )
      GET_SIND = ASIGN * SINTAB( IA )
C
C
      RETURN
C
C
      END
C
C
C
      FUNCTION RANDOM_GAUSSIAN(  NRANDOM_0_1, IRANSEED )
C     ==================================================
C
C
      SUM = 0.
C
C
      DO 999 I = 1, NRANDOM_0_1
        SUM = SUM + RAN( IRANSEED ) -.5
999    CONTINUE
C
C
      RANDOM_GAUSSIAN = SUM
C
C
      RETURN
C
C
      END
C
C
C
      FUNCTION ARCTAN_DEGREES( Y, X)
C     ===============================
C
C
      IF ( X .EQ. 0.0 )  THEN
        IF ( Y .EQ. 0.0 )  A = 0.
        IF ( Y .GT. 0.0 )  A = 90.
        IF ( Y .LT. 0.0 )  A = -90.
      ELSE
        IF ( Y .EQ. 0.0 )  THEN
          IF ( X .GE. 0.0 )  A = 0.
          IF ( X .LT. 0.0 )  A = 180.
        ELSE
          A = ATAN2D( Y, X)
        END IF
      END IF
C
C
      ARCTAN_DEGREES = A
C
C
      RETURN
C
C
      END
C
