C
C     parser.f: wrappers and utilies for C parser functions
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
C          CCP4 PARSER Routines
C          ====================
C
C  Original Author: Based on Mike Levitt's routine of the same name.
C  Modified By: Peter Brick, Phil Evans, Eleanor Dodson, Dave Love
C
C  Martyn Winn: legacy routines retained here. PARSER, PARSE and PARSDL
C    are now wrappers to C functions.
C
C     Library parser.f contains the following subroutines and functions,
C     some of which are currently unsued and commented out.
C
C  SUBROUTINES 
C
C    KEYNUM(N,NSTART,LINE,IBEG,IEND,ITYP,NTOK)
C  [ KEYERR(I,MODE,LINE,IBEG,IEND,ITYP) ] - internal subroutine (KEYNUM)
C    GTNREA(N,M,X,NTOK,ITYP,FVALUE)
C    GTNINT(N,M,J,NTOK,ITYP,FVALUE)
C    GTPREA(N,X,NTOK,ITYP,FVALUE)
C    GTPINT(N,I,NTOK,ITYP,FVALUE)
C    SBLANK(ARRAY,N1,N2)
C    CHKKEY(KEY,WORDS,NWORDS,IKEY)
C    PUTLIN(STROUT,OUTWIN)
C    BLANK(OUTWIN,NLINES)
C    LERROR(ERRFLG,IFAIL,ERRMSG)
C    RDSYMM(JTOK,LINE,IBEG,IEND,ITYP,FVALUE,NTOK,SPGNAM,NUMSGP,PGNAME,
C           NSYM,NSYMP,RSYM)
C    RDHEAD(JTOK,LINE,IBEG,IEND,ITYP,FVALUE,NTOK,MTZPRT,MTZBPR)
C    RDCELL(ITOK,ITYPE,FVALUE,NTOK,CELL)
C    RDRESO(ITOK,ITYPE,FVALUE,NTOK,RESMIN,RESMAX,SMIN,SMAX)
C    RDSCAL(ITOK,LINE,IBEG,IEND,ITYP,FVALUE,NTOK,NLPRGI,LSPRGI,ILPRGI,SCAL,BB)
C    RDRESL(ITOK,ITYPE,FVALUE,CVALUE,NTOK,RESMIN,RESMAX,SMIN,SMAX,ISTAT)
C    RDATOMSELECT(JTOK,INAT0,INAT1,IRES0,IRES1,CHNAM,IMODE,NTOK,LINE,IBEG,
C           IEND,ITYP,IDEC,FVALUE,IFAIL)
C    GTTREA(N,X,LFLAG,NTOK,ITYP,FVALUE)
C    GTTINT(N,I,LFLAG,NTOK,ITYP,FVALUE)
C
C FUNCTIONS
C
C    LOGICAL FUNCTION CMATCH(STRING1,STRING2,NCHAR)
C
C_BEGIN_INTRO
C          CCP4 PARSER Routines
C          ====================
C
C The PARSER module of the CCP4 library contains routines which are
C mainly used for `free-format' `keyworded' input of control data for
C programs.  Most programs have a loop over input records which are
C initially fed to the routine PARSER to tokenise them and extract the
C initial keyword.  PARSER can cope with continued, commented input
C lines and included files.  It calls PARSE to tokenise individual
C records and PARSE is sometimes useful itself to compensate for the
C lack of free-format internal READs in the fortran77 standard.  See
C the entries below for details.
C
C The library also contains routines to decode the parameters
C following the `standard' program keywords SYMMETRY, RESOLUTION,
C SCALE and CELL and to extract real and integer numbers from fields.
C 
C_END_INTRO
C
C
C_BEGIN_KEYNUM
C     ====================================================
      SUBROUTINE KEYNUM(N,NSTART,LINE,IBEG,IEND,ITYP,NTOK)
C     ====================================================
C  Check that correct number of numbers (numeric fields) are present
C
C--- Arguments:
C
C  N      (I) INTEGER        Number of consecutive numeric fields expected
C
C  NSTART (I) INTEGER        Number of first field to check
C
C  LINE   (I) CHARACTER*(*)  Array containing the fields
C
C  IBEG   (I) INTEGER(*)     First column number of fields (from PARSER)
C
C  IEND   (I) INTEGER(*)     Last column number of fields (from PARSER)
C
C  ITYP   (I) INTEGER(*)     =0  null field
C                            =1  character string
C                            =2  number
C                            (from PARSER)
C
C  NTOK   (I) INTEGER        Number of fields (from PARSER)
C
C_END_KEYNUM
C
C     .. Scalar Arguments ..
      INTEGER           N,NSTART,NTOK
      CHARACTER LINE*(*)
C     ..
C     .. Array Arguments ..
      INTEGER           IBEG(*),IEND(*),ITYP(*)
C     ..
C     .. Local Scalars ..
      INTEGER           I
      CHARACTER   LINERR*200
C     ..
C     .. External Subroutines ..
      EXTERNAL          KEYERR, CCPERR
C     ..
C
      DO 10 I = NSTART,NSTART + N - 1
          IF (I.GT.NTOK) THEN
              GO TO 30
          ELSE IF (ITYP(I).NE.2) THEN
              GO TO 20
          END IF
   10 CONTINUE
C
      RETURN
C
C          *******************************
   20 CALL KEYERR(I,2,LINE,IBEG,IEND,ITYP)
C          *******************************
C
      CALL CCPERR(1, 'Keyword error')
   30 CONTINUE
C
          WRITE (LINERR,FMT='(A,I4,A,I4,A)') 
     +     ' *** TOO FEW NUMBERS - ', (I - NSTART),
     +     ' FOUND WHEN ',N,' EXPECTED'
          CALL CCPERR(1, LINERR)
C
      END
C
      SUBROUTINE KEYERR(I,MODE,LINE,IBEG,IEND,ITYP)
C     =============================================
C  Print warning when token not of correct type.
C  Internal subroutine, called from KEYNUM.
C
C     .. Scalar Arguments ..
      INTEGER           I,MODE
      CHARACTER LINE*(*)
C     ..
C     .. Array Arguments ..
      INTEGER           IBEG(*),IEND(*),ITYP(*)
C     ..
C     .. Local Arrays ..
      CHARACTER         TYPE(3)*12
C     ..
C     .. Local Scalars ..
      CHARACTER LINERR*150
C     ..
C     .. External Subroutines ..
      EXTERNAL LERROR
C     ..
C     .. Data statements ..
      DATA              TYPE/'alphanumeric','numeric     ',
     +                  'quoted      '/
C     ..
C
C
      IF (MODE.EQ.0) THEN
          WRITE (LINERR,FMT='(A,A,A)') 
     +  ' ** ERROR : Key word < ',
     +  LINE(IBEG(I) : IEND(I)),
     +  ' > not recognized and has therefore been ignored'
C
C              ****************************
          CALL LERROR(1,0,LINERR)
C              ****************************
C
      ELSE
          WRITE (LINERR,FMT='(A,A,A,A,A,A,A)') 
     + ' ** ERROR: Token < ',
     +  LINE(IBEG(I) : IEND(I)),
     + ' > is ',
     + TYPE(ITYP(I)),
     + ' while a ',TYPE(I),' token was expected'
C
C              ****************************
          CALL LERROR(1,0,LINERR)
C              ****************************
      END IF
C
      END
C
C_BEGIN_GTNREA
C     ========================================
      SUBROUTINE GTNREA(N,M,X,NTOK,ITYP,FVALUE)
C     ========================================
C  Extract M real numbers X starting from N'th value of Parser
C  array FVALUE, if possible. If no value, X = 0.0 .
C  If illegal, write message.
C
C--- Arguments:
C
C N      (I) INTEGER    Number of 1st element of FVALUE to be extracted
C
C M      (I) INTEGER    Number of elements to be extracted
C
C X      (O) REAL(M)    Put extracted elements into this array
C
C NTOK   (I) INTEGER    Total number of fields (from PARSER)
C
C ITYP   (I) INTEGER(*)  =0  null field
C                        =1  character string
C                        =2  number
C
C FVALUE (I) REAL(*)     Array of numbers to be extracted (from PARSER)
C
C_END_GTNREA
C
C     .. Scalar Arguments ..
      INTEGER           M,N,NTOK
C     ..
C     .. Array Arguments ..
      INTEGER           ITYP(*)
      REAL              X(M),FVALUE(*)
C     ..
C     .. Local Scalars ..
      INTEGER           I,K
      CHARACTER LINERR*100
C     ..
C     .. External Subroutines ..
      EXTERNAL LERROR
C     ..
C
      DO 10 I = 1,M
          K = I + N - 1
          X(I) = 0.0
C
          IF (K.LE.NTOK) THEN
              IF (ITYP(K).EQ.2) THEN
                  X(I) = FVALUE(K)
              ELSE IF (ITYP(K).EQ.1) THEN
           WRITE (LINERR,FMT='(A,I4)') 
     +    ' Illegal number in field ',K
C
C              ****************************
          CALL LERROR(1,0,LINERR)
C              ****************************
C
              END IF
          END IF
   10 CONTINUE
      END
C
C_BEGIN_GTNINT
C     ========================================
      SUBROUTINE GTNINT(N,M,J,NTOK,ITYP,FVALUE)
C     ========================================
C Extract M integers J starting from N'th value of Parser array FVALUE,
C if possible. If no value, J = 0 . If illegal, write message
C
C--- Arguments:
C
C N      (I) INTEGER     Number of 1st element of FVALUE to be extracted
C
C M      (I) INTEGER     Number of elements to be extracted
C
C J      (O) INTEGER(M)  Put extracted elements into this array
C
C NTOK   (I) INTEGER     Total number of fields (from PARSER)
C
C ITYP   (I) INTEGER(*)  =0  null field
C                        =1  character string
C                        =2  number
C
C FVALUE (I) REAL(*)     Array of numbers to be extracted (from PARSER)
C
C_END_GTNINT
C
C     .. Scalar Arguments ..
      INTEGER           M,N,NTOK
C     ..
C     .. Array Arguments ..
      INTEGER           J(M),ITYP(*)
      REAL              FVALUE(*)
C     ..
C     .. Local Scalars ..
      INTEGER           I,K
      CHARACTER LINERR*200
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC         NINT
C     ..
C     .. External Subroutines ..
      EXTERNAL LERROR
C     ..
C
      DO 10 I = 1,M
          K = I + N - 1
          J(I) = 0
          IF (K.LE.NTOK) THEN
              IF (ITYP(K).EQ.2) THEN
                  J(I) = NINT(FVALUE(K))
              ELSE IF (ITYP(K).EQ.1) THEN
           WRITE (LINERR,FMT='(A,I4)') 
     +    ' Illegal number in field ',K
C
C              ****************************
          CALL LERROR(1,0,LINERR)
C              ****************************
C
              END IF
          END IF
   10 CONTINUE
C
      END
C
C_BEGIN_GTPREA
C     ======================================
      SUBROUTINE GTPREA(N,X,NTOK,ITYP,FVALUE)
C     ======================================
C Extract real number X from N'th value Parser array FVALUE, if possible
C If no value, leave X unchanged. If illegal, write message
C
C--- Arguments:
C
C N      (I) INTEGER    Number of 1st element of FVALUE to be extracted
C
C X      (O) REAL       Extracted number put here
C
C NTOK   (I) INTEGER    Total number of fields (from PARSER)
C
C ITYP   (I) INTEGER(*)  =0  null field
C                        =1  character string
C                        =2  number
C
C FVALUE (I) REAL(*)     Array of numbers to be extracted (from PARSER)
C
C_END_GTPREA
C
C     .. Scalar Arguments ..
      REAL X
      INTEGER N,NTOK
C     ..
C     .. Array arguments ..
      REAL FVALUE(*)
      INTEGER ITYP(*)
C     ..
C     .. Local Scalars ..
      CHARACTER LINERR*200
C     ..
C     .. External Subroutines ..
      EXTERNAL LERROR
C     ..
C
      IF (N.LE.NTOK) THEN
        IF (ITYP(N).EQ.2) THEN
          X = FVALUE(N)
        ELSE IF (ITYP(N).EQ.1) THEN
           WRITE (LINERR,FMT='(A,I4)') 
     +    ' Illegal number in field ',N
C
C              ****************************
          CALL LERROR(1,0,LINERR)
C              ****************************
C
        END IF
      ELSE
        CALL LERROR (1, 0, 'Real number expected at end of line')
      END IF
C
      END
C
C_BEGIN_GTPINT
C     ======================================
      SUBROUTINE GTPINT(N,I,NTOK,ITYP,FVALUE)
C     ======================================
C Extract integer I from N'th value Parser array FVALUE, if possible
C If no value, leave I unchanged. If illegal, write message
C
C--- Arguments:
C
C N      (I) INTEGER    Number of 1st element of FVALUE to be extracted
C
C I      (O) INTEGER    Extracted number put here
C
C NTOK   (I) INTEGER    Total number of fields (from PARSER)
C
C ITYP   (I) INTEGER(*)  =0  null field
C                        =1  character string
C                        =2  number
C
C FVALUE (I) REAL(*)     Array of numbers to be extracted (from PARSER)
C
C_END_GTPINT
C
C     .. Scalar Arguments ..
      INTEGER I,N,NTOK
C     ..
C     .. Arrays arguments ..
      REAL FVALUE(*)
      INTEGER ITYP(*)
C     ..
C     .. Local Scalars ..
      INTEGER ISTERR,IFGERR
      CHARACTER LINERR*100
C     ..
C     .. External Subroutines ..
      EXTERNAL LERROR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC NINT
C     ..
C
      IF (N.LE.NTOK) THEN
        IF (ITYP(N).EQ.2) THEN
          I = NINT(FVALUE(N))
        ELSE IF (ITYP(N).EQ.1) THEN
           WRITE (LINERR,FMT='(A,I4)') 
     +    ' Illegal number in field ',N
          ISTERR = 1
          IFGERR = 0
C
C              ****************************
          CALL LERROR(ISTERR,IFGERR,LINERR)
C              ****************************
C
        END IF
      END IF
C
      END
C
C_BEGIN_SBLANK
C     ==============================
      SUBROUTINE SBLANK(ARRAY,N1,N2)
C     ==============================
C Blank characters N1 to N2 of ARRAY
C
C--- Arguments:
C
C ARRAY (I/O)  CHARACTER(*)
C
C N1    (I)    INTEGER
C
C N2    (I)    INTEGER
C
C_END_SBLANK
C
      CHARACTER*1 ARRAY(*)
      INTEGER I,N1,N2
C
      DO 10 I=N1,N2
         ARRAY(I)=' '
10     CONTINUE
C
      RETURN
      END
C
C_BEGIN_CMATCH
C     ==============================================
      LOGICAL FUNCTION CMATCH(STRING1,STRING2,NCHAR)
C     ==============================================
C
C---- Compare nchar character in string1 and string2
C     return cmatch .true. if all match, else .false.
C
C---- Arguments:
C
C STRING1 (I) CHARACTER*(*)  1st string to compare
C
C STRING2 (I) CHARACTER*(*)  2nd string to compare
C
C NCHAR   (I) INTEGER        number of characters to compare
C
C_END_CMATCH
C
      CHARACTER*(*) STRING1,STRING2
      INTEGER NCHAR
C
      IF(STRING1(1:NCHAR).EQ.STRING2(1:NCHAR)) THEN
          CMATCH=.TRUE.
      ELSE
          CMATCH=.FALSE.
      ENDIF
      END
C
C_BEGIN_CHKKEY
C     ========================================
      SUBROUTINE CHKKEY(KEY,WORDS,NWORDS,IKEY)
C     ========================================
C Check keyword KEY against list of NWORDS possible keywords in WORDS.
C Allows abbreviated or extended keys provided they are not ambiguous.
C
C---- Arguments:
C
C KEY    (I) CHARACTER*(*)         Keyword for checking
C
C WORDS  (I) CHARACTER(NWORDS)*(*) List of possible keywords
C
C NWORDS (I) INTEGER               Number of keywords in WORDS
C
C IKEY (I/O) INTEGER               = '?', list all words
C                                  Returns:
C                                  = keyword number found (.gt.0)
C                                  = 0 if not found or null
C                                  = -1 if ambiguous
C
C_END_CHKKEY
C
      INTEGER NFMAX
      PARAMETER (NFMAX=20)
C     .. Scalar Arguments ..
      INTEGER NWORDS, IKEY
      CHARACTER KEY*(*)
C     ..
C     .. Array Arguments ..
      CHARACTER WORDS(NWORDS)*(*)
C     ..
C     .. Local Scalars ..
      INTEGER LK,I,L,NFOUND,JDO
      CHARACTER LINERR*200
C     ..
C     .. Local Arrays ..
      INTEGER LFOUND(20)
C     ..
C     .. External Functions ..
      INTEGER LENSTR
      EXTERNAL LENSTR
C     ..
C     .. External Subroutines ..
      EXTERNAL PUTLIN,LERROR
C     ..
C
C---- Get minimum significant length of KEY 
C     ( function LENSTR returns the length
C     of the character string excluding trailing blanks)
C
      IKEY=0
C
C        ***********
      LK=LENSTR(KEY)
C        ***********
C
C---- Ignore null string
C
      IF(LK.LE.0) RETURN
C
      IF(KEY(1:1).EQ.'?') THEN
C  
C           ****************
       CALL PUTLIN(' Possible keywords are:','HLPWIN')
C           ****************
C
       DO 10 JDO = 1,NWORDS
C
C            ****************
        CALL PUTLIN(WORDS(JDO),'HLPWIN')
C            ****************
C
10      CONTINUE
C
            IKEY=0
            RETURN
      ENDIF
C
      NFOUND=0
C
C---- Check all possible words in case of ambiguities
C
      DO 20 I=1,NWORDS
C
C----  Key may be longer than word in list
C
C              ****************
      L=MIN(LK,LENSTR(WORDS(I)))
C              ****************
C
      IF(L.LE.0) GO TO 20
C
C---- Find out if KEY is an initial substring of this option word
C
      IF(INDEX(WORDS(I),KEY(1:L)).EQ.1) THEN
            NFOUND=NFOUND+1
C
            IF(NFOUND.GT.NFMAX) THEN
          WRITE (LINERR,FMT='(A,I5)') 
     +  ' CHKKEY: too many ambiguities : ',NFMAX
C
C              ****************************
          CALL LERROR(1,0,LINERR)
C              ****************************
C
                  NFOUND=NFMAX
            ELSE
                  LFOUND(NFOUND)=I
            ENDIF
       ENDIF
20     CONTINUE
C
C---- If keyword is ambiguous, list possibilities
C
      IF(NFOUND.GT.1) THEN
          WRITE (LINERR,FMT='(A,A,A)') 
     +   ' Keyword ',
     +   KEY(1:LK),
     +  ' is ambiguous: possibilities are -'
C
C              ****************************
          CALL LERROR(1,0,LINERR)
C              ****************************
C
       DO 30 JDO = 1,NWORDS
C
C            ****************
        CALL PUTLIN(WORDS(JDO),'HLPWIN')
C            ****************
C
30      CONTINUE
            IKEY=-1
      ELSEIF (NFOUND.EQ.1) THEN
C
C---- Success if only 1 found
C
            IKEY=LFOUND(1)
      ENDIF
      END
C
C_BEGIN_PUTLIN
C     ================================
      SUBROUTINE PUTLIN(STROUT,OUTWIN)
C     ================================
C---- This is a dummy PUTLIN to link with the MTZ routines mark 1 -
C     all it does is write the line in STROUT to lun 6. Later the
C     routines will be linked with the Compose-Parser etc. from Kim
C     where PUTLIN does a few more things !
C
C---- Arguments:
C
C STROUT (I) CHARACTER*(*)  Input line
C
C OUTWIN (O) CHARACTER*(*)  Not used
C
C_END_PUTLIN
C
C     .. Scalar Arguments ..
      CHARACTER OUTWIN* (*)
      CHARACTER STROUT* (*)
C     ..
C     .. Local Scalars ..
      INTEGER LUNOUT,LL,LS,LX
C     ..
C     .. External Functions ..
      INTEGER LENSTR
      EXTERNAL LENSTR
C     ..
C     .. Data statements ..
      DATA LUNOUT/6/
C     ..
C
      LL = LENSTR(STROUT)
      IF (LL.GE.132) THEN
      LX = 1
      LS = 131
10    CONTINUE
      WRITE (LUNOUT,FMT=6000) STROUT(LX:LS)
      IF (LS.EQ.LL) GOTO 20
      LX = LS  + 1
      LS = LS + 130
      IF (LS.GT.LL) LS = LL
      GO TO 10
      ELSE
       IF (LL.EQ.0) THEN
           WRITE(LUNOUT,FMT=6000)
            ELSE
      WRITE (LUNOUT,FMT=6000) STROUT(1:LL)
           END IF
      END IF
20    CONTINUE
C
C---- Format statements
C
 6000 FORMAT (' ',A)
C
      END
C
C_BEGIN_BLANK
C     ===============================
      SUBROUTINE BLANK(OUTWIN,NLINES)
C     ===============================
C---- This subroutine calls PUTLIN to output NLINES blank lines to the
C     window OUTWIN
C
C---- Arguments:
C
C     OUTWIN  (I)   CHARACTER*6     output window
C
C     NLINES  (I)   INTEGER         number of blank lines to output
C
C_END_BLANK
C
C     .. Scalar Arguments ..
      INTEGER NLINES
      CHARACTER OUTWIN*(*)
C     ..
C     .. Local Scalars ..
      INTEGER JDO10
C     ..
C     .. External Subroutines ..
      EXTERNAL PUTLIN
C     ..
C
      DO 10 JDO10 = 1,MAX(NLINES,1)
C
C            **************
        CALL PUTLIN(' ',OUTWIN)
C            **************
C
   10 CONTINUE
      END
C
C_BEGIN_LERROR
C     =======================================
      SUBROUTINE LERROR(ERRFLG,IFAIL,ERRMSG)
C     =======================================
C---- General error reporting subroutine, for the MTZ routines, etc
C
C---- Arguments:
C
C     ERRFLG  (I)  INTEGER         =1 output meesage as warning
C                                  =2 output message as fatal
C
C     IFAIL   (I)  INTEGER         =0 return after fatal error
C                                  =-1 STOP after reporting fatal error
C
C     ERRMSG  (I)  CHARACTER*(*)   character string containing error
C                                  message to output
C
C_END_LERROR
C
C     .. Scalar Arguments ..
      INTEGER ERRFLG,IFAIL
      CHARACTER ERRMSG* (*)
C     ..
C     ..
C     .. External Subroutines ..
      EXTERNAL BLANK,PUTLIN
C     ..
C
      IF (ERRFLG.EQ.1) THEN
C
C---- Output a warning message and return
C
        CALL BLANK('ERRWIN',1)
        CALL PUTLIN('***  Warning','ERRWIN')
        CALL PUTLIN(ERRMSG,'ERRWIN')
        CALL BLANK('ERRWIN',1)
C
      ELSE IF (ERRFLG.EQ.2) THEN
C
C---- Output a fatal message, and quit or return depending on IFAIL
C
        CALL BLANK('ERRWIN',1)
        CALL PUTLIN('***  Error','ERRWIN')
        CALL PUTLIN(ERRMSG,'ERRWIN')
        IF (IFAIL.LT.0) THEN
          call ccperr(1,'*** Program Terminated ')
        ELSE
          CALL BLANK('ERRWIN',1)
        END IF
        RETURN
      ELSE
C
C---- Bad errflg, output message and continue
C
        CALL BLANK('ERRWIN',1)
        CALL PUTLIN('*** Unrecognised  error','ERRWIN')
        CALL PUTLIN(ERRMSG,'ERRWIN')
        CALL PUTLIN('Program continuing ...','ERRWIN')
        CALL BLANK('ERRWIN',1)
C
      END IF
      END
C
C
C_BEGIN_RDSYMM
C     =======================================================
      SUBROUTINE RDSYMM(JTOK,LINE,IBEG,IEND,ITYP,FVALUE,NTOK,
     .    SPGNAM,NUMSGP,PGNAME,NSYM,NSYMP,RSYM)
C     =======================================================
C---- Read and decode symmetry specification
C
C---- Arguments:
C
C   JTOK    (I)  INTEGER        Number of first field to interpret
C
C   LINE    (I)  CHARACTER*(*)  Input string (from PARSER)
C
C   IBEG    (I)  INTEGER(*)     1st column number of tokens in field 
C                               (from PARSER)
C
C   IEND    (I)  INTEGER(*)     Last column number of tokens in field
C                               (from PARSER)
C
C   ITYP    (I)  INTEGER(*)     =0  null field
C                               =1  character string
C                               =2  number
C                               (from PARSER)
C
C   FVALUE  (I)  REAL(*)        Array of numbers. (from PARSER)
C
C   NTOK    (I)  INTEGER        The number of fields parsed. (from PARSER)
C
C     
C   NSYM  (I/O)  INTEGER        Number of symmetry operations already read,
C                               including non-primitive.
C                               (should be cleared to 0 at beginning)
C
C   SPGNAM  (I/O) CHARACTER*(*)   Space group name
C
C   NUMSGP  (O) INTEGER         Space group number
C
C   PGNAME  (O) CHARACTER*(*)   Point group name
C
C   NSYMP   (O) INTEGER         Number of primitive symmetry operations
C
C   RSYM    (O) REAL(4,4,*)     Symmetry matrices. * should be at least =NSYM
C
C_END_RDSYMM
C     
      INTEGER JTOK,NTOK
      INTEGER IBEG(*),IEND(*),ITYP(*)
      REAL FVALUE(*)
      CHARACTER*(*)LINE,SPGNAM,PGNAME
      CHARACTER*20 SPGNAMS
      INTEGER NUMSGP,NSYM,NSYMP,LENSTR
      EXTERNAL LENSTR

      REAL RSYM(4,4,*)
C     
C---- Look at next field on line: this can be
C     (a) a space-group number
C     (b) a space-group name, ie a string beginning P,I,R,F,A,B or C
C     (c) a symmetry operation (anything else)
C     
C---- for cases (a) & (b), this is a single field:
C     case (c) is more than 1 field
C     
      SPGNAMS = ' '
      SPGNAM = ' '
      IF (JTOK.GT.NTOK) THEN
         CALL  PUTLIN(' No symmetry data !!!','CURWIN')
      ELSE
         IF (JTOK.EQ.NTOK) THEN
            IF (NSYM.GT.0) THEN
               CALL  PUTLIN('Warning: symmetry already given','CURWIN')
            ENDIF
C     
C---- A single field, see if it is a number or a string
C
            IF (ITYP(JTOK).EQ.2) THEN
C     
C---- it's a number, treat as space-group number
C     
               NUMSGP = NINT(FVALUE(JTOK))
            ELSE
C     
C---- it's a string, treat as space-group name
C     
               SPGNAM = LINE(IBEG(JTOK) :IEND(JTOK))
               CALL  CCPUPC(SPGNAM)
               IF (SPGNAM(1:1).EQ.'R' .AND. 
     +             INDEX(SPGNAM,':H').LE.0) THEN
                 WRITE(6,'(A,A)') 'Warning: rhombohedral axes implied',
     +                  ' - if you have hexagonal axes then use H'
               ENDIF
               NUMSGP = 0
            END IF
C     
C---- Read symmetry (all operations) from SYMOP
C     open symop on channel 24 - closed at end of reading
C     NSYMP returns number of primitive operations
C     
C           CALL  MSYMLB(24,NUMSGP,SPGNAM,PGNAME,NSYMP,NSYM,RSYM)
            nsymp = 0
            nsym = 0 
            CALL  MSYMLB3(24,NUMSGP,SPGNAM,SPGNAMS,PGNAME,NSYMP,NSYM,
     +                    RSYM)
         ELSE
C     
C     
C---- Read symmetry operations
C    
            NSYM = NSYM + 1
            NSYMP = NSYM
            CALL  CCPUPC(LINE)
            CALL  SYMFR2(LINE,IBEG(JTOK),NSYM,RSYM)
            NUMSGP = 0
            PGNAME = ' '
C     
         END IF
      END IF
      END     
C     
C_BEGIN_RDHEAD
C     ======================================================
      SUBROUTINE RDHEAD(JTOK,LINE,IBEG,IEND,ITYP,FVALUE,NTOK,
     .    MTZPRT,MTZBPR)
C     ======================================================
C---- Read and decode HEADER command, to set print flags for MTZ headers
C
C---- Arguments:
C 
C   JTOK   (I) INTEGER       Number of first field to interpret
C
C   LINE   (I) CHARACTER*(*) Input string (from PARSER)
C
C   IBEG   (I) INTEGER(*)    1st column number of tokens in field 
C                            (from PARSER)
C
C   IEND   (I) INTEGER(*)    Last column number of tokens in field
C                            (from PARSER)
C
C   ITYP   (I) INTEGER(*)    =0  null field
C                            =1  character string
C                            =2  number
C                            (from PARSER)
C
C   FVALUE (I) REAL(*)       Array of numbers. (from PARSER)
C
C   NTOK   (I) INTEGER       The number of fields parsed. (from PARSER)
C
C     
C   MTZPRT (O) INTEGER       Flag to control printout from MTZ file header
C                            NONE    sets MTZPRT = 0
C                             no header o/p
C                            BRIEF   sets MTZPRT = 1 (default)
C                             brief header o/p
C                            HISTORY sets MTZPRT = 2
C                             brief + mtz history
C                            ALL     sets MTZPRT = 3
C                             full header o/p from mtz reads
C     
C   MTZBPR (O) INTEGER       Controls printout from BATCH HEADERS
C                            NOBATCH     sets MTZBPR = 0
C                             no batch header o/p
C                            BATCH       sets MTZBPR = 1  (default)
C                             batch titles o/p
C                            ORIENTATION sets MTZBPR = 2
C                             batch orientation also
C
C_END_RDHEAD
C     
      INTEGER JTOK,NTOK
      INTEGER IBEG(*),IEND(*),ITYP(*)
      REAL FVALUE(*)
      CHARACTER*(*) LINE
      INTEGER MTZPRT,MTZBPR
C     
C     Locals
      INTEGER I,IKEY
      CHARACTER KEY*12
C     
      INTEGER NKEYS
      PARAMETER (NKEYS=7)
      CHARACTER*12 KEYS(NKEYS)
      DATA KEYS/'NONE','BRIEF','HISTORY','ALL',
     $     'NOBATCH','BATCH','ORIENTATION'/
C     
C     Set defaults
      MTZPRT = 1
      MTZBPR = 1
C     
C     Loop keywords
      IF (NTOK .GE. JTOK) THEN
         DO 10, I=JTOK,NTOK
            KEY = LINE(IBEG(I):IEND(I))
            CALL CCPUPC(KEY)
            CALL CHKKEY(KEY,KEYS,NKEYS,IKEY)
            IF (IKEY .LE. 0) THEN
              CALL PUTLIN
     +             ('Unrecognized or ambiguous subkeyword to HEADER: '
     +             // KEY,'CURWIN')
            ELSE
               IF (IKEY .EQ. 1) MTZPRT = 0
               IF (IKEY .EQ. 2) MTZPRT = 1
               IF (IKEY .EQ. 3) MTZPRT = 2
               IF (IKEY .EQ. 4) MTZPRT = 3
               IF (IKEY .EQ. 5) MTZBPR = 0
               IF (IKEY .EQ. 6) MTZBPR = 1
               IF (IKEY .EQ. 7) MTZBPR = 2
            ENDIF
 10      CONTINUE 
      ENDIF
      END
C     
C_BEGIN_RDCELL
C     ==============================================
      SUBROUTINE RDCELL(ITOK,ITYPE,FVALUE,NTOK,CELL)
C     ==============================================     
C---- Read and decode cell parameters 
C     
C---- Arguments:
C
C   ITOK   (I) INTEGER     Number of first field to interpret
C
C   ITYPE  (I) INTEGER(*)  =0  null field
C                          =1  character string
C                          =2  number
C                          (from PARSER)
C
C   FVALUE (I) REAL(*)     Array of numbers. (from PARSER)
C
C   NTOK   (I) INTEGER     The number of fields parsed. (from PARSER)
C
C   CELL   (O) REAL(6)     Cell parameters a, b, c, alpha, beta, gamma.
C
C_END_RDCELL
C     
C     .. Scalar Arguments ..
      INTEGER           ITOK,NTOK
C     ..
C     .. Array Arguments ..
      REAL              CELL(6),FVALUE(*)
      INTEGER           ITYPE(*)
C     ..
C     .. External Subroutines ..
      EXTERNAL          GTPREA
C     ..
C     
      IF (NTOK .LT. ITOK+2) THEN
        CALL LERROR (1,0,'Cell a, b and c not given -- ignored')
        RETURN
      END IF
      CELL(4) = 90.0
      CELL(5) = 90.0
      CELL(6) = 90.0
C     
C     ***************************************
      CALL GTPREA(ITOK,CELL(1),NTOK,ITYPE,FVALUE)
      CALL GTPREA(ITOK+1,CELL(2),NTOK,ITYPE,FVALUE)
      CALL GTPREA(ITOK+2,CELL(3),NTOK,ITYPE,FVALUE)
C     ***************************************
C     
C     *********************************************
      IF (ITOK+3.LE.NTOK) CALL GTPREA(ITOK+3,CELL(4),NTOK,ITYPE,FVALUE)
      IF (ITOK+4.LE.NTOK) CALL GTPREA(ITOK+4,CELL(5),NTOK,ITYPE,FVALUE)
      IF (ITOK+5.LE.NTOK) CALL GTPREA(ITOK+5,CELL(6),NTOK,ITYPE,FVALUE)
C     *********************************************
      END
C     
C     
C_BEGIN_RDRESO
C     ================================================
      SUBROUTINE RDRESO(ITOK,ITYPE,FVALUE,NTOK,RESMIN,
     +                  RESMAX,SMIN,SMAX)
C     ================================================
C---- Read and decode resolution limits.
C     
C---- Arguments:
C
C     ITOK    (I) INTEGER     Number of first field to interpret
C     
C     ITYPE   (I) INTEGER(*)  =0  null field
C                             =1  character string
C                             =2  number
C                             (from PARSER)
C
C     FVALUE  (I) REAL(*)     Array of numbers. (from PARSER)
C
C     NTOK    (I) INTEGER     The number of fields parsed. (from PARSER)
C
C     
C     RESMIN  (O) REAL        Minimum resolution (in As)
C
C     RESMAX  (O) REAL        Maximum resolution (in As)
C
C     SMIN    (O) REAL        Minimum resolution ( 4sin**2/lambda**2)
C
C     SMAX    (O) REAL        Maximum resolution ( 4sin**2/lambda**2)
C
C_END_RDRESO
C     .. Scalar Arguments ..
      REAL              RESMAX,RESMIN,SMAX,SMIN
      INTEGER           ITOK,NTOK
C     ..
C     .. Array Arguments ..
      REAL              FVALUE(*)
      INTEGER           ITYPE(*)
C     ..
C     .. Local Scalars ..
      REAL              RESTEM,STEM
C     ..
C     .. External Subroutines ..
      EXTERNAL          GTPREA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     ..
C
C---- Global defaults set here
C
        RESMAX = 10000.0
        RESMIN = .1
        IF (NTOK.LT.ITOK) THEN
          CALL LERROR(1,0,'Resolution missing')
          RETURN
        END IF
C     
C---- Look at next field on line: this can be
C     read resolution limits in A, if only one treat as high
C     resolution limit
C     
C     *************************************
      CALL GTPREA(ITOK,RESMIN,NTOK,ITYPE,FVALUE)
C     *************************************
C     
      IF (ABS(RESMIN).LE.0.000001) RESMIN = 0.00001
C     
C     ***************************************
      IF (NTOK.GT.ITOK)
     +      CALL GTPREA(ITOK+1,RESMAX,NTOK,ITYPE,FVALUE)
C     ***************************************
C     
      IF (ABS(RESMAX).LE.0.0000001) RESMAX = 100.0
C     
C     
      IF (RESMIN.LE.RESMAX) THEN
         RESTEM = RESMAX
         RESMAX = RESMIN
         RESMIN = RESTEM
      END IF
C     
C---- option to read 4sin**2/lamda**2
C     
      IF (RESMIN.LE.1.0 .AND. RESMAX.LE.1.0) THEN
C     
C---- swap over smin and resmin etc
C     
         SMIN = RESMIN
         SMAX = RESMAX
         RESMAX = SQRT(1.0/SMIN)
         RESMIN = SQRT(1.0/SMAX)
      ELSE
         SMIN = 1.0/RESMAX**2
         SMAX = 1.0/RESMIN**2
      END IF
C     
C     
      IF (SMIN.GT.SMAX) THEN
         STEM = SMAX
         SMAX = SMIN
         SMIN = STEM
      END IF
C     
C     
      RETURN
      END
C     
C     
C_BEGIN_RDSCAL     
C     ======================================================
      SUBROUTINE RDSCAL(ITOK,LINE,IBEG,IEND,ITYP,FVALUE,NTOK,
     .    NLPRGI,LSPRGI,ILPRGI,SCAL,BB)
C     ======================================================
C---- Read and decode SCALE .
C     
C---- Arguments:
C
C  ITOK   (I/O) INTEGER     Input: number of first field to interpret
C                           Output: number of next token to interpret (.gt. 0)
C                                  =  0 if line exhausted (SCAL & BB OK)
C                                  = -1 if no scale given
C                                  = -2 unrecognized label
C
C  LINE   (I) CHARACTER*(*) Input string (from PARSER)
C
C  IBEG   (I) INTEGER(*)    1st column number of tokens in field 
C                           (from PARSER)
C
C  IEND   (I) INTEGER(*)    Last column number of tokens in field
C                           (from PARSER)
C
C  ITYP   (I) INTEGER(*)    =0  null field
C                           =1  character string
C                           =2  number
C                           (from PARSER)
C
C  FVALUE (I) REAL(*)       Array of numbers. (from PARSER)
C
C  NTOK   (I) INTEGER       The number of fields parsed. (from PARSER)
C
C  LSPRGI (I) CHARACTER(*)*30  Program label strings.
C                                  L(abel) S(tring) PRG(ram) I(nput)
C
C  NLPRGI (I) INTEGER        Number of label strings in LSPRGI
C
C  ILPRGI (O) INTEGER        Number in array of LSPRGI whose scale has been reset
C
C  SCAL   (O) REAL           Scale factor, no default
C
C  BB     (O) REAL           Temperature factor, default = 0.0
C
C_END_RDSCAL
C     
      INTEGER ITOK,NTOK,ILPRGI,NLPRGI,JDO
      INTEGER IBEG(*),IEND(*),ITYP(*)
      REAL FVALUE(*)
      CHARACTER*(*) LINE
      CHARACTER*30 LSPRGI(*),CWORK
      REAL SCAL,BB
C     
      CWORK = LINE(IBEG(ITOK) :IEND(ITOK))
      DO 10 JDO = 1,NLPRGI
C     
         IF (CWORK.EQ.LSPRGI(JDO)) GO TO 20
C     
 10   CONTINUE
C     
C     ***********************
      CALL PUTLIN('**** Error input assignment does not match'//
     +     ' program labels','ERRWIN')
C     ***********************
C     
      ITOK = -2
      RETURN
C     
 20   ILPRGI = JDO
      IF(ITOK+1.GT.NTOK) THEN
         ITOK = -1
         RETURN
      ELSE
         IF (ITYP(ITOK+1) .EQ. 2) THEN
            CALL GTPREA(ITOK+1,SCAL,NTOK,ITYP,FVALUE)
         ELSE
            ITOK = -1
            RETURN
         ENDIF
      ENDIF
C
      BB = 0
      IF(ITOK+2.LE.NTOK) THEN
         IF (ITYP(ITOK+2) .EQ. 2) THEN
            CALL GTPREA(ITOK+2,BB,NTOK,ITYP,FVALUE)
            ITOK = ITOK + 3
         ELSE
            ITOK = ITOK + 2
         ENDIF
         IF (ITOK .GT. NTOK) ITOK = 0
      ELSE
         ITOK = 0
      ENDIF
C     
      RETURN
      END 
C
C
C_BEGIN_RDRESL
C     ======================================================
      SUBROUTINE RDRESL(ITOK,ITYPE,FVALUE,CVALUE,NTOK,RESMIN,
     +                  RESMAX,SMIN,SMAX,ISTAT)
C     ======================================================     
C---- Read and decode resolution limits.
C     Subkeywords in CVALUE recognized:
C       LOW   read next number as low resolution limit
C       HIGH  read next number as high resolution limit
C
C     If LOW & HIGH are both present, the limits will still be swapped
C     to the correct order
C
C     If only LOW or HIGH are given, the unset limit (ie either RESMAX, SMAX
C     or RESMIN, SMIN) will be set to -1.0. If only one number is given,
C     it is treated as a high resolution limit
C
C     If both limits are given without keywords, and both are .lt. 1.0,
c     it is assumed that the limits are 4(sin theta/lambda)**2 rather than A
C
C---- Arguments:
C
C  ITOK   (I) INTEGER         Number of first field to interpret
C
C  ITYP   (I) INTEGER(*)      =0  null field
C                             =1  character string
C                             =2  number
C                             (from PARSER)
C
C  FVALUE (I) REAL(*)         Array of numbers. (from PARSER)
C
C  NTOK   (I) INTEGER         The number of fields parsed. (from PARSER)
C
C  CVALUE (I) CHARACTER(*)*4  Parsed tokens from program input. (from PARSER)
C
C  RESMIN  (O) REAL           Minimum resolution (in As) (ie low resolution)
C
C  RESMAX  (O) REAL           Maximum resolution (in As) (ie high resolution)
C
C  SMIN    (O) REAL           Minimum resolution ( 4sin**2/lambda**2)
C                                (ie low resolution)
C
C  SMAX    (O) REAL           Maximum resolution ( 4sin**2/lambda**2)
C                                (ie high resolution)
C
C  ISTAT   (O) INTEGER        =0  OK
C                             =-1 illegal subkeyword
C                             =+1 no limits set
C                             =+2 illegal number (probably can't happen)
C_END_RDRESL
C     
C     .. Scalar Arguments ..
      REAL              RESMAX,RESMIN,SMAX,SMIN
      INTEGER           ITOK,NTOK,ISTAT
C     ..
C     .. Array Arguments ..
      REAL              FVALUE(*)
      INTEGER           ITYPE(*)
      CHARACTER*4       CVALUE(*)
C     ..
C     .. Local Scalars ..
      REAL              RESTEM,STEM
      INTEGER           N, KMNMX, NSET, LFLAG, NKEYS, IKEY
      LOGICAL           BOTH, KEYWRD
      CHARACTER*4 SUBKEY(2)
C     ..
C     .. External Subroutines ..
      EXTERNAL          GTTREA, CCPUPC, CHKKEY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     ..
      DATA SUBKEY /'LOW', 'HIGH'/
      DATA NKEYS /2/
C
C---- Global defaults set here
C
      RESMAX = -1.0
      RESMIN = -1.0
      SMIN   = -1.0
      SMAX   = -1.0
      NSET   = 0
      KMNMX  = 1
      ISTAT  = 0
      BOTH = .TRUE.
      KEYWRD = .FALSE.
C     
      N  = ITOK
C
 1    IF (N .LE. NTOK) THEN
C
         IF (ITYPE(N) .EQ. 1) THEN
C String
            CALL CCPUPC(CVALUE(N))
            CALL CHKKEY(CVALUE(N),SUBKEY,NKEYS,IKEY)
            IF(IKEY.LE.0) THEN
               ISTAT = -1
               RETURN
C              ======
            ELSEIF (IKEY .EQ. 1) THEN
C----- subkey LOW
               KMNMX = 1
            ELSEIF (IKEY .EQ. 2) THEN
C----- subkey HIGH
               KMNMX = 2
            ENDIF
            BOTH = .NOT. BOTH
            KEYWRD = .TRUE.
         ELSE
C Number
            RESTEM = 0.0
C                ******************************************
            CALL GTTREA(N,RESTEM,LFLAG,NTOK,ITYPE,FVALUE)
C                ******************************************
            IF (LFLAG .EQ. 0) THEN
               IF (KMNMX .EQ. 1) THEN
                  RESMIN = RESTEM
                  NSET   = NSET+1
                  KMNMX  = 2
               ELSEIF (KMNMX .EQ. 2) THEN
                  RESMAX = RESTEM
                  NSET = NSET+1
                  KMNMX  = 1
               ENDIF
            ELSE
               ISTAT = +2
            ENDIF
         ENDIF
         N = N+1
         GO TO 1
      ENDIF
C
C  Have any numbers been set?
      IF (NSET .EQ. 0) THEN
         ISTAT = +1
         RETURN
C        ======
      ELSEIF (NSET .EQ. 1) THEN
C One only set, if no keywords have been defined, use single number as
C     high resolution limit
         IF (BOTH) THEN
            RESMAX = RESMIN
            RESMIN = -1.0
         ENDIF
      ENDIF
C     
C---- option to read 4sin**2/lamda**2
      IF (.NOT. KEYWRD .AND. NSET .EQ. 2) THEN
         IF (RESMIN .GT. 0.0 .AND. RESMIN .LE. 1.0 .AND.
     $       RESMAX .GT. 0.0 .AND. RESMAX .LE. 1.0) THEN
C---- swap over SMIN and RESMIN
C     
            SMIN = RESMIN
            RESMIN = SQRT(1.0/SMIN)
C---- swap over SMAX and RESMAX 
            SMAX = RESMAX
            RESMAX = SQRT(1.0/SMAX)
         END IF
      ENDIF
C     
      IF (RESMIN .GT. 0.0) THEN
         SMIN = 1.0/RESMIN**2
      END IF
C     
      IF (RESMAX .GT. 0.0) THEN
            SMAX = 1.0/RESMAX**2
      ENDIF
C     
C---- Check that they are in the correct order, if both limits read
C     
      IF (NSET .EQ. 2) THEN
         IF (RESMIN.LE.RESMAX) THEN
            RESTEM = RESMAX
            RESMAX = RESMIN
            RESMIN = RESTEM
         ENDIF
         IF (SMIN.GT.SMAX) THEN
            STEM = SMAX
            SMAX = SMIN
            SMIN = STEM
         ENDIF
      ENDIF
      END
C
C_BEGIN_GTTREA
C     =============================================
      SUBROUTINE GTTREA(N,X,LFLAG,NTOK,ITYP,FVALUE)
C     =============================================
C---- Extract real number X from N'th value of Parser array FVALUE,
C     if possible.
C
C     If no value, leave X unchanged. If illegal, write message
C
C---- Arguments:
C
C  N      (I) INTEGER     Number of 1st element of FVALUE to be extracted
C
C  X      (O) REAL        Put extracted number here
C
C  LFLAG  (O) INTEGER     =  0  OK (valid number or null field)
C                         = -1  beyond end of line
C                         = +1  illegal number
C
C  NTOK   (I) INTEGER     Total number of fields (from PARSER)
C
C  ITYP   (I) INTEGER(*)  =0  null field
C                         =1  character string
C                         =2  number
C                         (from PARSER)
C
C  FVALUE (I) REAL(*)     Array of numbers to be extracted (from PARSER)
C
C_END_GTTREA
C
C     .. Scalar Arguments ..
      REAL X
      INTEGER N,NTOK,LFLAG
C     ..
C     .. Array arguments ..
      REAL FVALUE(*)
      INTEGER ITYP(*)
C     ..
C     .. Local Scalars ..
      CHARACTER LINERR*200
C     ..
C     .. External Subroutines ..
      EXTERNAL LERROR
C     ..
C
      LFLAG = 0
      IF (N.LE.NTOK) THEN
        IF (ITYP(N).EQ.2) THEN
          X = FVALUE(N)
        ELSE IF (ITYP(N).EQ.1) THEN
           WRITE (LINERR,FMT='(A,I4)') 
     +    ' Illegal number in field ',N
C
C              ****************************
          CALL LERROR(1,0,LINERR)
C              ****************************
C
          LFLAG = +1
        END IF
      ELSE
         LFLAG = -1
      END IF
C
      END
C
C_BEGIN_GTTINT
C     =============================================
      SUBROUTINE GTTINT(N,I,LFLAG,NTOK,ITYP,FVALUE)
C     =============================================
C---- Extract integer I from N'th value of Parser array FVALUE,
C     if possible.
C
C     If no value, leave I unchanged. If illegal, write message.
C
C---- Arguments:
C
C  N      (I) INTEGER     Number of 1st element of FVALUE to be extracted
C
C  I      (O) INTEGER     Put extracted number here
C
C  LFLAG  (O) INTEGER     =  0  OK (valid number or null field)
C                         = -1  beyond end of line
C                         = +1  illegal number
C
C  NTOK   (I) INTEGER     Total number of fields (from PARSER)
C
C  ITYP   (I) INTEGER(*)  =0  null field
C                         =1  character string
C                         =2  number
C                         (from PARSER)
C
C  FVALUE (I) REAL(*)     Array of numbers to be extracted (from PARSER)
C
C_END_GTTINT
C
C      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER I,N,NTOK,LFLAG
C     ..
C     .. Arrays arguments ..
      REAL FVALUE(*)
      INTEGER ITYP(*)
C     ..
C     .. Local Scalars ..
      CHARACTER LINERR*100
C     ..
C     .. External Subroutines ..
      EXTERNAL LERROR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC NINT
C     ..
C
      LFLAG = 0
      IF (N.LE.NTOK) THEN
        IF (ITYP(N).EQ.2) THEN
          I = NINT(FVALUE(N))
        ELSE IF (ITYP(N).EQ.1) THEN
           WRITE (LINERR,FMT='(A,I4)') ' Illegal number in field ',N
C
C              ****************************
          CALL LERROR(1,0,LINERR)
C              ****************************
C
          LFLAG = +1
        END IF
      ELSE
         LFLAG = -1
      END IF
      END
C
C_BEGIN_RDATOMSELECT
C     ==================================================================
      SUBROUTINE RDATOMSELECT(JTOK,INAT0,INAT1,IRES0,IRES1,CHNAM,
     +                        IMODE,NTOK,LINE,IBEG,IEND,ITYP,IDEC,
     +                        FVALUE,IFAIL)
C     ==================================================================
C
C  Subroutine to process atom selection keyword with the following
C  general syntax:
C
C    <Keywords...> ATOM <inat0> [ [TO] <inat1> ] |
C                  RESIDUE [ALL | ONS | CA] [ CHAIN <chnam> ]
C                  <ires0> [ [TO] <ires1> ]
C
C     e.g. kywd atom 1 to 100
C          kywd residue chain A 20 to 30
C          kywd residue all 11 32    etc...
C
C  To be compatible with DISTANG, CONTACT etc the ordering of the
C  RESIDUE subarguments is flexible, eg RESIDUE 1 TO 9 CA CHAIN B
C  is the same as RESIDUE CA CHAIN B 1 TO 9...
C
C  The subroutine returns the selection entered by the user and expects the
C  calling program to deal with the results. The preceeding keywords are
C  relevant for this subroutine
C
C  ARGUMENTS
C  =========
C
C     JTOK    (I) INTEGER       Number of first field to interpret
C     NTOK    (I) INTEGER       The number of fields parsed, from PARSER
C     LINE    (I) CHARACTER*(*) Input string, from PARSER
C     IBEG    (I) INTEGER(*)    1st column number of tokens in field
C                               (from PARSER)
C     IEND    (I) INTEGER(*)    Last column number of tokens in field
C                               (from PARSER)
C     ITYP    (I) INTEGER(*)    =0  null field
C                               =1  character string
C                               =2  number   (from PARSER)
C     IDEC    (I) INTEGER(*)    Number of characters/digits in each token
C                               (from PARSER)
C     FVALUE  (I) REAL(*)       Array of numbers. (from PARSER)
C
C     INAT0   (O) INTEGER       Lower limit of atom range (-99 if not set)
C     INAT1   (O) INTEGER       Upper limit of atom range (-99 if not set)
C     IRES0   (O) INTEGER       Lower limit of residue range (-99 if not set)
C     IRES1   (O) INTEGER       Upper limit of residue range (-99 if not set)
C     CHNAM   (O) CHARACTER*(*) Chain identifier (' ' if not set)
C     IMODE (I/O) INTEGER       On entry: -1 = don't allow MODE
C                                         any other value = allow MODE
C                               On exit:  Type of atoms to include:
C                                          1=ALL   2=ONS   3=CA (see eg CONTACT)
C     IFAIL (I/O) INTEGER       On entry:  0 = suppress warnings
C                                         -1 = print warnings
C                               On exit:   0 = LINE parsed ok
C                                         >0 = error occured parsing line
C                                              (value of IFAIL is no. of bad token)
C
C  RETURNED VALUES
C  ===============
C
C  The subroutine returns either:
C
C  1. first/last atom numbers, defining a range of atoms, or
C  2. first/last residue numbers, defining a range of residues, plus
C        (optionally) chain identifier
C        (optionally) a MODE which specifies which type of atoms to
C        include: all = (default) all atoms in residue range
C                 ons = only oxygen and nitrogen atoms
C                 ca  = only CA atoms
C        (see CONTACT/DISTANG)
C
C  Unset atoms/residue numbers will be returned < 0 (i.e. -99)
C  Unset chain identifier will be returned as a blank, i.e. ' '
C  Mode defaults to 1 = include all types of atoms.
C
C_END_RDATOMSELECT
C
      IMPLICIT NONE
C
C     ..Parameters
      INTEGER MAXTOK
      PARAMETER (MAXTOK=20)
C
C     ..Scalar arguments
      INTEGER   NTOK,JTOK,INAT0,INAT1,IRES0,IRES1,IMODE,IFAIL
      CHARACTER LINE*80,CHNAM*(*)
C
C     ..Array arguments
      INTEGER IBEG(MAXTOK),IEND(MAXTOK),ITYP(MAXTOK),IDEC(MAXTOK)
      REAL    FVALUE(MAXTOK)
C
C     ..Local scalars
      INTEGER ITOK,NLIMIT,TEMP
      CHARACTER KEY*4,ERRLINE*80
      LOGICAL LATOM,LRESI,LCHAIN,LLIMIT
C
C     ..Local arrays
      INTEGER ILIMIT(2)
C
C     ..External subroutines/functions
      EXTERNAL CCPUPC,CCPERR
C
C---- Initial checks
C
      IF (NTOK.GT.MAXTOK) THEN
        ERRLINE = 'RD_KEY_SELECT: too many arguments'
        GO TO 9010
      END IF
C
      IF (JTOK.LT.1 .OR. JTOK.GT.NTOK) THEN
        ERRLINE = 'RD_KEY_SELECT: JTOK out of range'
        GO TO 9010
      END IF
C
C---- Initialise values
C
      INAT0 = -99
      INAT1 = -99
      LATOM = .FALSE.
      IRES0 = -99
      IRES1 = -99
      CHNAM = ' '
C
C---- IMODE
C
      IF (IMODE.NE.-1) IMODE = 0
C
C---- Flags for modes
C
      LATOM  = .FALSE.
      LRESI  = .FALSE.
      LCHAIN = .FALSE.
      LLIMIT = .FALSE.
      NLIMIT = 0
      ERRLINE = ' '
C
C---- Step through line token at a time
C
      ITOK = JTOK
C
 9000 CONTINUE
      KEY  = LINE(IBEG(ITOK):IEND(ITOK))
      CALL CCPUPC(KEY)
C
C---- ATOM
C     ====
      IF (KEY(1:4).EQ.'ATOM') THEN
        IF (LATOM.OR.LRESI) ERRLINE = 'Multiple ATOM/RESidue keywords'
        LATOM = .TRUE.
C
C---- RESIDUE
C     =======
      ELSE IF (KEY(1:3).EQ.'RES') THEN
        IF (LATOM.OR.LRESI) ERRLINE = 'Multiple ATOM/RESidue keywords'
        LRESI = .TRUE.
C
C---- MODE: CA / ONS / ALL
C     ====================
      ELSE IF (KEY(1:3).EQ.'ALL'.OR.
     +         KEY(1:3).EQ.'ONS'.OR.
     +         KEY(1:2).EQ.'CA') THEN
C
        IF (IMODE.EQ.-1)
     +  ERRLINE = 'ALL/ONS/CA: invalid specifiers'
        IF (.NOT.LRESI)
     +  ERRLINE = 'ALL/ONS/CA not allowed without RESidue'
        IF (IMODE.GT.0) ERRLINE = 'Only one of CA/ONS/ALL allowed'
C
        IF (KEY(1:3).EQ.'ALL') IMODE = 1
        IF (KEY(1:3).EQ.'ONS') IMODE = 2
        IF (KEY(1:2).EQ.'CA')  IMODE = 3
C
C---- CHAIN <chnam>
C     =============
      ELSE IF (KEY(1:4).EQ.'CHAI') THEN
        IF (.NOT.LRESI) ERRLINE = 'CHAIN only allowed after RESidue'
        IF (LCHAIN) ERRLINE = 'Only one CHAIN allowed per line'
        ITOK = ITOK + 1
        IF (ITYP(ITOK).EQ.1 .AND. IDEC(ITOK).LE.2) THEN
          CHNAM = LINE(IBEG(ITOK):IEND(ITOK))
          LCHAIN = .TRUE.
        ELSE
          ERRLINE = 'Chain name should be one or two characters'
        END IF
C
C---- Number ...
C     ==========
C
C     These are atom or residue limits ... Process them all together
C     The possibilities are: 1 number, 2 numbers or 2 numbers separated
C     by "TO"
      ELSE IF (ITYP(ITOK).EQ.2) THEN
        IF (.NOT.(LATOM .OR. LRESI))
     +  ERRLINE = ' Missing ATOM or RESIDUE keyword'
        IF (LLIMIT)
     +  ERRLINE = ' Already read a set of atom/residue limits'
        ILIMIT(1) = INT(FVALUE(ITOK))
        ITOK = ITOK + 1
        NLIMIT = NLIMIT + 1
C
C     Check the next argument - is it "TO"?
        IF (ITYP(ITOK).EQ.1) THEN
          KEY = LINE(IBEG(ITOK):IEND(ITOK))
          CALL CCPUPC(KEY)
          IF (KEY(1:2).EQ.'TO') THEN
            IF(ITYP(ITOK+1).NE.2 .OR. ITOK.EQ.NTOK)
     +      ERRLINE = 'TO should be followed by a number'
            ITOK = ITOK + 1
          END IF
        END IF
C
C     Check if the next argument is the second limit
        IF (ITYP(ITOK).EQ.2) THEN
          ILIMIT(2) = INT(FVALUE(ITOK))
          NLIMIT = NLIMIT + 1
        ELSE
C
C     Otherwise, need to step back one token for next round
          ITOK = ITOK - 1
        END IF
        LLIMIT = .TRUE.
C
C---- Check for TO out of place
C
      ELSE IF (KEY(1:2).EQ.'TO') THEN
        ERRLINE = 'TO out of place'
C
C---- Keyword unrecognised otherwise
C
      ELSE
        write(6,*)key
        ERRLINE = 'Unrecognised subkeyword'
      END IF
C
C---- Use ERRLINE to trap errors
C
      IF (ERRLINE.NE.' ') GO TO 9010
C
C---- Parse some more?
C
      ITOK = ITOK + 1
      IF (ITOK.LE.NTOK) GO TO 9000
C
C---- Sort out the input
C
      IF (NLIMIT.EQ.0) THEN
        ERRLINE = 'Atom/residue limits unspecified'
        GO TO 9010
      ELSE IF (NLIMIT.EQ.1) THEN
        ILIMIT(2) = ILIMIT(1)
      ELSE IF (ILIMIT(1).GT.ILIMIT(2)) THEN
        TEMP = ILIMIT(1)
        ILIMIT(1) = ILIMIT(2)
        ILIMIT(2) = TEMP
      END IF
C
      IF (LATOM) THEN
        INAT0 = ILIMIT(1)
        INAT1 = ILIMIT(2)
      ELSE
        IRES0 = ILIMIT(1)
        IRES1 = ILIMIT(2)
      END IF
C
      IF (IMODE.EQ.0) IMODE = 1
C
      IFAIL = 0
      RETURN
C
C---- Errors come to here
C
 9010 IF (IFAIL.LT.0) CALL CCPERR(2, ERRLINE)
      IFAIL = ITOK
      RETURN
C
      END
