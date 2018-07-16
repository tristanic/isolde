C
C     ccplib.f: Pseudo-machine-independent low-level routines
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
C_BEGIN_CCPLIB
C     These are supposedly-machine-independent low-level routines.
C     They're actually machine-dependent at least insofar as some
C     contain non-standard code, but they do compile with the compilers
C     tried on unix as well as VMS.
C
C     fixme: the bit-twiddling should be in library.c, not here.
C     amalgamate ccppsf and fdir/fext/froot.  also add tests of these
C     routines to testlib.
C
C     $Id$
C
C      CCFILL    Set specified number of elements of byte array
C      CCPALC    Call subroutine with allocated memory
C      CCPALE    Call subroutine with allocated memory set from environment
C      CCPBYI    Copy array of unsigned (or signed) bytes into integer array
C      CCPBYT    Indicate whether byte handling is available
C      CCPCPI    Copy array of BYTE or INTEGER*2 elements into integer array
C      CCPDEX    Periodicity reduction of 1024 (for PROLSQ)
C      CCPDPN    more friendly CCPOPN
C      CCPE2I    read integer from logical name value
C      CCPGI2    Get unsigned integer*2 value from 0 to 65535 from N'th
C                unsigned integer*2 element of array.
C      CCPGTB    Get unsigned byte value from 0 to 255 from N'th byte of
C                array.
C      CCPI2I    Copy an array of INTEGER*2 elements into an integer array
C      CCPIBY    Copy array of integers into array of bytes.
C      CCPII2    Copy array of integers into array of INTEGER*2 elements.
C      CCPMDE    If byte handling available return nos. of bytes for map
C                modes
C      CCPMVB    Move bytes from one non-character array to another if
C                byte handling is available
C      CCPMVI    Move words from one integer array to another
C      CCPMVR    Move words from one real array to another
C      CCPNUN    Return an unconnected i/o unit number
C      CCPONL    See if program is being run interactively
C      CCPPSF    Parse file name into components
C      CSETNV    Associate logical name with file name
C      CCPPAG    Set paging parameters if available
C      CCPSI2    Set integer value from 0 to 65535 into the N'th
C                unsigned integer*2 element of an array.
C      CCPSTB    Set integer value from 0 to 255 into N'th byte of array.
C      CCPSUM    Sum the elements of an array
C      CCPTOI    Convert n'th byte or I*2 in a non-character array to an
C                integer
C      CCPUFL    Supress underflow messages
C      CCPZBI    Sets an array of bytes to zero
C      CCPZI     Set 'n' words of an integer array to zero using a simple loop
C      CCPZR     Set 'n' words of a real array to zero using a simple loop
C      FDIR      Returns the directory part of a file name
C      FEXTN     Returns the extension of a file name
C      FROOT     Returns the root of a file name
C      LITEND    determine endianness
C      LENSTR    length of string to last non-space ( C equiv 
C                ccp4_utils_flength
C      LUNSTI    Get logical unit number for input
C      LUNSTO    Get logical unit number for output
C      NBITST    Return the (unsigned) integer value held within a bit
C                field in a word
C      NOCRLF    write line supressing cr/lf to standard output
C      STBITS    Set a bit field within a word to a given (unsigned)
C                integer value
C      HKLEQ     Are the reflection indices the equal
C      Subroutines for generating and accessing a hash table
C      CCP4_HASH_SETUP
C      CCP4_HASH_LOOKUP
C      CCP4_HASH_ZEROIT
C
C_END_CCPLIB
C
C
C
C_BEGIN_CCFILL
      SUBROUTINE CCFILL(ARR1,SCAL,NTIMES)
C     ===================================
C
C      CCFILL    Set NTIMES bytes array ARR1 to value SCAL
C
C Arguments:
C ==========
C
C        ARR1 (O)   BYTE ARRAY (*): WHERE BYTES ARE TO BE COPIED
C        SCAL (I)   BYTE: value to be copied into ARR1
C      NTIMES (I)   INTEGER: NUMBER OF BYTES TO BE COPIED
C_END_CCFILL
C
C     .. Scalar Arguments ..
      INTEGER NTIMES
      INTEGER*1 SCAL
C     ..
C     .. Array Arguments ..
      INTEGER*1 ARR1(*)
C     ..
C     .. Local Scalars ..
      INTEGER N
C     ..
      DO 10 N = 1,NTIMES
        ARR1(N) = SCAL
   10 CONTINUE
C
      END
C
C
C_BEGIN_CCPALC
      SUBROUTINE CCPALC(ROUTNE, N, TYPE, LENGTH)
C     ==========================================
C
C     Arrange to call subroutine ROUTNE with N array arguments each of
C     length LENGTH (i) and type indicated by TYPE (i): 'i' == integer,
C     'r' == real, 'd' == double precision, 'c' == complex, 'b' ==
C     "byte" (logical*1 or integer*1, unportable and deprecated) .  TYPE
C     elements may have either case.
C     Consider `call ccpalc (fred, 3, types, lens)' with types = (/'i',
C     'r', 'c'/)  and lens = (/1000, 2000, 3000/).  This effectively does
C        call fred (1000, arr1, 2000, arr2, 3000, arr3)
C     with
C        subroutine fred (n1, foo, n2, bar, n3, baz)
C        integer n1, n2, n3, foo (n1)
C        real bar (n2)
C        complex baz (n3)
C        ...
C     Obviously all communication with ROUTNE must be by COMMON (or,
C     possibly, extra ENTRYs).  The allocated memory is freed on return
C     from ROUTNE.  As a concession, it's initially filled with zeroed
C     bytes.

C
C Arguments:
C ==========
C
C      ROUTNE (I)   EXTERNAL: routine to call
C           N (I)   INTEGER: number of arguments to ROUTNE (<=12)
C        TYPE (I)   CHARACTER*1 (*): type of arguments to ROUTNE:
C                      'I': INTEGER; 'R': REAL; 'D': DOUBLE PRECISION;
C                      'C': COMPLEX; 'B': LOGICAL*1 or INTEGER*1
C      LENGTH (I)   INTEGER*(*): number of elements in each (array)
C                       argument of ROUTNE
C_END_CCPALC
C
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      CHARACTER TYPE (*)
      INTEGER LENGTH (*)
C     ..
      EXTERNAL ROUTNE, CCPAL1, CCPUPC
      INTEGER I, ITYPE (12)
      CHARACTER TTYPE (12)
C     ..
      IF (N.LT.1 .OR. N.GT.12)
     +     CALL CCPERR (1, 'CCPALC: bad number of arguments')
      DO 10 I=1,N
        TTYPE (I) = TYPE (I)
        CALL CCPUPC (TTYPE (I))
        ITYPE (I) = INDEX ('IRDCB', TTYPE (I))
        IF (ITYPE (I) .EQ. 0) CALL CCPERR (1, 'CCPALC: bad TYPE: '//
     +       TYPE (I))
        IF (LENGTH (I).LE.0) CALL CCPERR (1, 'CCPALC: length <=0')
 10   CONTINUE
      CALL CCPAL1 (ROUTNE, N, ITYPE, LENGTH)
      END
C
C
C_BEGIN_CCPALE
      SUBROUTINE CCPALE(ROUTNE, N, TYPE, LENGTH, LENDEF, PRINT)
C     =================================================
C
C     Arrange to call subroutine ROUTNE with N array arguments each of
C     length LENGTH (i) and type indicated by TYPE (i): 'i' == integer,
C     'r' == real, 'd' == double precision, 'c' == complex, 'b' == byte.
C     TYPE elements may have either case.  LENGTH points to an array of
C     environment variable (logical) names from which integer values are
C     read.  The lengths default to values from LENDEF.
C     This is a convenient interface to CCPALC to allow configuring of
C     the memory requirements on the command line where appropriate.
C     This may be useful if the memory requirements can't be determined
C     initially and it's necessary to guess.
C
C Arguments:
C ==========
C
C      ROUTNE (I)   EXTERNAL: routine to call
C           N (I)   INTEGER: number of arguments to ROUTNE (<=12)
C        TYPE (I)   CHARACTER*1 (*): type of arguments to ROUTNE:
C                      'I': INTEGER; 'R': REAL; 'D': DOUBLE PRECISION;
C                      'C': COMPLEX; 'B': LOGICAL*1 or INTEGER*1
C     LENGTH (I)   CHARACTER *(*): logical names representing the number
C                       of elements in each (array) argument of ROUTNE
C     LENDEF (I)   INTEGER (*): default lengths for the argument arrays
C     used if the appropriate LENGTH argument doesn't represent a
C     defined logical
C     PRINT  (I)   LOGICAL: whether or not to print the values of the
C     array lengths
C_END_CCPALE
C
C     .. Scalar Arguments ..
      INTEGER N
      LOGICAL PRINT
C     ..
C     .. Array Arguments ..
      CHARACTER TYPE (*),  LENGTH (*)*(*)
      INTEGER LENDEF (*)
C     ..
      EXTERNAL ROUTNE, CCPE2I, CCPALC, LUNSTO
      INTEGER I, LENG (12), CCPE2I, LUNSTO
C     ..
      DO 10 I=1,N
        LENG (I) = CCPE2I (LENGTH (I), LENDEF (I))
 10   CONTINUE
      IF (PRINT) THEN
        WRITE (LUNSTO(1), 
     +     '(/'' Memory allocation (logical name, type, elements):'')')
        WRITE (LUNSTO(1), '(3X, A, 1X, A, 3X, I10)')
     +       (LENGTH (I), TYPE (I), LENG (I), I=1,N)
      ENDIF
      CALL CCPALC (ROUTNE, N, TYPE, LENG)
      END
C
C
C
C SUBROUTINE 'CCPBYI'
C ===================
C
C_BEGIN_CCPBYI
      SUBROUTINE CCPBYI(IA,IBYT,NB)
C     =============================
C
C COPY AN ARRAY OF UNSIGNED (OR SIGNED) BYTES INTO AN INTEGER ARRAY
C
C (MUST BE IMPLEMENTED IF CCPBYT FUNCTION RETURNS .TRUE.)
C [added for LAUE]
C
C Arguments:
C ==========
C
C      IA (O)   INTEGER ARRAY(*): TO RETURN INTEGER VALUES
C    IBYT (I)   BYTE ARRAY(*): DATA (MAY BE AN INTEGER ARRAY FOR EXAMPLE
C               WITH DATA PACKED INTO ADJACENT BYTES
C      NB (I)   INTEGER: IF >0, THE NUMBER OF UNSIGNED BYTES TO BE COPIED
C                        IF <0, -THE NUMBER OF SIGNED BYTES TO BE COPIED
C_END_CCPBYI
C
C SPECIFICATION STATEMENTS
C ------------------------
C
      INTEGER IA(*)
      INTEGER*1 IBYT(*)
      INTEGER*1 JBYT(4)
      EQUIVALENCE (JA,JBYT(1))
      LOGICAL CALLED, LITEND
      INTEGER IND
      EXTERNAL LITEND
      SAVE CALLED, IND
      DATA CALLED/.FALSE./
C
      IF (.NOT.CALLED) THEN
        CALLED=.TRUE.
        IF (LITEND(1)) THEN
          IND = 1
        ELSE
          IND = 4
        ENDIF
      ENDIF
C
C COPY DATA
C ---------
C
      NE = NB
      IF (NE.GT.0) THEN
         JA=0
         DO 10 I=1,NE
           JBYT(IND)=IBYT(I)
           IA(I)=JA
 10      CONTINUE
      ELSE
         NE = -NE
         DO 20 I=1,NE
         IA(I) = IBYT(I)
20       CONTINUE
      END IF
      END
C
C
C
C_BEGIN_CCPBYT
      LOGICAL FUNCTION CCPBYT(NBW)
C     ============================
C
C---- This function indicates whether byte handling is available or not.
C      if a value of .true. is returned then the subroutines ccpmde and
C      ccpmvb must be fully implemented.
C
C Arguments:
C ==========
C
C         NBW (O)   INTEGER: RETURNS THE NUMBER OF BYTES PER WORD OR A VALUE
C                   OF 1 IF NO BYTE HANDLING IS AVAILABLE.
C
C  RETURNS   CCPBYT  = .TRUE.  BYTE HANDLING AND ASSOCIATED CCPLIB
C                              ROUTINES AVAILABLE.
C                    = .FALSE. NO BYTE HANDLING AVAILABLE.
C_END_CCPBYT
C
C     .. Scalar Arguments ..
      INTEGER NBW
C     ..
      CCPBYT = .TRUE.
      NBW = 4
      END
C
C
C SUBROUTINE 'CCPCPI'
C ===================
C_BEGIN_CCPCPI
      SUBROUTINE CCPCPI(IA,IB,MINEL,MAXEL,ITYP)
C     =========================================
C
C Copy an array of BYTE or INTEGER*2 elements into an integer array
C
C (Must be implemented if ccpbyt function returns .TRUE.)
C [for LAUE]
C
C Arguments:
C ==========
C
C      IA (O)   INTEGER Array(*): to return values
C      IB (I)   INTEGER Array(*): holding data with data packed into adjacant 
C                                 BYTE or INTEGER*2 elements
C   MINEL (I)   INTEGER: Minimum element to copy
C   MAXEL (I)   INTEGER: Maximum element to copy
C    ITYP (I)   INTEGER: Type =1 unsigned byte
C                             =2 signed byte
C                             =3 unsigned two byte integer
C                             =4 signed two byte integer
C
C               Note: if MINEL>MAXEL elements will be copied in reverse order
C_END_CCPCPI
C
C====== Specification statements
C
      INTEGER IA(*)
      INTEGER*1 IB(*)
      INTEGER*2 J2(2)
      INTEGER*1 JBYT(4)
      EQUIVALENCE (JA,J2(1),JBYT(1))
      LOGICAL CALLED, LITEND
      EXTERNAL LITEND
      INTEGER IND1, IND2, INDB
      SAVE CALLED, IND1, IND2, INDB
      DATA CALLED/.FALSE./
C
      IF (.NOT.CALLED) THEN
        IF (LITEND(1)) THEN
          IND1 = 1
          IND2 = 2
          INDB = 1
        ELSE
          IND1 = 3
          IND2 = 4
          INDB = 4
        ENDIF
        CALLED=.TRUE.
      ENDIF
C
C====== Copy data
C
      ISTEP = 1
      IF (MINEL.GT.MAXEL) ISTEP=-1
      IF (ITYP.EQ.1) THEN
         JA=0
         J=0
         DO 10 I=MINEL,MAXEL,ISTEP
            J=J+1
            JBYT(INDB)=IB(I)
            IA(J)=JA
10       CONTINUE
      ELSE IF (ITYP.EQ.2) THEN
         J=0
         DO 20 I=MINEL,MAXEL,ISTEP
            J=J+1
            IA(J)=IB(I)
20       CONTINUE
      ELSE IF (ITYP.EQ.3) THEN
         JA=0
         J=0
         DO 30 I=MINEL,MAXEL,ISTEP
            J=J+1
            JBYT(IND1)=IB(2*I-1)
            JBYT(IND2)=IB(2*I)
            IA(J)=JA
30       CONTINUE
      ELSE IF (ITYP.EQ.4) THEN
         J=0
         DO 40 I=MINEL,MAXEL,ISTEP
            J=J+1
            JBYT(1)=IB(2*I-1)
            JBYT(2)=IB(2*I)
            IA(J)=J2(1)
40       CONTINUE
      END IF
      RETURN
      END
C
C
C
C_BEGIN_CCPDEX
      SUBROUTINE CCPDEX(INDX,N)
C     =========================
C
C---- This subroutine performs a periodicity reduction for a period
C     of 1024 for the elements of an array. written particularly for
C     'prolsq' to allow for use of the 'and' function on the cray or
C     'moveb' on the m28oh(iap).
C      These are much faster than the mod function used in
C      the standard fortran77 version.
C
C Arguments:
C ==========
C
C        INDX (I/O) INTEGER ARRAY(*): NUMBERS FOR PERIODICITY REDUCTION
C           N (I)   INTEGER: NO. OF ELEMENTS IN INDX
C
C EXAMPLE OF FUNCTIONS:
C
C FORTRAN77     INDX(I)=MOD(INDX(I),1024)+1
C CRAY-1S       INDX(I)=AND(INDX(I),1023)+1
C M280H(IAP)    CALL MOVEB(INDX(I),1,0,1,22)
C               INDX(I)=INDX(I)+1
C_END_CCPDEX
C
C SPECIFICATION STATEMENTS AND CODE
C
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER INDX(N)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
      DO 10 I = 1,N
        INDX(I) = MOD(INDX(I),1024) + 1
   10 CONTINUE
      END
C
C_BEGIN_CCPDPN
      SUBROUTINE CCPDPN(IUN,LOGNAM,STATUS,TYPE,LREC,IFAIL)
C     ====================================================
C
C---- Calls CCPOPN to open a file, but with mnemonic arguments
C
C Arguments:
C ==========
C
C         IUN (I)   INTEGER: UNIT NUMBER
C
C      LOGNAM (I)   CHARACTER*(*): LOGICAL FILE NAME
C
C      STATUS (I)   CHARACTER*(*): FILE STATUS FLAG:
C                                     'UNKNOWN'
C                                     'SCRATCH'
C                                     'OLD'
C                                     'NEW'
C                                     'READONLY'
C                                     'PRINTER'
C
C        TYPE (I)   CHARACTER*(*): FILE TYPE FLAG:
C                                  ='F', 'SEQUENTIAL' 'FORMATTED'
C                                  ='U', 'SEQUENTIAL' 'UNFORMATTED'
C                                  ='DF', 'DIRECT'     'FORMATTED'
C                                  ='DU', 'DIRECT'     'UNFORMATTED'
C     [STATUS and TYPE are case-insensitive]
C
C        LREC (I)   INTEGER: RECORD LENGTH FOR DIRECT ACCESS FILE (NO. OF
C                   CHARACTERS FOR A FORMATTED FILE OR WORDS FOR
C                   AN UNFORMATTED FILE). NOT RELEVANT FOR A SEQUENTIAL
C                   FILE
C
C       IFAIL (I/O) INTEGER: ON INPUT     =0, STOP ON OPEN FAILURE
C                                         =1, CONTINUE AFTER OPEN FAILURE
C                                             (only on file not found)
C                                         =2, CONTINUE SILENTLY AFTER OPEN FAILURE
C                                         =-1, As 0, but silent on success
C                                             (equivalent to negative IUN)
C                            ON OUTPUT    UNCHANGED IF FILE OPEN OK
C                                         =-1, ERROR IN OPENING FILE
C_END_CCPDPN
C
C     .. Scalar Arguments ..
      INTEGER IFAIL,IUN,IUN1,LREC
      CHARACTER LOGNAM* (*),STATUS* (*),TYPE* (*)
C     ..
C     .. Local Scalars ..
      INTEGER ISTAT,ITYPE
      CHARACTER ERRSTR*80
C     ..
C     .. Local Arrays ..
      CHARACTER TYPES(4)*2,STATS(6)*8, STAT*8, TYP*2
C     ..
C     .. External Functions ..
      INTEGER CCPNUN,LENSTR
      EXTERNAL CCPNUN,LENSTR
C     ..
C     .. External Subroutines ..
      EXTERNAL CCPOPN
C     ..
C     .. Data statements ..
      DATA STATS/'UNKNOWN','SCRATCH','OLD','NEW','READONLY','PRINTER'/
      DATA TYPES/'F','U','DF','DU'/
C     ..
C
      IF (IUN .EQ. 0) IUN = CCPNUN()
      STAT = STATUS
      TYP = TYPE
      CALL CCPUPC(STAT)
      CALL CCPUPC(TYP)
      DO 10 ISTAT = 1,6
        IF (STAT.EQ.STATS(ISTAT)) GO TO 20
   10 CONTINUE
      ERRSTR = ' CCPDPN: illegal status : '
      ERRSTR(LENSTR(ERRSTR)+2:) = STATUS
      CALL CCPERR(1,ERRSTR)
C
   20 DO 30 ITYPE = 1,4
        IF (TYP.EQ.TYPES(ITYPE)) GO TO 40
   30 CONTINUE
      ERRSTR = ' CCPDPN: illegal type: '
      ERRSTR(LENSTR(ERRSTR)+2:) = TYPE
      CALL CCPERR(1,ERRSTR)
C
 40   CONTINUE
      IUN1 = IUN
C  If IFAIL lt 0 No open message from CCPOPN
      IF(IFAIL.LT.0 .AND. IUN.GT.0) THEN
        IUN1 = -IUN
        IFAIL = 0
      ENDIF
      CALL CCPOPN(IUN1,LOGNAM,ISTAT,ITYPE,LREC,IFAIL)
C
      END
C
C
C_BEGIN_CCPE2I
      INTEGER FUNCTION CCPE2I (NAME, DEFVAL)
C     ======================================
C
C     Return an integer extracted from enviroment variable NAME.  If
C     NAME isn't defined, use DEFVAL as the default.  If the value of
C     NAME isn't a representation of an integer, abort.
C
C     Arguments
C     =========
C
C     NAME (I)    CHARACTER *(*)
C     DEFVAL (I)  INTEGER
C_END_CCPE2I
      CHARACTER *(*) NAME
      CHARACTER BUFFER*80, EMESS*100
      INTEGER DEFVAL, LENSTR
      EXTERNAL UGTENV, LENSTR
      BUFFER = ' '
      CALL UGTENV (NAME, BUFFER)
      IF (BUFFER.EQ.' ') THEN
        CCPE2I = DEFVAL
        RETURN
      ENDIF
      READ (BUFFER, '(BN,I80)', ERR=99) CCPE2I
      RETURN 
 99   EMESS = ' Logical name '
      EMESS(LENSTR(EMESS)+2:) = NAME(1:LENSTR(NAME))
      IF(LENSTR(EMESS) .LE. 99) THEN
        EMESS(LENSTR(EMESS)+1:) =' should represent an integer and is: '
        IF(LENSTR(EMESS) .LE. 98) 
     .           EMESS(LENSTR(EMESS)+2:) = BUFFER(1:LENSTR(BUFFER))
      ENDIF
      CALL CCPERR (1, EMESS)
      END
C
C
C SUBROUTINE 'CCPGI2'
C ===================
C
C_BEGIN_CCPGI2
      SUBROUTINE CCPGI2(IVAL,IA,N)
C     ============================
C
C GET AN UNSIGNED INTEGER*2 VALUE FROM 0 TO 65535 FROM THE N'TH unsigned
C INTEGER*2 ELEMENT OF AN INTEGER (OR OTHER) ARRAY.
C
C (MUST BE IMPLEMENTED IF CCPBYT FUNCTION RETURNS .TRUE.)
C [added for LAUE]
C
C Arguments:
C ==========
C
C    IVAL (O)   INTEGER: THE RETURNED VALUE FROM 0 TO 65535
C      IA (I/O) INTEGER*2 ARRAY(*): FROM WHICH THE UNSIGNED INTEGER*2 VALUE
C               IS TO BE RETRIEVED
C       N (I)   INTEGER: POSITION IN 'IA' WHERE THE UNSIGNED INTEGER*2 VALUE
C               IS TO BE RETRIEVED
C_END_CCPGI2
C
C SPECIFICATION STATEMENTS
C ------------------------
C
      INTEGER*2 IA(*)
      INTEGER*2 JBYT(2)
      EQUIVALENCE (JA,JBYT(1))
      LOGICAL CALLED, LITEND
      EXTERNAL LITEND
      INTEGER IND
      SAVE CALLED, IND
      DATA CALLED/.FALSE./
C
      IF (.NOT.CALLED) THEN
        IF (LITEND(1)) THEN
          IND = 1
        ELSE
          IND = 2
        ENDIF
        CALLED=.TRUE.
      ENDIF
C
C GET UNSIGNED INTEGER*2
C ----------------------
C
      JA=0
      JBYT(IND)=IA(N)
      IVAL=JA
      END
C
C
C SUBROUTINE 'CCPGTB'
C ===================
C
C_BEGIN_CCPGTB
      SUBROUTINE CCPGTB(IVAL,IA,N)
C     ============================
C
C GET AN UNSIGNED BYTE VALUE FROM 0 TO 255 FROM THE N'TH BYTE OF AN INTEGER
C (OR OTHER) ARRAY.
C
C (MUST BE IMPLEMENTED IF CCPBYT FUNCTION RETURNS .TRUE.)
C [for LAUE]
C
C Arguments:
C ==========
C
C    IVAL (O)   INTEGER: THE RETURNED VALUE FROM 0 TO 255
C      IA (I/O) BYTE ARRAY(*): FROM WHICH THE BYTE VALUE IS TO BE RETRIEVED
C       N (I)   INTEGER: THE POSITION IN 'IA' WHERE THE BYTE VALUE IS
C               TO BE RETRIEVED
C_END_CCPGTB
C
C SPECIFICATION STATEMENTS
C ------------------------
C
      INTEGER*1 IA(*)
      INTEGER*1 JBYT(4)
      EQUIVALENCE (JA,JBYT(1))
      LOGICAL CALLED, LITEND
      EXTERNAL LITEND
      INTEGER IND
      SAVE CALLED, IND
      DATA CALLED/.FALSE./
C
      IF (.NOT.CALLED) THEN
        IF (LITEND(1)) THEN
          IND = 1
        ELSE
          IND = 4
        ENDIF
        CALLED=.TRUE.
      ENDIF
C
C GET BYTE
C --------
C
      JA=0
      JBYT(IND)=IA(N)
      IVAL=JA
      END
C
C
C SUBROUTINE 'CCPI2I'
C ===================
C_BEGIN_CCPI2I
      SUBROUTINE CCPI2I(IA,I2,NE,SIGNED,SWAPB)
C     ========================================
C
C Copy an array of INTEGER*2 elements into an integer array
C
C (Must be implemented if ccpbyt function returns .TRUE.)
C [for LAUE]
C
C Arguments:
C ==========
C
C      IA (O)   INTEGER Array(*): to return values
C      I2 (I)   INTEGER*2 Array(*): holding data (may be an INTEGER array for
C               example with data packed into adjacant INTEGER*2 elements
C      NE (I)   INTEGER: The number of elements to be copied
C  SIGNED (I)   LOGICAL: =.TRUE.  Copy as signed integer*2 values
C                        =.FALSE. Copy as unsigned integer*2 values
C   SWAPB (I)   LOGICAL: =.TRUE.  Swap bytes in the integer*2 elements
C                        =.FALSE. Do not swap bytes
C_END_CCPI2I
C
C====== Specification statements
C
      LOGICAL SIGNED, SWAPB
      INTEGER IA(*),JA
      INTEGER*2 I2(*)
      INTEGER*2 J2(2)
      INTEGER*2 IEIGHT, I255
      PARAMETER (I255=255)
      EQUIVALENCE (JA,J2(1))
      LOGICAL CALLED, LITEND
      EXTERNAL LITEND
      INTEGER IND,I,NE
      SAVE CALLED, IND
      DATA CALLED/.FALSE./
C
      IF (.NOT.CALLED) THEN
        IF (LITEND(1)) THEN
          IND = 1
        ELSE
          IND = 2
        ENDIF
        CALLED=.TRUE.
      ENDIF
C
C====== Swap bytes if required
C
      IEIGHT = 8
      IF (SWAPB) THEN
         DO 10 I = 1,NE
            I2(I) = IOR(IAND(ISHFT(I2(I),-IEIGHT),I255),
     +              ISHFT(I2(I),IEIGHT))
10       CONTINUE
      END IF
C
C====== Copy data
C
      IF (SIGNED) THEN
         DO 20 I=1,NE
            IA(I) = I2(I)
20       CONTINUE
      ELSE
         JA=0
         DO 30 I=1,NE
         J2(IND)=I2(I)
         IA(I)=JA
30       CONTINUE
      END IF
      END
C
C
C
C SUBROUTINE 'CCPIBY'
C ===================
C
C_BEGIN_CCPIBY
      SUBROUTINE CCPIBY(IBYT,IA,NB)
C     =============================
C
C COPY AN ARRAY OF INTEGERS INTO AN ARRAY OF UNSIGNED (OR UNSIGNED) BYTES.
C NOTE: NO OVERFLOW CHECKING IS DONE.
C
C (MUST BE IMPLEMENTED IF CCPBYT FUNCTION RETURNS .TRUE.)
C [for LAUE]
C
C Arguments:
C ==========
C
C    IBYT (O)   BYTE ARRAY(*) RETURNING DATA (MAY BE AN INTEGER ARRAY
C               FOR EXAMPLE WITH DATA PACKED INTO ADJACENT BYTES)
C      IA (I)   INTEGER ARRAY(*): Input values
C      NB (I)   INTEGER: IF >0, THE NUMBER OF ELEMENTS TO BE COPIED TO
C                        UNSIGNED BYTES
C                        IF <0, -THE NUMBER OF ELEMENTS TO BE COPIED TO
C                        SIGNED BYTES
C_END_CCPIBY
C
C SPECIFICATION STATEMENTS
C ------------------------
C
      INTEGER IA(*)
      INTEGER*1 IBYT(*)
      INTEGER*1 JBYT(4)
      EQUIVALENCE (JA,JBYT(1))
      LOGICAL CALLED, LITEND
      EXTERNAL LITEND
      INTEGER IND
      SAVE CALLED, IND
      DATA CALLED/.FALSE./
C
      IF (.NOT.CALLED) THEN
        IF (LITEND(1)) THEN
          IND = 1
        ELSE
          IND = 4
        ENDIF
        CALLED=.TRUE.
      ENDIF
C
C COPY DATA
C ---------
C
      NE = NB
      IF (NE.GT.0) THEN
         DO 10 I=1,NE
         JA=IA(I)
         IBYT(I)=JBYT(IND)
10       CONTINUE
      ELSE
         NE = -NE
         DO 20 I=1,NE
         IBYT(I) = IA(I)
20       CONTINUE
      END IF
      END
C
C
C
C SUBROUTINE 'CCPII2'
C ===================
C
C_BEGIN_CCPII2
      SUBROUTINE CCPII2(I2,IA,NE,SIGNED,SWAPB)
C     ========================================
C
C Copy an array of integers into an array of INTEGER*2 elements.
C NOTE: No overflow checking is done.
C
C (Must be implemented if ccpbyt function returns .TRUE.)
C [for LAUE]
C
C Arguments:
C ==========
C
C      I2 (O)   INTEGER*2 ARRAY(*): returning data (may be an INTEGER array for
C               example with data packed into adjacent INTEGER*2 elements)
C      IA (I)   INTEGER ARRAY(*): holding input values
C      NE (I)   INTEGER: The number of elements to be copied
C  SIGNED (I)   LOGICAL: =.TRUE.  Copy as signed integer*2 values
C                        =.FALSE. Copy as unsigned integer*2 values
C   SWAPB (I)   LOGICAL: =.TRUE.  Swap bytes in the integer*2 elements
C                        =.FALSE. Do not swap bytes
C_END_CCPII2
C
C====== Specification statements
C
      LOGICAL SIGNED, SWAPB
      INTEGER IA(*)
      INTEGER*2 I2(*)
      INTEGER*2 J2(2)
      INTEGER*2 IEIGHT, I255
      PARAMETER (I255=255)
      EQUIVALENCE (JA,J2(1))
      LOGICAL CALLED, LITEND
      EXTERNAL LITEND
      INTEGER IND
      SAVE CALLED, IND
      DATA CALLED/.FALSE./
C
      IF (.NOT.CALLED) THEN
        IF (LITEND(1)) THEN
          IND = 1
        ELSE
          IND = 2
        ENDIF
        CALLED=.TRUE.
      ENDIF
C
C====== Copy data
C
      IEIGHT = 8
      IF (SIGNED) THEN
         DO 10 I=1,NE
            I2(I) = IA(I)
10       CONTINUE
      ELSE
         DO 20 I=1,NE
            JA=IA(I)
            I2(I)=J2(IND)
20       CONTINUE
      ENDIF
C
C====== Swap bytes if required
C
      IF (SWAPB) THEN
         DO 30 I = 1,NE
            I2(I) = IOR(IAND(ISHFT(I2(I),-IEIGHT),I255),
     +              ISHFT(I2(I),IEIGHT))
30       CONTINUE
      END IF
      END
C
C
C
C_BEGIN_CCPMDE
      SUBROUTINE CCPMDE(MODE,NBYT)
C     ============================
C
C---- If byte handling is available (see ccpbyt) then this subroutine
C     returns the number of bytes per data item for the different modes
C     used, in particular, in the map handling subroutines.
C
C---- If byte handling is not available, then the number of words per
C     item is returned with zeros for the undefined items
C
C Arguments:
C ==========
C
C        MODE (I)   INTEGER:
C                   MODE = 0,   BYTES
C                        = 1,   SHORT (2 BYTE) INTEGERS
C                        = 2,   REAL/INTEGER (SINGLE WORD)
C                        = 3,   SHORT COMPLEX (2 * 2 BYTE INTEGERS)
C                        = 4,   COMPLEX (TWO WORDS)
C
C        NBYT (O)   INTEGER:
C                        > 0,   THE NUMBER OF BYTES FOR THE ITEM  IF
C                               CCPBYT RETURNS .TRUE. OR THE  NUMBER
C                               OF WORDS IF CCPBYT RETURNS .FALSE.
C                        = 0,   NO VALUE AVAILABLE FOR THIS MODE
C                        = -1,  INVALID MODE
C
C  TYPICAL VALUES:  1  2  4  4  8    IF BYTE HANDLING AVAILABLE WITH 4
C                                    BYTES/WORD
C                   0  0  1  0  2    IF BYTE HANDLING UNAVAILABLE
C_END_CCPMDE
C
C SPECIFICATION STATEMENTS
C ------------------------
C
C     .. Scalar Arguments ..
      INTEGER MODE,NBYT
C     ..
C     .. Local Arrays ..
      INTEGER MODES(0:4)
C     ..
C     .. Data statements ..
      DATA MODES/1,2,4,4,8/
C     ..
C
C---- Get number of bytes or words
C
      NBYT = -1
      IF (MODE.GE.0 .AND. MODE.LE.4) NBYT = MODES(MODE)
      END
C
C
C
C_BEGIN_CCPMVB
      SUBROUTINE CCPMVB(ARR1,I1,ARR2,I2,NTOMOV)
C     =========================================
C
C---- This subroutine moves bytes from one non-character array
C     to another. I must be implemented if ccpbyt returns .true.
C     but will otherwise be a dummy routine.
C
C Arguments:
C ==========
C
C        ARR1 (I/O) BYTE ARRAY(*): TO WHICH BYTES ARE TO BE COPIED
C          I1 (I)   INTEGER: THE START BYTE NUMBER IN ARR1 WHERE THE BYTES ARE
C                   TO BE COPIED
C        ARR2 (I)   BYTE ARRAY(*): FROM WHICH BYTES ARE TO BE COPIED
C          I2 (I)   THE START BYTE NUMBER IN ARR2 FROM WHICH THE BYTES
C                   ARE TO BE COPIED
C      NTOMOV (I)   INTEGER: THE NUMBER OF BYTES TO BE COPIED
C_END_CCPMVB
C
C     .. Scalar Arguments ..
      INTEGER I1,I2,NTOMOV
C     ..
C     .. Array Arguments ..
      INTEGER*1 ARR1(*),ARR2(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,N
C     ..
      I = I1 - 1
      J = I2 - 1
      DO 10 N = 1,NTOMOV
        I = I + 1
        J = J + 1
        ARR1(I) = ARR2(J)
   10 CONTINUE
C
      END
C
C
C_BEGIN_CCPMVI
      SUBROUTINE CCPMVI (IARR1,IARR2,NUM)
C     =================================
C
C  This routine assigns the first NUM words of IARR2 to IARR1
C
C Arguments:
C ==========
C
C    IARR1 (O)   INTEGER ARRAY(*)
C    IARR2 (O)   INTEGER ARRAY(*)
C     NUM (I)   Number of words to copy
C_END_CCPMVI
C
C  Arguments
      INTEGER NUM
      REAL IARR1(*),IARR2(*)
C
      INTEGER J
C
      DO 10 J=1,NUM
   10 IARR1(J)=IARR2(J)
      END
C
C
C_BEGIN_CCPMVR
      SUBROUTINE CCPMVR (ARR1,ARR2,NUM)
C     =================================
C
C  This routine assigns the first NUM elements of ARR2 to ARR1
C
C Arguments:
C ==========
C
C    ARR1 (O)   REAL ARRAY(*)
C    ARR2 (O)   REAL ARRAY(*)
C     NUM (I)   Number of words to copy
C_END_CCPMVR
C
C  Arguments
      INTEGER NUM
      REAL ARR1(*),ARR2(*)
C
      INTEGER J
C
      DO 10 J=1,NUM
   10 ARR1(J)=ARR2(J)
      END
C
C_BEGIN_CCPNUN
      INTEGER FUNCTION CCPNUN ()
C     ==========================
C
C     Return (the next) unused (not connected) i/o unit number.
C     Use this to select an arbitrary unit for i/o to avoid clashes with
C     other code.  (The value returned will be the same until the unit in
C     question is opened or a lower-numbered one is closed.)
C
C_END_CCPNUN
      LOGICAL OD, EX
      EXTERNAL CCPERR
      INTEGER IOS
C     The `standard' unit 5 and 6 may or may not be reported as open,
C     normally depending on whether an appropriate read or write has
C     happened, so we'll start at 7.  Lower-numbered ones might be used
C     for other things such as standard error.  99 seems a reasonable
C     place to stop.
      DO 10 CCPNUN=7,99
        INQUIRE (UNIT=CCPNUN, OPENED=OD, IOSTAT=IOS, EXIST=EX)
        IF (EX .AND. (.NOT.OD) .AND. IOS.EQ.0) RETURN
 10   CONTINUE
      CALL CCPERR (1, 'CCPNUN: Can''t find an unused unit')
      END
C
C
C
C_BEGIN_CCPONL
      LOGICAL FUNCTION CCPONL(IDUM)
C     =============================
C
C---- This function determines whether a program is being run on-line
C     if this information is available
C
C Arguments:
C ==========
C
C        IDUM (D)   DUMMY
C
C RETURNS .TRUE.  IF PROGRAM IS BEING RUN ON-LINE
C RETURNS .FALSE. IF BATCH MODE OR STATUS UNKNOWN
C_END_CCPONL
C
C     .. Scalar Arguments ..
      INTEGER IDUM
C     ..
C     .. Local Scalars ..
      INTEGER IYES,ITERM
C     ..
C     .. External Functions ..
      EXTERNAL UISATT
C     ..
C
C      test for fortran unit=6 o/p
C
      IYES = 0
      ITERM = 6
      CALL UISATT(ITERM,IYES)
      CCPONL = IYES.EQ.1
      END
C
C
C
C SUBROUTINE 'CCPPSF'
C ===================
C
C_BEGIN_CCPPSF
      SUBROUTINE CCPPSF(FILNAM,PATH,NAME,TYPE,VERS)
C     =============================================
C
C PARSE FILE NAME INTO COMPONENTS
C
C NOTE: THE ROUTINE  CONTAINS MACHINE DEPENDENT CODE
C
C
C Arguments:
C ==========
C
C      FILNAM (I)   CHARACTER*(*): FILE NAME STRING (NO EMBEDDED BLANKS ASSUMED)
C
C        PATH (O)   CHARACTER*(*): STRING RETURNING PATH OR, FOR VAX VMS,
C                   THE PART OF THE FILE SPECIFICATION UP TO THE
C                   END OF THE DIRECTORY SPECIFICATION (BLANK IF NONE)
C                   (INCLUDES TERMINATING ] or : or /)
C
C        NAME (O)   CHARACTER*(*): STRING RETURNING NAME.  (BLANK IF NONE)
C
C        TYPE (O)   CHARACTER*(*): STRING RETURNING FILE TYPE/EXTENSION
C                   (BLANK IF NONE)
C
C        VERS (O)   CHARACTER*(*): STRING RETURNING THE VERSION.
C                   (BLANK IF NONE)
C
C AFTER REMOVAL OF THE PATH PART OF THE STRING, IF PRESENT, THE VERSION ON
C A VAX IS TAKEN AS ANY TEXT FOLLOWING A SEMICOLON IN THE STRING OR, IF NO
C SEMICOLON IS PRESENT, ANY TEXT FOLLOWING THE LAST DOT IN THE STRING
C PROVIDED THAT AT LEAST TWO DOTS ARE PRESENT. ON A UNIX SYSTEM THE VERSION
C WILL ALWAYS BE RETURNED AS A BLANK.
C
C AFTER THE REMOVAL OF THE PATH AND VERSION PARTS OF THE STRING THEN, IF
C THERE IS AT LEAST ONE DOT, THE NAME IS THE STRING UP TO THE LAST DOT
C REMAINING AND THE TYPE IS THE PART OF THE STRING AFTER THE DOT. IF
C NO DOT IS PRESENT THEN THE REMAINING STRING IS THE NAME AND THE TYPE
C IS BLANK.
C_END_CCPPSF
C
C SPECIFICATION STATEMENTS
C ------------------------
C
      CHARACTER*(*) FILNAM,PATH,NAME,TYPE,VERS
      INTEGER LMAX,LMIN,L,LSC,LDOT,NDOT,LENSTR
      EXTERNAL VAXVMS, WINMVS, RTNBKS, LENSTR
      LOGICAL VAXVMS, WINMVS, VMS, MVS
      CHARACTER RTNBKS*1, BKS*1
C
C INITIALISATIONS
C ---------------
C
      PATH=' '
      NAME=' '
      TYPE=' '
      VERS=' '
      LMAX=LENSTR(FILNAM)
      IF (LMAX.EQ.0) RETURN
      LMIN=0
      VMS = VAXVMS()
      MVS = WINMVS()
      BKS = RTNBKS()      
10    LMIN=LMIN+1
      IF (FILNAM(LMIN:LMIN).EQ.' ') GO TO 10
C
C GET PATH
C --------
C
      IF (VMS) THEN
        DO 20 L=LMAX,LMIN,-1
          IF (FILNAM(L:L).EQ.':'.OR.FILNAM(L:L).EQ.']') GO TO 30
 20     CONTINUE
      ELSEIF (MVS) THEN
        DO 21 L=LMAX,LMIN,-1
          IF (FILNAM(L:L).EQ.BKS .OR. FILNAM(L:L).EQ.'/')GO TO 30
 21     CONTINUE
      ELSE
        DO 22 L=LMAX,LMIN,-1
          IF (FILNAM(L:L).EQ.'/')GO TO 30
 22     CONTINUE
      ENDIF
      GO TO 40
30    PATH=FILNAM(LMIN:L)
      LMIN=L+1
      IF (LMIN.GT.LMAX) RETURN
C
C GET VERSION IF PRESENT
C ----------------------
C
 40   CONTINUE
      IF (VMS) THEN
        LSC=INDEX(FILNAM(LMIN:LMAX),';')
        IF (LSC.GT.0) THEN
          LSC=LSC+LMIN-1
          IF (LSC.LT.LMAX) VERS=FILNAM(LSC+1:LMAX)
          LMAX=LSC-1
        ELSE
          LDOT=0
          NDOT=0
          DO 50 L=LMAX,LMIN,-1
            IF (FILNAM(L:L).EQ.'.') THEN
              NDOT=NDOT+1
              IF (LDOT.EQ.0) LDOT=L
            ENDIF
 50       CONTINUE
          IF (NDOT.GT.1) THEN
            IF (LDOT.LT.LMAX) VERS=FILNAM(LDOT+1:LMAX)
            LMAX=LDOT-1
          ENDIF
        ENDIF
      ENDIF
C
C GET NAME AND TYPE
C -----------------
C
      IF (LMAX.LT.LMIN) RETURN
      LDOT=0
      DO 60 L=LMAX,LMIN,-1
      IF (FILNAM(L:L).EQ.'.') THEN
         LDOT=L
         GO TO 70
      ENDIF
60    CONTINUE
70    IF (LDOT.EQ.0) THEN
         NAME=FILNAM(LMIN:LMAX)
         RETURN
      ELSE
         IF (LDOT.GT.LMIN) NAME=FILNAM(LMIN:LDOT-1)
         IF (LDOT.LT.LMAX) TYPE=FILNAM(LDOT+1:LMAX)
      ENDIF
      END
C
C
C_BEGIN_CSETNV
      SUBROUTINE CSETNV(LNAME,FILNAM,ENAME,ETYPE,EXTN,ICOUNT,LSKIP)
C     =============================================================
C
C     Associate `logical name' LNAME with value FILNAM.  It is passed
C     arrays of (name, type, extension) for ICOUNT number of name lines
C     read from environ.def.  Doesn't re-define existing name if LSKIP is true.
C
C Arguments:
C ==========
C
C   LNAME (I)   CHARACTER*(*): Logical name (environment variable).
C
C  FILNAM (I/O) CHARACTER*(*): File name, if extension is omitted it is appended
C
C   ENAME (I/O) CHARACTER(150)*20 ARRAY: containing list of environment
C               variables; if LNAME is not in list it is appended
C               (also to ETYPE & EXTN arrays).
C
C   ETYPE (I,O) CHARACTER(150)*5 ARRAY: containing list of in/out types.
C
C    EXTN (I/O) CHARACTER(150)*4 ARRAY: containing list of extensions.
C
C  ICOUNT (I/O) INTEGER: Length of arrays ENAME, ETYPE & EXTN.
C
C   LSKIP (I)   LOGICAL: If .TRUE. existing name not re-defined.
C
C_END_CSETNV
C
C     .. Parameters ..
      INTEGER ILIMIT,ISTRLN,IENV
      PARAMETER (ILIMIT=150,ISTRLN=200,IENV=20)
C     ..
C     .. Scalar Arguments ..
      INTEGER ICOUNT
      CHARACTER LNAME* (*),FILNAM* (*)
      LOGICAL LSKIP
C     ..
C     .. Array Arguments ..
      CHARACTER ENAME(ILIMIT)* (IENV),ETYPE(ILIMIT)* (5),
     +          EXTN(ILIMIT)* (4)
C     ..
C     .. Local Scalars ..
      INTEGER I,II,ISTAT,JJ,PROCID
      LOGICAL VAX,MVS,EXIST
      CHARACTER ERRSTR* (ISTRLN),LIBFIL* (ISTRLN),PROGNM* (ISTRLN),
     +          TMPNAM* (ISTRLN),LINE* (ISTRLN),SCRFIL* (ISTRLN),
     +          BKS*(1)
C     ..
C     .. External Functions ..
      INTEGER GETPID,LENSTR
      LOGICAL VAXVMS, WINMVS
      CHARACTER FDIR* (ISTRLN),FEXTN* (ISTRLN),FROOT* (ISTRLN), RTNBKS
      EXTERNAL LENSTR,VAXVMS,FDIR,FEXTN,FROOT
C     ..
C     .. External Subroutines ..
      EXTERNAL CCPERR,UGTARG,QPRINT,UGTENV,USTENV
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC INDEX
C     ..
      SAVE
      DATA PROGNM/' '/
C
C---- Check Logical Name does not already exist (unless processing
C     command line, in which case LSKIP will be true to override
C     environment)
C
      CALL UGTENV(LNAME,TMPNAM)
      IF (TMPNAM.NE.' ' .AND. LSKIP ) RETURN
      VAX = VAXVMS()
      MVS = WINMVS()
      BKS = RTNBKS()
C
C---- Get program name (argv[0]), but check if we have it already
C
      IF (PROGNM.EQ.' ') THEN
        CALL UGTARG(0,TMPNAM)
        PROGNM = FROOT(TMPNAM)
      ENDIF
C
C---- look through list for a match (possibly abbreviated) [is this
C     abbreviation possibility documented?]
C
      DO 10 JJ = 1,ICOUNT
        IF (ENAME(JJ).EQ.LNAME(1:LENSTR(ENAME(JJ)))) GO TO 20
   10 CONTINUE
C
C---- Unknown logical name add it to the list.
C
      TMPNAM = 'Non standard logical name '
      TMPNAM(27:) = LNAME
      CALL QPRINT(2,TMPNAM)
      ICOUNT = ICOUNT + 1
      IF (ICOUNT.GT.ILIMIT)
     +     CALL CCPERR(1,'Too many logical names')
      ENAME(ICOUNT) = LNAME
      ETYPE(ICOUNT) = 'undef'
      IF (LENSTR(FEXTN(FILNAM)).GT.4)
     +  CALL CCPERR(2,
     +  'Extension too long in s/r CSETNV: '//FEXTN(FILNAM))
      EXTN(ICOUNT) = FEXTN(FILNAM)
      JJ = ICOUNT
C
C---- Known logical name processing
C
   20 IF (FEXTN(FILNAM).EQ.' ') THEN
C
C---- Add extension
C
        IF (FILNAM.EQ.'/dev/null' .OR. FILNAM.EQ.'NL:') THEN
C          but not if FILNAM is /dev/null or NL:
          GOTO 333
        ELSE
          II = LENSTR(FILNAM) + 1
          FILNAM(II:) = EXTN(JJ)
        ENDIF
      ENDIF
      IF (FDIR(FILNAM).EQ.' ') THEN
CCC       this didn't agree with documentation:
CCC        IF (EXTN(JJ).EQ.'.lib' .OR. EXTN(JJ).EQ.'.prt' .OR.
CCC     +      EXTN(JJ).EQ.'.bes' .OR. EXTN(JJ).EQ.'.dic') THEN
        TMPNAM = FEXTN(FILNAM)
        IF (VAX) CALL CCPLWC(TMPNAM)
        IF (TMPNAM.EQ.'lib' .OR. TMPNAM.EQ.'prt' .OR.
     +      TMPNAM.EQ.'bes' .OR. TMPNAM.EQ.'dic') THEN
C         look for files without path but with standard extension in the
C         standard place
          CALL UGTENV('CLIBD',LIBFIL)
C         add the standard directory qualifier
          IF (VAX) THEN
C           fixme: should we insist that VMS defines CLIBD as well as un*x?
            IF (LIBFIL.NE.' ') THEN
              TMPNAM = 'CLIBD:'
              TMPNAM(7:) = FILNAM
            ELSE
              TMPNAM = FILNAM
            ENDIF
          ELSEIF (MVS) THEN
            IF (LIBFIL.EQ.' ') CALL CCPERR(1,'CLIBD not defined')
            II = LENSTR(LIBFIL)
            TMPNAM = LIBFIL(:II)//BKS
            II = II + 2
            TMPNAM(II:) = FILNAM
          ELSE
            IF (LIBFIL.EQ.' ') CALL CCPERR(1,'CLIBD not defined')
            II = LENSTR(LIBFIL)
            TMPNAM = LIBFIL(:II)//'/'
            II = II + 2
            TMPNAM(II:) = FILNAM
          END IF
          FILNAM = TMPNAM
        ELSE IF (EXTN(JJ).EQ.'.scr' .OR. FEXTN(FILNAM).EQ.'scr') THEN
C         scratch files in a special place
C         actually create <ccp4_scr>/<prognm>_.<pid>
          CALL UGTENV('CCP4_SCR',TMPNAM)
          IF (VAX) THEN
            IF (TMPNAM.EQ.' ') THEN
              TMPNAM = PROGNM
            ELSE
              TMPNAM = 'CCP4_SCR:' // PROGNM
            ENDIF
          ELSEIF (MVS) THEN
            IF (TMPNAM.EQ.' ') CALL CCPERR(1,'CCP4_SCR not defined')
            II = LENSTR(TMPNAM) + 1
            TMPNAM(II:) = BKS//PROGNM
          ELSE
            IF (TMPNAM.EQ.' ') CALL CCPERR(1,'CCP4_SCR not defined')
            II = LENSTR(TMPNAM) + 1
            TMPNAM(II:) = '/'//PROGNM
          END IF
          II = LENSTR(TMPNAM) + 1
          TMPNAM(II:II) = '_'
          II = II + 1
          I = INDEX(FILNAM,'.')
          TMPNAM(II:) = FILNAM(:I)
          IF (VAX) THEN
            WRITE (SCRFIL,'(Z8.8)') GETPID()
          ELSE
            PROCID = MOD(GETPID(),100000)
C     Windows98 does not return a pid so make a number up
            IF (PROCID.GT.99999 .OR. PROCID.LE.0) THEN
              CALL USTIME(PROCID)
              PROCID = MOD(PROCID,100000)
            ENDIF
            WRITE (SCRFIL,'(I5.5)') PROCID
          ENDIF
          FILNAM = TMPNAM(1:LENSTR(TMPNAM))//SCRFIL
        END IF
      END IF
333   CONTINUE
C
C---- Now test input files do exist (but not for defaults, to avoid
C     checking 40 or 50 files listed in default.def which the setup
C     should guarantee)
C
      IF (ETYPE(JJ).EQ.'in' .AND. .NOT.LSKIP) THEN
        INQUIRE(FILE=FILNAM,EXIST=EXIST)
        IF (.NOT.EXIST) THEN
          ERRSTR = 'Cannot find file '
          ERRSTR(18:) = FILNAM
          CALL CCPERR (-1,ERRSTR)
        END IF
      END IF
      II = LENSTR(LNAME) + 1
      LINE = LNAME
      LINE(II:II) = '='
      II = II + 1
      LINE(II:) = FILNAM
C     =======================================
      CALL USTENV(LINE(1:LENSTR(LINE)),ISTAT)
C     =======================================
      IF (ISTAT.NE.0) THEN
        IF (VAX) THEN
          ERRSTR = 'Cannot create environment variable '
        ELSE
          ERRSTR = 'Cannot create logical name '
        ENDIF
        ERRSTR(36:) = LNAME
        CALL CCPERR (-1,ERRSTR)
      END IF
      CALL QPRINT(3,LINE)
      END
C
C
C
C_BEGIN_CCPPAG
      SUBROUTINE CCPPAG(IUN,NCOL,NLIN)
C     ================================
C
C---- This subroutine returns the number of columns and lines
C     for a printer output page on a given fortran unit number
C     if the information is available
C
C Arguments:
C ==========
C
C         IUN (I)   INTEGER: FORTRAN UNIT NUMBER
C        NCOL (O)   INTEGER: NUMBER OF COLUMNS IN THE PAGE
C        NLIN (O)   INTEGER: NUMBER OF LINES IN THE PAGE
C
C Return 80,132 unless a terminal whence 0,80
C_END_CCPPAG
C
C     .. Scalar Arguments ..
      INTEGER IUN,NCOL,NLIN
C     ..
C     .. Local Scalars ..
      INTEGER IYES
C     ..
C     .. External Subroutines ..
      EXTERNAL UISATT
C     ..
      CALL UISATT(IUN,IYES)
      IF (IYES.EQ.1) THEN
        NLIN = 0
        NCOL = 80
      ELSE
        NLIN = 80
        NCOL = 132
      END IF
      END
C
C
C
C SUBROUTINE 'CCPSI2'
C ===================
C
C_BEGIN_CCPSI2
      SUBROUTINE CCPSI2(IVAL,IA,N)
C     ============================
C
C SET AN INTEGER VALUE FROM 0 TO 65535 INTO THE N'TH UNSIGNED INTEGER*2 ELEMENT
C OF AN INTEGER (OR OTHER) ARRAY.
C NOTE: NO OVERFLOW CHECKING IS DONE.
C
C (MUST BE IMPLEMENTED IF CCPBYT FUNCTION RETURNS .TRUE.)
C [for LAUE]
C
C Arguments:
C ==========
C
C    IVAL (I)   INTEGER: VALUE FROM 0 TO 65535
C
C      IA (I/O) INTEGER*2 ARRAY: WHERE THE UNSIGNED INTEGER*2 VALUE IS TO BE
C               INSERTED
C
C       N (I)   INTEGER: THE POSITION IN 'IA' WHERE THE UNSIGNED INTEGER*2
C               VALUE IS TO BE INSERTED
C_END_CCPSI2
C
C SPECIFICATION STATEMENTS
C ------------------------
C
      INTEGER*2 IA(*)
      INTEGER*2 JBYT(2)
      INTEGER   JA,IVAL,N
      EQUIVALENCE (JA,JBYT(1))
      LOGICAL CALLED, LITEND
      EXTERNAL LITEND
      INTEGER IND
      SAVE CALLED, IND
      DATA CALLED/.FALSE./
C
      IF (.NOT.CALLED) THEN
        IF (LITEND(1)) THEN
          IND = 1
        ELSE
          IND = 2
        ENDIF
        CALLED=.TRUE.
      ENDIF
C
C SET UNSIGNED INTEGER*2
C ----------------------
C
      JA=IVAL
      IA(N)=JBYT(IND)
      END
C
C
C
C SUBROUTINE 'CCPSTB'
C ===================
C
C_BEGIN_CCPSTB
      SUBROUTINE CCPSTB(IVAL,IA,N)
C     ============================
C
C SET AN INTEGER VALUE FROM 0 TO 255 INTO THE N'TH BYTE OF AN INTEGER
C (OR OTHER) ARRAY.
C NOTE: NO OVERFLOW CHECKING IS DONE.
C
C (MUST BE IMPLEMENTED IF CCPBYT FUNCTION RETURNS .TRUE.)
C [for LAUE]
C
C Arguments:
C ==========
C
C    IVAL (I)   INTEGER: VALUE FROM 0 TO 255
C      IA (I/O) BYTE ARRAY(*): WHERE THE BYTE VALUE IS TO BE INSERTED
C       N (I)   INTEGER: THE POSITION IN 'IA' WHERE THE BYTE VALUE IS TO
C               BE INSERTED
C_END_CCPSTB
C
C SPECIFICATION STATEMENTS
C ------------------------
C
      INTEGER*1 IA(*)
      INTEGER*1 JBYT(4)
      EQUIVALENCE (JA,JBYT(1))
      EQUIVALENCE (JA,JBYT(1))
      LOGICAL CALLED, LITEND
      EXTERNAL LITEND
      INTEGER IND
      SAVE CALLED, IND
      DATA CALLED/.FALSE./
C
      IF (.NOT.CALLED) THEN
        IF (LITEND(1)) THEN
          IND = 1
        ELSE
          IND = 4
        ENDIF
        CALLED=.TRUE.
      ENDIF
C
C SET BYTE
C --------
C
      JA=IVAL
      IA(N)=JBYT(IND)
      END
C
C
C
C_BEGIN_CCPSUM
      REAL FUNCTION CCPSUM(A,N,L)
C     ===========================
C
C---- This function sums the elements of an array. (for the cray this
C     function will call the cray 'ssum' function)
C
C Arguments:
C ==========
C
C           A (I)   REAL ARRAY(N): ARRAY TO BE SUMMED
C           N (I)   INTEGER: NO. OF ELEMENTS IN THE ARRAY
C           L (I)   INTEGER: SUM EVERY L'TH ELEMENT
C
C  CCPSUM RETURNS THE SUM
C_END_CCPSUM
C
C SPECIFICATION STATEMENTS AND CODE
C
C     .. Scalar Arguments ..
      INTEGER L,N
C     ..
C     .. Array Arguments ..
      REAL A(N)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      CCPSUM = 0.0
      DO 10 I = 1,N,L
        CCPSUM = A(I) + CCPSUM
   10 CONTINUE
      END
C
C
C
C_BEGIN_CCPTOI
      SUBROUTINE CCPTOI(ARRAY,N,II,ITYP,IFAIL)
C     ========================================
C
C---- This subroutine converts the n'th byte or integer*2 element in a
C     non-character array to an integer value. it is used by the
C     map file handling routines and must be implemented if map modes
C     0,1,3 or 5 are to be used.
C
C Arguments:
C ==========
C
C       ARRAY (I)   REAL ARRAY(*): CONTAINING THE ELEMENTS TO BE CONVERTED
C
C           N (I)   INTEGER: THE NUMBER OF THE ELEMENT TO BE CONVERTED
C
C          II (O)   INTEGER: THE CALCULATED INTEGER VALUE (FOR BYTES THIS WILL
C                   BE IN THE RANGE 0-255)
C
C        ITYP (I)   INTEGER: THE CONVERSION TYPE =1, BYTE TO INTEGER
C                                                =2, INTEGER*2 TO INTEGER
C
C       IFAIL (I/O) INTEGER: ON INPUT   =0, STOP IF CONVERSION NOT AVAILABLE
C                                       =1, RETURN FROM SUBROUTINE ALWAYS
C                            ON OUTPUT  UNCHANGED IF CONVERSION CARRIED OUT
C                                       =-1 IF CONVERSION NOT AVAILABLE
C_END_CCPTOI
C
C     .. Scalar Arguments ..
      INTEGER IFAIL,II,ITYP,N
C     ..
C     .. Array Arguments ..
      REAL ARRAY(*)
C     ..
C     .. Local Scalars ..
      REAL RR
      INTEGER IA,NB,NIH,NW
C     ..
C     .. Local Arrays ..
      INTEGER*1 IBYT(4),JBYT(4)
      INTEGER*2 JHALF(2)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
C     .. Equivalences ..
      EQUIVALENCE (IA,IBYT(1))
      EQUIVALENCE (RR,JHALF(1),JBYT(1))
C     ..
      LOGICAL CALLED, LITEND
      INTEGER IND
      EXTERNAL LITEND, CCPERR
      SAVE CALLED, IND
      DATA CALLED/.FALSE./
C
      IF (.NOT.CALLED) THEN
        CALLED=.TRUE.
        IF (LITEND(1)) THEN
          IND = 1
        ELSE
          IND = 4
        ENDIF
      ENDIF
C
      GO TO (10,20) ITYP
C
C---- Byte to integer value
C
   10 NW = (N-1)/4 + 1
      NB = MOD(N-1,4) + 1
      IA = 0
      RR = ARRAY(NW)
      IBYT(IND) = JBYT(NB)
      II = IA
      IF (II.LT.0 .OR. II.GT.255) THEN
        IF (IFAIL .EQ. 0) THEN
          CALL CCPERR(1,' *** Error in CCPTOI, bad convertion ***') 
        ELSE
          IFAIL = -1
        ENDIF
      ENDIF
      RETURN
C
C---- Integer*2 to integer value
C
   20 NW = (N-1)/2 + 1
      NIH = MOD(N-1,2) + 1
      RR = ARRAY(NW)
      II = JHALF(NIH)
      IF (II.LT.0 .OR. II.GT.65535) THEN
        IF (IFAIL .EQ. 0) THEN
          CALL CCPERR(1,' *** Error in CCPTOI, bad convertion ***')
        ELSE
          IFAIL = -1
        ENDIF
      ENDIF
      END
C
C
C
C_BEGIN_CCPUFL
      SUBROUTINE CCPUFL
C     =================
C
C---- This subroutine is called to suppress underflow error messages
C     if required and if the facility is available.
C
C Arguments:  NONE
C ==========
C
C----  Not implemented.
C_END_CCPUFL
      END
C
C
C
C_BEGIN_CCP4_PATCH_LEVEL
      SUBROUTINE CCP4_PATCH_LEVEL(PATCH_LEVEL)
C     =================================
C
C---- Return current CCP4 patch level as string
C
C Arguments:
C ==========
C
C      PATCH_LEVEL  (O)   CHARACTER*(*): current patch level of CCP4 suite
C_END_CCP4_PATCH_LEVEL
C
C     .. Scalar Arguments ..
      CHARACTER*(*) PATCH_LEVEL

      PATCH_LEVEL = '5.0f'

      END
C
C
C
C_BEGIN_CCPZBI
      SUBROUTINE CCPZBI (ARR1,NUM)
C     ============================
C
C  This routine zeros NUM bytes of the array ARR1
C
C Arguments:
C
C    ARR1 (O)   BYTE ARRAY(*): array to be zeroed
C     NUM (I)   INTEGER: Number of bytes
C_END_CCPZBI
C
C  Arguments ......
      INTEGER NUM
      INTEGER*1 ARR1(*)
C
      INTEGER J
C
      DO 10 J=1,NUM
   10 ARR1(J)=0
      END
C
C
C_BEGIN_CCPZI
      SUBROUTINE CCPZI (IARR1,NUM)
C     ===========================
C
C  This routine assigns zero to IARR1 using NUM words
C
C Arguments:
C
C    IARR1 (O)   INTEGER ARRAY(*): array to be zeroed
C     NUM (I)   INTEGER: Number of words
C_END_CCPZI
C
C  Arguments ..........
      INTEGER NUM, IARR1(*)
C
      INTEGER J
C
      DO 10 J=1,NUM
   10 IARR1(J)=0
      END
C
C
C_BEGIN_CCPZR
      SUBROUTINE CCPZR (ARR1,NUM)
C     ===========================
C
C  This routine assigns zero to ARR1 using NUM words
C
C Arguments:
C
C    ARR1 (O)   REAL ARRAY(*): array to be zeroed
C     NUM (I)   INTEGER: Number of words
C_END_CCPZR
C
C  Arguments ..........
      INTEGER NUM
      REAL ARR1(*)
C
      INTEGER J
C
      DO 10 J=1,NUM
   10 ARR1(J)=0.0
      END
C
C
C     ===================================
C_BEGIN_FDIR
      FUNCTION FDIR(FILNAM)
      CHARACTER*(*) FDIR
C     ===================================
C
C---- Returns the path (directory) of a file name or ' '
C
C Arguments:
C
C  FILNAM (I)   CHARACTER*(*): File name
C_END_FDIR
      CHARACTER FILNAM* (*)
      CHARACTER*1 NAME, TYPE, VERS
      EXTERNAL CCPPSF
C
      CALL CCPPSF(FILNAM, FDIR, NAME, TYPE, VERS)
      END
C
C
C     ====================================
C_BEGIN_FEXTN
      FUNCTION FEXTN(FILNAM)
      CHARACTER*(*) FEXTN
C     ====================================
C
C---- Returns the extension of a file name or ' '
C
C Arguments:
C
C  FILNAM (I)   CHARACTER*(*): File name
C_END_FEXTN
      CHARACTER FILNAM* (*)
      CHARACTER*1 PATH, NAME, VERS
      EXTERNAL CCPPSF
C
      CALL CCPPSF(FILNAM, PATH, NAME, FEXTN, VERS)
      END
C
C
C     ====================================
C_BEGIN_FROOT
      FUNCTION FROOT(FILNAM)
      CHARACTER*(*) FROOT
C     ====================================
C
C---- Returns a file name minus an extension.
C
C Arguments:
C
C  FILNAM (I)   CHARACTER*(*): File name
C_END_FROOT
C
      CHARACTER FILNAM* (*)
      CHARACTER*1 PATH, TYPE, VERS
      EXTERNAL CCPPSF
C
      CALL CCPPSF(FILNAM, PATH, FROOT, TYPE, VERS)
      END
C
C
C
C_BEGIN_LITEND
         LOGICAL FUNCTION LITEND(IDUM)
C        =============================
C
C---- Check endedness, Returns TRUE if little endian (VAX, FX2800,
C                                                   Ultrix)
C                              FALSE if big endian (IBM,IRIS,ESV)
C
C Arguments:
C ==========
C
C    IDUM (D)   DUMMY
C_END_LITEND
C
         INTEGER I, IDUM
         INTEGER*1 B(4)
         EQUIVALENCE (I,B(1))
C
C---- Initialise B
C
          DO 10 JDO=1,4
            B(JDO) = 0
 10       CONTINUE
C
          I = 1
C
          IF (B(1) .NE. 0) THEN
              LITEND = .TRUE.
          ELSE
              LITEND = .FALSE.
          END IF
C
        END
C      
C^L
C====================================================================== 
C                      
C_BEGIN_LENSTR
      INTEGER FUNCTION LENSTR(STRING)
C     ===============================
C      
C---- Returns significant string length excluding trailing spaces       
C                     
C     NB: LENSTR removes trailing spaces, plus trailing "null"
C     characters (ascii code 0) and carriage-returns (ascii code 13).
C     Carriage-returns may be inserted by editing on non-unix
C     systems? (PJX)
C      
C Arguments:
C ==========
C      
C  STRING (I)   CHARACTER*(*): Input string
C_END_LENSTR
C      
C     .. Scalar Arguments ..
      CHARACTER STRING* (*)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LEN
C     ..
      LENSTR = LEN(STRING)
 10   CONTINUE
      IF (LENSTR.NE.0) THEN
        IF(  ICHAR(STRING(LENSTR:LENSTR)).EQ.0   .OR.
     .       ICHAR(STRING(LENSTR:LENSTR)).EQ.13  .OR.
     .             STRING(LENSTR:LENSTR) .EQ.' '      ) THEN
          LENSTR = LENSTR - 1
          GO TO 10
        END IF
      END IF
C      
      END
C
C
C FUNCTION 'LUNSTI'
C =================
C
C_BEGIN_LUNSTI
      FUNCTION LUNSTI(IDUM)
C     =====================
C
C Returns the fortran standard input unit number
C
C Arguments:
C ==========
C
C       IDUM (D)   Dummy
C_END_LUNSTI
C
      INTEGER LUNSTI,IDUM
C
      LUNSTI = 5
      END
C
C
C FUNCTION 'LUNSTO'
C =================
C
C_BEGIN_LUNSTO
      FUNCTION LUNSTO(IDUM)
C     =====================
C
C Returns the fortran standard output unit number
C
C Arguments:
C ==========
C
C       IDUM (I)   Dummy argument
C_END_LUNSTO
C
      INTEGER LUNSTO,IDUM
C
      LUNSTO = 6
      END
C
C
C FUNCTION 'NBITST'
C =================
C
C_BEGIN_NBITST
      FUNCTION NBITST(IWORD,LSB,NBITS)
C     ================================
C
C Return the (unsigned) integer value held within a bit field in a word
C [for LAUE]
C
C Arguments:
C ==========
C
C      IWORD (I)    INTEGER: The word containing the bits to be examined
C
C        LSB (I)    INTEGER: The least significant bit offset for the bit field
C
C      NBITS (I)    INTEGER: The number of bits in the bit field (Must be less
C                   than the word length)
C_END_NBITST
C
C====== Get the bit value
C
      INTEGER IWORD,LSB,NBITS,IVAL,KMSK,KVAL
C
      KMSK = 2**NBITS - 1
      NBITST = IAND(ISHFT(IWORD,-LSB),KMSK)
      END
C
C
C     =======================
C_BEGIN_NOCRLF
      SUBROUTINE NOCRLF(LINE)
C     =======================
C
C---- Output a line supressing cr/lf.
C
C Arguments:
C ==========
C
C    LINE (I)   CHARACTER*(*): Line to output.
C_END_NOCRLF
C
      EXTERNAL LUNSTO,TTSEND
      INTEGER LUNSTO
      CHARACTER*(*) LINE
      CALL TTSEND(LUNSTO(1),LINE,0)
      END
C
C
C SUBROUTINE 'STBITS'
C ===================
C
C_BEGIN_STBITS
      SUBROUTINE STBITS (IWORD,LSB,NBITS,IVAL)
C     ========================================
C Set a bit field within a word to a given (unsigned) integer value
C [for LAUE]
C
C Arguments:
C ==========
C
C      IWORD (I/O)  INTEGER: The word in which the bits are to be set
C
C        LSB (I)    INTEGER: The least significant bit offset for the bit field
C
C      NBITS (I)    INTEGER: The number of bits in the bit field (must be less
C                   than the word length)
C
C       IVAL (I)    INTEGER: The unsigned integer value to be set in the bit
C                   field (The user should ensure that this value will
C                   fit within the requested bit field)
C_END_STBITS
C
C====== Set the bits
C
      KMSK = 2**NBITS - 1
      KVAL = IVAL
      KMSK = ISHFT(KMSK,LSB)
      KMSK = NOT(KMSK)
      KVAL = ISHFT(KVAL,LSB)
      IWORD = IOR(IAND(IWORD,KMSK),KVAL)
      END
C      
C^L
C     ============================
      LOGICAL FUNCTION HKLEQ(IH,KH)
C     =============================    
C    
C---- Returns true if indices ih = kh
C        
C     .. Array Arguments ..   
      INTEGER IH(3),KH(3)
C     ..
C
      HKLEQ = .FALSE.
C 
      IF (IH(1).EQ.KH(1) .AND. IH(2).EQ.KH(2) .AND. 
     +    IH(3).EQ.KH(3)) HKLEQ = .TRUE.
C 
      END
C  
C Group 7: Subroutines for generating and accessing a hash table:
C======================================================================
C 
C---- CCP4_HASH_SETUP CCP4_HASH_LOOKUP CCP4_HASH_ZEROIT
C   
C         Routines and functions used to initialise, set up and access
C         an internal look-up table. Not clear why these routines are
C         here in particular.
C    
C---- SUBROUTINE CCP4_HASH_SETUP(NSER,NFIND)
C 
C---- This subroutine sets up a value for the function ccp4_hash_lookup
C     when ccp4_hash_lookup(nser) is later evaluated it will return nfind
C     this function will allow the efficient retrieval of an identifier
C     for a large range variable (such as a crystal number).  the values
C     of the function ccp4_hash_lookup(nser) are stored in the array
C     it(2, kpri) where kpri is the prime number used to generate the
C     function.
C     The array it  lives in the common look which is shared by
C     ccp4_hash_setup and the function ccp4_hash_lookup
C    
C     NOTES: A hash table is a way of storing information so that it
C     easily be retrieved without the need for indexing or long searches.
C     NSER is referred to as the "key", which is "hashed" (computer-
C     science speak for "messed up") by the hashing function (in this
C     case MOD(NSER4,KPRI) + 1) to determine where the value pair will
C     be stored. The function LOOKUP can then search on the same basis
C     when supplied with the key, to retreive the pair in (at most) 3
C     calculations. Note that KPRI (the table size) MUST BE A PRIME in
C     order for this method to work.   
C                  
C     IT(1, NDX) = NSER,  IT(2, NDX) = NFIND
C    
C---- INTEGER FUNCTION CCP4_HASH_LOOKUP(NSER)
C
C---- The function ccp4_hash_lookup returns the value nfind (which was
C     input when setting up the function in the subroutine ccp4_hash_setup)
C     for the large range variable nser.  Uses hashing. (see comments for
C     CCP4_HASH_SETUP for description of hashing method).      
C                       
C---- SUBROUTINE CCP4_HASH_ZEROIT()      
C
C     Initialises elements of array it used in ccp4_hash_setup and
C     ccp4_hash_lookup to zero.

C^L
C     =======================================
      INTEGER FUNCTION CCP4_HASH_LOOKUP(NSER)
C     =======================================
C
C---- The function ccp4_hash_lookup returns the value nfind (which was
C     input when setting up the function in the subroutine ccp4_hash_setup)
C     for the large range variable nser.  Uses hashing. (see comments for
C     CCP4_HASH_SETUP for description of hashing method).
C 
      IMPLICIT NONE
C     .. Parameter (table size: MUST BE A PRIME NUMBER)
      INTEGER KPRI
      PARAMETER (KPRI=12007)
C
C     .. Scalar Arguments ..
      INTEGER NSER
C     ..
C     .. Arrays in Common ..
      INTEGER IT
C     ..
C     .. Local Scalars ..
      INTEGER NDX,NSER4
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
C     .. Common blocks ..  
      COMMON /LOOK/IT(2,KPRI)
      SAVE /LOOK/
C     .. 
C
      NSER4 = NSER
C     
   10 CONTINUE
C
      NDX = MOD(NSER4,KPRI) + 1
      IF (NSER.NE.IT(1,NDX)) THEN
        IF (IT(1,NDX).NE.0) THEN
          NSER4 = NSER4 + 3
          GO TO 10
        END IF
      END IF
C 
      CCP4_HASH_LOOKUP = IT(2,NDX)
C
      END
C
C^L
C     ======================================
      SUBROUTINE CCP4_HASH_SETUP(NSER,NFIND)
C     ======================================
C     
C---- This subroutine sets up a value for the function ccp4_hash_lookup
C     when ccp4_hash_lookup(nser) is later evaluated it will return nfind
C     this function will allow the efficient retrieval of an identifier
C     for a large range variable (such as a crystal number).  the values
C     of the function ccp4_hash_lookup(nser) are stored in the array
C     it(2, kpri) where kpri is the prime number used to generate the
C     function
C     The array it  lives in the common look which is shared by
C     ccp4_hash_setup and the function ccp4_hash_lookup
C
C     NOTES: A hash table is a way of storing information so that it
C     easily be retrieved without the need for indexing or long searches.
C     NSER is referred to as the "key", which is "hashed" (computer-
C     science speak for "messed up") by the hashing function (in this
C     case MOD(NSER4,KPRI) + 1) to determine where the value pair will
C     be stored. The function LOOKUP can then search on the same basis
C     when supplied with the key, to retreive the pair in (at most) 3      
C     calculations. Note that KPRI (the table size) MUST BE A PRIME in   
C     order for this method to work.
C
C     IT(1, NDX) = NSER,  IT(2, NDX) = NFIND
C 
      IMPLICIT NONE
C     .. Parameter (table size: MUST BE A PRIME NUMBER)
      INTEGER KPRI
      PARAMETER (KPRI=12007)
C
C     .. Scalar Arguments ..
      INTEGER NFIND,NSER
C     ..
C     .. Arrays in Common ..
      INTEGER IT
C     ..
C     .. Local Scalars ..
      INTEGER NDX,NSER4
      CHARACTER STROUT*140
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
C     .. Common blocks ..  
      COMMON /LOOK/IT(2,KPRI)
      SAVE /LOOK/
C     .. 
C
      NSER4 = NSER
   10 CONTINUE
      NDX = MOD(NSER4,KPRI) + 1  
      IF ((NSER4-NSER) .GE. 3*KPRI) THEN
         WRITE (STROUT, '(A,I8)')
     $     ' **** Error in SETUP: overflowed hash table, size ', KPRI
         CALL PUTLIN(STROUT,'CURWIN')
         CALL CCPERR(1,'*** Filled hash table in SETUP ***')
      ENDIF
      IF (IT(1,NDX).NE.0) THEN
        NSER4 = NSER4 + 3
        GO TO 10
      END IF
C 
      IT(1,NDX) = NSER
      IT(2,NDX) = NFIND
      RETURN
      END
C
C^L
C         
C     =============================
      SUBROUTINE CCP4_HASH_ZEROIT()
C     =============================
C 
      IMPLICIT NONE
C     .. Parameter (table size: MUST BE A PRIME NUMBER)
      INTEGER KPRI
      PARAMETER (KPRI=12007)
C     
C     .. Arrays in Common ..
      INTEGER IT
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     ..
C     .. Common blocks ..
      COMMON /LOOK/IT(2,KPRI)
      SAVE /LOOK/
C     
      DO 20 I = 1,KPRI
        IT(1,I) = 0
        IT(2,I) = 0
   20 CONTINUE
      RETURN
C     
      END
C
C^L
C
C     ==================================
      REAL FUNCTION LSTLSQ(MINDX, IH, IK, IL)
C     ==================================
C
      IMPLICIT NONE
C     
      EXTERNAL LSTLSQ1
C
C     PARAMETERS
      INTEGER MINDX, IH, IK, IL
C
C     reso RETURN value
      REAL RESO
C
      CALL LSTLSQ1(RESO, MINDX, IH, IK, IL)
      LSTLSQ=RESO
      RETURN
      END
C     
C^L   
C     
C     ==================================
      REAL FUNCTION STHLSQ(IH, IK, IL)
C     ==================================
C     
      IMPLICIT NONE
C     
      EXTERNAL STHLSQ1
C
C     PARAMETERS
      INTEGER IH, IK, IL
C
C     reso RETURN value
      REAL RESO
C
      CALL STHLSQ1(RESO, IH, IK, IL)
      STHLSQ=RESO
      RETURN
      END
C     
C^L   
C     
C     ==================================
      REAL FUNCTION STS3R4(IH, IK, IL)
C     ==================================
C     
      IMPLICIT NONE
C
      EXTERNAL STS3R41
C     
C     PARAMETERS
      INTEGER IH, IK, IL
C
C     reso RETURN value
      REAL RESO
C
      CALL STS3R41(RESO, IH, IK, IL)
      STS3R4=RESO
      RETURN
      END

