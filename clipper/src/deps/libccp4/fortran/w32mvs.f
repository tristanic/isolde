C
C     w32mvs.f: platform-dependent low-level functions for MVS
C     Copyright (C) 1999  Alun Ashton
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
C  *** this file is an equivelant to unix.f or vms.for but 
C  *** for the microsoft visual studio packages.
C
C ========
C w32mvs.f
C ========
C
C CCPOPN - open a file
C UBYTES - Returns number of bytes per word and 'words'/'bytes'
C          to indicate if byte handling is available
C UCPUTM - Get CPU time
C UGERR  - Get error explanation
C UGTENV - Get value of env. variable
C UGTUID - Get user id - it's name
C UIDATE - Get date in 3 integer 
C UISATT - Is file a terminal?
C USTIME - Get absolute time in seconds (-1 for VMS)
C UTIME  - Get current time
C VAXVMS - Logical function returns TRUE if VAX/VMS
C WINMVS - Logical function returns TRUE if W32/Microsoft dev studio
C TTSEND - Write string to terminal with various carriage control
C     options
C GETELAPSED - Print timing info for CCPERR - removed for now
C UGTARG - Get command-line argument
C GETREF - Abstracted from abscale since it has BYTE declaration.
C CCPSPW - Spawns a new process to run shell command
C RTNBKS - Returns a Backslash as unix compilers are fussy!
C
C
      SUBROUTINE CCPOPN(IIUN,LOGNAM,KSTAT,ITYPE,LREC,IFAIL)
C     ====================================================
C
C---- This subroutine is used to open a file
C
C     The requirement to specify that leading carriage control
C     characters in the output records should be obeyed (or not) can't
C     be implemented portably; likewise specifying readonly opening.
C     Some compilers accept VAXtran `carriagecontrol=' and `readonly'
C     specifiers; if so we use them.  Others have IOINIT, which can be
C     used to specify the carriage control.  The HPUX compiler is said
C     not to have any means of doing this and AIX seems to be likewise,
C     sigh; they both seem to obey the normal unix convention of
C     printing the  as-is rather than obeying the first character
C     as carriage control.  Concentrix does obey the first column a la
C     VMS and `traditional' Fortran; the MIPS compilers have a compile
C     (link?) option to do so.  Unfortunately, carriagecontrol
C     specification isn't even defined in Fortan90, although
C     `ACTION="READ"' can be used.
C
C PARAMETERS
C ==========
C
C        IIUN (I)   UNIT NUMBER
C      LOGNAM (I)   LOGICAL FILE NAME (UP TO 8 CHARACTERS)
C       KSTAT (I)   FILE STATUS FLAG =1, 'UNKNOWN'
C                                    =2, 'SCRATCH'
C                                    =3, 'OLD'
C                                    =4, 'NEW'
C                                    =5, 'READONLY'
C                                    =6, 'PRINTER'
C       ITYPE (I)   FILE TYPE FLAG =1, 'SEQUENTIAL' 'FORMATTED'
C                                  =2, 'SEQUENTIAL' 'UNFORMATTED'
C                                  =3, 'DIRECT'     'FORMATTED'
C                                  =4, 'DIRECT'     'UNFORMATTED'
C        LREC (I)   RECORD LENGTH FOR DIRECT ACCESS FILE (NO. OF
C                   CHARACTERS FOR A FORMATTED FILE OR WORDS FOR
C                   AN UNFORMATTED FILE). NOT RELEVANT FOR A SEQUENTIAL
C                   FILE
C       IFAIL (I/O) ON INPUT:     =0, STOP ON OPEN FAILURE
C                                 =1, CONTINUE AFTER OPEN FAILURE
C                                 =2, CONTINUE SILENTLY AFTER OPEN FAILURE
C                   ON OUTPUT:    UNCHANGED IF FILE OPEN OK
C                                 =-1, ERROR IN OPENING FILE
C
C     .. Scalar Arguments ..
      INTEGER IFAIL,KSTAT,ITYPE,IIUN,LREC
      CHARACTER LOGNAM* (*)
C     ..
C     .. Local Scalars ..
      INTEGER LLREC,IUN,IBYTES,ISTAT,L,IOS
      CHARACTER CCNTRL*7,ST*7,FRM*12,ERRSTR*500,
     +     NAMFIL*255,HANDLE*5,OPNVAR*20, access*10
      INTEGER UNKNWN, SCRTCH, OLD, NEW, RDONLY, PRINTR
      PARAMETER (UNKNWN=1, SCRTCH=2, OLD=3, NEW=4, RDONLY=5, PRINTR=6)
      LOGICAL CCPEXS, LNONAM
C     ..
C     .. Local Arrays ..
      CHARACTER STAT(6)*7, DISP*6
C     ..
C     .. External Functions ..
      INTEGER LENSTR, LUNSTO
      EXTERNAL LENSTR, LUNSTO
C     ..
C     .. External Subroutines ..
      EXTERNAL UGERR,UGTENV, CCPEXS
C     ..
C     .. Data statements ..
C     NB mustn't have SCRATCH in here, because result is system
C     -dependent

c     but as this is a system dependant file its ok!
      DATA STAT/'UNKNOWN','SCRATCH','OLD','NEW','OLD','UNKNOWN'/
C     ..
C     
      ISTAT = KSTAT
C     Negative unit number means don't give messages for successful open
      IUN = IIUN
      IF (IIUN.LT.0) IUN = -IIUN
C     Check args:
      IF (ISTAT.LT.1 .OR. ISTAT.GT.6 .OR. ITYPE.LT.1 .OR. ITYPE.GT.4)
     +     THEN 
        IF (IFAIL.EQ.0) THEN
          CALL CCPERR(1,
     +         '**CCPOPN ERROR** Invalid parameters in call')
        ELSE
          WRITE (LUNSTO(1),
     +         '('' **CCPOPN ERROR** Invalid parameters in call'',/)')
          IFAIL = -1
        END IF
        RETURN
      ENDIF 
C
C     Do nothing for pre-connected units (what's the significance of
C     `TERM...'?) 
      IF (LOGNAM.EQ.'DATA' .OR. LOGNAM.EQ.'PRINTER' .OR.
     +     LOGNAM(:4).EQ.'TERM') RETURN
C
C     if environment variable CCP4_OPEN has value `UNKNOWN', open files
C     with status UNKNOWN rather than new if they exist
      IF (ISTAT.EQ.NEW) THEN
        OPNVAR = ' '
        CALL UGTENV('CCP4_OPEN',OPNVAR)
        IF (OPNVAR.EQ.'UNKNOWN' .OR. OPNVAR.EQ.'unknown') ISTAT = 1
      END IF
C
C     check for `logical name' referencing real file
      CALL UGTENV(LOGNAM,NAMFIL)
      LNONAM = .FALSE.
      IF (NAMFIL.EQ.' ') THEN
        IF (.NOT. CCPEXS(LOGNAM)) LNONAM = .TRUE.
        NAMFIL = LOGNAM
      END IF
C     VMS null device (VMS code canonicalises /dev/null)
      IF (NAMFIL.EQ.'NL:' .OR. NAMFIL.EQ.'nl:') NAMFIL='/dev/null'
C     Special case:  /dev/null should be opened UNKNOWN
      IF ( NAMFIL.EQ.'/dev/null') ISTAT = 1
C
C     type of open
      ST = STAT(ISTAT)
      IF (ITYPE.EQ.2 .OR. ITYPE.EQ.4) THEN
        FRM = 'UNFORMATTED'
      ELSE
        FRM = 'FORMATTED'
      ENDIF 
      IF (ITYPE .EQ. 1 .OR. ITYPE.EQ.2) THEN
        ACCESS='SEQUENTIAL'
      ELSE
        ACCESS='DIRECT'
      ENDIF
C
      IF (ISTAT.EQ.SCRTCH) THEN
        DISP = 'DELETE'
      ELSE
        DISP = 'KEEP'
      ENDIF
C     
      IF (access.eq.'DIRECT') THEN
C       Need to check is record length in words or bytes and set LLREC
C       accordingly. 
        CALL UBYTES (IBYTES,HANDLE)
        LLREC = LREC*IBYTES
        IF (HANDLE.EQ.'WORDS'.AND.ITYPE.EQ.4) LLREC=LLREC/IBYTES
        IF (ISTAT.EQ.RDONLY) THEN
C          may be defined as null or as `READONLY,'
          OPEN(UNIT=IUN,STATUS='UNKNOWN',ACCESS='DIRECT',FORM=FRM,
     +         READONLY,
     +         FILE=NAMFIL,RECL=LLREC,IOSTAT=IOS,ERR=5)
        ELSE
          OPEN(UNIT=IUN,STATUS='UNKNOWN',ACCESS='DIRECT',FORM=FRM,
     +         DISPOSE=DISP,
     +         FILE=NAMFIL,RECL=LLREC,IOSTAT=IOS,ERR=5)
        ENDIF 
      ELSE
C       if available, carriagecontrol='fortran' for print file, else = 
C       'list'.  we can use ioinit instead where it's available (see e.g.
C       Sun manual). 
        IF (ISTAT.EQ.PRINTR) THEN
C         want to obey  characters in column 1
          CCNTRL = 'FORTRAN'
          FRM = 'FORMATTED'
        ELSE
C         no special significance to column 1
          CCNTRL = 'LIST'
        END IF
        IF (FRM .EQ. 'UNFORMATTED') THEN
C         (carriage control not relevant)
          IF (ISTAT.EQ.RDONLY) THEN
            OPEN(UNIT=IUN, FILE=NAMFIL, STATUS=ST, ACCESS='SEQUENTIAL',
     +           READONLY,
     +           FORM=FRM, ERR=5, IOSTAT=IOS)
          ELSE
            OPEN(UNIT=IUN, FILE=NAMFIL, STATUS=ST, ACCESS='SEQUENTIAL',
     +           DISPOSE=DISP,
     +           FORM=FRM, ERR=5, IOSTAT=IOS)
          ENDIF
        ELSE
          IF (ISTAT.EQ.RDONLY) THEN
            OPEN(UNIT=IUN, FILE=NAMFIL, STATUS=ST, ACCESS='SEQUENTIAL',
     +           READONLY,
     +           CARRIAGECONTROL=CCNTRL,
     +           FORM=FRM, ERR=5, IOSTAT=IOS)
          ELSE
            OPEN(UNIT=IUN, FILE=NAMFIL, STATUS=ST, ACCESS='SEQUENTIAL',
     +           CARRIAGECONTROL=CCNTRL,
     +           DISPOSE=DISP,
     +           FORM=FRM, ERR=5, IOSTAT=IOS)
          ENDIF
        ENDIF
      ENDIF
C
C     Scratch files are immediately unlinked from the directory; they
C     become inaccessible only when closed, but don't appear in the
C     directory and the name can be re-used.
C     NB this may break with REWIND if that is implemented as close +
C     reopen, sigh.  See also  above
c  *** this now removed for windows. we will try actually opeing scratch files instead
c	IF (ISTAT.EQ.SCRTCH) CALL CUNLINK(NAMFIL)
C
C     Error check
 5    CONTINUE
C     don't report UNKNOWN if actually SCRATCH
      IF (ISTAT.EQ.SCRTCH) ST = 'SCRATCH'
      IF (IOS.NE.0) THEN
        CALL UGERR(IOS,ERRSTR)
        IF (IFAIL.EQ.0) THEN
C         warning if there was no file associated with logical name
          IF (LNONAM) THEN
             ERRSTR = 'CCPOPN Logical name '//LOGNAM
             ERRSTR(LENSTR(ERRSTR)+2:) = 'has no associated file name'
             CALL CCPERR(2,ERRSTR)
          END IF
C         hard failure
          WRITE (LUNSTO (1),FMT=6002) IUN, NAMFIL(1:LENSTR(NAMFIL)),
     +         LOGNAM(1:LENSTR(LOGNAM))
 6002     FORMAT (' Open failed: Unit:',I4,', File: ',A, ' (logical: ',
     +         A, ')')
          ERRSTR = ' Open failed: File: ' // NAMFIL
          CALL CCPERR(-1, ERRSTR)
        else
C         soft failure
          IF (IFAIL.EQ.1) WRITE (lunsto (1),FMT=6004) FRM, ST, IUN, 
     +         LOGNAM(1:LENSTR(LOGNAM)), NAMFIL(1:LENSTR(NAMFIL)),
     +         ERRSTR(1:LENSTR(ERRSTR))
 6004     FORMAT (' **CCPOPN ERROR**  ',A,3X,A,
     +         ' file open failure on unit ',I3,/' Logical name: ',
     +         A,', ','File name: ',A/1X,A/)
          IFAIL = -1
          RETURN            
        ENDIF
      ELSE
        IF (IIUN.LE.0) RETURN 
        WRITE (ERRSTR,FMT=6000) FRM,ST,IUN
        CALL QPRINT (1, ' ')
        CALL QPRINT (1, ERRSTR)
        call ccp4h_summary_beg()
        ERRSTR = 'Logical name: '
        ERRSTR (15:) = LOGNAM
        L = MIN(LENSTR (ERRSTR) + 1, LEN (ERRSTR))
        ERRSTR (L:) = ', Filename: ' // NAMFIL
        CALL QPRINT (1, ERRSTR)
        call ccp4h_summary_end()
        CALL QPRINT (1, ' ')
 6000 FORMAT (A,3X,A,' file opened on unit ',I3)
      ENDIF 
      END
C
C
C     ==============================
      SUBROUTINE UBYTES(INUM,STRING)
C     ==============================
C
C UBYTES - Return statistics about byte handling
C
C Input:  none
C
C Output:    INUM - number of bytes per word
C            HANDLE - 'WORDS' or 'BYTES'
C            HANDLE - For unformatted files records are usually
C                     counted in 'BYTES', however both VAX and 
C                     SGI swap to 'WORDS' for this file type.
C
C Arguments: INTEGER     INUM
C            CHARACTER*5 HANDLE
C
C Usage:     CALL UBYTES (INUM,HANDLE)
C
C     .. Scalar Arguments ..
      INTEGER INUM
      CHARACTER STRING*5
C     ..
C
C
      INUM = 4
      STRING = 'BYTES'
C
      END
C
C
C     ======================
      SUBROUTINE UCPUTM(SEC)
C     ======================
C
C     Get CPU time in seconds
C
C     Parameter:
C     REAL SEC (i/o): If sec<=0.0, initialize timer and return current
C                     elapsed cpu time since start of execution, otherwise
C                     return elapsed cpu since timer was initialized.
C                     Time is in seconds.
C
C     .. Scalar Arguments ..
      REAL SEC
C     ..
C     .. Local Scalars ..
      REAL TLAST
C     ..
C     .. Local Arrays ..
      REAL TARRAY(2)
C     ..
C     .. Save statement ..
      SAVE TLAST
C     ..
      IF (SEC.LE.0.0) THEN
        TLAST = ETIME (TARRAY)
        SEC = TLAST
      ELSE
        SEC = ETIME (TARRAY) - TLAST
      ENDIF
      
      END
C
C
C     ===============================
      SUBROUTINE UGERR(STATUS,ERRSTR)
C     ===============================
cDEC$ IF DEFINED (__INTEL_COMPILER) 
      USE IFCORE 
cDEC$ ENDIF
C UGERR - Get error message string for error number in STATUS
C     (supposedly).  Actually it ignores STATUS and always uses the
C     *last* error that occurred.
C
C Input:     STATUS - Error number (if negative print error message)
C
C Output:    ERRSTR - Error message string
C
C Arguments: INTEGER       STATUS
C            CHARACTER*(*) ERRSTR
C
C Usage:     CALL UGERR(STATUS, ERRSTR)
C
C     .. Scalar Arguments ..
      INTEGER STATUS
      CHARACTER ERRSTR* (*)
C     ..
C     .. Local Scalars ..
      LOGICAL IPRINT
C     ..
C     .. External Subroutines ..
      INTEGER IERRNO, LUNSTO
      EXTERNAL IERRNO, LUNSTO
C     ..
      IPRINT = .FALSE.
      IF (STATUS.LT.0) THEN
        IPRINT = .TRUE.
        STATUS = -STATUS
      END IF
C
C---- Get error message from system
C
      IF (IERRNO().NE.0) THEN
        CALL GERROR(ERRSTR)
      ELSE
        ERRSTR = ' '
      ENDIF
      IF (IPRINT) WRITE (LUNSTO(1),FMT=6000) 'UGERR',ERRSTR
C
 6000 FORMAT (' ',A,': ',A)
      END
C
C     ================================
      SUBROUTINE UGTENV(NAMENV,VALENV)
C     ================================
C
C UGTENV - Get value of env. variable
C
C Input:     NAMENV - Logical Name (trailing blanks are stripped)
C
C Output:    VALENV - Its value
C
C Arguments: CHARACTER*(*) NAMENV, VALENV
C
C Usage:     CALL UGTENV(NAMENV, VALENV)
C
C     .. Scalar Arguments ..
      CHARACTER NAMENV* (*),VALENV* (*)
C     ..
C     .. External Subroutines ..
C     don't declare getenv
      INTEGER LENSTR
      EXTERNAL LENSTR
C     ..
      CALL GETENV(NAMENV(:LENSTR(NAMENV)),VALENV)
C
      END
C
C
C     =========================
      SUBROUTINE UGTUID(USRNAM)
C     =========================
C
C UGTUID - Get user ID
C
C Customised for win nt
C
C Input:     none
C
C Output:    UID - user ID string
C
C Arguments: CHARACTER*(*) UID
C
C Usage:     CALL UGTUID(UID)
C
C     .. Scalar Arguments ..
      CHARACTER USRNAM* (*)
C     ..
C     .. External Subroutines ..
C     don't declare getenv
C     ..
C      CALL GETENV('USER',USRNAM)
      CALL GETLOG(USRNAM)
      IF (USRNAM.EQ.' ') CALL GETENV('LOGNAME',USRNAM)
C
      END
C
C
C     ====================================
      SUBROUTINE UIDATE(IMONTH,IDAY,IYEAR)
C     ====================================
C
C UIDATE - Get date in 3 integer . Alliant uses INTEGER*4
C          and order is IDAY,IMONTH,IYEAR
C
C Input:     none
C
C Output:    MONTH,DAY,YEAR
C
C Arguments: INTEGER MONTH, DAY, YEAR
C
C Usage:     CALL UIDATE(MONTH, DAY, YEAR)
C
C     .. Scalar Arguments ..
      INTEGER IDAY,IMONTH,IYEAR
C     .. Local Arrays ..
      INTEGER IARRAY(3)
C     ..
C     don't declare IDATE -- it's often an intrinsic to avoid confusion
C     between the two possible calling sequences.  On SGI, it's only
C     documented in one style but seems to work with the canonical
C     Unix one too; however, the order of the arguments of the
C     documented version (used here) *isn't* the same as for VMS...
C
      CALL IDATE (IARRAY(1), IARRAY(2), IARRAY(3))
      IDAY = IARRAY(1)
      IMONTH = IARRAY(2)
      IYEAR = IARRAY(3)
C
      END
C
C
C     ==============================
      SUBROUTINE UISATT(FLUN,ANSWER)
C     ==============================
C
C UISATT - This function determines whether a program is being
C          run on-line if this information is available.
C
C Input:     FLUN - Fortran Unit Number
C
C Output:    ANS - 1 for on-line, 0 otherwise
C
C Arguments: INTEGER FLUN, ANS
C
C Usage:     CALL UISATT (FLUN,ANS)
C
C     .. Scalar Arguments ..
      INTEGER ANSWER,FLUN
C     ..
      LOGICAL ISATTY
      EXTERNAL ISATTY
      ANSWER = 0
      IF (ISATTY(FLUN)) ANSWER = 1
C
      END
C
C
C     =======================
      SUBROUTINE USTIME(ISEC)
C     =======================
C
C USTIME - Get absolute time in seconds.
C          Convex uses STIME (), others seem to use TIME ().
C
C Input:     none
C
C Output:    SEC
C
C Arguments: INTEGER SEC
C
C Usage:     CALL USTIME(SEC)
C
      INTEGER ISEC
C
      INTEGER TIME
C
      ISEC = TIME()
C
      END
C
C
C     =======================
      SUBROUTINE UTIME(CTIME)
C     =======================
C
C UTIME - Get current time hh:mm:ss
C
C Input:     none
C
C Output:    TIME - as ASCII string
C
C Arguments: CHARACTER*(*) CTIME
C
C Usage:     CALL UTIME(CTIME)
C
C     .. Scalar Arguments ..
      CHARACTER CTIME* (*)
C     ..
C     .. Local Arrays ..
      INTEGER IARRAY(3)
C     ..
      CALL ITIME(IARRAY)
      WRITE (CTIME,FMT=6000) IARRAY(1),IARRAY(2),IARRAY(3)
 6000 FORMAT (I2,2 (':',I2.2))
      END
C
C
C     =========================
      LOGICAL FUNCTION VAXVMS()
C     =========================
C
C VAXVMS - Operating Sytem in use returns .TRUE. if VAXVMS
C
C Input:     none
C
C Returns:   .TRUE. for VAXVMS, .FALSE. otherwise
C
C Arguments: none
C
C Usage:     VAXVMS ()
C
      VAXVMS = .FALSE.
C
      END
C
C
C     =========================
      LOGICAL FUNCTION WINMVS()
C     =========================
C
C WINMVS - Windows mircrosoft Visual Studio
C
C Input:     none
C
C Returns:   .TRUE. for WINMVS, .FALSE. otherwise
C
C Arguments: none
C
C Usage:     WINMVS ()
C
      WINMVS = .TRUE.
C
      END
C
C
C SUBROUTINE 'TTSEND'
C ===================
C
C Write a string to a terminal with various carriage control options
C for LAUE
C
      SUBROUTINE TTSEND (IUN, STR, ICC)
C
C Parameters:
C
C         IUN (I)   Unit number for the output
C         STR (I)   The string to be output
C         ICC (I)   = 0, no carriage control at the end of the string
C                        (for prompts)
C                        e.g. for routine TPROMP
C                   = 1, normal carriage control
C                        e.g. for routine TWRITE
C                   = 2, no carriage control (for sending escape/control
C                        character sequences to ANSI/T4014 terminals)
C                        e.g. for QSCREEN graphics routines
C                   = 3, Output line at current point on screen (no leading
C                        line feed or carriage return - trailing does not
C                        matter)
C
C Machine dependence examples: Convex   1000  FORMAT (A,$)
C                                       1001  FORMAT (A)
C                                       1002  FORMAT (A,$)
C                                       1003  FORMAT (A)
C                              
C                              Vax      1000  FORMAT (' ',A,$)
C                                       1001  FORMAT (' ',A)
C                                       1002  FORMAT ('+',A,$)
C                                       1003  FORMAT ('+',A)
C
C
C====== Specification statements
C
      CHARACTER*(*) STR
      CHARACTER*10 CCNTRL
C
C====== Write string
C
C     'LIST' is the equivalent of the normal unix state
      CCNTRL = 'LIST'
C     in the case of systems obeying the carriagecontrol specifier, 
C     we assume the stream has actually been opened, so that the
C     specifier is suitably defined -- on the Alliant, for instance,
C     it will be 'UNKNOWN' for an unopened stream (6 is pre-opened)
C
      IF (CCNTRL .EQ. 'FORTRAN') THEN
C       VMS-type
        IF (ICC.EQ.0) THEN
          WRITE (IUN,1004) STR
        ELSE IF (ICC.EQ.2) THEN
          WRITE (IUN,1006) STR
        ELSE IF (ICC.EQ.3) THEN
          WRITE (IUN,1007) STR
        ELSE
          WRITE (IUN,1005) STR
        ENDIF
      ELSE
        IF (ICC.EQ.0) THEN
          WRITE (IUN,1000) STR
        ELSE IF (ICC.EQ.2) THEN
          WRITE (IUN,1002) STR
        ELSE IF (ICC.EQ.3) THEN
          WRITE (IUN,1003) STR
        ELSE
          WRITE (IUN,1001) STR
        ENDIF
      ENDIF
C     these formats are mostly non-standard, of course...
1000  FORMAT (A,$)
1001  FORMAT (A)
1002  FORMAT (A,$)
1003  FORMAT (A)
 1004 FORMAT (' ',A,$)
 1005 FORMAT (' ',A)
 1006 FORMAT ('+',A,$)
 1007 FORMAT ('+',A)
      END
C
C
C     =====================
c      SUBROUTINE GETELAPSED
cC     =====================
cC
c      EXTERNAL LUNSTO, USTIME
c      INTEGER LUNSTO
c      REAL TARRAY(2), JUNK
c      INTEGER ELAPS, START
c      LOGICAL INITED
c      SAVE START, INITED
c      DATA INITED /.FALSE./
cC     
c      JUNK = ETIME(TARRAY)
c      CALL USTIME(ELAPS)
c      ELAPS = ELAPS - START
cC     don't print anything if it hasn't been initialised (by CCPFYP)
c      IF (INITED) WRITE(LUNSTO(1),6000) TARRAY(1), TARRAY(2), 
c     +     ELAPS/60, MOD(ELAPS, 60)
c 6000 FORMAT(' Times: User: ', F9.1, 's System: ', F6.1, 's Elapsed:',
c     +     I5 , ':',I2.2)
cC     
c      ENTRY INITFYP
c      CALL USTIME(START)
c      INITED = .TRUE.
cC     Machine-dependent startup, e.g. set FPE on SunOS

c      END
C
      SUBROUTINE UGTARG(I, ARG)
      INTEGER I
      CHARACTER *(*) ARG
      CALL GETARG(I, ARG)
      END
C     
C     =====================================================
      SUBROUTINE GETREF(KERFLG,NREAD,NSPOTS,DYNAM,MAXOMITL)
C     =====================================================
C
C     This has been abtracted from ABSCALE because of the BYTE
C     declaration.
C
C        implicit none
C     
C     
C     
C     
C     
C     Read one reflection into common /IND/, skipping unmeasured reflections
C     Return 1 if end of file or all N spots found
C     Both integrated and profile fitted I's and SD's are stored, one in
C     INTT,SD and the other in INTT2,SD2. The values in INTT,SD are used
C     in scaling, and this is chosen on input card 1 to be either the 
C     integrated or profile fitted value.
C     
C
C This routine is probably VAX specific in its unpacking of indices
C
C
C
C---- IC  generate file variables
C
C
C
C     .. Scalar Arguments ..
      INTEGER           NREAD,NSPOTS,KERFLG,MAXOMITL
      LOGICAL DYNAM
C     ..
C     .. Scalars in Common ..
      INTEGER           IREC,IX,IY,JGUNIT,JH,JK,JL,MND
      LOGICAL           PROFILE
C     ..
C     .. Arrays in Common ..
      REAL              SPACER(12)
      INTEGER           INTT(3),INTT2(3),ISD(3),ISD2(3),JUNK(2)
C     ..
C     .. Local Scalars ..
      INTEGER           I,ICOL,ICOL2,IER,I4INTS,I4INTP
      BYTE              IR,IM
C     ..
C     .. Local Arrays ..
cejd      INTEGER*2         IBUF(18)
      INTEGER*2         IBUF(19)
      BYTE              B(2)
C     ..
C     .. External Subroutines ..
      EXTERNAL          QREAD
C     ..
C     .. Common blocks ..
       LOGICAL BRIEF
       INTEGER IBRIEF
       COMMON /BRF/ BRIEF,IBRIEF
      COMMON      /IND/JH,JK,JL,MND,JUNK,IX,IY,SPACER,INTT,ISD,
     +            INTT2,ISD2
      COMMON      /INREC/JGUNIT,IREC
      COMMON      /INTTYP/PROFILE
C     ..
C     .. Equivalences ..
      EQUIVALENCE       (B(1),IBUF(4)), (B(1),IR), (B(2),IM)
      EQUIVALENCE       (I4INTS,IBUF(7)),(I4INTP,IBUF(13))
C     ..
      SAVE
C
C
          KERFLG = 0
C
C
   10 CONTINUE
      NREAD = NREAD + 1
C
C
      IF (NREAD.GT.NSPOTS) THEN
          GO TO 40
      ELSE
C
C              *************************
          CALL QREAD(JGUNIT,IBUF,36,IER)
C              *************************
C
          IREC = IREC + 1
          IF (IER.NE.0) THEN
              GO TO 30
C
C---- If rejected, skip to next refl
C
CAL ALLOW IR TO HAVE VALUES 5,6
          ELSE IF ((IR.NE.0).AND.(IR.LE.4)) THEN
              GO TO 10
          END IF
      END IF
C
C
      JH = IBUF(1)
      JK = IBUF(2)
      JL = IBUF(3)
      MND = IM
      IF (MND.LT.0) MND = 8
      IX = IBUF(5)
      IY = IBUF(6)
C
C---- A film intensity in ibuf(7) for integrated intensities or
C     ibuf(13) for profile fitted intensities
C
      IF (PROFILE) THEN
          ICOL = 13
          ICOL2 = 7
      ELSE
          ICOL = 7
          ICOL2 = 13
      END IF
C
C
      DO 20 I = 1,3
          IF (DYNAM) THEN
           ISD(I) = IBUF(ICOL+2)
           ISD2(I) = IBUF(ICOL2+2)
           IF (PROFILE) THEN
             INTT(I) = I4INTP
             INTT2(I) = I4INTS
           ELSE
             INTT(I) = I4INTS
             INTT2(I) = I4INTP
           END IF
          ELSE
           INTT(I) = IBUF(ICOL)
           ISD(I) = IBUF(ICOL+1)
           INTT2(I) = IBUF(ICOL2)
           ISD2(I) = IBUF(ICOL2+1)
          END IF
C
C---- Test for badspots (isd=-9999) change to unmeasured
C     this will also reject overloaded reflections
C-AL   Change this so overloads are rejected (and counted) in RDREF
C
       IF ( (ISD(I)   .EQ. -9999) .AND.
     +      (INTT(I)  .NE. MAXOMITL) )       INTT(I) = -9999
       IF ( (ISD2(I)  .EQ. -9999) .AND.
     +      (INTT2(I) .NE. MAXOMITL) ) 
     +                                     INTT2(I) = -9999
C
C
          ICOL = ICOL + 2
          ICOL2 = ICOL2 + 2
   20     CONTINUE
      RETURN
   30 KERFLG = -1
      RETURN 
   40 KERFLG = -1
      RETURN
C
C
      END
C_BEGIN_CCPSPW
      SUBROUTINE CCPSPW(STRING)
C     =========================
C
C     Spawns a new process to run shell command
C
C Arguments:
C ==========
C
C  STRING (I)   CHARACTER*(*): string containing command
C_END_CCPSPW
C
       CHARACTER STRING*(*)
       EXTERNAL SYSTEM
       CALL SYSTEM(STRING)
       END
C
      SUBROUTINE CCPAL1 (ROUTNE, N, TYPE, LENGTH)
C
C     Arrange to call ROUTNE with N TYPEd array arguments of given
C     LENGTH (see CCPALC)
C
      EXTERNAL ROUTNE
      INTEGER N, TYPE (*), LENGTH (*)
      INTEGER I, SIZES (5), POINTER (12), ISTAT
C     bytes per word (assuming 32 bit words...)
      DATA SIZES /4,4,8,8,1/

C     The calling routine, CCPALC, will have checked that the arguments
C     are in range
      DO I=1,N
        POINTER(I) = malloc(SIZES(TYPE (I))*LENGTH(I))
        IF (POINTER(I) .eq. 0)
     +       CALL CCPERR (-1, 'CCPALC: can''t allocate memory')
        CALL CCPZBI (%VAL(POINTER(I)), SIZES(TYPE (I))*LENGTH(I))
      ENDDO
      IF (N.EQ.1) THEN
        CALL ROUTNE (LENGTH (1), %VAL(POINTER(1)))
      ELSE IF (N.EQ.2) THEN
        CALL ROUTNE (
     +       LENGTH (1), %VAL(POINTER(1)), LENGTH (2), %VAL(POINTER(2)))
      ELSE IF (N.EQ.3) THEN
        CALL ROUTNE (
     +       LENGTH (1), %VAL(POINTER(1)), LENGTH (2), %VAL(POINTER(2)),
     +       LENGTH (3), %VAL(POINTER(3)))
      ELSE IF (N.EQ.4) THEN
        CALL ROUTNE (
     +       LENGTH (1), %VAL(POINTER(1)), LENGTH (2), %VAL(POINTER(2)),
     +       LENGTH (3), %VAL(POINTER(3)), LENGTH (4), %VAL(POINTER(4)))
      ELSE IF (N.EQ.5) THEN
        CALL ROUTNE (
     +       LENGTH (1), %VAL(POINTER(1)), LENGTH (2), %VAL(POINTER(2)),
     +       LENGTH (3), %VAL(POINTER(3)), LENGTH (4), %VAL(POINTER(4)),
     +       LENGTH (5), %VAL(POINTER(5)))
      ELSE IF (N.EQ.6) THEN
        CALL ROUTNE (
     +       LENGTH (1), %VAL(POINTER(1)), LENGTH (2), %VAL(POINTER(2)),
     +       LENGTH (3), %VAL(POINTER(3)), LENGTH (4), %VAL(POINTER(4)),
     +       LENGTH (5), %VAL(POINTER(5)), LENGTH (6), %VAL(POINTER(6)))
      ELSE IF (N.EQ.7) THEN
        CALL ROUTNE (
     +       LENGTH (1), %VAL(POINTER(1)), LENGTH (2), %VAL(POINTER(2)),
     +       LENGTH (3), %VAL(POINTER(3)), LENGTH (4), %VAL(POINTER(4)),
     +       LENGTH (5), %VAL(POINTER(5)), LENGTH (6), %VAL(POINTER(6)),
     +       LENGTH (7), %VAL(POINTER(7)))
      ELSE IF (N.EQ.8) THEN
        CALL ROUTNE (
     +       LENGTH (1), %VAL(POINTER(1)), LENGTH (2), %VAL(POINTER(2)),
     +       LENGTH (3), %VAL(POINTER(3)), LENGTH (4), %VAL(POINTER(4)),
     +       LENGTH (5), %VAL(POINTER(5)), LENGTH (6), %VAL(POINTER(6)),
     +       LENGTH (7), %VAL(POINTER(7)), LENGTH (8), %VAL(POINTER(8)))
      ELSE IF (N.EQ.9) THEN
        CALL ROUTNE (
     +       LENGTH (1), %VAL(POINTER(1)), LENGTH (2), %VAL(POINTER(2)),
     +       LENGTH (3), %VAL(POINTER(3)), LENGTH (4), %VAL(POINTER(4)),
     +       LENGTH (5), %VAL(POINTER(5)), LENGTH (6), %VAL(POINTER(6)),
     +       LENGTH (7), %VAL(POINTER(7)), LENGTH (8), %VAL(POINTER(8)),
     +       LENGTH (9), %VAL(POINTER(9)))
      ELSE IF (N.EQ.10) THEN
        CALL ROUTNE (
     +       LENGTH (1), %VAL(POINTER(1)), LENGTH (2), %VAL(POINTER(2)),
     +       LENGTH (3), %VAL(POINTER(3)), LENGTH (4), %VAL(POINTER(4)),
     +       LENGTH (5), %VAL(POINTER(5)), LENGTH (6), %VAL(POINTER(6)),
     +       LENGTH (7), %VAL(POINTER(7)), LENGTH (8), %VAL(POINTER(8)),
     +       LENGTH (9), %VAL(POINTER(9)),
     +       LENGTH (10), %VAL(POINTER(10)))
      ELSE IF (N.EQ.11) THEN
        CALL ROUTNE (
     +       LENGTH (1), %VAL(POINTER(1)), LENGTH (2), %VAL(POINTER(2)),
     +       LENGTH (3), %VAL(POINTER(3)), LENGTH (4), %VAL(POINTER(4)),
     +       LENGTH (5), %VAL(POINTER(5)), LENGTH (6), %VAL(POINTER(6)),
     +       LENGTH (7), %VAL(POINTER(7)), LENGTH (8), %VAL(POINTER(8)),
     +       LENGTH (9), %VAL(POINTER(9)),
     +       LENGTH (10), %VAL(POINTER(10)),
     +       LENGTH (11), %VAL(POINTER(11)))
      ELSE IF (N.EQ.12) THEN
        CALL ROUTNE (
     +       LENGTH (1), %VAL(POINTER(1)), LENGTH (2), %VAL(POINTER(2)),
     +       LENGTH (3), %VAL(POINTER(3)), LENGTH (4), %VAL(POINTER(4)),
     +       LENGTH (5), %VAL(POINTER(5)), LENGTH (6), %VAL(POINTER(6)),
     +       LENGTH (7), %VAL(POINTER(7)), LENGTH (8), %VAL(POINTER(8)),
     +       LENGTH (9), %VAL(POINTER(9)),
     +       LENGTH (10), %VAL(POINTER(10)),
     +       LENGTH (11), %VAL(POINTER(11)),
     +       LENGTH (12), %VAL(POINTER(12)))
      ENDIF
      DO I=1,N
        CALL free(POINTER(I))
      ENDDO
      END
C
      SUBROUTINE CEXIT (ICODE)
C     trivial interface to system-dependent EXIT routine
      INTEGER ICODE
      CALL EXIT (ICODE)
      END
 
      subroutine clear
      end




        subroutine gdummy
        character *(*) char_dummy
        entry  qreset
          return
        entry  reshap
          return
        entry  qdevic(keybd)
          return
        entry  winope(char_dummy,i0)
          return
        entry  keepas(i1,i2)
          return
        entry  draw2i(i3,i4)
          return
        entry  move2i(i5,i6)
          return
        entry  loadma(i7)
          return
        entry  gconfi
          return
        entry  mmode(i8)
          return
        entry  foregr
          return
        entry  getval(i9)
          return
        entry  color(i10)
          return
        entry  getsiz(r1,r2)
          return
cc        entry  clear this is in somewhere else in -ltermcap
cc          return
        entry  ortho2(r3,r4,r5,r6)
          return
        entry  getori(r7,r8)
          return
        end
C
C
C     =========================
       CHARACTER FUNCTION RTNBKS()
C     =========================
C
C RTNBKS - Returns a Backslash for nt as unix compilers are fussy!
C
C Input:     none
C
C Returns:   \ if WIN32 or space if unix or vms
C
C Arguments: none
C
C Usage:     RTNBKS ()
C
      RTNBKS='\'
C
      END

c     ============================
      subroutine hciftime(ciftime)
c     ============================

c     Uses f90 intrinsic Date_and_Time. Using f77:
c       works on VMS Fortran V7 but not earlier versions
c       works on Digital UNIX V4.0F 
c       doesn't work on IRIX 6.5
c
      implicit none
c
      character ciftime*(*)
c
      character cdate*8,ctime*10,czone*5
      integer ivalues(8)
c
c ... check if the argument can hold 25 characters
c     (better to return an error flag, of course ;-)
c
      if (len(ciftime) .lt. 25) then
        print *,'error --- hciftime: string too short'
        ciftime = ' '
        return
      end if
c
      CALL Date_and_Time(CDATE,CTIME,CZONE,IVALUES)
c
      write (ciftime,fmt=6000) IVALUES(1),IVALUES(2),IVALUES(3),
     +   IVALUES(5),IVALUES(6),IVALUES(7),CZONE(1:3),CZONE(4:5)
c
c ... NOTE: "i4" in the following format makes that this routine
c           is not Year-10,000-compliant !!!
c
 6000 format (i4,'-',i2.2,'-',i2.2,'T',i2.2,':',i2.2,':',i2.2,
     +   a3,':',a2)
c
      return
      end
