C
C     vms.for
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
C   ---------------- vms.for --------------------
C   replaces vmssupport.for and is replaced by unix.for 
C   for unix systems.
C
C  $Date$
C
C VMS.FOR
C =======
C
C UGTENV - Get value of env. variable
C UGERR  - Get error explanation
C UGTUID - Get user id - it's name
C UIDATE - Get date in 3 integer format
C UTIME  - Get current time
C USTIME - Get absolute time in seconds
C UCPUTM - Get CPU time
C UISATT - Is file a terminal?
C URENAM - Rename file
C VAXVMS - Logical function returns TRUE if VAX/VMS
C WINMVS - Logical function returns TRUE if windows
C UBYTES - Returns number of bytes per word and 'words'/'bytes'
C          to indicate if byte handling is available
C GETPID - Get unique process id.
C USTENV - Create logical name.
C TTSEND - Write string to terminal with various carriage control options
C GETELAPSED - Print timing info for CCPERR
C UGTARG - Get command-line argument 
C GETREF - Abstracted from abscale since it has BYTE declaration.
C CCPSPW - Spawns a new process to run DCL command
C CCPAL1 - Support for CCPALC interface
C CEXIT  - Trivial interface to EXIT
C
C
C     ================================                         
      SUBROUTINE UGTENV(NAMENV,VALENV)
C     ================================
C
C UGTENV - Get value of env. variable
C
C Input:  NAMENV - Logical Name
C
C Output: VALENV - It's value
C
C Arguments: CHARACTER*(*) NAMENV, VALENV
C
C Usage:     CALL UGTENV(NAMENV, VALENV)
C
      CHARACTER*(*) NAMENV,VALENV
C
      INCLUDE '($LNMDEF)'
      INCLUDE '($SSDEF)'
C
      INTEGER       LN,LENGTH
      INTEGER*4     ITEMLIST(4),SYS$TRNLNM
      INTEGER*2     NAME_LEN_CODE(2)
C
C---- Equivalences
C
      EQUIVALENCE (NAME_LEN_CODE(1),ITEMLIST(1))
C                                  
      VALENV = ' '
      LN = LENSTR(NAMENV)
      IF (LN.LE.0) RETURN
C
C---- Setup item list for routine
C
      NAME_LEN_CODE(1) = LEN(VALENV) ! Length of buffer
      NAME_LEN_CODE(2) = LNM$_STRING ! item code for returning equivalence name
      ITEMLIST(2) = %LOC(VALENV)     ! Address to return equivalence name
      ITEMLIST(3) = %LOC(LENGTH)     ! Address to return name length
      ITEMLIST(4) = 0                ! terminator
C
C
10    IERR=SYS$TRNLNM(LNM$M_CASE_BLIND,'LNM$DCL_LOGICAL',
     .    NAMENV(1:LN),,ITEMLIST)
C
      END
C
C
C     ===============================
      SUBROUTINE UGERR(STATUS,ERRSTR)
C     ===============================
C
C UGERR - Get error message string for error number in status
C
C Input:  STATUS - Error number (if negative then print error message,
C     if zero, then use last error)
C
C Output: ERRSTR - Error message string
C
C Arguments: INTEGER       STATUS
C            CHARACTER*(*) ERRSTR
C
C Usage:     CALL UGERR(STATUS, ERRSTR)
C
      INCLUDE       '($SSDEF)'
      INTEGER       STATUS,IPRINT,ISTAT
      CHARACTER*(*) ERRSTR
      CHARACTER*100 OUTLIN
      INTEGER       ISTART,IEND,COND,IFLAGS,IRET
      INTEGER LENSTR
      EXTERNAL LENSTR
C
C
C---- IFLAGS masks out irrelevant parts if the error message
C
      IFLAGS = 13
C
C---- Remember STATUS because a call to ERRSNS destroys it !
C
      ISTAT = STATUS
C
C---- Set up print option
C
      IPRINT = 0                                  
      IF (ISTAT .LT. 0) THEN
        IPRINT = 1
        ISTAT = -ISTAT
      ENDIF
C
C---- Get error message from system
C
      IF (ISTAT.NE.0) THEN
        CALL ERRSNS (ISTAT,,,,COND)
      ELSE
        CALL ERRSNS (,,,,COND)
      ENDIF
C
C---- Translate it
C                
      IRET = LIB$SYS_GETMSG(COND,ILEN,ERRSTR,IFLAGS)
C
C---- If not a fortran error then get system error instead
C
      IF (IRET .EQ. SS$_MSGNOTFND) 
     +    IRET = LIB$SYS_GETMSG(ISTAT,ILEN,ERRSTR,IFLAGS)
C      WRITE (6,*) 'UGERR: ',ERRSTR(:LENSTR(ERRSTR))
C
C---- Remove rubbish
C
      ISTART = INDEX(ERRSTR,' ') + 1
      IEND   = ISTART + INDEX(ERRSTR(ISTART:),'!') - 2
      IF (IEND.LT.ISTART) IEND = LEN(ERRSTR) 
      ERRSTR = ERRSTR(ISTART:IEND)
C
C---- Print result if appropriate
C
      IF (IPRINT.EQ.1) THEN
        WRITE (OUTLIN,100) ISTAT
        OUTLIN(LENSTR(OUTLIN)+2:) = ERRSTR
        WRITE(6,FMT='(A)') OUTLIN(1:LENSTR(OUTLIN))
100     FORMAT (' OS error: ',I5,' Message: ')
      ENDIF
C
      END
C
C
C     ===========================
      SUBROUTINE UGTUID(USERNAME)
C     ===========================
C
C UGTUID - Get user ID
C
C Input:  none
C
C Output: UID - user ID string
C               
C Arguments: CHARACTER*(*) UID
C
C Usage:     CALL UGTUID(UID)
C
      CHARACTER*(*) USERNAME
C
      INCLUDE '($JPIDEF)'
C
      CALL LIB$GETJPI(JPI$_USERNAME,,,,USERNAME)
C
      END
C
C
C     =================================
      SUBROUTINE UIDATE(MONTH,DAY,YEAR)
C     =================================
C
C UIDATE - Get date in 3 integer format
C
C Input:  none
C
C Output: MONTH,DAY,YEAR
C
C Arguments: INTEGER MONTH, DAY, YEAR
C
C Usage:     CALL UIDATE(MONTH, DAY, YEAR)
C
      INTEGER     MONTH,DAY,YEAR
C
      CALL IDATE(MONTH,DAY,YEAR)
C
      END
C
C
C     =======================
      SUBROUTINE UTIME(CTIME)
C     =======================
C
C UTIME - Get current time  as hh:mm:ss
C
C Input:  none
C         
C Output: TIME - as ASCII string
C 
C Arguments: CHARACTER*(*) TIME
C
C Usage:     CALL UTIME(TIME)
C
      CHARACTER*(*) CTIME
C
      CALL TIME(CTIME)
C
      END
C
C
C     ==============================
      SUBROUTINE UISATT(FLUN,ANSWER)
C     ==============================
C
C UISATT - This function determines whether a program is being run 
C          on-line if this information is available.
C
C Input:  FLUN - Fortran Unit Number
C
C Output: ANS - 1 for on-line, 0 otherwise
C
C Arguments: INTEGER FLUN, ANS
C
C Usage:     CALL UISATT (FLUN,ANS)
C
      INTEGER   FLUN,ANSWER
      INTEGER*2 LENGTH,CODE,ITEMLIST(8),RLN
      INTEGER   BUFADDR,RLNADDR,SYS$GETJPI
      CHARACTER ERRSTR*100
C
C---- Equivalences
C
      EQUIVALENCE (ITEMLIST(1),LENGTH),(ITEMLIST(2),CODE),
     .  (ITEMLIST(3),BUFADDR),(ITEMLIST(5),RLNADDR),(ITEMLIST(7),JEND)
      INTEGER BUF(5)
      INCLUDE '($JPIDEF)'
C
C---- Set up item list
C
      LENGTH=20                 ! Length of return buffer in bytes
      CODE=JPI$_TERMINAL        ! Code for information required ( = '31D'X)
      BUFADDR=%LOC(BUF)         ! Address of return buffer
      RLNADDR=%LOC(RLN)         ! Address to receive length of returned 
      JEND=0                    ! Terminator of item list
      I=SYS$GETJPI(,,,%REF(LENGTH),,,)
      ANSWER = 0
C
C---- Set mode. Length of information = 0 if in batch
C
      IF (RLN.NE.0) ANSWER = 1
C
      END
C
C
C     =========================
      LOGICAL FUNCTION VAXVMS()
C     =========================
C
C VAXVMS - Operating Sytem in use returns .TRUE. if VAXVMS
C
C Input:  none
C
C Returns: .TRUE. for VAXVMS, .FALSE. otherwise
C
C Arguments: none
C
C Usage:     VAXVMS ()
C
      VAXVMS = .TRUE.
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
      WINMVS = .FALSE.
C
      END
C
C
C     ===============================
      SUBROUTINE UBYTES (INUM,STRING)
C     ===============================
C
C UBYTES - Return statistics about byte handling
C
C Input:  none
C
C Output: INUM - number of bytes per word
C         HANDLE - 'WORDS' or 'BYTES'
C
C Arguments: INTEGER     INUM
C            CHARACTER*5 HANDLE
C
C Usage:     CALL UBYTES (INUM,HANDLE)
C
      INTEGER INUM
      CHARACTER STRING*5
C
      INUM = 4
      STRING = 'WORDS'
C
      END
C
C
C     =====================================
      SUBROUTINE URENAM(NAME1,NAME2,STATUS)
C     =====================================
C
C URENAM - Rename file assigned to NAME1 to NAME2. 
C
C Input:  NAME1, NAME2   the file names
C
C Output: STATUS       =0 if ok <>0 if an error
C
C Arguments: CHARACTER*(*) NAME1, NAME2
C            INTEGER       STATUS
C
C Usage:     CALL URENAM (NAME1,NAME2,STATUS)
C
      INCLUDE '($SSDEF)'
C
      INTEGER       STATUS
      CHARACTER*(*) NAME1,NAME2
C
C---- Rename file
C
      STATUS = LIB$RENAME_FILE (NAME1,NAME2)
      IF (STATUS .EQ. SS$_NORMAL) STATUS = 0
C
      END
C
C
C
C     ==============================
      INTEGER FUNCTION IARGC(IDUMMY)
C     ==============================
C
C---- Return number of Command line arguments
C
C     ..
C     .. Scalar Arguments ..
      INTEGER IDUMMY
C     ..
C     .. Scalars in Common ..
      INTEGER IARG
C     ..
C     .. External Subroutines ..
      EXTERNAL INITFYP
C     ..
C     .. Common blocks ..
      COMMON /ARGCOUNT/IARG
C     ..
C     .. Save statement ..
      SAVE
C     ..
      IF (IARG.EQ.-1) CALL INITFYP
      IARGC = IARG
C
      END
C
C
C
C     ============================
      SUBROUTINE GETARG(INUM,LINE)
C     ============================
C
C---- Get INUM'th command line argument or ' ' into LINE
C
C     .. Parameters ..
      INTEGER MAXLEN,MAXPAR
      PARAMETER (MAXLEN=70,MAXPAR=41)
C     ..
C     .. Scalar Arguments ..
      INTEGER INUM
      CHARACTER LINE* (*)
C     ..
C     .. Scalars in Common ..
      INTEGER IARG
C     ..
C     .. Arrays in Common ..
      CHARACTER ARGNAM* (MAXLEN)
C     ..
C     .. External Subroutines ..
      EXTERNAL INITFYP
C     ..
C     .. Common blocks ..
      COMMON /ARGCOUNT/IARG
      COMMON /ARGS/ARGNAM(MAXPAR)
C     ..
C     .. Save statement ..
      SAVE
C     ..
      IF (IARG.EQ.-1) CALL INITFYP
      IF (INUM.LT.0 .OR. INUM.GT.IARG) THEN
        LINE = ' '
      ELSE
        LINE = ARGNAM(INUM + 1)
      END IF
C
      END
C
C
C
C     ============================
      INTEGER FUNCTION GETPID(IDUMMY)
C     ============================
C
C     Get process ID
C
      IMPLICIT NONE
      INTEGER PID,IDUMMY
C
      INTEGER ISTAT
C
      INTEGER SYS$GETJPI
C
      INCLUDE '($SSDEF)'
C
      PID = 0
      ISTAT = SYS$GETJPI(,PID,,0,,,)
      IF (ISTAT.NE.SS$_NORMAL) CALL LIB$SIGNAL(%VAL(ISTAT))
      GETPID = PID
C
      END
C
C
C
C     ===============================
      SUBROUTINE USTENV(LINE,IRESULT)
C     ===============================
C
C     Logical name assignment
C     LINE is '<logical name>=<filename>'.  IRESULT.eq.0 iff successful
C
      IMPLICIT NONE
      INTEGER IRESULT
      CHARACTER LOGNAM*80,FILNAM*200,LINE*(*)
C
      STRUCTURE /ITMLST/
        UNION
          MAP
            INTEGER*2 LB,IC
            INTEGER*4 IA,IR
          ENDMAP
          MAP
            INTEGER*4 IE
          ENDMAP
        ENDUNION
      ENDSTRUCTURE
C
      RECORD /ITMLST/ ITM(2)
C
      INTEGER ISTAT,ILEN
C
      INTEGER LENSTR,SYS$CRELNM
C
      INCLUDE '($SSDEF)'
      INCLUDE '($LNMDEF)'
C
      IRESULT = 0
      ILEN = INDEX(LINE, '=')
      IF (ILEN.EQ.0) THEN
        IRESULT = 1 
        RETURN
      ENDIF
      LOGNAM = LINE(1: ILEN - 1)
      FILNAM = LINE(ILEN + 1:)
      ITM(1).LB = LENSTR(FILNAM)
      ITM(1).IC = LNM$_STRING
      ITM(1).IA = %LOC(FILNAM)
      ITM(2).IE = 0
      ISTAT = SYS$CRELNM(,'LNM$PROCESS',LOGNAM(:LENSTR(LOGNAM)),,ITM)
      IF (ISTAT.NE.SS$_NORMAL .AND. ISTAT.NE.SS$_SUPERSEDE) IRESULT = 2
C
      END
C
C
C
      INTEGER FUNCTION SRTBEG(NKEYS,KEYB,LRECL,MEMSIZE)
C     *************************************************
      IMPLICIT NONE
C     # of keys, key descriptions, fixed record length, memory (not used)
      INTEGER   NKEYS,KEYB(*),LRECL,MEMSIZE
C
      INTEGER      ISTAT,NORMAL,DATASIZ,LUNSTO,LUNOUT,NFILSZ
      INTEGER      I,J,JOLD
      INTEGER*2    KEYBUF(401)
C
C    .. External Functions ..
      INTEGER  SOR$BEGIN_SORT
C
C     Definition of data type = single-precision floating - only this one
C     is implemented here
      EXTERNAL DSC$K_DTYPE_F,LUNSTO
C
C   Things for descriptor of ADATA
      INTEGER*4    MDATA(2)
      INTEGER*2    M2DATA(4)
      EQUIVALENCE  (MDATA(1),M2DATA(1))
      COMMON /BISRTC/MDATA,NORMAL,DATASIZ,LUNOUT
C
      SAVE KEYBUF
C
C
C---- Set up key buffer for number of keys. A descriptor consists of:
C     2-byte length (in bytes)
C Each following "block" contains:
C     2-byte data type (here real)
C     2-byte ascending (0) or descending (1) order
C     2-byte relative address of key in record (from 0)
C     2-byte key length in bytes
C
C     NORMAL return value from VMS sort subroutines
      NORMAL = 1
C     Length of tada type i.e. 4 for REAL
      DATASIZ = 4
      LUNOUT = LUNSTO(1)
      KEYBUF(1) = NKEYS
      DO 10 I = 1,NKEYS
         J = (I-1)*4 + 2
         JOLD = (I-1)*5 + 1
C
C---- Sort Data Type
C
         IF (KEYB(JOLD).NE.7) THEN
            WRITE (LUNOUT,FMT=6010)
 6010       FORMAT (' SRTBEG only REAL data type implemented')
            STOP
         END IF
         KEYBUF(J)=%LOC(DSC$K_DTYPE_F)
C
C---- Sort Order ascending/descending
C
         KEYBUF(J+1) = KEYB(JOLD+1)
C
C---- position of 1st byte in key
C
         KEYBUF(J+2) = KEYB(JOLD+2)
C
C---- keylength in BYTES
C
         KEYBUF(J+3) = DATASIZ * KEYB(JOLD+3)
         IF (KEYB(JOLD+4).NE.0) THEN
            WRITE (LUNOUT,FMT=6011)
 6011       FORMAT(' SRTBEG - on VMS MASK fields must be .EQ. 0')
            STOP
         ENDIF
 10   CONTINUE
C
C
C     Make string descriptors for data record
C     A descriptor consists of:
C     2-byte length (in bytes)
C     2-byte class & type (probably not used)
C     4-byte address of array
C     Note MDATA is equivalenced to M2DATA
C     length of array
      M2DATA(1)=LRECL
C     class = 1, type = 0 - never used but must be present
      M2DATA(2)='100'X
C     address of array  - filled in SRTRLS & SRTRET
C
C---- Initialise sort, set parameters, etc
C
C    ******************************
      NFILSZ = 0
      ISTAT = SOR$BEGIN_SORT(KEYBUF,LRECL,,NFILSZ,,,,,)
C    ******************************
C
      IF (ISTAT.NE.NORMAL) THEN
         WRITE (LUNOUT,FMT=6008) ISTAT
 6008    FORMAT (' Sort fail : BEGIN, status=',Z9)
C
         CALL LIB$STOP(%VAL(ISTAT))
      END IF
      SRTBEG = 0
      RETURN
      END
C
C
C
      INTEGER FUNCTION SRTRLS(RECORD)
C     *******************************
      IMPLICIT NONE
      REAL  RECORD(*)
C
      INTEGER ISTAT,NORMAL,DATASIZ,LUNOUT
C
      INTEGER  SOR$RELEASE_REC
C
C     MDATA is descriptor (ie indirect address) of RECORD
C     Things for descriptor of ADATA
      INTEGER*4    MDATA(2)
      INTEGER*2    M2DATA(4)
      EQUIVALENCE  (MDATA(1),M2DATA(1))
      COMMON /BISRTC/MDATA,NORMAL,DATASIZ,LUNOUT
C
      MDATA(2)=%LOC(RECORD)
C     *************
      ISTAT = SOR$RELEASE_REC(MDATA)
C     *************
C
C     IF (ISTAT.EQ.0) THEN
      IF (ISTAT.EQ.NORMAL) THEN
         SRTRLS = 0
      ELSE
         WRITE (LUNOUT,FMT=6010) ISTAT
 6010    FORMAT (' Sort fail : SRTRLS, status=',Z9)
C
         CALL LIB$STOP(%VAL(ISTAT))
         STOP
      END IF
      RETURN
      END
C
C
C
      INTEGER FUNCTION SRTMRG()
C     *************************
      IMPLICIT NONE
C
      INTEGER  ISTAT,NORMAL,DATASIZ,LUNOUT
      INTEGER  SOR$SORT_MERGE
C
C     MDATA is descriptor (ie indirect address) of RECORD
C     Things for descriptor of ADATA
      INTEGER*4    MDATA(2)
      INTEGER*2    M2DATA(4)
      EQUIVALENCE  (MDATA(1),M2DATA(1))
      COMMON /BISRTC/MDATA,NORMAL,DATASIZ,LUNOUT
C
C     ********
      ISTAT = SOR$SORT_MERGE()
C     ********
C
C     IF (ISTAT.NE.0) THEN
      IF (ISTAT.NE.NORMAL) THEN
         WRITE (LUNOUT,FMT=6014) ISTAT
 6014    FORMAT (' Sort fail : MERGE, status=',Z9)
         CALL LIB$STOP(%VAL(ISTAT))
         STOP
      ENDIF
      SRTMRG = 0
      RETURN
      END
C
C
C
      INTEGER FUNCTION SRTRET(RECORD)
C     *******************************
      IMPLICIT NONE
C     record array
      REAL RECORD(*)
C
      INTEGER ISTAT,NORMAL,DATASIZ,LUNOUT
C
      INTEGER  SOR$RETURN_REC
      INTEGER  SOR$END_SORT
      EXTERNAL SS$_ENDOFFILE
C  NRL not used, but still present (the length, which is already known)
      INTEGER*2  NRLVMS
C
C     MDATA is descriptor (ie indirect address) of RECORD
C     Things for descriptor of ADATA
      INTEGER*4    MDATA(2)
      INTEGER*2    M2DATA(4)
      EQUIVALENCE  (MDATA(1),M2DATA(1))
      COMMON /BISRTC/MDATA,NORMAL,DATASIZ,LUNOUT
C
C     set record address
      MDATA(2)=%LOC(RECORD)
C     *************
      ISTAT = SOR$RETURN_REC(MDATA,NRLVMS)
C     *************
C
      IF (ISTAT.EQ.%LOC(SS$_ENDOFFILE)) THEN
         ISTAT = SOR$END_SORT()
         IF(ISTAT.NE.NORMAL) THEN
            WRITE(6,1005) ISTAT
 1005       FORMAT(' Sort fail : END ,status=',Z9)
            CALL LIB$STOP(%VAL(ISTAT))
         ENDIF
         SRTRET = -1
         RETURN
      ELSE IF (ISTAT.NE.NORMAL) THEN
         WRITE(6,1006) ISTAT
 1006    FORMAT(' Sort fail : RETURN, status=',Z9)
         CALL LIB$STOP(%VAL(ISTAT))
         STOP
      ENDIF
      SRTRET = 0
      RETURN
      END
C     
C     
C     
      SUBROUTINE CCPOPN(IIUN,LOGNAM,KSTAT,ITYPE,LREC,IFAIL)
C     ====================================================
C
C---- This subroutine is used to open a file
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
      INTEGER LLREC,IUN,IBYTES,ISTAT
C     ERRSTR should be big enough to hold more than 2 255-long paths
      CHARACTER CCNTRL*7,ST*7,FRM*12,ERRSTR*600,FULNAM*255,
     +     NAMFIL*255,HANDLE*5,OPNVAR*20, ACCESS*10, DISPOS*6
      INTEGER UNKNWN, SCRTCH, OLD, NEW, RDONLY, PRINTR
      PARAMETER (UNKNWN=1, SCRTCH=2, OLD=3, NEW=4, RDONLY=5, PRINTR=6)
      LOGICAL LNONAM
C     ..
C     .. Local Arrays ..
      CHARACTER STAT(6)*7, OUTLIN*100
C     ..
C     .. External Functions ..
      INTEGER LENSTR, LUNSTO
      LOGICAL VAXVMS,CCPEXS
      EXTERNAL LENSTR,VAXVMS,LUNSTO,CCPEXS
C     ..
C     .. External Subroutines ..
      EXTERNAL UGERR,UGTENV
C     ..
C     .. Data statements ..
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
          CALL CCPERR(LUNSTO (1),
     +         '**CCPOPN ERROR** Invalid parameters in call')
        ELSE
          WRITE (LUNSTO (1),
     +         '('' **CCPOPN ERROR** Invalid parameters in call'',/)')
          IFAIL = -1
        END IF
        RETURN
      ENDIF 
C
C     Do nothing for pre-connected units (what's the significance of
C     `TERM...'?) 
      IF (LOGNAM.EQ.'DATA' .OR. LOGNAM.EQ.'PRINTER' .OR.
     $     LOGNAM(:4).EQ.'TERM') RETURN
C
C     if environment variable CCP4_OPEN has value `UNKNOWN', open files
C     with status UNKNOWN rather than new if they exist
      IF (ISTAT.EQ.NEW) THEN
        OPNVAR = ' '
        CALL UGTENV('CCP4_OPEN',OPNVAR)
        IF (OPNVAR.EQ.'UNKNOWN') ISTAT = 1
      END IF
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
      IF (ISTAT.EQ.SCRTCH) THEN
C       scratch file
        DISPOS = 'DELETE'
      ELSE
        DISPOS = 'KEEP'
      ENDIF
C
C     check for `logical name' referencing real file
      CALL UGTENV(LOGNAM,NAMFIL)
      LNONAM = .FALSE.
      IF (NAMFIL.EQ.' ') THEN
        IF (.NOT. CCPEXS(LOGNAM)) LNONAM = .TRUE.
        NAMFIL = LOGNAM
      END IF
C     Unix null device (defined as canonical if programs need it)
      IF (NAMFIL.EQ.'/dev/null') NAMFIL = 'NL:'
C     Opening /dev/null is necessary in Unix; not sure if this is needed
C     but presumably does no harm.
      IF (NAMFIL.EQ.'NL:') ISTAT = 1
C       
      IF (ACCESS.EQ.'DIRECT') THEN
C       Need to check is record length in words or bytes and set LLREC
C       accordingly. 
        CALL UBYTES (IBYTES,HANDLE)
        LLREC = LREC*IBYTES
        IF (HANDLE.EQ.'WORDS'.AND.ITYPE.EQ.4) LLREC=LLREC/IBYTES
        IF (ISTAT.EQ.RDONLY) THEN
C         READONLY, may be defined as null or as `READONLY,'
          OPEN(UNIT=IUN,STATUS='UNKNOWN',ACCESS='DIRECT',FORM=FRM,
     +         READONLY,
     +         FILE=NAMFIL,RECL=LLREC,IOSTAT=IOS,ERR=5)
        ELSE
          OPEN(UNIT=IUN,STATUS='UNKNOWN',ACCESS='DIRECT',FORM=FRM,
     +         DISPOSE = DISPOS,
     +         FILE=NAMFIL,RECL=LLREC,IOSTAT=IOS,ERR=5)
        ENDIF 
      ELSE
C       carriagecontrol='fortran' for print file, else = 
C       'list'. 
        IF (ISTAT.EQ.PRINTR) THEN
C         want to obey format characters in column 1
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
     +           FORM=FRM, ERR=5, IOSTAT=IOS, DISPOSE=DISPOS)
          ENDIF
        ELSE
          IF (ISTAT.EQ.RDONLY) THEN
            OPEN(UNIT=IUN, FILE=NAMFIL, STATUS=ST, ACCESS='SEQUENTIAL',
     +           READONLY, DISPOSE=DISPOS,
     +           CARRIAGECONTROL=CCNTRL,
     +           FORM=FRM, ERR=5, IOSTAT=IOS)
          ELSE
            OPEN(UNIT=IUN, FILE=NAMFIL, STATUS=ST, ACCESS='SEQUENTIAL',
     +           CARRIAGECONTROL=CCNTRL, DISPOSE=DISPOS,
     +           FORM=FRM, ERR=5, IOSTAT=IOS)
          ENDIF
        ENDIF
      ENDIF
C     Error check
 5    CONTINUE
C     don't report UNKNOWN if actually SCRATCH
      IF (ISTAT.EQ.SCRTCH) ST = 'SCRATCH'
      IF (IOS.NE.0) THEN
        IF (IFAIL.EQ.0) THEN
C         warning if there was no file associated with logical name
          IF (LNONAM) THEN
             ERRSTR = 'CCPOPN Logical name '//LOGNAM
             ERRSTR(LENSTR(ERRSTR)+2:) = 'has no associated file name'
             CALL CCPERR(2,ERRSTR)
          END IF
C         hard failure
          WRITE (ERRSTR,FMT=6002) IUN, NAMFIL(1:LENSTR(NAMFIL)),
     +         LOGNAM(1:LENSTR(LOGNAM))
 6002     FORMAT ('Open failed: Unit:',I4,', File: ',A, ' (logical: ', A
     +         , ')')
          CALL QPRINT (0, ERRSTR)
          ERRSTR = ' Open failed: File: ' // NAMFIL
          CALL CCPERR(1, ERRSTR)
        ELSE
C         soft failure
          IF (IFAIL.EQ.1) WRITE (LUNSTO (1),FMT=6004) FRM, ST, IUN, 
     +         LOGNAM(1:LENSTR(LOGNAM)), NAMFIL(1:LENSTR(NAMFIL))
 6004     FORMAT (' **CCPOPN ERROR**  ',A,3X,A,
     +         ' file open failure on unit ',I3)
          ERRSTR = 'Logical name: ' //LOGNAM
          CALL QPRINT (0, ERRSTR)
          ERRSTR = 'File name: ' // NAMFIL
          CALL QPRINT (0, ERRSTR)
          CALL UGERR(IOS,ERRSTR)
          CALL QPRINT (0, ERRSTR)
          CALL QPRINT (0, ' ')
          IFAIL = -1
          RETURN            
        ENDIF
      ELSE
        IF (IIUN.LE.0) RETURN 
        INQUIRE (FILE=NAMFIL,NAME=FULNAM)
C       DJGL: why is this inquire necessary rather than using NAMFIL?
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
C SUBROUTINE 'TTSEND'
C ===================
C
C Write a string to a terminal with various carriage control options
C [for LAUE]
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
C====== Specification statements
C
      INTEGER IUN, ICC
      CHARACTER*(*) STR
      INTEGER LENSTR
      LOGICAL VAXVMS
      EXTERNAL CCPERR, LENSTR, VAXVMS
C
C====== Write string
C
      IF (VAXVMS()) THEN
        IF (LENSTR(STR) .GT. 132)
     .      CALL CCPERR(1,' TTSEND: Output string is greater than 132')
      ENDIF

      IF (ICC.EQ.0) THEN
         WRITE (IUN,1000) STR
      ELSE IF (ICC.EQ.2) THEN
         WRITE (IUN,1002) STR
      ELSE IF (ICC.EQ.3) THEN
         WRITE (IUN,1003) STR
      ELSE
         WRITE (IUN,1001) STR
      ENDIF
C
C====== Format statements
C
1000  FORMAT (' ',A,$)
1001  FORMAT (' ',A)
1002  FORMAT ('+',A,$)
1003  FORMAT ('+',A)
      END
C
      SUBROUTINE UGTARG(I, ARG)
      INTEGER I
      CHARACTER *(*) ARG
      CALL GETARG(I, ARG)
      END

C     ========================
      SUBROUTINE GETELAPSED
C     ========================
C
C GETELAPSED - print CPU and ELAPSED times since job started.
C
C**** NOTE - Code is not VAX/VMS specific but will only work correctly on VAX,
C**** because on other systems USTIME returns system time relative to arbitrary
C**** zero, whereas in VMS.FOR USTIME has been set up to return time relative
C**** to start of job.
C==== 8-NOV-93 Made s/r INITFYP entry in GETELAPSED to initialise elapsed time
C     to bring into line with unix.for.
C
C     .. Local Scalars .. (GETELAPSED)
      REAL       CPUTIM
      INTEGER    CPUMIN, CPUSEC, CPUTIC, JOBMIN, JOBSEC, JOBSAV
      LOGICAL    INITED
C     ..
C     .. Parameters .. (INITFYP)
      INTEGER MAXLEN,MAXPAR
      PARAMETER (MAXLEN=70,MAXPAR=41)
C     ..
C     .. Scalars in Common ..
      INTEGER IARG
C     ..
C     .. Arrays in Common ..
      CHARACTER ARGNAM* (MAXLEN)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,K,L,LENARG,ISTAT
      CHARACTER CLIARG*700,NAME*200
C     ..
C     .. External Subroutines ..
      EXTERNAL LIB$GET_FOREIGN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC INDEX
C     ..
C     .. Common blocks ..
      COMMON /ARGCOUNT/IARG
      COMMON /ARGS/ARGNAM(MAXPAR)
C     ..
C     .. Save statement ..
      SAVE CPUTIM, INITED, JOBSAV
C     ..
C     .. Data statements ..
      DATA INITED /.FALSE./
      DATA CPUTIM /0./
      DATA IARG/-1/
C     ..
      INCLUDE '($SSDEF)'
      INCLUDE '($JPIDEF)'
C
C     This should be included (either VAX C or DEC C version, determined
C     by the logical name) to link with the C version of the library.
C     However, it's not clear that it's actually necessary...
      INCLUDE 'CRTLINIT'
C     
C     Don't print anything if it hasn't been initialised (by CCPFYP/INITFYP).
C
      CALL UCPUTM(CPUTIM)
      IF (INITED) THEN
        CPUTIC = NINT(100.*CPUTIM)
        CPUMIN = CPUTIC/6000
        CPUSEC = MOD(CPUTIC/100,60)
        CPUTIC = MOD(CPUTIC,100)
C
        CALL USTIME(JOBSEC)
        JOBSEC = JOBSEC - JOBSAV
        JOBMIN = JOBSEC/60
        JOBSEC = MOD(JOBSEC,60)
C
        WRITE (6,'(/,A,I6,2(A,I2.2),8X,A,I6,A,I2.2)')
     &  ' CPU time used:',CPUMIN,':',CPUSEC,'.',CPUTIC,
     &  'Elapsed time:',JOBMIN,':',JOBSEC
      ENDIF
C
      CALL USTIME(JOBSAV)
      INITED = .TRUE.
C
      RETURN
C
C
C     =============
      ENTRY INITFYP
C     =============
C
      CPUTIM = 0.
      CALL UCPUTM(CPUTIM)
      CALL USTIME(JOBSAV)
      INITED = .TRUE.
C
C---- Parse CLI argument: get command line
C
      IARG = 0
      DO 5 J = 1, MAXPAR
        ARGNAM(J) = ' '
5     CONTINUE
C
C---- get user id and use for argv[0]
C
      CALL UGTUID(ARGNAM(1))
C
      CALL LIB$GET_FOREIGN(CLIARG,,LENARG)
C
C---- Split command line into arguments.
C
      IF (LENARG.GT.0) THEN
        J = 1
   10   CONTINUE
        K = INDEX(CLIARG(J:LENARG),' ')
        IF (K.EQ.0) THEN
          K = LENARG
        ELSE
          K = J + K - 2
        END IF
        IARG = IARG + 1
        IF (IARG.EQ.MAXPAR) RETURN
        ARGNAM(IARG + 1) = CLIARG(J:K)
        DO 20 J = K + 2,LENARG
          IF (CLIARG(J:J).NE.' ') GO TO 10
   20   CONTINUE
      END IF
C
      END

      SUBROUTINE UCPUTM(CPUTIM)
C     =========================
C
C     Get CPU time in seconds.
C
C     Argument:
C     REAL CPUTIM (i/o): If <= 0, initialize timer and return current
C                        elapsed cpu time since start of execution, otherwise
C                        return elapsed cpu since timer was initialized.
C
C******************  VAX/VMS SPECIFIC !  *********************
C
      INCLUDE        '($JPIDEF)'
      INTEGER        IOSB(2), IS, IT
      REAL           CPUSAV, CPUTIM
      INTEGER        SYS$GETJPIW
C
      STRUCTURE      /ITMLST/
        UNION
          MAP
            INTEGER*2    LB,IC
            INTEGER      IA,IR
          ENDMAP
          MAP
            INTEGER      IE
          ENDMAP
        ENDUNION
      ENDSTRUCTURE
C
      RECORD         /ITMLST/ JPI(2)
C
      DATA           IT /0/, CPUSAV /0./
      SAVE           IT, CPUSAV
C
      IF (IT.EQ.0) THEN
        JPI(1).LB = 4
        JPI(1).IC = JPI$_CPUTIM
        JPI(1).IA = %LOC(IT)
        JPI(2).IE = 0
      ENDIF
C
      IS = SYS$GETJPIW(,,,JPI,IOSB,,)
      IF (IS) IS = IOSB(1)
      IF (.NOT.IS) CALL LIB$SIGNAL(%VAL(IS))
C
      IF (CPUTIM.LE.0.) THEN
        CPUSAV = .01*IT
        CPUTIM = CPUSAV
      ELSE
        CPUTIM = .01*IT - CPUSAV
      ENDIF
      END

C     =========================
      SUBROUTINE USTIME(JOBTIM)         
C     =========================
C
C USTIME - Get elapsed job time in seconds to nearest second.
C
C Argument: INTEGER JOBTIM
C
C
C******************  VAX/VMS SPECIFIC !  *********************
C
      INCLUDE        '($JPIDEF)'
      LOGICAL        LT
      INTEGER        IOSB(2), IS, JOBTIM, T0(2), TN(2)
      REAL           DT
      INTEGER        SYS$GETJPIW, SYS$GETTIM
C
      STRUCTURE      /ITMLST/
        UNION
          MAP
            INTEGER*2    LB,IC
            INTEGER      IA,IR
          ENDMAP
          MAP
            INTEGER      IE
          ENDMAP
        ENDUNION
      ENDSTRUCTURE
C
      RECORD         /ITMLST/ JPI(2)
C
      DATA           LT /.TRUE./
      SAVE           LT, T0
C
      IF (LT) THEN
        JPI(1).LB = 8
        JPI(1).IC = JPI$_LOGINTIM
        JPI(1).IA = %LOC(T0)
        JPI(2).IE = 0
C
        IS = SYS$GETJPIW(,,,JPI,IOSB,,)
        IF (IS) IS = IOSB(1)
        IF (.NOT.IS) CALL LIB$SIGNAL(%VAL(IS))
        LT = .FALSE.
      ENDIF
C
      IS = SYS$GETTIM(TN)
      IF (.NOT.IS) CALL LIB$SIGNAL(%VAL(IS))
C
C==== T0 = Absolute time at job login in 10^-7 sec units (INTEGER*8).
C==== TN = Absolute time now in 10^-7 sec units (INTEGER*8).
C==== JOBTIM = Time difference in seconds (INTEGER*4).
C
      DT = REAL(TN(1)) - REAL(T0(1))
      IF (T0(1).LT.0) DT = DT - 2.**32
      IF (TN(1).LT.0) DT = DT + 2.**32
      JOBTIM = NINT(1.E-7*(DT + 2.**32*(REAL(TN(2)) - REAL(T0(2)))))
      END
C     
C     =====================================================
      SUBROUTINE GETREF(KERFLG,NREAD,NSPOTS,DYNAM,MAXOMITL)
C     =====================================================
C
C     [This has been abtracted from ABSCALE because of the BYTE
C     declaration.]
C
        implicit none
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
C---- IC format generate file variables
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
      END
C
C
C
C_BEGIN_CCPSPW
      SUBROUTINE CCPSPW(STRING)
C     =========================
C
C     Spawns a new process using shell command
C
C Arguments:
C ==========
C
C  STRING (I)   CHARACTER*(*): string containing command
C_END_CCPSPW
C
       CHARACTER STRING*(*)
       CALL LIB$SPAWN(STRING)
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
        ISTAT = LIB$GET_VM (SIZES(TYPE (I))*LENGTH(I), POINTER(I))
        IF (.NOT.ISTAT)
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
        ISTAT = LIB$FREE_VM (SIZES(TYPE (I))*LENGTH(I), POINTER(I))
        IF (.NOT.ISTAT) CALL CCPERR (-1, 'CCPALC: can''t free memory')
      ENDDO
      END
C
      SUBROUTINE CEXIT (ICODE)
C     trivial interface to system-dependent EXIT routine
      CALL EXIT (ICODE)
      END

C
CA dummy function for unix and vms
C     =========================
       CHARACTER FUNCTION RTNBKS()
C     =========================
C
C RTNBKS - Returns a Backslash for nt as unix compilers are fussy!
C
C Input:     none
C
C Returns:   \ if WIN32 or not if unix or vms
C
C Arguments: none
C
C Usage:     RTNBKS ()
C
      RTNBKS=' '
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
      character cdate*8,ctime*9,czone*5
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

