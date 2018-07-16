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
C**** VAXDISKIO.FOR ****************************************************
C
C       This file contains a set of FORTRAN subroutines for
C       doing either sequential or random access disk I/O
C       with variable record lengths.  The files are actually
C       written as fixed-record direct-access files, but this
C       is transparent to the user.
C
C       Note: Only 10 files may be opened simultaneously!!!!
C                  --
C       This is a VAX version, and will require changes for other machines.
C
C       Within these routines, all counting is done as bytes (== minimal
C       machine element).  However, the entries QREAD, QWRITE, QSEEK, QBACK,
C       QSKIP and QLOCATE count in items.  The conversion between bytes and
C       items may be changed by a call to QMODE, which defines the mode
C       (corresponding to map/image modes), and stores and returns the number
C       of bytes per item NCHITM. NCHITM defaults to 1 for all streams unless
C       changed by a call to QMODE.  This equivalencing of bytes and items
C       should make it possible to write higher level routines which are
C       portable between computers, although the routines here are NOT
C       portable to non-byte computers.
C
C       the calls provided are:
C
C       call QOPEN(IUNIT,NAME,ATTRIBUTE)		;open file
C       call QCLOSE(IUNIT)				;close file
C       call QMODE(IUNIT,MODE,NCHITM)			;set mode
C       call QREAD(IUNIT,ARRAY,NITEMS,IER)		;read nitems
C       call QWRITE(IUNIT,ARRAY,NITEMS)			;write nitems
C       call QSEEK(IUNIT,IRECORD,IELEMENT,IRECLENGTH)	;move to irec,iel
C       call QBACK(IUNIT,IRECLENGTH)			;backspace 1 record
C       call QSKIP(IUNIT,IRECLENGTH)			;skip 1 record
C       call QLOCATE(IUNIT,LOCATION)			;return disk location
C
C       For files created by QOPEN, the user may set numeric values for the
C       logical names CCP4_INITIAL and CCP4_EXTEND to specify the initial
C       (default = 100) and extend sizes (default = 30) for the file.
C
C       QSEEK calculates the location as follows:
C       LOCATION = (IRECORD - 1)*IRECLENGTH + IELEMENT
C
C       NOTE: as in FORTRAN, addressing begins at 1 for BOTH record & element
C
C       where:
C       IUNIT      = Variable used to contain channel # (SET by QOPEN!!!)
C       NAME       = either a LOGICAL or FILE name to be associated with unit
C       ATTRIBUTE  = one of the following: 'NEW' 'OLD' 'RO' 'SCRATCH' 'UNKNOWN'
C       ARRAY      = starting location for data storage in core
C       NITEMS     = number of items to transfer
C       IER        = error flag (0 = no error) else # of bytes transferred
C       IRECORD    = desired record #     (starts @ 1)
C       IELEMENT   = desired element # within record (items)   (starts @ 1)
C       IRECLENGTH = record length in items
C       LENGTH     = returned file size in bytes
C       LOCATION   = returned location of disk pointer in items (starts @ 1)
C       MODE       = file mode, to set correspondence between an item in the
C                    file, and one machine character (byte).
C                    mode 0 = bytes, 1 = short integers (*2), 2 = words
C                    (reals or integers), 3 = short complex (2 integer*2),
C                    4 = complex
C       NCHITM     = returned as number of characters (bytes) / item
C
C
C
C       Last Update:    01.10.82                DAA             VAX
C                       02.05.84                PRE     item counting
C                       24.09.93                IJT     remove macro calls
C
C    For CCP compatibilty, added entries
C         QQOPEN(IUNIT,NAME,ISTAT)     calls QOPEN, sets file status
C                                    from ISTAT = 1 'UNKNOWN'
C                                               = 2 'SCRATCH'
C                                               = 3 'OLD'
C                                               = 4 'NEW'
C                                               = 5 'READONLY'
C                                     Also sets default mode = 2
C
C         QQINQ(IUNIT,LFN,NAME,LENGTH)  get full name & length of file
C
C
      SUBROUTINE DISKIO
C     =================
C
C==== NRBLK = Number of physical blocks/logical block
C.... Note the user may try varying NRBLK in the range 1...127 (=MAXDO),
C     but the other parameters should be left alone!
C==== NCMAX = Maximum number of channels open simultaneously
C==== NPBSZ = Physical block size on disk
C==== MAXDO = Maximum number of physical blocks in one transfer
C==== MAXMO = Maximum mode value
C==== NBALQ = Default initial physical block allocation
C==== NBDEQ = Default extension physical block allocation
C
      IMPLICIT   NONE
      INTEGER    NRBLK, NCMAX, NPBSZ, MAXDO, MAXMO, NBALQ, NBDEQ, NSIZE,
     &           NSIZE1, I2P15, I2P16
      PARAMETER  (NRBLK=64, NCMAX=10, NPBSZ=512, MAXDO=127, MAXMO=6,
     &           NBALQ=100, NBDEQ=30, NSIZE=NRBLK*NPBSZ, NSIZE1=NSIZE+1,
     &           I2P15=2**15, I2P16=2**16)
C
C==== Declare system calls, parameters, record structures etc etc ...
C
      INCLUDE    '($SYSSRVNAM)'
      INCLUDE    '($RMSDEF)'
      INCLUDE    '($FABDEF)'
      INCLUDE    '($NAMDEF)'
      INCLUDE    '($RABDEF)'
      INCLUDE    '($XABDEF)'
      INCLUDE    '($XABFHCDEF)'
C
      RECORD     /FABDEF/ FAB(NCMAX)
      RECORD     /NAMDEF/ NAM(NCMAX)
      RECORD     /RABDEF/ RAB(NCMAX)
C
      STRUCTURE  /MAPXABDEF/
       UNION
        MAP
         RECORD  /XABDEF/ XAB
        ENDMAP
        MAP
         RECORD  /XABFHCDEF/ XABFHC
        ENDMAP
       ENDUNION
      ENDSTRUCTURE
C
      RECORD     /MAPXABDEF/ MAPXAB(NCMAX)
C
      BYTE       ARRAY(*), BUFFER(NSIZE,NCMAX), MS(4)
      CHARACTER  ATBUTE*(*), ATTRIB*1, FILNAM*(*), MSG*80, NAME*(*),
     &           UNK*7
      CHARACTER  FNAME(NCMAX)*127
      LOGICAL    EXIST
      LOGICAL    BUFACT(NCMAX), WRTACT(NCMAX), NEW(NCMAX), LUNIT(NCMAX)
      INTEGER    I, IEL, IER, IREC, IRET, IUNIT, L, LENGTH, LOCATION,
     &           LRECL, MCHITM, MODE, MTIT, MTRT, NBYTES, MODE1, MAXNC
      INTEGER    IBY, INDEX, ISAVE, ISTAT, IVBLK, J, JEL, KEL, KRECSZ,
     &           NBY, NDO, NLEFT, NOBLK, NTOMOV
      INTEGER    LRECSZ(NCMAX), MRECSZ(NCMAX), MODES(0:MAXMO),
     &           MSIT(NCMAX), MSRT(NCMAX), NCHITM(NCMAX), NEL(NCMAX),
     &           NMODE(NCMAX), NRECSZ(NCMAX), NVBLK(NCMAX)
      INTEGER    LENSTR
      REAL       BUF(NSIZE/4,NCMAX)
      EQUIVALENCE(BUFFER,BUF)
      COMMON /DISKIO_MODES_/ MAXNC, NMODE
      DATA MAXNC/NCMAX/
      DATA       UNK/'UNKNOWN'/, LUNIT/NCMAX*.TRUE./
      SAVE
C
C Item length for each mode: these correspond to map/image modes 0 to 6
C          byte, integer*2, real*4, complex*4, complex*8, unused, integer*4
      DATA MODES/ 1,      2,      4,      4,      8,     0,     4 /
C
C
C#######################################################################
C
C
C*QOPEN
      ENTRY QOPEN(IUNIT,NAME,ATBUTE)
C     ==============================
C
C   FIND NEXT AVAILABLE UNIT #
C
      DO 10 J = 1,NCMAX
10      IF (LUNIT(J)) GOTO 20
      WRITE(6,'(/A,I4)') ' Max available unit= ',NCMAX
C
C==== Maybe this should be a soft fail: IUNIT = -1 ?
C
      CALL CCPERR(1,'QOPEN error: No available units !!')
      IUNIT = -1
      RETURN
C
20    IUNIT = J
      LUNIT(IUNIT) = .FALSE.
      NMODE(IUNIT) = 2
      NCHITM(IUNIT) = 4
C
C==== Identify File Access Block.
C
      FAB(IUNIT).FAB$B_BID = FAB$C_BID
      FAB(IUNIT).FAB$B_BLN = FAB$C_BLN
C
C==== File name address & size, clear file-processing options.
C
      FAB(IUNIT).FAB$L_FNA = %LOC(NAME)
      L = MIN(LENSTR(NAME),255)
      IF (L.GT.127) L = L - 256
      FAB(IUNIT).FAB$B_FNS = L
      FAB(IUNIT).FAB$L_FOP = 0
C
C==== File organisation is sequential.
C
      FAB(IUNIT).FAB$B_ORG = FAB$C_SEQ
C
C==== Pretend record format is stream_lf for C's benefit.
C
      FAB(IUNIT).FAB$B_RFM = FAB$C_STMLF
C
C==== NAM and XAB addresses.
C
      FAB(IUNIT).FAB$L_NAM = %LOC(NAM(IUNIT))
      FAB(IUNIT).FAB$L_XAB = %LOC(MAPXAB(IUNIT))
C
C==== Identify NAMe block.
C
      NAM(IUNIT).NAM$B_BID = NAM$C_BID
      NAM(IUNIT).NAM$B_BLN = NAM$C_BLN
C
C==== Expanded string address and size.
C
      NAM(IUNIT).NAM$L_ESA = %LOC(FNAME(IUNIT))
      NAM(IUNIT).NAM$B_ESS = LEN(FNAME(1))
      FNAME(IUNIT) = ' '
C
C==== Identify as eXtended Attributes Block for File Header Characteristics.
C
      MAPXAB(IUNIT).XAB.XAB$B_COD = XAB$C_FHC
      MAPXAB(IUNIT).XAB.XAB$B_BLN = XAB$C_FHCLEN
C
      ATTRIB = ATBUTE
      IF (ATTRIB.EQ.'R') THEN
C
C==== File access mode is block i/o readonly.
C
        FAB(IUNIT).FAB$B_FAC = FAB$M_BIO .OR. FAB$M_GET
C
C==== Open existing file.
C
        ISTAT = SYS$OPEN(FAB(IUNIT))
        NEW(IUNIT) = .FALSE.
      ELSEIF (ATTRIB.EQ.'O') THEN
C
C==== File access mode is block i/o read/write.
C
        FAB(IUNIT).FAB$B_FAC = FAB$M_BIO .OR. FAB$M_GET .OR. FAB$M_PUT
C
C==== Open existing file.
C
        ISTAT = SYS$OPEN(FAB(IUNIT))
        NEW(IUNIT) = .FALSE.
      ELSE
C
C==== Default allocation and extension quantities.
C
        FAB(IUNIT).FAB$L_ALQ = NBALQ
        FAB(IUNIT).FAB$W_DEQ = NBDEQ
C
C==== Reset from user-defined logical names if present.
C
        CALL UGTENV('CCP4_INITIAL',MSG)
        IF (MSG.NE.' ') READ (MSG,*,ERR=30) FAB(IUNIT).FAB$L_ALQ
30      CALL UGTENV('CCP4_EXTEND',MSG)
        IF (MSG.NE.' ') READ (MSG,*,ERR=40) FAB(IUNIT).FAB$W_DEQ
C
C==== File access mode is block i/o read/write.
C
40      FAB(IUNIT).FAB$B_FAC = FAB$M_GET .OR. FAB$M_PUT .OR. FAB$M_BIO
C
        IF (ATTRIB.EQ.'S') THEN
C
C==== Scratch file is temporary/marked-for-delete.
C
          FAB(IUNIT).FAB$L_FOP = FAB$M_TMD
C
        ELSE
C
C==== Create max version & truncate new file at EOF.
C
          FAB(IUNIT).FAB$L_FOP = FAB$M_MXV .OR. FAB$M_TEF
C
          IF (ATTRIB.EQ.'N') THEN
            CALL UGTENV('CCP4_OPEN',MSG)
            L = LENSTR(MSG)
            IF (L.GT.0 .AND. L.LE.7) THEN
              IF (MSG(:L).EQ.UNK(:L)) ATTRIB = 'U'
            ENDIF
          ENDIF
C
C==== If UNKNOWN set Create-IF.
C
          IF (ATTRIB.EQ.'U')
     &    FAB(IUNIT).FAB$L_FOP = FAB(IUNIT).FAB$L_FOP .OR. FAB$M_CIF
        ENDIF
C
C==== Create new or open existing file, depending on CIF setting.
C
        ISTAT = SYS$CREATE(FAB(IUNIT))
        NEW(IUNIT) = .TRUE.
      ENDIF
C
      IF (ISTAT) THEN
C
C==== Identify Record Access Block.
C
        RAB(IUNIT).RAB$B_BID = RAB$C_BID
        RAB(IUNIT).RAB$B_BLN = RAB$C_BLN
C
C==== FAB address and internal stream identifier.
C
        RAB(IUNIT).RAB$L_FAB = %LOC(FAB(IUNIT))
        RAB(IUNIT).RAB$W_ISI = 0
C
C==== Record processing option is block i/o
C
        RAB(IUNIT).RAB$L_ROP = RAB$M_BIO
C
C==== Connect RAB.
C
        ISTAT = SYS$CONNECT(RAB(IUNIT))
      ENDIF
C
      IF (.NOT.ISTAT) THEN
        CALL LIB$SYS_GETMSG(ISTAT,L,MSG)
        WRITE (*,'(/1X,2A,I4,2A)')
     &  ATBUTE(1:MIN(30,LENSTR(ATBUTE))),' file open error on unit',
     &  IUNIT,',   Logical name: ',NAME(:MIN(50,LENSTR(NAME)))
        WRITE(*,FMT='(//,1X,A)') MSG(1:MIN(130,L))
C
C==== Maybe this should be a soft fail: IUNIT = -1 ?
C
        CALL CCPERR(1,'QOPEN: FATAL ERROR')
        IUNIT = -1
        RETURN
      ENDIF
C
      WRITE (*,'(/1X,2A,I4,2A)')
     &ATBUTE(:MIN(30,LENSTR(ATBUTE))),' file opened on unit',IUNIT,
     &',   Logical name: ',NAME(:MIN(50,LENSTR(NAME)))
      WRITE (*,'(/,2A)')
     &      ' File name: ',FNAME(IUNIT)(:MIN(100,LENSTR(FNAME(IUNIT))))
      IF (NEW(IUNIT)) THEN
        WRITE (*,'(A,2I8/)')
     &  ' Initial & extend sizes in physical blocks =',
     &  FAB(IUNIT).FAB$L_ALQ,FAB(IUNIT).FAB$W_DEQ
      ELSE
        WRITE (*,'(A,I10/)') ' File size in physical blocks =',
     &  MAPXAB(IUNIT).XABFHC.XAB$L_EBK - 1 +
     &  (MAPXAB(IUNIT).XABFHC.XAB$W_FFB + NPBSZ - 1) / NPBSZ
      ENDIF
C
      LRECSZ(IUNIT) = NSIZE
      MRECSZ(IUNIT) = 0
      NRECSZ(IUNIT) = 0
      NVBLK(IUNIT) = 1
      NEL(IUNIT) = 1
      BUFACT(IUNIT) = .FALSE.
      WRTACT(IUNIT) = .FALSE.
C
C==== Initialise source machine integer and real types.
C
      MSIT(IUNIT)=0
      MSRT(IUNIT)=0
C
C==== Fill buffer with 0's.
C
      IF (NSIZE.LT.I2P15) THEN
        CALL LIB$MOVC5(0,,0,NSIZE,BUFFER(1,IUNIT))
      ELSE
C==== Convert buffer size to unsigned integer.
        CALL LIB$MOVC5(0,,0,NSIZE-I2P16,BUFFER(1,IUNIT))
      ENDIF
      RETURN
C
C
C#######################################################################
C
C
C*QCLOSE
      ENTRY QCLOSE(IUNIT)
C     ===================
C
C==== Maybe this should be a soft fail: IUNIT = -1 ?
C
      IF (LUNIT(IUNIT)) THEN
        CALL CCPERR(1,'QCLOSE error: File not open.')
        IUNIT = -1
        RETURN
      ENDIF
C
      ISTAT = 1
      IF (WRTACT(IUNIT) .AND. BUFACT(IUNIT)) THEN
        IF (MRECSZ(IUNIT).GT.0) THEN
          NOBLK = (LRECSZ(IUNIT) - 1) / NPBSZ
          KEL = NPBSZ * NOBLK
          KRECSZ = MRECSZ(IUNIT)
          IF (KRECSZ.LT.NRECSZ(IUNIT)) KRECSZ =
     &    MIN(NPBSZ * ((KRECSZ - 1) / NPBSZ + 1), NRECSZ(IUNIT))
          KRECSZ = KRECSZ - KEL
          IF (KRECSZ.GE.I2P15) KRECSZ = KRECSZ - I2P16
C
C==== Virtual block number, record buffer address and size.
C
          RAB(IUNIT).RAB$L_BKT = NVBLK(IUNIT) + NOBLK
          RAB(IUNIT).RAB$L_RBF = %LOC(BUFFER(KEL+1,IUNIT))
          RAB(IUNIT).RAB$W_RSZ = KRECSZ
C
C==== Write out the block.
C
          ISTAT = SYS$WRITE(RAB(IUNIT))
D         WRITE (*,'(A,5I8)') ' SYS$WRITE 1:',
D    &    IUNIT,NVBLK(IUNIT)+NOBLK,KEL+1,KRECSZ,ISTAT
          IF (.NOT.ISTAT) THEN
            CALL LIB$SYS_GETMSG(ISTAT,L,MSG)
            WRITE (*,'(/1X,A/)') MSG(:MIN(130,L))
          ENDIF
        ENDIF
      ENDIF
      LUNIT(IUNIT) = .TRUE.                             !FREE THIS UNIT #
C
C==== Close the file.
C
      ISAVE = ISTAT
      ISTAT = SYS$CLOSE(FAB(IUNIT))
      IF (.NOT.ISTAT) THEN
        CALL LIB$SYS_GETMSG(ISTAT,L,MSG)
        WRITE (*,'(/1X,A/)') MSG(:MIN(130,L))
      ENDIF
C      IF (.NOT.ISAVE .OR. .NOT.ISTAT)
C     &CALL CCPERR(1,'QCLOSE: FATAL ERROR.')
      RETURN
C
C
C#######################################################################
C
C
C*QMODE
      ENTRY QMODE(IUNIT,MODE,MCHITM)
C     ==============================
C
C==== Maybe this should be a soft fail: IUNIT = -1 ?
C
      IF (LUNIT(IUNIT)) THEN
        CALL CCPERR(1,'QMODE error: File not open.')
        IUNIT = -1
        RETURN
      ENDIF
C Check if recognized mode
      IF (MODE.GE.0 .AND. MODE.LE.MAXMO) THEN
        MCHITM=MAX(MODES(MODE),1)
      ELSE
        MCHITM=1
      ENDIF
      NMODE(IUNIT)=MODE
      NCHITM(IUNIT)=MCHITM
      RETURN
C
C
C#######################################################################
C
C
C*QREAD
      ENTRY QREAD(IUNIT,ARRAY,NBYTES,IER)
      ENTRY QREADR(IUNIT,ARRAY,NBYTES,IER)
      ENTRY QREADI(IUNIT,ARRAY,NBYTES,IER)
      ENTRY QREADQ(IUNIT,ARRAY,NBYTES,IER)
C     ===================================
C
C==== Maybe this should be a soft fail: IUNIT = -1 ?
C
 45   IF (IUNIT.GT.MAXNC .OR. IUNIT.LT.1)
     +     CALL CCPERR (1,'QREAD: bad stream number')
      IF (LUNIT(IUNIT)) THEN
        CALL CCPERR(1,'QREAD error: File not open.')
        IUNIT = -1
        RETURN
      ENDIF
      NTOMOV = NBYTES * NCHITM(IUNIT)
      IER = 0
      INDEX = 1
      NEW(IUNIT) = .FALSE.
D     WRITE (*,'(A,2I8)') ' QREAD     1:',NTOMOV,NEL(IUNIT)
C
      IF (.NOT.BUFACT(IUNIT)) THEN
C
C==== Virtual block number, user buffer address and size.
C
        RAB(IUNIT).RAB$L_BKT = NVBLK(IUNIT)
        RAB(IUNIT).RAB$L_UBF = %LOC(BUFFER(1,IUNIT))
        RAB(IUNIT).RAB$W_USZ = NSIZE
C
C==== Read in the block.
C
        ISTAT = SYS$READ(RAB(IUNIT))
        NRECSZ(IUNIT) = RAB(IUNIT).RAB$W_RSZ
        IF (NRECSZ(IUNIT).LT.0) NRECSZ(IUNIT) = NRECSZ(IUNIT) + I2P16
D       WRITE (*,'(A,6I8)') ' SYS$READ  1:',
D    &  IUNIT,NVBLK(IUNIT),1,NSIZE,ISTAT,NRECSZ(IUNIT)
        IF (.NOT.ISTAT) GOTO 60
C        WRITE (*,'(1X,16F8.0)') (BUF(I,IUNIT),I=1,NSIZE/4)
        BUFACT(IUNIT) = .TRUE.
      ENDIF
C
50    NLEFT = NRECSZ(IUNIT) + 1 - NEL(IUNIT)
      IF (NTOMOV.LE.NLEFT) THEN
D       WRITE (*,'(A,3I8)') ' QREAD     2:',NTOMOV,NEL(IUNIT),INDEX
C
C==== Move NTOMOV bytes from record BUFFER to user's ARRAY.
C
        IF (NTOMOV.LT.I2P15) THEN
          CALL LIB$MOVC3(NTOMOV,BUFFER(NEL(IUNIT),IUNIT),ARRAY(INDEX))
        ELSE
C==== Convert buffer size to unsigned integer.
          CALL LIB$MOVC3(NTOMOV-I2P16,BUFFER(NEL(IUNIT),IUNIT),
     &    ARRAY(INDEX))
        ENDIF
        NEL(IUNIT) = NEL(IUNIT) + NTOMOV
        GOTO 100
      ELSE
        IF (NLEFT.GT.0) THEN
D         WRITE (*,'(A,3I8)') ' QREAD     3:',NLEFT,NEL(IUNIT),INDEX
C
C==== Move NLEFT bytes from record BUFFER to user's ARRAY.
C
          IF (NLEFT.LT.I2P15) THEN
            CALL LIB$MOVC3(NLEFT,BUFFER(NEL(IUNIT),IUNIT),ARRAY(INDEX))
          ELSE
C==== Convert buffer size to unsigned integer.
            CALL LIB$MOVC3(NLEFT-I2P16,BUFFER(NEL(IUNIT),IUNIT),
     &      ARRAY(INDEX))
          ENDIF
          NTOMOV = NTOMOV - NLEFT
          INDEX = INDEX + NLEFT
        ENDIF
        IF (WRTACT(IUNIT)) THEN
          IF (MRECSZ(IUNIT).GT.0) THEN
            NOBLK = (LRECSZ(IUNIT) - 1) / NPBSZ
            KEL = NPBSZ * NOBLK
            KRECSZ = MRECSZ(IUNIT)
            IF (KRECSZ.LT.NRECSZ(IUNIT)) KRECSZ =
     &      MIN(NPBSZ * ((KRECSZ - 1) / NPBSZ + 1), NRECSZ(IUNIT))
            KRECSZ = KRECSZ - KEL
            IF (KRECSZ.GE.I2P15) KRECSZ = KRECSZ - I2P16
            RAB(IUNIT).RAB$L_BKT = NVBLK(IUNIT) + NOBLK
            RAB(IUNIT).RAB$L_RBF = %LOC(BUFFER(KEL+1,IUNIT))
            RAB(IUNIT).RAB$W_RSZ = KRECSZ
            ISTAT = SYS$WRITE(RAB(IUNIT))
D           WRITE (*,'(A,5I8)') ' SYS$WRITE 2:',
D    &      IUNIT,NVBLK(IUNIT)+NOBLK,KEL+1,KRECSZ,ISTAT
            IF (.NOT.ISTAT) THEN
              CALL LIB$SYS_GETMSG(ISTAT,L,MSG)
              WRITE (*,'(/1X,A/)') MSG(:MIN(130,L))
              CALL CCPERR(1,'QREAD: FATAL ERROR.')
              IER = -1
              RETURN
            ENDIF
            LRECSZ(IUNIT) = NSIZE
            MRECSZ(IUNIT) = 0
            NRECSZ(IUNIT) = 0
          ENDIF
          WRTACT(IUNIT) = .FALSE.
        ENDIF
        NVBLK(IUNIT) = NVBLK(IUNIT) + NRBLK
        NEL(IUNIT) = 1
        IF (NTOMOV.GE.NSIZE) THEN
          NDO = MIN(NTOMOV/NPBSZ, MAXDO)
          NBY = NDO*NPBSZ
C
C==== Convert buffer size to unsigned integer.
C
          IF (NBY.LT.I2P15) THEN
            IBY = NBY
          ELSE
            IBY = NBY - I2P16
          ENDIF
          RAB(IUNIT).RAB$L_BKT = NVBLK(IUNIT)
          RAB(IUNIT).RAB$L_UBF = %LOC(ARRAY(INDEX))
          RAB(IUNIT).RAB$W_USZ = IBY
          ISTAT = SYS$READ(RAB(IUNIT))
D         NRECSZ(IUNIT) = RAB(IUNIT).RAB$W_RSZ
D         IF (NRECSZ(IUNIT).LT.0) NRECSZ(IUNIT) = NRECSZ(IUNIT) + I2P16
D         WRITE (*,'(A,6I8)') ' SYS$READ  2:',
D    &    IUNIT,NVBLK(IUNIT),INDEX,IBY,ISTAT,NRECSZ(IUNIT)
          IF (.NOT.ISTAT) GOTO 60
          NTOMOV = NTOMOV - NBY
          NVBLK(IUNIT) = NVBLK(IUNIT) + NDO
          IF (NTOMOV.EQ.0) THEN
            BUFACT(IUNIT) = .FALSE.
            GOTO 100
          ENDIF
          INDEX = INDEX + NBY
        ENDIF
        RAB(IUNIT).RAB$L_BKT = NVBLK(IUNIT)
        RAB(IUNIT).RAB$L_UBF = %LOC(BUFFER(1,IUNIT))
        RAB(IUNIT).RAB$W_USZ = NSIZE
        ISTAT = SYS$READ(RAB(IUNIT))
        NRECSZ(IUNIT) = RAB(IUNIT).RAB$W_RSZ
        IF (NRECSZ(IUNIT).LT.0) NRECSZ(IUNIT) = NRECSZ(IUNIT) + I2P16
D       WRITE (*,'(A,6I8)') ' SYS$READ  3:',
D    &  IUNIT,NVBLK(IUNIT),1,NSIZE,ISTAT,NRECSZ(IUNIT)
        IF (.NOT.ISTAT) GOTO 60
C        WRITE (*,'(1X,16F8.0)') (BUF(I,IUNIT),I=1,NSIZE/4)
      ENDIF
      GOTO 50
C
60    IF (ISTAT.NE.RMS$_EOF) THEN
        CALL LIB$SYS_GETMSG(ISTAT,L,MSG)
        WRITE (*,'(/1X,A/)') MSG(:MIN(130,L))
        CALL CCPERR(1,'QREAD: FATAL ERROR.')
      ENDIF
      BUFACT(IUNIT) = .FALSE.
      IER = NBYTES - NTOMOV/NCHITM(IUNIT)
      IF (IER.EQ.0) IER = -1
C
100	IF (IER.EQ.0) THEN
	  NTOMOV=NBYTES*NCHITM(IUNIT)
	ELSE
	  NTOMOV=MAX(IER,0)*NCHITM(IUNIT)
	ENDIF
C	WRITE(6,*)'NMODE =',NMODE(IUNIT)
C
C====	Check if source machine integer conversion required.
C
	IF (MSIT(IUNIT).NE.0) THEN
C
C====	INTEGER*2 or *4 mode ?
C
	  IF (NMODE(IUNIT).EQ.1 .OR. NMODE(IUNIT).EQ.3) THEN
	    CALL QSWAP1(NTOMOV,ARRAY)
	  ELSEIF (NMODE(IUNIT).EQ.6) THEN
	    CALL QSWAP3(NTOMOV,ARRAY)
	  ENDIF
	ENDIF
C
C====	Check if source machine real conversion required.
C
	IF (MSRT(IUNIT).NE.0) THEN
	  IF (NMODE(IUNIT).EQ.2 .OR. NMODE(IUNIT).EQ.4) THEN
	    IF (MTRT.EQ.2) THEN
C
C====	Target machine is VAX; check source machine type for byte swap.
C
	      IF (MSRT(IUNIT).NE.4) CALL QSWAP3(NTOMOV,ARRAY)
C
C====	Check for Convex native.  If IEEE convert from IEEE.
C
	      IF (MSRT(IUNIT).NE.5) CALL QFIEEE(NTOMOV/4,ARRAY)
C
C====	Swap to VAX byte order.
C
	      CALL QSWAP2(NTOMOV,ARRAY)
C
	    ELSEIF (MTRT.EQ.4) THEN
C
C====	Target machine is L-E/IEEE; check source machine type for byte swap.
C
	      IF (MSRT(IUNIT).EQ.2) THEN
	        CALL QSWAP2(NTOMOV,ARRAY)
	      ELSE
	        CALL QSWAP3(NTOMOV,ARRAY)
	      ENDIF
C
C====	If not IEEE convert to IEEE.
C
	      IF (MSRT(IUNIT).NE.1) CALL QTIEEE(NTOMOV/4,ARRAY)
C
	    ELSEIF (MTRT.EQ.1) THEN
C
C====	Target machine is B-E/IEEE; check source machine type for byte swap.
C
	      IF (MSRT(IUNIT).EQ.2) THEN
	        CALL QSWAP1(NTOMOV,ARRAY)
	      ELSEIF (MSRT(IUNIT).EQ.4) THEN
	        CALL QSWAP3(NTOMOV,ARRAY)
	      ENDIF
C
C====	If not IEEE convert to IEEE.
C
	      IF (MSRT(IUNIT).NE.4) CALL QTIEEE(NTOMOV/4,ARRAY)
	    ELSE
	      CALL CCPERR(1,'*** QREAD: Bad machine REAL type.')
	    ENDIF
	  ENDIF
	ENDIF
      RETURN
C
C
C#######################################################################
C
C
C*QWRITE
      ENTRY QWRITE(IUNIT,ARRAY,NBYTES)
      ENTRY QWRITI(IUNIT,ARRAY,NBYTES)
      ENTRY QWRITR(IUNIT,ARRAY,NBYTES)
      ENTRY QWRITQ(IUNIT,ARRAY,NBYTES)
C     ================================
C
C==== Maybe this should be a soft fail: IUNIT = -1 ?
C
      IF (IUNIT.GT.MAXNC .OR. IUNIT.LT.1)
     +     CALL CCPERR (1,'QREAD: bad stream number')
      IF (LUNIT(IUNIT)) THEN
        CALL CCPERR(1,'QWRITE error: File not open.')
        IUNIT = -1
        RETURN
      ENDIF
      INDEX = 1
      NTOMOV = NBYTES * NCHITM(IUNIT)
D     WRITE (*,'(A,2I8)') ' QWRITE    1:',NTOMOV,NEL(IUNIT)
      IF (NEL(IUNIT).GT.NSIZE) THEN
        IF (WRTACT(IUNIT)) THEN
          IF (MRECSZ(IUNIT).GT.0) THEN
            NOBLK = (LRECSZ(IUNIT) - 1) / NPBSZ
            KEL = NPBSZ * NOBLK
            KRECSZ = NPBSZ * ((MRECSZ(IUNIT) - 1) / NPBSZ + 1) - KEL
            IF (KRECSZ.GE.I2P15) KRECSZ = KRECSZ - I2P16
            RAB(IUNIT).RAB$L_BKT = NVBLK(IUNIT) + NOBLK
            RAB(IUNIT).RAB$L_RBF = %LOC(BUFFER(KEL+1,IUNIT))
            RAB(IUNIT).RAB$W_RSZ = KRECSZ
            ISTAT = SYS$WRITE(RAB(IUNIT))
D           WRITE (*,'(A,5I8)') ' SYS$WRITE 3:',
D    &      IUNIT,NVBLK(IUNIT)+NOBLK,KEL+1,KRECSZ,ISTAT
            IF (.NOT.ISTAT) GOTO 80
            LRECSZ(IUNIT) = NSIZE
            MRECSZ(IUNIT) = 0
            NRECSZ(IUNIT) = 0
          ENDIF
        ENDIF
        NEL(IUNIT) = 1
        NVBLK(IUNIT) = NVBLK(IUNIT) + NRBLK
        BUFACT(IUNIT) = .FALSE.
      ENDIF
      WRTACT(IUNIT) = .TRUE.
C
70    NLEFT = NSIZE1 - NEL(IUNIT)
      IF (.NOT.BUFACT(IUNIT)) THEN
        IF (NEL(IUNIT).GT.1 .OR. NTOMOV.LT.NSIZE) THEN
          IF (.NOT.NEW(IUNIT)) THEN
            RAB(IUNIT).RAB$L_BKT = NVBLK(IUNIT)
            RAB(IUNIT).RAB$L_UBF = %LOC(BUFFER(1,IUNIT))
            RAB(IUNIT).RAB$W_USZ = NSIZE
            ISTAT = SYS$READ(RAB(IUNIT))
            NRECSZ(IUNIT) = RAB(IUNIT).RAB$W_RSZ
            IF (NRECSZ(IUNIT).LT.0)
     &      NRECSZ(IUNIT) = NRECSZ(IUNIT) + I2P16
D           WRITE (*,'(A,6I8)') ' SYS$READ  4:',
D    &      IUNIT,NVBLK(IUNIT),1,NSIZE,ISTAT,NRECSZ(IUNIT)
C
C     ERROR=BLANK FILE
C
C==== Fill buffer with 0's.
C
            IF (.NOT.ISTAT) THEN
              IF (NSIZE.LT.I2P15) THEN
                CALL LIB$MOVC5(0,,0,NSIZE,BUFFER(1,IUNIT))
              ELSE
C==== Convert buffer size to unsigned integer.
                CALL LIB$MOVC5(0,,0,NSIZE-I2P16,BUFFER(1,IUNIT))
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        BUFACT(IUNIT) = .TRUE.
      ENDIF
C
      IF (NTOMOV.LE.NLEFT) THEN
D       WRITE (*,'(A,3I8)') ' QWRITE    2:',NTOMOV,NEL(IUNIT),INDEX
C
C==== Move NTOMOV bytes from user's ARRAY to record BUFFER.
C
        IF (NTOMOV.LT.I2P15) THEN
          CALL LIB$MOVC3(NTOMOV,ARRAY(INDEX),BUFFER(NEL(IUNIT),IUNIT))
        ELSE
C==== Convert buffer size to unsigned integer.
          CALL LIB$MOVC3(NTOMOV-I2P16,ARRAY(INDEX),
     &    BUFFER(NEL(IUNIT),IUNIT))
        ENDIF
        LRECSZ(IUNIT)=MIN(LRECSZ(IUNIT),NEL(IUNIT))
        NEL(IUNIT) = NEL(IUNIT) + NTOMOV
        MRECSZ(IUNIT)=MAX(MRECSZ(IUNIT),NEL(IUNIT)-1)
        NRECSZ(IUNIT)=MAX(NRECSZ(IUNIT),MRECSZ(IUNIT))
        RETURN
      ELSE
D       WRITE (*,'(A,3I8)') ' QWRITE    3:',NLEFT,NEL(IUNIT),INDEX
C
C==== Move NLEFT bytes from user's ARRAY to record BUFFER.
C
        IF (NLEFT.LT.I2P15) THEN
          CALL LIB$MOVC3(NLEFT,ARRAY(INDEX),BUFFER(NEL(IUNIT),IUNIT))
        ELSE
C==== Convert buffer size to unsigned integer.
          CALL LIB$MOVC3(NLEFT-I2P16,ARRAY(INDEX),
     &    BUFFER(NEL(IUNIT),IUNIT))
        ENDIF
        NOBLK = (MIN(LRECSZ(IUNIT),NEL(IUNIT)) - 1) / NPBSZ
        KEL = NPBSZ * NOBLK
        KRECSZ = NPBSZ * NRBLK - KEL
        IF (KRECSZ.GE.I2P15) KRECSZ = KRECSZ - I2P16
        RAB(IUNIT).RAB$L_BKT = NVBLK(IUNIT) + NOBLK
        RAB(IUNIT).RAB$L_RBF = %LOC(BUFFER(KEL+1,IUNIT))
        RAB(IUNIT).RAB$W_RSZ = KRECSZ
        ISTAT = SYS$WRITE(RAB(IUNIT))
D       WRITE (*,'(A,5I8)') ' SYS$WRITE 4:',
D    &  IUNIT,NVBLK(IUNIT)+NOBLK,KEL+1,KRECSZ,ISTAT
        IF (.NOT.ISTAT) GOTO 80
        LRECSZ(IUNIT) = NSIZE
        MRECSZ(IUNIT) = 0
        NRECSZ(IUNIT) = 0
        NVBLK(IUNIT) = NVBLK(IUNIT) + NRBLK
        NEL(IUNIT) = 1
        NTOMOV = NTOMOV - NLEFT
        INDEX = INDEX + NLEFT
        BUFACT(IUNIT) = .FALSE.
C
        IF (NTOMOV.GT.NSIZE) THEN
          NDO = MIN(NTOMOV/NPBSZ, MAXDO)
          NBY = NDO*NPBSZ
          IF (NBY.LT.I2P15) THEN
            IBY = NBY
          ELSE
            IBY = NBY - I2P16
          ENDIF
          RAB(IUNIT).RAB$L_BKT = NVBLK(IUNIT)
          RAB(IUNIT).RAB$L_RBF = %LOC(ARRAY(INDEX))
          RAB(IUNIT).RAB$W_RSZ = IBY
          ISTAT = SYS$WRITE(RAB(IUNIT))
D         WRITE (*,'(A,5I8)') ' SYS$WRITE 5:',
D    &    IUNIT,NVBLK(IUNIT),1,IBY,ISTAT
          IF (.NOT.ISTAT) GOTO 80
          NVBLK(IUNIT) = NVBLK(IUNIT) + NDO
          INDEX = INDEX + NBY
          NTOMOV = NTOMOV - NBY
        ENDIF
      ENDIF
      GOTO 70
C
80    CALL LIB$SYS_GETMSG(ISTAT,L,MSG)
      WRITE (*,'(/1X,A/)') MSG(:MIN(130,L))
      CALL CCPERR(1,'QWRITE: FATAL ERROR.')
      RETURN
C
C
C#######################################################################
C
C
C*QBACK
      ENTRY QBACK(IUNIT,LRECL)
C     ========================
C
C==== Maybe this should be a soft fail: IUNIT = -1 ?
C
      IF (LUNIT(IUNIT)) THEN
        CALL CCPERR(1,'QBACK error: File not open.')
        IUNIT = -1
        RETURN
      ENDIF
      JEL = NEL(IUNIT) - LRECL * NCHITM(IUNIT)
      IVBLK = (JEL - NSIZE)/NSIZE
      JEL = JEL - IVBLK*NSIZE
      IVBLK = NVBLK(IUNIT) + IVBLK*NRBLK
      IF (IVBLK.LT.1) THEN
        IVBLK = 1
        JEL = 1
      ENDIF
      GOTO 90
C
C
C#######################################################################
C
C
C*QSKIP
      ENTRY QSKIP(IUNIT,LRECL)
C     ========================
C
C==== Maybe this should be a soft fail: IUNIT = -1 ?
C
      IF (LUNIT(IUNIT)) THEN
        IF (LUNIT(IUNIT)) CALL CCPERR(1,'QSKIP error: File not open.')
        IUNIT = -1
        RETURN
      ENDIF
      JEL = NEL(IUNIT) + LRECL*NCHITM(IUNIT) - 1
      IVBLK = JEL/NSIZE
      JEL = JEL - IVBLK*NSIZE + 1
      IVBLK = NVBLK(IUNIT) + IVBLK*NRBLK
      GOTO 90
C
C
C#######################################################################
C
C
C*QSEEK
      ENTRY QSEEK(IUNIT,IREC,IEL,LRECL)
C     =================================
C
C==== Maybe this should be a soft fail: IUNIT = -1 ?
C
      IF (LUNIT(IUNIT)) THEN
        CALL CCPERR(1,'QSEEK error: File not open.')
        IUNIT = -1
        RETURN
      ENDIF
      JEL = MAX(((IREC - 1)*LRECL + IEL - 1)*NCHITM(IUNIT), 0)
      IVBLK = JEL/NSIZE
      JEL = JEL - IVBLK*NSIZE + 1
      IVBLK = IVBLK*NRBLK + 1
C
90    IF (IVBLK.NE.NVBLK(IUNIT)) THEN
        IF (WRTACT(IUNIT) .AND. NEL(IUNIT).NE.1) THEN
          IF (MRECSZ(IUNIT).GT.0) THEN
            NOBLK = (LRECSZ(IUNIT) - 1) / NPBSZ
            KEL = NPBSZ * NOBLK
            KRECSZ = MRECSZ(IUNIT)
            IF (KRECSZ.LT.NRECSZ(IUNIT)) KRECSZ =
     &      MIN(NPBSZ * ((KRECSZ - 1) / NPBSZ + 1), NRECSZ(IUNIT))
            KRECSZ = KRECSZ - KEL
            IF (KRECSZ.GE.I2P15) KRECSZ = KRECSZ - I2P16
            RAB(IUNIT).RAB$L_BKT = NVBLK(IUNIT) + NOBLK
            RAB(IUNIT).RAB$L_RBF = %LOC(BUFFER(KEL+1,IUNIT))
            RAB(IUNIT).RAB$W_RSZ = KRECSZ
            ISTAT = SYS$WRITE(RAB(IUNIT))
D           WRITE (*,'(A,5I8)') ' SYS$WRITE 6:',
D    &      IUNIT,NVBLK(IUNIT)+NOBLK,KEL+1,KRECSZ,ISTAT
C
            IF (.NOT.ISTAT) THEN
              CALL LIB$SYS_GETMSG(ISTAT,L,MSG)
              WRITE (*,'(/1X,A/)') MSG(:MIN(130,L))
              CALL CCPERR(1,'QSEEK: FATAL ERROR -- MAYBE CORRUPT FILE.')
              RETURN
            ENDIF
            LRECSZ(IUNIT) = NSIZE
            MRECSZ(IUNIT) = 0
            NRECSZ(IUNIT) = 0
          ENDIF
          WRTACT(IUNIT) = .FALSE.
        ENDIF
        NVBLK(IUNIT) = IVBLK
        BUFACT(IUNIT) = .FALSE.
        NEW(IUNIT) = .FALSE.
      ENDIF
C
      NEL(IUNIT) = JEL
D     WRITE (*,'(A,8X,I8)') ' QSEEK      :',NEL(IUNIT)
      RETURN
C
C
C#######################################################################
C
C
C*QLOCATE
      ENTRY QLOCATE(IUNIT,LOCATION)
C     =============================
C
C==== Maybe this should be a soft fail: IUNIT = -1 ?
C
      IF (LUNIT(IUNIT)) THEN
        CALL CCPERR(1,'QLOCATE error: File not open.')
        IUNIT = -1
        RETURN
      ENDIF
      LOCATION = (NEL(IUNIT)+(NVBLK(IUNIT)-1)*NPBSZ+NCHITM(IUNIT)-1)/
     &           NCHITM(IUNIT)
      RETURN
C
C
C#######################################################################
C
C
C*QQINQ
      ENTRY QQINQ(IUNIT,NAME,FILNAM,LENGTH)
C     =====================================
C
C RETURNS FILE NAME AND LENGTH 
C IF NO NAME RETURN NAME AS BLANK, IF NO LENGTH RETURN LENGTH AS -1
C
      IF (LUNIT(IUNIT)) THEN
        INQUIRE (FILE=NAME,NAME=FILNAM)
        LENGTH = -1
      ELSE
        FILNAM = FNAME(IUNIT)
        LENGTH = NPBSZ * (MAPXAB(IUNIT).XABFHC.XAB$L_EBK - 1) +
     &  MAPXAB(IUNIT).XABFHC.XAB$W_FFB
      ENDIF
      RETURN
C
C
C#######################################################################
C
C
	ENTRY QSETMS(IUNIT,MS,IRET)
C	===========================
C
C====	Get machine integer and real types.
C
	CALL QIRTYP(MTIT,MTRT)
C
C====	Get file integer and real types.
C
	MSIT(IUNIT)=MS(2)/16
	IF (MSIT(IUNIT).NE.0 .AND. MSIT(IUNIT).NE.1 .AND.
     &	MSIT(IUNIT).NE.4) THEN
	  WRITE (6,*) '*** BAD FILE INTEGER TYPE, ASSUME NATIVE ***'
	  MSIT(IUNIT)=0
	ENDIF
	MSRT(IUNIT)=MS(1)/16
	IF ((MSRT(IUNIT).LT.0 .OR. MSRT(IUNIT).GT.2) .AND.
     &	MSRT(IUNIT).NE.4 .AND. MSRT(IUNIT).NE.5) THEN
	  WRITE (6,*) '*** BAD FILE REAL TYPE, ASSUME NATIVE ***'
	  MSIT(IUNIT)=0
	ENDIF
	IF (MSIT(IUNIT).EQ.0) MSRT(IUNIT)=0
	IRET=16*MSIT(IUNIT)+MSRT(IUNIT)
	IF (MSIT(IUNIT).EQ.MTIT) MSIT(IUNIT)=0
	IF (MSRT(IUNIT).EQ.MTRT) MSRT(IUNIT)=0
C	WRITE (6,*) 'MSIT, MTIT, MSRT, MTRT =',
C     &	MSIT(IUNIT),MTIT,MSRT(IUNIT),MTRT
	END
C
C
C#######################################################################
C
C
      SUBROUTINE QQOPEN(IUNIT,NAME,ISTAT)
C     ===================================
C
C Open file with name NAME, returns IUNIT as unit number
C
C
      IMPLICIT       NONE
      CHARACTER*(*)  NAME
      INTEGER        IUNIT, ISTAT, NCH
C
C
      GOTO (10,20,30,40,50),ISTAT
10    CALL QOPEN(IUNIT,NAME,'UNKNOWN')
      GOTO 100
20    CALL QOPEN(IUNIT,NAME,'SCRATCH')
      GOTO 100
30    CALL QOPEN(IUNIT,NAME,'OLD')
      GOTO 100
40    CALL QOPEN(IUNIT,NAME,'NEW')
      GOTO 100
50    CALL QOPEN(IUNIT,NAME,'RO')
100   RETURN
      END
C
C
C#######################################################################
C
C
	SUBROUTINE QIRTYP(MI,MR)
C	========================
C====	Integer/Real types; assumes equal endedness of integers and reals.
	IMPLICIT NONE
	INTEGER*2 I(2)
	INTEGER J,K,MI,MR
	REAL A,B,R
	PARAMETER (A=513./256., B=513./1024.)
	EQUIVALENCE (I,J), (K,R)
	DATA I/1,0/, K/1073758208/
	IF (J.EQ.65536) THEN
C====	Big-endian.
	  MI=1
	  IF (R.EQ.A) THEN
C====	IEEE.
	    MR=1
	  ELSEIF (R.EQ.B) THEN
C====	Convex native.
	    MR=5
	  ELSE
	    MR=0
	  ENDIF
	ELSEIF (J.EQ.1) THEN
C====	Little-endian.
	  MI=4
	  IF (R.EQ.A) THEN
C====	IEEE.
	    MR=4
	  ELSEIF (R.EQ.B) THEN
C====	VAX native.
	    MR=2
	  ELSE
	    MR=0
	  ENDIF
	ELSE
	  MI=0
	  MR=0
	ENDIF
	END
C
C
C#######################################################################
C
C
	SUBROUTINE QSWAP1(N,A)
C	======================
C
C====	Swap bytes N/N+1 for N = 0,2,4...
C
	IMPLICIT NONE
	INTEGER I,N
	BYTE A(N),S
	DO 10 I=1,N,2
	  S=A(I)
	  A(I)=A(I+1)
10	  A(I+1)=S
	END 
C
C
C#######################################################################
C
C
	SUBROUTINE QSWAP2(N,A)
C	======================
C
C====	Swap bytes N/N+2 and N+1/N+3 for N = 0,4,8...
C
	IMPLICIT NONE
	INTEGER I,N
	BYTE A(N),S
	DO 10 I=1,N,4
	  S=A(I)
	  A(I)=A(I+2)
	  A(I+2)=S
	  S=A(I+1)
	  A(I+1)=A(I+3)
10	  A(I+3)=S
	END
C
C
C#######################################################################
C
C
	SUBROUTINE QSWAP3(N,A)
C	======================
C
C====	Swap bytes N/N+3 and N+1/N+2 for N = 0,4,8...
C
	IMPLICIT NONE
	INTEGER I,N
	BYTE A(N),S
	DO 10 I=1,N,4
	  S=A(I)
	  A(I)=A(I+3)
	  A(I+3)=S
	  S=A(I+1)
	  A(I+1)=A(I+2)
10	  A(I+2)=S
	END
C
C
C#######################################################################
C
C
	SUBROUTINE QFIEEE(N,A)
C	======================
C
C====	Convert from IEEE to VAX native, assuming IEEE byte order.
C====	NaN (not a number) becomes ROP (reserved operand).
C
	IMPLICIT NONE
	INTEGER IEXP,MANT,MDN1,MDN2,MEXP,MNAN,ROP
	PARAMETER (IEXP='01000000'X,MANT='003FFFFF'X,MDN1='00400000'X,
     &	           MDN2='00200000'X,MEXP='7E800000'X,MNAN='7F800000'X,
     &	            ROP='80000000'X)
	INTEGER E,I,N
	INTEGER A(N)
C
	DO 10 I=1,N
	  E=A(I).AND.MNAN
	  IF (E.GT.0) THEN
	    IF (E.LE.MEXP) THEN
	      A(I)=A(I)+IEXP
	    ELSEIF (E.LT.MNAN) THEN
	      A(I)=ROP
	    ELSE
	      A(I)=ROP.OR.(A(I).AND..NOT.MNAN)
	    ENDIF
	  ELSEIF ((A(I).AND.MDN1).NE.0) THEN
	    A(I)=(A(I).AND.ROP).OR.2*(A(I).AND.MANT).OR.IEXP
	  ELSEIF ((A(I).AND.MDN2).NE.0) THEN
	    A(I)=(A(I).AND.ROP).OR.4*(A(I).AND.MANT)
	  ELSE
	    A(I)=0
	  ENDIF
10	CONTINUE
	END
C
C
C#######################################################################
C
C
	SUBROUTINE QTIEEE(N,A)
C	======================
C
C====	Convert to IEEE from VAX/Convex native, assuming IEEE byte order.
C====	ROP (reserved operand) becomes NaN (not a number).
C
	IMPLICIT NONE
	INTEGER IEXP,MDN1,MDN2,MNAN
	PARAMETER (IEXP='01000000'X,MDN1='00400000'X,MDN2='00200000'X,
     &	           MNAN='7F800000'X)
	INTEGER E,I,N
	INTEGER A(N)
C
	DO 10 I=1,N
	  E=A(I).AND.MNAN
	  IF (E.GT.IEXP) THEN
	    A(I)=A(I)-IEXP
	  ELSEIF (E.EQ.0) THEN
	    IF (A(I).LT.0) THEN
	      A(I)=A(I).OR.MNAN
	    ELSE
	      A(I)=0
	    ENDIF
	  ELSE
	    IF (E.EQ.IEXP) THEN
	      A(I)=(A(I)/2.AND..NOT.MNAN).OR.MDN1
	    ELSE
	      A(I)=((A(I)+SIGN(2,A(I)))/4.AND..NOT.MNAN).OR.MDN2
	    ENDIF
	  ENDIF
10	CONTINUE
	END
C
C
C#######################################################################
C
C
	SUBROUTINE QRARCH(IUNIT,IEL,IRET)
C	=================================
C====	Read machine integer/real type from file.
	IMPLICIT NONE
	INTEGER I,IUNIT,IEL,IRET
	CHARACTER MSG*16
	BYTE MS(4)
C
	CALL UGTENV('CCP4_NATIVE',MSG)
	IF (MSG.NE.' ') THEN
	  IRET=0
	ELSE
	  CALL QMODE(IUNIT,0,IRET)
	  CALL QSEEK(IUNIT,1,4*IEL+1,0)
	  CALL QREAD(IUNIT,MS,2,IRET)
	  IF (IRET.NE.0) CALL CCPERR(1,'*** FATAL ERROR IN QRARCH ***')
C	  WRITE(6,*)'MS =',MS(1),MS(2)
	  CALL QSETMS(IUNIT,MS,IRET)
	ENDIF
	END
C
C
C#######################################################################
C
C
	SUBROUTINE QWARCH(IUNIT,IEL)
C	============================
C====	Write machine integer/real type to file.
	IMPLICIT NONE
	INTEGER IUNIT,IEL,IRET,MTIT,MTRT
	BYTE MS(4)
	DATA MS/4*0/
C
	CALL QIRTYP(MTIT,MTRT)
C	WRITE (6,*) 'MTIT, MTRT =',MTIT,MTRT
	MS(1)=17*MTRT
	MS(2)=16*MTIT+1
	CALL QMODE(IUNIT,0,IRET)
	CALL QSEEK(IUNIT,1,4*IEL+1,0)
	CALL QWRITE(IUNIT,MS,4)
	END
C
C
C#######################################################################
C
C
	SUBROUTINE QNAN(VALUE)
C====	Get undefined number value.
	IMPLICIT NONE
	INTEGER NAN,ROPC,ROPV
	PARAMETER (NAN='FFFA5A5A'X,ROPC='80000000'X,ROPV='8000'X)
	INTEGER INAN,MTIT,MTRT
	REAL RNAN, VALUE
	EQUIVALENCE(INAN,RNAN)
C
	CALL QIRTYP(MTIT,MTRT)
	IF (MTRT.EQ.1 .OR. MTRT.EQ.4) THEN
	  INAN=NAN
	ELSEIF (MTRT.EQ.5) THEN
	  INAN=ROPC
	ELSEIF (MTRT.EQ.2) THEN
	  INAN=ROPV
	ELSE
	  CALL CCPERR(1,'*** FATAL: QNAN machine type undefined.')
	ENDIF
	VALUE=RNAN
	END
C
C
C#######################################################################
C
C
	LOGICAL FUNCTION QISNAN(A)
C====	Test undefined number value.
	IMPLICIT NONE
	INTEGER MNAN,MROPC,MROPV,ROPC,ROPV
	PARAMETER (MNAN='7F800000'X,MROPC='FF800000'X,MROPV='FF80'X,
     &	ROPC='80000000'X,ROPV='8000'X)
	INTEGER N,MTIT,MTRT
	REAL A,R
	EQUIVALENCE(N,R)
C
	R=A
	CALL QIRTYP(MTIT,MTRT)
	IF (MTRT.EQ.2) THEN
	  QISNAN=(N.AND.MROPV).EQ.ROPV
	ELSEIF (MTRT.EQ.5) THEN
	  QISNAN=(N.AND.MROPC).EQ.ROPC
	ELSE IF (MTRT.EQ.1 .OR. MTRT.EQ.4) THEN
	  QISNAN=(N.AND.MNAN).EQ.MNAN
	ELSE
	  CALL CCPERR(1,'*** FATAL: QISNAN machine type undefined.')
	ENDIF
	END
C
C
C#######################################################################
C
C
      SUBROUTINE CCPBML (N, A)
C     Reset BIOMOL absence flags to zero in N elements of array A,
C     testing for Rops.  (We rely on the default optimisation to inline
C     the QISNAN, which otherwise will be called for each number in the
C     MTZ file.)
      INTEGER N, I
      REAL A (*)
      LOGICAL QISNAN
      EXTERNAL QISNAN
      DO I=1,N
        IF (.NOT.QISNAN (A(I))) THEN
          IF (A(I).LE.-0.99E10) A(I) = 0.0
        ENDIF
      ENDDO
      END
C
C
C#######################################################################
C
C
      SUBROUTINE CCPWRG (N, A, WRANGE)
C     This a an Rop-safe routine to update the column ranges needed by
C     mtzlib
      INTEGER I, N
      REAL A (*), WRANGE (2,N,*)
      LOGICAL QISNAN
      EXTERNAL QISNAN
      DO I=1,N
        IF (.NOT.QISNAN (A(I))) THEN
          IF (A(I).NE.-1E-10) THEN
            IF (A(I).LT.WRANGE(1,I,1)) WRANGE(1,I,1) = A(I)
            IF (A(I).GT.WRANGE(2,I,1)) WRANGE(2,I,1) = A(I)            
          ENDIF
        ENDIF
      ENDDO
      END
C     
C
C#######################################################################

C     Reading and writing CHARACTER BUFFER, len(buffer) bytes, mode 0.
C     (It may be possible to declare BUFFER as a structure and get the
C     data address directly, without calling STR$ANALYZE_SDESC.)

      SUBROUTINE QREADC (IUNIT,BUFFER,RESULT)
      PARAMETER (NCMAX=10)
      INTEGER IUNIT, NITEMS, RESULT, MODE, OLDMODE, NMODE(NCMAX), MITEM
      CHARACTER BUFFER*(*)
      COMMON /DISKIO_MODES_/ MAXNC, NMODE
      IF (IUNIT.GT.MAXNC .OR. IUNIT.LT.1)
     +     CALL CCPERR (1, 'QREADC: Bad unit number')
C     save the old mode and change to bytes
      OLDMODE = NMODE (IUNIT)
      CALL QMODE (IUNIT,0,MITEM)
      CALL STR$ANALYZE_SDESC (BUFFER, NITEMS, IPTR)
      CALL QREAD (IUNIT,%VAL(IPTR),NITEMS,RESULT)
C     restore mode
      CALL QMODE (IUNIT,OLDMODE,MITEM)
      END

      SUBROUTINE QWRITC (IUNIT,BUFFER,RESULT)
      PARAMETER (NCMAX=10)
      INTEGER IUNIT, NITEMS, RESULT, MODE, OLDMODE, NMODE(NCMAX), MITEM
      CHARACTER BUFFER*(*)
      COMMON /DISKIO_MODES_/ MAXNC, NMODE
      IF (IUNIT.GT.MAXNC .OR. IUNIT.LT.1)
     +     CALL CCPERR (1, 'QREADC: Bad unit number')
C     save the old mode and change to bytes
      OLDMODE = NMODE (IUNIT)
      CALL QMODE (IUNIT,0,MITEM)
      CALL STR$ANALYZE_SDESC (BUFFER, NITEMS, IPTR)
      CALL QWRITE (IUNIT,%VAL(IPTR),NITEMS,RESULT)
C     restore mode
      CALL QMODE (IUNIT,OLDMODE,MITEM)
      END

      SUBROUTINE CUNLINK (FILE)
C     Dummy version -- can't unlink VMS files
      CHARACTER *(*) FILE
      END
