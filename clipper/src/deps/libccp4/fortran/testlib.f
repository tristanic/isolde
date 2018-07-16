C
C     testlib.f: test program for Machine dependent routines
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
      PROGRAM MCTEST
C
C---- Test file for Machine dependent routines
C
C     .. Parameters ..
      INTEGER LBUF,LSTR
      PARAMETER (LBUF=500,LSTR=12)
C     ..
C     .. Local Scalars ..
      REAL ERSTR,SEC,NAN
      INTEGER I,IBYTE,IDAY,IER,ILENGTH,ILOOP,IMON,ISEC,ISEED,ISTAT,
     +        IYEAR,IYES,J,LDUM,LUN,LUNIN,LUNOUT,NREC
      CHARACTER ERRSTR*40,HANDLE* (LSTR),NAME1* (LSTR),NAME2* (LSTR),
     +          ENVNAM* (120),USRNAM* (LSTR),UDATE* (LSTR),
     +          USRTIM* (LSTR),REPLY* (LSTR),TSTNAM*(LSTR),FOO*3
C     ..
C     .. Local Arrays ..
      REAL BUFFER(LBUF)
C     ..
C     .. External Functions ..
      LOGICAL LITEND,VAXVMS, WINMVS, QISNAN
      REAL RANU
      INTEGER CCPNUN
      EXTERNAL LITEND,VAXVMS, WINMVS, QISNAN, RANU, CCPNUN
C     ..
C     .. External Subroutines ..
      EXTERNAL CCPERR,CCPFYP,NOCRLF,QCLOSE,QMODE,QOPEN,QQINQ,QREAD,
     +         QSEEK,QWRITE,UBYTES,UCPUTM,UGERR,UGTENV,UGTUID,
     +         UIDATE,UISATT,USTIME,UTIME,CCPRCS,CCPDPN, QNAN
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC NINT
C     ..
C     .. Data statements ..
      DATA LUNIN/5/,LUNOUT/6/,NAME1/'AAA.TST'/,NAME2/'ZZZ.TST'/,
     +     ISEED/0/,ILOOP/100/,TSTNAM/'TESTNAME'/
C     ..
C
C---- Initialise CPU timer
C
      SEC = -1.0
      CALL UCPUTM(SEC)
C
C---- Parse command line arguments, open printer stream and print
C     version 
C
      CALL CCPFYP
      I = 0
      CALL CCPDPN(LUNOUT,'PRINTER','PRINTER','F',0,I)
      IF (I.NE.0) CALL CCPERR(1,'Can''t open printer stream')
      CALL CCPRCS(6,'TESLIB','$Date$')
C
C---- Other initialisations
C
      I = 0
      CALL CCPDPN(10,NAME1,'NEW','F',0,I)
      IF (I.NE.0) CALL CCPERR(1,' Failed to open file '//NAME1)
      WRITE (10,FMT='(A)') ' I am the test file, delete me.'
      CLOSE (10,STATUS='DELETE')
C
C---- Start of Tests
C
      WRITE (LUNOUT,FMT=6000) ' Routine  Result      Comments'
C
C---- UGTENV
C
      CALL UGTENV('TESTENV',ENVNAM)
      WRITE (LUNOUT,FMT=6002) ' UGTENV',ENVNAM(1:12),
     .                                    'Get value of TESTENV'
C
C---- UGTUID
C
      CALL UGTUID(USRNAM)
      WRITE (LUNOUT,FMT=6002) ' UGTUID',USRNAM,'Get Users name'
C
C---- UIDATE
C
      CALL UIDATE(IMON,IDAY,IYEAR)
      WRITE (UDATE,FMT=6006) IDAY,IMON,IYEAR
      WRITE (LUNOUT,FMT=6002) ' UIDATE',UDATE,'Todays date'
C
C---- UTIME
C
      CALL UTIME(USRTIM)
      WRITE (LUNOUT,FMT=6002) ' UTIME ',USRTIM,'Current time'
C
C---- USTIME
C
      CALL USTIME(ISEC)
      WRITE (UDATE,FMT=6020) ISEC
      WRITE (LUNOUT,FMT=6002) ' USTIME',UDATE,'Absolute time (VMS=-1)'
C
C---- UISATT
C
      CALL UISATT(LUNOUT,IYES)
      IF (IYES.EQ.0) THEN
        REPLY = 'No'
      ELSE
        WRITE (REPLY,FMT=6004) 'Yes, ',LUNOUT
      END IF
      WRITE (LUNOUT,FMT=6002) ' UISATT',REPLY,
     +  'are we attached to at tty? Unit number?'
C
C---- VAXVMS
C
      REPLY = 'No'
      IF (VAXVMS()) REPLY = 'Yes'
      WRITE (LUNOUT,FMT=6002) ' VAXVMS',REPLY,'Is this VMS?'
C
C---- WINMVS
C
      REPLY = 'No'
      IF (WINMVS()) REPLY = 'Yes'
      WRITE (LUNOUT,FMT=6002) ' WINMVS',REPLY,'Is this Win NT et al?'
C
C---- UBYTES
C
      CALL UBYTES(IBYTE,HANDLE)
      WRITE (REPLY,FMT=6008) HANDLE,IBYTE
      WRITE (LUNOUT,FMT=6002) ' UBYTES',REPLY,
     +  'Get BYTE/WORD Handling and number of bytes per word'
C
C---- LITEND
C
      REPLY = 'Big'
      IF (LITEND(IDUM)) REPLY = 'Little'
      WRITE (LUNOUT,FMT=6002) ' LITEND',REPLY,'Big/Little end machine'
CCCC
CCCC---- URENAM
CCCC
CCC      CALL URENAM(NAME1,NAME2,ISTAT)
CCC      ERRSTR = 'OK'
CCC      IF (ISTAT.NE.0) CALL UGERR(ISTAT,ERSTR)
CCC      WRITE (LUNOUT,FMT=6002) ' URENAM',ERRSTR,'Check rename status'
CCC      CALL CUNLINK (NAME2)
C
C---- UGERR
C
      OPEN (21,FILE='TESTFILE',STATUS='OLD',IOSTAT=I)
      CALL UGERR(I,ERRSTR)
      WRITE (REPLY,FMT=6020) I
      WRITE (LUNOUT,FMT=6002) ' UGERR ',REPLY,ERRSTR
C
C---- UCPUTM
C
      SEC = 99.99
      CALL UCPUTM(SEC)
      WRITE (REPLY,FMT=6016) SEC
      WRITE (LUNOUT,FMT=6002) ' UCPUTM',REPLY,'Show elapsed CPU time'
C
C---- NOCRLF
C     
      CALL NOCRLF('NOCRLF')
      WRITE(LUNOUT,'(''+'',14X,A)') 'Should be on same line'
C
C --- CCPNUN
C
      I = CCPNUN ()
      WRITE (REPLY,FMT=6020) I
      WRITE (LUNOUT,FMT=6002) ' CCPNUN',REPLY,'Next free unit'
C
C --- QNAN/QISNAN (`magic numbers')
C
      CALL QNAN (NAN)
      IF ((.NOT.QISNAN (NAN)) .OR. QISNAN (1.0)) THEN
        WRITE (LUNOUT,'(/'' *** QNAN/QISNAN test failed''/)')
      ELSE
        WRITE (LUNOUT,'('' QNAN/QISNAN test OK'')') 
      ENDIF
C
C---- End of tests
C
      WRITE (LUNOUT,FMT=6000) ' Now test diskio routines'
C
C---- Now test the diskio stuff
C
      CALL QOPEN(LUN,'DISKIO','UNKNOWN')
      CALL QMODE(LUN,2,LDUM)
C
C---- Write a file of size LBUF x LBUF x WORDSIZE
C
      DO 20 I = 1,LBUF
        DO 10 J = 1,LBUF
10        BUFFER(J) = 0.
        BUFFER(I) = I
        CALL QWRITE(LUN,BUFFER,LBUF)
   20 CONTINUE
C
C---- Close the file
C
      CALL QCLOSE(LUN)
C
C---- reset the array buffer(*)
C
      DO 30 I = 1,LBUF
        BUFFER(I) = 0.0
   30 CONTINUE
C
C---- Now do some reads on the file just created
C
      CALL QOPEN(LUN,'DISKIO','OLD')
      CALL QMODE(LUN,2,LDUM)
C
C---- test file size
C
      CALL QQINQ(LUN,'DISKIO',REPLY,ILENGTH)
      ISTAT = LDUM*LBUF*LBUF
      WRITE (6,'(A,2(I8,A)//)')
     &' DISKIO should be',ISTAT,' bytes; is',ILENGTH,' bytes.'
      IF (ILENGTH.NE.ISTAT) CALL CCPERR(1, '*** FILE SIZE ERROR ***')
C
C---- Seed random Number Generator
C
      REPLY = ' '
      CALL UGTENV('SEED',REPLY)
      IF (REPLY.NE.' ') READ (REPLY,*) ISEED
C
C---- Get number of reads to perform
C
      CALL UGTENV('READS',REPLY)
      IF (REPLY.NE.' ') READ (REPLY,*) ILOOP
C
C---- Do random reads & writes on the file
C
      CALL UGTENV('NUMRECORD',ENVNAM)
      IF (ENVNAM.EQ.' ') THEN
        DO 40 I = 1,ILOOP
          NREC = NINT(100.*RANU(ISEED) + 1.)
          IF (NREC.LE.0 .OR. NREC.GT.LBUF)
     &    CALL CCPERR(1,'*** RECORD ERROR ***')
          CALL QSEEK(LUN,NREC,1,LBUF)
          IF (RANU(ISEED).LT..5) THEN
            CALL QREAD(LUN,BUFFER,LBUF,IER)
            WRITE (LUNOUT,FMT=6014) NREC,BUFFER(NREC),IER
            IF (BUFFER(NREC).NE.NREC)
     &      CALL CCPERR(1,'*** VERIFY ERROR ***')
          ELSE
            DO 70 J = 1,LBUF
70            BUFFER(J) = 0.
            BUFFER(NREC) = NREC
            CALL QWRITE(LUN,BUFFER,LBUF)
            WRITE (LUNOUT,FMT=6015) NREC,BUFFER(NREC)
          ENDIF
   40   CONTINUE
      ELSE
   50   CONTINUE
          CALL NOCRLF('Record to seek (-ve to write, Ctrl/Z to stop)> ')
          READ (LUNIN,6020,END=60) NREC
          I=IABS(NREC)
          IF (I.EQ.0 .OR. I.GT.LBUF) THEN
            WRITE (LUNOUT,*) '*** RECORD ERROR ***'
            GOTO 50
          ENDIF
          CALL QSEEK(LUN,I,1,LBUF)
          IF (NREC.GT.0) THEN
            CALL QREAD(LUN,BUFFER,LBUF,IER)
            WRITE (LUNOUT,FMT=6014) NREC,BUFFER(NREC),IER
            IF (BUFFER(NREC).NE.NREC)
     &      WRITE (LUNOUT,*) '*** VERIFY ERROR ***'
          ELSE
            DO 80 J = 1,LBUF
80            BUFFER(J) = 0.
            BUFFER(I) = I
            CALL QWRITE(LUN,BUFFER,LBUF)
            WRITE (LUNOUT,FMT=6015) I,BUFFER(I)
          ENDIF
        GOTO 50
      ENDIF
   60 CALL QCLOSE (LUN)
      CALL CUNLINK ('DISKIO')
C     Now check we can open and close a scratch file
      CALL QOPEN (LUN, 'foo.bar', 'SCRATCH')
      CALL QCLOSE (LUN)
C     and can we rewind a scratch file?  (make sure something's been
C     written to it first)
      I = 0
      CALL CCPDPN (LUN,'FOOEY','SCRATCH','F',0,I)
      WRITE (LUN,'(A)') 'foo'
      REWIND (LUN,ERR=170)
      READ (LUN,'(A)') FOO
      WRITE (LUNOUT, *) 'contents of temp file: ', FOO
      CALL CCPERR(0,'Normal Termination')
 170  CALL CCPERR (1,'Can''t rewind scratch file')
90    CALL CCPERR(1,'*** EOF ERROR ***')
C
C---- Format Statements
C
 6000 FORMAT (//A,/)
 6002 FORMAT (A7,3X,A12,A)
 6004 FORMAT (A,I3)
 6006 FORMAT (I2.2,'/',I2.2,'/',I4)
 6008 FORMAT (A5,I3)
 6014 FORMAT (' Seek Record:',I5,'  Read:  ',F8.2,'  Status: ',I4)
 6015 FORMAT (' Seek Record:',I5,'  Write: ',F8.2)
 6016 FORMAT (F8.2)
 6020 FORMAT (I10)
      END
C
C
      REAL FUNCTION RANU(K)
C==== UNIFORM PSEUDO-RANDOM NUMBER IN THE RANGE >= 0 AND < 1.
C
C     Set the seed K zero or negative to start or restart the sequence.
C     K must be a variable since a new value is returned each time.
C
      IMPLICIT           NONE
      INTEGER            M, IA, IC
      REAL               RM
      PARAMETER         (M=714025, IA=1366, IC=150889, RM=1./M)
      INTEGER            J, K, IY, IR(97)
      LOGICAL            FF
      SAVE               IY, IR, FF
      DATA               FF /.TRUE./
C
      IF (K.LE.0 .OR. FF) THEN
        FF = .FALSE.
        K = MOD(IC-K,M)
        DO 11 J = 1, 97
          K = MOD(IA*K + IC, M)
11        IR(J) = K
        K = MOD(IA*K + IC, M)
        IY = K
      ENDIF
C
      J = 1 + 97*IY/M
      IY = IR(J)
      RANU = IY*RM
      K = MOD(IA*K + IC, M)
      IR(J) = K
      END
