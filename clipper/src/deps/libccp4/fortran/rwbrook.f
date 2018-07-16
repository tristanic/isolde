C
C     rwbrook.f: Fortran interface to MMDB for handling coordinates
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
C**************************************************************************
C
C
C**************************************************************************
C
C
C Modifications
C
C     CCB 23/4/01
C     test to determine if H or R spacegroup, XYZOPEN2
C
      SUBROUTINE XYZINIT()
C     ====================
C
C_BEGIN_XYZINIT
C
C	This subroutine initialises the common block RBRKAA ready for reading
C and writing coordinate files. Also, the common blocks associated with 
C storing the header information are initialised. It should be called only 
C once from the top of the program.
C
C Parameters:
C
C       NONE
C
C COMMONS:
C
C         /RBRKAA/ FILESOPEN,LOGUNIT(MAXFILESOPEN),UNIT(MAXFILESOPEN),
C                  TYPE(MAXFILESOPEN)
C
C           FILESOPEN       no. of current coordinate files open.
C           LOGUNIT          logical name of file
C           UNIT            if the file is PDB then the unit is the physical
C                           channel opened. If CIF then is related to blocks.
C           TYPE            indicates whether PDB (1,-1) or CIF (2,-2). If
C                           negative then file is output file.
C
C_END_XYZINIT
C
C     .. Parameters ..
      INTEGER MAXFILESOPEN
      PARAMETER (MAXFILESOPEN=90)
C     ..
C     .. Variables in Common ..
      REAL CELL,CELLAS,RF,RO,RR,VOL
      INTEGER FILESOPEN,ITYP,NCODE,TYPE,UNIT
      CHARACTER LOGUNIT*80,BRKSPGRP*15
      LOGICAL IFCRYS,IFHDOUT,IFNEWCRYS,IFSCAL,MATRIX
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      EXTERNAL MMDB_F_INIT,SIMRWBROOK
C     ..
C     .. Common Blocks ..
      COMMON /RBRKAA/ FILESOPEN,LOGUNIT(MAXFILESOPEN),
     +                UNIT(MAXFILESOPEN),TYPE(MAXFILESOPEN)
      COMMON /RBRKXX/ IFCRYS,IFSCAL,ITYP,MATRIX,IFHDOUT,IFNEWCRYS
      COMMON /RBRKZZ/ CELL(6),RR(3,3,6),VOL,CELLAS(6)
      COMMON /ORTHOG/ RO(4,4),RF(4,4),NCODE
      COMMON /ORTHOGU/ ROU(4,4),RFU(4,4)
      COMMON /RBRKSPGRP/BRKSPGRP
C     ..
C     .. Save ..
      SAVE /RBRKAA/,/RBRKXX/,/RBRKZZ/,/ORTHOG/,/ORTHOGU/
C     ..

      FILESOPEN = 0

      DO 10 I=1,MAXFILESOPEN
        LOGUNIT(I) = ' '
        UNIT(I) = 0
        TYPE(I) = 0
   10 CONTINUE

      DO 20 I=1,6
        CELL(I) = 0.0
        CELLAS(I) = 0.0
   20 CONTINUE

      IFCRYS=.FALSE.
      IFSCAL=.FALSE.
      MATRIX=.FALSE.
      IFHDOUT=.FALSE.
      IFNEWCRYS=.FALSE.
      NCODE=0
      ITYP=0
C
C
      DO 30 I=1,3
        DO 40 J=I+1,4
          RO(I,J)=0.0
          RO(J,I)=0.0
          RF(I,J)=0.0
          RF(J,I)=0.0
          ROU(I,J)=0.0
          ROU(J,I)=0.0
          RFU(I,J)=0.0
          RFU(J,I)=0.0
40      CONTINUE
        RO(I,I)=1.0
        RF(I,I)=1.0
        ROU(I,I)=1.0
        RFU(I,I)=1.0
30    CONTINUE
C
C
      RO(4,4)=1.0
      RF(4,4)=1.0
      ROU(4,4)=1.0
      RFU(4,4)=1.0

      BRKSPGRP = ' '

C... mmdb fortran call to initialise
      CALL MMDB_F_INIT
C... mmdb mimic rwbrook
      CALL SIMRWBROOK(1)

      RETURN
      END
C
C
C
      SUBROUTINE RBINIT(IUNIT)
C     ========================
C
C_BEGIN_RBINIT
C
C      This routine is obsolete and should be removed.
C
C_END_RBINIT
C
C     .. Parameters ..
      INTEGER MAXFILESOPEN
      PARAMETER (MAXFILESOPEN=90)
C     ..
C     .. Scalar Arguments ..
      INTEGER IUNIT
C     ..
C     .. Local Scalars ..
      INTEGER I,J,IRET
C     ..
C     .. Scalars in Common ..
      INTEGER FILESOPEN,NCODE,ITYP
      LOGICAL IFHDOUT,IFNEWCRYS,IFCRYS,IFSCAL,MATRIX
C     ..
C     .. Arrays in Common ..
      REAL RF,RO
      INTEGER UNIT,TYPE
      CHARACTER LOGUNIT*80
C     ..
C     .. External Routines ..
      EXTERNAL MMDB_F_INIT,MMDB_F_REWD,SIMRWBROOK
C     ..
C     .. Common Blocks ..
      COMMON /RBRKAA/ FILESOPEN,LOGUNIT(MAXFILESOPEN),
     +                UNIT(MAXFILESOPEN),TYPE(MAXFILESOPEN)
      COMMON /RBRKXX/IFCRYS,IFSCAL,ITYP,MATRIX,IFHDOUT,IFNEWCRYS
      COMMON /ORTHOG/RO(4,4),RF(4,4),NCODE
      COMMON /ORTHOGU/ ROU(4,4),RFU(4,4)
C     ..
C     .. Save Statement ..
      SAVE /RBRKAA/,/RBRKXX/,/ORTHOG/
C     ..
C
C---- Make sure that MMDB_F_INIT is only called once. However,
C     this is not fool proof
C
      DO 10 I=1,10
        IF (IUNIT .EQ. UNIT(I)) GOTO 20
   10 CONTINUE

C...  mmdb, initialise and simulate rwbrook
      CALL MMDB_F_INIT 
      CALL SIMRWBROOK(1)

C...  mmdb, rewind channel
   20 CALL MMDB_F_REWD(IUNIT,IRET)
      IFCRYS=.FALSE.
      IFSCAL=.FALSE.
      MATRIX=.FALSE.
      NCODE=0
      ITYP=0
      IBRKFL=0
C
C
      DO 30 I=1,3
        DO 40 J=I+1,4
          RO(I,J)=0.0
          RO(J,I)=0.0
          RF(I,J)=0.0
          RF(J,I)=0.0
          ROU(I,J)=0.0
          ROU(J,I)=0.0
          RFU(I,J)=0.0
          RFU(J,I)=0.0
   40   CONTINUE
        RO(I,I)=1.0
        RF(I,I)=1.0
        ROU(I,I)=1.0
        RFU(I,I)=1.0
   30 CONTINUE
      RO(4,4)=1.0
      RF(4,4)=1.0
      ROU(4,4)=1.0
      RFU(4,4)=1.0


      RETURN
      END
C
C
C
      SUBROUTINE XYZOPEN(LOGNAM,RWSTAT,FILTYP,IUNIT,IFAIL)
C     =====================================================
C
C_BEGIN_XYZOPEN
C
C      Calls XYZOPEN2 which has an extra paramater for ignoring
C      Symmetry and Cryst cards
C
C Parameters:
C
C         LOGNAM (I)   CHARACTER*(*) but maximum of eight? The logical unit 
C                                    to which the file is assigned
C         RWSTAT (I)   CHARACTER*(*) if 'INPUT' then file is readonly
C                                    if 'OUTPUT' then file is an output file.
C         FILTYP (I)   CHARACTER*(*) if 'CIF' then the file type is treated as
C                                    CIF. If 'PDB' then the file type is 
C                                    treated as PDB. If blank then file type is
C                                    automatically determined for input files 
C                                    and for output file the file type will be
C                                    the same as the first file opened or 
C                                    defaulted to PDB.
C         IUNIT  (I/O) INTEGER       If zero then unit is decided else
C		  	             file opened on that unit
C                                    checked against previous data if 
C                                    applicable. NOT used with output files.
C         IFAIL  (I/O) INTEGER       On input    = 0 stop on failure 
C                                                = 1 continue on failure
C
C                                    On output   unchanged if OK
C                                                = -1  if error
C
C
C_END_XYZOPEN
C
C   
      implicit none
C     ..
C     .. Arguments ..
      INTEGER IUNIT, IFAIL,ICRYST
      CHARACTER*(*) FILTYP,LOGNAM,RWSTAT
      ICRYST = 0

      CALL XYZOPEN2(LOGNAM,RWSTAT,FILTYP,IUNIT,IFAIL,ICRYST)

      RETURN
      END
C
C
C
      SUBROUTINE XYZOPEN2(LOGNAM,RWSTAT,FILTYP,IUNIT,IFAIL,ICRYST)
C     =====================================================
C
C_BEGIN_XYZOPEN2
C
C      Opens a coordinate file for input or output. The channel number can
C be determined automatically or set on input. The header info.
C is also read: cell, orthog. matrix and symmetry.
C This is a version of XYZOPEN with an extra argument to flag whether or
C not the CRYST and SCALE cards are required.
C
C Parameters:
C
C         LOGNAM (I)   CHARACTER*(*) but maximum of eight? The logical unit 
C                                    to which the file is assigned
C         RWSTAT (I)   CHARACTER*(*) if 'INPUT' then file is readonly
C                                    if 'OUTPUT' then file is an output file.
C         FILTYP (I)   CHARACTER*(*) if 'CIF' then the file type is treated as
C                                    CIF. If 'PDB' then the file type is 
C                                    treated as PDB. If blank then file type is
C                                    automatically determined for input files 
C                                    and for output file the file type will be
C                                    the same as the first file opened or 
C                                    defaulted to PDB.
C         IUNIT  (I/O) INTEGER       If zero then unit is decided else
C		  	             file opened on that unit
C                                    checked against previous data if 
C                                    applicable. NOT used with output files.
C         IFAIL  (I/O) INTEGER       On input    = 0 stop on failure 
C                                                = 1 continue on failure
C
C                                    On output   unchanged if OK
C                                                = -1  if error
C
C         ICRYST (I/O) INTEGER       If zero, then check for and use CRYST
C                                    and SCALE cards in input PDB.
C                                    If one, ignore these cards even if
C                                    present.
C
C
C_END_XYZOPEN2
C
C   
      implicit none
C     .. Parameters ..
      INTEGER MAXFILESOPEN,MAXSYM
      PARAMETER (MAXFILESOPEN=90,MAXSYM=96)
      INTEGER RWBERR_Ok,RWBERR_NoMatrices
      PARAMETER (RWBERR_Ok=0,RWBERR_NoMatrices=-16)
C     ..
C     .. Arguments ..
      INTEGER IFAIL,IUNIT,ICRYST,II,III,JJ,K,ISYM
      REAL AM,BM,RCHK1,RCHK2,FAC
      CHARACTER*(*) FILTYP,LOGNAM,RWSTAT
C     ..
C     .. Variables in Common ..
      REAL CELL,CELLAS,RF,RO,RR,VOL,ROU,RFU
      INTEGER FILESOPEN,ITYP,NCODE,TYPE,UNIT
      CHARACTER*80 LOGUNIT,BRKSPGRP*15
      LOGICAL IFCRYS,IFHDOUT,IFNEWCRYS,IFSCAL,MATRIX
c  Check symmetry stuff
C
      integer  nsymchk,lspgrp,nsymppdbs,nsympdbs,ist
      real rsymchk(4,4,maxsym), rsympdbs(4,4,maxsym)
      character NAMSPG_CIFS*20,nampg*10
C     ..
C     .. Local Scalars ..
      REAL ERROR,VOLCHK
      INTEGER IRET
      INTEGER I,IORTH,IFILTYP,J,LL,ILEN,KLEN
      CHARACTER ERRLIN*600,FILNAM*255
      CHARACTER LFILTYP*3,LRWSTAT*5, SPGCHK*30
      CHARACTER*40 ORTH(6)
C     ..
C     .. Local Arrays ..
      REAL P(4,4)
      CHARACTER IEC(3)*2
C     ..
C     .. External Functions ..
      INTEGER LENSTR,CCPNUN
      LOGICAL CCPEXS
      EXTERNAL CCPEXS,LENSTR,CCPNUN
C     ..
C     .. External Routines ..
      EXTERNAL CCPDPN,CCPERR,CCPUPC,RBFROR,RBRINV,UGTENV,
     *          MMDB_F_OPENL,MMDB_F_RBSPGRP,RBERRSTOP,
     *          MMDB_F_RBCELL,MMDB_F_WBSPGRP,MMDB_F_RBORF

C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common Blocks ..
      COMMON /RBRKAA/ FILESOPEN,LOGUNIT(MAXFILESOPEN),
     +                UNIT(MAXFILESOPEN),TYPE(MAXFILESOPEN)
      COMMON /RBRKXX/IFCRYS,IFSCAL,ITYP,MATRIX,IFHDOUT,IFNEWCRYS
      COMMON /RBRKZZ/CELL(6),RR(3,3,6),VOL,CELLAS(6)
      COMMON /ORTHOG/RO(4,4),RF(4,4),NCODE
      COMMON /ORTHOGU/ ROU(4,4),RFU(4,4)
      COMMON /RBRKSPGRP/BRKSPGRP
C     ..
C     .. Save ..
      SAVE /RBRKAA/,/RBRKXX/,/RBRKZZ/,/ORTHOG/
C     ..
C     .. Data Statement ..
      DATA IEC /'E1','E2','E3'/
      DATA ORTH/'A // XO, C* // ZO (STANDARD PDB)',
     *          'B // XO, A* // ZO',
     *          'C // XO, B* // ZO',
     *          'HEX A+B // XO, C* // ZO',
     *          'A* // XO, C // ZO (ROLLETT)',
     *          'A // XO, B* // YO'/
C     ..
      I = 0
      II = 0
      IFILTYP = 1
      LL = 0
      LRWSTAT = RWSTAT
      LFILTYP = FILTYP
      CALL CCPUPC(LFILTYP)
      CALL CCPUPC(LRWSTAT)
C
C---- If too many files opened
C
      IF (FILESOPEN .EQ. MAXFILESOPEN) THEN
        CALL CCPERR(1,' *** ERROR: too many coordinate files open. ***')
      ENDIF
C...
C...  MMDB open channel (need to test for failures)
      FILNAM = ' '
      CALL UGTENV(LOGNAM,FILNAM)
      IF (LRWSTAT(1:5) .EQ. 'INPUT' .AND. FILNAM.EQ.' ') THEN
        CALL MMDB_F_OPEN(LOGNAM,LRWSTAT,LFILTYP,IUNIT,IRET)
      ELSEIF (LRWSTAT(1:5) .EQ. 'INPUT') THEN
        CALL MMDB_F_OPENL(LOGNAM,LRWSTAT,LFILTYP,IUNIT,IRET)
      ELSE
        IF (LFILTYP(1:1).EQ.' ' .AND. FILESOPEN.GE.1) THEN
          IF (ABS(TYPE(1)).EQ.1) THEN
            IFILTYP = 1
            LFILTYP = 'PDB'
          ELSEIF (ABS(TYPE(1)).EQ.2) THEN
            IFILTYP = 2
            LFILTYP = 'CIF'
          ENDIF
        ENDIF
        CALL MMDB_F_OPENL(LOGNAM,LRWSTAT,LFILTYP,IUNIT,IRET)
      ENDIF
      IF (IRET.NE.RWBERR_Ok) THEN
        IF (IRET.NE.RWBERR_NoMatrices) THEN
C...  MMDB error information
          CALL RBERRSTOP(1,IRET,IUNIT,1)
          ERRLIN = ' XYZOPEN: Error opening logical name '//LOGNAM
          IF (IFAIL .EQ. 0) THEN
            CALL CCPERR(1,ERRLIN(1:LENSTR(ERRLIN)))
          ELSE
            CALL CCPERR(2,ERRLIN(1:LENSTR(ERRLIN)))
            IFAIL = -1
            RETURN
          ENDIF
        ENDIF
      ENDIF
C
C==== If the file is an INPUT file
C
      IF (LRWSTAT(1:5) .EQ. 'INPUT') THEN
C
C---- Determine whether CIF or PDB
C
        IF (LFILTYP(1:1) .EQ. ' ') THEN
          FILESOPEN = FILESOPEN + 1
          LOGUNIT(FILESOPEN) = LOGNAM
          UNIT(FILESOPEN) = IUNIT
          TYPE(FILESOPEN) = 1
        ENDIF  
C
C---- If known as CIF
C
       IF (LFILTYP(1:3).EQ.'CIF' .OR. IFILTYP.EQ.2) THEN
          FILESOPEN = FILESOPEN + 1
          LOGUNIT(FILESOPEN) = LOGNAM
          UNIT(FILESOPEN) = IUNIT
          TYPE(FILESOPEN) = 2
        ENDIF
C
C---- If known as a PDB file
C
        IF (LFILTYP(1:3).EQ.'PDB') THEN
          FILESOPEN = FILESOPEN + 1
          LOGUNIT(FILESOPEN) = LOGNAM
          UNIT(FILESOPEN) = IUNIT
          TYPE(FILESOPEN) = 1
        ENDIF
C
C---- Cell card found - calculate standard orthogonalising matrix
C     Check if you already have a cell which is inconsistent with 
C     this one
C
        IF (ICRYST.EQ.0) THEN
            ITYP=1
            BRKSPGRP = ' '
C...  MMDB get spacegroup and cell (cache)
            CALL MMDB_F_RBSPGRP(IUNIT,SPGCHK,IRET)
            ILEN = LENSTR(SPGCHK)
            IF (ILEN.LE.1) THEN 
               CALL CCPERR(2,' No Space group given on PDB CRYST1 line')
            ELSE
               BRKSPGRP = SPGCHK
            ENDIF
C...
C...  MMDB get CELL and VOL for cache
            CALL MMDB_F_RBCELL(IUNIT,CELL,VOL,IRET)
            IF (IRET.EQ.0) THEN
              IFCRYS=.TRUE.
            ELSE
              IFCRYS=.FALSE.
            ENDIF
C...  MMDB get the orthogonalisation and fractional matrices
            RO(1,1) = 0.0
            CALL MMDB_F_RBORF(IUNIT,RO,RF,NCODE,IRET)
            IF (NCODE.LT.0) THEN
               ERRLIN = 
     +        'XYZOPEN2: Orthogonalisation code not determined.  '//
     +        'Possible disagreement between CRYST1 and SCALEx cards.'
               CALL CCPERR(2,ERRLIN(1:LENSTR(ERRLIN)))
            ENDIF
C IFSCAL indicates SCALEx cards found, MATRIX indicates RO,RF in /ORTHOG/ set up
            IF (IRET.EQ.0) THEN
              IFSCAL=.TRUE.
              MATRIX=.TRUE.
            ELSE
              IFSCAL=.FALSE.
              MATRIX=.FALSE.
            ENDIF
C
C If BRKSPGRP contains "/" it is probably a Patterson group and may
C occupy the full 15 characters. Else it may be from the PDB and
C may contain Z value which must be removed.

             IF (INDEX(BRKSPGRP,'/').EQ.0) BRKSPGRP(12:15) = ' '

C  Read symmetry operators.
             IF (BRKSPGRP.NE.' ') THEN
C  Consistency check for R and H space groups.
C      (using angles only)
C---- H name associated with a=b; Gamma = 120
                IF( BRKSPGRP(1:1) .EQ. 'R' ) THEN
                   IF(abs(cell(4)-90.0).LT.0.001 .AND.
     +                  abs(cell(5)-90.0).LT.0.001 .AND.
     +                  abs(cell(6)-120.0).LT.0.001) THEN
                      BRKSPGRP(1:1)='H'
                      CALL CCPERR(2,
     +                     ' Changing "rhombohedral" to "hexagonal"')
                      CALL MMDB_F_WBSPGRP(IUNIT,BRKSPGRP,IRET)
                   END IF
C---- R name associated with a=b=c; Alpha=Beta=Gamma 
                ELSE IF( BRKSPGRP(1:1) .EQ. 'H' ) THEN
                   IF(abs(cell(4)-cell(5)).LT.0.001.AND.
     +                  abs(cell(5)-cell(6)).LT.0.001.AND.
     +                  abs(cell(6)-cell(4)).LT.0.001) THEN
                      BRKSPGRP(1:1)='R'
                      CALL CCPERR(2,
     +                     ' Changing "hexagonal" to "rhombohedral"')
                      CALL MMDB_F_WBSPGRP(IUNIT,BRKSPGRP,IRET)
                   END IF
                END IF
C
             END IF
C...   set up standard orth matrices (cache)
            IF (IFCRYS) CALL RBFROR
          ENDIF
C
C==== If the file is an OUTPUT file
C
      ELSE

        IF (LFILTYP(1:3) .EQ. 'CIF') IFILTYP = 2 
        IF (LFILTYP(1:3) .EQ. 'PDB') IFILTYP = 1
C
C---- Open output PDB file
C
        IF (IFILTYP .EQ. 1) THEN
          FILESOPEN = FILESOPEN + 1
          LOGUNIT(FILESOPEN) = LOGNAM
          UNIT(FILESOPEN) = IUNIT
          TYPE(FILESOPEN) = -1
        ENDIF

        IF (IFILTYP .EQ. 2) THEN
          FILESOPEN = FILESOPEN + 1
          LOGUNIT(FILESOPEN) = LOGNAM
          UNIT(FILESOPEN) = IUNIT
          TYPE(FILESOPEN) = -2
        ENDIF
      ENDIF

C     Generate ROU and RFU for AnisoU stuff
       IF (MATRIX) THEN
         RFU(4,4) = 1.0
         DO I = 1,3
         FAC= SQRT(RF(I,1)*RF(I,1) + RF(I,2)*RF(I,2) +
     +              RF(I,3)*RF(I,3))
         RFU(I,1) = RF(I,1)/FAC
         RFU(I,2) = RF(I,2)/FAC
         RFU(I,3) = RF(I,3)/FAC
         RFU(I,4) = 0.0
         RFU(4,I) = 0.0
         END DO 
         CALL RBRINV(RFU,ROU)
       END IF 
C     If reading in: check SCAL and CRYST1 cards
      IF (IFILTYP.EQ.1 .AND. LRWSTAT(1:5).EQ.'INPUT' .AND. ICRYST.EQ.0)
     +              THEN
        IF (.NOT.IFCRYS) THEN
          WRITE(ERRLIN,FMT='(A,A)') ' NO CRYST CARDS READ FROM ',LOGNAM
          CALL CCPERR (2,ERRLIN(1:LENSTR(ERRLIN)))
        END IF
        IF (.NOT.IFSCAL) THEN
          WRITE(ERRLIN,FMT='(A,A)') ' NO SCALE CARDS READ FROM ',LOGNAM
          CALL CCPERR (2,ERRLIN(1:LENSTR(ERRLIN)))
        END IF
      END IF
      RETURN
      END
C
C
C
      SUBROUTINE XYZCLOSE(IUNIT)
C     ==========================
C
C_BEGIN_XYZCLOSE
C
C	This subroutine closes a coordinate file. 
C
C Parameters:
C
C         IUNIT  (I)   INTEGER    Unit number for file
C
C_END_XYZCLOSE
C
C     .. Parameters ..
      INTEGER MAXFILESOPEN
      PARAMETER (MAXFILESOPEN=90)
C     ..
C     .. Arguments ..
      INTEGER IUNIT
C     ..
C     .. Variables in Common ..
      INTEGER FILESOPEN, UNIT, TYPE
      CHARACTER*80 LOGUNIT
C     ..
C     .. Local Scalars ..
      INTEGER I,II,IRET
C     ..
      EXTERNAL MMDB_F_CLOSE
C     ..
C     .. Common Blocks ..
      COMMON /RBRKAA/ FILESOPEN,LOGUNIT(MAXFILESOPEN),
     +                UNIT(MAXFILESOPEN),TYPE(MAXFILESOPEN)
C     ..
C     .. Save ..
      SAVE /RBRKAA/
C     ..
C
      II = 0

      DO 10 I=1,FILESOPEN
        IF (UNIT(I) .EQ. IUNIT) THEN
          II = I
          GOTO  20
        ENDIF
   10 CONTINUE

   20 IF (II .NE. 0) THEN
        IF (FILESOPEN.NE.1 .AND. II.NE.FILESOPEN) THEN
          LOGUNIT(II) = LOGUNIT(FILESOPEN)
          UNIT(II) = UNIT(FILESOPEN)
          TYPE(II) = TYPE(FILESOPEN)
        ENDIF
        FILESOPEN = FILESOPEN - 1
      ENDIF

C... MMDB, close the channel (some info in cache)
      CALL MMDB_F_CLOSE(IUNIT,IRET)

      RETURN
      END
C
C
C
      SUBROUTINE XYZADVANCE(IUNIT,IOUT,ITER,*,*)
C     ==========================================
C
C_BEGIN_XYZADVANCE
C
C When IUNIT represents an input coordinate file, this subroutine reads 
C recognised data lines into a buffer BROOK, from which XYZATOM and 
C XYZCOORD can extract useful information. Optionally, if the card is 
C unrecognised (eg REMARK) then the line can be echoed to an output file.
C
C When IUNIT represents an output coordinate file, this subroutine writes 
C out the contents of the buffer WBROOK or WBROOK1. This buffer is filled
C from an input file, or by calls to XYZATOM and XYZCOORD.
C
C Parameters:
C
C      IUNIT  (I) Channel number of the coordinate file
C
C      These arguments are not relevant for output files.
C        IOUT (I) Logical unit number to which non-atom/hetatm/anisou records 
C                 are to be written (may be blank if reading only)
C        ITER (I) FLAG =1, return if 'ter' card found (via return 1)
C                      =0, do not return when 'ter' card found
C      RETURN 1   Return on TER card if ITER=1
C      RETURN 2   Return on end of file.
C
C_END_XYZADVANCE
C
C     .. Paramters ..
      INTEGER MAXFILESOPEN
      PARAMETER (MAXFILESOPEN=90)
      INTEGER RWBERR_Ok     
      PARAMETER (RWBERR_Ok=0)
C     ..
C     .. Arguments ..
      INTEGER IOUT,ITER,IUNIT
C     ..
C     .. Variables in Common ..
      REAL CELL,CELLAS,RF,RO,RR,VOL
      INTEGER FILESOPEN,NCODE,TYPE,UNIT
      CHARACTER LOGUNIT*80,BRKSPGRP*15
      LOGICAL IFCRYS,IFHDOUT,IFNEWCRYS,IFSCAL,MATRIX
C     ..
C     .. Local Variables ..
      INTEGER I,II,IRET,IRE,ITE
      CHARACTER*80 ERRLIN
      CHARACTER*6 ITYPE(7)
C     ..
C     .. External Routines ..
      EXTERNAL CCPERR,MMDB_F_COPY,MMDB_F_ADVANCE
C     ..
C     .. Common Blocks ..
      COMMON /RBRKAA/ FILESOPEN,LOGUNIT(MAXFILESOPEN),
     +                UNIT(MAXFILESOPEN),TYPE(MAXFILESOPEN)
      COMMON /RBRKXX/ IFCRYS,IFSCAL,ITYP,MATRIX,IFHDOUT,IFNEWCRYS
      COMMON /RBRKZZ/ CELL(6),RR(3,3,6),VOL,CELLAS(6)
      COMMON /ORTHOG/ RO(4,4),RF(4,4),NCODE
      COMMON /ORTHOGU/ ROU(4,4),RFU(4,4)
      COMMON /RBRKSPGRP/BRKSPGRP
C     ..
C     .. Save Statement ..
      SAVE /RBRKAA/,/RBRKXX/,/RBRKZZ/,/ORTHOG/

      II = 0
      DO 10 I=1,FILESOPEN
        IF (IUNIT .EQ. UNIT(I)) THEN
          II = I
          GOTO 15
        ENDIF
   10 CONTINUE

      ERRLIN = ' ERROR: in XYZADVANCE file has not been opened'
      CALL CCPERR(1,ERRLIN)
C
C...  opened read
 15   IF (TYPE(II).GT.0) THEN
C...    copying across
        IF (IOUT.NE.0) THEN
          IF(.NOT.IFHDOUT)THEN
            IF(IFCRYS .AND. .NOT.IFNEWCRYS) THEN
C..   MMDB, copy header information
              CALL MMDB_F_COPY(IOUT,IUNIT,2,IRE)
              IFHDOUT=.TRUE.
            ELSE 
C...  MMDB, copy header less CRYST1 cards
              CALL MMDB_F_COPY(IOUT,IUNIT,3,IRE)
            ENDIF
          ENDIF 
C...  MMDB, advance pointer allowing for TER cards
          CALL MMDB_F_ADVANCE(IUNIT,IOUT,ITER,IRET)
        ELSE
C...   MMDB, standard advance pointer
          CALL MMDB_F_ADVANCE(IUNIT,IOUT,ITER,IRET)
        ENDIF
      ELSE
        IF(.NOT.IFHDOUT .AND. IFCRYS) THEN
          CALL MMDB_F_WBSPGRP(IOUT,BRKSPGRP,IRE)
          CALL MMDB_F_WBCELL(IOUT,CELL,NCODE,IRE)
          IFHDOUT=.TRUE.
        ENDIF
        CALL MMDB_F_ADVANCE(IUNIT,IOUT,ITER,IRET)
      ENDIF

      IF (IRET.EQ.1) RETURN 1
      IF (IRET.EQ.2) RETURN 2
C
      END
C
C
C
      SUBROUTINE XYZATOM(IUNIT,ISER,ATNAM,RESNAM,CHNNAM,IRESN,
     *RESNO,INSCOD,ALTCOD,SEGID,IZ,ID)
C     ========================================================
C
C_BEGIN_XYZATOM
C
C	This subroutine reads or writes the atom name, residue name, 
C chain name etc. into the buffer. XYZADVANCE actually advances a line 
C or atom. The character strings have undefined length in order to make 
C change easier. However, these data items will be strictly defined in 
C the working format.
C
C Parameters:
C
C       IUNIT  (I)  Logical unit of the input coordinate file
C        ISER (I/O) Atom serial number
C       ATNAM (I/O) Atom name        (left justified)
C      RESNAM (I/O) Residue name     
C      CHNNAM (I/O) Chain name       
C       IRESN (I/O) Residue number as an integer
C       RESNO  (O)  Residue number as character, NOT used for output file
C      INSCOD (I/O) The insertion code
C      ALTCOD (I/O) The alternate conformation code.
C       SEGID (I/O) Segid is here to complete PDB standard.
C          IZ  (O)  Atomic number (returned as 7 from ambiguous atoms),
C                   NOT used for output file
C          ID (I/O) Atomic ID related to atomic number (element symbol
C                   right justified), plus the ionic state +2, +3 etc..
C
C_END_XYZATOM
C     
C     .. Paramters ..
      INTEGER MAXFILESOPEN
      PARAMETER (MAXFILESOPEN=90)
      INTEGER RWBERR_Ok,RWBWAR_MASK 
      PARAMETER (RWBERR_Ok=0,RWBWAR_MASK=16384)
      INTEGER RWBWAR_WrongSerial
      PARAMETER (RWBWAR_WrongSerial=16448)
C     ..
C     .. Arguments ..
      INTEGER IRESN,ISER,IUNIT,IZ
      CHARACTER*(*) RESNO,ATNAM,RESNAM,CHNNAM
      CHARACTER*(*) ID,INSCOD,ALTCOD,SEGID
      CHARACTER*6 ITYPE(7)
C     ..
C     .. Variables in Common ..
      INTEGER FILESOPEN,ITYP,UNIT,TYPE
      CHARACTER LOGUNIT*80
      LOGICAL IFCRYS,IFHDOUT,IFNEWCRYS,IFSCAL,MATRIX
C     ..
C     .. Local Scalars ..
      REAL U(6),OCC,X,Y,Z
      INTEGER I,II,IRET
      CHARACTER*100 ERRLIN
C     ..
C     .. External Routines/Functions ..
      EXTERNAL CCPERR,MMDB_F_ATOM,RBERRSTOP
C     ..
C     .. Common Blocks ..
      COMMON /RBRKAA/ FILESOPEN,LOGUNIT(MAXFILESOPEN),
     +                UNIT(MAXFILESOPEN),TYPE(MAXFILESOPEN)
      COMMON /RBRKXX/ IFCRYS,IFSCAL,ITYP,MATRIX,IFHDOUT,IFNEWCRYS
C
C     .. Save ..
      SAVE /RBRKAA/,/RBRKXX/
C     ..
c
      II = 0

      DO 10 I=1,FILESOPEN
        IF (IUNIT .EQ. UNIT(I)) THEN
          II = I
          GOTO 20
        ENDIF
   10 CONTINUE

      ERRLIN = ' ERROR: in XYZATOM file has not been opened'
      CALL CCPERR(1,ERRLIN)
C
   20 CALL MMDB_F_ATOM(IUNIT,ISER,ATNAM,RESNAM,CHNNAM,IRESN,
     +    RESNO,INSCOD,ALTCOD,SEGID,IZ,ID,IRET)
      IF (IRET.NE.RWBERR_Ok.AND.IRET.NE.RWBWAR_WrongSerial) 
     +  THEN
C...  MMDB error information
        CALL RBERRSTOP(1,IRET,IUNIT,0)
        IF (IAND(IRET,RWBWAR_MASK).EQ.0) THEN
          ERRLIN = ' ERROR: XYZATOM'
          CALL CCPERR(1,ERRLIN) 
        ENDIF
      ENDIF

      RETURN
      END
C
C
C
      SUBROUTINE XYZCOORD(IUNIT,XFLAG,BFLAG,X,Y,Z,OCC,BISO,U)
C     =======================================================
C
C_BEGIN_XYZCOORD
C
C	This subroutine reads or writes x, y, z, occupancy and b from/to 
C the internal buffer. The buffer is updated from the file by 
C XYZADVANCE. The coordinates can be input/output (to the subroutine) 
C as orthogonal or fractional. 
C
C  PDB files contain anisotropic temperature factors as orthogonal Us.
C The anisotropic temperature factors can be input/output to this routine 
C  as orthogonal or as crystallographic Us. 
C  
C  Shelx defines Uf to calculate temperature factor as:
C T(aniso_Uf) = exp (-2PI**2 ( (h*ast)**2 Uf_11 + (k*bst)**2 Uf_22 + ... 
C                            + 2hk*ast*bst*Uf_12 +..)
C
C   Note:   Uo_ji == Uo_ij and  Uf_ji == Uf_ij.
C
C  [Uo_ij] listed on ANISOU card satisfy  the relationship:
C  [Uo_ij] =   [RFu]-1 [Uf_ij] {[RFu]-1}T   C
C        where [Rfu] is the normalised [Rf] matrix read from the SCALEi cards.
C        see code.   [ROu] ==  [RFu]-1
C
C T(aniso_Uo) = U(11)*H**2 + U(22)*K**2 + 2*U(12)*H*K + ...
C where H,K,L are orthogonal reciprocal lattice indecies. ( EJD: I think????)
C
C Biso     = 8*PI**2 (Uo_11 + Uo_22 + Uo_33) / 3.0
C
C   [Uf(symm_j)] = [Symm_j] [Uf] [Symm_j]T
C
C Parameters:
C
C       IUNIT  (I)  Channel number of the input coordinate file
C       XFLAG  (I)  For input file
C                     ='F' will get fractional coords. 
C                     ='O' will get orthogonal coords.
C                   For output file
C                     ='F' passed coordinates are fractional
C                     ='O' passed coordinates are orthogonal
C       BFLAG  (I)  For input file
C                     ='F' will get fractional us
C                     ='O' will get orthogonal Us.
C                   For output file
C                     ='F' have fractional us
C                     ='O' have othogonal Us
C           X (I/O) Coordinates (orthogonal angstrom coordinates as
C           Y (I/O)     "        stored)
C           Z (I/O)     "
C         OCC (I/O) Occupancy
C        BISO  (O)  Isotropic temperature factor, NOT used for output file.
C        U(6) (I/O) Orthogonal Anisotropic temperature factor, unless only U(1) defined.
C
C_END_XYZCOORD
C     
C     .. Paramters ..
      INTEGER MAXFILESOPEN
      PARAMETER (MAXFILESOPEN=90)
      INTEGER RWBERR_Ok,RWBWAR_MASK
      PARAMETER (RWBERR_Ok=0,RWBWAR_MASK=16384)
C     ..
C     .. Arguments ..
      REAL U(6),BISO,X,Y,Z
      INTEGER IUNIT
      CHARACTER*1 BFLAG,XFLAG
C     ..
C     .. Variables in Common ..
      INTEGER FILESOPEN,ITYP,UNIT,TYPE
      CHARACTER LOGUNIT*80
      LOGICAL IFCRYS,IFHDOUT,IFNEWCRYS,IFSCAL,MATRIX
C     ..
C     .. Local Scalars ..
      REAL EIGHTPI2,XX,YY,ZZ
      INTEGER I,II,IRET
      INTEGER IRESN,ISER,IZ
      CHARACTER*100 ERRLIN
      CHARACTER ATNAM*4,RESNAM*4,RESNO*4,ID*4,CHNNAM*1,SEGID*4
C     ..
C     .. External Routines/Functions ..
      EXTERNAL CCPERR,MMDB_F_COORD,RBERRSTOP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,NINT
C     ..
C     .. Common Blocks ..
      COMMON /RBRKAA/ FILESOPEN,LOGUNIT(MAXFILESOPEN),
     +                UNIT(MAXFILESOPEN),TYPE(MAXFILESOPEN)
      COMMON /RBRKXX/ IFCRYS,IFSCAL,ITYP,MATRIX,IFHDOUT,IFNEWCRYS
C     ..
C     .. Save ..
      SAVE /RBRKAA/,/RBRKXX/
C     ..
C     .. Data Statement ..
      DATA EIGHTPI2 /78.956835/

      II = 0
      DO 10 I=1,FILESOPEN
        IF (IUNIT .EQ. UNIT(I)) THEN
          II = I
          GOTO 20
        ENDIF
   10 CONTINUE

      ERRLIN = ' ERROR: in XYZCOORD has not been opened'
      CALL CCPERR(1,ERRLIN)

  20  CALL MMDB_F_COORD(IUNIT,XFLAG,BFLAG,X,Y,Z,OCC,BISO,U,IRET)
      IF (IRET.NE.RWBERR_Ok) THEN
C...  MMDB error information
        CALL RBERRSTOP(1,IRET,IUNIT,0)
        IF (IOR(IRET,RWBWAR_MASK).EQ.0) THEN
          ERRLIN = ' ERROR: XYZATOM'
          CALL CCPERR(1,ERRLIN)
        ENDIF
      ENDIF

      RETURN
      END
C
C
C
      SUBROUTINE XYZREWD(IUNIT)
C     =========================
C
C_BEGIN_XYZREWD
C
C	This routine is resets pointer to the begining of the file ie.
C rewind the file.
C
C Parameters:
C
C      IUNIT  (I) INTEGER  Channel number for file.
C
C_END_XYZREWD
C
C     .. Parameters ..
      INTEGER MAXFILESOPEN
      PARAMETER (MAXFILESOPEN=90)
      INTEGER RWBWAR_RewOutput
      PARAMETER (RWBWAR_RewOutput=34)
C     ..
C     .. Arguments ..
      INTEGER IUNIT
C     ..
C     .. Variables in Common ..
      INTEGER FILESOPEN,TYPE,UNIT
      CHARACTER*80 LOGUNIT
C     ..
C     .. Local Scalars ..
      INTEGER I,II
      CHARACTER*100 ERRLIN
C     ..
C     .. External Functions/Routines ..
      EXTERNAL CCPERR,MMDB_F_REWD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common Blocks ..
      COMMON /RBRKAA/ FILESOPEN,LOGUNIT(MAXFILESOPEN),
     +                UNIT(MAXFILESOPEN),TYPE(MAXFILESOPEN)
C     ..
C     .. Save ..
      SAVE /RBRKAA/
C     ..
      II = 0
      DO 10 I=1,FILESOPEN
        IF (IUNIT .EQ. UNIT(I)) THEN
          II = I
          GOTO 20
        ENDIF
   10 CONTINUE

      ERRLIN = ' ERROR: in XYZREWD file has not been opened'
      CALL CCPERR(1,ERRLIN)

   20 CALL MMDB_F_REWD(IUNIT,IRET)
      IF (IRET.EQ. RWBWAR_RewOutput) THEN
        CALL CCPERR(2,
     +    ' WARNING: you are rewinding an output file!!')
      ENDIF

      RETURN
      END
C
C
C
      SUBROUTINE XYZBKSP(IUNIT)
C     =========================
C
C_BEGIN_XYZBKSP
C
C	This routine is the opposite to XYZADVANCE in that it retreats 
C one atom ie. backspacing.
C
C Parameters:
C
C      IUNIT  (I) INTEGER  Channel number for file.
C
C_END_XYZBKSP
C
C     .. Parameters ..
      INTEGER MAXFILESOPEN
      PARAMETER (MAXFILESOPEN=90)
      INTEGER RWBWAR_RewOutput 
      PARAMETER (RWBWAR_RewOutput=34)
C     ..
C     .. Arguments ..
      INTEGER IUNIT
C     ..
C     .. Variables in Common ..
      INTEGER FILESOPEN,TYPE,UNIT
      CHARACTER*80 LOGUNIT
C     ..
C     .. Local Scalars ..
      INTEGER I,II,IRET
      CHARACTER*100 ERRLIN
C     ..
C     .. External Functions/Routines ..
      EXTERNAL CCPERR,MMDB_F_BKSP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Common Blocks ..
      COMMON /RBRKAA/ FILESOPEN,LOGUNIT(MAXFILESOPEN),
     +                UNIT(MAXFILESOPEN),TYPE(MAXFILESOPEN)
C     ..
C     .. Save ..
      SAVE /RBRKAA/
C     ..
      II = 0
      DO 10 I=1,FILESOPEN
        IF (IUNIT .EQ. UNIT(I)) THEN
          II = I
          GOTO 20
        ENDIF
   10 CONTINUE

      ERRLIN = ' ERROR: in XYZBKSP file has not been opened'
      CALL CCPERR(1,ERRLIN)

   20 CALL MMDB_F_BKSP(IUNIT,IRET) 
      IF (IRET .EQ. RWBWAR_RewOutput) THEN
        CALL CCPERR(2,
     +    ' WARNING: you are backspacing an output file!!')
      ENDIF

      RETURN
      END
C
C
C
      SUBROUTINE RBROOK(IUNIT,ISER,ATNAM,RESNAM,CHNNAM,IRESN,RESNO,IS,
     *X,Y,Z,OCC,B,IZ,IOUT,MSG1,MSG2,ITER,*,*)
C     ================================================================
C
C_BEGIN_RBROOK
C
C      This subroutine is obsolete and should be removed. May be 
C PROBLEM in that routine returns orthogonal coordinates and not fractional
C ones.
C
C_END_RBROOK
C
C     .. Scalar Arguments ..
      REAL X,Y,Z,OCC,B
      INTEGER IUNIT,ISER,IRESN,IS,IOUT,IZ,MSG1,MSG2,ITER
      CHARACTER*(*) RESNO
      CHARACTER ATNAM*4,RESNAM*4,CHNNAM*1
C     ..
C     .. Local Scalars ..
      CHARACTER*4 ID,SEGID
      CHARACTER*1 INSCOD,ALTCOD
C     ..
C     .. Local Arrays ..
      REAL U(6)
C     ..
C     .. External Routines ..
      EXTERNAL XYZADVANCE,XYZATOM,XYZCOORD
C     ..
      IS = 0

   10 CALL XYZADVANCE(IUNIT,IOUT,ITER,*1000,*2000)
      CALL XYZCOORD(IUNIT,'O','U',X,Y,Z,OCC,B,U)
C
C---- Skip ANISOU card in PDB. Test on X, Y and Z are not strictly necessary
C     as routines can only read PDB currently.
C
      IF (U(2).NE.0.0 .AND. U(3).NE.0.0) THEN
        IF (X.EQ.0.0 .AND. Y.EQ.0.0 .AND. Z.EQ.0.0) GOTO 10
      ENDIF
      CALL XYZATOM(IUNIT,ISER,ATNAM,RESNAM,CHNNAM,IRESN,RESNO,INSCOD,
     +             ALTCOD,SEGID,IZ,ID)


      RETURN
 1000 RETURN 1
 2000 RETURN 2
      END
C
C
C
      SUBROUTINE WBROOK(IUNIT,ISER,ATNAM,RESNAM,CHNNAM,IRESN,IS,
     *X,Y,Z,OCC,B,IZ)
C     ================================================================
C
C_BEGIN_RBROOK
C
C      This subroutine is obsolete and should be removed. May be 
C PROBLEM in that routine does not cater for IS.
C
C_END_RBROOK
C
C     .. Scalar Arguments ..
      REAL X,Y,Z,OCC,B
      INTEGER IUNIT,ISER,IRESN,IS,IZ
      CHARACTER ATNAM*4,RESNAM*4,CHNNAM*1
C     ..
C     .. Local Scalars ..
      CHARACTER*4 ID,SEGID,RESNO
      CHARACTER*1 INSCOD,ALTCOD
C     ..
C     .. Local Arrays ..
      REAL U(6)
C     ..
C     .. External Routines ..
      EXTERNAL XYZADVANCE,XYZATOM,XYZCOORD
C     ..
      SEGID = ' '
      ID = ' '
      INSCOD = ' '
      ALTCOD = ' '
      RESNO = ' '
      DO 10 I=1,6
        U(I) = 0.0
   10 CONTINUE


      CALL XYZADVANCE(IUNIT,0,0,*1000,*1000)
      CALL XYZCOORD(IUNIT,'O','U',X,Y,Z,OCC,B,U)
      CALL XYZATOM(IUNIT,ISER,ATNAM,RESNAM,CHNNAM,IRESN,RESNO,INSCOD,
     +             ALTCOD,SEGID,IZ,ID)
C
C---- This label is here for completeness but is not used (see XYZADVANCE).
C
 1000 CONTINUE
      RETURN
      END
C
C
C
C     SUBROUTINE PDBREAD(ISER,ATNAM,RESNAM,CHNNAM,IRESN,RESNO,
C    *X,Y,Z,OCC,U,IZ,SEGID,ID)
C     ========================================================
C
C_BEGIN_PDBREAD
C
C      The subroutine PDBREAD is used to read coordinates from a PDB
C format coordinate file. This routine should not be used stand alone 
C but through XYZADVANCE.
C 
C Parameters
C
C        ISER (O) Atom serial number
C       ATNAM (O) Atom name        (character*4 left justified)
C      RESNAM (O) Residue name     (character*4)
C      CHNNAM (O) Chain name       (character*1)
C       IRESN (O) Residue number as an integer
C       RESNO (O) Residue number   (character*4 or character*5)
C                 If character*5 then the 5th character will be the
C                 insertion code.
C           X (O) Coordinates (orthogonal angstrom coordinates as
C           Y (O)     "        stored)
C           Z (O)     "
C         OCC (O) Occupancy
C        U(6) (O) Temperature factor
C          IZ (O) Atomic number (returned as 7 from ambiguous atoms)
C          ID (O) Atomic ID related to atomic number + ionic state. 
C                 (character*4)
C
C  COMMON BLOCKS
C
C  COMMON /RBRKXX/IFCRYS,IFSCAL,ITYP,MATRIX
C
C      IFCRYS   .TRUE. IF 'CRYST1' CARD READ,  OTHERWISE .FALSE.
C      IFSCAL   .TRUE. IF 'SCALE' CARDS READ, OTHERWISE .FALSE.
C       ITYP    TYPE OF LAST CARD READ =1, 'CRYST1'
C                                      =2, 'SCALE'
C                                      =3, 'TER'
C                                      =4, 'ATOM'
C                                      =5, 'HETATM'
C     MATRIX    .TRUE. IF FRACT/ORTHOG MATRICES CALCULATED
C               .FALSE. IF NOT
C
C
C      BROOK    CHARACTER*1 ARRAY WHICH IS THE BUFFER FOR PDB FILES
C
C      COMMON/RBRKZZ/CELL(6),RR(3,3,6),VOL,CELLAS(6)
C
C       CELL    CELL DIMENSIONS FROM 'CRYST1' CARD IF READ
C               (CHECK IFCRYS)
C         RR    STANDARD ORTHOGONALISING MATRICES CALCULATED IF THE
C               'CRYST1' CARD WAS READ (CHECK IFCRYS)
C
C  COMMON /ORTHOG/RO(4,4),RF(4,4),NCODE
C
C        RO    ORTHOGONALISING MATRIX (ONLY SET IF 'CRYST1' OR 'SCALE'
C              CARDS PRESENT - CHECK 'MATRIX' FLAG)
C        RF    FRACTIONALISING MATRIX (ONLY SET IF 'CRYST1' OR 'SCALE'
C              CARDS PRESENT - CHECK 'MATRIX' FLAG)
C     NCODE    FLAG INDICATING SETTING FOUND, 0 IF NOT ONE THAT WAS
C              RECOGNISED
C
C_END_PDBREAD
C
C
C     .. Parameters ..
C     INTEGER MAXIATM, MAXIHATM
C     PARAMETER (MAXIATM=102,MAXIHATM=14)
C     ..
C     .. Arguments ..
C     REAL U(6),OCC,X,Y,Z
C     INTEGER IRESN,ISER,IZ
C     CHARACTER*(*) RESNO
C     CHARACTER ATNAM*4,CHNNAM*1,ID*4,RESNAM*4,SEGID*4
C     ..
C     .. Variables in Common ..
C     REAL CELL,CELLAS,RF,RO,RR,VOL
C     INTEGER ITYP,NCODE
C     CHARACTER BROOK*1,WBROOK*1,WBROOK1*1
C     LOGICAL IFCRYS,IFSCAL,IFTER,MATRIX,IFHDOUT,IFNEWCRYS
C     ..
C     .. Local Scalars ..
C     INTEGER I,II,J
C     CHARACTER*100 ERRLIN
C     CHARACTER*80 BROOKA
C     CHARACTER*4 IRTYPE
C     CHARACTER*2 IAA,IAT,IE
C     CHARACTER*1 ISP
C     ..
C     .. Local Arrays ..
C     INTEGER IU(6)
C     CHARACTER*40 ORTH(5)
C     CHARACTER*2 IATM(MAXIATM),IHATM(MAXIHATM)
C     ..
C     .. External Routines/Functions ..
C     INTEGER LENSTR
C     EXTERNAL CCPERR,CCPUPC,LENSTR
C     ..
C     .. Intrinsic Functions ..
C     INTRINSIC ABS
C     ..
C     .. Common Blocks ..
C     COMMON /RBRKXX/IFCRYS,IFSCAL,ITYP,MATRIX,IFHDOUT,IFNEWCRYS
C     COMMON/RBRKZZ/CELL(6),RR(3,3,6),VOL,CELLAS(6)
C     COMMON /ORTHOG/RO(4,4),RF(4,4),NCODE
C     COMMON /ORTHOGU/ ROU(4,4),RFU(4,4)
C     ..
C     .. Equivalences ..
C     EQUIVALENCE (IRTYPE,BROOK(1)),(IE,BROOK(5)),(BROOKA,BROOK(1))
C     ..
C     .. Save ..
C     SAVE /RBRKXX/,/RBRKZZ/,/ORTHOG/
C     ..
C     .. Data Statement ..
C     DATA IATM/' H','HE','LI','BE',' B',' C',' N',' O',' F','NE',
C    *          'NA','MG','AL','SI',' P',' S','CL','AR',' K','CA',
C    *          'SC','TI',' V','CR','MN','FE','CO','NI','CU','ZN',
C    *          'GA','GE','AS','SE','BR','KR','RB','SR',' Y','ZR',
C    *          'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN',
C    *          'SB','TE',' I','XE','CS','BA','LA','CE','PR','ND',
C    *          'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB',
C    *          'LU','HF','TA',' W','RE','OS','IR','PT','AU','HG',
C    *          'TL','PB','BI','PO','AT','RN','FR','RA','AC','TH',
C    *          'PA',' U','NP','PU','AM','CM','BK','CF','ES','FM',
C    *          ' D','AN'/
C     DATA IHATM/'0H','1H','2H','3H','4H','5H','6H','7H','8H','9H',
C    +           'HH','*H',"'H",'"H'/
C     DATA IAA/' A'/,ISP/' '/
C     DATA ORTH/'A // XO, C* // ZO (STANDARD PDB)',
C    *          'B // XO, A* // ZO',
C    *          'C // XO, B* // ZO',
C    *          'HEX A+B // XO, C* // ZO',
C    *          'A* // XO, C // ZO (ROLLETT)'/
C      DATA ITYPE/'CRYS','SCAL','TER ','ATOM','HETA','ANIS','END'/
C
C
C     IFTER=.FALSE.
C
C---- Atom/hetatm card processing
C
C     IF (IRTYPE.EQ.'ATOM' .OR. IRTYPE.EQ.'HETA' .OR.
C    +    IRTYPE.EQ.'ANIS' .OR. IRTYPE.EQ.'TER ') THEN
C       IF (IRTYPE.EQ.'TER ') THEN
C
C---- 'ter' card found
C
C         ITYP=3
C         IFTER=.TRUE.
C         GO TO 450
C       ENDIF
C       IF(BROOK(13).EQ.ISP) then
C
C----   ATNAM (O) Atom name  (character*4 left justified)
C
C         ATNAM=BROOK(14)//BROOK(15)//BROOK(16)//' '
C
C---- atnam should never had ALTCODE added to it 
C
C        else
C         ATNAM=BROOK(13)//BROOK(14)//BROOK(15)//BROOK(16)
C        end if
C
C
C450    READ(BROOKA,1006)ISER,IRESN
C       RESNAM=BROOK(18)//BROOK(19)//BROOK(20)//ISP
C       RESNO=BROOK(23)//BROOK(24)//BROOK(25)//BROOK(26)
C       IF(LEN(RESNO).GT.4)RESNO(5:5)=BROOK(27)
C       CHNNAM=BROOK(22)
C       IF(IFTER)GO TO 500
C       SEGID = BROOK(73)//BROOK(74)//BROOK(75)//BROOK(76)
C
C---- Find atomic number and ID, ID can be kept in columns 77-80
C
C       ID = BROOK(77)//BROOK(78)//BROOK(79)//BROOK(80)
C       CALL CCPUPC(ID)
C       IAT=BROOK(13)//BROOK(14)
C       CALL CCPUPC(IAT)
C
C---- Fast initial check for C, O, N or H
C       II = 0
C       IF (ID(1:4) .NE. '    ') THEN
C         IF (ID(1:2) .EQ. IATM(6)) THEN
C           II = 6
C           GOTO 480
C         ENDIF
C         IF (ID(1:2) .EQ. IATM(7)) THEN
C           II = 7
C           GOTO 480
C         ENDIF
C         IF (ID(1:2) .EQ. IATM(8)) THEN
C           II = 8
C           GOTO 480
C         ENDIF
C         IF (ID(1:2) .EQ. IATM(1)) THEN
C           II = 1
C           GOTO 480
C         ENDIF
C
C---- Must be a different element - check against all
C     possibilities, which is slower
C         DO 452 I=1,MAXIATM
C           IF (ID(1:2) .EQ. IATM(I)) THEN
C             II = I
C             GOTO 480
C           ENDIF
C452      CONTINUE
C       END IF 
C
C     If no ID match then make sure it is reset to be empty
C
C       ID = '    '
C
C     Check against first characters of atom name:
C
C---- Initial fast check against C, O, N or H
C       IF (IAT.EQ.IATM(6)) THEN
C         II = 6
C         GO TO 480
C       ENDIF
C       IF (IAT.EQ.IATM(7)) THEN
C         II = 7
C         GO TO 480
C       ENDIF
C       IF (IAT.EQ.IATM(8)) THEN
C         II = 8
C         GO TO 480
C       ENDIF
C       IF (IAT.EQ.IATM(1)) THEN
C         II = 1
C         GO TO 480
C       ENDIF
C
C---- Could be a hydrogen? Check for things like 0H, HH, etc
C       II=1
C       DO 454 J=1,MAXIHATM
C         IF (IAT.EQ.IHATM(J)) GO TO 480
C454    CONTINUE
C
C---- Must be a different element - check everything else
C       DO 456 I=1,MAXIATM
C         IF (IAT.EQ.IATM(I)) THEN
C           II = I
C           Should issue a warning if AC or AN to the effect that this
C           is based on ambigious input and should be checked
C           IF (II.EQ.89 .OR. II.EQ.102) THEN
C             WRITE(ERRLIN,2002)ATNAM,RESNAM,RESNO(1:4),IATM(II)
C             CALL CCPERR(2,ERRLIN)
C           END IF
C           GO TO 480
C         ENDIF
C456    CONTINUE
C
C---- No match for ID to anything using the first 2 characters
C     of the atom name
C
C  If the atom name begins with " A" set the atom_type to N
C  " A" is an ambigious atom name so presumably using N is
C  just a default?
C  Otherwise it's completely unknown
C
C       II=0
C       IF(IAT.EQ.IAA) II=7
C
C  Ambigious name...
C       IF (II .EQ. 7) THEN
C         WRITE(ERRLIN,2002)ATNAM,RESNAM,RESNO(1:4),IATM(II)
C         CALL CCPERR(2,ERRLIN)
C       ELSE
C
C---- II is zero so try some other tricks to get a match
C
C     Desperate measure: try to match second character only - 
C     This will deal with the NO7* horrors..
C         IF (IAT(1:1).NE.' ') THEN
C           IAT = ' '//BROOK(14)
C           DO I=1,MAXIATM
C             IF (IAT.EQ.IATM(I)) THEN
C               II = I
C               Issue a warning about this match since it
C               is based on incomplete data and assumptions
C               WRITE(ERRLIN,2002)ATNAM,RESNAM,RESNO(1:4),IATM(II)
C               CALL CCPERR(2,ERRLIN)
C               GO TO 480
C             ENDIF
C           END DO
C         ENDIF
C       END IF
C
C Still completely unrecognised... give up
C
C       IF (II .EQ. 0) THEN
C         WRITE(ERRLIN,2001)ATNAM,RESNAM,RESNO(1:4)
C         CALL CCPERR(2,ERRLIN)
C       END IF
C
C---- Atom number decided
C480    IZ=II
C
C        IF (IZ .EQ. 0) THEN
C          ID = ' '
C       ELSE 
C
C---- Keep the ionic state if valid, OR from atom name.
C
C         IF (ID(1:1) .EQ. ' ') THEN
C           IF (ATNAM(3:3).EQ.'+' .OR. ATNAM(3:3).EQ.'-') 
C    +                                       ID(3:4) = ATNAM(3:4)
C         ELSE
C           IF (ID(3:3).NE.'+' .AND. ID(3:3).NE.'-') ID(3:4) = '  '
C         ENDIF
C         ID(1:2) = IATM(IZ)
C       ENDIF
C
C---- Put elment ID into output buffer.
C
C       DO 485 J=1,4
C         WBROOK(76+J) = ID(J:J)
C         WBROOK1(76+J) = ID(J:J)
C485    CONTINUE
C        IF (IRTYPE .EQ. 'ATOM'.or.IRTYPE.EQ.'HETA') THEN
C  This is the ONLY flag that you have read a ATOM or HETATOM card..
C         DO 40 I=1,6
C           U(I) = 0.0
C  40     CONTINUE
C           IF (IRTYPE.EQ.'ATOM') ITYP=4
C           IF (IRTYPE.EQ.'HETA') ITYP=5
C         READ(BROOKA,1005)X,Y,Z,OCC,U(1)
C
C---- AnisoU cards
C
C       ELSE IF (IRTYPE .EQ. 'ANIS') THEN
C
C       READ(BROOKA,1010)IU(1),IU(2),IU(3),IU(4),IU(5),IU(6)
C       DO 510 I=1,6
C         U(I) = IU(I)/1.0E+04
C510    CONTINUE
C  Get rid of this, sometimes useful to know xyz 
C   use values of U(i) to check for ANISOU 
C       X = 0.0
C       Y = 0.0
C       Z = 0.0
C       ENDIF
C
C       RETURN        
C     ENDIF

C500  RETURN
C
C---- Format statements
C
C1005  FORMAT(30X,3F8.3,2F6.2)
C1006  FORMAT(6X,I5,11X,I4)
C1010  FORMAT(28X,6I7)
C2001  FORMAT(' *UNKNOWN ATOMIC FORMFACTOR ',A4,' IN ',A4,1X,A4,'*')
C2002  FORMAT(' *AMBIGUOUS ATOMIC FORMFACTOR ',A4,' IN ',A4,1X,A4,
C    +     ' ASSIGNED AS ',A2,' *')
C     END
C
C
C
      SUBROUTINE RBFRAC2(A,B,C,AL,BE,GA,ARGNCODE)
C     ===========================================
C
C_BEGIN_RBFRAC2
C
C
C This subroutine is used to calculate the default transformation
C matrices between orthogonal angstrom and fractional coordinates
C The sensible name is for Phil, as RBFRAC2 was changed from the original.
C
C
C PARAMETERS
C
C    A,B,C,AL,BE,GA (I) (REAL)     CELL PARAMETERS IN ANGSTROMS AND DEGREES
C    ARGNCODE       (I) (INTEGER)  ORTHOGONALISATION CODE 1-6
C
C_END_RBFRAC2
C
C     .. Arguments ..
      REAL A,B,C,AL,BE,GA
      INTEGER ARGNCODE
C     ..
C     .. Variables in Common ..
      REAL CELL,CELLAS,RF,RO,RR,VOL
      INTEGER ITYP,NCODE
      LOGICAL IFCRYS,IFSCAL,MATRIX,IFHDOUT,IFNEWCRYS
C     ..
C     .. External Routines ..
      EXTERNAL CCPERR,RBRINV
C     ..
C     .. Common Blocks ..
      COMMON /RBRKXX/IFCRYS,IFSCAL,ITYP,MATRIX,IFHDOUT,IFNEWCRYS
      COMMON /ORTHOG/RO(4,4),RF(4,4),NCODE
      COMMON /RBRKZZ/CELL(6),RR(3,3,6),VOL,CELLAS(6)
C     ..
C     .. Save Statement ..
      SAVE /RBRKXX/,/ORTHOG/,/RBRKZZ/
C     ..
      IF (ARGNCODE.NE.0) NCODE = ARGNCODE
C
C---- Calculate matrices
C

      DO 10 I=1,4
      DO 10 J=1,4
      RO(I,J)=0
      RF(I,J)=0
10    CONTINUE
C
C
      RO(4,4)=1.0
      RF(4,4)=1.0
      CELL(1)=A
      CELL(2)=B
      CELL(3)=C
      CELL(4)=AL
      CELL(5)=BE
      CELL(6)=GA
      CALL RBFROR
C
C
      DO 20 I=1,3
      DO 20 J=1,3
      RO(I,J)=RR(I,J,1)
20    CONTINUE
C
C
      CALL RBRINV(RO,RF)
      MATRIX=.TRUE.
      CALL CCPERR(4,' ')
      CALL CCPERR(4,
     +     ' STANDARD PDB COORDINATE SETTING WILL BE ASSUMED')
      CALL CCPERR(4,
     +     ' IF NO SCALE CARDS PRESENT  IN  INPUT  COORDINATE  FILE')
      CALL CCPERR(4,' ')

      RETURN
      END
C
C
C
      SUBROUTINE RBFRAC(A,B,C,AL,BE,GA,MSG)
C     =====================================
C
C_BEGIN_RBFRAC
C
C	This routine is obsolete and should be removed.
C
C_END_RBFRAC
C
C     .. Scalar Arguments ..
      REAL A,B,C,AL,BE,GA
      INTEGER MSG
C     ..
C     .. External Routines ..
      EXTERNAL RBFRAC2
C     ..
      CALL RBFRAC2(A,B,C,AL,BE,GA,1)

      RETURN
      END
C
C
C
      SUBROUTINE RBRORF(ROO,RFF)
C     ==========================
C
C_BEGIN_RBRORF
C
C	This routine is obsolete and should be removed.
C
C      SUBROUTINE RBRORF(ROO,RFF)
C
C     Subroutine to  fill or return RF (fractionalising) and Ro
C     (orthogonalising) matrices. 
C
C PARAMETERS
C
C          ROO (I) (REAL(4,4))  4*4 MATRIX TO BE INVERTED
C          RFF (O) (REAL(4,4))  INVERSE MATRIX
C
C common blocks
C
C
C
      DIMENSION ROO(4,4),RFF(4,4)
C
      LCODE = 0
      CALL RBRORF2(ROO,RFF,LCODE)
      END
C
C
C
      SUBROUTINE RBRORF2(ROO,RFF,LCODE)
C     ==========================
C
C_BEGIN_RBRORF
C
C      SUBROUTINE RBRORF2(ROO,RFF,LCODE)
C
C     Subroutine to  fill or return RF (fractionalising) and Ro
C     (orthogonalising) matrices. 
C
C PARAMETERS
C
C          ROO (I) (REAL(4,4))  4*4 MATRIX TO BE INVERTED
C          RFF (O) (REAL(4,4))  INVERSE MATRIX
C
C common blocks
C
C      COMMON /RBRKXX/IFCRYS,IFSCAL,ITYP,MATRIX
C      COMMON /ORTHOG/RO(4,4),RF(4,4),NCODE
C
C_END_RBRORF

      LOGICAL IFCRYS,IFSCAL,MATRIX,IFHDOUT,IFNEWCRYS
      COMMON /RBRKXX/IFCRYS,IFSCAL,ITYP,MATRIX,IFHDOUT,IFNEWCRYS
      COMMON /ORTHOG/RO(4,4),RF(4,4),NCODE
      SAVE /ORTHOG/, /RBRKXX/
C
C
      DIMENSION ROO(4,4),RFF(4,4)
C
C---- Get cofactors of 'a' in array 'c'
C
      IF(ROO(1,1) .LE.0.0000000001) THEN 
        DO 40 II=1,4
        DO 30 JJ=1,4
        RFF(II,JJ)=RF(II,JJ)
        ROO(II,JJ)=RO(II,JJ)
30      CONTINUE
40      CONTINUE
        LCODE = NCODE
        RETURN
      END IF
C  FILL...
      IF(ROO(1,1) .GT.0.0000000001) THEN 
        DO 50 II=1,4
        DO 60 JJ=1,4
        RF(II,JJ)=RFF(II,JJ)
        RO(II,JJ)=ROO(II,JJ)
60      CONTINUE
50      CONTINUE
        NCODE = LCODE
        MATRIX = .TRUE.
        RETURN
      END IF
      END
C
C
      SUBROUTINE RBRINV(A,AI)
C     =======================
C
C_BEGIN_RBRINV
C
C      SUBROUTINE RBRINV(A,AI)
C
C
C Subroutine to invert 4*4 matrices for conversion between
C fractional and orthogonal axes. 
C
C
C PARAMETERS
C
C           A (I) (REAL(4,4))  4*4 MATRIX TO BE INVERTED
C          AI (O) (REAL(4,4))  INVERSE MATRIX
C
C_END_RBRINV
C
      REAL A(4,4),AI(4,4),C(4,4),X(3,3)
C
C---- Get cofactors of 'a' in array 'c'
C
      DO 40 II=1,4
      DO 30 JJ=1,4
      I=0
      DO 20 I1=1,4
      IF(I1.EQ.II)GO TO 20
      I=I+1
      J=0
      DO 10 J1=1,4
      IF(J1.EQ.JJ)GO TO 10
      J=J+1
      X(I,J)=A(I1,J1)
10    CONTINUE
20    CONTINUE
      AM=X(1,1)*X(2,2)*X(3,3)-X(1,1)*X(2,3)*X(3,2)+X(1,2)*X(2,3)*X(3,1)
     *  -X(1,2)*X(2,1)*X(3,3)+X(1,3)*X(2,1)*X(3,2)-X(1,3)*X(2,2)*X(3,1)
      C(II,JJ)=(-1)**(II+JJ)*AM
30    CONTINUE
40    CONTINUE
C
C---- Calculate determinant
C
      D=0
      DO 50 I=1,4
      D=D+A(I,1)*C(I,1)
50    CONTINUE
C
C---- Get inverse matrix
C
      DO 70 I=1,4
      DO 60 J=1,4
      AI(I,J)=C(J,I)/D
60    CONTINUE
70    CONTINUE
      RETURN
      END
C
C
C
C
      SUBROUTINE RBFROR
C     =================
C
C_BEGIN_RBFROR
C
C      SUBROUTINE RBFROR
C
C THIS SUBROUTINE CALCULATES MATRICES FOR STANDARD ORTHOGONALISATIONS
c   and cell volume
C
C  this generates the various orthogonalising matrices
C     ' NCODE =1 -  ORTHOG AXES ARE DEFINED TO HAVE'
C                    A PARALLEL TO XO   CSTAR PARALLEL TO ZO'
C     ' NCODE =2 -  ORTHOG AXES ARE DEFINED TO HAVE'
C     '               B PARALLEL TO XO   ASTAR PARALLEL TO ZO'
C     ' NCODE =3 -  ORTHOG AXES ARE DEFINED TO HAVE'
C     '               C PARALLEL TO XO   BSTAR PARALLEL TO ZO'
C     ' NCODE =4 -  ORTHOG AXES ARE DEFINED TO HAVE'
C     '         HEX A+B PARALLEL TO XO   CSTAR PARALLEL TO ZO'
C     ' NCODE =5 -  ORTHOG AXES ARE DEFINED TO HAVE'
C     '           ASTAR PARALLEL TO XO       C PARALLEL TO ZO'
C     ' NCODE =6 -  ORTHOG AXES ARE DEFINED TO HAVE'
C                    A  PARALLEL TO XO   BSTAR PARALLEL TO YO'
C
C   SET UP MATRICES TO ORTHOGONALISE H K L AND X Y Z FOR THIS CELL.
C
C Common Blocks
C
C     .. Scalar Arguments ..
      REAL VOLL
C     ..
C     .. Array Arguments ..
      REAL CEL(6),RRR(3,3,6)
C
C      COMMON/RBRKZZ/CELL(6),RR(3,3,6),VOL,CELLAS(6)
C      COMMON /RBREC/AC(6)
C
C_END_RBFROR
C
      COMMON/RBRKZZ/CELL(6),RR(3,3,6),VOL,CELLAS(6)
      COMMON /RBREC/AC(6)
      SAVE /RBRKZZ/, /RBREC/
C
C---- Initialisations
C
        DO I = 1,6
         CEL(I)=CELL(I)
        END DO
        VOLL = VOL
      CALL RBFRO1(CEL,VOLL,RRR)
C
      RETURN
      END
C
C
C     ===============================
      SUBROUTINE RBFRO1(CEL,VOLL,RRR)
C     ===============================
C
C_BEGIN_RBFRO1
C
C      SUBROUTINE RBFRO1(CEL,VOLL,RRR)
C
C---- This subroutine is a duplicate of rbfror with a different call.
C
C PARAMETERS
C
C   CEL  (I) (REAL(6))     Cell dimensions
C   VOLL (O) (REAL)        Cell volume
C   RRR  (O) (REAL(3,3,6)) Standard orthogonisational matrices
C
C_END_RBFRO1
C
C THIS SUBROUTINE CALCULATES MATRICES FOR STANDARD ORTHOGONALISATIONS
c   and cell volume
C
C  this generates the various orthogonalising matrices
C     ' NCODE =1 -  ORTHOG AXES ARE DEFINED TO HAVE'
C                    A PARALLEL TO XO   CSTAR PARALLEL TO ZO'
C     ' NCODE =2 -  ORTHOG AXES ARE DEFINED TO HAVE'
C     '               B PARALLEL TO XO   ASTAR PARALLEL TO ZO'
C     ' NCODE =3 -  ORTHOG AXES ARE DEFINED TO HAVE'
C     '               C PARALLEL TO XO   BSTAR PARALLEL TO ZO'
C     ' NCODE =4 -  ORTHOG AXES ARE DEFINED TO HAVE'
C     '         HEX A+B PARALLEL TO XO   CSTAR PARALLEL TO ZO'
C     ' NCODE =5 -  ORTHOG AXES ARE DEFINED TO HAVE'
C     '           ASTAR PARALLEL TO XO       C PARALLEL TO ZO'
C     ' NCODE =6 -  ORTHOG AXES ARE DEFINED TO HAVE'
C                    A  PARALLEL TO XO   BSTAR PARALLEL TO YO'
C
C   SET UP MATRICES TO ORTHOGONALISE H K L AND X Y Z FOR THIS CELL.
C
C     .. Scalar Arguments ..
      REAL VOLL
C     ..
C     .. Array Arguments ..
      REAL CEL(6),RRR(3,3,6)
C     ..
C     .. Scalars in Common ..
      REAL VOL
C     ..
C     .. Arrays in Common ..
      REAL AC,CELL,CELLAS,RR
C     ..
C     .. Local Scalars ..
      REAL A,ALPH,ALPHAS,AS,B,BET,BETAS,BS,C,CONV,COSA,COSAS,COSB,COSBS,
     +     COSG,COSGS,CS,FCT,GAMM,GAMMAS,SINA,SINAS,SINB,SINBS,SING,
     +     SINGS,SUM,V
      INTEGER I,J,K,N,NCODE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN2,COS,SIN,SQRT
C     ..
C     .. Common blocks ..
      COMMON /RBREC/AC(6)
      COMMON /RBRKZZ/CELL(6),RR(3,3,6),VOL,CELLAS(6)
C     ..
C     .. Save statement ..
      SAVE /RBRKZZ/, /RBREC/
C     ..
C
C---- Initialisations
C
      CELDEL = 0.0
      IF (CEL(1).GT.0.0) THEN
        IF (CELL(1).GT.0.0) THEN
          IWARN = 0
          DO 101 I = 1,6
            CELDEL = ABS(CEL(I)-CELL(I))/CEL(I)
            IF (CELDEL.GT.0.01) IWARN = 1
 101      CONTINUE
          IF(IWARN.NE.0) WRITE(6,9876)CELL,CEL
9876      FORMAT(' Inconsistency in Cell Dimensions - replacing old:',
     +      /,' Old cell:',3X,6F10.5,/,' New cell:',3X,6F10.5)
        ENDIF
        DO 10 I = 1,6
          CELL(I) = CEL(I)
          IF (CELL(I).EQ.0.0) call ccperr(1,
     +' **** Incorrect (0.0) CELL element in  subroutine rbfro1?? ****')
   10   CONTINUE
      ENDIF
C
C
      CONV = 3.14159/180.0
      FCT = 8.0*3.14159*3.14159
      ALPH = CELL(4)*CONV
      BET = CELL(5)*CONV
      GAMM = CELL(6)*CONV
      SUM = (ALPH+BET+GAMM)*0.5
      V = SQRT(SIN(SUM-ALPH)*SIN(SUM-BET)*SIN(SUM-GAMM)*SIN(SUM))
      VOL = 2.0*CELL(1)*CELL(2)*CELL(3)*V
      SINA = SIN(ALPH)
      COSA = COS(ALPH)
      SINB = SIN(BET)
      COSB = COS(BET)
      SING = SIN(GAMM)
      COSG = COS(GAMM)
      COSAS = (COSG*COSB-COSA)/ (SINB*SING)
      SINAS = SQRT(1.0-COSAS*COSAS)
      COSBS = (COSA*COSG-COSB)/ (SINA*SING)
      SINBS = SQRT(1.0-COSBS*COSBS)
      COSGS = (COSA*COSB-COSG)/ (SINA*SINB)
      SINGS = SQRT(1.0-COSGS*COSGS)
      A = CELL(1)
      B = CELL(2)
      C = CELL(3)
      AS = B*C*SINA/VOL
      BS = C*A*SINB/VOL
      CS = A*B*SING/VOL
      ALPHAS = ATAN2(SINAS,COSAS)/CONV
      BETAS  = ATAN2(SINBS,COSBS)/CONV
      GAMMAS = ATAN2(SINGS,COSGS)/CONV
      CELLAS(1) = AS
      CELLAS(2) = BS
      CELLAS(3) = CS
      CELLAS(4) = ALPHAS
      CELLAS(5) = BETAS
      CELLAS(6) = GAMMAS
C
C---- Set useful things for calculating dstar
C
      AC(1) = AS*AS
      AC(2) = BS*BS
      AC(3) = CS*CS
      AC(4) = 2.0*BS*CS*COSAS
      AC(5) = 2.0*CS*AS*COSBS
      AC(6) = 2.0*AS*BS*COSGS
C
C---- Zero matrices
C
      DO 40 N = 1,6
        DO 30 I = 1,3
          DO 20 J = 1,3
            RR(I,J,N) = 0.0
   20     CONTINUE
   30   CONTINUE
   40 CONTINUE
C
C---- Calculate matrices
C
C---- XO along a  Zo along c*
C
      NCODE = 1
      RR(1,1,NCODE) = A
      RR(1,2,NCODE) = B*COSG
      RR(1,3,NCODE) = C*COSB
      RR(2,2,NCODE) = B*SING
      RR(2,3,NCODE) = -C*SINB*COSAS
      RR(3,3,NCODE) = C*SINB*SINAS
C
C---- XO along b  Zo along a*
C
      NCODE = 2
      RR(1,1,NCODE) = A*COSG
      RR(1,2,NCODE) = B
      RR(1,3,NCODE) = C*COSA
      RR(2,1,NCODE) = -A*SING*COSBS
      RR(2,3,NCODE) = C*SINA
      RR(3,1,NCODE) = A*SING*SINBS
C
C---- XO along c  Zo along b*
C
      NCODE = 3
      RR(1,1,NCODE) = A*COSB
      RR(1,2,NCODE) = B*COSA
      RR(1,3,NCODE) = C
      RR(2,1,NCODE) = A*SINB
      RR(2,2,NCODE) = -B*SINA*COSGS
      RR(3,2,NCODE) = B*SINA*SINGS
C
C---- trigonal only - XO along a+b  YO alon a-b  Zo along c*
C
      NCODE = 4
      RR(1,1,NCODE) = A/2.0
      RR(1,2,NCODE) = A/2.0
      RR(2,1,NCODE) = -A*SING
      RR(2,2,NCODE) = A*SING
      RR(3,3,NCODE) = C
C
C---- XO along a*   ZO along c
C
      NCODE = 5
      RR(1,1,NCODE) = A*SINB*SINGS
      RR(2,1,NCODE) = -A*SINB*COSGS
      RR(2,2,NCODE) = B*SINA
      RR(3,1,NCODE) = A*COSB
      RR(3,2,NCODE) = B*COSA
      RR(3,3,NCODE) = C
C
C---- Grr*! to  Gerard Bricogne - his setting for P1 in SKEW.
C     XO along a  Yo along b*
C
      NCODE = 6
      RR(1,1,NCODE) = A
      RR(1,2,NCODE) = B*COSG
      RR(1,3,NCODE) = C*COSB
      RR(2,2,NCODE) = B*SING*SINAS
      RR(3,2,NCODE) = -B*SING*COSAS
      RR(3,3,NCODE) = C*SINB
C
C---- copy rr(...) into rrr(...)
C
      DO 11 I=1,6
        CEL(I) = CELL(I)
11    CONTINUE
      DO 70 I = 1,6
        DO 60 J = 1,3
          DO 50 K = 1,3
            RRR(K,J,I) = RR(K,J,I)
   50     CONTINUE
   60   CONTINUE
   70 CONTINUE
C
C
      VOLL = VOL
      IF (CELDEL.GT.0.01) VOLL = -VOL
C
C
      END
C
C
      SUBROUTINE CVFRAC2(X,Y,Z,XX,YY,ZZ,IFLAG)
C     ============================================
C
C_BEGIN_CVFRAC2
C
C      This subroutine is used to convert between the stored  orthogonal  and
C fractional coordinates using the  matrices  set  up  in  the  common  block
C /ORTHOG/ by the subroutine XYZOPEN. If no matrices have been set up then the
C program will stop with an error message.
C                                         
C Call:  CALL CVFRAC2(X,Y,Z,XX,YY,ZZ,IFLAG)
C                                             
C Arguments:
C            
C       X (I)   (REAL)  Input coordinates.
C       Y (I)   (REAL)       "
C       Z (I)   (REAL)       "
C      XX (O)   (REAL)  Output coordinates.
C      YY (O)   (REAL)       "
C      ZZ (O)   (REAL)       "
C   IFLAG (I) (INTEGER)  Flag =0, Convert coordinates from fractional to orthogonal
C                             =1, Convert coordinates from orthogonal to fractional
C
C_END_CVFRAC2
C
      LOGICAL IFCRYS,IFSCAL,MATRIX,IFHDOUT,IFNEWCRYS
      COMMON /RBRKXX/IFCRYS,IFSCAL,ITYP,MATRIX,IFHDOUT,IFNEWCRYS
      COMMON /ORTHOG/RO(4,4),RF(4,4),NCODE
      SAVE /RBRKXX/,/ORTHOG/
C
C---- Check that matrices set up
C
      IF(.NOT.MATRIX)GO TO 800
C
C---- Perform transformation
C
      IF(IFLAG.NE.0)GO TO 10
      XX=RO(1,1)*X + RO(1,2)*Y +RO(1,3)*Z +RO(1,4)
      YY=RO(2,1)*X + RO(2,2)*Y +RO(2,3)*Z +RO(2,4)
      ZZ=RO(3,1)*X + RO(3,2)*Y +RO(3,3)*Z +RO(3,4)
      RETURN
10    XX=RF(1,1)*X + RF(1,2)*Y +RF(1,3)*Z +RF(1,4)
      YY=RF(2,1)*X + RF(2,2)*Y +RF(2,3)*Z +RF(2,4)
      ZZ=RF(3,1)*X + RF(3,2)*Y +RF(3,3)*Z +RF(3,4)
      RETURN
C
C---- Error condition
C
800   CALL CCPERR(4,' **FRACTIONAL/ORTHOGONAL MATRICES NOT SET UP**')
      call ccperr(1,' No knowledge of input orthogonalisation')
C
C---- Format statements
C
      END
C
C
C
      SUBROUTINE CVFRAC(X,Y,Z,XX,YY,ZZ,IFLAG,MSG)
C     ===========================================
C
C_BEGIN_CVFRAC
C
C	Another silly obsolete routine that really should be deleted.
C MSG value is in fact useless. Library output controlled by CCPERR.
C
C_END_CVFRAC
C
C     .. Scalar Arguments ..
      REAL X,Y,Z,XX,YY,ZZ
      INTEGER IFLAG,MSG
C     ..
C     .. External Routines ..
      EXTERNAL CVFRAC2
C     ..
      CALL CVFRAC2(X,Y,Z,XX,YY,ZZ,IFLAG)

      RETURN
      END
C
C
C
      SUBROUTINE CVANISOB(B,IFLAG)
C     THIS SUBROUTINE SHOULD NOT BE USED
C     SEE CVANISOU
      REAL B(6)
      INTEGER IFLAG
      EXTERNAL CVANISOu, CCPERR
      WRITE (6,*) 'ERR: THIS PROGRAM USES S/R CVANISOB'
      WRITE (6,*) 'ERR: IT SHOULD NOT USE THIS ROUTINE'
      WRITE (6,*) 'ERR: CVANISOU IS CALLED AUTOMATICALLY'
      CALL CCPERR(2, 'CHANGE YOUR CODE')
      CALL CVANISOU(B, IFLAG)
      RETURN
      END
C
C
C
C
      SUBROUTINE CVANISOU(U,IFLAG)
C     ============================
C
C_BEGIN_CVANISOU
C
C      This subroutine is used to convert between crystallographic bs and 
C orthogonal Us or the other way round. The orthogonal matrices are 
C required, if no matrices have been set up then the program will stop 
C with an error message. The temperature factors are defined below;
C
C  PDB files contain anisotropic temperature factors as orthogonal Us.
C The anisotropic temperature factors can be input/output to this routine 
C  as orthogonal or as crystallographic Us. 
C  
C  Shelx defines Uf to calculate temperature factor as:
C T(aniso_Uf) = exp (-2PI**2 ( (h*ast)**2 Uf_11 + (k*bst)**2 Uf_22 + ... 
C                            + 2hk*ast*bst*Uf_12 +..)
C
C   Note:   Uo_ji == Uo_ij and  Uf_ji == Uf_ij.
C
C  [Uo_ij] listed on ANISOU card satisfy  the relationship:
C  [Uo_ij] =   [RFu]-1 [Uf_ij] {[RFu]-1}T   
C
C        where [Rfu] is the normalised [Rf] matrix read from the SCALEi cards.
C        see code.   [ROu] ==  [RFu]-1
C  Hence:
C  [Uf_ij] =   [RFu]   [Uo_ij] {[RFu]  }T   
C
C T(aniso_Uo) = U(11)*H**2 + U(22)*K**2 + 2*U(12)*H*K + ...
C where H,K,L are orthogonal reciprocal lattice indecies. ( EJD: I think????)
C
C Biso     = 8*PI**2 (Uo_11 + Uo_22 + Uo_33) / 3.0
C
C   [Uf(symm_j)] = [Symm_j] [Uf] [Symm_j]T
C
C 
C Arguments:
C            
C    U(6) (I/O  (REAL)  Input coordinates.
C   IFLAG (I) (INTEGER)  Flag =0, Convert coordinates from fract. to orthog.
C                             =1, Convert coordinates from orthog. to fract.
C
C_END_CVANISOU
C
C     .. Arguments ..
      REAL U(6)
      INTEGER IFLAG
C     ..
C     .. Variables in Common ..
      REAL RF,RO
      INTEGER NCODE
      LOGICAL IFCRYS,IFSCAL,MATRIX,IFHDOUT,IFNEWCRYS
C     ..
C     .. Local Variables ..
      INTEGER I,J
      REAL TWOPI2
      REAL A(3,3),AT(3,3),TMP(3,3),TMPMAT(3,3)
C     ..
C     .. External Routines ..
      EXTERNAL CCPERR,MATMUL
C     ..
C     .. Common Blocks ..
      COMMON /RBRKXX/IFCRYS,IFSCAL,ITYP,MATRIX,IFHDOUT,IFNEWCRYS
      COMMON /ORTHOG/RO(4,4),RF(4,4),NCODE
      COMMON /ORTHOGU/ ROU(4,4),RFU(4,4)
C     ..
C     .. Save Statement ..
      SAVE /RBRKXX/,/ORTHOG/,/ORTHOGU/
C     ..
C     .. Data Statements ...
      DATA TWOPI2 /19.739209/
C     ..
C
C---- Check that matrices set up
C
      IF(.NOT.MATRIX)GO TO 800
C
C---- Perform transformation
C
        TMP(1,1)=U(1)
        TMP(2,2)=U(2)
        TMP(3,3)=U(3)
        TMP(1,2)=U(4)
        TMP(2,1)=U(4)
        TMP(1,3)=U(5)
        TMP(3,1)=U(5)
        TMP(2,3)=U(6)
        TMP(3,2)=U(6)
C
C   IFLAG (I) (INTEGER)  Flag =0, Convert coordinates from fract. to orthog.
C   IFLAG (I) (INTEGER)  Flag =1, Convert coordinates from orthog. to fract.
C
      IF (IFLAG .EQ. 0) THEN
        DO 10 I=1,3
          DO 10 J=1,3
            A(J,I)=ROU(J,I)
            AT(I,J)=ROU(J,I)
   10   CONTINUE
      ELSE
        DO 20 I=1,3
          DO 20 J=1,3
            A(J,I) = RFU(J,I)
            AT(I,J) = RFU(J,I)
   20   CONTINUE
      ENDIF
C
        CALL MATMUL(TMPMAT,TMP,AT)
        CALL MATMUL(TMP,A,TMPMAT)
        U(1) = TMP(1,1)
        U(2) = TMP(2,2)
        U(3) = TMP(3,3)
        U(4) = TMP(1,2)
        U(5) = TMP(1,3)
        U(6) = TMP(2,3)

      RETURN
C
C---- Error condition
C
800   CALL CCPERR(4,' **FRACTIONAL/ORTHOGONAL MATRICES NOT SET UP**')
      call ccperr(1,' No knowledge of input orthogonalisation')
C
C---- Format statements
C
      END
C
C
C
      SUBROUTINE RBCELL(CELLD,CVOL)
C     ============================
C
C_BEGIN_RBCELL
C
C      SUBROUTINE RBCELL(CELLD,CVOL)
C
C Returns cell dimensions and unit cell  volume.
C
C PARAMETERS
C     CELLD (O)  (REAL(6))  cell dimensions
C     CVOL (O)  (REAL)     cell volume
C
C Common blocks
C
C      COMMON/RBRKZZ/CELL(6),RR(3,3,6),VOL,CELLAS(6)
C
C_END_RBCELL
C
      COMMON/RBRKZZ/CELL(6),RR(3,3,6),VOL,CELLAS(6)
      REAL CELL,RR,VOL,CELLAS
      REAL CELLD(6), CVOL
      SAVE /RBRKZZ/
C
      CVOL = VOL
      DO 1 I=1,6
1       CELLD(I) = CELL(I) 
      END
C
C
C
      SUBROUTINE RBRCEL(RCEL,RVOL)
C     ============================
C
C_BEGIN_RBRCEL
C
C      SUBROUTINE RBRCEL(RCEL,RVOL)
C
C THIS SUBROUTINE RETURNS Reciprocal cell dimensions, and reciprocal
C                       unit cell  volume.
C
C PARAMETERS
C     RCEL (O)  (REAL(6)) reciprocal cell dimensions
C     RVOL (O)  (REAL)    reciprocal cell volume
C
C Common blocks
C
C      COMMON/RBRKZZ/CELL(6),RR(3,3,6),VOL,CELLAS(6)
C
C_END_RBRCEL
C
      COMMON/RBRKZZ/CELL(6),RR(3,3,6),VOL,CELLAS(6)
      REAL RCEL(6),RVOL
      REAL CELL,CELLAS,RR,VOL
      SAVE /RBRKZZ/
C
      IF (VOL.EQ.0.0) THEN
        RVOL = 0.0
      ELSE
        RVOL = 1.0/VOL
      ENDIF
      DO 1 I=1,6
1       RCEL(I) = CELLAS(I) 
      RETURN
      END
C
C
C
      SUBROUTINE RBSPGRP(SPGRP)
C     ============================
C
C_BEGIN_RBSPGRP
C
C      SUBROUTINE SUBROUTINE RBSPGRP(SPGRP)
C
C Returns spacegrpup from pdb
C
C PARAMETERS
C     SPGRP (O) (CHARACTER*15) 
C
C Common blocks
C
C      COMMON /RBRKSPGRP/BRKSPGRP
C
C_END_RBSPGRP
C
      CHARACTER SPGRP*(*)
      CHARACTER BRKSPGRP*15
      INTEGER ILEN,KLEN,J
      COMMON /RBRKSPGRP/BRKSPGRP
      SPGRP = BRKSPGRP
C
C Make sure that the returned name is left-justified
C
      ILEN = LENSTR(SPGRP)
      KLEN = ILEN
C
      DO J = 1,ILEN-1
        IF (SPGRP(1:1) .EQ. ' ') THEN
          SPGRP = SPGRP(2:KLEN)
          KLEN  = KLEN - 1
        END IF
      END DO
C
      END
C
C
C
      SUBROUTINE WBSPGRP(SPGRP)
C     =============================
C
C_BEGIN_WBSPGRP
C
C      SUBROUTINE WBSPGRP(SPGRP)
C
C Sets the internal spacegroup of a pdb file
C
C PARAMETERS
C     SPGRP (I) (CHARACTER*15)
C
C Common Blocks
C
C      COMMON /RBRKSPGRP/BRKSPGRP
C
C_END_WBSPGRP
C
C...  Local scalars
      CHARACTER SPGRP*(*)
C...  Common block
      CHARACTER BRKSPGRP*15
      COMMON /RBRKSPGRP/BRKSPGRP
      BRKSPGRP = SPGRP
C
      END
C
C
C
      SUBROUTINE RES3TO1(RESNM3,RESNM1)
C     ================================
C
C_BEGIN_RES3TO1
C
C      SUBROUTINE RES3TO1(RESNM3,RESNM1)
C
C       FIND 3 CHARACTER RESIDUE NAME FROM 1 CHARACTER CODE OR
C       FIND 1 CHARACTER RESIDUE NAME FROM 3 CHARACTER CODE.
C       SUBROUTINE IS CALLED WITH EITHER RESNM3 OR RESNM1 PREVIOUSLY 
C       ASSIGNED, AND THE OTHER IS ASSIGNED  HERE.
C 
C Parameters
C
C   RESNM3 (I/O)  CHAR*4    3 character residue name
C   RESNM1 (I/O)  CHAR*1    1 character residue name
C
C_END_RES3TO1
C
      CHARACTER*4 RESNM3
      CHARACTER*1 RESNM1
      CHARACTER*4 MAACD3(26)
      CHARACTER*1 MAACD1(26)
      DATA NAACID/26/
      DATA MAACD3/'ALA ','ARG ','ASN ','ASP ','CYS ','CYH ','GLN ',
     1 'GLU ','GLY ','HIS ','ILE ','LEU ','LYS ','MET ','PHE ','PRO ',
     2 'SER ','THR ','TRP ','TYR ','VAL ','HEM ','WAT ','SUL ','END ',
     3 'DUM '/
      DATA MAACD1/'A','R','N','D','C','C','Q',
     1 'E','G','H','I','L','K','M','F','P',
     2 'S','T','W','Y','V','X','O','U','Z','Z'/
C
C---- Routine to find one character amino acid name
C
        IF(RESNM3.NE.' ')THEN
          DO 1 I=1,NAACID
          IF(RESNM3.EQ.MAACD3(I)) GO TO 2
1         CONTINUE
        I=NAACID
2       RESNM1=MAACD1(I)
        RETURN
        ENDIF
C
C---- Routine to find three character amino acid name
C
        IF(RESNM1.NE.' ')THEN
          DO 11 I=1,NAACID
          IF(RESNM1.EQ.MAACD1(I)) GO TO 12
11        CONTINUE
        I=NAACID
12      RESNM3=MAACD3(I)
        RETURN
        ENDIF
      END
C
C
C
        SUBROUTINE RBRECIP(IH,IK,IL,S)
C       ==============================
C
C_BEGIN_BRECIP
C
C        SUBROUTINE RBRECIP(IH,IK,IL,S)
C
C---- This subroutine calculates 4SIN**2/L**2
C
C PARAMETERS
C         IH,IK,IL (I) (INTEGER)  reflection indices
C                S (O) (REAL)     4SIN**2/L**2
C
C_END_BRECIP
C
      COMMON /RBREC/AC(6)
      SAVE /RBREC/
C
      S = 
     .(AC(1)*IH*IH+AC(2)*IK*IK+AC(3)*IL*IL
     .+AC(4)*IK*IL+AC(5)*IL*IH+AC(6)*IH*IK)
      RETURN
      END


C     =====================================================
      SUBROUTINE SFREAD2(ID,NG,A,B,C,IWT,IELEC,CU,MO,Ifail)
C     =====================================================
C
C  Inputs: ID     atom identifier
C          This should match an atom type in the atomsf.lib
C          If an atom is identified as NE2+ say, characters are 
C          subtracted from the ID till a match is found, or there are 
C          no characters left. 
C          EG: Routine tests first NE2+, then NE2, then NE, then N.
C            All matching checks UPPER CASE strings.
C
C          NG     num. of gaussian approximations (2 or 5 (default))
C          IFILFF  .TRUE. if want to open the library file assigned
C                 to ATOMSF (default `atomsf.lib')
C
C  Output: A(4)   coefficient for structure factor calculation
C          B(4)   coefficient for structure factor calculation
C          C      coefficient for structure factor calculation
C          IWT    atomic weight
C          IELEC  number of electrons
C          CU(2)  delta F' and delta F'' for Cu
C          MO(2)  delta F' and delta F'' for Mo
C          Ifail  = -1 if atom not found at all
C                 =  0 OK
C                 =  1 for two gaussian case that does not exist
C
C 20/11/2000 C. Vonrhein
C     fixed problem with one character atom name matching
C     (S, I and B affected)
C
C     .. Scalar Arguments ..
      REAL C
      INTEGER IELEC,Ifail,IWT,NG
      CHARACTER ID*4,IDCHK*4,STRING*200
C     ..
C     .. Array Arguments ..
      REAL A(4),B(4),CU(2),MO(2)
C     ..
C     .. Local Scalars ..
      INTEGER NGauss,IOS
      CHARACTER ID2*6,IDIN*6
      LOGICAL OP
C     ..
C     .. External Subroutines ..
      EXTERNAL CCPDPN,CCPERR
C     ..

      IDCHK = ID
      CALL CCPUPC(IDCHK)
      ID2    = IDCHK//'  '
      NGauss = NG
C
      IF (NGauss.EQ.2) THEN
        ID2(6:6) = '2'
      ELSE
        NGauss = 5
      END IF
C
C---- Check to open file
C
      INQUIRE (UNIT=45, OPENED=OP, IOSTAT=IOS)
      IF (IOS .NE. 0) CALL CCPERR(1,'Error opening ATOMSF file')

      IF (.NOT.OP) THEN
        Ifail  = 1
        CALL CCPDPN(45,'ATOMSF','READONLY','F',0,Ifail)
        IF (Ifail.LT.0) CALL CCPERR(1,'Error opening library file')
      ELSE
        REWIND 45
      END IF
C
C---- Search for atom identifier ID
      IFAIL = -1
      LID = LENSTR(ID)
C  Big loop over possible sub-strings of ID
      DO 25 NID = LID,1,-1
        REWIND 45

C  Small loop over lines in library file
   10   CONTINUE
        READ (45,FMT=6002,END=50,ERR=40) IDIN
C
        CALL CCPUPC(IDIN)
        IF (ID2(1:NID).EQ.IDIN(1:NID)) THEN
c
c       20/11/2000 C. Vonrhein
c
c       special precautions for single character atom types:
c       skip this atom if second character is NOT '+', '-' or ' '.
c
          IF ((NID.NE.1).OR.( (NID.EQ.1).AND.(
     .                      (IDIN(2:2).EQ."+").OR.
     .                      (IDIN(2:2).EQ."-").OR.
     .                      (IDIN(2:2).EQ." ")    )
     .                    )) THEN
            Ifail = 1
            IF (NGauss.NE.2 .OR. IDIN(6:6).NE.' ') GO TO 60
          END IF
        END IF

C  Up for next line of library file
        GO TO 10
C
C---- Error reading library file
C
   40   CALL CCPERR(1,'Error reading library file')
C
C---- No match
C
   50   CONTINUE
        IF (NID.GT.1) THEN
          WRITE(STRING,'(A,A,A)') ' No match for atom ID ',ID2(1:NID),
     +     ' subtracting one character '
          CALL CCPERR(4,STRING(1:LENSTR(STRING)))
          ID2 = ID2(1:NID-1)//' '//ID2(NID+1:6)
          ID  = ID (1:NID-1)//'    '
        ENDIF
C  End of big loop. Back for smaller substring of ID.
  25  CONTINUE

      WRITE(STRING,'(A,A,A)') ' No match for atom ID ',ID2(1:1),
     +   ' giving up! '
      CALL CCPERR(4,STRING(1:LENSTR(STRING)))
      IFAIL = -1
      RETURN
C
C---- Matched atom
C
   60 READ (45,FMT=6006) IWT,IELEC,C
      READ (45,FMT=6008) A(1),A(2),A(3),A(4)
      READ (45,FMT=6008) B(1),B(2),B(3),B(4)
      READ (45,FMT=6008) CU(1),CU(2),MO(1),MO(2)
      Ifail = 0
C
C---- Format statements
C
 6002 FORMAT (A6)
 6006 FORMAT (2X,I8,2X,I8,2X,F14.6)
 6008 FORMAT (4 (2X,F14.6))
      END
C
C
C
      SUBROUTINE SFREAD(ID,NG,A,B,C,IWT,IELEC,CU,MO,IFAIL,IFILFF)
C     ===========================================================
C
C_BEGIN_SFREAD
C
C	Obsolete routine should be deleted. IFILFF not used.
C
C_END_SFREAD
C
C     .. Scalar Arguments ..
      REAL C
      INTEGER IELEC,Ifail,IWT,NG
      LOGICAL IFILFF
      CHARACTER ID*4
C     ..
C     .. Array Arguments ..
      REAL A(4),B(4),CU(2),MO(2)
C     ..
C     .. External Routines ..
      EXTERNAL SFREAD2
C     ..
      CALL SFREAD2(ID,NG,A,B,C,IWT,IELEC,CU,MO,IFAIL)
 
      RETURN
      END
C
C
C
        SUBROUTINE WBCELL(IUNIT,ARGCELL,ARGNCODE)
C       =========================================
C
C_BEGIN_WBCELL
C
C   This subroutine writes out the cell and orthogonalisation matrices, to 
C the output file. If the input parameters are null then the cell etc. are
C taken from the COMMON blocks.
C
C PARAMETERS
C
C            IUNIT (I) (INTEGER)   Channel number for output file.
C
C       ARGCELL(6) (I) (REAL)      crystallographic cell taken from COMMON
C                                  if cell = 0
C         ARGNCODE (I) (INTEGER)   NCODE number taken from COMMON if NCODE=0
C
C_END_WBCELL
C
C     .. Parameters ..
      INTEGER MAXFILESOPEN
      PARAMETER (MAXFILESOPEN=90)
C     ..
C     .. Agruments ..
      REAL ARGCELL(6)
      INTEGER ARGNCODE,IUNIT
C     ..
C     .. Variables in Common ..
      REAL CELL, RO, RF, RR
      INTEGER FILESOPEN, NCODE, TYPE, UNIT
      CHARACTER*80 LOGUNIT
      CHARACTER BRKSPGRP*15
      LOGICAL IFCRYS,IFSCAL,MATRIX,IFHDOUT,IFNEWCRYS
C     ..
C     .. Local Scalars ..
      INTEGER I, II, J, IRET
      CHARACTER*80 ERRLIN
C     ..
C     .. External Routines/Functions ..
      EXTERNAL CCPERR,RBFROR,RBRINV,MMDB_F_WBSPGRP,
     +         MMDB_F_WBCELL
C     ..
C     .. Common Blocks ..
      COMMON /RBRKAA/ FILESOPEN,LOGUNIT(MAXFILESOPEN),
     +                UNIT(MAXFILESOPEN),TYPE(MAXFILESOPEN)
      COMMON /ORTHOG/RO(4,4),RF(4,4),NCODE
      COMMON /RBRKZZ/CELL(6),RR(3,3,6),VOL,CELLAS(6)
      COMMON /RBRKXX/ IFCRYS,IFSCAL,ITYP,MATRIX,IFHDOUT,IFNEWCRYS
      COMMON /RBRKSPGRP/BRKSPGRP
C     ..
      SAVE /ORTHOG/,/RBRKAA/,/RBRKXX/,/RBRKZZ/

      II = 0
      DO 10 I=1,FILESOPEN
        IF (IUNIT .EQ. UNIT(I)) THEN
          II = I
          GOTO 20
        ENDIF
   10 CONTINUE

      ERRLIN = ' ERROR: in WBCELL file has not been opened'
      CALL CCPERR(1,ERRLIN)

   20 IF (ARGCELL(1) .EQ. 0.0) THEN
        IF (IFCRYS) CALL MMDB_F_WBCELL(IUNIT, CELL, ARGNCODE, IRET)
      ELSE
        CALL MMDB_F_WBCELL(IUNIT, ARGCELL, ARGNCODE, IRET)
      ENDIF
C...  update spacegroup information from cache
      CALL MMDB_F_WBSPGRP(IUNIT,BRKSPGRP,IRET)

      IF (ARGNCODE .NE. 0) THEN
        DO 30 I = 1,6
          CELL(I) = ARGCELL(I)
   30   CONTINUE

        CALL RBFROR
        DO 40 I = 1,3
          DO 40 J = 1,3
            RO(J,I) = RR(J,I,ARGNCODE)
   40   CONTINUE

        RO(4,4) = 1.0     
        DO 50 I=1,3
          RO(I,4) = 0.0
   50   CONTINUE

        CALL RBRINV(RO,RF)
      ENDIF

      IFNEWCRYS = .TRUE.
      RETURN
      END
C
C
C
        SUBROUTINE WREMARK(IUNIT,LINE)
C       ==============================
C
C_BEGIN_WREMARK
C
C	This subroutine writes a line to the output file. Its main use is for 
C REMARK statements in PDB.
C
C Parameters:
C
C            IUNIT (I) (CHARACTER*(*))  Channel number
C             LINE (I) (CHARACTER*(*))  line to be written, best
C                                       if declared as *80
C
C_END_WREMARK
C
C     .. Parameters ..
      INTEGER MAXFILESOPEN
      PARAMETER (MAXFILESOPEN=90)
C     ..
C     .. Arguments ..
      INTEGER IUNIT
      CHARACTER*(*) LINE
C     ..
C     .. Variables in Common ..
      INTEGER FILESOPEN,TYPE,UNIT
      CHARACTER*80 LOGUNIT
C     ..
C     .. Locals ..
      INTEGER II, IRET
      CHARACTER OUTLIN*80,ERRLIN*80
C     ..
C     .. Common Blocks ..
      COMMON /RBRKAA/ FILESOPEN,LOGUNIT(MAXFILESOPEN),
     +                UNIT(MAXFILESOPEN),TYPE(MAXFILESOPEN)
C     ..
C     .. Save Statement ..
      SAVE /RBRKAA/
C     ..
C
C---- The remark line will be truncated if it > 80 characters
C     this fits with PDB and CIF syntax.
C
      II = 0
      DO 10 I=1,FILESOPEN
        IF (IUNIT .EQ. UNIT(I)) THEN
          II = I
          GOTO 20
        ENDIF
   10 CONTINUE

      ERRLIN = ' ERROR: in WREMARK file has not been opened'
      CALL CCPERR(1,ERRLIN)

   20 CALL MMDB_F_WREMARK(IUNIT, LINE, IRET)

      RETURN
      END
C
C<FF>
C
      SUBROUTINE RWBFIN(IUN,IOUT)
C     ===========================
C
C_BEGIN_RWBFIN
C
C	This subroutine copies remaining lines straight from input to 
C output. 
C
C_END_RWBFIN
C
C     .. Scalar Arguments ..
      INTEGER IUN,IOUT
C     ..
C     .. External Routines ..
      EXTERNAL XYZADVANCE

   10 CALL XYZADVANCE(IUN,IOUT,0,*1000,*1000)
      CALL XYZADVANCE(IOUT,0,0,*1000,*1000)
      GOTO 10

 1000 RETURN
      END
C
C
C
      SUBROUTINE RWNOHEAD()
C     =====================
C
C_BEGIN_RWNOHEAD
C
C	This subroutine resets the logical variable IFHDOUT in the RWBROOK
C     common block RBRKXX, and should be called once before either
C     XYZADVANCE or WBCELL in order to prevent those routines from writing
C     headers to an output pdb file.
C     Effectively we are fooling the library that the header has already
C     been written.
C
C_END_RWNOHEAD
C
C     .. Arguments in common ..
      INTEGER ITYP
      LOGICAL IFCRYS,IFHDOUT,IFNEWCRYS,IFSCAL,MATRIX
C
C     .. Common blocks ..
      COMMON /RBRKXX/ IFCRYS,IFSCAL,ITYP,MATRIX,IFHDOUT,IFNEWCRYS
C
      IFHDOUT = .TRUE.
C
      RETURN
      END
C
C
      SUBROUTINE NUM_EXPECTED_WATERS(RESO,TEMP,FHOH,SEFHOH)
C     =====================================================
C
C_BEGIN_NUM_EXPECTED_WATERS
C
C Arguments:
C
C         RESO    (I)   REAL           Resolution in A.
C         TEMP    (I)   CHARACTER*4    'ROOM' or 'LOWT' for temperature.
C         FHOH    (O)   REAL           Prediction of N_HOH/N_at
C         SEFHOH  (O)   REAL           Standard error of FHOH
C
C     Given the resolution, this routine returns the number of water
C     molecules expected to be determined by PX (and the associated
C     standard error) as a fraction of the number of protein atoms,
C     as estimated by the statistical analysis of Carugo & Bordo, 
C     Acta Cryst D, 55, 479 (1999). Two expressions are given, one
C     for room temperature structures and one for low temperature
C     structures.
C
C_END_NUM_EXPECTED_WATERS
C
      CHARACTER*4 TEMP,TTEMP
      REAL RESO,FHOH,SEFHOH

      TTEMP = TEMP
      CALL CCPUPC(TTEMP)
      IF (TTEMP.EQ.'ROOM') THEN

        FHOH = 0.301 - 0.095*RESO
        SEFHOH = 0.092 * SQRT(0.00114 + 0.005*(RESO - 2.3)**2)

      ELSEIF (TTEMP.EQ.'LOWT') THEN

        FHOH = 0.334 - 0.110*RESO
        SEFHOH = 0.043 * SQRT(0.030 + 0.167*(RESO - 2.2)**2)

      ENDIF

      RETURN
      END
C
C

C
        SUBROUTINE WBCELLS(IUNIT,ARGCELL,ARGNCODE,NAMSPG_CIF)
C       =========================================
C
C_BEGIN_WBCELL
C
C   This subroutine writes out the cell and orthogonalisation matrices, to 
C the output file. If the input parameters are null then the cell etc. are
C taken from the COMMON blocks.
C
C PARAMETERS
C
C            IUNIT (I) (INTEGER)   Channel number for output file.
C
C       ARGCELL(6) (I) (REAL)      crystallographic cell taken from COMMON
C                                  if cell = 0
C         ARGNCODE (I) (INTEGER)   NCODE number taken from COMMON if NCODE=0
C
C_END_WBCELL
C
C     .. Parameters ..
      INTEGER MAXFILESOPEN
      PARAMETER (MAXFILESOPEN=90)
C     ..
C     .. Agruments ..
      REAL ARGCELL(6)
      INTEGER ARGNCODE,IUNIT
C     ..
C     .. Variables in Common ..
      REAL CELL, RO, RF, RR
      INTEGER FILESOPEN, NCODE, TYPE, UNIT
      CHARACTER*80 LOGUNIT
      CHARACTER BRKSPGRP*15,NAMSPG_CIF*(*)
      LOGICAL IFCRYS,IFSCAL,MATRIX,IFHDOUT,IFNEWCRYS
C     ..
C     .. Local Scalars ..
      INTEGER I, II, J, IRET
      CHARACTER*80 ERRLIN
      CHARACTER*15 SPGRP
C     ..
C     .. External Routines/Functions ..
      EXTERNAL CCPERR,RBFROR,RBRINV,MMDB_F_RBSPGRP
C     ..
C     .. Common Blocks ..
      COMMON /RBRKAA/ FILESOPEN,LOGUNIT(MAXFILESOPEN),
     +                UNIT(MAXFILESOPEN),TYPE(MAXFILESOPEN)
      COMMON /ORTHOG/RO(4,4),RF(4,4),NCODE
      COMMON /RBRKZZ/CELL(6),RR(3,3,6),VOL,CELLAS(6)
      COMMON /RBRKXX/ IFCRYS,IFSCAL,ITYP,MATRIX,IFHDOUT,IFNEWCRYS
      COMMON /RBRKSPGRP/BRKSPGRP
C     ..
      SAVE /ORTHOG/,/RBRKAA/,/RBRKXX/,/RBRKZZ/
C
      II = 0
      DO 10 I=1,FILESOPEN
        IF (IUNIT .EQ. UNIT(I)) THEN
          II = I
          GOTO 20
        ENDIF
   10 CONTINUE

      ERRLIN = ' ERROR: in WBCELL file has not been opened'
      CALL CCPERR(1,ERRLIN)
 
   20 CALL MMDB_F_RBSPGRP(IUNIT,SPGRP,IRET)
C    Has the spacegroup name been set yet?
      ILEN1 = LENSTR(SPGRP)
      ILEN2 = LENSTR(NAMSPG_CIF)
 
       IF(SPGRP.NE.' '.AND.
     +    SPGRP(1:ILEN1).NE.NAMSPG_CIF(1:ILEN2)) then
         WRITE(6,'(A,/,4(A,5x))') '*** Incompatible space group names:',
     +   ' From PDB input:         ',SPGRP,
     +   ' To be written to output:',NAMSPG_CIF
         CALL CCPERR(2,' Incompatible space group names')
      END IF
      BRKSPGRP = NAMSPG_CIF(1:ILEN2)

      CALL WBCELL(IUNIT,ARGCELL,ARGNCODE)
      RETURN
      END
C
