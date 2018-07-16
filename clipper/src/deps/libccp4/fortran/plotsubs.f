C
C     plotsubs.f
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
      SUBROUTINE PLTAXS(VWPORT,XLOW,XHIGH,NXTICK,NXSTCK,LXTICK,XLABEL,
     .                         YLOW,YHIGH,NYTICK,NYSTCK,LYTICK,YLABEL)
C     ================================================================
C
C
C  Draw box of axes, and set up user coordinate transformation
C
C
C  VWPORT(2,2)  viewport, ie coordinates (in range 0 to 1) of
C               bottom left corner (i,1) and top right corner (i,2)
C               of box
C  XLOW, XHIGH  range of values on x
C  NXTICK       if .gt.0 approximate number of tick marks on x (=0 if none)
C               If .lt.0, use exactly the range (XLOW - XHIGH) & number
C               of ticks specified
C  NXSTCK       number of subsidiary tick marks / division on x
C  LXTICK       |LXTICK| = 1 ticks inside box
C                        = 2 ticks outside box
C               .gt. 0     draw tick & label at end of axis
C               .lt. 0     omit last tick & label at end of axis
C  XLABEL       label for x axis
C
C    similarly,  YLOW,YHIGH,NYTICK,NYSTCK,LYTICK,YLABEL
C
C
      REAL VWPORT(2,2),XLOW,XHIGH,YLOW,YHIGH
      INTEGER NXTICK,NXSTCK,LXTICK,NYTICK,NYSTCK,LYTICK
      CHARACTER*(*) XLABEL,YLABEL
C
      REAL X,Y,TSIZE,WINDOW(2,2),WIN1,WIN2
      INTEGER L,NLAB,IFLAG,JTICK,LENSTR
      EXTERNAL LINSUB,TXTSUB
C
C
      TSIZE=1.2
C
C x axis
      IF(NXTICK.GT.0) THEN
C   Find "best" tick positions
            CALL UGLNDX(XLOW,XHIGH,NXTICK-2,NXTICK+2,
     .          WIN1,WIN2,NLAB)
C           avoid aliasing array elements
            WINDOW (1,1) = WIN1
            WINDOW (1,2) = WIN2
      ELSE
C   Use given values
            NLAB=IABS(NXTICK)
            WINDOW(1,1)=XLOW
            WINDOW(1,2)=XHIGH
      ENDIF
C   Get best format for numbers
      CALL PLTNDG(WINDOW(1,1),WINDOW(1,2),NLAB)
C   bottom axis
      IFLAG=1
      JTICK=LXTICK
      CALL UGLNAX(LINSUB,TXTSUB,IFLAG,
     .   VWPORT(1,1),VWPORT(2,1),VWPORT(1,2),VWPORT(2,1),
     .   WINDOW(1,1),WINDOW(1,2),NLAB,NXSTCK,JTICK)
C   top
      IFLAG=0
      JTICK=(3-IABS(LXTICK))*ISIGN(1,LXTICK)
      CALL UGLNAX(LINSUB,TXTSUB,IFLAG,
     .   VWPORT(1,1),VWPORT(2,2),VWPORT(1,2),VWPORT(2,2),
     .   WINDOW(1,1),WINDOW(1,2),NLAB,NXSTCK,JTICK)
C
C  Plot label
      L=LENSTR(XLABEL)
      X=0.5*(VWPORT(1,1)+VWPORT(1,2))
      Y=VWPORT(2,1)-0.07
      CALL PLTCTX(XLABEL(1:L),X,Y,TSIZE,2,0.0)
C
C y axis
      IF(NYTICK.GT.0) THEN
C   Find "best" tick positions
            CALL UGLNDX(YLOW,YHIGH,NYTICK-2,NYTICK+2,
     .              WIN1,WIN2,NLAB)
C           avoid aliasing array elements
            WINDOW (2,1) = WIN1
            WINDOW (2,2) = WIN2
      ELSE
C   Use given values
            NLAB=IABS(NYTICK)
            WINDOW(2,1)=YLOW
            WINDOW(2,2)=YHIGH
      ENDIF
C   Get best format for numbers
      CALL PLTNDG(WINDOW(2,1),WINDOW(2,2),NLAB)
C   left
      IFLAG=2
      JTICK=(3-IABS(LYTICK))*ISIGN(1,LYTICK)
      CALL UGLNAX(LINSUB,TXTSUB,IFLAG,
     .   VWPORT(1,1),VWPORT(2,1),VWPORT(1,1),VWPORT(2,2),
     .   WINDOW(2,1),WINDOW(2,2),NLAB,NYSTCK,JTICK)
C   right
      IFLAG=0
      JTICK=LYTICK
      CALL UGLNAX(LINSUB,TXTSUB,IFLAG,
     .   VWPORT(1,2),VWPORT(2,1),VWPORT(1,2),VWPORT(2,2),
     .   WINDOW(2,1),WINDOW(2,2),NLAB,NYSTCK,JTICK)
C
C  Plot label
      L=LENSTR(YLABEL)
      X=VWPORT(1,1)-0.11
      Y=0.5*(VWPORT(2,1)+VWPORT(2,2))
      CALL PLTCTX(YLABEL(1:L),X,Y,TSIZE,3,90.0)
C
C Set up user transformation
      CALL PLTWIN(VWPORT,WINDOW)
C
      END
      SUBROUTINE PLTCHS(SCALE)
C     ========================
C
C Set character scale as fraction of space 0-1
C
      REAL SCALE
C
      CALL GSSCLC(SCALE,SCALE)
      END
      SUBROUTINE PLTCTX(STRING,X,Y,SIZE,MODE,ANGLE)
C     =============================================
C
C*** Plot84 version
C  Plot STRING at X,Y, rotated by ANGLE, relative size SIZE
C    MODE = 1, 2, 3 justify left, middle, right
C
      CHARACTER*(*) STRING
      REAL X,Y,ANGLE,DTOR,SIZE
      INTEGER MODE
C
C
      DTOR=ATAN(1.)/45.
      CALL GSANCU(X,Y)
      CALL GSCROT(ANGLE*DTOR,(ANGLE+90.)*DTOR)
      CALL GSCETS(STRING,SIZE,SIZE,MODE)
      END
      SUBROUTINE PLTDBU(X,Y)
C     ======================
C
C Draw by point X,Y (user units)
C
      REAL X,Y,XV,YV
C
C  transform coordinates
      CALL PLTTNS(X,Y,XV,YV)
      CALL PLTDBY(XV,YV)
      END
      SUBROUTINE PLTDBY(X,Y)
C     ======================
C
C*** Plot84 version
C   Draw to X,Y
C
      REAL X,Y
      CALL GSDWBY(X,Y)
      END
      SUBROUTINE PLTDRW(X,Y)
C     ======================
C
C*** Plot84 version
C   Draw to X,Y
C
      REAL X,Y
      CALL GSDWTO(X,Y)
      END
      SUBROUTINE PLTDSH(M,R,D,T)
C     ==========================
C
C  Set dash parameters
C
C   M =0 SOLID LINE
C   M=1 DASHED LINE (DOT DUMMY PARAMETER)
C   M=2  CHAINED LINE
C   R == REPEAT repeat length
C   D == DASH  dash length
C   T == DOT   not used
C  REPEAT DASH DOT  are in basic picture units (ie range 0 to 1 )
C
      INTEGER M
      REAL R,D,T
C
      COMMON /DASH/ LPT,MODE,REPEAT,DASH,DOT,A1,B1
      REAL REPEAT,DASH,DOT,A1,B1
      INTEGER LPT,MODE
      SAVE /DASH/
C
C
      MODE=M
      REPEAT=R
      DASH=D
      DOT=T
      LPT=0
      END
      SUBROUTINE PLTDWU(X,Y)
C     ======================
C
C Draw to point X,Y (user units)
C
      REAL X,Y,XV,YV
C
C  transform coordinates
      CALL PLTTNF(X,Y,XV,YV)
      CALL PLTDRW(XV,YV)
      END
      SUBROUTINE PLTFNM(FNUM,NDIG,NAFTER,X,Y,SIZE,MODE,ANGLE)
C     ========================================================
C
C*** Plot84 version
C  Plot real FNUM at X,Y, rotated by ANGLE
C   NDIG digits, NAFTER digits after decimal point
C    MODE = 1, 2, 3 justify left, middle, right
C
      REAL FNUM,X,Y,ANGLE,DTOR,SIZE
      INTEGER NDIG,NAFTER,MODE
C
C
      DTOR=ATAN(1.)/45.
      CALL GSANCU(X,Y)
      CALL GSCROT(ANGLE*DTOR,(ANGLE+90.)*DTOR)
      CALL GSFNUM(FNUM,NDIG,NAFTER,SIZE,SIZE,MODE)
      END
      SUBROUTINE PLTINI
C     =================
C
C*** Plot84 version
C  Initialize plot on logical name PLOT
C
      REAL XSIZE,YSIZE,CSIZE
      DATA XSIZE,YSIZE,CSIZE/200.,200.,0.015/
C
      CALL GSINIT('PLOT')
      CALL GSBSIZ(XSIZE,YSIZE)
C
10    CALL GSPICT
C  User coordinates in range 0 - 1 on x & y
      CALL GSORGD(0.01*XSIZE,0.01*YSIZE)
      CALL GSSCLU(XSIZE*0.99,YSIZE*0.99)
C  Scale characters to CSIZE (user units)
      CALL GSTLNK(1)
      CALL GSSCLC(CSIZE,CSIZE)
C  Centred characters
      CALL GSCENC(1)
      RETURN
C
      ENTRY PLTPIC
C       =============
C
C Start new picture
C
      CALL GSENDP
      GO TO 10
      END
      SUBROUTINE PLTINM(INUM,NDIG,X,Y,SIZE,MODE,ANGLE)
C     ============================================
C
C*** Plot84 version
C  Plot integer INUM at X,Y, rotated by ANGLE
C   NDIG digits
C    MODE = 1, 2, 3 justify left, middle, right
C
      REAL X,Y,ANGLE,DTOR,SIZE
      INTEGER INUM,NDIG,MODE
C
C
      DTOR=ATAN(1.)/45.
      CALL GSANCU(X,Y)
      CALL GSCROT(ANGLE*DTOR,(ANGLE+90.)*DTOR)
      CALL GSINUM(INUM,NDIG,SIZE,SIZE,MODE)
      END
      SUBROUTINE PLTLIN(X,Y,NPT)
C     =========================
C
C Draw lines to join NPT points X,Y
C Coordinates are in user units
C
      INTEGER NPT
      REAL X(NPT),Y(NPT)
C
      COMMON /DASH/ LPT,MODE,REPEAT,DASH,DOT,A1,B1
      REAL REPEAT,DASH,DOT,A1,B1
      INTEGER LPT,MODE
      SAVE /DASH/
C
      REAL SMALL,GAP,A,B,X1,Y1,X2,Y2,DX,DY,R,XD,YD,XG,YG,D
      INTEGER I,J
      DATA SMALL/0.001/
C
C MODE =0 FOR SOLID LINES
      IF(NPT.LE.1) RETURN
      IF(MODE.GT.0) GO TO 10
C
C SOLID
      CALL PLTMVU(X(1),Y(1))
      DO 1 I=2,NPT
1     CALL PLTDWU(X(I),Y(I))
      RETURN
C
C
C DASHED LINES, REPEAT IS REPEAT LENGTH, DASH IS DASH LENGTH
C DOT IS DUMMY (NO CHAIN LINES)
10    GAP=REPEAT-DASH
      A=DASH
      B=GAP
C LPT=0 IF 1ST CALL SINCE CALL TO S/R DASHED, OTHERWISE KEEP DASHES IN PHASE
      IF(LPT.EQ.0) GO TO 11
      A=A1
      B=B1
C
11    LPT=1
C  Convert user units to picture units before working out dashing
      CALL PLTTNF(X(1),Y(1),X1,Y1)
CC      X1=X(1)
CC      Y1=Y(1)
      I=2
      CALL PLTMOV(X1,Y1)
C
C COME HERE FOR NEW LINE
15    CALL PLTTNF(X(I),Y(I),X2,Y2)
CC15    X2=X(I)
CC      Y2=Y(I)
      DX=X2-X1
      DY=Y2-Y1
      D=SQRT(DX*DX+DY*DY)
      IF(D.LT.SMALL) GO TO 30
      DX=DX/D
      DY=DY/D
      R=D
      A1=0.
      B1=0.
      J=-1
C
C R IS REMAINING LINE LENGTH
C COME HERE FOR EACH DASH REPEAT
20    IF(R.LT.(A+B)) GO TO 25
      IF(J) 26,26,27
C
C LAST BIT, GET REMAINING DASH LENGTH
25    IF(R.LT.A) GO TO 21
C LAST PART IS IN GAP
      A1=0.
      B1=B-R+A
      B=R-A
      GO TO 22
C LAST PART IS IN DASH
21    A1=A-R
      B1=GAP
      A=R
      B=0.
22    J=-2
C
26    J=J+1
      XD=A*DX
      YD=A*DY
      XG=B*DX
      YG=B*DY
C
27    IF(A.GT.SMALL) CALL PLTDBY(XD,YD)
      IF(B.GT.SMALL) CALL PLTMBY(XG,YG)
      R=R-A-B
      IF(J) 30,28,20
C RESET DASH LENGTH AFTER 1ST REPEAT
28    A=DASH
      B=GAP
      J=0
      GO TO 20
C
C END OF LINE, RESTORE UNUSED DASH LENGTH
30    I=I+1
      IF(I.GT.NPT) RETURN
      X1=X2
      Y1=Y2
      A=A1
      B=B1
      GO TO 15
C
      END
      SUBROUTINE PLTMBU(X,Y)
C     ======================
C
C Move by point X,Y (user units)
C
      REAL X,Y,XV,YV
C
C  transform coordinates
      CALL PLTTNS(X,Y,XV,YV)
      CALL PLTMBY(XV,YV)
      END
      SUBROUTINE PLTMBY(X,Y)
C     ======================
C
C*** Plot84 version
C   Move by X,Y
C
      REAL X,Y
      CALL GSMVBY(X,Y)
      END
      SUBROUTINE PLTMOV(X,Y)
C     ======================
C
C*** Plot84 version
C   Move to X,Y
C
      REAL X,Y
      CALL GSMVTO(X,Y)
      END
      SUBROUTINE PLTMVU(X,Y)
C     ======================
C
C Move to point X,Y (user units)
C
      REAL X,Y,XV,YV
C
C  transform coordinates
      CALL PLTTNF(X,Y,XV,YV)
      CALL PLTMOV(XV,YV)
      END
      SUBROUTINE PLTNDG(LOW,HIGH,NLAB)
C     ================================
C
C Try to determine "best" number format for axis labelling
C   Input:
C       low, high       low & high values for tick marks
C       nlab            number of tick marks
C
C   Output in COMMON /PLTDIG/
C       ndec    number of digits in number
C       idec    number of digits after decimal point (0 for integer)
C
      REAL RANGE,LOW,HIGH,DVAL,BIG
      INTEGER NLAB,N
      LOGICAL LINT
C
      COMMON /PLTDIG/ NDEC,IDEC
      INTEGER NDEC,IDEC
      SAVE /PLTDIG/
C
      RANGE=HIGH-LOW
      BIG=MAX(ABS(LOW),ABS(HIGH))
C
      NDEC=5
      IDEC=1
      LINT=.FALSE.
C
      DVAL=1.0
      IF(NLAB.GT.1) DVAL=RANGE/(NLAB-1)
      IDEC=1-ALOG10(ABS(DVAL))
C
      IF(ABS(RANGE).GT.9.99) THEN
            LINT=.TRUE.
            NDEC=2+ALOG10(BIG)
      ENDIF
C
      IF(AMOD(DVAL,1.0).EQ.0.0 .AND.
     .     AMOD(LOW,1.0).EQ.0.0) LINT=.TRUE.
C
      IF(LINT) IDEC=0
      IF((LOW.LT.0.0.OR.HIGH.LT.0.0)) THEN
            N=0
            IF(IDEC.GT.0) N=2+IDEC
            IF(NDEC.LT.N) NDEC=NDEC+1
      ENDIF
C
      END
      SUBROUTINE PLTSTP
C     =================
C
C*** Plot84 version
C  Stop plotting
C
      CALL GSENDP
      CALL GSSTOP
      END
      SUBROUTINE PLTTNF(X,Y,XV,YV)
C     ============================
C
C Transform user X,Y to plotter XV,YV (range 0 to 1)
C
C
      REAL X,Y,XV,YV
C
      COMMON /PLTTRN/ ORIGX,ORIGY,VSCALX,VSCALY
      REAL ORIGX,ORIGY,VSCALX,VSCALY
      SAVE /PLTTRN/
C
      XV=X*VSCALX+ORIGX
      YV=Y*VSCALY+ORIGY
      END
C
      SUBROUTINE PLTTNS(X,Y,XV,YV)
C     ============================
C
C Transform user X,Y to plotter XV,YV (range 0 to 1)
C    without translation (scaling only)
C
C
      REAL X,Y,XV,YV
C
      COMMON /PLTTRN/ ORIGX,ORIGY,VSCALX,VSCALY
      SAVE /PLTTRN/
      REAL ORIGX,ORIGY,VSCALX,VSCALY
C
      XV=X*VSCALX
      YV=Y*VSCALY
      RETURN
      END
C
      SUBROUTINE PLTWIN(VWPORT,WINDOW)
C     ================================
C
C Set up mapping of viewport VWPORT to user window WINDOW
C
      REAL VWPORT(2,2),WINDOW(2,2)
C
      COMMON /PLTTRN/ ORIGX,ORIGY,VSCALX,VSCALY
      REAL ORIGX,ORIGY,VSCALX,VSCALY
      SAVE /PLTTRN/
C
      VSCALX=(VWPORT(1,2)-VWPORT(1,1))/(WINDOW(1,2)-WINDOW(1,1))
      VSCALY=(VWPORT(2,2)-VWPORT(2,1))/(WINDOW(2,2)-WINDOW(2,1))
      ORIGX=VWPORT(1,1)-WINDOW(1,1)*VSCALX
      ORIGY=VWPORT(2,1)-WINDOW(2,1)*VSCALY
      END
      SUBROUTINE UGLNAX(LSUB,TSUB,TFLG,XCLO,YCLO,XCHI,YCHI,
     X                  LOLB,HILB,NLAB,NSUBTK,LLTICK)
C
C *******************  THE UNIFIED GRAPHICS SYSTEM  *******************
C *                    AXIS GENERATION SUBROUTINE                     *
C *                                                                   *
C *  THIS SUBROUTINE MAY BE USED TO GENERATE A DESCRIPTION OF AN      *
C *  AXIS WITH LINEAR LABELING.  A PAIR OF USER SUPPLIED SUBROUTINES  *
C *  ARE CALLED TO PROCESS THE LINE SEGMENT END POINTS AND THE TEXT.  *
C *                                                                   *
C *  THE CALLING SEQUENCE IS:                                         *
C *    CALL UGLNAX(LSUB,TSUB,TFLG,XCLO,YCLO,XCHI,YCHI,                *
C *                LOLB,HILB,NLAB,NSUBTK,LLTICK)                      *
C *                                                                   *
C *  THE PARAMETERS IN THE CALLING SEQUENCE ARE:                      *
C *    LSUB  THE LINE SEGMENT END POINT SUBROUTINE.                   *
C *    TSUB  THE LABEL SUBROUTINE.                                    *
C *    TFLG  A FLAG THAT IS PASSED TO TSUB.                           *
C *    XCLO  X COORDINATE OF THE LOW END OF THE AXIS.                 *
C *    YCLO  Y COORDINATE OF THE LOW END OF THE AXIS.                 *
C *    XCHI  X COORDINATE OF THE HIGH END OF THE AXIS.                *
C *    YCHI  Y COORDINATE OF THE HIGH END OF THE AXIS.                *
C *    LOLB  DATA VALUE AT THE LOW END OF THE AXIS.                   *
C *    HILB  DATA VALUE AT THE HIGH END OF THE AXIS.                  *
C *    NLAB  NUMBER OF LABELS ON THE AXIS. if =0, no tick marks       *
C *    NSUBTK Number of subsidiary tick marks between main ones.      *
C *    LLTICK =1 tick on left, =2 tick on right =3 both               *
C *          if .lt.0,  don't plot tick or label at end of axis       *
C *                                                                   *
C *                          ROBERT C. BEACH                          *
C *                    COMPUTATION RESEARCH GROUP                     *
C *                STANFORD LINEAR ACCELERATOR CENTER                 *
C *                                                                   *
C *********************************************************************
C  PRE 28/7/88 bodged
C
      EXTERNAL      LSUB,TSUB
      INTEGER       TFLG
      REAL          XCLO,YCLO,XCHI,YCHI,LOLB,HILB
      INTEGER       NLAB,NSUBTK,LLTICK
C
C
      INTEGER     EXNT
      REAL        EXLS,EXRS
C
      REAL          DLAX,DLAY
      REAL          DLTX,DLTY
      REAL          DLBX,DLBY,DLBL
      REAL          DLCX,DLCY
      REAL          XCRD,YCRD,VALU
      REAL          XCDS,YCDS
C
      REAL          FLT1
      INTEGER       INT1,INT2,LTICK,J
C
      LTICK=IABS(LLTICK)
C
      EXNT=MAX(0,NSUBTK)
      EXLS=0.0
      EXRS=EXLS
      FLT1=0.01*MAX(ABS(XCHI-XCLO),ABS(YCHI-YCLO))
      IF(MOD(LTICK,2).EQ.1) EXLS=FLT1
      IF(LTICK.GE.2) EXRS=FLT1
C
C  DRAW THE AXIS STARTING AT THE LOW END.
      IF (NLAB.LT.2) THEN
C  Just draw line if .le.2 ticks
            CALL LSUB(XCLO,YCLO,0)
            CALL LSUB(XCHI,YCHI,1)
            GO TO 201
      ENDIF
C
      CALL LSUB(XCLO,YCLO,0)
      CALL LSUB(XCHI,YCHI,1)
C
C  INITIALIZE THE TIC MARK AND LABEL LOOP.
      DLAX=XCHI-XCLO
      DLAY=YCHI-YCLO
      FLT1=SQRT(DLAX*DLAX+DLAY*DLAY)
      DLTX=DLAY/FLT1
      DLTY=-DLAX/FLT1
      FLT1=REAL(NLAB-1)
      DLBX=DLAX/FLT1
      DLBY=DLAY/FLT1
      DLBL=(HILB-LOLB)/FLT1
      IF (EXNT.GE.1) THEN
        FLT1=REAL(EXNT+1)
        DLCX=DLBX/FLT1
        DLCY=DLBY/FLT1
      END IF
C
C  LOOP TO GENERATE THE TIC MARKS AND LABELS.
C  Skip first tick (end of axis) if LLTICK .lt. 0
      J=1
      IF(LLTICK.LT.0) J=2
      DO 102 INT1=J,NLAB
C  GENERATE TIC MARK POSITION.
        FLT1=REAL(INT1-1)
        XCRD=XCHI-FLT1*DLBX
        YCRD=YCHI-FLT1*DLBY
C  GENERATE A LABEL.
        VALU=HILB-FLT1*DLBL
        CALL TSUB(XCRD,YCRD,VALU,TFLG)
C  GENERATE A PRIMARY TIC MARK.
        CALL LSUB(XCRD+EXRS*DLTX,YCRD+EXRS*DLTY,0)
        CALL LSUB(XCRD-EXLS*DLTX,YCRD-EXLS*DLTY,1)
C  GENERATE SECONDARY TIC MARKS.
        IF ((INT1.NE.NLAB).AND.(EXNT.GE.1)) THEN
          DO 101 INT2=1,EXNT
            FLT1=REAL(INT2)
            XCDS=XCRD-FLT1*DLCX
            YCDS=YCRD-FLT1*DLCY
            CALL LSUB(XCDS+0.5*EXRS*DLTX,YCDS+0.5*EXRS*DLTY,2)
            CALL LSUB(XCDS-0.5*EXLS*DLTX,YCDS-0.5*EXLS*DLTY,3)
  101     CONTINUE
        END IF
  102 CONTINUE
C
  201 RETURN
C
C
      END
      SUBROUTINE UGLNDX(LODA,HIDA,MINL,MAXL,LOLB,HILB,NLAB)
C
C *******************  THE UNIFIED GRAPHICS SYSTEM  *******************
C *               AUXILIARY AXIS GENERATION PROGRAM                   *
C *                                                                   *
C *  THIS SUBROUTINE IS AN AID IN USING SUBROUTINE UGLNAX.  IF A      *
C *  PROGRAMMER DETERMINES THE EXTENT OF THE DATA AT EXECUTION TIME,  *
C *  IT CAN BE A PROBLEM TO SUPPLY VALUES OF THE ARGUMENTS LOLB,      *
C *  HILB, AND NLAB WHICH RESULT IN "ROUNDED NUMBERS" FOR THE         *
C *  LABELS.  THIS SUBROUTINE ACCEPTS AS ITS INPUT THE EXTENT OF THE  *
C *  DATA AND THE APPROXIMATE NUMBER OF LABELS TO BE USED.  ITS       *
C *  OUTPUT IS THE THREE ARGUMENTS FOR UGLNAX.                        *
C *                                                                   *
C *  THE CALLING SEQUENCE IS:                                         *
C *    CALL UGLNDX(LODA,HIDA,MINL,MAXL,LOLB,HILB,NLAB)                *
C *                                                                   *
C *  THE PARAMETERS IN THE CALLING SEQUENCE ARE:                      *
C *    LODA  LOW LIMIT OF ACTUAL DATA.                                *
C *    HIDA  HIGH LIMIT OF ACTUAL DATA.                               *
C *    MINL  MINIMUM NUMBER OF LABELS.                                *
C *    MAXL  MAXIMUM NUMBER OF LABELS.                                *
C *    LOLB  COMPUTED LOW LIMIT OF AXIS.                              *
C *    HILB  COMPUTED HIGH LIMIT OF AXIS.                             *
C *    NLAB  COMPUTED NUMBER OF LABELS.                               *
C *                                                                   *
C *                          ROBERT C. BEACH                          *
C *                    COMPUTATION RESEARCH GROUP                     *
C *                STANFORD LINEAR ACCELERATOR CENTER                 *
C *                                                                   *
C *********************************************************************
C
      REAL          LODA,HIDA
      INTEGER       MINL,MAXL
      REAL          LOLB,HILB
      INTEGER       NLAB
C
CCCC      INCLUDE       'UGSYSTEM:UGERRCBK.FOR/LIST'
C
      REAL          PVAL(3)
      INTEGER       NPVL
C
      REAL          DALO,DAHI,SCAL,DLOX,DHIX,SCLX,OVLP
      REAL          TVR1,TVR2
      INTEGER       IMIN,IMAX
C
      INTEGER       INT1,INT2,INT3,INT4
C
      DATA          PVAL/1.0,2.0,5.0/
      DATA          NPVL/3/
C
C  ADJUST THE GIVEN DATA LIMITS.
      IF (LODA.LT.HIDA) THEN
        DALO=LODA
        DAHI=HIDA
      ELSE
        DALO=HIDA-0.5
        DAHI=LODA+0.5
      END IF
      TVR1=0.0005*(DAHI-DALO)
      DALO=DALO+TVR1
      DAHI=DAHI-TVR1
C  INITIALIZE THE LOOP TO FIND THE BEST VALUES.
      OVLP=1E20
      IMIN=MAX(2,MINL)
      IMAX=MAX(IMIN,MAXL)
C  LOOP TO FIND THE BEST VALUES.
      DO 104 INT1=IMIN,IMAX
        SCAL=INT1-1
        TVR2=(DAHI-DALO)/SCAL
        TVR1=LOG10(TVR2)
        INT2=INT(TVR1)
        IF (TVR1.LT.0.0) INT2=INT2-1
        TVR2=TVR2/(10.0**INT2)
        IF (TVR2.GT.PVAL(NPVL)) THEN
          INT2=INT2+1
          TVR2=TVR2/10.0
        END IF
        DO 101 INT4=1,NPVL
          INT3=INT4
          IF (PVAL(INT3).GE.TVR2) GO TO 102
  101   CONTINUE
  102   SCLX=PVAL(INT3)*(10.0**INT2)
        TVR1=(DAHI+DALO-SCLX*SCAL)/(2.0*SCLX)
        DLOX=AINT(TVR1)*SCLX
        IF (TVR1.LT.0.0) DLOX=DLOX-SCLX
  103   DHIX=DLOX+SCLX*SCAL
        IF (DAHI.GT.DHIX) THEN
          DLOX=DLOX+SCLX
          IF ((DLOX-DALO).LE.(0.005*SCLX)) GO TO 103
          IF (NPVL.NE.INT3) THEN
            INT3=INT3+1
          ELSE
            INT3=1
            INT2=INT2+1
          END IF
          GO TO 102
        END IF
        IF ((DHIX-DLOX).LT.OVLP) THEN
          OVLP=DHIX-DLOX
          LOLB=DLOX
          HILB=DHIX
          NLAB=INT1
        END IF
  104 CONTINUE
C
      RETURN
C
      END
      SUBROUTINE LINSUB(X,Y,IFLAG)
C     ============================
C
C  Line drawing subroutine called by UGLNAX
C   Move (IFLAG = 0 or 2) or draw (IFLAG = 1 or 3) to X,Y
C
      REAL X,Y
      INTEGER IFLAG
C
      IF(MOD(IFLAG,2).EQ.0) THEN
            CALL PLTMOV(X,Y)
      ELSE
            CALL PLTDRW(X,Y)
      ENDIF
      END
      SUBROUTINE TXTSUB(X,Y,VALU,IFLAG)
C     =================================
C
C Text drawing subroutine called by UGLNAX
C
C   IFLAG = 1 horizontal axis
C         = 2 vertical axis
C
      REAL X,Y,VALU,XX,YY,ANGLE,SIZE
      INTEGER MODE,IFLAG
C
      COMMON /PLTDIG/ NDEC,IDEC
      INTEGER NDEC,IDEC
      SAVE /PLTDIG/
C
      SIZE=1.0
      XX=X
      YY=Y
      IF(IFLAG.EQ.1) THEN
            YY=YY-0.03
            MODE=2
            ANGLE=0.0
      ELSEIF(IFLAG.EQ.2) THEN
            XX=XX-0.02
            MODE=3
            ANGLE=0.0
      ELSE
            RETURN
      ENDIF
C
      IF(IDEC.EQ.0) THEN
            CALL PLTINM(NINT(VALU),NDEC,XX,YY,SIZE,MODE,ANGLE)
      ELSE
            CALL PLTFNM(VALU,NDEC,IDEC,XX,YY,SIZE,MODE,ANGLE)
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE PL84CLEAR
C     =====================
C
C
C---- Start a new picture
C
C
C     .. External Subroutines ..
      EXTERNAL GSPICT
C     ..
      CALL GSPICT
      END
C
C
C
      SUBROUTINE PL84CLOSE
C     =====================
C
C
C     .. External Subroutines ..
      EXTERNAL GSSTOP
C     ..
      CALL GSSTOP
      END
C
C
C
      SUBROUTINE PL84DRAW(IX,IY)
C     ===========================
C
C
C---- Convert from 10 micron units to mm
C
C
C     .. Scalar Arguments ..
      INTEGER IX,IY
C     ..
C     .. Local Scalars ..
      REAL X,Y
C     ..
C     .. External Subroutines ..
      EXTERNAL GSDWTO
C     ..
      X = 0.01*IX
      Y = 0.01*IY
      CALL GSDWTO(X,Y)
      END
C
C
C
      SUBROUTINE PL84END
C     ===================
C
C
C     .. External Subroutines ..
      EXTERNAL GSENDP
C     ..
      CALL GSENDP
      END
C
C
C
C    ****  PLOT84 ROUTINES *****
C          ===============
C
C
      SUBROUTINE PL84INIT
C     ====================
C
C
C---- Initialise PLOT84 file
C
C
C     .. External Subroutines ..
      EXTERNAL GSINIT
C     ..
      CALL GSINIT('PLOT')
      END
C
C
C
      SUBROUTINE PL84INTEG(INUM,NDIG,SIZX,SIZY,NJUST)
C     ===============================================
C
C
C     .. Scalar Arguments ..
      REAL SIZX,SIZY
      INTEGER INUM,NDIG,NJUST
C     ..
C     .. External Subroutines ..
      EXTERNAL GSINUM
C     ..
      CALL GSINUM(INUM,NDIG,SIZX,SIZY,NJUST)
      END
C
C
C
      SUBROUTINE PL84MOVE(IX,IY)
C     ===========================
C
C
C
C---- Convert from 10 micron units to mm
C
C
C     .. Scalar Arguments ..
      INTEGER IX,IY
C     ..
C     .. Local Scalars ..
      REAL X,Y
C     ..
C     .. External Subroutines ..
      EXTERNAL GSMVTO
C     ..
      X = 0.01*IX
      Y = 0.01*IY
      CALL GSMVTO(X,Y)
      END
C
C
C
      SUBROUTINE PL84NEWLINE(X)
C     ==========================
C
C
C     .. Scalar Arguments ..
      REAL X
C     ..
C     .. External Subroutines ..
      EXTERNAL GSLNFD
C     ..
      CALL GSLNFD(X)
      END
C
C
C
      SUBROUTINE PL84ORIGIN(X,Y)
C     ===========================
C
C
C     .. Scalar Arguments ..
      REAL X,Y
C     ..
C     .. External Subroutines ..
      EXTERNAL GSANCD
C     ..
      CALL GSANCD(X,Y)
      END
C
C
C
      SUBROUTINE PL84REAL(X,NDIG,NAFTER,SIZX,SIZY,NJUST)
C     ===================================================
C
C
C     .. Scalar Arguments ..
      REAL SIZX,SIZY,X
      INTEGER NAFTER,NDIG,NJUST
C     ..
C     .. External Subroutines ..
      EXTERNAL GSFNUM
C     ..
      CALL GSFNUM(X,NDIG,NAFTER,SIZX,SIZY,NJUST)
      END
C
C
C
      SUBROUTINE PL84STRING(STR)
C     ===========================
C
C
C
C
C     .. Scalar Arguments ..
      CHARACTER STR* (*)
C     ..
C     .. External Subroutines ..
      EXTERNAL GSSTRC
C     ..
      CALL GSSTRC(STR)
      END
C
C
C
      SUBROUTINE PL84XCUR(X)
C     =======================
C
C
C     .. Scalar Arguments ..
      REAL X
C     ..
C     .. External Subroutines ..
      EXTERNAL GSSCUR
C     ..
      CALL GSSCUR(X)
      END
