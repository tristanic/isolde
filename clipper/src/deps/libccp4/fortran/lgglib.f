C
C     lgglib.f: functions used by Guoguang Lu's programs
C     Copyright (C) 1999  Guoguang Lu
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
      function angle(v1,v2)
c calculate the angle between two vectors v1, and v2
c v1, v2 are the input 3*vectors
c angle is the output angle between them (in degrees)
      real*4 v1(3),v2(3)
c
      external acosd
c	write(7,*) 'v1,v2',v1,v2
c	write(7,*) 'vem v1,v2',vem(3,v1),vem(3,v2)
c	write(7,*) 'poimut',poimult(3,3,v1,v2)
      arc = poimult(3,3,v1,v2)/(vem(3,v1)*vem(3,v2))
c	write(7,*) 'arc=',arc
      if (abs(arc).gt.1.) then
       if (abs(arc)-1.gt.1e-5) write(7,*)'Warning arccosd > 1'
       if (arc.gt.1.) arc = 1.
       if (arc.lt.-1.) arc = -1.
      end if
      angle = acosd(arc)
      end
c
c
      SUBROUTINE ANGLZ(CELL,V,A)
      DIMENSION CELL(6),V(3),A(3)
      DO 10 I=1,3
10	A(I)=CELL(I)*V(I)
      END
c
c
      SUBROUTINE ANTIARR(IM,IN,A1,A)
c    Transpose an array
      DIMENSION A1(IM,IN),A(IN,IM)
      DO 10 II=1,IM
      DO 10 IJ=1,IN
10	A(IJ,II)=A1(II,IJ)
      END
c
c
      SUBROUTINE ARRAD(IM,IN,A1,A2,OUT)
c    Add two arrays      
      DIMENSION A1(IM,IN),A2(IM,IN),OUT(IM,IN)
      DO 10 II=1,IM
      DO 10 IJ=1,IN
10	OUT(II,IJ)=A1(II,IJ)+A2(II,IJ)
      END
c
c
      SUBROUTINE ARRGIVE(N,XYZ0,XYZ)
C    Copy an N-dimensional real array XYZ0 to another one XYZ.
      REAL*4 XYZ(N),XYZ0(N)
      DO I = 1, N
       XYZ(I) = XYZ0(I)
      END DO
      END
c
c
      SUBROUTINE ARRI(IM,IN,A,V,I)
c    Extract the i'th row of an array       
      INTEGER IM,IN
      DIMENSION A(IM,IN)
      DIMENSION V(IN)
      DO 20 I1=1,IN
20	V(I1)=A(I,I1)
      RETURN
      END
c
c
      SUBROUTINE ARRJ(IM,IN,A,VJ,IJ)
c    Extract the j'th column of an array  
      DIMENSION A(IM,IN)
      DIMENSION VJ(IM)
      DO 20 I1=1,IM
20	VJ(I1)=A(I1,IJ)
      RETURN
      END
c
c
C------MULTIPLY A REAL ARRAY BY A CONSTANT
      SUBROUTINE ARRMC(IM,IN,A1,C,A)
      DIMENSION A1(IM,IN),A(IM,IN)
      DO 10 I1=1,IM
      DO 10 I2=1,IN
10	A(I1,I2)=A1(I1,I2)*C
      END
c
c
C-------A STANDARD SUBPROGRAM TO MULTIPLY TWO ARRAYS
C	A1(IM1,IN1)*A2(IM2,IN2)=RSLT(IM1,IN2)
      SUBROUTINE ARRMULT(IM1,IN1,IM2,IN2,A1,A2,RSLT,V1,V2)
      INTEGER IM1,IM2,IN1,IN2
      DIMENSION A1(IM1,IN1),A2(IM2,IN2),RSLT(IM1,IN2)
      DIMENSION V1(IN1),V2(IM2)
      IF (IN1.NE.IM2) WRITE(6,*) 'The two arrays cannot be multiplied'
      DO 50 I=1,IM1
      CALL ARRI(IM1,IN1,A1,V1,I)
      DO 40 IJ=1,IN2
      CALL ARRJ(IM2,IN2,A2,V2,IJ)
40	RSLT(I,IJ)=POIMULT(IN1,IM2,V1,V2)
50	CONTINUE
      RETURN
      END

      SUBROUTINE ARRPS(IM,IN,A1,A2,OUT)
      DIMENSION OUT(IM,IN),A1(IM,IN),A2(IM,IN)
      DO 10 II=1,IM
      DO 10 IJ=1,IN
10	OUT(II,IJ)=A1(II,IJ)-A2(II,IJ)
      END
c
c
      subroutine arrrtoi(im,in,rmat,imat,err)
c    convert a matrix from real to integer
c      IM and IN are the dimensions.
c     rmat is input real IM*IN-dimension matrix
c     imat is output integer IM*IN-dimension matrix
c     err: if difference between real and integer values is bigger than
c           this parameter program will warn you
      real*4 rmat(im,in)
      integer imat(im,in)
c
      do j = 1, in
      do i = 1, im
      if (rmat(i,j).ge.0) then
       imat(i,j) = int(rmat(i,j)+0.5)
      else
       imat(i,j) = int(rmat(i,j)-0.5)
      end if
      if (abs(rmat(i,j)-float(imat(i,j))).gt.err) then
c	 WRITE(6,*) 'Warning: r-i > ',err
c	 WRITE(6,*)  'i,j,real,integer',i,j,rmat(i,j),imat(i,j)
      end if
      end do
      end do
      end
c      
c
      subroutine arrvalue(n,array,value)
c    initialise an n-dimensional array
      real array(n)
      do i = 1, n
        array(i) = value
      end do
      end
c
c
      SUBROUTINE AVERG(M,N,XIN,AVE)
C A routine to calculate the row averages of an M*N array XIN.
C The mean M-dimensional array is AVE which average the N data in XIN.
      DIMENSION XIN(M,N),AVE(M)
      XN = N
      DO I = 1, M
      AVE(I) = 0.
      DO J = 1, N
      AVE(I) = AVE(I) + XIN(I,J)
      END DO
      AVE(I) = AVE(I)/ XN
      END DO
      END 
c
c
      function bondangle(xyz1,xyz2,xyz3)
c subroutine to calculate a bond angle by giving 3 atoms. the output angle
c will be vector 2-1 and 2-3
c  1     3
c   \   /
c    \2/
c   bondangle
c xyz1,xyz2,xyz3 are the 3 atom positions and bondangle is shown in the fig
c	real xyz1(3),xyz2(3),xyz3(3)
      real xyz1(3),xyz2(3),xyz3(3)
      real b1(3),b2(3)
      call arrps(3,1,xyz1,xyz2,b1)
      call arrps(3,1,xyz3,xyz2,b2)
      bondangle = angle(b1,b2)
      end 
c
c
      function bonddihed(xyz1,xyz2,xyz3,xyz4)
c subroutine to calculate a dihedral angle by giving 4 atoms. 
c the output angle will be between plane of 1-2-3 and 2-3-4
c  1         4                      1        
c   \       /                        \       
c    \2---3/                          \2---3\  bonddihed=180 in this case
c   bonddihed= 0 in this case                \
c                                             4
c xyz1,xyz2,xyz3 xyz4 are the 4 atom positions and bonddihed is 
c shown in the fig
c	real xyz1(3),xyz2(3),xyz3(3),xyz4(3)
      real xyz1(3),xyz2(3),xyz3(3),xyz4(3)
c	real b1(3),b2(3)
      real b321(3),b432(3)
      real b23(3),b32(3),b34(3),b21(3)
      real by(3),bdih(3),bd(2)
      real view(3,3),views(3,3)
      external acosd, asind
c
      call arrps(3,1,xyz3,xyz2,b23)
      call arrps(3,1,xyz1,xyz2,b21)
      call veccrsmlt(b23,b21,b321)
      call arrps(3,1,xyz2,xyz3,b32)
      call arrps(3,1,xyz4,xyz3,b34)
      call veccrsmlt(b34,b32,b432)
c	call arrgive(3,b321,view(1,1))
c	call arrgive(3,b23,view(1,3))
      call veccrsmlt(b32,b321,by)
      call arrmc(3,1,b321,1./vem(3,b321),view(1,1))
      call arrmc(3,1,b32,1./vem(3,b32),view(1,3))
      call arrmc(3,1,by,1./vem(3,by),view(1,2))
      call antiarr(3,3,view,views)
      call matmult(3,3,3,1,views,b432,bdih)
      if (bdih(3).gt.1.0e-5) then
       WRITE(6,*) 'Warning in dihedral',bdih
      end if
c
      call arrmc(2,1,bdih,1./vem(2,bdih),bd)
c
c dump here
c	WRITE(6,*) 'views',views
c	WRITE(6,*) 'bdih',bdih
c	WRITE(6,*) 'bd',bd
c
      if (bd(1).ge.0. .and. bd(2).ge. 0.) then
       bonddihed = acosd(bd(1))
      else if (bd(1).ge.0. .and. bd(2).lt.0) then
       bonddihed = asind(bd(2))
      else if (bd(1).lt.0. .and. bd(2).ge.0) then
       bonddihed = acosd(bd(1))
      else if (bd(1).lt.0. .and. bd(2).lt.0) then
       bonddihed = -acosd(bd(1))
c	else
c	 stop ' what could be?'
      end if
      end
c
c
      subroutine caseres(res1,res2)
c A subroutine to change the case of amino acid e.g. PHE to Phe
c res1 input residue name
c res2 output residue name
c
      character*3 res1,res2,nam1(20),nam2(20)
c...  data statements.  Separate declaration and init required for f2c
      data nam1 /
     1 'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
     2 'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'/
      data nam2 /
     1 'Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile',
     2 'Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp','Tyr','Val'/
c
      do i = 1,20
       if (res1.eq.nam1(i)) then
        res2 = nam2(i)
        return
       end if
      end do
      res2 = res1
      end
c
c
      SUBROUTINE CLLZ(NSYM,SYM,CELL)
      DIMENSION SYM(3,4,NSYM),CELL(6)
      DO 20 IS=1,NSYM
      DO 20 IC=1,3
20	SYM(IC,4,IS)=SYM(IC,4,IS)*CELL(IC)
      END	
c
c
      SUBROUTINE COMPARE ( NATM, XYZ1, XYZ2, A, T )
C	DIMENSION XYZ1(3,NATM), XYZ2(3,NATM), A(3,3), T(3)
C A subroutine to compare two structure without superposition.
c MOL2 i.e. XYZ2(3*NATM) = A(3*3) * XYZ1(3*NATM) + T(3)
C Where    NATM represents the number of the atoms in each molecule
c                which will be superimposed.
c          XYZ1(3,NATM) represents the coordinates of MOL1
c          XYZ2(3,NATM) represents the coordinates of MOL2
C          A(3*3)       represents the output matrix.
c          T(3*3)       represents the output vector.
c
C NATM can not be larger than NATM0 (currently 50000) or the parameter NATM0
C should be modified in this routine. 
c
C	COMMON/SHIFS/ SCALR,SCALT
c 	DATA SCALR/1./,SCALT/1./
C Note: There should be a scale of the shifts. i.e. Shift= scale * shifts
c Here SCALR is the rotation shiftscale
c Here SCALT is the translation shiftscale
c The initial and suggested SCALR is 1.
c The initial and suggested SCALt is 1.
C
C	COMMON/RMS/ RMS,SMEAN,NREF,NREF1
C RMS is the final R.M.S between two molecules.
C SMEAN is the mean distance.
c NREF is the number of cycles of rotation refinement
c NREF1 is not used by this routine.
c They could be used to judge the SCALR and SCALS
c
C	COMMON/IAT/ IAT1,IAT2,IAT3,IAT
C	DATA SCALR/1./,IAT/0/
c                                           by Guoguang Lu
C                                               01/31/1995 at Purdue
C
      COMMON/RMS/ RMS,SMEAN,NREF,NREF1
c	COMMON/SHIFS/ SCALR,SCALS
c	COMMON/IAT/ IAT1,IAT2,IAT3,IAT
c	DATA SCALR/1./,SCALS/1./,IAT/0/,IAT1/1/,IAT2/2/,IAT3/3/
      PARAMETER (NATM0 = 50000)
c   NATM0 is the largest number of atoms.
      DIMENSION XYZ1(3,NATM), XYZ2(3,NATM), A(3,3), T(3)
      DIMENSION X1(3,NATM0),X2(3,NATM0),DIS(3,NATM0)
      DIMENSION CEN1(3),CEN2(3)
      DIMENSION B1(3),B2(3),B3(3),T0(3)
c
c	DEG=180./3.14159265359
c
C Number of atoms cannot be more than NATM0 or less than 3.
c
C
      CALL RTMOV(NATM,XYZ1,A,T,DIS)
      CALL ARRPS(3,NATM,DIS,XYZ2,DIS)
C
      ATM = NATM
      SMEAN = 0
      DO I=1, NATM
      D = VEM(3,DIS(1,I))
      SMEAN = SMEAN + D
      END DO
      SMEAN = SMEAN /ATM
C
c SMEAN = SIGMA (DISTANCE) / N
C           i        i
      RMS = SQRT ( DOSQ(NATM*3,DIS) / ATM )
C
c RMS   = SQRT( SIGMA (DISTANCE^2) / N )
C           i              i
c 
      WRITE(6,*)
      WRITE(6,*) 'R.M.S.'
      WRITE(6,*) '       natm'
      WRITE(6,*) 'SQRT( SIGMA (Distance(i)^2)/natm ) = ',RMS
      WRITE(6,*) '       i=1'
C
      WRITE(6,*)
      WRITE(6,*) 'Mean Distance:'
      WRITE(6,*) ' natm'
      WRITE(6,*) 'SIGMA (Distance(i)) /natm = ',SMEAN
      WRITE(6,*) ' i=1'
C
      WRITE(6,*)
      WRITE(6,*) 'Mol1 is superposed to Mol2.'
      WRITE(6,*) 'The matrix and the vector are:'
      WRITE(6,*) 
      WRITE(6,1) (A(1,J),J=1,3),CEN1(1),T0(1)
      WRITE(6,2) (A(2,J),J=1,3),CEN1(2),T0(2)
      WRITE(6,1) (A(3,J),J=1,3),CEN1(3),T0(3)
1	FORMAT(1X,'      (',3F10.6,' )   (     ',f10.5,' )   ('
     1 ,f10.5,' )')
2	FORMAT(1X,' X2 = (',3F10.6,' ) * ( X1 -',f10.5,' ) + ('
     1 ,f10.5,' )')
C
      CALL MATMULT(3,3,3,1,A,CEN1,B3)
      CALL ARRPS(3,1,T0,B3,B3)
      CALL ARRAD(3,1,CEN1,B3,T)
C
      WRITE(6,*) 
      WRITE(6,*) 
      WRITE(6,4) (A(1,J),J=1,3),T(1)
      WRITE(6,5) (A(2,J),J=1,3),T(2)
      WRITE(6,4) (A(3,J),J=1,3),T(3)
4	FORMAT(1X,'      (',3F10.6,' )   (    )   (',f10.5,' )')
5	FORMAT(1X,' X2 = (',3F10.6,' ) * ( X1 ) + (',f10.5,' )')
C
      END 
c
c
      SUBROUTINE CRSTLARR(CELL,CRST,CRIS)
      DIMENSION CRST(3,3),CRIS(3,3),BUF1(3),BUF2(3),
     1IUF(3),CELL(6)
      INTEGER IUF
      external cosd, sind
      AP=CELL(4)
      BA=CELL(5)
      GA=CELL(6)
      CALL ARRMC(3,1,CRIS,0.,CRIS)
      COASTAR=(COSD(BA)*COSD(GA)-COSD(AP))/(SIND(BA)
     1*SIND(GA))
      SIASTAR=SQRT(1.-COASTAR*COASTAR)
      CRIS(1,1)=1.	
      CRIS(1,2)=COSD(GA)
      CRIS(1,3)=COSD(BA)
      CRIS(2,2)=SIND(GA)
      CRIS(2,3)=-SIND(BA)*COASTAR
      CRIS(3,3)=SIND(BA)*SIASTAR
      CALL ARRMC(3,3,CRIS,1.,CRST)
      CALL IVSN(3,CRST,BUF1,BUF2,IUF,VAL,1E-6)
      RETURN
      END
c
c
      SUBROUTINE CRTMOVE(NATM,XIN,A,T0,T,XOUT)
C A routine to compute: XOUT = A * (XIN-T0) + T
C        where   XIN is the input 3*NATM array
C                XOUT is the output 3*NATM array
c                A(3,3) is a 3*3 rotation matrix
c                T0(3) is a 3-dimensional translational vector.
c                T(3) is a 3-dimensional translational vector.
C Note: XIN has been changed into XIN-T0 after this subroutine.
c 89/11/06, BMC,UU,SE
c	DIMENSION XIN(3,NATM),XOUT(3,NATM),A(3,3),T(3)
C
      DIMENSION XIN(3,NATM),XOUT(3,NATM),A(3,3),T(3),T0(3)
      DO I = 1, NATM
      CALL ARRPS(3,1,XIN(1,I),T0,XIN(1,I))
      END DO
      CALL MATMULT(3,3,3,NATM,A,XIN,XOUT)
      DO I = 1, NATM
      CALL ARRAD(3,1,XOUT(1,I),T,XOUT(1,I))
      END DO
      END
c
c
      SUBROUTINE LGG_CRYSTAL(CELL,CELLS,DEOR,ORTH,DEORS,ORTHS)
C PJX - renamed from CRYSTAL to avoid clashes with similarly named
C common block in mapmask, ncsmask and dm
C
C this is a routine which inputs a cell parameter, CELL and calculates
c CELLS*6 --- cell dimensions in reciprocal space
c DEOR3*3 --- Deorthogonalization matrix in real space 
c ORTH3*3 --- Orthogonalization matrix in recepical space.
c DEORS3*3 --- Deorthogonalization matrix in reciprocal space 
c ORTHS3*3 --- Orthogonalization matrix in reciprocal space.
c All the matrices have cell dimension information in them.
c Guoguang 931118
c
      REAL DEOR(3,3),ORTH(3,3),CELL(6)
      REAL DEORS(3,3),ORTHS(3,3),CELLS(6)
      REAL BUFF(3,3)
      DIMENSION B1(3),B2(3),IUF(3)
      external cosd, sind
      AP=CELL(4)
      BA=CELL(5)
      GA=CELL(6)
      COASTAR=(COSD(BA)*COSD(GA)-COSD(AP))/(SIND(BA)
     1*SIND(GA))
      SIASTAR=SQRT(1.-COASTAR*COASTAR)
      CALL ARRVALUE(9,ORTH,0.)
      ORTH(1,1)=CELL(1)	
      ORTH(1,2)=CELL(2)*COSD(GA)
      ORTH(2,2)=CELL(2)*SIND(GA)
      ORTH(1,3)=CELL(3)*COSD(BA)
      ORTH(2,3)=-CELL(3)*SIND(BA)*COASTAR
      ORTH(3,3)=CELL(3)*SIND(BA)*SIASTAR
C
      CALL ARRGIVE(9,ORTH,DEOR)
      CALL IVSN(3,DEOR,B1,B2,IUF,DE,1.0E-6)
      CALL ANTIARR(3,3,ORTH,DEORS)
      CALL ANTIARR(3,3,DEOR,ORTHS)
C
      DO I = 1, 3
      CELLS(I) = VEM(3,orths(1,I))
      END DO
C
      CELLS(4) = ANGLE(ORTHS(1,2),ORTHS(1,3))	
      CELLS(5) = ANGLE(ORTHS(1,3),ORTHS(1,1))	
      CELLS(6) = ANGLE(ORTHS(1,1),ORTHS(1,2))	
      RETURN
      END
c
c
      SUBROUTINE CRYSTREC(CELL,CELLS,DEOR,ORTH,DEORS,ORTHS)
C this is a routine which inputs a cell parameter, CELL and calculates
c matrices for deorthogonalization and orthgonalization
c The convention in this routine is x=a* y in plane of a*-b*, z is C direction
c CELLS*6 --- cell dimensions in reciprocal space
c DEOR3*3 --- Deorthogonalization matrix in real space 
c ORTH3*3 --- Orthogonalization matrix in reciprocal space.
c DEORS3*3 --- Deorthogonalization matrix in reciprocal space 
c ORTHS3*3 --- Orthogonalization matrix in reciprocal space.
c All the matrices have cell dimension information in them.
c Guoguang 950914
c
      REAL DEOR(3,3),ORTH(3,3),CELL(6)
      REAL DEORS(3,3),ORTHS(3,3),CELLS(6)
      REAL BUFF(3,3),CELL1(6)
      DIMENSION B1(3),B2(3),IUF(3)
c Get a PDB convention 
      CALL LGG_CRYSTAL(CELL,CELLS,DEOR,ORTH,DEORS,ORTHS)
      CALL LGG_CRYSTAL(CELLS,CELL1,DEORS,ORTHS,DEOR,ORTH)
c
      end
c
c
c This is a subroutine to split one line into multiple lines by detecting
c a separation character for example ; | and so on
c the multiple lines will be written to unit iout
c iout -- unit to be written to
c il   -- length of the line
c line -- the line of characters to be split
c sep  -- "separation character"
c iline -- output number of lines after splitting.
c guoguang 940616
      subroutine cutline(iout,line,sep,iline)
      character*(*) line
      character*1 sep
      iline = 0
      is = 1
c	WRITE(6,*)  line
c	WRITE(6,*)  'sep ',sep
      il = lenstr(line)
      do it=il,1,-1
       if (ichar(line(it:it)).le.31) then
        line(it:it) = ' '
       else
        goto 9
       end if
      end do
9	il = it 
c	WRITE(6,*)  'il',il,line(1:il)
      rewind (iout)
      if (il.eq.0) return
10	ip = index(line(is:il),sep)
      if (ip.eq.0) then
       il1 = lenstr(line(is:il))
       if (il1.gt.0) then
        do isps = is, is+il1-1
         if (line(isps:isps).ne.' ') then
          write(iout,'(a)') line(isps:is+il1-1) 
          iline = iline + 1
          rewind (iout)
          return
         end if
        end do
        stop 'Strange'
       end if
      else 
       if (ip.gt.1) then
        il1 = lenstr(line(is:ip-2+is))
        if (il1.gt.0) then
         do isps = is, is+il1-1
          if (line(isps:isps).ne.' ') then
           write(iout,'(a)') line(isps:is+il1-1)
           iline = iline + 1
           is = ip + is
           if (is.le.il) goto 10
           rewind (10)
           return
          end if
         end do
         stop 'Strange'
        end if
       end if
       is = ip + is
       if (is.le.il) goto 10
       rewind (10)
       return
      end if
C	stop 'CUTLINE'
      end
c
c
      function dedx2(x1,x2,x3,y1,y2,y3,abc)
c   this is a program which calculates dy/dx(x2) by fitting 
c   a quadratic through 3 points. 
c   y = a*x^2 + b*x + c
c   dy/dx = 2*a*x + b
c  Here dedv = dy/dx(x2) = 2*a*x2 + b
c  When there are 3 points x1,y1, x2,y2, x3,y3, the coefficients a,b,c
c  can be obtained from a linear equation and dy/dx(x2) can thus be 
c  obtained.
c
c 
      dimension abc(3),y(3),arr(3,3)
      dimension b1(3),b2(3),m(3)
c
      y(1) = y1
      y(2) = y2
      y(3) = y3
      arr(1,1) = x1 * x1
      arr(2,1) = x2 * x2
      arr(3,1) = x3 * x3
      arr(1,2) = x1 
      arr(2,2) = x2 
      arr(3,2) = x3
      arr(1,3) = 1.
      arr(2,3) = 1.
      arr(3,3) = 1.
      call ivsn(3,arr,b1,b2,m,de,1e-5)
      call matmult(3,3,3,1,arr,y,abc)
      dedx2 = 2.*abc(1)*x2 + abc(2)
      return
      end
c
c
c a subroutine to convert denzo missetting angle to a matrix
c definition [rot]=[rotx]*[roty]*[rotz]
c 
      subroutine denmis(phi,rot)
      real phi(3),rot(3,3)
      external cosd, sind
      sinx = sind(phi(1))
      cosx = cosd(phi(1))
      siny = sind(phi(2))
      cosy = cosd(phi(2))
      sinz = sind(phi(3))
      cosz = cosd(phi(3))
      rot(1,1) =  cosy*cosz
      rot(2,1) = -cosx*sinz+sinx*siny*cosz
      rot(3,1) =  sinx*sinz+cosx*siny*cosz
      rot(1,2) =  cosy*sinz
      rot(2,2) =  cosx*cosz+sinx*siny*sinz
      rot(3,2) = -sinx*cosz+cosx*siny*sinz
      rot(1,3) = -siny
      rot(2,3) =  sinx*cosy
      rot(3,3) =  cosx*cosy
      end	
c
c
      function dihedral(x1,x2,x3,x4)
c calculate the dihedral angle of four input positions
c x1, x2, x3, x4 are the four input positions
c dihedral is the output dihedral angle in degrees
      real*4 x1(3),x2(3),x3(3),x4(3)
      real*4 b1(3),b2(3),b3(3),b4(3)
c
      call arrps(3,1,x1,x2,b1)	
      call arrps(3,1,x3,x2,b2)	
      call veccrsmlt(b1,b2,b3)
      call arrps(3,1,x2,x3,b1)	
      call arrps(3,1,x4,x3,b2)
      call veccrsmlt(b1,b2,b4)
      dihedral = angle(b3,b4)
      end
c
c
      FUNCTION DIST(XYZ1,XYZ2)
c Subroutine to calculate distance between two positions
      REAL xyz1(3),xyz2(3)
      DX = XYZ1(1)-XYZ2(1)
      DY = XYZ1(2)-XYZ2(2)
      DZ = XYZ1(3)-XYZ2(3)
      DIST = SQRT(DX*DX + DY*DY + DZ*DZ)
      END 
c
c
      FUNCTION DOSQ(N,V)
c	DIMENSION V(N)
C DOSQ = SIGMA ( V(i)^2 )
C          i
      DIMENSION V(N)
      DOSQ = 0.
      DO 10 I=1,N
10	DOSQ = DOSQ + V(I)*V(I)
      RETURN
      END
c
c
c------------------------------------------------------------
      subroutine down(txt,len)
c
c convert character string to lower case
c
      character*(*) txt
      character*80 save
      character*26 ualpha,lalpha
      data lalpha /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data ualpha /'abcdefghijklmnopqrstuvwxyz'/

      do 9000 i=1,len
        if((txt(i:i).ge.'A').and.(txt(i:i).le.'Z')) then
          match = index(lalpha,txt(i:i))
          save(i:i) = ualpha(match:match)
        else
          save(i:i) = txt(i:i)
        end if
c      end do
9000	continue

      txt = save
      return
      end
c
c
      SUBROUTINE DRVRBD(ITH,ROH,A)
C This is a subroutine to calculate the matrix of deriv(E)/deriv(theta i)
c where E is a 3*3 rotation matrix and theta i (i=1,2,3) is Eulerian Angle
c in ROH system.
C	DIMENSION A(3,3),ROH(3)
c parameter description:
c ITH:  when ITH = 1
c                     the output matrix is deriv(E) /deriv(theta1)
c       when ITH = 2
c                     the output matrix is deriv(E) /deriv(theta2)
c       when ITH = 3
c                     the output matrix is deriv(E) /deriv(theta3)
c ROH are the three Eulerian angles
c A  is the output deriv matrix. i.e. A = deriv(E) /deriv)(theta i)
c                               written by Guoguang Lu
c
      DIMENSION A(3,3),ROH(3)
      DIMENSION A1(3,3)
      ARC=3.14159265359/180.0
      SIPS = SIN(ROH(1)*ARC)
      COPS = COS(ROH(1)*ARC)
      SITH = SIN(ROH(2)*ARC)
      COTH = COS(ROH(2)*ARC)
      SIPH = SIN(ROH(3)*ARC)
      COPH = COS(ROH(3)*ARC)
c
c	Rossmann & Blow Matrix
c
c	E(1,1) =-SIPS*COTH*SIPH + COPS*COPH
c	E(1,2) =-SIPS*COTH*COPH - COPS*SIPH
c	E(1,3) =                  SIPS*SITH
c	E(2,1) = COPS*COTH*SIPH + SIPS*COPH
c	E(2,2) = COPS*COTH*COPH - SIPS*SIPH
c	E(2,3) =                - COPS*SITH
c	E(3,1) =                  SITH*SIPH
c	E(3,2) =                  SITH*COPH
c	E(3,3) = COTH


      IF (ITH.EQ.1) THEN
      A(1,1) =-coPS*COTH*SIPH - siPS*COPH
      A(1,2) =-coPS*COTH*COPH + siPS*SIPH
      A(1,3) =                  coPS*SITH
      A(2,1) =-siPS*COTH*SIPH + coPS*COPH
      A(2,2) =-siPS*COTH*COPH - coPS*SIPH
      A(2,3) =                + siPS*SITH
      A(3,1) = 0.
      A(3,2) = 0.
      A(3,3) = 0.
      ELSE IF (ITH.EQ.2) THEN
      A(1,1) =+SIPS*siTH*SIPH
      A(1,2) =+SIPS*siTH*COPH
      A(1,3) =                  SIPS*coTH
      A(2,1) =-COPS*siTH*SIPH
      A(2,2) =-COPS*siTH*COPH
      A(2,3) =                - COPS*coTH
      A(3,1) =                  coTH*SIPH
      A(3,2) =                  coTH*COPH
      A(3,3) = -siTH
      ELSE IF (ITH.EQ.3) THEN
      A(1,1) =-SIPS*COTH*coPH - COPS*siPH
      A(1,2) =+SIPS*COTH*siPH - COPS*coPH
      A(1,3) = 0.
      A(2,1) = COPS*COTH*coPH - SIPS*siPH
      A(2,2) =-COPS*COTH*siPH - SIPS*coPH
      A(2,3) = 0.
      A(3,1) =                  SITH*coPH
      A(3,2) =                 -SITH*siPH
      A(3,3) = 0.
      ELSE
      STOP 'invalid parameter ITH '
      END IF
C SINCE ROH IS IN DEGREES, THE MATRIX HAS TO BE MULTIPLIED BY ARC
      CALL ARRMC(3,3,A,ARC,A1)
      call arrgive(9,a1,a)
      END
c
c
      SUBROUTINE DRVROHTH(ITH,ROH,A)
C This is a subroutine to calculate the matrix of deriv(E)/deriv(theta i)
c where E is a 3*3 rotation matrix and theta i (i=1,2,3) is Eulerian Angle
c in ROH system.
C	DIMENSION A(3,3),ROH(3)
c parameter description:
c ITH:  when ITH = 1
c                     the output matrix is deriv(E) /deriv(theta1)
c       when ITH = 2
c                     the output matrix is deriv(E) /deriv(theta2)
c       when ITH = 3
c                     the output matrix is deriv(E) /deriv(theta3)
c ROH are the three Eulerian angles
c A  is the output deriv matrix. i.e. A = deriv(E) /deriv)(theta i)
c                               written by Guoguang Lu
c
      DIMENSION A(3,3),ROH(3)
      ARC=3.14159265359/180.0
      SIPS = SIN(ROH(1)*ARC)
      COPS = COS(ROH(1)*ARC)
      SITH = SIN(ROH(2)*ARC)
      COTH = COS(ROH(2)*ARC)
      SIPH = SIN(ROH(3)*ARC)
      COPH = COS(ROH(3)*ARC)
C
C
C *  ROH - MATRIX
C
C  110 E(1,1) = COPS*COPH - SIPS*SIPH*SITH
C      E(1,2) =-SIPS*COTH
C      E(1,3) = COPS*SIPH + SIPS*SITH*COPH
C      E(2,1) = SIPS*COPH + COPS*SITH*SIPH
C      E(2,2) = COPS*COTH
C      E(2,3) = SIPS*SIPH - COPS*SITH*COPH
C      E(3,1) =-COTH*SIPH
C      E(3,2) = SITH
C      E(3,3) = COTH*COPH

      IF (ITH.EQ.1) THEN
      A(1,1) = -SIPS*COPH - COPS*SIPH*SITH
      A(1,2) = -COPS*COTH
      A(1,3) = -SIPS*SIPH + COPS*SITH*COPH
      A(2,1) = COPS*COPH - SIPS*SITH*SIPH
      A(2,2) = -SIPS*COTH
      A(2,3) = COPS*SIPH + SIPS*SITH*COPH
      A(3,1) = 0.
      A(3,2) = 0.
      A(3,3) = 0.
      ELSE IF (ITH.EQ.2) THEN
      A(1,1) = - SIPS*SIPH*COTH
      A(1,2) = SIPS*SITH
      A(1,3) = SIPS*COTH*COPH
      A(2,1) = COPS*COTH*SIPH
      A(2,2) = -COPS*SITH
      A(2,3) = - COPS*COTH*COPH
      A(3,1) = SITH*SIPH
      A(3,2) = COTH
      A(3,3) = -SITH*COPH
      ELSE IF (ITH.EQ.3) THEN
      A(1,1) = -COPS*SIPH - SIPS*COPH*SITH
      A(1,2) = 0.
      A(1,3) = COPS*COPH - SIPS*SITH*SIPH
      A(2,1) =- SIPS*SIPH + COPS*SITH*COPH
      A(2,2) = 0.
      A(2,3) = SIPS*COPH + COPS*SITH*SIPH
      A(3,1) =-COTH*COPH
      A(3,2) = 0.
      A(3,3) = - COTH*SIPH
      ELSE
      STOP 'invalid parameter ITH '
      END IF
      END
c
c
      SUBROUTINE DRVROHTHARC(ITH,ROH,A)
C This is a subroutine to calculate the matrix of deriv(E)/deriv(theta i)
c where E is a 3*3 rotation matrix and theta i (i=1,2,3) is Eulerian Angle
c in ROH system.
C	DIMENSION A(3,3),ROH(3)
c parameter description:
c ITH:  when ITH = 1
c                     the output matrix is deriv(E) /deriv(theta1)
c       when ITH = 2
c                     the output matrix is deriv(E) /deriv(theta2)
c       when ITH = 3
c                     the output matrix is deriv(E) /deriv(theta3)
c ROH are the three Eulerian angles
c A  is the output deriv matrix. i.e. A = deriv(E) /deriv)(theta i)
c                               written by Guoguang Lu
c
      DIMENSION A(3,3),ROH(3)
C	ARC=3.14159265359/180.0
C	ARC= 1.0
      SIPS = SIN(ROH(1))
      COPS = COS(ROH(1))
      SITH = SIN(ROH(2))
      COTH = COS(ROH(2))
      SIPH = SIN(ROH(3))
      COPH = COS(ROH(3))
C
C
C *  ROH - MATRIX
C
C  110 E(1,1) = COPS*COPH - SIPS*SIPH*SITH
C      E(1,2) =-SIPS*COTH
C      E(1,3) = COPS*SIPH + SIPS*SITH*COPH
C      E(2,1) = SIPS*COPH + COPS*SITH*SIPH
C      E(2,2) = COPS*COTH
C      E(2,3) = SIPS*SIPH - COPS*SITH*COPH
C      E(3,1) =-COTH*SIPH
C      E(3,2) = SITH
C      E(3,3) = COTH*COPH

      IF (ITH.EQ.1) THEN
      A(1,1) = -SIPS*COPH - COPS*SIPH*SITH
      A(1,2) = -COPS*COTH
      A(1,3) = -SIPS*SIPH + COPS*SITH*COPH
      A(2,1) = COPS*COPH - SIPS*SITH*SIPH
      A(2,2) = -SIPS*COTH
      A(2,3) = COPS*SIPH + SIPS*SITH*COPH
      A(3,1) = 0.
      A(3,2) = 0.
      A(3,3) = 0.
      ELSE IF (ITH.EQ.2) THEN
      A(1,1) = - SIPS*SIPH*COTH
      A(1,2) = SIPS*SITH
      A(1,3) = SIPS*COTH*COPH
      A(2,1) = COPS*COTH*SIPH
      A(2,2) = -COPS*SITH
      A(2,3) = - COPS*COTH*COPH
      A(3,1) = SITH*SIPH
      A(3,2) = COTH
      A(3,3) = -SITH*COPH
      ELSE IF (ITH.EQ.3) THEN
      A(1,1) = -COPS*SIPH - SIPS*COPH*SITH
      A(1,2) = 0.
      A(1,3) = COPS*COPH - SIPS*SITH*SIPH
      A(2,1) =- SIPS*SIPH + COPS*SITH*COPH
      A(2,2) = 0.
      A(2,3) = SIPS*COPH + COPS*SITH*SIPH
      A(3,1) =-COTH*COPH
      A(3,2) = 0.
      A(3,3) = - COTH*SIPH
      ELSE
      STOP 'invalid parameter ITH '
      END IF
      END
c
c
      SUBROUTINE DRVROHTHD(ITH,ROH,A)
C This is a subroutine to calculate the matrix of deriv(E)/deriv(theta i)
c where E is a 3*3 rotation matrix and theta i (i=1,2,3) is Eulerian Angle
c in ROH system.
C	DIMENSION A(3,3),ROH(3)
c parameter description:
c ITH:  when ITH = 1
c                     the output matrix is deriv(E) /deriv(theta1)
c       when ITH = 2
c                     the output matrix is deriv(E) /deriv(theta2)
c       when ITH = 3
c                     the output matrix is deriv(E) /deriv(theta3)
c ROH are the three Eulerian angles
c A  is the output deriv matrix. i.e. A = deriv(E) /deriv)(theta i)
c                               written by Guoguang Lu
c
      DIMENSION A(3,3),ROH(3)
      DIMENSION A1(3,3)
      ARC=3.14159265359/180.0
      SIPS = SIN(ROH(1)*ARC)
      COPS = COS(ROH(1)*ARC)
      SITH = SIN(ROH(2)*ARC)
      COTH = COS(ROH(2)*ARC)
      SIPH = SIN(ROH(3)*ARC)
      COPH = COS(ROH(3)*ARC)
C
C
C *  ROH - MATRIX
C
C  110 E(1,1) = COPS*COPH - SIPS*SIPH*SITH
C      E(1,2) =-SIPS*COTH
C      E(1,3) = COPS*SIPH + SIPS*SITH*COPH
C      E(2,1) = SIPS*COPH + COPS*SITH*SIPH
C      E(2,2) = COPS*COTH
C      E(2,3) = SIPS*SIPH - COPS*SITH*COPH
C      E(3,1) =-COTH*SIPH
C      E(3,2) = SITH
C      E(3,3) = COTH*COPH

      IF (ITH.EQ.1) THEN
      A(1,1) = -SIPS*COPH - COPS*SIPH*SITH
      A(1,2) = -COPS*COTH
      A(1,3) = -SIPS*SIPH + COPS*SITH*COPH
      A(2,1) = COPS*COPH - SIPS*SITH*SIPH
      A(2,2) = -SIPS*COTH
      A(2,3) = COPS*SIPH + SIPS*SITH*COPH
      A(3,1) = 0.
      A(3,2) = 0.
      A(3,3) = 0.
      ELSE IF (ITH.EQ.2) THEN
      A(1,1) = - SIPS*SIPH*COTH
      A(1,2) = SIPS*SITH
      A(1,3) = SIPS*COTH*COPH
      A(2,1) = COPS*COTH*SIPH
      A(2,2) = -COPS*SITH
      A(2,3) = - COPS*COTH*COPH
      A(3,1) = SITH*SIPH
      A(3,2) = COTH
      A(3,3) = -SITH*COPH
      ELSE IF (ITH.EQ.3) THEN
      A(1,1) = -COPS*SIPH - SIPS*COPH*SITH
      A(1,2) = 0.
      A(1,3) = COPS*COPH - SIPS*SITH*SIPH
      A(2,1) =- SIPS*SIPH + COPS*SITH*COPH
      A(2,2) = 0.
      A(2,3) = SIPS*COPH + COPS*SITH*SIPH
      A(3,1) =-COTH*COPH
      A(3,2) = 0.
      A(3,3) = - COTH*SIPH
      ELSE
      STOP 'invalid parameter ITH '
      END IF
C SINCE ROH IS IN DEGREES, THE MATRIX HAS TO MULTIPLIED BY ARC
      CALL ARRMC(3,3,A,ARC,A1)
      CALL arrgive(9,A1,A1)
      END
c
c
      function dstll2(x1,x2,y1,y2)
c a routine to calculate the line-line distance.
c input:
c x1 and x2 are two positions on line X
c y1 and y2 are two positions on line Y
c output
c  dstll1 is the distance between x and y 
      real x1(3),x2(3)
      real y1(3),y2(3)
      real b1(3),b2(3),b3(3),b4(3),b5(3),b6(3)
c
      call arrps(3,1,x2,x1,b1)
      call arrps(3,1,y2,y1,b2)
      call veccrsmlt(b1,b2,b3)
      if (vem(3,b3).ge.1e-6) then
       dstll2 = dstps1(b3,x1,y1)
      else
       call arrps(3,1,y1,x1,b5)
       call arrmc(3,1,b1,1./vem(3,b1),b4)
       d1 = poimult(3,3,b5,b4)
       d2 = vem(3,b5)
       dstll2 = sqrt(d2*d2-d1*d1)
      end if
      end 
c
c
      FUNCTION DSTPL1(VL,XYZ0,XYZ,VPL)
C-----------DISTANCE FROM A POINT TO LINE	
c	FUNCTION DSTPL1(VL,XYZ0,XYZ,VPL)
c  VL -3- is the input direction vector of the line
c  XYZ0 -3- is the input position of a point on the line
c  XYZ -3-  is the input coordinates of the point.
c  VPL -3-  is the output vector from the point to the line which is
c           perpendicular to the direction vector of the line.
c  DISTPL1 is the output distance from the point to the line i.e. /VPL(3)/
C  so the equation of the line should be 
C         x-XYZ0(1)      y-XYZ0(2)       z-XYZ0(3)
C        ----------- =  -----------  =  -----------
C          VL(1)          VL(2)           VL(3)
C        where x,y,z is any position on the line
c

      DIMENSION VL(3),XYZ0(3),XYZ(3),VPL(3)
      DIMENSION B(3),BUF(3)
      CALL ARRPS(3,1,XYZ,XYZ0,BUF)
      T=POIMULT(3,3,BUF,VL)/POIMULT(3,3,VL,VL)
      CALL ARRMC(3,1,VL,T,BUF)
      CALL ARRAD(3,1,BUF,XYZ0,B)
      CALL ARRPS(3,1,XYZ,B,VPL)
      DSTPL1=VEM(3,VPL)
      RETURN
      END
c
c
      FUNCTION DSTPL2(VL,XYZ0,XYZ,VP0,VPL)
C-----------DISTANCE FROM A POINT TO LINE	
c	FUNCTION DSTPL1(VL,XYZ0,XYZ,VP0,VPL)
c  VL -3- is the input directional vector of the line
c  XYZ0 -3- is the a input position which belong to the line
c  XYZ -3-  is the input coordinate of the point.
c  VP0 -3-  is the output coordinate whis is the image of the point in the
c  line 
C  VPL -3- is the vector from VP0 to XYZ
c  DISTPL1 is the output distance from the point to the line i.e. /VPL(3)/
C  so the equation of the line should be 
C         x-XYZ0(1)      y-XYZ0(2)       z-XYZ0(3)
C        ----------- =  -----------  =  -----------
C          VL(1)          VL(2)           VL(3)
C        where the x,y,z is any position of the line
c

      DIMENSION VL(3),XYZ0(3),XYZ(3),VP0(3)
      DIMENSION VPL(3),BUF(3)
      CALL ARRPS(3,1,XYZ,XYZ0,BUF)
      DVL = VEM(3,VL)
      T= POIMULT(3,3,BUF,VL)/DVL*DVL
      CALL ARRMC(3,1,VL,T,BUF)
      CALL ARRAD(3,1,BUF,XYZ0,VP0)
      CALL ARRPS(3,1,XYZ,VP0,VPL)
      DSTPL2=VEM(3,VPL)
      RETURN
      END
c
c
      FUNCTION DSTPLROTR(A,T,XYZ)
C-----------DISTANCE FROM A POINT TO A ROTATION axis defined by ROT,T
c	FUNCTION DSTPL1(VL,XYZ0,XYZ,VPL)
c  A 3*3 is the rotation matrix which defines the axis
c  T  3 is the translation vector which defines the axis
c When there is a non-crystallographic rotation and translation
c x2 = A*x1 + t
c where A(3,3) is the 3*3 matrix
c 	t(3) is translation vector.
c  XYZ -3-  is the input coordinate of the point.
c  DISTPLROTR is the output distance from the point to the line i.e. /VPL(3)/
C  so the equation of the line should be 
C        where the x,y,z is any position of the line
c		Guoguang 970921
c
      real a(3,3),t(3),xyz(3),POL(3),xyz0(3),vpl(3)
      real B(3),BUF(3),vl(3),vl0(3)
      real b1(3),b2(3),x1(3),x2(3)
      external sind
c
      call arrvalue(3,xyz0,0.)
      call mtovec(a,vl,vkapa)
      vlm = vem(3,vl)
      if (vlm.lt.1.e-4) then
       WRITE(6,*) 'a',a
       WRITE(6,*) 'vlm',vlm
       stop 'something strange in dstplrotr'
      end if
      if (abs(vlm-1.).gt.1.e-5) then
       call arrmc(3,1,vl,1./vlm,vl0)
       call arrgive(3,vl0,vl)
      end if
      if (abs(vkapa).lt.1e-4) then
       go = vem(3,t)
       dstplrotr=9999.
      else 
       go = poimult(3,3,t,vl)
      end if
      call arrmc(3,1,vl,go,b1)
      call rtmov(3,xyz,a,t,x1)
      call arrps(3,1,x1,b1,x2)
      dstplrotr = dist(x2,xyz)/(2.*sind(vkapa/2.))
      end
c
c
      FUNCTION DSTPS1(VS,XYZ0,XYZ)
C------DISTANCE FROM A POINT TO A SECTION
c	FUNCTION DSTPS1(VL,XYZ0,XYZ,VPL)
c  VL -3- is the input directional vector of the plane section
c  XYZ0 -3- is the a input position which belong to the section
c  XYZ -3-  is the input coordinate of the point.
c  DISTPL1 is the output distance from the point to the section
C  so the equation of the section should be 
C
c VS(1) * ( x - XYZ0(1) ) + VS(2) * ( y - XYZ0(2) ) + VS(3) ( z - XYZ0(3) ) = 0
C
C        where the x,y,z is any position of the section
c
      DIMENSION VS(3),XYZ0(3),XYZ(3),BUF(3)
      CALL ARRPS(3,1,XYZ,XYZ0,BUF)
      DSTPS1=POIMULT(3,3,VS,BUF)/VEM(3,VS)
      RETURN
      END
c
c
c----------------------------------------------------------
      subroutine elb(txt,len)
c
c eliminate leading blanks from a character string
c
      character*(*) txt
      character*80 save

      do 9000 i=1,len
        if(txt(i:i).ne.' ') then
          nonb = i
          go to 100
        end if
9000  continue
      return
100   continue
      save = txt(nonb:len)
      txt = save
      return
      end
c
c
      SUBROUTINE ELIZE(N,E)
      DIMENSION E(N,N)
      DO 10 I=1,N
      E(I,I)=1.
c	WRITE(6,*) I,E
      DO 5 IJ=1,N
      IF (I.NE.IJ) E(I,IJ)=0.
5	CONTINUE
10	CONTINUE
      RETURN
      END
c
c
      function error2d(xy1,xy2)
      real xy1(2),xy2(2)
      e1 = xy2(1)-xy1(1)
      e2 = xy2(2)-xy1(2)
      error2d = sqrt(e1*e1+e2*e2)
      end
c
c
c read in a file and copy to another file
c nin --- unit of input file
c nout -- unit of output file
c file -- name of the input file
      subroutine filein(nin,nout,file)
      character*(*) file
      character*132 line
      call ccpdpn(nin,file,'readonly','F',0,0)
5	read(nin,'(a)',end=10) line
      il = lenstr(line)
      write(nout,'(a)') line(1:il)
      goto 5
10	continue
      end
c
c
      SUBROUTINE FMATIN(N,X)
      REAL X(N)
      DO I = 1, N
       CALL FRCINSIDE(X(i))
      END DO
      END
c
c
      SUBROUTINE FRAC2ANG(XYZ)
C A subroutine to convert fractional coordinates xyz into orthogonal 
c coordinates in Angstrom units.
C 
C	SUBROUTINE FRAC2ANG(XYZ)
C	COMMON/CELL/ CELL,DEOR,ORTH
C	DIMENSION CELL(6),DEOR(3,3),ORTH(3,3),XYZ(3)
C XYZ is an input fractional coordinate 3-dim array.
c CELL: 6-dim array to save the cell constants
c DEOR is 3*3 array, representing the deorthogonalizing array.
c ORTH is 3*3 array, representing the orthogonalizing array.
c
      COMMON/CELL/ CELL,DEOR,ORTH
      DIMENSION CELL(6),DEOR(3,3),ORTH(3,3),XYZ(3)
      DIMENSION B1(3),B2(3),B3(3)
      B3(1) = XYZ(1) * CELL(1)
      B3(2) = XYZ(2) * CELL(2)
      B3(3) = XYZ(3) * CELL(3)
C
      CALL matmult(3,3,3,1,ORTH,B3,XYZ,B1)
      END
c
c
      SUBROUTINE FRCINSIDE(X)
10	CONTINUE
      IF (X.GE.1.) THEN
      X = MOD(X,1.) 
      ELSE IF (X.LT.0.) THEN
      X = MOD(X,1.) + 1.
      END IF
      END
c
c
      SUBROUTINE FRCTOANG(XYZ)
C A subroutine to convert fractional coordinates xyz into orthogonal 
c coordinates in Angstrom units.
C 
C	SUBROUTINE FRAC2ANG(XYZ)
C	COMMON/CELL/ CELL,DEOR,ORTH
C	DIMENSION CELL(6),DEOR(3,3),ORTH(3,3),XYZ(3)
C XYZ is an input fractional coordinate 3-dim array.
c CELL: 6-dim array to save the cell constants
c DEOR is 3*3 array, representing the deorthogonalizing array.
c ORTH is 3*3 array, representing the orthogonalizing array.
c obs: different from frac2ang, the deor and orth must be from 
c 	subroutine crystal
c
      COMMON/CELL/ CELL,DEOR,ORTH
      DIMENSION CELL(6),DEOR(3,3),ORTH(3,3),XYZ(3)
      DIMENSION B1(3),B2(3),B3(3)
c	B3(1) = XYZ(1) * CELL(1)
c	B3(2) = XYZ(2) * CELL(2)
c	B3(3) = XYZ(3) * CELL(3)
C
      B3(1) = XYZ(1)
      B3(2) = XYZ(2)
      B3(3) = XYZ(3)
      CALL matmult(3,3,3,1,ORTH,B3,XYZ,B1)
      END
      character*80 function getnam(filnam)
      character*(*) filnam
c	character*80 cmd
c	call spstrunct(filnam)
c	il = lenstr(filnam)
c	write(cmd,'(3a)') 'printenv ',filnam(1:il),' > GETNAM'
cc	call system(cmd)
c       call ccpdpn(31,'GETNAM','old','F',0,0)
c	read(31,'(a)') getnam
c	call spstrunct(getnam)
c	close (31)
      getnam = ' '
      call ugtenv(filnam,getnam)
      end
c
c
      SUBROUTINE GETPDB(X,ATOM,RES,SEQ,NATM,FILNAM)
C A subroutine to read the Protein data bank format coordinate file
c X 3*21000 arrays, output xyz coordinates
c ATOM 21000 character*4 arrays, output atom names
c RES 21000 character*4 arrays, output residue number e.g. A101 A103...
C  If residue number is like A   1 then it will become A1
c SEQ 21000 character*4 arrays, output peptide name. e.g. PHE TRY ....
C NATM is the number of atoms in the coordinate file
c FILNAT is the file name of the coordinate.
c
      DIMENSION X(3,21000)
      CHARACTER*4 ATOM(21000),RES(21000),SEQ(21000),RES2
      CHARACTER FILNAM*80,HEAD*6,RES1*5
C	DATA FILNAM/'DIAMOND'/
C
      CALL SPSTRUNCT(FILNAM)
        call ccpdpn(3,FILNAM,'READONLY','F',0,0)
C
      NATM = 1
5	READ(3,10,ERR=5,END=20) HEAD, ATOM(NATM), SEQ(NATM), RES1
     1 , (X(I,NATM),I=1,3)
10	FORMAT(A6,7X,2A4,A5,4X,3F8.3)
      IF (HEAD.NE.'ATOM  '.AND.HEAD.NE.'HEDATM') GOTO 5
C
      IF (RES1(1:1).NE.' '.AND.RES1(2:2).EQ.' ') THEN
      RES2(1:1) = RES1(1:1)
      RES2(2:4) = RES1(3:5)
      ELSE
      RES2 = RES1(2:5)
      END IF
C
      call spstrunct(atom(natm))
      call spstrunct(seq(natm))
      call spstrunct(RES2)
      res(natm)=res2
C
12	CONTINUE
C
      NATM = NATM + 1
      IF (NATM.GT.21000) STOP 'Atoms can not be more than 21000.'
      GOTO 5
20	CONTINUE
      NATM = NATM - 1
      CLOSE (3)
C
      END
c
c
      SUBROUTINE GETPDB1(X,ATOM,SEQ,RES,BFAC,OCC,NATM,FILNAM)
C A subroutine to read the Protein data bank format coordinate file
c X 3*15000 arrays, output xyz coordinates
c ATOM 15000 character*4 arrays, output atom names
c RES 15000 character*4 arrays, output residue number e.g. A101 A103...
C  If residue number is like A   1 then it will become A1
c SEQ 15000 character*4 arrays, output peptide name. e.g. PHE TRY ....
C NATM is the number of atoms in the coordinate file
c FILNAT is the file name of the coordinate.
c
c in this subroutine, hydrogen atoms were excluded.
c
      PARAMETER (MAXATM = 25000)
      DIMENSION X(3,MAXATM)
      real bfac(MAXATM),occ(MAXATM)
      CHARACTER*4 ATOM(MAXATM),RES(MAXATM),SEQ(MAXATM),RES2
      CHARACTER FILNAM*80,HEAD*6,RES1*5
C	DATA FILNAM/'DIAMOND'/
C
      CALL SPSTRUNCT(FILNAM)
c	WRITE(6,*)  'open file'
c:	WRITE(6,*)  filnam
        call ccpdpn(3,FILNAM,'READONLY','F',0,0)
C
      NATM = 1
5	READ(3,10,ERR=5,END=20) HEAD, ATOM(NATM), SEQ(NATM), RES1
     1 , (X(I,NATM),I=1,3),occ(natm),bfac(natm)
c normal format
c10	FORMAT(A6,7X,2A4,A5,4X,3F8.3,2f6.2)
c xplor format
10	FORMAT(A6,6X,2A4,1X,A5,4X,3F8.3,2f6.2)
      IF (HEAD.NE.'ATOM  '.AND.HEAD.NE.'HETATM') GOTO 5
      IF (ATOM(NATM)(1:1).EQ.'H'.or.ATOM(NATM)(1:2).EQ.' H') goto 5
C
      IF (RES1(1:1).NE.' '.AND.RES1(2:2).EQ.' ') THEN
      RES2(1:1) = RES1(1:1)
      RES2(2:4) = RES1(3:5)
      ELSE
      RES2 = RES1(2:5)
      END IF
C
      call spstrunct(atom(natm))
      call spstrunct(seq(natm))
      call spstrunct(RES2)
      res(natm)=res2
C
12	CONTINUE
C
      NATM = NATM + 1
      IF (NATM.GT.MAXATM) STOP 'Atoms can not be more than MAXATM.'
      GOTO 5
20	CONTINUE
      NATM = NATM - 1
      CLOSE (3)
C
      END
c
c
      SUBROUTINE GETPDB2(X,ATOM,SEQ,RES,BFAC,OCC,NATM,FILNAM)
C A subroutine to read the Protein data bank format coordinate file
c X 3*maxatom arrays, output xyz coordinates
c ATOM maxatom character*4 arrays, output atom names
c RES maxatom character*4 arrays, output residue number e.g. A101 A103...
C  If residue number is like A   1 then it will become A1
c SEQ maxatom character*4 arrays, output peptide name. e.g. PHE TRY ....
C NATM is the number of atoms in the coordinate file
c FILNAT is the file name of the coordinate.
c
c in this subroutine, hydrogen atoms were included.
c
c	getpdb1 is for xplor file and exclude Hydrogen
c	getpdb2 is for any file
      parameter (maxatom = 50000)
      DIMENSION X(3,maxatom)
      real bfac(maxatom),occ(maxatom)
      CHARACTER*4 ATOM(maxatom),RES(maxatom),SEQ(maxatom),RES2
      CHARACTER FILNAM*80,HEAD*6,RES1*5
C	DATA FILNAM/'DIAMOND'/
C
      CALL SPSTRUNCT(FILNAM)
c	WRITE(6,*)  'open file'
c:	WRITE(6,*)  filnam
        call ccpdpn(3,FILNAM,'READONLY','F',0,0)
C
      NATM = 1
5	READ(3,10,ERR=5,END=20) HEAD, ATOM(NATM), SEQ(NATM), RES1
     1 , (X(I,NATM),I=1,3),occ(natm),bfac(natm)
c normal format
c10	FORMAT(A6,7X,2A4,A5,4X,3F8.3,2f6.2)
c xplor format
10	FORMAT(A6,6X,2A4,1X,A5,4X,3F8.3,2f6.2)
      IF (HEAD.NE.'ATOM  '.AND.HEAD.NE.'HETATM') GOTO 5
c	IF (ATOM(NATM)(1:1).EQ.'H'.or.ATOM(NATM)(1:2).EQ.' H') goto 5
C
      IF (RES1(1:1).NE.' '.AND.RES1(2:2).EQ.' ') THEN
      RES2(1:1) = RES1(1:1)
      RES2(2:4) = RES1(3:5)
      ELSE
      RES2 = RES1(2:5)
      END IF
C
      call spstrunct(atom(natm))
      call spstrunct(seq(natm))
      call spstrunct(RES2)
      res(natm)=res2
C
12	CONTINUE
C
      NATM = NATM + 1
      IF (NATM.GT.maxatom) STOP 'Atoms can not be more than maxatom.'
      GOTO 5
20	CONTINUE
      NATM = NATM - 1
      CLOSE (3)
C
      END
c
c
      SUBROUTINE GETPDB3(X,ATOM,SEQ,CH,RES,BFAC,OCC,NATM,FILNAM)
C A subroutine to read the Protein data bank format coordinate file
c X 3*15000 arrays, output xyz coordinates
c ATOM 15000 character*4 arrays, output atom names
c RES 15000 character*4 arrays, output residue number e.g. A101 A103...
C  If residue number is like A   1 then it will become A1
c SEQ 15000 character*4 arrays, output peptide name. e.g. PHE TRY ....
C NATM is the number of atoms in the coordinate file
c FILNAT is the file name of the coordinate.
c
c in this subroutine, hydrogen atoms were included.
c
c	getpdb1 is for xplor file and exclude Hydrogen
c	getpdb2 is for any file
      parameter (maxpro = 50000)
      DIMENSION X(3,maxpro)
      real bfac(maxpro),occ(maxpro)
      CHARACTER*4 ATOM(maxpro),RES(maxpro),SEQ(maxpro),RES2
      character*1 ch(maxpro)
      CHARACTER FILNAM*80,HEAD*6,RES1*5
C
      CALL SPSTRUNCT(FILNAM)
c	WRITE(6,*)  'open file'
c:	WRITE(6,*)  filnam
        call ccpdpn(3,FILNAM,'READONLY','F',0,0)
C
      NATM = 1
5	READ(3,10,ERR=5,END=20) HEAD, ATOM(NATM), SEQ(NATM), ch(natm),
     1  res(natm), (X(I,NATM),I=1,3),occ(natm),bfac(natm)
c normal format
c10	FORMAT(A6,7X,2A4,A5,4X,3F8.3,2f6.2)
c xplor format
10	FORMAT(A6,6X,2A4,1X,a1,A4,4X,3F8.3,2f6.2)
      IF (HEAD.NE.'ATOM  '.AND.HEAD.NE.'HETATM') GOTO 5
c	IF (ATOM(NATM)(1:1).EQ.'H'.or.ATOM(NATM)(1:2).EQ.' H') goto 5
C
c	IF (RES1(1:1).NE.' '.AND.RES1(2:2).EQ.' ') THEN
c	RES2(1:1) = RES1(1:1)
c	RES2(2:4) = RES1(3:5)
c	ELSE
c	RES2 = RES1(2:5)
c	END IF
C
      call spstrunct(atom(natm))
      call spstrunct(seq(natm))
c	call spstrunct(RES2)
c	res(natm)=res2
C
12	CONTINUE
C
      NATM = NATM + 1
      IF (NATM.GT.maxpro) then
       WRITE(6,*)  'Error: atom number is',natm
       WRITE(6,*)  'Error.. Atoms cannot be more than',maxpro
       stop
      end if
      GOTO 5
20	CONTINUE
      NATM = NATM - 1
      CLOSE (3)
C
      END
c
c
      SUBROUTINE GETPDBCA(X,ATOM,RES,SEQ,NATM,FILNAM)
C A subroutine to read the Protein data bank format coordinate file
c In this version only Ca atoms are read. -- Guoguang 950922
c X 3*50000 arrays, output xyz coordinates
c ATOM 50000 character*4 arrays, output atom names
c RES 50000 character*4 arrays, output residue number e.g. A101 A103...
C  If residue number is like A   1 then it will become A1
c SEQ 50000 character*4 arrays, output peptide name. e.g. PHE TRY ....
C NATM is the number of atoms in the coordinate file
c FILNAT is the file name of the coordinate.
c
      DIMENSION X(3,50000)
      CHARACTER*4 ATOM(50000),RES(50000),SEQ(50000),RES2
      CHARACTER FILNAM*80,HEAD*6,RES1*5
C	DATA FILNAM/'DIAMOND'/
C
      CALL SPSTRUNCT(FILNAM)
        call ccpdpn(3,FILNAM,'READONLY','F',0,0)
C
      NATM = 1
5	READ(3,10,ERR=5,END=20) HEAD, ATOM(NATM), SEQ(NATM), RES1
     1 , (X(I,NATM),I=1,3)
10	FORMAT(A6,7X,2A4,A5,4X,3F8.3)
      IF (HEAD(1:3).EQ.'END') GOTO 20
      IF (HEAD.NE.'ATOM  '.AND.HEAD.NE.'HEDATM') GOTO 5
      IF (ATOM(NATM).NE.'CA  ') GOTO 5
C
      IF (RES1(1:1).NE.' '.AND.RES1(2:2).EQ.' ') THEN
      RES2(1:1) = RES1(1:1)
      RES2(2:4) = RES1(3:5)
      ELSE
      RES2 = RES1(2:5)
      END IF
C
      call spstrunct(atom(natm))
      call spstrunct(seq(natm))
      call spstrunct(RES2)
      res(natm)=res2
C
12	CONTINUE
C
      NATM = NATM + 1
      IF (NATM.GT.50000) STOP 'Atoms cannot be more than 50000.'
      GOTO 5
20	CONTINUE
      NATM = NATM - 1
      CLOSE (3)
C
      END
c
c
      SUBROUTINE GETRDI(X,ATOM,RES,SEQ,NATM,FILNAM)
C A subroutine to read the Diamond format coordinate file
c X 3*5500 arrays, output xyz coordinates
c ATOM 5500 character*4 arrays, output atom names
c RES 5500 character*4 arrays, output residue number e.g. A101 A103...
C  If res is 'A  1'  then it will be truncated to 'A1  '
c SEQ 5500 character*4 arrays, output peptide name. e.g. PHE TRY ....
C NATM is the number of atoms in the coordinate file
c FILNAT is the file name of the coordinate.
c
      DIMENSION X(3,5500)
      CHARACTER*4 ATOM(5500),RES(5500),SEQ(5500),RES1*4
      CHARACTER FILNAM*80
C	DATA FILNAM/'DIAMOND'/
C
      CALL SPSTRUNCT(FILNAM)
        call ccpdpn(2,FILNAM,'READONLY','F',0,0)
C
      READ(2,'(A)') ATOM(1),ATOM(1),ATOM(1)	
C
      NATM = 1
5	READ(2,10,ERR=5,END=20) (X(I,NATM),I=1,3),
     1 SEQ(NATM), RES(NATM), ATOM(NATM)
10	FORMAT(3F10.4,35X,A3,A4,3X,A4)
      IF (RES(NATM).EQ.'END '.OR.RES(natm).EQ.'CAST') GOTO 5
C
      CALL SPSTRUNCT(ATOM(NATM))
      CALL SPSTRUNCT(SEQ(NATM))
      CALL SPSTRUNCT(RES(NATM))
C
12	NATM = NATM + 1
      GOTO 5
20	CONTINUE
      NATM = NATM - 1
      I = 1
      CLOSE (2)
C
25	IF (I.LE.NATM.AND.SEQ(I).NE.'    ') THEN
C
      DO J = I-1, 1, -1
      IF (J.LT.1) GOTO 30
      IF (RES(I).EQ.RES(J)) THEN
      SEQ(J) = SEQ(I)
      ELSE
      GOTO 30
      END IF
      END DO
C
30	DO J = I+1, NATM
      IF (RES(I).EQ.RES(J)) THEN
      SEQ(J) = SEQ(I)
      ELSE 
      GOTO 40
      END IF
      END DO
C
40	CONTINUE
      I = J 
      GOTO 25
      ELSE IF (I.LT.NATM) THEN
      I = I + 1
      GOTO 25
      ELSE
      RETURN
      END IF

      END
c
c
      SUBROUTINE GETSHELX ( X, ATOM, NATM, FILNAM )
c A routine to get coordinates and atom names from a SHELXX format
c file.
c X is 3*500 array to save the coordinates.
c ATOM is 4 bytes 500-dim array to represent the atom names
c NATM represents the number of atoms
c FILNAM is the SHELXX format file name.
c Sample coordinate format:
c FVAR   0.90195
c ND      5   0.59828   0.24644  10.24980  11.00000   0.01625   0.00675 =
c         0.03496   0.00372   0.00326   0.00215
c O1      4   0.58560   0.24368   0.55481  11.00000   0.02699   0.04010 =
c         0.01761   0.00117   0.00328   0.01509
c   ......
c   ......
c END
c 
      PARAMETER (MAXATM = 500 )
      CHARACTER FILNAM*80,ATOM(MAXATM)*4,KEY*80
      DIMENSION X(3,MAXATM)
C
C Open the SHELXX format coordinate file.
C
        call ccpdpn(1,FILNAM,'READONLY','F',0,0)
C	
      NATM = 0
      READ(1,*)
C
5	NATM = NATM + 1
      IF (NATM.GT.MAXATM) 
     +   STOP 'Atoms of Shelxx cannot be more than 500.'
      READ(1,10) ATOM(NATM),(X(I,NATM),I=1,3)
C	WRITE(6,*) NATM,ATOM(NATM),(X(I,NATM),I=1,3)
10	FORMAT(A4,5X,3F10.5)
      IF (ATOM(NATM).EQ.'END '.OR. ATOM(NATM).EQ.'    ') GOTO 20
      READ(1,*)
      GOTO 5
20	NATM = NATM -1
      CLOSE (1)
      RETURN
      END
C
C *  ROH - MATRIX
C
c
      SUBROUTINE HUBER(ROT,E)

      DIMENSION E(3,3),ROT(3)
c	EQUIVALENCE  (PS,AL,PA), (TH,BE,PB), (PH,GA,PC)
c	EQUIVALENCE  (S1,SIPS),  (S2,SITH),  (S3,SIPH)
c	1   ,(C1,COPS),  (C2,COTH),  (C3,COPH)

      PS=ROT(1)
      TH=ROT(2)
      PH=ROT(3)

      ARC=3.14159265359/180.0
      SIPS = SIN(PS*ARC)
      COPS = COS(PS*ARC)
      SITH = SIN(TH*ARC)
      COTH = COS(TH*ARC)
      SIPH = SIN(PH*ARC)
      COPH = COS(PH*ARC)
      E(1,1) = COPS*COPH - SIPS*SIPH*SITH
      E(1,2) =-SIPS*COTH
      E(1,3) = COPS*SIPH + SIPS*SITH*COPH
      E(2,1) = SIPS*COPH + COPS*SITH*SIPH
      E(2,2) = COPS*COTH
      E(2,3) = SIPS*SIPH - COPS*SITH*COPH
      E(3,1) =-COTH*SIPH
      E(3,2) = SITH
      E(3,3) = COTH*COPH
      END
C
C
C *  ROH - MATRIX
C
      SUBROUTINE HUBERARC(ROT,E)

      DIMENSION E(3,3),ROT(3)
      EQUIVALENCE  (PS,AL,PA), (TH,BE,PB), (PH,GA,PC)
      EQUIVALENCE  (S1,SIPS),  (S2,SITH),  (S3,SIPH)
     1   ,(C1,COPS),  (C2,COTH),  (C3,COPH)

      PS=ROT(1)
      TH=ROT(2)
      PH=ROT(3)

C	ARC=3.14159265359/180.0
C	ARC= 1.0
      SIPS = SIN(PS)
      COPS = COS(PS)
      SITH = SIN(TH)
      COTH = COS(TH)
      SIPH = SIN(PH)
      COPH = COS(PH)
      E(1,1) = COPS*COPH - SIPS*SIPH*SITH
      E(1,2) =-SIPS*COTH
      E(1,3) = COPS*SIPH + SIPS*SITH*COPH
      E(2,1) = SIPS*COPH + COPS*SITH*SIPH
      E(2,2) = COPS*COTH
      E(2,3) = SIPS*SIPH - COPS*SITH*COPH
      E(3,1) =-COTH*SIPH
      E(3,2) = SITH
      E(3,3) = COTH*COPH
      END
c
c
      SUBROUTINE IARRGIVE(N,XYZ0,XYZ)
C    Copy an N-dimensional integer array XYZ0 to another one XYZ.
      INTEGER*4 XYZ(N),XYZ0(N)
      DO I = 1, N
       XYZ(I) = XYZ0(I)
      END DO
      END
c
c
C------MULTIPLY AN INTEGER ARRAY BY A CONSTANT
      SUBROUTINE IARRMC(IM,IN,A1,C,A)
      INTEGER A1(IM,IN),A(IM,IN),C
      DO 10 I1=1,IM
      DO 10 I2=1,IN
10	A(I1,I2)=A1(I1,I2)*C
      END
c
c
      subroutine imatext(isym0,sout,ss1,ss2,ss3,trans)
c A subroutine to transfer a symmetric matrix into a string.
c isym0 is 3*3 symmetric matrix.
c in is 3-dimensional strings of the input like 'x''y''z' or 'h''k''l'
c trans is true. recognize the matrix as a transpose one.
c sout is a 3-dimensional output string
c
c Guoguang 91-02-06
      character sout(3)*12,ss1*1,ss2*1,ss3*1
      logical*1 trans
      integer isym0(3,3)
      integer isym(3,3)
      character out(3,3)*4,temp*4,in(3)*1
c	equivalence (sout,out)
c
      in(1) = ss1
      in(2) = ss2
      in(3) = ss3
C
      if (trans) then
       do j = 1,3
       do i = 1,3
        isym(i,j) = isym0(j,i)
       end do
       end do
      else 
       do j = 1, 3
       do i = 1, 3
        isym(i,j) = isym0(i,j)
       end do
       end do
      end if
c
      do i = 1,3
      do j = 1,3
       out(i,j)(4:4) = in(j)
       if (isym(i,j).eq.0) then
        out(i,j)(1:4) = '    '
       else if(abs(isym(i,j)).ge.10) then
        write(6,*) 'isym(',i,j,')=',isym(i,j),' > 10'
        stop 'error'
       else
        if (isym(i,j).gt.0) out(i,j)(1:1) = '+'
        if (isym(i,j).lt.0) out(i,j)(1:1) = '-'
        if (abs(isym(i,j)).ne.1) then
         write(out(i,j)(2:3),'(i1,a1)') abs(isym(i,j)),'*'
        else
         out(i,j)(2:3) = '  '
        end if
       end if
      end do
      end do
c
      do i = 1, 3
       write(sout(i),'(3a4)') (out(i,j),j=1,3)
c	 write(6,*) sout(i)
c	 if (sout(i)(1:1).eq.'+') sout(i)(1:1) = ' '
c make the text shorter
       il = lenstr(sout(i))
       if (il.gt.1) then
        iln = il
        do j = il-1, 1, -1
         if (sout(i)(j:j).eq.' ') then
         lg = iln - j
         sout(i)(j:iln-1) = sout(i)(j+1:iln)
         sout(i)(iln:iln) = ' '
         iln = iln - 1
         end if
        end do
       end if
       if (sout(i)(1:1).eq.'+') sout(i)(1:1) = ' '
      end do
c
      end
c
c
C-------A STANDARD SUBPROGRAM TO MULTIPLY TWO ARRAYS (a new quick version
C --- this is for I4 version
C) BMC 91/02/11/ Guoguang
C	A1(IM1,IN1)*A2(IM2,IN2)=RSLT(IM1,IN2)
      SUBROUTINE IMATMULT(IM1,IN1,IM2,IN2,A1,A2,RSLT)
      INTEGER A1(IM1,IN1),A2(IM2,IN2),RSLT(IM1,IN2)
      IF (IN1.NE.IM2) THEN
      stop 'The two arrays cannot be multiplied'
      end if
      DO I = 1, IM1
      DO J = 1, IN2
      RSLT(I,J) = 0.
      DO k = 1, in1
      RSLT(I,J) = RSLT(I,J) + A1(I,K)*A2(K,J)
      END DO
      END DO
      END DO
      RETURN
      END
c
c
c if you have a rigid movement x2=R*x1+1
c where R(3,3) is the rotation and t(3) is the translation
c then x1 =R2*x1-t2
c where R2=R-1
c       t2 =R-1*t
c The program reads in R,T matrix (3,4) and returns a inversed RT2
c nmat is the number of matrix
c
      subroutine invrt(nmat,rt,rt2)
      real rt(3,4,nmat),rt2(3,4,nmat)
      real b1(3),b2(3)
      integer me(3)
c
      do i = 1, nmat
       call arrgive(9,rt(1,1,i),rt2(1,1,i))
       call ivsn(3,rt2(1,1,i),b1,b2,me,de,1.e-4)
       call arrmc(3,1,rt(1,4,i),-1.,b1)
       call matmult(3,3,3,1,rt2(1,1,i),b1,rt2(1,4,i))
      end do
c
      end
c
c
C--------- IN PLACE MATRIX INVERSION ---------

      SUBROUTINE IVSN(N,A,B,C,ME,DE,EP)
      
C----  N is the dimension of the matrix A with N rows and N columns
c----  A-N,N- is the input matrix with N rows and N columns.
c---   ME(N), B(N) AND C(N) are the buffer of N-dimensional array.
c---   DE could any value
c      EP is the least value of the pivot i.e. when |A| < EP, the program
c      will stop.

c  on output, DE is the determinant of the original matrix
c  on output, A is the inverse of the original input A
c  ME is an index array which keeps track of row swaps

      DIMENSION A(N,N),ME(N),B(N),C(N)
      INTEGER I,J,K,ME
      
      DE = 1.0
      
      DO 10 J=1,N
10	    ME(J)=J

      DO 20 I=1,N
          Y=0.
          DO 30 J=I,N
              IF (ABS(A(I,J)).LE.ABS(Y)) GOTO 30
              K=J
              Y=A(I,J)
30	    CONTINUE
          DE = DE*Y
          IF (ABS(Y).LT.EP) then
              write(6,*) 'Ill conditioned matrix:'
              write(6,'(3f12.6)') a
              WRITE(6,*)  'the pivot of the matrix is ', y, ' N = ',N
              STOP 4444
          end if
          Y = 1.0/Y
          DO 40 J=1,N
              C(J) = A(J,K)
              A(J,K) = A(J,I)
              A(J,I) = -C(J)*Y
              B(J) = A(I,J)*Y
40	        A(I,J) = A(I,J)*Y
          A(I,I) = Y
          J = ME(I)
          ME(I) = ME(K)
          ME(K) = J
          DO 11 K=1,N
              IF(K.EQ.I) GOTO 11
              DO 12 J=1,N
                  IF (J.EQ.I) GOTO 12
                  A(K,J) = A(K,J) - B(J)*C(K)
12	        CONTINUE
11	    CONTINUE

20	CONTINUE

      DO 33 I=1,N
          DO 44 K=1,N
              IF(ME(K).EQ.I) GOTO 55
44	    CONTINUE
55	    IF(K.EQ.I) GOTO 33
          DO 66 J=1,N
              W = A(I,J)        ! swap rows k and i
              A(I,J) = A(K,J)
66	        A(K,J) = W
          IW = ME(I)            ! keep track of index changes 
          ME(I) = ME(K)
          ME(K) = IW
          DE = -DE  ! determinant changes sign when two rows swapped
33	CONTINUE
      RETURN
      END
c
c
      SUBROUTINE KABMOD(TH1,TH2,TH3,TH11,TH12,TH13)
      DIMENSION TH(6)
C
C --- ERRECHNET EULERWINKEL IM BEREICH -180<TH<180
C
      TH(1) = TH1
      TH(2) = TH2
      TH(3) = TH3
      TH(4) = TH11
      TH(5) = TH12
      TH(6) = TH13
      DO 1 I=1,6
      IF (TH(I).LE.-180.) TH(I) = TH(I) + 360.
    1 IF (TH(I).GT.180.) TH(I) = TH(I) - 360.
      TH1 = TH(1)
      TH2 = TH(2)
      TH3 = TH(3)
      TH11 = TH(4)
      TH12 = TH(5)
      TH13 = TH(6)
      RETURN
      END
c
c
      SUBROUTINE LATTIC(NSYM,SYM,LAT)
C A subroutine to add the translation vector to the symmetry arrays 
c and modify the number of symmetry operations according to the
c Bravais Lattice.
C	CHARACTER*1 LAT
C	DIMENSION SYM(3,4,192)
c Parameter description: 
c NSYM are the number of the symmetric operation. The input value do
c   not include the non-1 translation operation. But output includes it.
c SYM(3,4,NSYM) is the symmetry op arrays. 
c LAT is the Bravais lattice. It could be "P", "A", "B", "C", "I",
c "F",or "R"
C                                  Written by      Guoguang Lu
c                                                 27/02/88
      CHARACTER*1 LAT
      DIMENSION SYM(3,4,192),TR(3,4)
c
      IF (LAT.EQ.'P') RETURN
      CALL ARRMC(3,4,TR,0.,TR)
      IF (LAT.EQ.'A') THEN
      NTR = 2
      TR(2,2) = .5
      TR(3,2) = .5
      ELSE IF (LAT.EQ.'B') THEN
      NTR = 2
      TR(1,2) = .5
      TR(3,2) = .5
      ELSE IF (LAT.EQ.'C') THEN
      NTR = 2
      TR(2,2) = .5
      TR(1,2) = .5
      ELSE IF (LAT.EQ.'I') THEN
      NTR = 2
      TR(1,2) = .5
      TR(2,2) = .5
      TR(3,2) = .5
      ELSE IF (LAT.EQ.'F') THEN
      NTR = 4
      TR(1,2) = .5
      TR(2,2) = .5
      TR(2,3) = .5
      TR(3,3) = .5
      TR(3,4) = .5
      TR(1,4) = .5
      ELSE IF (LAT.EQ.'R') THEN
      NTR = 3
      TR(1,2) = 1./3.
      TR(2,2) = 2./3.
      TR(3,2) = 2./3.
      TR(1,3) = 2./3.
      TR(2,3) = 1./3.
      TR(3,3) = 1./3.
      ELSE 
      WRITE(6,*)  'Error Lattice>> No lattice ',LAT
      WRITE(6,*)  'Only P,A,B,C,I,F and R are allowed.'
      END IF
      DO ISYM = 1, NSYM
      DO ITR = 2, NTR
      I = (ITR-1)*NSYM + ISYM
      CALL ARRMC(3,3, SYM(1,1,ISYM), 1., SYM(1,1,I) )
      CALL ARRAD(3,1, TR(1,ITR), SYM(1,4,ISYM), SYM(1,4,I) )
      DO J = 1, 3
5	IF (SYM(J,4,I).GE.1.) THEN
      SYM(J,4,I) = SYM(J,4,I) - 1.
      GOTO 5
      END IF
6	IF (SYM(J,4,I).LT.-1.) THEN
      SYM(J,4,I) = SYM(J,4,I) + 1.
      GOTO 6
      END IF
      END DO		! J
      END DO		! ITR
      END DO		! ISYM
      NSYM = NSYM * NTR
      RETURN	
      END
c
c
      SUBROUTINE LSQEQ(M,N,A,PHI,DETA,ATA,B1)
C	DIMENSION A(M,N),PHI(M),DETA(N),ATA(N,N)
C	DIMENSION B1(N)
C This is a subroutine to compute the solution of least square equation
C of a normal matrix. i.e. AT*A*DETA = AT*PHI ==> DETA = (AT*A)^(-1)*AT*PHI
C It could be get the deta the make the SIGMA (phi)^2 = minimum
c        where  there are M equations
c                   .......
c                  phi i = 0
c                   .......
c
c  Note:  In the program parameter PHI(I) = - phi0 i  in the equation.
c
C                                             written by Guoguang Lu
c                                               23/01/1989
c parameter description:
c M is the number of equations. i.e. PHIi = 0 (i=1,...M)
C N is the number of the variables.
C           (N should be less than 3*N0 given in parameter statement.)
C           (At this moment, N0 is 6000.)
c (A,PHI) is a M*(N+1) normal  matrix
c A is a M*N array
c PHI is a M-dimensional vector. Here, PHI should be -phi0 in the equation.
c DETA is the output shift.
c ATA(N,N ) is a N*N buffer array.
c B1(N)     is a N-dimensional buffer array.
c Note :   M must not be less than N.
      PARAMETER (N0 = 6000 )
      DIMENSION A(M,N),PHI(M),DETA(N),ATA(N,N)
      DIMENSION B1(N),ME(3*N0)
      IF (M.LT.N) THEN
      STOP 'Equation number is not enough'
      else if (M.EQ.N) THEN
      CALL IVSN(M,A,B1,DETA,ME,DE,1.E-5)
c	CALL ARRMC(N,1,DETA,0.,DETA)
      call arrvalue(n,deta,0.)
      DO 50 I = 1, N
      DO 50 J = 1, N
50	DETA(I) = DETA(I) + ATA(I,J)*PHI(J)
      ELSE
C
C        T
C ATA = A * A
C
      call arrvalue(n*n,ata,0.)
c	CALL ARRMC(N,N,ATA,0.,ATA)
      DO 10 I = 1,N
      DO 10 J = I,N
      DO 10 K = 1,M
10	ATA(I,J) = ATA(I,J) + A(K,I)*A(K,J)
      IF (N.GE.2) THEN
      DO 20 I = 1, N
      DO 20 J = 1, I-1
20	ATA(I,J) = ATA(J,I)
      END IF
C
C
C ATA = ATA^-1
C
      CALL IVSN (N,ATA,B1,DETA,ME,DE,1.E-5)
C
C       T
C B1 = A * PHI
C
c	CALL ARRMC(N,1,B1,0.,B1)
      call arrvalue(n,b1,0.)
      DO 30 I = 1, N
30	B1(I) = POIMULT(M, M, A(1,I), PHI)
C
C  DETA = ATA^-1 * B1
C
C                T        T
C i.e. DETA  = (A * A) * A * PHI
C
c	CALL ARRMC(N,1,DETA,0.,DETA)
      call arrvalue(n,deta,0.)
      DO 40 I = 1, N
      DO 40 J = 1, N
40	DETA(I) = DETA(I) + ATA(I,J)*B1(J)
C
      END IF
C
      END
c
c
      subroutine matitor(im,in,imat,rmat)
c    Convert a matrix from integer*2 to real*4
c
      real*4 rmat(im,in)
      integer imat(im,in)
c
      do j = 1, in
      do i = 1, im
      rmat(i,j) = imat(i,j)
      end do
      end do
      end
c
c
C-------A STANDARD SUBPROGRAM TO MULTIPLY TWO ARRAYS (a new quick version
C) BMC 89/11/05/ Guoguang
C	A1(IM1,IN1)*A2(IM2,IN2)=RSLT(IM1,IN2)
      SUBROUTINE MATMULT(IM1,IN1,IM2,IN2,A1,A2,RSLT)
      DIMENSION A1(IM1,IN1),A2(IM2,IN2),RSLT(IM1,IN2)
      IF (IN1.NE.IM2) THEN
      stop 'The two arrays cannot be multiplied'
      end if
      DO I = 1, IM1
      DO J = 1, IN2
      RSLT(I,J) = 0.
      DO k = 1, in1
      RSLT(I,J) = RSLT(I,J) + A1(I,K)*A2(K,J)
      END DO
      END DO
      END DO
      RETURN
      END
c
c
      subroutine matrtoi(im,in,rmat,imat)
c    Convert a matrix from real to integer
c
      real*4 rmat(im,in)
      integer imat(im,in)
c
      do j = 1, in
      do i = 1, im
      imat(i,j) = rmat(i,j)
      end do
      end do
      end
c
c
      FUNCTION MATSYM(S,TEXT,ICOL) 
C 
C A FORTRAN FUNCTION SUBPROGRAM TO SCAN A SYMMETRY CARD AND BUILD 
C THE CORRESPONDING SYMMETRY-OPERATION MATRIX.  THE CONTENTS OF THE 
C CARD AND ERRORS ENCOUNTERED ARE WRITTEN OFF-LINE IN THE CALLING 
C PROGRAM.  THE FUNCTIONAL VALUE IS NORMALLY ZERO, WHEN AN ERROR HAS
C BEEN DETECTED THE NUMBERS 1-4 ARE GIVEN BACK: 
C ERROR-NUMBER 1:  INVALID CHARACTER AT POSITION .... 
C              2:  SYNTAX ERROR AT NONBLANK POSITION .... 
C              3:  BAD COMMAS 
C              4:  DETERMINANT OF ROTATION MATRIX NOT +1 OR -1
C 
C EXPLANATION OF ARGUMENTS
C S IS THE 1ST ELEMENT OF THE 3X4 ARRAY WHERE THE SYM. OP. IS TO BE 
C STORED. 
C TEXT  IS THE ARRAY OF THE TEXT RECORD WHICH IS TO BE SCANNED  (LENGTH:
C 36 CHARACTERS, CONTIGUOUSLY)
C ICOL POINTS - IN CASE OF ERRORS - TO THE POSITION OF AN INVALID 
C CHARACTER OR AN SYNTAX ERROR WITHIN THE SYMMETRY CARD 
C 
C     MODIFIED FOR PDP-11 06/JAN/78 S.J. REMINGTON
C ***** LAST UPDATE 07/10/76   MRS                WRS:PR2.MATSYM
C ***** PREVIOUS UPDATES                              07/07/75  06/11/72
C 
C 
      DIMENSION S(3,4), ICHAR(36)
      CHARACTER*1 CHAR,VALID(14),BLANK
      CHARACTER*36 TEXT
      EQUIVALENCE (BLANK,VALID(14))
      DATA VALID/'1','2','3','4','5','6','X','Y','Z','-','+','/',',',
     . ' '/
C 
C----INITIALIZE SIGNALS 
      ICOL = 0
      IFIELD = 0
C----CLEAR SYMOP ARRAY
      DO 4 J=1,4
      DO 4 I=1,3
 4    S(I,J) = 0.0
C 
C----TEST THE SYMMETRY TEXT (TEXT) FOR ILLEGAL CHARACTERS
C 
C      THE STATEMENT IS SCANNED AND ARRAY ICHAR IS FILLED WITH THE INDICES
C      OF EACH ELEMENT OF TEXT INTO THE ARRAY 'VALID' (I.E. IF THE SEVENTH
C      CHAR OF 'TEXT' -IGNORING BLANKS- IS 'Y', ICHAR(7)=8)
C
      ICOLMX = 0
      DO 20 I=1,36
      CHAR = TEXT(I:I)
      IF(CHAR .EQ. BLANK) GO TO 20
      ICOLMX = ICOLMX + 1
      DO 10 J=1,13
      IF( VALID(J) .NE. CHAR) GO TO 10
      ICHAR(ICOLMX) = J
      GO TO 20
   10 CONTINUE
C--IF GOT TO HERE, THE CHAR ISN'T IN 'VALID'
      JTXTSC = I
      GO TO 9975
   20 CONTINUE
C
C----BEGIN FIELD LOOP 
C 
 101  IFIELD = IFIELD + 1 
      IF(IFIELD-3) 104,104,103
C  HAVE ALL CHARACTERS  BEEN ANALYSED IN 3 FIELDS?
 103  IF(ICOL-ICOLMX) 9987,1000,1000
C---NO MORE THAN THREE FIELDS PERMITTED 
 104  SIGN = 1.0
      ICOL = ICOL+1 
C----MARCH ON, FIND WHAT'S NEXT
      IFNUM=ICHAR(ICOL) 
      IF (11-IFNUM) 9980,110,110
C----TEST FOR SIGNS 
 110  IF (10-IFNUM) 118,115,112 
C----TEST TO DISTINGUISH BETWEEN NUMBERS AND LETTERS
 112  IF (6-IFNUM) 125,130,130
C----FOUND A SIGN 
C     COME HERE FOR MINUS 
 115  SIGN = -1.0 
C     COME HERE FOR PLUS
 118  ICOL = ICOL + 1 
C----NEXT CHARACTER MUST BE A NUMBER OR LETTER
      IFNUM=ICHAR(ICOL) 
      IF (IFNUM-9) 122,122,9980 
 122  IF (IFNUM-6) 130,130,125
C----A LETTER 
 125  IFNUM = IFNUM - 6 
      S(IFIELD,IFNUM) = SIGN
      GO TO 150 
C----A NUMBER, FLOAT IT INTO DIGIT
 130  DIGIT = IFNUM 
      ICOL = ICOL + 1 
C----IS NEXT CHARACTER A SLASH
      IF(ICHAR(ICOL).EQ.12)  GO TO 135
C----NO SLASH, TRANSLATION TERM COMPLETE
      ICOL=ICOL-1 
      GO TO  140
C----A SLASH IS NEXT, IS THERE A NUMBER AFTER IT
 135  ICOL = ICOL + 1 
      IFNUM=ICHAR(ICOL) 
      IF (IFNUM-6) 138,138,9980 
C----MAKE FRACTION;   '5' NOT ALLOWED IN DENOMINATOR
  138 IF(IFNUM.EQ.5) GO TO 9980 
      DIGIT = DIGIT/FLOAT(IFNUM)
C----ACCUMULATE THE VECTOR COMPONENT
 140  S(IFIELD,4) = S(IFIELD,4) + SIGN*DIGIT
C----TEST TO SEE IF THERE ARE MORE CHARACTERS IN THIS SUBFIELD
C    OR A NEW SUBFIELD
 150  ICOL = ICOL + 1 
      SIGN = 1.0
C----TEST IF ALL DONE 
      IF(ICOL - ICOLMX) 151,151,1000
 151  CONTINUE
C----A SUBFIELD BEGINS WITH A PLUS OR MINUS 
      IF(ICHAR(ICOL).EQ.11)  GO TO 118
      IF(ICHAR(ICOL).EQ.10)  GO TO 115
C----A COMMA MUST BE NEXT UNLESS ICOL IS ICOLMX+1 WHICH MEANS THE END 
      IF(ICHAR(ICOL).EQ.13)   GO TO 101 
 163  IF (ICOLMX+1-ICOL) 9980,1000,9980 
C----EVERYTHING SEEMS OK SEE IF DETERMINANT IS A NICE + OR - 1. 
 1000 CONTINUE
      D=S(1,1)*S(2,2)*S(3,3)
     $ -S(1,2)*S(2,1)*S(3,3)
     $ -S(1,1)*S(2,3)*S(3,2)
     $ -S(1,3)*S(2,2)*S(3,1)
     $ +S(1,2)*S(2,3)*S(3,1)
     $ +S(1,3)*S(2,1)*S(3,2)
      IF(ABS(ABS(D) - 1.0) -.0001) 1001,9985,9985 
 1001 CONTINUE
      IF(IFIELD - 3) 9987,1002,9987 
 1002 CONTINUE
C----LEGAL RETURN 
      MATSYM = 0
 2005 RETURN
C 
C 
C----ILLEGAL CHARACTER ON SYMMETRY CARD 
 9975 MATSYM = 1
      ICOL = JTXTSC 
      GO TO 2005
C 
C----SYNTAX ERROR AT NONBLANK CHARACTER ICOL
 9980 MATSYM = 2
      GO TO 2005
C 
C----DETERMINANT IS NOT + OR - 1. 
 9985 MATSYM = 4
      GO TO 2005
C 
C----INCORRECT NUMBER OF COMMAS 
 9987 MATSYM = 3
      GO TO 2005
      END 
c
c
c this is subroutine to transfer missetting angle to matrix
c  ROT = (phiz)*(phiy)*(phix)
      subroutine misseting(phi,rot)
      real phi(3),rot(3,3)
      external cosd, sind
      sinx = sind(phi(1))
      cosx = cosd(phi(1))
      siny = sind(phi(2))
      cosy = cosd(phi(2))
      sinz = sind(phi(3))
      cosz = cosd(phi(3))
      rot(1,1) = cosz*cosy
      rot(2,1) = sinz*cosy
      rot(3,1) = -siny
      rot(1,2) = cosz*siny*sinx - sinz*cosx
      rot(2,2) = sinz*siny*sinx + cosz*cosx
      rot(3,2) = cosy*sinx
      rot(1,3) = cosz*siny*cosx + sinz*sinx
      rot(2,3) = sinz*siny*cosx - cosz*sinx
      rot(3,3) = cosy*cosx
      end	
c if you have a rigid movement x2=R*x1+1
c where R(3,3) is the rotation and t(3) is the translation
c then x1 =R2*x1-t2
c where R2=R-1
c       t2 =R-1*t
c The program reads in R,T matrix (3,4) 1 and 2, then returns an RT
c which is rt = rt1*rt2
c
c 
c
      subroutine morert(rt1,rt2,rt)
      real rt(3,4),rt2(3,4),rt1(3,4)
      real gt(4,4),gt2(4,4),gt1(4,4)
      real extr(4)
      real b1(3),b2(3)
      integer me(3)
c...  data statements.  Separate declaration and init required for f2c
      data extr /0.,0.,0.,1./
c
      do i = 1, 4
       gt1(i,4) = extr(i)
       gt2(i,4) = extr(i)
      end do
c
      do j = 1, 3
       do i = 1, 3
        gt1(i,j) = rt1(i,j)
        gt2(i,j) = rt2(i,j)
       end do
      end do
c
      do j = 1, 3
       gt1(4,j) = rt1(j,4)
       gt2(4,j) = rt2(j,4)
      end do
c
      call matmult(4,4,4,4,gt1,gt2,gt)
c
      do j = 1, 3
       do i = 1, 3
        rt(i,j) = rt(i,j)
       end do
      end do
c
      do j = 1, 4
       rt(j,4) = gt(4,j)
      end do
c
c	WRITE(6,*) 'gt'
c	write(6,'(4f10.5)') gt1
c	write(6,'(4f10.5)') gt2
c	write(6,'(4f10.5)') gt
c
      end
c
c
c a subroutine to transfer a normalized matrix to denzo missetting angle 
c definition [rot]=[rotx]*[roty]*[rotz]
c  GuoGuang 930706
      subroutine mtodenmis(rot0,phi)
      real rot0(3,3)
      real rot(3,3),phi(3)
      external asind, atand, cosd
c 	sinx = sind(phi(1))
c 	cosx = cosd(phi(1))
c 	siny = sind(phi(2))
c 	cosy = cosd(phi(2))
c 	sinz = sind(phi(3))
c 	cosz = cosd(phi(3))
c	rot0(1,1) =  cosy*cosz
c	rot0(2,1) = -cosx*sinz+sinx*siny*cosz
c	rot0(3,1) =  sinx*sinz+cosx*siny*cosz
c	rot0(1,2) =  cosy*sinz
c	rot0(2,2) =  cosx*cosz+sinx*siny*sinz
c	rot0(3,2) = -sinx*cosz+cosx*siny*sinz
c	rot0(1,3) = -siny
c	rot0(2,3) =  sinx*cosy
c	rot0(3,3) =  cosx*cosy
c
      call antiarr(3,3,rot0,rot)
c 	sinx = sind(phi(1))
c 	cosx = cosd(phi(1))
c 	siny = sind(phi(2))
c 	cosy = cosd(phi(2))
c 	sinz = sind(phi(3))
c 	cosz = cosd(phi(3))
c	rot(1,1) = cosz*cosy
c	rot(2,1) = sinz*cosy
c	rot(3,1) = -siny
c	rot(1,2) = cosz*siny*sinx - sinz*cosx
c	rot(2,2) = sinz*siny*sinx + cosz*cosx
c	rot(3,2) = cosy*sinx
c	rot(1,3) = cosz*siny*cosx + sinz*sinx
c	rot(2,3) = sinz*siny*cosx - cosz*sinx
c	rot(3,3) = cosy*cosx
      phi(2) = asind(-rot(3,1))
      if (abs(phi(2)).ne.90.) then
       if (rot(1,1).eq.0.) then
        phi(3) = 90.
       else
        phi(3) = atand(rot(2,1)/rot(1,1))
        if (phi(3).gt.0..and.rot(1,1)/cosd(phi(2)).lt.0.)
     1  phi(3) = 180. + phi(3)
        if (phi(3).lt.0..and.rot(2,1)/cosd(phi(2)).gt.0.) 
     1 phi(3) = 180. + phi(3)
       end if
       if (rot(3,3).eq.0.) then
        phi(1) = 90.
       else
        phi(1) = atand(rot(3,2)/rot(3,3))
        if (phi(1).gt.0..and.rot(3,3)/cosd(phi(2)).lt.0.)
     1  phi(1) = 180. + phi(1)
        if (phi(1).lt.0..and.rot(3,2)/cosd(phi(2)).gt.0.) 
     1 phi(1) = 180. + phi(1)
       end if
      else
       WRITE(6,*)  'not implemented yet'
      end if
      end 
c
c
      SUBROUTINE MTOHUBER(E,HB,HB1)
      DIMENSION E(3,3),HB(3),HB1(3)
C
C * ROH ANGLES:  -180 < PSI < 180
C                 -90 < TH  <  90
C                -180 < PHI < 180
C
      ARC=3.14159265359/180.0
      if (e(3,2).gt.1.0) e(3,2) = 1.0
      if (e(3,2).lt.-1.0) e(3,2) = -1.0
210	TH = ASIN(E(3,2))
      COSTH = COS(TH)
c	IF (COSTH.GE.0.01) GO TO 212
      if (costh.lt.1.0e-6) then
       ph = 0.
       ps = acos(e(1,1))
       if (e(2,1).lt.0) ps = -ps
       goto 1500
      end if
c	WRITE (6,1212) COSTH
c1212	FORMAT('-*** WARNING :  COS(THETA) = ', E8.2, ' ***')
212	PH = E(3,3)/COSTH
C	IF (PH.GT.1.0.OR.PH.LT.-1.0) WRITE(6,1291) PH
c	IF (ABS(PH).GT.(1.00001)) WRITE(6,1291) PH
1291	FORMAT(' ARCOS(E33/COSTH)=',F12.5)
1020	FORMAT(' *** ERROR IN ROH *** ',2F10.4/)
      IF (PH.GT.1.0) PH=1.0
      IF (PH.LT.-1.0) PH=-1.0
      PH=ACOS(PH)
      TST=-COSTH*SIN(PH)
      IF (TST*E(3,1).LT.0.0) PH=-PH
      PS = E(2,2)/COSTH
c	IF (PS.GT.1.0.OR.PS.LT.-1.0) WRITE(6,1292) PS
c	if (abs(ps).gt.1.0001) write(6,1292) ps
1292 	FORMAT(' ARCOS(E22/COSTH)=',F12.5)
      IF (PS.GT.1.0) PS=1.0
      IF (PS.LT.-1.0) PS=-1.0
      PS=ACOS(PS)
      TST=-COSTH*SIN(PS)
      IF (TST*E(1,2).LT.0.0) PS=-PS
      TST = -SIN(PS)*SIN(PH)*SIN(TH) + COS(PS)*COS(PH)
c	IF (ABS(TST-E(1,1)) .GT. 0.1) WRITE (6,1020) TST,E(1,1)
      TST = -COS(PS)*SIN(TH)*COS(PH) + SIN(PS)*SIN(PH)
c	IF (ABS(TST-E(2,3)) .GT. 0.1) WRITE (6,1020) TST,E(2,3)
1500	continue
      PS = PS/ARC
      TH = TH/ARC
      PH = PH/ARC
      PS1 = 180. + PS
      TH1 = 180. - TH
      PH1 = 180. + PH
      CALL KABMOD(PS,TH,PH,PS1,TH1,PH1)
      HB(1)=PS
      HB(2)=TH
      HB(3)=PH
      HB1(1)=PS1
      HB1(2)=TH1
      HB1(3)=PH1
      END
c
c
      SUBROUTINE MTOHUBERARC(E,HB,HB1)
      DIMENSION E(3,3),HB(3),HB1(3)
C
C * ROH ANGLES:  -180 < PSI < 180
C                 -90 < TH  <  90
C                -180 < PHI < 180
C
C	ARC=3.14159265359/180.0
C	ARC=1.0
210	TH = ASIN(E(3,2))
      COSTH = COS(TH)
      IF (COSTH.GE.0.01) GO TO 212
      WRITE (6,1212) COSTH
1212	FORMAT('-*** WARNING :  COS(THETA) = ', E8.2, ' ***')
212	PH = E(3,3)/COSTH
c	IF (PH.GT.1.0.OR.PH.LT.-1.0) WRITE(6,1291) PH
1291	FORMAT(' ARCOS(E33/COSTH)=',F12.5)
1020	FORMAT(' *** ERROR IN ROH *** ',2F10.4/)
      IF (PH.GT.1.0) PH=1.0
      IF (PH.LT.-1.0) PH=-1.0
      PH=ACOS(PH)
      TST=-COSTH*SIN(PH)
      IF (TST*E(3,1).LT.0.0) PH=-PH
      PS = E(2,2)/COSTH
c	IF (PS.GT.1.0.OR.PS.LT.-1.0) WRITE(6,1292) PS
1292 	FORMAT(' ARCOS(E22/COSTH)=',F12.5)
      IF (PS.GT.1.0) PS=1.0
      IF (PS.LT.-1.0) PS=-1.0
      PS=ACOS(PS)
      TST=-COSTH*SIN(PS)
      IF (TST*E(1,2).LT.0.0) PS=-PS
      TST = -SIN(PS)*SIN(PH)*SIN(TH) + COS(PS)*COS(PH)
c	IF (ABS(TST-E(1,1)) .GT. 0.1) WRITE (6,1020) TST,E(1,1)
      TST = -COS(PS)*SIN(TH)*COS(PH) + SIN(PS)*SIN(PH)
c	IF (ABS(TST-E(2,3)) .GT. 0.1) WRITE (6,1020) TST,E(2,3)
C	PS = PS/ARC
C	TH = TH/ARC
C	PH = PH/ARC
      PS1 = 180. + PS
      TH1 = 180. - TH
      PH1 = 180. + PH
      CALL KABMOD(PS,TH,PH,PS1,TH1,PH1)
      HB(1)=PS
      HB(2)=TH
      HB(3)=PH
      HB1(1)=PS1
      HB1(2)=TH1
      HB1(3)=PH1
      END
c
c
c this is subroutine to transfer missetting angle to matrix
c  ROT = (phiz)*(phiy)*(phix)
      subroutine mtomisset(rot,phi)
      real phi(3),rot(3,3)
      external asind, atand, cosd
c 	sinx = sind(phi(1))
c 	cosx = cosd(phi(1))
c 	siny = sind(phi(2))
c 	cosy = cosd(phi(2))
c 	sinz = sind(phi(3))
c 	cosz = cosd(phi(3))
c	rot(1,1) = cosz*cosy
c	rot(2,1) = sinz*cosy
c	rot(3,1) = -siny
c	rot(1,2) = cosz*siny*sinx - sinz*cosx
c	rot(2,2) = sinz*siny*sinx + cosz*cosx
c	rot(3,2) = cosy*sinx
c	rot(1,3) = cosz*siny*cosx + sinz*sinx
c	rot(2,3) = sinz*siny*cosx - cosz*sinx
c	rot(3,3) = cosy*cosx
      phi(2) = asind(-rot(3,1))
      if (abs(phi(2)).ne.90.) then
       if (rot(1,1).eq.0.) then
        phi(3) = 90.
       else
        phi(3) = atand(rot(2,1)/rot(1,1))
        if (phi(3).gt.0..and.rot(1,1)/cosd(phi(2)).lt.0.)
     1  phi(3) = 180. + phi(3)
        if (phi(3).lt.0..and.rot(2,1)/cosd(phi(2)).gt.0.) 
     1 phi(3) = 180. + phi(3)
       end if
       if (rot(3,3).eq.0.) then
        phi(1) = 90.
       else
        phi(1) = atand(rot(3,2)/rot(3,3))
        if (phi(1).gt.0..and.rot(3,3)/cosd(phi(2)).lt.0.)
     1  phi(1) = 180. + phi(1)
        if (phi(1).lt.0..and.rot(3,2)/cosd(phi(2)).gt.0.) 
     1 phi(1) = 180. + phi(1)
       end if
      else
       WRITE(6,*)  'not implemented yet'
      end if
      end 
c
C
C * POLAR ANGLES
C
      SUBROUTINE MTOPOLOR(E,POL,POL1)
      DIMENSION E(3,3), POL(3),POL1(3)
      ARC=3.14159265359/180.0
230	COSPC = (E(1,1)+E(2,2)+E(3,3)-1.0)/2.0
      IF (COSPC .LT.-1.0)
     1 WRITE(6,1290)
1290	FORMAT(1X,' SUM OF DIAGONAL-ELEMENTS LESS
     1 THAN -1./IS SET TO -1.')
      IF (COSPC .LT.-1.0) COSPC=-1.0
C --- PC == KAPPA (ALIAS CHI)
      PC = ACOS(COSPC)
      ARG = (E(2,2)-COSPC)/(1.0-COSPC)
      IF (ARG .GE. 0.0) GO TO 235
      WRITE(6,1235) ARG
1235	FORMAT('-*** ARG = ',E10.4, ' * ARG=0 ASSUMED ***')
      ARG = 0.0
235	COSPA = SQRT(ARG)
C --- PA == THETA (ALIAS PSI)
      PA = ACOS(COSPA)
      SINPA = SIN(PA)
      TST = 2*COSPA*SIN(PC)
      IF ((E(1,3)-E(3,1))*TST .LT. 0.0) PC=-PC
      SINPC = SIN(PC)
      IF (SINPA.NE.0.0) GO TO 237
      WRITE (6,1237)
1237	FORMAT('-*** SIN(THETA) = 0;  PHI UNDETERMINED ***')
      PB = 0.
      GO TO 250
C --- PB == PHI
237	IF(ABS(SINPC) .LT. 0.0001 .AND. 
     1 ABS(COSPA) .LT. 0.0001) GOTO 245
      IF ( ABS(SINPA*SINPC) .LT. 0.001) GO TO 240
      PB = ACOS( (E(3,2)-E(2,3)) / (2*SINPA*SINPC) )
      TST = 2.0*SINPA*SIN(PB)*SINPC
      IF (TST*(E(1,2)-E(2,1)) .LT. 0.0) PB=-PB
      GO TO 250
240	FAC = 2.0*SINPA*COSPA*(1.0-COSPC)
      PB  = ACOS((E(1,2)+E(2,1))/FAC)
      TST = FAC*SIN(PB)
      IF (-(E(2,3)+E(3,2))*TST .LT. 0.0) PB = -PB
      GO TO 250
245	DENOM = (1.0-ARG)*(COSPC-1.0)
      PB = 0.5*ASIN((E(1,3)+E(3,1))/DENOM)
      SINSPB = (COSPC-E(3,3))/DENOM
      IF (SIN(PB)**2 .LT. SINSPB) PB = 90.*ARC - PB
250	PA = PA/ARC
      PB = PB/ARC
      PC = PC/ARC
      TH1 = 180. - PA
      PH1 = 180. + PB
      PK1 =      - PC
      CALL KABMOD(PA,PB,PC,TH1,PH1,PK1)
      POL(1)=PA
      POL(2)=PB
      POL(3)=PC
      POL1(1)=TH1
      POL1(2)=PH1
      POL1(3)=PK1

      END
c
c
      subroutine mtopolors(e,p1,p2)
      dimension e(3,3),p1(3),p2(3)
      dimension b1(3),b2(3)
      external acosd

      call mtovec(e,b1,p1(3))
C
      if (p1(3).eq.0.) then
      p1(1) = 0.
      p1(2) = 0.
      p2(1) = 0.
      p2(2) = 0.
      p2(3) = 0.
      return
      end if
C
      p1(1) = acosd(b1(2))
      dos = sqrt(b1(1)*b1(1)+b1(3)*b1(3))
      if (dos.eq.0.) then
       if (vem(3,b1).eq.0.and.p1(3).ne.0.) then
       WRITE(6,*)  'vec:',b1,'kapa:',p1(3)
       end if
      p1(2) = 0.
      p2(2) = 0.
      p2(1) = 180. - p1(1)
      if (p2(1).ge.360.) p2(1) = 0.
      return
      end if
c
c
      p1(2) = acosd( b1(1)/dos )
c	b22 = asind( b1(3)/sqrt(b1(1)*b1(1)+b1(3)*b1(3)) )
      if (b1(3).gt.0.) p1(2) = -p1(2)
c
      p2(3) = -p1(3)
      p2(2) = 180. + p1(2)
      if (p2(2).gt.180) p2(2) = p2(2) -360.
      p2(1) = 180. - p1(1)
      end
c
c
      subroutine mtopolorz(e,p1,p2)
c 
      real e(3,3),p1(3),p2(3)
      real b1(3),b2(3)
      external acosd
      call mtovec(e,b1,p1(3))
c
      if (abs(p1(3)).lt.1e-3) then
      p1(1) = 0
      p2(1) = 0
      p1(2) = 0
      p2(2) = 0
      p2(3) = 0
      return
      end if

c	p1(1) = acosd(b1(2))
c	p1(2) = acosd( b1(1)/sqrt(b1(1)*b1(1)+b1(3)*b1(3)) )
c	b22 = asind( b1(3)/sqrt(b1(1)*b1(1)+b1(3)*b1(3)) )
      p1(1) = acosd(b1(3))
c
      if (p1(1).eq.0.or.p1(1).eq.180.) then
      p1(2) = 0.
      p2(2) = 0.
      p2(1) = p1(1) + 180.
      if (p2(1).ge.360) p2(1) = p2(1) - 360.
      p2(3) = -p1(3)
      return
      end if
c
      p1(2) = acosd( b1(1)/sqrt(b1(1)*b1(1)+b1(2)*b1(2)) )
c	b22 = asind( b1(2)/sqrt(b1(1)*b1(1)+b1(3)*b1(3)) )
c	if (b22.gt.0.) p1(2) = -p1(2)
      if (b1(2).lt.0.) p1(2) = -p1(2)
c
      p2(3) = -p1(3)
      p2(1) = 180. - p1(1)
      p2(2) = 180. + p1(2)
      if (p2(2).gt.180) p2(2) = p2(2) -360.
      end
c
c
      SUBROUTINE MTOR_B(IOPT2,E,ROT1,ROT2)
C
C * R&B ANGLES:  -180 < AL < 180
C                   0 < BE < 180
C                -180 < GA < 180
C
      EQUIVALENCE  (PS,AL,PA), (TH,BE,PB), (PH,GA,PC)
      EQUIVALENCE  (S1,SIPS),  (S2,SITH),  (S3,SIPH)
     *            ,(C1,COPS),  (C2,COTH),  (C3,COPH)
      DIMENSION   E(3,3),E1(3,3),E2(3,3),X(3),X1(3)
      DIMENSION ROT1(3),ROT2(3)
      DATA        E2 / 1., 3*0., 1., 3*0., 1./
      ARC=3.14159265359/180.0
  220 BE = ACOS(E(3,3))
      GA = E(3,2)/SIN(BE)
c      IF (GA.GT.1.0.OR.GA.LT.-1.0) WRITE(6,1293) GA
1293  FORMAT(' ARCOS(E32/SINBE)=',F12.5)
      IF (GA.GT.1.0) GA=1.0
      IF (GA.LT.-1.0) GA=-1.0
      GA=ACOS(GA)
      TST= SIN(BE)*SIN(GA)
      IF (TST*E(3,1).LT.0.0) GA=-GA
      AL = -E(2,3)/SIN(BE)
c      IF (AL.GT.1.0.OR.AL.LT.-1.0) WRITE(6,1294) AL
1294  FORMAT(' ARCOS(-E23/SINBE)=',F12.5)
      IF (AL.GT.1.0) AL=1.0
      IF (AL.LT.-1.0) AL=-1.0
      AL=ACOS(AL)
      TST= SIN(BE)*SIN(AL)
      IF (TST*E(1,3).LT.0.0) AL=-AL
      TST = -SIN(AL)*COS(BE)*SIN(GA) + COS(AL)*COS(GA)
      IF (ABS(TST-E(1,1)) .GT. 0.1) WRITE (6,1010) TST,E(1,1)
      TST = COS(AL)*COS(BE)*COS(GA) - SIN(AL)*SIN(GA)
      IF (ABS(TST-E(2,2)) .GT. 0.1) WRITE (6,1010) TST,E(2,2)
 1010 FORMAT(' *** ERROR IN R&B *** ',2F10.4/)
      AL = AL/ARC
      BE = BE/ARC
      GA = GA/ARC
      IF (IOPT2.NE.4) GO TO 225
      AL=AL-90.0
      GA=GA+90.0
  225 CONTINUE
      AL1 = 180. + AL
      BE1 =      - BE
      GA1 = 180. + GA
      CALL KABMOD(AL,BE,GA,AL1,BE1,GA1)
C      IF (IOPT2.EQ.2) WRITE (6,1220) AL,BE,GA,AL1,BE1,GA1
C      IF (IOPT2.EQ.4) WRITE (6,1221) AL,BE,GA,AL1,BE1,GA1
C      GO TO 100
      ROT1(1)=AL
      ROT1(2)=BE
      ROT1(3)=GA
      ROT2(1)=AL1
      ROT2(2)=BE1
      ROT2(3)=GA1
      END
c
c
      subroutine mtovec(a,vec,vkapa)
c this is a subroutine which give a transfer a orthogonal matrix
c to a polar angle rotate a direction of a vector which is 1.
c a is the input 3*3  orthogonal matrix.
c vec is the output 3-dim vector which indicate the axis direction
c kappa is the output polar angle.
c Guoguang 901001 bmc
      dimension a(3,3),vec(3),a1(3,3),a2(9),vec1(3)
      integer*2 sgn(3,8)
      external acosd
      data sgn/ 1,1,1,  1,1,-1,  1,-1,1,  1,-1,-1,
     1 	 -1,1,1, -1,1,-1, -1,-1,1, -1,-1,-1/

c
      trace = a(1,1) + a(2,2) + a(3,3)
c	write(6,*),' trace ', trace
      if (trace .gt.3.) then
        if ( (trace-3.) .gt. 1e-4) write(6,*) 'The trace of the '//
     1 'matrix is',trace,' more than 3. I set to 3 here.'
         trace = 3.
      else if (trace .lt.-1.) then
         if (( trace+1 ) .lt. -1e-4) write(6,*) 'The trace of the '//
     1 'matrix is',trace,' less than -1. I set to -1 here.'
         trace = -1.
      end if
      
      coskapa = (trace - 1.) / 2.
      vkapa = acosd(coskapa)
c	write(6,*),' coskapa, kapa ', coskapa, kapa
c
      if (abs(vkapa).lt.0.03) then
       vec(1) = 0.
       vec(2) = 0.
       vec(3) = 0.
       return
      else
       do i = 1, 3
        temp = (a(i,i)-coskapa)/(1.-coskapa)
        if (temp.lt.0.) then
         if (temp.lt.-1e-4.and.(coskapa.ne.-1.)) then
          write(6,*) 'Warning: your matrix might be not a orthogonal'//
     1 ' one'
          write(6,'(a,i1,a,i1,a)') 'vec(',i,')**2 ='
          write(6,*) temp,' coskapa,kapa,trace',coskapa,vkapa,trace
         end if
         temp = 0.
        end if
        vec1(i) = sqrt(temp)
       end do
       vv = vem(3,vec1)
cc	 if (abs(vv-1.).gt.1e-4) then
c	  write(6,*) 'Warning: modulus of vector is not 0, but',vv
c	 end if
       call arrmc(3,1,vec1,1./vv,vec1)
      end if
c
c	if (coskapa.eq.-1.) return
      errmin = 1.0
      ireal = 0
      do is = 1, 8
        do i =1, 3
         vec(i) = vec1(i) * sgn(i,is)
        end do
c
        call rotvsvec(vec,vkapa,a1)
        call arrps(3,3,a1,a,a2)
c
        err = 0.
        do ip = 1,9
         err = abs(a2(ip)) + err
        end do
      if (err.lt.errmin) ireal = is
      errmin = min(errmin,err)
c
10	end do
c
      if (ireal. eq. 0) then
c
       write(6,*) 'something wrong in this program a,a1'
       do is = 1, 8
        do i =1, 3
         vec(i) = vec1(i) * sgn(i,is)
        end do
c
        call rotvsvec(vec,vkapa,a1)
        call arrps(3,3,a1,a,a2)
        write(6,*) 'a',a
        write(6,*) 'a1',a1
        write(6,*) 'a2',a2
        err = 0.
        do ip = 1,9
         err = abs(a2(ip)) + err
        end do
        write(6,*) vec,vkapa,'vk'
        write(6,*) 'Error',err
c
       end do
       stop ' '
      else 
        do i =1, 3
         vec(i) = vec1(i) * sgn(i,ireal)
        end do
        if (errmin.gt.0.01) then
         write(6,*) 'ireal =', ireal
         call rotvsvec(vec,vkapa,a1)
         call arrps(3,3,a1,a,a2)
         write(6,*) a
         write(6,*) a1
         write(6,*) a2
         write(6,*) vec
        end if
      end if
c check
c	v2 = asind(sinkapa)
c	if (coskapa.lt.0) v2 = -180. - v2
c
c	WRITE(6,*)  'kapa1,kapa2,error',v2,vkapa,abs(v2-vkapa)
c
      end
c
c
      subroutine nodir(name1,name2,len2)
c a subroutine which extracts the file name from a complete filename.
c e.g. input name1 = /nfs/pdb/full/1cnd.pdb
c      output name2 = 1cnd.pdb
c      len2 is the length of the second word
      character*(*) name1,name2
c	i1 = index(name1,'.pdb')
      i1 = lenstr(name1)
      do i = i1, 1, -1
       if (name1(i:i).eq.'/') goto 20
      end do
      i = 0
20	continue
      len2 = i1 - i
      name2(1:len2) = name1(i+1:i1)
      end
c
c
      SUBROUTINE ORIEN ( NATM, XYZ1, XYZ2, A )
c The subroutine selects three atoms to calculate the initial superposing
c matrix A from two molecules in 3*4 array XYZ1 and XYZ2 individually.
c Parameter description:
c NATM is the number of the atoms
c	DIMENSION XYZ1(3,NATM), XYZ2(3,NATM), A(3,3)
c XYZ1(3,NATM) represent the XYZ of molecule 1
c XYZ2(3,NATM) represent the XYZ of molecule 2
c A is the output superposing matrix.
c	COMMON/IAT/ IAT1,IAT2,IAT3,IAT
c	DATA IAT/0/
c Atom IAT1, IAT2 and IAT3 is used to calculate the matrix. It could be
c decided outside the routine if IAT not eq 0.
c
c 
      COMMON/IAT/ IAT1,IAT2,IAT3,IAT
c	DATA IAT/0/,IAT1/1/,IAT2/2/,IAT3/3/
      DIMENSION XYZ1(3,NATM), XYZ2(3,NATM), A(3,3)
      DIMENSION V11(3),V12(3)
      DIMENSION V21(3),V22(3)
      DIMENSION C1(3,3),C2(3,3),C1P(3,3)
      DIMENSION B1(3),B2(3)
c local elements
      REAL T(3)
c
c If IAT = 0, define the three atoms used to calculate the 
c  initial orientation matrix
c
      IF (IAT.EQ.0) THEN
       IAT1 = 1
       IF (NATM.GE.6) THEN
        IAT2 = INT(NATM/3) + 1
        IAT3 = 2*INT(NATM/3) + 1
       ELSE
        IAT2 = 2
c	  IAT3 = 3
        IAT3 = NATM
       END IF
      END IF
c
15	if (dist(xyz1(1,iat1),xyz1(1,iat2)).lt.1.e-3.or.
     +      dist(xyz2(1,iat1),xyz2(1,iat2)).lt.1.e-3) then
       iat1=iat1+1
       if (iat1.lt.iat2) goto 15
       call elize(3,a)
       print*, 'Warning: not give any matrix: -orien'
       return
      end if
16	if (dist(xyz1(1,iat2),xyz1(1,iat3)).lt.1.e-3.or.
     +      dist(xyz2(1,iat2),xyz2(1,iat3)).lt.1.e-3) then
       iat2=iat2+1
       if (iat2.lt.iat3) goto 15
       call elize(3,a)
       print*, 'Warning: not give any matrix: -orien'
      end if
C
C Calculate the initial matrix and the Eulerian angle.
c 
C
C  calculate the matrix 1
20	continue
      CALL ARRPS(3,1,XYZ1(1,IAT2),XYZ1(1,IAT1),V11)  
      CALL ARRPS(3,1,XYZ1(1,IAT2),XYZ1(1,IAT3),V12)
      IF (VEM(3,V12).LT.1E-2) GOTO 50
      CALL ARRMC (3,1,V12,1./VEM(3,V12),C1(1,1))
      CALL VECCRSMLT (V12,V11,B1)
      IF (VEM(3,B1).LT.1E-2) GOTO 50
      CALL ARRMC(3,1,B1,1./VEM(3,B1),C1(1,3))
      CALL VECCRSMLT (C1(1,3),C1(1,1),C1(1,2))

C  calculate the matrix 2
      CALL ARRPS(3,1,XYZ2(1,IAT2),XYZ2(1,IAT1),V21)
      CALL ARRPS(3,1,XYZ2(1,IAT2),XYZ2(1,IAT3),V22)
      IF (VEM(3,V12).LT.1E-2) GOTO 50
      CALL ARRMC (3,1,V22,1./VEM(3,V22),C2(1,1))
      CALL VECCRSMLT (V22,V21,B1)
      IF (VEM(3,B1).LT.1E-2) GOTO 50
      CALL ARRMC(3,1,B1,1./VEM(3,B1),C2(1,3))
      CALL VECCRSMLT (C2(1,3),C2(1,1),C2(1,2))

c calculate the matrix A = c2 * c1^-1  (x2= A * x1)
      call ANTIARR(3,3,C1,C1P)
c	print*, 'iat1,iat2,iat3',iat1,iat2,iat3
      CALL MATMULT(3,3,3,3,C2,C1P,A)
C
      return
50	iat2 = iat2 + 1
      if (iat2.ne.iat3.and.iat2.le.natm) goto 20
      WRITE(6,*) 'warning: cannot give the orientation here'
      WRITE(6,*) 'natm,iat,iat1,iat2,iat3',natm,iat,iat1,iat2,iat3
      call elize(3,a)
      call arrvalue(3,t,0.)
      END
c
c
C-----------APPENDIX PROGRAM FOR THE NEWGIN SYSTEM ANGLE
C              WHEN INOPT1=8
      SUBROUTINE OXFORD(ROT,E)
      DIMENSION E(3,3),ROT(3)
      external cosd, sind
2000	A1=ROT(1)
        A2=ROT(2)
        A3=ROT(3)
      E(1,1)=COSD(A1)*COSD(A2)*COSD(A3)-SIND(A1)
     1*SIND(A3)
      E(2,1)=SIND(A1)*COSD(A2)*COSD(A3)+COSD(A1)
     1*SIND(A3)
      E(3,1)=-(SIND(A2)*COSD(A3))
      E(1,2)=-COSD(A1)*COSD(A2)*SIND(A3)-SIND(A1)
     1*COSD(A3)
      E(2,2)=-SIND(A1)*COSD(A2)*SIND(A3)+COSD(A1)
     1*COSD(A3)
      E(3,2)=SIND(A2)*SIND(A3)
      E(1,3)=COSD(A1)*SIND(A2)
      E(2,3)=SIND(A1)*SIND(A2)
      E(3,3)=COSD(A2)
      END
c
c
      subroutine packexpnd(cell,natom,xyz,nsym,sym)
      parameter(maxatm=25000)
      real xyz(3,maxatm),sym(3,4,160)
      real cell(6),deor(3,3),orth(3,3)
      real cells(6),deors(3,3),orths(3,3)
      real sym1(3,4,160)
      real tmove0(3),ave(3)
      real b1(3),b2(3),b3(3),b4(3),b5(3)
      real tmove(3,27)
c...  data statements.  Separate declaration and init required for f2c
      data tmove
     1 /-1.,-1.,-1.,   -1.,-1.,0.,   -1.,-1.,1.,
     2  -1., 0.,-1.,   -1., 0.,0.,   -1., 0.,1.,
     3  -1., 1.,-1.,   -1., 1.,0.,   -1., 1.,1.,
     4   0.,-1.,-1.,   -0.,-1.,0.,   -0.,-1.,1.,
     5   0., 0.,-1.,    0. ,0.,0.,    0., 0.,1.,
     6   0., 1.,-1.,    0., 1.,0.,    0., 1.,1.,
     7   1.,-1.,-1.,    1.,-1.,0.,    1.,-1.,1.,
     8   1., 0.,-1.,    1., 0.,0.,    1., 0.,1.,
     9   1., 1.,-1.,    1., 1.,0.,    1., 1.,1./
c
      nmove = 1
      call lgg_crystal(CELL,CELLS,DEOR,ORTH,deors,orths)
      call averg(3,natom,xyz,ave)
c	WRITE(6,*)  'Average xyz',ave
      do isym = 1, nsym
       call matmult(3,3,3,1,deor,ave,b1)
       call matmult(3,3,3,1,sym(1,1,isym),b1,b2)
       call arrad(3,1,b2,sym(1,4,isym),b3)
       call arrgive(3,b3,b4)
c	 WRITE(6,*) 'b3',b4
       do i = 1, 3
        call frcinside(b4(i))
        b5(i) = b4(i) - b3(i)
       end do
c	 WRITE(6,*) 'b4-b5',b4,b5
       call arrad(3,1,b5,sym(1,4,isym),tmove0)
c	 WRITE(6,*)  'tmove0',tmove0
       call arrgive(9,sym(1,1,isym),sym1(1,1,nmove))
       call arrgive(3,tmove0,sym1(1,4,nmove))
c	 write(6,*) 'Operation:',nmove
c	 write(6,11) ((sym1(i,j,nmove),j=1,4),i=1,3)
       nmove = nmove + 1
       call arrgive(12,sym(1,1,isym),sym1(1,1,nmove))
       do jmove = 1,27
        if (jmove.eq.14) goto 40
        call arrad(3,1,tmove(1,jmove),tmove0,sym1(1,4,nmove))
        do iatom = 1, natom
         call matmult(3,3,3,1,deor,xyz(1,iatom),b1)
         call matmult(3,3,3,1,sym1(1,1,nmove),b1,b2)
         call arrad(3,1,b2,sym1(1,4,nmove),b3)
         if (b3(1).ge.0. .and. b3(1).lt.1. .and.
     1      b3(2).ge.0. .and. b3(2).lt.1. .and.
     2      b3(3).ge.0. .and. b3(3).lt.1. ) then
c	    write(6,*) 'move-iele-iatom',jmove,iele,iatom
c	    WRITE(6,*) 'xyz-iatom',xyz(1,iatom),xyz(3,iatom),xyz(3,iatom)
c	    WRITE(6,*) 'b1-deor',b1
c	    WRITE(6,*) 'b2-symrot',b2
c	    WRITE(6,*) 'b3-jmove',b3
c	    write(6,*) 'Operation:',nmove
c	    write(6,11) ((sym1(i,j,nmove),j=1,4),i=1,3)
          nmove = nmove + 1
          call arrgive(12,sym(1,1,isym),sym1(1,1,nmove))
          goto 40
         end if
30	  end do
40	 end do
      end do
11	format(4f10.4)
c
      nsym = nmove - 1
c	write(6,*) 'Total',nsym,' operations can put the molecule into
c	1 the unit cell'
      call arrgive(nsym*12,sym1,sym)
c
      end
c
c
      CHARACTER*16 FUNCTION PLTNAM(IDUM3)
c	CALL SYSTEM('printenv USER > PLTNAM')
c       call ccpdpn(55,'PLTNAM','UNKNOWN','F',0,0)
c	read(55,'(a)') pltnam
c	do i = 1, 16
c	if (pltnam(i:i).eq.' ') goto 10
c	end do
cc10	continue
c	do j = i,16
c	 pltnam(j:j) = ' '
c	end do
      call getlog(pltnam)
      call up(pltnam,16)
      END
c
c
C------------STANDARD FUNCTION SUBPROGRAM FOR
C----------------DOT PRODUCT OF TWO VECTORS
      FUNCTION POIMULT(N1,N2,V1,V2)
      DIMENSION V1(N1),V2(N2)
      POIMULT=0.
      IF (N1.NE.N2) GOTO 90
      DO 10 I=1,N1
10	POIMULT=V1(I)*V2(I)+POIMULT
      RETURN
90	PRINT*,'POIMULT> Two vectors do not have same dimension'
      END
c
c
      SUBROUTINE POLOR(POL,E)
      DIMENSION E(3,3),POL(3)
      EQUIVALENCE  (S1,SIPS),  (S2,SITH),  (S3,SIPH)
     1  ,(C1,COPS),  (C2,COTH),  (C3,COPH)

C
C * POLAR ANGLE MATRIX  (TRANSPOSED ACCORDING TO TOLLIN, ERRATA OF R&B)
C
      ARC=3.14159265359/180.0
      PS=POL(1)
      TH=POL(2)
      PH=POL(3)

      SIPS = SIN(PS*ARC)
      COPS = COS(PS*ARC)
      SITH = SIN(TH*ARC)
      COTH = COS(TH*ARC)
      SIPH = SIN(PH*ARC)
      COPH = COS(PH*ARC)

130	S1SQ = S1*S1
      C1SQ = C1*C1
      CPC  = 1.0 - C3
      E(1,1) = C3             + S1SQ*C2*C2*CPC
      E(2,1) = S1*C1*C2*CPC   - S1*S2*S3
      E(3,1) =-S1SQ*C2*S2*CPC - C1*S3
      E(1,2) = S1*C1*C2*CPC   + S1*S2*S3
      E(2,2) = C3             + C1SQ*CPC
      E(3,2) =-S1*C1*S2*CPC   + S1*C2*S3
      E(1,3) =-S1SQ*S2*C2*CPC + C1*S3
      E(2,3) =-S1*C1*S2*CPC   - S1*C2*S3
      E(3,3) = C3             + S1SQ*S2*S2*CPC
      END
c
c
      subroutine polors(pol,e)
c my program to calculate matrix from polar angle
c pol is psi phi kappa
c where psi is angle between axis and Y
c 	phi is angle between image of axis in XZ and X
c	kappa is rotating angle

      dimension e(3,3),pol(3)
      dimension b1(3)
      external cosd, sind
      b1(2) = cosd(pol(1))	
      b1(1) = sind(pol(1))*cosd(pol(2))	
      b1(3) = - sind(pol(1))*sind(pol(2))
      call rotvsvec(b1,pol(3),e)
      end
c
c
      subroutine polorz(pol,e)
c my program to calculate matrix from polar angle (z system)
c pol is psi phi kappa
c where psi is angle between axis and Z
c 	phi is angle between image of axis in XY and X
c	kappa is rotating angle

      real e(3,3),pol(3)
      real b1(3)
      external cosd, sind
c	b1(2) = cosd(pol(1))	
c	b1(1) = sind(pol(1))*cosd(pol(2))	
c	b1(3) = - sind(pol(1))*sind(pol(2))
      b1(3) = cosd(pol(1))	
      b1(1) = sind(pol(1))*cosd(pol(2))	
      b1(2) = sind(pol(1))*sind(pol(2))
      call rotvsvec(b1,pol(3),e)
      end
c
c
      SUBROUTINE POS2VEC(NATM,POS,VEC)
C This subroutine is to calculate the center-atom vector from the 
c coordinate. NATM vector will be produced from NATM atoms when NATM >= 3. 
c The NATM vectors is CEN->atom2, CEN->atom3,....CEN->atomN, where
C CEN is the centre of the NATM atoms It is used to superimpose the 
c orientation of two molecules.
c                                  written by Guoguang Lu
c 
C	DIMENSION POS(3,NATM), VEC(3,NATM)
C Parameter description:
C       NATM is the number of atoms
c       POS(1->3,i) is the coordinate of the i'th atom
c       VEC  is the output vectors
c
      DIMENSION POS(3,NATM), VEC(3,NATM),C(3)
      IF (NATM.LT.3) STOP 'Number of the atoms is less than 3'
      CALL AVERG(3,NATM,POS,C)
C  C is the center coordinate of the atoms
      DO I=1,NATM
      CALL ARRPS(3,1,POS(1,I),C,VEC(1,I))
C	CALL ARRMC (3,1,VEC(1,I),1./VEM(3,VEC(1,I)),VEC(1,I))
      END DO
      END
c
C
C *  R&B - MATRIX
C
      SUBROUTINE R_B(ROT,E,IOTP)
      DIMENSION E(3,3),ROT(3)
      EQUIVALENCE  (PS,AL,PA), (TH,BE,PB), (PH,GA,PC)
      EQUIVALENCE  (S1,SIPS),  (S2,SITH),  (S3,SIPH)
     1   ,(C1,COPS),  (C2,COTH),  (C3,COPH)

      PS=ROT(1)
      TH=ROT(2)
      PH=ROT(3)

      IF (IOTP.EQ.4) THEN
      PS=PS+90.0
      PH=PH-90.0
      END IF

      ARC=3.14159265359/180.0
      SIPS = SIN(PS*ARC)
      COPS = COS(PS*ARC)
      SITH = SIN(TH*ARC)
      COTH = COS(TH*ARC)
      SIPH = SIN(PH*ARC)
      COPH = COS(PH*ARC)

120	E(1,1) =-S1*C2*S3 + C1*C3
      E(2,1) = C1*C2*S3 + S1*C3
      E(3,1) =            S2*S3
      E(1,2) =-S1*C2*C3 - C1*S3
      E(2,2) = C1*C2*C3 - S1*S3
      E(3,2) =            S2*C3
      E(1,3) =            S1*S2
      E(2,3) =          - C1*S2
      E(3,3) = C2
      END 
c
c
      real function lgg_radii(natm,xyz)
c pjx: renamed from radii to avoid clash with common block
c of same name in refmac4
c Explicitly typed to real
c
c--calculate average radii of a group of atoms.
c first calculate the center of gravity
c then the mean distance from each atom to this center.
c
      real xyz(3,natm)
      real b1(3),b2(3)
c
      lgg_radii = 0.
      call averg(3,natm,xyz,b1)
      do i = 1, natm
       call arrps(3,1,xyz(1,i),b1,b2)
       lgg_radii = lgg_radii + vem(3,b2)
      end do
      lgg_radii = lgg_radii/float(natm)
      end
c
c
c this is a subroutine to change orthogonalization matrix according to
c cell dimension raxis spindle and xray axis 
c input 
c character cell(6)   ------ cell dimensions
c character aspin*3   ------ raxis spindle axis
c character axray*2   ------ raxis xray axis
c output
c CELLS*6 --- cell dimensions in reciprocal space
c DEOR3*3 --- Deorthogonalization matrix in real space
c ORTH3*3 --- Orthogonalization matrix in receprocal space.
c DEORS3*3 --- Deorthogonalization matrix in reciprocal space
c ORTHS3*3 --- Orthogonalization matrix in reciprocal space.
c AND COMMON/DET/ DET 3*3 MATRIX transfer from statnd X-a orthogonal matrix
c
c Guoguang 930723
c
c PJX 010125 renamed common block det to lggdet to avoid
c clash with subroutine det in rsps.
c This common block doesn't appear to be used anywhere else?
c
      subroutine raxcrystl(cell,aspin,axray,cells,deor,orth,deors,orths)
      common/lggdet/ det
      real det(3,3),buf(3,3),buf2(3,3)
      real cell(6),cells(6)
      real deor(3,3),orth(3,3)
      real deors(3,3),orths(3,3)
      real pdir(3,6)
      INTEGER IC(6)
      real cof(6)
      character*3 aspin,pspin(6)
      character*2 axray,pxray(6)
      DIMENSION B1(3),B2(3),IUF(3)
c...  data statements.  Separate declaration and init required for f2c
      data pdir /1.,0.,0.,0.,1.,0.,0.,0.,1.,
     1	-1.,0.,0.,0.,-1.,0.,0.,0.,-1./
      data ic /1,2,3,1,2,3/
      data cof /1.,1.,1.,-1.,-1.,-1./
      data pspin /'+a*','+b*','+c*','-a*','-b*','-c*'/
      data pxray /'+a','+b','+c','-a','-b','-c'/
c
c standard orthogonalization matrix
c
      call lgg_crystal(cell,cells,deor,orth,deors,orths)
c
c transferring matrix from new system to old det
c
c xray axis
      do i = 1, 6
       if (axray.eq.pxray(i)) then
        call arrmc(3,1,orth(1,ic(i)),
     1  cof(i)/vem(3,orth(1,ic(i))),det(1,1))
        goto 30
       end if
      end do
        WRITE(6,*) 'X-ray axis',(pxray(i),' ',i=1,6)
      WRITE(6,*) 'but you have ',axray
        stop 'wrong xray  axis'
30	continue
c spindle axis
        do i = 1, 6
         if (aspin.eq.pspin(i)) then
        call arrmc(3,1,orths(1,ic(i)),
     1   cof(i)/vem(3,orths(1,ic(i))),det(1,3))
        goto 40
       end if
      end do
        WRITE(6,*) 'spindle axis',(pspin(i),' ',i=1,6)
      WRITE(6,*) 'but you have ',aspin
        stop 'wrong xray  axis'
c
40	continue
      call veccrsmlt(det(1,3),det(1,1),det(1,2))
c
      call arrgive(9,orth,buf)
      call antiarr(3,3,det,buf2)
      call matmult(3,3,3,3,buf2,buf,orth)
c
      CALL ARRGIVE(9,ORTH,DEOR)
      CALL IVSN(3,DEOR,B1,B2,IUF,DE,1.0E-6)
      CALL ANTIARR(3,3,ORTH,DEORS)
      CALL ANTIARR(3,3,DEOR,ORTHS)
C
c	DO I = 1, 3
c	CELLS(I) = VEM(3,orths(1,I))
c	END DO
C
c	CELLS(4) = ANGLE(ORTHS(1,2),ORTHS(1,3))	
c	CELLS(5) = ANGLE(ORTHS(1,3),ORTHS(1,1))	
c	CELLS(6) = ANGLE(ORTHS(1,1),ORTHS(1,2))	
c
      end
c
c
c this is a subroutine to change raxis spindle and xray axis into U matrix
c input 
c character aspin*3   ------ raxis spindle axis
c character axray*2   ------ raxis xray axis
c real*4 umat(3,3)     ------ raxis U matrix
c Guoguang 930723
c
      subroutine raxumat(aspin,axray,umat)
      real umat(3,3)
      real rmat(3,3)
      real pdir(3,6)
      character*3 aspin,pspin(6)
      character*2 axray,pxray(6)
c...  data statements.  Separate declaration and init required for f2c
      data pdir /1.,0.,0.,0.,1.,0.,0.,0.,1.,
     1	-1.,0.,0.,0.,-1.,0.,0.,0.,-1./
      data pspin /'+a*','+b*','+c*','-a*','-b*','-c*'/
      data pxray /'+a','+b','+c','-a','-b','-c'/
c
      do i = 1, 6
       if (aspin.eq.pspin(i)) then
        call arrgive(3,pdir(1,i),umat(1,3))
        goto 20
       end if 
      end do
      WRITE(6,*) ' allowed Spindle axis',(pspin(i),' ',i=1,6)
      stop 'wrong spindle axis'
20	continue
c
      do i = 1, 6
       if (axray.eq.pxray(i)) then
        call arrgive(3,pdir(1,i),umat(1,1))
        goto 30
       end if 
      end do
      WRITE(6,*) 'X-ray axis',(pxray(i),' ',i=1,6)
      stop 'wrong xray  axis'
30	continue
      call veccrsmlt(umat(1,3),umat(1,1),umat(1,2))
      call antiarr(3,3,umat,rmat)
      call arrgive(9,rmat,umat)
c
      end
c
c
      SUBROUTINE RECEPICAL(CELL,DEOR,ORTH)
C This is a routine which inputs cell parameters, CELL and converts 
c them into reciprocal cell parameters. It then gives a deorthogonalization
c and orthogonalization matrix in reciprocal space.
c
      DIMENSION DEOR(3,3),ORTH(3,3),CELL(6)
      DIMENSION BUFF(3,3)
      DIMENSION B1(3),B2(3),IUF(3)
      external acosd, cosd, sind
      AP=CELL(4)
      BA=CELL(5)
      GA=CELL(6)
      CALL ARRMC(3,3,BUFF,0.,BUFF)
      COASTAR=(COSD(BA)*COSD(GA)-COSD(AP))/(SIND(BA)
     1*SIND(GA))
      SIASTAR=SQRT(1.-COASTAR*COASTAR)
      BUFF(1,1)=1.	
      BUFF(1,2)=COSD(GA)
      BUFF(2,2)=SIND(GA)
      BUFF(1,3)=COSD(BA)
      BUFF(2,3)=-SIND(BA)*COASTAR
      BUFF(3,3)=SIND(BA)*SIASTAR
C
      DO I =1, 3
      CALL ARRMC(3,1,BUFF(1,I),CELL(I),BUFF(1,I))
      END DO
c
c	WRITE(6,*)  buff
C
      VOLUM = VLDIM3(BUFF)
C
      CALL VECCRSMLT(BUFF(1,2),BUFF(1,3),DEOR(1,1))
      CALL VECCRSMLT(BUFF(1,3),BUFF(1,1),DEOR(1,2))
      CALL VECCRSMLT(BUFF(1,1),BUFF(1,2),DEOR(1,3))
      CALL ARRMC(3,3,DEOR,1./VOLUM,DEOR)
C
c	WRITE(6,*) volum
c	WRITE(6,*) i,cell(i),vem(3,deor(1,i))
      DO I = 1, 3
      CELL(I) = VEM(3,DEOR(1,I))
c	WRITE(6,*) CELL(I)
c	WRITE(6,*) DEOR
      CALL ARRMC(3,1,DEOR(1,I),1./CELL(I),DEOR(1,i))
      END DO
C
      CELL(4) = ACOSD(POIMULT(3,3,DEOR(1,2),DEOR(1,3)))
      CELL(5) = ACOSD(POIMULT(3,3,DEOR(1,3),DEOR(1,1)))
      CELL(6) = ACOSD(POIMULT(3,3,DEOR(1,1),DEOR(1,2)))
      
C
      CALL ARRMC(3,3,DEOR,1.,ORTH)
      CALL IVSN(3,ORTH,B1,B2,IUF,VAL,1E-6)
      RETURN
      END
c
c
      SUBROUTINE RECEPICAL0(CELL,DEOR,ORTH)
C This is a routine which inputs cell parameters, CELL and converts 
c them into reciprocal cell parameters. It then gives a deorthogonalization
c and orthogonalization matrix in reciprocal space.
c
      DIMENSION DEOR(3,3),ORTH(3,3),CELL(6)
      DIMENSION BUFF(3,3)
      DIMENSION B1(3),B2(3),IUF(3)
      external acosd, cosd, sind
      AP=CELL(4)
      BA=CELL(5)
      GA=CELL(6)
      CALL ARRVALUE(3*3,BUFF,0.)
      COASTAR=(COSD(BA)*COSD(GA)-COSD(AP))/(SIND(BA)
     1*SIND(GA))
      SIASTAR=SQRT(1.-COASTAR*COASTAR)
      BUFF(1,1)=1.	
      BUFF(1,2)=COSD(GA)
      BUFF(2,2)=SIND(GA)
      BUFF(1,3)=COSD(BA)
      BUFF(2,3)=-SIND(BA)*COASTAR
      BUFF(3,3)=SIND(BA)*SIASTAR
C
      DO I =1, 3
      CALL ARRMC(3,1,BUFF(1,I),CELL(I),BUFF(1,I))
      END DO
c
c	WRITE(6,*)  buff
C
      VOLUM = VLDIM3(BUFF)
C
      CALL VECCRSMLT(BUFF(1,2),BUFF(1,3),ORTH(1,1))
      CALL VECCRSMLT(BUFF(1,3),BUFF(1,1),ORTH(1,2))
      CALL VECCRSMLT(BUFF(1,1),BUFF(1,2),ORTH(1,3))
      CALL ARRMC(3,3,ORTH,1./VOLUM,ORTH)
C
c	WRITE(6,*) volum
c	WRITE(6,*) i,cell(i),vem(3,ORTH(1,i))
      DO I = 1, 3
      CELL(I) = VEM(3,ORTH(1,I))
c	WRITE(6,*) CELL(I)
c	WRITE(6,*) ORTH
      CALL ARRMC(3,1,ORTH(1,I),1./CELL(I),DEOR(1,i))
      END DO
C
      CELL(4) = ACOSD(POIMULT(3,3,DEOR(1,2),DEOR(1,3)))
      CELL(5) = ACOSD(POIMULT(3,3,DEOR(1,3),DEOR(1,1)))
      CELL(6) = ACOSD(POIMULT(3,3,DEOR(1,1),DEOR(1,2)))
      
C
      CALL ARRGIVE(3*3,ORTH,DEOR)
      CALL IVSN(3,DEOR,B1,B2,IUF,VAL,1E-6)
      RETURN
      END
c
c
      SUBROUTINE REDSTRIN(NUNIT,NCHA,TXT,NPAR)
C A subroutine to write a string of characters in unit NUNIT according 
c to the spaces in the string. It is used to change a character into 
c parameters under help of a file NUNIT.
c NUNIT is the unit number.
c NCHA is length of character string.
c TXT is the NCHA bytes character.
C NPAR is how many parameters in this string.
      CHARACTER*120 TXT
      NPAR = 0
      REWIND (NUNIT)
      I = 1
10	IF (I.LE.NCHA) THEN
         IF (TXT(I:I).NE.' ') THEN
            J = I
20	      IF (J.GE.NCHA) THEN
               WRITE(NUNIT,'(A)') TXT(I:J)
               NPAR = NPAR + 1
            ELSE IF (TXT(J+1:J+1).EQ.' ') THEN
               WRITE(NUNIT,'(A)') TXT(I:J)
               NPAR = NPAR + 1         
            ELSE	
               J = J + 1
               GOTO 20
            END IF
            I = J + 1
            GOTO 10
         ELSE IF (I.LT.NCHA) THEN
            I = I + 1
            GOTO 10
         END IF
      END IF
C
      REWIND(NUNIT)
      END
c
c
      subroutine refmtoroh(amat,roh,roh2,rms)
c This is a subroutine to calculate from a matrix to Eulerian angles in
c Rober O Huber system and refine ROH angles to minimize 
c Sigma((Aroh(i,j)-Amat(i,j)**2). The output will include two
c symmetric ROH angles a rms value of input matrix and the matrix
c from ROH angles.
c Input:	amat(3,3) --- Input matrix
c Output:	roh(3,3) --- output ROH angle
c		roh2(3,3) --- output ROH angle
c		rms       --- rms value between amat and mat from ROH angle.
c Guoguang Lu, 950831. Purdue Univ
c
      real amat(3,3),roh(3),roh2(3)
      real bmat(3,3),b1(3),b2(3),cmat(9)
      real rohnew(3)
      real vt(9,3),vvv(3,3),delta(3),det(3)
c	real droh(3,3,3)
c
      ncyc = 0
      deg = 180./3.141592654
      call mtohuber(amat,roh,roh2)
      call huber(roh,bmat)
      call arrps(9,1,bmat,amat,cmat)
c	WRITE(6,*) 'cmat',cmat
      rms = vem(9,cmat)
c11	format(3f10.7)
c	write(6,*) 'bmat'
c	write(6,11) bmat
c
10	continue
c	WRITE(6,*) 'cycle,rms',rms,ncyc
c	WRITE(6,*) 'roh',roh
      ncyc = ncyc + 1
c
      call drvrohthd(1,roh,vt(1,1))
      call drvrohthd(2,roh,vt(1,2))
      call drvrohthd(3,roh,vt(1,3))
c
c	write(6,*) 'vt'
c	write(6,'(3f12.7,2x,f12.7)') ((vt(i,j),j=1,3),cmat(i),i=1,9)
      call lsqeq(9,3,vt,cmat,det,vvv,b1)
c	WRITE(6,*) 'det',det
      call arrps(3,1,roh,det,rohnew)
c	WRITE(6,*) 'rohnew',rohnew
      call huber(rohnew,bmat)
c	write(6,*) 'bmatnew'
c	write(6,11) bmat
      call arrps(9,1,bmat,amat,cmat)
c	WRITE(6,*) 'cmatnew',cmat
      rms1 = vem(9,cmat)
      if (rms1.lt.rms) then
       call arrgive(3,rohnew,roh)
       rms = rms1
       goto 10
      else
c	 WRITE(6,*) 'Refinement finished with rms,rms1:',rms,rms1
       roh2(1) = 180. + roh(1)
       roh2(2) = 180. - roh(2)
       roh2(3) = 180. + roh(3)
      end if
c
      end
c
c
      SUBROUTINE REFORN ( NATM, XYZ1, XYZ2, A, VT, DIS )
C	DIMENSION XYZ1(3,NATM), XYZ2(3,NATM), A(3,3)
C	DIMENSION VT(3*NATM,3), DIS(NATM,3)
C A subroutine to give the superimpose matrix and vector from MOL1 to 
c MOL2 i.e. XYZ2(3*NATM) = A(3*3) * XYZ1(3*NATM) + T(3)
C Where    NATM represents the number of the atoms in each molecule
c                which will be superimposed.
c          XYZ1(3,NATM) represents the coordinates of MOL1
c          XYZ2(3,NATM) represents the coordinates of MOL2
C          A(3*3)       represents the input initial matrix and 
c                       output matrix.
C          VT(3*6*NATM)   is a 6*NATM-dimensional buffer array
C          DIS(3*NATM)  is a 3*NATM-dimensional buffer array
c
C	COMMON/SHIFS/ SCALR,SCALT
C Note: There should be a scale of the shifts. i.e. Shift= scale * shifts
c Here SCALR is the rotation shiftscale
c Here SCALT is the translation shiftscale
c The initial and suggested SCALR is 1.
c The initial and suggested SCALt is 1.
C
C	COMMON/RMS/ RMS,SMEAN,NREF,NREF1
C RMS is the final R.M.S between two molecules.
C SMEAN is the mean distance.
c NREF is the number of cycles of rotation refinement
c NREF1 is the number of cycles of translation refinement
c They could be used to judge the SCALR and SCALS
c
c The routine uses the atom-atom vector to perform the Eulerian angle least 
c square refinement. Testing shows that that this method might be more
c accurate than others at least in orientation.
c
C SUPOS1 is almost same with SUPOS, but used as a subroutine as a first
c step of the refinement by subroutine SUPRIMP.
c
c                                           by Guoguang Lu
cv
c     added an upper limit for number of refinement cycles to do.
c     C. Vonrhein
cv
C                                               30/01/1989
C
C
      COMMON/RMS/ RMS,SMEAN,NREF,NREF1
c	COMMON/SHIFS/ SCALR,SCALS
      PARAMETER (NATM0 = 50000)
c   NATM0 is the largest number of atoms.
      DIMENSION XYZ1(3,NATM), XYZ2(3,NATM), A(3,3)
      DIMENSION VT(3*NATM,3), DIS(3,NATM)
      DIMENSION VEC1(3,NATM0),VEC2(3,NATM0),VEC1P(3,NATM0)
C
      DIMENSION DA(3,3,3),VVV(3,3)
      DIMENSION B1(3),B2(3),ROHOLD(3),ROH(3),DETAROH(3),TOLD(3)
C MAXREF = maximum number of refinements to do
C
      INTEGER MAXREF
      PARAMETER (MAXREF=100)
C
c...  data statements.  Separate declaration and init required for f2c
      DATA SCALR/1./,SCALT/1./
C
      DEG=180./3.14159265359
      NREF=0
c	NREF1=0
C
      CALL MTOHUBER(A,ROH,B1)
c
c calculate the atom-atom vector
c
      CALL POS2VEC(NATM,XYZ1,VEC1)
      CALL POS2VEC(NATM,XYZ2,VEC2)
C
      CALL MATMULT(3,3,3,NATM,A,VEC1,VEC1P)
      CALL ARRPS(3,NATM,VEC1P,VEC2,DIS)
C                                  A*V1 - V2
      RMS1= DOSQ( 3*NATM, DIS )
      RMS = sqrt(RMS1)/float(NATM)
C
50	CONTINUE
C 
C  Refine the superimposing Eulerian Angle'
c
C	WRITE(6,'(A)') '0'
C	WRITE(6,*) 'Cycle of Eulerian angle refinement ',NREF
C	WRITE(6,*) 'Superimposing matrix'
C	WRITE(6,10) ((A(I,J),J=1,3),I=1,3)
C	WRITE(6,*) 'Eulerian Angle', ROH
C10	format (3F10.5)
C	WRITE(6,*) 'SIGMA (Deta vector)^2 =', RMS
C
C
c
C COMPUTE THE DETATH1,DETATH2,DETATH3 than the normal matrix
C
      DO I = 1, 3
      CALL DRVROHTH(I,ROH,DA(1,1,I))
      CALL MATMULT(3,3,3,NATM,DA(1,1,I),VEC1,VT(1,I))
C                                (DERIV(A)/DERIV(THETAi))*VEC1
      END DO
c
      CALL LSQEQ(NATM*3,3,VT,DIS,DETAROH,VVV,B1)
C To change it into degree from arc unit
      CALL ARRMC(3,1,DETAROH,DEG,DETAROH)
C
      CALL ARRMC(3,1,DETAROH,-SCALR,DETAROH)
C	WRITE(6,*) 'DetaROH = ',DETAROH	
      CALL ARRMC(3,1,ROH,1.,ROHOLD)
C                            SAVE THE PREVIOUS ROH
      CALL ARRAD(3,1,ROHOLD,DETAROH,ROH)
C                                     ROH=ROHOLD+DETAROH
      CALL HUBER(ROH,A)
C                       NEW MATRIX
      CALL MATMULT(3,3,3,NATM,A,VEC1,VEC1P)
      CALL ARRPS(3,NATM,VEC1P,VEC2,DIS)
C                                  A*V1 - V2
      RMS=dosq(3*NATM,DIS)
C
C     do only a maximum of MAXREF refinement cycles
C
      IF ((NREF.LE.MAXREF).AND.(RMS.LT.RMS1)) THEN
c
        RMS1=RMS
        NREF=NREF+1
        GOTO 50
      END IF

C	WRITE(6,'(A)') '0'
C	WRITE(6,*) 'That is the final rotation result.'
C	WRITE(6,'(A)') '0'
      CALL ARRMC(3,1,ROHOLD,1.,ROH)
      CALL HUBER(ROH,A)

      END 
c
c
      SUBROUTINE REFORNFIN ( NATM, XYZ1, XYZ2, A, VT, DIS )
C	DIMENSION XYZ1(3,NATM), XYZ2(3,NATM), A(3,3)
C	DIMENSION VT(3*NATM,3), DIS(NATM,3)
C A subroutine to give the superimpose matrix and vector from MOL1 to 
c MOL2 i.e. XYZ2(3*NATM) = A(3*3) * XYZ1(3*NATM) + T(3)
C Where    NATM represents the number of the atoms in each molecule
c                which will be superimposed.
c          XYZ1(3,NATM) represents the coordinates of MOL1
c          XYZ2(3,NATM) represents the coordinates of MOL2
C          A(3*3)       represents the input initial matrix and 
c                       output matrix.
C          VT(3*6*NATM)   is a 6*NATM-dimensional buffer array
C          DIS(3*NATM)  is a 3*NATM-dimensional buffer array
c
C	COMMON/SHIFS/ SCALR,SCALT
c 	DATA SCALR/1./,SCALT/1./
C Note: There should be a scale of the shifts. i.e. Shift= scale * shifts
c Here SCALR is the rotation shiftscale
c Here SCALT is the translation shiftscale
c The initial and suggested SCALR is 1.
c The initial and suggested SCALt is 1.
C
C	COMMON/RMS/ RMS,SMEAN,NREF,NREF1
C RMS is the final R.M.S between two molecules.
C SMEAN is the mean distance.
c NREF is the number of cycles of rotation refinement
c NREF1 is the number of cycles of translation refinement
c They could be used to judge the SCALR and SCALS
c
c The program use the atom-atom vector to perform the Eulerian angle least 
c square refinement. Testing shows that this method might be more
c accurate than others at least in orientation.
c
C SUPOS1 is almost the same as SUPOS, but used as a subroutine as a first
c step of the refinement by subroutine SUPRIMP.
c
c                                           by Guoguang Lu
C                                               30/01/1989
C
C
      COMMON/RMS/ RMS,SMEAN,NREF,NREF1
C	COMMON/SHIFS/ SCALR,SCALS
      PARAMETER (NATM0 = 50000)
c   NATM0 is the largest number of atoms.
      DIMENSION XYZ1(3,NATM), XYZ2(3,NATM), A(3,3)
      DIMENSION VT(3*NATM,3), DIS(3,NATM)
      DIMENSION VEC1(3,NATM0),VEC2(3,NATM0),VEC1P(3,NATM0)
C
      DIMENSION DA(3,3,3),VVV(3,3)
      DIMENSION B1(3),B2(3),ROHOLD(3),ROH(3),DETAROH(3),TOLD(3)
C
C
C	DEG=180./3.14159265359
C
      SCALR = 1.
      NREF=0
c	NREF1=0
C
      CALL MTOHUBERARC(A,ROH,B1)
c
c calculate the atom-atom vector
c
      CALL POS2VEC(NATM,XYZ1,VEC1)
      CALL POS2VEC(NATM,XYZ2,VEC2)
C
      CALL MATMULT(3,3,3,NATM,A,VEC1,VEC1P)
      CALL ARRPS(3,NATM,VEC1P,VEC2,DIS)
C                                  A*V1 - V2
      RMS1= DOSQ( 3*NATM, DIS )
      RMS = RMS1
C
50	CONTINUE
C 
C  Refine the superimposing Eulerian Angle'
c
C	WRITE(6,'(A)') '0'
C	WRITE(6,*) 'Cycle of Eulerian angle refinement ',NREF
C	WRITE(6,*) 'Superimposing matrix'
C	WRITE(6,10) ((A(I,J),J=1,3),I=1,3)
C	WRITE(6,*) 'Eulerian Angle', ROH
C10	format (3F10.5)
C	WRITE(6,*) 'SIGMA (Deta vector)^2 =', RMS
C
C
c
C COMPUTE THE DETATH1,DETATH2,DETATH3 than the normal matrix
C
      DO I = 1, 3
      CALL DRVROHTHARC(I,ROH,DA(1,1,I))
      CALL MATMULT(3,3,3,NATM,DA(1,1,I),VEC1,VT(1,I))
C                                (DERIV(A)/DERIV(THETAi))*VEC1
      END DO
c
      CALL LSQEQ(NATM*3,3,VT,DIS,DETAROH,VVV,B1)
C To change it into degree from arc unit
C	CALL ARRMC(3,1,DETAROH,DEG,DETAROH)
C
30	CONTINUE
      CALL ARRMC(3,1,DETAROH,-SCALR,DETAROH)
C	WRITE(6,*) 'SCALR ',SCALR
C	WRITE(6,*) 'DetaROH = ',DETAROH	
      CALL ARRMC(3,1,ROH,1.,ROHOLD)
C                            SAVE THE PREVIOUS ROH
      CALL ARRAD(3,1,ROHOLD,DETAROH,ROH)
C                                     ROH=ROHOLD+DETAROH
      CALL HUBERARC(ROH,A)
C                       NEW MATRIX
      CALL MATMULT(3,3,3,NATM,A,VEC1,VEC1P)
      CALL ARRPS(3,NATM,VEC1P,VEC2,DIS)
C                                  A*V1 - V2
      RMS=dosq(3*NATM,DIS)
C
C	WRITE(6,*)  'RMS NEW AND OLD', RMS,RMS1
      IF (RMS.LT.RMS1) THEN
      RMS1=RMS
      NREF=NREF+1
      GOTO 50
      ELSE IF(SCALR.GT.1E-4) THEN
      CALL ARRMC(3,1,DETAROH,-1./SCALR,DETAROH)
      SCALR = SCALR*0.5
      CALL ARRMC(3,1,ROHOLD,1.,ROH)
      CALL HUBERARC(ROH,A)
      GOTO 30
      END IF

C	WRITE(6,'(A)') '0'
C	WRITE(6,*) 'That is the final rotation result.'
C	WRITE(6,'(A)') '0'
      CALL ARRMC(3,1,ROHOLD,1.,ROH)
      CALL HUBERARC(ROH,A)

      END 
c
c
      SUBROUTINE REFRT ( NATM, X1, X2, A, T, VT, DIS )
C	DIMENSION X1(3,NATM), X2(3,NATM), A(3,3), T(3)
C	DIMENSION VT(3*NATM,6),DIS(3,NATM)
C A subroutine to refine the superimpose matrix and vector from MOL1 to 
c MOL2 i.e. XYZ2(3*NATM) = A(3*3) * XYZ1(3*NATM) + T(3)
C Where    NATM represents the number of the atoms in each molecule
c                which will be superimposed.
c          XYZ1(3,NATM) represents the coordinates of MOL1
c          XYZ2(3,NATM) represents the coordinates of MOL2
C          A(3*3)       represents the input initial matrix and 
c                       output matrix.
c          T(3*3)       represents the input initial translation vector
c                       output vector.
C          VT(3*6*NATM)   is a 6*NATM-dimensional buffer array
C          DIS(3*NATM)  is a 3*NATM-dimensional buffer array
c
C	COMMON/SHIFS/ SCALR,SCALT
c 	DATA SCALR/1./,SCALT/1./
C Note: There should be a scale of the shifts. i.e. Shift= scale * shifts
c Here SCALR is the rotation shiftscale
c Here SCALT is the translation shiftscale
c The initial and suggested SCALR is 60.
c The initial and suggested SCALt is 1.
C
C	COMMON/RMS/ RMS,SMEAN,NREF,NREF1
C RMS is the final R.M.S between two molecules.
C SMEAN is the mean distance.
c NREF is the number of cycles of refinement
c NREF1 is not useful.
c They could be used to judge the SCALR and SCALS
c
c The program use translation linear to least square refinement. Testing 
c shows that it works quite well. The function of this routine is similar
c to SUPOS but the r.m.s. is less and the orientation might be more dangerous.
c
c
c                                           by Guoguang Lu
C                                               27/01/1989
C
      COMMON/RMS/ RMS,SMEAN,NREF1,NREF
c	COMMON/SHIFS/ SCALR,SCALT
c	DATA SCALR/1./,SCALT/1./,IAT/0/
C
      DIMENSION X1(3,NATM), X2(3,NATM), A(3,3), T(3)
      DIMENSION VT(NATM*3,6),DIS(3,NATM)
      DIMENSION VVV(6,6), VV1(6)
C
      DIMENSION DA(3,3,3),DROH(3),DT(3),DETA(6)
      DIMENSION B1(3),B2(3),ROHOLD(3),ROH(3),TOLD(3)
      REAL EMT(3,3)
      EQUIVALENCE (DETA,DROH)
      EQUIVALENCE (DETA(4),DT)
      data emt /1.,0.,0.,0.,1.,0.,0.,0.,1./
C
      DEG=180./3.14159265359
      SCALR = 1.
      SCALT = 1.
      iat = 0.
      NREF=0
      epsilon = 1.
C
      CALL MTOHUBER(A,ROH,B1)
c
      CALL RTMOV(NATM,X1,A,T,DIS)
      CALL ARRPS(3,NATM,DIS,X2,DIS)
      RMS = DOSQ(NATM*3, DIS)
      RMS1 = RMS
        DO 5 I = 1, 3
         DO 5 J = 1, NATM
 5	CALL ARRGIVE(3,EMT(1,I),VT((J-1)*3+1,I+3))
C
C   Refine th1, th2, th3, tx, ty and tz
c
10	CONTINUE
c
      if ( (nref.gt.99) .or. (abs(epsilon).le.1.0e-5) ) goto 25
C	WRITE(6,*)
c	WRITE(6,*) 'Cycle of refinement:',NREF
C	WRITE(6,*) 'Matrix'
C	WRITE(6,20) A
C20	FORMAT(3F10.5)
C	WRITE(6,*) 'Eulerian Angle:',ROH
C	WRITE(6,*) 'Translation:   ',T
C	WRITE(6,*)  'SIGMA (distance^2)=',RMS
C
C compute the normal matrix
      DO I =1, 3
      CALL DRVROHTH(I,ROH,DA(1,1,I))
      CALL MATMULT(3,3,3,NATM,DA(1,1,I),X1,VT(1,I))
C
c	CALL ARRMC(3,1,B1,0.,B1)
c	B1(I) = 1.
c	DO J = 1, NATM
c	CALL ARRMC(3,1,B1, 1., VT((J-1)*3+1,I+3) )
c	END DO
C
      END DO
C
      CALL LSQEQ(3*NATM,6,VT,DIS,DETA,VVV,VV1)
C
C SHIFT = SHIFT * SCALE
      CALL ARRMC(3,1,DROH,-SCALR*deg,DROH)
      CALL ARRMC(3,1,DT,-SCALT,DT)
C	TYPE *, 'Shift ROH',DROH
C	TYPE *, 'Shift T  ',DT
C SAVE THE OLD ANGLE AND TRANSLATION VECTOR
      CALL ARRMC(3,1,ROH,1.,ROHOLD)
      CALL ARRMC(3,1,T,1.,TOLD)
C
      CALL ARRAD(3,1,ROH,DROH,ROH)
      CALL ARRAD(3,1,T,DT,T)
C
      CALL HUBER(ROH,A)
C
      CALL RTMOV(NATM,X1,A,T,DIS)
      CALL ARRPS(3,NATM,DIS,X2,DIS)
C
      RMS = DOSQ(NATM*3, DIS)
c	TYPE *, 'rms new and old ',rms, rms1 
C
      IF (RMS.LT.RMS1) THEN
      epsilon = (rms1-rms)/rms1
      NREF = NREF+1
      RMS1=RMS
      GOTO 10
      END IF
C
   25 continue
      CALL ARRMC(3,1,ROHOLD,1.,ROH)
      CALL ARRMC(3,1,TOLD,1.,T)
      CALL HUBER(ROH,A)
C
      CALL RTMOV(NATM,X1,A,T,DIS)
      CALL ARRPS(3,NATM,DIS,X2,DIS)
C
      ATM = NATM
      SMEAN = 0
      DO I=1, NATM
      D = VEM(3,DIS(1,I))
      SMEAN = SMEAN + D
      END DO
      SMEAN = SMEAN /ATM
C
c SMEAN = SIGMA (DISTANCE) / N
C           i        i
      RMS = SQRT( RMS1/ATM )
C
c RMS   = SQRT( SIGMA (DISTANCE^2) / N )
C           i              i
c 
c
      END
c
c
      SUBROUTINE REFRTFIN ( NATM, X1, X2, A, T, VT, DIS )
C	DIMENSION X1(3,NATM), X2(3,NATM), A(3,3), T(3)
C	DIMENSION VT(3*NATM,6),DIS(3,NATM)
C A subroutine to refine the superimpose matrix and vector from MOL1 to 
c MOL2 i.e. XYZ2(3*NATM) = A(3*3) * XYZ1(3*NATM) + T(3)
C Where    NATM represents the number of the atoms in each molecule
c                which will be superimposed.
c          XYZ1(3,NATM) represents the coordinates of MOL1
c          XYZ2(3,NATM) represents the coordinates of MOL2
C          A(3*3)       represents the input initial matrix and 
c                       output matrix.
c          T(3*3)       represents the input initial translation vector
c                       output vector.
C          VT(3*6*NATM)   is a 6*NATM-dimensional buffer array
C          DIS(3*NATM)  is a 3*NATM-dimensional buffer array
c
C	COMMON/SHIFS/ SCALR,SCALT
c 	DATA SCALR/1./,SCALT/1./
C Note: There should be a scale of the shifts. i.e. Shift= scale * shifts
c Here SCALR is the rotation shiftscale
c Here SCALT is the translation shiftscale
c The initial and suggested SCALR is 60.
c The initial and suggested SCALt is 1.
C
C	COMMON/RMS/ RMS,SMEAN,NREF,NREF1
C RMS is the final R.M.S between two molecules.
C SMEAN is the mean distance.
c NREF is the number of cycles of refinement
c NREF1 is not useful.
c They could be used to judge the SCALR and SCALS
c
c The routine uses translation linear to least square refinement. Testing
c shows that it works quite well. The function of this routine is similar
c to SUPOS but the r.m.s. is less and the orientation might be more dangerous.
c
c
c                                           by Guoguang Lu
C                                               27/01/1989
C
      COMMON/RMS/ RMS,SMEAN,NREF1,NREF
c	COMMON/SHIFS/ SCALR,SCALT
C
      DIMENSION X1(3,NATM), X2(3,NATM), A(3,3), T(3)
      DIMENSION VT(NATM*3,6),DIS(3,NATM)
      DIMENSION VVV(6,6), VV1(6)
C
      DIMENSION DA(3,3,3),DROH(3),DT(3),DETA(6)
      DIMENSION B1(3),B2(3),ROHOLD(3),ROH(3),TOLD(3)
      EQUIVALENCE (DETA,DROH)
      EQUIVALENCE (DETA(4),DT)
c...  data statements.  Separate declaration and init required for f2c
      DATA SCALR/1./,SCALT/1./,IAT/0/
C
      DEG=180./3.14159265359
      NREF=0
C
      CALL MTOHUBERARC(A,ROH,B1)
c
      CALL RTMOV(NATM,X1,A,T,DIS)
      CALL ARRPS(3,NATM,DIS,X2,DIS)
      RMS = DOSQ(NATM*3, DIS)
      RMS1 = RMS
C
C   Refine th1, th2, th3, tx, ty and tz
c
10	CONTINUE
c
C	WRITE(6,*)
c	WRITE(6,*) 'Cycle of refinement:',NREF
C	WRITE(6,*) 'Matrix'
C	WRITE(6,20) A
C20	FORMAT(3F10.5)
C	WRITE(6,*) 'Eulerian Angle:',ROH
C	WRITE(6,*) 'Translation:   ',T
C	WRITE(6,*)  'SIGMA (distance^2)=',RMS
C
C compute the normal matrix
      DO I =1, 3
      CALL DRVROHTHARC(I,ROH,DA(1,1,I))
      CALL MATMULT(3,3,3,NATM,DA(1,1,I),X1,VT(1,I))
C
      CALL ARRMC(3,1,B1,0.,B1)
      B1(I) = 1.
      DO J = 1, NATM
      CALL ARRMC(3,1,B1, 1., VT((J-1)*3+1,I+3) )
      END DO
C
      END DO
C
      CALL LSQEQ(3*NATM,6,VT,DIS,DETA,VVV,VV1)
C
C SHIFT = SHIFT * SCALE
C	CALL ARRMC(3,1,DROH,-SCALR*deg,DROH)
30	continue
      CALL ARRMC(3,1,DROH,-SCALR,DROH)
      CALL ARRMC(3,1,DT,-SCALT,DT)
C	TYPE *, 'scalr,scalt',scalr,scalt
C	TYPE *, 'Shift ROH',DROH
C	TYPE *, 'Shift T  ',DT
C SAVE THE OLD ANGLE AND TRANSLATION VECTOR
      CALL ARRMC(3,1,ROH,1.,ROHOLD)
      CALL ARRMC(3,1,T,1.,TOLD)
C
      CALL ARRAD(3,1,ROH,DROH,ROH)
      CALL ARRAD(3,1,T,DT,T)
C
      CALL HUBERARC(ROH,A)
C
      CALL RTMOV(NATM,X1,A,T,DIS)
      CALL ARRPS(3,NATM,DIS,X2,DIS)
C
      RMS = DOSQ(NATM*3, DIS)
C	TYPE *, 'refrt rms new and old ',rms, rms1
C
      IF (RMS.LT.RMS1) THEN
      NREF = NREF+1
      RMS1=RMS
      GOTO 10
      else if (scalr.gt.1.e-4) then
      CALL ARRMC(3,1,DROH,-1./SCALR,DROH)
      CALL ARRMC(3,1,DT,-1./SCALT,DT)
      scalr = scalr*0.5
      scalt = scalt*0.5
      CALL ARRMC(3,1,ROHOLD,1.,ROH)
      CALL ARRMC(3,1,TOLD,1.,T)
      CALL HUBERARC(ROH,A)
      goto 30
c
      END IF
C
      CALL ARRMC(3,1,ROHOLD,1.,ROH)
      CALL ARRMC(3,1,TOLD,1.,T)
      CALL HUBERARC(ROH,A)
C
      CALL RTMOV(NATM,X1,A,T,DIS)
      CALL ARRPS(3,NATM,DIS,X2,DIS)
C
      ATM = NATM
      SMEAN = 0
      DO I=1, NATM
      D = VEM(3,DIS(1,I))
      SMEAN = SMEAN + D
      END DO
      SMEAN = SMEAN /ATM
C
c SMEAN = SIGMA (DISTANCE) / N
C           i        i
      RMS = SQRT( RMS1/ATM )
C
c RMS   = SQRT( SIGMA (DISTANCE^2) / N )
C           i              i
c 
c
      END
c
c
      SUBROUTINE REFRTFIN1 ( NATM, X1, X2, A, T, VT, DIS, DIST )
C	DIMENSION X1(3,NATM), X2(3,NATM), A(3,3), T(3)
C	DIMENSION VT(NATM,6),DIS(3,NATM)
C A subroutine to refine the superimpose matrix and vector from MOL1 to 
c MOL2 i.e. XYZ2(3*NATM) = A(3*3) * XYZ1(3*NATM) + T(3)
C Where    NATM represents the number of the atoms in each molecule
c                which will be superimposed.
c          XYZ1(3,NATM) represents the coordinates of MOL1
c          XYZ2(3,NATM) represents the coordinates of MOL2
C          A(3*3)       represents the input initial matrix and 
c                       output matrix.
c          T(3*3)       represents the input initial translation vector
c                       output vector.
C          VT(6*NATM)   is a 6*NATM-dimensional buffer array
C          DIS(3*NATM)  is a 3*NATM-dimensional buffer array
C          DIST(NATM)  is a 3*NATM-dimensional buffer array
c
C	COMMON/SHIFS/ SCALR,SCALT
c 	DATA SCALR/1./,SCALT/1./
C Note: There should be a scale of the shifts. i.e. Shift= scale * shifts
c Here SCALR is the rotation shiftscale
c Here SCALT is the translation shiftscale
c The initial and suggested SCALR is 60.
c The initial and suggested SCALt is 1.
C
C	COMMON/RMS/ RMS,SMEAN,NREF,NREF1
C RMS is the final R.M.S between two molecules.
C SMEAN is the mean distance.
c NREF is the number of cycles of refinement
c NREF1 is not useful.
c They could be used to judge the SCALR and SCALS
c
c The routine uses translation linear to least square refinement. Testing
c shows that it works quite well. The function of this routine is similar
c to SUPOS but the r.m.s. is less and the orientation might be more dangerous.
c
c
c                                           by Guoguang Lu
C                                               27/01/1989
C
      COMMON/RMS/ RMS,SMEAN,NREF1,NREF
c	COMMON/SHIFS/ SCALR,SCALT
c	DATA SCALR/1./,SCALT/1./,IAT/0/
C
      DIMENSION X1(3,NATM), X2(3,NATM), A(3,3), T(3)
      DIMENSION VT(NATM,6),DIS(3,NATM),DIST(NATM)
      DIMENSION VVV(6,6), VV1(6)
C
      DIMENSION DA(3,3,3),DROH(3),DT(3),DETA(6)
      DIMENSION B1(3),B2(3),ROHOLD(3),ROH(3),TOLD(3)
      EQUIVALENCE (DETA,DROH)
      EQUIVALENCE (DETA(4),DT)
C
c	DEG=180./3.14159265359
      scalr = 1.0
      scalt = 1.0
      NREF=0
C
      CALL MTOHUBERarc(A,ROH,B1)
c
      CALL RTMOV(NATM,X1,A,T,DIS)
      CALL ARRPS(3,NATM,DIS,X2,DIS)
      DO I =1, NATM
      DIST(I) = VEM(3,DIS(1,I))
      END DO
      RMS = DOSQ(NATM*3, DIS)
      RMS1 = RMS
C
C   Refine th1, th2, th3, tx, ty and tz
c
10	CONTINUE
c
c	WRITE(6,*)
c	WRITE(6,*) 'Cycle of refinement:',NREF
c	WRITE(6,*) 'Matrix'
c	WRITE(6,20) A
c20	FORMAT(3F10.5)
c	WRITE(6,*) 'Eulerian Angle:',ROH
c	WRITE(6,*) 'Translation:   ',T
c	WRITE(6,*)  'SIGMA (distance^2)=',RMS
C
C compute the normal matrix
      DO I =1, 3
      CALL DRVROHTHarc(I,ROH,DA(1,1,I))
C
      DO J = 1, NATM
      CALL MATMULT(3,3,3,1,DA(1,1,I),X1(1,J),B1)
      VT(J,I) = POIMULT(3,3,DIS(1,J),B1) / DIST(J)
      VT(J,I+3) = DIS(I,J) / DIST(J)
      END DO
C
      END DO
      CALL LSQEQ(NATM,6,VT,DIST,DETA,VVV,VV1)
C	TYPE *, 'Shift ', DETA
C SHIFT = SHIFT * SCALe
c	CALL ARRMC(3,1,DROH,-SCALR*deg,DROH)
30	continue
C	TYPE *, 'SCALR, SCALT ', SCALR, SCALT
      CALL ARRMC(3,1,DROH,-SCALR,DROH)
      CALL ARRMC(3,1,DT,-SCALT,DT)
C	TYPE *, 'Shift ROH',DROH
C	TYPE *, 'Shift T  ',DT
C SAVE THE OLD ANGLE AND TRANSLATION VECTOR
      CALL ARRMC(3,1,ROH,1.,ROHOLD)
      CALL ARRMC(3,1,T,1.,TOLD)
C
      CALL ARRAD(3,1,ROH,DROH,ROH)
      CALL ARRAD(3,1,T,DT,T)
C
      CALL HUBERarc(ROH,A)
C
      CALL RTMOV(NATM,X1,A,T,DIS)
      CALL ARRPS(3,NATM,DIS,X2,DIS)
      DO I =1, NATM
      DIST(I) = VEM(3,DIS(1,I))
      END DO
C
      RMS = DOSQ(NATM*3, DIS)
C	TYPE *, 'rms new and old ',rms, rms1
C
      IF (RMS.LT.RMS1) THEN
      NREF = NREF+1
      RMS1=RMS
      GOTO 10
c	else IF (RMS.ge.0) THEN
      else IF (scalr.Gt.1.e-4) THEN
      CALL ARRMC(3,1,DROH,-1./SCALR,DROH)
      CALL ARRMC(3,1,DT,-1./SCALT,DT)
      scalr = scalr*0.5
      scalt = scalt*0.5
C
      CALL ARRMC(3,1,ROHOLD,1.,ROH)
      CALL ARRMC(3,1,TOLD,1.,T)
      CALL HUBERarc(ROH,A)
      goto 30
      END IF
C
      CALL ARRMC(3,1,ROHOLD,1.,ROH)
      CALL ARRMC(3,1,TOLD,1.,T)
      CALL HUBERarc(ROH,A)
      CALL RTMOV(NATM,X1,A,T,DIS)
      CALL ARRPS(3,NATM,DIS,X2,DIS)
      DO I =1, NATM
      DIST(I) = VEM(3,DIS(1,I))
      END DO

c	RMS = DOSQ(NATM*3, DIS)
c	TYPE *, 'rms repeat and old ',rms, rms1


c	TYPE *, 'ending rms ',rms
      ATM = NATM
      SMEAN = 0
      DO I=1, NATM
      D = VEM(3,DIS(1,I))
      SMEAN = SMEAN + D
      END DO
      SMEAN = SMEAN /ATM
c
c SMEAN = SIGMA (DISTANCE) / N
C           i        i
      RMS = SQRT( RMS1/ATM )
C
c RMS   = SQRT( SIGMA (DISTANCE^2) / N )
C           i              i
c 
c
      END
c
c
      function gg_res3to1(rin)
      character gg_res3to1*1,rin*3
      character*3 res(20),abr(20)
c...  data statements.  Separate declaration and init required for f2c
      data RES /'ALA','GLU','GLN','ASP','ASN','LEU','GLY'
     1 ,'LYS','SER','VAL','ARG','THR'
     2 ,'PRO','ILE','MET','PHE','TYR','CYS','TRP','HIS'/
      data ABR /'A  ','E  ','Q  ','D  ','N  ','L  ','G  '
     1 ,'K  ','S  ','V  ','R  ','T  '
     2 ,'P  ','I  ','M  ','F  ','Y  ','C  ','W  ','H  '/
c
      do i = 1, 20
       if (rin.eq.res(i)) then
        gg_res3to1 = abr(i)(1:1)
        return
       end if
      end do
      if (rin.eq.'   ') then
       gg_res3to1 = ' '
       return
      end if
c	write(6,*) 'Warning: no such residues ',rin
      gg_res3to1 = ' '
      end
c
c
      subroutine rotvsvec(vec,vkapa,a)
c a subroutine to give a matrix when you know there is a rotation
c which angle is vkapa along the vector vec as axis.
c vec is the input 3-dimensional vector where the rotation axis is
c vkapa is the rotating angle
c  a is the output 3*3 rotating matrix.
c
      real vec(3),a(3,3)
      external sind
c	WRITE(6,*)  'vec,vkapa'
C	WRITE(6,*)  vec,vkapa
c There is a strange bug on Alpha, ---- cosd(vkapa) has the wrong sign!!!!
c Instead I have to have cos(vkapa*pai/180.)
      pai = 3.141592654
c
      sinkapa = sind(vkapa)
      coskapa = cos(vkapa*pai/180.)
c
      a(1,1) = (1.-vec(1)*vec(1))*coskapa+vec(1)*vec(1)
      a(2,2) = (1.-vec(2)*vec(2))*coskapa+vec(2)*vec(2)
      a(3,3) = (1.-vec(3)*vec(3))*coskapa+vec(3)*vec(3)
c
      a(1,2) = vec(1)*vec(2)*(1.-coskapa) -vec(3)*sinkapa
      a(2,1) = vec(1)*vec(2)*(1.-coskapa) +vec(3)*sinkapa
c
      a(1,3) = vec(3)*vec(1)*(1.-coskapa) +vec(2)*sinkapa
      a(3,1) = vec(3)*vec(1)*(1.-coskapa) -vec(2)*sinkapa
c
      a(2,3) = vec(2)*vec(3)*(1.-coskapa) -vec(1)*sinkapa
      a(3,2) = vec(2)*vec(3)*(1.-coskapa) +vec(1)*sinkapa
c
      end
c
c
      SUBROUTINE RTMOV(NATM,XIN,A,T,XOUT)
C A routine to compute: XOUT = A * XIN + T
C        where   XIN is the input 3*NATM array
C                XOUT is the Output 3*NATM array
c                A(3,3) is a 3*3 rotation matrix
c                T(3) is a 3-dimensional translation vector.
c	DIMENSION XIN(3,NATM),XOUT(3,NATM),A(3,3),T(3)
C
      DIMENSION XIN(3,NATM),XOUT(3,NATM),A(3,3),T(3)
      CALL MATMULT(3,3,3,NATM,A,XIN,XOUT)
      DO I = 1, NATM
      CALL ARRAD(3,1,XOUT(1,I),T,XOUT(1,I))
      END DO
      END
c
c
      function screw(a,t)
c	dimension a(3,3),t(3),POL(3),xyz0(3)
C this is a routine to calculate how much the movement is along the
c the rotation axis.
c When there is a non-crystallographic rotation and translation
c x2 = A*x1 + t
c where A(3,3) is the 3*3 matrix
c 	t(3) is translation vector.
c the routine output POL(3) (psi, phi and kappa), go and xyz0(3)
c where psi (POL(1))  is the angle between the rotation axis and y axis
c	phi is the angle between X axis and the image of the rotation axis
c	kappa is the rotating angle.
c 	go is the screw amount along the axis direction
c		Guoguang 900927
c

      dimension a(3,3),t(3)
      dimension vec(3)
c
      call mtovec(a,vec,vkapa)
      if (vkapa.lt.1e-6) then
       screw = vem(3,t)
      end if
c
      if (vem(3,t).lt.(1.e-6)) then
       screw = 0.
       return
      end if
c vec is the direction of the axis vector.
      screw = poimult(3,3,t,vec)
      end
c
c
      subroutine sectime(sec,nhour,nmin,sres)
c input seconds and output the hour and minute and residue second
      nmin = 0
      nhour = 0
      if (sec.gt.60.) then
       smin = sec/60.
       nmin = int(smin)
       sres = sec - nmin*60
      end if
      if (nmin.gt.60) then
       nmin0 = nmin
       nhour = int(float(nmin)/60.)
       nmin = nmin0 - nhour*60 
      end if
      end
c
c
      subroutine seekaxis(a,t,cen1,cen2,POL,go,xyz0)
c	dimension a(3,3),t(3),POL(3),xyz0(3)
C this is a routine to calculate how much the movement is along the
c the rotation axis and where is the rotation axis from a rigid movement
c matrix.
c When there is a non-crystallographic rotation and translation
c x2 = A*x1 + t
c where A(3,3) is the 3*3 matrix
c 	t(3) is translation vector.
c cen(3), center of the two mass molecule
c the routine output POL(3) (psi, phi and kappa), go and xyz0(3)
c where psi (POL(1))  is the angle between the rotation axis and y axis
c	phi is the angle between X axis and the image of the rotation axis
c	kappa is the rotating angle.
c 	go is the screw amount along the axis direction
c		Guoguang 900927
c	xyz0(3) is one point where the axis pass by. This point should be 
c        a project point of CEN on the axis line, where CEN is center of
c	two input molecules (CEN1+CEN2)/2.
c

      dimension a(3,3),t(3),POL(3),xyz0(3),cen(3),cen1(3),cen2(3)
      dimension b1(3),b2(3),b3(3)
      dimension rot(3,3), rot1(3,3)
      dimension b4(3),b5(3),cenp1(3),cenp2(3),xyzp(3)
      dimension vec(3)
      logical*1 zax
      external tand
c
      call arrmc(3,1,t,0.,xyz0)
      call mtopolorz(a,POL,b1)
      if (abs(POL(3)).lt.1e-6) then
       go = vem(3,t)
      end if
c
      if (vem(3,t).lt.(1.e-6)) then
       go = 0.
       return
      end if
c
      call arrad(3,1,cen1,cen2,b4)
      call arrmc(3,1,b4,0.5,cen)
c
c	b1(2) = cosd(POL(1))	
c	b1(1) = sind(POL(1))*cosd(POL(2))	
c	b1(3) = - sind(POL(1))*sind(POL(2))
      call mtovec(a,vec,vkapa)
      call arrgive(3,vec,b3)
c b1 and vec is the direction of the axis vector.
      go = poimult(3,3,t,b3)
c now define a new coordinate system 
c X0 = rot * X'
c where Z' is the rotating axis.
c and X'is in the Z'/cen2-cen1 plane 
c Y' is perpendicular to X' and Z'
c so b1 is Z' vector.
      call arrps(3,1,cen2,cen1,b4)
      call veccrsmlt(b3,b4,b5)
      temp = vem(3,b5)
      if (temp.gt.(1.e-6)) then
       zax = .false.
 	 call arrmc(3,1,b5,1./temp,b2)
c here b2 is X' vector
       call veccrsmlt(b2,b3,b1)
c b3 is Y' vector
c
       call arrgive(3,b2,rot(1,2))
       call arrgive(3,b3,rot(1,3))
       call arrgive(3,b1,rot(1,1))
      else
       call elize(3,rot)
       zax = .true.
      end if
c
      call antiarr(3,3,rot,rot1)
c rot1 = rot**-1
c
c b1 is the translation in new coordinate system
11	call matmult(3,3,3,1,rot1,t,b1)
      if (abs(b1(3)-go).gt.1e-6) then
       if (zax) then
        rot(3,3) = -1.
        rot(2,2) = -1.
        rot1(3,3) = -1.
        rot1(2,2) = -1.
        zax = .false.
        goto 11
       else
        write(6,*) 'go ', go
        write(6,*) 't ',t
        write(6,*) 'b1 ',b1
        write(6,*) 'rot1'
        write(6,'(3f12.5)') rot1
        stop 'guoguang, you are wrong to get go'
       end if
      end if
c
c take cen as origin. take cen1 cen2 into new coordinates as cenp1,cenp2
      call arrps(3,1,cen1,cen,b1)
      call matmult(3,3,3,1,rot1,b1,cenp1)
      call arrps(3,1,cen2,cen,b2)
      call matmult(3,3,3,1,rot1,b2,cenp2)
c
      if (cenp1(2).gt.1.e-4) print *, 'Warning cenp1',cenp1
      if (cenp2(2).gt.1.e-4) print *, 'Warning cenp2',cenp2
c the passing point in new system
      xyzp(1) = 0.
      xyzp(3) = 0.
      if (abs(180.-vkapa).lt.0.001) then
       xyzp(2) = 0.
      else
       xyzp(2) = -cenp1(1)/tand(vkapa/2.)
      end if
c the passing point in old system
      call matmult(3,3,3,1,rot,xyzp,b4)
      call arrad(3,1,b4,cen,xyz0)
c
      end
c
c
      subroutine seekscrew(a,t,POL,go,xyz0)
c	dimension a(3,3),t(3),POL(3),xyz0(3)
C this is a routine to calculate how much the movement is along the
c the rotation axis.
c When there is a non-crystallographic rotation and translation
c x2 = A*x1 + t
c where A(3,3) is the 3*3 matrix
c 	t(3) is translation vector.
c the routine output POL(3) (psi, phi and kappa), go and xyz0(3)
c where psi (POL(1))  is the angle between the rotation axis and y axis
c	phi is the angle between X axis and the image of the rotation axis
c	kappa is the rotating angle.
c 	go is the screw amount along the axis direction
c		Guoguang 900927
c	xyz0(3) is one point where the axis pass by.
c

      dimension a(3,3),t(3),POL(3),xyz0(3)
      dimension b1(3),b2(3),b3(3)
      dimension rot(3,3), rot1(3,3),rot3(3,3),rot2(3,3)
      dimension b4(3)
      dimension vec(3)
      logical*1 zax
      external acosd, cosd
c
       call arrmc(3,1,t,0.,xyz0)
      call mtopolors(a,POL,b1)
      if (abs(POL(3)).lt.1e-6) then
       go = vem(3,t)
      end if
c
      if (vem(3,t).lt.(1.e-6)) then
       go = 0.
       return
      end if
c
c	b1(2) = cosd(POL(1))	
c	b1(1) = sind(POL(1))*cosd(POL(2))	
c	b1(3) = - sind(POL(1))*sind(POL(2))
      call mtovec(a,b1,vkapa)
      call arrmc(3,1,b1,1.,vec)
c b1 and vec is the direction of the axis vector.
      go = poimult(3,3,t,b1)
c
      b2(1) = 0.
      b2(2) = 0.
      b2(3) = 1.
c now decide a new coordinate system 
c X0 = rot * X'
c where Z' is the rotating axis.
c and X'is in the XY plane.
c so b1 is Z' vector.
      call veccrsmlt(b1,b2,b3)
      temp = vem(3,b3)
      if (temp.gt.(1.e-6)) then
       zax = .false.
 	 call arrmc(3,1,b3,1./temp,b2)
c here b2 is X' vector
       call veccrsmlt(b1,b2,b3)
c b3 is Y' vector
c
       call arrmc(3,1,b2,1.,rot(1,1))
       call arrmc(3,1,b3,1.,rot(1,2))
       call arrmc(3,1,b1,1.,rot(1,3))
      else
       call elize(3,rot)
       zax = .true.
      end if
c
      call antiarr(3,3,rot,rot1)
c rot1 = rot**-1
c
c b1 is the translation in new coordinate system
11	call matmult(3,3,3,1,rot1,t,b1)
      if (abs(b1(3)-go).gt.1e-6) then
       if (zax) then
        rot(3,3) = -1.
        rot(2,2) = -1.
        rot1(3,3) = -1.
        rot1(2,2) = -1.
        zax = .false.
       goto 11
       else
        write(6,*) 'go ', go
        write(6,*) 't ',t
        write(6,*) 'b1 ',b1
        write(6,*) 'rot1'
        write(6,'(3f12.5)') rot1
        stop 'guoguang, you are wrong to get go'
       end if
      end if
c
      if (vem(2,b1).lt.1e-6) then
c the axes pass by the origin.
       call arrmc(3,1,t,0.,xyz0)
       return
      end if
c
c set a new coordinate system to make axis x'' is the image of t in new
c x'y' plane. and X' = rot2 * X''
c 
c
      b1(3) = 0.
      call elize(3,rot2)
      call arrmc(3,1,b1, 1./vem(3,b1), rot2(1,1) )
      call veccrsmlt(rot2(1,3),rot2(1,1),rot2(1,2))
      call antiarr(3,3,rot2,rot3)
c x'' = rot3 * x' = rot3 * rot1 * x0
      call matmult(3,3,3,1,rot3,b1,b4)
      if ( abs(b4(1)-vem(2,b1)).gt.1e-3 .or. abs(b4(2)).gt.1e-3 .or
     1 . abs(b4(3)).gt.1e-3) THEN
       WRITE(6,*)  'b1',b1
       WRITE(6,*)  'b4',b4
       write(6,*) 'rot2'
       write(6,'(3f12.5)') rot2
       write(6,*) 'rot3'
       write(6,'(3f12.5)') rot3
       stop 'something wrong in rot2 or rot3'
      end if
      x0 = b4(1)
c
      b2(3) = 0.
      b2(1) = x0/2
      vkapa = POL(3)
      b2(2) = sqrt( (1+cosd(vkapa))/(1-cosd(vkapa)) ) * b2(1)
      if (vkapa.lt.0.) b2(2) = - b2(2)
c b2 is what you want in x'' system
c the mathematics is very simple, isn't it?
c
c check result
      call arrps(3,1,b2,b4,b3)
      ang = acosd( poimult(3,3,b2,b3)/(vem(3,b2)*vem(3,b3)) )
c	if (abs(temp).gt.1) WRITE(6,*)  temp
c	ang = asind(temp)
      if (abs(ang-abs(vkapa)).gt.1e-2) then
       WRITE(6,*)  'something wrong in seekscrew'
       WRITE(6,*) 'angle ',ang
       WRITE(6,*)  'b1 ',b1
       WRITE(6,*)  'b2 ',b2
       WRITE(6,*)  'b3 ',b3
       WRITE(6,*)  'b4 ',b4
      end if	
c now x0 = rot * x' = rot * rot2 * x'' = rot1 * x''
      call matmult(3,3,3,3,rot,rot2,rot1)
      call matmult(3,3,3,1,rot1,b2,xyz0)
c to check the result.
      ds1 = dstpl2(vec,xyz0,t,b1,b4)
      call arrmc(3,1,t,0.,b3)
      ds2 = dstpl2(vec,xyz0,b3,b2,b4)
      df = abs(ds2-ds1)
c
      if (ds1.gt.1e-4) then
        if (df/ds1.gt.1e-4) 
     + stop 'distance between 0-axis and t to axis is not same'
      
      else
        if (ds2.gt.2e-4) 
     + stop 'distance between 0-axis  and t to axis is not same'
      end if

      call arrps(3,1,b1,b2,b3)
      go1 = vem(3,b3)
      if (abs(abs(go)-go1) .gt. 1e-4) then
      WRITE(6,*)  'problem in go'
      WRITE(6,*)  'go ', go
      WRITE(6,*)  'screw', go1
      end if
c go is the real screw value along axis
      end
c
c
      subroutine seekscrewz(a,t,POL,go,xyz0)
c	dimension a(3,3),t(3),POL(3),xyz0(3)
C this is a routine to calculate how much the movement is along the
c the rotation axis.
c When there is a non-crystallographic rotation and translation
c x2 = A*x1 + t
c where A(3,3) is the 3*3 matrix
c 	t(3) is translation vector.
c the routine output POL(3) (psi, phi and kappa), go and xyz0(3)
c where psi (POL(1))  is the angle between the rotation axis and y axis
c	phi is the angle between X axis and the image of the rotation axis
c	kappa is the rotating angle.
c 	go is the screw amount along the axis direction
c		Guoguang 900927
c	xyz0(3) is one point where the axis pass by.
c

      dimension a(3,3),t(3),POL(3),xyz0(3)
      dimension b1(3),b2(3),b3(3)
      dimension rot(3,3), rot1(3,3),rot3(3,3),rot2(3,3)
      dimension b4(3)
      dimension vec(3)
      logical*1 zax
      external acosd, cosd
c
       call arrmc(3,1,t,0.,xyz0)
      call mtopolorz(a,POL,b1)
      if (abs(POL(3)).lt.1e-6) then
       go = vem(3,t)
      end if
c
      if (vem(3,t).lt.(1.e-6)) then
       go = 0.
       return
      end if
c
c	b1(2) = cosd(POL(1))	
c	b1(1) = sind(POL(1))*cosd(POL(2))	
c	b1(3) = - sind(POL(1))*sind(POL(2))
      call mtovec(a,b1,vkapa)
      call arrmc(3,1,b1,1.,vec)
c b1 and vec is the direction of the axis vector.
      go = poimult(3,3,t,b1)
c
      b2(1) = 0.
      b2(2) = 0.
      b2(3) = 1.
c now decide a new coordinate system 
c X0 = rot * X'
c where Z' is the rotating axis.
c and X'is in the XY plane.
c so b1 is Z' vector.
      call veccrsmlt(b1,b2,b3)
      temp = vem(3,b3)
      if (temp.gt.(1.e-6)) then
       zax = .false.
 	 call arrmc(3,1,b3,1./temp,b2)
c here b2 is X' vector
       call veccrsmlt(b1,b2,b3)
c b3 is Y' vector
c
       call arrmc(3,1,b2,1.,rot(1,1))
       call arrmc(3,1,b3,1.,rot(1,2))
       call arrmc(3,1,b1,1.,rot(1,3))
      else
       call elize(3,rot)
       zax = .true.
      end if
c
      call antiarr(3,3,rot,rot1)
c rot1 = rot**-1
c
c b1 is the translation in new coordinate system
11	call matmult(3,3,3,1,rot1,t,b1)
      if (abs(b1(3)-go).gt.1e-6) then
       if (zax) then
        rot(3,3) = -1.
        rot(2,2) = -1.
        rot1(3,3) = -1.
        rot1(2,2) = -1.
        zax = .false.
       goto 11
       else
        write(6,*) 'go ', go
        write(6,*) 't ',t
        write(6,*) 'b1 ',b1
        write(6,*) 'rot1'
        write(6,'(3f12.5)') rot1
        stop 'guoguang, you are wrong to get go'
       end if
      end if
c
      if (vem(2,b1).lt.1e-6) then
c the axes pass by the origin.
       call arrmc(3,1,t,0.,xyz0)
       return
      end if
c
c set a new coordinate system to make axis x'' is the image of t in new
c x'y' plane. and X' = rot2 * X''
c 
c
      b1(3) = 0.
      call elize(3,rot2)
      call arrmc(3,1,b1, 1./vem(3,b1), rot2(1,1) )
      call veccrsmlt(rot2(1,3),rot2(1,1),rot2(1,2))
      call antiarr(3,3,rot2,rot3)
c x'' = rot3 * x' = rot3 * rot1 * x0
      call matmult(3,3,3,1,rot3,b1,b4)
      if ( abs(b4(1)-vem(2,b1)).gt.1e-3 .or. abs(b4(2)).gt.1e-3 .or
     1 . abs(b4(3)).gt.1e-3) THEN
       WRITE(6,*)  'b1',b1
       WRITE(6,*)  'b4',b4
       write(6,*) 'rot2'
       write(6,'(3f12.5)') rot2
       write(6,*) 'rot3'
       write(6,'(3f12.5)') rot3
       stop 'something wrong in rot2 or rot3'
      end if
      x0 = b4(1)
c
      b2(3) = 0.
      b2(1) = x0/2
      vkapa = POL(3)
      b2(2) = sqrt( (1+cosd(vkapa))/(1-cosd(vkapa)) ) * b2(1)
      if (vkapa.lt.0.) b2(2) = - b2(2)
c b2 is what you want in x'' system
c the mathematics is very simple, isn't it?
c
c check result
      call arrps(3,1,b2,b4,b3)
      ang = acosd( poimult(3,3,b2,b3)/(vem(3,b2)*vem(3,b3)) )
c	if (abs(temp).gt.1) WRITE(6,*)  temp
c	ang = asind(temp)
      if (abs(ang-abs(vkapa)).gt.1e-2) then
       WRITE(6,*)  'something wrong in seekscrew'
       WRITE(6,*) 'angle ',ang
       WRITE(6,*)  'b1 ',b1
       WRITE(6,*)  'b2 ',b2
       WRITE(6,*)  'b3 ',b3
       WRITE(6,*)  'b4 ',b4
      end if	
c now x0 = rot * x' = rot * rot2 * x'' = rot1 * x''
      call matmult(3,3,3,3,rot,rot2,rot1)
      call matmult(3,3,3,1,rot1,b2,xyz0)
c to check the result.
      ds1 = dstpl2(vec,xyz0,t,b1,b4)
      call arrmc(3,1,t,0.,b3)
      ds2 = dstpl2(vec,xyz0,b3,b2,b4)
      df = abs(ds2-ds1)
c
      if (ds1.gt.1e-4) then
        if (df/ds1.gt.1e-4) 
     +  stop 'distance between 0-axis and t to axis is not same'
      else
        if (ds2.gt.2e-4) 
     +   stop 'distance between 0-axis  and t to axis is not same'
      end if

      call arrps(3,1,b1,b2,b3)
      go1 = vem(3,b3)
      if (abs(abs(go)-go1) .gt. 1e-4) then
      WRITE(6,*)  'problem in go'
      WRITE(6,*)  'go ', go
      WRITE(6,*)  'screw', go1
      end if
c go is the real screw value along axis
      end
c
c
      SUBROUTINE SQSTLZ(M,P,PS,U,V)
C a subroutine to calculate the orientation matrix and the two axes
c  from a symmetric matrix of an ellipse.         by LGG
c   M(2,2)  is the input symmetric matrix of the ellipse.
C   P(2,2)  is one of the two orientation matrices of the ellipse
C   PS(2,2) is the inverse matrix of the matrix PS(2,2)
C   U, V  are the lengths of the two axes of the ellipse.
C
      DIMENSION M(2,2),P(2,2),PS(2,2)
      REAL M,P,A,B,C,U,V,X1,X2
10	FORMAT(1X,4F10.3)
      A=1
      B=-(M(1,1)+M(2,2))
      C=M(1,1)*M(2,2)-M(2,1)*M(1,2)
      CALL SQUROOT(A,B,C,X1,X2)
      U=SQRT(1./X1)
      V=SQRT(1./X2)
      P(2,1)=(X1-M(1,1))/M(1,2)
      P(2,2)=(X2-M(1,1))/M(1,2)
      S1=SQRT(P(2,1)*P(2,1)+1)
      S2=SQRT(P(2,2)*P(2,2)+1)
      P(2,1)=P(2,1)/S1
      P(2,2)=P(2,2)/S2
      P(1,1)=1./S1
      P(1,2)=1./S2
      PS(1,2)=P(2,1)
      PS(2,1)=P(1,2)
      PS(1,1)=P(1,1)
      PS(2,2)=P(2,2)
      RETURN
      END
c
c
      subroutine spacegp(nunit,file,latnam,nsym,nrot,sym,nsp)
c In this subroutine the user gives a space group name, the subroutine reads the 
c ccp4 symmetry library, then it outputs the symmetry operations
c Input
c       nunit --- unit number of library file
c	file--- symmetry libary file name
c	latnam  ---  *14 character of a spacegroup name
c Output
c 	nsymm --- Total number of symmetry operations in this spacegroup
c	nrot ---- only rotation symmetric operation in this spacegroup
c	symm  ---- a 3*4*nsym symmmetric operation c matrix.
c	nsp --- space group number
c in common symdump,  if dump is true, dump all the information
c maxop is the maximum number of symmetry operations
c library format
c155 6  2    R32
c X,Y,Z * -Y,X-Y,Z * Y-X,-X,Z
c Y,X,-Z * -X,Y-X,-Z * X-Y,-Y,-Z
c X+1/3,Y+2/3,Z+2/3 * -Y+1/3,X-Y+2/3,Z+2/3 * Y-X+1/3,-X+2/3,Z+2/3  
c Y+1/3,X+2/3,-Z+2/3 * -X+1/3,Y-X+2/3,-Z+2/3 * X-Y+1/3,-Y+2/3,-Z+2/3
c X+2/3,Y+1/3,Z+1/3 * -Y+2/3,X-Y+1/3,Z+1/3 * Y-X+2/3,-X+1/3,Z+1/3
c Y+2/3,X+1/3,-Z+1/3 * -X+2/3,Y-X+1/3,-Z+1/3 * X-Y+2/3,-Y+1/3,-Z+1/3
      parameter (maxop=96)
      character file*80,latnam*14
      real sym(3,4,maxop)
      integer*4 nunit,nsp,nsym,nrot
      logical*1 dump
c	COMMON/SYMDUMP/ DUMP
c
      character*80 key,str*40
c...  data statements.  Separate declaration and init required for f2c
      data dump /.false./
c
c	WRITE(6,*)  'nsp',nsp
      nsym=0
      if (nunit.eq.0) nunit=25
      if (file(1:1).eq.' ') call spstrunct(file)
      if (file(1:1).eq.' ') file='SYMOP'
        call ccpdpn(nunit,file(1:lenstr(file)),'READONLY','F',0,0)
10	read(nunit,'(a)',end=220) key
      im = lenstr(latnam)
      in = lenstr(key)
      is = index(key(1:in),latnam(1:im))
      if (is.eq.0) goto 10
c	call redstrin(12,i1-1,key,npar)
      read(key(1:12),*,err=10,end=10) isp,iall,irot
c	read(key(12:80),'(a)') latnam
c	if (isp.ne.nsp) goto 10
      write(6,*) 'Space Group  >>> ',Latnam,isp
      do i = 1, iall
       read(nunit,'(a)') key
       il = lenstr(key)
       key(il+1:il+1)='*'
       iend = 0
       if (dump)  write(6,'(1x,a)') key(1:il)
20	 ist = iend+1
        iend = index(key(ist:il+1),'*')+IST-1
        str = '                                        '
        str(1:iend-ist) = key(ist:iend-1)
        NSYM=NSYM+1
        IF (MATSYM(SYM(1,1,nsym),STR,ICOL).EQ.0) GOTO 502
        WRITE(6,*) 'Error in symop after column ',ICOL
        write(6,*) key
        write(6,*) str
        stop 'Check you file SYMOP'
502	continue
       if (dump) then
        write(6,'(4f8.4)') ((sym(j1,j2,nsym),j2=1,4),j1=1,3)
       end if
       if (iend.lt.il+1) goto 20
       if (i.eq.irot) nrot = nsym
      end do
      write(6,45) nsym,nrot
45	format(1x,'Symmetric operation ----',6x,
     1 'Total: ',i3,6x,'Rotation:',i3)
      close (unit=nunit)
      return
210	WRITE(6,*)  'Error opening the SYMOP file ',FILE(1:LENSTR(FILE))
220	WRITE(6,*)  'Space group',latnam(1:im),
     +              ' was not found in SYMOP file'
      nsym = 0
      nrot = 0
      latnam = '      '
      STOP 'Check your SYMOP file'
      END
c
c
      SUBROUTINE SPSTRUNCT(STRING)
      CHARACTER*(*) STRING
      LENS = LEN(STRING)
      IL = LENSTR(STRING)
C
5	CONTINUE
      ISP = INDEX(STRING(1:IL),' ')
      IF (ISP.EQ.0.OR.ISP.GE.IL) RETURN
      STRING(ISP:IL-1) = STRING(ISP+1:IL)
      STRING(IL:IL) = ' '
      IL = IL - 1
      GOTO 5
      END
c
c
      SUBROUTINE SQUROOT(A,B,C,X1,X2)
c    Finds roots of a quadratic equation
      REAL A,B,C,X1,X2
      X1=(-B+SQRT(B*B-4*A*C))/(2*A)
      X2=(-B-SQRT(B*B-4*A*C))/(2*A)
      RETURN
      END
c
c
C a subroutine to calculate the statistics of a group of numbers
c	subroutine statis(n,value,vmean,rms,dmean)
c	dimension value(n)
c  n is the number of value in the group
c value is the n-dimensional array whose statistics are required.
c vmean is the output mean value of the array value(n)
c rms is the output r.m.s in value(n)
c dmean is the output mean deviation of value(n)
      subroutine statis(n,value,vmean,rms,dmean)
      dimension value(n)
      vmean = 0
      do i =1, n
      vmean = vmean + value(i)
      end do
      vmean = vmean/n
c
      rms = 0.
      dmean = 0.
      do i =1, n
      temp = value(i) - vmean
      dmean = dmean + abs(temp)
      rms = rms + temp*temp
      end do
      rms = sqrt(rms/n)
      dmean = dmean / n
      end 
c
c
      SUBROUTINE LGG_STRING(NUNIT,NCHA,TXT,NPAR,cont)
C A subroutine to write a character string to unit NUNIT using 
c the spaces in the string to split the string into a list of
c parameters.
c NUNIT is the unit number.
c NCHA is the length of the character string.
c TXT is the character string.
C NPAR is the number of parameters in this string.
c cont is a flag which represent if this is a continue string or start
c	.true. = continue
c	.false. = start
c
c
      CHARACTER*(*) TXT
      logical*1 cont 
      
      if (.not.cont) then
         NPAR = 0
         REWIND (NUNIT)
      end if
      
      jcha = 0
      if (ncha.gt.0) then
         if (txt(ncha:ncha).eq.'-') then
            jcha = 1
            txt(ncha:ncha)=' '
         end if
      end if
      
      I = 1
 10   IF (I.LE.NCHA-jcha) then 
         if (TXT(I:I).NE.' ') THEN
            J = I
 20         CONTINUE
            IF ((J+1).GT.NCHA) THEN
               WRITE(NUNIT,'(A)') TXT(I:J)
               NPAR = NPAR + 1
            else if (TXT(J+1:J+1).EQ.' ') then
               WRITE(NUNIT,'(A)') TXT(I:J)
               NPAR = NPAR + 1            
            ELSE  
               J = J + 1
               GOTO 20
            END IF
            I = J + 1
            GOTO 10
         ELSE IF (I.LT.NCHA-jcha) THEN
            I = I + 1
            GOTO 10
         END IF
      endif
      
      if (jcha.eq.1) TXT(NCHA:NCHA)='-'
C
c	if (jcha.eq.0) REWIND(NUNIT)
C
      END
c
c
      SUBROUTINE SUPIM ( NATM, XYZ1, XYZ2, A, T )
C	DIMENSION XYZ1(3,NATM), XYZ2(3,NATM), A(3,3), T(3)
C A subroutine to give the superimpose matrix and vector from MOL1 to 
c MOL2 i.e. XYZ2(3*NATM) = A(3*3) * XYZ1(3*NATM) + T(3)
C Where    NATM represents the number of the atoms in each molecule
c                which will be superimposed.
c          XYZ1(3,NATM) represents the coordinates of MOL1
c          XYZ2(3,NATM) represents the coordinates of MOL2
C          A(3*3)       represents the output matrix.
c          T(3*3)       represents the output vector.
c
C	COMMON/SHIFS/ SCALR,SCALT
c 	DATA SCALR/60./,SCALT/1./
C Note: There should be a scale of the shifts. i.e. Shift= scale * shifts
c Here SCALR is the rotation shiftscale
c Here SCALT is the translation shiftscale
c The initial and suggested SCALR is 60.
c The initial and suggested SCALt is 1.
C
C	COMMON/RMS/ RMS,SMEAN,NREF,NREF1
C RMS is the final R.M.S between two molecules.
C SMEAN is the mean distance.
c NREF is the number of cycles of rotation refinement
c NREF1 is not useful 
c They could be used to judge the SCALR and SCALS
c
C	COMMON/IAT/ IAT1,IAT2,IAT3,IAT
C	DATA SCALR/60./,IAT/0/
C      The routine uses the first three atoms to get the initial Eulerian 
C angle. IAT1, IAT2 and IAT3 are the indexes of the three selected atoms.
c If IAT = 0, these atoms will be selected inside this routine. If not,
c the three atoms will be defined outside the routine through the common
c block. If the program does not work or you want special atoms, change
c IAT1,IAT2,IAT3 in COMMON/IAT/ and set IAT to 0 in order to select 
c the three atoms in the main routine.
c
C NATM cannot be larger than NATM0 (currently 50000) or the parameter NATM0
C should be modified in this routine. 
c
c The program use translation linear to least square refinement. Testing 
c shows that it works quite well. The function of this routine is similar
c to SUPOS but the r.m.s. is less and the orientation might be more dangerous.
c
c
c                                           by Guoguang Lu
C                                               27/01/1989
C
C
      COMMON/RMS/ RMS,SMEAN,NREF,NREF1
      COMMON/SHIFS/ SCALR,SCALS
      COMMON/IAT/ IAT1,IAT2,IAT3,IAT
c	DATA SCALR/60./,SCALS/1./
C	DATA IAT/0/,IAT1/1/,IAT2/2/,IAT3/3/
      PARAMETER (NATM0 = 50000)
c   NATM0 is the largest number of atoms.
      DIMENSION XYZ1(3,NATM), XYZ2(3,NATM), A(3,3), T(3)
      DIMENSION X1(3,NATM0),X2(3,NATM0),VT(3*6*NATM0),DIS(3*NATM0)
      DIMENSION CEN1(3),CEN2(3)
      DIMENSION B1(3),B2(3),B3(3),T0(3)
c
C Number of atoms cannot be more than NATM0 or less than 3.
c
      IF (NATM.GT.NATM0) THEN
      WRITE(6,*) 'ERROR> Atom is more than ',NATM0
      STOP 
      ELSE IF (NATM.LT.3) THEN
      STOP 'ERROR> Atom is less than 3.'
      END IF

c Compute the initial matrix.
      CALL ORIEN(NATM,XYZ1,XYZ2,A)
C compute the  gravity of mol1 and mol2 to CEN1 and CEN2
      CALL AVERG(3,NATM,XYZ1,CEN1)
      CALL AVERG(3,NATM,XYZ2,CEN2)
C T is the initial translation vector
      CALL ARRPS(3,1,CEN2,CEN1,T0)
c Change the origin to CEN1.
      CALL TMOVE(3,NATM,XYZ1,CEN1,-1.,X1)
      CALL TMOVE(3,NATM,XYZ2,CEN1,-1.,X2)
c
C Refine th1, th2, th3, tx, ty and tz
c
      CALL REFRT( NATM, X1, X2, A, T0, VT, DIS)
c
      WRITE(6,*)
      WRITE(6,*) 'R.M.S.'
      WRITE(6,*) '       natm'
      WRITE(6,*) 'SQRT( SIGMA (Distance(i)^2)/natm ) = ',RMS
      WRITE(6,*) '       i=1'
C
      WRITE(6,*)
      WRITE(6,*) 'Mean Distance:'
      WRITE(6,*) ' natm'
      WRITE(6,*) 'SIGMA (Distance(i)) /natm = ',SMEAN
      WRITE(6,*) ' i=1'
C
      WRITE(6,'(A)') '0'
      WRITE(6,'(A)') '0'
      WRITE(6,*) 'Mol1 is superposed to Mol2.'
      WRITE(6,*) 'The matrix and the vector are:'
      WRITE(6,*) 
      WRITE(6,1) (A(1,J),J=1,3),CEN1(1),T0(1)
      WRITE(6,2) (A(2,J),J=1,3),CEN1(2),T0(2)
      WRITE(6,1) (A(3,J),J=1,3),CEN1(3),T0(3)
1	FORMAT(1X,'      (',3F10.6,' )   (     ',F8.3,' )   ('
     1 ,F8.3,' )')
2	FORMAT(1X,' X2 = (',3F10.6,' ) * ( X1 -',F8.3,' ) + ('
     1 ,F8.3,' )')
C
      CALL MATMULT(3,3,3,1,A,CEN1,B3)
      CALL ARRPS(3,1,T0,B3,B3)
      CALL ARRAD(3,1,CEN1,B3,T)
C
      WRITE(6,'(A)') '0'
      WRITE(6,'(A)') '0'
      WRITE(6,*) 
      WRITE(6,4) (A(1,J),J=1,3),T(1)
      WRITE(6,5) (A(2,J),J=1,3),T(2)
      WRITE(6,4) (A(3,J),J=1,3),T(3)
4	FORMAT(1X,'      (',3F10.6,' )   (    )   (',F8.3,' )')
5	FORMAT(1X,' X2 = (',3F10.6,' ) * ( X1 ) + (',F8.3,' )')
C
      END 
c
c
      SUBROUTINE SUPRIMP ( NATM, XYZ1, XYZ2, A, T )
C	DIMENSION XYZ1(3,NATM), XYZ2(3,NATM), A(3,3), T(3)
C A subroutine to give the superimpose matrix and vector from MOL1 to 
c MOL2 i.e. XYZ2(3*NATM) = A(3*3) * XYZ1(3*NATM) + T(3)
C Where    NATM represents the number of the atoms in each molecule
c                which will be superimposed.
c          XYZ1(3,NATM) represents the coordinates of MOL1
c          XYZ2(3,NATM) represents the coordinates of MOL2
C          A(3*3)       represents the output matrix.
c          T(3*3)       represents the output vector.
c
C NATM cannot be larger than NATM0 (currently 50000) or the parameter NATM0
C should be modified in this routine. 
c
C	COMMON/SHIFS/ SCALR,SCALT
c 	DATA SCALR/1./,SCALT/1./
C Note: There should be a scale of the shifts. i.e. Shift= scale * shifts
c Here SCALR is the rotation shiftscale
c Here SCALT is the translation shiftscale
c The initial and suggested SCALR is 1.
c The initial and suggested SCALt is 1.
C
C	COMMON/RMS/ RMS,SMEAN,NREF,NREF1
C RMS is the final R.M.S between two molecules.
C SMEAN is the mean distance.
c NREF is the number of cycles of rotation refinement
c NREF1 is not used by this routine.
c They could be used to judge the SCALR and SCALS
c
C	COMMON/IAT/ IAT1,IAT2,IAT3,IAT
C	DATA SCALR/1./,IAT/0/
C      The routine use the first three atoms to get the initial Eulerian 
C angle. IAT1, IAT2 and IAT3 are the indexes of the three select atoms.
c If IAT = 0, these atoms will be selected inside this routine. If not,
c the three atoms will be defined outside the routine through the common
c block. If the program does not work or you want special atoms, change 
c IAT1,IAT2,IAT3 in COMMON/IAT/ and set IAT to 0 in order to select 
c the three atoms in the main routine. Subroutine ORIEN performs this
c computation.
c 
c     Then the program use the subroutine REFORN to refine the orientation
c by vector method. The refinement equation is same with subroutine SUPOS.
C
c The program use translation linear and Eulerian non-linear least square 
c refinement to refine both th1,th2,th3 and tx,ty,tz. Testing shows
c that it can give very low r.m.s., much less than any other method 
c in most cases, however it is a little expensive in computer time.
c
c                                           by Guoguang Lu
C                                               27/01/1989
C
      COMMON/RMS/ RMS,SMEAN,NREF,NREF1
      COMMON/SHIFS/ SCALR,SCALS
      COMMON/IAT/ IAT1,IAT2,IAT3,IAT
c	DATA SCALR/1./,SCALS/1./,IAT/0/,IAT1/1/,IAT2/2/,IAT3/3/
      PARAMETER (NATM0 = 50000)
c   NATM0 is the largest number of atoms.
      DIMENSION XYZ1(3,NATM), XYZ2(3,NATM), A(3,3), T(3)
      DIMENSION X1(3,NATM0),X2(3,NATM0),VT(3*6*NATM0),DIS(3*NATM0)
      DIMENSION CEN1(3),CEN2(3)
      DIMENSION B1(3),B2(3),B3(3),T0(3)
c
c	DEG=180./3.14159265359
c
C Number of atoms cannot be more than NATM0 or less than 3.
c
      IF (NATM.GT.50000) THEN
      WRITE(6,*) 'ERROR> Atom is more than ',NATM0
      STOP 
      ELSE IF (NATM.LT.3) THEN
      STOP 'ERROR> Atom is less than 3.'
      END IF

c Compute the initial matrix.
      CALL ORIEN(NATM,XYZ1,XYZ2,A)
c	NREF1 = NREF
C
C Refine the orientation by vector method.
C
      CALL REFORN ( NATM, XYZ1, XYZ2, A, VT, DIS)
C
C compute the  gravity of mol1 and mol2 to CEN1 and CEN2
      CALL AVERG(3,NATM,XYZ1,CEN1)
      CALL AVERG(3,NATM,XYZ2,CEN2)
C T is the initial translation vector
      CALL ARRPS(3,1,CEN2,CEN1,T0)
c Change the origin to CEN1.
      CALL TMOVE(3,NATM,XYZ1,CEN1,-1.,X1)
      CALL TMOVE(3,NATM,XYZ2,CEN1,-1.,X2)
c
C Refine th1, th2, th3, tx, ty and tz
c
      CALL REFRT( NATM, X1, X2, A, T0, VT, DIS)
c
      WRITE(6,*)
      WRITE(6,*) 'R.M.S.'
      WRITE(6,*) '       natm'
      WRITE(6,*) 'SQRT( SIGMA (Distance(i)^2)/natm ) = ',RMS
      WRITE(6,*) '       i=1'
C
      WRITE(6,*)
      WRITE(6,*) 'Mean Distance:'
      WRITE(6,*) ' natm'
      WRITE(6,*) 'SIGMA (Distance(i)) /natm = ',SMEAN
      WRITE(6,*) ' i=1'
C
      WRITE(6,*)
      WRITE(6,*) 'Mol1 is superposed to Mol2.'
      WRITE(6,*) 'The matrix and the vector are:'
      WRITE(6,*) 
      WRITE(6,1) (A(1,J),J=1,3),CEN1(1),T0(1)
      WRITE(6,2) (A(2,J),J=1,3),CEN1(2),T0(2)
      WRITE(6,1) (A(3,J),J=1,3),CEN1(3),T0(3)
1	FORMAT(1X,'      (',3F10.6,' )   (     ',f10.5,' )   ('
     1 ,f10.5,' )')
2	FORMAT(1X,' X2 = (',3F10.6,' ) * ( X1 -',f10.5,' ) + ('
     1 ,f10.5,' )')
C
      CALL MATMULT(3,3,3,1,A,CEN1,B3)
      CALL ARRPS(3,1,T0,B3,B3)
      CALL ARRAD(3,1,CEN1,B3,T)
C
      WRITE(6,*) 
      WRITE(6,*) 
      WRITE(6,4) (A(1,J),J=1,3),T(1)
      WRITE(6,5) (A(2,J),J=1,3),T(2)
      WRITE(6,4) (A(3,J),J=1,3),T(3)
4	FORMAT(1X,'      (',3F10.6,' )   (    )   (',f10.5,' )')
5	FORMAT(1X,' X2 = (',3F10.6,' ) * ( X1 ) + (',f10.5,' )')
C
      END 
c
c
      SUBROUTINE SUPRIMPFIN ( NATM, XYZ1, XYZ2, A, T )
C	DIMENSION XYZ1(3,NATM), XYZ2(3,NATM), A(3,3), T(3)
C A subroutine to give the superimpose matrix and vector from MOL1 to 
c MOL2 i.e. XYZ2(3*NATM) = A(3*3) * XYZ1(3*NATM) + T(3)
C Where    NATM represents the number of the atoms in each molecule
c                which will be superimposed.
c          XYZ1(3,NATM) represents the coordinates of MOL1
c          XYZ2(3,NATM) represents the coordinates of MOL2
C          A(3*3)       represents the output matrix.
c          T(3*3)       represents the output vector.
c
C NATM cannot be larger than NATM0 (currently 50000) or the parameter NATM0
C should be modified in this routine. 
c
C	COMMON/SHIFS/ SCALR,SCALT
c 	DATA SCALR/1./,SCALT/1./
C Note: There should be a scale of the shifts. i.e. Shift= scale * shifts
c Here SCALR is the rotation shiftscale
c Here SCALT is the translation shiftscale
c The initial and suggested SCALR is 1.
c The initial and suggested SCALt is 1.
C
C	COMMON/RMS/ RMS,SMEAN,NREF,NREF1
C RMS is the final R.M.S between two molecules.
C SMEAN is the mean distance.
c NREF is the number of the cycles of the rotation refinement
c NREF1 is not used by this routine.
c They could be used to judge the SCALR,and SCALS
c
C	COMMON/IAT/ IAT1,IAT2,IAT3,IAT
C	DATA SCALR/1./,IAT/0/
C      The routine use the first three atoms to get the initial Eulerian 
C angle. IAT1, IAT2 and IAT3 are the indexes of the three select atoms.
c If IAT = 0, these atoms will be selected inside this routine. If not,
c the three atoms will be defined outside the routine through the common
c block. If the program does not work or you want special atoms, change
c IAT1,IAT2,IAT3 in COMMON/IAT/ and give the IAT to 0 in order to select 
c the three atoms in the main routine. Subroutine ORIEN is to proceed this
c computing
c 
c     Then the program uses the subroutine REFORN to refine the orientation
c by vector method. The refinement equation is same with subroutine SUPOS.
C
c The program use translation linear and Eulerian non-linear least square 
c refinement to refine both th1,th2,th3 and tx,ty,tz. Testing shows
c that it can give very low r.m.s., much less than any other method 
c in most cases, however it is a little expensive in computer time.
c
c                                           by Guoguang Lu
C                                               27/01/1989
C
      COMMON/RMS/ RMS,SMEAN,NREF,NREF1
      COMMON/SHIFS/ SCALR,SCALS
      COMMON/IAT/ IAT1,IAT2,IAT3,IAT
c	DATA SCALR/1./,SCALS/1./,IAT/0/,IAT1/1/,IAT2/2/,IAT3/3/
      PARAMETER (NATM0 = 50000)
c   NATM0 is the maximum number of atoms.
      DIMENSION XYZ1(3,NATM), XYZ2(3,NATM), A(3,3), T(3)
      DIMENSION X1(3,NATM0),X2(3,NATM0),VT(3*6*NATM0),DIS(3*NATM0)
      DIMENSION CEN1(3),CEN2(3)
      DIMENSION B1(3),B2(3),B3(3),T0(3)
c
c	DEG=180./3.14159265359
c
C Number of atoms cannot be more than NATM0 or less than 3.
c
      IF (NATM.GT.50000) THEN
      WRITE(6,*) 'ERROR> Atom is more than ',NATM0
      STOP 
      ELSE IF (NATM.LT.3) THEN
      STOP 'ERROR> Atom is less than 3.'
      END IF

c Compute the initial matrix.
      CALL ORIEN(NATM,XYZ1,XYZ2,A)
c	NREF1 = NREF
C
C Refine the orientation by vector method.
cc
      CALL REFORNFIN ( NATM, XYZ1, XYZ2, A, VT, DIS)
      NREFORN = NREF
C
C compute the  gravity of mol1 and mol2 to CEN1 and CEN2
      CALL AVERG(3,NATM,XYZ1,CEN1)
      CALL AVERG(3,NATM,XYZ2,CEN2)
C T is the initial translation vector
      CALL ARRPS(3,1,CEN2,CEN1,T0)
c Change the origin to CEN1.
      CALL TMOVE(3,NATM,XYZ1,CEN1,-1.,X1)
      CALL TMOVE(3,NATM,XYZ2,CEN1,-1.,X2)
c
C Refine th1, th2, th3, tx, ty and tz
c
      CALL REFRTFIN ( NATM, X1, X2, A, T0, VT, DIS )
      NREFRT = NREF1
      CALL REFRTFIN1 ( NATM, X1, X2, A, T0, VT, DIS, VT(1+2*6*NATM) )
      NREFRT1 = NREF1
c
      WRITE(6,*)  'Final fit cycle>>>', NREFORN, NREFRT, NREFRT1
      WRITE(6,*)
      WRITE(6,*) 'R.M.S.'
      WRITE(6,*) '       natm'
      WRITE(6,*) 'SQRT( SIGMA (Distance(i)^2)/natm ) = ',RMS
      WRITE(6,*) '       i=1'
C
      WRITE(6,*)
      WRITE(6,*) 'Mean Distance:'
      WRITE(6,*) ' natm'
      WRITE(6,*) 'SIGMA (Distance(i)) /natm = ',SMEAN
      WRITE(6,*) ' i=1'
C
      WRITE(6,*)
      WRITE(6,*) 'Mol1 is superposed to Mol2.'
      WRITE(6,*) 'The matrix and the vector are:'
      WRITE(6,*) 
      WRITE(6,1) (A(1,J),J=1,3),CEN1(1),T0(1)
      WRITE(6,2) (A(2,J),J=1,3),CEN1(2),T0(2)
      WRITE(6,1) (A(3,J),J=1,3),CEN1(3),T0(3)
1	FORMAT(1X,'      (',3F10.6,' )   (     ',f10.5,' )   ('
     1 ,f10.5,' )')
2	FORMAT(1X,' X2 = (',3F10.6,' ) * ( X1 -',f10.5,' ) + ('
     1 ,f10.5,' )')
C
      CALL MATMULT(3,3,3,1,A,CEN1,B3)
      CALL ARRPS(3,1,T0,B3,B3)
      CALL ARRAD(3,1,CEN1,B3,T)
C
      WRITE(6,*) 
      WRITE(6,*) 
      WRITE(6,4) (A(1,J),J=1,3),T(1)
      WRITE(6,5) (A(2,J),J=1,3),T(2)
      WRITE(6,4) (A(3,J),J=1,3),T(3)
4	FORMAT(1X,'      (',3F10.6,' )   (    )   (',f10.5,' )')
5	FORMAT(1X,' X2 = (',3F10.6,' ) * ( X1 ) + (',f10.5,' )')
C
      END 
c
c
      subroutine lgg_symlib(nunit,file,nsp,nsym,nrot,sym,latnam)
c In this subroutine the user gives a space group number, the subroutine reads the 
c ccp4 symmetry library, then it outputs the symmetry operations
c Input
c	file--- symmetry library file name
c	nsp --- space group number
c Output
c 	nsymm --- Total number of symmetry operations in this spacegroup
c	nrot ---- only rotation symmetric operations in this spacegroup
c	symm  ---- a 3*4*nsym symmetric operation c matrix.
c	latnam  ---  *14 character spacegroup name
c in common symdump,  if dump is true, dump all the information
c maxop is the maximum number of symmetry operations
c libary format
c155 6  2    R32
c X,Y,Z * -Y,X-Y,Z * Y-X,-X,Z
c Y,X,-Z * -X,Y-X,-Z * X-Y,-Y,-Z
c X+1/3,Y+2/3,Z+2/3 * -Y+1/3,X-Y+2/3,Z+2/3 * Y-X+1/3,-X+2/3,Z+2/3  
c Y+1/3,X+2/3,-Z+2/3 * -X+1/3,Y-X+2/3,-Z+2/3 * X-Y+1/3,-Y+2/3,-Z+2/3
c X+2/3,Y+1/3,Z+1/3 * -Y+2/3,X-Y+1/3,Z+1/3 * Y-X+2/3,-X+1/3,Z+1/3
c Y+2/3,X+1/3,-Z+1/3 * -X+2/3,Y-X+1/3,-Z+1/3 * X-Y+2/3,-Y+1/3,-Z+1/3
      parameter (maxop=96)
      character file*80,latnam*14
      real sym(3,4,maxop)
      integer*4 nunit,nsp,nsym,nrot
      logical*1 dump
c	COMMON/SYMDUMP/ DUMP
c
      character*80 key,str*40
c
c...  data statements.  Separate declaration and init required for f2c
      data dump /.false./
c
c	WRITE(6,*)  'nsp',nsp
      nsym=0
      if (nunit.eq.0) nunit=25
      if (file(1:1).eq.' ') call spstrunct(file)
      if (file(1:1).eq.' ') file='SYMOP'
        call ccpdpn(nunit,file(1:lenstr(file)),'READONLY','F',0,0)
10	read(nunit,'(a)',end=220) key
c	il = lenstr(key)
c	i1 = index(key(1:il),'PG')
c	if (i1.eq.0) goto 10
c	call redstrin(12,i1-1,key,npar)
      read(key(1:12),*,err=10,end=10) isp,iall,irot
c	WRITE(6,*) key(1:lenstr(key)),isp,nsp,latnam
      if (isp.ne.nsp) goto 10
c	i1 = index(key(1:il),'PG')
      read(key(12:80),'(a)') latnam
c	call spstrunct(latnam)
      if (latnam(1:1).eq.' ') latnam(1:13) = latnam(2:14)
      if (isp.ne.nsp) goto 10
      write(6,*) 'Space Group  >>> ',Latnam,isp
      do i = 1, iall
       read(nunit,'(a)') key
       il = lenstr(key)
       key(il+1:il+1)='*'
       iend = 0
       if (dump)  write(6,'(1x,a)') key(1:il)
20	 ist = iend+1
        iend = index(key(ist:il+1),'*')+IST-1
        str = '                                        '
        str(1:iend-ist) = key(ist:iend-1)
        NSYM=NSYM+1
        IF (MATSYM(SYM(1,1,nsym),STR,ICOL).EQ.0) GOTO 502
        WRITE(6,*) 'Error in symop after colunm ',ICOL
        write(6,*) key
        write(6,*) str
        stop 'Check you file SYMOP'
502	continue
       if (dump) then
        write(6,'(4f8.4)') ((sym(j1,j2,nsym),j2=1,4),j1=1,3)
       end if
       if (iend.lt.il+1) goto 20
       if (i.eq.irot) nrot = nsym
      end do
      write(6,45) nsym,nrot
45	format(1x,'Symmetric operation ----',6x,
     1 'Total: ',i3,6x,'Rotation:',i3)
      close (unit=nunit)
      return
210	WRITE(6,*)  'Error opening the SYMOP file ',FILE(1:LENSTR(FILE))
220	WRITE(6,*)  'Space group',nsp,' was not found in SYMOP file'
      nsym = 0
      nrot = 0
      latnam = '      '
      STOP 'Check your SYMOP file'
      END
c
c
      SUBROUTINE TMOVE(M,NATM,XIN,T,CONST,XOUT)
C  XOUT = XIN + CONST * T
C where    XIN is  input M*NATM-dimensional array.
C          XOUT N is  output M*NATM-dimensional array.
c          T is a M-dimensional vector.
c          CONST is a constant.
C	DIMENSION XIN(M,NATM),XOUT(M,NATM),T(M),B(100)
c
      DIMENSION XIN(M,NATM),XOUT(M,NATM),T(M),B(100)
      CALL ARRMC(M,1,T,CONST,B)
      DO I = 1, NATM
      CALL ARRAD(M,1,XIN(1,I),B,XOUT(1,I))
      END DO
      END
c
c
c------------------------------------------------------------
      subroutine up(txt,len)
c
c convert character string to upper case
c
      character*(*) txt
      character*80 save
      character*26 ualpha,lalpha
      data ualpha /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data lalpha /'abcdefghijklmnopqrstuvwxyz'/

      do 9000 i=1,len
        if((txt(i:i).ge.'a').and.(txt(i:i).le.'z')) then
          match = index(lalpha,txt(i:i))
          save(i:i) = ualpha(match:match)
        else
          save(i:i) = txt(i:i)
        end if
c      end do
9000	continue

      txt = save
      return
      end
c
c
      SUBROUTINE VECCRSMLT(A1,A2,A)
      DIMENSION A1(3),A2(3),A(3)
      A(1)=VLDIM2(A1(2),A1(3),A2(2),A2(3))
      A(2)=VLDIM2(A1(3),A1(1),A2(3),A2(1))
      A(3)=VLDIM2(A1(1),A1(2),A2(1),A2(2))
      END
c
c
      FUNCTION VEM(N,V)
      DIMENSION V(N)
      C=0
      DO 10 I=1,N
10	C=C+V(I)*V(I)
      VEM=SQRT(C)
      RETURN
      END
c
c
      FUNCTION VLDIM3(AT)
C A function to calculate the modulus of a 3*3-dimension matrix.
c VLDIM3 = | AT3*3 |
      dimension at(3,3)
      dimension B1(3)
C
      call veccrsmlt(at(1,1),at(1,2),b1)
      vldim3 = poimult(3,3,b1,at(1,3))
      end
c
c
      FUNCTION VLDIM2(A11,A12,A21,A22)
      VLDIM2=A11*A22-A21*A12
      END
c
c
        REAL FUNCTION COSD(ANGINDEG)
        REAL ANGINDEG
        COSD = COS(ANGINDEG*3.1415926/180.)
        END
c
c
        REAL FUNCTION SIND(ANGINDEG)
        REAL ANGINDEG
        SIND = SIN(ANGINDEG*3.1415926/180.)
        END
c
c
        REAL FUNCTION TAND(ANGINDEG)
        REAL ANGINDEG
        TAND = TAN(ANGINDEG*3.1415926/180.)
        END
c
c
        REAL FUNCTION ACOSD(VAL)
        REAL VAL
        ACOSD = ACOS(VAL)*180./3.1415926
        END
c
c
        REAL FUNCTION ASIND(VAL)
        REAL VAL
        ASIND = ASIN(VAL)*180./3.1415926
        END
c
c
        REAL FUNCTION ATAND(VAL)
        REAL VAL
        ATAND = ATAN(VAL)*180./3.1415926
        END
