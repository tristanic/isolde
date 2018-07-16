C
C     sorting_main.f: sorting function library
C     Copyright (C) 2001  Garib Murshudov
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
                                !
                                !   A module to arrange memory for sorting
                                !   records (allocate, reallocate, deallocate)
      module sorting_commons
      implicit none
      integer :: nmax_rec = 2000000
      INTEGER NKEYS_S,NRECORD_S,KEY_ADDRESS,RECORD_ADDRESS
      integer  NRECORD_NOW,NRECORD_RETURN,INDEX_ADDRESS
      integer NRECORD_IN_RUNS,NRECORD_IN_THIS_RUN
      integer  NMAX_RECORD
      integer nkeys_l,nrec_l
      INTEGER SAVE_KEYS(5),SAVE_RECORD(200),ASCEND_DESCEND(5)
      REAL, allocatable :: keys_mem(:)
      real, allocatable :: array_mem(:)
      real, allocatable :: index_mem(:)


      contains

      subroutine sorting_allocate_this()


      if(nmax_rec.le.0)nmax_rec = 2000000

      allocate(keys_mem(nkeys_l*nmax_rec))
      allocate(array_mem(max(1,(nrec_l-nkeys_l)*nmax_rec)))
      allocate(index_mem(nmax_rec))

      end subroutine

      subroutine sorting_reallocate_this()

      integer nk1,nr1,ni1
      real, allocatable :: temp_keys(:)
      real, allocatable :: temp_recs(:)
      real, allocatable :: temp_inds(:)

      nk1 = nkeys_l*nmax_rec
      nr1 = max(1,(nrec_l-nkeys_l)*nmax_rec)
      ni1 = nmax_rec

      allocate(temp_keys(nk1))
      allocate(temp_recs(nr1))
      allocate(temp_inds(ni1))

      temp_keys(1:nk1) = keys_mem(1:nk1)
      temp_recs(1:nr1) = array_mem(1:nr1)
      temp_inds(1:ni1) = index_mem(1:ni1)

      nmax_rec = nint(nmax_rec*1.5)
      deallocate(keys_mem)
      deallocate(array_mem)
      deallocate(index_mem)

      allocate(keys_mem(nkeys_l*nmax_rec))
      allocate(array_mem(max(1,(nrec_l-nkeys_l)*nmax_rec)))
      allocate(index_mem(nmax_rec))
      
      keys_mem(1:nk1) = temp_keys(1:nk1)
      array_mem(1:nr1) = temp_recs(1:nr1)
      index_mem(1:ni1) = temp_inds(1:ni1)

      deallocate(temp_keys)
      deallocate(temp_recs)
      deallocate(temp_inds)


      return
      end subroutine

      subroutine sorting_deallocate_this()

      deallocate(keys_mem)
      deallocate(array_mem)
      deallocate(index_mem)

      end subroutine

      end module sorting_commons
                                !
                                !  now sorting bit
      INTEGER FUNCTION SRTBEG(NKEYS,KEYBUF,NRECL,ISS)
      use sorting_commons
C
      IMPLICIT NONE
      INTEGER NKEYS,NRECL,KEYBUF(*),ISS
C
      INTEGER IK,IREC,J,IRECORD_POSITIONS
      logical found
C
C----Find out how many records could be sorted. We need memory for
C----NRECL = NKEYS + NRECORDS number of fields and one field for
C----index of records which is going to be sorted.
      nkeys_l = nkeys
      nrec_l = nrecl/4
C
C---allocate enough memory (if it is not enough then we will reallocate)
      call sorting_allocate_this()

C
C----Save positions of keys and remaining records in an array
C----It would be easier if keybuf would have only info about 
C----ascending or descending and position of key. As it stands
C----it is a bit more comlicated to derive. To minimise 
C----change in the calling subroutine use previous binsort style
C----keybuf. 
      DO   IK=1,NKEYS
         J = (IK-1)*5+1
         SAVE_KEYS(IK) = KEYBUF(J+2)/4 + 1
         ASCEND_DESCEND(IK) = 1
         IF(KEYBUF(J+1).NE.0) ASCEND_DESCEND(IK) = -1
      ENDDO
      NKEYS_S = NKEYS
      IRECORD_POSITIONS = 0
      NRECORD_S = 0
      DO    IREC=1,NREC_L
         found = .FALSE.
         DO   IK=1,NKEYS
            IF(SAVE_KEYS(IK).EQ.IREC) then 
               found = .TRUE.
               exit
            endif
         ENDDO
         if(.not.found) then
            IRECORD_POSITIONS = IRECORD_POSITIONS + 1
            SAVE_RECORD(IRECORD_POSITIONS) = IREC
            NRECORD_S = NRECORD_S + 1
         endif
      ENDDO
c
C---Now we are ready to accept records. It will be done by another 
C---routine.
C---If external merging is necessary it should be initialised here.
C
      NRECORD_NOW = 0
      NRECORD_RETURN = 0
      SRTBEG = 0
      RETURN
      END
C
      INTEGER FUNCTION SRTRLS(ADATA)
      use sorting_commons
      IMPLICIT NONE

C
C----Recieve one record and save it in the list.
C----For large number of records if number of current records
C----is equal to maximum number of records then they should be 
C----sorted and then written to external files taken care of
C----distribution. (polyphase merging is possible option to use)
      REAL ADATA(*)

C
      INTEGER IKEY_NOW,IREC_NOW,IK,IR
C
C---First save keys.
      NRECORD_NOW = NRECORD_NOW + 1
      IF(NRECORD_NOW.GT.NMAX_REC) THEN
                                ! Memory is not sufficient. Reallocate
         call sorting_reallocate_this()
      ELSEIF(NRECORD_NOW.EQ.NMAX_REC) THEN
C
C---Memory is not enough for internal sorting. External sorting
C---part should be written. In that case available records 
C---should be sorted and written to file with distribution
c---which is going to be used. (polyphase distribution is one possibility)
C
Cmdw  binsort uses 4 scratch files (T=4)
Cmdw  Reads from stdin (see fillworka).
Cmdw  Each time buffer is full, do sort (see top of readrun; later
Cmdw  fwrite is only if all data fits into first buffer and disk not needed).
Cmdw  Write sorted buffer to scratch file (see writerun).
Cmdw  Read more, and write to next scratch file - see D2-D4 in merge().
Cmdw  When finished, rewind scratch files and go to merge step D5.

Cmdw  Open scratch files here (with QOPEN) iff necessary.
Cmdw  When finished, need to rewind with QSEEK.

      ENDIF

      IKEY_NOW = (NRECORD_NOW-1)*NKEYS_S
      DO   IK=1,NKEYS_S
         keys_mem(ikey_now+ik) = adata(save_keys(ik))*ascend_descend(ik)
      ENDDO
      IREC_NOW = (NRECORD_NOW-1)*NRECORD_S

      DO  IR=1,NRECORD_S
         array_mem(irec_now + ir) = adata(save_record(ir))
      ENDDO
      index_mem(nrecord_now) = nrecord_now
C
C---Normal return. No disaster.
      SRTRLS = 0
      RETURN
      END
C
      INTEGER FUNCTION SRTMRG()
      use sorting_commons
C
      IMPLICIT NONE
C
C---This function should do merging. But here we use only sorting
C---It will have to be expanded for merging for large number of records
C
      CALL HEAP_SORT(NRECORD_NOW,NKEYS_S,keys_mem,index_mem)
C
C---Records have been sorted. They should be distributed. 
C---But it is next stage. 
      SRTMRG = 0
      RETURN
      END
C
      INTEGER FUNCTION SRTRET(ADATA)
      use sorting_commons
      IMPLICIT NONE
C
C----Retrieve next record from the sorted list of the records.
      REAL ADATA(*)
C
C---This function should do merging. But here we use only sorting
C---It will have to be expanded for merging for large number of records
C
      INTEGER IK,IR,IKEY_REC,IREC_R
      integer kk
C
C----Take keys first.
      NRECORD_RETURN = NRECORD_RETURN + 1
      IF(NRECORD_RETURN.GT.NRECORD_NOW) THEN
C
C----In case of external search read from file new set of records
C----and then check if everything is o.k
C
C---Distaster. Calling subroutine wants more than it has given.
         SRTRET = -1
         RETURN
      ENDIF
C
C---Take keys.
      IKEY_REC = (NRECORD_RETURN-1)*NKEYS_S
      DO   IK=1,NKEYS_S
         kk = save_keys(ik)
        ADATA(kk) = keys_mem(ikey_rec + kk)*ascend_descend(ik)
      ENDDO
C
C--Now take records.
      IREC_R = (nint(index_mem(NRECORD_RETURN))-1)*nrecord_s
      DO   IR=1,NRECORD_S
         ADATA(SAVE_RECORD(IR)) = array_mem(IREC_R +  IR)
      ENDDO
C
C---Succesful retrieval
      SRTRET = 0
      RETURN
      END
C
C---Internal sorting part. It could be improved
C-----------
      SUBROUTINE  HEAP_SORT(N,NKEYS,A_KEY,INDEX_R)
      IMPLICIT NONE
C
C----Sorting using heapsort. 
C----A        contains keys
C----INDEX_R  index of records to be sorted
C----N        number of records
C----NKEYS    number of keys
C
C---Reference
C
c---Knuth, The art of computer programming Volume 3
C---1998
C
      INTEGER    N,NKEYS
      REAL       INDEX_R(*)
      REAL       A_KEY(*)
C
      REAL    T_KEY(5),T1_KEY(5)
C
      REAL      INDEX_C
      INTEGER   L,M,I_TEMP,L1,MKEY1_ADD,MKEY2_ADD,I_KEY,N1
C
C--------------------------------------------------------------
      L1 = N/2
C
C---Create heap. everybody tends to reach his level of incomptence
c      WRITE(*,*)'In heap_sort'

      DO   L=L1,1,-1
        CALL SIFT_UP(N,L,NKEYS,A_KEY,INDEX_R)
      ENDDO
C
c--Remove heaps  one after another
      N1 = N-1
      DO   M=N1,1,-1
        INDEX_C = INDEX_R(1)
        INDEX_R(1) = INDEX_R(M+1)
        INDEX_R(M+1) = INDEX_C
        MKEY1_ADD = M*NKEYS
        DO   I_KEY=1,NKEYS
          MKEY2_ADD = MKEY1_ADD + I_KEY
          T_KEY(I_KEY) = A_KEY(I_KEY)
          A_KEY(I_KEY) = A_KEY(MKEY2_ADD)
          A_KEY(MKEY2_ADD) = T_KEY(I_KEY)
        ENDDO

        CALL SIFT_UP(M,1,NKEYS,A_KEY,INDEX_R)
      ENDDO
      END
C
      SUBROUTINE SIFT_UP(M,L,NKEYS,A_KEY,INDEX_R)
c
c---Sift ip process in heap sort. This implementation is 
C---intended to work for multykey cases. I am not sure this treatment is best
C---for multy key cases. Other idea will have to be tested.
      IMPLICIT NONE
      INTEGER M,L,NKEYS
      REAL INDEX_R(*)
      REAL A_KEY(*)
C
      INTEGER I,J,IKEY_ADD,NKEYS1,JKEY_ADD,JKEY1_ADD,LKEY_ADD,
     &        I_KEY,J_KEY,I1_KEY
      REAL KEY_T(5),INDEX_C
C
C---Intitialise sift up (in Knuth's terminology H2)
      INDEX_C = INDEX_R(L)
      NKEYS1  = -NKEYS
      LKEY_ADD = L*NKEYS+NKEYS1
      DO  I_KEY=1,NKEYS
        KEY_T(I_KEY) = A_KEY(LKEY_ADD + I_KEY)
      ENDDO
c
C---H3
      J = L
 40   CONTINUE
C
c---Start sift up. Compare I with 2*I or 2*I+1
C---H4
      I = J
      J = J + J
C
C---If(J.GT.M) then return time
      IKEY_ADD = I*NKEYS + NKEYS1
      JKEY_ADD = J*NKEYS + NKEYS1
      IF(J.LE.M) THEN
         IF(J.LT.M) THEN
C
c---J.lt.M find maximum of J and J + 1 (which is 2*I, 2*I+1)
           JKEY1_ADD = JKEY_ADD + NKEYS
C
C---H5 
           DO   I_KEY=1,NKEYS
             IF(A_KEY(JKEY_ADD+I_KEY).LT.A_KEY(JKEY1_ADD+I_KEY)) THEN
C
c---Next record is smaller, Take it
               J = J + 1
               JKEY_ADD = JKEY1_ADD
               GOTO 55
             ELSE IF(A_KEY(JKEY_ADD+I_KEY).GT.
     &                                A_KEY(JKEY1_ADD+I_KEY))THEN
C
c--Next record is greater. No need to proceed
               GOTO 55
             ENDIF
C
c----Cannot decide yet. Check another key
           ENDDO
         ENDIF
 55      CONTINUE
         DO   I_KEY = 1,NKEYS
C
C---compare saved record in initialisation with record J. If
C---saved is smaller than current then then put for i (2*I or 2*I+1)
           IF(KEY_T(I_KEY).LT.A_KEY(JKEY_ADD + I_KEY)) THEN
C
C---H6 above H7 below
             INDEX_R(I) = INDEX_R(J)
             DO   I1_KEY=1,NKEYS
               A_KEY(IKEY_ADD+I1_KEY) = A_KEY(JKEY_ADD+I1_KEY)
             ENDDO
             GOTO 40
           ELSE IF(KEY_T(I_KEY).GT.A_KEY(JKEY_ADD + I_KEY)) THEN
             GOTO 80
           ENDIF
         ENDDO
       ENDIF
 80    CONTINUE
C
c---Put saved record to i-th and exit from sift up
C---H8
       INDEX_R(I) = INDEX_C
       DO   I_KEY=1,NKEYS
         A_KEY(IKEY_ADD+I_KEY) = KEY_T(I_KEY)
       ENDDO
       RETURN
       END

