C
C     keyparse.f: high-level interface to CCP4 parser functions
C     Copyright (C) 1995  Dave Love, Kevin Cowtan
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
C     A simplified interface to the %$!$#@$ (so-called) parser:
C
C     The idea is to call a routine (MEMOPARSE) to read a keyworded line
C     into a hidden internal buffer and then have a simple chain of
C     calls to a set of routines which check the appropriate keyword,
C     look for any extra arguments associated with it, and set values in
C     their arguments as appropriate.  No if ... then ... else (or, much
C     worse, assigned goto) is necessary, nor maintaining the parser
C     arrays -- the relevant checking is done internally.  At the end of
C     the checks, call PARSEDIAGNOSE to print any appropriate messages
C     and loop back to MEMOPARSE if PARSEDIAGNOSE's argument is true
C     (otherwise continue and process the data read in).  You don't need
C     to check for `END' or end-of-file.  Geddit?
C
C     Escape hatch: use PARSEKEYARG to get all the tokens after a random
C     keyword and call PARSE (or whatever) to deal with them as
C     necessary.  This is usually deprecated, however -- go for
C     simply-structured input.
C
C     10 call memoparse (.true.)  ! setup and echo i/p
C        call parseint ('IVAL', ival)  ! ival.eq.3 after `IVAL 3'
C        call parsecell (cell)   ! assign cell from `CELL 20 30 40' etc.
C        call parsekeyarg ('FRED', rest)  ! the full horror in REST
C        call parse (rest, ....)
C         [firkle with the `parse'd arrays]
C         ...
C        call parsediagnose (cont) ! check
C        if (cont) goto 10
C   c     now do something useful...
C
C     Fixme: the set of routines below might need extending...
C     Fixme: consider whether more obscure names worthwhile to avoid
C     possible clashes.
C
C     Dave Love $Date$
C     Additions by Kevin Cowtan
C     
      subroutine memoparse (echo)
C
C     Call PARSER and stash the returned values away for later testing
C     when the other entrypoints are called. (OK, so it's not really
C     memoisation...).
C
C     ECHO is set to echo the parser i/p.

      implicit none
C     Args
      character*(*) key, subkey, spgnam, pgname, rest
      character*30 prglab(*)
      character*1 chnam
      logical echo, flag, cont
      integer ival, nsym, numsgp, nsymp, n, ivals (n), nth, mtznum,
     +  nprglab, toks, inat0, inat1, ires0, ires1, imode, ifail
      real rval, cell (6), rsym (4,4,*), resmin, resmax, smin, smax,
     +  rvals(n)

C     Stores for parser stuff
      integer maxtoks, maxline
C     Let's not mess around...
      parameter (maxtoks = 500, maxline=2000)
      integer ibeg(maxtoks), iend(maxtoks), ityp(maxtoks), idec(maxtoks)
      real fvalue (maxtoks)
      character*4 cvalue (maxtoks), memokey, memosubkey
      character line*(maxline)
      integer ntok

C     locals 
      logical someerr, eof, argerr, success(maxtoks)
      integer i, k

      save
      data someerr, eof /2*.false./
      
      ntok=maxtoks
      argerr = .false.
      line = ' '
C     in case of immediate EOF:
      success(1) = .true.
      call parser(memokey, line, ibeg, iend, ityp, fvalue, cvalue, idec,
     +     ntok, eof, echo)
C     END == EOF always
      if (memokey.eq.'END') eof = .true.
C     not sure if necessary:
      if (eof) memokey = ' '
      do i=1,ntok
       success(i) = .false.
      enddo
      return
      
      entry parsekey (key, flag)
C     bare KEY -- set FLAG if found
      if (memokey.eq.key) then
C       matched key
        if (ntok.eq.1) then
          success(1) = .true.
          flag = .true.
        else
          argerr = .true.
          call lerror (1, 0, 'No argument expected')
        end if
      end if
      return

      entry parsekeyarg (key, rest)
C     KEY + rest of line -- rest of line returned in REST
      if (memokey.eq.key) then
C       matched key
        if (ntok.gt.1) then
          do i=1,ntok
           success(i) = .true.
          enddo
          rest = line (ibeg(2):iend(ntok))
        else
          rest = ' '
          success(1) = .true.
C          argerr = .true.
C          call lerror (1, 0, 'Argument expected')
        end if
      end if
      return

      entry parseint (key, ival)
C     KEY + integer -- returned in IVAL
      if (memokey.eq.key) then
C       matched key
        if (ntok.eq.2 .and. ityp (2).eq.2) then
          ival = nint (fvalue(2))
          do i=1,2
           success(i) = .true.
          enddo
        else 
          argerr = .true.
          call lerror (1, 0, 'Integer argument expected')
        end if
      end if
      return
      
      entry parsereal (key, rval)
C     KEY + real -- returned in RVAL
      if (memokey.eq.key) then
C       matched key
        if (ntok.eq.2 .and. ityp (2).eq.2) then
          rval = fvalue(2)
          do i=1,2
           success(i) = .true.
          enddo
        else 
          argerr = .true.
          call lerror (1, 0, 'Real argument expected')
        end if
      end if
      return
      
      entry parsenargs(key, toks)
      toks = 0
      if (memokey .eq. key) toks = ntok
      return

      entry parsenints (key, n, ivals)
C     KEY + upto N integers -- N reset to number found, returned in IVALS
      if (memokey.eq.key) then
       success(1)=.true.
        if (ntok.ge.2 .and. ntok.le. n+1) then
          do i = 1, min (n, ntok-1)
            if (ityp (i+1).ne.2) return
            ivals (i) = nint (fvalue (i+1))
            n = i
            success(i+1) = .true.
          end do
        else
          argerr = .true.
          n = ntok
          call lerror (1, 0, 'Incorrect number of integer arguments')
        end if
      end if
      return
      
      entry parsenreals (key, n, rvals)
C     KEY + upto N reals -- N reset to number found, returned in RVALS
      if (memokey.eq.key) then
       success(1)=.true.
        if (ntok.ge.2 .and. ntok.le. n+1) then
          do i = 1, min (n, ntok-1)
            if (ityp (i+1).ne.2) return
            rvals (i) = fvalue (i+1)
            n = i
            success(i+1) = .true.
          end do
        else
          argerr = .true.
          n = ntok
          call lerror (1, 0, 'Incorrect number of real arguments')
        end if
      end if
      return
      
      entry parsesubkey (key, subkey, flag)
C     KEY + subkeyword SUB -- set FLAG if found
      if (memokey.eq.key) then
C       matched key
       success(1) = .true.
       if (subkey.ne.' ') then
        do k=2,ntok
         memosubkey=cvalue(k)
         call ccpupc(memosubkey)
         if (memosubkey.eq.subkey) then
          flag=.true.
          success(k)=.true.
         endif
        enddo
       else
        flag=.true.
       end if
      end if
      return
      
      entry parsesubkeyarg (key, subkey, nth, rest)
C     KEY + subkey + rest of line -- rest of line returned in REST
      if (memokey.eq.key) then
C       matched key
        success(1) = .true.
        k = 1
        if (subkey .ne. ' ') then
          k = 9999
          do i=2,ntok
            memosubkey=cvalue(i)
            call ccpupc(memosubkey)
            if (memosubkey.eq.subkey) then
             k = i
             goto 90
            endif
          enddo
        endif
   90   if (k.le.ntok) then
          if (ntok.ge.k+nth) then
            rest = line (ibeg(k+nth):iend(ntok))
            do i=2,ntok
              success(i) = .true.
            enddo
          else
            argerr = .true.
            call lerror (1, 0, 'Argument expected after sub-keyword')
          end if
        end if
      end if
      return

      entry parsesubint (key, subkey, nth, flag, ival)
C     KEY + n'th integer after subkey -- returned in IVAL
C     ERROR only if flag=true
      if (memokey.eq.key) then
C  ... matched key
       success(1) = .true.
       k=1
       if (subkey.ne.' ') then
        k=9999
        do i=2,ntok
         memosubkey=cvalue(i)
         call ccpupc(memosubkey)
         if (memosubkey.eq.subkey) then
          k=i
          goto 100
         endif
        enddo
       endif
  100  if (k.le.ntok) then
C  .... matched subkey (if set)
        success(k) = .true.
        if (ntok.ge.nth+k .and. ityp(nth+k).eq.2) then
         ival = nint (fvalue(nth+k))
         success(nth+k) = .true.
        else if (flag) then
         argerr = .true.
         call lerror (1, 0, 'Integer sub-argument expected')
        endif
       endif
      endif
      return
      
      entry parsesubreal (key, subkey, nth, flag, rval)
C     KEY + n'th real after subkey -- returned in RVAL
C     ERROR only if flag=true
      if (memokey.eq.key) then
C  ... matched key
       success(1) = .true.
       k=1
       if (subkey.ne.' ') then
        k=9999
        do i=2,ntok
         memosubkey=cvalue(i)
         call ccpupc(memosubkey)
         if (memosubkey.eq.subkey) then
          k=i
          goto 110
         endif
        end do
       endif
  110  if (k.le.ntok) then
C  .... matched subkey (if set)
        success(k) = .true.
        if (ntok.ge.nth+k .and. ityp(nth+k).eq.2) then
         rval = fvalue(nth+k)
         success(nth+k) = .true.
        else if (flag) then
         argerr = .true.
         call lerror (1, 0, 'Real sub-argument expected')
        endif
       endif
      endif
      return
      
      entry parsesubchar (key, subkey, nth, flag, rest)
C     KEY + n'th string after subkey -- returned in REST
      if (memokey.eq.key) then
C  ... matched key
       success(1) = .true.
       k=1
       if (subkey.ne.' ') then
        k=9999
        do i=2,ntok
         memosubkey=cvalue(i)
         call ccpupc(memosubkey)
         if (memosubkey.eq.subkey) then
          k=i
          goto 120
        endif
        end do
       endif
  120  if (k.le.ntok) then
C  .... matched subkey (if set)
        success(k) = .true.
        if (ntok.ge.nth+k .and. (ityp(nth+k).ne.0 .or. .not.flag) ) then
         rest = line(ibeg(nth+k):iend(nth+k))
         success(nth+k) = .true.
        else if (flag) then
         argerr = .true.
         call lerror (1, 0, 'Character sub-argument expected')
        endif
       endif
      endif
      return
      
      entry parsecell (cell)
C     CELL -- returned in CELL
      if (memokey.eq.'CELL') then
        call rdcell (2, ityp, fvalue, ntok, cell)
        do i=1,ntok
         success(i) = .true.
        enddo
      end if
      return
      
      entry parsesymm (spgnam, numsgp, pgname, nsym, nsymp, rsym)
C     SYMMetry -- usual values returned
      if (memokey.eq.'SYMM') then
        nsym = 0
        call rdsymm(2, line, ibeg, iend, ityp, fvalue, ntok, spgnam,
     +       numsgp, pgname, nsym, nsymp, rsym)
        do i=1,ntok
         success(i) = .true.
        enddo
      end if
      return

      entry parseatomselect(key, inat0, inat1, ires0, ires1, chnam,
     +                      imode)
C     KEYword followed by atom/residue selection syntax
C     key atom <ai> [[to] <aj>] |
C         residue [all|ons|ca] [chain <chn>] <ri> [[to] <rj>]
C     Returns values of <ai>,<aj>,imode...
      if (memokey.eq.key) then
C       matched key
        ifail = -1
        call rdatomselect(2, inat0, inat1, ires0, ires1, chnam, imode,
     +                     ntok, line, ibeg, iend, ityp, idec, fvalue,
     +                     ifail)
        do i=1,ntok
          if (i.ne.ifail) success(i)=.true.
        enddo
      endif
      return

      entry parsereso (resmin, resmax, smin, smax)
C     RESOlution -- usual values returned
      if (memokey.eq.'RESO') then
        call rdreso (2, ityp, fvalue, ntok, resmin, resmax, smin, smax)
        do i=1,ntok
         success(i) = .true.
        enddo
      end if
      return
      
      entry parselabin(mtznum,prglab,nprglab)
C     LABIn -- requires mtz file number, program labels, and number of.
      if (memokey.eq.'LABI') then
        call lkyin(mtznum,prglab,nprglab,ntok,line,ibeg,iend)
        do i=1,ntok
         success(i) = .true.
        enddo
      end if
      return
      
      entry parselabout(mtznum,prglab,nprglab)
C     LABOut -- requires mtz file number, program labels, and number of.
      if (memokey.eq.'LABO') then
        call lkyout(mtznum,prglab,nprglab,ntok,line,ibeg,iend)
        do i=1,ntok
         success(i) = .true.
        enddo
      end if
      return

      entry parsefail()
      do i=1,ntok
       success(i) = .false.
      enddo
      return

      entry parseend(key)
C     KEY -- is treated as input terminator as same as END
      if (memokey .eq. key) eof = .true.
      return
      
      entry parsediagnose (cont)
C     Call at end of tests for possible 'Invalid keyword' diagnostic or
C     abort if at EOF and had an error.  Continue processing (no EOF)
C     if (cont).
      if (.not.argerr .and. .not.eof) then
       if (.not.success(1)) then
        call lerror (1, 0, 'Invalid keyword')
        someerr = .true.
       endif
       do i=2,ntok
        if (.not.success(i)) then
         write (line,950)i
 950     format ('Invalid sub-keyword in position ',i2)
         call lerror (1, 0, line)
         someerr = .true.
        endif
       enddo
      endif
      if (argerr) then
        someerr = .true.
      end if
      argerr = .false.
      if (eof) then
        cont = .false.
        if(someerr) call ccperr (1, 'Input error (see above)')
      else
        cont = .true.
      end if
      do i=1,ntok
       success(i) = .false.
      enddo
      return
C
      end
