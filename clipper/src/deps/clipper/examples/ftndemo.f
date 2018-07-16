      program ftndemo
C crystal definitions
      integer nspgr
      real cell(6)
      real resolution
C atom data definitions
      integer natoms
      real x(10),y(10),z(10),occ(10),b(10)
      integer atno(10)
C reflection data definitions
      integer nrefls
      integer h(100000),k(100000),l(100000)
      real f(100000),phi(100000)
      integer h0, k0, l0

C crystal data
      data nspgr/19/
      data cell/30.0,40.0,50.0,90.0,90.0,90.0/

C atom data
      data natoms/10/
      data x/12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0/
      data y/1.0,2.0,3.0,1.0,2.0,3.0,1.0,2.0,3.0,1.0/
      data z/6.0,8.0,7.0,9.0,6.0,8.0,7.0,9.0,6.0,8.0/
      data occ/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/
      data b/20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0,20.0/
      data atno/6,7,6,8,6,7,6,8,6,7/

c make a crude list of reflections
      nrefls = 0
      do h0 = 0,10
         do k0 = 0,10
            do l0 = 0,10
               if ( h0*h0+k0*k0+l0*l0 .lt. 100 ) then
                  nrefls = nrefls + 1
                  h(nrefls) = h0
                  k(nrefls) = k0
                  l(nrefls) = l0
               endif
            enddo
         enddo
      enddo

c do structure factor calc
      call sfcalc( nspgr, cell,
     +             natoms, x, y, z, occ, b, atno,
     +             nrefls, h, k, l, f, phi )

c output results
      do i = 1, nrefls
         write (*,*)h(i),k(i),l(i),f(i),phi(i)
      enddo

      stop
      end
