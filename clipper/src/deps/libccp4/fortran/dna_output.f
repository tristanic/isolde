C
C     dna_output.f: write output for the DNA project
C     Copyright (C) 2004  CCLRC, Graeme Winter
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
c     dna_output.f
c     maintained by G.Winter
c     "standard" DNA output XML subroutines for 
c     fortran programs - these will create named tables, lists and items
c     
c     tables containeth lists, which in turn containeth items.
c     
c     These may be applied to programs other than Mosflm - so long
c     as I can code up a clever way of handling this output!
c     
c     
c     
c     
c     
c     $Id$
c     

c     set no output flag - this is the alternative to starting the output
c     and should be called for good form - an alternative would be to
c     call dna_start with filename = ' ' - which I could interpret
c     appropriately

      subroutine dna_set_no_output
      implicit none
      include 'dna_header.fh'
      dnaout = .false.
      return
      end

c     
c     start writing the dna output - this takes a program name, 
c     for instance 'mosflm' - this should be called as soon as
c     the file is opened
c     

      subroutine dna_start(filename, progname)
      implicit none
      include 'dna_header.fh'
      character*(*) progname, filename
      integer ifail, lenstr
      external lenstr
      dnainlist = .false.
      dnaintable = .false.
      dna_image = ' '
      ifail = 1
      call ccpdpn(dnafd, filename, 'UNKNOWN', 'F' ,0 , ifail)
      if (ifail.ne.1) then
        call ccperr(2, 'dna_start: error opening file')
        dnaout = .false.
      else
        dnaout = .true.
      endif
 1    format('<?xml version="1.0"?><!DOCTYPE dna_tables>')
 2    format('<dna_tables program="', a, '">')
      if(dnaout) write(dnafd, 1)
      if(dnaout) write(dnafd, 2) progname(1:lenstr(progname))
      return
      end

c     
c     finish writing dna output - this should go at the end of the 
c     program
c     

      subroutine dna_end
      implicit none
      include 'dna_header.fh'
 1    format('</dna_tables>')

c     close any ongoing tables
      if (dnainlist) call dna_list_end
      if (dnaintable) call dna_table_end

      if(dnaout) write(dnafd, 1)
      if(dnaout) close(dnafd)
      dnaout = .false.
      return
      end

c     
c     write an item containing a "real"
c     

      subroutine dna_real_item(name, value)      
      implicit none
      include 'dna_header.fh'
      character *(*) name
      real value
      integer lenstr
      external lenstr
 1    format('      <item name="', a, '">', e15.6, '</item>')
      if(dnaout) write(dnafd, 1) name(1:lenstr(name)), value
      return
      end

      subroutine dna_double_item(name, value)      
      implicit none
      include 'dna_header.fh'
      character *(*) name
      double precision value
      integer lenstr
      external lenstr
 1    format('      <item name="', a, '">', e15.6, '</item>')
      if(dnaout) write(dnafd, 1) name(1:lenstr(name)), value
      return
      end

c     
c     as above sed 's/real/integer/'
c     


      subroutine dna_integer_item(name, value)
      implicit none
      include 'dna_header.fh'
      character *(*) name
      integer value, lenstr
      external lenstr
 1    format('      <item name="', a, '">', i15, '</item>')
      if(dnaout) write(dnafd, 1) name(1:lenstr(name)), value
      return
      end

c     
c     as above sed 's/real/character*(*)/'
c     

      subroutine dna_character_item(name, value)
      implicit none
      include 'dna_header.fh'
      character *(*) name
      character *(*) value
      integer lenstr
      external lenstr
 1    format('      <item name="', a, '">', a, '</item>')
      if(dnaout) write(dnafd, 1) name(1:lenstr(name)), value
      return
      end

c     
c     start a named list
c     

      subroutine dna_list_start(name)
      implicit none
      include 'dna_header.fh'
      character *(*) name
 1    format('    <list name="', a, '">')

c     check that we are not already in a list, and if we are
c     close it!

      if (dnainlist) call dna_list_end
      dnainlist = .true.

      if(dnaout) write(dnafd, 1) name
      return
      end

      subroutine dna_ilist_start(name, index)
c     this is the same as the above but with an integer index
c     so that you can have any lists with the same name - very 
c     important for tabular output.
      implicit none
      include 'dna_header.fh'
      character *(*) name
      integer index
 1    format('    <list name="', a, '" index="', i5, '">')

c     check that we are not already in a list, and if we are
c     close it!

      if (dnainlist) call dna_list_end
      dnainlist = .true.

      if(dnaout) write(dnafd, 1) name, index
      return
      end

c     
c     finish a named list
c     

      subroutine dna_list_end
      implicit none
      include 'dna_header.fh'

 1    format('    </list>')

c     check that we are in a list
      if (dnainlist) then
         if(dnaout) write(dnafd, 1)
      end if
      dnainlist = .false.
      return
      end

c     
c     start a named table
c     

      subroutine dna_table_start(name)
      implicit none
      include 'dna_header.fh'
      character *(*) name
      integer lenstr
      external lenstr
 1    format('  <table name="', a, '" image="', a, '">')
 2    format('  <table name="', a, '">')

c     check that we are not already inside a table
      if (dnaintable) call dna_table_end
      dnaintable = .true.

      if (dna_image .eq. ' ') then
         if(dnaout) write(dnafd, 2) name(1:lenstr(name))
      else
         if(dnaout) write(dnafd, 1) name(1:lenstr(name)), 
     +        dna_image(1:lenstr(dna_image))
      end if         
      return
      end

c     
c     finish a named table
c     

      subroutine dna_table_end
      implicit none
      include 'dna_header.fh'
 1    format('  </table>')

c     check we are not still inside a list

      if (dnainlist) then
         call dna_list_end
      end if

      if (dnaintable) then
         if(dnaout) write(dnafd, 1)
      end if
      dnaintable = .false.
      return
      end

      subroutine dna_error(message)
      implicit none
      character*(*) message
      integer lenstr
      external lenstr

      call dna_table_start('error')
      call dna_list_start('error')
      call dna_character_item('message', message(1:lenstr(message)))
      call dna_list_end
      call dna_table_end

      return
      end

      subroutine dna_warning(message)
      implicit none
      character*(*) message
      integer lenstr
      external lenstr

      call dna_table_start('warning')
      call dna_list_start('warning')
      call dna_character_item('message', message(1:lenstr(message)))
      call dna_list_end
      call dna_table_end

      return
      end

