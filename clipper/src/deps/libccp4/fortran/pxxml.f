C
C     pxxml.f: write XML tags into program outout
C     Copyright (C) 2001  Alun Ashton
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
C
C
C    =========
C    pxxml.f
C    =========
C
C     Author  : a.w.ashton@ccp4.ac.uk
C     Version : 1.0  
C     future  : translate to C or C++
C  
C    Assumption 1 : Exsisting/New file
C     XML files are only to be writen from new. The design must take into consideration 
C     the possability of future extension so that XML file might be reopned. 
C
C    Assumption 2 : Sizes
C     The maximum size for an Element value or name is 80 characters.
C     The maximum size for an Attribute value or name is 80 characters. 
C     The maximum size for program name is 20 characters.
C     The maximum size for XML filename – including path is 70 characters.
C     The maximum number of elements is 80.
C
C    Assumption 3: Language
C     In the first instance the libs will be writen in fortran to test the other 
C     assumptions!

C     SUBROUTINE XMLOpen (XMLFileUnit, rootvalue, ifail)
C
C     First subroutine to be called.
C     Arguments: 	XMLFIleUnit : Input - No value
C                           Output returns the file Unit number
C               root value = eg html
C		Ifail : same as usual!
C
C     Purpose: open an XML file.
C
C     See if XMLFILE has been assigned e.g. on command line.
C     If XMLFILE not assigned use program name.xml.
C     If file exists overwrite it.
C     Write out the header.
C     No provision for existing XML files to be reopened

      subroutine XMLOpen (XMLFileUnit, rootvalue, ifail)

      logical ccpexs
      
      character*(*) rootvalue

      integer XMLFileUnit, ifail, elementstatus(80), elementnumber 
      character*80 openelements(80), XMLFileName
      character*1024 XMLOutputLine

      character version*10, date*8, progname*32
      
      common /ccp4xmld/ elementnumber, elementstatus, openelements, 
     $     XMLFileName,XMLOutputLine

      integer lenstr
      external ccpexs,lenstr

      XMLFileName=' '
      elementnumber=0
      elementstatus(1)=0

      call ccp4_version(version)
      call ccpdat(date)
      call ccppnm(progname)

      call ugtenv('XMLFILE',XMLFileName)

      if (XMLFileUnit.ne.6) then
        if (XMLFileName.eq.' ') 
     +       XMLFileName=progname(1:lenstr(progname))//'.xml'
        
        if (ccpexs(XMLFileName)) call ccperr(2,'XMLFILE OVERWRITTEN')
        
        call ccpdpn(XMLFileUnit,XMLFileName,'UNKNOWN','F',0,ifail)
      endif

C     check to see if the rootvalue has a value or default to html
      if (lenstr(rootvalue).le.1) rootvalue='html'

      if (rootvalue(1:lenstr(rootvalue)).ne.'html') 
     $     write (XMLFileUnit,'(A)') '<?xml version="1.0"?>'

      call XMLOpenElement(XMLFileUnit, rootvalue(1:lenstr(rootvalue)), 
     $     ifail) 

      if (rootvalue(1:lenstr(rootvalue)).ne.'html') then
C     First element should be information about the program
        call XMLOpenElement(XMLFileUnit, 
     $       progname(1:lenstr(progname)), ifail)
        call XMLWriteAttribute(XMLFileUnit, 
     $       progname(1:lenstr(progname)),
     $       'ccp4_version', version(1:lenstr(version)), ifail)
        call XMLWriteAttribute(XMLFileUnit, 
     $       progname(1:lenstr(progname)),
     $   'date', date(1:lenstr(date)), ifail)
        call XMLCloseElement(XMLFileUnit, progname(1:lenstr(progname)), 
     $       ifail)
      endif

      return
      end

C      subroutine XMLWriteJournalCitation(XMLFileUnit, AuthorsString, 
C     +  JournalString, VolumeString, PageStart, PageEnd, YearString, 
C     +  ifail)
C      For example:
C        call XMLWriteJournalCitation(XMLFILEUNIT,'B.W.Matthews',
C     +   'J.Mol.Biol.','33','491','497','1968',ifail)
C      If using the Bioxhit schema, this should only be called once 
C      for each processing_stage.

      subroutine XMLWriteJournalCitation(XMLFileUnit, AuthorsString, 
     +  JournalString, VolumeString, IssueString, PageStart, PageEnd, 
     +  YearString, ifail)

      integer XMLFileUnit, ifail
      character*(*) AuthorsString, JournalString, VolumeString, 
     +  IssueString, PageStart, PageEnd, YearString

      call XMLOpenElement(XMLFileUnit,'citation',ifail)

      call XMLWriteElement(XMLFileUnit, 
     +    'citation_type', 'journal', ifail)
      if (AuthorsString.ne.' ') call XMLWriteElement(XMLFileUnit, 
     +    'citation_id', AuthorsString, ifail)
      if (JournalString.ne.' ') call XMLWriteElement(XMLFileUnit, 
     +    'journal_abbreviation', JournalString, ifail)
      if (VolumeString.ne.' ') call XMLWriteElement(XMLFileUnit, 
     +     'journal_volume', VolumeString, ifail)
      if (IssueString.ne.' ') call XMLWriteElement(XMLFileUnit, 
     +     'journal_issue', IssueString, ifail)
      if (PageStart.ne.' ') call XMLWriteElement(XMLFileUnit, 
     +     'page_from', PageStart, ifail)
      if (PageEnd.ne.' ') call XMLWriteElement(XMLFileUnit, 
     +     'page_to', PageEnd, ifail)
      if (YearString.ne.' ') call XMLWriteElement(XMLFileUnit, 
     +     'publication_year', YearString, ifail)

      call XMLCloseElement(XMLFileUnit,'citation',ifail)

      return
      end

C      subroutine XMLWriteMessages(XMLFileUnit, ccperr_level, 
C     +  MessageString, ifail)
C      For example: 
C        call XMLWriteMessage(XMLFileUnit,0,'Normal Termination',ifail)
C      If using the Bioxhit schema, this should only be called once 
C      for each processing_stage.

      subroutine XMLWriteMessages(XMLFileUnit, ccperr_level, 
     +  MessageString, ifail)

      integer XMLFileUnit, ifail, ccperr_level
      character*(*) MessageString
      character MessageCode*4,MessageSeverity*12

      write(MessageCode,'(I4)') ccperr_level
      if (ccperr_level.eq.0) MessageSeverity='information'
      if (ccperr_level.eq.1) MessageSeverity='fatal error'
      if (ccperr_level.eq.2) MessageSeverity='warning'
      if (ccperr_level.eq.3) MessageSeverity='information'
      if (ccperr_level.eq.4) MessageSeverity='information'

      call XMLOpenElement(XMLFileUnit,'messages',ifail)
      call XMLWriteElement(XMLFileUnit, 
     +    'message_code', MessageCode, ifail)
      if (MessageString.ne.' ') call XMLWriteElement(XMLFileUnit, 
     +    'message_text', MessageString, ifail)
      call XMLWriteElement(XMLFileUnit, 
     +    'message_severity', MessageSeverity, ifail)        
      call XMLCloseElement(XMLFileUnit,'messages',ifail)

      return
      end

C subroutine XMLOpenElement(XMLFileUnit, ElementName, ifail)
C     XMLFileUnit - integer - if file not already opened it will be
C     element name = value is char string.
C        completes a tag and reates a sub element but leaves tag open
C        e.g. <oldtag
C             >
C              <leftopen

      subroutine XMLOpenElement(XMLFileUnit, ElementName, ifail)

      logical ccpexs

      character*(*) ElementName
      character*32 progname

      character*80 indentline
      integer indent

      integer XMLFileUnit, ifail, elementstatus(80), elementnumber 
      character*80 openelements(80), XMLFileName
      character*1024 XMLOutputLine
      
      common /ccp4xmld/ elementnumber, elementstatus, openelements,
     $     XMLFileName,XMLOutputLine
c      save /ccp4xml/

      integer lenstr
      external ccpexs,lenstr

      indentline=' ' 

      if (.not.ccpexs(XMLFileName) .and. .not.(XMLFileUnit.eq.6)) then
        call ccppnm(progname)
        call ccplwc(progname)
        call XMLOpen (XMLFileUnit,progname(1:lenstr(progname)), ifail)
      endif

      elementnumber=elementnumber+1

      openelements(elementnumber)=ElementName(1:lenstr(ElementName))

      if (elementnumber.eq.1) then
        if (elementstatus(1).ne.2) then       
          write (XMLFileUnit, 95) ElementName(1:lenstr(ElementName))
          elementstatus(elementnumber) = 2
        endif
        return
      endif

c     close previous element if needed
      if (elementnumber.gt.1) then
        if (elementstatus(elementnumber-1).eq.1) then
c     add the closing bracket and print
          write (XMLOutputLine(lenstr(XMLOutputLine)+1:),100)
          write (XMLFileUnit,'(a)') 
     $              XMLOutputLine(1:lenstr(XMLOutputLine))
          XMLOutputLine = ' '
c     change status to complete and open
          elementstatus(elementnumber-1)=2
        endif
      endif
c
c     open current element
c     firstly indent line
      do 10 indent = 1, elementnumber
        write (indentline(indent:indent), 110)
 10   continue
c     write indent and tag
      write (XMLOutputLine,130) indentline(1:elementnumber),
     $     ElementName(1:lenstr(ElementName))

      elementstatus(elementnumber)=1

c unindented tag
 95   format(' <',a,'>')
c tag end and start...
 100  format('>')
 105  format('<')
c this is the indent
 110  format(' ')
c simply a string...
 120  format(a)
c indented open tag
 130  format(a,'<',a)
      return
      end

C     subroutine XMLWriteAttribute 
C     ElementName - if element not already opened then it will be
C     these are strings: AttributeName, AttributeValue,
C
      subroutine XMLWriteAttribute(XMLFileUnit, ElementName,
     $     AttributeName, AttributeValue, ifail)

      integer lenstr

      character*(*) ElementName, AttributeName, AttributeValue

      integer indent

      integer XMLFileUnit, ifail, elementstatus(80), elementnumber 
      character*80 openelements(80), XMLFileName
      character*1024 XMLOutputLine
      
      common /ccp4xmld/ elementnumber, elementstatus, openelements,
     $     XMLFileName,XMLOutputLine
c      save /ccp4xml/

      external lenstr

c     firstly check element is open otherwise open!
      if (elementstatus(elementnumber).ne.1) then
        call XMLOpenElement(XMLFileUnit, ElementName, ifail)
      endif 

c     secondly write the attribute 
      write (XMLOutputLine(lenstr(XMLOutputLine)+1:), 120) 
     $     AttributeName(1:lenstr(AttributeName)),
     $     AttributeValue(1:lenstr(AttributeValue))


c this is the indent
 110  format(' ')
c simply a string...
 120  format('  ',a,'="',a,'" ')

      return
      end

C     subroutine XMLWriteElement
C     for inputing the value of the element. if its not already opened
C     it will be. Closes the element as well
C
      subroutine XMLWriteElement(XMLFileUnit, ElementName,
     $     ElementValue, ifail)
      
      character*(*) ElementName, ElementValue

      integer XMLFileUnit, ifail, elementstatus(80), elementnumber 
      character*80 openelements(80), XMLFileName
      character*1024 XMLOutputLine
      
      common /ccp4xmld/ elementnumber, elementstatus, openelements,
     $     XMLFileName,XMLOutputLine
c      save /ccp4xml/

      integer lenstr
      external lenstr

      if (ElementName(1:lenstr(ElementName)).ne.
     $     openelements(elementnumber)) then
        call XMLOpenElement(XMLFileUnit, ElementName, ifail)
      endif

      if (elementnumber.ne.1) then
c     complete element and give the value
        write (XMLOutputLine(lenstr(XMLOutputLine)+1:), 125) 
     $               ElementValue(1:lenstr(ElementValue))
        elementstatus(elementnumber)=2
      endif

      call XMLCloseElement(XMLFileUnit, ElementName, ifail)
      

c an indented string
 125  format('>',a)

      return
      end

c subroutine XMLCloseElement
C the subroutine will close the element if it is the current element
C or go up the tree closing all elements till it reaches teh one specified.
C if its the root element being specified then it will close the file as well.
C
      subroutine XMLCloseElement(XMLFileUnit, ElementName, ifail)

      character*(*) ElementName

      character*80 indentline
      integer indent, elementlength

      integer XMLFileUnit, ifail, elementstatus(80), elementnumber 
      character*80 openelements(80), XMLFileName
      
      integer lenstr
      external lenstr
      character*1024 XMLOutputLine

      common /ccp4xmld/ elementnumber, elementstatus, openelements,
     $     XMLFileName,XMLOutputLine
c      save /ccp4xml/

c  sanity check to see if element is open
      do i = elementnumber,1,-1
        if (ElementName(1:lenstr(ElementName)).eq.
     $     openelements(i)) goto 1
      enddo
      call ccperr(2,'XMLCloseElement: trying to close wrong element')
      return

 1    continue

      elementlength=lenstr(openelements(elementnumber))

      indentline=' '
      indent=0

      do 10 indent = 1, elementnumber
        write (indentline(indent:indent), 110)
 10   continue

      if (ElementName(1:lenstr(ElementName)).ne.
     $     openelements(elementnumber)) then
C     if the element to close is not the open one then close 
C     the open one - can assume its just e.g. <wibble and not
C     <wibble> so only need a />
        write (XMLOutputLine(lenstr(XMLOutputLine)+1:), 155)
        write (XMLFileUnit,'(a)') 
     $              XMLOutputLine(1:lenstr(XMLOutputLine))
        XMLOutputLine = ' '
        openelements(elementnumber)=' '
        elementstatus(elementnumber)=0
        elementnumber=elementnumber-1
        goto 1
      else
        if (elementstatus(elementnumber).eq.1) then
          write (XMLOutputLine(lenstr(XMLOutputLine)+1:), 155)
          elementstatus(elementnumber)=2
          write (XMLFileUnit,'(a)') 
     $              XMLOutputLine(1:lenstr(XMLOutputLine))
          XMLOutputLine = ' '
        else
          if (lenstr(XMLOutputLine).eq.0) then
            write (XMLOutputLine,145) indentline(1:elementnumber),
     $            openelements(elementnumber)(1:elementlength)
          else
            write (XMLOutputLine(lenstr(XMLOutputLine)+1:), 150)
     $         openelements(elementnumber)(1:elementlength)
          endif
          write (XMLFileUnit,'(a)') 
     $              XMLOutputLine(1:lenstr(XMLOutputLine))
          XMLOutputLine = ' '
        endif
        openelements(elementnumber)=' '
        elementstatus(elementnumber)=0
        elementnumber=elementnumber-1
      endif

      if (elementnumber.eq.0) then
        close (XMLFileUnit)
      endif

c this is the indent
 110  format(' ')
c indented complete end tag
 145  format(a,'</',a,'>')
c indented complete end tag
 150  format('</',a,'>')
c indented complete end tag
 155  format(' />')

      return
      end 






