#include <fem.hpp> // Fortran EMulation library of fable module

namespace placeholder_please_replace {

using namespace fem::major_types;

void
ccp4h_init_clib(...)
{
  throw std::runtime_error(
    "Missing function implementation: ccp4h_init_clib");
}

void
ccppnm(...)
{
  throw std::runtime_error(
    "Missing function implementation: ccppnm");
}

int
lenstr(...)
{
  throw std::runtime_error(
    "Missing function implementation: lenstr");
}

void
lunsto(...)
{
  throw std::runtime_error(
    "Missing function implementation: lunsto");
}

void
parsekeyarg(...)
{
  throw std::runtime_error(
    "Missing function implementation: parsekeyarg");
}

void
parsesubkey(...)
{
  throw std::runtime_error(
    "Missing function implementation: parsesubkey");
}

void
ugtenv(...)
{
  throw std::runtime_error(
    "Missing function implementation: ugtenv");
}

struct common_ccp4hdat
{
  int lpt;
  bool html;
  bool logsumm;
  fem::str<160> cbin;
  fem::str<160> chtml;
  fem::str<160> cpid;
  int htmlinit;
  bool htmlopen;
  bool summopen;
  int summlevel;

  common_ccp4hdat() :
    lpt(fem::int0),
    html(fem::bool0),
    logsumm(fem::bool0),
    cbin(fem::char0),
    chtml(fem::char0),
    cpid(fem::char0),
    htmlinit(fem::int0),
    htmlopen(fem::bool0),
    summopen(fem::bool0),
    summlevel(fem::int0)
  {}
};

struct common :
  fem::common,
  common_ccp4hdat
{
  fem::cmn_sve ccp4h_init_lib_sve;

  common(
    int argc,
    char const* argv[])
  :
    fem::common(argc, argv)
  {}
};

struct ccp4h_init_lib_save
{
};

void
ccp4h_init_lib(
  common& cmn,
  int const& ihtml,
  int const& isumm)
{
  FEM_CMN_SVE(ccp4h_init_lib);
  // COMMON ccp4hdat
  bool& html = cmn.html;
  bool& logsumm = cmn.logsumm;
  fem::str<160>& cbin = cmn.cbin;
  fem::str<160>& chtml = cmn.chtml;
  fem::str<160>& cpid = cmn.cpid;
  int& htmlinit = cmn.htmlinit;
  //
  if (is_called_first_time) {
    htmlinit = -1;
    cmn.htmlopen = false;
  }
  //C   ccp4h_init_lib(ihtml,isum) - initialise variables in common block
  //C
  //C   Arguments: ihtml = 0 (use default settings for html tags)
  //C                     -1 (switch off html tags)
  //C              isumm = 0 (use default settings for summary tags)
  //C                     -1 (switch off summary tags)
  //C
  int idum = fem::int0;
  fem::str<160> dummy = fem::char0;
  int ihtml_level = fem::int0;
  int isumm_level = fem::int0;
  if (htmlinit < 0) {
    cmn.lpt = lunsto(idum);
    cbin = " ";
    chtml = " ";
    cpid = " ";
    dummy = " ";
    ugtenv("CBIN", cbin);
    ugtenv("CHTML", chtml);
    ugtenv("CCP_PROGRAM_ID", cpid);
    ugtenv("CCP_SUPPRESS_HTML", dummy);
    if (ihtml < 0) {
      html = false;
    }
    else {
      html = (dummy == " ");
    }
    ugtenv("CCP_SUPPRESS_SUMMARY", dummy);
    if (isumm < 0) {
      logsumm = false;
    }
    else {
      logsumm = (dummy == " ");
    }
    ihtml_level = 0;
    if (html) {
      ihtml_level = 1;
    }
    isumm_level = 0;
    if (logsumm) {
      isumm_level = 1;
    }
    ccp4h_init_clib(ihtml_level, isumm_level);
    htmlinit = 0;
    cmn.summopen = false;
    cmn.summlevel = 0;
  }
}

void
ccp4h_summary_beg(
  common& cmn)
{
  common_write write(cmn);
  // COMMON ccp4hdat
  int& lpt = cmn.lpt;
  bool& summopen = cmn.summopen;
  int& summlevel = cmn.summlevel;
  //
  //C   ccp4h_summary_beg() - begin summary section
  ccp4h_init_lib(cmn, 0, 0);
  //C   logsumm should be true unless CCP_SUPPRESS_SUMMARY set
  if (!cmn.logsumm) {
    return;
  }
  //C
  //C   if summopen is true then there is already an "unclosed" BEGIN tag
  //C   in this case don't write another
  //C   keep track of nested level of calls to ccp4h_summary_beg and
  //C   ccp4h_summary_end
  if (summopen) {
    summlevel++;
    return;
  }
  else {
    summlevel = 1;
  }
  if (cmn.html) {
    write(lpt, "('<B><FONT COLOR=\"#FF0000\"><!--SUMMARY_BEGIN-->')");
  }
  else {
    write(lpt, "('<!--SUMMARY_BEGIN-->')");
  }
  //C   set summopen to indicate there is an "unclosed" BEGIN
  summopen = true;
}

void
ccp4h_summary_end(
  common& cmn)
{
  common_write write(cmn);
  // COMMON ccp4hdat
  int& lpt = cmn.lpt;
  bool& summopen = cmn.summopen;
  int& summlevel = cmn.summlevel;
  //
  //C   ccp4h_summary_end() - end summary section
  if (!cmn.logsumm) {
    return;
  }
  //C   if summopen is not true then there is no matching BEGIN tag
  //C   in this case don't write an END
  //C   keep track of nested level of calls to ccp4h_summary_beg and
  //C   ccp4h_summary_end
  if (!summopen) {
    summlevel = 0;
    return;
  }
  if (summlevel > 1) {
    summlevel = summlevel - 1;
    return;
  }
  if (cmn.html) {
    write(lpt, "('<!--SUMMARY_END--></FONT></B>')");
  }
  else {
    write(lpt, "('<!--SUMMARY_END-->')");
  }
  //C   set summopen and summlevel to indicate there is no unmatched
  //C   BEGIN tag
  summopen = false;
  summlevel = 0;
}

void
ccp4h_rule(
  common& cmn)
{
  common_write write(cmn);
  // COMMON ccp4hdat
  int& lpt = cmn.lpt;
  //
  //C   ccp4h_rule() - rule (html <hr> tag)
  if (cmn.html) {
    write(lpt, "('<hr>')");
  }
  else {
    write(lpt, "('-------------------------------------------------------')");
  }
}

//C
//C     libhtml.f: write HTML tags into program output
//C     Copyright (C) 1998  Kevin Cowtan et al
//C
//C     This library is free software: you can redistribute it and/or
//C     modify it under the terms of the GNU Lesser General Public License
//C     version 3, modified in accordance with the provisions of the
//C     license to address the requirements of UK law.
//C
//C     You should have received a copy of the modified GNU Lesser General
//C     Public License along with this library.  If not, copies may be
//C     downloaded from http://www.ccp4.ac.uk/ccp4license.php
//C
//C     This program is distributed in the hope that it will be useful,
//C     but WITHOUT ANY WARRANTY; without even the implied warranty of
//C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//C     GNU Lesser General Public License for more details.
//C
//C     November 1998 KDC
//C     WROTE IT!
//C     Autumn 1999 SPE, AWA, MDW
//C     New routines ccp4h_init_lib, ccp4h_html_close, ccp4h_summary_beg,
//C                  ccp4h_summary_end
//C
//C    =========
//C    libhtml.f
//C    =========
//C
//C Link destinations can be of the for "#name" for a link in the same
//C document. Any other link will be prepended with $CHTML to link to
//C the documentation area
//C
//C The contents of the environment variable $CCP_PROGRAM_ID is concatenated
//C to the end of each link name. This allows multiple runs of a single
//C program from a script to be combined in one file without the links
//C conflicting.
//C
//C The environment variable $CCP_SUPPRESS_HTML if set suppresses the output
//C of most HTML tags. This is useful for devotees of plain text, and
//C especially for the command input section.
//C
//C   ccp4h_init_lib(int ihtml, int isumm) - initialise the routines for
//C     writing tags (html and summary)
//C
//C   ccp4h_init() - write initial comment to identify file
//C
//C   ccp4h_html_close() - write html close tag
//C
//C   ccp4h_toc_beg() - write starting Contents section tags
//C
//C   ccp4h_toc_ent(char*(*) text, char*(*) dest) - write Contents entry
//C     text= the link text
//C     dest= the link destination
//C
//C   ccp4h_toc_end() - write ending Contents section tags
//C
//C   ccp4h_graph_beg(x,y) - write starting JLogGraph applet tag
//C     x = width of graph in pixels
//C     y = height of graph in pixels
//C      both can be zero, then they default to 400,300.
//C
//C   ccp4h_graph_end() - write ending JLogGraph applet tag
//C
//C   ccp4h_summary_beg() - begin summary section
//C
//C   ccp4h_summary_end() - end summary section
//C
//C   ccp4h_pre_beg() - begin preformatted (html <pre> tag)
//C
//C   ccp4h_pre_end() - end preformatted (html <pre> tag)
//C
//C   ccp4h_rule() - rule (html <hr> tag)
//C
//C   ccp4h_link(char*(*) text, char*(*) dest) - link (html <a> tag)
//C     text= the link text
//C     dest= the link destination
//C       if dest is not an anchor name (i.e. begins with '#'), then
//C       the $CHTML path is prefixed automatically.
//C
//C   ccp4h_link_key(char*(*) text, char*(*) dest) - link to CCP4 documentation
//C     text= the keyparser keyword
//C     dest= the link destination filename and anchor,
//C       the $CHTML path is prefixed automatically.
//C     Between call memparse(.false.) and the main keyparser processing,
//C     call this routine with any keywords and links which should be
//C     linked to the program doc. Then call parsefail() to restart
//C     syntax checking.... e.g.
//C      call memoparse(.false)
//C      call ccp4h_link_key('LABIN   ','fft.html#labin')
//C      call ccp4h_link_key('GRID    ','fft.html#grid')
//C      call ccp4h_link_key('FFTSPGRP','fft.html#fftspgrp')
//C      call parsefail()
//C
//C   ccp4h_header(char*(*) text, char*(*) name, int level) - header
//C     text= the link text
//C     name= header name to link to
//C     level=0-6. Header size. 0= plain text
//C
void
ccp4h_init(
  common& cmn)
{
  common_write write(cmn);
  // COMMON ccp4hdat
  bool& htmlopen = cmn.htmlopen;
  //
  //C   ccp4h_init() - write initial comment to identify file
  ccp4h_init_lib(cmn, 0, 0);
  if (cmn.html && !htmlopen) {
    ccp4h_summary_beg(cmn);
    write(cmn.lpt, "('<html> <!-- CCP4 HTML LOGFILE -->')");
    htmlopen = true;
    ccp4h_rule(cmn);
    ccp4h_summary_end(cmn);
  }
}

void
ccp4h_html_close(
  common& cmn)
{
  common_write write(cmn);
  //C   ccp4h_html_close() - write html close tag
  if (cmn.html) {
    ccp4h_summary_beg(cmn);
    write(cmn.lpt, "('</html>')");
    cmn.htmlopen = false;
    ccp4h_summary_end(cmn);
  }
}

void
ccp4h_header(
  common& cmn,
  str_cref text,
  str_cref name,
  int const& level)
{
  common_write write(cmn);
  // COMMON ccp4hdat
  int& lpt = cmn.lpt;
  fem::str<160>& cpid = cmn.cpid;
  //
  //C   ccp4h_header(char*(*) text, char*(*) name, int level) - header
  //C     text= the link text
  //C     name= header name to link to
  //C     level=0-6. Header size. 0= plain text
  fem::str<160> pn = fem::char0;
  ccppnm(pn);
  fem::str<60> underline = fem::char0;
  if (cmn.html) {
    if (level > 0) {
      write(lpt,
        "(/,'<a name=\"',a,a,a,'\"><h',i1,'>',a,'</h',i1,'></a>')"),
        name, pn(1, lenstr(pn)), cpid(1, lenstr(cpid)), level, text,
        level;
    }
    else {
      write(lpt, "(/,'<a name=\"',a,a,a,'\">',a,'</a>')"), name, pn(1,
        lenstr(pn)), cpid(1, lenstr(cpid)), text;
    }
  }
  else {
    if (level <= 0) {
      write(lpt, "(a)"), text(1, lenstr(text));
    }
    else if (level == 1 || level == 2 || level == 3) {
      underline =
        "------------------------------------------------------------";
      write(lpt, "(/,/,a,/,a,/)"), text(1, lenstr(text)), underline(1,
        fem::min(lenstr(text), 60));
    }
    else {
      write(lpt, "(/,a,/)"), text(1, lenstr(text));
    }
  }
}

void
ccp4h_toc_beg(
  common& cmn)
{
  common_write write(cmn);
  //C   ccp4h_toc_beg() - write starting Contents section tags
  ccp4h_header(cmn, "Contents", "toc", 2);
  if (cmn.html) {
    write(cmn.lpt, "('<ul>')");
  }
}

void
ccp4h_toc_ent(
  common& cmn,
  str_cref text,
  str_cref dest)
{
  common_write write(cmn);
  // COMMON ccp4hdat
  int& lpt = cmn.lpt;
  fem::str<160>& cpid = cmn.cpid;
  //
  //C   ccp4h_toc_ent(char*(*) text, char*(*) dest) - write Contents entry
  //C     text= the link text
  //C     dest= the link destination
  fem::str<160> pn = fem::char0;
  ccppnm(pn);
  if (cmn.html) {
    write(lpt, "('<li><a href=\"',a,a,a,'\">',a,'</a>')"), dest, pn(1,
      lenstr(pn)), cpid(1, lenstr(cpid)), text;
  }
  else {
    write(lpt, "(a)"), text;
  }
}

void
ccp4h_toc_end(
  common& cmn)
{
  common_write write(cmn);
  //C   ccp4h_toc_end() - write ending Contents section tags
  if (cmn.html) {
    write(cmn.lpt, "('</ul>')");
  }
}

void
ccp4h_graph_beg(
  common& cmn,
  int const& x,
  int const& y)
{
  common_write write(cmn);
  // COMMON ccp4hdat
  fem::str<160>& cbin = cmn.cbin;
  //
  //C   ccp4h_graph_beg() - write starting JLogGraph applet tag
  //C     x = width of graph in pixels
  //C     y = height of graph in pixels
  //C      both can be zero, then they default to 400,300.
  int x1 = x;
  int y1 = y;
  if (x1 <= 0) {
    x1 = 400;
  }
  if (y1 <= 0) {
    y1 = 300;
  }
  if (cmn.html) {
    write(cmn.lpt,
      "('<applet width=\"',i4,'\" height=\"',i4,"
      "'\" code=\"JLogGraph.class\" ',/,'codebase=\"',a,"
      "'\"><param name=\"table\" value=\"')"),
      x1, y1, cbin(1, lenstr(cbin));
  }
}

void
ccp4h_graph_end(
  common& cmn)
{
  common_write write(cmn);
  //C   ccp4h_graph_end() - write ending JLogGraph applet tag
  if (cmn.html) {
    write(cmn.lpt,
      "('\"><b>For inline graphs use a Java browser</b></applet>')");
  }
}

void
ccp4h_pre_beg(
  common& cmn)
{
  common_write write(cmn);
  //C   ccp4h_pre_beg() - begin preformatted (html <pre> tag)
  if (cmn.html) {
    write(cmn.lpt, "('<pre>')");
  }
}

void
ccp4h_pre_end(
  common& cmn)
{
  common_write write(cmn);
  //C   ccp4h_pre_end() - end preformatted (html <pre> tag)
  if (cmn.html) {
    write(cmn.lpt, "('</pre>')");
  }
}

void
ccp4h_link(
  common& cmn,
  str_cref text,
  str_cref dest)
{
  common_write write(cmn);
  // COMMON ccp4hdat
  int& lpt = cmn.lpt;
  fem::str<160>& chtml = cmn.chtml;
  fem::str<160>& cpid = cmn.cpid;
  //
  //C   ccp4h_link(char*(*) text, char*(*) dest) - link (html <a> tag)
  //C     text= the link text
  //C     dest= the link destination
  fem::str<160> pn = fem::char0;
  ccppnm(pn);
  if (cmn.html) {
    if (dest(1, 1) == "#") {
      write(lpt, "('<a href=\"',a,a,a,'\">',a,'</a>')"), dest, pn(1,
        lenstr(pn)), cpid(1, lenstr(cpid)), text;
    }
    else {
      write(lpt, "('<a href=\"',a,'/',a,'\">',a,'</a>')"), chtml(1,
        lenstr(chtml)), dest, text;
    }
  }
  else {
    write(lpt, "(a)"), text;
  }
}

void
ccp4h_link_key(
  common& cmn,
  str_cref key,
  str_cref dest)
{
  common_write write(cmn);
  // COMMON ccp4hdat
  int& lpt = cmn.lpt;
  fem::str<160>& chtml = cmn.chtml;
  //
  //C   ccp4h_link(char*(*) text, char*(*) dest) - link keyword
  //C     text= the link text
  //C     dest= the link destination
  bool flag = false;
  fem::str<4> kw = key;
  parsesubkey(kw, "    ", flag);
  fem::str<120> rest = fem::char0;
  if (flag) {
    parsekeyarg(kw, rest);
    if (cmn.html) {
      write(lpt,
        "(' Data line--- <a href=\"',a,'/',a,'\">',a,'</a> ',a)"),
        chtml(1, lenstr(chtml)), dest(1, lenstr(dest)), key, rest(1,
        lenstr(rest));
    }
    else {
      write(lpt, "(' Data line--- ',a,a)"), key, rest(1, lenstr(rest));
    }
  }
}

} // namespace placeholder_please_replace
