#include <fem.hpp> // Fortran EMulation library of fable module

namespace placeholder_please_replace {

using namespace fem::major_types;

void
ccpdpn(...)
{
  throw std::runtime_error(
    "Missing function implementation: ccpdpn");
}

void
ccperr(...)
{
  throw std::runtime_error(
    "Missing function implementation: ccperr");
}

void
ccpfyp(...)
{
  throw std::runtime_error(
    "Missing function implementation: ccpfyp");
}

int
ccpnun(...)
{
  throw std::runtime_error(
    "Missing function implementation: ccpnun");
}

void
ccprcs(...)
{
  throw std::runtime_error(
    "Missing function implementation: ccprcs");
}

void
cunlink(...)
{
  throw std::runtime_error(
    "Missing function implementation: cunlink");
}

bool
litend(...)
{
  throw std::runtime_error(
    "Missing function implementation: litend");
}

void
nocrlf(...)
{
  throw std::runtime_error(
    "Missing function implementation: nocrlf");
}

void
qclose(...)
{
  throw std::runtime_error(
    "Missing function implementation: qclose");
}

bool
qisnan(...)
{
  throw std::runtime_error(
    "Missing function implementation: qisnan");
}

void
qmode(...)
{
  throw std::runtime_error(
    "Missing function implementation: qmode");
}

void
qnan(...)
{
  throw std::runtime_error(
    "Missing function implementation: qnan");
}

void
qopen(...)
{
  throw std::runtime_error(
    "Missing function implementation: qopen");
}

void
qqinq(...)
{
  throw std::runtime_error(
    "Missing function implementation: qqinq");
}

void
qread(...)
{
  throw std::runtime_error(
    "Missing function implementation: qread");
}

void
qseek(...)
{
  throw std::runtime_error(
    "Missing function implementation: qseek");
}

void
qwrite(...)
{
  throw std::runtime_error(
    "Missing function implementation: qwrite");
}

void
ubytes(...)
{
  throw std::runtime_error(
    "Missing function implementation: ubytes");
}

void
ucputm(...)
{
  throw std::runtime_error(
    "Missing function implementation: ucputm");
}

void
ugerr(...)
{
  throw std::runtime_error(
    "Missing function implementation: ugerr");
}

void
ugtenv(...)
{
  throw std::runtime_error(
    "Missing function implementation: ugtenv");
}

void
ugtuid(...)
{
  throw std::runtime_error(
    "Missing function implementation: ugtuid");
}

void
uidate(...)
{
  throw std::runtime_error(
    "Missing function implementation: uidate");
}

void
uisatt(...)
{
  throw std::runtime_error(
    "Missing function implementation: uisatt");
}

void
ustime(...)
{
  throw std::runtime_error(
    "Missing function implementation: ustime");
}

void
utime(...)
{
  throw std::runtime_error(
    "Missing function implementation: utime");
}

bool
vaxvms(...)
{
  throw std::runtime_error(
    "Missing function implementation: vaxvms");
}

bool
winmvs(...)
{
  throw std::runtime_error(
    "Missing function implementation: winmvs");
}

struct common :
  fem::common
{
  fem::cmn_sve ranu_sve;
  fem::cmn_sve program_mctest_sve;

  common(
    int argc,
    char const* argv[])
  :
    fem::common(argc, argv)
  {}
};

struct ranu_save
{
  bool ff;
  arr<int> ir;
  int iy;

  ranu_save() :
    ff(fem::bool0),
    ir(dimension(97), fem::fill0),
    iy(fem::int0)
  {}
};

float
ranu(
  common& cmn,
  int& k)
{
  float return_value = fem::float0;
  FEM_CMN_SVE(ranu);
  // SAVE
  bool& ff = sve.ff;
  arr_ref<int> ir(sve.ir, dimension(97));
  int& iy = sve.iy;
  //
  if (is_called_first_time) {
    ff = true;
  }
  //C==== UNIFORM PSEUDO-RANDOM NUMBER IN THE RANGE >= 0 AND < 1.
  //C
  //C     Set the seed K zero or negative to start or restart the sequence.
  //C     K must be a variable since a new value is returned each time.
  //C
  const int ic = 150889;
  const int m = 714025;
  int j = fem::int0;
  const int ia = 1366;
  if (k <= 0 || ff) {
    ff = false;
    k = fem::mod(ic - k, m);
    FEM_DO_SAFE(j, 1, 97) {
      k = fem::mod(ia * k + ic, m);
      ir(j) = k;
    }
    k = fem::mod(ia * k + ic, m);
    iy = k;
  }
  //C
  j = 1 + 97 * iy / m;
  iy = ir(j);
  const float rm = 1.f / m;
  return_value = iy * rm;
  k = fem::mod(ia * k + ic, m);
  ir(j) = k;
  return return_value;
}

struct program_mctest_save
{
  static const int lstr = 12;

  int iloop;
  int iseed;
  int lunin;
  int lunout;
  fem::str<lstr> name1;
  fem::str<lstr> name2;
  fem::str<lstr> tstnam;

  program_mctest_save() :
    iloop(fem::int0),
    iseed(fem::int0),
    lunin(fem::int0),
    lunout(fem::int0),
    name1(fem::char0),
    name2(fem::char0),
    tstnam(fem::char0)
  {}
};

const int program_mctest_save::lstr;

//C
//C     testlib.f: test program for Machine dependent routines
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
void
program_mctest(
  int argc,
  char const* argv[])
{
  common cmn(argc, argv);
  FEM_CMN_SVE(program_mctest);
  common_read read(cmn);
  common_write write(cmn);
  int& iloop = sve.iloop;
  int& iseed = sve.iseed;
  int& lunin = sve.lunin;
  int& lunout = sve.lunout;
  const int lstr = 12;
  fem::str<lstr>& name1 = sve.name1;
  if (is_called_first_time) {
    lunin = 5;
    lunout = 6;
    name1 = "AAA.TST";
    sve.name2 = "ZZZ.TST";
    iseed = 0;
    iloop = 100;
    sve.tstnam = "TESTNAME";
  }
  float sec = fem::float0;
  int i = fem::int0;
  fem::str<120> envnam = fem::char0;
  fem::str<lstr> usrnam = fem::char0;
  int imon = fem::int0;
  int iday = fem::int0;
  int iyear = fem::int0;
  fem::str<lstr> udate = fem::char0;
  fem::str<lstr> usrtim = fem::char0;
  int isec = fem::int0;
  int iyes = fem::int0;
  fem::str<lstr> reply = fem::char0;
  int ibyte = fem::int0;
  fem::str<lstr> handle = fem::char0;
  int idum = fem::int0;
  fem::str<40> errstr = fem::char0;
  float nan = fem::float0;
  int lun = fem::int0;
  int ldum = fem::int0;
  const int lbuf = 500;
  int j = fem::int0;
  arr<float> buffer(dimension(lbuf), fem::fill0);
  int ilength = fem::int0;
  int istat = fem::int0;
  int nrec = fem::int0;
  int ier = fem::int0;
  fem::str<3> foo = fem::char0;
  static const char* format_6000 = "(/,/,a,/)";
  static const char* format_6002 = "(a7,3x,a12,a)";
  static const char* format_6014 =
    "(' Seek Record:',i5,'  Read:  ',f8.2,'  Status: ',i4)";
  static const char* format_6015 = "(' Seek Record:',i5,'  Write: ',f8.2)";
  static const char* format_6020 = "(i10)";
  //C
  //C---- Test file for Machine dependent routines
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Local Arrays ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Data statements ..
  //C     ..
  //C
  //C---- Initialise CPU timer
  //C
  sec = -1.0f;
  ucputm(sec);
  //C
  //C---- Parse command line arguments, open printer stream and print
  //C     version
  //C
  ccpfyp();
  i = 0;
  ccpdpn(lunout, "PRINTER", "PRINTER", "F", 0, i);
  if (i != 0) {
    ccperr(1, "Can't open printer stream");
  }
  ccprcs(6, "TESLIB", "$Date$");
  //C
  //C---- Other initialisations
  //C
  i = 0;
  ccpdpn(10, name1, "NEW", "F", 0, i);
  if (i != 0) {
    ccperr(1, " Failed to open file " + name1);
  }
  write(10, "(a)"), " I am the test file, delete me.";
  cmn.io.close(10)
    .status("DELETE");
  //C
  //C---- Start of Tests
  //C
  write(lunout, format_6000), " Routine  Result      Comments";
  //C
  //C---- UGTENV
  //C
  ugtenv("TESTENV", envnam);
  write(lunout, format_6002), " UGTENV", envnam(1, 12), "Get value of TESTENV";
  //C
  //C---- UGTUID
  //C
  ugtuid(usrnam);
  write(lunout, format_6002), " UGTUID", usrnam, "Get Users name";
  //C
  //C---- UIDATE
  //C
  uidate(imon, iday, iyear);
  write(udate, "(i2.2,'/',i2.2,'/',i4)"), iday, imon, iyear;
  write(lunout, format_6002), " UIDATE", udate, "Todays date";
  //C
  //C---- UTIME
  //C
  utime(usrtim);
  write(lunout, format_6002), " UTIME ", usrtim, "Current time";
  //C
  //C---- USTIME
  //C
  ustime(isec);
  write(udate, format_6020), isec;
  write(lunout, format_6002), " USTIME", udate, "Absolute time (VMS=-1)";
  //C
  //C---- UISATT
  //C
  uisatt(lunout, iyes);
  if (iyes == 0) {
    reply = "No";
  }
  else {
    write(reply, "(a,i3)"), "Yes, ", lunout;
  }
  write(lunout, format_6002), " UISATT", reply,
    "are we attached to at tty? Unit number?";
  //C
  //C---- VAXVMS
  //C
  reply = "No";
  if (vaxvms()) {
    reply = "Yes";
  }
  write(lunout, format_6002), " VAXVMS", reply, "Is this VMS?";
  //C
  //C---- WINMVS
  //C
  reply = "No";
  if (winmvs()) {
    reply = "Yes";
  }
  write(lunout, format_6002), " WINMVS", reply, "Is this Win NT et al?";
  //C
  //C---- UBYTES
  //C
  ubytes(ibyte, handle);
  write(reply, "(a5,i3)"), handle, ibyte;
  write(lunout, format_6002), " UBYTES", reply,
    "Get BYTE/WORD Handling and number of bytes per word";
  //C
  //C---- LITEND
  //C
  reply = "Big";
  if (litend(idum)) {
    reply = "Little";
  }
  write(lunout, format_6002), " LITEND", reply, "Big/Little end machine";
  //CCCC
  //CCCC---- URENAM
  //CCCC
  //CCC      CALL URENAM(NAME1,NAME2,ISTAT)
  //CCC      ERRSTR = 'OK'
  //CCC      IF (ISTAT.NE.0) CALL UGERR(ISTAT,ERSTR)
  //CCC      WRITE (LUNOUT,FMT=6002) ' URENAM',ERRSTR,'Check rename status'
  //CCC      CALL CUNLINK (NAME2)
  //C
  //C---- UGERR
  //C
  cmn.io.open(21, "TESTFILE")
    .status("OLD")
    .iostat(i);
  ugerr(i, errstr);
  write(reply, format_6020), i;
  write(lunout, format_6002), " UGERR ", reply, errstr;
  //C
  //C---- UCPUTM
  //C
  sec = 99.99f;
  ucputm(sec);
  write(reply, "(f8.2)"), sec;
  write(lunout, format_6002), " UCPUTM", reply, "Show elapsed CPU time";
  //C
  //C---- NOCRLF
  //C
  nocrlf("NOCRLF");
  write(lunout, "('+',14x,a)"), "Should be on same line";
  //C
  //C --- CCPNUN
  //C
  i = ccpnun();
  write(reply, format_6020), i;
  write(lunout, format_6002), " CCPNUN", reply, "Next free unit";
  //C
  //C --- QNAN/QISNAN (`magic numbers')
  //C
  qnan(nan);
  if ((!qisnan(nan)) || qisnan(1.0f)) {
    write(lunout, "(/,' *** QNAN/QISNAN test failed',/)");
  }
  else {
    write(lunout, "(' QNAN/QISNAN test OK')");
  }
  //C
  //C---- End of tests
  //C
  write(lunout, format_6000), " Now test diskio routines";
  //C
  //C---- Now test the diskio stuff
  //C
  qopen(lun, "DISKIO", "UNKNOWN");
  qmode(lun, 2, ldum);
  //C
  //C---- Write a file of size LBUF x LBUF x WORDSIZE
  //C
  FEM_DO_SAFE(i, 1, lbuf) {
    FEM_DO_SAFE(j, 1, lbuf) {
      buffer(j) = 0.f;
    }
    buffer(i) = i;
    qwrite(lun, buffer.begin(), lbuf);
  }
  //C
  //C---- Close the file
  //C
  qclose(lun);
  //C
  //C---- reset the array buffer(*)
  //C
  FEM_DO_SAFE(i, 1, lbuf) {
    buffer(i) = 0.0f;
  }
  //C
  //C---- Now do some reads on the file just created
  //C
  qopen(lun, "DISKIO", "OLD");
  qmode(lun, 2, ldum);
  //C
  //C---- test file size
  //C
  qqinq(lun, "DISKIO", reply, ilength);
  istat = ldum * lbuf * lbuf;
  write(6, "(a,2(i8,a),/,/)"), " DISKIO should be", istat,
    " bytes; is", ilength, " bytes.";
  if (ilength != istat) {
    ccperr(1, "*** FILE SIZE ERROR ***");
  }
  //C
  //C---- Seed random Number Generator
  //C
  reply = " ";
  ugtenv("SEED", reply);
  if (reply != " ") {
    read(reply, star), iseed;
  }
  //C
  //C---- Get number of reads to perform
  //C
  ugtenv("READS", reply);
  if (reply != " ") {
    read(reply, star), iloop;
  }
  //C
  //C---- Do random reads & writes on the file
  //C
  ugtenv("NUMRECORD", envnam);
  if (envnam == " ") {
    FEM_DO_SAFE(i, 1, iloop) {
      nrec = fem::nint(100.f * ranu(cmn, iseed) + 1.f);
      if (nrec <= 0 || nrec > lbuf) {
        ccperr(1, "*** RECORD ERROR ***");
      }
      qseek(lun, nrec, 1, lbuf);
      if (ranu(cmn, iseed) < .5f) {
        qread(lun, buffer.begin(), lbuf, ier);
        write(lunout, format_6014), nrec, buffer(nrec), ier;
        if (buffer(nrec) != nrec) {
          ccperr(1, "*** VERIFY ERROR ***");
        }
      }
      else {
        FEM_DO_SAFE(j, 1, lbuf) {
          buffer(j) = 0.f;
        }
        buffer(nrec) = nrec;
        qwrite(lun, buffer.begin(), lbuf);
        write(lunout, format_6015), nrec, buffer(nrec);
      }
    }
  }
  else {
    statement_50:
    nocrlf("Record to seek (-ve to write, Ctrl/Z to stop)> ");
    try {
      read(lunin, format_6020), nrec;
    }
    catch (fem::read_end const&) {
      goto statement_60;
    }
    i = fem::iabs(nrec);
    if (i == 0 || i > lbuf) {
      write(lunout, star), "*** RECORD ERROR ***";
      goto statement_50;
    }
    qseek(lun, i, 1, lbuf);
    if (nrec > 0) {
      qread(lun, buffer.begin(), lbuf, ier);
      write(lunout, format_6014), nrec, buffer(nrec), ier;
      if (buffer(nrec) != nrec) {
        write(lunout, star), "*** VERIFY ERROR ***";
      }
    }
    else {
      FEM_DO_SAFE(j, 1, lbuf) {
        buffer(j) = 0.f;
      }
      buffer(i) = i;
      qwrite(lun, buffer.begin(), lbuf);
      write(lunout, format_6015), i, buffer(i);
    }
    goto statement_50;
  }
  statement_60:
  qclose(lun);
  cunlink("DISKIO");
  //C     Now check we can open and close a scratch file
  qopen(lun, "foo.bar", "SCRATCH");
  qclose(lun);
  //C     and can we rewind a scratch file?  (make sure something's been
  //C     written to it first)
  i = 0;
  ccpdpn(lun, "FOOEY", "SCRATCH", "F", 0, i);
  write(lun, "(a)"), "foo";
  try {
    cmn.io.rewind(lun);
  }
  catch (fem::io_err const&) {
    goto statement_170;
  }
  read(lun, "(a)"), foo;
  write(lunout, star), "contents of temp file: ", foo;
  ccperr(0, "Normal Termination");
  statement_170:
  ccperr(1, "Can't rewind scratch file");
  ccperr(1, "*** EOF ERROR ***");
  //C
  //C---- Format Statements
  //C
}

} // namespace placeholder_please_replace

int
main(
  int argc,
  char const* argv[])
{
  return fem::main_with_catch(
    argc, argv,
    placeholder_please_replace::program_mctest);
}
