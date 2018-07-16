check_PROGRAMS = rfltest sgtest symtest

if BUILD_CCP4
check_PROGRAMS += ccad diffmap maptest mtztest
if BUILD_CONTRIB
check_PROGRAMS += cmaplocal ffttest sktest
endif
if BUILD_CNS
check_PROGRAMS += cnstest
endif
if BUILD_MMDB
check_PROGRAMS += cns2mtz maketestdata sftest
endif
if BUILD_PHS
check_PROGRAMS += phstest
endif
endif

if BUILD_CIF
if BUILD_PHS
check_PROGRAMS += ciftophs
endif
endif

if BUILD_FORTRAN
check_PROGRAMS += ftndemo
endif

if BUILD_CCTBX
check_PROGRAMS += cctbxtest
endif

ccad_SOURCES = ccad.cpp
cctbxtest_SOURCES = cctbxtest.cpp
ciftophs_SOURCES = ciftophs.cpp
cnstest_SOURCES = cnstest.cpp
cmaplocal_SOURCES = cmaplocal.cpp
cns2mtz_SOURCES = cns2mtz.cpp cns2mtz.h
diffmap_SOURCES = diffmap.cpp
ffttest_SOURCES = ffttest.cpp
ftndemo_SOURCES = ftndemo.f
ftndemo_LINK = $(CXXLINK)
maketestdata_SOURCES = maketestdata.cpp
maptest_SOURCES = maptest.cpp
mtztest_SOURCES = mtztest.cpp
phstest_SOURCES = phstest.cpp
rfltest_SOURCES = rfltest.cpp
sftest_SOURCES = sftest.cpp
sgtest_SOURCES = sgtest.cpp
sktest_SOURCES = sktest.cpp
symtest_SOURCES = symtest.cpp

AM_CPPFLAGS = -I$(top_srcdir)
LIBS = -L$(libdir)

# LDADD is used for programs for which program_LDADD is not set
LDADD = ../clipper/core/libclipper-core.la $(FFTW2_LIBS)

ADD_CLIPPER_CCP4 = ../clipper/ccp4/libclipper-ccp4.la
ADD_CLIPPER_MMDB = ../clipper/mmdb/libclipper-mmdb.la
ADD_CLIPPER_CIF = ../clipper/cif/libclipper-cif.la
ADD_CLIPPER_CONTRIB = ../clipper/contrib/libclipper-contrib.la

ccad_LDADD = $(ADD_CLIPPER_CCP4) $(LDADD)
cctbxtest_LDADD = ../clipper/cctbx/libclipper-cctbx.la $(ADD_CLIPPER_CONTRIB) $(LDADD)
ciftophs_LDADD = ../clipper/phs/libclipper-phs.la $(ADD_CLIPPER_CIF) $(LDADD)
cmaplocal_LDADD = $(ADD_CLIPPER_CCP4) $(ADD_CLIPPER_CONTRIB) $(LDADD)
cns2mtz_LDADD = $(ADD_CLIPPER_CCP4) $(ADD_CLIPPER_MMDB) $(LDADD)
cnstest_LDADD = ../clipper/cns/libclipper-cns.la $(ADD_CLIPPER_CCP4) $(LDADD)
diffmap_LDADD = $(ADD_CLIPPER_CCP4) $(LDADD)
ffttest_LDADD = $(ADD_CLIPPER_CCP4) $(ADD_CLIPPER_CONTRIB) $(LDADD)
ftndemo_LDADD = ../clipper/fortran/libclipper-fortran.la $(LDADD) $(FLIBS)
maketestdata_LDADD = $(ADD_CLIPPER_MMDB) $(ADD_CLIPPER_CCP4) $(ADD_CLIPPER_CONTRIB) $(LDADD)
maptest_LDADD = $(ADD_CLIPPER_CCP4) $(LDADD)
mtztest_LDADD = $(ADD_CLIPPER_CCP4) $(LDADD)
phstest_LDADD = ../clipper/phs/libclipper-phs.la $(ADD_CLIPPER_CCP4) $(LDADD)
sftest_LDADD = $(ADD_CLIPPER_CCP4) $(ADD_CLIPPER_MMDB) $(ADD_CLIPPER_CONTRIB) $(LDADD)
sktest_LDADD = $(ADD_CLIPPER_CCP4) $(ADD_CLIPPER_CONTRIB) $(LDADD)

EXTRA_DIST = test.csh test.log testfile.mtz
