ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = ccp4/libccp4c.la

if FORTRAN_LIB
lib_LTLIBRARIES += fortran/libccp4f.la
endif

if GNU_FORTRAN
AM_FFLAGS = -fno-second-underscore
endif

ccp4datadir = $(datarootdir)/ccp4
ccp4includedir = $(includedir)/ccp4
pkgconfigdir = $(libdir)/pkgconfig

ccp4_libccp4c_la_SOURCES = \
ccp4/ccp4_array.c ccp4/cmap_data.c ccp4/cmtzlib.c ccp4/ccp4_general.c \
ccp4/cmap_header.c ccp4/csymlib.c ccp4/ccp4_parser.c ccp4/cmap_labels.c \
ccp4/cvecmat.c ccp4/ccp4_program.c ccp4/cmap_open.c ccp4/library_err.c \
ccp4/ccp4_unitcell.c ccp4/cmap_skew.c ccp4/library_file.c \
ccp4/cmap_accessor.c ccp4/cmap_stats.c ccp4/library_utils.c ccp4/cmap_close.c \
ccp4/cmap_symop.c ccp4/pack_c.c

# FIXME: is ccp4/w32mvs.h needed here?
#        is cmaplib_f.h a public header?
ccp4include_HEADERS = \
ccp4/ccp4_file_err.h  ccp4/ccp4_program.h   ccp4/ccp4_unitcell.h  \
ccp4/cmap_errno.h     ccp4/cmap_stats.h     ccp4/csymlib.h        \
ccp4/library_file.h   ccp4/w32mvs.h         ccp4/ccp4_fortran.h   \
ccp4/ccp4_spg.h       ccp4/ccp4_utils.h     ccp4/cmap_header.h    \
ccp4/cmaplib.h        ccp4/cvecmat.h        ccp4/mtzdata.h        \
ccp4/ccp4_array.h     ccp4/ccp4_general.h   ccp4/ccp4_sysdep.h    \
ccp4/ccp4_vars.h      ccp4/cmap_labels.h    ccp4/cmaplib_f.h      \
ccp4/overview.h       ccp4/ccp4_errno.h     ccp4/ccp4_parser.h    \
ccp4/ccp4_types.h     ccp4/cmap_data.h      ccp4/cmap_skew.h      \
ccp4/cmtzlib.h        ccp4/pack_c.h


fortran_libccp4f_la_SOURCES = \
fortran/ccp4_diskio_f.c fortran/ccp4_parser_f.c fortran/cmaplib_f.c \
fortran/csymlib_f.c fortran/ccp4_general_f.c fortran/ccp4_unitcell_f.c \
fortran/cmtzlib_f.c fortran/library_f.c \
 \
fortran/ccplib.f fortran/fftlib.f fortran/modlib.f fortran/parser.f \
fortran/rwbrook.f fortran/libhtml.f fortran/keyparse.f \
fortran/dna_output.f fortran/pack_f.f fortran/pxxml.f fortran/lgglib.f \
fortran/plot84lib.f fortran/plotsubs.f fortran/sorting_main.f \
fortran/dna_header.fh

nodist_fortran_libccp4f_la_SOURCES = fortran/unix.f

fortran_libccp4f_la_CPPFLAGS = -I$(srcdir)/ccp4
fortran_libccp4f_la_LIBADD = ccp4/libccp4c.la -lmmdb2

BUILT_SOURCES = fortran/unix.f
CLEANFILES = fortran/unix.f
fortran/unix.f: fortran/unix.m4 Makefile
	m4 $(M4_ARG) $< > $@

ccp4_libccp4c_la_LDFLAGS = -no-undefined
fortran_libccp4f_la_LDFLAGS = -no-undefined

# Files default.def and environ.def are not part of libccp4.
# They are included in CCP4 and may be needed for running programs
# compiled with libccp4 (see http://www.ccp4.ac.uk/dist/html/ccp4.html)

dist_ccp4data_DATA = ccp4-config.cmake
dist_data_DATA = data/atomsf_neutron.lib data/atomsf_electron.lib \
		     data/atomsf.lib data/syminfo.lib data/symop.lib
pkgconfig_DATA = libccp4c.pc
if FORTRAN_LIB
  dist_data_DATA += data/font84.ascii
  pkgconfig_DATA += libccp4f.pc
endif


# FIXME: all files but the first are unused now, should they be distributed?
EXTRA_DIST = fortran/unix.m4 \
ccp4/w32mvs.c ccp4/vmslibrary.c \
fortran/vmsdiskio.for fortran/vms.for \
fortran/testlib.f fortran/mthlib.for fortran/w32mvs.f

# Clean up Fortran modules.
clean-local: clean-modules
clean-modules:
	$(RM) *.mod
