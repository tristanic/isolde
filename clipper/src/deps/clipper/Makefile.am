
pkginclude_HEADERS = \
clipper-ccp4.h          clipper-contrib.h       clipper-phs.h \
clipper-cctbx.h         clipper-minimol.h       clipper.h \
clipper-cif.h           clipper-mmdb.h          minimal-clipper-hkl.h \
clipper-cns.h           minimal-clipper-map.h   imex.h

SUBDIRS = core

AM_CPPFLAGS=-I$(top_srcdir)

if BUILD_CONTRIB
SUBDIRS += contrib
endif

lib_LTLIBRARIES =

if BUILD_MMDB
lib_LTLIBRARIES += mmdb/libclipper-mmdb.la
mmdb_libclipper_mmdb_la_LDFLAGS = @VERSION_INFO@
mmdb_libclipper_mmdb_la_LIBADD = core/libclipper-core.la -lmmdb2
mmdb_libclipper_mmdb_la_SOURCES = mmdb/clipper_mmdb.cpp
mmdb_includedir=$(pkgincludedir)/mmdb
mmdb_include_HEADERS = mmdb/clipper_mmdb.h
endif

if BUILD_MINIMOL
lib_LTLIBRARIES += minimol/libclipper-minimol.la
minimol_libclipper_minimol_la_LDFLAGS = @VERSION_INFO@
# minimol uses libmmdb directly (although <mmdb2/mmdb_manager.h> is included
# indirectly, via clipper/clipper-mmdb.h)
minimol_libclipper_minimol_la_LIBADD = core/libclipper-core.la \
				       mmdb/libclipper-mmdb.la -lmmdb2
minimol_libclipper_minimol_la_SOURCES = \
 minimol/minimol.cpp    minimol/container_minimol.cpp minimol/minimol_data.cpp \
 minimol/minimol_io.cpp minimol/minimol_seq.cpp       minimol/minimol_utils.cpp
minimol_includedir=$(pkgincludedir)/minimol
minimol_include_HEADERS = \
 minimol/container_minimol.h minimol/minimol_data.h  minimol/minimol_seq.h   \
 minimol/minimol.h           minimol/minimol_io.h    minimol/minimol_utils.h
endif

if BUILD_PHS
lib_LTLIBRARIES += phs/libclipper-phs.la
phs_libclipper_phs_la_LDFLAGS = @VERSION_INFO@
phs_libclipper_phs_la_LIBADD = core/libclipper-core.la
phs_libclipper_phs_la_SOURCES = phs/phs_io.cpp
phs_includedir=$(pkgincludedir)/phs
phs_include_HEADERS = phs/phs_io.h
endif

if BUILD_CIF
lib_LTLIBRARIES += cif/libclipper-cif.la
cif_libclipper_cif_la_LDFLAGS = @VERSION_INFO@
cif_libclipper_cif_la_LIBADD = core/libclipper-core.la -lmmdb2
cif_libclipper_cif_la_SOURCES = cif/cif_data_io.cpp
cif_includedir=$(pkgincludedir)/cif
cif_include_HEADERS = cif/cif_data_io.h
endif

if BUILD_CCP4
lib_LTLIBRARIES += ccp4/libclipper-ccp4.la
ccp4_libclipper_ccp4_la_LDFLAGS = @VERSION_INFO@
ccp4_libclipper_ccp4_la_LIBADD = core/libclipper-core.la -lccp4c
ccp4_libclipper_ccp4_la_SOURCES = \
        ccp4/ccp4_mtz_types.cpp         ccp4/ccp4_mtz_io.cpp \
        ccp4/ccp4_map_io.cpp            ccp4/ccp4_utils.cpp
ccp4_includedir=$(pkgincludedir)/ccp4
ccp4_include_HEADERS = \
        ccp4/ccp4_map_io.h              ccp4/ccp4_mtz_io.h \
        ccp4/ccp4_mtz_types.h           ccp4/ccp4_utils.h
endif

if BUILD_CCTBX
lib_LTLIBRARIES += cctbx/libclipper-cctbx.la
cctbx_libclipper_cctbx_la_LDFLAGS = @VERSION_INFO@
cctbx_libclipper_cctbx_la_LIBADD = core/libclipper-core.la -lsgtbx -luctbx
cctbx_libclipper_cctbx_la_SOURCES = cctbx/clipper_cctbx.cpp
cctbx_includedir=$(pkgincludedir)/cctbx
cctbx_include_HEADERS = cctbx/clipper_cctbx.h
endif

if BUILD_CNS
lib_LTLIBRARIES += cns/libclipper-cns.la
cns_libclipper_cns_la_LDFLAGS = @VERSION_INFO@
cns_libclipper_cns_la_LIBADD = core/libclipper-core.la
cns_libclipper_cns_la_SOURCES = cns/cns_hkl_io.cpp cns/cns_map_io.cpp
cns_includedir=$(pkgincludedir)/cns
cns_include_HEADERS = cns/cns_hkl_io.h cns/cns_map_io.h
endif

if BUILD_FORTRAN
lib_LTLIBRARIES += fortran/libclipper-fortran.la
fortran_libclipper_fortran_la_LDFLAGS = @VERSION_INFO@
fortran_libclipper_fortran_la_LIBADD = core/libclipper-core.la \
				       contrib/libclipper-contrib.la
fortran_libclipper_fortran_la_SOURCES = fortran/f_interface.cpp
endif
