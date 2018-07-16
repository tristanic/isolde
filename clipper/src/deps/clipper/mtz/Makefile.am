#   Copyright 2002 
#   For licence conditions, see the file COPYING in the distribution.

libclipper_mtz_la_SOURCES = mtz_types.cpp mtz_io.cpp map_io.cpp

libclipper_mtz_a_SOURCES = $(libclipper_mtz_la_SOURCES)

library_include_HEADERS = map_io.h     mtz_io.h     mtz_types.h

library_includedir = $(pkgincludedir)/mtz

#libtool
lib_LTLIBRARIES = libclipper-mtz.la

#non libtool
#noinst_LIBRARIES = libclipper-mtz.a

libclipper_mtz_la_LDFLAGS = -version-info 1:0

INCLUDES = $(MCCP4_CXXFLAGS)

libclipper_mtz_a_AR = $(AR) $(AR_FLAGS)
