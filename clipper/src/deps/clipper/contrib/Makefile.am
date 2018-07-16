#   Copyright 2002 
#   For licence conditions, see the file COPYING in the distribution.

DIST_EXTRA = clipper_contrib.dox

lib_LTLIBRARIES = libclipper-contrib.la
libclipper_contrib_la_LDFLAGS = @VERSION_INFO@
libclipper_contrib_la_LIBADD = ../core/libclipper-core.la

libclipper_contrib_la_SOURCES = \
convolution_search.cpp          sfcalc.cpp \
edcalc.cpp                      sfcalc_obs.cpp \
fffear.cpp                      sfscale.cpp \
function_object_bases.cpp       sfweight.cpp \
mapfilter.cpp                   skeleton.cpp \
originmatch.cpp                 test_contrib.cpp

library_includedir = $(pkgincludedir)/contrib
library_include_HEADERS = \
convolution_search.h    mapfilter.h             sfscale.h \
edcalc.h                originmatch.h           sfweight.h \
fffear.h                sfcalc.h                skeleton.h \
function_object_bases.h sfcalc_obs.h            test_contrib.h

