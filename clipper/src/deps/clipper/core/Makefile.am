# Copyright 2003
# For licence conditions, see the file COPYING in the distribution

lib_LTLIBRARIES = libclipper-core.la
libclipper_core_la_LDFLAGS = @VERSION_INFO@
libclipper_core_la_LIBADD = $(FFTW2_LIBS) $(PTHREAD_LIBS)
libclipper_core_la_CPPFLAGS = $(PTHREAD_CFLAGS)

libclipper_core_la_SOURCES = \
atomsf.cpp              coords.cpp              nxmap_operator.cpp \
cell.cpp                derivs.cpp              ramachandran.cpp \
clipper_instance.cpp    fftmap.cpp              resol_basisfn.cpp \
clipper_memory.cpp      fftmap_sparse.cpp       resol_fn.cpp \
clipper_message.cpp     hkl_compute.cpp         resol_targetfn.cpp \
clipper_stats.cpp       hkl_data.cpp            rotation.cpp \
clipper_test.cpp        hkl_datatypes.cpp       spacegroup.cpp \
clipper_types.cpp       hkl_info.cpp            spacegroup_data.cpp \
clipper_util.cpp        hkl_lookup.cpp          symop.cpp \
container.cpp           hkl_operators.cpp       container_hkl.cpp \
map_interp.cpp          container_map.cpp       map_utils.cpp \
xmap.cpp                container_types.cpp     nxmap.cpp \
clipper_thread.cpp      test_core.cpp           test_data.cpp

library_includedir=$(pkgincludedir)/core
library_include_HEADERS = \
atomsf.h                container_types.h       nxmap_operator.h \
cell.h                  coords.h                ramachandran.h \
clipper_instance.h      derivs.h                resol_basisfn.h \
clipper_memory.h        fftmap.h                resol_fn.h \
clipper_message.h       fftmap_sparse.h         resol_targetfn.h \
clipper_precision.h     hkl_compute.h           rotation.h \
clipper_stats.h         hkl_data.h              spacegroup.h \
clipper_sysdep.h        hkl_datatypes.h         spacegroup_data.h \
clipper_test.h          hkl_info.h              symop.h \
clipper_types.h         hkl_lookup.h            clipper_util.h \
hkl_operators.h         container.h             map_interp.h \
xmap.h                  container_hkl.h         map_utils.h \
container_map.h         nxmap.h                 clipper_thread.h \
test_core.h             test_data.h

