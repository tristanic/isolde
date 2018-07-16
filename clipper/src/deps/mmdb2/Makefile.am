
ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = mmdb2/libmmdb2.la

pkginclude_HEADERS = \
mmdb2/hybrid_36.h		mmdb2/mmdb_io_file.h		mmdb2/mmdb_math_graph.h		mmdb2/mmdb_seqsuperpose.h \
mmdb2/mmdb_atom.h		mmdb2/mmdb_io_stream.h		mmdb2/mmdb_math_linalg.h	mmdb2/mmdb_symop.h        \
mmdb2/mmdb_bondmngr.h		mmdb2/mmdb_machine_.h		mmdb2/mmdb_math_rand.h		mmdb2/mmdb_tables.h       \
mmdb2/mmdb_chain.h		mmdb2/mmdb_manager.h		mmdb2/mmdb_mattype.h		mmdb2/mmdb_title.h        \
mmdb2/mmdb_cifdefs.h		mmdb2/mmdb_mask.h		mmdb2/mmdb_mmcif_.h		mmdb2/mmdb_uddata.h       \
mmdb2/mmdb_coormngr.h		mmdb2/mmdb_math_.h		mmdb2/mmdb_model.h		mmdb2/mmdb_utils.h        \
mmdb2/mmdb_cryst.h		mmdb2/mmdb_math_align.h		mmdb2/mmdb_root.h		mmdb2/mmdb_xml_.h         \
mmdb2/mmdb_defs.h		mmdb2/mmdb_math_bfgsmin.h	mmdb2/mmdb_rwbrook.h		mmdb2/mmdb_ficif.h	  \
mmdb2/mmdb_math_fft.h		mmdb2/mmdb_selmngr.h


mmdb2_libmmdb2_la_SOURCES = \
mmdb2/hybrid_36.cpp		mmdb2/mmdb_io_stream.cpp	mmdb2/mmdb_math_linalg.cpp	mmdb2/mmdb_symop.cpp      \
mmdb2/mmdb_atom.cpp		mmdb2/mmdb_machine_.cpp		mmdb2/mmdb_math_rand.cpp	mmdb2/mmdb_tables.cpp     \
mmdb2/mmdb_bondmngr.cpp		mmdb2/mmdb_manager.cpp		mmdb2/mmdb_mattype.cpp		mmdb2/mmdb_title.cpp      \
mmdb2/mmdb_chain.cpp		mmdb2/mmdb_mask.cpp		mmdb2/mmdb_mmcif_.cpp		mmdb2/mmdb_uddata.cpp     \
mmdb2/mmdb_cifdefs.cpp		mmdb2/mmdb_math_.cpp		mmdb2/mmdb_model.cpp		mmdb2/mmdb_utils.cpp      \
mmdb2/mmdb_coormngr.cpp		mmdb2/mmdb_math_align.cpp	mmdb2/mmdb_root.cpp		mmdb2/mmdb_xml_.cpp       \
mmdb2/mmdb_cryst.cpp		mmdb2/mmdb_math_bfgsmin.cpp	mmdb2/mmdb_rwbrook.cpp		mmdb2/mmdb_ficif.cpp	  \
mmdb2/mmdb_math_fft.cpp		mmdb2/mmdb_selmngr.cpp		mmdb2/mmdb_io_file.cpp		mmdb2/mmdb_math_graph.cpp \
mmdb2/mmdb_seqsuperpose.cpp

mmdb2_libmmdb2_la_LDFLAGS = -no-undefined

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = mmdb2.pc


