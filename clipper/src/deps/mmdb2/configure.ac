
AC_INIT(mmdb2, 2.0.12) # keep in sync with mmdb2/mmdb_defs.h

AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_MACRO_DIR(m4)
AC_CONFIG_SRCDIR(mmdb2/mmdb_math_align.cpp)
AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects silent-rules -Wall])
AM_PROG_AR
LT_INIT([disable-shared win32-dll])
AC_PROG_CXX
# AC_PROG_CXX doesn't fail if compiler is not found, test it explicitely
AC_LANG_PUSH(C++)
AC_LINK_IFELSE([AC_LANG_PROGRAM([], [])], [],
               [AC_MSG_ERROR([C++ compiler not found.])])
AC_LANG_POP(C++)

AM_MAINTAINER_MODE dnl disable (by default) maintainer mode

AC_SEARCH_LIBS([cos], [m], [], [AC_MSG_ERROR([need math library])])
AC_OUTPUT(Makefile mmdb2.pc)
