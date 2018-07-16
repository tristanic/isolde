
AC_INIT(clipper, 2.1)

AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_MACRO_DIR(m4)
AC_CONFIG_SRCDIR(clipper/clipper.h)
AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects silent-rules -Wall])

AC_PROG_CC
AC_PROG_CXX
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_INSTALL

LT_INIT([disable-shared])

AM_MAINTAINER_MODE dnl disable (by default) maintainer mode

# optional libraries built by default
AC_ARG_ENABLE(contrib, AS_HELP_STRING([--disable-contrib],
                                      [disable contrib library]))
AC_ARG_ENABLE(phs, AS_HELP_STRING([--disable-phs],
                                  [disable PHASEs file interface library]))
AC_ARG_ENABLE(cns, AS_HELP_STRING([--disable-cns],
                                  [enable cns-hkl-interface library]))

# optional libraries not built by default
AC_ARG_ENABLE(mmdb, AS_HELP_STRING([--enable-mmdb],
                     [enable mmdb-interface library (requires mmdb library)]))
AC_ARG_ENABLE(minimol, AS_HELP_STRING([--enable-minimol],
                  [enable minimol library (requires mmdb library)]))
AC_ARG_ENABLE(cif, AS_HELP_STRING([--enable-cif],
                     [enable cif-interface library (requires mmdb library)]))
AC_ARG_ENABLE(ccp4, AS_HELP_STRING([--enable-ccp4],
                     [enable ccp4-interface library (requires ccp4c library)]))
AC_ARG_ENABLE(cctbx, AS_HELP_STRING([--enable-cctbx],
                   [enable cctbx-interface library (requires cctbx package)]))
AC_ARG_ENABLE(fortran, AS_HELP_STRING([--enable-fortran],
                            [enable fortran interface routines]))

AM_CONDITIONAL([BUILD_CONTRIB], [test "x$enable_contrib" != xno])
AM_CONDITIONAL([BUILD_PHS],     [test "x$enable_phs" != xno])
AM_CONDITIONAL([BUILD_CNS],     [test "x$enable_cns" != xno])
AM_CONDITIONAL([BUILD_MMDB],    [test "x$enable_mmdb" = xyes])
AM_CONDITIONAL([BUILD_MMDBOLD], [test "x$enable_mmdbold" = xyes])
AM_CONDITIONAL([BUILD_MINIMOL], [test "x$enable_minimol" = xyes])
AM_CONDITIONAL([BUILD_CIF],     [test "x$enable_cif" = xyes])
AM_CONDITIONAL([BUILD_CCP4],    [test "x$enable_ccp4" = xyes])
AM_CONDITIONAL([BUILD_CCTBX],   [test "x$enable_cctbx" = xyes])
AM_CONDITIONAL([BUILD_FORTRAN], [test "x$enable_fortran" = xyes])

if test "x$enable_fortran" = "xyes"; then
    AC_PROG_F77
    AC_F77_WRAPPERS
fi

AC_SEARCH_LIBS(cos, m, , AC_MSG_ERROR([math library not found.]))
SINGLE_FFTW2

case $host_os in
  cygwin* | mingw* | pw32* | cegcc*)
    AC_CHECK_HEADER(windows.h, :, AC_MSG_ERROR(Windows headers not found))
    ;;
  *)
    AX_PTHREAD
    ;;
esac

CLIPPER_LIBS="-lclipper-core"
test "x$enable_contrib" != xno && CLIPPER_LIBS="-lclipper-contrib $CLIPPER_LIBS"
test "x$enable_phs" != xno     && CLIPPER_LIBS="-lclipper-phs $CLIPPER_LIBS"
test "x$enable_cns" != xno     && CLIPPER_LIBS="-lclipper-cns $CLIPPER_LIBS"
test "x$enable_mmdb" = xyes    && CLIPPER_LIBS="-lclipper-mmdb $CLIPPER_LIBS"
test "x$enable_minimol" = xyes && CLIPPER_LIBS="-lclipper-minimol $CLIPPER_LIBS"
test "x$enable_cif" = xyes     && CLIPPER_LIBS="-lclipper-cif $CLIPPER_LIBS"
test "x$enable_ccp4" = xyes    && CLIPPER_LIBS="-lclipper-ccp4 $CLIPPER_LIBS"
#test "x$enable_cctbx" = xyes   && CLIPPER_LIBS="-lclipper-cctbx $CLIPPER_LIBS"
test "x$enable_fortran" = xyes && CLIPPER_LIBS="-lclipper-fortran $CLIPPER_LIBS"

if test "x$enable_mmdb" = "xyes" -o "x$enable_mmdbold" = "xyes" -o \
        "x$enable_minimol" = "xyes" -o "x$enable_cif" = "xyes"; then
    AC_LANG_PUSH(C++)
    AC_CHECK_LIB(mmdb2, MMDB_CCIF_Init, :,
                 AC_MSG_ERROR(mmdb2 library not found))
    AC_CHECK_HEADER(mmdb2/mmdb_manager.h, :,
                    AC_MSG_ERROR(mmdb2 headers not found))
    AC_LANG_POP(C++)
    # CLIPPER_REQUIRES is only for clipper.pc
    CLIPPER_REQUIRES="mmdb2"
fi

if test "x$enable_ccp4" = "xyes"; then
    AC_CHECK_LIB(ccp4c, ccp4_banner, :, AC_MSG_ERROR(ccp4c library not found))
    AC_CHECK_HEADER(ccp4/ccp4_general.h, :,
                                       AC_MSG_ERROR(libccp4 headers not found))
    CLIPPER_REQUIRES="libccp4c $CLIPPER_REQUIRES"
fi

if test "x$enable_cctbx" = "xyes"; then
    AC_LANG_PUSH(C++)
    AC_CHECK_HEADER(cctbx/miller.h, :, AC_MSG_ERROR(cctbx headers not found))
    AC_MSG_CHECKING(for cctbx)
    AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <cctbx/miller.h>],
                                    [cctbx::Miller::Index a])],
                   AC_MSG_RESULT(yes),
                   [AC_MSG_RESULT(no)
                    AC_MSG_ERROR(cctbx library not found)])
    AC_LANG_POP(C++)
fi

VERSION_INFO="-version-info 2:1"
AC_SUBST(VERSION_INFO)
AC_SUBST(CLIPPER_LIBS)
AC_SUBST(CLIPPER_REQUIRES)

AC_OUTPUT([Makefile
           clipper/Makefile
           clipper/core/Makefile
           clipper/contrib/Makefile
           examples/Makefile
           clipper.pc])


echo "       "
echo "       " Configuration Summary
echo "" ------------------------------
echo core: "        "  yes
echo contrib: "     "  ${enable_contrib:-yes}
echo phs: "         "  ${enable_phs:-yes}
echo cns: "         "  ${enable_cns:-yes}
echo mmdb: "        "  ${enable_mmdb:-no}
echo minimol: "     "  ${enable_minimol:-no}
echo cif: "         "  ${enable_cif:-no}
echo ccp4: "        "  ${enable_ccp4:-no}
echo cctbx: "       "  ${enable_cctbx:-no}
echo fortran: "     "  ${enable_fortran:-no}
echo " "

