
AC_INIT(libccp4, 7.0.0) # keep version in sync with ccp4/ccp4_program.h

AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_MACRO_DIR(m4)
AC_CONFIG_SRCDIR(ccp4/cmtzlib.c)
AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects silent-rules -Wall])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
LT_INIT([disable-shared win32-dll])

AC_PROG_CC
dnl fortran/rwbrook.f provides interface to fortran mmdb interface
dnl (mmdb/mmdb_rwbrook.h/cpp). CXX is used only to test linking with mmdb.
AC_PROG_CXX
AC_PROG_INSTALL

AM_MAINTAINER_MODE dnl disable (by default) maintainer mode

AC_ARG_ENABLE(fortran, AS_HELP_STRING([--disable-fortran],
                                      [do not compile fortran API]))

AM_CONDITIONAL([FORTRAN_LIB], [ test "x$enable_fortran" != xno ])

if test "x$enable_fortran" != "xno"; then
  # sorting_main.f (mix of F77 and F90) cannot be compiled with g77 
  AC_PROG_F77([_AC_F90_FC _AC_F95_FC _AC_F77_FC])
  # AC_PROG_F77 doesn't fail if compiler is not found, test it explicitely
  AC_LANG_PUSH(Fortran 77)
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], [      write(*,*) 'hello' ])], [],
                 [AC_MSG_ERROR([Fortran compiler not found.])])
  AC_LANG_POP(Fortran 77)
  AC_PATH_PROG(M4,m4)
  if test x"$M4" = "x"; then
      AC_MSG_ERROR([m4 program not found.])
  fi
  # TODO: it'd be better to check for features rather than compilers,
  # but I'm afraid of changing this and unix.m4.
  cc_basename=`$ECHO "X$compiler_F77" | $Xsed -e 's%.*/%%' -e "s%^$host_alias-%%"`
  case $host_os in
    linux*|mingw* )
      case $cc_basename,$host_cpu in
        ifort*)
          M4_ARG=${M4_ARG:-"-D_ifc=8"}
          AC_DEFINE([IFC],[8],[using intel fortran compiler])
          ;;
        ifc*)
          M4_ARG=${M4_ARG:-"-D_efc=1"}
          AC_DEFINE([IFC],[1],[using intel fortran compiler])
          ;;
        g77*)
          M4_ARG=${M4_ARG:-"-D_g77=1"}
          AC_DEFINE([G77],[1],[using g77 compiler])
          ;;
        gfortran*)
          M4_ARG=${M4_ARG-"-D_gfort=1"}
          AC_DEFINE([GFORTRAN],[1],[using gfortran compiler])
          ;;
        xlf*)
          M4_ARG=${M4_ARG-"-D_xlf=1"} ;;
        pgf*) ;;
        esac
    ;;
    darwin* )
      case $cc_basename,$host_cpu in
        ifort*)
          M4_ARG=${M4_ARG:-"-D_ifc=8"}
          AC_DEFINE([IFC],[8],[using intel fortran compiler])
          ;;
        ifc*)
          M4_ARG=${M4_ARG:-"-D_ifc=1"}
          AC_DEFINE([IFC],[1],[using intel fortran compiler])
          ;;
        g77*)
          M4_ARG=${M4_ARG:-"-D_g77=1"}
          AC_DEFINE([G77],[1],[using g77 compiler])
          ;;
        gfortran*)
          M4_ARG=${M4_ARG-"-D_gfort=1"}
          AC_DEFINE([GFORTRAN],[1],[using gfortran compiler])
          ;;
        xlf*)
          M4_ARG=${M4_ARG-"-D_xlf=1"} ;;
        esac
    ;;
    irix* )
      case $cc_basename,$host_cpu in
        f77* )
          M4_ARG=${M4_ARG:-"-D_sgi=1"}
          ;;
        f90*)
          M4_ARG=${M4_ARG:-"-D_sgi=1 -D_f90=1"}
          ;;
        esac
  esac

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [return _WIN32;])],
                    [M4_ARG="$M4_ARG -D_windows=1"], [])

  AC_MSG_NOTICE([fortran/unix.m4 will be processed with m4 $M4_ARG])

  AC_SEARCH_LIBS(cos, m, , AC_MSG_ERROR([math library not found.]))
  AC_LANG_PUSH([C++])
  AC_CHECK_LIB(mmdb2, MMDB_CCIF_Init, :,
               AC_MSG_WARN([mmdb2 library not found.]))
  AC_LANG_POP([C++])

fi #dnl fortran API elements

AM_CONDITIONAL([GNU_FORTRAN], [test "x$G77" = "xyes"])

# TODO: prefix or pkgdatadir can be changed after configure is run,
# setting PACKAGE_ROOT here can fail, it would be better to set:
# ccp4_libccp4c_la_CPPFLAGS= -DDATA_DIR=\"$(pkgdatadir)\"
PACKAGE_ROOT=$prefix
test "x$PACKAGE_ROOT" = xNONE && PACKAGE_ROOT=$ac_default_prefix
AC_DEFINE_UNQUOTED([PACKAGE_ROOT],["$PACKAGE_ROOT"],[package install directory])

AC_SUBST(M4_ARG)

AC_OUTPUT(Makefile libccp4c.pc libccp4f.pc)

