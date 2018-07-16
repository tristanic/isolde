

m4_define([_AM_PATH_CCTBX_EXTRA],
[

AC_ARG_VAR(BOOST, [boost top dir -optional])

if test "x$ac_cv_env_BOOST_set" != xset; then
  if test "x$cctbx_prefix" != x; then
    BOOST="$cctbx_prefix/../boost"
  fi
fi

ac_CCTBX_CXXFLAGS="$ac_CCTBX_CXXFLAGS -I$BOOST"
#extend for systems that need it
case "$host_os" in
  *osf* | *64* | *irix* )
    ac_CCTBX_CXXFLAGS="$ac_CCTBX_CXXFLAGS -I$BOOST/boost/compatibility/ccp_c_headers"
esac

ac_CCTBX_LDOPTS="$ac_CCTBX_LDOPTS -lboost_python"

])
 
# AM_PATH_CCTBX([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
AC_DEFUN([AM_PATH_CCTBX],
[
AC_PROVIDE([AM_PATH_CCTBX])

AC_ARG_WITH(cctbx,
  AC_HELP_STRING( [--with-cctbx=PFX], [use cctbx package (default is NO) and set prefix] ),
  [
    test "$withval" = no || with_cctbx=yes 
    test "$withval" = yes || cctbx_prefix="$withval" ],
  [ with_cctbx="$enable_cctbx" ] ) #dnl default is no for now

if test "x${with_cctbx}" = xyes ; then  
AS_IF([test "x$CCTBX_LIBS" != x && test "x$CCTBX_CXXFLAGS" != x ],
[
  have_cctbx=yes
],
[
saved_LIBS="$LIBS"
saved_CXXFLAGS="$CXXFLAGS"
CCTBX_LIBS=""
CCTBX_CXXFLAGS=""

if test "x$cctbx_prefix" != x; then
ac_cctbx_dirs='
.
lib
include
build/cctbx/lib'
for ac_dir in $ac_cctbx_dirs; do
  if test -r "$cctbx_prefix/$ac_dir/cctbx/miller.h"; then
    ac_CCTBX_CXXFLAGS="-I$cctbx_prefix/$ac_dir"
    break
    fi
  done
for ac_dir in $ac_cctbx_dirs; do
  for ac_extension in a so sl dylib; do
  if test -r "$cctbx_prefix/$ac_dir/libsgtbx.$ac_extension"; then
    ac_CCTBX_LDOPTS="-L$cctbx_prefix/$ac_dir -lsgtbx -luctbx"
    break 2
    fi
  done
  done
else
 ac_CCTBX_CXXFLAGS=""
 ac_CCTBX_LDOPTS="-lsgtbx -luctbx"
fi

_AM_PATH_CCTBX_EXTRA

AC_MSG_CHECKING([for CCTBX and BOOST])

LIBS="$ac_CCTBX_LDOPTS $saved_LIBS"
CXXFLAGS="$ac_CCTBX_CXXFLAGS $saved_CXXFLAGS"
#
# AC_TRY_LINK uses the c compiler (set by AC_LANG), so we will
# temporarily reassign $CC to the c++ compiler.
#
AC_LANG_PUSH(C++)
AC_TRY_LINK([#include "cctbx/miller.h"] ,[  cctbx::Miller::Index a;  ], have_cctbx=yes, have_cctbx=no)
AC_LANG_POP(C++)  # the language we have just quit
AC_MSG_RESULT($have_cctbx)

 LIBS="$saved_LIBS"
 CXXFLAGS="$saved_CXXFLAGS"
]) # user override

AS_IF([test x$have_cctbx = xyes],
 [
   test "x$CCTBX_CXXFLAGS" = x && CCTBX_CXXFLAGS="$ac_CCTBX_CXXFLAGS"
   test "x$CCTBX_LIBS" = x && CCTBX_LIBS="$ac_CCTBX_LDOPTS"
   ifelse([$1], , :, [$1]) ],
 [
   ifelse([$2], , :, [$2]) ]
)
fi #dnl --with-cctbx

AC_SUBST(CCTBX_CXXFLAGS)
AC_SUBST(CCTBX_LIBS)
])
