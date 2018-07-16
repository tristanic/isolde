/*
     library_f.c: Fortran API to library.c
     Copyright (C) 2001  CCLRC, Charles Ballard

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the
     license to address the requirements of UK law.

     You should have received a copy of the modified GNU Lesser General
     Public License along with this library.  If not, copies may be
     downloaded from http://www.ccp4.ac.uk/ccp4license.php

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
*/

/** @file library_f.c
 *  FORTRAN API for library.c.
 *  Charles Ballard
 */

/* FORTRAN API for library.c                                                */
/* also include missing routines and wrappers for C commands                */
/*                                                                          */
/* revisions:                                                               */
/*           (4/5/01) C.Ballard                                             */
/*                    respond to first steps to make library.c              */
/*                    more "C" like.                                        */
/*           (21/8/01) C.Ballard                                            */
/*                     error catching from library.c                        */
/*                                                                          */
/* % Copyright Daresbury Laboratory 1992--2001                              */
/* % This is a CCP4 `part (i)' file for the purposes of copyright.          */
/* % See the CCP4 distribution conditions for explanation.                  */
/*                                                                          */
/* % \documentstyle[a4wide,times,noweb,makeidx]{article}                    */
/*                                                                          */
/* % \newcommand{\ac}[1]{{\rm\normalshape\sc #1}}   % acronym               */
/*                                                                          */
/* \documentclass{article}                                                  */
/* \usepackage{a4wide,times,noweb,makeidx}                                  */
/* \newcommand{\ac}[1]{\textsc{#1}}   % acronym                             */
/* \newcommand{\meta}[1]{\mbox{$\langle$\sl #1\/$\rangle$}}                 */
/* \newcommand{\ft}{\idx{Fortran}}                                          */
/* \newcommand{\idx}[1]{#1\index{#1}}                                       */
/* \newcommand{\fixme}[1]{\index{Fixme!}[{\bf Fixme!:} #1\@.]}              */
/*                                                                          */
/* \title{FORTRAN wrapper library routines}                                 */
/* \date{$ $Date$ $}                                  */
/* \author{This version: Martyn Winn, Charles Ballard @ Daresbury}          */
/*                                                                          */
/* \makeindex                                                               */
/*                                                                          */
/* \noweboptions{longchunks,smallcode}                                      */
/*                                                                          */
/* \begin{document}                                                         */
/*                                                                          */
/* \maketitle                                                               */
/*                                                                          */
/* \noindent                                                                */
/* This file contains the wrappers for calling library.c from FORTRAN and   */
/* some "missing" routines.                                                 */
/* \bigskip                                                                 */
/*                                                                          */
/* \tableofcontents                                                         */
/*                                                                          */
/*                                                                          */
/* \section{Summary}                                                        */
/*                                                                          */
/* The following routines are defined:                                      */
/* \bigskip                                                                 */
/*                                                                          */
/* \noindent                                                                */
/* \begin{tabular}{ll}                                                      */
/*                                                                          */
/* Routine and  arguments &                      Purpose \\                 */
/* \hline                                                                   */
/* [[ustenv(string, result)]]         & set an environment variable  \\     */
/* [[ccpal1(routine,number,type,length)]] & allocate arrays of type and     */
/*                                    & call routine \\                     */
/* \end{tabular}                                                            */
/*                                                                          */
/*                                                                          */
/* \section{Portability and Code}                                           */
/*                                                                          */
/* System dependent names are handled in the FORTRAN_SUBR,                  */
/* FORTRAN_FUN, FORTRAN_CALL macros defined in the header file.             */
/* fpstr is a typedef which masks the intricacies of FORTRAN string         */
/* passing.                                                                 */
/*                                                                          */
/* <*>=                                                                     */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include "../ccp4/ccp4_utils.h"
#include "../ccp4/ccp4_errno.h"
#include "../ccp4/ccp4_fortran.h"

#if defined(GFORTRAN) || defined (G95)
#include <time.h>
#endif

/* rcsid[] = "$Id$" */

/** Creates a null-terminated C string from an input
 * string obtained from a Fortran call. Trailing blanks are
 * removed. If input string is blank then return string "\0".
 * Memory assigned by malloc, so can be freed.
 * @param str1 pointer to string
 * @param str1_len Fortran length of string
 */
char *ccp4_FtoCString(fpstr str1, int str1_len)
{
  char *str2;

  size_t length = ccp4_utils_flength(FTN_STR(str1),str1_len);
  str2 = (char *) ccp4_utils_malloc((length+1)*sizeof(char));
  if(length) strncpy(str2, FTN_STR(str1), length);
  str2[length] = '\0';

  return str2;
}

/** Creates a Fortran string from an input C string for passing back to
 * Fortran call. Characters after null-terminator may be junk, so pad
 * with spaces. If input cstring is NULL, return blank string.
 * @param str1 pointer Fortran to string
 * @param str1_len Fortran length of string
 * @param cstring input C string
 */
void ccp4_CtoFString(fpstr str1, int str1_len, const char *cstring)
{
  int i;

  if (!cstring) {
    for (i = 0; i < str1_len; ++i)
      str1[i] = ' ';
  } else if (str1_len > strlen(cstring)) {
    strcpy(FTN_STR(str1),cstring);
    for (i = strlen(cstring); i < str1_len; ++i)
      str1[i] = ' ';
  } else {
    strncpy(FTN_STR(str1),cstring,str1_len);
  }
}

/* \section{Miscellaneous routines}                                         */
/* \subsection{{\tt subroutine ustenv(\meta{string}, \meta{result})}}       */
/*                                                                          */
/* This sets an environment variable \meta{var} to \meta{val}, where the    */
/* argument \meta{string}[[==']]\meta{var}[['//'='//']]\meta{val}[[']].     */
/* This is for use by the `\idx{logical name}' mechanism for specifying     */
/* file connexions.  Note that a \idx{VMS} varsion is supplied in {\tt      */
/*   vms.for} and that there is no standard way of setting and              */
/* environment variable.  In a minimal \ac{posix} system it might be        */
/* necessary to twiddle the environment strings explicitly.                 */
/* Upon exit result contains [[0]] on Success, [[-1]] on Failure.           */
/*                                                                          */
/*                                                                          */
/* <miscellaneous routines>=                                                */
/* <ustenv code>=                                                           */
#if ! defined (VMS)
FORTRAN_SUBR ( USTENV, ustenv,
         (fpstr str, int *result, int str_len),
         (fpstr str, int *result),
         (fpstr str, int str_len, int *result))
{
  char *temp_name;

  temp_name = ccp4_FtoCString(FTN_STR(str), FTN_LEN(str));

  if ((*result = ccp4_utils_setenv (temp_name)) != 0)
    ccp4_fatal("USTENV/CCP4_SETENV: Memory allocation failure");
  free(temp_name);
}
#endif

#if ! defined (_MSC_VER)
FORTRAN_SUBR ( USTIME, ustime,
         (int *isec),
         (int *isec),
         (int *isec))
{
  *isec = time(NULL);
}
#endif

/* \section{Miscellaneous routines}                                         */
/* \subsection{{\tt outbuf()}}                                              */
/*                                                                          */
/* This sets stdout to line buffering (error not fatal)                     */
/*                                                                          */
/* <miscellaneous routines>=                                                */
/* <outbuf code>=                                                           */
FORTRAN_SUBR ( OUTBUF, outbuf, (), (), ())
{
#if defined (__APPLE__) && defined (_CALL_SYSV)
   char *s = "buffering=disable_preconn";
   int s_len = strlen(s);
   FORTRAN_CALL (SETRTEOPTS,setrteopts,(s,s_len),(s,s_len),(s,s_len));
#endif
  if(ccp4_utils_outbuf())
    ccp4_utils_print("OUTBUF:Can't turn off buffering");
}

/* \subsection{{\tt subroutine cunlink (\meta{filename})}}                  */
/* This unlinks \meta{filename} from the directory.  It's intended for      */
/* use with scratch files, so that they can be hidden when opened but       */
/* still be available as long as they remain connected (see [[CCPOPN]]).    */
/* This functionality doesn't seem to exist in \idx{VMS}\@.  Failure to     */
/* unlink isn't fatal (it's been observed, apparently spuriously).          */
/*                                                                          */
/* <miscellaneous routines>=                                                */
FORTRAN_SUBR ( CUNLINK, cunlink,
      (fpstr filename, int filename_len),
      (fpstr filename),
      (fpstr filename, int filename_len))
{
#ifdef VMS
  return;                       /* can't do it */
#else
  char *temp_name;

  temp_name = ccp4_FtoCString(FTN_STR(filename), FTN_LEN(filename));

  if( unlink(temp_name) )
    ccp4_utils_print("CUNLINK: Can't unlink");
  free(temp_name);
#endif /* VMS */
}

/* \section{Dynamic memory allocation}                                      */
/* It's nice to be able to determine array sizes at run time to avoid       */
/* messy recompilation.  The only way effectively to get dynamic            */
/* allocation in Fortran77 reasonably portably is to do the allocation,     */
/* e.g.\ in C, and invoke the Fortran routine passed as a parameter with    */
/* pointers to the allocated memory which it will treat as arrays.  If we   */
/* want to allow more than one array, it's more tricky.                     */
/*                                                                          */
/* \subsection{{\tt subroutine ccpal1 (\meta{routne}, \meta{n}.             */
/*     \meta{type}, \meta{length})}}                                        */
/* Arranges to call subroutine \meta{routne} with \meta{n} array            */
/* arguments.  Each has a type indicated by \meta{type}$(i)$ and a length   */
/* given by \meta{length}($i$).  \meta{type} is an integer array with       */
/* values 1, 2, 3, 4 inidcating {\tt                                        */
/*   INTEGER}, {\tt REAL}, {\tt DOUBLE PRECISION} and {\tt COMPLEX}         */
/* respectively.                                                            */
/* It's not immediately clear what all the Fortran/C                        */
/* conventions are for passing [[CHARACTER]] arrays, so we'll arrange a     */
/* higher-level interface and have [[types]] here just numeric.  The        */
/* Fortran ([[CCPALC]]) will also do argument validation.  Also the rules   */
/* for passing external routines as arguments aren't clear---assume         */
/* the obvious way.                                                         */
/*                                                                          */
/* There's a \idx{VMS} Fortran version of this, although the code here      */
/* does work fine in VMS\@.                                                 */
/*                                                                          */
/* NB: there's a possibility of a hook here to use memory-mapped files on   */
/* systems with the capability and insufficient VM\@.                       */
/*                                                                          */
/* Under protest, this now allocates zeroed storage for where programs      */
/* make bad assumptions.                                                    */
/*                                                                          */
/* <miscellaneous routines>=                                                */
#ifndef VMS                     /* we'll use the Fortran version in VMS*/
#ifndef _MSC_VER
FORTRAN_SUBR ( CCPAL1, ccpal1,
     (void (* routne) (), int *n, int type[], int length[]),
     (void (* routne) (), int *n, int type[], int length[]),
     (void (* routne) (), int *n, int type[], int length[]))
{
  static int item_sizes[] = {
    (int) sizeof (char),           /* 0: bytes */
    (int) sizeof (short int),      /* 1: (integer) half words */
    (int) sizeof (float),          /* 2: reals/words */
    (int) sizeof (int),            /* 3: `short complex' (pairs of half words).
                                         NB int rather than 2*short since must fit
                                          into fortran integer */
    (int) 2*sizeof (float),        /* 4: complex (pairs of words) */
    (int) sizeof (int),            /* 5: not used */
    (int) sizeof (int)             /* 6: integers */
  };
  int i, size, *leng[13];
  void *pointer[13];

  for (i=0; i<*n; i++) {
    switch (type[i]) {
    case 1:
      size = item_sizes[6]; break; /* integer */
    case 2:
      size = item_sizes[2]; break; /* real */
    case 3:
      size = 2*item_sizes[2]; break; /* double */
    case 4:
      size = 2*item_sizes[2]; break; /* complex */
    case 5:
      size = item_sizes[1]; break; /* bytes (logical or integer *1) */
    }
    pointer[i+1] = calloc ((size_t) length[i], (size_t) size);
    if (pointer[i+1] == NULL) ccp4_fatal ("CCPALC: can't allocate memory");
    leng[i+1] = &(length[i]);   /* convenience */
  }
  switch (*n) {
  case 1:
    (* routne) (leng[1], pointer[1]);
    break;
  case 2:
    (* routne) (leng[1], pointer[1], leng[2], pointer[2]);
    break;
  case 3:
    (* routne) (leng[1], pointer[1], leng[2], pointer[2],
                leng[3], pointer[3]);
    break;
  case 4:
    (* routne) (leng[1], pointer[1], leng[2], pointer[2],
                leng[3], pointer[3], leng[4], pointer[4]);
    break;
  case 5:
    (* routne) (leng[1], pointer[1], leng[2], pointer[2],
                leng[3], pointer[3], leng[4], pointer[4],
                leng[5], pointer[5]);
    break;
  case 6:
    (* routne) (leng[1], pointer[1], leng[2], pointer[2],
                leng[3], pointer[3], leng[4], pointer[4],
                leng[5], pointer[5], leng[6], pointer[6]);
    break;
  case 7:
    (* routne) (leng[1], pointer[1], leng[2], pointer[2],
                leng[3], pointer[3], leng[4], pointer[4],
                leng[5], pointer[5], leng[6], pointer[6],
                leng[7], pointer[7]);
    break;
  case 8:
    (* routne) (leng[1], pointer[1], leng[2], pointer[2],
                leng[3], pointer[3], leng[4], pointer[4],
                leng[5], pointer[5], leng[6], pointer[6],
                leng[7], pointer[7], leng[8], pointer[8]);
    break;
  case 9:
    (* routne) (leng[1], pointer[1], leng[2], pointer[2],
                leng[3], pointer[3], leng[4], pointer[4],
                leng[5], pointer[5], leng[6], pointer[6],
                leng[7], pointer[7], leng[8], pointer[8],
                leng[9], pointer[9]);
    break;
  case 10:
    (* routne) (leng[1], pointer[1], leng[2], pointer[2],
                leng[3], pointer[3], leng[4], pointer[4],
                leng[5], pointer[5], leng[6], pointer[6],
                leng[7], pointer[7], leng[8], pointer[8],
                leng[9], pointer[9], leng[10], pointer[10]);
    break;
  case 11:
    (* routne) (leng[1], pointer[1], leng[2], pointer[2],
                leng[3], pointer[3], leng[4], pointer[4],
                leng[5], pointer[5], leng[6], pointer[6],
                leng[7], pointer[7], leng[8], pointer[8],
                leng[9], pointer[9], leng[10], pointer[10],
                leng[11], pointer[11]);
    break;
  case 12:
    (* routne) (leng[1], pointer[1], leng[2], pointer[2],
                leng[3], pointer[3], leng[4], pointer[4],
                leng[5], pointer[5], leng[6], pointer[6],
                leng[7], pointer[7], leng[8], pointer[8],
                leng[9], pointer[9], leng[10], pointer[10],
                leng[11], pointer[11], leng[12], pointer[12]);
    break;
  }
  for (i=0; i<*n; i++)
    free (pointer[i+1]);
}
#endif /* VMS */
#endif

/* \section{`Magic' numbers}                                                */
/*                                                                          */
/* When, for instance, an $F$ is unobserved in a derivative, we might       */
/* want to give it a special value---a `\idx{magic number}'---possibly in   */
/* addition to a special value of the $\sigma$, like a negative one.        */
/* Using such a number in a calculation (by mistake, through ignoring the   */
/* value of $\sigma$, say) should not allow one to get half-sensible        */
/* results as one might if this number was $-9999$ or some such.  (There    */
/* is non-enforced connexion between the $F$ and its $\sigma$ in the MTZ    */
/* file, although one could think of adding extra columns to the file       */
/* with bit-encoded flags telling whether the $F$ in a given column was     */
/* observed.)                                                               */
/*                                                                          */
/* The obvious tactic with \ac{ieee} arithmetic is to use a \idx{NaN}       */
/* value in such situations.  Things may be set up so that we either get    */
/* an exception on using it in arithmetic or it silently propagates to all  */
/* values using it and its presence is indicated by a NaN in the output.    */
/* On a \idx{VAX} architecture we can't use NaN, but there is the           */
/* possibility of using a                                                   */
/* `reserved operand'\index{reserved operand|see{Rop}}                      */
/* (`\idx{Rop}') value,                                                     */
/* which will cause an exception (by experiment: when used for              */
/* floating-point arithmetic {\em or\/} printed, but not when assigned).    */
/* The \idx{Convex} native mode is similar, except that the Rop may be      */
/* printed (in the form {\tt Rop0x}\meta{fraction part}).                   */
/*                                                                          */
/* On, say, the \idx{IBM 370 architecture}---which we don't currently       */
/* support---anything's a valid floating point number, and the best ploy    */
/* is probably to use the largest representable number as the `magic'       */
/* value.  This would stand a good chance of raising an overflow            */
/* exception if used.  Anyhow, if such bad use of an undefined value is     */
/* made in a program due to insufficient checking by the code, it should    */
/* be spotted on the \ac{ieee} systems and the bug fixed---it's not         */
/* strictly necessary that it should cause a fatal error on all             */
/* architectures.                                                           */
/*                                                                          */
/* We need to provide a means of setting the magic number and checking      */
/* whether a given value is such.  These are architecture-dependent         */
/* bit-level operations, hence their presence in the C code.                */
/*                                                                          */
/* The suite doesn't currently use these routines, but should do soon.      */
/* \subsection{Setting a value: {\tt subroutine qnan(value)}}               */
/*                                                                          */
/* [[qnan]] was originally a \ft{} [[real function]] returning the value    */
/* (and actually done in 2 stages) with a subroutine implementation like    */
/* this called by the \ft{} function to avoid problems under \idx{VMS}      */
/* and native \idx{Convex}.  However, the \idx{f2c} calling convention      */
/* for a function loses in that case since it assumes a [[double]] value    */
/* returned which is cast to [[float]] with a SIGFPE, sigh.                 */
/*                                                                          */
/* <magic numbers>=                                                         */
FORTRAN_SUBR ( QNAN, qnan,
    (union float_uint_uchar *realnum),
    (union float_uint_uchar *realnum),
    (union float_uint_uchar *realnum))
{
  *realnum = ccp4_nan ();
}
/* \subsection{Testing a value: {\tt int qisnan(\meta{real})}}              */
/*                                                                          */
/* We want a \ft{} logical function [[qisnan]] to test whether its argument */
/* is a \idx{NaN} or \idx{Rop}.  We have to do this by writing a C          */
/* [[int]]-valued procedure and testing the returned value in the \ft{}     */
/* so that we don't have to assume how it represents logical values.  The   */
/* {\tt diskio}\index{diskio} library module provides the                   */
/* trivial interface [[QISNAN]].                                            */
/*                                                                          */
/* <magic numbers>=                                                         */
FORTRAN_FUN (int, QISNAN, qisnan,
	     (union float_uint_uchar *realnum),
	     (union float_uint_uchar *realnum),
	     (union float_uint_uchar *realnum))
{
  return (_BTOLV(ccp4_utils_isnan (realnum)));
}

/* \subsection{Absent data test for {\tt mtzlib}: {\tt subroutine           */
/*     ccpbml (\meta{ncols}, \meta{cols})}}                                 */
/* In {\tt mtzlib} there's a fudge for \idx{BIOMOL}-convention absence      */
/* flags, which are re-written to zeroes.  To do the real number            */
/* comparison, though, it's necessary to do a [[qnan]]-type test first.     */
/* We don't want to call [[qnan]] (which calls [[cisnan]]) on every         */
/* number in the data file, so the tests are amortised in this routine      */
/* which deals with a whole array \meta{cols} of length \meta{ncols}.       */
/*                                                                          */
/* <magic numbers>=                                                         */
FORTRAN_SUBR ( CCPBML, ccpbml,
    (int *ncols, union float_uint_uchar cols[]),
    (int *ncols, union float_uint_uchar cols[]),
    (int *ncols, union float_uint_uchar cols[]))
{
  ccp4_utils_bml (*ncols, cols) ;
}

/* \subsection{Updating MTZ column ranges: {\tt subroutine ccpwrg           */
/*     (\meta{ncols}, \meta{rcols}, \meta{wmin}, \meta{wmax})}}             */
/* This is a similar fudge to [[ccpbml]] to avoid [[QISNAN]] calls in       */
/* updating the MTZ column ranges in {\tt mtzlib}.  Note that [[wminmax]]   */
/* actually indexes a 3-D Fortran array with the first                      */
/* dimension range of 2, indicating minimum and maximum values respectively. */
/*                                                                          */
/* <magic numbers>=                                                         */
FORTRAN_SUBR ( CCPWRG, ccpwrg,
    (int *ncols, union float_uint_uchar cols[], float wminmax[]),
    (int *ncols, union float_uint_uchar cols[], float wminmax[]),
    (int *ncols, union float_uint_uchar cols[], float wminmax[]))
{
  ccp4_utils_wrg (*ncols, cols, wminmax) ;
}

/* \subsection{Routines for Data Harvesting: {\tt subroutine hgetlimits}}    */
/* Returns largest int and largest float as defined in <limits.h> and       */
/* <float.h>                                                                 */
FORTRAN_SUBR ( HGETLIMITS, hgetlimits,
    (int *IValueNotDet, float *ValueNotDet),
    (int *IValueNotDet, float *ValueNotDet),
    (int *IValueNotDet, float *ValueNotDet))
{
  ccp4_utils_hgetlimits (IValueNotDet, ValueNotDet);
}

/* Wrap-around for mkdir function. Returns 0 if successful, 1 if directory  */
/* already exists, and -1 if other error.                                   */
FORTRAN_SUBR ( CMKDIR, cmkdir,
    (const fpstr path, const fpstr cmode, int *result, int path_len, int cmode_len),
    (const fpstr path, const fpstr cmode, int *result),
    (const fpstr path, int path_len, const fpstr cmode, int cmode_len, int *result))
{
  char *temp_path, *temp_cmode;

  temp_path = ccp4_FtoCString(FTN_STR(path), FTN_LEN(path));
  temp_cmode = ccp4_FtoCString(FTN_STR(cmode), FTN_LEN(cmode));

  *result = ccp4_utils_mkdir (temp_path, temp_cmode);
  free(temp_path);
  free(temp_cmode);
}

/* Wrap-around for mkdir function. Returns 0 if successful, 1 if directory     */
/* already exists, and -1 if other error.                                      */
FORTRAN_SUBR ( CCHMOD, cchmod,
    (const fpstr path, const fpstr cmode, int *result, int path_len, int cmode_len),
    (const fpstr path, const fpstr cmode, int *result),
    (const fpstr path, int path_len, const fpstr cmode, int cmode_len, int *result))
{
  char *temp_path, *temp_cmode;

  temp_path = ccp4_FtoCString(FTN_STR(path), FTN_LEN(path));
  temp_cmode = ccp4_FtoCString(FTN_STR(cmode), FTN_LEN(cmode));

  *result = ccp4_utils_chmod (temp_path, temp_cmode);
  free(temp_path);
  free(temp_cmode);
}

/* isatty doesnt seem to be in Mircrosoft Visual Studdio so this is a fudge */
#if defined (CALL_LIKE_MVS)
# if CALL_LIKE_MVS == 1
int __stdcall ISATTY (int *lunit)
{
  return 0;
}

/* erfc doesnt seem to be in Mircrosoft Visual Studdio so this is a fudge */
float __stdcall ERFC(float *value)
{
  return (float) ccp4_erfc( (double) *value);
}

#else

int isatty_ (int *lunit)
{
  return 0;
}

float erfc_ (float *value)
{
  return (float) ccp4_erfc( (double) *value);
}

# endif
#endif

#if defined(F2C)
/* <f2c support>=                                                           */
int exit_ (status)
     int *status;
{
  f_exit ();                    /* may or may not be registered with
                                   exit, depending on the C libraries
                                   capabilities, but is idempotent */
  exit (*status);
}

int time_ ()
{
  return (int) time (NULL);
}

int getpid_ ()
{
  return (int) getpid ();
}

/* following are from libI77/fio.h */
#define MXUNIT 100
typedef struct
{       FILE *ufd;      /*0=unconnected*/
        char *ufnm;
        long uinode;
        int udev;
        int url;        /*0=sequential*/
        flag useek;     /*true=can backspace, use dir, ...*/
        flag ufmt;
        flag uprnt;
        flag ublnk;
        flag uend;
        flag uwrt;      /*last io was write*/
        flag uscrtch;
} unit;
extern unit f__units[];
#define TRUE_ (1)
#define FALSE_ (0)
#define err(f,m,s) {if(f) errno= m; else f__fatal(m,s); return(m);}
/* end of fio.h extract */

int isatty_ (lunit)
     int *lunit;
{
  if (*lunit>=MXUNIT || *lunit<0)
    err(1,101,"isatty");
  /* f__units is a table of descriptions for the unit numbers (defined
     in io.h) with file descriptors rather than streams */
  return (isatty(fileno((f__units[*lunit]).ufd)) ? TRUE_ : FALSE_);
}

/* FORTRAN gerror intrinsic */
int gerror_ (str, Lstr)
char *str;
int  Lstr;
{
  int i;

  if (errno == 0) {             /* Avoid `Error 0' or some such message */
    for (i=1; Lstr; i++)
      str[i] = ' ';
  } else {
    (void) strncpy (str, strerror (errno), Lstr);
    for (i = strlen (str); i < Lstr; i++) str[i] = ' ';  /* pad with spaces */
  }
  return 0;
}

/* FORTRAN IErrNo intrinsic */
int ierrno_ () {
  return errno;
}

int itime_ (array)
     int array[3];
{
     struct tm *lt;
     time_t tim;
     tim = time(NULL);
     lt = localtime(&tim);
     array[0] = lt->tm_hour; array[1] = lt->tm_min; array[2] = lt->tm_sec;
}
/* These ought to be intrinsic, but they should only be applied to          */
/* [[INTEGER]] arguments.  The types [[integer]] and [[logical]] are both   */
/* assumed to be [[int]].                                                   */
/*                                                                          */
/* <f2c support>=                                                           */
int /* integer */ ibset_ (a, b)
     int /* integer */ *a, *b;
{
  return (*a) | 1<<(*b);
}

int /* integer */ ibclr_ (a, b)
     int /* integer */ *a, *b;
{
  return (*a) & ~(1<<(*b));
}

int /* logical */ btest_ (a, b)
     int /* integer */ *a, *b;
{
  return ((((unsigned long) *a)>>(*b)))&1 ? TRUE_ : FALSE_;
}
#endif              /* F2C support  */

#if defined (__hpux) || defined (_AIX)
/* <AIX and HPUX support>=                                                  */

#ifdef _AIX
int isatty_ (int *fd) {
  return(isatty(*fd));
}
#endif

void gerror  (str, Lstr)
char *str;
int  Lstr;
{
  int i;

  if (errno == 0) {             /* Avoid `Error 0' or some such message */
    for (i=1; Lstr; i++)
      str[i] = ' ';
  } else {
    (void) strncpy (str, strerror (errno), Lstr);
    for (i = strlen (str); i < Lstr; i++) str[i] = ' ';  /* pad with spaces */
  }
} /* End of gerror (str, Lstr) */

int ierrno () {
  return errno;
}

#endif             /*  HPUX and AIX support */

#if !( defined(G95) || defined(GFORTRAN) || defined(F2C) ) && \
    ( ( defined(__linux__) && defined(_CALL_SYSV) ) || defined(__APPLE__) )
/* linuxppc xlf and apple xlf support */
void gerror_ (str, Lstr)
char *str;
int  Lstr;
{
  int i;

  if (errno == 0) {             /* Avoid `Error 0' or some such message */
    for (i=1; Lstr; i++)
      str[i] = ' ';
  } else {
    (void) strncpy (str, strerror (errno), Lstr);
    for (i = strlen (str); i < Lstr; i++) str[i] = ' ';  /* pad with spaces */
  }
} /* End of gerror (str, Lstr) */

int isatty_(int *iunit)
{
  return isatty(*iunit);
}
#endif /* end of linuxppc/apple xlf support */

#if defined (sun)
int isatty_(int *iunit)
{
  return isatty(*iunit);
}
#endif

/* neither gfortran or g95 have isatty */
/* not true, since August 05 this has been added to gfortran */

/* G95 support */
#if defined(G95)
int isatty_(int *iunit)
{
  return isatty(*iunit);
}
#endif
#if defined(G95) || defined (GFORTRAN)

/* FORTRAN gerror intrinsic */
int gerror_(str, Lstr)
char *str;
int  Lstr;
{
  int i;

  if (errno == 0) {             /* Avoid `Error 0' or some such message */
    for (i=1; Lstr; i++)
      str[i] = ' ';
  } else {
    (void) strncpy (str, strerror (errno), Lstr);
    for (i = strlen (str); i < Lstr; i++) str[i] = ' ';  /* pad with spaces */
  }
  return 0;
}

/* FORTRAN IErrNo intrinsic */
int ierrno_() {
  return errno;
}

void ltime_(int *stime, int tarray[9])
{
  int i;
  struct tm ldatim;
  time_t t = *stime;

#ifdef __MINGW32__ // no localtime_r in MinGW
  struct tm* lt = localtime(&t);
  if (lt != NULL) {
    ldatim = *lt;
#else
  if (localtime_r(&t, &ldatim) != NULL) {
#endif
    tarray[0] = ldatim.tm_sec;
    tarray[1] = ldatim.tm_min;
    tarray[2] = ldatim.tm_hour;
    tarray[3] = ldatim.tm_mday;
    tarray[4] = ldatim.tm_mon;
    tarray[5] = ldatim.tm_year;
    tarray[6] = ldatim.tm_wday;
    tarray[7] = ldatim.tm_yday;
    tarray[8] = ldatim.tm_isdst;
  } else {
    for (i=0; i<9; i++)
      tarray[i] = 0;
  }

}

void idate_ (int *day, int *month, int *year)
{
     struct tm *lt=NULL;
     time_t tim;
     tim = time(NULL);
     lt = localtime(&tim);
     *day = lt->tm_mday;
     *month = lt->tm_mon+1;  /* need range 1-12 */
     *year = lt->tm_year + 1900;
}

void gmtime_(int *stime, int gmarray[9])
{
  int i;
  struct tm udatim;
  time_t t = *stime;

#ifdef __MINGW32__ // no gmtime_r in MinGW
  struct tm *p = gmtime(&t);
  if (p != NULL) {
    udatim = *p;
#else
  if (gmtime_r(&t, &udatim) != NULL) {
#endif
    gmarray[0] = udatim.tm_sec;
    gmarray[1] = udatim.tm_min;
    gmarray[2] = udatim.tm_hour;
    gmarray[3] = udatim.tm_mday;
    gmarray[4] = udatim.tm_mon;
    gmarray[5] = udatim.tm_year;
    gmarray[6] = udatim.tm_wday;
    gmarray[7] = udatim.tm_yday;
    gmarray[8] = udatim.tm_isdst;
  } else {
    for (i=0; i<9; i++)
      gmarray[i] = 0;
  }

}

void system_(int *status, char *cmd, int cmd_len)
{
  char *str = calloc( cmd_len+1, sizeof(char));
  str = strncpy( str, cmd, cmd_len);

  if ( (*status = system( str)) == -1 )
     printf(" Forked command %s failed\n",cmd);

  free( str);
  return;
}

#endif

#if defined (G95)
int time_()
{
  int ltim;
  time_t t_ltim;

  t_ltim = time(NULL);
  ltim = (int) t_ltim;

  return ltim;
}

#endif /* G95 support */
