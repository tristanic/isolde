//  $Id: mmdb_machine.h $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2013.
//
//    This library is free software: you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License version 3, modified in accordance with the provisions
//    of the license to address the requirements of UK law.
//
//    You should have received a copy of the modified GNU Lesser
//    General Public License along with this library. If not, copies
//    may be downloaded from http://www.ccp4.ac.uk/ccp4license.php
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//  =================================================================
//
//    12.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :   Machine  <interface>
//       ~~~~~~~~~
//  **** Functions : mmdb::machine::GetMachineID   - returns ID code
//       ~~~~~~~~~~~                                 for the machine
//                   mmdb::machine::GetMachineName - returns name of
//                                                   the machine
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#ifndef  __MMDB_Machine__
#define  __MMDB_Machine__

#include "mmdb_mattype.h"

/*
// Programs written in plain C, should define __PlainC each time
// this header is invoked.

#ifdef __PlainC
# undef __CPlusPlus
#else
# define __CPlusPlus
#endif
*/

namespace mmdb  {

  namespace machine  {


    //  ==================  List of known machines  =================

    enum MACHINE  {
      MACHINE_SGI       = 1,
      MACHINE_RS6000    = 2,
      MACHINE_ALLIANT   = 3,
      MACHINE_ARDENT    = 4,
      MACHINE_TITAN     = 5,
      MACHINE_STARDENT  = 6,
      MACHINE_CONVEX    = 7,
      MACHINE_ESV       = 8,
      MACHINE_HP9000    = 9,
      MACHINE_SOLBOURNE = 10,
      MACHINE_SOLARIS   = 11,
      MACHINE_ALPHA     = 12,
      MACHINE_VMS       = 13,
      MACHINE_MVS       = 14,
      MACHINE_F2C_G77   = 15,
      MACHINE_LINUX     = 16,
      MACHINE_UNKNOWN   = 100
    };

    //  =============  Identification of the machine  ===============

// IBM Unix RS/6000
#if defined(_AIX) || defined(___AIX)
# define CALL_LIKE_HPUX 1

// Alliant
#elif defined(alliant)
# define CALL_LIKE_SUN  1

// Ardent, Stardent/Titan
#elif defined(ardent)
# define CALL_LIKE_STARDENT 1
#elif defined(titan)
# define CALL_LIKE_STARDENT 2
#elif defined(stardent)
# define CALL_LIKE_STARDENT 3

// Convex
#elif defined(__convex__) || defined(__convexc__)
# define CALL_LIKE_SUN  2

// Evans and Sutherland
#elif defined(ESV)
# define CALL_LIKE_SUN  3

// Hewlett Packard 9000/750 (RISC) models
#elif defined(__hpux)
# define CALL_LIKE_HPUX 2

// Silicon Graphics IRIX systems, Iris'es, Indigo's, Crimson's etc.
#elif defined(__sgi) || defined(sgi)
# define CALL_LIKE_SUN  4

// Solbourne's are Sun clones.
#elif defined(solbourne)
# define CALL_LIKE_SUN  5

// Solaris 1 and 2
#elif defined(sun) || defined(__sun)
# define CALL_LIKE_SUN  6

//  DEC, OSF/1, Alpha and Ultrix
#elif defined(ultrix) || defined(__OSF1__) || defined(__osf__)
# define CALL_LIKE_SUN  7

// VMS
#elif defined(vms) || defined(__vms) || defined(__VMS)
# define CALL_LIKE_VMS  1

// MVS stands for Microsoft Visual Studio
#elif defined(_MVS)
# define CALL_LIKE_MVS  1

#elif defined(F2C) || defined(G77)
# define CALL_LIKE_SUN  8

#elif defined(linux)
# define CALL_LIKE_SUN  9

#else
//# error System type is not known -- see the Installation Guide
# define CALL_LIKE_SUN  100

#endif



//  =================  Machine-dependent definitions  ==================

#ifdef CALL_LIKE_STARDENT
    // StrPar is used in Ardent-like machines' fortran calls
    // for passing a string parameter
    DefineStructure(StrPar)
    struct StrPar  {
      pstr S;
      int  len;
      int  id;
    };
#endif


//
//   Macro  FORTRAN_SUBR(NAME,name,p_send,p_struct,p_sflw)
// makes function header statements that allow for linking with
// programs written in FORTRAN.
//
//   Parameters:
//
//   NAME      name of the FORTRAN subroutine in capital letters
//   name      name of the FORTRAN subroutine in small letters
//   p_send    parameter list (in brackets) with string lengths
//             attached to the end of it (see below)
//   p_struct  parameter list (in brackets) with strings passed
//             as complex parameters, or structures
//   p_sflw    parameter list (in brackets) with string lengths
//             following immediately the string parameters
//             (see below)
//
//   All non-string parameters must be passed as pointers, in
// the same order as they enter the FORTRAN call. Rules for
// the string parameters are as follows.
//
//   1. All strings should be specified as of 'fpstr' type.
//      The 'fpstr' type is defined below and depends on the
//      platform:
//
//        a) whenever length of string is passed as a separate
//           parameter ( CALL_LIKE_SUN, CALL_LIKE_HPUX,
//           CALL_LIKE_MVS )  'fpstr' is identical to 'pstr'.
//           You may choose arbitrary name for the string,
//           but you MUST use the same name, appended with
//           suffix '_len', for its length (see example below).
//
//        b) whenever string and its length are passed as
//           complex parameter, 'fpstr' is identical to the
//           pointer on the corresponding structure:
//             CALL_LIKE_STARDENT :
//                 'fpstr' is identical to 'PStrPar'
//             CALL_LIKE_VMS      :
//                 'fpstr' is identical to 'dsc$descriptor_s *'
//
//      With 'fpstr' type, two important macro definition come:
//
//        i)  FTN_STR(s)  - returns pointer to fortran-passed
//                          string s. This pointer is always
//                          of 'pstr' type
//        ii) FTN_LEN(s)  - returns integer length of fortran-
//                          passed string s. For this macro to
//                          work properly with SUN- and MVS-like
//                          machines, always use suffix '_len'
//                          for the string length parameters as
//                          described in a) above.
//
//   2. Three parameter lists, each enclosed in brackets, should
//      be given. These lists retain the general order of
//      parameters in the corresponding fortran call. Non-string
//      parameters are passed as pointers. String parameters
//      and their lengths are passed differently in different
//      lists:
//
//       p_send    strings enter their place in the list as in
//                 the corresponding FORTRAN call, having 'fpstr'
//                 parameter type. Their lengths are appended as
//                 'int' to the end of the list. They should
//                 retain the order in which the strings appear
//                 in the list.
//
//       p_struct strings enter their place in the list as in
//                 the corresponding FORTRAN call, having 'fpstr'
//                 parameter type.
//
//       p_sflw    strings enter their place in the list as in
//                 the corresponding FORTRAN call, having 'fpstr'
//                 type and being immediately followed by their
//                 lengths as 'int' parameters.
//
//
//
// Example:
//
//   FORTRAN statement
//
//     subroutine  SomeSub ( k,s1,a,s2,m )
//     integer       k,m
//     real          a
//     character*(*) s1,s2
//
//   is translated to
//
//     FORTRAN_SUBR ( SOMESUB, somesub,
//       ( int * k, fpstr s1, float * a, fpstr s2, int * m,
//         int s1_len, int s2_len ),
//       ( int * k, fpstr s1, float * a, fpstr s2, int * m ),
//       ( int * k, fpstr s1, int s1_len, float * a,
//         fpstr s2, int s2_len, int * m ) )
//
//
//   The macro should replace ordinary function header
// statements to assure compatibility with FORTRAN links.
// In header files, do not forget to add semicolumn:
//
//   FORTRAN_SUBR ( .... );
//
// while in source files use simply
//
//   FORTRAN_SUBR ( .... )  {
//    <source body, operators>
//   }
//
//
//
//   Macro  FORTRAN_CALL(NAME,name,p_send,p_struct,p_sflw)
// calls function defined with macro FORTRAN_SUBR(...), from
// a C/C++ application. Its parameters and their meaning are
// exactly identical to those of FORTRAN_SUBR(...).
// FORTRAN_CALL(...) should be followed by semicolon.
//


//  **** type of real numbers in the API functions
//       comment or uncomment the proper string

    typedef  float      apireal;   // FORTRAN  real*4
/*
    typedef  double     apireal;    // FORTRAN  real*8
*/


#if defined(CALL_LIKE_SUN)

    typedef pstr fpstr;

# define FTN_STR(s)  s
# define FTN_LEN(s)  s##_len

# define char_struct(s)           \
    mmdb::pstr  s;                      \
    int   s##_len;
# define fill_char_struct(s,str)  \
    s  = str;                     \
    s##_len = strlen(str);

# ifdef __cplusplus
#   define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    extern "C" void name##_ p_sun
# else
#   define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    void name##_ p_sun
# endif

# define FORTRAN_EXTERN(NAME,name,p_sun,p_stardent,p_mvs) \
    extern "C" void name##_ p_sun

# define FORTRAN_CALL(NAME,name,p_sun,p_stardent,p_mvs) \
    name##_ p_sun

# elif defined(CALL_LIKE_HPUX)

    typedef pstr fpstr;

# define FTN_STR(s)  s
# define FTN_LEN(s)  s##_len

# define char_struct(s)  \
    pstr  s;             \
    int   s##_len;
# define fill_char_struct(s,str)  \
    s  = str;                     \
    s##_len = strlen(str);

# ifdef __cplusplus
#   define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    extern "C" void name p_sun
# else
#   define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    void name p_sun
# endif

# define FORTRAN_EXTERN(NAME,name,p_sun,p_stardent,p_mvs) \
    extern "C" void name p_sun

# define FORTRAN_CALL(NAME,name,p_sun,p_stardent,p_mvs) \
    name p_sun

#elif defined(CALL_LIKE_STARDENT)

    typedef PStrPar fpstr;

# define FTN_STR(s)  s->S
# define FTN_LEN(s)  s->len

# define char_struct(s)           \
    StrPar s;
# define fill_char_struct(s,str)  \
    s.S   = str;                  \
    s.len = strlen(FName);        \
    s.id  = 0;

# ifdef __cplusplus
#   define FORTRAN_SUBR(NAME,name,p_send,p_struct,p_sflw) \
    extern "C" void NAME p_stardent
# else
#   define FORTRAN_SUBR(NAME,name,p_send,p_struct,p_sflw) \
    void NAME p_stardent
# endif

# define FORTRAN_EXTERN(NAME,name,p_sun,p_stardent,p_mvs) \
    extern "C" void NAME p_stardent

# define FORTRAN_CALL(NAME,name,p_send,p_struct,p_sflw) \
    NAME p_stardent

#elif defined(CALL_LIKE_VMS)

    typedef dsc$descriptor_s * fpstr;

# define FTN_STR(s)  s->dsc$a_pointer;
# define FTN_LEN(s)  s->dsc$w_length;

# define character(s)                \
    dsc$descriptor_s s;
# define fill_char_struct(s,str)     \
    s.dsc$a_pointer = str;           \
    s.dsc$w_length  = strlen(str);   \
    s.dsc$b_dtype   = DSC$K_DTYPE_T; \
    s.dsc$b_class   = DSC$K_CLASS_S;

# ifdef __cplusplus
#   define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    extern "C" void NAME p_stardent
# else
#   define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    void NAME p_stardent
# endif

# define FORTRAN_EXTERN(NAME,name,p_sun,p_stardent,p_mvs) \
    extern "C" void NAME p_stardent

# define FORTRAN_CALL(NAME,name,p_sun,p_stardent,p_mvs) \
    NAME p_stardent

#elif defined(CALL_LIKE_MVS)

    typedef pstr fpstr;

# define FTN_STR(s)  s
# define FTN_LEN(s)  s##_len

# define char_struct(s)  \
    pstr  s;             \
    int   s##_len;
# define fill_char_struct(s,str)  \
    s  = str;                     \
    s##_len = strlen(str);

# ifdef __cplusplus
#   define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    extern "C" void __stdcall NAME p_mvs
# else
#   define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    void __stdcall NAME p_mvs
# endif

# define FORTRAN_EXTERN(NAME,name,p_sun,p_stardent,p_mvs) \
    extern "C" void NAME p_mvs

# define FORTRAN_CALL(NAME,name,p_sun,p_stardent,p_mvs) \
    NAME p_mvs

#else

# error  Unknown machine!!!

    typedef pstr fpstr;

# define FTN_STR(s)  s
# define FTN_LEN(s)  s##_len

# define char_struct(s)  \
    pstr  s;             \
    int   s##_len;
# define fill_char_struct(s,str)  \
    s  = str;                     \
    s##_len = strlen(str);

# ifdef __cplusplus
#   define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    extern "C" void name##_ p_sun
# else
#   define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    void name##_ p_sun
# endif

# define FORTRAN_EXTERN(NAME,name,p_sun,p_stardent,p_mvs) \
    extern "C" void name##_ p_sun

# define FORTRAN_CALL(NAME,name,p_sun,p_stardent,p_mvs) \
    name##_ p_sun

#endif


    //  ==============  Machine-dependent functions  ===============

    extern int         GetMachineID   ();
    extern mmdb::cpstr GetMachineName ();
    extern mmdb::cpstr GetMachineName ( int MachineID );

  }

}

#endif
