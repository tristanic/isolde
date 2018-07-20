/*! \file lib/clipper_sysdep.h
    Header file for clipper system dependencies
*/
//C Copyright (C) 2000-2006 Kevin Cowtan and University of York
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA


#ifndef CLIPPER_SYSDEP
#define CLIPPER_SYSDEP


/* cmath is so variable, that we ignore it and use the C version:
 - include the math library
 - mirror into the 'std' namespace as below:
*/

#if defined(sun) || defined(sgi) || defined(__osf__) || defined(_MSC_VER)
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979
#endif
//#define rint(x) (floor((x)+0.5))
namespace std { using ::floor; using ::ceil; using ::fabs; using ::fmod; using ::sqrt; using ::sin; using ::cos; using ::tan; using ::asin; using ::acos; using ::atan; using ::sinh; using ::cosh; using ::tanh; using ::atan2; using ::exp; using ::log; using ::pow; }
#else
#include <cmath>
#endif

/* fixes for unexpected macros */

#if defined(isnan)
#undef isnan
#endif


/* numeric types, for where they are critical */

namespace clipper {
  typedef float     ftype32;
  typedef double    ftype64;
  typedef int       itype32;
  typedef unsigned int uitype32;
# define CLIPPER_NAN_MASK_A_32 0x7f800000U
# define CLIPPER_NAN_MASK_B_32 0x007fffffU
# define CLIPPER_NULL_MASK_32  0x7fc00000U
#if defined(__osf__) || defined(__amd64__)
  typedef long itype64;
  typedef unsigned long uitype64;
# define CLIPPER_NAN_MASK_A_64 0x7ff0000000000000UL
# define CLIPPER_NAN_MASK_B_64 0x000fffffffffffffUL
# define CLIPPER_NULL_MASK_64  0x7ff8000000000000UL
#else
  typedef long long itype64;
  typedef unsigned long long uitype64;
# define CLIPPER_NAN_MASK_A_64 0x7ff0000000000000ULL
# define CLIPPER_NAN_MASK_B_64 0x000fffffffffffffULL
# define CLIPPER_NULL_MASK_64  0x7ff8000000000000ULL
#endif
}


/* threading libraries and definitions */

#ifndef CLIPPER_DISABLE_THREADS

#if defined(__WIN32__) || defined(_WIN32)
// Through this header file, all windows macros are included in all programs
// that use clipper. Avoid some conflicts.
#define HKL HKL_RENAMED

// adding WIN32_LEAN_AND_MEAN makes some or all #undefs not necessary,
// but for now I'm leaving older workarounds unchanged
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN 1
#endif

#include <windows.h>

#undef small
#undef CompareString
#undef NO_ERROR
#undef HKL

#undef AddAtom
#undef GetAtomName
#undef DeleteAtom
#undef FindAtom

#define CLIPPER_MUTEX_INIT(MUTEX) InitializeCriticalSection(MUTEX)
#define CLIPPER_MUTEX_FREE(MUTEX) DeleteCriticalSection(MUTEX)
#define CLIPPER_MUTEX_LOCK(MUTEX) EnterCriticalSection(MUTEX)
#define CLIPPER_MUTEX_UNLK(MUTEX) LeaveCriticalSection(MUTEX)
#define CLIPPER_MUTEX_TYPE        CRITICAL_SECTION
#define CLIPPER_THREAD_EXEC(THREAD,ENTRY,ARG) (THREAD=CreateThread(0,0,(LPTHREAD_START_ROUTINE)ENTRY,(void*)ARG,0,NULL),THREAD!=NULL)
#define CLIPPER_THREAD_JOIN(THREAD)           (WaitForSingleObject(THREAD,INFINITE)>=0)
#define CLIPPER_THREAD_TYPE       HANDLE
#define CLIPPER_THREAD_ARGTYPE    LPVOID
#define CLIPPER_THREAD_RETTYPE    DWORD
#else

#include <pthread.h>
#define CLIPPER_MUTEX_INIT(MUTEX) pthread_mutex_init(MUTEX,NULL)
#define CLIPPER_MUTEX_FREE(MUTEX) pthread_mutex_destroy(MUTEX)
#define CLIPPER_MUTEX_LOCK(MUTEX) pthread_mutex_lock(MUTEX)
#define CLIPPER_MUTEX_UNLK(MUTEX) pthread_mutex_unlock(MUTEX)
#define CLIPPER_MUTEX_TYPE        pthread_mutex_t
#define CLIPPER_THREAD_EXEC(THREAD,ENTRY,ARG) (pthread_create(&THREAD,NULL,ENTRY,(void*)ARG)==0)
#define CLIPPER_THREAD_JOIN(THREAD)           (pthread_join(THREAD,NULL)==0)
#define CLIPPER_THREAD_TYPE       pthread_t
#define CLIPPER_THREAD_ARGTYPE    void*
#define CLIPPER_THREAD_RETTYPE    void*

#endif
#endif

#endif
