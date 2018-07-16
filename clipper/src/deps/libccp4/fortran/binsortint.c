/*
     binsortint.c: binary sorting functions

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
/*************************************************************************
	binsortint.c
	Z131093


IMPORTANT

This version is compatible with binsort version Z130891 and later.
Agains version Z080389 arguments of SRTBEG have been changed.
Modify and recompile your programs.

WHAT IS binsortint

binsortint is a set of routines used as an interface to
binsort mainly from FORTRAN programs.
Link this module in form of object module with your program.

CALLS:
  SRTBEG		sort initialisation
  SRTRLS		puts one record into sort
  SRTMRG		finishes sequence of input records
  SRTRET		gets one record from sort

For information about key type values see binsortkey.h

				Good Luck	J. Zelinka
***************************************************************************/

/* Probably, these function have never been used on Windows. */
#ifndef _WIN32

#include "binsort.h"

#include <stdio.h>
#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#ifndef NOUNISTD		/* ESV, for instance doesn't have it */
#  include <unistd.h>
#endif
#include <stddef.h>
#include <errno.h>
#include <sys/wait.h>
#include <sys/types.h>		/* necessary on Titan, at least, and
				   for POSIX wait */

#if __STDC__ && !defined(PROTOTYPE)
#define PROTOTYPE 1
#endif

#define	TRUE                    1
#define	FALSE                   0

#ifndef EIO
#  define EIO 5
#endif

static int	fildesout[2],	/* file descriptors for "r" & "w" access */
		fildesin[2];
static FILE	*filout,        /* file pointers for fildesout, fildesin */
		*filin;
static int	pid;		/* child process id */
static int	recl;		/* record length (fixed) */
static int      f_verbose = FALSE;      /* verbose flag */

/*=======================================================
SRTBEG:	Sort initialisation

	INTEGER		SRTBEG		(function) returns 0 if O.K.
                                        errno ot -1 otherwise
	INTEGER		NKYES		number of keys
	INTEGER 	KEYBUF(at least # of keys * 5)	keys description
	INTEGER		LRECL		(FIXED) length of record (BYTES)
	INTEGER		MEMSIZE		size of memory used by sort (BYTES)
						if .EQ. 0, default size

	I = SRTBEG(NKEYS,KEYBUF,LRECL,MEMSIZE)  I.EQ.0 status O.K.

KEYBUF consist of NKEYS entries, each of the form:
        INTEGER 	KEYTYP		key type (see binsortproc.h)
	INTEGER 	ORDER		sort order ( - || - )
	INTEGER 	KEYPOS		position within record (BYTES, 1st is
						position 0)
	INTEGER 	KEYLNG		length (data units chars, shorts ...)
	INTEGER         MASK            mask applied to data element before
	                                comparison (if 0 no mask applied)
=======================================================*/

/* F2C isn't disjoint with the rest, so be careful (also for other
   routines). */
#if defined (PROTOTYPE)

#if defined (sgi) || \
    defined(__OSF1__) || defined(__osf__) || defined(F2C) || \
    defined(G77) || defined(linux) || defined(__linux__) || defined (sun) || \
    defined(__APPLE__)
  int srtbeg_ (int *nkeys, int *keybuf, int *lrecl, int *memsize)

#else

#if defined (__hpux) || defined (_AIX) || defined (___AIX) 
  int srtbeg (int *nkeys, int *keybuf, int *lrecl, int *memsize)
#endif

#endif

/* no PROTOTYPE */
#else

#if defined (sgi) || \
    defined(__OSF1__) || defined(__osf__) || defined(F2C) || \
    defined(G77) || defined(linux) || defined(__linux__) || defined (sun) || \
    defined(__APPLE__)
  int srtbeg_ (nkeys, keybuf, lrecl, memsize)

#else

#if defined (__hpux) || defined (_AIX) || defined (___AIX)
  int srtbeg (nkeys, keybuf, lrecl, memsize)
#endif

#endif

int     	*keybuf;	/* keys description */
int	        *lrecl;		/* length of record */
int	        *nkeys;		/* number of keys */
int	        *memsize;       /* size of memory (BYTES) used by sort */

#endif
{
  char         **argv, **pargv;
  int            argc;
  char          *binsortname;
  static char   *rclswitch = "-r";
  char           lengthbuf [10];
  static char   *memswitch = "-m";
  char           memsizbuf [10];
  int            numkeys, i;
  static char   *keyswitch = "-k";
  char           keydatatype;
  char           sortorder;
  char           *charvalue;

  if (charvalue = (char *) getenv ("BINSORT_VERB"))
    f_verbose = charvalue ? TRUE : FALSE;
  binsortname = "binsort";
  if (charvalue = (char *) getenv ("BINSORT_NAME"))
    binsortname = charvalue;

  pipe(fildesin);
  pipe(fildesout);

  recl = *lrecl;

  if ((pid = fork()) == 0 ) {	/* the child process */

    /* create input and output pipes */

    if ((close(0) != 0) ||		/* close stdin */
	(dup(fildesout[0]) < 0) ||
	(close(1) != 0) ||
	(dup(fildesin[1]) < 0) ||
	(close(fildesout[1]) != 0) ||
	(close(fildesin[0]) != 0)) {
      perror("Binsort -- failed to handle I/O pipes to the parent process");
      _exit(1);
    }

    /* prepare binsort command line */

    argc = (*nkeys * 2) + 2 + 1 + 1;
    if (*memsize)
      argc += 2;
    argv = (char **)malloc(argc * sizeof(char *));
    if (argv == NULL) {
      fprintf(stderr,"malloc failed in SRTBEG\n");
      _exit(1);
    }
    argv [argc-1] = (char *)NULL;
    pargv = argv;
    *pargv++ = binsortname;
    *pargv++ = rclswitch;                  /* -r length */
    sprintf(lengthbuf, "%d", recl);
    *pargv++ = lengthbuf;
    if (*memsize) {
      *pargv++ = memswitch;                /* -m size */
      sprintf(memsizbuf, "%d", *memsize);
      *pargv++ = memsizbuf;
    }
    for (numkeys = *nkeys; numkeys; --numkeys, keybuf += 5) {
      *pargv++ = keyswitch;                /* -k keydescription */
      *pargv = (char *)malloc(256);
      if (*pargv == NULL) {
	fprintf(stderr,"malloc failed in SRTBEG\n");
	_exit(1);
      }
      switch (keybuf [0]) {
	case CHAR:     keydatatype = 'c';
	               break;
        case UCHAR:    keydatatype = 'C';
	               break;
        case SHORT:    keydatatype = 's';
	               break;
	case USHORT:   keydatatype = 'S';
	               break;
	case LONG:     keydatatype = 'l';
	               break;
	case ULONG:    keydatatype = 'L';
	               break;
	case FLOAT:    keydatatype = 'f';
	               break;
	case DOUBLE:   keydatatype = 'd';
	               break;
	}
      sortorder = (keybuf [1] == ASC) ? 'a' : 'd';
      switch (keybuf [0]) {
      case UCHAR:
      case USHORT:
      case ULONG:  sprintf(*pargv, "%c:%c:%d:%d:%x", keydatatype, sortorder,
			   keybuf [2], keybuf [3], keybuf [4]);
	           break;
      default:     sprintf(*pargv, "%c:%c:%d:%d", keydatatype, sortorder,
			   keybuf [2], keybuf [3]);
      }
      ++pargv;
    }

    if (f_verbose)
      fprintf(stderr, " binsortint -- calling program \"%s\"\n",binsortname);
    if (f_verbose) {
        for (i=0; i < argc ;++i) {
	  fprintf(stderr, " binsortint -- argument #%d = \"%s\"\n",i,argv[i]);
	}
    }
    execvp(binsortname, argv);
    perror("Trying to execute binsort");
    _exit(errno);
  }
  else if (pid == -1) {		/* fork failed */
    perror("Trying to fork for binsort");
    return(errno);
  }
  else {                         /* THE PARENT */
    close(fildesout[0]);
    close(fildesin[1]);
    if (!(filout = fdopen(fildesout[1], "w")))
      return(EIO);
  }
  return(0);
}

/*=======================================================
SRTRLS:	Release one record into Sort
	INTEGER		                             SRTRLS  (function)
	CHARACTER*(at least length of the record)    RECORD   pointer to record

	I = SRTRLS(RECORD)		I.EQ.0 status O.K.
                                        errno otherwise
=======================================================*/

#if defined (PROTOTYPE)

#if defined (sgi) || \
    defined(__OSF1__) || defined(__osf__) || defined(F2C) || \
    defined(G77) || defined(linux) || defined(__linux__) || defined (sun) || \
    defined(__APPLE__)
  int srtrls_ (char *record)

#else

#if defined (__hpux) || defined (_AIX) || defined (___AIX)
  int srtrls (char *record)
#endif

#endif

/* no PROTOTYPE */
#else

#if defined (sgi) || \
    defined(__OSF1__) || defined(__osf__) || defined(F2C) || \
    defined(G77) || defined(linux) || defined(__linux__) || defined (sun) || \
    defined(__APPLE__)
  int srtrls_ (record)

#else

#if defined (__hpux) || defined (_AIX) || defined (___AIX)
  int srtrls (record)
#endif

#endif

char            *record;
#endif
{
  register size_t ret;

  ret = fwrite(record, sizeof(char), recl, filout);
  return(ret == recl ? 0 : ferror(filout));
}

/*=======================================================
SRTMRG:	Merge - finish release phase
	INTEGER		SRTMRG          (function)

	I = SRTMRG()			I.EQ.0 status O.K.
                                        errno otherwise
=======================================================*/


#if defined (sgi) || \
    defined(__OSF1__) || defined(__osf__) || defined(F2C) || \
    defined(G77) || defined(linux) || defined(__linux__) || defined (sun) || \
    defined(__APPLE__)
  int srtmrg_ ()

#else

#if defined (__hpux) || defined (_AIX) || defined (___AIX)
  int srtmrg ()
#endif

#endif
{
    fclose(filout);
    if (!(filin = fdopen(fildesin[0], "r")))
      return(EIO);
    return(0);
}


/*=======================================================
SRTRET:	Return 1 record from sort
	INTEGER		                             SRTRET  (function)
	CHARACTER*(at least length of the record)    RECORD   pointer to record

	I = SRTRET(RECORD)		I.EQ.0 status O.K.
					I.EQ. -1 End of records
					errno otherwise
=======================================================*/

#if defined (PROTOTYPE)
#if defined (sgi) || \
    defined(__OSF1__) || defined(__osf__) || defined(F2C) || \
    defined(G77) || defined(linux) || defined(__linux__) || defined (sun) || \
    defined(__APPLE__)
  int srtret_ (char *record)

#else

#if defined (__hpux) || defined (_AIX) || defined (___AIX)
  int srtret (char *record)
#endif

#endif

/* no PROTOTYPE */
#else

#if defined (sgi) || \
    defined(__OSF1__) || defined(__osf__) || defined(F2C) || \
    defined(G77) || defined(linux) || defined(__linux__) || defined (sun) || \
    defined(__APPLE__)
  int srtret_ (record)

#else

#if defined (__hpux) || defined (_AIX) || defined (___AIX)
  int srtret (record)
#endif

#endif

char            *record;
#endif
{
    register size_t	ret;
    int reterr;

    int status=0;

    if ((ret = fread(record, sizeof(char), recl, filin)) == recl) {
      return(0);
    }
    /* else EOF or read error */
    if ((int) wait (&status) < 0) { /* some error with sub-process */
      fclose(filin);
      return (EIO);
    }
    if (feof(filin) && status == 0
	&& ret == 0) {		/* ensure record not truncated */
      fclose(filin);
      return(-1);		/* normal EOF */
    }
    reterr=ferror(filin);
    fclose(filin);
    if (status != 0)
      return (status);		/* sub-process abended */
    if (reterr != 0) {
      return (reterr);
    } else			/* e.g. premature EOF */
      return (EIO);
}
#endif /* _WIN32 */
