/*
     ccp4_diskio_f.c FORTRAN API for file i/o.
     Copyright (C) 2002  CCLRC, Charles Ballard and Martyn Winn

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

/** @page diskio_f_page Fortran API for low level input/output. 
 *
 *  @section diskio_f_file_list File list

<ul>
<li>ccp4_diskio_f.c
</ul>

 *  @section diskio_f_overview Overview

This library consists of a set of wrappers to C functions which
perform random access input/output (I/O) on various data items, including 
bytes, to stream-mode files.

 */

/** @file ccp4_diskio_f.c
 *  FORTRAN API for file i/o.
 *  Charles Ballard and Martyn Winn
 */

/* FORTRAN API for library.c                                              
 *                                                                        
 * revisions:                                                             
 *           (4/5/01) C.Ballard                                          
 *                    respond to first steps to make library.c         
 *                    more "C" like.                                   
 *           (21/8/01) C.Ballard                                       
 *                     error catching from library.c                   
 *                                                                     
 *                                                                      
 * Portability and Code
 *                                                                      
 * System dependent names are handled in the FORTRAN_SUBR,              
 * FORTRAN_FUN, FORTRAN_CALL macros defined in the header file.         
 * fpstr is a typedef which masks the intricacies of FORTRAN string     
 * passing.                                                             
 */                                                                      

#include <string.h>
#include "../ccp4/ccp4_utils.h"
#include "../ccp4/ccp4_errno.h"
#include "../ccp4/ccp4_fortran.h"
#include "../ccp4/ccp4_file_err.h"
/* rcsid[] = "$Id$" */


/**
 * _ioChannels:
 * structure to hold files 
 */
typedef struct _CCP4IObj CCP4IObj;

enum FILE_KLASS {NONE,CCP4_FILE,CCP4_MAP};

struct _CCP4IObj {
  enum FILE_KLASS klass;
  CCP4File *iobj;
};

static CCP4IObj *_ioChannels[MAXFILES];

int _get_channel()
{
  int i;
  for ( i = 1; i < MAXFILES ; i++) 
    if (!_ioChannels[i]) return i;
  return -1;
}

CCP4IObj *_iobj_init()
{
  return (CCP4IObj *) malloc(sizeof(CCP4IObj));
}

static int file_attribute[] = { /* DISKIO file modes */
  O_RDWR | O_TRUNC,   /* 'UNKNOWN'   open as 'OLD'/'NEW' check existence */
  O_TMP | O_RDWR | O_TRUNC,   /* 'SCRATCH'   open as 'OLD' and delete on closing */
  O_RDWR,   /* 'OLD'       file MUST exist or program halts */
  O_RDWR | O_TRUNC,   /* 'NEW'       create (overwrite) new file */
  O_RDONLY     /* 'READONLY'  self explanatory */
};

FORTRAN_SUBR ( QOPEN, qopen,
    (int *iunit, fpstr lognam, fpstr atbuta, int lognam_len, int atbuta_len),
    (int *iunit, fpstr lognam, fpstr atbuta),
    (int *iunit, fpstr lognam, int lognam_len, fpstr atbuta, int atbuta_len))
{
  char *atbut2, *temp_lognam, *fname;
  int istat;

  atbut2 = ccp4_FtoCString(FTN_STR(atbuta), FTN_LEN(atbuta));
  
  switch (atbut2[0]) {
  case 'U':
  case 'u':
    istat = 0; 
    break;
  case 'S':
  case 's':
    istat = 1; 
    break;
  case 'O':
  case 'o':
    istat = 2; 
    break;
  case 'N':
  case 'n':
    istat = 3; 
#ifndef _MSC_VER
    if (strcasecmp(getenv("CCP4_OPEN"),"UNKNOWN"))
#else
    if (_stricmp(getenv("CCP4_OPEN"),"UNKNOWN"))
#endif
      istat = 0;
    break;
  case 'R':
  case 'r':
    istat = 4; 
    break;
  default:
    ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_BadMode),
                "QOPEN", NULL);
  }
  if (atbut2) free(atbut2);

  if ((*iunit = _get_channel()) == -1)
   ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_MaxFile),
                "COPEN1", NULL);

  _ioChannels[*iunit] = _iobj_init();

  temp_lognam = ccp4_FtoCString(FTN_STR(lognam), FTN_LEN(lognam));
  if (!(fname = getenv(temp_lognam))) 
    fname = temp_lognam;

  if (!(_ioChannels[*iunit]->iobj = ccp4_file_open (fname, 
                 file_attribute[istat]) ) ) {
    printf("  Can't open file %s\n",fname);
    ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_CantOpenFile),
                "COPEN2", NULL);
  }

  _ioChannels[*iunit]->klass = CCP4_FILE;

  if(temp_lognam) free(temp_lognam);
}

FORTRAN_SUBR ( QQOPEN, qqopen,
    (int *iunit, fpstr lognam, const int *istat, int lognam_len),
    (int *iunit, fpstr lognam, const int *istat),
    (int *iunit, fpstr lognam, int lognam_len, const int *istat))
{
  char *fname, *temp_lognam;
  int jstat;
  temp_lognam = ccp4_FtoCString(FTN_STR(lognam), FTN_LEN(lognam));

  if (*istat < 1 || *istat > 5) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_BadMode),
                "QQOPEN (mode)", NULL);

  jstat = *istat;

  if (jstat == 4) 
#ifndef _MSC_VER
    if (strcasecmp(getenv("CCP4_OPEN"),"UNKNOWN"))
#else
    if (_stricmp(getenv("CCP4_OPEN"),"UNKNOWN"))
#endif
      jstat = 1;

  if (!(fname = getenv(temp_lognam))) 
    fname = temp_lognam;

  if ((*iunit = _get_channel()) == -1)
   ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_MaxFile),
                "QQOPEN", NULL);

  _ioChannels[*iunit] = _iobj_init();

  if (!(_ioChannels[*iunit]->iobj = ccp4_file_open (fname, 
		   file_attribute[jstat-1]) ) ) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_MaxFile),
                "QQOPEN", NULL);

  _ioChannels[*iunit]->klass = CCP4_FILE;

  if(temp_lognam) free(temp_lognam);
}

/**
 * Opens filename on io stream iunit. istat corresponds to the open mode. 
 * @param iunit iochannel number
 * @param filename fortran character array giving filename
 * @param istat file mode 
 */
FORTRAN_SUBR ( COPEN, copen,
    (int *iunit, fpstr filename, int *istat, int filename_len),
    (int *iunit, fpstr filename, int *istat),
    (int *iunit, fpstr filename, int filename_len, int *istat))
{
  char *tempfile;

  tempfile = ccp4_FtoCString(FTN_STR(filename), FTN_LEN(filename));

  if ((*iunit = _get_channel()) == -1)
   ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_MaxFile),
                "COPEN", NULL);

  _ioChannels[*iunit] = _iobj_init();

  if (!(_ioChannels[*iunit]->iobj = ccp4_file_open (tempfile, 
		   file_attribute[*istat-1]) ) ) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_MaxFile),
                "COPEN", NULL);

  _ioChannels[*iunit]->klass = CCP4_FILE;

  free(tempfile);
}

/** 
 * qrarch: 
 * @iunit: iochannel
 * @ipos: position in file
 * @ireslt: return value
 *                                                                          
 * For binary files with a well-determined structure in terms of            
 * [[float]]s and [[int]]s we may want to set up the connected stream to    
 * do transparent reading of files written on a machine with a different    
 * architecture.  This is currently the case for map files and   
 * MTZ files and this routine is called from mtzlib and maplib. 
 *                                                                          
 * qrarch reads the machine stamp at word ipos     
 * for the diskio file on stream iunit and sets up the appropriate   
 * bit-twiddling for subsequent qreads on that stream.  The             
 * information read from the file is returned in \meta{ireslt} in the       
 * form fileFT+16fileIT.  If the stamp is zero      
 * (as it would be for files written with a previous version of the         
 * library) we assume the file is in native format and needs no             
 * conversion in qread; in this case ireslt will be zero and     
 * the caller can issue a warning.  Iconvert and Fconvert are       
 * used by qread to determine the type of conversion (if any) to be     
 * applied to integers and reals.                                           
 *                                                                          
 * Fudge:fudge Ian Tickle reports old VAX files which have a machine
 * stamp which is byte-flipped from the correct VAX value,although it should
 * always have been zero as far as I can see.  To accommodate this, set the 
 * logical NATIVEMTZ and the machine stamp won't be read for any      
 * input files for which qrarch is called.                              
 *                                                                          
 * Extra feature: logical/environment variable CONVERT_FROM may be set  
 * to one of BEIEEE, LEIEEE, VAX or CONVEXNATIVE to avoid   
 * reading the machine stamp and assume the file is from the stipulated     
 * archictecture for all input MTZ and map files for which qrarch is
 * called.                                                                  
 *                                                                          
 * N.B.: leaves the stream positioned just after the machine stamp.         
 *                                                                          
 */
FORTRAN_SUBR ( QRARCH, qrarch,
    (int *iunit, int *ipos, int *ireslt),
    (int *iunit, int *ipos, int *ireslt),
    (int *iunit, int *ipos, int *ireslt))
{
  if (ccp4_file_setstamp(_ioChannels[*iunit]->iobj, *ipos)) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_BadMode),
                "QRARCH", NULL);
  if ((*ireslt = ccp4_file_rarch (_ioChannels[*iunit]->iobj)) == -1)
    ccp4_signal(CCP4_ERRLEVEL(4), "QRARCH", NULL);
}

/**
 * qwarch
 * @iunit: io channel
 * @ipos: position
 *
 * This is the complement of qrarch, writing the native machine        
 * architecture information machine stamp to diskio stream      
 * iunit at word ipos.  Currently called from mtzlib and maplib.
 *
 * The machine stamp in mtstring is four nibbles in order, indicating 
 * complex and real format (must both be the same), integer format and
 * character format (currently irrelevant).  The last two bytes of    
 * mtstring are currently unused and always zero.                 
 *                                                                    
 * N.B.: leaves the stream positioned just after the machine stamp.  
 * 
 */
FORTRAN_SUBR ( QWARCH, qwarch,
    (int *iunit, int *ipos),
    (int *iunit, int *ipos),
    (int *iunit, int *ipos))
{
  if (ccp4_file_setstamp(_ioChannels[*iunit]->iobj, *ipos))  
    ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_BadMode),
                "QWARCH", NULL); 
  if (ccp4_file_warch (_ioChannels[*iunit]->iobj) == -1)
    ccp4_signal(CCP4_ERRLEVEL(4), "QWARCH", NULL);
}

/**
 * qclose:
 * @iunit: io channel
 *
 * Closes the file open on diskio stream iunit
 */
FORTRAN_SUBR ( QCLOSE, qclose,
    (int *iunit),
    (int *iunit),
    (int *iunit))
{
  if (ccp4_file_close (_ioChannels[*iunit]->iobj))
    ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_CloseFail),
                "QCLOSE", NULL);
  free(_ioChannels[*iunit]);
  _ioChannels[*iunit]=NULL;
}

/**
 * qmode: 
 * @iunit: io channel
 * @mode: access mode
 * @size: item size
 *
 * Changes the diskio access mode for stream @iunit to  
 * @mode.  The resulting size in bytes of items for transfer is  
 * returned as @size
 *                   
 */
FORTRAN_SUBR ( QMODE, qmode,
    (int *iunit, int *mode, int *size),
    (int *iunit, int *mode, int *size),
    (int *iunit, int *mode, int *size))
{
  if ( (*size = ccp4_file_itemsize(_ioChannels[*iunit]->iobj)) == -1 ||
       ccp4_file_setmode(_ioChannels[*iunit]->iobj,*mode) == -1) 
   ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_BadMode),
                "QMODE", NULL);
}

/**
 * qread:
 * @iunit: io channel
 * @buffer:
 * @nitems: number of items
 * @result: return value
 *                                                                         
 * Reads @nitems in the current mode qmode() from diskio  
 * stream @iunit previously opened by qopen() and      
 * returns @result which is %0 on success or %-1 at EOF.     
 * It aborts on an i/o error.                                              
 * Numbers written in a foreign format will be translated if necessary if  
 * the stream is connected to an MTZ or map file.                          
 *                                                                        
 */
FORTRAN_SUBR ( QREAD, qread,
    (int *iunit, uint8 *buffer, int *nitems, int *result),
    (int *iunit, uint8 *buffer, int *nitems, int *result),
    (int *iunit, uint8 *buffer, int *nitems, int *result))
{
  *result = 0;
  if ( ccp4_file_read (_ioChannels[*iunit]->iobj, buffer, *nitems) != *nitems){
    if ( ccp4_file_feof(_ioChannels[*iunit]->iobj) ) 
      *result = -1;
    else 
      ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_ReadFail),
		  "QREAD", NULL);
  }
}

/**
 * qreadi:
 * @iunit: io channel
 * @buffer:
 * @result:
 *                                              
 * Fills INT buffer in int mode from diskio stream                      
 * @iunit previously opened by qopen() and returns  
 * @result} which %0 on success or %-1 on EOF.               
 * It aborts on an i/o failure.                                         
 * Call it with a character substring if necessary to control the number
 * of bytes read.                                                       
 *                                                                      
 */
FORTRAN_SUBR ( QREADI, qreadi,
    (int *iunit, uint8* buffer, int *nitems, int *result),
    (int *iunit, uint8* buffer, int *nitems, int *result),
    (int *iunit, uint8* buffer, int *nitems, int *result))
    {
  *result = 0;
  if ( ccp4_file_readint (_ioChannels[*iunit]->iobj, buffer, *nitems) != *nitems) {
    if ( ccp4_file_feof(_ioChannels[*iunit]->iobj) ) 
      *result = -1;
    else 
      ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_ReadFail),
		  "QREADI", NULL);
  }
}

/**
 * qreadi2:
 * @iunit: io channel
 * @buffer:
 * @result:
 *                                              
 * Fills INT*2 buffer in int mode from diskio stream                      
 * @iunit previously opened by qopen() and returns  
 * @result} which %0 on success or %-1 on EOF.               
 * It aborts on an i/o failure.                                         
 * Call it with a character substring if necessary to control the number
 * of bytes read.                                                       
 *                                                                      
 */
FORTRAN_SUBR ( QREADI2, qreadi2,
    (int *iunit, uint8* buffer, int *nitems, int *result),
    (int *iunit, uint8* buffer, int *nitems, int *result),
    (int *iunit, uint8* buffer, int *nitems, int *result))
    {
  *result = 0;
  if ( ccp4_file_readshort (_ioChannels[*iunit]->iobj, buffer, *nitems) != *nitems) {
    if ( ccp4_file_feof(_ioChannels[*iunit]->iobj) ) 
      *result = -1;
    else 
      ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_ReadFail),
		  "QREADI2", NULL);
  }
}

/**
 * qreadr:
 * @iunit:
 * @buffer:
 * @result:
 *
 * Fills REAL] buffer in int mode from diskio stream                   
 * @iunit previously opened by qopen() and returns  
 * @result which 0 on success or -1 on EOF.                
 * It aborts on an i/o failure.                                          
 * Call it with a character substring if necessary to control the number 
 * of bytes read.                                                        
 *                                                                       
 */
FORTRAN_SUBR ( QREADR, qreadr,
    (int *iunit, uint8* buffer, int *nitems, int *result),
    (int *iunit, uint8* buffer, int *nitems, int *result),
    (int *iunit, uint8* buffer, int *nitems, int *result))
    {
  *result = 0;
  if ( ccp4_file_readfloat (_ioChannels[*iunit]->iobj, buffer, *nitems) != *nitems) { 
    if ( ccp4_file_feof(_ioChannels[*iunit]->iobj) ) 
      *result = -1;
    else 
      ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_ReadFail),
		  "QREADR", NULL);
  }
}

/**
 * qreadq:
 * @iunit:
 * @buffer:
 * i@result:
 *
 * Fills COMPLEX buffer in int mode from diskio stream                 
 * @iunit previously opened by qopen() and returns 
 * @result which 0 on success or -1 on EOF.              
 * It aborts on an i/o failure.                                           
 * Call it with a character substring if necessary to control the number  
 * of bytes read.                                                         
 *
 */
FORTRAN_SUBR ( QREADQ, qreadq,
    (int *iunit, uint8* buffer, int *nitems, int *result),
    (int *iunit, uint8* buffer, int *nitems, int *result),
    (int *iunit, uint8* buffer, int *nitems, int *result))
    {
  *result = 0;
  if ( ccp4_file_readcomp (_ioChannels[*iunit]->iobj, buffer, *nitems) != *nitems) {
    if ( ccp4_file_feof(_ioChannels[*iunit]->iobj) ) 
      *result = -1;
    else 
      ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_ReadFail),
		  "QREADQ", NULL);
  }
}

/**
 * qreadc::
 * @iunit:
 * @buffer:
 * @result:
 *                                                                      
 * Fills CHARACTER buffer in byte mode from diskio stream           
 * @iunit previously opened by qopen() and returns  
 * @result which is 0 on success or -1 on EOF.               
 * It aborts on an i/o failure.                                         
 * Call it with a character substring if necessary to control the number
 * of bytes read.                                                       
 *
 */
FORTRAN_SUBR ( QREADC, qreadc,
    (int *iunit, fpstr buffer, int *result, int buffer_len),
    (int *iunit, fpstr buffer, int *result),
    (int *iunit, fpstr buffer, int buffer_len, int *result))
{
  int n;

  n = FTN_LEN(buffer);

  if (ccp4_file_readchar (_ioChannels[*iunit]->iobj, (uint8 *) FTN_STR(buffer), (size_t) n) != n)
      ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_ReadFail),
		  "QREADC", NULL);
  *result = 0;
}

/**
 * qwrite:
 * @iunit:
 * @buffer:
 * @meta:
 *
 * This write @nitems items from @buffer to qopen() 
 * stream \meta{iunit} using the current mode.                      
 *                                                                  
 */
FORTRAN_SUBR ( QWRITE, qwrite,
    (int *iunit, uint8 * buffer, int *nitems),
    (int *iunit, uint8 * buffer, int *nitems),
    (int *iunit, uint8 * buffer, int *nitems))
{
  if ( ccp4_file_write (_ioChannels[*iunit]->iobj, buffer, *nitems) != *nitems)
      ccp4_signal(CCP4_ERRLEVEL(4), "QWRITE", NULL);
}

/**
 * qwriti:
 * @iunit:
 * @buffer:
 * @nitems:
 *
 * This write @nitems items from @buffer to qopen()
 * stream @iunit using the INT32 mode.   
 *
 */
FORTRAN_SUBR ( QWRITI, qwriti,
    (int *iunit, uint8 * buffer, int *nitems),
    (int *iunit, uint8 * buffer, int *nitems),
    (int *iunit, uint8 * buffer, int *nitems))
{
  if ( ccp4_file_writeint (_ioChannels[*iunit]->iobj, buffer, *nitems) != *nitems)
      ccp4_signal(CCP4_ERRLEVEL(4), "QWRITI", NULL);
}

/**
 * qwritr:
 * @iunit:
 * @buffer:
 * @nitems:
 *
 * This write @nitems items from @buffer to qopen()  
 * stream @iunit using the FLOAT32 mode.                    
 *                                                                    
 */
FORTRAN_SUBR ( QWRITR, qwritr,
    (int *iunit, uint8 * buffer, int *nitems),
    (int *iunit, uint8 * buffer, int *nitems),
    (int *iunit, uint8 * buffer, int *nitems))
{
  if ( ccp4_file_writefloat (_ioChannels[*iunit]->iobj, buffer, *nitems) != *nitems)
      ccp4_signal(CCP4_ERRLEVEL(4), "QWRITR", NULL);
}

/**
 * qwrite:
 * @iunit:
 * @buffer:
 * @nitems:
 *
 * This write @nitems items from @buffer to qopen()
 * stream @iunit using the COMP64 mode. 
 *                                                
 */
FORTRAN_SUBR ( QWRITQ, qwritq,
    (int *iunit, uint8 * buffer, int *nitems),
    (int *iunit, uint8 * buffer, int *nitems),
    (int *iunit, uint8 * buffer, int *nitems))
{
  if ( ccp4_file_writecomp (_ioChannels[*iunit]->iobj, buffer, *nitems) != *nitems)
      ccp4_signal(CCP4_ERRLEVEL(4), "QWRITQ", NULL);
}

/* \subsection{{\tt subroutine qwritc (\meta{iunit}, \meta{buffer})}}       */
/*                                                                          */
/* Writes [[CHARACTER*(*)]] \meta{buffer} to [[qopen]]ed                    */
/* stream \meta{iunit} in byte mode.                                        */
/*                                                                          */
/* <diskio routines>=                                                       */
FORTRAN_SUBR ( QWRITC, qwritc,
    (int *iunit, fpstr buffer, int buffer_len),
    (int *iunit, fpstr buffer),
    (int *iunit, fpstr buffer, int buffer_len))
{
  int n;

  n = FTN_LEN(buffer);

  if (ccp4_file_writechar (_ioChannels[*iunit]->iobj, (uint8 *) FTN_STR(buffer), 
			   (size_t) n) != n)
      ccp4_signal(CCP4_ERRLEVEL(4), "WWRITC", NULL);
}

/**
 * qseek:
 * @iunit:
 * @irec:
 * @iel:
 * @lrecl:
 *
 * Seeks to element @iel in record @irec in diskio stream 
 * @iunit whose record length is @lrecl.                  
 *                                                                    
 */
FORTRAN_SUBR ( QSEEK, qseek,
    (int *iunit, int *irec, int *iel, int *lrecl),
    (int *iunit, int *irec, int *iel, int *lrecl),
    (int *iunit, int *irec, int *iel, int *lrecl))
{
  /*switch from FORTRAN offset to C offset */
  if (ccp4_file_seek (_ioChannels[*iunit]->iobj, (*irec-1)*(*lrecl)+(*iel-1),SEEK_SET) )
      ccp4_signal(CCP4_ERRLEVEL(4), "QSEEK", NULL);
}

/**
 * qback:
 * @iunit:
 * @lrecl:
 *
 * Backspaces one record, of length @lrecl on diskio stream @iunit.
 *                                                                 
 */
FORTRAN_SUBR ( QBACK, qback,
    (int *iunit, int *lrecl),
    (int *iunit, int *lrecl),
    (int *iunit, int *lrecl))
{
  if (ccp4_file_seek (_ioChannels[*iunit]->iobj, -(*lrecl), SEEK_CUR) )
      ccp4_signal(CCP4_ERRLEVEL(4), "QBACK", NULL);
}

/**
 * qskip:
 * @iunit:
 * @lrecl:
 *
 * Skip forward 1 record of length @lrecl on diskio stream @iunit.
 *       
 */
FORTRAN_SUBR ( QSKIP, qskip,
    (int *iunit, int *lrecl),
    (int *iunit, int *lrecl),
    (int *iunit, int *lrecl))
{
  if (ccp4_file_seek (_ioChannels[*iunit]->iobj, *lrecl, SEEK_CUR) )
      ccp4_signal(CCP4_ERRLEVEL(4), "QSKIP", NULL);
}

/**
 * qqinq:
 * @istrm:
 * @filnam:
 * @length:
 *
 * Returns the name @filnam and @length of the file (if any)
 * open on diskio stream @istrm.             
 *                                                 
 */
FORTRAN_SUBR ( QQINQ, qqinq,
    (int *istrm, fpstr logname, fpstr filename, int *length, 
     int logname_len, int filename_len),
    (int *istrm, fpstr logname, fpstr filename, int *length),
    (int *istrm, fpstr logname, int logname_len, fpstr filename, 
     int filename_len, int *length))
{
  char *log_name = NULL, *file_name;

  if ( *istrm < 1 || *istrm >= MAXFILES || !_ioChannels[*istrm]->iobj) {
    *length = -1;
    if (!(log_name = ccp4_FtoCString(FTN_STR(logname), 
				     FTN_LEN(logname))))
      log_name = strdup("diskio.dft"); 
    if (!(file_name = getenv(log_name)))
      file_name = log_name;
    for ( *istrm = 1; *istrm != MAXFILES; (*istrm)++)
      if (!strcmp(file_name,_ioChannels[*istrm]->iobj->name)) break;
  }
  if (*istrm != MAXFILES) {
    *length = ccp4_file_length(_ioChannels[*istrm]->iobj);
    strncpy(FTN_STR(filename), _ioChannels[*istrm]->iobj->name,
                    MIN(strlen(_ioChannels[*istrm]->iobj->name), 
			FTN_LEN(filename))); }
    
  if ( *length  == -1) 
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_SeekFail),
		  "QINQ", NULL);
  if (log_name != NULL) free(log_name);
}

/**
 * qlocate:
 * @iunit:
 * @locate:
 *
 * Returns the current position \meta{locate} in the diskio stream @iunit.
 *                                                                        
 */
FORTRAN_SUBR ( QLOCATE, qlocate,
    (int *iunit, int *locate),
    (int *iunit, int *locate),
    (int *iunit, int *locate))
{
  if ( (*locate = (int) ccp4_file_tell (_ioChannels[*iunit]->iobj) ) 
       == -1) 
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_SeekFail),
		  "QLOCATE", NULL);
}

