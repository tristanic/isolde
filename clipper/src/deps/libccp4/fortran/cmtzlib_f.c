/*
     cmtzlib_f.c: Fortran API to CMTZlib
     Copyright (C) 2001  CCLRC, Martyn Winn

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

/** @page cmtz_f_page Fortran API to CMTZ 
 *
 *  @section cmtz_f_file_list File list

<ul>
<li>cmtzlib_f.c
</ul>

 *  @section cmtz_f_overview Overview

This library consists of a set of wrappers to the CMTZ library
giving the same API as the original mtzlib.f For details of the
API, see the original <a href="../mtzlib.html">documentation</a>.
This document covers some peculiarities of the C implementation.

*   @section cmtz_f_batches Batches

Batch headers are held as a linked list in the MTZ data structure.
Upon reading an MTZ file, a list is created holding the all the batch
headers from the input file (mtz->n_orig_bat holds the number read
in). LRBATS will return a list of batch numbers. LRBAT can be used
to read the batch headers one at a time.
<p>
LWBAT can be used to write new or updated batch headers. Updated
batch headers are in fact added as new batch headers at the end of the
list of input batch headers. This mimics the use of RBATW/CBATW in
the Fortran routines. There are problems with this approach, and
this may be re-implemented.

 */
 
/** @file cmtzlib_f.c
 *
 *  @brief Fortran API for input, output and manipulation of MTZ files.
 *
 *  @author Martyn Winn 
 */

/*#define FORTRAN_CALL_DEBUG 1*/

#if defined (FORTRAN_CALL_DEBUG)
#  define CMTZLIB_DEBUG(x) x
#else
#  define CMTZLIB_DEBUG(x)
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../ccp4/ccp4_fortran.h"
#include "../ccp4/ccp4_utils.h"
#include "../ccp4/cmtzlib.h"
#include "../ccp4/csymlib.h"
#include "../ccp4/ccp4_program.h"
#include "../ccp4/ccp4_general.h"
/* rcsid[] = "$Id$" */

#define MFILES 9
#define MAXSYM 192

static int cmtz_in_memory = 0;
static int rlun[MFILES] = {0};
static int wlun[MFILES] = {0};
static MTZ *mtzdata[MFILES] = {NULL};         /* cf. Eugene's channel for rwbrook */
static char fileout[MFILES][MAXFLEN];
static int irref[MFILES] = {0};
static int iwref[MFILES] = {0};
static char user_label_in[MFILES][MCOLUMNS][2][31];
static char user_label_out[MFILES][MCOLUMNS][2][31];
static MTZCOL *collookup[MFILES][MCOLUMNS];
static MTZCOL *collookup_out[MFILES][MCOLUMNS];
static int sortorder[MFILES][5] = {{0}};
static int logmss[MFILES][MCOLUMNS];
static int ndatmss[MFILES];
static MTZBAT *rbat[MFILES];
static int nbatw[MFILES] = {0};
static double coefhkl[MFILES][6] = {{0}};

/* MVS defaults to int and doesn't like it */
#ifdef _MSC_VER
#if (CALL_LIKE_MVS==2)
void CCP4H_INIT_LIB(int *ihtml, int *isumm);
#else
 void __stdcall CCP4H_INIT_LIB(int *ihtml, int *isumm);
#endif
#endif

void MtzMemTidy(void) {

  int i;

  /* free any existing Mtz struct */
  for (i = 0; i < MFILES; ++i) 
    if (mtzdata[i]) MtzFree(mtzdata[i]);
}

/* Utility function for checking subroutine input 
   rwmode = 1 for read file, 2 for write file ... */
int MtzCheckSubInput(const int mindx, const char *subname, const int rwmode) {

 if (mindx <= 0 || mindx > MFILES) {
   printf("Error in %s: mindx %d out of range!\n",subname,mindx);
   return 1;
 }

 if (rwmode == 1 && rlun[mindx-1] == 0) {
   printf("Error in %s: mindx %d not open for read!\n",subname,mindx);
   return 1;
 }

 if (rwmode == 2 && wlun[mindx-1] == 0) {
   printf("Error in %s: mindx %d not open for write!\n",subname,mindx);
   return 1;
 }

 return 0;
}

/** Mainly for backwards compatibility. The only thing it does now is
 * initialise the html/summary stuff, which is usually done by CCPFYP
 * anyway.
 */
FORTRAN_SUBR ( MTZINI, mtzini,
               ( ),
               ( ),
               ( ))
{ int ihtml, isumm;

  /* initialise html/summary stuff */
  ihtml = 0;
  isumm = 0;
  FORTRAN_CALL (CCP4H_INIT_LIB, ccp4h_init_lib, (&ihtml,&isumm), (&ihtml,&isumm), (&ihtml,&isumm));
  return;
}

/** Open an MTZ file for reading. This is a wrapper to MtzGet.
 * @param mindx (I) MTZ file index
 * @param filename (I) Filename to open (real or logical)
 * @param iprint (I) Specifies how much header information to print out.
 * @param ifail (O) Returns 0 if successful, non-zero otherwise.
 */
FORTRAN_SUBR ( LROPEN, lropen,
               (int *mindx, fpstr filename, int *iprint, int *ifail, int filename_len),
               (int *mindx, fpstr filename, int *iprint, int *ifail),
               (int *mindx, fpstr filename, int filename_len, int *iprint, int *ifail))

{ 
  int i;
  char *temp_name, *fullfilename;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LROPEN");)

 *ifail = 0;

 if (*mindx <= 0 || *mindx > MFILES) {
   printf("Error: mindx out of range!\n");
   *ifail = 1;
   return;
 }

 if (rlun[*mindx-1] != 0) {
   printf("Error: mindx already used for read!\n");
   *ifail = 1;
   return;
 }

 if (wlun[*mindx-1] != 0) {
   printf("Error: mindx already used for write!\n");
   *ifail = 1;
   return;
 }

 temp_name = ccp4_FtoCString(FTN_STR(filename), FTN_LEN(filename));
 if (getenv(temp_name) != NULL) {
   fullfilename = strdup(getenv(temp_name));
 } else {
   fullfilename = strdup(temp_name);
 }

 /* In-memory mode not tested so well. Think it might be broken for
    some Fortran programs, e.g. when generating extra reflections. */
 if (getenv("CMTZ_IN_MEMORY")) cmtz_in_memory = 1;

 /* open file and read into memory */
 mtzdata[*mindx-1] = MtzGet(temp_name,cmtz_in_memory);
 if (mtzdata[*mindx-1] == NULL) {
   printf("Error: failed to open file for read!\n");
   *ifail = -1;
   free(fullfilename); 
   free(temp_name);
   return;
 }
 rlun[*mindx-1] = 1;

 /* Make sure HKL are in base dataset */
 MtzAssignHKLtoBase(mtzdata[*mindx-1]);

 ccp4printf(1,"\n OPENED INPUT MTZ FILE \n");
 ccp4printf(1," Logical Name: %s   Filename: %s \n\n",temp_name,fullfilename);
 if (*iprint > 0) ccp4_lhprt(mtzdata[*mindx-1], *iprint);

 /* set some control variables for Fortran calls */
 irref[*mindx-1] = 0;
 if (mtzdata[*mindx-1]->n_orig_bat > 0)
   rbat[*mindx-1] = mtzdata[*mindx-1]->batch;

 /* calculate hkl coefficients and store in coefhkl */
 for (i = 0; i < mtzdata[*mindx-1]->nxtal; ++i)
   if (mtzdata[*mindx-1]->xtal[i]->cell[0] > 0.001) {
     MtzHklcoeffs(mtzdata[*mindx-1]->xtal[i]->cell, coefhkl[*mindx-1]);
     break;
   }

 free(fullfilename); 
 free(temp_name);

}

/** Get title from MTZ file opened for read.
 * @param mindx (I) MTZ file index
 * @param ftitle (O) Title.
 * @param len (O) Length of returned title.
 */
FORTRAN_SUBR ( LRTITL, lrtitl,
	       (int *mindx, fpstr ftitle, int *len, int ftitle_len),
	       (int *mindx, fpstr ftitle, int *len),
	       (int *mindx, fpstr ftitle, int ftitle_len, int *len))
{
  char temp_title[71];

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRTITL");)
  if (MtzCheckSubInput(*mindx,"LRTITL",1)) return;

  *len = ccp4_lrtitl(mtzdata[*mindx-1], temp_title);
  ccp4_CtoFString(FTN_STR(ftitle),FTN_LEN(ftitle),temp_title);

}

/** Get history lines from MTZ file opened for read.
 * @param mindx (I) MTZ file index
 * @param hstrng (O) Array of history lines.
 * @param nlines (I/O) On input, dimension of hstrng, i.e. the maximum
 *   number of history lines to be returned. On output, actual number of 
 *   history lines returned.
 */
FORTRAN_SUBR ( LRHIST, lrhist,
	       (int *mindx, fpstr hstrng, int *nlines, int hstrng_len),
	       (int *mindx, fpstr hstrng, int *nlines),
	       (int *mindx, fpstr hstrng, int hstrng_len, int *nlines))

{

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRHIST");)

 if (MtzCheckSubInput(*mindx,"LRHIST",1)) return;

 *nlines = ccp4_lrhist(mtzdata[*mindx-1], hstrng, *nlines);

}

/** Fortran wrapper to functions returning number of reflections, columns
 * and ranges. In fact, it returns current values on MINDX rather than 
 * those of the input file.
 * @param mindx MTZ file index
 * @param versnx (O) MTZ version. This is the version of the current
    library rather than that of the input file. If these differ, a
    warning will have been issued by LROPEN/MtzGet.
 * @param ncolx Number of columns.
 * @param nreflx Number of reflections.
 * @param ranges Array of column ranges.
 */
FORTRAN_SUBR ( LRINFO, lrinfo,
	       (int *mindx, fpstr versnx, int *ncolx, int *nreflx, float *ranges, int versnx_len),
	       (int *mindx, fpstr versnx, int *ncolx, int *nreflx, float *ranges),
	       (int *mindx, fpstr versnx, int versnx_len, int *ncolx, int *nreflx, float *ranges))

{
  int i,j,k,iarray;
  char mtzvers[20];

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRINFO");)

 if (MtzCheckSubInput(*mindx,"LRINFO",1)) return;

  sprintf(mtzvers,"MTZ:V%d.%d",MTZ_MAJOR_VERSN,MTZ_MINOR_VERSN);
  ccp4_CtoFString(FTN_STR(versnx),FTN_LEN(versnx),mtzvers);
  *ncolx = MtzNumActiveCol(mtzdata[*mindx-1]);
  *nreflx = MtzNref(mtzdata[*mindx-1]);

  /* get ranges of active columns */
  iarray = 0;
  for (i = 0; i < mtzdata[*mindx-1]->nxtal; ++i) 
   for (j = 0; j < mtzdata[*mindx-1]->xtal[i]->nset; ++j) 
    for (k = 0; k < mtzdata[*mindx-1]->xtal[i]->set[j]->ncol; ++k) 
      if (mtzdata[*mindx-1]->xtal[i]->set[j]->col[k]->active) {
        ranges[iarray*2] = mtzdata[*mindx-1]->xtal[i]->set[j]->col[k]->min;
        ranges[iarray*2+1] = mtzdata[*mindx-1]->xtal[i]->set[j]->col[k]->max;
        iarray++;
      }
}

/** Fortran wrapper to function returning number of columns read in
 * from input file.
 * @param mindx MTZ file index
 * @param ncolx Number of columns.
 */
FORTRAN_SUBR ( LRNCOL, lrncol,
	       (int *mindx, int *ncolx),
	       (int *mindx, int *ncolx),
	       (int *mindx, int *ncolx))

{
  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRNCOL");)

 if (MtzCheckSubInput(*mindx,"LRNCOL",1)) return;

 *ncolx = MtzNumSourceCol(mtzdata[*mindx-1]);

}

/** Returns the current reflection number from an
 *  MTZ file opened for read. Files are normally read sequentially,
 *  the number returned is the number of the *NEXT* reflection
 *  record to be read. If you are going to jump about the file
 *  with LRSEEK then use this to record the current position before
 *  you start, so that it can be restored afterwards, if required.
 * @param mindx MTZ file index
 * @param nreflx the reflection record number of the next reflection to be read
 */
FORTRAN_SUBR ( LRNREF, lrnref,
	       (int *mindx, int *nreflx),
	       (int *mindx, int *nreflx),
	       (int *mindx, int *nreflx))

{
  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRNREF");)

 if (MtzCheckSubInput(*mindx,"LRNREF",1)) return;

  *nreflx = irref[*mindx-1] + 1;
}

/* Fortran wrapper for ccp4_lrsort */
FORTRAN_SUBR ( LRSORT, lrsort,
	       (int *mindx, int sortx[5]),
	       (int *mindx, int sortx[5]),
	       (int *mindx, int sortx[5]))

{
  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRSORT");)

 if (MtzCheckSubInput(*mindx,"LRSORT",1)) return;

  ccp4_lrsort(mtzdata[*mindx-1], sortx);

}

/** Fortran wrapper for ccp4_lrbats.
 * May be called for non-multirecord file, just to check there are 
 * no batches.
 * @param mindx MTZ file index
 * @param nbatx Number of batches found.
 * @param batchx Array of batch numbers.
 */
FORTRAN_SUBR ( LRBATS, lrbats,
	       (int *mindx, int *nbatx, int batchx[]),
	       (int *mindx, int *nbatx, int batchx[]),
	       (int *mindx, int *nbatx, int batchx[]))

{
  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRBATS");)

 if (MtzCheckSubInput(*mindx,"LRBATS",1)) return;

  ccp4_lrbats(mtzdata[*mindx-1], nbatx, batchx);

}

/* Fortran wrapper for MtzListColumn */
FORTRAN_SUBR ( LRCLAB, lrclab,
	       (int *mindx, fpstr clabs, fpstr ctyps, int *ncol, int clabs_len, int ctyps_len),
	       (int *mindx, fpstr clabs, fpstr ctyps, int *ncol),
	       (int *mindx, fpstr clabs, int clabs_len, fpstr ctyps, int ctyps_len, int *ncol))

{
 int i,j,k;
 char cclabs[MCOLUMNS][31]={{0}}, cctyps[MCOLUMNS][3]={{0}};
 int ccsetid[MCOLUMNS];

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRCLAB");)

 if (MtzCheckSubInput(*mindx,"LRCLAB",1)) return;

  *ncol = MtzListInputColumn(mtzdata[*mindx-1], cclabs, cctyps, ccsetid);

  for (i = 0; i < *ncol; ++i) {
    for (j = 0; j < clabs_len; ++j) {
      if (cclabs[i][j] == '\0') {
        for (k = j; k < clabs_len; ++k) {
          clabs[clabs_len*i+k] = ' ';
        }
        break;
      } else {
        clabs[clabs_len*i+j] = cclabs[i][j];
      }
    }
  }

  for (i = 0; i < *ncol; ++i) {
    for (j = 0; j < ctyps_len; ++j) {
      if (cctyps[i][j] == '\0') {
        for (k = j; k < ctyps_len; ++k) {
          ctyps[ctyps_len*i+k] = ' ';
        }
        break;
      } else {
        ctyps[ctyps_len*i+j] = cctyps[i][j];
      }
    }
  }
}

/* Fortran wrapper for MtzListColumn */
FORTRAN_SUBR ( LRCLID, lrclid,
	       (int *mindx, int csetid[], int *ncol),
	       (int *mindx, int csetid[], int *ncol),
	       (int *mindx, int csetid[], int *ncol))

{
 char cclabs[MCOLUMNS][31]={{0}}, cctyps[MCOLUMNS][3]={{0}};

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRCLID");)

 if (MtzCheckSubInput(*mindx,"LRCLID",1)) return;

  *ncol = MtzListInputColumn(mtzdata[*mindx-1], cclabs, cctyps, csetid);
}

/** Return the nominal cell dimensions of the MTZ structure. In
 * fact, these are the cell dimensions of the 1st crystal in the
 * MTZ structure. It is better to use the cell dimensions of the
 * correct crystal, as obtained from LRIDX.
 * @param mindx (I) MTZ file index
 * @param cell (O) Cell dimensions.
 */
FORTRAN_SUBR ( LRCELL, lrcell,
	       (int *mindx, float cell[]),
	       (int *mindx, float cell[]),
	       (int *mindx, float cell[]))

{
  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRCELL");)

 if (MtzCheckSubInput(*mindx,"LRCELL",1)) return;

  ccp4_lrcell(mtzdata[*mindx-1]->xtal[0], cell);

}

/** Return the overall resolution limits of the MTZ structure.
 * These are the widest limits over all crystals present.
 * @param mindx (I) MTZ file index
 * @param minres (O) minimum resolution
 * @param maxres (O) maximum resolution
 */
FORTRAN_SUBR ( LRRSOL, lrrsol,
	       (int *mindx, float *minres, float *maxres),
	       (int *mindx, float *minres, float *maxres),
	       (int *mindx, float *minres, float *maxres))

{
  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRRSOL");)

 if (MtzCheckSubInput(*mindx,"LRRSOL",1)) return;

  MtzResLimits(mtzdata[*mindx-1], minres, maxres);

}

/* Fortran wrapper for ccp4_lrsymi_c 
   As for LRSYMI, but with confidence flag added. */
FORTRAN_SUBR ( LRSYMI_C, lrsymi_c,
	       (int *mindx, int *nsympx, fpstr ltypex, int *nspgrx, fpstr spgrnx,
                  fpstr pgnamx, fpstr spgconf, 
                  int ltypex_len, int spgrnx_len, int pgnamx_len, int spgconf_len),
	       (int *mindx, int *nsympx, fpstr ltypex, int *nspgrx, fpstr spgrnx,
                  fpstr pgnamx, fpstr spgconf),
	       (int *mindx, int *nsympx, fpstr ltypex, int ltypex_len, int *nspgrx, 
                  fpstr spgrnx, int spgrnx_len, fpstr pgnamx, int pgnamx_len, 
                  fpstr spgconf, int spgconf_len))

{ 
  char ltypex_temp[2],spgrnx_temp[MAXSPGNAMELENGTH+1], pgnamx_temp[11], spgconf_temp[2];

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRSYMI_C");)

 if (MtzCheckSubInput(*mindx,"LRSYMI_C",1)) return;

  ccp4_lrsymi_c(mtzdata[*mindx-1], nsympx, ltypex_temp, nspgrx, spgrnx_temp, 
                pgnamx_temp, spgconf_temp);

  ccp4_CtoFString(FTN_STR(ltypex),FTN_LEN(ltypex),ltypex_temp);
  ccp4_CtoFString(FTN_STR(spgrnx),FTN_LEN(spgrnx),spgrnx_temp);
  ccp4_CtoFString(FTN_STR(pgnamx),FTN_LEN(pgnamx),pgnamx_temp);
  ccp4_CtoFString(FTN_STR(spgconf),FTN_LEN(spgconf),spgconf_temp);
}

/* Fortran wrapper for ccp4_lrsymi */
FORTRAN_SUBR ( LRSYMI, lrsymi,
	       (int *mindx, int *nsympx, fpstr ltypex, int *nspgrx, fpstr spgrnx,
                  fpstr pgnamx, int ltypex_len, int spgrnx_len, int pgnamx_len),
	       (int *mindx, int *nsympx, fpstr ltypex, int *nspgrx, fpstr spgrnx,
                  fpstr pgnamx),
	       (int *mindx, int *nsympx, fpstr ltypex, int ltypex_len, int *nspgrx, 
                  fpstr spgrnx, int spgrnx_len, fpstr pgnamx, int pgnamx_len))

{ 
 char ltypex_temp[2],spgrnx_temp[MAXSPGNAMELENGTH+1], pgnamx_temp[11];

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRSYMI");)

 if (MtzCheckSubInput(*mindx,"LRSYMI",1)) return;

  ccp4_lrsymi(mtzdata[*mindx-1], nsympx, ltypex_temp, nspgrx, spgrnx_temp, pgnamx_temp);

  ccp4_CtoFString(FTN_STR(ltypex),FTN_LEN(ltypex),ltypex_temp);
  ccp4_CtoFString(FTN_STR(spgrnx),FTN_LEN(spgrnx),spgrnx_temp);
  ccp4_CtoFString(FTN_STR(pgnamx),FTN_LEN(pgnamx),pgnamx_temp);
}

/** Get symmetry matrices from MTZ structure.
 * @param mindx (I) MTZ file index.
 * @param nsymx (O) Number of symmetry operators held in MTZ header.
 * @param rsymx (O) Symmetry operators as 4 x 4 matrices, in the order they
 *   are held in the MTZ header. Each matrix has translations in
 *   elements [3][*]. Note that a Fortran application will reverse
 *   the order of indices.
 */
FORTRAN_SUBR ( LRSYMM, lrsymm,
	       (int *mindx, int *nsymx, float rsymx[MAXSYM][4][4]),
	       (int *mindx, int *nsymx, float rsymx[MAXSYM][4][4]),
	       (int *mindx, int *nsymx, float rsymx[MAXSYM][4][4]))

{
  int i,j,k;
  float rsym[MAXSYM][4][4];

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRSYMM");)

  if (MtzCheckSubInput(*mindx,"LRSYMM",1)) return;

  ccp4_lrsymm(mtzdata[*mindx-1], nsymx, rsym);

  /* ccp4_lrsymm returns matrices in C order, so swap here */
  for (i = 0; i < *nsymx; ++i)
    for (j = 0; j < 4; ++j)
      for (k = 0; k < 4; ++k) 
        rsymx[i][j][k] = rsym[i][k][j];

  /* register this spacegroup with csymlib_f */
  ccp4spg_register_by_symops(mtzdata[*mindx-1]->mtzsymm.nsym,mtzdata[*mindx-1]->mtzsymm.sym);

}

/** Fortran wrapper to MtzParseLabin. This matches tokens from the LABIN
 * line against the program labels supplied by the program.
 * @param mindx (I) MTZ file index.
 * @param lsprgi (I) Array of program column labels. These are the column
 *  labels used by the program, as opposed to the column labels in the file.
 * @param nlprgi (I) Number of input program labels.
 * @param ntok (I) From Parser: number of tokens on line
 * @param labin_line (I) From Parser: input line
 * @param ibeg (I) From Parser: array of starting delimiters for each token
 * @param iend (I) From Parser: array of ending delimiters for each token
 */
FORTRAN_SUBR ( LKYIN, lkyin,
	       (const int *mindx, const fpstr lsprgi, const int *nlprgi, 
                   const int *ntok, const fpstr labin_line, const int ibeg[], 
                   const int iend[], int lsprgi_len, int labin_line_len),
	       (const int *mindx, const fpstr lsprgi, const int *nlprgi, 
                   const int *ntok, const fpstr labin_line, const int ibeg[], 
                   const int iend[]),
	       (const int *mindx, const fpstr lsprgi, int lsprgi_len, 
                   const int *nlprgi, const int *ntok, const fpstr labin_line, 
                   int labin_line_len, const int ibeg[], const int iend[] ))

{ int i,j;
  char *label;
  char *temp_name;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LKYIN");)

 if (*mindx <= 0 || *mindx > MFILES) {
   printf("Error: mindx out of range!\n");
   return;
 }
  label = (char *) ccp4_utils_malloc((*nlprgi)*31*sizeof(char));

 temp_name = ccp4_FtoCString(FTN_STR(labin_line), FTN_LEN(labin_line));

  for (i = 0; i < *nlprgi; ++i) {
    for (j = 0; j < lsprgi_len; ++j) {
      if (lsprgi[lsprgi_len*i+j] == ' ') {
        break;
      } else {
        label[i*31+j] = lsprgi[lsprgi_len*i+j];
      }
    }
    label[i*31+j] = '\0';
  }
 
  if (MtzParseLabin(temp_name,label,*nlprgi,user_label_in[*mindx-1]) == -1) 
    ccperror(1,"Error in label assignments in LABIN");

  free(temp_name);
  free(label);
}

/* Fortran wrapper for MtzParseLabin */
FORTRAN_SUBR ( LKYOUT, lkyout,
	       (const int *mindx, const fpstr lsprgo, const int *nlprgo, 
                   const int *ntok, const fpstr labin_line, const int ibeg[], 
                   const int iend[], int lsprgo_len, int labin_line_len),
	       (const int *mindx, const fpstr lsprgo, const int *nlprgo, 
                   const int *ntok, const fpstr labin_line, const int ibeg[], 
                   const int iend[]),
	       (const int *mindx, const fpstr lsprgo, int lsprgo_len, 
                   const int *nlprgo, const int *ntok, const fpstr labin_line, 
                   int labin_line_len, const int ibeg[], const int iend[] ))

{ int i,j;
  char *label;
  char *temp_name;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LKYOUT");)

 if (*mindx <= 0 || *mindx > MFILES) {
   printf("Error: mindx out of range!\n");
   return;
 }

 if (*nlprgo <= 0) {
   ccp4printf(1,"Warning from LKYOUT: no program output labels !\n");
   return;
 }

 if (*ntok <= 1) {
   ccp4printf(1,"Warning from LKYOUT: no arguments to LABOUT !\n");
   return;
 }
  label = (char *) ccp4_utils_malloc((*nlprgo)*31*sizeof(char));

 temp_name = ccp4_FtoCString(FTN_STR(labin_line), FTN_LEN(labin_line));

  for (i = 0; i < *nlprgo; ++i) {
    for (j = 0; j < lsprgo_len; ++j) {
      if (lsprgo[lsprgo_len*i+j] == ' ') {
        break;
      } else {
        label[i*31+j] = lsprgo[lsprgo_len*i+j];
      }
    }
    label[i*31+j] = '\0';
  }
 
  if (MtzParseLabin(temp_name,label,*nlprgo,user_label_out[*mindx-1]) == -1) 
    ccperror(1,"Error in label assignments in LABOUT");

  free(temp_name);
  free(label);
}

/* Fortran wrapper for MtzParseLabin */
FORTRAN_SUBR ( LKYSET, lkyset,
	       (const fpstr lsprgi, const int *nlprgi, 
                   fpstr lsusrj, int kpoint[],
                   const int *itok, const int *ntok, const fpstr labin_line, 
                   const int ibeg[], const int iend[], 
                   int lsprgi_len, int lsusrj_len, int labin_line_len),
	       (const fpstr lsprgi, const int *nlprgi, 
                   fpstr lsusrj, int kpoint[],
                   const int *itok, const int *ntok, const fpstr labin_line, 
                   const int ibeg[], const int iend[]),
	       (const fpstr lsprgi, int lsprgi_len, const int *nlprgi, 
                   fpstr lsusrj, int lsusrj_len, int kpoint[],
                   const int *itok, const int *ntok, const fpstr labin_line, 
                   int labin_line_len, const int ibeg[], const int iend[] ))

{ int i,j,k;
  char *label;
  char *user_lab;
  char *temp_name;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LKYSET");)

 temp_name = ccp4_FtoCString(FTN_STR(labin_line), FTN_LEN(labin_line));
  label = (char *) ccp4_utils_malloc((*nlprgi)*31*sizeof(char));
  user_lab = (char *) ccp4_utils_malloc((*nlprgi)*62*sizeof(char));

  for (i = 0; i < *nlprgi; ++i) {
    for (j = 0; j < 30; ++j) {
      if (lsprgi[lsprgi_len*i+j] == ' ') {
        break;
      } else {
        label[i*31+j] = lsprgi[lsprgi_len*i+j];
      }
    }
    label[i*31+j] = '\0';
  }
 
  if (MtzParseLabin(temp_name,label,*nlprgi,user_lab) == -1) 
    ccperror(1,"Error in label assignments in LKYSET");

  for (i = 0; i < *nlprgi; ++i) {
    /* leave kpoint unchanged unless user label exists */
    if (strcmp(user_lab+31+i*62,"") != 0)
      kpoint[i] = -1;
    for (j = 0; j < lsprgi_len; ++j) {
      if (user_lab[31+i*62+j] == '\0') {
        for (k = j; k < lsprgi_len; ++k) {
          lsusrj[lsprgi_len*i+k] = ' ';
        }
        break;
      } else {
        lsusrj[lsprgi_len*i+j] = user_lab[31+i*62+j];
      }
    }
  }

 free(temp_name);
 free(label);
 free(user_lab);
}

/** Fortran wrapper to ccp4_lrassn. First this updates labels from 
 * user_label_in if set by lkyin, then sets collookup array of pointers 
 * to columns.
 * @param mindx (I) MTZ file index.
 * @param lsprgi (I) Array of program column labels. These are the column
 *  labels used by the program, as opposed to the column labels in the file.
 * @param nlprgi (I) Number of input program labels.
 * @param lookup (I/O) On input, indicates whether a column is compulsory
 *  or not (-1 = hard compulsory - program will fail if column not found,
 *  1 = soft compulsory - program will attempt to find column even if not
 *  assigned on LABIN, 0 = optional). On output, gives the index of the
 *  column in the input file.
 * @param ctprgi (I) Array of column types.
 */
FORTRAN_SUBR ( LRASSN, lrassn,
	       (const int *mindx, fpstr lsprgi, int *nlprgi, int lookup[], fpstr ctprgi, 
                      int lsprgi_len, int ctprgi_len),
	       (const int *mindx, fpstr lsprgi, int *nlprgi, int lookup[], fpstr ctprgi),
	       (const int *mindx, fpstr lsprgi, int lsprgi_len, int *nlprgi, 
                      int lookup[], fpstr ctprgi, int ctprgi_len))

{ int i,j,k,l,err;
  char *label;
  char *type;
  MTZCOL **colarray;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRASSN");)

  err = 0;
  if (MtzCheckSubInput(*mindx,"LRASSN",1)) return;

  label = (char *) ccp4_utils_malloc((*nlprgi)*31*sizeof(char));
  type = (char *) ccp4_utils_malloc((*nlprgi)*3*sizeof(char));

 /* if lkyin has been called, used those labels */
 /* for soft and hard compulsory labels, use program labels */
  for (i = 0; i < *nlprgi; ++i) {
    if (strcmp(user_label_in[*mindx-1][i][1],"") != 0) {

      if (!MtzColLookup(mtzdata[*mindx-1],user_label_in[*mindx-1][i][1])) {
       printf("Error in LABIN: label %s not found in file!\n",user_label_in[*mindx-1][i][1]);
       err++;
      }
      strcpy(label+i*31,user_label_in[*mindx-1][i][1]);

    } else if (lookup[i] != 0) {
      for (j = 0; j < 30; ++j) {
        if (j == lsprgi_len || lsprgi[lsprgi_len*i+j] == ' ') {
          break;
        } else {
          label[i*31+j] = lsprgi[lsprgi_len*i+j];
        }
      }
      label[i*31+j] = '\0';
    } else {
      label[i*31] = '\0';
    }
  }

  /* Exit on error */
  if (err) {
    ccperror(1,"Error in label assignments in LABIN");
  }

  for (i = 0; i < *nlprgi; ++i) {
    for (j = 0; j < 2; ++j) {
      if (j == ctprgi_len || ctprgi[ctprgi_len*i+j] == ' ') {
        break;
      } else {
        type[i*3+j] = ctprgi[ctprgi_len*i+j];
      }
    }
    type[i*3+j] = '\0';
  }

  colarray = ccp4_lrassn(mtzdata[*mindx-1],label,*nlprgi,type);
  for (l = 0; l < *nlprgi; ++l) {
    collookup[*mindx-1][l] = colarray[l];
  }
  /* Set maximum number of columns to read. This is used in LRREFF and LRREFM */
  ndatmss[*mindx-1] = *nlprgi;

  for (l = 0; l < *nlprgi; ++l) {
 /* Loop over all columns */
   for (i = 0; i < mtzdata[*mindx-1]->nxtal; ++i) {
    for (j = 0; j < mtzdata[*mindx-1]->xtal[i]->nset; ++j) {
     for (k = 0; k < mtzdata[*mindx-1]->xtal[i]->set[j]->ncol; ++k) {
       if (mtzdata[*mindx-1]->xtal[i]->set[j]->col[k] == collookup[*mindx-1][l]) {
         lookup[l] = mtzdata[*mindx-1]->xtal[i]->set[j]->col[k]->source;
         goto next_label;
       } 
     }
    }
   }
   /* label not assigned */ 
   /* Have compulsory labels been found? */
   if (lookup[l] == -1) {
     printf("Error: label %s not found in file!\n",label+l*31);
     err++;
   }
   lookup[l] = 0;
   next_label:
   ;
  }

  free(colarray);
  free(label);
  free(type);

  /* Exit on error */
  if (err) {
    ccperror(1,"Error in label assignments");
  }

  return;
}

/** Fortran wrapper for ccp4_lridx.
 * Return dataset information. Note requirement to input how much
 * memory allocated in calling routine.
 * @param mindx MTZ file index
 * @param project_name
 * @param crystal_name
 * @param dataset_name
 * @param isets
 * @param datcell
 * @param datwave
 * @param ndatasets On input: space reserved for dataset information.
 *  On output: number of datasets found.
 */
FORTRAN_SUBR ( LRIDX, lridx,
	       (const int *mindx, fpstr project_name, 
                  fpstr crystal_name, fpstr dataset_name,
		  int *isets, float *datcell, float *datwave,
                  int *ndatasets, int project_name_len, 
                  int crystal_name_len, int dataset_name_len),
	       (const int *mindx, fpstr project_name, 
                  fpstr crystal_name, fpstr dataset_name,
		  int *isets, float *datcell, float *datwave,
                  int *ndatasets),
	       (const int *mindx, fpstr project_name, int project_name_len, 
	          fpstr crystal_name, int crystal_name_len,
		  fpstr dataset_name, int dataset_name_len,
		  int *isets, float *datcell, float *datwave,
                  int *ndatasets))
{
  int j,k,x,d,iset;
  char cxname[64], cdname[64], cpname[64];
  int cisets;
  float cdatcell[6], cdatwave;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRIDX");)

  if (MtzCheckSubInput(*mindx,"LRIDX",1)) return;

  /* Loop over crystals */
  for (iset = 0, x = 0; x < mtzdata[*mindx-1]->nxtal; ++x) {
   /* Loop over datasets for each crystal */
   for (d = 0; d < mtzdata[*mindx-1]->xtal[x]->nset; ++d ) {
      ccp4_lridx(mtzdata[*mindx-1], mtzdata[*mindx-1]->xtal[x]->set[d], 
		 cxname, cdname, cpname, &cisets, cdatcell, &cdatwave);

      /* check calling Fortran has assigned enough memory */
      if (iset+1 > *ndatasets) {
       printf("Warning in LRIDX. You have only reserved enough memory for %d datasets but there are more in the MTZ file. \n",*ndatasets);
       printf("Only returning partial dataset information. \n");
       *ndatasets = -1;
       return;
      }

      for (j = 0; j < project_name_len; ++j) {
       if (cpname[j] == '\0') {
        for (k = j; k < project_name_len; ++k) {
          project_name[project_name_len*iset+k] = ' ';
        }
        break;
       } else {
         project_name[project_name_len*iset+j] = cpname[j];
       }
      }

      for (j = 0; j < crystal_name_len; ++j) {
       if (cxname[j] == '\0') {
        for (k = j; k < crystal_name_len; ++k) {
          crystal_name[crystal_name_len*iset+k] = ' ';
        }
        break;
       } else {
         crystal_name[crystal_name_len*iset+j] = cxname[j];
       }
      }

      for (j = 0; j < dataset_name_len; ++j) {
       if (cdname[j] == '\0') {
        for (k = j; k < dataset_name_len; ++k) {
          dataset_name[dataset_name_len*iset+k] = ' ';
        }
        break;
       } else {
         dataset_name[dataset_name_len*iset+j] = cdname[j];
       }
      }

      isets[iset] = cisets;
      for (j = 0; j < 6; ++j) 
        datcell[iset*6+j] = cdatcell[j];
      datwave[iset] = cdatwave;

      ++iset;
   }
  }
  *ndatasets = iset;

}

/* Return cell parameters associated with the specified dataset
 * @param mindx On input: MTZ file index
 * @param iset On input: integer specifying the setid for which the
 * cell parameters are required.
 * @param mtzcell Array of 6 reals. On output: populated with the 
 * required cell parameters.
 */
FORTRAN_SUBR ( LRCELX, lrcelx,
	       (const int *mindx, const int *iset, float *mtzcell),
	       (const int *mindx, const int *iset, float *mtzcell),
	       (const int *mindx, const int *iset, float *mtzcell))
{
  MTZ *mtz;
  MTZXTAL *xtal;
  MTZSET *set;
  int i,j,k;

  /* NB This function interacts directly with the mtz data structure stored
     in memory */
  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRCELX");)

  if (MtzCheckSubInput(*mindx,"LRCELX",1)) return;

  mtz = mtzdata[*mindx-1];
  for (i = 0; i < mtz->nxtal; ++i) {
    xtal = mtz->xtal[i];
    for (j = 0; j < xtal->nset; ++j) {
      set = xtal->set[j];
      /* Check if we have found the set in question */
      if (set->setid == *iset) {
	/* Copy the cell for the crystal that the dataset belongs to */
	for (k = 0; k < 6; ++k) {
	  mtzcell[k] = xtal->cell[k];
	}
	/* Job done - escape! */
	return;
      }
    }
  }
  printf("LRCELX: error, dataset %d not found\n",*iset);
  return;
}

/* Fortran wrapper for ccp4_lridx */
FORTRAN_SUBR ( LRIDC, lridc,
	       (const int *mindx, fpstr project_name, fpstr dataset_name,
		  int *isets, float *datcell, float *datwave,
                  int *ndatasets, int project_name_len, int dataset_name_len),
	       (const int *mindx, fpstr project_name, fpstr dataset_name,
		  int *isets, float *datcell, float *datwave,
                  int *ndatasets),
	       (const int *mindx, fpstr project_name, int project_name_len,
		  fpstr dataset_name, int dataset_name_len,
		  int *isets, float *datcell, float *datwave,
                  int *ndatasets))
{
  int j,k,x,d,iset;
  char cxname[64], cdname[64], cpname[64];
  int cisets;
  float cdatcell[6], cdatwave;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRIDC");)

  if (MtzCheckSubInput(*mindx,"LRIDC",1)) return;

  /* Loop over crystals */
  for (iset = 0, x = 0; x < mtzdata[*mindx-1]->nxtal; ++x) {
   /* Loop over datasets for each crystal */
   for (d = 0; d < mtzdata[*mindx-1]->xtal[x]->nset; ++d ) {
      ccp4_lridx(mtzdata[*mindx-1], mtzdata[*mindx-1]->xtal[x]->set[d], 
		 cxname, cdname, cpname, &cisets, cdatcell, &cdatwave);

      for (j = 0; j < project_name_len; ++j) {
       if (cpname[j] == '\0') {
        for (k = j; k < project_name_len; ++k) {
          project_name[project_name_len*iset+k] = ' ';
        }
        break;
       } else {
         project_name[project_name_len*iset+j] = cpname[j];
       }
      }

      for (j = 0; j < dataset_name_len; ++j) {
       if (cdname[j] == '\0') {
        for (k = j; k < dataset_name_len; ++k) {
          dataset_name[dataset_name_len*iset+k] = ' ';
        }
        break;
       } else {
         dataset_name[dataset_name_len*iset+j] = cdname[j];
       }
      }

      isets[iset] = cisets;
      for (j = 0; j < 6; ++j) 
        datcell[iset*6+j] = cdatcell[j];
      datwave[iset] = cdatwave;

      ++iset;
   }
  }
  *ndatasets = iset;
}

/* Fortran wrapper for ccp4_lridx */
FORTRAN_SUBR ( LRID, lrid,
	       (const int *mindx, fpstr project_name, fpstr dataset_name,
		  int *isets, int *ndatasets, 
                  int project_name_len, int dataset_name_len),
	       (const int *mindx, fpstr project_name, fpstr dataset_name,
		  int *isets, int *ndatasets),
	       (const int *mindx, fpstr project_name, int project_name_len,
		  fpstr dataset_name, int dataset_name_len,
		  int *isets, int *ndatasets))
{
  int j,k,x,d,iset;
  char cxname[64], cdname[64], cpname[64];
  int cisets;
  float cdatcell[6], cdatwave;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRID");)

  if (MtzCheckSubInput(*mindx,"LRID",1)) return;

  /* Loop over crystals */
  for (iset = 0, x = 0; x < mtzdata[*mindx-1]->nxtal; ++x) {
   /* Loop over datasets for each crystal */
   for (d = 0; d < mtzdata[*mindx-1]->xtal[x]->nset; ++d ) {
      ccp4_lridx(mtzdata[*mindx-1], mtzdata[*mindx-1]->xtal[x]->set[d], 
		 cxname, cdname, cpname, &cisets, cdatcell, &cdatwave);

      for (j = 0; j < project_name_len; ++j) {
       if (cpname[j] == '\0') {
        for (k = j; k < project_name_len; ++k) {
          project_name[project_name_len*iset+k] = ' ';
        }
        break;
       } else {
         project_name[project_name_len*iset+j] = cpname[j];
       }
      }

      for (j = 0; j < dataset_name_len; ++j) {
       if (cdname[j] == '\0') {
        for (k = j; k < dataset_name_len; ++k) {
          dataset_name[dataset_name_len*iset+k] = ' ';
        }
        break;
       } else {
         dataset_name[dataset_name_len*iset+j] = cdname[j];
       }
      }

      isets[iset] = cisets;

      ++iset;
   }
  }
  *ndatasets = iset;
}

/* If we are reading from memory rather than file, this simply sets irref */
FORTRAN_SUBR ( LRSEEK, lrseek,
	       (const int *mindx, int *nrefl),
	       (const int *mindx, int *nrefl),
	       (const int *mindx, int *nrefl))
{
  int respos;

 CMTZLIB_DEBUG(puts("CMTZLIB_F: LRSEEK");)

 if (MtzCheckSubInput(*mindx,"LRSEEK",1)) return;

 /* lrrefl / lrreff will increment this */
 irref[*mindx-1] = *nrefl - 1;
 if (!cmtz_in_memory) {
   respos = (*nrefl - 1) * MtzNumSourceCol(mtzdata[*mindx-1]) + SIZE1;
   ccp4_file_seek(mtzdata[*mindx-1]->filein, respos, SEEK_SET); 
 }
}

/* Fortran wrapper for ccp4_lrrefl */
/* adata returned in file order */
FORTRAN_SUBR ( LRREFL, lrrefl,
	       (const int *mindx, float *resol, float adata[], ftn_logical *eof),
	       (const int *mindx, float *resol, float adata[], ftn_logical *eof),
	       (const int *mindx, float *resol, float adata[], ftn_logical *eof))
{
 int ieof,biomol,mindex;

 /*CMTZLIB_DEBUG(puts("CMTZLIB_F: LRREFL");)*/

 /* BIOMOL compatibility */
 if (*mindx > 1000) {
   biomol = 1;
   mindex = *mindx - 1000;
 } else {
   biomol = 0;
   mindex = *mindx;
 }

 if (MtzCheckSubInput(mindex,"LRREFL",1)) return;

 ndatmss[mindex-1] = MtzNumSourceCol(mtzdata[mindex-1]);

 ++irref[mindex-1];
 ieof = ccp4_lrrefl(mtzdata[mindex-1], resol, adata, logmss[mindex-1], irref[mindex-1]);

 if (ieof) {
   *eof = FORTRAN_LOGICAL_TRUE;
 } else {
   *eof = FORTRAN_LOGICAL_FALSE;
 }
}

/* Fortran wrapper for ccp4_lrreff */
/* adata returned in order of requested columns */
FORTRAN_SUBR ( LRREFF, lrreff,
	       (const int *mindx, float *resol, float adata[], ftn_logical *eof),
	       (const int *mindx, float *resol, float adata[], ftn_logical *eof),
	       (const int *mindx, float *resol, float adata[], ftn_logical *eof))
{
 int ieof,biomol,mindex;

  /*   CMTZLIB_DEBUG(puts("CMTZLIB_F: LRREFF");) */

 /* BIOMOL compatibility */
 if (*mindx > 1000) {
   biomol = 1;
   mindex = *mindx - 1000;
 } else {
   biomol = 0;
   mindex = *mindx;
 }

 if (MtzCheckSubInput(mindex,"LRREFF",1)) return;

 ++irref[mindex-1];
 ieof = ccp4_lrreff(mtzdata[mindex-1], resol, adata, logmss[mindex-1], 
             collookup[mindex-1], ndatmss[mindex-1], irref[mindex-1]);
 if (ieof) {
   *eof = FORTRAN_LOGICAL_TRUE;
 } else {
   *eof = FORTRAN_LOGICAL_FALSE;
 }
}

/* Fortran wrapper for ccp4_lrrefm */
FORTRAN_SUBR ( LRREFM, lrrefm,
	       (const int *mindx, ftn_logical logmiss[]),
	       (const int *mindx, ftn_logical logmiss[]),
	       (const int *mindx, ftn_logical logmiss[]))
{
 int i;

 if (MtzCheckSubInput(*mindx,"LRREFM",1)) return;

  for (i = 0; i < ndatmss[*mindx-1]; ++i) {
   if (logmss[*mindx-1][i]) {
     logmiss[i] = FORTRAN_LOGICAL_TRUE;
   } else {
     logmiss[i] = FORTRAN_LOGICAL_FALSE;
   }
  }

}

/* Check for MNFs in array of ndata reals. Results returned in logmiss. */
FORTRAN_SUBR ( MTZ_CHECK_FOR_MNF, mtz_check_for_mnf,
	       (const int *mindx, const int *ndata, float adata[], ftn_logical logmiss[]),
	       (const int *mindx, const int *ndata, float adata[], ftn_logical logmiss[]),
	       (const int *mindx, const int *ndata, float adata[], ftn_logical logmiss[]))
{
 int i;

 if (MtzCheckSubInput(*mindx,"MTZ_CHECK_FOR_MNF",1)) return;

 for (i = 0; i < *ndata; ++i) {
   if (ccp4_ismnf(mtzdata[*mindx-1], adata[i])) {
     logmiss[i] = FORTRAN_LOGICAL_TRUE;
   } else {
     logmiss[i] = FORTRAN_LOGICAL_FALSE;
   }
 }
}
/* print out header information */
FORTRAN_SUBR ( LHPRT, lhprt,
	       (const int *mindx, const int *iprint), 
	       (const int *mindx, const int *iprint), 
	       (const int *mindx, const int *iprint))
{
  CMTZLIB_DEBUG(puts("CMTZLIB_F: LHPRT");)

 if (MtzCheckSubInput(*mindx,"LHPRT",1)) return;
 
 ccp4_lhprt(mtzdata[*mindx-1], *iprint);

}

/* print out header information */
FORTRAN_SUBR ( LHPRT_ADV, lhprt_adv,
	       (const int *mindx, const int *iprint), 
	       (const int *mindx, const int *iprint), 
	       (const int *mindx, const int *iprint))
{
  CMTZLIB_DEBUG(puts("CMTZLIB_F: LHPRT_ADV");)

 if (MtzCheckSubInput(*mindx,"LHPRT_ADV",1)) return;
 
 ccp4_lhprt_adv(mtzdata[*mindx-1], *iprint);

}

/** Fortran wrapper for ccp4_lrbat. Returns the header info for the next 
 * batch from the multi-record MTZ file open on index MINDX, as
 * the two arrays RBATCH (for numbers) and CBATCH (for characters).
 * @param mindx MTZ file index
 * @param batno On return, batch number
 * @param rbatch On return, real and integer batch data.
 * @param cbatch On return, character batch data (title and axes names).
 * @param iprint =0 no printing, =1 print title only, >1 print full header.
 */
FORTRAN_SUBR ( LRBAT, lrbat,
	       (const int *mindx, int *batno, float rbatch[], fpstr cbatch, 
                        const int *iprint, int cbatch_len),
	       (const int *mindx, int *batno, float rbatch[], fpstr cbatch, 
                        const int *iprint),
	       (const int *mindx, int *batno, float rbatch[], fpstr cbatch, 
                        int cbatch_len, const int *iprint))

{
  MTZBAT *batch;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRBAT");)

 if (MtzCheckSubInput(*mindx,"LRBAT",1)) return;

 if (mtzdata[*mindx-1]->n_orig_bat <= 0) {
   printf("Error: file on mindx is not a multi-record file! \n");
   return;
 }

 if (rbat[*mindx-1] == NULL) {
   *batno = -1;
   return;
 }

 batch = rbat[*mindx-1];
 *batno = batch->num;

 ccp4_lrbat(batch, rbatch, cbatch, *iprint); 

 /* advance to next batch (may be NULL) */
 rbat[*mindx-1] = batch->next;

}

/* Print batch header information. In principle, this need not be
   connected to existing data structure */
FORTRAN_SUBR ( LBPRT, lbprt,
	       (const int *ibatch, const int *iprint, float rbatch[], fpstr cbatch, 
                        int cbatch_len),
               (const int *ibatch, const int *iprint, float rbatch[], fpstr cbatch),
               (const int *ibatch, const int *iprint, float rbatch[], fpstr cbatch, 
                        int cbatch_len))
{
  char btitle[71];

  /* For batches */
  int *intbuf = (int *) rbatch;
  float *fltbuf = rbatch + NBATCHINTEGERS;
  MTZBAT *batch;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LBPRT");)

  if (! *iprint) return;

  if (*iprint == 1) {
    strncpy(btitle,cbatch,70);
    btitle[70] = '\0';
    printf(" Batch number: \n");
    printf(" %6d    %s\n",*ibatch,btitle);
  }

  batch = MtzMallocBatch();

  strncpy(batch->title,cbatch,70); 
  strncpy(batch->gonlab[0],cbatch+70,8); 
  strncpy(batch->gonlab[1],cbatch+78,8); 
  strncpy(batch->gonlab[2],cbatch+86,8); 
  batch->num = *ibatch;
  MtzArrayToBatch(intbuf, fltbuf, batch);

  MtzPrintBatchHeader(batch);

  MtzFreeBatch(batch);

  return;
}

/* Set batch pointer to batch number batno. This is
   held in the static rbat[]. If batno == 0 then
   set rbat[] to mtz->batch */
FORTRAN_SUBR ( LRBRES, lrbres,
	       (const int *mindx, const int *batno),
	       (const int *mindx, const int *batno),
	       (const int *mindx, const int *batno))

{ int istat=-1;
  MTZBAT *batch;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRBRES");)

 if (MtzCheckSubInput(*mindx,"LRBRES",1)) return;

 if (mtzdata[*mindx-1]->n_orig_bat <= 0) {
   printf("Error: file on mindx is not a multi-record file! \n");
   return;
 }

 if (*batno == 0) {
   rbat[*mindx-1] = mtzdata[*mindx-1]->batch;
 } else {
   batch = mtzdata[*mindx-1]->batch;
   while (batch != NULL) {
     if (*batno == batch->num) {
       rbat[*mindx-1] = batch;
       istat = 0;
       break;
     }
     batch = batch->next;
   }
   if (istat == -1) 
    printf("Error: file on %d has no batch %d ! \n",*mindx,*batno);
 }

}

/* Fortran wrapper for ccp4_lrbat */
FORTRAN_SUBR ( LRBTIT, lrbtit,
	       (const int *mindx, const int *batno, fpstr tbatch, 
                        const int *iprint, int tbatch_len),
	       (const int *mindx, const int *batno, fpstr tbatch, 
                        const int *iprint),
	       (const int *mindx, const int *batno, fpstr tbatch, 
                        int tbatch_len, const int *iprint))

{ int istat=-1;
  MTZBAT *batch;
  float rbatch[NBATCHWORDS];
  char cbatch[94];

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRBTIT");)

 if (MtzCheckSubInput(*mindx,"LRBTIT",1)) return;

 if (mtzdata[*mindx-1]->n_orig_bat <= 0) {
   printf("Error: file on mindx is not a multi-record file! \n");
   return;
 }

 batch = mtzdata[*mindx-1]->batch;
 while (batch != NULL) {
   if (*batno == batch->num) {
     rbat[*mindx-1] = batch;
     istat = 0;
     ccp4_lrbat(batch, rbatch, cbatch, *iprint); 
     break;
   }
   batch = batch->next;
 }
 if (istat == -1) 
    printf("Error: file on %d has no batch %d ! \n",*mindx,*batno);

 /* advance to next batch (may be NULL) */
 rbat[*mindx-1] = batch->next;

 strncpy(tbatch,cbatch,70);

}

/* Fortran wrapper for ccp4_lrbat */
FORTRAN_SUBR ( LRBSCL, lrbscl,
	       (const int *mindx, const int *batno, float batscl[], int *nbatsc),
	       (const int *mindx, const int *batno, float batscl[], int *nbatsc),
	       (const int *mindx, const int *batno, float batscl[], int *nbatsc))

{ int istat=-1, nbscal, iprint=0;
  MTZBAT *batch;
  float rbatch[NBATCHWORDS];
  int *intbatch = (int *) rbatch;
  char cbatch[94];

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRBSCL");)

 if (MtzCheckSubInput(*mindx,"LRBSCL",1)) return;

 if (mtzdata[*mindx-1]->n_orig_bat <= 0) {
   printf("Error in lrbscl: file on mindx is not a multi-record file! \n");
   return;
 }

 batch = mtzdata[*mindx-1]->batch;
 while (batch != NULL) {
   if (*batno == batch->num) {
     rbat[*mindx-1] = batch;
     istat = 0;
     ccp4_lrbat(batch, rbatch, cbatch, iprint); 
     break;
   }
   batch = batch->next;
 }
 if (istat == -1) {
    printf("Error: file on %d has no batch %d ! \n",*mindx,*batno);
   return;
 }

 /* advance to next batch (may be NULL) */
 rbat[*mindx-1] = batch->next;

 nbscal = intbatch[16];
 if (nbscal > *nbatsc) {
   printf("From LRBSCL : %d too many batch scales in orientation block for batch %d, maximum %d \n",nbscal,*batno,*nbatsc);
   return;
 }
 *nbatsc = nbscal;
 batscl[0] = rbatch[72];
 batscl[1] = rbatch[73];
 batscl[2] = rbatch[74];
 batscl[3] = rbatch[75];

}

/* Fortran wrapper for ccp4_lrbat */
FORTRAN_SUBR ( LRBSETID, lrbsetid,
	       (const int *mindx, const int *batno, int *bsetid),
	       (const int *mindx, const int *batno, int *bsetid),
	       (const int *mindx, const int *batno, int *bsetid))

{ int istat=-1, iprint=0;
  MTZBAT *batch;
  float rbatch[NBATCHWORDS];
  int *intbatch = (int *) rbatch;
  char cbatch[94];

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LRBSETID");)

 if (MtzCheckSubInput(*mindx,"LRBSETID",1)) return;

 if (mtzdata[*mindx-1]->n_orig_bat <= 0) {
   printf("Error: file on mindx is not a multi-record file! \n");
   return;
 }

 batch = mtzdata[*mindx-1]->batch;
 while (batch != NULL) {
   if (*batno == batch->num) {
     istat = 0;
     ccp4_lrbat(batch, rbatch, cbatch, iprint); 
     break;
   }
   batch = batch->next;
 }
 if (istat == -1) 
    printf("Error: file on %d has no batch %d ! \n",*mindx,*batno);

 *bsetid = intbatch[20];

}

/* Rewind to first reflection */
FORTRAN_SUBR ( LRREWD, lrrewd,
	       (const int *mindx), 
	       (const int *mindx), 
	       (const int *mindx))

{ 

 if (MtzCheckSubInput(*mindx,"LRREWD",1)) return;
 
 irref[*mindx-1] = 0;
 /* Position file at start of reflections */
 if (!cmtz_in_memory)
   ccp4_file_seek(mtzdata[*mindx-1]->filein, SIZE1, SEEK_SET); 

}

/* Fortran wrapper for MtzHklcoeffs */
FORTRAN_SUBR ( LSTRSL, lstrsl,
	       (const int *mindx, const float *a, const float *b, const float *c,
                    const float *alpha, const float *beta, const float *gamma ), 
	       (const int *mindx, const float *a, const float *b, const float *c,
                    const float *alpha, const float *beta, const float *gamma ), 
	       (const int *mindx, const float *a, const float *b, const float *c,
                    const float *alpha, const float *beta, const float *gamma ))

{ float cell[6];

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LSTRSL");)

 cell[0] = *a;
 cell[1] = *b;
 cell[2] = *c;
 cell[3] = *alpha;
 cell[4] = *beta;
 cell[5] = *gamma;

 /* calculate hkl coefficients and store in coefhkl */
 MtzHklcoeffs(cell, coefhkl[*mindx-1]);

}

/* Fortran wrapper for MtzInd2reso */
FORTRAN_SUBR (LSTLSQ1, lstlsq1,
	       (float *reso, const int *mindx, const int *ih, const int *ik, const int *il),
               (float *reso, const int *mindx, const int *ih, const int *ik, const int *il),
               (float *reso, const int *mindx, const int *ih, const int *ik, const int *il))
{
  int in[3];

/*   CMTZLIB_DEBUG(puts("CMTZLIB_F: LSTLSQ");) */

  in[0] = *ih;
  in[1] = *ik;
  in[2] = *il;

  (*reso) = 0.25*MtzInd2reso(in, coefhkl[*mindx-1]);

  return;

}

/* In fact, already closed, so free memory */
FORTRAN_SUBR ( LRCLOS, lrclos,
	       (const int *mindx), 
	       (const int *mindx), 
	       (const int *mindx))

{ 

 if (MtzCheckSubInput(*mindx,"LRCLOS",1)) return;
 
 rlun[*mindx-1] = 0;
 if (wlun[*mindx-1] == 0) {
   MtzFree(mtzdata[*mindx-1]);
   mtzdata[*mindx-1] = NULL;
 }

}

/** Fortran wrapper to open output MTZ file. In fact, if 
 * reflection data is being held in memory, defer opening
 * until MtzPut call. But if reflections are written
 * immediately to file, need to open now.
 * @param mindx MTZ file index.
 * @param filename Output file name.
 * @param ifail (O) Returns 0 if successful, non-zero otherwise.
 */
FORTRAN_SUBR ( LWOPEN_NOEXIT, lwopen_noexit,
	       (const int *mindx, fpstr filename, int *ifail, int filename_len),
	       (const int *mindx, fpstr filename, int *ifail),
	       (const int *mindx, fpstr filename, int filename_len, int *ifail))

{ 
  char *temp_name, *fullfilename, err_str[300];
  int nxtal=1, nset[1]={1};
  int i,j,k,icol;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWOPEN_NOEXIT");)

 *ifail = 0;

 if (*mindx <= 0 || *mindx > MFILES) {
   printf("Error: mindx out of range!\n");
   *ifail = 1;
   return;
 }

 temp_name = ccp4_FtoCString(FTN_STR(filename), FTN_LEN(filename));

 if (getenv("CMTZ_IN_MEMORY")) cmtz_in_memory = 1;

 /* no corresponding read file, so set up empty structure with
    HKL_base crystal and dataset. Label set ID 0, so added ones
    start at 1. */
 if (rlun[*mindx-1] == 0) {
  mtzdata[*mindx-1] = MtzMalloc(nxtal,nset);
  strcpy(mtzdata[*mindx-1]->xtal[0]->xname,"HKL_base");
  strcpy(mtzdata[*mindx-1]->xtal[0]->pname,"HKL_base");
  mtzdata[*mindx-1]->xtal[0]->set[0]->setid = 0;
  strcpy(mtzdata[*mindx-1]->xtal[0]->set[0]->dname,"HKL_base");
  mtzdata[*mindx-1]->refs_in_memory = cmtz_in_memory;
 }

 wlun[*mindx-1] = 1;
 strcpy(fileout[*mindx-1],temp_name);
 iwref[*mindx-1] = 0;

 if (!cmtz_in_memory) {

   if (mtzdata[*mindx-1]->filein) {
     if (getenv(temp_name) != NULL) {
       fullfilename = strdup(getenv(temp_name));
     } else {
       fullfilename = strdup(temp_name);
     }
     if (!strcmp(mtzdata[*mindx-1]->filein->name,fullfilename)) {
       strcpy(err_str,"LWOPEN_NOEXIT: output file is same as open input file: ");
       strncat(err_str,fullfilename,245);
       ccperror(2,err_str);
       free(temp_name);
       free(fullfilename); 
       *ifail = 1;
       return;
     }

     free(fullfilename); 
   }

   if ( !(mtzdata[*mindx-1]->fileout = MtzOpenForWrite(temp_name)) ) {
     strcpy(err_str,"LWOPEN_NOEXIT: failed to open output file ");
     strncat(err_str,temp_name,260);
     ccperror(2,err_str);
     free(temp_name);
     *ifail = 1;
     return;
   }

   /* assign existing columns for output */
   /* if lwclab/lwassn are called with iappnd=0 then these are overwritten */
   /* needs to be here for those rare progs which don't call lwclab/lwassn */
   for (i = 0; i < mtzdata[*mindx-1]->nxtal; ++i)
    for (j = 0; j < mtzdata[*mindx-1]->xtal[i]->nset; ++j)
     for (k = 0; k < mtzdata[*mindx-1]->xtal[i]->set[j]->ncol; ++k)
       if ((icol = mtzdata[*mindx-1]->xtal[i]->set[j]->col[k]->source) != 0)
         collookup_out[*mindx-1][icol-1] = mtzdata[*mindx-1]->xtal[i]->set[j]->col[k];
 }

 free(temp_name);
}

/** Fortran wrapper to open output MTZ file. In fact, if 
 * reflection data is being held in memory, defer opening
 * until MtzPut call. But if reflections are written
 * immediately to file, need to open now.
 * @param mindx MTZ file index.
 * @param filename Output file name.
 */
FORTRAN_SUBR ( LWOPEN, lwopen,
	       (const int *mindx, fpstr filename, int filename_len),
	       (const int *mindx, fpstr filename),
	       (const int *mindx, fpstr filename, int filename_len))

{ 
  int ifail;

  FORTRAN_CALL ( LWOPEN_NOEXIT, lwopen_noexit,
	       (mindx, filename, &ifail, filename_len),
	       (mindx, filename, &ifail),
	       (mindx, filename, filename_len, &ifail));

 if ( ifail )
   ccperror(1,"LWOPEN: failed to open output file");
}

/** Set title for output file.
 * @param mindx MTZ file index.
 * @param ftitle Title to be added to output MTZ file.
 * @param flag =0 replace old title with new one, or
 *             =1 append new one to old, with one space
 */
FORTRAN_SUBR ( LWTITL, lwtitl,
	       (const int *mindx, const fpstr ftitle, const int *flag, int ftitle_len),
	       (const int *mindx, const fpstr ftitle, const int *flag),
	       (const int *mindx, const fpstr ftitle, int ftitle_len, const int *flag))

{
 char *temp_title;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWTITL");)

 if (MtzCheckSubInput(*mindx,"LWTITL",2)) return;

 temp_title = ccp4_FtoCString(FTN_STR(ftitle), FTN_LEN(ftitle));

 ccp4_lwtitl(mtzdata[*mindx-1], temp_title, *flag);
 free(temp_title);

}

/** Set sort order for output file. The integer array is stored
 * as static. Try to set sort order now, but may not be possible
 * if LWCLAB/LWASSN not yet called.
 * @param mindx MTZ file index.
 * @param sortx Sort order as integer array.
 */
FORTRAN_SUBR ( LWSORT, lwsort,
	       (const int *mindx, int sortx[5]),
	       (const int *mindx, int sortx[5]),
	       (const int *mindx, int sortx[5]))

{ int i,j,k,l,icol=0;
  MTZCOL *colsort[5];

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWSORT");)

 if (MtzCheckSubInput(*mindx,"LWSORT",2)) return;

 for (i = 0; i < 5; ++i) {
   sortorder[*mindx-1][i] = sortx[i];
   colsort[i] = NULL;
 }
 for (i = 0; i < mtzdata[*mindx-1]->nxtal; ++i) 
   for (j = 0; j < mtzdata[*mindx-1]->xtal[i]->nset; ++j) 
     for (k = 0; k < mtzdata[*mindx-1]->xtal[i]->set[j]->ncol; ++k) {
       ++icol;
       for (l = 0; l < 5; ++l)
         if (sortx[l] == icol)
           colsort[l] = mtzdata[*mindx-1]->xtal[i]->set[j]->col[k];
     }

  MtzSetSortOrder(mtzdata[*mindx-1],colsort);

}

/* Fortran wrapper for MtzAddHistory */
FORTRAN_SUBR ( LWHIST, lwhist,
	       (int *mindx, fpstr hstrng, int *nlines, int hstrng_len),
	       (int *mindx, fpstr hstrng, int *nlines),
	       (int *mindx, fpstr hstrng, int hstrng_len, int *nlines))

{
  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWHIST");)

 if (MtzCheckSubInput(*mindx,"LWHIST",2)) return;

 MtzAddHistory(mtzdata[*mindx-1], hstrng, *nlines);

}

/* Fortran wrapper for MtzAddHistory */
/* Also includes progname and date */
FORTRAN_SUBR ( LWHSTL, lwhstl,
	       (int *mindx, const fpstr hstrng, int hstrng_len),
	       (int *mindx, const fpstr hstrng),
	       (int *mindx, const fpstr hstrng, int hstrng_len))

{
 char hline[MTZRECORDLENGTH],date[11],time[9],*temp_hstrng;
 size_t len=0, Length;
 int i;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWHSTL");)

 if (*mindx <= 0 || *mindx > MFILES) {
   printf("Error in lwhstl: mindx out of range!\n");
   return;
 }

 if (!mtzdata[*mindx-1]) {
   printf("Error in lwhstl: mindx not open for write (or read)!\n");
   return;
 }

 strcpy(hline,"From ");
 if (ccp4ProgramName(NULL)) {
   len = strlen(strcpy(hline+5,ccp4ProgramName(NULL))) + 5;
   hline[len++] = ' ';
 }
 if (ccp4_utils_date(date)) {
   len = len + strlen(strcpy(hline+len,date));
   hline[len++] = ' ';
 }
 if (ccp4_utils_time(time)) {
   len = len + strlen(strcpy(hline+len,time));
   hline[len++] = ' ';
 }

/* append hstrng to hline - hline is not necessarily null-terminated */
/* MtzAddHistory assumes that hline is either null-terminated, or 
 * full - otherwise, garbage may be added 
 */
 temp_hstrng = ccp4_FtoCString(FTN_STR(hstrng), FTN_LEN(hstrng));
 Length = strlen(temp_hstrng);
 if (Length > MTZRECORDLENGTH-len) Length = MTZRECORDLENGTH-len;
 strncpy(hline+len,temp_hstrng,Length);
 for (i = len+Length; i < MTZRECORDLENGTH; ++i) hline[i] = '\0';

 MtzAddHistory(mtzdata[*mindx-1], hline, 1);

 free(temp_hstrng);
}

/** Fortran wrapper to ccp4_lwidx for writing dataset header information.
 * As for LWIDX except crystal name is not provided, and defaults to
 * supplied project name. Also cell and wavelength are not provided and
 * default to zero. This exists for backwards-compatibility - use LWIDX
 * instead.
 * @param mindx MTZ file index.
 * @param project_name Name of project that parent crystal belongs to.
 * @param dataset_name Name of dataset.
 */
FORTRAN_SUBR ( LWID, lwid,
	       (const int *mindx, const fpstr project_name, const fpstr dataset_name, 
                  int project_name_len, int dataset_name_len),
	       (const int *mindx, const fpstr project_name, const fpstr dataset_name),
	       (const int *mindx, const fpstr project_name, int project_name_len,
                  const fpstr dataset_name, int dataset_name_len))

{ char *temp_xname, *temp_pname, *temp_dname;
  float datcell[6]={0.0},datwave=0.0;
 MTZ *mtz;
 int i;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWID");)

 if (MtzCheckSubInput(*mindx,"LWID",2)) return;

 temp_pname = ccp4_FtoCString(FTN_STR(project_name), FTN_LEN(project_name));
 temp_dname = ccp4_FtoCString(FTN_STR(dataset_name), FTN_LEN(dataset_name));

 /* default crystal name to project name */
 temp_xname = strdup(temp_pname);

 mtz = mtzdata[*mindx-1];

 /* if it's a new crystal, then default cell dimensions to existing ones */
 if (MtzXtalLookup(mtz,temp_xname) == NULL && 
     mtz->nxtal > 0 && mtz->xtal[0]->cell[0] != 0.0) 
  for (i = 0; i < 6; ++i) 
    datcell[i] = mtz->xtal[0]->cell[i];

  ccp4_lwidx(mtz, temp_xname, temp_dname, temp_pname, datcell, &datwave);
  free(temp_xname); 
  free(temp_pname); 
  free(temp_dname); 
}

/** Fortran wrapper to ccp4_lwidx for writing dataset header information.
 * As for LWIDX except crystal name is not provided, and defaults to
 * supplied project name. This exists for backwards-compatibility - use LWIDX
 * instead.
 * @param mindx MTZ file index.
 * @param project_name Name of project that parent crystal belongs to.
 * @param dataset_name Name of dataset.
 * @param datcell Cell dimensions of parent crystal.
 * @param datwave Wavelength of dataset.
 */
FORTRAN_SUBR ( LWIDC, lwidc,
	       (const int *mindx, const fpstr project_name, const fpstr dataset_name,
		  float datcell[6], float *datwave,
                  int project_name_len, int dataset_name_len),
	       (const int *mindx, const fpstr project_name, const fpstr dataset_name,
		  float datcell[6], float *datwave),
	       (const int *mindx, const fpstr project_name, int project_name_len,
		  const fpstr dataset_name, int dataset_name_len,
		  float datcell[6], float *datwave))

{ char *temp_xname, *temp_pname, *temp_dname;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWIDC");)

 if (MtzCheckSubInput(*mindx,"LWIDC",2)) return;

 temp_pname = ccp4_FtoCString(FTN_STR(project_name), FTN_LEN(project_name));
 temp_dname = ccp4_FtoCString(FTN_STR(dataset_name), FTN_LEN(dataset_name));

  /* default crystal name to project name */
 temp_xname = strdup(temp_pname);

  ccp4_lwidx(mtzdata[*mindx-1], temp_xname, temp_dname, temp_pname, datcell, datwave);
  free(temp_xname); 
  free(temp_pname); 
  free(temp_dname); 
}

/** Fortran wrapper to ccp4_lwidx for writing dataset header information.
 * @param mindx MTZ file index.
 * @param project_name Name of project that parent crystal belongs to.
 * @param crystal_name Name of parent crystal.
 * @param dataset_name Name of dataset.
 * @param datcell Cell dimensions of parent crystal.
 * @param datwave Wavelength of dataset.
 */
FORTRAN_SUBR ( LWIDX, lwidx,
	       (const int *mindx, const fpstr project_name, const fpstr crystal_name,
		  const fpstr dataset_name, float datcell[6], float *datwave,
                  int project_name_len, int crystal_name_len, int dataset_name_len),
	       (const int *mindx, const fpstr project_name, const fpstr crystal_name,
		  const fpstr dataset_name, float datcell[6], float *datwave),
	       (const int *mindx, const fpstr project_name, int project_name_len,
		  const fpstr crystal_name, int crystal_name_len,
		  const fpstr dataset_name, int dataset_name_len,
		  float datcell[6], float *datwave))

{ char *temp_xname, *temp_pname, *temp_dname;

 CMTZLIB_DEBUG(puts("CMTZLIB_F: LWIDX");)

 if (MtzCheckSubInput(*mindx,"LWIDX",2)) return;

 temp_pname = ccp4_FtoCString(FTN_STR(project_name), FTN_LEN(project_name));
 temp_xname = ccp4_FtoCString(FTN_STR(crystal_name), FTN_LEN(crystal_name));
 temp_dname = ccp4_FtoCString(FTN_STR(dataset_name), FTN_LEN(dataset_name));

  ccp4_lwidx(mtzdata[*mindx-1], temp_xname, temp_dname, temp_pname, datcell, datwave);
  free(temp_xname); 
  free(temp_pname); 
  free(temp_dname); 
}

/** Fortran wrapper to update cell of output MTZ file. Overall
 * cell is obsolete - we only store crystal cell dimensions.
 * Therefore this simply writes the cell dimensions for any
 * crystal which has not yet been set. Crystal cell dimensions
 * should be set directly with lwidc.
 * @param mindx MTZ file index.
 * @param cell Output cell dimensions.
 */
FORTRAN_SUBR ( LWCELL, lwcell,
	       (const int *mindx, float cell[6]),
	       (const int *mindx, float cell[6]),
	       (const int *mindx, float cell[6]))

{
 MTZ *mtz;
 int i,j;
 char xname[200]="unknown", pname[200]="unknown", dname[200]="unknown";
 float datwave = 0.0;

 CMTZLIB_DEBUG(puts("CMTZLIB_F: LWCELL");)

 if (MtzCheckSubInput(*mindx,"LWCELL",2)) return;

 /* nothing to do if input cell is zero */
 if (cell[0] < 0.001) return;

 mtz = mtzdata[*mindx-1];
 if (mtz->nxtal == 0) {
  ccp4_lwidx(mtz, xname, dname, pname, cell, &datwave);
 } else {
  for (i = 0; i < mtz->nxtal; ++i) {
   if (mtz->xtal[i]->cell[0] == 0.0) 
     for (j = 0; j < 6; ++j) 
       mtz->xtal[i]->cell[j] = cell[j];
  }
 }

 /* calculate hkl coefficients and store in coefhkl */
 MtzHklcoeffs(cell, coefhkl[*mindx-1]);
}

/** Obsolete - use LWIDASX.
 * @param mindx MTZ file index.
 * @param nlprgo Number of output columns.
 * @param pname Array of project names.
 * @param dname Array of dataset names.
 * @param iappnd If 0 then assign all columns, if 1 then assign appended columns.
 */
FORTRAN_SUBR ( LWIDAS, lwidas,
	       (const int *mindx, int *nlprgo, fpstr pname, fpstr dname, int *iappnd,
                      int pname_len, int dname_len),
	       (const int *mindx, int *nlprgo, fpstr pname, fpstr dname, int *iappnd),
	       (const int *mindx, int *nlprgo, fpstr pname, int pname_len, 
                      fpstr dname, int dname_len, int *iappnd))
{int i,j,k,istart;
  char *project_name;
  char *crystal_name;
  char *dataset_name;
  MTZSET *baseset=NULL;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWIDAS");)

  if (MtzCheckSubInput(*mindx,"LWIDAS",2)) return;

  project_name = (char *) ccp4_utils_malloc((*nlprgo)*(pname_len+1)*sizeof(char));
  crystal_name = (char *) ccp4_utils_malloc((*nlprgo)*(pname_len+1)*sizeof(char));
  dataset_name = (char *) ccp4_utils_malloc((*nlprgo)*(dname_len+1)*sizeof(char));

  for (i = 0; i < *nlprgo; ++i) {
    for (j = 0; j < pname_len; ++j) {
      if (pname[pname_len*i+j] == ' ') {
        break;
      } else {
       project_name[i*(pname_len+1)+j] = pname[pname_len*i+j];
      }
    }
    project_name[i*(pname_len+1)+j] = '\0';
  /* default crystal name to project name */
    strcpy(crystal_name+i*(pname_len+1),project_name+i*(pname_len+1));
  }

  for (i = 0; i < *nlprgo; ++i) {
    for (j = 0; j < dname_len; ++j) {
      if (dname[dname_len*i+j] == ' ') {
        break;
      } else {
        dataset_name[i*(dname_len+1)+j] = dname[dname_len*i+j];
      }
    }
    dataset_name[i*(dname_len+1)+j] = '\0';
  }

  /* assignment request is in terms of pname/dname but data structure
     is in terms of xname/dname. We need to find appropriate xname.
     Use first crystal of correct pname/dname. If none found, above
     default is used. */
  for (i = 0; i < *nlprgo; ++i) 
   for (j = 0; j < mtzdata[*mindx-1]->nxtal; ++j) 
    if (!strcmp(mtzdata[*mindx-1]->xtal[j]->pname,
                project_name+i*(pname_len+1))) 
     for (k = 0; k < mtzdata[*mindx-1]->xtal[j]->nset; ++k) 
      if (!strcmp(mtzdata[*mindx-1]->xtal[j]->set[k]->dname,
                  dataset_name+i*(dname_len+1))) {
        strncpy(crystal_name+i*(pname_len+1),
                mtzdata[*mindx-1]->xtal[j]->xname,pname_len);
        *(crystal_name+i*(pname_len+1)+pname_len) = '\0';
      }

  /* if we are appending columns, shift collookup_out */
  istart = 0;
  if (*iappnd == 1) istart = MtzNumSourceCol(mtzdata[*mindx-1]);

  /* get base dataset if it exists */
  baseset = MtzSetLookup(mtzdata[*mindx-1],"HKL_base/HKL_base");

  for (i = 0; i < *nlprgo; ++i) {

    /* sanity check */
    if (!collookup_out[*mindx-1][i+istart]) {
      printf(" LWIDAS: severe warning - column %d does not exist. \n",i);
      continue;
    }

    if (baseset && (strcmp(collookup_out[*mindx-1][i+istart]->type,"H") == 0) ) {
      MtzAssignColumn(mtzdata[*mindx-1], collookup_out[*mindx-1][i+istart], 
          "HKL_base","HKL_base");
    } else if (strlen(crystal_name+i*(pname_len+1)) && 
               strlen(dataset_name+i*(dname_len+1))) {
      MtzAssignColumn(mtzdata[*mindx-1], collookup_out[*mindx-1][i+istart], 
          crystal_name+i*(pname_len+1), dataset_name+i*(dname_len+1));
    }
  }

  free(project_name);
  free(crystal_name);
  free(dataset_name);

}

/** Assign output columns to crystal/datasets.
 * @param mindx MTZ file index.
 * @param nlprgo Number of output columns.
 * @param xname Array of crystal names for columns.
 * @param dname Array of dataset names for columns.
 * @param iappnd If 0 then assign all columns, if 1 then assign appended columns.
 */
FORTRAN_SUBR ( LWIDASX, lwidasx,
	       (const int *mindx, int *nlprgo, fpstr xname, fpstr dname, int *iappnd,
                      int xname_len, int dname_len),
	       (const int *mindx, int *nlprgo, fpstr xname, fpstr dname, int *iappnd),
	       (const int *mindx, int *nlprgo, fpstr xname, int xname_len, 
                      fpstr dname, int dname_len, int *iappnd))
{int i,j,istart;
  char *crystal_name;
  char *dataset_name;
  MTZSET *baseset=NULL;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWIDASX");)

  if (MtzCheckSubInput(*mindx,"LWIDASX",2)) return;

  crystal_name = (char *) ccp4_utils_malloc((*nlprgo)*(xname_len+1)*sizeof(char));
  dataset_name = (char *) ccp4_utils_malloc((*nlprgo)*(dname_len+1)*sizeof(char));

  for (i = 0; i < *nlprgo; ++i) {
    for (j = 0; j < xname_len; ++j) {
      if (xname[xname_len*i+j] == ' ') {
        break;
      } else {
       crystal_name[i*(xname_len+1)+j] = xname[xname_len*i+j];
      }
    }
    crystal_name[i*(xname_len+1)+j] = '\0';
  }

  for (i = 0; i < *nlprgo; ++i) {
    for (j = 0; j < dname_len; ++j) {
      if (dname[dname_len*i+j] == ' ') {
        break;
      } else {
        dataset_name[i*(dname_len+1)+j] = dname[dname_len*i+j];
      }
    }
    dataset_name[i*(dname_len+1)+j] = '\0';
  }

  /* if we are appending columns, shift collookup_out */
  istart = 0;
  if (*iappnd == 1) istart = MtzNumSourceCol(mtzdata[*mindx-1]);

  /* get base dataset if it exists */
  baseset = MtzSetLookup(mtzdata[*mindx-1],"HKL_base/HKL_base");

  for (i = 0; i < *nlprgo; ++i) {

    /* sanity check */
    if (!collookup_out[*mindx-1][i+istart]) {
      printf(" LWIDASX: severe warning - column %d does not exist. \n",i);
      continue;
    }

    if (baseset && (strcmp(collookup_out[*mindx-1][i+istart]->type,"H") == 0) ) {
      MtzAssignColumn(mtzdata[*mindx-1], collookup_out[*mindx-1][i+istart], 
          "HKL_base","HKL_base");
    } else if (strlen(crystal_name+i*(xname_len+1)) && 
               strlen(dataset_name+i*(dname_len+1))) {
      MtzAssignColumn(mtzdata[*mindx-1], collookup_out[*mindx-1][i+istart], 
          crystal_name+i*(xname_len+1), dataset_name+i*(dname_len+1));
    }
  }

  free(crystal_name);
  free(dataset_name);
}

/** Assign output columns to crystal/datasets.
 * This is a simpler version of LWIDASX to assign all columns to one dataset 
 * (except for HKL which are assigned to base dataset). 
 * @param mindx MTZ file index.
 * @param xname Crystal name for all columns.
 * @param dname Dataset name for all columns.
 */
FORTRAN_SUBR ( LWIDALL, lwidall,
	       (const int *mindx, fpstr xname, fpstr dname,
                      int xname_len, int dname_len),
	       (const int *mindx, fpstr xname, fpstr dname),
	       (const int *mindx, fpstr xname, int xname_len, 
                      fpstr dname, int dname_len))
{ int i,j,k,l=0;
  MTZCOL **colarray;
  char *crystal_name;
  char *dataset_name;
  MTZSET *baseset=NULL;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWIDALL");)
  if (MtzCheckSubInput(*mindx,"LWIDALL",2)) return;

  crystal_name = ccp4_FtoCString(FTN_STR(xname), FTN_LEN(xname));
  dataset_name = ccp4_FtoCString(FTN_STR(dname), FTN_LEN(dname));

  /* get base dataset if it exists */
  baseset = MtzSetLookup(mtzdata[*mindx-1],"HKL_base/HKL_base");

  if (strlen(crystal_name) && strlen(dataset_name)) {

   /* need direct pointers to columns as we are going to mess with xtal/set hierarchy */
   colarray = (MTZCOL **) ccp4_utils_malloc(MtzNcol(mtzdata[*mindx-1])*sizeof(MTZCOL *));

   for (i = 0; i < mtzdata[*mindx-1]->nxtal; ++i)
    for (j = 0; j < mtzdata[*mindx-1]->xtal[i]->nset; ++j)
      for (k = 0; k < mtzdata[*mindx-1]->xtal[i]->set[j]->ncol; ++k)
        colarray[l++] = mtzdata[*mindx-1]->xtal[i]->set[j]->col[k];

   for (l = 0; l < MtzNcol(mtzdata[*mindx-1]); ++l) {
     if (baseset && (strcmp(colarray[l]->type,"H") == 0) ) {
       MtzAssignColumn(mtzdata[*mindx-1], colarray[l], "HKL_base","HKL_base");
     } else {
       MtzAssignColumn(mtzdata[*mindx-1], colarray[l], crystal_name, dataset_name);
     }
   }

   free(colarray);

  }

  free(crystal_name);
  free(dataset_name);
}

/** Write or update symmetry information for MTZ header. 
 * @param mindx (I) MTZ file index.
 * @param nsymx (I) number of symmetry operators
 * @param nsympx (I) number of primitive symmetry operators
 * @param rsymx (I) Array of symmetry operators as 4 x 4 matrices. Each matrix 
 *   is input with translations in elements [3][*] (i.e. reversed with respect
 *   to the way the Fortran application sees it). This function reverses
 *   the order before passing to ccp4_lwsymm.
 * @param ltypex (I) lattice type
 * @param nspgrx (I) spacegroup number
 * @param spgrnx (I) spacegroup name
 * @param pgnamx (I) point group name
 */
FORTRAN_SUBR ( LWSYMM, lwsymm,
	       (int *mindx, int *nsymx, int *nsympx, float rsymx[MAXSYM][4][4],
                  fpstr ltypex, int *nspgrx, fpstr spgrnx, fpstr pgnamx, 
                  int ltypex_len, int spgrnx_len, int pgnamx_len),
	       (int *mindx, int *nsymx, int *nsympx, float rsymx[MAXSYM][4][4],
                  fpstr ltypex, int *nspgrx, fpstr spgrnx, fpstr pgnamx),
	       (int *mindx, int *nsymx, int *nsympx, float rsymx[MAXSYM][4][4],
                  fpstr ltypex, int ltypex_len, int *nspgrx, fpstr spgrnx, 
		int spgrnx_len, fpstr pgnamx, int pgnamx_len))
{
 char *temp_ltypex, *temp_spgrnx, *temp_pgnamx;
 int i,ii,j,k,nsym;
 float rsym[MAXSYM][4][4];

 CMTZLIB_DEBUG(puts("CMTZLIB_F: LWSYMM");)

 if (MtzCheckSubInput(*mindx,"LWSYMM",2)) return;

 temp_ltypex = ccp4_FtoCString(FTN_STR(ltypex), FTN_LEN(ltypex));
 temp_spgrnx = ccp4_FtoCString(FTN_STR(spgrnx), FTN_LEN(spgrnx));
 temp_pgnamx = ccp4_FtoCString(FTN_STR(pgnamx), FTN_LEN(pgnamx));

 /* deal with array conventions */
 nsym = *nsymx;
 if (nsym > MAXSYM) nsym = MAXSYM;
 for (i = 0; i < nsym; ++i)
   for (j = 0; j < 4; ++j)
     for (k = 0; k < 4; ++k) 
       rsym[i][j][k] = rsymx[i][k][j];

 /* if there is a cell, check specified symmetry is consistent with it */
 if (nsym > 0)
	{
	for(ii=0;ii<mtzdata[*mindx-1]->nxtal;ii++)
		{
		if(mtzdata[*mindx-1]->xtal[ii]->cell[0] != 0.0 && 
		   strcmp(mtzdata[*mindx-1]->xtal[ii]->xname,"HKL_base"))
		if (!ccp4spg_check_symm_cell(nsym,rsym,mtzdata[*mindx-1]->xtal[ii]->cell)) 
			{
			printf(" LWSYMM: severe warning - specified symmetry not consistent with cell dimensions! \n");
			printf(" Spacegroup %s \n",temp_spgrnx);
			printf(" Cell dimensions %f %f %f %f %f %f \n",mtzdata[*mindx-1]->xtal[ii]->cell[0],
			mtzdata[*mindx-1]->xtal[ii]->cell[1],mtzdata[*mindx-1]->xtal[ii]->cell[2],
			mtzdata[*mindx-1]->xtal[ii]->cell[3],mtzdata[*mindx-1]->xtal[ii]->cell[4],
				mtzdata[*mindx-1]->xtal[ii]->cell[5]);
			ccperror(1,"Error in spacegroup or cell dimensions.");
			}
		}
	}
  ccp4_lwsymm(mtzdata[*mindx-1], nsym, *nsympx, rsym, temp_ltypex, 
            *nspgrx, temp_spgrnx, temp_pgnamx);

  free(temp_ltypex);
  free(temp_spgrnx);
  free(temp_pgnamx);

}
/** Write or update symmetry information for MTZ header. 
 * @param mindx (I) MTZ file index.
 * @param spgconf (I) one-character spacegroup confidence flag
 */
FORTRAN_SUBR ( LWSYMCONF, lwsymconf,
	       (int *mindx, fpstr spgconf, int spgconf_len),
	       (int *mindx, fpstr spgconf),
	       (int *mindx, fpstr spgconf, int spgconf_len))
{
 char *temp_spgconf;

 CMTZLIB_DEBUG(puts("CMTZLIB_F: LWSYMCONF");)

 if (MtzCheckSubInput(*mindx,"LWSYMCONF",2)) return;

 temp_spgconf = ccp4_FtoCString(FTN_STR(spgconf), FTN_LEN(spgconf));

  ccp4_lwsymconf(mtzdata[*mindx-1], temp_spgconf);

  free(temp_spgconf);
}

/** Fortran wrapper to assign columns of output MTZ file. 
 * First this updates labels from user_label_out if set by lkyout,
 * then sets collookup_out array of pointers to columns.
 * @param mindx MTZ file index.
 * @param lsprgo array of output labels
 * @param nlprgo number of output labels
 * @param ctprgo array of output column types
 * @param iappnd if = 0 replace all existing columns, else if = 1 "append" to 
 *   existing columns. Note that columns are appended to the relevant datasets
 *   and are not therefore necessarily at the end of the list of columns.
 */
FORTRAN_SUBR ( LWASSN, lwassn,
	       (const int *mindx, fpstr lsprgo, const int *nlprgo, fpstr ctprgo, int *iappnd,
                      int lsprgo_len, int ctprgo_len),
	       (const int *mindx, fpstr lsprgo, const int *nlprgo, fpstr ctprgo, int *iappnd),
	       (const int *mindx, fpstr lsprgo, int lsprgo_len, const int *nlprgo, 
                      fpstr ctprgo, int ctprgo_len, int *iappnd))
{   int i,j,istart;
    char *label;
    char *type;
    MTZCOL **colarray, *colsort[5];

    CMTZLIB_DEBUG(puts("CMTZLIB_F: LWASSN");)

    if (MtzCheckSubInput(*mindx,"LWASSN",2)) return;

    label = (char *) ccp4_utils_malloc((*nlprgo)*31*sizeof(char));
    type = (char *) ccp4_utils_malloc((*nlprgo)*3*sizeof(char));

    /* Construct label array to pass to ccp4_lwassn.
       This has undergone several iterations. Watch out for
       special cases: 
         output assignment is taken from LABIN e.g. revise
           (same prog label in lsprgi and lsprgo, and iappnd = 0)
         dm has same prog label FDM in lsprgi and lsprgo, and 
           iappnd = 1 so these are distinct columns
         sfall has "labout allin" which should be sorted out by
           sfall itself
     */
    for (i = 0; i < *nlprgo; ++i) {

      /* set label to program default */
      for (j = 0; j != 30 && j != lsprgo_len && lsprgo[lsprgo_len*i+j] != ' ' ; ++j) {
        label[i*31+j] = lsprgo[lsprgo_len*i+j];
      }
      label[i*31+j] = '\0';

      /* if there was a user output assignment (LABOUT), overwrite
         with this */
      if (strcmp(user_label_out[*mindx-1][i][1],"")) {
        strcpy(label+i*31,user_label_out[*mindx-1][i][1]);
      } else {

	/* else check for input assignment (see e.g. revise.exam) */
	for (j = 0; j < MCOLUMNS; ++j) {
	  if (strcmp(user_label_in[*mindx-1][j][1],"") && !strcmp(user_label_in[*mindx-1][j][0],label+i*31)) {
	    strcpy(label+i*31,user_label_in[*mindx-1][j][1]);
	    break;
	  }
	}

      }

      /* check for duplicated labels */
      for (j = 0; j < i; ++j)
        if ( ! strncmp(label+i*31,label+j*31,31) ) {
          printf("LWASSN: duplicate column labels in output file, columns %d and %d both have the label %s \n",j,i,label+i*31);
          ccperror(1,"Duplicate column labels in output file");
        }

    } /* end of labels loop */

  /* types loop */
  for (i = 0; i < *nlprgo; ++i) {
      for (j = 0; j < 2; ++j) {
          if (j == ctprgo_len || ctprgo[ctprgo_len*i+j] == ' ') {
              break;
          } else {
              type[i*3+j] = ctprgo[ctprgo_len*i+j];
          }
      }
      type[i*3+j] = '\0';
  } /* type loop */

  /* if we are appending columns, shift collookup_out */
    istart = 0;
    if (*iappnd == 1) istart = MtzNumSourceCol(mtzdata[*mindx-1]);
    
    /* assign new columns for output */
    colarray = ccp4_lwassn(mtzdata[*mindx-1],label,*nlprgo,type,*iappnd); 

    for (j = 0; j < 5; ++j)
        colsort[j] = NULL;
    for (i = 0; i < *nlprgo; ++i) {
        collookup_out[*mindx-1][i+istart] = colarray[i];
        /* register sort order */
        for (j = 0; j < 5; ++j)
            if (sortorder[*mindx-1][j] == (i + 1))
                colsort[j] = colarray[i];
    }

    /* If LWSORT has been called, then sortorder[][] should be set, and we make
       sure sort order is up-to-date with output columns. If sortorder[][] not
       set, we retain sort order of input file. */
    if (sortorder[*mindx-1][0]) 
      MtzSetSortOrder(mtzdata[*mindx-1],colsort);
    
    free(colarray);
    free(label);
    free(type);
    
    }

/* Fortran wrapper for ccp4_lwassn */
/* As lwassn except doesn't check user_label_out */
FORTRAN_SUBR ( LWCLAB, lwclab,
	       (const int *mindx, fpstr lsprgo, const int *nlprgo, fpstr ctprgo, int *iappnd,
                      int lsprgo_len, int ctprgo_len),
	       (const int *mindx, fpstr lsprgo, const int *nlprgo, fpstr ctprgo, int *iappnd),
	       (const int *mindx, fpstr lsprgo, int lsprgo_len, const int *nlprgo, 
                      fpstr ctprgo, int ctprgo_len, int *iappnd))
{int i,j,istart;
  char *label;
  char *type;
  MTZCOL **colarray, *colsort[5];

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWCLAB");)

 if (MtzCheckSubInput(*mindx,"LWCLAB",2)) return;

  label = (char *) ccp4_utils_malloc((*nlprgo)*31*sizeof(char));
  type = (char *) ccp4_utils_malloc((*nlprgo)*3*sizeof(char));

  for (i = 0; i < *nlprgo; ++i) {
    for (j = 0; j < 30; ++j) {
      if (j == lsprgo_len || lsprgo[lsprgo_len*i+j] == ' ') {
        break;
      } else {
       label[i*31+j] = lsprgo[lsprgo_len*i+j];
      }
    }
    label[i*31+j] = '\0';

    /* check for duplicated labels */
    for (j = 0; j < i; ++j)
      if ( ! strncmp(label+i*31,label+j*31,31) ) {
        printf("LWCLAB: duplicate column labels in output file, columns %d and %d both have the label %s \n",j,i,label+i*31);
        ccperror(1,"Duplicate column labels in output file");
      }
  }

  for (i = 0; i < *nlprgo; ++i) {
    for (j = 0; j < 3; ++j) {
      if (j == ctprgo_len || ctprgo[ctprgo_len*i+j] == ' ') {
        break;
      } else {
        type[i*3+j] = ctprgo[ctprgo_len*i+j];
      }
    }
    type[i*3+j] = '\0';
  }

  /* if we are appending columns, shift collookup_out */
  istart = 0;
  if (*iappnd == 1) istart = MtzNumSourceCol(mtzdata[*mindx-1]);

  /* assign new columns for output */
  colarray = ccp4_lwassn(mtzdata[*mindx-1],label,*nlprgo,type,*iappnd); 

  for (j = 0; j < 5; ++j)
    colsort[j] = NULL;
  for (i = 0; i < *nlprgo; ++i) {
    collookup_out[*mindx-1][i+istart] = colarray[i];
    /* register sort order */
    for (j = 0; j < 5; ++j)
      if (sortorder[*mindx-1][j] == (i + 1)) 
        colsort[j] = colarray[i];
  }

  /* If LWSORT has been called, then sortorder[][] should be set, and we make
     sure sort order is up-to-date with output columns. If sortorder[][] not
     set, we retain sort order of input file. */
  if (sortorder[*mindx-1][0]) 
    MtzSetSortOrder(mtzdata[*mindx-1],colsort);

  free(colarray);
  free(label);
  free(type);

}

/** Write batch header for batch number batno.
 * @param mindx MTZ file index
 * @param batno Serial number of batch.
 * @param rbatch Real/integer batch information.
 * @param cbatch Character batch information.
 */
FORTRAN_SUBR ( LWBAT, lwbat,
	       (const int *mindx, int *batno, float rbatch[], fpstr cbatch, 
                        int cbatch_len),
	       (const int *mindx, int *batno, float rbatch[], fpstr cbatch),
	       (const int *mindx, int *batno, float rbatch[], fpstr cbatch, 
                        int cbatch_len))

{
  MTZBAT *batch;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWBAT");)

 if (MtzCheckSubInput(*mindx,"LWBAT",2)) return;

 /* No check on mtzdata[*mindx-1]->n_orig_bat  It might be 0 if this is the
    first batch written. */

 /* batno = 0 is special flag to remove batch information */
 /* Used for example in SCALA to write merged file from unmerged input */

 if (*batno == 0) {
   mtzdata[*mindx-1]->n_orig_bat = 0;
   MtzFreeBatch(mtzdata[*mindx-1]->batch);
   mtzdata[*mindx-1]->batch = NULL;
   return;
 }

 /* add as new batch */
 batch = NULL;

 if (!ccp4_lwbat(mtzdata[*mindx-1], batch, *batno, rbatch, cbatch) )
   ccperror(1,"LWBAT: error in ccp4_lwbat, see messages above");

 /* record number of batch headers for output */
 ++nbatw[*mindx-1];
}

/** Write batch header for batch number batno. Only the batch title
 * is provided, so dummy header is written.
 * @param mindx MTZ file index
 * @param batno Serial number of batch.
 * @param tbatch Batch title.
 */
FORTRAN_SUBR ( LWBTIT, lwbtit,
	       (const int *mindx, int *batno, fpstr tbatch, int tbatch_len),
	       (const int *mindx, int *batno, fpstr tbatch),
	       (const int *mindx, int *batno, fpstr tbatch, int tbatch_len))

{
  char cbatch[95]=" ";
  int length;
  float rbatch[NBATCHWORDS]={0.0};
  int *intbuf = (int *) rbatch;
  MTZBAT *batch;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWBTIT");)

 if (MtzCheckSubInput(*mindx,"LWBTIT",2)) return;

 /* No check on mtzdata[*mindx-1]->n_orig_bat  It might be 0 if this is the
    first batch written. */

 /* batno = 0 is special flag to remove batch information */
 /* Used for example in SCALA to write merged file from unmerged input */

 if (*batno == 0) {
   mtzdata[*mindx-1]->n_orig_bat = 0;
   MtzFreeBatch(mtzdata[*mindx-1]->batch);
   mtzdata[*mindx-1]->batch = NULL;
   return;
 }

 length = (FTN_LEN(tbatch) < 70) ? FTN_LEN(tbatch) : 70 ;
 strncpy(cbatch,FTN_STR(tbatch),length);
 cbatch[length] = '\0';
 intbuf[0] = NBATCHWORDS;
 intbuf[1] = NBATCHINTEGERS;
 intbuf[2] = NBATCHREALS;

 /* add as new batch */
 batch = NULL;

 if (!ccp4_lwbat(mtzdata[*mindx-1], batch, *batno, rbatch, cbatch) )
   ccperror(1,"LWBTIT: error in ccp4_lwbat, see messages above"); 

 /* record number of batch headers for output */
 ++nbatw[*mindx-1];
}

/** Write batch header for batch number batno. New batch scales
 * are set. batno must correspond to pre-existing batch.
 * @param mindx MTZ file index
 * @param batno Serial number of batch.
 * @param batscl Array of batch scales.
 * @param nbatsc Number of batch scales.
 */
FORTRAN_SUBR ( LWBSCL, lwbscl,
	       (const int *mindx, int *batno, float batscl[], int *nbatsc),
	       (const int *mindx, int *batno, float batscl[], int *nbatsc),
	       (const int *mindx, int *batno, float batscl[], int *nbatsc))

{
  int i,istat=-1, iprint=0;
  MTZBAT *batch;
  float rbatch[NBATCHWORDS];
  int *intbatch = (int *) rbatch;
  char cbatch[94];

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWBSCL");)

 if (MtzCheckSubInput(*mindx,"LWBSCL",2)) return;

 /* No check on mtzdata[*mindx-1]->n_orig_bat  It might be 0 if this is the
    first batch written. */

 batch = mtzdata[*mindx-1]->batch;
 while (batch != NULL) {
   if (*batno == batch->num) {
     istat = 0;
     ccp4_lrbat(batch, rbatch, cbatch, iprint); 
     break;
   }
   batch = batch->next;
 }
 if (istat == -1) {
   printf("Error: file on %d has no batch %d ! \n",*mindx,*batno);
   return;
 }

 /* add as new batch */
 batch = NULL;

 intbatch[16] = *nbatsc;
 for (i = 0; i < *nbatsc; ++i) 
   rbatch[72+i] = batscl[i];

 if (!ccp4_lwbat(mtzdata[*mindx-1], batch, *batno, rbatch, cbatch) )
   ccperror(1,"LWBSCL: error in ccp4_lwbat, see messages above"); 
}

/** Obsolete. Use LWBSETIDX
 * @param mindx MTZ file index
 * @param batno Serial number of batch.
 * @param project_name Project Name
 * @param dataset_name Dataset Name
 */
FORTRAN_SUBR ( LWBSETID, lwbsetid,
	       (const int *mindx, const int *batno, const fpstr project_name, 
                  const fpstr dataset_name, 
                  int project_name_len, int dataset_name_len),
	       (const int *mindx, const int *batno, const fpstr project_name, 
                  const fpstr dataset_name),
	       (const int *mindx, const int *batno, const fpstr project_name, 
                  int project_name_len,
                  const fpstr dataset_name, int dataset_name_len))

{ char *temp_xname, *temp_pname, *temp_dname;
  int i,istat=-1;
  MTZBAT *batch;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWBSETID");)
  if (MtzCheckSubInput(*mindx,"LWBSETID",2)) return;

  /* No check on mtzdata[*mindx-1]->n_orig_bat  It might be 0 if this is the
    first batch written. */

  temp_pname = ccp4_FtoCString(FTN_STR(project_name), FTN_LEN(project_name));
  temp_dname = ccp4_FtoCString(FTN_STR(dataset_name), FTN_LEN(dataset_name));

  /* default crystal name to project name */
  temp_xname = strdup(temp_pname);

  batch = mtzdata[*mindx-1]->batch;
  /* if new batch headers have been written, skip the old ones */
  if (MtzNbat(mtzdata[*mindx-1]) > mtzdata[*mindx-1]->n_orig_bat) {
     for (i=0; i < mtzdata[*mindx-1]->n_orig_bat; ++i)
       batch = batch->next;
  }

  while (batch != NULL) {
    if (*batno == batch->num) {
     istat = 0;
     break;
    }
    batch = batch->next;
  }
  if (istat == -1) {
    printf("Error in lwbsetid: file on %d has no batch %d ! \n",*mindx,*batno);
  } else {
    ccp4_lwbsetid(mtzdata[*mindx-1], batch, temp_xname, temp_dname);
  }

  free(temp_xname); 
  free(temp_pname); 
  free(temp_dname); 
}

/** Assign a batch to a particular dataset, identified by crystal name
 * and dataset name.
 * @param mindx MTZ file index
 * @param batno Serial number of batch.
 * @param crystal_name Crystal Name
 * @param dataset_name Dataset Name
 */
FORTRAN_SUBR ( LWBSETIDX, lwbsetidx,
	       (const int *mindx, const int *batno, const fpstr crystal_name, 
                  const fpstr dataset_name, 
                  int crystal_name_len, int dataset_name_len),
	       (const int *mindx, const int *batno, const fpstr crystal_name, 
                  const fpstr dataset_name),
	       (const int *mindx, const int *batno, const fpstr crystal_name, 
                  int crystal_name_len,
                  const fpstr dataset_name, int dataset_name_len))

{ char *temp_xname, *temp_dname;
  int i,istat=-1;
  MTZBAT *batch;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWBSETIDX");)
  if (MtzCheckSubInput(*mindx,"LWBSETIDX",2)) return;

  /* No check on mtzdata[*mindx-1]->n_orig_bat  It might be 0 if this is the
    first batch written. */

  temp_xname = ccp4_FtoCString(FTN_STR(crystal_name), FTN_LEN(crystal_name));
  temp_dname = ccp4_FtoCString(FTN_STR(dataset_name), FTN_LEN(dataset_name));

  batch = mtzdata[*mindx-1]->batch;
  /* if new batch headers have been written, skip the old ones */
  if (MtzNbat(mtzdata[*mindx-1]) > mtzdata[*mindx-1]->n_orig_bat) {
     for (i=0; i < mtzdata[*mindx-1]->n_orig_bat; ++i)
       batch = batch->next;
  }

  while (batch != NULL) {
    if (*batno == batch->num) {
     istat = 0;
     break;
    }
    batch = batch->next;
  }
  if (istat == -1) {
    printf("Error in lwbsetidx: file on %d has no batch %d ! \n",*mindx,*batno);
  } else {
    ccp4_lwbsetid(mtzdata[*mindx-1], batch, temp_xname, temp_dname);
  }

  free(temp_xname); 
  free(temp_dname); 
}

/** Set whole array to MNF. The value of the MNF is taken from
 * the MTZ struct on unit mindx.
 * @param mindx MTZ file index
 * @param adata Array of reflection data to be initialised.
 * @param ncol Number of columns in the array to be initialised.
 */
FORTRAN_SUBR ( EQUAL_MAGIC, equal_magic,
	       (const int *mindx, float adata[], const int *ncol),
	       (const int *mindx, float adata[], const int *ncol),
	       (const int *mindx, float adata[], const int *ncol))
{ int i;
  union float_uint_uchar uf;

  /* CMTZLIB_DEBUG(puts("CMTZLIB_F: EQUAL_MAGIC");)*/

 if (*mindx <= 0 || *mindx > MFILES) {
   printf("Error in equal_magic: mindx out of range!\n");
   return;
 }

 if ( ! mtzdata[*mindx-1] ) {
   printf("Error in equal_magic: mindx %d not open yet!\n",*mindx);
   return;
 }

 if (strncmp (mtzdata[*mindx-1]->mnf.amnf,"NAN",3) == 0) {
   uf = ccp4_nan();
 } else {
   uf.f = mtzdata[*mindx-1]->mnf.fmnf;
 }
 for (i = 0; i < *ncol; ++i) {
   adata[i] = uf.f;
 }

}

/** Set or get MNF of file.
 * @param mindx MTZ file index
 * @param val_magic Value of MNF
 * @param setval If true, set the MNF with the value in val_magic.
 * If false, return value of MNF in val_magic. Returned as true, unless
 * there is an error.
 */
FORTRAN_SUBR ( SET_MAGIC, set_magic,
	       (const int *mindx, float *val_magic, ftn_logical *setval),
	       (const int *mindx, float *val_magic, ftn_logical *setval),
	       (const int *mindx, float *val_magic, ftn_logical *setval))
{
  union float_uint_uchar uf;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: SET_MAGIC");)

 if (*mindx <= 0 || *mindx > MFILES) {
   printf("Error in set_magic: mindx out of range!\n");
   return;
 }

 if ( ! mtzdata[*mindx-1] ) {
   printf("Warning in set_magic: mindx %d not open yet! MNF not set.\n",*mindx);
   return;
 }

 if (*setval == FORTRAN_LOGICAL_FALSE) {
   /* return current MNF */
   if (strncmp (mtzdata[*mindx-1]->mnf.amnf,"NAN",3) == 0) {
     uf = ccp4_nan();
   } else {
     uf.f = mtzdata[*mindx-1]->mnf.fmnf;
   }
   *val_magic = uf.f;

 } else {
   /* set the current MNF */
   if (ccp4_utils_isnan((union float_uint_uchar *) val_magic)) {
     sprintf(mtzdata[*mindx-1]->mnf.amnf,"NAN");
   } else {
     mtzdata[*mindx-1]->mnf.fmnf = *val_magic;
   }
 }

 *setval = FORTRAN_LOGICAL_TRUE;
}

/* Change value of MNF in an array  */
FORTRAN_SUBR ( RESET_MAGIC, reset_magic,
	       (const int *mindx, const float adata[], float bdata[], 
                 const int *ncol, const float *val_magica, const float *val_magicb),
	       (const int *mindx, const float adata[], float bdata[], 
                 const int *ncol, const float *val_magica, const float *val_magicb),
	       (const int *mindx, const float adata[], float bdata[], 
                 const int *ncol, const float *val_magica, const float *val_magicb))
{ int i;
  float val_magic;
  union float_uint_uchar uf;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: RESET_MAGIC");)

  val_magic = *val_magica;
  /* replace val_magica by file value if appropriate */
 if (*mindx > 0) {
  if (rlun[*mindx-1] > 0 || wlun[*mindx-1] > 0) {
   if (strncmp (mtzdata[*mindx-1]->mnf.amnf,"NAN",3) == 0) {
     uf = ccp4_nan();
     val_magic = uf.f;
   } else {
     val_magic = mtzdata[*mindx-1]->mnf.fmnf;
   }
  }
 }

 /* if adata[i] is a MNF replace it's value */
 for (i = 0; i < *ncol; ++i) {
   bdata[i] = adata[i];
   if (ccp4_utils_isnan((union float_uint_uchar *) &val_magic)) {
     if (ccp4_utils_isnan((union float_uint_uchar *) &adata[i])) bdata[i] = *val_magicb;
   } else {
     if (adata[i] == val_magic) bdata[i] = *val_magicb;
   }
 }

}

/** Write a one array of reflection data to output file. 
 * This is a wrapper for ccp4_lwrefl.
 * @param mindx (I) MTZ file index
 * @param adata (I) Array of reflection data to write.
 * @param ifail (O) Returns 0 if successful, non-zero otherwise.
 */
FORTRAN_SUBR ( LWREFL_NOEXIT, lwrefl_noexit,
	       (const int *mindx, const float adata[], int *ifail),
	       (const int *mindx, const float adata[], int *ifail),
	       (const int *mindx, const float adata[], int *ifail))
{
 /*  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWREFL_NOEXIT");) */

 *ifail = 0;

 if (MtzCheckSubInput(*mindx,"LWREFL_NOEXIT",2)) {
   *ifail = 1;
   return;
 }

 ++iwref[*mindx-1];
 if (!ccp4_lwrefl(mtzdata[*mindx-1],adata,collookup_out[*mindx-1],
		  MtzNumActiveCol(mtzdata[*mindx-1]),iwref[*mindx-1]) ) {
   *ifail = 1;
   return;
 }
}

/** Write a one array of reflection data to output file. 
 * This is a wrapper for ccp4_lwrefl. This routine exits upon
 * failure.
 * @param mindx (I) MTZ file index
 * @param adata (I) Array of reflection data to write.
 */
FORTRAN_SUBR ( LWREFL, lwrefl,
	       (const int *mindx, const float adata[]),
	       (const int *mindx, const float adata[]),
	       (const int *mindx, const float adata[]))
{
  int ifail;

  FORTRAN_CALL ( LWREFL_NOEXIT, lwrefl_noexit,
	       (mindx, adata, &ifail),
	       (mindx, adata, &ifail),
	       (mindx, adata, &ifail));

 if ( ifail )
   ccperror(1,"LWREFL: failed to write reflection");
}

/** Write MTZ file header and close output file. Wrapper for MtzPut.
 * @param mindx (I) MTZ file index
 * @param iprint (I) Specify whether to write output file header to log.
 * @param ifail (O) Returns 0 if successful, non-zero otherwise.
 */
FORTRAN_SUBR ( LWCLOS_NOEXIT, lwclos_noexit,
	       (const int *mindx, int *iprint, int *ifail),
	       (const int *mindx, int *iprint, int *ifail),
	       (const int *mindx, int *iprint, int *ifail))

{
  char *fullfilename;

  CMTZLIB_DEBUG(puts("CMTZLIB_F: LWCLOS_NOEXIT");)

 *ifail = 0;

 if (MtzCheckSubInput(*mindx,"LWCLOS_NOEXIT",2)) {
   *ifail = 1;
   return;
 }

 /* fix number of reflections at the number "written out" */
 mtzdata[*mindx-1]->nref = iwref[*mindx-1];
 if ( !MtzPut(mtzdata[*mindx-1],fileout[*mindx-1]) ) {
   *ifail = 1;
   return;
 }

 if (getenv(fileout[*mindx-1]) != NULL) {
   fullfilename = strdup(getenv(fileout[*mindx-1]));
 } else {
   fullfilename = strdup(fileout[*mindx-1]);
 }
 ccp4printf(1,"\n WRITTEN OUTPUT MTZ FILE \n");
 ccp4printf(1," Logical Name: %s   Filename: %s \n\n",fileout[*mindx-1],fullfilename);
 if (*iprint > 0) ccp4_lhprt(mtzdata[*mindx-1], *iprint);

 wlun[*mindx-1] = 0;
 if (rlun[*mindx-1] == 0) {
   MtzFree(mtzdata[*mindx-1]);
   mtzdata[*mindx-1] = NULL;
 }

 free(fullfilename); 
}

/** Write MTZ file header and close output file. Wrapper for MtzPut.
 * @param mindx (I) MTZ file index
 * @param iprint (I) Specify whether to write output file header to log.
 */
FORTRAN_SUBR ( LWCLOS, lwclos,
	       (const int *mindx, int *iprint),
	       (const int *mindx, int *iprint),
	       (const int *mindx, int *iprint))

{
  int ifail;

  FORTRAN_CALL ( LWCLOS_NOEXIT, lwclos_noexit,
	       (mindx, iprint, &ifail),
	       (mindx, iprint, &ifail),
	       (mindx, iprint, &ifail));

 if ( ifail )
   ccperror(1,"LWCLOS: failed to write output file");
}

/* old internal routines - obsolete! */

FORTRAN_SUBR ( RBATHD, rbathd,
	       (),(),())
{
   printf("RBATHD: Obsolete internal routine: you should not be calling this!\n");
   return;
}

FORTRAN_SUBR ( WBATHD, wbathd,
	       (),(),())
{
   printf("WBATHD: Obsolete internal routine: you should not be calling this!\n");
   return;
}

FORTRAN_SUBR ( LRHDRL, lrhdrl,
	       (),(),())
{
   printf("LRHDRL: Obsolete internal routine: you should not be calling this!\n");
   return;
}

FORTRAN_SUBR ( LABPRT, labprt,
	       (),(),())
{
   printf("LABPRT: Obsolete internal routine: you should not be calling this!\n");
   return;
}

FORTRAN_SUBR ( LBPRTH, lbprth,
	       (),(),())
{
   printf("LBPRTH: Obsolete internal routine: you should not be calling this!\n");
   return;
}

FORTRAN_SUBR ( SORTUP, sortup,
	       (),(),())
{
   printf("SORTUP: Obsolete internal routine: you should not be calling this!\n");
   return;
}

FORTRAN_SUBR ( ADDLIN, addlin,
	       (),(),())
{
   printf("ADDLIN: Obsolete internal routine: you should not be calling this!\n");
   return;
}

FORTRAN_FUN (int, NEXTLN, nextln,
	       (),(),())
{
   printf("NEXTLN: Obsolete internal function: you should not be using this!\n");
   return -1;
}

FORTRAN_SUBR ( IS_MAGIC, is_magic,
	       (const float *val_magic, const float *valtst, ftn_logical *lvalms),
	       (const float *val_magic, const float *valtst, ftn_logical *lvalms),
	       (const float *val_magic, const float *valtst, ftn_logical *lvalms))
{
   printf("IS_MAGIC: Obsolete internal routine: you should not be calling this!\n");
   return;
}
