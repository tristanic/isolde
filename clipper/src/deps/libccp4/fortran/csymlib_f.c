/*
     csymlib_f.c: Fortran API to CCP4 symmetry handling functions
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

/** @page csym_f_page Fortran API to CSYM 
 *
 *  @section csym_f_file_list File list

<ul>
<li>csymlib_f.c
</ul>
 *
 *  @section csym_f_overview Overview

This library consists of a set of wrappers to the CSYM library
giving the same API as the original symlib.f For details of the
API, see the original <a href="../symlib.html">documentation</a>.
This document covers some peculiarities of the C implementation.

 *   @section csym_f_multiple Multiple Spacegroups

The set of Fortran calls which mimic the original symlib.f assume
you are working within a single spacegroup. All calls access the
same spacegroup data structure, in analogy with the COMMON blocks
of symlib.f For cases where you wish to work with multiple
spacegroups (e.g. in the program <a href="../reindex.html">REINDEX</a>,
a different set of calls is provided (the names of which generally
start with "CCP4SPG_F_"). These identify the spacegroup of interest
via an index "sindx" (by analogy with the "mindx" of mtzlib).

 *   @section csym_f_mtz Symmetry information from MTZ files

MTZ file headers contain 2 types of symmetry records:
<dl>
<dt>SYMINF
<dd>Contains number of symmetry operators, number of primitive symmetry
operators, lattice type, spacegroup number, spacegroup name and point
group name.
<dt>SYMM
<dd>A series of records holding the symmetry operators.
</dl>
Note that the spacegroup name is likely to be ambiguous at best, with
no indication of the particular setting used. The primary source of
symmetry information is therefore taken to be the list of symmetry
operators. Note also that the order of operators is important if an
ISYM column is present.

 */
 
/** @file csymlib_f.c
 *
 *  @brief Fortran API for symmetry information.
 *
 *  @author Martyn Winn 
 */

/*#define FORTRAN_CALL_DEBUG 1*/

#if defined (FORTRAN_CALL_DEBUG)
#  define CSYMLIB_DEBUG(x) x
#else
#  define CSYMLIB_DEBUG(x)
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "../ccp4/ccp4_fortran.h"
#include "../ccp4/ccp4_general.h"
#include "../ccp4/ccp4_parser.h"
#include "../ccp4/csymlib.h"
#include "../ccp4/cmtzlib.h"
#include "../ccp4/cvecmat.h"
/* rcsid[] = "$Id$" */

#define MSPAC 4
#define MAXSYM 192
/* the two constants below are also defined in ccp4/csymlib.c, keep in sync */
#define MAXSYMOPS 20
#define MAXLENSYMOPSTR 80

static CCP4SPG *spacegroup = NULL;          /* allow more than one spacegroup ?? */
static CCP4SPG *spacegrp[MSPAC] = {NULL};   /* cf. Eugene's channel for rwbrook */

void ccp4spg_mem_tidy(void) {

  CSYMLIB_DEBUG(puts("CSYMLIB_F: ccp4spg_mem_tidy");)

  /* free any existing spacegroup */
  if ( spacegroup ) ccp4spg_free(&spacegroup);

}

FORTRAN_SUBR ( INVSYM, invsym,
               (const float a[4][4], float ai[4][4]),
               (const float a[4][4], float ai[4][4]),
               (const float a[4][4], float ai[4][4]))
{
  CSYMLIB_DEBUG(puts("CSYMLIB_F: INVSYM");)

  invert4matrix(a,ai);
}

FORTRAN_SUBR ( SYMFR3, symfr3,
               (const fpstr icol, const int *i1, int *nsym, float rot[MAXSYM][4][4],
                     int *eflag, int icol_len),
               (const fpstr icol, const int *i1, int *nsym, float rot[MAXSYM][4][4],
                     int *eflag),
               (const fpstr icol, int icol_len, const int *i1, int *nsym, 
                     float rot[MAXSYM][4][4], int *eflag))
/* symfr3   ---- Read and interpret symmetry operations

   This is the same as symfr2 except that it doesn't abort on error
   Instead the error status is returned in eflag (0=success, otherwise
   indicates an error occured).
*/
{ 
  char *temp_name;
  int i,j,k,ns;
  float tmp_rot[MAXSYMOPS][4][4];

  CSYMLIB_DEBUG(puts("CSYMLIB_F: SYMFR3");)

  /* nsym is the position to store the first symop in
     Convert from Fortran (starts at 1) to C (starts at 0) */
  *nsym = *nsym - 1;
  if (*nsym < 0) *nsym = 0;

  /* Get the input string to interpret */
  temp_name = ccp4_FtoCString(FTN_STR(icol)+(*i1-1), FTN_LEN(icol)-(*i1-1));
  /* Fetch the matrices */
  if ((ns = symfr_driver(temp_name,tmp_rot)) >= 0) {
    /* Store the matrices in Fortran ordering
       i.e. reverse of that normally used in C */
    for (i = 0; i < ns; ++i)
      for (j = 0; j < 4; ++j) 
	for (k = 0; k < 4; ++k) 
	  rot[*nsym+i][j][k] = tmp_rot[i][k][j];
    *nsym = *nsym + ns;
    *eflag = 0;
  } else {
    /* Error occured in symfr_driver - return error*/
    *eflag = 1;
  }
  /* Tidy up */
  if (temp_name) free(temp_name);
  return;
}

FORTRAN_SUBR( SYMFR2, symfr2,
	      (fpstr symchs, int *icol, int *nsym, float rot[MAXSYM][4][4], int symchs_len),
	      (fpstr symchs, int *icol, int *nsym, float rot[MAXSYM][4][4]),
	      (fpstr symchs, int symchs_len, int *icol, int *nsym, float rot[MAXSYM][4][4]))
/* symfr2   ---- Read and interpret symmetry operations

   SYMFR2 recognises the following types of input:
      real space symmetry operations,  e.g. X+1/2,Y-X,Z
      reciprocal space operations,     e.g. h,l-h,-k
      reciprocal axis vectors,         e.g. a*+c*,c*,-b*
      real space axis vectors,         e.g. a,c-a,-b

   The subroutine returns the appropriate 4x4 transformation
   matrix for each operation.  The calling program must 
   interpret the resutling matrix(ces) correctly.

   Multiple symmetry operations can be specified in a single
   input line, and must be separated by * (with spaces either
   side).

   On entry, icol is the first character to look at 
             nsym is the number of the first symmetry
	     operation to be read, and returns with the last
	     one read
*/
{
  char *temp_name;
  int i,j,k,ns;
  float tmp_rot[MAXSYMOPS][4][4];

  CSYMLIB_DEBUG(puts("CSYMLIB_F: SYMFR2");)

  /* nsym is the position to store the first symop in
     Convert from Fortran (starts at 1) to C (starts at 0) */
  *nsym = *nsym - 1;
  if (*nsym < 0) *nsym = 0;

  /* Get the input string to interpret */
  temp_name = ccp4_FtoCString(FTN_STR(symchs)+(*icol-1), FTN_LEN(symchs)-(*icol-1));
  /* Fetch the matrices */
  if ((ns = symfr_driver(temp_name,tmp_rot)) >= 0) {
    /* Store the matrices in Fortran ordering
       i.e. reverse of that normally used in C */
    for (i = 0; i < ns; ++i)
      for (j = 0; j < 4; ++j) 
	for (k = 0; k < 4; ++k) 
	  rot[*nsym+i][j][k] = tmp_rot[i][k][j];
    *nsym = *nsym + ns;
  } else {
    /* Error occured in symfr_driver - abort */
    ccperror(1," **SYMMETRY OPERATOR ERROR**");
    return;
  }
  /* Tidy up */
  if (temp_name) free(temp_name);
  return;
}

/** Fortran wrapper for mat4_to_symop.
 * @param nsm number of symmetry matrices passed.
 * @param rsm symmetry matrices.
 * @param symchs symmetry strings returned.
 * @param iprint print flag.
 */
FORTRAN_SUBR ( SYMTR3, symtr3,
               (const int *nsm, const float rsm[MAXSYM][4][4], 
                     fpstr symchs, const int *iprint, int symchs_len),
               (const int *nsm, const float rsm[MAXSYM][4][4], 
                     fpstr symchs, const int *iprint),
               (const int *nsm, const float rsm[MAXSYM][4][4], 
                     fpstr symchs, int symchs_len, const int *iprint))

{ char temp_symch[80];
 int i,j,k;
 float rsym[4][4];

  CSYMLIB_DEBUG(puts("CSYMLIB_F: SYMTR3");)

  for (i = 0; i < *nsm; ++i) {  
    /* need to transpose F to C */
    for (j = 0; j < 4; ++j) 
      for (k = 0; k < 4; ++k) 
        rsym[j][k] = rsm[i][k][j];
    mat4_to_symop(temp_symch,temp_symch+79,(const float (*)[4])rsym);
    /* mat4_to_symop fills temp_symch with spaces */
    /* ccp4_CtoFString will perform strlen(temp_symch) */
    temp_symch[79] = '\0';
    ccp4_CtoFString(FTN_STR(symchs+i*FTN_LEN(symchs)),FTN_LEN(symchs),temp_symch);

    if (*iprint) {
      printf("Symmetry %d %s \n",i+1,temp_symch);
    }
  }
}

/** Fortran wrapper for mat4_to_symop.
 * @param nsm number of symmetry matrices passed.
 * @param rsm symmetry matrices.
 * @param symchs symmetry strings returned.
 */
FORTRAN_SUBR ( SYMTR4, symtr4,
               (const int *nsm, const float rsm[MAXSYM][4][4], 
                     fpstr symchs, int symchs_len),
               (const int *nsm, const float rsm[MAXSYM][4][4], 
                     fpstr symchs),
               (const int *nsm, const float rsm[MAXSYM][4][4], 
                     fpstr symchs, int symchs_len))

{ char temp_symch[80];
 int i,j,k;
 float rsym[4][4];

  CSYMLIB_DEBUG(puts("CSYMLIB_F: SYMTR4");)

  for (i = 0; i < *nsm; ++i) {  
    /* need to transpose F to C */
    for (j = 0; j < 4; ++j) 
      for (k = 0; k < 4; ++k) 
        rsym[j][k] = rsm[i][k][j];
    mat4_to_symop(temp_symch,temp_symch+80,(const float (*)[4])rsym);
    /* mat4_to_symop will pad with spaces, but ccp4_CtoFString needs 
     * null-terminated 
     */
    temp_symch[79] = '\0';
    ccp4_CtoFString(FTN_STR(symchs+i*FTN_LEN(symchs)),FTN_LEN(symchs),temp_symch);
  }
}

FORTRAN_SUBR ( PGMDF, pgmdf,
               (int *jlass, int*jcentr, int jscrew[3]),
               (int *jlass, int*jcentr, int jscrew[3]),
               (int *jlass, int*jcentr, int jscrew[3]))
{
  static int klass, icentr, iscrew[3];

  CSYMLIB_DEBUG(puts("CSYMLIB_F: PGMDF");)

  if (*jlass==0) {
    /* need to set these variables */
    *jlass = klass;
    *jcentr = icentr;
    jscrew[0] = iscrew[0];
    jscrew[1] = iscrew[1];
    jscrew[2] = iscrew[2];
  } else {
    klass = *jlass;
    icentr = *jcentr;
    iscrew[0] = jscrew[0];
    iscrew[1] = jscrew[1];
    iscrew[2] = jscrew[2];
  }
  /* sorry, too lazy to do write statements! */
}

FORTRAN_SUBR ( PGDEFN, pgdefn,
               (fpstr nampg, int *nsymp, const int *nsym, float rsmt[192][4][4],
                const ftn_logical *lprint, int nampg_len),
               (fpstr nampg, int *nsymp, const int *nsym, float rsmt[192][4][4],
                const ftn_logical *lprint),
               (fpstr nampg, int nampg_len, int *nsymp, const int *nsym, 
                float rsmt[192][4][4], const ftn_logical *lprint))
{
  int i,j,k,l,nsym1;
  ccp4_symop *op1;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: PGDEFN");)

  /* free any existing spacegroup and start again */
  if ( spacegroup ) ccp4spg_free(&spacegroup);

  op1 = (ccp4_symop *) ccp4_utils_malloc(*nsym*sizeof(ccp4_symop));
  for (i = 0; i < *nsym; ++i) {
    for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l)
	op1[i].rot[k][l] = rsmt[i][l][k];
      /* Discard any translational component - it's not required
	 anyway when looking up the point group */
      op1[i].trn[k] = 0.0;
    }
  }

  /* Throw away symops that are duplicated once the
     translations have been removed */
  nsym1 = *nsym;
  i = 0; 
  while ( i < nsym1 ) {
    j = i + 1;
    while ( j < nsym1 ) {
      if (ccp4_symop_code( op1[i] ) == ccp4_symop_code( op1[j] )) {
	/* Duplication - overwrite this with the symop
	   at the end of the list */
	--nsym1;
	for (k = 0; k < 3; ++k) {
	  for (l = 0; l < 3; ++l) {
	    op1[j].rot[k][l] = op1[nsym1].rot[k][l];
	  }
	  /* Nb don't increment j as we need to test the 'new'
	     symop for duplication before stepping on */
	}
      } else {
	/* Look at next symop */
	++j;
      }
    }
    /* Look at next symop */
    ++i;
  }

  /* first, identify a spacegroup from supplied symops */
  spacegroup = ccp4_spgrp_reverse_lookup(nsym1,op1);
  free(op1);

  if (!spacegroup) ccperror(1,"Fatal error in PGDEFN");

  ccp4_CtoFString(FTN_STR(nampg),FTN_LEN(nampg),spacegroup->point_group);
  *nsymp = spacegroup->nsymop_prim;

}


/** Return Laue number and name for current spacegroup. 
 * @param nampg Point group name (unused in this implementation)
 * @param nlaue Laue number
 * @param launam Laue name
 */
FORTRAN_SUBR ( PGNLAU, pgnlau,
               (const fpstr nampg, int *nlaue, fpstr launam,
                int nampg_len, int launam_len),
               (const fpstr nampg, int *nlaue, fpstr launam),
               (const fpstr nampg, int nampg_len, int *nlaue, 
                fpstr launam, int launam_len))
{
  char *temp_pgname;   

  CSYMLIB_DEBUG(puts("CSYMLIB_F: PGNLAU");)

  temp_pgname = ccp4_FtoCString(FTN_STR(nampg), FTN_LEN(nampg));
  if (!spacegroup || !ccp4spg_pgname_equal(spacegroup->point_group,temp_pgname)) {
    printf("PGNLAU: No spacegroup or incorrect spacegroup loaded! \n");
    free(temp_pgname);
    return;
  }

  /* We should check we have the right spacegroup! However,
     nampg is typically in the format of the MTZ header record,
     which is different from that recorded in syminfo.lib */

  *nlaue = spacegroup->nlaue;
  ccp4_CtoFString(FTN_STR(launam),FTN_LEN(launam),spacegroup->laue_name);

  free(temp_pgname);
}

/** Return Laue number and name for a spacegroup onto index "sindx". 
 * @param sindx index of this spacegroup.
 * @param nlaue Laue number
 * @param launam Laue name
 */
FORTRAN_SUBR ( CCP4SPG_F_GET_LAUE, ccp4spg_f_get_laue,
               (const int *sindx, int *nlaue, fpstr launam, int launam_len),
               (const int *sindx, int *nlaue, fpstr launam),
               (const int *sindx, int *nlaue, fpstr launam, int launam_len))
{
  CSYMLIB_DEBUG(puts("CSYMLIB_F: CCP4SPG_F_GET_LAUE");)

  if (*sindx <= 0 || *sindx > MSPAC) {
    printf("Error in CCP4SPG_F_GET_LAUE: sindx %d out of range!\n",*sindx);
    return;
  }

  if ( ! spacegrp[*sindx-1] ) {
    printf("CCP4SPG_F_GET_LAUE: No spacegroup loaded on channel %d ! \n",*sindx);
    return;
  }

  *nlaue = spacegrp[*sindx-1]->nlaue;
  ccp4_CtoFString(FTN_STR(launam),FTN_LEN(launam),spacegrp[*sindx-1]->laue_name);

}
/** Return ranges on H K L appropriate to spacegroup.
 * @param sindx index of this spacegroup.
 * @param nlaue Laue number
 * @param launam Laue name
 */
FORTRAN_SUBR ( HKLRANGE, hklrange,
               (int *ihrng0, int *ihrng1, int *ikrng0, int *ikrng1, int *ilrng0, int *ilrng1),
               (int *ihrng0, int *ihrng1, int *ikrng0, int *ikrng1, int *ilrng0, int *ilrng1),
               (int *ihrng0, int *ihrng1, int *ikrng0, int *ikrng1, int *ilrng0, int *ilrng1))
{
  int i,j,itest;
  int test[8],max;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: HKLRANGE");)

  if (!spacegroup) {
    ccperror(2,"HKLRANGE: No spacegroup loaded yet! \n");
    return;
  }

  /* set up maximum ranges */
  *ihrng0 = - (*ihrng1);
  *ikrng0 = - (*ikrng1);
  *ilrng0 = - (*ilrng1);

  max = *ihrng1;
  if (*ikrng1 > max) max = *ikrng1;
  if (*ilrng1 > max) max = *ilrng1;
  test[0] = -max-2;
  test[1] = -max-1;
  test[2] = -max+1;
  test[3] = -1;
  test[4] = 1;
  test[5] = max-1;
  test[6] = max+1;
  test[7] = max+2;

  /* now try to cut it down by testing points */
  /* this is overkill but should be safe */
  /* update: not so simple. Didn't work for R32, see bugzilla 4149 */
  /* should be fixed now, but if further problems then consider not cutting down at all */ 
  itest = 0;
  for (i = 0; i < 8; ++i)
    for (j = 0; j < 8; ++j)
      if (ccp4spg_is_in_asu(spacegroup,*ihrng0,test[i],test[j])) itest = 1;
  if (!itest) *ihrng0 = 0;
  itest = 0;
  for (i = 0; i < 8; ++i)
    for (j = 0; j < 8; ++j)
      if (ccp4spg_is_in_asu(spacegroup,*ihrng1,test[i],test[j])) itest = 1;
  if (!itest) *ihrng1 = 0;
  itest = 0;
  for (i = 0; i < 8; ++i)
    for (j = 0; j < 8; ++j)
      if (ccp4spg_is_in_asu(spacegroup,test[i],*ikrng0,test[j])) itest = 1;
  if (!itest) *ikrng0 = 0;
  itest = 0;
  for (i = 0; i < 8; ++i)
    for (j = 0; j < 8; ++j)
      if (ccp4spg_is_in_asu(spacegroup,test[i],*ikrng1,test[j])) itest = 1;
  if (!itest) *ikrng1 = 0;
  itest = 0;
  for (i = 0; i < 8; ++i)
    for (j = 0; j < 8; ++j)
      if (ccp4spg_is_in_asu(spacegroup,test[i],test[j],*ilrng0)) itest = 1;
  if (!itest) *ilrng0 = 0;
  itest = 0;
  for (i = 0; i < 8; ++i)
    for (j = 0; j < 8; ++j)
      if (ccp4spg_is_in_asu(spacegroup,test[i],test[j],*ilrng1)) itest = 1;
  if (!itest) *ilrng1 = 0;

}

/** Return the Patterson group name and number corresponding to a spacegroup
 * identified by spacegroup name and point group name.
 * @param spgnam On input, spacegroup name.
 * @param pgname On input, point group name.
 * @param patnam On return, Patterson spacegroup name.
 * @param lpatsg On return, Patterson spacegroup number.
 */
FORTRAN_SUBR ( PATSGP, patsgp,
               (const fpstr spgnam, const fpstr pgname, fpstr patnam, int *lpatsg, 
                int spgnam_len, int pgname_len, int patnam_len),
               (const fpstr spgnam, const fpstr pgname, fpstr patnam, int *lpatsg),
               (const fpstr spgnam, int spgnam_len, const fpstr pgname, 
                int pgname_len, fpstr patnam, int patnam_len, int *lpatsg))
{
  CCP4SPG *tmp_spacegroup;
  char *temp_spgnam, *temp_pgname;   

  CSYMLIB_DEBUG(puts("CSYMLIB_F: PATSGP");)

  temp_spgnam = ccp4_FtoCString(FTN_STR(spgnam), FTN_LEN(spgnam));
  temp_pgname = ccp4_FtoCString(FTN_STR(pgname), FTN_LEN(pgname));
  if ( !spacegroup || !ccp4spg_name_equal_to_lib(spacegroup->symbol_xHM,temp_spgnam) ||
              !ccp4spg_pgname_equal(spacegroup->point_group,temp_pgname) ) {

    /* load temporary spacegroup */
    if ( ! (tmp_spacegroup = ccp4spg_load_by_ccp4_spgname(temp_spgnam)) ) {
      printf("PATSGP: failed to load spacegroup info from SYMINFO! \n");
      free(temp_spgnam);
      free(temp_pgname);
      return;
    }
    *lpatsg = tmp_spacegroup->npatt;
    ccp4_CtoFString(FTN_STR(patnam),FTN_LEN(patnam),tmp_spacegroup->patt_name);
    free(tmp_spacegroup); 

  } else {

    *lpatsg = spacegroup->npatt;
    ccp4_CtoFString(FTN_STR(patnam),FTN_LEN(patnam),spacegroup->patt_name);

  }
  free(temp_spgnam);
  free(temp_pgname);
}

/** Set spacegroup for subsequent calls to ASUPUT, ASUGET, ASUSYM and ASUPHP.
 * @param spgnam spacegroup name
 * @param numsgp spacegroup number
 * @param pgname On return, point group name
 * @param msym number of symmetry matrices passed.
 * @param rrsym symmetry matrices (preferred method of identifying spacegroup).
 * @param msymp On return, number of primitive symmetry operators
 * @param mlaue On return, number of Laue group.
 * @param lprint If true, print symmetry information.
 */
FORTRAN_SUBR ( ASUSET, asuset,
	       (fpstr spgnam, int *numsgp, fpstr pgname,
                int *msym, float rrsym[192][4][4], int *msymp,
                int *mlaue, ftn_logical *lprint, int spgnam_len, int pgname_len),
	       (fpstr spgnam, int *numsgp, fpstr pgname,
                int *msym, float rrsym[192][4][4], int *msymp,
                int *mlaue, ftn_logical *lprint),
	       (fpstr spgnam, int spgnam_len, int *numsgp, 
                fpstr pgname,int pgname_len,
                int *msym, float rrsym[192][4][4], int *msymp,
                int *mlaue, ftn_logical *lprint))
{
  int i,k,l;
  ccp4_symop *op1;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: ASUSET");)

  /* free any existing spacegroup and start again */
  if ( spacegroup ) ccp4spg_free(&spacegroup);

  op1 = (ccp4_symop *) ccp4_utils_malloc(*msym*sizeof(ccp4_symop));
  for (i = 0; i < *msym; ++i) {
    for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l)
	op1[i].rot[k][l] = rrsym[i][l][k];
      op1[i].trn[k] = rrsym[i][3][k];
    }
  }

  /* Loading by symops ensures spacegroup has desired ordering of symops.
     This is important for ASUGET which may use ISYM stored in MTZ file. */
  spacegroup = ccp4_spgrp_reverse_lookup(*msym,op1);

  /* If we fail to find match for symops, fall back on spacegroup number. */
  if (!spacegroup ) {
    if (*numsgp > 0) {
      if ( ! (spacegroup = ccp4spg_load_by_ccp4_num(*numsgp)) ) {
        printf("ASUSET: failed to load spacegroup info from SYMINFO! \n");
        ccperror(1,"Fatal error in ASUSET.");
        return;
      }
    } else {
      printf("ASUSET: no spacegroup info! \n");
      ccperror(1,"Fatal error in ASUSET.");
      return;
    }
  }

  ccp4_CtoFString(FTN_STR(pgname),FTN_LEN(pgname),spacegroup->point_group);
  *msymp = spacegroup->nsymop_prim;
  *mlaue = spacegroup->nlaue;

  if (*lprint != FORTRAN_LOGICAL_FALSE) ccp4spg_print_recip_spgrp(spacegroup);

  free(op1);
}

/** Return symmetry operators and inverses, set up by ASUSET.
 * @param rassym symmetry operators.
 * @param rinsym inverse symmetry operators.
 * @param nisym number of symmetry operators returned.
 */
FORTRAN_SUBR ( ASUSYM, asusym,
	       (float rassym[384][4][4], float rinsym[384][4][4], int *nisym),
	       (float rassym[384][4][4], float rinsym[384][4][4], int *nisym),
	       (float rassym[384][4][4], float rinsym[384][4][4], int *nisym))
{
  int i,j,k,l;
  float sgn;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: ASUSYM");)

  if (spacegroup) {
    *nisym = 0;
    for (i = 0; i < spacegroup->nsymop_prim; ++i) {
      sgn = +1.0;
      for (j = 0; j < 2; ++j) {
        for (k = 0; k < 3; ++k) {
          for (l = 0; l < 3; ++l) {
            rassym[*nisym][l][k] = sgn * spacegroup->symop[i].rot[k][l];
            rinsym[*nisym][l][k] = sgn * spacegroup->invsymop[i].rot[k][l];
          }
          rassym[*nisym][3][k] = sgn * spacegroup->symop[i].trn[k];
          rinsym[*nisym][3][k] = sgn * spacegroup->invsymop[i].trn[k];
          rassym[*nisym][k][3] = 0.0;
          rinsym[*nisym][k][3] = 0.0;
        }
        rassym[*nisym][3][3] = 1.0;
        rinsym[*nisym][3][3] = 1.0;
        ++(*nisym);
        sgn = -1.0;
      }
    }
  } else {
    ccperror(2,"ASUSYM: No spacegroup loaded yet! \n");
  }

}

/** Put reflection in asymmetric unit, as set up by ASUSET.
 * @param ihkl input indices.
 * @param jhkl output indices.
 * @param isym symmetry operation applied (ISYM number).
 */
FORTRAN_SUBR ( ASUPUT, asuput,
               (const int ihkl[3], int jhkl[3], int *isym),
               (const int ihkl[3], int jhkl[3], int *isym),
               (const int ihkl[3], int jhkl[3], int *isym))
{
  int hin,kin,lin,hout,kout,lout;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: ASUPUT");)

  hin = ihkl[0]; kin = ihkl[1]; lin = ihkl[2];

  *isym = ccp4spg_put_in_asu(spacegroup, hin, kin, lin, &hout, &kout, &lout);

  jhkl[0] = hout; jhkl[1] = kout; jhkl[2] = lout; 
}

/** Get the original indices jkhl from input indices ihkl generated
 * under symmetry operation isym.
 * @param ihkl input indices.
 * @param jhkl output indices (recovered original indices).
 * @param isym symmetry operation to be applied (ISYM number).
 */
FORTRAN_SUBR ( ASUGET, asuget,
               (const int ihkl[3], int jhkl[3], const int *isym),
               (const int ihkl[3], int jhkl[3], const int *isym),
               (const int ihkl[3], int jhkl[3], const int *isym))
{
  int hin,kin,lin,hout,kout,lout;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: ASUGET");)

  hin = ihkl[0]; kin = ihkl[1]; lin = ihkl[2];

  ccp4spg_generate_indices(spacegroup, *isym, hin, kin, lin, &hout, &kout, &lout);

  jhkl[0] = hout; jhkl[1] = kout; jhkl[2] = lout; 
}

/** Generate phase of symmetry equivalent JHKL from that of IHKL.
 * @param jhkl indices hkl generated in ASUPUT
 * @param lsym symmetry number for generating JHKL
 * @param isign 1   for I+ , -1   for I-
 * @param phasin phase for reflection IHKL
 * @param phasout phase for reflection JHKL
 */
FORTRAN_SUBR ( ASUPHP, asuphp,
               (const int jhkl[3], const int *lsym, const int *isign, 
                const float *phasin, float *phasout),
               (const int jhkl[3], const int *lsym, const int *isign, 
                const float *phasin, float *phasout),
               (const int jhkl[3], const int *lsym, const int *isign, 
                const float *phasin, float *phasout))
{
  int hin,kin,lin;
  float trans[3];

  CSYMLIB_DEBUG(puts("CSYMLIB_F: ASUPHP");)

  trans[0] = spacegroup->symop[*lsym-1].trn[0];
  trans[1] = spacegroup->symop[*lsym-1].trn[1];
  trans[2] = spacegroup->symop[*lsym-1].trn[2];

  hin = jhkl[0]; kin = jhkl[1]; lin = jhkl[2];

  *phasout = ccp4spg_phase_shift(hin, kin, lin, *phasin, trans, *isign);
}

/** Loads a spacegroup onto index "sindx". The spacegroup is
 * identified by the spacegroup name.
 * @param sindx index of this spacegroup.
 * @param namspg spacegroup name.
 */
FORTRAN_SUBR ( CCP4SPG_F_LOAD_BY_NAME, ccp4spg_f_load_by_name,
	       (const int *sindx, fpstr namspg, int namspg_len),
	       (const int *sindx, fpstr namspg),
	       (const int *sindx, fpstr namspg, int namspg_len))
{ 
  char *temp_name;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: CCP4SPG_F_LOAD_BY_NAME");)

  if (*sindx <= 0 || *sindx > MSPAC) {
    printf("Error in CCP4SPG_F_LOAD_BY_NAME: sindx %d out of range!\n",*sindx);
    return;
  }

  /* free any existing spacegroup and start again */
  if ( spacegrp[*sindx-1] ) ccp4spg_free(&spacegrp[*sindx-1]);

  temp_name = ccp4_FtoCString(FTN_STR(namspg), FTN_LEN(namspg));
  if (strlen(temp_name)) {
    spacegrp[*sindx-1] = ccp4spg_load_by_ccp4_spgname(temp_name);
  }
  free (temp_name);
}

/** Loads a spacegroup onto index "sindx". The spacegroup is
 * identified by the set of symmetry matrices.
 * @param sindx index of this spacegroup.
 * @param msym number of symmetry matrices passed.
 * @param rrsym symmetry matrices.
 */
FORTRAN_SUBR ( CCP4SPG_F_LOAD_BY_OPS, ccp4spg_f_load_by_ops,
	       (const int *sindx, int *msym, float rrsym[192][4][4]),
	       (const int *sindx, int *msym, float rrsym[192][4][4]),
	       (const int *sindx, int *msym, float rrsym[192][4][4]))
{
  int i,k,l;
  ccp4_symop *op1;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: CCP4SPG_F_LOAD_BY_OPS");)

  if (*sindx <= 0 || *sindx > MSPAC) {
    printf("Error in CCP4SPG_F_LOAD_BY_OPS: sindx %d out of range!\n",*sindx);
    return;
  }

  /* free any existing spacegroup and start again */
  if ( spacegrp[*sindx-1] ) ccp4spg_free(&spacegrp[*sindx-1]);

  op1 = (ccp4_symop *) ccp4_utils_malloc(*msym*sizeof(ccp4_symop));
  for (i = 0; i < *msym; ++i) {
    for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l)
	op1[i].rot[k][l] = rrsym[i][l][k];
      op1[i].trn[k] = rrsym[i][3][k];
    }
  }

  /* Loading by symops ensures spacegroup has desired ordering of symops.
     This is important for ASUGET which may use ISYM stored in MTZ file. */
  spacegrp[*sindx-1] = ccp4_spgrp_reverse_lookup(*msym,op1);

  if (!spacegroup ) {
    printf("CCP4SPG_F_LOAD_BY_OPS: no spacegroup info! \n");
    ccperror(1,"Fatal error in CCP4SPG_F_LOAD_BY_OPS.");
    return;
  }

  ccp4spg_print_recip_spgrp(spacegrp[*sindx-1]);

  free(op1);
}

/** Compare two sets of symmetry operators to see if they are
 * in the same order. This is important for the consistent use
 * of ISYM which encodes the operator position in the list.
 * @param msym1 number of symmetry matrices passed in first list.
 * @param rrsym1 first list of symmetry matrices.
 * @param msym2 number of symmetry matrices passed in second list.
 * @param rrsym2 second list of symmetry matrices.
 * @return 1 if operator lists are equal and in the same order, 0 otherwise
 */
FORTRAN_FUN (int, CCP4SPG_F_EQUAL_OPS_ORDER, ccp4spg_f_equal_ops_order,
	       (int *msym1, float rrsym1[192][4][4],int *msym2, float rrsym2[192][4][4]),
	       (int *msym1, float rrsym1[192][4][4],int *msym2, float rrsym2[192][4][4]),
	       (int *msym1, float rrsym1[192][4][4],int *msym2, float rrsym2[192][4][4]))
{
  int i,k,l,ret;
  ccp4_symop *op1, *op2;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: CCP4SPG_F_EQUAL_OPS_ORDER");)

  op1 = (ccp4_symop *) ccp4_utils_malloc(*msym1*sizeof(ccp4_symop));
  for (i = 0; i < *msym1; ++i) {
    for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l)
	op1[i].rot[k][l] = rrsym1[i][l][k];
      op1[i].trn[k] = rrsym1[i][3][k];
    }
  }

  op2 = (ccp4_symop *) ccp4_utils_malloc(*msym2*sizeof(ccp4_symop));
  for (i = 0; i < *msym2; ++i) {
    for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l)
	op2[i].rot[k][l] = rrsym2[i][l][k];
      op2[i].trn[k] = rrsym2[i][3][k];
    }
  }

  ret = ccp4_spgrp_equal_order(*msym1, op1, *msym2, op2);

  free(op1);
  free(op2);

  return ret;
}

/** Put reflection in asymmetric unit of spacegroup on index sindx.
 * @param sindx index of this spacegroup.
 * @param ihkl input indices.
 * @param jhkl output indices.
 * @param isym symmetry operation applied (ISYM number).
 */
FORTRAN_SUBR ( CCP4SPG_F_ASUPUT, ccp4spg_f_asuput,
               (const int *sindx, const int ihkl[3], int jhkl[3], int *isym),
               (const int *sindx, const int ihkl[3], int jhkl[3], int *isym),
               (const int *sindx, const int ihkl[3], int jhkl[3], int *isym))
{
  int hin,kin,lin,hout,kout,lout;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: CCP4SPG_F_ASUPUT");)

  if (*sindx <= 0 || *sindx > MSPAC) {
    printf("Error in CCP4SPG_F_ASUPUT: sindx %d out of range!\n",*sindx);
    return;
  }

  if ( ! spacegrp[*sindx-1] ) {
    printf("CCP4SPG_F_ASUPUT: No spacegroup loaded on channel %d ! \n",*sindx);
    return;
  }

  hin = ihkl[0]; kin = ihkl[1]; lin = ihkl[2];

  *isym = ccp4spg_put_in_asu(spacegrp[*sindx-1], hin, kin, lin, &hout, &kout, &lout);

  jhkl[0] = hout; jhkl[1] = kout; jhkl[2] = lout; 
}

/** Test whether reflection or it's Friedel mate is in asu.
 * The argument nlaue is checked against the value for the current
 * spacegroup: if it differs then spacegroup->nlaue is updated temporarily.
 * @param ihkl reflection indices.
 * @param nlaue Laue group number.
 * @return 1 if in asu, -1 if -h -k -l is in asu, 0 otherwise
 */
FORTRAN_FUN (int, INASU, inasu,
	       (const int ihkl[3], const int *nlaue),
               (const int ihkl[3], const int *nlaue),
               (const int ihkl[3], const int *nlaue))
{
  int ih, ik, il, nlaue_save = -1, retval;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: INASU");)

  if (!spacegroup) {
    ccperror(2,"INASU: No spacegroup loaded yet! \n");
    return 999;
  }

  if (spacegroup->nlaue != *nlaue) {
    /* The requested Laue number is different to that for the
       current spacegroup
       Save the current Laue code and load the data for the requested code */
    nlaue_save = spacegroup->nlaue;
    if (ccp4spg_load_laue(spacegroup,*nlaue)) {
      printf("INASU: unrecognised CCP4 Laue code\n");
      return 999;
    }
  }
  ih = ihkl[0];
  ik = ihkl[1];
  il = ihkl[2];
  retval = ccp4spg_is_in_pm_asu(spacegroup,ih,ik,il);
  if (nlaue_save > -1) {
    /* Restore previous settings */
    ccp4spg_load_laue(spacegroup,nlaue_save);
  }

  return retval;
}

/** Test whether reflection or it's Friedel mate is in the asymmetric
 * unit of the spacegroup on index "sindx".
 * @param sindx index of this spacegroup.
 * @param ihkl reflection indices.
 * @return 1 if in asu, -1 if -h -k -l is in asu, 0 otherwise
 */
FORTRAN_FUN (int, CCP4SPG_F_INASU, ccp4spg_f_inasu,
	       (const int *sindx, const int ihkl[3]),
               (const int *sindx, const int ihkl[3]),
               (const int *sindx, const int ihkl[3]))
{
  int ih, ik, il, retval;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: CCP4SPG_F_INASU");)

  if (*sindx <= 0 || *sindx > MSPAC) {
    printf("Error in CCP4SPG_F_INASU: sindx %d out of range!\n",*sindx);
    return 999;
  }

  if ( ! spacegrp[*sindx-1] ) {
    printf("CCP4SPG_F_INASU: No spacegroup loaded on channel %d ! \n",*sindx);
    return 999;
  }
  ih = ihkl[0];
  ik = ihkl[1];
  il = ihkl[2];
  retval = ccp4spg_is_in_pm_asu(spacegrp[*sindx-1],ih,ik,il);

  return retval;
}

FORTRAN_SUBR ( PRTRSM, prtrsm,
	       (const fpstr pgname, const int *nsymp, 
                const float rsymiv[192][4][4], int pgname_len),
	       (const fpstr pgname, const int *nsymp, 
                const float rsymiv[192][4][4]),
	       (const fpstr pgname, int pgname_len, const int *nsymp, 
                const float rsymiv[192][4][4]))
{

  CSYMLIB_DEBUG(puts("CSYMLIB_F: PRTRSM");)

  ccp4spg_print_recip_ops(spacegroup);

}

void ccp4spg_register_by_ccp4_num(int numspg) {

  CSYMLIB_DEBUG(puts("CSYMLIB_F: ccp4spg_register_by_ccp4_num");)

   /* free any existing spacegroup and start again */
   if ( spacegroup ) ccp4spg_free(&spacegroup);

   spacegroup = ccp4spg_load_by_ccp4_num(numspg);

   if (!spacegroup) ccperror(1,"Fatal error in ccp4spg_register_by_ccp4_num");

}

void ccp4spg_register_by_symops(int nops, float rsm[][4][4]) {

  int i,k,l;
  ccp4_symop *op1;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: ccp4spg_register_by_symops");)

  /* free any existing spacegroup and start again */
  if ( spacegroup ) ccp4spg_free(&spacegroup);

  /* identify spacegroup from supplied symops */
  op1 = (ccp4_symop *) ccp4_utils_malloc(nops*sizeof(ccp4_symop));
  for (i = 0; i < nops; ++i) {
    for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l) {
	op1[i].rot[k][l] = rsm[i][k][l];
      }
      op1[i].trn[k] = rsm[i][k][3];
    }
  }
  spacegroup = ccp4_spgrp_reverse_lookup(nops,op1);

  free(op1);

  if (!spacegroup) ccperror(1,"Fatal error in ccp4spg_register_by_symops");
}

/** Fortran wrapper for ccp4spg_load_by_* functions.
 * @param ist Obsolete parameter.
 * @param lspgrp Spacegroup number in CCP4 convention. If set on
 * entry, used to search for spacegroup. Returned value is that found.
 * @param namspg_cif Spacegroup name. If set on
 * entry, used to search for spacegroup. Returned value is the full
 * extended Hermann Mauguin symbol, with one slight alteration. Symbols
 * such as 'R 3 :H' are converted to 'H 3'. This is for backwards compatibility.
 * @param namspg_cifs On output, contains the spacegroup name without
 * any spaces.
 * @param nampg On output, the point group name.
 * @param nsymp On output, the number of primitive symmetry operators.
 * @param nsym On output, the total number of symmetry operators.
 * @param rlsymmmatrx On output, the symmetry operators.
 */
FORTRAN_SUBR ( MSYMLB3, msymlb3,
	       (const int *ist, int *lspgrp, fpstr namspg_cif,
		fpstr namspg_cifs, fpstr nampg, int *nsymp, int *nsym, 
                float rlsymmmatrx[192][4][4], int namspg_cif_len,
                int namspg_cifs_len, int nampg_len),
	       (const int *ist, int *lspgrp, fpstr namspg_cif,
		fpstr namspg_cifs, fpstr nampg, int *nsymp, int *nsym, 
                float rlsymmmatrx[192][4][4]),
	       (const int *ist, int *lspgrp, fpstr namspg_cif, int namspg_cif_len,
		fpstr namspg_cifs, int namspg_cifs_len, fpstr nampg, int nampg_len, 
                int *nsymp, int *nsym, float rlsymmmatrx[192][4][4]))
{
  int i,j,k;
  char *temp_name, *shortname=NULL, *no_colon_name=NULL;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: MSYMLB3");)

  /* search by number first
     we assume that lspgrp is in CCP4 convention */
  if (*lspgrp > 0) {

    /* free any existing spacegroup and start again */
    if ( spacegroup ) ccp4spg_free(&spacegroup);

    spacegroup = ccp4spg_load_by_ccp4_num(*lspgrp);

  } else {

    /* else try to search by name */
    temp_name = ccp4_FtoCString(FTN_STR(namspg_cif), FTN_LEN(namspg_cif));
    if (strlen(temp_name)) {

      /* free any existing spacegroup and start again */
      if ( spacegroup ) ccp4spg_free(&spacegroup);

      spacegroup = ccp4spg_load_by_ccp4_spgname(temp_name);

    }
    free (temp_name);
  }

  if (spacegroup) {
    if (spacegroup->spg_ccp4_num > 0) {
      *lspgrp = spacegroup->spg_ccp4_num;
    } else {
      *lspgrp = spacegroup->spg_num;
    }
    /* produce de-coloned version of xHM symbol */
    if (strlen(spacegroup->symbol_xHM) > 0) {
      no_colon_name = (char *) ccp4_utils_malloc((strlen(spacegroup->symbol_xHM)+1)*sizeof(char));
      strcpy(no_colon_name,spacegroup->symbol_xHM);
    } else {
      /* If no _xHM try _old. This should only happen in exceptional circumstances! */
      no_colon_name = (char *) ccp4_utils_malloc((strlen(spacegroup->symbol_old)+1)*sizeof(char));
      strcpy(no_colon_name,spacegroup->symbol_old);
    }
    ccp4spg_name_de_colon(no_colon_name);
    ccp4_CtoFString(FTN_STR(namspg_cif),FTN_LEN(namspg_cif),no_colon_name);
    if (spacegroup->symbol_old) {
     if (strlen(spacegroup->symbol_old) > 0) {
      shortname = (char *) ccp4_utils_malloc((strlen(spacegroup->symbol_old)+1)*sizeof(char));
      ccp4spg_to_shortname(shortname,spacegroup->symbol_old);
     }
    } 
    if (!shortname) {
     if (strlen(no_colon_name) > 0) {
      shortname = (char *) ccp4_utils_malloc((strlen(no_colon_name)+1)*sizeof(char));
      ccp4spg_to_shortname(shortname,no_colon_name);
     }
    }
    ccp4_CtoFString(FTN_STR(namspg_cifs),FTN_LEN(namspg_cifs),shortname);
    free(shortname);
    ccp4_CtoFString(FTN_STR(nampg),FTN_LEN(nampg),spacegroup->point_group);
    *nsymp = spacegroup->nsymop_prim;
    *nsym = spacegroup->nsymop;
    for (i = 0; i < *nsym; ++i) {
      for (j = 0; j < 3; ++j) {
        for (k = 0; k < 3; ++k) 
          rlsymmmatrx[i][k][j] = spacegroup->symop[i].rot[j][k];
        rlsymmmatrx[i][3][j] = spacegroup->symop[i].trn[j];
        rlsymmmatrx[i][j][3] = 0.0;
      }
      rlsymmmatrx[i][3][3] = 1.0;
    }
  }
}

FORTRAN_SUBR ( MSYMLB, msymlb,
	       (const int *ist, int *lspgrp, fpstr namspg_cif,
		fpstr nampg, int *nsymp, int *nsym, 
                float rlsymmmatrx[192][4][4], int namspg_cif_len,
                int nampg_len),
	       (const int *ist, int *lspgrp, fpstr namspg_cif,
		fpstr nampg, int *nsymp, int *nsym, 
                float rlsymmmatrx[192][4][4]),
	       (const int *ist, int *lspgrp, fpstr namspg_cif, int namspg_cif_len,
		fpstr nampg, int nampg_len, 
                int *nsymp, int *nsym, float rlsymmmatrx[192][4][4]))
{
  char namspg_cifs;
  int namspg_cifs_len=0;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: MSYMLB");)

  FORTRAN_CALL ( MSYMLB3, msymlb3,
	       (ist, lspgrp, namspg_cif, &namspg_cifs, nampg, nsymp, nsym, 
                rlsymmmatrx, namspg_cif_len, namspg_cifs_len, nampg_len),
	       (ist, lspgrp, namspg_cif, &namspg_cifs, nampg, nsymp, nsym, 
                rlsymmmatrx),
	       (ist, lspgrp, namspg_cif, namspg_cif_len, &namspg_cifs, 
                namspg_cifs_len, nampg, nampg_len, nsymp, nsym, rlsymmmatrx));

}

FORTRAN_SUBR ( MSYMLB2, msymlb2,
	       (const int *ist, int *lspgrp, fpstr namspg_cif,
		fpstr nampg, int *nsymp, int *nsym, 
                float rlsymmmatrx[192][4][4], int namspg_cif_len,
                int nampg_len),
	       (const int *ist, int *lspgrp, fpstr namspg_cif,
		fpstr nampg, int *nsymp, int *nsym, 
                float rlsymmmatrx[192][4][4]),
	       (const int *ist, int *lspgrp, fpstr namspg_cif, int namspg_cif_len,
		fpstr nampg, int nampg_len, 
                int *nsymp, int *nsym, float rlsymmmatrx[192][4][4]))
{
  char namspg_cifs;
  int namspg_cifs_len=0;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: MSYMLB2");)

  FORTRAN_CALL ( MSYMLB3, msymlb3,
	       (ist, lspgrp, namspg_cif, &namspg_cifs, nampg, nsymp, nsym, 
                rlsymmmatrx, namspg_cif_len, namspg_cifs_len, nampg_len),
	       (ist, lspgrp, namspg_cif, &namspg_cifs, nampg, nsymp, nsym, 
                rlsymmmatrx),
	       (ist, lspgrp, namspg_cif, namspg_cif_len, &namspg_cifs, 
                namspg_cifs_len, nampg, nampg_len, nsymp, nsym, rlsymmmatrx));

}

FORTRAN_SUBR ( MSYGET, msyget,
	       (const int *ist, int *lspgrp, int *nsym, 
                float rlsymmmatrx[192][4][4]),
	       (const int *ist, int *lspgrp, int *nsym, 
                float rlsymmmatrx[192][4][4]),
	       (const int *ist, int *lspgrp, int *nsym, 
                float rlsymmmatrx[192][4][4]))
{
  char namspg_cif, namspg_cifs, nampg;
  int namspg_cif_len=0, namspg_cifs_len=0, nampg_len=0, nsymp=0;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: MSYGET");)

  FORTRAN_CALL ( MSYMLB3, msymlb3,
	       (ist, lspgrp, &namspg_cif, &namspg_cifs, &nampg, &nsymp, nsym, 
                rlsymmmatrx, namspg_cif_len, namspg_cifs_len, nampg_len),
	       (ist, lspgrp, &namspg_cif, &namspg_cifs, &nampg, &nsymp, nsym, 
                rlsymmmatrx),
	       (ist, lspgrp, &namspg_cif, namspg_cif_len, &namspg_cifs, 
                namspg_cifs_len, &nampg, nampg_len, &nsymp, nsym, rlsymmmatrx));

}

/** Epsilon zones currently set up in ccp4spg_load_spacegroup
 * If these are not available, use lookup by symops.
 * @param nsm number of symmetry operators.
 * @param nsmp number of primitive symmetry operators.
 * @param rsm symmetry matrices.
 * @param iprint If iprint > 0 then a summary of epsilon zones is printed.
 */
FORTRAN_SUBR ( EPSLN, epsln,
	       (const int *nsm, const int *nsmp, const float rsm[192][4][4],
		const int *iprint),
	       (const int *nsm, const int *nsmp, const float rsm[192][4][4],
		const int *iprint),
	       (const int *nsm, const int *nsmp, const float rsm[192][4][4],
		const int *iprint))
{
  int i,k,l;
  ccp4_symop *op1;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: EPSLN");)

  /* identify spacegroup from supplied symops */
  op1 = (ccp4_symop *) ccp4_utils_malloc(*nsm*sizeof(ccp4_symop));
  for (i = 0; i < *nsm; ++i) {
    for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l) {
	op1[i].rot[k][l] = rsm[i][l][k];
      }
      op1[i].trn[k] = rsm[i][3][k];
    }
  }
  spacegroup = ccp4_spgrp_reverse_lookup(*nsm,op1);

  if (spacegroup && *iprint > 0) ccp4spg_print_epsilon_zones(spacegroup);

  free(op1);
}

FORTRAN_SUBR ( EPSLON, epslon,
	       (const int ih[3], float *epsi, int *isysab),
	       (const int ih[3], float *epsi, int *isysab),
	       (const int ih[3], float *epsi, int *isysab))
{
  int h,k,l;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: EPSLON");)

  if (!spacegroup) {
    ccperror(2,"EPSLON: No spacegroup loaded yet! \n");
    return;
  }

  h = ih[0]; k = ih[1]; l = ih[2]; 

  *epsi = (float) ccp4spg_get_multiplicity(spacegroup, h, k, l);
  *isysab = ccp4spg_is_sysabs(spacegroup, h, k, l);
}

FORTRAN_SUBR ( CCP4SPG_F_EPSLON, ccp4spg_f_epslon,
	       (const int *sindx, const int ih[3], float *epsi, int *isysab),
	       (const int *sindx, const int ih[3], float *epsi, int *isysab),
	       (const int *sindx, const int ih[3], float *epsi, int *isysab))
{
  int h,k,l;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: CCP4SPG_F_EPSLON");)

  if (*sindx <= 0 || *sindx > MSPAC) {
    printf("Error in CCP4SPG_F_EPSLON: sindx %d out of range!\n",*sindx);
    return;
  }

  if ( ! spacegrp[*sindx-1] ) {
    printf("CCP4SPG_F_EPSLON: No spacegroup loaded on channel %d ! \n",*sindx);
    return;
  }

  h = ih[0]; k = ih[1]; l = ih[2]; 

  *epsi = (float) ccp4spg_get_multiplicity(spacegrp[*sindx-1], h, k, l);
  *isysab = ccp4spg_is_sysabs(spacegrp[*sindx-1], h, k, l);
}

FORTRAN_SUBR ( SYSAB, sysab,
	       (const int in[3], int *isysab),
	       (const int in[3], int *isysab),
	       (const int in[3], int *isysab))
{
  int h,k,l;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: SYSAB");)

  h = in[0]; k = in[1]; l = in[2]; 

  *isysab = ccp4spg_is_sysabs(spacegroup, h, k, l);
}

FORTRAN_SUBR ( CCP4SPG_F_IS_SYSABS, ccp4spg_f_is_sysabs,
	       (const int *sindx, const int in[3], int *isysab),
	       (const int *sindx, const int in[3], int *isysab),
	       (const int *sindx, const int in[3], int *isysab))
{
  int h,k,l;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: CCP4SPG_F_IS_SYSABS");)

  if (*sindx <= 0 || *sindx > MSPAC) {
    printf("Error in CCP4SPG_F_IS_SYSABS: sindx %d out of range!\n",*sindx);
    return;
  }

  if ( ! spacegrp[*sindx-1] ) {
    printf("CCP4SPG_F_IS_SYSABS: No spacegroup loaded on channel %d ! \n",*sindx);
    return;
  }

  h = in[0]; k = in[1]; l = in[2]; 

  *isysab = ccp4spg_is_sysabs(spacegrp[*sindx-1], h, k, l);
}

/** Set up centric zones based on symmetry operators.
 * Convention: translations are in rsm[isym][3][*]
 * @param nsm number of symmetry matrices passed.
 * @param rsm symmetry matrices.
 * @param iprint If iprint > 0 then a summary of centric zones is printed.
 */
FORTRAN_SUBR ( CENTRIC, centric,
	       (const int *nsm, const float rsm[192][4][4],
		const int *iprint),
	       (const int *nsm, const float rsm[192][4][4],
		const int *iprint),
	       (const int *nsm, const float rsm[192][4][4],
		const int *iprint))
{
  int i,k,l;
  ccp4_symop *op1=NULL;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: CENTRIC");)

  /* identify spacegroup from supplied symops */
  op1 = (ccp4_symop *) ccp4_utils_malloc(*nsm*sizeof(ccp4_symop));
  for (i = 0; i < *nsm; ++i) {
    for (k = 0; k < 3; ++k) {
      for (l = 0; l < 3; ++l) {
	op1[i].rot[k][l] = rsm[i][l][k];
      }
      op1[i].trn[k] = rsm[i][3][k];
    }
  }
  spacegroup = ccp4_spgrp_reverse_lookup(*nsm,op1);

  if (spacegroup && *iprint > 0) ccp4spg_print_centric_zones(spacegroup);

  free(op1);
}

FORTRAN_SUBR ( CENTR, centr,
	       (const int ih[3], int *ic),
	       (const int ih[3], int *ic),
	       (const int ih[3], int *ic))
{
  int h,k,l;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: CENTR");)

  h = ih[0]; k = ih[1]; l = ih[2]; 

  *ic = ccp4spg_is_centric(spacegroup, h, k, l);

  if (*ic == -1) ccperror(1,"Fatal error in CENTR.");
}

FORTRAN_SUBR ( CCP4SPG_F_IS_CENTRIC, ccp4spg_f_is_centric,
	       (const int *sindx, const int ih[3], int *ic),
	       (const int *sindx, const int ih[3], int *ic),
	       (const int *sindx, const int ih[3], int *ic))
{
  int h,k,l;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: CCP4SPG_F_IS_CENTRIC");)

  if (*sindx <= 0 || *sindx > MSPAC) {
    printf("Error in CCP4SPG_F_IS_CENTRIC: sindx %d out of range!\n",*sindx);
    return;
  }

  if ( ! spacegrp[*sindx-1] ) {
    printf("CCP4SPG_F_IS_CENTRIC: No spacegroup loaded on channel %d ! \n",*sindx);
    return;
  }

  h = ih[0]; k = ih[1]; l = ih[2]; 

  *ic = ccp4spg_is_centric(spacegrp[*sindx-1], h, k, l);

  if (*ic == -1) ccperror(1,"Fatal error in CCP4SPG_F_IS_CENTRIC.");
}

FORTRAN_SUBR ( CENTPHASE, centphase,
	       (const int ih[3], float *cenphs),
	       (const int ih[3], float *cenphs),
	       (const int ih[3], float *cenphs))
{
  int h,k,l;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: CENTPHASE");)

  h = ih[0]; k = ih[1]; l = ih[2]; 

  if (! ccp4spg_is_centric(spacegroup, h, k, l) ) {
    printf("CENTPHASE: This is not a centric reflection!\n");
    return;
  }

  *cenphs = ccp4spg_centric_phase(spacegroup, h, k, l);
}

FORTRAN_SUBR ( CCP4SPG_F_CENTPHASE, ccp4spg_f_centphase,
	       (const int *sindx, const int ih[3], float *cenphs),
	       (const int *sindx, const int ih[3], float *cenphs),
	       (const int *sindx, const int ih[3], float *cenphs))
{
  int h,k,l;

  CSYMLIB_DEBUG(puts("CSYMLIB_F: CCP4SPG_F_CENTPHASE");)

  if (*sindx <= 0 || *sindx > MSPAC) {
    printf("Error in CCP4SPG_F_IS_CENTPHASE: sindx %d out of range!\n",*sindx);
    return;
  }

  if ( ! spacegrp[*sindx-1] ) {
    printf("CCP4SPG_F_IS_CENTPHASE: No spacegroup loaded on channel %d ! \n",*sindx);
    return;
  }

  h = ih[0]; k = ih[1]; l = ih[2]; 

  if (! ccp4spg_is_centric(spacegrp[*sindx-1], h, k, l) ) {
    printf("CCP4SPG_F_CENTPHASE: This is not a centric reflection!\n");
    return;
  }

  *cenphs = ccp4spg_centric_phase(spacegrp[*sindx-1], h, k, l);
}

/* Returns map limits of a.s.u. in fractional units.
   These are rounded up or down to mimic <= or < respectively.
   In fact, these limits may be larger than 1 a.s.u. but always
   have one corner at the origin */
FORTRAN_SUBR ( SETLIM, setlim,
	       (const int *lspgrp, float xyzlim[3][2]),
	       (const int *lspgrp, float xyzlim[3][2]),
	       (const int *lspgrp, float xyzlim[3][2]))
{ 
  CCP4SPG *tmp_spacegroup;      

  CSYMLIB_DEBUG(puts("CSYMLIB_F: SETLIM");)

  if (!spacegroup || spacegroup->spg_ccp4_num != *lspgrp) {
    /* load spacegroup if necessary */
    /* spacegroup only temporary, as setlim is not expected to
       interact with other calls */
    if ( ! (tmp_spacegroup = ccp4spg_load_by_ccp4_num(*lspgrp)) ) {
      printf("SETLIM: failed to load spacegroup info from SYMINFO! \n");
      return;
    }
    xyzlim[0][1] = tmp_spacegroup->mapasu_ccp4[0];
    xyzlim[1][1] = tmp_spacegroup->mapasu_ccp4[1];
    xyzlim[2][1] = tmp_spacegroup->mapasu_ccp4[2];
    free(tmp_spacegroup); 
  } else {
    xyzlim[0][1] = spacegroup->mapasu_ccp4[0];
    xyzlim[1][1] = spacegroup->mapasu_ccp4[1];
    xyzlim[2][1] = spacegroup->mapasu_ccp4[2];
  }
  xyzlim[0][0] = 0.0;
  xyzlim[1][0] = 0.0;
  xyzlim[2][0] = 0.0;
}

/* Returns map limits of a.s.u. in fractional units.
   These are rounded up or down to mimic <= or < respectively.
   In fact, these limits may be larger than 1 a.s.u. but always
   have one corner at the origin.
   This version uses mapasu_zero limits from sgtbx */
FORTRAN_SUBR ( SETLIM_ZERO, setlim_zero,
	       (const int *lspgrp, float xyzlim[3][2]),
	       (const int *lspgrp, float xyzlim[3][2]),
	       (const int *lspgrp, float xyzlim[3][2]))
{ 
  CCP4SPG *tmp_spacegroup;      

  CSYMLIB_DEBUG(puts("CSYMLIB_F: SETLIM_ZERO");)

  if (!spacegroup || spacegroup->spg_ccp4_num != *lspgrp) {
    /* load spacegroup if necessary */
    /* spacegroup only temporary, as setlim is not expected to
       interact with other calls */
    if ( ! (tmp_spacegroup = ccp4spg_load_by_ccp4_num(*lspgrp)) ) {
      printf("SETLIM_ZERO: failed to load spacegroup info from SYMINFO! \n");
      return;
    }
    xyzlim[0][1] = tmp_spacegroup->mapasu_zero[0];
    xyzlim[1][1] = tmp_spacegroup->mapasu_zero[1];
    xyzlim[2][1] = tmp_spacegroup->mapasu_zero[2];
    free(tmp_spacegroup); 
  } else {
    xyzlim[0][1] = spacegroup->mapasu_zero[0];
    xyzlim[1][1] = spacegroup->mapasu_zero[1];
    xyzlim[2][1] = spacegroup->mapasu_zero[2];
  }
  xyzlim[0][0] = 0.0;
  xyzlim[1][0] = 0.0;
  xyzlim[2][0] = 0.0;
}

FORTRAN_SUBR ( SETGRD, setgrd,
	       (const int *nlaue, const float *sample, const int *nxmin,
                const int *nymin, const int *nzmin, int *nx, int *ny, int *nz),
	       (const int *nlaue, const float *sample, const int *nxmin,
                const int *nymin, const int *nzmin, int *nx, int *ny, int *nz),
	       (const int *nlaue, const float *sample, const int *nxmin,
                const int *nymin, const int *nzmin, int *nx, int *ny, int *nz))
{
  int nlaue_save = -1;

  if (!spacegroup) {
    ccperror(2,"SETGRD: No spacegroup loaded yet! \n");
    return;
  }

  if (spacegroup->nlaue != *nlaue) {
    printf("SETGRD: supplied CCP4 Laue code is different from that currently stored\n");
    printf("NLAUE (supplied) = %d\n",*nlaue);
    printf("NLAUE (library)  = %d\n",spacegroup->nlaue);
    printf("(For program FFT and certain spacegroups, this is OK.)\n");
    /* The requested Laue number is different to that for the
       current spacegroup
       Save the current Laue code and load the data for the requested code */
    nlaue_save = spacegroup->nlaue;
    if (ccp4spg_load_laue(spacegroup,*nlaue)) {
      printf("SETGRD: unrecognised CCP4 Laue code, couldn't set FFT grid\n");
      return;
    }
  }
  set_fft_grid(spacegroup, *nxmin, *nymin, *nzmin, *sample, nx, ny, nz);
  if (nlaue_save > -1) {
    /* Restore previous settings */
    ccp4spg_load_laue(spacegroup,nlaue_save);
  }
  return;
}

FORTRAN_SUBR ( FNDSMP, fndsmp,
	       (const int *minsmp, const int *nmul, const float *sample, int *nsampl),
	       (const int *minsmp, const int *nmul, const float *sample, int *nsampl),
	       (const int *minsmp, const int *nmul, const float *sample, int *nsampl))
{

  *nsampl = get_grid_sample(*minsmp, *nmul, *sample);

}

FORTRAN_SUBR ( CALC_ORIG_PS, calc_orig_ps,
	       (fpstr namspg_cif, int *nsym, float rsym[192][4][4], int *norig,
		float orig[96][3], ftn_logical *lpaxisx, ftn_logical *lpaxisy,
		ftn_logical *lpaxisz, int namspg_cif_len),
	       (fpstr namspg_cif, int *nsym, float rsym[192][4][4], int *norig,
		float orig[96][3], ftn_logical *lpaxisx, ftn_logical *lpaxisy,
		ftn_logical *lpaxisz, int namspg_cif_len),
	       (fpstr namspg_cif, int namspg_cif_len, int *nsym, float rsym[192][4][4], 
                int *norig, float orig[96][3], ftn_logical *lpaxisx, 
                ftn_logical *lpaxisy, ftn_logical *lpaxisz))
{
  char *temp_namspg;
  int i,j,k;
  int polarx, polary, polarz;
  float crsym[192][4][4];

  CSYMLIB_DEBUG(puts("CSYMLIB_F: CALC_ORIG_PS");)

  temp_namspg = ccp4_FtoCString(FTN_STR(namspg_cif), FTN_LEN(namspg_cif));
  for (i = 0; i < *nsym; ++i) {
    for (j = 0; j < 4; ++j) {
      for (k = 0; k < 4; ++k) {
        crsym[i][k][j] = rsym[i][j][k];
      }
    }
  }
  *norig = ccp4spg_generate_origins(temp_namspg, *nsym, (const float (*)[4][4])crsym,
                                    orig, &polarx, &polary, &polarz, 1);
  *lpaxisx = polarx ? FORTRAN_LOGICAL_TRUE : FORTRAN_LOGICAL_FALSE;
  *lpaxisy = polary ? FORTRAN_LOGICAL_TRUE : FORTRAN_LOGICAL_FALSE;
  *lpaxisz = polarz ? FORTRAN_LOGICAL_TRUE : FORTRAN_LOGICAL_FALSE;

  free(temp_namspg);
}

static double coefhkl[6];

FORTRAN_SUBR ( SETRSL, setrsl,
	       (const float *a, const float *b, const float *c,
                const float *alpha, const float *beta, const float *gamma),
	       (const float *a, const float *b, const float *c,
                const float *alpha, const float *beta, const float *gamma),
	       (const float *a, const float *b, const float *c,
                const float *alpha, const float *beta, const float *gamma))
{
  float cell[6];

  CSYMLIB_DEBUG(puts("CSYMLIB_F: SETRSL");)

  cell[0] = *a;
  cell[1] = *b;
  cell[2] = *c;
  cell[3] = *alpha;
  cell[4] = *beta;
  cell[5] = *gamma;

  MtzHklcoeffs(cell, coefhkl);

}

FORTRAN_SUBR (STHLSQ1, sthlsq1,
	     (float *reso, const int *ih, const int *ik, const int *il),
	     (float *reso, const int *ih, const int *ik, const int *il),
	     (float *reso, const int *ih, const int *ik, const int *il))
{
  int in[3];

  CSYMLIB_DEBUG(puts("CSYMLIB_F: STHLSQ");)

  in[0] = *ih;
  in[1] = *ik;
  in[2] = *il;

  (*reso) = 0.25*MtzInd2reso(in, coefhkl);
  return;
}

FORTRAN_SUBR (STS3R41, sts3r41,
	     (float *reso, const int *ih, const int *ik, const int *il),
	     (float *reso, const int *ih, const int *ik, const int *il),
	     (float *reso, const int *ih, const int *ik, const int *il))
{
  int in[3];

  CSYMLIB_DEBUG(puts("CSYMLIB_F: STS3R4");)

  in[0] = *ih;
  in[1] = *ik;
  in[2] = *il;

  (*reso) = 0.25*MtzInd2reso(in, coefhkl);

  return;
}

/* straight translation, needs to be done properly, used in phistats */

FORTRAN_SUBR ( HANDCHANGE, handchange,
	       (const int *lspgrp, float *cx, float *cy, float *cz),
	       (const int *lspgrp, float *cx, float *cy, float *cz),
	       (const int *lspgrp, float *cx, float *cy, float *cz))
{
  CSYMLIB_DEBUG(puts("CSYMLIB_F: HANDCHANGE");)

  switch (*lspgrp) {
  case 80:
    *cx=0.0;
    *cy=0.5;
    *cz=0.0;
    break;
  case 98:
    *cx=0.0;
    *cy=0.5;
    *cz=0.25;
    break;
  case 210:
    *cx=0.75;
    *cy=0.25;
    *cz=0.75;
    break;
  case 214:
    *cx=0.25;
    *cy=0.25;
    *cz=0.25;
    break;
  }

}
