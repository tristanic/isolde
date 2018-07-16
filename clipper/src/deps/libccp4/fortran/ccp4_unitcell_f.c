/*
     ccp4_unitcell_f.c: Fortran API to ccp4_unitcell.c
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

/** @file ccp4_unitcell_f.c
 *  Fortran API to ccp4_unitcell.c.
 *  Martyn Winn
 */

#include "../ccp4/ccp4_fortran.h"
#include "../ccp4/ccp4_unitcell.h"
/* rcsid[] = "$Id$" */

/* from input cell and orthogonalisation code, find orthogonalisation
   and fractionalisation matrices. Returns cell volume. */

FORTRAN_SUBR ( CCP4UC_F_FRAC_ORTH_MAT, ccp4uc_f_frac_orth_mat,
          (const float cell[6], const int *ncode, 
	   float ro[3][3], float rf[3][3], float *volume),
          (const float cell[6], const int *ncode, 
	   float ro[3][3], float rf[3][3], float *volume),
          (const float cell[6], const int *ncode, 
	   float ro[3][3], float rf[3][3], float *volume))
{
  int i,j;
  double ro_cmat[3][3], rf_cmat[3][3], dcell[6];

  for (i = 0; i < 6; ++i) 
    dcell[i] = (double) cell[i];

  *volume =  (float) ccp4uc_frac_orth_mat(dcell, *ncode, ro_cmat, rf_cmat);
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      ro[i][j] = (float) ro_cmat[j][i];
      rf[i][j] = (float) rf_cmat[j][i];
    }
  }
}

FORTRAN_SUBR ( CCP4UC_F_CALC_RCELL, ccp4uc_f_calc_rcell,
          (const float cell[6], float rcell[6], float *rvolume),
          (const float cell[6], float rcell[6], float *rvolume),
          (const float cell[6], float rcell[6], float *rvolume))
{
  int i;
  double dcell[6],drcell[6];

  for (i = 0; i < 6; ++i) 
    dcell[i] = (double) cell[i];

  *rvolume = (float) ccp4uc_calc_rcell(dcell, drcell);

  for (i = 0; i < 6; ++i) 
    rcell[i] = (float) drcell[i];

}

FORTRAN_SUBR ( CCP4UC_F_ORTH_TO_FRAC, ccp4uc_f_orth_to_frac,
          (const float rf[3][3], const float xo[3], float xf[3]),
          (const float rf[3][3], const float xo[3], float xf[3]),
	  (const float rf[3][3], const float xo[3], float xf[3]))
{
  int i,j;
  double rf_cmat[3][3], dxo[3], dxf[3];

  for (i = 0; i < 3; ++i) {
    dxo[i] = (double) xo[i];
    for (j = 0; j < 3; ++j) 
      rf_cmat[i][j] = (double) rf[j][i];
  }
  ccp4uc_orth_to_frac((const double (*)[3])rf_cmat, dxo, dxf);
  for (i = 0; i < 3; ++i) 
    xf[i] = (float) dxf[i];

}

FORTRAN_SUBR ( CCP4UC_F_FRAC_TO_ORTH, ccp4uc_f_frac_to_orth,
          (const float ro[3][3], const float xf[3], float xo[3]),
          (const float ro[3][3], const float xf[3], float xo[3]),
	  (const float ro[3][3], const float xf[3], float xo[3]))
{
  int i,j;
  double ro_cmat[3][3], dxf[3], dxo[3];

  for (i = 0; i < 3; ++i) {
    dxf[i] = (double) xf[i];
    for (j = 0; j < 3; ++j) 
      ro_cmat[i][j] = (double) ro[j][i];
  }
  ccp4uc_orth_to_frac((const double (*)[3])ro_cmat, dxf, dxo);
  for (i = 0; i < 3; ++i) 
    xo[i] = (float) dxo[i];

}

FORTRAN_SUBR ( CCP4UC_F_ORTHU_TO_FRACU, ccp4uc_f_orthu_to_fracu,
          (const float rf[3][3], const float uo[3], float uf[3]),
          (const float rf[3][3], const float uo[3], float uf[3]),
	  (const float rf[3][3], const float uo[3], float uf[3]))
{
  int i,j;
  double rf_cmat[3][3], duo[3], duf[3];

  for (i = 0; i < 3; ++i) {
    duo[i] = (double) uo[i];
    for (j = 0; j < 3; ++j) 
      rf_cmat[i][j] = (double) rf[j][i];
  }
  ccp4uc_orthu_to_fracu((const double (*)[3])rf_cmat, duo, duf);
  for (i = 0; i < 3; ++i) 
    uf[i] = (float) duf[i];

}

FORTRAN_SUBR ( CCP4UC_F_FRACU_TO_ORTHU, ccp4uc_f_fracu_to_orthu,
          (const float ro[3][3], const float uf[3], float uo[3]),
          (const float ro[3][3], const float uf[3], float uo[3]),
	  (const float ro[3][3], const float uf[3], float uo[3]))
{
  int i,j;
  double ro_cmat[3][3], duf[3], duo[3];

  for (i = 0; i < 3; ++i) {
    duf[i] = (double) uf[i];
    for (j = 0; j < 3; ++j) 
      ro_cmat[i][j] = (double) ro[j][i];
  }
  ccp4uc_orthu_to_fracu((const double (*)[3])ro_cmat, duf, duo);
  for (i = 0; i < 3; ++i) 
    uo[i] = (float) duo[i];

}

FORTRAN_SUBR ( CELLCHK, cellchk,
	       (const float cell1[6], const float cell2[6], const float *errfrc, int *ierr),
	       (const float cell1[6], const float cell2[6], const float *errfrc, int *ierr),
	       (const float cell1[6], const float cell2[6], const float *errfrc, int *ierr))
{
  int i;
  double dcell1[6], dcell2[6];

  for (i = 0; i < 6; ++i) {
    dcell1[i] = (double) cell1[i];
    dcell2[i] = (double) cell2[i];
  }

  *ierr = ccp4uc_cells_differ(dcell1, dcell2, (double) *errfrc);

}

