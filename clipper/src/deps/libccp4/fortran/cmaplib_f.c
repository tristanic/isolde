/*
     cmaplib_f.c: Fortran API for accessing CCP4 map files
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
/*
*/

/** @page cmap_f_page Fortran API to CMAP 
 *
 *  @section cmap_f_file_list File list

<ul>
<li>cmaplib_f.c
</ul>

 *  @section cmap_f_overview Overview

This library consists of a set of wrappers to the CMAP library
giving the same API as the original maplib.f For details of the
API, see the original <a href="../maplib.html">documentation</a>.
This document covers some peculiarities of the C implementation.

 */
 
/** @file cmaplib_f.c
 *
 *  @brief Fortran API for input, output and manipulation of CCP4 map files.
 *
 *  @author Charles Ballard
 */

#include<string.h>
#include<stdio.h>
#include<fcntl.h>
#include"../ccp4/cmaplib_f.h"
#include"../ccp4/cmap_errno.h"
#include"../ccp4/ccp4_fortran.h"
#include"../ccp4/ccp4_parser.h"
#include"../ccp4/csymlib.h"
#include"../ccp4/ccp4_general.h"
#include"../ccp4/ccp4_program.h"
/* rcsid[] = "$Id$" */

static struct _IOConvMap *ioArray[MAXFILES];
static int last_Read = -1;
static int last_Write = -1;

static int SetChannel()
{
  register int i;
  for (i=0;i!=MAXFILES;i++) 
    if (ioArray[i] == NULL)
      break; 
  return i;
}

static int GetChannel(int iunit)
{
  register int i;
  for (i=0;i!=MAXFILES;i++)
    if (ioArray[i] && ioArray[i]->ipc == iunit)
      break;
  return i;
}

static int ioArrayPrint(IOConvMap *ioMap)
{
  char *filename = ccp4_file_name((*ioMap->mapfile).stream);
  long length = ccp4_file_length((*ioMap->mapfile).stream);
  unsigned rw_mode = ccp4_file_is_read((*ioMap->mapfile).stream);

  if (ccp4VerbosityLevel(-1) < 1) goto out;

  printf("\n Logical Name: %s   Filename: %s \n",
	 ioMap->logname,filename);
  
  if (rw_mode == 1) {
    fprintf(stdout,"\nFile name for input map file on unit %3d : %s\n",
      ioMap->ipc,filename);
    fprintf(stdout,"file size %ld ; logical name %s\n\n",length,ioMap->logname);
  } else {
    fprintf(stdout,"\nFile name for output map file on unit %3d : %s\n",
      ioMap->ipc,filename);
    fprintf(stdout,"logical name %s\n\n",ioMap->logname);
  }
  out:
  free(filename);      /* we strdup it in ccp4_file_name */
  return 1;
}

static int HeaderInit(CMMFile *mfile, const char *title, int mode,
                    const int *order, const int *grid,
                    int nc_start, int nc_end, int nr_start, int nr_end,
                    int ns_start, int nsecs, int lspgrp, const float *cell,
                    float min, float max, double mean, double rms)
{
  int tmp_array[3];

  ccp4_cmap_set_datamode(mfile,mode);
  ccp4_cmap_set_cell(mfile,cell);
  ccp4_cmap_set_label(mfile,title,0);
  ccp4_cmap_set_grid(mfile,grid);
  ccp4_cmap_set_order(mfile,order);
  ccp4_cmap_set_mapstats(mfile,min,max,mean,rms);
  ccp4_cmap_set_spacegroup(mfile,lspgrp);
  
  tmp_array[0] = nc_end - nc_start +1;
  tmp_array[1] = nr_end - nr_start +1;
  tmp_array[2] = nsecs;
  ccp4_cmap_set_dim(mfile,tmp_array);

  tmp_array[0] = nc_start;
  tmp_array[1] = nr_start;
  tmp_array[2] = ns_start;
  ccp4_cmap_set_origin(mfile,tmp_array);

  return (1);
}

static int HeaderReturn(const CMMFile *mfile, char *title, int *mode,
                int *order, int *grid,
                int *nc_start, int *nc_end, int *nr_start,
                int *nr_end, int *ns_start, int *nsecs, int *spacegroup,
                float *cell, float *min, float *max, double *mean,
                double *rms)
{
  int tmp_array[3];
  char *tmp_label;

  if (ccp4_cmap_number_label(mfile)) 
    if ((tmp_label = ccp4_cmap_get_label(mfile,0)) != NULL)
      strcpy(title,tmp_label);

  *mode = ccp4_cmap_get_datamode(mfile);

  ccp4_cmap_get_origin(mfile,tmp_array);
  *nc_start = tmp_array[0];
  *nr_start = tmp_array[1];
  *ns_start = tmp_array[2];
  
  ccp4_cmap_get_dim(mfile,tmp_array);
  *nc_end = *nc_start + tmp_array[0] -1;
  *nr_end = *nr_start + tmp_array[1] -1;
  *nsecs = tmp_array[2];

  ccp4_cmap_get_mapstats(mfile,min,max,mean,rms);

  ccp4_cmap_get_order(mfile,order);
  ccp4_cmap_get_grid(mfile,grid);
  ccp4_cmap_get_cell(mfile,cell);
  *spacegroup = ccp4_cmap_get_spacegroup(mfile);

  return (1);
}

static int HeaderPrint(const CMMFile *mfile)
{
  static const char axes[]={' ','X','Y','Z'};

  if(!mfile) {
    if (ccp4VerbosityLevel(-1) > 0)
      fprintf(stderr,"WARNING: no header information to print.\n");
    return (0);}

  /* C. Vonrhein: make the output identical to old (pre-5.0) Fortran */
  /*              libraries                                          */
  fprintf(stdout,"\n\n");
  fprintf(stdout,"           Number of columns, rows, sections ...............%5d%5d%5d\n",
          mfile->map_dim[0],mfile->map_dim[1],mfile->map_dim[2]);
  fprintf(stdout,"           Map mode ........................................%5d\n",
	  mfile->data_mode);
  fprintf(stdout,"           Start and stop points on columns, rows, sections %5d%5d%5d%5d%5d%5d\n",
          mfile->origin[0],
          mfile->origin[0]+mfile->map_dim[0] - 1,
          mfile->origin[1],
          mfile->origin[1]+mfile->map_dim[1] - 1,
          mfile->origin[2],
          mfile->origin[2]+mfile->map_dim[2] - 1);
  fprintf(stdout,"           Grid sampling on x, y, z ........................%5d%5d%5d\n",
          mfile->cell_grid[0], mfile->cell_grid[1], mfile->cell_grid[2]);
  fprintf(stdout,"           Cell dimensions .................................%10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
          mfile->cell[0], mfile->cell[1], mfile->cell[2],
          mfile->cell[3], mfile->cell[4], mfile->cell[5]);
  fprintf(stdout,"           Fast, medium, slow axes .........................    %c    %c    %c\n",
          axes[mfile->axes_order[0]],
          axes[mfile->axes_order[1]],
          axes[mfile->axes_order[2]]);
  fprintf(stdout,"           Minimum density .................................%12.5f\n",
          mfile->stats.min);
  fprintf(stdout,"           Maximum density .................................%12.5f\n",
          mfile->stats.max);
  fprintf(stdout,"           Mean density ....................................%12.5f\n",
          (float) mfile->stats.mean);
  fprintf(stdout,"           Rms deviation from mean density .................%12.5f\n",
          (float) mfile->stats.rms);
  fprintf(stdout,"           Space-group .....................................%5d\n",
          mfile->spacegroup);
  fprintf(stdout,"           Number of titles ................................%5d\n",
          ccp4_cmap_number_label(mfile));
  fprintf(stdout,"\n\n");

  fprintf(stdout,"     %-35s\n", "Labels: ");
  {
    int i, nlabels = ccp4_cmap_number_label(mfile);
    for (i=0 ; i != nlabels ; i++)
      fprintf(stdout,"  %s\n",ccp4_cmap_get_label(mfile,i));
  }
  fprintf(stdout,"\n\n");
  return (0);
} 

/* belongs to fortran interface , also provide translators for others?
   or should be in csym lib */
static int SymopToFMat4(char *symchs_begin, char *symchs_end, float *rot)
{
  int no_real =0, no_recip = 0, no_axis = 0;          /* counters */
  int col = 3, nops = 0;
  int nsym = 0;
  float sign = 1.0f, value = 0.0f, value2;
  char *cp, *ptr_symchs = symchs_begin, ch;
  int j,k;                                 /* loop variables */
  unsigned char Isep = FALSE;                             /* parsed seperator? */
 
  /* initialise and clear the relevant */
  /* use knowledge that are using [4][4] for rot */
  for (j = 0 ; j !=4 ; ++j)
    for (k = 0; k !=4 ; ++k)
      rot[(k << 2) +j] = 0.0f;
  rot[(3 << 2) +3] = 1.0f;

  while (ptr_symchs <= symchs_end) {
    ch = *ptr_symchs;

    /* parse it */
    if (isspace(ch) || (ch == '*' && no_axis != 0) ) {
      ++ptr_symchs;
      continue;
    } else if (ch == ',' || ch == '*') {
      ++ptr_symchs;
      if (value == 0.0f && col == 3) {
        /* nothing set, this is a problem */
        fprintf(stderr, "Problem \n");
        return (no_recip)? -(nsym) : nsym ;
        /* salvage what we can */
      } else {
        Isep = TRUE;     /* drop through to evaluation*/
      }
    } else if (ch == 'X' || ch == 'x') {
      no_real++, col = 0;
      if (value == 0.0f) value = sign * 1.0f;
      ++ptr_symchs;
      continue;
    } else if (ch == 'Y' || ch == 'y') {
      no_real++, col = 1;
      if (value == 0.0f) value = sign * 1.0f;
      ++ptr_symchs;
      continue;
    } else if (ch == 'Z' || ch == 'z') {
      no_real++, col = 2;
      if (value == 0.0f) value = sign * 1.0f;
      ++ptr_symchs;
      continue;
    } else if (ch == 'H' || ch == 'h') {
      no_recip++, col = 0;
      if (value == 0.0f) value = sign * 1.0f;
      ++ptr_symchs;
      continue;
    } else if (ch == 'K' || ch == 'k') {
      no_recip++, col = 1;
      if (value == 0.0f) value = sign * 1.0f;
      ++ptr_symchs;
      continue;
    } else if (ch == 'L' || ch == 'l') {
      no_recip++, col = 2; 
      if (value == 0.0f) value = sign * 1.0f;
      ++ptr_symchs;
      continue;
    } else if (ch == 'A' || ch == 'a') {
      no_axis++, col = 0;                  
      if (value == 0.0f) value = sign * 1.0f;             
      ++ptr_symchs;
      continue;
    } else if (ch == 'B' || ch == 'b') { 
      no_axis++, col = 1;
      if (value == 0.0f) value = sign * 1.0f;
      ++ptr_symchs;
      continue; 
    } else if (ch == 'C' || ch == 'c') {
      no_axis++, col = 2;
      if (value == 0.0f) value = sign * 1.0f;
      ++ptr_symchs;  
      continue;
    } else if (ch == '+' || ch == '-') {
      sign = ((ch == '+')? 1.0f : -1.0f) ;
      ++ptr_symchs;
      if ( value == 0.0f && col == 3)
        continue; 
      /* drop through to evaluation */
    } else if ( ch == '/') {
      ++ptr_symchs;
      if (value == 0.0f) {
        fprintf(stderr,"ooops\n"); 
        /* error */
      } 
      value2 = strtod(ptr_symchs, &cp);
      if (!value2) return (1);
      value = (float) value/value2;
      ptr_symchs = cp; 
      continue; 
    } else if ( isdigit(ch) || ch == '.') {
      value = strtod(ptr_symchs, &cp);
      ptr_symchs = cp;
      continue;
    } else {
      ++ptr_symchs;
      continue;
    } 
      
    /* value to be entered in rot */
      
    rot[(((nsym << 2) + col) << 2) + nops] = value;
    
    /* have we passed a operator seperator */
    if (Isep) { 
      Isep = FALSE;
      ++nops;  
      sign = 1.0f;
      if (nops == 3) {
        nops = 0;
        ++nsym;
        /* init for next operator */
        for (j = 0 ; j !=4 ; ++j) 
          for (k = 0; k !=4 ; ++k)
            rot[(((nsym << 2) +k) << 2) +j] = 0.0f;
        rot[(((nsym << 2) +3) << 2) +3] = 1.0f;
      
      }
    } 
    /* reset for next cycle */ 
    col = 3;
    value = 0.0f;
    
  }   
      
  /* tidy up */
  if (value)
    rot[(((nsym << 2) + col) << 2) + nops] = value;
      
  return (no_recip)? -(++nsym) : ++nsym ;
      
}     

static int NumberToSymop(char **symchs_begin, int spacegroup)
{ 
  int i,nsym;
  const int n_byt_symop = 80;
  CCP4SPG *spg;            
  
  spg=ccp4spg_load_by_ccp4_num(spacegroup);
  nsym = spg->nsymop;

  /* 
   * We allocate an extra byte for a null.
   * ccp4_cmap_set_symop assumes that the buffer is null terminated for strlen.
   * We really need a strnlen function
   */

  if ( !(*symchs_begin = (char *) malloc( 1 + nsym * n_byt_symop * sizeof(char)) )) {
     ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_AllocFail),
                 "Write_Sym_Gp",NULL);
     return (-1); }

  memset(*symchs_begin,' ',n_byt_symop*nsym);
  (*symchs_begin)[n_byt_symop*nsym] = 0;

  for (i = 0; i < nsym; ++i)
    rotandtrn_to_symop(*symchs_begin+n_byt_symop*i, (*symchs_begin)+n_byt_symop*(i+1), spg->symop[i]);

  ccp4spg_free(&spg);
  
  return nsym;
} 

/**
 * mwrhdl:
 * @iunit:     (I)      map stream number, internal to program
 * @mapnam:    (I)      logical file name
 * @title:     (I)      map title
 * @nsecs:     (I)      number of sections in the map
 * @iuvw:[3]   (I)      fast, medium and slow axes storage
 * @mxyz:[3]   (I)      sampling intervals along X, Y, Z
 * @nw1:       (I)      number of first section
 * @nu1:       (I)      start of section on fast axis
 * @nu2:       (I)      end of section on fast axis
 * @nv1:       (I)      start of section on medium axis
 * @nv2:       (I)      end of section on medium axis
 * @cell:[6]   (I)      cell dimensions in Angstroms and degrees
 * @lspgrp:    (I)      space group number
 * @lmode:     (I)      map data mode =0, FORTRAN logical*1 (char)
 *		                      =1, FORTRAN integer*2 (short?)
 *				      =2, FORTRAN real*4 (float)
 *				      =3, complex integer*2
 *				      =4, complex real
 *				      =5, ==0
 *				      =10, bricked byte map
 *
 * description: Function used to open an output map file and
 * set up the header information.  The actual header 
 * is only written to the file when the file is closed
 */
FORTRAN_SUBR ( MWRHDL, mwrhdl,
    (int *iunit, const fpstr mapnam, const fpstr title, int *nsecs, int iuvw[3], 
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2, 
     float cell[6], int *lspgrp, int *lmode, int mapnam_len, int title_len),
    (int *iunit, const fpstr mapnam, const fpstr title, int *nsecs, int iuvw[3],
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2, 
     float cell[6], int *lspgrp, int *lmode),
    (int *iunit, const fpstr mapnam, int mapnam_len, const fpstr title, 
     int title_len, int *nsecs, int iuvw[3], int mxyz[3], int *nw1, int *nu1, 
     int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode))
{
  char *temp_title, *temp_map, *file;
  int ii;

  temp_map = ccp4_FtoCString(FTN_STR(mapnam), FTN_LEN(mapnam) );
  temp_title = ccp4_FtoCString(FTN_STR(title), FTN_LEN(title) );
  
  if (!(file = getenv(temp_map)))
    file = temp_map;

  if ( (ii = SetChannel()) == MAXFILES )
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MWRHDL", NULL);

  ioArray[ii] = (IOConvMap *) malloc(sizeof(IOConvMap));

  if ((ioArray[ii]->mapfile = (CMMFile *) ccp4_cmap_open(file, O_WRONLY)) == NULL) 
      ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_CantOpenFile),
		  "MWRHDL", NULL);
  ioArray[ii]->ipc = *iunit;
  ioArray[ii]->logname = strdup(temp_map);

  ioArrayPrint(ioArray[ii]);

  HeaderInit(ioArray[ii]->mapfile, temp_title, *lmode, iuvw, mxyz,
		       *nu1, *nu2, *nv1, *nv2, *nw1, *nsecs,
		       *lspgrp, cell, 0.0f, 0.0f, 0.0, 0.0);

  /* Record for FORTRAN API */
  last_Write = ii;
  free(temp_title);
  free(temp_map);
}

FORTRAN_SUBR ( CCP4_MAP_WRITE_OPEN_HEADER_BY_NAME, 
	       ccp4_map_write_open_header_by_name,
    (int *iunit, const fpstr mapnam, const fpstr title, int *nsecs, int iuvw[3], 
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2, 
     float cell[6], int *lspgrp, int *lmode, int mapnam_len, int title_len),
    (int *iunit, const fpstr mapnam, const fpstr title, int *nsecs, int iuvw[3], 
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2, 
     float cell[6], int *lspgrp, int *lmode),
    (int *iunit, const fpstr mapnam, int mapnam_len, const fpstr title, 
     int title_len, int *nsecs, int iuvw[3], int mxyz[3], int *nw1, int *nu1, 
     int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode))
{
  char *temp_title, *temp_map, *file;
  int ii;

  temp_map = ccp4_FtoCString(FTN_STR(mapnam), FTN_LEN(mapnam) );
  temp_title = ccp4_FtoCString(FTN_STR(title), FTN_LEN(title) );
  
  if (!(file = getenv(temp_map)))
    file = temp_map;

  if ( (ii = SetChannel()) == MAXFILES )
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MWRHDL", NULL);

  ioArray[ii] = (IOConvMap *) malloc(sizeof(IOConvMap));

  if ((ioArray[ii]->mapfile = (CMMFile *) ccp4_cmap_open(file, O_WRONLY)) == NULL) 
      ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_CantOpenFile),
		  "MWRHDL", NULL);
  ioArray[ii]->ipc = *iunit;
  ioArray[ii]->logname = strdup(temp_map);

  ioArrayPrint(ioArray[ii]);

  HeaderInit(ioArray[ii]->mapfile, temp_title, *lmode, iuvw, mxyz,
		       *nu1, *nu2, *nv1, *nv2, *nw1, *nsecs,
		       *lspgrp, cell, 0.0f, 0.0f, 0.0, 0.0);

  /* Record for FORTRAN API */
  last_Write = ii;
  free(temp_title);
  free(temp_map);
}

/**
 * mwrhdr:
 * @iunit:     (I)      map stream number, internal to program
 * @title:     (I)      map title
 * @nsecs:     (I)      number of sections in the map
 * @iuvw:[3]   (I)      fast, medium and slow axes storage
 * @mxyz:[3]   (I)      sampling intervals along X, Y, Z
 * @nw1:       (I)      number of first section
 * @nu1:       (I)      start of section on fast axis
 * @nu2:       (I)      end of section on fast axis
 * @nv1:       (I)      start of section on medium axis
 * @nv2:       (I)      end of section on medium axis
 * @cell:[6]   (I)      cell dimensions in Angstroms and degrees
 * @lspgrp:    (I)      space group number
 * @lmode:     (I)      map data mode =0, FORTRAN logical*1 (char)
 *                                    =1, FORTRAN integer*2 (short?)
 *                                    =2, FORTRAN real*4 (float)
 *              	      	      =3, complex integer*2
 *				      =4, complex real
 *				      =5, ==0
 *				      =10, bricked byte map
 *
 * description: similar to mwrhdl().  Difference is that the logical filename for
 * the map is assumed to be %MAPOUT.
 */
FORTRAN_SUBR ( MWRHDR, mwrhdr,
    (int *iunit, const fpstr title, int *nsecs, int iuvw[3], 
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2, 
     float cell[6], int *lspgrp, int *lmode, int title_len),
    (int *iunit, const fpstr title, int *nsecs, int iuvw[3], 
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2, 
     float cell[6], int *lspgrp, int *lmode),
    (int *iunit, const fpstr title, 
     int title_len, int *nsecs, int iuvw[3], int mxyz[3], int *nw1, int *nu1, 
     int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode))
{
  char *temp_title, *file;
  int ii;

  temp_title = ccp4_FtoCString(FTN_STR(title), FTN_LEN(title) );
  
  if (!(file = getenv("MAPOUT")))
    file = "MAPOUT";

  if ( (ii = SetChannel()) == MAXFILES )
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		 "MWRHDR", NULL);
  ioArray[ii] = (IOConvMap *) malloc(sizeof(IOConvMap));

  if ((ioArray[ii]->mapfile = ccp4_cmap_open(file, O_WRONLY)) == NULL) 
      ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_CantOpenFile),
		  "MWRHDR", NULL);
  ioArray[ii]->ipc = *iunit;
  ioArray[ii]->logname = strdup("MAPOUT");

  ioArrayPrint(ioArray[ii]);

  HeaderInit(ioArray[ii]->mapfile, temp_title, *lmode, iuvw, mxyz,
	     *nu1, *nu2, *nv1, *nv2, *nw1, *nsecs,
	     *lspgrp, cell, 0.0f, 0.0f, 0.0, 0.0);

  free(temp_title);

  /* record for FORTRAN API */
  last_Write = ii;
}

FORTRAN_SUBR ( CCP4_MAP_WRITE_OPEN_HEADER_BY_ID, 
	       ccp4_map_write_open_header_by_id,
    (int *iunit, const fpstr title, int *nsecs, int iuvw[3], 
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2, 
     float cell[6], int *lspgrp, int *lmode, int title_len),
    (int *iunit, const fpstr title, int *nsecs, int iuvw[3], 
     int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2, 
     float cell[6], int *lspgrp, int *lmode),
    (int *iunit, const fpstr title, 
     int title_len, int *nsecs, int iuvw[3], int mxyz[3], int *nw1, int *nu1, 
     int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode))
     /*  see MWRHDR */
{
  char *temp_title, *file;
  int ii;

  temp_title = ccp4_FtoCString(FTN_STR(title), FTN_LEN(title) );
  
  if (!(file = getenv("MAPOUT")))
    file = "MAPOUT";

  if ( (ii = SetChannel()) == MAXFILES )
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		 "MWRHDR", NULL);
  ioArray[ii] = (IOConvMap *) malloc(sizeof(IOConvMap));

  if ((ioArray[ii]->mapfile = ccp4_cmap_open(file, O_WRONLY)) == NULL) 
      ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_CantOpenFile),
		  "MWRHDL", NULL);
  ioArray[ii]->ipc = *iunit;
  ioArray[ii]->logname = strdup("MAPOUT");

  ioArrayPrint(ioArray[ii]);

  HeaderInit(ioArray[ii]->mapfile, temp_title, *lmode, iuvw, mxyz,
	     *nu1, *nu2, *nv1, *nv2, *nw1, *nsecs,
	     *lspgrp, cell, 0.0f, 0.0f, 0.0, 0.0);

  free(temp_title);

  /* record for FORTRAN API */
  last_Write = ii;
}

/**
 * mrdhds:
 * @iunit:     (I)      stream number
 * @mapnam:    (I)      logical (file) name
 * @title:     (O)      title of map (char[80])
 * @nsec:      (O)      number of sections
 * @iuvw:[3]   (O)      fast, medium, slow axes (1,2,3 for X,Y,Z)
 * @mxyz:[3]   (O)      sampling intervals along x,y,z
 * @nw1:       (O)      first section number
 * @nu1:       (O)      limits on fast axis
 * @nu2:       (O)
 * @nv1:       (O)      limits on medium axis
 * @nv2:       (O)
 * @cell:[6]   (O)      cell dimensions in Angstroms and degrees
 * @lspgrp:    (O)      space-group number
 * @lmode:     (O)      map data mode
 * @rhmin:     (O)      minimum density
 * @rhmax:     (O)      maximum density
 * @rhmean:    (O)      mean densitry
 * @rhrms:     (O)      rms deviation form mean density
 * @ifail:    (I/O)     On input:   =0, stop on error
 *                                  =1, return on error
 *                      On output:  =-1, error
 *                                  =unchanged, no error
 * @iprint:    (I)      =0, silent
 *                      =other, print info
 *
 * description: Read map header from stream @iunit, logical name in @mapnam
 */
FORTRAN_SUBR( MRDHDS, mrdhds,
              (int *iunit, const fpstr mapnam, fpstr title, 
	       int *nsec, int iuvw[3], int mxyz[3], int *nw1, int *nu1, 
	       int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp, 
	       int *lmode, float *rhmin, float *rhmax, float *rhmean, 
	       float *rhrms, int *ifail, int *iprint, int mapnam_len, 
	       int title_len),
              (int *iunit, const fpstr mapnam, fpstr title, int *nsec, 
	       int iuvw[3], int mxyz[3], int *nw1, int *nu1, int *nu2, 
	       int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode, 
	       float *rhmin, float *rhmax, float *rhmean, float * rhrms, 
	       int *ifail, int *iprint),
              (int *iunit, const fpstr mapnam, int mapnam_len,  
	       fpstr title, int title_len, int *nsec, int iuvw[3], 
	       int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2, 
	       float cell[6], int *lspgrp, int *lmode, float *rhmin, 
	       float *rhmax, float *rhmean, float * rhrms, int *ifail, 
	       int *iprint))
{
  char temp_title[81], *temp_map, *file;
  int ii;
  double drhmean,drhrms;

  temp_map = ccp4_FtoCString(FTN_STR(mapnam), FTN_LEN(mapnam) );

  if (!(file = getenv(temp_map)))
    file = temp_map;

  if ( (ii = SetChannel()) == MAXFILES )
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MRDHDS", NULL);

  ioArray[ii] = (IOConvMap *) malloc(sizeof(IOConvMap));

  if ((ioArray[ii]->mapfile = ccp4_cmap_open(file, O_RDONLY)) == NULL) {
    if (!(*ifail)) {
      ccp4_signal(CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_CantOpenFile),
		  "MRDHDS", NULL);
      ccperror(1,"Error in opening input map file.");
    } else {
      ccp4_signal(CMAP_ERRNO(CMERR_CantOpenFile),
		  "MRDHDS", NULL);
      *ifail = -1; 
      return ;
    }
  }
  ioArray[ii]->ipc = *iunit;
  ioArray[ii]->logname = strdup(temp_map);

  HeaderReturn(ioArray[ii]->mapfile, temp_title, lmode, iuvw, 
	       mxyz, nu1, nu2, nv1, nv2, nw1, 
	       nsec, lspgrp, cell, rhmin, rhmax, &drhmean, &drhrms);
  *rhmean = (float) drhmean;
  *rhrms = (float) drhrms;
  if (*iprint != 0) {
    ioArrayPrint(ioArray[ii]);
    HeaderPrint(ioArray[ii]->mapfile);
  }

  strncpy(title,temp_title,MIN(strlen(temp_title),FTN_LEN(title)));
  /* free malloc'ed memory */
  free(temp_map);
  /* record for FORTRAN API */
  last_Read = ii;
}

FORTRAN_SUBR( CCP4_MAP_READ_OPEN_HEADER_CHECK, 
	      ccp4_map_read_open_header_check,
              (int *iunit, const fpstr mapnam, fpstr title, 
	       int *nsec, int iuvw[3], int mxyz[3], int *nw1, int *nu1, 
	       int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp, 
	       int *lmode, float *rhmin, float *rhmax, float *rhmean, 
	       float *rhrms, int *ifail, int *iprint, int mapnam_len, 
	       int title_len),
              (int *iunit, const fpstr mapnam, fpstr title, int *nsec, 
	       int iuvw[3], int mxyz[3], int *nw1, int *nu1, int *nu2, 
	       int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode, 
	       float *rhmin, float *rhmax, float *rhmean, float * rhrms, 
	       int *ifail, int *iprint),
              (int *iunit, const fpstr mapnam, int mapnam_len,  
	       fpstr title, int title_len, int *nsec, int iuvw[3], 
	       int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2, 
	       float cell[6], int *lspgrp, int *lmode, float *rhmin, 
	       float *rhmax, float *rhmean, float * rhrms, int *ifail, 
	       int *iprint))
     /* see MRDHDS */
{
  char temp_title[81], *temp_map, *file;
  int ii;
  double drhmean,drhrms;

  temp_map = ccp4_FtoCString(FTN_STR(mapnam), FTN_LEN(mapnam) );

  if (!(file = getenv(temp_map)))
    file = temp_map;

  if ( (ii = SetChannel()) == MAXFILES )
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MRDHDS", NULL);

  ioArray[ii] = (IOConvMap *) malloc(sizeof(IOConvMap));

  if ((ioArray[ii]->mapfile = ccp4_cmap_open(file, O_RDONLY)) == NULL) {
    if (!(*ifail)) {
      ccp4_signal(CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_CantOpenFile),
		  "MRDHDS", NULL);
      ccperror(1,"Error in opening input map file.");
    } else {
      ccp4_signal(CMAP_ERRNO(CMERR_CantOpenFile),
		  "MRDHDS", NULL);
    *ifail = -1; 
    return;
    }
  }
  ioArray[ii]->ipc = *iunit;
  ioArray[ii]->logname = strdup(temp_map);

  HeaderReturn(ioArray[ii]->mapfile, temp_title, lmode, iuvw, 
	       mxyz, nu1, nu2, nv1, nv2, nw1, 
	       nsec, lspgrp, cell, rhmin, rhmax, &drhmean, &drhrms);
  *rhmean = (float) drhmean;
  *rhrms = (float) drhrms;
  if (*iprint != 0) {
    ioArrayPrint(ioArray[ii]);
    HeaderPrint(ioArray[ii]->mapfile);
  }

  strncpy(title,temp_title,MIN(strlen(temp_title),FTN_LEN(title)));
  /* free malloc'ed memory */
  free(temp_map);
  /* record for FORTRAN API */
  last_Read = ii;
}

/**
 * mrdhdr:
 * @iunit:     (I)      stream number
 * @mapnam:    (I)      logical (file) name
 * @title:     (O)      title of map (char[80])
 * @nsec:      (O)      number of sections
 * @iuvw:[3]   (O)      fast, medium, slow axes (1,2,3 for X,Y,Z)
 * @mxyz:[3]   (O)      sampling intervals along x,y,z
 * @nw1:       (O)      first section number
 * @nu1:       (O)      limits on fast axis
 * @nu2:       (O)
 * @nv1:       (O)      limits on medium axis
 * @nv2:       (0)      
 * @cell:[6]   (O)      cell dimensions in Angstroms and degrees
 * @lspgrp:    (O)      space-group number
 * @lmode:     (O)      map data mode
 * @rhmin:     (O)      minimum density
 * @rhmax:     (O)      maximum density
 * @rhmean:    (O)      mean density deviation
 * @rhrms:     (O)      rms deviation
 *
 * description: Read map header from stream @iunit, logical name in @mapnam.
 * Differs from mrdhds() in that there are no IFAIL nor IPRINT
 */
FORTRAN_SUBR( MRDHDR, mrdhdr,
              (int *iunit, const fpstr mapnam, fpstr title, 
	       int *nsec, int iuvw[3], int mxyz[3], int *nw1, int *nu1, 
	       int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp, 
	       int *lmode, float *rhmin, float *rhmax, float *rhmean, 
	       float *rhrms, int mapnam_len, int title_len),
              (int *iunit, const fpstr mapnam, fpstr title, int *nsec, 
	       int iuvw[3], int mxyz[3], int *nw1, int *nu1, int *nu2, 
	       int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode, 
	       float *rhmin, float *rhmax, float *rhmean, float * rhrms),
              (int *iunit, const fpstr mapnam, int mapnam_len,  
	       fpstr title, int title_len, int *nsec, int iuvw[3], 
	       int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2, 
	       float cell[6], int *lspgrp, int *lmode, float *rhmin, 
	       float *rhmax, float *rhmean, float * rhrms))
{  
  char temp_title[81], *temp_map, *file;
  int ii;
  double drhmean,drhrms;

  temp_map = ccp4_FtoCString(FTN_STR(mapnam), FTN_LEN(mapnam) );

  if (!(file = getenv(temp_map)))
    file = temp_map;

  if ( (ii = SetChannel()) == MAXFILES )
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MRDHDR", NULL);

  ioArray[ii] = (IOConvMap *) malloc(sizeof(IOConvMap));

  if ((ioArray[ii]->mapfile = ccp4_cmap_open(file, O_RDONLY)) == NULL) {
       ccp4_signal(CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_CantOpenFile),
		  "MRDHDR", NULL);
       ccperror(1,"Error in opening input map file.");
  }
  
  ioArray[ii]->ipc = *iunit;
  ioArray[ii]->logname = strdup(temp_map);

  HeaderReturn(ioArray[ii]->mapfile, temp_title, lmode, iuvw, 
	       mxyz, nu1, nu2, nv1, nv2, nw1, 
	       nsec, lspgrp, cell, rhmin, rhmax, &drhmean, &drhrms);
  *rhmean = (float) drhmean;
  *rhrms = (float) drhrms;
  ioArrayPrint(ioArray[ii]);
  HeaderPrint(ioArray[ii]->mapfile);

  ccp4_CtoFString(FTN_STR(title),FTN_LEN(title),temp_title);
  /* free malloc'ed memory */
  free(temp_map);
  /* record for FORTRAN API */
  last_Read = ii;
}

FORTRAN_SUBR( CCP4_MAP_READ_OPEN_HEADER, 
	      ccp4_map_read_open_header,
              (int *iunit, const fpstr mapnam, fpstr title, 
	       int *nsec, int iuvw[3], int mxyz[3], int *nw1, int *nu1, 
	       int *nu2, int *nv1, int *nv2, float cell[6], int *lspgrp, 
	       int *lmode, float *rhmin, float *rhmax, float *rhmean, 
	       float *rhrms, int mapnam_len, int title_len),
              (int *iunit, const fpstr mapnam, fpstr title, int *nsec, 
	       int iuvw[3], int mxyz[3], int *nw1, int *nu1, int *nu2, 
	       int *nv1, int *nv2, float cell[6], int *lspgrp, int *lmode, 
	       float *rhmin, float *rhmax, float *rhmean, float * rhrms),
              (int *iunit, const fpstr mapnam, int mapnam_len,  
	       fpstr title, int title_len, int *nsec, int iuvw[3], 
	       int mxyz[3], int *nw1, int *nu1, int *nu2, int *nv1, int *nv2, 
	       float cell[6], int *lspgrp, int *lmode, float *rhmin, 
	       float *rhmax, float *rhmean, float * rhrms))
     /* see MRDHDR */
{
  char temp_title[81], *temp_map, *file;
  int ii;
  double drhmean,drhrms;

  temp_map = ccp4_FtoCString(FTN_STR(mapnam), FTN_LEN(mapnam) );

  if (!(file = getenv(temp_map)))
    file = temp_map;

  if ( (ii = SetChannel()) == MAXFILES )
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MRDHDR", NULL);

  ioArray[ii] = (IOConvMap *) malloc(sizeof(IOConvMap));

  if ((ioArray[ii]->mapfile = ccp4_cmap_open(file, O_RDONLY)) == NULL) {
       ccp4_signal(CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_CantOpenFile),
		  "MRDHDR", NULL);
      ccperror(1,"Error in opening input map file.");
  }
  
  ioArray[ii]->ipc = *iunit;
  ioArray[ii]->logname = strdup(temp_map);

  HeaderReturn(ioArray[ii]->mapfile, temp_title, lmode, iuvw, 
	       mxyz, nu1, nu2, nv1, nv2, nw1, 
	       nsec, lspgrp, cell, rhmin, rhmax, &drhmean, &drhrms);
  *rhmean = (float) drhmean;
  *rhrms = (float) drhrms;
  ioArrayPrint(ioArray[ii]);
  HeaderPrint(ioArray[ii]->mapfile);

  strncpy(title,temp_title,MIN(strlen(temp_title),FTN_LEN(title)));
  /* free malloc'ed memory */
  free(temp_map);
  /* record for FORTRAN API */
  last_Read = ii;
}

/**
 * @iunit:     (I)      stream number
 *
 * description: Write out header to mapfile on stream @iunit, and close it.
 */
FORTRAN_SUBR( MWCLOSE, mwclose, (int *iunit), (int *iunit), (int *iunit))
{
  int ii;

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
      ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		  "MWCLOSE", NULL);

  ccp4_cmap_close(ioArray[ii]->mapfile);

  /* record for FORTRAN API */
  free(ioArray[ii]->logname);
  free(ioArray[ii]);
  ioArray[ii] = NULL;
  last_Write = -1;
}

FORTRAN_SUBR( CCP4_MAP_WRITE_CLOSE_AUTO,
	      ccp4_map_write_close_auto, 
	      (int *iunit), (int *iunit), (int *iunit))
     /* see MWCLOSE */
{
  int ii;

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
      ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		  "MWCLOSE", NULL);

  ccp4_cmap_close(ioArray[ii]->mapfile);

  /* record for FORTRAN API */
  free(ioArray[ii]->logname);
  free(ioArray[ii]);
  ioArray[ii] = NULL;
  last_Write = -1;
}

/**
 * mclose:
 * @iunit:     (I)      map stream number
 * @min:       (I)      minimum density in map
 * @max:       (I)      maximum density in map
 * @mean:      (I)      the sum of all the densities in the map.
 *                      (this will be divided internally by the number
 *                       of points in the map to give the mean density
 *                       which is then stored)
 * @rms        (I)      the sum of the squares of the density values in the
 *                      map (this will be used internally to calculate the
 *                      rms deviation from the mean value, which is then
 *                      stored)
 *
 * description:  Write out header to map file on stream @iunit, and close it.
 */
FORTRAN_SUBR( MCLOSE, mclose,
	      (int *iunit, float *min, float *max, float *mean, float *rms),
	      (int *iunit, float *min, float *max, float *mean, float *rms),
	      (int *iunit, float *min, float *max, float *mean, float *rms))
{
  int ii;

  if ( (ii = GetChannel(*iunit)) == MAXFILES ||
       !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MCLOSE", NULL);

  ccp4_cmap_closemode(ioArray[ii]->mapfile, 2);
  ccp4_cmap_set_mapstats(ioArray[ii]->mapfile,*min, *max, (double) *mean, (double) *rms);
  ccp4_cmap_close(ioArray[ii]->mapfile);

  /* record for FORTRAN API */
  if (ioArray[ii]->logname) free(ioArray[ii]->logname);
  free(ioArray[ii]);
  ioArray[ii] = NULL;
  last_Write = -1;
}

FORTRAN_SUBR( CCP4_MAP_WRITE_CLOSE_USER_SUM, 
	      ccp4_map_write_close_user_sum,
	      (int *iunit, float *min, float *max, float *mean, float *rms),
	      (int *iunit, float *min, float *max, float *mean, float *rms),
	      (int *iunit, float *min, float *max, float *mean, float *rms))
     /* see MCLOSE */
{
  int ii;

  if ( (ii = GetChannel(*iunit)) == MAXFILES ||
       !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MCLOSE", NULL);

  ccp4_cmap_closemode(ioArray[ii]->mapfile, 2);
  ccp4_cmap_set_mapstats(ioArray[ii]->mapfile,*min, *max, (double) *mean, (double) *rms);
  ccp4_cmap_close(ioArray[ii]->mapfile);

  /* record for FORTRAN API */
  if (ioArray[ii]->logname) free(ioArray[ii]->logname);
  free(ioArray[ii]);
  ioArray[ii] = NULL;
  last_Write = -1;
}

/**
 * mclosc:
 * @iunit:     (I)      map stream number
 * @min:       (I)      minimum density in map
 * @max:       (I)      maximum density in map
 * @mean:      (I)      the sum of all the densities in the map.
 * @rms:       (I)      the sum of the squares of the density values in the map
 *
 * description: Write out header to map file on stream @iunit, and close it.
 * This routine differs from mwclose() in that the @mean and @rms
 * values are written directly.
 */
FORTRAN_SUBR( MCLOSC, mclosc,
	      (int *iunit, float *min, float *max, float *mean, float *rms),
	      (int *iunit, float *min, float *max, float *mean, float *rms),
	      (int *iunit, float *min, float *max, float *mean, float *rms))
{
  int ii;

  if ( (ii = GetChannel(*iunit)) == MAXFILES ||
       !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MCLOSC", NULL);

  ccp4_cmap_closemode(ioArray[ii]->mapfile, 1);
  ccp4_cmap_set_mapstats(ioArray[ii]->mapfile,*min, *max, (double) *mean, (double) *rms);
  ccp4_cmap_close(ioArray[ii]->mapfile);

  /* record for FORTRAN API */
  if (ioArray[ii]->logname) free(ioArray[ii]->logname);
  free(ioArray[ii]);
  ioArray[ii] = NULL;
  last_Write = -1;
}

FORTRAN_SUBR( CCP4_MAP_WRITE_CLOSE_USER_MEAN, 
	      ccp4_map_write_close_user_mean,
	      (int *iunit, float *min, float *max, float *mean, float *rms),
	      (int *iunit, float *min, float *max, float *mean, float *rms),
	      (int *iunit, float *min, float *max, float *mean, float *rms))
     /* see MCLOSC */
{
  int ii;

  if ( (ii = GetChannel(*iunit)) == MAXFILES ||
       !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MCLOSC", NULL);

  ccp4_cmap_closemode(ioArray[ii]->mapfile, 1);
  ccp4_cmap_set_mapstats(ioArray[ii]->mapfile,*min, *max, (double) *mean, (double) *rms);
  ccp4_cmap_close(ioArray[ii]->mapfile);

  /* record for FORTRAN API */
  if (ioArray[ii]->logname) free(ioArray[ii]->logname);
  free(ioArray[ii]);
  ioArray[ii] = NULL;
  last_Write = -1;
}

/**
 * mrclos:
 * @iunit:     (I)      map stream number
 */
FORTRAN_SUBR( MRCLOS, mrclos, (int *iunit), (int *iunit), (int *iunit))
{
  int ii;

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
	       "MRCLOS", NULL);

  ccp4_cmap_close(ioArray[ii]->mapfile);

  /* record for FORTRAN API */
  if (ioArray[ii]->logname) free(ioArray[ii]->logname);
  free(ioArray[ii]);
  ioArray[ii] = NULL;
  last_Read = -1;
}

FORTRAN_SUBR( CCP4_MAP_READ_CLOSE, 
	      ccp4_map_read_close, 
	      (int *iunit), (int *iunit), (int *iunit))
     /* see MRCLOS */
{
  int ii;

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
	       "MRCLOS", NULL);

  ccp4_cmap_close(ioArray[ii]->mapfile);

  /* record for FORTRAN API */
  if (ioArray[ii]->logname) free(ioArray[ii]->logname);
  free(ioArray[ii]);
  ioArray[ii] = NULL;
  last_Read = -1;
}

/**
 * msywrt:
 * @iunit:     (I)      map stream number
 * @nsym:      (I)      number of symmetry operators
 * @rot:       (I)      rotation/translation matrices (4,4,nsym)
 *
 * description: Write @nsym symmetry operators to iostream @iunit.
 * note: the symmetry operators are written to the file one per
 * line, and may have a different format to those in the 
 * SYMOP file (uses symtr3()).
 */
FORTRAN_SUBR( MSYWRT, msywrt, 
	      (int *iunit, int *nsym, float *rot),
	      (int *iunit, int *nsym, float *rot),
	      (int *iunit, int *nsym, float *rot))
{
  int ii,i,j,k;
  float rsm[4][4];
  char buffer[80];

  if ( (ii = GetChannel(*iunit)) == MAXFILES ||
       !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MSYWRT", NULL);

  for (i=0; i != *nsym ; ++i) {
    memset(buffer,' ',80U);
    for (j=0; j != 4 ; ++j) 
      for (k=0; k != 4 ; ++k) 
        rsm[j][k] = *(rot+16*i+j+4*k);
    mat4_to_symop(buffer,&buffer[80],(const float (*)[4])rsm);
    ccp4_cmap_set_symop(ioArray[ii]->mapfile,buffer);
  }
  /* record for FORTRAN API */
  last_Write = ii;
}

FORTRAN_SUBR( CCP4_MAP_WRITE_SYMM_MATRIX, ccp4_map_write_symm_matrix, 
	      (int *iunit, int *nsym, float *rot),
	      (int *iunit, int *nsym, float *rot),
	      (int *iunit, int *nsym, float *rot))
     /* see MSYWRT */  
{
  int ii,i,j,k;
  float rsm[4][4];
  char buffer[80];

  if ( (ii = GetChannel(*iunit)) == MAXFILES ||
       !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MSYWRT", NULL);

  for (i=0; i != *nsym ; ++i) {
    memset(buffer,' ',80U);
    for (j=0; j != 4 ; ++j) 
      for (k=0; k != 4 ; ++k) 
        rsm[j][k] = *(rot+16*i+j+4*k);
    mat4_to_symop(buffer,&buffer[80],(const float (*)[4])rsm);
    ccp4_cmap_set_symop(ioArray[ii]->mapfile,buffer);
  }
  /* record for FORTRAN API */
  last_Write = ii;
}

/**
 * msymop:
 * @iunit:     (I)      map stream number
 * @nsym:      (O)      number of symmetry operations
 * @rot:       (O)      rotation/translation matrices (4,4,nsym)
 *
 * description: Read symmetry operations from the mapfile on iostream @iunit
 * (after call to mrdhdr() to initialise header).
 * Convert to matrices and vectors.
 * Note: when no symops are found the I matrix is return (P1 spacegroup)
 */
FORTRAN_SUBR( MSYMOP, msymop,
	      (int *iunit, int *nsym, float *rot),
	      (int *iunit, int *nsym, float *rot),
	      (int *iunit, int *nsym, float *rot))
{
  char buffer[81];
  int nsymop, ii, i, result= 0;

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MSYMOP", NULL);

  nsymop = ccp4_cmap_num_symop(ioArray[ii]->mapfile);

  ccp4_cmap_seek_symop(ioArray[ii]->mapfile, 0, SEEK_SET);

  for (i=0 ; i != nsymop ; ++i) {
    ccp4_cmap_get_symop(ioArray[ii]->mapfile,buffer);
    ccp4printf(1,"  Symmetry operations:  %s\n",buffer);
    result += SymopToFMat4(buffer,&buffer[80],rot+16*i);
  }

  if ( result != nsymop) 
   ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_ReadFail),
	       "MSYMOP", NULL);
  else if ( nsymop == 0) { 
    /* if no symops are found, what to return? here the I matrix */
    int k;
    ccp4printf(0,"WARNING: 0 dimension symops, returning P1\n");
    memset(rot, '\0', 16*sizeof(float));
    for (k = 0; k !=4 ; k++)
      rot[(k<<2)+k] = 1.0f; 
      *nsym = 1; }
  else
    *nsym = result;

  /* record for FORTRAN API */
  last_Read = ii;
}

FORTRAN_SUBR( CCP4_MAP_READ_SYMM_MATRIX,
	      ccp4_map_read_symm_matrix,
	      (int *iunit, int *nsym, float *rot),
	      (int *iunit, int *nsym, float *rot),
	      (int *iunit, int *nsym, float *rot))
     /* see MSYMOP */
{
  char buffer[81];
  int ii, i, nsymop, result= 0;

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MSYMOP", NULL);

  nsymop = ccp4_cmap_num_symop(ioArray[ii]->mapfile);

  ccp4_cmap_seek_symop(ioArray[ii]->mapfile, 0, SEEK_SET);

  for (i=0 ; i != nsymop ; ++i) {
    ccp4_cmap_get_symop(ioArray[ii]->mapfile,buffer);
    ccp4printf(1,"  Symmetry operations:  %s",buffer);
    result += SymopToFMat4(buffer,&buffer[80],rot+16*i);
  }

  if ( result != nsymop) 
   ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_ReadFail),
	       "MSYMOP", NULL);
  else if ( nsymop == 0) { 
    /* if no symops are found, what to return? here the I matrix */
    int k;
    ccp4printf(0,"WARNING: 0 dimension symops, returning P1\n");
    memset(rot, '\0', 16*sizeof(float));
    for (k = 0; k !=4 ; k++)
      rot[(k<<2)+k] = 1.0f; 
      *nsym = 1; }
  else
    *nsym = result;

  /* record for FORTRAN API */
  last_Read = ii;
}

/** 
 * msyput:
 * @sunit:     (I)      iostream to use for SYMOP library (ignored)
 * @spacegroup:(I)      spacegroup number
 * @iunit:     (I)      map stream number
 *
 * description:   Copy symmetry operators to map
 * stream @iunit, leaving space at head of file for n_byt_hdr items
 * of header records.
 */
FORTRAN_SUBR( MSYPUT, msyput,
	      (int *sunit, int *spacegroup, int *iunit),
	      (int *sunit, int *spacegroup, int *iunit),
	      (int *sunit, int *spacegroup, int *iunit))
{
  const int n_byt_symop=80;
  int i,ii, nsymop;
  char *symops, *buffer;

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MSYPUT", NULL);

  nsymop = NumberToSymop(&symops,*spacegroup);
  
  for (i=0, buffer=symops ; i != nsymop ; ++i, buffer += n_byt_symop) {
    ccp4_cmap_set_symop(ioArray[ii]->mapfile,buffer);
  }
  free(symops);
  /* record for FORTRAN API */
  last_Write = ii;
}

FORTRAN_SUBR( CCP4_MAP_WRITE_SPGNAME, 
	      ccp4_map_write_spgname,
	      (int *sunit, int *spacegroup, int *iunit),
	      (int *sunit, int *spacegroup, int *iunit),
	      (int *sunit, int *spacegroup, int *iunit))
     /*  see MSYPUT */
{
  const int n_byt_symop=80;
  int i,ii, nsymop;
  char *symops, *buffer;

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MSYPUT", NULL);

  nsymop = NumberToSymop(&symops,*spacegroup);
  
  for (i=0, buffer=symops ; i != nsymop ; ++i, buffer += n_byt_symop) {
    ccp4_cmap_set_symop(ioArray[ii]->mapfile,buffer);
  }

  /* record for FORTRAN API */
  last_Write = ii;
}

/**
 * mspew:
 * @iunit:     (I)      stream number
 * @section:   (I)      array holding the map section.
 *
 * description: Write whole map section to stream @iunit.
 * This routine is only suitable if the whole map section is to
 * be written.
 */
FORTRAN_SUBR( MSPEW, mspew,
              (int *iunit, float *section),
              (int *iunit, float *section),
              (int *iunit, float *section))
{
  int ii;

  if ( (ii = GetChannel(*iunit)) == MAXFILES ||
       !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MSPEW", NULL);

  if (ccp4_cmap_write_section(ioArray[ii]->mapfile, section) == 0) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_WriteFail),
		"MSPEW", NULL);

  /* record for FORTRAN API */
  last_Write = ii;
}

FORTRAN_SUBR( CCP4_MAP_WRITE_ALL_SECTION, 
	      ccp4_map_write_all_section,
              (int *iunit, float *section),
              (int *iunit, float *section),
              (int *iunit, float *section))
     /* see MSPEW */    
{
  int ii;

  if ( (ii = GetChannel(*iunit)) == MAXFILES ||
       !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MSPEW", NULL);

  if (ccp4_cmap_write_section(ioArray[ii]->mapfile, section) == 0) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_WriteFail),
		"MSPEW", NULL);

  /* record for FORTRAN API */
  last_Write = ii;
}

/**
 * mgulp:
 * @iunit:     (I)      map stream number
 * @section:   (O)      array containing section of map data read from file
 * @ier:       (O)      error code =0, OK
 *                                 =non-zero, error or end-of-file
 *
 * description: Read next whole map section from stream @iunit to array @section.
 *     @section is returned in the same mode as in the file, no data
 *     conversion is performed.
 */
FORTRAN_SUBR( MGULP, mgulp,
	      (int *iunit, float *section, int *ier),
	      (int *iunit, float *section, int *ier),
	      (int *iunit, float *section, int *ier))
{
  int ii;

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile)  
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MGULP", NULL);

  if ( (*ier = ccp4_cmap_read_section(ioArray[ii]->mapfile,section)) == 0) 
    ccp4_signal(CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_ReadFail),
	       "MGULP", NULL);

  *ier = (*ier == 0) ? -1 : 0;

  /* record for FORTRAN API */
  last_Read = ii;
}

FORTRAN_SUBR( CCP4_MAP_READ_WHOLE_SECTION_AS_MODE,
	      ccp4_map_read_whole_section_as_mode,
	      (int *iunit, float *section, int *ier),
	      (int *iunit, float *section, int *ier),
	      (int *iunit, float *section, int *ier))
     /* see MGULP */
{
  int ii;

 if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
   ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
	       "MGULP", NULL);

  if ((*ier = ccp4_cmap_read_section(ioArray[ii]->mapfile, section)) == 0) 
   ccp4_signal(CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_ReadFail),
	       "MGULP", NULL);

  *ier = (*ier >= 0) ? 0 : -1;

  /* record for FORTRAN API */
  last_Read = ii;
}

/**
 * mwrsec:
 * @iunit:     (I)      map stream number
 * @section:   (I)      array holding the map section
 * @mu:        (I)      the number of points along the whole fast axis
 * @mv:        (I)      the number of points along the whole medium axis
 * @iu1:       (I)      the starting array index along the fast axis
 * @iu2:       (I)      the final array index along the fast axis
 * @iv1:       (I)      the starting array index along the medium axis
 * @iv2:       (I)      the final array index along the medium axis
 *
 * description: Write part of map section X(MU,MV) to stream @iunit
 *
 * Note: the C routine this calls uses the C/C++ convention of
 * [) rather than [].  Therefore the the @iu1 and @iv1 values
 * are decremented by 1, and the @iu2 and @iv2 are unchanged
 * as paramenters (C uses 0 offset for arrays).  Also note
 * that the dimensions are fall the storage array, not the map.
 */
FORTRAN_SUBR( MWRSEC, mwrsec,
	      (int *iunit, float *section, int *mu, int *mv, int *iu1,
	       int *iu2, int *iv1, int *iv2),
	      (int *iunit, float *section, int *mu, int *mv, int *iu1,
	       int *iu2, int *iv1, int *iv2),
	      (int *iunit, float *section, int *mu, int *mv, int *iu1,
	       int *iu2, int *iv1, int *iv2))
{
  int ii, jj;

  if ( (ii = GetChannel(*iunit)) == MAXFILES ||
       !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MWRSEC", NULL);

  for (jj=*iv1-1 ; jj != *iv2; jj++) 
    if (ccp4_cmap_write_data(ioArray[ii]->mapfile, section+jj*(*mu)+(*iu1)-1, 
                        (*iu2)-(*iu1)+1) == -1) 
      ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_WriteFail),
		  "MWRSEC", NULL);
      
  /* record for FORTRAN API */
  last_Write = ii;
}

FORTRAN_SUBR( CCP4_MAP_WRITE_PART_SECTION, 
	      ccp4_map_write_part_section,
	      (int *iunit, float *section, int *mu, int *mv, int *iu1,
	       int *iu2, int *iv1, int *iv2),
	      (int *iunit, float *section, int *mu, int *mv, int *iu1,
	       int *iu2, int *iv1, int *iv2),
	      (int *iunit, float *section, int *mu, int *mv, int *iu1,
	       int *iu2, int *iv1, int *iv2))
     /* see MWRSEC */
{
  int ii, jj;

  if ( (ii = GetChannel(*iunit)) == MAXFILES ||
       !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MWRSEC", NULL);

  for (jj=*iv1-1 ; jj != *iv2; jj++) 
    if (ccp4_cmap_write_data(ioArray[ii]->mapfile, section+jj*(*mu)+(*iu1)-1, 
			     (*iu2)-(*iu1)+1) == -1) 
      ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_WriteFail),
	    	  "MWRSEC", NULL);
      
  /* record for FORTRAN API */
  last_Write = ii;
}

/**
 * mrdlin:
 * @iunit:     (I)      stream number
 * @line:      (O)      array to contain the line of data from the map
 * @ier:       (O)      error flag =0, OK
 *                                 =non-zero, error or EOF.
 *
 * description: Read the next line of the map from stream @iunit into array
 * @line.  @line is returned in the same mode as in file, 
 * no data conversion is performed.
 */
FORTRAN_SUBR( MRDLIN, mrdlin,
	      (int *iunit, float *line, int *ier),
	      (int *iunit, float *line, int *ier),
	      (int *iunit, float *line, int *ier))
{
  int ii;

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MRDLIN", NULL);

  if ((*ier = ccp4_cmap_read_row(ioArray[ii]->mapfile, line)) == 0) 
   ccp4_signal(CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_ReadFail),
	       "MRDLIN", NULL);

  *ier = (*ier == 0) ? -1 : 0;

  /* record for FORTRAN API */
  last_Read = ii;
}

FORTRAN_SUBR( CCP4_MAP_READ_LINE_AS_MODE, 
	      ccp4_map_read_line_as_mode,
	      (int *iunit, float *line, int *ier),
	      (int *iunit, float *line, int *ier),
	      (int *iunit, float *line, int *ier))
     /* see MRDLIN */
{
  int ii;

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MRDLIN", NULL);

  if ((*ier = ccp4_cmap_read_row(ioArray[ii]->mapfile, line)) == 0) 
   ccp4_signal(CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_ReadFail),
	       "MRDLIN", NULL);

  *ier = (*ier == 0) ? -1 : 0;

  /* record for FORTRAN API */
  last_Read = ii;
}

/**
 * mgulpr:
 * @iunit:     (I)      stream number
 * @section:   (O)      array to contain the section read from the file
 * @ier:       (O)      error flag =0, OK
 *                                 =non-zero, error or EOF.
 *
 * description: Read next whole map section from stream @iunit to array @section.
 * For map modes other than 2, the array is converted to real; for
 * complex maps (MODE = 3 or 4) the complex amplitude is evaluated.
 */
FORTRAN_SUBR( MGULPR, mgulpr,
	      (int *iunit, float *section, int *ier),
	      (int *iunit, float *section, int *ier),
	      (int *iunit, float *section, int *ier))
{
  uint8 *buffer;
  int ii, read_dim, mode;
  int nstore[3];

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MGULPR", NULL);

  ccp4_cmap_get_dim(ioArray[ii]->mapfile,nstore);
  read_dim = nstore[0]*nstore[1];
  
  if ((mode = ccp4_cmap_get_datamode( ioArray[ii]->mapfile)) == FLOAT32)
    buffer = (uint8 *) section;
  else 
    buffer = (uint8 *) calloc(read_dim , ccp4_file_itemsize((*ioArray[ii]->mapfile).stream));
    
  if (ccp4_cmap_read_section(ioArray[ii]->mapfile, buffer) == 0) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_ReadFail),
	       "MGULPR", NULL);
    *ier = -1;
  } else {
    ccp4_utils_translate_mode_float(section, buffer, read_dim, mode);
    *ier = 0;
  }
  /* record for FORTRAN API */
  last_Read = ii;
}

FORTRAN_SUBR( CCP4_MAP_READ_WHOLE_SECT_AS_REAL, 
	      ccp4_map_read_whole_sect_as_real,
	      (int *iunit, float *section, int *ier),
	      (int *iunit, float *section, int *ier),
	      (int *iunit, float *section, int *ier))
     /* see MGULPR */
{
  uint8 *buffer;
  int ii, read_dim, mode;
  int nstore[3];

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MGULPR", NULL);

  ccp4_cmap_get_dim(ioArray[ii]->mapfile,nstore);
  read_dim = nstore[0]*nstore[1];
  
  if ((mode = ccp4_cmap_get_datamode( ioArray[ii]->mapfile)) == FLOAT32)
    buffer = (uint8 *) section;
  else 
    buffer = (uint8 *) calloc(read_dim , ccp4_file_itemsize((*ioArray[ii]->mapfile).stream));

  if (ccp4_cmap_read_section(ioArray[ii]->mapfile, buffer) == 0) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_ReadFail),
	       "MGULPR", NULL);
    *ier = -1;
  } else {
    ccp4_utils_translate_mode_float(section, buffer, read_dim, mode);
    *ier = 0;
  }
    
  /* record for FORTRAN API */
  last_Read = ii;
}

/**
 * mposn:
 * @iunit:     (I)      stream number
 * @sec:       (I)      position the input map before @sec
 *
 * description: Sets the position in the map file so that the next section
 * to be read is @sec.  Note: the section number is input.
 */
FORTRAN_SUBR( MPOSN, mposn,
	      (int *iunit, int *sec),
	      (int *iunit, int *sec),
	      (int *iunit, int *sec))
{
  int ii, origin[3], section;

  if ( (ii = GetChannel(*iunit)) == MAXFILES ||
       !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MPOSN", NULL);

  ccp4_cmap_get_origin(ioArray[ii]->mapfile,origin);
  section = *sec - origin[2];
  
  if ( ccp4_cmap_seek_section(ioArray[ii]->mapfile, 
			      section, SEEK_SET) == -1) 
      ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_SeekFail),
		"MPOSN", NULL);

  /* record for FORTRAN API */
  last_Read = ii;
}

FORTRAN_SUBR( CCP4_MAP_READ_POSITION_SELECTION, 
	      ccp4_map_read_position_selection,
	      (int *iunit, int *sec),
	      (int *iunit, int *sec),
	      (int *iunit, int *sec))
     /* see MPOSN */
{
  int ii, origin[3], section;

  if ( (ii = GetChannel(*iunit)) == MAXFILES ||
       !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MPOSN", NULL);
                
  ccp4_cmap_get_origin(ioArray[ii]->mapfile,origin);
  section = *sec - origin[2];

  if ( ccp4_cmap_seek_section(ioArray[ii]->mapfile, 
			      section, SEEK_SET) == -1) 
      ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_SeekFail),
		"MPOSN", NULL);

  /* record for FORTRAN API */
  last_Read = ii;
}

/**
 * mposnw:
 * @iunit:     (I)      stream number
 * @sec:       (I)      position the input map before @sec
 *
 * description: Sets the position in the map file so that the next section
 * to be written is @sec (used to overwrite)
 */
FORTRAN_SUBR( MPOSNW, mposnw,
	      (int *iunit, int *sec),
	      (int *iunit, int *sec),
	      (int *iunit, int *sec))
{
  int ii, origin[3], section;

  if ( (ii = GetChannel(*iunit)) == MAXFILES ||
       !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MPOSNW", NULL);
                
  ccp4_cmap_get_origin(ioArray[ii]->mapfile,origin);
  section = *sec - origin[2];

  if ( ccp4_cmap_seek_section(ioArray[ii]->mapfile, 
			      section, SEEK_SET) == -1) 
      ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_SeekFail),
		"MPOSNW", NULL);

  /* record for FORTRAN API */
  last_Write = ii;
}

FORTRAN_SUBR( CCP4_MAP_WRITE_POSITION_SECTION, 
	      ccp4_map_write_position_section,
	      (int *iunit, int *sec),
	      (int *iunit, int *sec),
	      (int *iunit, int *sec))
     /* see MPOSNW */
{
  int ii, origin[3], section;

  if ( (ii = GetChannel(*iunit)) == MAXFILES ||
       !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"MPOSNW", NULL);
                
  ccp4_cmap_get_origin(ioArray[ii]->mapfile,origin);
  section = *sec - origin[2];

  if ( ccp4_cmap_seek_section(ioArray[ii]->mapfile, 
			      section, SEEK_SET) == -1) 
      ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_SeekFail),
		"MPOSNW", NULL);

  /* record for FORTRAN API */
  last_Write = ii;
}

/**
 * msycpy:
 * @iunit:     (I)      input map stream number
 * @ounit:     (I)      output map stream number
 *
 * description: Copy symmetry data from file on iostream @iunit to the file
 * on iostream @iunit (after header initialisation calls).
 */
FORTRAN_SUBR( MSYCPY, msycpy,
	      (int *iunit, int *ounit),
	      (int *iunit, int *ounit),
	      (int *iunit, int *ounit))
{  
  int ii,jj, nsym, i;
  char symop[81];

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile
       || (jj = GetChannel(*ounit)) == MAXFILES || !ioArray[jj]->mapfile) 
     ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		 "MSYCPY", NULL);

  nsym = ccp4_cmap_num_symop(ioArray[ii]->mapfile);
  ccp4_cmap_seek_symop(ioArray[ii]->mapfile,0,SEEK_SET);

  for (i=0; i != nsym ; ++i) {
    if (ccp4_cmap_get_symop(ioArray[ii]->mapfile, symop) == 1) 
      ccp4_cmap_set_symop(ioArray[jj]->mapfile, symop);
  }

  last_Read = ii;
  last_Write = jj;
}

FORTRAN_SUBR( CCP4_MAP_COPY_SYMMETRY, 
	      ccp4_map_copy_symmetry,
	      (int *iunit, int *ounit),
	      (int *iunit, int *ounit),
	      (int *iunit, int *ounit))
     /* see MSYCPY */
{  
  int ii,jj, nsym, i;
  char symop[81];

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile
       || (jj = GetChannel(*ounit)) == MAXFILES || !ioArray[jj]->mapfile) 
     ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		 "MSYCPY", NULL);

  nsym = ccp4_cmap_num_symop(ioArray[ii]->mapfile);
  ccp4_cmap_seek_symop(ioArray[ii]->mapfile,0,SEEK_SET);

  for (i=0; i != nsym ; ++i) {
    if (ccp4_cmap_get_symop(ioArray[ii]->mapfile, symop) == 1) 
      ccp4_cmap_set_symop(ioArray[jj]->mapfile, symop);
  }

  last_Read = ii;
  last_Write = jj;
}

/**
 * mttrep:
 * @label:     (I)      new label (max. length 80 chars)
 * @posn:      (I)      position for label in label array.
 *
 * description: Replace label at @posn in header of current output mapfile.
 * Maximum of 10 labels, is exceeded overwrites final label.
 */
FORTRAN_SUBR( MTTREP, mttrep,
	      (const fpstr label, int *posn, int label_len),
	      (const fpstr label, int *posn),
	      (const fpstr label, int label_len, int *posn))
{
  char *temp_label;

  /* no output map, nothing to do */
  if (last_Write == -1) return;

  temp_label = ccp4_FtoCString(FTN_STR(label), FTN_LEN(label));

  ccp4_cmap_set_label(ioArray[last_Write]->mapfile, temp_label,*posn-1 );

  free(temp_label);
}

FORTRAN_SUBR( CCP4_MAP_WRITE_REPLACE_TITLE, 
	      ccp4_map_write_replace_title,
	      (const fpstr label, int *posn, int label_len),
	      (const fpstr label, int *posn),
	      (const fpstr label, int label_len, int *posn))
     /* see MTTREP */
{
  char *temp_label;

  /* no output map, nothing to do */
  if (last_Write == -1) return;

  temp_label = ccp4_FtoCString(FTN_STR(label), FTN_LEN(label));

  ccp4_cmap_set_label(ioArray[last_Write]->mapfile, temp_label,*posn-1 );

  free(temp_label);
}

/**
 * mttcpy:
 * @label:     (I)      new label
 *
 * description: Copy all labels from the current input map to the header (struct)
 * of the current output mapfile, adding @label to the end.
 */
FORTRAN_SUBR( MTTCPY, mttcpy,
	      (const fpstr label, int label_len),
	      (const fpstr label),
	      (const fpstr label, int label_len))
{
  char *temp_label;
  int nlabel=0,i;

  /* no output map, nothing to do */
  if (last_Write == -1) return;

  temp_label = ccp4_FtoCString(FTN_STR(label), FTN_LEN(label));

  if (last_Read != -1) {
    nlabel = ccp4_cmap_number_label(ioArray[last_Read]->mapfile);
    for (i=0 ; i != nlabel ; ++i) 
      ccp4_cmap_set_label(ioArray[last_Write]->mapfile,
			ccp4_cmap_get_label(ioArray[last_Read]->mapfile,i),
			i);
  }

  ccp4_cmap_set_label(ioArray[last_Write]->mapfile, temp_label,nlabel);

  free(temp_label);
}

FORTRAN_SUBR( CCP4_MAP_COPY_TITLE, 
	      ccp4_map_copy_title,
	      (const fpstr label, int label_len),
	      (const fpstr label),
	      (const fpstr label, int label_len))
     /* see MTTCPY */
{
  char *temp_label;
  int nlabel=0,i;

  /* no output map, nothing to do */
  if (last_Write == -1) return;

  temp_label = ccp4_FtoCString(FTN_STR(label), FTN_LEN(label));
  
  if (last_Read != -1) {
    nlabel = ccp4_cmap_number_label(ioArray[last_Read]->mapfile);
    for (i=0 ; i != nlabel ; ++i) 
      ccp4_cmap_set_label(ioArray[last_Write]->mapfile,
			ccp4_cmap_get_label(ioArray[last_Read]->mapfile,i),
			i);
  }

  ccp4_cmap_set_label(ioArray[last_Write]->mapfile, temp_label,nlabel);

  free(temp_label);
}


/**
 * mskput:
 * @skew_rot:  (I)      skew rotation matrix S[3][3];S11, S12, etc
 * @skew_trans (I)     skew translation vector T[3]
 *
 * description: Put skew transformation in header struct of current output
 * mapfile.
 * Note: Phil Evans (9/93), should probably not use this routine.
 */
FORTRAN_SUBR( MSKPUT, mskput,
	      (float *skew_rot, float *skew_trans),
	      (float *skew_rot, float *skew_trans),
	      (float *skew_rot, float *skew_trans))

{
  ccp4_cmap_set_mask(ioArray[last_Write]->mapfile, skew_rot, skew_trans);
}

FORTRAN_SUBR( CCP4_MAP_WRITE_SKEW_INFO, 
	      ccp4_map_write_skew_info,
	      (float *skew_rot, float *skew_trans),
	      (float *skew_rot, float *skew_trans),
	      (float *skew_rot, float *skew_trans))
     /* see MSKPUT */
{
  ccp4_cmap_set_mask(ioArray[last_Write]->mapfile, skew_rot, skew_trans);
}

/**
 * mskget:
 * @skew_rot:  (I)      skew rotation matrix S[3][3]; S11, S12, etc
 * @skew_trans (I)     skew translation vector T[3]
 * returns: value of skew_flag.
 *
 * description: Get skew translation from input header struct.
 * Skew transformation form orthogonal atom fram to orthogonal map frame:
 *       Xo(map) = S*(Xo(atoms)-T)
 */
FORTRAN_FUN(int, MSKGET, mskget,
	    (float *mask_rot, float *mask_trans),
	    (float *mask_rot, float *mask_trans),
	    (float *mask_rot, float *mask_trans))

{
  return ccp4_cmap_get_mask(ioArray[last_Read]->mapfile, mask_rot, mask_trans);
}

/**
 * ccp4_map_write_section_header:
 * @iunit:     (I)      stream number
 * @section:   (I)      array holding the map section.
 * @local_hdr: (I)      local header information
 *
 * description: Write whole map section to stream @iunit.
 * This routine is only suitable if the whole map section is to
 * be written.
 */
FORTRAN_SUBR( CCP4_MAP_WRITE_SECTION_HEADER,
              ccp4_map_write_section_header,
              (int *iunit, float *section, const fpstr local_hdr, int local_hdr_len),
              (int *iunit, float *section, const fpstr local_hdr),
              (int *iunit, float *section, const fpstr local_hdr, int local_hdr_len))
{
  int ii;
  char *tmp_label;

  if ( (ii = GetChannel(*iunit)) == MAXFILES ||
       !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"CCP4_MAP_WRITE_SECTION_HEADER", NULL);

  if (ccp4_cmap_write_section(ioArray[ii]->mapfile, section) == EOF) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_WriteFail),
		"CCP4_MAP_WRITE_SECTION_HEADER", NULL);

  tmp_label = ccp4_FtoCString(FTN_STR(local_hdr), FTN_LEN(local_hdr) );

  if (ccp4_cmap_write_section_header(ioArray[ii]->mapfile, tmp_label) == EOF) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_WriteFail),
		"CCP4_MAP_WRITE_SECTION_HEADER", NULL);
  free(tmp_label);
  /* record for FORTRAN API */
  last_Write = ii;
}

/**
 * ccp4_map_read_section_header:
 * @iunit:     (I)      map stream number
 * @section:   (O)      array containing section of map data read from file
 * @local_hdr: (O)      local header information as CHAR array
 * @ier:       (O)      error code =0, OK
 *                                 =non-zero, error or end-of-file
 *
 * description: Read next whole map section from stream @iunit to array @section.
 *     @section is returned in the same mode as in the file, no data
 *     conversion is performed.
 */
FORTRAN_SUBR( CCP4_MAP_READ_SECTION_HEADER, 
              ccp4_map_read_section_header,
	      (int *iunit, float *section, const fpstr local_hdr, 
               int *ier, int local_hdr_len),
	      (int *iunit, float *section, const fpstr local_hdr, 
               int *ier),
	      (int *iunit, float *section, const fpstr local_hdr, 
               int local_hdr_len, int *ier))
{
  int ii;

  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"CCP4_MAP_READ_SECTION_HEADER", NULL);

  if ((*ier = ccp4_cmap_read_section(ioArray[ii]->mapfile,section)) == -1) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_ReadFail),
	       "CCP4_MAP_READ_SECTION_HEADER", NULL);
    return; }

  if ((*ier = ccp4_cmap_read_section_header(ioArray[ii]->mapfile, FTN_STR(local_hdr)))
      == EOF)
    ccp4_signal(CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_ReadFail),
	       "CCP4_MAP_READ_SECTION_HEADER", NULL);

  *ier = (*ier >= 0) ? 0 : -1;

  /* record for FORTRAN API */
  last_Read = ii;
}


FORTRAN_SUBR( CCP4_MAP_SET_LOCAL_HEADER,
              ccp4_map_set_local_header,
              (int *iunit, int *size),
              (int *iunit, int *size),
              (int *iunit, int *size))
{
  int ii;
  
  if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"CCP4_MAP_SET_LOCAL_HEADER", NULL);

   ccp4_cmap_set_local_header(ioArray[ii]->mapfile, (size_t)*size);
}

FORTRAN_SUBR( CCP4_MAP_GET_LOCAL_HEADER,
              ccp4_map_get_local_header,
              (int *iunit, int *size),
              (int *iunit, int *size),
              (int *iunit, int *size))
{
  int ii;
  
    if ( (ii = GetChannel(*iunit)) == MAXFILES || !ioArray[ii]->mapfile) 
    ccp4_signal(CCP4_ERRLEVEL(4) | CMAP_ERRNO(CMERR_NoChannel),
		"CCP4_MAP_GET_LOCAL_HEADER", NULL);

   *size = (int) ccp4_cmap_get_local_header(ioArray[ii]->mapfile);
}

/**
 * mwfnam:
 * @logname:   (O)      logical name of last open file
 *
 * description: Returns the logical name for the last written to file
 */
FORTRAN_SUBR( MWFNAM, mwfnam, 
	      (fpstr logname, int logname_len), 
	      (fpstr logname), 
	      (fpstr logname, int logname_len))
{
  strncpy(FTN_STR(logname), ioArray[last_Write]->logname, 
	  MIN(strlen(ioArray[last_Write]->logname),
	      FTN_LEN(logname)));
}

FORTRAN_SUBR( CCP4_MAP_GET_LAST_WRITE_FILENAME, 
	      ccp4_map_get_last_write_filename, 
	      (fpstr logname, int logname_len), 
	      (fpstr logname), 
	      (fpstr logname, int logname_len))
     /*see MWFNAM */
{
  strncpy(FTN_STR(logname), ioArray[last_Write]->logname, 
	  MIN(strlen(ioArray[last_Write]->logname),
	      FTN_LEN(logname)));
}

/**
 * mwfnam:
 * @logname:   (O)      logical name of last open file
 *
 * description: Returns the logical name for the last read from file
 */
FORTRAN_SUBR( MRFNAM, mrfnam, 
	      (fpstr logname, int logname_len), 
	      (fpstr logname), 
	      (fpstr logname, int logname_len))
{
  strncpy(FTN_STR(logname), ioArray[last_Read]->logname, 
	  MIN(strlen(ioArray[last_Read]->logname),
	      FTN_LEN(logname)));
}

FORTRAN_SUBR( CCP4_MAP_GET_LAST_READ_FILENAME,
	      ccp4_map_get_last_read_filename,
	      (fpstr logname, int logname_len), 
	      (fpstr logname), 
	      (fpstr logname, int logname_len))
     /* see MRFNAM */
{
  strncpy(FTN_STR(logname), ioArray[last_Read]->logname, 
	  MIN(strlen(ioArray[last_Read]->logname),
	      FTN_LEN(logname)));
}

/**
 * mstmst:
 * @map:       (O)      "MAP "
 *
 * description: Set integer MAP to string "MAP "
 */
FORTRAN_SUBR( MSTMST, mstmst, (int *map), (int *map), (int *map))
{
  strncpy((char *)map, "MAP ", strlen("MAP "));
}

/**
 * modecv:
 * @output:    (O)        output array contain translated reals
 * @input:     (I)        pre-translation array (read from file)
 * @number:    (I)        number of elements to be translated
 * @mode:      (I)        storage mode, = 0 FORTRAN character (uint8)
 *                                      = 1 FORTRAN I*2 (uint16) 
 *					= 3 complex short
 *					= 4 complex real (uint32)
 *
 * description: Convert @number items form @input (mode @mode) to float in @output.
 */
FORTRAN_SUBR( MODECV, modecv,
	      (float *output, float *input, int *number, int *mode),
	      (float *output, float *input, int *number, int *mode),
	      (float *output, float *input, int *number, int *mode))
{
  if( (ccp4_utils_translate_mode_float(output, (void *)input, *number, *mode)
       != *number))
     ccp4_signal(CMAP_ERRNO(CMERR_ParamError),"MODECV", NULL);
}

FORTRAN_SUBR( CCP4_MAP_MODE_TO_REAL, 
	      ccp4_map_mode_to_real,
	      (float *output, float *input, int *number, int *mode),
	      (float *output, float *input, int *number, int *mode),
	      (float *output, float *input, int *number, int *mode))
     /* see MODECV */
{
  if( (ccp4_utils_translate_mode_float(output, (void *)input, *number, *mode)
       != *number))
     ccp4_signal(CMAP_ERRNO(CMERR_ParamError),"MODECV", NULL);
}

FORTRAN_FUN(int, NBYTXX, nbytxx,  (int *nwords), (int *nwords),(int *nwords)) 
/*
      Returns the number of machine items in NWORDS words (floats)
*/
{
  return ((int) sizeof(float)) * (*nwords);
}
