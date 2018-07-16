/*
     ccp4_parser_f.c: Fortran API to ccp4_parser.c
     Copyright (C) 2001  CCLRC, Peter Briggs

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

/** @page cparser_f_page Fortran API to CParser
 *
 *  @section cparser_f_file_list File list

<ul>
<li>ccp4_parser_f.c
</ul>
 *
 *  @section cparser_f_overview Overview

This library consists of a set of wrappers to the CParser library
giving the same API as the original parser.f

 */

/** @file ccp4_parser_f.c
 *
 *  @brief Fortran API to ccp4_parser.c.
 *
 *  @author Peter Briggs
 */

/*   ccp4_parser_f.c
     Peter Briggs CCP4 August 2001

     This file contains Fortran API functions which provide the following
     functionality originally found in the parser.f file:

     PARSER  read and interpret data from the input stream
     PARSE   free format read routine
     PARSEDL change delimiters

     plus internal utility functions which should not be accessed
     directly from an application.
*/

/*#define FORTRAN_CALL_DEBUG 1*/

#if defined (FORTRAN_CALL_DEBUG)
#  define PARSER_DEBUG(x) x
#else
#  define PARSER_DEBUG(x)
#endif

/* CCP4 header file for machine dependency */
#include "../ccp4/ccp4_sysdep.h"

/* CCP4 header file for Fortran interfaces */
#include "../ccp4/ccp4_fortran.h"

/* C parser header file */
#include "../ccp4/ccp4_parser.h"
#include "../ccp4/ccp4_general.h"

/* rcsid[] = "$Id$" */

/*------------------------------------------------------------------*/

/* Utility functions */     

/*------------------------------------------------------------------*/

/* fparse_isblank

   Given a string and its length, returns 1 if the only characters are
   whitespace or newlines, and 0 if there are other characters.
*/
int fparse_isblank(const char *line, const int line_len)
{
  int i;
  if (line_len > 0) {
    for (i=0; i<line_len; i++) {
      if (!charmatch(line[i]," \t\n")) return 0;
    }
    return 1;
  }
  return 0;
}     
/*------------------------------------------------------------------*/

/* fparse_strncpypad

   Copy lfstr characters from C-style string cstr into Fortran-style string
   fstr. lfstr is the maximum length of the Fortran-style string.

   If cstr contains fewer characters than lfstr then fstr is padded with
   spaces. Note that fstr will not contain an end-of-string null character,
   so printf using %s may also show trailing garbage.
 */
int fparse_strncpypad(char *fstr, const char *cstr, const int lfstr)
{
  int i;

  if (!cstr || !fstr || !lfstr) return 0;

  strncpy(fstr,cstr,lfstr);
  if (strlen(cstr) < lfstr) {
    for (i=strlen(cstr); i<lfstr; i++) fstr[i] = ' ';
  }
  return 1;
}     

/*------------------------------------------------------------------*/

/* fparse_populate_arrays

   Given a pointer to a CCP4PARSERARRAY, will populate the equivalent
   arrays supplied by the Fortran application.

   Should only be called from PARSER or PARSE.
 */
int fparse_populate_arrays(CCP4PARSERARRAY *parser, int *ibeg, int *iend,
			   int *ityp, float *fvalue, fpstr cvalue, int cvalue_len,
			   int *idec)
{
  int i,iindex;
  char *pstr;
  CCP4PARSERTOKEN *token = NULL;

  /* Check we have something to work with */
  if (!parser) return 0;

  /* Set a convenient pointer to the array of tokens */
  token = parser->token;

  /* Loop over all tokens */
  for (i=0; i<parser->ntokens; i++) {
    PARSER_DEBUG(printf("fparse_populate_arrays: Token %d\n",i);)

    /* Positions for first and last characters
       Add one on because Fortran numbers array elements from one */
    ibeg[i] = token[i].ibeg+1;
    iend[i] = token[i].iend+1;
    PARSER_DEBUG(printf("fparse_populate_arrays: ibeg = %d, iend = %d\n",ibeg[i],iend[i]);)

    /* Token types
       0 = null field
       1 = character string
       2 = number
    */
    if (!token[i].isnull) {
      /* Either string or number */
      
      if (token[i].isstring) {
	/* String token */
	PARSER_DEBUG(printf("fparse_populate_arrays: string token\n");)
	ityp[i] = 1;
	PARSER_DEBUG(printf("fparse_populate_arrays: 4-letter version is \"%s\"\n",token[i].word);)
	if (token[i].strlength > 4) {
	  idec[i] = 4;
	} else {
	  idec[i] = token[i].strlength;
	}
	PARSER_DEBUG(printf(" ityp = %d\n",ityp[i]);)
	
      } else if (token[i].isnumber) {
	/* Numerical token */
	PARSER_DEBUG(printf("fparse_populate_arrays: number token\n");)
	ityp[i] = 2;
	fvalue[i] = (float) token[i].value;
	if (token[i].frcdigits) {
	  /* Real number */
	  idec[i] = (token[i].intdigits+1)*100 + token[i].frcdigits;
	} else {
	  /* Integer number */
	  idec[i] = token[i].intdigits;
	}
	PARSER_DEBUG(printf(" ityp = %d\n",ityp[i]);)
      } else {
	/* Unrecognised token type */
	PARSER_DEBUG(printf("fparse_populate_arrays: unrecognised token type - setting to null");)
	ityp[i] = 0;
      }
      
      /* cvalue is set for both string and numerical tokens
	 It appears that Fortran allocates a block of memory for
	 string arrays which is basically string(1) followed by
	 string(2) etc, and then the length string_len is the
	 number of characters in a single element.
      */
      iindex = i*cvalue_len;
      pstr   = &(cvalue[iindex]);
      PARSER_DEBUG({
	printf("fparse_populate_arrays: token string stored as \"%s\"\n",token[i].word);
	printf("fparse_populate_arrays: initial value of cvalue[%d] is \"%s\"\n",i,pstr);
      })
      /* Store the value and pad if necessary */
      fparse_strncpypad(pstr,token[i].word,cvalue_len);
      PARSER_DEBUG(printf("fparse_populate_arrays: cvalue[%d] is now \"%s\"\n",i,pstr);)
      
    } else {
      /* Null token */
      PARSER_DEBUG(printf("fparse_populate_arrays: null field");)
      ityp[i] = 0;
      PARSER_DEBUG(printf(" ityp = %d\n",ityp[i]);)
    }
    /* All assignments complete for this token */
  }
  return 1;
}     


/*------------------------------------------------------------------*/

/* fparse_delimiters

   Set and get delimiters for use in ccp4_parse(r). This should only
   be called internally (from PARSDL).

   If called with parser set to NULL, the lists of characters in
   new_delimiters and new_nulldelimiters are stored.

   When called with a valid pointer to a parser structure, the stored
   lists of characters form the new delimiter sets via a call to
   ccp4_parse_delimiters.

   NB: if delimiter strings are set with a call to fparse_delimiters
   then they persist in memory unless a subsequent call is made with
   all arguments set to NULL. This represents a potential memory leak.
*/
int fparse_delimiters(CCP4PARSERARRAY *parser, char *new_delimiters,
		      char *new_nulldelimiters)
{
  int ndelim;
  static char *delimiters=NULL, *nulldelimiters=NULL;

  PARSER_DEBUG(puts("fparse_delimiters: starting");)

  if (!parser) {

    PARSER_DEBUG({
      puts("fparse_delimiters: setting local lists");
      if (new_delimiters) printf("fparse_delimiters: new_delimiters = \"%s\"\n",new_delimiters);
      if (new_nulldelimiters) printf("fparse_delimiters: new_nulldelimiters = \"%s\"\n",new_nulldelimiters);
    })
    /* Set the local lists of delimiters
       depending on the input */
    /* Standard delimiters */
    if (!new_delimiters) {
      if (delimiters) free(delimiters);
      delimiters = NULL;
    } else {
      ndelim = strlen(new_delimiters) + 1;
      delimiters = (char *) realloc(delimiters,ndelim*sizeof(char));
      if (delimiters) {
	strncpy(delimiters,new_delimiters,ndelim);
	PARSER_DEBUG(printf("fparse_delimiters: delimiters set to \"%s\"\n",delimiters);)
      } else {
	ccperror(4,"fparse_delimiters: couldn't reallocate delimiters");
      }
    }
    /* Null delimiters */
    if (!new_nulldelimiters) {
      if (nulldelimiters) free(nulldelimiters);
      nulldelimiters = NULL;
    } else {
      ndelim = strlen(new_nulldelimiters) + 1;
      nulldelimiters = (char *) realloc(nulldelimiters,ndelim*sizeof(char));
      if (nulldelimiters) {
	strncpy(nulldelimiters,new_nulldelimiters,ndelim);
	PARSER_DEBUG(printf("fparse_delimiters: nulldelimiters set to \"%s\"\n",nulldelimiters);)
      } else {
	ccperror(4,"fparse_delimiters: couldn't reallocate null delimiters");
      }
    }
  } else {
    /* Set the parser array so that it uses the delimiters
       already stored */

    PARSER_DEBUG({
      puts("fparse_delimiters: setting parser to use stored delimiters");
    })
    if (!ccp4_parse_delimiters(parser,delimiters,nulldelimiters))
        ccperror(4,"fparse_delimiters: couldn't reset delimiters");
    PARSER_DEBUG(printf("fparse_delimiters: now set to \"%s\" and \"%s\"\n",parser->delim, parser->nulldelim);)
    return 0;
  }
  return 1;
}     

/*------------------------------------------------------------------*/

/* PARSER

   This function implements the Fortran equivalent of:

   SUBROUTINE PARSER(KEY,LINE,IBEG,IEND,ITYP,FVALUE,CVALUE,IDEC,NTOK,LEND,PRINT)

   which is the original CCP4 keyword parser.
   The function wraps the ccp4_parser routines to mimic the behaviour of the
   original subroutine.
 */

FORTRAN_SUBR(PARSER,parser,
             (fpstr key, fpstr line, int *ibeg, int *iend, int *ityp, float *fvalue,
              fpstr cvalue, int *idec, int *ntok, ftn_logical *lend, 
              const ftn_logical *print, int key_len, int line_len, int cvalue_len),
             (fpstr key, fpstr line, int *ibeg, int *iend, int *ityp, float *fvalue,
              fpstr cvalue, int *idec, int *ntok, ftn_logical *lend, 
              const ftn_logical *print),
             (fpstr key, int key_len, fpstr line, int line_len, int *ibeg, int *iend,
              int *ityp, float *fvalue, fpstr cvalue, int cvalue_len, int *idec,
              int *ntok, ftn_logical *lend, const ftn_logical *print))
{
  static FILE *fparse_fp = NULL;
  int  max_line_len,lline,cprint = 0;
  char *cline;
  CCP4PARSERARRAY *parser = NULL;

  PARSER_DEBUG(puts("PARSER: starting");)

  /* On input ntok is the maximum number of fields to be parsed
     If ntok < 20 then set to 20 */
  if (*ntok < 20) *ntok = 20;
  /* also reset if ntok is silly - this is a trap for uninitialised ntok */
  if (*ntok > 10000) *ntok = 10000;

  PARSER_DEBUG({
    printf("PARSER: set maximum number of tokens to %d\n",*ntok);
    printf("PARSER: line is initially \"%s\"\n",line);
  })

  /* copy to C string, but also allocate enough memory for return value */
  cline  = (char *) ccp4_utils_malloc((FTN_LEN(line)+1)*sizeof(char));
  lline = ccp4_utils_flength(FTN_STR(line),FTN_LEN(line));
  strncpy(cline,FTN_STR(line),lline);
  cline[lline] = '\0';
  /* Get the maximum line length
     Use the FTN_LEN macro for this since we can't rely on
     getting the line_len argument on all systems */
  max_line_len = FTN_LEN(line);
  PARSER_DEBUG(printf("PARSER: line length is %d\n",max_line_len);)

  /* Set up a parser array to handle ntok tokens */
  parser = (CCP4PARSERARRAY *) ccp4_parse_start(*ntok);
  if (!parser) {
    PARSER_DEBUG(printf("PARSER: failed to allocate memory");)
    *lend = FORTRAN_LOGICAL_TRUE;
    return;
  }
  PARSER_DEBUG(puts("PARSER: parser array initialised");)

  /* Set up the delimiters */
    PARSER_DEBUG(puts("PARSER: fetching delimiters");)
  fparse_delimiters(parser,NULL,NULL);
  PARSER_DEBUG(puts("PARSER: delimiters set");)

  /* Set up the maximum and minimum exponents used to avoid
     under/overflow
     Since Fortran REALs are equivalent to C floats, make sure
     these limits are appropriate for floats */
  ccp4_parse_maxmin(parser,FLT_MAX_10_EXP,FLT_MIN_10_EXP);

  /* Silent or verbose output? */
    PARSER_DEBUG(printf("PARSER: print set to %d\n",*print);)
  if (*print != FORTRAN_LOGICAL_FALSE) cprint = 1;

  PARSER_DEBUG({ 
    if (cprint) {
      puts("PARSER: verbose output ON");
    } else {
      puts("PARSER: verbose output OFF");
    }
  })

  /* Was ccp4_parser reading from an external file last time? */
  if (fparse_fp) {
    PARSER_DEBUG(printf("PARSER: we were reading from an external file\n");)
    parser->fp = fparse_fp;
  }

  /* Call ccp4_parser to do the work */
  PARSER_DEBUG({
    printf("PARSER: line sent as \"%s\"\n",cline);
  })
  *ntok = ccp4_parser(cline,max_line_len,parser,cprint);
  PARSER_DEBUG({
    printf("PARSER: returned %d tokens from ccp4_parser\n",*ntok);
    printf("PARSER: line returned as \"%s\"\n",cline);
  })

  /* Check for end-of-file */
  if (!*ntok) {
    PARSER_DEBUG(puts("PARSER: reached end of file");)
    *lend = FORTRAN_LOGICAL_TRUE;
  } else {
    PARSER_DEBUG(puts("PARSER: end of file not reached yet");)
    *lend = FORTRAN_LOGICAL_FALSE;

    /* Keyword
       NB You need to pad the string with spaces before sending back to
       Fortran
    */
    PARSER_DEBUG(printf("PARSER: ccp4_parser keyword is \"%s\"\n",parser->keyword);)
    ccp4_CtoFString(FTN_STR(key),FTN_LEN(key),parser->keyword);
    PARSER_DEBUG(printf("PARSER: PARSER keyword is \"%s\"\n",key);)

    PARSER_DEBUG(printf("PARSER: ccp4_parser line is \"%s\"\n",cline);)
    ccp4_CtoFString(FTN_STR(line),FTN_LEN(line),cline);
    PARSER_DEBUG(printf("PARSER: PARSER line is \"%s\"\n",line);)

    /* Populate the Fortranic arrays */
    PARSER_DEBUG(printf("PARSER: about to populate the arrays\n");)
    fparse_populate_arrays(parser,ibeg,iend,ityp,fvalue,cvalue,FTN_LEN(cvalue),idec);
    PARSER_DEBUG(printf("PARSER: arrays populated\n");)
  }
  PARSER_DEBUG(puts("PARSER: finished assignments");)

  /* Check if ccp4_parser was reading from an external file
     which is still open */
  if (parser->fp) {
    PARSER_DEBUG(printf("PARSER: reading from external file\n");)
    fparse_fp = parser->fp;
  } else {
    fparse_fp = NULL;
  }

  /* Free the parser array */
  ccp4_parse_end(parser);
  PARSER_DEBUG(puts("PARSER: freed the memory, returning");)

  free(cline);

  return;
}     

/*------------------------------------------------------------------*/

/* PARSE

   This function implements the Fortran equivalent of:

   SUBROUTINE PARSE(LINE,IBEG,IEND,ITYP,FVALUE,CVALUE,IDEC,N)
*/
FORTRAN_SUBR(PARSE,parse,
	     (fpstr line, int *ibeg, int *iend, int *ityp, float *fvalue,
	      fpstr cvalue, int *idec, int *n, int line_len, int cvalue_len),
	     (fpstr line, int *ibeg, int *iend, int *ityp, float *fvalue,
	      fpstr cvalue, int *idec, int *n),
	     (fpstr line, int line_len, int *ibeg, int *iend,
	      int *ityp, float *fvalue, fpstr cvalue, int cvalue_len, int *idec,
	      int *n))
{
  char *temp_line;
  static int maxtok = 0;
  CCP4PARSERARRAY *parser = NULL;

  temp_line = ccp4_FtoCString(FTN_STR(line), FTN_LEN(line));

  PARSER_DEBUG({
    puts("PARSE: starting");
    printf("PARSE: line is initially \"%s\"\n",temp_line);
  })

  /* Set a value for the maximum number of tokens

     When called for the first time n should be less than zero, and
     ABS(n) is the maximum number of tokens.
     When called subsequently, if n > 0 then this is the number of
     tokens already read in.

     To handle this, set maxtok only when a -ve value of n is input.
  */
  if (*n < 0) {
    maxtok = -1*(*n);
  }
  /* If ntok is still zero then this is probably an error so return
     without any action */
  if (maxtok == 0) {
    printf("PARSE: zero number of tokens specified - aborting\n");
    return;
  }

  /* Set up a parser array to handle maxtok tokens */
  parser = (CCP4PARSERARRAY *) ccp4_parse_start(maxtok);
  if (!parser) {
    PARSER_DEBUG(printf("PARSE: failed to allocate memory");)
    return;
  }
  PARSER_DEBUG(puts("PARSE: parser array initialised");)

  /* Set up the delimiters */
  PARSER_DEBUG(puts("PARSER: fetching delimiters");)
  fparse_delimiters(parser,NULL,NULL);
  PARSER_DEBUG(puts("PARSER: delimiters set");)
  
  /* Call ccp4_parse to do the work */
  *n = ccp4_parse(temp_line,parser);

  PARSER_DEBUG(printf("PARSE: returned %d tokens from ccp4_parse\n",*n);)

  /* Populate the Fortranic arrays */
  PARSER_DEBUG(printf("PARSE: about to populate the arrays\n");)
  fparse_populate_arrays(parser,ibeg,iend,ityp,fvalue,cvalue,FTN_LEN(cvalue),idec);
  PARSER_DEBUG(printf("PARSE: arrays populated\n");)
  
  /* Free the parser array */
  ccp4_parse_end(parser);
  PARSER_DEBUG(puts("PARSE: freed the memory, returning");)

  free(temp_line);

  return;
}

/*------------------------------------------------------------------*/

/* PARSDL

   Set delimiters for PARSE, PARSER

   SUBROUTINE PARSDL(NEWDLM,NNEWDL,NSPECD)
*/
FORTRAN_SUBR(PARSDL,parsdl,
	     (fpstr newdlm, int *nnewdl, int *nspecd, int newdlm_len),
	     (fpstr newdlm, int *nnewdl, int *nspecd),
	     (fpstr newdlm, int newdlm_len, int *nnewdl, int *nspecd))
{
  int nnulldl;
  char *delim=NULL, *nulldelim=NULL;

  PARSER_DEBUG(puts("PARSDL: starting");)

  if (*nnewdl <= 0) {

    /* Reset the delimiters to the defaults */
    PARSER_DEBUG(puts("PARSDL: resetting to default delimiters");)
    fparse_delimiters(NULL,NULL,NULL);
  } else {
    /* Check the number of new delimiters is no more than the
       number of characters actually supplied */
    if (FTN_LEN(newdlm) < *nnewdl) {
      ccperror(4,"PARSDL: too few delimiter characters supplied");
      return;
    }
    /* Extract the list of delimiter characters */
    delim = (char *) malloc((*nnewdl+1)*sizeof(char));
    if (delim) {
      strncpy(delim,newdlm,*nnewdl);
      delim[*nnewdl] = '\0';
    }
    PARSER_DEBUG(printf("PARSDL: delimiters are \"%s\"\n",delim);)
    /* Extract the list of null delimiter characters (if any) */
    nnulldl = *nnewdl - *nspecd;
    if (nnulldl > 0) {
      nulldelim = (char *) malloc((nnulldl+1)*sizeof(char));
      if (nulldelim) {
	strncpy(nulldelim,&(newdlm[*nspecd]),nnulldl);
	nulldelim[nnulldl] = '\0';
      }
      PARSER_DEBUG(printf("PARSDL: null delimiters are \"%s\"\n",nulldelim);)
    } else {
      nulldelim = (char *) malloc(sizeof(char));
      if (nulldelim) {
	nulldelim[0] = '\0';
      }
    }
    /* Store the delimiters for later reference */
    fparse_delimiters(NULL,delim,nulldelim);
  }
  /* Free locally allocated memory */
  if (delim) free(delim);
  if (nulldelim) free(nulldelim);

  PARSER_DEBUG(puts("PARSDL: finished");)
  return;
}

