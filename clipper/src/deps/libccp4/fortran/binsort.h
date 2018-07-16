/*
     binsort.h: header for binary sorting functions

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
/****************************************************************************
  binsort.h
  Z130891

binsort - key data types definition.
****************************************************************************/

/*** Data types definition name begining with U means unsigned... ***/

#define CHAR        (int)1
#define UCHAR       (int)2
#define SHORT       (int)3
#define USHORT      (int)4
#define LONG        (int)5
#define ULONG       (int)6
#define FLOAT       (int)7
#define DOUBLE      (int)8

/*** Sorting order ***/
#define ASC         (int)0      /* ascending */
#define DESC        (int)1      /* descending */
