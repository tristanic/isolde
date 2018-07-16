//  $Id: mmdb_tables.cpp $
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
//    24.07.15   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :   MMDB_Tables <implementation>
//       ~~~~~~~~~
//  **** Project :   MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//
//  **** Namespace :  mmdb::
//
//  **** Functions :
//       ~~~~~~~~~~~
//
//  **** Constants :  AName  ( array of 2-character atom names       )
//       ~~~~~~~~~~~  HAName ( array of 2=character heteroatom names )
//                    RName  ( 3-characters amino acid names         )
//                    RName1 ( 1-characters amino acid names         )
//
//
//  (C) E. Krissinel  2000-2015
//
//  =================================================================
//

#include <string.h>

#include "mmdb_tables.h"
#include "mmdb_defs.h"

namespace mmdb  {

  //  ===============================================================

  cpstr const ElementName[nElementNames] = {
    " H", "HE",
    "LI", "BE", " B", " C", " N", " O", " F", "NE",
    "NA", "MG", "AL", "SI", " P", " S", "CL", "AR",
    " K", "CA",
          "SC", "TI", " V", "CR", "MN", "FE",
          "CO", "NI", "CU", "ZN",
                "GA", "GE", "AS", "SE", "BR", "KR",
    "RB", "SR",
          " Y", "ZR", "NB", "MO", "TC", "RU",
          "RH", "PD", "AG", "CD",
                "IN", "SN", "SB", "TE", " I", "XE",
    "CS", "BA",
          "LA", "CE", "PR", "ND", "PM", "SM", "EU",
          "GD", "TB", "DY", "HO", "ER", "TM", "YB",
                "LU", "HF", "TA", " W", "RE", "OS",
                "IR", "PT", "AU", "HG",
                      "TL", "PB", "BI", "PO", "AT", "RN",
    "FR", "RA",
          "AC", "TH", "PA", " U", "NP", "PU", "AM",
          "CM", "BK", "CF", "ES", "FM", "MD", "NO",
                "LR", "RF", "DB", "SG", "BH", "HS",
                "MT", "UN", "UU", "UB",
                            "UQ",       "UH",       "UO",
    " D", "AN"
  };

  realtype const MolecWeight[nElementNames] = {
    1.0079, 4.0026,
    6.9410, 9.0122, 10.811, 12.011, 14.007, 15.999, 18.998, 20.180,
    22.990, 24.305, 26.982, 28.086, 30.974, 32.066, 35.453, 39.948,
    39.098, 40.078,
            44.956, 47.867, 50.942, 51.996, 54.938, 55.845,
            58.993, 58.693, 63.546, 65.390,
                    69.723, 72.610, 74.922, 78.960, 79.904, 83.800,
    85.468, 87.620,
            88.906, 91.224, 92.906, 95.940, 97.907, 101.07,
            102.91, 106.42, 107.87, 112.41,
                    114.82, 118.71, 121.76, 127.60, 126.90, 131.29,
    132.91, 137.33,
            138.91, 140.12, 140.91, 144.24, 144.91, 150.36, 151.96,
            157.25, 158.93, 162.50, 164.93, 167.26, 168.93, 173.04,
                    174.97, 178.49, 180.95, 183.84, 186.21, 190.23,
                    192.22, 195.08, 196.97, 200.59,
                            204.38, 207.20, 208.98, 208.98, 209.99, 222.02,
    232.02, 226.03,
            227.03, 232.04, 231.04, 238.03, 237.05, 244.06, 243.06,
            247.07, 247.07, 251.08, 252.08, 257.10, 258.10, 259.10,
                    262.11, 263.11, 262.11, 266.12, 264.12, 269.13,
                    268.14, 272.15, 272.15, 277.00,
                            289.00, 289.00, 293.00,
    2.0200, 3.0300
  };


  realtype const CovalentRadius[nElementNames] = {
    0.32, 0.93,
    1.23, 0.90, 0.82, 0.77, 0.75, 0.73, 0.72, 0.71,
    1.54, 1.36, 1.18, 1.11, 1.06, 1.02, 0.99, 0.98,
    2.03, 1.91,
          1.62, 1.45, 1.34, 1.18, 1.17, 1.17,
          1.16, 1.15, 1.17, 1.25,
                1.26, 1.22, 1.20, 1.16, 1.14, 1.12,
    2.16, 1.91,
          1.62, 1.45, 1.34, 1.30, 1.27, 1.25,
          1.25, 1.28, 1.34, 1.48,
                1.44, 1.41, 1.40, 1.36, 1.33, 1.31,
    2.35, 1.98,
          1.69, 1.44, 1.34, 1.30, 1.28, 1.26, 1.27,
          1.30, 1.34, 1.49, 1.48, 1.47, 1.46, 1.46,
                1.45, 1.43, 2.50, 2.40, 2.20, 1.65,
                1.65, 1.64, 1.63, 1.62,
                      1.85, 1.61, 1.59, 1.59, 1.58, 1.57,
    1.56, 1.74,
          1.56, 1.65, 1.65, 1.42, 1.65, 1.65, 1.65,
          1.65, 1.65, 1.65, 1.65, 1.65, 1.65, 1.65,
                1.65, 0.32, 0.10, /**/
                                  0.20, 0.20, 0.20,
                0.20, 0.20, 0.20, 0.20,
                      0.20, 0.20, 0.20,
    0.32, 0.32
  };


  realtype const VdWaalsRadius[nElementNames] = {
    1.20, 1.40,
    1.82, 1.78, 1.74, 1.70, 1.55, 1.52, 1.47, 1.54,
  //      ^^^^  ^^^^   <- only a guess
    2.27, 1.73, 1.80, 2.10, 1.80, 1.80, 1.75, 1.88,
  //            ^^^^
    2.75, 2.65,
  //      ^^^^
          2.55, 2.45, 2.35, 2.20, 1.73, 1.90,
  //      ^^^^  ^^^^  ^^^^  ^^^^  ^^^^  ^^^^
          1.75, 1.63, 1.40, 1.39,
  //      ^^^^
                1.87, 1.86, 1.85, 1.90, 1.85, 2.02,
  //                  ^^^^
    2.75, 2.65,
  //^^^^  ^^^^
          2.55, 2.45, 2.35, 2.20, 2.05, 1.90,
  //      ^^^^  ^^^^  ^^^^  ^^^^  ^^^^  ^^^^
          1.75, 1.63, 1.72, 1.58,
  //      ^^^^
                1.93, 2.17, 2.10, 2.06, 1.98, 2.16,
  //                        ^^^^
    2.75, 2.75,
  //^^^^  ^^^^
          2.75, 2.75, 2.75, 2.75, 2.75, 2.75, 2.75,
  //      ^^^^  ^^^^  ^^^^  ^^^^  ^^^^  ^^^^  ^^^^
          2.75, 2.75, 2.75, 2.75, 2.75, 2.65, 2.55,
  //      ^^^^  ^^^^  ^^^^  ^^^^  ^^^^  ^^^^  ^^^^
                2.45, 2.35, 2.25, 2.15, 2.05, 1.95,
  //            ^^^^  ^^^^  ^^^^  ^^^^  ^^^^  ^^^^
                1.85, 1.75, 1.66, 1.55,
  //            ^^^^
                      1.96, 2.02, 2.00, 2.00, 2.00, 2.00,
  //                              ^^^^  ^^^^  ^^^^  ^^^^
    2.75, 2.75,
  //^^^^  ^^^^
          2.50, 2.25, 1.95, 1.86, 1.80, 1.80, 1.80,
  //      ^^^^  ^^^^  ^^^^        ^^^^  ^^^^  ^^^^
          1.80, 1.80, 1.80, 1.80, 1.80, 1.80, 1.80,
  //      ^^^^  ^^^^  ^^^^  ^^^^  ^^^^  ^^^^  ^^^^
                1.80, 1.80, 1.80, 1.80, 1.80, 1.80,
  //            ^^^^  ^^^^  ^^^^  ^^^^  ^^^^  ^^^^
                1.80, 1.80, 1.80, 1.80,
  //            ^^^^  ^^^^  ^^^^  ^^^^
                            1.80,       1.80,       1.80,
  //                        ^^^^        ^^^^        ^^^^
    1.30, 1.50
  //^^^^  ^^^^
  };

  realtype const IonicRadius[nElementNames] = {
       0.79, 0.49, 2.05, 1.40, 1.17, 0.91, 0.75, 0.65, 0.57, 0.51,
       2.23, 1.72, 1.82, 1.46, 1.23, 1.09, 0.97, 0.88, 2.77, 2.23,
       2.09, 2.00, 1.92, 1.85, 1.79, 1.72, 1.67, 1.62, 1.57, 1.53,
       1.81, 1.52, 1.33, 1.22, 1.12, 1.03, 2.98, 2.45, 2.27, 2.16,
       2.09, 2.01, 1.95, 1.89, 1.83, 1.79, 1.75, 1.71, 2.00, 1.72,
       1.53, 1.42, 1.32, 1.24, 3.34, 2.78, 2.74, 2.16, 2.09, 2.02,
       1.97, 1.92, 1.87, 1.83, 1.79, 1.76, 2.08, 1.81, 1.63, 1.53,
       1.43, 1.34, 3.50, 3.00, 3.20, 2.70, 2.67, 2.64, 2.62, 2.59,
       2.56, 2.54, 2.51, 2.49, 2.47, 2.45, 2.42, 2.40, 2.25, 3.16,
       3.14, 3.11, 3.08, 3.05, 3.02, 2.99, 2.97, 2.95, 2.92, 2.90,
       2.87, 2.85
  };

  cpstr const ElementMetal[nElementMetals] = {
    "LI", "BE", "NA", "MG", "AL", " K", "CA", "SC", "TI", " V",
    "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "RB", "SR", " Y",
    "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN",
    "SN", "SB", "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM",
    "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF",
    "TA", " W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB",
    "BI", "PO", "FR", "RA", "AC", "TH", "PA", " U", "NP", "PU",
    "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF",
    "DB", "SG", "BH", "HS", "MT", "UN", "UU", "UB", "UQ", "UH",
    "UO"
  };


  cpstr const HydAtomName[nHydAtomNames] = {
    "0H", "1H", "2H", "3H", "4H", "5H", "6H", "7H", "8H", "9H",
    "HH", "*H", "'H", """H"
  };



  bool isMetal ( cpstr element )  {
  char name[3];
  bool isThere;
  int  i;
    if (!element[1])  {
      name[0] = ' ';
      name[1] = element[0];
    } else
      strncpy ( name,element,2 );
    name[2] = char(0);
    isThere = false;
    for (i=0;(i<nElementMetals) && (!isThere);i++)
      isThere = (!strcmp(ElementMetal[i],name));
    return isThere;
  }

  int  getElementNo ( cpstr element )  {
  int  type=0;
  char El[3];
    if ((!element[1]) || (element[1]==' '))  {
      El[0] = ' ';
      El[1] = element[0];
    } else  {
      El[0] = element[0];
      El[1] = element[1];
    }
    El[2] = char(0);
    UpperCase ( El );
    while (type<nElementNames)
      if (!strcmp(El,ElementName[type]))  break;
                                    else  type++;
    if (type>=nElementNames)  return ELEMENT_UNKNOWN;
    return type+1;  // so that hydrogen is 1
  }

  realtype  getMolecWeight ( cpstr element )  {
  int  type=0;
  char El[3];
    if ((!element[1]) || (element[1]==' '))  {
      El[0] = ' ';
      El[1] = element[0];
    } else  {
      El[0] = element[0];
      El[1] = element[1];
    }
    El[2] = char(0);
    UpperCase ( El );
    while (type<nElementNames)
      if (!strcmp(El,ElementName[type]))  break;
                                    else  type++;
    if (type>=nElementNames)  return 1.0;
    return MolecWeight[type];
  }

  realtype  getCovalentRadius ( cpstr element )  {
  int  type=0;
  char El[3];
    if ((!element[1]) || (element[1]==' '))  {
      El[0] = ' ';
      El[1] = element[0];
    } else  {
      El[0] = element[0];
      El[1] = element[1];
    }
    El[2] = char(0);
    UpperCase ( El );
    while (type<nElementNames)
      if (!strcmp(El,ElementName[type]))  break;
                                    else  type++;
    if (type>=nElementNames)  return 2.2*CovalentRadius[0];
    return CovalentRadius[type];
  }

  realtype  getVdWaalsRadius ( cpstr element )  {
  int  type=0;
  char El[3];
    if ((!element[1]) || (element[1]==' '))  {
      El[0] = ' ';
      El[1] = element[0];
    } else  {
      El[0] = element[0];
      El[1] = element[1];
    }
    El[2] = char(0);
    UpperCase ( El );
    while (type<nElementNames)
      if (!strcmp(El,ElementName[type]))  break;
                                    else  type++;
    if (type>=nElementNames)  return 1.8;
    return VdWaalsRadius[type];
  }

  cpstr const ResidueName[nResNames] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYH", "GLN",
    "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HEM", "WAT", "SUL", "END", "DUM"
  };

  int getResidueNo ( cpstr resName )  {
  int i,m;
    m = -1;
    for (i=0;(i<nResNames) && (m<0);i++)
      if (!strcmp(resName,ResidueName[i]))
        m = i;
    return m;
  }

  char const ResidueName1[nResNames] = {
    'A', 'R', 'N', 'D', 'C', 'C', 'Q',
    'E', 'G', 'H', 'I', 'L', 'K', 'M',
    'F', 'P', 'S', 'T', 'W', 'Y', 'V',
    'X', 'O', 'U', 'Z', 'Z'
  };


  cpstr const StdSolventName[nSolventNames] = {
    "ADE", "CYT", "GUA", "INO", "THY", "URA",
    "WAT", "HOH", "TIP", "H2O", "DOD", "MOH"
  };


  AAProperty const AAProperties[nAminoacidNames] = {
    { "ALA",  1.8,  0.0,  0.42 },
    { "ARG", -4.5,  1.0, -1.37 },
    { "ASN", -3.5,  0.0, -0.82 },
    { "ASP", -3.5, -1.0, -1.05 },
    { "ASX", -3.5, -0.5, -1.05 },  // by analogy
    { "CYS",  2.5,  0.0,  1.34 },
    { "CYH",  2.5,  0.0,  1.34 },  // by analogy
    { "GLN", -3.5,  0.0, -0.30 },
    { "GLU", -3.5, -1.0, -0.87 },
    { "GLX", -3.5, -0.5,  0.0  },  // by analogy
    { "GLY", -0.4,  0.0,  0.0  },
    { "HIS", -3.2,  0.0,  0.18 },
    { "ILE",  4.5,  0.0,  2.46 },
    { "LEU",  3.8,  0.0,  2.32 },
    { "LYS", -3.9,  1.0, -1.35 },
    { "MET",  1.9,  0.0,  1.68 },
    { "PHE",  2.8,  0.0,  2.44 },
    { "PRO", -1.6,  0.0,  0.98 },
    { "SER", -0.8,  0.0, -0.05 },
    { "THR", -0.7,  0.0,  0.35 },
    { "TRP", -0.9,  0.0,  3.07 },
    { "TYR", -1.3,  0.0,  1.31 },
    { "VAL",  4.2,  0.0,  1.66 }
  };

  realtype GetAAHydropathy ( cpstr resName )  {
  int  i1,j;
    i1 = -1;
    j  = 0;
    while (j<nAminoacidNames)
      if (!strcasecmp(resName,AAProperties[j].name))  {
        i1 = j;
        break;
      } else
        j++;
    if (i1<0)  return -MaxReal;
    return AAProperties[i1].hydropathy;
  }

  realtype GetAACharge ( cpstr resName )  {
  int  i1,j;
    i1 = -1;
    j  = 0;
    while (j<nAminoacidNames)
      if (!strcasecmp(resName,AAProperties[j].name))  {
        i1 = j;
        break;
      } else
        j++;
    if (i1<0)  return 0.0;
    return AAProperties[i1].charge;
  }

  realtype GetAASolvationEnergy ( cpstr resName )  {
  int  i1,j;
    i1 = -1;
    j  = 0;
    while (j<nAminoacidNames)
      if (!strcasecmp(resName,AAProperties[j].name))  {
        i1 = j;
        break;
      } else
        j++;
    if (i1<0)  return 0.0;
    return AAProperties[i1].relSolvEnergy;
  }


  int const AASimilarity[nAminoacidNames][nAminoacidNames] = {
  /*  A A A A A C C G G G G H I L L M P P S T T T V            */
  /*  L R S S S Y Y L L L L I L E Y E H R E H R Y A            */
  /*  A G N P X S H N U X Y S E U S T E O R R P R L            */
    { 5,0,0,0,0,0,0,0,0,2,2,0,2,2,0,2,2,2,1,1,2,2,2 },  /* ALA */
    { 0,5,2,2,2,0,0,2,2,0,0,2,0,0,4,0,0,0,0,0,0,0,0 },  /* ARG */
    { 0,2,5,3,5,0,0,3,3,0,0,0,0,0,2,0,0,0,0,0,0,0,0 },  /* ASN */
    { 0,2,3,5,5,0,0,3,3,0,0,0,0,0,2,0,0,0,0,0,0,0,0 },  /* ASP */
    { 0,2,5,5,5,0,0,3,3,0,0,0,0,0,2,0,0,0,0,0,0,0,0 },  /* ASX */
    { 0,0,0,0,0,5,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },  /* CYS */
    { 0,0,0,0,0,4,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },  /* CYH */
    { 0,2,3,3,3,0,0,5,3,0,0,0,0,0,2,0,0,0,0,0,0,0,0 },  /* GLN */
    { 0,2,3,3,3,0,0,3,5,0,0,0,0,0,2,0,0,0,0,0,0,0,0 },  /* GLU */
    { 2,0,0,0,0,0,0,0,0,5,5,0,1,1,0,1,1,2,0,0,1,1,1 },  /* GLX */
    { 2,0,0,0,0,0,0,0,0,5,5,0,1,1,0,1,1,2,0,0,1,1,1 },  /* GLY */
    { 0,2,0,0,0,0,0,0,0,0,0,5,1,1,2,2,2,0,1,1,2,2,1 },  /* HIS */
    { 2,0,0,0,0,0,0,0,0,1,1,1,5,4,0,4,4,0,0,0,4,4,4 },  /* ILE */
    { 2,0,0,0,0,0,0,0,0,1,1,1,4,5,0,4,4,0,0,0,4,4,4 },  /* LEU */
    { 0,4,2,2,2,0,0,2,2,0,0,2,0,0,5,0,0,0,0,0,0,0,0 },  /* LYS */
    { 2,0,0,0,0,0,0,0,0,1,1,2,4,4,0,5,4,0,0,0,4,4,4 },  /* MET */
    { 2,0,0,0,0,0,0,0,0,1,1,2,4,4,0,4,5,0,0,0,5,5,4 },  /* PHE */
    { 2,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0,0,5,0,0,0,0,0 },  /* PRO */
    { 1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,5,3,0,0,0 },  /* SER */
    { 1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,3,5,0,0,0 },  /* THR */
    { 2,0,0,0,0,0,0,0,0,1,1,2,4,4,0,4,5,0,0,0,5,5,4 },  /* TRP */
    { 2,0,0,0,0,0,0,0,0,1,1,2,4,4,0,4,5,0,0,0,5,5,4 },  /* TYR */
    { 2,0,0,0,0,0,0,0,0,1,1,1,4,4,0,4,4,0,0,0,4,4,5 }   /* VAL */
  };

  int  GetAAPIndex ( cpstr resName )  {
  // 0..nAminoacidNames-1
  int i,k;
    k = -1;
    for (i=0;(i<nAminoacidNames) && (k<0);i++)
      if (!strcasecmp(resName,AAProperties[i].name))
        k = i;
    return k;
  }


  int  GetAASimilarity ( cpstr resName1, cpstr resName2 )  {
  int  i1,i2,j;

    i1 = -1;
    j  = 0;
    while (j<nAminoacidNames)
      if (!strcasecmp(resName1,AAProperties[j].name))  {
        i1 = j;
        break;
      } else
        j++;
    if (i1<0)  return -1;

    i2 = -1;
    j  = 0;
    while (j<nAminoacidNames)
      if (!strcasecmp(resName2,AAProperties[j].name))  {
        i2 = j;
        break;
      } else
        j++;
    if (i2<0)  return -2;

    return AASimilarity[i1][i2];

  }





  /*


  pstr const AminoacidName[nAminoacidNames] = {
    "ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "CYH", "GLN",
    "GLU", "GLX", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
  };


  //    The higher is the hydropathy scale value, the more
  //  hydrophobic the residue is.
  //    Some sources suggest that sure hydrophilic residues
  //  are those with hydropathy scale value less than -1.5,
  //  while sure hydropoobic residues have hydropathy scale
  //  value greater than 0.0.
  //    The scale is after J. Kyte and R. F. Doolittle,
  //  A simple method for displaying the hydropathic character
  //  of a protein, J. Mol. Biol. (1982) 157, 105-132
  realtype const AAHydropathyScale[nAminoacidNames] = {
      1.8,  -4.5,  -3.5,  -3.5,  -3.5,   2.5,   2.5,  -3.5,
     -3.5,  -3.5,  -0.4,  -3.2,   4.5,   3.8,  -3.9,   1.9,
      2.8,  -1.6,  -0.8,  -0.7,  -0.9,  -1.3,   4.2
  };

  realtype GetAAHydropathy ( cpstr resName )  {
  int  i1,j;
    i1 = -1;
    j  = 0;
    while (j<nAminoacidNames)
      if (!strcasecmp(resName,AminoacidName[j]))  {
        i1 = j;
        break;
      } else
        j++;
    if (i1<0)  return -MaxReal;
    return AAHydropathyScale[i1];
  }

  //  Acid residues are those with negative charge
  //  Base residues are those with positive charge
  realtype const AACharge[nAminoacidNames] = {
      0.0,   1.0,   0.0,  -1.0,  -0.5,   0.0,   0.0,   0.0,
     -1.0,  -0.5,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0,
      0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0
  };

  int  GetAASimilarity ( cpstr resName1, cpstr resName2 )  {
  int  i1,i2,j;

    i1 = -1;
    j  = 0;
    while (j<nAminoacidNames)
      if (!strcasecmp(resName1,AminoacidName[j]))  {
        i1 = j;
        break;
      } else
        j++;
    if (i1<0)  return -1;

    i2 = -1;
    j  = 0;
    while (j<nAminoacidNames)
      if (!strcasecmp(resName2,AminoacidName[j]))  {
        i2 = j;
        break;
      } else
        j++;
    if (i2<0)  return -2;

    return AASimilarity[i1][i2];

  }

  bool isAminoacid ( cpstr resName )  {
  bool isThere;
  int     i;
    isThere = false;
    for (i=0;(i<nAminoacidNames) && (!isThere);i++)
      isThere = (!strcmp(AminoacidName[i],resName));
    return isThere;
  }

  */


  cpstr const NucleotideName[nNucleotideNames] = {
     "A",  "C",  "G",  "I",   "T",   "U",
    "+A", "+C", "+G", "+I",  "+T",  "+U",
    "DA", "DC", "DG", "DI",  "DT",  "DU",
    "RA", "RC", "RG", "RU", "5NC", "TYD"
  };


  bool isSolvent ( cpstr resName )  {
  bool isThere;
  int     i;
    isThere = false;
    for (i=0;(i<nSolventNames) && (!isThere);i++)
      isThere = (!strcmp(StdSolventName[i],resName));
    return isThere;
  }

  bool isAminoacid ( cpstr resName )  {
  bool isThere;
  int     i;
    isThere = false;
    for (i=0;(i<nAminoacidNames) && (!isThere);i++)
      isThere = (!strcmp(AAProperties[i].name,resName));
    return isThere;
  }

  bool isNucleotide ( cpstr resName )  {
  bool isThere;
  int     i;
    isThere = false;
    for (i=0;(i<nNucleotideNames) && (!isThere);i++)
      isThere = (!strcmp(NucleotideName[i],resName));
    return isThere;
  }

  int isDNARNA ( cpstr resName )  {
  bool isThere;
  int     i;
    isThere = false;
    for (i=0;(i<nNucleotideNames) && (!isThere);i++)
      isThere = (!strcmp(NucleotideName[i],resName));
    if (!isThere)         return 0;  // neither
    if (resName[0]=='D')  return 1;  // DNA
    return 2;                        // RNA
  }

  bool isSugar ( cpstr resName )  {
  UNUSED_ARGUMENT(resName);
    return false;
  }


  cpstr const Res1Code[] = {

    // standard aminoacids

    "ALA A",  // Alanine
    "ARG R",  // Arginine
    "ASN N",  // Asparagine
    "ASP D",  // Aspartic acid (Aspartate)
    "CYS C",  // Cysteine
    "GLN Q",  // Glutamine
    "GLU E",  // Glutamic acid (Glutamate)
    "GLY G",  // Glycine
    "HIS H",  // Histidine
    "ILE I",  // Isoleucine
    "LEU L",  // Leucine
    "LYS K",  // Lysine
    "MET M",  // Methionine
    "PHE F",  // Phenylalanine
    "PRO P",  // Proline
    "SER S",  // Serine
    "THR T",  // Threonine
    "TRP W",  // Tryptophan
    "TYR Y",  // Tyrosine
    "VAL V",  // Valine
    "ASX B",  // Aspartic acid or Asparagine
    "GLX Z",  // Glutamine or Glutamic acid.
    //  ???     X       Any amino acid.

    // other

    "1PA A",   "1PI A",   "2AS D",   "2ML L",   "2MR R",   "3GA A",
    "5HP E",   "ACB D",   "ACL R",   "AGM R",   "AHB D",   "ALM A",
    "ALN A",   "ALO T",   "ALT A",   "ALY K",   "APH A",   "APM A",
    "AR2 R",   "ARM R",   "ARO R",   "ASA D",   "ASB D",   "ASI D",
    "ASK D",   "ASL D",   "ASQ D",   "AYA A",   "B1F A",   "B2A A",
    "B2F A",   "B2I I",   "B2V V",   "BAL A",   "BCS C",   "BFD D",
    "BHD D",   "BLE L",   "BLY K",   "BNN F",   "BNO L",   "BTA L",
    "BTC C",   "BTR W",   "BUC C",   "BUG L",   "C5C C",   "C6C C",
    "CAF C",   "CAS C",   "CAY C",   "CCS C",   "CEA C",   "CGU E",
    "CHG G",   "CHP G",   "CLB A",   "CLD A",   "CLE L",   "CME C",
    "CMT C",   "CSB C",   "CSD A",   "CSE C",   "CSO C",   "CSP C",
    "CSR C",   "CSS C",   "CSW C",   "CSX C",   "CSY C",   "CSZ C",
    "CTH T",   "CXM M",   "CY1 C",   "CYM C",   "CZZ C",   "DAH F",
    "DAL A",   "DAM A",   "DAR R",   "DAS D",   "DBY Y",   "DCY C",
    "DGL E",   "DGN Q",   "DHI H",   "DHN V",   "DIL I",   "DIV V",
    "DLE L",   "DLY K",   "DNP A",   "DOH D",   "DPH F",   "DPN F",
    "DPR P",   "DSE S",   "DSN S",   "DSP D",   "DTH T",   "DTR W",
    "DTY Y",   "DVA V",   "EFC C",   "EHP F",   "EYS C",   "FLA A",
    "FLE L",   "FME M",   "FTY Y",   "GGL E",   "GHP G",   "GSC G",
    "GT9 C",   "H5M P",   "HAC A",   "HAR R",   "HIC H",   "HIP H",
    "HMR R",   "HPH F",   "HPQ F",   "HTR W",   "HV5 A",   "HYP P",
    "IAS N",   "IIL I",   "ILG Q",   "IML I",   "IN2 K",   "ISO A",
    "IVA V",   "IYR Y",   "KCX K",   "KPH K",   "LLY K",   "LOL L",
    "LPL L",   "LTR W",   "LYM K",   "LYZ K",   "M3L K",   "MAA A",
    "MAI R",   "MEN N",   "MGN Q",   "MGY G",   "MHL L",   "MHO M",
    "MHS H",   "MIS S",   "MLE L",   "MLY K",   "MLZ K",   "MME M",
    "MNL L",   "MNV V",   "MPQ G",   "MSE M",   "MSO M",   "MTY Y",
    "MVA V",   "NAL A",   "NAM A",   "NCY C",   "NEM H",   "NEP H",
    "NFA F",   "NIT A",   "NLE L",   "NLN L",   "NNH R",   "NPH C",
    "NVA V",   "OAS S",   "OCS C",   "OCY C",   "OMT M",   "OPR R",
    "PAQ F",   "PBB C",   "PCA E",   "PEC C",   "PGY G",   "PHA F",
    "PHD N",   "PHI F",   "PHL F",   "PHM F",   "PLE L",   "POM P",
    "PPH F",   "PPN F",   "PR3 C",   "PRR A",   "PRS P",   "PTH Y",
    "PTR Y",   "PYA A",   "RON V",   "S1H S",   "SAC S",   "SAH C",
    "SAM M",   "SBD A",   "SBL A",   "SCH C",   "SCS C",   "SCY C",
    "SEB S",   "SEG A",   "SEP S",   "SET S",   "SHC C",   "SHP G",
    "SLZ K",   "SMC C",   "SME M",   "SNC C",   "SOC C",   "STY Y",
    "SVA S",   "TBG G",   "TCR W",   "THC T",   "THO T",   "TIH A",
    "TNB C",   "TPL W",   "TPO T",   "TPQ F",   "TRF W",   "TRG K",
    "TRN W",   "TRO W",   "TYB Y",   "TYI Y",   "TYN Y",   "TYQ Y",
    "TYS Y",   "TYY A",   "VAD V",   "VAF V",   "YOF Y",   ""

  };


  void Get1LetterCode ( cpstr res3name, pstr res1code )  {
  char r[4];
  int  i;

    strncpy ( r,res3name,3 );
    r[3] = char(0);
    UpperCase ( r );
    i = 0;
    res1code[0] = char(1);
    while (Res1Code[i][0])  {
      if (Res1Code[i][0]==r[0])  {
        if (Res1Code[i][1]==r[1])  {
          if (Res1Code[i][2]==r[2])  {
            res1code[0] = Res1Code[i][4];
            break;
          }
        }
      }
      i++;
    }

    if (res1code[0]!=char(1))  res1code[1] = char(0);
    else if (isNucleotide(r))  strcpy ( res1code,r   );
                         else  strcpy ( res1code,"X" );

  }

  void Get1LetterCode ( cpstr res3name, char & res1code )  {
  char r[4];
  int  i;

    strncpy ( r,res3name,3 );
    r[3] = char(0);
    UpperCase ( r );
    i = 0;
    res1code = char(1);
    while (Res1Code[i][0])  {
      if (Res1Code[i][0]==r[0])  {
        if (Res1Code[i][1]==r[1])  {
          if (Res1Code[i][2]==r[2])  {
            res1code = Res1Code[i][4];
            break;
          }
        }
      }
      i++;
    }

    if (res1code==char(1))  {
      if (isNucleotide(r))  res1code = r[0];
                      else  res1code = 'X';
    }

  }


  void Get3LetterCode ( cpstr res1name, pstr res3code )  {
  int  i;

    strcpy ( res3code,"XXX" );
    i = 0;
    while (Res1Code[i][0])  {
      if (Res1Code[i][4]==res1name[0])  {
        res3code[0] = Res1Code[i][0];
        res3code[1] = Res1Code[i][1];
        res3code[2] = Res1Code[i][2];
        break;
      }
      i++;
    }

  }

} // namespace mmdb
