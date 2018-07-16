//  $Id: mmdb_ficif.cpp,v 1.19 2012/01/26 17:52:20 ekr Exp $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2008.
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
//    24.04.03   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  MMDB_FICIF  <implementation>
//       ~~~~~~~~~
//  **** Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2000-2008
//
//  =================================================================
//

#include "mmdb_mmcif_.h"
#include "mmdb_ficif.h"


//  ==================================================================

mmcif::PData mmCIFData = NULL;


FORTRAN_SUBR ( MMDB_FCIF_INIT, mmdb_fcif_init,(),(),() )  {
  InitMatType();
  mmCIFData = NULL;
}

void MMDB_CCIF_Init()  {
  InitMatType();
  mmCIFData = NULL;
}

FORTRAN_SUBR ( MMDB_FCIF_QUIT, mmdb_fcif_quit,(),(),() )  {
  if (mmCIFData)  delete mmCIFData;
  mmCIFData = NULL;
}

void MMDB_CCIF_Quit()  {
  if (mmCIFData)  delete mmCIFData;
  mmCIFData = NULL;
}

pstr makeString ( pstr S, int SLen, pstr FS, int FSLen )  {
  GetStrTer ( S,FS,FSLen,SLen,FSLen );
  CutSpaces ( S,SCUTKEY_END );
  return S;
}

FORTRAN_SUBR ( MMDB_FCIF_CREATE, mmdb_fcif_create,
               (    // lengths-at-end list
                fpstr DataName,     // file name
                int   DataName_len  // fortran-hidden length of DataName
               ), ( // lengths-in-structure list
                fpstr DataName
               ), ( // lengths-follow list
                fpstr DataName, int DataName_len
               ) )  {
char S[500];

  if (mmCIFData)  delete mmCIFData;
  mmCIFData = new mmcif::Data ( makeString(S,sizeof(S),
                               FTN_STR(DataName),FTN_LEN(DataName)) );

}

void MMDB_CCIF_Create ( pstr DataName )  {
  if (mmCIFData)  delete mmCIFData;
  mmCIFData = new mmcif::Data ( DataName );
}


FORTRAN_SUBR ( MMDB_FCIF_WRITE, mmdb_fcif_write,
               (    // lengths-at-end list
                fpstr FileName,     // file name
                int * iRet,         // return code
                int   FileName_len  // fortran-hidden length of FileName
               ), ( // lengths-in-structure list
    fpstr FileName, int *iRet
               ), ( // lengths-follow list
    fpstr FileName, int FileName_len, int * iRet
               ) )  {
pstr S;

  if (!mmCIFData)
    *iRet = -1000;
  else  {
    S = new char[FTN_LEN(FileName)+10];
    if (mmCIFData->WriteMMCIFData(makeString(S,FTN_LEN(FileName)+5,
                                  FTN_STR(FileName),FTN_LEN(FileName))))
         *iRet = 0;
    else *iRet = 1;
    delete[] S;
  }

}


int MMDB_CCIF_Write ( pstr FileName )  {
  if (!mmCIFData)  return -1000;
  else if (mmCIFData->WriteMMCIFData(FileName))  return 0;
                                           else  return 1;
}


FORTRAN_SUBR ( MMDB_FCIF_PUTDATE, mmdb_fcif_putdate,
               (    // lengths-at-end list
                fpstr CatName,     // category name
                fpstr Tag,         // tag
                int * iRet,        // return code
                int   CatName_len, // fortran-hidden length of CatName
                int   Tag_len      // fortran-hidden length of Tag
               ), ( // lengths-in-structure list
                fpstr CatName, fpstr Tag, int * iRet
               ), ( // lengths-follow list
                fpstr CatName, int CatName_len,
                fpstr Tag,     int Tag_len,
                int * iRet
               ) )  {
char CN[200],TN[200];

  if (!mmCIFData)
       *iRet = -1000;
  else *iRet = mmCIFData->PutDate ( makeString(CN,sizeof(CN),
                                    FTN_STR(CatName),FTN_LEN(CatName)),
                                    makeString(TN,sizeof(TN),
                                    FTN_STR(Tag),FTN_LEN(Tag)) );
}

int MMDB_CCIF_PutDate ( pstr CatName, pstr Tag )  {
  if (!mmCIFData)  return -1000;
             else  return  mmCIFData->PutDate ( CatName,Tag );
}



FORTRAN_SUBR ( MMDB_FCIF_PUTDOT, mmdb_fcif_putdot,
               (    // lengths-at-end list
                fpstr CatName,     // category name
                fpstr Tag,         // tag
                int * iRet,        // return code
                int   CatName_len, // fortran-hidden length of CatName
                int   Tag_len      // fortran-hidden length of Tag
               ), ( // lengths-in-structure list
                fpstr CatName, fpstr Tag, int * iRet
               ), ( // lengths-follow list
                fpstr CatName, int CatName_len,
                fpstr Tag,     int Tag_len,
                int * iRet
               ) )  {
char CN[200],TN[200];
  if (!mmCIFData)
       *iRet = -1000;
  else *iRet = mmCIFData->PutNoData ( mmcif::CIF_NODATA_DOT,
                                      makeString(CN,sizeof(CN),
                                      FTN_STR(CatName),FTN_LEN(CatName)),
                                      makeString(TN,sizeof(TN),
                                      FTN_STR(Tag),FTN_LEN(Tag)) );
}


int MMDB_CCIF_PutDot ( pstr CatName, pstr Tag )  {
  if (!mmCIFData)  return  -1000;
             else  return  mmCIFData->PutNoData ( mmcif::CIF_NODATA_DOT,
                                                  CatName,Tag );
}


FORTRAN_SUBR ( MMDB_FCIF_PUTQUESTION, mmdb_fcif_putquestion,
               (    // lengths-at-end list
                fpstr CatName,     // category name
                fpstr Tag,         // tag
                int * iRet,        // return code
                int   CatName_len, // fortran-hidden length of CatName
                int   Tag_len      // fortran-hidden length of Tag
               ), ( // lengths-in-structure list
                fpstr CatName, fpstr Tag, int * iRet
               ), ( // lengths-follow list
                fpstr CatName, int CatName_len,
                fpstr Tag,     int Tag_len,
                int * iRet
               ) )  {
char CN[200],TN[200];
  if (!mmCIFData)
       *iRet = -1000;
  else *iRet = mmCIFData->PutNoData ( mmcif::CIF_NODATA_QUESTION,
                                      makeString(CN,sizeof(CN),
                                      FTN_STR(CatName),FTN_LEN(CatName)),
                                      makeString(TN,sizeof(TN),
                                      FTN_STR(Tag),FTN_LEN(Tag)) );
}


int MMDB_CCIF_PutQuestion ( pstr CatName, pstr Tag )  {
  if (!mmCIFData)
    return  -1000;
  else
    return  mmCIFData->PutNoData ( mmcif::CIF_NODATA_QUESTION,
                                                  CatName,Tag );
}


FORTRAN_SUBR ( MMDB_FCIF_PUTSTRING, mmdb_fcif_putstring,
               (    // lengths-at-end list
                fpstr Data,        // data string to store
                fpstr CatName,     // category name
                fpstr Tag,         // tag
                int * iRet,        // return code
                int   Data_len,    // fortran-hidden length of Data
                int   CatName_len, // fortran-hidden length of CatName
                int   Tag_len      // fortran-hidden length of Tag
               ), ( // lengths-in-structure list
                fpstr Data, fpstr CatName, fpstr Tag, int * iRet
               ), ( // lengths-follow list
    fpstr Data,    int Data_len,
                fpstr CatName, int CatName_len,
                fpstr Tag,     int Tag_len,
                int * iRet
               ) )  {
char CN[200],TN[200];
pstr S;

  if (!mmCIFData)
    *iRet = -1000;
  else  {
    S = new char[FTN_LEN(Data)+10];
    *iRet = mmCIFData->PutString ( makeString(S,FTN_LEN(Data)+5,
           FTN_STR(Data),FTN_LEN(Data)),
                                   makeString(CN,sizeof(CN),
                                   FTN_STR(CatName),FTN_LEN(CatName)),
                                   makeString(TN,sizeof(TN),
                                   FTN_STR(Tag),FTN_LEN(Tag)) );
    delete[] S;
  }

}

int MMDB_CCIF_PutString ( pstr Data, pstr CatName, pstr Tag )  {
  if (!mmCIFData)  return -1000;
             else  return mmCIFData->PutString ( Data,CatName,Tag );
}


FORTRAN_SUBR ( MMDB_FCIF_PUTREAL, mmdb_fcif_putreal,
               (    // lengths-at-end list
                apireal * V,       // real value to store
                fpstr CatName,     // category name
                fpstr Tag,         // tag
                int * iRet,        // return code
                int   CatName_len, // fortran-hidden length of CatName
                int   Tag_len      // fortran-hidden length of Tag
               ), ( // lengths-in-structure list
                apireal * V, fpstr CatName, fpstr Tag, int * iRet
               ), ( // lengths-follow list
    apireal * V,
                fpstr CatName, int CatName_len,
                fpstr Tag,     int Tag_len,
                int * iRet
               ) )  {
char CN[200],TN[200];

  if (!mmCIFData) *iRet = -1000;
             else *iRet = mmCIFData->PutReal ( *V,
                                   makeString(CN,sizeof(CN),
                                   FTN_STR(CatName),FTN_LEN(CatName)),
                                   makeString(TN,sizeof(TN),
                                   FTN_STR(Tag),FTN_LEN(Tag)) );

}


int MMDB_CCIF_PutReal ( realtype V, pstr CatName, pstr Tag )  {
  if (!mmCIFData)  return  -1000;
             else  return  mmCIFData->PutReal ( V,CatName,Tag );
}


FORTRAN_SUBR ( MMDB_FCIF_PUTINTEGER, mmdb_fcif_putinteger,
               (    // lengths-at-end list
                int * I,           // integer value to store
                fpstr CatName,     // category name
                fpstr Tag,         // tag
                int * iRet,        // return code
                int   CatName_len, // fortran-hidden length of CatName
                int   Tag_len      // fortran-hidden length of Tag
               ), ( // lengths-in-structure list
                int * I, fpstr CatName, fpstr Tag, int * iRet
               ), ( // lengths-follow list
    int * I,
                fpstr CatName, int CatName_len,
                fpstr Tag,     int Tag_len,
                int * iRet
               ) )  {
char CN[200],TN[200];

  if (!mmCIFData) *iRet = -1000;
             else *iRet = mmCIFData->PutInteger ( *I,
                                   makeString(CN,sizeof(CN),
                                   FTN_STR(CatName),FTN_LEN(CatName)),
                                   makeString(TN,sizeof(TN),
                                   FTN_STR(Tag),FTN_LEN(Tag)) );

}

int MMDB_CCIF_PutInteger ( int I, pstr CatName, pstr Tag )  {
  if (!mmCIFData)  return  -1000;
             else  return  mmCIFData->PutInteger ( I,CatName,Tag );
}


FORTRAN_SUBR ( MMDB_FCIF_PUTLOOPDOT, mmdb_fcif_putloopdot,
               (    // lengths-at-end list
                fpstr CatName,     // category name
                fpstr Tag,         // tag
                int * nrow,        // row number
                int * iRet,        // return code
                int   CatName_len, // fortran-hidden length of CatName
                int   Tag_len      // fortran-hidden length of Tag
               ), ( // lengths-in-structure list
                fpstr CatName, fpstr Tag, int * nrow, int * iRet
               ), ( // lengths-follow list
                fpstr CatName, int CatName_len,
                fpstr Tag,     int Tag_len,
                int * nrow,    int * iRet
               ) )  {
char CN[200],TN[200];
  if (!mmCIFData)
       *iRet = -1000;
  else *iRet = mmCIFData->PutLoopNoData ( mmcif::CIF_NODATA_DOT,
                                      makeString(CN,sizeof(CN),
                                      FTN_STR(CatName),FTN_LEN(CatName)),
                                      makeString(TN,sizeof(TN),
                                      FTN_STR(Tag),FTN_LEN(Tag)),*nrow );
}


int MMDB_CCIF_PutLoopDot ( pstr CatName, pstr Tag, int nrow )  {
  if (!mmCIFData)
    return  -1000;
  else
    return  mmCIFData->PutLoopNoData ( mmcif::CIF_NODATA_DOT,
                                                  CatName,Tag,nrow );
}


FORTRAN_SUBR ( MMDB_FCIF_PUTLOOPQUESTION, mmdb_fcif_putloopquestion,
               (    // lengths-at-end list
                fpstr CatName,     // category name
                fpstr Tag,         // tag
                int * nrow,        // row number
                int * iRet,        // return code
                int   CatName_len, // fortran-hidden length of CatName
                int   Tag_len      // fortran-hidden length of Tag
               ), ( // lengths-in-structure list
                fpstr CatName, fpstr Tag, int * nrow, int * iRet
               ), ( // lengths-follow list
                fpstr CatName, int CatName_len,
                fpstr Tag,     int Tag_len,
                int * nrow,    int * iRet
               ) )  {
char CN[200],TN[200];
  if (!mmCIFData)
       *iRet = -1000;
  else *iRet = mmCIFData->PutLoopNoData ( mmcif::CIF_NODATA_QUESTION,
                                      makeString(CN,sizeof(CN),
                                      FTN_STR(CatName),FTN_LEN(CatName)),
                                      makeString(TN,sizeof(TN),
                                      FTN_STR(Tag),FTN_LEN(Tag)),*nrow );
}

int MMDB_CCIF_PutLoopQuestion ( pstr CatName, pstr Tag, int nrow )  {
  if (!mmCIFData)
    return  -1000;
  else
    return  mmCIFData->PutLoopNoData ( mmcif::CIF_NODATA_QUESTION,
                                                   CatName,Tag,nrow );
}


FORTRAN_SUBR ( MMDB_FCIF_PUTLOOPSTRING, mmdb_fcif_putloopstring,
               (    // lengths-at-end list
                fpstr Data,        // data string to store
                fpstr CatName,     // category name
                fpstr Tag,         // tag
                int * nrow,        // row number
                int * iRet,        // return code
                int   Data_len,    // fortran-hidden length of Data
                int   CatName_len, // fortran-hidden length of CatName
                int   Tag_len      // fortran-hidden length of Tag
               ), ( // lengths-in-structure list
                fpstr Data, fpstr CatName, fpstr Tag,
                int * nrow, int * iRet
               ), ( // lengths-follow list
    fpstr Data,    int Data_len,
                fpstr CatName, int CatName_len,
                fpstr Tag,     int Tag_len,
                int * nrow,    int * iRet
               ) )  {
char CN[200],TN[200];
pstr S;

  if (!mmCIFData)
    *iRet = -1000;
  else  {
    S = new char[FTN_LEN(Data)+10];
    *iRet = mmCIFData->PutLoopString ( makeString(S,FTN_LEN(Data)+5,
               FTN_STR(Data),FTN_LEN(Data)),
                                       makeString(CN,sizeof(CN),
                                       FTN_STR(CatName),FTN_LEN(CatName)),
                                       makeString(TN,sizeof(TN),
                                       FTN_STR(Tag),FTN_LEN(Tag)),*nrow );
    delete[] S;
  }

}


int MMDB_CCIF_PutLoopString ( pstr Data, pstr CatName, pstr Tag, int nrow )  {
  if (!mmCIFData)  return -1000;
             else  return mmCIFData->PutLoopString ( Data,CatName,Tag,nrow );
}


FORTRAN_SUBR ( MMDB_FCIF_PUTLOOPREAL, mmdb_fcif_putloopreal,
               (    // lengths-at-end list
                apireal * V,       // real value to store
                fpstr CatName,     // category name
                fpstr Tag,         // tag
                int * nrow,        // row number
                int * iRet,        // return code
                int   CatName_len, // fortran-hidden length of CatName
                int   Tag_len      // fortran-hidden length of Tag
               ), ( // lengths-in-structure list
                apireal * V, fpstr CatName, fpstr Tag,
                int * nrow, int * iRet
               ), ( // lengths-follow list
    apireal * V,
                fpstr CatName, int CatName_len,
                fpstr Tag,     int Tag_len,
                int * nrow,    int * iRet
               ) )  {
char CN[200],TN[200];

  if (!mmCIFData) *iRet = -1000;
             else *iRet = mmCIFData->PutLoopReal ( *V,
                                   makeString(CN,sizeof(CN),
                                   FTN_STR(CatName),FTN_LEN(CatName)),
                                   makeString(TN,sizeof(TN),
                                   FTN_STR(Tag),FTN_LEN(Tag)),*nrow );

}


int MMDB_CCIF_PutLoopReal ( realtype V, pstr CatName, pstr Tag, int nrow )  {
  if (!mmCIFData)  return -1000;
             else  return mmCIFData->PutLoopReal ( V,CatName,Tag,nrow );
}


FORTRAN_SUBR ( MMDB_FCIF_PUTLOOPINTEGER, mmdb_fcif_putloopinteger,
               (    // lengths-at-end list
                int * I,           // integer value to store
                fpstr CatName,     // category name
                fpstr Tag,         // tag
                int * nrow,        // row number
                int * iRet,        // return code
                int   CatName_len, // fortran-hidden length of CatName
                int   Tag_len      // fortran-hidden length of Tag
               ), ( // lengths-in-structure list
                int * I, fpstr CatName, fpstr Tag,
                int * nrow, int * iRet
               ), ( // lengths-follow list
    int * I,
                fpstr CatName, int CatName_len,
                fpstr Tag,     int Tag_len,
                int * nrow,    int * iRet
               ) )  {
char CN[200],TN[200];

  if (!mmCIFData) *iRet = -1000;
             else *iRet = mmCIFData->PutLoopInteger ( *I,
                                   makeString(CN,sizeof(CN),
                                   FTN_STR(CatName),FTN_LEN(CatName)),
                                   makeString(TN,sizeof(TN),
                                   FTN_STR(Tag),FTN_LEN(Tag)),*nrow );

}


int MMDB_CCIF_PutLoopInteger ( int I, pstr CatName, pstr Tag, int nrow )  {
  if (!mmCIFData)  return -1000;
             else  return mmCIFData->PutLoopInteger ( I,CatName,Tag,nrow );
}
