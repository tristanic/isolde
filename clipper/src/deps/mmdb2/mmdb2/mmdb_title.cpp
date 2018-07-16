//  $Id: mmdb_title.cpp $
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
//    28.09.15   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  MMDB_Title <implementation>
//       ~~~~~~~~~
//  **** Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::TitleContainer  (container of title classes)
//       ~~~~~~~~~  mmdb::ObsLine
//                  mmdb::TitleLine
//                  mmdb::Caveat
//                  mmdb::Compound
//                  mmdb::Source
//                  mmdb::KeyWords
//                  mmdb::ExpData
//                  mmdb::MdlType
//                  mmdb::Author
//                  mmdb::RevData
//                  mmdb::Supersede
//                  mmdb::Journal
//                  mmdb::Remark
//                  mmdb::Biomolecule
//                  mmdb::Title       ( MMDB title section )
//
//   (C) E. Krissinel 2000-2015
//
//  =================================================================
//

#include <string.h>
#include <stdlib.h>

#include "mmdb_title.h"
#include "mmdb_cifdefs.h"

namespace mmdb  {

  //  ==============  TitleContainer  ====================

  PContainerClass TitleContainer::MakeContainerClass ( int ClassID )  {
    switch (ClassID)  {
      default :
      case ClassID_Template  : return
                           ClassContainer::MakeContainerClass(ClassID);
      case ClassID_ObsLine   : return new ObsLine  ();
      case ClassID_CAVEAT    : return new Caveat   ();
      case ClassID_TitleLine : return new TitleLine();
      case ClassID_Compound  : return new Compound ();
      case ClassID_Source    : return new Source   ();
      case ClassID_ExpData   : return new ExpData  ();
      case ClassID_Author    : return new Author   ();
      case ClassID_RevData   : return new RevData  ();
      case ClassID_Supersede : return new Supersede();
      case ClassID_Journal   : return new Journal  ();
      case ClassID_Remark    : return new Remark   ();
    }
  }

  MakeStreamFunctions(TitleContainer)


  //  ================  ObsLine  ===================

  ObsLine::ObsLine() : ContainerClass()  {
    InitObsLine();
  }

  ObsLine::ObsLine ( cpstr S ) : ContainerClass()  {
    InitObsLine();
    ConvertPDBASCII ( S );
  }

  ObsLine::ObsLine ( io::RPStream Object ) : ContainerClass(Object)  {
    InitObsLine();
  }

  ObsLine::~ObsLine() {}

  void  ObsLine::InitObsLine()  {
  int i;
    strcpy ( repDate,"DD-MMM-YYYY" );
    strcpy ( idCode, "----" );
    for (i=0;i<8;i++)
      strcpy ( rIdCode[i],"    " );
  }

  void  ObsLine::PDBASCIIDump ( pstr S, int N )  {
  //  makes the ASCII PDB OBSLTE line number N
  //  from the class' data
  int i;
    if (N==0)  strcpy  ( S,"OBSLTE    " );
         else  sprintf ( S,"OBSLTE  %2i",N+1 );
    PadSpaces ( S,80 );
    Date11to9 ( repDate,&(S[11]) );
    strncpy   ( &(S[21]),idCode,4 );
    for (i=0;i<8;i++)
      strncpy ( &(S[31+5*i]),rIdCode[i],4 );
  }

  void  ObsLine::MakeCIF ( mmcif::PData CIF, int )  {
  mmcif::PLoop Loop;
  int          RC,i,j;
  char         DateCIF[20];
    RC = CIF->AddLoop ( CIFCAT_OBSLTE,Loop );
    if (RC!=mmcif::CIFRC_Ok)  {
      // the category was (re)created, provide tags
      Loop->AddLoopTag ( CIFTAG_ID             );
      Loop->AddLoopTag ( CIFTAG_DATE           );
      Loop->AddLoopTag ( CIFTAG_REPLACE_PDB_ID );
      Loop->AddLoopTag ( CIFTAG_PDB_ID         );
    }
    Date11toCIF ( repDate,DateCIF );
    for (i=0;i<8;i++)  {
      j = 0;
      while (rIdCode[i][j] && (rIdCode[i][j]==' '))  j++;
      if (rIdCode[i][j])  {
        Loop->AddString ( pstr("OBSLTE") );
        Loop->AddString ( DateCIF    );
        Loop->AddString ( idCode     );
        Loop->AddString ( rIdCode[i] );
      }
    }
  }

  ERROR_CODE ObsLine::ConvertPDBASCII ( cpstr S )  {
  int i;
    Date9to11 ( &(S[11]),repDate );
    strncpy   ( idCode,&(S[21]),4 );
    idCode[4] = char(0);
    for (i=0;i<8;i++)  {
      strncpy ( rIdCode[i],&(S[31+i*5]),4 );
      rIdCode[i][4] = char(0);
    }
    return Error_NoError;
  }

  ERROR_CODE  ObsLine::GetCIF ( mmcif::PData CIF, int & n )  {
  mmcif::PLoop Loop;
  int         i,RC;
  pstr        F,FDate,FID;
  char        DateCIF [20];
  char        DateCIF0[20];
  IDCode      idCode1;

    Loop = CIF->GetLoop ( CIFCAT_OBSLTE );
    if (!Loop)  {
      n = -1;  // signal to finish processing of this structure
      return Error_EmptyCIF;
    }
    i = 0;
    do  {
      F = Loop->GetString ( CIFTAG_ID,n,RC );
      if (RC)  {
        if (i==0)  n = -1;
        return Error_MissingCIFField;
      }
      if (F)  {
        if (!strcmp(F,"OBSLTE"))  {
          FDate = Loop->GetString ( CIFTAG_DATE,n,RC );
          if ((!RC) && FDate)
                strncpy ( DateCIF,FDate,15 );
          else  strcpy  ( DateCIF,"YYYY-MMM-DD" );
          FID = Loop->GetString ( CIFTAG_REPLACE_PDB_ID,n,RC );
          if ((!RC) && FID)
                strncpy ( idCode1,FID,sizeof(IDCode)-1 );
          else  idCode1[0] = char(0);
          if (i==0)  {
            DateCIFto11 ( DateCIF,repDate );
            DateCIF[11] = char(0);
            strcpy ( idCode  ,idCode1 );
            strcpy ( DateCIF0,DateCIF );
          } else if ((strcmp(DateCIF0,DateCIF)) ||
                     (strcmp(idCode,idCode1)))
            return Error_MissingCIFField;
          FID = Loop->GetString ( CIFTAG_PDB_ID,n,RC );
          if ((!RC) && FID)
               strncpy ( rIdCode[i],FID,sizeof(IDCode)-1 );
          else rIdCode[i][0] = char(0);
          Loop->DeleteField ( CIFTAG_ID            ,n );
          Loop->DeleteField ( CIFTAG_DATE          ,n );
          Loop->DeleteField ( CIFTAG_REPLACE_PDB_ID,n );
          Loop->DeleteField ( CIFTAG_PDB_ID        ,n );
          i++;
        }
      }
      n++;
    } while (i<8);

    return Error_NoError;

  }

  void  ObsLine::Copy ( PContainerClass ObsLine )  {
  int i;
    strcpy ( repDate,PObsLine(ObsLine)->repDate );
    strcpy ( idCode ,PObsLine(ObsLine)->idCode  );
    for (i=0;i<8;i++)
      strcpy ( rIdCode[i],PObsLine(ObsLine)->rIdCode[i] );
  }

  void  ObsLine::write ( io::RFile f )  {
  int  i;
  byte Version=1;
   f.WriteByte    ( &Version      );
   f.WriteTerLine ( repDate,false );
   f.WriteTerLine ( idCode ,false );
   for (i=0;i<8;i++)
     f.WriteTerLine ( rIdCode[i],false );
  }

  void  ObsLine::read  ( io::RFile f ) {
  int  i;
  byte Version;
   f.ReadByte    ( &Version      );
   f.ReadTerLine ( repDate,false );
   f.ReadTerLine ( idCode ,false );
   for (i=0;i<8;i++)
     f.ReadTerLine ( rIdCode[i],false );
  }

  MakeStreamFunctions(ObsLine)


  //  ===================  TitleLine  ======================

  TitleLine::TitleLine() : ContString()  {
    InitTitleLine();
  }

  TitleLine::TitleLine ( cpstr S ) : ContString()  {
    InitTitleLine();
    ConvertPDBASCII ( S );
  }

  TitleLine::TitleLine ( io::RPStream Object ) : ContString(Object)  {
    InitTitleLine();
  }

  TitleLine::~TitleLine() {
  }

  void  TitleLine::InitTitleLine()  {
    CreateCopy ( CIFCategory,CIFCAT_STRUCT );
    CreateCopy ( CIFTag,     CIFTAG_TITLE  );
  }

  ERROR_CODE TitleLine::ConvertPDBASCII ( cpstr S )  {
    if (strlen(S)>10)
         CreateCopy ( Line,&(S[10]) );
    else CreateCopy ( Line,pstr(" ") );
    return Error_NoError;
  }

  void  TitleLine::PDBASCIIDump ( pstr S, int N )  {
    if (N==0)  strcpy  ( S,"TITLE     " );
         else  sprintf ( S,"TITLE   %2i",N+1 );
    strcat ( S,Line );
  }

  void  TitleLine::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte ( &Version );
    ContString::write ( f );
  }

  void  TitleLine::read ( io::RFile f )  {
  byte Version;
    f.ReadByte ( &Version );
    ContString::read ( f );
  }

  MakeStreamFunctions(TitleLine)



  //  ===================  Caveat  ======================

  Caveat::Caveat() : ContString()  {
    InitCaveat();
  }

  Caveat::Caveat ( cpstr S ) : ContString()  {
    InitCaveat();
    ConvertPDBASCII ( S );
  }

  Caveat::Caveat ( io::RPStream Object ) : ContString(Object)  {
    InitCaveat();
  }

  Caveat::~Caveat() {}

  void  Caveat::InitCaveat()  {
    strcpy ( idCode,"----" );
    CreateCopy ( CIFCategory,CIFCAT_DATABASE_PDB_CAVEAT );
    CreateCopy ( CIFTag     ,CIFTAG_TEXT                );
  }

  ERROR_CODE Caveat::ConvertPDBASCII ( cpstr S )  {
    if (strlen(S)>12)  {
      strncpy ( idCode,&(S[11]),4 );
      idCode[4] = char(0);
      if (strlen(S)>19)
            CreateCopy ( Line,&(S[19])  );
      else  CreateCopy ( Line,pstr(" ") );
    } else
      CreateCopy ( Line,pstr(" ") );
    return Error_NoError;
  }

  void  Caveat::PDBASCIIDump ( pstr S, int N )  {
    if (N==0)  strcpy  ( S,"CAVEAT     " );
         else  sprintf ( S,"CAVEAT  %2i ",N+1 );
    strcat ( S,idCode );
    strcat ( S,"    " );
    strcat ( S,Line   );
  }

  void  Caveat::MakeCIF ( mmcif::PData CIF, int N )  {
  char S[500];
    CIF->PutString ( idCode,CIFCAT_DATABASE_PDB_CAVEAT,CIFTAG_ID,false );
    strcpy  ( S,"\n" );
    strncat ( S,Line,sizeof(S)-2 );
    S[sizeof(S)-1] = char(0);
    CIF->PutString ( S,CIFCAT_DATABASE_PDB_CAVEAT,CIFTAG_TEXT,(N!=0) );
  }

/*
  void  Caveat::GetCIF1 ( mmcif::PData CIF, ERROR_CODE & Signal,
                          int & pos )  {
  pstr F;
  int  RC;
    F = CIF->GetString ( CIFCAT_DATABASE_PDB_CAVEAT,CIFTAG_ID,RC );
    if ((!RC) && F)  {
      strncpy ( idCode,F,sizeof(IDCode) );
      idCode[sizeof(IDCode)-1] = char(0);
    }
    ContString::GetCIF1 ( CIF,Signal,pos );
    if (Signal<0)
      CIF->DeleteField ( CIFCAT_DATABASE_PDB_CAVEAT,CIFTAG_ID );
  }
*/

  void  Caveat::Copy ( PContainerClass Caveat )  {
    strcpy ( idCode,PCaveat(Caveat)->idCode );
    ContString::Copy ( Caveat );
  }

  void  Caveat::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte    ( &Version     );
    f.WriteTerLine ( idCode,false );
    ContString::write ( f );
  }

  void  Caveat::read ( io::RFile f )  {
  byte Version;
    f.ReadByte    ( &Version     );
    f.ReadTerLine ( idCode,false );
    ContString::read ( f );
  }

  MakeStreamFunctions(Caveat)



  //  ===================  Compound  ======================

  Compound::Compound() : ContString()  {
    InitCompound();
  }

  Compound::Compound ( cpstr S ) : ContString()  {
    InitCompound();
    ConvertPDBASCII ( S );
  }

  Compound::Compound ( io::RPStream Object ) : ContString(Object)  {
    InitCompound();
  }

  Compound::~Compound() {}

  void  Compound::InitCompound()  {
    CreateCopy ( CIFCategory,CIFCAT_STRUCT         );
    CreateCopy ( CIFTag     ,CIFTAG_NDB_DESCRIPTOR );
  }

  ERROR_CODE Compound::ConvertPDBASCII ( cpstr S )  {
    if (strlen(S)>10)
         CreateCopy ( Line,&(S[10])  );
    else CreateCopy ( Line,pstr(" ") );
    return Error_NoError;
  }

  void  Compound::PDBASCIIDump ( pstr S, int N )  {
    if (N==0)  strcpy  ( S,"COMPND    " );
         else  sprintf ( S,"COMPND  %2i",N+1 );
    strcat ( S,Line );
  }

  void  Compound::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte ( &Version );
    ContString::write ( f );
  }

  void  Compound::read ( io::RFile f )  {
  byte Version;
    f.ReadByte ( &Version );
    ContString::read ( f );
  }

  MakeStreamFunctions(Compound)



  //  ===================  Source  ======================

  Source::Source() : ContString()  {
    InitSource();
  }

  Source::Source ( cpstr S ) : ContString()  {
    InitSource();
    ConvertPDBASCII ( S );
  }

  Source::Source ( io::RPStream Object ) : ContString(Object)  {
    InitSource();
  }

  Source::~Source() {}

  void  Source::InitSource()  {
    CreateCopy ( CIFCategory,CIFCAT_STRUCT );
    CreateCopy ( CIFTag     ,CIFTAG_SOURCE );
  }

  ERROR_CODE Source::ConvertPDBASCII ( cpstr S )  {
    if (strlen(S)>10)
         CreateCopy ( Line,&(S[10])  );
    else CreateCopy ( Line,pstr(" ") );
    return Error_NoError;
  }

  void  Source::PDBASCIIDump ( pstr S, int N )  {
    if (N==0)  strcpy  ( S,"SOURCE    " );
         else  sprintf ( S,"SOURCE  %2i",N+1 );
    strcat ( S,Line );
  }

  void  Source::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte ( &Version );
    ContString::write ( f );
  }

  void  Source::read ( io::RFile f )  {
  byte Version;
    f.ReadByte ( &Version );
    ContString::read ( f );
  }

  MakeStreamFunctions(Source)


  //  ===================  KeyWords  ======================

  KeyWords::KeyWords() : io::Stream()  {
    Init();
  }

  KeyWords::KeyWords ( cpstr S ) : io::Stream()  {
    Init();
    ConvertPDBASCII ( S );
  }

  KeyWords::KeyWords ( io::RPStream Object ) : io::Stream(Object)  {
    Init();
  }

  KeyWords::~KeyWords()  {
    Delete();
  }

  void  KeyWords::Init()  {
    nKeyWords = 0;
    KeyWord   = NULL;
    Cont      = false;
  }

  void  KeyWords::Delete()  {
  int i;
    if (KeyWord)  {
      for (i=0;i<nKeyWords;i++)
        if (KeyWord[i])
          delete[] KeyWord[i];
      delete[] KeyWord;
    }
    nKeyWords = 0;
    KeyWord   = NULL;
    Cont      = false;
  }

  int  KeyWords::ConvertPDBASCII ( cpstr S )  {
  //  we anticipate that length of S is 80 characters
  //  -- pad with spaces if necessary
  char   L[85];
  int    i,k,m;
  pstr * KW;

    i = 10;  // scan PDB line from ith character

    k = nKeyWords;
    if (!Cont)  k++;  // 1st keyword does not continue from previous line
    m = 0;
    while (S[i] && (i<70))  {
      if (S[i]==',')  k++;  // count keywords
      if (S[i]!=' ')  m++;  // count non-spaces to see if the line is empty
      i++;
    }

    if (m==0)  return 0;    //  empty line

    KW = new pstr[k];
    if (KeyWord)  {
      for (i=0;i<nKeyWords;i++)
        KW[i] = KeyWord[i];
      delete[] KeyWord;
    }
    for (i=nKeyWords;i<k;i++)
      KW[i] = NULL;       // null new pointers
    KeyWord = KW;

    i = 10;
    if (Cont)  nKeyWords--;
    while (S[i] && (i<70))  {
      while ((S[i]==' ') && (i<70))  i++;  // skip leading spaces
      if (Cont)  {
        strcpy ( L," " );
        m = 1;
      } else
        m = 0;
      while (S[i] && (S[i]!=',') && (i<70))
        L[m++] = S[i++];
      m--;
      while ((m>0) && (L[m]==' '))  m--;  // remove padding spaces
      m++;
      L[m] = char(0);
      if (Cont)  CreateConcat ( KeyWord[nKeyWords],L );
           else  CreateCopy   ( KeyWord[nKeyWords],L );
      if (S[i]==',')  {
        i++;
        Cont = false;
      } else
        Cont = true;
      nKeyWords++;
    }

    return 0;

  }

  void  KeyWords::PDBASCIIDump ( io::RFile f )  {
  int  N,i,k,m1,m2,ms;
  char S[85];
  char c;
    if (KeyWord)  {
      N = 0;
      i = 0;
      while (i<nKeyWords)  {
        if (N==0)  strcpy  ( S,"KEYWDS    " );
             else  sprintf ( S,"KEYWDS  %2i ",N+1 );
        do  {
          while ((i<nKeyWords) && (!KeyWord[i]))  i++;
          if (i<nKeyWords) {
            m1 = 0;
            while (KeyWord[i][m1])  {
              while (KeyWord[i][m1]==' ')  m1++;
              m2 = m1;
              ms = -1;
              while ((KeyWord[i][m2]) && ((m2-m1)<58))  {
                if (KeyWord[i][m2]==' ')  ms = m2;
                m2++;
              }
              if ((ms<0) || ((m2-m1)<58))  ms = m2;
              c = KeyWord[i][ms];
              KeyWord[i][ms] = char(0);
              strcat ( S,&(KeyWord[i][m1]) );
              KeyWord[i][ms] = c;
              m1 = ms;
              if (c)  {
                PadSpaces   ( S,80 );
                f.WriteLine ( S );
                N++;
                sprintf ( S,"KEYWDS  %2i ",N+1 );
              }
            }
            i++;
            if (i<nKeyWords)  {
              k = strlen(S) + strlen(KeyWord[i]) + 2;
              if (i<nKeyWords)
                strcat ( S,", " );
            } else
              k = 80;
          } else
            k = 80;
        } while (k<70);
        PadSpaces   ( S,80 );
        f.WriteLine ( S );
        N++;
      }
    }
  }

  void  KeyWords::MakeCIF ( mmcif::PData CIF )  {
  int  i,k;
  char S[500];
    strcpy ( S,"\n" );
    for (i=0;i<nKeyWords;i++)
      if (KeyWord[i])  {
        k = strlen(KeyWord[i]);
        if (strlen(S)+k>70)  {
          CIF->PutString ( S,CIFCAT_STRUCT_KEYWORDS,CIFTAG_TEXT,true );
          if (k>(int)sizeof(S))  {
            CIF->PutString ( KeyWord[i],CIFCAT_STRUCT_KEYWORDS,
                             CIFTAG_TEXT,true );
            k = 0;
          }
          strcpy ( S,"\n" );
        }
        if (k>0) {
          strcat ( S,KeyWord[i] );
          if (i<nKeyWords-1)  strcat ( S,", " );
        }
      }
    if (strlen(S)>1)
      CIF->PutString ( S,CIFCAT_STRUCT_KEYWORDS,CIFTAG_TEXT,true );
  }

  void  KeyWords::GetCIF ( mmcif::PData CIF )  {
  pstr  F;
  int   i,j,k;
  bool  NB;
  char  c;
    Delete();
    F = CIF->GetString ( CIFCAT_STRUCT_KEYWORDS,CIFTAG_TEXT,i );
    k = 0;
    if ((!i) && F)  {
      i  = 0;
      NB = false;
      while (F[i])  {
        if (F[i]==',')  {
          nKeyWords++;
          NB = false;
        } else if (F[i]!=' ')
          NB = true;
        i++;
      }
      if (NB)  nKeyWords++;
      KeyWord = new pstr[nKeyWords];
      i = 0;
      while (F[i] && (k<nKeyWords))  {
        while ((F[i]==' ') || (F[i]=='\n') || (F[i]=='\r'))  i++;
        j = i;
        while (F[i] && (F[i]!=','))  i++;
        c    = F[i];
        F[i] = char(0);
        KeyWord[k] = NULL;
        CreateCopy ( KeyWord[k],&(F[j]) );
        j = 0;
        while (KeyWord[k][j])  {
          if ((KeyWord[k][j]=='\n') || (KeyWord[k][j]=='\r'))
            KeyWord[k][j] = ' ';
          j++;
        }
        F[i] = c;
        k++;
        if (F[i])  i++;
      }
    }
    while (k<nKeyWords)  KeyWord[k++] = NULL;
    CIF->DeleteField ( CIFCAT_STRUCT_KEYWORDS,CIFTAG_TEXT );
  }


  void  KeyWords::Copy ( PKeyWords KeyWords )  {
  int i;
    Delete();
    nKeyWords = KeyWords->nKeyWords;
    if (nKeyWords>0)  {
      KeyWord = new pstr[nKeyWords];
      for (i=0;i<nKeyWords;i++)  {
        KeyWord[i] = NULL;
        CreateCopy ( KeyWord[i],KeyWords->KeyWord[i] );
      }
    }
  }

  void  KeyWords::write ( io::RFile f )  {
  int i;
  byte Version=1;
    f.WriteByte ( &Version   );
    f.WriteInt  ( &nKeyWords );
    for (i=0;i<nKeyWords;i++)
      f.CreateWrite ( KeyWord[i] );
  }

  void  KeyWords::read ( io::RFile f )  {
  int  i;
  byte Version;
    Delete();
    f.ReadByte ( &Version   );
    f.ReadInt  ( &nKeyWords );
    if (nKeyWords>0)  {
      KeyWord = new pstr[nKeyWords];
      for (i=0;i<nKeyWords;i++)  {
        KeyWord[i] = NULL;
        f.CreateRead ( KeyWord[i] );
      }
    }
  }

  MakeStreamFunctions(KeyWords)


  //  ===================  ExpData  ======================

  ExpData::ExpData() : ContString()  {
    InitExpData();
  }

  ExpData::ExpData ( cpstr S ) : ContString()  {
    InitExpData();
    ConvertPDBASCII ( S );
  }

  ExpData::ExpData ( io::RPStream Object ) : ContString(Object)  {
    InitExpData();
  }

  ExpData::~ExpData() {}

  void  ExpData::InitExpData()  {
    CreateCopy ( CIFCategory,CIFCAT_EXPTL  );
    CreateCopy ( CIFTag     ,CIFTAG_METHOD );
  }

  ERROR_CODE ExpData::ConvertPDBASCII ( cpstr S )  {
    if (strlen(S)>10)
         CreateCopy ( Line,&(S[10])  );
    else CreateCopy ( Line,pstr(" ") );
    return Error_NoError;
  }

  void  ExpData::PDBASCIIDump ( pstr S, int N )  {
    if (N==0)  strcpy  ( S,"EXPDTA    " );
         else  sprintf ( S,"EXPDTA  %2i",N+1 );
    strcat ( S,Line );
  }

  void  ExpData::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte ( &Version );
    ContString::write ( f );
  }

  void  ExpData::read ( io::RFile f )  {
  byte Version;
    f.ReadByte ( &Version );
    ContString::read ( f );
  }

  MakeStreamFunctions(ExpData)




  //  ===================  MdlType  ======================

  MdlType::MdlType() : ContString()  {
    InitMdlType();
  }

  MdlType::MdlType ( cpstr S ) : ContString()  {
    InitMdlType();
    ConvertPDBASCII ( S );
  }

  MdlType::MdlType ( io::RPStream Object ) : ContString(Object)  {
    InitMdlType();
  }

  MdlType::~MdlType() {}

  void  MdlType::InitMdlType()  {
    CreateCopy ( CIFCategory,CIFCAT_EXPTL  );
    CreateCopy ( CIFTag     ,CIFTAG_METHOD );
  }

  ERROR_CODE MdlType::ConvertPDBASCII ( cpstr S )  {
    if (strlen(S)>10)
         CreateCopy ( Line,&(S[10])  );
    else CreateCopy ( Line,pstr(" ") );
    return Error_NoError;
  }

  void  MdlType::PDBASCIIDump ( pstr S, int N )  {
    if (N==0)  strcpy  ( S,"MDLTYP    " );
         else  sprintf ( S,"MDLTYP  %2i",N+1 );
    strcat ( S,Line );
  }

  void  MdlType::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte ( &Version );
    ContString::write ( f );
  }

  void  MdlType::read ( io::RFile f )  {
  byte Version;
    f.ReadByte ( &Version );
    ContString::read ( f );
  }

  MakeStreamFunctions(MdlType)


  //  ===================  Author  ======================

  Author::Author() : ContString()  {
    InitAuthor();
  }

  Author::Author ( cpstr S ) : ContString()  {
    InitAuthor();
    ConvertPDBASCII ( S );
  }

  Author::Author ( io::RPStream Object ) : ContString(Object)  {
    InitAuthor();
  }

  Author::~Author() {}

  void  Author::InitAuthor()  {
    CreateCopy ( CIFCategory,CIFCAT_AUDIT_AUTHOR );
    CreateCopy ( CIFTag     ,CIFTAG_NAME         );
  }

  ERROR_CODE Author::ConvertPDBASCII ( cpstr S )  {
    if (strlen(S)>10)
         CreateCopy ( Line,&(S[10])  );
    else CreateCopy ( Line,pstr(" ") );
    return Error_NoError;
  }

  void  Author::PDBASCIIDump ( pstr S, int N )  {
    if (N==0)  strcpy  ( S,"AUTHOR    " );
         else  sprintf ( S,"AUTHOR  %2i",N+1 );
    strcat ( S,Line );
  }

  void  Author::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte ( &Version );
    ContString::write ( f );
  }

  void  Author::read ( io::RFile f )  {
  byte Version;
    f.ReadByte ( &Version );
    ContString::read ( f );
  }

  MakeStreamFunctions(Author)



  //  ================  RevData  ===================

  RevData::RevData() : ContainerClass()  {
    InitRevData();
  }

  RevData::RevData ( cpstr S ) : ContainerClass()  {
    InitRevData();
    ConvertPDBASCII ( S );
  }

  RevData::RevData ( io::RPStream Object ) : ContainerClass(Object)  {
    InitRevData();
  }

  RevData::~RevData() {}

  void  RevData::InitRevData()  {
  int i;
    modNum  = 0;
    strcpy ( modDate,"DD-MMM-YYYY" );
    strcpy ( modId  , "----" );
    modType = -1;
    for (i=0;i<4;i++)
      strcpy ( record[i],"      " );
    Warning = 0;
  }

  void  RevData::PDBASCIIDump ( pstr S, int N )  {
  //  makes the ASCII PDB REVDATA line number N
  //  from the class' data
  int i;
    if (N==0)  sprintf ( S,"REVDAT %3i  " ,modNum     );
         else  sprintf ( S,"REVDAT %3i%2i",modNum,N+1 );
    i = strlen(S);
    while (i<80)
      S[i++] = ' ';
    S[i] = char(0);
    Date11to9 ( modDate,&(S[13]) );
    strncpy   ( &(S[23]),modId,5 );
    S[31] = char(modType+int('0'));
    for (i=0;i<4;i++)
      strncpy ( &(S[39+i*7]),record[i],6 );
  }

  void RevData::MakeCIF ( mmcif::PData CIF, int N )  {
  mmcif::PLoop Loop;
  int         RC,i,j;
  char        DateCIF[20];
    RC = CIF->AddLoop ( CIFCAT_DATABASE_PDB_REV,Loop );
    if ((RC!=mmcif::CIFRC_Ok) || (N==0))  {
      // the category was (re)created, privide tags
      Loop->AddLoopTag ( CIFTAG_NUM                   );
      Loop->AddLoopTag ( CIFTAG_DATE                  );
      Loop->AddLoopTag ( CIFTAG_REPLACES              );
      Loop->AddLoopTag ( CIFTAG_MOD_TYPE              );
      Loop->AddLoopTag ( CIFTAG_RCSB_RECORD_REVISED_1 );
      Loop->AddLoopTag ( CIFTAG_RCSB_RECORD_REVISED_2 );
      Loop->AddLoopTag ( CIFTAG_RCSB_RECORD_REVISED_3 );
      Loop->AddLoopTag ( CIFTAG_RCSB_RECORD_REVISED_4 );
    }
    Date11toCIF ( modDate,DateCIF );
    Loop->AddInteger ( modNum  );
    Loop->AddString  ( DateCIF );
    Loop->AddString  ( modId   );
    Loop->AddInteger ( modType );
    for (i=0;i<4;i++)  {
      j = 0;
      while (record[i][j] && (record[i][j]==' '))  j++;
      if (record[i][j])  Loop->AddString ( record[i] );
                   else  Loop->AddString ( NULL      );
    }

  }

  ERROR_CODE RevData::GetCIF ( mmcif::PData CIF, int & n )  {
  mmcif::PLoop Loop;
  int          RC;
  pstr         F;

    Loop = CIF->GetLoop ( CIFCAT_DATABASE_PDB_REV );
    if (!Loop)  {
      n = -1;
      return Error_EmptyCIF;
    }

    RC = Loop->GetInteger ( modNum,CIFTAG_NUM,n,true );
    if (RC==mmcif::CIFRC_WrongIndex)  {
      n = -1;
      return Error_EmptyCIF;
    }
    if (RC==mmcif::CIFRC_WrongFormat)  {
      sprintf ( CIFErrorLocation,"loop %s.%s row %i",
                CIFCAT_DATABASE_PDB_REV,CIFTAG_NUM,n );
      n = -Error_UnrecognizedInteger-1;
      return Error_UnrecognizedInteger;
    }

    F = Loop->GetString ( CIFTAG_DATE,n,RC );
    if ((!RC) && F)  DateCIFto11 ( F,modDate );
    F = Loop->GetString ( CIFTAG_REPLACES,n,RC );
    if ((!RC) && F)  strcpy ( modId,F );
    RC = Loop->GetInteger ( modType,CIFTAG_MOD_TYPE,n,true );
    if (RC==mmcif::CIFRC_WrongFormat)  {
      sprintf ( CIFErrorLocation,"loop %s.%s row %i",
                CIFCAT_DATABASE_PDB_REV,CIFTAG_MOD_TYPE,n );
      n = -Error_UnrecognizedInteger-1;
      return Error_UnrecognizedInteger;
    }

    F = Loop->GetString ( CIFTAG_RCSB_RECORD_REVISED_1,n,RC );
    if ((!RC) && F)  strcpy ( record[0],F );
    F = Loop->GetString ( CIFTAG_RCSB_RECORD_REVISED_2,n,RC );
    if ((!RC) && F)  strcpy ( record[1],F );
    F = Loop->GetString ( CIFTAG_RCSB_RECORD_REVISED_3,n,RC );
    if ((!RC) && F)  strcpy ( record[2],F );
    F = Loop->GetString ( CIFTAG_RCSB_RECORD_REVISED_4,n,RC );
    if ((!RC) && F)  strcpy ( record[3],F );

    Loop->DeleteField ( CIFTAG_DATE                 ,n );
    Loop->DeleteField ( CIFTAG_REPLACES             ,n );
    Loop->DeleteField ( CIFTAG_RCSB_RECORD_REVISED_1,n );
    Loop->DeleteField ( CIFTAG_RCSB_RECORD_REVISED_2,n );
    Loop->DeleteField ( CIFTAG_RCSB_RECORD_REVISED_3,n );
    Loop->DeleteField ( CIFTAG_RCSB_RECORD_REVISED_4,n );

    n++;

    return Error_NoError;

  }

  ERROR_CODE RevData::ConvertPDBASCII ( cpstr S )  {
  int  i;
  pstr endptr;
  char N[20];
    Warning = 0;
    strncpy   ( N,&(S[7]),3 );
    N[3]    = char(0);
    modNum  = mround(strtod(N,&endptr));
    if (endptr==N)  Warning |= REVDAT_WARN_MODNUM;
    Date9to11 ( &(S[13]),modDate );
    strncpy   ( modId,&(S[23]),5 );
    modId[5] = char(0);
    modType  = int(S[31]) - int('0');
    if (modType>9)  Warning |= REVDAT_WARN_MODTYPE;
    for (i=0;i<4;i++)  {
      strncpy ( record[i],&(S[39+i*7]),6 );
      record[i][6] = char(0);
    }
    return Error_NoError;
  }

  void  RevData::Copy ( PContainerClass RevData )  {
  int i;
    modNum  = PRevData(RevData)->modNum;
    modType = PRevData(RevData)->modType;
    strcpy ( modDate,PRevData(RevData)->modDate );
    strcpy ( modId  ,PRevData(RevData)->modId   );
    for (i=0;i<4;i++)
      strcpy ( record[i],PRevData(RevData)->record[i] );
  }

  void  RevData::write ( io::RFile f )  {
  int  i;
  byte Version=1;
    f.WriteByte  ( &Version );
    f.WriteInt   ( &modNum  );
    f.WriteInt   ( &modType );
    f.WriteWord  ( &Warning );
    f.WriteTerLine ( modDate,false );
    f.WriteTerLine ( modId  ,false );
    for (i=0;i<4;i++)
      f.WriteTerLine ( record[i],false );
  }

  void  RevData::read  ( io::RFile f ) {
  int  i;
  byte Version;
    f.ReadByte  ( &Version );
    f.ReadInt   ( &modNum  );
    f.ReadInt   ( &modType );
    f.ReadWord  ( &Warning );
    f.ReadTerLine ( modDate,false );
    f.ReadTerLine ( modId  ,false );
    for (i=0;i<4;i++)
      f.ReadTerLine ( record[i],false );
  }

  MakeStreamFunctions(RevData)



  //  ================  Supersede  ===================

  Supersede::Supersede() : ContainerClass()  {
    InitSupersede();
  }

  Supersede::Supersede ( cpstr S ) : ContainerClass()  {
    InitSupersede();
    ConvertPDBASCII ( S );
  }

  Supersede::Supersede ( io::RPStream Object ) : ContainerClass(Object)  {
    InitSupersede();
  }

  Supersede::~Supersede() {}

  void  Supersede::InitSupersede()  {
  int i;
    strcpy ( sprsdeDate,"DD-MMM-YYYY" );
    strcpy ( idCode, "----" );
    for (i=0;i<8;i++)
      strcpy ( sIdCode[i],"    " );
  }

  void  Supersede::PDBASCIIDump ( pstr S, int N )  {
  //  makes the ASCII PDB OBSLTE line number N
  //  from the class' data
  int i;
    if (N==0)  strcpy  ( S,"SPRSDE    " );
         else  sprintf ( S,"SPRSDE  %2i",N+1 );
    i = strlen(S);
    while (i<80)
      S[i++] = ' ';
    S[i] = char(0);
    if (N==0)  {
      Date11to9 ( sprsdeDate,&(S[11]) );
      strncpy   ( &(S[21]),idCode,4 );
    }
    for (i=0;i<8;i++)
      strncpy ( &(S[31+5*i]),sIdCode[i],4 );
  }

  void  Supersede::MakeCIF ( mmcif::PData CIF, int )  {
  mmcif::PLoop Loop;
  int         RC,i,j;
  char        DateCIF[20];
    RC = CIF->AddLoop ( CIFCAT_SPRSDE,Loop );
    if (RC!=mmcif::CIFRC_Ok)  {
      // the category was (re)created, privide tags
      Loop->AddLoopTag ( CIFTAG_ID             );
      Loop->AddLoopTag ( CIFTAG_DATE           );
      Loop->AddLoopTag ( CIFTAG_REPLACE_PDB_ID );
      Loop->AddLoopTag ( CIFTAG_PDB_ID         );
    }
    Date11toCIF ( sprsdeDate,DateCIF );
    for (i=0;i<8;i++)  {
      j = 0;
      while (sIdCode[i][j] && (sIdCode[i][j]==' '))  j++;
      if (sIdCode[i][j])  {
        Loop->AddString ( pstr("SPRSDE") );
        Loop->AddString ( DateCIF    );
        Loop->AddString ( idCode     );
        Loop->AddString ( sIdCode[i] );
      }
    }
  }

  ERROR_CODE Supersede::ConvertPDBASCII ( cpstr S )  {
  int i;
    if (S[9]==' ')  {
      Date9to11 ( &(S[11]),sprsdeDate );
      strncpy   ( idCode,&(S[21]),4 );
      idCode[4] = char(0);
    }
    for (i=0;i<8;i++)  {
      strncpy ( sIdCode[i],&(S[31+i*5]),4 );
      sIdCode[i][4] = char(0);
    }
    return Error_NoError;
  }

  ERROR_CODE Supersede::GetCIF ( mmcif::PData CIF, int & n )  {
  mmcif::PLoop Loop;
  int         i,RC;
  pstr        F,FDate,FID;
  char        DateCIF [20];
  char        DateCIF0[20];
  IDCode      idCode1;

    Loop = CIF->GetLoop ( CIFCAT_SPRSDE );
    if (!Loop)  {
      n = -1;  // signal to finish processing of this structure
      return Error_EmptyCIF;
    }
    i = 0;
    do  {
      F = Loop->GetString ( CIFTAG_ID,n,RC );
      if (RC)  {
        if (i==0)  {
          n = -1;
          return Error_EmptyCIF;
        } else
          return Error_NoError;
      }
      if (F)  {
        if (!strcmp(F,"SPRSDE"))  {
          FDate = Loop->GetString ( CIFTAG_DATE,n,RC );
          if ((!RC) && FDate)
                strncpy ( DateCIF,FDate,15 );
          else  strcpy  ( DateCIF,"YYYY-MMM-DD" );
          FID = Loop->GetString ( CIFTAG_REPLACE_PDB_ID,n,RC );
          if ((!RC) && FID)
                strncpy ( idCode1,FID,sizeof(IDCode)-1 );
          else  idCode1[0] = char(0);
          if (i==0)  {
            DateCIFto11 ( DateCIF,sprsdeDate );
            DateCIF[11] = char(0);
            strcpy ( idCode  ,idCode1 );
            strcpy ( DateCIF0,DateCIF );
          } else if ((strcmp(DateCIF0,DateCIF)) ||
                     (strcmp(idCode,idCode1)))
            return Error_NoError;
          FID = Loop->GetString ( CIFTAG_PDB_ID,n,RC );
          if ((!RC) && FID)
               strncpy ( sIdCode[i],FID,sizeof(IDCode)-1 );
          else sIdCode[i][0] = char(0);
          Loop->DeleteField ( CIFTAG_ID            ,n );
          Loop->DeleteField ( CIFTAG_DATE          ,n );
          Loop->DeleteField ( CIFTAG_REPLACE_PDB_ID,n );
          Loop->DeleteField ( CIFTAG_PDB_ID        ,n );
          i++;
        }
      }
      n++;
    } while (i<8);

    return Error_NoError;

  }

  void  Supersede::Copy ( PContainerClass Supersede )  {
  int i;
    strcpy ( sprsdeDate,PSupersede(Supersede)->sprsdeDate );
    strcpy ( idCode    ,PSupersede(Supersede)->idCode     );
    for (i=0;i<8;i++)
      strcpy ( sIdCode[i],PSupersede(Supersede)->sIdCode[i] );
  }

  void  Supersede::write ( io::RFile f )  {
  int  i;
  byte Version=1;
    f.WriteByte  ( &Version );
    f.WriteTerLine ( sprsdeDate,false );
    f.WriteTerLine ( idCode    ,false );
    for (i=0;i<8;i++)
      f.WriteTerLine ( sIdCode[i],false );
  }

  void  Supersede::read  ( io::RFile f ) {
  int  i;
  byte Version;
    f.ReadByte  ( &Version );
    f.ReadTerLine ( sprsdeDate,false );
    f.ReadTerLine ( idCode    ,false );
    for (i=0;i<8;i++)
      f.ReadTerLine ( sIdCode[i],false );
  }

  MakeStreamFunctions(Supersede)


  //  ===================  Journal  ======================

  Journal::Journal() : ContString()  {
    InitJournal();
  }

  Journal::Journal ( cpstr S ) : ContString()  {
    InitJournal();
    ConvertPDBASCII ( S );
  }

  Journal::Journal ( io::RPStream Object ) : ContString(Object)  {
    InitJournal();
  }

  Journal::~Journal() {}

  void  Journal::InitJournal()  {
    CreateCopy ( CIFCategory,CIFCAT_CITATION );
    CreateCopy ( CIFTag     ,CIFTAG_TEXT     );
  }

  ERROR_CODE Journal::ConvertPDBASCII ( cpstr S )  {
    if (strlen(S)>10)
         CreateCopy ( Line,&(S[10]) );
    else CreateCopy ( Line,pstr(" ") );
    return Error_NoError;
  }

  void  Journal::PDBASCIIDump ( pstr S, int )  {
    strcpy ( S,"JRNL      " );
    strcat ( S,Line         );
  }

  void  Journal::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte ( &Version );
    ContString::write ( f );
  }

  void  Journal::read ( io::RFile f )  {
  byte Version;
    f.ReadByte ( &Version );
    ContString::read ( f );
  }

  MakeStreamFunctions(Journal)



  //  ===================  Remark  ======================

  Remark::Remark() : ContainerClass()  {
    InitRemark();
  }

  Remark::Remark ( cpstr S ) : ContainerClass()  {
    InitRemark();
    ConvertPDBASCII ( S );
  }

  Remark::Remark ( io::RPStream Object ) : ContainerClass(Object)  {
    InitRemark();
  }

  Remark::~Remark() {
    if (remark)  delete[] remark;
  }

  void  Remark::InitRemark()  {
    remarkNum = 0;
    remark    = NULL;
  }

  ERROR_CODE Remark::ConvertPDBASCII ( cpstr S )  {
  int i;
    GetInteger ( remarkNum,&(S[7]),3 );
    if (remarkNum==MinInt4)  CreateCopy ( remark,S );
    else if (strlen(S)>11)   CreateCopy ( remark,&(S[11])  );
                       else  CreateCopy ( remark,pstr(" ") );
    i = strlen(remark)-1;
    while ((i>0) && (remark[i]==' '))  i--;
    remark[i+1] = char(0);
    return Error_NoError;
  }

  void  Remark::PDBASCIIDump ( pstr S, int )  {
    if (remarkNum==MinInt4)
      strcpy ( S,remark );
    else  {
      strcpy     ( S,"REMARK" );
      PadSpaces  ( S,80 );
      PutInteger ( &(S[7]) ,remarkNum,3 );
      strncpy    ( &(S[11]),remark,IMin(68,strlen(remark)) );
    }
  }

  void  Remark::MakeCIF ( mmcif::PData CIF, int N )  {
  mmcif::PLoop Loop;
  int         RC;
    RC = CIF->AddLoop ( CIFCAT_NDB_DATABASE_REMARK,Loop );
    if ((RC!=mmcif::CIFRC_Ok) || (N==0))  {
      // the category was (re)created, privide tags
      Loop->AddLoopTag ( CIFTAG_ID   );
      Loop->AddLoopTag ( CIFTAG_TEXT );
    }
    if (remarkNum==MinInt4)  Loop->AddString  ( NULL      );
                       else  Loop->AddInteger ( remarkNum );
    Loop->AddString ( remark );
  }

  ERROR_CODE  Remark::GetCIF ( mmcif::PData CIF, int & n )  {
  mmcif::PLoop Loop;
  int          RC;

    Loop = CIF->GetLoop ( CIFCAT_NDB_DATABASE_REMARK );
    if (!Loop)  {
      n = -1;
      return Error_EmptyCIF;
    }
    if (n>=Loop->GetLoopLength() )  {
      n = -1;
      return Error_EmptyCIF;
    }

    RC = Loop->GetInteger ( remarkNum,CIFTAG_ID,n,true );
    if (RC==mmcif::CIFRC_WrongFormat)  {
      sprintf ( CIFErrorLocation,"loop %s.%s row %i",
                CIFCAT_NDB_DATABASE_REMARK,CIFTAG_ID,n );
      n = -Error_UnrecognizedInteger-1;
      return Error_UnrecognizedInteger;
    } else if (RC)
      remarkNum = MinInt4;
    Loop->GetString ( remark,CIFTAG_TEXT,n,true );

    n++;

    return Error_NoError;

  }

  void  Remark::Copy ( PContainerClass RemarkClass )  {
    remarkNum = PRemark(RemarkClass)->remarkNum;
    CreateCopy ( remark,PRemark(RemarkClass)->remark );
  }

  void  Remark::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte   ( &Version   );
    f.WriteInt    ( &remarkNum );
    f.CreateWrite ( remark     );
  }

  void  Remark::read ( io::RFile f )  {
  byte Version;
    f.ReadByte   ( &Version   );
    f.ReadInt    ( &remarkNum );
    f.CreateRead ( remark     );
  }

  MakeStreamFunctions(Remark)


  //  =================  Biomolecule  =====================

  #define  R350_ERRBIOMT     (-3)
  #define  R350_ERROR        (-2)
  #define  R350_END          (-1)
  #define  R350_NONE           0
  #define  R350_BIOMOLECULE    1
  #define  R350_CHAINS         2
  #define  R350_BIOMT          3

  void getRemarkKey ( RPRemark rem, int & lkey )  {
    if (rem)  {
      if (rem->remarkNum!=350)  lkey = R350_END;
      else if (rem->remark)  {
        if (strcasestr(rem->remark,"BIOMOLECULE:"))
          lkey = R350_BIOMOLECULE;
        else if (strcasestr(rem->remark,"CHAINS:"))
          lkey = R350_CHAINS;
        else if (strcasestr(rem->remark,"BIOMT1") ||
                 strcasestr(rem->remark,"BIOMT2") ||
                 strcasestr(rem->remark,"BIOMT3"))
          lkey = R350_BIOMT;
        else
          lkey = R350_NONE;
      }
    }
  }

  int lookupRemarks ( int & i, RPRemark rem,
                      RTitleContainer Remark )  {
  int l,lkey;

    l    = Remark.Length();
    lkey = R350_NONE;
    while ((i<l) && (lkey==R350_NONE))  {
      getRemarkKey ( rem,lkey );
      if (lkey==R350_NONE)  {
        i++;
        rem = (PRemark)Remark.GetContainerClass ( i );
      }
    }

    return lkey;

  }



  BMApply::BMApply() : io::Stream()  {
    InitBMApply();
  }

  BMApply::BMApply ( io::RPStream Object ) : io::Stream ( Object )  {
    InitBMApply();
  }

  BMApply::~BMApply()  {
    FreeMemory();
  }

  void  BMApply::InitBMApply()  {
    chain     = NULL;
    nChains   = 0;
    tm        = NULL;
    nMatrices = 0;
  }

  void  BMApply::FreeMemory()  {
    if (chain)  delete[] chain;
    if (tm)     delete[] tm;
    chain     = NULL;
    nChains   = 0;
    tm        = NULL;
    nMatrices = 0;
  }

  int  BMApply::addChains ( int & i, RPRemark rem,
                             RTitleContainer Remark )  {
  PChainID ch1;
  pstr     p;
  int      l,lkey,nAdd,j;

    l    = Remark.Length();
    lkey = R350_NONE;

    while ((i<l) && (lkey==R350_NONE))  {

      p = strcasestr ( rem->remark,"CHAINS:" );
      if (p)  p += 7;
      else  {
        p = rem->remark;
        while (*p==' ')  p++;
        if ((p[1]!=',') && (p[1]!=' '))  p = NULL;
      }

      if (p)  {
        nAdd  = strlen(p)/2 + 3;
        ch1   = chain;
        chain = new ChainID[nChains+nAdd];
        for (j=0;j<nChains;j++)
          strcpy ( chain[j],ch1[j] );
        if (ch1)  delete[] ch1;

        while (*p)  {
          while ((*p==' ') || (*p==','))  p++;
          if (*p)  {
            if ((p[1]==',') || (p[1]==' ') || (p[1]==char(0)))  {
              chain[nChains][0] = *p;
              chain[nChains][1] = char(0);
              nChains++;
              p++;
            } else
              break;
          }
        }
      }

      do  {
        i++;
        if (i<l)  {
          rem = (PRemark)Remark.GetContainerClass ( i );
          if (rem)  {
            if (rem->remarkNum!=350)  lkey = R350_END;
            else getRemarkKey ( rem,lkey );
          }
        } else
          lkey = R350_END;
      } while ((!rem) && (lkey==R350_NONE));

    }

    return lkey;

  }

  int getBIOMT ( RPRemark rem, int biomtNo, mat44 & t,
                 RTitleContainer Remark, int & i )  {
  char PN[20];
  pstr p1,p2;
  int  l,j,lkey;

    sprintf ( PN,"BIOMT%1i",biomtNo );
    p1 = strcasestr ( rem->remark,PN );
    if (!p1)  return R350_ERRBIOMT;

    p1 += 6;
    while (*p1==' ')  p1++;
    while (*p1 && (*p1!=' '))  p1++;

    l = biomtNo - 1;
    t[l][0] = strtod ( p1,&p2 );
    if (p1==p2)  return R350_ERRBIOMT;
    t[l][1] = strtod ( p2,&p1 );
    if (p1==p2)  return R350_ERRBIOMT;
    t[l][2] = strtod ( p1,&p2 );
    if (p1==p2)  return R350_ERRBIOMT;
    t[l][3] = strtod ( p2,&p1 );
    if (p1==p2)  return R350_ERRBIOMT;

    if (biomtNo==3)  {
      for (j=0;j<3;j++)
        t[3][j] = 0.0;
      t[3][3] = 1.0;
    }

    l    = Remark.Length();
    lkey = R350_BIOMT;
    do  {
      i++;
      if (i<l)  {
        rem = (PRemark)Remark.GetContainerClass ( i );
        if (rem)  {
          if (rem->remarkNum!=350)  lkey = R350_END;
                              else  getRemarkKey ( rem,lkey );
        }
      } else
        lkey = R350_END;
    } while ((lkey==R350_NONE) || ((!rem) && (lkey==R350_BIOMT)));

    return lkey;

  }

  int  BMApply::addMatrices ( int & i, RPRemark rem,
                               RTitleContainer Remark )  {
  pmat44 tm1;
  int    l,lkey,j,k1,k2,nAlloc;

    l      = Remark.Length();
    lkey   = R350_BIOMT;
    nAlloc = nMatrices;

    while ((i<l) && (lkey==R350_BIOMT))  {

      if (nMatrices>=nAlloc)  {
        nAlloc = nMatrices + 10;
        tm1    = tm;
        tm     = new mat44[nAlloc];
        for (j=0;j<nMatrices;j++)
          for (k1=0;k1<4;k1++)
            for (k2=0;k2<4;k2++)
              tm[j][k1][k2] = tm1[j][k1][k2];
        if (tm1)  delete[] tm1;
      }

      lkey = getBIOMT ( rem,1,tm[nMatrices],Remark,i );
      if (lkey==R350_BIOMT)
        lkey = getBIOMT ( rem,2,tm[nMatrices],Remark,i );
      if (lkey==R350_BIOMT)
        lkey = getBIOMT ( rem,3,tm[nMatrices],Remark,i );
      nMatrices++;

    }

    return lkey;

  }

  void  BMApply::Copy ( PBMApply BMA )  {
  // if BMA is NULL, then empties the class
  int  i,j,k;

    FreeMemory();

    if (BMA)  {

      nChains = BMA->nChains;
      if (nChains>0)  {
        chain = new ChainID[nChains];
        for (i=0;i<nChains;i++)
          strcpy ( chain[i],BMA->chain[i] );
      }

      nMatrices = BMA->nMatrices;
      if (nMatrices>0)  {
        tm = new mat44[nMatrices];
        for (i=0;i<nMatrices;i++)
          for (j=0;j<4;j++)
            for (k=0;k<4;k++)
              tm[i][j][k] = BMA->tm[i][j][k];
       }
    }

  }

  void  BMApply::write ( io::RFile f )  {
  int i,j,k;
    f.WriteInt ( &nChains );
    for (i=0;i<nChains;i++)
      f.WriteTerLine ( chain[i],false );
    f.WriteInt ( &nMatrices );
    for (i=0;i<nMatrices;i++)
      for (j=0;j<3;j++)
        for (k=0;k<4;k++)
          f.WriteReal ( &(tm[i][j][k]) );
  }

  void  BMApply::read ( io::RFile f )  {
  int i,j,k;
    FreeMemory();
    f.ReadInt ( &nChains );
    if (nChains>0)  {
      chain = new ChainID[nChains];
      for (i=0;i<nChains;i++)
        f.ReadTerLine ( chain[i],false );
    }
    f.ReadInt ( &nMatrices );
    if (nMatrices>0)  {
      tm = new mat44[nMatrices];
      for (i=0;i<nMatrices;i++)  {
        for (j=0;j<3;j++)  {
          for (k=0;k<4;k++)
            f.ReadReal ( &(tm[i][j][k]) );
          tm[i][3][j] = 0.0;
        }
        tm[i][3][3] = 1.0;
      }
    }
  }

  MakeStreamFunctions(BMApply)


  Biomolecule::Biomolecule() : io::Stream()  {
    InitBiomolecule();
  }

  Biomolecule::Biomolecule ( io::RPStream Object )
              : io::Stream ( Object )  {
    InitBiomolecule();
  }

  Biomolecule::~Biomolecule()  {
    FreeMemory();
  }

  void  Biomolecule::InitBiomolecule()  {
    bmApply = NULL;
    nBMAs   = 0;
  }

  void  Biomolecule::FreeMemory()  {
  int i;
    if (bmApply)  {
      for (i=0;i<nBMAs;i++)
        if (bmApply[i])  delete bmApply[i];
      delete[] bmApply;
      bmApply = NULL;
    }
    nBMAs = 0;
  }


  PBMApply Biomolecule::addBMApply()  {
  PPBMApply bmA1;
  int       i;
    bmA1 = bmApply;
    bmApply = new PBMApply[nBMAs+1];
    for (i=0;i<nBMAs;i++)
      bmApply[i] = bmA1[i];
    if (bmA1)  delete[] bmA1;
    bmApply[nBMAs] = new BMApply();
    nBMAs++;
    return bmApply[nBMAs-1];
  }

  int Biomolecule::Size()  {
  int i,k;
    k = 0;
    for (i=0;i<nBMAs;i++)
      k += bmApply[i]->nChains*bmApply[i]->nMatrices;
    return k;
  }

  bool Biomolecule::checkComposition ( PChainID chID, ivector occ,
                                           ivector  wocc, int n )  {
  // chID[n] is list of chain IDs
  // occ[n]  is list of chain occurencies
  // wocc[n] is working array
  int     i,j,k,k1;
  bool cmp;

    for (i=0;i<n;i++)
      wocc[i] = 0;

    cmp = true;

    for (i=0;(i<nBMAs) && cmp;i++)
      for (j=0;(j<bmApply[i]->nChains) && cmp;j++)  {
        k1 = -1;
        for (k=0;(k<n) && (k1<0);k++)
          if (!strcmp(chID[k],bmApply[i]->chain[j]))
            k1 = k;
        if (k1<0)  cmp = false;  // chain not found in the list
             else  wocc[k1] += bmApply[i]->nMatrices;
      }

    for (i=0;(i<n) && cmp;i++)
      if (occ[i]!=wocc[i])  cmp = false;

    return cmp;

  }

  void  Biomolecule::Copy ( PBiomolecule B )  {
  // if B is NULL, then empties the class
  int  i;

    FreeMemory();

    if (B)  {

      nBMAs = B->nBMAs;
      if (nBMAs>0)  {
        bmApply = new PBMApply[nBMAs];
        for (i=0;i<nBMAs;i++)
          if (B->bmApply[i])  {
            bmApply[i] = new BMApply();
            bmApply[i]->Copy ( B->bmApply[i] );
          } else
            bmApply[i] = NULL;
      }

    }

  }

  void  Biomolecule::write ( io::RFile f )  {
  int i;
    f.WriteInt ( &nBMAs );
    for (i=0;i<nBMAs;i++)
      StreamWrite ( f,bmApply[i] );
  }

  void  Biomolecule::read ( io::RFile f )  {
  int i;
    FreeMemory();
    f.ReadInt ( &nBMAs );
    if (nBMAs>0)  {
      bmApply = new PBMApply[nBMAs];
      for (i=0;i<nBMAs;i++)  {
        bmApply[i] = NULL;
        StreamRead ( f,bmApply[i] );
      }
    }
  }

  MakeStreamFunctions(Biomolecule)


  //  =====================   Title   =======================

  Title::Title() : io::Stream() {
    Init();
  }

  Title::Title ( io::RPStream Object )  : io::Stream(Object)  {
    Init();
  }

  void  Title::Init()  {

    //  Header data
    classification = NULL;
    depDate[0]     = char(0);
    idCode [0]     = char(0);
    resolution     = -2.0;
    col73          = false;

    biomolecule    = NULL;
    nBiomolecules  = 0;

  }

  Title::~Title() {
    FreeMemory ( false );
  }

  void  Title::FreeMemory ( bool keepBiomolecules )  {

    if (classification)  delete[] classification;
    classification = NULL;
    resolution     = -2.0;

    obsData  .FreeContainer();
    title    .FreeContainer();
    caveat   .FreeContainer();
    compound .FreeContainer();
    source   .FreeContainer();
    keyWords .Delete       ();
    expData  .FreeContainer();
    mdlType  .FreeContainer();
    author   .FreeContainer();
    revData  .FreeContainer();
    supersede.FreeContainer();
    journal  .FreeContainer();
    remark   .FreeContainer();

    col73 = false;

    if (!keepBiomolecules)
      FreeBiomolecules();

  }

  void  Title::FreeBiomolecules()  {
  int  i;
    if (biomolecule)  {
      for (i=0;i<nBiomolecules;i++)
        if (biomolecule[i])  delete biomolecule[i];
      delete[] biomolecule;
      biomolecule = NULL;
    }
    nBiomolecules = 0;
  }


  void  Title::SetHeader ( cpstr Classification,
                                cpstr DepDate,
                                cpstr IDCode )  {
  // fills the PDB file header
    CreateCopy ( classification ,Classification  );
    strncpy    ( depDate,DepDate,sizeof(depDate) );
    strncpy    ( idCode ,IDCode ,sizeof(idCode)  );
    depDate[sizeof(depDate)-1] = char(0);
    idCode [sizeof(idCode) -1] = char(0);
  }

  ERROR_CODE Title::ConvertPDBString ( pstr PDBString ) {
  // Interprets the ASCII PDB line belonging to the title section
  // and fills the corresponding fields.
  //   Returns zero if the line was converted, otherwise returns a
  // non-negative value of Error_XXXX.
  //   PDBString must be not shorter than 81 characters.
  int             i;
  char            c;
  PContainerClass ContainerClass;

    //  pad input line with spaces, if necessary
    PadSpaces ( PDBString,80 );

    if (!strncmp(PDBString,"HEADER",6))  {

      i = 49;
      while ((i>=10) && (PDBString[i]==' '))  i--;
      i++;
      c = PDBString[i];
      PDBString[i] = char(0);
      CreateCopy ( classification,&(PDBString[10]) );
      PDBString[i] = c;

      Date9to11 ( &(PDBString[50]),depDate );

      strncpy ( idCode,&(PDBString[62]),4 );
      idCode[4] = char(0);

    } else if (!strncmp(PDBString,"OBSLTE",6))  {

      ContainerClass = new ObsLine(PDBString);
      obsData.AddData ( ContainerClass );

    } else if (!strncmp(PDBString,"TITLE ",6))  {

      ContainerClass = new TitleLine(PDBString);
      title.AddData ( ContainerClass );

    } else if (!strncmp(PDBString,"CAVEAT",6))  {

      ContainerClass = new Caveat(PDBString);
      caveat.AddData ( ContainerClass );

    } else if (!strncmp(PDBString,"COMPND",6))  {

      ContainerClass = new Compound(PDBString);
      compound.AddData ( ContainerClass );

    } else if (!strncmp(PDBString,"SOURCE",6))  {

      ContainerClass = new Source(PDBString);
      source.AddData ( ContainerClass );

    } else if (!strncmp(PDBString,"KEYWDS",6))  {

      keyWords.ConvertPDBASCII ( PDBString );

    } else if (!strncmp(PDBString,"EXPDTA",6))  {

      ContainerClass = new ExpData(PDBString);
      expData.AddData ( ContainerClass );

    } else if (!strncmp(PDBString,"MDLTYPE",6))  {

      ContainerClass = new MdlType(PDBString);
      mdlType.AddData ( ContainerClass );

    } else if (!strncmp(PDBString,"AUTHOR",6))  {

      ContainerClass = new Author(PDBString);
      author.AddData ( ContainerClass );

    } else if (!strncmp(PDBString,"REVDAT",6))  {

      ContainerClass = new RevData(PDBString);
      revData.AddData ( ContainerClass );

    } else if (!strncmp(PDBString,"SPRSDE",6))  {

      ContainerClass = new Supersede(PDBString);
      supersede.AddData ( ContainerClass );

    } else if (!strncmp(PDBString,"JRNL  ",6))  {

      ContainerClass = new Journal(PDBString);
      journal.AddData ( ContainerClass );

    } else if (!strncmp(PDBString,"REMARK",6))  {

      ContainerClass = new Remark(PDBString);
      remark.AddData ( ContainerClass );

    } else if (!strncmp(PDBString,"SPLIT ",6))  {
      // do nothing at the moment
    } else
      return Error_WrongSection;

    //  check for ID code in columns 73-80

    if (!col73)  {
      if (('0'<=idCode[0]) && (idCode[0]<='9'))  {
        if (!strncasecmp(idCode,&(PDBString[72]),4))
          col73 = true;
      }
    }

    return  Error_NoError;

  }

  realtype Title::GetResolution()  {
  //  returns -1.0 if there is no resolution record in the file
  PRemark rem;
  pstr     p,eptr;
  int      i,l;
    if (resolution>-1.5)  return resolution;
    l = remark.Length();
    for (i=0;(i<l) && (resolution<-1.5);i++)  {
      rem = (PRemark)remark.GetContainerClass ( i );
      if (rem)  {
        if (rem->remarkNum==2)  {
          if (rem->remark)  {
            p = strcasestr ( rem->remark,"RESOLUTION" );
            if (p)  {
              while ((*p) && (*p!=' '))  p++;
              if (*p)  {
                resolution = strtod ( p,&eptr );
                if ((resolution<0.0) || (eptr==p))
                  resolution = -1.0;
              }
            }
          }
        } else if (rem->remarkNum>2)
          resolution = -1.0;
      }
    }
    return resolution;
  }

  PBiomolecule Title::addBiomolecule()  {
  PPBiomolecule  BM1;
  int             i;
    BM1 = biomolecule;
    biomolecule = new PBiomolecule[nBiomolecules+1];
    for (i=0;i<nBiomolecules;i++)
      biomolecule[i] = BM1[i];
    if (BM1)  delete[] BM1;
    biomolecule[nBiomolecules] = new Biomolecule();
    nBiomolecules++;
    return biomolecule[nBiomolecules-1];
  }

  int Title::ParseBiomolecules()  {
  PRemark       rem;
  PBiomolecule  BMol;
  PBMApply      BMA;
  int            i,l, lkey;

    FreeBiomolecules();

    l    = remark.Length();
    i    = 0;
    lkey = 0;
    while ((i<l) && (!lkey))  {
      rem = (PRemark)remark.GetContainerClass ( i );
      if (rem)  {
        if (rem->remarkNum==350)      lkey = 1;
        else if (rem->remarkNum>350)  lkey = -1;
      }
      if (!lkey) i++;
    }

    BMol = NULL;
    BMA  = NULL;

    while (lkey>0)  {

      rem = (PRemark)remark.GetContainerClass ( i );
      lkey = lookupRemarks ( i,rem,remark );

      switch (lkey)  {
        case R350_BIOMOLECULE : BMol = addBiomolecule();
                                i++;
                              break;
        case R350_CHAINS      : if (BMol)  {
                                  BMA = BMol->addBMApply();
                                  while (lkey==R350_CHAINS)
                                    lkey = BMA->addChains(i,rem,remark);
                                } else
                                  lkey = R350_ERROR;
                              break;
        case R350_BIOMT       : if (BMA)
                                  lkey = BMA->addMatrices(i,rem,remark);
                                else
                                  lkey = R350_ERROR;
                              break;
        default : i++;
      }

    }

    if (lkey<=R350_ERROR)  {
      FreeBiomolecules();
      return lkey;
    }

    return nBiomolecules;

  }

  int Title::GetNofBiomolecules()  {
    return nBiomolecules;
  }

  void Title::GetBiomolecules ( PPBiomolecule & BM, int & nBMs ) {
    BM   = biomolecule;
    nBMs = nBiomolecules;
  }

  PBiomolecule Title::GetBiomolecule ( int bmNo )  { // bmno=0,1,..
    if ((0<=bmNo) && (bmNo<nBiomolecules))
      return biomolecule[bmNo];
    return NULL;
  }


  ERROR_CODE Title::GetCIF ( mmcif::PData CIF )  {
  pstr       S;
  ERROR_CODE RC;

    S = NULL;
    CIF->GetDataName ( S,true );
    if (!S) CIF->GetString ( S,CIFCAT_DATABASE,CIFTAG_ENTRY_ID,true );
    if (!S) CIF->GetString ( S,CIFCAT_DATABASE,CIFTAG_CODE_NDB,true );
    if (!S) CIF->GetString ( S,CIFCAT_DATABASE,CIFTAG_CODE_PDB,true );
    if (S)  {
      strncpy ( idCode,S,sizeof(IDCode)-1 );
      idCode[sizeof(IDCode)-1] = char(0);
      delete[] S;
      S = NULL;
      CIF->DeleteField ( CIFCAT_DATABASE,CIFTAG_ENTRY_ID );
      CIF->DeleteField ( CIFCAT_DATABASE,CIFTAG_CODE_NDB );
      CIF->DeleteField ( CIFCAT_DATABASE,CIFTAG_CODE_PDB );
    } else
      idCode[0] = char(0);
    CIF->GetString ( classification,CIFCAT_STRUCT_KEYWORDS,
                                    CIFTAG_NDB_KEYWORDS,true );
    CIF->GetString ( S,CIFCAT_DATABASE,CIFTAG_DATE_ORIGINAL,true );
    if (S)  {
      DateCIFto11 ( S,depDate );
      delete[] S;
      S = NULL;
    } else
      depDate[0] = char(0);

    if (CIF->GetReal(resolution,CIFCAT_REFINE,
                     CIFTAG_LS_D_RES_HIGH,false)!=mmcif::CIFRC_Ok)
      resolution = -2.0;

    obsData .GetCIF ( CIF,ClassID_ObsLine   );
    title   .GetCIF ( CIF,ClassID_TitleLine );
    caveat  .GetCIF ( CIF,ClassID_CAVEAT    );
    compound.GetCIF ( CIF,ClassID_Compound  );
    source  .GetCIF ( CIF,ClassID_Source    );
    keyWords.GetCIF ( CIF );
    expData .GetCIF ( CIF,ClassID_ExpData   );
    mdlType .GetCIF ( CIF,ClassID_MdlType   );
    author  .GetCIF ( CIF,ClassID_Author    );
    RC = revData.GetCIF ( CIF,ClassID_RevData );
    if (RC!=Error_NoError)  {
      supersede.GetCIF ( CIF,ClassID_Supersede );
      journal  .GetCIF ( CIF,ClassID_Journal   );
      RC = remark.GetCIF ( CIF,ClassID_Remark );
    }
    return RC;

  }

  void  Title::MakePDBHeaderString ( pstr PDBString )  {
  //  makes the ASCII PDB HEADER line from the class' data
  int i;

    if (classification)  {

      strcpy ( PDBString,"HEADER    " );
      strcat ( PDBString,classification );
      i = strlen(PDBString);
      while (i<80)
        PDBString[i++] = ' ';
      PDBString[IMin(i,80)] = char(0);
      Date11to9 ( depDate,&(PDBString[50]) );
      strncpy   ( &(PDBString[62]),idCode,4 );

    } else
      strcpy ( PDBString,
        "HEADER    XXXXXXXXXXXXXXXXXXXXXXXXXXXX            XX-XXX-XX   ----" );

  }

  pstr  Title::GetStructureTitle ( pstr & S )  {
  // GetStructureTitle() returns the contents of TITLE record
  // unfolded into single line. If Title is missing, returns
  // contents of COMPND(:MOLECULE). If COMPND is missing, returns
  // HEADER. If Header is missing, returns PDB code. If no PDB
  // code is there, returns "Not available".
  PTitleLine TLine;
  PCompound  CLine;
  pstr        p;
  int         i,cl,l;
  bool     B;

    if (S)  delete[] S;
    S  = NULL;

    cl = title.Length();
    if (cl>0)  {
      l = 0;
      for (i=0;i<cl;i++)  {
        TLine = PTitleLine(title.GetContainerClass(i));
        if (TLine)  l += strlen_des(TLine->Line)+5;
      }
      S = new char[l];
      S[0] = char(0);
      for (i=0;i<cl;i++)  {
        TLine = PTitleLine(title.GetContainerClass(i));
        if (TLine)  {
          if (i>0)  strcat ( S," " );
          strcat_des ( S,TLine->Line );
        }
      }
    } else  {
      cl = compound.Length();
      if (cl>0)  {
        l = 0;
        p = NULL;
        B = true;
        for (i=0;(i<cl) && B;i++)  {
          CLine = PCompound(compound.GetContainerClass(i));
          if (CLine)  {
            if (!p)  {
              p = strstr(CLine->Line,"MOLECULE:");
              if (p)  l += strlen_des(&(p[9]))+5;
            } else  {
              p = strstr(CLine->Line,"MOLECULE:");
              if (p)
                l += strlen_des(&(p[9]))+5;
              else {
                p = FirstOccurence(CLine->Line,':');
                if (!p)  {
                  l += strlen_des(CLine->Line)+5;
                  p = CLine->Line;
                } else
                  B = false;
              }
            }
          }
        }
        if (l>0)  {
          S = new char[l];
          S[0] = char(0);
          p = NULL;
          B = true;
          for (i=0;(i<cl) && B;i++)  {
            CLine = PCompound(compound.GetContainerClass(i));
            if (CLine)  {
              if (!p)  {
                p = strstr(CLine->Line,"MOLECULE:");
                if (p)  strcat_des ( S,&(p[9]) );
              } else  {
                p = strstr(CLine->Line,"MOLECULE:");
                if (p)
                  strcat_des ( S,&(p[9]) );
                else {
                  p = FirstOccurence(CLine->Line,':');
                  if (!p)  {
                    strcat_des ( S,CLine->Line );
                    p = CLine->Line;
                  } else
                    B = false;
                }
              }
              l = strlen(S)-1;
              if (S[l]==';')  S[l] = char(0);
            }
          }
        } else  {
          l = 0;
          for (i=0;i<cl;i++)  {
            CLine = PCompound(compound.GetContainerClass(i));
            if (CLine)  l += strlen_des(CLine->Line)+5;
          }
          S = new char[l];
          S[0] = char(0);
          for (i=0;i<cl;i++)  {
            CLine = PCompound(compound.GetContainerClass(i));
            if (CLine)  {
              if (i>0)  strcat ( S," " );
              strcat_des ( S,CLine->Line );
            }
          }
        }
      } else if (classification)
        CreateCopy ( S,classification );
      else if (idCode[0])
        CreateCopy ( S,idCode );
      else
        CreateCopy ( S,pstr("Not available") );
    }

    if (!S[0])  CreateCopy ( S,pstr("Not available") );

    return S;

  }

  void  Title::PDBASCIIDump ( io::RFile f )  {
  char  PDBString[100];
    if (classification)  {
      MakePDBHeaderString ( PDBString );
      f.WriteLine ( PDBString );
    }
    obsData  .PDBASCIIDump ( f );
    title    .PDBASCIIDump ( f );
    caveat   .PDBASCIIDump ( f );
    compound .PDBASCIIDump ( f );
    source   .PDBASCIIDump ( f );
    keyWords .PDBASCIIDump ( f );
    expData  .PDBASCIIDump ( f );
    mdlType  .PDBASCIIDump ( f );
    author   .PDBASCIIDump ( f );
    revData  .PDBASCIIDump ( f );
    supersede.PDBASCIIDump ( f );
    journal  .PDBASCIIDump ( f );
    remark   .PDBASCIIDump ( f );
  }


  void  Title::MakeCIF ( mmcif::PData CIF )  {
  realtype res;
  char     DateCIF[20];

    if (idCode[0])  {
      CIF->PutDataName ( idCode );
      CIF->PutString   ( idCode, CIFCAT_DATABASE,CIFTAG_ENTRY_ID );
      CIF->PutString   ( idCode, CIFCAT_DATABASE,CIFTAG_CODE_NDB );
      CIF->PutString   ( idCode, CIFCAT_DATABASE,CIFTAG_CODE_PDB );
    } else  {
      CIF->PutDataName ( pstr("")                              );
      CIF->PutString   ( NULL, CIFCAT_DATABASE,CIFTAG_ENTRY_ID );
      CIF->PutString   ( NULL, CIFCAT_DATABASE,CIFTAG_CODE_NDB );
      CIF->PutString   ( NULL, CIFCAT_DATABASE,CIFTAG_CODE_PDB );
    }
    CIF->PutString   ( classification, CIFCAT_STRUCT_KEYWORDS,
                                       CIFTAG_NDB_KEYWORDS );
    if (depDate[0])  {
      Date11toCIF ( depDate,DateCIF );
      CIF->PutString ( DateCIF,CIFCAT_DATABASE,CIFTAG_DATE_ORIGINAL );
    } else
      CIF->PutString ( NULL,CIFCAT_DATABASE,CIFTAG_DATE_ORIGINAL );

    res = GetResolution();
    if (res>=0.0)
         CIF->PutReal   ( res,CIFCAT_REFINE,CIFTAG_LS_D_RES_HIGH,3 );
    else CIF->PutNoData ( mmcif::CIF_NODATA_QUESTION,
                          CIFCAT_REFINE,CIFTAG_LS_D_RES_HIGH );

    obsData  .MakeCIF ( CIF );
    title    .MakeCIF ( CIF );
    caveat   .MakeCIF ( CIF );
    compound .MakeCIF ( CIF );
    source   .MakeCIF ( CIF );
    keyWords .MakeCIF ( CIF );
    expData  .MakeCIF ( CIF );
    mdlType  .MakeCIF ( CIF );
    author   .MakeCIF ( CIF );
    revData  .MakeCIF ( CIF );
    supersede.MakeCIF ( CIF );
    journal  .MakeCIF ( CIF );
    remark   .MakeCIF ( CIF );

  }

  void  Title::Copy ( PTitle TS )  {
  int  i;

    FreeBiomolecules();

    if (TS)  {

      CreateCopy ( classification,TS->classification );
      strcpy     ( depDate       ,TS->depDate        );
      strcpy     ( idCode        ,TS->idCode         );
      resolution = TS->resolution;

      obsData  .Copy ( &(TS->obsData)   );
      title    .Copy ( &(TS->title)     );
      caveat   .Copy ( &(TS->caveat)    );
      compound .Copy ( &(TS->compound)  );
      source   .Copy ( &(TS->source)    );
      keyWords .Copy ( &(TS->keyWords)  );
      expData  .Copy ( &(TS->expData)   );
      mdlType  .Copy ( &(TS->mdlType)   );
      author   .Copy ( &(TS->author)    );
      revData  .Copy ( &(TS->revData)   );
      supersede.Copy ( &(TS->supersede) );
      journal  .Copy ( &(TS->journal)   );
      remark   .Copy ( &(TS->remark)    );

      nBiomolecules = TS->nBiomolecules;
      if (nBiomolecules>0)  {
        biomolecule = new PBiomolecule[nBiomolecules];
        for (i=0;i<nBiomolecules;i++)
          if (TS->biomolecule[i])  {
            biomolecule[i] = new Biomolecule();
            biomolecule[i]->Copy ( TS->biomolecule[i] );
          } else
            biomolecule[i] = NULL;
      }

    } else  {

      if (classification)  delete[] classification;
      classification = NULL;
      resolution     = -2.0;
      obsData  .FreeContainer();
      title    .FreeContainer();
      caveat   .FreeContainer();
      compound .FreeContainer();
      source   .FreeContainer();
      keyWords .Delete       ();
      expData  .FreeContainer();
      mdlType  .FreeContainer();
      author   .FreeContainer();
      revData  .FreeContainer();
      supersede.FreeContainer();
      journal  .FreeContainer();
      remark   .FreeContainer();

    }

  }

  void  Title::TrimInput ( pstr PDBString )  {
    if (col73)  {
      if (!strncasecmp(idCode,&(PDBString[72]),4))
        PDBString[72] = char(0);
    }
    PadSpaces ( PDBString,80 );
  }

  void  Title::write ( io::RFile f )  {
  // writes header to PDB binary file
  int  i;
  byte Version=3;

    f.WriteByte    ( &Version       );

    //  Header data
    f.CreateWrite  ( classification );
    f.WriteTerLine ( depDate,false  );
    f.WriteTerLine ( idCode ,false  );
    f.WriteReal    ( &resolution    );

    obsData  .write ( f );  //  Obsoletion data
    title    .write ( f );  //  Title
    caveat   .write ( f );  //  Error data
    compound .write ( f );  //  Compound
    source   .write ( f );  //  Source
    keyWords .write ( f );  //  Key words
    expData  .write ( f );  //  Experimental data
    mdlType  .write ( f );  //  Model descriptions
    author   .write ( f );  //  Author data
    revData  .write ( f );  //  Revision data
    supersede.write ( f );  //  Supersede records
    journal  .write ( f );  //  Journal records
    remark   .write ( f );  //  Remarks

    f.WriteInt ( &nBiomolecules );
    for (i=0;i<nBiomolecules;i++)
      StreamWrite ( f,biomolecule[i] );

  }

  void  Title::read ( io::RFile f )  {
  // reads header from PDB binary file
  int  i;
  byte Version;

    f.ReadByte    ( &Version );

    //  Header data
    f.CreateRead  ( classification );
    f.ReadTerLine ( depDate,false  );
    f.ReadTerLine ( idCode ,false  );
    if (Version>1)
      f.ReadReal  ( &resolution    );
    else
      resolution = -2.0;

    obsData  .read ( f );   //  Obsoletion data
    title    .read ( f );   //  Title
    caveat   .read ( f );   //  Error data
    compound .read ( f );   //  Compound
    source   .read ( f );   //  Source
    keyWords .read ( f );   //  Key words
    expData  .read ( f );   //  Experimental data
    if (Version>2)
      mdlType.read ( f );   //  Model descriptions
    author   .read ( f );   //  Author data
    revData  .read ( f );   //  Revision data
    supersede.read ( f );   //  Supersede records
    journal  .read ( f );   //  Journal records
    remark   .read ( f );   //  Remarks

    FreeBiomolecules();
    if (Version>1)  {
      f.ReadInt ( &nBiomolecules );
      if (nBiomolecules>0)  {
        biomolecule = new PBiomolecule[nBiomolecules];
        for (i=0;i<nBiomolecules;i++)  {
          biomolecule[i] = NULL;
          StreamRead ( f,biomolecule[i] );
        }
      }
    }

  }

  MakeStreamFunctions(Title)

}  // namespace mmdb


// ===================================================================

/*
void  TestHeader()  {
PTitle  Hdr;
char         S[81],S1[81];

  Hdr = new Title();

  Hdr->SetHeader ( pstr("MUSCLE PROTEIN"),pstr("02-JUN-1993"),pstr("1MYS") );
  Hdr->MakePDBHeaderString ( S );
  printf ( "1234567890123456789012345678901234567890"
           "1234567890123456789012345678901234567890\n" );
  printf ( S );
  printf ( "\n" );

  strcpy ( S,
// 1234567890123456789012345678901234567890123456789012345678901234567890
  "HEADER    HYDROLASE (CARBOXYLIC ESTER)            07-APR-01   2PHI" );

  Hdr->ConvertPDBString ( S );
  Hdr->MakePDBHeaderString    ( S1 );
  printf ( "1234567890123456789012345678901234567890"
           "1234567890123456789012345678901234567890\n" );
  printf ( S1 );
  printf ( "\n" );

  Hdr->SetHeader (
     pstr("MUSCLE PROTEIN;**A VERY LONG TITLE TEST;**ARBITRARY LENGTH"),
     pstr("02-JUN-1993"),pstr("1MYS") );
  Hdr->MakePDBHeaderString ( S );
  printf ( "1234567890123456789012345678901234567890"
           "1234567890123456789012345678901234567890\n" );
  printf ( S );
  printf ( "\n" );

  delete Hdr;

  printf ( " header deleted \n" );

}

void  TestTitle() {
// reads PDB title from file 'in.title'
// and rewrites it into 'out.title' and 'abin.title'
CFile        f;
char         S[81];
PTitle  Title;

  Title = new Title();

  f.assign ( pstr("in.title"),true );
  if (f.reset()) {
    while (!f.FileEnd()) {
      f.ReadLine ( S,sizeof(S) );
      Title->ConvertPDBString ( S );
    }
    f.shut();
  } else {
    printf ( " Can't open input file 'in.title' \n" );
    delete Title;
    return;
  }

  f.assign ( pstr("out.title"),true );
  if (f.rewrite()) {
    Title->PDBASCIIDump ( f );
    f.shut();
  } else {
    printf ( " Can't open output file 'out.title' \n" );
    delete Title;
    return;
  }



  f.assign ( pstr("mmdb.title.bin"),false );
  if (f.rewrite()) {
    Title->write ( f );
    f.shut();
  } else {
    printf ( "  Can't open binary file for writing.\n" );
    delete Title;
    return;
  }

  delete Title;
  printf ( "   Title deleted.\n" );

  Title = new Title();
  if (f.reset()) {
    Title->read ( f );
    f.shut();
  } else {
    printf ( "  Can't open binary file for reading.\n" );
    delete Title;
    return;
  }

  f.assign ( pstr("abin.title"),true );
  if (f.rewrite()) {
    Title->PDBASCIIDump ( f );
    f.shut();
  } else {
    printf ( " Can't open output file 'abin.title' \n" );
  }

  delete Title;

}


*/
