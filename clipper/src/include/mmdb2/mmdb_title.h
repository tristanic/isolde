//  $Id: mmdb_title.h $
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
//    12.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  MMDB_Title <interface>
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
//   (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#ifndef __MMDB_Title__
#define __MMDB_Title__

#include "mmdb_io_stream.h"
#include "mmdb_defs.h"
#include "mmdb_utils.h"
#include "mmdb_mmcif_.h"

namespace mmdb  {

  //  ======================  TitleContainer  =======================

  DefineClass(TitleContainer);
  DefineStreamFunctions(TitleContainer);

  class TitleContainer : public ClassContainer  {

    public :

      TitleContainer () : ClassContainer() {}
      TitleContainer ( io::RPStream Object )
                        : ClassContainer ( Object ) {}
      ~TitleContainer() {}

      PContainerClass MakeContainerClass ( int ClassID );

  };


  //  ==================  ObsLine  ========================

  DefineClass(ObsLine);
  DefineStreamFunctions(ObsLine);

  class ObsLine : public ContainerClass  {

    public :

      Date   repDate;    //  date of replacement
      IDCode idCode;     //  ID code of replaced entry
      IDCode rIdCode[8]; //  ID codes of entries that replaced this one

      ObsLine ();
      ObsLine ( cpstr S );
      ObsLine ( io::RPStream Object );
      ~ObsLine();

      void       PDBASCIIDump    ( pstr S, int N   );
      void       MakeCIF         ( mmcif::PData CIF, int N );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      ERROR_CODE GetCIF          ( mmcif::PData CIF, int & n );
      CLASS_ID   GetClassID      () { return ClassID_ObsLine; }

      void  Copy  ( PContainerClass ObsLine );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitObsLine();

  };


  //  ====================  TitleLine  =====================

  DefineClass(TitleLine);
  DefineStreamFunctions(TitleLine);

  class TitleLine : public ContString  {

    public :

      TitleLine ();
      TitleLine ( cpstr S );
      TitleLine ( io::RPStream Object );
      ~TitleLine();

      ERROR_CODE ConvertPDBASCII ( cpstr S );
      void       PDBASCIIDump    ( pstr S, int N );
      bool       PDBASCIIDump1   ( io::RFile ) { return false; }
      CLASS_ID   GetClassID      () { return ClassID_TitleLine; }

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitTitleLine();

  };


  //  ====================  Caveat  =====================

  DefineClass(Caveat);
  DefineStreamFunctions(Caveat);

  class Caveat : public ContString  {

    public :

      IDCode idCode;   //  ID code of the entry

      Caveat ();
      Caveat ( cpstr S );
      Caveat ( io::RPStream Object );
      ~Caveat();

      void       PDBASCIIDump    ( pstr S, int N );
      bool       PDBASCIIDump1   ( io::RFile ) { return false; }
      void       MakeCIF         ( mmcif::PData CIF, int N );
      ERROR_CODE ConvertPDBASCII ( cpstr S );

//      virtual void  GetCIF1      ( mmcif::PData CIF, ERROR_CODE & Signal,
//                                   int & pos );

      CLASS_ID   GetClassID      () { return ClassID_CAVEAT; }

      void  Copy  ( PContainerClass Caveat );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitCaveat();

  };


  //  ====================  Compound  =====================

  DefineClass(Compound);
  DefineStreamFunctions(Compound);

  class Compound : public ContString  {

    public :

      Compound ();
      Compound ( cpstr S );
      Compound ( io::RPStream Object );
      ~Compound();

      void       PDBASCIIDump    ( pstr S, int N );
      bool       PDBASCIIDump1   ( io::RFile  ) { return false; }
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      CLASS_ID   GetClassID      () { return ClassID_Compound; }

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitCompound();

  };


  //  ====================  Source  =====================

  DefineClass(Source);
  DefineStreamFunctions(Source);

  class Source : public ContString  {

    public :

      Source ();
      Source ( cpstr S );
      Source ( io::RPStream Object );
      ~Source();

      void       PDBASCIIDump    ( pstr S, int N );
      bool       PDBASCIIDump1   ( io::RFile  ) { return false; }
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      CLASS_ID   GetClassID      () { return ClassID_Source; }

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitSource();

  };


  //  ====================  KeyWords  =====================

  DefineClass(KeyWords);
  DefineStreamFunctions(KeyWords);

  class KeyWords : public io::Stream  {

    public :
      int      nKeyWords;     // number of key words
      psvector KeyWord;       // key word array

      KeyWords ();
      KeyWords ( cpstr S );
      KeyWords ( io::RPStream Object );
      ~KeyWords();

      void  Delete          ();

      void  PDBASCIIDump    ( io::RFile f );
      void  MakeCIF         ( mmcif::PData CIF );

      int   ConvertPDBASCII ( cpstr S );
      void  GetCIF          ( mmcif::PData CIF );

      void  Copy  ( PKeyWords KeyWords );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :
      bool Cont;

      void  Init();

  };


  //  ====================  ExpData  =====================

  DefineClass(ExpData);
  DefineStreamFunctions(ExpData);

  class ExpData : public ContString  {

    public :

      ExpData ();
      ExpData ( cpstr S );
      ExpData ( io::RPStream Object );
      ~ExpData();

      void       PDBASCIIDump  ( pstr S, int N );
      bool       PDBASCIIDump1 ( io::RFile ) { return false; }
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      CLASS_ID   GetClassID      () { return ClassID_ExpData; }

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitExpData();

  };


  //  ====================  MdlType  =====================

  DefineClass(MdlType);
  DefineStreamFunctions(MdlType);

  class MdlType : public ContString  {

    public :

      MdlType ();
      MdlType ( cpstr S );
      MdlType ( io::RPStream Object );
      ~MdlType();

      void       PDBASCIIDump  ( pstr S, int N );
      bool       PDBASCIIDump1 ( io::RFile ) { return false; }
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      CLASS_ID   GetClassID      () { return ClassID_MdlType; }

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitMdlType();

  };


  //  ====================  Author  =====================

  DefineClass(Author);
  DefineStreamFunctions(Author);

  class Author : public ContString  {

    public :

      Author ();
      Author ( cpstr S );
      Author ( io::RPStream Object );
      ~Author();

      void       PDBASCIIDump  ( pstr S, int N   );
      bool       PDBASCIIDump1 ( io::RFile ) { return false; }
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      CLASS_ID   GetClassID      () { return ClassID_Author; }

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitAuthor();

  };


  //  ====================  RevData  =====================

  DefineClass(RevData);
  DefineStreamFunctions(RevData);

  enum REVDAT_WARNING  {
    REVDAT_WARN_MODNUM  = 0x00000001,
    REVDAT_WARN_MODTYPE = 0x00000002
  };

  class RevData : public ContainerClass  {

    public :
      int     modNum;
      Date    modDate;
      char    modId[13];
      int     modType;
      RecName record[4];
      word    Warning;

      RevData ();
      RevData ( cpstr S );
      RevData ( io::RPStream Object );
      ~RevData();

      void       PDBASCIIDump    ( pstr S, int N );
      void       MakeCIF         ( mmcif::PData CIF, int N );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      ERROR_CODE GetCIF          ( mmcif::PData CIF, int & n );
      CLASS_ID   GetClassID      () { return ClassID_RevData; }

      void  Copy  ( PContainerClass RevData );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitRevData();

  };


  //  ==================  Supersede  ========================

  DefineClass(Supersede);
  DefineStreamFunctions(Supersede);

  class Supersede : public ContainerClass  {

    public :
      Date   sprsdeDate;  //  date of supersede
      IDCode idCode;      //  ID code of the entry
      IDCode sIdCode[8];  //  ID codes of superseded entries

      Supersede ();
      Supersede ( cpstr S );
      Supersede ( io::RPStream Object );
      ~Supersede();

      void       PDBASCIIDump    ( pstr S, int N );
      void       MakeCIF         ( mmcif::PData CIF, int N );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      ERROR_CODE GetCIF          ( mmcif::PData CIF, int & n );
      CLASS_ID   GetClassID      () { return ClassID_Supersede; }

      void  Copy  ( PContainerClass Supersede );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitSupersede();

  };


  //  ====================  Journal  =====================

  DefineClass(Journal);
  DefineStreamFunctions(Journal);

  class Journal : public ContString  {

    public :

      Journal ();
      Journal ( cpstr S );
      Journal ( io::RPStream Object );
      ~Journal();

      void       PDBASCIIDump    ( pstr S, int N );
      bool       PDBASCIIDump1   ( io::RFile ) { return false; }
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      CLASS_ID   GetClassID      () { return ClassID_Journal; }

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitJournal();

  };


  //  ====================  Remark  =====================

  DefineClass(Remark);
  DefineStreamFunctions(Remark);

  class Remark : public ContainerClass  {

    public :

      int  remarkNum;  // remark id
      pstr remark;     // remark line

      Remark ();
      Remark ( cpstr S );
      Remark ( io::RPStream Object );
      ~Remark();

      void       PDBASCIIDump    ( pstr S, int N );
      void       MakeCIF         ( mmcif::PData CIF, int N );
      ERROR_CODE ConvertPDBASCII ( cpstr S );
      ERROR_CODE GetCIF          ( mmcif::PData CIF, int & n );
      CLASS_ID   GetClassID      () { return ClassID_Remark; }

      void  Copy  ( PContainerClass RemarkClass );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      void InitRemark();

  };

  //  =================  Biomolecule  =====================

  DefineClass(BMApply);
  DefineStreamFunctions(BMApply);

  class BMApply : public io::Stream  {

    public :
      PChainID  chain;
      int       nChains;
      pmat44    tm;
      int       nMatrices;

      BMApply ();
      BMApply ( io::RPStream Object );
      ~BMApply();

      void  FreeMemory();

      int   addChains ( int & i, RPRemark rem, RTitleContainer Remark );
      int addMatrices ( int & i, RPRemark rem, RTitleContainer Remark );

      void  Copy  ( PBMApply BMA );  // if BMA is NULL, then empties
                                      // the class
      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :
      void  InitBMApply();

  };


  DefineClass(Biomolecule);
  DefineStreamFunctions(Biomolecule);

  class Biomolecule : public io::Stream  {

    public :
      PPBMApply bmApply;
      int       nBMAs;

      Biomolecule ();
      Biomolecule ( io::RPStream Object );
      ~Biomolecule();

      void  FreeMemory();

      PBMApply addBMApply();

      int   Size();
      bool  checkComposition ( PChainID chID, ivector occ,
                               ivector  wocc, int n );

      void  Copy  ( PBiomolecule B );  // if B is NULL, then empties
                                        // the class
      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :
      void  InitBiomolecule();

  };

  //  =================  Title  =======================

  DefineClass(Title);
  DefineStreamFunctions(Title);

  class Title : public io::Stream  {

    friend class Model;
    friend class Chain;
    friend class Root;

    public :

      Title ();
      Title ( io::RPStream Object );
      ~Title();

      void  FreeMemory ( bool keepBiomolecules );

      // Fills the PDB file header
      void  SetHeader ( cpstr Classification, // any length is Ok
                        cpstr DepDate,    // DD-MMM-YYYY
                        cpstr ID_Code );  // not more than 11 chars

      // Interprets the ASCII PDB line belonging to the title section
      // and fills the corresponding fields.
      //   Returns zero if the line was converted, otherwise returns a
      // non-negative value of Error_XXXX.
      //   PDBString must be not shorter than 81 characters.
      ERROR_CODE ConvertPDBString ( pstr PDBString );

      // MakePDBString() makes the ASCII PDB HEADER line from the
      // class data. PDBString must be not shorter than 81 characters.
      void  MakePDBHeaderString ( pstr PDBString );

      // GetStructureTitle() returns the contents of TITLE record
      // unfolded into single line. If Title is missing, returns
      // contents of COMPND(:MOLECULE). If COMPND is missing, returns
      // HEADER. If Header is missing, returns PDB code. If no PDB
      // code is there, returns "Not available".
      pstr  GetStructureTitle ( pstr & S );

      PTitleContainer GetObsData () { return &obsData;  }
      PTitleContainer GetCaveat  () { return &caveat;   }
      PTitleContainer GetCompound() { return &compound; }
      PTitleContainer GetSource  () { return &source;   }
      PKeyWords       GetKeyWords() { return &keyWords; }
      PTitleContainer GetExpData () { return &expData;  }
      PTitleContainer GetMdlType () { return &mdlType;  }
      PTitleContainer GetRemarks () { return &remark;   }
      PTitleContainer GetJournal () { return &journal;  }

      realtype GetResolution(); // -1.0 mean no resolution record in file

      int   ParseBiomolecules(); // returns the number of biomolecules,
                                 // -2 for general format error
                                 // -3 for errors in BIOMT records

      int   GetNofBiomolecules();
      void  GetBiomolecules   ( PPBiomolecule & BM, int & nBMs );
      PBiomolecule GetBiomolecule ( int bmNo ); // bmno=0,1,..
                                 // returns NULL if bmNo is incorrect

      void  PDBASCIIDump ( io::RFile      f   );
      void  MakeCIF      ( mmcif::PData CIF );

      //   GetCIF(..) returns the same code as ConvertPDBString(..)
      // save for Error_WrongSection
      ERROR_CODE  GetCIF ( mmcif::PData CIF );

      inline pstr  GetIDCode() { return idCode; }
      inline bool  GetCol73 () { return col73;  }
      void  TrimInput ( pstr PDBString );

      void  Copy  ( PTitle TS );  // if TS is NULL, then empties
                                       // the class

      void  write ( io::RFile f );    // writes header to PDB binary file
      void  read  ( io::RFile f );    // reads header from PDB binary file

    protected :

      //   Header data
      pstr     classification;  // classification of the molecule
      Date     depDate;         // deposition date DD-MMM-YYYY
      IDCode   idCode;          // unique PDB identifier
      realtype resolution;      // resolution
      bool     col73;           // True if columns 73-80 contain PDB ID

      TitleContainer obsData;     // obsoletion data
      TitleContainer title;       // title data
      TitleContainer caveat;      // error data
      TitleContainer compound;    // compound data
      TitleContainer source;      // source
      KeyWords       keyWords;    // key words
      TitleContainer expData;     // experimental data
      TitleContainer mdlType;     // model descriptions
      TitleContainer author;      // author data
      TitleContainer revData;     // revision data
      TitleContainer supersede;   // supersede records
      TitleContainer journal;     // journal records
      TitleContainer remark;      // remark records

      PPBiomolecule  biomolecule;
      int            nBiomolecules;

      void  Init();
      void  FreeBiomolecules();

      PBiomolecule addBiomolecule();

  };

  extern void  TestHeader();
  extern void  TestTitle (); // reads PDB title from file 'in.title'
                             // and rewrites it into 'out.title' and
                             // 'abin.title'

}  // namespace mmdb


#endif

