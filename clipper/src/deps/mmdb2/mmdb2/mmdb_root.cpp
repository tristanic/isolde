//  $Id: mmdb_root.h $
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
//    13.07.17   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  MMDB_Root <implementation>
//       ~~~~~~~~~
//       Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::Root
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2000-2016
//
//  =================================================================
//

#include "string.h"
#include "stdlib.h"

#include "mmdb_root.h"
#include "mmdb_atom.h"
#include "mmdb_mmcif_.h"
#include "mmdb_cifdefs.h"
#include "mmdb_tables.h"
#include "mmdb_defs.h"


namespace mmdb  {

  //  =====================   Root   =======================

  Root::Root() : UDData()  {
    InitMMDBRoot();
  }

  Root::Root ( io::RPStream Object ) : UDData(Object)  {
    InitMMDBRoot();
  }

  Root::~Root()  {
    FreeFileMemory();
  }

  void  Root::InitMMDBRoot()  {
    nModels = 0;
    model   = NULL;
    nAtoms  = 0;
    atmLen  = 0;
    atom    = NULL;
    CIF     = NULL;
    crModel = NULL;
    crChain = NULL;
    crRes   = NULL;
    lcount  = 0;
    strcpy ( S,"" );
  //  Flags   = 0x00000000;           // no special effects
    Flags   = MMDBF_IgnoreElement;  // done at request for default
    FType   = MMDB_FILE_Undefined;  // undefined file operation
    Exclude = true;
    ignoreRemarks     = false;  // used temporarily
    allowDuplChID     = false;  // used temporarily
    enforceUniqueChID = false;  // used temporarily
    modelCnt          = 0;      // used only at reading files
  }


  void  Root::FreeCoordMemory()  {
    //int i;

  /*
    //   All atoms are kept in array Atom. Models, chains
    // and residues have only references to Atom and
    // they do not dispose Atoms when disposed themselves.
    //   It is important, however, to dispose Atom at
    // still alive residues, because each atom wipes out
    // reference to itself from the corresponding residue
    // before it dies.
    if (Atom)  {
      for (i=0;i<atmLen;i++)
        if (atom[i]) delete atom[i];
      delete Atom;
    }
    Atom    = NULL;
    atmLen  = 0;
    nAtoms  = 0;
  */
    DeleteAllModels();
    if (model)  delete[] model;
    model   = NULL;
    nModels = 0;

    crModel = NULL;
    crChain = NULL;
    crRes   = NULL;

    if (atom)  delete[] atom;

    atom    = NULL;
    atmLen  = 0;
    nAtoms  = 0;

    modelCnt = 0;

  }

  void  Root::FreeFileMemory()  {

    FreeCoordMemory  ();
    title.FreeMemory ( false );
    cryst.FreeMemory ();

    SA      .FreeContainer();
    Footnote.FreeContainer();
    SB      .FreeContainer();
    SC      .FreeContainer();

    if (CIF)  delete CIF;
    CIF = NULL;

    lcount = 0;
    S[0]   = char(0);

  }

  // virtual to be served by MMDB manager classes
  void Root::ResetManager() {
    cryst.Reset();
  }

  void Root::SetFlag ( word Flag )  {
    Flags |= Flag;
    ignoreSegID            = (Flags & MMDBF_IgnoreSegID            ) != 0;
    ignoreElement          = (Flags & MMDBF_IgnoreElement          ) != 0;
    ignoreCharge           = (Flags & MMDBF_IgnoreCharge           ) != 0;
    ignoreNonCoorPDBErrors = (Flags & MMDBF_IgnoreNonCoorPDBErrors ) != 0;
    ignoreUnmatch          = (Flags & MMDBF_IgnoreUnmatch          ) != 0;
    allowDuplChID          = (Flags & MMDBF_AllowDuplChainID       ) != 0;
    enforceUniqueChID      = (Flags & MMDBF_EnforceUniqueChainID   ) != 0;
    cryst.processSG        = (Flags & MMDBF_DoNotProcessSpaceGroup ) == 0;
    cryst.fixSpaceGroup    = (Flags & MMDBF_FixSpaceGroup          ) != 0;
  }

  void Root::RemoveFlag ( word Flag )  {
    Flags &= ~Flag;
    ignoreSegID            = (Flags & MMDBF_IgnoreSegID            ) != 0;
    ignoreElement          = (Flags & MMDBF_IgnoreElement          ) != 0;
    ignoreCharge           = (Flags & MMDBF_IgnoreCharge           ) != 0;
    ignoreNonCoorPDBErrors = (Flags & MMDBF_IgnoreNonCoorPDBErrors ) != 0;
    ignoreUnmatch          = (Flags & MMDBF_IgnoreUnmatch          ) != 0;
    allowDuplChID          = (Flags & MMDBF_AllowDuplChainID       ) != 0;
    enforceUniqueChID      = (Flags & MMDBF_EnforceUniqueChainID   ) != 0;
    cryst.processSG        = (Flags & MMDBF_DoNotProcessSpaceGroup ) == 0;
    cryst.fixSpaceGroup    = (Flags & MMDBF_FixSpaceGroup          ) != 0;
  }


  ERROR_CODE Root::ReadPDBASCII1 ( cpstr PDBLFName, io::GZ_MODE gzipMode )  {
  pstr FName;
    FName = getenv ( PDBLFName );
    if (FName)  return ReadPDBASCII ( FName,gzipMode );
          else  return Error_NoLogicalName;
  }

  void Root::ReadPDBLine ( io::RFile f, pstr L, int maxlen )  {
  int     i;
  bool Done;
    do {
      f.ReadLine ( L,maxlen );
      Done = true;
      if (ignoreRemarks)  {
        if (!strncasecmp(L,"REMARK",6))  Done = false;
      }
      if (Flags & MMDBF_IgnoreBlankLines)  {
        i = 0;
        while (L[i] && (L[i]==' '))  i++;
        if (!L[i])  Done = false;
      }
      if ((Flags & MMDBF_IgnoreHash) && (L[0]=='#'))
        Done = false;
    } while ((!f.FileEnd()) && (!Done));
    PadSpaces  ( L,80 );
  }

  ERROR_CODE Root::ReadPDBASCII ( cpstr PDBFileName, io::GZ_MODE gzipMode )  {
  io::File   f;
  ERROR_CODE RC;

    //  open the file as ASCII for reading
    //  opening it in pseudo-binary mode helps reading various
    //  line terminators for files coming from different platforms
    f.assign ( PDBFileName,false,false,gzipMode );

    if (f.reset(true)) {

      RC = ReadPDBASCII ( f );
      f.shut();

    } else  {

      RC =  Error_CantOpenFile;
      ResetManager  ();
      FreeFileMemory();
      FType = MMDB_FILE_PDB;

    }

    return RC;

  }


  ERROR_CODE Root::ReadPDBASCII ( io::RFile f )  {
  PContString contString;
  word        cleanKey;
  int         modNum;
  bool        fend;
  ERROR_CODE  RC;

    //  remove previous data
    ResetManager  ();
    FreeFileMemory();

    FType = MMDB_FILE_PDB;
    SetFlag ( 0 );

    if (f.FileEnd())  return Error_EmptyFile;

    lcount = 1;  // line counter

    // read title section
    RC = Error_NoError;
    ReadPDBLine ( f,S,sizeof(S) );
    if (Flags & MMDBF_EnforceSpaces)  EnforceSpaces ( S );
    do  {
      if (!strncmp(S,"FTNOTE",6))  {
        contString = new ContString(S);
        Footnote.AddData ( contString );
      } else  {
        RC = title.ConvertPDBString(S);
        if ((RC!=Error_WrongSection) && ignoreNonCoorPDBErrors)
          RC = Error_NoError;
        if (RC)  break;
      }
      fend = f.FileEnd();
      if (!fend)  {
        ReadPDBLine ( f,S,sizeof(S) );
        lcount++;
      }
    } while (!fend);

    title.GetResolution(); // only to fetch resolution from remarks

    if (RC!=Error_WrongSection)  return RC;

    ignoreRemarks = (Flags & MMDBF_IgnoreRemarks)!=0;

    // read primary structure section
    SwitchModel ( 1 );
    if (!crModel)  return Error_GeneralError1;
    do {
      if (!strncmp(S,"FTNOTE",6))  {
        contString = new ContString(S);
        Footnote.AddData ( contString );
      } else  {
        RC = crModel->ConvertPDBString ( S );
        if ((RC!=Error_WrongSection) && ignoreNonCoorPDBErrors)
          RC = Error_NoError;
        if (RC)  break;
      }
      fend = f.FileEnd();
      if (!fend)  {
        ReadPDBLine ( f,S,sizeof(S) );
        title.TrimInput ( S );
        lcount++;
      }
    } while (!fend);

    if (RC!=Error_WrongSection)  return RC;

    // temporary solution: the rest of file is stored
    // in the form of strings
    while (!f.FileEnd()          &&
           strncmp(S,"CRYST" ,5) &&
           strncmp(S,"ORIGX" ,5) &&
           strncmp(S,"SCALE" ,5) &&
           strncmp(S,"MTRIX" ,5) &&
           strncmp(S,"TVECT" ,5) &&
           strncmp(S,"MODEL ",6) &&
           strncmp(S,"ATOM  ",6) &&
           strncmp(S,"SIGATM",6) &&
           strncmp(S,"ANISOU",6) &&
           strncmp(S,"SIGUIJ",6) &&
           strncmp(S,"TER   ",6) &&
           strncmp(S,"HETATM",6) &&
           strncmp(S,"ENDMDL",6))  {
      if (!strncmp(S,"LINK  ",6))
        crModel->ConvertPDBString ( S );
      else if (!strncmp(S,"LINKR ",6))
        crModel->ConvertPDBString ( S );
      else if (!strncmp(S,"CISPEP",6)) {
        GetInteger ( modNum,&(S[43]),3 );
        if (modNum<=0)  modNum = 1;
        if (modNum!=1)  SwitchModel ( modNum );
        crModel->ConvertPDBString ( S );
        if (modNum!=1)  SwitchModel ( 1 );
      } else  {
        contString = new ContString(S);
        SA.AddData ( contString );
      }
      ReadPDBLine ( f,S,sizeof(S) );
      title.TrimInput ( S );
      lcount++;
    }

    // read crystallographic information section
    do {
      RC = cryst.ConvertPDBString ( S );
      if ((RC!=Error_WrongSection) && ignoreNonCoorPDBErrors)
        RC = Error_NoError;
      if (RC)  break;
      fend = f.FileEnd();
      if (!fend)  {
        ReadPDBLine ( f,S,sizeof(S) );
        title.TrimInput ( S );
        lcount++;
      }
    } while (!fend);

    if (!RC)  {
      RC = cryst.ConvertPDBString ( S );
      if ((RC!=Error_WrongSection) && ignoreNonCoorPDBErrors)
        RC = Error_WrongSection;
    }

    cryst.CalcCoordTransforms();
    if (Flags & MMDBF_SimRWBROOK)
      cryst.RWBROOKReadPrintout();

    if (RC!=Error_WrongSection)  return RC;

    // temporary solution: the rest of file is stored
    // in the form of strings
    while (!f.FileEnd()          &&
           strncmp(S,"MODEL ",6) &&
           strncmp(S,"ATOM  ",6) &&
           strncmp(S,"SIGATM",6) &&
           strncmp(S,"ANISOU",6) &&
           strncmp(S,"SIGUIJ",6) &&
           strncmp(S,"TER   ",6) &&
           strncmp(S,"HETATM",6) &&
           strncmp(S,"ENDMDL",6))  {
      contString = new ContString(S);
      SB.AddData ( contString );
      ReadPDBLine ( f,S,sizeof(S) );
      title.TrimInput ( S );
      lcount++;
    }

    if (Flags & MMDBF_NoCoordRead)  return Error_NoError;

    // read coordinate section
    RC = Error_NoError;
    do {
      RC = ReadPDBAtom ( S );
      if (RC)  break;
      fend = f.FileEnd();
      if (!fend)  {
        ReadPDBLine ( f,S,sizeof(S) );
        title.TrimInput ( S );
        lcount++;
      }
    } while (!fend);
  //  if (!RC)
  //    RC = ReadPDBAtom(S);
  //  commented on 28.05.2004, it appears that "CHAIN_ORDER" should not
  //  be enforced here
  //  cleanKey = PDBCLEAN_ATNAME | PDBCLEAN_CHAIN_ORDER;
    cleanKey = 0x00000000;
    if (Flags & MMDBF_EnforceAtomNames)
      cleanKey = PDBCLEAN_ATNAME;
    if (Flags & MMDBF_AutoSerials)
      cleanKey |= PDBCLEAN_SERIAL;

    if (cleanKey)
      PDBCleanup ( cleanKey );

    if ((!f.FileEnd()) && (RC!=Error_WrongSection))  return RC;

    // temporary solution: the rest of file is stored
    // in the form of strings
    while (!f.FileEnd())  {
      if (strncmp(S,"END   ",6))  {  // END is added automatically
        contString = new ContString(S);
        SC.AddData ( contString );
      }
      ReadPDBLine ( f,S,sizeof(S) );
      title.TrimInput ( S );
      lcount++;
    }
    lcount--;  // last line was not read

    return Error_NoError;

  }


  ERROR_CODE Root::ReadCIFASCII1 ( cpstr CIFLFName, io::GZ_MODE gzipMode )  {
  pstr FName;
    FName = getenv ( CIFLFName );
    if (FName)  return ReadCIFASCII ( FName,gzipMode );
          else  return Error_NoLogicalName;
  }

  ERROR_CODE Root::ReadCIFASCII ( cpstr CIFFileName, io::GZ_MODE gzipMode )  {
  io::File   f;
  ERROR_CODE rc;

    //  open the file as ASCII for reading
    //  opening it in pseudo-binary mode helps reading various
    //  line terminators for files coming from different platforms
    f.assign ( CIFFileName,false,false,gzipMode );

    if (f.reset(true)) {
      rc = ReadCIFASCII ( f );
      f.shut();
    } else
      rc = Error_CantOpenFile;

    return rc;

  }

  ERROR_CODE Root::ReadCIFASCII ( io::RFile f )  {
  int        W;
  ERROR_CODE RC;

    //  remove previous data
    ResetManager  ();
    FreeFileMemory();
    FType = MMDB_FILE_CIF;

    SetFlag ( 0 );

    CIFErrorLocation[0] = char(0);  // CIF reading phase

    lcount = 0;  // line counter
    S[0]   = char(0);

    if (f.FileEnd())
      return Error_EmptyFile;

    if (!CIF)  CIF = new mmcif::Data();
    CIF->SetStopOnWarning  ( true );
    CIF->SetPrintWarnings  ( (Flags & MMDBF_PrintCIFWarnings)!=0 );
    W = CIF->ReadMMCIFData ( f,S,lcount );

    if (W)  {
      if (W == mmcif::CIFRC_NoDataLine)      return Error_NotACIFFile;
      if (W & mmcif::CIFW_UnrecognizedItems) return Error_UnrecognCIFItems;
      if (W & mmcif::CIFW_MissingField)      return Error_MissingCIFField;
      if (W & mmcif::CIFW_EmptyLoop)         return Error_EmptyCIFLoop;
      if (W & mmcif::CIFW_UnexpectedEOF)     return Error_UnexpEndOfCIF;
      if (W & mmcif::CIFW_LoopFieldMissing)  return Error_MissgCIFLoopField;
      if (W & mmcif::CIFW_NotAStructure)     return Error_NotACIFStructure;
      if (W & mmcif::CIFW_NotALoop)          return Error_NotACIFLoop;
      return Error_Unknown;
    }

    RC = ReadFromCIF ( CIF );
    if (CIF)  {
      delete CIF;
      CIF = NULL;
    }

    return RC;

  }


  ERROR_CODE Root::ReadFromCIF ( mmcif::PData CIFD )  {
  mmcif::PLoop  Loop1,Loop2;
  pstr          F,FC;
  word          cleanKey;
  int           i,l,j,n,retc;
  ERROR_CODE    RC;

    RC = title.GetCIF ( CIFD );

    if (RC!=Error_NoError)  {
      CIFD->Optimize();
      return RC;
    }

    SwitchModel ( 1 );
    if (!crModel)  return Error_GeneralError1;
    RC = crModel->GetCIF ( CIFD );
    if (RC!=Error_NoError)  {
      CIFD->Optimize();
      return RC;
    }

    RC = cryst.GetCIF ( CIFD );
    if (RC!=Error_NoError)  {
      CIFD->Optimize();
      return RC;
    }
    cryst.CalcCoordTransforms();
    if (Flags & MMDBF_SimRWBROOK)
      cryst.RWBROOKReadPrintout();

    RC = ReadCIFAtom ( CIFD );

    Loop1 = CIFD->GetLoop ( CIFCAT_ENTITY      );
    Loop2 = CIFD->GetLoop ( CIFCAT_STRUCT_ASYM );
    if (Loop1 && Loop2)  {
      // make 'Het' atoms
      l = Loop1->GetLoopLength();
      n = Loop2->GetLoopLength();
      for (i=0;i<l;i++)  {
        F = Loop1->GetString ( CIFTAG_TYPE,i,retc );
        if (F && (!retc))  {
          if (!strcasecmp(F,"non-polymer"))  {
            F = Loop1->GetString ( CIFTAG_ID,i,retc );
            if (F && (!retc))
              for (j=0;j<n;j++)  {
                FC = Loop2->GetString ( CIFTAG_ENTITY_ID,j,retc );
                if (FC && (!retc))  {
                  if (!strcasecmp(FC,F))  {
                    FC = Loop2->GetString ( CIFTAG_ID,j,retc );
                    if (FC && (!retc))
                      MakeHetAtoms ( FC,true );
                  }
                }
              }
          }
        }
      }
    }

    if (RC==Error_NoError)  {
      //  deleting these CIF loops here is a temporary solution
      // taken in order to avoid mess at rewriting the CIF file.
      CIFD->DeleteLoop ( CIFCAT_ATOM_SITE           );
      CIFD->DeleteLoop ( CIFCAT_ATOM_SITE_ANISOTROP );
      CIFD->Optimize   ();
    }

    cleanKey = 0x00000000;
    if (Flags & MMDBF_EnforceAtomNames)
      cleanKey = PDBCLEAN_ATNAME;
    if (Flags & MMDBF_AutoSerials)
      cleanKey |= PDBCLEAN_SERIAL;
    if (cleanKey)
      PDBCleanup ( cleanKey );

    return RC;

  }

  ERROR_CODE Root::ReadCoorFile1 ( cpstr LFName, io::GZ_MODE gzipMode )  {
  pstr FName;
    FName = getenv ( LFName );
    if (FName)  return ReadCoorFile ( FName,gzipMode );
          else  return Error_NoLogicalName;
  }

  ERROR_CODE Root::ReadCoorFile ( cpstr CFName, io::GZ_MODE gzipMode )  {
  // auto format recognition
  int  kin;
  bool IBL;

    kin = isMMDBBIN ( CFName,gzipMode );
    if (kin==Error_EmptyFile)
                return Error_EmptyFile;
    if (kin<0)  return Error_CantOpenFile;

    if (kin==0) return  ReadMMDBF ( CFName,gzipMode );

    IBL = ((Flags & MMDBF_IgnoreBlankLines)!=0);
    if (isPDB(CFName,gzipMode,IBL)==0)
      return ReadPDBASCII ( CFName,gzipMode );

    if (mmcif::isCIF(CFName,gzipMode)==0)
      return ReadCIFASCII ( CFName,gzipMode );

    return Error_ForeignFile;

  }


  ERROR_CODE Root::ReadCoorFile ( io::RFile f )  {
  // auto format recognition
  int  kin;
  bool IBL;

    kin = isMMDBBIN ( f );
    f.reset ( true );
    if (kin==Error_EmptyFile)
                return Error_EmptyFile;
    if (kin<0)  return Error_CantOpenFile;

    if (kin==0) return  ReadMMDBF ( f );

    IBL = ((Flags & MMDBF_IgnoreBlankLines)!=0);
    kin = isPDB ( f,IBL );
    f.reset ( true );
    if (kin==0)
      return ReadPDBASCII ( f );

    kin = mmcif::isCIF ( f );
    f.reset ( true );
    if (kin==0)
      return ReadCIFASCII ( f );

    return Error_ForeignFile;

  }


  word  Root::PDBCleanup ( word CleanKey )  {
  //  cleans coordinate part to comply with PDB standards:
  //
  //    CleanKey          Action
  //  PDBCLEAN_ATNAME  pads atom names with spaces to form 4-symbol names
  //  PDBCLEAN_TER     inserts TER cards in the end of each chain
  //  PDBCLEAN_CHAIN   generates 1-character chain ids instead of
  //                   those many-character
  //  PDBCLEAN_CHAIN_STRONG generates 1-character chain ids starting
  //                   from 'A' on for all ids, including single-char
  //  PDBCLEAN_ALTCODE generates 1-character alternative codes instead
  //                   of those many-character
  //  PDBCLEAN_ALTCODE_STRONG generates 1-character alternative codes
  //                   from 'A' on for all codes, including
  //                   single-character ones
  //  PDBCLEAN_SERIAL  puts serial numbers in due order
  //  PDBCLEAN_INDEX   reorders the internal index of atoms such that
  //                   it follows the actual order of atoms in
  //                   the object hierarchy
  //  PDBCLEAN_SEQNUM  renumbers all residues so that they go
  //                   incrementally-by-one without insertion codes
  //  PDBCLEAN_CHAIN_ORDER puts chains in order of atom's serial numbers
  //  PDBCLEAN_SOLVENT moves solvent chains at the end of each model
  //  PDBCLEAN_ELEMENT calculates PDB element names where they are not
  //                   found in the chemical element table
  //  PDBCLEAN_ELEMENT_STRONG  calculates all chemical element names
  //
  //  Return codes (as bits):
  //  0                Ok
  //  PDBCLEAN_CHAIN   too many chains for assigning them 1-letter codes
  //  PDBCLEAN_ATNAME  element names were not available
  //  PDBCLEAN_ALTCODE too many alternative codes encountered.
  //
  PPAtom    atom1;
  PPChain   Chain1,Chain2;
  PModel    crModel0;
  PChain    crChain0;
  PResidue  crRes0;
  PPResidue resTable;
  PAtom     A;
  AltLoc  * altLoc;
  ChainID * chain_ID;
//  ChainID   chainID;
//  pstr      chID;
  char      aLoc [257];
  char      chnID[257];
  word      RC;
  int       i,j,k,nal,nch,nr, nch1,nch2;
//  int       modN,modl;
  char      c;
  bool      Done,Solvent;
//  bool      NewChain;

    RC = 0;
    if (nAtoms<=0)  return RC;

    if (CleanKey & PDBCLEAN_ATNAME)
      for (i=0;i<nAtoms;i++)
        if (atom[i])
          if (!atom[i]->MakePDBAtomName())
            RC |= PDBCLEAN_ATNAME;

    k = -1;

    if (CleanKey & PDBCLEAN_TER)  {
      for (i=0;i<nModels;i++)
        if (model[i])  {
          model[i]->GetChainTable ( Chain1,nch1 );
          for (j=0;j<nch1;j++)
            if (Chain1[j])  {
              Chain1[j]->TrimResidueTable();
              Chain1[j]->GetResidueTable ( resTable,nr );
              if (nr>0)  {
                nr--;
                Done = false;
                while ((nr>=0) && (!Done))  {
                  Done = resTable[nr]->isAminoacid();
                  if (!Done)
                    nr--;
                }
                if (nr>=0)  {
                  resTable[nr]->TrimAtomTable();
                  resTable[nr]->GetAtomTable ( atom1,nal );
                  if (nal>0)  {
                    nal--;
                    if (!atom1[nal]->isTer())  {
                      A = newAtom();
                      A->Copy ( atom1[nal] );
                      A->MakeTer();
                      resTable[nr]->AddAtom ( A );
                    }
                  }
                }
              }
            }
        }
    
    /*
      modN     = -1;
      crModel0 = crModel;
      for (i=0;i<nAtoms;i++)
        if (atom[i])  {
          modl = atom[i]->GetModelNum();
          chID = atom[i]->GetChainID ();
          if (modN<0)  {
            modN = modl;
            SwitchModel ( modN );
            if (chID)  strcpy ( chainID,chID );
                 else  chainID[0] = char(0);
          } else  {
            if (modN!=modl)   NewChain = true;
            else if (chID)    NewChain = strcmp(chID,chainID)!=0;
                        else  NewChain = chainID[0]!=char(0);
            if (NewChain)  {
              if (k>=0)  {
                if ((!atom[k]->Ter) && (!atom[k]->Het))  {
                  // insert 'Ter' before atom in position 'i'
                  PutAtom ( -(i+1),atom[k]->serNum+1,pstr("TER"),
                            atom[k]->GetResName(),atom[k]->GetChainID(),
                            atom[k]->GetSeqNum (),atom[k]->GetInsCode(),
                            pstr(" "),pstr(" "),pstr(" ") );
                  atom[i]->MakeTer();
                }
              }
              modN = modl;
              SwitchModel ( modN );
              if (chID)  strcpy ( chainID,chID );
                   else  chainID[0] = char(0);
            }
          }
          k = i;
        }

      if (k>=0)  {
        if ((!atom[k]->Ter) && (!atom[k]->Het))  {  // add last TER
          i = nAtoms;
          SwitchModel ( atom[k]->GetModelNum() );
          PutAtom ( 0,nAtoms+1,pstr("TER"),atom[k]->GetResName(),
                    atom[k]->GetChainID(),atom[k]->GetSeqNum(),
                    atom[k]->GetInsCode(),pstr(" "),pstr(" "),
                    pstr(" ") );
          atom[i]->MakeTer();
        }
      }

      crModel = crModel0;
      
      */

    }


    if (CleanKey & (PDBCLEAN_CHAIN | PDBCLEAN_CHAIN_STRONG))  {
      chain_ID = new ChainID[256];
      for (i=0;i<nModels;i++)
        if (model[i])  {
          for (j=0;j<256;j++)  {
            strcpy ( chain_ID[j]," " );
            chnID[j] = char(0);
          }
          chnID[256] = char(0);
          nch = 0;
          for (j=0;j<model[i]->nChains;j++)  {
            crChain0 = model[i]->chain[j];
            if (crChain0)  {
              if (!crChain0->chainID[0])
                strcpy ( crChain0->chainID," " );
              k = 0;
              while ((k<nch) && (strcmp(chain_ID[k],crChain0->chainID)))
                k++;
              if (k>=nch)  {
                if (nch>=255)  RC |= PDBCLEAN_CHAIN;
                else  {
                  strcpy ( chain_ID[nch],crChain0->chainID );
                  if (!chain_ID[nch][1])
                    chnID[nch] = chain_ID[nch][0];
                  nch++;
                }
              }
            }
          }
          c = 'A';
          if (CleanKey & PDBCLEAN_CHAIN_STRONG)  {
            // rename all chains through from A to Z
            for (k=0;k<nch;k++)  {
              chnID[k] = c;
              c = char(int(c)+1);
            }
          } else  {
            // rename only multi-character chain IDs
            for (j=0;(j<nch) && (k<256);j++)  {
              k = 0;
              do  {
                while ((k<nch) && (chnID[k]!=c))  k++;
                if (k<nch)  c = char(int(c)+1);
              } while (k<nch);
              k = 0;
              while ((k<256) && (chnID[k]))  k++;
              if (k<256)  {
                chnID[k] = c;
                c = char(int(c)+1);
              }
            }
          }
          // assign new chain IDs
          for (j=0;j<model[i]->nChains;j++)  {
            crChain0 = model[i]->chain[j];
            if (crChain0)  {
              k = 0;
              while ((k<nch) && (strcmp(chain_ID[k],crChain0->chainID)))
                k++;
              strcpy ( crChain0->prevChainID,crChain0->chainID );
              crChain0->chainID[0] = chnID[k];
              crChain0->chainID[1] = char(0);
            }
          }
        }
      delete[] chain_ID;
    }


    if (CleanKey & (PDBCLEAN_ALTCODE | PDBCLEAN_ALTCODE_STRONG))  {
      altLoc = new AltLoc[256];
      for (i=0;i<256;i++)  {
        strcpy ( altLoc[i]," " );
        aLoc[i] = char(0);
      }
      aLoc[0]   = ' ';
      aLoc[256] = char(0);
      nal = 1;
      for (i=0;i<nAtoms;i++)
        if (atom[i])  {
          if (!atom[i]->altLoc[0])  strcpy ( atom[i]->altLoc," " );
          else  {
            k = 0;
            while ((k<nal) && (strcmp(altLoc[k],atom[i]->altLoc)))  k++;
            if (k>=nal)  {
              if (nal>=255)  RC |= PDBCLEAN_ALTCODE;
              else  {
                strcpy ( altLoc[nal],atom[i]->altLoc );
                if (!altLoc[nal][1])  aLoc[nal] = altLoc[nal][0];
                nal++;
              }
            }
          }
        }
      c = 'A';
      if (CleanKey & PDBCLEAN_ALTCODE_STRONG)
        for (i=1;i<nal;i++)  {
          aLoc[i] = c;
          c = char(int(c)+1);
        }
      else
        for (i=1;(i<nal) && (k<256);i++)  {
          k = 0;
          do  {
            while ((k<nal) && (aLoc[k]!=c))  k++;
            if (k<nal)  c = char(int(c)+1);
          } while (k<nal);
          k = 0;
          while ((k<256) && (aLoc[k]))  k++;
          if (k<256)  {
            aLoc[k] = c;
            c = char(int(c)+1);
          }
        }
      for (i=0;i<nAtoms;i++)
        if (atom[i])  {
          k = 0;
          while ((k<nal) && (strcmp(altLoc[k],atom[i]->altLoc)))  k++;
          atom[i]->altLoc[0] = aLoc[k];
          atom[i]->altLoc[1] = char(0);
        }
      delete[] altLoc;
    }


    if (CleanKey & PDBCLEAN_SEQNUM)
      for (i=0;i<nModels;i++)  {
        crModel0 = model[i];
        if (crModel0)
          for (j=0;j<crModel0->nChains;j++)  {
            crChain0 = crModel0->chain[j];
            if (crChain0)  {
              nr = 0;
              for (k=0;k<crChain0->nResidues;k++)  {
                crRes0 = crChain0->residue[k];
                if (crRes0)  {
                  nr++;
                  crRes0->seqNum     = nr;
                  crRes0->insCode[0] = char(0);
                }
              }
            }
          }
      }

    if (CleanKey & PDBCLEAN_SOLVENT)  {
      atom1 = new PAtom[nAtoms];
      k = 1;
      for (i=0;i<nModels;i++)
        if (model[i])  {
          if (model[i]->nChains>k)  k = model[i]->nChains;
        }
      Chain1 = new PChain[k];
      Chain2 = new PChain[k];
      k = 0;
      for (i=0;i<nModels;i++)  {
        crModel0 = model[i];
        if (crModel0)  {
          nch1 = 0;
          nch2 = 0;
          for (nch=0;nch<crModel0->nChains;nch++)  {
            crChain0 = crModel0->chain[nch];
            if (crChain0)  {
              Solvent = false;
              for (nr=0;(nr<crChain0->nResidues) && (!Solvent);nr++)  {
                crRes0 = crChain0->residue[nr];
                if (crRes0)
                  for (j=0;(j<nSolventNames) && (!Solvent);j++)
                    Solvent = !strcmp ( StdSolventName[j],crRes0->name );
              }
              if (Solvent)  Chain2[nch2++] = crChain0;
                      else  Chain1[nch1++] = crChain0;
            }
          }
          for (nch=0;nch<nch1;nch++)  {
            crChain0 = Chain1[nch];
            for (nr=0;nr<crChain0->nResidues;nr++)  {
              crRes0 = crChain0->residue[nr];
              if (crRes0)
                for (j=0;j<crRes0->nAtoms;j++)
                  if (crRes0->atom[j])  {
                    atom1[k] = crRes0->atom[j];
                    atom1[k]->index = k+1;
                    k++;
                  }
            }
            crModel0->chain[nch] = Chain1[nch];
          }
          for (nch=0;nch<nch2;nch++)  {
            crChain0 = Chain2[nch];
            for (nr=0;nr<crChain0->nResidues;nr++)  {
              crRes0 = crChain0->residue[nr];
              if (crRes0)
                for (j=0;j<crRes0->nAtoms;j++)
                  if (crRes0->atom[j])  {
                    atom1[k] = crRes0->atom[j];
                    atom1[k]->index = k+1;
                    k++;
                  }
            }
            crModel0->chain[nch1++] = Chain2[nch];
          }
          crModel0->nChains = nch1;
        }
      }
      delete[] Chain1;
      delete[] Chain2;
      if (atom)  delete[] atom;
      atom   = atom1;
      atmLen = nAtoms;
      nAtoms = k;
    }

    if (CleanKey & (PDBCLEAN_CHAIN_ORDER | PDBCLEAN_CHAIN_ORDER_IX))  {
      for (i=0;i<nModels;i++)  {
        crModel0 = model[i];
        if (crModel0)  {
          k = 0;
          for (j=0;j<crModel0->nChains;j++)  {
            crChain0 = crModel0->chain[j];
            if (crChain0)  {
              crChain0->nWeights = 0;
              crChain0->Weight   = 0.0;
              if (k<j)  {
                crModel0->chain[k] = crModel0->chain[j];
                crModel0->chain[j] = NULL;
              }
              k++;
            }
          }
          crModel0->nChains = k;
        }
      }
      if (CleanKey & PDBCLEAN_CHAIN_ORDER)  {
        for (i=0;i<nAtoms;i++)
          if (atom[i])  {
            crChain0 = atom[i]->GetChain();
            crChain0->nWeights++;
            crChain0->Weight += atom[i]->serNum;
          }
      } else  {
        for (i=0;i<nAtoms;i++)
          if (atom[i])  {
            crChain0 = atom[i]->GetChain();
            crChain0->nWeights++;
            crChain0->Weight += atom[i]->GetIndex();
          }
      }
      for (i=0;i<nModels;i++)  {
        crModel0 = model[i];
        if (crModel0)  {
          for (j=0;j<crModel0->nChains;j++)  {
            crChain0 = crModel0->chain[j];
            if (crChain0->nWeights)
              crChain0->Weight /= crChain0->nWeights;
          }
          //  bubble sorting
          do {
            Done = true;
            for (j=1;j<crModel0->nChains;j++)
              if (crModel0->chain[j-1]->Weight >
                  crModel0->chain[j]->Weight)  {
                crChain0             = crModel0->chain[j-1];
                crModel0->chain[j-1] = crModel0->chain[j];
                crModel0->chain[j]   = crChain0;
                Done = false;
              }
          } while (!Done);
        }
      }
    }

    if (CleanKey & PDBCLEAN_INDEX)  {
      k = 0;
      for (i=0;i<nModels;i++)  {
        crModel0 = model[i];
        if (crModel0)  {
          for (nch=0;nch<crModel0->nChains;nch++)  {
            crChain0 = crModel0->chain[nch];
            if (crChain0)  {
              for (nr=0;nr<crChain0->nResidues;nr++)  {
                crRes0 = crChain0->residue[nr];
                if (crRes0)  {
                  for (j=0;j<crRes0->nAtoms;j++)  {
                    A = crRes0->atom[j];
                    if (A)  {
                      atom[A->index-1] = atom[k];
                      if (atom[k])
                        atom[k]->index = A->index;
                      atom[k] = A;
                      k++;
                      A->index = k;
                    }
                  }
                }
              }
            }
          }
        }
      }
      nAtoms = k;
    }

    if (CleanKey & PDBCLEAN_SERIAL)  {
      k = 0;
      for (i=0;i<nAtoms;i++)
        if (atom[i])  {
          if (k<i)  {
            atom[k] = atom[i];
            atom[i] = NULL;
          }
          atom[k]->index  = k+1;
          atom[k]->serNum = atom[k]->index;
          k++;
        }
      nAtoms = k;
    }

    if (CleanKey & PDBCLEAN_ELEMENT)  {
      for (i=0;i<nAtoms;i++)
        if (atom[i] && (!atom[i]->Ter))  {
          if (getElementNo(atom[i]->element)==ELEMENT_UNKNOWN)  {
            strcpy ( atom[i]->element,"  " );
            //atom[i]->MakePDBAtomName();
            atom[i]->RestoreElementName();
          }
        }
    }

    if (CleanKey & PDBCLEAN_ELEMENT_STRONG)  {
      for (i=0;i<nAtoms;i++)
        if (atom[i] && (!atom[i]->Ter))  {
          strcpy ( atom[i]->element,"  " );
          //atom[i]->MakePDBAtomName();
          atom[i]->RestoreElementName();
        }
    }

    return RC;

  }

  void  Root::MakeHetAtoms ( cpstr chainID, bool Make )  {
  //  Makes all atoms in chain 'chainID', in all models, as 'Het' atoms
  //  if Make is set true, and makes them 'ordinary' atoms otherwise.
  //  'Ter' is automatically removed when converting to 'Het' atoms,
  //  and is automatically added when converting to 'ordinary' atoms.
  int       i,j,k,l,n;
  PModel   crModel0;
  PChain   crChain0;
  PResidue crRes0;
    crModel0 = crModel;
    for (i=0;i<nModels;i++)
      if (model[i])
        for (j=0;j<model[i]->nChains;j++)  {
          crChain0 = model[i]->chain[j];
          if (crChain0)  {
            if (!strcmp(crChain0->chainID,chainID))  {
              n = 0;
              for (k=0;k<crChain0->nResidues;k++)  {
                crRes0 = crChain0->residue[k];
                if (crRes0)
                  for (l=0;l<crRes0->nAtoms;l++)
                    if (crRes0->atom[l])  {
                      crRes0->atom[l]->Het = Make;
                      n = crRes0->atom[l]->index;
                    }
              }
              if (n>0)  {
                n--;
                if (atom[n]->Het && atom[n]->Ter)  RemoveAtom ( n+1 );
                else if ((!atom[n]->Het) && (!atom[n]->Ter))  {
                  SwitchModel ( model[i]->GetSerNum() );
                  if (n<nAtoms-1)
                    PutAtom ( -(n+2),atom[n]->serNum+1,pstr("TER"),
                            atom[n]->GetResName(),atom[n]->GetChainID(),
                            atom[n]->GetSeqNum (),atom[n]->GetInsCode(),
                            pstr(" "),pstr(" "),pstr(" ") );
                  else
                    PutAtom ( 0,nAtoms+1,pstr("TER"),
                            atom[n]->GetResName(),atom[n]->GetChainID(),
                            atom[n]->GetSeqNum (),atom[n]->GetInsCode(),
                            pstr(" "),pstr(" "),pstr(" ") );
                  atom[n+1]->MakeTer();
                }
              }
            }
          }
        }
    crModel = crModel0;
  }


  void Root::RemoveAtom ( int index )  {
  //    Removes atom at the specified index in the Atom array.
  // This index is always accessible as atom[index]->index.
  // If this leaves a residue empty, the residue is removed.
  // If this leaves an empty chain, the chain is removed as well;
  // the same happens to the model.
  PResidue crRes0;
  PChain   crChain0;
  PModel   crModel0;
  int      i,j;

    if ((index>0) && (index<=nAtoms))  {
      if (atom[index-1])  {
        crRes0 = atom[index-1]->residue;
        if (crRes0)  {
          if (crRes0->_ExcludeAtom(index))  {
            // the residue appears empty after the exclusion
            if (crRes)  {
              if ((crRes->seqNum==crRes0->seqNum) &&
                  (!strcmp(crRes->insCode,crRes0->insCode)))
                crRes = NULL;
            }
            crChain0 = crRes0->chain;
            if (crChain0)  {
              if (crChain0->_ExcludeResidue(crRes0->name,crRes0->seqNum,
                                            crRes0->insCode))  {
                // the chain appears empty after the exclusion
                if (crChain)  {
                  if (!strcmp(crChain->chainID,crChain0->chainID))
                    crChain = NULL;
                }
                crModel0 = PModel(crChain0->model);
                if (crModel0)  {
                  if (crModel0->_ExcludeChain(crChain0->chainID))  {
                    // the model appears ampty after the exclusion
                    if (crModel)  {
                      if (crModel->serNum==crModel0->serNum)
                        crModel = NULL;
                    }
                    i = crModel0->serNum-1;
                    delete model[i];
                    model[i] = NULL;
                  }
                }
                delete crChain0;  // it is already excluded from the hierarchy!
              }
            }
            delete crRes0;  // it is already excluded from the hierarchy!
          }
        }
        delete atom[index-1];  // it is already excluded from the hierarchy!
        atom[index-1] = NULL;
        // now rearrange and re-index atoms.
        j = 0;
        for (i=0;i<nAtoms;i++)
          if (atom[i])  {
            if (j<i)  {
              atom[j] = atom[i];
              atom[i] = NULL;
            }
            atom[j]->index = j+1;
            j++;
          }
        nAtoms = j;
      }
    }
  }


  int  Root::_ExcludeModel ( int serNum )  {
  //   _ExcludeModel(..) excludes (but does not dispose!) a model
  // from the file. Returns 1 if the file gets empty and 0 otherwise.
  int  i,k;

    if (!Exclude)  return 0;

    if ((0<serNum) && (serNum<=nModels))
      model[serNum-1] = NULL;

    k = 0;
    for (i=0;i<nModels;i++)
      if (model[i])  {
        if (k<i)  {
          model[k] = model[i];
          model[i] = NULL;
        }
        model[k]->serNum = k+1;
        k++;
      }

    nModels = k;

    if (nModels<=0)  return 1;
               else  return 0;

  }


  int  Root::FinishStructEdit()  {
  // Makes a new atom index after insertion or deletion of atoms.
  // This function may change atoms' positions in the index and
  // correspondingly the Atom::index field.
  PResidue res;
  PChain   chain;
  PModel   Model1;
  PPAtom   Atom1;
  int       i,j,k,l,n,index,nAtoms1;

    //  calculate new number of atoms
    nAtoms1 = 0;
    for (i=0;i<nModels;i++)  {
      Model1 = model[i];
      if (Model1)  {
        for (j=0;j<Model1->nChains;j++)  {
          chain = Model1->chain[j];
          if (chain)  {
            for (k=0;k<chain->nResidues;k++)  {
              res = chain->residue[k];
              if (res)  {
                res->TrimAtomTable();
                nAtoms1 += res->nAtoms;
              }
            }
            chain->TrimResidueTable();
          }
        }
        Model1->TrimChainTable();
      }
    }
    TrimModelTable();

    // compile a new index and null the old one

    if (nAtoms1>0)  Atom1 = new PAtom[nAtoms1];
              else  Atom1 = NULL;

    n = 0;
    for (i=0;i<nModels;i++)  {
      Model1 = model[i];
      for (j=0;j<Model1->nChains;j++)  {
        chain = Model1->chain[j];
        for (k=0;k<chain->nResidues;k++)  {
          res = chain->residue[k];
          for (l=0;l<res->nAtoms;l++)  {
            Atom1[n] = res->atom[l];
            index    = Atom1[n]->index;
            if ((index>0) && (index<=atmLen))
              atom[index-1] = NULL;
            Atom1[n]->index = n+1;
            n++;
          }
        }
      }
    }

  //  if (n!=nAtoms1)  {
  //    printf ( " **** PROGRAM ERROR IN Root::FinishStructEdit\n" );
  //    exit ( 1 );
  //  }


    // check if there are dead atoms in the old index
    for (i=0;i<atmLen;i++)
      if (atom[i])  delete atom[i];

    // dispose old index and replace it with the new one
    if (atom)  delete[] atom;

    atom   = Atom1;
    atmLen = n;
    nAtoms = n;

    if (n==nAtoms1)  return 0;  // Ok
               else  return 1;  // not Ok; should never happen

  }

  void Root::TrimModelTable()  {
  int i,j;
    j = 0;
    for (i=0;i<nModels;i++)
      if (model[i])  {
        if (j<i)  {
          model[j] = model[i];
          model[i] = NULL;
        }
        model[j]->serNum = j+1;
        j++;
      }
    nModels = j;
  }


  int  Root::GenerateNCSMates()  {
  //
  //   Generates NCS mates according to NCS matrices given
  // in cryst. This will result in generating many-character
  // chain names, composed as 'x_n' where 'x' is the original
  // name and 'n' is a unique number, which will coincide with
  // the symmetry operation (order) number. Another side
  // effect will be a disorder in atoms' serial numbers.
  //   The hierarchy should therefore be cleaned after
  // generating the NCS mates. An appropriate way to do that
  // is to issue the following call:
  //
  //   PDBCleanup ( PDBCLEAN_TER | PDBCLEAN_ALTCODE_STRONG |
  //                PDBCLEAN_CHAIN_STRONG | PDBCLEAN_SERIAL );
  //
  PPChain chainTable,chain;
  PChain  chn;
  mat44    ncs_m;
  ChainID  chainID;
  int      i,j,k,nNCSOps,nChains,iGiven;

    nNCSOps = cryst.GetNumberOfNCSMatrices();
    if (nNCSOps<=0)  return 1;

    for (i=0;i<nModels;i++)
      if (model[i])  {
        model[i]->GetChainTable ( chainTable,nChains );
        if (nChains>0)  {
          chain = new PChain[nChains];
          for (j=0;j<nChains;j++)
            chain[j] = chainTable[j];
          for (j=0;j<nChains;j++)
            if (chain[j])  {
              for (k=0;k<nNCSOps;k++)
                if (cryst.GetNCSMatrix(k,ncs_m,iGiven))  {
                  if (!iGiven)  {
                    chn = newChain();
                    chn->Copy ( chain[j] );
                    sprintf ( chainID,"%s_%i",
                              chain[j]->GetChainID(),k+1 );
                    chn->SetChainID     ( chainID );
                    chn->ApplyTransform ( ncs_m   );
                    model[i]->AddChain  ( chn     );
                  }
                }
            }
          delete[] chain;
        }
      }

    return 0;

  }


  void  Root::ApplyNCSTransform ( int NCSMatrixNo )  {
  mat33 t;
  vect3 v;
  int   i;
    if (!cryst.GetNCSMatrix(NCSMatrixNo,t,v))  return;
    for (i=0;i<nAtoms;i++)
      if (atom[i])  atom[i]->Transform ( t,v );
  }


  ERROR_CODE Root::PutPDBString ( cpstr PDBString )  {
  PContString contString;
  ERROR_CODE  RC;

    strcpy    ( S,PDBString );  // maintain the buffer!
    PadSpaces ( S,80 );
    lcount++;

    // belongs to title?
    RC = title.ConvertPDBString ( S );
    if (RC!=Error_WrongSection)  return RC;

    // belongs to primary structure section?
    SwitchModel ( 1 );
    RC = crModel->ConvertPDBString ( S );
    if (RC!=Error_WrongSection)  return RC;

    // belongs to the crystallographic information section?
    RC = cryst.ConvertPDBString ( S );
    if (RC!=Error_WrongSection)  {
  //    if (RC==0)  cryst.CalcCoordTransforms();
      return RC;
    }

    // belongs to the coordinate section?
    RC = ReadPDBAtom ( S );
    if (RC!=Error_WrongSection)  return RC;

    // temporary solution: the rest of file is stored
    // in the form of strings
    if ((S[0]) && (S[0]!=' ') && (strncmp(S,"END   ",6)))  {
      // END is added automatically
      contString = new ContString(S);
      SC.AddData ( contString );
    }

    return Error_NoError;

  }


  ERROR_CODE Root::AddPDBASCII1 ( cpstr PDBLFName, io::GZ_MODE gzipMode )  {
  pstr FName;
    FName = getenv ( PDBLFName );
    if (FName)  return AddPDBASCII ( FName,gzipMode );
          else  return Error_NoLogicalName;
  }

  ERROR_CODE Root::AddPDBASCII ( cpstr PDBFileName, io::GZ_MODE gzipMode ) {
  io::File   f;
  ERROR_CODE RC;
    //  open the file as ASCII for reading
    //  opening it in pseudo-binary mode helps reading various
    //  line terminators for files coming from different platforms
    f.assign ( PDBFileName,false,false,gzipMode );
    if (f.reset(true)) {
      lcount = 1;  // line counter
      RC     = Error_NoError;
      while ((!f.FileEnd()) && (!RC))  {
        ReadPDBLine ( f,S,sizeof(S) );
        RC = PutPDBString ( S );
      }
      f.shut();
    } else
      RC = Error_CantOpenFile;
    return RC;
  }


  void Root::GetInputBuffer ( pstr Line, int & count )  {
    if (FType==MMDB_FILE_PDB)  {  // PDB File
      strcpy ( Line,S );
      count = lcount;
    } else if (FType==MMDB_FILE_CIF)  {
      if (!CIFErrorLocation[0])  {  // CIF reading phase
        strcpy ( Line,S );
        count = lcount;
      } else  {
        strcpy ( Line,CIFErrorLocation );
        count = -1;  // CIF interpretation phase
      }
    } else {
      Line[0] = char(0);
      count = -2;
    }
  }

  bool Root::isCompactBinary()  {
    if (Flags & MMDBF_MakeCompactBinary)  return true;
    return false;
  }

  
  int  Root::CrystReady()  {
  //    Returns flags:
  // CRRDY_Complete       if crystallographic information is complete
  // CRRDY_NotPrecise     if cryst. inf-n is not precise
  // CRRDY_isTranslation  if cryst. inf-n contains translation
  // CRRDY_NoOrthCode     no orthogonalization code
  //    Fatal:
  // CRRDY_NoTransfMatrices  if transform. matrices were not calculated
  // CRRDY_Unchecked         if cryst. inf-n was not checked
  // CRRDY_Ambiguous         if cryst. inf-n is ambiguous
  // CRRDY_NoCell            if cryst. inf-n is unusable
  // CRRDY_NoSpaceGroup      if space group is not set
  int k;

    if (!(cryst.WhatIsSet & CSET_Transforms))
      return CRRDY_NoTransfMatrices;

    if ((cryst.WhatIsSet & CSET_CellParams)!=CSET_CellParams)
      return CRRDY_NoCell;

    if (!(cryst.WhatIsSet & CSET_SpaceGroup))
      return CRRDY_NoSpaceGroup;

    if (cryst.CellCheck & CCHK_Unchecked)
      return CRRDY_Unchecked;

    if (cryst.CellCheck & CCHK_Disagreement)
      return CRRDY_Ambiguous;

    k = 0x0000;
    if (cryst.CellCheck & CCHK_Error)        k |= CRRDY_NotPrecise;
    if (cryst.CellCheck & CCHK_Translations) k |= CRRDY_isTranslation;
    if (cryst.CellCheck & CCHK_NoOrthCode)   k |= CRRDY_NoOrthCode;

    return k;

  }


  bool Root::isCrystInfo()  {
    return (((cryst.WhatIsSet & CSET_CellParams)==CSET_CellParams) &&
             (cryst.WhatIsSet & CSET_SpaceGroup));
  }

  bool Root::isCellInfo()  {
    return ((cryst.WhatIsSet & CSET_CellParams)==CSET_CellParams);
  }

  bool Root::isSpaceGroup()  {
    return (cryst.WhatIsSet & CSET_SpaceGroup);
  }

  bool Root::isTransfMatrix()  {
    return cryst.areMatrices();
  }

  bool Root::isScaleMatrix()  {
    return ((cryst.WhatIsSet & CSET_ScaleMatrix)==CSET_ScaleMatrix);
  }

  bool Root::isNCSMatrix()  {
    return cryst.isNCSMatrix();
  }

  int  Root::AddNCSMatrix ( mat33 & ncs_m, vect3 & ncs_v,
                                 int iGiven )  {
    return cryst.AddNCSMatrix ( ncs_m,ncs_v,iGiven );
  }

  int  Root::GetNumberOfNCSMatrices()  {
    return cryst.GetNumberOfNCSMatrices();
  }

  int  Root::GetNumberOfNCSMates()  {
  // Returns the number of NCS mates not given in the file (iGiven==0)
    return cryst.GetNumberOfNCSMates();
  }

  bool  Root::GetNCSMatrix ( int NCSMatrixNo, // 0..N-1
                                     mat44 & ncs_m, int & iGiven )  {
    return cryst.GetNCSMatrix ( NCSMatrixNo,ncs_m,iGiven );
  }

  ERROR_CODE Root::ReadPDBAtom ( cpstr L )  {
  //   If string L belongs to the coordinate section
  // (records ATOM, SIGATM, ANISOU, SIGUIJ, TER, HETATM),
  // the correspondent information is retrieved and
  // stored in the dynamic Atom array. In parallel, the
  // structures of Model/Chain/Residue are generated and
  // referenced to the corresponding Atom.
  //   If processing of L was successful, the return is 0,
  // otherwise it returns the corresponding Error_XXX
  // code.
  //   If L does not belong to the coordinate section,
  // Error_WrongSection is returned.
  int        index,i;
  ERROR_CODE RC;

    if (!strncmp(L,"ATOM  ",6)) {

      index = nAtoms+1;  // index for the next atom in Atom array
      RC    = CheckAtomPlace ( index,L );
      if (!RC)  RC = atom[index-1]->ConvertPDBATOM ( index,L );

    } else if (!strncmp(L,"SIGATM",6)) {

      index = nAtoms;    // keep index!
      RC    = CheckAtomPlace ( index,L );
      if (!RC)  RC = atom[index-1]->ConvertPDBSIGATM ( index,L );

    } else if (!strncmp(L,"ANISOU",6)) {

      index = nAtoms;    // keep index
      RC    = CheckAtomPlace ( index,L );
      if (!RC)  RC = atom[index-1]->ConvertPDBANISOU ( index,L );

    } else if (!strncmp(L,"SIGUIJ",6)) {

      index = nAtoms;    // keep index
      RC    = CheckAtomPlace ( index,L );
      if (!RC)  RC = atom[index-1]->ConvertPDBSIGUIJ ( index,L );

    } else if (!strncmp(L,"TER   ",6)) {

      index = nAtoms+1;  // new place in Atom array
      RC    = CheckAtomPlace ( index,L );
      if (!RC)  RC = atom[index-1]->ConvertPDBTER ( index,L );

    } else if (!strncmp(L,"HETATM",6)) {

      index = nAtoms+1;  // new place in Atom array
      RC    = CheckAtomPlace ( index,L );
      if (!RC)  RC = atom[index-1]->ConvertPDBHETATM ( index,L );

    } else if (!strncmp(L,"MODEL ",6)) {

      modelCnt++;
      RC = SwitchModel ( L );
      for (i=0;(i<nModels) && (!RC);i++)
        if (model[i] && (model[i]!=crModel))  {
          if (crModel->serNum==model[i]->serNum)
            RC = Error_DuplicatedModel;
        }
  //    if (!RC)  {
  //      if (crModel->serNum!=modelCnt)
  //        RC = Error_DuplicatedModel;
  //    }

    } else if (!strncmp(L,"ENDMDL",6)) {

      crModel = NULL;
      crChain = NULL;
      crRes   = NULL;

      RC      = Error_NoError;

    } else
      return Error_WrongSection;

    return RC;

  }


  ERROR_CODE Root::ReadCIFAtom ( mmcif::PData CIFD )  {
  mmcif::PLoop Loop,LoopAnis;
  int          i,index,nATS;
  ERROR_CODE   RC;

    Loop = CIFD->GetLoop ( CIFCAT_ATOM_SITE );
    if (!Loop)  return Error_NoError;  // no atom coordinates in the file

    LoopAnis = CIFD->GetLoop ( CIFCAT_ATOM_SITE_ANISOTROP );
    nATS     = Loop->GetLoopLength();

    for (i=1;i<=nATS;i++)  {
      // nAtoms and i should always coincide at this point. This piece
      // of code was however left in order to reach identity with
      // ReadPDBAtom(..).
      index = nAtoms+1;  // index for the next atom in Atom array
      RC    = CheckAtomPlace ( index,Loop );
      if (!RC)  RC = atom[index-1]->GetCIF ( i,Loop,LoopAnis );
      if (RC && (RC!=Error_CIF_EmptyRow))  return RC;
    }
    if (Flags & MMDBF_AutoSerials)
      PDBCleanup ( PDBCLEAN_SERIAL );

    return Error_NoError;

  }

  int  Root::PutAtom ( int            index,
                       int            serNum,
                       const AtomName atomName,
                       const ResName  resName,
                       const ChainID  chainID,
                       int            seqNum,
                       const InsCode  insCode,
                       const AltLoc   altLoc,
                       const SegID    segID,
                       const Element  element )  {

  //   An atom with the specified properties is put into the
  // structure. The current model is used; if no model is
  // set (crModel==NULL), one is created. Coordinates and
  // other parameters of the atom need to be set separately.
  //
  //   If index is positive and there is already an atom at
  // this position in the system, the new atom will REPLACE
  // it. The corresponding residues are automatically
  // updated.
  //
  //   If index is null (=0), the new atom will be put on
  // the top of the structure, i.e. it will be put into
  // (index=nAtoms+1)-th position.
  //
  //   If index is negative, then the new atom is INSERTED
  // BEFORE the atom in the (-index)th position. For
  // saving the computational efforts, this WILL NOT cause
  // the recalculation of all atoms' serial numbers
  // according to their actual positions. It will be needed
  // however to put the things in order by calling
  // Root::OrderAtoms() at a certain point, especially
  // before writing an output ASCII file. NOTE that this
  // ordering is never done automatically.
  //
  //   Limitation: if PutAtom implies creating new
  // chains/residues, these are always created on the top
  // of existing chains/residues.


  int i,kndex,RC;

    kndex = index;

    if (kndex<0)  {  // the new atom is to be inserted

      kndex = -kndex;
      if (kndex>atmLen)
        ExpandAtomArray ( kndex+1000-atmLen );

      if (atom[kndex-1]!=NULL)  { // the position is occupied

        // expand the array if necessary
        if (nAtoms>=atmLen)
          ExpandAtomArray ( IMax(kndex,nAtoms)+1000-atmLen );

        // now shift all atoms from (kndex-1)th to the end of array.
        // note that this does not affect residues as they keep only
        // pointers on atoms
        for (i=nAtoms;i>=kndex;i--)  {
          atom[i] = atom[i-1];
          atom[i]->index = i+1;  // this is Ok because residues keep
                                 // POINTERS rather than indices!
        }
        atom[kndex-1] = NULL;
        nAtoms++;

      }

    }

    if (kndex==0)  kndex = nAtoms+1;

    if (!crModel)  SwitchModel ( 1 );


    RC = AllocateAtom ( kndex,chainID,chainID,resName,resName,
                        seqNum,seqNum,1,insCode,true );
    if (!RC)
      atom[kndex-1]->SetAtomName ( kndex,serNum,atomName,altLoc,
                                   segID,element );
    return RC;

  }


  int Root::PutAtom ( int   index,  // same meaning as above
                      PAtom A,      // pointer to completed atom class
                      int   serNum  // 0 means that the serial
                                    // number will be set equal
                                    // to "index". Otherwise,
                                    // the serial number is set
                                    // to the specified value
                    )  {
  int i,kndex,RC,sn;

    if (!A)  return -1;

    kndex = index;

    if (kndex<0)  {  // the new atom is to be inserted

      kndex = -kndex;

      if (kndex>atmLen)
        ExpandAtomArray ( kndex+1000-atmLen );

      if (atom[kndex-1]!=NULL)  { // the position is occupied

        // expand the array if necessary
        if (nAtoms>=atmLen)
          ExpandAtomArray ( IMax(kndex,nAtoms)+1000-atmLen );
        // now shift all atoms from (kndex-1)th to the end of array.
        // note that this does not affect residues as they keep only
        // pointers on atoms

        for (i=nAtoms;i>=kndex;i--)  {
          atom[i] = atom[i-1];
          atom[i]->index = i+1;  // this is Ok because residues keep
                                 // POINTERS rather than indices!
        }

        atom[kndex-1] = NULL;
        nAtoms++;

      }

    }

    if (kndex==0)  kndex = nAtoms+1;


    RC = AllocateAtom ( kndex,A->GetChainID(),A->GetLabelAsymID(),
                              A->GetResName(),A->GetLabelCompID(),
                              A->GetSeqNum (),A->GetLabelSeqID (),
                              A->GetLabelEntityID(),A->GetInsCode(),
                              true );

    if (serNum<=0)  sn = kndex;
              else  sn = serNum;
    if (!RC)  {
      atom[kndex-1]->Copy ( A );
      atom[kndex-1]->serNum = sn;
    }

    return RC;

  }

  int Root::CheckInAtom ( int index, // same meaning as above
                               PAtom  A  // pointer to completed
                                          // atom class
                             )  {
  int i,kndex;

    if (!A)  return -1;

    kndex = index;

    if (kndex<0)  {  // the new atom is to be inserted

      kndex = -kndex;

      if (kndex>atmLen)
        ExpandAtomArray ( kndex+1000-atmLen );

      if (atom[kndex-1]!=NULL)  { // the position is occupied

        // expand the array if necessary
        if (nAtoms>=atmLen)
          ExpandAtomArray ( IMax(kndex,nAtoms)+1000-atmLen );
        // now shift all atoms from (kndex-1)th to the end of array.
        // note that this does not affect residues as they keep only
        // pointers on atoms

        for (i=nAtoms;i>=kndex;i--)  {
          atom[i] = atom[i-1];
          if (atom[i])
            atom[i]->index = i+1;  // this is Ok because residues keep
                                   // POINTERS rather than indices!
        }

      }

      nAtoms++;

    } else  {
      if (kndex==0)      kndex = nAtoms + 1;  // add atom on the very top
      if (kndex>atmLen)  ExpandAtomArray ( kndex+1000-atmLen );
      if (kndex>nAtoms)  nAtoms = kndex;
      if (atom[kndex-1]) delete atom[kndex-1];
    }

    atom[kndex-1] = A;
    A->index = kndex;

    return 0;

  }

  int Root::CheckInAtoms ( int index, // same meaning as above
                                PPAtom A, // array of atoms to check in
                                int natms  // number of atoms to check in
                              )  {
  PPAtom A1;
  int     i,j,k,k1,kndex;

    if (!A)  return -1;

    A1    = NULL;
    kndex = index;

    if (kndex<0)  {  // the new atoms are to be inserted

      kndex = -kndex;

      if (nAtoms+natms>=atmLen)
        ExpandAtomArray ( IMax(kndex,nAtoms)+1000+natms-atmLen );

      if (kndex<nAtoms)
      A1 = new PAtom[natms];
      k = kndex-1;
      j = 0;
      for (i=0;i<natms;i++)
        if (A[i])  {
          if (atom[k])  A1[j++] = atom[k];
          atom[k] = A[i];
          atom[k]->index = k+1;
          k++;
        }

      if (j>0)  {
        // insert removed atoms into the gap
        nAtoms += j;
        k1      = k+j;
        for (i=nAtoms-1;i>=k1;i--)  {
          atom[i] = atom[i-j];
          if (atom[i])
            atom[i]->index = i+1;  // this is Ok because residues keep
                                   // POINTERS rather than indices!
        }
        for (i=0;i<j;i++)  {
          atom[k] = A1[i];
          atom[k]->index = k+1;
          k++;
        }
      }

      delete[] A1;

    } else  {

      if (kndex==0)      kndex = nAtoms + 1;  // add atom on the very top
      k = kndex + natms;
      if (k>atmLen)  ExpandAtomArray ( k+1000-atmLen );
      kndex--;
      for (i=0;i<natms;i++)
        if (A[i])  {
          if (atom[kndex]) delete atom[kndex];
          atom[kndex] = A[i];
          atom[kndex]->index = kndex+1;
          kndex++;
        }
      nAtoms = IMax(nAtoms,kndex);

    }

    return 0;

  }


  ERROR_CODE Root::SwitchModel ( cpstr L )  {
  int nM;

    if (!GetInteger(nM,&(L[10]),4))
      return Error_UnrecognizedInteger;

    return SwitchModel ( nM );

  }

  ERROR_CODE Root::SwitchModel ( int nM )  {
  PPModel Mdl;
  int      i;
  bool  Transfer;

    if (nM<=0)
      return Error_WrongModelNo;

    if (nM>nModels)  {
      if ((nModels==1) && model[0])  Transfer = (nAtoms<=0);
                               else  Transfer = false;
      Mdl = new PModel[nM];
      for (i=0;i<nModels;i++)
        Mdl[i] = model[i];
      for (i=nModels;i<nM;i++)
        Mdl[i] = NULL;
      if (model) delete[] model;
      model   = Mdl;
      nModels = nM;
      if (Transfer)  {
        model[nM-1] = model[0];
        model[0]    = NULL;
      }
    }

    if (!model[nM-1])
      model[nM-1] = newModel();
    model[nM-1]->SetMMDBManager ( PManager(this),nM );

    crModel = model[nM-1];
    crChain = NULL;  // new model - new chain
    crRes   = NULL;  // new chain - new residue

    return Error_NoError;

  }

  ERROR_CODE Root::CheckAtomPlace ( int index, cpstr L )  {
  //   This function gets the residue/chain information stored
  // in PDB string L (the records should start with the
  // keywords ATOM, SIGATM, ANISOU, SIGUIJ, TER, HETATM) and
  // sets the pointers crChain and crRes to the respective.
  // chain and residue. If there is no chain/residue to place
  // the atom in, these will be created.
  //   The function prepares place for the atom in the index-th
  // cell of the Atom array, expanding it as necessary. If the
  // corresponding element in the Atom array was not initialized,
  // a Atom class is created with reference to the current
  // residue.
  //   This function DOES NOT check the PDB string L for
  // atom keywords.
  ResName  resName;
  int      seqNum;
  ChainID  chainID;
  InsCode  insCode;

    // get the residue sequence number/ insert code
    if (!GetIntIns(seqNum,insCode,&(L[22]),4))  {
      if (strncmp(L,"TER   ",6))
            return Error_UnrecognizedInteger;
      else  { // we allow for empty TER card here
        seqNum  = 0;
        insCode[0] = char(1);  // unprintable symbol! used as
                               // flag that TER card does not
                               // have serial number
        insCode[1] = char(0);
      }
    }

    // get chain ID
    if (L[20]!=' ')  {
      chainID[0] = L[20];
      chainID[1] = L[21];
      chainID[2] = char(0);
    } else if (L[21]!=' ')  {
      chainID[0] = L[21];
      chainID[1] = char(0);
    } else
      chainID[0] = char(0);

    // get residue name
    strcpy_ncss ( resName,&(L[17]),3 );
    if ((!resName[0]) && (!strncmp(L,"TER   ",6)))  {
      insCode[0] = char(1);
      insCode[1] = char(0);
    }

    return AllocateAtom ( index ,chainID,chainID,resName,resName,
                          seqNum,seqNum,1,insCode,false );

  }

  ERROR_CODE Root::CheckAtomPlace ( int index, mmcif::PLoop Loop )  {
  //   Version of CheckAtomPlace(..) for reading from CIF file.
  ResName  resName,label_comp_id;
  int      seqNum ,label_seq_id,label_entity_id,RC,k,nM;
  ChainID  chainID,label_asym_id;
  InsCode  insCode;
  pstr     F;

    // Get the residue sequence number/insert code. They are
    // removed from the file after reading.
    k = index-1;
  //  if (!CIFGetInteger1(seqNum,Loop,CIFTAG_LABEL_SEQ_ID,k))
    if (!CIFGetInteger1(seqNum,Loop,CIFTAG_AUTH_SEQ_ID,k))
      CIFGetString  ( insCode,Loop,CIFTAG_NDB_HELIX_CLASS_PDB,k,
                      sizeof(InsCode),pstr("") );
    else  {
      F = Loop->GetString ( CIFTAG_GROUP_PDB,k,RC );
      if ((!F) || (RC)) return  Error_CIF_EmptyRow;
      if (strcmp(F,"TER"))  {
        seqNum = MinInt4;  // only at reading CIF we allow this
        CIFGetString ( insCode,Loop,CIFTAG_NDB_HELIX_CLASS_PDB,k,
                       sizeof(InsCode),pstr("") );
      } else  { // we allow for empty TER card here
        seqNum     = 0;
        insCode[0] = char(1);  // unprintable symbol! used as
                               // flag that TER card does not
                               // have serial number
        insCode[1] = char(0);
      }
    }

    CIFGetInteger1 ( label_seq_id   ,Loop,CIFTAG_LABEL_SEQ_ID   ,k );
    CIFGetInteger1 ( label_entity_id,Loop,CIFTAG_LABEL_ENTITY_ID,k );

    // get chain/residue ID
    CIFGetString ( chainID,Loop,CIFTAG_AUTH_ASYM_ID,k,
                   sizeof(ChainID),pstr("") );
    CIFGetString ( resName,Loop,CIFTAG_AUTH_COMP_ID,k,
                   sizeof(ResName),pstr("") );

    CIFGetString ( label_asym_id,Loop,CIFTAG_LABEL_ASYM_ID,k,
                   sizeof(ChainID),pstr("") );
    CIFGetString ( label_comp_id,Loop,CIFTAG_LABEL_COMP_ID,k,
                   sizeof(ResName),pstr("") );

    if (!resName[0])  strcpy ( resName,label_comp_id );

    if (!CIFGetInteger1(nM,Loop,CIFTAG_PDBX_PDB_MODEL_NUM,k))  {
      if (crModel)  {
        if (nM!=crModel->serNum)  SwitchModel ( nM );
      } else
        SwitchModel ( nM );
    }

    return AllocateAtom ( index ,chainID,label_asym_id,resName,
                          label_comp_id,seqNum,label_seq_id,
                          label_entity_id,insCode,false );

  }


  ERROR_CODE Root::AllocateAtom ( int           index,
                                  const ChainID chainID,
                                  const ChainID label_asym_id,
                                  const ResName resName,
                                  const ResName label_comp_id,
                                  int           seqNum,
                                  int           label_seq_id,
                                  int           label_entity_id,
                                  const InsCode insCode,
                                  bool          Replace )  {

    if ((!resName[0]) && (insCode[0]!=char(1)))
      return Error_EmptyResidueName;

    // check if there is a pointer to model
    if (!crModel)  {
      // the model pointer was not set. Check if there are
      // models already defined
      if (!model)
           SwitchModel ( 1 );  // creates a model
      else return Error_NoModel;
    }

    if (crChain && (insCode[0]!=char(1)))  {
      //   If crChain is not NULL, the model pointer was not
      // changed and we may try to keep using crChain as
      // pointer to the being-read chain. However, we must
      // check that the record still belongs to the same chain.
      //   All this does not work if insCode[0] is set to 1
      // which indicates a special case of 'TER' card without
      // parameters.
      if (enforceUniqueChID)  {
        // enforcing unique chain IDs should be used only in case
        // of multi-chain complexes where 1-letter chain IDs are
        // not enough to accomodate all chains. Then chains are
        // dynamically renamed like A0,A1,A2,.. etc. Therefore, we
        // check only first symbol here.
        if (chainID[0]!=crChain->chainID[0])
          crChain = NULL;  // the chain has to be changed
      } else if (strcmp(chainID,crChain->chainID))
        crChain = NULL;  // the chain has to be changed
    }
    if (!crChain) {
      // either the model or chain was changed  -- get a new chain
      if (allowDuplChID)
            crChain = crModel->CreateChain    ( chainID );
      else  crChain = crModel->GetChainCreate ( chainID,
                                                enforceUniqueChID );
      crRes = NULL;  // new chain - new residue
    }

    if (crRes && (insCode[0]!=char(1)))  {
      //   If crRes is not NULL, neither the model nor chain were
      // changed. Check if this record still belongs to the
      // same residue.
      //   All this does not work if insCode[0] is set to 1
      // which indicates a special case of 'TER' card without
      // parameters.
      if ((seqNum!=crRes->seqNum)         ||
           strcmp(insCode,crRes->insCode) ||
           strcmp(resName,crRes->name))
        crRes = NULL;  // the residue has to be changed
    }
    if (!crRes)  {
      // either the chain or residue was changed -- get a new residue
      crRes = crChain->GetResidueCreate ( resName,seqNum,insCode,
                                        Flags & MMDBF_IgnoreDuplSeqNum );
      if (!crRes)  return  Error_DuplicateSeqNum;
    }

    strcpy ( crRes->label_asym_id,label_asym_id );
    strcpy ( crRes->label_comp_id,label_comp_id );
    crRes->label_seq_id    = label_seq_id;
    crRes->label_entity_id = label_entity_id;

    // now check if there is place in the Atom array
    if (index>atmLen)
      // there is no place, expand Atom by 1000 atom places at once
      ExpandAtomArray ( index+1000-atmLen );
    nAtoms = IMax(nAtoms,index);

    // delete the to-be-replaced atom if there is any
    if (Replace && atom[index-1])  {
      delete atom[index-1];
      atom[index-1] = NULL;
    }
    if (!atom[index-1])  {
      atom[index-1] = newAtom();
      crRes->_AddAtom ( atom[index-1] );
      atom[index-1]->index = index;
    }

    return Error_NoError;

  }

  void Root::ExpandAtomArray ( int inc )  {
  // Expands the Atom array by adding more inc positions.
  // The length of Atom array is increased unconditionally.
  PPAtom Atom1;
  int     i;
    atmLen += inc;
    Atom1   = new PAtom[atmLen];
    for (i=0;i<nAtoms;i++)
      Atom1[i] = atom[i];
    for (i=nAtoms;i<atmLen;i++)
      Atom1[i] = NULL;
    if (atom) delete[] atom;
    atom = Atom1;
  }

  void Root::AddAtomArray ( int inc )  {
  // Checks if 'inc' atoms may be added into Atom array,
  // and if not, expands the Atom array such that to
  // allocate exactly 'inc' atoms more than is currently
  // contained.
  PPAtom Atom1;
  int     i;
    if (nAtoms+inc>atmLen)  {
      atmLen = nAtoms+inc;
      Atom1  = new PAtom[atmLen];
      for (i=0;i<nAtoms;i++)
        Atom1[i] = atom[i];
      for (i=nAtoms;i<atmLen;i++)
        Atom1[i] = NULL;
      if (atom) delete[] atom;
      atom = Atom1;
    }
  }


  ERROR_CODE Root::WritePDBASCII1 ( cpstr PDBLFName,
                                    io::GZ_MODE gzipMode )  {
  pstr FName;
    FName = getenv ( PDBLFName );
    if (FName)  return WritePDBASCII ( FName,gzipMode );
          else  return Error_NoLogicalName;
  }

  ERROR_CODE Root::WritePDBASCII ( cpstr PDBFileName,
                                   io::GZ_MODE gzipMode )  {
  io::File f;

    //  opening it in pseudo-text mode ensures that the line
    //  endings will correspond to the system MMDB is running on
    f.assign ( PDBFileName,true,false,gzipMode );
    FType = MMDB_FILE_PDB;

    if (f.rewrite())  {
      WritePDBASCII ( f );
      f.shut();
    } else
      return Error_CantOpenFile;

    return Error_NoError;

  }


  void  Root::WritePDBASCII ( io::RFile f )  {
  int  i;

    FType = MMDB_FILE_PDB;

    title.PDBASCIIDump ( f );

    i = 0;
    while (i<nModels)
      if (model[i])  break;
               else  i++;
    if (i<nModels)
      model[i]->PDBASCIIDumpPS ( f );

    // output cispep records
    for (i=0;i<nModels;i++)
      if (model[i])
        model[i]->PDBASCIIDumpCP ( f );

    SA      .PDBASCIIDump ( f );
    Footnote.PDBASCIIDump ( f );
    cryst   .PDBASCIIDump ( f );
    SB      .PDBASCIIDump ( f );

    for (i=0;i<nModels;i++)
      if (model[i])
        model[i]->PDBASCIIDump ( f );

    SC.PDBASCIIDump ( f );

    f.WriteLine ( pstr("END") );

  }


  ERROR_CODE Root::WriteCIFASCII1 ( cpstr CIFLFName,
                                    io::GZ_MODE gzipMode )  {
  pstr FName;
    FName = getenv ( CIFLFName );
    if (FName)  return WriteCIFASCII ( FName,gzipMode );
          else  return Error_NoLogicalName;
  }

  ERROR_CODE Root::WriteCIFASCII ( cpstr CIFFileName,
                                   io::GZ_MODE gzipMode )  {
  int  i;

    if (!CIF)  CIF = new mmcif::Data();
    CIF->SetStopOnWarning ( true );
    CIF->SetPrintWarnings ( (Flags & MMDBF_PrintCIFWarnings)!=0 );
    FType = MMDB_FILE_CIF;

    title.MakeCIF ( CIF );

    i = 0;
    while (i<nModels)
      if (model[i])  break;
               else  i++;
    if (i<nModels)
      model[i]->MakePSCIF ( CIF );

    cryst.MakeCIF ( CIF );

    for (i=0;i<nModels;i++)
      if (model[i])
        model[i]->MakeAtomCIF ( CIF );

    CIF->Optimize();
    CIF->WriteMMCIFData ( CIFFileName,gzipMode );

    return Error_NoError;

  }


  PAtom  Root::GetAtomI ( int index )  {
    if (index>nAtoms)  return NULL;
    if (index<1)       return NULL;
    if (!atom)         return NULL;
    return atom[index-1];
  }


  #define MMDBFLabel  "**** This is MMDB binary file ****"
  #define Edition     1

  ERROR_CODE Root::ReadMMDBF1 ( cpstr MMDBLFName,
                                io::GZ_MODE gzipMode )  {
  pstr FName;
    FName = getenv ( MMDBLFName );
    if (FName)  return ReadCoorFile ( FName,gzipMode );
          else  return Error_NoLogicalName;
  }

  ERROR_CODE Root::ReadMMDBF ( cpstr MMDBRootName,
                               io::GZ_MODE gzipMode )  {
  io::File   f;
  ERROR_CODE rc;

    f.assign ( MMDBRootName,false,true,gzipMode );
    FType = MMDB_FILE_Binary;
    if (f.reset(true))  {
      rc = ReadMMDBF ( f );
      f.shut();
    } else
      rc = Error_CantOpenFile;

    return rc;

  }

  ERROR_CODE Root::ReadMMDBF ( io::File & f )  {
  char  Label[100];
  byte  Version;

    FType = MMDB_FILE_Binary;
    f.ReadFile ( Label,sizeof(MMDBFLabel) );
    if (strncmp(Label,MMDBFLabel,sizeof(MMDBFLabel)))
      return Error_ForeignFile;

    f.ReadByte ( &Version );
    if (Version>Edition)
      return Error_WrongEdition;

    read ( f );

    return Error_NoError;

  }


  ERROR_CODE Root::WriteMMDBF1 ( cpstr MMDBLFName, io::GZ_MODE gzipMode )  {
  pstr FName;
    FName = getenv ( MMDBLFName );
    if (FName)  return WriteMMDBF ( FName,gzipMode );
          else  return Error_NoLogicalName;
  }

  /*
  ERROR_CODE Root::WriteMMDBF ( cpstr MMDBRootName, io::GZ_MODE gzipMode )  {
  io::File f;
  char     Label[100];
  byte     Version=Edition;

    f.assign ( MMDBRootName,false,true,gzipMode );
    FType = MMDB_FILE_Binary;
    if (f.rewrite())  {
      strcpy ( Label,MMDBFLabel );
      f.WriteFile ( Label,sizeof(MMDBFLabel) );
      f.WriteByte ( &Version );
      write ( f );
      f.shut();
    } else
      return Error_CantOpenFile;

    return Error_NoError;

  }
*/

  ERROR_CODE Root::WriteMMDBF ( cpstr       MMDBRootName,
                                io::GZ_MODE gzipMode )  {
  io::File f;
    f.assign ( MMDBRootName,false,true,gzipMode );
    if (f.rewrite())  {
      WriteMMDBF ( f );
      f.shut();
    } else
      return Error_CantOpenFile;
    return Error_NoError;
  }
  
  void Root::WriteMMDBF ( io::File & f )  {
  char  Label[100];
  byte  Version=Edition;
    FType = MMDB_FILE_Binary;
    strcpy ( Label,MMDBFLabel );
    f.WriteFile ( Label,sizeof(MMDBFLabel) );
    f.WriteByte ( &Version );
    write ( f );
  }
  

  pstr  Root::GetEntryID()  {
    return title.idCode;
  }

  void  Root::SetEntryID ( const IDCode idCode )  {
    strcpy ( title.idCode,idCode );
  }

  void Root::SetSyminfoLib ( cpstr syminfo_lib )  {
    cryst.SetSyminfoLib ( syminfo_lib );
  }

  pstr Root::GetSyminfoLib()  {
    return cryst.GetSyminfoLib();
  }

  int Root::SetSpaceGroup ( cpstr spGroup )  {
    return cryst.SetSpaceGroup ( spGroup );
  }

  pstr Root::GetSpaceGroup()  {
    return cryst.GetSpaceGroup();
  }

  pstr Root::GetSpaceGroupFix()  {
    return cryst.GetSpaceGroupFix();
  }

  void  Root::GetAtomStatistics ( RAtomStat AS )  {
  int i;
    AS.Init();
    for (i=0;i<nModels;i++)
      if (model[i])  model[i]->CalAtomStatistics ( AS );
    AS.Finish();
  }

  void Root::SetIgnoreSCALEi ( bool ignoreScalei )  {
    cryst.ignoreScalei = ignoreScalei;
  }

  void Root::SetCell ( realtype cell_a,
                       realtype cell_b,
                       realtype cell_c,
                       realtype cell_alpha,
                       realtype cell_beta,
                       realtype cell_gamma,
                       int      OrthCode )  {
    cryst.SetCell ( cell_a,cell_b,cell_c,cell_alpha,cell_beta,
                    cell_gamma,OrthCode );
  }

  void Root::PutCell ( realtype cell_a,
                       realtype cell_b,
                       realtype cell_c,
                       realtype cell_alpha,
                       realtype cell_beta,
                       realtype cell_gamma,
                       int      OrthCode )  {
    cryst.PutCell ( cell_a,cell_b,cell_c,cell_alpha,cell_beta,
                    cell_gamma,OrthCode );
  }

  int  Root::GetCell ( realtype & cell_a,
                            realtype & cell_b,
                            realtype & cell_c,
                            realtype & cell_alpha,
                            realtype & cell_beta,
                            realtype & cell_gamma,
                            realtype & vol,
                            int      & OrthCode )  {
    if (cryst.WhatIsSet & CSET_CellParams)  {
      cryst.GetCell ( cell_a,cell_b,cell_c,cell_alpha,cell_beta,
                      cell_gamma,vol );
      OrthCode = cryst.NCode + 1;
      return 1;
    } else {
      cell_a     = 0.0;    cell_b    = 0.0;    cell_c     = 0.0;
      cell_alpha = 0.0;    cell_beta = 0.0;    cell_gamma = 0.0;
      vol        = 0.0;    OrthCode  = 0;
      return 0;
    }
  }

  int  Root::GetRCell ( realtype & cell_as,
                             realtype & cell_bs,
                             realtype & cell_cs,
                             realtype & cell_alphas,
                             realtype & cell_betas,
                             realtype & cell_gammas,
                             realtype & vols,
                             int      & OrthCode )  {
    if (cryst.WhatIsSet & CSET_CellParams)  {
      cryst.GetRCell ( cell_as,cell_bs,cell_cs,cell_alphas,cell_betas,
                       cell_gammas,vols );
      OrthCode = cryst.NCode + 1;
      return 1;
    } else {
      cell_as     = 0.0;    cell_bs    = 0.0;    cell_cs     = 0.0;
      cell_alphas = 0.0;    cell_betas = 0.0;    cell_gammas = 0.0;
      vols        = 0.0;    OrthCode   = 0;
      return 0;
    }
  }

  int Root::GetNumberOfSymOps()  {
    if (cryst.WhatIsSet & CSET_SpaceGroup)
          return cryst.GetNumberOfSymOps();
    else  return 0;
  }

  pstr Root::GetSymOp ( int Nop )  {
    return cryst.GetSymOp ( Nop );
  }


  void Root::GetROMatrix ( mat44 & RO )  {
    Mat4Copy ( cryst.RO,RO );
  }

  int Root::GetTMatrix ( mat44 & TMatrix, int Nop,
                         int cellshift_a, int cellshift_b,
                         int cellshift_c )  {
  //  GetTMatrix(..) calculates and returns the coordinate transformation
  //  matrix, which converts orthogonal coordinates according to
  //  the symmetry operation number Nop and places them into unit cell
  //  shifted by cellshift_a a's, cellshift_b b's and cellshift_c c's.
  //
  //  Return 0 means everything's fine,
  //         1 there's no symmetry operation Nop defined
  //         2 fractionalizing/orthogonalizing matrices were not
  //           calculated
  //         3 cell parameters were not set up.
    return cryst.GetTMatrix ( TMatrix,Nop,cellshift_a,cellshift_b,
                              cellshift_c,NULL );
  }


  int Root::GetUCTMatrix ( mat44 & TMatrix, int Nop,
                           realtype x, realtype y, realtype z,
                           int cellshift_a, int cellshift_b,
                           int cellshift_c )  {
  //  GetUCTMatrix(..) calculates and returns the coordinate
  //  transformation matrix, which converts orthogonal coordinates
  //  according to the symmetry operation number Nop. Translation
  //  part of the resulting matrix is being chosen such that point
  //  (x,y,z) has least distance to the center of primary (333)
  //  unit cell, and then it is shifted by cellshift_a a's,
  //  cellshift_b b's and cellshift_c c's.
  //
  //  Return 0 means everything's fine,
  //         1 there's no symmetry operation Nop defined
  //         2 fractionalizing/orthogonalizing matrices were not
  //           calculated
  //         3 cell parameters were not set up.
    return cryst.GetUCTMatrix ( TMatrix,Nop,x,y,z,
                                cellshift_a,cellshift_b,cellshift_c,
                                NULL );
  }


  int Root::GetFractMatrix ( mat44 & TMatrix, int Nop,
                             int cellshift_a, int cellshift_b,
                             int cellshift_c )  {
  //  GetFractMatrix(..) calculates and returns the coordinate
  //  transformation matrix, which converts fractional coordinates
  //  according to the symmetry operation number Nop and places them
  //  into unit cell shifted by cellshift_a a's, cellshift_b b's and
  //  cellshift_c c's.
  //
  //  Return 0 means everything's fine,
  //         1 there's no symmetry operation Nop defined
  //         2 fractionalizing/orthogonalizing matrices were not
  //           calculated
  //         3 cell parameters were not set up.
    return cryst.GetFractMatrix ( TMatrix,Nop,cellshift_a,cellshift_b,
                                  cellshift_c,NULL );
  }

  int  Root::GetSymOpMatrix ( mat44 & TMatrix, int Nop )  {
  //
  //  GetSymOpMatrix(..) returns the transformation matrix for
  //  Nop-th symmetry operator in the space group
  //
  //  Return 0 means everything's fine,
  //         1 there's no symmetry operation Nop defined
  //         2 fractionalizing/orthogonalizing matrices were not
  //           calculated
  //         3 cell parameters were not set up.
  //
    return cryst.GetSymOpMatrix ( TMatrix,Nop );
  }

  //  -------------  User-Defined Data  ------------------------

  int  Root::RegisterUDInteger ( UDR_TYPE udr_type, cpstr UDDataID )  {
    return udRegister.RegisterUDInteger ( udr_type,UDDataID );
  }

  int  Root::RegisterUDReal ( UDR_TYPE udr_type, cpstr UDDataID )  {
    return udRegister.RegisterUDReal ( udr_type,UDDataID );
  }

  int  Root::RegisterUDString ( UDR_TYPE udr_type, cpstr UDDataID )  {
    return udRegister.RegisterUDString ( udr_type,UDDataID );
  }

  int  Root::GetUDDHandle ( UDR_TYPE udr_type, cpstr UDDataID )  {
    return udRegister.GetUDDHandle ( udr_type,UDDataID );
  }



  //  ----------------------------------------------------------

  int Root::DeleteAllModels()  {
  int i,k;
    Exclude = false;
    k = 0;
    for (i=0;i<nModels;i++)  {
      if (model[i])  {
        delete model[i];
        model[i] = NULL;
        k++;
      }
    }
    Exclude = true;
    FinishStructEdit();
    return k;
  }

  bool Root::GetNewChainID ( int modelNo, ChainID chID,
                                     int length )  {
    if ((modelNo>=1) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return  model[modelNo-1]->GetNewChainID ( chID,length );
    }
    return false;
  }

  //  -------------------------------------------------------------

  PMask Root::GetSelMask ( int selHnd )  {
  UNUSED_ARGUMENT(selHnd);
    return NULL;
  }

  //  -------------------------------------------------------------

  int  Root::GetNofExpDataRecs()  {
    return title.expData.Length();
  }

  pstr  Root::GetExpDataRec ( int recNo )  {
  PExpData  expData;
    expData = PExpData(title.expData.GetContainerClass(recNo));
    if (expData)  return expData->Line;
    return NULL;
  }


  //  -------------------------------------------------------------

  int  Root::GetNofMdlTypeRecs()  {
    return title.mdlType.Length();
  }

  pstr  Root::GetMdlTypeRec ( int recNo )  {
  PMdlType  mdlType;
    mdlType = PMdlType(title.mdlType.GetContainerClass(recNo));
    if (mdlType)  return mdlType->Line;
    return NULL;
  }


  //  -------------------  Stream functions  ----------------------

  void  Root::Copy ( PRoot MMDBRoot )  {
  int i;

    title.Copy ( &MMDBRoot->title );
    cryst.Copy ( &MMDBRoot->cryst );

    //   It is important to copy atoms _before_ models,
    // residues and chains!
    Flags  = MMDBRoot->Flags;
    nAtoms = MMDBRoot->nAtoms;
    atmLen = nAtoms;
    if (nAtoms>0)  {
      atom = new PAtom[atmLen];
      for (i=0;i<nAtoms;i++)
        if (MMDBRoot->atom[i])  {
          atom[i] = newAtom();
          atom[i]->Copy ( MMDBRoot->atom[i] );
          atom[i]->index = i+1;
          // the internal atom references are installed
          // by residue classes when they are copied in
          // model->chain below
        } else
          atom[i] = NULL;
    }

    nModels = MMDBRoot->nModels;
    if (nModels>0)  {
      model = new PModel[nModels];
      for (i=0;i<nModels;i++)  {
        if (MMDBRoot->model[i])  {
          model[i] = newModel();
          model[i]->SetMMDBManager ( PManager(this),i+1 );
          model[i]->_copy ( MMDBRoot->model[i] );
        } else
          model[i] = NULL;
      }
    }

    SA      .Copy ( &MMDBRoot->SA       );
    Footnote.Copy ( &MMDBRoot->Footnote );
    SB      .Copy ( &MMDBRoot->SB       );
    SC      .Copy ( &MMDBRoot->SC       );

    if (MMDBRoot->CIF)  {
      CIF = new mmcif::Data;
      CIF->Copy ( MMDBRoot->CIF );
    }

  }



  // -------  user-defined data handlers

  int  Root::PutUDData ( int UDDhandle, int iudd )  {
    if (UDDhandle & UDRF_HIERARCHY)
          return  UDData::putUDData ( UDDhandle,iudd );
    else  return  UDDATA_WrongUDRType;
  }

  int  Root::PutUDData ( int UDDhandle, realtype rudd )  {
    if (UDDhandle & UDRF_HIERARCHY)
          return  UDData::putUDData ( UDDhandle,rudd );
    else  return  UDDATA_WrongUDRType;
  }

  int  Root::PutUDData ( int UDDhandle, cpstr sudd )  {
    if (UDDhandle & UDRF_HIERARCHY)
          return  UDData::putUDData ( UDDhandle,sudd );
    else  return  UDDATA_WrongUDRType;
  }

  int  Root::GetUDData ( int UDDhandle, int & iudd )  {
    if (UDDhandle & UDRF_HIERARCHY)
          return  UDData::getUDData ( UDDhandle,iudd );
    else  return  UDDATA_WrongUDRType;
  }

  int  Root::GetUDData ( int UDDhandle, realtype & rudd )  {
    if (UDDhandle & UDRF_HIERARCHY)
          return  UDData::getUDData ( UDDhandle,rudd );
    else  return  UDDATA_WrongUDRType;
  }

  int  Root::GetUDData ( int UDDhandle, pstr sudd, int maxLen )  {
    if (UDDhandle & UDRF_HIERARCHY)
          return  UDData::getUDData ( UDDhandle,sudd,maxLen );
    else  return  UDDATA_WrongUDRType;
  }

  int  Root::GetUDData ( int UDDhandle, pstr & sudd )  {
    if (UDDhandle & UDRF_HIERARCHY)
          return  UDData::getUDData ( UDDhandle,sudd );
    else  return  UDDATA_WrongUDRType;
  }


  pstr Root::GetStructureTitle ( pstr & L )  {
    return title.GetStructureTitle ( L );
  }

  void  Root::SetCompactBinary()  {
  // leaves only coordinates in binary files
  int i;
    SetFlag ( MMDBF_MakeCompactBinary );
    for (i=0;i<nAtoms;i++)
      if (atom[i])  atom[i]->SetCompactBinary();
  }


  void  Root::write ( io::RFile f )  {
  int  i,k;
  byte Version=2;

    f.WriteByte ( &Version );
    f.WriteWord ( &Flags   );
    
    f.WriteInt  ( &nAtoms );
    for (i=0;i<nAtoms;i++)  {
      if (atom[i])  k = 1;
              else  k = 0;
      f.WriteInt ( &k );
      if (atom[i]) atom[i]->write ( f );
    }

    f.WriteInt ( &nModels );
    for (i=0;i<nModels;i++)  {
      if (model[i])  k = 1;
               else  k = 0;
      f.WriteInt ( &k );
      if (model[i]) model[i]->write ( f );
    }
    
    if (Flags & MMDBF_MakeCompactBinary)  {
      
      f.WriteTerLine ( title.idCode,false      );
      f.WriteReal    ( &title.resolution );
      title.title.write ( f );
      cryst      .write ( f );
      
    } else  {

      UDData::write ( f );
  
      title     .write ( f );
      cryst     .write ( f );
      udRegister.write ( f );
      DefPath   .write ( f );
  
      SA      .write ( f );
      Footnote.write ( f );
      SB      .write ( f );
      SC      .write ( f );
  
      StreamWrite ( f,CIF );
    
    }

  }

  void  Root::read ( io::RFile f )  {
  int  i,k;
  byte Version;

    ResetManager  ();
    FreeFileMemory();

    f.ReadByte ( &Version );
    f.ReadWord ( &Flags   );
    
    //   It is important to read atoms before models,
    // residues and chains!
    f.ReadInt  ( &nAtoms );
    atmLen = nAtoms;
    if (nAtoms>0)  {
      atom = new PAtom[atmLen];
      for (i=0;i<nAtoms;i++)  {
        f.ReadInt ( &k );
        if (k)  {
          atom[i] = newAtom();
          atom[i]->read ( f );
          // the internal atom references are installed
          // by residue classes when they are read in
          // model->chain below
        } else
          atom[i] = NULL;
      }
    }

    f.ReadInt ( &nModels );
    if (nModels>0)  {
      model = new PModel[nModels];
      for (i=0;i<nModels;i++)  {
        f.ReadInt ( &k );
        if (k)  {
          model[i] = newModel();
          model[i]->SetMMDBManager ( PManager(this),0 );
          model[i]->read ( f );
        } else
          model[i] = NULL;
      }
    }
    if (Flags & MMDBF_MakeCompactBinary)  {
      
      f.ReadTerLine ( title.idCode,false      );
      f.ReadReal    ( &title.resolution );
      title.title.read ( f );
      cryst      .read ( f );
      
    } else  {
  
      UDData::read ( f );
  
      title     .read ( f );
      cryst     .read ( f );
      udRegister.read ( f );
      DefPath   .read ( f );
  
      SA      .read ( f );
      Footnote.read ( f );
      SB      .read ( f );
      SC      .read ( f );
  
      StreamRead ( f,CIF );
      
    }

  }


  MakeStreamFunctions(Root)


  int isMMDBBIN ( cpstr FName, io::GZ_MODE gzipMode )  {
  io::File f;
  int   rc;

    f.assign ( FName,false,true,gzipMode );
    if (f.reset(true))  {
      rc = isMMDBBIN ( f );
      f.shut();
    } else
      rc = -1;

    return rc;

  }

  int isMMDBBIN ( io::RFile f )  {
  char  Label[100];
  byte  Version;

    if (f.FileEnd())
      return Error_EmptyFile;

    f.ReadFile ( Label,sizeof(MMDBFLabel) );
    if (strncmp(Label,MMDBFLabel,sizeof(MMDBFLabel)))
      return 1;

    f.ReadByte ( &Version );

    if (Version>Edition)  return 2;
                    else  return 0;

  }


  int isPDB ( cpstr FName, io::GZ_MODE gzipMode,
              bool IgnoreBlankLines )  {
  io::File f;
  int      rc;

    f.assign ( FName,false,false,gzipMode );
    if (f.reset(true))  {
      //  opening it in pseudo-binary mode helps reading various
      //  line terminators for files coming from different platforms
      rc = isPDB ( f,IgnoreBlankLines );
      f.shut();
    } else
      rc = -1;

    return rc;

  }

  int isPDB ( io::RFile f, bool IgnoreBlankLines )  {
  char    S[256];
  int     i;
  bool Done;

    if (f.FileEnd())
      return Error_EmptyFile;

    do {
      Done = true;
      f.ReadLine ( S,sizeof(S)-1 );
      if (IgnoreBlankLines)  {
        i = 0;
        while (S[i] && (S[i]==' '))  i++;
        if (!S[i])  Done = false;
      }
    } while ((!f.FileEnd()) && (!Done));

    PadSpaces  ( S,80 );
    if (!strncasecmp(S,"HEADER",6))  return 0;
    if (!strncasecmp(S,"OBSLTE",6))  return 0;
    if (!strncasecmp(S,"TITLE ",6))  return 0;
    if (!strncasecmp(S,"CAVEAT",6))  return 0;
    if (!strncasecmp(S,"COMPND",6))  return 0;
    if (!strncasecmp(S,"SOURCE",6))  return 0;
    if (!strncasecmp(S,"KEYWDS",6))  return 0;
    if (!strncasecmp(S,"EXPDTA",6))  return 0;
    if (!strncasecmp(S,"AUTHOR",6))  return 0;
    if (!strncasecmp(S,"REVDAT",6))  return 0;
    if (!strncasecmp(S,"SPRSDE",6))  return 0;
    if (!strncasecmp(S,"JRNL  ",6))  return 0;
    if (!strncasecmp(S,"REMARK",6))  return 0;
    if (!strncasecmp(S,"DBREF ",6))  return 0;
    if (!strncasecmp(S,"SEQADV",6))  return 0;
    if (!strncasecmp(S,"SEQRES",6))  return 0;
    if (!strncasecmp(S,"MODRES",6))  return 0;
    if (!strncasecmp(S,"HET   ",6))  return 0;
    if (!strncasecmp(S,"HETNAM",6))  return 0;
    if (!strncasecmp(S,"HETSYN",6))  return 0;
    if (!strncasecmp(S,"FORMUL",6))  return 0;
    if (!strncasecmp(S,"HELIX ",6))  return 0;
    if (!strncasecmp(S,"SHEET ",6))  return 0;
    if (!strncasecmp(S,"TURN  ",6))  return 0;
    if (!strncasecmp(S,"SSBOND",6))  return 0;
    if (!strncasecmp(S,"LINK  ",6))  return 0;
    if (!strncasecmp(S,"HYDBND",6))  return 0;
    if (!strncasecmp(S,"SLTBRG",6))  return 0;
    if (!strncasecmp(S,"CISPEP",6))  return 0;
    if (!strncasecmp(S,"SITE  ",6))  return 0;
    if (!strncasecmp(S,"CRYST1",6))  return 0;
    if (!strncasecmp(S,"CRYST ",6))  return 0;
    if (!strncasecmp(S,"ORIGX1",6))  return 0;
    if (!strncasecmp(S,"ORIGX2",6))  return 0;
    if (!strncasecmp(S,"ORIGX3",6))  return 0;
    if (!strncasecmp(S,"SCALE1",6))  return 0;
    if (!strncasecmp(S,"SCALE2",6))  return 0;
    if (!strncasecmp(S,"SCALE3",6))  return 0;
    if (!strncasecmp(S,"MTRIX1",6))  return 0;
    if (!strncasecmp(S,"MTRIX2",6))  return 0;
    if (!strncasecmp(S,"MTRIX3",6))  return 0;
    if (!strncasecmp(S,"TVECT ",6))  return 0;
    if (!strncasecmp(S,"MODEL ",6))  return 0;
    if (!strncasecmp(S,"ATOM  ",6))  return 0;
    if (!strncasecmp(S,"SIGATM",6))  return 0;
    if (!strncasecmp(S,"ANISOU",6))  return 0;
    if (!strncasecmp(S,"SIGUIJ",6))  return 0;
    if (!strncasecmp(S,"TER   ",6))  return 0;
    if (!strncasecmp(S,"HETATM",6))  return 0;
    if (!strncasecmp(S,"ENDMDL",6))  return 0;
    if (!strncasecmp(S,"CONECT",6))  return 0;
    if (!strncasecmp(S,"MASTER",6))  return 0;
    if (!strncasecmp(S,"END   ",6))  return 0;
    if (!strncasecmp(S,"USER  ",6))  return 0;

    return  1;

  }

}  // namespace mmdb
