//  $Id: mmdb_rwbrook.cpp $
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
//    16.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  MMDB_RWBrook  <implementation>
//       ~~~~~~~~~
//  **** Project :  MacroMolecular Data Base (MMDB), "new rwbrook"
//       ~~~~~~~~~
//  **** Functions : mmdb_f_init_     ( initializer                       )
//       ~~~~~~~~~~~ mmdb_f_quit_     ( disposer                          )
//                   autoserials_     ( switch to the autoserials regime  )
//                   setreadcoords_   ( switch for reading coordinates    )
//                   simrwbrook_      ( simulates old RWBROOK printout    )
//                   mmdb_f_openl_    ( associates a unit with a file     )
//                   mmdb_f_open_     ( associates a unit with a file     )
//                   mmdb_f_copy_     ( copies contents of units          )
//                   mmdb_f_delete_   ( deletes part of a unit            )
//                   mmdb_f_settype_  ( changes type of file and r/w mode )
//                   mmdb_f_setname_  ( changes file name                 )
//                   mmdb_f_write_    ( writes a data structure into file )
//                   mmdb_f_close_    ( closes and disposes a data str-re )
//                   mmdb_f_advance_  ( advances the internal pointer     )
//                   mmdb_f_rewd_     ( sets internal pointer on the top  )
//                   mmdb_f_bksp_     ( shifts int-l pointer 1 atom back  )
//                   mmdb_f_atom_     ( reads/writes atom properties      )
//                   mmdb_f_coord_    ( reads/writes atom coordinates     )
//                   mmdb_f_setcell_  ( sets the crystal cell parameters  )
//                   mmdb_f_wbspgrp_  ( sets the space group              )
//                   mmdb_f_rbspgrp_  ( gets the space group              )
//                   mmdb_f_wbcell_   ( sets the crystal cell parameters  )
//                   mmdb_f_rbcell_   ( gets the crystal cell parameters  )
//                   mmdb_f_rbcelln_  ( gets the crystal cell parameters  )
//                   mmdb_f_rbrcel_   ( gets the recipricol cell          )
//                   mmdb_f_rborf_    ( returns or fill transf. matrices  )
//                   mmdb_f_orthmat_  ( calc. standard othogonalisations  )
//                   mmdb_f_cvanisou_ ( converts between cryst-c units    )
//                   mmdb_f_wremark_  ( writes a remark statement         )
//                   mmdb_f_setter
//                   mmdb_f_sethet
//                   mmdb_f_getnofncsmates_
//                   rberrstop_       ( error messenger                   )
//                   rbcheckerr_      ( a simple error messenger          )
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#include "string.h"
#include "stdlib.h"
#include "math.h"

#include "mmdb_rwbrook.h"
#include "mmdb_manager.h"
#include "mmdb_tables.h"
#include "hybrid_36.h"

//  ==========================  Channel  ===========================

DefineClass(Channel)

class Channel {

  public :

    int            nUnit;       // unit number
    int            nType;       // unit type: 0- PDB; 1- CIF; 2- binary
    int            nRead;       // 0: input, 1: output
    mmdb::PManager MMDBManager; // MMDB manager
    mmdb::pstr     FName;       // file name
    int            fPos;        // "position" in the file
    int            ErrCode;     // error code of last operation
    bool           FAutoSer;    // autoserials flag for reading PDB
    bool           FReadCoords; // flag to read coordinate section
    bool           FSimRWBROOK; // flag to simulate old RWBROOK's printout

    Channel ();
    ~Channel();

    void    Dispose();
    void    Init   ();

    void    SetFileType ( mmdb::pstr FType );
    void    SetFileName ( mmdb::pstr FileName, int FNameLen );
    void    IdentifyFile( mmdb::pstr ExistingFName );

    bool    EndOfFile   ();
    mmdb::PAtom * GetAtomArray();
    mmdb::PAtom   GetAtomI    ( int index );

    mmdb::PCryst GetCryst ();

    bool   areCrystMatrices();
    void   Frac2Orth   (
                mmdb::realtype x,    mmdb::realtype y,    mmdb::realtype z,
                mmdb::realtype & xx, mmdb::realtype & yy, mmdb::realtype & zz );
    void   Orth2Frac   (
                mmdb::realtype x,    mmdb::realtype y,    mmdb::realtype z,
                mmdb::realtype & xx, mmdb::realtype & yy, mmdb::realtype & zz );
    void   Cryst2Orth  ( mmdb::rvector U );
    void   Orth2Cryst  ( mmdb::rvector U );
    int    SetCell     ( mmdb::realtype cell_a,
                         mmdb::realtype cell_b,
                         mmdb::realtype cell_c,
                         mmdb::realtype cell_alpha,
                         mmdb::realtype cell_beta,
                         mmdb::realtype cell_gamma,
                         int      OrthCode );
    int    PutCell     ( mmdb::realtype cell_a,
                         mmdb::realtype cell_b,
                         mmdb::realtype cell_c,
                         mmdb::realtype cell_alpha,
                         mmdb::realtype cell_beta,
                         mmdb::realtype cell_gamma,
                         int      OrthCode );
    int    SetSpGroup  ( mmdb::pstr     spGroup );
    int    GetSpGroup  ( mmdb::pstr     spGroup );
    int    GetCell     ( mmdb::realtype & cell_a,
                         mmdb::realtype & cell_b,
                         mmdb::realtype & cell_c,
                         mmdb::realtype & cell_alpha,
                         mmdb::realtype & cell_beta,
                         mmdb::realtype & cell_gamma,
                         mmdb::realtype & cell_v,
                         int      & OrthCode );
    int    GetRCell    ( mmdb::realtype & cell_as,
                         mmdb::realtype & cell_bs,
                         mmdb::realtype & cell_cs,
                         mmdb::realtype & cell_alphas,
                         mmdb::realtype & cell_betas,
                         mmdb::realtype & cell_gammas,
                         mmdb::realtype & cell_vs );

    int  GetNumberOfNCSMates();  // Returns the number of
                                 // NCS mates not given in
                                 // the file (iGiven==0)


    void MakeCoordStructure();
    void Read ();
    void Write();

    void GetInputBuffer ( mmdb::pstr Line, int & count );

  protected :

    void TranslateError();

};

Channel::Channel()  {
  Init();
}

Channel::~Channel()  {
  Dispose();
}

void Channel::Init()  {
  nUnit       = -1;
  nType       = -1;
  nRead       = 0;
  MMDBManager = NULL;
  FName       = NULL;
  ErrCode     = 0;
  fPos        = 0;
  FAutoSer    = false;
  FReadCoords = true;
  FSimRWBROOK = false;
}

void Channel::Dispose()  {
  if (MMDBManager)  delete MMDBManager;
  if (FName)        delete[] FName;
  MMDBManager = NULL;
  FName       = NULL;
  nUnit       = -1;
  nType       = -1;
  nRead       = 0;
  ErrCode     = 0;
  fPos        = 0;
}


void Channel::SetFileType ( mmdb::pstr FType )  {
  switch (FType[0])  {
    default  :
    case ' ' :  if (nRead==0)
                     nType = -1;                  // auto at reading
                else if (MMDBManager)
                     nType = MMDBManager->GetFileType(); // auto at writing
                else nType = -1;
              break;
    case 'P' :  nType = 0;   break;  // PDB
    case 'C' :  nType = 1;   break;  // CIF
    case 'B' :  nType = 2;   break;  // BIN
  }
}

void Channel::IdentifyFile ( mmdb::pstr ExistingFName )  {
  if (nType==-1)  {
    if (ExistingFName)  {
      if  (mmdb::isMMDBBIN(ExistingFName)==0)  nType = 2;
      else if (mmdb::isPDB(ExistingFName,mmdb::io::GZM_CHECK,true)==0)
                                         nType = 0;
      else if (mmdb::mmcif::isCIF(ExistingFName)==0)  nType = 1;
                                   else  nType = -2;  // unidentified
    } else  {
      if (MMDBManager)  {
        if (MMDBManager->GetFileType()<0)
              nType = 0;                  // PDB
        else  nType = MMDBManager->GetFileType(); // same as it was on last input
      } else  nType = 0;
    }
  }
}

void Channel::SetFileName ( mmdb::pstr FileName, int FNameLen )  {
  if (FName)  delete[] FName;
  FName = new char[FNameLen+1];
  strncpy ( FName,FileName,FNameLen );
  FName[FNameLen] = char(0);
}

void Channel::MakeCoordStructure()  {
  if (MMDBManager)
    MMDBManager->Delete ( mmdb::MMDBFCM_All );
  else  {
    MMDBManager = new mmdb::Manager();
    MMDBManager->SetFlag ( mmdb::MMDBF_AllowDuplChainID );
  }
}

int Channel::GetNumberOfNCSMates()  {
// Returns the number of NCS mates not given in the file (iGiven==0)
  if (!MMDBManager)  return RWBERR_NoData;
  return MMDBManager->GetNumberOfNCSMates();
}

void Channel::Read()  {
int RC;

  ErrCode = -2;
  if (!FName)  return;

  MakeCoordStructure();

  IdentifyFile ( FName );

  if (FAutoSer)     MMDBManager->SetFlag    ( mmdb::MMDBF_AutoSerials );
           else     MMDBManager->RemoveFlag ( mmdb::MMDBF_AutoSerials );
  if (FReadCoords)  MMDBManager->RemoveFlag ( mmdb::MMDBF_NoCoordRead );
              else  MMDBManager->SetFlag    ( mmdb::MMDBF_NoCoordRead );
  if (FSimRWBROOK)  MMDBManager->SetFlag    ( mmdb::MMDBF_SimRWBROOK  );
              else  MMDBManager->RemoveFlag ( mmdb::MMDBF_SimRWBROOK  );

  MMDBManager->SetFlag ( mmdb::MMDBF_IgnoreDuplSeqNum |
                         mmdb::MMDBF_IgnoreBlankLines |
                         mmdb::MMDBF_IgnoreRemarks    |
                         mmdb::MMDBF_IgnoreNonCoorPDBErrors |
                         mmdb::MMDBF_AllowDuplChainID );

  switch (nType)  {
    default : nType   = 0;  // nType=-2: unidentified: try PDB
    case  0 : ErrCode = MMDBManager->ReadPDBASCII ( FName );  break;
    case  1 : ErrCode = MMDBManager->ReadCIFASCII ( FName );  break;
    case  2 : ErrCode = MMDBManager->ReadMMDBF    ( FName );  break;
  }
  if (ErrCode==0)  {
    RC = MMDBManager->CrystReady();
    switch (RC)  {
      case mmdb::CRRDY_NoTransfMatrices : ErrCode = RWBERR_NoMatrices;   break;
      case mmdb::CRRDY_Unchecked        : ErrCode = RWBERR_NoCheck;      break;
      case mmdb::CRRDY_Ambiguous        : ErrCode = RWBERR_Disagreement; break;
      case mmdb::CRRDY_NoCell           : ErrCode = RWBERR_NoCellParams; break;
      default : ;
    }
  }
  fPos = 0;  // begining of the file
  TranslateError();
}

void Channel::Write()  {
  ErrCode = -3;
  if ((!MMDBManager) || (!FName))  return;
  IdentifyFile ( FName );
  switch (nType)  {
    default : nType   = 0;  // nType=-2: unidentified: make PDB
    case  0 : ErrCode = MMDBManager->WritePDBASCII ( FName );  break;
    case  1 : ErrCode = MMDBManager->WriteCIFASCII ( FName );  break;
    case  2 : ErrCode = MMDBManager->WriteMMDBF    ( FName );  break;
  }
  // we do not change fPos here!
  TranslateError();
}

void  Channel::TranslateError()  {

  switch (ErrCode)  {

    case mmdb::Error_CantOpenFile        : ErrCode = RWBERR_CantOpenFile;     break;
    case mmdb::Error_UnrecognizedInteger : ErrCode = RWBERR_WrongInteger;     break;
    case mmdb::Error_NoData             : ErrCode = RWBERR_NotACIFFile;       break;
    case mmdb::Error_WrongModelNo       : ErrCode = RWBERR_WrongModelNo;      break;
    case mmdb::Error_DuplicatedModel    : ErrCode = RWBERR_DuplicatedModel;   break;
    case mmdb::Error_ForeignFile        : ErrCode = RWBERR_ForeignFile;       break;
    case mmdb::Error_WrongEdition       : ErrCode = RWBERR_WrongEdition;      break;
    case mmdb::Error_ATOM_Unrecognized  : ErrCode = RWBERR_ATOM_Unrecognd;    break;
    case mmdb::Error_ATOM_AlreadySet    : ErrCode = RWBERR_ATOM_AlreadySet;   break;
    case mmdb::Error_ATOM_NoResidue     : ErrCode = RWBERR_ATOM_NoResidue;    break;
    case mmdb::Error_ATOM_Unmatch       : ErrCode = RWBERR_ATOM_Unmatch;      break;
    case mmdb::Error_NotACIFFile        : ErrCode = RWBERR_NotACIFFile;       break;
    case mmdb::Error_UnrecognCIFItems   : ErrCode = RWBERR_UnrecognCIFItems;  break;
    case mmdb::Error_MissingCIFField    : ErrCode = RWBERR_MissingCIFField;   break;
    case mmdb::Error_EmptyCIFLoop       : ErrCode = RWBERR_EmptyCIFLoop;      break;
    case mmdb::Error_UnexpEndOfCIF      : ErrCode = RWBERR_UnexpEndOfCIF;     break;
    case mmdb::Error_MissgCIFLoopField  : ErrCode = RWBERR_MissgCIFLoopField; break;
    case mmdb::Error_NotACIFStructure   : ErrCode = RWBERR_NotACIFStructure;  break;
    case mmdb::Error_NotACIFLoop        : ErrCode = RWBERR_NotACIFLoop;       break;
    case mmdb::Error_UnrecognizedReal   : ErrCode = RWBERR_WrongReal;         break;

    case mmdb::Error_Ok                 : ErrCode = RWBERR_Ok;                break;
    case mmdb::Error_WrongChainID       : ErrCode = RWBERR_WrongChainID;      break;
    case mmdb::Error_WrongEntryID       : ErrCode = RWBERR_WrongEntryID;      break;
    case mmdb::Error_SEQRES_serNum      : ErrCode = RWBERR_SEQRES_serNum;     break;
    case mmdb::Error_SEQRES_numRes      : ErrCode = RWBERR_SEQRES_numRes;     break;
    case mmdb::Error_SEQRES_extraRes    : ErrCode = RWBERR_SEQRES_exraRes;    break;
    case mmdb::Error_NCSM_Unrecognized  : ErrCode = RWBERR_NCSM_Unrecogn;     break;
    case mmdb::Error_NCSM_AlreadySet    : ErrCode = RWBERR_NCSM_AlreadySet;   break;
    case mmdb::Error_NCSM_WrongSerial   : ErrCode = RWBERR_NCSM_WrongSerial;  break;
    case mmdb::Error_NCSM_UnmatchIG     : ErrCode = RWBERR_NCSM_UnmatchIG;    break;
    case mmdb::Error_NoModel            : ErrCode = RWBERR_NoModel;           break;
    case mmdb::Error_NoSheetID          : ErrCode = RWBERR_NoSheetID;         break;
    case mmdb::Error_WrongSheetID       : ErrCode = RWBERR_WrongSheetID;      break;
    case mmdb::Error_WrongStrandNo      : ErrCode = RWBERR_WrongStrandNo;     break;
    case mmdb::Error_WrongNumberOfStrands : ErrCode = RWBERR_WrongNofStrands; break;
    case mmdb::Error_WrongSheetOrder    : ErrCode = RWBERR_WrongSheetOrder;   break;
    case mmdb::Error_HBondInconsistency : ErrCode = RWBERR_HBondInconsis;     break;
    case mmdb::Error_EmptyResidueName   : ErrCode = RWBERR_EmptyResidueName;  break;
    case mmdb::Error_DuplicateSeqNum    : ErrCode = RWBERR_DuplicateSeqNum;   break;
    case mmdb::Error_NoLogicalName      : ErrCode = RWBERR_NoLogicalName;     break;
    case mmdb::Error_GeneralError1      : ErrCode = RWBERR_GeneralError1;     break;

    default : ;
  }


}

bool Channel::EndOfFile()  {
int nA;
  if (MMDBManager)  {
    nA = MMDBManager->GetNumberOfAtoms();
    if (fPos>nA)  {
      fPos = nA+1;
      return true;
    }
  } else
    return  true;
  return false;
}

mmdb::PAtom * Channel::GetAtomArray()  {
  if (MMDBManager)  return MMDBManager->GetAtomArray();
              else  return NULL;
}

mmdb::PAtom Channel::GetAtomI ( int index )  {
// returns index-th atom, as counted from the
// top of file
  if (MMDBManager)  return MMDBManager->GetAtomI ( index );
              else  return NULL;
}

mmdb::PCryst Channel::GetCryst()  {
  if (MMDBManager)  return MMDBManager->GetCrystData();
              else  return NULL;
}

bool Channel::areCrystMatrices()  {
  if (MMDBManager)  return MMDBManager->isTransfMatrix();
              else  return false;
}

void  Channel::Frac2Orth (
                mmdb::realtype x,    mmdb::realtype y,    mmdb::realtype z,
                mmdb::realtype & xx, mmdb::realtype & yy, mmdb::realtype & zz )  {
  if (MMDBManager)
    MMDBManager->Frac2Orth ( x,y,z,xx,yy,zz );
  else  {
    xx = x;
    yy = y;
    zz = z;
  }
}

void  Channel::Orth2Frac (
                mmdb::realtype x,    mmdb::realtype y,    mmdb::realtype z,
                mmdb::realtype & xx, mmdb::realtype & yy, mmdb::realtype & zz )  {
  if (MMDBManager)
    MMDBManager->Orth2Frac ( x,y,z,xx,yy,zz );
  else  {
    xx = x;
    yy = y;
    zz = z;
  }
}

void  Channel::Cryst2Orth ( mmdb::rvector U )  {
  if (MMDBManager)
    MMDBManager->GetCrystData()->Cryst2Orth ( U );
}

void  Channel::Orth2Cryst ( mmdb::rvector U )  {
  if (MMDBManager)
    MMDBManager->GetCrystData()->Orth2Cryst ( U );
}


int  Channel::PutCell ( mmdb::realtype cell_a,
                         mmdb::realtype cell_b,
                         mmdb::realtype cell_c,
                         mmdb::realtype cell_alpha,
                         mmdb::realtype cell_beta,
                         mmdb::realtype cell_gamma,
                         int      OrthCode )  {

  if (MMDBManager)  {
    mmdb::PCryst cryst = MMDBManager->GetCrystData();

    cryst->PutCell ( cell_a,cell_b,cell_c,
                     cell_alpha,cell_beta,cell_gamma,
                     OrthCode );

    if ((cell_a!=0.0) || (OrthCode>0))  {
      if (cryst->CellCheck & mmdb::CCHK_Disagreement)
        return RWBERR_Disagreement;
      if (cryst->CellCheck & mmdb::CCHK_NoOrthCode)
        return RWBERR_NoOrthCode;
      if (cryst->CellCheck & mmdb::CCHK_Unchecked)
        return RWBERR_NoCheck;
    }

    return RWBERR_Ok;

  } else

    return RWBERR_NoFile;

}


int  Channel::SetCell ( mmdb::realtype cell_a,
                         mmdb::realtype cell_b,
                         mmdb::realtype cell_c,
                         mmdb::realtype cell_alpha,
                         mmdb::realtype cell_beta,
                         mmdb::realtype cell_gamma,
                         int      OrthCode )  {

  if (MMDBManager)  {
    mmdb::PCryst cryst = MMDBManager->GetCrystData();

    cryst->SetCell ( cell_a,cell_b,cell_c,
                     cell_alpha,cell_beta,cell_gamma,
                     OrthCode );

    if (cryst->CellCheck & mmdb::CCHK_Disagreement)
      return RWBERR_Disagreement;
    if (cryst->CellCheck & mmdb::CCHK_NoOrthCode)
      return RWBERR_NoOrthCode;
    if (cryst->CellCheck & mmdb::CCHK_Unchecked)
      return RWBERR_NoCheck;

    return RWBERR_Ok;

  } else

    return RWBERR_NoFile;

}


int  Channel::SetSpGroup ( mmdb::pstr spGroup )  {
  if (MMDBManager)  {
    MMDBManager->SetSpaceGroup(spGroup);
    return RWBERR_Ok;
  } else
    return RWBERR_NoFile;
}


int  Channel::GetSpGroup ( mmdb::pstr spGroup )  {
  if (MMDBManager)  {
    mmdb::PCryst cryst = MMDBManager->GetCrystData();
    if (cryst->WhatIsSet & mmdb::CSET_SpaceGroup)
          strcpy ( spGroup,cryst->spaceGroup );
    else  strcpy ( spGroup," " );
    return RWBERR_Ok;
  } else
    return RWBERR_NoFile;
}


int  Channel::GetCell ( mmdb::realtype & cell_a,
                         mmdb::realtype & cell_b,
                         mmdb::realtype & cell_c,
                         mmdb::realtype & cell_alpha,
                         mmdb::realtype & cell_beta,
                         mmdb::realtype & cell_gamma,
                         mmdb::realtype & cell_v,
                         int      & OrthCode )  {

  if (MMDBManager)  {
    mmdb::PCryst cryst = MMDBManager->GetCrystData();
    cell_a     = cryst->a;
    cell_b     = cryst->b;
    cell_c     = cryst->c;
    cell_alpha = cryst->alpha;
    cell_beta  = cryst->beta;
    cell_gamma = cryst->gamma;
    cell_v     = cryst->Vol;
    OrthCode   = cryst->NCode;
    if (!(cryst->WhatIsSet & mmdb::CSET_CellParams))
      return RWBERR_NoCellParams;
    if (!(cryst->WhatIsSet & mmdb::CSET_Transforms))
      return RWBERR_NoCheck;
//    if (MMDBManager->Cryst.CellCheck & mmdb::CCHK_NoOrthCode)
//      return RWBERR_NoOrthCode;

    return RWBERR_Ok;

  } else

    return RWBERR_NoFile;

}

int Channel::GetRCell ( mmdb::realtype & cell_as,
                         mmdb::realtype & cell_bs,
                         mmdb::realtype & cell_cs,
                         mmdb::realtype & cell_alphas,
                         mmdb::realtype & cell_betas,
                         mmdb::realtype & cell_gammas,
                         mmdb::realtype & cell_vs )  {
  if (MMDBManager)  {
    mmdb::PCryst cryst = MMDBManager->GetCrystData();
    cryst->GetRCell ( cell_as,cell_bs,cell_cs,
                      cell_alphas,cell_betas,cell_gammas,
                      cell_vs );
    if (!(cryst->WhatIsSet & mmdb::CSET_CellParams))
      return RWBERR_NoCellParams;
    if (!(cryst->WhatIsSet & mmdb::CSET_Transforms))
      return RWBERR_NoCheck;
    return RWBERR_Ok;
  } else
    return RWBERR_NoFile;
}

void Channel::GetInputBuffer ( mmdb::pstr Line, int & count )  {
  if (MMDBManager)
    MMDBManager->GetInputBuffer ( Line,count );
  else  {
    strcpy ( Line,"" );
    count = -1;
  }
}



//  ========================  static data  ===========================

static int        nChannels;    // number of channels in processing
static PChannel * channel;      // array of channels in processing

static bool       FAutoSer;     // flag to automatically generate
                                // serial numbers at reading PDB files
static bool       FReadCoords;  // flag to read coordinates; if set to
                                // false, only the header of PDB file
                                // is read
static bool       FSimRWBROOK;  // flag to simulate old RWBROOK printout
                                // as closely as possible

static char       LastFunc[80]; // name of the last called function
static int        LastUnit;     // number of the last unit called
static int        LastRC;       // last return code
static int        LastSer;      // last serial number kept for
                                // certain warnings


//  ========================  RWBrook API  ===========================


FORTRAN_SUBR ( MMDB_F_INIT, mmdb_f_init,(),(),() )  {
  mmdb::InitMatType();
  nChannels   = 0;
  channel     = NULL;
  strcpy ( LastFunc,"MMDB_F_Init" );
  LastUnit    = -1;
  LastRC      = 0;
  LastSer     = 0;
  FAutoSer    = false;
  FReadCoords = true;
  FSimRWBROOK = false;
}


FORTRAN_SUBR ( MMDB_F_QUIT, mmdb_f_quit,(),(),() )  {
int i;
  for (i=0;i<nChannels;i++)
    if (channel[i])  delete channel[i];
  if (channel) delete[] channel;
  channel   = NULL;
  nChannels = 0;
  strcpy ( LastFunc,"MMDB_F_Quit" );
  LastUnit  = -1;
  LastRC    = 0;
  LastSer   = 0;
  FAutoSer  = false;
}


FORTRAN_SUBR ( AUTOSERIALS, autoserials,
               ( int * iOnOff ),
               ( int * iOnOff ),
               ( int * iOnOff ) )  {
  FAutoSer = (*iOnOff!=0);
}


FORTRAN_SUBR ( SETREADCOORDS,setreadcoords,
               ( int * iOnOff ),
               ( int * iOnOff ),
               ( int * iOnOff ) )  {
  FReadCoords = (*iOnOff!=0);
}


FORTRAN_SUBR ( SIMRWBROOK,simrwbrook,
               ( int * iOnOff ),
               ( int * iOnOff ),
               ( int * iOnOff ) )  {
  FSimRWBROOK = (*iOnOff!=0);
}


int GetChannel ( int iUnit )  {
//   Returns serial number of the channle associated with
// unit iUnit.
//   If the channel is not found, returns -1
int i;
  for (i=0;i<nChannels;i++)
    if (channel[i])  {
      if (channel[i]->nUnit==iUnit)
        return i;
    }
  return -1;
}


int MakeChannel ( int iUnit )  {
//   If iUnit-th unit already exists, it is
// reinitialized. Otherwise the function looks
// for a not used channel, and if there is one,
// associates the new iUnit-th unit with it.
// If there is no unused channels, the new one
// is created and the new iUnit-th unit is
// associated with it.
//   Returns serial number of the channel
// associated with the newly reinitialized
// or created unit.
int        i,m;
PChannel * channel1;

  m = GetChannel ( iUnit );

  if (m>=0)  {  // such channel already exists
    channel[m]->Dispose();  // clear it first
    channel[m]->Init();     // reinitialize it
    channel[m]->nUnit = iUnit;
    return m;
  }

  for (i=0;i<nChannels;i++)  // look for free channel
    if (!channel[i])  {
      m = i;  // found!
      break;
    }

  if (m<0)  {  // no free channel
    // create new channel place
    channel1 = new PChannel[nChannels+1];
    for (i=0;i<nChannels;i++)
      channel1[i] = channel[i];
    if (channel) delete[] channel;
    channel = channel1;
    m = nChannels;
    nChannels++;  // increase number of channels
  }

  channel[m] = new Channel();  // create new channel
  channel[m]->nUnit = iUnit;

  return m;

}

FORTRAN_SUBR ( MMDB_F_OPEN, mmdb_f_open,
               (    // lengths-at-end list
                mmdb::machine::fpstr FName,      // file name
                mmdb::machine::fpstr RWStat,     // "INPUT" or "OUTPUT"
                mmdb::machine::fpstr FType,      // "PDB", "CIF", "BIN" or " "
                int * iUnit,      // channel number
                int * iRet,       // returns error code
                int   FName_len,  // fortran-hidden length of FName
                int   RWStat_len, // fortran-hidden length of RWStat
                int   FType_len   // fortran-hidden length of FType
               ), ( // lengths-in-structure list
                mmdb::machine::fpstr FName,  mmdb::machine::fpstr RWStat, mmdb::machine::fpstr FType,
                int * iUnit,  int * iRet
               ), ( // lengths-follow list
                mmdb::machine::fpstr FName,   int FName_len,
                mmdb::machine::fpstr RWStat,  int RWStat_len,
                mmdb::machine::fpstr FType,   int FType_len,
                int * iUnit,   int * iRet
               ) )  {

UNUSED_ARGUMENT(RWStat_len);
UNUSED_ARGUMENT(FType_len);

int k;
char  L[500];

#ifdef WIN32
 mmdb::GetStrTerWin32File ( L,FTN_STR(FName),0,sizeof(L),FTN_LEN(FName) );
#else
 mmdb::GetStrTer ( L,FTN_STR(FName),0,sizeof(L),FTN_LEN(FName) );
#endif

  strcpy ( LastFunc,"MMDB_F_Open" );
  LastUnit = *iUnit;

  if (*iUnit==0)  {  // generate unit number
    *iUnit = 1;
    do {
      k = GetChannel ( *iUnit );
      if (k>=0)  *iUnit = *iUnit+1;
    } while (k>=0);
  }

  // create channel
  k = MakeChannel ( *iUnit );

  if (k>=0)  {

    if (FTN_STR(RWStat)[0]=='I')  {
      channel[k]->nRead       = 0;
      channel[k]->FAutoSer    = FAutoSer;
      channel[k]->FReadCoords = FReadCoords;
      channel[k]->FSimRWBROOK = FSimRWBROOK;
    } else
      channel[k]->nRead = 1;

    // store file name
    channel[k]->SetFileName ( L,sizeof(L) );

    // store unit type
    channel[k]->SetFileType ( FTN_STR(FType) );
    channel[k]->IdentifyFile( L );

    if (FSimRWBROOK)  {
      switch (channel[k]->nType)  {
        default : printf ( "  unknown-format" );  break;
        case  0 : printf ( "  PDB"   );           break;
        case  1 : printf ( "  mmCIF" );           break;
        case  2 : printf ( "  MMDB BINARY" );
      }
      printf ( " file is being opened on unit %i",*iUnit );
      if (FTN_STR(RWStat)[0]=='I')  printf ( " for INPUT.\n\n" );
                              else  printf ( " for OUTPUT.\n\n" );
    }

    if (FTN_STR(RWStat)[0]=='I')  {
      channel[k]->Read();
      *iRet = channel[k]->ErrCode;
    } else  {
      channel[k]->MakeCoordStructure();
      channel[k]->fPos = 1;
      *iRet = RWBERR_Ok;
    }

  } else
    *iRet = RWBERR_NoChannel;

  LastRC = *iRet;

}

FORTRAN_SUBR ( MMDB_F_OPENL, mmdb_f_openl,
               (    // lengths-at-end list
                mmdb::machine::fpstr LName,      // logical name
                mmdb::machine::fpstr RWStat,     // "INPUT" or "OUTPUT"
                mmdb::machine::fpstr FType,      // "PDB", "CIF", "BIN" or " "
                int * iUnit,      // channel number
                int * iRet,       // returns error code
                int   LName_len,  // fortran-hidden length of LName
                int   RWStat_len, // fortran-hidden length of RWStat
                int   FType_len   // fortran-hidden length of FType
               ), ( // lengths-in-structure list
                mmdb::machine::fpstr LName,  mmdb::machine::fpstr RWStat, mmdb::machine::fpstr FType,
                int * iUnit,  int * iRet
               ), ( // lengths-follow list
                mmdb::machine::fpstr LName,   int LName_len,
                mmdb::machine::fpstr RWStat,  int RWStat_len,
                mmdb::machine::fpstr FType,   int FType_len,
                int * iUnit,   int * iRet
               ) )  {
char        L[200];
mmdb::pstr  S;
char_struct(FName)

  strcpy ( LastFunc,"MMDB_F_Openl" );

  mmdb::GetStrTer ( L,FTN_STR(LName),0,sizeof(L),FTN_LEN(LName) );

  S = getenv ( L );

  if (S)  {

    fill_char_struct(FName,S)

  } else if (FTN_STR(RWStat)[0]=='O') {

    // The user may not have assigned a logical
    // for output, so that the program should write file "XYZOUT". This
    // is allowed as a convenience when user is not really interested
    // in output file.
    fill_char_struct(FName,L)

  } else {
    *iRet = RWBERR_NoLogicalName;
    return;
  }

  printf ( "\n  Logical name: %s  File name: %s\n",L,FName );

  FORTRAN_CALL ( MMDB_F_OPEN, mmdb_f_open,
                   ( FName,RWStat,FType,iUnit,iRet,
                     FName_len,RWStat_len,FType_len ),
                   ( &FName,RWStat,FType,iUnit,iRet ),
                   ( FName,FName_len,RWStat,RWStat_len,
                     FType,FType_len,iUnit,iRet ) );

}

FORTRAN_SUBR ( MMDB_F_COPY, mmdb_f_copy,
               (    // lengths-at-end list
                int * iUnit1,    // destination unit
                int * iUnit2,    // source unit
                int * copyKey,   // copy key:
                                 //  = 1  copy all
                                 //  = 2  copy all except coordinates
                                 //  = 3  copy title section only
                                 //  = 4  copy crystallographic
                                 //       section only
                                 //  = 5  copy coordinate section only
                                 // any other value does not do anything
                int * iRet       // return code:
                                 //   =0 if success
                                 //   =RWBERR_NoChannel if a unit
                                 //                   does not exist
               ), ( // lengths-in-structure list
                int * iUnit1,  int * iUnit2,
                int * copyKey, int * iRet
               ), ( // lengths-follow list
                int * iUnit1,  int * iUnit2,
                int * copyKey, int * iRet
               ) )  {
int             k1,k2;
mmdb::COPY_MASK copyMask;

  strcpy ( LastFunc,"MMDB_F_Copy" );

  LastUnit = *iUnit1;
  k1 = GetChannel ( LastUnit );

  if (k1>=0)  {
    if (channel[k1]->MMDBManager)  {
      LastUnit = *iUnit2;
      k2 = GetChannel ( LastUnit );
      if (k2>=0)  {
        if (channel[k2]->MMDBManager)  {
          switch (*copyKey)  {
            case 1  :  copyMask = mmdb::MMDBFCM_All;    break;
            case 2  :  copyMask = mmdb::MMDBFCM_Top;    break;
            case 3  :  copyMask = mmdb::MMDBFCM_Title;  break;
            case 4  :  copyMask = mmdb::MMDBFCM_Cryst;  break;
            case 5  :  copyMask = mmdb::MMDBFCM_Coord;  break;
            default :  copyMask = mmdb::MMDBFCM_None;
          }
          channel[k1]->MMDBManager->Copy ( channel[k2]->MMDBManager,
                                           copyMask );
          *iRet = RWBERR_Ok;
       } else
          *iRet = RWBERR_NoFile;
      } else
        *iRet = RWBERR_NoChannel;
    } else
      *iRet = RWBERR_NoFile;
  } else
    *iRet = RWBERR_NoChannel;

  LastRC = *iRet;

}


FORTRAN_SUBR ( MMDB_F_DELETE, mmdb_f_delete,
               (    // lengths-at-end list
                int * iUnit,     // unit number; *iUnit<=0 means
                                 // "the last mentioned unit"
                int * delKey,    // delete key:
                                 //  = 1  delete all
                                 //  = 2  delete all except coordinates
                                 //  = 3  delete title section only
                                 //  = 4  delete crystallographic
                                 //       section only
                                 //  = 5  delete coordinate section only
                                 // any other value does not do anything
                int * iRet       // return code:
                                 //   =0 if success
                                 //   =RWBERR_NoChannel if a unit
                                 //                   does not exist
                                 //   =RWBERR_NoFile    if a unit
                                 //                   was not opened
               ), ( // lengths-in-structure list
                int * iUnit, int * delKey, int * iRet
               ), ( // lengths-follow list
                int * iUnit, int * delKey, int * iRet
               ) )  {
int        k;
mmdb::word delMask;

  strcpy ( LastFunc,"MMDB_F_Delete" );

  if (*iUnit>0)
    LastUnit = *iUnit;
  k = GetChannel ( LastUnit );

  if (k>=0)  {
    if (channel[k]->MMDBManager)  {
      switch (*delKey)  {
        case 1  :  delMask = mmdb::MMDBFCM_All;    break;
        case 2  :  delMask = mmdb::MMDBFCM_Top;    break;
        case 3  :  delMask = mmdb::MMDBFCM_Title;  break;
        case 4  :  delMask = mmdb::MMDBFCM_Cryst;  break;
        case 5  :  delMask = mmdb::MMDBFCM_Coord;  break;
        default :  delMask = 0x0000;
      }
      channel[k]->MMDBManager->Delete ( delMask );
      *iRet = RWBERR_Ok;
    } else
      *iRet = RWBERR_NoFile;
  } else
    *iRet = RWBERR_NoChannel;

  LastRC = *iRet;

}


FORTRAN_SUBR ( MMDB_F_SETTYPE, mmdb_f_settype,
               (    // lengths-at-end list
                int * iUnit,     // unit number
                mmdb::machine::fpstr FType,     // "PDB", "CIF", "BIN" or " "
                mmdb::machine::fpstr RWStat,    // "INPUT" or "OUTPUT"
                int * iRet,      // returns -1 if unit not found,
                                 // otherwise 0
                int   FType_len, // fortran-hidden length of FType
                int   RWStat_len // fortran-hidden length of RWStat
               ), ( // lengths-in-structure list
                int * iUnit,  mmdb::machine::fpstr FType,
                mmdb::machine::fpstr RWStat, int * iRet
               ), ( // length-follow list
                int * iUnit,
                mmdb::machine::fpstr FType,   int FType_len,
                mmdb::machine::fpstr RWStat,  int RWStat_len,
                int * iRet
               ) )  {
UNUSED_ARGUMENT(FType_len);
UNUSED_ARGUMENT(RWStat_len);

int k;

  strcpy ( LastFunc,"MMDB_F_SetType" );
  if (*iUnit>0)
    LastUnit = *iUnit;

  k = GetChannel ( LastUnit );

  if (k>=0)  {
    // store unit type
    channel[k]->SetFileType ( FTN_STR(FType) );
    // store unit mode
    if (FTN_STR(RWStat)[0]=='I')  channel[k]->nRead = 0;
                            else  channel[k]->nRead = 1;
    *iRet = RWBERR_Ok;
  } else
    *iRet = RWBERR_NoChannel;

  LastRC = *iRet;

}



FORTRAN_SUBR ( MMDB_F_SETNAME, mmdb_f_setname,
               (    // lengths-at-end list
                int * iUnit,    // unit number
                mmdb::machine::fpstr FName,    // file name
                int * iRet,     // returns -1 if unit not found,
                                // otherwise 0
                int   FName_len // fortran-hidden length of FName
               ), ( // lengths-in-structure list
                int * iUnit, mmdb::machine::fpstr FName, int * iRet
               ), ( // lengths-follow list
                int * iUnit,
                mmdb::machine::fpstr FName, int FName_len,
                int * iRet
               ) )  {

int k;

  strcpy ( LastFunc,"MMDB_F_SetName" );
  if (*iUnit>0)
    LastUnit = *iUnit;

  k = GetChannel ( LastUnit );

  if (k<0)
    *iRet = RWBERR_NoChannel;
  else  {
    // store file name
    channel[k]->SetFileName ( FTN_STR(FName),FTN_LEN(FName) );
    *iRet = RWBERR_Ok;
  }

  LastRC = *iRet;

}


FORTRAN_SUBR ( MMDB_F_WRITE,  mmdb_f_write,
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ) )  {
int k;

  strcpy ( LastFunc,"MMDB_F_Write" );
  if (*iUnit>0)
    LastUnit = *iUnit;

  k = GetChannel ( LastUnit );

  if (k<0)
    *iRet = RWBERR_NoChannel;
  else  {
    channel[k]->Write();
    *iRet = channel[k]->ErrCode;
  }

  LastRC = *iRet;

}


FORTRAN_SUBR ( MMDB_F_CLOSE, mmdb_f_close,
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ) )  {
int k;

  strcpy ( LastFunc,"MMDB_F_Close" );
  if (*iUnit>0)
    LastUnit = *iUnit;

  k = GetChannel ( LastUnit );

  if (k<0)
    *iRet = RWBERR_NoChannel;
  else if (channel[k]->nRead==1)  {
    channel[k]->Write();
    *iRet = channel[k]->ErrCode;
    if (!(*iRet))  {
      delete channel[k];
      channel[k] = NULL;
    }
  } else  {
    delete channel[k];
    channel[k] = NULL;
    *iRet = RWBERR_Ok;
  }

  LastRC = *iRet;

}



FORTRAN_SUBR ( MMDB_F_ADVANCE, mmdb_f_advance,
               (   // lengths-at-end list
                int * iUnit, // unit number
                int * iOut,  // output echo file
                int * iTer,  // FLAG =1, return iRet=1 if 'ter' card found
                             //      =0, do not return on 'ter' card
                int * iRet   // =0  if normal return
                             // =1  if return on 'ter' card (iTer=1)
                             // =2  if return on end of file
                             // =3  if return on 'hetatm' card
                             // =RWBERR_NoChannel if unit does not exist
                             // =RWBERR_NoAdvance if pointer was not
                             //                   advanced
               ), ( // lengths-in-structure list
                int * iUnit, int * iOut, int * iTer, int * iRet
               ), ( // lengths-follow list
                int * iUnit, int * iOut, int * iTer, int * iRet
               ) )  {

UNUSED_ARGUMENT(iOut);

int    k;
mmdb::PAtom atom;

  strcpy ( LastFunc,"mmdb_f_advance" );
  LastUnit = *iUnit;

  k = GetChannel ( *iUnit );

  if (k<0)

    *iRet = RWBERR_NoChannel;

  else if (channel[k]->nRead==0)  {

    // in the input file, try to get pointer on the next atom

    do {
      channel[k]->fPos++;  // advance the pointer on Atom array
      if (channel[k]->EndOfFile())  {
        atom = NULL;
        break;
      }
      atom = channel[k]->GetAtomI ( channel[k]->fPos );
      if (atom)  {
        if ((atom->Ter) && (*iTer==0))  {
          // ignore 'ter' card if iTer is set to 0
          atom = NULL;
        }
      }
    } while (!atom);

    if (!atom)  *iRet = 2; // no atom found == end of file
    else if (atom->Ter)  *iRet = 1; // 'ter' card encountered
    else if (atom->Het)  *iRet = 3; // 'hetatm' card encountered
                   else  *iRet = 0; // advance ok; normal return

  } else  {

    // in the output file, just advance the pointer

    if (channel[k]->fPos==0)  {
      channel[k]->fPos++;
      *iRet = 0;
    } else  {
      atom = channel[k]->GetAtomI ( channel[k]->fPos );
      if (atom)  {
        // the previous atom was set -- advance the pointer
        channel[k]->fPos++;
        *iRet = 0;
      } else
        // no atom was set; make no advancement
        *iRet = RWBERR_NoAdvance;
    }

  }

  LastRC = *iRet;

}



FORTRAN_SUBR ( MMDB_F_REWD, mmdb_f_rewd,
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ) )  {
int k;

  strcpy ( LastFunc,"MMDB_F_Rewd" );
  LastUnit = *iUnit;

  k = GetChannel ( *iUnit );
  if (k>=0)  {
    channel[k]->fPos = 0;
    if (channel[k]->nRead!=0)  *iRet = RWBWAR_RewOutput;
                         else  *iRet = RWBERR_Ok;
  } else
    *iRet = RWBERR_NoChannel;

  LastRC = *iRet;

}



FORTRAN_SUBR ( MMDB_F_BKSP, mmdb_f_bksp,
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ) )  {
int k;

  strcpy ( LastFunc,"MMDB_F_BkSp" );
  LastUnit = *iUnit;

  k = GetChannel ( *iUnit );
  if (k>=0)  {
    *iRet = RWBERR_Ok;
    if (channel[k]->fPos==0)  *iRet |= RWBWAR_FileTop;
                        else  channel[k]->fPos--;
    if (channel[k]->nRead!=0) *iRet |= RWBWAR_RewOutput;
  } else
    *iRet = RWBERR_NoChannel;

  LastRC = *iRet;

}



FORTRAN_SUBR ( MMDB_F_SEEK, mmdb_f_seek,
               (    // lengths-at-end list
                int * iUnit,  // unit number
                int * fPos,   // position to set
                int * iRet    // return code:
                              //  0   Ok
                              //  1   'ter' card met
                              //  2   end of file
                              //  3   'hetatm' card met
                              //  <0 error:
                              //  RWBERR_NoChannel
                              //      iUnit was not
                              //      initialized
                              //  RWBERR_EmptyPointer
                              //      fPos-th position
               ), ( // lengths-in-structure list
                int * iUnit, int * fPos, int * iRet
               ), ( // lengths-follow list
                int * iUnit, int * fPos, int * iRet
               ) )  {
int    k;
mmdb::PAtom atom;

  strcpy ( LastFunc,"MMDB_F_Seek" );
  LastUnit = *iUnit;

  k = GetChannel ( *iUnit );

  if (k<0)

    *iRet = RWBERR_NoChannel;

  else  {

    // set the pointer
    channel[k]->fPos = mmdb::IMax(0,*fPos);
    if (*fPos==0)  *iRet = RWBWAR_FileTop;
             else  *iRet = RWBERR_Ok;

    if (channel[k]->nRead==0)  {

      // in the input file, check the end-of-file state
      // and analyze the atom

      if (channel[k]->EndOfFile())  *iRet = 2;

      atom = channel[k]->GetAtomI ( channel[k]->fPos );

      if (!atom)          *iRet = RWBERR_EmptyPointer; // empty place
      else if (atom->Ter) *iRet = 1;  // 'ter' card encountered
      else if (atom->Het) *iRet = 3;  // 'hetatm' card encountered

    }

    // in the output file, there is nothing to do

  }

  LastRC = *iRet;

}


void  Make_AN_ID_IZ ( mmdb::PAtom atom, mmdb::pstr AtNam, int AtNam_L,
                      mmdb::pstr ID, int ID_L, int * IZ, int * iRet )  {
char chrg[10];
int  i,k;

  if (atom->Ter)  {

    mmdb::strcpy_ns ( AtNam,mmdb::pstr(" "),AtNam_L );
    mmdb::strcpy_ns ( ID   ,mmdb::pstr(" "),ID_L    );
    *IZ = 7;

  } else  {

    if (atom->name[0]==' ')  mmdb::strcpy_ns ( AtNam,&(atom->name[1]),4 );
                       else  mmdb::strcpy_ns ( AtNam,atom->name,4 );

    // first try to identify the atom with the element name
    mmdb::strcpy_ns ( ID,atom->element,ID_L );  // not more than ID_L symbols
                                // from element until but not including
                                // the terminated null are copied into
                                // ID, and the latter is padded with
                                // spaces up to the length of ID_L

    if (ID_L>3)  {  // if length permits, add ID with atom charge
                    // (always 2 symbols).
      atom->GetAtomCharge(chrg);
      ID[2] = chrg[0];
      ID[3] = chrg[1];
    }

    k = 0;
    while ((k<mmdb::nElementNames) &&
           ((atom->element[0]!=mmdb::ElementName[k][0]) ||
            (atom->element[1]!=mmdb::ElementName[k][1])))  k++;

    if (k>=mmdb::nElementNames)  {

      // no match for atom ID -- make sure to set it blank
      mmdb::strcpy_ns ( ID,mmdb::pstr(" "),ID_L );

      //  try to identify the atom using the atom name
      k = 0;
      while ((k<mmdb::nElementNames) &&
             ((atom->name[0]!=mmdb::ElementName[k][0]) ||
              (atom->name[1]!=mmdb::ElementName[k][1])))  k++;

      // try to identify a heteroatom
      i = 0;
      while ((i<mmdb::nHydAtomNames) && (k>=mmdb::nElementNames))  {
        if ((atom->name[0]==mmdb::HydAtomName[i][0]) &&
            (atom->name[1]==mmdb::HydAtomName[i][1]))
          k = 0;
        i++;
      }

      if (k>=mmdb::nElementNames)  {
        // unknown or ambiguous formfactor
        k = -1;
        if ((atom->name[0]==' ') &&
            (atom->name[1]=='A'))  k = 6;
        if (k==-1)  *iRet |= RWBWAR_UnkFormFactor;
              else  *iRet |= RWBWAR_AmbFormFactor;
      }

    }

    *IZ = k+1;
    if (*IZ==0)
      mmdb::strcpy_ns ( ID,mmdb::pstr(" "),ID_L );
    else  {
      if (ID_L>3)  {
        if (ID[0]==' ')  {
          if ((AtNam[2]=='+') ||
              (AtNam[2]=='-'))  {
            ID[2] = AtNam[2];
            ID[3] = AtNam[3];
          }
        } else if ((ID[2]!='+') && (ID[2]!='-'))  {
          ID[2] = ' ';
          ID[3] = ' ';
        }
      }
      mmdb::strcpy_ns ( ID,mmdb::ElementName[k],mmdb::IMin(2,ID_L) );
    }

  }

}



FORTRAN_SUBR ( MMDB_F_ATOM,  mmdb_f_atom,
           (    // lengths-at-end list
            int * iUnit,    // unit number
            int * iSer,     // atom serial number
            mmdb::machine::fpstr AtNam,    // atom name (left justified)
            mmdb::machine::fpstr ResNam,   // residue name
            mmdb::machine::fpstr ChnNam,   // chain name
            int * iResN,    // residue number as an integer
            mmdb::machine::fpstr ResNo,    // residue number as character (input only)
            mmdb::machine::fpstr InsCod,   // the insertion code
            mmdb::machine::fpstr AltCod,   // the alternate conformation code
            mmdb::machine::fpstr segID,    // segment ID
            int * IZ,       // atomic number (input only, returned as
                            // 7 from ambiguous atoms)
            mmdb::machine::fpstr ID,       // atomic ID related to atomic number
                            // (element symbol right justified), plus
                            // the ionic state +2, +3 etc..
                            //
            int * iRet,     // returns
                            //  RWBERR_NoChannel     if iUnit was not
                            //                       initialized
                            //  RWBERR_EmptyPointer  if atom was not
                            //                       advanced
                            //  RWBERR_Error1        internal error #1
                            //  RWBERR_Error2        internal error #2
                            //  RWBERR_Error3        internal error #3
                            //
                            //  >=0 : success, warning flags:
                            //  RWBWAR_WrongSerial   if serial number
                            //               differs from the position
                            //               number in the file
                            //  RWBWAR_UnkFormFactor unknown formfactor
                            //  RWBWAR_AmbFormFactor ambiguous formfactor
                            //
            int AtNam_len,  // fortran-hidden length of AtNam
            int ResNam_len, // fortran-hidden length of ResNam
            int ChnNam_len, // fortran-hidden length of ChnNam
            int ResNo_len,  // fortran-hidden length of ResNo
            int InsCod_len, // fortran-hidden length of InsCod
            int AltCod_len, // fortran-hidden length of AltCod
            int segID_len,  // fortran-hidden length of SegID
            int ID_len      // fortran-hidden length of ID
           ), ( // lengths-in-structure list
            int * iUnit,  int * iSer,  mmdb::machine::fpstr AtNam,  mmdb::machine::fpstr ResNam,
            mmdb::machine::fpstr ChnNam, int * iResN, mmdb::machine::fpstr ResNo,  mmdb::machine::fpstr InsCod,
            mmdb::machine::fpstr AltCod, mmdb::machine::fpstr segID, int * IZ,     mmdb::machine::fpstr ID,
            int * iRet
           ), ( // lengths-follow list
            int * iUnit,  int * iSer,
            mmdb::machine::fpstr AtNam,  int   AtNam_len,
            mmdb::machine::fpstr ResNam, int   ResNam_len,
            mmdb::machine::fpstr ChnNam, int   ChnNam_len,
            int * iResN,
            mmdb::machine::fpstr ResNo,  int   ResNo_len,
            mmdb::machine::fpstr InsCod, int   InsCod_len,
            mmdb::machine::fpstr AltCod, int   AltCod_len,
            mmdb::machine::fpstr segID,  int   segID_len,
            int * IZ,
            mmdb::machine::fpstr ID,     int   ID_len,
            int * iRet
           ) )  {
int            k,RC;
mmdb::ChainID  chainID;
mmdb::ResName  resName;
mmdb::InsCode  insCode;
mmdb::AtomName atomName;
mmdb::AltLoc   altLoc;
mmdb::SegID    sgID;
mmdb::Element  element;
mmdb::PAtom    atom;
char           charge[10];

  strcpy ( LastFunc,"MMDB_F_Atom" );
  LastUnit = *iUnit;

  k = GetChannel ( *iUnit );
  if (k<0)  {
    *iRet  = RWBERR_NoChannel;
    LastRC = *iRet;
    return;
  }

  *iRet = RWBERR_Ok;

  if (channel[k]->nRead==0)  {

    // reading the atom characteristics

    atom = channel[k]->GetAtomI ( channel[k]->fPos );
    if (!atom)  {
      // atom position was not advanced properly
      *iRet  = RWBERR_EmptyPointer;
      LastRC = *iRet;
      return;
    }

    *iSer = atom->serNum;
    if (*iSer!=channel[k]->fPos)  *iRet |= RWBWAR_WrongSerial;
    LastSer = *iSer;
    Make_AN_ID_IZ ( atom,FTN_STR(AtNam),FTN_LEN(AtNam),
                    FTN_STR(ID),FTN_LEN(ID),IZ,iRet );
    if (atom->residue)  {
      mmdb::strcpy_ns  ( FTN_STR(ResNam),atom->residue->name,FTN_LEN(ResNam) );
      *iResN = atom->residue->seqNum;
      mmdb::PutInteger ( FTN_STR(ResNo),*iResN,mmdb::IMin(4,FTN_LEN(ResNo)) );
      mmdb::strcpy_ns  ( FTN_STR(InsCod),atom->residue->insCode,FTN_LEN(InsCod) );
      mmdb::strcpy_ns  ( &(FTN_STR(ResNo)[4]),FTN_STR(InsCod),FTN_LEN(ResNo)-4 );
      mmdb::strcpy_ns  ( FTN_STR(ChnNam),atom->GetChainID(),FTN_LEN(ChnNam) );
    } else  {
      mmdb::strcpy_ns  ( FTN_STR(ResNam),mmdb::pstr("   "),FTN_LEN(ResNam) );
      mmdb::strcpy_ns  ( FTN_STR(ChnNam),mmdb::pstr(" ")  ,FTN_LEN(ChnNam) );
      *iResN = 0;
      mmdb::strcpy_ns  ( FTN_STR(ResNo) ,mmdb::pstr("0")  ,FTN_LEN(ResNo)  );
      mmdb::strcpy_ns  ( FTN_STR(InsCod),mmdb::pstr(" ")  ,FTN_LEN(InsCod) );
    }
    mmdb::strcpy_ns ( FTN_STR(AltCod),atom->altLoc,FTN_LEN(AltCod) );
    mmdb::strcpy_ns ( FTN_STR(segID) ,atom->segID ,FTN_LEN(segID)  );

  } else  {

    // storing the atom characteristics

    if (!channel[k]->MMDBManager)  {
      *iRet  = RWBERR_Error1;   // should never happen
      LastRC = *iRet;
      return;
    }

    mmdb::GetStrTer ( chainID,FTN_STR(ChnNam),1,sizeof(chainID),FTN_LEN(ChnNam) );
    mmdb::GetStrTer ( resName,FTN_STR(ResNam),3,sizeof(resName),FTN_LEN(ResNam) );
    mmdb::GetStrTer ( insCode,FTN_STR(InsCod),1,sizeof(insCode),FTN_LEN(InsCod) );
    mmdb::GetStrTer ( altLoc ,FTN_STR(AltCod),1,sizeof(altLoc) ,FTN_LEN(AltCod) );
    mmdb::GetStrTer ( sgID   ,FTN_STR(segID) ,4,sizeof(sgID)   ,FTN_LEN(segID)  );
    element[0] = FTN_STR(ID)[0];
    element[1] = FTN_STR(ID)[1];
    element[2] = char(0);
    if (FTN_LEN(ID)>3)  {
      charge [0] = FTN_STR(ID)[2];
      charge [1] = FTN_STR(ID)[3];
      charge [2] = char(0);
    } else
      charge [0] = char(0);

    if (FTN_STR(ID)[0]==' ')  {
      atomName[0] = char(0);
//      if ((FTN_STR(AtNam)[1]=='H') ||
//          ((FTN_STR(AtNam)[1]=='D') && (FTN_STR(ID)[2]=='D')))  {
//        int i = 0;
//        while ((i<nHydAtomNames) &&
//               (FTN_STR(AtNam)[0]!=HydAtomName[i][0])) i++;
//        if (i<nHydAtomNames)
//          GetStrTer ( atomName,FTN_STR(AtNam),4,5,FTN_LEN(AtNam) );
//      }
      if ((FTN_STR(AtNam)[0]=='H') && (FTN_STR(AtNam)[3]!=' '))
        mmdb::GetStrTer ( atomName,FTN_STR(AtNam),4,5,FTN_LEN(AtNam) );
      if (!atomName[0])  {
        atomName[0] = ' ';
        mmdb::GetStrTer ( &(atomName[1]),FTN_STR(AtNam),3,4,FTN_LEN(AtNam) );
      }
    } else
      mmdb::GetStrTer ( atomName,FTN_STR(AtNam),4,5,4 );

    RC = channel[k]->MMDBManager->PutAtom ( channel[k]->fPos,*iSer,
                              atomName,resName,chainID,*iResN,
                              insCode,altLoc,sgID,element );

    if (RC)  {
      *iRet  = RWBERR_Error2;  // should never happen
      LastRC = *iRet;
      return;
    }

    mmdb::DelSpaces ( charge );
    if (charge[0])  {
      atom  = channel[k]->GetAtomI ( channel[k]->fPos );
      if (!atom)  {
        *iRet  = RWBERR_EmptyPointer; // should never be so
        LastRC = *iRet;
        return;
      }
      atom->SetCharge ( charge );
    }

    if (*iSer!=channel[k]->fPos)  {
      *iRet |= RWBWAR_WrongSerial; // this is not the right thing at all
      atom  = channel[k]->GetAtomI ( channel[k]->fPos );
      if (!atom)  {
        *iRet  = RWBERR_EmptyPointer; // should never be so
        LastRC = *iRet;
        return;
      }
      //      atom->serNum = *iSer;        // - we allow for a mess in serials
    }

    LastSer = *iSer;

  }

  LastRC = *iRet;

}



FORTRAN_SUBR ( MMDB_F_SETTER, mmdb_f_setter,
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ) )  {
int    k;
mmdb::PAtom atom;

  strcpy ( LastFunc,"MMDB_F_SetTer" );
  LastUnit = *iUnit;

  k = GetChannel ( *iUnit );
  if (k<0)  {
    *iRet  = RWBERR_NoChannel;
    LastRC = *iRet;
    return;
  }

  atom  = channel[k]->GetAtomI ( channel[k]->fPos );
  *iRet = RWBERR_Ok;

  if (!atom)  {
    *iRet  = RWBERR_EmptyPointer;  // atom position was not advanced properly
    LastRC = *iRet;
    return;
  }

  atom->Ter       = true;
  atom->WhatIsSet |= mmdb::ASET_Coordinates;

}



FORTRAN_SUBR ( MMDB_F_SETHET, mmdb_f_sethet,
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ),
               ( int * iUnit, int * iRet ) )  {
int    k;
mmdb::PAtom atom;

  strcpy ( LastFunc,"MMDB_F_SetHet" );
  LastUnit = *iUnit;

  k = GetChannel ( *iUnit );
  if (k<0)  {
    *iRet  = RWBERR_NoChannel;
    LastRC = *iRet;
    return;
  }

  atom  = channel[k]->GetAtomI ( channel[k]->fPos );
  *iRet = RWBERR_Ok;

  if (!atom)  {
    *iRet  = RWBERR_EmptyPointer;  // atom position was not advanced properly
    LastRC = *iRet;
    return;
  }

  atom->Het = true;
  atom->WhatIsSet |= mmdb::ASET_Coordinates;

}

FORTRAN_SUBR ( MMDB_F_GETHET, mmdb_f_gethet,
               ( int * iUnit, int * isHet, int * iRet ),
               ( int * iUnit, int * isHet, int * iRet ),
               ( int * iUnit, int * isHet, int * iRet ) )  {
int    k;
mmdb::PAtom atom;

  strcpy ( LastFunc,"MMDB_F_GetHet" );
  LastUnit = *iUnit;

  *isHet = 0;  //  no HETATM record

  k = GetChannel ( *iUnit );
  if (k<0)  {
    *iRet  = RWBERR_NoChannel;
    LastRC = *iRet;
    return;
  }

  atom  = channel[k]->GetAtomI ( channel[k]->fPos );
  *iRet = RWBERR_Ok;

  if (!atom)  {
    *iRet  = RWBERR_EmptyPointer;  // atom position was not advance properly
    LastRC = *iRet;
    return;
  }

  if (atom->Het)  *isHet = 1;      // HETATM

}



// ------------------------------------------------------------------

//   mmdb_f_getnofncsmates_(..) returns the number of NCS mates not
//  given in the file (iGiven=0).
//
//  Negative returns N<0  mean an error.
//
//  FORTRAN equivalent:   subroutine MMDB_F_GetNofNCSMates ( iUnit,N )
//  ~~~~~~~~~~~~~~~~~~~   integer  iUnit,N

FORTRAN_SUBR ( MMDB_F_GETNOFNCSMATES, mmdb_f_getnofncsmates,
               ( int * iUnit, int * N ),
               ( int * iUnit, int * N ),
               ( int * iUnit, int * N ) )  {
int k;

  strcpy ( LastFunc,"mmdb_f_getnofncsmates" );
  LastUnit = *iUnit;

  k = GetChannel ( *iUnit );
  if (k<0)  {
    *N     = RWBERR_NoChannel;
    LastRC = *N;
    return;
  }

  *N = channel[k]->GetNumberOfNCSMates();
  return;

}


FORTRAN_SUBR ( MMDB_F_COPYATOM, mmdb_f_copyatom,
               ( int * iUnit1, int * iUnit2, int * iRet ),
               ( int * iUnit1, int * iUnit2, int * iRet ),
               ( int * iUnit1, int * iUnit2, int * iRet ) )  {
int    k1,k2,RC;
mmdb::PAtom atom;

  strcpy ( LastFunc,"mmdb_f_copyatom" );
  LastUnit = *iUnit1;

  k1 = GetChannel ( *iUnit1 );
  if (k1<0)  {
    *iRet  = RWBERR_NoChannel;
    LastRC = *iRet;
    return;
  }

  k2 = GetChannel ( *iUnit2 );
  if (k2<0)  {
    *iRet  = RWBERR_NoChannel;
    LastRC = *iRet;
    return;
  }

  atom = channel[k1]->GetAtomI ( channel[k1]->fPos );
  *iRet = RWBERR_Ok;

  if (!atom)  {
    *iRet  = RWBERR_EmptyPointer;  // atom position was not advanced
                                   // properly
    LastRC = *iRet;
    return;
  }

  RC = channel[k2]->MMDBManager->PutAtom ( channel[k2]->fPos,atom,
             atom->serNum );
  if (RC)  {
    *iRet  = RWBERR_Error2;  // should never happen
    LastRC = *iRet;
    return;
  }

  LastSer = atom->serNum;

}


FORTRAN_SUBR ( MMDB_F_COORD, mmdb_f_coord,
               (    // lengths-at-end list
                int * iUnit,    // unit number
                mmdb::machine::fpstr XFlag,    // "F" or "O" flag for the fractional
                                // or orthogonal coordinates x,y,z
                                // for output files XFlag may also be
                                // set to "HF" or "HO", where "F" and
                                // "O" have the same meaning as before
                                // and "H" indicates that the atom
                                // should be marked as heteroatom
                mmdb::machine::fpstr BFlag ,   // "F" or "O" flag for temperature
                                // factor in fractional or orthogonal
                                // Us
                mmdb::machine::apireal * x,    // x-coordinate
                mmdb::machine::apireal * y,    // y-coordinate
                mmdb::machine::apireal * z,    // z-coordinate
                mmdb::machine::apireal * occ,  // occupancy
                mmdb::machine::apireal * BIso, // isotropic temperature factor
                mmdb::machine::apireal * U,    // array(6) of the anisotr. t-factor
                int * iRet,     // returns
                                //  RWBERR_NoChannel     if iUnit was not
                                //                       initialized
                                //  RWBERR_EmptyPointer  if atom was not
                                //                       advanced
                                //  RWBERR_NoMatrices    if transformation
                                //                       matrices are
                                //                       undefined
                                //  RWBERR_NoCoordinates if coordinates were
                                //                       not set in the atom
                                //
                                //  >=0 : success, warning flags:
                                //  RWBERR_NoOccupancy   if occupancy was
                                //                       not set in the atom
                                //  RWBERR_NoTempFactor  if temp. factor was
                                //                       not set in the atom
                                //
                int XFlag_len,  // fortran-hidden length of XFlag
                int BFlag_len   // fortran-hidden length of BFlag
               ), ( // lengths-in-structure list
                int * iUnit,   mmdb::machine::fpstr XFlag,    mmdb::machine::fpstr BFlag,
                mmdb::machine::apireal * x,   mmdb::machine::apireal * y,    mmdb::machine::apireal * z,
                mmdb::machine::apireal * occ, mmdb::machine::apireal * BIso, mmdb::machine::apireal * U,
                int * iRet
               ), ( // lengths-follow list
                int * iUnit,
                mmdb::machine::fpstr XFlag,   int XFlag_len,
                mmdb::machine::fpstr BFlag,   int BFlag_len,
                mmdb::machine::apireal * x,   mmdb::machine::apireal * y,    mmdb::machine::apireal * z,
                mmdb::machine::apireal * occ, mmdb::machine::apireal * BIso, mmdb::machine::apireal * U,
                int * iRet
               ) )  {

UNUSED_ARGUMENT(XFlag_len);
UNUSED_ARGUMENT(BFlag_len);

mmdb::realtype AU[6];
mmdb::realtype xx,yy,zz;
int            k,i,m;
mmdb::PAtom    atom;

  strcpy ( LastFunc,"MMDB_F_Coord" );
  LastUnit = *iUnit;

  k = GetChannel ( *iUnit );
  if (k<0)  {
    *iRet  = RWBERR_NoChannel;
    LastRC = *iRet;
    return;
  }

  atom  = channel[k]->GetAtomI ( channel[k]->fPos );
  *iRet = RWBERR_Ok;

  if (!atom)  {
    *iRet  = RWBERR_EmptyPointer;  // atom position was not advanced properly
    LastRC = *iRet;
    return;
  }

  if ((FTN_STR(XFlag)[0]=='H') ||
      (FTN_STR(XFlag)[0]=='h'))  m = 1;
                           else  m = 0;

  if (channel[k]->nRead==0)  {

    // reading the atomic coordinates

    if (atom->Ter)  {
      *x    = 0.0;
      *y    = 0.0;
      *z    = 0.0;
      *occ  = 1.0;
      *BIso = 1.0;
      U[0]  = 1.0;
      U[1]  = 0.0;
      U[2]  = 0.0;
      U[3]  = 0.0;
      U[4]  = 0.0;
      U[5]  = 0.0;
    } else  {

      if (atom->WhatIsSet & mmdb::ASET_Coordinates)  {
        if ((FTN_STR(XFlag)[m]=='F') ||
            (FTN_STR(XFlag)[m]=='f'))  {
          //  receive fractional coordinates
          if (channel[k]->areCrystMatrices())  {
            channel[k]->Orth2Frac ( atom->x,atom->y,atom->z,xx,yy,zz );
            *x = (mmdb::machine::apireal)xx;
            *y = (mmdb::machine::apireal)yy;
            *z = (mmdb::machine::apireal)zz;
          } else  {
            *x = (mmdb::machine::apireal)atom->x;
            *y = (mmdb::machine::apireal)atom->y;
            *z = (mmdb::machine::apireal)atom->z;
            *iRet = RWBERR_NoMatrices;
          }
        } else  {
          // receive orthogonal coordinates
          *x = (mmdb::machine::apireal)atom->x;
          *y = (mmdb::machine::apireal)atom->y;
          *z = (mmdb::machine::apireal)atom->z;
        }
      } else  {
        *x = 0.0;
        *y = 0.0;
        *z = 0.0;
        *iRet = RWBERR_NoCoordinates;
      }

      // calculate isotropic Uf from Uo, and convert it
      // if necessary
      if (atom->WhatIsSet & mmdb::ASET_Anis_tFac)  {
        AU[0] = atom->u11;  // this intermediate array is
        AU[1] = atom->u22;  // required because of possible
        AU[2] = atom->u33;  // type difference between
        AU[3] = atom->u12;  // 'mmdb::machine::apireal' and 'realtype'
        AU[4] = atom->u13;
        AU[5] = atom->u23;
        *BIso = (mmdb::machine::apireal)(8.0*mmdb::Pi*mmdb::Pi*(AU[0]+AU[1]+AU[2])/3.0);
        if ((FTN_STR(BFlag)[0]=='F') ||
            (FTN_STR(BFlag)[0]=='f'))  {
          if (channel[k]->areCrystMatrices())
             channel[k]->Orth2Cryst ( AU );
          else if (*iRet==RWBERR_Ok)
             *iRet = RWBERR_NoMatrices;
        }
        for (i=0;i<6;i++)
          U[i] = (mmdb::machine::apireal)AU[i];
      } else  {
        for (i=0;i<6;i++)
          U[i] = 0.0;
        if (atom->WhatIsSet & mmdb::ASET_tempFactor)
          U[0] = (mmdb::machine::apireal)atom->tempFactor;
        else if (*iRet>=RWBERR_Ok)
          *iRet |= RWBWAR_NoTempFactor;
        *BIso = U[0];
      }

      // get occupancy now
      if (atom->WhatIsSet & mmdb::ASET_Occupancy)
        *occ = (mmdb::machine::apireal)atom->occupancy;
      else  {
        *occ = 0.0;
        if (*iRet>=RWBERR_Ok)  *iRet |= RWBWAR_NoOccupancy;
      }

    }

  } else  {

    // storing the atomic coordinates

    if (atom->Ter)  {
      atom->x = 0.0;
      atom->y = 0.0;
      atom->z = 0.0;
      atom->WhatIsSet |= mmdb::ASET_Coordinates;
      atom->occupancy  = 1.0;
      atom->tempFactor = 1.0;
      atom->u11 = 0.0;
      atom->u22 = 0.0;
      atom->u33 = 0.0;
      atom->u12 = 0.0;
      atom->u13 = 0.0;
      atom->u23 = 0.0;
    } else  {

      if ((FTN_STR(XFlag)[m]=='F') ||
          (FTN_STR(XFlag)[m]=='f'))  {
        //  convert fractional coordinates
        if (channel[k]->areCrystMatrices())  {
          xx = *x;
          yy = *y;
          zz = *z;
          channel[k]->Frac2Orth ( xx,yy,zz,atom->x,atom->y,atom->z );
          atom->WhatIsSet |= mmdb::ASET_Coordinates;
        } else  {
          atom->x = *x;
          atom->y = *y;
          atom->z = *z;
          *iRet   = RWBERR_NoMatrices;
          atom->WhatIsSet &= ~mmdb::ASET_Coordinates;
        }
      } else  {
        // store orthogonal coordinates
        atom->x = *x;
        atom->y = *y;
        atom->z = *z;
        atom->WhatIsSet |= mmdb::ASET_Coordinates;
      }

      atom->Het = (m>0);

      // calculate isotropic Uf from Uo, and convert it
      // if necessary
      if ((U[1]!=0.0) || (U[2]!=0.0))  {
        for (i=0;i<6;i++)
          AU[i] = U[i];
        if ((FTN_STR(BFlag)[0]=='F') ||
            (FTN_STR(BFlag)[0]=='f'))  {
          if (channel[k]->areCrystMatrices())
                channel[k]->Cryst2Orth ( AU );
          else  *iRet = RWBERR_NoMatrices;
        }
        *BIso = (mmdb::machine::apireal)(8.0*mmdb::Pi*mmdb::Pi*(AU[0]+AU[1]+AU[2])/3.0);
        atom->tempFactor = *BIso;
        atom->u11 = AU[0];
        atom->u22 = AU[1];
        atom->u33 = AU[2];
        atom->u12 = AU[3];
        atom->u13 = AU[4];
        atom->u23 = AU[5];
        atom->WhatIsSet |= mmdb::ASET_tempFactor | mmdb::ASET_Anis_tFac;
      } else  {
        *BIso = U[0];
        atom->tempFactor = *BIso;
        atom->u11 = 0.0;
        atom->u22 = 0.0;
        atom->u33 = 0.0;
        atom->u12 = 0.0;
        atom->u13 = 0.0;
        atom->u23 = 0.0;
        atom->WhatIsSet |= mmdb::ASET_tempFactor;
      }

      // store occupancy now
      atom->occupancy = *occ;
      atom->WhatIsSet |= mmdb::ASET_Occupancy;

    }

  }

  LastRC = *iRet;

}



FORTRAN_SUBR ( MMDB_F_SETCELL, mmdb_f_setcell,
               (   //   lengths-at-end list
                int     * iUnit,    // unit number
                mmdb::machine::apireal * a,        // cell parameter a, angstroms
                mmdb::machine::apireal * b,        // cell parameter b, angstroms
                mmdb::machine::apireal * c,        // cell parameter c, angstroms
                mmdb::machine::apireal * alpha,    // cell parameter alpha, degrees
                mmdb::machine::apireal * beta,     // cell parameter beta,  degrees
                mmdb::machine::apireal * gamma,    // cell parameter gamma, degrees
                int     * ArgNCode, // orthogonalization code, 1-6
                int     * iRet      // return code:
                                    //   RWBERR_Ok  - success
                                    //   RWBERR_NoChannel     if unit
                                    //              iUnit was not
                                    //              initialized
                                    //   RWBERR_NoFile        if unit
                                    //              has been disposed
                                    //   RWBERR_Disagreement  if a
                                    //              disagreement in
                        //              cell parameters
                        //              was found
                        //   RWBERR_NoOrthCode    if no
                                    //              orthogonalization
                        //              code was found
                        //   RWBERR_NoCheck       if check
                        //              of cell parameters
                        //              has failed.
                                    //   The last three returns would
                        // rather indicate a programming
                        // error in mmdb_rwbrook.cpp
               ), ( // lengths-in-structure list
                int     * iUnit,
                mmdb::machine::apireal * a,        mmdb::machine::apireal * b,    mmdb::machine::apireal * c,
                mmdb::machine::apireal * alpha,    mmdb::machine::apireal * beta, mmdb::machine::apireal * gamma,
                int     * ArgNCode, int     * iRet
               ), ( // lengths-follow list
                int     * iUnit,
                mmdb::machine::apireal * a,        mmdb::machine::apireal * b,    mmdb::machine::apireal * c,
                mmdb::machine::apireal * alpha,    mmdb::machine::apireal * beta, mmdb::machine::apireal * gamma,
                int     * ArgNCode, int     * iRet
               ) )  {
int  k;

  strcpy ( LastFunc,"MMDB_F_SetCell" );
  if (*iUnit>0)
    LastUnit = *iUnit;

  k = GetChannel ( LastUnit );
  if (k<0)
    *iRet = RWBERR_NoChannel;
  else
    *iRet = channel[k]->SetCell ( *a,*b,*c,*alpha,*beta,*gamma,
                                  *ArgNCode );

  LastRC = *iRet;

}


FORTRAN_SUBR ( MMDB_F_WBSPGRP, mmdb_f_wbspgrp,
               (   //   lengths-at-end list
                int   * iUnit,  // unit number; *iUnit<=0 means
                                // "the last mentioned unit"
                mmdb::machine::fpstr spGroup,  // space group
                int   * iRet,   // return code:
                                //   RWBERR_Ok  - success
                                //   RWBERR_NoChannel     if unit
                                //              iUnit was not
                                //              initialized
                                //   RWBERR_NoFile        if unit
                                //              has been disposed
                int spGroup_len // fortran-hidden length of spGroup
               ), ( // lengths-in-structure list
                int * iUnit, mmdb::machine::fpstr spGroup, int * iRet
               ), ( // lengths-follow list
                int * iUnit, mmdb::machine::fpstr spGroup, int spGroup_len,
                int * iRet
               )
       )  {
int            k;
mmdb::SymGroup spaceGroup;

  strcpy ( LastFunc,"MMDB_F_WBSpGrp" );
  if (*iUnit>0)
    LastUnit = *iUnit;

  k = GetChannel ( LastUnit );
  if (k<0)
    *iRet = RWBERR_NoChannel;
  else  {
//    GetStrTer ( spaceGroup,FTN_STR(spGroup),0,
//                sizeof(spaceGroup),FTN_LEN(spGroup) );
    mmdb::strcpy_ncss(spaceGroup,FTN_STR(spGroup),mmdb::IMin(FTN_LEN(spGroup),
                 sizeof(spaceGroup)-1) );
    *iRet = channel[k]->SetSpGroup ( spaceGroup );
  }

  LastRC = *iRet;

}



FORTRAN_SUBR ( MMDB_F_RBSPGRP, mmdb_f_rbspgrp,
               (   //   lengths-at-end list
                int   * iUnit,  // unit number; *iUnit<=0 means
                                // "the last mentioned unit"
                mmdb::machine::fpstr spGroup,  // space group
                int   * iRet,   // return code:
                                //   RWBERR_Ok  - success
                                //   RWBERR_NoChannel     if unit
                                //              iUnit was not
                                //              initialized
                                //   RWBERR_NoFile        if unit
                                //              has been disposed
                int spGroup_len // fortran-hidden length of spGroup
               ), ( // lengths-in-structure list
                int * iUnit, mmdb::machine::fpstr spGroup, int * iRet
               ), ( // lengths-follow list
                int * iUnit, mmdb::machine::fpstr spGroup, int spGroup_len,
                int * iRet
               )
       )  {
int  k;
char SpaceGroup[100];

  strcpy ( LastFunc,"MMDB_F_RBSpGrp" );
  if (*iUnit>0)
    LastUnit = *iUnit;

  SpaceGroup[0] = char(0);
  k = GetChannel ( LastUnit );
  if (k<0)  *iRet = RWBERR_NoChannel;
      else  *iRet = channel[k]->GetSpGroup ( SpaceGroup );

// all extra "superficial spaces" are killed in the following
  mmdb::CutSpaces ( SpaceGroup,mmdb::SCUTKEY_BEGEND );
  mmdb::strcpy_ns ( FTN_STR(spGroup),SpaceGroup,FTN_LEN(spGroup) );

  LastRC = *iRet;

}



FORTRAN_SUBR ( MMDB_F_WBCELL , mmdb_f_wbcell,
               (    // lengths-at-end list
                int     * iUnit,    // unit number; *iUnit<=0 means
                                    // "the last mentioned unit"
                mmdb::machine::apireal * ArgCell,  // array to accept the cell parameters
                                    // if ArgCell(1) is set to 0, then
                                    // the cell does not change
                int     * ArgNCode, // orthogonalisation code
                                    // if ArgNCode is set to 0, then
                                    // the orthogonalisation matrices
                                    // do not change
                int     * iRet      // return code
                                    //   RWBERR_Ok  - success
                                    //   RWBERR_NoChannel     if unit
                                    //              iUnit was not
                                    //              initialized
                                    //   RWBERR_NoFile        if unit
                                    //              has been disposed
               ), ( // lengths-in-structure list
                int * iUnit,    mmdb::machine::apireal * ArgCell,
                int * ArgNCode, int     * iRet
               ), ( // lengths-follow list
                int * iUnit,    mmdb::machine::apireal * ArgCell,
                int * ArgNCode, int     * iRet
               )
             )  {
int k;

  strcpy ( LastFunc,"MMDB_F_WBCell" );
  if (*iUnit>0)
    LastUnit = *iUnit;

  k = GetChannel ( LastUnit );
  if (k<0)
    *iRet = RWBERR_NoChannel;
  else
    *iRet = channel[k]->PutCell ( ArgCell[0],ArgCell[1],ArgCell[2],
                                  ArgCell[3],ArgCell[4],ArgCell[5],
                                  *ArgNCode );

  LastRC = *iRet;

}



FORTRAN_SUBR ( MMDB_F_RBCELL, mmdb_f_rbcell,
               (    // lengths-at-end list
                int     * iUnit,    // unit number
                mmdb::machine::apireal * celld,    // array to accept the cell parameters
                mmdb::machine::apireal * cvol,     // returns the cell volume
                int     * iRet      // return code
                        //   RWBERR_Ok  - success
                        //   RWBERR_NoChannel     if unit
                        //              iUnit was not
                        //              initialized
                        //   RWBERR_NoFile        if unit
                        //              has been disposed
                        //   RWBERR_Parameters    if the
                        //              cell parameters
                        //              were not set
                        //   RWBERR_NoOrthCode    if no
                                    //              orthogonalization
                        //              code was found
                        //   RWBERR_NoCheck       if check
                        //              of cell parameters
                        //              has failed.
                                    //   The last three returns would
                        // rather indicate a programming
                        // error in mmdb_rwbrook.cpp
               ), ( // lengths-in-structure list
                int     * iUnit,  mmdb::machine::apireal * celld,
                mmdb::machine::apireal * cvol,   int     * iRet
               ), ( // lengths-follow list
                int     * iUnit,  mmdb::machine::apireal * celld,
                mmdb::machine::apireal * cvol,   int     * iRet
               ) )  {
mmdb::realtype p[6];
mmdb::realtype v;
int            k,i,nc;

  strcpy ( LastFunc,"MMDB_F_RBCell" );
  if (*iUnit>0)
    LastUnit = *iUnit;

  k = GetChannel ( LastUnit );
  if (k<0)  {
    *iRet  = RWBERR_NoChannel;
    LastRC = *iRet;
    return;
  }

  *iRet = channel[k]->GetCell ( p[0],p[1],p[2],p[3],p[4],p[5],v,nc );

  if (*iRet==RWBERR_Ok)  {
    for (i=0;i<6;i++)
      celld[i] = (mmdb::machine::apireal)p[i];
    *cvol = (mmdb::machine::apireal)v;
  }

  LastRC = *iRet;

}



FORTRAN_SUBR ( MMDB_F_RBCELLN, mmdb_f_rbcelln,
               (    // lengths-at-end list
                int     * iUnit,    // unit number
                mmdb::machine::apireal * celld,    // array to accept the cell parameters
                mmdb::machine::apireal * cvol,     // returns the cell volume
                int     * ArgNCode, // returns the orthogonalization code, 1-6
                int     * iRet      // return code
                        //   RWBERR_Ok  - success
                        //   RWBERR_NoChannel     if unit
                        //              iUnit was not
                        //              initialized
                        //   RWBERR_NoFile        if unit
                        //              has been disposed
                        //   RWBERR_Parameters    if the
                        //              cell parameters
                        //              were not set
                        //   RWBERR_NoOrthCode    if no
                                    //              orthogonalization
                        //              code was found
                        //   RWBERR_NoCheck       if check
                        //              of cell parameters
                        //              has failed.
                                    //   The last three returns would
                        // rather indicate a programming
                        // error in mmdb_rwbrook.cpp
               ), ( // lengths-in-structure list
                int * iUnit,    mmdb::machine::apireal * celld, mmdb::machine::apireal * cvol,
                int * ArgNCode, int     * iRet
               ), ( // lengths-follow list
                int * iUnit,    mmdb::machine::apireal * celld, mmdb::machine::apireal * cvol,
                int * ArgNCode, int     * iRet
               ) )  {
mmdb::realtype p[6];
mmdb::realtype v;
int            k,i,nc;

  strcpy ( LastFunc,"MMDB_F_RBCellN" );
  if (*iUnit>0)
    LastUnit = *iUnit;

  k = GetChannel ( LastUnit );
  if (k<0)  {
    *iRet  = RWBERR_NoChannel;
    LastRC = *iRet;
    return;
  }

  *iRet = channel[k]->GetCell ( p[0],p[1],p[2],p[3],p[4],p[5],v,nc );
  if (*iRet==RWBERR_Ok)  {
    for (i=0;i<6;i++)
      celld[i] = (mmdb::machine::apireal)p[i];
    *cvol     = (mmdb::machine::apireal)v;
    *ArgNCode = nc;
  }

  LastRC = *iRet;

}



FORTRAN_SUBR ( MMDB_F_RBRCEL, mmdb_f_rbrcel,
               (    // lengths-at-end list
                int     * iUnit,    // unit number
                mmdb::machine::apireal * rcell,    // array to accept the reciprocal
                                    // cell parameters
                mmdb::machine::apireal * rvol,     // returns the reciprocal cell volume
                int     * iRet      // return code
                        //   RWBERR_Ok  - success
                        //   RWBERR_NoChannel     if unit
                        //              iUnit was not
                        //              initialized
                        //   RWBERR_NoFile        if unit
                        //              has been disposed
                        //   RWBERR_Parameters    if the
                        //              cell parameters
                        //              were not set
                        //   RWBERR_NoOrthCode    if no
                                    //              orthogonalization
                        //              code was found
                        //   RWBERR_NoCheck       if check
                        //              of cell parameters
                        //              has failed.
                                    //   The last three returns would
                        // rather indicate a programming
                        // error in mmdb_rwbrook.cpp
               ), ( // lengths-in-structure list
                int * iUnit,    mmdb::machine::apireal * rcell, mmdb::machine::apireal * rvol,
                int * iRet
               ), ( // lengths-follow list
                int * iUnit,    mmdb::machine::apireal * rcell, mmdb::machine::apireal * rvol,
                int * iRet
               ) )  {
mmdb::realtype p[6];
mmdb::realtype v;
int            k,i;

  strcpy ( LastFunc,"MMDB_F_RBRCel" );
  if (*iUnit>0)
    LastUnit = *iUnit;

  k = GetChannel ( LastUnit );
  if (k<0)  {
    *iRet  = RWBERR_NoChannel;
    LastRC = *iRet;
    return;
  }

  *iRet = channel[k]->GetRCell ( p[0],p[1],p[2],p[3],p[4],p[5],v );
  if (*iRet==RWBERR_Ok)  {
    for (i=0;i<6;i++)
      rcell[i] = (mmdb::machine::apireal)p[i];
    *rvol = (mmdb::machine::apireal)v;
  }

  LastRC = *iRet;

}



FORTRAN_SUBR ( MMDB_F_RBORF, mmdb_f_rborf,
               (     // lengths-at-end list
                int     * iUnit, // unit number
                mmdb::machine::apireal * RO,    // array for orthogonalising matrix
                mmdb::machine::apireal * RF,    // array for fractionalising matrix
                int     * LCode, // buffer for orthogonalisation code
                int     * iRet   // return code:
                                 //   RWBERR_Ok  - success
                                 //   RWBERR_NoChannel     if unit
                                 //              iUnit was not
                                 //              initialized
                                 //   RWBERR_NoFile        if unit
                                 //              has been disposed
                                 //   RWBERR_NoMatrices    if the
                                 //              orthogonalisation
                                 //              matrices were not
                                 //              calculated
               ), (  // lengths-in-structure list
                int * iUnit, mmdb::machine::apireal * RO, mmdb::machine::apireal * RF,
                int * LCode, int * iRet
               ), (  // lengths-follow list
                int * iUnit, mmdb::machine::apireal * RO, mmdb::machine::apireal * RF,
                int * LCode, int * iRet )
               )  {
int         i,j,k,l;
mmdb::PCryst Cryst;

  strcpy ( LastFunc,"MMDB_F_RBORF" );
  if (*iUnit>0)
    LastUnit = *iUnit;

  k = GetChannel ( LastUnit );
  if (k<0)  {
    *iRet  = RWBERR_NoChannel;
    LastRC = *iRet;
    return;
  }

  Cryst = channel[k]->GetCryst();
  if (Cryst==NULL)  {
    *iRet  = RWBERR_NoFile;
    LastRC = *iRet;
    return;
  }

  *iRet = RWBERR_Ok;

  l = 0;

  if (RO[0]<=0.0000000001)  {
    for (j=0;j<4;j++)
      for (i=0;i<4;i++)  {
        RF[l] = (mmdb::machine::apireal)Cryst->RF[i][j];
        RO[l] = (mmdb::machine::apireal)Cryst->RO[i][j];
        l++;
      }
    *LCode = Cryst->NCode;
    if (!(Cryst->WhatIsSet & mmdb::CSET_Transforms))
      *iRet = RWBERR_NoMatrices;
  } else  {
    for (j=0;j<4;j++)
      for (i=0;i<4;i++)  {
        Cryst->RF[i][j] = RF[l];
        Cryst->RO[i][j] = RO[l];
        l++;
      }
    Cryst->NCode = *LCode;
    Cryst->WhatIsSet |= mmdb::CSET_Transforms;
  }

  LastRC = *iRet;

}


FORTRAN_SUBR ( MMDB_F_ORTHMAT, mmdb_f_orthmat,
               (     // lengths-at-end list
                int     * iUnit, // unit number; *iUnit<=0 means
                                 // "the last mentioned unit"
                mmdb::machine::apireal * Cell,  // array of cell parameters:
                                 //  Cell(1) - a   Cell(4) - alpha
                                 //  Cell(2) - b   Cell(5) - beta
                                 //  Cell(3) - c   Cell(6) - gamma
                mmdb::machine::apireal * Vol,   // returns cell volume
                mmdb::machine::apireal * RRR,   // array (3,3,6), returns
                                 // orthogonalisation matrices
                int     * iRet   // return code:
                                 //   RWBERR_Ok  - success
                                 //   RWBERR_NoChannel     if unit
                                 //              iUnit was not
                                 //              initialized
                                 //   RWBERR_NoFile        if unit
                                 //              has been disposed
                                 //   RWBERR_NoMatrices    if the
                                 //              orthogonalisation
                                 //              matrices were not
                                 //              calculated
               ), ( // lengths-in-structure list
                int     * iUnit, mmdb::machine::apireal * Cell, mmdb::machine::apireal * Vol,
                mmdb::machine::apireal * RRR,   int * iRet
               ), ( // lengths-follow list
                int     * iUnit, mmdb::machine::apireal * Cell, mmdb::machine::apireal * Vol,
                mmdb::machine::apireal * RRR,   int * iRet
               )
             )  {
int            i,j,k,l,m;
mmdb::PCryst   Cryst;
mmdb::realtype CelDel;

  strcpy ( LastFunc,"MMDB_F_OrthMat" );
  if (*iUnit>0)
    LastUnit = *iUnit;

  k = GetChannel ( LastUnit );
  if (k<0)  {
    *iRet  = RWBERR_NoChannel;
    LastRC = *iRet;
    return;
  }

  Cryst = channel[k]->GetCryst();
  if (Cryst==NULL)  {
    *iRet  = RWBERR_NoFile;
    LastRC = *iRet;
    return;
  }

  CelDel = 0.0;
  if (Cell[0]>0.0)  {
    if ((Cryst->WhatIsSet & mmdb::CSET_CellParams)==mmdb::CSET_CellParams)  {
      CelDel = fabs((Cell[0]-Cryst->a)/Cell[0]);
      if (Cell[1]!=0.0)
        CelDel = mmdb::RMax(CelDel,fabs((Cell[1]-Cryst->b)/Cell[1]));
      if (Cell[2]!=0.0)
        CelDel = mmdb::RMax(CelDel,fabs((Cell[2]-Cryst->c)/Cell[2]));
      if (Cell[3]!=0.0)
        CelDel = mmdb::RMax(CelDel,fabs((Cell[3]-Cryst->alpha)/Cell[3]));
      if (Cell[4]!=0.0)
        CelDel = mmdb::RMax(CelDel,fabs((Cell[4]-Cryst->beta )/Cell[4]));
      if (Cell[5]!=0.0)
        CelDel = mmdb::RMax(CelDel,fabs((Cell[5]-Cryst->gamma)/Cell[5]));
      if (FSimRWBROOK && (CelDel>0.01))
        printf ( "\n Inconsistency in Cell Dimensions"
                 " - replacing old:\n"
                 " Old cell:   "
                 "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n"
                 " New cell:   "
                 "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n",
                 Cryst->a,Cryst->b,Cryst->c,
                 Cryst->alpha,Cryst->beta,Cryst->gamma,
                 Cell[0],Cell[1],Cell[2],Cell[3],Cell[4],Cell[5] );
    }
    Cryst->a     = Cell[0];
    Cryst->b     = Cell[1];
    Cryst->c     = Cell[2];
    Cryst->alpha = Cell[3];
    Cryst->beta  = Cell[4];
    Cryst->gamma = Cell[5];
    Cryst->WhatIsSet |= mmdb::CSET_CellParams;
  } else  {
    Cell[0] = (mmdb::machine::apireal)Cryst->a;
    Cell[1] = (mmdb::machine::apireal)Cryst->b;
    Cell[2] = (mmdb::machine::apireal)Cryst->c;
    Cell[3] = (mmdb::machine::apireal)Cryst->alpha;
    Cell[4] = (mmdb::machine::apireal)Cryst->beta;
    Cell[5] = (mmdb::machine::apireal)Cryst->gamma;
  }

  if ((Cryst->WhatIsSet & mmdb::CSET_CellParams)!=mmdb::CSET_CellParams)  {
    *iRet  = RWBERR_NoCellParams;
    LastRC = *iRet;
    return;
  }

  *iRet  = RWBERR_Ok;

  //  Cryst->CalcOrthMatrices();  <-- old version, changed 09.01.2004
  Cryst->CalcCoordTransforms();
  Cryst->WhatIsSet |= mmdb::CSET_Transforms;

  if (CelDel>0.01)  *Vol = -(mmdb::machine::apireal)Cryst->Vol;
              else  *Vol =  (mmdb::machine::apireal)Cryst->Vol;

  l = 0;
  for (j=0;j<3;j++)
    for (i=0;i<3;i++)
      for (m=0;m<6;m++)
        RRR[l++] = (mmdb::machine::apireal)Cryst->RR[m][j][i];

  LastRC = *iRet;

}


FORTRAN_SUBR ( MMDB_F_CVANISOU, mmdb_f_cvanisou,
               (     // lengths-at-end list
                int     * iUnit, // unit number; *iUnit<=0 means
                                 // "the last mentioned unit"
                mmdb::machine::apireal * U,     // array of coordinates to convert
                int     * iFlag, // =0: convert from fract. to orthog.
                                 // =1: convert from orthog. to fract.
                int     * iRet   // return code:
                                 //   RWBERR_Ok  - success
                                 //   RWBERR_NoChannel     if unit
                                 //              iUnit was not
                                 //              initialized
                                 //   RWBERR_NoFile        if unit
                                 //              has been disposed
                                 //   RWBERR_NoMatrices    if the
                                 //              orthogonalisation
                                 //              matrices were not
                                 //              calculated
               ), ( // lengths-in-structure list
                int * iUnit, mmdb::machine::apireal * U, int * iFlag, int * iRet
               ), ( // lengths-follow list
                int * iUnit, mmdb::machine::apireal * U, int * iFlag, int * iRet
               )
             )  {
int            k,i;
mmdb::PCryst   Cryst;
mmdb::realtype U1[6];

  strcpy ( LastFunc,"MMDB_F_CVAnisou" );
  if (*iUnit>0)
    LastUnit = *iUnit;

  k = GetChannel ( LastUnit );
  if (k<0)  {
    *iRet  = RWBERR_NoChannel;
    LastRC = *iRet;
    return;
  }

  Cryst = channel[k]->GetCryst();
  if (Cryst==NULL)  {
    *iRet  = RWBERR_NoFile;
    LastRC = *iRet;
    return;
  }

  *iRet = RWBERR_Ok;
  for (i=0;i<6;i++)
    U1[i] = U[i];

  if (iFlag==0)  {
    if (!Cryst->Cryst2Orth(U1))  *iRet = RWBERR_NoMatrices;
  } else  {
    if (!Cryst->Orth2Cryst(U1))  *iRet = RWBERR_NoMatrices;
  }

  if (*iRet==RWBERR_Ok)
    for (i=0;i<6;i++)
      U[i] = (mmdb::machine::apireal)U1[i];

  LastRC = *iRet;

}



FORTRAN_SUBR ( MMDB_F_WREMARK, mmdb_f_wremark,
               (    // lengths-at-end list
                int     * iUnit, // unit number; *iUnit<=0 means
                                 // "the last mentioned unit"
                mmdb::machine::fpstr     Line,  // line to be added
                int     * iRet,  // return code:
                                 //   RWBERR_Ok  - success
                                 //   RWBERR_NoChannel     if unit
                                 //              iUnit was not
                                 //              initialized
                                 //   RWBERR_NoFile        if unit
                                 //              has been disposed
                                 // other return codea are those
                                 // returned by xyzopen1_(..)
                int    Line_len  // fortran-hidden length of Line
               ), ( // lengths-in-structure list
                int * iUnit, mmdb::machine::fpstr Line, int * iRet
               ), ( // lengths-follow list
                int * iUnit, mmdb::machine::fpstr Line, int Line_len, int *iRet
               )
             )  {
int  k;
char S[500];

  strcpy ( LastFunc,"MMDB_F_WRemark" );
  if (*iUnit>0)
    LastUnit = *iUnit;

  k = GetChannel ( LastUnit );
  if (k<0)  {
    *iRet  = RWBERR_NoChannel;
    LastRC = *iRet;
    return;
  }

  if (channel[k]->MMDBManager)  {
    mmdb::GetStrTer ( S,FTN_STR(Line),FTN_LEN(Line),sizeof(S),FTN_LEN(Line) );
    *iRet =  channel[k]->MMDBManager->PutPDBString ( S );
  } else
    *iRet = RWBERR_NoFile;

  LastRC = *iRet;

}


/*
FORTRAN_SUBR ( RBRINV, rbrinv,
               ( mmdb::machine::apireal * A, mmdb::machine::apireal * AI ),
               ( mmdb::machine::apireal * A, mmdb::machine::apireal * AI ),
               ( mmdb::machine::apireal * A, mmdb::machine::apireal * AI ) )  {
mat44  A1,AI1;
int    i,j,k;

  k = 0;
  for (j=0;j<4;j++)
    for (i=0;i<4;i++)
      A1[j][i] = A[k++];

  Mat4Inverse ( A1,AI1 );

  k = 0;
  for (j=0;j<4;j++)
    for (i=0;i<4;i++)
      AI[k++] = AI1[j][i];

}
*/
/*
FORTRAN_SUBR ( RES3TO1, res3to1,
               (     // lengths-at-end list
                mmdb::machine::fpstr ResNm3,   // 3-char name, 4th char
                                // will be set blank
                mmdb::machine::fpstr ResNm1,   // 1-char name
                int ResNm3_len, // fortran-hidden length of ResNm3
                int ResNm1_len  // fortran-hidden length of ResNm3
               ), ( // lengths-in-structure list
                mmdb::machine::fpstr ResNm3, mmdb::machine::fpstr ResNm1
               ), ( // lengths-follow list
                mmdb::machine::fpstr ResNm3, int ResNm3_len,
                mmdb::machine::fpstr ResNm1, int ResNm1_len
               )
             )  {
int i;

  if (FTN_STR(ResNm3)[0]==' ')  {
    for (i=0;i<nResNames;i++)
      if ((FTN_STR(ResNm3)[0]==ResidueName[i][0]) &&
          (FTN_STR(ResNm3)[1]==ResidueName[i][1]) &&
          (FTN_STR(ResNm3)[2]==ResidueName[i][2]))  {
        FTN_STR(ResNm1)[0] = ResidueName1[i];
        return;
      }
    FTN_STR(ResNm1)[0] = ResidueName1[nResNames-1];
    return;
  }

  if (FTN_STR(ResNm1)[0]==' ')  {
    for (i=0;i<nResNames;i++)
      if (FTN_STR(ResNm1)[0]==ResidueName1[i])  {
        FTN_STR(ResNm3)[0] = ResidueName[i][0];
        FTN_STR(ResNm3)[1] = ResidueName[i][1];
        FTN_STR(ResNm3)[2] = ResidueName[i][2];
        FTN_STR(ResNm3)[3] = ' ';
        return;
      }
    FTN_STR(ResNm3)[0] = ResidueName[nResNames-1][0];
    FTN_STR(ResNm3)[1] = ResidueName[nResNames-1][1];
    FTN_STR(ResNm3)[2] = ResidueName[nResNames-1][2];
    FTN_STR(ResNm3)[3] = ' ';
    return;
  }

}
*/

static mmdb::pstr MSG_NoChannel       = mmdb::pstr("unassigned unit");
static mmdb::pstr MSG_NoFile          = mmdb::pstr("unassigned unit or disposed file");
static mmdb::pstr MSG_NoLogicalName   = mmdb::pstr("logical name does not exist");

static mmdb::pstr MSG_CantOpenFile    = mmdb::pstr("cannot open a file");
static mmdb::pstr MSG_WrongInteger    = mmdb::pstr("unrecognized integer at reading a file");
static mmdb::pstr MSG_WrongModelNo    = mmdb::pstr("wrong model number read from a file");
static mmdb::pstr MSG_DuplicatedModel = mmdb::pstr("duplicated model number");
static mmdb::pstr MSG_ForeignFile     = mmdb::pstr("unknown file format");
static mmdb::pstr MSG_WrongEdition    = mmdb::pstr("unknown file version");

static mmdb::pstr MSG_ATOM_Unrecognd  = mmdb::pstr("unrecognized data in coordinate section");
static mmdb::pstr MSG_ATOM_AlreadySet = mmdb::pstr("duplicate atom serial number");
static mmdb::pstr MSG_ATOM_NoResidue  = mmdb::pstr("residue for atom cannot be found");
static mmdb::pstr MSG_ATOM_Unmatch    = mmdb::pstr("ambiguous data in coordinate section");

static mmdb::pstr MSG_NoAdvance       = mmdb::pstr("atom position was not advanced");
static mmdb::pstr MSG_EmptyPointer    = mmdb::pstr("atom was not allocated");
static mmdb::pstr MSG_NoMatrices      = mmdb::pstr("no coordinate transformation matrices");

static mmdb::pstr MSG_NoCoordinates   = mmdb::pstr("no atom coordinates set");

static mmdb::pstr MSG_Disagreement    = mmdb::pstr("ambiguous cell parameters");
static mmdb::pstr MSG_NoOrthCode      = mmdb::pstr("no orthogonalization code");
static mmdb::pstr MSG_NoCheck         = mmdb::pstr("missing check of cell parameters");

static mmdb::pstr MSG_NoCellParams    = mmdb::pstr("no cell parameters");

static mmdb::pstr MSG_NotACIFFile     = mmdb::pstr("not a CIF file: 'data_' tag missing");
static mmdb::pstr MSG_NoData          = mmdb::pstr("expected data is not met at reading a file");
static mmdb::pstr MSG_UnrecognCIFItems  = mmdb::pstr("unrecognized CIF items (syntax error?)");
static mmdb::pstr MSG_MissingCIFField   = mmdb::pstr("missing CIF data field");
static mmdb::pstr MSG_EmptyCIFLoop      = mmdb::pstr("CIF loop does not contain any data");
static mmdb::pstr MSG_UnexpEndOfCIF     = mmdb::pstr("unexpected end of CIF file");
static mmdb::pstr MSG_MissgCIFLoopField = mmdb::pstr("CIF loop is incomplete");
static mmdb::pstr MSG_NotACIFStructure  = mmdb::pstr("wrong use of CIF structure (as a loop?)");
static mmdb::pstr MSG_NotACIFLoop       = mmdb::pstr("wrong use of CIF loop (as a structure?)");
static mmdb::pstr MSG_WrongReal         = mmdb::pstr("unrecognized real at reading a file");

static mmdb::pstr MSG_WrongChainID      = mmdb::pstr("Wrong or inconsistent chain ID");
static mmdb::pstr MSG_WrongEntryID      = mmdb::pstr("Wrong or insonsistent entry ID");
static mmdb::pstr MSG_SEQRES_serNum     = mmdb::pstr("Wrong serial number in SEQRES");
static mmdb::pstr MSG_SEQRES_numRes     = mmdb::pstr("Wrong number of residues in SEQRES");
static mmdb::pstr MSG_SEQRES_extraRes   = mmdb::pstr("Extra residues in SEQRES");
static mmdb::pstr MSG_NCSM_Unrecogn     = mmdb::pstr("Unrecognized item in NCSM cards");
static mmdb::pstr MSG_NCSM_AlreadySet   = mmdb::pstr("Attempt to reset NCSM");
static mmdb::pstr MSG_NCSM_WrongSerial  = mmdb::pstr("Wrong serial number in NCSM cards");
static mmdb::pstr MSG_NCSM_UnmatchIG    = mmdb::pstr("Unmatched IG parameter in NCSM cards");
static mmdb::pstr MSG_NoModel           = mmdb::pstr("MMDB's error in structuring models");
static mmdb::pstr MSG_NoSheetID         = mmdb::pstr("No sheet ID on SHEET card(s)");
static mmdb::pstr MSG_WrongSheetID      = mmdb::pstr("Wrong sheet ID on SHEET card(s)");
static mmdb::pstr MSG_WrongStrandNo     = mmdb::pstr("Wrong strand no. on SHEET card(s)");
static mmdb::pstr MSG_WrongNofStrands   = mmdb::pstr("Wrong number of strands in sheet");
static mmdb::pstr MSG_WrongSheetOrder   = mmdb::pstr("Wrong sheet ordering");
static mmdb::pstr MSG_HBondInconsistency = mmdb::pstr("Inconsistency in H-bonds");
static mmdb::pstr MSG_EmptyResidueName  = mmdb::pstr("No (blank) residue name");
static mmdb::pstr MSG_DuplicateSeqNum   = mmdb::pstr("Duplicated sequence number and insertion code");
static mmdb::pstr MSG_GeneralError1     = mmdb::pstr("MMDB's general error #1");


static mmdb::pstr MSG_Error1          = mmdb::pstr("internal error #1 -- report to developer");
static mmdb::pstr MSG_Error2          = mmdb::pstr("internal error #2 -- report to developer");
static mmdb::pstr MSG_Error3          = mmdb::pstr("internal error #3 -- report to developer");

static mmdb::pstr MSG_Unknown         = mmdb::pstr("unknown return code");


#define nWarnings  7

static int RWBWarCode[nWarnings] = {
  RWBWAR_RewOutput,
  RWBWAR_FileTop,
  RWBWAR_WrongSerial,
  RWBWAR_UnkFormFactor,
  RWBWAR_AmbFormFactor,
  RWBWAR_NoOccupancy,
  RWBWAR_NoTempFactor
};

static mmdb::pstr RWBWarning[nWarnings] = {
  mmdb::pstr("output file rewind"),
  mmdb::pstr("rewind or backspace at top of file"),
  mmdb::pstr("atom serial number does not match position"),
  mmdb::pstr("unknown form factor encountered"),
  mmdb::pstr("ambiguous form factor encountered"),
  mmdb::pstr("occupancy was not set"),
  mmdb::pstr("temperature factor was not set")
};



FORTRAN_SUBR ( RBERRSTOP, rberrstop,
               (   //    lengths-at-end list
                int * iPlace, // (unique) identificator inside an application
                int * iRet,   // return code to check
                int * iUnit,  // unit number
                int * iStop   // if 0 then stop if error
               ), ( // lengths-in-structure list
                int * iPlace, int * iRet,
                int * iUnit,  int * iStop
               ), ( // lengths-follow list
                int * iPlace, int * iRet,
                int * iUnit,  int * iStop
               ) )  {
int        i,k,lcount;
mmdb::pstr Msg;
char       ErrLine[500];

  strcpy ( ErrLine,"" );
  lcount = -11;
  k      = GetChannel(*iUnit);

  switch (*iRet)  {

    case RWBERR_Ok                : return;

    case RWBERR_NoChannel         : Msg = MSG_NoChannel;          break;
    case RWBERR_NoFile            : Msg = MSG_NoFile;             break;
    case RWBERR_NoLogicalName     : Msg = MSG_NoLogicalName;      break;
    case RWBERR_CantOpenFile      : Msg = MSG_CantOpenFile;       break;


    case RWBERR_WrongInteger      : Msg = MSG_WrongInteger;       break;
    case RWBERR_WrongModelNo      : Msg = MSG_WrongModelNo;       break;
    case RWBERR_DuplicatedModel   : Msg = MSG_DuplicatedModel;    break;

    case RWBERR_ForeignFile       : Msg = MSG_ForeignFile;        break;
    case RWBERR_WrongEdition      : Msg = MSG_WrongEdition;       break;

    case RWBERR_ATOM_Unrecognd    : Msg = MSG_ATOM_Unrecognd;     break;
    case RWBERR_ATOM_AlreadySet   : Msg = MSG_ATOM_AlreadySet;    break;
    case RWBERR_ATOM_NoResidue    : Msg = MSG_ATOM_NoResidue;     break;
    case RWBERR_ATOM_Unmatch      : Msg = MSG_ATOM_Unmatch;       break;

    case RWBERR_NoAdvance         : Msg = MSG_NoAdvance;          break;
    case RWBERR_EmptyPointer      : Msg = MSG_EmptyPointer;       break;
    case RWBERR_NoMatrices        : Msg = MSG_NoMatrices;         break;

    case RWBERR_NoCoordinates     : Msg = MSG_NoCoordinates;      break;

    case RWBERR_Disagreement      : Msg = MSG_Disagreement;       break;
    case RWBERR_NoOrthCode        : Msg = MSG_NoOrthCode;         break;
    case RWBERR_NoCheck           : Msg = MSG_NoCheck;            break;

    case RWBERR_NoCellParams      : Msg = MSG_NoCellParams;       break;

    case RWBERR_NotACIFFile       : Msg = MSG_NotACIFFile;        break;
    case RWBERR_NoData            : Msg = MSG_NoData;             break;
    case RWBERR_UnrecognCIFItems  : Msg = MSG_UnrecognCIFItems;   break;
    case RWBERR_MissingCIFField   : Msg = MSG_MissingCIFField;    break;
    case RWBERR_EmptyCIFLoop      : Msg = MSG_EmptyCIFLoop;       break;
    case RWBERR_UnexpEndOfCIF     : Msg = MSG_UnexpEndOfCIF;      break;
    case RWBERR_MissgCIFLoopField : Msg = MSG_MissgCIFLoopField;  break;
    case RWBERR_NotACIFStructure  : Msg = MSG_NotACIFStructure;   break;
    case RWBERR_NotACIFLoop       : Msg = MSG_NotACIFLoop;        break;
    case RWBERR_WrongReal         : Msg = MSG_WrongReal;          break;

    case RWBERR_WrongChainID      : Msg = MSG_WrongChainID;       break;
    case RWBERR_WrongEntryID      : Msg = MSG_WrongEntryID;       break;
    case RWBERR_SEQRES_serNum     : Msg = MSG_SEQRES_serNum;      break;
    case RWBERR_SEQRES_numRes     : Msg = MSG_SEQRES_numRes;      break;
    case RWBERR_SEQRES_exraRes    : Msg = MSG_SEQRES_extraRes;    break;
    case RWBERR_NCSM_Unrecogn     : Msg = MSG_NCSM_Unrecogn;      break;
    case RWBERR_NCSM_AlreadySet   : Msg = MSG_NCSM_AlreadySet;    break;
    case RWBERR_NCSM_WrongSerial  : Msg = MSG_NCSM_WrongSerial;   break;
    case RWBERR_NCSM_UnmatchIG    : Msg = MSG_NCSM_UnmatchIG;     break;
    case RWBERR_NoModel           : Msg = MSG_NoModel;            break;
    case RWBERR_NoSheetID         : Msg = MSG_NoSheetID;          break;
    case RWBERR_WrongSheetID      : Msg = MSG_WrongSheetID;       break;
    case RWBERR_WrongStrandNo     : Msg = MSG_WrongStrandNo;      break;
    case RWBERR_WrongNofStrands   : Msg = MSG_WrongNofStrands;    break;
    case RWBERR_WrongSheetOrder   : Msg = MSG_WrongSheetOrder;    break;
    case RWBERR_HBondInconsis     : Msg = MSG_HBondInconsistency; break;
    case RWBERR_EmptyResidueName  : Msg = MSG_EmptyResidueName;   break;
    case RWBERR_DuplicateSeqNum   : Msg = MSG_DuplicateSeqNum;    break;
    case RWBERR_GeneralError1     : Msg = MSG_GeneralError1;      break;

    case RWBERR_Error1            : Msg = MSG_Error1;             break;
    case RWBERR_Error2            : Msg = MSG_Error2;             break;
    case RWBERR_Error3            : Msg = MSG_Error3;             break;

    default                     :
      if ((*iRet & RWBWAR_Warning)==RWBWAR_Warning)  {
        Msg = NULL;
        printf ( "\n\n *** Warning(s): point code unit    function\n" );
        printf ( " ***             %5i %4i %4i    %s\n",
                           *iPlace,*iRet,*iUnit,LastFunc );
        if (k>=0)
          printf ( " *** file   : %s\n",channel[k]->FName );
        for (i=0;i<nWarnings;i++)
          if ((*iRet & RWBWarCode[i])==RWBWarCode[i])  {
            Msg = RWBWarning[i];
            printf ( " *** warning: %s\n",Msg );
            if ((*iRet & RWBWAR_WrongSerial)==RWBWAR_WrongSerial)  {
              if (k>0)
                printf ( " *** position %i, serial number %i\n",
                         channel[k]->fPos,LastSer );
              else
                printf ( " *** position unavailable, serial number %i\n",
                         LastSer );
            }
          }
        if (!Msg)
          printf ( " *** warning: unknown warning code" );
        return;
      } else
        Msg = MSG_Unknown;
  }

  if ((k>=0) && (
      ((*iRet<=RWBERR_WrongInteger)   && (*iRet>=RWBERR_DuplicatedModel)) ||
      ((*iRet<=RWBERR_ATOM_Unrecognd) && (*iRet>=RWBERR_ATOM_Unmatch))    ||
      ((*iRet<=RWBERR_NoData)         && (*iRet>=RWBERR_DuplicateSeqNum))
     ))
    channel[k]->GetInputBuffer ( ErrLine,lcount );

  printf ( " \n *** RWBROOK error: point code unit    function\n"    );
  printf ( " ***                %5i %4i %4i    %s\n",*iPlace,*iRet,
                                                     *iUnit,LastFunc );
  k = GetChannel(*iUnit);
  if (k>=0)
    printf ( " *** file   : %s\n",channel[k]->FName );

  printf ( " *** reason : %s\n",Msg );
  if (lcount>=0)
    printf ( " ***          at input line #%i:\n"
             " %s\n",lcount,ErrLine );
  else if (lcount==-1)
    printf ( " ***          at taking the following data from CIF:\n"
             "              %s\n",ErrLine );

  if (*iStop==0)  {  // will stop it
    printf ( " *** Execution stopped.\n \n" );
    FORTRAN_CALL ( MMDB_F_QUIT, mmdb_f_quit,(),(),() );
    // xyzquit_();
    exit(0);
  } else     // just warn, but no guarantee that it will not crash
    printf ( " *** continue running, may crash ...\n \n" );

}



FORTRAN_SUBR ( RBCHECKERR, rbcheckerr,
               (   //    lengths-at-end list
                int * iPlace, // (unique) identificator inside an application
                int * iStop   // if 0 then stop if error
               ), ( // lengths-in-structure list
                int * iPlace, int * iStop
               ), ( // lengths-follow list
                int * iPlace, int * iStop
               ) )  {
  FORTRAN_CALL ( RBERRSTOP, rberrstop,
                 ( iPlace,&LastRC,&LastUnit,iStop ),
                 ( iPlace,&LastRC,&LastUnit,iStop ),
                 ( iPlace,&LastRC,&LastUnit,iStop ) );
}


/* hybrid-36 encoder: converts integer value to string result

      iwidth: must be 4 (e.g. for residue sequence numbers)
                  or 5 (e.g. for atom serial numbers)

      value: integer value to be converted

      strval: char array containing string result
*/

FORTRAN_SUBR ( HY36ENCODE_F, hy36encode_f,
               (const int *iwidth, int *value,
                mmdb::machine::fpstr strval, int strval_len),
               (const int *iwidth, int *value,
                mmdb::machine::fpstr strval),
               (const int *iwidth, int *value,
                mmdb::machine::fpstr strval, int strval_len))
{
  unsigned width;
  char result[6];

  width = (unsigned) *iwidth;

  if (hy36encode(width, *value, result)) {
    printf("problem in hy36encode_f! \n");
  }
  mmdb::strcpy_ns(FTN_STR(strval),result,FTN_LEN(strval));

}

/*  hybrid-36 decoder: converts string s to integer result

      iwidth: must be 4 (e.g. for residue sequence numbers)
                  or 5 (e.g. for atom serial numbers)

      strval: string to be converted

      value: integer holding the conversion result
*/


FORTRAN_SUBR ( HY36DECODE_F, hy36decode_f,
               (const int *iwidth,
                mmdb::machine::fpstr strval, int *value, int strval_len),
               (const int *iwidth,
                mmdb::machine::fpstr strval, int *value),
               (const int *iwidth,
                mmdb::machine::fpstr strval, int strval_len, int *value))
{ UNUSED_ARGUMENT(strval);
  unsigned width;
  size_t length = FTN_LEN(strval);
  char* s;

  width = (unsigned) *iwidth;
  s = (char *) malloc((length+1)*sizeof(char));
  s[length] = '\0';

  if (hy36decode(width, s, width, value)) {
    printf("problem in hy36decode_f! \n");
  }

}

