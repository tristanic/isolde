//  $Id: mmdb_manager.cpp $
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
//    15.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  mmdb_manager <implementation>
//       ~~~~~~~~~
//       Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::Manager  ( MMDB file manager )
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//


#include <string.h>

#include "mmdb_manager.h"

namespace mmdb  {

//  =====================   Manager   =======================

  Manager::Manager() : BondManager()  {
  }

  Manager::Manager ( io::RPStream Object ) : BondManager(Object)  {
  }

  Manager::~Manager()  {}

  void  Manager::Copy ( PManager MMDB, COPY_MASK CopyMask )  {
  PModel  mdl;
  PPChain chain;
  PChain  ch;
  ChainID chID;
  int     i,j, nchains;

    if (CopyMask & MMDBFCM_Flags)  Flags = MMDB->Flags;

    if (CopyMask & MMDBFCM_Title)  title.Copy ( &(MMDB->title) );
    if (CopyMask & MMDBFCM_Cryst)  cryst.Copy ( &(MMDB->cryst) );

    if (CopyMask & MMDBFCM_Coord)  {

      FreeCoordMemory    ();
      DeleteAllSelections();

      nAtoms = MMDB->nAtoms;
      atmLen = nAtoms;
      if (nAtoms>0)  {
        atom = new PAtom[atmLen];
        for (i=0;i<nAtoms;i++)  {
          if (MMDB->atom[i])  {
            atom[i] = newAtom();
            atom[i]->Copy ( MMDB->atom[i] );
            // the internal atom references are installed
            // by residue classes when they are read in
            // model->chain below
            atom[i]->SetAtomIndex ( i+1 );
          } else
            atom[i] = NULL;
        }
      }

      nModels = MMDB->nModels;
      if (nModels>0)  {
        model = new PModel[nModels];
        for (i=0;i<nModels;i++)  {
          if (MMDB->model[i])  {
            model[i] = newModel();
            model[i]->SetMMDBManager ( this,0 );
            model[i]->_copy ( MMDB->model[i] );
          } else
            model[i] = NULL;
        }
      }

      crModel = NULL;
      crChain = NULL;
      crRes   = NULL;

      if (MMDB->crModel)  {

        for (i=0;i<nModels;i++)
          if (model[i])  {
            if (model[i]->serNum==MMDB->crModel->serNum)  {
              crModel = model[i];
              break;
            }
          }

        if (crModel && crModel->chain && MMDB->crChain)
          for (i=0;i<crModel->nChains;i++)
            if (crModel->chain[i])  {
              if (!strcmp(crModel->chain[i]->chainID,
                          MMDB->crModel->chain[i]->chainID))  {
                crChain = crModel->chain[i];
                break;
              }
            }

        if (crChain && crChain->residue && MMDB->crRes)
          for (i=0;i<crChain->nResidues;i++)
            if (crChain->residue[i])  {
              if ((!strcmp(crChain->residue[i]->name,
                           MMDB->crRes->name))                       &&
                  (crChain->residue[i]->seqNum==MMDB->crRes->seqNum) &&
                  (!strcmp(crChain->residue[i]->insCode,
                           MMDB->crRes->insCode)))  {
                crRes = crChain->residue[i];
                break;
              }
            }
      }

      /*
      if ((MMDB->nSelections>0) && MMDB->Mask)  {
        nSelections = MMDB->nSelections;
        if (nSelections>0)  {
          Mask      = new PCMask [nSelections];
          SelAtom   = new PPAtom[nSelections];
          nSelAtoms = new int    [nSelections];
          for (i=0;i<nSelections;i++)  {
            Mask[i] = new CMask();
            Mask[i]->CopyMask ( MMDB->Mask[i] );
            nSelAtoms[i] = MMDB->nSelAtoms[i];
            if (nSelAtoms[i]>0)  {
              SelAtom[i] = new PAtom[nSelAtoms[i]];
              for (j=0;j<nSelAtoms[i];j++)
                SelAtom[i][j] = Atom[MMDB->SelAtom[i][j]->index];
            } else
              SelAtom[i] = NULL;
          }
        }
      }
      */

    } else if (CopyMask & (MMDBFCM_HetInfo | MMDBFCM_SecStruct |
                            MMDBFCM_Links | MMDBFCM_CisPeps |
                            MMDBFCM_ChainAnnot))  {

      for (i=0;i<MMDB->nModels;i++)
        if (MMDB->model[i])  {

          mdl = GetModel ( i+1 );
          if (!mdl)  {
            mdl = new Model( NULL,i+1 );
            AddModel ( mdl );
          }

          if (CopyMask & MMDBFCM_HetInfo)
            mdl->CopyHets ( MMDB->model[i] );
          if (CopyMask & MMDBFCM_SecStruct)
            mdl->CopySecStructure ( MMDB->model[i] );
          if (CopyMask & MMDBFCM_Links)  {
            mdl->CopyLinks  ( MMDB->model[i] );
            mdl->CopyLinkRs ( MMDB->model[i] );
          }
          if (CopyMask & MMDBFCM_CisPeps)
            mdl->CopyCisPeps ( MMDB->model[i] );
          if (CopyMask & MMDBFCM_ChainAnnot)  {
            MMDB->GetChainTable ( i+1,chain,nchains );
            for (j=0;j<nchains;j++)
              if (chain[j])  {
                chain[j]->GetChainID ( chID );
                ch = mdl->GetChain ( chID );
                if (!ch)  {
                  ch = new Chain();
                  ch->SetChainID ( chID );
                  mdl->AddChain ( ch );
                }
                ch->CopyAnnotations ( chain[j] );
              }

          }

        }

    }

    if (CopyMask & MMDBFCM_SA)  SA.Copy ( &(MMDB->SA) );
    if (CopyMask & MMDBFCM_SB)  SB.Copy ( &(MMDB->SB) );
    if (CopyMask & MMDBFCM_SC)  SC.Copy ( &(MMDB->SC) );
    if (CopyMask & MMDBFCM_Footnotes)
                         Footnote.Copy ( &(MMDB->Footnote) );

    if (CopyMask & MMDBFCM_Buffer)  {
      lcount = MMDB->lcount;
      strncpy ( S,MMDB->S,sizeof(S) );
    }

  }

  void  Manager::Delete ( word DelMask )  {
  PPModel model;
  PPChain chain;
  int      i,j,nm, nchains;

    if (DelMask & MMDBFCM_Flags)  Flags = 0;

    if (DelMask & MMDBFCM_Title)        title.Copy ( NULL );
    if (DelMask & MMDBFCM_TitleKeepBM)  title.FreeMemory ( true );
    if (DelMask & MMDBFCM_Cryst)        cryst.Copy ( NULL );

    if (DelMask & MMDBFCM_Coord)  {
      FreeCoordMemory    ();
      DeleteAllSelections();
    }

    if (DelMask & MMDBFCM_SecStruct)  {
      GetModelTable ( model,nm );
      if (model)
        for (i=0;i<nm;i++)
          if (model[i])
            model[i]->RemoveSecStructure();
    }

    if (DelMask & MMDBFCM_HetInfo)  {
      GetModelTable ( model,nm );
      if (model)
        for (i=0;i<nm;i++)
          if (model[i])
            model[i]->RemoveHetInfo();
    }

    if (DelMask & MMDBFCM_Links)  {
      GetModelTable ( model,nm );
      if (model)
        for (i=0;i<nm;i++)
          if (model[i])  {
            model[i]->RemoveLinks ();
            model[i]->RemoveLinkRs();
          }
    }

    if (DelMask & MMDBFCM_CisPeps)  {
      GetModelTable ( model,nm );
      if (model)
        for (i=0;i<nm;i++)
          if (model[i])
            model[i]->RemoveCisPeps();
    }

    if (DelMask & MMDBFCM_ChainAnnot)  {
      nm = GetNumberOfModels();
      for (i=1;i<=nm;i++)  {
        GetChainTable ( i,chain,nchains );
        if (chain)
          for (j=0;j<nchains;j++)
            if (chain[j])
              chain[j]->FreeAnnotations();
      }
    }

    if (DelMask & MMDBFCM_SA)        SA.FreeContainer();
    if (DelMask & MMDBFCM_SB)        SB.FreeContainer();
    if (DelMask & MMDBFCM_SC)        SC.FreeContainer();
    if (DelMask & MMDBFCM_Footnotes) Footnote.FreeContainer();

    if (DelMask & MMDBFCM_Buffer)  {
      lcount = 0;
      S[0]   = char(0);
    }

  }

  PTitleContainer Manager::GetRemarks()  {
    return title.GetRemarks();
  }


  PTitleContainer Manager::GetJournal()  {
    return title.GetJournal();
  }

  realtype Manager::GetResolution()  {
    return title.GetResolution();
  }

  int Manager::ParseBiomolecules()  {
    return title.ParseBiomolecules();
  }

  int Manager::GetNofBiomolecules()  {
    return title.GetNofBiomolecules();
  }

  void Manager::GetBiomolecules ( PPBiomolecule & BM, int & nBMs )  {
    title.GetBiomolecules ( BM,nBMs );
  }

  PBiomolecule Manager::GetBiomolecule ( int bmNo )  {
    return title.GetBiomolecule ( bmNo );
  }

  PManager Manager::MakeBiomolecule ( int bmNo, int modelNo ) {
  PManager M;
  PPChain      ch;
  PChain       chain;
  PModel       model;
  PBiomolecule BM;
  int           i,j,k,n,n0,nChains;

    BM = title.GetBiomolecule ( bmNo );
    if (!BM)  return NULL;

    GetChainTable ( modelNo,ch,nChains );
    if ((!ch) || (nChains<=0))  return NULL;

    n0    = 0;
    model = new Model();

    for (i=0;(i<BM->nBMAs) && (n0>=0);i++)
      if (BM->bmApply[i])  {
        for (j=0;(j<BM->bmApply[i]->nMatrices) && (n0>=0);j++)
          for (k=0;(k<BM->bmApply[i]->nChains) && (n0>=0);k++)  {
            n0 = -1;
            for (n=0;(n<nChains) && (n0<0);n++)
              if (!strcmp(ch[n]->GetChainID(),BM->bmApply[i]->chain[k]))
                n0 = n;
            if (n0>=0)  {
              chain = new Chain();
              chain->Copy ( ch[n0] );
              chain->ApplyTransform ( BM->bmApply[i]->tm[j] );
              model->AddChain ( chain );
            }
          }
      }

    if (n0>=0)  {
      M = new Manager();
      M->AddModel ( model );
      M->PDBCleanup ( PDBCLEAN_SERIAL | PDBCLEAN_INDEX );
    } else  {
      delete model;
      M = NULL;
    }

    return M;

  }


  //  -------------------  Stream functions  ----------------------


  void  Manager::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte ( &Version );
    BondManager::write ( f );
  }

  void  Manager::read ( io::RFile f )  {
  byte Version;
    f.ReadByte ( &Version );
    BondManager::read ( f );
  }


  MakeStreamFunctions(Manager)

}  // namespace mmdb
