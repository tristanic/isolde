//  $Id: mmdb_selmngr.cpp $
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
//  **** Module  :  mmdb_selmngr <implementation>
//       ~~~~~~~~~
//       Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::Manager ( MMDB atom selection manager )
//       ~~~~~~~~~
//
//   (C) E. Krissinel 2000-2013
//
//  =================================================================
//


#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "mmdb_selmngr.h"


namespace mmdb  {

  const int ANY_RES = MinInt4;

  //  ====================   SelManager   =====================

  SelManager::SelManager() : CoorManager()  {
    InitSelManager();
  }

  SelManager::SelManager ( io::RPStream Object )
                 : CoorManager(Object)  {
    InitSelManager();
  }

  SelManager::~SelManager()  {
    DeleteAllSelections();
  }

  void  SelManager::ResetManager()  {
    CoorManager::ResetManager();
    DeleteAllSelections();
    InitSelManager ();
  }

  void  SelManager::InitSelManager()  {
    nSelections = 0;     // number of selections
    mask        = NULL;  // vector of selections
    selType     = NULL;  // vector of selection types
    nSelItems   = NULL;  // numbers of selected items
    selection   = NULL;  // vector of selected items
  }


  // ------------------------  Selection  -----------------------------

  int  SelManager::NewSelection()  {
  PMask    M;
  PPMask   Mask1;
  PPMask * Selection1;
  ivector  nSelItems1;
  SELECTION_TYPE * SelType1;
  int      i,l;

    M = new Mask();
    M->NewMask ( mask,nSelections );

    i = 0;
    while (i<nSelections)
      if (!mask[i])  break;
               else  i++;

    if (i>=nSelections)  {
      l          = nSelections+10;
      Mask1      = new PMask [l];
      Selection1 = new PPMask[l];
      nSelItems1 = new int[l];
      SelType1   = new SELECTION_TYPE[l];
      for (i=0;i<nSelections;i++)  {
        Mask1     [i] = mask     [i];
        Selection1[i] = selection[i];
        nSelItems1[i] = nSelItems[i];
        SelType1  [i] = selType  [i];
      }
      for (i=nSelections;i<l;i++)  {
        Mask1     [i] = NULL;
        Selection1[i] = NULL;
        nSelItems1[i] = 0;
        SelType1  [i] = STYPE_UNDEFINED;
      }
      if (mask)      delete[] mask;
      if (selection) delete[] selection;
      if (nSelItems) delete[] nSelItems;
      if (selType)   delete[] selType;
      mask        = Mask1;
      selection   = Selection1;
      nSelItems   = nSelItems1;
      selType     = SelType1;
      i           = nSelections;
      nSelections = l;
    }

    mask[i] = M;
    if (selection[i])  delete[] selection[i];
    selection[i] = NULL;
    nSelItems[i] = 0;
    selType  [i] = STYPE_UNDEFINED;

    return i+1;

  }

  int  SelManager::GetSelType ( int selHnd )  {
  int k;
    if ((selHnd>0) && (selHnd<=nSelections))  {
      k = selHnd-1;
      if (mask[k])  return selType[k];
    }
    return STYPE_INVALID;
  }

  void  SelManager::DeleteSelection ( int selHnd )  {
  int i,k;
    if ((selHnd>0) && (selHnd<=nSelections))  {
      k = selHnd-1;
      if (mask[k])  {
        for (i=0;i<nSelItems[k];i++)
          if (selection[k][i])
            selection[k][i]->RemoveMask ( mask[k] );

        //      for (i=0;i<nAtoms;i++)
        //        if (atom[i])
        //          atom[i]->RemoveMask ( mask[k] );

        delete mask[k];
      }
      mask[k] = NULL;
      if (selection[k])  delete[] selection[k];
      selection[k] = NULL;
      nSelItems[k] = 0;
      selType  [k] = STYPE_UNDEFINED;
    }
  }


  PMask SelManager::GetSelMask ( int selHnd )  {
    if ((selHnd>0) && (selHnd<=nSelections))
         return mask[selHnd-1];
    else return NULL;
  }

  void  SelManager::DeleteAllSelections()  {
  PResidue res  ,res1;
  PChain   chain,chain1;
  PModel   model,model1;
  int      i;

    if (mask)  {
      res   = NULL;
      chain = NULL;
      model = NULL;
      if (atom)
        for (i=0;i<nAtoms;i++)
          if (atom[i])  {
            atom[i]->ClearMask();
            res1 = atom[i]->GetResidue();
            if (res1!=res)  {
              res = res1;
              res->ClearMask();
              chain1 = res->GetChain();
              if (chain1!=chain)  {
                chain = chain1;
                chain->ClearMask();
                model1 = chain->GetModel();
                if (model1!=model)  {
                  model = model1;
                  model->ClearMask();
                }
              }
            }
          }
      for (i=0;i<nSelections;i++)  {
        if (mask     [i])  delete   mask[i];
        if (selection[i])  delete[] selection[i];
      }
      delete[] mask;
      if (selection) delete[] selection;
      if (nSelItems) delete[] nSelItems;
      if (selType)   delete[] selType;
    }

    nSelections = 0;
    mask        = NULL;
    selection   = NULL;
    nSelItems   = NULL;
    selType     = NULL;

  }

  void  SelManager::SelectAtoms ( int selHnd, int iSer1, int iSer2,
                                  SELECTION_KEY sKey )  {
  //   SelectAtoms(..) selects atoms in the serial number range
  // of iSer1 to iSer2 by adding them to the set of atoms
  // marked by the given mask. If iSer1=iSer2=0 then all atoms
  // are selected. Each atom may be selected by a number of masks
  // simultaneously
  int           i,s1,s2,k,nsel;
  SELECTION_KEY sk;

    if ((selHnd<=0) || (selHnd>nSelections) || (nAtoms<=0))  return;

    k  = selHnd-1;
    sk = sKey;

    if ((selType[k]==STYPE_UNDEFINED) || (sKey==SKEY_NEW))
             selType[k] = STYPE_ATOM;
    else if (selType[k]!=STYPE_ATOM)  return;

    switch (sKey)  {
      case SKEY_NEW : for (i=0;i<nAtoms;i++)
                        if (atom[i])  atom[i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : nsel = 0;             break;
      case SKEY_XOR : nsel = nSelItems[k];  break;
      case SKEY_CLR : nsel = nSelItems[k];
                      if (nsel<=0)  return;
                    break;
      case SKEY_XAND: nsel = 0;             break;
    }

    if ((iSer1==0) && (iSer2==0))  {
      for (i=0;i<nAtoms;i++)
        if (atom[i])  {
          if (!atom[i]->Ter)
            SelectAtom ( atom[i],k,sk,nsel );
        }
    } else  {
      if (iSer1<=iSer2)  {
        s1 = iSer1;
        s2 = iSer2;
      } else  {
        s1 = iSer2;
        s2 = iSer1;
      }
      // for a very general use, we allow the serial number
      // to differ from the atom's index, although this is
      // against PDB format. Therefore we apply here the most
      // primitive and less efficient way of selection
      for (i=0;i<nAtoms;i++)
        if (atom[i])  {
          if (!atom[i]->Ter)  {
            if ((s1<=atom[i]->serNum) && (atom[i]->serNum<=s2))
              SelectAtom ( atom[i],k,sk,nsel );
            else if (sk==SKEY_AND)
              atom[i]->RemoveMask ( mask[k] );
          }
        }
    }

    MakeSelIndex ( selHnd,STYPE_ATOM,nsel );

  }


  void  SelManager::SelectAtoms ( int selHnd, ivector asn, int nsn,
                                  SELECTION_KEY selKey )  {
  //   SelectAtoms(..) selects atoms with serial numbers given in
  // vector asn[0..nsn-1].
  QuickSort     QS;
  ivector       asn1;
  int           i,k,nsn1,j,j1,j2, sn,nsel;
  SELECTION_KEY sk;

    if ((selHnd<=0) || (selHnd>nSelections) || (nAtoms<=0))  return;

    k  = selHnd-1;
    sk = selKey;

    if ((selType[k]==STYPE_UNDEFINED) ||
        (selKey==SKEY_NEW))           selType[k] = STYPE_ATOM;
    else if (selType[k]!=STYPE_ATOM)  return;

    switch (selKey)  {
      case SKEY_NEW : for (i=0;i<nAtoms;i++)
                        if (atom[i])  atom[i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : nsel = 0;             break;
      case SKEY_XOR : nsel = nSelItems[k];  break;
      case SKEY_CLR : nsel = nSelItems[k];
                      if (nsel<=0)  return;
                    break;
      case SKEY_XAND: nsel = 0;             break;
    }

    GetVectorMemory ( asn1,nsn,0 );
    for (i=0;i<nsn;i++)
      asn1[i] = asn[i];

    QS.Sort ( asn1,nsn );
    nsn1 = nsn-1;

    for (i=0;i<nAtoms;i++)
      if (atom[i])  {
        if (!atom[i]->Ter)  {
          sn = atom[i]->serNum;
          if ((asn1[0]<=sn) && (sn<=asn1[nsn1]))  {
            // binary search
            j1 = 0;
            j2 = nsn1;
            do  {
              j = (j1+j2)/2;
              if (sn<asn1[j])      j2 = j;
              else if (sn>asn1[j]) j1 = j;
                              else j1 = j2;
            } while (j1<j2-1);
            if ((sn==asn1[j]) || (sn==asn1[j1]) || (sn==asn1[j2]))
              SelectAtom ( atom[i],k,sk,nsel );
            else if (sk==SKEY_AND)
              atom[i]->RemoveMask ( mask[k] );
          } else if (sk==SKEY_AND)
            atom[i]->RemoveMask ( mask[k] );
        }
      }

    FreeVectorMemory ( asn1,0 );

    MakeSelIndex ( selHnd,STYPE_ATOM,nsel );

  }


  void  SelManager::UnselectAtoms ( int selHnd, int iSer1, int iSer2 )  {
  //   UnselectAtoms(..) clears the specified mask for atoms in
  // the serial number range of iSer1 to iSer2. If iSer1=iSer2=0
  // then all atoms are cleared of the specified mask. If selHnd
  // is set to 0, then the atoms are cleared of any mask.
  int i,s1,s2,k;

    if ((selHnd<=nSelections) && (nAtoms>0))  {

      k = selHnd-1;

      if (selType[k]==STYPE_UNDEFINED)  selType[k] = STYPE_ATOM;
      else if (selType[k]!=STYPE_ATOM)  return;

      if ((iSer1==0) && (iSer2==0))  {
        if (k<0) {
          for (i=0;i<nAtoms;i++)
            if (atom[i]) atom[i]->ClearMask();
        } else  {
          for (i=0;i<nAtoms;i++)
            if (atom[i]) atom[i]->RemoveMask ( mask[k] );
        }
      } else  {
        if (iSer1<=iSer2)  {
          s1 = iSer1;
          s2 = iSer2;
        } else  {
          s1 = iSer2;
          s2 = iSer1;
        }
        // for a very general use, we allow the serial number
        // to differ from the atom's index, although this is
        // against PDB format. Therefore we apply here the most
        // primitive and less efficient way of selection
        if (k<0)  {
          for (i=0;i<nAtoms;i++)
            if (atom[i])  {
              if ((s1<=atom[i]->serNum) && (atom[i]->serNum<=s2))
                atom[i]->ClearMask();
            }
        } else  {
          for (i=0;i<nAtoms;i++)
            if (atom[i])  {
              if ((s1<=atom[i]->serNum) && (atom[i]->serNum<=s2))
                atom[i]->RemoveMask ( mask[k] );
            }
        }
      }

      MakeSelIndex ( selHnd,STYPE_ATOM,-1 );

    }

  }


  pstr MakeList ( cpstr S )  {
  // makes the list of selecting items:
  //   1st character - special use,
  //       then each item from S embraced by commas
  pstr L;
  int  i,j;
    i = 0;
    while (S[i]==' ')  i++;
    if (S[i]!='*')  {
      // compile a searchable list
      L = new char[strlen(S)+5];
      if (S[i]=='!')  {
        L[0] = '!';
        i++;
      } else
        L[0] = ' ';
      if (FirstOccurence(S,'['))  L[1] = '"';
                    else  L[1] = ' ';
      L[2] = ',';
      j    = 3;
      while (S[i])  {
        while (S[i]==' ')  i++;
        if (S[i]=='[')  {
          while (S[i] && (S[i]!=']'))
            L[j++] = S[i++];
          L[j++] = ']';
          if (S[i]==']')  i++;
        } else
          while (S[i] && (S[i]!=' ') && (S[i]!=','))
            L[j++] = S[i++];
        while (S[i]==' ')  i++;
        L[j++] = ',';
        if (S[i]==',')  {
          i++;
          if (!S[i])  L[j++] = ',';  // blank chain ID at the end assumed
        }
      }
      if (j==3)  L[j++] = ',';
      L[j] = char(0);
    } else
      L = NULL;
    return L;
  }

  bool MatchName ( pstr L, pstr N )  {
  char M[sizeof(maxMMDBName)+5];
  int  i,j;
    if (L)  {
      i    = 0;
      M[0] = ',';
      j    = 1;
      while (N[i])
        if (N[i]==' ')  i++;
                  else  M[j++] = N[i++];
      M[j++] = ',';
      M[j]   = char(0);
      if (strstr(&(L[2]),M))  return (L[0]!='!');
      else if (L[1]!='"')     return (L[0]=='!');
      else  {
        strcpy ( M,",[" );
        strcat ( M,N    );
        strcat ( M,"]," );
        if (strstr(&(L[2]),M))  return (L[0]!='!');
                          else  return (L[0]=='!');
      }
    } else
      return true;
  }

  bool MatchCharge ( pstr L, PAtom atom )  {
  char N[100];
    if (L)  {
      if (atom->WhatIsSet & ASET_Charge)  {
        sprintf ( N,"%+2i",mround(atom->charge) );
        return MatchName ( L,N );
      } else
        return false;
    } else
      return true;
  }


  void SelManager::SelectAtom ( int selHnd, PAtom A,
                                SELECTION_KEY selKey,
                                bool makeIndex )  {
  int           i, k, nsel;
  SELECTION_KEY sk;

    if ((selHnd<=0) || (selHnd>nSelections))  return;

    k  = selHnd-1;
    sk = selKey;

    if ((selType[k]==STYPE_UNDEFINED) ||
        (selKey==SKEY_NEW))           selType[k] = STYPE_ATOM;
    else if (selType[k]!=STYPE_ATOM)  return;

    switch (selKey)  {
      case SKEY_NEW : for (i=0;i<nSelItems[k];i++)
                        if (selection[k][i])
                            selection[k][i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : if (nSelItems[k]==0)  return;
                      nsel = 0;
                    break;
      case SKEY_XOR : nsel = nSelItems[k];
                    break;
      case SKEY_CLR : nsel = nSelItems[k];
                      if (nsel<=0)  return;
                    break;
      case SKEY_XAND: nsel = 0;             break;
    }

    SelectAtom ( A,k,sk,nsel);
    if (makeIndex)  MakeSelIndex ( selHnd,STYPE_ATOM,nsel );

  }


  void SelManager::SelectResidue ( int selHnd, PResidue Res,
                                   SELECTION_TYPE sType,
                                   SELECTION_KEY  sKey,
                                   bool makeIndex )  {
  //  Selects residue Res or all its atoms depending on selType
  PPAtom        A;
  int           i, k, nsel, nat;
  SELECTION_KEY sk;

    if ((selHnd<=0) || (selHnd>nSelections))  return;

    k  = selHnd-1;
    sk = sKey;

    if ((selType[k]==STYPE_UNDEFINED) || (sKey==SKEY_NEW))
       selType[k] = sType;
    else if (selType[k]!=sType)  return;

    switch (sKey)  {
      case SKEY_NEW : for (i=0;i<nSelItems[k];i++)
                        if (selection[k][i])
                            selection[k][i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : if (nSelItems[k]==0)  return;
                      nsel = 0;
                    break;
      case SKEY_XOR : nsel = nSelItems[k];
                    break;
      case SKEY_CLR : nsel = nSelItems[k];
                      if (nsel<=0)  return;
                    break;
      case SKEY_XAND: nsel = 0;             break;
    }

    switch (sType)  {
      case STYPE_ATOM    :  Res->GetAtomTable ( A,nat );
                            for (i=0;i<nat;i++)
                              if (A[i])  {
                                if (!A[i]->Ter)
                                  SelectAtom ( A[i],k,sk,nsel);
                              }
                          break ;
      case STYPE_RESIDUE :  SelectObject ( Res,k,sk,nsel );
                          break ;
      default : ;
    }

    if (makeIndex)  MakeSelIndex ( selHnd,sType,nsel );

  }


  void SelManager::SelectChain ( int selHnd, PChain Chain,
                                 SELECTION_TYPE sType,
                                 SELECTION_KEY  sKey,
                                 bool makeIndex )  {
  //  Selects chain Chain or all its residues or atoms depending on selType
  PPAtom    A;
  PPResidue Res;
  int       i,j, k, nsel, nat,nres;
  SELECTION_KEY sk;

    if ((selHnd<=0) || (selHnd>nSelections))  return;

    k  = selHnd-1;
    sk = sKey;

    if ((selType[k]==STYPE_UNDEFINED) || (sKey==SKEY_NEW))
       selType[k] = sType;
    else if (selType[k]!=sType)  return;

    switch (sKey)  {
      case SKEY_NEW : for (i=0;i<nSelItems[k];i++)
                        if (selection[k][i])
                            selection[k][i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : if (nSelItems[k]==0)  return;
                      nsel = 0;
                    break;
      case SKEY_XOR : nsel = nSelItems[k];
                    break;
      case SKEY_CLR : nsel = nSelItems[k];
                      if (nsel<=0)  return;
                    break;
      case SKEY_XAND: nsel = 0;             break;
    }

    switch (sType)  {
      case STYPE_ATOM    :  Chain->GetResidueTable ( Res,nres );
                            for (i=0;i<nres;i++)
                              if (Res[i])  {
                                Res[i]->GetAtomTable ( A,nat );
                                for (j=0;j<nat;j++)
                                  if (A[j])  {
                                    if (!A[j]->Ter)
                                      SelectAtom ( A[j],k,sk,nsel);
                                  }
                              }
                          break ;
      case STYPE_RESIDUE :  Chain->GetResidueTable ( Res,nres );
                            for (i=0;i<nres;i++)
                              if (Res[i])
                                SelectObject ( Res[i],k,sk,nsel );
                          break ;
      case STYPE_CHAIN   :  SelectObject ( Chain,k,sk,nsel );
                          break ;
      default : ;
    }

    if (makeIndex)  MakeSelIndex ( selHnd,sType,nsel );

  }


  void SelManager::SelectModel ( int selHnd, PModel model,
                                 SELECTION_TYPE sType,
                                 SELECTION_KEY  sKey,
                                 bool makeIndex )  {
  //  Selects model or all its chains or residues or atoms depending
  // on selType
  PPAtom    A;
  PPResidue Res;
  PPChain   Chain;
  int       i,j,n, k, nsel, nat,nres,nch;
  SELECTION_KEY sk;

    if ((selHnd<=0) || (selHnd>nSelections))  return;

    k  = selHnd-1;
    sk = sKey;

    if ((selType[k]==STYPE_UNDEFINED) || (sKey==SKEY_NEW))
             selType[k] = sType;
    else if (selType[k]!=sType)  return;

    switch (sKey)  {
      case SKEY_NEW : for (i=0;i<nSelItems[k];i++)
                        if (selection[k][i])
                            selection[k][i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : if (nSelItems[k]==0)  return;
                      nsel = 0;
                    break;
      case SKEY_XOR : nsel = nSelItems[k];
                    break;
      case SKEY_CLR : nsel = nSelItems[k];
                      if (nsel<=0)  return;
                    break;
      case SKEY_XAND: nsel = 0;             break;
    }

    switch (sType)  {
      case STYPE_ATOM    :  model->GetChainTable ( Chain,nch );
                            for (i=0;i<nch;i++)
                              if (Chain[i])  {
                                Chain[i]->GetResidueTable ( Res,nres );
                                for (j=0;j<nres;j++)
                                  if (Res[j])  {
                                    Res[j]->GetAtomTable ( A,nat );
                                    for (n=0;n<nat;n++)
                                      if (A[n])  {
                                        if (!A[n]->Ter)
                                          SelectAtom ( A[n],k,sk,nsel);
                                      }
                                  }
                              }
                          break ;
      case STYPE_RESIDUE :  model->GetChainTable ( Chain,nch );
                            for (i=0;i<nch;i++)
                              if (Chain[i])  {
                                Chain[i]->GetResidueTable ( Res,nres );
                                for (j=0;j<nres;j++)
                                  if (Res[j])
                                    SelectObject ( Res[j],k,sk,nsel );
                              }
                          break ;
      case STYPE_CHAIN   :  model->GetChainTable ( Chain,nch );
                            for (i=0;i<nch;i++)
                              if (Chain[i])
                                SelectObject ( Chain[i],k,sk,nsel );
                          break ;
      case STYPE_MODEL   :  SelectObject ( model,k,sk,nsel );
                          break ;
      default : ;
    }

    if (makeIndex)  MakeSelIndex ( selHnd,sType,nsel );

  }


  int SelManager::MakeSelIndex ( int selHnd )  {
  int k;
    if ((selHnd<=0) || (selHnd>nSelections))  return 0;
    k = selHnd-1;
    if (selType[k]==STYPE_UNDEFINED)  return 0;
    MakeSelIndex ( selHnd,selType[k],-1 );
    return nSelItems[k];
  }

  void SelManager::MakeAllSelIndexes()  {
  int k;
    for (k=0;k<nSelections;k++)
      if (selType[k]!=STYPE_UNDEFINED)
        MakeSelIndex ( k+1,selType[k],-1 );
  }

  void  SelManager::SelectAtoms (
               int   selHnd,   // must be obtained from NewSelection()
               int   iModel,   // model number; iModel=0 means
                               // 'any models'
               cpstr Chains,   // may be several chains "A,B,W";
                               // "*" means 'any chain' (in model)
               int   ResNo1,   // starting residue number
               cpstr Ins1,     // starting residue insertion code;
                               // "*" means 'any code'
               int   ResNo2,   // ending residue number.
                               // ResNo1=ResNo2=ANY_RES
                               // means 'any residue number'
                               // (in chain)
               cpstr Ins2,     // ending residue insertion code
                               // "*" means 'any code'
               cpstr RNames,   // may be several residue names
                               // "ALA,GLU,CIS"; "*" means
                               // 'any residue name'
               cpstr ANames,   // may be several names "CA,CB";
                               // "*" means 'any atom' (in residue)
               cpstr Elements, // may be several element types like
                               // "H,C,O,CU"; "*" means 'any element'
               cpstr altLocs,  // may be several alternative
                               // locations 'A,B'; "*" means 'any
                               // alternative location'
               cpstr Segments, // may be several segment IDs like
                               // "S1,S2,A234"; "*" means 'any
                               // segment'
               cpstr Charges,  // may be several charges like
                               // "+1,-2,  "; "*" means 'any charge'
               realtype occ1,  // lowest occupancy
               realtype occ2,  // highest occupancy;
                               // occ1=occ2<0.0 means
                               // "any occupancy"
               realtype x0,    // reference x-point
               realtype y0,    // reference y-point
               realtype z0,    // reference z-point
               realtype d0,    // selection distance from the reference
                               // reference point; d0<=0.0 means
                               // 'any distance" and values of
                               // x0, y0 and z0 are ignored
               SELECTION_KEY selKey     // selection key
                      )  {

    Select ( selHnd,STYPE_ATOM,iModel,Chains,ResNo1,Ins1,ResNo2,Ins2,
             RNames,ANames,Elements,altLocs,Segments,Charges,
             occ1,occ2,x0,y0,z0,d0,selKey );

  }


  #define  hetIndicator '@'

  void  SelManager::Select (
               int   selHnd,   // must be obtained from NewSelection()
               SELECTION_TYPE sType,  // selection type STYPE_XXXXX
               int   iModel,   // model number; iModel=0 means
                               // 'any models'
               cpstr Chains,   // may be several chains "A,B,W";
                               // "*" means 'any chain' (in model)
               int   ResNo1,   // starting residue number
               cpstr Ins1,     // starting residue insertion code;
                               // "*" means 'any code'
               int   ResNo2,   // ending residue number.
                               // ResNo1=ResNo2=ANY_RES means 'any
                               // residue number' (in chain)
               cpstr Ins2,     // ending residue insertion code
                               // "*" means 'any code'
               cpstr RNames,   // may be several residue names
                               // "ALA,GLU,CIS"; "*" means 'any
                               // residue name'
               cpstr ANames,   // may be several names "CA,CB";"*"
                               // means 'any atom' (in residue)
               cpstr Elements, // may be several element types like
                               // "H,C,O,CU"; "*" means 'any element'
               cpstr altLocs,  // may be several alternative
                               // locations 'A,B'; "*" means 'any
                               // alternative location'
               cpstr Segments, // may be several segment IDs like
                               // "S1,S2,A234"; "*" means 'any
                               // segment'
               cpstr Charges,  // may be several charges like
                               // "+1,-2,  "; "*" means 'any charge'
               realtype occ1,  // lowest occupancy
               realtype occ2,  // highest occupancy;
                               // occ1=occ2<0.0 means
                               // "any occupancy"
               realtype x0,    // reference x-point
               realtype y0,    // reference y-point
               realtype z0,    // reference z-point
               realtype d0,    // selection distance from the reference
                               // reference point; d0<=0.0 means
                               // 'any distance" and values of
                               // x0, y0 and z0 are ignored
               SELECTION_KEY sKey     // selection key
                      )  {
  PModel   mdl;
  PChain   chain;
  PResidue res;
  PAtom    atom;
  pstr     chain_l;
  pstr     res_l;
  pstr     atom_l;
  pstr     elem_l;
  pstr     altLocs1;
  pstr     aloc_l;
  pstr     segm_l;
  pstr     charge_l;
  realtype dx,dy,dz,d02;
  int      i,j,k,n,m1,m2,c, nsel;
  bool     noRes,Occ,Dist,Sel,selAND;
  bool     modelSel,chainSel,resSel;
  SELECTION_KEY sk;

    if ((selHnd<=0) || (selHnd>nSelections) || (nAtoms<=0))  return;

    modelSel = false;

    k  = selHnd-1;
    sk = sKey;

    if ((selType[k]==STYPE_UNDEFINED) || (sKey==SKEY_NEW))
             selType[k] = sType;
    else if (selType[k]!=sType)  return;

    // if something goes wrong, sk should be assigned SKEY_OR if
    // selKey is set to SKEY_NEW or SKEY_OR below
    switch (sKey)  {
      case SKEY_NEW : for (i=0;i<nSelItems[k];i++)
                        if (selection[k][i])
                            selection[k][i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : nsel = 0;             break;
      case SKEY_XOR : nsel = nSelItems[k];  break;
      case SKEY_CLR : nsel = nSelItems[k];
                      if (nsel<=0)  return;
                    break;
      default       : return;
    }

    selAND   = (sKey==SKEY_AND);

    altLocs1 = NULL;
    if (altLocs)  {
      if (FirstOccurence(altLocs,hetIndicator))  {
        CreateCopy ( altLocs1,altLocs );
        DelSpaces  ( altLocs1 );
        aloc_l = FirstOccurence ( altLocs1,hetIndicator );
        aloc_l[0] = ' ';
        if (aloc_l[1])  aloc_l[1] = ' ';  // instead of comma
        else if (aloc_l!=altLocs1)  {
          aloc_l--;
          aloc_l[0] = ' ';
        }
        DelSpaces  ( altLocs1 );
        aloc_l = MakeList ( altLocs1 );
      } else
        aloc_l = MakeList ( altLocs );
    } else
      aloc_l = MakeList ( altLocs );

    chain_l  = MakeList ( Chains   );
    res_l    = MakeList ( RNames   );
    atom_l   = MakeList ( ANames   );
    elem_l   = MakeList ( Elements );
    segm_l   = MakeList ( Segments );
    charge_l = MakeList ( Charges  );

    //  noRes==true means no residue restrictions
    noRes = (ResNo1==ResNo2)   && (ResNo1==ANY_RES) &&
            (Ins1[0]==Ins2[0]) && (Ins1[0]=='*');

    Occ  = (occ1>=0.0) || (occ2>=0.0);
    Dist = (d0>0.0);
    d02  = d0*d0;

    m1   = iModel-1;

    if (m1>=0)
      m2 = m1+1;     // will take only this model
    else  {
      m1 = 0;        // will take
      m2 = nModels;  //   all models
    }

    if (m1>=nModels)  return;

    for (n=0;n<nModels;n++)  {
      mdl = model[n];
      if (mdl)  {  // check for safety
        if ((m1<=n) && (n<m2))  {
          modelSel = false; // will be true on any selection in the model
          for (c=0;c<mdl->nChains;c++)  {
            chain = mdl->chain[c];
            if (chain)  {   // again check for safety
              if (MatchName(chain_l,chain->chainID))  {
                // the chain has to be taken
                i = 0;
                if (!noRes)
                  while (i<chain->nResidues)  {
                    res = chain->residue[i];
                    if (res)  {
                      if ((res->seqNum==ResNo1) &&
                          MatchName(res_l,res->name) &&
                          ((Ins1[0]=='*') ||
                           (!strcmp(res->insCode,Ins1))))
                        break;
                      else if (selAND)  {
                        if (sType==STYPE_ATOM)
                          res->UnmaskAtoms ( mask[k] );
                        else if (sType==STYPE_RESIDUE)
                          res->RemoveMask ( mask[k] );
                      }
                    }
                    i++;
                  }
                while (i<chain->nResidues)  {
                  res = chain->residue[i];
                  if (res)  {
                    resSel = false; // will be true on 1st sel-n in the res-e
                    if (MatchName(res_l,res->name))  {
                      for (j=0;j<res->nAtoms;j++)  {
                        atom = res->atom[j];
                        if (atom)  {
                          if ((!atom->Ter)                      &&
                              MatchName(atom_l  ,atom->name   ) &&
                              MatchName(elem_l  ,atom->element) &&
                              MatchName(aloc_l  ,atom->altLoc ) &&
                              MatchName(segm_l  ,atom->segID  ) &&
                              MatchCharge(charge_l,atom       ) &&
                              ((!altLocs1) || atom->Het))  {
                            Sel = true;
                            if (Occ)
                              Sel = ((occ1<=atom->occupancy) &&
                                     (atom->occupancy<=occ2));
                            if (Dist)  {
                              dx  = atom->x - x0;
                              dy  = atom->y - y0;
                              dz  = atom->z - z0;
                              Sel = Sel && ((dx*dx+dy*dy+dz*dz)<=d02);
                            }
                          } else
                            Sel = false;
                          if (Sel)  {
                            SelectObject ( sType,atom,k,sk,nsel );
                            resSel   = true;
                            chainSel = true;
                            modelSel = true;
                          } else if (selAND && (sType==STYPE_ATOM))
                            atom->RemoveMask ( mask[k] );
                        }
                        if (resSel && (sType!=STYPE_ATOM))  break;
                      }
                    } else if (selAND && (sType==STYPE_ATOM))
                        res->UnmaskAtoms ( mask[k] );
                    if ((!resSel) && selAND && (sType==STYPE_RESIDUE))
                      res->RemoveMask ( mask[k] );
                    if (chainSel && (sType>STYPE_RESIDUE))  break;
                    if (!noRes)  {
                      if ((res->seqNum==ResNo2) &&
                          ((Ins2[0]=='*') || (!strcmp(res->insCode,Ins2)))
                         )  break;
                    }
                  }
                  i++;
                }
                if (selAND)  {
                  if (sType==STYPE_ATOM)
                    while (i<chain->nResidues)  {
                      res = chain->residue[i];
                      if (res)  res->UnmaskAtoms ( mask[k] );
                      i++;
                    }
                  if (sType==STYPE_RESIDUE)
                    while (i<chain->nResidues)  {
                      res = chain->residue[i];
                      if (res)  res->RemoveMask ( mask[k] );
                      i++;
                    }
                }
              } else if (selAND)
                switch (sType)  {
                  case STYPE_ATOM    : chain->UnmaskAtoms    ( mask[k] ); break;
                  case STYPE_RESIDUE : chain->UnmaskResidues ( mask[k] ); break;
                  case STYPE_CHAIN   : chain->RemoveMask     ( mask[k] ); break;
                  default            : ;
                }
              if ((!chainSel) && selAND && (sType==STYPE_CHAIN))
                chain->RemoveMask ( mask[k] );
              if (modelSel && (sType>STYPE_CHAIN))  break;
            }
          }
        } else if (selAND)
          switch (sType)  {
            case STYPE_ATOM    : mdl->UnmaskAtoms    ( mask[k] ); break;
            case STYPE_RESIDUE : mdl->UnmaskResidues ( mask[k] ); break;
            case STYPE_CHAIN   : mdl->UnmaskChains   ( mask[k] ); break;
            default            : ;
          }
        if ((!modelSel) && selAND && (sType==STYPE_MODEL))
          mdl->RemoveMask ( mask[k] );
      }
    }

    // release dynamic memory
    if (chain_l)  delete[] chain_l;
    if (res_l)    delete[] res_l;
    if (atom_l)   delete[] atom_l;
    if (elem_l)   delete[] elem_l;
    if (altLocs1) delete[] altLocs1;
    if (aloc_l)   delete[] aloc_l;
    if (segm_l)   delete[] segm_l;
    if (charge_l) delete[] charge_l;

    MakeSelIndex ( selHnd,STYPE_ATOM,nsel );

  }


  void  SelManager::SelectAtoms (
               int   selHnd,   // must be obtained from NewSelection()
               int   iModel,   // model number; iModel=0 means
                               // 'any models'
               cpstr Chains,   // may be several chains "A,B,W";
                               // "*" means 'any chain' (in model)
               int   ResNo1,   // starting residue number
               cpstr Ins1,     // starting residue insertion code;
                               // "*" means 'any code'
               int   ResNo2,   // ending residue number.
                               // ResNo1=ResNo2=ANY_RES means 'any
                               // residue number' (in chain)
               cpstr Ins2,     // ending residue insertion code
                               // "*" means 'any code'
               cpstr RNames,   // may be several residue names
                               // "ALA,GLU,CIS"; "*" means 'any
                               // residue name'
               cpstr ANames,   // may be several names "CA,CB"; "*"
                               // means 'any atom' (in residue)
               cpstr Elements, // may be several element types like
                               // "H,C,O,CU"; "*" means 'any
                               // element'
               cpstr altLocs,  // may be several alternative
                               // locations 'A,B'; "*" means 'any
                               // alternative location'
               SELECTION_KEY sKey    // selection key
                      )  {
    Select ( selHnd,STYPE_ATOM,iModel,Chains,ResNo1,Ins1,ResNo2,Ins2,
             RNames,ANames,Elements,altLocs,sKey );
  }


  int  SelManager::Select (
               int           selHnd, // must be obtained from NewSelection()
               SELECTION_TYPE sType, // selection type STYPE_XXXXX
               cpstr          CID,   // coordinate ID
               SELECTION_KEY  sKey   // selection key
                          )  {
  InsCode insCode1,insCode2;
  pstr    RNames;
  pstr    ANames;
  pstr    Elements;
  pstr    altLocs;
  pstr    Chains;
  int     seqNum1 ,seqNum2;
  int     iModel,l,RC;

    l = IMax(10,strlen(CID))+1;
    Chains   = new char[l];
    RNames   = new char[l];
    ANames   = new char[l];
    Elements = new char[l];
    altLocs  = new char[l];

    if (strcmp(CID,"-all"))  {
      RC = ParseSelectionPath ( CID,iModel,Chains,seqNum1,insCode1,
                                seqNum2,insCode2,RNames,ANames,
                                Elements,altLocs );
    } else  {
      iModel = 0;
      strcpy ( Chains,"*" );
      seqNum1 = ANY_RES;
      seqNum2 = ANY_RES;
      strcpy ( insCode1,"*" );
      strcpy ( insCode2,"*" );
      strcpy ( RNames  ,"*" );
      strcpy ( ANames  ,"*" );
      strcpy ( Elements,"*" );
      strcpy ( altLocs ,""  ); // only main conformation by default
      RC = 0;
    }

    if (!RC)  {
      Select ( selHnd,sType,iModel,Chains,seqNum1,insCode1,
               seqNum2,insCode2,RNames,ANames,Elements,altLocs,sKey );
      RC = 0;
    }

    delete[] Chains;
    delete[] RNames;
    delete[] ANames;
    delete[] Elements;
    delete[] altLocs;

    return RC;

  }

  void  SelManager::Select (
               int   selHnd,   // must be obtained from NewSelection()
               SELECTION_TYPE sType,  // selection type STYPE_XXXXX
               int   iModel,   // model number; iModel=0 means
                               // 'any model'
               cpstr Chains,   // may be several chains "A,B,W";
                               // "*" means 'any chain' (in model)
               int   ResNo1,   // starting residue number
               cpstr Ins1,     // starting residue insertion code;
                               // "*" means 'any code'
               int   ResNo2,   // ending residue number.
                               // ResNo1=ResNo2=ANY_RES means 'any
                               // residue number' (in chain)
               cpstr Ins2,     // ending residue insertion code
                               // "*" means 'any code'
               cpstr RNames,   // may be several residue names
                               // "ALA,GLU,CIS"; "*" means 'any
                               // residue name'
               cpstr ANames,   // may be several names "CA,CB"; "*"
                               // means 'any atom' (in residue)
               cpstr Elements, // may be several element types like
                               // "H,C,O,CU"; "*" means 'any element'
               cpstr altLocs,  // may be several alternative
                               // locations 'A,B'; "*" means 'any
                               // alternative location'
               SELECTION_KEY sKey    // selection key
                      )  {
  PModel   mdl;
  PChain   chain;
  PResidue res;
  PAtom    atom;
  pstr     chain_l;
  pstr     res_l;
  pstr     atom_l;
  pstr     elem_l;
  pstr     altLocs1;
  pstr     aloc_l;
  int      i,j,k,n,m1,m2,c, nsel;
  bool     noRes,modelSel,chainSel,resSel,selAND;
  SELECTION_KEY sk;

    if ((selHnd<=0) || (selHnd>nSelections) || (nAtoms<=0))  return;

    modelSel = false;

    k  = selHnd-1;
    sk = sKey;

    if ((selType[k]==STYPE_UNDEFINED) || (sKey==SKEY_NEW))
             selType[k] = sType;
    else if (selType[k]!=sType)  return;

    // if something goes wrong, sk should be assigned SKEY_OR if
    // selKey is set to SKEY_NEW or SKEY_OR below
    switch (sKey)  {
      case SKEY_NEW : for (i=0;i<nSelItems[k];i++)
                        if (selection[k][i])
                            selection[k][i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : nsel = 0;             break;
      case SKEY_XOR : nsel = nSelItems[k];  break;
      case SKEY_CLR : nsel = nSelItems[k];
                      if (nsel<=0)  return;
                    break;
      default       : return;
    }

    selAND  = (sKey==SKEY_AND);

    altLocs1 = NULL;
    if (altLocs)  {
      if (FirstOccurence(altLocs,hetIndicator))  {
        CreateCopy ( altLocs1,altLocs );
        DelSpaces  ( altLocs1 );
        aloc_l = FirstOccurence ( altLocs1,hetIndicator );
        aloc_l[0] = ' ';
        if (aloc_l[1])  aloc_l[1] = ' ';  // instead of comma
        else if (aloc_l!=altLocs1)  {
          aloc_l--;
          aloc_l[0] = ' ';
        }
        DelSpaces  ( altLocs1 );
        aloc_l = MakeList ( altLocs1 );
      } else
        aloc_l = MakeList ( altLocs );
    } else
      aloc_l = MakeList ( altLocs );

    chain_l = MakeList ( Chains   );
    res_l   = MakeList ( RNames   );
    atom_l  = MakeList ( ANames   );
    elem_l  = MakeList ( Elements );

    //  noRes==true means no residue restrictions
    noRes   = (ResNo1==ResNo2) && (ResNo1==ANY_RES) &&
              (Ins1[0]=='*')   && (Ins2[0]=='*');

    m1      = iModel-1;
    if (m1>=0)
      m2 = m1+1;     // will take only this model
    else  {
      m1 = 0;        // will take
      m2 = nModels;  //   all models
    }

    if (m1>=nModels)  return;

    for (n=0;n<nModels;n++)  {
      mdl = model[n];
      if (mdl)  {  // check for safety
        if ((m1<=n) && (n<m2))  {
          modelSel = false; // will be true on any selection in the model
          for (c=0;c<mdl->nChains;c++)  {
            chain = mdl->chain[c];
            if (chain)  {  // again check for safety
              chainSel = false; // will be true on 1st sel-n in the chain
              if (MatchName(chain_l,chain->chainID))  {
                // the chain is to be taken
                i = 0;
                if (!noRes)  // skip "leading" residues
                  while (i<chain->nResidues)  {
                    res = chain->residue[i];
                    if (res)  {
                      if ((res->seqNum==ResNo1) &&
                          MatchName(res_l,res->name) &&
                          ((Ins1[0]=='*') ||
                           (!strcmp(res->insCode,Ins1))))
                        break;
                      else if (selAND)  {
                        if (sType==STYPE_ATOM)
                          res->UnmaskAtoms ( mask[k] );
                        else if (sType==STYPE_RESIDUE)
                          res->RemoveMask ( mask[k] );
                      }
                    }
                    i++;
                  }
                while (i<chain->nResidues)  {
                  res = chain->residue[i];
                  i++;
                  if (res)  {
                    resSel = false; // will be true on 1st selection
                                    // in the residue
                    if (MatchName(res_l,res->name))  {
                      for (j=0;j<res->nAtoms;j++)  {
                        atom = res->atom[j];
                        if (atom)  {
                          if ((!atom->Ter)                    &&
                              MatchName(atom_l,atom->name   ) &&
                              MatchName(elem_l,atom->element) &&
                              MatchName(aloc_l,atom->altLoc ) &&
                              ((!altLocs1) || atom->Het))  {
                            SelectObject ( sType,atom,k,sk,nsel );
                            resSel   = true;
                            chainSel = true;
                            modelSel = true;
                          } else if (selAND && (sType==STYPE_ATOM))
                            atom->RemoveMask ( mask[k] );
                        }
                        if (resSel && (sType!=STYPE_ATOM))  break;
                      }
                    } else if (selAND && (sType==STYPE_ATOM))
                        res->UnmaskAtoms ( mask[k] );
                    if ((!resSel) && selAND && (sType==STYPE_RESIDUE))
                      res->RemoveMask ( mask[k] );
                    if (chainSel && (sType>STYPE_RESIDUE))  break;
                    if (!noRes)  {
                      if ((res->seqNum==ResNo2) &&
                          ((Ins2[0]=='*') || (!strcmp(res->insCode,Ins2)))
                         ) break;
                    }
                  }
                }
                if (selAND)  {
                  if (sType==STYPE_ATOM)
                    while (i<chain->nResidues)  {
                      res = chain->residue[i];
                      if (res)  res->UnmaskAtoms ( mask[k] );
                      i++;
                    }
                  if (sType==STYPE_RESIDUE)
                    while (i<chain->nResidues)  {
                      res = chain->residue[i];
                      if (res)  res->RemoveMask ( mask[k] );
                      i++;
                    }
                }
              } else if (selAND)
                switch (sType)  {
                  case STYPE_ATOM    : chain->UnmaskAtoms    ( mask[k] ); break;
                  case STYPE_RESIDUE : chain->UnmaskResidues ( mask[k] ); break;
                  default            : ;
                }
              if ((!chainSel) && selAND && (sType==STYPE_CHAIN))
                chain->RemoveMask ( mask[k] );
              if (modelSel && (sType>STYPE_CHAIN))  break;
            }
          }
        } else if (selAND)
          switch (sType)  {
            case STYPE_ATOM    : mdl->UnmaskAtoms    ( mask[k] ); break;
            case STYPE_RESIDUE : mdl->UnmaskResidues ( mask[k] ); break;
            case STYPE_CHAIN   : mdl->UnmaskChains   ( mask[k] ); break;
            default            : ;
          }
        if ((!modelSel) && selAND && (sType==STYPE_MODEL))
          mdl->RemoveMask ( mask[k] );
      }
    }

    // release dynamic memory
    if (chain_l)  delete[] chain_l;
    if (res_l)    delete[] res_l;
    if (atom_l)   delete[] atom_l;
    if (elem_l)   delete[] elem_l;
    if (altLocs1) delete[] altLocs1;
    if (aloc_l)   delete[] aloc_l;

    MakeSelIndex ( selHnd,sType,nsel );

  }


  void SelManager::Select (
               int          selHnd1, // destination, must be obtained
                                     //   from NewSelection()
               SELECTION_TYPE sType, // selection type STYPE_XXXXX
               int          selHnd2, // source, must be obtained from
                                     // NewSelection() and have been
                                     //   used for selection
               SELECTION_KEY sKey    // selection key
                               )  {
  //  SKEY_XOR works only downward the hierarchy!
  PAtom    atom;
  PResidue res;
  PChain   chain;
  PModel   model;
  int      k1,k2,i,j,l,n,nsel;
  SELECTION_KEY sk;

    if ((selHnd1<=0) || (selHnd1>nSelections) ||
        (selHnd2<=0) || (selHnd2>nSelections) || (nAtoms<=0))  return;

    k1 = selHnd1-1;
    k2 = selHnd2-1;
    sk = sKey;

    if ((selType[k1]==STYPE_UNDEFINED) || (sKey==SKEY_NEW))
             selType[k1] = sType;
    else if (selType[k1]!=sType)  return;

    if (selType[k2]==STYPE_UNDEFINED)  return;

    switch (sKey)  {
      case SKEY_NEW : for (i=0;i<nSelItems[k1];i++)
                        if (selection[k1][i])
                            selection[k1][i]->RemoveMask ( mask[k1] );
                      nSelItems[k1] = 0;
                      sk   = SKEY_OR;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k1]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k1];
                    break;
      case SKEY_AND : if (nSelItems[k1]==0)  return;
                      sk   = SKEY_XAND;
                      nsel = 0;
                    break;
      case SKEY_XOR : nsel = nSelItems[k1];  break;
      case SKEY_CLR : nsel = nSelItems[k1];
                      if (nsel<=0)  return;
                    break;
      default       : return;
    }


    switch (selType[k2])  {

      case STYPE_ATOM    :
          for (i=0;i<nSelItems[k2];i++)  {
            atom = (PAtom)selection[k2][i];
            if (atom)  {
              if (!atom->Ter)
                SelectObject ( sType,atom,k1,sk,nsel );
            }
          }
        break;

      case STYPE_RESIDUE :
          for (i=0;i<nSelItems[k2];i++)  {
            res = (PResidue)selection[k2][i];
            if (res)
              switch (sType)  {
                case STYPE_ATOM  : for (j=0;j<res->nAtoms;j++)  {
                                     atom = res->atom[j];
                                     if (atom)  {
                                       if (!atom->Ter)
                                         SelectObject (atom,k1,sk,nsel);
                                     }
                                   }
                                 break;
                case STYPE_RESIDUE : //if (res->chain)
                                     SelectObject ( res,k1,sk,nsel );
                                 break;
                case STYPE_CHAIN : if (res->chain)
                                     SelectObject ( res->chain,k1,
                                                    sk,nsel );
                                 break;
                case STYPE_MODEL : if (res->chain)  {
                                     if (res->chain->model)
                                       SelectObject ( res->chain->model,
                                                      k1,sk,nsel );
                                   }
                default          : ;
              }
          }
        break;

      case STYPE_CHAIN   :
          for (i=0;i<nSelItems[k2];i++)  {
            chain = (PChain)selection[k2][i];
            if (chain)
              switch (sType)  {
                case STYPE_ATOM    : for (j=0;j<chain->nResidues;j++)  {
                                       res = chain->residue[j];
                                       if (res)
                                         for (l=0;l<res->nAtoms;l++)  {
                                           atom = res->atom[l];
                                           if (atom)  {
                                             if (!atom->Ter)
                                               SelectObject ( atom,k1,
                                                               sk,nsel );
                                           }
                                         }
                                     }
                                 break;
                case STYPE_RESIDUE : for (j=0;j<chain->nResidues;j++)  {
                                       res = chain->residue[j];
                                       if (res)
                                         SelectObject ( res,k1,sk,nsel );
                                     }
                                 break;
                case STYPE_CHAIN   : //if (chain->model)
                                       SelectObject ( chain,k1,sk,nsel );
                                 break;
                case STYPE_MODEL   : if (chain->model)
                                       SelectObject ( chain->model,k1,
                                                      sk,nsel );
                default            : ;
              }
          }
        break;

      case STYPE_MODEL   :
          for (i=0;i<nSelItems[k2];i++)  {
            model = (PModel)selection[k2][i];
            if (model)
              switch (sType)  {
                case STYPE_ATOM    :
                          for (j=0;j<model->nChains;j++)  {
                            chain = model->chain[j];
                            if (chain)
                              for (l=0;l<chain->nResidues;l++) {
                                res = chain->residue[l];
                                if (res)
                                  for (n=0;n<res->nAtoms;n++)  {
                                    atom = res->atom[n];
                                    if (atom)  {
                                      if (!atom->Ter)
                                        SelectObject ( atom,k1,sk,nsel );
                                    }
                                  }
                              }
                          }
                        break;
                case STYPE_RESIDUE :
                          for (j=0;j<model->nChains;j++)  {
                            chain = model->chain[j];
                            if (chain)
                              for (l=0;l<chain->nResidues;l++)  {
                                res = chain->residue[j];
                                if (res)
                                  SelectObject ( res,k1,sk,nsel );
                              }
                          }
                        break;
                case STYPE_CHAIN   : for (j=0;j<model->nChains;j++)  {
                                       chain = model->chain[j];
                                       if (chain)
                                         SelectObject (chain,k1,sk,nsel);
                                     }
                                 break;
                case STYPE_MODEL   : SelectObject ( model,k1,sk,nsel );
                default            : ;
              }
          }
        break;

      default : ;

    }

    if (sKey==SKEY_AND)
      for (i=0;i<nSelItems[k1];i++)
        if (selection[k1][i])
          selection[k1][i]->XadMask ( mask[k1] );

    MakeSelIndex ( selHnd1,sType,nsel );

  }

  void  SelManager::SelectProperty (
                    int  selHnd, // must be obtained from NewSelection()
                    SELECTION_PROPERTY propKey, // property key
                    SELECTION_TYPE     sType, // selection type STYPE_XXXXX
                    SELECTION_KEY      sKey  // selection key
                  )  {
  PModel   mdl;
  PChain   chain;
  PResidue res;
  int      i,k,selHnd1,nsel, m,c,r;
  bool     doSelect;
  SELECTION_KEY sk;

    if ((selHnd<=0) || (selHnd>nSelections) || (nAtoms<=0))  return;

    k  = selHnd-1;
    if ((selType[k]==STYPE_UNDEFINED) || (sKey==SKEY_NEW))
            selType[k] = sType;
    else if (selType[k]!=sType)  return;

    if (sType!=STYPE_RESIDUE)  {
      selHnd1 = NewSelection();
      if ((sKey==SKEY_AND) || (sKey==SKEY_CLR))
        Select ( selHnd1,STYPE_RESIDUE,selHnd,SKEY_NEW );
    } else
      selHnd1 = selHnd;

    k          = selHnd1-1;
    selType[k] = STYPE_RESIDUE;
    sk         = sKey;

    switch (sKey)  {
      case SKEY_NEW : for (i=0;i<nSelItems[k];i++)
                        if (selection[k][i])
                            selection[k][i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      sk   = SKEY_OR;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : if (nSelItems[k]==0)  return;
                      sk   = SKEY_XAND;
                      nsel = 0;
                    break;
      case SKEY_XOR : nsel = nSelItems[k];  break;
      case SKEY_CLR : nsel = nSelItems[k];
                      if (nsel<=0)  return;
                    break;
      default       : return;
    }

    if ((sKey==SKEY_AND) || (sKey==SKEY_CLR))  {

      for (i=0;i<nSelItems[k];i++)  {
        res = (PResidue)selection[k][i];
        if (res)  {
          switch (propKey)  {
            case SELPROP_Solvent    : doSelect = res->isSolvent();
                                    break;
            case SELPROP_Aminoacid  : doSelect = res->isAminoacid();
                                    break;
            case SELPROP_Nucleotide : doSelect = res->isNucleotide();
                                    break;
            case SELPROP_Sugar      : doSelect = res->isSugar();
                                    break;
            case SELPROP_ModRes     : doSelect = res->isModRes();
                                    break;
            default : doSelect = false;
          }
          if (doSelect)  SelectObject ( res,k,sk,nsel );
        }
      }

      if (sKey==SKEY_AND)
        for (i=0;i<nSelItems[k];i++)
          if (selection[k][i])
            selection[k][i]->XadMask ( mask[k] );

    } else  {

      for (m=0;m<nModels;m++)  {
        mdl = model[m];
        if (mdl)  {
          for (c=0;c<mdl->nChains;c++)  {
            chain = mdl->chain[c];
            if (chain)  {
              for (r=0;r<chain->nResidues;r++)  {
                res = chain->residue[r];
                if (res)  {
                  switch (propKey)  {
                    case SELPROP_Solvent    : doSelect = res->isSolvent();
                                           break;
                    case SELPROP_Aminoacid  : doSelect = res->isAminoacid();
                                           break;
                    case SELPROP_Nucleotide : doSelect = res->isNucleotide();
                                           break;
                    case SELPROP_Sugar      : doSelect = res->isSugar();
                                           break;
                    case SELPROP_ModRes     : doSelect = res->isModRes();
                                           break;
                    default : doSelect = false;
                  }
                  if (doSelect)  SelectObject ( res,k,sk,nsel );
                }
              }
            }
          }
        }
      }

    }


    MakeSelIndex ( selHnd1,STYPE_RESIDUE,nsel );

    if (sType!=STYPE_RESIDUE)  {
      Select ( selHnd,sType,selHnd1,SKEY_NEW );
      DeleteSelection ( selHnd1 );
    }

  }


  void SelManager::SelectUDD (
               int    selHnd, // must be obtained from NewSelection()
               SELECTION_TYPE sType, // selection type STYPE_XXXXX
               int UDDhandle, // UDD handle
               int    selMin, // lower selection boundary
               int    selMax, // upper selection boundary
               SELECTION_KEY sKey  // selection key
             )  {
  PModel   mdl;
  PChain   chain;
  PResidue res;
  PAtom    atom;
  int      i,k,nsel,iudd, n,c,r,a;
  bool     selAND;
  SELECTION_KEY sk;

    k  = selHnd-1;
    sk = sKey;

    if ((selType[k]==STYPE_UNDEFINED) || (sKey==SKEY_NEW))
            selType[k] = sType;
    else if (selType[k]!=sType)  return;


    if ((selHnd<=0) || (selHnd>nSelections))  return;

    switch (sType)  {
      case STYPE_ATOM    : if ((UDDhandle & UDRF_ATOM)==0)    return;
                        break;
      case STYPE_RESIDUE : if ((UDDhandle & UDRF_RESIDUE)==0) return;
                        break;
      case STYPE_CHAIN   : if ((UDDhandle & UDRF_CHAIN)==0)   return;
                        break;
      case STYPE_MODEL   : if ((UDDhandle & UDRF_MODEL)==0)   return;
                        break;
      default            : return;
    }


    // if something goes wrong, sk should be assigned SKEY_OR if
    // selKey is set to SKEY_NEW or SKEY_OR below
    switch (sKey)  {
      case SKEY_NEW : for (i=0;i<nSelItems[k];i++)
                        if (selection[k][i])
                            selection[k][i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : if (nSelItems[k]==0)  return;
                      nsel = 0;
                    break;
      case SKEY_XOR : nsel = nSelItems[k];  break;
      case SKEY_CLR : nsel = nSelItems[k];
                      if (nsel<=0)  return;
                    break;
      default       : return;
    }

    selAND = (sKey==SKEY_AND);


    for (n=0;n<nModels;n++)  {

      mdl = model[n];
      if (mdl)  {  // check for safety

        if (sType==STYPE_MODEL)  {

          mdl->getUDData ( UDDhandle,iudd );
          if ((selMin<=iudd) && (iudd<=selMax))
            SelectObject ( mdl,k,sk,nsel );
          else if (selAND)
            mdl->RemoveMask ( mask[k] );

        } else  {

          for (c=0;c<mdl->nChains;c++)  {

            chain = mdl->chain[c];
            if (chain)  {   // again check for safety

              if (sType==STYPE_CHAIN)  {
                chain->getUDData ( UDDhandle,iudd );
                if ((selMin<=iudd) && (iudd<=selMax))
                  SelectObject ( chain,k,sk,nsel );
                else if (selAND)
                  chain->RemoveMask ( mask[k] );

              } else  {

                for (r=0;r<chain->nResidues;r++)  {

                  res = chain->residue[r];
                  if (res)  {

                    if (sType==STYPE_RESIDUE)  {
                      res->getUDData ( UDDhandle,iudd );
                      if ((selMin<=iudd) && (iudd<=selMax))
                        SelectObject ( res,k,sk,nsel );
                      else if (selAND)
                        res->RemoveMask ( mask[k] );

                    } else  {

                      for (a=0;a<res->nAtoms;a++)  {
                        atom = res->atom[a];
                        if (atom)  {
                          if (!atom->Ter)  {
                            atom->getUDData ( UDDhandle,iudd );
                            if ((selMin<=iudd) && (iudd<=selMax))
                              SelectObject ( atom,k,sk,nsel );
                            else if (selAND)
                              atom->RemoveMask ( mask[k] );
                          }
                        }

                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    MakeSelIndex ( selHnd,sType,nsel );

  }


  void SelManager::SelectUDD (
               int      selHnd, // must be obtained from NewSelection()
               SELECTION_TYPE sType, // selection type STYPE_XXXXX
               int   UDDhandle, // UDD handle
               realtype selMin, // lower selection boundary
               realtype selMax, // upper selection boundary
               SELECTION_KEY sKey  // selection key
             )  {
  PModel   mdl;
  PChain   chain;
  PResidue res;
  PAtom    atom;
  realtype rudd;
  int      i,k,nsel, n,c,r,a;
  bool     selAND;
  SELECTION_KEY sk;

    k  = selHnd-1;
    sk = sKey;

    if ((selType[k]==STYPE_UNDEFINED) || (sKey==SKEY_NEW))
            selType[k] = sType;
    else if (selType[k]!=sType)  return;


    if ((selHnd<=0) || (selHnd>nSelections))  return;

    switch (sType)  {
      case STYPE_ATOM    : if ((UDDhandle & UDRF_ATOM)==0)    return;
                        break;
      case STYPE_RESIDUE : if ((UDDhandle & UDRF_RESIDUE)==0) return;
                        break;
      case STYPE_CHAIN   : if ((UDDhandle & UDRF_CHAIN)==0)   return;
                        break;
      case STYPE_MODEL   : if ((UDDhandle & UDRF_MODEL)==0)   return;
                        break;
      default            : return;
    }


    // if something goes wrong, sk should be assigned SKEY_OR if
    // selKey is set to SKEY_NEW or SKEY_OR below
    switch (sKey)  {
      case SKEY_NEW : for (i=0;i<nSelItems[k];i++)
                        if (selection[k][i])
                            selection[k][i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : if (nSelItems[k]==0)  return;
                      nsel = 0;
                    break;
      case SKEY_XOR : nsel = nSelItems[k];  break;
      case SKEY_CLR : nsel = nSelItems[k];
                      if (nsel<=0)  return;
                    break;
      default       : return;
    }

    selAND = (sKey==SKEY_AND);


    for (n=0;n<nModels;n++)  {

      mdl = model[n];
      if (mdl)  {  // check for safety

        if (sType==STYPE_MODEL)  {

          mdl->getUDData ( UDDhandle,rudd );
          if ((selMin<=rudd) && (rudd<=selMax))
            SelectObject ( mdl,k,sk,nsel );
          else if (selAND)
            mdl->RemoveMask ( mask[k] );

        } else  {

          for (c=0;c<mdl->nChains;c++)  {

            chain = mdl->chain[c];
            if (chain)  {   // again check for safety

              if (sType==STYPE_CHAIN)  {
                chain->getUDData ( UDDhandle,rudd );
                if ((selMin<=rudd) && (rudd<=selMax))
                  SelectObject ( chain,k,sk,nsel );
                else if (selAND)
                  chain->RemoveMask ( mask[k] );

              } else  {

                for (r=0;r<chain->nResidues;r++)  {

                  res = chain->residue[r];
                  if (res)  {

                    if (sType==STYPE_RESIDUE)  {
                      res->getUDData ( UDDhandle,rudd );
                      if ((selMin<=rudd) && (rudd<=selMax))
                        SelectObject ( res,k,sk,nsel );
                      else if (selAND)
                        res->RemoveMask ( mask[k] );

                    } else  {

                      for (a=0;a<res->nAtoms;a++)  {
                        atom = res->atom[a];
                        if (atom)  {
                          if (!atom->Ter)  {
                            atom->getUDData ( UDDhandle,rudd );
                            if ((selMin<=rudd) && (rudd<=selMax))
                              SelectObject ( atom,k,sk,nsel );
                            else if (selAND)
                              atom->RemoveMask ( mask[k] );
                          }
                        }

                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    MakeSelIndex ( selHnd,sType,nsel );

  }


  bool selSUDD ( cpstr sudd, cpstr selStr, int cmpRule, int ssLen )  {
    if (!sudd)  return false;
    switch (cmpRule)  {
      case UDSCR_LT        : return (strcmp(sudd,selStr)<0);
      case UDSCR_LE        : return (strcmp(sudd,selStr)<=0);
      case UDSCR_EQ        : return (strcmp(sudd,selStr)==0);
      case UDSCR_NE        : return (strcmp(sudd,selStr)!=0);
      case UDSCR_GE        : return (strcmp(sudd,selStr)>=0);
      case UDSCR_GT        : return (strcmp(sudd,selStr)>=0);
      case UDSCR_LTcase    : return (strcasecmp(sudd,selStr)<0);
      case UDSCR_LEcase    : return (strcasecmp(sudd,selStr)<=0);
      case UDSCR_EQcase    : return (strcasecmp(sudd,selStr)==0);
      case UDSCR_NEcase    : return (strcasecmp(sudd,selStr)!=0);
      case UDSCR_GEcase    : return (strcasecmp(sudd,selStr)>=0);
      case UDSCR_GTcase    : return (strcasecmp(sudd,selStr)>=0);
      case UDSCR_LTn       : return (strncmp(sudd,selStr,ssLen)<0);
      case UDSCR_LEn       : return (strncmp(sudd,selStr,ssLen)<=0);
      case UDSCR_EQn       : return (strncmp(sudd,selStr,ssLen)==0);
      case UDSCR_NEn       : return (strncmp(sudd,selStr,ssLen)!=0);
      case UDSCR_GEn       : return (strncmp(sudd,selStr,ssLen)>=0);
      case UDSCR_GTn       : return (strncmp(sudd,selStr,ssLen)>=0);
      case UDSCR_LTncase   : return (strncasecmp(sudd,selStr,ssLen)<0);
      case UDSCR_LEncase   : return (strncasecmp(sudd,selStr,ssLen)<=0);
      case UDSCR_EQncase   : return (strncasecmp(sudd,selStr,ssLen)==0);
      case UDSCR_NEncase   : return (strncasecmp(sudd,selStr,ssLen)!=0);
      case UDSCR_GEncase   : return (strncasecmp(sudd,selStr,ssLen)>=0);
      case UDSCR_GTncase   : return (strncasecmp(sudd,selStr,ssLen)>=0);
      case UDSCR_Substr    : return (strstr(sudd,selStr)!=NULL);
      case UDSCR_NoSubstr  : return (strstr(sudd,selStr)==NULL);
      case UDSCR_Substr1   : return (strstr(selStr,sudd)!=NULL);
      case UDSCR_NoSubstr1 : return (strstr(selStr,sudd)==NULL);
      default              : return false;
    }
  }


  void SelManager::SelectUDD (
               int   selHnd,    // must be obtained from NewSelection()
               SELECTION_TYPE sType,   // selection type STYPE_XXXXX
               int   UDDhandle, // UDD handle
               cpstr selStr,    // selection string
               int   cmpRule,   // comparison rule
               SELECTION_KEY sKey     // selection key
             )  {
  PModel   mdl;
  PChain   chain;
  PResidue res;
  PAtom    atom;
  int      i,k,nsel,ssLen, n,c,r,a;
  bool     selAND;
  SELECTION_KEY sk;

    k  = selHnd-1;
    sk = sKey;

    if ((selType[k]==STYPE_UNDEFINED) || (sKey==SKEY_NEW))
            selType[k] = sType;
    else if (selType[k]!=sType)  return;


    if ((selHnd<=0) || (selHnd>nSelections))  return;

    switch (sType)  {
      case STYPE_ATOM    : if ((UDDhandle & UDRF_ATOM)==0)    return;
                        break;
      case STYPE_RESIDUE : if ((UDDhandle & UDRF_RESIDUE)==0) return;
                        break;
      case STYPE_CHAIN   : if ((UDDhandle & UDRF_CHAIN)==0)   return;
                        break;
      case STYPE_MODEL   : if ((UDDhandle & UDRF_MODEL)==0)   return;
                        break;
      default            : return;
    }


    // if something goes wrong, sk should be assigned SKEY_OR if
    // selKey is set to SKEY_NEW or SKEY_OR below
    switch (sKey)  {
      case SKEY_NEW : for (i=0;i<nSelItems[k];i++)
                        if (selection[k][i])
                            selection[k][i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : if (nSelItems[k]==0)  return;
                      nsel = 0;
                    break;
      case SKEY_XOR : nsel = nSelItems[k];  break;
      case SKEY_CLR : nsel = nSelItems[k];
                      if (nsel<=0)  return;
                    break;
      default       : return;
    }

    selAND = (sKey==SKEY_AND);
    ssLen = strlen ( selStr );

    for (n=0;n<nModels;n++)  {

      mdl = model[n];
      if (mdl)  {  // check for safety

        if (sType==STYPE_MODEL)  {

          if (selSUDD(mdl->getUDData(UDDhandle),selStr,
                                            cmpRule,ssLen))
            SelectObject ( mdl,k,sk,nsel );
          else if (selAND)
            mdl->RemoveMask ( mask[k] );

        } else  {

          for (c=0;c<mdl->nChains;c++)  {

            chain = mdl->chain[c];
            if (chain)  {   // again check for safety

              if (sType==STYPE_CHAIN)  {
                if (selSUDD(chain->getUDData(UDDhandle),selStr,
                                                  cmpRule,ssLen))
                  SelectObject ( chain,k,sk,nsel );
                else if (selAND)
                  chain->RemoveMask ( mask[k] );

              } else  {

                for (r=0;r<chain->nResidues;r++)  {

                  res = chain->residue[r];
                  if (res)  {

                    if (sType==STYPE_RESIDUE)  {
                      if (selSUDD(res->getUDData(UDDhandle),selStr,
                                                      cmpRule,ssLen))
                        SelectObject ( res,k,sk,nsel );
                      else if (selAND)
                        res->RemoveMask ( mask[k] );

                    } else  {

                      for (a=0;a<res->nAtoms;a++)  {
                        atom = res->atom[a];
                        if (atom)  {
                          if (!atom->Ter)  {
                            if (selSUDD(atom->getUDData(UDDhandle),selStr,
                                                             cmpRule,ssLen))
                              SelectObject ( atom,k,sk,nsel );
                            else if (selAND)
                              atom->RemoveMask ( mask[k] );
                          }
                        }

                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    MakeSelIndex ( selHnd,sType,nsel );

  }


  void SelManager::SelectSphere (
               int  selHnd, // must be obtained from NewSelection()
               SELECTION_TYPE sType, // selection type STYPE_XXXXX
               realtype  x, // x-coordinate of the sphere's center
               realtype  y, // y-coordinate of the sphere's center
               realtype  z, // z-coordinate of the sphere's center
               realtype  r, // radius of the sphere
               SELECTION_KEY sKey  // selection key
             )  {
  //  Selecting a sphere
  PPAtom    A;
  PAtom     atm;
  PResidue  res;
  PChain    chain;
  PModel    mdl;
  realtype  dx,dy,dz, r2;
  int       i,k, nat,nsel, im,ic,ir;
  bool      ASel, resSel,chainSel,modelSel,selAND;
  SELECTION_KEY sk;

    if ((selHnd<=0) || (selHnd>nSelections) || (r<=0.0))  return;

    k   = selHnd-1;
    sk  = sKey;
    A   = atom;
    nat = nAtoms;

    if ((selType[k]==STYPE_UNDEFINED) || (sKey==SKEY_NEW))
            selType[k] = sType;
    else if (selType[k]!=sType)  return;

    // if something goes wrong, sk should be assigned SKEY_OR if
    // selKey is set to SKEY_NEW or SKEY_OR below
    switch (sKey)  {
      case SKEY_NEW : for (i=0;i<nSelItems[k];i++)
                        if (selection[k][i])
                            selection[k][i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : nsel = 0;
                      nat  = nSelItems[k];
                      A    = (PPAtom)selection[k];
                    break;
      case SKEY_XOR : nsel = nSelItems[k];
                    break;
      case SKEY_CLR : nsel = nSelItems[k];
                      nat  = nSelItems[k];
                      A    = (PPAtom)selection[k];
                    break;
      default       : return;
    }

    selAND = (sKey==SKEY_AND);

    if ((nat<=0) || (!A))  return;

    r2 = r*r;

    if (sType==STYPE_ATOM)  {

      for (i=0;i<nat;i++)
        if (A[i])  {
          ASel = (sk!=SKEY_AND);
          if ((!A[i]->Ter) && (A[i]->WhatIsSet & ASET_Coordinates))  {
            dx = fabs(A[i]->x-x);
            if (dx<=r)  {
              dy = fabs(A[i]->y-y);
              if (dy<=r)  {
                dz = fabs(A[i]->z-z);
                if (dz<=r)  {
                  if (dx*dx+dy*dy+dz*dz<=r2)  {
                    ASel = true;
                    SelectAtom ( A[i],k,sk,nsel );
                  }
                }
              }
            }
          }
          if (!ASel)  A[i]->RemoveMask ( mask[k] );
        }

    } else  {

      for (im=0;im<nModels;im++)  {
        mdl = model[im];
        if (mdl)  {
          modelSel = false;
          for (ic=0;ic<mdl->nChains;ic++)  {
            chain = mdl->chain[ic];
            if (chain)  {
              chainSel = false;
              for (ir=0;ir<chain->nResidues;ir++)  {
                res = chain->residue[ir];
                if (res)  {
                  resSel = false;
                  for (i=0;i<res->nAtoms;i++)  {
                    atm = res->atom[i];
                    if (atm) {
                      ASel = false;
                      if ((!atm->Ter) &&
                          (atm->WhatIsSet & ASET_Coordinates))  {
                        dx = fabs(atm->x-x);
                        if (dx<=r)  {
                          dy = fabs(atm->y-y);
                          if (dy<=r)  {
                            dz = fabs(atm->z-z);
                            if (dz<=r)  {
                              if (dx*dx+dy*dy+dz*dz<=r2)  {
                                SelectObject ( sType,atm,k,sk,nsel );
                                ASel     = true;
                                resSel   = true;
                                chainSel = true;
                                modelSel = true;
                              }
                            }
                          }
                        }
                      }
                      if (ASel)  break;  // selType>=STYPE_RESIDUE
                    }
                  }
                  if ((!resSel) && selAND && (sType==STYPE_RESIDUE))
                    res->RemoveMask ( mask[k] );
                  if (chainSel && (sType>STYPE_RESIDUE))  break;
                }
              }
              if ((!chainSel) && selAND && (sType==STYPE_CHAIN))
                chain->RemoveMask ( mask[k] );
              if (modelSel && (sType>STYPE_CHAIN))  break;
            }
          }
          if ((!modelSel) && selAND && (sType==STYPE_MODEL))
            mdl->RemoveMask ( mask[k] );
        }
      }

    }

    MakeSelIndex ( selHnd,sType,nsel );

  }


  void SelManager::SelectCylinder (
               int  selHnd, // must be obtained from NewSelection()
               SELECTION_TYPE sType, // selection type STYPE_XXXXX
               realtype x1, // x-coordinate of the cylinder axis' 1st end
               realtype y1, // y-coordinate of the cylinder axis' 1st end
               realtype z1, // z-coordinate of the cylinder axis' 1st end
               realtype x2, // x-coordinate of the cylinder axis' 2nd end
               realtype y2, // y-coordinate of the cylinder axis' 2nd end
               realtype z2, // z-coordinate of the cylinder axis' 2nd end
               realtype  r, // radius of the cylinder
               SELECTION_KEY sKey  // selection key
             )  {
  //
  //  Selecting a cylinder
  //
  //  Method : given a line running through (x1,y1,z1) to (x2,y2,z2) on,
  //  a point (x,y,z) is then projected on it at distance
  //
  //              c1 = (c^2-a^2+b^2)/(2c),
  //
  //  from (x1,y1,z1), where
  //      'a' is the distance between (x,y,z) and (x2,y2,z2)
  //      'b' is the distance between (x,y,z) and (x1,y1,z1)
  //      'c' is the distance between (x1,y1,z1) and (x2,y2,z2).
  //  The distance between point (x,y,z) and line is determined as
  //
  //              h^2 = b^2 - c1^2
  //
  //  If c1>=0 and c1<=c and h^2<=r^2  then point (x,y,z) is inside
  //  a cylinder of radius 'r' with axis running from point
  //  (x1,y1,z1) to (x2,y2,z2).
  //
  PPAtom   A;
  PAtom    atm;
  PResidue res;
  PChain   chain;
  PModel   mdl;
  realtype dx,dy,dz, c,dc,c1,c2, a2,b2, r2;
  int      i,k, nat,nsel, im,ic,ir;
  bool     resSel,chainSel,modelSel,selAND;
  SELECTION_KEY sk;

    if ((selHnd<=0) || (selHnd>nSelections) || (r<=0.0))  return;

    dx = x1-x2;
    dy = y1-y2;
    dz = z1-z2;
    c2 = dx*dx + dy*dy + dz*dz;
    if (c2<=0.0)  return;
    c  = sqrt(c2);
    dc = 2.0*c;
    r2 = r*r;

    k   = selHnd-1;
    sk  = sKey;
    A   = atom;
    nat = nAtoms;

    if ((selType[k]==STYPE_UNDEFINED) || (sKey==SKEY_NEW))
            selType[k] = sType;
    else if (selType[k]!=sType)  return;

    // if something goes wrong, sk should be assigned SKEY_OR if
    // selKey is set to SKEY_NEW or SKEY_OR below
    switch (sKey)  {
      case SKEY_NEW : for (i=0;i<nSelItems[k];i++)
                        if (selection[k][i])
                            selection[k][i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : nsel = 0;
                      nat  = nSelItems[k];
                      A    = (PPAtom)selection[k];
                    break;
      case SKEY_XOR : nsel = nSelItems[k];
                    break;
      case SKEY_CLR : nsel = nSelItems[k];
                      nat  = nSelItems[k];
                      A    = (PPAtom)selection[k];
                    break;
      default       : return;
    }

    selAND = (sKey==SKEY_AND);

    if ((nat<=0) || (!A))  return;

    if (sType==STYPE_ATOM)  {

      for (i=0;i<nat;i++)
        if (A[i])  {
          if ((!A[i]->Ter) && (A[i]->WhatIsSet & ASET_Coordinates))  {
            dx = fabs(A[i]->x-x1);
            dy = fabs(A[i]->y-y1);
            dz = fabs(A[i]->z-z1);
            a2 = dx*dx + dy*dy + dz*dz;
            dx = fabs(A[i]->x-x2);
            dy = fabs(A[i]->y-y2);
            dz = fabs(A[i]->z-z2);
            b2 = dx*dx + dy*dy + dz*dz;
            c1 = (c2-a2+b2)/dc;
            if ((0.0<=c1) && (c1<=c) && (b2-c1*c1<=r2))
              SelectAtom ( A[i],k,sk,nsel );
            else if (sk==SKEY_AND)
              A[i]->RemoveMask ( mask[k] );
          }
        }

    } else  {

      for (im=0;im<nModels;im++)  {
        mdl = model[im];
        if (mdl)  {
          modelSel = false;
          for (ic=0;ic<mdl->nChains;ic++)  {
            chain = mdl->chain[ic];
            if (chain)  {
              chainSel = false;
              for (ir=0;ir<chain->nResidues;ir++)  {
                res = chain->residue[ir];
                if (res)  {
                  resSel = false;
                  for (i=0;i<res->nAtoms;i++)  {
                    atm = res->atom[i];
                    if (atm) {
                      if ((!atm->Ter) &&
                          (atm->WhatIsSet & ASET_Coordinates))  {
                        dx = fabs(atm->x-x1);
                        dy = fabs(atm->y-y1);
                        dz = fabs(atm->z-z1);
                        a2 = dx*dx + dy*dy + dz*dz;
                        dx = fabs(atm->x-x2);
                        dy = fabs(atm->y-y2);
                        dz = fabs(atm->z-z2);
                        b2 = dx*dx + dy*dy + dz*dz;
                        c1 = (c2-a2+b2)/dc;
                        if ((0.0<=c1) && (c1<=c) && (b2-c1*c1<=r2))  {
                          SelectObject ( sType,atm,k,sk,nsel );
                          resSel   = true;
                          chainSel = true;
                          modelSel = true;
                          break;  // selType>=STYPE_RESIDUE
                        }
                      }
                    }
                  }
                  if ((!resSel) && selAND && (sType==STYPE_RESIDUE))
                    res->RemoveMask ( mask[k] );
                  if (chainSel && (sType>STYPE_RESIDUE))  break;
                }
              }
              if ((!chainSel) && selAND && (sType==STYPE_CHAIN))
                chain->RemoveMask ( mask[k] );
              if (modelSel && (sType>STYPE_CHAIN))  break;
            }
          }
          if ((!modelSel) && selAND && (sType==STYPE_MODEL))
            mdl->RemoveMask ( mask[k] );
        }
      }

    }

    MakeSelIndex ( selHnd,sType,nsel );

  }


  void SelManager::SelectSlab (
               int  selHnd, // must be obtained from NewSelection()
               SELECTION_TYPE sType, // selection type STYPE_XXXXX
               realtype  a, // a-parameter of the plane  ax+by+cz=d
               realtype  b, // b-parameter of the plane  ax+by+cz=d
               realtype  c, // c-parameter of the plane  ax+by+cz=d
               realtype  d, // d-parameter of the plane  ax+by+cz=d
               realtype  r, // distance to the plane
               SELECTION_KEY sKey  // selection key
             )  {
  //
  //  Selecting all atoms on a given distance from a plane
  //
  //  Method : the distance between a point (x0,y0,z0) and a plane
  //  defined by equation
  //
  //              a*x + b*y + c*z = d
  //
  //  is found as
  //
  //              h = (d-a*x0-b*y0-c*z0)/sqrt(a^2+b^2+c^2)
  //
  //  If |h|<d then point (x0,y0,z0) belongs to the slab.
  //
  PPAtom   A;
  PAtom    atm;
  PResidue res;
  PChain   chain;
  PModel   mdl;
  realtype v,h;
  int      i,k, nat,nsel, im,ic,ir;
  bool     resSel,chainSel,modelSel,selAND;
  SELECTION_KEY sk;

    if ((selHnd<=0) || (selHnd>nSelections) || (r<=0.0))  return;

    v   = sqrt(a*a + b*b + c*c);
    if (v<=0.0)  return;

    k   = selHnd-1;
    sk  = sKey;
    A   = atom;
    nat = nAtoms;

    if ((selType[k]==STYPE_UNDEFINED) || (sKey==SKEY_NEW))
            selType[k] = sType;
    else if (selType[k]!=sType)  return;

    // if something goes wrong, sk should be assigned SKEY_OR if
    // selKey is set to SKEY_NEW or SKEY_OR below
    switch (sKey)  {
      case SKEY_NEW : for (i=0;i<nSelItems[k];i++)
                        if (selection[k][i])
                            selection[k][i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : nsel = 0;
                      nat  = nSelItems[k];
                      A    = (PPAtom)selection[k];
                    break;
      case SKEY_XOR : nsel = nSelItems[k];
                    break;
      case SKEY_CLR : nsel = nSelItems[k];
                      nat  = nSelItems[k];
                      A    = (PPAtom)selection[k];
                    break;
      default       : return;
    }

    selAND = (sKey==SKEY_AND);

    if ((nat<=0) || (!A))  return;

    if (sType==STYPE_ATOM)  {

      for (i=0;i<nat;i++)
        if (A[i])  {
          if ((!A[i]->Ter) && (A[i]->WhatIsSet & ASET_Coordinates))  {
            h = fabs(d-a*A[i]->x-b*A[i]->y-c*A[i]->z)/v;
            if (h<=r)
              SelectAtom ( A[i],k,sk,nsel );
            else if (sk==SKEY_AND)
              A[i]->RemoveMask ( mask[k] );
          }
        }

    } else  {

      for (im=0;im<nModels;im++)  {
        mdl = model[im];
        if (mdl)  {
          modelSel = false;
          for (ic=0;ic<mdl->nChains;ic++)  {
            chain = mdl->chain[ic];
            if (chain)  {
              chainSel = false;
              for (ir=0;ir<chain->nResidues;ir++)  {
                res = chain->residue[ir];
                if (res)  {
                  resSel = false;
                  for (i=0;i<res->nAtoms;i++)  {
                    atm = res->atom[i];
                    if (atm) {
                      if ((!atm->Ter) &&
                          (atm->WhatIsSet & ASET_Coordinates))  {
                        h = fabs(d-a*A[i]->x-b*A[i]->y-c*A[i]->z)/v;
                        if (h<=r)  {
                          SelectObject ( sType,atm,k,sk,nsel );
                          resSel   = true;
                          chainSel = true;
                          modelSel = true;
                          break;  // selType>=STYPE_RESIDUE
                        }
                      }
                    }
                  }
                  if ((!resSel) && selAND && (sType==STYPE_RESIDUE))
                    res->RemoveMask ( mask[k] );
                  if (chainSel && (sType>STYPE_RESIDUE))  break;
                }
              }
              if ((!chainSel) && selAND && (sType==STYPE_CHAIN))
                chain->RemoveMask ( mask[k] );
              if (modelSel && (sType>STYPE_CHAIN))  break;
            }
          }
          if ((!modelSel) && selAND && (sType==STYPE_MODEL))
            mdl->RemoveMask ( mask[k] );
        }
      }

    }

    MakeSelIndex ( selHnd,sType,nsel );

  }


  void SelManager::SelectNeighbours (
               int  selHnd, // must be obtained from NewSelection()
               SELECTION_TYPE sType, // selection type STYPE_XXXXX
               PPAtom  sA, // array of already selected atoms
               int    alen, // length of A
               realtype d1, // minimal distance to already selected atoms
               realtype d2, // maximal distance to already selected atoms
               SELECTION_KEY sKey  // selection key
                                         )  {
  // Selecting all atoms on a given distance from already selected
  PPAtom   A;
  PBrick   B;
  PAtom    atm;
  PResidue res;
  PChain   chain;
  PModel   mdl;
  realtype x,y,z, dx,dy,dz, dst, d12,d22;
  int      i,j,k, dn, nx,ny,nz, nat,nsel, im,ic,ir;
  int      ix1,ix2,ix, iy1,iy2,iy, iz1,iz2,iz;
  bool     ASel,resSel,chainSel,modelSel,selAND;
  SELECTION_KEY sk;

    if ((selHnd<=0) || (selHnd>nSelections) ||
        (d2<=0.0)   || (d2<d1))  return;

    k   = selHnd-1;
    sk  = sKey;
    A   = atom;
    nat = nAtoms;
    d12 = d1*d1;
    d22 = d2*d2;

    if ((selType[k]==STYPE_UNDEFINED) || (sKey==SKEY_NEW))
            selType[k] = sType;
    else if (selType[k]!=sType)  return;

    if ((alen<1) || (!sA))  {
      if ((sKey==SKEY_NEW) || (sKey==SKEY_AND))  {
        for (i=0;i<nSelItems[k];i++)
          if (selection[k][i])
            selection[k][i]->RemoveMask ( mask[k] );
        nSelItems[k] = 0;
      }
      return;
    }

    // if something goes wrong, sk should be assigned SKEY_OR if
    // selKey is set to SKEY_NEW or SKEY_OR below
    switch (sKey)  {
      case SKEY_NEW : for (i=0;i<nSelItems[k];i++)
                        if (selection[k][i])
                            selection[k][i]->RemoveMask ( mask[k] );
                      nSelItems[k] = 0;
                      nsel = 0;
                    break;
      case SKEY_OR  : if (nSelItems[k]==0)  sk = SKEY_NEW;
                      nsel = nSelItems[k];
                    break;
      case SKEY_AND : nsel = 0;
                      nat  = nSelItems[k];
                      A    = (PPAtom)selection[k];
                    break;
      case SKEY_XOR : nsel = nSelItems[k];
                    break;
      case SKEY_CLR : nsel = nSelItems[k];
                      nat  = nSelItems[k];
                      A    = (PPAtom)selection[k];
                    break;
      default       : return;
    }

    selAND = (sk==SKEY_AND);

    if ((nat<=0) || (!A))  return;

    MakeBricks ( sA,alen,d2*1.5 );
    dn = mround(d2/brick_size)+1;

    if (brick && (sType==STYPE_ATOM))  {

      for (i=0;i<nat;i++)
        if (A[i])  {
          if (!A[i]->Ter)  {
            ASel = false;
            GetBrickCoor ( A[i],nx,ny,nz );
            if (nx<0) nx++;
            ix1 = IMax ( 0,nx-dn );
            iy1 = IMax ( 0,ny-dn );
            iz1 = IMax ( 0,nz-dn );
            ix2 = IMin ( nbrick_x,nx+dn+1 );
            iy2 = IMin ( nbrick_y,ny+dn+1 );
            iz2 = IMin ( nbrick_z,nz+dn+1 );
            x   = A[i]->x;
            y   = A[i]->y;
            z   = A[i]->z;
            for (ix=ix1;(ix<ix2) && (!ASel);ix++)
              if (brick[ix])
                for (iy=iy1;(iy<iy2) && (!ASel);iy++)
                  if (brick[ix][iy])
                    for (iz=iz1;(iz<iz2) && (!ASel);iz++)  {
                      B = brick[ix][iy][iz];
                      if (B)
                        for (j=0;(j<B->nAtoms) && (!ASel);j++)
                          if (B->atom[j]!=A[i])  {
                            dx = fabs(x-B->atom[j]->x);
                            if (dx<=d2)  {
                              dy = fabs(y-B->atom[j]->y);
                              if (dy<=d2)  {
                                dz = fabs(z-B->atom[j]->z);
                                if (dz<=d2)  {
                                  dst = dx*dx+dy*dy+dz*dz;
                                  if ((dst>=d12) && (dst<=d22))  {
                                    ASel = true;
                                    SelectAtom ( A[i],k,sk,nsel );
                                  }
                                }
                              }
                            }
                          }
                    }
            if ((!ASel) && selAND)  A[i]->RemoveMask ( mask[k] );
          }
        }

    } else if (brick)  {

      for (im=0;im<nModels;im++)  {
        mdl = model[im];
        if (mdl)  {
          modelSel = false;
          for (ic=0;ic<mdl->nChains;ic++)  {
            chain = mdl->chain[ic];
            if (chain)  {
              chainSel = false;
              for (ir=0;ir<chain->nResidues;ir++)  {
                res = chain->residue[ir];
                if (res)  {
                  resSel = false;
                  for (i=0;(i<res->nAtoms) && (!resSel);i++)  {
                    atm = res->atom[i];
                    if (atm) {
                      if ((!atm->Ter) &&
                          (atm->WhatIsSet & ASET_Coordinates))  {
                        GetBrickCoor ( atm,nx,ny,nz );
                        if (nx<0) nx++;
                        ix1 = IMax ( 0,nx-dn );
                        iy1 = IMax ( 0,ny-dn );
                        iz1 = IMax ( 0,nz-dn );
                        ix2 = IMin ( nbrick_x,nx+dn+1 );
                        iy2 = IMin ( nbrick_y,ny+dn+1 );
                        iz2 = IMin ( nbrick_z,nz+dn+1 );
                        x   = atm->x;
                        y   = atm->y;
                        z   = atm->z;
                        for (ix=ix1;(ix<ix2) && (!resSel);ix++)
                          if (brick[ix])
                            for (iy=iy1;(iy<iy2) && (!resSel);iy++)
                              if (brick[ix][iy])
                                for (iz=iz1;(iz<iz2) && (!resSel);iz++) {
                                  B = brick[ix][iy][iz];
                                  if (B)
                                    for (j=0;(j<B->nAtoms) &&
                                             (!resSel);j++)
                                      if (B->atom[j]!=atm)  {
                                        dx = fabs(x-B->atom[j]->x);
                                        if (dx<=d2)  {
                                          dy = fabs(y-B->atom[j]->y);
                                          if (dy<=d2)  {
                                            dz = fabs(z-B->atom[j]->z);
                                            if (dz<=d2)  {
                                              dst = dx*dx+dy*dy+dz*dz;
                                              if ((dst>=d12) &&
                                                  (dst<=d22))  {
                                                SelectObject ( sType,
                                                        atm,k,sk,nsel );
                                                resSel   = true;
                                                chainSel = true;
                                                modelSel = true;
                                              }
                                            }
                                          }
                                        }
                                      }
                                }

                      }
                    }
                  }
                  if ((!resSel) && selAND && (sType==STYPE_RESIDUE))
                    res->RemoveMask ( mask[k] );
                  if (chainSel && (sType>STYPE_RESIDUE))  break;
                }
              }
              if ((!chainSel) && selAND && (sType==STYPE_CHAIN))
                chain->RemoveMask ( mask[k] );
              if (modelSel && (sType>STYPE_CHAIN))  break;
            }
          }
          if ((!modelSel) && selAND && (sType==STYPE_MODEL))
            mdl->RemoveMask ( mask[k] );
        }
      }

    }

    MakeSelIndex ( selHnd,sType,nsel );

  }



  int TakeChainID ( pstr & p, pstr chainID )  {
  int RC,k;
    chainID[0] = char(0);
    if (*p)  {
      RC = 0;
      if (*p==':')  {
        // starts with colon <=> empty chain ID
        chainID[0] = char(0);
        p++;  // advance to residue number
      } else if (p[1]==':')  {
        // second symbol is colon <=> regular chain ID
        chainID[0] = *p;
        chainID[1] = char(0);
        p++;
        p++;  // advance to residue number
      } else if (*p=='\'')  {
        // starts with a strip <=> assume empty chain ID
        chainID[0] = char(0);
        p++;
        if (*p=='\'')  {
          // closing strip must be followed by colon
          p++;
          if (*p!=':')  RC = -1;
        } else  {
          // no concluding strip <=> could be a strip chain ID,
          // although this must be captured by 2nd "if"
          chainID[0] = '\'';
          chainID[1] = char(0);
          // assume that residue number is following the strip
        }
      } else if ((int(*p)>=int('0')) && (int(*p)<=int('9')))  {
        // a digit without following semicolon looks very much
        // like residue number with unspecified empty chain ID
        chainID[0] = char(0);
        // assume that p already points on residue number
      } else  {
        // assume a long chain ID
        k = 0;
        while (*p && (*p!=':') && (k<(int)sizeof(ChainID)-1)) {
          chainID[k++] = *p;
          p++;
        }
        if (*p==':')  {
          chainID[k] = char(0);
        } else  {
          // a mistake
          chainID[0] = char(0);
          RC = -1;
        }
      }
      while (*p==' ')  p++;
    } else
      RC = 1;
    return RC;
  }

  int TakeResID( pstr & p, int & seqNum, pstr inscode )  {
  char N[100];
  int  i,RC;
  pstr endptr;
    RC = 1;
    inscode[0] = '*';
    inscode[1] = char(0);
    seqNum     = ANY_RES;
    if (((*p) &&
         (int(*p)>=int('0')) && (int(*p)<=int('9'))) || (*p=='-'))  {
      N[0] = *p;
      p++;
      i = 1;
      while ((*p) && (int(*p)>=int('0')) && (int(*p)<=int('9')))  {
        N[i++] = *p;
        p++;
      }
      N[i] = char(0);
      seqNum = mround(strtod(N,&endptr));
      if ((seqNum==0) && (endptr==N))
        RC = -1;
      else  {
        RC = 0;
        if ((*p) && (*p!='-') && (*p!=',') && (*p!=' '))  {
          inscode[0] = *p;
          inscode[1] = char(0);
          p++;
        } else
          inscode[0] = char(0);
        if ((*p=='-') || (*p==','))  p++;
      }
      while (*p==' ')  p++;
    }
    return RC;
  }


  int  SelManager::SelectDomain ( int selHnd , cpstr domainRange,
                                  SELECTION_TYPE sType,
                                  SELECTION_KEY  sKey,
                                  int modelNo )  {
  // domainRange is of the following format:
  //    "*", "(all)"            - take all file
  //    "-"                     - take chain without chain ID
  //    "a:Ni-Mj,b:Kp-Lq,..."   - take chain a residue number N
  //                              insertion code i to residue number M
  //                              insertion code j plus chain b
  //                              residue number K insertion code p to
  //                              residue number L insertion code q and
  //                              so on.
  //    "a:,b:..."              - take whole chains a and b and so on
  //    "a:,b:Kp-Lq,..."        - any combination of the above.
  ChainID       chainID;
  InsCode       insCode1,insCode2;
  pstr          S,p;
  int           seqNum1,seqNum2,rc;
  SELECTION_KEY selKey1;

    if ((selHnd<=0) || (selHnd>nSelections))  return 1;

    // leave only required residues

    rc = 1;
    if (!domainRange)  rc = 0;
    else if ((!domainRange[0]) || (domainRange[0]=='*'))  rc = 0;
    else if (!strcasecmp(domainRange,"(all)"))  rc = 0;
    if (!rc)  {
      // select all
      Select ( selHnd,sType,modelNo,"*",ANY_RES,"*",ANY_RES,"*",
                                      "*","*","*","*",sKey );
      return 0;
    }
    if (!strcasecmp(domainRange,"-"))  {
      // select chain without chain ID
      Select ( selHnd,sType,modelNo,"",ANY_RES,"*",ANY_RES,"*",
                                      "*","*","*","*",sKey );
      return 0;
    }

    S = new char[strlen(domainRange)+10];
    strcpy    ( S,domainRange );
    DelSpaces ( S );
  //  UpperCase ( S );

    p  = S;
    rc = 0;

    selKey1 = sKey;

    while ((*p) && (!rc))  {

      if (TakeChainID(p,chainID)<0)             rc = -1;
      else if (TakeResID(p,seqNum1,insCode1)<0) rc = -2;
      else if (TakeResID(p,seqNum2,insCode2)<0) rc = -3;
      else  {
        Select ( selHnd,sType,modelNo,chainID,
                 seqNum1,insCode1,seqNum2,insCode2,
                 "*","*","*","*",selKey1 );
        if (*p==',')  p++;
        if (selKey1==SKEY_NEW)  selKey1 = SKEY_OR;
      }

    }

    delete[] S;

    return rc;

  }


  int SelManager::GetSelLength ( int selHnd )  {
    if ((selHnd>0) && (selHnd<=nSelections))
          return nSelItems[selHnd-1];
    else  return 0;
  }


  void SelManager::GetSelIndex ( int       selHnd,
                                      PPAtom & SelAtom,
                                      int &     nSelAtoms )  {
    if ((selHnd>0) && (selHnd<=nSelections))  {
      if (selType[selHnd-1]!=STYPE_ATOM)  {
        SelAtom   = NULL;
        nSelAtoms = 0;
      } else  {
        SelAtom   = (PPAtom)selection[selHnd-1];
        nSelAtoms = nSelItems[selHnd-1];
      }
    } else  {
      SelAtom   = NULL;
      nSelAtoms = 0;
    }
  }

  void SelManager::GetSelIndex ( int          selHnd,
                                      PPResidue & SelResidue,
                                      int &        nSelResidues )  {
    if ((selHnd>0) && (selHnd<=nSelections))  {
      if (selType[selHnd-1]!=STYPE_RESIDUE)  {
        SelResidue   = NULL;
        nSelResidues = 0;
      } else  {
        SelResidue   = (PPResidue)selection[selHnd-1];
        nSelResidues = nSelItems[selHnd-1];
      }
    } else  {
      SelResidue   = NULL;
      nSelResidues = 0;
    }
  }

  void SelManager::GetSelIndex ( int        selHnd,
                                      PPChain & SelChain,
                                      int &      nSelChains )  {
    if ((selHnd>0) && (selHnd<=nSelections))  {
      if (selType[selHnd-1]!=STYPE_CHAIN)  {
        SelChain   = NULL;
        nSelChains = 0;
      } else  {
        SelChain   = (PPChain)selection[selHnd-1];
        nSelChains = nSelItems[selHnd-1];
      }
    } else  {
      SelChain   = NULL;
      nSelChains = 0;
    }
  }

  void SelManager::GetSelIndex ( int        selHnd,
                                      PPModel & SelModel,
                                      int &      nSelModels )  {
    if ((selHnd>0) && (selHnd<=nSelections))  {
      if (selType[selHnd-1]!=STYPE_MODEL)  {
        SelModel   = NULL;
        nSelModels = 0;
      } else  {
        SelModel   = (PPModel)selection[selHnd-1];
        nSelModels = nSelItems[selHnd-1];
      }
    } else  {
      SelModel   = NULL;
      nSelModels = 0;
    }
  }


  void SelManager::GetAtomStatistics ( int selHnd, RAtomStat AS )  {
  int  i,k;
    AS.Init();
    if ((selHnd>0) && (selHnd<=nSelections))  {
      k = selHnd-1;
      switch (selType[k])  {
        case STYPE_MODEL   : if (selection[k])
                               for (i=0;i<nSelItems[k];i++)
                                 ((PModel)selection[k][i])->
                                   CalAtomStatistics ( AS );
                         break;
        case STYPE_CHAIN   : if (selection[k])
                               for (i=0;i<nSelItems[k];i++)
                                 ((PChain)selection[k][i])->
                                   CalAtomStatistics ( AS );
                         break;
        case STYPE_RESIDUE : if (selection[k])
                               for (i=0;i<nSelItems[k];i++)
                                 ((PResidue)selection[k][i])->
                                   CalAtomStatistics ( AS );
                         break;
        case STYPE_ATOM    : if (selection[k])
                               for (i=0;i<nSelItems[k];i++)
                                 ((PAtom)selection[k][i])->
                                   CalAtomStatistics ( AS );
                         break;
        default            : break;
      }
    }
    AS.Finish();
  }


  void SelManager::SelectAtom ( PAtom atom, int maskNo,
                                SELECTION_KEY sKey, int & nsel )  {
  bool ASel;
    ASel = atom->CheckMask ( mask[maskNo] );
    switch (sKey)  {
      default       :
      case SKEY_NEW :
      case SKEY_OR  : if (!ASel)  {
                        atom->SetMask ( mask[maskNo] );
                        nsel++;
                      }
                    break;
      case SKEY_AND : if (ASel)  nsel++;
                    break;
      case SKEY_XOR : if (ASel)  {
                        atom->RemoveMask ( mask[maskNo] );
                        nsel--;
                      } else  {
                        atom->SetMask ( mask[maskNo] );
                        nsel++;
                      }
                    break;
      case SKEY_CLR : if (ASel)  {
                        atom->RemoveMask ( mask[maskNo] );
                        nsel--;
                      }
    }
  }


  void SelManager::SelectObject ( SELECTION_TYPE sType,
                                  PAtom          atm,
                                  int            maskNo,
                                  SELECTION_KEY  sKey,
                                  int          & nsel )  {
  PMask  object;
    switch (sType)  {
      default :
      case STYPE_UNDEFINED : return;
      case STYPE_ATOM      : object = atm;                break;
      case STYPE_RESIDUE   : object = atm->GetResidue();  break;
      case STYPE_CHAIN     : object = atm->GetChain  ();  break;
      case STYPE_MODEL     : object = atm->GetModel  ();  break;
    }
    if (!object)  return;
    SelectObject ( object,maskNo,sKey,nsel );
  }


  void SelManager::SelectObject ( PMask object, int maskNo,
                                  SELECTION_KEY sKey, int & nsel )  {
  bool ASel;
    ASel = object->CheckMask ( mask[maskNo] );
    switch (sKey)  {
      default        :
      case SKEY_NEW  :
      case SKEY_OR   : if (!ASel)  {
                         object->SetMask ( mask[maskNo] );
                         nsel++;
                       }
                    break;
      case SKEY_AND  : if (ASel)  nsel++;
                    break;
      case SKEY_XOR  : if (ASel)  {
                         object->RemoveMask ( mask[maskNo] );
                         nsel--;
                       } else  {
                         object->SetMask ( mask[maskNo] );
                         nsel++;
                       }
                    break;
      case SKEY_CLR  : if (ASel)  {
                         object->RemoveMask ( mask[maskNo] );
                         nsel--;
                       }
                    break;
      case SKEY_XAND : if (ASel)  {
                         object->RemoveMask ( mask[maskNo] );
                         nsel++;
                       }
    }
  }


  void  SelManager::DeleteSelObjects ( int selHnd )  {
  PPModel   model;
  PPChain   chain;
  PPResidue res;
  PPAtom    atom;
  int        i,k,nSel;

    if ((selHnd>0) && (selHnd<=nSelections))  {

      k    = selHnd-1;
      nSel = nSelItems[k];
      switch (selType[k])  {

        case STYPE_MODEL   : model = (PPModel)selection[k];
                             for (i=0;i<nSel;i++)
                               delete model[i];
                          break;

        case STYPE_CHAIN   : chain = (PPChain)selection[k];
                             for (i=0;i<nSel;i++)
                               delete chain[i];
                          break;

        case STYPE_RESIDUE : res   = (PPResidue)selection[k];
                             for (i=0;i<nSel;i++)
                               delete res[i];
                          break;

        case STYPE_ATOM    : atom  = (PPAtom)selection[k];
                             for (i=0;i<nSel;i++)
                               delete atom[i];
                          break;

        default : ;

      }

      if (selection[k])  delete[] selection[k];
      selection[k] = NULL;
      nSelItems[k] = 0;

    }

  }

  // ------------------------------------------------------------------------

  void SelManager::MakeSelIndex ( int selHnd,
                                  SELECTION_TYPE sType, int nsel )  {
  // if nsel is less than 0, the number of selected atoms will
  // be calculated.
  PModel   mdl;
  PChain   chain;
  PResidue res;
  int      k,i,j,n,ns,k1,k2, nns;

    if ((selHnd>0) && (selHnd<=nSelections))  {
      k1 = selHnd-1;
      k2 = k1+1;
    } else  {
      k1 = 0;
      k2 = nSelections;
    }

    for (k=k1;k<k2;k++)  {
      if (nsel<0)  {
        ns = 0;
        switch (sType)  {
          case STYPE_ATOM    : for (i=0;i<nAtoms;i++)
                                 if (atom[i])
                                   if (atom[i]->CheckMask(mask[k]))  ns++;
                            break;
          case STYPE_RESIDUE : for (n=0;n<nModels;n++)  {
                                 mdl = model[n];
                                 if (mdl)
                                   for (i=0;i<mdl->nChains;i++) {
                                     chain = mdl->chain[i];
                                     if (chain)
                                       for (j=0;j<chain->nResidues;j++) {
                                         res = chain->residue[j];
                                         if (res)
                                           if (res->CheckMask(mask[k]))  ns++;
                                       }
                                   }
                               }
                            break;
          case STYPE_CHAIN   : for (i=0;i<nModels;i++)  {
                                 mdl = model[i];
                                 if (mdl)
                                   for (j=0;j<mdl->nChains;j++) {
                                     chain = mdl->chain[j];
                                     if (chain)
                                       if (chain->CheckMask(mask[k]))  ns++;
                                   }
                               }
                            break;
          case STYPE_MODEL   : for (i=0;i<nModels;i++)
                                 if (model[i])
                                   if (model[i]->CheckMask(mask[k]))  ns++;
                            break;
          default : ;
        }
      } else
        ns = nsel;
      if (selection[k])  delete[] selection[k];
      if (ns>0)  {
        selection[k] = new PMask[ns];
        nns = 0;
        switch (sType)  {
          case STYPE_ATOM    : for (i=0;i<nAtoms;i++)
                                 if (atom[i])  {
                                   if (atom[i]->CheckMask(mask[k]))  {
                                     selection[k][nns++] = atom[i];
                                     if (nns>=ns)  nns = ns-1;
                                   }
                                 }
                            break;
          case STYPE_RESIDUE : for (n=0;n<nModels;n++)  {
                                 mdl = model[n];
                                 if (mdl)
                                   for (i=0;i<mdl->nChains;i++) {
                                     chain = mdl->chain[i];
                                     if (chain)
                                       for (j=0;j<chain->nResidues;j++)  {
                                         res = chain->residue[j];
                                         if (res)
                                           if (res->CheckMask(mask[k]))  {
                                             selection[k][nns++] = res;
                                             if (nns>=ns)  nns = ns-1;
                                           }
                                       }
                                   }
                               }
                            break;
          case STYPE_CHAIN   : for (i=0;i<nModels;i++)  {
                                 mdl = model[i];
                                 if (mdl)
                                   for (j=0;j<mdl->nChains;j++) {
                                     chain = mdl->chain[j];
                                     if (chain)
                                       if (chain->CheckMask(mask[k]))  {
                                         selection[k][nns++] = chain;
                                         if (nns>=ns)  nns = ns-1;
                                       }
                                   }
                               }
                            break;
          case STYPE_MODEL   : for (i=0;i<nModels;i++)
                                 if (model[i])
                                   if (model[i]->CheckMask(mask[k]))  {
                                     selection[k][nns++] = model[i];
                                     if (nns>=ns)  nns = ns-1;
                                   }
                            break;
          default : ;
        }

      } else
        selection[k] = NULL;

      nSelItems[k] = ns;
    }

  }


  //  -------------------  Stream functions  ----------------------


  void  SelManager::write ( io::RFile f )  {
  int  i,sType;
  byte Version=1;

    f.WriteByte ( &Version );

    CoorManager::write ( f );
    
    if (!isCompactBinary())  {
      f.WriteInt ( &nSelections );
      for (i=0;i<nSelections;i++)  {
        StreamWrite ( f,mask[i]       );
        f.WriteInt  ( &(nSelItems[i]) );
        sType = selType[i];
        f.WriteInt  ( &(sType) );
      }
    }

  }

  void  SelManager::read ( io::RFile f )  {
  int  i,sType;
  byte Version;

    f.ReadByte ( &Version );

    DeleteAllSelections();

    CoorManager::read ( f );

    if (!isCompactBinary())  {
      f.ReadInt ( &nSelections );
      if (nSelections>0)  {
        mask      = new PMask [nSelections];
        selection = new PPMask[nSelections];
        nSelItems = new int    [nSelections];
        selType   = new SELECTION_TYPE[nSelections];
        for (i=0;i<nSelections;i++)  {
          mask[i] = NULL;
          StreamRead ( f,mask[i]       );
          f.ReadInt  ( &(nSelItems[i]) );
          f.ReadInt  ( &(sType)        );
          selType  [i] = (SELECTION_TYPE)sType;
          selection[i] = NULL;
          if (mask[i])
               MakeSelIndex ( i+1,selType[i],-1 );
          else nSelItems[i] = 0;
        }
      }
    }

  }


  MakeStreamFunctions(SelManager)

} // namespace mmdb

