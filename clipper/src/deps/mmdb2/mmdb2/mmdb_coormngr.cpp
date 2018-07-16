//  $Id: mmdb_coormngr.h $
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
//    23.10.15   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  mmdb_coormngr <implementation>
//       ~~~~~~~~~
//       Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::Brick       ( space brick                  )
//       ~~~~~~~~~  mmdb::CoorManager ( MMDB atom coordinate manager )
//
//  (C) E. Krissinel 2000-2015
//
//  =================================================================
//

#include <math.h>
#include <string.h>

#include "mmdb_coormngr.h"
#include "mmdb_math_linalg.h"
#include "mmdb_tables.h"

namespace mmdb  {

  // ==========================  Brick  =============================

  Brick::Brick()  {
    InitBrick();
  }

  Brick::~Brick()  {
    Clear();
  }

  void  Brick::InitBrick()  {
    atom        = NULL;
    id          = NULL;
    nAtoms      = 0;
    nAllocAtoms = 0;
  }

  void  Brick::Clear()  {
    if (atom)  delete[] atom;
    FreeVectorMemory ( id,0 );
    atom        = NULL;
    nAtoms      = 0;
    nAllocAtoms = 0;
  }

  void  Brick::AddAtom ( PAtom A, int atomid )  {
  int     i;
  PPAtom  atom1;
  ivector id1;
    if (nAtoms>=nAllocAtoms)  {
      nAllocAtoms = nAtoms+10;
      atom1       = new PAtom[nAllocAtoms];
      GetVectorMemory ( id1,nAllocAtoms,0 );
      for (i=0;i<nAtoms;i++)  {
        atom1[i] = atom[i];
        id1  [i] = id  [i];
      }
      for (i=nAtoms;i<nAllocAtoms;i++)  {
        atom1[i] = NULL;
        id1  [i] = -1;
      }
      if (atom)  delete[] atom;
      FreeVectorMemory ( id,0 );
      atom = atom1;
      id   = id1;
    }
    atom[nAtoms] = A;
    id  [nAtoms] = atomid;
    nAtoms++;
  }


  // ===========================  mbrick  =============================

  MBrick::MBrick ( int nStructures )  {
    InitMBrick ( nStructures );
  }

  MBrick::~MBrick()  {
    Clear();
  }

  void  MBrick::InitMBrick ( int nStructures )  {
  int i;
    nStruct = nStructures;
    atom   = new PPAtom[nStruct];
    id     = new ivector[nStruct];
    GetVectorMemory ( nAtoms,nStruct,0 );
    GetVectorMemory ( nAlloAtoms,nStruct,0 );
    for (i=0;i<nStruct;i++)  {
      atom       [i] = NULL;
      id         [i] = NULL;
      nAtoms     [i] = 0;
      nAlloAtoms[i] = 0;
    }
  }

  void  MBrick::Clear()  {
  int i;
    if (atom)  {
      for (i=0;i<nStruct;i++)
        if (atom[i])  delete[] atom[i];
      delete[] atom;
      atom = NULL;
    }
    FreeMatrixMemory ( id,nStruct,0,0 );
    FreeVectorMemory ( nAtoms,0 );
    FreeVectorMemory ( nAlloAtoms,0 );
    nStruct = 0;
  }

  void  MBrick::AddAtom ( PAtom A, int structNo, int atomid )  {
  int     i,natoms,nalloc;
  PPAtom  atom0,atom1;
  ivector id0,id1;
    natoms = nAtoms     [structNo];
    nalloc = nAlloAtoms[structNo];
    atom0  = atom       [structNo];
    id0    = id         [structNo];
    if (natoms>=nalloc)  {
      nalloc = natoms+10;
      atom1 = new PAtom[nalloc];
      GetVectorMemory ( id1,nalloc,0 );
      for (i=0;i<natoms;i++)  {
        atom1[i] = atom0[i];
        id1  [i] = id0  [i];
      }
      for (i=natoms;i<nalloc;i++)  {
        atom1[i] = NULL;
        id1  [i] = -1;
      }
      if (atom0)  delete[] atom0;
      FreeVectorMemory ( id0,0 );
      atom[structNo] = atom1;
      id  [structNo] = id1;
      nAlloAtoms[structNo] = nalloc;
      atom0 = atom1;
      id0   = id1;
    }
    atom0 [natoms]   = A;
    id0   [natoms]   = atomid;
    nAtoms[structNo] = natoms+1;
  }



  //  ====================  GenSym  ========================

  GenSym::GenSym() : SymOps()  {
    InitGenSym();
  }

  GenSym::GenSym ( io::RPStream Object ) : SymOps(Object)  {
    InitGenSym();
  }

  GenSym::~GenSym()  {}  // virtual FreeMmeory is called by ~SymOps()

  void GenSym::InitGenSym()  {
    chID1    = NULL;
    chID2    = NULL;
    nChains  = NULL;
    nOpAlloc = 0;
  }

  void GenSym::FreeMemory()  {
  int i;
    for (i=0;i<nOpAlloc;i++)  {
      if (chID1[i]) delete[] chID1[i];
      if (chID2[i]) delete[] chID2[i];
    }
    if (chID1) delete[] chID1;
    if (chID2) delete[] chID2;
    FreeVectorMemory ( nChains,0 );
    nOpAlloc = 0;
    SymOps::FreeMemory();
  }

  int GenSym::AddSymOp ( cpstr XYZOperation )  {
  int        RC,i;
  PChainID * ch1ID;
  PChainID * ch2ID;
  ivector    nChains1;

    RC = SymOps::AddSymOp ( XYZOperation );
    if (Nops>nOpAlloc)  {
      ch1ID = new PChainID[Nops];
      ch2ID = new PChainID[Nops];
      GetVectorMemory ( nChains1,Nops,0 );
      for (i=0;i<nOpAlloc;i++)  {
        ch1ID[i]    = chID1[i];
        ch2ID[i]    = chID2[i];
        nChains1[i] = nChains[i];
      }
      for (i=nOpAlloc;i<Nops;i++)  {
        ch1ID[i]    = NULL;
        ch2ID[i]    = NULL;
        nChains1[i] = 0;
      }
      if (chID1)  delete[] chID1;
      if (chID2)  delete[] chID2;
      FreeVectorMemory ( nChains,0 );
      chID1    = ch1ID;
      chID2    = ch2ID;
      nChains  = nChains1;
      nOpAlloc = Nops;
    }
    return RC;
  }

  int  GenSym::AddRenChain ( int Nop, const ChainID ch1,
                                       const ChainID ch2 )  {
  int      i;
  PChainID c1,c2;
    if ((0<=Nop) && (Nop<Nops))  {
      c1 = new ChainID[nChains[Nop]+1];
      c2 = new ChainID[nChains[Nop]+1];
      for (i=0;i<nChains[Nop];i++)  {
        strcpy ( c1[i],chID1[Nop][i] );
        strcpy ( c2[i],chID2[Nop][i] );
      }
      strcpy ( c1[nChains[Nop]],ch1 );
      strcpy ( c2[nChains[Nop]],ch2 );
      if (chID1[Nop])  delete[] chID1[Nop];
      if (chID2[Nop])  delete[] chID2[Nop];
      chID1[Nop] = c1;
      chID2[Nop] = c2;
      nChains[Nop]++;
      return SYMOP_Ok;
    } else
      return SYMOP_NoSymOps;
  }

  void GenSym::Copy ( PSymOps GenSym )  {
  int i,j;
    SymOps::Copy ( GenSym );
    if (Nops>0)  {
      nOpAlloc = Nops;
      chID1 = new PChainID[Nops];
      chID2 = new PChainID[Nops];
      GetVectorMemory ( nChains,Nops,0 );
      for (i=0;i<Nops;i++)  {
        nChains[i] = PGenSym(GenSym)->nChains[i];
        if (nChains[i]<=0)  {
          chID1[i] = NULL;
          chID2[i] = NULL;
        } else  {
          chID1[i] = new ChainID[nChains[i]];
          chID2[i] = new ChainID[nChains[i]];
          for (j=0;j<nChains[i];j++)  {
            strcpy ( chID1[i][j],PGenSym(GenSym)->chID1[i][j] );
            strcpy ( chID2[i][j],PGenSym(GenSym)->chID2[i][j] );
          }
        }
      }
    }
  }

  void  GenSym::write ( io::RFile f )  {
  int  i,j;
  byte Version=1;
    f.WriteByte ( &Version  );
    SymOps::write ( f );
    f.WriteInt ( &nOpAlloc );
    for (i=0;i<nOpAlloc;i++)  {
      f.WriteInt ( &(nChains[i]) );
      for (j=0;j<nChains[i];j++)  {
        f.WriteTerLine ( chID1[i][j],false );
        f.WriteTerLine ( chID2[i][j],false );
      }
    }
  }

  void  GenSym::read ( io::RFile f )  {
  int  i,j;
  byte Version;
    f.ReadByte ( &Version  );
    SymOps::read ( f );
    f.ReadInt ( &nOpAlloc );
    if (nOpAlloc>0)  {
      chID1 = new PChainID[nOpAlloc];
      chID2 = new PChainID[nOpAlloc];
      GetVectorMemory ( nChains,nOpAlloc,0 );
      for (i=0;i<nOpAlloc;i++)  {
        f.ReadInt ( &(nChains[i]) );
        if (nChains[i]>0)  {
          chID1[i] = new ChainID[nChains[i]];
          chID2[i] = new ChainID[nChains[i]];
          for (j=0;j<nChains[i];j++)  {
            f.ReadTerLine ( chID1[i][j],false );
            f.ReadTerLine ( chID2[i][j],false );
          }
        } else  {
          chID1[i] = NULL;
          chID2[i] = NULL;
        }
      }
    }
  }


  MakeStreamFunctions(GenSym)



  // =======================  ContactIndex  ==========================

  void Contact::Copy ( RContact c )  {
    id1   = c.id1;
    id2   = c.id2;
    group = c.group;
    dist  = c.dist;
  }

  void Contact::Swap ( RContact c )  {
  int      ib;
  long     lb;
  realtype rb;
    ib = id1;     id1   = c.id1;     c.id1   = ib;
    ib = id2;     id2   = c.id2;     c.id2   = ib;
    lb = group;   group = c.group;   c.group = lb;
    rb = dist;    dist  = c.dist;    c.dist  = rb;
  }

  DefineClass(ContactIndex)

  class ContactIndex  {

    friend class SelManager;

    public :

      ContactIndex ( PContact contact,
                      int       maxlen,
                      int       ncontacts,
                      int       max_alloc );
      ~ContactIndex();

      void AddContact ( int id1, int id2,   realtype dist, int group  );
      void GetIndex   ( RPContact contact, int & ncontacts );

    protected :

      PContact contact_index; // contact index
      int      max_index;     // if <=0 then dynamical index
                              // otherwise fixed by max_index
      int      n_contacts;    // number of contacts
      int      alloc_index;   // physical length of contact_index
                              // when dynamical
      int      alloc_max;     // physical limit on allocation

  };


  ContactIndex::ContactIndex ( PContact contact,
                               int       maxlen,
                               int       ncontacts,
                               int       max_alloc )  {
    contact_index = contact;
    max_index     = maxlen;
    if (!contact_index)  n_contacts = 0;
                   else  n_contacts = IMax(0,ncontacts);
    alloc_index = n_contacts;
    alloc_max   = n_contacts + max_alloc;
  }

  ContactIndex::~ContactIndex() {
    if (contact_index)  delete[] contact_index;
    contact_index = NULL;
    n_contacts    = 0;
    alloc_index   = 0;
  }

  void ContactIndex::AddContact ( int id1, int id2, realtype dist,
                                  int group )  {
  PContact cont1;
  int      i;

    if ((alloc_max<=0) || (n_contacts<alloc_max))  {
      if (max_index>0)  {
        if (n_contacts<max_index)  {
          contact_index[n_contacts].id1   = id1;
          contact_index[n_contacts].id2   = id2;
          contact_index[n_contacts].dist  = dist;
          contact_index[n_contacts].group = group;
        }
      } else  {
        if (n_contacts>=alloc_index)  {
          alloc_index = n_contacts+IMax(alloc_index/4+10,10);
          if ((alloc_max>0) && (alloc_index>alloc_max))
            alloc_index = alloc_max;
          cont1 = new Contact[alloc_index];
          for (i=0;i<n_contacts;i++)
            cont1[i].Copy ( contact_index[i] );
          if (contact_index)  delete[] contact_index;
          contact_index = cont1;
        }
        contact_index[n_contacts].id1   = id1;
        contact_index[n_contacts].id2   = id2;
        contact_index[n_contacts].dist  = dist;
        contact_index[n_contacts].group = group;
      }
      n_contacts++;
    }
  }

  void  ContactIndex::GetIndex ( RPContact contact, int & ncontacts )  {
    contact       = contact_index;
    ncontacts     = n_contacts;
    contact_index = NULL;
    n_contacts    = 0;
    alloc_index   = 0;
  }


  // ========================  MContact  =============================

  MContact::MContact ( int nStructures )  {
  int i;
    nStruct = nStructures;
    if (nStruct>0)  {
      atom = new PPAtom[nStruct];
      id   = new ivector[nStruct];
      GetVectorMemory ( nAtoms,nStruct,0 );
      GetVectorMemory ( nAlloc,nStruct,0 );
      for (i=0;i<nStruct;i++)  {
        atom  [i] = NULL;
        id    [i] = NULL;
        nAtoms[i] = 0;
        nAlloc[i] = 0;
      }
    } else  {
      atom   = NULL;
      nAtoms = NULL;
      nAlloc = NULL;
    }
  }

  MContact::~MContact()  {
  int i;
    if (atom)  {
      for (i=0;i<nStruct;i++)
        if (atom[i])  delete[] atom[i];
      delete[] atom;
      atom = NULL;
    }
    FreeMatrixMemory ( id,nStruct,0,0 );
    FreeVectorMemory ( nAtoms,0 );
    FreeVectorMemory ( nAlloc,0 );
    nStruct = 0;
  }

  void MContact::AddContact ( PAtom A, int structNo, int atomid )  {
  PPAtom  A1,A2;
  ivector id1,id2;
  int     nat,nal,i;
    A1  = atom  [structNo];
    id1 = id    [structNo];
    nat = nAtoms[structNo];
    nal = nAlloc[structNo];
    if (nat>=nal)  {
      nal = nat+10;
      A2  = new PAtom[nal];
      GetVectorMemory ( id2,nal,0 );
      for (i=0;i<nat;i++)  {
        A2 [i] = A1 [i];
        id2[i] = id1[i];
      }
      for (i=nat;i<nal;i++)  {
        A2 [i] = NULL;
        id2[i] = 0;
      }
      if (A1)  delete[] A1;
      FreeVectorMemory ( id1,0 );
      atom[structNo] = A2;
      id  [structNo] = id2;
      A1  = A2;
      id1 = id2;
      nAlloc[structNo] = nal;
    }
    A1 [nat] = A;
    id1[nat] = atomid;
    nAtoms[structNo] = nat+1;
  }


  void  DeleteMContacts ( PPMContact & mcontact, int nContacts )  {
  int i;
    if (mcontact)  {
      for (i=0;i<nContacts;i++)
        if (mcontact[i])  delete mcontact[i];
      delete[] mcontact;
      mcontact = NULL;
    }
  }


  //  ====================   CoorManager   =====================

  CoorManager::CoorManager() : Root()  {
    InitMMDBCoorManager();
  }

  CoorManager::CoorManager ( io::RPStream Object ) : Root(Object)  {
    InitMMDBCoorManager();
  }

  CoorManager::~CoorManager()  {
    RemoveBricks ();
    RemoveMBricks();
  }

  void  CoorManager::ResetManager()  {
    Root::ResetManager();
    RemoveBricks       ();
    RemoveMBricks      ();
    InitMMDBCoorManager();
  }

  void  CoorManager::InitMMDBCoorManager()  {

    CoorIDCode  = CID_Ok;

    brick_size  = 6.0;  // angstroms
    xbrick_0    = 0.0;
    ybrick_0    = 0.0;
    zbrick_0    = 0.0;
    nbrick_x    = 0;
    nbrick_y    = 0;
    nbrick_z    = 0;
    brick       = NULL;

    mbrick_size = 6.0;  // angstroms
    xmbrick_0   = 0.0;
    ymbrick_0   = 0.0;
    zmbrick_0   = 0.0;
    nmbrick_x   = 0;
    nmbrick_y   = 0;
    nmbrick_z   = 0;
    mbrick      = NULL;

  }


  int  CoorManager::SetDefaultCoorID ( cpstr CID )  {
    return DefPath.SetPath ( CID );
  }

  PModel CoorManager::GetFirstDefinedModel()  {
  PModel mdl;
  int    i;
    mdl = NULL;
    for (i=0;(i<nModels) && (!mdl);i++)
      mdl = model[i];
    return mdl;
  }

  int CoorManager::GetFirstModelNum()  {
  PModel mdl;
  int    i;
    mdl = NULL;
    for (i=0;(i<nModels) && (!mdl);i++)
      mdl = model[i];
    if (mdl)  return mdl->GetSerNum();
    return 1;
  }


  PModel CoorManager::GetModel ( int modelNo )  {
    if ((modelNo>=1) && (modelNo<=nModels))
          return model[modelNo-1];
    else  return NULL;
  }

  PModel CoorManager::GetModel ( cpstr CID )  {
  int      modno,sn,rc;
  ChainID  chname;
  InsCode  ic;
  ResName  resname;
  AtomName aname;
  Element  elname;
  AltLoc   aloc;

    CoorIDCode = CID_Ok;
    rc = ParseAtomPath ( CID,modno,chname,sn,ic,resname,
                         aname,elname,aloc,&DefPath );
    if ((rc<0) || (rc & APATH_WC_ModelNo))  {
      CoorIDCode = CID_WrongPath;
      return NULL;
    }

    if ((modno>=1) && (modno<=nModels))
          return model[modno-1];
    else  return NULL;

  }

  void CoorManager::GetModelTable ( PPModel & modelTable,
                                         int & NumberOfModels )  {
    NumberOfModels = nModels;
    modelTable     = model;
  }

  int CoorManager::DeleteModel ( int modelNo )  {
    if ((modelNo>=1) && (modelNo<=nModels))  {
      if (model[modelNo-1])  {
        Exclude = false;
        delete model[modelNo-1];
        model[modelNo-1] = NULL;
        Exclude = true;
        return 1;
      }
    }
    return 0;
  }

  int CoorManager::DeleteModel ( cpstr CID )  {
  int      modno,sn,rc;
  ChainID  chname;
  InsCode  ic;
  ResName  resname;
  AtomName aname;
  Element  elname;
  AltLoc   aloc;

    CoorIDCode = CID_Ok;
    rc = ParseAtomPath ( CID,modno,chname,sn,ic,resname,
                         aname,elname,aloc,&DefPath );
    if ((rc<0) || (rc & APATH_WC_ModelNo))  {
      CoorIDCode = CID_WrongPath;
      return 0;
    }

    if ((modno>=1) && (modno<=nModels))  {
      if (model[modno-1])  {
        Exclude = false;
        delete model[modno-1];
        model[modno-1] = NULL;
        Exclude = true;
        return 1;
      }
    }

    return 0;

  }


  int CoorManager::DeleteSolvent()  {
  int i,k;
    Exclude = false;
    k = 0;
    for (i=0;i<nModels;i++)
      if (model[i])  {
        k += model[i]->DeleteSolvent();
        model[i]->TrimChainTable();
        if (model[i]->nChains<=0)  {
          delete model[i];
          model[i] = NULL;
        }
      }
    Exclude = true;
    return k;
  }


  //  ----------------  Adding/Inserting models  ---------------

  int  CoorManager::AddModel ( PModel mdl )  {
  PPModel model1;
  int     i,nnat,nat1;

    if (!mdl)
      return  nModels;

    for (i=0;i<nModels;i++)
      if (model[i]==mdl)  return -i;

    nnat = mdl->GetNumberOfAtoms ( true );
    AddAtomArray ( nnat );         // get space for new atoms

    if (mdl->GetCoordHierarchy())  {
      SwitchModel ( nModels+1 ); // get one more model at the end
      nat1 = nAtoms;
      model[nModels-1]->_copy ( mdl,atom,nat1 );
      model[nModels-1]->serNum = nModels;
      nAtoms = nat1;
    } else  {
      model1 = new PModel[nModels+1];
      for (i=0;i<nModels;i++)
        model1[i] = model[i];
      if (model)  delete[] model;
      model = model1;
      model[nModels] = mdl;
      model[nModels]->SetMMDBManager ( PManager(this),nModels+1 );
      model[nModels]->CheckInAtoms();
      nModels++;
    }

    return nModels;

  }

  int  CoorManager::InsModel ( PModel mdl, int modelNo )  {
    AddModel     ( mdl );
    RotateModels ( modelNo,nModels,1 );
    return nModels;
  }

  void CoorManager::RotateModels ( int modelNo1, int modelNo2,
                                   int rotdir )  {
  PModel mdl;
  PPAtom A;
  int    m1,m2,i11,i12,i21,i22,nat,i,k;

    m1 = IMax ( 0,modelNo1-1 );
    m2 = IMin ( nModels,modelNo2) - 1;
    if (m1>m2)  ISwap ( m1,m2 );

    if (m1!=m2)  {

      if (model[m1] && model[m2])  {
        model[m1]->GetAIndexRange ( i11,i12 );
        model[m2]->GetAIndexRange ( i21,i22 );
        if ((i11<i12) && (i21<i22) && (i12<i22))  {
          i11--;    i12--;
          i21--;    i22--;
          if (rotdir<0)  {
            //  rotate anticlockwise
            nat = i12-i11+1;
            A = new PAtom[nat];
            k = 0;
            for (i=i11;i<=i12;i++)
              A[k++] = atom[i];
            k = i11;
            for (i=i12+1;i<=i22;i++)  {
              atom[k] = atom[i];
              if (atom[k])  atom[k]->index = k+1;
              k++;
            }
            for (i=0;i<nat;i++)  {
              atom[k] = A[i];
              if (atom[k])  atom[k]->index = k+1;
              k++;
            }
          } else  {
            //  rotate anticlockwise
            nat = i22-i21+1;
            A = new PAtom[nat];
            k = 0;
            for (i=i21;i<=i22;i++)
              A[k++] = atom[i];
            k = i22;
            for (i=i21-1;i>=i11;i--)  {
              atom[k] = atom[i];
              if (atom[k])  atom[k]->index = k+1;
              k--;
            }
            for (i=nat-1;i>=0;i--)  {
              atom[k] = A[i];
              if (atom[k])  atom[k]->index = k+1;
              k--;
            }
          }
          delete[] A;
        }
      }

      if (rotdir<0)  {
        //  rotate anticlockwise
        mdl = model[m1];
        for (i=m1;i<m2;i++)  {
          model[i] = model[i+1];
          model[i]->serNum = i+1;
        }
        model[m2] = mdl;
        model[m2]->serNum = m2+1;
      } else  {
        //  rotate clockwise
        mdl = model[m2];
        for (i=m2;i>m1;i--)  {
          model[i] = model[i-1];
          model[i]->serNum = i+1;
        }
        model[m1] = mdl;
        model[m1]->serNum = m1+1;
      }

    }

  }


  void CoorManager::SwapModels ( int modelNo1, int modelNo2 )  {
  PModel mdl;
  PPAtom A;
  int    m1,m2,i11,i12,i21,i22,i,k,n;

    n = 0;  // tp depress "uninitialized" warning

    m1 = IMax ( 0,modelNo1-1 );
    m2 = IMin ( nModels,modelNo2) - 1;
    if (m1>m2)  ISwap ( m1,m2 );

    if (m1!=m2)  {

      if (model[m1])
        model[m1]->GetAIndexRange ( i11,i12 );
      else  {
        n = m1;
        while ((!model[n]) && (n<m2))  n++;
        if (n<m2)  {
          model[n]->GetAIndexRange ( i11,i12 );
          i12 = i11-1;
        } else
          n = -1;
      }

      if (n>=0)  {
        if (model[m2])
          model[m2]->GetAIndexRange ( i21,i22 );
        else  {
          n = m2;
          while ((!model[n]) && (m1<n))  n--;
          if (m1<n)  {
            model[n]->GetAIndexRange ( i21,i22 );
            i22 = i21-1;
          } else
            n = -1;
        }
      }

      if (n>=0)  {

        i11--;    i12--;
        i21--;    i22--;

        A = new PAtom[atmLen];
        k = 0;

        for (i=0     ;i<i11   ;i++)  A[k++] = atom[i];
        for (i=i21   ;i<=i22  ;i++)  A[k++] = atom[i];
        for (i=i12+1 ;i<i21   ;i++)  A[k++] = atom[i];
        for (i=i11   ;i<=i12  ;i++)  A[k++] = atom[i];

        for (i=0     ;i<nAtoms;i++)  if (A[i]) A[i]->index = i+1;
        for (i=nAtoms;i<atmLen;i++)  A[i] = NULL;

        if (atom)  delete[] atom;
        atom = A;

      }

      mdl       = model[m2];
      model[m2] = model[m1];
      model[m1] = mdl;

      model[m1]->serNum = m1+1;
      model[m2]->serNum = m2+1;

    }

  }


  PChain CoorManager::GetChain ( int modelNo, const ChainID chainID )  {
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->GetChain ( chainID );
    }
    return NULL;
  }

  PChain CoorManager::GetChain ( int modelNo, int chainNo )  {
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->GetChain ( chainNo );
    }
    return NULL;
  }

  PChain CoorManager::GetChain ( cpstr CID )  {
  int      modno,sn,rc;
  ChainID  chname;
  InsCode  ic;
  ResName  resname;
  AtomName aname;
  Element  elname;
  AltLoc   aloc;

    CoorIDCode = CID_Ok;
    rc = ParseAtomPath ( CID,modno,chname,sn,ic,resname,
                         aname,elname,aloc,&DefPath );
    if ((rc<0) || (rc & (APATH_WC_ModelNo | APATH_WC_ChainID)))  {
      CoorIDCode = CID_WrongPath;
      return NULL;
    }
    return GetChain ( modno,chname );

  }

  void  CoorManager::GetChainTable ( int modelNo,
                                          PPChain & chainTable,
                                          int & NumberOfChains )  {
    chainTable     = NULL;
    NumberOfChains = 0;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])  {
        chainTable     = model[modelNo-1]->chain;
        NumberOfChains = model[modelNo-1]->nChains;
      }
    }
  }

  void  CoorManager::GetChainTable ( cpstr CID,
                                          PPChain & chainTable,
                                          int & NumberOfChains )  {
  int      modno,sn,rc;
  ChainID  chname;
  InsCode  ic;
  ResName  resname;
  AtomName aname;
  Element  elname;
  AltLoc   aloc;

    chainTable     = NULL;
    NumberOfChains = 0;
    CoorIDCode     = CID_Ok;

    rc = ParseAtomPath ( CID,modno,chname,sn,ic,resname,
                         aname,elname,aloc,&DefPath );
    if ((rc<0) || (rc & APATH_WC_ModelNo))  {
      CoorIDCode = CID_WrongPath;
      return;
    }

    if ((0<modno) && (modno<=nModels))  {
      if (model[modno-1])  {
        chainTable     = model[modno-1]->chain;
        NumberOfChains = model[modno-1]->nChains;
      }
    }
  }


  int CoorManager::DeleteChain ( int modelNo, const ChainID chID )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteChain ( chID );
    }
    return 0;
  }

  int CoorManager::DeleteChain ( int modelNo, int chainNo )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteChain ( chainNo );
    }
    return 0;
  }

  int CoorManager::DeleteAllChains ( int modelNo )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAllChains();
    }
    return 0;
  }

  int CoorManager::DeleteAllChains()  {
  int i,k;
    k = 0;
    for (i=0;i<nModels;i++)
      if (model[i])  k += model[i]->DeleteAllChains();
    return k;
  }

  int CoorManager::AddChain ( int modelNo, PChain chain )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->AddChain ( chain );
    }
    return 0;
  }


  PResidue CoorManager::GetResidue ( int           modelNo,
                                     const ChainID chainID,
                                     int           seqNo,
                                     const InsCode insCode )  {
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->GetResidue ( chainID,seqNo,insCode );
    }
    return NULL;
  }

  PResidue CoorManager::GetResidue ( int modelNo, int chainNo,
                                     int seqNo,
                                     const InsCode insCode )  {
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->GetResidue ( chainNo,seqNo,insCode );
    }
    return NULL;
  }

  PResidue CoorManager::GetResidue ( int modelNo,
                                     const ChainID chainID,
                                     int resNo )  {
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->GetResidue ( chainID,resNo );
    }
    return NULL;
  }

  PResidue CoorManager::GetResidue ( int modelNo, int chainNo,
                                     int resNo )  {
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->GetResidue ( chainNo,resNo );
    }
    return NULL;
  }

  PResidue CoorManager::GetResidue ( cpstr CID )  {
  int      modno,sn,rc;
  ChainID  chname;
  InsCode  ic;
  ResName  resname;
  AtomName aname;
  Element  elname;
  AltLoc   aloc;

    CoorIDCode = CID_Ok;
    rc = ParseAtomPath ( CID,modno,chname,sn,ic,resname,
                         aname,elname,aloc,&DefPath );
    if ((rc<0) || (rc & (APATH_WC_ModelNo | APATH_WC_ChainID |
                         APATH_WC_SeqNum  | APATH_WC_InsCode)))  {
      CoorIDCode = CID_WrongPath;
      return NULL;
    }
    return GetResidue ( modno,chname,sn,ic );

  }


  int CoorManager::GetResidueNo ( int           modelNo,
                                  const ChainID chainID,
                                  int           seqNo,
                                  const InsCode insCode )  {
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->GetResidueNo ( chainID,seqNo,insCode );
    }
    return -3;
  }

  int CoorManager::GetResidueNo ( int modelNo, int chainNo,
                                  int seqNo,
                                  const InsCode insCode )  {
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->GetResidueNo ( chainNo,seqNo,insCode );
    }
    return -3;
  }

  void CoorManager::GetResidueTable ( PPResidue & resTable,
                                      int & NumberOfResidues )  {
  //    resTable has to be NULL or it will be reallocated. The
  //  application is responsible for deallocating the resTable (but not
  //  of its residues!). This does not apply to other GetResidueTable
  //  functions.
  PPChain   chain;
  PPResidue res;
  int       i,j,k,n,nChains,nResidues;

    if (resTable)  {
      delete[] resTable;
      resTable = NULL;
    }

    NumberOfResidues = 0;
    for (i=0;i<nModels;i++)
      if (model[i])  {
        model[i]->GetChainTable ( chain,nChains );
        for (j=0;j<model[i]->nChains;j++)
          if (chain[j])  {
            chain[j]->GetResidueTable ( res,nResidues );
            NumberOfResidues += nResidues;
          }
      }

    if (NumberOfResidues>0)  {
      resTable = new PResidue[NumberOfResidues];
      k = 0;
      for (i=0;i<nModels;i++)
        if (model[i])  {
          model[i]->GetChainTable ( chain,nChains );
          for (j=0;j<model[i]->nChains;j++)
            if (chain[j])  {
              chain[j]->GetResidueTable ( res,nResidues );
              for (n=0;n<nResidues;n++)
                if (res[n])  resTable[k++] = res[n];
            }
        }
      NumberOfResidues = k;
    }

  }

  void CoorManager::GetResidueTable ( int modelNo,
                                           const ChainID chainID,
                                           PPResidue & resTable,
                                           int & NumberOfResidues )  {
  PChain chain;
    resTable         = NULL;
    NumberOfResidues = 0;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])  {
        chain = model[modelNo-1]->GetChain ( chainID );
        if (chain)  {
          resTable         = chain->residue;
          NumberOfResidues = chain->nResidues;
        }
      }
    }
  }

  void CoorManager::GetResidueTable ( int modelNo, int chainNo,
                                           PPResidue & resTable,
                                           int & NumberOfResidues )  {
  PChain chain;
    resTable         = NULL;
    NumberOfResidues = 0;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])  {
        chain = model[modelNo-1]->GetChain ( chainNo );
        if (chain)  {
          resTable         = chain->residue;
          NumberOfResidues = chain->nResidues;
        }
      }
    }
  }

  void CoorManager::GetResidueTable ( cpstr CID,
                                           PPResidue & resTable,
                                           int & NumberOfResidues )  {
  int      modno,sn,rc;
  ChainID  chname;
  InsCode  ic;
  ResName  resname;
  AtomName aname;
  Element  elname;
  AltLoc   aloc;
  PChain  chain;

    resTable         = NULL;
    NumberOfResidues = 0;
    CoorIDCode       = CID_Ok;

    rc = ParseAtomPath ( CID,modno,chname,sn,ic,resname,
                         aname,elname,aloc,&DefPath );
    if ((rc<0) || (rc & (APATH_WC_ModelNo | APATH_WC_ChainID)))  {
      CoorIDCode = CID_WrongPath;
      return;
    }

    if ((0<modno) && (modno<=nModels))  {
      if (model[modno-1])  {
        chain = model[modno-1]->GetChain ( chname );
        if (chain)  {
          resTable         = chain->residue;
          NumberOfResidues = chain->nResidues;
        }
      }
    }

  }


  int CoorManager::DeleteResidue ( int           modelNo,
                                        const ChainID chainID,
                                        int           seqNo,
                                        const InsCode insCode )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteResidue ( chainID,seqNo,insCode );
    }
    return 0;
  }

  int CoorManager::DeleteResidue ( int           modelNo,
                                        const ChainID chainID,
                                        int           resNo )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteResidue ( chainID,resNo );
    }
    return 0;
  }

  int CoorManager::DeleteResidue ( int modelNo, int chainNo,
                                        int seqNo,
                                        const InsCode insCode )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteResidue ( chainNo,seqNo,insCode );
    }
    return 0;
  }

  int CoorManager::DeleteResidue ( int modelNo, int chainNo,
                                        int resNo )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteResidue ( chainNo,resNo );
    }
    return 0;
  }

  int CoorManager::DeleteAllResidues ( int modelNo,
                                            const ChainID chainID )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAllResidues ( chainID );
    }
    return 0;
  }

  int CoorManager::DeleteAllResidues ( int modelNo, int chainNo )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAllResidues ( chainNo );
    }
    return 0;
  }

  int CoorManager::DeleteAllResidues ( int modelNo )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAllResidues();
    }
    return 0;
  }

  int CoorManager::DeleteAllResidues()  {
  int i,k;
    k = 0;
    for (i=0;i<nModels;i++)
      if (model[i])  k += model[i]->DeleteAllResidues();
    return k;
  }

  int CoorManager::AddResidue ( int modelNo, const ChainID chainID,
                                     PResidue res )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->AddResidue ( chainID,res );
    }
    return 0;
  }

  int CoorManager::AddResidue ( int modelNo, int chainNo,
                                     PResidue res )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->AddResidue ( chainNo,res );
    }
    return 0;
  }



  int  CoorManager::GetNumberOfChains ( int modelNo )  {
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return  model[modelNo-1]->nChains;
    }
    return 0;
  }

  int  CoorManager::GetNumberOfChains ( cpstr CID )  {
  int      modno,sn,rc;
  ChainID  chname;
  InsCode  ic;
  ResName  resname;
  AtomName aname;
  Element  elname;
  AltLoc   aloc;

    CoorIDCode = CID_Ok;
    rc = ParseAtomPath ( CID,modno,chname,sn,ic,resname,
                         aname,elname,aloc,&DefPath );
    if ((rc<0) || (rc & APATH_WC_ModelNo))  {
      CoorIDCode = CID_WrongPath;
      return 0;
    }

    if ((0<modno) && (modno<=nModels))  {
      if (model[modno-1])
        return  model[modno-1]->nChains;
    }

    return 0;

  }

  int  CoorManager::GetNumberOfResidues ( int modelNo,
                                               const ChainID chainID )  {
  PChain chain;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])  {
        chain = model[modelNo-1]->GetChain ( chainID );
        if (chain)  return chain->nResidues;
      }
    }
    return 0;
  }

  int  CoorManager::GetNumberOfResidues ( int modelNo, int chainNo )  {
  PChain chain;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])  {
        if ((0<=chainNo) && (chainNo<model[modelNo-1]->nChains))  {
          chain = model[modelNo-1]->chain[chainNo];
          if (chain)  return chain->nResidues;
        }
      }
    }
    return 0;
  }

  int  CoorManager::GetNumberOfResidues ( cpstr CID )  {
  int      modno,sn,rc;
  ChainID  chname;
  InsCode  ic;
  ResName  resname;
  AtomName aname;
  Element  elname;
  AltLoc   aloc;
  PChain  chain;

    CoorIDCode = CID_Ok;
    rc = ParseAtomPath ( CID,modno,chname,sn,ic,resname,
                         aname,elname,aloc,&DefPath );
    if ((rc<0) || (rc & (APATH_WC_ModelNo | APATH_WC_ChainID)))  {
      CoorIDCode = CID_WrongPath;
      return 0;
    }

    if ((0<modno) && (modno<=nModels))  {
      if (model[modno-1])  {
        chain = model[modno-1]->GetChain ( chname );
        if (chain)  return chain->nResidues;
      }
    }

    return 0;

  }


  int  CoorManager::GetNumberOfAtoms ( int           modelNo,
                                            const ChainID chainID,
                                            int           seqNo,
                                            const InsCode insCode )  {
  PChain   chain;
  PResidue res;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])  {
        chain = model[modelNo-1]->GetChain ( chainID );
        if (chain)  {
          res = chain->GetResidue ( seqNo,insCode );
          if (res)  return res->nAtoms;
        }
      }
    }
    return 0;
  }

  int  CoorManager::GetNumberOfAtoms ( int modelNo, int chainNo,
                                            int seqNo,
                                            const InsCode insCode )  {
  PChain   chain;
  PResidue res;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])  {
        if ((0<=chainNo) && (chainNo<model[modelNo-1]->nChains))  {
          chain = model[modelNo-1]->chain[chainNo];
          if (chain)  {
            res = chain->GetResidue ( seqNo,insCode );
            if (res)  return res->nAtoms;
          }
        }
      }
    }
    return 0;
  }

  int  CoorManager::GetNumberOfAtoms ( int           modelNo,
                                            const ChainID chainID,
                                            int           resNo )  {
  PChain   chain;
  PResidue res;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])  {
        chain = model[modelNo-1]->GetChain ( chainID );
        if (chain)  {
          if ((0<=resNo) && (resNo<chain->nResidues))  {
            res = chain->residue[resNo];
            if (res)  return res->nAtoms;
          }
        }
      }
    }
    return 0;
  }

  int  CoorManager::GetNumberOfAtoms ( int modelNo, int chainNo,
                                            int resNo )  {
  PChain   chain;
  PResidue res;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])  {
        if ((0<=chainNo) && (chainNo<model[modelNo-1]->nChains))  {
          chain = model[modelNo-1]->chain[chainNo];
          if (chain)  {
            if ((0<=resNo) && (resNo<chain->nResidues))  {
              res = chain->residue[resNo];
              if (res)  return res->nAtoms;
            }
          }
        }
      }
    }
    return 0;
  }

  int  CoorManager::GetNumberOfAtoms ( cpstr CID )  {
  // returns number of atoms in residues identified by CID
  int       modno,sn,rc;
  ChainID   chname;
  InsCode   ic;
  ResName   resname;
  AtomName  aname;
  Element   elname;
  AltLoc    aloc;
  PChain   chain;
  PResidue res;

    CoorIDCode = CID_Ok;
    rc = ParseAtomPath ( CID,modno,chname,sn,ic,resname,
                         aname,elname,aloc,&DefPath );
    if ((rc<0) || (rc & (APATH_WC_ModelNo | APATH_WC_ChainID |
                         APATH_WC_SeqNum  | APATH_WC_InsCode)))  {
      CoorIDCode = CID_WrongPath;
      return 0;
    }

    if ((0<modno) && (modno<=nModels))  {
      if (model[modno-1])  {
        chain = model[modno-1]->GetChain ( chname );
        if (chain)  {
          res = chain->GetResidue ( sn,ic );
          if (res)  return res->nAtoms;
        }
      }
    }

    return 0;

  }


  // --------------------  Extracting atoms  -----------------------

  PAtom  CoorManager::GetAtom (
                     int            modelNo, // model serial number 1...
                     const ChainID  chID,    // chain ID
                     int            seqNo,   // residue sequence number
                     const InsCode  insCode, // residue insertion code
                     const AtomName aname,   // atom name
                     const Element  elmnt,   // chemical element code or '*'
                     const AltLoc   aloc     // alternate location indicator
                   )  {
  PModel   mdl;
  PChain   chn;
  PResidue res;
  PAtom    atm;

    if ((1<=modelNo) && (modelNo<=nModels))
          mdl = model[modelNo-1];
    else  mdl = NULL;
    if (!mdl)  {
      CoorIDCode = CID_NoModel;
      return NULL;
    }

    chn = mdl->GetChain ( chID );
    if (!chn)  {
      CoorIDCode = CID_NoChain;
      return NULL;
    }

    res = chn->GetResidue ( seqNo,insCode );
    if (!res)  {
      CoorIDCode = CID_NoResidue;
      return NULL;
    }

    atm = res->GetAtom ( aname,elmnt,aloc );
    if (!atm)  CoorIDCode = CID_NoAtom;
         else  CoorIDCode = CID_Ok;

    return atm;

  }

  PAtom CoorManager::GetAtom (
                       int           modelNo, // model serial number 1...
                       const ChainID chID,    // chain ID
                       int           seqNo,   // residue sequence number
                       const InsCode insCode, // residue insertion code
                       int           atomNo   // atom number 0..
                   )  {
  PModel   mdl;
  PChain   chn;
  PResidue res;
  PAtom    atm;

    if ((1<=modelNo) && (modelNo<=nModels))
          mdl = model[modelNo-1];
    else  mdl = NULL;
    if (!mdl)  {
      CoorIDCode = CID_NoModel;
      return NULL;
    }

    chn = mdl->GetChain ( chID );
    if (!chn)  {
      CoorIDCode = CID_NoChain;
      return NULL;
    }

    res = chn->GetResidue ( seqNo,insCode );
    if (!res)  {
      CoorIDCode = CID_NoResidue;
      return NULL;
    }

    if ((0<=atomNo) && (atomNo<res->nAtoms))
          atm = res->atom[atomNo];
    else  atm = NULL;
    if (!atm)  CoorIDCode = CID_NoAtom;
         else  CoorIDCode = CID_Ok;

    return atm;

  }

  PAtom CoorManager::GetAtom (
                       int            modelNo, // model serial number 1...
                       const ChainID  chID,    // chain ID
                       int            resNo,   // residue number 0..
                       const AtomName aname,   // atom name
                       const Element  elmnt,   // chemical element code or '*'
                       const AltLoc   aloc     // alternate location indicator
                   )  {
  PModel   mdl;
  PChain   chn;
  PResidue res;
  PAtom    atm;

    if ((1<=modelNo) && (modelNo<=nModels))
          mdl = model[modelNo-1];
    else  mdl = NULL;
    if (!mdl)  {
      CoorIDCode = CID_NoModel;
      return NULL;
    }

    chn = mdl->GetChain ( chID );
    if (!chn)  {
      CoorIDCode = CID_NoChain;
      return NULL;
    }

    if ((0<=resNo) && (resNo<chn->nResidues))
          res = chn->residue[resNo];
    else  res = NULL;
    if (!res)  {
      CoorIDCode = CID_NoResidue;
      return NULL;
    }

    atm = res->GetAtom ( aname,elmnt,aloc );
    if (!atm)  CoorIDCode = CID_NoAtom;
         else  CoorIDCode = CID_Ok;

    return atm;

  }

  PAtom CoorManager::GetAtom (
                       int           modelNo, // model serial number 1...
                       const ChainID chID,    // chain ID
                       int           resNo,   // residue number 0..
                       int           atomNo   // atom number 0..
                   )  {
  PModel   mdl;
  PChain   chn;
  PResidue res;
  PAtom    atm;

    if ((1<=modelNo) && (modelNo<=nModels))
          mdl = model[modelNo-1];
    else  mdl = NULL;
    if (!mdl)  {
      CoorIDCode = CID_NoModel;
      return NULL;
    }

    chn = mdl->GetChain ( chID );
    if (!chn)  {
      CoorIDCode = CID_NoChain;
      return NULL;
    }

    if ((0<=resNo) && (resNo<chn->nResidues))
          res = chn->residue[resNo];
    else  res = NULL;
    if (!res)  {
      CoorIDCode = CID_NoResidue;
      return NULL;
    }

    if ((0<=atomNo) && (atomNo<res->nAtoms))
          atm = res->atom[atomNo];
    else  atm = NULL;
    if (!atm)  CoorIDCode = CID_NoAtom;
         else  CoorIDCode = CID_Ok;

    return atm;

  }

  PAtom CoorManager::GetAtom (
                       int            modelNo, // model serial number 1...
                       int            chNo,    // chain number 0..
                       int            seqNo,   // residue sequence number
                       const InsCode  insCode, // residue insertion code
                       const AtomName aname,   // atom name
                       const Element  elmnt,   // chemical element code or '*'
                       const AltLoc   aloc     // alternate location indicator
                   )  {
  PModel   mdl;
  PChain   chn;
  PResidue res;
  PAtom    atm;

    if ((1<=modelNo) && (modelNo<=nModels))
          mdl = model[modelNo-1];
    else  mdl = NULL;
    if (!mdl)  {
      CoorIDCode = CID_NoModel;
      return NULL;
    }

    if ((0<=chNo) && (chNo<mdl->nChains))
          chn = mdl->chain[chNo];
    else  chn = NULL;
    if (!chn)  {
      CoorIDCode = CID_NoChain;
      return NULL;
    }

    res = chn->GetResidue ( seqNo,insCode );
    if (!res)  {
      CoorIDCode = CID_NoResidue;
      return NULL;
    }

    atm = res->GetAtom ( aname,elmnt,aloc );
    if (!atm)  CoorIDCode = CID_NoAtom;
         else  CoorIDCode = CID_Ok;

    return atm;

  }

  PAtom CoorManager::GetAtom (
                       int           modelNo, // model serial number 1...
                       int           chNo,    // chain number 0...
                       int           seqNo,   // residue sequence number
                       const InsCode insCode, // residue insertion code
                       int           atomNo   // atom number 0...
                   )  {
  PModel   mdl;
  PChain   chn;
  PResidue res;
  PAtom    atm;

    if ((1<=modelNo) && (modelNo<=nModels))
          mdl = model[modelNo-1];
    else  mdl = NULL;
    if (!mdl)  {
      CoorIDCode = CID_NoModel;
      return NULL;
    }

    if ((0<=chNo) && (chNo<mdl->nChains))
          chn = mdl->chain[chNo];
    else  chn = NULL;
    if (!chn)  {
      CoorIDCode = CID_NoChain;
      return NULL;
    }

    res = chn->GetResidue ( seqNo,insCode );
    if (!res)  {
      CoorIDCode = CID_NoResidue;
      return NULL;
    }

    if ((0<=atomNo) && (atomNo<res->nAtoms))
          atm = res->atom[atomNo];
    else  atm = NULL;
    if (!atm)  CoorIDCode = CID_NoAtom;
         else  CoorIDCode = CID_Ok;

    return atm;

  }

  PAtom CoorManager::GetAtom (
                       int            modelNo, // model serial number 1...
                       int            chNo,    // chain number 0...
                       int            resNo,   // residue number 0...
                       const AtomName aname,   // atom name
                       const Element  elmnt,   // chemical element code or '*'
                       const AltLoc   aloc     // alternate location indicator
                   )  {
  PModel   mdl;
  PChain   chn;
  PResidue res;
  PAtom    atm;

    if ((1<=modelNo) && (modelNo<=nModels))
          mdl = model[modelNo-1];
    else  mdl = NULL;
    if (!mdl)  {
      CoorIDCode = CID_NoModel;
      return NULL;
    }

    if ((0<=chNo) && (chNo<mdl->nChains))
          chn = mdl->chain[chNo];
    else  chn = NULL;
    if (!chn)  {
      CoorIDCode = CID_NoChain;
      return NULL;
    }

    if ((0<=resNo) && (resNo<chn->nResidues))
          res = chn->residue[resNo];
    else  res = NULL;
    if (!res)  {
      CoorIDCode = CID_NoResidue;
      return NULL;
    }

    atm = res->GetAtom ( aname,elmnt,aloc );
    if (!atm)  CoorIDCode = CID_NoAtom;
         else  CoorIDCode = CID_Ok;

    return atm;

  }

  PAtom CoorManager::GetAtom (
                       int modelNo, // model serial number 1...
                       int chNo,    // chain number 0...
                       int resNo,   // residue number 0...
                       int atomNo   // atom number 0...
                   )  {
  PModel   mdl;
  PChain   chn;
  PResidue res;
  PAtom    atm;

    if ((1<=modelNo) && (modelNo<=nModels))
          mdl = model[modelNo-1];
    else  mdl = NULL;
    if (!mdl)  {
      CoorIDCode = CID_NoModel;
      return NULL;
    }

    if ((0<=chNo) && (chNo<mdl->nChains))
          chn = mdl->chain[chNo];
    else  chn = NULL;
    if (!chn)  {
      CoorIDCode = CID_NoChain;
      return NULL;
    }

    if ((0<=resNo) && (resNo<chn->nResidues))
          res = chn->residue[resNo];
    else  res = NULL;
    if (!res)  {
      CoorIDCode = CID_NoResidue;
      return NULL;
    }

    if ((0<=atomNo) && (atomNo<res->nAtoms))
          atm = res->atom[atomNo];
    else  atm = NULL;
    if (!atm)  CoorIDCode = CID_NoAtom;
         else  CoorIDCode = CID_Ok;

    return atm;

  }


  PAtom CoorManager::GetAtom ( cpstr CID )  {
  int      modno,sn,rc;
  ChainID  chname;
  InsCode  ic;
  ResName  resname;
  AtomName aname;
  Element  elname;
  AltLoc   aloc;

    CoorIDCode = CID_Ok;
    rc = ParseAtomPath ( CID,modno,chname,sn,ic,resname,
                         aname,elname,aloc,&DefPath );
    if ((rc<0) || (rc & APATH_Incomplete))  {
      CoorIDCode = CID_WrongPath;
      return NULL;
    }

    return GetAtom ( modno,chname,sn,ic,aname,elname,aloc );

  }


  void CoorManager::GetAtomTable ( PPAtom & atomTable,
                                        int & NumberOfAtoms )  {
    atomTable     = atom;
    NumberOfAtoms = nAtoms;
  }

  void CoorManager::GetAtomTable ( int           modelNo,
                                        const ChainID chainID,
                                        int           seqNo,
                                        const InsCode insCode,
                                        PPAtom &     atomTable,
                                        int &         NumberOfAtoms )  {
  PResidue res;
    atomTable     = NULL;
    NumberOfAtoms = 0;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])  {
        res = model[modelNo-1]->GetResidue ( chainID,seqNo,insCode );
        if (res)  {
          atomTable     = res->atom;
          NumberOfAtoms = res->nAtoms;
        }
      }
    }
  }

  void CoorManager::GetAtomTable ( int           modelNo,
                                        int           chainNo,
                                        int           seqNo,
                                        const InsCode insCode,
                                        PPAtom &     atomTable,
                                        int &         NumberOfAtoms )  {
  PResidue res;
    atomTable     = NULL;
    NumberOfAtoms = 0;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])  {
        res = model[modelNo-1]->GetResidue ( chainNo,seqNo,insCode );
        if (res)  {
          atomTable     = res->atom;
          NumberOfAtoms = res->nAtoms;
        }
      }
    }
  }

  void CoorManager::GetAtomTable ( int           modelNo,
                                        const ChainID chainID,
                                        int           resNo,
                                        PPAtom &     atomTable,
                                        int &         NumberOfAtoms )  {
  PResidue res;
    atomTable     = NULL;
    NumberOfAtoms = 0;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])  {
        res = model[modelNo-1]->GetResidue ( chainID,resNo );
        if (res)  {
          atomTable     = res->atom;
          NumberOfAtoms = res->nAtoms;
        }
      }
    }
  }

  void CoorManager::GetAtomTable ( int modelNo, int chainNo,
                                        int resNo, PPAtom & atomTable,
                                        int & NumberOfAtoms )  {
  PResidue res;
    atomTable     = NULL;
    NumberOfAtoms = 0;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])  {
        res = model[modelNo-1]->GetResidue ( chainNo,resNo );
        if (res)  {
          atomTable     = res->atom;
          NumberOfAtoms = res->nAtoms;
        }
      }
    }
  }

  void CoorManager::GetAtomTable ( cpstr CID,
                                        PPAtom & atomTable,
                                        int & NumberOfAtoms )  {
  int       modno,sn,rc;
  ChainID   chname;
  InsCode   ic;
  ResName   resname;
  AtomName  aname;
  Element   elname;
  AltLoc    aloc;
  PResidue res;

    atomTable     = NULL;
    NumberOfAtoms = 0;

    CoorIDCode = CID_Ok;
    rc = ParseAtomPath ( CID,modno,chname,sn,ic,resname,
                         aname,elname,aloc,&DefPath );
    if ((rc<0) || (rc & (APATH_WC_ModelNo | APATH_WC_ChainID |
                         APATH_WC_SeqNum  | APATH_WC_InsCode)))  {
      CoorIDCode = CID_WrongPath;
      return;
    }

    res = GetResidue ( modno,chname,sn,ic );
    if (res)  {
      atomTable     = res->atom;
      NumberOfAtoms = res->nAtoms;
    }

  }


  void CoorManager::GetAtomTable1 ( PPAtom & atomTable,
                                         int & NumberOfAtoms )  {
  int i,j;
    if (atomTable)  delete[] atomTable;
    if (nAtoms>0)  {
      atomTable = new PAtom[nAtoms];
      j = 0;
      for (i=0;i<nAtoms;i++)
        if (atom[i])  {
          if (!atom[i]->Ter)
            atomTable[j++] = atom[i];
        }
      NumberOfAtoms = j;
    } else  {
      atomTable     = NULL;
      NumberOfAtoms = 0;
    }
  }

  void CoorManager::GetAtomTable1 ( int           modelNo,
                                         const ChainID chainID,
                                         int           seqNo,
                                         const InsCode insCode,
                                         PPAtom &     atomTable,
                                         int &         NumberOfAtoms )  {
  PResidue res;
    res = NULL;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        res = model[modelNo-1]->GetResidue ( chainID,seqNo,insCode );
    }
    if (res)
      res->GetAtomTable1 ( atomTable,NumberOfAtoms );
    else  {
      if (atomTable)  delete[] atomTable;
      atomTable     = NULL;
      NumberOfAtoms = 0;
    }
  }

  void CoorManager::GetAtomTable1 ( int           modelNo,
                                         int           chainNo,
                                         int           seqNo,
                                         const InsCode insCode,
                                         PPAtom &     atomTable,
                                         int &         NumberOfAtoms )  {
  PResidue res;
    res = NULL;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        res = model[modelNo-1]->GetResidue ( chainNo,seqNo,insCode );
    }
    if (res)
      res->GetAtomTable1 ( atomTable,NumberOfAtoms );
    else  {
      if (atomTable)  delete[] atomTable;
      atomTable     = NULL;
      NumberOfAtoms = 0;
    }
  }

  void CoorManager::GetAtomTable1 ( int           modelNo,
                                         const ChainID chainID,
                                         int           resNo,
                                         PPAtom &     atomTable,
                                         int &         NumberOfAtoms )  {
  PResidue res;
    res = NULL;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        res = model[modelNo-1]->GetResidue ( chainID,resNo );
    }
    if (res)
      res->GetAtomTable1 ( atomTable,NumberOfAtoms );
    else  {
      if (atomTable)  delete[] atomTable;
      atomTable     = NULL;
      NumberOfAtoms = 0;
    }
  }

  void CoorManager::GetAtomTable1 ( int modelNo, int chainNo,
                                         int resNo,
                                         PPAtom & atomTable,
                                         int & NumberOfAtoms )  {
  PResidue res;
    res = NULL;
    if ((0<modelNo) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        res = model[modelNo-1]->GetResidue ( chainNo,resNo );
    }
    if (res)
      res->GetAtomTable1 ( atomTable,NumberOfAtoms );
    else  {
      if (atomTable)  delete[] atomTable;
      atomTable     = NULL;
      NumberOfAtoms = 0;
    }
  }

  void CoorManager::GetAtomTable1 ( cpstr CID, PPAtom & atomTable,
                                         int & NumberOfAtoms )  {
  int       modno,sn,rc;
  ChainID   chname;
  InsCode   ic;
  ResName   resname;
  AtomName  aname;
  Element   elname;
  AltLoc    aloc;
  PResidue res;

    atomTable     = NULL;
    NumberOfAtoms = 0;

    CoorIDCode = CID_Ok;
    rc = ParseAtomPath ( CID,modno,chname,sn,ic,resname,
                         aname,elname,aloc,&DefPath );
    if ((rc<0) || (rc & (APATH_WC_ModelNo | APATH_WC_ChainID |
                         APATH_WC_SeqNum  | APATH_WC_InsCode)))  {
      CoorIDCode = CID_WrongPath;
      return;
    }

    res = GetResidue ( modno,chname,sn,ic );
    if (res)
      res->GetAtomTable1 ( atomTable,NumberOfAtoms );
    else  {
      if (atomTable)  delete[] atomTable;
      atomTable     = NULL;
      NumberOfAtoms = 0;
    }

  }



  int CoorManager::DeleteAtom ( int            modelNo,
                                     const ChainID  chID,
                                     int            seqNo,
                                     const InsCode  insCode,
                                     const AtomName aname,
                                     const Element  elmnt,
                                     const AltLoc   aloc )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAtom ( chID,seqNo,insCode,
                                              aname,elmnt,aloc );
    }
    return 0;
  }

  int CoorManager::DeleteAtom ( int           modelNo,
                                     const ChainID chID,
                                     int           seqNo,
                                     const InsCode insCode,
                                     int           atomNo )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAtom ( chID,seqNo,insCode,atomNo );
    }
    return 0;
  }

  int CoorManager::DeleteAtom ( int            modelNo,
                                     const ChainID  chID,
                                     int            resNo,
                                     const AtomName aname,
                                     const Element  elmnt,
                                     const AltLoc   aloc )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAtom ( chID,resNo,
                                              aname,elmnt,aloc );
    }
    return 0;
  }

  int CoorManager::DeleteAtom ( int           modelNo,
                                     const ChainID chID,
                                     int           resNo,
                                     int           atomNo )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAtom ( chID,resNo,atomNo );
    }
    return 0;
  }

  int CoorManager::DeleteAtom ( int modelNo, int chNo, int seqNo,
                                     const InsCode  insCode,
                                     const AtomName aname,
                                     const Element  elmnt,
                                     const AltLoc   aloc )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAtom ( chNo,seqNo,insCode,
                                              aname,elmnt,aloc );
    }
    return 0;
  }

  int CoorManager::DeleteAtom ( int modelNo, int chNo, int seqNo,
                                     const InsCode insCode, int atomNo )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAtom ( chNo,seqNo,insCode,atomNo );
    }
    return 0;
  }

  int CoorManager::DeleteAtom ( int modelNo, int chNo, int resNo,
                                     const AtomName aname,
                                     const Element  elmnt,
                                     const AltLoc   aloc )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAtom ( chNo,resNo,
                                              aname,elmnt,aloc );
    }
    return 0;
  }

  int CoorManager::DeleteAtom ( int modelNo, int chNo, int resNo,
                                     int atomNo )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAtom ( chNo,resNo,atomNo );
    }
    return 0;
  }

  int CoorManager::DeleteAllAtoms ( int           modelNo,
                                         const ChainID chID,
                                         int           seqNo,
                                         const InsCode insCode )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAllAtoms ( chID,seqNo,insCode );
    }
    return 0;
  }

  int CoorManager::DeleteAllAtoms ( int modelNo, const ChainID chID,
                                         int resNo )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAllAtoms ( chID,resNo );
    }
    return 0;
  }

  int CoorManager::DeleteAllAtoms ( int modelNo, const ChainID chID )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAllAtoms ( chID );
    }
    return 0;
  }

  int CoorManager::DeleteAllAtoms ( int modelNo, int chNo, int seqNo,
                                         const InsCode insCode )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAllAtoms ( chNo,seqNo,insCode );
    }
    return 0;
  }

  int CoorManager::DeleteAllAtoms ( int modelNo, int chNo,
                                         int resNo )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAllAtoms ( chNo,resNo );
    }
    return 0;
  }

  int CoorManager::DeleteAllAtoms ( int modelNo, int chNo )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAllAtoms ( chNo );
    }
    return 0;
  }

  int CoorManager::DeleteAllAtoms ( int modelNo )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->DeleteAllAtoms();
    }
    return 0;
  }

  int CoorManager::DeleteAllAtoms()  {
  int i,k;
    k = 0;
    for (i=0;i<nModels;i++)
      if (model[i])  k += model[i]->DeleteAllAtoms();
    return k;
  }


  /*
  int CoorManager::DeleteAltLocs()  {
  //  This function leaves only alternative location with maximal
  // occupancy, if those are equal or unspecified, the one with
  // "least" alternative location indicator.
  //  The function returns the number of deleted atoms and optimizes
  // the atom index.
  ChainID  chID;
  ResName  rName;
  InsCode  iCode;
  AtomName aname;
  AltLoc   aLoc,aL;
  realtype occupancy,occ;
  int      seqNum;
  int      i,j,k,i1,i2,n;

    k = 0;
    n = 0;
    i = 0;
    while (i<nAtoms)  {

      if (atom[i])  {
        seqNum    = atom[i]->GetSeqNum   ();
        occupancy = atom[i]->GetOccupancy();
        strcpy ( chID ,atom[i]->GetChainID() );
        strcpy ( rName,atom[i]->GetResName() );
        strcpy ( iCode,atom[i]->GetInsCode() );
        strcpy ( aname,atom[i]->name   );
        strcpy ( aLoc ,atom[i]->altLoc );
        j  = i+1;
        i1 = -1;
        i2 = i;
        while (j<nAtoms)
          if (atom[j])  {
            if ((atom[j]->GetSeqNum()==seqNum)         &&
                (!strcmp(atom[j]->name,aname))         &&
                (!strcmp(atom[j]->GetInsCode(),iCode)) &&
                (!strcmp(atom[j]->GetResName(),rName)) &&
                (!strcmp(atom[j]->GetChainID(),chID )))  {
              occ = atom[j]->GetOccupancy();
              if (occ>occupancy)  {
                occupancy = occ;
                i1 = j;
              }
              if (aLoc[0])  {
                strcpy ( aL,atom[j]->altLoc );
                if (!aL[0])  {
                  aLoc[0] = char(0);
                  i2 = j;
                } else if (strcmp(aL,aLoc)<0)  {
                  strcpy ( aLoc,aL );
                  i2 = j;
                }
              }
              j++;
            } else
              break;
          } else
            j++;
        if (i1<0)  {
          if (atom[i]->WhatIsSet & ASET_Occupancy)  i1 = i;
                                              else  i1 = i2;
        }
        while (i<j)  {
          if (atom[i])  {
            if (i!=i1)  {
              delete atom[i];
              atom[i] = NULL;
              n++;
            } else  {
              if (k<i)  {
                atom[k] = atom[i];
                atom[k]->index = k+1;
              }
              k++;
            }
          }
          i++;
        }

      } else
        i++;

    }

    nAtoms = k;
    return n;

  }
  */

  int CoorManager::DeleteAltLocs()  {
  //  This function leaves only alternative location with maximal
  // occupancy, if those are equal or unspecified, the one with
  // "least" alternative location indicator.
  //  The function returns the number of deleted atoms. All tables
  // remain untrimmed, so that explicit trimming or calling
  // FinishStructEdit() at some point is required.
  int i,n;

    n = 0;
    for (i=0;i<nModels;i++)
      if (model[i])  n += model[i]->DeleteAltLocs();

    return n;

  }

  int CoorManager::AddAtom ( int           modelNo,
                                  const ChainID chID,
                                  int           seqNo,
                                  const InsCode insCode,
                                  PAtom atom )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->AddAtom ( chID,seqNo,insCode,atom );
    }
    return 0;
  }

  int CoorManager::AddAtom ( int modelNo, const ChainID chID,
                                  int resNo, PAtom atom )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->AddAtom ( chID,resNo,atom );
    }
    return 0;
  }

  int CoorManager::AddAtom ( int modelNo, int chNo,
                                  int seqNo, const InsCode insCode,
                                  PAtom atom )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->AddAtom ( chNo,seqNo,insCode,atom );
    }
    return 0;
  }

  int CoorManager::AddAtom ( int modelNo, int chNo,
                                  int resNo, PAtom atom )  {
    if ((modelNo>0) && (modelNo<=nModels))  {
      if (model[modelNo-1])
        return model[modelNo-1]->AddAtom ( chNo,resNo,atom );
    }
    return 0;
  }


  void  CoorManager::RemoveBricks()  {
  int i,j,k;
    if (brick)  {
      for (i=0;i<nbrick_x;i++)
        if (brick[i])  {
          for (j=0;j<nbrick_y;j++)
            if (brick[i][j])  {
              for (k=0;k<nbrick_z;k++)
                if (brick[i][j][k])  delete brick[i][j][k];
              delete[] brick[i][j];
            }
          delete[] brick[i];
        }
      delete[] brick;
    }
    brick    = NULL;
    nbrick_x = 0;
    nbrick_y = 0;
    nbrick_z = 0;
  }

  void  CoorManager::GetBrickCoor ( PAtom A,
                                    int & nx, int & ny, int & nz ) {
    nx = (int)floor((A->x-xbrick_0)/brick_size);
    ny = (int)floor((A->y-ybrick_0)/brick_size);
    nz = (int)floor((A->z-zbrick_0)/brick_size);
    if ((ny<0) || (nz<0) || (nx>=nbrick_x) ||
        (ny>=nbrick_y) || (nz>=nbrick_z))  nx = -1;
  }
  
  void  CoorManager::GetBrickCoor ( realtype x, realtype y,
                                    realtype z, int & nx,
                                    int & ny, int & nz ) {
    nx = (int)floor((x-xbrick_0)/brick_size);
    ny = (int)floor((y-ybrick_0)/brick_size);
    nz = (int)floor((z-zbrick_0)/brick_size);
    if ((ny<0) || (nz<0) || (nx>=nbrick_x) ||
        (ny>=nbrick_y) || (nz>=nbrick_z))  nx = -1;
  }
  
  void  CoorManager::GetBrickCoor ( vect3 & xyz,
                                    int & nx, int & ny, int & nz ) {
    nx = (int)floor((xyz[0]-xbrick_0)/brick_size);
    ny = (int)floor((xyz[1]-ybrick_0)/brick_size);
    nz = (int)floor((xyz[2]-zbrick_0)/brick_size);
    if ((ny<0) || (nz<0) || (nx>=nbrick_x) ||
        (ny>=nbrick_y) || (nz>=nbrick_z))  nx = -1;
  }

  void  CoorManager::GetBrickDimension (
                        int & nxmax, int & nymax, int & nzmax )  {
    if (!brick)  {
      nxmax = 0;  nymax = 0;  nzmax = 0;
    } else  {
      nxmax = nbrick_x;
      nymax = nbrick_y;
      nzmax = nbrick_z;
    }
  }

  PBrick CoorManager::GetBrick ( int nx, int ny, int nz )  {
    if (!brick)  return NULL;
    if ((nx>=0) && (nx<nbrick_x) &&
        (ny>=0) && (ny<nbrick_y) &&
        (nz>=0) && (nz<nbrick_z))  {
      if (!brick[nx])     return NULL;
      if (!brick[nx][ny]) return NULL;
      return brick[nx][ny][nz];
    }
    return NULL;
  }

  void  CoorManager::MakeBricks ( PPAtom   atmvec,
                                  int      avlen,
                                  realtype Margin,
                                  realtype BrickSize )  {
  //    Makes bricking for atoms contained in vector atmvec of length
  // avlen, with brick size BrickSize (in angstroms). The previous
  // bricking, if there was any, is removed.
  int      i,j, nx,ny,nz, alen;
  realtype x1,x2, y1,y2, z1,z2, dx,dy,dz;
  PPAtom   A;
  
    RemoveBricks();

    brick_size = BrickSize;

    if (atmvec)  {
      A    = atmvec;
      alen = avlen;
    } else  {
      A    = atom;
      alen = nAtoms;
    }

    if (alen>0)  {
      //  find the range of coordinates
      x1 = MaxReal;
      x2 = -x1;
      y1 = MaxReal;
      y2 = -y1;
      z1 = MaxReal;
      z2 = -z1;
      for (i=0;i<alen;i++)
        if (A[i])  {
          if ((!A[i]->Ter) && (A[i]->WhatIsSet & ASET_Coordinates))  {
            if (A[i]->x<x1)  x1 = A[i]->x;
            if (A[i]->x>x2)  x2 = A[i]->x;
            if (A[i]->y<y1)  y1 = A[i]->y;
            if (A[i]->y>y2)  y2 = A[i]->y;
            if (A[i]->z<z1)  z1 = A[i]->z;
            if (A[i]->z>z2)  z2 = A[i]->z;
          }
        }
      if (x1<MaxReal)  {
        x1 -= Margin; x2 += Margin;
        y1 -= Margin; y2 += Margin;
        z1 -= Margin; z2 += Margin;
        dx = x2-x1;  nbrick_x = mround(dx/brick_size+0.0001)+1;
        dy = y2-y1;  nbrick_y = mround(dy/brick_size+0.0001)+1;
        dz = z2-z1;  nbrick_z = mround(dz/brick_size+0.0001)+1;
        xbrick_0 = x1 - (nbrick_x*brick_size-dx)/2.0;
        ybrick_0 = y1 - (nbrick_y*brick_size-dy)/2.0;
        zbrick_0 = z1 - (nbrick_z*brick_size-dz)/2.0;
        for (i=0;i<alen;i++)
          if (A[i])  {
            if ((!A[i]->Ter) && (A[i]->WhatIsSet & ASET_Coordinates))  {
              GetBrickCoor ( A[i],nx,ny,nz );
              if (nx>=0)  {
                if (!brick)  {
                  brick = new PPPBrick[nbrick_x];
                  for (j=0;j<nbrick_x;j++)
                    brick[j] = NULL;
                }
                if (!brick[nx])  {
                  brick[nx] = new PPBrick[nbrick_y];
                  for (j=0;j<nbrick_y;j++)
                    brick[nx][j] = NULL;
                }
                if (!brick[nx][ny])  {
                  brick[nx][ny] = new PBrick[nbrick_z];
                  for (j=0;j<nbrick_z;j++)
                    brick[nx][ny][j] = NULL;
                }
                if (!brick[nx][ny][nz])
                  brick[nx][ny][nz] = new Brick();
                brick[nx][ny][nz]->AddAtom ( A[i],i );
              } else
                printf ( " error in "
                         "CoorManager::MakeBricks!!!\n" );
            }
          }
      }
    }

  }


  void  CoorManager::RemoveMBricks()  {
  int i,j,k;
    if (mbrick)  {
      for (i=0;i<nmbrick_x;i++)
        if (mbrick[i])  {
          for (j=0;j<nmbrick_y;j++)
            if (mbrick[i][j])  {
              for (k=0;k<nmbrick_z;k++)
                if (mbrick[i][j][k])  delete mbrick[i][j][k];
              delete[] mbrick[i][j];
            }
          delete[] mbrick[i];
        }
      delete[] mbrick;
    }
    mbrick    = NULL;
    nmbrick_x = 0;
    nmbrick_y = 0;
    nmbrick_z = 0;
  }

  void  CoorManager::GetMBrickCoor ( PAtom A,
                                       int & nx, int & ny, int & nz )  {
    nx = (int)floor((A->x-xmbrick_0)/mbrick_size);
    ny = (int)floor((A->y-ymbrick_0)/mbrick_size);
    nz = (int)floor((A->z-zmbrick_0)/mbrick_size);
    if ((ny<0) || (nz<0) || (nx>=nmbrick_x) ||
        (ny>=nmbrick_y)  || (nz>=nmbrick_z))  nx = -nx-1;
  }

  void  CoorManager::GetMBrickCoor (
                                  realtype x, realtype y, realtype z,
                                  int   & nx, int   & ny, int   & nz )  {
    nx = (int)floor((x-xmbrick_0)/mbrick_size);
    ny = (int)floor((y-ymbrick_0)/mbrick_size);
    nz = (int)floor((z-zmbrick_0)/mbrick_size);
    if ((ny<0) || (nz<0) || (nx>=nmbrick_x) ||
          (ny>=nmbrick_y) || (nz>=nmbrick_z))  nx = -nx-1;
  }

  void  CoorManager::GetMBrickDimension (
                        int & nxmax, int & nymax, int & nzmax )  {
    if (!brick)  {
      nxmax = 0;  nymax = 0;  nzmax = 0;
    } else  {
      nxmax = nmbrick_x;
      nymax = nmbrick_y;
      nzmax = nmbrick_z;
    }
  }

  PMBrick CoorManager::GetMBrick ( int nx, int ny, int nz )  {
    if (!mbrick)  return NULL;
    if ((nx>=0) && (nx<nmbrick_x) &&
        (ny>=0) && (ny<nmbrick_y) &&
        (nz>=0) && (nz<nmbrick_z))  {
      if (!mbrick[nx])     return NULL;
      if (!mbrick[nx][ny]) return NULL;
      return mbrick[nx][ny][nz];
    }
    return NULL;
  }

  void  CoorManager::MakeMBricks ( PPAtom * atmvec,  ivector avlen,
                                   int nStructures, realtype Margin,
                                   realtype BrickSize )  {
  //    Makes bricking for atoms contained in vectors atmvec with lengths
  // given in avlen, with brick size BrickSize (in angstroms).
  // The previous bricking, if there was any, is removed.
  int      i,j,k, nx,ny,nz;
  realtype x1,x2, y1,y2, z1,z2, dx,dy,dz;
  PAtom   A;

    RemoveMBricks();

    mbrick_size = BrickSize;

    //  find the range of coordinates
    x1 = MaxReal;
    x2 = -x1;
    y1 = MaxReal;
    y2 = -y1;
    z1 = MaxReal;
    z2 = -z1;
    for (i=0;i<nStructures;i++)
      for (j=0;j<avlen[i];j++)  {
        A = atmvec[i][j];
        if (A)  {
          if ((!A->Ter) && (A->WhatIsSet & ASET_Coordinates))  {
            if (A->x<x1)  x1 = A->x;
            if (A->x>x2)  x2 = A->x;
            if (A->y<y1)  y1 = A->y;
            if (A->y>y2)  y2 = A->y;
            if (A->z<z1)  z1 = A->z;
            if (A->z>z2)  z2 = A->z;
          }
        }
      }
    if (x1<MaxReal)  {
      x1 -= Margin; x2 += Margin;
      y1 -= Margin; y2 += Margin;
      z1 -= Margin; z2 += Margin;
      dx = x2-x1;  nmbrick_x = mround(dx/mbrick_size+0.0001)+1;
      dy = y2-y1;  nmbrick_y = mround(dy/mbrick_size+0.0001)+1;
      dz = z2-z1;  nmbrick_z = mround(dz/mbrick_size+0.0001)+1;
      xmbrick_0 = x1 - (nmbrick_x*mbrick_size-dx)/2.0;
      ymbrick_0 = y1 - (nmbrick_y*mbrick_size-dy)/2.0;
      zmbrick_0 = z1 - (nmbrick_z*mbrick_size-dz)/2.0;
      /*
      mbrick = new PPPMBrick[nmbrick_x];
      for (i=0;i<nmbrick_x;i++)  {
        mbrick[i] = new PPMBrick[nmbrick_y];
        for (j=0;j<nmbrick_y;j++)  {
          mbrick[i][j] = new PMBrick[nmbrick_z];
          for (k=0;k<nmbrick_z;k++)
            mbrick[i][j][k] = new mbrick(nStructures);
        }
      }
      */
      for (i=0;i<nStructures;i++)
        for (j=0;j<avlen[i];j++)  {
          A = atmvec[i][j];
          if (A)  {
            if ((!A->Ter) && (A->WhatIsSet & ASET_Coordinates))  {
              GetMBrickCoor ( A,nx,ny,nz );
              if (nx>=0)  {
                if (!mbrick)  {
                  mbrick = new PPPMBrick[nmbrick_x];
                  for (k=0;k<nmbrick_x;k++)
                    mbrick[k] = NULL;
                }
                if (!mbrick[nx])  {
                  mbrick[nx] = new PPMBrick[nmbrick_y];
                  for (k=0;k<nmbrick_y;k++)
                    mbrick[nx][k] = NULL;
                }
                if (!mbrick[nx][ny])  {
                  mbrick[nx][ny] = new PMBrick[nmbrick_z];
                  for (k=0;k<nmbrick_z;k++)
                    mbrick[nx][ny][k] = NULL;
                }
                if (!mbrick[nx][ny][nz])
                  mbrick[nx][ny][nz] = new MBrick(nStructures);
                mbrick[nx][ny][nz]->AddAtom ( A,i,j );
              } else
                printf ( " error in "
                         "CoorManager::MakeMBricks!!!\n" );
            }
          }
        }
    }

  }


  int  CoorManager::GenerateSymMates ( PGenSym genSym )  {
  //
  //   The function generates symmetry mates according to symmetry
  // operations found in GenSym. Results of first symmetry operation
  // (number 0) always replaces the existing set of atoms, others
  // are added as additional sets.
  //   If GenSym is set to NULL, the function generates all
  // symmetry mates for the unit cell taking the symmetry information
  // from cryst.SymOps.
  //   The newly generated chains are added to each model. These
  // chains have many-character chain names, composed as 'x_n',
  // where 'x' is the original name and 'n' is a unique number, which
  // coincides with the symmetry operation (order) number; number '_0'
  // (for the very first symmetry operatyion) is missing. Another
  // side effect is the disorder in atoms' serial numbers.
  //   The hierarchy should therefore be cleaned after
  // generating the symmetry mates. An appropriate way to do
  // that is to issue the following call:
  //
  //   PDBCleanup ( PDBCLEAN_TER | PDBCLEAN_ALTCODE_STRONG |
  //                PDBCLEAN_CHAIN_STRONG | PDBCLEAN_SERIAL );
  //
  PPCoorManager Mate;
  int           i,j,k,n,nMates,nMates1,nAtoms1;
  PPAtom        Atom1;
  PPModel       Model1;

    if (genSym)  nMates = genSym->GetNofSymOps();
           else  nMates = cryst.GetNumberOfSymOps();
    if (nMates<=0)  return GSM_NoSymOps;

    if (!cryst.areMatrices())       return GSM_NoTransfMatrices;
    if (!cryst.isCellParameters())  return GSM_NoCell;

    nMates1 = nMates-1;
    if (nMates1>0)  {

      //  Generate symmetry mates in parallel hierarchies
      Mate = new PCoorManager[nMates1];
      for (i=0;i<nMates1;i++)  {
        Mate[i] = new CoorManager();
        Mate[i]->Copy ( this );
        Mate[i]->ApplySymTransform ( i+1,genSym );
      }

      //  apply 1st symmetry operation:
      if (genSym)  ApplySymTransform ( 0,genSym );

      //  Gather all symmetry mates in 'this' hierarchy
      nAtoms1 = nMates*nAtoms;        // new number of atoms
      Atom1   = new PAtom[nAtoms1];  // new array of atoms

      if (nModels>0)  Model1 = new PModel[nModels];  // new array of
                else  Model1 = NULL;                  // models

      k = 0;  // index of collected atoms
      for (i=0;i<nModels;i++)
        if (model[i])  {
          Model1[i] = newModel();
          Model1[i]->SetMMDBManager ( PManager(this),i+1 );
          for (j=0;j<model[i]->nChains;j++)
            Model1[i]->MoveChain ( model[i]->chain[j],atom,Atom1,k,0 );
          for (n=0;n<nMates1;n++)
            for (j=0;j<model[i]->nChains;j++)
              Model1[i]->MoveChain ( Mate[n]->model[i]->chain[j],
                                     Mate[n]->atom,Atom1,k,n+1 );
        } else
          Model1[i] = NULL;

      if (model) delete[] model;
      model = Model1;

      for (i=0;i<nAtoms;i++)
        if (atom[i])  delete atom[i];  // should never happen
      if (atom)  delete[] atom;
      atom   = Atom1;
      atmLen = nAtoms1;   // length of Atom array
      nAtoms = k;

      //  Dispose parallel hierarchies
      for (i=0;i<nMates1;i++)
        delete Mate[i];
      delete[] Mate;

    } else  {
      //  just apply the only symmetry operation:
      if (genSym)  ApplySymTransform ( 0,genSym );
    }

    return GSM_Ok;

  }

  void  CoorManager::ApplyTransform ( const mat44 & TMatrix )  {
  // simply transforms all coordinates by multiplying with matrix TMatrix
  int i;
    for (i=0;i<nAtoms;i++)
      if (atom[i])  {
        if (!atom[i]->Ter)  atom[i]->Transform ( TMatrix );
      }
  }

  void  CoorManager::ApplySymTransform ( int SymOpNo, PGenSym genSym ) {
  //    This procedure applies the symmetry operation number SymOpNo
  // (starting from 0 on) and renames chains as specified in
  // GenSym.
  //    The chains don't have to be renamed. The number of chains
  // to be renamed is obtained as GenSym->nChains[SymOpNo], their
  // old names - as GenSym->chID1[SymOpNo][j], and their new names
  // - as GenSym->chID2[SymOpNo][j],  0<=j<GenSym->nChains[SymOpNo].
  mat44   tmat;
  int     i,j,k,nChn;
  PPChain chain;
    if (cryst.GetTMatrix(tmat,SymOpNo,0,0,0,PSymOps(genSym))
         ==SYMOP_Ok)  {
      for (i=0;i<nAtoms;i++)
        if (atom[i])  {
          if (!atom[i]->Ter)  atom[i]->Transform ( tmat );
        }
      if (genSym)
        for (i=0;i<nModels;i++)
          if (model[i])  {
            model[i]->GetChainTable ( chain,nChn );
            for (j=0;j<genSym->nChains[SymOpNo];j++)
              for (k=0;k<nChn;k++)
                if (!strcmp(chain[k]->chainID,genSym->chID1[SymOpNo][j]))
                  chain[k]->SetChainID ( genSym->chID2[SymOpNo][j] );
          }
    }
  }


  void  GetEulerRotMatrix ( mat33 & erm,
                            realtype alpha,
                            realtype beta,
                            realtype gamma )  {
  //  Calculates the Euler rotation matrix for rotation:
  //                   1) about z-axis by angle alpha (0..2*Pi)
  //                   2) about new y-axis by angle beta (0..Pi)
  //                   3) about new z-axis by angle gamma (0..2*Pi)
  realtype ca,cb,cg, sa,sb,sg;

    ca = cos(alpha);
    sa = sin(alpha);
    cb = cos(beta);
    sb = sin(beta);
    cg = cos(gamma);
    sg = sin(gamma);

    erm[0][0] =  ca*cb*cg - sa*sg;
    erm[0][1] =  cb*cg*sa + ca*sg;
    erm[0][2] = -cg*sb;

    erm[1][0] = -cg*sa - ca*cb*sg;
    erm[1][1] =  ca*cg - cb*sa*sg;
    erm[1][2] =  sb*sg;

    erm[2][0] =  ca*sb;
    erm[2][1] =  sa*sb;
    erm[2][2] =  cb;

  }



  void  GetEulerTMatrix ( mat44 & erm,
                          realtype alpha,
                          realtype beta,
                          realtype gamma,
                          realtype x0,
                          realtype y0,
                          realtype z0 )  {
  //  Calculates the Euler rotation-translation matrix for rotation:
  //                   1) about z-axis by angle alpha
  //                   2) about new y-axis by angle beta
  //                   3) about new z-axis by angle gamma
  //  Point (x0,y0,z0) is the center of rotation.
  mat33 m;

    m[0][0] = 1.0;
    GetEulerRotMatrix ( m,alpha,beta,gamma );

    erm[0][0] = m[0][0];  erm[0][1] = m[0][1];  erm[0][2] = m[0][2];
    erm[1][0] = m[1][0];  erm[1][1] = m[1][1];  erm[1][2] = m[1][2];
    erm[2][0] = m[2][0];  erm[2][1] = m[2][1];  erm[2][2] = m[2][2];
    erm[3][0] = 0.0;      erm[3][1] = 0.0;      erm[3][2] = 0.0;

    erm[3][3] = 1.0;

    erm[0][3] = x0 - m[0][0]*x0 - m[0][1]*y0 - m[0][2]*z0;
    erm[1][3] = y0 - m[1][0]*x0 - m[1][1]*y0 - m[1][2]*z0;
    erm[2][3] = z0 - m[2][0]*x0 - m[2][1]*y0 - m[2][2]*z0;

  }


  void  EulerRotation ( PPAtom  A,
                        int      nA,
                        realtype alpha,
                        realtype beta,
                        realtype gamma,
                        realtype x0,
                        realtype y0,
                        realtype z0 )  {
  //  Euler rotation:  1) about z-axis by angle alpha
  //                   2) about new y-axis by angle beta
  //                   3) about new z-axis by angle gamma
  //  Point (x0,y0,z0) is the center of rotation.
  mat33    m;
  realtype x,y,z;
  int      i;

    m[0][0] = 1.0;
    GetEulerRotMatrix ( m,alpha,beta,gamma );

    for (i=0;i<nA;i++)
      if (A[i])  {
        if ((!A[i]->Ter) && (A[i]->WhatIsSet & ASET_Coordinates))  {
          x = A[i]->x - x0;
          y = A[i]->y - y0;
          z = A[i]->z - z0;
          A[i]->x = m[0][0]*x + m[0][1]*y + m[0][2]*z + x0;
          A[i]->y = m[1][0]*x + m[1][1]*y + m[1][2]*z + y0;
          A[i]->z = m[2][0]*x + m[2][1]*y + m[2][2]*z + z0;
        }
      }

  }


  void  GetVecRotMatrix ( mat33 & vrm,
                          realtype alpha,
                          realtype vx,
                          realtype vy,
                          realtype vz )  {
  //   Calculates the rotation matrix for rotation by angle alpha about
  // arbitrary vector directed as (vx,vy,vz) = (vx2-vx1,vy2-vy1,vz2-vz1).
  realtype ca,sa, rx,ry,rz, vl;

    ca = cos(alpha);
    sa = sin(alpha);

    vl = sqrt ( vx*vx + vy*vy + vz*vz );
    if (vl<=0.0)  return;
    rx = vx/vl;
    ry = vy/vl;
    rz = vz/vl;

    vrm[0][0] = rx*rx*(1.0-ca) + ca;
    vrm[0][1] = rx*ry*(1.0-ca) - rz*sa;
    vrm[0][2] = rx*rz*(1.0-ca) + ry*sa;

    vrm[1][0] = ry*rx*(1.0-ca) + rz*sa;
    vrm[1][1] = ry*ry*(1.0-ca) + ca;
    vrm[1][2] = ry*rz*(1.0-ca) - rx*sa;

    vrm[2][0] = rz*rx*(1.0-ca) - ry*sa;
    vrm[2][1] = rz*ry*(1.0-ca) + rx*sa;
    vrm[2][2] = rz*rz*(1.0-ca) + ca;

  }


  void  GetRotParameters ( mat33    & vrm,
                           realtype & alpha,
                           realtype & vx,
                           realtype & vy,
                           realtype & vz )  {
  //    Given the rotation matrix vrm, GetRotParameters(..)
  //  returns the rotation angle alpha and the normalized
  //  rotation axis vector (vx,vy,vz).
  //    The rotation angle and vector are determined up to
  //  their sign (however correlated, so that being substituted
  //  into GetVecRotMatrix(..) they yield the same rotation
  //  matrix).
  //    The function does not check for vrm to be a valid
  //  rotation matrix.
  realtype ca,sa,vl;
    ca = (vrm[0][0]+vrm[1][1]+vrm[2][2]-1.0)/2.0;
    if (ca<-1.0) ca = -1.0;  // does not work if rotation
    if (ca>1.0)  ca =  1.0;  //   matrix is correct
    sa = sqrt(1.0-ca*ca);
    if (sa>0.0)  {
      alpha = acos(ca);
      // coefficient of 2 is corrected by normalization below
      vx    = (vrm[2][1]-vrm[1][2])/sa;
      vy    = (vrm[0][2]-vrm[2][0])/sa;
      vz    = (vrm[1][0]-vrm[0][1])/sa;
      // the following code is formally redundant if rotation
      // matrix is correct, however it eliminates the round-offs
      vl    = sqrt(vx*vx+vy*vy+vz*vz);
      vx   /= vl;
      vy   /= vl;
      vz   /= vl;
    } else  {
      // zero rotation, arbitrary axis would do
      alpha = 0.0;
      vx    = 1.0;
      vy    = 0.0;
      vz    = 0.0;
    }
  }


  void  GetVecTMatrix ( mat44 & vrm,
                        realtype alpha,
                        realtype vx,
                        realtype vy,
                        realtype vz,
                        realtype x0,
                        realtype y0,
                        realtype z0 )  {
  //   Calculates the rotation-translation matrix for rotation by angle
  // alpha about arbitrary vector directed as (vx,vy,vz) =
  // (vx2-vx1,vy2-vy1,vz2-vz1). Point (x0,y0,z0) is the center of
  // rotation -- actually a point belonging to the rotation axis.
  mat33 m;

    GetVecRotMatrix ( m,alpha,vx,vy,vz );

    vrm[0][0] = m[0][0];  vrm[0][1] = m[0][1];  vrm[0][2] = m[0][2];
    vrm[1][0] = m[1][0];  vrm[1][1] = m[1][1];  vrm[1][2] = m[1][2];
    vrm[2][0] = m[2][0];  vrm[2][1] = m[2][1];  vrm[2][2] = m[2][2];
    vrm[3][0] = 0.0;      vrm[3][1] = 0.0;      vrm[3][2] = 0.0;

    vrm[3][3] = 1.0;

    vrm[0][3] = x0 - m[0][0]*x0 - m[0][1]*y0 - m[0][2]*z0;
    vrm[1][3] = y0 - m[1][0]*x0 - m[1][1]*y0 - m[1][2]*z0;
    vrm[2][3] = z0 - m[2][0]*x0 - m[2][1]*y0 - m[2][2]*z0;

  }


  void  VectorRotation ( PPAtom  A,
                         int      nA,
                         realtype alpha,
                         realtype vx,
                         realtype vy,
                         realtype vz,
                         realtype x0,
                         realtype y0,
                         realtype z0 )  {
  //   Vector rotation is rotation by angle alpha about arbitrary
  // vector directed as (vx,vy,vz) = (vx2-vx1,vy2-vy1,vz2-vz1).
  // Point (x0,y0,z0) is the center of rotation -- actually
  // a point belonging to the rotation axis.
  mat33    m;
  realtype x,y,z;
  int      i;

    GetVecRotMatrix ( m, alpha,vx,vy,vz );

    for (i=0;i<nA;i++)
      if (A[i])  {
        if ((!A[i]->Ter) && (A[i]->WhatIsSet & ASET_Coordinates))  {
          x = A[i]->x - x0;
          y = A[i]->y - y0;
          z = A[i]->z - z0;
          A[i]->x = m[0][0]*x + m[0][1]*y + m[0][2]*z + x0;
          A[i]->y = m[1][0]*x + m[1][1]*y + m[1][2]*z + y0;
          A[i]->z = m[2][0]*x + m[2][1]*y + m[2][2]*z + z0;
        }
      }

  }


  void  GetMassCenter ( PPAtom    A,   int        nA,
                        realtype & xmc, realtype & ymc,
                        realtype & zmc )  {
  realtype w,mass;
  int      i,k;

    xmc  = 0.0;
    ymc  = 0.0;
    zmc  = 0.0;
    mass = 0.0;

    for (i=0;i<nA;i++)
      if (A[i])  {
        if ((!A[i]->Ter) && (A[i]->WhatIsSet & ASET_Coordinates))  {
          k = getElementNo ( A[i]->element );
          if (k>=0)  w = MolecWeight[k];
               else  w = 1.0;
          xmc  += w*A[i]->x;
          ymc  += w*A[i]->y;
          zmc  += w*A[i]->z;
          mass += w;
        }
      }

    if (mass>0.0)  {
      xmc /= mass;
      ymc /= mass;
      zmc /= mass;
    }

  }

  int CoorManager::BringToUnitCell()  {
  // brings all chains into 0th unit cell
  PChain    chain;
  PPAtom    atom;
  realtype  x0,y0,z0, x,y,z, xf,yf,zf, sx,sy,sz;
  realtype  dx,dy,dz, d,d0;
  int       nAtoms;
  int       i,j,k,n,m,nt, ic,jc,kc, is,js,ks;

    if (!cryst.areMatrices())  return -1;

    is = 0;  // this is only
    js = 0;  //   to depress
    ks = 0;  //      "uninitialized" worning

    cryst.Frac2Orth ( 0.5,0.5,0.5, x0,y0,z0 );

    nt = 0;
    for (i=0;i<nModels;i++)
      if (model[i])  {
        for (j=0;j<model[i]->nChains;j++)  {
          chain = model[i]->chain[j];
          if (chain)  {

            x = 0.0;
            y = 0.0;
            z = 0.0;
            m = 0;
            for (k=0;k<chain->nResidues;k++)
              if (chain->residue[k])  {
                chain->residue[k]->GetAtomTable ( atom,nAtoms );
                for (n=0;n<nAtoms;n++)
                  if (atom[n])  {
                    if (!atom[n]->Ter)  {
                      x += atom[n]->x;
                      y += atom[n]->y;
                      z += atom[n]->z;
                      m++;
                    }
                  }
              }
            x /= m;
            y /= m;
            z /= m;

            cryst.Orth2Frac ( x,y,z, xf,yf,zf );
            sx = frac ( xf );
            sy = frac ( yf );
            sz = frac ( zf );
            d0 = MaxReal;
            for (ic=-3;ic<3;ic++)
              for (jc=-3;jc<3;jc++)
                for (kc=-3;kc<3;kc++)  {
                  cryst.Frac2Orth ( sx+ic,sy+jc,sz+kc, dx,dy,dz );
                  dx -= x0;
                  dy -= y0;
                  dz -= z0;
                  d = dx*dx + dy*dy + dz*dz;
                  if (d<d0)  {
                    d0 = d;
                    is = ic;
                    js = jc;
                    ks = kc;
                  }
                }

            sx = xf - (sx+is);
            sy = yf - (sy+js);
            sz = zf - (sz+ks);

            if ((fabs(sx)>1.0e-10) || (fabs(sy)>1.0e-10)
                                   || (fabs(sz)>1.0e-10))  {
              nt++;
              for (k=0;k<chain->nResidues;k++)
                if (chain->residue[k])  {
                  chain->residue[k]->GetAtomTable ( atom,nAtoms );
                  for (n=0;n<nAtoms;n++)
                    if (atom[n])  {
                      if (!atom[n]->Ter)  {
                        cryst.Orth2Frac ( atom[n]->x,atom[n]->y,
                                          atom[n]->z,
                                          xf,yf,zf );
                        cryst.Frac2Orth ( xf-sx,yf-sy,zf-sz,
                                          atom[n]->x,atom[n]->y,
                                          atom[n]->z );
                      }
                    }
                }
            }

          }
        }
      }

    return nt;  // number of converted chains

  }


  bool CoorManager::Frac2Orth (
                realtype   xfrac, realtype   yfrac, realtype   zfrac,
                realtype & xorth, realtype & yorth, realtype & zorth )  {
    return cryst.Frac2Orth ( xfrac,yfrac,zfrac, xorth,yorth,zorth );
  }

  bool CoorManager::Orth2Frac (
                realtype   xorth, realtype   yorth, realtype   zorth,
                realtype & xfrac, realtype & yfrac, realtype & zfrac )  {
    return cryst.Orth2Frac ( xorth,yorth,zorth, xfrac,yfrac,zfrac );
  }


  bool CoorManager::Frac2Orth ( mat44 & F, mat44 & T )  {
    return cryst.Frac2Orth ( F,T );
  }

  bool CoorManager::Orth2Frac ( mat44 & T, mat44 & F )  {
    return cryst.Orth2Frac ( T,F );
  }



  //  ------------------------  Contacts  -------------------------------


  #define  CA_CA_Dist2  16.0

  void  CoorManager::FindSeqSection ( PAtom atom, int seqDist,
                                           int  & seq1, int & seq2 )  {
  PAtom    a;
  PResidue res;
  PChain   chain;
  realtype  x0,y0,z0, x,y,z, dx,dy,dz, d2;
  int       i1;
  bool   B0,B;

    x  = 0.0;
    y  = 0.0;
    z  = 0.0;
    x0 = 0.0;
    y0 = 0.0;
    z0 = 0.0;

    res = atom->residue;
    if ((!res) || (seqDist<=0))  {
      seq1 = MaxInt4;
      seq2 = MinInt4;
      return;
    }

    chain = res->chain;
    if (!chain)  {
      seq1 = MaxInt4;
      seq2 = MinInt4;
      return;
    }

    if (seqDist==1)  {
      seq1 = res->index;
      seq2 = seq1;
      return;
    }

    a = res->GetAtom ( pstr("CA"),pstr("C"),NULL );
    if (a)  {
      x0 = a->x;
      y0 = a->y;
      z0 = a->z;
      B0 = true;
    } else
      B0 = false;
    if (B0)  {
      x = x0;
      y = y0;
      z = z0;
    }

    B    = B0;
    seq2 = res->index;
    i1   = IMin(chain->nResidues,seq2+seqDist)-1;
    while (seq2<i1)  {
      seq2++;
      if (chain->residue[seq2])  {
        a = chain->residue[seq2]->GetAtom ( pstr("CA"),pstr("C"),NULL );
        if (a && B)  {
          dx = x-a->x;
          dy = y-a->y;
          dz = z-a->z;
          d2 = dx*dx + dy*dy + dz*dz;
          if (d2>CA_CA_Dist2)  {
            seq2--;
            break;
          }
        }
        if (a)  {
          x = a->x;
          y = a->y;
          z = a->z;
          B = true;
        } else
          B = false;
      }
    }

    if (B0)  {
      x = x0;
      y = y0;
      z = z0;
    }
    B    = B0;
    seq1 = res->index;
    i1   = IMax(0,seq1-seqDist+1);
    while (seq1>i1)  {
      seq1--;
      if (chain->residue[seq1])  {
        a = chain->residue[seq1]->GetAtom ( pstr("CA"),pstr("C"),NULL );
        if (a && B)  {
          dx = x-a->x;
          dy = y-a->y;
          dz = z-a->z;
          d2 = dx*dx + dy*dy + dz*dz;
          if (d2>CA_CA_Dist2)  {
            seq1++;
            break;
          }
        }
        if (a)  {
          x = a->x;
          y = a->y;
          z = a->z;
          B = true;
        } else
          B = false;
      }
    }

  }


  bool  CoorManager::iContact ( PAtom    a1, PAtom    a2,
                                         int     seq1, int     seq2,
                                         realtype  dd, realtype d12,
                                         realtype d22, realtype & d2 )  {
  //  seq1..seq2 is forbidden region for residue sequence numbers
  PResidue  res1,res2;
  PChain    chain1,chain2;
  realtype   dx,dy,dz;

    if (a2->Ter)  return false;

    dx = fabs(a2->x-a1->x);
    if (dx<=dd)  {
      dy = fabs(a2->y-a1->y);
      if (dy<=dd)  {
        dz = fabs(a2->z-a1->z);
        if (dz<=dd)  {
          d2 = dx*dx + dy*dy + dz*dz;
          if ((d12<=d2) && (d2<=d22))  {
            if (seq1<=seq2)  {
              res1 = a1->residue;
              res2 = a2->residue;
              if (res1 && res2)  {
                chain1 = res1->chain;
                chain2 = res2->chain;
                if (chain1 && chain2)  {
                  if (!strcmp(chain1->chainID,chain2->chainID))  {
                    if ((seq1<=res2->index) && (res2->index<=seq2))
                      return false;
                  }
                }
              }
            }
            return true;
          }
        }
      }
    }

    return false;

  }

  bool  CoorManager::iContact ( realtype   x, realtype   y,
                                         realtype   z, PAtom    a2,
                                         realtype  dd, realtype d12,
                                         realtype d22, realtype & d2 )  {
  realtype dx,dy,dz;

    if (a2->Ter)  return false;

    dx = fabs(a2->x-x);
    if (dx<=dd)  {
      dy = fabs(a2->y-y);
      if (dy<=dd)  {
        dz = fabs(a2->z-z);
        if (dz<=dd)  {
          d2 = dx*dx + dy*dy + dz*dz;
          if ((d12<=d2) && (d2<=d22))  return true;
        }
      }
    }

    return false;

  }


  void  CoorManager::SeekContacts ( PPAtom    AIndex,
                                    int       ilen,
                                    int       atomNum,
                                    realtype  dist1,
                                    realtype  dist2,
                                    int       seqDist,
                                    RPContact contact,
                                    int &     ncontacts,
                                    int       maxlen,
                                    long      group )  {
  PContactIndex contactIndex;
  realtype      d12,d22,d2;
  int           i,seq1,seq2;

    if (!AIndex)              return;
    if (dist2<dist1)          return;
    if (!AIndex[atomNum])     return;
    if (AIndex[atomNum]->Ter) return;

    contactIndex = new ContactIndex ( contact,maxlen,ncontacts,ilen );

    FindSeqSection ( AIndex[atomNum],seqDist,seq1,seq2 );

    d12 = dist1*dist1;
    d22 = dist2*dist2;

    for (i=0;i<ilen;i++)
      if ((i!=atomNum) && AIndex[i])  {
        if (iContact(AIndex[atomNum],AIndex[i],seq1,seq2,dist2,
                      d12,d22,d2))
           contactIndex->AddContact ( atomNum,i,sqrt(d2),group );
      }

    contactIndex->GetIndex ( contact,ncontacts );

    delete contactIndex;

  }


  void  CoorManager::SeekContacts ( PAtom     A,
                                    PPAtom    AIndex,
                                    int       ilen,
                                    realtype  dist1,
                                    realtype  dist2,
                                    int       seqDist,
                                    RPContact contact,
                                    int &     ncontacts,
                                    int       maxlen,
                                    long      group
                                   )  {
  PContactIndex contactIndex;
  realtype      d12,d22,d2;
  int           i,seq1,seq2;

    if (!AIndex)     return;
    if (dist2<dist1) return;
    if (!A)          return;
    if (A->Ter)      return;

    contactIndex = new ContactIndex ( contact,maxlen,ncontacts,ilen );

    FindSeqSection ( A,seqDist,seq1,seq2 );

    d12 = dist1*dist1;
    d22 = dist2*dist2;

    for (i=0;i<ilen;i++)
      if ((AIndex[i]!=A) && AIndex[i])  {
        if (iContact(A,AIndex[i],seq1,seq2,dist2,d12,d22,d2))
          contactIndex->AddContact ( -1,i,sqrt(d2),group );
      }

    contactIndex->GetIndex ( contact,ncontacts );

    delete contactIndex;

  }


  void  CoorManager::SeekContacts ( PPAtom    AIndex1,
                                    int       ilen1,
                                    PPAtom    AIndex2,
                                    int       ilen2,
                                    realtype  dist1,
                                    realtype  dist2,
                                    int       seqDist,
                                    RPContact contact,
                                    int &     ncontacts,
                                    int       maxlen,
                                    mat44 *   TMatrix,
                                    long      group,
                                    int       bricking,
                                    bool      doSqrt
                                  )  {
  //  It is Ok to have NULL pointers in AIndex1 and AIndex2
  PContactIndex contactIndex;
  PPAtom        A1,A2;
  rvector       sx0,sy0,sz0;
  rvector       dx0,dy0,dz0;
  realtype      d12,d22,d2, eps;
  int           l1,l2, i,j, nx,ny,nz, dn;
  int           ix1,ix2, iy1,iy2, iz1,iz2, ix,iy,iz;
  int           seq1,seq2;
  PBrick        B;
  bool          swap,UnitT;

    if ((dist2<dist1) || (!AIndex1) || (!AIndex2))  return;

    contactIndex = new ContactIndex ( contact,maxlen,ncontacts,
                                      ilen1*ilen2 );

    sx0   = NULL;
    sy0   = NULL;
    sz0   = NULL;
    dx0   = NULL;
    dy0   = NULL;
    dz0   = NULL;
    UnitT = true;
    if (TMatrix)  {
      // Transformation matrix is given. Check that that is not
      // the unit one.
      eps = 1.0e-6;
      for (i=0;(i<3) && UnitT;i++)
        for (j=0;(j<4) && UnitT;j++)
          if (i==j)  UnitT = fabs(1.0-(*TMatrix)[i][j])<eps;
               else  UnitT = fabs((*TMatrix)[i][j])<eps;
      if (!UnitT)  {
        // A non-unit transformation to AIndex2 is required.
        // As AIndex1 and AIndex2 may overlap, we have to save
        // the original AIndex1 coordinates
        GetVectorMemory ( sx0,ilen1,0 );
        GetVectorMemory ( sy0,ilen1,0 );
        GetVectorMemory ( sz0,ilen1,0 );
        for (i=0;i<ilen1;i++)
          if (AIndex1[i])  {
            sx0[i] = AIndex1[i]->x;
            sy0[i] = AIndex1[i]->y;
            sz0[i] = AIndex1[i]->z;
          }
        // Save original AIndex2 coordinates and modify the index
        GetVectorMemory ( dx0,ilen2,0 );
        GetVectorMemory ( dy0,ilen2,0 );
        GetVectorMemory ( dz0,ilen2,0 );
        for (i=0;i<ilen2;i++)
          if (AIndex2[i])  {
            dx0[i] = AIndex2[i]->x;
            dy0[i] = AIndex2[i]->y;
            dz0[i] = AIndex2[i]->z;
            AIndex2[i]->Transform ( *TMatrix );
          }
      }
    }

    // choose A2 as the largest atom set convenient for
    // bricking (bricking on larger set is more efficient)
    if (bricking & BRICK_ON_1)  {
      A1   = AIndex2;
      A2   = AIndex1;
      l1   = ilen2;
      l2   = ilen1;
      swap = true;
    } else if (bricking & BRICK_ON_2)  {
      A1   = AIndex1;
      A2   = AIndex2;
      l1   = ilen1;
      l2   = ilen2;
      swap = false;
    } else if (ilen1<=ilen2)  {
      A1   = AIndex1;
      A2   = AIndex2;
      l1   = ilen1;
      l2   = ilen2;
      swap = false;
    } else  {
      A1   = AIndex2;
      A2   = AIndex1;
      l1   = ilen2;
      l2   = ilen1;
      swap = true;
    }

    d12 = dist1*dist1;
    d22 = dist2*dist2;

    if (((bricking & BRICK_READY)==0) || (!brick))
      MakeBricks ( A2,l2,dist2*1.5 );

    dn = mround(dist2/brick_size)+1;

    if (brick)
      for (i=0;i<l1;i++)
        if (A1[i])  {
          if (!A1[i]->Ter)  {
            if (UnitT)  {
              // No transformation -- AIndex1 and AIndex2 are unmodified.
              // Calculate the forbidden sequence region
              FindSeqSection ( A1[i],seqDist,seq1,seq2 );
              // And the brick location
              GetBrickCoor   ( A1[i],nx,ny,nz );
            } else  {
              // AIndex2 and AIndex1 are modified, but the sequence
              // distance does not apply to physically different chains
              // (meaning that transformation of A2 effectively makes
              // a different chain). Use unmodified atom coordinates
              // of 1st set for calculating the brick location.
              if (swap) GetBrickCoor ( A1[i],nx,ny,nz ); // A1 is AIndex2
                   else GetBrickCoor ( sx0[i],sy0[i],sz0[i],nx,ny,nz );
          }
          if (nx>=0)  {
            ix1 = IMax ( 0,nx-dn );
            iy1 = IMax ( 0,ny-dn );
            iz1 = IMax ( 0,nz-dn );
            ix2 = IMin ( nbrick_x,nx+dn+1 );
            iy2 = IMin ( nbrick_y,ny+dn+1 );
            iz2 = IMin ( nbrick_z,nz+dn+1 );
            if (UnitT)  {
              // AIndex1 unmodified, use it
              for (ix=ix1;ix<ix2;ix++)
                if (brick[ix])
                  for (iy=iy1;iy<iy2;iy++)
                    if (brick[ix][iy])
                      for (iz=iz1;iz<iz2;iz++)  {
                        B = brick[ix][iy][iz];
                        if (B)
                          for (j=0;j<B->nAtoms;j++)
                            if (B->atom[j]!=A1[i])  {
                              if (iContact(A1[i],B->atom[j],seq1,seq2,
                                            dist2,d12,d22,d2))  {
                                if (doSqrt)  d2 = sqrt(d2);
                                if (swap)  contactIndex->AddContact (
                                             B->id[j],i,d2,group );
                                     else  contactIndex->AddContact (
                                             i,B->id[j],d2,group );
                              }
                            }
                      }
            } else if (swap)  {
              // A1 stands for AIndex2, it is modified and we need to use
              // the modified coordinates
              for (ix=ix1;ix<ix2;ix++)
                if (brick[ix])
                  for (iy=iy1;iy<iy2;iy++)
                    if (brick[ix][iy])
                      for (iz=iz1;iz<iz2;iz++)  {
                        B = brick[ix][iy][iz];
                        if (B)
                          for (j=0;j<B->nAtoms;j++)
                            if (iContact(A1[i]->x,A1[i]->y,A1[i]->z,
                                          B->atom[j], dist2,d12,d22,d2))  {
                              if (doSqrt)  d2 = sqrt(d2);
                              contactIndex->AddContact ( B->id[j],i,d2,group );
                            }
                      }
            } else  {
              // A1 stands for AIndex1, it may be modified (if AIndex1
              // and AIndex2 overlap) -- use its unmodified coordinates
              // instead.
              for (ix=ix1;ix<ix2;ix++)
                if (brick[ix])
                  for (iy=iy1;iy<iy2;iy++)
                    if (brick[ix][iy])
                      for (iz=iz1;iz<iz2;iz++)  {
                        B = brick[ix][iy][iz];
                        if (B)
                          for (j=0;j<B->nAtoms;j++)
                            if (iContact(sx0[i],sy0[i],sz0[i],
                                          B->atom[j],dist2,d12,d22,d2))  {
                              if (doSqrt)  d2 = sqrt(d2);
                              contactIndex->AddContact ( i,B->id[j],d2,group );
                            }
                      }
            }
          }
        }
      }


    if (!UnitT)  {
      // restore original coordinates
      for (i=0;i<ilen1;i++)
        if (AIndex1[i])  {
          AIndex1[i]->x = sx0[i];
          AIndex1[i]->y = sy0[i];
          AIndex1[i]->z = sz0[i];
        }
      for (i=0;i<ilen2;i++)
        if (AIndex2[i])  {
          AIndex2[i]->x = dx0[i];
          AIndex2[i]->y = dy0[i];
          AIndex2[i]->z = dz0[i];
        }
      FreeVectorMemory ( sx0,0 );
      FreeVectorMemory ( sy0,0 );
      FreeVectorMemory ( sz0,0 );
      FreeVectorMemory ( dx0,0 );
      FreeVectorMemory ( dy0,0 );
      FreeVectorMemory ( dz0,0 );
    }

    contactIndex->GetIndex ( contact,ncontacts );

    delete contactIndex;

  }


  void  CoorManager::SeekContacts ( PPAtom   AIndex1,
                                    int      ilen1,
                                    PPAtom   AIndex2,
                                    int      ilen2,
                                    realtype contDist,
                                    PContact contact,
                                    int &    ncontacts,
                                    int      bricking
                                   )  {
  //  Simplified optimized for speed version:
  //    - no NULL pointers and Ters in AIndex1 and AIndex2
  //    - no checks for identity atoms in AIndex1 and AIndex2
  //    - contact must be pre-allocated with at least ilen1*ilen2 elements
  //    - contact returns square distances
  //    - ncontacts is always reset
  PPAtom   A1,A2;
  realtype contDist2, dx,dy,dz, d2;
  int      l1,l2, i,j, nx,ny,nz, dn;
  int      ix1,ix2, iy1,iy2, iz1,iz2, ix,iy,iz;
  PBrick   B;
  bool     swap;

    // choose A2 as the largest atom set convenient for
    // bricking (bricking on larger set is more efficient)
    if (bricking & BRICK_ON_1)  {
      A1   = AIndex2;
      A2   = AIndex1;
      l1   = ilen2;
      l2   = ilen1;
      swap = true;
    } else if (bricking & BRICK_ON_2)  {
      A1   = AIndex1;
      A2   = AIndex2;
      l1   = ilen1;
      l2   = ilen2;
      swap = false;
    } else if (ilen1<=ilen2)  {
      A1   = AIndex1;
      A2   = AIndex2;
      l1   = ilen1;
      l2   = ilen2;
      swap = false;
    } else  {
      A1   = AIndex2;
      A2   = AIndex1;
      l1   = ilen2;
      l2   = ilen1;
      swap = true;
    }

    contDist2 = contDist*contDist;

    if (((bricking & BRICK_READY)==0) || (!brick))
      MakeBricks ( A2,l2,contDist*1.5 );

    ncontacts = 0;

    if (!brick)  return;

    dn = (int)floor(contDist/brick_size)+1;

    if (swap)  {

      for (i=0;i<l1;i++)
        if (A1[i])  {
          // Find brick location
          GetBrickCoor ( A1[i],nx,ny,nz );
          if (nx>=0)  {
            ix1 = IMax ( 0,nx-dn );
            iy1 = IMax ( 0,ny-dn );
            iz1 = IMax ( 0,nz-dn );
            ix2 = IMin ( nbrick_x,nx+dn+1 );
            iy2 = IMin ( nbrick_y,ny+dn+1 );
            iz2 = IMin ( nbrick_z,nz+dn+1 );
            for (ix=ix1;ix<ix2;ix++)
              if (brick[ix])
                for (iy=iy1;iy<iy2;iy++)
                  if (brick[ix][iy])
                    for (iz=iz1;iz<iz2;iz++)  {
                      B = brick[ix][iy][iz];
                      if (B)
                        for (j=0;j<B->nAtoms;j++)  {
                          dx = A1[i]->x - B->atom[j]->x;
                          dy = A1[i]->y - B->atom[j]->y;
                          dz = A1[i]->z - B->atom[j]->z;
                          d2 = dx*dx + dy*dy + dz*dz;
                          if (d2<=contDist2)  {
                            contact[ncontacts].id1  = B->id[j];
                            contact[ncontacts].id2  = i;
                            contact[ncontacts].dist = d2;
                            ncontacts++;
                          }
                        }
                    }
          }
        }

    } else  {

      for (i=0;i<l1;i++)
        if (A1[i])  {
          // Find brick location
          GetBrickCoor ( A1[i],nx,ny,nz );
          if (nx>=0)  {
            ix1 = IMax ( 0,nx-dn );
            iy1 = IMax ( 0,ny-dn );
            iz1 = IMax ( 0,nz-dn );
            ix2 = IMin ( nbrick_x,nx+dn+1 );
            iy2 = IMin ( nbrick_y,ny+dn+1 );
            iz2 = IMin ( nbrick_z,nz+dn+1 );
            for (ix=ix1;ix<ix2;ix++)
              if (brick[ix])
                for (iy=iy1;iy<iy2;iy++)
                  if (brick[ix][iy])
                    for (iz=iz1;iz<iz2;iz++)  {
                      B = brick[ix][iy][iz];
                      if (B)
                        for (j=0;j<B->nAtoms;j++)  {
                          dx = A1[i]->x - B->atom[j]->x;
                          dy = A1[i]->y - B->atom[j]->y;
                          dz = A1[i]->z - B->atom[j]->z;
                          d2 = dx*dx + dy*dy + dz*dz;
                          if (d2<=contDist2)  {
                            contact[ncontacts].id1  = i;
                            contact[ncontacts].id2  = B->id[j];
                            contact[ncontacts].dist = d2;
                            ncontacts++;
                          }
                        }
                    }
          }
        }

    }

  }

  
  void  CoorManager::SeekContacts ( vect3  * xyz,
                                    int      nxyz,
                                    realtype contDist,
                                    PContact contact,
                                    int &    ncontacts
                                   )  {
  //  Simplified optimized for speed and convenience version:
  //    - bricking is pre-done
  //    - contacting set of atoms is given as a bare vect3 (xyz)
  //      coordinate vector
  //    - no checks for identity atoms
  //    - contact must be pre-allocated with at least ilen1*ilen2
  //      elements
  //    - contact returns square distances
  //    - ncontacts is always reset
  realtype contDist2, dx,dy,dz, d2;
  int      i,j, nx,ny,nz, dn;
  int      ix1,ix2, iy1,iy2, iz1,iz2, ix,iy,iz;
  PBrick   B;

    contDist2 = contDist*contDist;

    ncontacts = 0;

    if (!brick)  return;
    
    dn = (int)floor(contDist/brick_size)+1;

    for (i=0;i<nxyz;i++)  {
      // Find brick location
      GetBrickCoor ( xyz[i],nx,ny,nz );
      if (nx>=0)  {
        ix1 = IMax ( 0,nx-dn );
        iy1 = IMax ( 0,ny-dn );
        iz1 = IMax ( 0,nz-dn );
        ix2 = IMin ( nbrick_x,nx+dn+1 );
        iy2 = IMin ( nbrick_y,ny+dn+1 );
        iz2 = IMin ( nbrick_z,nz+dn+1 );
        for (ix=ix1;ix<ix2;ix++)
          if (brick[ix])
            for (iy=iy1;iy<iy2;iy++)
              if (brick[ix][iy])
                for (iz=iz1;iz<iz2;iz++)  {
                  B = brick[ix][iy][iz];
                  if (B)
                    for (j=0;j<B->nAtoms;j++)  {
                      dx = xyz[i][0] - B->atom[j]->x;
                      dy = xyz[i][1] - B->atom[j]->y;
                      dz = xyz[i][2] - B->atom[j]->z;
                      d2 = dx*dx + dy*dy + dz*dz;
                      if (d2<=contDist2)  {
                        contact[ncontacts].id1  = B->id[j];
                        contact[ncontacts].id2  = i;
                        contact[ncontacts].dist = d2;
                        ncontacts++;
                      }
                  }
                }
      }    
    }

  }
  

  void  CoorManager::SeekContacts ( PPAtom       AIndex1,
                                    int          ilen1,
                                    PPAtom *     AIndex2,
                                    ivector      ilen2,
                                    int          nStructures,
                                    realtype     dist1,
                                    realtype     dist2,
                                    PPMContact & contact,
                                    int          bricking
                                   )  {

  //  It is Ok to have NULL pointers in AIndex1 and AIndex2
  PMBrick  B;
  PAtom    A;
  realtype d12,d22,d2;
  int      dn, i,j,k, nx,ny,nz, ix1,iy1,iz1, ix2,iy2,iz2;
  int      ix,iy,iz;

    if (dist2<dist1)              return;
    if ((!AIndex1) || (!AIndex2)) return;

    d12 = dist1*dist1;
    d22 = dist2*dist2;

    if (((bricking & BRICK_READY)==0) || (!mbrick))
      MakeMBricks ( AIndex2,ilen2,nStructures,dist2*1.5 );

    contact = new PMContact[ilen1];

    dn = mround(dist2/brick_size)+1;

    if (mbrick)
      for (i=0;i<ilen1;i++)  {
        A = AIndex1[i];
        contact[i] = NULL;
        if (A)  {
          if (!A->Ter)  {
            contact[i] = new MContact(nStructures);
            contact[i]->contactID = i;
            //  Calculate the brick location
            GetMBrickCoor ( A,nx,ny,nz );
            if (nx>=0)  {
              ix1 = IMax ( 0,nx-dn );
              iy1 = IMax ( 0,ny-dn );
              iz1 = IMax ( 0,nz-dn );
              ix2 = IMin ( nmbrick_x,nx+dn+1 );
              iy2 = IMin ( nmbrick_y,ny+dn+1 );
              iz2 = IMin ( nmbrick_z,nz+dn+1 );
              for (ix=ix1;ix<ix2;ix++)
                if (mbrick[ix])
                  for (iy=iy1;iy<iy2;iy++)
                    if (mbrick[ix][iy])
                      for (iz=iz1;iz<iz2;iz++)  {
                        B = mbrick[ix][iy][iz];
                        if (B)
                          for (j=0;j<nStructures;j++)
                            for (k=0;k<B->nAtoms[j];k++)
                              if (B->atom[j][k]!=A)  {
                                if (iContact(A,B->atom[j][k],
                                              MaxInt4,MinInt4,
                                              dist2,d12,d22,d2))
                                  contact[i]->AddContact (
                                         B->atom[j][k],j,B->id[j][k] );
                              }
                      }
            }
          }
        }
      }
    else
      for (i=0;i<ilen1;i++)
        contact[i] = NULL;

  }



  DefineClass(QSortContacts)

  class QSortContacts : public QuickSort  {
    public :
      QSortContacts() : QuickSort() {}
      int  Compare ( int i, int j );
      void Swap    ( int i, int j );
      void Sort    ( PContact contact, int ncontacts, int sortmode );
    protected :
      int  mode;
  };

  int QSortContacts::Compare ( int i, int j )  {
  bool gt,lt;
    switch (mode)  {
      default          :
      case CNSORT_1INC : gt = (((PContact)data)[i].id1 >
                               ((PContact)data)[j].id1);
                         lt = (((PContact)data)[i].id1 <
                               ((PContact)data)[j].id1);
                     break;
      case CNSORT_1DEC : gt = (((PContact)data)[j].id1 >
                               ((PContact)data)[i].id1);
                         lt = (((PContact)data)[j].id1 <
                               ((PContact)data)[i].id1);
                     break;
      case CNSORT_2INC : gt = (((PContact)data)[i].id2 >
                               ((PContact)data)[j].id2);
                         lt = (((PContact)data)[i].id2 <
                               ((PContact)data)[j].id2);
                     break;
      case CNSORT_2DEC : gt = (((PContact)data)[j].id2 >
                               ((PContact)data)[i].id2);
                         lt = (((PContact)data)[j].id2 <
                               ((PContact)data)[i].id2);
                     break;
      case CNSORT_DINC : gt = (((PContact)data)[i].dist >
                               ((PContact)data)[j].dist);
                         lt = (((PContact)data)[i].dist <
                               ((PContact)data)[j].dist);
                     break;
      case CNSORT_DDEC : gt = (((PContact)data)[j].dist >
                               ((PContact)data)[i].dist);
                         lt = (((PContact)data)[j].dist <
                               ((PContact)data)[i].dist);
                     break;
    }
    if (gt)  return  1;
    if (lt)  return -1;
    return 0;
  }

  void QSortContacts::Swap ( int i, int j )  {
    ((PContact)data)[i].Swap ( ((PContact)data)[j] );
  }


  void QSortContacts::Sort ( PContact contact, int ncontacts,
                             int sortmode )  {
    mode = sortmode;
    if (mode!=CNSORT_OFF)
      QuickSort::Sort ( &(contact[0]),ncontacts );
  }


  void  SortContacts ( PContact contact, int ncontacts,
                       CNSORT_DIR sortmode )  {
  QSortContacts SC;
    if (sortmode!=CNSORT_OFF)
      SC.Sort ( contact,ncontacts,sortmode );
  }


  //  -------------------  Stream functions  ----------------------

  void  CoorManager::write ( io::RFile f )  {
  byte Version=1;
    f.WriteByte ( &Version    );
    Root::write ( f );
    if (!isCompactBinary())  {
      f.WriteInt  ( &CoorIDCode );
      f.WriteReal ( &brick_size );
      f.WriteReal ( &xbrick_0   );
      f.WriteReal ( &ybrick_0   );
      f.WriteReal ( &zbrick_0   );
      f.WriteInt  ( &nbrick_x   );
      f.WriteInt  ( &nbrick_y   );
      f.WriteInt  ( &nbrick_z   );
    }
  }

  void  CoorManager::read ( io::RFile f )  {
  byte Version;
    f.ReadByte ( &Version    );
    Root::read ( f );
    if (!isCompactBinary())  {
      f.ReadInt  ( &CoorIDCode );
      f.ReadReal ( &brick_size );
      f.ReadReal ( &xbrick_0   );
      f.ReadReal ( &ybrick_0   );
      f.ReadReal ( &zbrick_0   );
      f.ReadInt  ( &nbrick_x   );
      f.ReadInt  ( &nbrick_y   );
      f.ReadInt  ( &nbrick_z   );
    }
  }


  MakeStreamFunctions(CoorManager);



  // ===================================================================

  int  SuperposeAtoms ( mat44 & T, PPAtom A1, int nA, PPAtom A2,
                        ivector C )  {
  realtype xc1,yc1,zc1, xc2,yc2,zc2, det,B;
  rmatrix  A,U,V;
  rvector  W,RV1;
  vect3    vc1,vc2;
  int      i,j,k,i1,i2,nat;


    //  1. Set unit matrix as "default" return

    for (i=0;i<4;i++)  {
      for (j=0;j<4;j++)
        T[i][j] = 0.0;
      T[i][i] = 1.0;
    }


    //  2. Calculate mass centers

    xc1 = 0.0;
    yc1 = 0.0;
    zc1 = 0.0;
    xc2 = 0.0;
    yc2 = 0.0;
    zc2 = 0.0;

    nat = 0;
    if (C)  {

      for (i1=0;i1<nA;i1++)
        if (!A1[i1]->Ter)  {
          i2 = C[i1];
          if (i2>=0)  {
            xc1 += A1[i1]->x;
            yc1 += A1[i1]->y;
            zc1 += A1[i1]->z;
            xc2 += A2[i2]->x;
            yc2 += A2[i2]->y;
            zc2 += A2[i2]->z;
            nat++;
          }
        }

    } else  {

      for (i=0;i<nA;i++)
        if ((!A1[i]->Ter) && (!A2[i]->Ter))  {
          xc1 += A1[i]->x;
          yc1 += A1[i]->y;
          zc1 += A1[i]->z;
          xc2 += A2[i]->x;
          yc2 += A2[i]->y;
          zc2 += A2[i]->z;
          nat++;
        }

    }

    if (nat>1)  {
      xc1 /= nat;
      yc1 /= nat;
      zc1 /= nat;
      xc2 /= nat;
      yc2 /= nat;
      zc2 /= nat;
    } else if (nat>0)  {
      T[0][3] = xc2 - xc1;
      T[1][3] = yc2 - yc1;
      T[2][3] = zc2 - zc1;
      return SPOSEAT_Ok;
    } else
      return SPOSEAT_NoAtoms;


    //  3.  Calculate the correlation matrix

    GetMatrixMemory ( A,3,3,1,1 );

    for (i=1;i<=3;i++)
      for (j=1;j<=3;j++)
        A[i][j] = 0.0;

    if (C)  {

      for (i1=0;i1<nA;i1++)
        if (!A1[i1]->Ter)  {
          i2 = C[i1];
          if (i2>=0)  {
            vc1[0] = A1[i1]->x - xc1;
            vc1[1] = A1[i1]->y - yc1;
            vc1[2] = A1[i1]->z - zc1;
            vc2[0] = A2[i2]->x - xc2;
            vc2[1] = A2[i2]->y - yc2;
            vc2[2] = A2[i2]->z - zc2;
            for (i=1;i<=3;i++)
              for (j=1;j<=3;j++)
                A[i][j] += vc1[j-1]*vc2[i-1];
          }
        }

    } else  {

      for (k=0;k<nA;k++)
        if ((!A1[k]->Ter) && (!A2[k]->Ter))  {
          vc1[0] = A1[k]->x - xc1;
          vc1[1] = A1[k]->y - yc1;
          vc1[2] = A1[k]->z - zc1;
          vc2[0] = A2[k]->x - xc2;
          vc2[1] = A2[k]->y - yc2;
          vc2[2] = A2[k]->z - zc2;
          for (i=1;i<=3;i++)
            for (j=1;j<=3;j++)
              A[i][j] += vc1[j-1]*vc2[i-1];
        }

    }


    //  4. Calculate transformation matrix (to be applied to A1)

    det = A[1][1]*A[2][2]*A[3][3] +
          A[1][2]*A[2][3]*A[3][1] +
          A[2][1]*A[3][2]*A[1][3] -
          A[1][3]*A[2][2]*A[3][1] -
          A[1][1]*A[2][3]*A[3][2] -
          A[3][3]*A[1][2]*A[2][1];

    //  4.1 SV-decompose the correlation matrix

    GetMatrixMemory ( U  ,3,3,1,1 );
    GetMatrixMemory ( V  ,3,3,1,1 );
    GetVectorMemory ( W  ,3,1 );
    GetVectorMemory ( RV1,3,1 );

    math::SVD ( 3,3,3,A,U,V,W,RV1,true,true,i );

    if (i!=0)  {
      FreeVectorMemory ( RV1,1 );
      FreeVectorMemory ( W  ,1 );
      FreeMatrixMemory ( V  ,3,1,1 );
      FreeMatrixMemory ( U  ,3,1,1 );
      FreeMatrixMemory ( A  ,3,1,1 );
      return SPOSEAT_SVD_Fail;
    }

    //  4.2 Check for parasite inversion and fix it if found

    if (det<=0.0)  {
      k = 0;
      B = MaxReal;
      for (j=1;j<=3;j++)
        if (W[j]<B)  {
          B = W[j];
          k = j;
        }
      for (j=1;j<=3;j++)
        V[j][k] = -V[j][k];
    }

    //  4.3 Calculate rotational part of T

    for (j=1;j<=3;j++)
      for (k=1;k<=3;k++)  {
        B = 0.0;
        for (i=1;i<=3;i++)
          B += U[j][i]*V[k][i];
        T[j-1][k-1] = B;
      }


    //  4.4 Add translational part to T

    T[0][3] = xc2 - T[0][0]*xc1 - T[0][1]*yc1 - T[0][2]*zc1;
    T[1][3] = yc2 - T[1][0]*xc1 - T[1][1]*yc1 - T[1][2]*zc1;
    T[2][3] = zc2 - T[2][0]*xc1 - T[2][1]*yc1 - T[2][2]*zc1;


    //  5. Release memory and quit

    FreeVectorMemory ( RV1,1 );
    FreeVectorMemory ( W  ,1 );
    FreeMatrixMemory ( V  ,3,1,1 );
    FreeMatrixMemory ( U  ,3,1,1 );
    FreeMatrixMemory ( A  ,3,1,1 );

    return SPOSEAT_Ok;

  }

  realtype getPhi ( PPAtom A )  {
  //
  //   A0    A1    A2    A3
  //   o-----o-----o-----o
  //            |
  //           Phi
  //
  //  -Pi <= Phi <= +Pi
  //
  vect3    U,W,V, a,b,c;
  realtype Wmag,S,T;

    U[0] = A[0]->x - A[1]->x;
    U[1] = A[0]->y - A[1]->y;
    U[2] = A[0]->z - A[1]->z;

    W[0] = A[2]->x - A[1]->x;
    W[1] = A[2]->y - A[1]->y;
    W[2] = A[2]->z - A[1]->z;

    V[0] = A[3]->x - A[2]->x;
    V[1] = A[3]->y - A[2]->y;
    V[2] = A[3]->z - A[2]->z;

    a[0] = U[1]*W[2] - W[1]*U[2];
    a[1] = U[2]*W[0] - W[2]*U[0];
    a[2] = U[0]*W[1] - W[0]*U[1];

    b[0] = V[1]*W[2] - W[1]*V[2];
    b[1] = V[2]*W[0] - W[2]*V[0];
    b[2] = V[0]*W[1] - W[0]*V[1];

    c[0] = a[1]*b[2] - b[1]*a[2];
    c[1] = a[2]*b[0] - b[2]*a[0];
    c[2] = a[0]*b[1] - b[0]*a[1];

    Wmag = sqrt(W[0]*W[0]+W[1]*W[1]+W[2]*W[2]);

    S    = c[0]*W[0] + c[1]*W[1] + c[2]*W[2];
    T    = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    T   *= Wmag;

    if ((S==0.0) && (T==0.0))  return NO_TORSION;
                         else  return atan2(S,T);

  }

  realtype getPsi ( PPAtom A )  {
  vect3    v1,v2;
  realtype l1,l2;

    v1[0] = A[0]->x - A[1]->x;
    v1[1] = A[0]->y - A[1]->y;
    v1[2] = A[0]->z - A[1]->z;

    v2[0] = A[2]->x - A[1]->x;
    v2[1] = A[2]->y - A[1]->y;
    v2[2] = A[2]->z - A[1]->z;

    l1 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
    if (l1==0.0)  l1 = 1.0;
    l2 = v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2];
    if (l2==0.0)  l2 = 1.0;

    return  acos((v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])/sqrt(l1*l2));

  }

  const realtype NO_TORSION = -MaxReal;


}  // namespace mmdb
