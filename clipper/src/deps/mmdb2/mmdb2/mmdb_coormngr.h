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
//  **** Module  :  mmdb_coormngr <interface>
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

#ifndef __MMDB_CoorMngr__
#define __MMDB_CoorMngr__

#include "mmdb_root.h"
#include "imex.h"
namespace mmdb  {

  // ===========================  Brick  ==============================

  //  bricking control
  enum BRICK_STATE  {
    BRICK_ON_1  = 0x00000001,
    BRICK_ON_2  = 0x00000002,
    BRICK_READY = 0x00000004
  };

  DefineClass(Brick);
  typedef  PPBrick * PPPBrick;

  class MMDB_IMEX Brick  {

    public :
      int     nAtoms;  // number of atoms hit into brick
      PPAtom  atom;  // pointers to atoms
      ivector id;  // atom ids (in present realization, these are
                       // indices of atoms from the bricked array)

      Brick ();
      ~Brick();

      void  Clear   ();
      void  AddAtom ( PAtom A, int atomid );

    protected :
      int  nAllocAtoms;
      void InitBrick();

  };


  // ===========================  MBrick  =============================

  //  Bricking multiple structures

  DefineClass(MBrick);
  typedef  PPMBrick * PPPMBrick;

  class MMDB_IMEX MBrick  {

    public :
      ivector  nAtoms;  // number of atoms in the brick
      PPAtom  *atom;    // pointers to atoms
      imatrix    id;    // atom ids (in present realization, these are
                        // indices of atoms from the bricked array)

      MBrick ( int nStructures );
      ~MBrick();

      void  Clear   ();
      void  AddAtom ( PAtom A, int structNo, int atomid );

    protected :
      ivector  nAlloAtoms;
      int      nStruct;
      void InitMBrick ( int nStructures );

  };



  //  ====================  GenSym  ========================

  DefineClass(GenSym);
  DefineStreamFunctions(GenSym);

  class MMDB_IMEX GenSym : public SymOps  {

    friend class CoorManager;

    public :

      GenSym ();
      GenSym ( io::RPStream Object );
      ~GenSym();

      void FreeMemory();

      int  AddSymOp    ( cpstr XYZOperation );
      //  the number of just added operation may be obtained as
      //  Nop = GenSym::GetNofSymOps()-1 .

      int  AddRenChain ( int Nop, const ChainID ch1, const ChainID ch2 );

      void Copy  ( PSymOps genSym );

      void write ( io::RFile f );
      void read  ( io::RFile f );

    protected :

      PChainID * chID1;   // pairs of chains to rename from chID1[n][i]
      PChainID * chID2;   // to chID2[n][i] for each operation n<Nops
      ivector    nChains; // number of chains to rename for each oper-n

      void InitGenSym();

    private :
      int        nOpAlloc;  // number of allocated operations

  };


  // =========================  Contact  =============================

  DefineStructure(Contact);

  struct MMDB_IMEX Contact  {
    int      id1,id2;
    long     group;
    realtype dist;
    void Copy ( RContact c );
    void Swap ( RContact c );
  };


  // ========================  MContact  =============================

  DefineClass(MContact);

  class MMDB_IMEX MContact : public io::Stream  {

    public :
      int       nStruct,contactID;
      ivector   nAtoms;
      PPAtom *  atom;
      imatrix   id;

      MContact ( int nStructures  );
      ~MContact();

      void    AddContact ( PAtom A, int structNo, int atomid );

    protected:
      ivector nAlloc;

  };

  extern void DeleteMContacts ( PPMContact & mcontact, int nContacts );


  // ======================  CoorManager  =========================

  DefineClass(CoorManager);
  DefineStreamFunctions(CoorManager);

  //  ----  Atom extraction return
  enum CID_RC  {
    CID_Ok        = 0,
    CID_NoModel   = 1,
    CID_NoChain   = 2,
    CID_NoResidue = 3,
    CID_NoAtom    = 4,
    CID_WrongPath = 5
  };

  //  ----  generate symmetry mates return codes
  enum GSM_RC  {
    GSM_Ok               = 0,
    GSM_NoSymOps         = 1,
    GSM_NoTransfMatrices = 2,
    GSM_NoCell           = 3
  };

  class MMDB_IMEX CoorManager : public Root  {

    public :

      int CoorIDCode; // last return from atom extraction procedure

      CoorManager ();
      CoorManager ( io::RPStream Object );
      ~CoorManager();


      //  ----------------------------------------------------------

      int  SetDefaultCoorID ( cpstr CID );


      //  ----------------  Bricking  ------------------------------

      void  RemoveBricks ();
      bool  areBricks    () { return (brick!=NULL); }
      void  MakeBricks   ( PPAtom atmvec, int avlen,
                           realtype Margin, realtype BrickSize=6.0 );
      void  GetBrickDimension (
                           int & nxmax, int & nymax, int & nzmax );
      void  GetBrickCoor ( PAtom A, int & nx, int & ny, int & nz );
      void  GetBrickCoor ( realtype x, realtype y, realtype z,
                           int & nx, int & ny, int & nz );
      void  GetBrickCoor ( vect3 & xyz, int & nx, int & ny, int & nz );
      PBrick GetBrick    ( int   nx, int   ny, int   nz );

      void  RemoveMBricks ();
      bool  areMBricks    ()  { return (mbrick!=NULL); }
      void  MakeMBricks   ( PPAtom * atmvec, ivector avlen,
                            int nStructures, realtype Margin,
                            realtype BrickSize=6.0 );
      void  GetMBrickDimension (
                           int & nxmax, int & nymax, int & nzmax );
      void  GetMBrickCoor ( PAtom A, int & nx, int & ny, int & nz );
      void  GetMBrickCoor ( realtype x, realtype y, realtype z,
                            int & nx, int & ny, int & nz );
      PMBrick GetMBrick  ( int   nx, int   ny, int   nz );

      //  ----------------  Extracting models  ---------------------

      int    GetNumberOfModels ()  { return nModels; }
      int    GetFirstModelNum ();
      PModel GetFirstDefinedModel();
      PModel GetModel ( int   modelNo  );  // 1<=modelNo<=nModels
      PModel GetModel ( cpstr CID );
      void   GetModelTable ( PPModel & modTable,
                              int & NumberOfModels );

      //  ----------------  Deleting models  -----------------------

      int  DeleteModel ( cpstr CID );
      int  DeleteModel ( int modelNo );  // 1<=modelNo<=nOfModels

      //  ----------------  Adding/Inserting models  ---------------

      int  AddModel     ( PModel mdl );
      int  InsModel     ( PModel mdl, int modelNo );
      void RotateModels ( int  modelNo1, int modelNo2, int rotdir );
      void SwapModels   ( int  modelNo1, int modelNo2 );

      //  ----------------  Extracting chains  ---------------------

      int     GetNumberOfChains ( int   modelNo  );
      int     GetNumberOfChains ( cpstr CID );
      PChain  GetChain ( int modelNo, const ChainID chainID );
      PChain  GetChain ( int modelNo, int   chainNo );
      PChain  GetChain ( cpstr CID );
      void    GetChainTable ( int modelNo, PPChain & chainTable,
                              int & NumberOfChains );
      void    GetChainTable ( cpstr CID, PPChain & chainTable,
                              int & NumberOfChains );

      //  -----------------  Deleting chains  ----------------------

      int  DeleteChain     ( int modelNo, const ChainID chID );
      int  DeleteChain     ( int modelNo, int chainNo  );
      int  DeleteAllChains ( int modelNo );
      int  DeleteAllChains ();

      //  ------------------  Adding chains  -----------------------

      int  AddChain ( int modelNo, PChain chain );

      //  ----------------  Extracting residues  -------------------

      int GetNumberOfResidues ( int modelNo, const ChainID chainID );
      int GetNumberOfResidues ( int modelNo, int   chainNo );
      int GetNumberOfResidues ( cpstr CID );
      PResidue GetResidue   ( int modelNo, const ChainID chainID,
                              int seqNo,   const InsCode insCode );
      PResidue GetResidue   ( int modelNo, int   chainNo,
                              int seqNo,   const InsCode insCode );
      PResidue GetResidue   ( int modelNo, const ChainID chainID,
                                                          int resNo );
      PResidue GetResidue   ( int modelNo, int   chainNo, int resNo );
      PResidue GetResidue   ( cpstr CID );
      int      GetResidueNo ( int modelNo, const ChainID chainID,
                              int seqNo,   const InsCode insCode );
      int      GetResidueNo ( int modelNo, int   chainNo,
                              int seqNo,   const InsCode insCode );
      void    GetResidueTable ( PPResidue & resTable,
                                int & NumberOfResidues );
      void    GetResidueTable ( int modelNo, const ChainID chainID,
                                PPResidue & resTable,
                                int & NumberOfResidues );
      void    GetResidueTable ( int modelNo, int chainNo,
                                PPResidue & resTable,
                                int & NumberOfResidues );
      void    GetResidueTable ( cpstr CID, PPResidue & resTable,
                                int & NumberOfResidues );


      //  -----------------  Deleting residues  -----------------------

      int DeleteResidue     ( int modelNo, const ChainID chainID,
                              int seqNo,   const InsCode insCode );
      int DeleteResidue     ( int modelNo, const ChainID chainID,
                                                     int resNo );
      int DeleteResidue     ( int modelNo, int   chainNo,
                              int seqNo,   const InsCode insCode );
      int DeleteResidue     ( int modelNo, int   chainNo, int resNo );
      int DeleteAllResidues ( int modelNo, const ChainID chainID );
      int DeleteAllResidues ( int modelNo, int   chainNo );
      int DeleteAllResidues ( int modelNo );
      int DeleteAllResidues ();
      int DeleteSolvent     ();

      //  -------------------  Adding residues  -----------------------

      int AddResidue ( int modelNo, const ChainID chainID,
                                                   PResidue res );
      int AddResidue ( int modelNo, int   chainNo, PResidue res );

      // --------------------  Extracting atoms  ----------------------

      int   GetNumberOfAtoms ()  { return nAtoms;  }
      int   GetNumberOfAtoms ( int modelNo, const ChainID chainID,
                               int seqNo,   const InsCode insCode );
      int   GetNumberOfAtoms ( int modelNo, int   chainNo,
                               int seqNo,   const InsCode insCode );
      int   GetNumberOfAtoms ( int modelNo, const ChainID chainID,
                                                      int resNo );
      int   GetNumberOfAtoms ( int modelNo, int   chainNo, int resNo );
      int   GetNumberOfAtoms ( cpstr CID );

      PAtom GetAtom (
                int            modelNo, // model serial number 1...
                const ChainID  chID,    // chain ID
                int            seqNo,   // residue sequence number
                const InsCode  insCode, // residue insertion code
                const AtomName aname,   // atom name
                const Element  elmnt,   // chemical element code or '*'
                const AltLoc   aloc     // alternate location indicator
                     );

      PAtom GetAtom (
                int            modelNo, // model serial number 1...
                const ChainID  chID,    // chain ID
                int            seqNo,   // residue sequence number
                const InsCode  insCode, // residue insertion code
                int            atomNo   // atom number 0..
                     );

      PAtom GetAtom (
                int            modelNo, // model serial number 1...
                const ChainID  chID,    // chain ID
                int            resNo,   // residue number 0..
                const AtomName aname,   // atom name
                const Element  elmnt,   // chemical element code or '*'
                const AltLoc   aloc     // alternate location indicator
                     );

      PAtom GetAtom (
                int            modelNo, // model serial number 1...
                const ChainID  chID,    // chain ID
                int            resNo,   // residue number 0..
                int            atomNo   // atom number 0..
                     );

      PAtom GetAtom (
                int            modelNo, // model serial number 1...
                int            chNo,    // chain number 0..
                int            seqNo,   // residue sequence number
                const InsCode  insCode, // residue insertion code
                const AtomName aname,   // atom name
                const Element  elmnt,   // chemical element code or '*'
                const AltLoc   aloc     // alternate location indicator
                     );

      PAtom GetAtom (
                int            modelNo, // model serial number 1...
                int            chNo,    // chain number 0...
                int            seqNo,   // residue sequence number
                const InsCode  insCode, // residue insertion code
                int            atomNo   // atom number 0...
                     );

      PAtom GetAtom (
                int            modelNo, // model serial number 1...
                int            chNo,    // chain number 0...
                int            resNo,   // residue number 0...
                const AtomName aname,   // atom name
                const Element  elmnt,   // chemical element code or '*'
                const AltLoc   aloc     // alternate location indicator
                     );

      PAtom GetAtom (
                int      modelNo, // model serial number 1...
                int      chNo,    // chain number 0...
                int      resNo,   // residue number 0...
                int      atomNo   // atom number 0...
                     );


      //   GetAtom(CID) returns atom answering to the following
      // CID pattern:
      //   /mdl/chn/seq(res).i/atm[elm]:a
      // where
      //   mdl   - model number (mandatory); at least model #1 is always
      //           present
      //   chn   - chain identifier ( mandatory)
      //   seq   - residue sequence number (mandatory)
      //   (res) - residue name in round brackets (may be omitted)
      //   .i    - insert code after a dot; if '.i' or 'i' is missing
      //           then residue without an insertion code is looked
      //           for
      //   atm   - atom name (mandatory)
      //   [elm] - chemical element code in square brackets; it may
      //           be omitted but could be helpful for e.g.
      //           distinguishing C_alpha and CA
      //   :a    - alternate location indicator after colon; if
      //           ':a' or 'a' is missing then an atom without
      //           alternate location indicator is looked for.
      // All spaces are ignored, all identifiers should be in capital
      // letters (comparisons are case-sensitive).
      PAtom GetAtom ( cpstr CID );


      void GetAtomTable ( PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable ( int modelNo, const ChainID chainID,
                          int seqNo,   const InsCode insCode,
                          PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable ( int modelNo, int   chainNo,
                          int seqNo,   const InsCode insCode,
                          PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable ( int modelNo, const ChainID chainID, int resNo,
                          PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable ( int modelNo, int     chainNo, int resNo,
                          PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable ( cpstr CID, PPAtom & atomTable,
                          int & NumberOfAtoms );


      //   GetAtomTable1(..) returns atom table without TER atoms and
      // without NULL atom pointers. NumberOfAtoms returns the actual
      // number of atom pointers in atomTable.
      //   atomTable is allocated within the function. If it was
      // not set to NULL before calling the function, the function will
      // attempt to deallocate it first.
      //   The application is responsible for deleting atomTable,
      // however it must not touch atom pointers, i.e. use simply
      // "delete atomTable;". Never pass atomTable from GetAtomTable(..)
      // into this function, unless you set it to NULL before doing that.
      void GetAtomTable1 ( PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable1 ( int modelNo, const ChainID chainID,
                           int seqNo,   const InsCode insCode,
                           PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable1 ( int modelNo, int   chainNo,
                           int seqNo,   const InsCode insCode,
                           PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable1 ( int modelNo, const ChainID chainID, int resNo,
                           PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable1 ( int modelNo, int     chainNo, int resNo,
                           PPAtom & atomTable, int & NumberOfAtoms );
      void GetAtomTable1 ( cpstr CID, PPAtom & atomTable,
                           int & NumberOfAtoms );


      //  --------------------  Deleting atoms  -----------------------

      int DeleteAtom ( int            modelNo,
                       const ChainID  chID,
                       int            seqNo,
                       const InsCode  insCode,
                       const AtomName aname,
                       const Element  elmnt,
                       const AltLoc   aloc );
      int DeleteAtom ( int           modelNo,
                       const ChainID chID,
                       int           seqNo,
                       const InsCode insCode,
                       int           atomNo );
      int DeleteAtom ( int            modelNo,
                       const ChainID  chID,
                       int            resNo,
                       const AtomName aname,
                       const Element  elmnt,
                       const AltLoc   aloc );
      int DeleteAtom ( int modelNo, const ChainID chID,
                       int resNo, int atomNo );
      int DeleteAtom ( int modelNo, int chNo, int seqNo,
                       const InsCode  insCode,
                       const AtomName aname,
                       const Element  elmnt,
                       const AltLoc   aloc );
      int DeleteAtom ( int modelNo, int chNo, int seqNo,
                       const InsCode insCode, int atomNo );
      int DeleteAtom ( int modelNo, int chNo, int resNo,
                       const AtomName aname,
                       const Element  elmnt,
                       const AltLoc   aloc );
      int DeleteAtom ( int modelNo, int chNo, int resNo, int atomNo );

      int DeleteAllAtoms ( int modelNo, const ChainID chID,
                           int seqNo,   const InsCode insCode );
      int DeleteAllAtoms ( int modelNo, const ChainID chID, int resNo );
      int DeleteAllAtoms ( int modelNo, const ChainID chID );
      int DeleteAllAtoms ( int modelNo, int chNo, int seqNo,
                           const InsCode insCode );
      int DeleteAllAtoms ( int modelNo, int chNo, int resNo );
      int DeleteAllAtoms ( int modelNo, int chNo );
      int DeleteAllAtoms ( int modelNo );
      int DeleteAllAtoms ();

      //  This function leaves only alternative location with maximal
      // occupancy, if those are equal or unspecified, the one with
      // "least" alternative location indicator.
      //  The function returns the number of deleted atoms and optimizes
      // the atom index.
      int DeleteAltLocs  ();


      //  ---------------------  Adding atoms  ------------------------

      int AddAtom ( int modelNo, const ChainID chID,
                    int seqNo,   const InsCode insCode, PAtom atom );
      int AddAtom ( int modelNo, const ChainID chID, int resNo,
                                                        PAtom atom );
      int AddAtom ( int modelNo, int chNo, int seqNo,
                                 const InsCode insCode, PAtom atom );
      int AddAtom ( int modelNo, int chNo, int resNo, PAtom  atom );


      // --------------------  Transformations  -----------------------

      int  GenerateSymMates ( PGenSym genSym=NULL );
                               // 1: no Sym operations,
                               // 2: no fract/orth matrices
                               // 3: no cell parameters
                               // 0: Ok

      void ApplyTransform   ( const mat44 & TMatrix ); // transforms all
                                            // coordinates by multiplying
                                            // with matrix TMatrix

      int  BringToUnitCell();  // brings all chains into 0th unit cell

      //   Frac2Orth(..) and Orth2Frac(..) transform between fractional
      // and orthogonal coordinates, if areMatrices() returns true.
      // If the transformation matrices were not set, the functions just
      // copy the coordinates.  Returns true if the transformation was
      // done; False return means that transformation matrices were not
      // calculated
      bool Frac2Orth (
                realtype   xfrac, realtype   yfrac, realtype   zfrac,
                realtype & xorth, realtype & yorth, realtype & zorth );
      bool Orth2Frac (
                realtype   xorth, realtype   yorth, realtype   zorth,
                realtype & xfrac, realtype & yfrac, realtype & zfrac );


      //   Below, F and T are transformation matrices in fractional and
      // orthogonal coordinates, respectively.
      bool Frac2Orth ( mat44 & F, mat44 & T );
      bool Orth2Frac ( mat44 & T, mat44 & F );

      // ====================  Seeking contacts  ======================

      void  SeekContacts (
               PPAtom     AIndex,    // index of atoms [0..ilen-1]
               int        ilen,      // length of index
               int        atomNum,   // number of 1st contact atom
                                     // in the index. All other atoms
                                     // are checked for contact with
                                     // 1st atom
               realtype   dist1,     // minimal contact distance
               realtype   dist2,     // maximal contact distance
               int        seqDist,   // the sequence distance to neglect.
                                     // If seqDist==0, all atoms are
                                     // checked for contact. If
                                     // seqDist==1, the atoms belonging
                                     // to the same residue as atom
                                     // AIndex[atomNum], are neglected.
                                     // If seqDist>1, all atoms belonging
                                     // to residues closer than
                                     // +/-(seqDist-1) around that of
                                     // atom AIndex[atomNum], are
                                     // neglected. If chain is broken
                                     // (has a gap) on section
                                     // [-(seqDist-1)..seqDist-1], the
                                     // section of neglection is
                                     // shortened to that gap.
               RPContact  contact,   // indices of contacting atoms
                                     // [0..ncontacts-1]. contact[i].id1
                                     // is set to atomNum and
                                     // contact[i].id2 is set to the
                                     // index of 2nd contacting atom
                                     // in vector AIndex
               int &      ncontacts, // number of contacts found. If
                                     // ncontacts>0 on input, it is
                                     // assumed that new contacts that
                                     // newly found contacts should be
                                     // appended to those already
                                     // existing
               int        maxlen=0,  // if <=0, then vector contact is
                                     // allocated dynamically. If
                                     // contact!=NULL, then it is
                                     // appended with new contacts.
                                     // The application is responsible
                                     // for deallocation of contact
                                     // after use.
                                     //   If maxlen>0 then vector contact
                                     // is prohibited of dynamical
                                     // allocation/deallocation. In this
                                     // case, not more than maxlen
                                     // contacts will be returned.
               long       group=0    // a contact group ID, which will be
                                     // simply stored in contact[i].group
                                     // fields. This ID may be useful
                                     // if contacts are obtained in
                                     // multiple calls of the function
                         );

      void  SeekContacts (
               PAtom      A,         // 1st atom in contact
               PPAtom     AIndex,    // index of atoms [0..ilen-1] to
                                     // check for contact with 1st atom
               int        ilen,      // length of index
               realtype   dist1,     // minimal contact distance
               realtype   dist2,     // maximal contact distance
               int        seqDist,   // the sequence distance to neglect.
                                     // If seqDist==0, all atoms are
                                     // checked for contact. If
                                     // seqDist==1, the atoms belonging
                                     // to the same residue as atom
                                     // A, are neglected. If seqDist>1,
                                     // all atoms belonging to residues
                                     // closer than +/-(seqDist-1) around
                                     // that of atom A, are neglected. If
                                     // chain is broken (has a gap) on
                                     // section
                                     // [-(seqDist-1)..seqDist-1], the
                                     // section of neglection is
                                     // shortened to that gap.
               RPContact  contact,   // indices of contacting atoms
                                     // [0..ncontacts-1]. contact[i].id1
                                     // is set to -1, and contact[i].id2
                                     // is set to the index of 2nd
                                     // contacting atom in vector AIndex
               int &      ncontacts, // number of contacts found. If
                                     // ncontacts>0 on input, it is
                                     // assumed that new contacts that
                                     // newly found contacts should be
                                     // appended those already existing
               int        maxlen=0,  // if <=0, then vector contact is
                                     // allocated dynamically. If
                                     // contact!=NULL, then it is
                                     // appended with new contacts.
                                     // The application is responsible
                                     // for deallocation of contact
                                     // after use.
                                     //   If maxlen>0 then vector contact
                                     // is prohibited of dynamical
                                     // allocation/deallocation. In this
                                     // case, not more than maxlen
                                     // contacts will be returned.
               long       group=0    // a contact group ID, which will be
                                     // simply stored in contact[i].group
                                     // fields. This ID may be useful
                                     // if contacts are obtained in
                                     // multiple calls of the function
                         );

      void  SeekContacts (
               PPAtom     AIndex1,   //  1st atom index [0..ilen1-1]
               int        ilen1,     //  length of 1st index
               PPAtom     AIndex2,   //  2nd atom index [0..ilen2-1] to
                                     // check for contact with 1st index
               int        ilen2,     //  length of 2nd index
               realtype   dist1,     //  minimal contact distance
               realtype   dist2,     //  maximal contact distance
               int        seqDist,   //  the sequence distance to
                                     // neglect.
                                     //  If seqDist==0, all atoms are
                                     // checked for contact.
                                     //  If seqDist==1, the atoms
                                     // belonging to the same residue
                                     // are neglected.
                                     //  If seqDist>1, all atoms
                                     // belonging to residues closer than
                                     // +/-(seqDist-1) to each other,
                                     // are neglected. If chain is broken
                                     // (has a gap) on section
                                     // [-(seqDist-1)..seqDist-1], the
                                     // section of neglection is
                                     // shortened to that gap.
               RPContact  contact,   //  indices of contacting atoms
                                     // [0..ncontacts-1]. contact[i].id1
                                     // contains number of atom from 1st
                                     // index, and contact[i].id2
                                     // contains number of atom from 2nd
                                     // index, contacting with the former
                                     // one
               int &      ncontacts, //  number of contacts found. If
                                     // ncontacts>0 on input, it is
                                     // assumed that newly found
                                     // contacts should be appended to
                                     // those already existing
               int        maxlen=0,  //  if <=0, then vector contact is
                                     // allocated dynamically. If
                                     // contact!=NULL, then it is
                                     // appended with new contacts.
                                     // The application is responsible
                                     // for deallocation of contact
                                     // after use.
                                     //   If maxlen>0 then vector contact
                                     // is prohibited of dynamical
                                     // allocation/deallocation. In this
                                     // case, not more than maxlen
                                     // contacts will be returned.
               mat44 * TMatrix=NULL, //  transformation matrix for 2nd
                                     // set of atoms (AIndex2)
               long       group=0,   //  a contact group ID, which will
                                     // be stored in contact[i].group
                                     // fields. This ID may be useful
                                     // if contacts are obtained in
                                     // multiple calls of the function
               int     bricking=0,   //  bricking control; may be a
                                     // combination of BRICK_ON_1 or
                                     // BRICK_ON_2 with BRICK_READY
               bool    doSqrt=true   // if False, then Contact contains
                                     // square distances
                         );

      //  Simplified optimized for speed version:
      //    - no NULL pointers and Ters in AIndex1 and AIndex2
      //    - no checks for identity atoms in AIndex1 and AIndex2
      //    - contact must be pre-allocated with at least ilen1*ilen2
      //      elements
      //    - contact returns square distances
      //    - ncontacts is always reset
      void  SeekContacts (
               PPAtom     AIndex1,   //  1st atom index [0..ilen1-1]
               int        ilen1,     //  length of 1st index
               PPAtom     AIndex2,   //  2nd atom index [0..ilen2-1] to
                                     // check for contact with 1st index
               int        ilen2,     //  length of 2nd index
               realtype   contDist,  //  maximal contact distance
               PContact   contact,   //  indices of contacting atoms
                                     // [0..ncontacts-1]. contact[i].id1
                                     // contains number of atom from 1st
                                     // index, and contact[i].id2
                                     // contains number of atom from 2nd
                                     // index, contacting with the former
                                     // one. Must be pre-allocated
               int &      ncontacts, //  number of contacts found
               int        bricking=0 //  bricking control; may be a
                                     // combination of BRICK_ON_1 or
                                     // BRICK_ON_2 with BRICK_READY
                         );
      
      //  Simplified optimized for speed and convenience version:
      //    - bricking is pre-done
      //    - contacting set of atoms is given as a bare vect3 (xyz)
      //      coordinate vector
      //    - no checks for identity atoms
      //    - contact must be pre-allocated with at least ilen1*ilen2
      //      elements
      //    - contact returns square distances
      //    - ncontacts is always reset
      void  SeekContacts (
               vect3    * xyz,       //  2nd atom index [0..ilen2-1] to
                                     // check for contact with 1st index
                                     // which was used for bricking 
               int        nxyz,      //  length of 2nd index
               realtype   contDist,  //  maximal contact distance
               PContact   contact,   //  indices of contacting atoms
                                     // [0..ncontacts-1]. contact[i].id1
                                     // contains number of atom from 1st
                                     // index, and contact[i].id2
                                     // contains number of atom from 2nd
                                     // index, contacting with the former
                                     // one. Must be pre-allocated
               int &      ncontacts  //  number of contacts found
                         );
      
      void  SeekContacts (
             PPAtom        AIndex1,  //  1st atom index [0..ilen1-1]
             int           ilen1,    //  length of 1st index
             PPAtom *      AIndex2,  //  indexes of atoms to be checked
                                     // for contact with each atom from
                                     // Aindex1; dimension
                                     // [0..nStructures-1][0..ilen2[i]-1]
             ivector       ilen2,    //  lengths of indexes AIndex2
             int           nStructures, //  number of indexes AIndex2
             realtype      dist1,    //  minimal contact distance
             realtype      dist2,    //  maximal contact distance
             PPMContact &  contact,  // resulting contacts, one structure
                                     // per each position in AIndex1. If
                                     // AIndex1[i] is NULL, contact[i] is
                                     // also NULL. "contact" is always
                                     // allocated, no re-use or
                                     // re-allocation is attempted.
             int            bricking=0  //  bricking control; may be
                                     // BRICK_READY if AIndex2 does not
                                     // change
                         );

    protected :

      //  bricks
      realtype    brick_size, xbrick_0,ybrick_0,zbrick_0;
      int         nbrick_x,nbrick_y,nbrick_z;
      PPPBrick  * brick;

      realtype    mbrick_size, xmbrick_0,ymbrick_0,zmbrick_0;
      int         nmbrick_x,nmbrick_y,nmbrick_z;
      PPPMBrick * mbrick;

      //  ---------------  Stream I/O  -----------------------------
      void  write ( io::RFile f );
      void  read  ( io::RFile f );

      void  InitMMDBCoorManager();

      void  ApplySymTransform ( int SymMatrixNo, PGenSym genSym=NULL );

      void  ResetManager ();

      void  FindSeqSection ( PAtom  atom, int  seqDist,
                             int  &  seq1, int  &  seq2 );
      bool  iContact  ( PAtom    a1, PAtom    a2,
                        int     seq1, int     seq2,
                        realtype  dd, realtype d12,
                        realtype d22, realtype & d2 );
      bool  iContact  ( realtype   x, realtype   y,
                        realtype   z, PAtom    a2,
                        realtype  dd, realtype d12,
                        realtype d22, realtype & d2 );

  };



  //  ===================================================================



  //   GetEulerRotMatrix(..) calculates the Euler rotation matrix
  // for rotation:
  //                   1) about z-axis by angle alpha
  //                   2) about new y-axis by angle beta
  //                   3) about new z-axis by angle gamma
  extern MMDB_IMEX void  GetEulerRotMatrix ( mat33 & erm,   realtype alpha,
                                   realtype beta, realtype gamma );

  //  GetEulerTMatrix(..) calculates the Euler rotation-translation
  // matrix for rotation:
  //                   1) about z-axis by angle alpha
  //                   2) about new y-axis by angle beta
  //                   3) about new z-axis by angle gamma
  //  Point (x0,y0,z0) is the center of rotation.
  extern MMDB_IMEX void  GetEulerTMatrix ( mat44 & erm,   realtype alpha,
                            realtype beta, realtype gamma,
                            realtype x0,   realtype y0,  realtype z0 );

  //  Euler rotation:  1) about z-axis by angle alpha
  //                   2) about new y-axis by angle beta
  //                   3) about new z-axis by angle gamma
  //  Point (x0,y0,z0) is the center of rotation.
  extern MMDB_IMEX void EulerRotation ( PPAtom A, int nA,
                           realtype alpha, realtype beta, realtype gamma,
                           realtype x0,    realtype y0,   realtype z0 );

  //   GetVecRotMatrix(..) calculates the rotation matrix for
  // rotation by angle alpha about arbitrary vector directed
  // as (vx,vy,vz) = (vx2-vx1,vy2-vy1,vz2-vz1).
  extern MMDB_IMEX void GetVecRotMatrix ( mat33 & vrm,  realtype alpha,
                             realtype vx,  realtype vy, realtype vz );


  //    Given the rotation matrix vrm, GetRotParameters(..)
  //  returns the rotation angle alpha and the normalized
  //  rotation axis vector (vx,vy,vz).
  //    The rotation angle and vector are determined up to
  //  their sign (however correlated, so that being substituted
  //  into GetVecRotMatrix(..) they yield the same rotation
  //  matrix).
  //    The function does not check for vrm to be a valid
  //  rotation matrix.
  extern MMDB_IMEX void GetRotParameters ( mat33 & vrm, realtype & alpha,
                         realtype & vx, realtype & vy, realtype & vz );


  //   GetVecTMatrix(..) calculates the rotation-translation matrix
  // for rotation by angle alpha about arbitrary vector directed as
  // (vx,vy,vz) = (vx2-vx1,vy2-vy1,vz2-vz1). Point (x0,y0,z0) is
  // the center of rotation -- actually a point belonging to the
  // rotation axis.
  extern MMDB_IMEX void GetVecTMatrix  ( mat44 & vrm, realtype alpha,
                               realtype vx, realtype vy, realtype vz,
                               realtype x0, realtype y0, realtype z0 );

  //   Vector rotation is rotation by angle alpha about arbitrary
  // vector directed as (vx,vy,vz) = (vx2-vx1,vy2-vy1,vz2-vz1).
  // Point (x0,y0,z0) is the center of rotation -- actually
  // a point belonging to the rotation axis.
  extern MMDB_IMEX void VectorRotation ( PPAtom A, int nA,  realtype alpha,
                               realtype vx, realtype vy, realtype vz,
                               realtype x0, realtype y0, realtype z0 );

  extern MMDB_IMEX void GetMassCenter  ( PPAtom A, int nA,
                        realtype & xmc, realtype & ymc, realtype & zmc );


  enum SPOSEAT_RC  {
    SPOSEAT_Ok       = 0,
    SPOSEAT_NoAtoms  = 1,
    SPOSEAT_SVD_Fail = 2
  };

  //   Given two sets of atoms, A1 and A2, SuperposeAtoms(...) calculates
  // the rotational-translational matrix T such that |T*A1 - A2| is
  // minimal in least-square terms.
  //   If vector C is not given (default), all nA atoms of set A1 are
  // considered as corresponding to nA first atoms of set A2,
  // A1[i] <-> A2[i], 0<=i<nA .
  //   If vector C is given, then the correspondence of atoms is
  // established as A1[i] <-> A2[C[i]] only for those i that C[i]>=0.
  // The default option (C==NULL) is thus identical to C[i]==i, 0<=i<nA.
  //   Upon normal completion, the procedure returns SPOSEAT_Ok.

  extern MMDB_IMEX int SuperposeAtoms ( mat44 & T, PPAtom A1, int nA, PPAtom A2,
                              ivector C=NULL );

  enum CNSORT_DIR  {
    CNSORT_OFF  = 0,
    CNSORT_1INC = 1,
    CNSORT_1DEC = 2,
    CNSORT_2INC = 3,
    CNSORT_2DEC = 4,
    CNSORT_DINC = 5,
    CNSORT_DDEC = 6
  };

  extern MMDB_IMEX void  SortContacts ( PContact contact, int ncontacts,
                              CNSORT_DIR sortmode );


  extern const realtype NO_TORSION;

  extern MMDB_IMEX realtype getPhi ( PPAtom A );  // A[0] - A[3] used
  extern MMDB_IMEX realtype getPsi ( PPAtom A );  // A[0] - A[2] used

}  // namespace mmdb

#endif

