//  $Id: mmdb_atom.h $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2015.
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
//    09.03.16   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  MMDB_Atom <interface>
//       ~~~~~~~~~
//  **** Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::Atom     ( atom class    )
//       ~~~~~~~~~  mmdb::Residue  ( residue class )
//  **** Functions: mmdb::BondAngle
//       ~~~~~~~~~~
//
//  Copyright (C) E. Krissinel 2000-2016
//
//  =================================================================
//

#ifndef __MMDB_Atom__
#define __MMDB_Atom__

#include "mmdb_io_stream.h"
#include "mmdb_uddata.h"
#include "mmdb_utils.h"
#include "mmdb_defs.h"

#include "imex.h"

namespace mmdb  {

  //  ======================  Atom  ==========================

  // constants for the WhatIsSet field
  enum ASET_FLAG  {
    ASET_Coordinates  = 0x00000001,
    ASET_Occupancy    = 0x00000002,
    ASET_tempFactor   = 0x00000004,
    ASET_CoordSigma   = 0x00000010,
    ASET_OccSigma     = 0x00000020,
    ASET_tFacSigma    = 0x00000040,
    ASET_Anis_tFac    = 0x00000100,
    ASET_Anis_tFSigma = 0x00001000,
    ASET_Charge       = 0x00000080,
    ASET_All          = 0x000FFFFF
  };

  const int ATOM_NoSeqNum = MinInt4;

  extern bool  ignoreSegID;
  extern bool  ignoreElement;
  extern bool  ignoreCharge;
  extern bool  ignoreNonCoorPDBErrors;
  extern bool  ignoreUnmatch;


  DefineStructure(AtomStat);

  struct MMDB_IMEX AtomStat  {

    public :
      int       nAtoms;          // number of atoms in statistics

      realtype  xmin,ymin,zmin;  // minimums of coordinates
      realtype  xmax,ymax,zmax;  // maximums of coordinates
      realtype  xm  ,ym  ,zm;    // mediums  of coordinates
      realtype  xm2 ,ym2 ,zm2;   // square mediums of coordinates

      realtype  occ_min,occ_max; // minimum/maximum occupancy
      realtype  occ_m  ,occ_m2;  // medium and square medium occupancy

      realtype  tFmin,tFmax;     // minimum/maximum temperature factor
      realtype  tFm  ,tFm2;      // medium and sq. med. temp. factor

      realtype  u11_min,u11_max; // minimums and
      realtype  u22_min,u22_max; //   maximums of
      realtype  u33_min,u33_max; //     anisotropic
      realtype  u12_min,u12_max; //       temperature
      realtype  u13_min,u13_max; //         factors
      realtype  u23_min,u23_max;

      realtype  u11_m,u11_m2;    // mediums and
      realtype  u22_m,u22_m2;    //   square mediums of
      realtype  u33_m,u33_m2;    //     anisotropic
      realtype  u12_m,u12_m2;    //       temperature
      realtype  u13_m,u13_m2;    //         factors
      realtype  u23_m,u23_m2;

      word      WhatIsSet;       //   mask field

      void  Init  ();
      void  Finish();

      realtype GetMaxSize();

    private :
      bool finished;

  };


  DefineStructure(AtomBondI);

  struct MMDB_IMEX AtomBondI  {
    int  index;  //!< bonded atom index
    byte order;  //!< bond order
  };


  DefineStructure(AtomBond);

  struct MMDB_IMEX AtomBond  {
    PAtom atom;  //!< bonded atom pointer
    byte  order;  //!< bond order
  };


  DefineFactoryFunctions(Atom);

  class MMDB_IMEX Atom : public UDData  {

    friend class Residue;
    friend class Model;
    friend class Root;
    friend class CoorManager;
    friend class SelManager;

    public :
      int        serNum;         //!< serial number
      AtomName   name;           //!< atom name (ALIGNED)
      AtomName   label_atom_id;  //!< assigned atom name (not aligned)
      AltLoc     altLoc; //!< alternative location indicator ("" for none)
      SegID      segID;          //!< segment identifier
      Element    element;        //!< element symbol (ALIGNED TO RIGHT)
      EnergyType energyType;     //!< energy type (without spaces)
      PResidue   residue;        //!< reference to residue
      realtype   x,y,z;          //!< orthogonal coordinates in angstroms
      realtype   occupancy;      //!< occupancy
      realtype   tempFactor;     //!< temperature factor
      realtype   charge;         //!< charge on the atom
      realtype   sigX,sigY,sigZ; //!< standard deviations of the coords
      realtype   sigOcc;         //!< standard deviation of occupancy
      realtype   sigTemp;        //!< standard deviation of temp. factor
      realtype   u11,u22,u33;    //!< anisotropic temperature
      realtype   u12,u13,u23;    ///    factors
      realtype   su11,su22,su33; //!< standard deviations of
      realtype   su12,su13,su23; ///    anisotropic temperature factors
      bool       Het;            //!< indicator of het atom
      bool       Ter;            //!< chain terminator

      word       WhatIsSet;      //!<   mask      field
                         ///  0x0001   atomic coordinates
                         ///  0x0002   occupancy
                         ///  0x0004   temperature factor
                         ///  0x0010   coordinate standard deviations
                         ///  0x0020   deviation of occupancy
                         ///  0x0040   deviation of temperature factor
                         ///  0x0100   anisotropic temperature factors
                         ///  0x1000   anis. temp. fact-s st-d deviations

      Atom ();
      Atom ( PResidue     res    );
      Atom ( io::RPStream Object );
      ~Atom();

      void  SetResidue   ( PResidue     res );
      void  PDBASCIIDump ( io::RFile    f   );
      void  MakeCIF      ( mmcif::PData CIF );

      //    AddBond(...) adds a bond to the atom, that is a pointer
      //  to the bonded atom and the bond order. nAdd_bonds allows
      //  one to minimize the memory reallocations, if number of
      //  bonds is known apriori: Atom adds space for nAdd_bonds
      //  if currently allocated space is exchausted.
      //    Return:  <=0  - error: bond_atom is already "bonded"
      //              >0  - Ok, returns current number of bonds
      int   AddBond  ( PAtom bond_atom, int bond_order,
                                        int nAdd_bonds=1 );
      int   GetNBonds();

      //    This GetBonds(..) returns pointer to the Atom's
      //  internal Bond structure, IT MUST NOT BE DISPOSED.
      void  GetBonds ( RPAtomBond atomBond, int & nAtomBonds );
      void  FreeBonds();

      //    This GetBonds(..) disposes AtomBondI, if it was not set
      //  to NULL, allocates AtomBondI[nAtomBonds] and returns its
      //  pointer. AtomBondI MUST BE DISPOSED BY APPLICATION.
      void  GetBonds ( RPAtomBondI atomBondI, int & nAtomBonds );

      //    This GetBonds(..) does not dispose or allocate AtomBondI.
      //  It is assumed that length of AtomBondI is sufficient to
      //  accomodate all bonded atoms.
      void  GetBonds ( PAtomBondI atomBondI, int & nAtomBonds,
                       int maxlength );


      //   ConvertPDBxxxxxx() gets data from the PDB ASCII xxxxxx
      // record (xxxxxx stands for ATOM, SIGATM, ANISOU, SIGUIJ,
      // TER or HETATM).
      //   These functions DO NOT check the xxxxxx keyword and
      // do not decode the chain and residue parameters! These
      // must be treated by the calling process, see
      // CMMDBFile::ReadPDBAtom().
      //   The atom reference is updated in the corresponding
      // residue.
      ERROR_CODE ConvertPDBATOM   ( int ix, cpstr S );
      ERROR_CODE ConvertPDBSIGATM ( int ix, cpstr S );
      ERROR_CODE ConvertPDBANISOU ( int ix, cpstr S );
      ERROR_CODE ConvertPDBSIGUIJ ( int ix, cpstr S );
      ERROR_CODE ConvertPDBTER    ( int ix, cpstr S );
      ERROR_CODE ConvertPDBHETATM ( int ix, cpstr S );

      ERROR_CODE GetCIF           ( int ix, mmcif::PLoop Loop,
                                     mmcif::PLoop LoopAnis );

      bool RestoreElementName();
      bool MakePDBAtomName();

      void  SetAtomName    ( int            ix,      // index
                             int            sN,      // serial number
                             const AtomName aName,   // atom name
                             const AltLoc   aLoc, // alternative location
                             const SegID    sID,     // segment ID
                             const Element  eName ); // element name

      //  This only renames the atom
      void  SetAtomName    ( const AtomName atomName );
      void  SetElementName ( const Element  elName   );
      void  SetCharge      ( cpstr          chrg     );
      void  SetCharge      ( realtype       chrg     );

      void  SetAtomIndex   ( int ix ); // don't use in your applications!

      void  MakeTer();  // converts atom into 'ter'

      void  SetCoordinates ( realtype xx,  realtype yy, realtype zz,
                             realtype occ, realtype tFac );

      int   GetModelNum       ();
      pstr  GetChainID        ();
      pstr  GetLabelAsymID    ();
      pstr  GetResName        ();
      pstr  GetLabelCompID    ();
      int   GetAASimilarity   ( const ResName resName );
      int   GetAASimilarity   ( PAtom  A );
      realtype GetAAHydropathy();
      realtype GetOccupancy   ();
      int   GetSeqNum         ();
      int   GetLabelSeqID     ();
      int   GetLabelEntityID  ();
      pstr  GetInsCode        ();
      int   GetSSEType        ();  // works only after SSE calculations
      pstr  GetAtomName       () { return name;    }
      pstr  GetElementName    () { return element; }
      pstr  GetAtomCharge     ( pstr chrg );

      //   GetChainCalphas(...) is a specialized function for quick
      // access to C-alphas of chain which includes given atom.
      // This function works faster than an equivalent implementation
      // through MMDB's selection procedures.
      //    Parameters:
      //       Calphas   - array to accept pointers on C-alpha atoms
      //                  If Calphas!=NULL, then the function will
      //                  delete and re-allocate it. When the array
      //                  is no longer needed, the application MUST
      //                  delete it:  delete[] Calphas; Deleting
      //                  Calphas does not delete atoms from MMDB.
      //       nCalphas   - integer to accept number of C-alpha atoms
      //                  and the length of Calphas array.
      //       altLoc     - alternative location indicator. By default
      //                  (""), maximum-occupancy locations are taken.
      void  GetChainCalphas ( PPAtom & Calphas, int & nCalphas,
                              cpstr altLoc = "" );

      bool isTer         () { return Ter; }
      bool isMetal       ();
      bool isSolvent     ();  // works only for atom in a residue!
      bool isInSelection ( int selHnd );
      bool isNTerminus   ();
      bool isCTerminus   ();

      void  CalAtomStatistics ( RAtomStat AS );

      realtype GetDist2 ( PAtom a );
      realtype GetDist2 ( PAtom a, mat44 & tm );  // tm applies to 'a'
      realtype GetDist2 ( PAtom a, mat33 & r, vect3 & t );// tm applies to A
      realtype GetDist2 ( realtype ax, realtype ay, realtype az );
      realtype GetDist2 ( mat44 & tm,  // applies to 'this'
                          realtype ax, realtype ay, realtype az  );
      realtype GetDist2 ( vect3 & xyz );

      // GetCosine(a1,a2) calculates cosine of angle a1-this-a2,
      // i.e. that between vectors [a1,this] and [this,a2].
      realtype GetCosine ( PAtom a1, PAtom a2 );

      PResidue GetResidue  ();
      PChain   GetChain    ();
      PModel   GetModel    ();
      int      GetResidueNo();
      void *   GetCoordHierarchy();  // PRoot

      //  GetAtomID(..) generates atom ID in the form
      //     /m/c/r(rn).i/n[e]:a
      //  where  m  - model number
      //         c  - chain ID
      //         r  - residue sequence number
      //         rn - residue name
      //         i  - insertion code
      //         n  - atom name
      //         e  - chemical element specification
      //         a  - alternate location indicator
      //  If any of the fields is undefined, it is replaced by
      //  hyphen  '-'.
      //    No checks on the sufficiency of string buffer AtomID
      //  is made.
      //    GetAtomID returns AtomID.
      pstr  GetAtomID ( pstr AtomID );

      pstr  GetAtomIDfmt ( pstr AtomID );

      // -------  checking atom ID
      // CheckID(..) returns 1 if atom is identified, and 0 otherwise.
      //   Parameters:
      //     aname   - atom name. It may or may not be aligned (as in
      //               a PDB file), only first word of the name will
      //               be taken ("CA", " CA" and " CA B" are all
      //               considered as "CA"). aname may be set to NULL
      //               or '*', then this parameter is ignored.
      //     elname  - element code. It will work only if element code
      //               is supplied (which might not be the case if
      //               the atom was created in a tricky way). elname
      //               should be used to distinguih between, e.g.
      //               "Ca" and "C_alpha"). elname may be set to NULL,
      //               or '*', then this parameter is ignored.
      //     aloc    - the alternate location code. aloc may be set to
      //               NULL or '*', then this parameter is ignored.
      //  IMPORTANT: comparison is case-sensitive.
      //  The atom is considered as identified, if all non-NULL
      //  parameters do match. If all parameters are set NULL, any atom
      //  is identified.
      //  DEFAULT values correspond to 'any element' and
      //                 'no alternate location code'
      //  NOTE that " " is not an empty item.
      int   CheckID ( const AtomName aname, const Element elname=NULL,
                      const AltLoc aloc=pstr("") );

      // CheckIDS(..) works exactly like CheckID(..), but it takes
      // the only parameter, the atom ID, which is of the form:
      //    {name} {[element]} {:altcode}
      // Here {} means that the item may be omitted. Any item may be
      // represented by a wildcard '*', which means 'any value'. Just
      // absence of an item means 'empty', which makes sense only for
      // alternate location code. Missing name or element therefore
      // mean 'any name' or 'any element', correspondingly (same as a
      // wildcard). There should be no spaces in ID except for leading
      // spaces; any following space will terminate parsing.
      // The followings are perfectly valid IDs:
      //   CA[C]:A     (carbon C_alpha in location A)
      //   CA[*]:A     (either C_alpha or Ca in location A)
      //   CA:A        (same as above)
      //   CA          (either C_alpha or Ca with no location indicator)
      //   CA[]        (same as above)
      //   CA[C]:      (C_alpha with no location indicator)
      //   [C]         (any carbon with no location indicator)
      //   [C]:*       (any carbon with any location indicator)
      //   *[C]:*      (same as above)
      //   :A          (any atom in location A)
      //   *[*]:A      (same as above)
      //   *[*]:*      (any atom)
      //   *           (any atom with no alternate location indicator)
      int   CheckIDS ( cpstr ID );


      // -------  transform coordinates: x := m*x + v
      void  Transform     ( const mat33 & tm, vect3 & v );
      void  Transform     ( const mat44 & tm );
      void  TransformCopy ( const mat44 & tm,
                            realtype & xx, realtype & yy, realtype & zz );
      void  TransformCopy ( const mat44 & tm, vect3 & xyz );
      void  TransformSet  ( const mat44 & tm,
                            realtype xx, realtype yy, realtype zz );


      // -------  user-defined data handlers
      int   PutUDData ( int UDDhandle, int      iudd );
      int   PutUDData ( int UDDhandle, realtype rudd );
      int   PutUDData ( int UDDhandle, cpstr    sudd );

      int   GetUDData ( int UDDhandle, int      & iudd );
      int   GetUDData ( int UDDhandle, realtype & rudd );
      int   GetUDData ( int UDDhandle, pstr sudd, int maxLen );
      int   GetUDData ( int UDDhandle, pstr     & sudd );


      int   GetIndex()  { return index; }

      virtual void Copy ( PAtom atom );  // without references in
                                          // residues

      void  SetCompactBinary();  // leaves only coordinates in binary files
      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :
      int       index;   // index in the file
      int       nBonds;  // number of bonds in the lowest byte (!)
      PAtomBond Bond;    // atom bonds

      void  InitAtom       ();
      void  FreeMemory     ();
      void  StandardPDBOut ( cpstr Record, pstr S );
      void  GetData        ( cpstr S );
      ERROR_CODE CheckData ( cpstr S );
      void  GetStat        ( realtype   v,
                             realtype & v_min, realtype & v_max,
                             realtype & v_m,   realtype & v_m2 );
      void  _setBonds      ( PPAtom A ); // used only in Residue

  };


  //  ======================  Residue  ==========================

  enum ALTLOC_FLAG  {
    ALF_NoAltCodes    = 0x00000000,
    ALF_EmptyAltLoc   = 0x00000001,
    ALF_NoEmptyAltLoc = 0x00000002,
    ALF_Mess          = 0x00000004,
    ALF_Occupancy     = 0x00000008
  };

  enum SSE_FLAG  {
    SSE_None   = 0,
    SSE_Strand = 1,
    SSE_Bulge  = 2,
    SSE_3Turn  = 3,
    SSE_4Turn  = 4,
    SSE_5Turn  = 5,
    SSE_Helix  = 6
  };

  DefineFactoryFunctions(Residue);

  class MMDB_IMEX Residue : public UDData  {

    friend class Atom;
    friend class Chain;
    friend class Root;

    public :

      ResName  name;            //!< residue name - all spaces cut
      ResName  label_comp_id;   //!< assigned residue name
      ChainID  label_asym_id;   //!< assigned chain Id
      InsCode  insCode;         //!< residue insertion code
      PChain   chain;           //!< reference to chain
      PPAtom   atom;            //!< array of atoms
      int      seqNum;          //!< residue sequence number
      int      label_seq_id;    //!< assigned residue sequence number
      int      label_entity_id; //!< assigned entity id
      int      index;           //!< index in the chain
      int      nAtoms;          //!< number of atoms in the residue
      byte     SSE;             //!< SSE type

      Residue ();
      Residue ( PChain Chain_Owner );
      Residue ( PChain Chain_Owner, const ResName resName,
                int    sqNum,       const InsCode ins );
      Residue ( io::RPStream Object    );
      ~Residue();

      void  SetChain ( PChain Chain_Owner );
      void  SetResID ( const ResName resName, int sqNum,
                       const InsCode ins );
      void  SetChainID ( const ChainID chID );

      void  PDBASCIIAtomDump ( io::RFile f      );
      void  MakeAtomCIF      ( mmcif::PData CIF );

      PChain GetChain();
      PModel GetModel();

      int   GetModelNum   ();
      pstr  GetChainID    ();
      pstr  GetLabelAsymID();
      pstr  GetResName    ();
//      inline pstr  GetResName() { return name; }
      pstr  GetLabelCompID();
      int   GetAASimilarity ( const ResName resName );
      int   GetAASimilarity ( PResidue res );
      realtype GetAAHydropathy();
      void  SetResName      ( const ResName resName );
      int   GetSeqNum       ();
      int   GetLabelSeqID   ();
      int   GetLabelEntityID();
      pstr  GetInsCode      ();
      int   GetResidueNo    ();
      int   GetCenter       ( realtype & x, realtype & y, realtype & z );
      void * GetCoordHierarchy();  // PCMMDBFile

      void  GetAtomStatistics ( RAtomStat AS );
      void  CalAtomStatistics ( RAtomStat AS );

      pstr  GetResidueID ( pstr ResidueID );

      //   GetAltLocations(..) returns the number of different
      // alternative locations in nAltLocs, the locations themselves
      // - in aLoc and the corresponding occupancies - in occupancy.
      //   aLoc and occupancy are allocated dynamically; it is
      // responsibility of the application to deallocate aLoc prior
      // calling GetAltLocations(..) if they were previously allocated.
      // Either, the application is responsible for deallocating aLoc and
      // occupancy after use.
      //   occupancy[i] may return -1.0 if occupancies were not read
      // from coordinate file.
      //   alflag returns ALF_NoAltCodes if no alt codes was found,
      // otherwise the output is decoded according to bits:
      //   ALF_EmptyAltLoc   alternative locations include the
      //                     "no alt loc indicator" ("" for
      //                     Atom::altLoc).
      //                     This means that each atom that has alt locs
      //                     different of "", also includes one marked as
      //                     "".
      //  ALF_NoEmptyAltLoc  alternative locations do not include the
      //                     "no alt loc indicator" ("" for
      //                     Atom::altLoc).
      //                     This means that each atom has either ""
      //                     alt loc or at least two alt locs different
      //                     of "".
      //  ALF_Mess           incorrect residue: it mixes both
      //                     ""-including and not-""-including schemes
      //  ALF_Occupancy      warning that sum of occupancies for alt
      //                     located atoms differ from 1.0 by more
      //                     than 0.01.
      void  GetAltLocations   ( int & nAltLocs, PAltLoc & aLoc,
                                rvector & occupancy, int & alflag );
      int   GetNofAltLocations();

      bool isAminoacid   ();
      bool isNucleotide  ();
      int  isDNARNA      (); // 0(neither),1(DNA),2(RNA)
      bool isSugar       ();
      bool isSolvent     ();
      bool isModRes      ();
      bool isInSelection ( int selHnd );
      bool isNTerminus   ();
      bool isCTerminus   ();

      // -------  checking residue ID
      // CheckID(..) returns 1 if residue is identified, and 0 otherwise.
      //   Parameters:
      //     sname   - pointer to sequence number; if NULL then ignored.
      //     inscode - insertion code; if NULL or '*' then ignored.
      //     resname - residue name; if NULL or '*' then ignored.
      //  IMPORTANT: comparison is case-sensitive.
      //  The residue is considered as identified, if all non-NULL
      //  parameters do match. If all parameters are set NULL, any
      //  residue is identified.
      //  DEFAULT values correspond to 'any residue name' and
      //                 'no insertion code'
      //  NOTE that " " is not an empty item.
      int   CheckID ( int * snum, const InsCode inscode=pstr(""),
                      const ResName resname=NULL );

      // CheckIDS(..) works exactly like CheckID(..), but it takes
      // the only parameter, the residue ID, which is of the form:
      //    {seqnum} {(name)} {.inscode}
      // Here {} means that the item may be omitted. Any item may be
      // represented by a wildcard '*', which means 'any value'. Just
      // absence of a value means 'empty', which is meaningful only for
      // the insertion code. Missing sequence number or residue name
      // therefore mean 'any sequence number' or 'any residue name',
      // correspondingly (same as a wildcard).  There should be no
      // spaces in ID except for leading spaces; any following space will
      // terminate parsing. The followings are perfectly valid IDs:
      //        27(ALA).A   (residue 27A ALA)
      //        27().A      (residue 27A)
      //        27(*).A     (same as above)
      //        27.A        (same as above)
      //        27          (residue 27)
      //        27().       (same as above)
      //        (ALA)       (any ALA without insertion code)
      //        (ALA).      (same as above)
      //        (ALA).*     (any ALA)
      //        *(ALA).*    (any ALA)
      //        .A          (any residue with insertion code A)
      //        *(*).A      (same as above)
      //        *(*).*      (any residue)
      //        *           (any residue with no insertion code)
      int  CheckIDS ( cpstr ID );


      //  --------------------  Extracting atoms  ----------------------

      int  GetNumberOfAtoms ();
      int  GetNumberOfAtoms ( bool countTers );

      PAtom GetAtom ( const AtomName aname, const Element elname=NULL,
                      const AltLoc aloc=cpstr("") );
      PAtom GetAtom ( int atomNo );

      void GetAtomTable  ( PPAtom & atomTable, int & NumberOfAtoms );

      //   GetAtomTable1(..) returns atom table without TER atoms and
      // without NULL atom pointers. NumberOfAtoms returns the actual
      // number of atom pointers in atomTable.
      //   atomTable is allocated withing the function. If it was
      // not set to NULL before calling the function, the latter will
      // attempt to deallocate it first.
      //   The application is responsible for deleting atomTable,
      // however it must not touch atom pointers, i.e. use simply
      // "delete[] atomTable;". Never pass atomTable from GetAtomTable()
      // into this function, unless you set it to NULL before doing that.
      void GetAtomTable1 ( PPAtom & atomTable, int & NumberOfAtoms );


      //  ---------------------  Deleting atoms  -----------------------

      int  DeleteAtom ( const AtomName aname, const Element elname=NULL,
                        const AltLoc aloc=cpstr("") );
      int  DeleteAtom ( int atomNo );
      int  DeleteAllAtoms();

      //   DeleteAltLocs() leaves only alternative location with maximal
      // occupancy, if those are equal or unspecified, the one with
      // "least" alternative location indicator.
      //   The function returns the number of deleted atoms. The atom
      // table remains untrimmed, so that nAtoms are wrong until that
      // is done. Tables are trimmed by FinishStructEdit() or
      // explicitely.
      int  DeleteAltLocs ();

      void TrimAtomTable ();

      //  ----------------------  Adding atoms  ------------------------

      //   AddAtom(..) adds atom to the residue. If residue is associated
      // with a coordinate hierarchy, and atom 'atm' is not, the latter
      // is checked in automatically. If atom 'atm' belongs to any
      // coordinate hierarchy (even though that of the residue), it is
      // *copied* rather than simply taken over, and is checked in.
      //   If residue is not associated with a coordinate hierarchy, all
      // added atoms will be checked in automatically once the residue
      // is checked in.
      int  AddAtom ( PAtom atm );

      //   InsertAtom(..) inserts atom into the specified position of
      // the residue. If residue is associated with a coordinate
      // hierarchy, and atom 'atm' is not, the latter is checked in
      // automatically. If atom 'atm' belongs to any coordinate
      // hierarchy (even though that of the residue), it is *copied*
      // rather than simply taken over, and is checked in.
      //   If residue is not associated with a coordinate hierarchy, all
      // added atoms will be checked in automatically once the residue
      // is checked in.
      int  InsertAtom ( PAtom atm, int position );

      //   This version inserts before the atom with given name. If such
      // name is not found, the atom is appended to the end.
      int  InsertAtom ( PAtom atm, const AtomName aname );

      //  --------------------------------------------------------------

      void  ApplyTransform ( const mat44 & TMatrix );  // transforms all
                                                 // coordinates by
                                                 // multiplying with
                                                 // matrix TMatrix

      void  MaskAtoms   ( PMask Mask );
      void  UnmaskAtoms ( PMask Mask );


      // -------  user-defined data handlers
      int   PutUDData ( int UDDhandle, int      iudd );
      int   PutUDData ( int UDDhandle, realtype rudd );
      int   PutUDData ( int UDDhandle, cpstr    sudd );

      int   GetUDData ( int UDDhandle, int      & iudd );
      int   GetUDData ( int UDDhandle, realtype & rudd );
      int   GetUDData ( int UDDhandle, pstr sudd, int maxLen );
      int   GetUDData ( int UDDhandle, pstr     & sudd );


      bool isMainchainHBond ( PResidue res );

      void  Copy  ( PResidue res );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

    protected :

      int   AtmLen;   // length of atom array
      bool  Exclude;  // used internally

      void  InitResidue  ();
      void  FreeMemory   ();
      int   _AddAtom     ( PAtom atm );
      int   _ExcludeAtom ( int  kndex );  // 1: residue gets empty,
                                          // 0 otherwise
      void  _copy ( PResidue res );
      void  _copy ( PResidue res, PPAtom atm, int & atom_index );
      void  ExpandAtomArray ( int nAdd );
      void  CheckInAtoms ();

  };


  extern realtype  BondAngle ( PAtom A, PAtom B, PAtom C );

}  // namespace mmdb

#endif

