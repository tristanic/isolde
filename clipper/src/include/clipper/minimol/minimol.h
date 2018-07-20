/*! \file minimol.h
 Header file for atomic model types */

//C Copyright (C) 2000-2006 Kevin Cowtan and University of York
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA


#ifndef CLIPPER_MINIMOL
#define CLIPPER_MINIMOL


#include "../core/coords.h"
#include "minimol_data.h"
#include "../imex.h"

namespace clipper {


 // Forward definitions
 class MAtomIndex;
 class MAtomIndexSymmetry;


 //! dummy namespace to hold search modes
 namespace MM {
   enum MODE { UNIQUE, ANY };
   enum COPY { COPY_NONE, COPY_M, COPY_P, COPY_MP, COPY_C, COPY_MC, COPY_PC, COPY_MPC, MEMBERS=COPY_M, PROPERTIES=COPY_P, CHILDREN=COPY_C };
 }


 //! MiniMol atom object
 /*! The MiniMol atom is derived from the basic clipper::Atom, with
   the addition of an 'id', which is a unique identifier within a
   monomer in accordance with the mmCIF definition.

   In addition, it is a clipper::PropertyManager, which means you can
   add labelled properties of any type to the object. These may be
   simple strings, or complex objects such as maps, function objects,
   or whatever.

   The most commonly used properties are:
   - "CID" The original CID of this atom in an MMDB heirarchy.

   The id() is the unique key which identifies an atom. */
 class CLIPPER_IMEX MAtom : public Atom, public PropertyManager
 {
 public:
   MAtom() {}  //!< null constructor
   MAtom( const clipper::Atom& atom );  //!< constructor: from clipper::Atom

   //! get atom ID, e.g. " N  ", " CA ", " CG1", " CA :A"
   const String& id() const { return id_; }
   void set_id( const String& s );      //!< set atom ID

   //! get atom name, i.e. the ID, omitting any alternate conformation code
   String name() const { return id_.substr(0,4); }
   void set_name( const String s, const String altconf="" );  //!< set full id

   //-- const String& element() const;  //!< get element e.g. H, C, Zn2+
   //-- const Coord_orth& coord_orth() const;      //!< get orth coordinate
   //-- const ftype& occupancy() const;            //!< get occupancy
   //-- const ftype& u_iso() const;                //!< get isotropic U
   //-- const U_aniso_orth& u_aniso_orth() const;  //!< get anisotropic U
   //-- void set_element( const String& s );             //!< set element
   //-- void set_coord_orth( const Coord_orth& s );      //!< set coord_orth
   //-- void set_occupancy( const ftype& s );            //!< set occupancy
   //-- void set_u_iso( const ftype& s );                //!< set u_iso
   //-- void set_u_aniso_orth( const U_aniso_orth& s );  //!< set u_aniso
   //-- void transform( const RTop_orth rt );  //!< apply transform to object
   const Atom& atom() const { return *this; }  //!< explicitly get atom
   Atom& atom() { return *this; }              //!< explicitly set atom

   //! configureable copy function
   MAtom& copy( const MAtom& other, const MM::COPY& mode );

   static String id_tidy( const String& id );  //!< convert ID to std format
   static bool id_match( const String& id1, const String& id2, const MM::MODE& mode );  //!< convert ID to std format
 private:
   String id_;
 };


 //! MiniMol monomer (e.g. residue) object
 /*! The MiniMol monomer object contains a list of clipper::MAtom.

   It has two properties: a sequence number and a type. The sequence
   number need not reflect the order in which the monomers are stored
   in a polymer. MResidue is an alias for MMonomer.

   In addition, it is a clipper::PropertyManager, which means you can
   add labelled properties of any type to the object. These may be
   simple strings, or complex objects such as maps, function objects,
   or whatever.

   The most commonly used properties are:
   - "CID" The original CID of this atom in an MMDB heirarchy.

   The id() is the unique key which identifies a monomer. */
 class CLIPPER_IMEX MMonomer : public PropertyManager
 {
 public:
   const String& id() const { return id_; }  //!< get monomer ID
   void set_id( const String& s );           //!< set monomer ID

   const String& type() const { return type_; }   //!< get monomer type
   void set_type( const String& s );  //!< set monomer type, e.g. LYS, VAL, G

   int seqnum() const { return id_.i(); }  //!< get monomer seq number
   void set_seqnum( const int s, const String inscode="" ); //!< set full id

   // the following methods are similar for all levels of the hierarchy
   Atom_list atom_list() const;           //!< return list of contained atoms
   void transform( const RTop_orth rt );  //!< apply transformation to object
   //! number of atoms in monomer
   int size() const { return children.size(); }
   //! get atom
   const MAtom& operator[] ( const int& i ) const { return children[i]; }
   //! set atom
   MAtom& operator[] ( const int& i ) { return children[i]; }
   //! get atom by id
   const MAtom& find( const String& n, const MM::MODE mode=MM::UNIQUE ) const;
   //! set atom by id
   MAtom& find( const String& n, const MM::MODE mode=MM::UNIQUE );
   //! create selection
   MMonomer select( const String& sel, const MM::MODE mode=MM::UNIQUE ) const;
   //! get child indices matching a selection criteria
   std::vector<int> select_index( const String& sel, const MM::MODE mode=MM::UNIQUE ) const;
   //! lookup atom by id
   int lookup( const String& str, const MM::MODE& mode ) const;
   void insert( const MAtom& add, int pos=-1 );  //!< add atom

   //! and operator
   friend CLIPPER_IMEX MMonomer operator& ( const MMonomer& m1, const MMonomer& m2 );
   //! or operator
   friend CLIPPER_IMEX MMonomer operator| ( const MMonomer& m1, const MMonomer& m2 );

   //! configureable copy function
   MMonomer& copy( const MMonomer& other, const MM::COPY& mode );

   static String id_tidy( const String& id );  //!< convert ID to std format
   static bool id_match( const String& id1, const String& id2, const MM::MODE& mode );  //!< convert ID to std format

   //! Rotamer library type
   enum TYPE { Default, Dunbrack, Richardson };
   //! UTILITY: Build carbonyl oxygen, given next residue in chain
   void protein_mainchain_build_carbonyl_oxygen( const MMonomer& next );
   //! UTILITY: Build carbonyl oxygen, without next residue in chain
   void protein_mainchain_build_carbonyl_oxygen();
   //! UTILITY: get number of rotamers for protein sidechain
   int protein_sidechain_number_of_rotamers( TYPE t = default_type_ ) const;
   int protein_sidechain_number_of_rotomers() const { return protein_sidechain_number_of_rotamers(); }
   //! UTILITY: build numbered rotamer for protein sidechain
   ftype protein_sidechain_build_rotamer( const int& n, TYPE t = default_type_ );
   ftype protein_sidechain_build_rotomer( const int& n ) { return protein_sidechain_build_rotamer( n ); }
   //! UTILITY: test if two peptide are adjacent
   static bool protein_peptide_bond( const MMonomer& m1, const MMonomer& m2, ftype r = 1.5 );
   //! UTILITY: return Ramachandran phi, or NaN if atoms missing
   static double protein_ramachandran_phi( const MMonomer& m1, const MMonomer& m2 );
   //! UTILITY: return Ramachandran psi, or NaN if atoms missing
   static double protein_ramachandran_psi( const MMonomer& m1, const MMonomer& m2 );

   static TYPE& default_type() { return default_type_; }
 private:
   typedef MAtom CHILDTYPE;
   std::vector<CHILDTYPE> children;
   String id_, type_;
   static TYPE default_type_;
   static int rotamer_find( String res, int rota, TYPE t );
 };


 //! MiniMol polymer (e.g. chain) object
 /*! The MiniMol polymer object has one property: an identifying name.

   It contains a list of clipper::MMonomer.

   In addition, it is a clipper::PropertyManager, which means you can
   add labelled properties of any type to the object. These may be
   simple strings, or complex objects such as maps, function objects,
   or whatever.

   The most commonly used properties are:
   - "CID" The original CID of this atom in an MMDB heirarchy.

   The id() is the unique key which identifies a polymer. */
 class CLIPPER_IMEX MPolymer : public PropertyManager
 {
 public:
   const String& id() const { return id_; }  //!< get polymer ID
   void set_id( const String& s );           //!< set polymer ID

   // the following methods are similar for all levels of the hierarchy
   Atom_list atom_list() const;           //!< return list of contained atoms
   void transform( const RTop_orth rt );  //!< apply transformation to object
   //! number of monomers in polymer
   int size() const { return children.size(); }
   //! get monomer
   const MMonomer& operator[] ( const int& i ) const { return children[i]; }
   //! set monomer
   MMonomer& operator[] ( const int& i ) { return children[i]; }
   //! get monomer by id
   const MMonomer& find( const String& n, const MM::MODE mode=MM::UNIQUE ) const;
   //! set monomer by id
   MMonomer& find( const String& n, const MM::MODE mode=MM::UNIQUE );
   //! create selection
   MPolymer select( const String& sel, const MM::MODE mode=MM::UNIQUE ) const;
   //! get child indices matching a selection criteria
   std::vector<int> select_index( const String& sel, const MM::MODE mode=MM::UNIQUE ) const;
   //! lookup monomer by id
   int lookup( const String& str, const MM::MODE& mode ) const;
   void insert( const MMonomer& add, int pos=-1 );  //!< add monomer

   //! and operator
   friend CLIPPER_IMEX MPolymer operator& ( const MPolymer& m1, const MPolymer& m2 );
   //! or operator
   friend CLIPPER_IMEX MPolymer operator| ( const MPolymer& m1, const MPolymer& m2 );

   //! configureable copy function
   MPolymer& copy( const MPolymer& other, const MM::COPY& mode );

   static String id_tidy( const String& id );  //!< convert ID to std format
   static bool id_match( const String& id1, const String& id2, const MM::MODE& mode );  //!< convert ID to std format
 private:
   typedef MMonomer CHILDTYPE;
   std::vector<CHILDTYPE> children;
   String id_;
 };


 //! MiniMol model object
 /*! The MiniMol model object constains a list of clipper::MPolymer.

   In addition, it is a clipper::PropertyManager, which means you can
   add labelled properties of any type to the object. These may be
   simple strings, or complex objects such as maps, function objects,
   or whatever.

   The most commonly used properties are:
   - "CID" The original CID of this atom in an MMDB heirarchy.
 */
 class CLIPPER_IMEX MModel : public PropertyManager
 {
 public:
   // the following methods are similar for all levels of the hierarchy
   Atom_list atom_list() const;           //!< return list of contained atoms
   void transform( const RTop_orth rt );  //!< apply transformation to object
   //! number of polymers in model
   int size() const { return children.size(); }
   //! get polymer
   const MPolymer& operator[] ( const int& i ) const { return children[i]; }
   //! set polymer
   MPolymer& operator[] ( const int& i ) { return children[i]; }
   //! get polymer by id
   const MPolymer& find( const String& n, const MM::MODE mode=MM::UNIQUE ) const;
   //! set polymer by id
   MPolymer& find( const String& n, const MM::MODE mode=MM::UNIQUE );
   //! create selection
   MModel select( const String& sel, const MM::MODE mode=MM::UNIQUE ) const;
   //! get child indices matching a selection criteria
   std::vector<int> select_index( const String& sel, const MM::MODE mode=MM::UNIQUE ) const;
   //! lookup polymer by id
   int lookup( const String& str, const MM::MODE& mode ) const;
   void insert( const MPolymer& add, int pos=-1 );  //!< add polymer

   //! and operator
   friend CLIPPER_IMEX MModel operator& ( const MModel& m1, const MModel& m2 );
   //! or operator
   friend CLIPPER_IMEX MModel operator| ( const MModel& m1, const MModel& m2 );

   //! configureable copy function
   MModel& copy( const MModel& other, const MM::COPY& mode );

   //! Return atom by MAtomIndex
   const MAtom& atom( const MAtomIndex& index ) const;
   //! Return atom by MAtomIndex
   MAtom& atom( const MAtomIndex& index );
   //! Select and return MAtomIndex values
   std::vector<MAtomIndex> select_atom_index( const String& sel, const MM::MODE mode=MM::UNIQUE ) const;
 private:
   typedef MPolymer CHILDTYPE;
   std::vector<CHILDTYPE> children;
 };


 //! MiniMol lightweight coordinate model object
 /*! A MiniMol object is a model (clipper::MModel) embedded in a
   crystal frame, i.e. with additional spacegroup and cell
   information.

   The design of this object was inspired by and contributed to by
   Paul Emsley. */
 class CLIPPER_IMEX MiniMol : public MModel
 {
 public:
   enum MODE { UNIQUE, ANY };  //!< mode to use when matching IDs

   //! null constructor
   MiniMol();
   //! constructor: from spacegroup and cell
   MiniMol( const Spacegroup& spacegroup, const Cell& cell );
   //! initialiser: from spacegroup and cell
   void init( const Spacegroup& spacegroup, const Cell& cell );
   //! get the cell
   const Cell& cell() const { return cell_; }
   //! get the spacegroup
   const Spacegroup& spacegroup() const { return spacegroup_; }
   const MModel& model() const { return *this; }  //!< explicitly get model
   MModel& model() { return *this; }              //!< explicitly set model

   //! Return symmetry atom by MAtomIndexSymmetry
   MAtom symmetry_atom( const MAtomIndexSymmetry& index );

   bool is_null() const;  //!< test for null model
 private:
   Spacegroup spacegroup_;
   Cell cell_;
 };

 //! MResidue: an alternative name for an MMonomer
 typedef MMonomer MResidue;

 //! MChain: an alternative name for an MPolymer
 typedef MPolymer MChain;


} // namespace clipper

#endif

