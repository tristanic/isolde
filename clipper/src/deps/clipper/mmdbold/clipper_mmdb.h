/*!  \file clipper_mmdb.h
  Header file for MMDB wrapper
  \ingroup g_mmdb
*/
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


#ifndef CLIPPER_MMDBOLD
#define CLIPPER_MMDBOLD


#include "clipper_mmdb_types.h"
#include "clipper_mmdb_wrapper.h"


namespace clipper {

  // forward definitions
  class DBAtom;
  class DBResidue;
  class DBChain;
  class DBModel;
  class DBManager;
  class DBAtom_selection_inv;
  class DBResidue_selection_inv;


  // ==== abstract base classes ====

  //! Abstract base class for all objects with lists of atoms
  class DBAtom_list
  {
  public:
    virtual ~DBAtom_list() {}
    //! number of atoms in object
    virtual int size() const = 0;
    //! access N'th atom in object 
    virtual DBAtom operator[] ( const int& i ) const = 0;
  };

  //! Abstract base class for all objects with lists of residues
  class DBResidue_list
  {
  public:
    virtual ~DBResidue_list() {}
    //! number of residues in object
    virtual int size() const = 0;
    //! access N'th residue in object 
    virtual DBResidue operator[] ( const int& i ) const = 0;
  };

  //! Abstract base class for all objects with lists of chains
  class DBChain_list
  {
  public:
    virtual ~DBChain_list() {}
    //! number of chains in object
    virtual int size() const = 0;
    //! access N'th chain in object 
    virtual DBChain operator[] ( const int& i ) const = 0;
  };

  //! Abstract base class for all objects with lists of models
  class DBModel_list
  {
  public:
    virtual ~DBModel_list() {}
    //! number of models in object
    virtual int size() const = 0;
    //! access N'th model in object 
    virtual DBModel operator[] ( const int& i ) const = 0;
  };


  // ==== MMDB type wrappers (smart pointers) ====

  //! MMDB atom object wrapper
  /*! This object can be viewed in two ways:
   - As an atom object in the hierarchy, or
   - As a smart pointer to by used in place of mmdb::PCAtom
   This object may be assigned from a PCAtom, and all the methods of
   the mmdb CAtom object are available through the '->' operator, just
   as they would be if you were using a PCAtom.

   But you may also consider it as an object in its own right, using
   the convenience methods provided with no penalty. These methods
   allow access to a subset of the MMDB functionality, but are often
   more concise than their MMDB equivalents, as well as being more
   consistent with the rest of the Clipper interface. */
  class DBAtom : public DBAtom_base
  {
  public:
    //! null constructor
    DBAtom() : ptr(NULL) {}
    //! constructor: from MMDB atom
    DBAtom( mmdb::PCAtom p ) : ptr(p) {}
    bool is_null() const;                //!< test for null object
    bool is_atom() const;                //!< test if this is an atom
    DBResidue residue() const;           //!< get parent
    void set( const DBAtom_base& a );    //!< set all properties
    int index() const;  //!< return the index of this object in parent

    // standard atom properties
    String type() const;
    String element() const;
    String altconf() const;
    Coord_orth coord_orth() const;
    ftype occupancy() const;
    ftype u_iso() const;
    U_aniso_orth u_aniso_orth() const;
    Sig_Coord_orth sig_coord_orth() const;
    ftype sig_occupancy() const;
    ftype sig_u_iso() const;
    Sig_U_aniso_orth sig_u_aniso_orth() const;
    void set_type( const String& n );
    void set_element( const String& n );
    void set_altconf( const String& n );
    void set_coord_orth( const Coord_orth& v );
    void set_occupancy( const ftype& v );
    void set_u_iso( const ftype& v );
    void set_u_aniso_orth( const U_aniso_orth& v );
    void set_sig_coord_orth( const Sig_Coord_orth& s );
    void set_sig_occupancy( const ftype& s );
    void set_sig_u_iso( const ftype& s );
    void set_sig_u_aniso_orth( const Sig_U_aniso_orth& s );
    // other atom properties
    int serial_num() const;         //!< get atom serial number
    String charge() const;          //!< get atom charge

    //! direct access to MMDB Atom
    mmdb::PCAtom pcatom() { return ptr; }
    //! direct access to MMDB Atom
    mmdb::PCAtom operator -> () { return pcatom(); }
  private:
    mmdb::PCAtom ptr;
  };

  //! MMDB residue object wrapper
  class DBResidue : public DBAtom_list, public DBResidue_base
  {
  public:
    //! null constructor
    DBResidue() : ptr(NULL) {}
    //! constructor: from MMDB res
    DBResidue( mmdb::PCResidue p ) : ptr(p) {}
    bool is_null() const;                 //!< test for null object
    DBChain chain() const;                //!< get parent
    void set( const DBResidue_base& r );  //!< set all properties
    int index() const;  //!< return the index of this object in parent

    // standard residue properties
    String type() const;
    int seqnum() const;
    String inscode() const;
    void set_type( const String& n );
    void set_seqnum( const int& n );
    void set_inscode( const String& n );

    DBAtom add_atom( const DBAtom_base& a );  //!< add an atom to the residue
    //! access atom within residue by type
    DBAtom atom( const String& type ) const;

    int size() const;
    DBAtom operator[] ( const int& i ) const;
    //! direct access to MMDB Residue
    mmdb::PCResidue pcresidue() { return ptr; }
    //! direct access to MMDB Residue
    mmdb::PCResidue operator -> () { return pcresidue(); }
  private:
    mmdb::PCResidue ptr;
  };

  //! MMDB chain object wrapper
  class DBChain : public DBResidue_list, public DBChain_base
  {
  public:
    //! null constructor
    DBChain() : ptr(NULL) {}
    //! constructor: from MMDB chain
    DBChain( mmdb::PCChain p ) : ptr(p) {}
    bool is_null() const;               //!< test for null object
    DBModel model() const;              //!< get parent
    void set( const DBChain_base& c );  //!< set all properties
    int index() const;  //!< return the index of this object in parent

    // chain properties
    String id() const;
    void set_id( const String& n );

    DBResidue add_residue( const DBResidue_base& r );  //!< add residue at end
    //! access residue within chain by seqnum
    DBResidue residue( const int& seqnum ) const;

    int size() const;
    DBResidue operator[] ( const int& i ) const;
    //! direct access to MMDB Chain
    mmdb::PCChain pcchain() { return ptr; }
    //! direct access to MMDB Chain
    mmdb::PCChain operator -> () { return pcchain(); }
  private:
    mmdb::PCChain ptr;
  };

  //! MMDB model object wrapper
  class DBModel : public DBChain_list, public DBModel_base
  {
  public:
    //! null constructor
    DBModel() : ptr(NULL) {}
    //! constructor: from MMDB model
    explicit DBModel( mmdb::PCModel p ) : ptr(p) {}
    bool is_null() const;               //!< test for null object
    DBManager manager() const;          //!< get parent
    void set( const DBModel_base& m );  //!< set all properties
    int index() const;  //!< return the index of this object in parent

    // model properties
    String id() const;
    void set_id( const String& n );

    DBChain add_chain( const DBChain_base& c );  //!< add a chain to the model
    //! access chain within residue by type
    DBChain chain( const String& id ) const;

    int size() const;
    DBChain operator[] ( const int& i ) const;
    //! direct access to MMDB Model
    mmdb::PCModel pcmodel() { return ptr; }
    //! direct access to MMDB Model
    mmdb::PCModel operator -> () { return pcmodel(); }
  private:
    mmdb::PCModel ptr;
  };


  // ==== List (selection) types ====

  //! Atom selection object
  /*! An atom selection may hold an abitrary list of atoms in any
    order, including duplicates. However, if constructed using a
    selection function, the selection will always be sorted and
    unique. Application of logical operators to a selection also
    reduces it to a sorted, unique list. */
  class DBAtom_selection : public DBAtom_list
  {
  public:
    //! null constructor
    DBAtom_selection() {}
    //! constructor: from MMDB selection table
    DBAtom_selection( mmdb::PPCAtom p, int n );
    int size() const { return list.size(); }
    DBAtom operator[] ( const int& i ) const { return DBAtom( list[i] ); }
    //! add an atom to the end of the selection
    void add_atom( DBAtom a );
    friend DBAtom_selection operator&
      ( const DBAtom_selection& a1, const DBAtom_selection& a2 );
    friend DBAtom_selection operator|
      ( const DBAtom_selection& a1, const DBAtom_selection& a2 );
    friend DBAtom_selection operator^
      ( const DBAtom_selection& a1, const DBAtom_selection& a2 );
    friend DBAtom_selection_inv operator! ( const DBAtom_selection& a );
    //! direct access to MMDB Atoms
    mmdb::PPCAtom ppcatom() { return &list[0]; }
  private:
    std::vector<mmdb::PCAtom> list;
  };

  //! Inverse atom selection object \internal
  /*! A residue selection may hold an abitrary list of residues in any
    order, including duplicates. However, if constructed using a
    selection function, the selection will always be sorted and
    unique. Application of logical operators to a selection also
    reduces it to a sorted, unique list. */
  class DBAtom_selection_inv : protected DBAtom_selection
  {
  public:
    explicit DBAtom_selection_inv( DBAtom_selection s ) : DBAtom_selection(s) {}
    friend DBAtom_selection_inv operator& ( const DBAtom_selection_inv& a1, const DBAtom_selection_inv& a2 ) { return !(!a1 | !a2); }
    friend DBAtom_selection_inv operator| ( const DBAtom_selection_inv& a1, const DBAtom_selection_inv& a2 ) { return !(!a1 & !a2); }
    friend DBAtom_selection operator^ ( const DBAtom_selection_inv& a1, const DBAtom_selection_inv& a2 ) { return (!a1 ^ !a2); }
    friend DBAtom_selection operator& ( const DBAtom_selection& a1, const DBAtom_selection_inv& a2 ) { return a1 ^ (a1 & !a2); }
    friend DBAtom_selection operator& ( const DBAtom_selection_inv& a1, const DBAtom_selection& a2 ) { return a2 ^ (a2 & !a1); }
    friend DBAtom_selection operator! ( const DBAtom_selection_inv& a ) { return a; }
  };

  //! Residue selection object
  class DBResidue_selection : public DBResidue_list
  {
  public:
    //! null constructor
    DBResidue_selection() {}
    //! constructor: from MMDB selection table
    DBResidue_selection( mmdb::PPCResidue p, int n );
    int size() const { return list.size(); }
    DBResidue operator[] ( const int& i ) const { return DBResidue( list[i] ); }
    //! add a residue to the end of the selection
    void add_residue( DBResidue a );
    friend DBResidue_selection operator&
      ( const DBResidue_selection& a1, const DBResidue_selection& a2 );
    friend DBResidue_selection operator|
      ( const DBResidue_selection& a1, const DBResidue_selection& a2 );
    friend DBResidue_selection operator^
      ( const DBResidue_selection& a1, const DBResidue_selection& a2 );
    friend DBResidue_selection_inv operator! ( const DBResidue_selection& a );
    //! direct access to MMDB Residues
    mmdb::PPCResidue ppcresidue() { return &list[0]; }
  private:
    std::vector<mmdb::PCResidue> list;
  };

  //! Inverse residue selection object \internal
  class DBResidue_selection_inv : protected DBResidue_selection
  {
  public:
    explicit DBResidue_selection_inv( DBResidue_selection s ) : DBResidue_selection(s) {}
    friend DBResidue_selection_inv operator& ( const DBResidue_selection_inv& a1, const DBResidue_selection_inv& a2 ) { return !(!a1 | !a2); }
    friend DBResidue_selection_inv operator| ( const DBResidue_selection_inv& a1, const DBResidue_selection_inv& a2 ) { return !(!a1 & !a2); }
    friend DBResidue_selection operator^ ( const DBResidue_selection_inv& a1, const DBResidue_selection_inv& a2 ) { return (!a1 ^ !a2); }
    friend DBResidue_selection operator& ( const DBResidue_selection& a1, const DBResidue_selection_inv& a2 ) { return a1 ^ (a1 & !a2); }
    friend DBResidue_selection operator& ( const DBResidue_selection_inv& a1, const DBResidue_selection& a2 ) { return a2 ^ (a2 & !a1); }
    friend DBResidue_selection operator! ( const DBResidue_selection_inv& a ) { return a; }
  };

  //! Chain selection object
  class DBChain_selection : public DBChain_list
  {
  public:
    //! null constructor
    DBChain_selection() {}
    //! constructor: from MMDB selection table
    DBChain_selection( mmdb::PPCChain p, int n );
    int size() const { return list.size(); }
    DBChain operator[] ( const int& i ) const { return DBChain( list[i] ); }
    //! add a chain to the end of the selection
    void add_chain( DBChain a );
    //! direct access to MMDB Chains
    mmdb::PPCChain ppcchain() { return &list[0]; }
  private:
    std::vector<mmdb::PCChain> list;
  };

  //! Model selection object
  class DBModel_selection : public DBModel_list
  {
  public:
    //! null constructor
    DBModel_selection() {}
    //! constructor: from MMDB selection table
    DBModel_selection( mmdb::PPCModel p, int n );
    int size() const { return list.size(); }
    DBModel operator[] ( const int& i ) const { return DBModel( list[i] ); }
    //! add a model to the end of the selection
    void add_model( DBModel a );
    //! direct access to MMDB Models
    mmdb::PPCModel ppcmodel() { return &list[0]; }
  private:
    std::vector<mmdb::PCModel> list;
  };


  // ==== MMDB manager wrapper ====


  //! MMDB object wrapper
  /*! This is a smart pointer to an MMDB-manager. It should not be
    used directly, instead use the derived clipper::MMDB class. This
    class is returned when you ask for a reference to an MMDB, e.g. as
    the parent of a DBModel, however it should not be instantiated
    directly.
  */

  class DBManager : public DBModel_list
  {
  public:
    //! null constructor
    DBManager() {}
    //! constructor: from MMDB selection table
    explicit DBManager( mmdb::XMMDBManager* p ) : ptr(p) {}

    //! access N'th model (indexed from 1)
    DBModel model( const int i = 1);
    //! add a model
    DBModel add_model( const DBModel_base& m );
    //! finalise structure edit
    void finalise_edit();
    //! Select atoms by string selection function
    DBAtom_selection select_atoms( const String& s );
    //! Select residues by string selection function
    DBResidue_selection select_residues( const String& s );
    //! Select chains by string selection function
    DBChain_selection select_chains( const String& s );
    //! Select models by string selection function
    DBModel_selection select_models( const String& s );
    //! Select atoms neighbouring some atom selection
    DBAtom_selection select_atoms_near( DBAtom_selection& s, const ftype& r1, const ftype& r2 );
    //! Select residues neighbouring some atom selection
    DBResidue_selection select_residues_near( DBAtom_selection& s, const ftype& r1, const ftype& r2 );
    //! Select atoms by sphere about some centre
    DBAtom_selection select_atoms_sphere( const Coord_orth& c, const ftype& r );
    //! Select residues by sphere about some centre
    DBResidue_selection select_residues_sphere( const Coord_orth& c, const ftype& r );
    //! Select atoms by cylinder
    DBAtom_selection select_atoms_cylinder( const Coord_orth& c1, const Coord_orth& c2, const ftype& r );
    //! Select atoms by serial number
    DBAtom_selection select_atoms_serial( const int i1 = 0, const int i2 = 0 );
    int size() const;
    DBModel operator[] ( const int& i ) const;
    //! direct access to MMDB Manager
    mmdb::PCMMDBManager pcmmdbmanager() const { return ptr; }
    //! direct access to MMDB Manager
    mmdb::PCMMDBManager operator -> () { return pcmmdbmanager(); }

  protected:
    mmdb::XMMDBManager* ptr;
  };


  //! MMDB object
  /*! This is the Clipper form of an MMDB-manager. Unlike the other
    data objects, it does not use external i/o objects, the i/o
    methods are built in. All MMDB functionality may be accessed
    through the -> operator, however the new read_file() and
    write_file() methods should be used in preference since they set
    the spacegroup and cell properties correctly.
  */
  class MMDB : public DBManager
  {
  public:
    enum TYPE { Default=-1, PDB, CIF, Binary };
    //! null constructor
    MMDB();
    //! constructor: from spacegroup and cell
    MMDB( const Spacegroup& spacegroup, const Cell& cell );
    //! override copy constructor
    MMDB( const MMDB& c );
    //! destructor
    ~MMDB();
    //! initialise: from spacegroup and cell
    void init( const Spacegroup& spacegroup, const Cell& cell );

    //! test if object has been initialised
    bool is_null() const;

    //! get spacegroup
    const Spacegroup& spacegroup() const;
    //! get cell
    const Cell& cell() const;
    //! import model from file
    void read_file( const String& file );
    //! export model to file
    void write_file( const String& file, TYPE type = Default );

    // inherited functions listed for documentation purposes
    //-- DBModel model( const int i = 1);
    //-- DBModel add_model( const DBModel_base& m );
    //-- void finalise_edit();
    //-- DBAtom_selection select_atoms( const String& s );
    //-- DBResidue_selection select_residues( const String& s );
    //-- DBChain_selection select_chains( const String& s );
    //-- DBModel_selection select_models( const String& s );
    //-- DBAtom_selection select_atoms_near( DBAtom_selection& s, const ftype& r1, const ftype& r2 );
    //-- DBResidue_selection select_residues_near( DBAtom_selection& s, const ftype& r1, const ftype& r2 );
    //-- DBAtom_selection select_atoms_sphere( const Coord_orth& c, const ftype& r );
    //-- DBResidue_selection select_residues_sphere( const Coord_orth& c, const ftype& r );
    //-- DBAtom_selection select_atoms_cylinder( const Coord_orth& c1, const Coord_orth& c2, const ftype& r );
    //-- DBAtom_selection select_atoms_serial( const int i1 = 0, const int i2 = 0 );
    //-- int size() const;
    //-- DBModel operator[] ( const int& i ) const;
    //-- mmdb::PCMMDBManager pcmmdbmanager() const { return ptr; }
    //-- mmdb::PCMMDBManager operator -> () { return pcmmdbmanager(); }

    void debug() const;
  private:
    Spacegroup spacegroup_;
    Cell cell_;
  };


} // namespace clipper

#endif
