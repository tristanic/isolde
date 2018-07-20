/*!  \file clipper_mmdb.h
  Header file for MMDB wrapper
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


#ifndef CLIPPER_MMDB
#define CLIPPER_MMDB


#include "../core/coords.h"

#include <mmdb2/mmdb_manager.h>
#include "../imex.h"


namespace clipper {

  namespace mmdb {
    typedef ::mmdb::Atom CAtom;
    typedef ::mmdb::Residue CResidue;
    typedef ::mmdb::Chain CChain;
    typedef ::mmdb::Model CModel;
    typedef ::mmdb::Manager CMMDBManager;

    typedef ::mmdb::PAtom PCAtom;
    typedef ::mmdb::PResidue PCResidue;
    typedef ::mmdb::PChain PCChain;
    typedef ::mmdb::PModel PCModel;
    typedef ::mmdb::PManager PCMMDBManager;

    typedef ::mmdb::PPAtom PPCAtom;
    typedef ::mmdb::PPResidue PPCResidue;
    typedef ::mmdb::PPChain PPCChain;
    typedef ::mmdb::PPModel PPCModel;
    typedef ::mmdb::PPManager PPCMMDBManager;
  }

  // PDB 'sigma' classes (should be depracated)

  //! Standard deviation of orthognal coordinates
  /*! see clipper::Coord_orth
    \note In my view this is a stupid definition, but its in the PDB.
    \ingroup g_mmdb */
  class CLIPPER_IMEX Sig_Coord_orth : public Vec3<>
  {
  public:
    Sig_Coord_orth() {}    //!< null constructor
    explicit Sig_Coord_orth( const Vec3<>& v ) :
      Vec3<>( v ) {}  //!< constructor: copy/convert
    Sig_Coord_orth( const ftype& sigx, const ftype& sigy, const ftype& sigz ) :
      Vec3<>( sigx, sigy, sigz ) {}    //!< constructor: from sig(x,y,z)
    const ftype& sigx() const { return (*this)[0]; }  //!< get sigx
    const ftype& sigy() const { return (*this)[1]; }  //!< get sigy
    const ftype& sigz() const { return (*this)[2]; }  //!< get sigz
  };

  //! Standard deviation of anisotropic atomic displacement parameters
  /*! see clipper::U_aniso_orth
    \note In my view this is a stupid definition, but its in the PDB.
    \ingroup g_mmdb */
  class CLIPPER_IMEX Sig_U_aniso_orth : public Mat33sym<>
  {
  public:
    //! null constructor
    Sig_U_aniso_orth() {};
    //! constructor: from Mat33sym
    explicit Sig_U_aniso_orth( const Mat33sym<>& m ) : Mat33sym<>(m) {}
    //! constructor: from sig_Uij
    Sig_U_aniso_orth( const ftype& su11, const ftype& su22, const ftype& su33,
                 const ftype& su12, const ftype& su13, const ftype& su23 ) :
      Mat33sym<>( su11, su22, su33, su12, su13, su23 ) {}
  };


  //! MMDB atom object
  /*! This class is a trivial derivation of the corresponding MMDB,
    providing access in terms of Clipper types. Thus, when you such
    access, simply cast you MMDB object reference to this type to
    access the additional functions.
    For full documentation see: http://www.ebi.ac.uk/~keb/ */
  class CLIPPER_IMEX MMDBAtom : public mmdb::CAtom
  {
  public:
    //! null constructor
    MMDBAtom() {}
    //! constructor: from MMDB atom
    MMDBAtom( const mmdb::CAtom& a ) : mmdb::CAtom(a) {}

    // standard atom properties
    String id() const;              //!< Atom id, e.g. CA, CB, CH3
    String element() const;         //!< Atom element, e.g. C, H, Zn2+
    Coord_orth coord_orth() const;  //!< Atom coordinate (orthogonal Angstroms)
    ftype occupancy() const;        //!< Atom occupancy (0...1)
    ftype u_iso() const;            //!< Atom isotropic U
    U_aniso_orth u_aniso_orth() const;  //!< Atom anisotropic U (orthogonal As)
    void set_id( const String& n );                  //!< set id
    void set_element( const String& n );             //!< set element
    void set_coord_orth( const Coord_orth& v );      //!< set coordinate
    void set_occupancy( const ftype& v );            //!< set occupancy
    void set_u_iso( const ftype& v );                //!< set iso U
    void set_u_aniso_orth( const U_aniso_orth& v );  //!< set aniso U
    Sig_Coord_orth sig_coord_orth() const;  // stupid sigmas
    ftype sig_occupancy() const;
    ftype sig_u_iso() const;
    Sig_U_aniso_orth sig_u_aniso_orth() const;
    void set_sig_coord_orth( const Sig_Coord_orth& s );
    void set_sig_occupancy( const ftype& s );
    void set_sig_u_iso( const ftype& s );
    void set_sig_u_aniso_orth( const Sig_U_aniso_orth& s );
    // other atom properties
    String altconf() const;         //!< get atom alternate conformation code
    int serial_num() const;         //!< get atom serial number
    String charge() const;          //!< get atom charge
  };

  //! MMDB residue object wrapper
  /*! This class is a trivial derivation of the corresponding MMDB,
    providing access in terms of Clipper types. Thus, when you such
    access, simply cast you MMDB object reference to this type to
    access the additional functions.
    For full documentation see: http://www.ebi.ac.uk/~keb/ */
  class CLIPPER_IMEX MMDBResidue : public mmdb::CResidue
  {
  public:
    //! null constructor
    MMDBResidue() {}
    //! constructor: from MMDB residue
    MMDBResidue( const mmdb::CResidue& a ) : mmdb::CResidue(a) {}
 
    // standard residue properties
    String type() const;
    int seqnum() const;
    String inscode() const;
    void set_type( const String& n );
    void set_seqnum( const int& n );
    void set_inscode( const String& n );
  };

  //! MMDB chain object wrapper
  /*! This class is a trivial derivation of the corresponding MMDB,
    providing access in terms of Clipper types. Thus, when you such
    access, simply cast you MMDB object reference to this type to
    access the additional functions.
    For full documentation see: http://www.ebi.ac.uk/~keb/ */
  class CLIPPER_IMEX MMDBChain : public mmdb::CChain
  {
  public:
    //! null constructor
    MMDBChain() {}
    //! constructor: from MMDB chain
    MMDBChain( const mmdb::CChain& a ) : mmdb::CChain(a) {}

    // chain properties
    String id() const;
    void set_id( const String& n );
  };

  //! MMDB model object wrapper
  /*! This class is a trivial derivation of the corresponding MMDB,
    providing access in terms of Clipper types. Thus, when you such
    access, simply cast you MMDB object reference to this type to
    access the additional functions.
    For full documentation see: http://www.ebi.ac.uk/~keb/ */
  class CLIPPER_IMEX MMDBModel : public mmdb::CModel
  {
  public:
    //! null constructor
    MMDBModel()  {}
    //! constructor: from MMDB model
    MMDBModel( const mmdb::CModel& a ) : mmdb::CModel(a) {}

    // model properties
    String id() const;
    void set_id( const String& n );
  };

  //! MMDB manager wrapper
  /*! This class is a trivial derivation of the corresponding MMDB,
    providing access in terms of Clipper types. Thus, when you such
    access, simply cast you MMDB object reference to this type to
    access the additional functions.
    For full documentation see: http://www.ebi.ac.uk/~keb/ */
  class CLIPPER_IMEX MMDBManager : public mmdb::CMMDBManager
  {
  public:
    enum TYPE { Default=-1, PDB, CIF, Binary };
    MMDBManager();   //! null constructor
    ~MMDBManager();  //! destructor
    Spacegroup spacegroup() const;                        //!< get spacegroup
    Cell cell() const;                                    //!< get cell
    void set_spacegroup( const Spacegroup& spacegroup );  //!< set spacegroup
    void set_cell( const Cell& cell );                    //!< set cell
    //-- int ReadCoorFile( char* name );                       //!< For file i/o see http://www.ebi.ac.uk/~keb/cldoc/object/cl_obj_rdwr.html
    //-- void select( int hnd, int typ, char* str, int key );  //!< For selection functions see http://www.ebi.ac.uk/~keb/cldoc/object/cl_obj_selfnc.html
    //-- PCModel GetModel ( int  modelNo );                    //!< For accessor functions see http://www.ebi.ac.uk/~keb/cldoc/object/cl_obj_surf.html
    //-- int AddModel ( PCModel  model );                      //!< For editing functions see http://www.ebi.ac.uk/~keb/cldoc/object/cl_obj_edit.html
  };

  //! MMDB atom list class
  /*! This class is used to convert an MMDB PPCAtom to a Clipper Atom_list.
   It is a trivial derivation of a clipper::Atom_list, and may be used
   wherever an Atom_list is required. */
  class CLIPPER_IMEX MMDBAtom_list : public Atom_list
  {
  public:
    //! constructor: from PPCAtom
    MMDBAtom_list( const mmdb::PPCAtom ppcatom, const int natom );
  };


} // namespace clipper

#endif
