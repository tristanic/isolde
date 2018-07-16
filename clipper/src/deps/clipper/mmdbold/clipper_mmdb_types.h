/*! \file clipper_mmdb_types.h
  Header file for atomic model types
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


#ifndef CLIPPER_MMDBOLD_TYPES
#define CLIPPER_MMDBOLD_TYPES


#include "../core/coords.h"


namespace clipper {

  // PDB 'sigma' classes (should be depracated)

  //! Standard deviation of orthognal coordinates
  /*! see clipper::Coord_orth
    \note In my view this is a stupid definition, but its in the PDB.
    \ingroup g_mmdb */
  class Sig_Coord_orth : public Vec3<>
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
  class Sig_U_aniso_orth : public Mat33sym<>
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


  // abstract base classes ========================================

  //! Abstract base class for a single atom
  /*! \ingroup g_mmdb */
  class DBAtom_base
  {
  public:
    virtual ~DBAtom_base() {}
    //! get atom name (e.g. CA, CZ1, OH)
    /*! \return The atom type name. */
    virtual String type() const = 0;
    //! get atom element name (e.g. C, O, ZN2+ )
    /*! \return The atom element name. */
    virtual String element() const = 0;
    //! get atom alternate conformation code
    /*! \return The atom alternate conformation code. */
    virtual String altconf() const = 0;
    //! get atom coordinate (in orthogonal Angstroms)
    /*! \return The atom coordinate. (result.is_null() if not set) */
    virtual Coord_orth coord_orth() const = 0;
    //! get atom occupancy
    /*! \return The atom occupancy. (NaN if not set) */
    virtual ftype occupancy() const = 0;
    //! get atom atomic displacement parameter (U iso)
    /*! \return The atom temperature factor. (NaN if not set) */
    virtual ftype u_iso() const = 0;
    //! get atom anisotropic U
    /*! \return The atom U_aniso. (result.is_null() if not set) */
    virtual U_aniso_orth u_aniso_orth() const = 0;
    //! set atom name (e.g. CA, CZ1, OH)
    /*! \param n The atom type name. */
    virtual void set_type( const String& n ) = 0;
    //! set atom element name (e.g. C, O, ZN2+  )
    /*! \param n The atom element name. */
    virtual void set_element( const String& n ) = 0;
    //! set atom alternate conformation code
    /*! \param n The atom alternate conformation code. */
    virtual void set_altconf( const String& n ) = 0;
    //! set void atom coordinate (in orthogonal Angstroms)
    /*! \param v The atom coordinate. (Use null coordinate to clear) */
    virtual void set_coord_orth( const Coord_orth& v ) = 0;
    //! set atom occupancy
    /*! \param v The atom occupancy. (Use NaN to clear) */
    virtual void set_occupancy( const ftype& v ) = 0;
    //! set atom atomic displacement parameter (U iso)
    /*! \param v The atom temperature factor. (Use NaN to clear) */
    virtual void set_u_iso( const ftype& v ) = 0;
    //! set atom anisotropic U
    /*! \param v The atom U_aniso. (Use null U_aniso to clear) */
    virtual void set_u_aniso_orth( const U_aniso_orth& v ) = 0;

    //! apply an RT operator to the atom
    void transform( const RTop_orth rt );
  };

  //! Abstract base class for a single monomer (residue or base)
  /*! \ingroup g_mmdb */
  class DBResidue_base
  {
  public:
    virtual ~DBResidue_base() {}
    //! get monomer name/type (e.g. TYR/ARG, A/T/C/G)
    virtual String type() const = 0;
    //! get monomer sequence number
    virtual int seqnum() const = 0;
    //! get monomer insertion code
    virtual String inscode() const = 0;
    //! set monomer name/type (e.g. TYR/ARG, A/T/C/G)
    virtual void set_type( const String& n ) = 0;
    //! set monomer sequence number
    virtual void set_seqnum( const int& n ) = 0;
    //! set monomer insertion code
    virtual void set_inscode( const String& n ) = 0;
  };

  //! Abstract base class for a single polymer (chain)
  /*! \ingroup g_mmdb */
  class DBChain_base
  {
  public:
    virtual ~DBChain_base() {}
    //! get polymer id (e.g. A, B)
    virtual String id() const = 0;
    //! get polymer id (e.g. A, B)
    virtual void set_id( const String& n ) = 0;
  };

  //! Abstract base class for a model
  /*! \ingroup g_mmdb */
  class DBModel_base
  {
  public:
    virtual ~DBModel_base() {}
    //! get polymer id (e.g. A, B)
    virtual String id() const = 0;
    //! get polymer id (e.g. A, B)
    virtual void set_id( const String& n ) = 0;
  };


  // trivial non-abstract examples ========================================

  //! A single non-database atom
  /*! \ingroup g_mmdb */
  class NDBAtom : public DBAtom_base
  {
  public:
    NDBAtom() {}                      //!< null constuctor
    NDBAtom( const DBAtom_base& a );  //!< constuctor: from Atom base
    String type() const;
    String element() const;
    String altconf() const;
    Coord_orth coord_orth() const;
    ftype occupancy() const;
    ftype u_iso() const;
    U_aniso_orth u_aniso_orth() const;
    void set_type( const String& n );
    void set_element( const String& n );
    void set_altconf( const String& n );
    void set_coord_orth( const Coord_orth& v );
    void set_occupancy( const ftype& v );
    void set_u_iso( const ftype& v );
    void set_u_aniso_orth( const U_aniso_orth& v );
    static NDBAtom null();
  private:
    String type_, element_, altconf_;
    Coord_orth xyz; Sig_Coord_orth sig_xyz;
    ftype occ,  sig_occ;
    ftype u, sig_u;
    U_aniso_orth uij; Sig_U_aniso_orth sig_uij;
  };

  //! A non-database single monomer with no contents (residue or base)
  /*! \ingroup g_mmdb */
  class NDBResidue : public DBResidue_base
  {
  public:
    NDBResidue() {}                         //!< null constuctor
    NDBResidue( const DBResidue_base& r );  //!< constuctor: from Residue base
    String type() const;
    int seqnum() const;
    String inscode() const;
    void set_type( const String& n );
    void set_seqnum( const int& n );
    void set_inscode( const String& n );
  private:
    int seqnum_;
    String type_, inscode_;
  };

  //! A single non-database polymer with no contents (chain)
  /*! \ingroup g_mmdb */
  class NDBChain : public DBChain_base
  {
  public:
    NDBChain() {}                       //!< null constuctor
    NDBChain( const DBChain_base& c );  //!< constuctor: from Chain base
    String id() const;
    void set_id( const String& n );
  private:
    String id_;
  };

  //! A single non-database model with no contents
  /*! \ingroup g_mmdb */
  class NDBModel : public DBModel_base
  {
  public:
    NDBModel() {}                       //!< null constuctor
    NDBModel( const DBModel_base& m );  //!< constuctor: from Model base
    String id() const;
    void set_id( const String& n );
  private:
    String id_;
  };


} // namespace clipper

#endif
