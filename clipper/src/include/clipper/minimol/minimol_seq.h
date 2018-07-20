/*! \file minimol_seq.h
  Header file for atomic sequence types */

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


#ifndef CLIPPER_MINIMOL_SEQ
#define CLIPPER_MINIMOL_SEQ


#include "minimol.h"
#include "../imex.h"

namespace clipper {

  //! Polymer sequence object
  /*! The polymer sequence object represents the named sequence of a
    single chain. */
  class CLIPPER_IMEX MPolymerSequence
  {
  public:
    MPolymerSequence() {}  //!< null constructor
    const String& id() const { return id_; }  //!< get sequence ID
    void set_id( const String& s );           //!< set sequence ID
    const String& sequence() const { return seq_; }  //!< get sequence ID
    void set_sequence( const String& s );            //!< set sequence ID
    static String id_tidy( const String& id ) { return id; }  //!< convert ID to std format
    static bool id_match( const String& id1, const String& id2, const MM::MODE& mode ) { return id1 == id2; }  //!< compare two ids
  private:
    String id_;
    String seq_;
  };


  //! Molecule sequence object
  /*! The molecule sequence object is a list of polymer sequence
    objects representing the named sequences of all the chains in a
    molecule. */
  class CLIPPER_IMEX MMoleculeSequence
  {
  public:
    //! number of polymer sequences in model
    int size() const { return children.size(); }
    //! get polymer sequence
    const MPolymerSequence& operator[] ( const int& i ) const { return children[i]; }
    //! set polymer sequence
    MPolymerSequence& operator[] ( const int& i ) { return children[i]; }
    //! get polymer sequence by id
    const MPolymerSequence& find( const String& n, const MM::MODE mode=MM::UNIQUE ) const;
    //! set polymer sequence by id
    MPolymerSequence& find( const String& n, const MM::MODE mode=MM::UNIQUE );
    //! lookup polymer sequence by id
    int lookup( const String& str, const MM::MODE& mode ) const;
    void insert( const MPolymerSequence& add, int pos=-1 );  //!< add polymer sequence

    bool is_null() const { return (size()==0); }  //!< test for null model
  private:
    typedef MPolymerSequence CHILDTYPE;
    std::vector<CHILDTYPE> children;
  };


  //! Sequence alignment obeject
  /*! Provides methods to find an optimal alignment between two sequences. */
  class CLIPPER_IMEX MSequenceAlign {
  public:
    enum TYPE { GLOBAL, LOCAL };
    MSequenceAlign( TYPE type = GLOBAL, ftype match_score = 1.0, ftype miss_score = -0.5, ftype gap_score = -1.0 ) : type_(type), scrmat(match_score), scrmis(miss_score), scrgap(gap_score) {}
    std::pair<std::vector<int>,std::vector<int> > operator() ( const String& seq1, const String& seq2 ) const;
  private:
    TYPE type_;
    ftype32 scrmat, scrmis, scrgap;
  };

} // namespace clipper

#endif
