/*! \file ramachandran.h
    Ramachandran Plot class
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


#ifndef CLIPPER_RAMACHANDRAN
#define CLIPPER_RAMACHANDRAN


#include "clipper_types.h"
#include "../imex.h"

namespace clipper {

  //! 2-d angular probability distibution class
  /*! Base for Ramachandran class (and other similar classes, such as
    a pseudo-ramachandran plot or the JPD of two phases ). */
  class CLIPPER_IMEX Prob_phi_2d
  {
  public:
    //! initialise: with sampling
    void init( const int& size );
    //! accumulate new table of samples to probability
    void accumulate( const ftype32 table[] );
    //! accumulate new sample to probability
    void accumulate( const ftype& phi1, const ftype& phi2, ftype wgt = 1.0 );
    //! normalise to integrate to 1/(2pi)^2
    void normalise();
    //! get probability for a particular pair of angles
    ftype probability( const ftype& phi1, const ftype& phi2 ) const;
    //! formatted string representation (as C++ code)
    String format() const;
    //! 2d read access
    const ftype& data( const int& i, const int& j ) const
      { return data_[n*i+j]; }
    //! 2d write access
    ftype& data( const int& i, const int& j )
      { return data_[n*i+j]; }
  private:
    int n;  //!< sampling
    std::vector<ftype> data_;
  };


  //! Ramachandran plot class
  /*! This class provides a reference Ramachandran plot for Gly, Pro,
    other, and combinations of those types of residues. The source
    data comes from the best residues from the 'top500' best-determined
    structures list of D. C. and J. S. Richardson,
    http://kinemage.biochem.duke.edu/index.html

    The Ramachandran plot is normalised in inverse radians squared,
    so the mean value of a probability is 1/(2 pi)<sup>2</sup>. */
  class CLIPPER_IMEX Ramachandran : private Prob_phi_2d
  {
  public:
    //! enumeration of built-in Ramachandran tables
    enum TYPE { Gly, Pro, NonGlyPro, NonGly, All, Gly5, Pro5, NonGlyPro5, NonGly5, All5 };
    //! null constructor
    Ramachandran() {}
    //! constructor: from standard plot
    Ramachandran( TYPE type );
    //! initialise: from standard plot
    void init( TYPE type );
    //! change threshholds to different values
    void set_thresholds( ftype prob_favored = 0.01,
			 ftype prob_allowed = 0.0005 );
    //! get probability for a particular pair of angles
    ftype probability( const ftype& phi, const ftype& psi ) const
      { return Prob_phi_2d::probability( phi, psi ); }
    //! test if a pair of angles are in the favored region
    bool favored( const ftype& phi, const ftype& psi ) const
      { return ( probability( phi, psi ) > p_favored ); }
    //! test if a pair of angles are in the allowed region
    bool allowed( const ftype& phi, const ftype& psi ) const
      { return ( probability( phi, psi ) > p_allowed ); }
  private:
    ftype p_favored, p_allowed;
  };


}  // namespace clipper


#endif
