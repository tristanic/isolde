/*! \file atomsf.h
  Header file for atomic scattering factors
  \ingroup g_fftmapbspl
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


#ifndef CLIPPER_ATOMSF
#define CLIPPER_ATOMSF


#include "coords.h"
#include "../imex.h"

namespace clipper
{

  //! Atomic shape function object
  /*! The atomic scattering factor object is instantiated for each
    atom in turn, giving the atom parameters: position, element,
    occupancy and the isotropic or anisotropic U-value. (See
    clipper::Util for conversion from B-factors.). The methods of the
    class may then be called to return the scattering in reciprocal
    space or density in real space using either isotropic or
    anistropic models as required.

    If the atom only has an isotropic U, the faster isotropic methods
    will be used where available.

    This implementation uses the coefficients from Waasmaier & Kirfel
    (1995), Acta Cryst. A51, 416-431. The source data can be found at:
    ftp://wrzx02.rz.uni-wuerzburg.de/pub/local/Crystallography/sfac.dat
  */
  class CLIPPER_IMEX AtomShapeFn
  {
  public:
    enum TYPE { X, Y, Z, Uiso, Occ, U11, U22, U33, U12, U13, U23 };
    //! null constructor
    AtomShapeFn() {}
    //! constructor: from atom object
    AtomShapeFn( const Atom& atom );
    //! constructor: from coord, element, isotropic U, occupancy
    AtomShapeFn( const Coord_orth& xyz, const String& element, const ftype u_iso = 0.0, const ftype occ = 1.0 );
    //! constructor: from coord, element, anisotropic U, occupancy
    AtomShapeFn( const Coord_orth& xyz, const String& element, const U_aniso_orth& u_aniso, const ftype occ = 1.0 );
    //! initialiser:  from atom object
    void init( const Atom& atom );
    //! initialiser: from coord, element, isotropic U, occupancy
    void init( const Coord_orth& xyz, const String& element, const ftype u_iso = 0.0, const ftype occ = 1.0 );
    //! initialiser: from coord, element, anisotropic U, occupancy
    void init( const Coord_orth& xyz, const String& element, const U_aniso_orth& u_aniso, const ftype occ = 1.0 );

    //! return scattering factor as a function of reflection posn
    ftype f( const Coord_reci_orth& rfl ) const;
    //! return electron density as a function of coordinate
    ftype rho( const Coord_orth& xyz ) const;

    //! return Agarwal density gradients as a function of coordinate
    bool rho_grad( const Coord_orth& xyz, ftype& rho, std::vector<ftype>& grad ) const;
    //! return Agarwal density gradient/curvature as a function of coordinate
    bool rho_curv( const Coord_orth& xyz, ftype& rho, std::vector<ftype>& grad, Matrix<ftype>& curv ) const;

    //! \deprecated return Agarwal density gradients as a function of coordinate
    bool rho_grad( const Coord_orth& xyz, std::vector<ftype>& grad ) const;

    //! return (isotropic) scattering factor as a function of resolution
    ftype f( const ftype& invresolsq ) const;
    //! return (isotropic) electron density as a function of radius
    ftype rho( const ftype& rsq ) const;

    //! define parameters for Agarwal gradient/curvature calcs
    std::vector<TYPE>& agarwal_params() { return params; }
  private:
    //! look up atom coeffs
    void init( const String& element, const ftype& u_iso );
    // members
    Coord_orth coord_;
    U_aniso_orth u_aniso_;
    ftype u_iso_, occ_;
    ftype a[6],  b[6];                //!< atom coeffs
    ftype aw[6], bw[6];               //!< intermediate results (iso)
    std::vector<Mat33sym<> > uaninv;  //!< intermediate results (aniso)
    bool is_iso;
    std::vector<TYPE> params;
  };


  //! Atomic scattering factor object
  /*! \deprecated This class has been replaced by AtomShapeFn, which
    is smaller, faster, and more capable. This class is now a wrapper
    for that class. */
  class CLIPPER_IMEX AtomSF : private AtomShapeFn
  {
  public:
    AtomSF( const String& type, const ftype u_iso = 0.0, const ftype occ = 1.0 );
    AtomSF( const String& type, const U_aniso_orth& u_aniso, const ftype occ = 1.0 );
    void init( const String& type, const ftype u_iso = 0.0, const ftype occ = 1.0 );
    void init( const String& type, const U_aniso_orth& u_aniso, const ftype occ = 1.0 );
    ftype f_iso( const ftype& s ) const;
    ftype f_aniso( const Coord_reci_orth& rfl ) const;
    ftype rho_iso( const ftype& d2 ) const;
    ftype rho_aniso( const Coord_orth& uvw ) const;
  };


} // namespace clipper

#endif
