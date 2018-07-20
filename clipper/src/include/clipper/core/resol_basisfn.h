/*! \file lib/resol_basisfn.h
    Header file for resolution basis function 
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


#ifndef CLIPPER_RESOL_BASISFN
#define CLIPPER_RESOL_BASISFN

#include "resol_fn.h"
#include "../imex.h"

namespace clipper {


  //! Resolution ordinal gernerator
  /*! This class is a helper class for functions which need to divide
    reflections up by resolution whilst guaranteeing a certain
    distribution of number of reflections per range. It takes a list
    of reflections, one at a time, and calculates a function to get
    the approximate ordinal number of a reflection in a list sorted by
    resolution.
  */
  class CLIPPER_IMEX Resolution_ordinal : public Generic_ordinal
  {
  public:
    //! initialiser: takes an HKL_info and uses all reflections.
    void init( const HKL_info& hklinfo, const ftype& power );
    //! initialiser: takes an HKL_data & uses non-missing reflections.
    void init( const HKL_data_base& hkldata, const ftype& power );
    //! initialiser: takes an HKL_data + Cell & uses non-missing reflections.
    void init( const HKL_data_base& hkldata, const Cell& cell, const ftype& power );
  };


  //! simple binning basis function
  /*! This class bins reflections on the basis of resolution, i.e. it
    generates a resolution function from spherical shells. */
  class CLIPPER_IMEX BasisFn_binner : public BasisFn_base
  {
  public:
    //! constructor: include whole reflection list in histogram
    BasisFn_binner( const HKL_info& hklinfo, const int& nbins_, const ftype power = 1.0 ) : BasisFn_base( nbins_ ) { s_ord.init( hklinfo, power ); }
    //! constructor: include only non-missing reflections in histogram
    BasisFn_binner( const HKL_data_base& hkldata, const int& nbins_, const ftype power = 1.0  ) : BasisFn_base( nbins_ ) { s_ord.init( hkldata, hkldata.base_cell(), power ); }
    //! the value of the resolution function (override for speed)
    ftype f_s( const ftype& s, const std::vector<ftype>& params ) const;
    //! the derivative of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv_s( const ftype& s, const std::vector<ftype>& params ) const;
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return LINEAR; }
    //! number of non-zero diagonals in the upper triangle of the curvatures
    int num_diagonals() const { return 1; }
    //! the value of the resolution function (override for speed)
    ftype f( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const { return f_s( hkl.invresolsq( cell ), params ); }
    //! the derivatives of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const { return fderiv_s( hkl.invresolsq( cell ), params ); }
  private:
    Resolution_ordinal s_ord; //<! resolution ordinal
  };


  //! simple linear basis function
  /*! This class fits a piecewise linear function through reflections
    on the basis of resolution. */
  class CLIPPER_IMEX BasisFn_linear : public BasisFn_base
  {
  public:
    //! constructor: include whole reflection list in histogram
    BasisFn_linear( const HKL_info& hklinfo, const int& nbins_, const ftype power = 1.0 ) : BasisFn_base( nbins_ ) { s_ord.init( hklinfo, power ); }
    //! constructor: include only non-missing reflections in histogram
    BasisFn_linear( const HKL_data_base& hkldata, const int& nbins_, const ftype power = 1.0  ) : BasisFn_base( nbins_ ) { s_ord.init( hkldata, hkldata.base_cell(), power ); }
    //! the value of the resolution function (override for speed)
    ftype f_s( const ftype& s, const std::vector<ftype>& params ) const;
    //! the derivative of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv_s( const ftype& s, const std::vector<ftype>& params ) const;
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return LINEAR; }
    //! number of non-zero diagonals in the upper triangle of the curvatures
    int num_diagonals() const { return 2; }
    //! the value of the resolution function (override for speed)
    ftype f( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const { return f_s( hkl.invresolsq( cell ), params ); }
    //! the derivatives of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const { return fderiv_s( hkl.invresolsq( cell ), params ); }
  private:
    Resolution_ordinal s_ord; //<! resolution ordinal
  };


  //! simple smooth basis function
  /*! This class fits a Bspline through reflections on the basis of
    resolution. */
  class CLIPPER_IMEX BasisFn_spline : public BasisFn_base
  {
  public:
    //! constructor: include whole reflection list in histogram
    BasisFn_spline( const HKL_info& hklinfo, const int& nbins_, const ftype power = 1.0 ) : BasisFn_base( nbins_ ) { s_ord.init( hklinfo, power ); }
    //! constructor: include only non-missing reflections in histogram
    BasisFn_spline( const HKL_data_base& hkldata, const int& nbins_, const ftype power = 1.0  ) : BasisFn_base( nbins_ ) { s_ord.init( hkldata, hkldata.base_cell(), power ); }
    //! the value of the resolution function (override for speed)
    ftype f_s( const ftype& s, const std::vector<ftype>& params ) const;
    //! the derivative of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv_s( const ftype& s, const std::vector<ftype>& params ) const;
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return LINEAR; }
    //! number of non-zero diagonals in the upper triangle of the curvatures
    int num_diagonals() const { return 3; }
    //! the value of the resolution function (override for speed)
    ftype f( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const { return f_s( hkl.invresolsq( cell ), params ); }
    //! the derivatives of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const { return fderiv_s( hkl.invresolsq( cell ), params ); }
  private:
    Resolution_ordinal s_ord; //<! resolution ordinal
  };


  //! simple Gaussian basis function
  /*! This class provides a Gaussian basis function. */
  class CLIPPER_IMEX BasisFn_gaussian : public BasisFn_base
  {
  public:
    //! constructor:
    BasisFn_gaussian() : BasisFn_base( 2 ) {}
    //! the value of the resolution function
    //ftype f_s( const ftype& s, const std::vector<ftype>& params ) const;
    //! the derivatives of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv_s( const ftype& s, const std::vector<ftype>& params ) const;
    //! the value of the resolution function (override for speed)
    //ftype f( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const { return f_s( hkl.invresolsq( cell ), params ); }
    //! the derivatives of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const { return fderiv_s( hkl.invresolsq( cell ), params ); }
    //! return the scale factor corresponding to the Gaussian parameters
    ftype scale( const std::vector<ftype>& params ) const;
    //! return the isotropic U corresponding to the Gaussian parameters
    ftype u_iso( const std::vector<ftype>& params ) const;
  };


  //! simple anisotropic Gaussian basis function
  /*! This class provides a anisotropic Gaussian basis function. */
  class CLIPPER_IMEX BasisFn_aniso_gaussian : public BasisFn_base
  {
  public:
    //! constructor:
    BasisFn_aniso_gaussian() : BasisFn_base( 7 ) {}
    //! the derivatives of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv_coord( const Coord_reci_orth& xs, const std::vector<ftype>& params ) const;
    //! the value of the resolution function (override for speed)
    //ftype f( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const { return f_coord( hkl.coord_reci_orth( cell ), params ); }
    //! the derivatives of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const { return fderiv_coord( hkl.coord_reci_orth( cell ), params ); }
    //! return the scale factor corresponding to the Gaussian parameters
    ftype scale( const std::vector<ftype>& params ) const;
    //! return the anisotropic U corresponding to the Gaussian parameters
    U_aniso_orth u_aniso_orth( const std::vector<ftype>& params ) const;
  };


  //! simple log Gaussian basis function
  /*! This class provides a Log Gaussian basis function.  i.e. a
   quadratic function of resolution.
   Use this in conjunction with a Log-target function to get a fast
   estimate to a Gaussian fit. */
  class CLIPPER_IMEX BasisFn_log_gaussian : public BasisFn_base
  {
  public:
    //! constructor:
    BasisFn_log_gaussian() : BasisFn_base( 2 ) {}
    //! the value of the resolution function
    //ftype f_s( const ftype& s, const std::vector<ftype>& params ) const;
    //! the derivatives of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv_s( const ftype& s, const std::vector<ftype>& params ) const;
    //! the value of the resolution function (override for speed)
    //ftype f( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const { return f_s( hkl.invresolsq( cell ), params ); }
    //! the derivatives of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv( const HKL& hkl, const Cell& cell, const
std::vector<ftype>& params ) const { return fderiv_s( hkl.invresolsq( cell ), params ); }
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return LINEAR; }
    //! return the scale factor corresponding to the Gaussian parameters
    ftype scale( const std::vector<ftype>& params ) const;
    //! return the isotropic U corresponding to the Gaussian parameters
    ftype u_iso( const std::vector<ftype>& params ) const;
  };


  //! simple anisotropic Gaussian basis function
  /*! This class provides a anisotropic Gaussian basis function.
   i.e. a general quadratic function of resolution.
   Use this in conjunction with a Log-target function to get a fast
   estimate to a Gaussian fit. */
  class CLIPPER_IMEX BasisFn_log_aniso_gaussian : public BasisFn_base
  {
  public:
    //! constructor:
    BasisFn_log_aniso_gaussian() : BasisFn_base( 7 ) {}
    //! the derivatives of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv_coord( const Coord_reci_orth& xs, const std::vector<ftype>& params ) const;
    //! the value of the resolution function (override for speed)
    //ftype f( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const { return f_coord( hkl.coord_reci_orth( cell ), params ); }
    //! the derivatives of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv( const HKL& hkl, const Cell& cell, const
std::vector<ftype>& params ) const { return fderiv_coord( hkl.coord_reci_orth( cell ), params ); }
    //! the type of the function: optionally used to improve convergence
    FNtype type() const { return LINEAR; }
    //! return the scale factor corresponding to the Gaussian parameters
    ftype scale( const std::vector<ftype>& params ) const;
    //! return the anisotropic U corresponding to the Gaussian parameters
    U_aniso_orth u_aniso_orth( const std::vector<ftype>& params ) const;
  };


  //! simple Expcubic basis function
  /*! This class provides a Expcubic basis function. */
  class CLIPPER_IMEX BasisFn_expcubic : public BasisFn_base
  {
  public:
    //! constructor
    BasisFn_expcubic() : BasisFn_base( 4 ) {}
    //! the derivatives of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv_s( const ftype& s, const std::vector<ftype>& params ) const;
    //! the value of the resolution function (override for speed)
    //ftype f( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const { return f_s( hkl.invresolsq( cell ), params ); }
    //! the derivatives of the resolution function w.r.t. the parameters
    const BasisFn_base::Fderiv& fderiv( const HKL& hkl, const Cell& cell, const std::vector<ftype>& params ) const { return fderiv_s( hkl.invresolsq( cell ), params ); }
  };


} // namespace clipper

#endif
