/* sfweight.cpp: structure factor calculation vs observed implementation */
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


#include "sfcalc_obs.h"
#include "edcalc.h"
#include "../core/hkl_compute.h"
#include "../core/resol_targetfn.h"


namespace clipper {


template<class T> bool SFcalc_obs_bulk<T>::operator() ( HKL_data<datatypes::F_phi<T> >& fphi, const HKL_data<datatypes::F_sigF<T> >& fsig, const Atom_list& atoms )
{
  // set U value constants
  double u_atom = Util::b2u( 20.0 );
  double u_mask = Util::b2u( 50.0 );

  // increase the U values
  Atom_list atomu = atoms;
  U_aniso_orth uadd( u_atom ), u;
  for ( int i = 0; i < atomu.size(); i++ ) if ( !atomu[i].is_null() ) {
    u = atomu[i].u_aniso_orth();
    if ( u.is_null() ) u = U_aniso_orth( atomu[i].u_iso() );
    atomu[i].set_u_aniso_orth( u + uadd );
  }

  // now make the map for ed calcs
  const HKL_info&   hkls = fsig.base_hkl_info();
  const Spacegroup& spgr = hkls.spacegroup();
  const Cell&       cell = fsig.base_cell();
  HKL_data<datatypes::F_phi<T> >
    fphi_atom( hkls, cell ), fphi_mask( hkls, cell );
  const Grid_sampling grid( spgr, cell, hkls.resolution() );
  Xmap<float> xmap( spgr, cell, grid );

  // do ed calc from atomu
  EDcalc_aniso<ftype32> edcalc;
  edcalc( xmap, atomu );
  xmap.fft_to( fphi_atom );
  fphi_atom.compute( fphi_atom, datatypes::Compute_scale_u_iso<datatypes::F_phi<T> >( 1.0, u_atom ) );

  // do density calc from mask
  EDcalc_mask<ftype32> emcalc;
  emcalc( xmap, atomu );
  for ( Xmap<ftype32>::Map_reference_index ix = xmap.first();
	!ix.last(); ix.next() )
    xmap[ix] = 1.0 - xmap[ix];
  xmap.fft_to( fphi_mask );
  fphi_mask.compute( fphi_mask, datatypes::Compute_scale_u_iso<datatypes::F_phi<T> >( 1.0, -u_mask ) );

  // try some different scale factors
  std::vector<double> params( nparams, 1.0 );
  BasisFn_spline basisfn( hkls, nparams, 1.0 );
  TargetFn_scaleF1F2<datatypes::F_phi<T>,datatypes::F_sigF<T> > targetfn( fphi, fsig );
  T x1 = 0.35, dx = 0.35, x;
  ftype y[3] = { 0.0, 0.0, 0.0 };
  for ( int i = 0; i < 8; i++ ) {
    // take 3 samples of function
    for ( int d = -1; d <= 1; d++ ) if ( y[d+1] == 0.0 ) {
      HKL_data<data32::F_phi>::HKL_reference_index ih;
      x = x1 + T(d)*dx;
      for ( ih = fphi.first(); !ih.last(); ih.next() )
	fphi[ih] = std::complex<T>(fphi_atom[ih]) +
	       x * std::complex<T>(fphi_mask[ih]);
      ResolutionFn rfn( hkls, basisfn, targetfn, params );
      double r = 0.0;
      for ( ih = fsig.first(); !ih.last(); ih.next() )
	if ( !fsig[ih].missing() ) {
	  double eps = ih.hkl_class().epsilon();
	  r += (2.0/eps) * fabs( sqrt(rfn.f(ih))*fphi[ih].f() - fsig[ih].f() );
	  // r += ( 2.0/eps ) * pow( rfn.f(ih) * pow(fphi[ih].f(),2)/eps - pow(fsig[ih].f(),2)/eps, 2 );
	}
      y[d+1] = r;
      //std::cout << d << "\t" << x << "\t" << r << "\n";
    }
    // find minimum of current 3 samples
    if      ( y[0] < y[1] && y[0] < y[2] ) { y[1] = y[0]; x1 -= dx; }
    else if ( y[2] < y[1] && y[2] < y[0] ) { y[1] = y[2]; x1 += dx; }
    // reduce step and search again
    y[0] = y[2] = 0.0;
    dx /= 2.0;
  }

  // adopt final scale
  for ( HKL_data<data32::F_phi>::HKL_reference_index ih = fphi.first();
	!ih.last(); ih.next() )
    fphi[ih] = std::complex<T>(fphi_atom[ih]) +
          x1 * std::complex<T>(fphi_mask[ih]);

  // store stats
  ftype64 w, s0 = 0.0, s1 = 0.0;
  for ( Xmap<ftype32>::Map_reference_index ix = xmap.first();
	!ix.last(); ix.next() ) {
    w = 1.0/ftype64( xmap.multiplicity( ix.coord() ) );
    s0 += w;
    s1 += w*xmap[ix];
  }
  bulkfrc = s1/s0;
  bulkscl = x1;

  return true;
}


// compile templates

template class CLIPPER_IMEX SFcalc_obs_bulk<ftype32>;

template class CLIPPER_IMEX SFcalc_obs_bulk<ftype64>;


} // namespace clipper
