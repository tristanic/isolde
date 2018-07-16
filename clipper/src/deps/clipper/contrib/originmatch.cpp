/* originmatch.cpp: Origin matching implementation */
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


#include "originmatch.h"


namespace clipper {

template<class T> bool OriginMatch<T>::operator() ( bool& invert, clipper::Coord_frac& shift, const HKL_data<datatypes::F_phi<T> >& fphi1, const HKL_data<datatypes::F_phi<T> >& fphi2 ) const
{
  const Spacegroup& oldspgr = fphi1.base_hkl_info().spacegroup();
  const Cell&       oldcell = fphi1.base_hkl_info().cell();
  const Resolution& oldreso = fphi1.base_hkl_info().resolution();

  // calc point group
  const Spgr_descr spgrdscr( oldspgr.generator_ops().pgrp_ops() );
  const Spacegroup newspgr( spgrdscr );

  // new cell and resolution
  const Cell       newcell( oldcell );
  const Resolution newreso( Util::max( limit, oldreso.limit() ) );

  // new reflection list
  HKL_info hkls( newspgr, newcell, newreso, true );
  HKL_data<datatypes::F_phi<ftype32> > fphi(hkls), fphiinv(hkls);
  HKL_info::HKL_reference_index ih;
  HKL hkl;
  datatypes::F_phi<T> fp1, fp2;
  for ( ih = hkls.first(); !ih.last(); ih.next() ) {
    hkl = ih.hkl();
    fp1 = fphi1[hkl];
    fp2 = fphi2[hkl];
    if ( !fp1.missing() && !fp2.missing() ) {
      fphi[ih].f() = fphiinv[ih].f() = fp1.f() * fp2.f();
      fphi[ih].phi()    = fp1.phi() - fp2.phi();
      fphiinv[ih].phi() = fp1.phi() + fp2.phi();
    } else {
      fphi[ih].f() = fphiinv[ih].f() = fphi[ih].phi() = fphiinv[ih].phi() = 0.0;
    }
  }

  // new map
  Grid_sampling oldgrid( newspgr, newcell, newreso );
  Grid_sampling newgrid( ((oldgrid.nu()+11)/12)*12, ((oldgrid.nv()+11)/12)*12,
			 ((oldgrid.nw()+11)/12)*12 );
  Xmap<ftype32> xmap( newspgr, newcell, newgrid );
  Xmap<ftype32>::Map_reference_index ix;

  // search for highest peak
  T mapmax = 0.0;
  invert = false;
  shift  = Coord_frac( 0.0, 0.0, 0.0 );
  {
    xmap.fft_from( fphi );
    for ( ix = xmap.first(); !ix.last(); ix.next() )
      if ( xmap[ix] > mapmax ) {
	mapmax = xmap[ix];
	shift  = ix.coord().coord_frac( newgrid );
	invert = false;
      }
  }
  if ( oldspgr.invariant_under_change_of_hand() ) {
    xmap.fft_from( fphiinv );
    for ( ix = xmap.first(); !ix.last(); ix.next() )
      if ( xmap[ix] > mapmax ) {
	mapmax = xmap[ix];
	shift  = ix.coord().coord_frac( newgrid );
	invert = true;
      }
  }

  // check against symmetry restrictions

  // return the shift
  return true;
}


// compile templates

template class CLIPPER_IMEX OriginMatch<ftype32>;

template class CLIPPER_IMEX OriginMatch<ftype64>;


} // namespace clipper
