/* sfcalc.cpp: Structure factor calculation implementation */
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


#include "sfcalc.h"

#include "../core/atomsf.h"


namespace clipper {


template<class T> bool SFcalc_iso_sum<T>::operator() ( HKL_data<datatypes::F_phi<T> >& fphidata, const Atom_list& atoms ) const
{
  const HKL_info& hkls   = fphidata.base_hkl_info();
  const Cell& cell       = fphidata.base_cell();
  const Spacegroup& spgr = hkls.spacegroup();

  fphidata = datatypes::F_phi<T>( std::complex<T>( 0.0, 0.0 ) );

  HKL_info::HKL_reference_index ih;
  for ( int i = 0; i < atoms.size(); i++ ) if ( !atoms[i].is_null() ) {
    AtomShapeFn sf( atoms[i].coord_orth(), atoms[i].element(),
		    atoms[i].u_iso(), atoms[i].occupancy() );
    for ( int j = 0; j < spgr.num_symops(); j++ ) {
      Coord_frac uvw =
        atoms[i].coord_orth().coord_frac( cell ).transform( spgr.symop(j) );
      for ( ih = hkls.first(); !ih.last(); ih.next() ) {
        T phi = Util::twopi() * (ih.hkl().coord_reci_frac()*uvw);
        fphidata[ih] = std::complex<T>(fphidata[ih]) +
	  T( sf.f(ih.invresolsq()) ) *
          std::complex<T>( cos(phi), sin(phi) );
      }
    }
  }
  return true;
}


template<class T> bool SFcalc_aniso_sum<T>::operator() ( HKL_data<datatypes::F_phi<T> >& fphidata, const Atom_list& atoms ) const
{
  const HKL_info& hkls   = fphidata.base_hkl_info();
  const Cell& cell       = fphidata.base_cell();
  const Spacegroup& spgr = hkls.spacegroup();

  fphidata = datatypes::F_phi<T>( std::complex<T>( 0.0, 0.0 ) );

  HKL_info::HKL_reference_index ih;
  for ( int i = 0; i < atoms.size(); i++ ) if ( !atoms[i].is_null() ) {
    for ( int j = 0; j < spgr.num_symops(); j++ ) {
      Atom atom( atoms[i] );
      atom.transform( spgr.symop(j).rtop_orth( cell ) );
      AtomShapeFn sf( atom );
      Coord_frac uvw =	atom.coord_orth().coord_frac( cell );
      for ( ih = hkls.first(); !ih.last(); ih.next() ) {
	T phi = Util::twopi() * (ih.hkl().coord_reci_frac()*uvw);
	fphidata[ih] = std::complex<T>(fphidata[ih]) +
	  T(sf.f(ih.hkl().coord_reci_orth(cell))) *
	  std::complex<T>( cos(phi), sin(phi) );
      }
    }
  }
  return true;
}


template<class T> bool SFcalc_iso_fft<T>::operator() ( HKL_data<datatypes::F_phi<T> >& fphidata, const Atom_list& atoms ) const
{
  const HKL_info& hkls   = fphidata.base_hkl_info();
  const Cell& cell       = fphidata.base_cell();
  const Spacegroup& spgr = hkls.spacegroup();
  const Grid_sampling grid( spgr, cell, hkls.resolution(), rate_ );
  Xmap<ftype32> xmap( spgr, cell, grid );

  Coord_frac uvw, duvw;
  Coord_grid g0, g1;
  Grid_range gd( cell, grid, radius_ );
  Xmap<ftype32>::Map_reference_coord i0, iu, iv, iw;
  for ( int i = 0; i < atoms.size(); i++ ) if ( !atoms[i].is_null() ) {
    AtomShapeFn sf( atoms[i].coord_orth(), atoms[i].element(),
		    atoms[i].u_iso() + uadd_, atoms[i].occupancy() );
    uvw = atoms[i].coord_orth().coord_frac( cell );
    g0 = uvw.coord_grid( grid ) + gd.min();
    g1 = uvw.coord_grid( grid ) + gd.max();
    i0 = Xmap<ftype32>::Map_reference_coord( xmap, g0 );
    for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
      for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
	for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() )
	  xmap[iw] += sf.rho( iw.coord_orth() );
  }

  for ( Xmap<ftype32>::Map_reference_index ix = xmap.first();
	!ix.last(); ix.next() )
    xmap[ix] *= xmap.multiplicity( ix.coord() );

  xmap.fft_to( fphidata );
  if ( uadd_ != 0.0 ) {
    ftype u = Util::twopi2()*uadd_;
    for ( HKL_info::HKL_reference_index ih = fphidata.first_data();
	  !ih.last(); fphidata.next_data( ih ) )
      fphidata[ih].scale( exp( u * ih.invresolsq() ) );
  }

  return true;
}


template<class T> bool SFcalc_aniso_fft<T>::operator() ( HKL_data<datatypes::F_phi<T> >& fphidata, const Atom_list& atoms ) const
{
  const HKL_info& hkls   = fphidata.base_hkl_info();
  const Cell& cell       = fphidata.base_cell();
  const Spacegroup& spgr = hkls.spacegroup();
  const Grid_sampling grid( spgr, cell, hkls.resolution(), rate_ );
  Xmap<ftype32> xmap( spgr, cell, grid );

  U_aniso_orth uadd( uadd_ );
  Coord_frac uvw, duvw;
  Coord_grid g0, g1;
  Grid_range gd( cell, grid, radius_ );
  Xmap<ftype32>::Map_reference_coord i0, iu, iv, iw;
  for ( int i = 0; i < atoms.size(); i++ ) if ( !atoms[i].is_null() ) {
    U_aniso_orth u( atoms[i].u_aniso_orth() );
    if ( u.is_null() ) u = U_aniso_orth( atoms[i].u_iso() );
    AtomShapeFn sf( atoms[i].coord_orth(), atoms[i].element(),
		    u + uadd, atoms[i].occupancy() );
    uvw = atoms[i].coord_orth().coord_frac( cell );
    g0 = uvw.coord_grid( grid ) + gd.min();
    g1 = uvw.coord_grid( grid ) + gd.max();
    i0 = Xmap<ftype32>::Map_reference_coord( xmap, g0 );
    for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
      for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
	for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() )
	  xmap[iw] += sf.rho( iw.coord_orth() );
  }

  for ( Xmap<ftype32>::Map_reference_index ix = xmap.first();
	!ix.last(); ix.next() )
    xmap[ix] *= xmap.multiplicity( ix.coord() );

  xmap.fft_to( fphidata );
  if ( uadd_ != 0.0 ) {
    ftype u = Util::twopi2()*uadd_;
    for ( HKL_info::HKL_reference_index ih = fphidata.first_data();
	  !ih.last(); fphidata.next_data( ih ) )
      fphidata[ih].scale( exp( u * ih.invresolsq() ) );
  }

  return true;
}


// compile templates

template class CLIPPER_IMEX SFcalc_iso_sum<ftype32>;
template class CLIPPER_IMEX SFcalc_aniso_sum<ftype32>;
template class CLIPPER_IMEX SFcalc_iso_fft<ftype32>;
template class CLIPPER_IMEX SFcalc_aniso_fft<ftype32>;

template class CLIPPER_IMEX SFcalc_iso_sum<ftype64>;
template class CLIPPER_IMEX SFcalc_aniso_sum<ftype64>;
template class CLIPPER_IMEX SFcalc_iso_fft<ftype64>;
template class CLIPPER_IMEX SFcalc_aniso_fft<ftype64>;



} // namespace clipper
