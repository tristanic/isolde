/* minimol_utils.cpp: minimol utils */
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


#include "minimol_utils.h"

#include <algorithm>


namespace clipper {


/*! Set up a non-bonding type search database to allow fast access to
  a list of atoms within a given radius of a given coordinate. Setting
  up the class takes a litle (but not a great deal) of time, returning
  coordinates is very fast indeed. You can also specify a radius to
  the constructor, which determines memory use as (24*cell
  volume/rad^3) bytes. For small radii, set this to the search radius
  or larger. For large radii, a smaller value may improve
  performance.
  \param mol The MiniMol object containing the atoms.
  \param rad A radius for sampling crystal space.
*/
MAtomNonBond::MAtomNonBond( const clipper::MiniMol& mol, double rad )
{
  // store params
  mol_ = &mol;
  rad_ = rad;
  cell = mol.cell();
  spgr = mol.spacegroup();
  // pick a grid for the cell
  grid = Grid_sampling( Util::intc(1.0/(rad_*cell.a_star())), Util::intc(1.0/(rad_*cell.b_star())), Util::intc(1.0/(rad_*cell.c_star())));
  // now create a list of atoms in the cell
  int gi;
  Coord_frac cf;
  std::vector<std::pair<int,MAtomIndexSymmetry> > iatoms;
  for ( int p = 0; p < mol.size(); p++ )
    for ( int m = 0; m < mol[p].size(); m++ )
      for ( int a = 0; a < mol[p][m].size(); a++ )
	if ( !mol[p][m][a].is_null() ) {
	  for ( int s = 0; s < spgr.num_symops(); s++ ) {
	    cf = spgr.symop(s) * mol[p][m][a].coord_orth().coord_frac(cell);
	    gi = cf.coord_grid(grid).unit(grid).index(grid);
	    iatoms.push_back( std::pair<int,MAtomIndexSymmetry>( gi, MAtomIndexSymmetry( p, m, a, s ) ) );
	  }
	}
  // sort them by grid coordinate index
  std::sort( iatoms.begin(), iatoms.end() );
  // turn the index into a list of atom indices, indexed by a list of cells
  lookup.resize( grid.size()+1, iatoms.size() );
  atoms.resize( iatoms.size() );
  // fill out the atom index
  for ( int i = 0; i < iatoms.size(); i++ ) {
    atoms[i] = iatoms[i].second;
    if ( lookup[iatoms[i].first] == iatoms.size() )
      lookup[iatoms[i].first] = i;
  }
  // and fill in any missing terms
  for ( int i = lookup.size()-1; i > 0; i-- )
    if ( lookup[i-1] > lookup[i] ) lookup[i-1] = lookup[i];
}


/*! Return a list of atom indices in the MiniMol objects, and symop
  numbers, of atoms within the previously specified radius of the
  given coordinate. The function always returns all atoms within the
  specified radius, but also some beyond it.
  This function is very fast.
  \param co The coordinate about which to search.
  \param rad The radius to search.
  \return A vector of atom indices near the given coordinate.
*/
std::vector<MAtomIndexSymmetry> MAtomNonBond::atoms_near( const clipper::Coord_orth& co, double rad ) const
{
  std::vector<MAtomIndexSymmetry> result;
  Coord_grid cg = co.coord_frac(cell).coord_grid(grid);
  Coord_grid c;
  int d = Util::intc( rad/rad_ - 1.0e-4 );
  int gi;
  for ( int u = cg.u()-d; u <= cg.u()+d; u++ )
    for ( int v = cg.v()-d; v <= cg.v()+d; v++ )
      for ( int w = cg.w()-d; w <= cg.w()+d; w++ ) {
	gi = Coord_grid(u,v,w).unit(grid).index(grid);
	for ( int i = lookup[gi]; i < lookup[gi+1]; i++ )
	  result.push_back( atoms[i] );
      }
  return result;
}


/*! Return a list of atom indices in the MiniMol objects, and symop
  numbers, of atoms within the previously specified radius of the
  given coordinate. The function only returns atoms within the
  specified radius.
  This function is very fast.
  \param co The coordinate about which to search.
  \param rad The radius to search.
  \return A vector of atom indices near the given coordinate.
*/
std::vector<MAtomIndexSymmetry> MAtomNonBond::operator() ( const clipper::Coord_orth& co, double rad ) const
{
  std::vector<MAtomIndexSymmetry> result;
  std::vector<MAtomIndexSymmetry> atoms = atoms_near( co, rad );
  const MiniMol& mol = *mol_;
  Coord_frac f1, f2;
  f1 = co.coord_frac(cell);
  for ( int i = 0; i < atoms.size(); i++ ) {
    f2 = spgr.symop(atoms[i].symmetry()) * mol[atoms[i].polymer()][atoms[i].monomer()][atoms[i].atom()].coord_orth().coord_frac(cell);
    double l2 = ( f2.lattice_copy_near( f1 ) - f1 ).lengthsq( cell );
    if ( l2 < rad*rad ) result.push_back( atoms[i] );
  }
  return result;
}


void MAtomNonBond::debug() const
{
  std::cout << grid.size() << "\t" << atoms.size() << "\n";
  for ( int i = 0; i < lookup.size(); i++ )
    std::cout << i << "\t" << lookup[i] << "\n";
}


} // namespace clipper
