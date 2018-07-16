/* hkl_lookup.cpp: class file for reflection data lookup */
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


#include "hkl_lookup.h"


namespace clipper {


void HKL_lookup::init(const std::vector<HKL>& hkl)
{
  // generate an efficient 3-d lookup table for finding a reflection
  // given its index. We could use hash_map, but this approach is
  // more efficient and uses far less memory
  int i, h, k;

  if ( hkl.size() == 0 ) return;

  // make h lookup table
  for (i=0; i<hkl.size(); i++) {
    const int& h = hkl[i].h();  // find range of h
    h_ptr.min = Util::min( h_ptr.min, h );
    h_ptr.max = Util::max( h_ptr.max, h );
  }
  if ( h_ptr.max >= h_ptr.min )
    h_ptr.p.resize( h_ptr.max - h_ptr.min + 1 );    // now add k lists

  // make k lookup tables
  for (i=0; i<hkl.size(); i++) {
    klookup& k_ptr = h_ptr.p[hkl[i].h() - h_ptr.min];
    const int& k = hkl[i].k();  // find range of k
    k_ptr.min = Util::min( k_ptr.min, k );
    k_ptr.max = Util::max( k_ptr.max, k );
  }
  for (h=h_ptr.min; h<=h_ptr.max; h++) {
    klookup& k_ptr = h_ptr.p[h - h_ptr.min];
    if ( k_ptr.max >= k_ptr.min )
      k_ptr.p.resize( k_ptr.max - k_ptr.min + 1 );  // now add l lists
  }

  // make l lookup tables
  for (i=0; i<hkl.size(); i++) {
    klookup& k_ptr = h_ptr.p[hkl[i].h() - h_ptr.min];
    llookup& l_ptr = k_ptr.p[hkl[i].k() - k_ptr.min];
    const int& l = hkl[i].l();  // find range of l
    l_ptr.min = Util::min( l_ptr.min, l );
    l_ptr.max = Util::max( l_ptr.max, l );
  }
  for (h=h_ptr.min; h<=h_ptr.max; h++) {
    klookup& k_ptr = h_ptr.p[h - h_ptr.min];
    for (k=k_ptr.min; k<=k_ptr.max; k++) {
      llookup& l_ptr = k_ptr.p[k - k_ptr.min];      // init all l_ptrs to -1
      if ( l_ptr.max >= l_ptr.min )
	l_ptr.p.resize( l_ptr.max - l_ptr.min + 1, -1 );
    }
  }

  // now fill in the data
  for (i=0; i<hkl.size(); i++) {
    klookup& k_ptr = h_ptr.p[hkl[i].h() - h_ptr.min];
    llookup& l_ptr = k_ptr.p[hkl[i].k() - k_ptr.min];
    l_ptr.p[ hkl[i].l() - l_ptr.min ] = i;
  }
}

int HKL_lookup::index_of(const HKL& rfl) const
{
  const int& h = rfl.h(), k = rfl.k(), l = rfl.l();
  // lookup a reflection:
  if ( h >= h_ptr.min && h <= h_ptr.max ) {
    const klookup& k_ptr = h_ptr.p[h - h_ptr.min];
    if ( k >= k_ptr.min && k <= k_ptr.max ) {
      const llookup& l_ptr = k_ptr.p[k - k_ptr.min];
      if ( l >= l_ptr.min && l <= l_ptr.max ) {
	return l_ptr.p[ l - l_ptr.min ];
      }
    }
  }
  // if not found, return -1
  return -1;
}

void HKL_lookup::debug()
{
  int size=0, h, k;
  for (h=h_ptr.min; h<=h_ptr.max; h++) {
    klookup& k_ptr = h_ptr.p[h - h_ptr.min];
    for (k=k_ptr.min; k<=k_ptr.max; k++) {
      llookup& l_ptr = k_ptr.p[k - k_ptr.min];
      size+=l_ptr.p.size();
    }
  }
  std::cout << "lookup: size " << size << "\n";
}


} // namespace clipper
