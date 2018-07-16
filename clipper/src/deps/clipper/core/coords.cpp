/* coords.cpp: fundamental data types for the clipper libraries */
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


#include "coords.h"

#include "spacegroup.h"
#include "rotation.h"


namespace clipper {


/*! \param resol_ The resolution limit in Angstroms. */
Resolution::Resolution( const ftype& resol_ )
{ init( resol_ ); }

/*! \param resol_ The resolution limit in Angstroms. */
void Resolution::init( const ftype& resol_ )
{ resol = resol_; }

/*! \return The resolution limit in Angstroms. */
const ftype& Resolution::limit() const
{ return resol; }

/*! \return The resolution limit in inverse squared Angstroms. */
ftype Resolution::invresolsq_limit() const
{ return 1.0/(resol*resol); }

/*! \return true if the object has not been initalised. */
bool Resolution::is_null() const
{ return ( resol <= 0.0 ); }


/*! Determine the class of a reflection for a give spacegroup.
  \param spgr The spacegroup.
  \param hkl The reflection HKL */
HKL_class::HKL_class( const Spacegroup& spgr, const HKL& hkl )
{
  HKL equiv;
  ftype shift;
  epsilon_ = 1;
  allowed_ = 255;
  for ( int i = 1; i < spgr.num_symops(); i++ ) {
    equiv = hkl.transform(spgr.symop(i));
    shift = hkl.sym_phase_shift(spgr.symop(i));
    if ( equiv == hkl ) {
      // if reflection is related to itself, then it is special or sysabs
      if ( cos(shift) > 0.999 ) {
	epsilon_++;  // increase multiplicity
      } else {
	epsilon_ = allowed_ = 0; // flag sysabs
	break;
      }
    } else if ( equiv == -hkl ) {
      // if reflection is related to opposite, it is centric
      allowed_ = Util::intr( Util::mod(-0.5*shift, Util::pi()) / (Util::pi()/12.0) );
    }
  }
  if ( hkl.h() == 0 && hkl.k() == 0 && hkl.l() == 0 ) allowed_ = 0;
}


/*! Construct the operator which give the least-squares fit of one set
  of coordinates onto another. The coodinates are stored as STL
  vectors of Coord_orth. The lists must be the same size, and each
  atom in the source list must correspond to the same atom in the
  target list.  The algorithm employed is that of Kearsley,
  S.K. (1989) 'On the orthogonal transformation used for structural
  comparisons'. Acta Cryst. A45, 208-210.
  \param src The source list (i.e. the atoms to be transformed).
  \param tgt The target list (i.e. the fixed atoms). */
RTop_orth::RTop_orth( const std::vector<Coord_orth>& src, const std::vector<Coord_orth>& tgt )
{
  if ( src.size() != tgt.size() )
    Message::message(Message_fatal("RTop_orth: coordinate list size mismatch"));
  // get centre of mass
  Coord_orth s, t, p, m;
  Coord_orth src_cen(0.0,0.0,0.0);
  Coord_orth tgt_cen(0.0,0.0,0.0);
  int n = src.size();
  for ( int i = 0; i < n; i++ ) {
    src_cen = src_cen + src[i];
    tgt_cen = tgt_cen + tgt[i];
  }
  src_cen = (1.0/ftype(n)) * src_cen;
  tgt_cen = (1.0/ftype(n)) * tgt_cen;
  // prepare cross-sums
  Matrix<> mat( 4, 4, 0.0 );
  for ( int i = 0; i < n; i++ ) {
    s = src[i] - src_cen;
    t = tgt[i] - tgt_cen;
    p = s + t;
    m = s - t;
    mat(0,0) = mat(0,0) + m[0]*m[0] + m[1]*m[1] + m[2]*m[2];
    mat(1,1) = mat(1,1) + m[0]*m[0] + p[1]*p[1] + p[2]*p[2];
    mat(2,2) = mat(2,2) + p[0]*p[0] + m[1]*m[1] + p[2]*p[2];
    mat(3,3) = mat(3,3) + p[0]*p[0] + p[1]*p[1] + m[2]*m[2];
    mat(0,1) = mat(0,1) + m[2]*p[1] - m[1]*p[2];
    mat(0,2) = mat(0,2) + m[0]*p[2] - m[2]*p[0];
    mat(0,3) = mat(0,3) + m[1]*p[0] - m[0]*p[1];
    mat(1,2) = mat(1,2) + m[0]*m[1] - p[0]*p[1];
    mat(1,3) = mat(1,3) + m[0]*m[2] - p[0]*p[2];
    mat(2,3) = mat(2,3) + m[1]*m[2] - p[1]*p[2];
  }
  mat(1,0) = mat(0,1);
  mat(2,0) = mat(0,2);
  mat(2,1) = mat(1,2);
  mat(3,0) = mat(0,3);
  mat(3,1) = mat(1,3);
  mat(3,2) = mat(2,3);
  // eigenvalue calc
  std::vector<ftype> v = mat.eigen( true );
  // get result
  Rotation r( mat(0,0), mat(1,0), mat(2,0), mat(3,0) );
  Mat33<> rot = r.norm().matrix();
  (*this) = RTop_orth( rot, tgt_cen - rot*src_cen );
}

/*! Construct the operator which give the least-squares fit of one set
  of coordinates onto another. The coodinates are stored as STL
  vectors of Coord_orth. The lists must be the same size, and each
  atom in the source list must correspond to the same atom in the
  target list.  The algorithm employed is that of Kearsley,
  S.K. (1989) 'On the orthogonal transformation used for structural
  comparisons'. Acta Cryst. A45, 208-210.
  \param src The source list (i.e. the atoms to be transformed).
  \param tgt The target list (i.e. the fixed atoms).
  \param wgt The weight to apply to each atom. */
  RTop_orth::RTop_orth( const std::vector<Coord_orth>& src, const std::vector<Coord_orth>& tgt, const std::vector<ftype>& wgt )
{
  if ( src.size() != wgt.size() || tgt.size() != wgt.size() )
    Message::message(Message_fatal("RTop_orth: coordinate list size mismatch"));
  // get centre of mass
  Coord_orth s, t, p, m;
  Coord_orth src_cen(0.0,0.0,0.0);
  Coord_orth tgt_cen(0.0,0.0,0.0);
  ftype w = 0.0;
  int n = src.size();
  for ( int i = 0; i < n; i++ ) {
    w += wgt[i];
    src_cen += wgt[i] * src[i];
    tgt_cen += wgt[i] * tgt[i];
  }
  src_cen = (1.0/w) * src_cen;
  tgt_cen = (1.0/w) * tgt_cen;
  // prepare cross-sums
  Matrix<> mat( 4, 4, 0.0 );
  for ( int i = 0; i < n; i++ ) {
    w = wgt[i];
    s = src[i] - src_cen;
    t = tgt[i] - tgt_cen;
    p = s + t;
    m = s - t;
    mat(0,0) = mat(0,0) + w * ( m[0]*m[0] + m[1]*m[1] + m[2]*m[2] );
    mat(1,1) = mat(1,1) + w * ( m[0]*m[0] + p[1]*p[1] + p[2]*p[2] );
    mat(2,2) = mat(2,2) + w * ( p[0]*p[0] + m[1]*m[1] + p[2]*p[2] );
    mat(3,3) = mat(3,3) + w * ( p[0]*p[0] + p[1]*p[1] + m[2]*m[2] );
    mat(0,1) = mat(0,1) + w * ( m[2]*p[1] - m[1]*p[2] );
    mat(0,2) = mat(0,2) + w * ( m[0]*p[2] - m[2]*p[0] );
    mat(0,3) = mat(0,3) + w * ( m[1]*p[0] - m[0]*p[1] );
    mat(1,2) = mat(1,2) + w * ( m[0]*m[1] - p[0]*p[1] );
    mat(1,3) = mat(1,3) + w * ( m[0]*m[2] - p[0]*p[2] );
    mat(2,3) = mat(2,3) + w * ( m[1]*m[2] - p[1]*p[2] );
  }
  mat(1,0) = mat(0,1);
  mat(2,0) = mat(0,2);
  mat(2,1) = mat(1,2);
  mat(3,0) = mat(0,3);
  mat(3,1) = mat(1,3);
  mat(3,2) = mat(2,3);
  // eigenvalue calc
  std::vector<ftype> v = mat.eigen( true );
  // get result
  Rotation r( mat(0,0), mat(1,0), mat(2,0), mat(3,0) );
  Mat33<> rot = r.norm().matrix();
  (*this) = RTop_orth( rot, tgt_cen - rot*src_cen );
}

/*! \param cell The cell concerned \return The transformed coordinate. */
RTop_frac RTop_orth::rtop_frac( const Cell& cell ) const
{
  return RTop_frac( RTop<>(cell.matrix_frac()) * (*this) * RTop<>(cell.matrix_orth()) );
}

/*! \return The inverse of the operator. */
RTop_orth RTop_orth::inverse() const
{ return RTop_orth( RTop<>::inverse() ); }

/*! \param centre An arbitrary point.
  \return point on axis near the specified coordinate, 000 if rotation is zero */
Coord_orth RTop_orth::axis_coordinate_near( const Coord_orth& centre ) const
{
  Rotation R(rot());     // Rotation part
  if ( R.abs_angle() < 0.001) return Coord_orth(0.0,0.0,0.0);
  Vec3<> t = trn();  // translation part
  // direction cosines of rotation axis n
  Vec3<> n = Vec3<>(R.x(),R.y(),R.z()).unit();
  // parallel component of translation r = n (n.t)
  Vec3<> r = n * Vec3<>::dot(n, t);
  // Perpendicular component s = t - r
  Vec3<> s = t - r;
  // Rotation angle kappa
  double kappa = R.norm().polar_ccp4().kappa();
  // Origin x0 = 0.5*[s + (n x s)/tan(kappa/2)]
  Vec3<> x0 = 0.5*(s + Vec3<>::cross(n, s) * (1./tan(0.5*kappa)));
  // Find closest position on axis to centre 
  //  projection on to axis through origin
  return Coord_orth( x0 + n * Vec3<>::dot((centre-x0), n) );
}

/*! \return screw translation, 000 if rotation is zero */
Coord_orth RTop_orth::screw_translation() const
{
  Rotation R(rot());     // Rotation part
  if ( R.abs_angle() < 0.001) return Coord_orth(0.0,0.0,0.0);
  Vec3<> t = trn();  // translation part
  // direction cosines of rotation axis n
  Vec3<> n = Vec3<>(R.x(),R.y(),R.z()).unit();
  // parallel component of translation r = n (n.t)
  return Coord_orth( n * Vec3<>::dot(n, t) );
}

/*! \return The identity operator. */
RTop_orth RTop_orth::identity()
{ return RTop_orth( RTop<>::identity() ); }

/*! \return The null (uninitialised) operator. */
RTop_orth RTop_orth::null()
{ return RTop_orth( RTop<>::null() ); }


/*! \return The formatted text string */
String HKL::format() const
{ return "HKL = ("+String(h())+","+String(k())+","+String(l())+")"; }


/*! \return The formatted text string */
String Coord_reci_orth::format() const
{ return "x*y*z* = ("+String(xs())+","+String(ys())+","+String(zs())+")"; }


/*! \return The formatted text string */
String Coord_reci_frac::format() const
{ return "u*v*w* = ("+String(us())+","+String(vs())+","+String(ws())+")"; }


/*! The coordinate is calculated which extends the sequence of
  coordinates x1, x2, x3 with the specified distance to x3, angle to
  x2,x3, and torsion to x1,x2,x3.
  \param x1 First coordinate.
  \param x2 Second coordinate.
  \param x3 Third coordinate.
  \param length x3-new bond length in Angstroms.
  \param angle x2-x3-new opening angle in Radians.
  \param torsion x1-x2-x3-new torsion angle in Radians. */
Coord_orth::Coord_orth( const Coord_orth& x1, const Coord_orth& x2, const Coord_orth& x3, const ftype& length, const ftype& angle, const ftype& torsion )
{
  const Coord_orth xa( (x3-x2).unit() );
  const Coord_orth xc( Vec3<>::cross( x2-x1, xa ).unit() );
  const Coord_orth xb( Vec3<>::cross( xa, xc ) );
  const ftype wa = -length*cos( angle );
  const ftype wb = -length*sin( angle )*cos( -torsion );
  const ftype wc = -length*sin( angle )*sin( -torsion );
  (*this) = x3 + wa*xa + wb*xb + wc*xc;
}

/*! \return The bond length x1-x2 in Angstroms. */
ftype Coord_orth::length( const Coord_orth& x1, const Coord_orth& x2)
{ return sqrt( (x2-x1).lengthsq() ); }

/*! \return The bond angle x1-x2-x3 in Radians. */
ftype Coord_orth::angle( const Coord_orth& x1, const Coord_orth& x2, 
			 const Coord_orth& x3)
{ return acos( (x3-x2).unit() * (x1-x2).unit() ); }

/*! \return The bond torsion x1-x2-x3-x4 in Radians. */
ftype Coord_orth::torsion( const Coord_orth& x1, const Coord_orth& x2,
			   const Coord_orth& x3, const Coord_orth& x4)
{
  const Vec3<> xu( (x3-x2).unit() );
  const Vec3<> xa( Vec3<>::cross( x2-x1, xu ) );
  const Vec3<> xb( Vec3<>::cross( xu, x4-x3 ) );
  const Vec3<> xc( Vec3<>::cross( xa, xb ) );
  return atan2( xc*xu, xa*xb );
}

/*! \return The formatted text string */
String Coord_orth::format() const
{ return "xyz = ("+String(x(),10,4)+","+String(y(),10,4)+","+String(z(),10,4)+")"; }


Coord_frac Coord_frac::symmetry_copy_near(const Spacegroup& spgr, const Cell& cell, const Coord_frac& n) const {
  Coord_frac c, cmin(*this);
  double d2, d2min(1.0e12);
  for ( int k = 0; k < spgr.num_symops(); k++ ) {
    c = spgr.symop(k) * (*this);
    c = c.lattice_copy_near( n );
    d2 = ( c - n ).lengthsq( cell );
    if ( d2 < d2min ) {
      d2min = d2;
      cmin = c;
    }
  }
  return cmin;
}


/*! \return The formatted text string */
String Coord_frac::format() const
{ return "uvw = ("+String(u(),10,4)+","+String(v(),10,4)+","+String(w(),10,4)+")"; }


/*! \return The formatted text string */
String Coord_grid::format() const
{ return "uvw = ("+String(u())+","+String(v())+","+String(w())+")"; }


/*! \return The formatted text string */
String Coord_map::format() const
{ return "uvw = ("+String(u(),10,4)+","+String(v(),10,4)+","+String(w(),10,4)+")"; }


/*! \return The formatted text string */
String Grid::format() const
{ return "Nuvw = ("+String(nu())+","+String(nv())+","+String(nw())+")"; }

void Grid::debug() const
{ std::cout << format() << "\n"; }

/*! Make a map grid with an oblong bounded by the coordinates min and max.
  \param min The lower bound coordinate in u,v,w.
  \param max The upper bound coordinate in u,v,w. */
Grid_range::Grid_range( const Coord_grid& min, const Coord_grid& max )
{
  min_ = min;
  max_ = max;
  (*this)[0] = max_.u()-min_.u()+1;
  (*this)[1] = max_.v()-min_.v()+1;
  (*this)[2] = max_.w()-min_.w()+1;
}

/*! Make a map grid with an oblong bounded by the fractional
  coordinates min and max, when the sampling of the cell is g
  \param g The grid sampling of the whole unit cell.
  \param min The lower bound coordinate in u,v,w. 
  \param max The upper bound coordinate in u,v,w. */
Grid_range::Grid_range( const Grid& g, const Coord_frac& min, const Coord_frac& max )
{
  min_ = Coord_grid( Util::intc( g.nu() * min.u() ),
		     Util::intc( g.nv() * min.v() ),
		     Util::intc( g.nw() * min.w() ) );
  max_ = Coord_grid( Util::intf( g.nu() * max.u() ),
		     Util::intf( g.nv() * max.v() ),
		     Util::intf( g.nw() * max.w() ) );
  (*this)[0] = max_.u()-min_.u()+1;
  (*this)[1] = max_.v()-min_.v()+1;
  (*this)[2] = max_.w()-min_.w()+1;
}

/*! Make a map grid large enough to fully enclose a sphere about the
  origin of a given radius with a given cell and grid sampling.
  \param cell The cell parameters.
  \param grid The grid sampling of the whole cell.
  \param radius The radius of the sphere in Angstroms. */
Grid_range::Grid_range( const Cell& cell, const Grid& grid, const ftype& radius )
{
  Coord_grid lim( Util::intc( radius * cell.a_star() * ftype(grid.nu()) ),
		  Util::intc( radius * cell.b_star() * ftype(grid.nv()) ),
		  Util::intc( radius * cell.c_star() * ftype(grid.nw()) ) );
  min_ = -lim;
  max_ =  lim;
  (*this)[0] = max_.u()-min_.u()+1;
  (*this)[1] = max_.v()-min_.v()+1;
  (*this)[2] = max_.w()-min_.w()+1;
}

/*! Enlarge the grid by adding \c b cells in every direction.
  Will shrink the grid if \c b is negative.
  \param b The number of cells by which to enlarge/shrink. */
void Grid_range::add_border( const int b )
{
  min_ = min_ - Coord_grid(b,b,b);
  max_ = max_ + Coord_grid(b,b,b);
  (*this)[0] = max_.u()-min_.u()+1;
  (*this)[1] = max_.v()-min_.v()+1;
  (*this)[2] = max_.w()-min_.w()+1;
}

/*! A grid is chosen to represent the specified cell at the given
  resolution, obeying any restrictions imposed by the spacegroup. A
  slightly finer grid may be chosen if doing so is liable to
  significantly increase the speed of FFTs on that grid.

  \param spacegroup The spacegroup which the grid must obey.
  \param cell The cell which the grid must contain.
  \param resol The resolution to which the grid must sample.
  \param rate The linear Shannon rate (oversampling) required. If rate
  = 1, the grid spaceing will be half the resolution (the the minimum
  required). For a grid spaceing of resol/3, use the default rate=1.5.
*/
Grid_sampling::Grid_sampling( const Spacegroup& spacegroup, const Cell& cell, const Resolution& resol, const ftype rate )
{
  init( spacegroup, cell, resol, rate );
}

/*! A grid is chosen to represent the specified cell at the given
  resolution, obeying any restrictions imposed by the spacegroup. A
  slightly finer grid may be chosen if doing so is liable to
  significantly increase the speed of FFTs on that grid.

  \param spacegroup The spacegroup which the grid must obey.
  \param cell The cell which the grid must contain.
  \param resol The resolution to which the grid must sample.
  \param rate The linear Shannon rate (oversampling) required. If rate
  = 1, the grid spaceing will be half the resolution (the the minimum
  required). For a grid spaceing of resol/3, use the default rate=1.5.
*/
void Grid_sampling::init( const Spacegroup& spacegroup, const Cell& cell, const Resolution& resol, const ftype rate )
{
  int i, j, l, n, m, nbest;
  bool eqxy, eqxz, eqyz;
  ftype t, tbest;

  // search symops to find grid factors and equalities
  Grid g( 48, 48, 48 );
  Grid factors( 1, 1, 1 );
  eqxy = eqxz = eqyz = false;
  for ( i = 0; i < spacegroup.num_symops(); i++ ) {
    Isymop isymop( spacegroup.symop(i), g );
    eqxy = eqxy || ( isymop.rot()(0,1) != 0 );
    eqxz = eqxz || ( isymop.rot()(0,2) != 0 );
    eqyz = eqyz || ( isymop.rot()(1,2) != 0 );
    for ( j = 0; j < 3; j++ )
      factors[j] = Util::max( factors[j], g[j]/(Util::mod(isymop.trn()[j]-1, g[j])+1) );
  }

  // now try a grid
  // first (lowest) estimate
  Grid_sampling
    nuvw( Util::intc( 2.0 * cell.descr().a() * rate / resol.limit() ),
	  Util::intc( 2.0 * cell.descr().b() * rate / resol.limit() ),
	  Util::intc( 2.0 * cell.descr().c() * rate / resol.limit() ) );
  // now check against restrictions, speed
  for ( i = 0; i < 3; i++ ) {
    nbest = 0;
    tbest = 1.0e12;
    for ( n = nuvw[i]; n < 2*nuvw[i]+16; n++ ) {
      if ( n % (2*factors[i]) == 0 ) {
        l = 0; // sum of factors (approx. log n)
        m = n; // what is left in factorisation
        for ( j = 2; j <= n; j++ )
          while ( m%j == 0 ) { m /= j; l += j; }
        // FFT time O( n * l ) - introduce an extra n^2 for 3D.
        t = pow( ftype(n), 3 ) * ftype(l);
        if ( t < tbest ) { nbest = n; tbest = t; }
      }
    }
    nuvw[i] = nbest;
  }

  // now check symmetry relationships
  if ( eqxy ) { n = Util::max( nuvw[0], nuvw[1] ); nuvw[0] = nuvw[1] = n; }
  if ( eqxz ) { n = Util::max( nuvw[0], nuvw[2] ); nuvw[0] = nuvw[2] = n; }
  if ( eqyz ) { n = Util::max( nuvw[1], nuvw[2] ); nuvw[1] = nuvw[2] = n; }
  (*this) = nuvw;
}

/*! The result is an RT operator. This is a redudent representation,
  but is handy for assembling compound operators.
  \return The operator */
Mat33<> Grid_sampling::matrix_grid_frac() const
{
  Mat33<> m( Mat33<>::identity() );
  m(0,0) = 1.0/ftype(nu());
  m(1,1) = 1.0/ftype(nv());
  m(2,2) = 1.0/ftype(nw());
  return m;
}

/*! The result is an RT operator. This is a redudent representation,
  but is handy for assembling compound operators.
  \return The operator */
Mat33<> Grid_sampling::matrix_frac_grid() const
{
  Mat33<> m( Mat33<>::identity() );
  m(0,0) = ftype(nu());
  m(1,1) = ftype(nv());
  m(2,2) = ftype(nw());
  return m;
}

/*! \return true if the object has not been initalised. */
bool Grid_sampling::is_null() const
{ return ( size() <= 0 ); }


/*! Threshold value for scaling HKL-sampling coefficients */
itype64 HKL_sampling::sqrt_limit_value = 0x100000;

/*! Initialise to 'null' */
HKL_sampling::HKL_sampling()
{
  m00 = -1;
}

/*! Initialise using cell and resolution. */
HKL_sampling::HKL_sampling( const Cell& cell, const Resolution& resolution )
{
  itype64 limit_value = sqrt_limit_value*sqrt_limit_value;
  ftype s = ftype( limit_value ) / resolution.invresolsq_limit();
  m00 = itype64( s * cell.a_star()*cell.a_star() );
  m11 = itype64( s * cell.b_star()*cell.b_star() );
  m22 = itype64( s * cell.c_star()*cell.c_star() );
  m01 = itype64( s * 2.0*cell.a_star()*cell.b_star()*cos(cell.gamma_star()) );
  m02 = itype64( s * 2.0*cell.a_star()*cell.c_star()*cos(cell.beta_star() ) );
  m12 = itype64( s * 2.0*cell.b_star()*cell.c_star()*cos(cell.alpha_star()) );
}

/*! Returned HKL contains maximum possible values of H, K, L respectively.
  \return Limiting h,k,l. */
HKL HKL_sampling::hkl_limit() const
{
  itype64 s00(m00/sqrt_limit_value), s11(m11/sqrt_limit_value),
          s22(m22/sqrt_limit_value), s01(m01/sqrt_limit_value),
          s02(m02/sqrt_limit_value), s12(m12/sqrt_limit_value);
  s01 /= 2; s02 /= 2; s12 /= 2;
  itype64 det = s00*s11*s22 + s01*s12*s02 + s02*s12*s01
              - s00*s12*s12 - s11*s02*s02 - s22*s01*s01;
  itype64 n00 = ( sqrt_limit_value * ( s11*s22 - s12*s12 + 1 ) ) / ( det - 3 );
  itype64 n11 = ( sqrt_limit_value * ( s00*s22 - s02*s02 + 1 ) ) / ( det - 3 );
  itype64 n22 = ( sqrt_limit_value * ( s00*s11 - s01*s01 + 1 ) ) / ( det - 3 );
  return HKL( Util::isqrt(n00), Util::isqrt(n11), Util::isqrt(n22) );
}

/*! Returned resolution is an estimate based on highest reflection in list.
  \return The resolution. */
Resolution HKL_sampling::resolution( const Cell& cell ) const
{
  HKL lim = hkl_limit();
  HKL rfl;
  // make a list of valid reflections
  ftype slim(0.0);
  for (rfl.h()=       0; rfl.h()<=lim.h(); rfl.h()++)
    for (rfl.k()=-lim.k(); rfl.k()<=lim.k(); rfl.k()++)
      for (rfl.l()=-lim.l(); rfl.l()<=lim.l(); rfl.l()++)
	if ( in_resolution( rfl ) )
	  slim = Util::max( slim, rfl.invresolsq( cell ) );
  return Resolution( 0.999999 / sqrt( slim ) );
}

bool HKL_sampling::is_null() const { return m00 < 0; }

String HKL_sampling::format() const
{
  return "CANNOT FORMAT HKL_sampling";
  //return "m00=" + String(m00,12) + " m11=" + String(m11,12) + " m22=" + String(m22,12) + " m01=" + String(m01,12) + " m02=" + String(m02,12) + " m12=" + String(m12,12);
}


/*! The best isotropic U is the cube root of the determinant of the
  matrix of anisotropic coefficients. NOTE: This is not the
  conventional definition, but the mathematically correct one, and
  gives a better approximation to the anisotropic U (i.e. lower
  R-factors).
  \return The nearest isotropic U. */
ftype U_aniso_orth::u_iso() const
{
  return pow( det(), 0.3333333333 );
}

/*! \param cell The cell concerned \return The transformed coordinate. */
U_aniso_frac U_aniso_orth::u_aniso_frac( const Cell& cell ) const
{ return U_aniso_frac( Mat33sym<>( cell.matrix_frac() * Mat33<>(*this) * cell.matrix_frac().transpose() ) ); }

/*! The aniso U is transformed by the given RT op.
  \param u The aniso U. */
U_aniso_orth U_aniso_orth::transform( const RTop_orth& op ) const
{
  Mat33<> r = op.rot().inverse();
  return U_aniso_orth( Mat33sym<> ( r.transpose() * Mat33<>( *this ) * r ) );
}

/*! \param cell The cell concerned \return The transformed coordinate. */
U_aniso_orth U_aniso_frac::u_aniso_orth( const Cell& cell ) const
{ return U_aniso_orth( Mat33sym<>( cell.matrix_orth() * Mat33<>(*this) * cell.matrix_orth().transpose() ) ); }

/*! The aniso U is transformed by the given RT op.
  \param u The aniso U. */
U_aniso_frac U_aniso_frac::transform( const RTop_frac& op ) const
{
  Mat33<> r = op.rot().inverse();
  return U_aniso_frac( Mat33sym<> ( r.transpose() * Mat33<>( *this ) * r ) );
}

void Atom::set_element( const String& s ) { element_ = s; }
void Atom::set_coord_orth( const Coord_orth& s ) { coord_orth_ = s; }
void Atom::set_occupancy( const ftype& s ) { occupancy_ = s; }
void Atom::set_u_iso( const ftype& s ) { u_iso_ = s; }
void Atom::set_u_aniso_orth( const U_aniso_orth& s ) { u_aniso_orth_ = s; }

/*! The coordinates and U_aniso_orth are transformed. The sigmas are not,
  since without the full variance-covariance matrix this transformation
  is impossible.
  \param rt The operator to apply.
*/
void Atom::transform( const RTop_orth rt ) {
  set_coord_orth( coord_orth().transform( rt ) );
  set_u_aniso_orth( u_aniso_orth().transform( rt ) );
}

Atom Atom::null() {
  Atom atom;
  atom.set_element( "" );
  atom.set_coord_orth( Coord_orth( Coord_orth::null() ) );
  atom.set_u_aniso_orth( U_aniso_orth( U_aniso_orth::null() ) );
  atom.set_occupancy( Util::nan() );
  atom.set_u_iso( Util::nan() );
  return atom;
}

} // namespace clipper
